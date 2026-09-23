!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module kdtree
!
! This module implements the k-d tree build
!    and associated tree walking routines
!
! :References:
!    Gafton & Rosswog (2011), MNRAS 418, 770-781
!    Benz, Bowers, Cameron & Press (1990), ApJ 348, 647-667
!    Dehnen (2000), ApJL 536, L39; Dehnen (2002), JCoPh 179, 27
!    Marcello (2017), AJ 154, 92 (angular momentum conservation)
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: allocutils, boundary, dim, dtypekdtree, io, kernel,
!   mpibalance, mpidomain, mpitree, mpiutils, part, sortutils, timing
!
 use dim,         only:maxp,ncellsmax,minpart,use_apr,use_sinktree,maxptmass,maxpsph
 use io,          only:nprocs
 use dtypekdtree, only:kdnode,lenfgrav
 use part,        only:ll,iphase,treecache,maxphase, &
                       apr_level,aprmassoftype

 implicit none

 integer, public,  allocatable :: inoderange(:,:)
 integer, public,  allocatable :: inodeparts(:)
 type(kdnode),     allocatable :: refinementnode(:)
 real,             allocatable :: fnode_branch(:,:)
 integer,          allocatable :: neighnodecount_branch(:)
 integer,          allocatable :: neighnode_branch(:,:)
 integer,          allocatable :: neighnodecache(:)
 integer,          allocatable :: neighnodecache_start(:)
 integer,          allocatable :: neighnodecache_count(:)
 real,             allocatable :: fnodecache(:,:)
!$omp threadprivate(fnode_branch,neighnode_branch,neighnodecount_branch)
!
!--tree parameters
!
 integer,          parameter, public :: irootnode    = 1
 character(len=1), parameter, public :: labelax(3)   = (/'x','y','z'/)
 integer,          parameter         :: maxdepth     = 64
 integer,          parameter         :: maxnodecache_local = 512
 integer,          parameter         :: maxneigh_per_node  = 16
 integer,          parameter         :: keydepth_max = 52
!
!--runtime options for this module
!
 real,    public  :: tree_accuracy    = 0.5
 logical, public  :: use_geosplit     = .true.
  logical, public  :: use_geosplit_fast = .false. ! key-sorted fast build for geosplit trees (experimental)
  logical, public  :: use_tree_renumber = .true. ! deterministic DFS-preorder node numbering
  logical, public  :: use_cache        = .true.
 integer, private :: kdpat(0:keydepth_max-1) ! split axis (0,1,2) per dyadic level
 integer(kind=8), private, allocatable :: kdkey(:)     ! dyadic sort key per particle slot
 integer(kind=8), private, allocatable :: kdkey_buf(:) ! radix sort scratch
 integer,         private, allocatable :: kdord(:)     ! permutation carried through the key sort
 integer,         private, allocatable :: kdord_buf(:) ! radix sort scratch
  real,            private, allocatable :: kdtmp(:,:)   ! scratch for parallel key-order permute
  integer(kind=8), private :: kdkeymax = 0              ! max key value (sets radix passes)
  integer,         private, allocatable :: kddone(:)    ! marks fallback-built nodes (skipped in Phase B)
  integer,         private, allocatable :: kdmap(:)     ! old->new ids for deterministic renumbering
  integer,         private, allocatable :: kdbkt_start(:) ! level bucket starts (Phase B)
  integer,         private, allocatable :: kdbkt_list(:)  ! node ids by level (Phase B)
 logical, private :: done_init_kdtree = .false.
 logical, private :: already_warned   = .false.
 integer, private :: numthreads

! Index of the last node in the local tree that has been copied to
! the global tree
 integer :: irefine
 integer :: itail_neigh = 0

 public :: allocate_kdtree, deallocate_kdtree
 public :: maketree, revtree, getneigh,getneigh_dual,kdnode,lenfgrav
 public :: maketreeglobal
 public :: empty_tree
 public :: compute_M2L,expand_fgrav_in_taylor_series
 integer, public :: maxlevel_indexed, maxlevel

 ! neighbour cache indices (xyzcache); imported with only: from dens/force
 integer, parameter, public :: ix=1, iy=2, iz=3, ih1=4, im=5, irho=6, izetaomega=7, isoftomega=8

 type kdbuildstack
    integer :: node
    integer :: parent
    integer :: level
    integer :: npnode
    real    :: xmin(3)
    real    :: xmax(3)
 end type kdbuildstack

 private

contains

subroutine allocate_kdtree
 use dim, only:mpi,maxp,ncellsmax
 use allocutils, only:allocate_array

 call allocate_array('inoderange', inoderange, 2, ncellsmax+1)
 call allocate_array('inodeparts', inodeparts, maxp)
 call allocate_array('kdkey', kdkey, maxp)
 call allocate_array('kdkey_buf', kdkey_buf, maxp)
 call allocate_array('kdord', kdord, maxp)
 call allocate_array('kdord_buf', kdord_buf, maxp)
 call allocate_array('kdtmp', kdtmp, 5, maxp)
 call allocate_array('kddone', kddone, ncellsmax+1)
 call allocate_array('kdmap', kdmap, ncellsmax+1)
 call allocate_array('kdbkt_start', kdbkt_start, maxdepth+3)
 call allocate_array('kdbkt_list', kdbkt_list, ncellsmax+1)
 if (mpi) call allocate_array('refinementnode', refinementnode, ncellsmax+1)
 call allocate_array('fnodecache', fnodecache, lenfgrav, ncellsmax+1)
 call allocate_array('neighnodecache',neighnodecache,ncellsmax*maxneigh_per_node)
 call allocate_array('neighnodecache_start',neighnodecache_start,ncellsmax+1)
 call allocate_array('neighnodecache_count',neighnodecache_count,ncellsmax+1)
!$omp parallel
 call allocate_array('neighnodecount_branch',neighnodecount_branch,maxdepth)
 call allocate_array('neighnode_branch',neighnode_branch,maxnodecache_local,maxdepth)
 call allocate_array('fnode_branch', fnode_branch, lenfgrav, maxdepth)
!$omp end parallel
 itail_neigh = 0

end subroutine allocate_kdtree

subroutine deallocate_kdtree
 use dim, only:mpi
 if (allocated(inoderange)) deallocate(inoderange)
 if (allocated(inodeparts)) deallocate(inodeparts)
 if (allocated(kdkey)) deallocate(kdkey)
 if (allocated(kdkey_buf)) deallocate(kdkey_buf)
 if (allocated(kdord)) deallocate(kdord)
 if (allocated(kdord_buf)) deallocate(kdord_buf)
 if (allocated(kdtmp)) deallocate(kdtmp)
 if (allocated(kddone)) deallocate(kddone)
 if (allocated(kdmap)) deallocate(kdmap)
 if (allocated(kdbkt_start)) deallocate(kdbkt_start)
 if (allocated(kdbkt_list)) deallocate(kdbkt_list)
 if (mpi .and. allocated(refinementnode)) deallocate(refinementnode)
 if (allocated(fnodecache)) deallocate(fnodecache)
 if (allocated(neighnodecache)) deallocate(neighnodecache)
 if (allocated(neighnodecache_start)) deallocate(neighnodecache_start)
 if (allocated(neighnodecache_count)) deallocate(neighnodecache_count)
!$omp parallel
 if (allocated(neighnode_branch)) deallocate(neighnode_branch)
 if (allocated(neighnodecount_branch)) deallocate(neighnodecount_branch)
 if (allocated(fnode_branch)) deallocate(fnode_branch)
!$omp end parallel

end subroutine deallocate_kdtree

!--------------------------------------------------------------------------------
!+
!  Routine to build the tree from scratch
!
!  Notes/To do:
!  -openMP parallelisation of maketree_stack (done - April 2013)
!  -test centre of mass vs. geometric centre based cell sizes (done - 2013)
!  -need analysis module that times (and checks) build_tree and tree walk
!   for a given dump
!  -test bottom-up vs. top-down neighbour search
!  -should we try to store tree structure with particle arrays?
!  -need to compute centre of mass and moments for each cell on the fly (c.f. revtree?)
!  -need to implement long-range gravitational interaction (done - May 2013)
!  -implement revtree routine to update tree w/out rebuilding (done - Sep 2015)
!+
!-------------------------------------------------------------------------------
subroutine maketree(node, xyzh, np, leaf_is_active, ncells, apr_tree, refinelevels,nptmass,xyzmh_ptmass)
 use io,   only:fatal,warning,iprint,iverbose
!$ use omp_lib
 type(kdnode),    intent(out)   :: node(:) !ncellsmax+1)
 integer,         intent(in)    :: np
 real,            intent(inout) :: xyzh(:,:)  ! inout because of boundary crossing
 integer,         intent(out)   :: leaf_is_active(:) !ncellsmax+1)
 integer(kind=8), intent(out)   :: ncells
 logical,         intent(in)    :: apr_tree
 integer,         intent(out),   optional :: refinelevels
 integer,         intent(in),    optional :: nptmass
 real,            intent(inout), optional :: xyzmh_ptmass(:,:)

 integer :: i,npnode,il,ir,istack,nl,nr,mymum
 integer :: nnode,minlevel,level,nqueue
 real :: xmini(3),xmaxi(3),xminl(3),xmaxl(3),xminr(3),xmaxr(3)
 integer, parameter :: istacksize = 512
 type(kdbuildstack), save :: stack(istacksize)
 !$omp threadprivate(stack)
 type(kdbuildstack) :: queue(istacksize)
!$ integer :: threadid
 integer :: npcounter
 logical :: wassplit,finished,sinktree
 character(len=10) :: string

 if (present(nptmass) .and. present(xyzmh_ptmass)) then
    sinktree = .true.
 endif

 itail_neigh = 0
 leaf_is_active = 0

 ir = 0
 il = 0
 nl = 0
 nr = 0
 wassplit = .false.
 finished = .false.

 ! construct root node, i.e. find bounds of all particles
 if (sinktree) then
    call construct_root_node(np,npcounter,irootnode,xmini,xmaxi,leaf_is_active,xyzh,xyzmh_ptmass,nptmass)
 else
    call construct_root_node(np,npcounter,irootnode,xmini,xmaxi,leaf_is_active,xyzh)
 endif

 if (inoderange(1,irootnode)==0 .or. inoderange(2,irootnode)==0 ) then
    call fatal('maketree','no particles or all particles dead/accreted')
 endif

! Fast key-sorted build for geosplit trees: dyadic cells, single
! property scan per node, no per-level partitioning. Properties are
! computed with the same routines as the standard build, so every
! stored node quantity is identical for a given particle set; only
! the split planes (topology) differ.
 if (use_geosplit .and. use_geosplit_fast .and. .not.apr_tree) then
    call maketree_fast(node,np,npcounter,xmini,xmaxi,leaf_is_active,ncells,refinelevels)
    return
 endif

! Put root node on top of stack
 ncells = 1
 maxlevel = 0
 minlevel = maxdepth - 1
 istack = 1

 ! maximum level where 2^k indexing can be used (thus avoiding critical sections)
 ! deeper than this we access cells via a stack as usual
 maxlevel_indexed = int(log(real(ncellsmax+1))/log(2.)) - 1

 ! default number of cells is the size of the `indexed' part of the tree
 ! this can be *increased* by building tree beyond indexed levels
 ! and is decreased afterwards according to the maximum depth actually reached
 if (.not. use_geosplit) ncells = 2**(maxlevel_indexed+1) - 1

 ! need to number of particles in node during build
 ! this is counted above to remove dead/accreted particles
 call push_onto_stack(queue(istack),irootnode,0,0,npcounter,xmini,xmaxi)

 if (.not.done_init_kdtree) then
    ! 1 thread for serial, overwritten when using OpenMP
    numthreads = 1

    ! get number of OpenMP threads
    !$omp parallel default(none) shared(numthreads)
!$  numthreads = omp_get_num_threads()
    !$omp end parallel
    done_init_kdtree = .true.
 endif

 nqueue = numthreads
 ! build using a queue to build level by level until number of nodes = number of threads
 over_queue: do while (istack  <  nqueue)
    ! if the tree finished while building the queue, then we should just return
    ! only happens for small particle numbers
    if (istack <= 0) then
       finished = .true.
       exit over_queue
    endif
    ! pop off front of queue
    call pop_off_stack(queue(1), istack, nnode, mymum, level, npnode, xmini, xmaxi)

    ! shuffle queue forward
    do i=1,istack
       queue(i) = queue(i+1)
    enddo

    ! construct node
    if (sinktree) then
       call construct_node(node(nnode), nnode, mymum, level, xmini, xmaxi, npnode, .true., &  ! construct in parallel
                           il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr, ncells, leaf_is_active, &
                           minlevel, maxlevel, wassplit, .false., apr_tree, xyzmh_ptmass)
    else
       call construct_node(node(nnode), nnode, mymum, level, xmini, xmaxi, npnode, .true., &  ! construct in parallel
                           il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr, ncells, leaf_is_active, &
                           minlevel, maxlevel, wassplit, .false., apr_tree)
    endif

    if (wassplit) then ! add children to back of queue
       if (istack+2 > istacksize) call fatal('maketree',&
                                       'queue size exceeded in tree build, increase istacksize and recompile')

       istack = istack + 1
       call push_onto_stack(queue(istack),il,nnode,level+1,nl,xminl,xmaxl)
       istack = istack + 1
       call push_onto_stack(queue(istack),ir,nnode,level+1,nr,xminr,xmaxr)
    endif

 enddo over_queue

 ! fix the indices

 done: if (.not.finished) then

    ! build using a stack which builds depth first
    ! each thread grabs a node from the queue and builds its own subtree

    !$omp parallel default(none) &
    !$omp shared(queue) &
    !$omp shared(ll, leaf_is_active) &
    !$omp shared(xyzmh_ptmass) &
    !$omp shared(np) &
    !$omp shared(node, ncells) &
    !$omp shared(nqueue,apr_tree,sinktree) &
    !$omp private(istack) &
    !$omp private(nnode, mymum, level, npnode, xmini, xmaxi) &
    !$omp private(ir, il, nl, nr) &
    !$omp private(xminr, xmaxr, xminl, xmaxl) &
    !$omp private(threadid) &
    !$omp private(wassplit) &
    !$omp reduction(min:minlevel) &
    !$omp reduction(max:maxlevel)
    !$omp do schedule(static)
    do i = 1, nqueue

       stack(1) = queue(i)
       istack = 1

       over_stack: do while(istack > 0)

          ! pop node off top of stack
          call pop_off_stack(stack(istack), istack, nnode, mymum, level, npnode, xmini, xmaxi)

          ! construct node
          if (sinktree) then
             call construct_node(node(nnode), nnode, mymum, level, xmini, xmaxi, npnode, .false., &  ! don't construct in parallel
                                 il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr, ncells, leaf_is_active, &
                                 minlevel, maxlevel, wassplit, .false., apr_tree, xyzmh_ptmass)
          else
             call construct_node(node(nnode), nnode, mymum, level, xmini, xmaxi, npnode, .false., &  ! don't construct in parallel
                                 il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr, ncells, leaf_is_active, &
                                 minlevel, maxlevel, wassplit, .false., apr_tree)
          endif

          if (wassplit) then ! add children to top of stack
             if (istack+2 > istacksize) call fatal('maketree',&
                                       'stack size exceeded in tree build, increase istacksize and recompile')

             istack = istack + 1
             call push_onto_stack(stack(istack),il,nnode,level+1,nl,xminl,xmaxl)
             istack = istack + 1
             call push_onto_stack(stack(istack),ir,nnode,level+1,nr,xminr,xmaxr)
          endif

       enddo over_stack
    enddo
    !$omp enddo
    !$omp end parallel

 endif done

 ! decrease number of cells if tree is entirely within 2^k indexing limit
 if ((maxlevel < maxlevel_indexed) .and. (.not. use_geosplit)) then
    ncells = 2**(maxlevel+1) - 1
 endif
 !-- if octree is used, we need to propagate information from leaf to root (hmax and quads)
 ! if (use_geosplit) call propagate_upward(int(ncells), node)

 if (maxlevel > maxlevel_indexed .and. .not.already_warned) then
    write(string,"(i10)") 2**(maxlevel-maxlevel_indexed)
    if (iverbose > 0) call warning('maketree','maxlevel > max_indexed: will run faster if recompiled with '// &
               'NCELLSMAX='//trim(adjustl(string))//'*maxp,')
 endif

  if (present(refinelevels)) refinelevels = minlevel

  ! deterministic DFS-preorder numbering (geosplit serial builds only:
  ! the 2^k indexed and MPI-refined paths rely on their own layouts)
  if (use_tree_renumber .and. use_geosplit .and. nprocs == 1) &
     call renumber_tree_dfs(node,ncells,leaf_is_active)

  if (iverbose >= 3) then
     write(iprint,"(a,i10,3(a,i2))") ' maketree: nodes = ',ncells,', max level = ',maxlevel,&
        ', min leaf level = ',minlevel,' max level indexed = ',maxlevel_indexed
  endif

end subroutine maketree

!--------------------------------------------------------------------------------
!+
!  Fast build for geosplit trees using a precomputed dyadic sort key.
!
!  The split cells are strictly dyadic (each child box is the geometric
!  half of its parent, no re-fitting to particle extents), so every
!  node's particles form a contiguous range in key order. The tree is
!  therefore assembled from key-range cuts with a single property scan
!  per node and no per-level partitioning. Node properties are computed
!  with compute_nodes_cofm/set_nodes_properties, i.e. identical values
!  to the standard build for a given particle set; only the split
!  planes (topology) differ.
!+
!-------------------------------------------------------------------------------
subroutine maketree_fast(node,np,nproot,xminroot,xmaxroot,leaf_is_active,ncells,refinelevels)
 use io,   only:fatal,warning,iprint,iverbose
!$ use omp_lib
 use dim, only:minpart
 use sortutils, only:radixsort_i8
 type(kdnode),    intent(out)   :: node(:)
 integer,         intent(in)    :: np,nproot
 real,            intent(in)    :: xminroot(3),xmaxroot(3)
 integer,         intent(out)   :: leaf_is_active(:)
 integer(kind=8), intent(out)   :: ncells
 integer,         intent(out),   optional :: refinelevels
 integer :: i,npnode,il,ir,istack,nl,nr,mymum
 integer :: nnode,minlevel,level,nqueue,npass,nlvl,i0,i1
 real :: xmini(3),xmaxi(3),xminl(3),xmaxl(3),xminr(3),xmaxr(3)
 real :: dumbox(3),dumbox2(3)
 integer, parameter :: istacksize = 512
 type(kdbuildstack), save :: stack(istacksize)
 !$omp threadprivate(stack)
 type(kdbuildstack) :: queue(istacksize)
!$ integer :: threadid
 logical :: wassplit,dofallback,finished,dobig
 character(len=10) :: string

 itail_neigh = 0
 leaf_is_active = 0
 kddone = 0
 dumbox = 0.
 dumbox2 = 0.

 ncells = 1
 maxlevel = 0
 minlevel = maxdepth - 1
 finished = .false.
 maxlevel_indexed = int(log(real(ncellsmax+1))/log(2.)) - 1

 if (.not.done_init_kdtree) then
    numthreads = 1
    !$omp parallel default(none) shared(numthreads)
!$  numthreads = omp_get_num_threads()
    !$omp end parallel
    done_init_kdtree = .true.
 endif

 ! dyadic split-axis pattern from the root box (data independent)
 call get_dyadic_axis_pattern(xminroot,xmaxroot)

 ! one sort replaces all per-level partitioning
 call compute_dyadic_keys(nproot,xminroot,xmaxroot)
 ! number of radix digits from the largest key actually present
 npass = 1
 do while (ishft(kdkeymax,-8*npass) > 0_8 .and. npass < 8)
    npass = npass + 1
 enddo
 call radixsort_i8(nproot,kdkey,kdord,kdkey_buf,kdord_buf,npass)
 call permute_to_key_order(nproot)

 ! root range (already set by construct_root_node, ensured here)
 inoderange(1,irootnode) = 1
 inoderange(2,irootnode) = nproot
 istack = 1
 call push_onto_stack(queue(istack),irootnode,0,0,nproot,xminroot,xmaxroot)

 nqueue = numthreads
 ! Phase A: structure only. Build level by level until number of nodes
 ! equals number of threads; no property scans, leaves just stop here
 over_queue: do while (istack < nqueue)
    if (istack <= 0) then
       finished = .true.
       exit over_queue
    endif
    call pop_off_stack(queue(1), istack, nnode, mymum, level, npnode, xmini, xmaxi)
    do i=1,istack
       queue(i) = queue(i+1)
    enddo
    call try_split_fast(node,nnode,level,xmini,xmaxi, &
                        il,ir,nl,nr,xminl,xmaxl,xminr,xmaxr,ncells,leaf_is_active, &
                        wassplit,dofallback)
    if (dofallback) then
       call build_old_subtree(node,nnode,mymum,level,npnode,xmini,xmaxi,ncells, &
                              leaf_is_active,minlevel,maxlevel,.true.)
    elseif (wassplit) then
       if (istack+2 > istacksize) call fatal('maketree_fast',&
           'queue size exceeded in tree build, increase istacksize and recompile')
       maxlevel = max(level+1,maxlevel)
       istack = istack + 1
       call push_onto_stack(queue(istack),il,nnode,level+1,nl,xminl,xmaxl)
       istack = istack + 1
       call push_onto_stack(queue(istack),ir,nnode,level+1,nr,xminr,xmaxr)
    else
       maxlevel = max(level,maxlevel)
       minlevel = min(level,minlevel)
    endif
 enddo over_queue

 done: if (.not.finished) then
    ! each thread splits its subtree fully (structure only)
    !$omp parallel default(none) &
    !$omp shared(queue) &
    !$omp shared(ll, leaf_is_active) &
    !$omp shared(node, ncells) &
    !$omp shared(nqueue) &
    !$omp private(istack) &
    !$omp private(nnode, mymum, level, npnode, xmini, xmaxi) &
    !$omp private(ir, il, nl, nr) &
    !$omp private(xminr, xmaxr, xminl, xmaxl) &
    !$omp private(threadid) &
    !$omp private(wassplit,dofallback) &
    !$omp reduction(min:minlevel) &
    !$omp reduction(max:maxlevel)
    !$omp do schedule(static)
    do i = 1, nqueue
       stack(1) = queue(i)
       istack = 1
       over_stack: do while(istack > 0)
          call pop_off_stack(stack(istack), istack, nnode, mymum, level, npnode, xmini, xmaxi)
          call try_split_fast(node,nnode,level,xmini,xmaxi, &
                              il,ir,nl,nr,xminl,xmaxl,xminr,xmaxr,ncells,leaf_is_active, &
                              wassplit,dofallback)
          if (dofallback) then
             call build_old_subtree(node,nnode,mymum,level,npnode,xmini,xmaxi,ncells, &
                                    leaf_is_active,minlevel,maxlevel,.false.)
          elseif (wassplit) then
             if (istack+2 > istacksize) call fatal('maketree_fast',&
                 'stack size exceeded in tree build, increase istacksize and recompile')
             maxlevel = max(level+1,maxlevel)
             istack = istack + 1
             call push_onto_stack(stack(istack),il,nnode,level+1,nl,xminl,xmaxl)
             istack = istack + 1
             call push_onto_stack(stack(istack),ir,nnode,level+1,nr,xminr,xmaxr)
          else
             maxlevel = max(level,maxlevel)
             minlevel = min(level,minlevel)
          endif
       enddo over_stack
    enddo
    !$omp enddo
    !$omp end parallel
 endif done

 ! bucket nodes by level for the bottom-up property pass
 call bucket_nodes_by_level(node,int(ncells))

 ! Phase B: properties bottom-up, one level at a time (deepest first).
 ! Nodes on the same level are independent. Levels with fewer nodes
 ! than threads hold the biggest nodes: run those serially so each
 ! scan gets a full thread team (as in the old queue build).
 over_levels: do nlvl = maxlevel,0,-1
    i0 = kdbkt_start(nlvl+1)
    i1 = kdbkt_start(nlvl+2)-1
    if (i1 < i0) cycle
    dobig = (i1-i0+1 < numthreads)
    ! leaf test by range size (exact restatement of the split rule);
    ! child/parent pointers of unscanned nodes are not valid yet
    if (dobig) then
       ! few (big) nodes: serial loop, each scan gets a full thread team
       do i = i0,i1
          nnode = kdbkt_list(i)
          if (kddone(nnode) == 1) cycle
          mymum = node(nnode)%parent
          level = node(nnode)%level
          npnode = inoderange(2,nnode)-inoderange(1,nnode)+1
          if (npnode <= minpart) then
             call scan_leaf_fast(node(nnode),nnode,mymum,level,dumbox,dumbox2,.true., &
                                 leaf_is_active,minlevel,maxlevel)
          else
             call scan_internal_fast(node,nnode,mymum,level,dumbox,dumbox2,.true.)
          endif
       enddo
    else
       !$omp parallel do default(none) schedule(guided) &
       !$omp shared(kdbkt_list,node,inoderange,kddone,leaf_is_active,i0,i1,dumbox,dumbox2) &
       !$omp private(i,nnode,mymum,level,npnode) &
       !$omp reduction(min:minlevel) &
       !$omp reduction(max:maxlevel)
       do i = i0,i1
          nnode = kdbkt_list(i)
          if (kddone(nnode) == 1) cycle
          mymum = node(nnode)%parent
          level = node(nnode)%level
          npnode = inoderange(2,nnode)-inoderange(1,nnode)+1
          if (npnode <= minpart) then
             call scan_leaf_fast(node(nnode),nnode,mymum,level,dumbox,dumbox2,.false., &
                                 leaf_is_active,minlevel,maxlevel)
          else
             call scan_internal_fast(node,nnode,mymum,level,dumbox,dumbox2,.false.)
          endif
       enddo
       !$omp end parallel do
    endif
 enddo over_levels

 ! deterministic DFS-preorder numbering (serial builds only: the MPI
 ! global-tree refinement indexes the local tree by level ranges)
 if (use_tree_renumber .and. nprocs == 1) call renumber_tree_dfs(node,ncells,leaf_is_active)

 if (maxlevel > maxlevel_indexed .and. .not.already_warned) then
    write(string,"(i10)") 2**(maxlevel-maxlevel_indexed)
    if (iverbose > 0) call warning('maketree_fast','maxlevel > max_indexed: will run faster if recompiled with '// &
               'NCELLSMAX='//trim(adjustl(string))//'*maxp,')
 endif

 if (present(refinelevels)) refinelevels = minlevel

 if (iverbose >= 3) then
    write(iprint,"(a,i10,3(a,i2))") ' maketree_fast: nodes = ',ncells,', max level = ',maxlevel,&
       ', min leaf level = ',minlevel,' max level indexed = ',maxlevel_indexed
 endif

end subroutine maketree_fast

!--------------------------------------------------------------------
!+
!  Axis pattern for the dyadic subdivision: at each level split the
!  longest axis of the current (geometric) box in half, exactly as
!  maxloc(xmax-xmin) would select it during the build. Pure geometry,
!  hence data independent and precomputable from the root box.
!+
!--------------------------------------------------------------------
subroutine get_dyadic_axis_pattern(xminroot,xmaxroot)
 real, intent(in) :: xminroot(3),xmaxroot(3)
 real :: blen(3)
 integer :: l

 blen(:) = xmaxroot(:) - xminroot(:)
 do l = 0,keydepth_max-1
    kdpat(l) = maxloc(blen,1) - 1
    blen(kdpat(l)+1) = 0.5*blen(kdpat(l)+1)
 enddo

end subroutine get_dyadic_axis_pattern

!--------------------------------------------------------------------
!+
!  Compute the dyadic (morton-like) sort key of each particle: bit
!  (keydepth_max-1-l) is 0/1 depending on which half of the level-l
!  dyadic cell the particle falls in along the pattern axis. The
!  sorted key order therefore groups every dyadic cell contiguously.
!+
!--------------------------------------------------------------------
subroutine compute_dyadic_keys(nproot,xminroot,xmaxroot)
 integer, intent(in) :: nproot
 real,    intent(in) :: xminroot(3),xmaxroot(3)
 integer :: i,a,l,nbits(0:2),cursor(0:2)
 integer(kind=8) :: ix(0:2),k,imax
 real :: sca(0:2),ext
 real :: x0(3)

 nbits(:) = 0
 do l = 0,keydepth_max-1
    nbits(kdpat(l)) = nbits(kdpat(l)) + 1
 enddo
 do a = 0,2
    ext = xmaxroot(a+1) - xminroot(a+1)
    if (ext > 0. .and. nbits(a) > 0) then
       sca(a) = (2.0**nbits(a))/ext
    else
       sca(a) = 0.
    endif
 enddo

 kdkeymax = 0_8
 !$omp parallel do default(none) schedule(static) &
 !$omp shared(nproot,xminroot,sca,nbits,treecache,kdkey,kdord,kdpat) &
 !$omp private(i,a,l,ix,cursor,k,imax,x0) &
 !$omp reduction(max:kdkeymax)
 do i = 1,nproot
    x0(1) = treecache(1,i)
    x0(2) = treecache(2,i)
    x0(3) = treecache(3,i)
    do a = 0,2
       if (sca(a) > 0.) then
          ix(a) = int((x0(a+1)-xminroot(a+1))*sca(a),kind=8)
          imax = ishft(1_8,nbits(a)) - 1_8
          if (ix(a) < 0_8) ix(a) = 0_8
          if (ix(a) > imax) ix(a) = imax
       else
          ix(a) = 0_8
       endif
       cursor(a) = nbits(a)
    enddo
    k = 0_8
    do l = 0,keydepth_max-1
       a = kdpat(l)
       cursor(a) = cursor(a) - 1
       if (btest(ix(a),cursor(a))) k = ibset(k,keydepth_max-1-l)
    enddo
    kdkey(i) = k
    if (k > kdkeymax) kdkeymax = k
    kdord(i) = i
 enddo
 !$omp end parallel do

end subroutine compute_dyadic_keys

!--------------------------------------------------------------------
!+
!  Permute inodeparts/treecache into key order. treecache rows go
!  through the kdtmp scratch buffer (parallel streaming copy, then
!  parallel gather); inodeparts follows by serial cycles (cheap).
!+
!--------------------------------------------------------------------
subroutine permute_to_key_order(nproot)
 integer, intent(in) :: nproot
 integer :: i,k

 ! kdord_buf is free (radix scratch): save inodeparts through it
 !$omp parallel do default(none) schedule(static) shared(nproot,treecache,kdtmp,inodeparts,kdord_buf) private(i)
 do i = 1,nproot
    kdtmp(:,i) = treecache(:,i)
    kdord_buf(i) = inodeparts(i)
 enddo
 !$omp parallel do default(none) schedule(static) &
 !$omp shared(nproot,treecache,kdtmp,inodeparts,kdord,kdord_buf) private(i,k)
 do i = 1,nproot
    k = kdord(i)
    treecache(:,i) = kdtmp(:,k)
    inodeparts(i) = kdord_buf(k)
 enddo

end subroutine permute_to_key_order

!--------------------------------------------------------------------
!+
!  First index in kdkey(lo:hi) with bitpos set (keys are sorted, so
!  the bit is monotone 0..1 within any dyadic cell range). Returns
!  hi+1 if no key in the range has the bit set.
!+
!--------------------------------------------------------------------
pure integer function key_lower_bound(lo,hi,bitpos)
 integer, intent(in) :: lo,hi,bitpos
 integer :: a,b,mid

 a = lo
 b = hi + 1
 do while (a < b)
    mid = (a + b)/2
    if (btest(kdkey(mid),bitpos)) then
       b = mid
    else
       a = mid + 1
    endif
 enddo
 key_lower_bound = a

end function key_lower_bound

!--------------------------------------------------------------------
!+
!  Fast version of construct_node for dyadic cells: identical node
!  properties (same routines), but children ranges come from key
!  cuts instead of partitioning, and child boxes are the geometric
!  halves of the parent box instead of particle extents. Sets
!  dofallback=.true. (no children created) if the split would exceed
!  the key depth; the caller then builds that subtree the old way.
!+
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
!  Split a dyadic node by key-range cut (no scans, no partitioning):
!  children ranges come from the level bit in the sorted keys, child
!  boxes are the geometric halves of the parent box. Sets
!  dofallback=.true. (no children created) if the split would exceed
!  the key depth; the caller then builds that subtree the old way.
!+
!--------------------------------------------------------------------
subroutine try_split_fast(node,nnode,level,xmini,xmaxi, &
                          il,ir,nl,nr,xminl,xmaxl,xminr,xmaxr,ncells,leaf_is_active, &
                          wassplit,dofallback)
 use dim, only:minpart
 use io,  only:fatal
 type(kdnode),    intent(inout) :: node(:)
 integer,         intent(in)    :: nnode,level
 real,            intent(in)    :: xmini(3),xmaxi(3)
 integer,         intent(out)   :: il,ir,nl,nr
 real,            intent(out)   :: xminl(3),xmaxl(3),xminr(3),xmaxr(3)
 integer(kind=8), intent(inout) :: ncells
 integer,         intent(inout) :: leaf_is_active(:)
 logical,         intent(out)   :: wassplit,dofallback
 integer(kind=8) :: myslot
 integer :: npnode,iaxis,bitpos,m

 ir = 0
 il = 0
 nl = 0
 nr = 0
 dofallback = .false.
 wassplit = .false.

 ! stamp the level now: Phase B buckets nodes before any scan
 node(nnode)%level = level

 if (inoderange(1,nnode) > 0) then
    npnode = inoderange(2,nnode) - inoderange(1,nnode) + 1
 else
    npnode = 0
 endif
 if (npnode < 1) return
 wassplit = (npnode > minpart)
 if (.not. wassplit) return

 ! no key bit left to split on: build this subtree the old way
 if (level+1 > keydepth_max-1) then
    dofallback = .true.
    return
 endif
 if (level+1 > maxdepth) call fatal('maketree_fast','maximum tree depth reached !!')
 iaxis  = kdpat(level) + 1
 bitpos = keydepth_max-1-level
 m = key_lower_bound(inoderange(1,nnode),inoderange(2,nnode),bitpos)
 if (m == inoderange(1,nnode) .or. m == inoderange(2,nnode)+1) then
    ! all particles on one side: balanced split as in construct_node
    m = inoderange(1,nnode) + npnode/2
 endif
 !$omp atomic capture
 ncells = ncells + 2
 myslot = ncells
 !$omp end atomic
 ir = int(myslot)
 il = int(myslot-1)
 if (ir > ncellsmax) call fatal('maketree_fast',&
    'number of nodes exceeds array dimensions, increase ncellsmax and recompile',ival=int(ncellsmax))
 node(nnode)%leftchild  = il
 node(nnode)%rightchild = ir
 ! stamp linkage now: Phase B reads parent/level straight from the nodes
 node(il)%parent = nnode
 node(ir)%parent = nnode
 node(il)%level = level+1
 node(ir)%level = level+1

 leaf_is_active(nnode) = 0

 inoderange(1,il) = inoderange(1,nnode)
 inoderange(2,il) = m-1
 inoderange(1,ir) = m
 inoderange(2,ir) = inoderange(2,nnode)
 nl = m - inoderange(1,nnode)
 nr = inoderange(2,nnode) - m + 1

 xminl(:) = xmini(:)
 xmaxl(:) = xmaxi(:)
 xmaxl(iaxis) = 0.5*(xmini(iaxis) + xmaxi(iaxis))
 xminr(:) = xmini(:)
 xmaxr(:) = xmaxi(:)
 xminr(iaxis) = 0.5*(xmini(iaxis) + xmaxi(iaxis))

end subroutine try_split_fast

!--------------------------------------------------------------------
!+
!  Full property scan for a leaf node (identical values to the
!  standard build for the same particle set).
!+
!--------------------------------------------------------------------
subroutine scan_leaf_fast(nodeentry,nnode,mymum,level,xmini,xmaxi,doparallel, &
                          leaf_is_active,minlevel,maxlevel)
 use dim, only:ind_timesteps
 use io,  only:fatal
 type(kdnode),    intent(out)   :: nodeentry
 integer,         intent(in)    :: nnode,mymum,level
 real,            intent(inout) :: xmini(3),xmaxi(3)
 logical,         intent(in)    :: doparallel
 integer,         intent(inout) :: leaf_is_active(:)
 integer,         intent(inout) :: minlevel,maxlevel
 real    :: xyzcofm(3),totmass_node
 integer :: npnode,i
 logical :: nodeisactive

 nodeisactive = .false.
 if (inoderange(1,nnode) > 0) then
    do i = inoderange(1,nnode),inoderange(2,nnode)
       if (inodeparts(i) > 0) then
          nodeisactive = .true.
          exit
       endif
    enddo
    npnode = inoderange(2,nnode) - inoderange(1,nnode) + 1
 else
    npnode = 0
 endif

 call compute_nodes_cofm(npnode,nnode,xyzcofm,totmass_node,doparallel)
 if (totmass_node <= 0.) call fatal('maketree_fast','totmass_node==0',val=totmass_node)

 ! geosplit keeps the centre of mass at the node centre (as in construct_node)
 call set_nodes_properties(npnode,nnode,xyzcofm,totmass_node,mymum,nodeentry,xmini,xmaxi, &
                           level,.false.,doparallel,.true.)

 nodeentry%leftchild  = 0
 nodeentry%rightchild = 0
 maxlevel = max(level,maxlevel)
 minlevel = min(level,minlevel)
 if (ind_timesteps) then
    if (nodeisactive) then
       leaf_is_active(nnode) = 1
    else
       leaf_is_active(nnode) = -1
    endif
 else
    leaf_is_active(nnode) = 1
 endif

end subroutine scan_leaf_fast

!--------------------------------------------------------------------
!+
!  Property scan for an internal node whose children are complete:
!  centre of mass in O(1) from the children's masses and centres,
!  then a single scan for moments/size/hmax. Children ranges and
!  child pointers must already be set.
!+
!--------------------------------------------------------------------
subroutine scan_internal_fast(node,nnode,mymum,level,xmini,xmaxi,doparallel)
 use io, only:fatal
 type(kdnode),    intent(inout) :: node(:)
 integer,         intent(in)    :: nnode,mymum,level
 real,            intent(inout) :: xmini(3),xmaxi(3)
 logical,         intent(in)    :: doparallel
 integer :: il,ir,npnode
 real :: xyzcofm(3),totmass_node,ml,mr

 il = node(nnode)%leftchild
 ir = node(nnode)%rightchild
 ml = node(il)%mass
 mr = node(ir)%mass
 totmass_node = ml + mr
 if (totmass_node <= 0.) call fatal('maketree_fast','totmass_node==0 in internal combine')
 xyzcofm(:) = (ml*node(il)%xcen(:) + mr*node(ir)%xcen(:))/totmass_node
 npnode = inoderange(2,nnode) - inoderange(1,nnode) + 1

 call set_nodes_properties(npnode,nnode,xyzcofm,totmass_node,mymum,node(nnode),xmini,xmaxi, &
                           level,.false.,doparallel,.true.)

end subroutine scan_internal_fast

!--------------------------------------------------------------------
!+
!  Build the subtree rooted at nnode with the standard (partitioning)
!  construct_node. Used for nodes that would split beyond the key
!  depth, operating on their private contiguous slice only.
!+
!--------------------------------------------------------------------
subroutine build_old_subtree(node,nnode,mymum,level,npnode,xmini,xmaxi,ncells, &
                             leaf_is_active,minlevel,maxlevel,doparallel)
 use io, only:fatal
 type(kdnode),    intent(inout) :: node(:)
 integer,         intent(in)    :: nnode,mymum,level,npnode
 real,            intent(in)    :: xmini(3),xmaxi(3)
 integer(kind=8), intent(inout) :: ncells
 integer,         intent(inout) :: leaf_is_active(:)
 integer,         intent(inout) :: minlevel,maxlevel
 logical,         intent(in)    :: doparallel
 integer, parameter :: istacksize = 512
 type(kdbuildstack) :: ostack(istacksize)
 integer :: iold,onode,omum,olev,onp,oil,oir,onl,onr
 real :: oxmini(3),oxmaxi(3),oxminl(3),oxmaxl(3),oxminr(3),oxmaxr(3)
 logical :: owas

 call push_onto_stack(ostack(1),nnode,mymum,level,npnode,xmini,xmaxi)
 kddone(nnode) = 1
 iold = 1
 over_old: do while (iold > 0)
    call pop_off_stack(ostack(iold),iold,onode,omum,olev,onp,oxmini,oxmaxi)
    call construct_node(node(onode),onode,omum,olev,oxmini,oxmaxi,onp,doparallel, &
                        oil,oir,onl,onr,oxminl,oxmaxl,oxminr,oxmaxr,ncells,leaf_is_active, &
                        minlevel,maxlevel,owas,.false.,.false.)
    if (owas) then
       if (iold+2 > istacksize) call fatal('maketree_fast',&
           'stack size exceeded in fallback tree build, increase istacksize and recompile')
       kddone(oil) = 1
       kddone(oir) = 1
       iold = iold + 1
       call push_onto_stack(ostack(iold),oil,onode,olev+1,onl,oxminl,oxmaxl)
       iold = iold + 1
       call push_onto_stack(ostack(iold),oir,onode,olev+1,onr,oxminr,oxmaxr)
    endif
 enddo over_old

end subroutine build_old_subtree

!--------------------------------------------------------------------
!+
!  Bucket node ids 1..ncells by tree level (serial counting sort,
!  trivial cost next to particle scans) for the bottom-up Phase B.
!  Slots of level L run from kdbkt_start(L+1) to kdbkt_start(L+2)-1.
!+
!--------------------------------------------------------------------
subroutine bucket_nodes_by_level(node,ncells)
 type(kdnode), intent(in) :: node(:)
 integer,      intent(in) :: ncells
 integer :: id,lev,pos,tmp

 kdbkt_start = 0
 do id = 1,ncells
    lev = node(id)%level
    if (lev < 0) lev = 0
    if (lev > maxdepth) lev = maxdepth
    kdbkt_start(lev+1) = kdbkt_start(lev+1) + 1
 enddo
 pos = 1
 do lev = 0,maxdepth
    tmp = kdbkt_start(lev+1)
    kdbkt_start(lev+1) = pos
    pos = pos + tmp
 enddo
 kdbkt_start(maxdepth+2) = pos
 do id = 1,ncells
    lev = node(id)%level
    if (lev < 0) lev = 0
    if (lev > maxdepth) lev = maxdepth
    kdbkt_list(kdbkt_start(lev+1)) = id
    kdbkt_start(lev+1) = kdbkt_start(lev+1) + 1
 enddo
 do lev = maxdepth,0,-1
    kdbkt_start(lev+2) = kdbkt_start(lev+1)
 enddo
 kdbkt_start(1) = 1

end subroutine bucket_nodes_by_level

!--------------------------------------------------------------------
!+
!  Deterministic DFS-preorder renumbering of the tree (serial builds
!  only): each subtree becomes contiguous in memory, which matches
!  the top-down DFS order of the tree walks. Root keeps id 1.
!  Permutes node/inoderange/leaf_is_active consistently and rewrites
!  all child/parent pointers. O(ncells), no extra memory beyond kdmap
!  (kdbkt_list is reused as the inverse map scratch).
!+
!--------------------------------------------------------------------
subroutine renumber_tree_dfs(node,ncells,leaf_is_active)
 use io, only:fatal
 type(kdnode), intent(inout) :: node(:)
 integer(kind=8), intent(in) :: ncells
 integer,         intent(inout) :: leaf_is_active(:)
 integer :: n,newid,sp,top,old,new,k,j
 integer :: st(1024)
 type(kdnode) :: ntmp
 integer :: rtmp(2),itmp

 ! 1. left-first traversal from the root: kdmap(old) = new
 newid = 0
 sp = 1
 st(1) = irootnode
 do while (sp > 0)
    n = st(sp)
    sp = sp - 1
    if (n <= 0 .or. n > ncells) call fatal('renumber_tree_dfs','invalid node id in traversal')
    newid = newid + 1
    kdmap(n) = newid
    if (node(n)%rightchild /= 0) then
       sp = sp + 1
       if (sp > 1024) call fatal('renumber_tree_dfs','traversal stack exceeded')
       st(sp) = node(n)%rightchild
    endif
    if (node(n)%leftchild /= 0) then
       sp = sp + 1
       if (sp > 1024) call fatal('renumber_tree_dfs','traversal stack exceeded')
       st(sp) = node(n)%leftchild
    endif
 enddo
 if (newid /= ncells) call fatal('renumber_tree_dfs','unreachable nodes: cannot renumber')

 ! 2. rewrite pointers to new ids (in place, map is complete)
 do old = 1,int(ncells)
    if (node(old)%leftchild /= 0) node(old)%leftchild = kdmap(node(old)%leftchild)
    if (node(old)%rightchild /= 0) node(old)%rightchild = kdmap(node(old)%rightchild)
    if (node(old)%parent /= 0) node(old)%parent = kdmap(node(old)%parent)
 enddo

 ! 3. inverse map into bucket scratch: result(new) = arrays(old=inv(new))
 do old = 1,int(ncells)
    kdbkt_list(kdmap(old)) = old
 enddo

 ! 4. in-place cycle permutation of node/inoderange/leaf_is_active
 !    (result(new) = arrays(invmap(new)); rotate each cycle)
 top = int(ncells)
 do new = 1,top
    if (kdbkt_list(new) <= 0) cycle
    j = new
    ntmp = node(j)
    rtmp = inoderange(:,j)
    itmp = leaf_is_active(j)
    do
       k = kdbkt_list(j)
       if (k == new) exit
       node(j) = node(k)
       inoderange(:,j) = inoderange(:,k)
       leaf_is_active(j) = leaf_is_active(k)
       kdbkt_list(j) = -k
       j = k
    enddo
    node(j) = ntmp
    inoderange(:,j) = rtmp
    leaf_is_active(j) = itmp
    kdbkt_list(j) = -abs(kdbkt_list(j))
    kdbkt_list(new) = -abs(kdbkt_list(new))
 enddo

end subroutine renumber_tree_dfs

!----------------------------
!+
! routine to empty the tree
!+
!----------------------------
subroutine empty_tree(node)
 type(kdnode), intent(out) :: node(:)
 integer :: i

!$omp parallel do private(i)
 do i=1,size(node)
    node(i)%xcen = 0.
    node(i)%size = 0.
    node(i)%hmax = 0.
    node(i)%leftchild = 0
    node(i)%rightchild = 0
    node(i)%parent = 0
    node(i)%level  = 0
#ifdef GRAVITY
    node(i)%mass  = 0.
    node(i)%quads = 0.
    node(i)%octs  = 0.
#endif
 enddo
!$omp end parallel do

end subroutine empty_tree

!---------------------------------
!+
! routine to construct root node
!+
!---------------------------------
subroutine construct_root_node(np,nproot,irootnode,xmini,xmaxi,leaf_is_active,xyzh,xyzmh_ptmass,nptmass)
 use boundary, only:cross_boundary
 use mpidomain,only:isperiodic
 use part, only:iphase,iactive
 use part, only:isdead_or_accreted,ibelong
 use io,   only:fatal,id
 use dim,  only:ind_timesteps,mpi,periodic
 use part, only:isink,massoftype,igas,iamtype,maxphase,maxp,aprmassoftype,apr_level,ihsoft
 integer, intent(in)    :: np,irootnode
 integer, intent(out)   :: nproot
 real,    intent(out)   :: xmini(3), xmaxi(3)
 integer, intent(inout) :: leaf_is_active(:)
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(inout), optional :: xyzmh_ptmass(:,:)
 integer, intent(in),    optional :: nptmass
 integer :: i,ncross
 real    :: xminpart,yminpart,zminpart,xmaxpart,ymaxpart,zmaxpart
 real    :: xi, yi, zi

 xminpart = xyzh(1,1)
 yminpart = xyzh(2,1)
 zminpart = xyzh(3,1)
 xmaxpart = xminpart
 ymaxpart = yminpart
 zmaxpart = zminpart

 ncross = 0
 nproot = 0
 !$omp parallel default(none) &
 !$omp shared(np,xyzh,nptmass,xyzmh_ptmass) &
 !$omp shared(inodeparts,iphase,treecache,nproot) &
 !$omp shared(id,use_sinktree) &
 !$omp shared(isperiodic) &
 !$omp private(i,xi,yi,zi) &
 !$omp reduction(min:xminpart,yminpart,zminpart) &
 !$omp reduction(max:xmaxpart,ymaxpart,zmaxpart) &
 !$omp reduction(+:ncross)
 !$omp do schedule(guided,1)
 do i=1,np
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       if (periodic) call cross_boundary(isperiodic,xyzh(:,i),ncross)
       xi = xyzh(1,i)
       yi = xyzh(2,i)
       zi = xyzh(3,i)
       if (isnan(xi) .or. isnan(yi) .or. isnan(zi)) then
          call fatal('maketree','NaN in particle position, likely caused by NaN in force',i,var='x',val=xi)
       endif
       xminpart = min(xminpart,xi)
       yminpart = min(yminpart,yi)
       zminpart = min(zminpart,zi)
       xmaxpart = max(xmaxpart,xi)
       ymaxpart = max(ymaxpart,yi)
       zmaxpart = max(zmaxpart,zi)
    endif
 enddo
 !$omp enddo
 !$omp barrier
 if (use_sinktree) then
    if (nptmass>0) then
       !$omp do schedule(guided,1)
       do i=1,nptmass
          if (xyzmh_ptmass(4,i)>0.) then
             if (periodic) call cross_boundary(isperiodic,xyzmh_ptmass(1:3,i),ncross)
             xi = xyzmh_ptmass(1,i)
             yi = xyzmh_ptmass(2,i)
             zi = xyzmh_ptmass(3,i)
             if (isnan(xi) .or. isnan(yi) .or. isnan(zi)) then
                call fatal('maketree','NaN in ptmass position, likely caused by NaN in force',i,var='x',val=xi)
             endif
             xminpart = min(xminpart,xi)
             yminpart = min(yminpart,yi)
             zminpart = min(zminpart,zi)
             xmaxpart = max(xmaxpart,xi)
             ymaxpart = max(ymaxpart,yi)
             zmaxpart = max(zmaxpart,zi)
          endif
       enddo
       !$omp enddo
    endif
 endif
 !$omp end parallel

 do i=1,np
    isnotdead: if (.not.isdead_or_accreted(xyzh(4,i))) then
       nproot = nproot + 1

       if (ind_timesteps) then
          if (iactive(iphase(i))) then
             inodeparts(nproot) = i  ! +ve if active
          else
             inodeparts(nproot) = -i ! -ve if inactive
          endif
          if (use_apr) inodeparts(nproot) = abs(inodeparts(nproot))
       else
          inodeparts(nproot) = i
       endif
       treecache(1:4,nproot) = xyzh(1:4,i)
       if (maxphase==maxp) then
          if (use_apr) then
             treecache(5,nproot) = aprmassoftype(iamtype(iphase(i)),apr_level(i))
          else
             treecache(5,nproot) = massoftype(iamtype(iphase(i)))
          endif
       elseif (use_apr) then
          treecache(5,nproot) = aprmassoftype(igas,apr_level(i))
       else
          treecache(5,nproot) = massoftype(igas)
       endif
    endif isnotdead
 enddo

 if (use_sinktree) then
    if (nptmass > 0) then
       do i=1,nptmass
          if (mpi) then
             if (ibelong(maxpsph+i) /= id) cycle
          endif
          if (xyzmh_ptmass(4,i)<0.) cycle
          nproot = nproot + 1
          inodeparts(nproot) = (maxpsph) + i
          treecache(1:3,nproot) = xyzmh_ptmass(1:3,i)
          treecache(4,nproot)   = xyzmh_ptmass(ihsoft,i)
          treecache(5,nproot)   = xyzmh_ptmass(4,i)
       enddo
    endif
 endif

 if (nproot /= 0) then
    inoderange(1,irootnode) = 1
    inoderange(2,irootnode) = nproot
 else
    inoderange(:,irootnode) = 0
 endif

 xmini(1) = xminpart
 xmini(2) = yminpart
 xmini(3) = zminpart
 xmaxi(1) = xmaxpart
 xmaxi(2) = ymaxpart
 xmaxi(3) = zmaxpart


end subroutine construct_root_node

! also used for queue push
pure subroutine push_onto_stack(stackentry,node,parent,level,npnode,xmin,xmax)
 type(kdbuildstack), intent(out) :: stackentry
 integer,            intent(in)  :: node,parent,level
 integer,            intent(in)  :: npnode
 real,               intent(in)  :: xmin(3),xmax(3)

 stackentry%node   = node
 stackentry%parent = parent
 stackentry%level  = level
 stackentry%npnode = npnode
 stackentry%xmin   = xmin
 stackentry%xmax   = xmax

end subroutine push_onto_stack

! also used for queue pop
pure subroutine pop_off_stack(stackentry, istack, nnode, mymum, level, npnode, xmini, xmaxi)
 type(kdbuildstack), intent(in)    :: stackentry
 integer,            intent(inout) :: istack
 integer,            intent(out)   :: nnode, mymum, level, npnode
 real,               intent(out)   :: xmini(3), xmaxi(3)

 nnode  = stackentry%node
 mymum  = stackentry%parent
 level  = stackentry%level
 npnode = stackentry%npnode
 xmini  = stackentry%xmin
 xmaxi  = stackentry%xmax
 istack = istack - 1

end subroutine pop_off_stack

subroutine compute_nodes_cofm(npnode,nnode,xyzcofm,totmass_node,doparallel)
 use dim,       only:maxtypes
 use part,      only:massoftype,igas,npartoftype
 integer, intent(in)  :: npnode,nnode
 real,    intent(out) :: xyzcofm(3),totmass_node
 logical, intent(in)  :: doparallel
 real    :: pmassi,fac,dfac
 real    :: xi,yi,zi,xcofm,ycofm,zcofm
 integer :: i1,i
!
! to avoid round off error from repeated multiplication by pmassi (which is small)
! we compute the centre of mass with a factor relative to gas particles
! but only if gas particles are present
!
 pmassi = massoftype(igas)
 fac    = 1.
 totmass_node = 0.
 if (pmassi > 0.) then
    dfac = 1./pmassi
 else
    pmassi = massoftype(maxloc(npartoftype(2:maxtypes),1)+1)
    if (pmassi > 0.) then
       dfac = 1./pmassi
    else
       dfac = 1.
    endif
 endif

 ! note that dfac can be a constant value across all particles even if APR is used
 i1=inoderange(1,nnode)
 xcofm = 0.
 ycofm = 0.
 zcofm = 0.

 ! during initial queue build which is serial, we can parallelise this loop
 if (npnode > 1000 .and. doparallel) then
    !$omp parallel do schedule(static) default(none) &
    !$omp shared(npnode,dfac) &
    !$omp shared(treecache,i1) &
    !$omp private(i,xi,yi,zi) &
    !$omp firstprivate(pmassi,fac) &
    !$omp reduction(+:xcofm,ycofm,zcofm,totmass_node)
    do i=i1,i1+npnode-1
       xi = treecache(1,i)
       yi = treecache(2,i)
       zi = treecache(3,i)
       pmassi = treecache(5,i)
       fac    = pmassi*dfac ! to avoid round-off error
       totmass_node = totmass_node + pmassi
       xcofm = xcofm + fac*xi
       ycofm = ycofm + fac*yi
       zcofm = zcofm + fac*zi
    enddo
    !$omp end parallel do
 else
    do i=i1,i1+npnode-1
       xi = treecache(1,i)
       yi = treecache(2,i)
       zi = treecache(3,i)
       pmassi = treecache(5,i)
       fac    = pmassi*dfac ! to avoid round-off error
       totmass_node = totmass_node + pmassi
       xcofm = xcofm + fac*xi
       ycofm = ycofm + fac*yi
       zcofm = zcofm + fac*zi
    enddo
 endif

 xyzcofm = (/xcofm,ycofm,zcofm/)

 ! if there are no particles in this node, then the cofm will
 ! remain at zero
 if (totmass_node > 0.) then
    xyzcofm(:)   = xyzcofm(:)/(totmass_node*dfac)
 endif

end subroutine compute_nodes_cofm

subroutine set_nodes_properties(npnode,nnode,x0,totmass_node,mymum,nodeentry,xmini,xmaxi,&
                                level,global_build,doparallel,comp_node)
 use mpitree,   only:reduce_group
 use dim,       only:mpi
 type(kdnode),    intent(out)   :: nodeentry
 integer,         intent(in)    :: npnode,mymum,level,nnode
 real,            intent(inout) :: xmini(3), xmaxi(3), totmass_node
 real,            intent(in)    :: x0(3)
 logical,         intent(in)    :: doparallel,global_build,comp_node
 real    :: pmassi
 real    :: dx,dy,dz,dr2,xi,yi,zi,hi
 real    :: hmax,r2max,totmass
 integer :: i1,i
#ifdef GRAVITY
 real    :: quads(9)
 real    :: octs(10)
#endif

 pmassi  = 0.
 totmass = 0.
 r2max = 0.
 hmax  = 0.
#ifdef GRAVITY
 quads(:) = 0.
 octs(:) = 0.
#endif

 i1=inoderange(1,nnode)

 if (comp_node) then
    !--compute size of node
    ! parallelise this loop if node is large enough
    ! use !$omp parallel do when doparallel=.true. (not in parallel region)
    ! when doparallel=.false., we're already in a parallel region but can't use nested reductions
    ! so we'll use thread-local accumulators and combine at the end
    if (npnode > 1000 .and. doparallel) then
       !$omp parallel do schedule(static) default(none) &
       !$omp shared(npnode,treecache,x0,i1,use_geosplit) &
       !$omp private(i,xi,yi,zi,hi,dx,dy,dz,dr2) &
       !$omp firstprivate(pmassi) &
#ifdef GRAVITY
       !$omp reduction(+:totmass,quads,octs) &
#endif
       !$omp reduction(max:r2max,hmax)
       do i=i1,i1+npnode-1
          xi = treecache(1,i)
          yi = treecache(2,i)
          zi = treecache(3,i)
          hi = treecache(4,i)
          dx    = xi - x0(1)
          dy    = yi - x0(2)
          dz    = zi - x0(3)
          ! if (.not.use_geosplit) then
          dr2   = dx*dx + dy*dy + dz*dz
          r2max = max(r2max,dr2)
          ! endif
          hmax  = max(hmax,hi)
#ifdef GRAVITY
          pmassi = treecache(5,i)
          totmass  = totmass  + pmassi
          call add_node_moments(pmassi,dx,dy,dz,quads,octs)
#endif
       enddo
       !$omp end parallel do
    else
       do i=i1,i1+npnode-1
          xi = treecache(1,i)
          yi = treecache(2,i)
          zi = treecache(3,i)
          hi = treecache(4,i)
          dx    = xi - x0(1)
          dy    = yi - x0(2)
          dz    = zi - x0(3)
          ! if (.not.use_geosplit) then
          dr2   = dx*dx + dy*dy + dz*dz
          r2max = max(r2max,dr2)
          ! endif
          hmax = max(hmax,hi)
#ifdef GRAVITY
          pmassi = treecache(5,i)
          totmass  = totmass  + pmassi
          call add_node_moments(pmassi,dx,dy,dz,quads,octs)
#endif
       enddo
    endif
 endif

 ! if (use_geosplit) then
 !    r2max = 0.25*sum((xmaxi-xmini)**2)
 !    totmass_node  = totmass
 ! endif
 ! reduce node limits and quads across MPI tasks belonging to this group
 if (mpi .and. global_build) then
    r2max     = reduce_group(r2max,'max',level)
    hmax      = reduce_group(hmax,'max',level)

    xmini(1)  = reduce_group(xmini(1),'min',level)
    xmini(2)  = reduce_group(xmini(2),'min',level)
    xmini(3)  = reduce_group(xmini(3),'min',level)

    xmaxi(1)  = reduce_group(xmaxi(1),'max',level)
    xmaxi(2)  = reduce_group(xmaxi(2),'max',level)
    xmaxi(3)  = reduce_group(xmaxi(3),'max',level)
#ifdef GRAVITY
    totmass_node = reduce_group(totmass_node, "+", level)
    quads(1)  = reduce_group(quads(1),'+',level)
    quads(2)  = reduce_group(quads(2),'+',level)
    quads(3)  = reduce_group(quads(3),'+',level)
    quads(4)  = reduce_group(quads(4),'+',level)
    quads(5)  = reduce_group(quads(5),'+',level)
    quads(6)  = reduce_group(quads(6),'+',level)
    quads(7)  = reduce_group(quads(7),'+',level)
    quads(8)  = reduce_group(quads(8),'+',level)
    quads(9)  = reduce_group(quads(9),'+',level)
    octs(1)   = reduce_group(octs(1),'+',level)
    octs(2)   = reduce_group(octs(2),'+',level)
    octs(3)   = reduce_group(octs(3),'+',level)
    octs(4)   = reduce_group(octs(4),'+',level)
    octs(5)   = reduce_group(octs(5),'+',level)
    octs(6)   = reduce_group(octs(6),'+',level)
    octs(7)   = reduce_group(octs(7),'+',level)
    octs(8)   = reduce_group(octs(8),'+',level)
    octs(9)   = reduce_group(octs(9),'+',level)
    octs(10)  = reduce_group(octs(10),'+',level)
#endif
 endif

 ! assign properties to node
 nodeentry%xcen       = x0(:)
 nodeentry%size       = sqrt(r2max) + epsilon(r2max)
 nodeentry%hmax       = hmax
 nodeentry%parent     = mymum
 nodeentry%level      = level
#ifdef GRAVITY
 nodeentry%mass       = totmass_node
 nodeentry%quads      = quads
 nodeentry%octs       = octs
 nodeentry%tobecached = 1
 nodeentry%fcached    = .false.
 nodeentry%ncached    = .false.
#endif

end subroutine set_nodes_properties


!--------------------------------------------------------------------
!+
!  create all the properties for a given node such as centre of mass,
!  size, max smoothing length, etc
!  will also split the node if necessary, setting wassplit=true
!  returns the left and right child information if split
!+
!--------------------------------------------------------------------
subroutine construct_node(nodeentry, nnode, mymum, level, xmini, xmaxi, npnode, doparallel,&
                          il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr,ncells, leaf_is_active, &
                          minlevel, maxlevel, wassplit, global_build,apr_tree, &
                          xyzmh_ptmass)
 use dim,       only:maxtypes,mpi,ind_timesteps
 use io,        only:fatal,error
 use mpitree,   only:get_group_cofm,reduce_group
 type(kdnode),    intent(out)   :: nodeentry
 integer,         intent(in)    :: nnode, mymum, level
 real,            intent(inout) :: xmini(3), xmaxi(3)
 integer,         intent(in)    :: npnode
 logical,         intent(in)    :: doparallel
 integer,         intent(out)   :: il, ir, nl, nr
 real,            intent(out)   :: xminl(3), xmaxl(3), xminr(3), xmaxr(3)
 integer(kind=8), intent(inout) :: ncells
 integer,         intent(out)   :: leaf_is_active(:)
 integer,         intent(inout) :: maxlevel, minlevel
 logical,         intent(out)   :: wassplit
 logical,         intent(in)    :: global_build
 logical,         intent(in)    :: apr_tree
 real,            intent(in), optional :: xyzmh_ptmass(:,:)

 integer(kind=8) :: myslot
 real    :: xyzcofm(3)
 real    :: totmass_node
 real    :: xyzcofmg(3)
 real    :: totmassg
 integer :: npnodetot
 logical :: nodeisactive,comp_node
 integer :: i,npcounter,ipart
 real    :: x0(3)
 integer :: iaxis
 real    :: xpivot

 nodeisactive = .false.
 if (inoderange(1,nnode) > 0) then
    checkactive: do i = inoderange(1,nnode),inoderange(2,nnode)
       if (inodeparts(i) > 0) then
          nodeisactive = .true.
          exit checkactive
       endif
    enddo checkactive
    npcounter = inoderange(2,nnode) - inoderange(1,nnode) + 1
 else
    npcounter = 0
 endif

 if (npcounter /= npnode) then
    print*,'constructing node ',nnode,': found ',npcounter,' particles, expected:',npnode,' particles for this node'
    call fatal('maketree', 'expected number of particles in node differed from actual number')
 endif

 if (mpi .and. global_build) then
    npnodetot = reduce_group(npnode,'+',level)
 else
    npnodetot = npnode
 endif

 ! following lines to avoid compiler warnings on intent(out) variables
 ir = 0
 il = 0
 nl = 0
 nr = 0
 wassplit    = (npnodetot > minpart)
 if ((.not. global_build) .and. (npnode  <  1)) return ! node has no particles, just quit

 xyzcofm(:) = 0.

 call compute_nodes_cofm(npnode,nnode,xyzcofm,totmass_node,doparallel)
 ! if this is global node construction, get the cofm and total mass
 ! of all particles in this node (some on other MPI tasks)
 if (mpi .and. global_build) then
    call get_group_cofm(xyzcofm,totmass_node,level,xyzcofmg,totmassg)
    xyzcofm = xyzcofmg
    totmass_node = totmassg
 endif
 ! checks the reduced mass in the case of global maketree
 if (totmass_node<=0. .and. use_apr) call fatal('mtree + apr', &
    'totmass_node==0, something almost certainly wrong with aprmassoftype')
 if (totmass_node<=0.) call fatal('mtree','totmass_node==0',val=totmass_node)

 if (use_geosplit) then !--for geotree we use the middle point to split the node and propagate properties after
    x0        = (xmaxi+xmini)*0.5       ! middle point of the node
    comp_node = .true. !.not.wassplit
 else  !--for gravity and default KDtree, we need the centre of the node to be the centre of mass
    x0 = xyzcofm
    comp_node = .true.
 endif


 call set_nodes_properties(npnode,nnode,xyzcofm,totmass_node,mymum,nodeentry,xmini,xmaxi,&
                           level,global_build,doparallel,comp_node)

 if (apr_tree)   wassplit = (npnode > 2)

 if (.not. wassplit) then
    nodeentry%leftchild  = 0
    nodeentry%rightchild = 0
    maxlevel = max(level,maxlevel)
    minlevel = min(level,minlevel)
    ! individual timesteps where we mark leaf node as active/inactive
    if (ind_timesteps) then
       !
       !--mark leaf node as active (contains some active particles)
       !  or inactive by setting the firstincell to +ve (active) or -ve (inactive)
       !
       if (nodeisactive) then
          leaf_is_active(nnode) = 1
       else
          leaf_is_active(nnode) = -1
       endif
    else
       leaf_is_active(nnode) = 1
    endif
 else ! split this node and add children to stack
    iaxis  = maxloc(xmaxi - xmini,1) ! split along longest axis
    xpivot = x0(iaxis)               ! split middle longest axis

    if (maxlevel > maxdepth) call fatal('maketree','maximum tree depth reached !!')
    ! create two children nodes and point to them from current node
    ! always use G&R indexing for global tree
    if (((level < maxlevel_indexed) .or. global_build) .and. (.not. use_geosplit)) then
       il = 2*nnode   ! indexing as per Gafton & Rosswog (2011)
       ir = il + 1
    else
       ! no need to lock, we could just atomic the update
       !$omp atomic capture
       ncells = ncells + 2
       myslot = ncells
       !$omp end atomic
       ir = int(myslot)
       il = int(myslot-1)
       if (ir > ncellsmax) call fatal('maketree',&
          'number of nodes exceeds array dimensions, increase ncellsmax and recompile',ival=int(ncellsmax))
    endif
    nodeentry%leftchild  = il
    nodeentry%rightchild = ir

    leaf_is_active(nnode) = 0

    if (npnode > 0) then
       if (apr_tree) then
          ! apr special sort - only used for merging particles
          call special_sort_particles_in_cell(iaxis,inoderange(1,nnode),inoderange(2,nnode),inoderange(1,il),inoderange(2,il),&
                                    inoderange(1,ir),inoderange(2,ir),nl,nr,xpivot,treecache,inodeparts,&
                                    npnode)
       else
          ! regular sort
          call sort_particles_in_cell(iaxis,inoderange(1,nnode),inoderange(2,nnode),inoderange(1,il),inoderange(2,il),&
                                  inoderange(1,ir),inoderange(2,ir),nl,nr,xpivot,treecache,inodeparts)
       endif

       if (nr + nl  /=  npnode) then
          call error('maketree','number of left + right != parent while splitting (likely cause: NaNs in position arrays)')
       endif

       ! see if all the particles ended up in one node, if so, arbitrarily build 2 cells. This should never happen
       if ( (.not. global_build) .and. ((nl==npnode) .or. (nr==npnode)) ) then
          ! no need to move particles because if they all ended up in one node,
          ! then they are still in the original order
          nl = npnode / 2
          inoderange(1,il) = inoderange(1,nnode)
          inoderange(2,il) = inoderange(1,nnode) + nl - 1
          inoderange(1,ir) = inoderange(1,nnode) + nl
          inoderange(2,ir) = inoderange(2,nnode)
          nr = npnode - nl
       endif

       ! compute min/max with explicit loops for better cache behavior
       xminl(1) = treecache(1,inoderange(1,il))
       xminl(2) = treecache(2,inoderange(1,il))
       xminl(3) = treecache(3,inoderange(1,il))
       xmaxl(1) = xminl(1)
       xmaxl(2) = xminl(2)
       xmaxl(3) = xminl(3)
       do ipart=inoderange(1,il)+1,inoderange(2,il)
          xminl(1) = min(xminl(1),treecache(1,ipart))
          xminl(2) = min(xminl(2),treecache(2,ipart))
          xminl(3) = min(xminl(3),treecache(3,ipart))
          xmaxl(1) = max(xmaxl(1),treecache(1,ipart))
          xmaxl(2) = max(xmaxl(2),treecache(2,ipart))
          xmaxl(3) = max(xmaxl(3),treecache(3,ipart))
       enddo

       xminr(1) = treecache(1,inoderange(1,ir))
       xminr(2) = treecache(2,inoderange(1,ir))
       xminr(3) = treecache(3,inoderange(1,ir))
       xmaxr(1) = xminr(1)
       xmaxr(2) = xminr(2)
       xmaxr(3) = xminr(3)
       do ipart=inoderange(1,ir)+1,inoderange(2,ir)
          xminr(1) = min(xminr(1),treecache(1,ipart))
          xminr(2) = min(xminr(2),treecache(2,ipart))
          xminr(3) = min(xminr(3),treecache(3,ipart))
          xmaxr(1) = max(xmaxr(1),treecache(1,ipart))
          xmaxr(2) = max(xmaxr(2),treecache(2,ipart))
          xmaxr(3) = max(xmaxr(3),treecache(3,ipart))
       enddo
    else
       nl = 0
       nr = 0
       xminl = 0.0
       xmaxl = 0.0
       xminr = 0.0
       xmaxr = 0.0
    endif

    ! Reduce node limits of children across MPI tasks belonging to this group.
    ! The synchronisation needs to happen here, not at the next level, because
    ! the groups will be independent by then.
    if (mpi .and. global_build) then
       xminl(1) = reduce_group(xminl(1),'min',level)
       xminl(2) = reduce_group(xminl(2),'min',level)
       xminl(3) = reduce_group(xminl(3),'min',level)

       xmaxl(1) = reduce_group(xmaxl(1),'max',level)
       xmaxl(2) = reduce_group(xmaxl(2),'max',level)
       xmaxl(3) = reduce_group(xmaxl(3),'max',level)

       xminr(1) = reduce_group(xminr(1),'min',level)
       xminr(2) = reduce_group(xminr(2),'min',level)
       xminr(3) = reduce_group(xminr(3),'min',level)

       xmaxr(1) = reduce_group(xmaxr(1),'max',level)
       xmaxr(2) = reduce_group(xmaxr(2),'max',level)
       xmaxr(3) = reduce_group(xmaxr(3),'max',level)
    endif

 endif

end subroutine construct_node

!----------------------------------------------------------------
!+
!  Categorise particles into daughter nodes by whether they
!  fall to the left or the right of the pivot axis
!+
!----------------------------------------------------------------
subroutine sort_particles_in_cell(iaxis,imin,imax,min_l,max_l,min_r,max_r,nl,nr,xpivot,&
                                   treecache,inodeparts)
 integer, intent(in)    :: iaxis,imin,imax
 integer, intent(out)   :: min_l,max_l,min_r,max_r,nl,nr
 real,    intent(inout) :: xpivot,treecache(:,:)
 integer, intent(inout) :: inodeparts(:)
 logical :: i_lt_pivot,j_lt_pivot
 integer :: inodeparts_swap,i,j
 real :: xyzh_swap(5)
 real :: xi_coord, xj_coord

 !print*,'nnode ',imin,imax,' pivot = ',iaxis,xpivot
 i = imin
 j = imax

 xi_coord = treecache(iaxis,i)
 xj_coord = treecache(iaxis,j)
 i_lt_pivot = xi_coord <= xpivot
 j_lt_pivot = xj_coord <= xpivot
 !  k = 0

 do while(i < j)
    if (i_lt_pivot) then
       i = i + 1
       xi_coord = treecache(iaxis,i)
       i_lt_pivot = xi_coord <= xpivot
    else
       if (.not.j_lt_pivot) then
          j = j - 1
          xj_coord = treecache(iaxis,j)
          j_lt_pivot = xj_coord <= xpivot
       else
          ! swap i and j positions in list
          inodeparts_swap = inodeparts(i)
          xyzh_swap(1:5)  = treecache(1:5,i)

          inodeparts(i)   = inodeparts(j)
          treecache(1:5,i) = treecache(1:5,j)

          inodeparts(j)   = inodeparts_swap
          treecache(1:5,j) = xyzh_swap(1:5)

          i = i + 1
          j = j - 1
          xi_coord = treecache(iaxis,i)
          xj_coord = treecache(iaxis,j)
          i_lt_pivot = xi_coord <= xpivot
          j_lt_pivot = xj_coord <= xpivot
          ! k = k + 1
       endif
    endif
 enddo
 if (.not.i_lt_pivot) i = i - 1
 if (j_lt_pivot)      j = j + 1

 min_l = imin
 max_l = i
 min_r = j
 max_r = imax

 if ( j /= i+1) print*,' ERROR ',i,j
 nl = max_l - min_l + 1
 nr = max_r - min_r + 1

end subroutine sort_particles_in_cell

!----------------------------------------------------------------
!+
!  Categorise particles into daughter nodes by whether they
!  fall to the left or the right of the pivot axis, but additionally
!  force the cells to have a certain minimum number of particles per cell
!+
!----------------------------------------------------------------
subroutine special_sort_particles_in_cell(iaxis,imin,imax,min_l,max_l,min_r,max_r,&
                                nl,nr,xpivot,treecache,inodeparts,npnode)
 use io, only:error
 integer, intent(in)    :: iaxis,imin,imax,npnode
 integer, intent(out)   :: min_l,max_l,min_r,max_r,nl,nr
 real,    intent(inout) :: xpivot,treecache(:,:)
 integer, intent(inout) :: inodeparts(:)
 logical :: i_lt_pivot,j_lt_pivot,slide_l,slide_r
 integer :: inodeparts_swap,i,j,nchild_in
 integer :: k,ii,rem_nr,rem_nl
 real :: xyzh_swap(5),dpivot(npnode)

 dpivot = 0.0
 nchild_in = 2

 if (modulo(npnode,nchild_in) > 0) then
    call error('apr sort','number of particles sent in to kdtree is not divisible by 2')
 endif

! print*,'nnode ',imin,imax,npnode,' pivot = ',iaxis,xpivot
 i = imin
 j = imax

 i_lt_pivot = treecache(iaxis,i) <= xpivot
 j_lt_pivot = treecache(iaxis,j) <= xpivot
 dpivot(i-imin+1) = xpivot - treecache(iaxis,i)
 dpivot(j-imin+1) = xpivot - treecache(iaxis,j)
 !k = 0
 do while(i < j)
    if (i_lt_pivot) then
       i = i + 1
       dpivot(i-imin+1) = xpivot - treecache(iaxis,i)
       i_lt_pivot = treecache(iaxis,i) <= xpivot
    else
       if (.not.j_lt_pivot) then
          j = j - 1
          dpivot(j-imin+1) = xpivot - treecache(iaxis,j)
          j_lt_pivot = treecache(iaxis,j) <= xpivot
       else
          ! swap i and j positions in list
          inodeparts_swap = inodeparts(i)
          xyzh_swap(1:5)  = treecache(1:5,i)

          inodeparts(i)   = inodeparts(j)
          treecache(1:5,i) = treecache(1:5,j)

          inodeparts(j)   = inodeparts_swap
          treecache(1:5,j) = xyzh_swap(1:5)

          i = i + 1
          j = j - 1

          dpivot(i-imin+1) = xpivot - treecache(iaxis,i)
          dpivot(j-imin+1) = xpivot - treecache(iaxis,j)

          i_lt_pivot = treecache(iaxis,i) <= xpivot
          j_lt_pivot = treecache(iaxis,j) <= xpivot
       endif
    endif
 enddo

 if (.not.i_lt_pivot) then
    i = i - 1
    dpivot(i-imin+1) = xpivot - treecache(iaxis,i)
 endif
 if (j_lt_pivot) then
    j = j + 1
    dpivot(j-imin+1) = xpivot - treecache(iaxis,j)
 endif

 min_l = imin
 max_l = i
 min_r = j
 max_r = imax

 if ( j /= i+1) print*,' ERROR ',i,j
 nl = max_l - min_l + 1
 nr = max_r - min_r + 1

 ! does the pivot need to be adjusted?
 rem_nl = modulo(nl,nchild_in)
 rem_nr = modulo(nr,nchild_in)
 if (rem_nl == 0 .and. rem_nr == 0) return

 ! Decide which direction the pivot needs to go
 if (rem_nl < rem_nr) then
    slide_l = .true.
    slide_r = .false.
 else
    slide_l = .false.
    slide_r = .true.
 endif
 ! Override this if there's less than nchild*2 in the cell
 if (nl < nchild_in) then
    slide_r = .true.
    slide_l = .false.
 elseif (nr < nchild_in) then
    slide_r = .false.
    slide_l = .true.
 endif

 ! Move across particles by distance from xpivot till we get
 ! the right number of particles in each cell
 if (slide_r) then
    do ii = 1,rem_nr
       ! next particle to shift across
       k = minloc(dpivot,dim=1,mask=dpivot > 0.) + imin - 1
       if (k-imin+1==0) k = maxloc(dpivot,dim=1,mask=dpivot < 0.) + imin - 1

       ! swap this with the first particle on the j side
       inodeparts_swap = inodeparts(k)
       xyzh_swap(1:5)  = treecache(1:5,k)

       inodeparts(k)   = inodeparts(j)
       treecache(1:5,k) = treecache(1:5,j)

       inodeparts(j)   = inodeparts_swap
       treecache(1:5,j) = xyzh_swap(1:5)

       ! and now shift to the right
       i = i + 1
       j = j + 1

       ! ditch it, go again
       dpivot(k-imin+1) = huge(k-imin+1)
    enddo
 else
    do ii = 1,rem_nl
       ! next particle to shift across
       k = maxloc(dpivot,dim=1,mask=dpivot < 0.) + imin - 1
       if (k-imin+1==0) k = minloc(dpivot,dim=1,mask=dpivot > 0.) + imin - 1

       ! swap this with the last particle on the i side
       inodeparts_swap = inodeparts(k)
       xyzh_swap(1:5)  = treecache(1:5,k)

       inodeparts(k)   = inodeparts(i)
       treecache(1:5,k) = treecache(1:5,i)

       inodeparts(i)   = inodeparts_swap
       treecache(1:5,i) = xyzh_swap(1:5)

       ! and now shift to the left
       i = i - 1
       j = j - 1

       ! ditch it, go again
       dpivot(k-imin+1) = huge(k-imin+1)

    enddo
 endif

 ! tidy up outputs
 max_l = i
 min_r = j
 nl = max_l - min_l + 1
 nr = max_r - min_r + 1

end subroutine special_sort_particles_in_cell

subroutine propagate_upward(ncells,node)
 use io, only: fatal
!$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
 integer,      intent(in)    :: ncells
 type(kdnode), intent(inout) :: node(:)
 integer, allocatable :: levcount(:), levstart(:), nodelist(:)
 integer, allocatable :: levcount_t(:,:)
 integer :: i,lvl,id,istart,iend,il,ir,npnode,nthreads,it,tid
 integer :: rcnt(0:maxlevel)
 real    :: mnode

 !--
 !  sort internal nodes (both children present) by tree level using a
 !  parallel counting sort. The count pass keeps a histogram per thread
 !  (levcount_t); prefix sums give the start of each level in nodelist;
 !  for the scatter each thread uses a private running offset (rcnt),
 !  initialised from the cumulative histogram of the threads that precede
 !  it, so that each thread writes a disjoint slice of nodelist, race-free.
 !  Count and scatter use the same schedule(static), which guarantees that
 !  thread tid counts exactly the nodes it later scatters.
 !--
 nthreads = 1
!$ nthreads = omp_get_max_threads()
 allocate(levcount(0:maxlevel), levstart(0:maxlevel+1), nodelist(ncells))
 allocate(levcount_t(0:maxlevel,1:nthreads))
 levcount_t = 0
 levcount = 0

 !$omp parallel default(none) &
 !$omp shared(node,ncells,maxlevel,nthreads,nodelist,levcount,levstart,levcount_t,inoderange) &
 !$omp private(i,lvl,id,istart,iend,il,ir,mnode,npnode) &
 !$omp private(it,tid,rcnt)
 tid = omp_get_thread_num()

 !$omp do schedule(static)
 do i = 1, ncells
    if (node(i)%leftchild > 0) then
       lvl = node(i)%level
       if (lvl < 0 .or. lvl > maxlevel) cycle
       levcount_t(lvl,tid+1) = levcount_t(lvl,tid+1) + 1
    endif
 enddo
 !$omp end do

 !--combine the per-thread histograms and build the per-level starts
 !$omp single
 do it=1,nthreads
    levcount(:) = levcount(:) + levcount_t(:,it)
 enddo
 levstart(0) = 1
 do lvl = 1, maxlevel+1
    levstart(lvl) = levstart(lvl-1) + levcount(lvl-1)
 enddo
 !$omp end single

 !--scatter: start each thread's offsets after the counts of the preceding
 !  threads, then fill nodelist within each thread's own level slices
 rcnt(:) = levstart(0:maxlevel)
 do it=1,tid
    rcnt(:) = rcnt(:) + levcount_t(:,it)
 enddo
 !$omp do schedule(static)
 do i = 1, ncells
    if (node(i)%leftchild > 0) then
       lvl = node(i)%level
       if (lvl < 0 .or. lvl > maxlevel) cycle
       nodelist(rcnt(lvl)) = i
       rcnt(lvl) = rcnt(lvl) + 1
    endif
 enddo
 !$omp end do

 !--
 !  propagate properties upward, one level at a time (deepest first).
 !  nodes in the same level are not ancestor/descendant of each other, so
 !  each level can be processed in parallel, in place: a node only reads
 !  its children (already final from the previous level) and writes itself.
 !--
 do lvl = maxlevel, 0, -1

    istart = levstart(lvl)
    iend = levstart(lvl) + levcount(lvl) - 1
    if (iend < istart) cycle

    !$omp do schedule(runtime)
    do id = istart, iend
       i  = nodelist(id)
       il = node(i)%leftchild
       ir = node(i)%rightchild
       call translate_node(node,i,il,ir)

       npnode = inoderange(2,i) - inoderange(1,i) + 1
       mnode  = node(i)%mass
       if (npnode > 1 .and. mnode < epsilon(mnode)) then
          call fatal('mtree','mnode==0',val=mnode)
       endif
    enddo
    !$omp end do

 enddo
 !$omp end parallel

 deallocate(levcount,levstart,nodelist,levcount_t)

end subroutine propagate_upward

subroutine translate_node(node,ip,il,ir)
 type(kdnode), intent(inout) :: node(:)
 integer,      intent(in)    :: ip,il,ir
 real    :: dx(3),massp,massc,quadsp(9),quadsc(9),dips(3),hmaxc,hmaxp
 integer :: j,ic(2),k

 ic = (/il,ir/)
 massp  = 0.
 hmaxp  = 0.
 quadsp = 0.
 dips   = 0.

 do j=1,2
    k = ic(j)
    massc  = node(k)%mass
    quadsc = node(k)%quads
    hmaxc  = node(k)%hmax
    dx     = node(k)%xcen - node(ip)%xcen

    massp       = massp + massc
    dips(1)     = quadsc(1) + dx(1)*massc
    dips(2)     = quadsc(2) + dx(2)*massc
    dips(3)     = quadsc(3) + dx(3)*massc
    quadsp(1)   = quadsp(1) + dips(1)
    quadsp(2)   = quadsp(2) + dips(2)
    quadsp(3)   = quadsp(3) + dips(3)
    quadsp(4)   = quadsp(4) + quadsc(4) + dx(1)*dips(1) + quadsc(1)*dx(1)
    quadsp(5)   = quadsp(5) + quadsc(5) + dx(1)*dips(2) + quadsc(1)*dx(2)
    quadsp(6)   = quadsp(6) + quadsc(6) + dx(1)*dips(3) + quadsc(1)*dx(3)
    quadsp(7)   = quadsp(7) + quadsc(7) + dx(2)*dips(2) + quadsc(2)*dx(2)
    quadsp(8)   = quadsp(8) + quadsc(8) + dx(2)*dips(3) + quadsc(2)*dx(3)
    quadsp(9)   = quadsp(9) + quadsc(9) + dx(3)*dips(3) + quadsc(3)*dx(3)

    hmaxp = max(hmaxp,hmaxc)

 enddo

 node(ip)%mass  = massp
 node(ip)%quads = quadsp
 node(ip)%hmax  = hmaxp

end subroutine translate_node


!----------------------------------------------------------------
!+
!  Cache particles within identified neighbour nodes
!+
!----------------------------------------------------------------
subroutine cache_neighbours(nneigh,isrc,ixyzcachesize,maxcache,listneigh,xyzcache,xoffset,yoffset,zoffset)
 use part, only:rho,gradh
 use dim,  only:igradomega,igradzeta
#ifdef GRAVITY
 use dim,  only:igradsoft
#endif
 integer, intent(in)    :: isrc,ixyzcachesize,maxcache
 real,    intent(in)    :: xoffset,yoffset,zoffset
 integer, intent(inout) :: nneigh
 integer, intent(out)   :: listneigh(:)
 real,    intent(out)   :: xyzcache(:,:)
 integer :: npnode,ipart,num_to_cache,ip,inode

 npnode = inoderange(2,isrc) - inoderange(1,isrc) + 1

 if (nneigh + npnode <= ixyzcachesize) then
    num_to_cache = npnode
 elseif (nneigh < ixyzcachesize) then
    num_to_cache = ixyzcachesize - nneigh
 else
    num_to_cache = 0
 endif

 if (num_to_cache > 0) then
    do ipart=1,num_to_cache
       inode = inoderange(1,isrc)+ipart-1
       ip = abs(inodeparts(inode))
       listneigh(nneigh+ipart)  = ip
       xyzcache(1,nneigh+ipart) = treecache(1,inode) + xoffset
       xyzcache(2,nneigh+ipart) = treecache(2,inode) + yoffset
       xyzcache(3,nneigh+ipart) = treecache(3,inode) + zoffset
       if (maxcache >= 4) then
          xyzcache(ih1,nneigh+ipart) = 1./treecache(4,inode)
       endif
       if (maxcache >= 5) then
          xyzcache(im,nneigh+ipart) = treecache(5,inode)
       endif
       if (maxcache >= 7) then
          if (ip <= maxpsph) then
             xyzcache(irho,nneigh+ipart)        = rho(ip)
             xyzcache(izetaomega,nneigh+ipart) = real(gradh(igradzeta,ip))*real(gradh(igradomega,ip))
          else
             xyzcache(irho,nneigh+ipart)        = 0.
             xyzcache(izetaomega,nneigh+ipart) = 0.
          endif
       endif
#ifdef GRAVITY
       if (maxcache >= 8) then
          if (ip <= maxpsph) then
             xyzcache(isoftomega,nneigh+ipart) = real(gradh(igradsoft,ip))*real(gradh(igradomega,ip))
          else
             xyzcache(isoftomega,nneigh+ipart) = 0.
          endif
       endif
#endif
    enddo
 endif

 if (num_to_cache < npnode) then
    do ipart=num_to_cache+1,npnode
       listneigh(nneigh+ipart) = abs(inodeparts(inoderange(1,isrc)+ipart-1))
    enddo
 endif

 nneigh = nneigh + npnode

end subroutine cache_neighbours

!----------------------------------------------------------------
!+
!  Routine to walk tree for neighbour search
!  (all particles within a given h_i and optionally within h_j)
!+
!----------------------------------------------------------------
subroutine getneigh(node,xpos,xsizei,rcuti,listneigh,nneigh,xyzcache,ixyzcachesize,leaf_is_active,&
                    get_hj,get_f,fnode,remote_export,nq)
 use io,       only:fatal,id
 use part,     only:gravity
 use kernel,   only:radkern
 type(kdnode), intent(in)  :: node(:) !ncellsmax+1)
 integer,      intent(in)  :: ixyzcachesize
 real,         intent(in)  :: xpos(3)
 real,         intent(in)  :: xsizei,rcuti
 integer,      intent(out) :: listneigh(:)
 integer,      intent(out) :: nneigh
 real,         intent(out) :: xyzcache(:,:)
 integer,      intent(in)  :: leaf_is_active(:)
 logical,      intent(in)  :: get_hj
 logical,      intent(in)  :: get_f
 real,         intent(out), optional :: fnode(lenfgrav)
 logical,      intent(out), optional :: remote_export(:)
 integer,      intent(in),  optional :: nq
 integer :: maxcache
 integer :: n,istack,il,ir
 integer :: nstack(maxdepth)
 real :: dx,dy,dz,xsizej,rcutj
 real :: rcut,rcut2,r2
 real :: xoffset,yoffset,zoffset,tree_acc2
 logical :: open_tree_node
 logical :: global_walk
#ifdef GRAVITY
 real :: quads(9)
 real :: dr,totmass_node
#endif
 tree_acc2 = tree_accuracy*tree_accuracy
 if (get_f .and. .not.present(fnode)) then
    call fatal('getneigh','get_f but fnode not passed...')
 endif
 if (present(fnode)) fnode(:) = 0.
 rcut     = rcuti

 if (ixyzcachesize > 0) then
    maxcache = size(xyzcache,1)
 else
    maxcache = 0
 endif

 if (present(remote_export)) then
    remote_export = .false.
    global_walk = .true.
 else
    global_walk = .false.
 endif

 nneigh = 0
 istack = 1
 nstack(istack) = irootnode
 open_tree_node = .false.

 over_stack: do while(istack /= 0)
    n = nstack(istack)
    istack = istack - 1
    call get_sep(xpos,node(n)%xcen,dx,dy,dz,xoffset,yoffset,zoffset,r2)
    xsizej  = node(n)%size
    il      = node(n)%leftchild
    ir      = node(n)%rightchild
#ifdef GRAVITY
    totmass_node = node(n)%mass
    quads        = node(n)%quads
#endif

    if (get_hj) then  ! find neighbours within both hi and hj
       rcutj = radkern*node(n)%hmax
       rcut  = max(rcuti,rcutj)
    endif
    rcut2 = (xsizei + xsizej + rcut)**2   ! node size + search radius
    if (gravity) open_tree_node = tree_acc2*r2 < (xsizei + xsizej)**2   ! tree opening criterion for self-gravity
    if_open_node: if ((r2 < rcut2) .or. open_tree_node) then
       if_leaf: if (leaf_is_active(n) /= 0) then ! once we hit a leaf node, retrieve contents into trial neighbour cache
          if_global_walk: if (global_walk) then
             ! id is stored in cellatid (passed through into leaf_is_active) as id + 1
             if (leaf_is_active(n) /= (id + 1)) then
                remote_export(leaf_is_active(n)) = .true.
             endif
          else
             call cache_neighbours(nneigh,n,ixyzcachesize,maxcache,listneigh,xyzcache,xoffset,yoffset,zoffset)
          endif if_global_walk
       else
          if (istack+2 > ncellsmax+1) call fatal('getneigh','stack overflow in getneigh')
          if (il /= 0) then
             istack = istack + 1
             nstack(istack) = il
          endif
          if (ir /= 0) then
             istack = istack + 1
             nstack(istack) = ir
          endif
       endif if_leaf
#ifdef GRAVITY
    elseif (get_f) then ! if_open_node
       ! When searching for neighbours of this node, the tree walk may encounter
       ! nodes on the global tree that it does not need to open, so it should
       ! just add the contribution to fnode. However, when walking a different
       ! part of the tree, it may then become necessary to export this node to
       ! a remote task. When it arrives at the remote task, it will then walk
       ! the remote tree.
       !
       ! The complication arises when tree refinment is enabled, which puts part
       ! of the remote tree onto the global tree. fnode will be double counted
       ! if a contribution is made on the global tree and a separate branch
       ! causes it to be sent to a remote task, where that contribution is
       ! counted again.
       !
       ! The solution is to not count the parts of the local tree that have been
       ! added onto the global tree.

       count_gravity: if ( global_walk .or. (n > irefine) ) then
          !
          !--long range force on node due to distant node, along node centres
          !  along with derivatives in order to perform series expansion
          !
          dr = 1./sqrt(r2)
          call compute_M2L(dx,dy,dz,dr,totmass_node,quads,fnode)

       endif count_gravity
#endif

    endif if_open_node
 enddo over_stack

end subroutine getneigh

!----------------------------------------------------------------
!+
!  Routine to walk tree for neighbour search (SFMM version)
!  (all particles within a given h_i and optionally within h_j)
!  A dual tree walk is used to compute
!  every node-node interactions
!+
!----------------------------------------------------------------
subroutine getneigh_dual(node,xpos,xsizei,rcuti,listneigh,nneigh,xyzcache,ixyzcachesize,leaf_is_active,&
                              get_hj,get_f,fnode,icell)
 type(kdnode), intent(inout) :: node(:) !ncellsmax+1)
 integer,      intent(in)    :: ixyzcachesize
 real,         intent(in)    :: xpos(3)
 real,         intent(in)    :: xsizei,rcuti
 integer,      intent(out)   :: listneigh(:)
 integer,      intent(out)   :: nneigh
 real,         intent(out)   :: xyzcache(:,:)
 integer,      intent(in)    :: leaf_is_active(:)
 logical,      intent(in)    :: get_hj
 logical,      intent(in)    :: get_f
 real,         intent(out)   :: fnode(lenfgrav)
 integer,      intent(in)    :: icell
 integer :: istack,i,iparent,idstbranch,idst,isrc,maxcache,tobecached,ibase
 integer :: branch(maxdepth),nparents,stack(3,2048),startwith(2)
 real    :: dx,dy,dz,xoffset,yoffset,zoffset
 real    :: tree_acc2
 real    :: fnode_acc(lenfgrav)
 logical :: stackit,fcached

 tree_acc2 = tree_accuracy*tree_accuracy

 if (ixyzcachesize > 0) then
    maxcache = size(xyzcache,1)
 else
    maxcache = 0
 endif

 call get_list_of_parent_nodes(icell,node,branch,nparents,startwith)

 neighnodecount_branch = 0
 neighnode_branch = 0
 fnode_branch = 0.
 fnode_acc    = 0.
 nneigh = 0
 istack = 0
 xoffset = 0.
 yoffset = 0.
 zoffset = 0.

 if (use_cache .and. startwith(2) > 0) then
    ! print*, real(nparents-startwith(2)+1)/nparents,nparents,neighnodecache_count(startwith(1))*2
    do i=1,neighnodecache_count(startwith(1))
       isrc = neighnodecache(neighnodecache_start(startwith(1)) + i)
       call open_nodes(stack,istack,node(isrc),isrc,branch,startwith(2),&
                       listneigh,xyzcache,ixyzcachesize,nneigh,leaf_is_active,&
                       maxcache,xoffset,yoffset,zoffset)
    enddo
 else
    istack = istack + 1
    stack(1,istack) = irootnode
    stack(2,istack) = irootnode
    stack(3,istack) = nparents
 endif
!
!-- parallel select algorithm to check every interactions between the tree and the selected branch
!
 do while(istack > 0)
    !-- pop the stack
    idst       = stack(1,istack) ! dest node id
    isrc       = stack(2,istack) ! src node id
    idstbranch = stack(3,istack) ! dest id in branch array
    istack     = istack - 1

    if (idst == isrc) then !-- self interaction ignored (directly push onto stack)
       stackit = .true.
       xoffset = 0.
       yoffset = 0.
       zoffset = 0.
    else
       call node_interaction(node(idst),node(isrc),tree_acc2,fnode_branch(:,idstbranch),stackit,xoffset,yoffset,zoffset)
    endif

    if (stackit) then
       neighnodecount_branch(idstbranch) = neighnodecount_branch(idstbranch) + 1
       !-- if count overflow, we will not cache it during the downward pass
       if (neighnodecount_branch(idstbranch) <= maxnodecache_local) then
          neighnode_branch(neighnodecount_branch(idstbranch),idstbranch) = isrc
       endif


       call open_nodes(stack,istack,node(isrc),isrc,branch,idstbranch,&
                       listneigh,xyzcache,ixyzcachesize,nneigh,leaf_is_active,&
                       maxcache,xoffset,yoffset,zoffset)
    endif
 enddo

 !
 !-- Downward pass to accumulate on each leaf / Cache and fetch fgrav optimisation
 !
 do i=nparents,2,-1 ! parents(1) is equal to icell
    iparent = branch(i)
    ! -- Cache node if first thread to reach it or fetch fnode in memory
#ifdef GRAVITY
    if (use_cache) then
       !$omp atomic capture
       tobecached = node(iparent)%tobecached
       node(iparent)%tobecached = min(node(iparent)%tobecached,0)
       !$omp end atomic
       if (tobecached==1) then
          !always cached fnode (fnode always first cuz ncached is use as the main barrier)
          fnodecache(1:lenfgrav,iparent) = fnode_branch(1:lenfgrav,i)
          !$omp atomic write
          node(iparent)%fcached = .true.
          !$omp end atomic

          !-- store interaction list in the cache array if it fits
          if (neighnodecount_branch(i)> 0 .and. neighnodecount_branch(i) <= maxnodecache_local) then
             !$omp atomic capture
             ibase = itail_neigh
             itail_neigh = itail_neigh + neighnodecount_branch(i)
             !$omp end atomic
             if (ibase+neighnodecount_branch(i) <= size(neighnodecache)) then
                neighnodecache(ibase+1:ibase+neighnodecount_branch(i)) = neighnode_branch(1:neighnodecount_branch(i),i)
                neighnodecache_start(iparent) = ibase
                neighnodecache_count(iparent) = neighnodecount_branch(i)
                !$omp atomic write
                node(iparent)%ncached = .true.
                !$omp end atomic
                ! else
                !    print*,"overflow global !!! "
             endif
             ! else
             !    print*,"overflow local !!!",neighnodecount_branch(i),i
          endif
       else
          !$omp atomic read
          fcached = node(iparent)%fcached
          !$omp end atomic
          if (fcached) then
             !-- fetch fnode from the cache array
             fnode_branch(1:lenfgrav,i) = fnodecache(1:lenfgrav,iparent)
          endif
       endif
    endif
#else
    fcached = .true.
    tobecached=1
#endif
    call get_sep(node(iparent)%xcen,node(branch(i-1))%xcen,dx,dy,dz,xoffset,yoffset,zoffset)
    fnode = fnode_acc + fnode_branch(:,i)
    call propagate_fnode_to_node(fnode_acc,fnode,dx,dy,dz)
 enddo

 fnode = fnode_acc + fnode_branch(:,1)

end subroutine getneigh_dual

!-----------------------------------------------------------
!+
!  get the separation in 3D between two nodes of the tree
!+
!-----------------------------------------------------------
pure subroutine get_sep(x1,x2,dx,dy,dz,xoffset,yoffset,zoffset,r2)
#ifdef PERIODIC
 use boundary, only:dxbound,dybound,dzbound,hdlx,hdly,hdlz
#endif
 real, intent(in)  :: x1(3),x2(3)
 real, intent(out) :: dx,dy,dz,xoffset,yoffset,zoffset
 real, intent(out), optional :: r2

 xoffset = 0.
 yoffset = 0.
 zoffset = 0.

 dx = x2(1) - x1(1)
 dy = x2(2) - x1(2)
 dz = x2(3) - x1(3)

#ifdef PERIODIC
 if (abs(dx) > hdlx) then ! mod distances across boundary if periodic BCs
    xoffset = dxbound*SIGN(1.0,dx)
    dx = dx - xoffset
 endif
 if (abs(dy) > hdly) then
    yoffset = dybound*SIGN(1.0,dy)
    dy = dy - yoffset
 endif
 if (abs(dz) > hdlz) then
    zoffset = dzbound*SIGN(1.0,dz)
    dz = dz - zoffset
 endif
#endif

 if (present(r2)) r2 = dx*dx+dy*dy+dz*dz

end subroutine get_sep

!-----------------------------------------------------------
!+
!  get the size and rcut of two interacting nodes
!+
!-----------------------------------------------------------
pure subroutine get_node_size(node_dst,node_src,size_dst,size_src,rcut_dst,rcut_src)
 use kernel,   only:radkern
 type(kdnode), intent(in)  :: node_dst,node_src
 real,         intent(out) :: size_src,size_dst
 real,         intent(out) :: rcut_src,rcut_dst

 rcut_src = node_src%hmax*radkern
 rcut_dst = node_dst%hmax*radkern
 size_src = node_src%size
 size_dst = node_dst%size

end subroutine get_node_size

!-----------------------------------------------------------
!+
!  Taylor expand the contribution from direct parent nodes
!  to the child node centre
!+
!-----------------------------------------------------------
pure subroutine propagate_fnode_to_node(fnode,fnode_sup,dx,dy,dz)
 real, intent(in)  :: fnode_sup(lenfgrav),dx,dy,dz
 real, intent(out) :: fnode(lenfgrav)

 fnode(1)  = fnode_sup(1) + dx*(fnode_sup(4) + 0.5*(dx*fnode_sup(10) + dy*fnode_sup(11) +dz*fnode_sup(12)))& ! xx +0.5(xxx+xxy+xxz)
                          + dy*(fnode_sup(5) + 0.5*(dx*fnode_sup(11) + dy*fnode_sup(13) +dz*fnode_sup(14)))& ! xy +0.5(xxy+xyy+xyz)
                          + dz*(fnode_sup(6) + 0.5*(dx*fnode_sup(12) + dy*fnode_sup(14) +dz*fnode_sup(15)))  ! xz +0.5(xxz+xyz+xzz)
 fnode(2)  = fnode_sup(2) + dx*(fnode_sup(5) + 0.5*(dx*fnode_sup(11) + dy*fnode_sup(13) +dz*fnode_sup(14)))& ! xy +0.5(xxy+xyy+xyz)
                          + dy*(fnode_sup(7) + 0.5*(dx*fnode_sup(13) + dy*fnode_sup(16) +dz*fnode_sup(17)))& ! yy +0.5(xyy+yyy+yyz)
                          + dz*(fnode_sup(8) + 0.5*(dx*fnode_sup(14) + dy*fnode_sup(17) +dz*fnode_sup(18)))  ! yz +0.5(xyz+yyz+yyz)
 fnode(3)  = fnode_sup(3) + dx*(fnode_sup(6) + 0.5*(dx*fnode_sup(12) + dy*fnode_sup(14) +dz*fnode_sup(15)))& ! xz +0.5(xxz+xyz+xzz)
                          + dy*(fnode_sup(8) + 0.5*(dx*fnode_sup(14) + dy*fnode_sup(17) +dz*fnode_sup(18)))& ! yz +0.5(xyz+yyz+yzz)
                          + dz*(fnode_sup(9) + 0.5*(dx*fnode_sup(15) + dy*fnode_sup(18) +dz*fnode_sup(19)))  ! zz +0.5(xzz+yzz+zzz)
 fnode(4)  = fnode_sup(4) + dx*fnode_sup(10) + dy*fnode_sup(11) + dz*fnode_sup(12)                           ! xxx + xxy + xxz
 fnode(5)  = fnode_sup(5) + dx*fnode_sup(11) + dy*fnode_sup(13) + dz*fnode_sup(14)                           ! xxy + xyy + xyz
 fnode(6)  = fnode_sup(6) + dx*fnode_sup(12) + dy*fnode_sup(14) + dz*fnode_sup(15)                           ! xxz + xyz + xzz
 fnode(7)  = fnode_sup(7) + dx*fnode_sup(13) + dy*fnode_sup(16) + dz*fnode_sup(17)                           ! xyy + yyy + yyz
 fnode(8)  = fnode_sup(8) + dx*fnode_sup(14) + dy*fnode_sup(17) + dz*fnode_sup(18)                           ! xyz + yyz + yzz
 fnode(9)  = fnode_sup(9) + dx*fnode_sup(15) + dy*fnode_sup(18) + dz*fnode_sup(19)                           ! xzz + yzz + zzz
 fnode(10) = fnode_sup(10)
 fnode(11) = fnode_sup(11)
 fnode(12) = fnode_sup(12)
 fnode(13) = fnode_sup(13)
 fnode(14) = fnode_sup(14)
 fnode(15) = fnode_sup(15)
 fnode(16) = fnode_sup(16)
 fnode(17) = fnode_sup(17)
 fnode(18) = fnode_sup(18)
 fnode(19) = fnode_sup(19)
 fnode(20) = fnode_sup(20) - dx*(fnode_sup(1)+0.5*(dx*(fnode_sup(4)+(1./3.)*(dx*fnode_sup(10)+dy*fnode_sup(11)+dz*fnode_sup(12)))+&
                                                   dy*(fnode_sup(5)+(1./3.)*(dx*fnode_sup(11)+dy*fnode_sup(13)+dz*fnode_sup(14)))+&
                                                   dz*(fnode_sup(6)+(1./3.)*(dx*fnode_sup(12)+dy*fnode_sup(14)+dz*fnode_sup(15)))))&
                           - dy*(fnode_sup(2)+0.5*(dx*(fnode_sup(5)+(1./3.)*(dx*fnode_sup(11)+dy*fnode_sup(13)+dz*fnode_sup(14)))+&
                                                   dy*(fnode_sup(7)+(1./3.)*(dx*fnode_sup(13)+dy*fnode_sup(16)+dz*fnode_sup(17)))+&
                                                   dz*(fnode_sup(8)+(1./3.)*(dx*fnode_sup(14)+dy*fnode_sup(17)+dz*fnode_sup(18)))))&
                           - dz*(fnode_sup(3)+0.5*(dx*(fnode_sup(6)+(1./3.)*(dx*fnode_sup(12)+dy*fnode_sup(14)+dz*fnode_sup(15)))+&
                                                   dy*(fnode_sup(8)+(1./3.)*(dx*fnode_sup(14)+dy*fnode_sup(17)+dz*fnode_sup(18)))+&
                                                   dz*(fnode_sup(9)+(1./3.)*(dx*fnode_sup(15)+dy*fnode_sup(18)+dz*fnode_sup(19)))))

end subroutine propagate_fnode_to_node

!-----------------------------------------------------------
!+
!  return list of parents of current node
!+
!-----------------------------------------------------------
pure subroutine get_list_of_parent_nodes(inode,node,parents,nparents,startwith)
 integer,      intent(in)  :: inode
 type(kdnode), intent(in)  :: node(:)
 integer,      intent(out) :: parents(:)
 integer,      intent(out) :: nparents
 integer,      intent(out) :: startwith(2)
 integer :: j
 logical :: notfound

 j = inode
 notfound  = .true.
 nparents  = 1
 parents   = 0
 startwith = 0
 parents(nparents) = j ! set first elem to inode to use parents for propagation
 do while (node(j)%parent  /=  0)
    j = node(j)%parent
    nparents = nparents + 1
    if (node(j)%ncached .and. notfound ) then
       startwith(1) = j
       startwith(2) = nparents
       notfound = .false.
    endif
    parents(nparents) = j
 enddo

end subroutine get_list_of_parent_nodes

!-----------------------------------------------------------
!+
!  Compute node node gravity interactions
!+
!-----------------------------------------------------------
subroutine open_nodes(stack,istack,srcnode,isrc,branch,idstbranch,&
                           listneigh,xyzcache,ixyzcachesize,nneigh,leaf_is_active,&
                           maxcache,xoffset,yoffset,zoffset)
 type(kdnode), intent(in)    :: srcnode
 integer,      intent(in)    :: isrc,idstbranch
 integer,      intent(in)    :: branch(:)
 integer,      intent(in)    :: ixyzcachesize,maxcache
 integer,      intent(in)    :: leaf_is_active(:)
 integer,      intent(inout) :: listneigh(:)
 integer,      intent(inout) :: nneigh
 integer,      intent(inout) :: stack(:,:),istack
 real,         intent(inout) :: xyzcache(:,:)
 real,         intent(in)    :: xoffset,yoffset,zoffset
 integer :: ir,il,ibranchnext,idstnext
 logical :: isdstleaf

 il = srcnode%leftchild
 ir = srcnode%rightchild

 !-- find the new dst id to push onto the stack
 if (idstbranch-1>0) then !-- if not leaf
    ibranchnext = idstbranch-1
    isdstleaf   = .false.
 else
    ibranchnext = idstbranch ! leaf lowering if upper leaf
    isdstleaf   = .true.
 endif

 idstnext = branch(ibranchnext) ! new dest node id

 is_src_leaf: if (leaf_is_active(isrc) /= 0) then
    is_P2P: if (isdstleaf) then !-- P2P detected should be cached and tagged as neighbours
       call cache_neighbours(nneigh,isrc,ixyzcachesize,maxcache,listneigh,xyzcache,xoffset,yoffset,zoffset)
    else ! then you're a leaf -> leaf lowering
       istack = istack + 1
       stack(1,istack) = idstnext
       stack(2,istack) = isrc
       stack(3,istack) = ibranchnext
    endif is_P2P
 else
    if (il /= 0) then
       istack = istack + 1
       stack(1,istack) = idstnext
       stack(2,istack) = il
       stack(3,istack) = ibranchnext
    endif
    if (ir /= 0) then
       istack = istack + 1
       stack(1,istack) = idstnext
       stack(2,istack) = ir
       stack(3,istack) = ibranchnext
    endif
 endif is_src_leaf

end subroutine open_nodes

!-----------------------------------------------------------
!+
!  Test the separation between the node pair and compute
!  the interaction if needed
!+
!-----------------------------------------------------------
subroutine node_interaction(node_dst,node_src,tree_acc2,fnode,stackit,xoffset,yoffset,zoffset)
 type(kdnode), intent(in)    :: node_dst,node_src
 real,         intent(in)    :: tree_acc2
 real,         intent(inout) :: fnode(lenfgrav)
 real,         intent(out)   :: xoffset,yoffset,zoffset
 logical,      intent(out)   :: stackit
 real    :: dx,dy,dz,r2
 real    :: rcut_dst,rcut_src,rcut,rcut2
 real    :: size_dst,size_src
 logical :: wellsep,fcached
#ifdef GRAVITY
 real    :: dr1
#endif

 call get_sep(node_dst%xcen,node_src%xcen,dx,dy,dz,xoffset,yoffset,zoffset,r2)
 call get_node_size(node_dst,node_src,size_dst,size_src,rcut_dst,rcut_src)
#ifdef GRAVITY
 if (use_cache) then
    !$omp atomic read
    fcached = node_dst%fcached
    !$omp end atomic
 else
    fcached = .false.
 endif
#else
 fcached = .false.
#endif
 rcut  = max(rcut_dst,rcut_src)
 rcut2 = (size_dst+size_src+rcut)**2
 wellsep = (tree_acc2*r2 > (size_dst+size_src)**2) .and. (r2 > rcut2)

 if (wellsep) then
#ifdef GRAVITY
    if (.not.fcached) then
       dr1 = 1./sqrt(r2)
       call compute_M2L(dx,dy,dz,dr1,node_src%mass,node_src%quads,fnode)
       call add_torque_correction(dx,dy,dz,dr1,node_dst%mass,node_src%mass, &
                                  node_dst%octs,node_src%octs,fnode)
    endif
#endif
    stackit = .false.
 else
    stackit = .true.
 endif

end subroutine node_interaction

!-----------------------------------------------------------
!+
!  Compute the Taylor expansion coeffs between the node
!  centres using the quadrupole moments (p=3) (Dehnen 2002)
!+
!-----------------------------------------------------------
pure subroutine compute_M2L(dx,dy,dz,dr1,q0,quads,fnode)
 real, intent(in)    :: dx,dy,dz,dr1,q0
 real, intent(in)    :: quads(9)
 real, intent(inout) :: fnode(lenfgrav)
 real :: qx,qy,qz,qxx,qxy,qxz,qyy,qyz,qzz,dx2,dx3,dy2,dy3,dz2,dz3
 real :: dr12,D3(10),D2(6),D1(3),g0,g1,g2,g3,g2dx,g2dy,g2dz

! note: dr == 1/sqrt(r2)
 dr12 = dr1*dr1
 dx2  = dx*dx
 dx3  = dx*dx2
 dy2  = dy*dy
 dy3  = dy*dy2
 dz2  = dz*dz
 dz3  = dz*dz2
 ! be careful with the sign of your Green's function, it can mess up everything.
 ! We switched multiple signs here to match the Phantom sign convention
 g0   =  dr1
 g1   = -1.*dr12*g0
 g2   = -3.*dr12*g1
 g3   = -5.*dr12*g2
 g2dx = g2 * dx
 g2dy = g2 * dy
 g2dz = g2 * dz

 !D1, D2, D3 verified and agree with shamrock to float precision
 D3(1)  = 3. * g2dx + g3 * dx3    ! xxx
 D3(2)  = g2dy + g3 * dx2 * dy    ! xxy
 D3(3)  = g2dz + g3 * dx2 * dz    ! xxz
 D3(4)  = g2dx + g3 * dy2 * dx    ! xyy
 D3(5)  = g3 * dx * dy * dz       ! xyz
 D3(6)  = g2dx + g3 * dz2 * dx    ! xzz
 D3(7)  = 3. * g2dy + g3 * dy3    ! yyy
 D3(8)  = g2dz + g3 * dy2 * dz    ! yyz
 D3(9)  = g2dy + g3 * dz2 * dy    ! yzz
 D3(10) = 3. * g2dz + g3 * dz3    ! zzz

 D2(1)  = g1 + g2 * dx2 ! xx
 D2(2)  = g2dx * dy     ! xy
 D2(3)  = g2dx * dz     ! xz
 D2(4)  = g1 + g2 * dy2 ! yy
 D2(5)  = g2dy * dz     ! yz
 D2(6)  = g1 + g2 * dz2 ! zz

 D1(1)  = g1*dx
 D1(2)  = g1*dy
 D1(3)  = g1*dz

 qx  = quads(1)
 qy  = quads(2)
 qz  = quads(3)
 qxx = quads(4)
 qxy = quads(5)
 qxz = quads(6)
 qyy = quads(7)
 qyz = quads(8)
 qzz = quads(9)

 fnode(1)  = fnode(1)  + (D1(1)*q0  + D2(1)*qx + D2(2)*qy + D2(3)*qz +&
                     0.5*(D3(1)*qxx + 2.*(D3(2)*qxy + D3(3)*qxz + D3(5)*qyz) + D3(4)*qyy + D3(6)*qzz ))    ! C¹_x
 fnode(2)  = fnode(2)  + (D1(2)*q0  + D2(2)*qx + D2(4)*qy + D2(5)*qz +&
                     0.5*(D3(2)*qxx + 2.*(D3(4)*qxy + D3(5)*qxz + D3(8)*qyz) + D3(7)*qyy + D3(9)*qzz ))    ! C¹_y
 fnode(3)  = fnode(3)  + (D1(3)*q0  + D2(3)*qx + D2(5)*qy + D2(6)*qz +&
                     0.5*(D3(3)*qxx + 2.*(D3(5)*qxy + D3(6)*qxz + D3(9)*qyz) + D3(8)*qyy + D3(10)*qzz))   ! C¹_z
 fnode(4)  = fnode(4)  - (D2(1) * q0 + D3(1)*qx + D3(2)*qy + D3(3)*qz)  ! C²_xx
 fnode(5)  = fnode(5)  - (D2(2) * q0 + D3(2)*qx + D3(4)*qy + D3(5)*qz)! C²_xy
 fnode(6)  = fnode(6)  - (D2(3) * q0 + D3(3)*qx + D3(5)*qy + D3(6)*qz)! C²_xz
 fnode(7)  = fnode(7)  - (D2(4) * q0 + D3(4)*qx + D3(7)*qy + D3(8)*qz)! C²_yy
 fnode(8)  = fnode(8)  - (D2(5) * q0 + D3(5)*qx + D3(8)*qy + D3(9)*qz)! C²_yz
 fnode(9)  = fnode(9)  - (D2(6) * q0 + D3(6)*qx + D3(9)*qy + D3(10)*qz)! C²_zz
 fnode(10) = fnode(10) + D3(1) * q0    ! C³_xxx
 fnode(11) = fnode(11) + D3(2) * q0    ! C³_xxy
 fnode(12) = fnode(12) + D3(3) * q0    ! C³_xxz
 fnode(13) = fnode(13) + D3(4) * q0    ! C³_xyy
 fnode(14) = fnode(14) + D3(5) * q0    ! C³_xyz
 fnode(15) = fnode(15) + D3(6) * q0    ! C³_xzz
 fnode(16) = fnode(16) + D3(7) * q0    ! C³_yyy
 fnode(17) = fnode(17) + D3(8) * q0    ! C³_yyz
 fnode(18) = fnode(18) + D3(9) * q0    ! C³_yzz
 fnode(19) = fnode(19) + D3(10)* q0    ! C³_zzz
 fnode(20) = fnode(20) + g0*q0 + (D1(1)*qx + D1(2)*qy + D1(3)*qz)  + &
                         0.5*(D2(1)*qxx + D2(4)*qyy + D2(6)*qzz + 2*(D2(2)*qxy + D2(3)*qxz + D2(5)*qyz))! C⁰ (potential)

end subroutine compute_M2L

#ifdef GRAVITY
!----------------------------------------------------------------
!+
!  Accumulate quadrupole and octupole moments of a particle
!  about the node centre of mass (extensive Cartesian form).
!+
!----------------------------------------------------------------
pure subroutine add_node_moments(pmassi,dx,dy,dz,quads,octs)
 real, intent(in)    :: pmassi,dx,dy,dz
 real, intent(inout) :: quads(9),octs(10)
 real :: dx2,dy2,dz2

 dx2 = dx*dx
 dy2 = dy*dy
 dz2 = dz*dz
 quads(1) = quads(1) + pmassi*dx  ! Q_x
 quads(2) = quads(2) + pmassi*dy  ! Q_y
 quads(3) = quads(3) + pmassi*dz  ! Q_z
 quads(4) = quads(4) + pmassi*dx2          ! Q_xx
 quads(5) = quads(5) + pmassi*dx*dy        ! Q_xy
 quads(6) = quads(6) + pmassi*dx*dz        ! Q_xz
 quads(7) = quads(7) + pmassi*dy2          ! Q_yy
 quads(8) = quads(8) + pmassi*dy*dz        ! Q_yz
 quads(9) = quads(9) + pmassi*dz2          ! Q_zz
 octs(1)  = octs(1)  + pmassi*dx2*dx       ! xxx
 octs(2)  = octs(2)  + pmassi*dx2*dy       ! xxy
 octs(3)  = octs(3)  + pmassi*dx2*dz       ! xxz
 octs(4)  = octs(4)  + pmassi*dx*dy2       ! xyy
 octs(5)  = octs(5)  + pmassi*dx*dy*dz     ! xyz
 octs(6)  = octs(6)  + pmassi*dx*dz2       ! xzz
 octs(7)  = octs(7)  + pmassi*dy2*dy       ! yyy
 octs(8)  = octs(8)  + pmassi*dy2*dz       ! yyz
 octs(9)  = octs(9)  + pmassi*dy*dz2       ! yzz
 octs(10) = octs(10) + pmassi*dz2*dz       ! zzz

end subroutine add_node_moments

!----------------------------------------------------------------
!+
!  Marcello (2017) TCO torque correction: add a constant
!  acceleration Fc/M_dst to the destination cell so the net
!  cell-cell torque vanishes, while keeping equal-and-opposite
!  forces. Uses the pruned D'_ijkl contraction (his Eq. 19).
!+
!----------------------------------------------------------------
pure subroutine add_torque_correction(dx,dy,dz,dr1,mass_dst,mass_src,octs_dst,octs_src,fnode)
 real, intent(in)    :: dx,dy,dz,dr1,mass_dst,mass_src
 real, intent(in)    :: octs_dst(10),octs_src(10)
 real, intent(inout) :: fnode(lenfgrav)
 real :: s(10)
 real :: sxkk,sykk,szkk,sxrr,syrr,szrr
 real :: r5i,r7i,fac

 if (mass_dst <= 0. .or. mass_src <= 0.) return

 ! S_jkl = M_dst,jkl * M_src - M_dst * M_src,jkl
 s(:) = octs_dst*mass_src - mass_dst*octs_src

 ! traces S_i,kk
 sxkk = s(1) + s(4) + s(6)
 sykk = s(2) + s(7) + s(9)
 szkk = s(3) + s(8) + s(10)

 ! S_iab R_a R_b  (R is the node separation; even in R so dest-src is fine)
 sxrr = s(1)*dx*dx + s(4)*dy*dy + s(6)*dz*dz + 2.*(s(2)*dx*dy + s(3)*dx*dz + s(5)*dy*dz)
 syrr = s(2)*dx*dx + s(7)*dy*dy + s(9)*dz*dz + 2.*(s(4)*dx*dy + s(5)*dx*dz + s(8)*dy*dz)
 szrr = s(3)*dx*dx + s(8)*dy*dy + s(10)*dz*dz + 2.*(s(5)*dx*dy + s(6)*dx*dz + s(9)*dy*dz)

 r5i = dr1**5
 r7i = r5i*dr1*dr1
 ! Appendix Eq. 38 at P=3 uses 1/(n!(P-n)!) = 1/3!, not the 1/2 of Eq. 15.
 ! Combined with the D' contraction this is 3/2 rather than 9/2.
 fac = 1.5

 ! Fc_i = (3/2) (S_ikk/R^5 - 5 S_iab R_a R_b / R^7); add Fc/M_dst
 fnode(1) = fnode(1) - fac*(sxkk*r5i - 5.*sxrr*r7i)/mass_dst
 fnode(2) = fnode(2) - fac*(sykk*r5i - 5.*syrr*r7i)/mass_dst
 fnode(3) = fnode(3) - fac*(szkk*r5i - 5.*szrr*r7i)/mass_dst

end subroutine add_torque_correction
#endif

!----------------------------------------------------------------
!+
!  Internal subroutine to compute the Taylor-series expansion
!  of the gravitational force, given the force acting on the
!  centre of the node and its derivatives
!
! INPUT:
!   fnode: array containing force on node due to distant nodes
!          and first derivatives of f (i.e. Jacobian matrix)
!          and second derivatives of f (i.e. Hessian matrix)
!   dx,dy,dz: offset of the particle from the node centre of mass
!
! OUTPUT:
!   fxi,fyi,fzi : gravitational force at the new position
!+
!----------------------------------------------------------------
pure subroutine expand_fgrav_in_taylor_series(fnode,dx,dy,dz,fxi,fyi,fzi,poti)
 real, intent(in)  :: fnode(lenfgrav)
 real, intent(in)  :: dx,dy,dz
 real, intent(out) :: fxi,fyi,fzi,poti
 real :: dfxx,dfxy,dfxz,dfyy,dfyz,dfzz
 real :: d2fxxx,d2fxxy,d2fxxz,d2fxyy,d2fxyz,d2fxzz,d2fyyy,d2fyyz,d2fyzz,d2fzzz

 fxi = fnode(1)
 fyi = fnode(2)
 fzi = fnode(3)
 dfxx = fnode(4)
 dfxy = fnode(5)
 dfxz = fnode(6)
 dfyy = fnode(7)
 dfyz = fnode(8)
 dfzz = fnode(9)
 d2fxxx = fnode(10)
 d2fxxy = fnode(11)
 d2fxxz = fnode(12)
 d2fxyy = fnode(13)
 d2fxyz = fnode(14)
 d2fxzz = fnode(15)
 d2fyyy = fnode(16)
 d2fyyz = fnode(17)
 d2fyzz = fnode(18)
 d2fzzz = fnode(19)
 poti = fnode(20)

 fxi = fxi   + dx*(dfxx + 0.5*(dx*d2fxxx + dy*d2fxxy + dz*d2fxxz)) &
             + dy*(dfxy + 0.5*(dx*d2fxxy + dy*d2fxyy + dz*d2fxyz)) &
             + dz*(dfxz + 0.5*(dx*d2fxxz + dy*d2fxyz + dz*d2fxzz))
 fyi = fyi   + dx*(dfxy + 0.5*(dx*d2fxxy + dy*d2fxyy + dz*d2fxyz)) &
             + dy*(dfyy + 0.5*(dx*d2fxyy + dy*d2fyyy + dz*d2fyyz)) &
             + dz*(dfyz + 0.5*(dx*d2fxyz + dy*d2fyyz + dz*d2fyzz))
 fzi = fzi   + dx*(dfxz + 0.5*(dx*d2fxxz + dy*d2fxyz + dz*d2fxzz)) &
             + dy*(dfyz + 0.5*(dx*d2fxyz + dy*d2fyyz + dz*d2fyzz)) &
             + dz*(dfzz + 0.5*(dx*d2fxzz + dy*d2fyzz + dz*d2fzzz))
 ! Minus sign here as we are shifted of 1 in the (-1)^k compared to force
 poti = poti - dx*(fxi - 0.5*(dx*dfxx + dy*dfxy + dz*dfxz)) &
             - dy*(fyi - 0.5*(dx*dfxy + dy*dfyy + dz*dfyz)) &
             - dz*(fzi - 0.5*(dx*dfxz + dy*dfyz + dz*dfzz))

end subroutine expand_fgrav_in_taylor_series

!-----------------------------------------------
!+
!  Routine to update a constructed tree
!  Note: current version ONLY works if
!  tree is built to < maxlevel_indexed
!  That is, it relies on the 2^n style tree
!  indexing to sweep out each level
!+
!-----------------------------------------------
subroutine revtree(node, xyzh, leaf_is_active, ncells)
 use dim,  only:maxp,use_apr,ind_timesteps
 use part, only:maxphase,iphase,igas,massoftype,iamtype,aprmassoftype,&
                apr_level,iactive,treecache,isdead_or_accreted
 use io,   only:fatal
 type(kdnode),    intent(inout) :: node(:) !ncellsmax+1)
 real,            intent(in)    :: xyzh(:,:)
 integer,         intent(inout) :: leaf_is_active(:) !ncellsmax+1)
 integer(kind=8), intent(in)    :: ncells
 real :: hmax, r2max
 real :: xi, yi, zi, hi
 real :: dx, dy, dz, dr2
#ifdef GRAVITY
 real :: quads(9)
 real :: octs(10)
#endif
 integer :: inode, ipart, ipartidx, i, nptot
 real :: pmassi, totmass
 real :: x0(3)
 real :: xcofm, ycofm, zcofm, fac, dfac
 logical :: nodeisactive
 itail_neigh = 0
 pmassi = massoftype(igas)

 ! find maximum index in inodeparts that we need to update in treecache
 nptot = 0
 do i=1,int(ncells)
    if (i > 1 .and. node(i)%parent == 0) cycle
    if (inoderange(1,i) > 0 .and. inoderange(2,i) >= inoderange(1,i)) then
       nptot = max(nptot, inoderange(2,i))
    endif
 enddo

 ! update treecache for particles in the tree only
 ! mark dead/accreted particles by setting treecache(4,i) negative
 !$omp parallel default(none) &
 !$omp shared(nptot,inodeparts,xyzh,iphase,apr_level) &
 !$omp shared(massoftype,aprmassoftype,treecache) &
 !$omp shared(maxphase,maxp) &
 !$omp private(i,ipartidx)
 !$omp do schedule(static)
 do i=1,nptot
    if (inodeparts(i) == 0) cycle
    ipartidx = abs(inodeparts(i))
    treecache(1:4,i) = xyzh(1:4,ipartidx)
    ! compute and store mass
    if (maxphase==maxp) then
       if (use_apr) then
          treecache(5,i) = aprmassoftype(iamtype(iphase(ipartidx)),apr_level(ipartidx))
       else
          treecache(5,i) = massoftype(iamtype(iphase(ipartidx)))
       endif
    elseif (use_apr) then
       treecache(5,i) = aprmassoftype(igas,apr_level(ipartidx))
    else
       treecache(5,i) = massoftype(igas)
    endif
 enddo
 !$omp enddo
 !$omp end parallel

!$omp parallel default(none) &
!$omp shared(maxp,maxphase) &
!$omp shared(ncells) &
!$omp shared(node,inoderange,inodeparts,treecache,leaf_is_active) &
!$omp private(hmax,r2max,xi,yi,zi,hi) &
!$omp private(dx,dy,dz,dr2,inode,ipart,x0) &
!$omp private(xcofm,ycofm,zcofm,fac,dfac,nodeisactive) &
#ifdef GRAVITY
!$omp private(quads,octs) &
#endif
!$omp firstprivate(pmassi) &
!$omp private(totmass)
!$omp do schedule(guided)
 over_nodes: do inode=1,int(ncells)
    if (inode > 1 .and. node(inode)%parent == 0) cycle
    ! initialize node properties
    node(inode)%xcen(:) = 0.
    node(inode)%size    = 0.
    node(inode)%hmax    = 0.
#ifdef GRAVITY
    node(inode)%mass    = 0.
    node(inode)%quads(:)= 0.
    node(inode)%octs(:) = 0.
#endif
    ! initialize leaf_is_active (will be set for leaf nodes below)
    leaf_is_active(inode) = 0

    ! check if node has particles
    if (inoderange(1,inode) <= 0 .or. inoderange(2,inode) < inoderange(1,inode)) cycle over_nodes

    ! find centre of mass from particle list using same algorithm as maketree
    ! also check for active particles and compute hmax during this loop
    xcofm = 0.
    ycofm = 0.
    zcofm = 0.
    totmass = 0.0
    hmax = 0.
    dfac = 1.
    if (pmassi > 0.) then
       dfac = 1./pmassi
    endif
    nodeisactive = .false.
    do ipart = inoderange(1,inode), inoderange(2,inode)
       if (inodeparts(ipart) == 0) cycle
       xi = treecache(1,ipart)
       yi = treecache(2,ipart)
       zi = treecache(3,ipart)
       hi = treecache(4,ipart)
       ! check condition after loading (dead/accreted particles have hi <= 0)
       if (hi <= 0.) cycle
       hi = abs(hi)
       pmassi = treecache(5,ipart)
       ! check for active particles (for leaf_is_active flag)
       if (ind_timesteps .and. .not. nodeisactive) then
          if (inodeparts(ipart) > 0) nodeisactive = .true.
       endif
       fac = pmassi*dfac
       xcofm = xcofm + fac*xi
       ycofm = ycofm + fac*yi
       zcofm = zcofm + fac*zi
       totmass = totmass + pmassi
       hmax = max(hi, hmax)
    enddo
    if (.not. ind_timesteps) nodeisactive = .true.

    if (totmass <= 0.0) cycle over_nodes

    x0(1) = xcofm/(totmass*dfac)
    x0(2) = ycofm/(totmass*dfac)
    x0(3) = zcofm/(totmass*dfac)

    ! update cell size and quads
    r2max = 0.
#ifdef GRAVITY
    quads = 0.
    octs  = 0.
#endif
    do ipart = inoderange(1,inode), inoderange(2,inode)
       ! load all treecache values sequentially (1,2,3,4,5) for cache efficiency
       xi = treecache(1,ipart)
       yi = treecache(2,ipart)
       zi = treecache(3,ipart)
       hi = treecache(4,ipart)
       ! check condition after loading (dead/accreted particles have hi <= 0)
       if (hi <= 0.) cycle
       pmassi = treecache(5,ipart)
       dx = xi - x0(1)
       dy = yi - x0(2)
       dz = zi - x0(3)
       dr2 = dx*dx + dy*dy + dz*dz
       r2max = max(dr2, r2max)
#ifdef GRAVITY
       call add_node_moments(pmassi,dx,dy,dz,quads,octs)
#endif
    enddo

    node(inode)%xcen(1) = x0(1)
    node(inode)%xcen(2) = x0(2)
    node(inode)%xcen(3) = x0(3)
    node(inode)%size = sqrt(r2max) + epsilon(r2max)
    node(inode)%hmax = hmax
#ifdef GRAVITY
    node(inode)%mass = totmass
    node(inode)%quads = quads
    node(inode)%octs  = octs
    node(inode)%tobecached = 1
    node(inode)%fcached = .false.
    node(inode)%ncached = .false.
#endif

    ! set leaf_is_active flag for leaf nodes (matching maketree behavior)
    if (node(inode)%leftchild == 0 .and. node(inode)%rightchild == 0) then
       if (ind_timesteps) then
          if (nodeisactive) then
             leaf_is_active(inode) = 1
          else
             leaf_is_active(inode) = -1
          endif
       else
          leaf_is_active(inode) = 1
       endif
    endif
 enddo over_nodes
!$omp enddo
!$omp end parallel

end subroutine revtree

!--------------------------------------------------------------------------------
!+
!  Routine to build the global level tree
!+
!-------------------------------------------------------------------------------
subroutine maketreeglobal(nodeglobal,node,nodemap,globallevel,refinelevels,xyzh,&
                          np,cellatid,leaf_is_active,ncells,apr_tree,nptmass,xyzmh_ptmass)
 use io,           only:fatal,warning,id,nprocs,master
 use mpiutils,     only:reduceall_mpi
 use mpibalance,   only:balancedomains
 use mpitree,      only:tree_sync,tree_bcast
 use part,         only:isdead_or_accreted,iactive,ibelong,isink,massoftype,igas,&
                        iamtype,maxphase,maxp,aprmassoftype,apr_level,ihsoft
 use timing,       only:increment_timer,get_timings,itimer_balance
 use dim,          only:ind_timesteps

 type(kdnode),    intent(out)   :: nodeglobal(:)    ! ncellsmax+1
 type(kdnode),    intent(out)   :: node(:)          ! ncellsmax+1
 integer,         intent(out)   :: nodemap(:)       ! ncellsmax+1
 integer,         intent(out)   :: globallevel
 integer,         intent(out)   :: refinelevels
 integer,         intent(inout) :: np
 real,            intent(inout) :: xyzh(:,:)
 integer,         intent(out)   :: cellatid(:)      ! ncellsmax+1
 integer,         intent(out)   :: leaf_is_active(:)  ! ncellsmax+1)
 integer(kind=8), intent(out)   :: ncells
 logical,         intent(in)    :: apr_tree
 integer,         intent(in),    optional :: nptmass
 real,            intent(inout), optional :: xyzmh_ptmass(:,:)
 real                              :: xmini(3),xmaxi(3)
 real                              :: xminl(3),xmaxl(3)
 real                              :: xminr(3),xmaxr(3)
 integer                           :: minlevel, maxlevel
 integer                           :: idleft, idright
 integer                           :: groupsize,ifirstingroup,groupsplit
 type(kdnode)                      :: mynode(1)
 integer                           :: nl, nr
 integer                           :: il, ir, iself, parent
 integer                           :: level
 integer                           :: nnodestart, nnodeend,locstart,locend
 integer                           :: npcounter
 integer                           :: i, k, offset, roffset, roffset_prev, coffset
 integer                           :: inode
 integer                           :: npnode
 logical                           :: wassplit,sinktree
 real(kind=4)                      :: t1,t2,tcpu1,tcpu2

 sinktree = .false.
 if (present(nptmass).and.present(xyzmh_ptmass)) sinktree=.true.
 parent = 0
 iself = irootnode
 leaf_is_active = 0

 ! root is level 0
 globallevel = int(ceiling(log(real(nprocs)) / log(2.0)))

 minlevel = maxdepth - 1
 maxlevel = 0

 levels: do level = 0, globallevel
    groupsize = 2**(globallevel - level)
    ifirstingroup = (id / groupsize) * groupsize
    if (level == 0) then
       if (sinktree) then
          call construct_root_node(np,npcounter,irootnode,xmini,xmaxi,leaf_is_active,xyzh,&
                                   xyzmh_ptmass,nptmass)
       else
          call construct_root_node(np,npcounter,irootnode,xmini,xmaxi,leaf_is_active,xyzh)
       endif
    else
       npcounter = npnode
    endif
    if (sinktree) then
       call construct_node(mynode(1), iself, parent, level, xmini, xmaxi, npcounter, .false., &
                           il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr,ncells, leaf_is_active, &
                           minlevel, maxlevel, wassplit,.true.,apr_tree,xyzmh_ptmass)
    else
       call construct_node(mynode(1), iself, parent, level, xmini, xmaxi, npcounter, .false., &
                        il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr,ncells, leaf_is_active, &
                        minlevel, maxlevel, wassplit,.true.,apr_tree)
    endif

    if (.not.wassplit) then
       call fatal('maketreeglobal','insufficient particles for splitting at the global level: '// &
            'use more particles or less MPI threads')
    endif

    ! set which tree child this proc will belong to next
    groupsplit = ifirstingroup + (groupsize / 2)

    ! record parent for next round
    parent = iself

    ! which half of the tree this task is on
    if (id < groupsplit) then
       ! i for the next node we construct
       iself = il
       ! the left and right task IDs
       idleft = id
       idright = id + 2**(globallevel - level - 1)
       xmini = xminl
       xmaxi = xmaxl
    else
       iself = ir
       idleft = id - 2**(globallevel - level - 1)
       idright = id
       xmini = xminr
       xmaxi = xmaxr
    endif
    if (sinktree) then
       if (nptmass>0) then
          ibelong(maxpsph+1:maxpsph+nptmass) = -1
       endif
    endif
    if (npcounter > 0) then
       do i = inoderange(1,il), inoderange(2,il)
          ibelong(abs(inodeparts(i))) = idleft
       enddo
       do i = inoderange(1,ir), inoderange(2,ir)
          ibelong(abs(inodeparts(i))) = idright
       enddo
    endif

    call get_timings(t1,tcpu1)
    ! move particles to where they belong
    call balancedomains(np)
    call get_timings(t2,tcpu2)
    if (sinktree) ibelong(maxpsph+1:maxpsph+nptmass) = int(reduceall_mpi("max", ibelong(maxpsph+1:maxpsph+nptmass)))
    call increment_timer(itimer_balance,t2-t1,tcpu2-tcpu1)
    ! move particles from old array
    ! this is a waste of time, but maintains compatibility
    npnode = 0
    do i=1,np
       npnode = npnode + 1
       !
       ! tag inactive particles with negative index
       ! in the particle list for the node
       !
       if (ind_timesteps) then
          if (iactive(iphase(i))) then
             inodeparts(npnode) = i
          else
             inodeparts(npnode) = -i
          endif
       else
          inodeparts(npnode) = i
       endif
       treecache(1:4,npnode) = xyzh(1:4,i)
       if (maxphase==maxp) then
          if (use_apr) then
             treecache(5,npnode) = aprmassoftype(iamtype(iphase(i)),apr_level(i))
          else
             treecache(5,npnode) = massoftype(iamtype(iphase(i)))
          endif
       elseif (use_apr) then
          treecache(5,npnode) = aprmassoftype(igas,apr_level(i))
       else
          treecache(5,npnode) = massoftype(igas)
       endif
    enddo
    if (sinktree) then
       if (nptmass > 0) then
          do i=1,nptmass
             if (ibelong(maxpsph + i) /= id) cycle
             if (xyzmh_ptmass(4,i) < 0.) cycle ! dead sink particle
             npnode = npnode + 1
             inodeparts(npnode) = maxpsph + i
             treecache(1:3,npnode) = xyzmh_ptmass(1:3,i)
             treecache(4,npnode)   = xyzmh_ptmass(ihsoft,i)
             treecache(5,npnode)   = xyzmh_ptmass(4,i)
          enddo
       endif
    endif

    ! set all particles to belong to this node
    inoderange(1,iself) = 1
    inoderange(2,iself) = npnode

    ! range of newly written tree
    nnodestart = 2**level
    nnodeend = 2**(level + 1) - 1

    ! synchronize tree with other owners if this proc is the first in group
    call tree_sync(mynode,1,nodeglobal(nnodestart:nnodeend),nprocs/groupsize,ifirstingroup,level)

    ! at level 0, tree_sync already 'broadcasts'
    if (level > 0) then
       ! tree broadcast to non-owners
       call tree_bcast(nodeglobal(nnodestart:nnodeend), nnodeend - nnodestart + 1, level)
    endif

 enddo levels

 ! local tree
 if (sinktree) then
    call maketree(node,xyzh,np,leaf_is_active,ncells,apr_tree,refinelevels,nptmass,xyzmh_ptmass)
 else
    call maketree(node,xyzh,np,leaf_is_active,ncells,apr_tree,refinelevels)
 endif

 ! tree refinement
 refinelevels = int(reduceall_mpi('min',refinelevels),kind=kind(refinelevels))
 roffset_prev = 1

 irefine = 0
 do i = 1,refinelevels
    offset = 2**(globallevel + i)
    roffset = 2**i

    nnodestart = offset
    nnodeend   = 2*nnodestart-1

    if (nnodeend > ncellsmax) call fatal('kdtree', 'global tree refinement has exceeded ncellsmax')

    locstart   = roffset
    locend     = 2*locstart-1

    ! index shift the node to the global level
    do k = roffset,2*roffset-1
       refinementnode(k) = node(k)
       coffset = refinementnode(k)%parent - roffset_prev

       refinementnode(k)%parent = 2**(globallevel + i - 1) + id * roffset_prev + coffset

       if (i /= refinelevels) then
          refinementnode(k)%leftchild  = 2**(globallevel + i + 1) + 2*id*roffset + 2*(k - roffset)
          refinementnode(k)%rightchild = refinementnode(k)%leftchild + 1
       else
          refinementnode(k)%leftchild = 0
          refinementnode(k)%rightchild = 0
       endif
    enddo

    roffset_prev = roffset
    ! sync, replacing level with globallevel, since all procs will get synced
    ! and deeper comms do not exist
    call tree_sync(refinementnode(locstart:locend),roffset, &
                   nodeglobal(nnodestart:nnodeend),nnodestart-nnodeend, &
                   id,globallevel)

    ! get the mapping from the local tree to the global tree, for future hmax updates
    do inode = locstart,locend
       nodemap(inode) = nnodestart + (id * roffset) + (inode - locstart)
    enddo
 enddo
!  The index up to which the local tree is copied to the global tree
 irefine = 2*roffset-1

 ! cellatid is zero by default
 cellatid = 0
 do i = 1,nprocs
    offset = 2**(globallevel+refinelevels)
    roffset = 2**refinelevels
    do k = 1,roffset
       cellatid(offset + (i - 1) * roffset + (k - 1)) = i
    enddo
 enddo

end subroutine maketreeglobal

end module kdtree
