!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module neighkdtree
!
! This module contains all routines required for
!  tree based neighbour-finding
!
!  THIS VERSION USES A K-D TREE
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - tree_accuracy : *tree opening criterion (0.0-1.0)*
!
! :Dependencies: allocutils, boundary, dim, dtypekdtree, infile_utils, io,
!   kdtree, kernel, mpiutils, part
!
 use dim,          only:ncellsmax,ncellsmaxglobal
 use dtypekdtree,  only:kdnode
 implicit none

 integer,               allocatable :: cellatid(:)
 integer,               allocatable :: nodemap(:)
 type(kdnode),          allocatable :: nodeglobal(:)
 type(kdnode), public,  allocatable :: node(:)
 integer,      public,  allocatable :: leaf_is_active(:) ! : 0 internal node or empty cell, : 1 active cell, :- inactive cell
 integer,      public,  allocatable :: active_leaves(:)  ! the cells with leaf_is_active > 0, in order (set by build_tree)
 integer,      public               :: nactive_leaves = 0
 integer,      public , allocatable :: listneigh(:)
 integer,      public , allocatable :: listneigh_global(:)
!$omp threadprivate(listneigh)
 integer(kind=8), public            :: ncells
 real, public                       :: dxcell
 real, public                       :: dcellx = 0.,dcelly = 0.,dcellz = 0.
 logical, public                    :: use_dualtree = .true.
 ! MPI: dual tree walk from the global tree, then in rounds of kslab_mpi levels on
 ! the local tree, the remote nodes being exchanged between the rounds
 integer, public                    :: kslab_mpi = 3
 ! check at each round that the remote pairs are mirrored (debug, costs 4 collectives)
 logical, public                    :: check_dualtree_mpi = .true.
 ! result of the walk in rounds for the leaves, as records (a leaf with remote src
 ! left to open is walked again in the next round, adding a record): latest record of
 ! each leaf (0 if not walked), previous record of the same leaf, its local src leaves
 ! leafsrc(leafsrc_start+1:leafsrc_start+leafsrc_count), its remote src leaves (owner,id)
 ! leafrem(:,leafrem_start+1:leafrem_start+leafrem_count), and the expansion
 integer, allocatable :: leafslot(:),leafrec_prev(:),leafrec_cell(:)
 integer, allocatable :: leafsrc(:),leafsrc_start(:),leafsrc_count(:)
 integer, allocatable :: leafrem(:,:),leafrem_start(:),leafrem_count(:)
 real,    allocatable :: fnode_leafcell(:,:)
 integer              :: nleafslots = 0, nleafsrc = 0, nleafrem = 0
 ! remote leaves received as ghost particles (after npart in the particle arrays):
 ! (owner,leaf) sorted, first ghost particle, number of particles, centre of the leaf
 integer(kind=8), allocatable :: ghostkey(:)
 integer,         allocatable :: ghostfirst(:),ghostcount(:)
 real,            allocatable :: ghostxcen(:,:),ghostmass(:)
 integer                      :: nghostleaves = 0, ighostbase = 0
 integer                            :: globallevel,refinelevels

 public :: allocate_neigh, deallocate_neigh
 public :: build_tree, get_neighbour_list, write_options_tree, read_options_tree
 public :: get_distance_from_centre_of_mass, getneigh_pos
 public :: set_hmaxcell,get_hmaxcell
 public :: get_cell_location
 public :: sync_hmax_mpi
 public :: get_global_pairs,dualwalk_rounds,start_local_rounds,get_leaf_walk
 public :: get_remote_leaves,set_ghost_leaves

 private

contains

!-----------------------------------------------------------------------
!+
!  allocate memory for the neighbour list
!+
!-----------------------------------------------------------------------
subroutine allocate_neigh
 use allocutils, only:allocate_array
 use kdtree,     only:allocate_kdtree
 use dim,        only:maxp

 call allocate_array('cellatid',       cellatid,       ncellsmaxglobal+1 )
 call allocate_array('leaf_is_active', leaf_is_active, ncellsmax+1       )
 call allocate_array('active_leaves',  active_leaves,  ncellsmax+1       )
 call allocate_array('nodeglobal',     nodeglobal,     ncellsmaxglobal+1 )
 call allocate_array('node',           node,           ncellsmax+1       )
 call allocate_array('nodemap',        nodemap,        ncellsmax+1       )
 call allocate_kdtree()
!$omp parallel
 call allocate_array('listneigh',listneigh,maxp)
!$omp end parallel
 call allocate_array('listneigh_global',listneigh_global,maxp)

end subroutine allocate_neigh

!-----------------------------------------------------------------------
!+
!  deallocate memory for the neighbour list
!+
!-----------------------------------------------------------------------
subroutine deallocate_neigh
 use kdtree,   only:deallocate_kdtree

 if (allocated(cellatid)) deallocate(cellatid)
 if (allocated(leaf_is_active)) deallocate(leaf_is_active)
 if (allocated(active_leaves)) deallocate(active_leaves)
 if (allocated(nodeglobal)) deallocate(nodeglobal)
 if (allocated(node)) deallocate(node)
 if (allocated(nodemap)) deallocate(nodemap)
 if (allocated(leafslot)) deallocate(leafslot)
 if (allocated(ghostkey)) deallocate(ghostkey,ghostfirst,ghostcount,ghostxcen)
 if (allocated(ghostmass)) deallocate(ghostmass)
 if (allocated(leafsrc)) deallocate(leafsrc,leafsrc_start,leafsrc_count,leafrec_prev,leafrec_cell,fnode_leafcell,&
                                    leafrem,leafrem_start,leafrem_count)
!$omp parallel
 if (allocated(listneigh)) deallocate(listneigh)
!$omp end parallel
 if (allocated(listneigh_global)) deallocate(listneigh_global)
 call deallocate_kdtree()

end subroutine deallocate_neigh

!-----------------------------------------------------------------------
!+
!  get the hmax value of a cell
!+
!-----------------------------------------------------------------------
subroutine get_hmaxcell(inode,hmaxcell)
 integer, intent(in)  :: inode
 real,    intent(out) :: hmaxcell

 hmaxcell = node(inode)%hmax

end subroutine get_hmaxcell

!-----------------------------------------------------------------------
!+
!  set the hmax value of a cell and propagate the value up the tree
!+
!-----------------------------------------------------------------------
subroutine set_hmaxcell(inode,hmaxcell)
 integer, intent(in) :: inode
 real,    intent(in) :: hmaxcell
 integer :: n
 real    :: hmaxn

 n = inode
 node(n)%hmax = hmaxcell

 ! walk tree up, stopping at the first ancestor whose hmax already covers hmaxcell:
 ! a node's hmax is never below its children's, so neither is any of its ancestors'.
 ! Other threads update the same ancestors, so both the test and the update are atomic
 do while (node(n)%parent /= 0)
    n = node(n)%parent
!$omp atomic read
    hmaxn = node(n)%hmax
    if (hmaxn >= hmaxcell) exit
!$omp atomic
    node(n)%hmax = max(node(n)%hmax, hmaxcell)
 enddo

end subroutine set_hmaxcell

!-----------------------------------------------------------------------
!+
!  get the distance from the centre of mass of a cell
!+
!-----------------------------------------------------------------------
subroutine get_distance_from_centre_of_mass(inode,xi,yi,zi,dx,dy,dz,xcen)
 integer, intent(in)  :: inode
 real,    intent(in)  :: xi,yi,zi
 real,    intent(out) :: dx,dy,dz
 real,    intent(in), optional :: xcen(3)

 if (present(xcen)) then
    dx = xi - xcen(1)
    dy = yi - xcen(2)
    dz = zi - xcen(3)
 else
    dx = xi - node(inode)%xcen(1)
    dy = yi - node(inode)%xcen(2)
    dz = zi - node(inode)%xcen(3)
 endif

end subroutine get_distance_from_centre_of_mass

!-----------------------------------------------------------------------
!+
!  build the tree
!+
!-----------------------------------------------------------------------
subroutine build_tree(npart,nactive,xyzh,vxyzu,for_apr)
 use io,           only:nprocs
 use kdtree,       only:maketree,maketreeglobal!,revtree
 use dim,          only:mpi,use_sinktree
 use part,         only:nptmass,xyzmh_ptmass,maxp
 use allocutils,   only:allocate_array
 integer, intent(inout) :: npart
 integer, intent(in)    :: nactive
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(in)    :: vxyzu(:,:)
 logical, intent(in), optional :: for_apr
 logical :: apr_tree

 apr_tree = .false.
 if (present(for_apr)) apr_tree = for_apr

 !
 ! the listneigh array is threadprivate, but if the thread numbers or ids are changed
 ! then the memory might be lost. So the following lines are a failsafe
 ! to ensure that the listneigh array is always allocated for each thread
 !
 !$omp parallel
 if (.not. allocated(listneigh)) call allocate_array('listneigh',listneigh,maxp)
 !$omp end parallel

 if (mpi .and. nprocs > 1) then
    if (use_sinktree) then
       call maketreeglobal(nodeglobal,node,nodemap,globallevel,refinelevels,xyzh,npart,cellatid,leaf_is_active,ncells,&
                           apr_tree,nptmass,xyzmh_ptmass)
    else
       call maketreeglobal(nodeglobal,node,nodemap,globallevel,refinelevels,xyzh,npart,cellatid,leaf_is_active,ncells,&
                           apr_tree)
    endif
 else
    if (use_sinktree) then
       call maketree(node,xyzh,npart,leaf_is_active,ncells,apr_tree,nptmass=nptmass,xyzmh_ptmass=xyzmh_ptmass)
    else
       ! use revtree for small numbers of active particles to avoid tree rebuild overhead
       ! threshold: use revtree if < 0.1% of total particles
       !if (npart > 0 .and. nactive < 0.001*npart) then
       !   call revtree(node,xyzh,leaf_is_active,ncells)
       !else
       call maketree(node,xyzh,npart,leaf_is_active,ncells,apr_tree)
       !endif
    endif
 endif
 call list_active_leaves()

end subroutine build_tree

!-----------------------------------------------------------------------
!+
!  list the cells with leaf_is_active > 0, so that loops over cells
!  (density, force) need not visit every node of the tree: with
!  individual timesteps only a few leaves may be active
!
!  In parallel, in cell order: the active leaves are counted per chunk of
!  cells, and a running total gives each chunk its place in the list.
!+
!-----------------------------------------------------------------------
subroutine list_active_leaves()
!$ use omp_lib, only:omp_get_max_threads
 integer, allocatable :: nlist(:)
 integer :: icell,ic,nchunk,n,k

 n = int(ncells)
 nchunk = 1
!$ nchunk = omp_get_max_threads()
 allocate(nlist(0:nchunk))
 nlist = 0

 !$omp parallel default(none) &
 !$omp shared(n,leaf_is_active,nchunk,nlist,active_leaves,nactive_leaves) &
 !$omp private(icell,ic,k)
 !$omp do schedule(static)
 do ic=1,nchunk
    k = 0
    do icell=int((int(ic-1,8)*n)/nchunk)+1,int((int(ic,8)*n)/nchunk)
       if (leaf_is_active(icell) > 0) k = k + 1
    enddo
    nlist(ic) = k
 enddo
 !$omp end do
 !$omp single
 do ic=1,nchunk
    nlist(ic) = nlist(ic) + nlist(ic-1)
 enddo
 nactive_leaves = nlist(nchunk)
 !$omp end single
 !$omp do schedule(static)
 do ic=1,nchunk
    k = nlist(ic-1)
    do icell=int((int(ic-1,8)*n)/nchunk)+1,int((int(ic,8)*n)/nchunk)
       if (leaf_is_active(icell) > 0) then
          k = k + 1
          active_leaves(k) = icell
       endif
    enddo
 enddo
 !$omp end do
 !$omp end parallel
 deallocate(nlist)

end subroutine list_active_leaves

!-----------------------------------------------------------------------
!+
! Using the k-d tree, compiles the neighbour list for the
! current cell (this list is common to all particles in the cell)
!
! the list is returned in 'listneigh' (length nneigh)
!+
!-----------------------------------------------------------------------
subroutine get_neighbour_list(inode,mylistneigh,nneigh,xyzh,xyzcache,ixyzcachesize, &
                              getj,f,remote_export,cell_xpos,cell_xsizei,cell_rcuti)
 use io,       only:nprocs,warning
 use dim,      only:mpi
 use kdtree,   only:getneigh,getneigh_dual,lenfgrav
 use kernel,   only:radkern
 use part,     only:gravity,periodic
 use boundary, only:dxbound,dybound,dzbound
 integer, intent(in)  :: inode,ixyzcachesize
 integer, intent(out) :: mylistneigh(:)
 integer, intent(out) :: nneigh
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(out) :: xyzcache(:,:)
 logical, intent(in),  optional :: getj
 real,    intent(out), optional :: f(lenfgrav)
 logical, intent(out), optional :: remote_export(:)
 real,    intent(in),  optional :: cell_xpos(3),cell_xsizei,cell_rcuti
 real :: xpos(3)
 real :: fgrav(lenfgrav),fgrav_global(lenfgrav)
 real :: xsizei,rcuti
 logical :: get_j,global_search,get_f
!
!--retrieve geometric centre of the node and the search radius (e.g. 2*hmax)
!
 if (present(cell_xpos)) then
    xpos = cell_xpos
    xsizei = cell_xsizei
    rcuti = cell_rcuti
 else
    call get_cell_location(inode,xpos,xsizei,rcuti)
 endif

 if (present(remote_export)) then
    if (nprocs > 1) global_search = .true.
    remote_export = .false.
 else
    global_search = .false.
 endif

 if (periodic) then
    if (rcuti > 0.5*min(dxbound,dybound,dzbound)) then
       call warning('get_neighbour_list', '2h > 0.5*L in periodic neighb. '//&
                'search: USE HIGHER RES, BIGGER BOX or LOWER MINPART IN TREE')
    endif
 endif
 !
 !--perform top-down tree walk to find all particles within radkern*h
 !  and force due to node-node interactions
 !
 get_j = .false.
 if (present(getj)) get_j = getj

 get_f = (gravity .and. present(f))

 if (mpi .and. nprocs > 1 .and. present(f)) then
    ! force with MPI: leaves walked in rounds from the global walk (local and ghost neighbours)
    call get_leaf_walk(inode,mylistneigh,nneigh,xyzcache,ixyzcachesize,f)
    return
 endif

 if (mpi .and. global_search) then ! no sym fmm for now...
    ! Find MPI tasks that have neighbours of this cell, output to remote_export
    call getneigh(nodeglobal,xpos,xsizei,rcuti,mylistneigh,nneigh,xyzcache,ixyzcachesize,&
                  cellatid,get_j,get_f,fgrav_global,remote_export)
 elseif (get_f) then
    ! Set fgrav to zero, which matters if gravity is enabled but global search is not
    fgrav_global = 0.0
 endif

 ! Find neighbours of this cell on this node
 if (get_f .and. .not.(mpi) .and. use_dualtree) then
    call getneigh_dual(node,xpos,xsizei,rcuti,mylistneigh,nneigh,xyzcache,ixyzcachesize,&
                          leaf_is_active,get_j,get_f,fgrav,inode)
 else
    call getneigh(node,xpos,xsizei,rcuti,mylistneigh,nneigh,xyzcache,ixyzcachesize,&
                     leaf_is_active,get_j,get_f,fgrav)
 endif

 if (get_f) f = fgrav + fgrav_global

end subroutine get_neighbour_list

!-----------------------------------------------------------------------
!+
!  dual tree walk on the global (refined) tree, one refined leaf of this
!  task at a time (same walk as getneigh_dual), stopping at the frontier
!  between refined and local nodes. Returns, for each refined leaf, the
!  M2L expansion inherited from all its ancestors, and the pairs of
!  refined leaves that are not well separated: this is the restart state
!  of the walk on the local trees, with local (owner = id) or remote sources
!
!  fnode_leaf(:,j) : expansion at refined leaf j (local node 2**refinelevels+j-1)
!  pairs(:,i)      : (dst,src) global indices of the refined leaves, owner of src
!+
!-----------------------------------------------------------------------
subroutine get_global_pairs(fnode_leaf,npairs,pairs)
 use io,     only:id,nprocs,fatal
 use kdtree, only:getneigh_dual_global,reset_cachestate_global,lenfgrav
 real,    allocatable, intent(inout) :: fnode_leaf(:,:)
 integer,              intent(out)   :: npairs
 integer, allocatable, intent(inout) :: pairs(:,:)
 integer :: nleaves,nnodes,ifirstleaf,j,k,icell,ibase,nneigh
 logical :: allpairs

 if (iand(nprocs,nprocs-1) /= 0) call fatal('get_global_pairs','number of MPI tasks must be a power of 2')

 nleaves    = 2**refinelevels
 ifirstleaf = 2**(globallevel+refinelevels) + id*nleaves
 nnodes     = 2**(globallevel+refinelevels+1) - 1
 ! the node caches of the tree walk are shared with the local tree
 if (nnodes > ncellsmax+1) call fatal('get_global_pairs','global tree larger than ncellsmax')

 ! nleaves = 2**refinelevels is the number of refined leaves of each task, not of the
 ! whole global level: this task only walks its own block of refined leaves
 if (any(cellatid(ifirstleaf:ifirstleaf+nleaves-1) /= id+1)) &
    call fatal('get_global_pairs','refined leaves walked are not all owned by this task')


 if (allocated(fnode_leaf)) then
    if (size(fnode_leaf,2) < nleaves) deallocate(fnode_leaf)
 endif
 if (.not.allocated(fnode_leaf)) allocate(fnode_leaf(lenfgrav,nleaves))

 npairs   = max(8*nleaves,512)
 allpairs = .false.
 ! each leaf reserves its slots in the pair list with one atomic update. The list
 ! cannot grow inside the loop: if it is too small, the walk is redone once it has
 ! been resized to the number of pairs found (the list is kept between calls)
 do while(.not.allpairs)

    if (allocated(pairs)) deallocate(pairs)
    allocate(pairs(3,2*npairs))

    call reset_cachestate_global(nnodes)
    npairs = 0

    !$omp parallel do default(none) schedule(dynamic) &
    !$omp shared(nleaves,ifirstleaf,nodeglobal,cellatid,fnode_leaf,npairs,pairs) &
    !$omp private(j,k,icell,ibase,nneigh)
    do j=1,nleaves
       icell = ifirstleaf + j - 1
       call getneigh_dual_global(nodeglobal,cellatid,icell,listneigh,nneigh,fnode_leaf(:,j))

       !$omp atomic capture
       ibase  = npairs
       npairs = npairs + nneigh
       !$omp end atomic

       if (ibase + nneigh <= size(pairs,2)) then
          do k=1,nneigh
             pairs(1,ibase+k) = icell
             pairs(2,ibase+k) = listneigh(k)
             pairs(3,ibase+k) = cellatid(listneigh(k)) - 1
          enddo
       endif
    enddo
    !$omp end parallel do
    allpairs = (npairs <= size(pairs,2))
 enddo

end subroutine get_global_pairs

!-----------------------------------------------------------------------
!+
!  walk of the local tree from the pairs of the global walk: for each
!  refined leaf (dst), the local src as local nodes, and the remote src
!  as (owner, node of the owner), then the walk in rounds of kslab_mpi
!+
!-----------------------------------------------------------------------
subroutine start_local_rounds(fnode_leaf,npairs,pairs)
 use io,     only:id
 use kdtree, only:global_to_local
 real,    intent(in) :: fnode_leaf(:,:)
 integer, intent(in) :: npairs
 integer, intent(in) :: pairs(:,:)
 integer, allocatable :: roots(:),istart(:),icount(:),srclist(:),rstart(:),rcount(:),srcrem(:,:)
 integer :: nleaves,ifirstleaf,i,j,irank,ksrc,nrounds

 nleaves    = 2**refinelevels
 ifirstleaf = 2**(globallevel+refinelevels) + id*nleaves

 allocate(roots(nleaves),istart(nleaves),icount(nleaves),rstart(nleaves),rcount(nleaves))
 allocate(srclist(max(npairs,1)),srcrem(2,max(npairs,1)))
 icount = 0
 rcount = 0
 do i=1,npairs
    j = pairs(1,i) - ifirstleaf + 1
    if (pairs(3,i) == id) then
       icount(j) = icount(j) + 1
    else
       rcount(j) = rcount(j) + 1
    endif
 enddo
 istart(1) = 0
 rstart(1) = 0
 do j=2,nleaves
    istart(j) = istart(j-1) + icount(j-1)
    rstart(j) = rstart(j-1) + rcount(j-1)
 enddo
 icount = 0
 rcount = 0
 do i=1,npairs
    j = pairs(1,i) - ifirstleaf + 1
    call global_to_local(pairs(2,i),globallevel,irank,ksrc)
    if (irank == id) then
       icount(j) = icount(j) + 1
       srclist(istart(j)+icount(j)) = ksrc
    else
       rcount(j) = rcount(j) + 1
       srcrem(:,rstart(j)+rcount(j)) = (/irank,ksrc/)
    endif
 enddo
 ! refined leaf j is the local node nleaves+j-1
 do j=1,nleaves
    roots(j) = nleaves + j - 1
 enddo

 call dualwalk_rounds(kslab_mpi,nleaves,roots,istart,icount,srclist,rstart,rcount,srcrem,fnode_leaf,nrounds)

end subroutine start_local_rounds

!-----------------------------------------------------------------------
!+
!  dual tree walk of the local tree in rounds, from a list of fake roots,
!  each with its expansion and the src nodes it still has to open, local
!  (srclist) or remote ((owner,node) in srcrem).
!  Before each round, each task sends the subtree (kslab levels) of its
!  fake roots with remote src to their owners, and receives the subtrees
!  of its remote src (as the pairs are mirrored, no request is needed).
!  In each round, the cells are the nodes kslab levels below each fake
!  root (or the leaves above): each cell is walked from its fake root,
!  over its local src then its remote src, and stops there. Its output
!  makes the fake roots of the next round. The output of the leaves (src
!  leaves + expansion) is kept for get_leaf_walk; a leaf with remote src
!  left to open (truncated) is also a fake root of the next round
!+
!-----------------------------------------------------------------------
subroutine dualwalk_rounds(kslab,nroots_in,roots_in,istart_in,icount_in,srclist_in,&
                           rstart_in,rcount_in,srcrem_in,fnode_in,nrounds)
 use io,       only:fatal
 use mpiutils, only:reduceall_mpi
 use mpiforce, only:check_pair_mirror
 use kdtree,   only:getneigh_dual_from,lenfgrav,kdnode
 integer, intent(in)  :: kslab,nroots_in
 integer, intent(in)  :: roots_in(:),istart_in(:),icount_in(:),srclist_in(:)
 integer, intent(in)  :: rstart_in(:),rcount_in(:),srcrem_in(:,:)
 real,    intent(in)  :: fnode_in(:,:)
 integer, intent(out) :: nrounds
 integer, allocatable :: roots(:),istart(:),icount(:),srclist(:),rstart(:),rcount(:),srcrem(:,:)
 integer, allocatable :: roots_new(:),istart_new(:),icount_new(:),srclist_new(:)
 integer, allocatable :: rstart_new(:),rcount_new(:),srcrem_new(:,:)
 integer, allocatable :: cells(:),cellroot(:),rslot(:),rleaf(:),rowner(:),rid(:),lrem(:),lpend(:)
 integer, allocatable :: pairs_chk(:,:)
 real,    allocatable :: fnode_root(:,:),fnode_new(:,:)
 type(kdnode), allocatable :: rnode(:)
 integer :: nroots,ncell,nslab,nrnode,ir,ic,ibase,icell,nneigh,nrem,npend,islot,i0,j
 integer :: nnew,nnewl,nnewr,nrec0,nleafsrc0,nleafrem0,nmismatch,nchk
 logical :: allfit,anyroots
 real    :: fnode(lenfgrav),fnode2(lenfgrav),fzero(lenfgrav)

 if (kslab < 1) call fatal('dualwalk_rounds','slab depth must be >= 1')
 nslab = 2**kslab
 fzero = 0.

 nroots = nroots_in
 allocate(roots(nroots),istart(nroots),icount(nroots),rstart(nroots),rcount(nroots),fnode_root(lenfgrav,nroots))
 roots      = roots_in(1:nroots)
 istart     = istart_in(1:nroots)
 icount     = icount_in(1:nroots)
 rstart     = rstart_in(1:nroots)
 rcount     = rcount_in(1:nroots)
 fnode_root = fnode_in(:,1:nroots)
 allocate(srclist(max(size(srclist_in),1)),srcrem(2,max(size(srcrem_in,2),1)))
 srclist(1:size(srclist_in)) = srclist_in
 srcrem(:,1:size(srcrem_in,2)) = srcrem_in

 ! records of the leaves
 if (.not.allocated(leafslot)) allocate(leafslot(size(node)))
 leafslot(1:ncells) = 0
 if (.not.allocated(leafsrc)) then
    allocate(leafsrc(1024),leafsrc_start(1024),leafsrc_count(1024),leafrec_prev(1024),leafrec_cell(1024),&
             fnode_leafcell(lenfgrav,1024))
    allocate(leafrem(2,1024),leafrem_start(1024),leafrem_count(1024))
 endif
 nleafslots = 0
 nleafsrc   = 0
 nleafrem   = 0

 ! nodes of the remote src of the first round (rnode(rslot(i)) is the remote src i), and
 ! whether a task still has fake roots: the rounds go on as long as one has
 nrounds = 0
 call exchange_round(kslab,nroots,roots,rstart,rcount,srcrem,rnode,rleaf,rowner,rid,nrnode,rslot,anyroots)
 do while(anyroots)
    nrounds = nrounds + 1

    ! cells of the round: each fake root reserves 2**kslab slots for its fake leaves
    allocate(cells(max(nroots*nslab,1)),cellroot(max(nroots*nslab,1)))
    ncell = 0
    !$omp parallel do default(none) schedule(dynamic) &
    !$omp shared(nroots,roots,cells,cellroot,ncell,nslab,kslab) private(ir,ibase)
    do ir=1,nroots
       !$omp atomic capture
       ibase = ncell
       ncell = ncell + nslab
       !$omp end atomic
       call get_fake_leaves(roots(ir),kslab,cells(ibase+1:ibase+nslab))
       cellroot(ibase+1:ibase+nslab) = ir
    enddo
    !$omp end parallel do

    ! walk each cell from its fake root; output reserved with one atomic per cell. The
    ! lists cannot grow inside the loop: if they are too small, grow them and redo
    allocate(roots_new(max(ncell,1)),istart_new(max(ncell,1)),icount_new(max(ncell,1)),&
             rstart_new(max(ncell,1)),rcount_new(max(ncell,1)),fnode_new(lenfgrav,max(ncell,1)))
    if (.not.allocated(srclist_new)) allocate(srclist_new(max(size(srclist),1024)))
    if (.not.allocated(srcrem_new))  allocate(srcrem_new(2,max(size(srcrem,2),1024)))
    nrec0     = nleafslots
    nleafsrc0 = nleafsrc
    nleafrem0 = nleafrem
    allfit    = .false.
    do while(.not.allfit)
       nnew       = 0
       nnewl      = 0
       nnewr      = 0
       nleafslots = nrec0
       nleafsrc   = nleafsrc0
       nleafrem   = nleafrem0
       !$omp parallel default(none) &
       !$omp shared(ncell,cells,cellroot,roots,istart,icount,srclist,rstart,rcount,rslot,fnode_root,fzero) &
       !$omp shared(node,leaf_is_active,rnode,rleaf,rowner,rid,nrnode) &
       !$omp shared(nnew,nnewl,nnewr,roots_new,istart_new,icount_new,srclist_new,rstart_new,rcount_new,srcrem_new) &
       !$omp shared(fnode_new,nleafslots,nleafsrc,nleafrem,leafslot,leafrec_prev,leafrec_cell) &
       !$omp shared(leafsrc,leafsrc_start,leafsrc_count) &
       !$omp shared(leafrem,leafrem_start,leafrem_count,fnode_leafcell) &
       !$omp private(ic,ir,icell,nneigh,nrem,npend,fnode,fnode2,islot,i0,j,lrem,lpend)
       allocate(lrem(max(nrnode,1)),lpend(max(nrnode,1)))
       !$omp do schedule(dynamic)
       do ic=1,ncell
          icell = cells(ic)
          if (icell == 0) cycle
          ir = cellroot(ic)
          !-- local src, then remote src
          call getneigh_dual_from(node,leaf_is_active,node,leaf_is_active,.true.,roots(ir),&
                                  srclist(istart(ir)+1:istart(ir)+icount(ir)),fnode_root(:,ir),&
                                  icell,listneigh,nneigh,lpend,npend,fnode)
          if (npend /= 0) call fatal('dualwalk_rounds','truncated src in the local tree')
          nrem  = 0
          if (rcount(ir) > 0) then
             call getneigh_dual_from(node,leaf_is_active,rnode,rleaf,.false.,roots(ir),&
                                     rslot(rstart(ir)+1:rstart(ir)+rcount(ir)),fzero,&
                                     icell,lrem,nrem,lpend,npend,fnode2)
             fnode = fnode + fnode2
          endif

          if (leaf_is_active(icell) /= 0) then
             !-- leaf: record its src leaves and expansion
             !$omp atomic capture
             islot = nleafslots
             nleafslots = nleafslots + 1
             !$omp end atomic
             !$omp atomic capture
             i0 = nleafsrc
             nleafsrc = nleafsrc + nneigh
             !$omp end atomic
             if (islot < size(leafsrc_start) .and. i0+nneigh <= size(leafsrc)) then
                leafrec_prev(islot+1)  = leafslot(icell)
                leafrec_cell(islot+1)  = icell
                leafslot(icell)        = islot + 1
                leafsrc_start(islot+1) = i0
                leafsrc_count(islot+1) = nneigh
                fnode_leafcell(:,islot+1) = fnode
                leafsrc(i0+1:i0+nneigh) = listneigh(1:nneigh)
             endif
             !$omp atomic capture
             i0 = nleafrem
             nleafrem = nleafrem + nrem
             !$omp end atomic
             if (islot < size(leafrem_start) .and. i0+nrem <= size(leafrem,2)) then
                leafrem_start(islot+1) = i0
                leafrem_count(islot+1) = nrem
                do j=1,nrem
                   leafrem(:,i0+j) = (/rowner(lrem(j)),rid(lrem(j))/)
                enddo
             endif
             !-- remote src left to open: the leaf is a fake root of the next round
             nneigh = 0
             nrem   = npend
             lrem(1:npend) = lpend(1:npend)
          elseif (npend > 0) then
             call fatal('dualwalk_rounds','truncated src left to open by an internal node')
          endif

          if (leaf_is_active(icell) == 0 .or. nrem > 0) then
             !$omp atomic capture
             islot = nnew
             nnew = nnew + 1
             !$omp end atomic
             roots_new(islot+1)   = icell
             fnode_new(:,islot+1) = fnode
             !$omp atomic capture
             i0 = nnewl
             nnewl = nnewl + nneigh
             !$omp end atomic
             istart_new(islot+1) = i0
             icount_new(islot+1) = nneigh
             if (i0+nneigh <= size(srclist_new)) srclist_new(i0+1:i0+nneigh) = listneigh(1:nneigh)
             !$omp atomic capture
             i0 = nnewr
             nnewr = nnewr + nrem
             !$omp end atomic
             rstart_new(islot+1) = i0
             rcount_new(islot+1) = nrem
             if (i0+nrem <= size(srcrem_new,2)) then
                do j=1,nrem
                   srcrem_new(:,i0+j) = (/rowner(lrem(j)),rid(lrem(j))/)
                enddo
             endif
          endif
       enddo
       !$omp end do
       deallocate(lrem,lpend)
       !$omp end parallel

       allfit = (nnewl <= size(srclist_new) .and. nnewr <= size(srcrem_new,2) .and. &
                 nleafslots <= size(leafsrc_start) .and. nleafsrc <= size(leafsrc) .and. &
                 nleafrem <= size(leafrem,2))
       ! a leaf walked in this round goes back to its previous record before the redo
       ! (and before the record arrays are grown, keeping only the previous rounds)
       if (.not.allfit) then
          do ic=1,ncell
             if (cells(ic) > 0) then
                if (leafslot(cells(ic)) > nrec0) leafslot(cells(ic)) = leafrec_prev(leafslot(cells(ic)))
             endif
          enddo
       endif
       if (nnewl > size(srclist_new)) then
          deallocate(srclist_new)
          allocate(srclist_new(2*nnewl))
       endif
       if (nnewr > size(srcrem_new,2)) then
          deallocate(srcrem_new)
          allocate(srcrem_new(2,2*nnewr))
       endif
       if (nleafslots > size(leafsrc_start)) then
          call grow_int(leafsrc_start,2*nleafslots,nrec0)
          call grow_int(leafsrc_count,2*nleafslots,nrec0)
          call grow_int(leafrec_prev,2*nleafslots,nrec0)
          call grow_int(leafrec_cell,2*nleafslots,nrec0)
          call grow_int(leafrem_start,2*nleafslots,nrec0)
          call grow_int(leafrem_count,2*nleafslots,nrec0)
          call grow_real2(fnode_leafcell,2*nleafslots,nrec0)
       endif
       if (nleafsrc > size(leafsrc)) call grow_int(leafsrc,2*nleafsrc,nleafsrc0)
       if (nleafrem > size(leafrem,2)) call grow_int2(leafrem,2*nleafrem,nleafrem0)
    enddo

    ! the remote pairs left to open and the remote leaf-leaf pairs must be mirrored
    if (check_dualtree_mpi) then
    nchk = max(nnewr,nleafrem-nleafrem0,1)
    allocate(pairs_chk(3,nchk))
    do ir=1,nnew
       do j=1,rcount_new(ir)
          pairs_chk(:,rstart_new(ir)+j) = (/roots_new(ir),srcrem_new(2,rstart_new(ir)+j),srcrem_new(1,rstart_new(ir)+j)/)
       enddo
    enddo
    call check_pair_mirror(nnewr,pairs_chk,nmismatch)
    if (reduceall_mpi('+',nmismatch) > 0) call fatal('dualwalk_rounds','remote pairs left to open are not mirrored')
    nchk = 0
    do islot=nrec0+1,nleafslots
       do j=1,leafrem_count(islot)
          nchk = nchk + 1
          pairs_chk(:,nchk) = (/leafrec_cell(islot),leafrem(2,leafrem_start(islot)+j),leafrem(1,leafrem_start(islot)+j)/)
       enddo
    enddo
    call check_pair_mirror(nchk,pairs_chk,nmismatch)
    if (reduceall_mpi('+',nmismatch) > 0) call fatal('dualwalk_rounds','remote leaf-leaf pairs are not mirrored')
    deallocate(pairs_chk)
    endif
    ! the output of this round is the input of the next one
    nroots = nnew
    call move_alloc(roots_new,roots)
    call move_alloc(istart_new,istart)
    call move_alloc(icount_new,icount)
    call move_alloc(rstart_new,rstart)
    call move_alloc(rcount_new,rcount)
    call move_alloc(fnode_new,fnode_root)
    call move_alloc(srclist_new,srclist)
    call move_alloc(srcrem_new,srcrem)
    deallocate(cells,cellroot)

    ! nodes of the remote src of the next round
    call exchange_round(kslab,nroots,roots,rstart,rcount,srcrem,rnode,rleaf,rowner,rid,nrnode,rslot,anyroots)
 enddo

end subroutine dualwalk_rounds

!-----------------------------------------------------------------------
!+
!  exchange of the nodes for a round: the subtree (kslab levels) of the
!  fake roots with remote src goes to the owners of these src, and the
!  subtrees of the remote src are received, in increasing node order for
!  each pair of tasks (by symmetry, both sides know what the other needs).
!  rnode(rslot(i)) is the received node of the remote src i of srcrem
!+
!-----------------------------------------------------------------------
subroutine exchange_round(kslab,nroots,roots,rstart,rcount,srcrem,rnode,rleaf,rowner,rid,nrnode,rslot,anyroots)
 use io,        only:fatal,nprocs
 use mpiforce,  only:exchange_slabs
 use kdtree,    only:kdnode
 integer, intent(in) :: kslab,nroots
 integer, intent(in) :: roots(:),rstart(:),rcount(:),srcrem(:,:)
 type(kdnode), allocatable, intent(inout) :: rnode(:)
 integer,      allocatable, intent(inout) :: rleaf(:),rowner(:),rid(:),rslot(:)
 integer,                   intent(out)   :: nrnode
 logical,                   intent(out)   :: anyroots
 integer(kind=8), allocatable :: keysend(:),keyreq(:)
 integer, allocatable :: indx(:),reqslot(:)
 real,    allocatable :: sendbuf(:),recvbuf(:)
 integer :: nremtot,nsendkey,nreq,ir,j,i,n,irank,ipos,nfield,nkd,nn,ireq
 integer :: nsend(nprocs),nrecv(nprocs),slab(2**(kslab+1)),slabflag(2**(kslab+1)),slabchild(2**(kslab+1))
 type(kdnode) :: dummy

 nkd    = size(transfer(node(1),(/0./)))
 nfield = nkd + 5
 ! the lists of the fake roots are reserved in any order: the total is the end of the last one
 nremtot = 0
 if (nroots > 0) nremtot = maxval(rstart(1:nroots) + rcount(1:nroots))

 ! subtrees to send: (owner of the remote src, fake root), and to receive: (owner, remote src)
 allocate(keysend(max(nremtot,1)),keyreq(max(nremtot,1)),indx(max(nremtot,1)))
 do ir=1,nroots
    do j=rstart(ir)+1,rstart(ir)+rcount(ir)
       keysend(j) = int(srcrem(1,j),8)*2_8**32 + int(roots(ir),8)
       keyreq(j)  = int(srcrem(1,j),8)*2_8**32 + int(srcrem(2,j),8)
    enddo
 enddo
 call sort_unique(nremtot,keysend,nsendkey,indx)
 call sort_unique(nremtot,keyreq,nreq,indx)

 ! pack: the slabs for each task are contiguous, in increasing (task,node) order
 nsend = 0
 do i=1,nsendkey
    irank = int(keysend(i)/2_8**32)
    nsend(irank+1) = nsend(irank+1) + count_slab(int(mod(keysend(i),2_8**32)),kslab)*nfield
 enddo
 allocate(sendbuf(max(sum(nsend),1)))
 ipos = 0
 do i=1,nsendkey
    call get_slab(int(mod(keysend(i),2_8**32)),kslab,slab,slabflag,slabchild,nn)
    call pack_slab(slab,slabflag,slabchild,nn,nkd,sendbuf(ipos+1:ipos+nn*nfield))
    ipos = ipos + nn*nfield
 enddo

 call exchange_slabs(sendbuf,nsend,recvbuf,nrecv,(nroots > 0),anyroots)

 ! unpack: the slabs from each task come in the order of the requested src
 nrnode = sum(nrecv)/nfield
 if (allocated(rnode)) deallocate(rnode,rleaf,rowner,rid)
 allocate(rnode(max(nrnode,1)),rleaf(max(nrnode,1)),rowner(max(nrnode,1)),rid(max(nrnode,1)),reqslot(max(nreq,1)))
 n     = 0
 ireq  = 0
 ipos  = 0
 do irank=0,nprocs-1
    do while(ipos < sum(nrecv(1:irank+1)))
       ireq = ireq + 1
       nn = nint(recvbuf(ipos+nkd+5))
       if (ireq > nreq) call fatal('exchange_round','more subtrees received than requested')
       if (keyreq(ireq) /= int(irank,8)*2_8**32 + nint(recvbuf(ipos+nkd+1),8)) &
          call fatal('exchange_round','subtree received is not the one requested')
       reqslot(ireq) = n + 1
       do j=1,nn
          rnode(n+j)  = transfer(recvbuf(ipos+1:ipos+nkd),dummy)
          rid(n+j)    = nint(recvbuf(ipos+nkd+1))
          rnode(n+j)%leftchild  = slot_of(nint(recvbuf(ipos+nkd+2)),n)
          rnode(n+j)%rightchild = slot_of(nint(recvbuf(ipos+nkd+3)),n)
          rleaf(n+j)  = nint(recvbuf(ipos+nkd+4))
          rowner(n+j) = irank
          ipos = ipos + nfield
       enddo
       n = n + nn
    enddo
 enddo
 if (ireq /= nreq) call fatal('exchange_round','fewer subtrees received than requested')

 ! slot of each remote src
 if (allocated(rslot)) deallocate(rslot)
 allocate(rslot(max(nremtot,1)))
 do j=1,nremtot
    rslot(j) = reqslot(find_key(nreq,keyreq,int(srcrem(1,j),8)*2_8**32 + int(srcrem(2,j),8)))
 enddo

contains

integer function slot_of(irel,nbase)
 integer, intent(in) :: irel,nbase

 slot_of = 0
 if (irel > 0) slot_of = nbase + irel

end function slot_of

end subroutine exchange_round

!-----------------------------------------------------------------------
!+
!  nodes of the subtree of iroot down to kslab levels below, in breadth
!  first order (slab(1) = iroot), with their leaf flag (isrc_truncated
!  for the internal nodes at the bottom, whose children are not in it)
!  and the position of their left child in the slab (0 if not in it)
!+
!-----------------------------------------------------------------------
subroutine get_slab(iroot,kslab,slab,slabflag,slabchild,nn)
 use kdtree, only:isrc_truncated
 integer, intent(in)  :: iroot,kslab
 integer, intent(out) :: slab(:),slabflag(:),slabchild(:),nn
 integer :: i,n

 nn = 1
 slab(1) = iroot
 i = 0
 do while(i < nn)
    i = i + 1
    n = slab(i)
    slabchild(i) = 0
    if (node(n)%leftchild == 0 .or. leaf_is_active(n) /= 0) then
       slabflag(i) = leaf_is_active(n)
    elseif (node(n)%level - node(iroot)%level >= kslab) then
       slabflag(i) = isrc_truncated
    else
       slabflag(i)  = 0
       slabchild(i) = nn + 1
       slab(nn+1) = node(n)%leftchild
       slab(nn+2) = node(n)%rightchild
       nn = nn + 2
    endif
 enddo

end subroutine get_slab

!-----------------------------------------------------------------------
!+
!  number of nodes of the slab of iroot
!+
!-----------------------------------------------------------------------
integer function count_slab(iroot,kslab)
 integer, intent(in) :: iroot,kslab
 integer :: slab(2**(kslab+1)),slabflag(2**(kslab+1)),slabchild(2**(kslab+1))

 call get_slab(iroot,kslab,slab,slabflag,slabchild,count_slab)

end function count_slab

!-----------------------------------------------------------------------
!+
!  pack the nodes of a slab: the node, its id, its children (position
!  in the slab, 0 if not in it), its leaf flag and the number of nodes
!  of the slab
!+
!-----------------------------------------------------------------------
subroutine pack_slab(slab,slabflag,slabchild,nn,nkd,buf)
 integer, intent(in)  :: slab(:),slabflag(:),slabchild(:),nn,nkd
 real,    intent(out) :: buf(:)
 integer :: i,ipos,nfield

 nfield = nkd + 5
 do i=1,nn
    ipos = (i-1)*nfield
    buf(ipos+1:ipos+nkd) = transfer(node(slab(i)),buf(1:nkd))
    buf(ipos+nkd+1) = real(slab(i))
    if (slabchild(i) > 0) then
       buf(ipos+nkd+2) = real(slabchild(i))
       buf(ipos+nkd+3) = real(slabchild(i)+1)
    else
       buf(ipos+nkd+2:ipos+nkd+3) = 0.
    endif
    buf(ipos+nkd+4) = real(slabflag(i))
    buf(ipos+nkd+5) = real(nn)
 enddo

end subroutine pack_slab

!-----------------------------------------------------------------------
!+
!  sort keys and remove the duplicates
!+
!-----------------------------------------------------------------------
subroutine sort_unique(n,key,nunique,indx)
 use sortutils, only:indexx
 integer,         intent(in)    :: n
 integer(kind=8), intent(inout) :: key(:)
 integer,         intent(out)   :: nunique
 integer,         intent(inout) :: indx(:)
 integer(kind=8), allocatable :: tmp(:)
 integer :: i

 nunique = 0
 if (n == 0) return
 call indexx(n,key,indx)
 allocate(tmp(n))
 tmp = key(indx(1:n))
 do i=1,n
    if (nunique > 0) then
       if (tmp(i) == key(nunique)) cycle
    endif
    nunique = nunique + 1
    key(nunique) = tmp(i)
 enddo

end subroutine sort_unique

!-----------------------------------------------------------------------
!+
!  position of a key in a sorted list (binary search)
!+
!-----------------------------------------------------------------------
integer function find_key(n,key,k)
 use io, only:fatal
 integer,         intent(in) :: n
 integer(kind=8), intent(in) :: key(:),k
 integer :: ilo,ihi,imid

 ilo = 1
 ihi = n
 find_key = 0
 do while(ilo <= ihi)
    imid = (ilo+ihi)/2
    if (key(imid) == k) then
       find_key = imid
       return
    elseif (key(imid) < k) then
       ilo = imid + 1
    else
       ihi = imid - 1
    endif
 enddo
 call fatal('find_key','remote src not received')

end function find_key

!-----------------------------------------------------------------------
!+
!  fake leaves of a fake root: its descendants kslab levels below, or the
!  active leaves above (unused slots are set to 0)
!+
!-----------------------------------------------------------------------
subroutine get_fake_leaves(iroot,kslab,cells)
 integer, intent(in)  :: iroot,kslab
 integer, intent(out) :: cells(:)
 integer :: stack(128),istack,n,ideproot,idep,ncell

 ideproot = node(iroot)%level
 cells  = 0
 ncell  = 0
 istack = 1
 stack(1) = iroot
 do while(istack > 0)
    n      = stack(istack)
    istack = istack - 1
    idep   = node(n)%level - ideproot
    if (node(n)%leftchild == 0 .or. leaf_is_active(n) /= 0) then
       if (leaf_is_active(n) > 0) then
          ncell = ncell + 1
          cells(ncell) = n
       endif
    elseif (idep == kslab) then
       ncell = ncell + 1
       cells(ncell) = n
    else
       stack(istack+1) = node(n)%leftchild
       stack(istack+2) = node(n)%rightchild
       istack = istack + 2
    endif
 enddo

end subroutine get_fake_leaves

!-----------------------------------------------------------------------
!+
!  neighbours and expansion of a leaf from the walk in rounds: the
!  particles of its src leaves, and its expansion
!+
!-----------------------------------------------------------------------
subroutine get_leaf_walk(inode,mylistneigh,nneigh,xyzcache,ixyzcachesize,fnode)
 use io,     only:fatal
 use kdtree, only:getneigh_srcleaves,lenfgrav
 integer, intent(in)  :: inode,ixyzcachesize
 integer, intent(out) :: mylistneigh(:)
 integer, intent(out) :: nneigh
 real,    intent(out) :: fnode(lenfgrav)
 real,    intent(out) :: xyzcache(:,:)
 integer, allocatable :: srcleaves(:)
 integer :: irec,nsrc,j,ighost

 irec = leafslot(inode)
 if (irec == 0) call fatal('get_leaf_walk','leaf not walked in rounds')
 ! the expansion of the last record is complete
 fnode = fnode_leafcell(:,irec)
 ! local src leaves of all the records of this leaf
 nsrc = 0
 do while(irec > 0)
    nsrc = nsrc + leafsrc_count(irec)
    irec = leafrec_prev(irec)
 enddo
 allocate(srcleaves(max(nsrc,1)))
 nsrc = 0
 irec = leafslot(inode)
 do while(irec > 0)
    srcleaves(nsrc+1:nsrc+leafsrc_count(irec)) = leafsrc(leafsrc_start(irec)+1:leafsrc_start(irec)+leafsrc_count(irec))
    nsrc = nsrc + leafsrc_count(irec)
    irec = leafrec_prev(irec)
 enddo
 call getneigh_srcleaves(node,inode,srcleaves,nsrc,mylistneigh,nneigh,xyzcache,ixyzcachesize)

 ! then the ghost particles of the remote src leaves
 irec = leafslot(inode)
 do while(irec > 0)
    do j=leafrem_start(irec)+1,leafrem_start(irec)+leafrem_count(irec)
       ighost = find_ghost(int(leafrem(1,j),8)*2_8**32 + int(leafrem(2,j),8))
       call cache_ghosts(inode,ighost,mylistneigh,nneigh,xyzcache,ixyzcachesize)
    enddo
    irec = leafrec_prev(irec)
 enddo

end subroutine get_leaf_walk

!-----------------------------------------------------------------------
!+
!  remote leaves to exchange as ghost particles, from the remote leaf-leaf
!  pairs of the walk in rounds: to send, (owner of the remote leaf, local
!  leaf), and to receive, (owner, remote leaf), sorted (as for the nodes,
!  the pairs being mirrored, both sides agree on what to exchange)
!+
!-----------------------------------------------------------------------
subroutine get_remote_leaves(nsend,keysend,nreq,keyreq)
 integer,                      intent(out) :: nsend,nreq
 integer(kind=8), allocatable, intent(out) :: keysend(:),keyreq(:)
 integer, allocatable :: indx(:)
 integer :: irec,j,n

 allocate(keysend(max(nleafrem,1)),keyreq(max(nleafrem,1)),indx(max(nleafrem,1)))
 n = 0
 do irec=1,nleafslots
    do j=leafrem_start(irec)+1,leafrem_start(irec)+leafrem_count(irec)
       n = n + 1
       keysend(n) = int(leafrem(1,j),8)*2_8**32 + int(leafrec_cell(irec),8)
       keyreq(n)  = int(leafrem(1,j),8)*2_8**32 + int(leafrem(2,j),8)
    enddo
 enddo
 call sort_unique(n,keysend,nsend,indx)
 call sort_unique(n,keyreq,nreq,indx)

end subroutine get_remote_leaves

!-----------------------------------------------------------------------
!+
!  store the table of the ghost leaves: for the remote leaf keyreq(i), its
!  first ghost particle, its number of particles and its centre, and the
!  mass of each ghost particle (ghost k is particle ibase+k)
!+
!-----------------------------------------------------------------------
subroutine set_ghost_leaves(nreq,keyreq,ifirst,icount,xcen,ibase,nghost,mass)
 integer,         intent(in) :: nreq,ibase,nghost
 integer(kind=8), intent(in) :: keyreq(:)
 integer,         intent(in) :: ifirst(:),icount(:)
 real,            intent(in) :: xcen(:,:),mass(:)

 if (allocated(ghostkey)) deallocate(ghostkey,ghostfirst,ghostcount,ghostxcen)
 if (allocated(ghostmass)) deallocate(ghostmass)
 allocate(ghostkey(max(nreq,1)),ghostfirst(max(nreq,1)),ghostcount(max(nreq,1)),ghostxcen(3,max(nreq,1)))
 allocate(ghostmass(max(nghost,1)))
 nghostleaves = nreq
 ighostbase   = ibase
 ghostkey(1:nreq)    = keyreq(1:nreq)
 ghostfirst(1:nreq)  = ifirst(1:nreq)
 ghostcount(1:nreq)  = icount(1:nreq)
 ghostxcen(:,1:nreq) = xcen(:,1:nreq)
 ghostmass(1:nghost) = mass(1:nghost)

end subroutine set_ghost_leaves

!-----------------------------------------------------------------------
!+
!  position of a remote leaf in the table of the ghost leaves
!+
!-----------------------------------------------------------------------
integer function find_ghost(k)
 integer(kind=8), intent(in) :: k

 find_ghost = find_key(nghostleaves,ghostkey,k)

end function find_ghost

!-----------------------------------------------------------------------
!+
!  add the ghost particles of a ghost leaf to the neighbours of inode,
!  as cache_neighbours does for a local leaf
!+
!-----------------------------------------------------------------------
subroutine cache_ghosts(inode,ighost,mylistneigh,nneigh,xyzcache,ixyzcachesize)
 use part,     only:xyzh,rho,gradh,periodic
 use dim,      only:igradomega,igradzeta,igradsoft
 use kdtree,   only:ih1,im,irho,izetaomega,isoftomega
 use boundary, only:dxbound,dybound,dzbound
 integer, intent(in)    :: inode,ighost,ixyzcachesize
 integer, intent(inout) :: mylistneigh(:),nneigh
 real,    intent(inout) :: xyzcache(:,:)
 integer :: k,ip,maxcache
 real    :: offset(3),dx(3)

 maxcache = 0
 if (ixyzcachesize > 0) maxcache = size(xyzcache,1)
 ! periodic offset between the two leaves, as get_sep
 offset = 0.
 if (periodic) then
    dx = ghostxcen(:,ighost) - node(inode)%xcen
    if (abs(dx(1)) > 0.5*dxbound) offset(1) = -dxbound*sign(1.0,dx(1))
    if (abs(dx(2)) > 0.5*dybound) offset(2) = -dybound*sign(1.0,dx(2))
    if (abs(dx(3)) > 0.5*dzbound) offset(3) = -dzbound*sign(1.0,dx(3))
 endif
 do k=1,ghostcount(ighost)
    ip = ghostfirst(ighost) + k - 1
    nneigh = nneigh + 1
    mylistneigh(nneigh) = ip
    if (nneigh <= ixyzcachesize) then
       xyzcache(1:3,nneigh) = xyzh(1:3,ip) + offset
       if (maxcache >= 4) xyzcache(ih1,nneigh) = 1./xyzh(4,ip)
       if (maxcache >= 5) xyzcache(im,nneigh)  = ghostmass(ip - ighostbase)
       if (maxcache >= 7) then
          xyzcache(irho,nneigh)       = rho(ip)
          xyzcache(izetaomega,nneigh) = real(gradh(igradzeta,ip))*real(gradh(igradomega,ip))
       endif
       if (maxcache >= 8) then
          if (size(gradh,1) >= igradsoft) then
             xyzcache(isoftomega,nneigh) = real(gradh(igradsoft,ip))*real(gradh(igradomega,ip))
          else
             xyzcache(isoftomega,nneigh) = 0.
          endif
       endif
    endif
 enddo


end subroutine cache_ghosts

!-----------------------------------------------------------------------
!+
!  grow an array, keeping its first nkeep elements
!+
!-----------------------------------------------------------------------
subroutine grow_int(a,n,nkeep)
 integer, allocatable, intent(inout) :: a(:)
 integer,              intent(in)    :: n,nkeep
 integer, allocatable :: tmp(:)

 allocate(tmp(n))
 tmp(1:nkeep) = a(1:nkeep)
 call move_alloc(tmp,a)

end subroutine grow_int

subroutine grow_int2(a,n,nkeep)
 integer, allocatable, intent(inout) :: a(:,:)
 integer,              intent(in)    :: n,nkeep
 integer, allocatable :: tmp(:,:)

 allocate(tmp(size(a,1),n))
 tmp(:,1:nkeep) = a(:,1:nkeep)
 call move_alloc(tmp,a)

end subroutine grow_int2

subroutine grow_real2(a,n,nkeep)
 real, allocatable, intent(inout) :: a(:,:)
 integer,           intent(in)    :: n,nkeep
 real, allocatable :: tmp(:,:)

 allocate(tmp(size(a,1),n))
 tmp(:,1:nkeep) = a(:,1:nkeep)
 call move_alloc(tmp,a)

end subroutine grow_real2

!-----------------------------------------------------------------------
!+
!  get neighbours around an arbitrary position in space
!+
!-----------------------------------------------------------------------
subroutine getneigh_pos(xpos,xsizei,rcuti,mylistneigh,nneigh,xyzcache,ixyzcachesize,leaf_is_active,get_j)
 use kdtree, only:getneigh
 integer, intent(in)  :: ixyzcachesize
 real,    intent(in)  :: xpos(3)
 real,    intent(in)  :: xsizei,rcuti
 integer, intent(out) :: mylistneigh(:)
 integer, intent(out) :: nneigh
 real,    intent(out) :: xyzcache(:,:)
 integer, intent(in)  :: leaf_is_active(:) !ncellsmax+1)
 logical, intent(in), optional :: get_j
 logical :: getj

 getj = .false.
 if (present(get_j)) getj=get_j
 call getneigh(node,xpos,xsizei,rcuti,mylistneigh,nneigh,xyzcache,ixyzcachesize, &
               leaf_is_active,getj,.false.)

end subroutine getneigh_pos

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_tree(iunit)
 use kdtree,       only:tree_accuracy
 use infile_utils, only:write_inopt
 use part,         only:gravity
 integer, intent(in) :: iunit

 if (gravity) call write_inopt(tree_accuracy,'tree_accuracy','tree opening criterion (0.0-1.0)',iunit)

end subroutine write_options_tree

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_tree(db,nerr)
 use part,         only:gravity
 use kdtree,       only:tree_accuracy
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 if (gravity) call read_inopt(tree_accuracy,'tree_accuracy',db,errcount=nerr,min=0.,max=1.)

end subroutine read_options_tree

!-----------------------------------------------------------------------
!+
!  find the position and size of a tree node
!+
!-----------------------------------------------------------------------
subroutine get_cell_location(inode,xpos,xsizei,rcuti)
 use kernel, only:radkern
 integer, intent(in)  :: inode
 real,    intent(out) :: xpos(3)
 real,    intent(out) :: xsizei
 real,    intent(out) :: rcuti

 xpos    = node(inode)%xcen(1:3)
 xsizei  = node(inode)%size
 rcuti   = radkern*node(inode)%hmax

end subroutine get_cell_location

!-----------------------------------------------------------------------
!+
!  sync the hmax values across all MPI tasks
!+
!-----------------------------------------------------------------------
subroutine sync_hmax_mpi
 use mpiutils,  only:reduceall_mpi
 use io,        only:nprocs
 integer :: i, n
 real    :: hmax(2**(globallevel+refinelevels+1)-1)

 hmax(:) = 0.0
 ! copy hmax values into contiguous array
 do i = 2,2**(refinelevels+1)-1
    hmax(nodemap(i)) = node(i)%hmax
 enddo

 ! reduce across threads
 hmax = reduceall_mpi('max', hmax)

 ! put values back into node
 do i = 2*nprocs,2**(globallevel+refinelevels+1)-1
    nodeglobal(i)%hmax = hmax(i)
 enddo

 ! walk tree up
 do i = 2*nprocs,4*nprocs
    n = i
    do while (nodeglobal(n)%parent /= 0)
       n = nodeglobal(n)%parent
       nodeglobal(n)%hmax = max(nodeglobal(n)%hmax, hmax(i))
    enddo
 enddo

end subroutine sync_hmax_mpi

end module neighkdtree
