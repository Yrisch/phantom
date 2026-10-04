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
!   mpibalance, mpidomain, mpitree, mpiutils, part, timing
!
 use dim,         only:maxp,ncellsmax,minpart,use_apr,use_sinktree,maxptmass,maxpsph,gravity
 use io,          only:nprocs
 use dtypekdtree, only:kdnode,lenfgrav
 use part,        only:ll,iphase,treecache,maxphase, &
                       apr_level,aprmassoftype

 implicit none

 integer, public,  allocatable :: inoderange(:,:)
 integer, public,  allocatable :: inodeparts(:)
 type(kdnode),     allocatable :: refinementnode(:)
 integer,          allocatable :: cachestate(:)
 integer,          allocatable :: neighnodecount_branch(:)
 integer,          allocatable :: neighnode_branch(:,:)
 integer,          allocatable :: neighnodecache(:)
 integer,          allocatable :: neighnodecache_start(:)
 integer,          allocatable :: neighnodecache_count(:)
 real,             allocatable :: fnodecache(:,:)
 real,             allocatable :: fnode_branch(:,:)
!$omp threadprivate(fnode_branch,neighnode_branch,neighnodecount_branch)
!
!--tree parameters
!
 integer,          parameter, public :: irootnode    = 1
 character(len=1), parameter, public :: labelax(3)   = (/'x','y','z'/)
 integer,          parameter         :: maxdepth     = 64
 integer,          parameter         :: maxnodecache_local = 512
 integer,          parameter         :: maxneigh_per_node  = 16
 integer,          parameter         :: maxstacksize = 2048
!
!--runtime options for this module
!
 real,    public  :: tree_accuracy    = 0.5
 logical, public  :: use_geosplit     = gravity ! only debug flag / should be on gravity
 logical, public  :: use_cache        = .true.  ! only debug flag
 ! scratch space for the parallel partition in build_top_parallel
 real,    allocatable, private :: tcbuf(:,:)
 integer, allocatable, private :: ipbuf(:)

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
 interface
  module subroutine maketree(node, xyzh, np, leaf_is_active, ncells, apr_tree, refinelevels,nptmass,xyzmh_ptmass)
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

  end subroutine maketree
 end interface


!--------------------------------------------------------------------------------
!+
!  Routine to build the global level tree
!+
!-------------------------------------------------------------------------------
 interface
  module subroutine maketreeglobal(nodeglobal,node,nodemap,globallevel,refinelevels,xyzh,&
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
  end subroutine maketreeglobal
 end interface
!-----------------------------------------------
!+
!  Routine to update a constructed tree
!  Note: current version ONLY works if
!  tree is built to < maxlevel_indexed
!  That is, it relies on the 2^n style tree
!  indexing to sweep out each level
!+
!-----------------------------------------------
 interface
  module subroutine revtree(node, xyzh, leaf_is_active, ncells)
   use dim,  only:maxp,use_apr,ind_timesteps
   use part, only:maxphase,iphase,igas,massoftype,iamtype,aprmassoftype,&
                apr_level,iactive,treecache,isdead_or_accreted
   use io,   only:fatal
   type(kdnode),    intent(inout) :: node(:) !ncellsmax+1)
   real,            intent(in)    :: xyzh(:,:)
   integer,         intent(inout) :: leaf_is_active(:) !ncellsmax+1)
   integer(kind=8), intent(in)    :: ncells
  end subroutine revtree
 end interface


!----------------------------------------------------------------
!+
!  Routine to walk tree for neighbour search
!  (all particles within a given h_i and optionally within h_j)
!+
!----------------------------------------------------------------
 interface
  module subroutine getneigh(node,xpos,xsizei,rcuti,listneigh,nneigh,xyzcache,ixyzcachesize,leaf_is_active,&
                    get_hj,get_f,fnode,remote_export,nq)
   use io,       only:fatal,id
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
  end subroutine getneigh
 end interface

!----------------------------------------------------------------
!+
!  Routine to walk tree for neighbour search (SFMM version)
!  (all particles within a given h_i and optionally within h_j)
!  A dual tree walk is used to compute
!  every node-node interactions
!+
!----------------------------------------------------------------
 interface
  module subroutine getneigh_dual(node,xpos,xsizei,rcuti,listneigh,nneigh,xyzcache,ixyzcachesize,leaf_is_active,&
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
  end subroutine getneigh_dual
 end interface

!-----------------------------------------------------------
!+
!  Compute the Taylor expansion coeffs between the node
!  centres using the quadrupole moments (p=3) (Dehnen 2002)
!+
!-----------------------------------------------------------
 interface
  pure module subroutine compute_M2L(dx,dy,dz,dr1,q0,quads,fnode)
   real, intent(in)    :: dx,dy,dz,dr1,q0
   real, intent(in)    :: quads(6)
   real, intent(inout) :: fnode(lenfgrav)
  end subroutine compute_M2L
 end interface

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
 interface
  pure module subroutine expand_fgrav_in_taylor_series(fnode,dx,dy,dz,fxi,fyi,fzi,poti)
   real, intent(in)  :: fnode(lenfgrav)
   real, intent(in)  :: dx,dy,dz
   real, intent(out) :: fxi,fyi,fzi,poti
  end subroutine expand_fgrav_in_taylor_series
 end interface


 private

contains

subroutine allocate_kdtree
 use dim, only:mpi
 use allocutils, only:allocate_array

 call allocate_array('inoderange', inoderange, 2, ncellsmax+1)
 call allocate_array('inodeparts', inodeparts, maxp)
 if (mpi) call allocate_array('refinementnode', refinementnode, ncellsmax+1)
 call allocate_array('cachestate', cachestate, ncellsmax+1)
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
 if (mpi .and. allocated(refinementnode)) deallocate(refinementnode)
 if (allocated(cachestate)) deallocate(cachestate)
 if (allocated(fnodecache)) deallocate(fnodecache)
 if (allocated(neighnodecache)) deallocate(neighnodecache)
 if (allocated(neighnodecache_start)) deallocate(neighnodecache_start)
 if (allocated(neighnodecache_count)) deallocate(neighnodecache_count)
!$omp parallel
 if (allocated(neighnode_branch)) deallocate(neighnode_branch)
 if (allocated(neighnodecount_branch)) deallocate(neighnodecount_branch)
 if (allocated(fnode_branch)) deallocate(fnode_branch)
!$omp end parallel
 if (allocated(tcbuf)) deallocate(tcbuf,ipbuf)

end subroutine deallocate_kdtree

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


end module kdtree
