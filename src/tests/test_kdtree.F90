!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testkdtree
!
! This module performs unit tests of the kdtree module
!   The tests here are specific to the tree, some general
!   tests of neighbour finding are done in test_neigh
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, io, kdtree, kernel, mpidomain, neighkdtree, part,
!   testutils, timing, unifdis
!
 implicit none
 public :: test_kdtree

 private

contains
!-----------------------------------------------------------------------
!+
!   Unit tests of tree code
!+
!-----------------------------------------------------------------------
subroutine test_kdtree(ntests,npass)
 use dim,         only:maxp,periodic,ind_timesteps
 use io,          only:id,master,iverbose
 use neighkdtree, only:leaf_is_active,ncells,node
 use part,        only:npart,xyzh,hfact,massoftype,igas,maxphase,iphase,isetphase,iactive
 use kernel,      only:hfact_default
 use kdtree,      only:maketree,revtree,kdnode,empty_tree
 use unifdis,     only:set_unifdis
 use testutils,   only:checkvalbuf,checkvalbuf_end,update_test_scores,checkval
 use timing,      only:print_time,getused
 use mpidomain,   only:i_belong
 integer, intent(inout) :: ntests,npass
 logical :: test_revtree, test_dualrounds, test_all
 integer :: i,nfailed(22),nchecked(22),nfailed_leaf(1),nchecked_leaf(1),ierrmax_leaf(1)
 real    :: psep,tol,tol_octs,errmax(22)
 real(4) :: t2,t1,tmaketree
 type(kdnode), allocatable :: old_tree(:)
 integer, allocatable :: leaf_is_active_saved(:)

 test_all = .true.
 test_revtree = .true.
 test_dualrounds = .true.
 iverbose = 2

 if (id==master) write(*,"(a,/)") '--> TESTING KDTREE'

 if (test_revtree .or. test_all) then
    if (id==master) write(*,"(/,a)") '--> testing revtree routine'
    !
    ! set up a random particle distribution
    !
    psep = 1./100.
    hfact = hfact_default
    npart = 0
    call set_unifdis('random',id,master,-0.5,0.5,-0.5,0.5,-0.5,0.5,&
                     psep,hfact,npart,xyzh,periodic,mask=i_belong)
    massoftype(igas) = 1000./npart
    if (maxphase==maxp) iphase(:) = isetphase(igas,iactive=.true.)

    !
    ! call maketree to build the tree
    !
    call empty_tree(node)
    call cpu_time(t1)
    call maketree(node,xyzh,npart,leaf_is_active,ncells,apr_tree=.false.)
    call cpu_time(t2)
    call print_time(t2-t1,'maketree completed in')
    !
    ! now save the tree structure and leaf_is_active
    !
    allocate(old_tree(int(ncells)))
    old_tree(1:ncells) = node(1:ncells)
    allocate(leaf_is_active_saved(int(ncells)))
    leaf_is_active_saved(1:int(ncells)) = leaf_is_active(1:int(ncells))

    !
    ! erase all information in the existing tree except the structure
    !
    do i=1,int(ncells)
       node(i)%xcen(:) = 0.
       node(i)%size    = 0.
       node(i)%hmax    = 0.
#ifdef GRAVITY
       node(i)%mass    = 0.
       node(i)%quads(:)= 0.
       node(i)%octs(:) = 0.
#endif
       leaf_is_active(i) = 0
    enddo

    !
    ! call revtree to rebuild
    !
    tmaketree = t2-t1
    call cpu_time(t1)
    call revtree(node,xyzh,leaf_is_active,ncells)
    call cpu_time(t2)
    call print_time(t2-t1,'revtree completed in')
    if (id==master) print*,' ratio of revtree/maketree: ',(t2-t1)/tmaketree

    !
    ! check that the revised tree matches the tree built
    !
    nfailed(:)  = 0
    nchecked(:) = 0
    errmax(:)   = 0.
    tol = 2.e-11
    ! use larger tolerance for some octupole moments due to variation from openMP loop ordering
    tol_octs = 1.e-9
    do i=1,int(ncells)
       if (i > 1 .and. node(i)%parent == 0) cycle
       ! if (leaf_is_active(i) /= 0) then
       call checkvalbuf(node(i)%xcen(1),old_tree(i)%xcen(1),tol,'x0',nfailed(1),nchecked(1),errmax(1))
       call checkvalbuf(node(i)%xcen(2),old_tree(i)%xcen(2),tol,'y0',nfailed(2),nchecked(2),errmax(2))
       call checkvalbuf(node(i)%xcen(3),old_tree(i)%xcen(3),tol,'z0',nfailed(3),nchecked(3),errmax(3))
!       call checkvalbuf(node(i)%size,old_tree(i)%size,tol,'size',nfailed(4),nchecked(4),errmax(4))
       call checkvalbuf((node(i)%size + tol >= old_tree(i)%size),.true.,'size',nfailed(4),nchecked(4))
       call checkvalbuf(node(i)%hmax,old_tree(i)%hmax,tol,'hmax',nfailed(5),nchecked(5),errmax(5))
#ifdef GRAVITY
       call checkvalbuf(node(i)%mass,old_tree(i)%mass,tol,'mass',nfailed(6),nchecked(6),errmax(6))
       call checkvalbuf(node(i)%quads(1),old_tree(i)%quads(1),tol,'qxx',nfailed(7),nchecked(7),errmax(7))
       call checkvalbuf(node(i)%quads(2),old_tree(i)%quads(2),2.*tol,'qxy',nfailed(8),nchecked(8),errmax(8))
       call checkvalbuf(node(i)%quads(3),old_tree(i)%quads(3),tol,'qxz',nfailed(9),nchecked(9),errmax(9))
       call checkvalbuf(node(i)%quads(4),old_tree(i)%quads(4),tol,'qyy',nfailed(10),nchecked(10),errmax(10))
       call checkvalbuf(node(i)%quads(5),old_tree(i)%quads(5),tol,'qyz',nfailed(11),nchecked(11),errmax(11))
       call checkvalbuf(node(i)%quads(6),old_tree(i)%quads(6),tol,'qzz',nfailed(12),nchecked(12),errmax(12))
       call checkvalbuf(node(i)%octs(1),old_tree(i)%octs(1),tol_octs,'oxxx',nfailed(13),nchecked(13),errmax(13))
       call checkvalbuf(node(i)%octs(2),old_tree(i)%octs(2),tol,'oxxy',nfailed(14),nchecked(14),errmax(14))
       call checkvalbuf(node(i)%octs(3),old_tree(i)%octs(3),tol,'oxxz',nfailed(15),nchecked(15),errmax(15))
       call checkvalbuf(node(i)%octs(4),old_tree(i)%octs(4),tol_octs,'oxyy',nfailed(16),nchecked(16),errmax(16))
       call checkvalbuf(node(i)%octs(5),old_tree(i)%octs(5),tol,'oxyz',nfailed(17),nchecked(17),errmax(17))
       call checkvalbuf(node(i)%octs(6),old_tree(i)%octs(6),tol_octs,'oxzz',nfailed(18),nchecked(18),errmax(18))
       call checkvalbuf(node(i)%octs(7),old_tree(i)%octs(7),tol,'oyyy',nfailed(19),nchecked(19),errmax(19))
       call checkvalbuf(node(i)%octs(8),old_tree(i)%octs(8),tol,'oyyz',nfailed(20),nchecked(20),errmax(20))
       call checkvalbuf(node(i)%octs(9),old_tree(i)%octs(9),tol,'oyzz',nfailed(21),nchecked(21),errmax(21))
       call checkvalbuf(node(i)%octs(10),old_tree(i)%octs(10),tol,'ozzz',nfailed(22),nchecked(22),errmax(22))
#endif
       ! endif
    enddo
    call checkvalbuf_end('x0',nchecked(1),nfailed(1),errmax(1),tol)
    call checkvalbuf_end('y0',nchecked(2),nfailed(2),errmax(2),tol)
    call checkvalbuf_end('z0',nchecked(3),nfailed(3),errmax(3),tol)
    call checkvalbuf_end('size',nchecked(4),nfailed(4),errmax(4),tol)
    call checkvalbuf_end('hmax',nchecked(5),nfailed(5),errmax(5),tol)
#ifdef GRAVITY
    call checkvalbuf_end('mass',nchecked(6),nfailed(6),errmax(6),tol)
    call checkvalbuf_end('qxx',nchecked(7),nfailed(7),errmax(7),tol)
    call checkvalbuf_end('qxy',nchecked(8),nfailed(8),errmax(8),tol)
    call checkvalbuf_end('qxz',nchecked(9),nfailed(9),errmax(9),tol)
    call checkvalbuf_end('qyy',nchecked(10),nfailed(10),errmax(10),tol)
    call checkvalbuf_end('qyz',nchecked(11),nfailed(11),errmax(11),tol)
    call checkvalbuf_end('qzz',nchecked(12),nfailed(12),errmax(12),tol)
    call checkvalbuf_end('oxxx',nchecked(13),nfailed(13),errmax(13),tol_octs)
    call checkvalbuf_end('oxxy',nchecked(14),nfailed(14),errmax(14),tol)
    call checkvalbuf_end('oxxz',nchecked(15),nfailed(15),errmax(15),tol)
    call checkvalbuf_end('oxyy',nchecked(16),nfailed(16),errmax(16),tol_octs)
    call checkvalbuf_end('oxyz',nchecked(17),nfailed(17),errmax(17),tol)
    call checkvalbuf_end('oxzz',nchecked(18),nfailed(18),errmax(18),tol_octs)
    call checkvalbuf_end('oyyy',nchecked(19),nfailed(19),errmax(19),tol)
    call checkvalbuf_end('oyyz',nchecked(20),nfailed(20),errmax(20),tol)
    call checkvalbuf_end('oyzz',nchecked(21),nfailed(21),errmax(21),tol)
    call checkvalbuf_end('ozzz',nchecked(22),nfailed(22),errmax(22),tol)
#endif
    call update_test_scores(ntests,nfailed,npass)

    !
    ! check that leaf_is_active matches what maketree set
    !
    nfailed_leaf(:) = 0
    nchecked_leaf(:) = 0
    ierrmax_leaf(:) = 0
    do i=1,int(ncells)
       ! only check leaf nodes (non-zero leaf_is_active)
       if (leaf_is_active_saved(i) /= 0) then
          call checkvalbuf(leaf_is_active(i),leaf_is_active_saved(i),0,'leaf_is_active', &
                          nfailed_leaf(1),nchecked_leaf(1),ierrmax_leaf(1))
       endif
    enddo
    if (nchecked_leaf(1) > 0) then
       call checkvalbuf_end('leaf_is_active',nchecked_leaf(1),nfailed_leaf(1),ierrmax_leaf(1),0)
    endif
    call update_test_scores(ntests,nfailed_leaf,npass)

    deallocate(old_tree)
    deallocate(leaf_is_active_saved)
 endif

 if (test_dualrounds .or. test_all) call test_dual_rounds(ntests,npass)

 if (id==master) write(*,"(/,a,/)") '<-- KDTREE TEST COMPLETE'

end subroutine test_kdtree

!-----------------------------------------------------------------------
!+
!   Test of the dual tree walk in rounds (dualwalk_rounds): the src
!   leaves and expansion of each leaf after the rounds must give the same
!   neighbours and expansion as a full walk from the root. Done from the
!   root, and from fake roots at depth R given by the walk of the tree
!   truncated at depth R (as the walk of the global tree under MPI)
!+
!-----------------------------------------------------------------------
subroutine test_dual_rounds(ntests,npass)
 use dim,         only:maxp,periodic
 use io,          only:id,master
 use neighkdtree, only:leaf_is_active,ncells,node,dualwalk_rounds
 use part,        only:npart,xyzh,hfact,massoftype,igas,maxphase,iphase,isetphase
 use kernel,      only:hfact_default
 use kdtree,      only:maketree,empty_tree,lenfgrav,irootnode,getneigh_dual_global,reset_cachestate_global
 use unifdis,     only:set_unifdis
 use testutils,   only:checkval,update_test_scores
 use mpidomain,   only:i_belong
 integer, intent(inout) :: ntests,npass
 integer, allocatable :: list(:),leaves(:),idepth_leaf(:),iflag(:)
 integer, allocatable :: roots(:),istart(:),icount(:),srclist(:),rstart(:),rcount(:),srcrem(:,:)
 real,    allocatable :: fnode_root(:,:)
 integer :: kslab,nrounds,nleaves,irootdepth,nneigh,iroots,k,nroots,nsrc
 integer :: nfailed(4)
 real    :: psep

 if (id==master) write(*,"(/,a)") '--> testing dual tree walk in rounds'
 psep  = 1./32.
 hfact = hfact_default
 npart = 0
 call set_unifdis('random',id,master,-0.5,0.5,-0.5,0.5,-0.5,0.5,&
                  psep,hfact,npart,xyzh,periodic,mask=i_belong)
 massoftype(igas) = 1000./npart
 if (maxphase==maxp) iphase(:) = isetphase(igas,iactive=.true.)

 allocate(list(maxp),leaves(maxp),idepth_leaf(maxp))

 do iroots=1,2
    do kslab=1,3
       call empty_tree(node)
       call maketree(node,xyzh,npart,leaf_is_active,ncells,apr_tree=.false.)
       call get_leaves(nleaves,leaves,idepth_leaf)

       if (iroots==1) then
          ! from the real root, which has itself left to open
          irootdepth = 0
          nroots = 1
          allocate(roots(1),istart(1),icount(1),srclist(1),fnode_root(lenfgrav,1))
          roots = irootnode
          istart = 0
          icount = 1
          srclist = irootnode
          fnode_root = 0.
       else
          ! fake roots at the minimum leaf depth, from the walk of the tree truncated
          ! at that depth (the global walk under MPI)
          irootdepth = minval(idepth_leaf(1:nleaves))
          allocate(iflag(ncells))
          iflag = 0
          call flag_depth(irootdepth,iflag)
          nroots = count(iflag /= 0)
          allocate(roots(nroots),istart(nroots),icount(nroots),srclist(nroots*nroots),fnode_root(lenfgrav,nroots))
          call reset_cachestate_global(int(ncells))
          nroots = 0
          nsrc   = 0
          do k=1,int(ncells)
             if (iflag(k) == 0) cycle
             nroots = nroots + 1
             call getneigh_dual_global(node,iflag,k,list,nneigh,fnode_root(:,nroots))
             roots(nroots)  = k
             istart(nroots) = nsrc
             icount(nroots) = nneigh
             srclist(nsrc+1:nsrc+nneigh) = list(1:nneigh)
             nsrc = nsrc + nneigh
          enddo
          deallocate(iflag)
       endif

       ! no remote src: the remote lists are empty
       allocate(rstart(nroots),rcount(nroots),srcrem(2,1))
       rstart = 0
       rcount = 0
       call dualwalk_rounds(kslab,nroots,roots,istart,icount,srclist,rstart,rcount,srcrem,fnode_root,nrounds)
       if (id==master) write(*,"(a,i2,a,i5,a,i2,a,i3,a)") ' root depth ',irootdepth,' (',nroots,&
                                                         ' roots), slab depth ',kslab,': ',nrounds,' rounds'
       nfailed = 0
       call checkval(nrounds > 0,.true.,nfailed(1),'rounds done')
       call check_leaf_walks(nleaves,leaves,nfailed(2:4))
       call update_test_scores(ntests,nfailed,npass)
       deallocate(roots,istart,icount,srclist,fnode_root,rstart,rcount,srcrem)
    enddo
 enddo

 deallocate(list,leaves,idepth_leaf)

end subroutine test_dual_rounds

!-----------------------------------------------------------------------
!+
!   leaves of the tree and their depth, found by walking down from the root
!+
!-----------------------------------------------------------------------
subroutine get_leaves(nleaves,leaves,idepth_leaf)
 use neighkdtree, only:leaf_is_active,node
 use kdtree,      only:irootnode
 integer, intent(out) :: nleaves,leaves(:),idepth_leaf(:)
 integer :: stack(2,128),istack,n,idep

 nleaves = 0
 istack  = 1
 stack(:,1) = (/irootnode,0/)
 do while(istack > 0)
    n    = stack(1,istack)
    idep = stack(2,istack)
    istack = istack - 1
    if (node(n)%leftchild == 0 .or. leaf_is_active(n) /= 0) then
       if (leaf_is_active(n) /= 0) then
          nleaves = nleaves + 1
          leaves(nleaves) = n
          idepth_leaf(nleaves) = idep
       endif
    else
       stack(:,istack+1) = (/node(n)%leftchild,idep+1/)
       stack(:,istack+2) = (/node(n)%rightchild,idep+1/)
       istack = istack + 2
    endif
 enddo

end subroutine get_leaves

!-----------------------------------------------------------------------
!+
!   flag the nodes at a given depth
!+
!-----------------------------------------------------------------------
subroutine flag_depth(idepth,iflag)
 use neighkdtree, only:leaf_is_active,node
 use kdtree,      only:irootnode
 integer, intent(in)    :: idepth
 integer, intent(inout) :: iflag(:)
 integer :: stack(2,128),istack,n,idep

 istack = 1
 stack(:,1) = (/irootnode,0/)
 do while(istack > 0)
    n    = stack(1,istack)
    idep = stack(2,istack)
    istack = istack - 1
    if (idep == idepth) then
       iflag(n) = 1
    elseif (node(n)%leftchild /= 0 .and. leaf_is_active(n) == 0) then
       stack(:,istack+1) = (/node(n)%leftchild,idep+1/)
       stack(:,istack+2) = (/node(n)%rightchild,idep+1/)
       istack = istack + 2
    endif
 enddo

end subroutine flag_depth

!-----------------------------------------------------------------------
!+
!   compare, for every leaf, the result of the walk in rounds with a full
!   walk from the root
!+
!-----------------------------------------------------------------------
subroutine check_leaf_walks(nleaves,leaves,nfailed)
 use dim,         only:maxp
 use neighkdtree, only:leaf_is_active,node,get_leaf_walk
 use kdtree,      only:getneigh_dual,use_cache,lenfgrav
 use testutils,   only:checkvalbuf,checkvalbuf_end
 integer, intent(in)    :: nleaves,leaves(:)
 integer, intent(inout) :: nfailed(3)
 integer, allocatable :: list(:),list_ref(:)
 integer :: i,j,icell,nneigh,nneigh_ref,nchecked(3),ierrmax(2)
 real    :: fnode(lenfgrav),fnode_ref(lenfgrav),errmax,xpos(3),xyzcache(1,1)
 real, parameter :: tol = 1.e-10
 logical :: use_cache_saved

 allocate(list(maxp),list_ref(maxp))
 use_cache_saved = use_cache
 use_cache = .false.
 xpos     = 0.
 nchecked = 0
 ierrmax  = 0
 errmax   = 0.
 do i=1,nleaves
    icell = leaves(i)
    call get_leaf_walk(icell,list,nneigh,xyzcache,0,fnode)
    call getneigh_dual(node,xpos,0.,0.,list_ref,nneigh_ref,xyzcache,0,leaf_is_active,&
                       .true.,.true.,fnode_ref,icell)
    call checkvalbuf(nneigh,nneigh_ref,0,'nneigh',nfailed(1),nchecked(1),ierrmax(1))
    if (nneigh == nneigh_ref) then
       call sort_int(list(1:nneigh))
       call sort_int(list_ref(1:nneigh))
       do j=1,nneigh
          call checkvalbuf(list(j),list_ref(j),0,'neighbours',nfailed(2),nchecked(2),ierrmax(2))
       enddo
    endif
    do j=1,lenfgrav
       call checkvalbuf(fnode(j),fnode_ref(j),tol,'fnode',nfailed(3),nchecked(3),errmax)
    enddo
 enddo
 use_cache = use_cache_saved
 call checkvalbuf_end('nneigh',nchecked(1),nfailed(1),ierrmax(1),0)
 call checkvalbuf_end('neighbours',nchecked(2),nfailed(2),ierrmax(2),0)
 call checkvalbuf_end('fnode',nchecked(3),nfailed(3),errmax,tol)
 deallocate(list,list_ref)

end subroutine check_leaf_walks

!-----------------------------------------------------------------------
!+
!   insertion sort of a short integer list
!+
!-----------------------------------------------------------------------
pure subroutine sort_int(a)
 integer, intent(inout) :: a(:)
 integer :: i,j,tmp

 do i=2,size(a)
    tmp = a(i)
    j = i - 1
    do while(j >= 1)
       if (a(j) <= tmp) exit
       a(j+1) = a(j)
       j = j - 1
    enddo
    a(j+1) = tmp
 enddo

end subroutine sort_int
end module testkdtree
