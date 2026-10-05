!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
submodule (kdtree) kdtree_walk
 use io,       only:fatal,id
 use part,     only:gravity
 use kernel,   only:radkern
 implicit none

contains

!----------------------------------------------------------------
!+
!  Routine to walk tree for neighbour search
!  (all particles within a given h_i and optionally within h_j)
!+
!----------------------------------------------------------------
module procedure getneigh
 integer :: maxcache
 integer :: n,istack,il,ir
 integer :: nstack(maxdepth)
 real :: dx,dy,dz,xsizej,rcutj
 real :: rcut,rcut2,r2
 real :: xoffset,yoffset,zoffset,tree_acc2
 logical :: open_tree_node
 logical :: global_walk
#ifdef GRAVITY
 real :: quads(6)
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
          if (istack+2 > maxdepth+1) call fatal('getneigh','stack overflow in getneigh')
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

end procedure getneigh

!----------------------------------------------------------------
!+
!  Routine to walk tree for neighbour search (SFMM version)
!  (all particles within a given h_i and optionally within h_j)
!  A dual tree walk is used to compute
!  every node-node interactions
!+
!----------------------------------------------------------------
module procedure getneigh_dual
 integer :: istack,i,iparent,idstbranch,idst,isrc,maxcache,ibase,nodestate
 integer :: branch(maxdepth),nparents,stack(3,maxstacksize),startwith(2)
 real    :: dx,dy,dz,xoffset,yoffset,zoffset
 real    :: tree_acc2
 real    :: fnode_acc(lenfgrav)
 logical :: stackit

 tree_acc2 = tree_accuracy*tree_accuracy

 if (ixyzcachesize > 0) then
    maxcache = size(xyzcache,1)
 else
    maxcache = 0
 endif

 call get_list_of_parent_nodes(icell,node,branch,nparents,startwith)

 neighnodecount_branch(1:nparents) = 0
 ! neighnode_branch(:,1:nparents) = 0 ! no need to reset neighnode_branch as neighnodecount_branch act as a switch
 fnode_branch(:,1:nparents) = 0.
 fnode_acc = 0.
 nneigh = 0
 istack = 0
 xoffset = 0.
 yoffset = 0.
 zoffset = 0.

 if (use_cache .and. startwith(2) > 0) then
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
       call node_interaction(idst,node(idst),node(isrc),tree_acc2,fnode_branch(:,idstbranch),stackit,xoffset,yoffset,zoffset)
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
    if (use_cache) then
       !$omp atomic read
       nodestate = cachestate(iparent)
       !$omp end atomic
       if (nodestate == 0) then ! first fence to avoid capture collision
          !$omp atomic capture
          nodestate = cachestate(iparent)
          cachestate(iparent) = max(cachestate(iparent),1)
          !$omp end atomic
          if (nodestate == 0) then ! if still the winner then cache
             !-- winner: publish fnode first ...
             fnodecache(1:lenfgrav,iparent) = fnode_branch(1:lenfgrav,i)
             !$omp atomic write
             cachestate(iparent) = 2
             !$omp end atomic

             !-- then store interaction list in the cache array if it fits
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
                   cachestate(iparent) = 3
                   !$omp end atomic
                endif
             endif
          endif
       elseif (nodestate>=2) then
          !-- fetch fnode from the cache array
          fnode_branch(1:lenfgrav,i) = fnodecache(1:lenfgrav,iparent)
       endif
    endif

    call get_sep(node(iparent)%xcen,node(branch(i-1))%xcen,dx,dy,dz,xoffset,yoffset,zoffset)
    fnode = fnode_acc + fnode_branch(:,i)
    call propagate_fnode_to_node(fnode_acc,fnode,dx,dy,dz)
 enddo

 fnode = fnode_acc + fnode_branch(:,1)

end procedure getneigh_dual

!-----------------------------------------------------------
!+
!  return list of parents of current node
!+
!-----------------------------------------------------------
subroutine get_list_of_parent_nodes(inode,node,parents,nparents,startwith)
 integer,      intent(in)  :: inode
 type(kdnode), intent(in)  :: node(:)
 integer,      intent(out) :: parents(:)
 integer,      intent(out) :: nparents
 integer,      intent(out) :: startwith(2)
 integer :: j,nodestate

 j = inode
 nparents  = 1
 parents   = 0
 startwith = 0
 parents(nparents) = j ! set first elem to inode to use parents for propagation
 do while (node(j)%parent  /=  0)
    j = node(j)%parent
    nparents = nparents + 1
    parents(nparents) = j
    !$omp atomic read
    nodestate = cachestate(j)
    !$omp end atomic
    if (nodestate==3 .and. startwith(2)==0) then
       ! deepest fully-cached ancestor: candidate pruned start
       startwith(1) = j
       startwith(2) = nparents
    elseif (nodestate<2 .and. startwith(2)/=0) then
       ! if non-cached ancestor node on the pruned branch reset the pruned start
       startwith = 0
    endif
 enddo

end subroutine get_list_of_parent_nodes
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
    xoffset = -dxbound*SIGN(1.0,dx)
    dx = dx + xoffset
 endif
 if (abs(dy) > hdly) then
    yoffset = -dybound*SIGN(1.0,dy)
    dy = dy + yoffset
 endif
 if (abs(dz) > hdlz) then
    zoffset = -dzbound*SIGN(1.0,dz)
    dz = dz + zoffset
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
!  Compute node node gravity interactions
!+
!-----------------------------------------------------------
subroutine open_nodes(stack,istack,srcnode,isrc,branch,idstbranch,&
                           listneigh,xyzcache,ixyzcachesize,nneigh,leaf_is_active,&
                           maxcache,xoffset,yoffset,zoffset)
 use io, only:fatal
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
       if (istack+1 > maxstacksize) call fatal('getneigh','stack overflow in getneigh')
       istack = istack + 1
       stack(1,istack) = idstnext
       stack(2,istack) = isrc
       stack(3,istack) = ibranchnext
    endif is_P2P
 else
    if (il /= 0) then
       if (istack+1 > maxstacksize) call fatal('getneigh','stack overflow in getneigh')
       istack = istack + 1
       stack(1,istack) = idstnext
       stack(2,istack) = il
       stack(3,istack) = ibranchnext
    endif
    if (ir /= 0) then
       if (istack+1 > maxstacksize) call fatal('getneigh','stack overflow in getneigh')
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
subroutine node_interaction(idst,node_dst,node_src,tree_acc2,fnode,stackit,xoffset,yoffset,zoffset)
 type(kdnode), intent(in)    :: node_dst,node_src
 integer,      intent(in)    :: idst
 real,         intent(in)    :: tree_acc2
 real,         intent(inout) :: fnode(lenfgrav)
 real,         intent(out)   :: xoffset,yoffset,zoffset
 logical,      intent(out)   :: stackit
 real    :: dx,dy,dz,r2
 real    :: rcut_dst,rcut_src,rcut,rcut2
 real    :: size_dst,size_src
 logical :: wellsep
 integer :: dststate
#ifdef GRAVITY
 real    :: dr1
#endif

 call get_sep(node_dst%xcen,node_src%xcen,dx,dy,dz,xoffset,yoffset,zoffset,r2)
 call get_node_size(node_dst,node_src,size_dst,size_src,rcut_dst,rcut_src)

 if (use_cache) then
    !$omp atomic read
    dststate = cachestate(idst)
    !$omp end atomic
 else
    dststate = 0
 endif

 rcut  = max(rcut_dst,rcut_src)
 rcut2 = (size_dst+size_src+rcut)**2
 wellsep = (tree_acc2*r2 > (size_dst+size_src)**2) .and. (r2 > rcut2)

 if (wellsep) then
#ifdef GRAVITY
    if (dststate<2) then
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

!-----------------------------------------------------------
!+
!  Compute the Taylor expansion coeffs between the node
!  centres using the quadrupole moments (p=3) (Dehnen 2002)
!+
!-----------------------------------------------------------
module procedure compute_M2L
 real :: qxx,qxy,qxz,qyy,qyz,qzz,dx2,dx3,dy2,dy3,dz2,dz3
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

 qxx = quads(1)
 qxy = quads(2)
 qxz = quads(3)
 qyy = quads(4)
 qyz = quads(5)
 qzz = quads(6)

 fnode(1)  = fnode(1)  + D1(1)*q0 + 0.5*(D3(1)*qxx + 2.*(D3(2)*qxy + D3(3)*qxz + D3(5)*qyz) + D3(4)*qyy + D3(6)*qzz)    ! C¹_x
 fnode(2)  = fnode(2)  + D1(2)*q0 + 0.5*(D3(2)*qxx + 2.*(D3(4)*qxy + D3(5)*qxz + D3(8)*qyz) + D3(7)*qyy + D3(9)*qzz)    ! C¹_y
 fnode(3)  = fnode(3)  + D1(3)*q0 + 0.5*(D3(3)*qxx + 2.*(D3(5)*qxy + D3(6)*qxz + D3(9)*qyz) + D3(8)*qyy + D3(10)*qzz)   ! C¹_z
 fnode(4)  = fnode(4)  - (D2(1) * q0)  ! C²_xx
 fnode(5)  = fnode(5)  - (D2(2) * q0)  ! C²_xy
 fnode(6)  = fnode(6)  - (D2(3) * q0)  ! C²_xz
 fnode(7)  = fnode(7)  - (D2(4) * q0)  ! C²_yy
 fnode(8)  = fnode(8)  - (D2(5) * q0)  ! C²_yz
 fnode(9)  = fnode(9)  - (D2(6) * q0)  ! C²_zz
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
 fnode(20) = fnode(20) + g0*q0 + 0.5*(D2(1)*qxx + D2(4)*qyy + D2(6)*qzz + 2*(D2(2)*qxy + D2(3)*qxz + D2(5)*qyz))! C⁰ (potential)

end procedure compute_M2L

#ifdef GRAVITY

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
module procedure expand_fgrav_in_taylor_series
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

end procedure expand_fgrav_in_taylor_series

end submodule kdtree_walk
