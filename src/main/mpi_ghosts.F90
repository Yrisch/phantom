!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module mpighosts
!
! Ghost particles of the neighbouring MPI domains. Each task sends the
! particles within reach of the domain (bounding box) of another task.
! They are written after npart in the particle arrays and built into the
! local tree as inactive particles, so the neighbour search finds them.
! Two sets of fields are sent: those read for a neighbour in density, and
! those read for a neighbour in force. The arrays are the ones passed to
! derivs (the predicted values in the step)
!
! :References: None
!
! :Owner: Yann Bernard
!
! :Runtime parameters: None
!
! :Dependencies: boundary, dim, io, kernel, mpi, mpiforce, mpiutils, part
!
 use io, only:nprocs,fatal
 implicit none
 private

 public :: copy_ghost,exchange_tree_ghosts,refresh_tree_ghosts,clear_tree_ghosts
 public :: check_pair_mirror,exchange_slabs

 ! fields sent: those read for a neighbour in density, or in force
 integer, parameter, public :: ighost_dens = 1, ighost_force = 2

 ! number of ghost particles after npart
 integer, public :: nghost_tree = 0
 ! h may grow in the density iterations: the ghosts are chosen with hfac_ghost*h
 real,    public :: hfac_ghost = 1.2
 ! the ghosts hold all the neighbours of the particles of this task with h <= hmax_ghost
 ! (density is a gather over radkern*h): beyond it they must be chosen again
 real,    public :: hmax_ghost = huge(1.)

 ! particles sent to each task, in task order: nsend_task(r) particles to task r-1
 integer, allocatable :: sendlist(:),nsend_task(:)

contains

!----------------------------------------------------------------
!+
!  select the ghost particles of each task and send them
!  (when the tree is built, after the domain decomposition)
!+
!----------------------------------------------------------------
subroutine exchange_tree_ghosts(iset,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,&
                                rad,radprop,dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level)
 integer,         intent(in)    :: iset,npart
 real,            intent(inout) :: xyzh(:,:),vxyzu(:,:),fxyzu(:,:),fext(:,:),Bevol(:,:),rad(:,:)
 real,            intent(inout) :: radprop(:,:),dustprop(:,:),dustfrac(:,:),filfac(:),eos_vars(:,:)
 real,            intent(inout) :: dens(:),metrics(:,:,:,:)
 real(kind=4),    intent(inout) :: divcurlv(:,:),divcurlB(:,:)
 integer(kind=1), intent(inout) :: apr_level(:)

 call select_ghosts(npart,xyzh)
 call send_ghosts(iset,.false.,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,&
                  rad,radprop,dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level)

end subroutine exchange_tree_ghosts

!----------------------------------------------------------------
!+
!  send the values of the same ghost particles again
!+
!----------------------------------------------------------------
subroutine refresh_tree_ghosts(iset,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,&
                               rad,radprop,dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level)
 integer,         intent(in)    :: iset,npart
 real,            intent(inout) :: xyzh(:,:),vxyzu(:,:),fxyzu(:,:),fext(:,:),Bevol(:,:),rad(:,:)
 real,            intent(inout) :: radprop(:,:),dustprop(:,:),dustfrac(:,:),filfac(:),eos_vars(:,:)
 real,            intent(inout) :: dens(:),metrics(:,:,:,:)
 real(kind=4),    intent(inout) :: divcurlv(:,:),divcurlB(:,:)
 integer(kind=1), intent(inout) :: apr_level(:)

 call send_ghosts(iset,.true.,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,&
                  rad,radprop,dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level)

end subroutine refresh_tree_ghosts

!----------------------------------------------------------------
!+
!  remove the ghost particles (before the particles are moved
!  between the tasks). iphase is left as it is: particles set
!  later in these slots may rely on it
!+
!----------------------------------------------------------------
subroutine clear_tree_ghosts(npart)
 use part, only:xyzh
 integer, intent(in) :: npart
 integer :: nclear

 nclear = min(nghost_tree,size(xyzh,2)-npart)
 if (nclear > 0) xyzh(:,npart+1:npart+nclear) = 0.
 nghost_tree = 0
 hmax_ghost  = huge(1.)

end subroutine clear_tree_ghosts

!----------------------------------------------------------------
!+
!  particles to send to each task: those within radkern*max(h,hmax)
!  of the bounding box of its particles. The tasks whose box is out
!  of reach of the box of this task are skipped
!+
!----------------------------------------------------------------
subroutine select_ghosts(npart,xyzh)
 use io,   only:id
 use part, only:isdead_or_accreted
!$ use omp_lib, only:omp_get_max_threads
 integer, intent(in) :: npart
 real,    intent(in) :: xyzh(:,:)
 real    :: box(7),boxes(7,nprocs),bmin(3),bmax(3),hmax
 integer :: i,irank,ic,nchunk,ntot
 integer, allocatable :: ncount(:)

 ! bounding box and hmax of the particles of this task, for all tasks
 bmin =  huge(1.)
 bmax = -huge(1.)
 hmax = 0.
 !$omp parallel do default(none) shared(npart,xyzh) private(i) &
 !$omp reduction(min:bmin) reduction(max:bmax,hmax)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       bmin = min(bmin,xyzh(1:3,i))
       bmax = max(bmax,xyzh(1:3,i))
       hmax = max(hmax,xyzh(4,i))
    endif
 enddo
 !$omp end parallel do
 box = (/bmin,bmax,hmax/)
 call allgather_box(box,boxes)
 ! a particle j of another task is sent if its distance to the box of this task is
 ! below radkern*hfac_ghost*max(hj,hmax): all the neighbours within radkern*h are
 ! sent for any particle of this task with h <= hfac_ghost*hmax
 hmax_ghost = hfac_ghost*hmax

 if (.not.allocated(nsend_task)) allocate(nsend_task(nprocs))
 if (.not.allocated(sendlist)) allocate(sendlist(1024))
 nsend_task = 0
 ntot   = 0
 nchunk = 1
!$ nchunk = omp_get_max_threads()
 allocate(ncount(0:nchunk))
 do irank=0,nprocs-1
    if (irank /= id .and. within_reach(bmin,bmax,hmax,boxes(:,irank+1))) then
       ! count per chunk of particles, then fill in particle order
       ncount = 0
       !$omp parallel do default(none) schedule(static) &
       !$omp shared(nchunk,npart,xyzh,boxes,irank,ncount) private(ic,i)
       do ic=1,nchunk
          do i=chunk_start(ic),chunk_start(ic+1)-1
             if (is_ghost(i)) ncount(ic) = ncount(ic) + 1
          enddo
       enddo
       !$omp end parallel do
       ncount(0) = ntot
       do ic=1,nchunk
          ncount(ic) = ncount(ic) + ncount(ic-1)
       enddo
       if (ncount(nchunk) > size(sendlist)) call grow(sendlist,2*ncount(nchunk),ntot)
       !$omp parallel do default(none) schedule(static) &
       !$omp shared(nchunk,npart,xyzh,boxes,irank,ncount,sendlist) private(ic,i)
       do ic=1,nchunk
          do i=chunk_start(ic),chunk_start(ic+1)-1
             if (is_ghost(i)) then
                ncount(ic-1) = ncount(ic-1) + 1
                sendlist(ncount(ic-1)) = i
             endif
          enddo
       enddo
       !$omp end parallel do
       nsend_task(irank+1) = ncount(nchunk-1) - ntot
       ntot = ncount(nchunk-1)
    endif
 enddo

contains

integer function chunk_start(ic)
 integer, intent(in) :: ic

 chunk_start = int((int(ic-1,8)*npart)/nchunk) + 1

end function chunk_start

logical function is_ghost(i)
 integer, intent(in) :: i

 is_ghost = .false.
 if (.not.isdead_or_accreted(xyzh(4,i))) &
    is_ghost = within_reach(xyzh(1:3,i),xyzh(1:3,i),xyzh(4,i),boxes(:,irank+1))

end function is_ghost

end subroutine select_ghosts

!----------------------------------------------------------------
!+
!  true if a box (or a particle, amin = amax) is within reach of
!  the box b (min, max, hmax): distance < radkern*hfac_ghost*max(h,hmax)
!+
!----------------------------------------------------------------
pure logical function within_reach(amin,amax,h,b)
 use kernel,   only:radkern
 use dim,      only:periodic
 use boundary, only:dxbound,dybound,dzbound
 real, intent(in) :: amin(3),amax(3),h,b(7)
 real :: d(3),l(3)
 integer :: k

 l = (/dxbound,dybound,dzbound/)
 do k=1,3
    d(k) = max(b(k)-amax(k),amin(k)-b(k+3),0.)
    if (periodic) d(k) = min(d(k),max(b(k)+l(k)-amax(k),amin(k)-b(k+3)-l(k),0.),&
                                  max(b(k)-l(k)-amax(k),amin(k)-b(k+3)+l(k),0.))
 enddo
 within_reach = dot_product(d,d) < (radkern*hfac_ghost*max(h,b(7)))**2

end function within_reach

!----------------------------------------------------------------
!+
!  send the particles of sendlist and write the particles received
!  after npart (refresh: the same particles as the last selection)
!+
!----------------------------------------------------------------
subroutine send_ghosts(iset,refresh,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,&
                       rad,radprop,dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level)
 integer,         intent(in)    :: iset,npart
 logical,         intent(in)    :: refresh
 real,            intent(inout) :: xyzh(:,:),vxyzu(:,:),fxyzu(:,:),fext(:,:),Bevol(:,:),rad(:,:)
 real,            intent(inout) :: radprop(:,:),dustprop(:,:),dustfrac(:,:),filfac(:),eos_vars(:,:)
 real,            intent(inout) :: dens(:),metrics(:,:,:,:)
 real(kind=4),    intent(inout) :: divcurlv(:,:),divcurlB(:,:)
 integer(kind=1), intent(inout) :: apr_level(:)
 real, allocatable :: sendbuf(:),recvbuf(:)
 integer :: nfield,ipos,k,n,nsend(nprocs),nrecv(nprocs)
 logical :: anyflag
 real    :: buf1(4096)

 ! number of fields of a particle
 ipos = 0
 call copy_ghost(iset,.true.,1,buf1,ipos,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,&
                 rad,radprop,dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level)
 nfield = ipos

 nsend = nsend_task*nfield
 allocate(sendbuf(max(sum(nsend),1)))
 !$omp parallel do default(none) schedule(static) &
 !$omp shared(iset,nsend_task,sendlist,sendbuf,nfield) &
 !$omp shared(xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,rad,radprop) &
 !$omp shared(dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level) &
 !$omp private(k,ipos)
 do k=1,sum(nsend_task)
    ipos = (k-1)*nfield
    call copy_ghost(iset,.true.,sendlist(k),sendbuf,ipos,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,&
                    rad,radprop,dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level)
 enddo
 !$omp end parallel do

 call exchange_slabs(sendbuf,nsend,recvbuf,nrecv,.true.,anyflag)

 n = sum(nrecv)/nfield
 if (refresh) then
    if (n /= nghost_tree) call fatal('mpighosts','ghost particles refreshed differ from those of the tree')
 elseif (npart + n > size(xyzh,2)) then
    call fatal('mpighosts','no room for the ghost particles: increase maxp')
 endif
 !$omp parallel do default(none) schedule(static) &
 !$omp shared(iset,n,npart,recvbuf,nfield) &
 !$omp shared(xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,rad,radprop) &
 !$omp shared(dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level) &
 !$omp private(k,ipos)
 do k=1,n
    ipos = (k-1)*nfield
    call copy_ghost(iset,.false.,npart+k,recvbuf,ipos,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,&
                    rad,radprop,dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level)
 enddo
 !$omp end parallel do
 nghost_tree = n

end subroutine send_ghosts

!----------------------------------------------------------------
!+
!  bounding boxes of all the tasks
!+
!----------------------------------------------------------------
subroutine allgather_box(box,boxes)
#ifdef MPI
 use mpi
 use mpiutils, only:mpierr
#endif
 real, intent(in)  :: box(7)
 real, intent(out) :: boxes(7,nprocs)

#ifdef MPI
 call MPI_ALLGATHER(box,7,MPI_REAL8,boxes,7,MPI_REAL8,MPI_COMM_WORLD,mpierr)
#else
 boxes(:,1) = box
#endif

end subroutine allgather_box

!----------------------------------------------------------------
!+
!  pack (or unpack) the fields of particle i read for a neighbour
!  in density (iset = ighost_dens) or in force (iset = ighost_force)
!  into buf, from position ipos. Only the arrays allocated for all
!  the particles are copied (the others are not used in this setup)
!+
!----------------------------------------------------------------
subroutine copy_ghost(iset,pack,i,buf,ipos,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,&
                      rad,radprop,dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level)
 use part, only:iphase,rho,gradh,alphaind,dvdx,eta_nimhd,fxyz_dragold,ibin_old
 integer,         intent(in)    :: iset
 logical,         intent(in)    :: pack
 integer,         intent(in)    :: i
 real,            intent(inout) :: buf(:)
 integer,         intent(inout) :: ipos
 real,            intent(inout) :: xyzh(:,:),vxyzu(:,:),fxyzu(:,:),fext(:,:),Bevol(:,:),rad(:,:)
 real,            intent(inout) :: radprop(:,:),dustprop(:,:),dustfrac(:,:),filfac(:),eos_vars(:,:)
 real,            intent(inout) :: dens(:),metrics(:,:,:,:)
 real(kind=4),    intent(inout) :: divcurlv(:,:),divcurlB(:,:)
 integer(kind=1), intent(inout) :: apr_level(:)
 integer :: nfull

 nfull = size(xyzh,2)
 ! read for a neighbour in both
 call copy_r8(xyzh)
 call copy_i1(iphase)
 call copy_i1(apr_level)
 call copy_r8(vxyzu)
 call copy_r8_1(rho)
 call copy_r8(Bevol)
 call copy_r8(rad)
 select case(iset)
 case(ighost_dens)
    ! acceleration for the Cullen & Dehnen switch
    call copy_r8(fxyzu)
    call copy_r8(fext)
 case(ighost_force)
    call copy_r4(gradh)
    call copy_r4(divcurlv)
    call copy_r4(divcurlB)
    call copy_r4(dvdx)
    call copy_r4(alphaind)
    call copy_r8(eos_vars)
    call copy_r8(radprop)
    call copy_r8(dustprop)
    call copy_r8(dustfrac)
    call copy_r8_1(filfac)
    call copy_r8(eta_nimhd)
    call copy_r8(fxyz_dragold)
    call copy_r8_1(dens)
    if (size(metrics,4) == nfull) then
       if (pack) then
          buf(ipos+1:ipos+size(metrics(:,:,:,i))) = reshape(metrics(:,:,:,i),(/size(metrics(:,:,:,i))/))
       else
          metrics(:,:,:,i) = reshape(buf(ipos+1:ipos+size(metrics(:,:,:,i))),shape(metrics(:,:,:,i)))
       endif
       ipos = ipos + size(metrics(:,:,:,i))
    endif
    call copy_i1(ibin_old)
 end select

contains

subroutine copy_r8(a)
 real, intent(inout) :: a(:,:)
 integer :: nf

 if (size(a,2) /= nfull) return
 nf = size(a,1)
 if (pack) then
    buf(ipos+1:ipos+nf) = a(:,i)
 else
    a(:,i) = buf(ipos+1:ipos+nf)
 endif
 ipos = ipos + nf

end subroutine copy_r8

subroutine copy_r4(a)
 real(kind=4), intent(inout) :: a(:,:)
 integer :: nf

 if (size(a,2) /= nfull) return
 nf = size(a,1)
 if (pack) then
    buf(ipos+1:ipos+nf) = real(a(:,i))
 else
    a(:,i) = real(buf(ipos+1:ipos+nf),kind=4)
 endif
 ipos = ipos + nf

end subroutine copy_r4

subroutine copy_r8_1(a)
 real, intent(inout) :: a(:)

 if (size(a) /= nfull) return
 if (pack) then
    buf(ipos+1) = a(i)
 else
    a(i) = buf(ipos+1)
 endif
 ipos = ipos + 1

end subroutine copy_r8_1

subroutine copy_i1(a)
 integer(kind=1), intent(inout) :: a(:)

 if (size(a) /= nfull) return
 if (pack) then
    buf(ipos+1) = real(a(i))
 else
    a(i) = int(nint(buf(ipos+1)),kind=1)
 endif
 ipos = ipos + 1

end subroutine copy_i1

end subroutine copy_ghost

!----------------------------------------------------------------
!+
!  grow an array, keeping its first nkeep elements
!+
!----------------------------------------------------------------
subroutine grow(a,n,nkeep)
 integer, allocatable, intent(inout) :: a(:)
 integer,              intent(in)    :: n,nkeep
 integer, allocatable :: tmp(:)

 allocate(tmp(n))
 tmp(1:nkeep) = a(1:nkeep)
 call move_alloc(tmp,a)

end subroutine grow

!----------------------------------------------------------------
!+
!  check that the remote pairs of the global dual tree walk are
!  mirrored: (a,b) on this task for remote r <=> (b,a) on task r
!  pairs(:,i) = (dst,src) global node indices, owner of src
!  (pairs with a local src are skipped)
!  returns the number of mismatched pairs found on this task
!+
!----------------------------------------------------------------
subroutine check_pair_mirror(npairs,pairs,nmismatch)
#ifdef MPI
 use mpi
 use mpiutils, only:mpierr
 use sortutils, only:indexx
 use io,        only:id
#endif
 integer, intent(in)  :: npairs
 integer, intent(in)  :: pairs(:,:)
 integer, intent(out) :: nmismatch
#ifdef MPI
 integer :: nsend(nprocs),nrecv(nprocs),isdispl(nprocs),irdispl(nprocs),ifill(nprocs)
 integer :: i,j,irank,n,nrecvtot
 integer, allocatable :: sendbuf(:),recvbuf(:),indx(:),jndx(:)
 integer(kind=8), allocatable :: key(:),keyrecv(:)

 nmismatch = 0

 ! group my pairs by owner of src
 nsend = 0
 do i=1,npairs
    if (pairs(3,i) == id) cycle
    nsend(pairs(3,i)+1) = nsend(pairs(3,i)+1) + 2
 enddo
 isdispl(1) = 0
 do irank=2,nprocs
    isdispl(irank) = isdispl(irank-1) + nsend(irank-1)
 enddo
 allocate(sendbuf(max(2*npairs,1)))
 ifill = isdispl
 do i=1,npairs
    if (pairs(3,i) == id) cycle
    irank = pairs(3,i)+1
    sendbuf(ifill(irank)+1) = pairs(1,i)
    sendbuf(ifill(irank)+2) = pairs(2,i)
    ifill(irank) = ifill(irank) + 2
 enddo

 call MPI_ALLTOALL(nsend,1,MPI_INTEGER,nrecv,1,MPI_INTEGER,MPI_COMM_WORLD,mpierr)
 irdispl(1) = 0
 do irank=2,nprocs
    irdispl(irank) = irdispl(irank-1) + nrecv(irank-1)
 enddo
 nrecvtot = sum(nrecv)
 allocate(recvbuf(max(nrecvtot,1)))
 call MPI_ALLTOALLV(sendbuf,nsend,isdispl,MPI_INTEGER,recvbuf,nrecv,irdispl,MPI_INTEGER,&
                    MPI_COMM_WORLD,mpierr)

 ! what task r sent as (a,b) must be (b,a) in my list for r: compare sorted (dst,src) keys
 n = max(maxval(nsend),1)/2 + 1
 allocate(key(n),keyrecv(n),indx(n),jndx(n))
 do irank=1,nprocs
    if (nrecv(irank) /= nsend(irank)) then
       nmismatch = nmismatch + abs(nrecv(irank) - nsend(irank))/2
       cycle
    endif
    n = nsend(irank)/2
    if (n == 0) cycle
    do j=1,n
       key(j)     = int(sendbuf(isdispl(irank)+2*j-1),8)*2_8**31 + int(sendbuf(isdispl(irank)+2*j),8)
       keyrecv(j) = int(recvbuf(irdispl(irank)+2*j),8)*2_8**31   + int(recvbuf(irdispl(irank)+2*j-1),8)
    enddo
    call indexx(n,key,indx)
    call indexx(n,keyrecv,jndx)
    do j=1,n
       if (keyrecv(jndx(j)) /= key(indx(j))) nmismatch = nmismatch + 1
    enddo
 enddo
#else
 nmismatch = 0
#endif

end subroutine check_pair_mirror

!----------------------------------------------------------------
!+
!  exchange of the nodes received between two rounds of the walk:
!  nsend(r) reals of sendbuf (grouped by destination task) go to
!  task r-1, and nrecv(r) reals are received from it in recvbuf.
!  The sizes carry a flag: anyactive if one task is active
!+
!----------------------------------------------------------------
subroutine exchange_slabs(sendbuf,nsend,recvbuf,nrecv,active,anyactive)
#ifdef MPI
 use mpi
 use mpiutils, only:mpierr
#endif
 real,              intent(in)    :: sendbuf(:)
 integer,           intent(in)    :: nsend(nprocs)
 real, allocatable, intent(inout) :: recvbuf(:)
 integer,           intent(out)   :: nrecv(nprocs)
 logical,           intent(in)    :: active
 logical,           intent(out)   :: anyactive
#ifdef MPI
 integer :: isdispl(nprocs),irdispl(nprocs),irank,isize(2,nprocs),irsize(2,nprocs)

 isize(1,:) = nsend
 isize(2,:) = merge(1,0,active)
 call MPI_ALLTOALL(isize,2,MPI_INTEGER,irsize,2,MPI_INTEGER,MPI_COMM_WORLD,mpierr)
 nrecv     = irsize(1,:)
 anyactive = any(irsize(2,:) > 0)
 isdispl(1) = 0
 irdispl(1) = 0
 do irank=2,nprocs
    isdispl(irank) = isdispl(irank-1) + nsend(irank-1)
    irdispl(irank) = irdispl(irank-1) + nrecv(irank-1)
 enddo
 if (allocated(recvbuf)) deallocate(recvbuf)
 allocate(recvbuf(max(sum(nrecv),1)))
 call MPI_ALLTOALLV(sendbuf,nsend,isdispl,MPI_REAL8,recvbuf,nrecv,irdispl,MPI_REAL8,&
                    MPI_COMM_WORLD,mpierr)
#else
 nrecv = nsend
 anyactive = active
 if (allocated(recvbuf)) deallocate(recvbuf)
 allocate(recvbuf(max(sum(nrecv),1)))
 recvbuf(1:sum(nrecv)) = sendbuf(1:sum(nsend))
#endif

end subroutine exchange_slabs

end module mpighosts
