!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module mpiforce
!
! None
!
! :References: None
!
! :Owner: Conrad Chan
!
! :Runtime parameters: None
!
! :Dependencies: dim, io, mpi
!
 use io,       only:nprocs,fatal
 use dim,      only:minpart,maxfsum,maxxpartveciforce

 implicit none
 private

 public :: cellforce
 public :: stackforce
 public :: get_mpitype_of_cellforce
 public :: free_mpitype_of_cellforce
 public :: check_pair_mirror

 integer, parameter :: ndata = 20 ! number of elements in the cell (including padding)
 integer, parameter :: nbytes_cellforce = 8 * maxxpartveciforce * minpart + &  !  xpartvec(maxxpartveciforce,minpart)
                                          8 * maxfsum * minpart           + &  !  fsums(maxfsum,minpart)
                                          8 * 20                          + &  !  fgrav(20)
                                          8 * 3                           + &  !  xpos(3)
                                          8                               + &  !  xsizei
                                          8                               + &  !  rcuti
                                          8 * minpart                     + &  !  tsmin(minpart)
                                          8 * minpart                     + &  !  vsigmax(minpart)
                                          4                               + &  !  icell
                                          4                               + &  !  npcell
                                          4 * minpart                     + &  !  arr_index(minpart)
                                          4                               + &  !  ndrag
                                          4                               + &  !  nstokes
                                          4                               + &  !  nsuper
                                          4                               + &  !  owner
                                          4                               + &  !  waiting_index
                                          1 * minpart                     + &  !  iphase(minpart)
                                          1 * minpart                     + &  !  ibinneigh(minpart)
                                          1 * minpart                          !  apr_level

 type cellforce
    sequence
    real             :: xpartvec(maxxpartveciforce,minpart)
    real             :: fsums(maxfsum,minpart)
    real             :: fgrav(20)
    real             :: xpos(3)
    real             :: xsizei
    real             :: rcuti
    real             :: tsmin(minpart)
    real             :: vsigmax(minpart)
    integer          :: icell
    integer          :: npcell                                 ! number of particles in here
    integer          :: arr_index(minpart)
    integer          :: ndrag
    integer          :: nstokes
    integer          :: nsuper
    integer          :: owner                                  ! id of the process that owns this
    integer          :: waiting_index
    integer(kind=1)  :: iphase(minpart)
    integer(kind=1)  :: ibinneigh(minpart)
    integer(kind=1)  :: apr(minpart)                           ! apr resolution level (not in xpartvec because integer)

    ! pad the array to 8-byte boundaries
    integer(kind=1)  :: pad(8 - mod(nbytes_cellforce, 8)) !padding to maintain alignment of elements
 end type cellforce

 type stackforce
    sequence
    type(cellforce), pointer  :: cells(:)
    integer                   :: maxlength = 0
    integer                   :: n = 0
    integer                   :: number
    integer                   :: ibuffer   ! to avoid ifort error
 end type stackforce

contains

subroutine get_mpitype_of_cellforce(dtype)
#ifdef MPI
 use mpi
#endif
 integer, intent(out) :: dtype
#ifdef MPI
 integer                         :: nblock, blens(ndata), mpitypes(ndata)
 integer(kind=MPI_ADDRESS_KIND)  :: disp(ndata)

 type(cellforce)                 :: cell
 integer(kind=MPI_ADDRESS_KIND)  :: addr,start,lb,extent
 integer                         :: mpierr

 nblock = 0

 call MPI_GET_ADDRESS(cell,start,mpierr)

 nblock = nblock + 1
 blens(nblock) = size(cell%xpartvec)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%xpartvec,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%fsums)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%fsums,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%fgrav)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%fgrav,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%xpos)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%xpos,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%xsizei,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%rcuti,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%tsmin)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%tsmin,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%vsigmax)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%vsigmax,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%icell,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%npcell,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%arr_index)
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%arr_index,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%ndrag,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%nstokes,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%nsuper,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%owner,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%waiting_index,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%iphase)
 mpitypes(nblock) = MPI_INTEGER1
 call MPI_GET_ADDRESS(cell%iphase,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%ibinneigh)
 mpitypes(nblock) = MPI_INTEGER1
 call MPI_GET_ADDRESS(cell%ibinneigh,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%apr)
 mpitypes(nblock) = MPI_INTEGER1
 call MPI_GET_ADDRESS(cell%apr,addr,mpierr)
 disp(nblock) = addr - start

 ! padding must come last
 nblock = nblock + 1
 blens(nblock) = 8 - mod(nbytes_cellforce, 8)
 mpitypes(nblock) = MPI_INTEGER1
 call MPI_GET_ADDRESS(cell%pad,addr,mpierr)
 disp(nblock) = addr - start

 call MPI_TYPE_CREATE_STRUCT(nblock,blens(1:nblock),disp(1:nblock),mpitypes(1:nblock),dtype,mpierr)
 call MPI_TYPE_COMMIT(dtype,mpierr)

 ! check extent okay
 call MPI_TYPE_GET_EXTENT(dtype,lb,extent,mpierr)
 if (extent /= sizeof(cell)) then
    call fatal('mpi_force','MPI_TYPE_GET_EXTENT has calculated the extent incorrectly')
 endif

#else
 dtype = 0
#endif

end subroutine get_mpitype_of_cellforce

subroutine free_mpitype_of_cellforce(dtype)
#ifdef MPI
 use mpi
#endif
 integer, intent(inout) :: dtype
#ifdef MPI
 integer                :: mpierr

 call MPI_Type_free(dtype,mpierr)
#else
 dtype = 0
#endif
end subroutine free_mpitype_of_cellforce

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

end module mpiforce
