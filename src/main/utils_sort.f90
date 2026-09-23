!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module sortutils
!
! contains low level sorting utilities
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
 public :: indexx,indexxfunc,Knnfunc,parqsort,find_rank
 public :: sort_by_radius
 public :: radixsort_i8
 public :: r2func,r2func_origin,set_r2func_origin
 interface indexx
  module procedure indexx_r4, indexx_i8
 end interface indexx

 private
 real, private :: x0,y0,z0

contains

!----------------------------------------------------------------
!+
!  function returning radius squared given xyzh
!+
!----------------------------------------------------------------
real function r2func(xyzh)
 real, intent(in) :: xyzh(4)

 r2func = xyzh(1)**2 + xyzh(2)**2 + xyzh(3)**2

end function r2func

!----------------------------------------------------------------
!+
!  function returning radius squared given xyzh and origin pos
!+
!----------------------------------------------------------------
subroutine set_r2func_origin(xorigin,yorigin,zorigin)
 real, intent(in) :: xorigin,yorigin,zorigin

 x0 = xorigin
 y0 = yorigin
 z0 = zorigin

end subroutine set_r2func_origin

real function r2func_origin(xyzh)
 real, intent(in) :: xyzh(4)

 r2func_origin = (xyzh(1)-x0)**2 + (xyzh(2)-y0)**2 + (xyzh(3)-z0)**2

end function r2func_origin

real function r2func_pos(xyzh,xpos)
 real, intent(in) :: xyzh(4),xpos(3)

 r2func_pos = (xyzh(1)-xpos(1))**2 + (xyzh(2)-xpos(2))**2 + (xyzh(3)-xpos(3))**2

end function r2func_pos
!----------------------------------------------------------------
!+
!  standard sorting routine using Quicksort for real*4
!+
!----------------------------------------------------------------
subroutine indexx_r4(n, arr, indx)
 integer, intent(in)  :: n
 real,    intent(in)  :: arr(n)
 integer, intent(out) :: indx(n)
 integer, parameter :: m=7, nstack=500

 integer :: i,j,k,l,ir,jstack,indxt,itemp
 integer :: istack(nstack)
 real :: a

 do j = 1, n
    indx(j) = j
 enddo
 jstack = 0
 l = 1
 ir = n

1 if (ir - l < m) then
    do j = l + 1, ir
       indxt = indx(j)
       a = arr(indxt)
       do i = j - 1, 1, -1
          if (arr(indx(i)) <= a) goto 2
          indx(i + 1) = indx(i)
       enddo
       i = 0
2      indx(i + 1) = indxt
    enddo
    if (jstack==0) return
    ir = istack(jstack)
    l = istack(jstack - 1)
    jstack = jstack - 2
 else
    k = (l + ir)/2
    itemp = indx(k)
    indx(k) = indx(l + 1)
    indx(l + 1) = itemp
    if (arr(indx(l + 1)) > arr(indx(ir))) then
       itemp = indx(l + 1)
       indx(l + 1) = indx(ir)
       indx(ir) = itemp
    endif
    if (arr(indx(l)) > arr(indx(ir))) then
       itemp = indx(l)
       indx(l) = indx(ir)
       indx(ir) = itemp
    endif
    if (arr(indx(l + 1)) > arr(indx(l))) then
       itemp = indx(l + 1)
       indx(l + 1) = indx(l)
       indx(l) = itemp
    endif
    i = l + 1
    j = ir
    indxt = indx(l)
    a = arr(indxt)

3   continue
    i = i + 1
    if (arr(indx(i)) < a) goto 3
4   continue
    j = j - 1
    if (arr(indx(j)) > a) goto 4
    if (j < i) goto 5
    itemp = indx(i)
    indx(i) = indx(j)
    indx(j) = itemp
    goto 3

5   indx(l) = indx(j)
    indx(j) = indxt
    jstack = jstack + 2
    if (jstack > nstack) then
       print*,'fatal error!!! stacksize exceeded in sort'
       print*,'need to set parameter nstack higher in subroutine indexx '
       stop
    endif
    if (ir - i + 1 >= j - l) then
       istack(jstack) = ir
       istack(jstack - 1) = i
       ir = j - 1
    else
       istack(jstack) = j - 1
       istack(jstack - 1) = l
       l = i
    endif
 endif

 goto 1

end subroutine indexx_r4

!----------------------------------------------------------------
!+
!  standard sorting routine using Quicksort for integer*8
!+
!----------------------------------------------------------------
subroutine indexx_i8(n, arr, indx)
 integer,         intent(in)  :: n
 integer(kind=8), intent(in)  :: arr(n)
 integer,         intent(out) :: indx(n)
 integer,            parameter :: m=7, nstack=500

 integer :: i,j,k,l,ir,jstack,indxt,itemp
 integer :: istack(nstack)
 integer(kind=8) :: a

 do j = 1, n
    indx(j) = j
 enddo
 jstack = 0
 l = 1
 ir = n

1 if (ir - l < m) then
    do j = l + 1, ir
       indxt = indx(j)
       a = arr(indxt)
       do i = j - 1, 1, -1
          if (arr(indx(i)) <= a) goto 2
          indx(i + 1) = indx(i)
       enddo
       i = 0
2      indx(i + 1) = indxt
    enddo
    if (jstack==0) return
    ir = istack(jstack)
    l = istack(jstack - 1)
    jstack = jstack - 2
 else
    k = (l + ir)/2
    itemp = indx(k)
    indx(k) = indx(l + 1)
    indx(l + 1) = itemp
    if (arr(indx(l + 1)) > arr(indx(ir))) then
       itemp = indx(l + 1)
       indx(l + 1) = indx(ir)
       indx(ir) = itemp
    endif
    if (arr(indx(l)) > arr(indx(ir))) then
       itemp = indx(l)
       indx(l) = indx(ir)
       indx(ir) = itemp
    endif
    if (arr(indx(l + 1)) > arr(indx(l))) then
       itemp = indx(l + 1)
       indx(l + 1) = indx(l)
       indx(l) = itemp
    endif
    i = l + 1
    j = ir
    indxt = indx(l)
    a = arr(indxt)

3   continue
    i = i + 1
    if (arr(indx(i)) < a) goto 3
4   continue
    j = j - 1
    if (arr(indx(j)) > a) goto 4
    if (j < i) goto 5
    itemp = indx(i)
    indx(i) = indx(j)
    indx(j) = itemp
    goto 3

5   indx(l) = indx(j)
    indx(j) = indxt
    jstack = jstack + 2
    if (jstack > nstack) then
       print*,'fatal error!!! stacksize exceeded in sort'
       print*,'need to set parameter nstack higher in subroutine indexx '
       stop
    endif
    if (ir - i + 1 >= j - l) then
       istack(jstack) = ir
       istack(jstack - 1) = i
       ir = j - 1
    else
       istack(jstack) = j - 1
       istack(jstack - 1) = l
       l = i
    endif
 endif

 goto 1
end subroutine indexx_i8

!----------------------------------------------------------------
!+
!  customised low-memory sorting routine using Quicksort
!  sort key value on-the-fly by calling the function func
!  which can be any function of the particle positions
!+
!----------------------------------------------------------------
subroutine indexxfunc(n, func, xyzh, indx)
 integer, intent(in)  :: n
 real,    intent(in)  :: xyzh(:,:)
 integer, intent(out) :: indx(n)
 integer, parameter :: m=7, nstack=500
 real, external :: func

 integer :: i,j,k,l,ir,jstack,indxt,itemp
 integer :: istack(nstack)
 real :: a

 do j = 1, n
    indx(j) = j
 enddo
 jstack = 0
 l = 1
 ir = n

1 if (ir - l < m) then
    do j = l + 1, ir
       indxt = indx(j)
       a = func(xyzh(:,indxt))
       do i = j - 1, 1, -1
          if (func(xyzh(:,indx(i))) <= a) goto 2
          indx(i + 1) = indx(i)
       enddo
       i = 0
2      indx(i + 1) = indxt
    enddo
    if (jstack==0) return
    ir = istack(jstack)
    l = istack(jstack - 1)
    jstack = jstack - 2
 else
    k = (l + ir)/2
    itemp = indx(k)
    indx(k) = indx(l + 1)
    indx(l + 1) = itemp
    if (func(xyzh(:,indx(l+1))) > func(xyzh(:,indx(ir)))) then
       itemp = indx(l + 1)
       indx(l + 1) = indx(ir)
       indx(ir) = itemp
    endif
    if (func(xyzh(:,indx(l))) > func(xyzh(:,indx(ir)))) then
       itemp = indx(l)
       indx(l) = indx(ir)
       indx(ir) = itemp
    endif
    if (func(xyzh(:,indx(l+1))) > func(xyzh(:,indx(l)))) then
       itemp = indx(l + 1)
       indx(l + 1) = indx(l)
       indx(l) = itemp
    endif
    i = l + 1
    j = ir
    indxt = indx(l)
    a = func(xyzh(:,indxt))

3   continue
    i = i + 1
    if (func(xyzh(:,indx(i))) < a) goto 3
4   continue
    j = j - 1
    if (func(xyzh(:,indx(j))) > a) goto 4
    if (j < i) goto 5
    itemp = indx(i)
    indx(i) = indx(j)
    indx(j) = itemp
    goto 3

5   indx(l) = indx(j)
    indx(j) = indxt
    jstack = jstack + 2
    if (jstack > nstack) then
       print*,'fatal error!!! stacksize exceeded in sort'
       print*,'need to set parameter nstack higher in subroutine indexxfunc '
       stop
    endif
    if (ir - i + 1 >= j - l) then
       istack(jstack) = ir
       istack(jstack - 1) = i
       ir = j - 1
    else
       istack(jstack) = j - 1
       istack(jstack - 1) = l
       l = i
    endif
 endif

 goto 1
end subroutine indexxfunc

!----------------------------------------------------------------
!+
!  customised low-memory sorting routine using Quicksort
!  sort key value on-the-fly by calling the function func
!  which can be any function of the particle positions.
!  (Tweaked version of the original one to sort a list of
!   neighbours founded using the KD tree)
!+
!----------------------------------------------------------------
subroutine Knnfunc(n, xpos, xyzh, indx)
 integer, intent(in)  :: n
 real,    intent(in)  :: xpos(3)
 real,    intent(in)  :: xyzh(:,:)
 integer, intent(out) :: indx(n)
 integer, parameter :: m=7, nstack=500
 real, external :: func

 integer :: i,j,k,l,ir,jstack,indxt,itemp
 integer :: istack(nstack)
 real :: a

 jstack = 0
 l = 1
 ir = n

1 if (ir - l < m) then
    do j = l + 1, ir
       indxt = indx(j)
       a = r2func_pos(xyzh(:,indxt),xpos)
       do i = j - 1, 1, -1
          if (r2func_pos(xyzh(:,indx(i)),xpos) <= a) goto 2
          indx(i + 1) = indx(i)
       enddo
       i = 0
2      indx(i + 1) = indxt
    enddo
    if (jstack==0) return
    ir = istack(jstack)
    l = istack(jstack - 1)
    jstack = jstack - 2
 else
    k = (l + ir)/2
    itemp = indx(k)
    indx(k) = indx(l + 1)
    indx(l + 1) = itemp
    if (r2func_pos(xyzh(:,indx(l+1)),xpos) > r2func_pos(xyzh(:,indx(ir)),xpos)) then
       itemp = indx(l + 1)
       indx(l + 1) = indx(ir)
       indx(ir) = itemp
    endif
    if (r2func_pos(xyzh(:,indx(l)),xpos) > r2func_pos(xyzh(:,indx(ir)),xpos)) then
       itemp = indx(l)
       indx(l) = indx(ir)
       indx(ir) = itemp
    endif
    if (r2func_pos(xyzh(:,indx(l+1)),xpos) > r2func_pos(xyzh(:,indx(l)),xpos)) then
       itemp = indx(l + 1)
       indx(l + 1) = indx(l)
       indx(l) = itemp
    endif
    i = l + 1
    j = ir
    indxt = indx(l)
    a = r2func_pos(xyzh(:,indxt),xpos)

3   continue
    i = i + 1
    if (r2func_pos(xyzh(:,indx(i)),xpos) < a) goto 3
4   continue
    j = j - 1
    if (r2func_pos(xyzh(:,indx(j)),xpos) > a) goto 4
    if (j < i) goto 5
    itemp = indx(i)
    indx(i) = indx(j)
    indx(j) = itemp
    goto 3

5   indx(l) = indx(j)
    indx(j) = indxt
    jstack = jstack + 2
    if (jstack > nstack) then
       print*,'fatal error!!! stacksize exceeded in sort'
       print*,'need to set parameter nstack higher in subroutine indexxfunc '
       stop
    endif
    if (ir - i + 1 >= j - l) then
       istack(jstack) = ir
       istack(jstack - 1) = i
       ir = j - 1
    else
       istack(jstack) = j - 1
       istack(jstack - 1) = l
       l = i
    endif
 endif

 goto 1
end subroutine Knnfunc

!----------------------------------------------------------------
!+
!  customised low-memory sorting routine using Quicksort
!  sort key value on-the-fly by calling the function func
!  which can be any function of the particle positions.
!  (Tweaked version of the original one to sort a list of
!   neighbours founded using the KD tree) (Parallel scheme, approx 2 times faster)
!+
!----------------------------------------------------------------
subroutine parqsort(n, arr,func, indx)
!$ use omp_lib,only:omp_get_num_threads
 implicit none
 integer, intent(in)    :: n
 real,    intent(in)    :: arr(n)
 integer, intent(inout) :: indx(n)
 integer, parameter :: m=8, nstack=500
 real, external :: func
 integer       :: i,j,k,il,ir,jstack,jqueue,indxt,itemp,nthreads,t,spt,nquick
 integer, save :: istack(nstack)
 !$omp threadprivate(istack)
 integer       :: iqueue(nstack)
 real :: a

 nthreads = 1

 !$omp parallel default(none) shared(nthreads)
!$ nthreads = omp_get_num_threads()
 !$omp end parallel

 spt = n/nthreads

 jstack = 0
 jqueue = 0
 iqueue = 0
 istack = 0
 il = 1
 ir = n

 do while (.true.)

    if (ir - il <= spt) then
       jqueue = jqueue + 2
       iqueue(jqueue) = ir
       iqueue(jqueue - 1) = il
       if (jstack==0) exit
       ir = istack(jstack)
       il = istack(jstack - 1)
       jstack = jstack - 2
    else
       k = (il + ir)/2
       i = il
       j = ir
       indxt = indx(k)
       a = func(arr(indxt))

       do while (j>i)
          do while(func(arr(indx(i))) < a)
             i = i + 1
          enddo
          do while (func(arr(indx(j))) > a)
             j = j - 1
          enddo
          if (j>i) then
             itemp = indx(i)
             indx(i) = indx(j)
             indx(j) = itemp
          endif
       enddo
       jstack = jstack + 2
       if (jstack > nstack) then
          print*,'fatal error!!! stacksize exceeded in sort'
          print*,'need to set parameter nstack higher in subroutine indexx '
          stop
       endif
       if (ir - i + 1 >= j - il) then
          istack(jstack) = ir
          istack(jstack - 1) = i
          ir = j - 1
       else
          istack(jstack) = j - 1
          istack(jstack - 1) = il
          il = i
       endif
    endif
 enddo

 istack = 0
 nquick = jqueue/2

 !$omp parallel do default(none) &
 !$omp shared(indx,arr,nquick,iqueue)&
 !$omp private(i,j,k,il,ir,a,jstack,indxt,itemp)
 do t=1,nquick
    ir = iqueue(2*t)
    il  = iqueue(2*t - 1)
    jstack = 0

    do while (.true.)
       if (ir - il < m) then
          !print*,il,ir
          do j = il , ir
             indxt = indx(j)
             a = func(arr(indxt))
             do i = j - 1, il, -1
                if (func(arr(indx(i))) <= a) goto 5
                indx(i + 1) = indx(i)
             enddo
             i = il-1
5            indx(i + 1) = indxt
          enddo
          if (jstack==0) exit
          ir = istack(jstack)
          il = istack(jstack - 1)
          jstack = jstack - 2
       else
          k = (il + ir)/2
          i = il
          j = ir
          indxt = indx(k)
          a = func(arr(indxt))

          do while (j>i)
             do while(func(arr(indx(i))) < a)
                i = i + 1
             enddo
             do while (func(arr(indx(j))) > a)
                j = j - 1
             enddo
             if (j>i) then
                itemp = indx(i)
                indx(i) = indx(j)
                indx(j) = itemp
             endif
          enddo
          jstack = jstack + 2
          if (jstack > nstack) then
             print*,'fatal error!!! stacksize exceeded in sort'
             print*,'need to set parameter nstack higher in subroutine indexx '
             stop
          endif
          if (ir - i + 1 >= j - il) then
             istack(jstack) = ir
             istack(jstack - 1) = i
             ir = j - 1
          else
             istack(jstack) = j - 1
             istack(jstack - 1) = il
             il = i
          endif
       endif
    enddo
 enddo

 end subroutine parqsort

!----------------------------------------------------------------
!+
!  Parallel LSD radix sort of 64-bit integer keys, carrying an
!  integer payload array (e.g. a permutation index).
!
!  Sorts key(1:n) in ascending order; indx(1:n) is permuted
!  identically. keybuf/indxbuf are caller-provided scratch arrays
!  of size >= n. The sort is stable and deterministic. npass 8-bit
!  digits are sorted (up to 8, covering keys < 2**64).
!+
!----------------------------------------------------------------
subroutine radixsort_i8(n,key,indx,keybuf,indxbuf,npass)
!$ use omp_lib, only: omp_get_num_threads
 integer,         intent(in)    :: n,npass
 integer(kind=8), intent(inout) :: key(n),keybuf(n)
 integer,         intent(inout) :: indx(n),indxbuf(n)
 integer, allocatable :: cnt(:,:),off(:,:)
 integer :: nthreads,ipass,i,np
 integer, parameter :: nradix = 256

 if (n <= 1) return
 np = min(max(npass,1),8)

 nthreads = 1
 !$omp parallel default(none) shared(nthreads)
 !$ nthreads = omp_get_num_threads()
 !$omp end parallel
 allocate(cnt(0:nradix-1,nthreads),off(0:nradix-1,nthreads))

 do ipass = 0,np-1
    if (mod(ipass,2) == 0) then
       call radix_pass_i8(n,key,indx,keybuf,indxbuf,ipass,nthreads,cnt,off)
    else
       call radix_pass_i8(n,keybuf,indxbuf,key,indx,ipass,nthreads,cnt,off)
    endif
 enddo

 ! odd number of passes: sorted data is in the scratch buffers
 if (mod(np,2) == 1) then
    !$omp parallel do default(none) shared(n,key,indx,keybuf,indxbuf) private(i)
    do i = 1,n
       key(i)  = keybuf(i)
       indx(i) = indxbuf(i)
    enddo
 endif
 deallocate(cnt,off)

end subroutine radixsort_i8

!----------------------------------------------------------------
!+
!  One digit pass of radixsort_i8: stable count-and-scatter of the
!  8-bit digit ipass of ksrc into kdst, carrying isrc into idst.
!+
!----------------------------------------------------------------
subroutine radix_pass_i8(n,ksrc,isrc,kdst,idst,ipass,nthreads,cnt,off)
!$ use omp_lib, only: omp_get_thread_num
 integer,         intent(in)    :: n,ipass,nthreads
 integer(kind=8), intent(in)    :: ksrc(n)
 integer,         intent(in)    :: isrc(n)
 integer(kind=8), intent(out)   :: kdst(n)
 integer,         intent(out)   :: idst(n)
 integer,         intent(inout) :: cnt(0:255,nthreads),off(0:255,nthreads)
 integer :: myoff(0:255)
 integer :: tid,t,istart,iend,pos,dg,i
 integer(kind=8) :: kval
 integer :: ival
 integer, parameter :: nradix = 256, nbits = 8

 !$omp parallel default(none) &
 !$omp shared(n,nthreads,ksrc,isrc,kdst,idst,ipass,cnt,off) &
 !$omp private(tid,t,istart,iend,dg,pos,i,kval,ival,myoff)
 tid = 1
 !$ tid = omp_get_thread_num() + 1
 istart = (tid-1)*n/nthreads + 1
 iend   = tid*n/nthreads

 ! count digit occurrences in this thread's chunk
 cnt(:,tid) = 0
 if (iend >= istart) then
    do i = istart,iend
       dg = int(ibits(ksrc(i),nbits*ipass,nbits))
       cnt(dg,tid) = cnt(dg,tid) + 1
    enddo
 endif
 !$omp barrier
 ! prefix sums over digits and threads -> global write offsets
 !$omp single
 pos = 1
 do dg = 0,nradix-1
    do t = 1,nthreads
       off(dg,t) = pos
       pos = pos + cnt(dg,t)
    enddo
 enddo
 !$omp end single
 ! stable scatter of this thread's chunk
 myoff(:) = off(:,tid)
 if (iend >= istart) then
    do i = istart,iend
       kval = ksrc(i)
       ival = isrc(i)
       dg = int(ibits(kval,nbits*ipass,nbits))
       pos = myoff(dg)
       kdst(pos) = kval
       idst(pos) = ival
       myoff(dg) = pos + 1
    enddo
 endif
 !$omp end parallel

end subroutine radix_pass_i8

!----------------------------------------------------------------
!+
!  Same as indexxfunc, except two particles can have the same
!  order/rank in the array ranki
!+
!----------------------------------------------------------------
subroutine find_rank(npart,func,xyzh,ranki)
 real,                 intent(in)  :: xyzh(:,:)
 integer,              intent(in)  :: npart
 integer, allocatable, intent(out) :: ranki(:)
 real, external :: func
 integer, allocatable :: iorder(:)
 real, parameter :: min_diff = tiny(1.)
 integer :: i,j,k

 ! First call indexxfunc
 allocate(iorder(npart),ranki(npart))
 call indexxfunc(npart,func,xyzh,iorder)
 ranki(iorder(1)) = 1

 do i=2,npart ! Loop over ranks sorted by indexxfunc
    j = iorder(i)
    k = iorder(i-1)
    if (abs(func(xyzh(:,j)) - func(xyzh(:,k))) > min_diff) then ! If particles have distinct radii
       ranki(j) = i
    else
       ranki(j) = ranki(k)       ! Else, give same ranks
    endif
 enddo

end subroutine find_rank

!----------------------------------------------------------------
!+
!  simplified interface to sort by radius given 3D cartesian
!  coordinates as input
!+
!----------------------------------------------------------------
subroutine sort_by_radius(n,xyzh,iorder,x0)
 integer, intent(in)  :: n
 real,    intent(in)  :: xyzh(4,n)
 integer, intent(out) :: iorder(n)
 real,    intent(in), optional :: x0(3)

 ! optional argument x0=[1,1,1] to set the origin
 if (present(x0)) then
    call set_r2func_origin(x0(1),x0(2),x0(3))
    call indexxfunc(n,r2func_origin,xyzh,iorder)
 else
    ! sort by r^2 using the r2func function
    call indexxfunc(n,r2func,xyzh,iorder)
 endif

end subroutine sort_by_radius

end module sortutils
