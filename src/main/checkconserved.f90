!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module checkconserved
!
! Utility routines to perform runtime checks that
! conservation laws are appropriately satisfied
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary_dyn, dim, externalforces, io, options, part
!
 use dim, only:maxdusttypes,use_apr
 implicit none
 real, public :: get_conserv = 1.0 ! to track when we have initial values for conservation laws
 real, public :: etot_in,angtot_in,totmom_in,mdust_in(maxdusttypes),mtot_in

 public :: init_conservation_checks, check_conservation_error
 public :: check_magnetic_stability

 private

contains

!----------------------------------------------------------------
!+
!  check if conservation of various properties *should* be
!  possible given the range of physics selected
!+
!----------------------------------------------------------------
subroutine init_conservation_checks(should_conserve_energy,should_conserve_momentum,&
                                    should_conserve_angmom,should_conserve_dustmass,&
                                    should_conserve_aprmass)
 use options,     only:icooling,ieos,ipdv_heating,ishock_heating,&
                       iresistive_heating,use_dustfrac,iexternalforce
 use dim,         only:mhd,maxvxyzu,periodic,inject_parts,use_apr
 use part,        only:iboundary,npartoftype
 use boundary_dyn,only:dynamic_bdy
 logical, intent(out) :: should_conserve_energy,should_conserve_momentum
 logical, intent(out) :: should_conserve_angmom,should_conserve_dustmass
 logical, intent(out) :: should_conserve_aprmass

 !
 ! should conserve energy if using adiabatic equation of state with no cooling
 ! as long as all heating terms are included
 !
 should_conserve_energy = (maxvxyzu==4 .and. ieos==2 .and. &
                          icooling==0 .and. ipdv_heating==1 .and. ishock_heating==1 &
                          .and. (.not.mhd .or. iresistive_heating==1))
 !
 ! code should conserve momentum unless boundary particles are employed
 !
 if (iexternalforce/=0) then
    should_conserve_momentum = .false.
 else
    should_conserve_momentum = (npartoftype(iboundary)==0)
 endif
 !
 ! code should conserve angular momentum as long as no boundaries (fixed or periodic)
 ! and as long as there are no non-radial forces (iexternalforce > 1)
 !
 should_conserve_angmom = (npartoftype(iboundary)==0 .and. .not.periodic &
                          .and. iexternalforce <= 1)
 !
 ! should always conserve dust mass
 !
 should_conserve_dustmass = use_dustfrac
 !
 ! Each injection routine will need to bookeep conserved quantities, but until then...
 !
 if (inject_parts .or. dynamic_bdy .or. use_apr) then
    should_conserve_energy   = .false.
    should_conserve_momentum = .false.
    should_conserve_angmom   = .false.
 endif

 ! This is to check that total mass is conserved when we use apr
 ! It can't be used if mass is accreted or injected
 should_conserve_aprmass = (iexternalforce==0 .and. use_apr .and. .not.inject_parts)

end subroutine init_conservation_checks

!----------------------------------------------------------------
!+
!  routine to check conservation errors during the calculation
!  and stop if it is too large
!+
!----------------------------------------------------------------
subroutine check_conservation_error(val,ref,tol,label,decrease)
 use io,             only:error,fatal,iverbose
 use options,        only:iexternalforce
 use externalforces, only:iext_corot_binary
 real, intent(in) :: val,ref,tol
 character(len=*), intent(in) :: label
 logical, intent(in), optional :: decrease
 real :: err

 if (abs(ref) > 1.e-3) then
    err = (val - ref)/abs(ref)
 else
    err = (val - ref)
 endif
 if (present(decrease)) then
    err = max(err,0.) ! allow decrease but not increase
 else
    err = abs(err)
 endif
 if (err > tol) then
    if ((trim(label) == 'angular momentum' .or. trim(label) == 'energy') &
        .and. iexternalforce == iext_corot_binary) then
       call error('evolve',trim(label)//' is not being conserved due to corotating frame',var='err',val=err)
    else
       call error('evolve','Large error in '//trim(label)//' conservation ',var='err',val=err)
       call do_not_publish_crap('evolve','Conservation errors too large to continue simulation')
    endif
 else
    if (iverbose >= 2) print "(a,es10.3)",trim(label)//' error is ',err
 endif

end subroutine check_conservation_error
!----------------------------------------------------------------
!+
!  routine to check the stability of the magnetic field based upon
!  the values of h |divB|/B
!  Although not at true conservation check, this is a stability check
!  so is related to the checks performed here
!+
!----------------------------------------------------------------
subroutine check_magnetic_stability(hdivBonB_ave,hdivBonB_max)
 use io,      only:fatal
 real, intent(in) :: hdivBonB_ave,hdivBonB_max

 if (hdivBonB_max > 100 .or. hdivBonB_ave > 0.1) then
    ! Tricco, Price & Bate (2016) suggest the average should remain lower than 0.01,
    ! but we will increase it here due to the nature of the exiting the code
    ! The suggestion of 512 was empirically determined in Dobbs & Wurster (2021)
    call do_not_publish_crap('evolve','h|divb|/b is too large; recommend to increase the overcleanfac')
 endif

end subroutine check_magnetic_stability

subroutine do_not_publish_crap(subr,msg)
 use io, only:fatal
 character(len=*), intent(in) :: subr,msg
 character(len=20) :: string

 call get_environment_variable('I_WILL_NOT_PUBLISH_CRAP',string)
 if (.not. (trim(string)=='yes')) then
    print "(2(/,a))",' You can ignore this error and continue by setting the ',&
                     ' environment variable I_WILL_NOT_PUBLISH_CRAP=yes to continue'
    call fatal(subr,msg)
 endif

end subroutine do_not_publish_crap
!----------------------------------------------------------------
end module checkconserved
