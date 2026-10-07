!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module deriv
!
! this module is a wrapper for the main derivative evaluation
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: HIIRegion, cons2prim, densityforce, derivutils, dim,
!   externalforces, forces, forcing, growth, io, metric_tools, neighkdtree,
!   options, part, porosity, ptmass, ptmass_radiation, radiation_implicit,
!   timestep, timestep_ind, timing
!
 implicit none

 public :: derivs, get_derivs_global, get_density_global
 real, private :: stressmax

 private

contains

!-------------------------------------------------------------
!+
!  calculates derivatives of all particle quantities
!  (wrapper for call to density and rates, calls neighbours etc first)
!+
!-------------------------------------------------------------
subroutine derivs(icall,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                  Bevol,dBevol,rad,drad,radprop,dustprop,ddustprop,&
                  dustevol,ddustevol,filfac,dustfrac,eos_vars,time,dt,dtnew,pxyzu,&
                  dens,metrics,apr_level)
 use dim,            only:mhd,fast_divcurlB,gr,periodic,do_radiation,driving,&
                          sink_radiation,use_dustgrowth,ind_timesteps,isothermal,mpi,gravity
 use io,             only:iprint,fatal,error
 use neighkdtree,    only:build_tree,refine_local_tree
 use mpighosts,      only:refresh_tree_ghosts,ighost_dens,ighost_force
 use densityforce,   only:densityiterate
 use ptmass,         only:ipart_rhomax,ptmass_calc_enclosed_mass,ptmass_boundary_crossing,get_pressure_on_sinks
 use externalforces, only:externalforce
 use part,           only:dustgasprop,Vrel_disp,dvdx,Bxyz,set_boundaries_to_active,&
                          nptmass,xyzmh_ptmass,sinks_have_heating,dust_temp,VrelVf,fxyz_drag,rho
 use timestep_ind,   only:nbinmax
 use timestep,       only:dtmax,dtcourant,dtforce,dtrad
 use forcing,        only:forceit
 use growth,           only:get_growth_rate
 use porosity,         only:get_disruption,get_probastick
 use ptmass_radiation, only:get_dust_temperature
 use timing,         only:get_timings
 use forces,         only:force,dualwalk_global_force,clear_ghosts
 use io,             only:nprocs
 use part,           only:mhd,gradh,alphaind,iradxi,ifluxx,ifluxy,ifluxz,ithick
 use derivutils,     only:do_timing
 use cons2prim,      only:cons2primall,cons2prim_everything
 use metric_tools,   only:init_metric
 use radiation_implicit, only:do_radiation_implicit,ierr_failed_to_converge
 use options,        only:implicit_radiation,implicit_radiation_store_drad,use_porosity,need_pressure_on_sinks
 use HIIRegion,      only:HIIupdateflag,iH2R,HII_feedback
 integer,         intent(in)    :: icall
 integer,         intent(inout) :: npart
 integer,         intent(in)    :: nactive
 real,            intent(inout) :: xyzh(:,:)
 real,            intent(inout) :: vxyzu(:,:)
 real,            intent(inout) :: fxyzu(:,:)
 real,            intent(inout) :: fext(:,:)       ! inout: ghost particles written after npart (MPI)
 real(kind=4),    intent(out)   :: divcurlv(:,:)
 real(kind=4),    intent(out)   :: divcurlB(:,:)
 real,            intent(inout) :: Bevol(:,:)     ! inout: ghost particles written after npart (MPI)
 real,            intent(out)   :: dBevol(:,:)
 real,            intent(inout) :: rad(:,:)
 real,            intent(out)   :: eos_vars(:,:)
 real,            intent(out)   :: drad(:,:)
 real,            intent(inout) :: radprop(:,:)
 real,            intent(in)    :: dustevol(:,:)
 real,            intent(inout) :: dustprop(:,:)
 real,            intent(out)   :: dustfrac(:,:)
 real,            intent(out)   :: ddustevol(:,:),ddustprop(:,:)
 real,            intent(inout) :: filfac(:)
 real,            intent(in)    :: time,dt
 real,            intent(out)   :: dtnew
 real,            intent(inout) :: pxyzu(:,:), dens(:)
 real,            intent(inout) :: metrics(:,:,:,:)
 integer(kind=1), intent(inout) :: apr_level(:)   ! inout: ghost particles written after npart (MPI)
 integer                     :: ierr,i
 real(kind=4)                :: t1,tcpu1,tlast,tcpulast
 logical                     :: redo_ghosts

 t1    = 0.
 tcpu1 = 0.
 call get_timings(t1,tcpu1)
 tlast    = t1
 tcpulast = tcpu1
!
!--check for errors in input options
!
 if (icall < 0 .or. icall > 2) call fatal('deriv','invalid icall on input')
!
! icall is a flag to say whether or not positions have changed
! since the last call to derivs.
!
! icall = 1 is the "standard" call to derivs: calculates all derivatives
! icall = 2 does not remake the tree build and does not recalculate density
!           (ie. only re-evaluates the SPH force term using updated values
!            of the input variables)
!
! build tree to prepare neighbour finding
!
 if (icall==1 .or. icall==0) then
    ! with MPI: the domains, then the local tree with the ghost particles read in density
    call build_tree(npart,nactive,xyzh,vxyzu,domains_only=.true.)
    if (mpi .and. nprocs > 1) &
       call build_ghosts_tree(ighost_dens,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,&
                              rad,radprop,dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level)

    if (gr) then
       ! update time-dependent metric (e.g. binary BH) and repack at particle positions
       call init_metric(npart,xyzh,metrics,time=time)
    endif

    if (nptmass > 0 .and. periodic) call ptmass_boundary_crossing(nptmass,xyzmh_ptmass)
 endif

 call do_timing('tree',tlast,tcpulast,start=.true.)

 !
 ! compute disruption of dust particles
 !
 if (use_dustgrowth .and. use_porosity) call get_disruption(npart,xyzh,filfac,dustprop,dustgasprop,rho)
!
! calculate density by direct summation
!

 if (icall==1) then
    redo_ghosts = .true.
    do while(redo_ghosts)
       call densityiterate(1,npart,nactive,xyzh,vxyzu,divcurlv,divcurlB,Bevol,&
                           stressmax,fxyzu,fext,alphaind,gradh,rad,radprop,dvdx,apr_level,redo_ghosts)
       ! MPI: h grew beyond the ghost particles, choose them again with the new h
       if (redo_ghosts) call build_ghosts_tree(ighost_dens,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,&
                                         rad,radprop,dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level)
    enddo
    if (.not. fast_divcurlB) then
       ! Repeat the call to calculate all the non-density-related quantities in densityiterate.
       ! This needs to be separate for an accurate calculation of divcurlB which requires an up-to-date rho.
       ! if fast_divcurlB = .false., then all additional quantities are calculated during the previous call
       ! (MPI: the ghost particles with their new rho)
       if (mpi .and. nprocs > 1) &
          call refresh_tree_ghosts(ighost_dens,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,&
                                   rad,radprop,dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level)
       call densityiterate(3,npart,nactive,xyzh,vxyzu,divcurlv,divcurlB,Bevol,&
                           stressmax,fxyzu,fext,alphaind,gradh,rad,radprop,dvdx,apr_level)
       ! put a similar flag for pressure calculation from dens: call cons2primall/everyhting and densityiterate(3. Import pressure from eos_vars in dens and use it to calculate delta_v
    endif
    set_boundaries_to_active = .false.     ! boundary particles are no longer treated as active
    call do_timing('dens',tlast,tcpulast)
 endif
!
!-- update ionising state of the particle if HII regions are used in cluster formation simulations
!
 if (iH2R >0) then
    if (HIIupdateflag) then
       call HII_feedback(nptmass,npart,xyzh,xyzmh_ptmass,vxyzu,rho,eos_vars)
       HIIupdateflag = .false.
    endif
    call do_timing('HII_region',tlast,tcpulast)
 endif

 if (gr) then
    call cons2primall(npart,xyzh,metrics,pxyzu,vxyzu,dens,eos_vars)
 else
    call cons2prim_everything(npart,xyzh,vxyzu,dvdx,rad,eos_vars,radprop,Bevol,Bxyz,dustevol,dustfrac,alphaind)
 endif
 call do_timing('cons2prim',tlast,tcpulast)

 !
 ! implicit radiation update
 !
 if (do_radiation .and. implicit_radiation .and. dt > 0.) then
    call do_radiation_implicit(dt,npart,rad,xyzh,vxyzu,radprop,drad,ierr)
    if (ierr /= 0 .and. ierr /= ierr_failed_to_converge) call fatal('radiation','Failed in radiation')
    call do_timing('radiation',tlast,tcpulast)
 endif

 !
 ! compute forces
 !
 if (driving) then
    ! forced turbulence -- call driving routine
    call forceit(time,npart,xyzh,vxyzu,fxyzu)
    call do_timing('driving',tlast,tcpulast)
 endif

 !
 ! compute SPH forces
 !
 stressmax = 0.
 if (sinks_have_heating(nptmass,xyzmh_ptmass)) call ptmass_calc_enclosed_mass(nptmass,npart,xyzh)
 if (mpi .and. nprocs > 1) then
    if (icall==0 .or. icall==1) then
       if (gravity) then
          ! local tree without the ghost particles, refined in the global tree
          call refine_local_tree(npart,xyzh)
       else
          ! ghost particles read in force, chosen with the new h, in the local tree
          call build_ghosts_tree(ighost_force,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,&
                                 rad,radprop,dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level)
       endif
    elseif (.not.gravity) then
       ! same ghost particles (neither positions nor h have changed): only their values
       call refresh_tree_ghosts(ighost_force,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,&
                                rad,radprop,dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level)
    endif
    ! dual tree walk over MPI: remote nodes and ghost particles for force
    if (gravity) call dualwalk_global_force(npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,&
                                            rad,radprop,dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level)
 endif
 call force(icall,npart,xyzh,vxyzu,fxyzu,divcurlv,divcurlB,Bevol,dBevol,&
            rad,drad,radprop,dustprop,dustgasprop,Vrel_disp,dustfrac,ddustevol,fext,fxyz_drag,&
            ipart_rhomax,dt,stressmax,eos_vars,dens,metrics,apr_level)
 if (mpi .and. nprocs > 1) call clear_ghosts(npart,xyzh)
 call do_timing('force',tlast,tcpulast)

 !
 ! compute growth rate of dust particles
 !
 if (use_dustgrowth) then
    call get_growth_rate(npart,xyzh,vxyzu,rho,dustgasprop,VrelVf,dustprop,filfac,ddustprop(1,:),Vrel_disp)!--we only get dm/dt (i.e 1st dimension of ddustprop)
    ! compute growth rate and probability of sticking/bouncing of porous dust
    if (use_porosity) call get_probastick(npart,xyzh,ddustprop(1,:),dustprop,dustgasprop,filfac)
 endif
!
! compute density and pressure at location of sink particles
!
 if (need_pressure_on_sinks) call get_pressure_on_sinks(nptmass,xyzmh_ptmass)
!
! compute dust temperature
!
 if (sink_radiation .and. .not.isothermal) then
    call get_dust_temperature(npart,xyzh,eos_vars,nptmass,xyzmh_ptmass,dust_temp)
 endif

 if (do_radiation .and. implicit_radiation .and. .not.implicit_radiation_store_drad) then
    !$omp parallel do shared(drad,fxyzu,npart) private(i)
    do i=1,npart
       drad(:,i) = 0.
       fxyzu(4,i) = 0.
    enddo
    !$omp end parallel do
 endif
!
! set new timestep from Courant/forces condition
!
 if (ind_timesteps) then
    dtnew = dtmax/2.**nbinmax  ! minimum timestep over all particles
 else
    dtnew = min(dtforce,dtcourant,dtrad,dtmax)
 endif

 call do_timing('total',t1,tcpu1,lunit=iprint)

end subroutine derivs

!-------------------------------------------------------------
!+
!  MPI: choose the ghost particles of the other tasks (the fields
!  read for a neighbour in density or in force, iset) and build
!  the local tree with them
!+
!-------------------------------------------------------------
subroutine build_ghosts_tree(iset,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,&
                             rad,radprop,dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level)
 use neighkdtree, only:build_ghost_tree
 use mpighosts,   only:exchange_tree_ghosts,nghost_tree
 integer,         intent(in)    :: iset,npart
 real,            intent(inout) :: xyzh(:,:),vxyzu(:,:),fxyzu(:,:),fext(:,:),Bevol(:,:),rad(:,:)
 real,            intent(inout) :: radprop(:,:),dustprop(:,:),dustfrac(:,:),filfac(:),eos_vars(:,:)
 real,            intent(inout) :: dens(:),metrics(:,:,:,:)
 real(kind=4),    intent(inout) :: divcurlv(:,:),divcurlB(:,:)
 integer(kind=1), intent(inout) :: apr_level(:)

 call exchange_tree_ghosts(iset,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,&
                           rad,radprop,dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level)
 call build_ghost_tree(npart,xyzh,nghost_tree)

end subroutine build_ghosts_tree

!--------------------------------------
!+
!  wrapper for the call to derivs
!  so only one line needs changing
!  if interface changes
!
!  this should NOT be called during timestepping, it is useful
!  for when one requires just a single call to evaluate derivatives
!  and store them in the global shared arrays
!+
!--------------------------------------
subroutine get_derivs_global(tused,dt_new,dt,icall)
 use part,         only:npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                        Bevol,dBevol,rad,drad,radprop,dustprop,ddustprop,filfac,&
                        dustfrac,ddustevol,eos_vars,pxyzu,dens,metrics,dustevol,gr,&
                        apr_level
 use timing,       only:printused,getused
 use io,           only:id,master
 use cons2prim,    only:prim2consall
 use metric_tools, only:init_metric
 real(kind=4), intent(out), optional :: tused
 real,         intent(out), optional :: dt_new
 real,         intent(in),  optional :: dt  ! optional argument needed to test implicit radiation routine
 integer,      intent(in),  optional :: icall
 real(kind=4) :: t1,t2
 real    :: dtnew,dti,time
 integer :: icalli

 time = 0.
 dti = 0.
 icalli = 1
 if (present(dt)) dti = dt
 if (present(icall)) icalli = icall
 call getused(t1)
 ! update conserved quantities in the GR code
 if (gr) then
    ! kernel-summed rho required for prim2consall (all MPI ranks must call this)
    call get_density_global(2,zero_fxyzu=.true.)
    call init_metric(npart,xyzh,metrics,time=time)
    call prim2consall(npart,xyzh,metrics,vxyzu,pxyzu,dens=dens)
 endif

 ! evaluate derivatives
 call derivs(icalli,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,&
             rad,drad,radprop,dustprop,ddustprop,dustevol,ddustevol,filfac,dustfrac,&
             eos_vars,time,dti,dtnew,pxyzu,dens,metrics,apr_level)

 call getused(t2)
 if (id==master .and. present(tused)) call printused(t1)
 if (present(tused)) tused = t2 - t1
 if (present(dt_new)) dt_new = dtnew

end subroutine get_derivs_global

!--------------------------------------
!+
!  wrapper for the call to densityiterate
!  so only one line needs changing
!  if interface changes
!
!  this should be used when one requires just a density calculation
!  and store results in the global shared arrays
!+
!--------------------------------------
subroutine get_density_global(icall,nactive,zero_fxyzu,make_tree)
 use part,         only:npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                        Bevol,alphaind,gradh,rad,radprop,dvdx,apr_level,&
                        dustprop,dustfrac,filfac,eos_vars,dens,metrics
 use dim,          only:mpi
 use io,           only:nprocs
 use densityforce, only:densityiterate
 use neighkdtree,  only:build_tree
 use mpighosts,    only:ighost_dens
 integer, intent(in) :: icall
 integer, intent(in), optional :: nactive
 logical, intent(in), optional :: zero_fxyzu
 logical, intent(in), optional :: make_tree
 integer :: nactivei
 logical :: do_tree,redo_ghosts
 real    :: stressmax

 nactivei = npart
 if (present(nactive)) nactivei = nactive

 do_tree = .true.
 if (present(make_tree)) do_tree = make_tree

 ! build tree to prepare neighbour finding (if requested)
 ! with MPI: the domains, then the local tree with the ghost particles read in density
 if (do_tree) then
    call build_tree(npart,nactivei,xyzh,vxyzu,domains_only=.true.)
    if (mpi .and. nprocs > 1) &
       call build_ghosts_tree(ighost_dens,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,&
                              rad,radprop,dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level)
 endif

 ! optionally zero fxyzu (useful for initialization)
 if (present(zero_fxyzu)) then
    if (zero_fxyzu) fxyzu = 0.
 endif

 ! evaluate density (MPI: again with new ghost particles if h grew beyond them)
 stressmax   = 0.
 redo_ghosts = .true.
 do while(redo_ghosts)
    call densityiterate(icall,npart,nactivei,xyzh,vxyzu,divcurlv,divcurlB,Bevol,stressmax,&
                        fxyzu,fext,alphaind,gradh,rad,radprop,dvdx,apr_level,redo_ghosts)
    if (redo_ghosts) call build_ghosts_tree(ighost_dens,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,&
                                        rad,radprop,dustprop,dustfrac,filfac,eos_vars,dens,metrics,apr_level)
 enddo

end subroutine get_density_global

end module deriv
