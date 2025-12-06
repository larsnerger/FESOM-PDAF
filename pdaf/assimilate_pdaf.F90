!>  Routine to call PDAF for analysis step
!!
!! This routine is called during the model integrations at each timestep.
!! It calls PDAF to check whether the forecast phase is completed and if
!! so, PDAF will perform the analysis step.
!!
!! __Revision history:__
!! * 2013-08 - Lars Nerger - Initial code for AWI-CM
!! * 2025-12 - Lars Nerger - Revision for PDAF3
!!
subroutine assimilate_pdaf(istep)

  use PDAF, &
       only: PDAF3_assimilate_local, PDAF_get_assim_flag
  use mod_parallel_pdaf, &        ! Parallelization variables
       only: mype_world, abort_parallel, task_id, mype_submodel, &
       COMM_COUPLE, filterpe
  use mod_assim_pdaf, &           ! Variables for assimilation
       only: filtertype, istep_asml, step_null, timemean, &
       dim_state_p, delt_obs_ocn, dim_ens, timemean_s, &
       monthly_state_sm, monthly_state_m, &
       compute_monthly_mm, compute_monthly_sm
  use fesom_pdaf, &
       only: timenew, daynew, yearnew, month, &
       num_day_in_month, fleapyear, daily_event, monthly_event
  use output_config_pdaf, &
       only: w_mm, w_sm, w_dayensm, w_monensm
  use mod_nc_out_routines, &
       only: netCDF_out
  use timer, &
       only: timeit

  implicit none
  include 'mpif.h'

! *** Arguments ***
  integer, intent(in) :: istep       !< current time step

! *** Local variables ***
  integer :: status_pdaf             ! PDAF status flag
  real, allocatable :: state_p(:)    ! Ensemble member / mean state
  real, allocatable :: stdev_p(:)    ! Standard deviation
  real, allocatable :: ensm_p(:)     ! Ensemble mean state
  integer :: mpierror
  real :: invdim_ens                 ! Inverse ensemble size
  real :: weights
  integer            :: assim_flag
  
  logical :: IsLastStepDay
  logical :: IsLastStepMonth
  
  integer, parameter :: int0 = 0

  ! External subroutines
  external :: collect_state_pdaf, &  ! Routine to collect a state vector from model fields
       distribute_state_pdaf, &      ! Routine to distribute a state vector to model fields
       next_observation_pdaf, &      ! Provide time step of next observation
       prepoststep_pdaf              ! User supplied pre/poststep routine
  ! Localization of state vector
  external :: init_n_domains_pdaf, & ! Provide number of local analysis domains
       init_dim_l_pdaf, &            ! Initialize state dimension for local analysis domain
       g2l_state_pdaf, &             ! Get state on local analysis domain from global state
       l2g_state_pdaf                ! Update global state from state on local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  external :: &
       init_dim_obs_pdafomi, &       ! Get dimension of full obs. vector for PE-local domain
       obs_op_pdafomi, &             ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_l_pdafomi        ! Get dimension of obs. vector for local analysis domain
  ! Subroutines used for generating observations
  external :: get_obs_f_pdaf         ! Get vector of synthetic observations from PDAF


! *********************************
! *** Call assimilation routine ***
! *********************************

  call timeit(6, 'new')

  istep_asml = istep + step_null  ! istep:       starting at 1 at each model (re)start
                                  ! istep_asml:  starting at 1 at beginning of each year

  if(mype_submodel==0 .and. task_id==1) write (*,'(a,1x,a,1x,a,1x,i5,1x,a,1x,i5,1x,a,1x,i3,1x,a,1x,i2,a,1x,i2,a)') &
          'FESOM-PDAF','assimilate_pdaf','istep', istep, 'istep_asml', istep_asml, 'day', daynew, &
          'time', floor(timenew/3600.0),'h',int(mod(timenew,3600.0)/60.0),'min'

  ! Call universal assimilate routine
  call PDAF3_assimilate_local(collect_state_pdaf, distribute_state_pdaf, &
       init_dim_obs_pdafomi, obs_op_pdafomi, &
       init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
       g2l_state_pdaf, l2g_state_pdaf, &
       prepoststep_pdaf, next_observation_pdaf, status_pdaf)

  ! Check for errors during execution of PDAF
  if (status_pdaf /= 0) then
     write (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in PDAF_put_state - stopping! (PE ', mype_world,')'
     call  abort_parallel()
  end if
  

  call timeit(6, 'old')
  call timeit(7, 'new')

  ! *********************************
  ! *** Compute daily mean        ***
  ! *********************************
  ! daily means ("m"-state) are averaged over 1 analysis step followed by step_per_day-minus-1 model forecast steps

  ! note: computing ensemble mean of state fields at every step is less efficient,
  ! but required to compute ensemble standard deviation at every step
     
  if (w_mm .or. w_sm) then

     call daily_event  (IsLastStepDay,  1)
     call monthly_event(IsLastStepMonth,1)

     if ( .not. allocated(state_p))            allocate(state_p(dim_state_p))
     if ( .not. allocated(ensm_p ))            allocate(ensm_p (dim_state_p))
     if ( .not. allocated(stdev_p) .and. w_sm) allocate(stdev_p(dim_state_p))

     call PDAF_get_assim_flag(assim_flag)
     
     ! *** in between assimilation steps, add forecast steps to m-fields ***
     ! note: assimilation step is at first time step of day
     if (assim_flag == 0) then
        ! collect instantenous model data
        call collect_state_pdaf(dim_state_p, state_p)
        ! compute ensemble mean
        invdim_ens = 1.0 / real(dim_ens)
        call MPI_ALLREDUCE((state_p*invdim_ens),ensm_p,dim_state_p,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_COUPLE,mpierror)
        ! compute ensemble mean of squared deviations on filterpe
        if (w_sm) call MPI_REDUCE(((ensm_p-state_p)*(ensm_p-state_p)*invdim_ens),stdev_p,dim_state_p,&
             MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM_COUPLE,mpierror)
        ! add to daily mean on filterpe
        if (filterpe) then
           timemean    = timemean    + ensm_p / delt_obs_ocn
           ! add daily mean of standard deviation
           if (w_sm) stdev_p     = sqrt(invdim_ens * stdev_p)
           if (w_sm) timemean_s  = timemean_s  + stdev_p / delt_obs_ocn
        endif ! filterpe
     endif ! assim_flag
     
     ! *** compute monthly means ***
     if (filterpe) then
     if (IsLastStepDay) then
        ! add m-fields to monthly mean
        if (compute_monthly_mm) monthly_state_m  = monthly_state_m  + timemean
        if (compute_monthly_sm) monthly_state_sm = monthly_state_sm + timemean_s
     endif
     if (IsLastStepMonth) then
        ! compute monthly mean
        weights = 1.0/real(num_day_in_month(fleapyear,month))
        if (compute_monthly_mm) monthly_state_m  = monthly_state_m  * weights
        if (compute_monthly_sm) monthly_state_sm = monthly_state_sm * weights
     endif
     
     ! *** write output and reset to zero ***
     if (IsLastStepMonth) then
        ! monthly and daily output
        if (w_dayensm .or. w_monensm) call netCDF_out('mm',timemean  , int0, IsLastStepMonth, m_state_p=monthly_state_m )
        if (w_dayensm .or. w_monensm) call netCDF_out('sm',timemean_s, int0, IsLastStepMonth, m_state_p=monthly_state_sm)
        ! reset monthly and daily
        timemean = 0
        if (w_sm) timemean_s = 0
        if (compute_monthly_mm) monthly_state_m = 0
        if (compute_monthly_sm) monthly_state_sm = 0
     elseif (IsLastStepDay) then
        ! daily output
        if (w_dayensm)            call netCDF_out('mm',timemean,   int0, IsLastStepMonth)
        if (w_dayensm .and. w_sm) call netCDF_out('sm',timemean_s, int0, IsLastStepMonth)
        ! reset daily
        timemean = 0
        if (w_sm) timemean_s = 0
     endif
     endif ! filterpe
  endif ! w_mm

  call timeit(7, 'old')

end subroutine assimilate_pdaf
