!>  Module for computing time means
!!
!! __Revision history:__
!! * 2025-12 - Lars Nerger - Revision for PDAF3, extracting code from other routines
!!
module means_pdaf

  implicit none
  public
  save

  real, allocatable :: timemean(:)     ! daily mean local state vector (mean of model forecast steps and analysis step)
  real, allocatable :: timemean_s(:)

  ! whether to compute monthly means:
  logical :: compute_monthly_ff = .false.
  logical :: compute_monthly_aa = .false.
  logical :: compute_monthly_mm = .false.
  logical :: compute_monthly_sf = .false.
  logical :: compute_monthly_sa = .false.
  logical :: compute_monthly_sm = .false.

  real, allocatable :: monthly_state_f(:)       ! forecasted monthly state
  real, allocatable :: monthly_state_a(:)       ! analyzed monthly state
  real, allocatable :: monthly_stddev_f(:)      ! forecasted monthly standard deviation
  real, allocatable :: monthly_stddev_a(:)      ! analyzed monthly standard deviation
  real, allocatable :: monthly_state_m(:)       ! monthly time-mean state
  real, allocatable :: monthly_stddev_m(:)      ! monthly time-mean standard deviation

contains

!-----------------------------------------------------------------------
!>  Routine to initialize time means 
!!
!! __Revision history:__
!! * 2025-12 - Lars Nerger - Revision for PDAF3 extracted from init_pdaf
!!
  subroutine init_means_pdaf(dim_p)

    use statevector_pdaf, &
         only: nfields, sfields
    use output_config_pdaf, &
         only: mm, aa, ff, ii, sa, sf, si, sm, oo, dd, ee

    implicit none

! *** Arguments ***
    integer, intent(in) :: dim_p

! *** Local variables ***
    integer :: s        ! counter

    allocate(timemean(dim_p))
    timemean = 0.0
  
    allocate(timemean_s(dim_p))
    timemean_s = 0.0

    do s=1, nfields
       !        -- output? -----------        -- ensemble mean? ----------            -- monthly? ----------------------
       ! have monthly forecast (ff) states to write?
       compute_monthly_ff = compute_monthly_ff .or. &
            ( sfields(s)%output(ff,oo) .and. (.not. sfields(s)%output(ff,ee)) .and. (.not. sfields(s)%output(ff,dd)))
       ! have monthly analysis (aa) states to write?
       compute_monthly_aa = compute_monthly_aa .or. &
            ( sfields(s)%output(aa,oo) .and. (.not. sfields(s)%output(aa,ee)) .and. (.not. sfields(s)%output(aa,dd)))
       ! have monthly daymean  (mm) states to write?
       compute_monthly_mm = compute_monthly_mm .or. &
            ( sfields(s)%output(mm,oo) .and. (.not. sfields(s)%output(mm,ee)) .and. (.not. sfields(s)%output(mm,dd)))
       ! have monthly forecast (sf) STD to write?
       compute_monthly_sf = compute_monthly_sf .or. &
            ( sfields(s)%output(sf,oo) .and. (.not. sfields(s)%output(sf,ee)) .and. (.not. sfields(s)%output(sf,dd)))
       ! have monthly analysis (sa) STD to write?
       compute_monthly_sa = compute_monthly_sa .or. &
            ( sfields(s)%output(sa,oo) .and. (.not. sfields(s)%output(sa,ee)) .and. (.not. sfields(s)%output(sa,dd)))
       ! have monthly daymean  (sm) STD to write?
       compute_monthly_sm = compute_monthly_sm .or. &
            ( sfields(s)%output(sm,oo) .and. (.not. sfields(s)%output(sm,ee)) .and. (.not. sfields(s)%output(sm,dd)))
    enddo

    if ((.not. allocated(monthly_state_a)))     allocate(monthly_state_a(dim_p))
    if ((.not. allocated(monthly_state_f)))     allocate(monthly_state_f(dim_p))
    if ((.not. allocated(monthly_state_m)))     allocate(monthly_state_m(dim_p))
    if ((.not. allocated(monthly_stddev_a)))    allocate(monthly_stddev_a(dim_p))
    if ((.not. allocated(monthly_stddev_f)))    allocate(monthly_stddev_f(dim_p))
    if ((.not. allocated(monthly_stddev_m)))    allocate(monthly_stddev_m(dim_p))

    ! init monthly state
    monthly_state_a=  0.0D0
    monthly_state_m=  0.0D0
    monthly_state_f=  0.0D0
    monthly_stddev_a= 0.0D0
    monthly_stddev_m= 0.0D0
    monthly_stddev_f= 0.0D0

  end subroutine init_means_pdaf


!-----------------------------------------------------------------------
!>  Routine to compute and output time means 
!!
!! This routine is called during the model integrations at each timestep
!! and enables to compute ensemble standard deviations at every step.
!!
!! __Revision history:__
!! * ~2022 - Frauke B - Initial code
!! * 2025-12 - Lars Nerger - Revision for PDAF3 extracted from assimilate_pdaf
!!
  subroutine do_means_pdaf()

    use mpi
    use PDAF, &
         only: PDAF_get_assim_flag
    use parallel_pdaf_mod, &
         only: COMM_COUPLE, filterpe
    use output_config_pdaf, &
         only: w_mm, w_sm, w_dayensm, w_monensm
    use mod_nc_out_routines, &
         only: netCDF_out
    use fesom_pdaf, &
         only: timenew, daynew, yearnew, month, &
         num_day_in_month, fleapyear, daily_event, monthly_event
    use assim_pdaf_mod, &           ! Variables for assimilation
         only: dim_state_p, delt_obs_ocn, dim_ens

    implicit none

! *** Local variables ***
    logical :: IsLastStepDay
    logical :: IsLastStepMonth
    integer, parameter :: int0 = 0
    integer            :: assim_flag
    real, allocatable :: state_p(:)    ! Ensemble member / mean state
    real, allocatable :: stdev_p(:)    ! Standard deviation
    real, allocatable :: ensm_p(:)     ! Ensemble mean state
    integer :: mpierror
    real :: invdim_ens                 ! Inverse ensemble size
    real :: weights


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
             timemean = timemean + ensm_p / delt_obs_ocn

             ! add daily mean of standard deviation
             if (w_sm) stdev_p = sqrt(invdim_ens * stdev_p)
             if (w_sm) timemean_s  = timemean_s  + stdev_p / delt_obs_ocn
          endif ! filterpe
       endif ! assim_flag

       ! *** compute monthly means ***
       if (filterpe) then

          if (IsLastStepDay) then
             ! add m-fields to monthly mean
             if (compute_monthly_mm) monthly_state_m  = monthly_state_m  + timemean
             if (compute_monthly_sm) monthly_stddev_m = monthly_stddev_m + timemean_s
          endif

          if (IsLastStepMonth) then
             ! compute monthly mean
             weights = 1.0/real(num_day_in_month(fleapyear,month))
             if (compute_monthly_mm) monthly_state_m  = monthly_state_m  * weights
             if (compute_monthly_sm) monthly_stddev_m = monthly_stddev_m * weights
          endif
     
          ! *** write output and reset to zero ***

          if (IsLastStepMonth) then
             ! monthly and daily output
             if (w_dayensm .or. w_monensm) call netCDF_out('mm',timemean  , int0, IsLastStepMonth, m_state_p=monthly_state_m )
             if (w_dayensm .or. w_monensm) call netCDF_out('sm',timemean_s, int0, IsLastStepMonth, m_state_p=monthly_stddev_m)

             ! reset monthly and daily
             timemean = 0
             if (w_sm) timemean_s = 0
             if (compute_monthly_mm) monthly_state_m = 0
             if (compute_monthly_sm) monthly_stddev_m = 0

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

  end subroutine do_means_pdaf


!-----------------------------------------------------------------------
!>  Routine to update time means 
!!
!! daily means ("m"-state) are averaged over one analysis step and the
!! consecutive model forecast steps of that day during analysis step, add to m-fields
!!
!! __Revision history:__
!! * 2025-12 - Lars Nerger - Revision for PDAF3 extracted from prepoststep_pdaf
!!
  subroutine update_means_pdaf(step, dim_p, state_p, stdev_p, delt_obs, write_now)

    use assim_pdaf_mod, &
         only: step_null
    use fesom_pdaf, &
         only: month, fleapyear, num_day_in_month
    use output_config_pdaf, &
         only: w_mm, w_sm

    implicit none

! *** Arguments ***
    integer, intent(in) :: step
    integer, intent(in) :: dim_p
    real, intent(in)    :: state_p(dim_p)
    real, intent(in)    :: stdev_p(dim_p)
    integer, intent(in) :: delt_obs
    logical, intent(in) :: write_now

! *** Local variables ***
    real :: weights


! **********************************  
! *** Add analysis step to means ***
! **********************************  

    ! Add analysis state to mean
    if (w_mm) then
       if ((step-step_null) > 0) then
          timemean = timemean + state_p / delt_obs
       endif ! step > 0
    endif ! w_mm
    if (w_sm) then
       if ((step-step_null) > 0) then
          timemean_s = timemean_s + stdev_p / delt_obs
       endif ! step > 0
    endif ! w_sm

    ! include state into monthly mean
    if ((step-step_null) > 0) then
       ! *** analyzed fields ***
       if (compute_monthly_aa) monthly_state_a  = monthly_state_a  + state_p
       if (compute_monthly_sa) monthly_stddev_a = monthly_stddev_a + stdev_p
    else if ((step-step_null) < 0) then
       ! *** forecasted fields ***
       if (compute_monthly_ff) monthly_state_f  = monthly_state_f  + state_p
       if (compute_monthly_sf) monthly_stddev_f = monthly_stddev_f + stdev_p
    end if

    if (write_now) then
       ! computing monthly mean at last day of month
       weights =  1.0/real(num_day_in_month(fleapyear,month))

       if ((step-step_null) > 0) then
          ! *** analyzed state fields ***
          if (compute_monthly_aa) monthly_state_a  = monthly_state_a  * weights
          if (compute_monthly_sa) monthly_stddev_a = monthly_stddev_a * weights
       else if ((step-step_null) < 0) then
          ! *** forecasted state fields ***
          if (compute_monthly_ff) monthly_state_f  = monthly_state_f  * weights
          if (compute_monthly_sf) monthly_stddev_f = monthly_stddev_f * weights
       end if
    endif

  end subroutine update_means_pdaf


!-----------------------------------------------------------------------
!>  Routine to reset time means 
!!
!! __Revision history:__
!! * 2025-12 - Lars Nerger - Revision for PDAF3 extracted from prepoststep_pdaf
!!
  subroutine reset_means_pdaf(step)

    use assim_pdaf_mod, &
         only: step_null

    implicit none

! *** Arguments ***
    integer, intent(in) :: step

    if ((step-step_null) > 0) then
       ! *** assimilated state fields ***
       if (compute_monthly_aa) monthly_state_a  = 0.0D0
       if (compute_monthly_sa) monthly_stddev_a = 0.0D0
    else if ((step-step_null) < 0) then
       ! *** forecasted state fields ***
       if (compute_monthly_ff) monthly_state_f  = 0.0D0
       if (compute_monthly_sf) monthly_stddev_f = 0.0D0
    end if

  end subroutine reset_means_pdaf

end module means_pdaf
