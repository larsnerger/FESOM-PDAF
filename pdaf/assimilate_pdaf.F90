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
module assimilate_pdaf_mod
contains

  subroutine assimilate_pdaf(istep)

    use PDAF, &
         only: PDAF3_assimilate_local, PDAF_get_assim_flag
    use timer, &
         only: timeit
    use parallel_pdaf_mod, &        ! Parallelization variables
         only: mype_world, abort_parallel, task_id, mype_submodel
    use assim_pdaf_mod, &           ! Variables for assimilation
         only: istep_asml, step_null
    use fesom_pdaf, &
         only: timenew, daynew
    use means_pdaf, &
         only: do_means_pdaf

    implicit none

! *** Arguments ***
    integer, intent(in) :: istep       !< current time step

! *** Local variables ***
    integer :: status_pdaf             ! PDAF status flag
    integer :: assim_flag              ! Flag whether assimilation step was just done

! *** External subroutines ***

    ! Interface between model and PDAF, and prepoststep
    external :: collect_state_pdaf, &  ! Collect a state vector from model fields
         distribute_state_pdaf, &      ! Distribute a state vector to model fields
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


! *********************************
! *** Call assimilation routine ***
! *********************************

    call timeit(6, 'new')

    istep_asml = istep + step_null  ! istep:       starting at 1 at each model (re)start
                                    ! istep_asml:  starting at 1 at beginning of each year

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

! *** Query whether analysis step was performed
    assim_flag = 0
    call PDAF_get_assim_flag(assim_flag)

    ! Output into NEMO's ocean.output file
    if(assim_flag==1 .and. mype_submodel==0 .and. task_id==1) then
       write (*,'(a,1x,a,1x,a,1x,i5,1x,a,1x,i5,1x,a,1x,i3,1x,a,1x,i2,a,1x,i2,a)') &
         'FESOM-PDAF','DA was applied at','istep', istep, 'istep_asml', istep_asml, 'day', daynew, &
         'time', floor(timenew/3600.0),'h',int(mod(timenew,3600.0)/60.0),'min'
    endif

    call timeit(6, 'old')


! *********************************
! *** Compute daily mean        ***
! *********************************

    call timeit(7, 'new')
    call do_means_pdaf()
    call timeit(7, 'old')

  end subroutine assimilate_pdaf

end module assimilate_pdaf_mod
