!> Initialize PDAF
!!
!! This routine collects the initialization of variables for PDAF.
!! In addition, the initialization routine PDAF_init is called
!! such that the internal initialization of PDAF is performed.
!!
!! __Revision history:__
!! * 2008-10 - Lars Nerger  - Initial code
!! * 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
!! * 2022-02 - Frauke B     - Adapted for FESOM2.1
!! * 2025-12 - Lars Nerger  - Update for PDAF3
!!
module init_pdaf_mod
contains

  subroutine init_pdaf(nsteps, mesh)

    use mpi
    use PDAF, only: &                                                       ! PDAF functions
         PDAF3_init, PDAF_set_iparam, PDAF_init_forecast, &
         PDAFomi_get_domain_limits_unstr, PDAF_reset_forget
    use timer, only: timeit
    use statevector_pdaf, only: setup_statevector, sfields
    use parallel_pdaf_mod, &                                                ! Parallelization variables for assimilation
         only: n_modeltasks, task_id, COMM_filter, COMM_couple, filterpe, &
         mype_world, COMM_model, abort_parallel, MPIerr, &
         mype_model, mype_filter, npes_filter, writepe, mype_submodel
    use assim_pdaf_mod, &                                                   ! Variables for assimilation
         only: dim_state, dim_state_p, dim_ens, dim_lag, &
         step_null, istep_asml, assim_time, screen, filtertype, subtype, &
         delt_obs_ocn, &
         type_forget, forget, locweight, cradius, sradius, &
         type_trans, type_sqrt, loc_ratio, loc_radius, &
         twin_experiment, dim_obs_max, use_global_obs,  &
         this_is_pdaf_restart, start_from_ENS_spinup, &
         resetforget, &
         proffiles_o, path_obs_rawprof, file_rawprof_prefix, file_rawprof_suffix, & ! EN4 profile data processing:
         state_p_init, ens_p_init                                           ! Initial state
    use coupled_da_mod, &
         only: cda_phy, cda_bio, cda_set_sweeps, DA_couple_type, cda_reset_filter_comm
    use fesom_pdaf, &
         only: mesh_fesom, topography_p, t_mesh, &
         myDim_nod2D, MPI_COMM_FESOM, myList_edge2D, myDim_edge2D, myDim_elem2D, &
         timeold, daynew, cyearold, yearnew, yearold
    use statevector_pdaf, &
         only: id, nfields, phymin, phymax, bgcmin, bgcmax
    use adaptive_lradius_pdaf, &
         only: init_adaptive_lradius_pdaf
    use mod_perturbation_pdaf, &                                            ! BGC parameter perturbation
         only: perturb_scale, perturb_params_bio, perturb_params_phy, &
         perturb_lognormal, perturb_scaleD, &
         do_perturb_param_bio, do_perturb_param_phy
    use obs_sss_smos_pdafomi, &
         only: assim_o_sss, rms_obs_sss, path_obs_sss, file_sss_prefix, file_sss_suffix, &
         sss_exclude_ice, sss_exclude_diff, bias_obs_sss, sss_fixed_rmse
    use obs_sss_cci_pdafomi, &
         only: assim_o_sss_cci, rms_obs_sss_cci, path_obs_sss_cci, file_sss_cci_prefix, file_sss_cci_suffix, &
         sss_cci_exclude_ice, sss_cci_exclude_diff, bias_obs_sss_cci, sss_cci_fixed_rmse
    use obs_ssh_cmems_pdafomi, &
         only: assim_o_ssh, rms_obs_ssh, path_obs_ssh, file_ssh_prefix, file_ssh_suffix, &
         ssh_exclude_ice, ssh_exclude_diff, bias_obs_ssh, ssh_fixed_rmse
    use obs_sst_pdafomi, &
         only: assim_o_sst, rms_obs_sst, path_obs_sst, file_sst_prefix, file_sst_suffix, &
         sst_exclude_ice, sst_exclude_diff, bias_obs_sst, sst_fixed_rmse
    use obs_TSprof_EN4_pdafomi, &
         only: assim_o_en4_s, assim_o_en4_t, &
         rms_obs_S, rms_obs_T, &
         path_obs_prof, file_prof_prefix, file_prof_suffix, &
         bias_obs_prof, prof_exclude_diff
    use mod_atmos_ens_stochasticity, &
         only: disturb_xwind, disturb_ywind, disturb_humi, &
         disturb_qlw, disturb_qsr, disturb_tair, &
         disturb_prec, disturb_snow, disturb_mslp, &
         init_atmos_ens_stochasticity, init_atmos_stochasticity_output,&
         atmos_stochasticity_ON, path_atm_cov
    use means_pdaf, &
         only: init_means_pdaf
    use mod_postprocess, &
         only: isPP, doPP
    use mod_nc_out_routines, &
         only: netCDF_init
    use output_config_pdaf, &
         only: configure_output, setoutput
    use cfluxes_diags_pdaf, &
         only: init_cfluxes_diags_out, init_cfluxes_diags_arrays

    implicit none

! *** Arguments ***
    integer, intent(inout) :: nsteps     ! number of model time steps
    type(t_mesh), intent(in), target :: mesh

! *** Local variables ***
    integer :: i                 ! Counter
    integer :: filter_param_i(2) ! Integer parameter array for filter
    real    :: filter_param_r(1) ! Real parameter array for filter
    integer :: status_pdaf       ! PDAF status flag
    integer :: doexit, steps     ! required arguments in call to PDAF; not defined here
    real    :: timenow           ! required arguments in call to PDAF; not defined here
    character(len=6) :: cdaval   ! Flag whether strongly-coupled DA is done

    ! External subroutines
    external :: init_ens_pdaf            ! Ensemble initialization
    external :: next_observation_pdaf, & ! Provide time step and model time  of next observation
         distribute_state_pdaf, &        ! Routine to distribute a state vector to model fields
         prepoststep_pdaf                ! User supplied pre/poststep routine


! ***************************
! ***   Initialize PDAF   ***
! ***************************

    call timeit(3, 'old')
    call timeit(4, 'new')

    ! Get process-ID in task of model compartment
    call MPI_Comm_Rank(MPI_COMM_FESOM, mype_submodel, MPIerr)

    ! In case of coupled Ocean-atmosphere DA reset communicator
    call cda_reset_filter_comm(MPI_COMM_FESOM)

    ! Set pointer to FESOM mesh variable
    mesh_fesom => mesh

    writepe = .false.
    if (filterpe) then
       if (mype_filter==0) writepe = .true.
    endif

    if (mype_submodel==0) then
       write (*,'(1x,a, i5)') 'FESOM-PDAF: INITIALIZE PDAF, task: ', task_id
    end if


! **********************************************************
! ***                  CONTROL OF PDAF                   ***
! ***              used in call to PDAF_init             ***
! **********************************************************

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! +++ Here some default values are specified mainly to   +++
    ! +++ show which configuration variable can be modified  +++
    ! +++                                                    +++ 
    ! +++ For an experiment, all configuration values should +++
    ! +++ be changed in the namelist file that can be        +++
    ! +++ archived with an assimialtion run.                 +++
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! *** IO options ***
    screen = 2             ! Write screen output (1) for output, (2) add timings

! *** Filter specific variables
    filtertype = 7         ! Type of filter
    dim_ens = n_modeltasks ! Size of ensemble
    dim_lag = 0            ! Size of lag in smoother
    subtype = 0            ! subtype of filter: 
    type_trans = 0         ! Type of ensemble transformation
    type_forget = 0        ! Type of forgetting factor
    forget  = 1.0          ! Forgetting factor
    resetforget = .false.  ! Whether to reset forgetting factor after initial phase
    type_sqrt = 0          ! Type of transform matrix square-root
 

! **********************************************************
! ***     CONTROL OF USER ROUTINES FOR ASSIMILATION      ***
! **********************************************************

! *** Forecast length (time interval between analysis steps) ***
    delt_obs_ocn = 32      ! Number of time steps between analysis/assimilation steps
    assim_time   = 2700    ! Time-of-day for assimilation step in seconds
  
    step_null = 0          ! read from namelist: 0 at beginning of each year
    istep_asml = step_null ! Assimilation step counter

! *** Whether to generate profile observation files
    proffiles_o = 0        ! (0) don't generate them; 
                           ! (1) generate distributed profile files
                           ! (2) generate global profile file

! *** Set assimilation variables
    assim_o_sst     = .false.
    assim_o_sss     = .false.
    assim_o_sss_cci = .false.
    assim_o_ssh     = .false.

! *** specifications for observations ***
    ! This error is the standard deviation for the Gaussian distribution 
    rms_obs_sst = 0.8         ! error for satellite SST observations
    rms_obs_sss = 0.5         ! error for satellite SSS observations (SMOS)
    rms_obs_sss_cci = 0.5     ! error for satellite SSS observations (CCI)
    rms_obs_ssh = 0.05        ! error for satellite SSH observations
    rms_obs_T = 0.8           ! error for temperature profile observations
    rms_obs_S = 0.5           ! error for salinity profile observations

    bias_obs_ssh = 0.0        ! observation bias  
    bias_obs_prof = 0.0       ! observation bias  
    sst_exclude_ice = .true.  ! Exclude SST observations at point with sea ice and T>0
    sst_exclude_diff = 0.0    ! Exclude SST observations if difference from ensemble mean is >sst_exclude_diff
    prof_exclude_diff = 0.0   ! Exclude profile T observations if difference from ensemble mean is >prof_exclude_diff
    use_global_obs = 1        ! Use global full obs. (1) or full obs. limited to process domains (0)

    twin_experiment = .false. ! Whether to run a twin experiment assimilating synthetic observations
    dim_obs_max = 80000       ! Expected maximum number of observations for synthetic obs.

! *** Coupled physics-BGC DA
    cda_phy = 'weak'          ! whether physics is updated from BGC assimilation
    cda_bio = 'weak'          ! whether BGC is updated from physics assimilation

! *** Localization settings
    locweight = 4             ! Type of localization weight
    cradius = 2.0e5           ! Cut-off radius
    sradius = cradius         ! Support radius

! *** Configuration of output frequency:
    setoutput( 1)=.false.     ! daily forecast and analysis ensemble members
    setoutput( 2)=.false.     ! daily forecast and analysis ensemble mean
    setoutput( 3)=.false.     ! monthly forecast and analysis ensemble mean of updated variables
    setoutput( 4)=.false.     ! daily m-fields ensemble mean
    setoutput( 5)=.false.     ! initial fields ensemble mean
    setoutput( 6)=.false.     ! initial fields ensemble members
    setoutput( 7)=.false.     ! monthly m-fields ensemble mean
    setoutput( 8)=.false.     ! daily m-fields of assimilated variables ensemble mean
    setoutput( 9)=.false.     ! daily m-fields of assimilate-able BGC variables ensemble mean
    setoutput(10)=.false.     ! daily m-fields of CO2 flux and pCO2 ensemble mean
    setoutput(11)=.false.     ! daily m-fields of variables defining the CO2 flux
    setoutput(12)=.false.     ! initial fields standard deviation
    setoutput(13)=.false.     ! daily forecast of standard deviation
    setoutput(14)=.false.     ! monthly forecast of standard deviation
    setoutput(15)=.false.     ! monthly forecast and analysis of standard deviation for updated fields
    setoutput(16)=.false.     ! monthly mm-fields of standard deviation
    setoutput(17)=.false.     ! if COMF-O2 assimilated: daily O2 forecast

! *** Read PDAF configuration from namelist ***
    call read_config_pdaf()


! *********************************************************
! *** For coupled DA in FESOM-REcoM:                    ***
! *** Define number and types of local analysis sweeps  ***
! *********************************************************

    call cda_set_sweeps()


! ***************************
! *** Define state vector ***
! ***************************

    call setup_statevector(dim_state, dim_state_p, screen)

    ! Set land mask and compute volumes
    call init_topography(dim_state_p, dim_state)

    ! Set configuration fro file output
    call configure_output()
  
! *** Initial Screen output PDAF ***

    if (mype_model==0 .and. task_id==1) call init_pdaf_info()

! *** Check ensemble size
    if (dim_ens /= n_modeltasks) then
       write (*,*) 'ERROR: Ensemble size (',dim_ens, &
            ') needs to be identical to the number of model tasks (',n_modeltasks,')'
       call abort_parallel()
    end if
  
! ***************************************
! *** Initializations for diagnostics ***
! ***************************************

    ! Initialize time-mean fields (over all time steps)
    call init_means_pdaf(dim_state_p)
  
    ! Initialize carbon flux diagnostics
    call init_cfluxes_diags_arrays()

  
! *************************************************
! *** Perturb model parameters for ensemble run ***
! *************************************************

    if (perturb_params_bio .and. dim_ens>1) then
       if (mype_model==0 .and. task_id==1) write(*,*) 'FESOM-PDAF', 'Perturbing BGC parameters'
       call do_perturb_param_bio()
    endif

    if (perturb_params_phy .and. dim_ens>1) then
       if (mype_model==0 .and. task_id==1) write(*,*) 'FESOM-PDAF', 'Perturbing PHY parameters'
       call do_perturb_param_phy()
    endif


! *********************************
! *** Compute velocaty at nodes ***
! *********************************

    call compute_vel_nodes(mesh_fesom)


! ******************************************************
! *** Initialize state from model in case of restart ***
! ******************************************************

    if (this_is_pdaf_restart .or. start_from_ENS_spinup) then

       allocate(state_p_init(dim_state_p))
       allocate(ens_p_init(dim_state_p,dim_ens))
   
       ! collect model initial fields
       ! note: diagnostic fields, e.g. pCO2, not available: computed at the end model time step
       call collect_state_pdaf(dim_state_p,state_p_init)
       state_p_init = topography_p * state_p_init
   
       call MPI_GATHER(state_p_init, dim_state_p, MPI_DOUBLE_PRECISION, &   ! send
            ens_p_init, dim_state_p, MPI_DOUBLE_PRECISION, &                ! receive
            0, COMM_COUPLE, mpierr)
    endif


! *****************************************************
! *** Call PDAF initialization routine on all PEs.  ***
! *****************************************************

    filter_param_i(1) = dim_state_p ! State dimension
    filter_param_i(2) = dim_ens     ! Size of ensemble
    filter_param_r(1) = forget      ! Forgetting factor

    call PDAF3_init(filtertype, subtype, step_null, &
         filter_param_i, 2,&
         filter_param_r, 1, &
         init_ens_pdaf, screen, status_pdaf)

    call PDAF_set_iparam(5, type_forget, status_pdaf) ! Type of forgetting factor
    call PDAF_set_iparam(6, type_trans, status_pdaf)  ! Type of ensemble transformation
    call PDAF_set_iparam(7, type_sqrt, status_pdaf)   ! Type of transform square-root (SEIK-sub2/ESTKF)

    ! *** Check whether initialization of PDAF was successful ***
    if (status_pdaf /= 0) then
       write (*,'(/1x,a6,i3,a43,i4,a1/)') &
            'ERROR ', status_pdaf, &
            ' in initialization of PDAF - stopping! (PE ', mype_world,')'
       call abort_parallel()
    end if

    if (this_is_pdaf_restart .or. start_from_ENS_spinup) then
       deallocate(state_p_init,ens_p_init)
    end if


! ***************************************
! *** Get domain limiting coordinates ***
! ***************************************

!~   IF (filterpe) CALL ignore_nod_pdaf() ! Seems to cause problems in FESOM2.0 (SigSegV)
                                          ! Not sure if needed

    call PDAFomi_get_domain_limits_unstr(myDim_nod2d, mesh_fesom%geo_coord_nod2D)

  

    if (.not. isPP) then
       
    ! *** Assimilation run ***

! ******************************
! *** Initialize file output ***
! ******************************

       ! Initialize PDAF netCDF output file:
       !    - at beginning of each new year
       !    - at start of assimililation experiment
       if ((.not. this_is_pdaf_restart) .or. (daynew==1)) then
          call netCDF_init()
       endif
     
       ! carbon flux diagnostics
       call init_cfluxes_diags_out()


! **********************************
! *** Prepare ensemble forecasts ***
! **********************************

       call MPI_BARRIER(MPI_COMM_WORLD, MPIerr)

       ! Initialize ensemble forecasting
       call PDAF_init_forecast(next_observation_pdaf, distribute_state_pdaf, &
            prepoststep_pdaf, status_pdaf)


! ***********************************************************************************
! *** Allocate arrays for effective observation dimension and localization radius ***
! ***********************************************************************************
  
       call init_adaptive_lradius_pdaf(mydim_nod2d, loc_ratio)

       allocate(loc_radius(mydim_nod2d))


! **************************************************
! *** Initialize file for synthetic observations ***
! **************************************************

       if (filtertype==100 .and. mype_world==0) then
          write(*,*) "Synthetic observations not yet implemented in this version - stopping!"
          call abort_parallel
       end if


! ***********************************
! **** Atmospheric stochasticity  ***
! ***********************************
  
       if (atmos_stochasticity_ON) then
  
          ! initialize atmospheric stochasticity at (re)start
          call init_atmos_ens_stochasticity()

          ! create stochasticity file
          !    - at beginning of every new year
          !    - at first start of assimilation experiment
          if ((.not. (this_is_pdaf_restart .or. start_from_ENS_spinup)) .or. (yearnew .ne. yearold)) then
             call init_atmos_stochasticity_output()
          endif
       endif
  
       ! resetforget scheme: in case of restart, reset forget to 0.99 or 1.00
       ! note: at restarts, forgetting factor is saved and read with atmospheric stochasticity
       if (resetforget .and. (this_is_pdaf_restart .or. start_from_ENS_spinup)) then
          call PDAF_reset_forget(forget)
       endif

    else

       ! *** Post processing setup ***

       call doPP(nsteps)

    endif ! isPP

    call timeit(4, 'old')
    call timeit(5, 'new')

  end subroutine init_pdaf

end module init_pdaf_mod
