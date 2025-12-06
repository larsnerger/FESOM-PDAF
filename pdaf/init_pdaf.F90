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

subroutine init_pdaf(nsteps, mesh)

  use mpi
  use PDAF                        ! PDAF interface definitions
  use statevector_pdaf, only: setup_statevector, sfields
  use mod_parallel_pdaf, &                                                ! Parallelization variables for assimilation
       only: n_modeltasks, task_id, COMM_filter, COMM_couple, filterpe, &
       mype_world, COMM_model, abort_parallel, MPIerr, &
       mype_model, mype_filter, npes_filter, writepe, mype_submodel
  use mod_assim_pdaf, &                                                   ! Variables for assimilation
       only: dim_state, dim_state_p, dim_ens, dim_lag, &
       step_null, istep_asml, assim_time, screen, filtertype, subtype, &
       delt_obs_ocn, &
       DA_couple_type, type_forget, &
       forget, locweight, cradius, sradius, &
       type_trans, type_sqrt, eff_dim_obs, loc_radius, loctype, &
       twin_experiment, dim_obs_max, use_global_obs,  &
       path_atm_cov, this_is_pdaf_restart, start_from_ENS_spinup, timemean, timemean_s, &
       resetforget, &
       debug_id_depth, debug_id_nod2, ens_member_debug, &   ! Debugging:
       proffiles_o, path_obs_rawprof, file_rawprof_prefix, file_rawprof_suffix, & ! EN4 profile data processing:
       n_sweeps, type_sweep, assimilateBGC, assimilatePHY, & ! Weak coupling of FESOM-REcoM:
       cda_phy, cda_bio, &
       state_p_init, ens_p_init  ! Initial state:
  use fesom_pdaf, &
       only: mesh_fesom, topography_p, t_mesh, &
       myDim_nod2D, MPI_COMM_FESOM, myList_edge2D, myDim_edge2D, myDim_elem2D, &
       timeold, daynew, cyearold, yearnew, yearold
  use statevector_pdaf, &
       only: id, nfields, phymin, phymax, bgcmin, bgcmax
  use mod_perturbation_pdaf, &           ! BGC parameter perturbation
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
       atmos_stochasticity_ON
  use mod_postprocess, &
       only: isPP, doPP
  use timer, only: timeit
  use mod_nc_out_routines, &
       only: netCDF_init
  use output_config_pdaf, &
       only: configure_output, setoutput
  use mod_carbon_fluxes_diags, &
       only: init_carbonfluxes_diags_out, init_carbonfluxes_diags_arrays

  implicit none

! *** Arguments ***
  integer, intent(inout) :: nsteps     ! number of model time steps
  type(t_mesh), intent(in), target :: mesh

! *** Local variables ***
  integer :: i,b,memb,n,k,s,nz ! Counters
  integer :: filter_param_i(7) ! Integer parameter array for filter
  real    :: filter_param_r(2) ! Real parameter array for filter
  integer :: status_pdaf       ! PDAF status flag
  integer :: doexit, steps     ! required arguments in call to PDAF; not defined here
  real    :: timenow           ! required arguments in call to PDAF; not defined here
  character(len=6) :: cdaval   ! Flag whether strongly-coupled DA is done
  integer, parameter :: int0=0 ! Zero
  real, allocatable :: aux(:)  ! Temporary array

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

  if (mype_submodel==0) then
     write (*,'(1x,a, i5)') 'FESOM-PDAF: INITIALIZE PDAF, task: ', task_id
  end if


! **********************************************************
! ***                  CONTROL OF PDAF                   ***
! ***              used in call to PDAF_init             ***
! **********************************************************

! *** IO options ***
  screen     = 2    ! Write screen output (1) for output, (2) add timings

! *** Filter specific variables
  filtertype = 7    ! Type of filter
                    !   (6) ESTKF
                    !   (7) LESTKF
                    !   (11) GENOBS: Generate synthetic observations
  dim_ens = n_modeltasks ! Size of ensemble for all ensemble filters
                    ! Number of EOFs to be used for SEEK
  dim_lag = 0       ! Size of lag in smoother
  subtype = 0       ! subtype of filter: 
                    !   ESTKF:
                    !     (0) Standard form of ESTKF
                    !   LESTKF:
                    !     (0) Standard form of LESTKF
  type_trans = 0    ! Type of ensemble transformation
                    !   SEIK/LSEIK and ESTKF/LESTKF:
                    !     (0) use deterministic omega
                    !     (1) use random orthonormal omega orthogonal to (1,...,1)^T
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
                    !   ETKF/LETKF:
                    !     (0) use deterministic symmetric transformation
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
  type_forget = 0   ! Type of forgetting factor in SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
                    !   (0) fixed
                    !   (1) global adaptive
                    !   (2) local adaptive for LSEIK/LETKF/LESTKF
  forget  = 1.0     ! Forgetting factor
  resetforget = .false. ! Whether to reset forgetting factor after initial phase
  type_sqrt = 0     ! Type of transform matrix square-root
                    !   (0) symmetric square root, (1) Cholesky decomposition
 

! **********************************************************
! ***     CONTROL OF USER ROUTINES FOR ASSIMILATION      ***
! **********************************************************

! *** Forecast length (time interval between analysis steps) ***
  delt_obs_ocn = 32   ! Number of time steps between analysis/assimilation steps
  assim_time   = 2700 ! Time-of-day for assimilation step in seconds
  
  step_null = 0 ! read from namelist: 0 at beginning of each year
  istep_asml = step_null

! *** Set weakly- or strongly-coupled DA of FESOM (ocean) and atmospheric component
  DA_couple_type = 0 ! (0) for weakly- (1) for strongly-coupled DA

! *** Whether to generate profile observation files
  proffiles_o = 0  ! (0) don't generate them; 
                   ! (1) generate distributed profile files
                   ! (2) generate global profile file
                   
! *** Set assimilation variables
  assim_o_sst     = .false.
  assim_o_sss     = .false.
  assim_o_sss_cci = .false.
  assim_o_ssh     = .false.

! *** specifications for observations ***
  ! This error is the standard deviation for the Gaussian distribution 
  rms_obs_sst = 0.8 ! error for satellite SST observations
  rms_obs_sss = 0.5 ! error for satellite SSS observations (SMOS)
  rms_obs_sss_cci = 0.5 ! error for satellite SSS observations (CCI)
  rms_obs_ssh = 0.05 ! error for satellite SSH observations
  rms_obs_T = 0.8    ! error for temperature profile observations
  rms_obs_S = 0.5    ! error for salinity profile observations
  bias_obs_ssh = 0.0    ! observation bias  
  bias_obs_prof = 0.0   ! observation bias  
  sst_exclude_ice = .true.  ! Exclude SST observations at point with sea ice and T>0
  sst_exclude_diff = 0.0     ! Exclude SST observations if difference from ensemble mean is >sst_exclude_diff
  prof_exclude_diff = 0.0    ! Exclude profile T observations if difference from ensemble mean is >prof_exclude_diff
  use_global_obs = 1 ! Use global full obs. (1) or full obs. limited to process domains (0)
  twin_experiment = .false.  ! Whether to run a twin experiment assimilating synthetic observations
  dim_obs_max = 80000        ! Expected maximum number of observations for synthetic obs.

! *** Localization settings
  locweight = 0     ! Type of localization weight
  cradius = 0.0     ! Cut-off radius
  sradius = cradius ! Support radius
  
! *** Configuration for atmospheric stochasticity:
  disturb_xwind=.true.
  disturb_ywind=.true.
  disturb_humi=.true.
  disturb_qlw=.true.
  disturb_qsr=.true.
  disturb_tair=.true.
  disturb_prec=.true.
  disturb_snow=.true.
  disturb_mslp=.true.
  
! *** Configuration of output frequency:
  setoutput( 1)=.false. ! daily forecast and analysis ensemble members
  setoutput( 2)=.false. ! daily forecast and analysis ensemble mean
  setoutput( 3)=.true.  ! monthly forecast and analysis ensemble mean of updated variables
  setoutput( 4)=.false. ! daily m-fields ensemble mean
  setoutput( 5)=.false. ! initial fields ensemble mean
  setoutput( 6)=.false. ! initial fields ensemble members
  setoutput( 7)=.true.  ! monthly m-fields ensemble mean
  setoutput( 8)=.false. ! daily m-fields of assimilated variables ensemble mean
  setoutput( 9)=.true.  ! daily m-fields of assimilate-able BGC variables ensemble mean
  setoutput(10)=.true.  ! daily m-fields of CO2 flux and pCO2 ensemble mean
  setoutput(11)=.true.  ! daily m-fields of variables defining the CO2 flux
  setoutput(12)=.false. ! initial fields standard deviation
  setoutput(13)=.false. ! daily forecast of standard deviation
  setoutput(14)=.true.  ! monthly forecast of standard deviation
  setoutput(15)=.true.  ! monthly forecast and analysis of standard deviation for updated fields
  setoutput(16)=.true.  ! monthly mm-fields of standard deviation
  setoutput(17)=.true.  ! if COMF-O2 assimilated: daily O2 forecast

! *** Read PDAF configuration from namelist ***
  call read_config_pdaf()


! **************************************************************************
! *** Configuration for FESOM-REcom coupling: Define local analysis loop ***
! **************************************************************************

  cda_phy = 'weak'  ! whether physics is updated from BGC assimilation
  cda_bio = 'weak'  ! whether BGC is updated from physics assimilation

  if ((assimilateBGC) .and. (assimilatePHY)) then
     ! Observations of both physics and BGC are assimilated
     n_sweeps = 2
     type_sweep(1) = 'phy'
     type_sweep(2) = 'bio'
  else
     ! Less than two observation categories
     n_sweeps = 1
     if (assimilatePHY) then
        ! Only observations of physics are assimilated
        type_sweep(1) = 'phy'
     elseif (assimilateBGC) then
        ! Only observations of BGC are assimilated
        type_sweep(1) = 'bio'
     else
        ! No observation active (free run); set sweep to physics
        type_sweep(1) = 'phy'
     end if
  end if

  if (mype_world == 0) then
     write (*,'(a,2x,a)') 'FESOM-PDAF', '*** Setup for coupled DA FESOM-REcoM ***'
     write (*, '(a,4x,a,i5)') 'FESOM-PDAF', 'Number of local analysis sweeps', n_sweeps
     write (*, '(a,4x,a)') 'FESOM-PDAF','Type of sweeps:'
     do i = 1, n_sweeps
        if (trim(type_sweep(i))=='phy') then
           cdaval = cda_phy
        else
           cdaval = cda_bio
        end if
        write (*, '(a,8x,a,i3,3x,a,a,3x,a,a)') &
             'FESOM-PDAF', 'sweep', i, ' observation type: ', trim(type_sweep(i)), 'CDA: ', trim(cdaval)
     end do
  end if


! ***********************************************************************************************
! ***   For weakly-coupled assimilation of FESOM and atmosphere re-define filter communicator ***
! ***********************************************************************************************

  if (DA_couple_type == 0) then

     ! Set filter communicator to the communicator of FESOM
     COMM_filter = MPI_COMM_FESOM

     if (filterpe) then
        call MPI_Comm_Size(COMM_filter, npes_filter, MPIerr)
        call MPI_Comm_Rank(COMM_filter, mype_filter, MPIerr)

        if (mype_filter==0) then
           write (*,'(a)') 'FESOM-Atmos-PDAF: Initialize weakly-coupled data assimilation'
        endif
     endif
  else
     if (filterpe) then
        if (mype_filter==0) then
           write (*,'(a)') 'FESOM-Atmpos-PDAF: Initialize strongly-coupled data assimilation'
        end if
     end if
  end if
  
  writepe = .false.
  if (filterpe) then
     if (mype_filter==0) writepe = .true.
  endif


! ******************************************
! *** Set pointer to FESOM mesh variable ***
! ******************************************

  mesh_fesom => mesh


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
  
! *** init m-fields and carbon diagnostic arrays

  allocate(timemean(dim_state_p))
  timemean = 0.0
  
  allocate(timemean_s(dim_state_p))
  timemean_s = 0.0
  
  call init_carbonfluxes_diags_arrays()

  
! *************************************************
! *** Perturb model parameters for ensemble run ***
! *************************************************

! Selected parameters to disturb chlorophyll and biomass:
! alfa       = 0.14   ! Initial slope of P I curve small phytoplankton
! alfa_d     = 0.19   ! Initial slope of P I curve Diatoms
! P_cm        = 3.0   ! Small phytoplankton maximum rate of phtotosynthesis
! P_cm_d      = 3.5   ! Diatom maximum rate of phtotosynthesis
! Chl2N_max   = 3.78  ! Small phytoplankton maximum Chlorophyll a to nitrogen ratio
! Chl2N_max_d = 4.2   ! Diatom maximum Chlorophyll a to nitrogen ratio
! deg_Chl     = 0.1   ! degradation rate constant
! deg_Chl_d   = 0.1   ! degradation rate constant
! graz_max, graz_max2 ! maximum grazing rates
! grazEff, grazEff2   ! grazing effeciency of ZooPlankton

! Selected parameters to disturb dissolved tracers, most importantly DIC:
! VDet, VDet_zoo2, Vdet_a   ! Sinking speed detritus
! k_din, k_din_d            ! Nitrate uptake
! res_phy, res_phy_d        ! Respiration
! res_het, res_zoo2
! biosynth
! calc_diss_rate, calc_diss_rate2 ! Dissolution during sinking
! rho_N, rho_C1             ! Remineralization of dissolved organic material
!                             (DON -> DIN; DOC --> DIC)
! lossN, lossN_d            ! Loss terms phytoplankton (PhyN --> DON)
! lossC, lossC_d
! reminN, reminC            ! Remineralization of detritus
! calc_prod_ratio           ! How much of small phytoplankton are calcifiers

  if (perturb_params_bio) then
     if (dim_ens <= 1) then
        if (mype_model==0 .and. task_id==1) write(*,*) 'FESOM-PDAF', 'Ensemble Size 1: Not perturbing BGC parameters'
     elseif (dim_ens>1) then
        if (mype_model==0 .and. task_id==1) write(*,*) 'FESOM-PDAF', 'Perturbing BGC parameters'
        call do_perturb_param_bio()
     endif
  endif

  if (perturb_params_phy) then
     if (dim_ens <= 1) then
        if (mype_model==0 .and. task_id==1) write(*,*) 'FESOM-PDAF', 'Ensemble Size 1: Not perturbing PHY parameters'
     elseif (dim_ens>1) then
        if (mype_model==0 .and. task_id==1) write(*,*) 'FESOM-PDAF', 'Perturbing PHY parameters'
        call do_perturb_param_phy()
     endif
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
          ens_p_init,   dim_state_p, MPI_DOUBLE_PRECISION, &   ! receive
          0, COMM_COUPLE, mpierr)
  endif


! *****************************************************
! *** Call PDAF initialization routine on all PEs.  ***
! ***                                               ***
! *** For all filters, first the arrays of integer  ***
! *** and real number parameters are initialized.   ***
! *** Subsequently, PDAF_init is called.            ***
! *****************************************************

  ! *** All other filters                       ***
  ! *** SEIK, LSEIK, ETKF, LETKF, ESTKF, LESTKF ***
  filter_param_i(1) = dim_state_p ! State dimension
  filter_param_i(2) = dim_ens     ! Size of ensemble
  filter_param_i(3) = 0           ! Smoother lag (not implemented here)
  filter_param_i(4) = 0           ! Not used
  filter_param_i(5) = type_forget ! Type of forgetting factor
  filter_param_i(6) = type_trans  ! Type of ensemble transformation
  filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
  filter_param_r(1) = forget      ! Forgetting factor
     
  call PDAF3_init(filtertype, subtype, step_null, &
       filter_param_i, 7,&
       filter_param_r, 1, &
       init_ens_pdaf, screen, status_pdaf)
       
  if (this_is_pdaf_restart .or. start_from_ENS_spinup) deallocate(state_p_init,ens_p_init)

! *** Check whether initialization of PDAF was successful ***
  if (status_pdaf /= 0) then
     write (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in initialization of PDAF - stopping! (PE ', mype_world,')'
     call abort_parallel()
  end if

  
! ***************************************
! *** Get domain limiting coordinates ***
! ***************************************

!~   IF (filterpe) CALL ignore_nod_pdaf() ! Seems to cause problems in FESOM2.0 (SigSegV)
                                          ! Not sure if needed

  call PDAFomi_get_domain_limits_unstr(myDim_nod2d, mesh_fesom%geo_coord_nod2D)

  
! *****************************
! *** Post processing setup ***
! *****************************

  if (isPP) then
     call doPP(nsteps)
  else


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
     call init_carbonfluxes_diags_out() ! if no file exists, file is created


! **********************************
! *** Prepare ensemble forecasts ***
! **********************************

     call MPI_BARRIER(MPI_COMM_WORLD, MPIerr)

     ! Initialize ensmeble forecasting
     call PDAF_init_forecast(next_observation_pdaf, distribute_state_pdaf, &
          prepoststep_pdaf, status_pdaf)


! ***********************************************************************************
! *** Allocate arrays for effective observation dimension and localization radius ***
! ***********************************************************************************
  
     allocate(eff_dim_obs(mydim_nod2d))
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
    
! *****************
! *** Debugging ***
! *****************

!   At this node and depth, we see a problem in the model output:
     debug_id_depth = 1
     debug_id_nod2  = 85615 !96487 !3833
     ens_member_debug = 0
   
  endif ! isPP

  call timeit(4, 'old')
  call timeit(5, 'new')

end subroutine init_pdaf
