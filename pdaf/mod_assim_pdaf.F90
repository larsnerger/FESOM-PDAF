!>  Module for assimilation variables
!!
!! This module provides variables needed for the 
!! assimilation within the routines of the dummy model.
!! For simplicity, all assimilation-related variables
!! are stored here, even if they are only used in
!! the main program for the filter initialization.
!! Most variables can be specified as a command line 
!! argument.
!!
!! Implementation for the 2D online example
!! with or without parallelization.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger  - Initial code
!! * 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
!! * 2022-03 - Frauke       - Removed sea ice, added ocean velocities for FESOM2.1
!! * 2025-12 - Lars Nerger  - update for PDAF3
!!
MODULE mod_assim_pdaf

!  USE statevector_pdaf
  USE fesom_pdaf

  IMPLICIT NONE
  SAVE

! *** Variables specific for FESOM-PDAF ***

  CHARACTER(len=100) :: nmlfile ='namelist.fesom.pdaf'    ! name of namelist file

! Variables for time stepping:
INTEGER :: step_null = 0          ! at initialization of PDAF, specify which time step we have, to assimilate the time-corresponding observations:
                                  !    - if model (re)starts on 1st-Jan, step_null must be zero
                                  !    - if model (re)starts later during year, step_null must be adapted (done in slurm-job-script)
INTEGER :: days_since_DAstart = 1 ! days since start of assimilation; counted onwards at model restarts
INTEGER :: delt_obs_ocn           ! time step interval between assimilation steps
INTEGER :: istep_asml             ! assimilation time step at end of an forecast phase (FESOM's "mstep" + step_null);
                                  ! FESOM's "mstep" starts counting from zero even at each restart (even mid-of year);
                                  ! istep_asml starts counting from zero on 1st-Jan
INTEGER :: assim_time             ! assimlation step time of day (UTC) in seconds

! Settings for observations:
INTEGER :: dim_obs          ! Number of observations
REAL    :: peak_obs_error   ! Peak value used to define the observation error
INTEGER :: proffiles_o      ! (0) don't generate profile observation files; 
                            ! (1) generate distributed profile files
                            ! (2) generate global profile file
INTEGER :: start_year_o, &  ! Which years to generate profile files
           end_year_o
INTEGER :: use_global_obs
INTEGER :: dim_obs_max      ! Expect max. number of observations for synthetic obs.

INTEGER :: DA_couple_type ! (0) for weakly-coupled, (1) for strongly-coupled assimilation

LOGICAL :: resetforget  ! Whether to reset forgetting factor after initial phase


! which fields to perturb
  LOGICAL :: perturb_ssh   = .TRUE.
  LOGICAL :: perturb_u     = .TRUE.
  LOGICAL :: perturb_v     = .TRUE.
  LOGICAL :: perturb_temp  = .TRUE.
  LOGICAL :: perturb_salt  = .TRUE.
  LOGICAL :: perturb_DIC   = .TRUE.
  LOGICAL :: perturb_Alk   = .TRUE.
  LOGICAL :: perturb_DIN   = .TRUE.
  LOGICAL :: perturb_O2    = .TRUE.

  CHARACTER(len=150) :: path_atm_cov



! Type variable holding field IDs in state vector
!TYPE(field_ids) :: id
!INTEGER :: nfields          ! Number of fields in state vector
! INTEGER :: phymin, phymax   ! First and last physics field in state vector
! INTEGER :: bgcmin, bgcmax   ! First and last biogeochemistry field in state vector

! Specific for local filters
INTEGER, ALLOCATABLE :: id_lstate_in_pstate(:) ! Indices of local state vector in PE-local global state vector
!REAL, ALLOCATABLE :: distance(:)               ! Distances of local observations

! Variables for adaptive localization radius
REAL, ALLOCATABLE :: eff_dim_obs(:)            ! Effective observation dimension
REAL, ALLOCATABLE :: loc_radius(:)             ! Effective observation dimension
INTEGER :: loctype       ! Type of localization
                         !   (0) Fixed radius defined by cradius
                         !   (1) Variable radius for constant effective observation dimension
REAL :: loc_ratio        ! Choose cradius so the effective observation dim. is loc_ratio times dim_ens

INTEGER, ALLOCATABLE :: id_nod2D_ice(:)        ! IDs of nodes with ice
INTEGER :: depth_excl_no
INTEGER, ALLOCATABLE :: depth_excl(:)          ! nodes excluded in each pe

! File output and input - available as as namelist read-in
LOGICAL :: read_inistate = .false.            ! Whether to read initial state from separate file
CHARACTER(len=150) :: DAoutput_path  = '.'    ! Path of DAoutput
CHARACTER(len=150) :: path_init = '.'         ! Path to initialization files
CHARACTER(len=150) :: file_init = 'covar_'    ! netcdf file holding distributed initial
                                              ! state and covariance matrix (added is _XX.nc)
CHARACTER(len=150) :: file_inistate = 'state_ini_' ! netcdf file holding distributed initial
                                                   ! state (added is _XX.nc)
CHARACTER(len=150) :: file_syntobs = 'syntobs.nc'  ! File name for synthetic observations
CHARACTER(len=150) :: path_obs_rawprof  = ''       ! Path to profile observations
CHARACTER(len=150) :: file_rawprof_prefix  = ''    ! file name prefix for profile observations 
CHARACTER(len=150) :: file_rawprof_suffix  = '.nc' ! file name suffix for profile observations 
LOGICAL :: ASIM_START_USE_CLIM_STATE = .true.

! Initial ensemble covariance
REAL    :: varscale=1.0           ! scaling factor for initial ensemble variance

! Restart information - set in slurm-job-script:
LOGICAL :: this_is_pdaf_restart = .false.            ! init_pdaf:        - at every start, initialize PDAF-netCDF-output
                                                     !                   - at restart, set forget from restart info
                                                     ! init_ens_pdaf:    - at restart, skip perturbation of initial fields
                                                     ! distribute_state: - at restart, skip distribution of initial fields
                                                     ! add_atmos_ens_st: - at restart, read perturbed atmospheric state
                                                     !                   - at restart, read forget and target value for RMSE
LOGICAL :: start_from_ENS_spinup = .false.           ! init_ens_pdaf:    - at start from perturbed ensemble, skip perturbation of initial fields
                                                     ! add_atmos_ens_st: - at start from perturbed ensemble, read perturbed atmospheric state
                                                     ! distribute_state: - at start from perturbed ensemble, skip distribution of initial fields


LOGICAL :: assimilateBGC = .false. ! whether to do a BGC assimilation step
LOGICAL :: assimilatePHY = .false. ! whether to do a physics assimilation step

! Other variables - NOT available as command line options / in the namelist:
INTEGER, ALLOCATABLE :: offset(:)          ! PE-local offsets of fields in state vector
INTEGER, ALLOCATABLE :: dim_fields(:)      ! PE-local dimensions of fields in state vector
INTEGER, ALLOCATABLE :: offset_glob(:)     ! Global offsets of fields in state vector
INTEGER, ALLOCATABLE :: dim_fields_glob(:) ! Global dimensions of fields in state vector
REAL                 :: coords_l(2)        ! Coordinates of local analysis domain
INTEGER, ALLOCATABLE :: dim_fields_l(:)    ! Field dimensions for local domain (i.e. field of vertical water column at 1 node)
INTEGER, ALLOCATABLE :: offset_l(:)        ! Field offsets for local domain

REAL, ALLOCATABLE :: state_fcst(:,:)       ! state prior to assimilation, saved to use for correction
REAL, ALLOCATABLE :: state_fcst_SSH_p(:,:) ! state prior to assimilation, saved to use for correction
REAL, ALLOCATABLE :: stdev_SSH_f_p(:)      ! forecast ensemble standard deviation at grid points for SSH field, saved to use for correction
REAL, ALLOCATABLE :: monthly_state_f(:)       ! forecasted monthly state
REAL, ALLOCATABLE :: monthly_state_a(:)       ! analyzed monthly state
REAL, ALLOCATABLE :: monthly_state_m(:)       ! monthly time-mean state
REAL, ALLOCATABLE :: monthly_state_ens_f(:,:)
REAL, ALLOCATABLE :: monthly_state_ens_a(:,:)
REAL, ALLOCATABLE :: timemean(:)     ! daily mean local state vector (mean of model forecast steps and analysis step)

REAL, ALLOCATABLE :: monthly_state_sf(:)       ! forecasted monthly standard deviation
REAL, ALLOCATABLE :: monthly_state_sa(:)       ! analyzed monthly standard deviation
REAL, ALLOCATABLE :: monthly_state_sm(:)       ! monthly time-mean standard deviation
REAL, ALLOCATABLE :: timemean_s(:)

! whether to compute monthly means:
LOGICAL :: compute_monthly_ff
LOGICAL :: compute_monthly_aa
LOGICAL :: compute_monthly_mm

LOGICAL :: compute_monthly_sf
LOGICAL :: compute_monthly_sa
LOGICAL :: compute_monthly_sm


! Julian-Gregorian date transformation of EN4 raw data 
INTEGER :: num_day_in_month(0:1,12), endday_of_month_in_year(0:1,12), startday_of_month_in_year(0:1,12)



! For weak coupling:
integer :: n_sweeps                 !< Number of sweeps in local analysis loop
character(len=3) :: type_sweep(2)   !< Type of sweep in local analysis loop
integer :: isweep                   !< Index of sweep during the local analysis loop
character(len=6) :: cda_phy   ! Flag whether strongly-coupled DA is done
character(len=6) :: cda_bio   ! Flag whether strongly-coupled DA is done

! Initial state in case of restarts:
real, allocatable :: state_p_init(:)
real, allocatable :: ens_p_init(:,:)

! For carbon diagnostics:
real, allocatable :: factor_mass(:,:)
real, allocatable :: factor_conc(:,:)

! Type variable for postprocessing:
type obs_PP
   real, allocatable :: isExclObs (:)! whether to exclude observation due to model topography
   real, allocatable :: isInnoOmit(:)! whether to exclude observation due to Inno Omit
   real, allocatable :: nod1_g(:)    ! observation indeces on global FESOM grid (nodes)
   real, allocatable :: nod2_g(:)    ! """
   real, allocatable :: nod3_g(:)    ! """
   real, allocatable :: elem_g(:)    ! observation indeces on global FESOM grid (elements)
   real, allocatable :: lon(:)       ! observation coordinates
   real, allocatable :: lat(:)       ! """
   real, allocatable :: nz(:)        ! observation layer indeces
   real, allocatable :: depth(:)     ! observation depth
   real, allocatable :: numrep(:)    ! number of observations on one single element
   real, allocatable :: volelem(:)   ! volume of FESOM element
end type obs_PP


! For debugging:
INTEGER :: debug_id_depth, & ! Location for debugging output
           debug_id_nod2           
INTEGER :: ens_member_debug
INTEGER :: mype_debug = 18
INTEGER :: node_debug = 1167
! DIC debugging: Have_obs on mype_debug=18; node_debug=1167 on 2010-01-05.
! Alk debugging: Have_obs on mype_debug=65; node_debug=915 .OR. 885 .OR. 268 on 2010-01-08.




! -----------------------------------------------------------------
! *** Below are the generic variables used for configuring PDAF ***
! *** Their values are set in init_PDAF                         ***

! Settings for state vector size
  INTEGER :: dim_state         !< Global model state dimension
  INTEGER :: dim_state_p       !< Model state dimension for PE-local domain

! Settings for time stepping - available as command line options
  LOGICAL :: model_error       !< Control application of model error
  REAL    :: model_err_amp     !< Amplitude for model error

! Settings for observations - available as command line options
  INTEGER :: delt_obs          !< time step interval between assimilation steps
  LOGICAL :: twin_experiment   !< Whether to run an twin experiment with synthetic observations
  INTEGER :: observe_ens=0     !< (0) apply H also to ensemble mean; (1) apply H only to ensemble states
  INTEGER :: type_obs_init=1   !< init obs. (0) before or (1) after call to prepostsstep
  LOGICAL :: do_omi_obsstats=.false. !< Whether to let OMI compute observation statistics

! General control of PDAF - available as command line options
  INTEGER :: screen       !< Control verbosity of PDAF
                          !< * (0) no outputs
                          !< * (1) progress info
                          !< * (2) add timings
                          !< * (3) debugging output
  INTEGER :: dim_ens      !< Size of ensemble
  INTEGER :: filtertype   !< Select filter algorithm:
                          !<   * SEIK (1), EnKF (2), LSEIK (3), ETKF (4)
                          !<   LETKF (5), ESTKF (6), LESTKF (7), NETF (9), LNETF (10)
                          !<   LKNETF (11), PF (12), GENOBS (100), 3DVAR (200)
  INTEGER :: subtype      !< Subtype of filter algorithm
                          !<   * SEIK:
                          !<     (0) ensemble forecast; new formulation
                          !<     (1) ensemble forecast; old formulation
                          !<     (2) SEIK with ensemble transformation
                          !<     (10) fixed error space basis
                          !<     (11) fixed state covariance matrix
                          !<   * LSEIK:
                          !<     (0) ensemble forecast;
                          !<     (2) LSEIK with ensemble transformation
                          !<     (10) fixed error space basis
                          !<     (11) fixed state covariance matrix
                          !<   * ETKF:
                          !<     (0) ETKF using T-matrix like SEIK
                          !<     (1) ETKF following Hunt et al. (2007)
                          !<       There are no fixed basis/covariance cases, as
                          !<       these are equivalent to SEIK subtypes 2/3
                          !<   * LETKF:
                          !<     (0) LETKF using T-matrix like SEIK
                          !<     (1) LETKF following Hunt et al. (2007)
                          !<       There are no fixed basis/covariance cases, as
                          !<       these are equivalent to LSEIK subtypes 2/3
                          !<   * EnKF:
                          !<     (0) analysis for large observation dimension
                          !<     (1) analysis for small observation dimension
                          !<   * LEnKF:
                          !<     (0) standard analysis
                          !<   * ESTKF:
                          !<     (0) standard ESTKF 
                          !<       There are no fixed basis/covariance cases, as
                          !<       these are equivalent to SEIK subtypes 2/3
                          !<   * LESTKF:
                          !<     (0) standard LESTKF 
                          !<       There are no fixed basis/covariance cases, as
                          !<       these are equivalent to LSEIK subtypes 2/3
                          !<   * NETF:
                          !<     (0) standard NETF 
                          !<   * LNETF:
                          !<     (0) standard LNETF 
                          !<   * LKNETF:
                          !<     (0) HNK: 2-step LKNETF with NETF before LETKF
                          !<     (1) HKN: 2-step LKNETF with LETKF before NETF
                          !<     (2) HSync: LKNETF synchronous
                          !<   * PF:
                          !<     (0) standard PF 
                          !<   * ENSRF/EAKF:
                          !<     (0) ENSRF with serial observation processing
                          !<     (1) EAKF with loca least square regression
                          !<   * 3D-Var:
                          !<     (0) parameterized 3D-Var
                          !<     (1) 3D Ensemble Var using LESTKF for ensemble update
                          !<     (2) 3D Ensemble Var using ESTKF for ensemble update
                          !<     (3) hybrid 3D-Var using LESTKF for ensemble update
                          !<     (4) hybrid 3D-Var using ESTKF for ensemble update
  INTEGER :: type_iau     !< Type of incremental updating:
                          !<     (0) no IAU
                          !<     (1) constant IAU weight
                          !<     (2) linear increase/decrease with maimum in middle of period
                          !<     (3) Null IAU: initialize increments arrays, but do not add increment
  INTEGER :: steps_iau    !< Number of time steps over which IAU is applied
  INTEGER :: dim_lag      !< Number of time instances for smoother

! Filter settings - available as command line options
!    ! General
  INTEGER :: type_forget  !< Type of forgetting factor
                          !<  SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
                          !<   (0) fixed
                          !<   (1) global adaptive
                          !<   (2) local adaptive for LSEIK/LETKF/LESTKF
                          !<  NETF/LNETF/PF
                          !<   (0) apply inflation on forecast ensemble
                          !<   (2) apply inflation on analysis ensemble
  REAL    :: forget       !< Forgetting factor for filter analysis
  INTEGER :: dim_bias     !< dimension of bias vector
!    ! All localized filters
  REAL    :: cradius       !< Cut-off radius for local observation domain
  INTEGER :: locweight     !< * Type of localizing weighting of observations
                           !<   (0) constant weight of 1
                           !<   (1) exponentially decreasing with SRADIUS
                           !<   (2) use 5th-order polynomial
                           !<   (3) regulated localization of R with mean error variance
                           !<   (4) regulated localization of R with single-point error variance
  REAL    :: sradius       !< Support radius for 5th order polynomial
                           !<   or radius for 1/e for exponential weighting
!    ! ENKF
  INTEGER :: rank_ana_enkf !< Rank to be considered for inversion of HPH in analysis of EnKF
                           !<  (0) for analysis w/o eigendecomposition
!    ! SEIK/ETKF/ESTKF/LSEIK/LETKF/LESTKF/NETF/LNETF/LKNETF
  INTEGER :: type_trans    !< Type of ensemble transformation 
                           !< * SEIK/LSEIK: 
                           !< (0) use deterministic omega
                           !< (1) use random orthonormal omega orthogonal to (1,...,1)^T
                           !< (2) use product of (0) with random orthonormal matrix with
                           !<     eigenvector (1,...,1)^T 
                           !< * ETKF/LETKF with subtype=4: 
                           !< (0) use deterministic symmetric transformation
                           !< (2) use product of (0) with random orthonormal matrix with
                           !<     eigenvector (1,...,1)^T 
                           !< * ESTKF/LESTKF:
                           !< (0) use deterministic omega
                           !< (1) use random orthonormal omega orthogonal to (1,...,1)^T
                           !< (2) use product of (0) with random orthonormal matrix with
                           !<     eigenvector (1,...,1)^T
                           !< * NETF/LNETF:
                           !< (0) use random orthonormal transformation orthogonal to (1,...,1)^T
                           !< (1) use identity transformation
                           !< * LKNETF:
                           !< (0) use random orthonormal transformation orthogonal to (1,...,1)^T
                           !< (1) use identity transformation
!    ! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
  INTEGER :: type_sqrt     !< * Type of the transform matrix square-root 
                           !<   (0) symmetric square root
                           !<   (1) Cholesky decomposition
!    ! NETF/LNETF/PF
  INTEGER :: type_winf     !< Set weights inflation: 
                           !<   (0) no weights inflation
                           !<   (1) use N_eff/N>limit_winf
  REAL    :: limit_winf    !< Limit for weights inflation: N_eff/N>limit_winf
!    ! hybrid LKNETF
  INTEGER :: type_hyb      !< * Type of hybrid weight:
                           !<   (0) use fixed hybrid weight hyb_gamma
                           !<   (1) use gamma_lin: (1 - N_eff/N_e)*hyb_gamma
                           !<   (2) use gamma_alpha: hybrid weight from N_eff/N>=hyb_gamma
                           !<   (3) use gamma_ska: 1 - min(s,k)/sqrt(hyb_kappa) with N_eff/N>=hyb_gamma
                           !<   (4) use gamma_sklin: 1 - min(s,k)/sqrt(hyb_kappa) >= 1-N_eff/N>=hyb_gamma
  REAL    :: hyb_gamma     !< Hybrid filter weight for state (1.0: LETKF, 0.0 LNETF)
  REAL    :: hyb_kappa     !< Hybrid norm for using skewness and kurtosis
!    ! Particle filter
  INTEGER :: pf_res_type   !< * Resampling type for PF
                           !<   (1) probabilistic resampling
                           !<   (2) stochastic universal resampling
                           !<   (3) residual resampling        
  INTEGER :: pf_noise_type !< * Resampling type for PF
                           !<   (0) no perturbations, (1) constant stddev, 
                           !<   (2) amplitude of stddev relative of ensemble variance
  REAL :: pf_noise_amp     !< Noise amplitude (>=0.0, only used if pf_noise_type>0)

!    ! 3D-Var
  INTEGER :: type_opt      !< * Type of minimizer for 3DVar
                           !<   (1) LBFGS (default)
                           !<   (2) CG+
                           !<   (3) plain CG
                           !<   (12) CG+ parallelized
                           !<   (13) plain CG parallelized
  INTEGER :: dim_cvec = 0  !< Size of control vector (parameterized part; for subtypes 0,1)
  INTEGER :: dim_cvec_ens = 0   !< Size of control vector (ensemble part; for subtypes 1,2)
  INTEGER :: mcols_cvec_ens = 1 !< Multiplication factor for number of columns for ensemble control vector
  REAL :: beta_3dvar = 0.5 !< Hybrid weight for hybrid 3D-Var
  INTEGER :: solver_iparam1 = 2 !< Solver specific parameter
                                !<  LBFGS: parameter m (default=5)
                                !<       Number of corrections used in limited memory matrix; 3<=m<=20
                                !<  CG+: parameter method (default=2)
                                !<       (1) Fletcher-Reeves, (2) Polak-Ribiere, (3) positive Polak-Ribiere
                                !<  CG: maximum number of iterations (default=200)
  INTEGER :: solver_iparam2 = 1 !< Solver specific parameter
                                !<  LBFGS: - not used - 
                                !<  CG+: parameter irest (default=1)
                                !<       (0) no restarts; (n>0) restart every n steps
                                !<  CG: - not used -
  REAL :: solver_rparam1 = 1.0e-6 !< Solver specific parameter
                                !<  LBFGS: limit for stopping iterations 'pgtol' (default=1.0e-5)
                                !<  CG+: convergence parameter 'eps' (default=1.0e-5)
                                !<  CG: conpergence parameter 'eps' (default=1.0e-6)
  REAL :: solver_rparam2 = 1.0e+7 !< Solver specific parameter
                                !<  LBFGS: tolerance in termination test 'factr' (default=1.0e+7) 
                                !<  CG+: - not used -
                                !<  CG: - not used -

!    ! Other variables - _NOT_ available as command line options!
  REAL    :: time               !< model time
  REAL, ALLOCATABLE :: coords_p(:,:)    !< Coordinates of process-local state vector entries
                                        !< needed to intiialize localization for LEnKF/ENSRF

END MODULE mod_assim_pdaf
