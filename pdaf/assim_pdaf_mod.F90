!>  Module for assimilation variables
!!
!! This module provides variables needed for the 
!! assimilation
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger  - Initial code
!! * 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
!! * 2022-03 - Frauke       - Removed sea ice, added ocean velocities for FESOM2.1
!! * 2025-12 - Lars Nerger  - update for PDAF3
!!
module assim_pdaf_mod

  implicit none
  save

! *** Variables specific for FESOM-PDAF ***

  character(len=100) :: nmlfile ='namelist.fesom.pdaf'    ! name of namelist file

  ! Variables for time stepping:
  integer :: step_null = 0 ! time step at initialization of PDAF
                           ! - if model (re)starts on 1st-Jan, step_null must be zero
                           ! - if model (re)starts later during year, step_null must be adapted (done in slurm-job-script)
  integer :: days_since_DAstart = 1 ! days since start of assimilation; counted onwards at model restarts
  integer :: delt_obs_ocn           ! time step interval between assimilation steps
  integer :: istep_asml             ! assimilation time step at end of an forecast phase (FESOM's "mstep" + step_null);
                                    ! FESOM's "mstep" starts counting from zero even at each restart (even mid-of year);
                                    ! istep_asml starts counting from zero on 1st-Jan
  integer :: assim_time             ! assimlation step time of day (UTC) in seconds

  ! Settings for observations:
  integer :: dim_obs          ! Number of observations
  real    :: peak_obs_error   ! Peak value used to define the observation error
  integer :: proffiles_o      ! (0) don't generate profile observation files; 
                              ! (1) generate distributed profile files
                              ! (2) generate global profile file
  integer :: start_year_o, &  ! Which years to generate profile files
       end_year_o
  integer :: use_global_obs   ! When to use global observations, or limited to sub-domain plus localiation halo
  integer :: dim_obs_max      ! Expect max. number of observations for synthetic obs.
  integer :: DA_couple_type   ! (0) for weakly-coupled, (1) for strongly-coupled assimilation
  logical :: resetforget      ! Whether to reset forgetting factor after initial phase

  ! which fields to perturb
  logical :: perturb_ssh   = .true.
  logical :: perturb_u     = .true.
  logical :: perturb_v     = .true.
  logical :: perturb_temp  = .true.
  logical :: perturb_salt  = .true.
  logical :: perturb_DIC   = .true.
  logical :: perturb_Alk   = .true.
  logical :: perturb_DIN   = .true.
  logical :: perturb_O2    = .true.

  ! Variables for adaptive localization radius
  real, allocatable :: loc_radius(:)  ! Varying localizatino radius
  integer :: loctype                  ! Type of localization
                                      ! (0) Fixed radius defined by cradius
                                      ! (1) Variable radius for constant effective observation dimension
  real :: loc_ratio                   ! Choose cradius so the effective observation dim. is loc_ratio times dim_ens

  ! Observation exclusions
  integer :: depth_excl_no
  integer, allocatable :: depth_excl(:)          ! nodes excluded in each pe

  ! File output and input - available as as namelist read-in
  logical :: read_inistate = .false.            ! Whether to read initial state from separate file
  character(len=150) :: DAoutput_path  = '.'    ! Path of DAoutput
  character(len=150) :: path_init = '.'         ! Path to initialization files
  character(len=150) :: file_init = 'covar_'    ! netcdf file holding distributed initial
                                              ! state and covariance matrix (added is _XX.nc)
  character(len=150) :: file_inistate = 'state_ini_' ! netcdf file holding distributed initial
                                                   ! state (added is _XX.nc)
  character(len=150) :: file_syntobs = 'syntobs.nc'  ! File name for synthetic observations
  character(len=150) :: path_obs_rawprof  = ''       ! Path to profile observations
  character(len=150) :: file_rawprof_prefix  = ''    ! file name prefix for profile observations 
  character(len=150) :: file_rawprof_suffix  = '.nc' ! file name suffix for profile observations 

  ! Initial ensemble covariance
  real    :: varscale=1.0           ! scaling factor for initial ensemble variance

  ! Restart information - set in slurm-job-script:
  logical :: this_is_pdaf_restart = .false.   ! init_pdaf:        - at every start, initialize PDAF-netCDF-output
                                              !                   - at restart, set forget from restart info
                                              ! init_ens_pdaf:    - at restart, skip perturbation of initial fields
                                              ! distribute_state: - at restart, skip distribution of initial fields
                                              ! add_atmos_ens_st: - at restart, read perturbed atmospheric state
                                              !                   - at restart, read forget and target value for RMSE
  logical :: start_from_ENS_spinup = .false.  ! init_ens_pdaf:    - at start from perturbed ensemble, skip perturbation of initial fields
                                              ! add_atmos_ens_st: - at start from perturbed ensemble, read perturbed atmospheric state
                                              ! distribute_state: - at start from perturbed ensemble, skip distribution of initial fields


  logical :: assimilateBGC = .false. ! whether to do a BGC assimilation step
  logical :: assimilatePHY = .false. ! whether to do a physics assimilation step

  ! Other variables - NOT available as command line options / in the namelist:

  ! Julian-Gregorian date transformation of EN4 raw data 
  integer :: num_day_in_month(0:1,12)

  ! For weakly coupled DA:
  integer :: n_sweeps                 !< Number of sweeps in local analysis loop
  character(len=3) :: type_sweep(2)   !< Type of sweep in local analysis loop
  integer :: isweep                   !< Index of sweep during the local analysis loop
  character(len=6) :: cda_phy         !< Flag whether strongly-coupled DA is done for physics data
  character(len=6) :: cda_bio         !< Flag whether strongly-coupled DA is done for bgc data

  ! Initial state in case of restarts:
  real, allocatable :: state_p_init(:)
  real, allocatable :: ens_p_init(:,:)

  ! Type variable for postprocessing:
  type obs_PP
     real, allocatable :: isExclObs (:)! whether to exclude observation due to model topography
     real, allocatable :: isInnoOmit(:)! whether to exclude observation due to Inno Omit
     real, allocatable :: nod1_g(:)    ! observation indices on global FESOM grid (nodes)
     real, allocatable :: nod2_g(:)    ! """
     real, allocatable :: nod3_g(:)    ! """
     real, allocatable :: elem_g(:)    ! observation indices on global FESOM grid (elements)
     real, allocatable :: lon(:)       ! observation coordinates
     real, allocatable :: lat(:)       ! """
     real, allocatable :: nz(:)        ! observation layer indices
     real, allocatable :: depth(:)     ! observation depth
     real, allocatable :: numrep(:)    ! number of observations on one single element
     real, allocatable :: volelem(:)   ! volume of FESOM element
  end type obs_PP



! -----------------------------------------------------------------
! *** Below are the generic variables used for configuring PDAF ***
! *** Their values are set in init_PDAF                         ***

! Settings for state vector size
  integer :: dim_state         !< Global model state dimension
  integer :: dim_state_p       !< Model state dimension for PE-local domain

! Settings for time stepping - available as command line options
  logical :: model_error       !< Control application of model error
  real    :: model_err_amp     !< Amplitude for model error

! Settings for observations - available as command line options
  integer :: delt_obs          !< time step interval between assimilation steps
  logical :: twin_experiment   !< Whether to run an twin experiment with synthetic observations
  integer :: observe_ens=0     !< (0) apply H also to ensemble mean; (1) apply H only to ensemble states
  integer :: type_obs_init=1   !< init obs. (0) before or (1) after call to prepostsstep
  logical :: do_omi_obsstats=.false. !< Whether to let OMI compute observation statistics

! General control of PDAF - available as command line options
  integer :: screen       !< Control verbosity of PDAF
                          !< * (0) no outputs
                          !< * (1) progress info
                          !< * (2) add timings
                          !< * (3) debugging output
  integer :: dim_ens      !< Size of ensemble
  integer :: filtertype   !< Select filter algorithm:
                          !<   * SEIK (1), EnKF (2), LSEIK (3), ETKF (4)
                          !<   LETKF (5), ESTKF (6), LESTKF (7), NETF (9), LNETF (10)
                          !<   LKNETF (11), PF (12), GENOBS (100), 3DVAR (200)
  integer :: subtype      !< Subtype of filter algorithm
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
  integer :: type_iau     !< Type of incremental updating:
                          !<     (0) no IAU
                          !<     (1) constant IAU weight
                          !<     (2) linear increase/decrease with maimum in middle of period
                          !<     (3) Null IAU: initialize increments arrays, but do not add increment
  integer :: steps_iau    !< Number of time steps over which IAU is applied
  integer :: dim_lag      !< Number of time instances for smoother

! Filter settings - available as command line options
!    ! General
  integer :: type_forget  !< Type of forgetting factor
                          !<  SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
                          !<   (0) fixed
                          !<   (1) global adaptive
                          !<   (2) local adaptive for LSEIK/LETKF/LESTKF
                          !<  NETF/LNETF/PF
                          !<   (0) apply inflation on forecast ensemble
                          !<   (2) apply inflation on analysis ensemble
  real    :: forget       !< Forgetting factor for filter analysis
  integer :: dim_bias     !< dimension of bias vector
!    ! All localized filters
  real    :: cradius       !< Cut-off radius for local observation domain
  integer :: locweight     !< * Type of localizing weighting of observations
                           !<   (0) constant weight of 1
                           !<   (1) exponentially decreasing with SRADIUS
                           !<   (2) use 5th-order polynomial
                           !<   (3) regulated localization of R with mean error variance
                           !<   (4) regulated localization of R with single-point error variance
  real    :: sradius       !< Support radius for 5th order polynomial
                           !<   or radius for 1/e for exponential weighting
!    ! ENKF
  integer :: rank_ana_enkf !< Rank to be considered for inversion of HPH in analysis of EnKF
                           !<  (0) for analysis w/o eigendecomposition
!    ! SEIK/ETKF/ESTKF/LSEIK/LETKF/LESTKF/NETF/LNETF/LKNETF
  integer :: type_trans    !< Type of ensemble transformation 
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
  integer :: type_sqrt     !< * Type of the transform matrix square-root 
                           !<   (0) symmetric square root
                           !<   (1) Cholesky decomposition
!    ! NETF/LNETF/PF
  integer :: type_winf     !< Set weights inflation: 
                           !<   (0) no weights inflation
                           !<   (1) use N_eff/N>limit_winf
  real    :: limit_winf    !< Limit for weights inflation: N_eff/N>limit_winf
!    ! hybrid LKNETF
  integer :: type_hyb      !< * Type of hybrid weight:
                           !<   (0) use fixed hybrid weight hyb_gamma
                           !<   (1) use gamma_lin: (1 - N_eff/N_e)*hyb_gamma
                           !<   (2) use gamma_alpha: hybrid weight from N_eff/N>=hyb_gamma
                           !<   (3) use gamma_ska: 1 - min(s,k)/sqrt(hyb_kappa) with N_eff/N>=hyb_gamma
                           !<   (4) use gamma_sklin: 1 - min(s,k)/sqrt(hyb_kappa) >= 1-N_eff/N>=hyb_gamma
  real    :: hyb_gamma     !< Hybrid filter weight for state (1.0: LETKF, 0.0 LNETF)
  real    :: hyb_kappa     !< Hybrid norm for using skewness and kurtosis
!    ! Particle filter
  integer :: pf_res_type   !< * Resampling type for PF
                           !<   (1) probabilistic resampling
                           !<   (2) stochastic universal resampling
                           !<   (3) residual resampling        
  integer :: pf_noise_type !< * Resampling type for PF
                           !<   (0) no perturbations, (1) constant stddev, 
                           !<   (2) amplitude of stddev relative of ensemble variance
  real :: pf_noise_amp     !< Noise amplitude (>=0.0, only used if pf_noise_type>0)

!    ! 3D-Var
  integer :: type_opt      !< * Type of minimizer for 3DVar
                           !<   (1) LBFGS (default)
                           !<   (2) CG+
                           !<   (3) plain CG
                           !<   (12) CG+ parallelized
                           !<   (13) plain CG parallelized
  integer :: dim_cvec = 0  !< Size of control vector (parameterized part; for subtypes 0,1)
  integer :: dim_cvec_ens = 0   !< Size of control vector (ensemble part; for subtypes 1,2)
  integer :: mcols_cvec_ens = 1 !< Multiplication factor for number of columns for ensemble control vector
  real :: beta_3dvar = 0.5 !< Hybrid weight for hybrid 3D-Var
  integer :: solver_iparam1 = 2 !< Solver specific parameter
                                !<  LBFGS: parameter m (default=5)
                                !<       Number of corrections used in limited memory matrix; 3<=m<=20
                                !<  CG+: parameter method (default=2)
                                !<       (1) Fletcher-Reeves, (2) Polak-Ribiere, (3) positive Polak-Ribiere
                                !<  CG: maximum number of iterations (default=200)
  integer :: solver_iparam2 = 1 !< Solver specific parameter
                                !<  LBFGS: - not used - 
                                !<  CG+: parameter irest (default=1)
                                !<       (0) no restarts; (n>0) restart every n steps
                                !<  CG: - not used -
  real :: solver_rparam1 = 1.0e-6 !< Solver specific parameter
                                !<  LBFGS: limit for stopping iterations 'pgtol' (default=1.0e-5)
                                !<  CG+: convergence parameter 'eps' (default=1.0e-5)
                                !<  CG: conpergence parameter 'eps' (default=1.0e-6)
  real :: solver_rparam2 = 1.0e+7 !< Solver specific parameter
                                !<  LBFGS: tolerance in termination test 'factr' (default=1.0e+7) 
                                !<  CG+: - not used -
                                !<  CG: - not used -

!    ! Other variables - _NOT_ available as command line options!
  real    :: time               !< model time

  ! Specific for local filters
  integer, allocatable :: id_lstate_in_pstate(:) !< Indices of local state vector in PE-local global state vector
  real                 :: coords_l(2)            !< Coordinates of local analysis domain

!$OMP THREADPRIVATE(coords_l, id_lstate_in_pstate)

end module assim_pdaf_mod
