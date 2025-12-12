!> PDAF-OMI template observation module 
!!
!! This module handles operations for one data type (called 'module-type' below).
!!
!! The subroutines in this module are for the particular handling of
!! a single observation type.
!! The routines are called by the different call-back routines of PDAF.
!! Most of the routines are generic so that in practice only 2 routines
!! need to be adapted for a particular data type. These are the routines
!! for the initialization of the observation information (init_dim_obs)
!! and for the observation operator (obs_op).
!!
!! The module and the routines are named according to the observation type.
!! This allows to distinguish the observation type and the routines in this
!! module from other observation types.
!!
!! The module uses two derived data type (obs_f and obs_l), which contain
!! all information about the full and local observations. Only variables
!! of the type obs_f need to be initialized in this module. The variables
!! in the type obs_l are initilized by the generic routines from PDAFomi.
!!
!!
!! These 2 routines need to be adapted for the particular observation type:
!! * init_dim_obs_TYPE \n
!!           Count number of process-local and full observations; 
!!           initialize vector of observations and their inverse variances;
!!           initialize coordinate array and index array for indices of
!!           observed elements of the state vector.
!! * obs_op_TYPE \n
!!           observation operator to get full observation vector of this type. Here
!!           one has to choose a proper observation operator or implement one.
!!
!! In addition, there are two optional routine, which are required if filters 
!! with localization are used:
!! * init_dim_obs_l_TYPE \n
!!           Only required if domain-localized filters (e.g. LESTKF, LETKF) are used:
!!           Count number of local observations of module-type according to
!!           their coordinates (distance from local analysis domain). Initialize
!!           module-internal distances and index arrays.
!! * localize_covar_TYPE \n
!!           Only required if the localized EnKF is used:
!!           Apply covariance localization in the LEnKF.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
module obs_n_merged_pdafomi

  use parallel_pdaf_mod, &
       only: mype_filter, writepe
  use PDAF, &                     ! Declaration of observation data types
       only: obs_f, obs_l
  use coupled_da_mod, &           ! Variables for coupled DA
       only: n_sweeps 
  use assim_pdaf_mod, &
       only: obs_PP
 
  implicit none
  save

  ! Variables which are inputs to the module (usually set in init_pdaf)
  logical :: assim_o_n_merged   !< Whether to assimilate

  ! Further variables specific for the profile observations
  character(len=110) :: path_obs_n_merged  = ''      !< Path to profile observations
  character(len=110) :: file_n_merged

  real    :: rms_obs_n_merged      ! Observation error
  real    :: bias_obs_n_merged     ! Observation bias

  real    :: lradius_n_merged      ! Localization radius
  real    :: sradius_n_merged      ! Support radius for localization function

  real    :: n_merged_exclude_diff ! limit difference beyond which observations are excluded
                                    ! using PDAF's inno_omit functionality
                                    ! set 0.0 to deactivate
  real    :: n_merged_excl_relative   ! relative difference to exclude observations based on forecast
  real    :: n_merged_excl_absolute   ! absolute difference to exclude observations based on forecast

  real, allocatable :: mean_n_p(:)              ! ensemble mean for observation exclusion
  real, allocatable :: loc_radius_n_merged(:)   ! localization radius array
  
  real, allocatable :: ivariance_obs_g(:)      ! global-earth inverse observation variances
  
  real, parameter   :: refdens  = 1.026        ! reference density of water for unit conversion
  real, parameter   :: irefdens = 0.975        ! inverse """
  real, parameter   :: third = 0.3333333333333333

  integer, parameter :: val1 =1
  integer, parameter :: val2 =2
  integer, parameter :: val3 =3
  integer, parameter :: dens1=4
  integer, parameter :: dens2=5
  integer, parameter :: dens3=6

! ***********************************************************************
! *** The following two data types are used in PDAFomi                ***
! *** They are declared in PDAFomi and only listed here for reference ***
! ***********************************************************************

! Data type to define the full observations by internally shared variables of the module
!   TYPE obs_f
!           Mandatory variables to be set in INIT_DIM_OBS
!      INTEGER :: doassim                   ! Whether to assimilate this observation type
!      INTEGER :: disttype                  ! Type of distance computation to use for localization
!                                           ! (0) Cartesian, (1) Cartesian periodic
!                                           ! (2) simplified geographic, (3) geographic haversine function
!      INTEGER :: ncoord                    ! Number of coordinates use for distance computation
!      INTEGER, ALLOCATABLE :: id_obs_p(:,:) ! Indices of observed field in state vector (process-local)
!           
!           Optional variables - they can be set in INIT_DIM_OBS
!      REAL, ALLOCATABLE :: icoeff_p(:,:)   ! Interpolation coefficients for obs. operator
!      REAL, ALLOCATABLE :: domainsize(:)   ! Size of domain for periodicity (<=0 for no periodicity)
!
!           Variables with predefined values - they can be changed in INIT_DIM_OBS
!      INTEGER :: obs_err_type=0            ! Type of observation error: (0) Gauss, (1) Laplace
!      INTEGER :: use_global_obs=1          ! Whether to use (1) global full obs. 
!                                           ! or (0) obs. restricted to those relevant for a process domain
!
!           The following variables are set in the routine PDAFomi_gather_obs
!      INTEGER :: dim_obs_p                 ! number of PE-local observations
!      INTEGER :: dim_obs_f                 ! number of full observations
!      INTEGER :: dim_obs_g                 ! global number of observations
!      INTEGER :: off_obs_f                 ! Offset of this observation in overall full obs. vector
!      INTEGER :: off_obs_g                 ! Offset of this observation in overall global obs. vector
!      INTEGER :: obsid                     ! Index of observation over all assimilated observations
!      REAL, ALLOCATABLE :: obs_f(:)        ! Full observed field
!      REAL, ALLOCATABLE :: ocoord_f(:,:)   ! Coordinates of full observation vector
!      REAL, ALLOCATABLE :: ivar_obs_f(:)   ! Inverse variance of full observations
!      INTEGER, ALLOCATABLE :: id_obs_f_lim(:) ! Indices of domain-relevant full obs. in global vector of obs.
!                                           ! (only if full obs. are restricted to process domain))
!   END TYPE obs_f

! Data type to define the local observations by internally shared variables of the module
!   TYPE obs_l
!      INTEGER :: dim_obs_l                 ! number of local observations
!      INTEGER :: off_obs_l                 ! Offset of this observation in overall local obs. vector
!      INTEGER, ALLOCATABLE :: id_obs_l(:)  ! Indices of local observations in full obs. vector 
!      REAL, ALLOCATABLE :: distance_l(:)   ! Distances of local observations
!      REAL, ALLOCATABLE :: ivar_obs_l(:)   ! Inverse variance of local observations
!      INTEGER :: locweight                 ! Specify localization function
!      REAL :: lradius                      ! localization radius
!      REAL :: sradius                      ! support radius for localization function
!   END TYPE obs_l
! ***********************************************************************

! Declare instances of observation data types used here
! We use generic names here, but one could renamed the variables
  type(obs_f), target, public :: thisobs      ! full observation
  type(obs_l), target, public :: thisobs_l    ! local observation
  
  type(obs_PP) :: thisobs_PP
  type(obs_PP) :: thisobs_PP_f
  
  logical :: isPP = .false.                 ! T: for postprocessing of simulation output
                                            ! F: for assimilation

!$OMP THREADPRIVATE(thisobs_l)


!-------------------------------------------------------------------------------

contains

!> Initialize information on the module-type observation
!!
!! The routine is called by each filter process.
!! at the beginning of the analysis step before 
!! the loop through all local analysis domains.
!! 
!! It has to count the number of observations of the
!! observation type handled in this module according
!! to the current time step for all observations 
!! required for the analyses in the loop over all local 
!! analysis domains on the PE-local state domain.
!!
!! The following four variables have to be initialized in this routine
!! * thisobs\%doassim     - Whether to assimilate this type of observations
!! * thisobs\%disttype    - type of distance computation for localization with this observaton
!! * thisobs\%ncoord      - number of coordinates used for distance computation
!! * thisobs\%id_obs_p    - index of module-type observation in PE-local state vector
!!
!! Optional is the use of
!! * thisobs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
!! * thisobs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
!! * thisobs\%obs_err_type - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
!! * thisobs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
!!                          (default: 1=use global full observations)
!!
!! Further variables are set when the routine PDAFomi_gather_obs is called.
!!
  subroutine init_dim_obs_n_merged(step, dim_obs)

    use PDAF, &
         only: PDAFomi_gather_obs, PDAFomi_set_debug_flag
    use assim_pdaf_mod, &
         only: use_global_obs, cradius, sradius
    use fesom_pdaf, &
         only: mesh_fesom, nlmax, mydim_nod2d, pi, &
         month, day_in_month, yearnew, timenew, daynew
    use statevector_pdaf, &
         only: id, sfields
    use parallel_pdaf_mod, &
         only: MPI_SUM, MPIerr, COMM_filter, MPI_INTEGER, MPI_MAX
         
    use netcdf

    implicit none

! *** Arguments ***
    integer, intent(in)    :: step      !< Current time step
    integer, intent(inout) :: dim_obs   !< Dimension of full observation vector

! *** Local variables ***
    integer :: i, s, k, j, &
               e, e_reps, e_found, n_found, &
               e_new, e1, n1, e2, n2, &
               esa, nsa                       ! Counters
    character(len=2) :: mype_string           ! String for process rank
    character(len=4) :: year_string           ! String for yearly observation data path
    character(len=2) :: mon_string            ! String for daily observation files
    character(len=2) :: day_string            ! String for daily observation files
    
    character(len=300) :: filename
    logical            :: FileExists
    
    ! observation data including duplicates at repeatedly sampled elements
    integer :: dim_obs_p_reps                         ! Number of process local observations
    integer :: dim_obs_reps                           ! Number of global observations
    real, allocatable :: obs_p_reps(:)                ! PE-local observation values
    real, allocatable :: lon_p_reps(:),lat_p_reps(:)  ! PE-local observed coords
    integer, allocatable :: nod1_p_reps(:), &
                            nod2_p_reps(:), &
                            nod3_p_reps(:)            ! PE-local observed node/element indices
    integer, allocatable :: elem_p_reps(:)
    integer, allocatable :: nod1_g_reps(:), &
                            nod2_g_reps(:), &
                            nod3_g_reps(:)            ! Global observed node/element indices
    integer, allocatable :: elem_g_reps(:)
    integer, allocatable :: nlay_p_reps(:)
    real, allocatable    :: dep_p_reps(:)
    
    real    :: excl_relative_upper, excl_relative_lower ! relative upper and lower limits for exclusion
    
    !  Step 1: sorting of observations and trivial exclusion criteria
        ! observation data at element
    real, allocatable    :: x_p1(:,:), y_p1(:,:), z_p1(:,:)
    real, allocatable    :: obs_p_sort1(:,:)
    real, allocatable    :: dep_p_sort1(:,:)
        ! observation statistics at element
    integer, allocatable :: num_obs_p_sort1(:)
    real, allocatable    :: var_obs_p_sort1(:)
    real, allocatable    :: avg_obs_p_sort1(:)
    integer              :: dim_obs_p_sort1      ! PE-local number of observed elements
    integer              :: dim_obs_sort1        ! Global number of observed elements
        ! element indices
    integer, allocatable :: nod1_p_sort1(:), &
                            nod2_p_sort1(:), &
                            nod3_p_sort1(:)
    integer, allocatable :: elem_p_sort1(:)
    integer, allocatable :: nod1_g_sort1(:), &
                            nod2_g_sort1(:), &
                            nod3_g_sort1(:)
    integer, allocatable :: elem_g_sort1(:)
    integer, allocatable :: nlay_p_sort1(:)
    
    !  Step 2: exclude outliers from observations
    !          exclude based on model topography
        ! observation data at element
    real, allocatable    :: x_p2(:,:), y_p2(:,:), z_p2(:,:)
    real, allocatable    :: obs_p_sort2(:,:)
    real, allocatable    :: dep_p_sort2(:,:)
        ! observation statistics at element
    integer, allocatable :: num_obs_p_sort2(:)
    real, allocatable    :: var_obs_p_sort2(:)
    real, allocatable    :: avg_obs_p_sort2(:)
    integer              :: dim_obs_p_sort2
        ! element indices
    integer, allocatable :: nod1_p_sort2(:), &
                            nod2_p_sort2(:), &
                            nod3_p_sort2(:)
    integer, allocatable :: elem_p_sort2(:)
    integer, allocatable :: nod1_g_sort2(:), &
                            nod2_g_sort2(:), &
                            nod3_g_sort2(:)
    integer, allocatable :: elem_g_sort2(:)
    integer, allocatable :: nlay_p_sort2(:)
    
    !  Step 3: define observation error from variance
    !          exclude observations based on model-observation difference
    !          compute average for grid element
        ! observation data at element
    real, allocatable    :: x_p_sortavg(:), y_p_sortavg(:), z_p_sortavg(:)
    real, allocatable    :: lon_p_sortavg(:), lat_p_sortavg(:)
    real, allocatable    :: obs_p_sortavg(:)
    real, allocatable    :: dep_p_sortavg(:)
        ! observation statistics at element
    integer, allocatable :: num_obs_p_sortavg(:)
    real, allocatable    :: var_obs_p_sortavg(:)
        ! element indices
    integer, allocatable :: nod1_p_sortavg(:), &
                            nod2_p_sortavg(:), &
                            nod3_p_sortavg(:)
    integer, allocatable :: elem_p_sortavg(:)
    integer, allocatable :: nod1_g_sortavg(:), &
                            nod2_g_sortavg(:), &
                            nod3_g_sortavg(:)
    integer, allocatable :: elem_g_sortavg(:)
    integer, allocatable :: nlay_p_sortavg(:)
    
    ! final unique observation data, trimmed arrays
    integer :: dim_obs_p                      ! Number of process local observations
    
    real, allocatable :: obs_p(:)             ! PE-local observation values
    real, allocatable :: dep_p(:)
    real, allocatable :: var_obs_p(:)         ! PE-local observation variance
    real, allocatable :: ocoord_p(:,:)        ! PE-local coordinates of observations
    real, allocatable :: lon_p(:),lat_p(:)
    real, allocatable :: ivariance_obs_p(:)   ! PE-local array of inverse observation errors
    integer, allocatable :: nod1_p(:), &
                            nod2_p(:), &
                            nod3_p(:)         ! Array of observation pe-local indices on FESOM grid
    integer, allocatable :: elem_p(:)
    integer, allocatable :: nod1_g(:), &
                            nod2_g(:), &
                            nod3_g(:)         ! Array of observation global indices on FESOM grid
    integer, allocatable :: elem_g(:)
    integer, allocatable :: nlay_p(:)         ! Array of observation layer indices
    
    ! netCDF handles
    integer :: ncstat                         ! Status for NetCDF functions
    integer :: ncid, dimid                    ! NetCDF IDs
    integer :: id_obs, id_nod1, id_nod2, &
               id_nod3, id_nl, id_lon, &
               id_lat, id_depth, id_elem, &
               id_err
               
    ! exclusion count
    integer :: ecntex_topo_p, ecntex_topo   ! Observed Element Count Exclusion Topography
    integer :: ecntex_diff_p, ecntex_diff   ! Observed Element Count Exclusion Model-Observation-Difference
    integer :: ncntex_topo_p, ncntex_topo   ! Observation Count Exclusion Topography
    integer :: ncntex_diff_p, ncntex_diff   ! Observation Count Exclusion Model-Observation-Difference
    integer :: ncntex_onan_p, ncntex_onan   ! Observation Count Exclusion: Observation is FillValue
    integer :: ncntex_oneg_p, ncntex_oneg   ! Observation Count Exclusion: Observation is zero / negative
    integer :: ncntex_dnan_p, ncntex_dnan   ! Observation Count Exclusion: Depth is FillValue
    integer :: ncntex_dneg_p, ncntex_dneg   ! Observation Count Exclusion: Depth is zero / negative
    integer :: ncntex_outl_p, ncntex_outl   ! Observation Count Exclusion: Value is outlier
    integer :: ncntex_halo_p, ncntex_halo   ! Observation Count Exclusion: Node not on PE
    
    
    ! temporary variables during loops
    logical :: is_included = .true. ! observation exclusion due to trivial criteria
    logical :: is_unique            ! not a duplicate of previous observation
    integer :: nzmin                ! number of wet vertical model layers at observation location considering model topography
    integer :: num_obs_p_max        ! maximum number of observations at same element
    integer :: num_obs              ! number of observations at element
    real    :: num_obs_inv
    real, allocatable :: diffobs(:)
    real    :: emean                ! mean of observations at element
    real    :: evar                 ! variance of observations at element
    real    :: emean_fcst           ! observed model forecast at sample location
    
    real :: deg2rad
    real :: rad2deg
    
    deg2rad = PI/180.0
    rad2deg = 180.0/PI
    

! *******************
! *** Initialize  ***
! *******************

    ! Specify type of distance computation
    thisobs%disttype = 2   ! 2=Geographic

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2

    ! set localization radius, everywhere the same
    lradius_n_merged = 1.0e6 ! 1000 km
    sradius_n_merged = 1.0e6 ! 1000 km
    
    ! set limit factor for omitted observation due to high innovation
    ! omit observation if innovation larger than this factor times
    ! observation error (only active for >0)
    thisobs%inno_omit =  n_merged_exclude_diff / rms_obs_n_merged
    
    ! relative observation exclusion limits
    if (n_merged_excl_relative > 0.0) then
       excl_relative_lower = 1.0 - n_merged_excl_relative
       excl_relative_upper = 1.0 / excl_relative_lower
    endif

! **********************************
! *** Read PE-local observations ***
! **********************************

   ! ------------------------------------
   ! init file and number of observations
   ! ------------------------------------


    ! Daily files: initialize complete file name
    write(mype_string,'(i2.2)') mype_filter
    write(year_string,'(i4.4)') yearnew
    write(mon_string, '(i2.2)') month
    write(day_string, '(i2.2)') day_in_month
    
    file_n_merged = trim(year_string//'-'//mon_string//'-'//day_string//'.'//mype_string//'.nc')
    
    ! Message:
    if (step < 0) then
        ! ---
        ! Postprocessing
        ! ---
        isPP = .true.
        if (writepe) write (*,'(a,5x,a,1x,i2,a1,i2,a1,i4,1x,f5.2,a,1x,a)') &
            'FESOM-PDAF', 'Postprocessing of merged DIN observations at', &
            day_in_month, '.', month, '.', yearnew, timenew/3600.0,&
            'h; read from file:', file_n_merged
        thisobs%use_global_obs = 1
        thisobs%doassim = 1
    else
       ! ---
       ! Assimilation
       ! ---
       if (writepe) write (*,'(a,5x,a,1x,i2,a1,i2,a1,i4,1x,f5.2,a,1x,a)') &
            'FESOM-PDAF', 'Assimilate merged DIN observations at', &
            day_in_month, '.', month, '.', yearnew, timenew/3600.0,&
            'h; read from file:', file_n_merged
       ! Store whether to use global observations
       thisobs%use_global_obs = use_global_obs
       ! Store whether to assimilate this observation type
       if (assim_o_n_merged) thisobs%doassim = 1
    end if
    
    ! complete filename
    filename=trim(path_obs_n_merged)//year_string//'/dist72/'//trim(file_n_merged)
    inquire(file=trim(filename), exist=FileExists)
    
    if (FileExists) then
       ! Open the pe-local NetCDF file for that day
       ncstat = nf90_open(trim(filename), nf90_nowrite, ncid)
       if (ncstat /= nf90_noerr) then
        print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error opening NetCDF file'
       end if
       
       ! Get the number of observations dim_obs_p_reps
       ncstat = nf90_inq_dimid(ncid, 'INDEX', dimid)
       if (ncstat /= nf90_noerr) then
          print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error getting dimension ID'
       end if
       ncstat = nf90_inquire_dimension(ncid,dimid,len=dim_obs_p_reps)
       if (ncstat /= nf90_noerr) then
          print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error getting number of observations from NetCDF file'
       end if
    else
       dim_obs_p_reps=0
    endif ! FileExists
    
    ! ------------------------------------------
    ! init variables with zero or dim_obs_p_reps
    ! ------------------------------------------
    
    ! start exclusion count
    ecntex_topo_p = 0
    ecntex_diff_p = 0
    ncntex_onan_p = 0
    ncntex_oneg_p = 0
    ncntex_dnan_p = 0
    ncntex_dneg_p = 0
    ncntex_topo_p = 0
    ncntex_outl_p = 0
    ncntex_diff_p = 0
    ncntex_halo_p = 0
    
    ! &&&&&&&&&&&&&&&&&&&&&&&&&
    ! &&&&&&&&&&&&&&&&&&&&&&&&&
    if (dim_obs_p_reps <= 0) then
    ! --> Have no observations
    
      dim_obs_p = 0
      dim_obs_p_sort1 = 0
      num_obs_p_max = 0
    
      allocate(obs_p(1))
      allocate(ivariance_obs_p(1))
      allocate(ocoord_p(2, 1))
      allocate(thisobs%id_obs_p(6,1))
      
      obs_p=0.0
      ivariance_obs_p=1e-12
      ocoord_p=0.0
      thisobs%id_obs_p=0
      
      allocate(lon_p(0),lat_p(0))
      
      if (isPP) then
        allocate(thisobs_PP%nod1_g(0),thisobs_PP%nod2_g(0),thisobs_PP%nod3_g(0))
        allocate(thisobs_PP%elem_g(0))
        allocate(thisobs_PP%depth(0),thisobs_PP%nz(0))
        allocate(thisobs_PP%isExclObs(0))
        allocate(thisobs_PP%lon(0),thisobs_PP%lat(0))
        allocate(thisobs_PP%numrep(0))
      endif
    
    ! &&&&&&&&&&&&&&&&&&&&&&&&&
    ! &&&&&&&&&&&&&&&&&&&&&&&&&
    else ! (i.e. dim_obs_p > 0)
    ! --> Have observations
    
      ! Allocate memory before reading from netCDF
      allocate(obs_p_reps(dim_obs_p_reps))
      allocate(nod1_p_reps(dim_obs_p_reps),nod2_p_reps(dim_obs_p_reps),nod3_p_reps(dim_obs_p_reps))
      allocate(nlay_p_reps(dim_obs_p_reps))
      allocate(lon_p_reps(dim_obs_p_reps),lat_p_reps(dim_obs_p_reps))
      allocate(elem_p_reps(dim_obs_p_reps))
      allocate(dep_p_reps(dim_obs_p_reps))
      
      if (isPP) then
         allocate(nod1_g_reps(dim_obs_p_reps),nod2_g_reps(dim_obs_p_reps),nod3_g_reps(dim_obs_p_reps))
         allocate(elem_g_reps(dim_obs_p_reps))
      endif
      
     ! ---------
     ! read data
     ! ---------
            
      ! Reading observations
      ncstat = nf90_inq_varid(ncid,'OBS', id_obs)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error getting id_obs from netCDF'
      ncstat = nf90_get_var(ncid, id_obs, obs_p_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error reading obs from netCDF'
   
      ! Reading nodes
      ncstat = nf90_inq_varid(ncid,'NOD1_P', id_nod1)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error getting id_nod1 from netCDF'
      ncstat = nf90_get_var(ncid, id_nod1, nod1_p_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error reading nod1 from netCDF'
   
      ncstat = nf90_inq_varid(ncid,'NOD2_P', id_nod2)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error getting id_nod2 from netCDF'
      ncstat = nf90_get_var(ncid, id_nod2, nod2_p_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error reading nod2 from netCDF'
   
      ncstat = nf90_inq_varid(ncid,'NOD3_P', id_nod3)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error getting id_nod3 from netCDF'
      ncstat = nf90_get_var(ncid, id_nod3, nod3_p_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error reading nod3 from netCDF'
      
      ncstat = nf90_inq_varid(ncid,'ELEM_P', id_elem)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error getting id_elem from netCDF'
      ncstat = nf90_get_var(ncid, id_elem, elem_p_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error reading elem from netCDF'
   
      ! Reading layer
      ncstat = nf90_inq_varid(ncid,'NZ1', id_nl)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error getting id_layer from netCDF'
      ncstat = nf90_get_var(ncid, id_nl, nlay_p_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error reading layer from netCDF'
      
      ! Reading depth
      ncstat = nf90_inq_varid(ncid,'DEPTH', id_depth)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error getting id_depth from netCDF'
      ncstat = nf90_get_var(ncid, id_depth, dep_p_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error reading depth from netCDF'
   
      ! Reading coordinates
      ncstat = nf90_inq_varid(ncid,'LON', id_lon)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error getting id_lon from netCDF'
      ncstat = nf90_get_var(ncid, id_lon, lon_p_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error reading lon from netCDF'
   
      ncstat = nf90_inq_varid(ncid,'LAT', id_lat)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error getting id_lat from netCDF'
      ncstat = nf90_get_var(ncid, id_lat, lat_p_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error reading lat from netCDF'
      
      if (isPP) then
      
      ! Reading nodes on globe
      ncstat = nf90_inq_varid(ncid,'NOD1_G', id_nod1)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error getting id_nod1 from netCDF'
      ncstat = nf90_get_var(ncid, id_nod1, nod1_g_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error reading nod1 from netCDF'
   
      ncstat = nf90_inq_varid(ncid,'NOD2_G', id_nod2)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error getting id_nod2 from netCDF'
      ncstat = nf90_get_var(ncid, id_nod2, nod2_g_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error reading nod2 from netCDF'
   
      ncstat = nf90_inq_varid(ncid,'NOD3_G', id_nod3)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error getting id_nod3 from netCDF'
      ncstat = nf90_get_var(ncid, id_nod3, nod3_g_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error reading nod3 from netCDF'
      
      ! Reading mesh element
      ncstat = nf90_inq_varid(ncid,'ELEM_G', id_elem)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error getting id_elem from netCDF'
      ncstat = nf90_get_var(ncid, id_elem, elem_g_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error reading elem from netCDF'
      endif ! isPP
      
      ! **********************
      ! * Finding duplicates *
      ! **********************
      
      ! Pre-Step: Dry run of first loop to determine how much memory is required.
      allocate(num_obs_p_sort1(dim_obs_p_reps),source=0)
      allocate(elem_p_sort1   (dim_obs_p_reps),source=0)
      allocate(nlay_p_sort1   (dim_obs_p_reps),source=0)
      dim_obs_p_sort1 = 0  ! number of observed elements
      do e_reps=1, dim_obs_p_reps ! go through repreated observations one-by-one
         is_included = .true. ! trivial observation exclusion
         if ((obs_p_reps(e_reps) == -999) .or. &
             (obs_p_reps(e_reps) <=    0) .or. &
             (dep_p_reps(e_reps) == -999) .or. &
             (dep_p_reps(e_reps) <=    0) .or. &
             (nod1_p_reps(e_reps) <=   0) .or. &
             (nod2_p_reps(e_reps) <=   0) .or. &
             (nod3_p_reps(e_reps) <=   0)) then
             is_included = .false. 
         endif ! trivial observation exclusion
         if (is_included) then
            is_unique = .true.
            ! searching for previous occurence of elem(e_reps) and nlay(e_reps) within sorted arrays
            do e=1,dim_obs_p_sort1
               if ((elem_p_reps(e_reps)==elem_p_sort1(e)) .and. &
                   (nlay_p_reps(e_reps)==nlay_p_sort1(e))) then
                  is_unique = .false.
                  e_found   = e
                  exit ! --> found previous occurence at e; stop searching.
               endif
            enddo ! e=1,dim_obs_p_sort1
            if (.not. is_unique) then
            ! --> already have observations of elem(e_reps) at e_found
                num_obs_p_sort1(e_found) = num_obs_p_sort1(e_found) + 1 ! number of observations at e, adding one.
            else
            ! --> is the first or only observation at elem(e_reps)
               dim_obs_p_sort1 = dim_obs_p_sort1 + 1 ! number of observed elements, adding one.
               num_obs_p_sort1(dim_obs_p_sort1) = 1  ! number of observations at e, adding first.
               elem_p_sort1   (dim_obs_p_sort1) = elem_p_reps (e_reps) ! index of observed element
               nlay_p_sort1   (dim_obs_p_sort1) = nlay_p_reps (e_reps) ! index of observed element
            endif ! is_unique
         endif ! is_included
      enddo! e_reps=1, dim_obs_p_reps
      num_obs_p_max = maxval(num_obs_p_sort1)
      deallocate(num_obs_p_sort1,elem_p_sort1,nlay_p_sort1)
      
      endif
      ! print pre-count
      call MPI_Allreduce( &
                   dim_obs_p_reps, dim_obs_reps, &
                   1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
      call MPI_Allreduce( &
                   num_obs_p_max, n_found, &
                   1, MPI_INTEGER, MPI_MAX, COMM_filter, MPIerr)
      call MPI_Allreduce( &
                   dim_obs_p_sort1, dim_obs_sort1, &
                   1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
      if (mype_filter == 0) then
        write (*,'(a,5x,a30,2x,i7)') 'FESOM-PDAF', &
        '- DIN merged; dim_obs_reps:    ', dim_obs_reps
        write (*,'(a,5x,a30,2x,i7)') 'FESOM-PDAF', &
        '- DIN merged; num_obs_max:     ', n_found
        write (*,'(a,5x,a30,2x,i7)') 'FESOM-PDAF', &
        '- DIN merged; dim_obs_sort:    ', dim_obs_sort1
        write (*,'(a,5x,a30,2x,i7)') 'FESOM-PDAF', &
        '- DIN merged; num x dim_sort:  ', n_found * dim_obs_sort1
      endif
      if (dim_obs_p_reps > 0) then
      
      ! -----------------------------------------------
      !  Step 1: sorting and trivial exclusion criteria
      ! -----------------------------------------------
      
      allocate(x_p1       (dim_obs_p_sort1,num_obs_p_max),source=0.0)
      allocate(y_p1       (dim_obs_p_sort1,num_obs_p_max),source=0.0)
      allocate(z_p1       (dim_obs_p_sort1,num_obs_p_max),source=0.0)
      allocate(obs_p_sort1(dim_obs_p_sort1,num_obs_p_max),source=0.0)
      
      allocate(elem_p_sort1(dim_obs_p_sort1),source=0)
      allocate(nlay_p_sort1(dim_obs_p_sort1),source=0)
      allocate(nod1_p_sort1(dim_obs_p_sort1),source=0)
      allocate(nod2_p_sort1(dim_obs_p_sort1),source=0)
      allocate(nod3_p_sort1(dim_obs_p_sort1),source=0)
      
      allocate(num_obs_p_sort1(dim_obs_p_sort1),source=0)
      
      if (isPP) then
         allocate(dep_p_sort1(dim_obs_p_sort1,num_obs_p_max),source=0.0)
         
         allocate(elem_g_sort1(dim_obs_p_sort1),source=0)
         allocate(nod1_g_sort1(dim_obs_p_sort1),source=0)
         allocate(nod2_g_sort1(dim_obs_p_sort1),source=0)
         allocate(nod3_g_sort1(dim_obs_p_sort1),source=0)
      endif
      
      dim_obs_p_sort1 = 0
      e_found = 0
      n_found = 0
      
      ! go through repreated observations one-by-one:
      do_elem_p_reps: do e_reps=1, dim_obs_p_reps
         
         ! trivial observation exclusion
         is_included = .true.
         ! FillValue for observation
         if     (obs_p_reps(e_reps) == -999) then
            is_included = .false.
            ncntex_onan_p = ncntex_onan_p +1
         ! Negative / zero observation value
         elseif (obs_p_reps(e_reps) <=    0) then
            is_included = .false.
            ncntex_oneg_p = ncntex_oneg_p +1
         ! FillValue for depth
         elseif (dep_p_reps(e_reps) == -999) then
            is_included = .false.
            ncntex_dnan_p = ncntex_dnan_p +1
         ! Negative / zero depth
         elseif (dep_p_reps(e_reps) <=    0) then
            is_included = .false.
            ncntex_dneg_p = ncntex_dneg_p +1
         ! Invalid Node
         elseif (nod1_p_reps(e_reps) <=   0) then
            is_included = .false.
            ncntex_halo_p = ncntex_halo_p +1
         elseif (nod2_p_reps(e_reps) <=   0) then
            is_included = .false.
            ncntex_halo_p = ncntex_halo_p +1
         elseif (nod3_p_reps(e_reps) <=   0) then
            is_included = .false.
            ncntex_halo_p = ncntex_halo_p +1 
         endif ! trivial observation exclusion
         
         if (is_included) then
            is_unique = .true.
            
            ! searching for previous occurence of elem(e_reps) and nlay(e_reps) within sorted arrays:
            do_elem_p: do e=1,dim_obs_p_sort1
               
               if ((elem_p_reps(e_reps)==elem_p_sort1(e)) .and. &
                   (nlay_p_reps(e_reps)==nlay_p_sort1(e))) then
                  is_unique = .false.
                  e_found   = e
                  exit ! --> found previous occurence at e; stop searching.
               endif
            enddo do_elem_p ! counter e
            
            if (.not. is_unique) then
            ! --> already have observations of elem(e_reps) at e_found:
                
                ! number of observations at e, adding one:
                n_found = num_obs_p_sort1(e_found) + 1
                num_obs_p_sort1(e_found) = n_found
                
                ! writing observation to array at e:
                obs_p_sort1(e_found,n_found) = obs_p_reps(e_reps)
                
                ! transforming to cartesian coordinates:
                x_p1(e_found,n_found) = cos(deg2rad*lat_p_reps(e_reps)) * cos(deg2rad*lon_p_reps(e_reps))
                y_p1(e_found,n_found) = cos(deg2rad*lat_p_reps(e_reps)) * sin(deg2rad*lon_p_reps(e_reps))
                z_p1(e_found,n_found) = sin(deg2rad*lat_p_reps(e_reps))
                            
                ! writing depth to array at e:
                if (isPP) dep_p_sort1(e_found,n_found) = dep_p_reps(e_reps)
               
            else
            ! --> is the first or only observation at elem(e_reps)
               
               ! number of observed elements, adding one:
               dim_obs_p_sort1 = dim_obs_p_sort1 + 1
               e_new = dim_obs_p_sort1
               
               ! one observation so far:
               num_obs_p_sort1(e_new) = 1
               
               ! indices of observed element:
               elem_p_sort1 (e_new) = elem_p_reps (e_reps)
               nlay_p_sort1 (e_new) = nlay_p_reps (e_reps)
               nod1_p_sort1 (e_new) = nod1_p_reps (e_reps)
               nod2_p_sort1 (e_new) = nod2_p_reps (e_reps)
               nod3_p_sort1 (e_new) = nod3_p_reps (e_reps)
               
               ! observation data:
               obs_p_sort1 (e_new,1) = obs_p_reps (e_reps)
               ! transforming to cartesian coordinates:
               x_p1 (e_new,1) = cos(deg2rad*lat_p_reps(e_reps)) * cos(deg2rad*lon_p_reps(e_reps))
               y_p1 (e_new,1) = cos(deg2rad*lat_p_reps(e_reps)) * sin(deg2rad*lon_p_reps(e_reps))
               z_p1 (e_new,1) = sin(deg2rad*lat_p_reps(e_reps))
               
               if (isPP) then
                  ! indices
                  elem_g_sort1 (e_new) = elem_g_reps(e_reps)
                  nod1_g_sort1 (e_new) = nod1_g_reps(e_reps)
                  nod2_g_sort1 (e_new) = nod2_g_reps(e_reps)
                  nod3_g_sort1 (e_new) = nod3_g_reps(e_reps)
                  ! data
                  dep_p_sort1 (e_new,1) = dep_p_reps(e_reps)
               endif ! isPP

            endif ! is_unique
         endif ! is_included
      enddo do_elem_p_reps ! counter e_reps
      
      ! clean up *_reps
      deallocate(nod1_p_reps,nod2_p_reps,nod3_p_reps,lon_p_reps,lat_p_reps,obs_p_reps,elem_p_reps,nlay_p_reps,dep_p_reps)
      if (isPP) then
         deallocate(nod1_g_reps,nod2_g_reps,nod3_g_reps)
         deallocate(elem_g_reps)
      endif
      
      
      ! -----------------------------------------------
      !  Step 2: exclude outliers from observations
      !          exclude based on model topography
      ! -----------------------------------------------
      
      ! maximum number of observations at same element
      num_obs_p_max = maxval(num_obs_p_sort1)
      
      allocate(x_p2       (dim_obs_p_sort1,num_obs_p_max),source=0.0)
      allocate(y_p2       (dim_obs_p_sort1,num_obs_p_max),source=0.0)
      allocate(z_p2       (dim_obs_p_sort1,num_obs_p_max),source=0.0)
      allocate(obs_p_sort2(dim_obs_p_sort1,num_obs_p_max),source=0.0)
      
      allocate(elem_p_sort2(dim_obs_p_sort1),source=0)
      allocate(nlay_p_sort2(dim_obs_p_sort1),source=0)
      allocate(nod1_p_sort2(dim_obs_p_sort1),source=0)
      allocate(nod2_p_sort2(dim_obs_p_sort1),source=0)
      allocate(nod3_p_sort2(dim_obs_p_sort1),source=0)
      
      allocate(num_obs_p_sort2(dim_obs_p_sort1),source=0)
      
      if (isPP) then
         allocate(dep_p_sort2(dim_obs_p_sort1,num_obs_p_max),source=0.0)
         allocate(elem_g_sort2(dim_obs_p_sort1),source=0)
         allocate(nod1_g_sort2(dim_obs_p_sort1),source=0)
         allocate(nod2_g_sort2(dim_obs_p_sort1),source=0)
         allocate(nod3_g_sort2(dim_obs_p_sort1),source=0)
      endif
      
      e2 = 0      
      ! loop observed elements
      do e1=1, dim_obs_p_sort1
      
         ! model topography
         ! nlevels <- index of first dry layer at node
         ! nzmin   <- index of first incomplete element
         nzmin = min (mesh_fesom% nlevels_nod2D( nod1_p_sort1(e1) ), &
                      mesh_fesom% nlevels_nod2D( nod2_p_sort1(e1) ), &
                      mesh_fesom% nlevels_nod2D( nod3_p_sort1(e1) ))
         if (nlay_p_sort1(e1) < nzmin) then
         ! --> included after model topography check
             
            ! observed elements, adding one:
            e2 = e2+1
            
            ! write to arrays
            elem_p_sort2 (e2) = elem_p_sort1 (e1)
            nlay_p_sort2 (e2) = nlay_p_sort1 (e1)
            nod1_p_sort2 (e2) = nod1_p_sort1 (e1)
            nod2_p_sort2 (e2) = nod2_p_sort1 (e1)
            nod3_p_sort2 (e2) = nod3_p_sort1 (e1)
            
            if (isPP) then
              elem_g_sort2 (e2) = elem_g_sort1 (e1)
              nod1_g_sort2 (e2) = nod1_g_sort1 (e1)
              nod2_g_sort2 (e2) = nod2_g_sort1 (e1)
              nod3_g_sort2 (e2) = nod3_g_sort1 (e1)
            endif
            
            ! compute mean and variance of observations at e1
            num_obs     = num_obs_p_sort1(e1)
            num_obs_inv = 1/real(num_obs)
            allocate(diffobs(num_obs))
            emean   = num_obs_inv * sum (obs_p_sort1(e1,1:num_obs))
            diffobs = obs_p_sort1(e1,1:num_obs) - emean
            evar    = num_obs_inv * sum (diffobs * diffobs)
            
            n2 = 0
            do n1=1, num_obs
               if ((evar <= 1.0) .or. (diffobs(n1)*diffobs(n1) <= (9.0 * evar))) then
               ! --> included after outlier check
               
                   ! number of observations at element, adding one:
                   n2 = n2+1
                   ! write to arrays
                   obs_p_sort2(e2,n2) = obs_p_sort1(e1,n1)
                   x_p2(e2,n2) = x_p1(e1,n1)
                   y_p2(e2,n2) = y_p1(e1,n1)
                   z_p2(e2,n2) = z_p1(e1,n1)
                   if (isPP) dep_p_sort2(e2,n2) = dep_p_sort1(e1,n1)
               
               else
                   ! exlusion counter
                   ncntex_outl_p = ncntex_outl_p +1
               endif ! outlier check
               ! number of observations after outlier exclusions
               num_obs_p_sort2(e2) = n2
            enddo
            
            deallocate(diffobs)
         else
            ! exclusion counter
            ecntex_topo_p = ecntex_topo_p +1
            ncntex_topo_p = ncntex_topo_p +num_obs_p_sort1(e1)
         endif ! model topography check
         ! number of observed elements after topography exclusions
         dim_obs_p_sort2 = e2
      enddo ! loop observed elements
      
      ! ecntex_topo_p = dim_obs_p_sort1 - dim_obs_p_sort2
      
      ! clean up *_sort1
      deallocate(x_p1       )
      deallocate(y_p1       )
      deallocate(z_p1       )
      deallocate(obs_p_sort1)
      deallocate(elem_p_sort1)
      deallocate(nlay_p_sort1)
      deallocate(nod1_p_sort1)
      deallocate(nod2_p_sort1)
      deallocate(nod3_p_sort1)
      deallocate(num_obs_p_sort1)
      if (isPP) then
         deallocate(dep_p_sort1)
         deallocate(elem_g_sort1)
         deallocate(nod1_g_sort1)
         deallocate(nod2_g_sort1)
         deallocate(nod3_g_sort1)
      endif
      
      ! -------------------------------------------------------------------
      !  Step 3: exclude observations based on model-observation difference
      !          define observation error from variance
      !          compute average for grid element
      ! -------------------------------------------------------------------
      
      allocate(obs_p_sortavg    (dim_obs_p_sort2))
      allocate(var_obs_p_sortavg(dim_obs_p_sort2))
      allocate(num_obs_p_sortavg(dim_obs_p_sort2))
      allocate(elem_p_sortavg(dim_obs_p_sort2))
      allocate(nlay_p_sortavg(dim_obs_p_sort2))
      allocate(nod1_p_sortavg(dim_obs_p_sort2),nod2_p_sortavg(dim_obs_p_sort2),nod3_p_sortavg(dim_obs_p_sort2))
      allocate(x_p_sortavg(dim_obs_p_sort2),y_p_sortavg(dim_obs_p_sort2),z_p_sortavg(dim_obs_p_sort2))
      allocate(lat_p_sortavg(dim_obs_p_sort2),lon_p_sortavg(dim_obs_p_sort2))

      obs_p_sortavg(:)     = 0
      var_obs_p_sortavg(:) = 0
      num_obs_p_sortavg(:) = 0
      elem_p_sortavg(:) = 0
      nlay_p_sortavg(:) = 0
      nod1_p_sortavg(:) = 0
      nod2_p_sortavg(:) = 0
      nod3_p_sortavg(:) = 0
      x_p_sortavg(:) = 0
      y_p_sortavg(:) = 0
      z_p_sortavg(:) = 0
      lat_p_sortavg(:) = 0
      lon_p_sortavg(:) = 0
      
      if (isPP) then
         allocate(elem_g_sortavg(dim_obs_p_sort2))
         allocate(nod1_g_sortavg(dim_obs_p_sort2),nod2_g_sortavg(dim_obs_p_sort2),nod3_g_sortavg(dim_obs_p_sort2))
         allocate(dep_p_sortavg(dim_obs_p_sort2))
         elem_g_sortavg(:)=0.0
         nod1_g_sortavg(:) = 0
         nod2_g_sortavg(:) = 0
         nod3_g_sortavg(:) = 0
         dep_p_sortavg(:) = 0
      endif ! isPP
      
      ! -------------------------------------------------------------------------
      if ((n_merged_excl_absolute > 0) .and. (n_merged_excl_relative > 0)) then
      ! perform loop, checking absolute and relative exclusion criteria
         
         ! loop observed elements
         esa=0
         do e2=1, dim_obs_p_sort2
           ! model forecast at observed element
           emean_fcst = irefdens * third * &
                        (  mean_n_p((nlmax) * (nod1_p_sort2(e2)-1) + nlay_p_sort2(e2)) &
                        +  mean_n_p((nlmax) * (nod2_p_sort2(e2)-1) + nlay_p_sort2(e2)) &
                        +  mean_n_p((nlmax) * (nod3_p_sort2(e2)-1) + nlay_p_sort2(e2)))
           ! loop observations at element
           num_obs = num_obs_p_sort2(e2)
           nsa=0
           do n2=1,num_obs
              if (abs(obs_p_sort2(e2,n2) - emean_fcst) < n_merged_excl_absolute) then
              ! observation passed absolute model difference check
                if (obs_p_sort2(e2,n2) > (excl_relative_lower * emean_fcst)) then
                ! observation passed relative difference check, lower limit
                  if (obs_p_sort2(e2,n2) < (excl_relative_upper * emean_fcst)) then
                  ! observation passed relative difference check, upper limit
                     
                     if (nsa==0) then
                     ! first observation at element
                       esa = esa+1
                       
                       ! write indices to arrays
                       elem_p_sortavg (esa) = elem_p_sort2 (e2)
                       nlay_p_sortavg (esa) = nlay_p_sort2 (e2)
                       nod1_p_sortavg (esa) = nod1_p_sort2 (e2)
                       nod2_p_sortavg (esa) = nod2_p_sort2 (e2)
                       nod3_p_sortavg (esa) = nod3_p_sort2 (e2)
                       
                       if (isPP) then
                         elem_g_sortavg (esa) = elem_g_sort2 (e2)
                         nod1_g_sortavg (esa) = nod1_g_sort2 (e2)
                         nod2_g_sortavg (esa) = nod2_g_sort2 (e2)
                         nod3_g_sortavg (esa) = nod3_g_sort2 (e2)
                       endif
                       
                       ! compute variance for observation error
                       num_obs_inv = 1/real(num_obs)
                       allocate(diffobs(num_obs))
                       emean   = num_obs_inv * sum (obs_p_sort2(e2,1:num_obs))
                       diffobs = obs_p_sort2(e2,1:num_obs) - emean
                       var_obs_p_sortavg(esa) = num_obs_inv * sum (diffobs * diffobs)
                       deallocate(diffobs)
                     endif ! nsa==0
                     
                     nsa=nsa+1
                     ! include this observation into sum of observations at element
                     obs_p_sortavg(esa) = obs_p_sortavg(esa) + obs_p_sort2(e2,n2)
                     x_p_sortavg  (esa) = x_p_sortavg  (esa) + x_p2(e2,n2)
                     y_p_sortavg  (esa) = y_p_sortavg  (esa) + y_p2(e2,n2)
                     z_p_sortavg  (esa) = z_p_sortavg  (esa) + z_p2(e2,n2)
                     
                     if (isPP) dep_p_sortavg(esa) = dep_p_sortavg(esa) + dep_p_sort2(e2,n2)
                  endif ! upper limit check
                endif ! lower limit check
              endif ! absolute difference check
           enddo ! loop observations
           
           ! exclusion count
           ncntex_diff_p = ncntex_diff_p +num_obs -nsa
           
           if (nsa>0) then
           ! after checks, at least one valid observation at element
           
             num_obs_p_sortavg(esa) = nsa
             num_obs_inv = 1 / real(nsa)
             ! mean of observations at element
             obs_p_sortavg(esa) = num_obs_inv * obs_p_sortavg(esa)
             ! mean location, in cartesian coordinates
             x_p_sortavg(esa) = num_obs_inv * x_p_sortavg(esa)
             y_p_sortavg(esa) = num_obs_inv * y_p_sortavg(esa)
             z_p_sortavg(esa) = num_obs_inv * z_p_sortavg(esa)
             ! mean depth
             if (isPP) dep_p_sortavg(esa) = num_obs_inv * dep_p_sortavg(esa)        
             ! transforming to spherical coordinates
             lat_p_sortavg(esa) = atan2( z_p_sortavg(esa), sqrt(x_p_sortavg(esa)*x_p_sortavg(esa) + y_p_sortavg(esa)*y_p_sortavg(esa)) )
             lon_p_sortavg(esa) = atan2( y_p_sortavg(esa), x_p_sortavg(esa))
           else
             ! exclusion count
             ecntex_diff_p = ecntex_diff_p +1
           endif ! nsa>0
         enddo ! loop observed elements
         dim_obs_p = esa
         
      ! ------------------------------------------------------------------------------
      elseif ((n_merged_excl_absolute > 0) .and. (n_merged_excl_relative == 0)) then
      ! perform loop, checking only absolute exclusion criterion
      
         ! loop observed elements
         esa=0
         do e2=1, dim_obs_p_sort2
           ! model forecast at observed element
           emean_fcst = irefdens * third * &
                        (  mean_n_p((nlmax) * (nod1_p_sort2(e2)-1) + nlay_p_sort2(e2)) &
                        +  mean_n_p((nlmax) * (nod2_p_sort2(e2)-1) + nlay_p_sort2(e2)) &
                        +  mean_n_p((nlmax) * (nod3_p_sort2(e2)-1) + nlay_p_sort2(e2)))
           ! loop observations at element
           num_obs = num_obs_p_sort2(e2)
           nsa=0
           do n2=1,num_obs
              if (abs(obs_p_sort2(e2,n2) - emean_fcst) < n_merged_excl_absolute) then
              ! observation passed absolute model difference check
                     
                     if (nsa==0) then
                     ! first observation at element
                       esa = esa+1
                       
                       ! write indices to arrays
                       elem_p_sortavg (esa) = elem_p_sort2 (e2)
                       nlay_p_sortavg (esa) = nlay_p_sort2 (e2)
                       nod1_p_sortavg (esa) = nod1_p_sort2 (e2)
                       nod2_p_sortavg (esa) = nod2_p_sort2 (e2)
                       nod3_p_sortavg (esa) = nod3_p_sort2 (e2)
                       
                       if (isPP) then
                         elem_g_sortavg (esa) = elem_g_sort2 (e2)
                         nod1_g_sortavg (esa) = nod1_g_sort2 (e2)
                         nod2_g_sortavg (esa) = nod2_g_sort2 (e2)
                         nod3_g_sortavg (esa) = nod3_g_sort2 (e2)
                       endif
                       
                       ! compute variance for observation error
                       num_obs_inv = 1/real(num_obs)
                       allocate(diffobs(num_obs))
                       emean   = num_obs_inv * sum (obs_p_sort2(e2,1:num_obs))
                       diffobs = obs_p_sort2(e2,1:num_obs) - emean
                       var_obs_p_sortavg(esa) = num_obs_inv * sum (diffobs * diffobs)
                       deallocate(diffobs)
                     endif ! nsa==0
                     
                     nsa=nsa+1
                     ! include this observation into sum of observations at element
                     obs_p_sortavg(esa) = obs_p_sortavg(esa) + obs_p_sort2(e2,n2)
                     x_p_sortavg  (esa) = x_p_sortavg  (esa) + x_p2(e2,n2)
                     y_p_sortavg  (esa) = y_p_sortavg  (esa) + y_p2(e2,n2)
                     z_p_sortavg  (esa) = z_p_sortavg  (esa) + z_p2(e2,n2)
                     
                     if (isPP) dep_p_sortavg(esa) = dep_p_sortavg(esa) + dep_p_sort2(e2,n2)
              endif ! absolute difference check
           enddo ! loop observations
           
           ! exclusion count
           ncntex_diff_p = ncntex_diff_p +num_obs -nsa
           
           if (nsa>0) then
           ! after checks, at least one valid observation at element
           
             num_obs_p_sortavg(esa) = nsa
             num_obs_inv = 1 / real(nsa)
             ! mean of observations at element
             obs_p_sortavg(esa) = num_obs_inv * obs_p_sortavg(esa)
             ! mean location, in cartesian coordinates
             x_p_sortavg(esa) = num_obs_inv * x_p_sortavg(esa)
             y_p_sortavg(esa) = num_obs_inv * y_p_sortavg(esa)
             z_p_sortavg(esa) = num_obs_inv * z_p_sortavg(esa)
             ! mean depth
             if (isPP) dep_p_sortavg(esa) = num_obs_inv * dep_p_sortavg(esa)        
             ! transforming to spherical coordinates
             lat_p_sortavg(esa) = atan2( z_p_sortavg(esa), sqrt(x_p_sortavg(esa)*x_p_sortavg(esa) + y_p_sortavg(esa)*y_p_sortavg(esa)) )
             lon_p_sortavg(esa) = atan2( y_p_sortavg(esa), x_p_sortavg(esa))
           else
             ! exclusion count
             ecntex_diff_p = ecntex_diff_p +1
           endif ! nsa>0
         enddo ! loop observed elements
         dim_obs_p = esa
      
      ! ------------------------------------------------------------------------------
      elseif ((n_merged_excl_absolute == 0) .and. (n_merged_excl_relative > 0)) then
      ! perform loop, checking only relative exclusion criterion
      
         ! loop observed elements
         esa=0
         do e2=1, dim_obs_p_sort2
           ! model forecast at observed element
           emean_fcst = irefdens * third * &
                        (  mean_n_p((nlmax) * (nod1_p_sort2(e2)-1) + nlay_p_sort2(e2)) &
                        +  mean_n_p((nlmax) * (nod2_p_sort2(e2)-1) + nlay_p_sort2(e2)) &
                        +  mean_n_p((nlmax) * (nod3_p_sort2(e2)-1) + nlay_p_sort2(e2)))
           ! loop observations at element
           num_obs = num_obs_p_sort2(e2)
           nsa=0
           do n2=1,num_obs
                if (obs_p_sort2(e2,n2) > (excl_relative_lower * emean_fcst)) then
                ! observation passed relative difference check, lower limit
                  if (obs_p_sort2(e2,n2) < (excl_relative_upper * emean_fcst)) then
                  ! observation passed relative difference check, upper limit
                     
                     if (nsa==0) then
                     ! first observation at element
                       esa = esa+1
                       
                       ! write indices to arrays
                       elem_p_sortavg (esa) = elem_p_sort2 (e2)
                       nlay_p_sortavg (esa) = nlay_p_sort2 (e2)
                       nod1_p_sortavg (esa) = nod1_p_sort2 (e2)
                       nod2_p_sortavg (esa) = nod2_p_sort2 (e2)
                       nod3_p_sortavg (esa) = nod3_p_sort2 (e2)
                       
                       if (isPP) then
                         elem_g_sortavg (esa) = elem_g_sort2 (e2)
                         nod1_g_sortavg (esa) = nod1_g_sort2 (e2)
                         nod2_g_sortavg (esa) = nod2_g_sort2 (e2)
                         nod3_g_sortavg (esa) = nod3_g_sort2 (e2)
                       endif
                       
                       ! compute variance for observation error
                       num_obs_inv = 1/real(num_obs)
                       allocate(diffobs(num_obs))
                       emean   = num_obs_inv * sum (obs_p_sort2(e2,1:num_obs))
                       diffobs = obs_p_sort2(e2,1:num_obs) - emean
                       var_obs_p_sortavg(esa) = num_obs_inv * sum (diffobs * diffobs)
                       deallocate(diffobs)
                     endif ! nsa==0
                     
                     nsa=nsa+1
                     ! include this observation into sum of observations at element
                     obs_p_sortavg(esa) = obs_p_sortavg(esa) + obs_p_sort2(e2,n2)
                     x_p_sortavg  (esa) = x_p_sortavg  (esa) + x_p2(e2,n2)
                     y_p_sortavg  (esa) = y_p_sortavg  (esa) + y_p2(e2,n2)
                     z_p_sortavg  (esa) = z_p_sortavg  (esa) + z_p2(e2,n2)
                     
                     if (isPP) dep_p_sortavg(esa) = dep_p_sortavg(esa) + dep_p_sort2(e2,n2)
                  endif ! upper limit check
                endif ! lower limit check
           enddo ! loop observations
           
           ! exclusion count
           ncntex_diff_p = ncntex_diff_p +num_obs -nsa
           
           if (nsa>0) then
           ! after checks, at least one valid observation at element
           
             num_obs_p_sortavg(esa) = nsa
             num_obs_inv = 1 / real(nsa)
             ! mean of observations at element
             obs_p_sortavg(esa) = num_obs_inv * obs_p_sortavg(esa)
             ! mean location, in cartesian coordinates
             x_p_sortavg(esa) = num_obs_inv * x_p_sortavg(esa)
             y_p_sortavg(esa) = num_obs_inv * y_p_sortavg(esa)
             z_p_sortavg(esa) = num_obs_inv * z_p_sortavg(esa)
             ! mean depth
             if (isPP) dep_p_sortavg(esa) = num_obs_inv * dep_p_sortavg(esa)        
             ! transforming to spherical coordinates
             lat_p_sortavg(esa) = atan2( z_p_sortavg(esa), sqrt(x_p_sortavg(esa)*x_p_sortavg(esa) + y_p_sortavg(esa)*y_p_sortavg(esa)) )
             lon_p_sortavg(esa) = atan2( y_p_sortavg(esa), x_p_sortavg(esa))
           else
             ! exclusion count
             ecntex_diff_p = ecntex_diff_p +1
           endif ! nsa>0
         enddo ! loop observed elements
         dim_obs_p = esa
      
      ! ------------------------------------------
      else
      ! perform loop, no exclusion checks required
      
         ! loop observed elements
         do e2=1, dim_obs_p_sort2
           ! no exclusions
           esa=e2
           num_obs_p_sortavg(esa) = num_obs_p_sort2(e2)
           num_obs = num_obs_p_sortavg(esa)
           num_obs_inv = 1/real(num_obs)
           ! write indices to arrays
           elem_p_sortavg (esa) = elem_p_sort2 (e2)
           nlay_p_sortavg (esa) = nlay_p_sort2 (e2)
           nod1_p_sortavg (esa) = nod1_p_sort2 (e2)
           nod2_p_sortavg (esa) = nod2_p_sort2 (e2)
           nod3_p_sortavg (esa) = nod3_p_sort2 (e2)
           if (isPP) then
             elem_g_sortavg (esa) = elem_g_sort2 (e2)
             nod1_g_sortavg (esa) = nod1_g_sort2 (e2)
             nod2_g_sortavg (esa) = nod2_g_sort2 (e2)
             nod3_g_sortavg (esa) = nod3_g_sort2 (e2)
           endif
           ! mean of observations at element
           emean   = num_obs_inv * sum (obs_p_sort2(e2,1:num_obs))
           obs_p_sortavg(esa) = emean
           ! variance for observation error
           allocate(diffobs(num_obs))
           diffobs = obs_p_sort2(e2,1:num_obs) - emean
           var_obs_p_sortavg(esa) = num_obs_inv * sum (diffobs * diffobs)
           deallocate(diffobs)
           ! mean depth
           if (isPP) dep_p_sortavg(esa) = num_obs_inv * sum (dep_p_sort2(e2,1:num_obs))
           ! mean location, in cartesian coordinates
           x_p_sortavg(esa) = num_obs_inv * sum (x_p2(e2,1:num_obs))
           y_p_sortavg(esa) = num_obs_inv * sum (y_p2(e2,1:num_obs))
           z_p_sortavg(esa) = num_obs_inv * sum (z_p2(e2,1:num_obs))
           ! transforming to spherical coordinates
           lat_p_sortavg(esa) = atan2( z_p_sortavg(esa), sqrt(x_p_sortavg(esa)*x_p_sortavg(esa) + y_p_sortavg(esa)*y_p_sortavg(esa)) )
           lon_p_sortavg(esa) = atan2( y_p_sortavg(esa), x_p_sortavg(esa))
         enddo ! loop observed elements
         dim_obs_p = dim_obs_p_sort2
      
      endif ! absolute / relative exclusion criteria
      
      ! clean up *_sort2
      deallocate(x_p2)
      deallocate(y_p2)
      deallocate(z_p2)
      deallocate(obs_p_sort2)
      deallocate(elem_p_sort2)
      deallocate(nlay_p_sort2)
      deallocate(nod1_p_sort2)
      deallocate(nod2_p_sort2)
      deallocate(nod3_p_sort2)
      deallocate(num_obs_p_sort2)
      if (isPP) then
         deallocate(dep_p_sort2)
         deallocate(elem_g_sort2)
         deallocate(nod1_g_sort2)
         deallocate(nod2_g_sort2)
         deallocate(nod3_g_sort2)
      endif

      ! -------------------------------------------------------------------
      !  Final: unique observation data, trimmed arrays
      ! -------------------------------------------------------------------
      
      ! copy data from *_sortavg to trimmed arrays
      
      allocate(obs_p(dim_obs_p))                                      ! PE-local observation values
      allocate(var_obs_p(dim_obs_p))                                  ! PE-local variance for observation error
      allocate(lon_p(dim_obs_p),lat_p(dim_obs_p))
      allocate(nod1_p(dim_obs_p),nod2_p(dim_obs_p),nod3_p(dim_obs_p)) ! PE-local indices on FESOM grid
      allocate(elem_p(dim_obs_p))
      allocate(nlay_p(dim_obs_p))                                     ! PE-local layer indices
      
      obs_p(:)     = obs_p_sortavg    (1:dim_obs_p)
      var_obs_p(:) = var_obs_p_sortavg(1:dim_obs_p)
      lon_p(:)     = lon_p_sortavg    (1:dim_obs_p)
      lat_p(:)     = lat_p_sortavg    (1:dim_obs_p)
      nod1_p(:)    = nod1_p_sortavg   (1:dim_obs_p)
      nod2_p(:)    = nod2_p_sortavg   (1:dim_obs_p)
      nod3_p(:)    = nod3_p_sortavg   (1:dim_obs_p)
      elem_p(:)    = elem_p_sortavg   (1:dim_obs_p)
      nlay_p(:)    = nlay_p_sortavg   (1:dim_obs_p)
      
      deallocate(obs_p_sortavg    )
      deallocate(var_obs_p_sortavg)
      deallocate(lon_p_sortavg    )
      deallocate(lat_p_sortavg    )
      deallocate(nod1_p_sortavg   )
      deallocate(nod2_p_sortavg   )
      deallocate(nod3_p_sortavg   )
      deallocate(elem_p_sortavg   )
      deallocate(nlay_p_sortavg   )
      
      ! copy data for postprocessing
      if (isPP) then
          allocate(elem_g(dim_obs_p))
          allocate(nod1_g(dim_obs_p),nod2_g(dim_obs_p),nod3_g(dim_obs_p))
          allocate(dep_p(dim_obs_p))
         
          nod1_g(:) = nod1_g_sortavg(1:dim_obs_p)
          nod2_g(:) = nod2_g_sortavg(1:dim_obs_p)
          nod3_g(:) = nod3_g_sortavg(1:dim_obs_p)
          elem_g(:) = elem_g_sortavg(1:dim_obs_p)
          dep_p(:)  = dep_p_sortavg (1:dim_obs_p)
          
          deallocate(nod1_g_sortavg)
          deallocate(nod2_g_sortavg)
          deallocate(nod3_g_sortavg)
          deallocate(elem_g_sortavg)
          deallocate(dep_p_sortavg )

          allocate(thisobs_PP%nod1_g(dim_obs_p),thisobs_PP%nod2_g(dim_obs_p),thisobs_PP%nod3_g(dim_obs_p))
          allocate(thisobs_PP%elem_g(dim_obs_p))
          allocate(thisobs_PP%depth(dim_obs_p),thisobs_PP%nz(dim_obs_p))
          allocate(thisobs_PP%lon(dim_obs_p),thisobs_PP%lat(dim_obs_p))
          allocate(thisobs_PP%numrep(dim_obs_p))
          
          thisobs_PP%nod1_g = nod1_g
          thisobs_PP%nod2_g = nod2_g
          thisobs_PP%nod3_g = nod3_g
          thisobs_PP%elem_g = elem_g
          thisobs_PP%nz     = nlay_p
          
          thisobs_PP%lon = lon_p / pi * 180.0
          thisobs_PP%lat = lat_p / pi * 180.0
          
          thisobs_PP%depth  = dep_p
          thisobs_PP%numrep = num_obs_p_sortavg(:dim_obs_p)
          
          allocate(thisobs_PP%isExclObs(dim_obs_p))
          thisobs_PP%isExclObs = 0 ! already trimmed: no more to exclude
      endif
      deallocate(x_p_sortavg,y_p_sortavg,z_p_sortavg,num_obs_p_sortavg)
      
      ! set inverse observation error variance
      allocate(ivariance_obs_p(dim_obs_p))
      do i = 1, dim_obs_p
         ! combine variance of duplicate observations (var_obs_p, computed above)
         ! and minimum observation error (rms_obs_n_merged, defined in namelist)
         ivariance_obs_p(i) = 1.0 / (var_obs_p(i) + (rms_obs_n_merged * rms_obs_n_merged))
      enddo

      ! set observations coordinates
      allocate(ocoord_p(2,dim_obs_p))
      ocoord_p(1,:) = lon_p
      ocoord_p(2,:) = lat_p
      
      ! *** Initialize index vector of observed surface nodes ***
      ! This array has as many rows as required for the observation operator
      ! 1 if observations are at grid points; >1 if interpolation is required
      allocate(thisobs%id_obs_p(6,dim_obs_p))
      do i = 1, dim_obs_p
        ! observed tracer
        thisobs%id_obs_p(val1,i) = (nlmax) * (nod1_p(i)-1) + nlay_p(i) + sfields(id%DIN)%off
        thisobs%id_obs_p(val2,i) = (nlmax) * (nod2_p(i)-1) + nlay_p(i) + sfields(id%DIN)%off
        thisobs%id_obs_p(val3,i) = (nlmax) * (nod3_p(i)-1) + nlay_p(i) + sfields(id%DIN)%off
        ! density for unit conversion
        thisobs%id_obs_p(dens1,i) = (nlmax) * (nod1_p(i)-1) + nlay_p(i) + sfields(id%sigma)%off
        thisobs%id_obs_p(dens2,i) = (nlmax) * (nod2_p(i)-1) + nlay_p(i) + sfields(id%sigma)%off
        thisobs%id_obs_p(dens3,i) = (nlmax) * (nod3_p(i)-1) + nlay_p(i) + sfields(id%sigma)%off
      end do
      
    endif ! IF (dim_obs_p = 0) ELSEIF (dim_obs_p > 0)
    ! &&&&&&&&&&&&&&&&&&&&&&&&&
    ! &&&&&&&&&&&&&&&&&&&&&&&&&
    ! regardless of whether we have observations or not,
    ! continue with what needs to be done 
    

    ! **************************************
    ! *** Gather full observation arrays ***
    ! **************************************

    call PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivariance_obs_p, ocoord_p, &
                            thisobs%ncoord, lradius_n_merged, dim_obs)
    
    ! Global inverse variance array (thisobs%ivar_obs_f)
    ! has been gathered during PDAFomi call.
    ! In case of coupled DA (multiple sweeps),
    ! we will reset ivar in the SUBROUTINE init_dim_obs_l_n_merged,
    ! depending on whether the observation type is to be assimilated during the sweep.
    ! Thus, save a copy, from which we can reset ivar.
    if (n_sweeps>1) then
       if (allocated(ivariance_obs_g)) deallocate(ivariance_obs_g)
       allocate(ivariance_obs_g(dim_obs))
       ivariance_obs_g = thisobs%ivar_obs_f
    end if
    
    ! Gather global observation exclusion statistics
    call MPI_Allreduce( &
                       dim_obs_p_reps, dim_obs_reps, &
                       1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
    call MPI_Allreduce( &
                       dim_obs_p_sort1, dim_obs_sort1, &
                       1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
    call MPI_Allreduce( &
                       ncntex_onan_p, ncntex_onan, &
                       1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
    call MPI_Allreduce( &
                       ncntex_oneg_p, ncntex_oneg, &
                       1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
    call MPI_Allreduce( &
                       ncntex_dnan_p, ncntex_dnan, &
                       1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
    call MPI_Allreduce( &
                       ncntex_dneg_p, ncntex_dneg, &
                       1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
    call MPI_Allreduce( &
                       ncntex_halo_p, ncntex_halo, &
                       1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
    call MPI_Allreduce( &
                       ncntex_topo_p, ncntex_topo, &
                       1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
    call MPI_Allreduce( &
                       ncntex_outl_p, ncntex_outl, &
                       1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
    call MPI_Allreduce( &
                       ncntex_diff_p, ncntex_diff, &
                       1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
    call MPI_Allreduce( &
                       ecntex_topo_p, ecntex_topo, &
                       1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
    call MPI_Allreduce( &
                       ecntex_diff_p, ecntex_diff, &
                       1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)

    ! Print observation exclusion statistics                     
    if (mype_filter == 0) then
    
        write (*,'(a,5x,a61,2x,i7)') 'FESOM-PDAF', &
        '- DIN merged. Number of samples, some of them duplicates      ', dim_obs_reps
        
        write (*,'(a,5x,a61,2x,i7)') 'FESOM-PDAF', &
        '- DIN merged. Samples excluded (observation is FillValue)     ', ncntex_onan
        
        write (*,'(a,5x,a61,2x,i7)') 'FESOM-PDAF', &
        '- DIN merged. Samples excluded (observation is zero/ negative)', ncntex_oneg
        
        write (*,'(a,5x,a61,2x,i7)') 'FESOM-PDAF', &
        '- DIN merged. Samples excluded (depth is FillValue)           ', ncntex_dnan
        
        write (*,'(a,5x,a61,2x,i7)') 'FESOM-PDAF', &
        '- DIN merged. Samples excluded (depth is zero/ negative)      ', ncntex_dneg
        
        write (*,'(a,5x,a61,2x,i7)') 'FESOM-PDAF', &
        '- DIN merged. Samples excluded (nod invalid)                  ', ncntex_halo
        
        write (*,'(a,5x,a61,2x,i7)') 'FESOM-PDAF', &
        '- DIN merged. Samples excluded (model topography)             ', ncntex_topo
        
        write (*,'(a,5x,a61,2x,i7)') 'FESOM-PDAF', &
        '- DIN merged. Samples excluded (outliers among samples)       ', ncntex_outl
        
        write (*,'(a,5x,a61,2x,i7)') 'FESOM-PDAF', &
        '- DIN merged. Samples excluded (model-observation difference) ', ncntex_diff
        
        write (*,'(a,5x,a61,2x,i7)') 'FESOM-PDAF', &
        '- DIN merged. Number of observed volumes:                     ', dim_obs_sort1
        
        write (*,'(a,5x,a61,2x,i7)') 'FESOM-PDAF', &
        '- DIN merged. Volumes excluded (model topography)             ', ecntex_topo
        
        write (*,'(a,5x,a61,2x,i7)') 'FESOM-PDAF', &
        '- DIN merged. Volumes excluded (model-observation difference) ', ecntex_diff
    endif

    ! *** Clean-up ***
    if (FileExists) then
       ncstat = nf90_close(ncid)
       if (ncstat /= nf90_noerr) then
          print *, 'FESOM-PDAF - obs_n_merged_pdafomi - Error closing NetCDF file'
       end if
    endif
    
    deallocate(ivariance_obs_p,ocoord_p,obs_p)
    deallocate(lat_p,lon_p)
    if (dim_obs_p_reps > 0) then
          deallocate(var_obs_p, nod1_p, nod2_p, nod3_p, elem_p, nlay_p)
          if (isPP) deallocate(nod1_g, nod2_g, nod3_g, elem_g, dep_p)
    endif

  end subroutine init_dim_obs_n_merged



!-------------------------------------------------------------------------------
!> Implementation of observation operator 
!!
!! This routine applies the full observation operator
!! for the type of observations handled in this module.
!!
!! One can choose a proper observation operator from
!! PDAFOMI_OBS_OP or add one to that module or 
!! implement another observation operator here.
!!
!! The routine is called by all filter processes.
!!
  subroutine obs_op_n_merged(dim_p, dim_obs, state_p, ostate)

    use PDAF, &
         only: PDAFomi_obs_op_gridavg, &
               PDAFomi_gather_obsstate

    implicit none

! *** Arguments ***
    integer, intent(in) :: dim_p                 !< PE-local state dimension
    integer, intent(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
    real, intent(in)    :: state_p(dim_p)        !< PE-local model state
    real, intent(inout) :: ostate(dim_obs)       !< Full observed state
    
    real, allocatable   :: ostate_p(:)           !< Pe-local observed state
    integer :: i                                 !< Counters

! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    if (thisobs%doassim == 1) then

       if (thisobs%dim_obs_p>0) then
       ! have obs
          allocate(ostate_p(thisobs%dim_obs_p))
          
          if (isPP) then
            do i = 1, thisobs%dim_obs_p
                ! -- unit conversion:
                !    from milli mol per m3 (model) --> micro mol per kg (observations)
                ! -- average values of 3 grid points
                ostate_p(i) =  ( state_p(thisobs%id_obs_p(val1,i)) * irefdens &
                               + state_p(thisobs%id_obs_p(val2,i)) * irefdens &
                               + state_p(thisobs%id_obs_p(val3,i)) * irefdens &
                               ) * third
            end do
          else
            ! initialize observed pe-local state vector
            do i = 1, thisobs%dim_obs_p
                ! -- unit conversion:
                !    from milli mol per m3 (model) --> micro mol per kg (observations)
                ! -- average values of 3 grid points
                   ostate_p(i) =  ( state_p(thisobs%id_obs_p(val1,i)) / state_p(thisobs%id_obs_p(dens1,i)) &
                                  + state_p(thisobs%id_obs_p(val2,i)) / state_p(thisobs%id_obs_p(dens2,i)) &
                                  + state_p(thisobs%id_obs_p(val3,i)) / state_p(thisobs%id_obs_p(dens3,i)) &
                                  ) * third
            end do
          endif ! isPP
                    
       else
       ! habe no obs
          allocate(ostate_p(1))
       end if

       ! *** Global: Gather full observed state vector
       call PDAFomi_gather_obsstate(thisobs, ostate_p, ostate)
       
       ! clean up
       deallocate(ostate_p)

       ! For profile observations handled here, the observation
       ! operator has to average the values of 3 grid points.
       ! For this the observation operator OBS_OP_F_GRIDAVG is used.
       ! CALL PDAFomi_obs_op_gridavg(thisobs, 3, state_p, ostate)
    
    end if ! (thisobs%doassim == 1)

  end subroutine obs_op_n_merged



!-------------------------------------------------------------------------------
!> Initialize local information on the module-type observation
!!
!! The routine is called during the loop over all local
!! analysis domains. It has to initialize the information
!! about local observations of the module type. It returns
!! number of local observations of the module type for the
!! current local analysis domain in DIM_OBS_L and the full
!! and local offsets of the observation in the overall
!! observation vector.
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type and  local analysis domain.
!!
  subroutine init_dim_obs_l_n_merged(domain_p, step, dim_obs, dim_obs_l)

    ! Include PDAFomi function
    use PDAF, only: PDAFomi_init_dim_obs_l

    ! Include localization radius and local coordinates
    use assim_pdaf_mod, only: coords_l, locweight, loctype
    
    ! Number of domains per sweep:
    use fesom_pdaf, only: myDim_nod2D

    implicit none

! *** Arguments ***
    integer, intent(in)  :: domain_p     !< Index of current local analysis domain: containing repititive sweeps
    integer, intent(in)  :: step         !< Current time step
    integer, intent(in)  :: dim_obs      !< Full dimension of observation vector
    integer, intent(inout) :: dim_obs_l  !< Local dimension of observation vector

    if (thisobs%doassim == 1) then
    
       ! ************************************************************
       ! *** Adapt observation error for coupled DA (multi sweep) ***
       ! ************************************************************
    
       if (n_sweeps>1) then
       
          ! Physics observations sweep.
          ! Set inverse observation error to small value
          if (domain_p==1) then

             if (mype_filter==0) &
                  write (*,'(a,4x,a)') 'FESOM-PDAF', &
                   '--- PHY sweep: set ivar_obs_f for merged DIN to 1.0e-12'
             thisobs%ivar_obs_f = 1.0e-12
             
          ! BGC observations sweep.
          elseif (domain_p==myDim_nod2D+1) then
             if (mype_filter==0) &
                  write (*,'(a,4x,a)') 'FESOM-PDAF', &
                  '--- BIO sweep: set ivar_obs_f for merged DIN to original'
             thisobs%ivar_obs_f(:) = ivariance_obs_g
          end if
       end if ! n_sweeps

       ! **********************************************
       ! *** Initialize local observation dimension ***
       ! **********************************************
       call PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
            locweight, lradius_n_merged, sradius_n_merged, dim_obs_l)
    
     end if ! thisobs%doassim
  end subroutine init_dim_obs_l_n_merged

end module obs_n_merged_pdafomi
