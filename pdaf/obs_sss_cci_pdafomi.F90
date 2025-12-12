!> PDAF-OMI observation module 
!!
!! This module handles operations for one data type (called 'module-type' below).
!! 
!! Module type here: Surface salinity observations from SMOS.
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
module obs_sss_cci_pdafomi

  use parallel_pdaf_mod, &
       only: mype_filter
  use PDAF, &                     ! Declaration of observation data types
       only: obs_f, obs_l
  use coupled_da_mod, &           ! Variables for coupled DA
       only: n_sweeps 

  implicit none
  save

  ! Variables which are inputs to the module (usually set in init_pdaf)
  logical :: assim_o_sss_cci      ! Whether to assimilate SSS data

  ! Further variables specific for the SMOS SSS observations
  character(len=100) :: path_obs_sss_cci  = ''      ! Path to observations
  character(len=110) :: file_sss_cci_prefix  = ''   ! file name prefix for observations 
  character(len=110) :: file_sss_cci_suffix  = '.nc'! file name suffix for observations 
!~   CHARACTER(len=110) :: file_syntobs_sss_cci = 'syntobs_sss_cci.nc' ! File name for synthetic observations

  real    :: rms_obs_sss_cci      ! Observation error standard deviation
  real    :: bias_obs_sss_cci     ! Observation bias

  real    :: lradius_sss_cci      ! Localization radius in the ocean
  real    :: sradius_sss_cci      ! Support radius for localization function

  logical :: sss_cci_exclude_ice  ! Whether to exclude observations at grid points with ice
  real    :: sss_cci_exclude_diff ! Limit difference beyond which observations are excluded (0.0 to deactivate)
  logical :: sss_cci_fixed_rmse   ! Whether to use a fixed RMS error or the error provided with the data

!~   REAL, ALLOCATABLE :: mean_ice_p (:)    ! Mean ice concentration for observation exclusion: IMPORT FROM SST obs. module!
  real, allocatable :: mean_sss_cci_p (:)    ! Mean value for observation exclusion
  real, allocatable :: loc_radius_sss_cci(:) ! Localization radius array
  real, allocatable :: ivariance_obs_g(:)    ! global-earth inverse observation variances

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
  subroutine init_dim_obs_sss_cci(step, dim_obs)

    use PDAF, &
         only: PDAFomi_gather_obs
    use assim_pdaf_mod, &
         only: twin_experiment, use_global_obs, delt_obs_ocn, &
         cradius, sradius
    use fesom_pdaf, &
         only: mesh_fesom, nlmax, r2g, mydim_nod2d, &
         month, day_in_month, yearnew, timenew, daynew
    use statevector_pdaf, &
         only: id, sfields
    use parallel_pdaf_mod, &
         only: MPI_SUM, MPIerr, COMM_filter, MPI_INTEGER
    use obs_sst_pdafomi, &
         only: mean_ice_p

    implicit none

    include 'netcdf.inc'

! *** Arguments ***
    integer, intent(in)    :: step       !< Current time step
    integer, intent(inout) :: dim_obs    !< Dimension of full observation vector

! *** Local variables ***
    integer :: i, iter_file, i_obs       ! Counters
    integer :: dim_obs_p                 ! number of PE-local observations
    integer :: fileid                    ! ID for NetCDF file
    integer :: id_state, id_std          ! ID for state
    integer :: stat(100)                 ! Status for NetCDF functions
    integer :: startv(2),countv(2)       ! Vectors for reading fields
    character(len=5)   :: mype_string    ! String for process rank
    character(len=100) :: obs_file = ''     ! Complete name of observation file without path
    real(4), allocatable :: all_obs_p(:)    ! PE-local complete observation field read from file
    real(4), allocatable :: all_std_p(:)    ! PE-local complete observation error field read from file
    real, allocatable :: obs_error_p(:)     ! PE-local observation error
    real, allocatable :: obs_p(:)           ! PE-local observed field
    real, allocatable :: ivariance_obs_p(:) ! PE-local inverse observation error variance
    real, allocatable :: ocoord_n2d_p(:,:)  ! PE-local coordinates of observations
    integer :: status                       ! Status flag for PDAF gather operation
    integer :: cnt_ex_ice_p, cnt_ex_diff_p  ! PE-local counts of excluded points due to ice and difference
    integer :: cnt_ex_ice, cnt_ex_diff      ! Global counts of excluded points due to ice and difference
    real, allocatable :: obs_g(:)           ! Global full observation vector (used in case of limited obs.)
    real, allocatable :: ivar_obs_g(:)      ! Global full inverse variances (used in case of limited obs.)
    real, allocatable :: ocoord_g(:,:)      ! Global full observation coordinates (used in case of limited obs.)
!~     INTEGER, ALLOCATABLE :: id_obs_p_all(:)        ! State vector index of all observations on process domain (including excluded observations)
    integer, allocatable :: obs_include_index(:)   ! Index of observed (not excluded!) surface nodes on process domain
    integer :: dim_obs_f                      ! Global full observation number


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    ! Store whether to assimilate this observation type
    if (assim_o_sss_cci) thisobs%doassim = 1

    ! Specify type of distance computation
    thisobs%disttype = 2   ! 2=Geographic

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2

    ! Initialize flag for type of full observations
    if (mype_filter==0) write(*,*) 'obs_sss_cci_pdafomi: use_global_obs', use_global_obs
    thisobs%use_global_obs = use_global_obs

    ! set localization radius
    lradius_sss_cci = cradius
    sradius_sss_cci = sradius
    if (allocated(loc_radius_sss_cci)) deallocate(loc_radius_sss_cci)
    allocate(loc_radius_sss_cci(mydim_nod2d))
    loc_radius_sss_cci(:) = lradius_sss_cci


! **********************************
! *** Read PE-local observations ***
! **********************************

    ! Initialize complete file name
    write(mype_string,'(i5.5)') mype_filter

    obs_file=trim(file_sss_cci_prefix)//trim(mype_string)//trim(file_sss_cci_suffix)

    ! Allocate array
    allocate(all_obs_p(myDim_nod2D))
    ! ALLOCATE(all_std_p(myDim_nod2D))

    ! Position to read from file
    iter_file = daynew
    
    ! Debugging message:
    if (mype_filter==0) then
       write (*,'(a,5x,a,i2,a,i2,a,i4,a,f5.2,a,i)') &
            'FESOM-PDAF', 'Assimilate SSS CCI observations - OBS_SSS_CCI at ', &
            day_in_month, '.', month, '.', yearnew, ' ', timenew/3600.0,&
            ' h; read at ', iter_file
    end if

    ! Read observation and standard deviation
         
    stat(1) = NF_OPEN(trim(path_obs_sss_cci)//trim(obs_file), NF_NOWRITE, fileid)        
 
    ! *** Read state estimate ***

    stat(2) = NF_INQ_VARID(fileid, 'obs', id_state)

    startv(2) = iter_file
    countv(2) = 1
    startv(1) = 1
    countv(1) = myDim_nod2D 

    stat(3) = NF_GET_VARA_REAL (fileid, id_state, startv, countv, all_obs_p)
 
    ! *** Read standard deviation ***

    ! stat(4) = NF_INQ_VARID(fileid, 'std', id_std)

    ! startv(2) = iter_file
    ! countv(2) = 1
    ! startv(1) = 1
    ! countv(1) = myDim_nod2D 

    ! stat(5) = NF_GET_VARA_REAL (fileid, id_std, startv, countv, all_std_p)

    ! *** close file  ***
    stat(6) = NF_CLOSE(fileid)
 
    ! check status flag
    do i=1,6
       if (stat(i).ne.NF_NOERR) write(*,*) &
            'NetCDF error in reading full SSS_CCI, no.',i, &
            ' file ',obs_file
       stat(i)=0
    end do


! ****************************
! *** Exclude observations ***
! ****************************

    ! *** Exclude observations if mean_ice is not zero ***

    exclude_ice: if (sss_cci_exclude_ice) then

       cnt_ex_ice_p = 0
       do i = 1, myDim_nod2D
          if ((mean_ice_p(i) > 0.0) &
             .and. &
              (abs(all_obs_p(i))<999.0)) &
          then
             all_obs_p(i) = 1.0e6
             cnt_ex_ice_p = cnt_ex_ice_p + 1
          end if
       end do
       
       ! *** Sum PE-local excluded nodes to global number of nodes excluded ***

       call MPI_Allreduce(cnt_ex_ice_p, cnt_ex_ice, 1, MPI_INTEGER, MPI_SUM, &
            COMM_filter, MPIerr)

       if (mype_filter == 0) &
            write (*,'(a,5x,a,i7)') 'FESOM-PDAF', &
            '--- Observations excluded because of ice', cnt_ex_ice

       ! *** Set localization radius to zero for grid points with ice ***

       if (mype_filter == 0) &
            write (*,'(a,5x,a,i7)') 'FESOM-PDAF', &
            '--- Set localization radius to zero for points with ice'
       do i = 1, myDim_nod2D
          if (mean_ice_p (i) > 0.0) then
             loc_radius_sss_cci(i) = 0.0
          end if
       end do

    end if exclude_ice


    ! *** Exclude observations if difference from ensemble mean is beyond limit SSS_EXCLUDE_DIFF ***

    exclude_diff: if (sss_cci_exclude_diff > 0.0) then

       cnt_ex_diff_p = 0
       do i = 1, myDim_nod2D
          if (abs(mean_sss_cci_p(i) - all_obs_p(i)) > sss_cci_exclude_diff .and. abs(all_obs_p(i))<=999.0) then
             all_obs_p(i) = 1.0e6
             cnt_ex_diff_p = cnt_ex_diff_p+1
          end if
       end do

       call MPI_Allreduce(cnt_ex_diff_p, cnt_ex_diff, 1, MPI_INTEGER, MPI_SUM, &
            COMM_filter, MPIerr)

       if (mype_filter == 0) &
            write (*,'(a,5x,a,f6.2,a,i7)') 'FESOM-PDAF', &
            '--- CCI Observations excluded due to difference >',sss_cci_exclude_diff,'psu:', cnt_ex_diff

    end if exclude_diff
  

! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

    ! *** Count PE-local number of observations ***
    dim_obs_p = 0
    do i = 1, myDim_nod2d
       if (abs(all_obs_p(i)) < 999.0) dim_obs_p=dim_obs_p+1
    enddo

    haveobs: if(dim_obs_p>0) then
    
		! *** Initialize index vector of observed surface nodes ***
		! This array has a many rows as required for the observation operator
		! 1 if observations are at grid points; >1 if interpolation is required
		
		allocate(thisobs%id_obs_p(1, dim_obs_p))
		allocate(obs_include_index(dim_obs_p))

		i_obs=0
		do i = 1, myDim_nod2d
		   if (abs(all_obs_p(i)) < 999.0) then
			  i_obs = i_obs + 1
			  
			  ! index for state vector
			  thisobs%id_obs_p(1, i_obs) = &
			  (i-1) * (nlmax) + 1 + sfields(id% salt)%off
			  
			  ! index for all_obs_p and surface nod2d vector, respectively. 
			  obs_include_index(i_obs) = i
		   endif
		enddo

		! *** Initialize PE-local vectors of observations and std error ***
		allocate(obs_p(dim_obs_p))
		allocate(obs_error_p(dim_obs_p))
		
		do i = 1, dim_obs_p
		   obs_p(i)       = real(all_obs_p(obs_include_index(i)), 8)
		   ! obs_error_p(i) = REAL(all_std_p(obs_include_index(i)), 8)
		enddo

		! *** Initialize coordinate arrays for PE-local observations
		allocate(ocoord_n2d_p(2, dim_obs_p))
		do i = 1, dim_obs_p
		   ! Rotate to geographic coordinates and store
		   call r2g(ocoord_n2d_p(1, i), ocoord_n2d_p(2, i), &
				mesh_fesom%coord_nod2d(1, obs_include_index(i)), mesh_fesom%coord_nod2d(2, obs_include_index(i)))
		enddo


	! ****************************************************************
	! *** Define observation errors for process-local observations ***
	! ****************************************************************

		if (sss_cci_fixed_rmse) then

		   ! *** Set constant observation error *** 
		   if (mype_filter == 0) &
				write (*, '(a, 5x, a, f12.3, a)') 'FESOM-PDAF', &
				'--- Use global SSS CCI observation error of ', rms_obs_sss_cci, 'psu'

		   obs_error_p(:) = rms_obs_sss_cci
		!  ELSE

		   !! *** Use variable error from file
		   !IF (mype_filter == 0) &
		   !		WRITE (*,'(a,5x,a,i7)') 'FESOM-PDAF', &
		   !		'--- Use variable SSS CCI observation error from file'
		end if

		! Set inverse observation error variance
		allocate(ivariance_obs_p(dim_obs_p))
		do i = 1, dim_obs_p
		   ivariance_obs_p(i) = 1.0 / obs_error_p(i)**2
		end do


	! ****************************
	! *** De-bias observations ***
	! ****************************

		if (bias_obs_sss_cci/=0.0 .and. mype_filter == 0) &
			 write (*, '(a, 5x, a, f12.3)') &
			 'FESOM-PDAF', '--- For SSS, use global observation bias of ', bias_obs_sss_cci

		obs_p = obs_p - bias_obs_sss_cci
		
	else ! (i.e. dim_obs_p==0)
	
		allocate(obs_p(1))
		allocate(ivariance_obs_p(1))
		allocate(ocoord_n2d_p(2, 1))
		allocate(thisobs%id_obs_p(1,1))
		thisobs%id_obs_p = 0
		ivariance_obs_p  = 1e-12
		
		allocate(obs_include_index(1))
		allocate(obs_error_p(1))
				
	endif haveobs

! **************************************
! *** Gather full observation arrays ***
! **************************************

    call PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivariance_obs_p, ocoord_n2d_p, &
         thisobs%ncoord, lradius_sss_cci, dim_obs)

    ! Global inverse variance array (thisobs%ivar_obs_f)
    ! has been gathered, but, in case of coupled DA / "double-sweep",
    ! it will be reset during each sweep. Thus, save a copy:
    if (n_sweeps>1) then
       if (allocated(ivariance_obs_g)) deallocate(ivariance_obs_g)
       allocate(ivariance_obs_g(dim_obs))
       ivariance_obs_g = thisobs%ivar_obs_f
    end if

! *********************************************************
! *** For twin experiment: Read synthetic observations  ***
! *********************************************************

    ! This only works with use_full_obs=.true.
!~     IF (twin_experiment) THEN
!~        CALL read_syn_obs(file_syntobs_sss_cci, dim_obs, thisobs%obs_f, 0, 1-mype_filter)
!~     END IF


! ********************
! *** Finishing up ***
! ********************

    ! Clean up arrays
    deallocate(all_obs_p, obs_error_p) ! all_std_p
    deallocate(obs_p, ocoord_n2d_p, ivariance_obs_p)
    if (allocated(obs_include_index)) deallocate(obs_include_index)

  end subroutine init_dim_obs_sss_cci



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
  subroutine obs_op_sss_cci(dim_p, dim_obs, state_p, ostate)

    use PDAF, &
         only: PDAFomi_obs_op_gridpoint

    implicit none

! *** Arguments ***
    integer, intent(in) :: dim_p                 !< PE-local state dimension
    integer, intent(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
    real, intent(in)    :: state_p(dim_p)        !< PE-local model state
    real, intent(inout) :: ostate(dim_obs)       !< Full observed state


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    if (thisobs%doassim == 1) then
       call PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)
    end if

  end subroutine obs_op_sss_cci



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
  subroutine init_dim_obs_l_sss_cci(domain_p, step, dim_obs, dim_obs_l)

    ! Include PDAFomi function
    use PDAF, only: PDAFomi_init_dim_obs_l
    ! Include routine for adaptive localization radius
    use adaptive_lradius_pdaf, only: get_adaptive_lradius_pdaf
    ! Include localization radius and local coordinates
    use assim_pdaf_mod, only: coords_l, locweight, loctype
    ! Number of domains per sweep:
    use fesom_pdaf, only: myDim_nod2D

    implicit none

! *** Arguments ***
    integer, intent(in)  :: domain_p     !< Index of current local analysis domain
    integer, intent(in)  :: step         !< Current time step
    integer, intent(in)  :: dim_obs      !< Full dimension of observation vector
    integer, intent(inout) :: dim_obs_l  !< Local dimension of observation vector


    if (thisobs%doassim == 1) then
       if (loctype == 1) then
          ! *** Variable localization radius for fixed effective observation dimension ***
          call get_adaptive_lradius_pdaf(thisobs, modulo(domain_p,myDim_nod2D), lradius_sss_cci, loc_radius_sss_cci)
       end if
       lradius_sss_cci = loc_radius_sss_cci(modulo(domain_p,myDim_nod2D))


       ! ************************************************************
       ! *** Adapt observation error for coupled DA (double loop) ***
       ! ************************************************************
    
       if (n_sweeps>1) then
       
          ! Physics observations sweep.
          if (domain_p==1) then
             if (mype_filter==0) &
                  write (*,'(a,4x,a)') 'FESOM-PDAF', &
                   '--- PHY sweep: leave ivar_obs_f for SSS-CCI as it is'
          ! BGC observations sweep.
          elseif (domain_p==myDim_nod2D+1) then
             if (mype_filter==0) &
                  write (*,'(a,4x,a)') 'FESOM-PDAF', &
                  '--- BIO sweep: set ivar_obs_f for SSS-CCI to 1.0e-12'
             thisobs%ivar_obs_f = 1.0e-12
          end if
       end if ! n_sweeps

       ! **********************************************
       ! *** Initialize local observation dimension ***
       ! **********************************************


       call PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
            locweight, lradius_sss_cci, sradius_sss_cci, dim_obs_l)
       
    end if

  end subroutine init_dim_obs_l_sss_cci

end module obs_sss_cci_pdafomi
