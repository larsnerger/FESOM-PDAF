!> Include FESOM-RECOM variables and initialize FESOM grid information
!!
!! This module includes variables from FESOM-RECOM. For PDAF in its
!! online coupling it is the single point which directly links
!! to FESOM-RECOM. All other routines include from this module.
!! For the offline case the variables are declared here and
!! separately initialized.
!!
MODULE fesom_pdaf

  USE MOD_MESH, ONLY: t_mesh
  USE g_parsup, ONLY: MPI_COMM_FESOM, &
       myDim_nod2D, &          !Process-local number of 2D nodes  
       myDim_nod2D, &          ! Process-local number of vertices
       myDim_elem2D, &         ! Process-local number of elements 
       myDim_edge2D, &
       myList_nod2D, myList_edge2D, &
       eDim_nod2D, eDim_elem2D 
  USE g_comm_auto, &
       ONLY: gather_nod
  USE o_arrays, &
       ONLY: eta_n, uv, wvel, tr_arr, unode, MLD1, MLD2, sigma0, &
       zbar_n_bot, zbar_n_srf, hnode_new, z_n, zbar_n
  USE REcoM_GloVar, &
       ONLY: GloPCO2surf, GloCO2flux, Diags3D, PAR3D, export, &
       PistonVelocity, alphaCO2
  USE i_arrays, &
       ONLY: a_ice
  use g_clock, &
       only: timenew, timeold, daynew, dayold, yearnew, yearold, &
       month, day_in_month, num_day_in_month, cyearnew, cyearold, &
       fleapyear, check_fleapyr, dt, clock
  use g_events, &
       only: daily_event, monthly_event
  USE g_comm_auto, &
       ONLY: exchange_nod, exchange_elem
  USE g_rotate_grid, &
       ONLY: r2g                 ! Transform from the mesh (rotated) coordinates to geographical coordinates  
  use o_param, &
       only: r_earth, rad
  use recom_config, &
       only: tiny_chl, tiny, chl2N_max, chl2N_max_d, NCmax, &      
       NCmax_d, SiCmax, Redfield, SecondsPerDay
  USE g_config, &
       ONLY: runid, step_per_day, ResultPath
  USE g_sbf, &
       ONLY: atmdata, &
       i_xwind, i_ywind, i_humi, &
       i_qsr, i_qlw, i_tair, i_prec, i_mslp, i_snow

  IMPLICIT NONE

  REAL, PARAMETER ::  pi=3.141592653589793   ! Pi

  ! FESOM mesh:
  TYPE(t_mesh), POINTER, SAVE :: mesh_fesom
  INTEGER, PARAMETER :: nlmax = 46            ! CORE2 mesh: deepest wet cells at mesh_fesom%nl-2
  REAL, ALLOCATABLE :: topography3D(:,:)      ! topography: 1 for wet nodes and 0 for dry nodes (array shape as in model)
  REAL, ALLOCATABLE :: topography_p(:)        ! """                                             (array shape as state_p)
  REAL, ALLOCATABLE :: topography3D_g(:,:)    ! """                                             (array shape as in model globally)
  REAL :: area_surf_glob(nlmax)               ! ocean area and standard volume to calculate area-/volume weighted means
  REAL :: inv_area_surf_glob(nlmax)
  REAL :: volo_full_glob, inv_volo_full_glob
  REAL, ALLOCATABLE :: cellvol(:,:)           ! standard volume of cells, NOT considering time-varying ALE layerwidth

END MODULE fesom_pdaf
