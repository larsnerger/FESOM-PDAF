!> Include FESOM-RECOM variables and initialize FESOM grid information
!!
!! This module includes variables from FESOM-RECOM. For PDAF in its
!! online coupling it is the single point which directly links
!! to FESOM-RECOM. All other routines include from this module.
!! For the offline case the variables are declared here and
!! separately initialized.
!!
!! __Revision history:__
!! 2025-12 - Lars Nerger - Initial code for PDAF3 revision
!!
module fesom_pdaf

  use MOD_MESH, only: t_mesh
  use g_parsup, only: MPI_COMM_FESOM, &
       myDim_nod2D, &          !Process-local number of 2D nodes  
       myDim_nod2D, &          ! Process-local number of vertices
       myDim_elem2D, &         ! Process-local number of elements 
       myDim_edge2D, &
       myList_nod2D, myList_edge2D, &
       eDim_nod2D, eDim_elem2D 
  use g_comm_auto, &
       only: gather_nod
  use o_arrays, &
       only: eta_n, uv, wvel, tr_arr, unode, MLD1, MLD2, sigma0, &
       zbar_n_bot, zbar_n_srf, hnode_new, z_n, zbar_n
  use REcoM_GloVar, &
       only: GloPCO2surf, GloCO2flux, Diags3D, PAR3D, export, &
       PistonVelocity, alphaCO2
  use i_arrays, &
       only: a_ice
  use g_clock, &
       only: timenew, timeold, daynew, dayold, yearnew, yearold, &
       month, day_in_month, num_day_in_month, cyearnew, cyearold, &
       fleapyear, check_fleapyr, dt, clock
  use g_events, &
       only: daily_event, monthly_event
  use g_comm_auto, &
       only: exchange_nod, exchange_elem, broadcast_nod
  use g_rotate_grid, &
       only: r2g                 ! Transform from the mesh (rotated) coordinates to geographical coordinates  
  use o_param, &
       only: r_earth, rad, num_tracers, pi
  use recom_config, &
       only: tiny_chl, tiny, chl2N_max, chl2N_max_d, NCmax, &      
       NCmax_d, SiCmax, Redfield, SecondsPerDay
  use g_config, &
       only: runid, step_per_day, ResultPath
  use g_sbf, &
       only: atmdata, i_xwind, i_ywind, i_humi, &
       i_qsr, i_qlw, i_tair, i_prec, i_mslp, i_snow

  implicit none

!  REAL, PARAMETER :: pi=3.14159265358979
!  REAL, PARAMETER :: pi=3.141592653589793   ! Pi

  ! FESOM mesh:
  type(t_mesh), pointer, save :: mesh_fesom
  integer, parameter :: nlmax = 46            ! CORE2 mesh: deepest wet cells at mesh_fesom%nl-2
  real, allocatable :: topography3D(:,:)      ! topography: 1 for wet nodes and 0 for dry nodes (array shape as in model)
  real, allocatable :: topography_p(:)        ! """                                             (array shape as state_p)
  real, allocatable :: topography3D_g(:,:)    ! """                                             (array shape as in model globally)
  real :: area_surf_glob(nlmax)               ! ocean area and standard volume to calculate area-/volume weighted means
  real :: inv_area_surf_glob(nlmax)
  real :: volo_full_glob, inv_volo_full_glob
  real, allocatable :: cellvol(:,:)           ! standard volume of cells, NOT considering time-varying ALE layerwidth

end module fesom_pdaf
