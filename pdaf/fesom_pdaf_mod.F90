!> Include FESOM-RECOM variables and initialize FESOM grid information
!!
!! This module includes variables from FESOM-RECOM. For PDAF in its
!! online coupling it is the single point which directly links
!! to FESOM-RECOM. All other routines include from this module.
!! For the offline case the variables are declared here and
!! separately initialized.
!!
!! __Revision history:__
!! * 2025-12 - Lars Nerger - Initial code for PDAF3 revision
!! * 2026-02 - Lars Nerger  - Adaption for FESOM2.6
!!
module fesom_pdaf

  use fesom_main_storage_module, &
       only: f, t_mesh, t_partit, t_dyn, t_tracer, t_ice, &
       timenew, timeold, daynew, dayold, yearnew, yearold, &
       month, day_in_month, num_day_in_month, cyearnew, cyearold, &
       fleapyear, check_fleapyr, dt, clock, step_per_day, &
       gather_nod, exchange_nod, broadcast_nod, exchange_elem, &
       MLD1, MLD2, rotated_grid, &
       runid, ResultPath, &
       r_earth, rad, pi
  use g_rotate_grid, &
       only: r2g                 ! Transform from the mesh (rotated) coordinates to geographical coordinates  
  use g_sbf, &
       only: atmdata, i_xwind, i_ywind, i_humi, &
       i_qsr, i_qlw, i_tair, i_prec, i_mslp, i_snow
  use recom_config, &
       only: tiny, tiny_chl, chl2N_max, chl2N_max_d, NCmax, &      
       NCmax_d, SiCmax, Redfield, SecondsPerDay
  use REcoM_GloVar, &
       only: GloPCO2surf, GloCO2flux, PAR3D
  !, Diags3D, export, PistonVelocity, alphaCO2

  implicit none

  ! Additional variables related to FESOM mesh
  type(t_mesh), pointer, save :: mesh_fesom
  integer, parameter :: nlmax = 46            ! CORE2 mesh: deepest wet cells at mesh_fesom%nl-2
  real, allocatable :: topography3D(:,:)      ! topography: 1 for wet, 0 for dry nodes (array shape as in model)
  real, allocatable :: topography_p(:)        ! """                                    (array shape as state_p)
  real, allocatable :: topography3D_g(:,:)    ! """                                    (array shape as in model globally)
  real, allocatable :: cellvol(:,:)           ! standard volume of cells, NOT considering time-varying ALE layerwidth
  real :: area_surf_glob(nlmax)               ! ocean area and standard volume to calculate area-/volume weighted means
  real :: inv_area_surf_glob(nlmax)
  real :: volo_full_glob, inv_volo_full_glob

  integer :: MPI_COMM_FESOM
  integer :: mydim_nod2d, edim_nod2d
  integer :: mydim_elem2d, edim_elem2d
  integer :: mydim_edge2D
  integer, pointer :: myList_nod2D(:)
  integer, pointer :: myList_edge2D(:)
  real, pointer :: zbar_n_srf(:)
  real, pointer :: zbar_n_bot(:)

  integer :: num_tracers
  type(t_partit), pointer :: partit
  type(t_mesh), pointer   :: mesh
  type(t_dyn), pointer    :: dynamics
  type(t_tracer), pointer :: tracers
  type(t_ice), pointer :: ice
  real, pointer :: hnode_new(:,:)
  real, pointer :: eta_n(:)
  real, pointer :: UV(:,:,:)
  real, pointer :: UVnode(:,:,:)
  real, pointer :: Wvel(:,:)
  real, pointer :: a_ice(:)
  real, pointer :: u_ice(:)
  real, pointer :: v_ice(:)

contains

!> Routine to set variables for use in PDAF user routines
!!
!! FESOM2.6 puts most variables into different Fortran type
!! variables and then combines those into the general 
!! type 'f'. To avoid the need to directly access 'f'
!! or different of the Fortran types, we here define
!! regular variables so that most of the user code for
!! PDAF can remain identical for FESOM2.6 and older
!! FESOM versions.
!!
  subroutine set_fesom_pdaf_vars()

    implicit none

    partit   => f%partit
    mesh     => f%mesh
    dynamics => f%dynamics
    tracers  => f%tracers
    ice      => f%ice

    MPI_COMM_FESOM = partit%MPI_COMM_FESOM

    mydim_nod2d  = partit%mydim_nod2d
    edim_nod2d   = partit%edim_nod2d
    mydim_elem2d = partit%mydim_elem2d
    edim_elem2d  = partit%edim_elem2d
    mydim_edge2d = partit%mydim_edge2d

    myList_nod2D  => partit%myList_nod2D
    myList_edge2D => partit%myList_edge2D
    hnode_new(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D) => mesh%hnode_new(:,:)

    zbar_n_srf => mesh%zbar_n_srf
    zbar_n_bot => mesh%zbar_n_bot

    num_tracers = tracers%num_tracers

    ! Model fields
    eta_n  => dynamics%eta_n(:)
    UV     => dynamics%uv(:,:,:)
    UVnode => dynamics%uvnode(:,:,:)
    Wvel   => dynamics%w(:,:)
    a_ice  => ice%data(1)%values(:)
    u_ice  => ice%uice(:)
    v_ice  => ice%vice(:)

  end subroutine set_fesom_pdaf_vars

end module fesom_pdaf
