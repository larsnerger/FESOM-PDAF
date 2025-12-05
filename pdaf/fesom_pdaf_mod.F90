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
  USE g_parsup, ONLY: &
       myDim_nod2D, &          !Process-local number of 2D nodes  
       myDim_nod2D, &          ! Process-local number of vertices
       myDim_elem2D, &         ! Process-local number of elements 
       myList_nod2D

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
