!> Building the state vector
!!
!! This module provides variables & routines for
!! defining the state vector.
!!
!! The module contains three routines
!! - **init_id** - initialize the array `id`
!! - **init_sfields** - initialize the array `sfields`
!! - **setup_statevector** - generic routine controlling the initialization
!!
!! The declarations of **id** and **sfields** as well as the
!! routines **init_id** and **init_sfields** might need to be
!! adapted to a particular modeling case. However, for most
!! parts also the configruation using the namelist is possible.
!!
!! __Revision history__
!! ~2022 Frauke - initial functionality distributed over different routines
!! 2025-12 - Lars Nerger - restructuring code introducing module statevector_pdaf
module statevector_pdaf

  implicit none
  save

  !---- `field_ids` and `state_field` can be adapted for a DA case -----

  ! Declare Fortran type holding the indices of model fields in the state vector
  ! This can be extended to any number of fields - it serves to give each field a name
  type field_ids
     integer :: ssh        ! physics
     integer :: u 
     integer :: v 
     integer :: w 
     integer :: temp 
     integer :: salt
     integer :: a_ice
     integer :: MLD1
     integer :: MLD2
     integer :: PhyChl     ! chlorophyll
     integer :: DiaChl
     integer :: DIC        ! dissolved tracers
     integer :: DOC
     integer :: Alk
     integer :: DIN
     integer :: DON
     integer :: O2
     integer :: pCO2s      ! surface carbon diagnostics
     integer :: CO2f
     integer :: alphaCO2
     integer :: PistonVel
     integer :: PhyN       ! small phyto
     integer :: PhyC
     integer :: PhyCalc
     integer :: DiaN       ! diatoms
     integer :: DiaC
     integer :: DiaSi
     integer :: Zo1C       ! zooplankton
     integer :: Zo1N
     integer :: Zo2C
     integer :: Zo2N
     integer :: DetC       ! detritus
     integer :: DetCalc
     integer :: DetSi
     integer :: DetN
     integer :: Det2C
     integer :: Det2Calc
     integer :: Det2Si
     integer :: Det2N
     integer :: PAR        ! diags
     integer :: NPPn
     integer :: NPPd
     integer :: export
     integer :: sigma
!     INTEGER :: TChl   ! Total chlorophyll = PhyChl + DiaChl
!     INTEGER :: TDN    ! Total dissolved N = DIN + DON
!     INTEGER :: TOC    ! Total organic carbon: PhyC + DiaC + DetC + DOC + HetC
  end type field_ids

  ! Declare Fortran type holding the definitions for model fields
  type state_field
     integer :: ndims = 0                  !< Number of field dimensions (2 or 3)
     integer :: dim = 0                    !< Dimension of the field
     integer :: off = 0                    !< Offset of field in state vector
     logical :: nz1 = .true.               !< Vertical coordinates (on levels / on layers)
     character(len=10) :: variable = ''    !< Name of field
     character(len=50) :: long_name = ''   !< Long name of field
     character(len=20) :: units = ''       !< Unit of variable
     integer :: varid(9)                   !< To write to netCDF file
     logical :: updated = .true.           !< Whether variable is updated through assimilation
     logical :: output(8,3) = .false.      !< How frequently output is written
     logical :: bgc = .false.              !< Whether variable is biogeochemistry (or physics)
     integer :: trnumfesom = -1            !< Tracer index in FESOM-REcoM
     integer :: tridfesom = -1             !< Tracer ID in FESOM-REcoM
     integer :: id_tr = -1                 !< Field index in list of 3D model tracer fields
  end type state_field


  ! Declare Fortran type holding the definitions for local model fields
  ! This is separate from state_field to support OpenMP
  type state_field_l
     integer :: dim = 0                    !< Dimension of the field
     integer :: off = 0                    !< Offset of field in state vector
  end type state_field_l

  integer :: phymin, phymax   ! First and last physics field in state vector
  integer :: bgcmin, bgcmax   ! First and last biogeochemistry field in state vector


  !---- The next variables usually do not need editing -----

  ! Type variable holding field IDs in state vector
  type(field_ids) :: id

  ! Type variable holding the definitions of model fields
  type(state_field), allocatable :: sfields(:)

  ! Type variable holding the definitions of local model fields
  ! This is separate from sfields to support OpenMP
  type(state_field_l), allocatable :: sfields_l(:)

!$OMP THREADPRIVATE(sfields_l)

  ! Variables to handle multiple fields in the state vector
  integer :: nfields           !< number of fields in state vector

contains

!> This routine initializes the array `id`
!!
  subroutine init_id(nfields)

    implicit none

! *** Arguments ***
    integer, intent(out) :: nfields

! Total number of fields
    nfields = 44

! physics part of state vector, specify start and end:
    phymin = 1
    phymax = 10
  
! BGC part of state vector, specify start and end:
    bgcmin = 11
    bgcmax = nfields

! Set field IDs
    id% ssh    =  1 ! sea surface height
    id% u      =  2 ! zonal velocity
    id% v      =  3 ! meridional velocity
    id% w      =  4 ! vertical velocity
    id% temp   =  5 ! temperature
    id% salt   =  6 ! salinity
    id% a_ice  =  7 ! sea-ice concentration
    id% MLD1   =  8 ! boundary layer depth (criterion after Large et al., 1997)
    id% MLD2   =  9 ! mixed layer depth (density treshold)
    id% sigma  = 10

    id% PhyChl = 11 ! chlorophyll-a small phytoplankton
    id% DiaChl = 12 ! chlorophyll-a diatoms

    id% DIC    = 13 ! dissolved tracers
    id% DOC    = 14
    id% Alk    = 15
    id% DIN    = 16
    id% DON    = 17
    id% O2     = 18

    id% pCO2s     = 19 ! surface carbon diags
    id% CO2f      = 20
    id% alphaCO2  = 21
    id% PistonVel = 22

    id% PhyN   = 23 ! small phyto
    id% PhyC   = 24
    id% PhyCalc= 25

    id% DiaN   = 26 ! diatoms
    id% DiaC   = 27
    id% DiaSi  = 28

    id% Zo1N   = 29 ! zooplankton
    id% Zo2C   = 30
    id% Zo2N   = 31
    id% Zo1C   = 32

    id% DetC      = 33 ! detritus
    id% DetCalc   = 34
    id% DetSi     = 35
    id% DetN      = 36
    id% Det2C     = 37
    id% Det2Calc  = 38
    id% Det2Si    = 39
    id% Det2N     = 40

    id% PAR    = 41 ! diags
    id% NPPn   = 42
    id% NPPd   = 43
    id% export = 44

  end subroutine init_id
! ===================================================================================

!> This initializes the array sfields
!!
!! This routine initializes the sfields array with specifications
!! of the fields in the state vector.
!!
  subroutine init_sfields()

    use fesom_pdaf, &
         only: myDim_nod2D, nlmax

    implicit none

! *** Local variables ***
    integer :: i           ! Counter


! *** Allocate ***

    allocate(sfields(nfields))


! ****************
! *** Physics ****
! ****************

! SSH
    sfields(id% ssh) % ndims = 1
    sfields(id% ssh) % variable = 'SSH'
    sfields(id% ssh) % long_name = 'Sea surface height'
    sfields(id% ssh) % units = 'm'
    sfields(id% ssh) % updated = .false.
    sfields(id% ssh) % bgc = .false.

! u
    sfields(id% u) % ndims = 2
    sfields(id% u) % nz1 = .true.
    sfields(id% u) % variable = 'u'
    sfields(id% u) % long_name = 'Zonal velocity (interpolated on nodes)'
    sfields(id% u) % units = 'm/s'
    sfields(id% u) % updated = .false.
    sfields(id% u) % bgc = .false.

! v
    sfields(id% v) % ndims = 2
    sfields(id% v) % nz1 = .true.
    sfields(id% v) % variable = 'v'
    sfields(id% v) % long_name = 'Meridional velocity (interpolated on nodes)'
    sfields(id% v) % units = 'm/s'
    sfields(id% v) % updated = .false.
    sfields(id% v) % bgc = .false.

! w
    sfields(id% w) % ndims = 2
    sfields(id% w) % nz1 = .false.
    sfields(id% w) % variable = 'w'
    sfields(id% w) % long_name = 'Vertical velocity'
    sfields(id% w) % units = 'm/s'
    sfields(id% w) % updated = .false.
    sfields(id% w) % bgc = .false.

! temp
    sfields(id% temp) % ndims = 2
    sfields(id% temp) % nz1 = .true.
    sfields(id% temp) % variable = 'T'
    sfields(id% temp) % long_name = 'Temperature'
    sfields(id% temp) % units = 'degC'
    sfields(id% temp) % updated = .false.
    sfields(id% temp) % bgc = .false.

! salt
    sfields(id% salt) % ndims = 2
    sfields(id% salt) % nz1 = .true.
    sfields(id% salt) % variable = 'S'
    sfields(id% salt) % long_name = 'Salinity'
    sfields(id% salt) % units = 'psu'
    sfields(id% salt) % updated = .false.
    sfields(id% salt) % bgc = .false.


! **********************
! *** Physics Diags ****
! **********************

! ice
    sfields(id% a_ice) % ndims = 1
    sfields(id% a_ice) % variable = 'ice'
    sfields(id% a_ice) % long_name = 'Sea-ice concentration'
    sfields(id% a_ice) % units = '1'
    sfields(id% a_ice) % updated = .false.
    sfields(id% a_ice) % bgc = .false.

! MLD1
    sfields(id% MLD1) % ndims = 1
    sfields(id% MLD1) % variable = 'MLD1'
    sfields(id% MLD1) % long_name = 'Boundary layer depth (Large et al. 1997)'
    sfields(id% MLD1) % units = 'm'
    sfields(id% MLD1) % updated = .false.
    sfields(id% MLD1) % bgc = .false.

! MLD2
    sfields(id% MLD2) % ndims = 1
    sfields(id% MLD2) % variable = 'MLD2'
    sfields(id% MLD2) % long_name = 'Mixed layer depth (sigma 0.03; Boyer-Montegut et al. 2004)'
    sfields(id% MLD2) % units = 'm'
    sfields(id% MLD2) % updated = .false.
    sfields(id% MLD2) % bgc = .false.

! **********************
! *** Chlorophyll   ****
! **********************

! chlorophyll small phytoplankton
    sfields(id% PhyChl) % ndims = 2
    sfields(id% PhyChl) % nz1 = .true.
    sfields(id% PhyChl) % variable = 'PhyChl'
    sfields(id% PhyChl) % long_name = 'Chlorophyll-a small phytoplankton'
    sfields(id% PhyChl) % units = 'mg chl m-3'
    sfields(id% PhyChl) % updated = .false.
    sfields(id% PhyChl) % bgc = .true.

! chlorophyll diatoms
    sfields(id% DiaChl) % ndims = 2
    sfields(id% DiaChl) % nz1 = .true.
    sfields(id% DiaChl) % variable = 'DiaChl'
    sfields(id% DiaChl) % long_name = 'Chlorophyll-a diatoms'
    sfields(id% DiaChl) % units = 'mg chl m-3'
    sfields(id% DiaChl) % updated = .false.
    sfields(id% DiaChl) % bgc = .true.

! **************************
! *** Dissolved tracers ****
! **************************

! DIC
    sfields(id% DIC) % ndims = 2
    sfields(id% DIC) % nz1 = .true.
    sfields(id% DIC) % variable = 'DIC'
    sfields(id% DIC) % long_name = 'Dissolved inorganic carbon'
    sfields(id% DIC) % units = 'mmol C m-3'
    sfields(id% DIC) % updated = .false.
    sfields(id% DIC) % bgc = .true.

! DOC
    sfields(id% DOC) % ndims = 2
    sfields(id% DOC) % nz1 = .true.
    sfields(id% DOC) % variable = 'DOC'
    sfields(id% DOC) % long_name = 'Dissolved organic carbon'
    sfields(id% DOC) % units = 'mmol C m-3'
    sfields(id% DOC) % updated = .false.
    sfields(id% DOC) % bgc = .true.

! Alkalinity
    sfields(id% Alk) % ndims = 2
    sfields(id% Alk) % nz1 = .true.
    sfields(id% Alk) % variable = 'Alk'
    sfields(id% Alk) % long_name = 'Alkalinity'
    sfields(id% Alk) % units = 'mmol m-3'
    sfields(id% Alk) % updated = .false.
    sfields(id% Alk) % bgc = .true.

! DIN
    sfields(id% DIN) % ndims = 2
    sfields(id% DIN) % nz1 = .true.
    sfields(id% DIN) % variable = 'DIN'
    sfields(id% DIN) % long_name = 'Dissolved inorganic nitrogen'
    sfields(id% DIN) % units = 'mmol m-3'
    sfields(id% DIN) % updated = .false.
    sfields(id% DIN) % bgc = .true.

! DON
    sfields(id% DON) % ndims = 2
    sfields(id% DON) % nz1 = .true.
    sfields(id% DON) % variable = 'DON'
    sfields(id% DON) % long_name = 'Dissolved organic nitrogen'
    sfields(id% DON) % units = 'mmol m-3'
    sfields(id% DON) % updated = .false.
    sfields(id% DON) % bgc = .true.

! Oxygen
    sfields(id% O2) % ndims = 2
    sfields(id% O2) % nz1 = .true.
    sfields(id% O2) % variable = 'O2'
    sfields(id% O2) % long_name = 'Oxygen'
    sfields(id% O2) % units = 'mmol m-3'
    sfields(id% O2) % updated = .false.
    sfields(id% O2) % bgc = .true.

! *****************************
! *** Surface Carbon Diags ****
! *****************************

! pCO2
    sfields(id% pCO2s) % ndims = 1
    sfields(id% pCO2s) % variable = 'pCO2s'
    sfields(id% pCO2s) % long_name = 'Partial pressure CO2 surface ocean'
    sfields(id% pCO2s) % units = 'micro atm'
    sfields(id% pCO2s) % updated = .false.
    sfields(id% pCO2s) % bgc = .true.

! CO2f
    sfields(id% CO2f) % ndims = 1
    sfields(id% CO2f) % variable = 'CO2f'
    sfields(id% CO2f) % long_name = 'CO2 flux from atmosphere into ocean'
    sfields(id% CO2f) % units = 'mmol C m-2 d-1'
    sfields(id% CO2f) % updated = .false.
    sfields(id% CO2f) % bgc = .true.

! Solubility of CO2
    sfields(id% alphaCO2) % ndims = 1
    sfields(id% alphaCO2) % variable = 'alphaCO2'
    sfields(id% alphaCO2) % long_name = 'solubility of surface CO2'
    sfields(id% alphaCO2) % units = 'mol kg-1 atm-1'
    sfields(id% alphaCO2) % updated = .false.
    sfields(id% alphaCO2) % bgc = .true.

! Piston velocity
    sfields(id% PistonVel) % ndims = 1
    sfields(id% PistonVel) % variable = 'Kw660'
    sfields(id% PistonVel) % long_name = 'air-sea piston velocity'
    sfields(id% PistonVel) % units = 'm/s'
    sfields(id% PistonVel) % updated = .false.
    sfields(id% PistonVel) % bgc = .true.

! *****************************
! *** Small Phyto          ****
! *****************************

! PhyN
    sfields(id% PhyN) % ndims = 2
    sfields(id% PhyN) % variable = 'PhyN'
    sfields(id% PhyN) % long_name = 'intracell nitrogen small phytoplankton'
    sfields(id% PhyN) % units = 'mmol m-3'
    sfields(id% PhyN) % updated = .false.
    sfields(id% PhyN) % bgc = .true.

! PhyC
    sfields(id% PhyC) % ndims = 2
    sfields(id% PhyC) % variable = 'PhyC'
    sfields(id% PhyC) % long_name = 'intracell carbon small phytoplankton'
    sfields(id% PhyC) % units = 'mmol C m-3'
    sfields(id% PhyC) % updated = .false.
    sfields(id% PhyC) % bgc = .true.

! PhyCalc
    sfields(id% PhyCalc) % ndims = 2
    sfields(id% PhyCalc) % variable = 'PhyCalc'
    sfields(id% PhyCalc) % long_name = 'calcium carbonate small phytoplankton'
    sfields(id% PhyCalc) % units = 'mmol m-3'
    sfields(id% PhyCalc) % updated = .false.
    sfields(id% PhyCalc) % bgc = .true.

! *****************************
! *** diatoms              ****
! *****************************

! DiaN
    sfields(id% DiaN) % ndims = 2
    sfields(id% DiaN) % variable = 'DiaN'
    sfields(id% DiaN) % long_name = 'intracell nitrogen diatoms'
    sfields(id% DiaN) % units = 'mmol m-3'
    sfields(id% DiaN) % updated = .false.
    sfields(id% DiaN) % bgc = .true.

! DiaC
    sfields(id% DiaC) % ndims = 2
    sfields(id% DiaC) % variable = 'DiaC'
    sfields(id% DiaC) % long_name = 'intracell carbon diatom'
    sfields(id% DiaC) % units = 'mmol C m-3'
    sfields(id% DiaC) % updated = .false.
    sfields(id% DiaC) % bgc = .true.

! DiaSi
    sfields(id% DiaSi) % ndims = 2
    sfields(id% DiaSi) % variable = 'DiaSi'
    sfields(id% DiaSi) % long_name = 'intracell Si diatom'
    sfields(id% DiaSi) % units = 'mmol m-3'
    sfields(id% DiaSi) % updated = .false.
    sfields(id% DiaSi) % bgc = .true.

! *****************************
! *** zooplankton          ****
! *****************************

! Zo1C
    sfields(id% Zo1C) % ndims = 2
    sfields(id% Zo1C) % variable = 'Zo1C'
    sfields(id% Zo1C) % long_name = 'carbon in small zooplankton'
    sfields(id% Zo1C) % units = 'mmol C m-3'
    sfields(id% Zo1C) % updated = .false.
    sfields(id% Zo1C) % bgc = .true.

! Zo1N
    sfields(id% Zo1N) % ndims = 2
    sfields(id% Zo1N) % variable = 'Zo1N'
    sfields(id% Zo1N) % long_name = 'nitrogen in small zooplankton'
    sfields(id% Zo1N) % units = 'mmol C m-3'
    sfields(id% Zo1N) % updated = .false.
    sfields(id% Zo1N) % bgc = .true.

! Zo2C
    sfields(id% Zo2C) % ndims = 2
    sfields(id% Zo2C) % variable = 'Zo2C'
    sfields(id% Zo2C) % long_name = 'carbon in macrozooplankton'
    sfields(id% Zo2C) % units = 'mmol C m-3'
    sfields(id% Zo2C) % updated = .false.
    sfields(id% Zo2C) % bgc = .true.

! Zo2N
    sfields(id% Zo2N) % ndims = 2
    sfields(id% Zo2N) % variable = 'Zo2N'
    sfields(id% Zo2N) % long_name = 'nitrogen in macrozooplankton'
    sfields(id% Zo2N) % units = 'mmol C m-3'
    sfields(id% Zo2N) % updated = .false.
    sfields(id% Zo2N) % bgc = .true.

! *****************************
! *** detritus             ****
! *****************************

! DetC
    sfields(id% DetC) % ndims = 2
    sfields(id% DetC) % variable = 'DetC'
    sfields(id% DetC) % long_name = 'carbon in small detritus'
    sfields(id% DetC) % units = 'mmol C m-3'
    sfields(id% DetC) % updated = .false.
    sfields(id% DetC) % bgc = .true.

! DetCalc
    sfields(id% DetCalc) % ndims = 2
    sfields(id% DetCalc) % variable = 'DetCalc'
    sfields(id% DetCalc) % long_name = 'calcite in small detritus'
    sfields(id% DetCalc) % units = 'mmol C m-3'
    sfields(id% DetCalc) % updated = .false.
    sfields(id% DetCalc) % bgc = .true.

! DetN
    sfields(id% DetN) % ndims = 2
    sfields(id% DetN) % variable = 'DetN'
    sfields(id% DetN) % long_name = 'nitrogen in small detritus'
    sfields(id% DetN) % units = 'mmol C m-3'
    sfields(id% DetN) % updated = .false.
    sfields(id% DetN) % bgc = .true.

! DetSi
    sfields(id% DetSi) % ndims = 2
    sfields(id% DetSi) % variable = 'DetSi'
    sfields(id% DetSi) % long_name = 'silicate in small detritus'
    sfields(id% DetSi) % units = 'mmol C m-3'
    sfields(id% DetSi) % updated = .false.
    sfields(id% DetSi) % bgc = .true.

! Det2 C
    sfields(id% Det2C) % ndims = 2
    sfields(id% Det2C) % variable = 'Det2C'
    sfields(id% Det2C) % long_name = 'carbon in large detritus'
    sfields(id% Det2C) % units = 'mmol C m-3'
    sfields(id% Det2C) % updated = .false.
    sfields(id% Det2C) % bgc = .true.

! Det2 Calc
    sfields(id% Det2Calc) % ndims = 2
    sfields(id% Det2Calc) % variable = 'Det2Calc'
    sfields(id% Det2Calc) % long_name = 'calcite in large detritus'
    sfields(id% Det2Calc) % units = 'mmol C m-3'
    sfields(id% Det2Calc) % updated = .false.
    sfields(id% Det2Calc) % bgc = .true.

! Det2 N
    sfields(id% Det2N) % ndims = 2
    sfields(id% Det2N) % variable = 'Det2N'
    sfields(id% Det2N) % long_name = 'nitrogen in large detritus'
    sfields(id% Det2N) % units = 'mmol C m-3'
    sfields(id% Det2N) % updated = .false.
    sfields(id% Det2N) % bgc = .true.

! Det2 Si
    sfields(id% Det2Si) % ndims = 2
    sfields(id% Det2Si) % variable = 'Det2Si'
    sfields(id% Det2Si) % long_name = 'silicate in large detritus'
    sfields(id% Det2Si) % units = 'mmol C m-3'
    sfields(id% Det2Si) % updated = .false.
    sfields(id% Det2Si) % bgc = .true.

! *****************************
! *** diagnostics          ****
! *****************************

! PAR
    sfields(id% PAR) % ndims = 2
    sfields(id% PAR) % variable = 'PAR'
    sfields(id% PAR) % long_name = 'photosynthetically active radiation'
    sfields(id% PAR) % units = 'W m-2'
    sfields(id% PAR) % updated = .false.
    sfields(id% PAR) % bgc = .true.

! NPPn
    sfields(id% NPPn) % ndims = 2
    sfields(id% NPPn) % variable = 'NPPn'
    sfields(id% NPPn) % long_name = 'mean net primary production small phytoplankton'
    sfields(id% NPPn) % units = 'mmol C m-2 d-1'
    sfields(id% NPPn) % updated = .false.
    sfields(id% NPPn) % bgc = .true.

! NPPd
    sfields(id% NPPd) % ndims = 2
    sfields(id% NPPd) % variable = 'NPPd'
    sfields(id% NPPd) % long_name = 'mean net primary production diatoms'
    sfields(id% NPPd) % units = 'mmol C m-2 d-1'
    sfields(id% NPPd) % updated = .false.
    sfields(id% NPPd) % bgc = .true.

! Export production
    sfields(id% export) % ndims = 1
    sfields(id% export) % variable = 'export'
    sfields(id% export) % long_name = 'export through particle sinking at 190m'
    sfields(id% export) % units = 'mmol m-2 day-1'
    sfields(id% export) % updated = .false.
    sfields(id% export) % bgc = .true.

! Potential density
    sfields(id% sigma) % ndims = 2
    sfields(id% sigma) % variable = 'sigma'
    sfields(id% sigma) % long_name = 'potential density'
    sfields(id% sigma) % units = 'kg liter-1'
    sfields(id% sigma) % updated = .false.
    sfields(id% sigma) % bgc = .false.

    !~ ! TChl
    !~ sfields(id% TChl) % ndims = 2
    !~ sfields(id% TChl) % variable = 'TChl'
    !~ sfields(id% TChl) % long_name = 'Total chlorophyll (PhyChl+DiaChl)'
    !~ sfields(id% TChl) % units = 'mg chl m-3'
    !~ sfields(id% TChl) % updated = .false.
    !~ sfields(id% TChl) % bgc = .true.

    !~ ! TDN
    !~ sfields(id% TDN) % ndims = 2
    !~ sfields(id% TDN) % variable = 'TDN'
    !~ sfields(id% TDN) % long_name = 'Total dissolved nitrogen (DIN+DON)'
    !~ sfields(id% TDN) % units = 'mmol m-3'
    !~ sfields(id% TDN) % updated = .false.
    !~ sfields(id% TDN) % bgc = .true.

    !~ ! TOC
    !~ sfields(id% TOC) % ndims = 2
    !~ sfields(id% TOC) % variable = 'TOC'
    !~ sfields(id% TOC) % long_name = 'Total Organic Carbon (PhyC+DiaC+DetC+DOC+HetC)'
    !~ sfields(id% TOC) % units = 'mmol C m-3'
    !~ sfields(id% TOC) % updated = .false.
    !~ sfields(id% TOC) % bgc = .true.


! **************************************
! ***   Set dimensions and offsets   ***
! **************************************

    ! *** Dimensions ***
    do i = 1, nfields
       if (sfields(i)%ndims == 1) then
          sfields(i)%dim = myDim_nod2D
       else if (sfields(i)%ndims == 2) then
          sfields(i)%dim = myDim_nod2D*nlmax
       else
          write (*, '(a,i2,a)') 'FESOM-PDAF: cannot handle', sfields(i)%ndims, ' number of dimensions.'
       end if
    end do

! *** Specify offset of fields in pe-local state vector ***

!    . . . A . . . . . B . . . . . C . . . . . D
!         . .         / .         / .         . .
!        .   .   2   /   .   3   /   .   5   .   .
! .     .     .     /     .     /     .     .     .
!  .   .   1   .   /   3   .   /   4   .   .       
!   . .         . /         . /         . .        
!    A . . . . . B . . . . . C . . . . . D . . . . 
!
!  A:    Internal nodes of left PE
!  B:    Internal nodes of left PE, simultanesously external nodes of right PE
!  C:    External nodes of left PE, simultanesously internal nodes of right PE
!  D:    Internal nodes of right PE
!  1:    Internal element of left PE, simultanesously wide-halo element of right PE (shares node B with right PE)
!  2:    Internal element of left PE, simultanesously small-halo element of right PE (shares edge BB with right PE)
!  3:    Internal (shared) elements of both PEs
!  4:    Small-halo element of left PE (shares edge CC with left PE), simultanesously internal element of right PE
!  5:    Wide-halo element of left PE (shares node C with left PE), simultanesously internal element of right PE
!
!  myDim_nod2D:      Number of internal nodes (A+B)
!  eDim_nod2D:       Number of external nodes (C)
!
!  myDim_elem2D:     Number of internal elements (1+2+3)
!  eDim_elem2D:      Number of small-halo elements (4)
!  xDim_elem2D:      Number of wide-halo elements (5)
!
!  mesh_fesom%nl:    Maximum number of fesom levels (1 is air-sea interface)
!  mesh_fesom%nl-1:  Maximum number of fesom layers (1 is surface layer, 0-5m)

!  CORE2 mesh: deepest wet cells at mesh_fesom%nl-2 or nlmax=46

    ! *** Define offsets in state vector ***
    sfields(1)%off = 0
    do i = 2, nfields
       sfields(i)%off = sfields(i-1)%off + sfields(i-1)%dim
    end do


! ************************************************
! ***   Tracer index and ID from FESOM-REcoM   ***
! ************************************************

    sfields(id%temp    )%trnumfesom = 1  ! temperature
    sfields(id%salt    )%trnumfesom = 2  ! salinity 
    sfields(id%DIN     )%trnumfesom = 3  ! DIN
    sfields(id%DIC     )%trnumfesom = 4  ! DIC
    sfields(id%Alk     )%trnumfesom = 5  ! Alk
    sfields(id%PhyN    )%trnumfesom = 6  ! PhyN 
    sfields(id%PhyC    )%trnumfesom = 7  ! PhyC   
    sfields(id%PhyChl  )%trnumfesom = 8  ! PhyChl
    sfields(id%DetN    )%trnumfesom = 9  ! DetN
    sfields(id%DetC    )%trnumfesom = 10 ! DetC
    sfields(id%Zo1N    )%trnumfesom = 11 ! HetN
    sfields(id%Zo1C    )%trnumfesom = 12 ! HetC
    sfields(id%DON     )%trnumfesom = 13 ! DON
    sfields(id%DOC     )%trnumfesom = 14 ! DOC
    sfields(id%DiaN    )%trnumfesom = 15 ! DiaN
    sfields(id%DiaC    )%trnumfesom = 16 ! DiaC
    sfields(id%DiaChl  )%trnumfesom = 17 ! DiaChl
    sfields(id%DiaSi   )%trnumfesom = 18 ! DiaSi
    sfields(id%DetSi   )%trnumfesom = 19 ! DetSi 
    !sfields(id%        )%trnumfesom = 20 ! DSi     
    !sfields(id%        )%trnumfesom = 21 ! Fe
    sfields(id%PhyCalc )%trnumfesom = 22 ! PhyCalc
    sfields(id%DetCalc )%trnumfesom = 23 ! DetCalc
    sfields(id%O2      )%trnumfesom = 24 ! Oxy
    sfields(id%Zo2N    )%trnumfesom = 25 ! Zoo2N
    sfields(id%Zo2C    )%trnumfesom = 26 ! Zoo2C
    sfields(id%Det2N   )%trnumfesom = 27 ! DetZ2N                              
    sfields(id%Det2C   )%trnumfesom = 28 ! DetZ2C                                    
    sfields(id%Det2Si  )%trnumfesom = 29 ! DetZ2Si                            
    sfields(id%Det2Calc)%trnumfesom = 30 ! DetZ2Calc

    sfields(id%temp    )%tridfesom =    0 ! temperature
    sfields(id%salt    )%tridfesom =    1 ! salinity
    sfields(id%DIN     )%tridfesom = 1001 ! DIN
    sfields(id%DIC     )%tridfesom = 1002 ! DIC
    sfields(id%Alk     )%tridfesom = 1003 ! Alk
    sfields(id%PhyN    )%tridfesom = 1004 ! PhyN 
    sfields(id%PhyC    )%tridfesom = 1005 ! PhyC   
    sfields(id%PhyChl  )%tridfesom = 1006 ! PhyChl
    sfields(id%DetN    )%tridfesom = 1007 ! DetN
    sfields(id%DetC    )%tridfesom = 1008 ! DetC
    sfields(id%Zo1N    )%tridfesom = 1009 ! HetN
    sfields(id%Zo1C    )%tridfesom = 1010 ! HetC
    sfields(id%DON     )%tridfesom = 1011 ! DON
    sfields(id%DOC     )%tridfesom = 1012 ! DOC
    sfields(id%DiaN    )%tridfesom = 1013 ! DiaN
    sfields(id%DiaC    )%tridfesom = 1014 ! DiaC
    sfields(id%DiaChl  )%tridfesom = 1015 ! DiaChl
    sfields(id%DiaSi   )%tridfesom = 1016 ! DiaSi
    sfields(id%DetSi   )%tridfesom = 1017 ! DetSi 
    !sfields(id%        )%tridfesom = 1018 ! DSi     
    !sfields(id%        )%tridfesom = 1019 ! Fe
    sfields(id%PhyCalc )%tridfesom = 1020 ! PhyCalc
    sfields(id%DetCalc )%tridfesom = 1021 ! DetCalc
    sfields(id%O2      )%tridfesom = 1022 ! Oxy
    sfields(id%Zo2N    )%tridfesom = 1023 ! Zoo2N
    sfields(id%Zo2C    )%tridfesom = 1024 ! Zoo2C
    sfields(id%Det2N   )%tridfesom = 1025 ! DetZ2N                              
    sfields(id%Det2C   )%tridfesom = 1026 ! DetZ2C                                    
    sfields(id%Det2Si  )%tridfesom = 1027 ! DetZ2Si                            
    sfields(id%Det2Calc)%tridfesom = 1028 ! DetZ2Calc


  end subroutine init_sfields
! ===================================================================================

!> Count 2D and 3D fields and initialize index arrays
!!
!! This functioanlity of optional and not used elsewhere
!! in the code.
!!
  subroutine set_field_types(verbose)

    use parallel_pdaf_mod, &
         only: mype_world

    implicit none

! *** Arguments ***
    integer, intent(in) :: verbose     ! Verbosity level

! *** Local variables ***
    integer :: i, cnt      ! Counters
    integer, allocatable :: ids_3D(:)         ! List of 3D-field IDs
    integer, allocatable :: ids_2D(:)         ! """       2D fields """
    integer, allocatable :: ids_phy(:)        ! """       physics fields """
    integer, allocatable :: ids_bgc(:)        ! """       biogeochem. fields """
    integer, allocatable :: ids_tr3D(:)       ! List of 3D model tracer IDs
    integer :: nfields_3D                     ! Number of 3D fields in state vector
    integer :: nfields_2D                     ! """       2D fields """
    integer :: nfields_phy                    ! """       physics fields """
    integer :: nfields_bgc                    ! """       biogeochem. fields """
    integer :: nfields_tr3D                   ! Number of 3D model tracer fields in state vector


! **************************************
! ***  Indices of by type of field   ***
! **************************************
  
! ***2D/3D***

    ! count number of 3D and 2D fields (tracers and diagnostics)
    nfields_3D = 0
    nfields_2D = 0
    do i=1,nfields
       if (sfields(i)%ndims == 2) nfields_3D = nfields_3D + 1
       if (sfields(i)%ndims == 1) nfields_2D = nfields_2D + 1
    enddo
    
    ! Store indices of 3D and 2D fields
    allocate(ids_3D(nfields_3D))
    cnt = 1
    do i=1,nfields
       if (sfields(i)%ndims == 2) then
          ids_3D(cnt) = i
          cnt = cnt+1
       endif
    enddo
    if (mype_world==0 .and. verbose>2) then
       do cnt=1,nfields_3D
          write (*,'(a, 10x,3a,1x,7a)') &
               'FESOM-PDAF', '3D fields in state vector: ', sfields(ids_3D(cnt))%variable
       enddo
    endif
    
    allocate(ids_2D(nfields_2D))
    cnt = 1
    do i=1,nfields
       if (sfields(i)%ndims == 1) then
          ids_2D(cnt) = i
          cnt = cnt+1
       endif
    enddo
    if (mype_world==0 .and. verbose>2) then
       do cnt=1,nfields_2D
          write (*,'(a, 10x,3a,1x,7a)') &
               'FESOM-PDAF', '2D fields in state vector: ', sfields(ids_2D(cnt))%variable
       enddo
    endif
    
! *** model 3D tracers ***

    ! count number of 3D model tracer fields in state vector (only tracers)
    nfields_tr3D = 0
    do i=1,nfields
       if (sfields(i)%trnumfesom > 0) nfields_tr3D = nfields_tr3D + 1
    enddo
    
    ! Store indices of 3D model tracer fields in state vector
    allocate(ids_tr3D(nfields_tr3D))
    cnt = 1
    do i=1,nfields
       if (sfields(i)%trnumfesom > 0) then
          ids_tr3D(cnt)=i
          sfields(i)%id_tr = cnt
          cnt = cnt+1
       endif
    enddo
    if (mype_world==0 .and. verbose>2) then
       do cnt=1,nfields_tr3D
          write (*,'(a, 10x,3a,1x,7a)') &
               'FESOM-PDAF', '3D model tracer fields in state vector: ', sfields(ids_tr3D(cnt))%variable
       enddo
    endif
  
! *** phy/bgc fields  ***_

    ! count number of phy/bgc fields
    nfields_phy = 0
    nfields_bgc = 0
    do i=1,nfields
       if (      sfields(i)%bgc) nfields_bgc = nfields_bgc + 1
       if (.not. sfields(i)%bgc) nfields_phy = nfields_phy + 1
    enddo
    
    ! Store indices of BGC and physics fields in state vector
    allocate(ids_bgc(nfields_bgc))
    cnt = 1
    do i=1,nfields
       if (sfields(i)%bgc) then
          ids_bgc(cnt) = i
          cnt = cnt+1
       endif
    enddo
    if (mype_world==0 .and. verbose>2) then
       do cnt=1,nfields_bgc
          write (*,'(a, 10x,3a,1x,7a)') &
               'FESOM-PDAF', 'bgc fields in state vector: ', sfields(ids_bgc(cnt))%variable
       enddo
    endif
    
    allocate(ids_phy(nfields_phy))
    cnt = 1
    do i=1,nfields
       if (.not. sfields(i)%bgc) then
          ids_phy(cnt) = i
          cnt = cnt+1
       endif
    enddo
    if (mype_world==0 .and. verbose>2) then
       do cnt=1,nfields_phy
          write (*,'(a, 10x,3a,1x,7a)') &
               'FESOM-PDAF', 'phy fields in state vector: ', sfields(ids_phy(cnt))%variable
       enddo
    endif

  end subroutine set_field_types
! ===================================================================================

!> Set which fields are updated by the DA
!!
!! This routine read from the namelist which fields should be update
!! and set the updated flags in sfields.
!!
  subroutine set_updated()

    use assim_pdaf_mod, &
         only: nmlfile
  use coupled_da_mod, &                 ! Variables for coupled DA
       only: assimilatePHY, assimilateBGC, cda_phy, cda_bio
    use parallel_pdaf_mod, &
         only: mype_world


    implicit none

    integer :: i                                                ! Counter
    logical :: upd_ssh, upd_u, upd_v, upd_w, upd_temp, upd_salt, upd_ice, &  ! Physics
         upd_MLD1, upd_MLD2, &                                  ! physics diagnostics
         upd_PhyChl, upd_DiaChl, &                              ! chlorophyll
         upd_DIC, upd_DOC, upd_Alk, upd_DIN, upd_DON, upd_O2, & ! dissolved tracers
         upd_pCO2s, upd_CO2f, upd_alphaCO2, upd_PistonVel, &    ! surface carbon diags
         upd_DiaN, upd_DiaC, upd_DiaSi, &                       ! diatoms
         upd_PhyCalc, upd_PhyC, upd_PhyN, &                     ! small phyto
         upd_Zo1C, upd_Zo1N, &                                  ! het
         upd_Zo2C, upd_Zo2N, &                                  ! zoo 2
         upd_DetC, upd_DetCalc, upd_DetSi, upd_DetN     , &     ! small det
         upd_Det2C, upd_Det2N, upd_Det2Si, upd_Det2Calc , &     ! large det
         upd_export, upd_PAR, upd_NPPn, upd_NPPd, upd_sigma     ! diags


! ************************************************
! ***   Read updated variables from namelist   ***
! ************************************************

  ! The logical "updated" describes whether a variables is updated in at least one sweep

  ! In case of diagnostic variables, "updated" is False. Setting diagnostics variables to False
  ! and the others to True, is set in namelist.
  ! In case of weak coupling and only physics or BGC assimilation, "updated" is False for the
  ! other type of fields. This is done below.
  ! "updated" is used in init_dim_l_pdaf: only updated fields are included in local state
  ! "updated" is used in the output routine: option to write out only updated fields

! *** Read namelist file ***
    if (mype_world==0) write(*,*) 'Read namelist file for updated variables: ',nmlfile
  
    namelist /updated/ &
         upd_ssh, upd_u, upd_v, upd_w, upd_temp, upd_salt, upd_ice, &  ! Physics
         upd_MLD1, upd_MLD2, &                                  ! physics diagnostics
         upd_PhyChl, upd_DiaChl, &                              ! chlorophyll
         upd_DIC, upd_DOC, upd_Alk, upd_DIN, upd_DON, upd_O2, & ! dissolved tracers
         upd_pCO2s, upd_CO2f, upd_alphaCO2, upd_PistonVel, &    ! surface carbon diags
         upd_DiaN, upd_DiaC, upd_DiaSi, &                       ! diatoms
         upd_PhyCalc, upd_PhyC, upd_PhyN, &                     ! small phyto
         upd_Zo1C, upd_Zo1N, &                                  ! het
         upd_Zo2C, upd_Zo2N, &                                  ! zoo 2
         upd_DetC, upd_DetCalc, upd_DetSi, upd_DetN     , &     ! small det
         upd_Det2C, upd_Det2N, upd_Det2Si, upd_Det2Calc , &     ! large det
         upd_export, upd_PAR, upd_NPPn, upd_NPPd, upd_sigma     ! diags

    open  (20,file=nmlfile)
    read  (20,NML=updated)
    close (20)

    ! *** Set 'updated' in sfields ***

    sfields(id%ssh      )%updated = upd_ssh
    sfields(id%u        )%updated = upd_u
    sfields(id%v        )%updated = upd_v
    sfields(id%w        )%updated = upd_w
    sfields(id%temp     )%updated = upd_temp
    sfields(id%salt     )%updated = upd_salt
    sfields(id%a_ice    )%updated = upd_ice
    sfields(id%MLD1     )%updated = upd_MLD1
    sfields(id%MLD2     )%updated = upd_MLD2
    sfields(id%sigma    )%updated = upd_sigma

    sfields(id%PhyChl   )%updated = upd_PhyChl
    sfields(id%DiaChl   )%updated = upd_DiaChl

    sfields(id%DIC      )%updated = upd_DIC
    sfields(id%DOC      )%updated = upd_DOC
    sfields(id%Alk      )%updated = upd_Alk
    sfields(id%DIN      )%updated = upd_DIN
    sfields(id%DON      )%updated = upd_DON
    sfields(id%O2       )%updated = upd_O2

    sfields(id%pCO2s    )%updated = upd_pCO2s
    sfields(id%CO2f     )%updated = upd_CO2f
    sfields(id%alphaCO2 )%updated = upd_alphaCO2
    sfields(id%PistonVel)%updated = upd_PistonVel  

    sfields(id%DiaN     )%updated = upd_DiaN
    sfields(id%DiaC     )%updated = upd_DiaC
    sfields(id%DiaSi    )%updated = upd_DiaSi

    sfields(id%PhyCalc  )%updated = upd_PhyCalc
    sfields(id%PhyC     )%updated = upd_PhyC
    sfields(id%PhyN     )%updated = upd_PhyN

    sfields(id%Zo1C     )%updated = upd_Zo1C
    sfields(id%Zo1N     )%updated = upd_Zo1N
    sfields(id%Zo2C     )%updated = upd_Zo2C
    sfields(id%Zo2N     )%updated = upd_Zo2N

    sfields(id%DetC     )%updated = upd_DetC
    sfields(id%DetCalc  )%updated = upd_DetCalc
    sfields(id%DetSi    )%updated = upd_DetSi
    sfields(id%DetN     )%updated = upd_DetN

    sfields(id%Det2C     )%updated = upd_Det2C
    sfields(id%Det2Calc  )%updated = upd_Det2Calc
    sfields(id%Det2Si    )%updated = upd_Det2Si
    sfields(id%Det2N     )%updated = upd_Det2N

    sfields(id%PAR      )%updated = upd_PAR
    sfields(id%NPPn     )%updated = upd_NPPn
    sfields(id%NPPd     )%updated = upd_NPPd
    sfields(id%export   )%updated = upd_export
  

    ! *** General settings ***

    ! Physics not assimilated and coupling weak: No update to physics
    if ((.not. assimilatePHY) .and. (cda_bio=='weak')) then
       sfields(phymin: phymax)%updated = .false.
    endif

    ! BGC not assimilated and coupling weak: No update to BGC
    if ((.not. assimilateBGC) .and. (cda_phy=='weak')) then
       sfields(bgcmin: bgcmax)%updated = .false.
    endif

  end subroutine set_updated

! ===================================================================================

!> Calculate the dimension of the process-local statevector.
!!
!! This routine is generic. case-specific adaptions should only
!! by done in the routines init_id and init_sfields.
!!
  subroutine setup_statevector(dim_state, dim_state_p, screen)

    use parallel_pdaf_mod, &
         only: mype=>mype_ens, npes=>npes_ens, task_id, comm_ensemble, &
         comm_model, MPI_SUM, MPI_INTEGER, MPIerr

    implicit none

! *** Arguments ***
    integer, intent(out) :: dim_state    !< Global dimension of state vector
    integer, intent(out) :: dim_state_p  !< Local dimension of state vector
    integer, intent(in)  :: screen       !< Verbosity flag

! *** Local variables ***
    integer :: i                 ! Counters


! ***********************************
! *** Initialize the state vector ***
! ***********************************

! *** Initialize array `id` ***

    call init_id(nfields)

! *** Initialize array `sfields` ***

    call init_sfields()
    call set_updated()
    call set_field_types(screen)

! *** Set state vector dimension ***

    dim_state_p = sum(sfields(:)%dim)

! *** Write information about the state vector ***

    if (mype==0) then
       write (*,'(/a,2x,a)') 'FESOM-PDAF', '*** Setup of state vector ***'
       write (*,'(a,5x,a,i5)') 'FESOM-PDAF', '--- Number of fields in state vector:', nfields
       write (*,'(a,a4,3x,a2,2x,a8,4x,a5,6x,a3,7x,a6,4x,a6,2x,a3,1x,a6)') &
            'FESOM-PDAF','pe','ID', 'variable', 'ndims', 'dim', 'offset', 'update', 'BGC', 'tracer'
    end if

    if (mype==0 .or. (task_id==1 .and. screen>2)) then
       do i = 1, nfields
          write (*,'(a, i4, i5,3x,a10,2x,i3,2x,i10,3x,i10,4x,l,4x,l,2x,i4)') 'FESOM-PDAF', &
               mype, i, sfields(i)%variable, sfields(i)%ndims, sfields(i)%dim, sfields(i)%off, sfields(i)%updated, &
               sfields(i)%bgc, sfields(i)%trnumfesom
       end do
    end if

    if (npes==1) then
       write (*,'(a,2x,a,1x,i10)') 'FESOM-PDAF', 'Full state dimension: ',dim_state_p
    else
       if (task_id==1) then
          if (screen>2 .or. mype==0) &
               write (*,'(a,2x,a,1x,i4,2x,a,1x,i10)') &
               'FESOM-PDAF', 'PE', mype, 'PE-local full state dimension: ',dim_state_p

          call MPI_Reduce(dim_state_p, dim_state, 1, MPI_INTEGER, MPI_SUM, 0, COMM_model, MPIerr)
          if (mype==0) then
             write (*,'(a,2x,a,1x,i10)') 'FESOM-PDAF', 'Global state dimension: ',dim_state
          end if
       end if
    end if
    call MPI_Barrier(comm_ensemble, MPIerr)

  end subroutine setup_statevector

end module statevector_pdaf
