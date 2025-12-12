!! This file provides interface routines between the call-back routines
!! of PDAF and the observation-specific routines in PDAF-OMI. This structure
!! collects all calls to observation-specific routines in this single file
!! to make it easier to find the routines that need to be adapted.
!!
!! The routines here are mainly pure pass-through routines. Thus they
!! simply call one of the routines from PDAF-OMI. Partly some addtional
!! variable is required, e.g. to specify the offset of an observation
!! in the observation vector containing all observation types. These
!! cases are described in the routines.
!!
!! **Adding an observation type:**
!!   When adding an observation type, one has to add one module
!!   obs_TYPE_pdafomi (based on the template obs_TYPE_pdafomi_TEMPLATE.F90).
!!   In addition one has to add a call to the different routines include
!!   in this file. It is recommended to keep the order of the calls
!!   consistent over all files. 
!! 
!! __Revision history:__
!! * 2019-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!

!-------------------------------------------------------------------------------

!> Call-back routine for init_dim_obs
!!
!! This routine calls the observation-specific
!! routines init_dim_obs_TYPE.
!!
subroutine init_dim_obs_pdafomi(step, dim_obs)

  ! Include functions for different observations
  use assim_pdaf_mod, &
       only:  proffiles_o, start_year_o, end_year_o
  use parallel_pdaf_mod, &
       only: abort_parallel
       
  use obs_sst_pdafomi, &
       only: assim_o_sst, init_dim_obs_sst
  use obs_sss_smos_pdafomi, &
       only: assim_o_sss, init_dim_obs_sss
  use obs_sss_cci_pdafomi, &
       only: assim_o_sss_cci, init_dim_obs_sss_cci
  use obs_ssh_cmems_pdafomi, &
       only: assim_o_ssh, init_dim_obs_ssh
  use obs_TSprof_EN4_pdafomi, &
       only: assim_o_en4_t, assim_o_en4_s, init_dim_obs_prof
       
  use obs_chl_cci_pdafomi, &
       only: assim_o_chl_cci, init_dim_obs_chl_cci
  use obs_DIC_glodap_pdafomi, &
       only: assim_o_DIC_glodap, init_dim_obs_DIC_glodap
  use obs_Alk_glodap_pdafomi, &
       only: assim_o_Alk_glodap, init_dim_obs_Alk_glodap
  use obs_pco2_SOCAT_pdafomi, &
       only: assim_o_pCO2_SOCAT, init_dim_obs_pCO2_SOCAT
  use obs_O2_comf_pdafomi, &
       only: assim_o_O2_comf, init_dim_obs_O2_comf
  use obs_N_comf_pdafomi, &
       only: assim_o_N_comf, init_dim_obs_N_comf
  use obs_O2_argo_pdafomi, &
       only: assim_o_O2_argo, init_dim_obs_O2_argo
  use obs_N_argo_pdafomi, &
       only: assim_o_N_argo, init_dim_obs_N_argo
  use obs_O2_merged_pdafomi, &
       only: assim_o_O2_merged, init_dim_obs_O2_merged
  use obs_n_merged_pdafomi, &
       only: assim_o_n_merged, init_dim_obs_n_merged

  implicit none

! *** Arguments ***
  integer, intent(in)  :: step      !< Current time step
  integer, intent(out) :: dim_obs   !< Dimension of full observation vector

! *** Local variables ***
  integer :: dim_obs_type    ! Full number of observation of one type

!   integer :: dim_obs_sst     ! Full number of SST observations
!   integer :: dim_obs_sss_cci ! Full number of SSS (CCI) observations
!   integer :: dim_obs_sss     ! Full number of SSS (SMOS) observations
!   integer :: dim_obs_ssh     ! Full number of SSH observations
!   integer :: dim_obs_prof    ! Full number of subsurface profile observations
!   integer :: dim_obs_en4ana  ! Full number of EN4 analysis profile observations
!   
!   integer :: dim_obs_chl_cci ! Full number of biogeochem.-TYPE observations
!   integer :: dim_obs_DIC_glodap
!   integer :: dim_obs_Alk_glodap
!   integer :: dim_obs_pCO2_socat
!   integer :: dim_obs_O2_comf
!   integer :: dim_obs_N_comf
!   integer :: dim_obs_O2_argo
!   integer :: dim_obs_N_argo
!   integer :: dim_obs_O2_merged
!   integer :: dim_obs_n_merged


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

!   ! Initialize number of observations
!   dim_obs_sst = 0
!   dim_obs_sss = 0
!   dim_obs_sss_cci = 0
!   dim_obs_ssh = 0
!   dim_obs_prof = 0
!   dim_obs_en4ana = 0
!   dim_obs_chl_cci = 0
!   dim_obs_DIC_glodap = 0
!   dim_obs_Alk_glodap = 0
!   dim_obs_pCO2_socat = 0
!   dim_obs_O2_comf = 0
!   dim_obs_N_comf = 0
!   dim_obs_O2_argo = 0
!   dim_obs_N_argo = 0
!   dim_obs_O2_merged = 0
!   dim_obs_n_merged = 0


  ! Call observation specific routines
  ! The routines are independent, so it is not relevant
  ! in which order they are called

  dim_obs = 0
  if (assim_o_sst) then
     call init_dim_obs_sst(step, dim_obs_type)
     dim_obs = dim_obs + dim_obs_type
  end if
  if (assim_o_sss) then
     call init_dim_obs_sss(step, dim_obs_type)
     dim_obs = dim_obs + dim_obs_type
  end if
  if (assim_o_sss_cci) then
     call init_dim_obs_sss_cci(step, dim_obs_type)
     dim_obs = dim_obs + dim_obs_type
  end if
  if (assim_o_ssh) then
     call init_dim_obs_ssh(step, dim_obs_type)
     dim_obs = dim_obs + dim_obs_type
  end if
  if (assim_o_en4_t .or. assim_o_en4_s) then
     call init_dim_obs_prof(step, dim_obs_type)
     dim_obs = dim_obs + dim_obs_type
  end if
  
  if (assim_o_chl_cci) then
     call init_dim_obs_chl_cci(step, dim_obs_type)
     dim_obs = dim_obs + dim_obs_type
  end if
  if (assim_o_DIC_glodap) then
     call init_dim_obs_DIC_glodap(step, dim_obs_type)
     dim_obs = dim_obs + dim_obs_type
  end if
  if (assim_o_Alk_glodap) then
     call init_dim_obs_Alk_glodap(step, dim_obs_type)
     dim_obs = dim_obs + dim_obs_type
  end if
  if (assim_o_pCO2_SOCAT) then
     call init_dim_obs_pCO2_SOCAT(step, dim_obs_type)
     dim_obs = dim_obs + dim_obs_type
  end if
  if (assim_o_O2_comf) then
     call init_dim_obs_O2_comf(step, dim_obs_type)
     dim_obs = dim_obs + dim_obs_type
  end if
  if (assim_o_N_comf) then
     call init_dim_obs_N_comf(step, dim_obs_type)
     dim_obs = dim_obs + dim_obs_type
  end if
  if (assim_o_O2_argo) then
     call init_dim_obs_O2_argo(step, dim_obs_type)
     dim_obs = dim_obs + dim_obs_type
  end if
  if (assim_o_N_argo) then
     call init_dim_obs_N_argo(step, dim_obs_type)
     dim_obs = dim_obs + dim_obs_type
  end if
  if (assim_o_O2_merged) then
     call init_dim_obs_O2_merged(step, dim_obs_type)
     dim_obs = dim_obs + dim_obs_type
  end if
  if (assim_o_n_merged) then
     call init_dim_obs_n_merged(step, dim_obs_type)
     dim_obs = dim_obs + dim_obs_type
  end if

!   if (assim_o_sst)     call init_dim_obs_sst(step, dim_obs_sst)
!   if (assim_o_sss)     call init_dim_obs_sss(step, dim_obs_sss)
!   if (assim_o_sss_cci) call init_dim_obs_sss_cci(step, dim_obs_sss_cci)
!   if (assim_o_ssh)     call init_dim_obs_ssh(step, dim_obs_ssh)
!   if (assim_o_en4_t .or. assim_o_en4_s) call init_dim_obs_prof(step, dim_obs_prof)
!   
!   if (assim_o_chl_cci)      call init_dim_obs_chl_cci(step, dim_obs_chl_cci)
!   if (assim_o_DIC_glodap)   call init_dim_obs_DIC_glodap(step, dim_obs_DIC_glodap)
!   if (assim_o_Alk_glodap)   call init_dim_obs_Alk_glodap(step, dim_obs_Alk_glodap)
!   if (assim_o_pCO2_SOCAT)   call init_dim_obs_pCO2_SOCAT(step, dim_obs_pCO2_socat)
!   if (assim_o_O2_comf)      call init_dim_obs_O2_comf(step, dim_obs_O2_comf)
!   if (assim_o_N_comf)       call init_dim_obs_N_comf(step, dim_obs_N_comf)
!   if (assim_o_O2_argo)      call init_dim_obs_O2_argo(step, dim_obs_O2_argo)
!   if (assim_o_N_argo)       call init_dim_obs_N_argo(step, dim_obs_N_argo)
!   if (assim_o_O2_merged)    call init_dim_obs_O2_merged(step, dim_obs_O2_merged)
!   if (assim_o_n_merged)     call init_dim_obs_n_merged(step, dim_obs_n_merged)
! 
!   dim_obs =   dim_obs_sst + dim_obs_sss + dim_obs_sss_cci + dim_obs_ssh + dim_obs_prof + dim_obs_en4ana &
!             + dim_obs_chl_cci + dim_obs_DIC_glodap + dim_obs_Alk_glodap + dim_obs_pCO2_socat &
!             + dim_obs_O2_comf + dim_obs_N_comf + dim_obs_o2_argo + dim_obs_N_argo &
!             + dim_obs_O2_merged + dim_obs_n_merged


  ! *** Generate profile observation files ***
  if (proffiles_o == 1) then
     ! Generate distributed files
     call init_dim_obs_f_proffile_pdaf(start_year_o,end_year_o)
          
  elseif (proffiles_o == 2) then
     ! Generate one global file
     write(*,*) 'Generation of global file from EN4 raw profile data ', &
                'not yet implemented in FESOM Version 2.0 - stopping!'
     call abort_parallel
!~      CALL init_dim_obs_f_proffile_g_pdaf(step,dim_obs_prof)
  end if

end subroutine init_dim_obs_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for obs_op
!!
!! This routine calls the observation-specific
!! routines obs_op_TYPE.
!!
subroutine obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)

  ! Include functions for different observations
  use obs_sst_pdafomi, only: obs_op_sst
  use obs_sss_smos_pdafomi, only: obs_op_sss
  use obs_sss_cci_pdafomi, only: obs_op_sss_cci
  use obs_ssh_cmems_pdafomi, only: obs_op_ssh
  use obs_TSprof_EN4_pdafomi, only: obs_op_prof
  
  use obs_chl_cci_pdafomi, only: obs_op_chl_cci
  use obs_DIC_glodap_pdafomi, only: obs_op_DIC_glodap
  use obs_Alk_glodap_pdafomi, only: obs_op_Alk_glodap
  use obs_pco2_SOCAT_pdafomi, only: obs_op_pCO2_SOCAT
  use obs_O2_comf_pdafomi, only: obs_op_O2_comf
  use obs_N_comf_pdafomi, only: obs_op_N_comf
  use obs_O2_argo_pdafomi, only: obs_op_o2_argo
  use obs_N_argo_pdafomi, only: obs_op_N_argo
  use obs_O2_merged_pdafomi, only: obs_op_o2_merged
  use obs_n_merged_pdafomi, only: obs_op_n_merged

  implicit none

! *** Arguments ***
  integer, intent(in) :: step                 !< Current time step
  integer, intent(in) :: dim_p                !< PE-local state dimension
  integer, intent(in) :: dim_obs              !< Dimension of full observed state
  real, intent(in)    :: state_p(dim_p)       !< PE-local model state
  real, intent(inout) :: ostate(dim_obs)      !< PE-local full observed state


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

  ! The order of the calls determines how the different observations
  ! are ordered in the full state vector
  call obs_op_sst    (dim_p, dim_obs, state_p, ostate)
  call obs_op_sss    (dim_p, dim_obs, state_p, ostate)
  call obs_op_sss_cci(dim_p, dim_obs, state_p, ostate)
  call obs_op_ssh    (dim_p, dim_obs, state_p, ostate)
  call obs_op_prof   (dim_p, dim_obs, state_p, ostate)
  
  call obs_op_chl_cci   (dim_p, dim_obs, state_p, ostate)
  call obs_op_DIC_glodap(dim_p, dim_obs, state_p, ostate)
  call obs_op_Alk_glodap(dim_p, dim_obs, state_p, ostate)
  call obs_op_pCO2_SOCAT(dim_p, dim_obs, state_p, ostate)
  call obs_op_O2_comf   (dim_p, dim_obs, state_p, ostate)
  call obs_op_N_comf    (dim_p, dim_obs, state_p, ostate)
  call obs_op_O2_argo   (dim_p, dim_obs, state_p, ostate)
  call obs_op_N_argo    (dim_p, dim_obs, state_p, ostate)
  call obs_op_O2_merged (dim_p, dim_obs, state_p, ostate)
  call obs_op_n_merged  (dim_p, dim_obs, state_p, ostate)

end subroutine obs_op_pdafomi


!-------------------------------------------------------------------------------
!> Call-back routine for init_dim_obs_l
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type
!!
subroutine init_dim_obs_l_pdafomi(domain_p, step, dim_obs, dim_obs_l)

  ! Include observation types:
  use obs_sst_pdafomi, only: init_dim_obs_l_sst
  use obs_sss_smos_pdafomi, only: init_dim_obs_l_sss
  use obs_sss_cci_pdafomi, only: init_dim_obs_l_sss_cci
  use obs_ssh_cmems_pdafomi, only: init_dim_obs_l_ssh
  use obs_TSprof_EN4_pdafomi, only: init_dim_obs_l_prof
  
  use obs_chl_cci_pdafomi, only: init_dim_obs_l_chl_cci
  use obs_DIC_glodap_pdafomi, only: init_dim_obs_l_DIC_glodap
  use obs_Alk_glodap_pdafomi, only: init_dim_obs_l_Alk_glodap
  use obs_pco2_SOCAT_pdafomi, only: init_dim_obs_l_pCO2_SOCAT
  use obs_o2_comf_pdafomi, only: init_dim_obs_l_o2_comf
  use obs_n_comf_pdafomi, only: init_dim_obs_l_n_comf
  use obs_o2_argo_pdafomi, only: init_dim_obs_l_o2_argo
  use obs_n_argo_pdafomi, only: init_dim_obs_l_n_argo
  use obs_o2_merged_pdafomi, only: init_dim_obs_l_o2_merged
  use obs_n_merged_pdafomi, only: init_dim_obs_l_n_merged

  implicit none

! *** Arguments ***
  integer, intent(in)  :: domain_p   !< Index of current local analysis domain
  integer, intent(in)  :: step       !< Current time step
  integer, intent(in)  :: dim_obs    !< Full dimension of observation vector
  integer, intent(out) :: dim_obs_l  !< Local dimension of observation vector


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

   call init_dim_obs_l_sst    (domain_p, step, dim_obs, dim_obs_l)
   call init_dim_obs_l_sss    (domain_p, step, dim_obs, dim_obs_l)
   call init_dim_obs_l_sss_cci(domain_p, step, dim_obs, dim_obs_l)
   call init_dim_obs_l_ssh    (domain_p, step, dim_obs, dim_obs_l)
   call init_dim_obs_l_prof   (domain_p, step, dim_obs, dim_obs_l)
   
   call init_dim_obs_l_chl_cci   (domain_p, step, dim_obs, dim_obs_l)
   call init_dim_obs_l_DIC_glodap(domain_p, step, dim_obs, dim_obs_l)
   call init_dim_obs_l_Alk_glodap(domain_p, step, dim_obs, dim_obs_l)
   call init_dim_obs_l_pCO2_SOCAT(domain_p, step, dim_obs, dim_obs_l)
   call init_dim_obs_l_o2_comf   (domain_p, step, dim_obs, dim_obs_l)
   call init_dim_obs_l_n_comf    (domain_p, step, dim_obs, dim_obs_l)
   call init_dim_obs_l_o2_argo   (domain_p, step, dim_obs, dim_obs_l)
   call init_dim_obs_l_n_argo    (domain_p, step, dim_obs, dim_obs_l)
   call init_dim_obs_l_o2_merged (domain_p, step, dim_obs, dim_obs_l)
   call init_dim_obs_l_n_merged  (domain_p, step, dim_obs, dim_obs_l)

end subroutine init_dim_obs_l_pdafomi

!--------------------------------------------
!> Determine which type of observation are assimilated
!!
!! This routine sets the variables assimilatePHY 
!! and assimilateBGC indicating whether observations
!! of physics or biogeochemistry are assimilated.
!! 
!! This routine is an add-on. It is called from
!! coupled_da_mod/cda_set_sweeps and does not call an
!! observation module. However, it is places here because
!! it is important to adapt it if an observation type is added.
!!
subroutine check_coupled_da_pdafomi(assimilatePHY, assimilateBGC)

  ! physics observations
  use obs_sst_pdafomi, &
       only: assim_o_sst
  use obs_sss_smos_pdafomi, &
       only: assim_o_sss
  use obs_sss_cci_pdafomi, &
       only: assim_o_sss_cci
  use obs_ssh_cmems_pdafomi, &
       only: assim_o_ssh
  use obs_TSprof_EN4_pdafomi, &
       only: ASSIM_O_en4_t, ASSIM_O_en4_s 

  ! biogeochem observations    
  use obs_chl_cci_pdafomi, &
       only: assim_o_chl_cci
  use obs_DIC_glodap_pdafomi, &
       only: assim_o_DIC_glodap
  use obs_Alk_glodap_pdafomi, &
       only: assim_o_Alk_glodap
  use obs_pco2_SOCAT_pdafomi, &
       only: assim_o_pCO2_SOCAT
  use obs_o2_comf_pdafomi, &
       only: assim_o_o2_comf
  use obs_n_comf_pdafomi, &
       only: assim_o_n_comf
  use obs_o2_argo_pdafomi, &
       only: assim_o_o2_argo
  use obs_N_argo_pdafomi, &
       only: assim_o_N_argo
  use obs_o2_merged_pdafomi, &
       only: assim_o_o2_merged
  use obs_N_merged_pdafomi, &
       only: assim_o_N_merged

  implicit none

! *** Arguments ***
  logical, intent(out) :: assimilatePHY
  logical, intent(out) :: assimilateBGC


  ! Check whether any physics observation is assimilated
  if (assim_o_sst .or. assim_o_sss .or. assim_o_sss_cci .or. assim_o_ssh &
       .or. ASSIM_O_en4_t .or. ASSIM_O_en4_s) then
     assimilatePHY = .true.
  else
     assimilatePHY = .false.
  end if

  ! Check whether any BGC observation is assimilated
  if (assim_o_chl_cci .or. assim_o_DIC_glodap .or. assim_o_Alk_glodap &
       .or. assim_o_pCO2_SOCAT .or. assim_o_o2_comf .or. assim_o_n_comf &
       .or. assim_o_o2_argo .or. assim_o_N_argo .or. assim_o_o2_merged &
       .or. assim_o_N_merged) then
     assimilateBGC = .true.
  else
     assimilateBGC = .false.
  end if

end subroutine check_coupled_da_pdafomi
  
