!> Print information to screen
!!
!! This routine performs a model-sided screen output about
!! the coniguration of the data assimilation system.
!! Using this output is optional. Most of the information
!! is also displayed by PDAF itself when it is initialized
!! in PDAF_init. Not displayed by PDAF is the assimilation
!! interval (delt_obs_ocn), which is unknown to PDAF.
!!
!! __Revision history:__
!! * 2011-05 - Lars Nerger - Initial code extracted from init_pdaf
!! * 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
!!
subroutine init_pdaf_info()

  use assim_pdaf_mod, & ! Variables for assimilation
       only: filtertype, subtype, dim_ens,  &
       forget, delt_obs_ocn, dim_state, &
       file_syntobs, twin_experiment
    
  implicit none
    
    
! *****************************
! *** Initial Screen output ***
! *****************************
    
  if (filtertype == 1) then
     write (*, '(a,21x, a)') 'FESOM-PDAF', 'Filter: SEIK'
     if (subtype == 2) then
        write (*, '(a,6x, a)') 'FESOM-PDAF', '-- fixed error-space basis'
     else if (subtype == 3) then
        write (*, '(a,6x, a)') 'FESOM-PDAF', '-- fixed state covariance matrix'
     else if (subtype == 4) then
        write (*, '(a,6x, a)') 'FESOM-PDAF', '-- use ensemble transformation'
     else if (subtype == 5) then
        write (*, '(a,6x, a)') 'FESOM-PDAF', '-- Offline mode'
     end if
  else if (filtertype == 2) then
     write (*, '(a,21x, a)') 'FESOM-PDAF', 'Filter: EnKF'
     if (subtype == 5) then
        write (*, '(a,6x, a)') 'FESOM-PDAF', '-- Offline mode'
     end if
  else if (filtertype == 3) then
     write (*, '(a,21x, a)') 'FESOM-PDAF', 'Filter: LSEIK'
     if (subtype == 2) then
        write (*, '(a,6x, a)') 'FESOM-PDAF', '-- fixed error-space basis'
     else if (subtype == 3) then
        write (*, '(a,6x, a)') 'FESOM-PDAF', '-- fixed state covariance matrix'
     else if (subtype == 4) then
        write (*, '(a,6x, a)') 'FESOM-PDAF', '-- use ensemble transformation'
     else if (subtype == 5) then
        write (*, '(a,6x, a)') 'FESOM-PDAF', '-- Offline mode'
     end if
  else if (filtertype == 4) then
     write (*, '(a,21x, a)') 'FESOM-PDAF', 'Filter: ETKF'
     if (subtype == 0) then
        write (*, '(a,6x, a)') 'FESOM-PDAF', '-- Variant using T-matrix'
     else if (subtype == 1) then
        write (*, '(a,6x, a)') 'FESOM-PDAF', '-- Variant following Hunt et al. (2007)'
     else if (subtype == 5) then
        write (*, '(a,6x, a)') 'FESOM-PDAF', '-- Offline mode'
     end if
  else if (filtertype == 5) then
     write (*, '(a,21x, a)') 'FESOM-PDAF', 'Filter: LETKF'
     if (subtype == 0) then
        write (*, '(a,6x, a)') 'FESOM-PDAF', '-- Variant using T-matrix'
     else if (subtype == 1) then
        write (*, '(a,6x, a)') 'FESOM-PDAF', '-- Variant following Hunt et al. (2007)'
     else if (subtype == 5) then
        write (*, '(a,6x, a)') 'FESOM-PDAF', '-- Offline mode'
     end if
  else if (filtertype == 6) then
     write (*, '(a,21x, a)') 'FESOM-PDAF', 'Filter: ESTKF'
     if (subtype == 0) then
        write (*, '(a,6x, a)') 'FESOM-PDAF', '-- Standard mode'
     else if (subtype == 5) then
        write (*, '(a,6x, a)') 'FESOM-PDAF', '-- Offline mode'
     end if
  else if (filtertype == 7) then
     write (*, '(a,21x, a)') 'FESOM-PDAF', 'Filter: LESTKF'
     if (subtype == 0) then
        write (*, '(a,6x, a)') 'FESOM-PDAF', '-- Standard mode'
     else if (subtype == 5) then
        write (*, '(a,6x, a)') 'FESOM-PDAF', '-- Offline mode'
     end if
  end if
  write (*, '(a, 14x, a, i5)') 'FESOM-PDAF','ensemble size:', dim_ens
  if (subtype /= 5) write (*, '(a, 6x, a, i5)') 'FESOM-PDAF','Assimilation interval:', delt_obs_ocn
  write (*, '(a, 10x, a, f5.2)') 'FESOM-PDAF','forgetting factor:', forget
  if (twin_experiment) &
       write (*, '(/a,6x, a)') 'FESOM-PDAF','Run twin experiment with synthetic observations'
  if (filtertype==100 .or. twin_experiment) &
       write (*, '(a,11x, a, a)') 'FESOM-PDAF','File for synthetic observations: ', trim(file_syntobs)
  write (*, '(a, 8x, a, 1x, i9)') 'FESOM-PDAF','FESOM state dimension:',dim_state


end subroutine init_pdaf_info
