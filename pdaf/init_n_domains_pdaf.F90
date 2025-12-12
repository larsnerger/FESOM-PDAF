!>  Routine to set number of local analysis domains
!!
!! User-supplied call-back routine for PDAF.
!!
!! The routine is called in PDAF_X_update 
!! or PDAF_assimialte_X at the beginning of the 
!! analysis step before the loop through all local
!! analysis domains. It has to set the number of 
!! local analysis domains for the process-local domain.
!!
!! The routine is called by all filter processes.
!!
!! __Revision history:__
!! 2005-09 - Lars Nerger - Initial code for AWI-CM
!!
subroutine init_n_domains_pdaf(step, n_domains_p)

  use g_parsup, &
      only: myDim_nod2D, myDim_elem2D
  use parallel_pdaf_mod, &
      only: mype_filter
  use coupled_da_mod, &
      only: n_sweeps

  implicit none

! *** Arguments ***
  integer, intent(in)  :: step        ! Current time step
  integer, intent(out) :: n_domains_p ! Process-local number of analysis domains
  
! ************************************
! *** Initialize number of domains ***
! ************************************

  n_domains_p = n_sweeps * myDim_nod2D

end subroutine init_n_domains_pdaf
