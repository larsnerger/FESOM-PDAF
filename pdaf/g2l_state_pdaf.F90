!> Routine to map a model state to a local analysis domain
!!
!! The routine is called during the loop over all
!! local analysis domains in the domain local filters
!! before the analysis on a single local analysis 
!! domain. It has to project the full process-local 
!! model state onto the current local analysis 
!! domain.
!!
!! The routine is called by each filter process.
!!
!! __Revision history:__
!! 2005-09 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
subroutine g2l_state_pdaf(step, domain_p, dim_p, state_p, dim_l, state_l)

  use assim_pdaf_mod, &                 ! Variables for assimilation
       only: id_lstate_in_pstate

  implicit none

! *** Arguments ***
  integer, intent(in) :: step           !< Current time step
  integer, intent(in) :: domain_p       !< Current local analysis domain
  integer, intent(in) :: dim_p          !< PE-local full state dimension
  integer, intent(in) :: dim_l          !< Local state dimension
  real, intent(in)    :: state_p(dim_p) !< PE-local full state vector 
  real, intent(out)   :: state_l(dim_l) !< State vector on local analysis domain
  
! *** Local variables *** 
  integer :: i                          ! Counter


! *************************************
! *** Initialize local state vector ***
! *************************************
 
  do i = 1, dim_l
     state_l(i) = state_p(id_lstate_in_pstate(i))
  enddo

end subroutine g2l_state_pdaf
