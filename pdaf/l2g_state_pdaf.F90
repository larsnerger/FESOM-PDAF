!>  Routine to initialize full state from local analysis
!!
!! The routine is called during the loop over all
!! local analysis domains in the domain local filters
!! after the analysis and ensemble transformation 
!! on a single local analysis domain. It has to 
!! initialize elements of the PE-local full state 
!! vector from the provided analysis state vector 
!! on the local analysis domain.
!!
!! The routine is called by each filter process.
!!
!! __Revision history:__
!! 2005-09 - Lars Nerger - Initial code
!! ~2022   - Frauke B - adaption for REcoM
!!
subroutine l2g_state_pdaf(step, domain, dim_l, state_l, dim_p, state_p)

  use assim_pdaf_mod, &                 ! Variables for assimilation
       only: id_lstate_in_pstate, isweep, &
             cda_bio, cda_phy, type_sweep
  use statevector_pdaf, &               ! Statevector variables
       only: nfields, sfields, sfields_l

  implicit none

! *** Arguments ***
  integer, intent(in) :: step           !< Current time step
  integer, intent(in) :: domain         !< Current local analysis domain
  integer, intent(in) :: dim_l          !< Local state dimension
  integer, intent(in) :: dim_p          !< Process-local full state dimension
  real, intent(in)    :: state_l(dim_l) !< State vector on local analysis domain
  real, intent(inout) :: state_p(dim_p) !< Process-local full state vector 
  
! *** Local variables *** 
  integer :: i, ifield                  ! Counters
  logical :: update_cda                 ! Whether to apply DA update


! **************************************************
! *** Initialize elements of global state vector ***
! **************************************************
  
  do ifield = 1, nfields
  
     ! Determine whether to apply update according to coupled data assimilation settings
     
     if (.not.(sfields(ifield)%bgc) .and. (trim(type_sweep(isweep))=='phy')) then
        ! Physics field and physics sweep:
        update_cda = .true.
     elseif (sfields(ifield)%bgc .and. (trim(type_sweep(isweep))=='bio')) then
        ! BGC field and BGC sweep:
        update_cda = .true.
     else
        ! Strongly coupled DA configuration:
        if (type_sweep(isweep)=='phy' .and. trim(cda_phy)=='strong') then
           update_cda = .true.
        elseif (type_sweep(isweep)=='bio' .and. trim(cda_bio)=='strong') then
           update_cda = .true.
        else
           ! Weak coupling and unequal type of field and sweep:
           update_cda = .false.
        end if
     end if
  
     if (update_cda) then
        ! update field in global state vector from local state
        do i = sfields_l(ifield)%off + 1, sfields_l(ifield)%off + sfields_l(ifield)%dim
           state_p(id_lstate_in_pstate(i)) = state_l(i)
        end do
     endif
     
  end do

end subroutine l2g_state_pdaf
