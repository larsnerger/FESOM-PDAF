!$Id: l2g_state_pdaf.F90 2466 2021-02-25 12:46:34Z lnerger $
!>  Routine to initialize full state from local analysis
!!
!! User-supplied call-back routine for PDAF.
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
!! * Later revisions - see repository log
!!
SUBROUTINE l2g_state_pdaf(step, domain, dim_l, state_l, dim_p, state_p)

  USE mod_assim_pdaf, &           ! Variables for assimilation
       ONLY: id_lstate_in_pstate, isweep, &
             cda_bio, cda_phy, type_sweep
  USE statevector_pdaf, &
       ONLY: nfields, sfields, sfields_l

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step           !< Current time step
  INTEGER, INTENT(in) :: domain         !< Current local analysis domain
  INTEGER, INTENT(in) :: dim_l          !< Local state dimension
  INTEGER, INTENT(in) :: dim_p          !< Process-local full state dimension
  REAL, INTENT(in)    :: state_l(dim_l) !< State vector on local analysis domain
  REAL, INTENT(inout) :: state_p(dim_p) !< Process-local full state vector 
  
! *** Local variables *** 
  INTEGER :: i, ifield            !< Counters
  LOGICAL :: update_cda           !< Whether to perform DA update
  
  
! **************************************************
! *** Initialize elements of global state vector ***
! **************************************************
  
  DO ifield = 1, nfields
  
     ! Determine whether to apply update according to coupled data assimilation settings
     
     IF (.NOT.(sfields(ifield)%bgc) .AND. (TRIM(type_sweep(isweep))=='phy')) THEN
        ! Physics field and physics sweep:
        update_cda = .TRUE.
     ELSEIF (sfields(ifield)%bgc .AND. (TRIM(type_sweep(isweep))=='bio')) THEN
        ! BGC field and BGC sweep:
        update_cda = .TRUE.
     ELSE
        ! Strongly coupled DA configuration:
        IF (type_sweep(isweep)=='phy' .AND. TRIM(cda_phy)=='strong') THEN
           update_cda = .TRUE.
        ELSEIF (type_sweep(isweep)=='bio' .AND. TRIM(cda_bio)=='strong') THEN
           update_cda = .TRUE.
        ELSE
           ! Weak coupling and unequal type of field and sweep:
           update_cda = .FALSE.
        END IF
     END IF
  
     IF (update_cda) THEN
        ! update field in global state vector from local state
        DO i = sfields_l(ifield)%off + 1, sfields_l(ifield)%off + sfields_l(ifield)%dim
           state_p(id_lstate_in_pstate(i)) = state_l(i)
        END DO
     ENDIF
     
  END DO

END SUBROUTINE l2g_state_pdaf
