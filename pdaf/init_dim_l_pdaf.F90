!> Set dimension of local model state
!!
!! The routine is called by PDAF during the
!! analysis step in the loop over all local
!! analysis domains. It has to set the dimension
!! of the local model state on the current analysis
!! domain. In addition, the coordinates of this
!! domain are stored and the index arrays for the
!! local state vector and the mapping between global
!! to local state vectors are initialized.
!!
!! The routine is called by each filter process.
!!
!! __Revision history:__
!! * 2005-09 - Lars Nerger - Initial code
!! * 2022-03 - Frauke B    - Adapted for FESOM 2.1
!! * 2025-12 - Lars Nerger - Revision for PDAF3
!!
SUBROUTINE init_dim_l_pdaf(step, nsweeped_domain_p, dim_l)

  USE mod_assim_pdaf, &
       ONLY: id_lstate_in_pstate, coords_l, isweep
  USE mod_parallel_pdaf, &
       ONLY: abort_parallel
  USE fesom_pdaf, &
       ONLY: mesh_fesom, nlmax, r2g
  USE statevector_pdaf, &
       ONLY: id, nfields, sfields, sfields_l, &
       bgcmin, bgcmax, phymin, phymax
  USE fesom_pdaf, &
       ONLY: myDim_nod2D

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step              !< Current time step
  INTEGER, INTENT(in)  :: nsweeped_domain_p !< Current local analysis domain, containing repititive sweeps
  INTEGER, INTENT(out) :: dim_l             !< Local state dimension

! *** Local variables ***
  INTEGER :: i, b, id_var                        ! Counters
  INTEGER :: nlay                                ! Number of layers for current domain
  INTEGER :: domain_p                            ! Local analysis domain accounting for multiple sweeps
  
  
! ********************************************************
! ***  Account for multi sweeps in local analysis loop ***
! ********************************************************

  IF (nsweeped_domain_p <= myDim_nod2D) THEN
     domain_p = nsweeped_domain_p
     isweep = 1
  ELSE
     domain_p = nsweeped_domain_p - myDim_nod2D
     isweep = 2
  END IF


! ****************************************
! *** Initialize local state dimension ***
! ****************************************

  ! The local state vector only contains fields that are updated
  
  ! Allocate array
  IF (ALLOCATED(sfields_l)) DEALLOCATE(sfields_l)
  ALLOCATE(sfields_l(nfields))

  nlay = mesh_fesom%nlevels_nod2D(domain_p)-1
  
  IF (nlay > nlmax) THEN
     WRITE(*,*) 'FESOM-PDAF ', 'init_dim_l_pdaf ', 'domain_p ', domain_p, ' nlay exceeds layer bounds!'
     CALL abort_parallel()
  ENDIF
  
  ! Physics:
  DO i = phymin, phymax

     IF (sfields(i)%updated) THEN
        ! surface fields:
        IF (sfields(i)%ndims == 1) sfields_l(i)%dim = 1
        ! 3D fields:
        IF (sfields(i)%ndims == 2) sfields_l(i)%dim = nlay
     ELSE
        ! not updated:
        sfields_l(i)%dim = 0
     ENDIF
  ENDDO

  ! BGC:
  DO i = bgcmin, bgcmax
     IF (sfields(i)%updated) THEN
        ! surface fields:
        IF (sfields(i)%ndims == 1) sfields_l(i)%dim = 1
        ! 3D fields:
        IF (sfields(i)%ndims == 2) sfields_l(i)%dim = nlay
     ELSE
        ! not updated:
        sfields_l(i)%dim = 0
     ENDIF
  ENDDO

  ! Set local offsets
  sfields_l(1)%off = 0
  DO i = 2, nfields
     sfields_l(i)%off = sfields_l(i-1)%off + sfields_l(i-1)%dim
  END DO

  ! *** Local state dimension
  dim_l = SUM(sfields_l(:)%dim)


! **********************************************
! *** Initialize coordinates of local domain ***
! **********************************************

  ! Get location of current water column (basis point)
  CALL r2g(coords_l(1), coords_l(2), mesh_fesom%coord_nod2D(1, domain_p), mesh_fesom%coord_nod2D(2, domain_p))
  

! ****************************************************
! *** Initialize array of indices for local domain ***
! ****************************************************

  ! Allocate array
  IF (ALLOCATED(id_lstate_in_pstate)) DEALLOCATE(id_lstate_in_pstate)
  ALLOCATE(id_lstate_in_pstate(dim_l))

  ! *** indices for full state vector ***

  ! SSH
  IF (sfields(id%ssh)%updated) THEN
     id_lstate_in_pstate (sfields_l(id%ssh)%off+1) &
          = sfields(id%ssh)%off + domain_p
  ENDIF
  
  ! U
  id_var = id%u
  IF (sfields(id_var)%updated) THEN
     DO i = 1, sfields_l(id_var)%dim
        id_lstate_in_pstate (sfields_l(id_var)%off + i) = &
            sfields(id_var)%off + (domain_p-1)*(nlmax) + i 
     END DO
  ENDIF
  
  ! V
  id_var = id%v
  IF (sfields(id_var)%updated) THEN
     DO i = 1, sfields_l(id_var)%dim
        id_lstate_in_pstate (sfields_l(id_var)%off + i) = &
            sfields(id_var)%off + (domain_p-1)*(nlmax) + i 
     END DO
  ENDIF
        
  ! W
  ! id_lstate_in_pstate (sfields_l(id%w)%off+1 : sfields_l(id%w+1)%off) &
  !      = sfields(id%w)%off &
  !      + (domain_p-1)*(nlmax) &
  !      + (/(i, i=1,sfields_l(id%w)%dim)/)
  
  ! Temp
  id_var = id%temp
  IF (sfields(id_var)%updated) THEN
     DO i = 1, sfields_l(id_var)%dim
        id_lstate_in_pstate (sfields_l(id_var)%off + i) = &
            sfields(id_var)%off + (domain_p-1)*(nlmax) + i 
     END DO
  ENDIF
  
  ! Salt
  id_var = id%salt
  IF (sfields(id_var)%updated) THEN
     DO i = 1, sfields_l(id_var)%dim
        id_lstate_in_pstate (sfields_l(id_var)%off + i) = &
            sfields(id_var)%off + (domain_p-1)*(nlmax) + i 
     END DO
  ENDIF
        
  ! BGC:
  DO b = bgcmin, bgcmax
  
     ! only updated fields:
     IF ((sfields(b)% updated)) THEN
        IF (sfields(b)%ndims == 1)   THEN

           ! surface fields:
           id_lstate_in_pstate (sfields_l(b)%off+1) &
                = sfields(b)%off + domain_p

        ELSEIF (sfields(b)%ndims == 2)   THEN

           ! 3D fields:
           id_var = b
           IF (sfields(id_var)%updated) THEN
              DO i = 1, sfields_l(id_var)%dim
                 id_lstate_in_pstate (sfields_l(id_var)%off + i) = &
                      sfields(id_var)%off + (domain_p-1)*(nlmax) + i 
              END DO
           ENDIF
!            id_lstate_in_pstate (sfields_l(i)%off + 1 : sfields_l(i)%off + sfields_l(i)%dim)&
!                 = sfields(i)%off &
!                 + (domain_p-1)*(nlmax) &
!                 + (/(k, k=1, sfields_l(i)%dim)/)
        ENDIF
     ENDIF
  ENDDO

END SUBROUTINE init_dim_l_pdaf
