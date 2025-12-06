!>  Initialize model fields from state vector
!!
!! During the forecast phase of the filter this subroutine
!! is called from PDAF_init_forecast or PDAF3_assimilate.
!! supplying a model state, which has to be evolved. 
!! The routine has to initialize the fields of the 
!! model (typically available through a module) from 
!! the state vector of PDAF. With parallelization, 
!! MPI communication might be required to 
!! initialize all subdomains on the model PEs.
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! __Revision history:__
!! * 2004-10 - Lars Nerger - Initial code
!! * 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
!! * 2020-03 - Frauke B     - Added velocities, removed sea-ice (FESOM2.0)
!! * 2025-12 - Lars Nerger - Revision for PDAF3
!!
SUBROUTINE distribute_state_pdaf(dim_p, state_p)

  USE mod_parallel_pdaf, &
       ONLY: mype_submodel, mype_world, task_id, mype_model, npes_model, &
       COMM_model, MPIerr, mype_filter, writepe
  USE mod_assim_pdaf, &
       ONLY: this_is_pdaf_restart, start_from_ENS_spinup, &
       istep_asml, step_null
  USE statevector_pdaf, &
       ONLY: id, sfields, nfields
  USE fesom_pdaf, &
       ONLY: daynew, timenew, nlmax, mesh_fesom, topography_p, &
       mydim_nod2d, myDim_elem2D, eDim_nod2D, eDim_elem2D, &
       eta_n, uv, wvel, tr_arr, unode, a_ice, &
       exchange_nod, exchange_elem

  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  !< PE-local state vector

! *** Local variables ***
  INTEGER :: i, k, b, s, istate, ifesom  ! Counters
  INTEGER :: node                        ! Node index
  REAL, ALLOCATABLE :: U_node_upd(:,:,:) ! Velocity update on nodes
  REAL, ALLOCATABLE :: U_elem_upd(:,:,:) ! Velocity update on elements
  LOGICAL, SAVE :: first_call = .TRUE.   ! Indicator for first call

! Debugging:
  LOGICAL            :: debugmode
  INTEGER            :: fileID_debug
  CHARACTER(len=5)   :: mype_string, tim_string
  CHARACTER(len=3)   :: day_string

  
  ! Set debug output
  debugmode = .FALSE.

! **********************
! *** Initialization ***
! **********************

  ! no need to distibute the ensemble in case of restart, just skip this routine:

  do_dist: IF (this_is_pdaf_restart .AND. (istep_asml==step_null)) THEN
     IF (mype_world==0) WRITE(*,'(a,3x,a)') &
          'FESOM-PDAF', 'This is a restart: Skip distribute_state_pdaf at initial step'
    
  ELSEIF (start_from_ENS_spinup .AND. (istep_asml==step_null)) THEN do_dist
     IF (mype_world==0) WRITE(*,'(a,3x,a)') &
          'FESOM-PDAF', 'Start from perturbed ensemble: Skip distribute_state_pdaf at initial step'
    
  ELSE do_dist
     IF (mype_submodel==0) WRITE (*,'(a,3x,a,i5)') &
          'FESOM-PDAF', 'distribute_state_pdaf, task: ', task_id
    
     ! ensure to distribute fields with valid topography at first call
     IF (first_call) THEN
        IF (writepe) WRITE (*,'(a,8x,a)') 'FESOM-PDAF', 'Distribute_state: set topography'
        state_p = state_p * topography_p
        first_call = .FALSE.
     END IF


! *******************************************
! *** Initialize model fields from state  ***
! *** Each model PE knows its sub-state   ***
!********************************************

  ! *** Dimensions of FESOM arrays:
  ! * eta_n          (myDim_nod2D + eDim_nod2D)            ! Dynamic topography
  ! * UV(1,:,:)      (1, nl-1, myDim_elem2D + eDim_elem2D) ! Velocity u
  ! * UV(2,:,:)      (1, nl-1, myDim_elem2D + eDim_elem2D) ! Velocity v
  ! * wvel           (nl,   myDim_nod2D + eDim_nod2D)      ! Velocity w
  ! * tr_arr(:,:,1)  (nl-1, myDim_nod2D + eDim_nod2D, 1)   ! Temperature
  ! * tr_arr(:,:,2)  (nl-1, myDim_nod2D + eDim_nod2D, 1)   ! Salinity
  ! * unode(1,:,:)   (1, nl-1, myDim_nod2D + eDim_nod2D)   ! Velocity u interpolated on nodes
  ! * unode(2,:,:)   (1, nl-1, myDim_nod2D + eDim_nod2D)   ! Velocity v interpolated on nodes
  ! * a_ice          (myDim_nod2D + eDim_nod2D)            ! Sea-ice concentration
  ! ***
  
     ! SSH (1)
     DO i = 1, myDim_nod2D
        eta_n(i) = state_p(i + sfields(id%SSH)%off)
     END DO

     ! u (2) and v (3) velocities
     ! 1. calculate update on nodes, i.e. analysis state (state_p) minus not-yet-updated model state (Unode)
     ALLOCATE(U_node_upd(2, mesh_fesom%nl-1, myDim_nod2D+eDim_nod2D))
     U_node_upd = 0.0

     ! u
     s = sfields(id%u)%off
     DO i = 1, myDim_nod2D
        DO k = 1, nlmax
           s = s + 1
           U_node_upd(1, k, i) = state_p(s) - Unode(1, k, i)
        END DO
     END DO

     ! v
     s = sfields(id%v)%off
     DO i = 1, myDim_nod2D
        DO k = 1, nlmax
           s = s + 1
           U_node_upd(2, k, i) = state_p(s) - Unode(2, k, i)
        END DO
     END DO

      ! 2. interpolate update from nodes to elements
     ALLOCATE(U_elem_upd(2, mesh_fesom%nl-1, myDim_elem2D+eDim_elem2D))
     U_elem_upd = 0.0

     CALL compute_vel_elems(U_node_upd,U_elem_upd)

     ! 3. add update to model velocity on elements (UV)
     UV = UV + U_elem_upd

     ! 4. adjust diagnostic model velocity on nodes (Unode)
     CALL compute_vel_nodes(mesh_fesom)


  ! w (4) velocity: not updated and thus no need to distribute.
  ! DO i = 1, myDim_nod2D
  !  DO k = 1, nlmax
  !      wvel(k,i) = state_p((i-1)*nlmax + k + sfields(id%w)%off) ! w
  !  END DO
  ! END DO
  
  ! Temp and salt are included in tracer field loop.
  ! Sea-ice concentration is needed in PDAF to not assimilate SST at sea-ice locations.
  ! But sea-ice itself is not assimilated, thus sea-ice update is not distributed to the model.

  
! *********************************
! *** Initialize external nodes ***
! *********************************

     CALL exchange_nod(eta_n)            ! SSH
     CALL exchange_elem(UV(:,:,:))       ! u and v (element-wise)


! *********************************
! *** Model 3D tracers          ***
! *********************************

     ! tracer field loop
     DO istate = 1, nfields

        IF (sfields(istate)%trnumfesom>0) THEN

           ifesom = sfields(istate)%trnumfesom  ! index of field in model tracer array

           s = 0
           DO i = 1, myDim_nod2D
              DO k = 1, nlmax
                 s = s + 1
                 tr_arr(k, i, ifesom) = state_p(s + sfields(istate)%off)
              ENDDO
           ENDDO

           ! initialize external nodes
           CALL exchange_nod(tr_arr(:,:,ifesom))

        END IF
     ENDDO 


! ********************
! *** Debug output ***
! ********************

     writedebug: IF (debugmode .AND. mype_world==0) THEN

        ! print state vector
        WRITE(day_string, '(i3.3)') daynew
        WRITE(tim_string, '(i5.5)') INT(timenew)
        fileID_debug=20
        OPEN(unit=fileID_debug, file='distribute_state_pdaf_'//day_string//'_'//tim_string//'.txt', status='unknown')

        ! SSH (1)
        DO i = 1, myDim_nod2D
           WRITE(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') &
                sfields(id%SSH)%variable, i+sfields(id%SSH)%off, eta_n(i), state_p(s)
        END DO

        ! u (2) and v (3) velocities

        s = sfields(id%u)%off
        DO i = 1, myDim_nod2D
           DO k = 1, nlmax
              s = s + 1
              WRITE(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id%u)%variable, s, Unode(1, k, i), state_p(s)
           END DO
        END DO

        s = sfields(id%v)%off
        DO i = 1, myDim_nod2D
           DO k = 1, nlmax
              s = s + 1
              WRITE(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id%v)%variable, s, Unode(2, k, i), state_p(s)
           END DO
        END DO

        CLOSE(fileID_debug)

     END IF writedebug

     ! clean up:
     DEALLOCATE(U_node_upd,U_elem_upd)

  END IF do_dist

END SUBROUTINE distribute_state_pdaf
