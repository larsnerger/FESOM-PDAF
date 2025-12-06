!> Collecting the state vector variables from model fields
!!
!! This subroutine is called during the forecast 
!! phase propagation of each ensemble member. 
!! The supplied state vector has to be initialized
!! from the model fields (typically via a module). 
!! With parallelization, MPI communication might be 
!! required to initialize state vectors for all 
!! subdomains on the model PEs. 
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! __Revision history:__
!! * 2004-11 - Lars Nerger - Initial code
!! * 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
!! * 2022-03 - Frauke B     - Removed sea-ice, added velocities
!! * 2025-12 - Lars Nerger - Revision for PDAF3
!!
SUBROUTINE collect_state_pdaf(dim_p, state_p)

  USE mod_parallel_pdaf, &
       ONLY: mype_world
  USE statevector_pdaf, &
       ONLY: id, sfields, nfields !, ids_tr3D, nfields_tr3D
  USE fesom_pdaf, &
       ONLY: nlmax, daynew, timenew, mydim_nod2d, &
       eta_n, uv, wvel, tr_arr, unode, MLD1, MLD2, sigma0, &
       GloPCO2surf, GloCO2flux, Diags3D, PAR3D, export, &
       PistonVelocity, alphaCO2, a_ice
  USE PDAF, &
       ONLY: PDAF_get_assim_flag

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p            !< PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)   !< Local state vector

! *** Local variables
  INTEGER :: i, k, b, s, istate, ifesom   !< Counters
  
! Debugging:
  LOGICAL            :: debugmode
  INTEGER            :: fileID_debug
  CHARACTER(len=5)   :: mype_string, tim_string
  CHARACTER(len=3)   :: day_string
  INTEGER            :: assim_flag

  ! Set debug output
  debugmode = .false.


! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

   ! *** Dimensions of FESOM arrays:
   ! * eta_n          (myDim_nod2D + eDim_nod2D)            ! Sea Surface Height
   ! * UV(1,:,:)      (1, nl-1, myDim_elem2D + eDim_elem2D) ! Velocity u
   ! * UV(2,:,:)      (1, nl-1, myDim_elem2D + eDim_elem2D) ! Velocity v
   ! * wvel           (nl,   myDim_nod2D + eDim_nod2D)      ! Velocity w
   ! * tr_arr(:,:,1)  (nl-1, myDim_nod2D + eDim_nod2D, 1)   ! Temperature
   ! * tr_arr(:,:,2)  (nl-1, myDim_nod2D + eDim_nod2D, 1)   ! Salinity
   ! * unode(1,:,:)   (1, nl-1, myDim_elem2D + eDim_elem2D) ! Velocity u interpolated on nodes
   ! * unode(2,:,:)   (1, nl-1, myDim_elem2D + eDim_elem2D) ! Velocity v interpolated on nodes
   ! * a_ice          (myDim_nod2D + eDim_nod2D)            ! Sea-ice concentration
   ! ***

  ! SSH
  DO i = 1, myDim_nod2D
     state_p(i + sfields(id%SSH)%off) = eta_n(i)
  END DO

  ! u and v velocities (interpolated on nodes)
   
  ! U
  s = sfields(id%u)%off
  DO i = 1, myDim_nod2D
     DO k =1, nlmax
        s = s + 1
        state_p(s) = Unode(1, k, i)
     END DO
  END DO

  ! V
  s = sfields(id%v)%off
  DO i = 1, myDim_nod2D
     DO k =1, nlmax
        s = s + 1
        state_p(s) = Unode(2, k, i)
     END DO
  END DO

  ! w velocity
  s = sfields(id%w)%off
  DO i = 1, myDim_nod2D
     DO k = 1, nlmax
        s = s + 1
        state_p(s) = wvel(k, i)
     END DO
  END DO
   
  ! temp and salt are model 3D tracers.
   
  ! sea-ice concentration
  s = sfields(id%a_ice)%off
  DO i = 1, myDim_nod2D
     s = s + 1
     state_p(s) = a_ice(i)
  END DO

  ! MLD1
  DO i = 1, myDim_nod2D
     state_p(i + sfields(id%MLD1)%off) = MLD1(i)
  END DO

  ! MLD2
  DO i = 1, myDim_nod2D
     state_p(i + sfields(id%MLD2)%off) = MLD2(i)
  END DO
   
  ! potential density
  s = sfields(id%sigma)%off
  DO i = 1, myDim_nod2D
     DO k = 1, nlmax
        s = s + 1
        state_p(s) = 1e-3 * sigma0(k, i)
     ENDDO
  ENDDO

  ! model 3D tracers
  DO istate = 1, nfields
!  DO b = 1, nfields_tr3D
           
     IF (sfields(istate)%trnumfesom>0) THEN
!     istate = ids_tr3D(b)                 ! index of field in state vector
        ifesom = sfields(istate)%trnumfesom  ! index of field in model tracer array

        s = 0
        DO i = 1, myDim_nod2D
           DO k = 1, nlmax
              s = s + 1
              state_p(s + sfields(istate)%off) = tr_arr(k, i, ifesom)
           ENDDO
        ENDDO
     END IF
  ENDDO

  ! biogeochem 2D fields
  DO i = 1, myDim_nod2D
     state_p(i + sfields(id% pCO2s)%off)    = GloPCO2surf(i)    ! surface ocean partial pressure CO2
  ENDDO

  DO i = 1, myDim_nod2D
     state_p(i + sfields(id% CO2f)%off)     = GloCO2flux(i)     ! CO2 flux (from atmosphere into ocean)
  ENDDO

  DO i = 1, myDim_nod2D
     state_p(i + sfields(id% export)%off)    = export(i)         ! Export through particle sinking
  ENDDO

  DO i = 1, myDim_nod2D
     state_p(i + sfields(id% alphaCO2)%off)  = alphaCO2(i)       ! Solubility of CO2
  ENDDO

  DO i = 1, myDim_nod2D
     state_p(i + sfields(id% PistonVel)%off) = PistonVelocity(i) ! Air-sea gas transfer velocity
  ENDDO

  ! diagnostic biogeochemical 3D fields
  s = 0
  DO i = 1, myDim_nod2D
     DO k = 1, nlmax        
        s = s + 1
        state_p(s + sfields(id% PAR)%off)    = PAR3D(k, i)      ! photosynthetically active radiation
     ENDDO
  ENDDO

  s = 0
  DO i = 1, myDim_nod2D
     DO k = 1, nlmax        
        s = s + 1
        state_p(s + sfields(id% NPPn)%off)   = diags3D(k, i, 1) ! net primary production small phytoplankton
     ENDDO
  ENDDO

  s = 0
  DO i = 1, myDim_nod2D
     DO k = 1, nlmax        
        s = s + 1
        state_p(s + sfields(id% NPPd)%off)   = diags3D(k, i, 2) ! net primary production diatoms
     ENDDO
  ENDDO


! ********************
! *** Debug output ***
! ********************

  ! Check whether analysis step was done
  CALL PDAF_get_assim_flag(assim_flag)

  writedebug:  IF (debugmode .and. assim_flag>0 .and. mype_world==0) THEN

     ! print state vector
     WRITE(day_string, '(i3.3)') daynew
     WRITE(tim_string, '(i5.5)') int(timenew)
     fileID_debug=10
     open(unit=fileID_debug, file='collect_state_pdaf_'//day_string//'_'//tim_string//'.txt', status='unknown')

     ! SSH
     DO i = 1, myDim_nod2D
        WRITE(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%SSH)%variable, i+sfields(id%SSH)%off, eta_n(i)
     END DO

     ! U
     s = sfields(id%u)%off
     DO i = 1, myDim_nod2D
        DO k =1, nlmax
           s = s + 1
           WRITE(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%u)%variable, s, Unode(1, k, i)
        END DO
     END DO
   
     ! V
     s = sfields(id%v)%off
     DO i = 1, myDim_nod2D
        DO k =1, nlmax
           s = s + 1
           WRITE(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%v)%variable, s, Unode(2, k, i)
        END DO
     END DO
 
     ! w velocity
     s = sfields(id%w)%off
     DO i = 1, myDim_nod2D
        DO k = 1, nlmax
           s = s + 1
           WRITE(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%w)%variable, s, wvel(k, i)
        END DO
     END DO

     ! temp and salt are model 3D tracers.
   
     ! sea-ice concentration
     s = sfields(id%a_ice)%off
     DO i = 1, myDim_nod2D
        s = s + 1
        WRITE(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%a_ice)%variable, s, a_ice(i)
     END DO

     ! MLD1
     DO i = 1, myDim_nod2D
        WRITE(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%MLD1)%variable, i+sfields(id%MLD1)%off, MLD1(i)
     END DO

     ! MLD2
     DO i = 1, myDim_nod2D
        WRITE(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%MLD2)%variable, i+sfields(id%MLD2)%off, MLD2(i)
     END DO
      
     ! biogeochem 2D fields
     DO i = 1, myDim_nod2D   
        write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%pCO2s    )%variable, i + sfields(id% pCO2s)%off   , GloPCO2surf(i)    
        write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%CO2f     )%variable, i + sfields(id% CO2f)%off    , GloCO2flux(i)     
        write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%export   )%variable, i + sfields(id% export)%off   , export(i)         
        write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%alphaCO2 )%variable, i + sfields(id% alphaCO2)%off , alphaCO2(i)       
        write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%PistonVel)%variable, i + sfields(id% PistonVel)%off, PistonVelocity(i) 
     ENDDO

     ! diagnostic biogeochemical 3D fields
     s = 0
     DO i = 1, myDim_nod2D
        DO k = 1, nlmax
           s = s + 1
           write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%PAR  )%variable, s + sfields(id% PAR)%off , PAR3D(k, i)     
           write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%NPPn )%variable, s + sfields(id% NPPn)%off, diags3D(k, i, 1)
           write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%NPPd )%variable, s + sfields(id% NPPd)%off, diags3D(k, i, 2)
        ENDDO
     ENDDO

     ! Close debug-file
     close(fileID_debug)

  END IF writedebug
  
END SUBROUTINE collect_state_pdaf
