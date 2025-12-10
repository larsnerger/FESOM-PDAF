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
subroutine collect_state_pdaf(dim_p, state_p)

  use PDAF, &
       only: PDAF_get_assim_flag
  use parallel_pdaf_mod, &
       only: mype_world
  use statevector_pdaf, &
       only: id, sfields, nfields
  use fesom_pdaf, &
       only: nlmax, daynew, timenew, mydim_nod2d, &
       eta_n, uv, wvel, tr_arr, unode, MLD1, MLD2, sigma0, &
       GloPCO2surf, GloCO2flux, Diags3D, PAR3D, export, &
       PistonVelocity, alphaCO2, a_ice

  implicit none
  
! *** Arguments ***
  integer, intent(in) :: dim_p            !< PE-local state dimension
  real, intent(inout) :: state_p(dim_p)   !< Local state vector

! *** Local variables
  integer :: i, k, b, s, istate, ifesom   !< Counters
  
! Debugging:
  logical            :: debugmode
  integer            :: fileID_debug
  character(len=5)   :: mype_string, tim_string
  character(len=3)   :: day_string
  integer            :: assim_flag

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
  do i = 1, myDim_nod2D
     state_p(i + sfields(id%SSH)%off) = eta_n(i)
  end do

  ! u and v velocities (interpolated on nodes)
   
  ! U
  s = sfields(id%u)%off
  do i = 1, myDim_nod2D
     do k =1, nlmax
        s = s + 1
        state_p(s) = Unode(1, k, i)
     end do
  end do

  ! V
  s = sfields(id%v)%off
  do i = 1, myDim_nod2D
     do k =1, nlmax
        s = s + 1
        state_p(s) = Unode(2, k, i)
     end do
  end do

  ! w velocity
  s = sfields(id%w)%off
  do i = 1, myDim_nod2D
     do k = 1, nlmax
        s = s + 1
        state_p(s) = wvel(k, i)
     end do
  end do
   
  ! temp and salt are model 3D tracers.
   
  ! sea-ice concentration
  s = sfields(id%a_ice)%off
  do i = 1, myDim_nod2D
     s = s + 1
     state_p(s) = a_ice(i)
  end do

  ! MLD1
  do i = 1, myDim_nod2D
     state_p(i + sfields(id%MLD1)%off) = MLD1(i)
  end do

  ! MLD2
  do i = 1, myDim_nod2D
     state_p(i + sfields(id%MLD2)%off) = MLD2(i)
  end do
   
  ! potential density
  s = sfields(id%sigma)%off
  do i = 1, myDim_nod2D
     do k = 1, nlmax
        s = s + 1
        state_p(s) = 1e-3 * sigma0(k, i)
     enddo
  enddo

  ! model 3D tracers
  do istate = 1, nfields
     if (sfields(istate)%trnumfesom>0) then

        ifesom = sfields(istate)%trnumfesom  ! index of field in model tracer array

        s = 0
        do i = 1, myDim_nod2D
           do k = 1, nlmax
              s = s + 1
              state_p(s + sfields(istate)%off) = tr_arr(k, i, ifesom)
           enddo
        enddo
     end if
  enddo

  ! biogeochem 2D fields
  do i = 1, myDim_nod2D
     state_p(i + sfields(id% pCO2s)%off)    = GloPCO2surf(i)    ! surface ocean partial pressure CO2
  enddo

  do i = 1, myDim_nod2D
     state_p(i + sfields(id% CO2f)%off)     = GloCO2flux(i)     ! CO2 flux (from atmosphere into ocean)
  enddo

  do i = 1, myDim_nod2D
     state_p(i + sfields(id% export)%off)    = export(i)         ! Export through particle sinking
  enddo

  do i = 1, myDim_nod2D
     state_p(i + sfields(id% alphaCO2)%off)  = alphaCO2(i)       ! Solubility of CO2
  enddo

  do i = 1, myDim_nod2D
     state_p(i + sfields(id% PistonVel)%off) = PistonVelocity(i) ! Air-sea gas transfer velocity
  enddo

  ! diagnostic biogeochemical 3D fields
  s = 0
  do i = 1, myDim_nod2D
     do k = 1, nlmax        
        s = s + 1
        state_p(s + sfields(id% PAR)%off)    = PAR3D(k, i)      ! photosynthetically active radiation
     enddo
  enddo

  s = 0
  do i = 1, myDim_nod2D
     do k = 1, nlmax        
        s = s + 1
        state_p(s + sfields(id% NPPn)%off)   = diags3D(k, i, 1) ! net primary production small phytoplankton
     enddo
  enddo

  s = 0
  do i = 1, myDim_nod2D
     do k = 1, nlmax        
        s = s + 1
        state_p(s + sfields(id% NPPd)%off)   = diags3D(k, i, 2) ! net primary production diatoms
     enddo
  enddo


! ********************
! *** Debug output ***
! ********************

  ! Check whether analysis step was done
  call PDAF_get_assim_flag(assim_flag)

  writedebug:  if (debugmode .and. assim_flag>0 .and. mype_world==0) then

     ! print state vector
     write(day_string, '(i3.3)') daynew
     write(tim_string, '(i5.5)') int(timenew)
     fileID_debug=10
     open(unit=fileID_debug, file='collect_state_pdaf_'//day_string//'_'//tim_string//'.txt', status='unknown')

     ! SSH
     do i = 1, myDim_nod2D
        write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%SSH)%variable, i+sfields(id%SSH)%off, eta_n(i)
     end do

     ! U
     s = sfields(id%u)%off
     do i = 1, myDim_nod2D
        do k =1, nlmax
           s = s + 1
           write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%u)%variable, s, Unode(1, k, i)
        end do
     end do
   
     ! V
     s = sfields(id%v)%off
     do i = 1, myDim_nod2D
        do k =1, nlmax
           s = s + 1
           write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%v)%variable, s, Unode(2, k, i)
        end do
     end do
 
     ! w velocity
     s = sfields(id%w)%off
     do i = 1, myDim_nod2D
        do k = 1, nlmax
           s = s + 1
           write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%w)%variable, s, wvel(k, i)
        end do
     end do

     ! temp and salt are model 3D tracers.
   
     ! sea-ice concentration
     s = sfields(id%a_ice)%off
     do i = 1, myDim_nod2D
        s = s + 1
        write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%a_ice)%variable, s, a_ice(i)
     end do

     ! MLD1
     do i = 1, myDim_nod2D
        write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%MLD1)%variable, i+sfields(id%MLD1)%off, MLD1(i)
     end do

     ! MLD2
     do i = 1, myDim_nod2D
        write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%MLD2)%variable, i+sfields(id%MLD2)%off, MLD2(i)
     end do
      
     ! biogeochem 2D fields
     do i = 1, myDim_nod2D   
        write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%pCO2s    )%variable, i + sfields(id% pCO2s)%off   , GloPCO2surf(i)    
        write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%CO2f     )%variable, i + sfields(id% CO2f)%off    , GloCO2flux(i)     
        write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%export   )%variable, i + sfields(id% export)%off   , export(i)         
        write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%alphaCO2 )%variable, i + sfields(id% alphaCO2)%off , alphaCO2(i)       
        write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%PistonVel)%variable, i + sfields(id% PistonVel)%off, PistonVelocity(i) 
     enddo

     ! diagnostic biogeochemical 3D fields
     s = 0
     do i = 1, myDim_nod2D
        do k = 1, nlmax
           s = s + 1
           write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%PAR  )%variable, s + sfields(id% PAR)%off , PAR3D(k, i)     
           write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%NPPn )%variable, s + sfields(id% NPPn)%off, diags3D(k, i, 1)
           write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%NPPd )%variable, s + sfields(id% NPPd)%off, diags3D(k, i, 2)
        enddo
     enddo

     ! Close debug-file
     close(fileID_debug)

  end if writedebug
  
end subroutine collect_state_pdaf
