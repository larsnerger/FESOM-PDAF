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
subroutine init_dim_l_pdaf(step, domain_p_all, dim_l)

  use assim_pdaf_mod, &
       only: id_lstate_in_pstate, coords_l
  use coupled_da_mod, &
       only: isweep
  use parallel_pdaf_mod, &
       only: abort_parallel
  use fesom_pdaf, &
       only: mesh_fesom, nlmax, r2g
  use statevector_pdaf, &
       only: id, nfields, sfields, sfields_l, &
       bgcmin, bgcmax, phymin, phymax
  use fesom_pdaf, &
       only: myDim_nod2D

  implicit none

! *** Arguments ***
  integer, intent(in)  :: step              !< Current time step
  integer, intent(in)  :: domain_p_all      !< Current local analysis domain, containing repititive sweeps
  integer, intent(out) :: dim_l             !< Local state dimension

! *** Local variables ***
  integer :: i, b, id_var                   ! Counters
  integer :: nlay                           ! Number of layers for current domain
  integer :: domain_p                       ! Local analysis domain accounting for multiple sweeps
  

! ********************************************************
! ***  Account for multi sweeps in local analysis loop ***
! ********************************************************

  if (domain_p_all <= myDim_nod2D) then
     domain_p = domain_p_all
     isweep = 1
  else
     domain_p = domain_p_all - myDim_nod2D
     isweep = 2
  end if


! ****************************************
! *** Initialize local state dimension ***
! ****************************************

  ! The local state vector only contains fields that are updated
  
  ! Allocate array
  if (allocated(sfields_l)) deallocate(sfields_l)
  allocate(sfields_l(nfields))

  nlay = mesh_fesom%nlevels_nod2D(domain_p)-1
  
  if (nlay > nlmax) then
     write(*,*) 'FESOM-PDAF ', 'init_dim_l_pdaf ', 'domain_p ', domain_p, ' nlay exceeds layer bounds!'
     call abort_parallel()
  endif
  
  ! Physics:
  do i = phymin, phymax

     if (sfields(i)%updated) then
        ! surface fields:
        if (sfields(i)%ndims == 1) sfields_l(i)%dim = 1
        ! 3D fields:
        if (sfields(i)%ndims == 2) sfields_l(i)%dim = nlay
     else
        ! not updated:
        sfields_l(i)%dim = 0
     endif
  enddo

  ! BGC:
  do i = bgcmin, bgcmax
     if (sfields(i)%updated) then
        ! surface fields:
        if (sfields(i)%ndims == 1) sfields_l(i)%dim = 1
        ! 3D fields:
        if (sfields(i)%ndims == 2) sfields_l(i)%dim = nlay
     else
        ! not updated:
        sfields_l(i)%dim = 0
     endif
  enddo

  ! Set local offsets
  sfields_l(1)%off = 0
  do i = 2, nfields
     sfields_l(i)%off = sfields_l(i-1)%off + sfields_l(i-1)%dim
  end do

  ! *** Local state dimension
  dim_l = sum(sfields_l(:)%dim)


! **********************************************
! *** Initialize coordinates of local domain ***
! **********************************************

  ! Get location of current water column (basis point)
  call r2g(coords_l(1), coords_l(2), mesh_fesom%coord_nod2D(1, domain_p), mesh_fesom%coord_nod2D(2, domain_p))
  

! ****************************************************
! *** Initialize array of indices for local domain ***
! ****************************************************

  ! Allocate array
  if (allocated(id_lstate_in_pstate)) deallocate(id_lstate_in_pstate)
  allocate(id_lstate_in_pstate(dim_l))

  ! *** indices for full state vector ***

  ! SSH
  if (sfields(id%ssh)%updated) then
     id_lstate_in_pstate (sfields_l(id%ssh)%off+1) &
          = sfields(id%ssh)%off + domain_p
  endif
  
  ! U
  id_var = id%u
  if (sfields(id_var)%updated) then
     do i = 1, sfields_l(id_var)%dim
        id_lstate_in_pstate (sfields_l(id_var)%off + i) = &
            sfields(id_var)%off + (domain_p-1)*(nlmax) + i 
     end do
  endif
  
  ! V
  id_var = id%v
  if (sfields(id_var)%updated) then
     do i = 1, sfields_l(id_var)%dim
        id_lstate_in_pstate (sfields_l(id_var)%off + i) = &
            sfields(id_var)%off + (domain_p-1)*(nlmax) + i 
     end do
  endif
        
  ! W
  ! id_lstate_in_pstate (sfields_l(id%w)%off+1 : sfields_l(id%w+1)%off) &
  !      = sfields(id%w)%off &
  !      + (domain_p-1)*(nlmax) &
  !      + (/(i, i=1,sfields_l(id%w)%dim)/)
  
  ! Temp
  id_var = id%temp
  if (sfields(id_var)%updated) then
     do i = 1, sfields_l(id_var)%dim
        id_lstate_in_pstate (sfields_l(id_var)%off + i) = &
            sfields(id_var)%off + (domain_p-1)*(nlmax) + i 
     end do
  endif
  
  ! Salt
  id_var = id%salt
  if (sfields(id_var)%updated) then
     do i = 1, sfields_l(id_var)%dim
        id_lstate_in_pstate (sfields_l(id_var)%off + i) = &
            sfields(id_var)%off + (domain_p-1)*(nlmax) + i 
     end do
  endif
        
  ! BGC:
  do b = bgcmin, bgcmax
  
     ! only updated fields:
     if ((sfields(b)%updated)) then
        if (sfields(b)%ndims == 1)   then

           ! surface fields:
           id_lstate_in_pstate (sfields_l(b)%off+1) &
                = sfields(b)%off + domain_p

        elseif (sfields(b)%ndims == 2)   then

           ! 3D fields:
           id_var = b
           if (sfields(id_var)%updated) then
              do i = 1, sfields_l(id_var)%dim
                 id_lstate_in_pstate (sfields_l(id_var)%off + i) = &
                      sfields(id_var)%off + (domain_p-1)*(nlmax) + i 
              end do
           endif
        endif
     endif
  enddo

end subroutine init_dim_l_pdaf
