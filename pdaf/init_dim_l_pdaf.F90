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

  use PDAF, only: PDAFlocal_set_indices, PDAFlocal_set_increment_weights
  use assim_pdaf_mod, &
       only: id_lstate_in_pstate, coords_l
  use coupled_da_mod, &
       only: isweep, type_sweep, cda_bio, cda_phy
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
  integer :: i, ifield                      ! Counters
  integer :: nlay                           ! Number of layers for current domain
  integer :: domain_p                       ! Local analysis domain accounting for multiple sweeps
  logical :: update_cda                     ! Whether to apply DA update
  real, allocatable :: weights_l(:)


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
  
  ! Count local state dimension
  do ifield = 1, nfields

     if (sfields(ifield)%updated) then
        ! surface fields:
        if (sfields(ifield)%ndims == 1) sfields_l(ifield)%dim = 1
        ! 3D fields:
        if (sfields(ifield)%ndims == 2) sfields_l(ifield)%dim = nlay
     else
        ! not updated:
        sfields_l(ifield)%dim = 0
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

  do ifield = 1, nfields
  
     ! only updated fields:
     if ((sfields(ifield)%updated)) then
        if (sfields(ifield)%ndims == 1)   then

           ! surface fields:
           id_lstate_in_pstate (sfields_l(ifield)%off+1) &
                = sfields(ifield)%off + domain_p

        elseif (sfields(ifield)%ndims == 2)   then

           ! 3D fields:
           if (sfields(ifield)%updated) then
              do i = 1, sfields_l(ifield)%dim
                 id_lstate_in_pstate(sfields_l(ifield)%off + i) = &
                      sfields(ifield)%off + (domain_p-1)*(nlmax) + i 
              end do
           endif
        endif
     endif
  enddo

  call PDAFlocal_set_indices(dim_l, id_lstate_in_pstate)


! ****************************************************************************
! *** Initialize array of increment weights for mapping state_l to state_p ***
! ****************************************************************************

  ! Allocate array
  if (allocated(weights_l)) deallocate(weights_l)
  allocate(weights_l(dim_l))
  weights_l(:) = 0.0

  do ifield = 1, nfields

     ! Determine whether to apply update according to coupled data assimilation settings

     if ((sfields(ifield)%updated)) then
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

        if (sfields(ifield)%ndims == 1)   then

           ! surface fields:
           if (update_cda) then
              weights_l(sfields_l(ifield)%off+1) = 1.0
           end if

        elseif (sfields(ifield)%ndims == 2)   then

           ! 3D fields:
           if (update_cda) then
              do i = 1, sfields_l(ifield)%dim
                 weights_l(sfields_l(ifield)%off+i) =  1.0
              end do
           endif
        endif
     end if
  end do

  call PDAFlocal_set_increment_weights(dim_l, weights_l)
  
end subroutine init_dim_l_pdaf
