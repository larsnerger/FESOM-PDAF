!> Module providing some functionality for coupled DA
!!
!! This module provides routines for weakly and strongly
!! coupled DA, both for the coupled ocean-atmosphere case
!! of AWI-CM/AWI-ESM and for coupled physics-biogeochemistry
!! DA (FESOM-REcoM).
!!
!! __Revision history:__
!! * 2025-12 - Lars Nerger - initial code from restructuring AWI-CM DA code
!!
module coupled_da_mod

  implicit none
  public
  save

  ! For weakly or strongly coupled DA determine by parallelization (e.g .ocean & atmosphere)
  integer :: DA_couple_type=0         ! (0) for weakly-coupled, (1) for strongly-coupled assimilation

  ! For coupled physics-BGC DA (set by user)
  character(len=6) :: cda_phy='weak'  !< Flag whether strongly-coupled DA is done for physics data
  character(len=6) :: cda_bio='weak'  !< Flag whether strongly-coupled DA is done for bgc data

  ! Internally set variables
  logical :: assimilateBGC            !< whether BGC observations are assimilated (set by check_coupled_da)
  logical :: assimilatePHY            !< whether physics observations are assimilated (set by check_coupled_da)
  integer :: n_sweeps                 !< Number of sweeps in local analysis loop
  character(len=3) :: type_sweep(2)   !< Type of sweep in local analysis loop
  integer :: isweep                   !< Index of sweep during the local analysis loop

contains

!-----------------------------------------------------------
!> Define local analysis sweeps for physics-BGC DA
!!
!! this routine determines how many sweeps over the model grid 
!! have to be performed for the local analysis of 
!! weakly-coupled physics-BGC assimilation.
!! 
  subroutine cda_set_sweeps()

    use parallel_pdaf_mod, &
         only: mype_world

    implicit none

    integer :: i                 ! Counter
    character(len=6) :: cdaval   ! Flag whether strongly-coupled DA is done

    ! Set flags for assimilating physics and BGC observations
    call check_coupled_da_pdafomi(assimilatePHY, assimilateBGC)

    if (assimilateBGC .and. assimilatePHY) then
       ! Observations of both physics and BGC are assimilated
       n_sweeps = 2
       type_sweep(1) = 'phy'
       type_sweep(2) = 'bio'
    else
       ! Less than two observation categories
       n_sweeps = 1
       if (assimilatePHY) then
          ! Only observations of physics are assimilated
          type_sweep(1) = 'phy'
       elseif (assimilateBGC) then
          ! Only observations of BGC are assimilated
          type_sweep(1) = 'bio'
       else
          ! No observation active (free run); set sweep to physics
          type_sweep(1) = 'phy'
       end if
    end if

    if (mype_world == 0) then
       write (*,'(a,2x,a)') 'FESOM-PDAF', '*** Setup for coupled DA FESOM-REcoM ***'
       write (*, '(a,4x,a,i5)') 'FESOM-PDAF', 'Number of local analysis sweeps', n_sweeps
       write (*, '(a,4x,a)') 'FESOM-PDAF','Type of sweeps:'
       do i = 1, n_sweeps
          if (trim(type_sweep(i))=='phy') then
             cdaval = cda_phy
          else
             cdaval = cda_bio
          end if
          write (*, '(a,8x,a,i3,3x,a,a,3x,a,a)') &
               'FESOM-PDAF', 'sweep', i, ' observation type: ', trim(type_sweep(i)), 'CDA: ', trim(cdaval)
       end do
    end if

  end subroutine cda_set_sweeps

!-----------------------------------------------------------
!> Set FESOM communicator for weakly-coupled ocean-atmosphere DA
!!
!! This routine is usd for the coupled DA in the coupled
!! ocean-atmosphere model setup (e.g. FESOM-oIFS) in which
!! the ocean and atmosphere are executed on distinct sets
!! of MPI tasks. For weakly-coupled DA, the filter communicator
!! in FESOM then only contains the processes in the ocean
!!
  subroutine cda_reset_filter_comm(comm)

    use MPI
    use parallel_pdaf_mod, &
         only: COMM_filter, npes_filter, mype_filter, &
         MPIerr, filterpe

    implicit none

    integer, intent(in) :: comm      ! Input communicator

    if (DA_couple_type == 0) then
       ! Set filter communicator to the communicator of FESOM
       COMM_filter = comm

       if (filterpe) then
          call MPI_Comm_Size(COMM_filter, npes_filter, MPIerr)
          call MPI_Comm_Rank(COMM_filter, mype_filter, MPIerr)

          if (mype_filter==0) then
             write (*,'(a)') 'FESOM-Atmos-PDAF: Initialize weakly-coupled data assimilation'
          endif
       endif
    else
       if (filterpe) then
          if (mype_filter==0) then
             write (*,'(a)') 'FESOM-Atmpos-PDAF: Initialize strongly-coupled data assimilation'
          end if
       end if
    end if

  end subroutine cda_reset_filter_comm

end module coupled_da_mod
