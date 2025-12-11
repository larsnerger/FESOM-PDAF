!> Check validity of state vector values and correct if necessary
!!
!! To ensure the consistency of the ensemble states after the
!! analysis update their values are checked here and corrected if required.
!!
!! __Revision history__
!! * 2025-12 - Lars Nerger - initial version extracing code by Frauke from prepoststep_pdaf
!!
module corrections_pdaf

  implicit none
  public
  save

  real, allocatable :: stddev_SSH_f_p(:)     ! forecast ensemble standard deviation at grid points for
  real, allocatable :: state_fcst_SSH_p(:,:) ! state prior to assimilation, saved to use for correction

contains

  subroutine correct_state(step, dim_ens, ens_p)

    use mpi
    use assim_pdaf_mod, &
         only: step_null
    use fesom_pdaf, &
         only: mesh_fesom, nlmax, myDim_nod2d
    use statevector_pdaf, &
         only: id, nfields, sfields
    use parallel_pdaf_mod, &
         only: mype_filter, comm_filter, MPIerr

    implicit none

    integer, intent(in) :: step             !< Current time step, starting from 0 at beginning of year
    integer, intent(in) :: dim_ens          !< Size of state ensemble
    real, intent(inout) :: ens_p(:,:)       !< PE-local ensemble


    integer :: i, member                    ! Counters
    integer, allocatable :: count_lim(:,:)  ! Count corrected state values
                                            ! first index: 1 salt, 2 ssh, 3 tempM2, 2nd index: member
    real :: diffm                           ! temporary array

    ! variables allocated and saved during forecast; and deallocated after analysis
    if (.not.allocated(stddev_SSH_f_p)) allocate(stddev_SSH_f_p(sfields(id%SSH)%dim))
    if (.not.allocated(state_fcst_SSH_p)) allocate(state_fcst_SSH_p(sfields(id%SSH)%dim, dim_ens))


    Corrections: if ((step-step_null)<0) then

       ! *** store forecast state fields temporarily to compare with analysis afterwards ***    
       state_fcst_SSH_p(1 : sfields(id%SSH)%dim, 1:dim_ens) &
            = ens_p(1 + sfields(id%SSH)%off : sfields(id%SSH)%dim + sfields(id%SSH)%off, 1:dim_ens)

    else if ((step-step_null)>0) then Corrections

       allocate(count_lim(3, dim_ens))
       count_lim = 0

       ! *** correcting assimilated state fields ***

       ! *** salinity must be > 0 ***
       do member = 1, dim_ens
          do i = 1, sfields(id%salt)%dim
             if (ens_p(i+ sfields(id%salt)%off,member) < 0.0D0) then
                ens_p(i+ sfields(id%salt)%off,member) = 0.0D0
                count_lim(1, member) = count_lim(1, member)+1
             end if
          end do
       end do
   
       ! *** SSH state update must be <= 2*sigma

       do member = 1, dim_ens
          do i = 1, sfields(id%SSH)%dim

             diffm = ens_p(sfields(id%SSH)%off+i, member) - state_fcst_SSH_p(i, member)

             if (abs(diffm) > 2.0*stddev_SSH_f_p(i)) then
                ens_p(sfields(id% SSH)%off+i,member) = state_fcst_SSH_p(i,member) + sign(2.0*stddev_SSH_f_p(i),diffm)
                count_lim(2, member) = count_lim(2, member)+1
             end if
          end do
       end do

       ! *** temperature must be > -2 degC ***
   
       do member = 1, dim_ens
          do i = 1, sfields(id%temp)%dim
             if (ens_p(i+ sfields(id% temp)%off,member) < -2.0) then
                ens_p(i+ sfields(id% temp)%off,member) = -2.0
                count_lim(3, member) = count_lim(3, member)+1
             end if
          end do
       end do

       ! *** Get global counts and write them to screen
       call MPI_Allreduce(MPI_IN_PLACE, count_lim, 4*dim_ens, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)

       if (mype_filter == 0) then
          write(*, *) 'FESOM-PDAF', &
               '--- Updated salinity limited to zero: ', (count_lim(1, member), member = 1, dim_ens)
          write(*,*) 'FESOM-PDAF', &
               '--- SSH updates limited to 2x standard deviation: ', (count_lim(2, member), member = 1, dim_ens)
          write(*,*) 'FESOM-PDAF', &
               '--- Updated temperature limited to -2 degC: ', (count_lim(3, member), member = 1, dim_ens)
       end if

       ! Deallocate arrays
       if (allocated(stddev_SSH_f_p))    deallocate(stddev_SSH_f_p)
       if (allocated(state_fcst_SSH_p)) deallocate(state_fcst_SSH_p)
       deallocate(count_lim)

       ! *** Check recom fields ***

       call correct_state_recom(step, dim_ens, ens_p)

    end if Corrections

  
  end subroutine correct_state

!-------------------------------------------------------------
!> Store SSH part of standard deviation field
!!
!! At forecast time this routine stores the ensemble standard
!! deviation field computed in prepoststep_pdaf kto use it for
!! checks in correct_state().
!!
  subroutine store_stddev(step, stddev_p)

    use assim_pdaf_mod, &
         only: step_null
    use statevector_pdaf, &
         only: id, sfields

    implicit none

    integer, intent(in) :: step             !< Current time step, starting from 0 at beginning of year
    real, intent(inout) :: stddev_p(:)      !< PE-local standard deviation field

    if ((step-step_null) < 0) then
       stddev_SSH_f_p = stddev_p( sfields(id%SSH)%off+1 : sfields(id%SSH)%off+sfields(id%SSH)%dim)
    endif

  end subroutine store_stddev


!-------------------------------------------------------------
!> Check REcoM fields for consistency and correct if required
!!
!! This routine check the values of the REcoM fields to
!! ensure that their values are not too low.
!!
  subroutine correct_state_recom(step, dim_ens, ens_p)

    use assim_pdaf_mod, &
         only: step_null
    use fesom_pdaf, &
         only: mesh_fesom, nlmax, myDim_nod2d, &
         tiny_chl, tiny, chl2N_max, chl2N_max_d, NCmax, &
         NCmax_d, SiCmax, Redfield
    use statevector_pdaf, &
         only: id, nfields, sfields
    use parallel_pdaf_mod, &
         only: mype_filter

    implicit none

    integer, intent(in) :: step             !< Current time step, starting from 0 at beginning of year
    integer, intent(in) :: dim_ens          !< Size of state ensemble
    real, intent(inout) :: ens_p(:,:)       !< PE-local ensemble


    integer :: i, f, k, s, member           ! Counters
    real :: diffm                           ! temporary array
    real :: tiny_N, tiny_C                  ! Min Phy N, C
    real :: tiny_N_d, tiny_C_d, tiny_Si     ! Min Dia N, C, Si
    real :: tiny_R                          ! Min ZoC


    Corrections: if ((step-step_null)>0) then

       ! *** BGC fields must be larger than "tiny" ***
       if (mype_filter == 0) write(*, *) 'FESOM-PDAF', '--- reset BGC to tiny'
   
       tiny_N   = tiny_chl/chl2N_max      ! Chl2N_max   = 0.00001/ 3.15 [mg CHL/mmol N] Maximum CHL-a:N ratio = 0.3 gCHL gN^{-1}
       tiny_N_d = tiny_chl/chl2N_max_d    ! Chl2N_max_d = 0.00001/ 4.2
       tiny_C   = tiny_N  /NCmax          ! NCmax       = 0.2           [mmol N/mmol C] Maximum cell quota of nitrogen (N:C)
       tiny_C_d = tiny_N_d/NCmax_d        ! NCmax_d     = 0.2 
       tiny_Si  = tiny_C_d/SiCmax         ! SiCmax      = 0.8
       tiny_R   = tiny * Redfield
   
       do member = 1, dim_ens
          do i = 1, myDim_nod2D
             do k = 1, mesh_fesom%nlevels_nod2D(i)-1 ! loop through all wet nodes

                s = (i-1) * (nlmax) + k ! index in state vector

                do f=1,nfields
                   ! biogeochemical model tracers
                   if ((sfields(f)%bgc) .and. (sfields(f)% trnumfesom > 0)) then
                      ens_p(sfields(f)%off+s,member) = max(tiny,ens_p(sfields(f)%off+s,member))
                   endif
                enddo

                ! small phytoplankton
                ens_p(s+ sfields(id% PhyN   )%off,member) = max(tiny_N  ,ens_p(s+ sfields(id% PhyN   )%off,member))
                ens_p(s+ sfields(id% PhyC   )%off,member) = max(tiny_C  ,ens_p(s+ sfields(id% PhyC   )%off,member))
                ens_p(s+ sfields(id% PhyChl )%off,member) = max(tiny_chl,ens_p(s+ sfields(id% PhyChl )%off,member))
                ens_p(s+ sfields(id% PhyCalc)%off,member) = max(tiny    ,ens_p(s+ sfields(id% PhyCalc)%off,member))

                ! diatoms
                ens_p(s+ sfields(id% DiaN   )%off,member) = max(tiny_N_d,ens_p(s+ sfields(id% DiaN   )%off,member))
                ens_p(s+ sfields(id% DiaC   )%off,member) = max(tiny_C_d,ens_p(s+ sfields(id% DiaC   )%off,member))
                ens_p(s+ sfields(id% DiaChl )%off,member) = max(tiny_chl,ens_p(s+ sfields(id% DiaChl )%off,member))
                ens_p(s+ sfields(id% DiaSi  )%off,member) = max(tiny_Si ,ens_p(s+ sfields(id% DiaSi  )%off,member))

                ! zooplankton 1
                ens_p(s+ sfields(id% Zo1N   )%off,member) = max(tiny    ,ens_p(s+ sfields(id% Zo1N   )%off,member))
                ens_p(s+ sfields(id% Zo1C   )%off,member) = max(tiny_R  ,ens_p(s+ sfields(id% Zo1C   )%off,member))

                ! zooplankton 2
                ens_p(s+ sfields(id% Zo2N   )%off,member) = max(tiny    ,ens_p(s+ sfields(id% Zo2N   )%off,member))
                ens_p(s+ sfields(id% Zo2C   )%off,member) = max(tiny_R  ,ens_p(s+ sfields(id% Zo2C   )%off,member))

             enddo ! k=1,nlmax
          enddo ! i=1,my_Dim_nod2D
       enddo ! member=1,dim_ens

    end if Corrections
  
  end subroutine correct_state_recom

end module corrections_pdaf
