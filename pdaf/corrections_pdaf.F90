module corrections_pdaf

contains

  subroutine correct_state(step, dim_ens, ens_p)

    use mpi
    use assim_pdaf_mod, &
         only: step_null, stdev_SSH_f_p, state_fcst_SSH_p
    use fesom_pdaf, &
         only: mesh_fesom, nlmax, myDim_nod2d, &
         tiny_chl, tiny, chl2N_max, chl2N_max_d, NCmax, &
         NCmax_d, SiCmax, Redfield
    use statevector_pdaf, &
         only: id, nfields, sfields
    use parallel_pdaf_mod, &
         only: mype_filter, comm_filter, MPIerr

    implicit none

    integer, intent(in) :: step
    integer, intent(in) :: dim_ens
    real, intent(inout) :: ens_p(:,:)


    integer :: i, f, k, s
    integer :: member
    integer, allocatable :: count_lim(:,:)
!     integer, allocatable :: count_lim_salt0_g(:) , count_lim_salt0_p(:)       
!     integer, allocatable :: count_lim_ssh_g(:)   , count_lim_ssh_p(:)
!     integer, allocatable :: count_lim_tempM2_g(:), count_lim_tempM2_p(:) 
    real :: diffm                         ! temporary array
    real :: tiny_N, tiny_C               ! Min Phy N, C
    real :: tiny_N_d, tiny_C_d, tiny_Si  ! Min Dia N, C, Si
    real :: tiny_R                       ! Min ZoC


    ! variables allocated and saved during forecast; and deallocated after analysis
    if (.not.allocated(stdev_SSH_f_p)) allocate(stdev_SSH_f_p(sfields(id%SSH)%dim))
    if (.not.allocated(state_fcst_SSH_p)) allocate(state_fcst_SSH_p(sfields(id%SSH)%dim, dim_ens))

  ! allocate correction-counters at initial time; never de-allocated; reset to zero during each analysis
!   allocate(count_lim_salt0_g  (dim_ens))
!   allocate(count_lim_salt0_p  (dim_ens))
!   allocate(count_lim_ssh_g    (dim_ens))
!   allocate(count_lim_ssh_p    (dim_ens))
!   allocate(count_lim_tempM2_g (dim_ens))
!   allocate(count_lim_tempM2_p (dim_ens))

  allocate(count_lim(3, dim_ens))
  count_lim = 0

    Corrections: if ((step-step_null)<0) then

     ! *** store forecast state fields temporarily to compare with analysis afterwards ***    
!     do member = 1, dim_ens
       state_fcst_SSH_p(1 : sfields(id%SSH)%dim, 1:dim_ens) &
            = ens_p(1 + sfields(id%SSH)%off : sfields(id%SSH)%dim + sfields(id%SSH)%off, 1:dim_ens)
!     enddo

  else if ((step-step_null)>0) then Corrections
     ! *** correcting assimilated state fields ***

     ! *** salinity must be > 0 ***
!     count_lim_salt0_p = 0
   
     do member = 1, dim_ens
        do i = 1, sfields(id%salt)%dim
           if (ens_p(i+ sfields(id%salt)%off,member) < 0.0D0) then
              ens_p(i+ sfields(id%salt)%off,member) = 0.0D0
!              count_lim_salt0_p(member) = count_lim_salt0_p(member)+1
              count_lim(1, member) = count_lim(1, member)+1
           end if
        end do
     end do

!      call MPI_Allreduce(count_lim_salt0_p, count_lim_salt0_g, dim_ens, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
!      if (mype_filter == 0) &
!           write(*, *) 'FESOM-PDAF', &
!           '--- Updated salinity limited to zero: ', (count_lim_salt0_g(member), member = 1, dim_ens)
   
     ! *** SSH state update must be <= 2*sigma
!     count_lim_ssh_p = 0

     do member = 1, dim_ens
        do i = 1, sfields(id%SSH)%dim

           diffm = ens_p(sfields(id%SSH)%off+i, member) - state_fcst_SSH_p(i, member)

           if (abs(diffm) > 2.0*stdev_SSH_f_p(i)) then
              ens_p(sfields(id% SSH)%off+i,member) = state_fcst_SSH_p(i,member) + sign(2.0*stdev_SSH_f_p(i),diffm)
!              count_lim_ssh_p(member) = count_lim_ssh_p(member)+1
              count_lim(2, member) = count_lim(2, member)+1
           end if

        end do
     end do
   
!      call MPI_Allreduce(count_lim_ssh_p, count_lim_ssh_g, dim_ens, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
!      if (mype_filter == 0) &
!           write(*,*) 'FESOM-PDAF', &
!           '--- SSH updates limited to 2x standard deviation: ', (count_lim_ssh_g(member), member = 1, dim_ens)

     ! *** temperature must be > -2 degC ***
!     count_lim_tempM2_p = 0
   
     do member = 1, dim_ens
        do i = 1, sfields(id%temp)%dim
           if (ens_p(i+ sfields(id% temp)%off,member) < -2.0) then
              ens_p(i+ sfields(id% temp)%off,member) = -2.0

!              count_lim_tempM2_p(member) = count_lim_tempM2_p(member)+1
              count_lim(3, member) = count_lim(3, member)+1
           end if
        end do
     end do
   
!      call MPI_Allreduce(count_lim_tempM2_p, count_lim_tempM2_g, dim_ens, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
!      if (mype_filter == 0) &
!           write(*,*) 'FESOM-PDAF', &
!           '--- Updated temperature limited to -2 degC: ', (count_lim_tempM2_g(member), member = 1, dim_ens)

     call MPI_Allreduce(MPI_IN_PLACE, count_lim, 4*dim_ens, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
     if (mype_filter == 0) then
        write(*, *) 'FESOM-PDAF', &
             '--- Updated salinity limited to zero: ', (count_lim(1, member), member = 1, dim_ens)
        write(*,*) 'FESOM-PDAF', &
             '--- SSH updates limited to 2x standard deviation: ', (count_lim(2, member), member = 1, dim_ens)
        write(*,*) 'FESOM-PDAF', &
             '--- Updated temperature limited to -2 degC: ', (count_lim(3, member), member = 1, dim_ens)
     end if

     ! *** BGC fields must be larger than "tiny" ***
     if (mype_filter == 0) &
          write(*, *) 'FESOM-PDAF', '--- reset BGC to tiny'
   
     tiny_N   = tiny_chl/chl2N_max      ! Chl2N_max   = 0.00001/ 3.15d0 [mg CHL/mmol N] Maximum CHL-a:N ratio = 0.3 gCHL gN^{-1}
     tiny_N_d = tiny_chl/chl2N_max_d    ! Chl2N_max_d = 0.00001/ 4.2d0
     tiny_C   = tiny_N  /NCmax          ! NCmax       = 0.2d0           [mmol N/mmol C] Maximum cell quota of nitrogen (N:C)
     tiny_C_d = tiny_N_d/NCmax_d        ! NCmax_d     = 0.2d0 
     tiny_Si  = tiny_C_d/SiCmax         ! SiCmax      = 0.8d0
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

!   deallocate(count_lim_salt0_g, count_lim_salt0_p, &
!        count_lim_ssh_g, count_lim_ssh_p, count_lim_tempM2_g, count_lim_tempM2_p)
  deallocate(count_lim)
  
end subroutine correct_state
end module corrections_pdaf
