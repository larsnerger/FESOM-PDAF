!> Timing and clean-up of PDAF
!!
!! The routine prints timing and memory information
!! and deallocates PDAF internal arrays.
!!
module finalize_pdaf_mod

contains
  subroutine finalize_pdaf()

    use mpi
    use PDAF, &
         only: PDAF_print_info, PDAF_deallocate
    use parallel_pdaf_mod, &
         only: mype_ens, npes_ens, comm_ensemble, mpierr
    use timer, &
         only: timeit, time_tot

    call timeit(5, 'old')
    call timeit(1, 'old')

    ! Show allocated memory for PDAF
    if (mype_ens==0) call PDAF_print_info(10)
    if (npes_ens>1) call PDAF_print_info(11)

    ! Print PDAF timings onto screen
    if (mype_ens==0) call PDAF_print_info(3)

    ! Deallocate PDAF arrays
    call PDAF_deallocate()

    if (mype_ens==0) then
       write (*, '(/a,10x,a)') 'NEMO-PDAF', 'Model-sided timing overview'
       write (*, '(a,2x,a)') 'NEMO-PDAF', '-----------------------------------'
       write (*, '(a,8x,a,F11.3,1x,a)') 'NEMO-PDAF', 'initialize MPI:  ', time_tot(2), 's'
       write (*, '(a,8x,a,F11.3,1x,a)') 'NEMO-PDAF', 'initialize model:', time_tot(3), 's'
       write (*, '(a,8x,a,F11.3,1x,a)') 'NEMO-PDAF', 'initialize PDAF :', time_tot(4), 's'
       write (*, '(a,8x,a,F11.3,1x,a)') 'NEMO-PDAF', 'main part:       ', time_tot(5), 's'
       write (*, '(a,8x,a,F11.3,1x,a)') 'NEMO-PDAF', 'PDAF analysis:    ', time_tot(6), 's'
       write (*, '(a,8x,a,F11.3,1x,a)') 'NEMO-PDAF', 'means in assimilate_pdaf:    ', time_tot(7), 's'
       write (*, '(a,8x,a,F11.3,1x,a)') 'NEMO-PDAF', 'total:         ', time_tot(1), 's'
    end if

    call mpi_barrier(comm_ensemble, mpierr)

  end subroutine finalize_pdaf

end module finalize_pdaf_mod
