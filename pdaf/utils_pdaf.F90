module utils_pdaf
contains

!> Routine to decide about output timing
!!
!! This routine decides whether it's time to do output
!! during analysis and forecast phase. This is analogous
!! to routines in gen_events.F90 of FESOM
!!
  subroutine monthly_event_assimstep(do_output, N_opt)

    use g_clock
    use assim_pdaf_mod, only: assim_time

    implicit none

! *** Arguments ***
    logical :: do_output
    integer, intent(in), optional :: N_opt ! length of time-period N, optional

! *** local variable ***
    integer :: N_nec                       ! length of time-period N, necessary


    if (present(N_opt)) then
       N_nec = N_opt
    else
       N_nec = 1
    endif

    if (day_in_month==num_day_in_month(fleapyear,month) .and. &
         timenew==assim_time                             .and. &
         (mod(month, N_nec)==0)) then
       do_output=.true.
    else
       do_output=.false.
    end if

  end subroutine monthly_event_assimstep

end module utils_pdaf
