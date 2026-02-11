!> Module for timings
!!
!! This module provides methods to perform timings of 
!! parts of a program execution. It uses the intrinsic 
!! function SYSTEM_CLOCK.
!!
!! __Revision history:__
!! * 2000-11 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
module timer

  implicit none
  save
  
  public :: timeit, time_tot, time_temp

  private
  integer :: t_rate
  integer, allocatable :: t_start(:), t_end(:)
  real, allocatable    :: t_total(:), t_temp(:)

contains
!-------------------------------------------------------------------------------
!> Initialize Counters and time regions
!!
!! Subroutine to initialize counters and to perform timing of a region
!! specified by timerID.
!! Usage:\\
!!   CALL PDAF_timeit(N,'ini') - Allocates and initializes N counters\\
!!   CALL PDAF_timeit(M,'new') - Start timing region for counter M\\
!!   CALL PDAF_timeit(M,'old') - End timing region for counter M\\
!!   CALL PDAF_timeit(M,'fin') - Finalized and deallocates all counters\\
!!
  subroutine timeit(timerID, operation)

    implicit none

    ! *** Arguments ***
    integer, intent(in) :: timerID             ! ID of timer
    character(len=3), intent(in) :: operation  ! Requested operation 


    ! Initialize timers
    if (operation == 'ini') then
       if ( .not. (allocated(t_start))) then
          allocate(t_start(timerID), t_end(timerID))
          allocate(t_total(timerID), t_temp(timerID))
       end if
        
       t_total = 0.0
    end if
    
    ! Begin timing region
    if (operation == 'new') then
       call system_clock(t_start(timerID))
    end if

    ! End timing region
    if (operation == 'old') then
       call system_clock(t_end(timerID), t_rate)
       t_temp(timerID) = real(t_end(timerID) - t_start(timerID)) &
            / real(t_rate)
       t_total(timerID) = t_total(timerID) + real(t_end(timerID) - &
            t_start(timerID)) / real(t_rate)
    end if
    
    ! Finalize timers
    if (operation == 'fin') then
       deallocate(t_start, t_end)
       deallocate(t_total, t_temp)
    end if
    
  end subroutine timeit

!-------------------------------------------------------------------------------
!> Read out timers for last timing interval
!!
!! Read out the value of the timer in seconds for the last 
!! passage of the timing region defined by timerID.
!!
  real function time_temp(timerID)

    implicit none

! *** Arguments ***
    integer, intent(in) :: timerID             ! ID of timer

    time_temp = t_temp(timerID)

  end function time_temp

!-------------------------------------------------------------------------------
!> Read out total time of a timing region.
!!
!! Read out the accumulated value of the timer in seconds
!! for the timing region define by timerID.
!!
    real function time_tot(timerID)

    implicit none

! *** Arguments ***
    integer, intent(in) :: timerID             ! ID of timer

    time_tot = t_total(timerID)

  end function time_tot

end module timer
