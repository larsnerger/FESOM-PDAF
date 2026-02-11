SUBROUTINE U_init_obs(step, dim_obs_f, observation_f)


  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step        !< Current time step
  INTEGER, INTENT(in) :: dim_obs_f   !< Dimension of full observation vector
  REAL, INTENT(out)   :: observation_f(dim_obs_f) !< Full observation vector

! *** local variables ***
  INTEGER :: i                ! Loop counter
  INTEGER :: offset_obs_f     ! Count offset of an observation type in full obs. vector
  INTEGER :: idummy           ! Dummy to prevent compiler warning


! ******************************************
! *** Initialize full observation vector ***
! ******************************************

END SUBROUTINE u_init_obs
