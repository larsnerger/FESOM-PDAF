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

!-------------------------------------------------------------------------------
!> Call-back routine for init_obsvar
!!
!! This routine calls the routine PDAFomi_init_obsvar_f
!! for each observation type
!!
SUBROUTINE u_init_obsvar(step, dim_obs_p, obs_p, meanvar)

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step          !< Current time step
  INTEGER, INTENT(in) :: dim_obs_p     !< PE-local dimension of observation vector
  REAL, INTENT(in) :: obs_p(dim_obs_p) !< PE-local observation vector
  REAL, INTENT(out)   :: meanvar       !< Mean observation error variance


END SUBROUTINE U_init_obsvar

!-------------------------------------------------------------------------------
!> Call-back routine for init_obsvar_l
!!
!! This routine calls the routine PDAFomi_init_obsvar_l
!! for each observation type
!!
SUBROUTINE init_obsvar_l(domain_p, step, dim_obs_l, obs_l, meanvar_l)

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p      !< Index of current local analysis domain
  INTEGER, INTENT(in) :: step          !< Current time step
  INTEGER, INTENT(in) :: dim_obs_l     !< Local dimension of observation vector
  REAL, INTENT(in) :: obs_l(dim_obs_l) !< Local observation vector
  REAL, INTENT(out)   :: meanvar_l     !< Mean local observation error variance

END SUBROUTINE init_obsvar_l

!-------------------------------------------------------------------------------
!> Call-back routine for obs_op
!!
!! This routine calls the observation-specific
!! routines obs_op_TYPE.
!!
subroutine u_obs_op(step, dim_p, dim_obs, state_p, ostate)

  implicit none

! *** Arguments ***
  integer, intent(in) :: step                 !< Current time step
  integer, intent(in) :: dim_p                !< PE-local state dimension
  integer, intent(in) :: dim_obs              !< Dimension of full observed state
  real, intent(in)    :: state_p(dim_p)       !< PE-local model state
  real, intent(inout) :: ostate(dim_obs)      !< PE-local full observed state

end subroutine u_obs_op

SUBROUTINE u_init_dim_obs(step, dim_obs)
  IMPLICIT NONE
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(out) :: dim_obs  !< Dimension of full observation vector
END SUBROUTINE u_init_dim_obs

SUBROUTINE u_init_obsvar_l(domain_p, step, dim_obs_l, obs_l, meanvar_l)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: domain_p      !< Index of current local analysis domain
  INTEGER, INTENT(in) :: step          !< Current time step
  INTEGER, INTENT(in) :: dim_obs_l     !< Local dimension of observation vector
  REAL, INTENT(in) :: obs_l(dim_obs_l) !< Local observation vector
  REAL, INTENT(out)   :: meanvar_l     !< Mean local observation error variance
END SUBROUTINE u_init_obsvar_l
