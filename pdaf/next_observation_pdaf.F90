!> Determining the Next Analysis Step
!!
!! The subroutine is called before each forecast phase.
!! It has to initialize the number of time steps to 
!! compute until the next available observation (`nsteps`).
!! It also set the exit-flag, indicating whether the data
!! assimilation process is completed such that the
!! ensemble loop in the model routine can be exited.
!! 
!! The routine is called by all processes.
!! 
!! !REVISION HISTORY:
!! 2004-10 - Lars Nerger - Initial code
!! 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
!! 2025-12 - Lars Nerger - update for PDAF3
!!
SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)

  USE parallel_pdaf_mod, &
       ONLY: mype_model, task_id
  USE assim_pdaf_mod, &
       ONLY: delt_obs_ocn, step_null, assim_time, steps_first_fcst
  USE recom_config, &
       ONLY: secondsperday

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: stepnow  !< Number of the current time step
  INTEGER, INTENT(out) :: nsteps   !< Number of time steps until next obs
  INTEGER, INTENT(out) :: doexit   !< Whether to exit forecasting (1 for exit)
  REAL, INTENT(out)    :: time     !< Current model (physical) time


! *************************************************************
! *** Determine number of time steps until next observation ***
! *************************************************************

  IF (stepnow==step_null) THEN
     ! at start, one assimilation step right away
     nsteps = steps_first_fcst
     assim_time = INT( REAL(nsteps) / REAL(delt_obs_ocn) * REAL(secondsperday))
  ELSE
      ! daily assimilation steps during model time loop
      nsteps=delt_obs_ocn
  ENDIF
  
  IF (mype_model==0 .AND. task_id==1) THEN
     WRITE (*,'(a,i8,a)') 'FESOM-PDAF: Next observation after ', nsteps ,' time steps'
  ENDIF


! *********************************
! *** Set current physical time ***
! *********************************

  time = 1.0

! *********************
! *** Set exit flag ***
! *********************

  doexit = 0

END SUBROUTINE next_observation_pdaf
