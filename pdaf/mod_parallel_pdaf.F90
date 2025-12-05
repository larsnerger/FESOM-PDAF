! PDAF in AWI-CM2 / Fesom 2.0

MODULE mod_parallel_pdaf

! !DESCRIPTION:
! This modules provides variables for the MPI parallelization
! to be shared between model-related routines. The are variables
! that are used in the model, even without PDAF and additional
! variables that are only used, if data assimialtion with PDAF
! is performed.
! In addition methods to initialize and finalize MPI are provided.
! The initialization routine is only for the model itself, the 
! more complex initialization of communicators for xecution with
! PDAF is peformed in init\_parallel\_pdaf.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger  - Initial code
! 2019-11 - Longjiang Mu - Initial commit for AWI-CM3

  !USES:
  IMPLICIT NONE
  SAVE 

  INCLUDE 'mpif.h'

  !PUBLIC DATA MEMBERS:
  ! Basic variables for model state integrations
  integer :: COMM_model         ! MPI communicator for model tasks
  integer :: mype_model         ! Rank in COMM_model
  integer :: npes_model         ! Size of COMM_model

  integer :: COMM_ensemble      ! Communicator for entire ensemble
  integer :: mype_ens           ! Rank in COMM_ensemble
  integer :: npes_ens           ! Size of COMM_ensemble

  INTEGER :: mype_world         ! Number of PEs in MPI_COMM_WORLD
  INTEGER :: npes_world         ! PE rank in MPI_COMM_WORLD

  INTEGER :: mype_submodel      ! Number of PEs in model compartment task

  ! Additional variables for use with PDAF
  INTEGER :: n_modeltasks = 1   ! Number of parallel model tasks

  INTEGER :: COMM_filter        ! MPI communicator for filter PEs 
  integer :: mype_filter        ! Rank in COMM_filter
  integer :: npes_filter        ! Size of COMM_filter

  INTEGER :: COMM_couple ! MPI communicator for coupling filter and model
  integer :: mype_couple        ! Rank in COMM_couple
  integer :: npes_couple        ! Size in COMM_couple

  logical :: modelpe            ! Whether we are on a PE in a COMM_model
  logical :: filterpe           ! Whether we are on a PE in a COMM_filter
  integer :: task_id            ! Index of my model task (1,...,n_modeltasks)
  character(len=10) :: task_str ! Task ID as string
  INTEGER :: MPIerr             ! Error flag for MPI
  INTEGER :: MPIstatus(MPI_STATUS_SIZE)       ! Status array for MPI
  INTEGER, ALLOCATABLE :: local_npes_model(:) ! # PEs per ensemble
  
  ! Additional variables for use with AWI-CM
  LOGICAL :: pairs       ! Whether we use pairs of fesom.x and OpenIFS.x
  LOGICAL :: writepe     ! Whether the process does file writing

!-------------------------------------------------------------------------------

CONTAINS
!
! !INTERFACE:
  SUBROUTINE abort_parallel()

! !DESCRIPTION:
! Routine to abort MPI program

    IMPLICIT NONE
    
    CALL  MPI_Abort(MPI_COMM_WORLD, 1, MPIerr)

  END SUBROUTINE abort_parallel

END MODULE mod_parallel_pdaf
