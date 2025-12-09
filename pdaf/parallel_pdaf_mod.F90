!> Setup parallelisation
!!
!! This module provides variables for the MPI parallelization
!! to be shared between model-related routines. There are variables
!! that are used in the model even without PDAF, and additional variables
!! that are only used if data assimilaion with PDAF is performed.
!! The initialization of communicators for execution with PDAF is
!! performed in `init_parallel_pdaf`.
!!
!! __Revision history:__
!! * 2004-10 - Lars Nerger  - Initial code
!! * 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
!! * 2004-10 - Lars Nerger  - Revision for PDAF3
!! 
module mod_parallel_pdaf

  implicit none
  save 

  include 'mpif.h'
  
  ! Particular variables for use with AWI-CM
  logical :: pairs       !< Whether we use pairs of fesom.x and OpenIFS.x
  logical :: writepe     !< Whether the process does file writing

  ! Basic variables for model state integrations
  integer :: COMM_model         !< MPI communicator for model tasks
  integer :: mype_model         !< Rank in COMM_model
  integer :: npes_model         !< Size of COMM_model

  integer :: COMM_ensemble      !< Communicator for entire ensemble
  integer :: mype_ens           !< Rank in COMM_ensemble
  integer :: npes_ens           !< Size of COMM_ensemble

  integer :: mype_world         !< Number of PEs in MPI_COMM_WORLD
  integer :: npes_world         !< PE rank in MPI_COMM_WORLD

  integer :: mype_submodel      !< Number of PEs in model compartment task

  ! Additional variables for use with PDAF
  integer :: n_modeltasks = 1   !< Number of parallel model tasks

  integer :: COMM_filter        !< MPI communicator for filter PEs 
  integer :: mype_filter        !< Rank in COMM_filter
  integer :: npes_filter        !< Size of COMM_filter

  integer :: COMM_couple !< MPI communicator for coupling filter and model
  integer :: mype_couple        !< Rank in COMM_couple
  integer :: npes_couple        !< Size in COMM_couple

  logical :: modelpe            !< Whether we are on a PE in a COMM_model
  logical :: filterpe           !< Whether we are on a PE in a COMM_filter
  integer :: task_id            !< Index of my model task (1,...,n_modeltasks)
  character(len=10) :: task_str !< Task ID as string
  integer :: MPIerr             !< Error flag for MPI
  integer :: MPIstatus(MPI_STATUS_SIZE)       !< Status array for MPI
  integer, allocatable :: local_npes_model(:) !< # PEs per ensemble

!-------------------------------------------------------------------------------

contains
!-------------------------------------------------------------------------------
!> Initialize MPI communicators for PDAF
!!
!! Split the FESOM world MPI communicator into MODEL,
!! FILTER and COUPLE communicators. Provide the 
!! communicators to PDAF by calling PDAF3_set_parallel.
!! Set communicator MPI_COMM_FESOM as well as npes, mype
!! for FESOM.
!!
!! The routine supports AWI-CM coupled with OASIS.
!!
  subroutine init_parallel_pdaf(MPI_COMM_FESOM, mype_fesom, npes_fesom)

    use mpi                         ! MPI
    use timer, only: timeit
    use PDAF, &                     ! PDAF routines
         only: PDAF3_set_parallel
    use mod_assim_pdaf, &
         only: dim_ens
#if defined(__oasis)
    use mod_oasis_data, &
         only: comm_cplmod, cplmodid, strcplmodid ! For AWI-CM
#endif

    implicit none
  
! *** Arguments ***
    integer, intent(inout) :: MPI_COMM_FESOM
    integer, intent(inout) :: mype_fesom
    integer, intent(inout) :: npes_fesom

! *** Local variables ***
    integer :: i, j                     ! Counters
    integer :: mype_couple, npes_couple ! Rank and size in COMM_couple
    integer :: pe_index                 ! Index of PE
    integer :: my_color, color_couple   ! Variables for communicator-splitting 
    logical :: iniflag                  ! Flag whether MPI is initialized
    integer :: flag                     ! Status flag
    character(len=32) :: handle         ! handle for command line parser
    integer :: screen=0                 ! Whether information is shown on the screen
                                        ! 1 fro active display (set in namelist)

    integer provided_mpi_thread_support_level
    character(:), allocatable :: provided_mpi_thread_support_level_name

#if defined(__oasis)  
    integer :: ncpus_atm          ! number of processes for each atmosphere task
    integer :: ncpus_ocn          ! number of processes for each ocean task
#endif
  
    call timeit(7, 'ini')
    call timeit(1, 'new')
    call timeit(2, 'new')
  
    ! *** Initialize MPI if not yet initialized ***
    call MPI_Initialized(iniflag, MPIerr)
    if (.not.iniflag) then
       call MPI_Init(MPIerr)
    end if

    ! *** Initialize PE information on COMM_world ***
    call MPI_Comm_size(MPI_COMM_WORLD, npes_world, MPIerr)
    call MPI_Comm_rank(MPI_COMM_WORLD, mype_world, MPIerr)
  
#if defined(__oasis)
    namelist /pdaf_parallel/ dim_ens, ncpus_atm, ncpus_ocn, pairs, screen
    if (mype_world == 0) then
       write(*,*) 'init_parallel_pdaf: oasis defined (i.e. fesom coupling for Awi-Cm ON).'
    end if
#endif
#if !defined(__oasis)
    namelist /pdaf_parallel/ dim_ens, screen
    if (mype_world == 0) then
       write(*,*) 'init_parallel_pdaf: oasis NOT defined (i.e. fesom coupling for Awi-Cm OFF).'
    end if
#endif

    ! *** Get parallelization information for ensemble runs ***
    open (10, file='namelist.fesom.pdaf')
    read (10, NML=pdaf_parallel)
    close (10)
  
    n_modeltasks = dim_ens

#if defined(__oasis)
    ! *** Initialize communicators for ensemble evaluations ***
    if (mype_world == 0) &
         write (*, '(/1x, a)') 'FESOM-PDAF: Initialize communicators for assimilation with PDAF'

    if (mype_world == 0) then
       write (*,'(a,3x,a,i6)') 'FESOM-PDAF', 'model tasks: ', n_modeltasks
       write (*,'(a,3x,a,i6)') 'FESOM-PDAF', 'CPUs per ocean task:      ', ncpus_ocn
       write (*,'(a,3x,a,i6)') 'FESOM-PDAF', 'CPUs per atmosphere task: ', ncpus_atm
       write (*,'(a,3x,a,i6)') 'FESOM-PDAF', 'CPUs per runoff mappertask: ', 1 
    end if
#endif


    ! *** Check consistency of number of parallel ensemble tasks ***
    consist1: if (n_modeltasks > npes_world) then
       ! *** # parallel tasks is set larger than available PEs ***
       n_modeltasks = npes_world
       if (mype_world == 0) write (*, '(3x, a)') &
            '!!! Resetting number of parallel ensemble tasks to total number of PEs!'
    end if consist1
    if (dim_ens > 0) then
       ! Check consistency with ensemble size
       consist2: if (n_modeltasks > dim_ens) then
          ! # parallel ensemble tasks is set larger than ensemble size
          n_modeltasks = dim_ens
          if (mype_world == 0) write (*, '(5x, a)') &
               '!!! Resetting number of parallel ensemble tasks to number of ensemble states!'
       end if consist2
    end if

    ! ***              COMM_ENSEMBLE                ***
    ! *** Generate communicator for ensemble runs   ***
    ! *** only used to generate model communicators ***
    COMM_ensemble = MPI_COMM_WORLD

    npes_ens = npes_world
    mype_ens = mype_world


    ! *** Store # PEs per ensemble                 ***
    ! *** used for info on PE 0 and for generation ***
    ! *** of model communicators on other Pes      ***
    allocate(local_npes_model(n_modeltasks))
  
#if defined(__oasis)
    if (pairs) then
       local_npes_model = ncpus_ocn + 1 + ncpus_atm
    else
       local_npes_model = ncpus_ocn
    end if
    if (n_modeltasks*(ncpus_ocn+1+ncpus_atm) /= npes_world) then
       write (*,*) 'ERROR: Inconsistent setting of processes for ensemble integrations - check namelist.fesom.pdaf.'
       call abort_parallel()
    end if
#endif

#if !defined(__oasis)
    local_npes_model = floor(real(npes_world) / real(n_modeltasks))
    do i = 1, (npes_world - n_modeltasks * local_npes_model(1))
       local_npes_model(i) = local_npes_model(i) + 1
    end do
#endif

    ! ***              COMM_MODEL               ***
    ! *** Generate communicators for model runs ***
    ! *** (Split COMM_ENSEMBLE)                 ***
    pe_index = 0
    doens1: do i = 1, n_modeltasks
       do j = 1, local_npes_model(i)
          if (mype_ens == pe_index) then
             task_id = i
             exit doens1
          end if
          pe_index = pe_index + 1
       end do
    end do doens1
    call MPI_Comm_split(COMM_ensemble, task_id, mype_ens, &
         COMM_model, MPIerr)
  
    ! *** Re-initialize PE informations   ***
    ! *** according to model communicator ***
    call MPI_Comm_Size(COMM_model, npes_model, MPIerr)
    call MPI_Comm_Rank(COMM_model, mype_model, MPIerr)
    if (screen > 1) then
       write (*,'(a,i5,a,i5,a,i5,a,i5)') 'FESOM-PDAF: mype(w)= ', mype_world, &
            '; model task: ', task_id, '; mype(m)= ', mype_model, '; npes(m)= ', npes_model
    end if


    ! Init flag FILTERPE (all PEs of model task 1)
    if (task_id == 1) then
       filterpe = .true.
    else
       filterpe = .false.
    end if

    ! ***         COMM_FILTER                 ***
    ! *** Generate communicator for filter    ***
    ! *** For simplicity equal to COMM_couple ***
    if (filterpe) then
       my_color = task_id
    else
       my_color = MPI_UNDEFINED
    endif

    call MPI_Comm_split(MPI_COMM_WORLD, my_color, mype_world, &
         COMM_filter, MPIerr)

    ! *** Initialize PE informations         ***
    ! *** according to coupling communicator ***
    if (filterpe) then
       call MPI_Comm_Size(COMM_filter, npes_filter, MPIerr)
       call MPI_Comm_Rank(COMM_filter, mype_filter, MPIerr)
    endif


    ! ***              COMM_COUPLE                 ***
    ! *** Generate communicators for communication ***
    ! *** between model and filter PEs             ***
    ! *** (Split COMM_ENSEMBLE)                    ***

    color_couple = mype_model + 1

    call MPI_Comm_split(MPI_COMM_WORLD, color_couple, mype_world, &
         COMM_couple, MPIerr)

    ! *** Initialize PE informations         ***
    ! *** according to coupling communicator ***
    call MPI_Comm_Size(COMM_couple, npes_couple, MPIerr)
    call MPI_Comm_Rank(COMM_couple, mype_couple, MPIerr)

    if (screen > 0) then
       if (mype_world == 0) then
          write (*, '(/a, 18x,a)') 'FESOM-PDAF Pconf','PE configuration:'
          write (*, '(a, 2x, a6, a9, a10, a14, a13, /a, 2x, a5, a9, a7, a7, a7, a7, a7, /a, 2x, a)') &
               'FESOM-PDAF Pconf', 'world', 'filter', 'model', 'couple', 'filterPE', &
               'FESOM-PDAF Pconf', 'rank', 'rank', 'task', 'rank', 'task', 'rank', 'T/F', &
               'FESOM-PDAF Pconf', '----------------------------------------------------------'
       end if

       call MPI_Barrier(MPI_COMM_WORLD, MPIerr)
       if (task_id == 1) then
          write (*, '(a, 2x, i6, 3x, i4, 4x, i3, 4x, i3, 4x, i3, 4x, i3, 4x, l3)') &
               'FESOM-PDAF Pconf', mype_world, mype_filter, task_id, mype_model, color_couple, &
               mype_couple, filterpe
       endif
       if (task_id > 1) then
          write (*,'(a, 2x, i6, 11x, i3, 4x, i3, 4x, i3, 4x, i3, 4x, l3)') &
               'FESOM-PDAF Pconf', mype_world, task_id, mype_model, color_couple, mype_couple, filterpe
       end if

       call MPI_Barrier(MPI_COMM_WORLD, MPIerr)

       if (mype_world == 0) write (*, '(/a)') ''

    end if


! ***************************************************
! *** Provide parallelization information to PDAF ***
! ***************************************************

    call PDAF3_set_parallel(COMM_ensemble, COMM_model, COMM_filter, COMM_couple, &
         task_id, n_modeltasks, filterpe, flag)


! ******************************************************************************
! *** Initialize model equivalents to COMM_model, npes_model, and mype_model ***
! ******************************************************************************

#if defined(__oasis)
    comm_cplmod = comm_model
    cplmodid = task_id ! Store task ID 
    write (strcplmodid,'(a1,i5.5)') '_',cplmodid ! Store string of coupled task ID
#endif
  
    MPI_COMM_FESOM = comm_model
    npes_fesom = npes_model
    mype_fesom = mype_model

    call timeit(2,'old')
    call timeit(3,'new')

  end subroutine init_parallel_pdaf

!-------------------------------------------------------------------------------
!> Initialize the MPI execution environment.
!!
  subroutine init_parallel()

    implicit none

    integer :: i
  
    call MPI_INIT(i);
    call MPI_Comm_Size(MPI_COMM_WORLD,npes_world,i)
    call MPI_Comm_Rank(MPI_COMM_WORLD,mype_world,i)

    ! Initialize model communicator, its size and the process rank
    ! Here the same as for MPI_COMM_WORLD
    comm_model = MPI_COMM_WORLD
    npes_model = npes_world
    mype_model = mype_world
   
  end subroutine init_parallel

!-------------------------------------------------------------------------------
!> Finalize the MPI execution environment.
!!
  subroutine finalize_parallel()

    implicit none
    
    call  MPI_Barrier(MPI_COMM_WORLD,MPIerr)
    call  MPI_Finalize(MPIerr)

  end subroutine finalize_parallel

!-------------------------------------------------------------------------------
!> Terminate the MPI execution environment.
!!
  subroutine abort_parallel()

    call MPI_Abort(MPI_COMM_WORLD, 1, MPIerr)

  end subroutine abort_parallel

end module mod_parallel_pdaf
