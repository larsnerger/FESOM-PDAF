!>  Initialize communicators for PDAF
!!
!! Parallelization routine for a model with attached PDAF. The subroutine is
!! called in the main program subsequently to the initialization of MPI. It
!! initializes MPI communicators for the model tasks, filter task and the
!! coupling between model and filter tasks. In addition some other variables 
!! for the parallelization are initialized.
!! The communicators and variables are handed over to PDAF in the call to 
!! PDAF_set_parallel toward the end of this routine.
!!
!! 3 Communicators are generated:
!! * _COMM_filter_: Communicator in which the filter itself operates
!! * _COMM_model_: Communicators for parallel model forecasts
!! * _COMM_couple_: Communicator for coupling between models and filter
!!
!! In addition there is the main communicator
!! * _COMM_ensemble_: The main communicator in which PDAF operates
!! COMM_ensemble is set to the communicator in which all model integration 
!! are computed. Typically, this is MPI_COMM_WORLD, but it can be defined
!! differently if the model only operators on a subset to MPI_COMM_WORLD.
!! This happens, e.g. if some processes are separated to operate an
!! I/O server or a model coupler for coupled model systems.
!!
!! Other variables that have to be initialized are:
!! * _filterpe_ - Logical: Does the PE execute the filter?
!! * _task_id_ - Integer: Index identifying the model task
!! * _my_ensemble_ - Integer: The index of the PE's model task
!! * _local_npes_model_ - Integer array holding numbers of PEs per model task
!!
!! For COMM_filter and COMM_model also the size of the communicators
!! (npes_filter and npes_model) and the rank of each process  (mype_filter,
!! mype_model) are initialized. These variables can be used in the model part 
!! of the program, but are not handed over to PDAF.
!!
!! __Revision history:__
!! * 2004-11 - Lars Nerger - Initial code
!! * 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
!! * 2025-12 - Lars Nerger - Revision for PDAF3
!!
subroutine init_parallel_pdaf(dim_ens_in, screen2)

  use mpi                         ! MPI
  use timer, only: timeit
  use PDAF, &                     ! PDAF routines
       only: PDAF3_set_parallel
  use mod_parallel_pdaf, &        ! PDAF parallelization variables
       only: mype_world, npes_world, mype_model, npes_model, &
       COMM_model, mype_filter, npes_filter, COMM_filter, filterpe, &
       n_modeltasks, local_npes_model, task_id, COMM_couple, MPIerr, &
       COMM_ensemble, mype_ens, npes_ens, pairs, abort_parallel
  use mod_assim_pdaf, &
       only: dim_ens
#if defined(__oasis)
  use mod_oasis_data, &
       only: comm_cplmod, cplmodid, strcplmodid ! For AWI-CM
#endif
  use g_PARSUP, &
       only: MPI_COMM_FESOM, npes, mype

  implicit none
  
! *** Arguments ***
  integer, intent(inout) :: dim_ens_in   !< Ensemble size
  !< Often the routine is called with dim_ens=0, because the real ensemble size
  !< is initialized later in the program. For dim_ens=0 no consistency check
  !< for the ensemble size with the number of model tasks is performed.
  integer, intent(in)    :: screen2       !< Whether screen information is shown (replced by namelist)

! *** Local variables ***
  integer :: i, j                     ! Counters
  integer :: mype_couple, npes_couple ! Rank and size in COMM_couple
  integer :: pe_index                 ! Index of PE
  integer :: my_color, color_couple   ! Variables for communicator-splitting 
  logical :: iniflag                  ! Flag whether MPI is initialized
  integer :: flag                     ! Status flag
  character(len=32) :: handle         ! handle for command line parser
  integer :: screen=0                 ! Whether information is shown o the screen
  
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
  namelist /pdaf_parallel/ dim_ens
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
  npes = npes_model
  mype = mype_model

  call timeit(2,'old')
  call timeit(3,'new')

end subroutine init_parallel_pdaf
