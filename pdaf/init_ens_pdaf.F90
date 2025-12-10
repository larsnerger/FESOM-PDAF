!> Initialize ensemble
!!
!! __Revision history:__
!! * 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
!! * 2022-02 - Frauke B     - Adapted for FESOM2.1
!!
subroutine init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

  use mpi
  use PDAF, &
       only: PDAF_seik_omega, PDAF_sampleens
  use assim_pdaf_mod, &
       only: file_init, path_init, read_inistate, file_inistate, varscale, &
       this_is_pdaf_restart, start_from_ENS_spinup, &
       perturb_ssh, perturb_u, &
       perturb_v, perturb_temp, perturb_salt, &
       perturb_DIC, perturb_Alk, perturb_DIN, perturb_O2, &
       ens_p_init
  use fesom_pdaf, &
       only: mesh_fesom, nlmax, topography_p, edim_nod2d, mydim_nod2d, &
       daynew, timenew, eta_n, tr_arr, uv, wvel, a_ice, tiny
  use statevector_pdaf, &
       only: id, sfields
  use parallel_pdaf_mod, &
       only: mype_filter, COMM_filter, abort_parallel, mype_world, MPIerr

  implicit none

  include 'netcdf.inc'

! *** Arguments ***
  integer, intent(in) :: filtertype                !< Type of filter to initialize
  integer, intent(in) :: dim_p                     !< PE-local state dimension
  integer, intent(in) :: dim_ens                   !< Size of ensemble
  real, intent(inout) :: state_p(dim_p)            !< PE-local model state
  real, intent(inout) :: Uinv(dim_ens-1,dim_ens-1) !< Array not referenced for SEIK
  real, intent(out)   :: ens_p(dim_p, dim_ens)     !< PE-local state ensemble
  integer, intent(inout) :: flag                   !< PDAF status flag

! *** local variables ***
  integer :: i, j, member, s, row, col, k, n ! Counters
  integer :: rank                      ! Rank stored in init file
  integer :: dim_p_cov                 ! Local state dimension in covariance file
  integer,allocatable :: offsets_cov(:)! Local field offsets in covariance file
  integer :: fileid                    ! ID for NetCDF file
  integer :: id_dim                    ! ID for dimension
  integer :: id_state,id_svals,id_eof  ! IDs for fields
  integer :: startv(2),countv(2)       ! Vectors for reading fields
  integer :: stat(7)                   ! Status flag for NetCDF commands
  real :: fac                          ! Square-root of dim_ens or dim_ens-1
  real,allocatable :: eof_p(:,:)       ! Matrix of eigenvectors of covariance matrix
  real,allocatable :: svals(:)         ! Singular values
  real,allocatable :: omega(:,:)       ! Transformation matrix Omega
  character(len=5)   :: mype_string    ! String for process rank
  character(len=5)   :: col_string     ! String for ensemble member
  character(len=150) :: infile         ! File holding initial state estimate
  integer :: dim_p_read
  logical :: runningmean               ! True: Initialize state vector from
                                       ! nc-file running mean
  real, allocatable :: ens_p_per(:,:)  ! Ensemble of field perturbations
  integer :: o1,o2,f1,f2               ! Indeces
  integer, allocatable :: id_per_mod(:)! Field indices of to-be-perturbed fields in model state vector
  integer, allocatable :: id_per_cov(:)! Field indices of to-be-perturbed fields in covariance file
  integer, allocatable :: id_cov_mod(:)! Field indices of covariance fields in model state vector
  integer :: nfields_per               ! Number of to-be-perturbed fields
  integer :: nfields_cov               ! Number of fields in covariance file

  integer :: type_generate_ens=0       ! 1 to use PDAF_sample_ens, 0 old init using PDAF_seik_omega

  type field_ids_cov                   ! Field IDs for covariance matrix
     integer :: ssh
     integer :: u
     integer :: v
     integer :: temp
     integer :: salt
     integer :: DIC
     integer :: Alk
     integer :: DIN
     integer :: O2
  end type field_ids_cov
  
  type(field_ids_cov) :: id_cov        ! Type variable holding field IDs in covariance matrix
  
  logical, allocatable :: perturb_this(:)
  
  integer :: n_treshold_ssh_p,  &
             n_treshold_temp_p, &
             n_treshold_salt_p, &
             n_treshold_sic_p          ! Local counters for treshold-based corrections
             
  integer :: n_treshold_ssh_g,  &
             n_treshold_temp_g, &
             n_treshold_salt_g, &
             n_treshold_sic_g          ! Global counters for treshold-based corrections
             
  ! Debugging:
  logical            :: debugmode
  logical            :: write_debug
  integer            :: fileID_debug
  character(len=3)   :: day_string
  character(len=5)   :: tim_string
  
  ! Set debug output
  debugmode    = .false.
  write_debug = .false.
  if (debugmode .and. mype_world==0) write_debug = .true.
  
  if (write_debug) then
     ! print state vector
     write(day_string, '(i3.3)') daynew
     write(tim_string, '(i5.5)') int(timenew)
     fileID_debug=10
     open(unit=fileID_debug, file='init_ens_pdaf_'//day_string//'_'//tim_string//'.txt', status='unknown')
  endif


! **********************
! *** INITIALIZATION ***
! **********************

  if (mype_filter == 0)  write(*,'(a,2x,a)') 'FESOM-PDAF','*** INIT_ENS_PDAF ***'

  type_init: if (this_is_pdaf_restart) then

     if (mype_filter == 0)  write(*,'(a,2x,a)') 'FESOM-PDAF','This is a restart, skipping init_ens_pdaf'
     ens_p   = ens_p_init
     state_p = sum(ens_p, dim=2) / real(dim_ens)

  elseif (start_from_ENS_spinup) then type_init

     if (mype_filter == 0)  write(*,'(a,2x,a)') 'FESOM-PDAF','Starting from a perturbed ensemble, skipping init_ens_pdaf'
     ens_p   = ens_p_init
     state_p = sum(ens_p, dim=2) / real(dim_ens)

  else type_init

     ! ************************************************
     ! *** Generate ensemble from covariance matrix ***
     ! ************************************************

     ! Fields in covariance input file:
     nfields_cov = 9

     id_cov% ssh  = 1
     id_cov% u    = 2
     id_cov% v    = 3
     id_cov% temp = 4
     id_cov% salt = 5
     id_cov% DIC  = 6
     id_cov% Alk  = 7
     id_cov% DIN  = 8
     id_cov% O2   = 9

     ! Positions of covariance fields in the full model state vector:
     allocate(id_cov_mod(nfields_cov))
     id_cov_mod = (/ id%ssh,  &
          id%u,    &
          id%v,    &
          id%temp, &
          id%salt, &
          id%DIC,  &
          id%Alk,  &
          id%DIN,  &
          id%O2      /)

     ! Dimension of covariance matrix and EOF must be sum of field dimensions:
     dim_p_cov =  0
     do j=1, nfields_cov
        dim_p_cov = dim_p_cov + sfields(id_cov_mod(j))%dim
     end do
  
     ! Offsets in covariance matrix:
     allocate(offsets_cov(nfields_cov))
     offsets_cov(1) = 0
     do j=2, nfields_cov
        offsets_cov(j) = offsets_cov(j-1) + sfields(id_cov_mod(j-1))%dim
     enddo
  
     ! Positions of to-be-perturbed state fields in the full model state vector:
     ! note: you might take out some fields here, code will still work and not perturb these
     !       more simple: you might just set perturb_field=.false. in namelist.fesom.pdaf
  
     nfields_per = 9
     allocate(id_per_mod(nfields_per))
     id_per_mod = (/ id%ssh,  &
          id%u,    &
          id%v,    &
          id%temp, &
          id%salt, &
          id%DIC,  &
          id%Alk,  &
          id%DIN,  &
          id%O2      /)

     ! Positions of to-be-perturbed fields in covariance vector:
     allocate(id_per_cov(nfields_per))
     id_per_cov = (/ id_cov%ssh,  &
          id_cov%u,    &
          id_cov%v,    &
          id_cov%temp, &
          id_cov%salt, &
          id_cov%DIC,  &
          id_cov%Alk,  &
          id_cov%DIN,  &
          id_cov%O2      /)

     ! Not-to-be perturbed fields according to namelist settings:
     allocate(perturb_this(nfields_per))
     perturb_this(:) = .false.
  
     do j = 1,nfields_per 
        if ((id_per_mod(j)==id%ssh ) .and. perturb_ssh ) perturb_this(j) = .true.
        if ((id_per_mod(j)==id%u   ) .and. perturb_u   ) perturb_this(j) = .true.
        if ((id_per_mod(j)==id%v   ) .and. perturb_v   ) perturb_this(j) = .true.
        if ((id_per_mod(j)==id%temp) .and. perturb_temp) perturb_this(j) = .true.
        if ((id_per_mod(j)==id%salt) .and. perturb_salt) perturb_this(j) = .true.
        if ((id_per_mod(j)==id%DIC ) .and. perturb_DIC ) perturb_this(j) = .true.
        if ((id_per_mod(j)==id%Alk ) .and. perturb_Alk ) perturb_this(j) = .true.
        if ((id_per_mod(j)==id%DIN ) .and. perturb_DIN ) perturb_this(j) = .true.
        if ((id_per_mod(j)==id%O2  ) .and. perturb_O2  ) perturb_this(j) = .true.
     enddo
  
     ! Print out message
     if (mype_filter==0) then
        do j=1,nfields_cov
           write (*,'(1X,A40,1X,A5,1X,I7,1X,I7,L2)') 'PE0-dim and offset of covariance field ', &
                trim(sfields(id_cov_mod(j))%variable), sfields(id_cov_mod(j))%dim, &
                offsets_cov(j), perturb_this(j)
        enddo
        do j=1,nfields_per
           write (*,'(1X,A40,1X,A5,1X,I7,1X,I7,L2)') 'PE0-dim and offset of to-be-perturbed fields ', &
                trim(sfields(id_per_mod(j))%variable), &
                sfields(id_per_mod(j))%dim, sfields(id_per_mod(j))%off, perturb_this(j)
        enddo
     endif
  
     ! Print initialization message:
     write(mype_string,'(i4.4)') mype_filter
     mype0: if (mype_filter == 0) then
        write (*, '(/a, 8x,a)') 'FESOM-PDAF', 'Generate state ensemble from covariance matrix'
        write (*, '(a, 8x,a)') &
             'FESOM-PDAF', '--- use 2nd order exact sampling (SEIK type)'
        write (*, '(a, 8x,a,i5)') 'FESOM-PDAF', '--- number of EOFs:',dim_ens-1
     end if mype0

     ! Allocate memory for temporary fields
     allocate(eof_p(dim_p_cov, dim_ens-1))
     allocate(svals(dim_ens-1))
     allocate(omega(dim_ens, dim_ens-1))


! *************************************************
! *** Initialize covar matrix                   ***
! *************************************************
 
     infile=trim(path_init)//trim(file_init)//trim(mype_string)//'.nc'
  
     if (mype_filter == 0) write (*,'(a, 8x,a)') 'FESOM-PDAF', '--- Read initial state from file ', infile
     s = 1
     stat(s) = NF_OPEN(infile, NF_NOWRITE, fileid)
     do i = 1,  s
        if (stat(i) /= NF_NOERR) then
           write(*, *) 'NetCDF error in opening initialization file, no.', i
           stop
        end if
     end do

     ! Read size of state vector
     s = 1
     stat(s) = NF_INQ_DIMID(fileid, 'dim_state', id_dim)
     s = s + 1
     stat(s) = NF_INQ_DIMLEN(fileid, id_dim, dim_p_read)

     ! Read rank stored in file
     s = s + 1
     stat(s) = NF_INQ_DIMID(fileid, 'rank', id_dim)
     s = s + 1
     stat(s) = NF_INQ_DIMLEN(fileid, id_dim, rank)

     do i = 1,  s
        if (stat(i) /= NF_NOERR) &
             write(*, *) 'NetCDF error in reading dimensions from file, no.', i
     end do

     checkdim: if (dim_p_read == dim_p_cov .and. rank >= dim_ens-1) then

        if (mype_filter == 0) write (*,'(8x,a)') '--- Read covariance matrix'

        ! Inquire IDs for mean state, singular vectors and values
        s = 1
        stat(s) = NF_INQ_VARID(fileid, 'running_meanstate', id_state)
        s = s + 1
        stat(s) = NF_INQ_VARID(fileid, 'V', id_eof)
        s = s + 1
        stat(s) = NF_INQ_VARID(fileid, 'sigma', id_svals)

        ! EOF and singular values
        startv(2) = 1
        countv(2) = dim_ens-1
        startv(1) = 1
        countv(1) = dim_p_read
        s = s + 1
        stat(s) = NF_GET_VARA_DOUBLE(fileid, id_eof, startv, countv, eof_p)

        s = s + 1
        stat(s) = NF_GET_VARA_DOUBLE(fileid, id_svals, 1, dim_ens-1, svals)

        s = s + 1
        stat(s) = nf_close(fileid)

        do i = 1,  s
           if (stat(i) /= NF_NOERR) &
                write(*, *) 'NetCDF error in reading initialization file, no.', i
        end do

        if (mype_filter==0) then
           write(*,*) 'svals', svals
        end if


! *******************************
! *** Generate ensemble       ***
! *******************************

        ! *** Generate perturbations ***

        allocate(ens_p_per(dim_p_cov,dim_ens))
        ens_p_per   = 0.0

        if (dim_ens>1) then

           if (type_generate_ens==1) then

              ! Use PDAF routine to generate perturbations from covariance matrix

              call PDAF_SampleEns(dim_p_cov, dim_ens, eof_p, svals, state_p, ens_p_per, 1, flag)

           else

              ! Old initialization type explicitly using Omega

              if (mype_filter==0) then

                 write (*,'(a,8x,a)') 'FESOM-PDAF','--- generate state ensemble'

                 ! *** Generate uniform orthogonal matrix OMEGA ***
                 call PDAF_seik_omega(dim_ens-1, Omega, 1, 1)

                 ! ***      Generate ensemble of states         ***
                 ! *** x_i = x + sqrt(FAC) eofV (Omega C^(-1))t ***

                 ! A = Omega C^(-1)
                 do col = 1, dim_ens-1
                    do row = 1, dim_ens
                       Omega(row, col) = Omega(row,col) * svals(col)
                    end do
                 end do
              end if

              call MPI_Bcast(Omega, dim_ens*(dim_ens-1), MPI_DOUBLE_PRECISION, 0, &
                   COMM_filter, MPIerr)

              ! *** state_ens = state + sqrt(dim_ens-1) eofV A^T ***
     
              fac = varscale * sqrt(real(dim_ens-1))
        
              ! =========          =====             =====        =============
              ! ens_p_per := fac * eof_p * transpose(Omega) + 0 * ens_p_per (0)

              call DGEMM('n', 't', dim_p_cov, dim_ens, dim_ens-1, &
                   fac, eof_p, dim_p_cov, Omega, dim_ens, 0.0, ens_p_per, dim_p_cov) ! matrix operation
           end if
        else
           if (mype_filter==0) write (*,'(a,8x,a)') 'FESOM-PDAF','--- ensemble size 1: no state ensemble generated'
        end if ! if (dim_ens>1)


        ! *** Set ensemble mean state to state collected from filter PEs ***

        call collect_state_PDAF(dim_p, state_p)
        do col = 1,dim_ens
           ens_p(1:dim_p,col) = state_p(1:dim_p)
        end do

   
        ! *** Add perturbation to to-be-perturbed part of full state vector ***

        do j = 1,nfields_per   

           ! field offset in covariance matrix
           o1 = offsets_cov(id_per_cov(j)) + 1
           o2 = offsets_cov(id_per_cov(j)) + sfields(id_per_mod(j))%dim

           ! field offset in full state vector
           f1 = sfields(id_per_mod(j))%off + 1
           f2 = sfields(id_per_mod(j))%off + sfields(id_per_mod(j))%dim
      
           if (write_debug) then
              do n=f1,f2
                 write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6,G15.6,G15.6)') &
                      trim(sfields(id_per_mod(j))%variable), n, ens_p(n,:)
              enddo
           endif

           ! add perturbation, if specified so in namelist
           if (perturb_this(j)) then
              ! add perturbation
              ens_p(f1:f2,:) = ens_p(f1:f2,:) + ens_p_per(o1:o2,:)
           endif

           if (write_debug) then
              do n=f1,f2
                 write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6,G15.6,G15.6)') &
                      trim(sfields(id_per_mod(j))%variable), n, ens_p(n,:)
              enddo
           endif

           if (mype_filter==0) then
              write (*,'(1X,A15,1X,A5,1X,I7,1X,I7,1X,A15,1X,I7,1X,I7)') 'Index in covar: ', &
                   trim(sfields(id_per_mod(j))%variable), o1, o2, &
                   'Index in state: ', f1, f2
           endif

        enddo
   
        deallocate(ens_p_per)


        ! *** Enforce treshold values ***

        treshold: do col= 1,dim_ens
   
           ! consider model topography
           ens_p(:,col) = ens_p(:,col) * topography_p

           ! surface fields

           n_treshold_ssh_p = 0
           fields_surf: do i= 1,myDim_nod2D

              ! SSH: set to +/- 1.7m where larger than that
              ! (note: control simulation has minimum of -1.34 in Jan 2016; and -1.67 to 1.56 full-year max.)
              if ( ens_p(i+sfields(id% ssh)%off,col) < -1.7 ) n_treshold_ssh_p = n_treshold_ssh_p + 1 ! counter
              if ( ens_p(i+sfields(id% ssh)%off,col) > +1.7 ) n_treshold_ssh_p = n_treshold_ssh_p + 1 ! counter

              ens_p(i+sfields(id% ssh)%off,col)=min(max(ens_p(i+sfields(id% ssh)%off,col),-1.7),1.7)       ! apply correction
           end do fields_surf

           ! write out correction counters for SSH and sea ice:
           call MPI_Allreduce(n_treshold_ssh_p, n_treshold_ssh_g, 1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
           if (mype_filter == 0) write(*,*) 'FESOM-PDAF init_ens_pdaf - SSH  - Number of field corrections: ', &
                n_treshold_ssh_g, ' on member ', col
      
           ! 3D fields

           n_treshold_temp_p = 0
           n_treshold_salt_p = 0
           fields_3D: do i = 1, (nlmax) * myDim_nod2D 

              ! temp: set to -1.895 Celsius where smaller than that
              if ( ens_p(i+sfields(id% temp)%off,col) < -1.895D0) n_treshold_temp_p = n_treshold_temp_p + 1 ! counter
              ens_p(i+sfields(id% temp)%off,col)= max(ens_p(i+sfields(id% temp)%off,col),-1.895D0)          ! apply correction

              ! salt: set to null where negative
              if ( ens_p(i+sfields(id% salt)%off,col) < 0.0) n_treshold_salt_p = n_treshold_salt_p + 1      ! counter
              ens_p(i+sfields(id% salt)%off,col)= max(ens_p(i+sfields(id% salt)%off,col),0.)                ! apply correction

              ! BGC: set to "tiny" where negative
              ens_p(i+ sfields(id%O2 )%off,col) = max(tiny     ,ens_p(i+ sfields(id%O2 )%off,col))
              ens_p(i+ sfields(id%DIN)%off,col) = max(tiny*1e-3,ens_p(i+ sfields(id%DIN)%off,col))
              ens_p(i+ sfields(id%DIC)%off,col) = max(tiny*1e-3,ens_p(i+ sfields(id%DIC)%off,col))
              ens_p(i+ sfields(id%Alk)%off,col) = max(tiny*1e-3,ens_p(i+ sfields(id%Alk)%off,col))
           end do fields_3D
      
           ! write out correction counters for temperature and salinity:
           call MPI_Allreduce(n_treshold_temp_p, n_treshold_temp_g, 1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
           call MPI_Allreduce(n_treshold_salt_p, n_treshold_salt_g, 1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
           if (mype_filter == 0) then
              write(*,*) 'FESOM-PDAF init_ens_pdaf - Temp - Number of field corrections: ', n_treshold_temp_g, ' on member ', col
              write(*,*) 'FESOM-PDAF init_ens_pdaf - Salt - Number of field corrections: ', n_treshold_salt_g, ' on member ', col
           end if
        end do treshold

     else checkdim

        ! *** Rank stored in file is smaller than requested EOF rank ***
        write(*,*) 'FESOM-PDAF: ','ERROR: Rank stored in file is smaller than requested EOF rank ...'
        write(*,*) 'FESOM-PDAF: ','... or dim_p not equal--------> init_ens_pdaf'
        call abort_parallel()

     end if checkdim


! ****************
! *** clean up ***
! ****************

     deallocate(svals, eof_p, omega, id_per_cov, id_per_mod, id_cov_mod, offsets_cov)

  endif type_init

  if (write_debug) close(fileID_debug)

end subroutine init_ens_pdaf
  
