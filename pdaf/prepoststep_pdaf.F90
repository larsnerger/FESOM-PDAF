!>  Used-defined Pre/Poststep routine for PDAF
!!
!! The routine is called for global filters (e.g. ESTKF)
!! before the analysis and after the ensemble transformation.
!! For local filters (e.g. LESTKF) the routine is called
!! before and after the loop over all local analysis
!! domains.
!!
!! The routine provides full access to the state 
!! estimate and the state ensemble to the user.
!! Thus, user-controlled pre- and poststep 
!! operations can be performed here. For example 
!! the forecast and the analysis states and ensemble
!! covariance matrix can be analyzed, e.g. by 
!! computing the estimated variances.
!!
!! If a user considers to perform adjustments to the 
!! estimates (e.g. for balances), this routine is 
!! the right place for it.
!!
!! __Revision history:__
!! * 2010-07 - Lars Nerger  - Initial code
!! * 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
!! * 2022    - Frauke       - Adapted for FESOM2.0
!! * 2025-12 - Lars Nerger  - revision for PDAF3
!!
subroutine prepoststep_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

  use mpi
  use PDAF, &
       only: PDAF_reset_forget, PDAF_diag_stddev, PDAF_diag_variance
  use parallel_pdaf_mod, &
       only: mype_filter, npes_filter, COMM_filter, mype_world, &
       writepe, MPIerr
  use assim_pdaf_mod, &
       only: step_null, filtertype, dim_lag, loctype, &
       forget, resetforget, DAoutput_path, proffiles_o, &
       this_is_pdaf_restart, days_since_DAstart, delt_obs_ocn, &
       depth_excl, depth_excl_no
  use cfluxes_diags_pdaf, &
       only: factor_mass, factor_conc, cffields, cffieldsasml, &
       id_s_asml_alk, id_s_asml_dic, id_s_asml_deadmatter, id_s_asml_livingmatter, &
       cfluxes_diags_output_tmean_asml
  use fesom_pdaf, &
       only: mesh_fesom, nlmax, MPI_COMM_FESOM, &
       area_surf_glob, inv_area_surf_glob, &
       volo_full_glob, inv_volo_full_glob, &
       cellvol, mydim_nod2d, &
       myList_edge2D, myDim_edge2D, myList_nod2D, &
       gather_nod, hnode_new, SecondsPerDay, & 
       yearnew, num_day_in_month, fleapyear, month, &      ! clock
       day_in_month, timenew                               ! clock
  use statevector_pdaf, &
       only: id, nfields, sfields
  use utils_pdaf, &
       only: monthly_event_assimstep
  use means_pdaf, &
       only: update_means_pdaf, reset_means_pdaf, &
       monthly_state_f, monthly_state_a, monthly_state_m, &
       monthly_stddev_f, monthly_stddev_a, monthly_stddev_m
  use adaptive_lradius_pdaf, &
       only: adaptive_lradius_stats_pdaf
  use mod_atmos_ens_stochasticity, &
      only: stable_rmse
  use mod_nc_out_routines, &
       only: netCDF_out
  use output_config_pdaf, &
       only: w_dayensm, w_daymemb, w_monensm, w_monmemb, w_mm, w_sm, &
       mm, aa, ff, ii, sa, sf, si, sm, oo, dd, ee
  ! mean state forecast for observation exclusion criteria
  use obs_TSprof_EN4_pdafomi, &
       only: assim_o_en4_t, assim_o_en4_s, prof_exclude_diff, mean_temp_p
  use obs_sst_pdafomi, &
       only: assim_o_sst, sst_exclude_ice, sst_exclude_diff, &
             mean_ice_p, mean_sst_p
  use obs_sss_smos_pdafomi, &
        only: assim_o_sss, sss_exclude_ice, sss_exclude_diff, &
              mean_sss_p
  use obs_sss_cci_pdafomi, &
        only: assim_o_sss_cci, sss_cci_exclude_ice, sss_cci_exclude_diff, &
              mean_sss_cci_p
  use obs_chl_cci_pdafomi, &
        only: assim_o_chl_cci, chl_cci_exclude_ice, chl_cci_exclude_diff, &
              mean_chl_cci_p
  use obs_o2_merged_pdafomi, &
        only: assim_o_o2_merged, o2_merged_excl_absolute, o2_merged_excl_relative, &
              mean_O2_p
  use obs_n_merged_pdafomi, &
        only: assim_o_n_merged, n_merged_excl_absolute, n_merged_excl_relative, &
              mean_n_p
  use corrections_pdaf, &
       only: correct_state, store_stddev
              
  use netcdf

  implicit none
  save

! *** Arguments ***
  integer, intent(in) :: step           !< Current time step, starting from 0 at beginning of year
                                        !< (When the routine is called before
                                        !< the analysis, -step is provided.)
  integer, intent(in) :: dim_p          !< PE-local state dimension
  integer, intent(in) :: dim_ens        !< Size of state ensemble
  integer, intent(in) :: dim_ens_p      !< PE-local size of ensemble
  integer, intent(in) :: dim_obs_p      !< PE-local dimension of observation vector
  real, intent(inout) :: state_p(dim_p) !< PE-local forecast/analysis state
  !< The array 'state_p' is not generally not initialized in the case of SEIK.
  !< It can be used freely here.
  real, intent(inout) :: Uinv(dim_ens-1, dim_ens-1) !< Inverse of matrix U
  real, intent(inout) :: ens_p(dim_p, dim_ens)      !< PE-local state ensemble
  integer, intent(in) :: flag           !< PDAF status flag

! *** Local variables ***
  integer :: f,i,k,n,nz,s,member        ! Counters
  real :: invdim_ens                    ! Inverse ensemble size

  real :: diffm                         ! temporary array
  character(len=1) :: typestr           ! Character indicating call type (intial, forecast, analysis)
  character(len=26) :: timestr          ! String with step value and date
  logical :: write_now

  character(len=3) :: forana           ! String indicating forecast or analysis
  real :: stddev_g
  real, allocatable :: ens_stddev(:) ! estimated RMS errors
  integer :: istart, iend             ! stard and end index of a field in state vector
  integer :: pdaf_status              ! status flag
  
  real, allocatable :: stdev_p(:)  ! ensemble standard deviation at grid proints
  real :: stdevglob_temp           ! global full ocean average of ensemble standard deviation for temperature field
  ! regional average of local ensemble standard deviation
  real :: stdev_surf_p(nfields,nlmax), stdev_surf_g(nfields,nlmax)  ! on surfaces
  real :: stdev_volo_p(nfields),       stdev_volo_g(nfields)        ! on full ocean volume
  
  integer, parameter :: int0 = 0

  ! variables for debugging:
  logical :: debug
  logical :: write_debug


! **********************
! *** INITIALIZATION ***
! **********************

  if (mype_filter==0) then

     ! Create step-date string
     write (timestr,'(i7,2x,i4,a1,i2.2,a1,i2.2,2x,i2.2,a1,i2.2)') &
          step, yearnew,'-',month,'-',day_in_month,floor(timenew/3600.0),':',int(mod(timenew,3600.0)/60.0)

     ! Write info
     if ((step-step_null)==0) then
        if (.not.(this_is_pdaf_restart)) then
           write (*,'(a, 8x,a,1x,a26)') 'FESOM-PDAF', 'PrePost: Analyze initial state ensemble at step', timestr
           write (typestr,'(a1)') 'i'
        else
           write (*,'(a, 8x,a,1x,a26)') 'FESOM-PDAF', 'PrePost: This is a PDAF restart. No initial fields at step', timestr
           write (typestr,'(a1)') 'i'
        end if
        forana = 'ini'
     else if ((step-step_null)>0) then
        write (*,'(a, 8x,a,1x,a26)') 'FESOM-PDAF', 'PrePost: Analyze assimilated state ensemble at step', timestr
        write (typestr,'(a1)') 'a'
        forana = 'ana'
     else if ((step-step_null)<0) then
        write (*,'(a, 8x,a,1x,a26)') 'FESOM-PDAF', 'PrePost: Analyze forecast state ensemble at step', timestr
        write (typestr,'(a1)') 'f'
        forana = 'for'
     end if
  end if

  ! initialize numbers
  invdim_ens = 1.0 / real(dim_ens)


! ****************************
! *** Perform pre/poststep ***
! ****************************

  ! monthly event
  call monthly_event_assimstep(write_now)


! ***************************************************
! *** Corrections to state vector fields          ***
! ***************************************************

  call correct_state(step, dim_ens, ens_p)


! ************************************************************
! *** Compute ensemble mean and standard deviation         ***
! *** (=RMS errors according to sampled covar matrix)      ***
! ************************************************************

  if (mype_filter==0) write (*,'(a, 8x,a)') 'FESOM-PDAF', '--- compute ensemble mean and standard deviations'

  ! Allocate fields
  allocate(ens_stddev(nfields))

  ! Compute ensemble deviation and mean separately
  ! for each field in the state vector
  do i = 1, nfields
     ! Start and end index
     istart = 1 + sfields(i)%off
     iend = sfields(i)%dim + sfields(i)%off

     call PDAF_diag_stddev(sfields(i)%dim, dim_ens, &
          state_p(istart:iend), ens_p(istart:iend,:), &
          ens_stddev(i), 1, COMM_filter, pdaf_status)
  end do

  ! Output ensemble standard deviations
  if (mype_world == 0) then
     write (*, '(a,6x,a)') 'FESOM-PDAF', 'Ensemble standard deviation (estimated RMS error)'
     do i = 1, nfields
        write (*,'(a,4x,a10,4x,a10,2x,es12.4)') &
             'FEOSM-PDAF', 'STDDEV-'//forana, trim(sfields(i)%variable), ens_stddev(i)
     end do
  end if
  deallocate(ens_stddev)


! *********************************************************************
! *** Store ensemble mean values for observation exclusion criteria ***
! *********************************************************************
! save values at forecast phase

  store_fcst: if ((step-step_null)<0) then

     ! -- Sea-ice concentration --
     ! save mean_ice_p
     if (    (sst_exclude_ice     .and. assim_o_sst     )&
          .or. (sss_exclude_ice     .and. assim_o_sss     )&
          .or. (sss_cci_exclude_ice .and. assim_o_sss_cci )&
          .or. (chl_cci_exclude_ice .and. assim_o_chl_cci )) then 
        if (mype_filter==0) write (*,'(a, 8x,a)') &
             'FESOM-PDAF', '--- save ensemble mean forecast SEA-ICE for observation exclusion'
        if (allocated(mean_ice_p)) deallocate(mean_ice_p)
        allocate (mean_ice_p(sfields(id%a_ice)%dim))
        mean_ice_p = state_p(sfields(id%a_ice)%off + 1 : &
             sfields(id%a_ice)%off + sfields(id% a_ice)%dim)
     end if

     ! -- SST --
     ! save mean_sst_p
     if ((sst_exclude_diff > 0.0) .and. assim_o_sst) then
        if (mype_filter==0) write (*,'(a, 8x,a)') &
             'FESOM-PDAF', '--- save ensemble mean forecast SST for observation exclusion'
        if (allocated(mean_sst_p)) deallocate(mean_sst_p)
        allocate (mean_sst_p(myDim_nod2D))
        do i = 1, myDim_nod2D
           mean_sst_p(i) = state_p(sfields(id% temp)%off + (i-1) * (nlmax) + 1)
        end do
     end if

     ! -- SSS (CASE SMOS) --
     ! save mean_sss_p
     if ((sss_exclude_diff > 0.0) .and. assim_o_sss) then
        if (mype_filter==0) write (*,'(a, 8x,a)') &
             'FESOM-PDAF', '--- save ensemble mean forecast SSS for observation exclusion'
        if (allocated(mean_sss_p)) deallocate(mean_sss_p)
        allocate (mean_sss_p(myDim_nod2D))
        do i = 1, myDim_nod2D
           mean_sss_p(i) = state_p(sfields(id% salt)%off + (i-1) * (nlmax) + 1)
        end do
     end if

     ! -- SSS (CASE CCI) --
     ! save mean_sss_cci_p
     if ((sss_cci_exclude_diff > 0.0) .and. assim_o_sss_cci) then
        if (mype_filter==0) write (*,'(a, 8x,a)') &
             'FESOM-PDAF', '--- save ensemble mean forecast SSS for observation exclusion'
        if (allocated(mean_sss_cci_p)) deallocate(mean_sss_cci_p)
        allocate (mean_sss_cci_p(myDim_nod2D))
        do i = 1, myDim_nod2D
           mean_sss_cci_p(i) = state_p(sfields(id% salt)%off + (i-1) * (nlmax) + 1)
        end do
     end if
     
     ! -- Chlorophyll --
     ! save mean_chl_cci_p
     if ((chl_cci_exclude_diff > 0.0) .and. assim_o_chl_cci) then
        if (mype_filter==0) write (*,'(a, 8x,a)') &
             'FESOM-PDAF', '--- save ensemble mean forecast CHL for observation exclusion'
        if (allocated(mean_chl_cci_p)) deallocate(mean_chl_cci_p)
        allocate (mean_chl_cci_p(myDim_nod2D))
        do i = 1, myDim_nod2D
           mean_chl_cci_p(i) = state_p(sfields(id% PhyChl)%off + (i-1) * (nlmax) + 1) &
                + state_p(sfields(id% DiaChl)%off + (i-1) * (nlmax) + 1)
        end do
     end if
     
     ! -- 3D temperature field --
     ! save mean_temp_p
     if ((assim_o_en4_t .or. assim_o_en4_s) .and. (prof_exclude_diff > 0.0)) then
        if (mype_filter==0) write (*,'(a, 8x,a)') &
             'FESOM-PDAF', '--- save ensemble mean temperature (3D) for observation exclusion'
        ! Store mean temperature for profile assimilation
        if (allocated(mean_temp_p)) deallocate(mean_temp_p)
        allocate (mean_temp_p(sfields(id%temp)%dim))
        mean_temp_p = state_p(sfields(id%temp)%off+1 : sfields(id%temp)%off+sfields(id%temp)%dim)
     end if
     
     ! -- oxygen --
     ! save mean_o2_p
     if (((o2_merged_excl_absolute > 0.0) .or. (o2_merged_excl_relative > 0.0)) &
          .and. assim_o_o2_merged) then
        if (mype_filter==0) write (*,'(a, 8x,a)') &
             'FESOM-PDAF', '--- save ensemble mean forecast oxygen for observation exclusion'
        if (allocated(mean_O2_p)) deallocate(mean_O2_p)
        allocate (mean_O2_p(sfields(id%O2)%dim))
        mean_O2_p = state_p(sfields(id%O2)%off+1 : sfields(id%O2)%off+sfields(id%O2)%dim)
     end if
     
     ! -- nitrate --
     ! save mean_n_p
     if (((n_merged_excl_absolute > 0.0) .or. (n_merged_excl_relative > 0.0)) &
          .and. assim_o_n_merged) then
        if (mype_filter==0) write (*,'(a, 8x,a)') &
             'FESOM-PDAF', '--- save ensemble mean forecast DIN for observation exclusion'
        if (allocated(mean_n_p)) deallocate(mean_n_p)
        allocate (mean_n_p(sfields(id%DIN)%dim))
        mean_n_p = state_p(sfields(id%DIN)%off+1 : sfields(id%DIN)%off+sfields(id%DIN)%dim)
     end if

  end if store_fcst


! ************************************************
! *** Carbon sources minus sinks diagnostics   ***
! ************************************************
  
  ! factor to convert concentration to mass
  factor_mass = mesh_fesom%areasvol(:nlmax,:myDim_nod2D) * hnode_new(:nlmax,:myDim_nod2D) / SecondsPerDay
  factor_conc = 1.0 / SecondsPerDay

  if ((step-step_null)<0) then
  ! forecast phase
  ! get fmass and fconc before analysis step
  
     if (mype_filter == 0) &
          write(*, *) 'FESOM-PDAF', '--- compute carbon diagnostics at forecast'

     do i = 1, myDim_nod2D
        do k = 1, nlmax
           s = (i-1) * (nlmax) + k ! index in state vector
           ! DIC
           cffields(id_s_asml_dic)%fconc(k, i) = state_p(s + sfields(id%DIC)%off)
           ! Alk
           cffields(id_s_asml_alk)%fconc(k, i) = state_p(s + sfields(id%Alk)%off)
           ! Living carbon biomass
           cffields(id_s_asml_livingmatter)%fconc(k, i) = &
                (state_p(s + sfields(id%PhyC)%off) &
                + state_p(s + sfields(id%DiaC)%off) &
                + state_p(s + sfields(id%Zo1C)%off) &
                + state_p(s + sfields(id%Zo2C)%off) &
                + state_p(s + sfields(id%PhyCalc)%off))
           ! Dead organic carbon
           cffields(id_s_asml_deadmatter)%fconc(k, i) = &
                (state_p(s + sfields(id%DOC)%off)     &
                + state_p(s + sfields(id%DetC)%off)    &
                + state_p(s + sfields(id%DetCalc)%off) &
                + state_p(s + sfields(id%Det2C)%off)   &
                + state_p(s + sfields(id%Det2Calc)%off))
        enddo ! k=1,nlmax
     enddo ! i=1,my_Dim_nod2D
    
     ! convert concentration to mass
     do s=1, size(cffieldsasml)
        i = cffieldsasml(s)
        cffields(i)%fmass = cffields(i)%fconc * factor_mass
        cffields(i)%fconc = cffields(i)%fconc * factor_conc
     enddo

  endif ! (forecast phase)
  
  if ((step-step_null)>0) then
     ! analysis phase
     ! get amass and aconc after analysis step
  
     if (mype_filter == 0) &
          write(*, *) 'FESOM-PDAF', '--- compute carbon diagnostics at analysis'
  
     do i = 1, myDim_nod2D
        do k = 1, nlmax
           s = (i-1) * (nlmax) + k ! index in state vector
           ! DIC
           cffields(id_s_asml_dic)%aconc (k, i) = state_p(s + sfields(id%DIC)%off)
           ! Alk
           cffields(id_s_asml_alk)%aconc (k, i) = state_p(s + sfields(id%Alk)%off)
           ! Living carbon biomass
           cffields(id_s_asml_livingmatter)%aconc (k, i) = &
                (state_p(s + sfields(id%PhyC)%off) &
                + state_p(s + sfields(id%DiaC)%off) &
                + state_p(s + sfields(id%Zo1C)%off) &
                + state_p(s + sfields(id%Zo2C)%off) &
                + state_p(s + sfields(id%PhyCalc)%off))
           ! Dead organic carbon
           cffields(id_s_asml_deadmatter)%aconc (k, i) = &
                (state_p(s + sfields(id%DOC)%off)     &
                + state_p(s + sfields(id%DetC)%off)    &
                + state_p(s + sfields(id%DetCalc)%off) &
                + state_p(s + sfields(id%Det2C)%off)   &
                + state_p(s + sfields(id%Det2Calc)%off))
        enddo ! k=1,nlmax
     enddo ! i=1,my_Dim_nod2D

     ! convert concentration to mass
     do s=1, size(cffieldsasml)
        i = cffieldsasml(s)
        cffields(i)%amass = cffields(i)%aconc * factor_mass
        cffields(i)%aconc = cffields(i)%aconc * factor_conc
     enddo

     ! save the difference of forecast and analysis
     call cfluxes_diags_output_tmean_asml()
  
  endif ! (analysis phase)
  



! *****************************************************************
! *** Compute ensemble spread (STD) for different fields        ***
! *****************************************************************
  
  ! Set debug output
  debug = .false.
  write_debug = .false.
  if (debug .and. mype_world==0) write_debug = .true.

  
  ! *** Compute ensemble standard deviation field at grid points ***

  if (mype_filter == 0) &
      write(*, *) 'FESOM-PDAF', '--- compute ensemble standard deviation'

  allocate(stdev_p(dim_p))

  ! Note: This stddev is different from Frauke's code: we normalize with
  ! dim_ens-1 while Frauke uses dim_ens
  CALL PDAF_diag_variance(dim_p, dim_ens, state_p, ens_p, stdev_p, &
     stddev_g, 0, 0, COMM_filter, pdaf_status)
  stdev_p = sqrt(stdev_p)

  ! Store standard deviation of SSH used for corrections after analysis step
  CALL store_stddev(step, stdev_p)

  
  ! -----------------------------------------------------------------------------------------------------
  ! --- stdev_surf_g (nfields)    | layerwise surface mean of grid-point ensemble STD for each field area-weighted
  ! -----------------------------------------------------------------------------------------------------
  ! Compute pe-local surface mean of ensemble STD for each field
  if (mype_filter == 0) &
      write(*, *) 'FESOM-PDAF', '--- compute ensemble standard deviation surface mean'
  stdev_surf_p = 0.0
  stdev_surf_g = 0.0
  
  do f=1,nfields
     do n=1,myDim_nod2D
        if (sfields(f)%ndims == 2) then
           ! 3D fields
           do nz=1,nlmax
              stdev_surf_p(f,nz) =  stdev_surf_p(f,nz) &
                                   +  mesh_fesom%areasvol(nz,n) * stdev_p( sfields(f)%off + (n-1)*(nlmax) + nz )
           enddo ! nz,nlmax
        else
           ! surface fields
           stdev_surf_p(f,1) =  stdev_surf_p(f,1) &
                               +  mesh_fesom%areasvol(1, n) * stdev_p( sfields(f)%off + n)
        endif
     enddo ! n, myDim_nod2D
  enddo ! f, nfields
  
  ! Reduce to global mean
  call MPI_Allreduce (stdev_surf_p, stdev_surf_g, nfields*nlmax, MPI_DOUBLE_PRECISION, MPI_SUM, &
                      MPI_COMM_FESOM, MPIerr)
  do nz=1,nlmax
     stdev_surf_g(:,nz) = stdev_surf_g(:,nz) * inv_area_surf_glob(nz)
  enddo
  
  ! -----------------------------------------------------------------------------------------------------
  ! --- stdev_volo_g (nfields)    | global mean of grid-point ensemble STD for each field volume-weighted
  ! -----------------------------------------------------------------------------------------------------
  ! Compute pe-local mean of ensemble STD for each field
  if (mype_filter == 0) &
      write(*, *) 'FESOM-PDAF', '--- compute ensemble standard deviation global ocean mean'
  stdev_volo_p = 0.0
  stdev_volo_g = 0.0
  
  do f=1,nfields
     do n=1,myDim_nod2D
           if (sfields(f)%ndims == 2) then
           ! 3D fields
             do nz=1,nlmax
                stdev_volo_p(f) =  stdev_volo_p(f) &
                                +  cellvol(nz,n) * stdev_p( sfields(f)%off + (n-1)*(nlmax) + nz )
             enddo ! nz,nlmax
           endif
     enddo ! n,myDim_nod2D
  enddo ! f,nfields
  
  ! Reduce to global mean
  call MPI_Allreduce (stdev_volo_p, stdev_volo_g, nfields, MPI_DOUBLE_PRECISION, MPI_SUM, &
                      MPI_COMM_FESOM, MPIerr)
  do f=1,nfields
     if (sfields(f)%ndims == 2) then
     ! 3D fields
       stdev_volo_g(f) = stdev_volo_g(f) * inv_volo_full_glob
     else
     ! surface fields
       stdev_volo_g(f) = stdev_surf_g(f,1)
     endif
  enddo
  
  ! Global 3D mean of temperature field used to tune ensemble inflation
  stdevglob_temp = stdev_volo_g(id%temp)

  ! Display RMS errors
  if (mype_filter==0) then
     write (*,'(a, 10x,a)') &
          'FESOM-PDAF', 'Ensemble standard deviation:'
     write (*,'(a,7x,    a14,   a14,   a14,   a14,  a14, /a, 10x,70a)') &
          'FESOM-PDAF', 'CO2f','pCO2','temp','DIC','Alk', &
          'FESOM-PDAF', ('-',i=1,70)
     write (*,'(a,10x,  5es14.4, 3x,a13,a1,/a, 10x,70a)')  &
          'FESOM-PDAF', stdev_surf_g(id%CO2f  ,1), &
                        stdev_surf_g(id%pCO2s ,1), &
                        stdev_surf_g(id%temp  ,1), &
                        stdev_surf_g(id%DIC   ,1), &
                        stdev_surf_g(id%Alk   ,1), &
                       'surface STDEV', typestr, 'FESOM-PDAF', ('-',i=1,70)
     write (*,'(a,10x,  5es14.4, 3x,a13,a1,/a, 10x,70a)')  &
          'FESOM-PDAF', stdev_surf_g(id%CO2f  ,11), &
                        stdev_surf_g(id%pCO2s ,11), &
                        stdev_surf_g(id%temp  ,11), &
                        stdev_surf_g(id%DIC   ,11), &
                        stdev_surf_g(id%Alk   ,11), &
                       '90-100m STDEV', typestr, 'FESOM-PDAF', ('-',i=1,70)
     write (*,'(a,10x,  5es14.4, 3x,a13,a1,/a, 10x,70a)')  &
          'FESOM-PDAF', stdev_volo_g(id%CO2f ), &
                        stdev_volo_g(id%pCO2s), &
                        stdev_volo_g(id%temp ), &
                        stdev_volo_g(id%DIC  ), &
                        stdev_volo_g(id%Alk  ), &
                       'vol oce STDEV', typestr, 'FESOM-PDAF', ('-',i=1,70)
  end if
  
! *******************************
! *** Reset forgetting factor ***
! *******************************
! Forgetting factor increases after the start of the assimilation (days_since_DAstart).
! In case of model restarts:
!  -  days_since_DAstart is set by slurm-job-script
!  -  current forgetting factor and target temperature ensemble standard deviation are read from atmos-perturbation file
  
  if ((step-step_null)<0) then
  ! forecast phase   
     
     if (mype_filter==0) write(*,*) 'FESOM-PDAF', ' days_since_DAstart ', days_since_DAstart
     
     if (resetforget) then
     ! reset forgetting factor:
     ! set value for 1st half-month of assimilation
     if     (days_since_DAstart==  1) then
       forget=0.95
       call PDAF_reset_forget(forget)
       if (mype_filter==0) write(*,*) 'FESOM-PDAF', ' Resetting forget to ', forget, ' at day ', days_since_DAstart
     
      ! set value for 2nd half of 1st month of assimilation
     elseif (days_since_DAstart== 16) then
       forget=0.96
       call PDAF_reset_forget(forget)
       if (mype_filter==0) write(*,*) 'FESOM-PDAF', ' Resetting forget to ', forget, ' at day ', days_since_DAstart
     
      ! set value for 2nd month of assimilation 
     elseif (days_since_DAstart== 32) then
       forget=0.97
       call PDAF_reset_forget(forget)
       if (mype_filter==0) write(*,*) 'FESOM-PDAF', ' Resetting forget to ', forget, ' at day ', days_since_DAstart
     
     ! set value for 3rd month of assimilation 
     elseif (days_since_DAstart== 60) then
       forget=0.98
       call PDAF_reset_forget(forget)
       if (mype_filter==0) write(*,*) 'FESOM-PDAF', ' Resetting forget to ', forget, ' at day ', days_since_DAstart
     
     ! set value for 4th-17th months of assimilation
     elseif (days_since_DAstart== 90) then
       forget=0.99
       call PDAF_reset_forget(forget)
       if (mype_filter==0) write(*,*) 'FESOM-PDAF', ' Resetting forget to ', forget, ' at day ', days_since_DAstart
     
     ! during month 17, save temperature ensemble standard deviation ("stable_rmse")
     ! month 17
     elseif ((days_since_DAstart >= 485) .and. (days_since_DAstart <= 516)) then
       stable_rmse = stable_rmse + stdevglob_temp/31
       if (mype_filter==0) write(*,*) 'FESOM-PDAF', ' Saving ensemble spread to adapt forget at day ', days_since_DAstart
     
     ! after month 17, tune ensemble inflation based on the temperature ensemble standard deviation
     ! reset forgetting factor if ensemble standard deviation becomes larger/smaller than target value
     
     ! after month 17, higher ensemble standard deviation
     elseif ((days_since_DAstart >= 516) .and. (stdevglob_temp > 1.1*stable_rmse)) then
       forget=1.00
       call PDAF_reset_forget(forget)
       if (mype_filter==0) write(*,*) 'FESOM-PDAF '   , 'Resetting Forget ',&
                                      'Current RMSE: ', stdevglob_temp,&
                                      ' Target RMSE: ', stable_rmse, &
                                      ' New forget: ' , forget
                                      
     ! after month 17, lower ensemble standard deviation
     elseif ((days_since_DAstart >= 516) .and. (stdevglob_temp < 0.9*stable_rmse)) then
       forget=0.99
       call PDAF_reset_forget(forget)
       if (mype_filter==0) write(*,*) 'FESOM-PDAF '   , 'Resetting Forget ',&
                                      'Current RMSE: ', stdevglob_temp,&
                                      ' Target RMSE: ', stable_rmse, &
                                      ' New forget: ' , forget
     
     ! after month 17, ensemble standard deviation around target value
     else
       if (mype_filter==0) write(*,*) 'FESOM-PDAF '   , 'Keeping Forget ',&
                                      'Current RMSE: ', stdevglob_temp,&
                                      ' Target RMSE: ', stable_rmse, &
                                      ' Forget:      ', forget
     
     endif ! days_since_DAstart == ... (whether to reset forgetting factor according to scheme)
     endif ! resetforget               (whether to use reset scheme for forgetting factor)
     
     ! day count during daily forecast phase:
     days_since_DAstart=days_since_DAstart+1
     
  endif ! (forecast phase)


! ***************************************************************
! *** Compute statistics for effective observation dimensions ***
! ***************************************************************

  if (loctype==1 .and. (step-step_null)>0)  call adaptive_lradius_stats_pdaf()

  
! ***************************
! *** Compute daily means ***
! ***************************

  call update_means_pdaf(step, dim_p, state_p, stdev_p, delt_obs_ocn, write_now)


! **************************
! *** Write output files ***
! **************************
! note: after monthly output is written, reset monthly fields to zero

  ! *** write initial state fields ***
  if ((step - step_null)==0 .and. ( .not. this_is_pdaf_restart)) then
      ! ensemble mean
      if (w_dayensm) call netCDF_out('ii',state_p, int0, write_now, stdev_surf_g=stdev_surf_g, stdev_volo_g=stdev_volo_g)
      if (w_dayensm) call netCDF_out('si',stdev_p, int0, write_now)
      ! ensemble members
      if (w_daymemb) then
        do member = 1, dim_ens
           call netCDF_out('ii',ens_p(:,member), member, write_now)
        enddo
      endif
  endif

  ! daily output
  if (.not. write_now) then
     if ((step-step_null) < 0) then
        ! *** write forecast state fields ***
        ! during forecast phase, additionally, datetime and forgetting factor are written to output file
        ! ensemble mean
        if (w_dayensm) call netCDF_out('ff',state_p,int0, write_now, stdev_surf_g=stdev_surf_g, stdev_volo_g=stdev_volo_g, forget=forget)
        if (w_dayensm) call netCDF_out('sf',stdev_p,int0, write_now)
        ! ensemble members
        if (w_daymemb) then
           do member = 1, dim_ens
              call netCDF_out('ff',ens_p(:,member), member, write_now)
           enddo
        endif
     else if ((step-step_null) > 0) then
        ! *** write analysis ***
        ! ensemble mean
        if (w_dayensm) call netCDF_out('aa',state_p   , int0, write_now, stdev_surf_g=stdev_surf_g, stdev_volo_g=stdev_volo_g)
        if (w_dayensm) call netCDF_out('sa',stdev_p   , int0, write_now)
        ! ensemble members
        if (w_daymemb) then
           do member = 1, dim_ens
              call netCDF_out('aa',ens_p(:,member), member, write_now)
           enddo
        endif
     end if
  endif
  
  ! monthly output
  if (write_now) then
  ! end of month: pass monthly output in addition to daily output
     if ((step-step_null) < 0) then
        ! *** write forecast state fields ***
        ! ensemble mean, adding monthly mean of forecast states
        if (w_dayensm .or. w_monensm) call netCDF_out('ff',state_p, int0, write_now, stdev_surf_g=stdev_surf_g, stdev_volo_g=stdev_volo_g, forget=forget, m_state_p=monthly_state_f )
        if (w_dayensm .or. w_monensm) call netCDF_out('sf',stdev_p, int0, write_now,                                                                      m_state_p=monthly_stddev_f)
        ! ensemble members, adding snapshot
        if (w_daymemb .or. w_monmemb) then
           do member = 1, dim_ens
              call netCDF_out('ff',ens_p(:,member), member, write_now, m_state_p=ens_p(:,member))
           enddo
        endif
     else if ((step-step_null) > 0) then
        ! *** write analysis and "m"-state fields ***
        ! ensemble mean, adding monthly mean of analysis
        if (w_dayensm .or. w_monensm) call netCDF_out('aa',state_p   , int0, write_now, stdev_surf_g=stdev_surf_g, stdev_volo_g=stdev_volo_g, m_state_p=monthly_state_a )
        if (w_dayensm .or. w_monensm) call netCDF_out('sa',stdev_p   , int0, write_now,                                                       m_state_p=monthly_stddev_a)
! ensemble members, adding snapshot
        if (w_daymemb .or. w_monmemb) then
           do member = 1, dim_ens
              call netCDF_out('aa',ens_p(:,member), member, write_now, m_state_p=ens_p(:,member))
           enddo
        endif
     end if
  end if

  ! at last day of month, reset monthly_state to zero (has been written)

  if (write_now) call reset_means_pdaf(step)


! ********************
! *** finishing up ***
! ********************

  ! variables deallocated after analysis step
  if ((step-step_null) >= 0) then
     if (allocated(mean_O2_p ))       deallocate(mean_O2_p )
     if (allocated(mean_n_p))         deallocate(mean_n_p)
     if (allocated(mean_chl_cci_p))   deallocate(mean_chl_cci_p)
     if (allocated(mean_temp_p))      deallocate(mean_temp_p)
     if (allocated(mean_sss_cci_p))   deallocate(mean_sss_cci_p)
     if (allocated(mean_sss_p))       deallocate(mean_sss_p)     
     if (allocated(mean_sst_p))       deallocate(mean_sst_p)     
     if (allocated(mean_ice_p))       deallocate(mean_ice_p)     
  endif

  deallocate(stdev_p)

end subroutine prepoststep_pdaf
