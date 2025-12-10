! PDAF in AWI-CM2 / Fesom 2.0

subroutine read_config_pdaf()

! !DESCRIPTION:
! This routine read the namelist file with
! parameters controlling data assimilation with
! PDAF.
! 2019-11 - Longjiang Mu - Initial commit for AWI-CM3

! !USES:


! general assimilation settings
  use parallel_pdaf_mod, &
       only: mype_model, n_modeltasks, task_id
  use assim_pdaf_mod, &
       only: dim_state, dim_state_p, dim_ens, dim_lag, &
       screen, filtertype, subtype, &
       delt_obs_ocn, &
       dim_bias, DA_couple_type, &
       type_forget, &
       forget, locweight, cradius, sradius, &
       type_trans, type_sqrt, step_null, &
       eff_dim_obs, loc_radius, loctype, loc_ratio, &
       path_init, file_init, file_inistate, read_inistate, &
       twin_experiment, dim_obs_max, use_global_obs, DAoutput_path, &
       this_is_pdaf_restart, &
       path_atm_cov, days_since_DAstart, assimilateBGC, assimilatePHY, &
       start_from_ENS_spinup, resetforget, &
       ! initial ensemble perturbation
       varscale, perturb_ssh, perturb_u, &
       perturb_v, perturb_temp, perturb_salt, &
       perturb_DIC, perturb_Alk, perturb_DIN, perturb_O2, &
       ! Temp-Salt-Profiles:
       path_obs_rawprof, file_rawprof_prefix, file_rawprof_suffix, &
       proffiles_o, start_year_o, end_year_o
  use output_config_pdaf, &       ! Output
       only: setoutput
  use mod_postprocess, &     ! Postprocessing
       only: isPP, yearPP, pathsim
       ! clock
  use g_clock, &
       only: yearold, yearnew
  
  
  ! physics observations
  use obs_sst_pdafomi, &
       only: assim_o_sst, rms_obs_sst, path_obs_sst, file_sst_prefix, file_sst_suffix, &
       sst_exclude_ice, sst_exclude_diff, bias_obs_sst, sst_fixed_rmse
  use obs_sss_smos_pdafomi, &
       only: assim_o_sss, rms_obs_sss, path_obs_sss, file_sss_prefix, file_sss_suffix, &
       sss_exclude_ice, sss_exclude_diff, bias_obs_sss, sss_fixed_rmse
  use obs_sss_cci_pdafomi, &
       only: assim_o_sss_cci, rms_obs_sss_cci, path_obs_sss_cci, file_sss_cci_prefix, file_sss_cci_suffix, &
       sss_cci_exclude_ice, sss_cci_exclude_diff, bias_obs_sss_cci, sss_cci_fixed_rmse
  use obs_ssh_cmems_pdafomi, &
       only: assim_o_ssh, rms_obs_ssh, path_obs_ssh, file_ssh_prefix, file_ssh_suffix, &
       ssh_exclude_ice, ssh_exclude_diff, bias_obs_ssh, ssh_fixed_rmse, lradius_ssh
  use obs_TSprof_EN4_pdafomi, &
       only: ASSIM_O_en4_t, ASSIM_O_en4_s, & 
       path_obs_prof, file_prof_prefix, file_prof_suffix, &
       rms_obs_S, rms_obs_T, &
       file_syntobs_prof, prof_exclude_diff, bias_obs_prof

  ! biogeochem observations    
  use obs_chl_cci_pdafomi, &
       only: assim_o_chl_cci, rms_obs_chl_cci, path_obs_chl_cci, file_chl_cci_prefix, file_chl_cci_suffix, &
       chl_cci_exclude_ice, chl_cci_exclude_diff, bias_obs_chl_cci, chl_cci_fixed_rmse, chl_logarithmic, &
       path_bias, file_bias_prefix
  use obs_DIC_glodap_pdafomi, &
       only: assim_o_DIC_glodap, path_obs_DIC_glodap, &
       rms_obs_DIC_glodap, DIC_glodap_exclude_diff
  use obs_Alk_glodap_pdafomi, &
       only: assim_o_Alk_glodap, path_obs_Alk_glodap, &
       rms_obs_Alk_glodap, Alk_glodap_exclude_diff
  use obs_pco2_SOCAT_pdafomi, &
       only: assim_o_pCO2_SOCAT, path_obs_pCO2_SOCAT, &
       rms_obs_pCO2_SOCAT, pCO2_SOCAT_exclude_diff
  use obs_o2_comf_pdafomi, &
       only: assim_o_o2_comf, path_obs_o2_comf, &
       rms_obs_o2_comf, o2_comf_exclude_diff
  use obs_n_comf_pdafomi, &
       only: assim_o_n_comf, path_obs_n_comf, &
       rms_obs_n_comf, n_comf_exclude_diff
  use obs_o2_argo_pdafomi, &
       only: assim_o_o2_argo, path_obs_o2_argo, &
       rms_obs_o2_argo, o2_argo_exclude_diff
  use obs_N_argo_pdafomi, &
       only: assim_o_N_argo, path_obs_N_argo, &
       rms_obs_N_argo, N_argo_exclude_diff
  use obs_o2_merged_pdafomi, &
       only: assim_o_o2_merged, path_obs_o2_merged, &
       rms_obs_o2_merged, o2_merged_exclude_diff, &
       o2_merged_excl_absolute, o2_merged_excl_relative
  use obs_N_merged_pdafomi, &
       only: assim_o_N_merged, path_obs_N_merged, &
       rms_obs_N_merged, N_merged_exclude_diff, &
       n_merged_excl_absolute, n_merged_excl_relative
  
  ! ensemble initialization    
  use mod_atmos_ens_stochasticity, &
       only: disturb_xwind, disturb_ywind, disturb_humi, &
       disturb_qlw, disturb_qsr, disturb_tair, &
       disturb_prec, disturb_snow, disturb_mslp, &
       atmos_stochasticity_ON, write_atmos_st, &
       varscale_wind, varscale_tair, &
       varscale_humi, varscale_qlw
  use mod_perturbation_pdaf, &
       only: perturb_scale, &
       perturb_params_bio, perturb_params_phy
  


  implicit none


! Local variables
  character(len=100) :: nmlfile ='namelist.fesom.pdaf'    ! name of namelist file
  character(len=32)  :: handle             ! Handle for command line parser
  logical :: printconfig = .true.          ! Print information on all configuration parameters
  character(len=4) :: year_string          ! Current year as string
  
       
  namelist /pdaf/ filtertype, subtype, screen, &
       type_forget, forget, resetforget, dim_bias, &
       cradius, locweight, sradius, DA_couple_type, &
       n_modeltasks, use_global_obs, &
       path_init, file_init, step_null, printconfig, &
       file_inistate, read_inistate, &
       type_trans, type_sqrt, dim_lag, &
       loctype, loc_ratio, delt_obs_ocn, &     
       dim_obs_max, &
       twin_experiment, &
       DAoutput_path, &
       this_is_pdaf_restart, &
       days_since_DAstart, start_from_ENS_spinup, &
       assimilatePHY, &
       ! initial ensemble perturbation:
       varscale, perturb_ssh, perturb_u, &
       perturb_v, perturb_temp, perturb_salt, &
       perturb_DIC, perturb_Alk, perturb_DIN, perturb_O2, &
       ! Salt SMOS:
       ASSIM_o_sss, path_obs_sss, file_sss_prefix, file_sss_suffix, &
       rms_obs_sss, sss_fixed_rmse, &
       sss_exclude_ice, sss_exclude_diff, &
       ! Salt CCI:
       ASSIM_o_sss_cci, path_obs_sss_cci, file_sss_cci_prefix, file_sss_cci_suffix, &
       rms_obs_sss_cci, sss_cci_fixed_rmse, &
       sss_cci_exclude_ice, sss_cci_exclude_diff, &
       ! SSH:
       ASSIM_o_ssh, path_obs_ssh, file_ssh_prefix, file_ssh_suffix, &
       rms_obs_ssh, ssh_fixed_rmse, bias_obs_ssh, lradius_ssh, &
       ! SST:
       ASSIM_o_sst, path_obs_sst, file_sst_prefix, file_sst_suffix, &
       rms_obs_sst, sst_fixed_rmse, bias_obs_sst, &
       sst_exclude_ice, sst_exclude_diff, &
       ! Profiles:
       ASSIM_o_en4_t, ASSIM_o_en4_S, &
       path_obs_prof, file_prof_prefix, file_prof_suffix, &
       rms_obs_S, rms_obs_T, &
       path_obs_rawprof, file_rawprof_prefix, proffiles_o, &
       start_year_o, end_year_o, &
       ! parameter perturbation:
       perturb_scale, perturb_params_bio, perturb_params_phy, &
       ! BioGeoChemistry:
       assimilateBGC, &
       ! Chl-a CCI:
       assim_o_chl_cci, rms_obs_chl_cci, path_obs_chl_cci, &
       file_chl_cci_prefix, file_chl_cci_suffix, &
       chl_cci_exclude_ice, chl_cci_exclude_diff, &
       bias_obs_chl_cci, chl_cci_fixed_rmse, chl_logarithmic, &
       path_bias, file_bias_prefix, &
       ! DIC GLODAP:
       assim_o_DIC_glodap, path_obs_DIC_glodap, &
       rms_obs_DIC_glodap, DIC_glodap_exclude_diff, &
       ! Alk GLODAP:
       assim_o_Alk_glodap, path_obs_Alk_glodap, &
       rms_obs_Alk_glodap, Alk_glodap_exclude_diff, &
       ! pCO2 SOCAT:
       assim_o_pCO2_SOCAT, path_obs_pCO2_SOCAT, &
       rms_obs_pCO2_SOCAT, pCO2_SOCAT_exclude_diff, &
       ! O2 COMFORT:
       assim_o_o2_comf, path_obs_o2_comf, &
       rms_obs_o2_comf, o2_comf_exclude_diff, &
       ! DIN COMFORT:
       assim_o_n_comf, path_obs_n_comf, &
       rms_obs_n_comf, n_comf_exclude_diff, &
       ! O2 ARGO:
       assim_o_o2_argo, path_obs_o2_argo, &
       rms_obs_o2_argo, o2_argo_exclude_diff, &
       ! DIN ARGO:
       assim_o_N_argo, path_obs_N_argo, &
       rms_obs_N_argo, N_argo_exclude_diff, &
       ! O2 MERGED:
       assim_o_o2_merged, path_obs_o2_merged, &
       rms_obs_o2_merged, o2_merged_exclude_diff, &
       o2_merged_excl_absolute, o2_merged_excl_relative, &
       ! DIN MERGED:
       assim_o_n_merged, path_obs_n_merged, &
       rms_obs_n_merged, n_merged_exclude_diff, &
       n_merged_excl_absolute, n_merged_excl_relative

  namelist /atmos_stoch/ &
       path_atm_cov, &
       disturb_xwind, disturb_ywind, disturb_humi, &
       disturb_qlw, disturb_qsr, disturb_tair, &
       disturb_prec, disturb_snow, disturb_mslp, &
       varscale_wind, varscale_tair, varscale_humi, varscale_qlw, &
       write_atmos_st
       
  namelist /pp/ &
       isPP, yearPP, pathsim
       
  namelist /pdafoutput/ &
       setoutput
              
! ****************************************************
! ***   Initialize PDAF parameters from namelist   ***
! ****************************************************

! *** Read namelist file ***
#ifdef DEBUG
  write(*,*) 'Read PDAF namelist file: ',nmlfile
#endif

  open (20,file=nmlfile)
  read (20,NML =pdaf)
  close(20)

  open (30,file=nmlfile)
  read (30,NML =atmos_stoch)
  close(30)
  
  open (40,file=nmlfile)
  read (40,NML =pp)
  close(40)
  
  open (50,file=nmlfile)
  read (50,NML =pdafoutput)
  close(50)

! *** Add trailing slash to paths ***
  call add_slash(path_obs_sst)
  call add_slash(path_obs_sss)
  call add_slash(path_obs_ssh)
  
  call add_slash(path_obs_prof)
  call add_slash(path_init)
  
  call add_slash(path_obs_chl_cci)
  call add_slash(path_obs_DIC_glodap)
  call add_slash(path_obs_Alk_glodap)
  call add_slash(path_obs_pCO2_SOCAT)
  call add_slash(path_obs_O2_comf)
  call add_slash(path_obs_N_comf)
  
! *** Is atmospheric stochasticity used at all?
if (disturb_humi   .or. &
    disturb_mslp   .or. &
    disturb_xwind  .or. &
    disturb_ywind  .or. &
    disturb_qlw    .or. &
    disturb_qsr    .or. &
    disturb_tair   .or. &
    disturb_prec   .or. &
    disturb_snow   ) then
    
    atmos_stochasticity_ON = .true.
else
    atmos_stochasticity_ON = .false.
endif

! Observation file prefixes:
write(year_string,'(i4.4)') yearnew

file_sst_prefix = 'OSTIA_SST_'//trim(year_string)//'0101_'//trim(year_string)//'1231_daily_dist72_'
file_sss_prefix = 'SMOS_SSS_'//trim(year_string)//'_dist72_'
file_sss_cci_prefix = 'CCI_SSS_'//trim(year_string)//'_dist72_'
file_chl_cci_prefix = 'CCI_OC_'//trim(year_string)//'_dist72_'

! *** Print configuration variables ***
  showconf: if (printconfig .and. mype_model==0 .and. task_id==1) then

     ! Overview of PDAF configuration
     write (*,'(/a,1x,a)')           'FESOM-PDAF',   '-- Overview of PDAF configuration --'
     write (*,'(a,3x,a)')            'FESOM-PDAF',   'PDAF [namelist: pdaf]:'
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'filtertype  ',         filtertype
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'subtype     ',         subtype
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'n_modeltasks',         n_modeltasks
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'dim_ens     ',         dim_ens
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'delt_obs_ocn',         delt_obs_ocn
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'step_null   ',         step_null
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'days_since_DAstart',   days_since_DAstart
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'screen      ',         screen
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'type_forget ',         type_forget
     write (*,'(a,5x,a20,1x,f10.4)') 'FESOM-PDAF',   'forget      ',         forget
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'resetforget ',         resetforget
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'dim_bias    ',         dim_bias
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'type_trans  ',         type_trans
     write (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'cradius     ',         cradius
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'locweight   ',         locweight
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'loctype     ',         loctype
     write (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'sradius     ',         sradius
     write (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'loc_ratio   ',         loc_ratio
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'use_global_obs',       use_global_obs
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'dim_lag     ',         dim_lag
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'DA_couple_type  ',     DA_couple_type
     
     ! Pre-processing of TS-profile data
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'proffiles_o  ',        proffiles_o
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'start_year_o ',        start_year_o
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'end_year_o   ',        end_year_o
     
     ! Physics observation-type settings
     write (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'bias_obs_sst  ',       bias_obs_sst
     write (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'bias_obs_prof ',       bias_obs_prof
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'sst_exclude_ice',      sst_exclude_ice
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'sst_exclude_diff',     sst_exclude_diff
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'prof_exclude_diff',    prof_exclude_diff
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_sst     ',    trim(path_obs_sst)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_sst_prefix  ',    trim(file_sst_prefix)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_sst_suffix  ',    trim(file_sst_suffix)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_prof    ',    trim(path_obs_prof)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_prof_prefix ',    trim(file_prof_prefix)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_prof_suffix ',    trim(file_prof_suffix)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_rawprof    ', trim(path_obs_rawprof)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_rawprof_prefix ', trim(file_rawprof_prefix)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_rawprof_suffix ', trim(file_rawprof_suffix)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'DAoutput_path ',       trim(DAoutput_path)
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_sss   ',       assim_o_sss
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_sss_cci',      assim_o_sss_cci
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_ssh   ',       assim_o_ssh
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_sst   ',       assim_o_sst
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_en4_t ',       assim_o_en4_t
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_en4_s ',       assim_o_en4_s
     write (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'rms_obs_sst ',         rms_obs_sst
     write (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'rms_obs_sss ',         rms_obs_sss
     write (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'rms_obs_sss_cci ',     rms_obs_sss_cci
     write (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'rms_obs_ssh ',         rms_obs_ssh
     write (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'lradius_ssh ',         lradius_ssh
     write (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'rms_obs_T   ',         rms_obs_T
     write (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'rms_obs_S   ',         rms_obs_S
     write (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'bias_obs_ssh',         bias_obs_ssh
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'ssh_fixed_rmse',       ssh_fixed_rmse
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'sst_fixed_rmse',       sst_fixed_rmse
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'sss_fixed_rmse',       sss_fixed_rmse
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'sss_cci_fixed_rmse',   sss_cci_fixed_rmse
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_sss     ',    trim(path_obs_sss)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_sss_prefix  ',    trim(file_sss_prefix)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_sss_suffix  ',    trim(file_sss_suffix)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_sss_cci ',    trim(path_obs_sss_cci)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_sss_cci_prefix',  trim(file_sss_cci_prefix)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_sss_cci_suffix',  trim(file_sss_cci_suffix)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_ssh     ',    trim(path_obs_ssh)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_ssh_prefix  ',    trim(file_ssh_prefix)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_ssh_suffix  ',    trim(file_ssh_suffix)
     
     ! BGC observation-type settings
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_chl_cci',      assim_o_chl_cci
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_chl_cci_prefix',  trim(file_chl_cci_prefix)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_chl_cci_suffix',  trim(file_chl_cci_suffix)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_chl_cci ',    trim(path_obs_chl_cci)
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'chl_cci_exclude_ice',  chl_cci_exclude_ice
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'chl_cci_fixed_rmse',   chl_cci_fixed_rmse
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_chl_cci',      rms_obs_chl_cci
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'chl_cci_exclude_diff', chl_cci_exclude_diff
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'bias_obs_chl_cci',     bias_obs_chl_cci
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'chl_logarithmic',      chl_logarithmic
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_bias',            trim(path_bias)
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_bias_prefix',     trim(file_bias_prefix)
     
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_DIC_glodap',      assim_o_DIC_glodap
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_DIC_glodap',     path_obs_DIC_glodap
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_DIC_glodap',      rms_obs_DIC_glodap
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'DIC_glodap_exclude_diff', DIC_glodap_exclude_diff
     
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_Alk_glodap',      assim_o_Alk_glodap
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_Alk_glodap',     path_obs_Alk_glodap
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_Alk_glodap',      rms_obs_Alk_glodap
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'Alk_glodap_exclude_diff', Alk_glodap_exclude_diff
     
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_pCO2_SOCAT',      assim_o_pCO2_SOCAT
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_pCO2_SOCAT',     path_obs_pCO2_SOCAT
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_pCO2_SOCAT',      rms_obs_pCO2_SOCAT
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'pCO2_SOCAT_exclude_diff', pCO2_SOCAT_exclude_diff
     
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_n_comf',       assim_o_n_comf
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_n_comf',      path_obs_n_comf
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_n_comf',       rms_obs_n_comf
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'N_comf_exclude_diff',  N_comf_exclude_diff
     
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_o2_comf',      assim_o_o2_comf
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_o2_comf',     path_obs_o2_comf
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_o2_comf',      rms_obs_o2_comf
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'o2_comf_exclude_diff', o2_comf_exclude_diff
     
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_o2_argo',      assim_o_o2_argo
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_o2_argo',     path_obs_o2_argo
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_o2_argo',      rms_obs_o2_argo
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'o2_argo_exclude_diff', o2_argo_exclude_diff
     
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_N_argo',       assim_o_N_argo
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_N_argo',      path_obs_N_argo
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_N_argo',       rms_obs_N_argo
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'N_argo_exclude_diff',  N_argo_exclude_diff
     
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_o2_merged',       assim_o_o2_merged
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_o2_merged',      path_obs_o2_merged
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_o2_merged',       rms_obs_o2_merged
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'o2_merged_exclude_diff',  o2_merged_exclude_diff
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'o2_merged_excl_absolute', o2_merged_excl_absolute
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'o2_merged_excl_relative', o2_merged_excl_relative
     
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_N_merged',       assim_o_n_merged
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_N_merged',      path_obs_n_merged
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_N_merged',       rms_obs_n_merged
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'N_merged_exclude_diff',  n_merged_exclude_diff
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'N_merged_excl_absolute', n_merged_excl_absolute
     write (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'N_merged_excl_relative', n_merged_excl_relative
          
     ! initial ensemble covariance
     write (*,'(a,5x,a20,1x,1x,a)')     'FESOM-PDAF',   'path_init   ',         trim(path_init)
     write (*,'(a,5x,a20,1x,1x,a)')     'FESOM-PDAF',   'file_init   ',         trim(file_init)
     write (*,'(a,5x,a20,1x,1x,l)')     'FESOM-PDAF',   'this_is_pdaf_restart', this_is_pdaf_restart
     write (*,'(a,5x,a20,1x,1x,l)')     'FESOM-PDAF',   'start_from_ENS_spinup',start_from_ENS_spinup
     write (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'varscale    ', varscale
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_ssh',  perturb_ssh
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_u',    perturb_u
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_v',    perturb_v
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_temp', perturb_temp
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_salt', perturb_salt
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_DIC',  perturb_DIC
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_Alk',  perturb_Alk
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_DIN',  perturb_DIN
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_O2',   perturb_O2
     
     ! Twin experiment
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'twin_experiment', twin_experiment
     if (filtertype==100 .or. twin_experiment) then
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'dim_obs_max ',    dim_obs_max
     end if
     if (read_inistate) then
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_inistate ',  trim(file_inistate)
     endif

     ! Atmospheric perturbation
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_atm_cov  ',       trim(path_atm_cov)
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'atmos_stochasticity_ON',atmos_stochasticity_ON
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'write_atmos_st',       write_atmos_st
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'disturb_xwind',        disturb_xwind
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'disturb_ywind',        disturb_ywind
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'disturb_humi',         disturb_humi
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'disturb_qlw',          disturb_qlw
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'disturb_qsr',          disturb_qsr
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'disturb_tair',         disturb_tair
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'disturb_prec',         disturb_prec
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'disturb_snow',         disturb_snow
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'disturb_mslp',         disturb_mslp
     write (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'varscale_wind ',       varscale_wind
     write (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'varscale_tair ',       varscale_tair
     write (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'varscale_humi ',       varscale_humi
     write (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'varscale_qlw  ',       varscale_qlw
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_params_bio',   perturb_params_bio
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_params_phy',   perturb_params_phy
     write (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'perturb_scale',        perturb_scale
     
     ! Updated state variables
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assimilatePHY',        assimilatePHY
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assimilateBGC',        assimilateBGC
     
     ! Postprocessing
     write (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'isPP',                 isPP
     write (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'yearPP',               yearPP
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'pathsim(1)',           trim(pathsim(1))
     write (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'pathsim(2)',           trim(pathsim(2))

     write (*,'(a,1x,a)') 'FESOM-PDAF','-- End of PDAF configuration overview --'
     

  end if showconf

end subroutine read_config_pdaf



! ==============================================================================
!
! !ROUTINE: add_slash --- Add trailing slash to path string
!
! !INTERFACE:
subroutine add_slash(path)

! !DESCRIPTION:
! This routine ensures that a string defining a path
! has a trailing slash.
!
! !USES:
  implicit none

! !ARGUMENTS:
  character(len=100) :: path  ! String holding the path

! *** Local variables ***
  integer :: strlength

! *** Add trailing slash ***
  strlength = len_trim(path)

  if (path(strlength:strlength) /= '/') then
     path = trim(path) // '/'
  end if

end subroutine add_slash
