MODULE output_config

! USES:

  IMPLICIT NONE

! *** Local Variables ***

  character(len=200) :: filename_std = ''       ! Full name of output file

  LOGICAL :: w_daymemb    = .false.       ! whether to write any daily ensemble member states
  LOGICAL :: w_dayensm    = .false.       ! whether to write any daily ensemble mean states
  LOGICAL :: w_monmemb    = .false.       ! whether to write any monthly ensemble member states
  LOGICAL :: w_monensm    = .false.       ! whether to write any monthly ensemble mean states
  LOGICAL :: w_mm         = .false.       ! whether to write any m-fields (day-averages)
  LOGICAL :: w_sm         = .false.       ! whether to write any m-fields (day-averages) of standard deviation

! Field description:

! Indices of whether to write: 
  INTEGER, PARAMETER :: ff=1, aa=2, mm=3, ii=4 ! forecast (ff), analysis (aa), mean (mm) and initial (ii)
  INTEGER, PARAMETER :: sf=5, sa=6, si=7, sm=8 ! ensemble standard deviation snapshots: forecast (sf), analysis (sa), initial (si) and mean (sm)
  INTEGER, PARAMETER :: oo=1, ee=2, dd=3       ! any output (oo), ensemble members (ee) and daily values (dd)

  LOGICAL :: setoutput(18)
                                             
CONTAINS

  SUBROUTINE configure_output()

    USE mod_assim_pdaf, &
         ONLY: assimilateBGC, assimilatePHY, &
         cda_phy, cda_bio, nlmax, mydim_nod2d
    USE statevector_pdaf, &
         only: id, sfields, nfields, &
         phymin, phymax, bgcmin, bgcmax
    USE mod_parallel_pdaf, &
         ONLY: writepe, mype_world
    USE obs_chl_cci_pdafomi,    ONLY: assim_o_chl_cci
    USE obs_DIC_glodap_pdafomi, ONLY: assim_o_DIC_glodap
    USE obs_Alk_glodap_pdafomi, ONLY: assim_o_Alk_glodap
    USE obs_pCO2_SOCAT_pdafomi, ONLY: assim_o_pCO2_SOCAT
    USE obs_o2_comf_pdafomi,    ONLY: assim_o_o2_comf
    USE obs_n_comf_pdafomi,     ONLY: assim_o_n_comf
    USE obs_o2_argo_pdafomi,    ONLY: assim_o_o2_argo
    USE obs_n_argo_pdafomi,     ONLY: assim_o_n_argo
    USE obs_o2_merged_pdafomi,  ONLY: assim_o_o2_merged
    USE obs_sss_smos_pdafomi,   ONLY: assim_o_sss
    USE obs_sss_cci_pdafomi,    ONLY: assim_o_sss_cci
    USE obs_ssh_cmems_pdafomi,  ONLY: assim_o_ssh 
    USE obs_sst_pdafomi,        ONLY: assim_o_sst 
    USE obs_TSprof_EN4_pdafomi, ONLY: assim_o_en4_s
    USE obs_TSprof_EN4_pdafomi, ONLY: assim_o_en4_t

    IMPLICIT NONE

! *** Local variables: ***
    INTEGER :: i, b,p,s,j                      ! Counter
    CHARACTER(len=12) :: outputmessage(8,3)


! *********************************************************
! ***   Set field-specific output type and frequency    ***
! *********************************************************

! [oo] True   - write output
!      False  - no output
! [ee] True   - write ensemble members
!      False  - write ensemble mean
! [dd] True   - daily
!      False  - monthly

! Defaults: False
! activate one (or multiple if not contradictory - think!) of the following predefined output schemes:


! *** write daily forecast and analysis ensemble members  ***
    IF (setoutput(1)) THEN
       DO s=1, nfields
          ! forecast
          sfields(s)% output(ff,oo) = .True.
          sfields(s)% output(ff,ee) = .True.
          sfields(s)% output(ff,dd) = .True.
          ! analysis
          sfields(s)% output(aa,oo) = .True.
          sfields(s)% output(aa,ee) = .True.
          sfields(s)% output(aa,dd) = .True.
       ENDDO
    ENDIF

! *** write daily forecast and analysis ensemble mean  ***
    IF (setoutput(2)) THEN
       DO s=1, nfields
          ! forecast
          sfields(s)% output(ff,oo) = .True.
          sfields(s)% output(ff,dd) = .True.
          ! analysis
          sfields(s)% output(aa,oo) = .True.
          sfields(s)% output(aa,dd) = .True.
       ENDDO
    ENDIF

! *** write monthly forecast and analysis ensemble mean of updated variables ***
    IF (setoutput(3)) THEN
       DO s=1, nfields
          ! forecast
          IF (sfields(s)%updated)   sfields(s)% output(ff,oo) = .True.
          ! analysis
          IF (sfields(s)%updated)   sfields(s)% output(aa,oo) = .True.
       ENDDO
    ENDIF

! *** write daily m-fields ensemble mean  ***
    IF (setoutput(4)) THEN
       DO s=1, nfields
          sfields(s)% output(mm,oo) = .True.
          sfields(s)% output(mm,dd) = .True.
       ENDDO
    ENDIF

! *** write initial fields ensemble mean  ***
    IF (setoutput(5)) THEN
       DO s=1, nfields
          sfields(s)% output(ii,oo) = .True.
       ENDDO
    ENDIF

! *** write initial fields ensemble members  ***
    IF (setoutput(6)) THEN
       DO s=1, nfields
          sfields(s)% output(ii,oo) = .True.
          sfields(s)% output(ii,ee) = .True.
       ENDDO
    ENDIF

! *** write monthly m-fields ensemble mean ***
    IF (setoutput(7)) THEN
       DO s=1, nfields
          sfields(s)% output(mm,oo) = .True.
       ENDDO
    ENDIF

! *** write daily m-fields of assimilated variables ensemble mean ***
    IF (setoutput(8)) THEN
       ! activate m-field output
       IF (assim_o_sst)          sfields(id% temp)   % output(mm,oo) = .True.
       IF (assim_o_sss)          sfields(id% salt)   % output(mm,oo) = .True.
       IF (assim_o_sss_cci)      sfields(id% salt)   % output(mm,oo) = .True.
       IF (assim_o_en4_t)        sfields(id% temp)   % output(mm,oo) = .True.
       IF (assim_o_en4_s)        sfields(id% salt)   % output(mm,oo) = .True.
       IF (assim_o_ssh)          sfields(id% SSH)    % output(mm,oo) = .True.
       IF (assim_o_chl_cci)      sfields(id% PhyChl) % output(mm,oo) = .True.
       IF (assim_o_chl_cci)      sfields(id% DiaChl) % output(mm,oo) = .True.
       IF (assim_o_DIC_glodap)   sfields(id% DIC)    % output(mm,oo) = .True.
       IF (assim_o_Alk_glodap)   sfields(id% Alk)    % output(mm,oo) = .True.
       IF (assim_o_pCO2_SOCAT)   sfields(id% pCO2s)  % output(mm,oo) = .True.
       IF (assim_o_o2_comf)      sfields(id% O2)     % output(mm,oo) = .True.
       IF (assim_o_n_comf)       sfields(id% DIN)    % output(mm,oo) = .True.
       IF (assim_o_o2_argo)      sfields(id% O2)     % output(mm,oo) = .True.
       IF (assim_o_n_argo)       sfields(id% DIN)    % output(mm,oo) = .True.
       IF (assim_o_o2_merged)    sfields(id% O2)     % output(mm,oo) = .True.

       ! set to daily
       IF (assim_o_sst)          sfields(id% temp)   % output(mm,dd) = .True.
       IF (assim_o_sss)          sfields(id% salt)   % output(mm,dd) = .True.
       IF (assim_o_sss_cci)      sfields(id% salt)   % output(mm,dd) = .True.
       IF (assim_o_en4_t)        sfields(id% temp)   % output(mm,dd) = .True.
       IF (assim_o_en4_s)        sfields(id% salt)   % output(mm,dd) = .True.
       IF (assim_o_ssh)          sfields(id% SSH)    % output(mm,dd) = .True.
       IF (assim_o_chl_cci)      sfields(id% PhyChl) % output(mm,dd) = .True.
       IF (assim_o_chl_cci)      sfields(id% DiaChl) % output(mm,dd) = .True.
       IF (assim_o_DIC_glodap)   sfields(id% DIC)    % output(mm,dd) = .True.
       IF (assim_o_Alk_glodap)   sfields(id% Alk)    % output(mm,dd) = .True.
       IF (assim_o_pCO2_SOCAT)   sfields(id% pCO2s)  % output(mm,dd) = .True.
       IF (assim_o_o2_comf)      sfields(id% O2)     % output(mm,dd) = .True.
       IF (assim_o_n_comf)       sfields(id% DIN)    % output(mm,dd) = .True.
       IF (assim_o_o2_argo)      sfields(id% O2)     % output(mm,dd) = .True.
       IF (assim_o_n_argo)       sfields(id% DIN)    % output(mm,dd) = .True.
       IF (assim_o_o2_merged)    sfields(id% O2)     % output(mm,dd) = .True.
    ENDIF

! *** write daily m-fields of assimilate-able BGC variables ensemble mean ***
    IF (setoutput(9)) THEN
       ! activate m-field output
       sfields(id% PhyChl) % output(mm,oo) = .True.
       sfields(id% DiaChl) % output(mm,oo) = .True.
       sfields(id% DIC)    % output(mm,oo) = .True.
       sfields(id% Alk)    % output(mm,oo) = .True.
       sfields(id% pCO2s)  % output(mm,oo) = .True.
       sfields(id% O2)     % output(mm,oo) = .True.
       sfields(id% DIN)    % output(mm,oo) = .True.

       ! set to daily
       sfields(id% PhyChl) % output(mm,dd) = .True.
       sfields(id% DiaChl) % output(mm,dd) = .True.
       sfields(id% DIC)    % output(mm,dd) = .True.
       sfields(id% Alk)    % output(mm,dd) = .True.
       sfields(id% pCO2s)  % output(mm,dd) = .True.
       sfields(id% O2)     % output(mm,dd) = .True.
       sfields(id% DIN)    % output(mm,dd) = .True.
    ENDIF
  
! *** write daily m-fields of CO2 flux and pCO2 ensemble mean ***
    IF (setoutput(10)) THEN
       sfields(id% CO2f)  % output(mm,oo) = .True.
       sfields(id% pCO2s) % output(mm,oo) = .True.
       sfields(id% CO2f)  % output(mm,dd) = .True.
       sfields(id% pCO2s) % output(mm,dd) = .True.
    ENDIF

! *** write daily m-fields of variables defining the CO2 flux ***
    IF (setoutput(11)) THEN
       sfields(id% CO2f     ) % output(mm,oo) = .True.
       sfields(id% pCO2s    ) % output(mm,oo) = .True.
       sfields(id% a_ice    ) % output(mm,oo) = .True.
       sfields(id% alphaCO2 ) % output(mm,oo) = .True.
       sfields(id% PistonVel) % output(mm,oo) = .True.

       sfields(id% CO2f     ) % output(mm,dd) = .True.
       sfields(id% pCO2s    ) % output(mm,dd) = .True.
       sfields(id% a_ice    ) % output(mm,dd) = .True.
       sfields(id% alphaCO2 ) % output(mm,dd) = .True.
       sfields(id% PistonVel) % output(mm,dd) = .True.
    ENDIF

! *** write initial fields standard deviation  ***
    IF (setoutput(12)) THEN
       DO s=1, nfields
          sfields(s)% output(si,oo) = .True.
       ENDDO
    ENDIF

! *** write daily forecast of standard deviation ***
    IF (setoutput(13)) THEN
       DO s=1, nfields
          sfields(s)% output(sf,oo) = .True.
          sfields(s)% output(sf,dd) = .True.
       ENDDO
    ENDIF

! *** write monthly forecast of standard deviation ***
    IF (setoutput(14)) THEN
       DO s=1, nfields
          sfields(s)% output(sf,oo) = .True.
       ENDDO
    ENDIF

! *** write monthly forecast and analysis of standard deviation for updated fields ***
    IF (setoutput(15)) THEN
       DO s=1, nfields
          IF (sfields(s)%updated)   sfields(s)% output(sf,oo) = .True.
          IF (sfields(s)%updated)   sfields(s)% output(sa,oo) = .True.
       ENDDO
    ENDIF

! *** write monthly mm-fields of standard deviation ***
    IF (setoutput(16)) THEN
       DO s=1, nfields
          sfields(s)% output(sm,oo) = .True.
       ENDDO
    ENDIF

! *** write daily forecast fields for oxygen and DIN         ***
! (to reconstruct inno_omit)
    IF (setoutput(17)) THEN
       ! Oxygen
       sfields(id% O2) % output(ff,oo) = .True. ! activate
       sfields(id% O2) % output(ff,dd) = .True. ! daily
       ! DIN
       sfields(id% DIN) % output(ff,oo) = .True. ! activate
       sfields(id% DIN) % output(ff,dd) = .True. ! daily
    ENDIF

! *** write daily analysis of assimilate-able BGC variables ensemble mean ***
    IF (setoutput(18)) THEN
       ! activate m-field output
       sfields(id% PhyChl) % output(aa,oo) = .True.
       sfields(id% DiaChl) % output(aa,oo) = .True.
       sfields(id% DIC)    % output(aa,oo) = .True.
       sfields(id% Alk)    % output(aa,oo) = .True.
       sfields(id% pCO2s)  % output(aa,oo) = .True.
       sfields(id% O2)     % output(aa,oo) = .True.
       sfields(id% DIN)    % output(aa,oo) = .True.

       ! set to daily
       sfields(id% PhyChl) % output(aa,dd) = .True.
       sfields(id% DiaChl) % output(aa,dd) = .True.
       sfields(id% DIC)    % output(aa,dd) = .True.
       sfields(id% Alk)    % output(aa,dd) = .True.
       sfields(id% pCO2s)  % output(aa,dd) = .True.
       sfields(id% O2)     % output(aa,dd) = .True.
       sfields(id% DIN)    % output(aa,dd) = .True.
    ENDIF


! *** FINALIZE        ***

    DO s=1, nfields
       ! no monthly initial fields. monthly fields are written at last day of month, but initial fields at first day.
       ! because monthly initial fields make no sense, we set "daily" (dd=True) for all initial fields.
       sfields(s)% output(ii,dd) = .True.
       sfields(s)% output(si,dd) = .True.
  
       ! do not compute full ensemble state for m-fields! this takes memory and time. simply use fesom-output instead.
       ! whatever settings made be before, we reset m-fields ens-member output to False, in the end.
       sfields(s)% output(mm,ee) = .False.
  
       ! the standard deviation is written to ensemble-mean file
       sfields(s)% output(si,ee) = .False.
       sfields(s)% output(sf,ee) = .False.
       sfields(s)% output(sa,ee) = .False.
       sfields(s)% output(sm,ee) = .False.
    ENDDO


! *** which files are needed?           ***_

    w_daymemb = .false.          
    w_dayensm = .false.          
    w_monmemb = .false.          
    w_monensm = .false.          

    DO s=1, nfields
  !                               -- output? -----------        -- ensemble members? ----------       -- daily? ----------------------
       w_daymemb = w_daymemb .or. any( sfields(s)%output(:,oo) .and. (      sfields(s)%output(:,ee)) .and. (      sfields(s)%output(:,dd)))  ! have any daily member states to write?
       w_dayensm = w_dayensm .or. any( sfields(s)%output(:,oo) .and. (.not. sfields(s)%output(:,ee)) .and. (      sfields(s)%output(:,dd)))  ! have any daily ensemble mean states to write?
       w_monmemb = w_monmemb .or. any( sfields(s)%output(:,oo) .and. (      sfields(s)%output(:,ee)) .and. (.not. sfields(s)%output(:,dd)))  ! have any monthly member states to write?
       w_monensm = w_monensm .or. any( sfields(s)%output(:,oo) .and. (.not. sfields(s)%output(:,ee)) .and. (.not. sfields(s)%output(:,dd)))  ! have any monthly ensemble mean states to write?
    ENDDO

    w_dayensm = .true. ! daily file to protocol forgetting factor and global standard deviation


    w_mm = any(sfields(:)%output(mm,oo))    ! any m-fields (day-mean) to be written?
    w_sm = any(sfields(:)%output(sm,oo))    ! any m-fields (day-mean) of standard deviation to be written?

    ! output message
    IF (mype_world==0) THEN

       DO s=1, nfields
          outputmessage(aa,1) = 'Analysis'
          outputmessage(ff,1) = 'Forecast'
          outputmessage(mm,1) = 'M-Field'
          outputmessage(ii,1) = 'Initial'
          outputmessage(sf,1) = 'STD Forcast'
          outputmessage(sa,1) = 'STD Analysis'
          outputmessage(si,1) = 'STD Initial'
          outputmessage(sm,1) = 'STD M-Field'

          DO j=1,8
             IF (sfields(s)%output(j,ee)) then
                outputmessage(j,ee) = 'Ens-Memb'
             ELSE
                outputmessage(j,ee) = 'Ens-Mean'
             ENDIF
             IF (sfields(s)%output(j,dd)) then
                outputmessage(j,dd) = 'Daily'
             ELSE
                outputmessage(j,dd) = 'Monthly'
             ENDIF
          ENDDO ! j=1,4

          outputmessage(ii,dd) = '-'
          outputmessage(si,dd) = '-'

          if (.not. any(sfields(s)%output(:,oo))) then ! this field any output?
             write (*, '(a,4x,a,1x,a10,1x,a9)') 'FESOM-PDAF', 'Field', sfields(s)%variable, 'No Output'
          else
             DO j=1,4
                if (sfields(s)%output(j,oo)) write (*, '(a,4x,a,1x,a10,1x,a12,1x,a10,1x,a10,1x,a10)') &
                     'FESOM-PDAF', 'Field', &
                     sfields(s)%variable, &
                     outputmessage(j, 1), &
                     outputmessage(j,ee), &
                     outputmessage(j,dd)
             ENDDO ! j=1,8
          endif ! this field any output?
       ENDDO ! s=1, nfields
    ENDIF ! writepe

  END SUBROUTINE configure_output
  
END MODULE output_config
