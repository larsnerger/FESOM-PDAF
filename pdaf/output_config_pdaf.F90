module output_config_pdaf

! USES:

  implicit none

! *** Local Variables ***

  character(len=200) :: filename_std = ''       ! Full name of output file

  logical :: w_daymemb    = .false.       ! whether to write any daily ensemble member states
  logical :: w_dayensm    = .false.       ! whether to write any daily ensemble mean states
  logical :: w_monmemb    = .false.       ! whether to write any monthly ensemble member states
  logical :: w_monensm    = .false.       ! whether to write any monthly ensemble mean states
  logical :: w_mm         = .false.       ! whether to write any m-fields (day-averages)
  logical :: w_sm         = .false.       ! whether to write any m-fields (day-averages) of standard deviation

! Field description:

! Indices of whether to write: 
  integer, parameter :: ff=1, aa=2, mm=3, ii=4 ! forecast (ff), analysis (aa), mean (mm) and initial (ii)
  integer, parameter :: sf=5, sa=6, si=7, sm=8 ! ensemble standard deviation snapshots: forecast (sf), analysis (sa), initial (si) and mean (sm)
  integer, parameter :: oo=1, ee=2, dd=3       ! any output (oo), ensemble members (ee) and daily values (dd)

  logical :: setoutput(18) = .false.
  ! Setoutput will be set ny the namelist 'pdafoutput' (or in init_pdaf)
  !  setoutput( 1)     ! daily forecast and analysis ensemble members
  !  setoutput( 2)     ! daily forecast and analysis ensemble mean
  !  setoutput( 3)     ! monthly forecast and analysis ensemble mean of updated variables
  !  setoutput( 4)     ! daily m-fields ensemble mean
  !  setoutput( 5)     ! initial fields ensemble mean
  !  setoutput( 6)     ! initial fields ensemble members
  !  setoutput( 7)     ! monthly m-fields ensemble mean
  !  setoutput( 8)     ! daily m-fields of assimilated variables ensemble mean
  !  setoutput( 9)     ! daily m-fields of assimilate-able BGC variables ensemble mean
  !  setoutput(10)     ! daily m-fields of CO2 flux and pCO2 ensemble mean
  !  setoutput(11)     ! daily m-fields of variables defining the CO2 flux
  !  setoutput(12)     ! initial fields standard deviation
  !  setoutput(13)     ! daily forecast of standard deviation
  !  setoutput(14)     ! monthly forecast of standard deviation
  !  setoutput(15)     ! monthly forecast and analysis of standard deviation for updated fields
  !  setoutput(16)     ! monthly mm-fields of standard deviation
  !  setoutput(17)     ! if COMF-O2 assimilated: daily O2 forecast
                                             
contains

  subroutine configure_output()

    use fesom_pdaf, &
         only: nlmax, mydim_nod2d
    use statevector_pdaf, &
         only: id, sfields, nfields, &
         phymin, phymax, bgcmin, bgcmax
    use parallel_pdaf_mod, &
         only: writepe, mype_world
    use obs_chl_cci_pdafomi,    only: assim_o_chl_cci
    use obs_DIC_glodap_pdafomi, only: assim_o_DIC_glodap
    use obs_Alk_glodap_pdafomi, only: assim_o_Alk_glodap
    use obs_pCO2_SOCAT_pdafomi, only: assim_o_pCO2_SOCAT
    use obs_o2_comf_pdafomi,    only: assim_o_o2_comf
    use obs_n_comf_pdafomi,     only: assim_o_n_comf
    use obs_o2_argo_pdafomi,    only: assim_o_o2_argo
    use obs_n_argo_pdafomi,     only: assim_o_n_argo
    use obs_o2_merged_pdafomi,  only: assim_o_o2_merged
    use obs_sss_smos_pdafomi,   only: assim_o_sss
    use obs_sss_cci_pdafomi,    only: assim_o_sss_cci
    use obs_ssh_cmems_pdafomi,  only: assim_o_ssh 
    use obs_sst_pdafomi,        only: assim_o_sst 
    use obs_TSprof_EN4_pdafomi, only: assim_o_en4_s
    use obs_TSprof_EN4_pdafomi, only: assim_o_en4_t

    implicit none

! *** Local variables: ***
    integer :: i, b,p,s,j                      ! Counter
    character(len=12) :: outputmessage(8,3)


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
    if (setoutput(1)) then
       do s=1, nfields
          ! forecast
          sfields(s)% output(ff,oo) = .true.
          sfields(s)% output(ff,ee) = .true.
          sfields(s)% output(ff,dd) = .true.
          ! analysis
          sfields(s)% output(aa,oo) = .true.
          sfields(s)% output(aa,ee) = .true.
          sfields(s)% output(aa,dd) = .true.
       enddo
    endif

! *** write daily forecast and analysis ensemble mean  ***
    if (setoutput(2)) then
       do s=1, nfields
          ! forecast
          sfields(s)% output(ff,oo) = .true.
          sfields(s)% output(ff,dd) = .true.
          ! analysis
          sfields(s)% output(aa,oo) = .true.
          sfields(s)% output(aa,dd) = .true.
       enddo
    endif

! *** write monthly forecast and analysis ensemble mean of updated variables ***
    if (setoutput(3)) then
       do s=1, nfields
          ! forecast
          if (sfields(s)%updated)   sfields(s)% output(ff,oo) = .true.
          ! analysis
          if (sfields(s)%updated)   sfields(s)% output(aa,oo) = .true.
       enddo
    endif

! *** write daily m-fields ensemble mean  ***
    if (setoutput(4)) then
       do s=1, nfields
          sfields(s)% output(mm,oo) = .true.
          sfields(s)% output(mm,dd) = .true.
       enddo
    endif

! *** write initial fields ensemble mean  ***
    if (setoutput(5)) then
       do s=1, nfields
          sfields(s)% output(ii,oo) = .true.
       enddo
    endif

! *** write initial fields ensemble members  ***
    if (setoutput(6)) then
       do s=1, nfields
          sfields(s)% output(ii,oo) = .true.
          sfields(s)% output(ii,ee) = .true.
       enddo
    endif

! *** write monthly m-fields ensemble mean ***
    if (setoutput(7)) then
       do s=1, nfields
          sfields(s)% output(mm,oo) = .true.
       enddo
    endif

! *** write daily m-fields of assimilated variables ensemble mean ***
    if (setoutput(8)) then
       ! activate m-field output
       if (assim_o_sst)          sfields(id% temp)   % output(mm,oo) = .true.
       if (assim_o_sss)          sfields(id% salt)   % output(mm,oo) = .true.
       if (assim_o_sss_cci)      sfields(id% salt)   % output(mm,oo) = .true.
       if (assim_o_en4_t)        sfields(id% temp)   % output(mm,oo) = .true.
       if (assim_o_en4_s)        sfields(id% salt)   % output(mm,oo) = .true.
       if (assim_o_ssh)          sfields(id% SSH)    % output(mm,oo) = .true.
       if (assim_o_chl_cci)      sfields(id% PhyChl) % output(mm,oo) = .true.
       if (assim_o_chl_cci)      sfields(id% DiaChl) % output(mm,oo) = .true.
       if (assim_o_DIC_glodap)   sfields(id% DIC)    % output(mm,oo) = .true.
       if (assim_o_Alk_glodap)   sfields(id% Alk)    % output(mm,oo) = .true.
       if (assim_o_pCO2_SOCAT)   sfields(id% pCO2s)  % output(mm,oo) = .true.
       if (assim_o_o2_comf)      sfields(id% O2)     % output(mm,oo) = .true.
       if (assim_o_n_comf)       sfields(id% DIN)    % output(mm,oo) = .true.
       if (assim_o_o2_argo)      sfields(id% O2)     % output(mm,oo) = .true.
       if (assim_o_n_argo)       sfields(id% DIN)    % output(mm,oo) = .true.
       if (assim_o_o2_merged)    sfields(id% O2)     % output(mm,oo) = .true.

       ! set to daily
       if (assim_o_sst)          sfields(id% temp)   % output(mm,dd) = .true.
       if (assim_o_sss)          sfields(id% salt)   % output(mm,dd) = .true.
       if (assim_o_sss_cci)      sfields(id% salt)   % output(mm,dd) = .true.
       if (assim_o_en4_t)        sfields(id% temp)   % output(mm,dd) = .true.
       if (assim_o_en4_s)        sfields(id% salt)   % output(mm,dd) = .true.
       if (assim_o_ssh)          sfields(id% SSH)    % output(mm,dd) = .true.
       if (assim_o_chl_cci)      sfields(id% PhyChl) % output(mm,dd) = .true.
       if (assim_o_chl_cci)      sfields(id% DiaChl) % output(mm,dd) = .true.
       if (assim_o_DIC_glodap)   sfields(id% DIC)    % output(mm,dd) = .true.
       if (assim_o_Alk_glodap)   sfields(id% Alk)    % output(mm,dd) = .true.
       if (assim_o_pCO2_SOCAT)   sfields(id% pCO2s)  % output(mm,dd) = .true.
       if (assim_o_o2_comf)      sfields(id% O2)     % output(mm,dd) = .true.
       if (assim_o_n_comf)       sfields(id% DIN)    % output(mm,dd) = .true.
       if (assim_o_o2_argo)      sfields(id% O2)     % output(mm,dd) = .true.
       if (assim_o_n_argo)       sfields(id% DIN)    % output(mm,dd) = .true.
       if (assim_o_o2_merged)    sfields(id% O2)     % output(mm,dd) = .true.
    endif

! *** write daily m-fields of assimilate-able BGC variables ensemble mean ***
    if (setoutput(9)) then
       ! activate m-field output
       sfields(id% PhyChl) % output(mm,oo) = .true.
       sfields(id% DiaChl) % output(mm,oo) = .true.
       sfields(id% DIC)    % output(mm,oo) = .true.
       sfields(id% Alk)    % output(mm,oo) = .true.
       sfields(id% pCO2s)  % output(mm,oo) = .true.
       sfields(id% O2)     % output(mm,oo) = .true.
       sfields(id% DIN)    % output(mm,oo) = .true.

       ! set to daily
       sfields(id% PhyChl) % output(mm,dd) = .true.
       sfields(id% DiaChl) % output(mm,dd) = .true.
       sfields(id% DIC)    % output(mm,dd) = .true.
       sfields(id% Alk)    % output(mm,dd) = .true.
       sfields(id% pCO2s)  % output(mm,dd) = .true.
       sfields(id% O2)     % output(mm,dd) = .true.
       sfields(id% DIN)    % output(mm,dd) = .true.
    endif
  
! *** write daily m-fields of CO2 flux and pCO2 ensemble mean ***
    if (setoutput(10)) then
       sfields(id% CO2f)  % output(mm,oo) = .true.
       sfields(id% pCO2s) % output(mm,oo) = .true.
       sfields(id% CO2f)  % output(mm,dd) = .true.
       sfields(id% pCO2s) % output(mm,dd) = .true.
    endif

! *** write daily m-fields of variables defining the CO2 flux ***
    if (setoutput(11)) then
       sfields(id% CO2f     ) % output(mm,oo) = .true.
       sfields(id% pCO2s    ) % output(mm,oo) = .true.
       sfields(id% a_ice    ) % output(mm,oo) = .true.
       sfields(id% alphaCO2 ) % output(mm,oo) = .true.
       sfields(id% PistonVel) % output(mm,oo) = .true.

       sfields(id% CO2f     ) % output(mm,dd) = .true.
       sfields(id% pCO2s    ) % output(mm,dd) = .true.
       sfields(id% a_ice    ) % output(mm,dd) = .true.
       sfields(id% alphaCO2 ) % output(mm,dd) = .true.
       sfields(id% PistonVel) % output(mm,dd) = .true.
    endif

! *** write initial fields standard deviation  ***
    if (setoutput(12)) then
       do s=1, nfields
          sfields(s)% output(si,oo) = .true.
       enddo
    endif

! *** write daily forecast of standard deviation ***
    if (setoutput(13)) then
       do s=1, nfields
          sfields(s)% output(sf,oo) = .true.
          sfields(s)% output(sf,dd) = .true.
       enddo
    endif

! *** write monthly forecast of standard deviation ***
    if (setoutput(14)) then
       do s=1, nfields
          sfields(s)% output(sf,oo) = .true.
       enddo
    endif

! *** write monthly forecast and analysis of standard deviation for updated fields ***
    if (setoutput(15)) then
       do s=1, nfields
          if (sfields(s)%updated)   sfields(s)% output(sf,oo) = .true.
          if (sfields(s)%updated)   sfields(s)% output(sa,oo) = .true.
       enddo
    endif

! *** write monthly mm-fields of standard deviation ***
    if (setoutput(16)) then
       do s=1, nfields
          sfields(s)% output(sm,oo) = .true.
       enddo
    endif

! *** write daily forecast fields for oxygen and DIN         ***
! (to reconstruct inno_omit)
    if (setoutput(17)) then
       ! Oxygen
       sfields(id% O2) % output(ff,oo) = .true. ! activate
       sfields(id% O2) % output(ff,dd) = .true. ! daily
       ! DIN
       sfields(id% DIN) % output(ff,oo) = .true. ! activate
       sfields(id% DIN) % output(ff,dd) = .true. ! daily
    endif

! *** write daily analysis of assimilate-able BGC variables ensemble mean ***
    if (setoutput(18)) then
       ! activate m-field output
       sfields(id% PhyChl) % output(aa,oo) = .true.
       sfields(id% DiaChl) % output(aa,oo) = .true.
       sfields(id% DIC)    % output(aa,oo) = .true.
       sfields(id% Alk)    % output(aa,oo) = .true.
       sfields(id% pCO2s)  % output(aa,oo) = .true.
       sfields(id% O2)     % output(aa,oo) = .true.
       sfields(id% DIN)    % output(aa,oo) = .true.

       ! set to daily
       sfields(id% PhyChl) % output(aa,dd) = .true.
       sfields(id% DiaChl) % output(aa,dd) = .true.
       sfields(id% DIC)    % output(aa,dd) = .true.
       sfields(id% Alk)    % output(aa,dd) = .true.
       sfields(id% pCO2s)  % output(aa,dd) = .true.
       sfields(id% O2)     % output(aa,dd) = .true.
       sfields(id% DIN)    % output(aa,dd) = .true.
    endif


! *** FINALIZE        ***

    do s=1, nfields
       ! no monthly initial fields. monthly fields are written at last day of month, but initial fields at first day.
       ! because monthly initial fields make no sense, we set "daily" (dd=True) for all initial fields.
       sfields(s)% output(ii,dd) = .true.
       sfields(s)% output(si,dd) = .true.
  
       ! do not compute full ensemble state for m-fields! this takes memory and time. simply use fesom-output instead.
       ! whatever settings made be before, we reset m-fields ens-member output to False, in the end.
       sfields(s)% output(mm,ee) = .false.
  
       ! the standard deviation is written to ensemble-mean file
       sfields(s)% output(si,ee) = .false.
       sfields(s)% output(sf,ee) = .false.
       sfields(s)% output(sa,ee) = .false.
       sfields(s)% output(sm,ee) = .false.
    enddo


! *** which files are needed?           ***_

    w_daymemb = .false.          
    w_dayensm = .false.          
    w_monmemb = .false.          
    w_monensm = .false.          

    do s=1, nfields
  !                               -- output? -----------        -- ensemble members? ----------       -- daily? ----------------------
       w_daymemb = w_daymemb .or. any( sfields(s)%output(:,oo) .and. (      sfields(s)%output(:,ee)) .and. (      sfields(s)%output(:,dd)))  ! have any daily member states to write?
       w_dayensm = w_dayensm .or. any( sfields(s)%output(:,oo) .and. (.not. sfields(s)%output(:,ee)) .and. (      sfields(s)%output(:,dd)))  ! have any daily ensemble mean states to write?
       w_monmemb = w_monmemb .or. any( sfields(s)%output(:,oo) .and. (      sfields(s)%output(:,ee)) .and. (.not. sfields(s)%output(:,dd)))  ! have any monthly member states to write?
       w_monensm = w_monensm .or. any( sfields(s)%output(:,oo) .and. (.not. sfields(s)%output(:,ee)) .and. (.not. sfields(s)%output(:,dd)))  ! have any monthly ensemble mean states to write?
    enddo

    w_dayensm = .true. ! daily file to protocol forgetting factor and global standard deviation


    w_mm = any(sfields(:)%output(mm,oo))    ! any m-fields (day-mean) to be written?
    w_sm = any(sfields(:)%output(sm,oo))    ! any m-fields (day-mean) of standard deviation to be written?

    ! output message
    if (mype_world==0) then

       write (*,'(/a,2x,a)') 'FESOM-PDAF', '*** Configuration of file output ***'

       do s=1, nfields
          outputmessage(aa,1) = 'Analysis'
          outputmessage(ff,1) = 'Forecast'
          outputmessage(mm,1) = 'M-Field'
          outputmessage(ii,1) = 'Initial'
          outputmessage(sf,1) = 'STD Forcast'
          outputmessage(sa,1) = 'STD Analysis'
          outputmessage(si,1) = 'STD Initial'
          outputmessage(sm,1) = 'STD M-Field'

          do j=1,8
             if (sfields(s)%output(j,ee)) then
                outputmessage(j,ee) = 'Ens-Memb'
             else
                outputmessage(j,ee) = 'Ens-Mean'
             endif
             if (sfields(s)%output(j,dd)) then
                outputmessage(j,dd) = 'Daily'
             else
                outputmessage(j,dd) = 'Monthly'
             endif
          enddo ! j=1,4

          outputmessage(ii,dd) = '-'
          outputmessage(si,dd) = '-'

          if (.not. any(sfields(s)%output(:,oo))) then ! this field any output?
             write (*, '(a,4x,a,1x,a10,1x,a9)') 'FESOM-PDAF', 'Field', sfields(s)%variable, 'No Output'
          else
             do j=1,4
                if (sfields(s)%output(j,oo)) write (*, '(a,4x,a,1x,a10,1x,a12,1x,a10,1x,a10,1x,a10)') &
                     'FESOM-PDAF', 'Field', &
                     sfields(s)%variable, &
                     outputmessage(j, 1), &
                     outputmessage(j,ee), &
                     outputmessage(j,dd)
             enddo ! j=1,8
          endif ! this field any output?
       enddo ! s=1, nfields
    endif ! writepe

  end subroutine configure_output
  
end module output_config_pdaf
