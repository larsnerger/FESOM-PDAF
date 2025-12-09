!$Id: adaptive_lradius_pdaf.F90 2466 2021-02-25 12:46:34Z lnerger $
!>  Routine to compute adaptive localization radius
!!
!! The routine computes an adaptive localization
!! radius following Kirchgessner et al., MWR (2014).
!!
!! The routine is called by each filter process.
!!
!! __Revision history:__
!! * 2012    - Lars Nerger - Initial code for FESOM
!! * 2025-12 - Lars Nerger - Reviison for PDAF3
!! * other revisions - see repository log
!!
module adaptive_lradius_pdaf
contains

  subroutine get_adaptive_lradius_pdaf(thisobs, domain_p, lradius, loc_radius)

    use PDAF, &
         only: PDAF_local_weight, obs_f
    use assim_pdaf_mod, &           ! Variables for assimilation
         only: locweight, loctype, loc_ratio, dim_ens, eff_dim_obs
    use fesom_pdaf, &
         only: mesh_fesom, pi, mydim_nod2d, r_earth, r2g
    use parallel_pdaf_mod, &        ! Parallelization variables
         only: mype_filter

    implicit none

! *** Arguments ***
    type(obs_f), intent(in) :: thisobs   !< Observation type variable
    integer, intent(in)  :: domain_p     !< Current local analysis domain
    real, intent(out) :: lradius         !< Uniform localization radius
    real, intent(inout) :: loc_radius(mydim_nod2d)  !< Variable localization radius

! *** Local variables ***
    integer :: node, i, j         ! Counter
    real :: wc_coord(2)           ! Coordinates of current water column
    real :: o_coord(2)            ! Coordinates of observation
    real :: dist2d(2)             ! Distance vector between water column and observation
    real :: dist                  ! distance for adaptive localization
    real :: distance2, lradius2   ! squared distance and localization radius
    integer,save :: allocflag=0   ! Allocation flag
    integer,save :: first=1       ! Flag for first call
    integer,save :: domain_save=1 ! Stored domain index
    integer :: wtype              ! Type of weight function
    integer :: rtype              ! Type of weight regulation
    real :: l_range               ! Localization radius
    real :: weight                ! Weight for observation localization
    logical, save :: firstround = .true.   ! Flag for first analysis time
    real :: l_ranges(2)           ! Minimum and maximum allowed localization radii
    real :: l_step                ! Step size for search of localization radius
    real :: matA(1,1)             ! Temporary array for calling PDAF_local_weight

! *** Variable localization radius for fixed effective observation dimension ***

    ! Set parameters for localization radius search
    l_ranges(1) = 1.0e4
    l_ranges(2) = 1.0e7
    l_step     = 5.0e3

    if (mype_filter==0 .and. (domain_p<=domain_save .or. first==1)) then
       write (*,'(8x,a)') '--- Choose localization radius according to effective obs. dimension'
       write (*,'(8x,a,es10.2)') '--- ratio to effective obs. dimension', loc_ratio
       write (*,'(8x,a,2es10.2)') '--- minimum and maximum allowed radii', l_ranges
       first=0
    end if
    domain_save = domain_p

    if (locweight == 0) then
       ! Uniform (unit) weighting
       wtype = 0
       rtype = 0
    else if (locweight == 1 .or. locweight == 3) then
       ! Exponential weighting
       wtype = 1
       rtype = 0
    else if (locweight == 2 .or. locweight == 4 .or. locweight == 5 &
         .or. locweight == 6 .or. locweight == 7) then
       ! 5th-order polynomial (Gaspari&Cohn, 1999)
       wtype = 2
       rtype = 0 ! Always non-regulated here
    end if

    if (firstround) then

       ! *** At the first analysis time, compute the loclaization radius ***

       ! Get location of current water column (basis point)
       ! wc_coord is the rotated coordinates in FESOM
       call r2g(wc_coord(1), wc_coord(2), mesh_fesom%coord_nod2d(1, domain_p), mesh_fesom%coord_nod2d(2, domain_p))

       lrange: do l_range = l_ranges(1), l_ranges(2), l_step

          eff_dim_obs(domain_p) = 0

          ! Scan through full domain to initialize dimension and index array
          scanpointsB: do node = 1, thisobs%dim_obs_f

             ! location of observation point
             o_coord(1) = thisobs%ocoord_f(1, node)
             o_coord(2) = thisobs%ocoord_f(2, node)
              
             ! approximate distances in longitude and latitude
             dist2d(1) = r_earth * min( abs(wc_coord(1) - o_coord(1))* cos(wc_coord(2)), &
                  abs(abs(wc_coord(1) - o_coord(1)) - 2.0*pi) * cos(wc_coord(2)))
             dist2d(2) = r_earth * abs(wc_coord(2) - o_coord(2))

             ! full distance in meters
             dist = sqrt(dist2d(1)**2 + dist2d(2)**2)
           
             ! If distance below limit increment effective observation dimension
             if (dist <= l_range) then
                call PDAF_local_weight(wtype, rtype, l_range, l_range, dist, &
                     1, 1, matA, 1.0, weight, 0)

                eff_dim_obs(domain_p) = eff_dim_obs(domain_p) + weight
             end if

             if (eff_dim_obs(domain_p) >= loc_ratio * dim_ens) then
                loc_radius(domain_p) = l_range
                exit lrange
             end if

          end do scanpointsB

       end do lrange

       if (domain_p==mydim_nod2d) firstround = .false.

    end if

    lradius = loc_radius(domain_p)

  end subroutine get_adaptive_lradius_pdaf

!-------------------------------------------------------
!>  Routine to compute statistics about adaptive localization radius
!!
!! The routine computes the minimum, maximum,
!! and mean values for the adaptive localization
!! radius following Kirchgessner et al., MWR (2014).
!!
!! The routine is called by each filter process.
!!
!! __Revision history:__
!! 2012 - Lars Nerger - Initial code for FESOM
!! * Later revisions - see repository log
!!
  subroutine adaptive_lradius_stats_pdaf() 

    use parallel_pdaf_mod, &
         only: mype_filter, npes=>npes_filter, &
         comm => COMM_filter, &
         MPI_REAL8, MPIerr, &
         MPI_INTEGER, MPI_SUM, MPI_MAX, MPI_MIN 
    use assim_pdaf_mod, &       ! Variables for assimilation
         only: eff_dim_obs, loctype
    use fesom_pdaf, &
         only: mesh_fesom, dim_p => mydim_nod2d

    implicit none

! *** Local variables ***
    integer :: i                                   ! Counters
    real :: min_eff_dim_obs, max_eff_dim_obs       ! Stats on effective observation dimensions
    real :: min_eff_dim_obs_g, max_eff_dim_obs_g   ! Stats on effective observation dimensions
    real :: sum_eff_dim_obs, avg_eff_dim_obs_g     ! Stats on effective observation dimensions


! ***************************************************************
! *** Compute statistics for effective observation dimensions ***
! ***************************************************************

    if (loctype==1) then

       max_eff_dim_obs = 0.0
       min_eff_dim_obs = 1.0e16
       sum_eff_dim_obs = 0.0

       do i = 1, dim_p
          if (eff_dim_obs(i) > max_eff_dim_obs) max_eff_dim_obs = eff_dim_obs(i)
          if (eff_dim_obs(i) < min_eff_dim_obs) min_eff_dim_obs = eff_dim_obs(i)
          sum_eff_dim_obs = sum_eff_dim_obs + eff_dim_obs(i)
       end do
       if (npes>1) then
          call MPI_Reduce(sum_eff_dim_obs, avg_eff_dim_obs_g, 1, MPI_REAL8, MPI_SUM, &
               0, COMM, MPIerr)
          call MPI_Reduce(max_eff_dim_obs, max_eff_dim_obs_g, 1, MPI_REAL8, MPI_MAX, &
               0, COMM, MPIerr)
          call MPI_Reduce(min_eff_dim_obs, min_eff_dim_obs_g, 1, MPI_REAL8, MPI_MIN, &
               0, COMM, MPIerr)
       else
          ! This is a work around for working with nullmpi.F90
          avg_eff_dim_obs_g = sum_eff_dim_obs
          min_eff_dim_obs_g = min_eff_dim_obs
          max_eff_dim_obs_g = max_eff_dim_obs
       end if

       if (mype_filter==0) then
          avg_eff_dim_obs_g = avg_eff_dim_obs_g / real(mesh_fesom%nod2D)

          write (*, '(8x, a)') &
               '--- Effective observation dimensions for local analysis:'
          write (*, '(12x, a, f12.2)') &
               'min. effective observation dimension:       ', min_eff_dim_obs_g
          write (*, '(12x, a, f12.2)') &
               'max. effective observation dimension:       ', max_eff_dim_obs_g
          write (*, '(12x, a, f12.2)') &
               'avg. effective observation dimension:       ', avg_eff_dim_obs_g
       end if
    end if

  end subroutine adaptive_lradius_stats_pdaf

end module adaptive_lradius_pdaf
