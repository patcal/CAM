module CAPT_mod
!
! CAPT_mod: Module containing CAPT processing routines.
!           The entire set of routines in the original cade are
!           included. Some are currently non-functional for use 
!           here, but are included as a starting point for future 
!           development. (e.g.the mass_fixer here requires global 
!           summs and needs to be modified for the MPI environment 
!           before it can be used).
!======================================================================
  ! Useful modules
  !------------------
  use shr_kind_mod,  only: r8=>SHR_KIND_R8, cs=>SHR_KIND_CS, cl=>SHR_KIND_CL
  use cam_abortutils,only: endrun
  use cam_logfile,   only: iulog

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public :: mask_ice
  public :: mask_sst
  public :: jra_25_press_full_levels
  public :: binning
  public :: cubic_opt1
  public :: cubic_slav
  public :: extx
  public :: extys
  public :: extyv
  public :: lcdbas
  public :: lcdint
  public :: tsadj
  private:: tsadj_rect
  private:: tsadj_se
  public :: vert_quad_opt1
  private:: vert_quad_opt1_rect
  private:: vert_quad_opt1_se
  public :: vert_int_opt1
  private:: vert_int_opt1_rect
  private:: vert_int_opt1_se
  public :: vert_int_opt2
  private:: vert_int_opt2_rect
  private:: vert_int_opt2_se
  public :: myminmax
  private:: myminmax_rect
  private:: myminmax_se
  public :: psadj
  private:: psadj_rect
  private:: psadj_se
  public :: q2rh
  private:: q2rh_rect
  private:: q2rh_se
  public :: rh2q
  private:: rh2q_rect
  private:: rh2q_se
  public :: esinti
  public :: gestbl
  public :: gffgch
  public :: mass_fixer
  private:: mass_fixer_rect
  private:: mass_fixer_se

  ! Interfaces
  !-------------
  interface tsadj
    procedure:: tsadj_rect
    procedure:: tsadj_se
  end interface

  interface vert_quad_opt1
    procedure:: vert_quad_opt1_rect
    procedure:: vert_quad_opt1_se
  end interface

  interface vert_int_opt1
    procedure:: vert_int_opt1_rect
    procedure:: vert_int_opt1_se
  end interface

  interface vert_int_opt2
    procedure:: vert_int_opt2_rect
    procedure:: vert_int_opt2_se
  end interface

  interface myminmax
    procedure:: myminmax_rect
    procedure:: myminmax_se
  end interface

  interface psadj
    procedure:: psadj_rect
    procedure:: psadj_se
  end interface

  interface q2rh
    procedure:: q2rh_rect
    procedure:: q2rh_se
  end interface

  interface rh2q
    procedure:: rh2q_rect
    procedure:: rh2q_se
  end interface

  interface mass_fixer
    procedure:: mass_fixer_rect
    procedure:: mass_fixer_se
  end interface

  ! Global Values
  !---------------
  integer,parameter:: plenest=250
  type comes_t
    real(r8):: estbl(plenest)
    real(r8):: tmin
    real(r8):: tmax
    real(r8):: ttrice
    real(r8):: pcf(6)
    logical :: icephs
  end type comes_t

  type(comes_t):: COM

contains
  !=====================================================================
  subroutine mask_ice(plat, plon, ice, landfrac)
    !
    ! Use land mask to overwrite land fraction with ice fraction to 
    ! prepare for interpolation of ice fraction to new resolution
    !===================================================================
    implicit none
    ! 
    ! Passed variables 
    !------------------
    integer :: plat                 ! latitude  dimension of input/output fields
    integer :: plon                 ! longitude dimension of input/output fields
    real(r8):: ice      (plon,plat) ! input/output ice fraction
    real(r8):: landfrac (plon,plat) ! land mask
    !
    ! Local values
    !--------------
    integer :: i,ii,j
    integer :: icount
    real(r8):: tw                   ! on lat circle, ocn point west of land
    integer :: iw                   ! index of western ocn point
    integer :: iw_tmp               ! index of western ocn point
    real(r8):: te                   ! on lat circle, ocn point east of land
    integer :: ie                   ! index of eastern ocn point
    integer :: ie_tmp               ! index of eastern ocn point
    integer :: offshore_pts         ! how many points offshore to sample ice fraction
    real(r8):: land_threshold       ! minimum fraction to be considered land
    real(r8):: wgt_w, wgt_e         ! interpolation weights
    logical :: found_ocn_w          ! true when ocn point is found on lat 
                                    ! circle west of current land point
    logical :: found_ocn_e          ! true when ocn point is found on lat 
                                    ! circle east of current land point
    logical :: found_lnd            ! true when lnd point is found

    land_threshold = 0._r8
    offshore_pts   = 3

    if((offshore_pts .lt.1).or.(offshore_pts.gt.10)) then
      write(iulog,*) 'number of points offshore (in east/west direction)'
      write(iulog,*) 'to sample ice must be at least 1 but not more than 10'
      write(iulog,*) 'offshore_pts = ', offshore_pts
      call endrun
    endif

    ! Overwrite land fraction with interpolated ice fraction
    !-------------------------------------------------------
    do j = 1, plat
      found_ocn_w = .false.
      found_ocn_e = .false.
      iw          = -999
      ie          = -999
      do i = 1,plon
        if(landfrac(i,j).le.land_threshold) then
          ! If found ocean, set tw/iw and move on to next point eastward
          !-------------------------------------------------------------
          found_ocn_w = .true.
          iw          = i
        else
          ! current point is land and will be overwritten
          !
          ! If "i=1", search westward over entire latitude circle for first
          ! ocean point west of current land point.
          !-----------------------------------------------------------------
          if(i .eq. 1) then
            found_ocn_w = .false.
            ii = plon + 1
            do while((ii .gt.1).and.(.not.found_ocn_w))
              ii = ii - 1
              if(landfrac(ii,j).le.land_threshold) then
                iw          = ii-plon
                found_ocn_w = .true.
              endif
            end do
          endif

          ! Now search eastward over entire latitude circle for first 
          ! ocean point encountered east of current land point
          !-------------------------------------------------------------
          if((ie.lt.i).and.(found_ocn_w)) then
            found_ocn_e = .false.
            ii = i
            ie = i
            do while((ie.le.plon+i).and.(.not.found_ocn_e))
              ie = ie + 1
              ii = ii + 1
              if(ii .gt. plon) ii = 1
              if(landfrac(ii,j).le.land_threshold) then
                found_ocn_e = .true.
              endif
            end do
          endif

          if((.not.found_ocn_w).and.(.not.found_ocn_e)) then
            ! If no ocean points were found in this latitude circle, assume 
            ! we are at Antarctica and set all ice fractions to 1.
            !--------------------------------------------------------------
            ice(i,j) = 1.
          elseif((found_ocn_w).and.(found_ocn_e)) then
            ! Else, interpolate ice fraction over land point
            !------------------------------------------------
            iw_tmp = iw
            ie_tmp = ie

            ! First, decide how many points offshore (in both directions) 
            ! to sample ice fraction for interpolation over land
            !--------------------------------------------------------------
            if(offshore_pts .gt. 1) then
              ! Eastern shore:  look east to make sure we are staying within 
              !                 half-way of the next land mass to the east
              !----------------------------------------------------------------
              ii        = ie_tmp
              icount    = 1
              found_lnd = .false.
              do while((icount.lt.2*offshore_pts).and.(.not.found_lnd))
                icount = icount + 1
                ii     = ii     + 1
                if(ii.gt.plon) ii = ii - plon
                if(landfrac(ii,j).gt.land_threshold) then
                  found_lnd = .true.
                endif
              end do
              ie_tmp = ie_tmp + icount/2 - 1

              ! Western shore:  look west to make sure we are staying within 
              !                 half-way of the next land mass to the west
              !----------------------------------------------------------------
              ii        = iw_tmp
              icount    = 1
              found_lnd = .false.
              do while((icount.lt.2*offshore_pts).and.(.not.found_lnd))
                icount = icount + 1
                ii     = ii     - 1
                if(ii.lt.1) ii = ii + plon
                if(landfrac(ii,j).gt.land_threshold) then
                  found_lnd = .true.
                endif
              end do
              iw_tmp = iw_tmp - icount/2 + 1
            endif

            if((i.le.iw_tmp).or.(i.ge.ie_tmp)) then
              write(iulog,*) 'Error:  Should never reach this branch.  Current land ' 
              write(iulog,*) 'point should be between the east and west ocn points.'
              write(iulog,*) 'iw, i, ie = ', iw_tmp, i, ie_tmp
              call endrun
            endif

            ! get E/W ice fraction from offshore points
            !---------------------------------------------
            ii = iw_tmp
            if(ii.lt.1) ii = ii+plon
            tw = ice(ii,j)
            ii = ie_tmp
            if(ii.gt.plon) ii = ii-plon
            te = ice(ii,j)

            wgt_w = float(ie_tmp -      i)/float(ie_tmp-iw_tmp)
            wgt_e = float(i      - iw_tmp)/float(ie_tmp-iw_tmp)
            if(wgt_w.gt.wgt_e) then
              ice(i,j) = tw
            else
              ice(i,j) = te
            endif
!             ice(i,j) = tw*wgt_w + te*wgt_e
          elseif((.not.found_ocn_w).and.(found_ocn_e)) then
            write(iulog,*) 'Error:  should never reach this branch. If even just one ' 
            write(iulog,*) 'ocean point exists on a latitude circle, it should be found in' 
            write(iulog,*) 'both directions'
            call endrun
          endif
        endif
      end do ! i = 1,plon
    end do ! j = 1, plat

    ! Sanity checks
    !---------------
    do j = 1,plat
    do i = 1,plon
      if((ice(i,j).gt.1.).or.(ice(i,j).lt.0.)) then
        write(iulog,*) 'ice fraction out of bounds'
        write(iulog,*) 'ice fraction  = ', ice (i,j)
        write(iulog,*) 'at i,j        = ', i,j
        call endrun
      end if
    end do
    end do

    ! End Routine
    !------------
    return
  end subroutine mask_ice
  !=====================================================================


  !=====================================================================
  subroutine mask_sst(plat, plon, sst, landfrac, icefrac)
    !
    ! Use land and ice masks to overwrite land fraction with SSTs to 
    ! prepare for interpolation of SSTs to new resolution
    !===================================================================
    implicit none
    !
    ! Passed variables
    !---------------------
    integer :: plat                 ! latitude  dimension of input/output fields
    integer :: plon                 ! longitude dimension of input/output fields
    real(r8):: sst     (plon,plat ) ! input: sst (K); output: masked sst in degrees C
    real(r8):: landfrac(plon,plat ) ! land mask
    real(r8):: icefrac (plon,plat ) ! ice  mask
    !
    ! Local Values
    !--------------
    integer :: i,ii,j
    integer :: icount
    real(r8):: tfreez             ! freezing point of fresh water(K)
                                  ! (C/K conversion factor)
    real(r8):: tfreezsw           ! freezing point of salt water (K)
    real(r8):: tw                 ! on lat circle, sst point west of land
    integer :: iw                 ! index of western sst point
    integer :: iw_tmp             ! index of western sst point
    real(r8):: te                 ! on lat circle, sst point east of land
    integer :: ie                 ! index of eastern sst point
    integer :: ie_tmp             ! index of eastern sst point
    integer :: offshore_pts       ! how many points offshore to sample SST's
    real(r8):: land_threshold     ! minimum fraction to be considered land
    real(r8):: wgt_w, wgt_e       ! interpolation weights
    logical :: found_ocn_w        ! true when ocn point is found on lat 
                                  ! circle west of current land point
    logical :: found_ocn_e        ! true when ocn point is found on lat 
                                  ! circle east of current land point
    logical :: found_lnd          ! true when lnd point is found

    tfreez         = 273.15_r8
    tfreezsw       = tfreez - 1.8_r8
    land_threshold = 0.0_r8
    offshore_pts   = 1

    if((offshore_pts.lt.1).or.(offshore_pts.gt.10)) then
      write(iulog,*) 'number of points offshore (in east/west direction)'
      write(iulog,*) 'to sample SSTs must be at least 1 but not more'
      write(iulog,*) 'than 10'
      write(iulog,*) 'offshore_pts = ', offshore_pts
      call endrun
    endif

    ! For each grid box that has seaice, set sst to "tfreezsw"
    !-----------------------------------------------------------
    do j = 1,plat
    do i = 1,plon
      if(icefrac(i,j).gt.0.) then
        sst(i,j) = tfreezsw
      endif
    end do
    end do

    ! Now, set all T's to a minimum of "tfreezsw"
    ! (no, this is not being redundant)
    !------------------------------------------------
    do j = 1,plat
    do i = 1,plon
      if(sst(i,j).lt.tfreezsw) then
        sst(i,j) = tfreezsw
      endif
    end do
    end do

    ! Overwrite land Ts with interpolated SST's
    !---------------------------------------------
    do j = 1,plat
      found_ocn_w = .false.
      found_ocn_e = .false.
      iw          = -999
      ie          = -999
      do i = 1,plon
        if(landfrac(i,j).le.land_threshold) then
          ! If found ocean, set tw/iw and move on to next point eastward
          !---------------------------------------------------------------
          found_ocn_w = .true.
          iw          = i
        else
          ! else, current point is land and will be overwritten
          !
          ! If "i=1", search westward over entire latitude circle for first
          ! ocean point west of current land point.
          !-------------------------------------------------------------------
          if(i.eq.1) then
            found_ocn_w = .false.
            ii = plon + 1
            do while((ii.gt.1).and.(.not.found_ocn_w))
              ii = ii - 1
              if(landfrac(ii,j).le.land_threshold) then
                iw          = ii-plon
                found_ocn_w = .true.
              endif
            end do
          endif

          ! Now search eastward over entire latitude circle for first ocean point
          ! encountered east of current land point
          !-----------------------------------------------------------------------
          if((ie.lt.i).and.(found_ocn_w)) then
            found_ocn_e = .false.
            ii = i
            ie = i
            do while((ie.le.plon+i).and.(.not.found_ocn_e))
              ie = ie + 1
              ii = ii + 1
              if(ii .gt. plon) ii = 1
              if(landfrac(ii,j) .le. land_threshold) then
                found_ocn_e = .true.
              endif
            end do
          endif

          if((.not.found_ocn_w).and.(.not.found_ocn_e)) then
            ! If no ocean points were found in this latitude circle, assume we
            ! are at Antarctica and set all SST's to seaice values
            !------------------------------------------------------------------
            sst(i,j) = tfreezsw
          elseif((found_ocn_w).and.(found_ocn_e)) then
            ! Else, interpolate SST's over land point
            !-----------------------------------------
            iw_tmp = iw
            ie_tmp = ie

            ! First, decide how many points offshore (in both directions) 
            ! to sample SST's for interpolation over land
            !--------------------------------------------------------------
            if(offshore_pts.gt.1) then
              ! Eastern shore:  look east to make sure we are staying within 
              !                 half-way of the next land mass to the east
              !---------------------------------------------------------------
              ii        = ie_tmp
              icount    = 1
              found_lnd = .false.
              do while((icount.lt.2*offshore_pts).and.(.not.found_lnd))
                icount = icount + 1
                ii     = ii     + 1
                if(ii .gt. plon) ii = ii - plon
                if(landfrac(ii,j).gt.land_threshold) then
                  found_lnd = .true.
                endif
              end do
              ie_tmp = ie_tmp + icount/2 - 1

              ! Western shore:  look west to make sure we are staying within 
              !                 half-way of the next land mass to the west
              !--------------------------------------------------------------
              ii        = iw_tmp
              icount    = 1
              found_lnd = .false.
              do while((icount.lt.2*offshore_pts).and.(.not.found_lnd))
                icount = icount + 1
                ii     = ii     - 1
                if(ii.lt.1) ii = ii + plon
                if(landfrac(ii,j) .gt. land_threshold) then
                  found_lnd = .true.
                endif
              end do
              iw_tmp = iw_tmp - icount/2 + 1
            endif

            if(i .le. iw_tmp .or. i .ge. ie_tmp) then
              write(iulog,*) 'Error:  Should never reach this'
              write(iulog,*) 'branch.  Current land point should'
              write(iulog,*) 'be between the east and west SST'
              write(iulog,*) 'points.'
              write(iulog,*) 'iw, i, ie = ', iw_tmp, i, ie_tmp
              call endrun
            endif

            ! get E/W SSTs from offshore points
            !------------------------------------
            ii = iw_tmp
            if(ii.lt.1) ii = ii + plon
            tw = sst(ii,j)
            ii = ie_tmp
            if(ii.gt.plon) ii = ii - plon
            te = sst(ii,j)

            wgt_w = float(ie_tmp -      i)/float(ie_tmp-iw_tmp)
            wgt_e = float(i      - iw_tmp)/float(ie_tmp-iw_tmp)
            if(wgt_w.gt.wgt_e) then
              sst(i,j) = tw
            else
              sst(i,j) = te
            endif
!            sst(i,j) = tw*wgt_w + te*wgt_e
          elseif((.not.found_ocn_w).and.(found_ocn_e)) then
            write(iulog,*) 'Error:  should never reach this branch.'
            write(iulog,*) 'If even just one ocean point exists on '
            write(iulog,*) 'a latitude circle, it should be found in'
            write(iulog,*) 'both directions'
            call endrun
          endif
        endif
      end do ! i = 1,plon
    end do ! j = 1, plat

    ! Sanity checks
    !----------------
    do j = 1, plat
    do i = 1,plon
      if((sst(i,j).gt.340.).or.(sst(i,j).lt.160.)) then
        write(iulog,*) 'sst out of bounds'
        write(iulog,*) 'sst           = ', sst     (i,j)
        write(iulog,*) 'at i,j        = ', i,j
        call endrun
      endif
    end do
    end do

    ! Output degrees C
    !-------------------
    do j = 1, plat
    do i = 1,plon
      sst(i,j) = sst(i,j) - tfreez
    end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine mask_sst
  !=====================================================================


  !=====================================================================
  subroutine jra_25_press_full_levels(plevp1, plev, plat, plon, psi, psm)
    !
    ! Interpolate JRA-25 full-level pressures from interface pressures
    !
    !===================================================================
    implicit none
    !
    ! Passed variables
    !-------------------
    integer:: plevp1                ! vertical  dim of interface  pressure field (input )
    integer:: plev                  ! vertical  dim of full-level pressure field (output)
    integer:: plat                  ! latitude  dim of input/output fields
    integer:: plon                  ! longitude dim of input/output fields
    real*8 :: psi(plon,plat,plevp1) ! input  interface  pressure field
    real*8 :: psm(plon,plat,plev  ) ! output full-level pressure field
    !
    ! Local values
    !-----------------
    integer:: i,j,k
    real*8 :: one, two

    one = 1._r8
    two = 2._r8

    do k = 2,plev
    do j = 1,plat
    do i = 1,plon
      psm(i,j,k) = (psi(i,j,k+1)*log(psi(i,j,k+1))) - (psi(i,j,k)*log(psi(i,j,k)))
      psm(i,j,k) = psm(i,j,k)/(psi(i,j,k+1)-psi(i,j,k))
      psm(i,j,k) = psm(i,j,k) - one
      psm(i,j,k) = exp(psm(i,j,k))
    end do
    end do
    end do

    do j = 1,plat
    do i = 1,plon
      psm(i,j,1) = psi(i,j,2)/two
    end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine jra_25_press_full_levels
  !=====================================================================


  !=====================================================================
  subroutine binning(plev, plato, plono, plat, plon, xx, yy, clat, clon, gw, &
                     clato, clono, gwo, bin_factor, dyn_flag, dyn_flago      )
    !
    ! Grid-Box Binning
    !===================================================================
    implicit none
    !
    ! Passed variables
    !------------------
    integer :: plev                    ! vertical dimension of input/output field
    integer :: plato                   ! latitude dimension of output field
    integer :: plono                   ! longitude dimension of output field
    integer :: plat                    ! latitude dimension of input field
    integer :: plon                    ! longitude dimension of input field
    real(r8):: xx   (plon ,plat ,plev) ! input analysis field
    real(r8):: yy   (plono,plato,plev) ! horizontally interpolated (output) field
    real(r8):: clat (plat )            ! Input lat in degrees starting from southern-most lat
    real(r8):: clon (plon )            ! Input lon in degrees from 0 deg and moving eastward
    real(r8):: gw   (plat )            ! Input Gaussian wgts (if relevant grid)
    real(r8):: clato(plato)            ! Output lat in degrees starting from southern-most lat
    real(r8):: clono(plono)            ! Output lon in degrees from 0 deg and moving eastward
    real(r8):: gwo  (plato)            ! Output Gaussian wgts (if relevant grid)
    integer :: dyn_flag                ! Dynamics flag of input grid:   Eul=1, FV=0
    integer :: dyn_flago               ! Dynamics flag of output grid:  Eul=1, FV=0
    real(r8):: bin_factor              ! bin-box area expansion/contraction factor 
                                    ! relative to output grid-box area.
    !
    ! Local values
    !----------------
    integer,parameter:: max_segs = 10000      ! Max # of box segments
    integer :: i, j, ii, jj, k, platp1, platop1  ! Indices
    integer :: nx, ny, nx_max, ny_max            ! Indices
    integer :: plon2, plonhalf
    integer :: i_in(max_segs),i_out(max_segs)
    integer :: j_in(max_segs),j_out(max_segs)
    real(r8):: sum
    real(r8):: xx_loc(plon*2, plat, plev)
    real(r8):: pi, pio180, pio2, one, factor
    real(r8):: flat  (plat   ),flon  (plon*2  )
    real(r8):: flato (plato  ),flono (plono   )
    real(r8):: flati (plat +1), floni(plon*2+1)
    real(r8):: flatoi(plato+1),flonoi(plono+1 )
    real(r8):: tmps, tmpn, tmp(plono,plato)
    real(r8):: edge_w (plon*2), edge_e (plon*2), edge_s (plat ), edge_n (plat )
    real(r8):: edgeo_w(plono ), edgeo_e(plono ), edgeo_s(plato), edgeo_n(plato)
    real(r8):: sin_s (plat ),sin_n (plat )
    real(r8):: sino_s(plato),sino_n(plato)
    real(r8):: dx(max_segs), dy(max_segs)
    real(r8):: distmin, dist, two, zero

    zero     = 0._r8
    platp1   = plat  + 1
    platop1  = plato + 1
    plon2    = plon*2
    plonhalf = plon/2
    one      = 1._r8
    pi       = 4._r8*atan(one)
    pio180   = pi/180._r8
    pio2     = pi/2._r8
    two      = 2._r8

    ! Sanity checks
    !------------------
    if(bin_factor.lt.0.05_r8) then
      write(iulog,*) 'ERROR ("BINNING"):  binning factor out of range'
      write(iulog,*) 'bin_factor = ', bin_factor
      call endrun
    endif

    ! Copy input data to wrap-around array 
    ! (wrap half-way around globe at each end of x-direction)
    !-------------------------------------------------------------
    do k = 1,plev
    do j = 1,plat
      ii = plonhalf
      do i = 1,plon2
        ii = ii + 1
        if(ii.gt.plon) ii = 1
        xx_loc(i,j,k) = xx(ii,j,k)
      end do
    end do
    end do

    ! Convert to radians 
    ! (wrap half-way around globe at each end of x-direction of input grid)
    !----------------------------------------------------------------------
    ii = plonhalf
    do i = 1,plon2
      ii = ii + 1
      if(ii.gt.plon) ii = 1
      if(i .le. plonhalf      ) flon(i) = clon(ii)*pio180-4*pio2
      if(i .gt.(plonhalf+plon)) flon(i) = clon(ii)*pio180+4*pio2
      if((i.gt.plonhalf).and.(i.le.plonhalf+plon)) then
        flon(i) = clon(ii)*pio180
      endif
    end do
    do j = 1,plat
      flat (j) = clat (j)*pio180
    end do
    do i = 1,plono
      flono(i) = clono(i)*pio180
    end do
    do j = 1,plato
      flato(j) = clato(j)*pio180
    end do

    ! Compute box edges
    !--------------------
    floni(       1) = 0.5_r8*(flon(1) + flon(plon*2) - pio2*8)
    floni(plon*2+1) = 0.5_r8*(flon(1) + flon(plon*2) + pio2*8)
    do i = 2,plon*2
      floni(i) = 0.5_r8*(flon(i-1) + flon(i))
    end do

    sum = 0._r8
    flati(1     ) = -pio2
    flati(platp1) =  pio2
    do j = 1,plat
      if(dyn_flag.eq.1) then
        sum = sum + gw (j)
        if((sum.gt.2._r8).and.(j.lt.plat)) then
          write(iulog,*) 'ERROR ("BINNING"): Something wrong with'
          write(iulog,*) 'input Gaussian weights'
          call endrun
        elseif((sum.gt.2._r8).and.(j.eq.plat)) then
          sum = min (sum, two)
        endif
        flati(j+1) = asin( sum-one )
      elseif(dyn_flag.eq.0) then
        if(j.lt.plat) flati(j+1) = 0.5_r8*(flat(j) + flat(j+1))
      else
        write(iulog,*) 'ERROR:  "BINNING" should not reach this branch'
        call endrun
      endif
    end do
      
    flonoi(      1) = 0.5_r8*(flono(1) + flono(plono) - pio2*4)
    flonoi(plono+1) = 0.5_r8*(flono(1) + flono(plono) + pio2*4)
    do i = 2,plono
      flonoi(i) = 0.5_r8*(flono(i-1) + flono(i))
    end do

    sum = 0._r8
    flatoi(1      ) = -pio2
    flatoi(platop1) =  pio2
    sum = 0._r8
    do j = 1,plato
      if(dyn_flago.eq.1) then
        sum = sum + gwo (j)
        if((sum.gt.2._r8).and.(j.lt. plato)) then
          write(iulog,*) 'ERROR ("BINNING"): Something wrong with'
          write(iulog,*) 'output Gaussian weights'
          call endrun
        elseif((sum.gt.2._r8).and.(j.eq.plato)) then
          sum = min (sum, two)
        endif
        flatoi(j+1) = asin( sum-one )
      elseif(dyn_flago.eq.0) then
        if(j.lt.plato) flatoi(j+1) = 0.5_r8*(flato(j) + flato(j+1))
      else
        write(iulog,*) 'ERROR:  "BINNING" should not reach this branch'
        call endrun
      endif
    end do

    ! Copy grid interfaces to "edge" arrays
    !---------------------------------------
    do i = 1,plon*2
      edge_w(i) = floni(i  )
      edge_e(i) = floni(i+1)
    end do

    do j = 1,plat
      edge_s(j) = flati(j  )
      edge_n(j) = flati(j+1)
      sin_s (j) = sin(edge_s(j))
      sin_n (j) = sin(edge_n(j))
    end do

    ! Expand/contract bin box area for each output grid box by "bin_factor"
    !----------------------------------------------------------------------
    factor = sqrt(bin_factor)

    do i = 1,plono
      edgeo_w(i) = flono(i) - ( flono (i  ) - flonoi(i) )*factor
      edgeo_e(i) = flono(i) + ( flonoi(i+1) - flono (i) )*factor
    end do

    do j = 1,plato
      tmps       = flato(j) - ( flato (j  ) - flatoi(j) )*factor
      tmpn       = flato(j) + ( flatoi(j+1) - flato (j) )*factor
      edgeo_s(j) = max( tmps, -pio2) - max( ( tmpn - pio2), zero)
      edgeo_n(j) = min( tmpn,  pio2) + max( (-pio2 - tmps), zero)
      sino_s (j) = sin(edgeo_s(j))
      sino_n (j) = sin(edgeo_n(j))
    end do

    ! Make vector of box segments in x-direction
    !-------------------------------------------
    nx = 0
    do i = 1,plono
    do ii = 1,plon*2
      if((edge_e(ii).gt.edgeo_w(i)).and.(edgeo_e(i).gt.edge_w(ii))) then
        nx = nx + 1
        if(nx.gt.max_segs) then
          write(iulog,*) 'ERROR  ("BINNING"):  number of box'
          write(iulog,*) 'segments greater than "max_segs"'
          call endrun
        endif
        i_in (nx) = ii
        i_out(nx) = i
        dx   (nx) = min(min(min(edge_e(ii)- edge_w(ii),   &
                               edgeo_e(i )-edgeo_w(i ) ), &
                                edge_e(ii)-edgeo_w(i ) ), &
                               edgeo_e(i )- edge_w(ii) )
      endif
      if(edge_w (ii).ge.edgeo_e( i)) exit
    end do
    end do

    ! Make vector of box segments in y-direction
    !----------------------------------------------
    ny = 0
    do j = 1,plato
    do jj = 1,plat
      if((edge_n(jj).gt.edgeo_s(j)).and.(edgeo_n(j).gt.edge_s(jj))) then
        ny = ny + 1
        if(ny.gt.max_segs) then
          write(iulog,*) 'ERROR  ("BINNING"):  number of box'
          write(iulog,*) 'segments greater than "max_segs"'
          call endrun
        endif
        j_in (ny) = jj
        j_out(ny) = j
        distmin   = edge_n(jj)-edge_s(jj)
        dy(ny)    = sin_n (jj)-sin_s (jj)
        dist      = edgeo_n(j)-edgeo_s(j)
        if(dist.lt.distmin) then
          distmin = dist
          dy(ny)  = sino_n(j)-sino_s(j)
        endif
        dist = edge_n(jj)-edgeo_s(j)
        if(dist.lt.distmin) then
          distmin = dist
          dy(ny)  = sin_n(jj)-sino_s(j)
        endif
        dist = edgeo_n(j)-edge_s(jj)
        if(dist.lt.distmin) then
          distmin = dist
          dy(ny)  = sino_n(j)-sin_s(jj)
        endif
      endif
      if(edge_s (jj) .ge. edgeo_n( j)) exit
    end do
    end do

    nx_max = nx
    ny_max = ny

    ! Begin weighted binning
    !-----------------------
    do k = 1,plev
    do j = 1,plato
    do i = 1,plono
      yy(i,j,k) = 0._r8
    end do
    end do
    end do

    do k = 1,plev
    do ny = 1,ny_max
      j  = j_out(ny)
      jj = j_in (ny)
      do nx = 1,nx_max
        i  = i_out(nx)
        ii = i_in (nx)
        yy(i,j,k) = yy(i,j,k) + xx_loc(ii,jj,k)*dx(nx)*dy(ny)
      end do
    end do
    end do

    ! Normalize
    !--------------
    do j = 1,plato
    do i = 1,plono
      tmp(i,j) = (edgeo_e(i) - edgeo_w(i))*(sino_n(j) - sino_s(j))
    end do
    end do

    do k = 1,plev
    do j = 1,plato
    do i = 1,plono
      yy(i,j,k) = yy(i,j,k)/tmp(i,j)
    end do
    end do
    end do

    ! End Routine
    !--------------
    return
  end subroutine binning
  !=====================================================================


  !=====================================================================
  subroutine cubic_opt1(plev, plato, plono, plat, plon, platm2, xx, yy,  &
                        pext, xx_exts, clat, clon, clato, clono, limdr)
    !
    ! Horizontal Cubic Interpolation Driver
    !===================================================================
    implicit none
    !
    ! Passed Variables
    !----------------------
    integer :: plev                      ! vertical dimension of input/output field
    integer :: plato                     ! latitude dimension of output field
    integer :: plono                     ! longitude dimension of output field
    integer :: plat                      ! latitude dimension of input field
    integer :: plon                      ! longitude dimension of input field
    integer :: platm2                    ! input latitude dimension(less possible pole points)
    integer :: pext                      ! # of latitude extensions
    real(r8):: xx     (plon ,plat ,plev) ! input analysis field
    real(r8):: yy     (plono,plato,plev) ! horizontally interpolated (output) field
    real(r8):: xx_exts(plon ,pext ,plev) ! input latitude extensions
    real(r8):: clat   (plat )            ! Input lat in degrees starting from southern-most lat
    real(r8):: clon   (plon )            ! Input lon in degrees from 0 deg and moving eastward
    real(r8):: clato  (plato)            ! Output lat in degrees starting from southern-most lat
    real(r8):: clono  (plono)            ! Output lon in degrees from 0 deg and moving eastward
    logical :: limdr                     ! Flag to use SCM0 derivative estimate limiter
    !
    ! Local Values
    !---------------
    integer,parameter:: nxpt   = 1 ! Extra point in local 4-pt interpolation stencil
    integer,parameter:: jintmx = 2 ! number of extra latitudes at southern and northern
                                   ! extremes of input latitude grid (for interp purposes)

    integer:: i, j, k      ! Indices
    integer:: platd        ! Total latitude dimension of input grid including extensions
    integer:: plond        ! Total longitude dimension of input grid including extensions
    logical:: npole        ! flag set if northern-most input latitude is a pole point
    logical:: spole        ! flag set if southern-most input latitude is a pole point
    logical:: lpole        ! flag set if both input latitudes are pole points

    real(r8):: xxm2(plon ,platm2 ,plev) ! input analysis field (less possible pole points)

    npole  = .false.
    spole  = .false.
    lpole  = .false.
    if(clat(   1) .lt. -89.9999_r8) spole = .true.
    if(clat(plat) .gt.  89.9999_r8) npole = .true.
    if(platm2 .ne. plat-2 ) then
      write(iulog,*) 'Error in CUBIC_OPT1:  '
      write(iulog,*) '"platm2" must be same as "plat-2" '
      call endrun
    endif
    if((spole.and.(.not.npole)).or.((.not.spole).and.npole)) then
      write(iulog,*) 'Error in CUBIC_OPT1:  '
      write(iulog,*) 'input data must either include BOTH pole points or none'
      call endrun
    endif

    platd = plat + 2*nxpt + 2*jintmx
    plond = plon + 1 + 2*nxpt

    ! If input data has pole points, strip them out for now
    !-------------------------------------------------------
    if(spole .and. npole) then
      platd = platm2 + 2*nxpt + 2*jintmx
      lpole = .true.
      do k = 1,plev
      do j = 1,platm2
      do i = 1,plon
        xxm2(i,j,k) = xx(i,j+1,k)
      end do
      end do
      end do
    endif

    ! Call cubic interpolator
    !-------------------------
    if(lpole) then
      call cubic_slav(nxpt, jintmx, platd, plond  ,          &
                      plev, plato , plono, platm2 , plon   , &
                      xxm2, yy    , pext , xx_exts, clat(2), &
                      clon, clato , clono, limdr             )
    else
      call cubic_slav(nxpt, jintmx, platd, plond  ,       &
                      plev, plato , plono, plat   , plon, &
                      xx  , yy    , pext , xx_exts, clat, &
                      clon, clato , clono, limdr          )
    endif

    ! End Routine
    !---------------
    return
  end subroutine cubic_opt1
  !=====================================================================


  !=====================================================================
  subroutine cubic_slav(nxpt, jintmx, platd, plond   ,      &
                        plev, plato , plono, plat    ,plon, &
                        xx  , yy    , pext , xx_exts ,clat, &
                        clon, clato , clono, limdr          )
    !
    ! Horizontal Cubic Interpolation
    !===================================================================
    implicit none
    !
    ! Passed Variables 
    !---------------------
    integer:: nxpt       ! Extra point in local 4-pt interpolation stencil
    integer:: jintmx     ! number of extra latitudes at southern and northern extremes of  
                         ! input latitude grid (created for interpolation purposes)
    integer:: platd      ! Total latitude dimension of input grid including extensions
    integer:: plond      ! Total longitude dimension of input grid including extensions
    integer:: plev       ! vertical dimension of input/output field
    integer:: plato      ! latitude dimension of output field
    integer:: plono      ! longitude dimension of output field
    integer:: plat       ! latitude dimension of input field
    integer:: plon       ! longitude dimension of input field
    integer:: pext       ! # of latitude extensions

    real(r8):: xx     (plon ,plat ,plev) ! input analysis field
    real(r8):: yy     (plono,plato,plev) ! horizontally interpolated (output) field
    real(r8):: xx_exts(plon ,pext ,plev) ! input latitude extensions
    real(r8):: clat   (plat )            ! Input lat in degrees starting from southern-most lat
    real(r8):: clon   (plon )            ! Input lon in degrees from 0 deg and moving eastward
    real(r8):: clato  (plato)            ! Output lat in degrees starting from southern-most lat
    real(r8):: clono  (plono)            ! Input lon in degrees from 0 deg and moving eastward
    logical :: limdr                     ! Flag to use SCM0 derivative estimate limiter
    !
    ! Local values
    !-----------------
    real(r8):: xxx   (plond,platd,plev) ! input field on extended grid
    real(r8):: dphi  (platd)            ! latitude intervals (radians)
    real(r8):: rdphi                    ! reciprocal of dphi for interpolant y-interval
    real(r8):: rdx                      ! recip of del-x
    real(r8):: phi   (platd)            ! latitudes on extended grid (radians)
    real(r8):: lam   (plond)            ! longitudes on extended grid radians
    real(r8):: dlam  (plond)            ! delta-lam
    real(r8):: phidp (plato)            ! latitudes of output field (radians)
    real(r8):: lamdp (plono)            ! longitudes of output field (radians)
    real(r8):: lbasdx(4,2,plono)        ! longitudinal basis functions for interpolation
    real(r8):: lbasdy(4,2,plato)        ! latitudinal basis functions for interpolation
    real(r8):: pi                       ! 3.1415926...
    real(r8):: pio180                   ! pi/180.
    real(r8):: pi2                      ! pi*2
    real(r8):: fintx(4)
    real(r8):: hs (plato)
    real(r8):: hn (plato)
    real(r8):: dhs(plato)
    real(r8):: dhn(plato)
    real(r8):: hl (plono)
    real(r8):: hr (plono)
    real(r8):: dhl(plono)
    real(r8):: dhr(plono)
    integer :: jdp(plato)
    integer :: idp(plono)
    integer :: istart                  ! index for first analysis long.
    integer :: istop                   ! index for last  analysis long.
    integer :: jstart                  ! index for first analysis lat.
    integer :: jstop                   ! index for last  analysis lat.
    integer :: i, j, k, n              ! Indices
    integer :: ii, jj

    ! Define constants
    !-------------------
    pi     = 4._r8*atan(1._r8)
    pio180 = pi/180._r8
    pi2    = pi*2._r8
    istart = nxpt   + 1
    jstart = nxpt   + 1 + jintmx
    istop  = istart - 1 + plon
    jstop  = jstart - 1 + plat

    ! Make sure latitude dimension of "xx_ext" is same as
    ! total number of extended latitudes
    !-----------------------------------------------------
    if((jstart-1)*2.ne.pext) then
      write(iulog,*) 'Error in CUBIC_SLAV:  '
      write(iulog,*) '   dimension of "xx_ext not same as total number'
      write(iulog,*) '   of extened latitudes'
      call endrun
    endif

    ! Define lats/lons (in radians) for input and output grids
    !----------------------------------------------------------
    do j = 1,plat
      phi(jstart-1+j) = clat(j)*pio180
    end do
    do i = 1,plon
      lam(istart-1+i) = clon(i)*pio180
    end do
    do j = 1,plato
      phidp(j) = clato(j)*pio180
    end do
    do i = 1,plono
      lamdp(i) = clono(i)*pio180
      if(lamdp(i).lt.0._r8) lamdp(i) = lamdp(i) + pi2
    end do

    ! North and south poles.
    !---------------------------
    phi(jstart-1) = -pi/2.0_r8
    phi(jstop +1) =  pi/2.0_r8

    ! Extend Gauss latitudes below south pole so that the spacing above
    ! the pole is symmetric, and phi is decreasing, i.e., phi < -pi/2
    !-------------------------------------------------------------------
    if(jstart.gt.2) then
      do j = 1,jstart-2
        phi(j) = -pi - phi(2*jstart-2-j)
      end do
    endif

    ! Analogously for Northern Hemisphere
    !--------------------------------------
    if(platd.gt.jstop+1)then
      do j = jstop+2,platd
        phi(j) = pi - phi(2*jstop+2-j)
      end do
    endif

    ! Fill East/West extensions of lam array
    !----------------------------------------
    do i = 1,istart-1
      lam(i) = lam(i+istop-istart+1) - pi2
    end do
    do i = istop+1,plond
      lam(i) = lam(i-istop+istart-1) + pi2
    end do

    ! Compute delta values in X/Y
    !-------------------------------
    do i = 1,plond-1
      dlam(i) = lam(i+1) - lam(i)
    end do
    do j = 1,platd-1
      dphi(j) = phi(j+1) - phi(j)
    end do

    ! Copy input to extended grid and fill extensions
    !-------------------------------------------------
    do k = 1,plev
      do j = 1,plat
      do i = 1,plon
        xxx(istart-1+i,jstart-1+j,k) = xx(i,j,k)
      end do
      end do

      ! fill N/S poles
      !----------------
      do i = 1,plon
        xxx(istart-1+i,jstart-1,k) = xx_exts(i,pext/2  ,k)
        xxx(istart-1+i,jstop +1,k) = xx_exts(i,pext/2+1,k)
      end do

      ! fill beyond N/S pole latitudes
      !-------------------------------
      if(jstart-2.ge.1)then
        do j=1,jstart-2
        do i = 1,plon
          xxx(istart-1+i,j,k) = xx_exts(i,j,k)
        end do
        end do
      endif

      if(jstop+2.le.platd) then
        jj = 1
        do j=jstop+2,platd
          jj = jj + 1
          do i = 1,plon
            xxx(istart-1+i,j,k) = xx_exts(i,pext/2+jj,k)
          end do
        end do
      endif
    end do ! k = 1,plev

    ! Fill E/W points
    !---------------------
    call extx(nxpt, platd , plond, plev, plon, istart, istop, jstart, jstop, xxx)

    ! Interpolation weights for x-direction
    !--------------------------------------
    do i = 1,plono
      call lcdbas(plond, lam, dlam, lamdp(i), idp(i),                        &
                  lbasdx(1,1,i), lbasdx(1,2,i), hl(i), hr(i) , dhl(i), dhr(i))
    end do

    ! Interpolation weights for y-direction
    !-----------------------------------------
    do j = 1,plato
      call lcdbas(platd, phi, dphi, phidp(j), jdp(j),                        &
                  lbasdy(1,1,j), lbasdy(1,2,j), hs(j), hn(j) , dhs(j), dhn(j))
    end do

    ! Hermite Cubic interpolation 
    ! (using 4X4 input stencil for each output grid point)
    !--------------------------------------------------------
    do k = 1,plev
    do jj = 1,plato
      rdphi = 1./dphi(jdp(jj))
      do ii = 1,plono
        rdx = 1./dlam(idp(ii))
        ! x-interpolation (over 4 adjacent latitudes in stencil)
        !--------------------------------------------------------
        do n = 1,4
          j = jdp(jj)-2+n
          call lcdint(xxx(idp(ii)-1,j,k),lbasdx(1,1,ii), rdx, limdr, &
                      hl(ii), hr(ii), dhl(ii), dhr(ii), fintx(n)     )
        end do

        ! y-interpolation
        !----------------
        call lcdint(fintx, lbasdy(1,1,jj), rdphi, limdr,         &
                    hs(jj), hn(jj), dhs(jj), dhn(jj), yy(ii,jj,k))
      end do
    end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine cubic_slav
  !=====================================================================


  !=====================================================================
  subroutine extx(nxpt, platd , plond, plev  , plon, istart, istop, jstart, jstop, fb)
    ! 
    ! Copy data to the longitude extensions of the extended array
    !===================================================================
    implicit none
    !
    ! Passed variables
    !---------------------
    integer :: nxpt                 ! Extra point in local 4-pt interpolation stencil
    integer :: platd                ! Total lat dimension of input grid including extensions
    integer :: plond                ! Total lon dimension of input grid including extensions
    integer :: plev                 ! vertical dimension of input/output field
    integer :: plon                 ! longitude dimension of input field
    integer :: istart               ! index for first analysis long.
    integer :: istop                ! index for last  analysis long.
    integer :: jstart               ! index for first analysis lat.
    integer :: jstop                ! index for last  analysis lat.
    real(r8):: fb(plond,platd,plev) ! input field on extended grid
    !
    ! Local values
    !----------------
    integer:: i              ! longitude index
    integer:: j              ! latitude  index
    integer:: k              ! vertical  index
    integer:: i2pi           ! start of eastern long. extension

    ! Fill west edge points.
    !------------------------
    if(nxpt .ge. 1) then
      do k=1,plev
      do j=1,platd
      do i=1,nxpt
        fb(i,j,k) = fb(i+plon,j,k)
      end do
      end do
      end do
    endif

    ! Fill east edge points
    !-----------------------
    i2pi = nxpt + plon + 1
    do k=1,plev
    do j=1,platd
    do i=i2pi,plond
      fb(i,j,k) = fb(i-plon,j,k)
    end do
    end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine extx
  !=====================================================================


  !=====================================================================
  subroutine extys(plon, plev, plat, extent_dim, clat, fb, fb_extents)

    ! 
    ! Fill latitude extensions of a scalar extended array
    ! 
    ! Method: 
    ! This is done in 2 steps:
    !   1) interpolate to the pole points; use the mean field value on the
    !      Gaussian latitude closest to the pole.
    !   2) add latitude lines beyond the poles.
    !===================================================================
    implicit none
    !
    ! Passed Variables
    !------------------
    integer :: plon                             ! longitude dimension of input field
    integer :: plev                             ! vertical dimension of input/output
    integer :: plat                             ! Total latitude dimension of input grid
    integer :: extent_dim                       ! # of latitude extensions
    real(r8)::  clat      (plat )                ! Input latitudes in degrees
    real(r8)::  fb        (plon,plat,plev)       ! input field on extended grid
    real(r8)::  fb_extents(plon,extent_dim,plev) ! latitude extensions
    !
    ! Local Values
    !----------------
    real(r8):: zave           ! accumulator for zonal averaging
    integer :: i,j,k          ! indices
    integer :: plon2          ! half the number of real longitudes
    integer :: platn, plats   ! indices
    logical :: npole          ! flag set if northern-most input latitude is a pole point
    logical :: spole          ! flag set if southern-most input latitude is a pole point

    ! This code Hardwired to 6 latitude extensions
    !----------------------------------------------
    if(extent_dim .ne. 6) then
      write(iulog,*) 'Error in EXTYS:  latitude extensions hardwired'
      write(iulog,*) 'to 6'
      call endrun
    endif

    ! Check if input dataset has pole points
    !-----------------------------------------
    npole  = .false.
    spole  = .false.
    if(clat(   1) .lt. -89.9999_r8) spole = .true.
    if(clat(plat) .gt.  89.9999_r8) npole = .true.

    plon2 = plon/2

    ! Fill north pole line.
    !----------------------
    if(npole) then
      do k=1,plev
      do i=1,plon
        fb_extents(i,4,k) = fb(i,plat,k)
      end do
      end do
    else
      do k=1,plev
        zave = 0.0_r8
        do i = 1,plon
          zave = zave + fb(i,plat,k)
        end do
        zave = zave/plon
        do i=1,plon
          fb_extents(i,4,k) = zave
        end do
      end do
    endif

    ! Fill northern lines beyond pole line.
    !--------------------------------------
    if(npole) then
      platn = plat-1
      plats = plat-2
    else
      platn = plat
      plats = plat-1
    endif

    do k=1,plev
    do i=1,plon2
      fb_extents(      i,5,k) = fb(plon2+i,platn,k)
      fb_extents(plon2+i,5,k) = fb(      i,platn,k)
      fb_extents(      i,6,k) = fb(plon2+i,plats,k)
      fb_extents(plon2+i,6,k) = fb(      i,plats,k)
    end do
    end do

    ! Fill south pole line.
    !------------------------
    if(spole) then
      do k=1,plev
      do i=1,plon
        fb_extents(i,3,k) = fb(i,1,k)
      end do
      end do
    else
      do k=1,plev
        zave = 0.0_r8
        do i = 1,plon
          zave = zave + fb(i,1,k)
        end do
        zave = zave/plon
        do i=1,plon
          fb_extents(i,3,k) = zave
        end do
      end do
    endif

    ! Fill southern lines beyond pole line.
    !--------------------------------------
    if(spole) then
      platn = 3
      plats = 2
    else
      platn = 2
      plats = 1
    endif

    do k=1,plev
    do i=1,plon2
      fb_extents(      i,1,k) = fb(plon2+i,platn,k)
      fb_extents(plon2+i,1,k) = fb(      i,platn,k)
      fb_extents(      i,2,k) = fb(plon2+i,plats,k)
      fb_extents(plon2+i,2,k) = fb(      i,plats,k)
    end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine extys
  !=====================================================================


  !=====================================================================
  subroutine extyv(plon, plev, plat, extent_dim,                 &
                   clon, clat, fbu, fbv, fbu_extents, fbv_extents)
    ! 
    ! Fill latitude extensions of a vector quantity
    ! 
    ! Method: 
    ! This is done in 2 steps:
    !   1) interpolate to the pole points; project the orthogonal wave 1
    !      of U and V to the pole; use the Gaussian latitude closest to the pole.
    !   2) add latitude lines beyond the poles.
    !===================================================================
    implicit none
    !
    ! Passed Variables
    !------------------
    integer :: plon                              ! longitude dimension of input field
    integer :: plev                              ! vertical dimension of input/output
    integer :: plat                              ! Total latitude dimension of input grid
    integer :: extent_dim                        ! # of latitude extensions
    real(r8):: clon       (plon )                ! Input longitude in degrees
    real(r8):: clat       (plat )                ! Input latitude  in degrees
    real(r8):: fbu        (plon,plat,plev)       ! input field on extended grid
    real(r8):: fbv        (plon,plat,plev)       ! input field on extended grid
    real(r8):: fbu_extents(plon,extent_dim,plev) ! latitude extensions
    real(r8):: fbv_extents(plon,extent_dim,plev) ! latitude extensions
    !
    ! Local Values
    !----------------
    real(r8):: zave              ! accumulator for zonal averaging
    integer :: i,j,k             ! indices
    integer :: plon2             ! half the number of real longitudes
    real(r8):: pi
    real(r8):: zavecv
    real(r8):: zavesv
    real(r8):: zavecu
    real(r8):: zavesu
    real(r8):: zavcus  
    real(r8):: zaucvs
    real(r8):: sinlam(plon)      ! Sin of lamda at all points in lat circle
    real(r8):: coslam(plon)      ! Cos of lamda at all points in lat circle
    logical :: npole             ! flag set if northern-most input latitude is a pole point
    logical :: spole             ! flag set if southern-most input latitude is a pole point

    ! This code Hardwired to 6 latitude extensions
    !---------------------------------------------
    if(extent_dim.ne.6) then
      write(iulog,*) 'Error in EXTYV:  latitude extensions hardwired to 6'
      call endrun
    endif

    ! Check if input dataset has pole points
    !---------------------------------------
    npole = .false.
    spole = .false.
    if(clat(   1).lt.-89.9999_r8) spole = .true.
    if(clat(plat).gt. 89.9999_r8) npole = .true.

    plon2 = plon/2
    pi = 4._r8*atan(1._r8)
    do i = 1,plon
      sinlam(i) = sin( clon(i)*pi/180._r8 )
      coslam(i) = cos( clon(i)*pi/180._r8 )
    end do

    ! Fill north pole line.
    !-----------------------
    if(npole) then
      do k = 1,plev
      do i = 1,plon
        fbv_extents(i,4,k) = fbv(i,plat,k)
        fbu_extents(i,4,k) = fbu(i,plat,k)
      end do
      end do
    else
      do k = 1,plev
        zavecv = 0.0_r8
        zavesv = 0.0_r8
        zavecu = 0.0_r8
        zavesu = 0.0_r8
        do i = 1,plon
          zavecv = zavecv + fbv(i,plat,k  )*coslam(i)
          zavesv = zavesv + fbv(i,plat,k  )*sinlam(i)
          zavecu = zavecu + fbu(i,plat,k  )*coslam(i)
          zavesu = zavesu + fbu(i,plat,k  )*sinlam(i)
        end do
        zavcus = (zavecv + zavesu)/plon
        zaucvs = (zavecu - zavesv)/plon
        do i = 1,plon
          fbv_extents(i,4,k) = zavcus*coslam(i) - zaucvs*sinlam(i)
          fbu_extents(i,4,k) = zaucvs*coslam(i) + zavcus*sinlam(i)
        end do
      end do
    endif

    ! Fill northern lines beyond pole line.
    !--------------------------------------
    do k=1,plev
    do i=1,plon2
      fbu_extents(      i,5,k) = fbu(plon2+i,plat  ,k)
      fbu_extents(plon2+i,5,k) = fbu(      i,plat  ,k)
      fbv_extents(      i,5,k) = fbv(plon2+i,plat  ,k)
      fbv_extents(plon2+i,5,k) = fbv(      i,plat  ,k)
      fbu_extents(      i,6,k) = fbu(plon2+i,plat-1,k)
      fbu_extents(plon2+i,6,k) = fbu(      i,plat-1,k)
      fbv_extents(      i,6,k) = fbv(plon2+i,plat-1,k)
      fbv_extents(plon2+i,6,k) = fbv(      i,plat-1,k)
    end do
    end do

    ! Fill south pole line.
    !----------------------
    if(spole) then
      do k = 1,plev
      do i = 1,plon
        fbv_extents(i,3,k) = fbv(i,1,k)
        fbu_extents(i,3,k) = fbu(i,1,k)
      end do
      end do
    else
      do k = 1,plev
        zavecv = 0.0_r8
        zavesv = 0.0_r8
        zavecu = 0.0_r8
        zavesu = 0.0_r8
        do i = 1,plon
          zavecv = zavecv + fbv(i,1,k  )*coslam(i)
          zavesv = zavesv + fbv(i,1,k  )*sinlam(i)
          zavecu = zavecu + fbu(i,1,k  )*coslam(i)
          zavesu = zavesu + fbu(i,1,k  )*sinlam(i)
        end do
        zavcus = (zavecv - zavesu)/plon
        zaucvs = (zavecu + zavesv)/plon
        do i = 1,plon
          fbv_extents(i,3,k) = zavcus*coslam(i) + zaucvs*sinlam(i)
          fbu_extents(i,3,k) = zaucvs*coslam(i) - zavcus*sinlam(i)
        end do
      end do
    endif

    ! Fill southern lines beyond pole line.
    !---------------------------------------
    do k=1,plev
    do i=1,plon2
      fbu_extents(      i,2,k) = fbu(plon2+i,1,k)
      fbu_extents(plon2+i,2,k) = fbu(      i,1,k)
      fbv_extents(      i,2,k) = fbv(plon2+i,1,k)
      fbv_extents(plon2+i,2,k) = fbv(      i,1,k)
      fbu_extents(      i,1,k) = fbu(plon2+i,2,k)
      fbu_extents(plon2+i,1,k) = fbu(      i,2,k)
      fbv_extents(      i,1,k) = fbv(plon2+i,2,k)
      fbv_extents(plon2+i,1,k) = fbv(      i,2,k)
    end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine extyv
  !=====================================================================


  !=====================================================================
  subroutine lcdbas(idim, grd, dgrd, dp, idp, dbas2, dbas3, hb, ht, dhb, dht)
    !
    ! Evaluate cubic basis functions for a 4-point stencil
    !===================================================================
    implicit none
    !
    ! Passed variables 
    !--------------------
    integer :: idim       ! dimension of input grid
    real(r8):: grd (idim) ! input grid
    real(r8):: dgrd(idim) ! input grid intervals
    real(r8):: dp         ! target grid point
    integer :: idp        ! pointer into grid for target grid pt
    real(r8):: dbas2(4)   ! derivatives at grid point 2.
    real(r8):: dbas3(4)   ! derivatives at grid point 3.
    real(r8):: hb         ! Interpolation weight for bot field value
    real(r8):: ht         ! Interpolation weight for top field value
    real(r8):: dhb        ! Interpolation weight for bot derivative
    real(r8):: dht        ! Interpolation weight for top derivative
    !
    ! Local Values
    !-----------------
    integer :: i
    real(r8):: dx
    real(r8):: xb
    real(r8):: xt
    real(r8):: x1                   !  |
    real(r8):: x2                   !  |- grid values
    real(r8):: x3                   !  |
    real(r8):: x4                   !  |
    real(r8):: x1mx2                !  |
    real(r8):: x1mx3                !  |
    real(r8):: x1mx4                !  |- differences of grid values
    real(r8):: x2mx3                !  |
    real(r8):: x2mx4                !  |
    real(r8):: x3mx4                !  |

    idp = -1
    do i = 2,idim-2
      if((dp.ge.grd(i)).and.(dp.lt.grd(i+1))) then
        idp = i
      endif
    end do
    if(idp.lt.0) then
      write(iulog,*) 'Error in LCDBAS:  target grid point outside of input grid'
      call endrun
    endif

    x1 = grd(idp-1)
    x2 = grd(idp  )
    x3 = grd(idp+1)
    x4 = grd(idp+2)

    x1mx2 = x1 - x2
    x1mx3 = x1 - x3
    x1mx4 = x1 - x4
    x2mx3 = x2 - x3
    x2mx4 = x2 - x4
    x3mx4 = x3 - x4

    dbas2(1) =   x2mx3 * x2mx4 / ( x1mx2 * x1mx3 * x1mx4 )
    dbas2(2) =   -1./x1mx2 + 1./x2mx3 + 1./x2mx4
    dbas2(3) = - x1mx2 * x2mx4 / ( x1mx3 * x2mx3 * x3mx4 )
    dbas2(4) =   x1mx2 * x2mx3 / ( x1mx4 * x2mx4 * x3mx4 )

    dbas3(1) = - x2mx3 * x3mx4 / ( x1mx2 * x1mx3 * x1mx4 )
    dbas3(2) =   x1mx3 * x3mx4 / ( x1mx2 * x2mx3 * x2mx4 )
    dbas3(3) =   -1./x1mx3 - 1./x2mx3 + 1./x3mx4
    dbas3(4) = - x1mx3 * x2mx3 / ( x1mx4 * x2mx4 * x3mx4 )

    dx       = dgrd(idp)
    xb       = ( grd(idp+1) - dp )/dx
    xt       = 1.0_r8 - xb
    hb       = ( 3.0_r8 - 2.0_r8*xb )*xb**2
    ht       = ( 3.0_r8 - 2.0_r8*xt )*xt**2
    dhb      = -dx*( xb - 1.0_r8  )*xb**2
    dht      =  dx*( xt - 1.0_r8  )*xt**2

    ! End Routine
    !------------
    return
  end subroutine lcdbas
  !=====================================================================


  !=====================================================================
  subroutine lcdint(x, lbasdx, rdx, limdr, hb, ht, dhb, dht, y)
    !
    ! Interpolate data on 4-pt stencil using cubic basis functions
    !===================================================================
    implicit none
    !
    ! Passed variables
    !---------------------
    real(r8)::  x(4)        ! input field on 4-pt stencil
    real(r8)::  lbasdx(4,2) ! basis functions
    real(r8)::  rdx         ! grid interval
    logical :: limdr       ! Flag to invoke SCM0 derivative limiting
    real(r8)::  hb          ! Interpolation weight for bot field value
    real(r8)::  ht          ! Interpolation weight for top field value
    real(r8)::  dhb         ! Interpolation weight for bot derivative
    real(r8)::  dht         ! Interpolation weight for top derivative
    real(r8)::  y           ! Interpolated value
    !
    ! Local values
    !-----------------
    real(r8):: fbot, ftop
    real(r8):: deli, tmp1, tmp2
    real(r8):: fac

    fbot = lbasdx(1,1)*x(1) &
          +lbasdx(2,1)*x(2) &
          +lbasdx(3,1)*x(3) &
          +lbasdx(4,1)*x(4)
    ftop = lbasdx(1,2)*x(1) &
          +lbasdx(2,2)*x(2) &
          +lbasdx(3,2)*x(3) &
          +lbasdx(4,2)*x(4)

    ! Apply SCM0 limiter to derivative estimates.
    !--------------------------------------------
    if(limdr) then
      fac  = 3._r8*(1._r8 - 10._r8*epsilon(fac))
      deli = ( x(3) - x(2) )*rdx
      tmp1 = fac*deli
      tmp2 = abs( tmp1 )
      if( deli*fbot   .le. 0.0_r8) fbot = 0._r8
      if( deli*ftop   .le. 0.0_r8) ftop = 0._r8
      if( abs( fbot ) .gt. tmp2  ) fbot = tmp1
      if( abs( ftop ) .gt. tmp2  ) ftop = tmp1
    endif

    y = x(2)*hb + fbot*dhb + x(3)*ht + ftop*dht

    ! End Routine
    !------------
    return
  end subroutine lcdint
  !=====================================================================


  !=====================================================================
  subroutine tsadj_rect(plat, plon, phis_old, phis_new, ts)
    !
    ! Adjust Ts based on difference between old and new phis.
    !===================================================================
    implicit none
    !
    ! Passed variables
    !--------------------
    integer :: plat                !  latitude dimension
    integer :: plon                !  longitude dimension
    real(r8):: phis_old(plon,plat) ! analysis phis (e.g., ECMWF)
    real(r8):: phis_new(plon,plat) ! model phis
    real(r8):: ts      (plon,plat) ! Surface Temp
    !
    ! Local values
    !---------------
    real(r8):: dtdz
    real(r8):: gravit
    real(r8):: del_z
    integer :: i, j, k           ! Indices

    dtdz    = -0.0065_r8           ! -6.5 deg/km
    gravit  = 9.80616_r8           ! acceleration of gravity ~ m/s^2

    do j = 1,plat
    do i = 1,plon
      del_z = (phis_new(i,j) - phis_old(i,j))/gravit
      ts(i,j) = ts(i,j) + dtdz*del_z
    end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine tsadj_rect
  !=====================================================================


  !=====================================================================
  subroutine tsadj_se(pcol, phis_old, phis_new, ts)
    !
    ! Adjust Ts based on difference between old and new phis.
    !===================================================================
    implicit none
    !
    ! Passed variables
    !--------------------
    integer :: pcol
    real(r8):: phis_old(pcol)  ! analysis phis (e.g., ECMWF)
    real(r8):: phis_new(pcol)  ! model phis
    real(r8):: ts      (pcol)  ! Surface Temp
    !
    ! Local values
    !----------------
    real(r8):: dtdz
    real(r8):: gravit
    real(r8):: del_z
    integer :: n, k             ! Indices

    dtdz    = -0.0065_r8           ! -6.5 deg/km
    gravit  = 9.80616_r8           ! acceleration of gravity ~ m/s^2

    do n = 1,pcol
      del_z = ( phis_new(n) - phis_old(n) )/gravit
      ts(n) = ts(n) + dtdz*del_z
    end do

    return
  end subroutine tsadj_se
  !=====================================================================


  !=====================================================================
  subroutine vert_quad_opt1_rect(plevi, plevip1, plev, plat, plon, t_old, pressi_m, &
                                 pressi_i, presso_m, phis_old, ps_old, t_new, loglin)
    !
    ! Quadratic interpolation (designed for Temperature interpolation)
    !
    !                        (if "loglin" == 1: in P
    !                         if "loglin" /= 1: in ln(P) )
    !
    !  Above input top                 :  quadratic using top, levels 1 and 2 
    !                                     (top defined as 1.e-10 Pa for now).  
    !                                     Top value set to value at level 1
    !  Between levels 1 and "bot"      :  quadratic interp using 3 closest levels
    !  Between levels "bot" and surface:  linear interpolation using "Tbot"and "Tsurf".
    !  Below surface                   :  You don"t wanna know
    !===================================================================
    implicit none
    !
    ! Passed variables
    !--------------------
    integer :: plevi                        ! vertical dimension of input analysis fields
    integer :: plevip1                      ! "plevi+1 (vert dimension of input interfaces)
    integer :: plev                         ! vertical dimension of model fields
    integer :: plat                         ! latitude dimension
    integer :: plon                         ! longitude dimension
    real(r8)::  t_old   (plevi  ,plon,plat)  ! analysis tempertatures
    real(r8)::  pressi_m(plevi  ,plon,plat)  ! analysis pressures at all levels
    real(r8)::  pressi_i(plevip1,plon,plat)  ! analysis interface pressures
    real(r8)::  presso_m(plev   ,plon,plat)  ! model pressures (based on adjusted PS)
    real(r8)::  phis_old(plon,        plat)  ! analysis phis
    real(r8)::  ps_old  (plon,        plat)  ! analysis surface pressure
    real(r8)::  ps_new  (plon,        plat)  ! "adjusted" model surface pressure
    real(r8)::  t_new   (plev   ,plon,plat)  ! Interpolated Temperatures
    integer :: loglin                       ! interpolation flag
    !
    ! Local Values
    !----------------
    real(r8)::  tsurf
    real(r8)::  t0
    real(r8)::  t_ref1
    real(r8)::  t_ref2
    real(r8)::  t_ref3
    real(r8)::  t_ref3_top
    real(r8)::  t_ref3_bot
    real(r8)::  z_ref_top
    real(r8)::  z_ref_bot
    real(r8)::  tbot
    real(r8)::  pbot
    real(r8)::  psurf
    real(r8)::  dtdz
    real(r8)::  lapse
    real(r8)::  boltz
    real(r8)::  avogad
    real(r8)::  mwdair
    real(r8)::  rgas
    real(r8)::  rdair
    real(r8)::  gravit
    real(r8)::  x
    real(r8)::  threshold
    real(r8)::  tmp
    real(r8)::  z
    real(r8)::  z_min
    real(r8)::  z_incr
    real(r8)::  hkk
    real(r8)::  x1
    real(r8)::  x2
    real(r8)::  x3
    real(r8)::  p1
    real(r8)::  p2
    real(r8)::  p3
    real(r8)::  px
    real(r8)::  pt                       ! top input pressure (ghost)
    real(r8)::  xt                       ! top input value (linear extrap)
    real(r8)::  tmp1
    real(r8)::  tmp2
    real(r8)::  tmp3
    real(r8)::  tmpt
    real(r8)::  beta
    real(r8)::  zero
    integer :: i, j, k, kk, kkp1, kkp2 ! Indices
    integer :: k_bot

    zero   = 0._r8
    dtdz   = -0.0065_r8          ! -6.5 deg/km
    gravit = 9.80616_r8          ! acceleration of gravity ~ m/s^2
    boltz  = 1.38065d-23         ! boltzmann's constant ~ J/k/molecule
    avogad = 6.02214d26          ! avogadro's number ~ molecules/kmole
    mwdair = 28.966_r8           ! molecular weight dry air ~ kg/kmole
    rgas   = avogad*boltz        ! universal gas constant ~ J/k/kmole
    rdair  = rgas/mwdair         ! constant for dry air   ~ J/k/kg

    t_ref1    = 290.5_r8
    t_ref2    = 255.0_r8
    t_ref3    = 298.0_r8
    z_ref_bot = 2000._r8
    z_ref_top = 2500._r8

    threshold = 0.001_r8
    pt        = 1.d-10
    if(loglin.ne.1) pt = log(pt)
      
    do j = 1,plat
    do i = 1,plon

      ! Tbot and Pbot are determined from the first model level that is at
      ! least 150m above the surface  ... will be used later
      !--------------------------------------------------------------------
      z_min = 150._r8
      z     = 0._r8

      do k = plevi,1,-1
        k_bot  = k
        hkk    = 0.5_r8*(pressi_i(k+1,i,j)-pressi_i(k,i,j))/pressi_m(k,i,j)
        z_incr = (rdair/gravit)*t_old(k,i,j)*hkk
        z      = z + z_incr
        if(z.gt.z_min) go to 10
        z = z + z_incr
      end do
      write(iulog,*) 'Error:  could not find model level above ',z_min
      call endrun
 10   continue

      lapse = -dtdz
      tbot  = t_old  (k_bot,i,j)
      pbot  = pressi_m(k_bot,i,j)
      tmp   = lapse*(rdair/gravit)*(ps_old(i,j)/pbot - 1._r8)
      tsurf = tbot*(1._r8 + tmp)

      ! Find bracketting input levels if possible
      !-------------------------------------------
      kk   = 1
      kkp1 = kk + 1
      kkp2 = kk + 2
      do k = 1,plev
 20     continue
          if(pressi_m(kkp1,i,j) .le. presso_m(k,i,j) ) then
            beta = presso_m(k,i,j) - pressi_m(kkp1,i,j)
            beta = beta/( pressi_m(kkp2,i,j) - pressi_m(kkp1,i,j) )
            if(beta .ge. 0.5_r8) then
              kk   = kk + 1
              kkp1 = kk + 1
              kkp2 = kk + 2
              if(kkp2 .ge. plevi+1) then
                kk   = plevi - 2
                kkp1 = kk + 1
                kkp2 = kk + 2
                go to 30
              endif
              go to 20
            endif
          endif
 30     continue
            
        x1    = t_old   (kk   ,i,j)
        x2    = t_old   (kkp1 ,i,j)
        x3    = t_old   (kkp2 ,i,j)
        p1    = pressi_m(kk   ,i,j)
        p2    = pressi_m(kkp1 ,i,j)
        p3    = pressi_m(kkp2 ,i,j)
        pbot  = pressi_m(k_bot,i,j)
        psurf = ps_old  (      i,j)
        px    = presso_m(k    ,i,j)

        ! Convert to log(P) for log(P) interpolation
        !--------------------------------------------
        if(loglin .ne. 1) then
          p1    = log(p1)
          p2    = log(p2)
          p3    = log(p3)
          pbot  = log(pbot)
          psurf = log(psurf)
          px    = log(px)
        endif

        if(px.lt.p1) then
          ! If above 1st analysis level: quadratic interp
          !-----------------------------------------------
          xt           = x1
          tmpt         = ((px-p1)*(px-p2))/((pt-p1)*(pt-p2))
          tmp1         = ((px-pt)*(px-p2))/((p1-pt)*(p1-p2))
          tmp2         = ((px-pt)*(px-p1))/((p2-pt)*(p2-p1))
          t_new(k,i,j) = xt*tmpt + x1*tmp1 + x2*tmp2
        elseif((px.ge.pbot).and.(px.le.psurf)) then
          ! Elseif between "pbot" and analysis surface: linear interp
          !------------------------------------------------------------
          t_new(k,i,j) = (tbot*(psurf-px) + tsurf*(px-pbot))/(psurf-pbot)
        elseif(px .gt. psurf ) then
          ! Elseif below analysis surface: special case - heuristic equations
          !-------------------------------------------------------------------
          t_ref3_bot =      tsurf - z_ref_bot*dtdz
          t_ref3_top = min( tsurf - z_ref_top*dtdz, t_ref3 )
          z = phis_old(i,j)/gravit

          if(z .ge. z_ref_bot) then
            if(z .ge. z_ref_top) then
              t0 = min((tsurf-z*dtdz),t_ref3)
            else
              t0 = ( t_ref3_bot*(z_ref_top - z        ) &
                    +t_ref3_top*(z         - z_ref_bot) )/(z_ref_top-z_ref_bot)
            endif
            lapse = max (((t0-tsurf)/z),zero)
          else
            lapse = -dtdz
          endif

          x   = lapse*rdair/gravit*log(presso_m(k,i,j)/ps_old(i,j))
          t_new(k,i,j) = tsurf*(1._r8 + x + x*x/2._r8 + x*x*x/6._r8)
        elseif(px.ge.p3) then
          ! Should never happen
          !----------------------
          write(iulog,*) 'Error:  extrapolation below input levels not allowed'
          call endrun
        else
          ! Else between 1st analysis level and "pbot":  quadratic interp
          !---------------------------------------------------------------
          tmp1         = ((px-p2)*(px-p3))/((p1-p2)*(p1-p3))
          tmp2         = ((px-p1)*(px-p3))/((p2-p1)*(p2-p3))
          tmp3         = ((px-p1)*(px-p2))/((p3-p1)*(p3-p2))
          t_new(k,i,j) = x1*tmp1 + x2*tmp2 + x3*tmp3
        endif
      end do ! k = 1,plev
    end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine vert_quad_opt1_rect
  !=====================================================================


  !=====================================================================
  subroutine vert_quad_opt1_se(plevi, plevip1, plev, pcol, t_old,                           &
                               pressi_m, pressi_i, presso_m, phis_old, ps_old, t_new, loglin)
    !
    ! Quadratic interpolation (designed for Temperature interpolation)
    !
    !                        (if "loglin" == 1: in P
    !                         if "loglin" /= 1: in ln(P) )
    !
    !  Above input top                 :  quadratic using top, levels 1 and 2 
    !                                     (top defined as 1.e-10 Pa for now).  
    !                                     Top value set to value at level 1
    !  Between levels 1 and "bot"      :  quadratic interp using 3 closest levels
    !  Between levels "bot" and surface:  linear interpolation using "Tbot"and "Tsurf".
    !  Below surface                   :  You don"t wanna know
    !===================================================================
    implicit none
    !
    ! Passed Variables
    !------------------
    integer :: plevi                  ! vertical dimension of input analysis fields
    integer :: plevip1                ! "plevi+1 (vert dimension of input interfaces)
    integer :: plev                   ! vertical dimension of model fields
    integer :: pcol
    real(r8):: t_old   (plevi  ,pcol)  ! analysis tempertatures
    real(r8):: pressi_m(plevi  ,pcol)  ! analysis pressures at all levels
    real(r8):: pressi_i(plevip1,pcol)  ! analysis interface pressures
    real(r8):: presso_m(plev   ,pcol)  ! model pressures (based on adjusted PS)
    real(r8):: phis_old(pcol)          ! analysis phis
    real(r8):: ps_old  (pcol)          ! analysis surface pressure
    real(r8):: ps_new  (pcol)          ! "adjusted" model surface pressure
    real(r8):: t_new   (plev,pcol)     ! Interpolated Temperatures
    integer :: loglin                 ! interpolation flag
    !
    ! Local Values
    !---------------
    real(r8):: tsurf
    real(r8):: t0
    real(r8):: t_ref1
    real(r8):: t_ref2
    real(r8):: t_ref3
    real(r8):: t_ref3_top
    real(r8):: t_ref3_bot
    real(r8):: z_ref_top
    real(r8):: z_ref_bot
    real(r8):: tbot
    real(r8):: pbot
    real(r8):: psurf
    real(r8):: dtdz
    real(r8):: lapse
    real(r8):: boltz
    real(r8):: avogad
    real(r8):: mwdair
    real(r8):: rgas
    real(r8):: rdair
    real(r8):: gravit
    real(r8):: x
    real(r8):: threshold
    real(r8):: tmp
    real(r8):: z
    real(r8):: z_min
    real(r8):: z_incr
    real(r8):: hkk
    real(r8):: x1
    real(r8):: x2
    real(r8):: x3
    real(r8):: p1
    real(r8):: p2
    real(r8):: p3
    real(r8):: px
    real(r8):: pt                       ! top input pressure (ghost)
    real(r8):: xt                       ! top input value (linear extrap)
    real(r8):: tmp1
    real(r8):: tmp2
    real(r8):: tmp3
    real(r8):: tmpt
    real(r8):: beta
    real(r8):: zero
    integer :: n, k, kk, kkp1, kkp2 ! Indices
    integer :: k_bot

    zero   = 0._r8
    dtdz   = -0.0065_r8         ! -6.5 deg/km
    gravit = 9.80616_r8         ! acceleration of gravity ~ m/s^2
    boltz  = 1.38065d-23        ! boltzmann's constant ~ J/k/molecule
    avogad = 6.02214d26         ! avogadro's number ~ molecules/kmole
    mwdair = 28.966_r8          ! molecular weight dry air ~ kg/kmole
    rgas   = avogad*boltz       ! universal gas constant ~ J/k/kmole
    rdair  = rgas/mwdair        ! constant for dry air   ~ J/k/kg

    t_ref1    = 290.5_r8
    t_ref2    = 255.0_r8
    t_ref3    = 298.0_r8
    z_ref_bot = 2000._r8
    z_ref_top = 2500._r8

    threshold = 0.001_r8
    pt        = 1.d-10
    if(loglin.ne.1) pt = log(pt)
      
    do n = 1,pcol

      ! Tbot and Pbot are determined from the first model level that is at
      ! least 150m above the surface ... will be used later
      !---------------------------------------------------------------------
      z_min = 150._r8
      z     = 0._r8

      do k = plevi,1,-1
        k_bot  = k
        hkk    = 0.5_r8*(pressi_i(k+1,n)-pressi_i(k,n))/pressi_m(k,n)
        z_incr = (rdair/gravit)*t_old(k,n)*hkk
        z      = z + z_incr
        if(z.gt.z_min) go to 10
        z = z + z_incr
      end do
      write(iulog,*) 'Error:  could not find model level above ',z_min
      call endrun
 10   continue

      lapse = -dtdz
      tbot  = t_old  (k_bot,n)
      pbot  = pressi_m(k_bot,n)
      tmp   = lapse*(rdair/gravit)*(ps_old(n)/pbot - 1._r8)
      tsurf = tbot*(1._r8 + tmp)

      ! Find bracketting input levels if possible
      !-------------------------------------------
      kk   = 1
      kkp1 = kk + 1
      kkp2 = kk + 2
      do k = 1,plev
 20     continue
          if(pressi_m(kkp1,n) .le. presso_m(k,n) ) then
            beta = presso_m(k,n) - pressi_m(kkp1,n)
            beta = beta/( pressi_m(kkp2,n) - pressi_m(kkp1,n) )
            if(beta .ge. 0.5) then
              kk   = kk + 1
              kkp1 = kk + 1
              kkp2 = kk + 2
              if(kkp2 .ge. plevi+1) then
                kk   = plevi - 2
                kkp1 = kk + 1
                kkp2 = kk + 2
                go to 30
              endif
              go to 20
            endif
          endif
 30     continue
            
        x1    = t_old   (kk   ,n)
        x2    = t_old   (kkp1 ,n)
        x3    = t_old   (kkp2 ,n)
        p1    = pressi_m(kk   ,n)
        p2    = pressi_m(kkp1 ,n)
        p3    = pressi_m(kkp2 ,n)
        pbot  = pressi_m(k_bot,n)
        psurf = ps_old  (      n)
        px    = presso_m(k    ,n)

        ! Convert to log(P) for log(P) interpolation
        !---------------------------------------------
        if(loglin.ne.1) then
          p1    = log(p1)
          p2    = log(p2)
          p3    = log(p3)
          pbot  = log(pbot)
          psurf = log(psurf)
          px    = log(px)
        endif

        if(px.lt.p1 ) then
          ! If above 1st analysis level: quadratic interp
          !------------------------------------------------
          xt         = x1
          tmpt       = ((px-p1)*(px-p2))/((pt-p1)*(pt-p2))
          tmp1       = ((px-pt)*(px-p2))/((p1-pt)*(p1-p2))
          tmp2       = ((px-pt)*(px-p1))/((p2-pt)*(p2-p1))
          t_new(k,n) = xt*tmpt + x1*tmp1 + x2*tmp2
        elseif((px.ge.pbot).and.(px.le.psurf)) then
          ! Elseif between "pbot" and analysis surface: linear interp
          !-----------------------------------------------------------
          t_new(k,n) = (tbot*(psurf-px) + tsurf*(px-pbot))/(psurf-pbot)
        elseif(px.gt.psurf) then
          ! Elseif below analysis surface: special case - heuristic equations
          !------------------------------------------------------------------
          t_ref3_bot =      tsurf - z_ref_bot*dtdz
          t_ref3_top = min( tsurf - z_ref_top*dtdz, t_ref3 )
          z = phis_old(n)/gravit

          if(z .ge. z_ref_bot) then
            if(z .ge. z_ref_top) then
              t0 = min( tsurf - z*dtdz, t_ref3 )
            else
              t0 = ( t_ref3_bot*(z_ref_top - z        ) &
                    +t_ref3_top*(z         - z_ref_bot) )/(z_ref_top-z_ref_bot)
            endif
            lapse = max ( (t0 - tsurf)/z , zero)
          else
            lapse = -dtdz
          endif

          x = lapse*rdair/gravit*log(presso_m(k,n)/ps_old(n))
          t_new(k,n) = tsurf*(1._r8 + x + x*x/2._r8 + x*x*x/6._r8)
        elseif(px.ge.p3) then
          ! Should never happen
          !--------------------
          write(iulog,*) 'Error:  extrapolation below input levels not allowed'
          call endrun
        else
          ! Else between 1st analysis level and "pbot":  quadratic interp
          !---------------------------------------------------------------
          tmp1       = ( (px-p2)*(px-p3) )/( (p1-p2)*(p1-p3) )
          tmp2       = ( (px-p1)*(px-p3) )/( (p2-p1)*(p2-p3) )
          tmp3       = ( (px-p1)*(px-p2) )/( (p3-p1)*(p3-p2) )
          t_new(k,n) = x1*tmp1 + x2*tmp2 + x3*tmp3
        endif
      end do ! k = 1,plev
    end do ! n = 1,pcol

    ! End Routine
    !-------------
    return
  end subroutine vert_quad_opt1_se
  !=====================================================================


  !=====================================================================
  subroutine vert_int_opt2_rect(plat, plon, plevi, plevo, pressi, presso, xxi, xxo, loglin)
    !
    ! Designed for vertical interpolation of U/V
    !
    ! Linear and Quadratic interpolation (if "loglin" == 1: in P
    !                                     if "loglin" /= 1: in ln(P) )
    !
    !  Above input top        :  quadratic using top, levels 1 and 2 (top
    !                            defined as 1.e-10 Pa for now).  Top value
    !                            determined from linear extrapolation from
    !                            levels 1 and 2
    !  Between levels 1 and 2 :  quadratic interp using levels 1,2, & 3
    !  Between levels 2 and K :  linear interpolation using adjacent levels
    !  Below level K          :  set equal to level K
    !===================================================================
    implicit none
    !
    ! Passed variables
    !-------------------
    integer :: plat                    ! latitude dimension
    integer :: plon                    ! longitude dimension
    integer :: plevi                   ! vertical dimension of analysis fields
    integer :: plevo                   ! vertical dimension of model fields
    real(r8):: pressi(plevi,plon,plat) ! analysis pressures
    real(r8):: presso(plevo,plon,plat) ! model pressures (based on adjusted PS)
    real(r8):: xxi   (plevi,plon,plat) ! input analysis field
    real(r8):: xxo   (plevo,plon,plat) ! model field
    integer :: loglin                  ! interpolation flag
    !
    ! Local values
    !----------------
    real(r8):: x1
    real(r8):: x2
    real(r8):: x3
    real(r8):: p1
    real(r8):: p2
    real(r8):: p3
    real(r8):: px
    real(r8):: pt                       ! top input pressure (ghost)
    real(r8):: xt                       ! top input value (linear extrap)
    real(r8):: tmp1
    real(r8):: tmp2
    real(r8):: tmp3
    real(r8):: tmpt
    integer :: i, j, k, kk, kkp1, kkp2 ! Indices

    pt = 1.d-10
    if(loglin.ne.1) pt = log(pt)

    do j = 1,plat
    do i = 1,plon

      ! Find bracketting analysis pressure levels
      !-------------------------------------------
      kk   = 1
      kkp1 = kk + 1
      kkp2 = kk + 2
      do k = 1,plevo
 10     continue
          if(pressi(kkp1,i,j) .le. presso(k,i,j) ) then
            kk   = kk + 1
            kkp1 = kk + 1
            kkp2 = kk + 2
            if(kkp1 .eq. plevi+1) then
              kk   = plevi - 1
              kkp1 = kk + 1
              kkp2 = kk + 2
              go to 20
            endif
            go to 10
          endif
 20     continue
            
        x1 = xxi   (kk  ,i,j)
        x2 = xxi   (kkp1,i,j)
        p1 = pressi(kk  ,i,j)
        p2 = pressi(kkp1,i,j)
        px = presso(k   ,i,j)

        if(kkp2 .le. plevi) then
          x3 = xxi   (kkp2,i,j)
          p3 = pressi(kkp2,i,j)
        else
          x3 = 1.d+36
          p3 = 1.d+36
        endif

        ! Convert to log(P) for log(P) interpolation
        !--------------------------------------------
        if(loglin .ne. 1) then
          p1 = log(p1)
          p2 = log(p2)
          p3 = log(p3)
          px = log(px)
        endif

        if(px.lt.p1) then
          ! If above 1st analysis level:  quadratic interp
          !------------------------------------------------
          xt         = ( x1*(p2 - pt) - x2*(p1 - pt) )/(p2 - p1)
          tmpt       = ( (px-p1)*(px-p2) )/( (pt-p1)*(pt-p2) )
          tmp1       = ( (px-pt)*(px-p2) )/( (p1-pt)*(p1-p2) )
          tmp2       = ( (px-pt)*(px-p1) )/( (p2-pt)*(p2-p1) )
          xxo(k,i,j) = xt*tmpt + x1*tmp1 + x2*tmp2
        elseif(px.ge.p2) then
          ! Elseif below bottom analysis level:  output = bottome analysis field value
          !--------------------------------------------------------------------------
          xxo(k,i,j) = x2
        elseif(kk.eq.1) then
          ! Elseif between 1st and 2nd analysis levels:  quadratic interp
          !---------------------------------------------------------------
          tmp1       = ( (px-p2)*(px-p3) )/( (p1-p2)*(p1-p3) )
          tmp2       = ( (px-p1)*(px-p3) )/( (p2-p1)*(p2-p3) )
          tmp3       = ( (px-p1)*(px-p2) )/( (p3-p1)*(p3-p2) )
          xxo(k,i,j) = x1*tmp1 + x2*tmp2 + x3*tmp3
        else
          ! Else, Linear interpolation
          !-----------------------------
          xxo(k,i,j) = ( x1*(p2 - px) + x2*(px - p1) )/(p2 - p1)
        endif
      end do ! k = 1,plevo

    end do
    end do

    ! End Rotine
    !-------------
    return
  end subroutine vert_int_opt2_rect
  !=====================================================================


  !=====================================================================
  subroutine vert_int_opt1_rect(plat, plon, plevi, plevo, pressi, presso, xxi, xxo, loglin)
    !
    ! Designed for moisture fields like q, cloud water, cloud ice, cloud frac, etc.
    !
    ! Linearly interpolate (if "loglin" == 1: in P
    !                       if "loglin" /= 1: in ln(P) )
    !
    !  Above input top        :  set equal to level 1
    !  Between levels 1 and K :  linear interpolation using adjacent levels
    !  Below level K          :  set equal to level K
    !===================================================================
    implicit none
    !
    ! Passed variables 
    !-------------------
    integer :: plat                     ! latitude dimension
    integer :: plon                     ! longitude dimension
    integer :: plevi                    ! vertical dimension of analysis fields
    integer :: plevo                    ! vertical dimension of model fields
    real(r8):: pressi(plevi,plon,plat)  ! analysis pressures
    real(r8):: presso(plevo,plon,plat)  ! model pressures (based on adjusted PS)
    real(r8):: xxi   (plevi,plon,plat)  ! analysis field
    real(r8):: xxo   (plevo,plon,plat)  ! model field
    integer :: loglin                   ! interpolation flag
    !
    ! Local values
    !---------------
    real(r8):: p1
    real(r8):: p2
    real(r8):: px
    integer :: i, j, k, kk, kkp1         ! Indices

    do j = 1,plat
    do i = 1,plon

      ! Find bracketting analysis pressure levels
      !--------------------------------------------
      kk   = 1
      kkp1 = kk + 1
      do k = 1,plevo
 10     continue
          if(pressi(kkp1,i,j) .le. presso(k,i,j) ) then
            kk   = kk + 1
            kkp1 = kk + 1
            if(kkp1 .eq. plevi+1) then
              kk   = plevi - 1
              kkp1 = kk + 1
              go to 20
            endif
            go to 10
          endif
 20     continue

        if(presso(k,i,j).lt.pressi(kk,i,j)) then
          ! If above 1st analysis level:  output = top analysis field value
          !----------------------------------------------------------------
          xxo(k,i,j) = xxi(kk  ,i,j)
        elseif(presso(k,i,j).ge.pressi(kkp1,i,j)) then
          ! If below bottom analysis level:  output = bottom analysis field value
          !-----------------------------------------------------------------------
          xxo(k,i,j) = xxi(kkp1,i,j)
        else
          ! Else, Linear interpolation
          !-----------------------------
          p1 = pressi(kk  ,i,j)
          p2 = pressi(kkp1,i,j)
          px = presso(k   ,i,j)
          if(loglin .ne. 1) then
            p1 = log(p1)
            p2 = log(p2)
            px = log(px)
          endif
          xxo(k,i,j) = xxi(kk,i,j)*(p2-px) + xxi(kkp1,i,j)*(px-p1)
          xxo(k,i,j) = xxo(k ,i,j)/(p2-p1)
        endif
      end do

    end do
    end do

    ! End Routine
    !------------
    return
  end subroutine vert_int_opt1_rect
  !=====================================================================


  !=====================================================================
  subroutine vert_int_opt2_se(pcol, plevi, plevo, pressi, presso, xxi, xxo, loglin)
    !
    ! Designed for vertical interpolation of U/V
    !
    ! Linear and Quadratic interpolation (if "loglin" == 1: in P
    !                                     if "loglin" /= 1: in ln(P) )
    !
    !  Above input top        :  quadratic using top, levels 1 and 2 (top
    !                            defined as 1.e-10 Pa for now).  Top value
    !                            determined from linear extrapolation from
    !                            levels 1 and 2
    !  Between levels 1 and 2 :  quadratic interp using levels 1,2, & 3
    !  Between levels 2 and K :  linear interpolation using adjacent levels
    !  Below level K          :  set equal to level K
    !===================================================================
    implicit none
    !
    ! Passed variables 
    !--------------------
    integer :: pcol
    integer :: plevi              ! vertical dimension of analysis fields
    integer :: plevo              ! vertical dimension of model fields
    real(r8)::  pressi(plevi,pcol) ! analysis pressures
    real(r8)::  presso(plevo,pcol) ! model pressures (based on adjusted PS)
    real(r8)::  xxi   (plevi,pcol) ! input analysis field
    real(r8)::  xxo   (plevo,pcol) ! model field
    integer :: loglin             ! interpolation flag
    !
    ! Local values
    !----------------
    real(r8):: x1
    real(r8):: x2
    real(r8):: x3
    real(r8):: p1
    real(r8):: p2
    real(r8):: p3
    real(r8):: px
    real(r8):: pt                       ! top input pressure (ghost)
    real(r8):: xt                       ! top input value (linear extrap)
    real(r8):: tmp1
    real(r8):: tmp2
    real(r8):: tmp3
    real(r8):: tmpt
    integer :: n, k, kk, kkp1, kkp2 ! Indices

    pt = 1.d-10
    if(loglin.ne.1) pt = log(pt)

    do n = 1,pcol

      ! Find bracketting analysis pressure levels
      !-------------------------------------------
      kk   = 1
      kkp1 = kk + 1
      kkp2 = kk + 2
      do k = 1,plevo
 10     continue
          if(pressi(kkp1,n) .le. presso(k,n) ) then
            kk   = kk + 1
            kkp1 = kk + 1
            kkp2 = kk + 2
            if(kkp1 .eq. plevi+1) then
              kk   = plevi - 1
              kkp1 = kk + 1
              kkp2 = kk + 2
              go to 20
            endif
            go to 10
          endif
 20     continue
            
        x1 = xxi   (kk  ,n)
        x2 = xxi   (kkp1,n)
        p1 = pressi(kk  ,n)
        p2 = pressi(kkp1,n)
        px = presso(k   ,n)

        if(kkp2 .le. plevi) then
          x3 = xxi   (kkp2,n)
          p3 = pressi(kkp2,n)
        else
          x3 = 1.d+36
          p3 = 1.d+36
        endif

        ! Convert to log(P) for log(P) interpolation
        !--------------------------------------------
        if(loglin .ne. 1) then
          p1 = log(p1)
          p2 = log(p2)
          p3 = log(p3)
          px = log(px)
        endif

        if(px.lt.p1) then
          ! If above 1st analysis level:  quadratic interp
          !--------------------------------------------------
          xt       = ( x1*(p2 - pt) - x2*(p1 - pt) )/(p2 - p1)
          tmpt     = ( (px-p1)*(px-p2) )/( (pt-p1)*(pt-p2) )
          tmp1     = ( (px-pt)*(px-p2) )/( (p1-pt)*(p1-p2) )
          tmp2     = ( (px-pt)*(px-p1) )/( (p2-pt)*(p2-p1) )
          xxo(k,n) = xt*tmpt + x1*tmp1 + x2*tmp2
        elseif(px.ge.p2) then
          ! Elseif below bottom analysis level:  output = bottome analysis field value
          !----------------------------------------------------------------------------
          xxo(k,n) = x2
        elseif(kk.eq.1) then
          ! Elseif between 1st and 2nd analysis levels:  quadratic interp
          !---------------------------------------------------------------
          tmp1     = ( (px-p2)*(px-p3) )/( (p1-p2)*(p1-p3) )
          tmp2     = ( (px-p1)*(px-p3) )/( (p2-p1)*(p2-p3) )
          tmp3     = ( (px-p1)*(px-p2) )/( (p3-p1)*(p3-p2) )
          xxo(k,n) = x1*tmp1 + x2*tmp2 + x3*tmp3
        else
          ! Else, Linear interpolation
          !-----------------------------
          xxo(k,n) = ( x1*(p2 - px) + x2*(px - p1) )/(p2 - p1)
        endif
      end do ! k = 1,plevo

    end do

    ! End Rotine
    !-----------
    return
  end subroutine vert_int_opt2_se
  !=====================================================================


  !=====================================================================
  subroutine vert_int_opt1_se(pcol, plevi, plevo, pressi, presso, xxi, xxo, loglin)
    !
    ! Designed for moisture fields like q, cloud water, cloud ice, cloud frac, etc.
    !
    ! Linearly interpolate (if "loglin" == 1: in P
    !                       if "loglin" /= 1: in ln(P) )
    !
    !  Above input top        :  set equal to level 1
    !  Between levels 1 and K :  linear interpolation using adjacent levels
    !  Below level K          :  set equal to level K
    !===================================================================
    implicit none
    !
    ! Passed variables 
    !--------------------
    integer :: pcol
    integer :: plevi               ! vertical dimension of analysis fields
    integer :: plevo               ! vertical dimension of model fields
    real(r8)::  pressi(plevi,pcol)  ! analysis pressures
    real(r8)::  presso(plevo,pcol)  ! model pressures (based on adjusted PS)
    real(r8)::  xxi   (plevi,pcol)  ! analysis field
    real(r8)::  xxo   (plevo,pcol)  ! model field
    integer :: loglin              ! interpolation flag
    !
    ! Local values
    !---------------
    real(r8):: p1
    real(r8):: p2
    real(r8):: px
    integer :: n, k, kk, kkp1         ! Indices

    do n = 1,pcol

      ! Find bracketting analysis pressure levels
      !-------------------------------------------
      kk   = 1
      kkp1 = kk + 1
      do k = 1,plevo
 10     continue
          if(pressi(kkp1,n) .le. presso(k,n) ) then
            kk   = kk + 1
            kkp1 = kk + 1
            if(kkp1 .eq. plevi+1) then
              kk   = plevi - 1
              kkp1 = kk + 1
              go to 20
            endif
            go to 10
          endif
 20     continue

        if(presso(k,n).lt.pressi(kk,n)) then
          ! If above 1st analysis level:  output = top analysis field value
          !-----------------------------------------------------------------
          xxo(k,n) = xxi(kk,n)
        elseif(presso(k,n) .ge. pressi(kkp1,n) ) then
          ! If below bottom analysis level:  output = bottom analysis field value
          !---------------------------------------------------------------------
          xxo(k,n) = xxi(kkp1,n)
        else
          ! Else, Linear interpolation
          !----------------------------
          p1 = pressi(kk  ,n)
          p2 = pressi(kkp1,n)
          px = presso(k   ,n)
          if(loglin .ne. 1) then
            p1 = log(p1)
            p2 = log(p2)
            px = log(px)
          endif
          xxo(k,n) = xxi(kk,n)*(p2-px) + xxi(kkp1,n)*(px-p1)
          xxo(k,n) = xxo(k ,n)/(p2-p1)
        endif
      end do ! k = 1,plevo

    end do

    ! End Routine
    !-------------
    return
  end subroutine vert_int_opt1_se
  !=====================================================================


  !=====================================================================
  subroutine myminmax_rect(plev, plat, plon, x, fmin, fmax)
    !
    ! Bracket "x" between fmin and fmax
    !===================================================================
    implicit none
    !
    ! Passed Variables
    !-------------------
    integer :: plev
    integer :: plat
    integer :: plon
    real(r8)::  x(plon,plat,plev)
    real(r8)::  fmin
    real(r8)::  fmax
    !
    ! Local values
    !--------------
    integer:: i, j, k           ! Indices

    do k = 1,plev
    do j = 1,plat
    do i = 1,plon
      x(i,j,k) = min(max(x(i,j,k),fmin),fmax)
    end do
    end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine myminmax_rect
  !=====================================================================


  !=====================================================================
  subroutine myminmax_se(plev, pcol, x, fmin, fmax)
    !
    ! Bracket "x" between fmin and fmax
    !===================================================================
    implicit none
    !
    ! Passed variables
    !-------------------
    integer :: plev
    integer :: pcol
    real(r8)::  x(pcol,plev)
    real(r8)::  fmin
    real(r8)::  fmax
    !
    ! Local values 
    !--------------
    integer:: n, k           ! Indices

    do k = 1,plev
    do n = 1,pcol
      x(n,k) = min(max(x(n,k),fmin),fmax)
    end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine myminmax_se
  !=====================================================================


  !=====================================================================
  subroutine psadj_rect(plev, plevp1, plat, plon, t,                        &
                        press_m, press_i, phis_old, phis_new, ps_old, ps_new)
    !
    ! Adjust Ps based on difference between "analysis" phis and model phis.
    ! Also uses T and P arrays
    !===================================================================
    implicit none
    !
    ! Passed variables 
    !---------------------
    integer :: plev                       ! vertical dimension
    integer :: plevp1                     ! "plev+1"
    integer :: plat                       ! latitude dimension
    integer :: plon                       ! longitude dimension
    real(r8):: t       (plon,plat,plev)   ! analysis Temperatures
    real(r8):: press_m (plon,plat,plev)   ! analysis pressures
    real(r8):: press_i (plon,plat,plevp1) ! analysis pressures (interfaces)
    real(r8):: phis_old(plon,     plat)   ! analysis phis
    real(r8):: phis_new(plon,     plat)   ! model phis
    real(r8):: ps_old  (plon,     plat)   ! analysis Ps (horizontal interpolated to model grid)
    real(r8):: ps_new  (plon,     plat)   ! adjusted model Ps
    !
    ! Local Values 
    !---------------
    real(r8):: tsurf
    real(r8):: t_ref1
    real(r8):: t_ref2
    real(r8):: t0
    real(r8):: tbot
    real(r8):: pbot
    real(r8):: dtdz
    real(r8):: lapse
    real(r8):: boltz
    real(r8):: avogad
    real(r8):: mwdair
    real(r8):: rgas
    real(r8):: rdair
    real(r8):: gravit
    real(r8):: del_phis
    real(r8):: x
    real(r8):: threshold
    real(r8):: tmp
    real(r8):: z
    real(r8):: z_min
    real(r8):: z_incr
    real(r8):: hkk
    integer :: i, j, k, kk         ! Indices

    dtdz   = -0.0065_r8        ! -6.5 deg/km
    gravit = 9.80616_r8        ! acceleration of gravity ~ m/s^2
    boltz  = 1.38065d-23       ! boltzmann's constant ~ J/k/molecule
    avogad = 6.02214d26        ! avogadro's number ~ molecules/kmole
    mwdair = 28.966_r8         ! molecular weight dry air ~ kg/kmole
    rgas   = avogad*boltz      ! universal gas constant ~ J/k/kmole
    rdair  = rgas/mwdair       ! constant for dry air   ~ J/k/kg

    t_ref1    = 290.5_r8
    t_ref2    = 255.0_r8
    threshold = 0.001_r8
      
    do j = 1,plat
    do i = 1,plon

      del_phis = phis_old(i,j) - phis_new(i,j)

      if(abs(del_phis) .le. threshold) then
        ! If difference between analysis and model phis is negligible,
        ! then set model Ps = analysis
        !-------------------------------------------------------------
        ps_new(i,j) = ps_old(i,j)
      else
        ! Else, go nuts...
        !
        ! Tbot and Pbot are determined from the first model level 
        ! that is at least 150m above the surface
        !------------------------------------------------------------
        z_min = 150._r8
        z     = 0._r8
        do k = plev,1,-1
          kk     = k
          hkk    = 0.5_r8*(press_i(i,j,k+1)-press_i(i,j,k))/press_m(i,j,k)
          z_incr = (rdair/gravit)*t(i,j,k)*hkk
          z      = z + z_incr
          if(z.gt.z_min) go to 10
          z = z + z_incr
        end do
        write(iulog,*) 'Error:  could not find model level above ',z_min
        call endrun
 10     continue
        lapse = -dtdz
        k     = kk

        ! Define Tbot & Pbot
        !---------------------
        tbot  = t      (i,j,k)
        pbot  = press_m(i,j,k)
        tmp   = lapse*(rdair/gravit)*(ps_old(i,j)/pbot - 1.)
        tsurf = tbot*(1._r8 + tmp)
        t0    = tsurf + lapse*phis_old(i,j)/gravit

        ! Based on heuristic equations:
        !------------------------------
        if((t0.gt.t_ref1).and.(tsurf.le.t_ref1)) then
          lapse = (t_ref1 - tsurf)*gravit/phis_old(i,j)
        elseif((t0.gt.t_ref1).and.(tsurf.gt.t_ref1)) then
          lapse = 0._r8
          tsurf = (t_ref1+tsurf)*0.5_r8
        endif

        if(tsurf .lt. t_ref2) then
          lapse = -dtdz
          tsurf = (t_ref2+tsurf)*0.5_r8
        endif              

        x   = lapse*del_phis/(gravit*tsurf)
        tmp = 1._r8 - x/2._r8 + x**2._r8/3._r8
        tmp = del_phis/(rdair*tsurf)*tmp
        ps_new(i,j) = ps_old(i,j)*exp(tmp)
      endif

    end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine psadj_rect
  !=====================================================================


  !=====================================================================
  subroutine psadj_se(plev, plevp1, pcol, t,                              &
                      press_m, press_i, phis_old, phis_new, ps_old, ps_new)
    !
    ! Adjust Ps based on difference between "analysis" phis and model phis.
    ! Also uses T and P arrays
    !===================================================================
    implicit none
    !
    ! Passed Variables 
    !--------------------
    integer :: plev                  ! vertical dimension
    integer :: plevp1                ! "plev+1"
    integer :: pcol      
    real(r8):: t       (pcol,plev)   ! analysis Temperatures
    real(r8):: press_m (pcol,plev)   ! analysis pressures
    real(r8):: press_i (pcol,plevp1) ! analysis pressures (interfaces)
    real(r8):: phis_old(pcol)        ! analysis phis
    real(r8):: phis_new(pcol)        ! model phis
    real(r8):: ps_old  (pcol)        ! analysis Ps (horizontal interpolated to model grid)
    real(r8):: ps_new  (pcol)        ! adjusted model Ps
    !
    ! Local Values 
    !------------------
    real(r8):: tsurf
    real(r8):: t_ref1
    real(r8):: t_ref2
    real(r8):: t0
    real(r8):: tbot
    real(r8):: pbot
    real(r8):: dtdz
    real(r8):: lapse
    real(r8):: boltz
    real(r8):: avogad
    real(r8):: mwdair
    real(r8):: rgas
    real(r8):: rdair
    real(r8):: gravit
    real(r8):: del_phis
    real(r8):: x
    real(r8):: threshold
    real(r8):: tmp
    real(r8):: z
    real(r8):: z_min
    real(r8):: z_incr
    real(r8):: hkk
    integer :: n, k, kk         ! Indices

    dtdz   = -0.0065_r8        ! -6.5 deg/km
    gravit = 9.80616_r8        ! acceleration of gravity ~ m/s^2
    boltz  = 1.38065d-23       ! boltzmann's constant ~ J/k/molecule
    avogad = 6.02214d26        ! avogadro's number ~ molecules/kmole
    mwdair = 28.966_r8         ! molecular weight dry air ~ kg/kmole
    rgas   = avogad*boltz      ! universal gas constant ~ J/k/kmole
    rdair  = rgas/mwdair       ! constant for dry air   ~ J/k/kg

    t_ref1    = 290.5_r8
    t_ref2    = 255.0_r8
    threshold = 0.001_r8
      
    do n = 1,pcol

      del_phis = phis_old(n) - phis_new(n)

      if(abs(del_phis).le.threshold) then
        ! If difference between analysis and model phis is negligible,
        ! then set model Ps = analysis
        !-------------------------------------------------------------
        ps_new(n) = ps_old(n)
      else
        ! Else, go nuts...
        !
        ! Tbot and Pbot are determined from the first model level that is at
        ! least 150m above the surface
        !--------------------------------------------------------------------
        z_min = 150._r8
        z     = 0._r8
        do k = plev,1,-1
          kk     = k
          hkk    = 0.5_r8*( press_i(n,k+1) - press_i(n,k))/press_m(n,k)
          z_incr = (rdair/gravit)*t(n,k)*hkk
          z      = z + z_incr
          if(z .gt. z_min) go to 10
          z      = z + z_incr
        end do
        write(iulog,*) 'Error:  could not find model level above ',z_min
        call endrun
 10     continue
        lapse = -dtdz
        k     = kk

        ! Define Tbot & Pbot
        !--------------------
        tbot  = t      (n,k)
        pbot  = press_m(n,k)
        tmp   = lapse*(rdair/gravit)*(ps_old(n)/pbot - 1._r8)
        tsurf = tbot*(1. + tmp)
        t0    = tsurf + lapse*phis_old(n)/gravit

        ! Based on heuristic equations:
        !-------------------------------
        if((t0.gt.t_ref1).and.(tsurf.le.t_ref1)) then
          lapse = (t_ref1 - tsurf)*gravit/phis_old(n)
        elseif((t0.gt.t_ref1).and.(tsurf.gt.t_ref1)) then
          lapse = 0.
          tsurf = (t_ref1 + tsurf)*0.5_r8
        endif

        if(tsurf .lt. t_ref2) then
          lapse = -dtdz
          tsurf = (t_ref2 + tsurf)*0.5_r8
        endif              

        x   = lapse*del_phis/(gravit*tsurf)
        tmp = 1._r8 - x/2._r8 + x**2._r8/3._r8
        tmp = del_phis/(rdair*tsurf)*tmp
        ps_new(n) = ps_old(n)*exp(tmp)
      endif

    end do

    ! End Routine
    !-------------
    return
  end subroutine psadj_se
  !=====================================================================


  !=====================================================================
  subroutine q2rh_se(plev, pcol, q, t, press)
    !
    ! Compute RH from T, Pressure, and Specific Humidity
    !===================================================================
    implicit none
    !
    ! Passed variables
    !------------------
    integer :: plev
    integer :: pcol
    real(r8):: q    (pcol,plev)
    real(r8):: t    (pcol,plev)
    real(r8):: press(pcol,plev)
    !
    !---------------------------- Commons ----------------------------------
    !
    ! Common block and statement functions for saturation vapor pressure
    ! look-up procedure, J. J. Hack, February 1990
    !
!    integer,parameter:: plenest=250 ! length of saturation vapor pressure table
    !
    ! Table of saturation vapor pressure values es from tmin degrees
    ! to tmax+1 degrees k in one degree increments.  ttrice defines the
    ! transition region where es is a combination of ice & water values
    !
!    common/comes/ estbl(plenest), tmin, tmax, ttrice, pcf(6), icephs
!    real(r8):: estbl      ! table values of saturation vapor pressure
!    real(r8):: tmin       ! min temperature (K) for table
!    real(r8):: tmax       ! max temperature (K) for table
!    real(r8):: ttrice     ! transition range from es over H2O to over ice
!    real(r8):: pcf        ! polynomial coeffs -> es transition h2o to ice
!    logical :: icephs     ! false => saturation vapor press over water only
    !
    ! Dummy variables for statement functions
    !-----------------------------------------
    real(r8):: td         ! dummy variable for function evaluation
    real(r8):: tlim       ! intermediate var for es look-up with estbl4
    real(r8):: estblf     ! statement function es look-up
    real(r8):: estbl4     ! statement function es look-up
    !
    ! Local values
    !--------------
    real(r8):: es
    real(r8):: qs
    real(r8):: epsilo 
    real(r8):: latvap 
    real(r8):: latice 
    real(r8):: rh2o   
    real(r8):: cpair  
    real(r8):: omeps              ! 1 - 0.622
    real(r8):: qmin,qmax
    real(r8):: one
    integer :: countneg,countgt1
    integer :: n, k           ! Indices
    integer :: nmin, kmin
    integer :: nmax, kmax
    integer :: itype
    !
    !-------------------------- Statement Functions ------------------------
    !
    ! Statement functions used in saturation vapor pressure table lookup
    ! there are two ways to use these three statement functions.
    ! For compilers that do a simple in-line expansion:
    ! => ttemp = tlim(t)
    !    es    = estbl4(ttemp)
    !
    ! For compilers that provide real optimization:
    ! => es    = estblf(t)
    !
    tlim(td) = max(min(td,COM%tmax),COM%tmin)
    estblf(td) = (COM%tmin + int(tlim(td)-COM%tmin) - tlim(td)+1._r8) &
                  *COM%estbl(int(tlim(td)-COM%tmin)+1)                &
                -(COM%tmin + int(tlim(td)-COM%tmin) - tlim(td)      ) &
                  *COM%estbl(int(tlim(td)-COM%tmin)+2)
    estbl4(td) = (COM%tmin+int(td-COM%tmin)+1._r8-td)*COM%estbl(int(td-COM%tmin)+1) &
                +( td-(COM%tmin+int(td-COM%tmin))   )*COM%estbl(int(td-COM%tmin)+2)
    !-----------------------------------------------------------------------

    one    = 1.0_r8
    epsilo = 0.622_r8
    latvap = 2.5104d06
    latice = 3.336d5
    rh2o   = 4.61d2
    cpair  = 1004.64_r8
    omeps  = 1.0_r8 - epsilo

    ! Build es table
    !-------------------
    call esinti(epsilo, latvap, latice, rh2o, cpair)

    qmin  =  1.d+36
    qmax  = -1.d+36
    countneg = 0
    countgt1 = 0

    do k = 1,plev
    do n = 1,pcol

      ! Saturation specific humidity
      !------------------------------
      es = estblf( t(n,k) )
      qs = epsilo*es/(press(n,k) - omeps*es)

      ! The following check is to avoid the generation of negative values
      ! that can occur in the upper stratosphere and mesosphere
      !-------------------------------------------------------------------
      qs = min(one,qs)
      if(qs.lt.0.0_r8) qs = 1.0_r8

      ! Compute RH
      !-----------
      q(n,k) = q(n,k)/qs

      if(q(n,k) .lt. qmin) then
        qmin = q(n,k)
        nmin = n
        kmin = k
      endif
      if(q(n,k) .gt. qmax) then
        qmax = q(n,k)
        nmax = n
        kmax = k
      endif

      if(q(n,k) .lt. 0.0_r8) then
        countneg = countneg + 1
        q(n,k) = 0._r8
      endif
      if(q(n,k) .gt. 1._r8) then
        countgt1 = countgt1 + 1
        q(n,k) = 1._r8
      endif

    end do
    end do

    write(iulog,*) ' '
    write(iulog,*) ' '
    write(iulog,1000) qmax*100.,nmax,kmax
    write(iulog,2000) qmin*100.,nmin,kmin
    if(countneg.gt.0) write(iulog,3000) countneg, pcol*plev,float(countneg)/(pcol*plev)*100.
    if(countgt1.gt.0) write(iulog,4000) countgt1, pcol*plev,float(countgt1)/(pcol*plev)*100.
    write(iulog,*) ' '
    write(iulog,*) ' '

 1000 format(' Maximum RH = ',f12.5,'% at n,k = ',2i5)
 2000 format(' Minimum RH = ',f12.5,'% at n,k = ',2i5)
 3000 format(' WARNING in Q2RH:  ',i10,' points (out of ',i10, &
             ' (',f12.5,'%) ) were negative.  All were set to 0')
 4000 format(' WARNING in Q2RH:  ',i10,' points (out of ',i10, &
             ' (',f12.5,'%) ) were  .gt. 1.   All were set to 1')

    ! End Routine
    !------------
    return
  end subroutine q2rh_se
  !=====================================================================


  !=====================================================================
  subroutine rh2q_se(plev, pcol, q, t, press)
    !
    ! Compute Q from T, Pressure, and RH
    !===================================================================
    implicit none
    !
    ! Passed variables
    !-------------------
    integer :: plev
    integer :: pcol
    real(r8):: q    (pcol,plev)
    real(r8):: t    (pcol,plev)
    real(r8):: press(pcol,plev)
    !
    !---------------------------- Commons ----------------------------------
    !
    ! Common block and statement functions for saturation vapor pressure
    ! look-up procedure, J. J. Hack, February 1990
    !
!    integer,parameter:: plenest=250 ! length of saturation vapor pressure table
    !
    ! Table of saturation vapor pressure values es from tmin degrees
    ! to tmax+1 degrees k in one degree increments.  ttrice defines the
    ! transition region where es is a combination of ice & water values
    !
!    common/comes/ estbl(plenest), tmin, tmax, ttrice, pcf(6), icephs
!    real(r8):: estbl      ! table values of saturation vapor pressure
!    real(r8):: tmin       ! min temperature (K) for table
!    real(r8):: tmax       ! max temperature (K) for table
!    real(r8):: ttrice     ! transition range from es over H2O to over ice
!    real(r8):: pcf        ! polynomial coeffs -> es transition h2o to ice
!    logical :: icephs     ! false => saturation vapor press over water only
    !
    ! Dummy variables for statement functions
    !-------------------------------------------
    real(r8):: td         ! dummy variable for function evaluation
    real(r8):: tlim       ! intermediate var for es look-up with estbl4
    real(r8):: estblf     ! statement function es look-up
    real(r8):: estbl4     ! statement function es look-up
    !
    ! Local Values 
    !-----------------
    real(r8):: es
    real(r8):: qs
    real(r8):: epsilo 
    real(r8):: latvap 
    real(r8):: latice 
    real(r8):: rh2o   
    real(r8):: cpair  
    real(r8):: omeps              ! 1 - 0.622
    real(r8):: one
    integer :: n, k           ! Indices
    integer :: itype
    !
    !-------------------------- Statement Functions ------------------------
    !
    ! Statement functions used in saturation vapor pressure table lookup
    ! there are two ways to use these three statement functions.
    ! For compilers that do a simple in-line expansion:
    ! => ttemp = tlim(t)
    !    es    = estbl4(ttemp)
    !
    ! For compilers that provide real optimization:
    ! => es    = estblf(t)
    !
    tlim(td) = max(min(td,COM%tmax),COM%tmin)
    estblf(td) = (COM%tmin + int(tlim(td)-COM%tmin) - tlim(td)+1._r8) &
                  *COM%estbl(int(tlim(td)-COM%tmin)+1)                &
                -(COM%tmin + int(tlim(td)-COM%tmin) - tlim(td)      ) &
                  *COM%estbl(int(tlim(td)-COM%tmin)+2)
    estbl4(td) = (COM%tmin+int(td-COM%tmin)+1.0_r8-td)*COM%estbl(int(td-COM%tmin)+1) &
                +( td-(COM%tmin+int(td-COM%tmin))    )*COM%estbl(int(td-COM%tmin)+2)
    !-----------------------------------------------------------------------

    one    = 1.0_r8
    epsilo = 0.622_r8
    latvap = 2.5104d06
    latice = 3.336d5
    rh2o   = 4.61d2
    cpair  = 1004.64_r8
    omeps  = 1.0_r8 - epsilo

    ! Build es table
    !------------------
    call esinti(epsilo, latvap, latice, rh2o, cpair)

    do k = 1,plev
    do n = 1,pcol

      ! Saturation specific humidity
      !------------------------------
      es = estblf( t(n,k) )
      qs = epsilo*es/(press(n,k) - omeps*es)

      ! The following check is to avoid the generation of negative values
      ! that can occur in the upper stratosphere and mesosphere
      !------------------------------------------------------------------
      qs = min(one,qs)
      if(qs.lt.0.0_r8) qs = 1.0_r8

      ! Compute Q from RH and Qs
      !--------------------------
      q(n,k) = q(n,k)*qs

    end do
    end do

    ! End Routine
    !--------------
    return
  end subroutine rh2q_se
  !=====================================================================


  !=====================================================================
  subroutine q2rh_rect(plat, plev, plon, q, t, press)
    !
    ! Compute RH from T, Pressure, and Specific Humidity
    !===================================================================
    implicit none
    !
    ! Passed variables
    !---------------------
    integer :: plat
    integer :: plev
    integer :: plon
    real(r8):: q    (plon,plat,plev)
    real(r8):: t    (plon,plat,plev)
    real(r8):: press(plon,plat,plev)
    !
    !---------------------------- Commons ----------------------------------
    !
    ! Common block and statement functions for saturation vapor pressure
    ! look-up procedure, J. J. Hack, February 1990
    !
!    integer,parameter:: plenest=250 ! length of saturation vapor pressure table
    !
    ! Table of saturation vapor pressure values es from tmin degrees
    ! to tmax+1 degrees k in one degree increments.  ttrice defines the
    ! transition region where es is a combination of ice & water values
    !
!    common/comes/ estbl(plenest), tmin, tmax, ttrice, pcf(6), icephs
!    real(r8):: estbl      ! table values of saturation vapor pressure
!    real(r8):: tmin       ! min temperature (K) for table
!    real(r8):: tmax       ! max temperature (K) for table
!    real(r8):: ttrice     ! transition range from es over H2O to over ice
!    real(r8):: pcf        ! polynomial coeffs -> es transition h2o to ice
!    logical :: icephs     ! false => saturation vapor press over water only
    !
    ! Dummy variables for statement functions
    !-----------------------------------------
    real(r8):: td         ! dummy variable for function evaluation
    real(r8):: tlim       ! intermediate var for es look-up with estbl4
    real(r8):: estblf     ! statement function es look-up
    real(r8):: estbl4     ! statement function es look-up
    !
    ! Local values
    !----------------
    real(r8):: es
    real(r8):: qs
    real(r8):: epsilo 
    real(r8):: latvap 
    real(r8):: latice 
    real(r8):: rh2o   
    real(r8):: cpair  
    real(r8):: omeps              ! 1 - 0.622
    real(r8):: qmin,qmax
    real(r8):: one
    integer :: countneg,countgt1
    integer :: i, j, k           ! Indices
    integer :: imin, jmin, kmin
    integer :: imax, jmax, kmax
    integer :: itype
    !
    !-------------------------- Statement Functions ------------------------
    !
    ! Statement functions used in saturation vapor pressure table lookup
    ! there are two ways to use these three statement functions.
    ! For compilers that do a simple in-line expansion:
    ! => ttemp = tlim(t)
    !    es    = estbl4(ttemp)
    !
    ! For compilers that provide real optimization:
    ! => es    = estblf(t)
    !
    tlim(td) = max(min(td,COM%tmax),COM%tmin)
    estblf(td) = (COM%tmin + int(tlim(td)-COM%tmin) - tlim(td)+1._r8) &
                  *COM%estbl(int(tlim(td)-COM%tmin)+1)                &
                -(COM%tmin + int(tlim(td)-COM%tmin) - tlim(td)      ) &
                  *COM%estbl(int(tlim(td)-COM%tmin)+2)
    estbl4(td) = (COM%tmin+int(td-COM%tmin)+1.0_r8-td)*COM%estbl(int(td-COM%tmin)+1) &
                +( td-(COM%tmin+int(td-COM%tmin))    )*COM%estbl(int(td-COM%tmin)+2)
    !-----------------------------------------------------------------------

    one    = 1.0_r8
    epsilo = 0.622_r8
    latvap = 2.5104d06
    latice = 3.336d5
    rh2o   = 4.61d2
    cpair  = 1004.64_r8
    omeps  = 1.0_r8 - epsilo

    ! Build es table
    !-----------------
    call esinti(epsilo, latvap, latice, rh2o, cpair)

    qmin  =  1.d+36
    qmax  = -1.d+36
    countneg = 0
    countgt1 = 0

    do k = 1,plev
    do j = 1,plat
    do i = 1,plon

      ! Saturation specific humidity
      !-------------------------------
      es = estblf( t(i,j,k) )
      qs = epsilo*es/(press(i,j,k) - omeps*es)

      ! The following check is to avoid the generation of negative values
      ! that can occur in the upper stratosphere and mesosphere
      !----------------------------------------------------------------------
      qs = min(one,qs)
      if(qs.lt.0.0_r8) qs = 1.0_r8

      ! Compute RH
      !-------------
      q(i,j,k) = q(i,j,k)/qs

      if(q(i,j,k) .lt. qmin) then
        qmin = q(i,j,k)
        imin = i
        kmin = k
        jmin = j
      endif
      if(q(i,j,k) .gt. qmax) then
        qmax = q(i,j,k)
        imax = i
        kmax = k
        jmax = j
      endif

      if(q(i,j,k) .lt. 0.0_r8) then
        countneg = countneg + 1
        q(i,j,k) = 0._r8
      endif
      if(q(i,j,k) .gt. 1._r8) then
        countgt1 = countgt1 + 1
        q(i,j,k) = 1._r8
      endif

    end do
    end do
    end do

    write(iulog,*) ' '
    write(iulog,*) ' '
    write(iulog,1000) qmax*100.,imax,jmax,kmax
    write(iulog,2000) qmin*100.,imin,jmin,kmin
    if(countneg.gt.0) write(iulog,3000) countneg,plon*plev*plat, &
                                        float(countneg)/(plon*plev*plat)*100.
    if(countgt1.gt.0) write(iulog,4000) countgt1,plon*plev*plat, &
                                        float(countgt1)/(plon*plev*plat)*100.
    write(iulog,*) ' '
    write(iulog,*) ' '

 1000 format(' Maximum RH = ',f12.5,'% at i,j,k = ',3i5)
 2000 format(' Minimum RH = ',f12.5,'% at i,j,k = ',3i5)
 3000 format(' WARNING in Q2RH:  ',i10,' points (out of ',i10, &
             ' (',f12.5,'%) ) were negative.  All were set to 0')
 4000 format(' WARNING in Q2RH:  ',i10,' points (out of ',i10, &
             ' (',f12.5,'%) ) were  .gt. 1.   All were set to 1')

    ! End Routine
    !------------
    return
  end subroutine q2rh_rect
  !=====================================================================


  !=====================================================================
  subroutine rh2q_rect(plat, plev, plon, q, t, press)
    !
    ! Compute Q from T, Pressure, and RH
    !===================================================================
    implicit none
    !
    ! Passed variables 
    !---------------------
    integer :: plat
    integer :: plev
    integer :: plon
    real(r8):: q    (plon,plat,plev)
    real(r8):: t    (plon,plat,plev)
    real(r8):: press(plon,plat,plev)
    !
    !---------------------------- Commons ----------------------------------
    !
    ! Common block and statement functions for saturation vapor pressure
    ! look-up procedure, J. J. Hack, February 1990
    !
!    integer,parameter:: plenest=250 ! length of saturation vapor pressure table
    !
    ! Table of saturation vapor pressure values es from tmin degrees
    ! to tmax+1 degrees k in one degree increments.  ttrice defines the
    ! transition region where es is a combination of ice & water values
    !
!    common/comes/ estbl(plenest), tmin, tmax, ttrice, pcf(6), icephs
!    real(r8):: estbl      ! table values of saturation vapor pressure
!    real(r8):: tmin       ! min temperature (K) for table
!    real(r8):: tmax       ! max temperature (K) for table
!    real(r8):: ttrice     ! transition range from es over H2O to over ice
!    real(r8):: pcf        ! polynomial coeffs -> es transition h2o to ice
!    logical :: icephs     ! false => saturation vapor press over water only
    !
    ! Dummy variables for statement functions
    !------------------------------------------
    real(r8):: td         ! dummy variable for function evaluation
    real(r8):: tlim       ! intermediate var for es look-up with estbl4
    real(r8):: estblf     ! statement function es look-up
    real(r8):: estbl4     ! statement function es look-up
    !
    ! Local Values
    !-----------------
    real(r8):: es
    real(r8):: qs
    real(r8):: epsilo 
    real(r8):: latvap 
    real(r8):: latice 
    real(r8):: rh2o   
    real(r8):: cpair  
    real(r8):: omeps              ! 1 - 0.622
    real(r8):: one
    integer :: i, j, k           ! Indices
    integer :: itype
    !
    !-------------------------- Statement Functions ------------------------
    !
    ! Statement functions used in saturation vapor pressure table lookup
    ! there are two ways to use these three statement functions.
    ! For compilers that do a simple in-line expansion:
    ! => ttemp = tlim(t)
    !    es    = estbl4(ttemp)
    !
    ! For compilers that provide real optimization:
    ! => es    = estblf(t)
    !
    tlim(td) = max(min(td,COM%tmax),COM%tmin)
    estblf(td) = (COM%tmin + int(tlim(td)-COM%tmin) - tlim(td)+1._r8) &
                  *COM%estbl(int(tlim(td)-COM%tmin)+1)                &
                -(COM%tmin + int(tlim(td)-COM%tmin) - tlim(td)      ) &
                  *COM%estbl(int(tlim(td)-COM%tmin)+2)
    estbl4(td) = (COM%tmin+int(td-COM%tmin)+1.0_r8-td)*COM%estbl(int(td-COM%tmin)+1) &
                +( td-(COM%tmin+int(td-COM%tmin))    )*COM%estbl(int(td-COM%tmin)+2)
    !-----------------------------------------------------------------------

    one    = 1.0_r8
    epsilo = 0.622_r8
    latvap = 2.5104d06
    latice = 3.336d5
    rh2o   = 4.61d2
    cpair  = 1004.64_r8
    omeps  = 1.0_r8 - epsilo

    ! Build es table
    !------------------
    call esinti(epsilo, latvap, latice, rh2o, cpair)

    do k = 1,plev
    do j = 1,plat
    do i = 1,plon

      ! Saturation specific humidity
      !---------------------------------
      es = estblf( t(i,j,k) )
      qs = epsilo*es/(press(i,j,k) - omeps*es)

      ! The following check is to avoid the generation of negative values
      ! that can occur in the upper stratosphere and mesosphere
      !---------------------------------------------------------------------
      qs = min(one,qs)
      if(qs.lt.0.0_r8) qs = 1.0_r8

      ! Compute Q from RH and Qs
      !----------------------------
      q(i,j,k) = q(i,j,k)*qs

    end do
    end do
    end do

    ! End Routine
    !--------------
    return
  end subroutine rh2q_rect
  !=====================================================================


  !=====================================================================
  subroutine esinti(epslon, latvap, latice, rh2o, cpair)
    !
    ! Initialize es lookup tables 
    !===================================================================
    implicit none
    !
    ! Passed Variables
    !-------------------
    real(r8):: epslon        ! Ratio of h2o to dry air molecular weights 
    real(r8):: latvap        ! Latent heat of vaporization
    real(r8):: latice        ! Latent heat of fusion
    real(r8):: rh2o          ! Gas constant for water vapor
    real(r8):: cpair         ! Specific heat of dry air
    !
    ! Local values
    !-----------------
    real(r8):: tmn           ! Minimum temperature entry in table
    real(r8):: tmx           ! Maximum temperature entry in table
    real(r8):: trice         ! Trans range from es over h2o to es over ice
    logical :: ip           ! Ice phase (true or false)

    ! Specify control parameters first
    !------------------------------------
    tmn   = 173.16_r8
    tmx   = 375.16_r8
    trice =  20.00_r8
    ip    = .true.

    ! Call gestbl to build saturation vapor pressure table.
    !---------------------------------------------------------
    call gestbl(tmn, tmx, trice, ip, epslon, latvap, latice, rh2o, cpair)

    ! End Routine
    !--------------
    return
  end subroutine esinti
  !=====================================================================


  !=====================================================================
  subroutine gestbl(tmn, tmx, trice, ip, epsil, latvap, latice, rh2o, cpair)
    !
    ! Builds saturation vapor pressure table for later lookup procedure.
    ! Uses Goff & Gratch (1946) relationships to generate the table
    ! according to a set of free parameters defined below.  Auxiliary
    ! routines are also included for making rapid estimates (well with 1%)
    ! of both es and d(es)/dt for the particular table configuration.
    !
    !===================================================================
    implicit none
    !
    ! Passed varaibles
    !--------------------
    real(r8):: tmn        ! Minimum temperature entry in es lookup table
    real(r8):: tmx        ! Maximum temperature entry in es lookup table
    real(r8):: epsil      ! Ratio of h2o to dry air molecular weights
    real(r8):: trice      ! Transition range from es over range to es over ice
    real(r8):: latvap     ! Latent heat of vaporization
    real(r8):: latice     ! Latent heat of fusion
    real(r8):: rh2o       ! Gas constant for water vapor
    real(r8):: cpair      ! Specific heat of dry air
    !
    ! Local values
    !---------------
    real(r8):: t          ! Temperature
    integer :: n          ! Increment counter
    integer :: lentbl     ! Calculated length of lookup table
    integer :: itype      ! Ice phase: 0 -> no ice phase
                          !            1 -> ice phase, no transition
                          !           -x -> ice phase, x degree transitn
    logical :: ip         ! Ice phase logical flag
    !
    !---------------------------- Commons ----------------------------------
    !
    ! Common block and statement functions for saturation vapor pressure
    ! look-up procedure, J. J. Hack, February 1990
    !
!    integer,parameter:: plenest=250 ! length of saturation vapor pressure table
    !
    ! Table of saturation vapor pressure values es from tmin degrees
    ! to tmax+1 degrees k in one degree increments.  ttrice defines the
    ! transition region where es is a combination of ice & water values
    !
!    common/comes/ estbl(plenest), tmin, tmax, ttrice, pcf(6), icephs
!    real(r8):: estbl      ! table values of saturation vapor pressure
!    real(r8):: tmin       ! min temperature (K) for table
!    real(r8):: tmax       ! max temperature (K) for table
!    real(r8):: ttrice     ! transition range from es over H2O to over ice
!    real(r8):: pcf        ! polynomial coeffs -> es transition h2o to ice
!    logical :: icephs     ! false => saturation vapor press over water only
    !-----------------------------------------------------------------------

    ! Set es table parameters
    !----------------------------
    COM%tmin   = tmn       ! Minimum temperature entry in table
    COM%tmax   = tmx       ! Maximum temperature entry in table
    COM%ttrice = trice     ! Trans. range from es over h2o to es over ice
    COM%icephs = ip        ! Ice phase (true or false)

    ! Set physical constants required for es calculation
    !--------------------------------------------------
    lentbl = COM%tmax-COM%tmin+2.000001_r8
    if(lentbl .gt. plenest) then
      write(iulog,9000) COM%tmax, COM%tmin, plenest
      call endrun     ! Abnormal termination
    endif

    ! Begin building es table.
    ! Check whether ice phase requested.
    ! If so, set appropriate transition range for temperature
    !----------------------------------------------------------
    if(COM%icephs) then
      if(COM%ttrice.ne.0.0_r8) then
        itype = -COM%ttrice
      else
        itype = 1
      endif
    else
      itype = 0
    endif

    t = COM%tmin - 1._r8
    do n=1,lentbl
      t = t + 1._r8
      call gffgch(t,COM%estbl(n),itype)
    end do

    do n=lentbl+1,plenest
      COM%estbl(n) = -99999.0_r8
    end do

    ! Table complete -- Set coefficients for polynomial approximation of
    ! difference between saturation vapor press over water and saturation
    ! pressure over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial
    ! is valid in the range -40 < t < 0 (degrees C).
    !
    !  Degree 5 approximation
    !---------------------------------
    COM%pcf(1) =  5.04469588506d-01
    COM%pcf(2) = -5.47288442819d+00
    COM%pcf(3) = -3.67471858735d-01
    COM%pcf(4) = -8.95963532403d-03
    COM%pcf(5) = -7.78053686625d-05
    !
    !  Degree 6 approximation ---
    !----------------------------------
!   COM%pcf(1) =  7.63285250063d-02
!   COM%pcf(2) = -5.86048427932d+00
!   COM%pcf(3) = -4.38660831780d-01
!   COM%pcf(4) = -1.37898276415d-02
!   COM%pcf(5) = -2.14444472424d-04
!   COM%pcf(6) = -1.36639103771d-06

 9000 format('GESTBL: FATAL ERROR *********************************',/, &
           ' TMAX AND TMIN REQUIRE A LARGER DIMENSION ON THE LENGTH', &
           ' OF THE SATURATION VAPOR PRESSURE TABLE ESTBL(PLENEST)',/, &
           ' TMAX, TMIN, AND PLENEST => ', 2f7.2, i3)

    ! End Routine
    !---------------
    return
  end subroutine gestbl
  !=====================================================================


  !=====================================================================
  subroutine gffgch(t, es, itype)
    !
    ! Computes saturation vapor pressure over water and/or over ice using
    ! Goff & Gratch (1946) relationships.  T (temperature), and itype are
    ! input parameters, while es (saturation vapor pressure) is an output
    ! parameter.  The input parameter itype serves two purposes: a value of
    ! zero indicates that saturation vapor pressures over water are to be
    ! returned (regardless of temperature), while a value of one indicates
    ! that saturation vapor pressures over ice should be returned when t is
    ! less than 273.16 degrees k.  If itype is negative, its absolute value
    ! is interpreted to define a temperature transition region below 273.16
    ! degrees k in which the returned saturation vapor pressure is a
    ! weighted average of the respective ice and water value.  That is, in
    ! the temperature range 0 => -itype degrees c, the saturation vapor
    ! pressures are assumed to be a weighted average of the vapor pressure
    ! over supercooled water and ice (all water at 0 c; all ice at -itype
    ! c).  Maximum transition range => 40 c
    !
    !===================================================================
    implicit none
    !
    ! Passed variables
    !-------------------
    real(r8):: t         ! Temperature
    integer :: itype     ! Flag for ice phase and associated transition
    real(r8):: es        ! Saturation vapor pressure
    !
    ! Local Values 
    !------------------
    real(r8):: e1        ! Intermediate scratch variable for es over water
    real(r8):: e2        ! Intermediate scratch variable for es over water
    real(r8):: eswtr     ! Saturation vapor pressure over water
    real(r8):: f         ! Intermediate scratch variable for es over water
    real(r8):: f1        ! Intermediate scratch variable for es over water
    real(r8):: f2        ! Intermediate scratch variable for es over water
    real(r8):: f3        ! Intermediate scratch variable for es over water
    real(r8):: f4        ! Intermediate scratch variable for es over water
    real(r8):: f5        ! Intermediate scratch variable for es over water
    real(r8):: ps        ! Reference pressure (mb)
    real(r8):: t0        ! Reference temperature (freezing point of water)
    real(r8):: term1     ! Intermediate scratch variable for es over ice
    real(r8):: term2     ! Intermediate scratch variable for es over ice
    real(r8):: term3     ! Intermediate scratch variable for es over ice
    real(r8):: tr        ! Transition range for es over liq to es over ice
    real(r8):: ts        ! Reference temperature (boiling point of water)
    real(r8):: weight    ! Intermediate scratch variable for es transition
    real(r8):: tfreez, one
    integer :: itypo    ! Intermediate scratch variable for holding itype

    tfreez = 273.16_r8
    one    = 1.0_r8

    ! Check on whether there is to be a transition region for es
    !------------------------------------------------------------
    if(itype.lt.0) then
      tr    = abs(float(itype))
      itypo = itype
      itype = 1
    else
      tr    = 0.0_r8
      itypo = itype
    endif
    if(tr .gt. 40.0_r8) then
      write(iulog,900) tr
      call endrun                 ! Abnormal termination
    endif

    if((t.lt.(tfreez - tr)).and.(itype.eq.1)) go to 10
    ! Water
    !--------
    ps = 1013.246_r8
    ts = 373.16_r8
    e1 =   11.344_r8*(1.0_r8 - t/ts)
    e2 = -3.49149_r8*(ts/t - 1.0_r8)
    f1 = -7.90298_r8*(ts/t - 1.0_r8)
    f2 = 5.02808_r8*log10(ts/t)
    f3 = -1.3816_r8*(10.0_r8**e1 - 1.0_r8)/10000000.0_r8
    f4 =  8.1328_r8*(10.0_r8**e2 - 1.0_r8)/1000.0_r8
    f5 = log10(ps)
    f  = f1 + f2 + f3 + f4 + f5
    es = (10.0_r8**f)*100.0_r8
    eswtr = es

    if(t.ge.tfreez .or. itype.eq.0) go to 20
 10 continue
    ! Ice
    !-----------
    t0    = tfreez
    term1 = 2.01889049_r8/(t0/t)
    term2 = 3.56654_r8*log(t0/t)
    term3 = 20.947031_r8*(t0/t)
    es    = 575.185606d10*exp(-(term1 + term2 + term3))

    if(t.lt.(tfreez-tr)) go to 20
      ! Weighted transition between water and ice
      !------------------------------------------
      weight = min((tfreez - t)/tr,one)
      es = weight*es + (1.0_r8 - weight)*eswtr

 20 continue
    itype = itypo

  900 format('GFFGCH: FATAL ERROR ******************************',/, &
             'TRANSITION RANGE FOR WATER TO ICE SATURATION VAPOR', &
             ' PRESSURE, TR, EXCEEDS MAXIMUM ALLOWABLE VALUE OF', &
             ' 40.0 DEGREES C',/, ' TR = ',f7.2)

    ! End Routine
    !-------------
    return
  end subroutine gffgch
  !=====================================================================
 

  !=====================================================================
  subroutine mass_fixer_rect(plev, plevp1, plat, plon, q,           &
                             hyai, hybi, gw, gravit, ps0, tmass0, ps)
    !
    ! Adjust atmospheric mass based upon Q.
    !===================================================================
    implicit none
    !
    ! Passed variables 
    !--------------------
    integer :: plev                   !  level dimension
    integer :: plevp1                 !  level dimension + 1
    integer :: plat                   !  latitude dimension
    integer :: plon                   !  longitude dimension
    real(r8):: hyai(plevp1)
    real(r8):: hybi(plevp1)
    real(r8):: gw  (plat)
    real(r8):: q   (plon,plat,plev)   ! Specific humidity
    real(r8):: gravit                 ! acceleration of gravity ~ m/s^2
    real(r8):: ps0                    ! Ref. Surface pressure (10**5 Pa)
    real(r8):: tmass0                 ! Dry mass of Ref. atmosphere
    real(r8):: ps(plon,plat)          ! Surface Pressure
    !
    ! Local Values
    !-----------------
    real(r8):: tmassf               ! unfixed mass of atmosphere
    real(r8):: qmass1               ! "a" contribution to mass of q
    real(r8):: qmass2               ! "b" contribution to mass of q
    real(r8):: pdela(plon,plev)     ! "a" contribution to del-P
    real(r8):: pdelb(plon,plev)     ! "b" contribution to del-P
    real(r8):: pssum                ! !
    real(r8):: dotproda             ! ! - accumulators
    real(r8):: dotprodb             ! !
    real(r8):: fixmas               ! Mass fix coefficient
    integer :: i, lat, k            ! Indices

    tmassf = 0._r8
    qmass1 = 0._r8
    qmass2 = 0._r8

    ! Compute pdel from "A" portion of hybrid vertical 
    ! grid for later use in global integrals
    !-----------------------------------------------------------
    do k = 1,plev
    do i = 1,plon
      pdela(i,k) = (hyai(k+1) - hyai(k))*ps0
    end do
    end do

    ! Compute integrals of mass, moisture, and geopotential height
    !-------------------------------------------------------------
    do lat = 1,plat
      ! Accumulate average mass of atmosphere
      !---------------------------------------
      do k = 1,plev
      do i = 1,plon
        pdelb(i,k) = (hybi(k+1) - hybi(k))*ps(i,lat)
      end do
      end do

      pssum  = 0._r8
      do i = 1,plon
        pssum  = pssum  + ps(i,lat)
      end do
      tmassf = tmassf + gw(lat)*pssum/plon

      ! Calculate global integrals needed for water vapor adjustment
      !--------------------------------------------------------------
      do k = 1,plev
        dotproda = 0._r8
        dotprodb = 0._r8
        do i = 1,plon
          dotproda = dotproda + q(i,lat,k)*pdela(i,k)
          dotprodb = dotprodb + q(i,lat,k)*pdelb(i,k)
        end do
        qmass1 = qmass1 + gw(lat)*dotproda/plon
        qmass2 = qmass2 + gw(lat)*dotprodb/plon
      end do
    end do ! lat = 1,plat

    ! Normalize average mass, height
    !---------------------------------
    tmassf = tmassf*.5_r8/gravit
    qmass1 = qmass1*.5_r8/gravit
    qmass2 = qmass2*.5_r8/gravit

    ! Compute and apply an initial mass fix factor which 
    ! preserves horizontal gradients of ln(ps).
    !---------------------------------------------------------
    fixmas = (tmass0 + qmass1)/(tmassf - qmass2)

    do lat = 1,plat
    do i = 1,plon
      ps(i,lat) = ps(i,lat)*fixmas
    end do
    end do

    write(iulog,1000) fixmas
1000  format("             FIXMAS = ", f18.14/)

    ! End Routine
    !--------------
    return
  end subroutine mass_fixer_rect
  !=====================================================================


  !=====================================================================
  subroutine mass_fixer_se(plev, plevp1, pcol, q,                   &
                           hyai, hybi, area, gravit, ps0, tmass0, ps)
    !
    ! Adjust atmospheric mass based upon Q.
    !===================================================================
    implicit none
    !
    ! Passed variables 
    !--------------------
    integer :: plev              !  level dimension
    integer :: plevp1            !  level dimension + 1
    integer :: pcol
    real(r8):: hyai(plevp1)
    real(r8):: hybi(plevp1)
    real(r8):: area(pcol)
    real(r8):: q   (pcol,plev)   ! Specific humidity
    real(r8):: gravit            ! acceleration of gravity ~ m/s^2
    real(r8):: ps0               ! Ref. Surface pressure (10**5 Pa)
    real(r8):: tmass0            ! Dry mass of Ref. atmosphere
    real(r8):: ps(pcol)          ! Surface Pressure
    !
    ! Local values
    !---------------
    real(r8):: tmassf         ! unfixed mass of atmosphere
    real(r8):: qmass1         ! "a" contribution to mass of q
    real(r8):: qmass2         ! "b" contribution to mass of q
    real(r8):: fixmas         ! Mass fix coefficient
    real(r8):: anorm          
    integer :: n, k          ! Indices


    ! Compute integrals of mass, moisture, and geopotential height
    ! Accumulate average mass of atmosphere
    ! Calculate global integrals needed for water vapor adjustment
    !---------------------------------------------------------------
    anorm=0._r8
    do n = 1,pcol
      anorm = anorm + area(n)
    end do

    tmassf = 0._r8
    qmass1 = 0._r8
    qmass2 = 0._r8
    do n = 1,pcol
      tmassf = tmassf + (area(n)*ps(n))
      do k = 1,plev
        qmass1 = qmass1 + (area(n)*q(n,k)*((hyai(k+1)-hyai(k))*ps0  ))
        qmass2 = qmass2 + (area(n)*q(n,k)*((hybi(k+1)-hybi(k))*ps(n)))
      end do
    end do ! n = 1,pcol

    ! Normalize average mass, height
    !--------------------------------
    tmassf = tmassf/(gravit*anorm)
    qmass1 = qmass1/(gravit*anorm)
    qmass2 = qmass2/(gravit*anorm)

    ! Compute and apply an initial mass fix factor which 
    ! preserves horizontal gradients of ln(ps).
    !-----------------------------------------------------
    fixmas = (tmass0+qmass1)/(tmassf-qmass2)
    do n = 1,pcol
      ps(n) = ps(n)*fixmas
    end do

    write(iulog,1000) fixmas
1000  format("             FIXMAS = ", f18.14/)

    ! End Routine
    !------------
    return
  end subroutine mass_fixer_se
  !=====================================================================

end module CAPT_mod
