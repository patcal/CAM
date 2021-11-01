module nudging
!=====================================================================
!
! Purpose: Implement Nudging of the model state of U,V,T,Q, and/or PS
!          toward specified values from analyses.
!
! Author: Patrick Callaghan
!
! Description:
!
!    This module assumes that the user has {U,V,T,Q,PS} values from analyses
!    which have been preprocessed onto the current model grid and adjusted
!    for differences in topography. It is also assumed that these resulting
!    values and are stored in individual files which are indexed with respect
!    to year, month, day, and second of the day. When the model is inbetween
!    the given begining and ending times, a relaxation forcing is added to
!    nudge the model toward the analyses values determined from the forcing
!    option specified. After the model passes the ending analyses time, the
!    forcing discontinues.
!
!    Some analyses products can have gaps in the available data, where values
!    are missing for some interval of time. When files are missing, the nudging
!    force is switched off for that interval of time, so we effectively 'coast'
!    thru the gap.
!
!    Currently, the nudging module is set up to accomodate nudging of PS
!    values, however that functionality requires forcing that is applied in
!    the selected dycore and is not yet implemented.
!
!    The nudging of the model toward the analyses data is controlled by
!    the 'nudging_nl' namelist in 'user_nl_cam'; whose variables control the
!    time interval over which nudging is applied, the strength of the nudging
!    tendencies, and its spatial distribution.
!
!    FORCING:
!    --------
!    Nudging tendencies are applied as a relaxation force between the current
!    model state values and target state values derived from the avalilable
!    analyses. The form of the target values is selected by the 'Nudge_Force_Opt'
!    option, the timescale of the forcing is determined from the given
!    'Nudge_TimeScale_Opt', and the nudging strength Alpha=[0.,1.] for each
!    variable is specified by the 'Nudge_Xcoef' values. Where X={U,V,T,Q,PS}
!
!           F_nudge = Alpha*((Target-Model(t_curr))/TimeScale
!
!
!    WINDOWING:
!    ----------
!    The region of applied nudging can be limited using Horizontal/Vertical
!    window functions that are constructed using a parameterization of the
!    Heaviside step function.
!
!    The Heaviside window function is the product of separate horizonal and vertical
!    windows that are controled via 12 parameters:
!
!        Nudge_Hwin_lat0:     Specify the horizontal center of the window in degrees.
!        Nudge_Hwin_lon0:     The longitude must be in the range [0,360] and the
!                             latitude should be [-90,+90].
!        Nudge_Hwin_latWidth: Specify the lat and lon widths of the window as positive
!        Nudge_Hwin_lonWidth: values in degrees.Setting a width to a large value (e.g. 999)
!                             renders the window a constant in that direction.
!        Nudge_Hwin_latDelta: Controls the sharpness of the window transition with a
!        Nudge_Hwin_lonDelta: length in degrees. Small non-zero values yeild a step
!                             function while a large value yeilds a smoother transition.
!        Nudge_Hwin_Invert  : A logical flag used to invert the horizontal window function
!                             to get its compliment.(e.g. to nudge outside a given window).
!
!        Nudge_Vwin_Lindex:   In the vertical, the window is specified in terms of model
!        Nudge_Vwin_Ldelta:   level indcies. The High and Low transition levels should
!        Nudge_Vwin_Hindex:   range from [0,(NLEV+1)]. The transition lengths are also
!        Nudge_Vwin_Hdelta:   specified in terms of model indices. For a window function
!                             constant in the vertical, the Low index should be set to 0,
!                             the High index should be set to (NLEV+1), and the transition
!                             lengths should be set to 0.001
!        Nudge_Vwin_Invert  : A logical flag used to invert the vertical window function
!                             to get its compliment.
!
!        EXAMPLE: For a channel window function centered at the equator and independent
!                 of the vertical (30 levels):
!                        Nudge_Hwin_lat0     = 0.         Nudge_Vwin_Lindex = 0.
!                        Nudge_Hwin_latWidth = 30.        Nudge_Vwin_Ldelta = 0.001
!                        Nudge_Hwin_latDelta = 5.0        Nudge_Vwin_Hindex = 31.
!                        Nudge_Hwin_lon0     = 180.       Nudge_Vwin_Hdelta = 0.001
!                        Nudge_Hwin_lonWidth = 999.       Nudge_Vwin_Invert = .false.
!                        Nudge_Hwin_lonDelta = 1.0
!                        Nudge_Hwin_Invert   = .false.
!
!                 If on the other hand one wanted to apply nudging at the poles and
!                 not at the equator, the settings would be similar but with:
!                        Nudge_Hwin_Invert = .true.
!
!    A user can preview the window resulting from a given set of namelist values before
!    running the model. Lookat_NudgeWindow.ncl is a script avalable in the tools directory
!    which will read in the values for a given namelist and display the resulting window.
!
!    The module is currently configured for only 1 window function. It can readily be
!    extended for multiple windows if the need arises.
!
!
! Input/Output Values:
!    Forcing contributions are available for history file output by
!    the names:    {'Nudge_U','Nudge_V','Nudge_T',and 'Nudge_Q'}
!    The target values that the model state is nudged toward are available for history
!    file output via the variables:  {'Target_U','Target_V','Target_T',and 'Target_Q'}
!
!    &nudging_nl
!      Nudge_Model         - LOGICAL toggle to activate nudging.
!                              TRUE  -> Nudging is on.
!                              FALSE -> Nudging is off.                            [DEFAULT]
!
!      Nudge_Path          - CHAR path to the analyses files.
!                              (e.g. '/glade/scratch/USER/inputdata/nudging/ERAI-Data/')
!
!      Nudge_File_Template - CHAR Analyses filename with year, month, day, and second
!                                 values replaced by %y, %m, %d, and %s respectively.
!                              (e.g. '%y/ERAI_ne30np4_L30.cam2.i.%y-%m-%d-%s.nc')
!
!      Nudge_Times_Per_Day - INT Number of analyses files available per day.
!                              1 --> daily analyses.
!                              4 --> 6 hourly analyses.
!                              8 --> 3 hourly.
!
!      Model_Times_Per_Day - INT Number of times to update the model state (used for nudging)
!                                each day. The value is restricted to be longer than the
!                                current model timestep and shorter than the analyses
!                                timestep. As this number is increased, the nudging
!                                force has the form of newtonian cooling.
!                              48 --> 1800 Second timestep.
!                              96 -->  900 Second timestep.
!
!      Nudge_Beg_Year      - INT nudging begining year.  [1979- ]
!      Nudge_Beg_Month     - INT nudging begining month. [1-12]
!      Nudge_Beg_Day       - INT nudging begining day.   [1-31]
!      Nudge_End_Year      - INT nudging ending year.    [1979-]
!      Nudge_End_Month     - INT nudging ending month.   [1-12]
!      Nudge_End_Day       - INT nudging ending day.     [1-31]
!
!      Nudge_Force_Opt     - INT Index to select the nudging Target for a relaxation
!                                forcing of the form:
!                                where (t'==Analysis times ; t==Model Times)
!
!                              0 -> NEXT-OBS: Target=Anal(t'_next)                 [DEFAULT]
!                              1 -> LINEAR:   Target=(F*Anal(t'_curr) +(1-F)*Anal(t'_next))
!                                                 F =(t'_next - t_curr )/Tdlt_Anal
!
!      Nudge_TimeScale_Opt - INT Index to select the timescale for nudging.
!                                where (t'==Analysis times ; t==Model Times)
!
!                              0 -->  TimeScale = 1/Tdlt_Anal                      [DEFAULT]
!                              1 -->  TimeScale = 1/(t'_next - t_curr )
!
!      Nudge_Uprof         - INT index of profile structure to use for U.  [0,1,2]
!      Nudge_Vprof         - INT index of profile structure to use for V.  [0,1,2]
!      Nudge_Tprof         - INT index of profile structure to use for T.  [0,1,2]
!      Nudge_Qprof         - INT index of profile structure to use for Q.  [0,1,2]
!      Nudge_PSprof        - INT index of profile structure to use for PS. [0,N/A]
!
!                                The spatial distribution is specified with a profile index.
!                                 Where:  0 == OFF      (No Nudging of this variable)
!                                         1 == CONSTANT (Spatially Uniform Nudging)
!                                         2 == HEAVISIDE WINDOW FUNCTION
!
!      Nudge_Ucoef         - REAL fractional nudging coeffcient for U.
!      Nudge_Vcoef         - REAL fractional nudging coeffcient for V.
!      Nudge_Tcoef         - REAL fractional nudging coeffcient for T.
!      Nudge_Qcoef         - REAL fractional nudging coeffcient for Q.
!      Nudge_PScoef        - REAL fractional nudging coeffcient for PS.
!
!                                 The strength of the nudging is specified as a fractional
!                                 coeffcient between [0,1].
!
!      Nudge_Hwin_lat0     - REAL latitudinal center of window in degrees.
!      Nudge_Hwin_lon0     - REAL longitudinal center of window in degrees.
!      Nudge_Hwin_latWidth - REAL latitudinal width of window in degrees.
!      Nudge_Hwin_lonWidth - REAL longitudinal width of window in degrees.
!      Nudge_Hwin_latDelta - REAL latitudinal transition length of window in degrees.
!      Nudge_Hwin_lonDelta - REAL longitudinal transition length of window in degrees.
!      Nudge_Hwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
!                                    TRUE = value=0 inside the specified window, 1 outside
!      Nudge_Vwin_Lindex   - REAL LO model index of transition
!      Nudge_Vwin_Hindex   - REAL HI model index of transition
!      Nudge_Vwin_Ldelta   - REAL LO transition length
!      Nudge_Vwin_Hdelta   - REAL HI transition length
!      Nudge_Vwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
!                                    TRUE = value=0 inside the specified window, 1 outside
!    /
!
!================
!
! TO DO:
! -----------
!    ** Implement Ps Nudging????
!
!=====================================================================
  ! Useful modules
  !------------------
  use ppgrid,         only: pver,pcols,begchunk,endchunk
  use shr_kind_mod,   only: r8=>SHR_KIND_R8, cs=>SHR_KIND_CS, cl=>SHR_KIND_CL
  use time_manager,   only: timemgr_time_ge, timemgr_time_inc, get_curr_date
  use time_manager,   only: get_step_size
  use cam_abortutils, only: endrun
  use cam_pio_utils,  only: cam_pio_openfile,cam_pio_closefile,cam_pio_find_var
  use pio,            only: file_desc_t, var_desc_t, pio_inq_varid, PIO_NOERR,PIO_NOWRITE
  use spmd_utils,     only: masterproc, mstrid=>masterprocid, mpicom
  use spmd_utils,     only: mpi_integer, mpi_real8, mpi_logical, mpi_character
  use cam_logfile,    only: iulog
  use Dataset_mod,    only: open_dataset, close_dataset, print_dataset_metadata
  use Hmaps_mod,      only: init_Hmaps_mod, print_known_Hmaps
  use Vmaps_mod,      only: init_Vmaps_mod, print_known_Vmaps
  use ncdio_atm,      only: register_infld_grid
  use phys_grid,      only: get_rlat_p,get_rlon_p,get_ncols_p
  use netcdf
  use OBSdata_mod,only: OBS_Data_t, INV_Data_t, DataFile_Template_t
  use OBSdata_mod,only: OBS_VarConfig_t, INV_VarConfig_t
  use shr_sys_mod,    only: shr_sys_flush      ! Standardized system subroutines

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public :: Nudge_Model,Nudge_ON
  public :: nudging_readnl
  public :: nudging_init
  public :: nudging_timestep_init
  public :: nudging_timestep_tend
  private:: set_horiz_profile
  private:: set_nudging_profile
  private:: calc_DryStaticEnergy
  private:: initialize_boundary_tau
  private:: interp_Bwin
  private:: lookup_DataType
  private:: replace_string
  private:: days_per_month

  ! Types
  !---------
  type NudgeVar_t
    character(len=8 )   :: Name
    character(len=32)   :: UseOpt
    logical             :: HORIZ_ONLY
    integer             :: Prof
    real(r8)            :: Coef
    real(r8),allocatable:: tau   (:,:,:)  !(pcols,pver,begchunk:endchunk)
    real(r8),allocatable:: target(:,:,:)  !(pcols,pver,begchunk:endchunk)
    real(r8),allocatable:: tend  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  end type NudgeVar_t

  type Nudge_t
    logical          :: Initialized =.false.
    integer          :: Force_Opt
    integer          :: TimeScale_Opt
    integer          :: Times_Per_Day
    integer          :: Beg_Year,Beg_Month
    integer          :: Beg_Day ,Beg_Sec
    integer          :: End_Year,End_Month
    integer          :: End_Day ,End_Sec
    real(r8)         :: Hwin_lat0
    real(r8)         :: Hwin_latWidth
    real(r8)         :: Hwin_latDelta
    real(r8)         :: Hwin_lon0
    real(r8)         :: Hwin_lonWidth
    real(r8)         :: Hwin_lonDelta
    logical          :: Hwin_Invert = .false.
    real(r8)         :: Hwin_lo
    real(r8)         :: Hwin_hi
    real(r8)         :: Vwin_Hindex
    real(r8)         :: Vwin_Hdelta
    real(r8)         :: Vwin_Lindex
    real(r8)         :: Vwin_Ldelta
    logical          :: Vwin_Invert =.false.
    real(r8)         :: Vwin_lo
    real(r8)         :: Vwin_hi
    real(r8)         :: Hwin_latWidthH
    real(r8)         :: Hwin_lonWidthH
    real(r8)         :: Hwin_max
    real(r8)         :: Hwin_min
    integer          :: Step
    logical                     :: use_BDY_file
    logical                     :: BDY_Invert
    character(len=cl)           :: BDY_File
    real(r8)        ,allocatable:: BDY_tau(:,:)    !(pcols,begchunk:endchunk)
    integer                     :: NumVar
    type(NudgeVar_t),allocatable:: var(:)
  end type Nudge_t

  type ModelState_t
    integer             :: Times_Per_Day
    integer             :: Curr_Year,Curr_Month
    integer             :: Curr_Day ,Curr_Sec
    integer             :: Next_Year,Next_Month
    integer             :: Next_Day ,Next_Sec
    integer             :: Step
    real(r8),allocatable:: U (:,:,:)  !(pcols,pver,begchunk:endchunk)
    real(r8),allocatable:: V (:,:,:)  !(pcols,pver,begchunk:endchunk)
    real(r8),allocatable:: T (:,:,:)  !(pcols,pver,begchunk:endchunk)
    real(r8),allocatable:: S (:,:,:)  !(pcols,pver,begchunk:endchunk)
    real(r8),allocatable:: Q (:,:,:)  !(pcols,pver,begchunk:endchunk)
    real(r8),allocatable:: PS(:,:)    !(pcols,begchunk:endchunk)
  end type ModelState_t

  ! Nudging Parameters
  !--------------------
  logical           :: Nudge_Model =.false.
  logical           :: Nudge_ON    =.false.
  type(Nudge_t)     :: Nudge
  type(ModelState_t):: Model

  ! Global OBS/INV data structures
  !-------------------------------
  type(DataFile_Template_t):: INV_FileTemplates
  type(DataFile_Template_t):: OBS_FileTemplates
  type(INV_VarConfig_t)    :: INV_VarConfig
  type(OBS_VarConfig_t)    :: OBS_VarConfig
  type(OBS_Data_t)         :: OBS
  type(INV_Data_t)         :: INV
  integer                  :: idxU,idxV,idxT,idxQ,idxPS,idxT2
  integer                  :: idxPHIS

contains
  !================================================================
  subroutine nudging_readnl(nlfile)
    !
    ! NUDGING_READNL: Initialize default values controlling the Nudging
    !                 process. Then read namelist values to override
    !                 them.
    !===============================================================
    use ppgrid,        only: pver
    use namelist_utils,only:find_group_name
    !
    ! Arguments
    !-------------
    character(len=*),intent(in):: nlfile
    !
    ! Local Values
    !---------------
    real(r8):: lonp,lon0,lonn,latp,lat0,latn
    real(r8):: Val1_p,Val2_p,Val3_p,Val4_p
    real(r8):: Val1_0,Val2_0,Val3_0,Val4_0
    real(r8):: Val1_n,Val2_n,Val3_n,Val4_n
    real(r8):: Nudge_Ucoef,Nudge_Vcoef
    integer :: Nudge_Uprof,Nudge_Vprof
    real(r8):: Nudge_Qcoef,Nudge_Tcoef
    integer :: Nudge_Qprof,Nudge_Tprof
    real(r8):: Nudge_PScoef
    integer :: Nudge_PSprof
    logical :: Nudge_Initialized
    integer :: Nudge_Force_Opt
    integer :: Nudge_TimeScale_Opt
    integer :: Nudge_Times_Per_Day
    integer :: Model_Times_Per_Day
    integer :: Nudge_Beg_Year ,Nudge_Beg_Month
    integer :: Nudge_Beg_Day  ,Nudge_Beg_Sec
    integer :: Nudge_End_Year ,Nudge_End_Month
    integer :: Nudge_End_Day  ,Nudge_End_Sec
    real(r8):: Nudge_Hwin_lat0
    real(r8):: Nudge_Hwin_latWidth
    real(r8):: Nudge_Hwin_latDelta
    real(r8):: Nudge_Hwin_lon0
    real(r8):: Nudge_Hwin_lonWidth
    real(r8):: Nudge_Hwin_lonDelta
    logical :: Nudge_Hwin_Invert
    real(r8):: Nudge_Hwin_lo
    real(r8):: Nudge_Hwin_hi
    real(r8):: Nudge_Vwin_Hindex
    real(r8):: Nudge_Vwin_Hdelta
    real(r8):: Nudge_Vwin_Lindex
    real(r8):: Nudge_Vwin_Ldelta
    logical :: Nudge_Vwin_Invert
    real(r8):: Nudge_Vwin_lo
    real(r8):: Nudge_Vwin_hi
    logical :: Nudge_Use_BDY_File
    logical :: Nudge_BDY_Invert
    character(len=cl):: Nudge_BDY_File
    integer :: ierr, unitn

  character(len=cl):: Nudge_Path
  character(len=cs):: Nudge_File,Nudge_File_Template

    namelist /nudging_nl/ Nudge_Model, Nudge_Path,                        &
                          Nudge_File_Template, Nudge_Force_Opt,           &
                          Nudge_TimeScale_Opt,                            &
                          Nudge_Times_Per_Day, Model_Times_Per_Day,       &
                          Nudge_Ucoef , Nudge_Uprof,                      &
                          Nudge_Vcoef , Nudge_Vprof,                      &
                          Nudge_Qcoef , Nudge_Qprof,                      &
                          Nudge_Tcoef , Nudge_Tprof,                      &
                          Nudge_PScoef, Nudge_PSprof,                     &
                          Nudge_Beg_Year, Nudge_Beg_Month, Nudge_Beg_Day, &
                          Nudge_End_Year, Nudge_End_Month, Nudge_End_Day, &
                          Nudge_Hwin_lat0, Nudge_Hwin_lon0,               &
                          Nudge_Hwin_latWidth, Nudge_Hwin_lonWidth,       &
                          Nudge_Hwin_latDelta, Nudge_Hwin_lonDelta,       &
                          Nudge_Hwin_Invert,                              &
                          Nudge_Vwin_Lindex, Nudge_Vwin_Hindex,           &
                          Nudge_Vwin_Ldelta, Nudge_Vwin_Hdelta,           &
                          Nudge_Vwin_Invert
 
    ! Nudging is NOT initialized yet, For now
    ! Nudging will always begin/end at midnight.
    !--------------------------------------------
    Nudge_Initialized =.false.
    Nudge_ON          =.false.
    Nudge_Beg_Sec=0
    Nudge_End_Sec=0

    ! Set Default Namelist values
    !-----------------------------
    Nudge_Model         = .false.
    Nudge_Path          = './Data/YOTC_ne30np4_001/'
    Nudge_File_Template = 'YOTC_ne30np4_L30.cam2.i.%y-%m-%d-%s.nc'
    Nudge_Force_Opt     = 0
    Nudge_TimeScale_Opt = 0
    Nudge_Times_Per_Day = 4
    Model_Times_Per_Day = 4
    Nudge_Ucoef         = 0._r8
    Nudge_Vcoef         = 0._r8
    Nudge_Qcoef         = 0._r8
    Nudge_Tcoef         = 0._r8
    Nudge_PScoef        = 0._r8
    Nudge_Uprof         = 0
    Nudge_Vprof         = 0
    Nudge_Qprof         = 0
    Nudge_Tprof         = 0
    Nudge_PSprof        = 0
    Nudge_Beg_Year      = 2010
    Nudge_Beg_Month     = 5
    Nudge_Beg_Day       = 1
    Nudge_End_Year      = 2010
    Nudge_End_Month     = 9
    Nudge_End_Day       = 1
    Nudge_Hwin_lat0     = 0._r8
    Nudge_Hwin_latWidth = 9999._r8
    Nudge_Hwin_latDelta = 1.0_r8
    Nudge_Hwin_lon0     = 180._r8
    Nudge_Hwin_lonWidth = 9999._r8
    Nudge_Hwin_lonDelta = 1.0_r8
    Nudge_Hwin_Invert   = .false.
    Nudge_Hwin_lo       = 0.0_r8
    Nudge_Hwin_hi       = 1.0_r8
    Nudge_Vwin_Hindex   = float(pver+1)
    Nudge_Vwin_Hdelta   = 0.001_r8
    Nudge_Vwin_Lindex   = 0.0_r8
    Nudge_Vwin_Ldelta   = 0.001_r8
    Nudge_Vwin_Invert   = .false.
    Nudge_Vwin_lo       = 0.0_r8
    Nudge_Vwin_hi       = 1.0_r8

    ! Read in namelist values
    !------------------------
    open(newunit=unitn, file=trim(nlfile), status='old')
    call find_group_name(unitn,'nudging_nl',status=ierr)
    if(ierr.eq.0) then
      read(unitn,nudging_nl,iostat=ierr)
      if(ierr.ne.0) then
        call endrun('nudging_readnl:: ERROR reading namelist')
      endif
    endif
    close(unitn)

    ! Initailize Boundary File and usage
    !-----------------------------------------------
    Nudge_BDY_File     = '/glade/scratch/patc/AMP_REPO/Boundary_NudgingWindows/' &
                         //'SAMwrf01_NudgingWindow_5X_5_5_0.999.nc'
!PFC    Nudge_Use_BDY_File = .true.
    Nudge_Use_BDY_File = .false.
    Nudge_BDY_Invert   = .true.
!xx    Nudge_BDY_File     = '/home/user/GRID_REPO/Boundary_NudgingWindows/' &
!xx                         //'SAMwrf01_NudgingWindow_5X_5_5_0.999.nc'
    
    ! Initialize a set of Dataset templates to use   
    !      (HARDCODED - convert this to use namelist entries)
    !------------------------------------------------
    OBS_FileTemplates%Num_Templates = 2
    allocate(OBS_FileTemplates%Path   (OBS_FileTemplates%Num_Templates))
    allocate(OBS_FileTemplates%Name   (OBS_FileTemplates%Num_Templates))
    allocate(OBS_FileTemplates%Type   (OBS_FileTemplates%Num_Templates))
    allocate(OBS_FileTemplates%TmapOpt(OBS_FileTemplates%Num_Templates))
!PFC NETCDF ERROR    OBS_FileTemplates%Path   (1) = '/glade/collections/rda/data/ds633.0/e5.oper.an.pl/'
!PFC NETCDF ERROR    OBS_FileTemplates%Name   (1) = '%y%m/e5.oper.an.pl.%v.%y%m%d00_%y%m%d23.nc'
    OBS_FileTemplates%Path   (1) = '/glade/scratch/patc/ERA5/e5.oper.an.pl/'
    OBS_FileTemplates%Name   (1) = '%y%m/e5.oper.an.pl.%v.%y%m%d00_%y%m%d23.nc'
!xx    OBS_FileTemplates%Path   (1) = '/home/user/nudging/ERA5/rda/e5.oper.an.pl/'
!xx    OBS_FileTemplates%Name   (1) = 'e5.oper.an.pl.%v.%y%m%d00_%y%m%d23.nc'
    OBS_FileTemplates%Type   (1) = 'ERA5'
    OBS_FileTemplates%TmapOpt(1) = 'LINEAR'
!PFC NETCDF ERROR    OBS_FileTemplates%Path   (2) = '/glade/collections/rda/data/ds633.0/e5.oper.an.sfc/'
!PFC NETCDF ERROR    OBS_FileTemplates%Name   (2) = '%y%m/e5.oper.an.sfc.%v.%y%m0100_%y%m%n23.nc'
    OBS_FileTemplates%Path   (2) = '/glade/scratch/patc/ERA5/e5.oper.an.sfc/'
    OBS_FileTemplates%Name   (2) = '%y%m/e5.oper.an.sfc.%v.%y%m0100_%y%m%n23.nc'
!xx    OBS_FileTemplates%Path   (2) = '/home/user/nudging/ERA5/rda/e5.oper.an.sfc/'
!xx    OBS_FileTemplates%Name   (2) = 'e5.oper.an.sfc.%v.%y%m0100_%y%m%n23.nc'
    OBS_FileTemplates%Type   (2) = 'ERA5'
    OBS_FileTemplates%TmapOpt(2) = 'LINEAR'

    ! Initialize the set of nudgined variables
    ! Set indicies for the DYN variables   (FIX THIS to set the indices from 'Name' values
    !      (HARDCODED - convert this to use namelist entries)
    !------------------------------------
    OBS_VarConfig%Num_Vars = 6
    allocate(OBS_VarConfig%Name        (OBS_VarConfig%Num_Vars))
    allocate(OBS_VarConfig%VmapOpt     (OBS_VarConfig%Num_Vars))
    allocate(OBS_VarConfig%HmapOpt     (OBS_VarConfig%Num_Vars))
    allocate(OBS_VarConfig%UseOpt      (OBS_VarConfig%Num_Vars))
    allocate(OBS_VarConfig%HORIZ_ONLY  (OBS_VarConfig%Num_Vars))
    allocate(OBS_VarConfig%SpaceProfile(OBS_VarConfig%Num_Vars))
    allocate(OBS_VarConfig%ScalingCoef (OBS_VarConfig%Num_Vars))

    idxU  = 1
    idxV  = 2
    idxT  = 3
    idxQ  = 4
    idxPS = 5
    idxT2 = 6
    OBS_VarConfig%Name        (idxU ) = 'U'
    OBS_VarConfig%VmapOpt     (idxU ) = 'INDEX_SHIFT' ! 'LOGP'
    OBS_VarConfig%HmapOpt     (idxU ) = 'BILINEAR_INTERP'
    OBS_VarConfig%UseOpt      (idxU ) = 'NUDGING'
    OBS_VarConfig%HORIZ_ONLY  (idxU ) = .false.
    OBS_VarConfig%SpaceProfile(idxU ) = Nudge_Uprof
    OBS_VarConfig%ScalingCoef (idxU ) = Nudge_Ucoef
    OBS_VarConfig%Name        (idxV ) = 'V'
    OBS_VarConfig%VmapOpt     (idxV ) = 'INDEX_SHIFT' ! 'LOGP'
    OBS_VarConfig%HmapOpt     (idxV ) = 'BILINEAR_INTERP'
    OBS_VarConfig%UseOpt      (idxV ) = 'NUDGING'
    OBS_VarConfig%HORIZ_ONLY  (idxV ) = .false.
    OBS_VarConfig%SpaceProfile(idxV ) = Nudge_Vprof
    OBS_VarConfig%ScalingCoef (idxV ) = Nudge_Vcoef
    OBS_VarConfig%Name        (idxT ) = 'T'
    OBS_VarConfig%VmapOpt     (idxT ) = 'INDEX_SHIFT' ! 'LOGP'
    OBS_VarConfig%HmapOpt     (idxT ) = 'BILINEAR_INTERP'
    OBS_VarConfig%UseOpt      (idxT ) = 'NUDGING'
    OBS_VarConfig%HORIZ_ONLY  (idxT ) = .false.
    OBS_VarConfig%SpaceProfile(idxT ) = Nudge_Tprof
    OBS_VarConfig%ScalingCoef (idxT ) = Nudge_Tcoef
    OBS_VarConfig%Name        (idxQ ) = 'Q'
    OBS_VarConfig%VmapOpt     (idxQ ) = 'INDEX_SHIFT' ! 'LOGP'
    OBS_VarConfig%HmapOpt     (idxQ ) = 'BILINEAR_INTERP'
    OBS_VarConfig%UseOpt      (idxQ ) = 'NUDGING'
    OBS_VarConfig%HORIZ_ONLY  (idxQ ) = .false.
    OBS_VarConfig%SpaceProfile(idxQ ) = Nudge_Qprof
    OBS_VarConfig%ScalingCoef (idxQ ) = Nudge_Qcoef
    OBS_VarConfig%Name        (idxPS) = 'PS'
    OBS_VarConfig%VmapOpt     (idxPS) = 'INDEX' 
    OBS_VarConfig%HmapOpt     (idxPS) = 'BILINEAR_INTERP'
    OBS_VarConfig%UseOpt      (idxPS) = 'NUDGING'
    OBS_VarConfig%HORIZ_ONLY  (idxPS) = .true.
    OBS_VarConfig%SpaceProfile(idxPS) = Nudge_PSprof
    OBS_VarConfig%ScalingCoef (idxPS) = Nudge_PScoef
    OBS_VarConfig%Name        (idxT2) = 'T2M'
    OBS_VarConfig%VmapOpt     (idxT2) = 'INDEX'
    OBS_VarConfig%HmapOpt     (idxT2) = 'BILINEAR_INTERP'
    OBS_VarConfig%UseOpt      (idxT2) = 'NUDGING'
    OBS_VarConfig%HORIZ_ONLY  (idxT2) = .true.
    OBS_VarConfig%SpaceProfile(idxT2) = 0._r8
    OBS_VarConfig%ScalingCoef (idxT2) = 0._r8

    ! Identify Mapping options that require inter-dependent Variables 
    !-----------------------------------------------------------------
    OBS_VarConfig%AdjOpt     = 'CAPT' ! 'CAPT'   ! 'NONE'
    OBS_VarConfig%Num_AdjVar = 6
    allocate(OBS_VarConfig%AdjVar(OBS_VarConfig%Num_AdjVar))
    OBS_VarConfig%AdjVar(1)  = 'U'
    OBS_VarConfig%AdjVar(2)  = 'V'
    OBS_VarConfig%AdjVar(3)  = 'PS'
    OBS_VarConfig%AdjVar(4)  = 'T'
    OBS_VarConfig%AdjVar(5)  = 'Q'
    OBS_VarConfig%AdjVar(6)  = 'T2M'

    ! Initialize a set of Filename templates for Invariant values (e.g.PHIS)
    !      (HARDCODED - convert this to use namelist entries)
    !------------------------------------------------------------------
    INV_FileTemplates%Num_Templates = 1
    allocate(INV_FileTemplates%Path(INV_FileTemplates%Num_Templates))
    allocate(INV_FileTemplates%Name(INV_FileTemplates%Num_Templates))
    allocate(INV_FileTemplates%Type(INV_FileTemplates%Num_Templates))
!PFC NETCDF ERROR    INV_FileTemplates%Path(1) = '/glade/collections/rda/data/ds633.0/e5.oper.invariant/'
!PFC NETCDF ERROR    INV_FileTemplates%Name(1) = '197901/e5.oper.invariant.%v.1979010100_1979010100.nc'
    INV_FileTemplates%Path(1) = '/glade/scratch/patc/ERA5/e5.oper.invariant/'
    INV_FileTemplates%Name(1) = '197901/e5.oper.invariant.%v.1979010100_1979010100.nc'
!xx    INV_FileTemplates%Path(1) = '/home/user/nudging/ERA5/rda/e5.oper.invariant/'
!xx    INV_FileTemplates%Name(1) = '197901/e5.oper.invariant.%v.1979010100_1979010100.nc'
    INV_FileTemplates%Type(1) = 'ERA5'

    ! Initialize the set of INVARIANT variables
    ! Set indicies for the variables
    !      (HARDCODED - convert this to use namelist entries)
    !------------------------------------
    INV_VarConfig%Num_Vars = 1
    allocate(INV_VarConfig%Name      (INV_VarConfig%Num_Vars))
    allocate(INV_VarConfig%VmapOpt   (INV_VarConfig%Num_Vars))
    allocate(INV_VarConfig%HmapOpt   (INV_VarConfig%Num_Vars))
    allocate(INV_VarConfig%HORIZ_ONLY(INV_VarConfig%Num_Vars))

    idxPHIS = 1
    INV_VarConfig%Name      (idxPHIS) = 'PHIS'
    INV_VarConfig%VmapOpt   (idxPHIS) = 'INDEX' 
    INV_VarConfig%HmapOpt   (idxPHIS) = 'BILINEAR_INTERP'
    INV_VarConfig%HORIZ_ONLY(idxPHIS) = .true.

    ! Set hi/lo values according to the given '_Invert' parameters
    !--------------------------------------------------------------
    if(Nudge_Hwin_Invert) then
      Nudge_Hwin_lo = 1.0_r8
      Nudge_Hwin_hi = 0.0_r8
    else
      Nudge_Hwin_lo = 0.0_r8
      Nudge_Hwin_hi = 1.0_r8
    endif

    if(Nudge_Vwin_Invert) then
      Nudge_Vwin_lo = 1.0_r8
      Nudge_Vwin_hi = 0.0_r8
    else
      Nudge_Vwin_lo = 0.0_r8
      Nudge_Vwin_hi = 1.0_r8
    endif

    ! Check for valid namelist values
    !----------------------------------
    if((Nudge_Hwin_lat0.lt.-90._r8).or.(Nudge_Hwin_lat0.gt.+90._r8)) then
      write(iulog,*) 'NUDGING: Window lat0 must be in [-90,+90]'
      write(iulog,*) 'NUDGING:  Nudge_Hwin_lat0=',Nudge_Hwin_lat0
      call endrun('nudging_readnl:: ERROR in namelist')
    endif

    if((Nudge_Hwin_lon0.lt.0._r8).or.(Nudge_Hwin_lon0.ge.360._r8)) then
      write(iulog,*) 'NUDGING: Window lon0 must be in [0,+360)'
      write(iulog,*) 'NUDGING:  Nudge_Hwin_lon0=',Nudge_Hwin_lon0
      call endrun('nudging_readnl:: ERROR in namelist')
    endif

    if((Nudge_Vwin_Lindex.gt.Nudge_Vwin_Hindex)                            .or. &
       (Nudge_Vwin_Hindex.gt.float(pver+1)).or.(Nudge_Vwin_Hindex.lt.0._r8).or. &
       (Nudge_Vwin_Lindex.gt.float(pver+1)).or.(Nudge_Vwin_Lindex.lt.0._r8)     ) then
      write(iulog,*) 'NUDGING: Window Lindex must be in [0,pver+1]'
      write(iulog,*) 'NUDGING: Window Hindex must be in [0,pver+1]'
      write(iulog,*) 'NUDGING: Lindex must be LE than Hindex'
      write(iulog,*) 'NUDGING:  Nudge_Vwin_Lindex=',Nudge_Vwin_Lindex
      write(iulog,*) 'NUDGING:  Nudge_Vwin_Hindex=',Nudge_Vwin_Hindex
      call endrun('nudging_readnl:: ERROR in namelist')
    endif

    if((Nudge_Hwin_latDelta.le.0._r8).or.(Nudge_Hwin_lonDelta.le.0._r8).or. &
       (Nudge_Vwin_Hdelta  .le.0._r8).or.(Nudge_Vwin_Ldelta  .le.0._r8)     ) then
      write(iulog,*) 'NUDGING: Window Deltas must be positive'
      write(iulog,*) 'NUDGING:  Nudge_Hwin_latDelta=',Nudge_Hwin_latDelta
      write(iulog,*) 'NUDGING:  Nudge_Hwin_lonDelta=',Nudge_Hwin_lonDelta
      write(iulog,*) 'NUDGING:  Nudge_Vwin_Hdelta=',Nudge_Vwin_Hdelta
      write(iulog,*) 'NUDGING:  Nudge_Vwin_Ldelta=',Nudge_Vwin_Ldelta
      call endrun('nudging_readnl:: ERROR in namelist')
    endif

    if((Nudge_Hwin_latWidth.le.0._r8).or.(Nudge_Hwin_lonWidth.le.0._r8)) then
      write(iulog,*) 'NUDGING: Window widths must be positive'
      write(iulog,*) 'NUDGING:  Nudge_Hwin_latWidth=',Nudge_Hwin_latWidth
      write(iulog,*) 'NUDGING:  Nudge_Hwin_lonWidth=',Nudge_Hwin_lonWidth
      call endrun('nudging_readnl:: ERROR in namelist')
    endif

    ! Initializae Nudging parameters 
    ! datastructure with namelist values
    !-----------------------------------------
    Nudge%Initialized   = Nudge_Initialized
    Nudge%Force_Opt     = Nudge_Force_Opt
    Nudge%TimeScale_Opt = Nudge_TimeScale_Opt
    Nudge%Beg_Year      = Nudge_Beg_Year
    Nudge%Beg_Month     = Nudge_Beg_Month
    Nudge%Beg_Day       = Nudge_Beg_Day
    Nudge%Beg_Sec       = Nudge_Beg_Sec
    Nudge%End_Year      = Nudge_End_Year
    Nudge%End_Month     = Nudge_End_Month
    Nudge%End_Day       = Nudge_End_Day
    Nudge%End_Sec       = Nudge_End_Sec
    Nudge%Times_Per_Day = Nudge_Times_Per_Day
    Nudge%Hwin_lat0     = Nudge_Hwin_lat0
    Nudge%Hwin_latWidth = Nudge_Hwin_latWidth
    Nudge%Hwin_latDelta = Nudge_Hwin_latDelta
    Nudge%Hwin_lon0     = Nudge_Hwin_lon0
    Nudge%Hwin_lonWidth = Nudge_Hwin_lonWidth
    Nudge%Hwin_lonDelta = Nudge_Hwin_lonDelta
    Nudge%Hwin_Invert   = Nudge_Hwin_Invert
    Nudge%Hwin_lo       = Nudge_Hwin_lo
    Nudge%Hwin_hi       = Nudge_Hwin_hi
    Nudge%Vwin_Hindex   = Nudge_Vwin_Hindex
    Nudge%Vwin_Hdelta   = Nudge_Vwin_Hdelta
    Nudge%Vwin_Lindex   = Nudge_Vwin_Lindex
    Nudge%Vwin_Ldelta   = Nudge_Vwin_Ldelta
    Nudge%Vwin_Invert   = Nudge_Vwin_Invert
    Nudge%Vwin_lo       = Nudge_Vwin_lo
    Nudge%Vwin_hi       = Nudge_Vwin_hi

    Nudge%use_BDY_file =      Nudge_Use_BDY_File
    Nudge%BDY_Invert   =      Nudge_BDY_Invert
    Nudge%BDY_File     = trim(Nudge_BDY_File)

    ! Initialize values for Horizontal window function
    !-------------------------------------------------
    lonp= 180._r8
    lon0=   0._r8
    lonn=-180._r8
    latp=  90._r8-Nudge%Hwin_lat0
    lat0=   0._r8
    latn= -90._r8-Nudge%Hwin_lat0

    Nudge%Hwin_lonWidthH=Nudge%Hwin_lonWidth/2._r8
    Nudge%Hwin_latWidthH=Nudge%Hwin_latWidth/2._r8

    Val1_p=(1._r8+tanh((Nudge%Hwin_lonWidthH+lonp)/Nudge%Hwin_lonDelta))/2._r8
    Val2_p=(1._r8+tanh((Nudge%Hwin_lonWidthH-lonp)/Nudge%Hwin_lonDelta))/2._r8
    Val3_p=(1._r8+tanh((Nudge%Hwin_latWidthH+latp)/Nudge%Hwin_latDelta))/2._r8
    Val4_p=(1._r8+tanh((Nudge%Hwin_latWidthH-latp)/Nudge%Hwin_latDelta))/2_r8
    Val1_0=(1._r8+tanh((Nudge%Hwin_lonWidthH+lon0)/Nudge%Hwin_lonDelta))/2._r8
    Val2_0=(1._r8+tanh((Nudge%Hwin_lonWidthH-lon0)/Nudge%Hwin_lonDelta))/2._r8
    Val3_0=(1._r8+tanh((Nudge%Hwin_latWidthH+lat0)/Nudge%Hwin_latDelta))/2._r8
    Val4_0=(1._r8+tanh((Nudge%Hwin_latWidthH-lat0)/Nudge%Hwin_latDelta))/2._r8

    Val1_n=(1._r8+tanh((Nudge%Hwin_lonWidthH+lonn)/Nudge%Hwin_lonDelta))/2._r8
    Val2_n=(1._r8+tanh((Nudge%Hwin_lonWidthH-lonn)/Nudge%Hwin_lonDelta))/2._r8
    Val3_n=(1._r8+tanh((Nudge%Hwin_latWidthH+latn)/Nudge%Hwin_latDelta))/2._r8
    Val4_n=(1._r8+tanh((Nudge%Hwin_latWidthH-latn)/Nudge%Hwin_latDelta))/2._r8

    Nudge%Hwin_max=     Val1_0*Val2_0*Val3_0*Val4_0
    Nudge%Hwin_min=min((Val1_p*Val2_p*Val3_n*Val4_n), &
                       (Val1_p*Val2_p*Val3_p*Val4_p), &
                       (Val1_n*Val2_n*Val3_n*Val4_n), &
                       (Val1_n*Val2_n*Val3_p*Val4_p)  )

    ! Initializae the time interval 
    ! for model state updates
    !--------------------------------------
    Model%Times_Per_Day = Model_Times_Per_Day

    ! End Routine
    !------------
    return
  end subroutine nudging_readnl
  !================================================================


  !================================================================
  subroutine nudging_init
    !
    ! NUDGING_INIT: Allocate space and initialize Nudging values
    !===============================================================
    use ppgrid          ,only: pver,pcols,begchunk,endchunk
    use error_messages  ,only: alloc_err
    use dycore          ,only: dycore_is
    use dyn_grid        ,only: get_horiz_grid_dim_d
    use cam_grid_support,only: cam_grid_id,cam_grid_get_dim_names,DLEN=>max_hcoordname_len
    use cam_history     ,only: addfld
    use shr_const_mod   ,only: SHR_CONST_PI
    use filenames       ,only: interpret_filename_spec
    use netcdf

    ! Local values
    !----------------
    integer :: Year,Month,Day,Sec
    integer :: YMD1,YMD
    logical :: After_Beg,Before_End
    integer :: istat,lchnk,ncol,icol,ilev
    integer :: ierr
    integer :: dtime
    real(r8):: rlat,rlon
    real(r8):: Wprof(pver)
    integer :: nn,tt

    integer             :: idxvar,ncid,kk
    integer             :: nlon,nlat
    real(r8),allocatable:: Bwindow(:,:)
    real(r8),allocatable:: B_lon(:)
    real(r8),allocatable:: B_lat(:)

    ! Get the time step size
    !-------------------------
    dtime = get_step_size()

    ! Check the time relative to the nudging window
    !------------------------------------------------
    call get_curr_date(Year,Month,Day,Sec)
    YMD=(Year*10000) + (Month*100) + Day
    YMD1=(Nudge%Beg_Year*10000) + (Nudge%Beg_Month*100) + Nudge%Beg_Day
    call timemgr_time_ge(YMD1,Nudge%Beg_Sec,         &
                         YMD ,Sec          ,After_Beg)
    YMD1=(Nudge%End_Year*10000) + (Nudge%End_Month*100) + Nudge%End_Day
    call timemgr_time_ge(YMD ,Sec          ,          &
                         YMD1,Nudge%End_Sec,Before_End)
 
    ! Allocate Space for Model State data arrays
    !-----------------------------------------
    allocate(Model%U(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'nudging_init','Model%U',pcols*pver*((endchunk-begchunk)+1))
    allocate(Model%V(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'nudging_init','Model%V',pcols*pver*((endchunk-begchunk)+1))
    allocate(Model%T(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'nudging_init','Model%T',pcols*pver*((endchunk-begchunk)+1))
    allocate(Model%S(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'nudging_init','Model%S',pcols*pver*((endchunk-begchunk)+1))
    allocate(Model%Q(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'nudging_init','Model%Q',pcols*pver*((endchunk-begchunk)+1))
    allocate(Model%PS(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'nudging_init','Model%PS',pcols*((endchunk-begchunk)+1))

    ! Set the Stepping intervals for Model and Nudging values
    ! Ensure that the Model_Step is not smaller then one timestep
    !  and not larger then the Nudge_Step.
    !--------------------------------------------------------
    Nudge%Step=86400/Nudge%Times_Per_Day   !>>>>>>       ?????? WHAT TO DO ABOUT THIS????
    Model%Step=86400/Model%Times_Per_Day
    if(Model%Step.lt.dtime) then
      write(iulog,*) ' '
      write(iulog,*) 'NUDGING: Model%Step cannot be less than a model timestep'
      write(iulog,*) 'NUDGING:  Setting Model%Step=dtime , dtime=',dtime
      write(iulog,*) ' '
      Model%Step=dtime
    endif
    if(Model%Step.gt.Nudge%Step) then
      write(iulog,*) ' '
      write(iulog,*) 'NUDGING: Model%Step cannot be more than Nudge_Step'
      write(iulog,*) 'NUDGING:  Setting Model%Step=Nudge_Step, Nudge_Step=',Nudge%Step
      write(iulog,*) ' '
      Model%Step=Nudge%Step
    endif
 
    if((After_Beg).and.(Before_End)) then
      ! Set Time indicies so that the next call to
      ! timestep_init will initialize the data arrays.
      !--------------------------------------------
      Model%Next_Year  = Year
      Model%Next_Month = Month
      Model%Next_Day   = Day
      Model%Next_Sec   = (Sec/Model%Step)*Model%Step
!>>>>      Nudge%Next_Year  = Year
!>>>>      Nudge%Next_Month = Month
!>>>>      Nudge%Next_Day   = Day
!>>>>      Nudge%Next_Sec   = (Sec/Nudge%Step)*Nudge%Step
    elseif(.not.After_Beg) then
      ! Set Time indicies to Nudging start,
      ! timestep_init will initialize the data arrays.
      !--------------------------------------------
      Model%Next_Year  = Nudge%Beg_Year
      Model%Next_Month = Nudge%Beg_Month
      Model%Next_Day   = Nudge%Beg_Day
      Model%Next_Sec   = Nudge%Beg_Sec
!>>>>      Nudge%Next_Year  = Nudge%Beg_Year
!>>>>      Nudge%Next_Month = Nudge%Beg_Month
!>>>>      Nudge%Next_Day   = Nudge%Beg_Day
!>>>>      Nudge%Next_Sec   = Nudge%Beg_Sec
    elseif(.not.Before_End) then
      ! Nudging will never occur, so switch it off
      !--------------------------------------------
      Nudge_Model = .false.
      Nudge_ON    = .false.
      write(iulog,*) ' '
      write(iulog,*) 'NUDGING: WARNING - Nudging has been requested by it will'
      write(iulog,*) 'NUDGING:           never occur for the given time values'
      write(iulog,*) ' '
    endif

    ! Initialize nudging observation OBS values to keep track of
    !-----------------------------------------------------------------
!x    call OBS%init(OBS_FileTemplates,OBS_VarConfig,Nudge%Beg_Year,  &
!x                                                  Nudge%Beg_Month, &
!x                                                  Nudge%Beg_Day    )
    call OBS%init(OBS_FileTemplates,OBS_VarConfig,2010,06,01)
    call OBS%print()

    ! Initialize Invariant Data values
    !-----------------------------------
    call INV%init(INV_FileTemplates, INV_VarConfig)
    call INV%print()

    ! Initialize boundary Coef values
    !---------------------------------
    if(Nudge%use_BDY_file) then
      if(allocated(Nudge%BDY_tau)) deallocate(Nudge%BDY_tau) 
      allocate(Nudge%BDY_tau(pcols,begchunk:endchunk))
      call initialize_boundary_tau(Nudge%use_BDY_file)
    endif

    ! Loop over OBS variables and create the associated
    ! Nudging values for each.
    !------------------------------------------------
    if(allocated(Nudge%var)) then
      do nn=1,Nudge%NumVar
        if(allocated(Nudge%var(nn)%Tau   )) deallocate(Nudge%var(nn)%Tau   )
        if(allocated(Nudge%var(nn)%Target)) deallocate(Nudge%var(nn)%Target)
        if(allocated(Nudge%var(nn)%Tend  )) deallocate(Nudge%var(nn)%Tend  )
      end do
      deallocate(Nudge%var)
    endif

    Nudge%NumVar = OBS_VarConfig%Num_Vars
    allocate(Nudge%var(Nudge%NumVar))
    do nn=1,OBS_VarConfig%Num_Vars
      Nudge%var(nn)%Name       = trim(OBS_VarConfig%Name        (nn))
      Nudge%var(nn)%UseOpt     = trim(OBS_VarConfig%UseOpt      (nn))
      Nudge%var(nn)%HORIZ_ONLY =      OBS_VarConfig%HORIZ_ONLY  (nn)
      Nudge%var(nn)%Prof       =      OBS_VarConfig%SpaceProfile(nn)
      Nudge%var(nn)%Coef       =      OBS_VarConfig%ScalingCoef (nn)
      if(Nudge%var(nn)%HORIZ_ONLY) then
        allocate(Nudge%var(nn)%Tau   (pcols,  1 ,begchunk:endchunk))
        allocate(Nudge%var(nn)%Target(pcols,  1 ,begchunk:endchunk))
        allocate(Nudge%var(nn)%Tend  (pcols,  1 ,begchunk:endchunk))
        Nudge%var(nn)%Tau   (:ncol,1,begchunk:endchunk) = 0._r8
        Nudge%var(nn)%Target(:ncol,1,begchunk:endchunk) = 0._r8
        Nudge%var(nn)%Tend  (:ncol,1,begchunk:endchunk) = 0._r8
      else
        allocate(Nudge%var(nn)%Tau   (pcols,pver,begchunk:endchunk))
        allocate(Nudge%var(nn)%Target(pcols,pver,begchunk:endchunk))
        allocate(Nudge%var(nn)%Tend  (pcols,pver,begchunk:endchunk))
        Nudge%var(nn)%Tau   (:ncol,:pver,begchunk:endchunk) = 0._r8
        Nudge%var(nn)%Target(:ncol,:pver,begchunk:endchunk) = 0._r8
        Nudge%var(nn)%Tend  (:ncol,:pver,begchunk:endchunk) = 0._r8
      endif
    end do ! nn=1,OBS%NumVar

    ! Initialize Nudging Coeffcient profiles in local arrays
    ! Load zeros into nudging arrays
    !------------------------------------------------------
    do lchnk=begchunk,endchunk
      ncol=get_ncols_p(lchnk)
      do icol=1,ncol
        rlat=get_rlat_p(lchnk,icol)*180._r8/SHR_CONST_PI
        rlon=get_rlon_p(lchnk,icol)*180._r8/SHR_CONST_PI
        do nn=1,Nudge%NumVar
          if(Nudge%var(nn)%HORIZ_ONLY) then
            Nudge%var(nn)%Tau(icol,1,lchnk)=set_horiz_profile(rlat,rlon,Nudge%var(nn)%Prof)
          else
            call set_nudging_profile(rlat,rlon,Nudge%var(nn)%Prof,Wprof,pver)
            Nudge%var(nn)%Tau(icol,:,lchnk)=Wprof(:)
          endif
        end do
      end do

      ! Scale the profiles with the given Coef values
      !-----------------------------------------------
      do nn=1,Nudge%NumVar
        if(Nudge%var(nn)%HORIZ_ONLY) then
          Nudge%var(nn)%Tau(:ncol,   1 ,lchnk) =                                    &
          Nudge%var(nn)%Tau(:ncol,   1 ,lchnk) * Nudge%var(nn)%Coef/float(Nudge%Step)
        else
          Nudge%var(nn)%Tau(:ncol,:pver,lchnk) =                                    &
          Nudge%var(nn)%Tau(:ncol,:pver,lchnk) * Nudge%var(nn)%Coef/float(Nudge%Step)
        endif
      end do

      ! Apply Boundary nudging conditions to all variables - merge with
      ! user-specified nudging by taking the MAX value.
      ! Initialize Tend/Target values to 0.
      !
      !  TODO: Make this more flexable for selective BC application
      !---------------------------------------------------------------
      if(Nudge%use_BDY_file) then
        do nn = 1,Nudge%NumVar
          if(Nudge%var(nn)%HORIZ_ONLY) then
            Nudge%var(nn)%Tau(:ncol,1,lchnk) = MAX(Nudge%var(nn)%Tau(:ncol,1,lchnk), &
                                                       Nudge%BDY_tau(:ncol,lchnk)    )
          else
            do kk = 1,pver
              Nudge%var(nn)%Tau(:ncol,kk,lchnk) = MAX(Nudge%var(nn)%Tau(:ncol,kk,lchnk), &
                                                          Nudge%BDY_tau(:ncol,lchnk)     )
            end do
          endif
        end do
      endif
    end do ! lchnk=begchunk,endchunk

    ! Register output fields with the cam history module
    !  TODO: Go Thru Nudge%var() and genetate addfld() calls
    !        only for active variables.
    !-----------------------------------------------------
    call addfld( 'Nudge_U',(/ 'lev' /),'A','m/s/s'  ,'U Nudging Tendency')
    call addfld( 'Nudge_V',(/ 'lev' /),'A','m/s/s'  ,'V Nudging Tendency')
    call addfld( 'Nudge_T',(/ 'lev' /),'A','K/s'    ,'T Nudging Tendency')
    call addfld( 'Nudge_Q',(/ 'lev' /),'A','kg/kg/s','Q Nudging Tendency')
    call addfld('Target_U',(/ 'lev' /),'A','m/s'    ,'U Nudging Target'  )
    call addfld('Target_V',(/ 'lev' /),'A','m/s'    ,'V Nudging Target'  )
    call addfld('Target_T',(/ 'lev' /),'A','K'      ,'T Nudging Target'  )
    call addfld('Target_Q',(/ 'lev' /),'A','kg/kg'  ,'Q Nudging Target  ')

    ! Initialization is done,
    !--------------------------
    Nudge%Initialized = .true.

    ! Informational Output
    !---------------------------
    if(masterproc) then
      write(iulog,*) ' '
      write(iulog,*) '---------------------------------------------------------'
      write(iulog,*) '  MODEL NUDGING INITIALIZED WITH THE FOLLOWING SETTINGS: '
      write(iulog,*) '---------------------------------------------------------'
      write(iulog,*) 'NUDGING: Nudge_Model='         ,Nudge_Model
!xx      write(iulog,*) 'NUDGING: Nudge_Path='          ,Nudge_Path
!xx      write(iulog,*) 'NUDGING: Nudge_File_Template =',Nudge_File_Template
      write(iulog,*) 'NUDGING: Nudge_Force_Opt='     ,Nudge%Force_Opt
      write(iulog,*) 'NUDGING: Nudge_TimeScale_Opt=' ,Nudge%TimeScale_Opt
      write(iulog,*) 'NUDGING: Nudge_Times_Per_Day=' ,Nudge%Times_Per_Day
      write(iulog,*) 'NUDGING: Model_Times_Per_Day=' ,Model%Times_Per_Day
      write(iulog,*) 'NUDGING: Nudge_Step='          ,Nudge%Step
      write(iulog,*) 'NUDGING: Model_Step='          ,Model%Step
      do nn=1,Nudge%NumVar
        write(iulog,*) 'NUDGING: ' &
                         //trim(Nudge%var(nn)%Name)//' Nudge_coef=',Nudge%var(nn)%Coef
        write(iulog,*) 'NUDGING: ' &
                         //trim(Nudge%var(nn)%Name)//' Nudge_prof=',Nudge%var(nn)%Prof
      end do
      write(iulog,*) 'NUDGING: Nudge_Beg_Year ='     ,Nudge%Beg_Year
      write(iulog,*) 'NUDGING: Nudge_Beg_Month='     ,Nudge%Beg_Month
      write(iulog,*) 'NUDGING: Nudge_Beg_Day  ='     ,Nudge%Beg_Day
      write(iulog,*) 'NUDGING: Nudge_End_Year ='     ,Nudge%End_Year
      write(iulog,*) 'NUDGING: Nudge_End_Month='     ,Nudge%End_Month
      write(iulog,*) 'NUDGING: Nudge_End_Day  ='     ,Nudge%End_Day
      write(iulog,*) 'NUDGING: Nudge_Hwin_lat0     =',Nudge%Hwin_lat0
      write(iulog,*) 'NUDGING: Nudge_Hwin_latWidth =',Nudge%Hwin_latWidth
      write(iulog,*) 'NUDGING: Nudge_Hwin_latDelta =',Nudge%Hwin_latDelta
      write(iulog,*) 'NUDGING: Nudge_Hwin_lon0     =',Nudge%Hwin_lon0
      write(iulog,*) 'NUDGING: Nudge_Hwin_lonWidth =',Nudge%Hwin_lonWidth
      write(iulog,*) 'NUDGING: Nudge_Hwin_lonDelta =',Nudge%Hwin_lonDelta
      write(iulog,*) 'NUDGING: Nudge_Hwin_Invert   =',Nudge%Hwin_Invert
      write(iulog,*) 'NUDGING: Nudge_Hwin_lo       =',Nudge%Hwin_lo
      write(iulog,*) 'NUDGING: Nudge_Hwin_hi       =',Nudge%Hwin_hi
      write(iulog,*) 'NUDGING: Nudge_Vwin_Hindex   =',Nudge%Vwin_Hindex
      write(iulog,*) 'NUDGING: Nudge_Vwin_Hdelta   =',Nudge%Vwin_Hdelta
      write(iulog,*) 'NUDGING: Nudge_Vwin_Lindex   =',Nudge%Vwin_Lindex
      write(iulog,*) 'NUDGING: Nudge_Vwin_Ldelta   =',Nudge%Vwin_Ldelta
      write(iulog,*) 'NUDGING: Nudge_Vwin_Invert   =',Nudge%Vwin_Invert
      write(iulog,*) 'NUDGING: Nudge_Vwin_lo       =',Nudge%Vwin_lo
      write(iulog,*) 'NUDGING: Nudge_Vwin_hi       =',Nudge%Vwin_hi
      write(iulog,*) 'NUDGING: Nudge_Hwin_latWidthH=',Nudge%Hwin_latWidthH
      write(iulog,*) 'NUDGING: Nudge_Hwin_lonWidthH=',Nudge%Hwin_lonWidthH
      write(iulog,*) 'NUDGING: Nudge_Hwin_max      =',Nudge%Hwin_max
      write(iulog,*) 'NUDGING: Nudge_Hwin_min      =',Nudge%Hwin_min
      write(iulog,*) 'NUDGING: Nudge_Initialized   =',Nudge%Initialized
      write(iulog,*) ' '
      write(iulog,*) ' '
    endif ! (masterproc) then

    ! End Routine
    !------------
    return
  end subroutine nudging_init
  !================================================================


  !================================================================
  subroutine nudging_timestep_init(phys_state)
    !
    ! NUDGING_TIMESTEP_INIT:
    !                 Check the current time and update Model/Nudging
    !                 arrays when necessary. Toggle the Nudging flag
    !                 when the time is withing the nudging window.
    !===============================================================
    use physconst    ,only: cpair
    use physics_types,only: physics_state
    use constituents ,only: cnst_get_ind
    use dycore       ,only: dycore_is
    use ppgrid       ,only: pver,pcols,begchunk,endchunk
    use filenames    ,only: interpret_filename_spec
    use ESMF

    ! Arguments
    !-----------
    type(physics_state),intent(in):: phys_state(begchunk:endchunk)

    ! Local values
    !----------------
    integer:: Year,Month,Day,Sec
    integer:: YMD1,YMD2,YMD
    logical:: Update_Model,Sync_Error
    logical:: After_Beg   ,Before_End
    integer:: lchnk,ncol,indw

    type(ESMF_Time)        :: Date1,Date2
    type(ESMF_TimeInterval):: DateDiff
    integer                :: DeltaT
    real(r8)               :: Tscale
    real(r8)               :: Tfrac
    integer                :: rc
    integer                :: nn
    integer                :: kk
    real(r8)               :: Sbar,Qbar,Wsum
    integer                :: dtime

    ! Check if Nudging is initialized
    !---------------------------------
    if(.not.Nudge%Initialized) then
      call endrun('nudging_timestep_init:: Nudging NOT Initialized')
    endif

    ! Get time step size
    !--------------------
    dtime = get_step_size()

    ! Get Current time
    !--------------------
    call get_curr_date(Year,Month,Day,Sec)
    YMD=(Year*10000) + (Month*100) + Day

    !-------------------------------------------------------
    ! Determine if the current time is AFTER the begining time
    ! and if it is BEFORE the ending time.
    !-------------------------------------------------------
    YMD1=(Nudge%Beg_Year*10000) + (Nudge%Beg_Month*100) + Nudge%Beg_Day
    call timemgr_time_ge(YMD1,Nudge%Beg_Sec,         &
                         YMD ,Sec          ,After_Beg)
 
    YMD1=(Nudge%End_Year*10000) + (Nudge%End_Month*100) + Nudge%End_Day
    call timemgr_time_ge(YMD ,Sec,                    &
                         YMD1,Nudge%End_Sec,Before_End)
 
    !--------------------------------------------------------------
    ! When past the NEXT time, Update Model Arrays and time indices
    !--------------------------------------------------------------
    YMD1=(Model%Next_Year*10000) + (Model%Next_Month*100) + Model%Next_Day
    call timemgr_time_ge(YMD1,Model%Next_Sec,            &
                         YMD ,Sec           ,Update_Model)

    if((Before_End).and.(Update_Model)) then
      ! Increment the Model times by the current interval
      !---------------------------------------------------
      Model%Curr_Year =Model%Next_Year
      Model%Curr_Month=Model%Next_Month
      Model%Curr_Day  =Model%Next_Day
      Model%Curr_Sec  =Model%Next_Sec
      YMD1=(Model%Curr_Year*10000) + (Model%Curr_Month*100) + Model%Curr_Day
      call timemgr_time_inc(YMD1,Model%Curr_Sec,              &
                            YMD2,Model%Next_Sec,Model%Step,0,0)
 
      ! Check for Sync Error where NEXT model time after the update
      ! is before the current time. If so, reset the next model
      ! time to a Model%Step after the current time.
      !--------------------------------------------------------------
      call timemgr_time_ge(YMD2,Model%Next_Sec,          &
                           YMD ,Sec           ,Sync_Error)
      if(Sync_Error) then
        Model%Curr_Year =Year
        Model%Curr_Month=Month
        Model%Curr_Day  =Day
        Model%Curr_Sec  =Sec
        call timemgr_time_inc(YMD ,Model%Curr_Sec,              &
                              YMD2,Model%Next_Sec,Model%Step,0,0)
        write(iulog,*) 'NUDGING: WARNING - Model%Time Sync ERROR... CORRECTED'
      endif
      Model%Next_Year =(YMD2/10000)
      YMD2            = YMD2-(Model%Next_Year*10000)
      Model%Next_Month=(YMD2/100)
      Model%Next_Day  = YMD2-(Model%Next_Month*100)

      ! Load values at Current into the Model arrays
      !-----------------------------------------------
      call cnst_get_ind('Q',indw)
      do lchnk=begchunk,endchunk
        ncol=phys_state(lchnk)%ncol
        Model%U(:ncol,:pver,lchnk)=phys_state(lchnk)%u(:ncol,:pver)
        Model%V(:ncol,:pver,lchnk)=phys_state(lchnk)%v(:ncol,:pver)
        Model%T(:ncol,:pver,lchnk)=phys_state(lchnk)%t(:ncol,:pver)
        Model%Q(:ncol,:pver,lchnk)=phys_state(lchnk)%q(:ncol,:pver,indw)
        Model%PS(:ncol,lchnk)=phys_state(lchnk)%ps(:ncol)
      end do

      ! Load Dry Static Energy values for Model
      !-----------------------------------------
      do lchnk=begchunk,endchunk
        ncol=phys_state(lchnk)%ncol
        Model%S(:ncol,:pver,lchnk)=cpair*Model%T(:ncol,:pver,lchnk)
      end do
    endif ! ((Before_End).and.(Update_Model)) then
 
    ! Update the time for OBS values
    !---------------------------------
    call OBS%setTime(Year,Month,Day,Sec)

    !----------------------------------------------------------------
    ! Toggle Nudging flag when the time interval is between
    ! beginning and ending times, and all of the analyses files exist.
    !----------------------------------------------------------------
    if((After_Beg).and.(Before_End)) then
      Nudge_ON=.true.
    else
      Nudge_ON=.false.
    endif

    !-------------------------------------------------------
    ! HERE Implement time dependence of Nudging Coefs HERE
    !-------------------------------------------------------


    !---------------------------------------------------
    ! If Data arrays have changed update stepping arrays
    !---------------------------------------------------
    if((Before_End).and.(Update_Model)) then

      ! Loop over Nudge Variables and update Target values
      ! with OBS data
      !--------------------------------------------------
      do nn = 1,Nudge%NumVar
        if(Nudge%var(nn)%HORIZ_ONLY) then
          call OBS%getData2D(trim(Nudge%var(nn)%Name),Nudge%var(nn)%Target)
        else
          call OBS%getData3D(trim(Nudge%var(nn)%Name),Nudge%var(nn)%Target)
        endif

        ! Internally we have DSE rather than Temperature, 
        ! so convert the input T values.
        !-------------------------------------------------
        if(trim(Nudge%var(nn)%Name).eq.'T') then
          Nudge%var(nn)%Target(:,:,:) = cpair*Nudge%var(nn)%Target(:,:,:)
        endif
      end do

      ! Set Tscale for the specified Forcing Option
      !-----------------------------------------------
      if(Nudge%TimeScale_Opt.eq.0) then
        Tscale=1._r8
      elseif(Nudge%TimeScale_Opt.eq.1) then
!z!FIX        call ESMF_TimeSet(Date1,YY=Year,MM=Month,DD=Day,S=Sec)
!z!FIX        call ESMF_TimeSet(Date2,YY=Nudge%Next_Year,MM=Nudge%Next_Month, &
!z!FIX                                DD=Nudge%Next_Day , S=Nudge%Next_Sec    )
!z!FIX        DateDiff =Date2-Date1
!z!FIX        call ESMF_TimeIntervalGet(DateDiff,S=DeltaT,rc=rc)
!z!FIX        Tscale=float(Nudge%Step)/float(DeltaT)
        call endrun('nudging_timestep_init:: NEED To Figure out Hoe to Implement this')
      else
        write(iulog,*) 'NUDGING: Unknown Nudge_TimeScale_Opt=',Nudge%TimeScale_Opt
        call endrun('nudging_timestep_init:: ERROR unknown Nudging_TimeScale_Opt')
      endif
 
      ! Update the nudging tendencies
      !--------------------------------
      do lchnk=begchunk,endchunk
        ncol=phys_state(lchnk)%ncol
        do nn = 1,Nudge%NumVar
          if(Nudge%var(nn)%HORIZ_ONLY) then
            if(trim(Nudge%var(nn)%Name).eq.'PS') then
               Nudge%var(nn)%Tend(:ncol,1,lchnk)=(Nudge%var(nn)%Target(:ncol,1,lchnk)  &
                                                             -Model%PS(:ncol,  lchnk)) &
                                             *Tscale*Nudge%var(nn)%Tau(:ncol,1,lchnk)
            endif
          else
            if(trim(Nudge%var(nn)%Name).eq.'U') then
              Nudge%var(nn)%Tend(:ncol,:pver,lchnk)=(Nudge%var(nn)%Target(:ncol,:pver,lchnk)  &
                                                                 -Model%U(:ncol,:pver,lchnk)) &
                                                *Tscale*Nudge%var(nn)%Tau(:ncol,:pver,lchnk)
            elseif(trim(Nudge%var(nn)%Name).eq.'V') then
              Nudge%var(nn)%Tend(:ncol,:pver,lchnk)=(Nudge%var(nn)%Target(:ncol,:pver,lchnk)  &
                                                                 -Model%V(:ncol,:pver,lchnk)) &
                                                *Tscale*Nudge%var(nn)%Tau(:ncol,:pver,lchnk)
            elseif(trim(Nudge%var(nn)%Name).eq.'T') then
              Nudge%var(nn)%Tend(:ncol,:pver,lchnk)=(Nudge%var(nn)%Target(:ncol,:pver,lchnk)  &
                                                                 -Model%T(:ncol,:pver,lchnk)) &
                                                *Tscale*Nudge%var(nn)%Tau(:ncol,:pver,lchnk)
            elseif(trim(Nudge%var(nn)%Name).eq.'Q') then
              Nudge%var(nn)%Tend(:ncol,:pver,lchnk)=(Nudge%var(nn)%Target(:ncol,:pver,lchnk)  &
                                                                 -Model%Q(:ncol,:pver,lchnk)) &
                                                *Tscale*Nudge%var(nn)%Tau(:ncol,:pver,lchnk)
            endif
          endif
        end do ! nn = 1,Nudge%NumVar
      end do ! lchnk=begchunk,endchunk
    endif ! ((Before_End).and.(Update_Model)) then

    ! End Routine
    !------------
    return
  end subroutine nudging_timestep_init
  !================================================================


  !================================================================
  subroutine nudging_timestep_tend(phys_state,phys_tend)
   !
   ! NUDGING_TIMESTEP_TEND:
   !                If Nudging is ON, return the Nudging contributions
   !                to forcing using the current contents of the Nudge
   !                arrays. Send output to the cam history module as well.
   !===============================================================
   use physconst    ,only: cpair
   use physics_types,only: physics_state,physics_ptend,physics_ptend_init
   use constituents ,only: cnst_get_ind,pcnst
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use cam_history  ,only: outfld

   ! Arguments
   !-------------
   type(physics_state),intent(in) :: phys_state
   type(physics_ptend),intent(out):: phys_tend

   ! Local values
   !--------------------
   integer:: indw,ncol,lchnk,nn
   logical:: lq(pcnst)

   call cnst_get_ind('Q',indw)
   lq(:)   =.false.
   lq(indw)=.true.
   call physics_ptend_init(phys_tend,phys_state%psetcols,               &
                           'nudging',lu=.true.,lv=.true.,ls=.true.,lq=lq)

   if(Nudge_ON) then
     lchnk=phys_state%lchnk
     ncol =phys_state%ncol

     do nn = 1,Nudge%NumVar
       if(trim(Nudge%var(nn)%Name).eq.'U') then
         phys_tend%u(:ncol,:pver) = Nudge%var(nn)%Tend(:ncol,:pver,lchnk)
       elseif(trim(Nudge%var(nn)%Name).eq.'V') then
         phys_tend%v(:ncol,:pver) = Nudge%var(nn)%Tend(:ncol,:pver,lchnk)
       elseif(trim(Nudge%var(nn)%Name).eq.'T') then
         phys_tend%s(:ncol,:pver) = Nudge%var(nn)%Tend(:ncol,:pver,lchnk)
       elseif(trim(Nudge%var(nn)%Name).eq.'Q') then
         phys_tend%q(:ncol,:pver,indw) = Nudge%var(nn)%Tend(:ncol,:pver,lchnk)
       endif
     end do ! nn = 1,Nudge%NumVar

     call outfld( 'Nudge_U',phys_tend%u                      ,pcols,lchnk)
     call outfld( 'Nudge_V',phys_tend%v                      ,pcols,lchnk)
     call outfld( 'Nudge_T',phys_tend%s/cpair                ,pcols,lchnk)
     call outfld( 'Nudge_Q',phys_tend%q(1,1,indw)            ,pcols,lchnk)
     call outfld('Target_U',Nudge%var(idxU)%Target(:,:,lchnk),pcols,lchnk)
     call outfld('Target_V',Nudge%var(idxV)%Target(:,:,lchnk),pcols,lchnk)
     call outfld('Target_T',Nudge%var(idxT)%Target(:,:,lchnk),pcols,lchnk)
     call outfld('Target_Q',Nudge%var(idxQ)%Target(:,:,lchnk),pcols,lchnk)
   endif

   ! End Routine
   !------------
   return
  end subroutine nudging_timestep_tend
  !================================================================


  !================================================================
  subroutine set_nudging_profile(rlat,rlon,Nudge_prof,Wprof,nlev)
   !
   ! SET_NUDGING_PROFILE: for the given lat,lon, and Nudging_prof, set
   !                      the verical profile of window coeffcients.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on nudging strength.
   !===============================================================

   ! Arguments
   !--------------
   integer :: nlev,Nudge_prof
   real(r8):: rlat,rlon
   real(r8):: Wprof(nlev)

   ! Local values
   !----------------
   integer :: ilev
   real(r8):: Hcoef,latx,lonx,Vmax,Vmin
   real(r8):: lon_lo,lon_hi,lat_lo,lat_hi,lev_lo,lev_hi

   !---------------
   ! set coeffcient
   !---------------
   if(Nudge_prof.eq.0) then
     ! No Nudging
     !-------------
     Wprof(:)=0.0_r8
   elseif(Nudge_prof.eq.1) then
     ! Uniform Nudging
     !-----------------
     Wprof(:)=1.0_r8
   elseif(Nudge_prof.eq.2) then
     ! Localized Nudging with specified Heaviside window function
     !------------------------------------------------------------
     if(Nudge%Hwin_max.le.Nudge%Hwin_min) then
       ! For a constant Horizontal window function,
       ! just set Hcoef to the maximum of Hlo/Hhi.
       !--------------------------------------------
       Hcoef=max(Nudge%Hwin_lo,Nudge%Hwin_hi)
     else
       ! get lat/lon relative to window center
       !------------------------------------------
       latx=rlat-Nudge%Hwin_lat0
       lonx=rlon-Nudge%Hwin_lon0
       if(lonx.gt. 180._r8) lonx=lonx-360._r8
       if(lonx.le.-180._r8) lonx=lonx+360._r8

       ! Calcualte RAW window value
       !-------------------------------
       lon_lo=(Nudge%Hwin_lonWidthH+lonx)/Nudge%Hwin_lonDelta
       lon_hi=(Nudge%Hwin_lonWidthH-lonx)/Nudge%Hwin_lonDelta
       lat_lo=(Nudge%Hwin_latWidthH+latx)/Nudge%Hwin_latDelta
       lat_hi=(Nudge%Hwin_latWidthH-latx)/Nudge%Hwin_latDelta
       Hcoef=((1._r8+tanh(lon_lo))/2._r8)*((1._r8+tanh(lon_hi))/2._r8) &
            *((1._r8+tanh(lat_lo))/2._r8)*((1._r8+tanh(lat_hi))/2._r8)

       ! Scale the horizontal window coef for specfied range of values.
       !--------------------------------------------------------
       Hcoef=(Hcoef-Nudge%Hwin_min)/(Nudge%Hwin_max-Nudge%Hwin_min)
       Hcoef=(1._r8-Hcoef)*Nudge%Hwin_lo + Hcoef*Nudge%Hwin_hi
     endif

     ! Load the RAW vertical window
     !------------------------------
     do ilev=1,nlev
       lev_lo=(float(ilev)-Nudge%Vwin_Lindex)/Nudge%Vwin_Ldelta
       lev_hi=(Nudge%Vwin_Hindex-float(ilev))/Nudge%Vwin_Hdelta
       Wprof(ilev)=((1._r8+tanh(lev_lo))/2._r8)*((1._r8+tanh(lev_hi))/2._r8)
     end do

     ! Scale the Window function to span the values between Vlo and Vhi:
     !-----------------------------------------------------------------
     Vmax=maxval(Wprof)
     Vmin=minval(Wprof)
     if((Vmax.le.Vmin).or.((Nudge%Vwin_Hindex.ge.(nlev+1)).and. &
                           (Nudge%Vwin_Lindex.le. 0      )     )) then
       ! For a constant Vertical window function,
       ! load maximum of Vlo/Vhi into Wprof()
       !--------------------------------------------
       Vmax=max(Nudge%Vwin_lo,Nudge%Vwin_hi)
       Wprof(:)=Vmax
     else
       ! Scale the RAW vertical window for specfied range of values.
       !--------------------------------------------------------
       Wprof(:)=(Wprof(:)-Vmin)/(Vmax-Vmin)
       Wprof(:)=Nudge%Vwin_lo + Wprof(:)*(Nudge%Vwin_hi-Nudge%Vwin_lo)
     endif

     ! The desired result is the product of the vertical profile
     ! and the horizontal window coeffcient.
     !----------------------------------------------------
     Wprof(:)=Hcoef*Wprof(:)
   else
     call endrun('set_nudging_profile:: Unknown Nudge_prof value')
   endif

   ! End Routine
   !------------
   return
  end subroutine set_nudging_profile
  !================================================================


  !================================================================
  real(r8) function set_horiz_profile(rlat,rlon,Nudge_prof)
   !
   ! SET_HORIZ_PROFILE: for the given lat and lon set the surface
   !                    pressure profile value for the specified index.
   !                    Values range from 0. to 1. to affect spatial
   !                    variations on nudging strength.
   !===============================================================

   ! Arguments
   !--------------
   real(r8):: rlat,rlon
   integer :: Nudge_prof

   ! Local values
   !----------------

   !---------------
   ! set coeffcient
   !---------------
   if(Nudge_prof.eq.0) then
     ! No Nudging
     !-------------
     set_horiz_profile=0.0_r8
   elseif(Nudge_prof.eq.1) then
     ! Uniform Nudging
     !-----------------
     set_horiz_profile=1.0_r8
   else
     call endrun('set_horiz_profile:: Unknown Nudge_prof value')
   endif

   ! End Routine
   !------------
   return
  end function set_horiz_profile
  !================================================================


  !================================================================
  subroutine calc_DryStaticEnergy(t, q, phis, ps, dse, ncol)
   !
   ! calc_DryStaticEnergy: Given the temperature, specific humidity, surface pressure,
   !                       and surface geopotential for a chunk containing 'ncol' columns,
   !                       calculate and return the corresponding dry static energy values.
   !--------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,       only: pver, pverp
   use dycore,       only: dycore_is
   use hycoef,       only: hyai, hybi, ps0, hyam, hybm
   use physconst,    only: zvir, gravit, cpair, rair
   !
   ! Input/Output arguments
   !-----------------------
   integer , intent(in) :: ncol      ! Number of columns in chunk
   real(r8), intent(in) :: t(:,:)    ! (pcols,pver) - temperature
   real(r8), intent(in) :: q(:,:)    ! (pcols,pver) - specific humidity
   real(r8), intent(in) :: ps(:)     ! (pcols)      - surface pressure
   real(r8), intent(in) :: phis(:)   ! (pcols)      - surface geopotential
   real(r8), intent(out):: dse(:,:)  ! (pcols,pver)  - dry static energy
   !
   ! Local variables
   !------------------
   logical :: fvdyn                 ! finite volume dynamics
   integer :: ii,kk                 ! Lon, level, level indices
   real(r8):: tvfac                 ! Virtual temperature factor
   real(r8):: hkk(ncol)             ! diagonal element of hydrostatic matrix
   real(r8):: hkl(ncol)             ! off-diagonal element
   real(r8):: pint(ncol,pverp)      ! Interface pressures
   real(r8):: pmid(ncol,pver )      ! Midpoint pressures
   real(r8):: zi(ncol,pverp)        ! Height above surface at interfaces
   real(r8):: zm(ncol,pver )        ! Geopotential height at mid level

   ! Set dynamics flag
   !-------------------
   fvdyn = dycore_is ('LR')

   ! Load Pressure values and midpoint pressures
   !----------------------------------------------
   do kk=1,pverp
   do ii=1,ncol
     pint(ii,kk)=(hyai(kk)*ps0)+(hybi(kk)*ps(ii))
   end do
   end do

   do kk=1,pver
   do ii=1,ncol
     pmid(ii,kk)=(hyam(kk)*ps0)+(hybm(kk)*ps(ii))
   end do
   end do

   ! The surface height is zero by definition.
   !-------------------------------------------
   do ii = 1,ncol
     zi(ii,pverp) = 0.0_r8
   end do

   ! Compute the dry static energy, zi, zm from bottom up
   ! Note, zi(i,k) is the interface above zm(i,k)
   !---------------------------------------------------------
   do kk=pver,1,-1

     ! First set hydrostatic elements consistent with dynamics
     !--------------------------------------------------------
     if(fvdyn) then
       do ii=1,ncol
         hkl(ii)=log(pint(ii,kk+1))-log(pint(ii,kk))
         hkk(ii)=1._r8-(hkl(ii)*pint(ii,kk)/(pint(ii,kk+1)-pint(ii,kk)))
       end do
     else
       do ii=1,ncol
         hkl(ii)=(pint(ii,kk+1)-pint(ii,kk))/pmid(ii,kk)
         hkk(ii)=0.5_r8*hkl(ii)
       end do
     endif

     ! Now compute zm, zi, and dse  (WACCM-X vars rairv/zairv/cpairv not used!)
     !------------------------------------------------------------------------
     do ii=1,ncol
       tvfac=t(ii,kk)*rair*(1._r8+(zvir*q(ii,kk)))/gravit
       zm (ii,kk)=zi(ii,kk+1) + (tvfac*hkk(ii))
       zi (ii,kk)=zi(ii,kk+1) + (tvfac*hkl(ii))
       dse(ii,kk)=(t(ii,kk)*cpair) + phis(ii) + (gravit*zm(ii,kk))
     end do

   end do ! kk=pver,1,-1

   ! End Routine
   !-----------
   return
  end subroutine calc_DryStaticEnergy
  !================================================================


  !================================================================
  subroutine initialize_boundary_tau(Use_BDY_File)
    ! 
    !==============================================================
    use shr_const_mod ,only: SHR_CONST_PI
    !
    ! Passed Variables
    !--------------------
    logical,intent(in):: Use_BDY_File
    !
    ! Local Values
    !-------------
    integer :: istat, ncid, idxvar, nlon, nlat
    integer :: lchnk, ncol, icol
    real(r8):: rlat, rlon
    real(r8),allocatable:: Bwindow(:,:)
    real(r8),allocatable:: B_lon(:)
    real(r8),allocatable:: B_lat(:)

    ! If not using a boundary file, then BDY_tau remains 
    ! unallocated and un used. Just return.
    !----------------------------------------------------
 write(iulog,*) ' PFCDIAG: Use_BDY_File=',Use_BDY_File
    if(.not.Use_BDY_File) then
      return
    endif

    ! Read in Boundary Nudging Window, 
    ! interpolate values from lat/lon grid, 
    !----------------------------------
    istat = nf90_open(trim(Nudge%BDY_File),NF90_NOWRITE,ncid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      write(iulog,*) 'NF90_OPEN: failed for file:'//trim(Nudge%BDY_File)
      call endrun('subroutine nudging_init')
    else
      if(masterproc) then
        write(iulog,*) 'NUDGING: Opened Boundary Window File:'//trim(Nudge%BDY_File) 
      endif
    endif

    istat=nf90_inq_dimid(ncid,'longitude',idxvar)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun('subroutine nudging_init')
    endif
    istat=nf90_inquire_dimension(ncid,idxvar,len=nlon)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun('subroutine nudging_init')
    endif

    istat=nf90_inq_dimid(ncid,'latitude',idxvar)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun('subroutine nudging_init')
    endif
    istat=nf90_inquire_dimension(ncid,idxvar,len=nlat)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun('subroutine nudging_init')
    endif

    allocate(Bwindow(nlon,nlat))
    allocate(B_lon(nlon))
    allocate(B_lat(nlat))

    istat=nf90_inq_varid(ncid,'longitude',idxvar)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun('subroutine nudging_init')
    endif
    istat=nf90_get_var(ncid,idxvar,B_lon)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun('subroutine nudging_init')
    endif

    istat=nf90_inq_varid(ncid,'latitude',idxvar)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun('subroutine nudging_init')
    endif
    istat=nf90_get_var(ncid,idxvar,B_lat)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun('subroutine nudging_init')
    endif

    istat=nf90_inq_varid(ncid,'refineMap',idxvar)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun('subroutine nudging_init')
    endif
    istat=nf90_get_var(ncid,idxvar,Bwindow)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun('subroutine nudging_init')
    endif
    istat=nf90_close(ncid)

    do lchnk=begchunk,endchunk
      ncol=get_ncols_p(lchnk)
      do icol=1,ncol
        rlat=get_rlat_p(lchnk,icol)*180._r8/SHR_CONST_PI
        rlon=get_rlon_p(lchnk,icol)*180._r8/SHR_CONST_PI
        if(rlon.lt.  0._r8) rlon = rlon + 360._r8
        if(rlon.gt.360._r8) rlon = rlon - 360._r8
 
        Nudge%BDY_tau(icol,lchnk) = interp_Bwin(rlon,rlat,B_lon,B_lat,Bwindow,nlon,nlat)
        if(Nudge%BDY_Invert) then
          ! INVERT WINDOW
          !---------------
          Nudge%BDY_tau(icol,lchnk) = 1._r8 - Nudge%BDY_tau(icol,lchnk)
        endif
  
        ! CONSISTENT SCALING
        !-------------------
        Nudge%BDY_tau(icol,lchnk) = Nudge%BDY_tau(icol,lchnk)/float(Nudge%Step) 
      end do
    end do

    deallocate(Bwindow)
    deallocate(B_lon)
    deallocate(B_lat)

    ! End Routine
    !------------
    return
  end subroutine initialize_boundary_tau
  !================================================================


  !================================================================
  real(r8) function interp_Bwin(rlon,rlat,B_lon,B_lat,Bwindow,nlon,nlat)
   ! 
   ! interp_Bwin: for the given lat and lon return the interpolated
   !              Bounary Nudinging coef from the given lat,lon 
   !              array of coeffcients.
   !===============================================================

   ! Arguments
   !--------------
   integer  nlat,nlon
   real(r8) rlat,rlon
   real(r8) B_lon  (nlon)
   real(r8) B_lat  (nlat)
   real(r8) Bwindow(nlon,nlat)

   ! Local values
   !----------------
   integer  ii
   integer  LonLo,LonHi,LatLo,LatHi
   integer  LonLx,LonHx,LatLx,LatHx
   real(r8) LonFrac,LatFrac

   ! Get the Bounding indices for longitude
   !----------------------------------------
   if(rlon.lt.B_lon(1)) then
     LonLo   = nlon
     LonHi   = 1
     LonFrac = ( rlon   +360._r8 - B_lon(nlon)) &
              /(B_lon(1)+360._r8 - B_lon(nlon))
   elseif(rlon.gt.B_lon(nlon)) then
     LonLo   = nlon
     LonHi   = 1
     LonFrac = ( rlon            - B_lon(nlon)) &
              /(B_lon(1)+360._r8 - B_lon(nlon))
   else
     do ii=1,(nlon-1)
       if((rlon.ge.B_lon(ii)).and.(rlon.lt.B_lon(ii+1))) then
         LonLo   = ii
         LonHi   = ii+1
         LonFrac = ( rlon       - B_lon(ii)) &
                  /(B_lon(ii+1) - B_lon(ii))
         exit
       endif
     end do
   endif

   ! Get the Bounding indices for latitude
   !---------------------------------------
   if(rlat.lt.B_lat(1)) then
     LatLo   = 0
     LatHi   = 1
     LatFrac = (rlat + B_lat(1))/(2._r8*B_lat(1))
   elseif(rlat.gt.B_lat(nlat)) then
     LatLo   = nlat
     LatHi   = nlat + 1
     LatFrac =  ( rlat         - B_lat(nlat)) &
               /(180._r8 - 2._r8*B_lat(nlat))
   else
     do ii=1,(nlat-1)
       if((rlat.ge.B_lat(ii)).and.(rlat.lt.B_lat(ii+1))) then
         LatLo   = ii
         LatHi   = ii+1
         LatFrac = ( rlat       - B_lat(ii)) &
                  /(B_lat(ii+1) - B_lat(ii))
         exit
       endif
     end do
   endif

   ! Interpolate
   !-------------
   if(LatHi.gt.nlat) then
     LatHx = LatLo
     LonLx = LonHi + (nlon/2)
     LonHx = LonLo + (nlon/2)
     if(LonLx.gt.nlon)  LonLx = LonLx - nlon
     if(LonHx.gt.nlon)  LonHx = LonHx - nlon
     interp_Bwin = Bwindow(LonLo,LatLo)*(1._r8-LonFrac)*(1._r8-LatFrac) &
                  +Bwindow(LonHi,LatLo)*       LonFrac *(1._r8-LatFrac) &
                  +Bwindow(LonHx,LatHx)*       LonFrac *       LatFrac  &
                  +Bwindow(LonLx,LatHx)*(1._r8-LonFrac)*       LatFrac
   elseif(LatLo.lt. 1  ) then
     LatLx = LatHi
     LonLx = LonHi + (nlon/2)
     LonHx = LonLo + (nlon/2)
     if(LonLx.gt.nlon)  LonLx = LonLx - nlon
     if(LonHx.gt.nlon)  LonHx = LonHx - nlon
     interp_Bwin = Bwindow(LonLx,LatLx)*(1._r8-LonFrac)*(1._r8-LatFrac) &
                  +Bwindow(LonHx,LatLx)*       LonFrac *(1._r8-LatFrac) &
                  +Bwindow(LonHi,LatHi)*       LonFrac *       LatFrac  &
                  +Bwindow(LonLo,LatHi)*(1._r8-LonFrac)*       LatFrac
   else
     interp_Bwin = Bwindow(LonLo,LatLo)*(1._r8-LonFrac)*(1._r8-LatFrac) &
                  +Bwindow(LonHi,LatLo)*       LonFrac *(1._r8-LatFrac) &
                  +Bwindow(LonHi,LatHi)*       LonFrac *       LatFrac  &
                  +Bwindow(LonLo,LatHi)*(1._r8-LonFrac)*       LatFrac
   endif

   ! End Routine
   !------------
   return
  end function interp_Bwin
  !================================================================


  !================================================================
  subroutine lookup_DataType(FileType,VarName,FileString,DSname)
    !
    !==============================================================
    ! 
    ! Passed variables 
    !-------------------
    character(len=*),intent(in ):: FileType
    character(len=*),intent(in ):: VarName
    character(len=cl),intent(out):: FileString
    character(len=16),intent(out):: DSname

    ! Get datastructure values for the given FileType
    !-------------------------------------------------------
    if(trim(FileType).eq.'ERA5') then
      if(trim(VarName).eq.'U') then
        FileString = '128_131_u.ll025uv'
        DSname     = trim(VarName)
      elseif(trim(VarName).eq.'V') then
        FileString = '128_132_v.ll025uv'
        DSname     = trim(VarName)
      elseif(trim(VarName).eq.'T') then
        FileString = '128_130_t.ll025sc'
        DSname     = trim(VarName)
      elseif(trim(VarName).eq.'Q') then
        FileString = '128_133_q.ll025sc'
        DSname     = trim(VarName)
      elseif(trim(VarName).eq.'PS') then
        FileString = '128_134_sp.ll025sc'
        DSname     = 'SP'
      elseif(trim(VarName).eq.'T2M') then
        FileString = '128_167_2t.ll025sc'
        DSname     = 'VAR_2T'
      elseif(trim(VarName).eq.'PHIS') then
        FileString = '128_129_z.ll025sc'
        DSname     = 'Z'
      endif
    elseif(trim(FileType).eq.'ERAI') then
      if(trim(VarName).eq.'U') then
        FileString = 'NONE'
        DSname     = trim(VarName)
      elseif(trim(VarName).eq.'V') then
        FileString = 'NONE'
        DSname     = trim(VarName)
      elseif(trim(VarName).eq.'T') then
        FileString = 'NONE'
        DSname     = trim(VarName)
      elseif(trim(VarName).eq.'Q') then
        FileString = 'NONE'
        DSname     = trim(VarName)
      elseif(trim(VarName).eq.'PS') then
        FileString = 'NONE'
        DSname     = trim(VarName)
      elseif(trim(VarName).eq.'T2M') then
        FileString = 'NONE'
        DSname     = trim(VarName)
      elseif(trim(VarName).eq.'PHIS') then
        FileString = 'NONE'
        DSname     = 'Z'
      endif
    elseif(trim(FileType).eq.'CESM') then
      if(trim(VarName).eq.'U') then
        FileString = 'NONE'
        DSname     = trim(VarName)
      elseif(trim(VarName).eq.'V') then
        FileString = 'NONE'
        DSname     = trim(VarName)
      elseif(trim(VarName).eq.'T') then
        FileString = 'NONE'
        DSname     = trim(VarName)
      elseif(trim(VarName).eq.'Q') then
        FileString = 'NONE'
        DSname     = trim(VarName)
      elseif(trim(VarName).eq.'PS') then
        FileString = 'NONE'
        DSname     = trim(VarName)
      elseif(trim(VarName).eq.'T2M') then
        FileString = 'NONE'
        DSname     = trim(VarName)
      elseif(trim(VarName).eq.'PHIS') then
        FileString = 'NONE'
        DSname     = 'PHIS'
      endif
    endif

    ! End Routine
    !------------
    return
  end subroutine lookup_DataType
  !================================================================


  !================================================================
  character(len=cl) function replace_string(Template,StringCode,String,StringLen)
    !
    ! replace_string: Scan the Template for the given string code, 
    !                 if found, replace the code with the given String
    !
    !       "%v" = code for Variable-dependent String
    !       "%n" = code for NN the number of days in the month
    !       "%y" = code for a YYYY  Year string
    !       "%m" = code for a MM    Month string.
    !       "%d" = code for a DD    Day string.
    !       "%h" = code for a HH    Hour string.
    !       "%s" = code for a SSSSS Second string.
    !==============================================================
    ! 
    ! Passed variables
    !-----------------
    character(len=*),intent(in):: Template
    character(len=*),intent(in):: StringCode
    character(len=*),intent(in):: String
    integer         ,intent(in):: StringLen
    !
    ! Local values
    !--------------
    character(len=cl):: tmp_string
    character(len=cl):: new_string
    integer          :: len_string
    logical          :: scan_not_done
    logical          :: code_found
    integer          :: ii
    
    ! Initialize the new string with the tmeplate
    !---------------------------------------------
    len_string = len_trim(Template)
    new_string = trim(Template)

    ! Replace all occurances of the StingCode with 
    ! the given String(1:StingLen)
    !---------------------------------------------
    scan_not_done = .true.
    do while (scan_not_done) 
      code_found = .false.
      do ii=1,(len_string-1)
        if(new_string(ii:ii+1).eq.StringCode(1:2)) then
          if(ii.eq.1) then
            tmp_string = String(1:StringLen)//new_string(3:len_string)
          elseif(ii.eq.len_string-1) then
            tmp_string = new_string(1:len_string-2)//String(1:StringLen)
          else
            tmp_string = new_string(1:ii-1)//String(1:StringLen)//new_string(ii+2:len_string)
          endif
          code_found = .true.
          exit
        endif
        if(ii.eq.(len_string-1)) scan_not_done = .false.
      end do
  
      if(code_found) then
        new_string = trim(tmp_string)
        len_string = len_trim(new_string)
      endif
    end do ! while (scan_not_done) 

    replace_string = trim(new_string)

    ! End Routine
    !------------
    return
  end function replace_string
  !================================================================


  !================================================================
  integer function days_per_month(Year,Month)
    !
    ! days_per_month: Simple function to return the number of days
    !                 for the given month year.
    !
    !==============================================================
    ! 
    ! Passed variables
    !-----------------
    integer,intent(in):: Year
    integer,intent(in):: Month
    !
    ! Local values
    !--------------
    integer:: month_days(12)
    data month_days /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

    days_per_month = month_days(Month)
    if((Year.eq.(4*(Year/4))).and.(Month.eq.2)) days_per_month=29

    ! End Routine
    !------------
    return
  end function days_per_month
  !================================================================

end module nudging
