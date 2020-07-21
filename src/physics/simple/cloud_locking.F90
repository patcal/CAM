module cloud_locking
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
!    analyses. The form of the target values is selected by the 'Clock_Force_Opt'
!    option, the timescale of the forcing is determined from the given 
!    'Clock_TimeScale_Opt', and the nudging strength Alpha=[0.,1.] for each 
!    variable is specified by the 'Clock_Xcoef' values. Where X={U,V,T,Q,PS}
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
!        Clock_Hwin_lat0:     Specify the horizontal center of the window in degrees. 
!        Clock_Hwin_lon0:     The longitude must be in the range [0,360] and the 
!                             latitude should be [-90,+90].
!        Clock_Hwin_latWidth: Specify the lat and lon widths of the window as positive 
!        Clock_Hwin_lonWidth: values in degrees.Setting a width to a large value (e.g. 999) 
!                             renders the window a constant in that direction.
!        Clock_Hwin_latDelta: Controls the sharpness of the window transition with a 
!        Clock_Hwin_lonDelta: length in degrees. Small non-zero values yeild a step 
!                             function while a large value yeilds a smoother transition.
!        Clock_Hwin_Invert  : A logical flag used to invert the horizontal window function 
!                             to get its compliment.(e.g. to nudge outside a given window).
!
!        Clock_Vwin_Lindex:   In the vertical, the window is specified in terms of model 
!        Clock_Vwin_Ldelta:   level indcies. The High and Low transition levels should 
!        Clock_Vwin_Hindex:   range from [0,(NLEV+1)]. The transition lengths are also 
!        Clock_Vwin_Hdelta:   specified in terms of model indices. For a window function 
!                             constant in the vertical, the Low index should be set to 0,
!                             the High index should be set to (NLEV+1), and the transition 
!                             lengths should be set to 0.001 
!        Clock_Vwin_Invert  : A logical flag used to invert the vertical window function 
!                             to get its compliment.
!
!        EXAMPLE: For a channel window function centered at the equator and independent 
!                 of the vertical (30 levels):
!                        Clock_Hwin_lat0     = 0.         Clock_Vwin_Lindex = 0.
!                        Clock_Hwin_latWidth = 30.        Clock_Vwin_Ldelta = 0.001
!                        Clock_Hwin_latDelta = 5.0        Clock_Vwin_Hindex = 31.
!                        Clock_Hwin_lon0     = 180.       Clock_Vwin_Hdelta = 0.001 
!                        Clock_Hwin_lonWidth = 999.       Clock_Vwin_Invert = .false.
!                        Clock_Hwin_lonDelta = 1.0
!                        Clock_Hwin_Invert   = .false.
!
!                 If on the other hand one wanted to apply nudging at the poles and
!                 not at the equator, the settings would be similar but with:
!                        Clock_Hwin_Invert = .true.
!
!    A user can preview the window resulting from a given set of namelist values before 
!    running the model. Lookat_ClockWindow.ncl is a script avalable in the tools directory 
!    which will read in the values for a given namelist and display the resulting window.
!
!    The module is currently configured for only 1 window function. It can readily be 
!    extended for multiple windows if the need arises.
!
!
! Input/Output Values:
!    Forcing contributions are available for history file output by 
!    the names:    {'Clock_U','Clock_V','Clock_T',and 'Clock_Q'}
!    The target values that the model state is nudged toward are available for history 
!    file output via the variables:  {'Target_MU'   ,'Target_LAMBDAC','Target_ICSWP','Target_ICLWP',
!                                     'Target_ICIWP','Target_DES'    ,'Target_DEI'  ,'Target_CLD'  ,
!                                     'Target_CLDFSNOW', 'Target_Q'   }
!
!    &clocking_nl
!      Clock_Model         - LOGICAL toggle to activate cloud locking.
!                              TRUE  -> Nudging is on.
!                              FALSE -> Nudging is off.                            [DEFAULT]
!
!      Clock_Path          - CHAR path to the analyses files.
!                              (e.g. '/glade/scratch/USER/inputdata/ERAI-Data/')
!
!      Clock_File_Template - CHAR Analyses filename with year, month, day, and second
!                                 values replaced by %y, %m, %d, and %s respectively.
!                              (e.g. '%y/ERAI_ne30np4_L30.cam2.i.%y-%m-%d-%s.nc')
!
!      Clock_Times_Per_Day - INT Number of analyses files available per day.
!                              1 --> daily analyses.
!                              4 --> 6 hourly analyses.
!                              8 --> 3 hourly.
!
!      Clock_Model_Times_Per_Day - INT Number of times to update the model state (used for clocking) 
!                                each day. The value is restricted to be longer than the 
!                                current model timestep and shorter than the analyses 
!                                timestep. As this number is increased, the nudging
!                                force has the form of newtonian cooling.
!                              48 --> 1800 Second timestep.
!                              96 -->  900 Second timestep.
!
!      Clock_Beg_Year      - INT clocking begining year.  [1979- ]
!      Clock_Beg_Month     - INT clocking begining month. [1-12]
!      Clock_Beg_Day       - INT clocking begining day.   [1-31]
!      Clock_End_Year      - INT clocking ending year.    [1979-]
!      Clock_End_Month     - INT clocking ending month.   [1-12]
!      Clock_End_Day       - INT clocking ending day.     [1-31]
!
!      Clock_Force_Opt     - INT Index to select the clocking Target for a relaxation 
!                                forcing of the form: 
!                                where (t'==Analysis times ; t==Model Times)
!
!                              0 -> NEXT-OBS: Target=Anal(t'_next)                 [DEFAULT]
!                              1 -> LINEAR:   Target=(F*Anal(t'_curr) +(1-F)*Anal(t'_next))
!                                                 F =(t'_next - t_curr )/Tdlt_Anal
!
!      Clock_TimeScale_Opt - INT Index to select the timescale for clocking.
!                                where (t'==Analysis times ; t==Model Times) 
!
!                              0 -->  TimeScale = 1/Tdlt_Anal                      [DEFAULT]
!                              1 -->  TimeScale = 1/(t'_next - t_curr )
!
!      Clock_MUprof        - INT index of profile structure to use for MU.       [0,1,2]
!      Clock_LAMBDACprof   - INT index of profile structure to use for LAMBDAC.  [0,1,2]
!      Clock_ICSWPprof     - INT index of profile structure to use for ICSWP.    [0,1,2]
!      Clock_ICLWPprof     - INT index of profile structure to use for ICLWP.    [0,1,2]
!      Clock_ICIWPprof     - INT index of profile structure to use for ICIWP.    [0,1,2]
!      Clock_DESprof       - INT index of profile structure to use for DES.      [0,1,2]
!      Clock_DEIprof       - INT index of profile structure to use for DEI.      [0,1,2]
!      Clock_CLDprof       - INT index of profile structure to use for CLD.      [0,1,2]
!      Clock_CLDFSNOWprof  - INT index of profile structure to use for CLDFSNOW. [0,1,2]
!      Clock_Qprof         - INT index of profile structure to use for Q.        [0,1,2]
!      Clock_PSprof        - INT index of profile structure to use for PS.       [0,N/A]
!
!                                The spatial distribution is specified with a profile index.
!                                 Where:  0 == OFF      (No Nudging of this variable)
!                                         1 == CONSTANT (Spatially Uniform Nudging)
!                                         2 == HEAVISIDE WINDOW FUNCTION
!
!      Clock_MUcoef        - REAL fractional clocking coeffcient for MU. 
!      Clock_LAMBDACcoef   - REAL fractional clocking coeffcient for LAMBDAC. 
!      Clock_ICSWPcoef     - REAL fractional clocking coeffcient for ICSWP. 
!      Clock_ICLWPcoef     - REAL fractional clocking coeffcient for ICLWP. 
!      Clock_ICIWPcoef     - REAL fractional clocking coeffcient for ICIWP. 
!      Clock_DEScoef       - REAL fractional clocking coeffcient for DES. 
!      Clock_DEIcoef       - REAL fractional clocking coeffcient for DEI. 
!      Clock_CLDcoef       - REAL fractional clocking coeffcient for CLD. 
!      Clock_CLDFSNOWcoef  - REAL fractional clocking coeffcient for CLDFSNOW. 
!      Clock_Qcoef         - REAL fractional clocking coeffcient for Q. 
!      Clock_PScoef        - REAL fractional clocking coeffcient for PS. 
!
!                                 The strength of the clocking is specified as a fractional 
!                                 coeffcient between [0,1].
!           
!      Clock_Hwin_lat0     - REAL latitudinal center of window in degrees.
!      Clock_Hwin_lon0     - REAL longitudinal center of window in degrees.
!      Clock_Hwin_latWidth - REAL latitudinal width of window in degrees.
!      Clock_Hwin_lonWidth - REAL longitudinal width of window in degrees.
!      Clock_Hwin_latDelta - REAL latitudinal transition length of window in degrees.
!      Clock_Hwin_lonDelta - REAL longitudinal transition length of window in degrees.
!      Clock_Hwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
!                                    TRUE = value=0 inside the specified window, 1 outside
!      Clock_Vwin_Lindex   - REAL LO model index of transition
!      Clock_Vwin_Hindex   - REAL HI model index of transition
!      Clock_Vwin_Ldelta   - REAL LO transition length 
!      Clock_Vwin_Hdelta   - REAL HI transition length 
!      Clock_Vwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
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
  use shr_kind_mod,   only:r8=>SHR_KIND_R8,cs=>SHR_KIND_CS,cl=>SHR_KIND_CL
  use time_manager,   only:timemgr_time_ge,timemgr_time_inc,get_curr_date,get_step_size
  use phys_grid   ,   only:scatter_field_to_chunk
  use cam_abortutils, only:endrun
  use spmd_utils  ,   only:masterproc
  use cam_logfile ,   only:iulog
#ifdef SPMD
  use mpishorthand
#endif

  ! Set all Global values and routines to private by default 
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public:: Clock_Model,Clock_ON
  public:: clocking_readnl
  public:: clocking_init
  public:: clocking_timestep_init
  public:: clocking_timestep_tend
  private::clocking_update_analyses_se
  private::clocking_update_analyses_fv
  private::clocking_set_PSprofile
  private::clocking_set_profile

  ! Nudging Parameters
  !--------------------
  logical          :: Clock_Model       =.false.
  logical          :: Clock_ON          =.false.
  logical          :: Clock_Initialized =.false.
  character(len=cl):: Clock_Path
  character(len=cs):: Clock_File,Clock_File_Template
  integer          :: Clock_Force_Opt
  integer          :: Clock_TimeScale_Opt
  integer          :: Clock_Times_Per_Day
  integer          :: Clock_Model_Times_Per_Day
  real(r8)         :: Clock_MUcoef   ,Clock_LAMBDACcoef
  integer          :: Clock_MUprof   ,Clock_LAMBDACprof
  real(r8)         :: Clock_Qcoef    ,Clock_ICSWPcoef
  integer          :: Clock_Qprof    ,Clock_ICSWPprof
  real(r8)         :: Clock_ICLWPcoef,Clock_ICIWPcoef
  integer          :: Clock_ICLWPprof,Clock_ICIWPprof
  real(r8)         :: Clock_DEScoef  ,Clock_DEIcoef  
  integer          :: Clock_DESprof  ,Clock_DEIprof  
  real(r8)         :: Clock_CLDcoef  ,Clock_CLDFSNOWcoef  
  integer          :: Clock_CLDprof  ,Clock_CLDFSNOWprof  
  real(r8)         :: Clock_PScoef
  integer          :: Clock_PSprof
  integer          :: Clock_Beg_Year ,Clock_Beg_Month
  integer          :: Clock_Beg_Day  ,Clock_Beg_Sec
  integer          :: Clock_End_Year ,Clock_End_Month
  integer          :: Clock_End_Day  ,Clock_End_Sec
  integer          :: Clock_Curr_Year,Clock_Curr_Month
  integer          :: Clock_Curr_Day ,Clock_Curr_Sec
  integer          :: Clock_Next_Year,Clock_Next_Month
  integer          :: Clock_Next_Day ,Clock_Next_Sec
  integer          :: Clock_Step
  integer          :: Model_Curr_Year,Model_Curr_Month
  integer          :: Model_Curr_Day ,Model_Curr_Sec
  integer          :: Model_Next_Year,Model_Next_Month
  integer          :: Model_Next_Day ,Model_Next_Sec
  integer          :: Model_Step
  real(r8)         :: Clock_Hwin_lat0
  real(r8)         :: Clock_Hwin_latWidth
  real(r8)         :: Clock_Hwin_latDelta
  real(r8)         :: Clock_Hwin_lon0
  real(r8)         :: Clock_Hwin_lonWidth
  real(r8)         :: Clock_Hwin_lonDelta
  logical          :: Clock_Hwin_Invert = .false.
  real(r8)         :: Clock_Hwin_lo
  real(r8)         :: Clock_Hwin_hi
  real(r8)         :: Clock_Vwin_Hindex
  real(r8)         :: Clock_Vwin_Hdelta
  real(r8)         :: Clock_Vwin_Lindex
  real(r8)         :: Clock_Vwin_Ldelta
  logical          :: Clock_Vwin_Invert =.false.
  real(r8)         :: Clock_Vwin_lo
  real(r8)         :: Clock_Vwin_hi
  real(r8)         :: Clock_Hwin_latWidthH
  real(r8)         :: Clock_Hwin_lonWidthH
  real(r8)         :: Clock_Hwin_max
  real(r8)         :: Clock_Hwin_min

  ! Nudging State Arrays
  !-----------------------
  integer Clock_nlon,Clock_nlat,Clock_ncol,Clock_nlev
  real(r8),allocatable::Target_MU      (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_LAMBDAC (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_ICSWP   (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_ICLWP   (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_ICIWP   (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_DES     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_DEI     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_CLD     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_CLDFSNOW(:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_Q       (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_PS      (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable::Model_MU       (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_LAMBDAC  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_ICSWP    (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_ICLWP    (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_ICIWP    (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_DES      (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_DEI      (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_CLD      (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_CLDFSNOW (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_Q        (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_PS       (:,:)    !(pcols,begchunk:endchunk)

  real(r8),allocatable:: Clock_MUtau      (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_LAMBDACtau (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_ICSWPtau   (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_ICLWPtau   (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_ICIWPtau   (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_DEStau     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_DEItau     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_CLDtau     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_CLDFSNOWtau(:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_Qtau       (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_PStau      (:,:)    !(pcols,begchunk:endchunk)

  real(r8),allocatable:: Clock_MUstep      (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_LAMBDACstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_ICSWPstep   (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_ICLWPstep   (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_ICIWPstep   (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_DESstep     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_DEIstep     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_CLDstep     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_CLDFSNOWstep(:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_Qstep       (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Clock_PSstep      (:,:)    !(pcols,begchunk:endchunk)

  ! Nudging Observation Arrays
  !-----------------------------
  integer               Clock_NumObs
  integer,allocatable:: Clock_ObsInd(:)
  logical ,allocatable::Clock_File_Present(:)
  real(r8),allocatable::Nobs_MU      (:,:,:,:) !(pcols,pver,begchunk:endchunk,Clock_NumObs)
  real(r8),allocatable::Nobs_LAMBDAC (:,:,:,:) !(pcols,pver,begchunk:endchunk,Clock_NumObs)
  real(r8),allocatable::Nobs_ICSWP   (:,:,:,:) !(pcols,pver,begchunk:endchunk,Clock_NumObs)
  real(r8),allocatable::Nobs_ICLWP   (:,:,:,:) !(pcols,pver,begchunk:endchunk,Clock_NumObs)
  real(r8),allocatable::Nobs_ICIWP   (:,:,:,:) !(pcols,pver,begchunk:endchunk,Clock_NumObs)
  real(r8),allocatable::Nobs_DES     (:,:,:,:) !(pcols,pver,begchunk:endchunk,Clock_NumObs)
  real(r8),allocatable::Nobs_DEI     (:,:,:,:) !(pcols,pver,begchunk:endchunk,Clock_NumObs)
  real(r8),allocatable::Nobs_CLD     (:,:,:,:) !(pcols,pver,begchunk:endchunk,Clock_NumObs)
  real(r8),allocatable::Nobs_CLDFSNOW(:,:,:,:) !(pcols,pver,begchunk:endchunk,Clock_NumObs)
  real(r8),allocatable::Nobs_Q       (:,:,:,:) !(pcols,pver,begchunk:endchunk,Clock_NumObs)
  real(r8),allocatable::Nobs_PS      (:,:,:)   !(pcols,begchunk:endchunk,Clock_NumObs)

contains
  !================================================================
  subroutine clocking_readnl(nlfile)
   ! 
   ! CLOCKING_READNL: Initialize default values controlling the Nudging 
   !                 process. Then read namelist values to override 
   !                 them.
   !===============================================================
   use ppgrid        ,only: pver
   use namelist_utils,only:find_group_name
   use units         ,only:getunit,freeunit
   !
   ! Arguments
   !-------------
   character(len=*),intent(in)::nlfile
   !
   ! Local Values
   !---------------
   integer ierr,unitn

   namelist /clocking_nl/ Clock_Model,Clock_Path,                       &
                          Clock_File_Template,Clock_Force_Opt,          &
                          Clock_TimeScale_Opt,                          &
                          Clock_Times_Per_Day,Clock_Model_Times_Per_Day,&
                          Clock_MUcoef      ,Clock_MUprof,              &
                          Clock_LAMBDACcoef ,Clock_LAMBDACprof,         &
                          Clock_ICSWPcoef   ,Clock_ICSWPprof,           &
                          Clock_ICLWPcoef   ,Clock_ICLWPprof,           &
                          Clock_ICIWPcoef   ,Clock_ICIWPprof,           &
                          Clock_DEScoef     ,Clock_DESprof  ,           &
                          Clock_DEIcoef     ,Clock_DEIprof  ,           &
                          Clock_CLDcoef     ,Clock_CLDprof  ,           &
                          Clock_CLDFSNOWcoef,Clock_CLDFSNOWprof  ,      &
                          Clock_Qcoef ,Clock_Qprof,                     &
                          Clock_PScoef,Clock_PSprof,                    &
                          Clock_Beg_Year,Clock_Beg_Month,Clock_Beg_Day, &
                          Clock_End_Year,Clock_End_Month,Clock_End_Day, &
                          Clock_Hwin_lat0,Clock_Hwin_lon0,              &
                          Clock_Hwin_latWidth,Clock_Hwin_lonWidth,      &
                          Clock_Hwin_latDelta,Clock_Hwin_lonDelta,      &
                          Clock_Hwin_Invert,                            &
                          Clock_Vwin_Lindex,Clock_Vwin_Hindex,          &
                          Clock_Vwin_Ldelta,Clock_Vwin_Hdelta,          &
                          Clock_Vwin_Invert                            

   ! Nudging is NOT initialized yet, For now
   ! Nudging will always begin/end at midnight.
   !--------------------------------------------
   Clock_Initialized =.false.
   Clock_ON          =.false.
   Clock_Beg_Sec=0
   Clock_End_Sec=0

   ! Set Default Namelist values
   !-----------------------------
   Clock_Model         = .false.
   Clock_Path          = './Data/YOTC_ne30np4_001/'
   Clock_File_Template = 'YOTC_ne30np4_L30.cam2.i.%y-%m-%d-%s.nc'
   Clock_Force_Opt     = 0
   Clock_TimeScale_Opt = 0
   Clock_Times_Per_Day = 4
   Clock_Model_Times_Per_Day = 4
   Clock_MUcoef        = 0._r8
   Clock_LAMBDACcoef   = 0._r8
   Clock_ICSWPcoef     = 0._r8
   Clock_ICLWPcoef     = 0._r8
   Clock_ICIWPcoef     = 0._r8
   Clock_DEScoef       = 0._r8
   Clock_DEIcoef       = 0._r8
   Clock_CLDcoef       = 0._r8
   Clock_CLDFSNOWcoef  = 0._r8
   Clock_Qcoef         = 0._r8
   Clock_PScoef        = 0._r8
   Clock_MUprof        = 0
   Clock_LAMBDACprof   = 0
   Clock_ICSWPprof     = 0
   Clock_ICLWPprof     = 0
   Clock_ICIWPprof     = 0
   Clock_DESprof       = 0
   Clock_DEIprof       = 0
   Clock_CLDprof       = 0
   Clock_CLDFSNOWprof  = 0
   Clock_Qprof         = 0
   Clock_PSprof        = 0
   Clock_Beg_Year      = 2008
   Clock_Beg_Month     = 5
   Clock_Beg_Day       = 1
   Clock_End_Year      = 2008
   Clock_End_Month     = 9
   Clock_End_Day       = 1
   Clock_Hwin_lat0     = 0._r8
   Clock_Hwin_latWidth = 9999._r8
   Clock_Hwin_latDelta = 1.0_r8
   Clock_Hwin_lon0     = 180._r8
   Clock_Hwin_lonWidth = 9999._r8
   Clock_Hwin_lonDelta = 1.0_r8
   Clock_Hwin_Invert   = .false.
   Clock_Hwin_lo       = 0.0_r8
   Clock_Hwin_hi       = 1.0_r8
   Clock_Vwin_Hindex   = float(pver+1)
   Clock_Vwin_Hdelta   = 0.001_r8
   Clock_Vwin_Lindex   = 0.0_r8
   Clock_Vwin_Ldelta   = 0.001_r8
   Clock_Vwin_Invert   = .false.
   Clock_Vwin_lo       = 0.0_r8
   Clock_Vwin_hi       = 1.0_r8

   ! Read in namelist values
   !------------------------
   if(masterproc) then
     unitn = getunit()
     open(unitn,file=trim(nlfile),status='old')
     call find_group_name(unitn,'clocking_nl',status=ierr)
     if(ierr.eq.0) then
       read(unitn,clocking_nl,iostat=ierr)
       if(ierr.ne.0) then
         call endrun('clocking_readnl:: ERROR reading namelist')
       endif
     endif
     close(unitn)
     call freeunit(unitn)
   endif

   ! Set hi/lo values according to the given '_Invert' parameters
   !--------------------------------------------------------------
   if(Clock_Hwin_Invert) then
     Clock_Hwin_lo = 1.0_r8
     Clock_Hwin_hi = 0.0_r8
   else
     Clock_Hwin_lo = 0.0_r8
     Clock_Hwin_hi = 1.0_r8
   endif

   if(Clock_Vwin_Invert) then
     Clock_Vwin_lo = 1.0_r8
     Clock_Vwin_hi = 0.0_r8
   else
     Clock_Vwin_lo = 0.0_r8
     Clock_Vwin_hi = 1.0_r8
   endif

   ! Check for valid namelist values 
   !----------------------------------
   if((Clock_Hwin_lat0.lt.-90._r8).or.(Clock_Hwin_lat0.gt.+90._r8)) then
     write(iulog,*) 'CLOCKING: Window lat0 must be in [-90,+90]'
     write(iulog,*) 'CLOCKING:  Clock_Hwin_lat0=',Clock_Hwin_lat0
     call endrun('clocking_readnl:: ERROR in namelist')
   endif

   if((Clock_Hwin_lon0.lt.0._r8).or.(Clock_Hwin_lon0.ge.360._r8)) then
     write(iulog,*) 'CLOCKING: Window lon0 must be in [0,+360)'
     write(iulog,*) 'CLOCKING:  Clock_Hwin_lon0=',Clock_Hwin_lon0
     call endrun('clocking_readnl:: ERROR in namelist')
   endif

   if((Clock_Vwin_Lindex.gt.Clock_Vwin_Hindex)                         .or. &
      (Clock_Vwin_Hindex.gt.float(pver+1)).or.(Clock_Vwin_Hindex.lt.0._r8).or. &
      (Clock_Vwin_Lindex.gt.float(pver+1)).or.(Clock_Vwin_Lindex.lt.0._r8)   ) then
     write(iulog,*) 'CLOCKING: Window Lindex must be in [0,pver+1]'
     write(iulog,*) 'CLOCKING: Window Hindex must be in [0,pver+1]'
     write(iulog,*) 'CLOCKING: Lindex must be LE than Hindex'
     write(iulog,*) 'CLOCKING:  Clock_Vwin_Lindex=',Clock_Vwin_Lindex
     write(iulog,*) 'CLOCKING:  Clock_Vwin_Hindex=',Clock_Vwin_Hindex
     call endrun('clocking_readnl:: ERROR in namelist')
   endif

   if((Clock_Hwin_latDelta.le.0._r8).or.(Clock_Hwin_lonDelta.le.0._r8).or. &
      (Clock_Vwin_Hdelta  .le.0._r8).or.(Clock_Vwin_Ldelta  .le.0._r8)    ) then
     write(iulog,*) 'CLOCKING: Window Deltas must be positive'
     write(iulog,*) 'CLOCKING:  Clock_Hwin_latDelta=',Clock_Hwin_latDelta
     write(iulog,*) 'CLOCKING:  Clock_Hwin_lonDelta=',Clock_Hwin_lonDelta
     write(iulog,*) 'CLOCKING:  Clock_Vwin_Hdelta=',Clock_Vwin_Hdelta
     write(iulog,*) 'CLOCKING:  Clock_Vwin_Ldelta=',Clock_Vwin_Ldelta
     call endrun('clocking_readnl:: ERROR in namelist')

   endif

   if((Clock_Hwin_latWidth.le.0._r8).or.(Clock_Hwin_lonWidth.le.0._r8)) then
     write(iulog,*) 'CLOCKING: Window widths must be positive'
     write(iulog,*) 'CLOCKING:  Clock_Hwin_latWidth=',Clock_Hwin_latWidth
     write(iulog,*) 'CLOCKING:  Clock_Hwin_lonWidth=',Clock_Hwin_lonWidth
     call endrun('clocking_readnl:: ERROR in namelist')
   endif

   ! Broadcast namelist variables
   !------------------------------
#ifdef SPMD
   call mpibcast(Clock_Path         ,len(Clock_Path)         ,mpichar,0,mpicom)
   call mpibcast(Clock_File_Template,len(Clock_File_Template),mpichar,0,mpicom)
   call mpibcast(Clock_Model        , 1, mpilog, 0, mpicom)
   call mpibcast(Clock_Initialized  , 1, mpilog, 0, mpicom)
   call mpibcast(Clock_ON           , 1, mpilog, 0, mpicom)
   call mpibcast(Clock_Force_Opt    , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_TimeScale_Opt, 1, mpiint, 0, mpicom)
   call mpibcast(Clock_Times_Per_Day, 1, mpiint, 0, mpicom)
   call mpibcast(Clock_Model_Times_Per_Day, 1, mpiint, 0, mpicom)
   call mpibcast(Clock_MUcoef       , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_LAMBDACcoef  , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_ICSWPcoef    , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_ICLWPcoef    , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_ICIWPcoef    , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_DEScoef      , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_DEIcoef      , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_CLDcoef      , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_CLDFSNOWcoef , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_Qcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_PScoef       , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_MUprof       , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_LAMBDACprof  , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_ICSWPprof    , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_ICLWPprof    , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_ICIWPprof    , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_DESprof      , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_DEIprof      , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_CLDprof      , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_CLDFSNOWprof , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_Qprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_PSprof       , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_Beg_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_Beg_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_Beg_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_Beg_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_End_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_End_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_End_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_End_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(Clock_Hwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_Hwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_Hwin_lat0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_Hwin_latWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_Hwin_latDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_Hwin_lon0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_Hwin_lonWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_Hwin_lonDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_Hwin_Invert,   1, mpilog, 0, mpicom)
   call mpibcast(Clock_Vwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_Vwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_Vwin_Hindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_Vwin_Hdelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_Vwin_Lindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_Vwin_Ldelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(Clock_Vwin_Invert,   1, mpilog, 0, mpicom)
#endif

   ! End Routine
   !------------
   return
  end subroutine ! clocking_readnl
  !================================================================


  !================================================================
  subroutine clocking_init
   ! 
   ! CLOCKING_INIT: Allocate space and initialize Nudging values
   !===============================================================
   use ppgrid        ,only: pver,pcols,begchunk,endchunk
   use error_messages,only: alloc_err
   use dycore        ,only: dycore_is
   use dyn_grid      ,only: get_horiz_grid_dim_d
   use phys_grid     ,only: get_rlat_p,get_rlon_p,get_ncols_p
   use cam_history   ,only: addfld
   use shr_const_mod ,only: SHR_CONST_PI
   use filenames     ,only: interpret_filename_spec

   ! Local values
   !----------------
   integer  Year,Month,Day,Sec
   integer  YMD1,YMD
   logical  After_Beg,Before_End
   integer  istat,lchnk,ncol,icol,ilev
   integer  hdim1_d,hdim2_d
   integer  dtime
   real(r8) rlat,rlon
   real(r8) Wprof(pver)
   real(r8) lonp,lon0,lonn,latp,lat0,latn
   real(r8) Val1_p,Val2_p,Val3_p,Val4_p
   real(r8) Val1_0,Val2_0,Val3_0,Val4_0
   real(r8) Val1_n,Val2_n,Val3_n,Val4_n
   integer               nn

   ! Get the time step size
   !------------------------
   dtime = get_step_size()

   ! Allocate Space for Nudging data arrays
   !-----------------------------------------
   allocate(Target_MU(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Target_MU',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_LAMBDAC(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Target_LAMBDAC',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_ICSWP(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Target_ICSWP',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_ICLWP(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Target_ICLWP',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_ICIWP(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Target_ICIWP',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_DES(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Target_DES',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_DEI(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Target_DEI',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_CLD(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Target_CLD',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_CLDFSNOW(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Target_CLDFSNOW',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_Q(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Target_Q',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_PS(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Target_PS',pcols*((endchunk-begchunk)+1))

   allocate(Model_MU(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Model_MU',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_LAMBDAC(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Model_LAMBDAC',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_ICSWP(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Model_ICSWP',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_ICLWP(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Model_ICLWP',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_ICIWP(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Model_ICIWP',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_DES(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Model_DES',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_DEI(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Model_DEI',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_CLD(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Model_CLD',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_CLDFSNOW(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Model_CLDFSNOW',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_Q(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Model_Q',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_PS(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Model_PS',pcols*((endchunk-begchunk)+1))

   ! Allocate Space for spatial dependence of 
   ! Nudging Coefs and Nudging Forcing.
   !-------------------------------------------
   allocate(Clock_MUtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_MUtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_LAMBDACtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_LAMBDACtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_ICSWPtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_ICSWPtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_ICLWPtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_ICLWPtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_ICIWPtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_ICIWPtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_DEStau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_DEStau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_DEItau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_DEItau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_CLDtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_CLDtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_CLDFSNOWtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_CLDFSNOWtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_Qtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_Qtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_PStau(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_PStau',pcols*((endchunk-begchunk)+1))

   allocate(Clock_MUstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_MUstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_LAMBDACstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_LAMBDACstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_ICSWPstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_ICSWPstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_ICLWPstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_ICLWPstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_ICIWPstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_ICIWPstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_DESstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_DESstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_DEIstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_DEIstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_CLDstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_CLDstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_CLDFSNOWstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_CLDFSNOWstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_Qstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_Qstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Clock_PSstep(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'clocking_init','Clock_PSstep',pcols*((endchunk-begchunk)+1))

   ! Register output fields with the cam history module
   !-----------------------------------------------------
!IMPLEMENT?   call addfld( 'Nudge_U'       ,(/ 'lev' /),'A','m/s/s'  ,'U Nudging Tendency')
!IMPLEMENT?   call addfld( 'Nudge_V'       ,(/ 'lev' /),'A','m/s/s'  ,'V Nudging Tendency')
!IMPLEMENT?   call addfld( 'Nudge_T'       ,(/ 'lev' /),'A','K/s'    ,'T Nudging Tendency')
!IMPLEMENT?   call addfld( 'Nudge_Q'       ,(/ 'lev' /),'A','kg/kg/s','Q Nudging Tendency')
   call addfld('Target_MU'      ,(/ 'lev' /),'A','m/s'    ,'MU        Target'  )
   call addfld('Target_LAMBDAC' ,(/ 'lev' /),'A','m/s'    ,'LAMBDAC   Target'  )
   call addfld('Target_ICSWP'   ,(/ 'lev' /),'A','K'      ,'ICSWP     Target'  )
   call addfld('Target_ICLWP'   ,(/ 'lev' /),'A','K'      ,'ICLWP     Target'  )
   call addfld('Target_ICIWP'   ,(/ 'lev' /),'A','K'      ,'ICIWP     Target'  )
   call addfld('Target_DES'     ,(/ 'lev' /),'A','K'      ,'DES       Target'  )
   call addfld('Target_DEI'     ,(/ 'lev' /),'A','K'      ,'DEI       Target'  )
   call addfld('Target_CLD'     ,(/ 'lev' /),'A','K'      ,'CLD       Target'  )
   call addfld('Target_CLDFSNOW',(/ 'lev' /),'A','K'      ,'CLDFSNOW  Target'  )
   call addfld('Target_Q'       ,(/ 'lev' /),'A','kg/kg'  ,'Q         Target  ')

   !-----------------------------------------
   ! Values initialized only by masterproc
   !-----------------------------------------
   if(masterproc) then

     ! Set the Stepping intervals for Model and Nudging values
     ! Ensure that the Model_Step is not smaller then one timestep
     !  and not larger then the Clock_Step.
     !--------------------------------------------------------
     Model_Step=86400/Clock_Model_Times_Per_Day
     Clock_Step=86400/Clock_Times_Per_Day
     if(Model_Step.lt.dtime) then
       write(iulog,*) ' '
       write(iulog,*) 'CLOCKING: Model_Step cannot be less than a model timestep'
       write(iulog,*) 'CLOCKING:  Setting Model_Step=dtime , dtime=',dtime
       write(iulog,*) ' '
       Model_Step=dtime
     endif
     if(Model_Step.gt.Clock_Step) then
       write(iulog,*) ' '
       write(iulog,*) 'CLOCKING: Model_Step cannot be more than Clock_Step'
       write(iulog,*) 'CLOCKING:  Setting Model_Step=Clock_Step, Clock_Step=',Clock_Step
       write(iulog,*) ' '
       Model_Step=Clock_Step
     endif

     ! Initialize column and level dimensions
     !--------------------------------------------------------
     call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
     Clock_nlon=hdim1_d
     Clock_nlat=hdim2_d
     Clock_ncol=hdim1_d*hdim2_d
     Clock_nlev=pver

     ! Check the time relative to the clocking window
     !------------------------------------------------
     call get_curr_date(Year,Month,Day,Sec)
     YMD=(Year*10000) + (Month*100) + Day
     YMD1=(Clock_Beg_Year*10000) + (Clock_Beg_Month*100) + Clock_Beg_Day
     call timemgr_time_ge(YMD1,Clock_Beg_Sec,         &
                          YMD ,Sec          ,After_Beg)
     YMD1=(Clock_End_Year*10000) + (Clock_End_Month*100) + Clock_End_Day
     call timemgr_time_ge(YMD ,Sec          ,          &
                          YMD1,Clock_End_Sec,Before_End)
  
     if((After_Beg).and.(Before_End)) then
       ! Set Time indicies so that the next call to 
       ! timestep_init will initialize the data arrays.
       !--------------------------------------------
       Model_Next_Year =Year
       Model_Next_Month=Month
       Model_Next_Day  =Day
       Model_Next_Sec  =(Sec/Model_Step)*Model_Step
       Clock_Next_Year =Year
       Clock_Next_Month=Month
       Clock_Next_Day  =Day
       Clock_Next_Sec  =(Sec/Clock_Step)*Clock_Step
     elseif(.not.After_Beg) then
       ! Set Time indicies to Nudging start,
       ! timestep_init will initialize the data arrays.
       !--------------------------------------------
       Model_Next_Year =Clock_Beg_Year
       Model_Next_Month=Clock_Beg_Month
       Model_Next_Day  =Clock_Beg_Day
       Model_Next_Sec  =Clock_Beg_Sec
       Clock_Next_Year =Clock_Beg_Year
       Clock_Next_Month=Clock_Beg_Month
       Clock_Next_Day  =Clock_Beg_Day
       Clock_Next_Sec  =Clock_Beg_Sec
     elseif(.not.Before_End) then
       ! Nudging will never occur, so switch it off
       !--------------------------------------------
       Clock_Model=.false.
       Clock_ON   =.false.
       write(iulog,*) ' '
       write(iulog,*) 'CLOCKING: WARNING - Cloud Locking has been requested by it will'
       write(iulog,*) 'CLOCKING:           never occur for the given time values'
       write(iulog,*) ' '
     endif

     ! Initialize values for window function  
     !----------------------------------------
     lonp= 180._r8
     lon0=   0._r8
     lonn=-180._r8
     latp=  90._r8-Clock_Hwin_lat0
     lat0=   0._r8
     latn= -90._r8-Clock_Hwin_lat0
    
     Clock_Hwin_lonWidthH=Clock_Hwin_lonWidth/2._r8
     Clock_Hwin_latWidthH=Clock_Hwin_latWidth/2._r8

     Val1_p=(1._r8+tanh((Clock_Hwin_lonWidthH+lonp)/Clock_Hwin_lonDelta))/2._r8
     Val2_p=(1._r8+tanh((Clock_Hwin_lonWidthH-lonp)/Clock_Hwin_lonDelta))/2._r8
     Val3_p=(1._r8+tanh((Clock_Hwin_latWidthH+latp)/Clock_Hwin_latDelta))/2._r8
     Val4_p=(1._r8+tanh((Clock_Hwin_latWidthH-latp)/Clock_Hwin_latDelta))/2_r8
     Val1_0=(1._r8+tanh((Clock_Hwin_lonWidthH+lon0)/Clock_Hwin_lonDelta))/2._r8
     Val2_0=(1._r8+tanh((Clock_Hwin_lonWidthH-lon0)/Clock_Hwin_lonDelta))/2._r8
     Val3_0=(1._r8+tanh((Clock_Hwin_latWidthH+lat0)/Clock_Hwin_latDelta))/2._r8
     Val4_0=(1._r8+tanh((Clock_Hwin_latWidthH-lat0)/Clock_Hwin_latDelta))/2._r8

     Val1_n=(1._r8+tanh((Clock_Hwin_lonWidthH+lonn)/Clock_Hwin_lonDelta))/2._r8
     Val2_n=(1._r8+tanh((Clock_Hwin_lonWidthH-lonn)/Clock_Hwin_lonDelta))/2._r8
     Val3_n=(1._r8+tanh((Clock_Hwin_latWidthH+latn)/Clock_Hwin_latDelta))/2._r8
     Val4_n=(1._r8+tanh((Clock_Hwin_latWidthH-latn)/Clock_Hwin_latDelta))/2._r8

     Clock_Hwin_max=     Val1_0*Val2_0*Val3_0*Val4_0
     Clock_Hwin_min=min((Val1_p*Val2_p*Val3_n*Val4_n), &
                        (Val1_p*Val2_p*Val3_p*Val4_p), &
                        (Val1_n*Val2_n*Val3_n*Val4_n), &
                        (Val1_n*Val2_n*Val3_p*Val4_p))

     ! Initialize number of clocking observation values to keep track of.
     ! Allocate and initialize observation indices 
     !-----------------------------------------------------------------
     if((Clock_Force_Opt.ge.0).and.(Clock_Force_Opt.le.1)) then
       Clock_NumObs=2
     else
       ! Additional Options may need OBS values at more times.
       !------------------------------------------------------
       Clock_NumObs=2
       write(iulog,*) 'CLOCKING: Setting Clock_NumObs=2'
       write(iulog,*) 'CLOCKING: WARNING: Unknown Clock_Force_Opt=',Clock_Force_Opt
       call endrun('CLOCKING: Unknown Forcing Option')
     endif
     allocate(Clock_ObsInd(Clock_NumObs),stat=istat)
     call alloc_err(istat,'clocking_init','Clock_ObsInd',Clock_NumObs)
     allocate(Clock_File_Present(Clock_NumObs),stat=istat)
     call alloc_err(istat,'clocking_init','Clock_File_Present',Clock_NumObs)
     do nn=1,Clock_NumObs
       Clock_ObsInd(nn) = Clock_NumObs+1-nn
     end do
     Clock_File_Present(:)=.false.

     ! Initialization is done, 
     !--------------------------
     Clock_Initialized=.true.

     ! Check that this is a valid DYCORE model
     !------------------------------------------
     if((.not.dycore_is('UNSTRUCTURED')).and. &
        (.not.dycore_is('EUL')         ).and. &
        (.not.dycore_is('LR')          )      ) then
       call endrun('CLOCKING IS CURRENTLY ONLY CONFIGURED FOR CAM-SE, FV, or EUL')
     endif

     ! Informational Output
     !---------------------------
     write(iulog,*) ' '
     write(iulog,*) '---------------------------------------------------------'
     write(iulog,*) '  MODEL CLOCKING INITIALIZED WITH THE FOLLOWING SETTINGS: '
     write(iulog,*) '---------------------------------------------------------'
     write(iulog,*) 'CLOCKING: Clock_Model=',Clock_Model
     write(iulog,*) 'CLOCKING: Clock_Path=',Clock_Path
     write(iulog,*) 'CLOCKING: Clock_File_Template =',Clock_File_Template
     write(iulog,*) 'CLOCKING: Clock_Force_Opt=',Clock_Force_Opt    
     write(iulog,*) 'CLOCKING: Clock_TimeScale_Opt=',Clock_TimeScale_Opt    
     write(iulog,*) 'CLOCKING: Clock_Times_Per_Day=',Clock_Times_Per_Day
     write(iulog,*) 'CLOCKING: Clock_Model_Times_Per_Day=',Clock_Model_Times_Per_Day
     write(iulog,*) 'CLOCKING: Clock_Step=',Clock_Step
     write(iulog,*) 'CLOCKING: Model_Step=',Model_Step
     write(iulog,*) 'CLOCKING: Clock_MUcoef  =',Clock_MUcoef
     write(iulog,*) 'CLOCKING: Clock_LAMBDACcoef  =',Clock_LAMBDACcoef
     write(iulog,*) 'CLOCKING: Clock_ICSWPcoef  =',Clock_ICSWPcoef
     write(iulog,*) 'CLOCKING: Clock_ICLWPcoef  =',Clock_ICLWPcoef
     write(iulog,*) 'CLOCKING: Clock_ICIWPcoef  =',Clock_ICIWPcoef
     write(iulog,*) 'CLOCKING: Clock_DEScoef  =',Clock_DEScoef
     write(iulog,*) 'CLOCKING: Clock_DEIcoef  =',Clock_DEIcoef
     write(iulog,*) 'CLOCKING: Clock_CLDcoef  =',Clock_CLDcoef
     write(iulog,*) 'CLOCKING: Clock_CLDFSNOWcoef  =',Clock_CLDFSNOWcoef
     write(iulog,*) 'CLOCKING: Clock_Qcoef  =',Clock_Qcoef
     write(iulog,*) 'CLOCKING: Clock_PScoef =',Clock_PScoef
     write(iulog,*) 'CLOCKING: Clock_MUprof  =',Clock_MUprof
     write(iulog,*) 'CLOCKING: Clock_LAMBDACprof  =',Clock_LAMBDACprof
     write(iulog,*) 'CLOCKING: Clock_ICSWPprof  =',Clock_ICSWPprof
     write(iulog,*) 'CLOCKING: Clock_ICLWPprof  =',Clock_ICLWPprof
     write(iulog,*) 'CLOCKING: Clock_ICIWPprof  =',Clock_ICIWPprof
     write(iulog,*) 'CLOCKING: Clock_DESprof  =',Clock_DESprof
     write(iulog,*) 'CLOCKING: Clock_DEIprof  =',Clock_DEIprof
     write(iulog,*) 'CLOCKING: Clock_CLDprof  =',Clock_CLDprof
     write(iulog,*) 'CLOCKING: Clock_CLDFSNOWprof  =',Clock_CLDFSNOWprof
     write(iulog,*) 'CLOCKING: Clock_Qprof  =',Clock_Qprof
     write(iulog,*) 'CLOCKING: Clock_PSprof =',Clock_PSprof
     write(iulog,*) 'CLOCKING: Clock_Beg_Year =',Clock_Beg_Year
     write(iulog,*) 'CLOCKING: Clock_Beg_Month=',Clock_Beg_Month
     write(iulog,*) 'CLOCKING: Clock_Beg_Day  =',Clock_Beg_Day
     write(iulog,*) 'CLOCKING: Clock_End_Year =',Clock_End_Year
     write(iulog,*) 'CLOCKING: Clock_End_Month=',Clock_End_Month
     write(iulog,*) 'CLOCKING: Clock_End_Day  =',Clock_End_Day
     write(iulog,*) 'CLOCKING: Clock_Hwin_lat0     =',Clock_Hwin_lat0
     write(iulog,*) 'CLOCKING: Clock_Hwin_latWidth =',Clock_Hwin_latWidth
     write(iulog,*) 'CLOCKING: Clock_Hwin_latDelta =',Clock_Hwin_latDelta
     write(iulog,*) 'CLOCKING: Clock_Hwin_lon0     =',Clock_Hwin_lon0
     write(iulog,*) 'CLOCKING: Clock_Hwin_lonWidth =',Clock_Hwin_lonWidth
     write(iulog,*) 'CLOCKING: Clock_Hwin_lonDelta =',Clock_Hwin_lonDelta
     write(iulog,*) 'CLOCKING: Clock_Hwin_Invert   =',Clock_Hwin_Invert  
     write(iulog,*) 'CLOCKING: Clock_Hwin_lo       =',Clock_Hwin_lo
     write(iulog,*) 'CLOCKING: Clock_Hwin_hi       =',Clock_Hwin_hi
     write(iulog,*) 'CLOCKING: Clock_Vwin_Hindex   =',Clock_Vwin_Hindex
     write(iulog,*) 'CLOCKING: Clock_Vwin_Hdelta   =',Clock_Vwin_Hdelta
     write(iulog,*) 'CLOCKING: Clock_Vwin_Lindex   =',Clock_Vwin_Lindex
     write(iulog,*) 'CLOCKING: Clock_Vwin_Ldelta   =',Clock_Vwin_Ldelta
     write(iulog,*) 'CLOCKING: Clock_Vwin_Invert   =',Clock_Vwin_Invert  
     write(iulog,*) 'CLOCKING: Clock_Vwin_lo       =',Clock_Vwin_lo
     write(iulog,*) 'CLOCKING: Clock_Vwin_hi       =',Clock_Vwin_hi
     write(iulog,*) 'CLOCKING: Clock_Hwin_latWidthH=',Clock_Hwin_latWidthH
     write(iulog,*) 'CLOCKING: Clock_Hwin_lonWidthH=',Clock_Hwin_lonWidthH
     write(iulog,*) 'CLOCKING: Clock_Hwin_max      =',Clock_Hwin_max
     write(iulog,*) 'CLOCKING: Clock_Hwin_min      =',Clock_Hwin_min
     write(iulog,*) 'CLOCKING: Clock_Initialized   =',Clock_Initialized
     write(iulog,*) ' '
     write(iulog,*) 'CLOCKING: Clock_NumObs=',Clock_NumObs
     write(iulog,*) ' '

   endif ! (masterproc) then

   ! Broadcast other variables that have changed
   !---------------------------------------------
#ifdef SPMD
   call mpibcast(Model_Step          ,            1, mpir8 , 0, mpicom)
   call mpibcast(Clock_Step          ,            1, mpir8 , 0, mpicom)
   call mpibcast(Model_Next_Year     ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Month    ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Day      ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Sec      ,            1, mpiint, 0, mpicom)
   call mpibcast(Clock_Next_Year     ,            1, mpiint, 0, mpicom)
   call mpibcast(Clock_Next_Month    ,            1, mpiint, 0, mpicom)
   call mpibcast(Clock_Next_Day      ,            1, mpiint, 0, mpicom)
   call mpibcast(Clock_Next_Sec      ,            1, mpiint, 0, mpicom)
   call mpibcast(Clock_Model         ,            1, mpilog, 0, mpicom)
   call mpibcast(Clock_ON            ,            1, mpilog, 0, mpicom)
   call mpibcast(Clock_Initialized   ,            1, mpilog, 0, mpicom)
   call mpibcast(Clock_ncol          ,            1, mpiint, 0, mpicom)
   call mpibcast(Clock_nlev          ,            1, mpiint, 0, mpicom)
   call mpibcast(Clock_nlon          ,            1, mpiint, 0, mpicom)
   call mpibcast(Clock_nlat          ,            1, mpiint, 0, mpicom)
   call mpibcast(Clock_Hwin_max      ,            1, mpir8 , 0, mpicom)
   call mpibcast(Clock_Hwin_min      ,            1, mpir8 , 0, mpicom)
   call mpibcast(Clock_Hwin_lonWidthH,            1, mpir8 , 0, mpicom)
   call mpibcast(Clock_Hwin_latWidthH,            1, mpir8 , 0, mpicom)
   call mpibcast(Clock_NumObs        ,            1, mpiint, 0, mpicom)
#endif

   ! All non-masterproc processes also need to allocate space
   ! before the broadcast of Clock_NumObs dependent data.
   !------------------------------------------------------------
   if(.not.masterproc) then
     allocate(Clock_ObsInd(Clock_NumObs),stat=istat)
     call alloc_err(istat,'clocking_init','Clock_ObsInd',Clock_NumObs)
     allocate(Clock_File_Present(Clock_NumObs),stat=istat)
     call alloc_err(istat,'clocking_init','Clock_File_Present',Clock_NumObs)
   endif
#ifdef SPMD
   call mpibcast(Clock_ObsInd        , Clock_NumObs, mpiint, 0, mpicom)
   call mpibcast(Clock_File_Present  , Clock_NumObs, mpilog, 0, mpicom)
#endif

   ! Allocate Space for Nudging observation arrays, initialize with 0's
   !---------------------------------------------------------------------
   allocate(Nobs_MU(pcols,pver,begchunk:endchunk,Clock_NumObs),stat=istat)
   call alloc_err(istat,'clocking_init','Nobs_MU',pcols*pver*((endchunk-begchunk)+1)*Clock_NumObs)
   allocate(Nobs_LAMBDAC(pcols,pver,begchunk:endchunk,Clock_NumObs),stat=istat)
   call alloc_err(istat,'clocking_init','Nobs_LAMBDAC',pcols*pver*((endchunk-begchunk)+1)*Clock_NumObs)
   allocate(Nobs_ICSWP(pcols,pver,begchunk:endchunk,Clock_NumObs),stat=istat)
   call alloc_err(istat,'clocking_init','Nobs_ICSWP',pcols*pver*((endchunk-begchunk)+1)*Clock_NumObs)
   allocate(Nobs_ICLWP(pcols,pver,begchunk:endchunk,Clock_NumObs),stat=istat)
   call alloc_err(istat,'clocking_init','Nobs_ICLWP',pcols*pver*((endchunk-begchunk)+1)*Clock_NumObs)
   allocate(Nobs_ICIWP(pcols,pver,begchunk:endchunk,Clock_NumObs),stat=istat)
   call alloc_err(istat,'clocking_init','Nobs_ICIWP',pcols*pver*((endchunk-begchunk)+1)*Clock_NumObs)
   allocate(Nobs_DES(pcols,pver,begchunk:endchunk,Clock_NumObs),stat=istat)
   call alloc_err(istat,'clocking_init','Nobs_DES',pcols*pver*((endchunk-begchunk)+1)*Clock_NumObs)
   allocate(Nobs_DEI(pcols,pver,begchunk:endchunk,Clock_NumObs),stat=istat)
   call alloc_err(istat,'clocking_init','Nobs_DEI',pcols*pver*((endchunk-begchunk)+1)*Clock_NumObs)
   allocate(Nobs_CLD(pcols,pver,begchunk:endchunk,Clock_NumObs),stat=istat)
   call alloc_err(istat,'clocking_init','Nobs_CLD',pcols*pver*((endchunk-begchunk)+1)*Clock_NumObs)
   allocate(Nobs_CLDFSNOW(pcols,pver,begchunk:endchunk,Clock_NumObs),stat=istat)
   call alloc_err(istat,'clocking_init','Nobs_CLDFSNOW',pcols*pver*((endchunk-begchunk)+1)*Clock_NumObs)
   allocate(Nobs_Q(pcols,pver,begchunk:endchunk,Clock_NumObs),stat=istat)
   call alloc_err(istat,'clocking_init','Nobs_Q',pcols*pver*((endchunk-begchunk)+1)*Clock_NumObs)
   allocate(Nobs_PS(pcols,begchunk:endchunk,Clock_NumObs),stat=istat)
   call alloc_err(istat,'clocking_init','Nobs_PS',pcols*((endchunk-begchunk)+1)*Clock_NumObs)

   Nobs_MU      (:pcols,:pver,begchunk:endchunk,:Clock_NumObs)=0._r8
   Nobs_LAMBDAC (:pcols,:pver,begchunk:endchunk,:Clock_NumObs)=0._r8
   Nobs_ICSWP   (:pcols,:pver,begchunk:endchunk,:Clock_NumObs)=0._r8
   Nobs_ICLWP   (:pcols,:pver,begchunk:endchunk,:Clock_NumObs)=0._r8
   Nobs_ICIWP   (:pcols,:pver,begchunk:endchunk,:Clock_NumObs)=0._r8
   Nobs_DES     (:pcols,:pver,begchunk:endchunk,:Clock_NumObs)=0._r8
   Nobs_DEI     (:pcols,:pver,begchunk:endchunk,:Clock_NumObs)=0._r8
   Nobs_CLD     (:pcols,:pver,begchunk:endchunk,:Clock_NumObs)=0._r8
   Nobs_CLDFSNOW(:pcols,:pver,begchunk:endchunk,:Clock_NumObs)=0._r8
   Nobs_Q       (:pcols,:pver,begchunk:endchunk,:Clock_NumObs)=0._r8
   Nobs_PS      (:pcols      ,begchunk:endchunk,:Clock_NumObs)=0._r8

!!DIAG
   if(masterproc) then
     write(iulog,*) 'CLOCKING: clocking_init() OBS arrays allocated and initialized'
     write(iulog,*) 'CLOCKING: clocking_init() SIZE#',(9*pcols*pver*((endchunk-begchunk)+1)*Clock_NumObs)
     write(iulog,*) 'CLOCKING: clocking_init() MB:',float(8*9*pcols*pver*((endchunk-begchunk)+1)*Clock_NumObs)/(1024._r8*1024._r8)
     write(iulog,*) 'CLOCKING: clocking_init() pcols=',pcols,' pver=',pver
     write(iulog,*) 'CLOCKING: clocking_init() begchunk:',begchunk,' endchunk=',endchunk
     write(iulog,*) 'CLOCKING: clocking_init() chunk:',(endchunk-begchunk+1),' Clock_NumObs=',Clock_NumObs
     write(iulog,*) 'CLOCKING: clocking_init() Clock_ObsInd=',Clock_ObsInd
     write(iulog,*) 'CLOCKING: clocking_init() Clock_File_Present=',Clock_File_Present
   endif
!!DIAG

   ! Initialize the analysis filename at the NEXT time for startup.
   !---------------------------------------------------------------
   Clock_File=interpret_filename_spec(Clock_File_Template      , &
                                       yr_spec=Clock_Next_Year , &
                                      mon_spec=Clock_Next_Month, &
                                      day_spec=Clock_Next_Day  , &
                                      sec_spec=Clock_Next_Sec    )
   if(masterproc) then
    write(iulog,*) 'CLOCKING: Reading analyses:',trim(Clock_Path)//trim(Clock_File)
   endif

   ! Rotate Clock_ObsInd() indices for new data, then update 
   ! the Nudge observation arrays with analysis data at the 
   ! NEXT==Clock_ObsInd(1) time.
   !----------------------------------------------------------
   if(dycore_is('UNSTRUCTURED')) then
     call clocking_update_analyses_se (trim(Clock_Path)//trim(Clock_File))
   else !if(dycore_is('LR')) then
     call clocking_update_analyses_fv (trim(Clock_Path)//trim(Clock_File))
   endif

   ! Initialize Nudging Coeffcient profiles in local arrays
   ! Load zeros into clocking arrays
   !------------------------------------------------------
   do lchnk=begchunk,endchunk
     ncol=get_ncols_p(lchnk)
     do icol=1,ncol
       rlat=get_rlat_p(lchnk,icol)*180._r8/SHR_CONST_PI
       rlon=get_rlon_p(lchnk,icol)*180._r8/SHR_CONST_PI

       call clocking_set_profile(rlat,rlon,Clock_MUprof,Wprof,pver)
       Clock_MUtau(icol,:,lchnk)=Wprof(:)
       call clocking_set_profile(rlat,rlon,Clock_LAMBDACprof,Wprof,pver)
       Clock_LAMBDACtau(icol,:,lchnk)=Wprof(:)
       call clocking_set_profile(rlat,rlon,Clock_ICSWPprof,Wprof,pver)
       Clock_ICSWPtau(icol,:,lchnk)=Wprof(:)
       call clocking_set_profile(rlat,rlon,Clock_ICLWPprof,Wprof,pver)
       Clock_ICLWPtau(icol,:,lchnk)=Wprof(:)
       call clocking_set_profile(rlat,rlon,Clock_ICIWPprof,Wprof,pver)
       Clock_ICIWPtau(icol,:,lchnk)=Wprof(:)
       call clocking_set_profile(rlat,rlon,Clock_DESprof,Wprof,pver)
       Clock_DEStau(icol,:,lchnk)=Wprof(:)
       call clocking_set_profile(rlat,rlon,Clock_DEIprof,Wprof,pver)
       Clock_DEItau(icol,:,lchnk)=Wprof(:)
       call clocking_set_profile(rlat,rlon,Clock_CLDprof,Wprof,pver)
       Clock_CLDtau(icol,:,lchnk)=Wprof(:)
       call clocking_set_profile(rlat,rlon,Clock_CLDFSNOWprof,Wprof,pver)
       Clock_CLDFSNOWtau(icol,:,lchnk)=Wprof(:)
       call clocking_set_profile(rlat,rlon,Clock_Qprof,Wprof,pver)
       Clock_Qtau(icol,:,lchnk)=Wprof(:)

       Clock_PStau(icol,lchnk)=clocking_set_PSprofile(rlat,rlon,Clock_PSprof)
     end do
     Clock_MUtau      (:ncol,:pver,lchnk) =                                    &
     Clock_MUtau      (:ncol,:pver,lchnk) * Clock_MUcoef      /float(Clock_Step)
     Clock_LAMBDACtau (:ncol,:pver,lchnk) =                                    &
     Clock_LAMBDACtau (:ncol,:pver,lchnk) * Clock_LAMBDACcoef /float(Clock_Step)
     Clock_ICSWPtau   (:ncol,:pver,lchnk) =                                    &
     Clock_ICSWPtau   (:ncol,:pver,lchnk) * Clock_ICSWPcoef   /float(Clock_Step)
     Clock_ICLWPtau   (:ncol,:pver,lchnk) =                                    &
     Clock_ICLWPtau   (:ncol,:pver,lchnk) * Clock_ICLWPcoef   /float(Clock_Step)
     Clock_ICIWPtau   (:ncol,:pver,lchnk) =                                    &
     Clock_ICIWPtau   (:ncol,:pver,lchnk) * Clock_ICIWPcoef   /float(Clock_Step)
     Clock_DEStau     (:ncol,:pver,lchnk) =                                    &
     Clock_DEStau     (:ncol,:pver,lchnk) * Clock_DEScoef     /float(Clock_Step)
     Clock_DEItau     (:ncol,:pver,lchnk) =                                    &
     Clock_DEItau     (:ncol,:pver,lchnk) * Clock_DEIcoef     /float(Clock_Step)
     Clock_CLDtau     (:ncol,:pver,lchnk) =                                    &
     Clock_CLDtau     (:ncol,:pver,lchnk) * Clock_CLDcoef     /float(Clock_Step)
     Clock_CLDFSNOWtau(:ncol,:pver,lchnk) =                                    &
     Clock_CLDFSNOWtau(:ncol,:pver,lchnk) * Clock_CLDFSNOWcoef/float(Clock_Step)
     Clock_Qtau(:ncol,:pver,lchnk) =                             &
     Clock_Qtau(:ncol,:pver,lchnk) * Clock_Qcoef/float(Clock_Step)
     Clock_PStau(:ncol,lchnk)=                             &
     Clock_PStau(:ncol,lchnk)* Clock_PScoef/float(Clock_Step)

     Clock_MUstep      (:pcols,:pver,lchnk)=0._r8
     Clock_LAMBDACstep (:pcols,:pver,lchnk)=0._r8
     Clock_ICSWPstep   (:pcols,:pver,lchnk)=0._r8
     Clock_ICLWPstep   (:pcols,:pver,lchnk)=0._r8
     Clock_ICIWPstep   (:pcols,:pver,lchnk)=0._r8
     Clock_DESstep     (:pcols,:pver,lchnk)=0._r8
     Clock_DEIstep     (:pcols,:pver,lchnk)=0._r8
     Clock_CLDstep     (:pcols,:pver,lchnk)=0._r8
     Clock_CLDFSNOWstep(:pcols,:pver,lchnk)=0._r8
     Clock_Qstep       (:pcols,:pver,lchnk)=0._r8
     Clock_PSstep      (:pcols,lchnk)=0._r8
     Target_MU      (:pcols,:pver,lchnk)=0._r8
     Target_LAMBDAC (:pcols,:pver,lchnk)=0._r8
     Target_ICSWP   (:pcols,:pver,lchnk)=0._r8
     Target_ICLWP   (:pcols,:pver,lchnk)=0._r8
     Target_ICIWP   (:pcols,:pver,lchnk)=0._r8
     Target_DES     (:pcols,:pver,lchnk)=0._r8
     Target_DEI     (:pcols,:pver,lchnk)=0._r8
     Target_CLD     (:pcols,:pver,lchnk)=0._r8
     Target_CLDFSNOW(:pcols,:pver,lchnk)=0._r8
     Target_Q       (:pcols,:pver,lchnk)=0._r8
     Target_PS      (:pcols,lchnk)=0._r8
   end do

   ! End Routine
   !------------
   return
  end subroutine ! clocking_init
  !================================================================


  !================================================================
  subroutine clocking_timestep_init(phys_state)
   ! 
   ! CLOCKING_TIMESTEP_INIT: 
   !                 Check the current time and update Model/Nudging 
   !                 arrays when necessary. Toggle the Nudging flag
   !                 when the time is withing the clocking window.
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
   integer Year,Month,Day,Sec
   integer YMD1,YMD2,YMD
   logical Update_Model,Update_Clock,Sync_Error
   logical After_Beg   ,Before_End
   integer lchnk,ncol,indw

   type(ESMF_Time)         Date1,Date2
   type(ESMF_TimeInterval) DateDiff
   integer                 DeltaT
   real(r8)                Tscale
   real(r8)                Tfrac
   integer                 rc
   integer                 nn
   integer                 kk
   real(r8)                Sbar,Qbar,Wsum
   integer                 dtime

   ! Check if Cloud Locking is initialized
   !---------------------------------
   if(.not.Clock_Initialized) then
     call endrun('clocking_timestep_init:: Cloud Locking NOT Initialized')
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
   YMD1=(Clock_Beg_Year*10000) + (Clock_Beg_Month*100) + Clock_Beg_Day
   call timemgr_time_ge(YMD1,Clock_Beg_Sec,         &
                        YMD ,Sec          ,After_Beg)

   YMD1=(Clock_End_Year*10000) + (Clock_End_Month*100) + Clock_End_Day
   call timemgr_time_ge(YMD ,Sec,                    &
                        YMD1,Clock_End_Sec,Before_End)

   !--------------------------------------------------------------
   ! When past the NEXT time, Update Model Arrays and time indices
   !--------------------------------------------------------------
   YMD1=(Model_Next_Year*10000) + (Model_Next_Month*100) + Model_Next_Day
   call timemgr_time_ge(YMD1,Model_Next_Sec,            &
                        YMD ,Sec           ,Update_Model)

   if((Before_End).and.(Update_Model)) then
     ! Increment the Model times by the current interval
     !---------------------------------------------------
     Model_Curr_Year =Model_Next_Year
     Model_Curr_Month=Model_Next_Month
     Model_Curr_Day  =Model_Next_Day
     Model_Curr_Sec  =Model_Next_Sec
     YMD1=(Model_Curr_Year*10000) + (Model_Curr_Month*100) + Model_Curr_Day
     call timemgr_time_inc(YMD1,Model_Curr_Sec,              &
                           YMD2,Model_Next_Sec,Model_Step,0,0)

     ! Check for Sync Error where NEXT model time after the update
     ! is before the current time. If so, reset the next model 
     ! time to a Model_Step after the current time.
     !--------------------------------------------------------------
     call timemgr_time_ge(YMD2,Model_Next_Sec,            &
                          YMD ,Sec           ,Sync_Error)
     if(Sync_Error) then
       Model_Curr_Year =Year
       Model_Curr_Month=Month
       Model_Curr_Day  =Day
       Model_Curr_Sec  =Sec
       call timemgr_time_inc(YMD ,Model_Curr_Sec,              &
                             YMD2,Model_Next_Sec,Model_Step,0,0)
       write(iulog,*) 'CLOCKING: WARNING - Model_Time Sync ERROR... CORRECTED'
     endif
     Model_Next_Year =(YMD2/10000)
     YMD2            = YMD2-(Model_Next_Year*10000)
     Model_Next_Month=(YMD2/100)
     Model_Next_Day  = YMD2-(Model_Next_Month*100)

     ! Load values at Current into the Model arrays
     !-----------------------------------------------
!IMPLEMENT?     call cnst_get_ind('Q',indw)
!IMPLEMENT?     do lchnk=begchunk,endchunk
!IMPLEMENT?       ncol=phys_state(lchnk)%ncol
!IMPLEMENT?       Model_MU      (:ncol,:pver,lchnk)=phys_state(lchnk)%u(:ncol,:pver)
!IMPLEMENT?       Model_LAMBDAC (:ncol,:pver,lchnk)=phys_state(lchnk)%v(:ncol,:pver)
!IMPLEMENT?       Model_ICSWP   (:ncol,:pver,lchnk)=phys_state(lchnk)%t(:ncol,:pver)
!IMPLEMENT?       Model_ICLWP   (:ncol,:pver,lchnk)=phys_state(lchnk)%t(:ncol,:pver)
!IMPLEMENT?       Model_ICIWP   (:ncol,:pver,lchnk)=phys_state(lchnk)%t(:ncol,:pver)
!IMPLEMENT?       Model_DES     (:ncol,:pver,lchnk)=phys_state(lchnk)%t(:ncol,:pver)
!IMPLEMENT?       Model_DEI     (:ncol,:pver,lchnk)=phys_state(lchnk)%t(:ncol,:pver)
!IMPLEMENT?       Model_CLD     (:ncol,:pver,lchnk)=phys_state(lchnk)%t(:ncol,:pver)
!IMPLEMENT?       Model_CLDFSNOW(:ncol,:pver,lchnk)=phys_state(lchnk)%t(:ncol,:pver)
!IMPLEMENT?       Model_Q       (:ncol,:pver,lchnk)=phys_state(lchnk)%q(:ncol,:pver,indw)
!IMPLEMENT?       Model_PS      (:ncol,lchnk)=phys_state(lchnk)%ps(:ncol)
!IMPLEMENT?     end do

   endif ! ((Before_End).and.(Update_Model)) then

   !----------------------------------------------------------------
   ! When past the NEXT time, Update Clocking Arrays and time indices
   !----------------------------------------------------------------
   YMD1=(Clock_Next_Year*10000) + (Clock_Next_Month*100) + Clock_Next_Day
   call timemgr_time_ge(YMD1,Clock_Next_Sec,            &
                        YMD ,Sec           ,Update_Clock)

   if((Before_End).and.(Update_Clock)) then
     ! Increment the Clock times by the current interval
     !---------------------------------------------------
     Clock_Curr_Year =Clock_Next_Year
     Clock_Curr_Month=Clock_Next_Month
     Clock_Curr_Day  =Clock_Next_Day
     Clock_Curr_Sec  =Clock_Next_Sec
     YMD1=(Clock_Curr_Year*10000) + (Clock_Curr_Month*100) + Clock_Curr_Day
     call timemgr_time_inc(YMD1,Clock_Curr_Sec,              &
                           YMD2,Clock_Next_Sec,Clock_Step,0,0)
     Clock_Next_Year =(YMD2/10000)
     YMD2            = YMD2-(Clock_Next_Year*10000)
     Clock_Next_Month=(YMD2/100)
     Clock_Next_Day  = YMD2-(Clock_Next_Month*100)

     ! Set the analysis filename at the NEXT time.
     !---------------------------------------------------------------
     Clock_File=interpret_filename_spec(Clock_File_Template      , &
                                         yr_spec=Clock_Next_Year , &
                                        mon_spec=Clock_Next_Month, &
                                        day_spec=Clock_Next_Day  , &
                                        sec_spec=Clock_Next_Sec    )
     if(masterproc) then
      write(iulog,*) 'CLOCKING: Reading analyses:',trim(Clock_Path)//trim(Clock_File)
     endif

     ! Rotate Clock_ObsInd() indices for new data, then update 
     ! the Clock observation arrays with analysis data at the 
     ! NEXT==Clock_ObsInd(1) time.
     !----------------------------------------------------------
     if(dycore_is('UNSTRUCTURED')) then
       call clocking_update_analyses_se (trim(Clock_Path)//trim(Clock_File))
     else !if(dycore_is('LR')) then
       call clocking_update_analyses_fv (trim(Clock_Path)//trim(Clock_File))
     endif
   endif ! ((Before_End).and.(Update_Clock)) then

   !----------------------------------------------------------------
   ! Toggle Clock flag when the time interval is between 
   ! beginning and ending times, and all of the analyses files exist.
   !----------------------------------------------------------------
   if((After_Beg).and.(Before_End)) then
     if(Clock_Force_Opt.eq.0) then
       ! Verify that the NEXT analyses are available
       !---------------------------------------------
       Clock_ON=Clock_File_Present(Clock_ObsInd(1))
     elseif(Clock_Force_Opt.eq.1) then
       ! Verify that the CURR and NEXT analyses are available
       !-----------------------------------------------------
       Clock_ON=(Clock_File_Present(Clock_ObsInd(1)).and. &
                 Clock_File_Present(Clock_ObsInd(2))      )
     else
       ! Verify that the ALL analyses are available
       !---------------------------------------------
       Clock_ON=.true.
       do nn=1,Clock_NumObs
         if(.not.Clock_File_Present(nn)) Clock_ON=.false.
       end do
     endif
     if(.not.Clock_ON) then
       if(masterproc) then
         write(iulog,*) 'CLOCKING: WARNING - analyses file NOT FOUND. Switching '
         write(iulog,*) 'CLOCKING:           clocking OFF to coast thru the gap. '
       endif
     endif
   else
     Clock_ON=.false.
   endif

   !-------------------------------------------------------
   ! HERE Implement time dependence of Clock Coefs HERE
   !-------------------------------------------------------


   !---------------------------------------------------
   ! If Data arrays have changed update stepping arrays
   !---------------------------------------------------
   if((Before_End).and.((Update_Clock).or.(Update_Model))) then

     ! Now Load the Target values for clocking tendencies
     !---------------------------------------------------
     if(Clock_Force_Opt.eq.0) then
       ! Target is OBS data at NEXT time
       !----------------------------------
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         Target_MU      (:ncol,:pver,lchnk)=Nobs_MU      (:ncol,:pver,lchnk,Clock_ObsInd(1))
         Target_LAMBDAC (:ncol,:pver,lchnk)=Nobs_LAMBDAC (:ncol,:pver,lchnk,Clock_ObsInd(1))
         Target_ICSWP   (:ncol,:pver,lchnk)=Nobs_ICSWP   (:ncol,:pver,lchnk,Clock_ObsInd(1))
         Target_ICLWP   (:ncol,:pver,lchnk)=Nobs_ICLWP   (:ncol,:pver,lchnk,Clock_ObsInd(1))
         Target_ICIWP   (:ncol,:pver,lchnk)=Nobs_ICIWP   (:ncol,:pver,lchnk,Clock_ObsInd(1))
         Target_DES     (:ncol,:pver,lchnk)=Nobs_DES     (:ncol,:pver,lchnk,Clock_ObsInd(1))
         Target_DEI     (:ncol,:pver,lchnk)=Nobs_DEI     (:ncol,:pver,lchnk,Clock_ObsInd(1))
         Target_CLD     (:ncol,:pver,lchnk)=Nobs_CLD     (:ncol,:pver,lchnk,Clock_ObsInd(1))
         Target_CLDFSNOW(:ncol,:pver,lchnk)=Nobs_CLDFSNOW(:ncol,:pver,lchnk,Clock_ObsInd(1))
         Target_Q       (:ncol,:pver,lchnk)=Nobs_Q       (:ncol,:pver,lchnk,Clock_ObsInd(1))
         Target_PS      (:ncol      ,lchnk)=Nobs_PS      (:ncol      ,lchnk,Clock_ObsInd(1))
       end do
     elseif(Clock_Force_Opt.eq.1) then
       ! Target is linear interpolation of OBS data CURR<-->NEXT time    
       !---------------------------------------------------------------
       call ESMF_TimeSet(Date1,YY=Year,MM=Month,DD=Day,S=Sec)
       call ESMF_TimeSet(Date2,YY=Clock_Next_Year,MM=Clock_Next_Month, &
                               DD=Clock_Next_Day , S=Clock_Next_Sec    )
       DateDiff =Date2-Date1
       call ESMF_TimeIntervalGet(DateDiff,S=DeltaT,rc=rc)
       Tfrac= float(DeltaT)/float(Clock_Step)
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         Target_MU      (:ncol,:pver,lchnk)=(1._r8-Tfrac)*Nobs_MU      (:ncol,:pver,lchnk,Clock_ObsInd(1)) &
                                                  +Tfrac *Nobs_MU      (:ncol,:pver,lchnk,Clock_ObsInd(2))
         Target_LAMBDAC (:ncol,:pver,lchnk)=(1._r8-Tfrac)*Nobs_LAMBDAC (:ncol,:pver,lchnk,Clock_ObsInd(1)) &
                                                  +Tfrac *Nobs_LAMBDAC (:ncol,:pver,lchnk,Clock_ObsInd(2))
         Target_ICSWP   (:ncol,:pver,lchnk)=(1._r8-Tfrac)*Nobs_ICSWP   (:ncol,:pver,lchnk,Clock_ObsInd(1)) &
                                                  +Tfrac *Nobs_ICSWP   (:ncol,:pver,lchnk,Clock_ObsInd(2))
         Target_ICLWP   (:ncol,:pver,lchnk)=(1._r8-Tfrac)*Nobs_ICLWP   (:ncol,:pver,lchnk,Clock_ObsInd(1)) &
                                                  +Tfrac *Nobs_ICLWP   (:ncol,:pver,lchnk,Clock_ObsInd(2))
         Target_ICIWP   (:ncol,:pver,lchnk)=(1._r8-Tfrac)*Nobs_ICIWP   (:ncol,:pver,lchnk,Clock_ObsInd(1)) &
                                                  +Tfrac *Nobs_ICIWP   (:ncol,:pver,lchnk,Clock_ObsInd(2))
         Target_DES     (:ncol,:pver,lchnk)=(1._r8-Tfrac)*Nobs_DES     (:ncol,:pver,lchnk,Clock_ObsInd(1)) &
                                                  +Tfrac *Nobs_DES     (:ncol,:pver,lchnk,Clock_ObsInd(2))
         Target_DEI     (:ncol,:pver,lchnk)=(1._r8-Tfrac)*Nobs_DEI     (:ncol,:pver,lchnk,Clock_ObsInd(1)) &
                                                  +Tfrac *Nobs_DEI     (:ncol,:pver,lchnk,Clock_ObsInd(2))
         Target_CLD     (:ncol,:pver,lchnk)=(1._r8-Tfrac)*Nobs_CLD     (:ncol,:pver,lchnk,Clock_ObsInd(1)) &
                                                  +Tfrac *Nobs_CLD     (:ncol,:pver,lchnk,Clock_ObsInd(2))
         Target_CLDFSNOW(:ncol,:pver,lchnk)=(1._r8-Tfrac)*Nobs_CLDFSNOW(:ncol,:pver,lchnk,Clock_ObsInd(1)) &
                                                  +Tfrac *Nobs_CLDFSNOW(:ncol,:pver,lchnk,Clock_ObsInd(2))
         Target_Q       (:ncol,:pver,lchnk)=(1._r8-Tfrac)*Nobs_Q       (:ncol,:pver,lchnk,Clock_ObsInd(1)) &
                                                  +Tfrac *Nobs_Q       (:ncol,:pver,lchnk,Clock_ObsInd(2))
         Target_PS      (:ncol      ,lchnk)=(1._r8-Tfrac)*Nobs_PS      (:ncol      ,lchnk,Clock_ObsInd(1)) &
                                                  +Tfrac *Nobs_PS      (:ncol      ,lchnk,Clock_ObsInd(2))
       end do
     else
       write(iulog,*) 'CLOCKING: Unknown Clock_Force_Opt=',Clock_Force_Opt
       call endrun('clocking_timestep_init:: ERROR unknown Clocking_Force_Opt')
     endif

     ! Set Tscale for the specified Forcing Option 
     !-----------------------------------------------
     if(Clock_TimeScale_Opt.eq.0) then
       Tscale=1._r8
     elseif(Clock_TimeScale_Opt.eq.1) then
       call ESMF_TimeSet(Date1,YY=Year,MM=Month,DD=Day,S=Sec)
       call ESMF_TimeSet(Date2,YY=Clock_Next_Year,MM=Clock_Next_Month, &
                               DD=Clock_Next_Day , S=Clock_Next_Sec    )
       DateDiff =Date2-Date1
       call ESMF_TimeIntervalGet(DateDiff,S=DeltaT,rc=rc)
       Tscale=float(Clock_Step)/float(DeltaT)
     else
       write(iulog,*) 'CLOCKING: Unknown Clock_TimeScale_Opt=',Clock_TimeScale_Opt
       call endrun('clocking_timestep_init:: ERROR unknown Clocking_TimeScale_Opt')
     endif

     ! Update the clocking tendencies
     !--------------------------------
     do lchnk=begchunk,endchunk
       ncol=phys_state(lchnk)%ncol
       Clock_MUstep      (:ncol,:pver,lchnk)=(Target_MU         (:ncol,:pver,lchnk)  &
                                              -Model_MU         (:ncol,:pver,lchnk)) &
                                       *Tscale*Clock_MUtau      (:ncol,:pver,lchnk)
       Clock_LAMBDACstep (:ncol,:pver,lchnk)=(Target_LAMBDAC    (:ncol,:pver,lchnk)  &
                                              -Model_LAMBDAC    (:ncol,:pver,lchnk)) &
                                       *Tscale*Clock_LAMBDACtau (:ncol,:pver,lchnk)
       Clock_ICSWPstep   (:ncol,:pver,lchnk)=(Target_ICSWP      (:ncol,:pver,lchnk)  &
                                              -Model_ICSWP      (:ncol,:pver,lchnk)) &
                                       *Tscale*Clock_ICSWPtau   (:ncol,:pver,lchnk)
       Clock_ICLWPstep   (:ncol,:pver,lchnk)=(Target_ICLWP      (:ncol,:pver,lchnk)  &
                                              -Model_ICLWP      (:ncol,:pver,lchnk)) &
                                       *Tscale*Clock_ICLWPtau   (:ncol,:pver,lchnk)
       Clock_ICIWPstep   (:ncol,:pver,lchnk)=(Target_ICIWP      (:ncol,:pver,lchnk)  &
                                              -Model_ICIWP      (:ncol,:pver,lchnk)) &
                                       *Tscale*Clock_ICIWPtau   (:ncol,:pver,lchnk)
       Clock_DESstep     (:ncol,:pver,lchnk)=(Target_DES        (:ncol,:pver,lchnk)  &
                                              -Model_DES        (:ncol,:pver,lchnk)) &
                                       *Tscale*Clock_DEStau     (:ncol,:pver,lchnk)
       Clock_DEIstep     (:ncol,:pver,lchnk)=(Target_DEI        (:ncol,:pver,lchnk)  &
                                              -Model_DEI        (:ncol,:pver,lchnk)) &
                                       *Tscale*Clock_DEItau     (:ncol,:pver,lchnk)
       Clock_CLDstep     (:ncol,:pver,lchnk)=(Target_CLD        (:ncol,:pver,lchnk)  &
                                              -Model_CLD        (:ncol,:pver,lchnk)) &
                                       *Tscale*Clock_CLDtau     (:ncol,:pver,lchnk)
       Clock_CLDFSNOWstep(:ncol,:pver,lchnk)=(Target_CLDFSNOW   (:ncol,:pver,lchnk)  &
                                              -Model_CLDFSNOW   (:ncol,:pver,lchnk)) &
                                       *Tscale*Clock_CLDFSNOWtau(:ncol,:pver,lchnk)
       Clock_Qstep       (:ncol,:pver,lchnk)=(Target_Q          (:ncol,:pver,lchnk)      &
                                              -Model_Q          (:ncol,:pver,lchnk))     &
                                       *Tscale*Clock_Qtau       (:ncol,:pver,lchnk)
       Clock_PSstep      (:ncol,      lchnk)=(Target_PS         (:ncol,lchnk)      &
                                              -Model_PS         (:ncol,lchnk))     &
                                       *Tscale*Clock_PStau      (:ncol,lchnk)
     end do

   endif ! ((Before_End).and.((Update_Clock).or.(Update_Model))) then

   ! End Routine
   !------------
   return
  end subroutine ! clocking_timestep_init
  !================================================================


  !================================================================
  subroutine clocking_timestep_tend(phys_state,phys_tend)
   ! 
   ! CLOCKING_TIMESTEP_TEND: 
   !                If Cloud Locking is ON, return the Clocking contributions 
   !                to forcing using the current contents of the Clocking 
   !                arrays. Send output to the cam history module as well.
   !===============================================================
   use physconst    ,only: cpair
   use physics_types,only: physics_state,physics_ptend,physics_ptend_init
   use constituents ,only: cnst_get_ind,pcnst
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use cam_history  ,only: outfld

   ! Arguments
   !-------------
   type(physics_state), intent(in) :: phys_state
   type(physics_ptend), intent(out):: phys_tend

   ! Local values
   !--------------------
   integer indw,ncol,lchnk
   logical lq(pcnst)

   call cnst_get_ind('Q',indw)
   lq(:)   =.false.
   lq(indw)=.true.
   call physics_ptend_init(phys_tend,phys_state%psetcols,'clocking',lu=.true.,lv=.true.,ls=.true.,lq=lq)

   if(Clock_ON) then
     lchnk=phys_state%lchnk
     ncol =phys_state%ncol
!IMPLEMET?     phys_tend%u(:ncol,:pver)     =Clock_MUstep      (:ncol,:pver,lchnk)
!IMPLEMET?     phys_tend%v(:ncol,:pver)     =Clock_LAMBDACstep (:ncol,:pver,lchnk)
!IMPLEMET?     phys_tend%s(:ncol,:pver)     =Clock_ICSWPstep   (:ncol,:pver,lchnk)
!IMPLEMET?     phys_tend%s(:ncol,:pver)     =Clock_ICLWPstep   (:ncol,:pver,lchnk)
!IMPLEMET?     phys_tend%s(:ncol,:pver)     =Clock_ICIWPstep   (:ncol,:pver,lchnk)
!IMPLEMET?     phys_tend%s(:ncol,:pver)     =Clock_DESstep     (:ncol,:pver,lchnk)
!IMPLEMET?     phys_tend%s(:ncol,:pver)     =Clock_DEIstep     (:ncol,:pver,lchnk)
!IMPLEMET?     phys_tend%s(:ncol,:pver)     =Clock_CLDstep     (:ncol,:pver,lchnk)
!IMPLEMET?     phys_tend%s(:ncol,:pver)     =Clock_CLDFSNOWstep(:ncol,:pver,lchnk)
!IMPLEMET?     phys_tend%q(:ncol,:pver,indw)=Clock_Qstep       (:ncol,:pver,lchnk)

!IMPLEMENT?     call outfld( 'Nudge_U',phys_tend%u                ,pcols,lchnk)
!IMPLEMENT?     call outfld( 'Nudge_V',phys_tend%v                ,pcols,lchnk)
!IMPLEMENT?     call outfld( 'Nudge_T',phys_tend%s/cpair          ,pcols,lchnk)
!IMPLEMENT?     call outfld( 'Nudge_Q',phys_tend%q(1,1,indw)      ,pcols,lchnk)
     call outfld('Target_MU'      ,Target_MU      (:,:,lchnk),pcols,lchnk)
     call outfld('Target_LAMBDAC' ,Target_LAMBDAC (:,:,lchnk),pcols,lchnk)
     call outfld('Target_ICSWP'   ,Target_ICSWP   (:,:,lchnk),pcols,lchnk)
     call outfld('Target_ICLWP'   ,Target_ICLWP   (:,:,lchnk),pcols,lchnk)
     call outfld('Target_ICIWP'   ,Target_ICIWP   (:,:,lchnk),pcols,lchnk)
     call outfld('Target_DES'     ,Target_DES     (:,:,lchnk),pcols,lchnk)
     call outfld('Target_DEI'     ,Target_DEI     (:,:,lchnk),pcols,lchnk)
     call outfld('Target_CLD'     ,Target_CLD     (:,:,lchnk),pcols,lchnk)
     call outfld('Target_CLDFSNOW',Target_CLDFSNOW(:,:,lchnk),pcols,lchnk)
     call outfld('Target_Q'       ,Target_Q       (:,:,lchnk),pcols,lchnk)
   endif

   ! End Routine
   !------------
   return
  end subroutine ! clocking_timestep_tend
  !================================================================


  !================================================================
  subroutine clocking_update_analyses_se(anal_file)
   ! 
   ! CLOCKING_UPDATE_ANALYSES_SE: 
   !                 Open the given analyses data file, read in 
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid ,only: pver,begchunk
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   integer lev
   integer ncol,plev,istat
   integer ncid,varid
   real(r8) Xanal(Clock_ncol,Clock_nlev)
   real(r8) PSanal(Clock_ncol)
   real(r8) Lat_anal(Clock_ncol)
   real(r8) Lon_anal(Clock_ncol)
   integer  nn,Nindex

   ! Rotate Clock_ObsInd() indices, then check the existence of the analyses 
   ! file; broadcast the updated indices and file status to all the other MPI nodes. 
   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     Nindex=Clock_ObsInd(Clock_NumObs)
     do nn=Clock_NumObs,2,-1
       Clock_ObsInd(nn)=Clock_ObsInd(nn-1)
     end do
     Clock_ObsInd(1)=Nindex
     inquire(FILE=trim(anal_file),EXIST=Clock_File_Present(Clock_ObsInd(1)))
     write(iulog,*)'CLOCKING: Clock_ObsInd=',Clock_ObsInd
     write(iulog,*)'CLOCKING: Clock_File_Present=',Clock_File_Present
   endif
#ifdef SPMD
   call mpibcast(Clock_File_Present, Clock_NumObs, mpilog, 0, mpicom)
   call mpibcast(Clock_ObsInd      , Clock_NumObs, mpiint, 0, mpicom)
#endif
   if(.not.Clock_File_Present(Clock_ObsInd(1))) return

   ! masterporc does all of the work here
   !-----------------------------------------
   if(masterproc) then
   
     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(anal_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(anal_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

     ! Read in Dimensions
     !--------------------
     istat=nf90_inq_dimid(ncid,'ncol',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=ncol)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

     istat=nf90_inq_dimid(ncid,'lev',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=plev)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

     istat=nf90_inq_varid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Lon_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

     istat=nf90_inq_varid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Lat_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

     if((Clock_ncol.ne.ncol).or.(plev.ne.pver)) then
      write(iulog,*) 'ERROR: clocking_update_analyses_se: ncol=',ncol,' Clock_ncol=',Clock_ncol
      write(iulog,*) 'ERROR: clocking_update_analyses_se: plev=',plev,' pver=',pver
      call endrun('clocking_update_analyses_se: analyses dimension mismatch')
     endif

     ! Read in and scatter data arrays
     !----------------------------------
     istat=nf90_inq_varid(ncid,'MU_rad',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_ncol,Xanal,    &
                               Nobs_MU(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'LAMBDAC_rad',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_ncol,Xanal,    &
                               Nobs_LAMBDAC(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'ICSWP_rad',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_ncol,Xanal,    &
                               Nobs_ICSWP(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'ICLWP_rad',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_ncol,Xanal,    &
                               Nobs_ICLWP(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'ICIWP_rad',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_ncol,Xanal,    &
                               Nobs_ICIWP(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'DES_rad',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_ncol,Xanal,    &
                               Nobs_DES(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'DEI_rad',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_ncol,Xanal,    &
                               Nobs_DEI(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'CLD_rad',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_ncol,Xanal,    &
                               Nobs_CLD(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'CLDFSNOW_rad',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_ncol,Xanal,    &
                               Nobs_CLDFSNOW(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'Q',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_ncol,Xanal,    &
                               Nobs_Q(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
    istat=nf90_inq_varid(ncid,'PS',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_SE')
    endif
    istat=nf90_get_var(ncid,varid,PSanal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_SE')
    endif

     ! Close the analyses file
     !-----------------------
     istat=nf90_close(ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,1,1,Clock_ncol,PSanal,           &
                               Nobs_PS(1,begchunk,Clock_ObsInd(1)))

   ! End Routine
   !------------
   return
  end subroutine ! clocking_update_analyses_se
  !================================================================


  !================================================================
  subroutine clocking_update_analyses_fv(anal_file)
   ! 
   ! CLOCKING_UPDATE_ANALYSES_FV: 
   !                 Open the given analyses data file, read in 
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid ,only: pver,begchunk
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   integer lev
   integer nlon,nlat,plev,istat
   integer ncid,varid
   integer ilat,ilon,ilev
   real(r8) Xanal(Clock_nlon,Clock_nlat,Clock_nlev)
   real(r8) PSanal(Clock_nlon,Clock_nlat)
   real(r8) Lat_anal(Clock_nlat)
   real(r8) Lon_anal(Clock_nlon)
   real(r8) Xtrans(Clock_nlon,Clock_nlev,Clock_nlat)
   integer  nn,Nindex

   ! Rotate Clock_ObsInd() indices, then check the existence of the analyses 
   ! file; broadcast the updated indices and file status to all the other MPI nodes. 
   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     Nindex=Clock_ObsInd(Clock_NumObs)
     do nn=Clock_NumObs,2,-1
       Clock_ObsInd(nn)=Clock_ObsInd(nn-1)
     end do
     Clock_ObsInd(1)=Nindex
     inquire(FILE=trim(anal_file),EXIST=Clock_File_Present(Clock_ObsInd(1)))
     write(iulog,*)'CLOCKING: Clock_ObsInd=',Clock_ObsInd
     write(iulog,*)'CLOCKING: Clock_File_Present=',Clock_File_Present
   endif
#ifdef SPMD
   call mpibcast(Clock_File_Present, Clock_NumObs, mpilog, 0, mpicom)
   call mpibcast(Clock_ObsInd      , Clock_NumObs, mpiint, 0, mpicom)
#endif
   if(.not.Clock_File_Present(Clock_ObsInd(1))) return

   ! masterporc does all of the work here
   !-----------------------------------------
   if(masterproc) then
   
     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(anal_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(anal_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     ! Read in Dimensions
     !--------------------
     istat=nf90_inq_dimid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlon)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_dimid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlat)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_dimid(ncid,'lev',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=plev)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_varid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Lon_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_varid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Lat_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     if((Clock_nlon.ne.nlon).or.(Clock_nlat.ne.nlat).or.(plev.ne.pver)) then
      write(iulog,*) 'ERROR: clocking_update_analyses_fv: nlon=',nlon,' Clock_nlon=',Clock_nlon
      write(iulog,*) 'ERROR: clocking_update_analyses_fv: nlat=',nlat,' Clock_nlat=',Clock_nlat
      write(iulog,*) 'ERROR: clocking_update_analyses_fv: plev=',plev,' pver=',pver
      call endrun('clocking_update_analyses_fv: analyses dimension mismatch')
     endif

     ! Read in, transpose lat/lev indices, 
     ! and scatter data arrays
     !----------------------------------
     istat=nf90_inq_varid(ncid,'MU_rad',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_nlon,Xtrans,   &
                               Nobs_MU(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'LAMBDAC_rad',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_nlon,Xtrans,   &
                               Nobs_LAMBDAC(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'ICSWP_rad',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_nlon,Xtrans,   &
                               Nobs_ICSWP(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'ICLWP_rad',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_nlon,Xtrans,   &
                               Nobs_ICLWP(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'ICIWP_rad',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_nlon,Xtrans,   &
                               Nobs_ICIWP(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'DES_rad',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_nlon,Xtrans,   &
                               Nobs_DES(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'DEI_rad',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_nlon,Xtrans,   &
                               Nobs_DEI(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'CLD_rad',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_nlon,Xtrans,   &
                               Nobs_CLD(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'CLDFSNOW_rad',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_nlon,Xtrans,   &
                               Nobs_CLDFSNOW(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'Q',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Clock_nlev,1,Clock_nlon,Xtrans,   &
                               Nobs_Q(1,1,begchunk,Clock_ObsInd(1)))

   if(masterproc) then
    istat=nf90_inq_varid(ncid,'PS',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_SE')
    endif
    istat=nf90_get_var(ncid,varid,PSanal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_SE')
    endif

     ! Close the analyses file
     !-----------------------
     istat=nf90_close(ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,1,1,Clock_nlon,PSanal,           &
                               Nobs_PS(1,begchunk,Clock_ObsInd(1)))

   ! End Routine
   !------------
   return
  end subroutine ! clocking_update_analyses_fv
  !================================================================


  !================================================================
  subroutine clocking_set_profile(rlat,rlon,Clock_prof,Wprof,nlev)
   ! 
   ! CLOCKING_SET_PROFILE: for the given lat,lon, and Clock_prof, set
   !                      the verical profile of window coeffcients.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on clocking strength.
   !===============================================================

   ! Arguments
   !--------------
   integer  nlev,Clock_prof
   real(r8) rlat,rlon
   real(r8) Wprof(nlev)

   ! Local values
   !----------------
   integer  ilev
   real(r8) Hcoef,latx,lonx,Vmax,Vmin
   real(r8) lon_lo,lon_hi,lat_lo,lat_hi,lev_lo,lev_hi

   !---------------
   ! set coeffcient
   !---------------
   if(Clock_prof.eq.0) then
     ! No Clocking
     !-------------
     Wprof(:)=0.0_r8
   elseif(Clock_prof.eq.1) then
     ! Uniform Clocking
     !-----------------
     Wprof(:)=1.0_r8
   elseif(Clock_prof.eq.2) then
     ! Localized Clocking with specified Heaviside window function
     !------------------------------------------------------------
     if(Clock_Hwin_max.le.Clock_Hwin_min) then
       ! For a constant Horizontal window function, 
       ! just set Hcoef to the maximum of Hlo/Hhi.
       !--------------------------------------------
       Hcoef=max(Clock_Hwin_lo,Clock_Hwin_hi)
     else
       ! get lat/lon relative to window center
       !------------------------------------------
       latx=rlat-Clock_Hwin_lat0
       lonx=rlon-Clock_Hwin_lon0
       if(lonx.gt. 180._r8) lonx=lonx-360._r8
       if(lonx.le.-180._r8) lonx=lonx+360._r8

       ! Calcualte RAW window value
       !-------------------------------
       lon_lo=(Clock_Hwin_lonWidthH+lonx)/Clock_Hwin_lonDelta
       lon_hi=(Clock_Hwin_lonWidthH-lonx)/Clock_Hwin_lonDelta
       lat_lo=(Clock_Hwin_latWidthH+latx)/Clock_Hwin_latDelta
       lat_hi=(Clock_Hwin_latWidthH-latx)/Clock_Hwin_latDelta
       Hcoef=((1._r8+tanh(lon_lo))/2._r8)*((1._r8+tanh(lon_hi))/2._r8) &
            *((1._r8+tanh(lat_lo))/2._r8)*((1._r8+tanh(lat_hi))/2._r8)

       ! Scale the horizontal window coef for specfied range of values.
       !--------------------------------------------------------
       Hcoef=(Hcoef-Clock_Hwin_min)/(Clock_Hwin_max-Clock_Hwin_min)
       Hcoef=(1._r8-Hcoef)*Clock_Hwin_lo + Hcoef*Clock_Hwin_hi
     endif

     ! Load the RAW vertical window
     !------------------------------
     do ilev=1,nlev
       lev_lo=(float(ilev)-Clock_Vwin_Lindex)/Clock_Vwin_Ldelta
       lev_hi=(Clock_Vwin_Hindex-float(ilev))/Clock_Vwin_Hdelta
       Wprof(ilev)=((1._r8+tanh(lev_lo))/2._r8)*((1._r8+tanh(lev_hi))/2._r8)
     end do 

     ! Scale the Window function to span the values between Vlo and Vhi:
     !-----------------------------------------------------------------
     Vmax=maxval(Wprof)
     Vmin=minval(Wprof)
     if((Vmax.le.Vmin).or.((Clock_Vwin_Hindex.ge.(nlev+1)).and. &
                           (Clock_Vwin_Lindex.le. 0      )     )) then
       ! For a constant Vertical window function, 
       ! load maximum of Vlo/Vhi into Wprof()
       !--------------------------------------------
       Vmax=max(Clock_Vwin_lo,Clock_Vwin_hi)
       Wprof(:)=Vmax
     else
       ! Scale the RAW vertical window for specfied range of values.
       !--------------------------------------------------------
       Wprof(:)=(Wprof(:)-Vmin)/(Vmax-Vmin)
       Wprof(:)=Clock_Vwin_lo + Wprof(:)*(Clock_Vwin_hi-Clock_Vwin_lo)
     endif

     ! The desired result is the product of the vertical profile 
     ! and the horizontal window coeffcient.
     !----------------------------------------------------
     Wprof(:)=Hcoef*Wprof(:)
   else
     call endrun('clocking_set_profile:: Unknown Clock_prof value')
   endif

   ! End Routine
   !------------
   return
  end subroutine ! clocking_set_profile
  !================================================================


  !================================================================
  real(r8) function clocking_set_PSprofile(rlat,rlon,Clock_PSprof)
   ! 
   ! CLOCKING_SET_PSPROFILE: for the given lat and lon set the surface
   !                      pressure profile value for the specified index.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on clocking strength.
   !===============================================================

   ! Arguments
   !--------------
   real(r8) rlat,rlon
   integer  Clock_PSprof

   ! Local values
   !----------------

   !---------------
   ! set coeffcient
   !---------------
   if(Clock_PSprof.eq.0) then
     ! No Clocking
     !-------------
     clocking_set_PSprofile=0.0_r8
   elseif(Clock_PSprof.eq.1) then
     ! Uniform Clocking
     !-----------------
     clocking_set_PSprofile=1.0_r8
   else
     call endrun('clocking_set_PSprofile:: Unknown Clock_prof value')
   endif

   ! End Routine
   !------------
   return
  end function ! clocking_set_PSprofile
  !================================================================


end module cloud_locking
