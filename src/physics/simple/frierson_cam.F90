module frierson_cam
!-----------------------------------------------------------------------
!
! Purpose: Implement idealized forcings described in 
!          Frierson, et al. (2006), " A Gray-Radiation Aquaplanet
!          Moist GCM, Part I. Static Stability and Eddy Scale"
!          J. Atmos. Sci, Vol 63, 2548-2566.
!
!      NOTE: This is a DIAGNOSTIC/DEVELOPMENT VERSION of this module
!            with options to toggle the use of previously existing 
!            simple formulations for PRCIP, PBL, RATIATION, and SURFACE
!            temperatures. 
!            This module serves as a template for a more general 
!            Simple Model fromulation in which that user will be able 
!            to select processes a-la-carte to establish an atmosphereic 
!            environment suited to thier needs.
!
!============================================================================
  ! Useful modules
  !-------------------
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use shr_const_mod,  only: pi => shr_const_pi
  use physconst,      only: gravit, cappa, rair, cpair, latvap, rh2o, epsilo, rhoh2o, zvir
  use ppgrid,         only: pcols, pver, begchunk, endchunk
  use constituents,   only: pcnst
  use physics_buffer, only: dtype_r8, pbuf_add_field, physics_buffer_desc, &
                            pbuf_set_field, pbuf_get_field
  use camsrfexch,     only: cam_in_t,cam_out_t
  use cam_history,    only: outfld
  use time_manager,   only: is_first_step
  use cam_abortutils, only: endrun
  use spmd_utils,     only: masterproc
  use cam_logfile,    only: iulog
  use hycoef,         only: ps0, etamid
#ifdef SPMD
  use mpishorthand
#endif

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure. 
  !---------------------------------------------------------
  implicit none
  private
  save

  public :: frierson_register
  public :: frierson_readnl
  public :: frierson_init
!  public :: frierson_timestep_init
  public :: frierson_convection_tend
  public :: frierson_condensate_tend
  public :: frierson_pbl_init
  public :: frierson_pbl_tend
  public :: frierson_radiative_tend
  private:: frierson_surface_init

  ! PBL Configuatons
  !------------------
  integer,parameter:: FRIERSON_IMPL_DIAG          = 1
  integer,parameter:: FRIERSON_IMPL_DIAGC         = 2
  integer,parameter:: FRIERSON_IMPL_CONFIG1       = 3
  integer,parameter:: FRIERSON_EXPSFC_DIAG        = 4
  integer,parameter:: FRIERSON_EXPSFC_DIAGC       = 5
  integer,parameter:: FRIERSON_EXPSFC_CONFIG1     = 6
  integer,parameter:: FRIERSON_EXPSFC_CONFIG2     = 7
  integer,parameter:: FRIERSON_EXPSFC             = 8

  ! Tags to identify optional model formulations
  !------------------------------------------------
  integer,parameter:: CONVECTION_NONE          = 0
  integer,parameter:: CONVECTION_FRIERSON      = 1
  integer,parameter:: CONDENSATE_NONE          = 0
  integer,parameter:: CONDENSATE_FRIERSON      = 1
  integer,parameter:: CONDENSATE_FRIERSON_TJ16 = 2
  integer,parameter:: RADIATION_FRIERSON       = 1
  integer,parameter:: SURFACE_INIT_FRIERSON    = 1
  integer,parameter:: SURFACE_UPDATE_NONE      = 0
  integer,parameter:: SURFACE_UPDATE_FRIERSON  = 1

  ! Options selecting which PRECIP, PBL, RADIATION, and SURFACE formulations to use.
  !---------------------------------------------------------------------------------
!  integer,parameter:: PBL_OPT            = FRIERSON_EXPSFC_CONFIG2 
  integer,parameter:: PBL_OPT            = FRIERSON_EXPSFC
  integer,parameter:: CONVECTION_OPT     = CONVECTION_NONE
  integer,parameter:: CONDENSATE_OPT     = CONDENSATE_FRIERSON
  integer,parameter:: RADIATION_OPT      = RADIATION_FRIERSON
  integer,parameter:: SURFACE_INIT_OPT   = SURFACE_INIT_FRIERSON
  integer,parameter:: SURFACE_UPDATE_OPT = SURFACE_UPDATE_FRIERSON

  ! Some Physics buffer indicies
  !-------------------------------
  integer:: prec_pcw_idx = 0
  integer:: prec_dp_idx  = 0
  integer:: relhum_idx   = 0

  ! Global values for Surface Temp, surface fluxes, and radiative heating
  !
  !  NOTE: These values are typically passed around with the cam_in/cam_out
  !        containers. However, the cam_in%sst() values are being munged 
  !        by the coupler between tphyscbc() and tphysac(). Scattered values
  !        are set to 0, I suspect by the application of a land mask. There was also 
  !        a problem with initiailizing the surface temperatures so that they
  !        were not 0 on the first step. Rather than fighting 'the system' to 
  !        communicate values, it is easier and more efficient to just 
  !        maintain the values here to keep Simple Physics simple.
  !----------------------------------------------------------------------
  integer ,allocatable:: LAST(:)
  integer ,allocatable:: CURR(:)        ! Time level indices.
  integer ,allocatable:: NEXT(:)
  real(r8),allocatable:: Tsurf (:,:,:)   ! Surface Temp at 3 time levels
  real(r8),allocatable:: Qsurf (:,:)     ! Surface Q
  real(r8),allocatable:: Fsolar(:,:)     ! Net Solar Heating
  real(r8),allocatable:: Fup   (:,:)     ! Upward Longwave heating
  real(r8),allocatable:: Fdown (:,:)     ! Downward Longwave heating
  real(r8),allocatable:: SHflux(:,:)     ! Sensible Heat flux
  real(r8),allocatable:: LHflux(:,:)     ! Latent Heat Flux
  real(r8),allocatable:: SWflux(:,:)     ! Surface Water flux
  real(r8),allocatable:: TUflux(:,:)     ! U momentum flux
  real(r8),allocatable:: TVflux(:,:)     ! V momentum flux
  real(r8),allocatable:: Evap  (:,:)     ! U momentum flux
  real(r8),allocatable:: Cd    (:,:)     ! V momentum flux
  real(r8),allocatable:: clat  (:,:)     ! latitudes(radians) for columns
  real(r8),allocatable:: Fnet  (:,:)     ! Net Radiative Surface Heating

  ! Global Tuning values
  !------------------------
  real(r8):: frierson_T0         = 273.16_r8
  real(r8):: frierson_E0         = 610.78_r8
  real(r8):: frierson_Erad       = 6.376d6
  real(r8):: frierson_Wind_min   = 1.0d-5
  real(r8):: frierson_Z0         = 3.21d-5
  real(r8):: frierson_Ri_c       = 1.0_r8
  real(r8):: frierson_Karman     = 0.4_r8
  real(r8):: frierson_Fb         = 0.1_r8
  real(r8):: frierson_Rs0        = 938.4_r8
  real(r8):: frierson_DeltaS     = 1.4_r8
  real(r8):: frierson_Tau_eqtr   = 6.0_r8
  real(r8):: frierson_Tau_pole   = 1.5_r8
  real(r8):: frierson_LinFrac    = 0.1_r8
  real(r8):: frierson_Boltz      = 5.6734d-8
  real(r8):: frierson_C0         = 1.e7_R8
  real(r8):: frierson_Tmin       = 271._R8
  real(r8):: frierson_Tdlt       = 39._R8
  real(r8):: frierson_Twidth     = 26._R8
  real(r8):: frierson_WetDryCoef = 1._R8

  ! Global values stored for upward sweep
  !----------------------------------------
  real(r8),allocatable:: Tsfc_bc(:,:)
  real(r8),allocatable:: Qsfc_bc(:,:)
  real(r8),allocatable:: Eval_m (:,:,:)
  real(r8),allocatable:: Eval_e (:,:,:)
  real(r8),allocatable:: Fval_t (:,:,:)
  real(r8),allocatable:: Fval_q (:,:,:)
  real(r8),allocatable:: Fval_u (:,:,:)
  real(r8),allocatable:: Fval_v (:,:,:)

  ! Values that should be passed to surface for implicit calc.
  !-------------------------------------------------------------
  real(r8),allocatable:: Cstar  (:,:)
  real(r8),allocatable:: MU_a   (:,:) 
  real(r8),allocatable:: Estar_t(:,:)
  real(r8),allocatable:: Estar_q(:,:)
  real(r8),allocatable:: Estar_u(:,:)
  real(r8),allocatable:: Estar_v(:,:)
  real(r8),allocatable:: dFa_dTa(:,:)
  real(r8),allocatable:: dFa_dQa(:,:)
  real(r8),allocatable:: dFa_dUa(:,:)
  real(r8),allocatable:: dFa_dVa(:,:)

  ! Redundant surface values that should be calculated 
  ! by FLUX and passed to SOM (calculated after downward sweep)
  !-------------------------------------------------------------
  real(r8),allocatable:: Ft      (:,:)
  real(r8),allocatable:: Fq      (:,:)
  real(r8),allocatable:: Fu      (:,:)
  real(r8),allocatable:: Fv      (:,:)
  real(r8),allocatable:: dFt_dTs (:,:)
  real(r8),allocatable:: dFq_dTs (:,:)
  real(r8),allocatable:: dFup_dTs(:,:)

  ! Redundant surface values that should be 
  ! passed back from Flux calculation
  !-------------------------------------------------------------------
  real(r8),allocatable:: FN_t(:,:)
  real(r8),allocatable:: FN_q(:,:)
  real(r8),allocatable:: FN_u(:,:)
  real(r8),allocatable:: FN_v(:,:)
  real(r8),allocatable:: EN_t(:,:)
  real(r8),allocatable:: EN_q(:,:)

contains
  !==============================================================================
  subroutine frierson_register()
    ! 
    ! frierson_register: Register physics buffer values
    !=====================================================================

    call pbuf_add_field('PREC_PCW','physpkg',dtype_r8, (/pcols/),     prec_pcw_idx)
    call pbuf_add_field('PREC_DP' ,'physpkg',dtype_r8, (/pcols/),     prec_dp_idx )
    call pbuf_add_field('RELHUM'  ,'physpkg',dtype_r8, (/pcols,pver/),relhum_idx  )

    ! End Routine
    !-------------
    return
  end subroutine frierson_register
  !==============================================================================


  !==============================================================================
  subroutine frierson_readnl(nlfile)
    ! 
    ! frierson_readnl: Read in parameters controlling Frierson parameterizations.
    !=====================================================================
    use namelist_utils,only: find_group_name
    use units         ,only: getunit, freeunit
    ! 
    ! Passed Variables
    !------------------
    character(len=*),intent(in):: nlfile
    ! 
    ! Local Values
    !--------------
    integer:: ierr,unitn

    namelist /frierson_nl/ frierson_T0 , frierson_E0    , frierson_Erad    , frierson_Wind_min, &
                           frierson_Z0 , frierson_Ri_c  , frierson_Karman  , frierson_Fb      , &
                           frierson_Rs0, frierson_DeltaS, frierson_Tau_eqtr, frierson_Tau_pole, &
                           frierson_C0 , frierson_Tmin  , frierson_Tdlt    , frierson_Twidth  , &
                           frierson_LinFrac, frierson_Boltz, frierson_WetDryCoef

    ! Set default namelist values
    !-----------------------------
    frierson_T0         = 273.16_r8
    frierson_E0         = 610.78_r8
    frierson_Erad       = 6.376d6
    frierson_Wind_min   = 1.0d-5
    frierson_Z0         = 3.21d-5
    frierson_Ri_c       = 1.0_r8
    frierson_Karman     = 0.4_r8
    frierson_Fb         = 0.1_r8
    frierson_Rs0        = 938.4_r8
    frierson_DeltaS     = 1.4_r8
    frierson_Tau_eqtr   = 6.0_r8
    frierson_Tau_pole   = 1.5_r8
    frierson_LinFrac    = 0.1_r8
    frierson_Boltz      = 5.6734d-8
    frierson_C0         = 1.e7_R8
    frierson_Tmin       = 271._R8
    frierson_Tdlt       = 39._R8
    frierson_Twidth     = 26._R8
    frierson_WetDryCoef = 1._R8

    ! Read in namelist values
    !-------------------------
    if(masterproc) then
      unitn = getunit()
      open(unitn,file=trim(nlfile),status='old')
      call find_group_name(unitn,'frierson_nl',status=ierr)
      if(ierr.eq.0) then
        read(unitn,frierson_nl,iostat=ierr)
        if(ierr.ne.0) then
          call endrun('frierson_readnl:: ERROR reading namelist')
        endif
      endif
      close(unitn)
      call freeunit(unitn)
    endif

    ! Sanity Check namelist values
    !--------------------------------


    ! Broadcast namelist values
    !---------------------------
#ifdef SPMD
    call mpibcast(frierson_T0        , 1, mpir8 , 0, mpicom)
    call mpibcast(frierson_E0        , 1, mpir8 , 0, mpicom)
    call mpibcast(frierson_Erad      , 1, mpir8 , 0, mpicom)
    call mpibcast(frierson_Wind_min  , 1, mpir8 , 0, mpicom)
    call mpibcast(frierson_Z0        , 1, mpir8 , 0, mpicom)
    call mpibcast(frierson_Ri_c      , 1, mpir8 , 0, mpicom)
    call mpibcast(frierson_Karman    , 1, mpir8 , 0, mpicom)
    call mpibcast(frierson_Fb        , 1, mpir8 , 0, mpicom)
    call mpibcast(frierson_Rs0       , 1, mpir8 , 0, mpicom)
    call mpibcast(frierson_DeltaS    , 1, mpir8 , 0, mpicom)
    call mpibcast(frierson_Tau_eqtr  , 1, mpir8 , 0, mpicom)
    call mpibcast(frierson_Tau_pole  , 1, mpir8 , 0, mpicom)
    call mpibcast(frierson_LinFrac   , 1, mpir8 , 0, mpicom)
    call mpibcast(frierson_Boltz     , 1, mpir8 , 0, mpicom)
    call mpibcast(frierson_C0        , 1, mpir8 , 0, mpicom)
    call mpibcast(frierson_Tmin      , 1, mpir8 , 0, mpicom)
    call mpibcast(frierson_Tdlt      , 1, mpir8 , 0, mpicom)
    call mpibcast(frierson_Twidth    , 1, mpir8 , 0, mpicom)
    call mpibcast(frierson_WetDryCoef, 1, mpir8 , 0, mpicom)

#endif

    ! End Routine
    !-------------
    return
  end subroutine frierson_readnl
  !==============================================================================


  !==============================================================================
  subroutine frierson_init(phys_state,pbuf2d)
    !
    ! frierson_init: allocate space for global arrays and initialize values. 
    !                Add variables to history outputs
    !=====================================================================
    use physics_types, only: physics_state
    use error_messages,only: alloc_err
    use cam_history,   only: addfld, add_default,horiz_only
    use phys_grid,     only: get_ncols_p, get_rlat_p
    use frierson,      only: frierson_set_const
    ! 
    ! Passed Variables
    !------------------
    type(physics_state)      ,pointer:: phys_state(:)
    type(physics_buffer_desc),pointer:: pbuf2d    (:,:)
    !
    ! Local Values
    !---------------
    integer :: istat,lchnk,icol,ncol,ll
    real(r8):: adjusted_E0

    ! Initialize constants in Frierson module
    !------------------------------------------
    adjusted_E0 = frierson_WetDryCoef*frierson_E0
    call frierson_set_const(gravit,cappa,rair,cpair,latvap,rh2o,epsilo,rhoh2o,zvir,ps0,etamid, &
                            frierson_T0 ,adjusted_E0    ,frierson_Erad    ,frierson_Wind_min,  &
                            frierson_Z0 ,frierson_Ri_c  ,frierson_Karman  ,frierson_Fb      ,  &
                            frierson_Rs0,frierson_DeltaS,frierson_Tau_eqtr,frierson_Tau_pole,  &
                            frierson_LinFrac,frierson_Boltz)

    ! Add values for history output 
    !---------------------------------
    call addfld('QRS',(/'lev'/),'A','K/s','Temperature tendency associated with the '//           &
                                          'relaxation toward the equilibrium temperature profile' )
    call addfld('KVH' ,(/'ilev'/),'A','m2/s'   ,'Vertical diffusion diffusivities (heat/moisture)')
    call addfld('KVM' ,(/'ilev'/),'A','m2/s'   ,'Vertical diffusion diffusivities (momentum)'     )
    call addfld('VSE' ,(/'lev' /),'A','K'      ,'VSE: (Tv + gZ/Cp)'                                )
    call addfld('Zm'  ,(/'lev' /),'A','m'      ,'ATM Layer Heights use in PBL'                    )
    call addfld('DTV' ,(/'lev' /),'A','K/s'    ,'T vertical diffusion'                            )
    call addfld('DUV' ,(/'lev' /),'A','m/s2'   ,'U vertical diffusion'                            )
    call addfld('DVV' ,(/'lev' /),'A','m/s2'   ,'V vertical diffusion'                            )
    call addfld('VD01',(/'lev' /),'A','kg/kg/s','Q tendency (vertical diffusion)'                 )
    call addfld('Cdrag',horiz_only,'A','n/a'    ,'Surface Drag'                                   )
    call addfld('Z_pbl',horiz_only,'I','m'      ,'PBL Height'                                     )
    call addfld('Rf'  ,(/'lev' /),'I','n/a'     ,'Another Richardson number / Ri_c'                )

    call addfld('R_Fsolar', horiz_only, 'I','W/m2', 'SW Solar Flux'             )
    call addfld('R_Fup'   , horiz_only, 'I','W/m2', 'LW Upward Radiative Flux'  )
    call addfld('R_Fdown' , horiz_only, 'I','W/m2', 'LW Downward Radiative Flux')
    call addfld('R_SHflux', horiz_only, 'I','W/m2', 'Sensible Heat Flux'        )
    call addfld('R_LHflux', horiz_only, 'I','W/m2', 'Latent Heat Flux'          )
    call addfld('R_Fnet  ', horiz_only, 'I','W/m2', 'Net Radiative Flux'        )
    call addfld('R_Tsurf ', horiz_only, 'I','K'   , 'Surface Temperature'       )
    call addfld('R_Qsurf ', horiz_only, 'I','kg/kg', 'Surface Water Vapor'      )
    call addfld('R_Cdrag' , horiz_only, 'I','n/a'  , 'Surface Drag'             )

    call add_default('QRS'  ,1,' ')
    call add_default('KVH'  ,1,' ')
    call add_default('KVM'  ,1,' ')
    call add_default('VSE'  ,1,' ')
    call add_default('Zm'   ,1,' ')
    call add_default('DTV'  ,1,' ')
    call add_default('DUV'  ,1,' ')
    call add_default('DVV'  ,1,' ')
    call add_default('VD01' ,1,' ')
    call add_default('Cdrag',1,' ')
    call add_default('Z_pbl',1,' ')
    call add_default('Rf'   ,1,' ')

    ! Initialize physics buffer values
    !----------------------------------
    if(is_first_step()) then
       call pbuf_set_field(pbuf2d, prec_pcw_idx, 0._r8)
       call pbuf_set_field(pbuf2d, prec_dp_idx , 0._r8)
    endif

    ! Allocate Global arrays
    !-------------------------
    allocate(Tsurf (pcols,3,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'Frierson INIT','Tsurf' ,pcols*3*(endchunk-begchunk+1))
    allocate(Qsurf (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Qsurf' ,pcols*(endchunk-begchunk+1))
    allocate(Fsolar(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fsolar',pcols*(endchunk-begchunk+1))
    allocate(Fup   (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fup'   ,pcols*(endchunk-begchunk+1))
    allocate(Fdown (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fdown' ,pcols*(endchunk-begchunk+1))
    allocate(SHflux(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','SHflux',pcols*(endchunk-begchunk+1))
    allocate(LHflux(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','LHflux',pcols*(endchunk-begchunk+1))
    allocate(SWflux(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','SWflux',pcols*(endchunk-begchunk+1))
    allocate(TUflux(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','TUflux',pcols*(endchunk-begchunk+1))
    allocate(TVflux(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','TVflux',pcols*(endchunk-begchunk+1))
    allocate(Evap  (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Evap'  ,pcols*(endchunk-begchunk+1))
    allocate(Cd    (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Cd'    ,pcols*(endchunk-begchunk+1))
    allocate(clat  (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','clat'  ,pcols*(endchunk-begchunk+1))
    allocate(LAST  (begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','LAST'  ,(endchunk-begchunk+1))
    allocate(CURR  (begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','CURR'  ,(endchunk-begchunk+1))
    allocate(NEXT  (begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','NEXT'  ,(endchunk-begchunk+1))

    allocate(Fnet(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fnet',pcols*(endchunk-begchunk+1))

    ! Global values stored for upward sweep
    !----------------------------------------
    allocate(Tsfc_bc(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Tsfc_bc' ,pcols*(endchunk-begchunk+1))
    allocate(Qsfc_bc(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Qsfc_bc' ,pcols*(endchunk-begchunk+1))
    allocate(Eval_m (pcols,pver,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Eval_m' ,pcols*pver*(endchunk-begchunk+1))
    allocate(Eval_e (pcols,pver,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Eval_e' ,pcols*pver*(endchunk-begchunk+1))
    allocate(Fval_t (pcols,pver,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fval_t' ,pcols*pver*(endchunk-begchunk+1))
    allocate(Fval_q (pcols,pver,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fval_q' ,pcols*pver*(endchunk-begchunk+1))
    allocate(Fval_u (pcols,pver,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fval_u' ,pcols*pver*(endchunk-begchunk+1))
    allocate(Fval_v (pcols,pver,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fval_v' ,pcols*pver*(endchunk-begchunk+1))

    ! Values that should be passed to surface for implicit calc.
    !-------------------------------------------------------------
    allocate(Cstar  (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Cstar' ,pcols*(endchunk-begchunk+1))
    allocate(MU_a   (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','MU_a' ,pcols*(endchunk-begchunk+1))
    allocate(Estar_t(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Estar_t' ,pcols*(endchunk-begchunk+1))
    allocate(Estar_q(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Estar_q' ,pcols*(endchunk-begchunk+1))
    allocate(Estar_u(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Estar_u' ,pcols*(endchunk-begchunk+1))
    allocate(Estar_v(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Estar_v' ,pcols*(endchunk-begchunk+1))
    allocate(dFa_dTa(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','dFa_dTa' ,pcols*(endchunk-begchunk+1))
    allocate(dFa_dQa(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','dFa_dQa' ,pcols*(endchunk-begchunk+1))
    allocate(dFa_dUa(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','dFa_dUa' ,pcols*(endchunk-begchunk+1))
    allocate(dFa_dVa(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','dFa_dVa' ,pcols*(endchunk-begchunk+1))

    ! Redundant surface values that should be calculated 
    ! by FLUX and passed to SOM (calculated after downward sweep)
    !-------------------------------------------------------------
    allocate(Ft      (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Ft' ,pcols*(endchunk-begchunk+1))
    allocate(Fq      (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fq' ,pcols*(endchunk-begchunk+1))
    allocate(Fu      (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fu' ,pcols*(endchunk-begchunk+1))
    allocate(Fv      (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fv' ,pcols*(endchunk-begchunk+1))
    allocate(dFt_dTs (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','dFt_dTs' ,pcols*(endchunk-begchunk+1))
    allocate(dFq_dTs (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','dFq_dTs' ,pcols*(endchunk-begchunk+1))
    allocate(dFup_dTs(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','dFup_dTs' ,pcols*(endchunk-begchunk+1))

    ! Redundant surface values that should be 
    ! passed back from Flux calculation
    !-------------------------------------------------------------------
    allocate(FN_t(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','FN_t' ,pcols*(endchunk-begchunk+1))
    allocate(FN_q(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','FN_q' ,pcols*(endchunk-begchunk+1))
    allocate(FN_u(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','FN_u' ,pcols*(endchunk-begchunk+1))
    allocate(FN_v(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','FN_v' ,pcols*(endchunk-begchunk+1))
    allocate(EN_t(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','EN_t' ,pcols*(endchunk-begchunk+1))
    allocate(EN_q(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','EN_q' ,pcols*(endchunk-begchunk+1))

    ! Initialize time indices and latitudes
    !----------------------------------------
    do lchnk = begchunk,endchunk
      ncol = get_ncols_p(lchnk)
      do icol = 1,ncol
        clat(icol,lchnk) = get_rlat_p(lchnk,icol)
      end do
      LAST(lchnk) = 1  
      CURR(lchnk) = 2
      NEXT(lchnk) = 3
    end do

    ! Initialize (DIAGNOSTIC) Surface temperatures
    !-----------------------------------------------------------------------
    do lchnk = begchunk,endchunk
      ncol = get_ncols_p(lchnk)

      ! phys_state PS values are Inf at this point.... HUH???
      ! Set to reference values for initialization
      !------------------------------------------------------------
      phys_state(lchnk)%ps(:ncol) = ps0

      call frierson_surface_init(ncol,         clat(:ncol,lchnk),             &
                               phys_state(lchnk)%ps(:ncol),                   &
                                              Tsurf(:ncol,LAST(lchnk),lchnk), &
                                              Qsurf(:ncol,lchnk)              )
      Tsurf(:,CURR(lchnk),lchnk) = Tsurf(:,LAST(lchnk),lchnk)
      Tsurf(:,NEXT(lchnk),lchnk) = Tsurf(:,LAST(lchnk),lchnk)
    end do
      
    ! Initialize radition and flux values to 0.0  (Add Init from restart file???)
    !---------------------------------------------------------------------------
    do lchnk = begchunk,endchunk
      Fsolar(:,lchnk) = 0._r8
      Fup   (:,lchnk) = 0._r8
      Fdown (:,lchnk) = 0._r8
      SHflux(:,lchnk) = 0._r8
      LHflux(:,lchnk) = 0._r8
      SWflux(:,lchnk) = 0._r8
      TUflux(:,lchnk) = 0._r8
      TVflux(:,lchnk) = 0._r8
      Evap  (:,lchnk) = 0._r8
      Cd    (:,lchnk) = 0._r8
      Fnet  (:,lchnk) = 0._r8
    end do

    ! Informational Output
    !----------------------
    write(iulog,*) ' '
    write(iulog,*) '-----------------------------------------------------------'
    write(iulog,*) '  FRIERSON MODULE INITIALIZED WITH THE FOLLOWING SETTINGS: '
    write(iulog,*) '-----------------------------------------------------------'
    write(iulog,*) 'FRIERSON: gravit='    , gravit
    write(iulog,*) 'FRIERSON: cappa='     , cappa
    write(iulog,*) 'FRIERSON: rair ='     , rair
    write(iulog,*) 'FRIERSON: cpair='     , cpair
    write(iulog,*) 'FRIERSON: latvap='    , latvap
    write(iulog,*) 'FRIERSON: rh2o='      , rh2o
    write(iulog,*) 'FRIERSON: epsilo='    , epsilo
    write(iulog,*) 'FRIERSON: rhoh2o='    , rhoh2o
    write(iulog,*) 'FRIERSON: zvir='      , zvir
    write(iulog,*) 'FRIERSON: ps0='       , ps0
    write(iulog,*) 'FRIERSON: etamid='    , etamid
    write(iulog,*) 'FRIERSON: T0='        , frierson_T0
    write(iulog,*) 'FRIERSON: E0='        , frierson_E0
    write(iulog,*) 'FRIERSON: Erad='      , frierson_Erad
    write(iulog,*) 'FRIERSON: Wind_min='  , frierson_Wind_min
    write(iulog,*) 'FRIERSON: Z0='        , frierson_Z0
    write(iulog,*) 'FRIERSON: Ri_c='      , frierson_Ri_c
    write(iulog,*) 'FRIERSON: Karman='    , frierson_Karman
    write(iulog,*) 'FRIERSON: Fb='        , frierson_Fb
    write(iulog,*) 'FRIERSON: Rs0='       , frierson_Rs0
    write(iulog,*) 'FRIERSON: DeltaS='    , frierson_DeltaS
    write(iulog,*) 'FRIERSON: Tau_eqtr='  , frierson_Tau_eqtr
    write(iulog,*) 'FRIERSON: Tau_pole='  , frierson_Tau_pole
    write(iulog,*) 'FRIERSON: LinFrac='   , frierson_LinFrac
    write(iulog,*) 'FRIERSON: Boltz='     , frierson_Boltz
    write(iulog,*) 'FRIERSON: C0='        , frierson_C0
    write(iulog,*) 'FRIERSON: Tmin='      , frierson_Tmin
    write(iulog,*) 'FRIERSON: Tdlt='      , frierson_Tdlt
    write(iulog,*) 'FRIERSON: Twidth='    , frierson_Twidth
    write(iulog,*) 'FRIERSON: WetDryCoef=', frierson_WetDryCoef
    write(iulog,*) ' '

    ! End Routine
    !--------------
    return
  end subroutine frierson_init
  !==============================================================================


  !==============================================================================
  subroutine frierson_convection_tend(state, ptend, ztodt, pbuf)
    !
    ! frierson_convection_tend: Run the selected process to compute precipitation 
    !                           due to convection.
    !=====================================================================
    use physics_types,only: physics_state, physics_ptend
    use physics_types,only: physics_ptend_init
    use physconst,    only: cpair
    use frierson,     only: frierson_convection_NONE,frierson_convection
    !
    ! Passed Variables
    !------------------
    type(physics_state)      ,intent(inout):: state
    real(r8)                 ,intent(in)   :: ztodt 
    type(physics_ptend)      ,intent(out)  :: ptend 
    type(physics_buffer_desc),pointer      :: pbuf(:)
    !
    ! Local Values
    !-----------------
    real(r8),pointer:: prec_dp (:)          ! convective precip
    real(r8),pointer:: relhum  (:,:)
    real(r8)        :: T (state%ncol, pver) ! T temporary
    real(r8)        :: qv(state%ncol, pver) ! Q temporary
    logical         :: lq(pcnst)            ! Calc tendencies?
    integer         :: lchnk                ! chunk identifier
    integer         :: ncol                 ! number of atmospheric columns
    integer         :: k

    ! Set local copies of values
    !---------------------------------
    lchnk       = state%lchnk
    ncol        = state%ncol
    T (:ncol,:) = state%T(:ncol,:)
    qv(:ncol,:) = state%Q(:ncol,:,1)

    ! initialize individual parameterization tendencies
    !---------------------------------------------------
    lq    = .false.
    lq(1) = .true.
    call physics_ptend_init(ptend, state%psetcols, 'Frierson convection', &
                                ls=.true., lu=.true., lv=.true., lq=lq)

    ! Get values from the physics buffer
    !------------------------------------
    call pbuf_get_field(pbuf,prec_dp_idx ,prec_dp )
    call pbuf_get_field(pbuf,  relhum_idx,relhum  )

    ! Call the Selected convection routine
    !--------------------------------------------------------
    if(CONVECTION_OPT.eq.CONVECTION_FRIERSON) then
      call frierson_convection(ncol,pver,ztodt,state%pmid(:ncol,:), &
                                               state%pdel(:ncol,:), &
                                                        T(:ncol,:), &
                                                       qv(:ncol,:), &
                                                   relhum(:ncol,:), &
                                                  prec_dp(:ncol)    )
    else ! (CONVECTION_OPT.eq.CONVECTION_NONE) then
      call frierson_convection_NONE(ncol,pver,ztodt,state%pmid(:ncol,:), &
                                                    state%pdel(:ncol,:), &
                                                             T(:ncol,:), &
                                                            qv(:ncol,:), &
                                                        relhum(:ncol,:), &
                                                       prec_dp(:ncol)    )
    endif

    ! Back out temperature and specific humidity 
    ! tendencies from updated fields
    !--------------------------------------------
    do k = 1, pver
      ptend%s(:ncol,k)   = (T (:,k)-state%T(:ncol,k)  )/ztodt*cpair
      ptend%q(:ncol,k,1) = (qv(:,k)-state%q(:ncol,k,1))/ztodt
    end do

    ! End Routine
    !--------------
    return
  end subroutine frierson_convection_tend
  !==============================================================================


  !==============================================================================
  subroutine frierson_condensate_tend(state, ptend, ztodt, pbuf)
    !
    ! frierson_condensate_tend: Run the selected process to compute precipitation 
    !                           due to large scale condensation.
    !=====================================================================
    use physics_types,only: physics_state, physics_ptend
    use physics_types,only: physics_ptend_init
    use physconst,    only: cpair
    use frierson,     only: frierson_condensate_NONE,frierson_condensate, &
                            frierson_condensate_TJ16
    !
    ! Passed Variables
    !------------------
    type(physics_state)      ,intent(inout):: state
    real(r8)                 ,intent(in)   :: ztodt 
    type(physics_ptend)      ,intent(out)  :: ptend 
    type(physics_buffer_desc),pointer      :: pbuf(:)
    !
    ! Local Values
    !-----------------
    real(r8),pointer:: prec_pcw(:)          ! large scale precip
    real(r8),pointer:: relhum  (:,:)
    real(r8)        :: T (state%ncol, pver) ! T temporary
    real(r8)        :: qv(state%ncol, pver) ! Q temporary
    logical         :: lq(pcnst)            ! Calc tendencies?
    integer         :: lchnk                ! chunk identifier
    integer         :: ncol                 ! number of atmospheric columns
    integer         :: k

    ! Set local copies of values
    !---------------------------------
    lchnk       = state%lchnk
    ncol        = state%ncol
    T (:ncol,:) = state%T(:ncol,:)
    qv(:ncol,:) = state%Q(:ncol,:,1)

    ! initialize individual parameterization tendencies
    !---------------------------------------------------
    lq    = .false.
    lq(1) = .true.
    call physics_ptend_init(ptend, state%psetcols, 'Frierson condensate', &
                                ls=.true., lu=.true., lv=.true., lq=lq)

    ! Get values from the physics buffer
    !------------------------------------
    call pbuf_get_field(pbuf,prec_pcw_idx,prec_pcw)
    call pbuf_get_field(pbuf,  relhum_idx,relhum  )

    ! Call the Selected condensation routine  ~~DEVO style~~
    !--------------------------------------------------------
    if(CONDENSATE_OPT.eq.CONDENSATE_FRIERSON) then
      call frierson_condensate(ncol,pver,ztodt,state%pmid(:ncol,:), &
                                               state%pdel(:ncol,:), &
                                                        T(:ncol,:), &
                                                       qv(:ncol,:), &
                                                   relhum(:ncol,:), &
                                                 prec_pcw(:ncol)    )
    elseif(CONDENSATE_OPT.eq.CONDENSATE_FRIERSON_TJ16) then
      call frierson_condensate_TJ16(ncol,pver,ztodt,state%pmid(:ncol,:), &
                                                    state%pdel(:ncol,:), &
                                                             T(:ncol,:), &
                                                            qv(:ncol,:), &
                                                        relhum(:ncol,:), &
                                                      prec_pcw(:ncol)    )
    else ! (CONDENSATE_OPT.eq.CONDENSATE_NONE) then
      call frierson_condensate_NONE(ncol,pver,ztodt,state%pmid(:ncol,:), &
                                                    state%pdel(:ncol,:), &
                                                             T(:ncol,:), &
                                                            qv(:ncol,:), &
                                                        relhum(:ncol,:), &
                                                      prec_pcw(:ncol)    )
    endif

    ! Back out temperature and specific humidity 
    ! tendencies from updated fields
    !--------------------------------------------
    do k = 1, pver
      ptend%s(:ncol,k)   = (T (:,k)-state%T(:ncol,k)  )/ztodt*cpair
      ptend%q(:ncol,k,1) = (qv(:,k)-state%q(:ncol,k,1))/ztodt
    end do

    ! End Routine
    !--------------
    return
  end subroutine frierson_condensate_tend
  !==============================================================================


  !============================================================================
  subroutine frierson_pbl_init(state, ptend, ztodt, cam_in)
    !
    ! frierson_pbl_init: Run the selected PBL process.
    !=========================================================================
    use physics_types,only: physics_state, physics_ptend
    use physics_types,only: physics_ptend_init
    use physconst,    only: cpair
    use phys_grid,    only: get_rlat_all_p
    use frierson,     only: frierson_pbl_down_impl, frierson_pbl_flux_impl, frierson_pbl_surface_impl, &
                            frierson_pbl_flux_expsfc
    !
    ! Passed Variables
    !-------------------
    type(physics_state),intent(in)   :: state
    real(r8),           intent(in)   :: ztodt
    type(physics_ptend),intent(out)  :: ptend
    type(cam_in_t),     intent(inout):: cam_in
    !
    ! Local Values
    !----------------
    real(r8) :: Km   (state%ncol,pver+1) ! Eddy diffusivity at layer interfaces (m2/s)
    real(r8) :: Ke   (state%ncol,pver+1) ! Eddy diffusivity at layer interfaces (m2/s)
    real(r8) :: VSE  (state%ncol,pver)   ! Dry Static Energy divided by Cp (K)
    real(r8) :: Rf   (state%ncol,pver)   ! 
    real(r8) :: Z_pbl(state%ncol)        ! 
    real(r8) :: Cdrag(state%ncol)        ! Cdrag coef from surface calculation
    real(r8) :: Ws_a (state%ncol)        ! 
    real(r8) :: rho_a(state%ncol)        ! 
    !
    ! Local work values needed for redundant calcs
    !----------------------------------------------
    real(r8):: dFt_dTa(state%ncol)
    real(r8):: dFq_dQa(state%ncol)
    real(r8):: dFu_dUa(state%ncol)
    real(r8):: dFv_dVa(state%ncol)

    integer:: lchnk   ! chunk identifier
    integer:: ncol    ! number of atmospheric columns

    real(r8) :: Tsfc      (state%ncol)        ! Surface T 
    real(r8) :: Qsfc      (state%ncol)        ! Surface Q (saturated)

    integer :: ii_diag
    real(r8):: rtd_diag

    ! Set local copies of values
    !---------------------------------
    lchnk = state%lchnk
    ncol  = state%ncol

    ! Initialize PBL with downward sweep of implicit method.
    !-----------------------------------------------------
    if((PBL_OPT.eq.FRIERSON_IMPL_DIAG   ).or. &
       (PBL_OPT.eq.FRIERSON_IMPL_DIAGC  ).or. &
       (PBL_OPT.eq.FRIERSON_EXPSFC_DIAG ).or. &
       (PBL_OPT.eq.FRIERSON_EXPSFC_DIAGC)     ) then
      continue  ! Nothig to do here
    elseif(PBL_OPT.eq.FRIERSON_IMPL_CONFIG1) then
      Tsfc_bc(:ncol,lchnk) = Tsurf(:ncol,CURR(lchnk),lchnk)
      Qsfc_bc(:ncol,lchnk) = Qsurf(:ncol,lchnk)
      call frierson_pbl_down_impl(ncol, pver, ztodt, state%pmid(:ncol,:),       &
                                                     state%pint(:ncol,:),       &
                                                     state%zm  (:ncol,:),       &
                                                     state%zi  (:ncol,1:pver),  &
                                                     state%ps  (:ncol),         &
                                                        Tsfc_bc(:ncol,lchnk),   &
                                                        Qsfc_bc(:ncol,lchnk),   &
                                                     state%T   (:ncol,:),       &
                                                     state%U   (:ncol,:),       &
                                                     state%V   (:ncol,:),       &
                                                     state%Q   (:ncol,:,1),     &
                                                         Eval_m(:ncol,:,lchnk), &
                                                         Eval_e(:ncol,:,lchnk), &
                                                         Fval_t(:ncol,:,lchnk), &
                                                         Fval_q(:ncol,:,lchnk), &
                                                         Fval_u(:ncol,:,lchnk), &
                                                         Fval_v(:ncol,:,lchnk), &
                                                          Cdrag(:ncol),         &
                                                           MU_a(:ncol,lchnk),   &
                                                        Estar_t(:ncol,lchnk),   &
                                                        Estar_q(:ncol,lchnk),   &
                                                        Estar_u(:ncol,lchnk),   &
                                                        Estar_v(:ncol,lchnk),   &
                                                        dFa_dTa(:ncol,lchnk),   &
                                                        dFa_dQa(:ncol,lchnk),   &
                                                        dFa_dUa(:ncol,lchnk),   &
                                                        dFa_dVa(:ncol,lchnk),   &
                                                             Km(:ncol,:),       &
                                                             Ke(:ncol,:),       &
                                                            VSE(:ncol,:),       &
                                                          Z_pbl(:ncol),         &
                                                           Ws_a(:ncol),         &
                                                          rho_a(:ncol),         &
                                                             Rf(:ncol,:)        )
      Cstar(:ncol,lchnk) = rho_a(:ncol)*Cdrag(:ncol)*Ws_a(:ncol)
      call frierson_pbl_flux_impl(ncol, ztodt, state%pmid(:ncol,pver),         &
                                                        state%T(:ncol,pver),   &
                                                        state%U(:ncol,pver),   &
                                                        state%V(:ncol,pver),   &
                                                        state%Q(:ncol,pver,1), &
                                                       state%ps(:ncol),        &
                                                        Tsfc_bc(:ncol,lchnk),  &
                                                        Qsfc_bc(:ncol,lchnk),  &
                                                          Cstar(:ncol,lchnk),  &
                                                           MU_a(:ncol,lchnk),  &
                                                        Estar_t(:ncol,lchnk),  &
                                                        Estar_q(:ncol,lchnk),  &
                                                        Estar_u(:ncol,lchnk),  &
                                                        Estar_v(:ncol,lchnk),  &
                                                        dFa_dTa(:ncol,lchnk),  &
                                                        dFa_dQa(:ncol,lchnk),  &
                                                        dFa_dUa(:ncol,lchnk),  &
                                                        dFa_dVa(:ncol,lchnk),  &
                                                             Ft(:ncol,lchnk),  &
                                                             Fq(:ncol,lchnk),  &
                                                             Fu(:ncol,lchnk),  &
                                                             Fv(:ncol,lchnk),  &
                                                            Fup(:ncol,lchnk),  &
                                                        dFt_dTs(:ncol,lchnk),  &
                                                        dFq_dTs(:ncol,lchnk),  &
                                                       dFup_dTs(:ncol,lchnk),  &
                                                           FN_t(:ncol,lchnk),  &
                                                           FN_q(:ncol,lchnk),  &
                                                           FN_u(:ncol,lchnk),  &
                                                           FN_v(:ncol,lchnk),  &
                                                           EN_t(:ncol,lchnk),  &
                                                           EN_q(:ncol,lchnk)   )
    elseif(PBL_OPT.eq.FRIERSON_EXPSFC_CONFIG1) then
      Tsfc_bc(:ncol,lchnk) = Tsurf(:ncol,CURR(lchnk),lchnk)
      Qsfc_bc(:ncol,lchnk) = Qsurf(:ncol,lchnk)
      call frierson_pbl_down_impl(ncol, pver, ztodt, state%pmid(:ncol,:),       &
                                                     state%pint(:ncol,:),       &
                                                     state%zm  (:ncol,:),       &
                                                     state%zi  (:ncol,1:pver),  &
                                                     state%ps  (:ncol),         &
                                                        Tsfc_bc(:ncol,lchnk),   &
                                                        Qsfc_bc(:ncol,lchnk),   &
                                                     state%T   (:ncol,:),       &
                                                     state%U   (:ncol,:),       &
                                                     state%V   (:ncol,:),       &
                                                     state%Q   (:ncol,:,1),     &
                                                         Eval_m(:ncol,:,lchnk), &
                                                         Eval_e(:ncol,:,lchnk), &
                                                         Fval_t(:ncol,:,lchnk), &
                                                         Fval_q(:ncol,:,lchnk), &
                                                         Fval_u(:ncol,:,lchnk), &
                                                         Fval_v(:ncol,:,lchnk), &
                                                          Cdrag(:ncol),         &
                                                           MU_a(:ncol,lchnk),   &
                                                        Estar_t(:ncol,lchnk),   &
                                                        Estar_q(:ncol,lchnk),   &
                                                        Estar_u(:ncol,lchnk),   &
                                                        Estar_v(:ncol,lchnk),   &
                                                        dFa_dTa(:ncol,lchnk),   &
                                                        dFa_dQa(:ncol,lchnk),   &
                                                        dFa_dUa(:ncol,lchnk),   &
                                                        dFa_dVa(:ncol,lchnk),   &
                                                             Km(:ncol,:),       &
                                                             Ke(:ncol,:),       &
                                                            VSE(:ncol,:),       &
                                                          Z_pbl(:ncol),         &
                                                           Ws_a(:ncol),         &
                                                          rho_a(:ncol),         &
                                                             Rf(:ncol,:)        )

      call frierson_pbl_flux_expsfc(ncol, ztodt, state%zm(:ncol,pver),   &
                                               state%pmid(:ncol,pver),   &
                                                  state%T(:ncol,pver),   &
                                                  state%U(:ncol,pver),   &
                                                  state%V(:ncol,pver),   &
                                                  state%Q(:ncol,pver,1), &
                                                 state%ps(:ncol),        &
                                                  Tsfc_bc(:ncol,lchnk),  &
                                                  Qsfc_bc(:ncol,lchnk),  &
                                                       Ft(:ncol,lchnk),  &
                                                       Fq(:ncol,lchnk),  &
                                                       Fu(:ncol,lchnk),  &
                                                       Fv(:ncol,lchnk),  &
                                                      Fup(:ncol,lchnk)   )

      ! Incorporate surface fluxes into implicit scheme, then
      ! update flux values and derivatives
      !------------------------------------------
      dFt_dTa(:ncol) = 0._r8
      dFq_dQa(:ncol) = 0._r8
      dFu_dUa(:ncol) = 0._r8
      dFv_dVa(:ncol) = 0._r8

      FN_u(:ncol,lchnk) =  (Estar_u(:ncol,lchnk) + MU_a(:ncol,lchnk)*Fu(:ncol,lchnk))    &
                          /(1._r8-MU_a(:ncol,lchnk)*(dFa_dUa(:ncol,lchnk)+dFu_dUa(:ncol)))
      FN_v(:ncol,lchnk) =  (Estar_v(:ncol,lchnk) + MU_a(:ncol,lchnk)*Fv(:ncol,lchnk))    &
                          /(1._r8-MU_a(:ncol,lchnk)*(dFa_dVa(:ncol,lchnk)+dFv_dVa(:ncol)))
      FN_t(:ncol,lchnk) =  (Estar_t(:ncol,lchnk) + MU_a(:ncol,lchnk)*Ft(:ncol,lchnk))    &
                          /(1._r8-MU_a(:ncol,lchnk)*(dFa_dTa(:ncol,lchnk)+dFt_dTa(:ncol)))
      FN_q(:ncol,lchnk) =  (Estar_q(:ncol,lchnk) + MU_a(:ncol,lchnk)*Fq(:ncol,lchnk))    &
                          /(1._r8-MU_a(:ncol,lchnk)*(dFa_dQa(:ncol,lchnk)+dFq_dQa(:ncol)))

    elseif(PBL_OPT.eq.FRIERSON_EXPSFC_CONFIG2) then
      Tsfc_bc(:ncol,lchnk)     = cam_in%ts (:ncol)
      Qsfc_bc(:ncol,lchnk)     = cam_in%ssq(:ncol)
      call frierson_pbl_down_impl(ncol, pver, ztodt, state%pmid(:ncol,:),       &
                                                     state%pint(:ncol,:),       &
                                                     state%zm  (:ncol,:),       &
                                                     state%zi  (:ncol,1:pver),  &
                                                     state%ps  (:ncol),         &
                                                        Tsfc_bc(:ncol,lchnk),   &
                                                        Qsfc_bc(:ncol,lchnk),   &
                                                     state%T   (:ncol,:),       &
                                                     state%U   (:ncol,:),       &
                                                     state%V   (:ncol,:),       &
                                                     state%Q   (:ncol,:,1),     &
                                                         Eval_m(:ncol,:,lchnk), &
                                                         Eval_e(:ncol,:,lchnk), &
                                                         Fval_t(:ncol,:,lchnk), &
                                                         Fval_q(:ncol,:,lchnk), &
                                                         Fval_u(:ncol,:,lchnk), &
                                                         Fval_v(:ncol,:,lchnk), &
                                                          Cdrag(:ncol),         &
                                                           MU_a(:ncol,lchnk),   &
                                                        Estar_t(:ncol,lchnk),   &
                                                        Estar_q(:ncol,lchnk),   &
                                                        Estar_u(:ncol,lchnk),   &
                                                        Estar_v(:ncol,lchnk),   &
                                                        dFa_dTa(:ncol,lchnk),   &
                                                        dFa_dQa(:ncol,lchnk),   &
                                                        dFa_dUa(:ncol,lchnk),   &
                                                        dFa_dVa(:ncol,lchnk),   &
                                                             Km(:ncol,:),       &
                                                             Ke(:ncol,:),       &
                                                            VSE(:ncol,:),       &
                                                          Z_pbl(:ncol),         &
                                                           Ws_a(:ncol),         &
                                                          rho_a(:ncol),         &
                                                             Rf(:ncol,:)        )

      call frierson_pbl_flux_expsfc(ncol, ztodt, state%zm(:ncol,pver),   &
                                               state%pmid(:ncol,pver),   &
                                                  state%T(:ncol,pver),   &
                                                  state%U(:ncol,pver),   &
                                                  state%V(:ncol,pver),   &
                                                  state%Q(:ncol,pver,1), &
                                                 state%ps(:ncol),        &
                                                  Tsfc_bc(:ncol,lchnk),  &
                                                  Qsfc_bc(:ncol,lchnk),  &
                                                       Ft(:ncol,lchnk),  &
                                                       Fq(:ncol,lchnk),  &
                                                       Fu(:ncol,lchnk),  &
                                                       Fv(:ncol,lchnk),  &
                                                      Fup(:ncol,lchnk)   )

      ! Incorporate surface fluxes into implicit scheme, then
      ! update flux values and derivatives
      !------------------------------------------
      dFt_dTa(:ncol) = 0._r8
      dFq_dQa(:ncol) = 0._r8
      dFu_dUa(:ncol) = 0._r8
      dFv_dVa(:ncol) = 0._r8

      FN_u(:ncol,lchnk) =  (Estar_u(:ncol,lchnk) + MU_a(:ncol,lchnk)*Fu(:ncol,lchnk))    &
                          /(1._r8-MU_a(:ncol,lchnk)*(dFa_dUa(:ncol,lchnk)+dFu_dUa(:ncol)))
      FN_v(:ncol,lchnk) =  (Estar_v(:ncol,lchnk) + MU_a(:ncol,lchnk)*Fv(:ncol,lchnk))    &
                          /(1._r8-MU_a(:ncol,lchnk)*(dFa_dVa(:ncol,lchnk)+dFv_dVa(:ncol)))
      FN_t(:ncol,lchnk) =  (Estar_t(:ncol,lchnk) + MU_a(:ncol,lchnk)*Ft(:ncol,lchnk))    &
                          /(1._r8-MU_a(:ncol,lchnk)*(dFa_dTa(:ncol,lchnk)+dFt_dTa(:ncol)))
      FN_q(:ncol,lchnk) =  (Estar_q(:ncol,lchnk) + MU_a(:ncol,lchnk)*Fq(:ncol,lchnk))    &
                          /(1._r8-MU_a(:ncol,lchnk)*(dFa_dQa(:ncol,lchnk)+dFq_dQa(:ncol)))

    elseif(PBL_OPT.eq.FRIERSON_EXPSFC        ) then
      Tsfc_bc(:ncol,lchnk)     = cam_in%ts (:ncol)
      Qsfc_bc(:ncol,lchnk)     = cam_in%ssq(:ncol)
      call frierson_pbl_down_impl(ncol, pver, ztodt, state%pmid(:ncol,:),       &
                                                     state%pint(:ncol,:),       &
                                                     state%zm  (:ncol,:),       &
                                                     state%zi  (:ncol,1:pver),  &
                                                     state%ps  (:ncol),         &
                                                        Tsfc_bc(:ncol,lchnk),   &
                                                        Qsfc_bc(:ncol,lchnk),   &
                                                     state%T   (:ncol,:),       &
                                                     state%U   (:ncol,:),       &
                                                     state%V   (:ncol,:),       &
                                                     state%Q   (:ncol,:,1),     &
                                                         Eval_m(:ncol,:,lchnk), &
                                                         Eval_e(:ncol,:,lchnk), &
                                                         Fval_t(:ncol,:,lchnk), &
                                                         Fval_q(:ncol,:,lchnk), &
                                                         Fval_u(:ncol,:,lchnk), &
                                                         Fval_v(:ncol,:,lchnk), &
                                                          Cdrag(:ncol),         &
                                                           MU_a(:ncol,lchnk),   &
                                                        Estar_t(:ncol,lchnk),   &
                                                        Estar_q(:ncol,lchnk),   &
                                                        Estar_u(:ncol,lchnk),   &
                                                        Estar_v(:ncol,lchnk),   &
                                                        dFa_dTa(:ncol,lchnk),   &
                                                        dFa_dQa(:ncol,lchnk),   &
                                                        dFa_dUa(:ncol,lchnk),   &
                                                        dFa_dVa(:ncol,lchnk),   &
                                                             Km(:ncol,:),       &
                                                             Ke(:ncol,:),       &
                                                            VSE(:ncol,:),       &
                                                          Z_pbl(:ncol),         &
                                                           Ws_a(:ncol),         &
                                                          rho_a(:ncol),         &
                                                             Rf(:ncol,:)        )

      call frierson_pbl_flux_expsfc(ncol, ztodt, state%zm(:ncol,pver),   &
                                               state%pmid(:ncol,pver),   &
                                                  state%T(:ncol,pver),   &
                                                  state%U(:ncol,pver),   &
                                                  state%V(:ncol,pver),   &
                                                  state%Q(:ncol,pver,1), &
                                                 state%ps(:ncol),        &
                                                  Tsfc_bc(:ncol,lchnk),  &
                                                  Qsfc_bc(:ncol,lchnk),  &
                                                       Ft(:ncol,lchnk),  &
                                                       Fq(:ncol,lchnk),  &
                                                       Fu(:ncol,lchnk),  &
                                                       Fv(:ncol,lchnk),  &
                                                      Fup(:ncol,lchnk)   )

      ! Incorporate surface fluxes into implicit scheme, then
      ! update flux values and derivatives
      !------------------------------------------
      dFt_dTa(:ncol) = 0._r8
      dFq_dQa(:ncol) = 0._r8
      dFu_dUa(:ncol) = 0._r8
      dFv_dVa(:ncol) = 0._r8

      FN_u(:ncol,lchnk) =  (Estar_u(:ncol,lchnk) + MU_a(:ncol,lchnk)*Fu(:ncol,lchnk))    &
                          /(1._r8-MU_a(:ncol,lchnk)*(dFa_dUa(:ncol,lchnk)+dFu_dUa(:ncol)))
      FN_v(:ncol,lchnk) =  (Estar_v(:ncol,lchnk) + MU_a(:ncol,lchnk)*Fv(:ncol,lchnk))    &
                          /(1._r8-MU_a(:ncol,lchnk)*(dFa_dVa(:ncol,lchnk)+dFv_dVa(:ncol)))
      FN_t(:ncol,lchnk) =  (Estar_t(:ncol,lchnk) + MU_a(:ncol,lchnk)*Ft(:ncol,lchnk))    &
                          /(1._r8-MU_a(:ncol,lchnk)*(dFa_dTa(:ncol,lchnk)+dFt_dTa(:ncol)))
      FN_q(:ncol,lchnk) =  (Estar_q(:ncol,lchnk) + MU_a(:ncol,lchnk)*Fq(:ncol,lchnk))    &
                          /(1._r8-MU_a(:ncol,lchnk)*(dFa_dQa(:ncol,lchnk)+dFq_dQa(:ncol)))

    else

      !** ERROR STOP

    endif
 
    ! Some values that need to be passed to surface for implicit calc.
    !----------------------------------------------------------------
!    cam_out%Cstar(ncol)
!    cam_out%MU_a (ncol) 
!    cam_out%Estar_t(ncol)
!    cam_out%Estar_q(ncol)
!    cam_out%Estar_u(ncol)
!    cam_out%Estar_v(ncol)
!    cam_out%dFa_dTa(ncol)
!    cam_out%dFa_dQa(ncol)
!    cam_out%dFa_dUa(ncol)
!    cam_out%dFa_dVa(ncol)

   !=================================================================
   ! Compute values that should be passed back from surface routines:
   !=================================================================

    ! End Routine
    !--------------
    return
  end subroutine frierson_pbl_init
  !============================================================================


  !============================================================================
  subroutine frierson_pbl_tend(state, ptend, ztodt, cam_in)
    !
    ! frierson_pbl_tend: Run the selected PBL process.
    !=========================================================================
    use physics_types,only: physics_state, physics_ptend
    use physics_types,only: physics_ptend_init
    use physconst,    only: cpair
    use phys_grid,    only: get_rlat_all_p
    use frierson,     only: frierson_pbl_impl_diag, frierson_pbl_expsfc_diag,    &
                            frierson_pbl_flux_impl, frierson_pbl_surface_impl,   &
                            frierson_pbl_up_impl,   frierson_pbl_impl_diagC,     &
                            frierson_pbl_expsfc_diagC, frierson_pbl_flux_expsfc, &
                            frierson_pbl_surface_expsfc, frierson_pbl_up_expsfc
    !
    ! Passed Variables
    !-------------------
    type(physics_state),intent(in)   :: state
    real(r8),           intent(in)   :: ztodt
    type(physics_ptend),intent(out)  :: ptend
    type(cam_in_t),     intent(inout):: cam_in
    !
    ! Local Values
    !----------------
!    real(r8) :: T         (state%ncol,pver)   ! T temporary
!    real(r8) :: qv        (state%ncol,pver)   ! Q temporary (specific humidity)
!    real(r8) :: U         (state%ncol,pver)   ! U temporary
!    real(r8) :: V         (state%ncol,pver)   ! V temporary
    logical  :: lq        (pcnst)             ! Calc tendencies?
!    real(r8) :: dqdt_vdiff(state%ncol,pver)   ! PBL Q vertical diffusion tend kg/kg/s
!    real(r8) :: dtdt_vdiff(state%ncol,pver)   ! PBL T vertical diffusion tend  K/s
!    real(r8) :: Km        (state%ncol,pver+1) ! Eddy diffusivity at layer interfaces (m2/s)
!    real(r8) :: Ke        (state%ncol,pver+1) ! Eddy diffusivity at layer interfaces (m2/s)
!    real(r8) :: VSE       (state%ncol,pver)   ! Dry Static Energy divided by Cp (K)
!    real(r8) :: Zm        (state%ncol,pver)   ! 
!    real(r8) :: Zi        (state%ncol,pver)   ! 
!    real(r8) :: Z_pbl     (state%ncol)        ! 
!    real(r8) :: Rf        (state%ncol,pver)   ! 
!    real(r8) :: Tsfc      (state%ncol)        ! Surface T 
!    real(r8) :: Qsfc      (state%ncol)        ! Surface Q (saturated)
!    real(r8) :: Cdrag     (state%ncol)        ! Cdrag coef from surface calculation
    real(r8) :: dTs       (state%ncol)   
    real(r8) :: dUa       (state%ncol,pver)   
    real(r8) :: dVa       (state%ncol,pver)   
    real(r8) :: dTa       (state%ncol,pver)   
    real(r8) :: dQa       (state%ncol,pver)   
    integer  :: lchnk                        ! chunk identifier
    integer  :: ncol                         ! number of atmospheric columns
    integer  :: k                            ! loop index

    real(r8),allocatable :: Zm        (:,:)  ! 
    real(r8),allocatable :: Zi        (:,:)  ! 
    real(r8),allocatable :: Tsfc      (:)    ! Surface T 
    real(r8),allocatable :: Qsfc      (:)    ! Surface Q (saturated)
    real(r8),allocatable :: T         (:,:)  ! T temporary
    real(r8),allocatable :: qv        (:,:)  ! Q temporary (specific humidity)
    real(r8),allocatable :: U         (:,:)  ! U temporary
    real(r8),allocatable :: V         (:,:)  ! V temporary
    real(r8),allocatable :: dqdt_vdiff(:,:)  ! PBL Q vertical diffusion tend kg/kg/s
    real(r8),allocatable :: dtdt_vdiff(:,:)  ! PBL T vertical diffusion tend  K/s
    real(r8),allocatable :: Km        (:,:)  ! Eddy diffusivity at layer interfaces (m2/s)
    real(r8),allocatable :: Ke        (:,:)  ! Eddy diffusivity at layer interfaces (m2/s)
    real(r8),allocatable :: VSE       (:,:)  ! Dry Static Energy divided by Cp (K)
    real(r8),allocatable :: Cdrag     (:)    ! Cdrag coef from surface calculation
    real(r8),allocatable :: Z_pbl     (:)    ! 
    real(r8),allocatable :: Rf        (:,:)  ! 

    integer :: ii_diag
    real(r8):: rtd_diag

    allocate(Zm        (state%ncol,pver))
    allocate(Zi        (state%ncol,pver))
    allocate(Tsfc      (state%ncol)     )
    allocate(Qsfc      (state%ncol)     )
    allocate(T         (state%ncol,pver))
    allocate(qv        (state%ncol,pver))
    allocate(U         (state%ncol,pver))
    allocate(V         (state%ncol,pver))
    allocate(dqdt_vdiff(state%ncol,pver))
    allocate(dtdt_vdiff(state%ncol,pver))
    allocate(Km        (state%ncol,pver+1))
    allocate(Ke        (state%ncol,pver+1))
    allocate(VSE       (state%ncol,pver))
    allocate(Cdrag     (state%ncol)     )
    allocate(Z_pbl     (state%ncol)     )
    allocate(Rf        (state%ncol,pver))

    ! Set local copies of values
    !---------------------------------
    lchnk             = state%lchnk
    ncol              = state%ncol
    Zm    (:ncol,:)   = state%zm    (:ncol,:)
    Zi    (:ncol,1:pver)   = state%zi    (:ncol,1:pver)
    T     (:ncol,:)   = state%T     (:ncol,:)
    U     (:ncol,:)   = state%U     (:ncol,:)
    V     (:ncol,:)   = state%V     (:ncol,:)
    qv    (:ncol,:)   = state%Q     (:ncol,:,1)

    ! Initialize individual parameterization tendencies
    !-----------------------------------------------------
    lq           = .false.
    lq(1)        = .true.
    call physics_ptend_init(ptend,state%psetcols,'Frierson pbl_tend', &
                                       ls=.true., lu=.true., lv=.true., lq=lq)

    ! Call the Selected PBL routine  
    !--------------------------------------------------------
    if(PBL_OPT.eq.FRIERSON_IMPL_DIAG) then
      Tsfc(:ncol) = Tsurf(:ncol,CURR(lchnk),lchnk)
      Qsfc(:ncol) = Qsurf(:ncol,lchnk)
      call frierson_pbl_impl_diag(ncol, pver, ztodt,state%pmid (:ncol,:),  &
                                                    state%pint (:ncol,:),  &
                                                    state%rpdel(:ncol,:),  &
                                                             Zm(:ncol,:),  &
                                                             Zi(:ncol,:),  &
                                                    state%ps   (:ncol)  ,  &
                                                          Tsfc (:ncol)  ,  &
                                                          Qsfc (:ncol)  ,  &
                                                              T(:ncol,:),  &
                                                              U(:ncol,:),  &
                                                              V(:ncol,:),  &
                                                             qv(:ncol,:),  &
                                                     dqdt_vdiff(:ncol,:),  &
                                                     dtdt_vdiff(:ncol,:),  &
                                                             Km(:ncol,:),  &
                                                             Ke(:ncol,:),  &
                                                            VSE(:ncol,:),  &
                                                          Cdrag(:ncol)  ,  &
                                                          Z_pbl(:ncol),Rf(:ncol,:), &
                                                         Fsolar(:ncol,lchnk),       &
                                                          Fdown(:ncol,lchnk)        )
    elseif(PBL_OPT.eq.FRIERSON_IMPL_DIAGC) then
      Tsfc(:ncol) = Tsurf(:ncol,CURR(lchnk),lchnk)
      Qsfc(:ncol) = Qsurf(:ncol,lchnk)
      call frierson_pbl_impl_diagC(ncol, pver, ztodt,state%pmid (:ncol,:),  &
                                                    state%pint (:ncol,:),  &
                                                    state%rpdel(:ncol,:),  &
                                                             Zm(:ncol,:),  &
                                                             Zi(:ncol,:),  &
                                                    state%ps   (:ncol)  ,  &
                                                          Tsfc (:ncol)  ,  &
                                                          Qsfc (:ncol)  ,  &
                                                              T(:ncol,:),  &
                                                              U(:ncol,:),  &
                                                              V(:ncol,:),  &
                                                             qv(:ncol,:),  &
                                                     dqdt_vdiff(:ncol,:),  &
                                                     dtdt_vdiff(:ncol,:),  &
                                                             Km(:ncol,:),  &
                                                             Ke(:ncol,:),  &
                                                            VSE(:ncol,:),  &
                                                          Cdrag(:ncol)  ,  &
                                                          Z_pbl(:ncol),Rf(:ncol,:), &
                                                         Fsolar(:ncol,lchnk),       &
                                                          Fdown(:ncol,lchnk)        )
    elseif(PBL_OPT.eq.FRIERSON_IMPL_CONFIG1) then
      Tsfc(:ncol) = Tsurf(:ncol,CURR(lchnk),lchnk)
      Qsfc(:ncol) = Qsurf(:ncol,lchnk)
      call frierson_pbl_surface_impl(ncol, ztodt, Ft      (:ncol,lchnk), &
                                                  Fq      (:ncol,lchnk), &
                                                  Fu      (:ncol,lchnk), &
                                                  Fv      (:ncol,lchnk), &
                                                  Fup     (:ncol,lchnk), &
                                                  Fdown   (:ncol,lchnk), &
                                                  Fsolar  (:ncol,lchnk), &
                                                   dFt_dTs(:ncol,lchnk), &
                                                   dFq_dTs(:ncol,lchnk), &
                                                  dFup_dTs(:ncol,lchnk), &
                                                  state%ps(:ncol),       &
                                                      Tsfc(:ncol),       &
                                                      Qsfc(:ncol)        )
      dTs(:) = Tsfc(:) - Tsfc_bc(:,lchnk)
      call frierson_pbl_up_impl(ncol, pver, ztodt, dTs   (:ncol),          &
                                                   Fval_t(:ncol,:,lchnk),  &
                                                   Fval_q(:ncol,:,lchnk),  &
                                                   Fval_u(:ncol,:,lchnk),  &
                                                   Fval_v(:ncol,:,lchnk),  &
                                                   Eval_m(:ncol,:,lchnk),  &
                                                   Eval_e(:ncol,:,lchnk),  &
                                                   FN_t  (:ncol,lchnk),    &
                                                   FN_q  (:ncol,lchnk),    &
                                                   FN_u  (:ncol,lchnk),    &
                                                   FN_v  (:ncol,lchnk),    &
                                                   EN_t  (:ncol,lchnk),    &
                                                   EN_q  (:ncol,lchnk),    &
                                                   T     (:ncol,:),        &
                                                   U     (:ncol,:),        &
                                                   V     (:ncol,:),        &
                                                   qv    (:ncol,:),        &
                                                   dTa   (:ncol,:),        &
                                                   dQa   (:ncol,:),        &
                                                   dUa   (:ncol,:),        &
                                                   dVa   (:ncol,:)         )
    elseif(PBL_OPT.eq.FRIERSON_EXPSFC_DIAG) then
      Tsfc(:ncol)   = Tsurf(:ncol,CURR(lchnk),lchnk)
      Qsfc(:ncol)   = Qsurf(:ncol,lchnk)
      call frierson_pbl_expsfc_diag(ncol, pver, ztodt, state%pmid (:ncol,:),  &
                                                       state%pint (:ncol,:),  &
                                                       state%rpdel(:ncol,:),  &
                                                                Zm(:ncol,:),  &
                                                                Zi(:ncol,:),  &
                                                       state%ps   (:ncol)  ,  &
                                                             Tsfc (:ncol)  ,  &
                                                             Qsfc (:ncol)  ,  &
                                                                 T(:ncol,:),  &
                                                                 U(:ncol,:),  &
                                                                 V(:ncol,:),  &
                                                                qv(:ncol,:),  &
                                                        dqdt_vdiff(:ncol,:),  &
                                                        dtdt_vdiff(:ncol,:),  &
                                                                Km(:ncol,:),  &
                                                                Ke(:ncol,:),  &
                                                               VSE(:ncol,:),  &
                                                             Cdrag(:ncol)  ,  &
                                                             Z_pbl(:ncol),Rf(:ncol,:), &
                                                            Fsolar(:ncol,lchnk),       &
                                                             Fdown(:ncol,lchnk)        )
    elseif(PBL_OPT.eq.FRIERSON_EXPSFC_DIAGC) then
      Tsfc(:ncol)   = Tsurf(:ncol,CURR(lchnk),lchnk)
      Qsfc(:ncol)   = Qsurf(:ncol,lchnk)
      call frierson_pbl_expsfc_diagC(ncol, pver, ztodt, state%pmid (:ncol,:),  &
                                                        state%pint (:ncol,:),  &
                                                        state%rpdel(:ncol,:),  &
                                                                 Zm(:ncol,:),  &
                                                                 Zi(:ncol,:),  &
                                                        state%ps   (:ncol)  ,  &
                                                              Tsfc (:ncol)  ,  &
                                                              Qsfc (:ncol)  ,  &
                                                                  T(:ncol,:),  &
                                                                  U(:ncol,:),  &
                                                                  V(:ncol,:),  &
                                                                 qv(:ncol,:),  &
                                                         dqdt_vdiff(:ncol,:),  &
                                                         dtdt_vdiff(:ncol,:),  &
                                                                 Km(:ncol,:),  &
                                                                 Ke(:ncol,:),  &
                                                                VSE(:ncol,:),  &
                                                              Cdrag(:ncol)  ,  &
                                                              Z_pbl(:ncol),Rf(:ncol,:), &
                                                             Fsolar(:ncol,lchnk),       &
                                                              Fdown(:ncol,lchnk)        )
    elseif(PBL_OPT.eq.FRIERSON_EXPSFC_CONFIG1) then
      Tsfc(:ncol)   = Tsurf(:ncol,CURR(lchnk),lchnk)
      Qsfc(:ncol)   = Qsurf(:ncol,lchnk)
      call frierson_pbl_flux_expsfc(ncol, ztodt, Zm(:ncol,pver),  &
                                         state%pmid(:ncol,pver),  &
                                                  T(:ncol,pver),  &
                                                  U(:ncol,pver),  &
                                                  V(:ncol,pver),  &
                                                 qv(:ncol,pver),  &
                                           state%ps(:ncol),       &
                                               Tsfc(:ncol),       &
                                               Qsfc(:ncol),       &
                                                 Ft(:ncol,lchnk), &
                                                 Fq(:ncol,lchnk), &
                                                 Fu(:ncol,lchnk), &
                                                 Fv(:ncol,lchnk), &
                                                Fup(:ncol,lchnk)  )

      call frierson_pbl_surface_expsfc(ncol, ztodt, Ft      (:ncol,lchnk), &
                                                    Fq      (:ncol,lchnk), &
                                                    Fu      (:ncol,lchnk), &
                                                    Fv      (:ncol,lchnk), &
                                                    Fup     (:ncol,lchnk), &
                                                    Fdown   (:ncol,lchnk), &
                                                    Fsolar  (:ncol,lchnk), &
                                                    state%ps(:ncol),       &
                                                        Tsfc(:ncol),       &
                                                        Qsfc(:ncol)        )
  
      call frierson_pbl_up_expsfc(ncol, pver, ztodt, Fval_t(:ncol,:,lchnk),  &
                                                     Fval_q(:ncol,:,lchnk),  &
                                                     Fval_u(:ncol,:,lchnk),  &
                                                     Fval_v(:ncol,:,lchnk),  &
                                                     Eval_m(:ncol,:,lchnk),  &
                                                     Eval_e(:ncol,:,lchnk),  &
                                                     FN_t  (:ncol,lchnk),    &
                                                     FN_q  (:ncol,lchnk),    &
                                                     FN_u  (:ncol,lchnk),    &
                                                     FN_v  (:ncol,lchnk),    &
                                                     T     (:ncol,:),        &
                                                     U     (:ncol,:),        &
                                                     V     (:ncol,:),        &
                                                     qv    (:ncol,:),        &
                                                     dTa   (:ncol,:),        &
                                                     dQa   (:ncol,:),        &
                                                     dUa   (:ncol,:),        &
                                                     dVa   (:ncol,:)         )

    elseif(PBL_OPT.eq.FRIERSON_EXPSFC_CONFIG2) then
      Tsfc  (:ncol)     = Tsfc_bc(:ncol,lchnk)
      Qsfc  (:ncol)     = Qsfc_bc(:ncol,lchnk)
      call frierson_pbl_up_expsfc(ncol, pver, ztodt, Fval_t(:ncol,:,lchnk),  &
                                                     Fval_q(:ncol,:,lchnk),  &
                                                     Fval_u(:ncol,:,lchnk),  &
                                                     Fval_v(:ncol,:,lchnk),  &
                                                     Eval_m(:ncol,:,lchnk),  &
                                                     Eval_e(:ncol,:,lchnk),  &
                                                     FN_t  (:ncol,lchnk),    &
                                                     FN_q  (:ncol,lchnk),    &
                                                     FN_u  (:ncol,lchnk),    &
                                                     FN_v  (:ncol,lchnk),    &
                                                     T     (:ncol,:),        &
                                                     U     (:ncol,:),        &
                                                     V     (:ncol,:),        &
                                                     qv    (:ncol,:),        &
                                                     dTa   (:ncol,:),        &
                                                     dQa   (:ncol,:),        &
                                                     dUa   (:ncol,:),        &
                                                     dVa   (:ncol,:)         )

    elseif(PBL_OPT.eq.FRIERSON_EXPSFC        ) then
      Tsfc  (:ncol)     = cam_in%ts (:ncol)
      Qsfc  (:ncol)     = cam_in%ssq(:ncol)
!      Tsfc  (:ncol)     = Tsfc_bc(:ncol,lchnk)
!      Qsfc  (:ncol)     = Qsfc_bc(:ncol,lchnk)
      call frierson_pbl_up_expsfc(ncol, pver, ztodt, Fval_t(:ncol,:,lchnk),  &
                                                     Fval_q(:ncol,:,lchnk),  &
                                                     Fval_u(:ncol,:,lchnk),  &
                                                     Fval_v(:ncol,:,lchnk),  &
                                                     Eval_m(:ncol,:,lchnk),  &
                                                     Eval_e(:ncol,:,lchnk),  &
                                                     FN_t  (:ncol,lchnk),    &
                                                     FN_q  (:ncol,lchnk),    &
                                                     FN_u  (:ncol,lchnk),    &
                                                     FN_v  (:ncol,lchnk),    &
                                                     T     (:ncol,:),        &
                                                     U     (:ncol,:),        &
                                                     V     (:ncol,:),        &
                                                     qv    (:ncol,:),        &
                                                     dTa   (:ncol,:),        &
                                                     dQa   (:ncol,:),        &
                                                     dUa   (:ncol,:),        &
                                                     dVa   (:ncol,:)         )

    else
      ! ** ERROR STOP
    endif


!DIAGtest    
    Tsurf(:ncol,CURR(lchnk),lchnk) = Tsfc(:ncol)
    Qsurf(:ncol,lchnk)             = Qsfc(:ncol)
!DIAGtest   

    ! Back out tendencies from updated fields
    !-----------------------------------------
    do k = 1, pver
      ptend%s(:ncol,k  ) = (T (:,k)-state%T(:ncol,k  ))/ztodt*cpair
      ptend%u(:ncol,k  ) = (U (:,k)-state%U(:ncol,k  ))/ztodt
      ptend%v(:ncol,k  ) = (V (:,k)-state%V(:ncol,k  ))/ztodt
      ptend%q(:ncol,k,1) = (qv(:,k)-state%q(:ncol,k,1))/ztodt
    end do

    ! Archive diagnostic fields
    !----------------------------
!    call outfld('KVH' ,Ke        ,ncol ,lchnk) ! Eddy diffusivity (heat and moisture,m2/s)
!    call outfld('KVM' ,Km        ,ncol ,lchnk) ! Eddy diffusivity (momentum, m2/s)
!    call outfld('VSE' ,VSE       ,ncol ,lchnk) ! Virtual Dry Static Energy divided by Cp (K)
!    call outfld('Zm'  ,Zm        ,ncol ,lchnk) ! 
!    call outfld('Z_pbl',Z_pbl    ,ncol ,lchnk) ! 
!    call outfld('Rf'  ,Rf        ,ncol ,lchnk) ! 
!    call outfld('DUV' ,ptend%u   ,pcols,lchnk) ! PBL u tendency (m/s2)
!    call outfld('DVV' ,ptend%v   ,pcols,lchnk) ! PBL v tendency (m/s2)
!    call outfld('DTV' ,dtdt_vdiff,ncol ,lchnk) ! PBL + surface flux T tendency (K/s)
!    call outfld('VD01',dqdt_vdiff,ncol ,lchnk) ! PBL + surface flux Q tendency (kg/kg/s)
!    call outfld('Cdrag',Cdrag    ,ncol ,lchnk) ! 

    call outfld('R_Tsurf' , Tsurf (:ncol,CURR(lchnk),lchnk),ncol,lchnk)  !DIAG
    call outfld('R_Qsurf' , Qsurf (:ncol,lchnk),ncol,lchnk)              !DIAG

    deallocate(Zm        )
    deallocate(Zi        )
    deallocate(Tsfc      )
    deallocate(Qsfc      )
    deallocate(T         )
    deallocate(qv        )
    deallocate(U         )
    deallocate(V         )
    deallocate(dqdt_vdiff)
    deallocate(dtdt_vdiff)
    deallocate(Km        )
    deallocate(Ke        )
    deallocate(VSE       )
    deallocate(Cdrag     )
    deallocate(Z_pbl     )
    deallocate(Rf        )

    ! End Routine
    !--------------
    return
  end subroutine frierson_pbl_tend
  !============================================================================


  !============================================================================
  subroutine frierson_radiative_tend(state, ptend, ztodt,cam_in,cam_out)
    !
    ! frierson_radiative_tend: Run the radiatvie process
    !=========================================================================
    use physics_types,only: physics_state, physics_ptend
    use physics_types,only: physics_ptend_init
    use physconst,    only: cpair
    use phys_grid,    only: get_rlat_all_p
    use frierson,     only: frierson_radiation
    !
    ! Passed Variables
    !------------------
    type(physics_state),intent(in)   :: state
    real(r8)           ,intent(in)   :: ztodt
    type(physics_ptend),intent(out)  :: ptend
    type(cam_in_t),     intent(inout):: cam_in
    type(cam_out_t),    intent(inout):: cam_out
    !
    ! Local Values
    !---------------
    real(r8):: T           (state%ncol,pver) ! T temporary
    real(r8):: qv          (state%ncol,pver) ! Q temporary
    real(r8):: dtdt_heating(state%ncol,pver) ! temperature tendency from relaxation in K/s
    real(r8):: Tsfc        (state%ncol)      ! Surface T 
    real(r8):: Qsfc        (state%ncol)      ! Surface Q (saturated)
    logical :: lq(pcnst)                     ! Calc tendencies?
    integer :: lchnk                         ! chunk identifier
    integer :: ncol                          ! number of atmospheric columns
    integer :: k                             ! loop index

    ! Copy to local values 
    !-------------------------------------------------
    lchnk         = state%lchnk
    ncol          = state%ncol
    T   (:ncol,:) = state%T(:ncol,:)
    qv  (:ncol,:) = state%Q(:ncol,:,1)

    if((PBL_OPT.eq.FRIERSON_IMPL_DIAG     ).or.(PBL_OPT.eq.FRIERSON_IMPL_DIAGC    ).or. &
       (PBL_OPT.eq.FRIERSON_IMPL_CONFIG1  ).or.(PBL_OPT.eq.FRIERSON_EXPSFC_DIAG   ).or. &
       (PBL_OPT.eq.FRIERSON_EXPSFC_DIAGC  ).or.(PBL_OPT.eq.FRIERSON_EXPSFC_CONFIG1)     ) then
      Tsfc(:ncol)   = Tsurf(:ncol,CURR(lchnk),lchnk)
      Qsfc(:ncol)   = Qsurf(:ncol,lchnk)
    else !  ((PBL_OPT.eq.FRIERSON_EXPSFC_CONFIG2).or.(PBL_OPT.eq.FRIERSON_EXPSFC))
      Tsfc(:ncol)   = cam_in%ts (:ncol)
      Qsfc(:ncol)   = cam_in%ssq(:ncol)
    endif

    ! initialize individual parameterization tendencies
    !---------------------------------------------------
    lq(:) = .false.
    call physics_ptend_init(ptend, state%psetcols, 'Frierson radiative_tend',   &
                                        ls=.true., lu=.false., lv=.false., lq=lq)

    ! Call the Selected radiative routine 
    !--------------------------------------------------------
    call frierson_radiation(ncol,pver,ztodt,clat(:ncol,lchnk), &
                                      state%pint(:ncol,:),     &
                                      state%pmid(:ncol,:),     &
                                        state%ps(:ncol),       &
                                            Tsfc(:ncol),       &
                                            Qsfc(:ncol),       &
                                               T(:ncol,:),     &
                                              qv(:ncol,:),     &
                                    dtdt_heating(:ncol,:),     &
                                          Fsolar(:ncol,lchnk), &
                                             Fup(:ncol,lchnk), &
                                           Fdown(:ncol,lchnk)  )
    Fnet  (:ncol,lchnk) = Fup(:ncol,lchnk) - Fdown (:ncol,lchnk)

    ! Copy downward LW radiative heating values to cam_out%
    !---------------------------------------------------------
    cam_out%flwds(:ncol) = Fdown (:ncol,lchnk)
    cam_out%netsw(:ncol) = Fsolar(:ncol,lchnk)
    cam_out%sols (:ncol) = Fsolar(:ncol,lchnk)
    cam_out%solsd(:ncol) = Fsolar(:ncol,lchnk)
    cam_out%soll (:ncol) = Fsolar(:ncol,lchnk)
    cam_out%solld(:ncol) = Fsolar(:ncol,lchnk)

    ! Back out tendencies from updated T field
    !--------------------------------------------
    do k = 1, pver
      ptend%s(:ncol,k) = (T(:,k)-state%T(:ncol,k))/ztodt*cpair
    end do

    ! Archive T tendency from temperature relaxation (mimics radiation, K/s)
    !-----------------------------------------------------------------------
    call outfld('QRS',dtdt_heating, ncol,lchnk) 

    ! End Routine
    !------------
    return
  end subroutine frierson_radiative_tend
  !============================================================================


  !=======================================================================
  subroutine frierson_surface_init(ncol, clat, PS, Tsfc, Qsfc)
    !
    !
    !==========================================================================
    ! Tuning parameters
    !---------------------
    real(r8),parameter:: T_min   = 271._r8             ! Minimum sst (K)
    real(r8),parameter:: del_T   = 39._r8  ! 29._r8    ! eq-polar difference sst (K)
    real(r8),parameter:: T_width = 26.0_r8*pi/180.0_r8 ! width parameter for sst (C)
    !
    ! Passed variables
    !--------------------
    integer ,intent(in) :: ncol
    real(r8),intent(in) :: clat (ncol)
    real(r8),intent(in) :: PS   (ncol)
    real(r8),intent(out):: Tsfc(ncol)
    real(r8),intent(out):: Qsfc(ncol)
    !
    ! Local values
    !--------------
    integer :: i

    ! set SST profile
    !------------------
    do i = 1, ncol
      Tsfc(i) = T_min + del_T*exp(-((clat(i)/T_width)**2)/2.0_r8)
      Qsfc(i) = epsilo*frierson_E0/PS(i)                             &
                *exp(-latvap/rh2o*((1._r8/Tsfc(i))-1._r8/frierson_T0))
    end do

    ! End Routine
    !--------------
    return
  end subroutine frierson_surface_init
  !=======================================================================

end module frierson_cam

