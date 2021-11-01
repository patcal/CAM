module Vmaps_mod
!===========================================================================
!
! MODULE: Vmap_mod
!
! DESCRIPTION: Module contaiing the known mapping of data in 
!              the vertical dimension
!
!===========================================================================
  ! Useful Modules
  !---------------------
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use spmd_utils,    only: masterproc
  use cam_abortutils,only: endrun
  use cam_logfile,   only: iulog
  use Dataset_mod,   only: dataset_var_t
  use PathDesc_type ,only: PathDesc_t
  use CAPT_mod 

  ! Set all Global values and routines to private by default 
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private
  save

  public :: vgrid_t
  public :: vmap_data_t 
  public :: vmap_desc_t 
  public :: vmap_t
  public :: known_vgrids_t
  public :: known_vmaps_t

  public :: init_Vmaps_mod
  public :: inq_Vmaps
  public :: print_known_Vmaps
  public :: copy_Vgrid
  public :: deallocate_Vgrid

  public ::   chk_Vmap
  public ::  init_Vmap
  public :: apply_Vmap

  public :: vmap_index
  private:: vmap_index_chk
  private:: vmap_index_init
  private:: vmap_index_apply

  public :: vmap_pressure
  private:: vmap_pressure_chk
  private:: vmap_pressure_init
  private:: vmap_pressure_apply

  public :: vmap_USER
  private:: vmap_USER_chk
  private:: vmap_USER_init
  private:: vmap_USER_apply

  ! Interfaces
  !-----------------
  interface vmap_index
    procedure:: vmap_index_chk
    procedure:: vmap_index_init
    procedure:: vmap_index_apply
  end interface

  interface vmap_pressure
    procedure:: vmap_pressure_chk
    procedure:: vmap_pressure_init
    procedure:: vmap_pressure_apply
  end interface

  interface vmap_USER
    procedure:: vmap_USER_chk
    procedure:: vmap_USER_init
    procedure:: vmap_USER_apply
  end interface

  ! Type Definitions
  !-------------------
  type vgrid_t
    character(len=80)            :: name
    character(len=32)            :: Type
    character(len=32)            :: Vdimname
    integer                      :: nlev  = 0
    integer                      :: nreq_var
    character(len=32),allocatable:: req_var (:)
    character(len=32),allocatable:: req_type(:)
    integer                      :: ngrid_var
    character(len=32),allocatable:: grid_var(:)
    real(r8),allocatable         :: HYA (:)
    real(r8),allocatable         :: HYB (:)
    real(r8),allocatable         :: PRES(:)
    real(r8),allocatable         :: HGHT(:)
    real(r8)                     :: Pref     = 100000._r8
  end type vgrid_t

  type preq_data_t
    character(len=32)   :: Name
    character(len=32)   :: Shape
    logical             :: CONST
    integer             :: ncol
    integer             :: nlev
    real(r8),allocatable:: val(:,:)
  end type preq_data_t

  type vmap_data_t
    integer                      :: npreq_in
    integer                      :: npreq_out
    type(preq_data_t),allocatable:: Vpreq_in (:)
    type(preq_data_t),allocatable:: Vpreq_out(:)
  end type vmap_data_t

  type vmap_t
    integer              :: idxVgrid_d
    integer              :: idxVgrid_m
    type(vgrid_t),pointer:: Vgrid_in    => null()
    type(vgrid_t),pointer:: Vgrid_out   => null()
    character(len=32)    :: VmapOpt
    character(len=32)    :: VmapFunction
    type(vmap_data_t)    :: VmapData
  end type vmap_t


  type vmap_desc_t
    character(len=32)            :: name
    character(len=32)            :: FunctionName
    character(len=32)            :: InType
    character(len=32)            :: OutType
    integer                      :: numMapOpts
    character(len=32),allocatable:: Opts      (:)
    integer                      :: nreq_in
    integer                      :: nreq_out
    character(len=32),allocatable:: Vreq_in (:)
    character(len=32),allocatable:: Vreq_out(:)
  end type vmap_desc_t

  type known_vgrids_t
    integer                  :: numVgrids
    type(vgrid_t),allocatable:: Vgrid(:)
  end type known_vgrids_t

  type known_vmaps_t
    integer                      :: numVmaps
    type(vmap_desc_t),allocatable:: Vmap(:)
  end type known_vmaps_t

  ! Global Values
  !------------------
  type(known_vgrids_t):: My_Vgrids
  type(known_vmaps_t ):: My_Vmaps

contains

  !================================================================
  subroutine init_Vmaps_mod
    ! 
    ! init_Vmaps_mod:  Initialize Vmaps for use.
    !
    !==============================================================
    !
    ! Passed variables
    !-------------------

    !
    ! Local Values
    !-------------
    integer:: Idx_HYP_LEV, Idx_HYP_ILEV, Idx_PRESSURE, Idx_HEIGHT, Idx_INDEX, Idx_USER

    ! Define Known Vertical Grids
    !--------------------------------
    My_Vgrids%numVgrids = 6
    allocate(My_Vgrids%Vgrid(My_Vgrids%numVgrids))
    Idx_HYP_LEV = 1
    My_Vgrids%Vgrid(1)%name        = "HYP_LEV"
    My_Vgrids%Vgrid(1)%Type        = "HYBRIDP"
    My_Vgrids%Vgrid(1)%Vdimname    = "lev"
    My_Vgrids%Vgrid(1)%ngrid_var   = 2
    My_Vgrids%Vgrid(1)%nreq_var    = 1
    allocate(My_Vgrids%Vgrid(1)%grid_var(My_Vgrids%Vgrid(1)%ngrid_var))
    allocate(My_Vgrids%Vgrid(1)%req_var (My_Vgrids%Vgrid(1)%nreq_var))
    allocate(My_Vgrids%Vgrid(1)%req_type(My_Vgrids%Vgrid(1)%nreq_var))
    My_Vgrids%Vgrid(1)%grid_var(1) = "hyam"
    My_Vgrids%Vgrid(1)%grid_var(2) = "hybm"
    My_Vgrids%Vgrid(1)%req_var (1) = "PS"
    My_Vgrids%Vgrid(1)%req_type(1) = "HORIZ_ONLY"
    Idx_HYP_ILEV = 2
    My_Vgrids%Vgrid(2)%name        = "HYP_ILEV"
    My_Vgrids%Vgrid(2)%Type        = "HYBRIDP"
    My_Vgrids%Vgrid(2)%Vdimname    = "ilev"
    My_Vgrids%Vgrid(2)%ngrid_var   = 2
    My_Vgrids%Vgrid(2)%nreq_var    = 1
    allocate(My_Vgrids%Vgrid(2)%grid_var(My_Vgrids%Vgrid(2)%ngrid_var))
    allocate(My_Vgrids%Vgrid(2)%req_var (My_Vgrids%Vgrid(2)%nreq_var))
    allocate(My_Vgrids%Vgrid(2)%req_type(My_Vgrids%Vgrid(2)%nreq_var))
    My_Vgrids%Vgrid(2)%grid_var(1) = "hyai"
    My_Vgrids%Vgrid(2)%grid_var(2) = "hybi"
    My_Vgrids%Vgrid(2)%req_var (1) = "PS"
    My_Vgrids%Vgrid(2)%req_type(1) = "HORIZ_ONLY"
    Idx_PRESSURE = 3
    My_Vgrids%Vgrid(3)%name        = "PRESSURE"
    My_Vgrids%Vgrid(3)%Type        = "PRESSURE"
    My_Vgrids%Vgrid(3)%Vdimname    = "pressure"
    My_Vgrids%Vgrid(3)%ngrid_var   = 1
    My_Vgrids%Vgrid(3)%nreq_var    = 0
    allocate(My_Vgrids%Vgrid(3)%grid_var(My_Vgrids%Vgrid(3)%ngrid_var))
    My_Vgrids%Vgrid(3)%grid_var(1) = "PRESSURE"
    Idx_HEIGHT = 4
    My_Vgrids%Vgrid(4)%name        = "HEIGHT"
    My_Vgrids%Vgrid(4)%Type        = "HEIGHT"
    My_Vgrids%Vgrid(4)%Vdimname    = "height"
    My_Vgrids%Vgrid(4)%ngrid_var   = 1
    My_Vgrids%Vgrid(4)%nreq_var    = 0
    allocate(My_Vgrids%Vgrid(4)%grid_var(My_Vgrids%Vgrid(4)%ngrid_var))
    My_Vgrids%Vgrid(4)%grid_var(1) = "HEIGHT"
    Idx_INDEX = 5
    My_Vgrids%Vgrid(5)%name        = "ANY"
    My_Vgrids%Vgrid(5)%Type        = "INDEX"
    My_Vgrids%Vgrid(5)%Vdimname    = "any"
    My_Vgrids%Vgrid(5)%ngrid_var   = 1
    My_Vgrids%Vgrid(5)%nreq_var    = 0
    allocate(My_Vgrids%Vgrid(5)%grid_var(My_Vgrids%Vgrid(5)%ngrid_var))
    My_Vgrids%Vgrid(5)%grid_var(1) = "ANYINDEX"
    Idx_USER = 6
    My_Vgrids%Vgrid(6)%name        = "ANY"
    My_Vgrids%Vgrid(6)%Type        = "USER"
    My_Vgrids%Vgrid(6)%Vdimname    = "any"
    My_Vgrids%Vgrid(6)%ngrid_var   = 1
    My_Vgrids%Vgrid(6)%nreq_var    = 0
    allocate(My_Vgrids%Vgrid(6)%grid_var(My_Vgrids%Vgrid(6)%ngrid_var))
    My_Vgrids%Vgrid(6)%grid_var(1) = "ANYINDEX"

    ! Define Known Vertical Maps
    !-----------------------------
    My_Vmaps%numVmaps = 5
    allocate(My_Vmaps%Vmap(My_Vmaps%numVmaps))

    ! Vmap for direct index copy
    !------------------------------
    My_Vmaps%Vmap(1)%name         = 'vmap_index'
    My_Vmaps%Vmap(1)%FunctionName = 'vmap_index'
    My_Vmaps%Vmap(1)%InType       = 'INDEX'
    My_Vmaps%Vmap(1)%OutType      = 'INDEX'
    My_Vmaps%Vmap(1)%numMapOpts   = 2
    allocate(My_Vmaps%Vmap(1)%Opts(My_Vmaps%Vmap(1)%numMapOpts ))
    My_Vmaps%Vmap(1)%Opts(1)      = 'INDEX'
    My_Vmaps%Vmap(1)%Opts(2)      = 'INDEX_SHIFT'
    My_Vmaps%Vmap(1)%nreq_in      = 0
    My_Vmaps%Vmap(1)%nreq_out     = 0

    ! Log-Pressure Interpolation Map
    !---------------------------------
    My_Vmaps%Vmap(2)%name         = 'vmap_HYP2HYP'
    My_Vmaps%Vmap(2)%FunctionName = 'vmap_Pressure'
    My_Vmaps%Vmap(2)%InType       = "HYBRIDP"
    My_Vmaps%Vmap(2)%OutType      = "HYBRIDP"
    My_Vmaps%Vmap(2)%numMapOpts   = 2
    allocate(My_Vmaps%Vmap(2)%Opts(My_Vmaps%Vmap(2)%numMapOpts))
    My_Vmaps%Vmap(2)%Opts(1)      = "LOG_PRESSURE"
    My_Vmaps%Vmap(2)%Opts(2)      = "PRESSURE"
    My_Vmaps%Vmap(2)%nreq_in      = 1
    My_Vmaps%Vmap(2)%nreq_out     = 1
    allocate(My_Vmaps%Vmap(2)%Vreq_in (My_Vmaps%Vmap(2)%nreq_in ))
    allocate(My_Vmaps%Vmap(2)%Vreq_out(My_Vmaps%Vmap(2)%nreq_out))
    My_Vmaps%Vmap(2)%Vreq_in (1)  = 'PS'
    My_Vmaps%Vmap(2)%Vreq_out(1)  = 'PS'

    My_Vmaps%Vmap(3)%name         = 'vmap_PRES2HYP'
    My_Vmaps%Vmap(3)%FunctionName = 'vmap_Pressure'
    My_Vmaps%Vmap(3)%InType       = "PRESSURE"
    My_Vmaps%Vmap(3)%OutType      = "HYBRIDP"
    My_Vmaps%Vmap(3)%numMapOpts   = 2
    allocate(My_Vmaps%Vmap(3)%Opts(My_Vmaps%Vmap(3)%numMapOpts))
    My_Vmaps%Vmap(3)%Opts(1)      = "LOG_PRESSURE"
    My_Vmaps%Vmap(3)%Opts(2)      = "PRESSURE"
    My_Vmaps%Vmap(3)%nreq_in      = 0
    My_Vmaps%Vmap(3)%nreq_out     = 1
    allocate(My_Vmaps%Vmap(3)%Vreq_out(My_Vmaps%Vmap(3)%nreq_out))
    My_Vmaps%Vmap(3)%Vreq_out(1)  = 'PS'

    My_Vmaps%Vmap(4)%name         = 'vmap_HYP2PRES'
    My_Vmaps%Vmap(4)%FunctionName = 'vmap_Pressure'
    My_Vmaps%Vmap(4)%InType       = "HYBRIDP"
    My_Vmaps%Vmap(4)%OutType      = "PRESSURE"
    My_Vmaps%Vmap(4)%numMapOpts   = 2
    allocate(My_Vmaps%Vmap(4)%Opts(My_Vmaps%Vmap(4)%numMapOpts))
    My_Vmaps%Vmap(4)%Opts(1)      = "LOG_PRESSURE"
    My_Vmaps%Vmap(4)%Opts(2)      = "PRESSURE"
    My_Vmaps%Vmap(4)%nreq_in      = 1
    My_Vmaps%Vmap(4)%nreq_out     = 0
    allocate(My_Vmaps%Vmap(4)%Vreq_in (My_Vmaps%Vmap(4)%nreq_in ))
    My_Vmaps%Vmap(4)%Vreq_in (1)  = 'PS'

    ! USER Defined Map
    !---------------------------------
    My_Vmaps%Vmap(5)%name         = 'vmap_USER'
    My_Vmaps%Vmap(5)%FunctionName = 'vmap_USER'
    My_Vmaps%Vmap(5)%InType       = "USER"
    My_Vmaps%Vmap(5)%OutType      = "USER"
    My_Vmaps%Vmap(5)%numMapOpts   = 1
    allocate(My_Vmaps%Vmap(5)%Opts(My_Vmaps%Vmap(5)%numMapOpts))
    My_Vmaps%Vmap(5)%Opts   (1)   = "USER"
    My_Vmaps%Vmap(5)%nreq_in      = 0
    My_Vmaps%Vmap(5)%nreq_out     = 0

    ! End Routine
    !-----------
    return
  end subroutine init_Vmaps_mod
  !================================================================


  !================================================================
  subroutine inq_Vmaps(Vgrids,Vmaps)
    ! 
    ! inq_Vmaps:  Return a datastructure contaiing the number of 
    !             maps in the vertical and the metadata for each 
    !             map. The metadata includes the grid templates 
    !             for input/output data and the mapping method.
    !
    !==============================================================
    !
    ! Passed variables
    !-------------------
    type(known_vgrids_t):: Vgrids
    type(known_vmaps_t ):: Vmaps
    !
    ! Local values
    !--------------
    integer:: nn,kk

    ! Return the set of Known Vertical Grids
    !----------------------------------------
    if(allocated(Vgrids%Vgrid)) then 
      do nn=1,Vgrids%numVgrids
        if(allocated(Vgrids%Vgrid(nn)%grid_var)) deallocate(Vgrids%Vgrid(nn)%grid_var)
        if(allocated(Vgrids%Vgrid(nn)%req_var )) deallocate(Vgrids%Vgrid(nn)%req_var )
      end do
      deallocate(Vgrids%Vgrid)
    endif

    Vgrids%numVgrids = My_Vgrids%numVgrids
    allocate(Vgrids%Vgrid(Vgrids%numVgrids))
    do nn=1,Vgrids%numVgrids
      Vgrids%Vgrid(nn)%name      = trim(My_Vgrids%Vgrid(nn)%name)
      Vgrids%Vgrid(nn)%Type      = trim(My_Vgrids%Vgrid(nn)%Type)
      Vgrids%Vgrid(nn)%Vdimname  = trim(My_Vgrids%Vgrid(nn)%Vdimname)
      Vgrids%Vgrid(nn)%ngrid_var = My_Vgrids%Vgrid(nn)%ngrid_var
      Vgrids%Vgrid(nn)%nreq_var  = My_Vgrids%Vgrid(nn)%nreq_var
      if(Vgrids%Vgrid(nn)%ngrid_var.ge.1) then
        allocate(Vgrids%Vgrid(nn)%grid_var(Vgrids%Vgrid(nn)%ngrid_var))
        do kk=1,Vgrids%Vgrid(nn)%ngrid_var
          Vgrids%Vgrid(nn)%grid_var(kk) = trim(My_Vgrids%Vgrid(nn)%grid_var(kk))
        end do
      endif
      if(Vgrids%Vgrid(nn)%nreq_var.ge.1) then
        allocate(Vgrids%Vgrid(nn)%req_var(Vgrids%Vgrid(nn)%nreq_var))
        do kk=1,Vgrids%Vgrid(nn)%nreq_var
          Vgrids%Vgrid(nn)%req_var(kk) = trim(My_Vgrids%Vgrid(nn)%req_var(kk))
        end do
      endif
    end do

    ! Return the set of Known Vertical Maps
    !----------------------------------------
    if(allocated(Vmaps%Vmap)) then
      do nn=1,Vmaps%numVmaps
        if(allocated(Vmaps%Vmap(nn)%Opts   )) deallocate(Vmaps%Vmap(nn)%Opts   )
      end do
      deallocate(Vmaps%Vmap)
    endif

    Vmaps%numVmaps = My_Vmaps%numVmaps
    allocate(Vmaps%Vmap(Vmaps%numVmaps))
    do nn=1,Vmaps%numVmaps
      Vmaps%Vmap(nn)%name         = trim(My_Vmaps%Vmap(nn)%name   )
      Vmaps%Vmap(nn)%FunctionName = trim(My_Vmaps%Vmap(nn)%FunctionName)
      Vmaps%Vmap(nn)%InType       = trim(My_Vmaps%Vmap(nn)%InType )
      Vmaps%Vmap(nn)%OutType      = trim(My_Vmaps%Vmap(nn)%OutType)
      Vmaps%Vmap(nn)%numMapOpts   =      My_Vmaps%Vmap(nn)%numMapOpts
      Vmaps%Vmap(nn)%nreq_in      =      My_Vmaps%Vmap(nn)%nreq_in
      Vmaps%Vmap(nn)%nreq_out     =      My_Vmaps%Vmap(nn)%nreq_out
      allocate(Vmaps%Vmap(nn)%Opts(Vmaps%Vmap(nn)%numMapOpts))
      do kk=1,Vmaps%Vmap(nn)%numMapOpts
        Vmaps%Vmap(nn)%Opts(kk) = My_Vmaps%Vmap(nn)%Opts(kk)
      end do
      allocate(Vmaps%Vmap(nn)%Vreq_in(My_Vmaps%Vmap(nn)%nreq_in))
      do kk=1,My_Vmaps%Vmap(nn)%nreq_in
        Vmaps%Vmap(nn)%Vreq_in(kk) = trim(My_Vmaps%Vmap(nn)%Vreq_in(kk))
      end do
      allocate(Vmaps%Vmap(nn)%Vreq_out(My_Vmaps%Vmap(nn)%nreq_out))
      do kk=1,My_Vmaps%Vmap(nn)%nreq_out
        Vmaps%Vmap(nn)%Vreq_out(kk) = trim(My_Vmaps%Vmap(nn)%Vreq_out(kk))
      end do
    end do

    ! End Routine
    !-----------
    return
  end subroutine inq_Vmaps
  !================================================================


  !================================================================
  subroutine print_known_Vmaps
    ! 
    ! print_known_Vmaps:  Print the metadata for the known
    !                     Vgrids and Vmaps between them.
    !
    !==============================================================
    !
    ! Local values
    !--------------
    integer:: nn,ii

    if(masterproc) then
      ! print out the set of known grids
      !----------------------------------------
      write(iulog,*) ' '
      write(iulog,*) '============================================================== '
      write(iulog,*) 'METADATA for Vmap Module'
      write(iulog,*) '    Number of Known Vgrids=',My_Vgrids%numVgrids
      write(iulog,*) '    Number of Known  Vmaps=',My_Vmaps%numVmaps
      write(iulog,*) ' '
      write(iulog,*) '  KNOWN VERTICAL GRIDS:'
      do nn=1,My_Vgrids%numVgrids
        write(iulog,*) '    nn=',nn,' -------------------------'
        write(iulog,*) '           name=',trim(My_Vgrids%Vgrid(nn)%name)
        write(iulog,*) '           type=',trim(My_Vgrids%Vgrid(nn)%Type)
        write(iulog,*) '        dimname=',trim(My_Vgrids%Vgrid(nn)%Vdimname)
        write(iulog,*) '           nlev=',My_Vgrids%Vgrid(nn)%nlev
        write(iulog,*) '      ngrid_var=',My_Vgrids%Vgrid(nn)%ngrid_var
        do ii=1,My_Vgrids%Vgrid(nn)%ngrid_var
          write(iulog,*) '        Grid Variable',ii,                                  &
                                        ' name=',trim(My_Vgrids%Vgrid(nn)%grid_var(ii))
        end do
        write(iulog,*) '       nreq_var=',My_Vgrids%Vgrid(nn)%nreq_var
        do ii=1,My_Vgrids%Vgrid(nn)%nreq_var
          write(iulog,*) '        Required Variable',ii,                             &
                                        ' name=',trim(My_Vgrids%Vgrid(nn)%req_var(ii))
        end do
        write(iulog,*) '            Pref = ',My_Vgrids%Vgrid(nn)%Pref
        if(allocated(My_Vgrids%Vgrid(nn)%HYA)) then
          write(iulog,*) '           HYA = ALLOCATED'
        else
          write(iulog,*) '           HYA = NOT ALLOCATED'
        endif
        if(allocated(My_Vgrids%Vgrid(nn)%HYB)) then
          write(iulog,*) '           HYB = ALLOCATED'
        else
          write(iulog,*) '           HYB = NOT ALLOCATED'
        endif
        if(allocated(My_Vgrids%Vgrid(nn)%PRES)) then
          write(iulog,*) '          PRES = ALLOCATED'
        else
          write(iulog,*) '          PRES = NOT ALLOCATED'
        endif
        if(allocated(My_Vgrids%Vgrid(nn)%HGHT)) then
          write(iulog,*) '          HGHT = ALLOCATED'
        else
          write(iulog,*) '          HGHT = NOT ALLOCATED'
        endif
      end do ! nn=1,My_Vgrids%numVgrids
      write(iulog,*) ' '
      write(iulog,*) '  KNOWN VERTICAL MAPS:'
      do nn=1,My_Vmaps%numVmaps
        write(iulog,*) '    nn=',nn,' -------------------------'
        write(iulog,*) '      function name     = ',trim(My_Vmaps%Vmap(nn)%name   )
        write(iulog,*) '      Input  Grid Type  = ',trim(My_Vmaps%Vmap(nn)%InType )
        write(iulog,*) '      Output Grid Type  = ',trim(My_Vmaps%Vmap(nn)%OutType)
        write(iulog,*) '      number of Map Opts= ',My_Vmaps%Vmap(nn)%numMapOpts
        do ii=1,My_Vmaps%Vmap(nn)%numMapOpts
          write(iulog,*) '          MapOpt=',ii,' Name= ',trim(My_Vmaps%Vmap(nn)%Opts(ii))
        end do
        write(iulog,*) '         nreq_in= ',My_Vmaps%Vmap(nn)%nreq_in
        do ii=1,My_Vmaps%Vmap(nn)%nreq_in
          write(iulog,*) '        Required Input Variable= ',ii,             &
                                        trim(My_Vmaps%Vmap(nn)%Vreq_in(ii))
        end do
        write(iulog,*) '         nreq_out= ',My_Vmaps%Vmap(nn)%nreq_out
        do ii=1,My_Vmaps%Vmap(nn)%nreq_out
          write(iulog,*) '        Required Output Variable= ',ii,             &
                                        trim(My_Vmaps%Vmap(nn)%Vreq_out(ii))
        end do
      end do ! nn=1,My_Vmaps%numVmaps
      write(iulog,*) '============================================================== '
      write(iulog,*) ' '
    endif

    ! End Routine
    !-----------
    return
  end subroutine print_known_Vmaps
  !================================================================


  !================================================================
  subroutine copy_Vgrid(OutVgrid,InVgrid)
    ! 
    ! copy_Vgrid: Copy the contents of the Input model grid
    !               to the Output model grid.
    !==============================================================
    !
    ! Passed variables
    !-------------------
    type(vgrid_t),intent(in )::  InVgrid
    type(vgrid_t),intent(out):: OutVgrid
    !
    ! Local Values
    !---------------
    integer:: ii

    OutVgrid%name      = trim(InVgrid%name)
    OutVgrid%Type      = trim(InVgrid%Type)
    OutVgrid%Vdimname  = trim(InVgrid%Vdimname)
    OutVgrid%nlev      =      InVgrid%nlev
    OutVgrid%Pref      =      InVgrid%Pref
    OutVgrid%ngrid_var =      InVgrid%ngrid_var
    OutVgrid%nreq_var  =      InVgrid%nreq_var
    if(allocated(OutVgrid%grid_var)) deallocate(OutVgrid%grid_var)
    if(allocated(OutVgrid%req_var) ) deallocate(OutVgrid%req_var )
    if(allocated(OutVgrid%HYA)     ) deallocate(OutVgrid%HYA     )
    if(allocated(OutVgrid%HYB)     ) deallocate(OutVgrid%HYB     )
    if(allocated(OutVgrid%PRES)    ) deallocate(OutVgrid%PRES    )
    if(allocated(OutVgrid%HGHT)    ) deallocate(OutVgrid%HGHT    )
    if(OutVgrid%ngrid_var.gt.0) then
      allocate(OutVgrid%grid_var(OutVgrid%ngrid_var))
      do ii=1,OutVgrid%ngrid_var
        OutVgrid%grid_var(ii) = trim(InVgrid%grid_var(ii))
      end do
    endif
    if(OutVgrid%nreq_var.gt.0) then
      allocate(OutVgrid%req_var(OutVgrid%nreq_var))
      do ii=1,OutVgrid%nreq_var
        OutVgrid%req_var(ii) = trim(InVgrid%req_var(ii))
      end do
    endif
    if(allocated(InVgrid%HYA)) then
      allocate(OutVgrid%HYA (OutVgrid%nlev))
      OutVgrid%HYA (:) = InVgrid%HYA (:)
    endif
    if(allocated(InVgrid%HYB)) then
      allocate(OutVgrid%HYB (OutVgrid%nlev))
      OutVgrid%HYB (:) = InVgrid%HYB (:)
    endif
    if(allocated(InVgrid%PRES)) then
      allocate(OutVgrid%PRES(OutVgrid%nlev))
      OutVgrid%PRES(:) = InVgrid%PRES(:)
    endif
    if(allocated(InVgrid%HGHT)) then
      allocate(OutVgrid%HGHT(OutVgrid%nlev))
      OutVgrid%HGHT(:) = InVgrid%HGHT(:)
    endif

    ! End Routine
    !-----------
    return
  end subroutine copy_Vgrid
  !================================================================


  !================================================================
  subroutine deallocate_Vgrid(Vgrid)
    ! 
    ! deallocate_Vgrid: Deallocate Hgrid arrays
    !==============================================================
    !
    ! Passed variables
    !-------------------
    type(vgrid_t),intent(inout)::  Vgrid

    if(allocated(Vgrid%grid_var)) deallocate(Vgrid%grid_var)
    if(allocated(Vgrid%req_var) ) deallocate(Vgrid%req_var )
    if(allocated(Vgrid%HYA)     ) deallocate(Vgrid%HYA     )
    if(allocated(Vgrid%HYB)     ) deallocate(Vgrid%HYB     )
    if(allocated(Vgrid%PRES)    ) deallocate(Vgrid%PRES    )
    if(allocated(Vgrid%HGHT)    ) deallocate(Vgrid%HGHT    )

    ! End Routine
    !-----------
    return
  end subroutine deallocate_Vgrid
  !================================================================


  !================================================================
  logical function chk_Vmap(Vgrid_Out,Vgrid_In,I_VmapOpt,known_Vmap)
    !
    ! chk_Vmap: This is a logical function to test if the input/output
    !           vertical grids satisfiy the requirements of the 
    !           given Map name for the given known Vmap template.
    !
    !   This is an interface function which calls the Function
    !   specified in the known_Vmap structure to check the output 
    !   grid and the input grid for the given map option name.
    !==============================================================
    !
    ! Passed Variables
    !------------------
    type(vgrid_t)    ,intent(in):: Vgrid_Out
    type(vgrid_t)    ,intent(in):: Vgrid_In
    character(len=*) ,intent(in):: I_VmapOpt
    type(vmap_desc_t),intent(in):: known_Vmap
    ! 
    ! Local Values
    !--------------
    character(len=32):: VmapOpt

    VmapOpt = '                                '
    VmapOpt = trim(I_VmapOpt)

    ! This should be cleaned up ASAP, to use function pointers
    !-------------------------------------------------------------
    if(trim(known_Vmap%FunctionName).eq.'vmap_index') then
      ! Maps that use vertical indices for mapping
      !----------------------------------------------------
      chk_Vmap = vmap_index(Vgrid_Out,Vgrid_In,trim(VmapOpt),known_Vmap)
    elseif(trim(known_Vmap%FunctionName).eq.'vmap_pressure') then
      ! Maps that use P or ln(P)  for mapping
      !----------------------------------------------------
      chk_Vmap = vmap_pressure(Vgrid_Out,Vgrid_In,trim(VmapOpt),known_Vmap)
    elseif(trim(known_Vmap%FunctionName).eq.'vmap_USER') then
      ! USER defined Maps
      !----------------------------------------------------
      chk_Vmap = vmap_USER(Vgrid_Out,Vgrid_In,trim(VmapOpt),known_Vmap)
    else
      ! ERROR: Unknown FunctionName
      !----------------------------------------------------
      chk_Vmap = .false.
    endif

    ! End Routine
    !------------------
    return
  end function chk_Vmap
  !================================================================


  !================================================================
  logical function init_Vmap(Vmap)
    ! 
    ! init_Vmap:  This is a logical function to fill-in/initialize the 
    !             VmapData data structure with info for reading and 
    !             mapping data from the input Vgrid to the output
    !             Vgrid for the given mapping option. 
    !             (e.g. to hold static PHIS values that may be used by a map)
    !
    !   This is an interface function which calls the specified Function
    !   with the input/output grids and mapping option. 
    !
    !   DEVNOTE: For now the mapping data is passed thru a specific 
    !            data structure (vmap_data_t). For more general options 
    !            this will be something that needs to be sorted out. 
    !            ?? Can I pass VmapData as a generic class [ class(*) ] so 
    !            that the called VmapFunctions can initialze a datastructure
    !            tailored to the function needs? 
    !==============================================================
    !
    ! Passed Variables
    !------------------
    type(vmap_t),intent(inout):: Vmap
    !
    ! Local values
    !-------------
    character(len=32):: VmapFunction
    character(len=32):: VmapOpt

    VmapFunction = trim(Vmap%VmapFunction)
    VmapOpt      = trim(Vmap%VmapOpt)

    ! This should be cleaned up ASAP, to use function pointers
    !-------------------------------------------------------------
    if(trim(VmapFunction).eq.'vmap_index') then
      ! Maps that use vertical indices for mapping
      !----------------------------------------------------
      init_Vmap = vmap_index(Vmap)
    elseif(trim(VmapFunction).eq.'vmap_pressure') then
      ! Maps that use P or ln(P)  for mapping
      !----------------------------------------------------
      init_Vmap = vmap_pressure(Vmap)
    elseif(trim(VmapFunction).eq.'vmap_USER') then
      ! USER defined Maps
      !----------------------------------------------------
      init_Vmap = vmap_USER(Vmap)
    else
      ! ERROR: Unknown FunctionName
      !----------------------------------------------------
      init_Vmap = .false.
    endif

    ! End Routine
    !------------------
    return
  end function init_Vmap
  !================================================================


  !================================================================
  logical function apply_Vmap(Vdat_Out,Vdat_In,DPdesc,I_VmapFunction,        &
                                                      I_VmapOpt     ,VmapData)
!PFC!   logical function apply_Vmap(VmapData,Vgrid_in,Vgrid_out,Vdata_in,Vdata_out) 

    ! 
    ! apply_Vmap:  This is a logical function to use the given Vmap
    !              to Vertically process/map the data from the input 
    !              source to the output grid. The data values are all 
    !              assumed to be real*8.
    !
    !   This is an interface function which calls the specified Function
    !   with the input/output grids and mapping option. 
    !
    !   DEVNOTE: For now the mapping data is passed thru a specific 
    !            data structure (vmap_data_t). For more general options 
    !            this will be something that needs to be sorted out. 
    !            ?? Can I pass VmapData as a generic class [ class(*) ] so 
    !            that the called VmapFunctions can initialze a datastructure
    !            tailored to the function needs? 
    !==============================================================
    !
    ! Passed Variables
    !------------------
    real(r8)         ,intent(inout):: Vdat_Out(:,:)
    real(r8)         ,intent(inout):: Vdat_In (:,:)
    type(PathDesc_t) ,intent(inout):: DPdesc
    character(len=*) ,intent(in   ):: I_VmapFunction
    character(len=*) ,intent(in   ):: I_VmapOpt
    type(vmap_data_t),intent(inout):: VmapData
    !
    ! Local Values
    !--------------
    character(len=32):: VmapFunction
    character(len=32):: VmapOpt

    VmapFunction = trim(I_VmapFunction)
    VmapOpt      = trim(I_VmapOpt)

    ! This should be cleaned up ASAP, to use function pointers
    !-------------------------------------------------------------
    if(trim(VmapFunction).eq.'vmap_index') then
      ! Maps that use vertical indices for mapping
      !----------------------------------------------------
      apply_Vmap = vmap_index(Vdat_Out,Vdat_In,DPdesc,trim(VmapOpt),VmapData)
    elseif(trim(VmapFunction).eq.'vmap_pressure') then
      ! Maps that use P or ln(P)  for mapping
      !----------------------------------------------------
      apply_Vmap = vmap_pressure(Vdat_Out,Vdat_In,DPdesc,trim(VmapOpt),VmapData)
    elseif(trim(VmapFunction).eq.'vmap_USER') then
      ! USER defined Maps
      !----------------------------------------------------
      apply_Vmap = vmap_USER(Vdat_Out,Vdat_In,DPdesc,trim(VmapOpt),VmapData)
    else
      ! ERROR: Unknown FunctionName
      !----------------------------------------------------
      apply_Vmap = .false.
    endif

    ! End Routine
    !------------------
    return
  end function apply_Vmap
  !================================================================


  !================================================================
  logical function vmap_index_chk(Vgrid_out,Vgrid_in,I_VmapOpt,Vmap)
    !
    ! vmap_index_chk: Given input Vgrid_in and output Vgrid_out grids
    !                 check them against the given VmapOpt and 
    !                 its requirements.
    !==============================================================
    !
    ! Passed variables  
    !-----------------
    type(vgrid_t)    ,intent(in):: Vgrid_out
    type(vgrid_t)    ,intent(in):: Vgrid_in
    character(len=*) ,intent(in):: I_VmapOpt
    type(vmap_desc_t),intent(in):: Vmap
    !
    ! Local values
    !---------------
    integer:: ii,jj
    logical:: map_found,grid_found
    character(len=32):: VmapOpt

    VmapOpt = trim(I_VmapOpt)

    ! First check that the in/out grid type match the map
    !------------------------------------------------------
    vmap_index_chk = .true.

    if((trim(Vmap%InType ).ne.'INDEX').or. &
       (trim(Vmap%OutType).ne.'INDEX')     ) then
      vmap_index_chk = .false.
      return
    endif

    ! Verify that the VmapOpt is in the Opts list
    !----------------------------------------------
    map_found = .false.
    do ii=1,Vmap%numMapOpts
      if(trim(VmapOpt).eq.trim(Vmap%Opts(ii))) then
        map_found = .true.
      endif
    end do
    if(.not.map_found) then
      vmap_index_chk = .false.
      return
    endif

    ! See if the Input grid requirements are satisfied
    !--------------------------------------------------
    if(Vmap%nreq_in.gt.0) then
      grid_found = .false.
      do ii=1,Vmap%nreq_in
      do jj=1,Vgrid_in%nreq_var
        if(trim(Vmap%Vreq_in(ii)).eq.trim(Vgrid_in%req_var(jj))) then
          grid_found = .true.
        endif
      end do
      end do
      if(.not.grid_found) then
        vmap_index_chk = .false.
        return
      endif
    endif
    
    ! See if the Output grid requirements are satisfied
    !--------------------------------------------------
    if(Vmap%nreq_out.gt.0) then
      grid_found = .false.
      do ii=1,Vmap%nreq_out
      do jj=1,Vgrid_out%nreq_var
        if(trim(Vmap%Vreq_out(ii)).eq.trim(Vgrid_out%req_var(jj))) then
          grid_found = .true.
        endif
      end do
      end do
      if(.not.grid_found) then
        vmap_index_chk = .false.
        return
      endif
    endif

    ! Check In/Out Vgrid Constraints 
    !---------------------------------
    if(trim(VmapOpt).eq.'INDEX') then
      ! The number of in/out levels must match
      !----------------------------------------
      if(Vgrid_out%nlev.ne.Vgrid_in%nlev) then
        vmap_index_chk = .false.
      endif
    elseif(trim(VmapOpt).eq.'INDEX_SHIFT') then
      ! The number of in/out levels don;t need to match
      !-------------------------------------------------
      vmap_index_chk = .true.
    endif

    ! End Routine
    !-----------
    return
  end function vmap_index_chk
  !================================================================


  !================================================================
  logical function vmap_index_init(Vmap)
    ! 
    ! vmap_index_init: Initialize the VMAP values needed to implement vertical 
    !                  interpolation of input data for the gridpoints contained 
    !                  in this PE.
    !
    !       Initialize the (Vmap%) datastructure values needed for the specified
    !       vertical interpolation method.
    !
    !==============================================================
    !
    ! Passed Variables
    !------------------
    type(vmap_t),intent(inout):: Vmap
    !
    ! Local Values
    !--------------
    character(len=32):: VmapOpt

    ! Set the Vertical Map option in VmapData
    !-----------------------------------------
    VmapOpt = trim(Vmap%VmapOpt)
!PFC???    VmapData%VmapOpt = VmapOpt
    vmap_index_init = .true.

    ! End Routine
    !-----------
    return
  end function vmap_index_init
  !================================================================


  !================================================================
  logical function vmap_index_apply(Vdat_Out,Vdat_In,DPdesc,I_VmapOpt,VmapData)
    ! 
    ! vmap_index_apply:  This is an index map in which values are copied 
    !                    thru level by level without any modification.
    !
    !   Constraints: For the INDEX_SHIFT option, the number of in/out levels 
    !                can differ. Only a sub-set of the levels are used using the 
    !                logic from the SD nudging method. 
    !                For the INDEX option, the number of in/out vertical levels 
    !                must be equal.
    !==============================================================
    !
    ! Passed variables
    !-------------------
    real(r8)         ,intent(inout):: Vdat_Out(:,:)
    real(r8)         ,intent(inout):: Vdat_In (:,:)
    type(PathDesc_t) ,intent(inout):: DPdesc
    character(len=*) ,intent(in   ):: I_VmapOpt
    type(vmap_data_t),intent(inout):: VmapData
    !
    ! Local Values
    !---------------
    character(len=32):: VmapOpt
    integer:: ind_lo,ind_hi
    integer:: nmaps,ii,jj

    ! Set local values for readability
    !--------------------------------------------
    VmapOpt = trim(I_VmapOpt)

    ! Select the mapping option
    !----------------------------
    if(trim(VmapOpt).eq.'INDEX_SHIFT') then
      ! Just a straight pass thru, but with index values shifted
      ! in a mannor consistent with the existing SD logic.
      !----------------------------------------------------------
      Vdat_Out(1:DPdesc%ncol,1:DPdesc%nlev_m) = 0._r8
      if(DPdesc%nlev_d .gt. DPdesc%nlev_m) then
        ind_lo = DPdesc%nlev_d-DPdesc%nlev_m + 1
        ind_hi = DPdesc%nlev_d
        Vdat_Out(1:DPdesc%ncol,1:DPdesc%nlev_m) = Vdat_In(1:DPdesc%ncol,ind_lo:ind_hi)
      else
        ind_lo = DPdesc%nlev_m-DPdesc%nlev_d + 1
        ind_hi = DPdesc%nlev_m
        Vdat_Out(1:DPdesc%ncol,ind_lo:ind_hi) = Vdat_In(1:DPdesc%ncol,1:DPdesc%nlev_d)
      endif
    elseif(trim(VmapOpt).eq.'INDEX') then
      ! Just a straight pass thru
      !---------------------------
      Vdat_Out(1:DPdesc%ncol,1:DPdesc%nlev_m) = Vdat_In(1:DPdesc%ncol,1:DPdesc%nlev_d)
    else
      ! ERROR: Unknown VmapOpt
      !------------------------
      write(iulog,*) '  ERROR: vmap_index_apply() unknown VmapOpt=',trim(VmapOpt)
      call endrun(' vmap_index_apply() ERROR')
    endif

    vmap_index_apply = .true.
  
    ! End Routine
    !-----------
    return
  end function vmap_index_apply
  !================================================================




  !================================================================
  logical function vmap_pressure_chk(Vgrid_out,Vgrid_in,I_VmapOpt,Vmap)
    !
    ! vmap_pressure_chk: Given input Vgrid_in and output Vgrid_out grids
    !                    check them against the given VmapOpt and 
    !                    its requirements.
    !==============================================================
    !
    ! Passed variables  
    !-----------------
    type(vgrid_t)    ,intent(in):: Vgrid_out
    type(vgrid_t)    ,intent(in):: Vgrid_in
    character(len=*) ,intent(in):: I_VmapOpt
    type(vmap_desc_t),intent(in):: Vmap
    !
    ! Local values
    !---------------
    integer:: ii,jj
    logical:: map_found,grid_found
    character(len=32):: VmapOpt

    VmapOpt = trim(I_VmapOpt)

    ! First check that the in/out grid type match the map
    !------------------------------------------------------
    vmap_pressure_chk = .true.

    if((trim( Vgrid_in%Type).ne.trim(Vmap%InType )).or. &
       (trim(Vgrid_out%Type).ne.trim(Vmap%OutType))     ) then
      vmap_pressure_chk = .false.
      return
    endif

    ! Verify that the VmapOpt is in the Opts list
    !----------------------------------------------
    map_found = .false.
    do ii=1,Vmap%numMapOpts
      if(trim(VmapOpt).eq.trim(Vmap%Opts(ii))) then
        map_found = .true.
      endif
    end do
    if(.not.map_found) then
      vmap_pressure_chk = .false.
      return
    endif

    ! See if the Input grid requirements are satisfied
    !--------------------------------------------------
    if(Vmap%nreq_in.gt.0) then
      grid_found = .false.
      do ii=1,Vmap%nreq_in
      do jj=1,Vgrid_in%nreq_var
        if(trim(Vmap%Vreq_in(ii)).eq.trim(Vgrid_in%req_var(jj))) then
          grid_found = .true.
        endif
      end do
      end do
      if(.not.grid_found) then
        vmap_pressure_chk = .false.
        return
      endif
    endif
    
    ! See if the Output grid requirements are satisfied
    !--------------------------------------------------
    if(Vmap%nreq_out.gt.0) then
      grid_found = .false.
      do ii=1,Vmap%nreq_out
      do jj=1,Vgrid_out%nreq_var
        if(trim(Vmap%Vreq_out(ii)).eq.trim(Vgrid_out%req_var(jj))) then
          grid_found = .true.
        endif
      end do
      end do
      if(.not.grid_found) then
        vmap_pressure_chk = .false.
        return
      endif
    endif

    ! Check In/Out Vgrid Constraints 
    !---------------------------------




    ! End Routine
    !-----------
    return
  end function vmap_pressure_chk
  !================================================================


  !================================================================
!PFC?  logical function vmap_pressure_init(DPdesc,Vmap,Vgrid_in,Vgrid_out, ) 
  logical function vmap_pressure_init(Vmap)
    ! 
    ! vmap_pressure_init: Initialize the VMAP values needed to implement vertical 
    !                     interpolation of input data for the gridpoints contained 
    !                     in this PE.
    !
    !       Initialize the (Vmap%) datastructure values needed for the specified
    !       vertical interpolation method.
    !
    !==============================================================
    !
    ! Passed Variables
    !------------------
    type(vmap_t),intent(inout):: Vmap
    !
    ! Local Values
    !--------------
    character(len=32):: VmapOpt

    VmapOpt = trim(Vmap%VmapOpt)


    ! Given an input/output grid and a matching Vmap, allocate and 
    ! initialize all needed values in the map datastructure
    !--------------------------------------------------------------


    ! 




    vmap_pressure_init = .true.

    ! End Routine
    !-----------
    return
  end function vmap_pressure_init
  !================================================================


  !================================================================
  logical function vmap_pressure_apply(Vdata_out,Vdata_in,DPdesc,I_VmapOpt,VmapData)
!  logical function vmap_pressure_apply(DPdesc,Vmap,Vdata_in,Vdata_out) 
    ! 
    ! vmap_pressure_apply:  This is an pressure based map in which values are 
    !                       interpolated using P or log(P) as the independent 
    !                       variable. This map accomodates both Hybrid-P as 
    !                       well as Std-Pressure vertical grids.
    !
    !   Constraints: 
    !
    !==============================================================
    !
    ! Passed variables
    !-------------------
    type(PathDesc_t),intent(inout):: DPdesc
    real(r8)        ,intent(inout):: Vdata_in (:,:)
    real(r8)        ,intent(inout):: Vdata_out(:,:)
    character(len=*) ,intent(in   ):: I_VmapOpt
    type(vmap_data_t),intent(inout):: VmapData
    !
    ! Local Values
    !---------------
    integer:: nn,ii,levo,levi

    ! Loop over the number of columns
    !--------------------------------
!    do nn=1,DPdesc%ncol
!
!      ! Loop over output levels 
!      !-------------------------
!      do levo=1,Vmap%Vgrid_out%nlev
!
!        ! Apply interp wgts to input to Vdata_in()
!        !-----------------------------------------
!        Vdata_out(nn,levo) = 0._r8
!        do ii=1,Vmap%VmapData%NumWgts(levo)
!          levi= ii + Vmap%VmapData%Ilevlo(levo) - 1
!          Vdata_out(nn,levo) = Vdata_out(nn,levo)                              &
!                             + Vdata_in (nn,levi)*Vmap%VmapData%Wgts(ii,levo,nn)
!        end do
!
!      end do 
!    end do ! nn=1,DPdesc%ncol

    vmap_pressure_apply = .true.
  
    ! End Routine
    !-----------
    return
  end function vmap_pressure_apply
  !================================================================


  !================================================================
  logical function vmap_USER_chk(Vgrid_out,Vgrid_in,I_VmapOpt,Vmap)
    !
    ! vmap_USER_chk: Given input Vgrid_in and output Vgrid_out grids
    !                check them against the given VmapOpt and 
    !                its requirements.
    !
    !       This is a TEMPLATE/STUB to help users wishing to implement 
    !       their own mapping method.
    !
    !==============================================================
    !
    ! Passed variables  
    !-----------------
    type(vgrid_t)    ,intent(in):: Vgrid_out
    type(vgrid_t)    ,intent(in):: Vgrid_in
    character(len=*) ,intent(in):: I_VmapOpt
    type(vmap_desc_t),intent(in):: Vmap
    !
    ! Local values
    !---------------
    integer:: ii,jj
    logical:: map_found,grid_found
    character(len=32):: VmapOpt

    VmapOpt = trim(I_VmapOpt)

    ! First check that the in/out grid type match the map
    !------------------------------------------------------
    vmap_USER_chk = .true.

    if((trim( Vgrid_in%Type).ne.trim(Vmap%InType )).or. &
       (trim(Vgrid_out%Type).ne.trim(Vmap%OutType))     ) then
      vmap_USER_chk = .false.
      return
    endif

    ! Verify that the VmapOpt is in the Opts list
    !----------------------------------------------
    map_found = .false.
    do ii=1,Vmap%numMapOpts
      if(trim(VmapOpt).eq.trim(Vmap%Opts(ii))) then
        map_found = .true.
      endif
    end do
    if(.not.map_found) then
      vmap_USER_chk = .false.
      return
    endif

    ! See if the Input grid requirements are satisfied
    !--------------------------------------------------
    if(Vmap%nreq_in.gt.0) then
      grid_found = .false.
      do ii=1,Vmap%nreq_in
      do jj=1,Vgrid_in%nreq_var
        if(trim(Vmap%Vreq_in(ii)).eq.trim(Vgrid_in%req_var(jj))) then
          grid_found = .true.
        endif
      end do
      end do
      if(.not.grid_found) then
        vmap_USER_chk = .false.
        return
      endif
    endif
    
    ! See if the Output grid requirements are satisfied
    !--------------------------------------------------
    if(Vmap%nreq_out.gt.0) then
      grid_found = .false.
      do ii=1,Vmap%nreq_out
      do jj=1,Vgrid_out%nreq_var
        if(trim(Vmap%Vreq_out(ii)).eq.trim(Vgrid_out%req_var(jj))) then
          grid_found = .true.
        endif
      end do
      end do
      if(.not.grid_found) then
        vmap_USER_chk = .false.
        return
      endif
    endif

    ! Check In/Out Vgrid Constraints 
    !---------------------------------




    ! End Routine
    !-----------
    return
  end function vmap_USER_chk
  !================================================================


  !================================================================
  logical function vmap_USER_init(Vmap)
    ! 
    ! vmap_USER_init: Initialize the VMAP values needed to implement 
    !                 some USER defined vertical interpolation of input 
    !                 data for the gridpoints contained in this PE.
    !
    !       This is a TEMPLATE/STUB to help users wishing to implement 
    !       their own mapping method.
    !
    !       Initialize the (Vmap%) datastructure values needed for the 
    !       specified vertical interpolation method.
    !
    !==============================================================
    !
    ! Passed Variables
    !------------------
    type(vmap_t),intent(inout):: Vmap
    !
    ! Local Values
    !--------------
    character(len=32):: VmapOpt

    VmapOpt = trim(Vmap%VmapOpt)

       ! ** DEVNOTE: NEED To add prereq info and figure out how to initialize 
       ! **          when PHIS values need to be mapped and save in the VmapData
       ! **          structure  

    vmap_USER_init = .true.

    ! End Routine
    !-----------
    return
  end function vmap_USER_init
  !================================================================


  !================================================================
  logical function vmap_USER_apply(Vdat_Out,Vdat_In,DPdesc,I_VmapOpt,VmapData)
    ! 
    ! vmap_USER_apply:  This is an USER based map in which values are 
    !                   interpolated using a user defined mapping.
    !
    !       This is a TEMPLATE/STUB to help users wishing to implement 
    !       their own mapping method.
    !
    !   Constraints: 
    !
    !==============================================================
    !
    ! Passed variables
    !-------------------
    real(r8)         ,intent(inout):: Vdat_Out(:,:)
    real(r8)         ,intent(inout):: Vdat_In (:,:)
    type(PathDesc_t) ,intent(inout):: DPdesc
    character(len=*) ,intent(in   ):: I_VmapOpt
    type(vmap_data_t),intent(inout):: VmapData
    !
    ! Local Values
    !---------------
    character(len=32):: VmapOpt

    ! Set local values for readability
    !--------------------------------------------
    VmapOpt = trim(I_VmapOpt)



    vmap_USER_apply = .true.
  
    ! End Routine
    !-----------
    return
  end function vmap_USER_apply
  !================================================================

end module Vmaps_mod
