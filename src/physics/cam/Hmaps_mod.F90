module Hmaps_mod
!===========================================================================
!
! MODULE: Hmaps_mod
!
! DESCRIPTION: Manage the mapping of data values from a given dataset to 
!              model data arrays.
!
!===========================================================================
  ! Useful Modules
  !---------------------
  use pio,           only: pio_offset_kind, file_desc_t, var_desc_t, pio_double,     &
                           pio_inq_dimid, pio_max_var_dims, io_desc_t, pio_setframe, &
                           iosystem_desc_t, pio_initdecomp, pio_inq_varid, pio_read_darray
  use shr_pio_mod,   only: shr_pio_getiosys
  use cam_instance,  only: atm_id
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use spmd_utils,    only: masterproc
  use cam_abortutils,only: endrun
  use cam_logfile,   only: iulog
  use Dataset_mod,   only: dataset_t, dataset_var_t, inq_dataset_ncid
  use Dataset_mod,   only: inq_dataset_ncid_nf90  ! PFC DIAG Temp function to test PIO ERROR with nf90
  use PathDesc_type ,only: PathDesc_t

  ! Set all Global values and routines to private by default 
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private
  save

  public :: packed_hgrid_t
  public :: dataset_hgrid_t 
  public :: hmap_rect_data_t          
  public :: hmaps_t          
  public :: known_hgrids_t
  public :: known_hmaps_t

  public :: init_Hmaps_mod
  public :: inq_Hmaps
  public :: print_known_Hmaps
  public :: copy_Hgrid
  private:: copy_Hgrid_m
  private:: copy_Hgrid_d
  public :: deallocate_Hgrid
  private:: deallocate_Hgrid_m
  private:: deallocate_Hgrid_d
  private:: ind_sort

  public ::   chk_Hmap
  public ::  init_Hmap
  public :: apply_Hmap

  public :: hmap_rect
  private:: hmap_rect_chk
  private:: hmap_rect_init
  private:: hmap_rect_apply

  public :: hmap_native
  private:: hmap_native_chk
  private:: hmap_native_init
  private:: hmap_native_apply

  public :: hmap_USER
  private:: hmap_USER_chk
  private:: hmap_USER_init
  private:: hmap_USER_apply

  ! Interfaces
  !-----------------
  interface hmap_rect
    procedure:: hmap_rect_chk
    procedure:: hmap_rect_init
    procedure:: hmap_rect_apply
  end interface

  interface hmap_native
    procedure:: hmap_native_chk
    procedure:: hmap_native_init
    procedure:: hmap_native_apply
  end interface

  interface hmap_USER
    procedure:: hmap_USER_chk
    procedure:: hmap_USER_init
    procedure:: hmap_USER_apply
  end interface

  interface copy_Hgrid
    procedure:: copy_Hgrid_m
    procedure:: copy_Hgrid_d
  end interface

  interface deallocate_Hgrid
    procedure:: deallocate_Hgrid_m
    procedure:: deallocate_Hgrid_d
  end interface

  ! Type Definitions
  !-------------------
  type packed_hgrid_t
    character(len=32)   :: name
    character(len=32)   :: Type
    character(len=32)   :: PackType
    character(len=32)   :: model_gridname
    integer             :: model_gridid
    integer             :: ncol
    real(r8),allocatable:: Lon (:)
    real(r8),allocatable:: Lat (:)
    integer ,allocatable:: Xmap(:)
    integer ,allocatable:: Ymap(:)
    character(len=8)    :: Hdim1name,Hdim2name
    integer             :: Hdim1b,Hdim1e,Hdim2b,Hdim2e
    integer             :: Hdim1len,Hdim2len
  end type packed_hgrid_t

  type dataset_hgrid_t 
    character(len=32)   :: name
    character(len=32)   :: Type
    integer             :: datasetID
    integer             :: nhdim
    character(len=32)   :: dimName(2)
    integer             :: dimID  (2)
    integer             :: dimLen (2)
    integer             :: nlon
    integer             :: nlat
    real(r8),allocatable:: Lon (:)
    real(r8),allocatable:: Lat (:)
  end type dataset_hgrid_t 

  type hmap_rect_data_t
    integer             :: nlon
    integer             :: nlat
    real(r8),allocatable:: lon(:)
    real(r8),allocatable:: lat(:)
    character(len=32)   :: HmapOpt
    integer             :: ncol
    integer             :: NbrCount
    integer             :: npts
    integer ,allocatable:: NbrList(:,:)
    real(r8),allocatable:: Vparity(:)
    integer ,allocatable:: lonInd (:)
    integer ,allocatable:: latInd (:)
    real(r8),allocatable:: NbrIwgt(:,:)
    integer ,allocatable:: dof2D  (:)     ! PFC?? Has this outlived its usefulness?
  end type hmap_rect_data_t

  type hmaps_t
    character(len=32)       :: name
    character(len=32)       :: FunctionName
    character(len=32)       :: InType
    character(len=32)       :: OutType
    integer                 :: numMapOpts
    character(len=32),allocatable:: Opts   (:)
  end type hmaps_t

  type known_hgrids_t
    integer:: numDataGrids
    integer:: numModelGrids
    type(dataset_hgrid_t),allocatable:: DataGrid (:)
    type(packed_hgrid_t) ,allocatable:: ModelGrid(:)
  end type known_hgrids_t

  type known_hmaps_t
    integer                  :: numDataMaps
    integer                  :: numModelMaps
    type(hmaps_t),allocatable:: DataMap (:)
    type(hmaps_t),allocatable:: ModelMap(:)
  end type known_hmaps_t

  ! Global Values
  !-----------------
  type(known_hgrids_t):: My_Hgrids
  type(known_hmaps_t ):: My_Hmaps

contains
  !================================================================
  subroutine init_Hmaps_mod
    ! 
    ! init_Hmaps_mod: initialization the internal datastucture
    !                 containing the set of horizontal mappings
    !                 supported by the Hmap_mod module.
    !
    !  We need to initializae a datastructure that identifies the
    !  number of Dataset mapping routines and the number of Model
    !  mapping routines in the module and for each mapping routine 
    !  provide:
    !              Name of routine
    !              Set of Mapping options
    !              Set of Input grids
    !              Set of Output grids
    !
    ! NOTE: The set of identifiers: (INPUT,OUTPUT,MAP_OPT) must be 
    !       unique (e.g. 2 mapping routines cannot duplicate this 
    !                    set of these identifiers)
    !==============================================================
    !
    ! Local Values
    !--------------
    integer:: Didx_LON_LAT,Didx_SLON_LAT,Didx_LON_SLAT,Didx_NCOL
    integer:: Midx_ANY,Midx_PHYSGRID

    ! Define Known DataSet Hgrids
    !-------------------------------
    My_Hgrids%numDataGrids = 4
    allocate(My_Hgrids%DataGrid(My_Hgrids%numDataGrids))
    Didx_LON_LAT = 1
    My_Hgrids%DataGrid(1)%name      = '(lon,lat)'
    My_Hgrids%DataGrid(1)%Type      = 'RECTILINEAR'
    My_Hgrids%DataGrid(1)%nhdim     = 2
    My_Hgrids%DataGrid(1)%dimName(1)= 'lon'
    My_Hgrids%DataGrid(1)%dimName(2)= 'lat'
    My_Hgrids%DataGrid(1)%dimID  (1)= -1
    My_Hgrids%DataGrid(1)%dimID  (2)= -1
    My_Hgrids%DataGrid(1)%dimLen (1)=  0
    My_Hgrids%DataGrid(1)%dimLen (2)=  0
    Didx_SLON_LAT = 2
    My_Hgrids%DataGrid(2)%name      = '(slon,lat)'
    My_Hgrids%DataGrid(2)%Type      = 'RECTILINEAR'
    My_Hgrids%DataGrid(2)%nhdim     = 2
    My_Hgrids%DataGrid(2)%dimName(1)= 'slon'
    My_Hgrids%DataGrid(2)%dimName(2)= 'lat'
    My_Hgrids%DataGrid(2)%dimID  (1)= -1
    My_Hgrids%DataGrid(2)%dimID  (2)= -1
    My_Hgrids%DataGrid(2)%dimLen (1)=  0
    My_Hgrids%DataGrid(2)%dimLen (2)=  0
    Didx_LON_SLAT = 3
    My_Hgrids%DataGrid(3)%name      = '(lon,slat)'
    My_Hgrids%DataGrid(3)%Type      = 'RECTILINEAR'
    My_Hgrids%DataGrid(3)%nhdim     = 2
    My_Hgrids%DataGrid(3)%dimName(1)= 'lon'
    My_Hgrids%DataGrid(3)%dimName(2)= 'slat'
    My_Hgrids%DataGrid(3)%dimID  (1)= -1
    My_Hgrids%DataGrid(3)%dimID  (2)= -1
    My_Hgrids%DataGrid(3)%dimLen (1)=  0
    My_Hgrids%DataGrid(3)%dimLen (2)=  0
    Didx_NCOL = 4
    My_Hgrids%DataGrid(4)%name      = '(ncol)'
    My_Hgrids%DataGrid(4)%Type      = 'SE_GRID'
    My_Hgrids%DataGrid(4)%nhdim     = 1
    My_Hgrids%DataGrid(4)%dimName(1)= 'ncol'
    My_Hgrids%DataGrid(4)%dimName(2)= 'NULL'
    My_Hgrids%DataGrid(4)%dimID  (1)= -1
    My_Hgrids%DataGrid(4)%dimID  (2)= -1
    My_Hgrids%DataGrid(4)%dimLen (1)=  0
    My_Hgrids%DataGrid(4)%dimLen (2)=  1

    My_Hgrids%numModelGrids = 2
    allocate(My_Hgrids%ModelGrid(My_Hgrids%numModelGrids))
    Midx_ANY = 1
    My_Hgrids%ModelGrid(1)%name          = 'ANY'
    My_Hgrids%ModelGrid(1)%model_gridname= 'ANY'
    My_Hgrids%ModelGrid(1)%Type          = 'PACKED_MODEL'
    My_Hgrids%ModelGrid(1)%PackType      = 'ANY'
    My_Hgrids%ModelGrid(1)%ncol          = -1
    Midx_PHYSGRID = 2
    My_Hgrids%ModelGrid(2)%name          = 'physgrid'
    My_Hgrids%ModelGrid(2)%model_gridname= 'physgrid'
    My_Hgrids%ModelGrid(2)%Type          = 'PACKED_MODEL'
    My_Hgrids%ModelGrid(2)%PackType      = 'UNSTRUCTURED'
    My_Hgrids%ModelGrid(2)%ncol          = -1

    ! Define Known DataSet Maps:
    !-------------------------------
    My_Hmaps%numDataMaps = 4
    allocate(My_Hmaps%DataMap(My_Hmaps%numDataMaps))

    ! Hmap for Rectilinear Dataset grids with embedded PIO reads
    !------------------------------------------------------------
    My_Hmaps%DataMap(1)%name         = 'hmap_rect'
    My_Hmaps%DataMap(1)%FunctionName = 'hmap_rect'
    My_Hmaps%DataMap(1)%InType       = 'RECTILINEAR'
    My_Hmaps%DataMap(1)%OutType      = 'PACKED_MODEL'
    My_Hmaps%DataMap(1)%numMapOpts   = 10
    allocate(My_Hmaps%DataMap(1)%Opts(My_Hmaps%DataMap(1)%numMapOpts))
    My_Hmaps%DataMap(1)%Opts(1 )     = 'BILINEAR_INTERP'
    My_Hmaps%DataMap(1)%Opts(2 )     = 'BICUBIC_INTERP'
    My_Hmaps%DataMap(1)%Opts(3 )     = 'BIQUINTIC_INTERP'
    My_Hmaps%DataMap(1)%Opts(4 )     = 'BICUBIC_FILTER'
    My_Hmaps%DataMap(1)%Opts(5 )     = 'BIQUINTIC_FILTER'
    My_Hmaps%DataMap(1)%Opts(6 )     = 'VECT_BILINEAR_INTERP'
    My_Hmaps%DataMap(1)%Opts(7 )     = 'VECT_BICUBIC_INTERP'
    My_Hmaps%DataMap(1)%Opts(8 )     = 'VECT_BIQUINTIC_INTERP'
    My_Hmaps%DataMap(1)%Opts(9 )     = 'VECT_BICUBIC_FILTER'
    My_Hmaps%DataMap(1)%Opts(10)     = 'VECT_BIQUINTIC_FILTER'

    ! Native Hmap from ANY Dataset grids 
    ! pure PIO reads onto model grids.
    !--------------------------------------------------------
    My_Hmaps%DataMap(2)%name         = 'hmap_native'
    My_Hmaps%DataMap(2)%FunctionName = 'hmap_native'
    My_Hmaps%DataMap(2)%InType       = 'ANY'    ! 'RECTILINEAR'
    My_Hmaps%DataMap(2)%OutType      = 'ANY'    ! 'PACKED_MODEL'
    My_Hmaps%DataMap(2)%numMapOpts   = 1
    allocate(My_Hmaps%DataMap(2)%Opts(My_Hmaps%DataMap(2)%numMapOpts))
    My_Hmaps%DataMap(2)%Opts(1)      = 'NATIVE' ! 'IDENTITY'

    ! DEVNOTE:  (NON-FUNCTIONAL ~Template~) for adding a SE grid option
    !-------------------------------------------------------------------
    My_Hmaps%DataMap(3)%name         = 'hmap_se'
    My_Hmaps%DataMap(3)%FunctionName = 'hmap_se'
    My_Hmaps%DataMap(3)%InType       = 'SE_GRID'
    My_Hmaps%DataMap(3)%OutType      = 'PACKED_MODEL'
    My_Hmaps%DataMap(3)%numMapOpts   = 1
    allocate(My_Hmaps%DataMap(3)%Opts(My_Hmaps%DataMap(3)%numMapOpts))
    My_Hmaps%DataMap(3)%Opts(1)      = 'ESMF_OPTS_HERE'

    ! USER Defined Map
    !--------------------------------------------------------
    My_Hmaps%DataMap(4)%name         = 'hmap_USER'
    My_Hmaps%DataMap(4)%FunctionName = 'hmap_USER'
    My_Hmaps%DataMap(4)%InType       = 'USER'
    My_Hmaps%DataMap(4)%OutType      = 'USER'
    My_Hmaps%DataMap(4)%numMapOpts   = 1
    allocate(My_Hmaps%DataMap(4)%Opts(My_Hmaps%DataMap(4)%numMapOpts))
    My_Hmaps%DataMap(4)%Opts(1)      = 'USER'

    ! DEVNOTE:  Model Maps: (NON-FUNCTIONAL ~Template) 
    !           e.g. Zonal Mean mapping from physgrid -> physgrid
    !--------------------------------------------------------------------------
    My_Hmaps%numModelMaps = 1
    allocate(My_Hmaps%ModelMap(My_Hmaps%numModelMaps))
    My_Hmaps%ModelMap(1)%name         = 'hmap_model'
    My_Hmaps%ModelMap(1)%FunctionName = 'hmap_model'
    My_Hmaps%ModelMap(1)%InType       = 'PACKED_MODEL'
    My_Hmaps%ModelMap(1)%OutType      = 'PACKED_MODEL'
    My_Hmaps%ModelMap(1)%numMapOpts   = 1
    allocate(My_Hmaps%ModelMap(1)%Opts(My_Hmaps%ModelMap(1)%numMapOpts))
    My_Hmaps%ModelMap(1)%Opts(1)      = 'ZONAL_MEAN'

    ! End Routine
    !------------------
    return
  end subroutine init_Hmaps_mod
  !=====================================================================


  !================================================================
  subroutine inq_Hmaps(Hgrids,Hmaps)
    ! 
    ! inq_Hmaps: Return Hgrids and Hmaps data strucures
    !            to provide information about this modules 
    !            functionality
    !==============================================================
    !
    ! Passed variables
    !------------------
    type(known_hgrids_t),intent(inout):: Hgrids
    type(known_hmaps_t ),intent(inout):: Hmaps
    !
    ! Local values
    !-------------
    integer:: nn,kk

    ! Return the set of known grids
    !------------------------------
    if(allocated(Hgrids%DataGrid)) deallocate(Hgrids%DataGrid)
    Hgrids%numDataGrids = My_Hgrids%numDataGrids
    allocate(Hgrids%DataGrid(Hgrids%numDataGrids))
    do nn=1,Hgrids%numDataGrids
      Hgrids%DataGrid(nn)%name      = trim(My_Hgrids%DataGrid(nn)%name)
      Hgrids%DataGrid(nn)%Type      = trim(My_Hgrids%DataGrid(nn)%Type)
      Hgrids%DataGrid(nn)%nhdim     =      My_Hgrids%DataGrid(nn)%nhdim
      Hgrids%DataGrid(nn)%dimName(1)= trim(My_Hgrids%DataGrid(nn)%dimName(1))
      Hgrids%DataGrid(nn)%dimName(2)= trim(My_Hgrids%DataGrid(nn)%dimName(2))
      Hgrids%DataGrid(nn)%dimID  (:)=      My_Hgrids%DataGrid(nn)%dimID  (:)
      Hgrids%DataGrid(nn)%dimLen (:)=      My_Hgrids%DataGrid(nn)%dimLen (:)
    end do

    if(allocated(Hgrids%ModelGrid)) deallocate(Hgrids%ModelGrid)
    Hgrids%numModelGrids = My_Hgrids%numModelGrids
    allocate(Hgrids%ModelGrid(Hgrids%numModelGrids))
    do nn=1,Hgrids%numModelGrids
      My_Hgrids%ModelGrid(nn)%name = trim(My_Hgrids%ModelGrid(nn)%name)
      My_Hgrids%ModelGrid(nn)%Type = trim(My_Hgrids%ModelGrid(nn)%Type)
      My_Hgrids%ModelGrid(nn)%ncol =      My_Hgrids%ModelGrid(nn)%ncol
    end do

    ! Return the set of known maps
    !------------------------------
    if(allocated(Hmaps%DataMap)) then
      do nn=1,Hmaps%numDataMaps
        if(allocated(Hmaps%DataMap(nn)%Opts)) deallocate(Hmaps%DataMap(nn)%Opts)
      end do
      deallocate(Hmaps%DataMap)
    endif
    if(allocated(Hmaps%ModelMap)) then
      do nn=1,Hmaps%numModelMaps
        if(allocated(Hmaps%ModelMap(nn)%Opts)) deallocate(Hmaps%ModelMap(nn)%Opts)
      end do
      deallocate(Hmaps%ModelMap)
    endif

    Hmaps%numDataMaps = My_Hmaps%numDataMaps
    allocate(Hmaps%DataMap(Hmaps%numDataMaps))
    do nn=1,Hmaps%numDataMaps
      Hmaps%DataMap(nn)%name         = trim(My_Hmaps%DataMap(nn)%name )
      Hmaps%DataMap(nn)%FunctionName = trim(My_Hmaps%DataMap(nn)%FunctionName)
      Hmaps%DataMap(nn)%InType       = trim(My_Hmaps%DataMap(nn)%InType )
      Hmaps%DataMap(nn)%OutType      = trim(My_Hmaps%DataMap(nn)%OutType)
      Hmaps%DataMap(nn)%numMapOpts   = My_Hmaps%DataMap(nn)%numMapOpts
      allocate(Hmaps%DataMap(nn)%Opts(Hmaps%DataMap(nn)%numMapOpts))
      do kk=1,Hmaps%DataMap(nn)%numMapOpts
        Hmaps%DataMap(nn)%Opts(kk) = trim(My_Hmaps%DataMap(nn)%Opts(kk))
      end do
    end do

    Hmaps%numModelMaps = My_Hmaps%numModelMaps
    allocate(Hmaps%ModelMap(Hmaps%numModelMaps))
    do nn=1,Hmaps%numModelMaps
      Hmaps%ModelMap(nn)%name         = trim(My_Hmaps%ModelMap(nn)%name )
      Hmaps%ModelMap(nn)%FunctionName = trim(My_Hmaps%ModelMap(nn)%FunctionName)
      Hmaps%ModelMap(nn)%InType       = trim(My_Hmaps%ModelMap(nn)%InType )
      Hmaps%ModelMap(nn)%OutType      = trim(My_Hmaps%ModelMap(nn)%OutType)
      Hmaps%ModelMap(nn)%numMapOpts   = My_Hmaps%ModelMap(nn)%numMapOpts
      allocate(Hmaps%ModelMap(nn)%Opts(Hmaps%ModelMap(nn)%numMapOpts))
      do kk=1,Hmaps%ModelMap(nn)%numMapOpts
        Hmaps%ModelMap(nn)%Opts(kk)   = trim(My_Hmaps%ModelMap(nn)%Opts(kk))
      end do
    end do

    ! End Routine
    !------------------
    return
  end subroutine inq_Hmaps
  !=====================================================================


  !=====================================================================
  subroutine print_known_Hmaps
    ! 
    ! print_known_Hmaps:  Print the metadata for the known
    !                     Hgrids and Hmaps between them.
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
      write(iulog,*) 'METADATA for Hmap Module'
      write(iulog,*) '    Number of Known DataSet Hgrids=',My_Hgrids%numDataGrids
      write(iulog,*) '    Number of Known Model   Hgrids=',My_Hgrids%numModelGrids
      write(iulog,*) '    Number of Known DataSet  Hmaps=',My_Hmaps%numDataMaps
      write(iulog,*) '    Number of Known Model    Hmaps=',My_Hmaps%numModelMaps
      write(iulog,*) ' '
      write(iulog,*) '  KNOWN HORIZONTAL DATASET GRIDS:'
      do nn=1,My_Hgrids%numDataGrids
        write(iulog,*) '    nn=',nn,' -------------------------'
        write(iulog,*) '       name= ',trim(My_Hgrids%DataGrid(nn)%name)
        write(iulog,*) '       type= ',trim(My_Hgrids%DataGrid(nn)%Type)
        write(iulog,*) '      nhdim= ',My_Hgrids%DataGrid(nn)%nhdim
        do ii=1,My_Hgrids%DataGrid(nn)%nhdim
          write(iulog,*) '      Dim =',ii,' Name= ',trim(My_Hgrids%DataGrid(nn)%dimName(ii))
        end do
      end do
      write(iulog,*) ' '
      write(iulog,*) '  KNOWN HORIZONTAL MODEL GRIDS:'
      do nn=1,My_Hgrids%numModelGrids
        write(iulog,*) '    nn=',nn,' -------------------------'
        write(iulog,*) '                name= ',trim(My_Hgrids%ModelGrid(nn)%name)
        write(iulog,*) '      model_gridname= ',trim(My_Hgrids%ModelGrid(nn)%model_gridname)
        write(iulog,*) '                type= ',trim(My_Hgrids%ModelGrid(nn)%Type)
        write(iulog,*) '                ncol= ',My_Hgrids%ModelGrid(nn)%ncol
      end do
      write(iulog,*) ' '
      write(iulog,*) '  KNOWN HORIZONTAL DATASET MAPS:'
      do nn=1,My_Hmaps%numDataMaps
        write(iulog,*) '    nn=',nn,' -------------------------'
        write(iulog,*) '      name              = ',trim(My_Hmaps%DataMap(nn)%name   )
        write(iulog,*) '      function_name     = ',trim(My_Hmaps%DataMap(nn)%FunctionName)
        write(iulog,*) '      Input  Grid Type  = ',trim(My_Hmaps%DataMap(nn)%InType )
        write(iulog,*) '      Output Grid Type  = ',trim(My_Hmaps%DataMap(nn)%OutType)
        write(iulog,*) '      number of Map Opts= ',My_Hmaps%DataMap(nn)%numMapOpts
        do ii=1,My_Hmaps%DataMap(nn)%numMapOpts
          write(iulog,*) '      MapOpt=',ii,' Name= ',trim(My_Hmaps%DataMap(nn)%Opts(ii))
        end do
      end do
      write(iulog,*) ' '
      write(iulog,*) '  KNOWN HORIZONTAL MODEL MAPS:'
      do nn=1,My_Hmaps%numModelMaps
        write(iulog,*) '    nn=',nn,' -------------------------'
        write(iulog,*) '      name              = ',trim(My_Hmaps%ModelMap(nn)%name   )
        write(iulog,*) '      function_name     = ',trim(My_Hmaps%ModelMap(nn)%FunctionName)
        write(iulog,*) '      Input  Grid Type  = ',trim(My_Hmaps%ModelMap(nn)%InType )
        write(iulog,*) '      Output Grid Type  = ',trim(My_Hmaps%ModelMap(nn)%OutType)
        write(iulog,*) '      number of Map Opts= ',My_Hmaps%ModelMap(nn)%numMapOpts
        do ii=1,My_Hmaps%ModelMap(nn)%numMapOpts
          write(iulog,*) '      MapOpt=',ii,' Name= ',trim(My_Hmaps%ModelMap(nn)%Opts(ii))
        end do
      end do
      write(iulog,*) '============================================================== '
      write(iulog,*) ' '
    endif

    ! End Routine
    !------------------
    return
  end subroutine print_known_Hmaps
  !=====================================================================


  !================================================================
  subroutine copy_Hgrid_m(OutHgrid_m,InHgrid_m)
    !
    ! copy_Hgrid_m: Copy the contents of the Input model grid
    !               to the Output model grid.
    !==============================================================
    !
    ! Passed variables  
    !-----------------
    type(packed_hgrid_t),intent(in )::  InHgrid_m
    type(packed_hgrid_t),intent(out):: OutHgrid_m

    ! Reallocate arrays as needed and Copy the values
    !-------------------------------------------------
    OutHgrid_m%name           = trim(InHgrid_m%name)
    OutHgrid_m%Type           = trim(InHgrid_m%Type)
    OutHgrid_m%PackType       = trim(InHgrid_m%PackType)
    OutHgrid_m%model_gridname = trim(InHgrid_m%model_gridname)
    OutHgrid_m%model_gridid   =      InHgrid_m%model_gridid
    OutHgrid_m%ncol           =      InHgrid_m%ncol
    OutHgrid_m%Hdim1name      = trim(InHgrid_m%Hdim1name)
    OutHgrid_m%Hdim2name      = trim(InHgrid_m%Hdim2name)
    OutHgrid_m%Hdim1b         =      InHgrid_m%Hdim1b
    OutHgrid_m%Hdim1e         =      InHgrid_m%Hdim1e
    OutHgrid_m%Hdim2b         =      InHgrid_m%Hdim2b
    OutHgrid_m%Hdim2e         =      InHgrid_m%Hdim2e
    OutHgrid_m%Hdim1len       =      InHgrid_m%Hdim1len
    OutHgrid_m%Hdim2len       =      InHgrid_m%Hdim2len
    if(allocated(OutHgrid_m%Lon) ) deallocate(OutHgrid_m%Lon) 
    if(allocated(OutHgrid_m%Lat) ) deallocate(OutHgrid_m%Lat) 
    if(allocated(OutHgrid_m%Xmap)) deallocate(OutHgrid_m%Xmap) 
    if(allocated(OutHgrid_m%Ymap)) deallocate(OutHgrid_m%Ymap) 
    allocate(OutHgrid_m%Lon (OutHgrid_m%ncol))
    allocate(OutHgrid_m%Lat (OutHgrid_m%ncol))
    allocate(OutHgrid_m%Xmap(OutHgrid_m%ncol))
    allocate(OutHgrid_m%Ymap(OutHgrid_m%ncol))
    OutHgrid_m%Lon (:) = InHgrid_m%Lon (:)
    OutHgrid_m%Lat (:) = InHgrid_m%Lat (:)
    OutHgrid_m%Xmap(:) = InHgrid_m%Xmap(:)
    OutHgrid_m%Ymap(:) = InHgrid_m%Ymap(:)

    ! End Routine
    !-----------
    return
  end subroutine copy_Hgrid_m
  !================================================================


  !================================================================
  subroutine copy_Hgrid_d(OutHgrid_d,InHgrid_d)
    !
    ! copy_Hgrid_d: Copy the contents of the Input model grid
    !               to the Output model grid.
    !==============================================================
    !
    ! Passed variables  
    !-----------------
    type(dataset_hgrid_t),intent(in )::  InHgrid_d
    type(dataset_hgrid_t),intent(out):: OutHgrid_d

    ! Reallocate arrays as needed and Copy the values
    !-------------------------------------------------
    OutHgrid_d%name       = trim(InHgrid_d%name)
    OutHgrid_d%Type       = trim(InHgrid_d%Type)
    OutHgrid_d%datasetID  =      InHgrid_d%datasetID
    OutHgrid_d%nhdim      =      InHgrid_d%nhdim
    OutHgrid_d%dimName(1) = trim(InHgrid_d%dimName(1))
    OutHgrid_d%dimName(2) = trim(InHgrid_d%dimName(2))
    OutHgrid_d%dimID  (1) =      InHgrid_d%dimID  (1)
    OutHgrid_d%dimID  (2) =      InHgrid_d%dimID  (2)
    OutHgrid_d%dimLen (1) =      InHgrid_d%dimLen (1)
    OutHgrid_d%dimLen (2) =      InHgrid_d%dimLen (2)
    OutHgrid_d%nlon       =      InHgrid_d%nlon   
    OutHgrid_d%nlat       =      InHgrid_d%nlat
    if(allocated(OutHgrid_d%Lon)) deallocate(OutHgrid_d%Lon) 
    if(allocated(OutHgrid_d%Lat)) deallocate(OutHgrid_d%Lat) 
    OutHgrid_d%Lon(:)     = InHgrid_d%Lon(:)
    OutHgrid_d%Lat(:)     = InHgrid_d%Lat(:)

    ! End Routine
    !-----------
    return
  end subroutine copy_Hgrid_d
  !================================================================


  !================================================================
  subroutine deallocate_Hgrid_m(Hgrid_m)
    !
    ! deallocate_Hgrid_m: Deallocate Hgrid_m arrays
    !==============================================================
    !
    ! Passed variables  
    !-----------------
    type(packed_hgrid_t),intent(inout):: Hgrid_m

    ! Reallocate arrays as needed and Copy the values
    !-------------------------------------------------
    if(allocated(Hgrid_m%Lon) ) deallocate(Hgrid_m%Lon) 
    if(allocated(Hgrid_m%Lat) ) deallocate(Hgrid_m%Lat) 
    if(allocated(Hgrid_m%Xmap)) deallocate(Hgrid_m%Xmap) 
    if(allocated(Hgrid_m%Ymap)) deallocate(Hgrid_m%Ymap) 
 
    ! End Routine
    !-----------
    return
  end subroutine deallocate_Hgrid_m
  !================================================================


  !================================================================
  subroutine deallocate_Hgrid_d(Hgrid_d)
    !
    ! deallocate_Hgrid_d: Deallocate Hgrid_d arrays
    !==============================================================
    !
    ! Passed variables  
    !-----------------
    type(dataset_hgrid_t),intent(inout):: Hgrid_d

    ! Deallocate arrays
    !--------------------
    if(allocated(Hgrid_d%Lon)) deallocate(Hgrid_d%Lon) 
    if(allocated(Hgrid_d%Lat)) deallocate(Hgrid_d%Lat) 

    ! End Routine
    !-----------
    return
  end subroutine deallocate_Hgrid_d
  !================================================================


  !================================================================
  subroutine ind_sort(ind_val,ind_map,ind_num)
    !
    ! ind_sort: A sort routine which sorts the given ind_val's 
    !           into ascending order.
    !==========================================================
    !
    ! Passed Values
    !--------------------
    integer ind_num
    integer ind_map(ind_num)
    integer ind_val(ind_num)
    !
    ! Local values
    !---------------
    integer ii,jj,ll,ir,ind,qq

    ! initialize the index array.
    !----------------------------
    do jj=1,ind_num
      ind_map(jj)=jj
    end do
    ll=ind_num/2+1
    ir=ind_num

    ! Do the sort
    !------------
 10 continue
      if(ll.gt.1) then
        ll=ll-1
        ind=ind_map(ll)
        qq=ind_val(ind)
      else
        ind=ind_map(ir)
        qq=ind_val(ind)
        ind_map(ir)=ind_map(1)
        ir=ir-1
        if(ir.eq.1) then
          ind_map(1)=ind
          goto 90
        endif
      endif
      ii=ll
      jj=ll+ll
 20 continue
      if(jj.le.ir) then
        if(jj.lt.ir) then
          if(ind_val(ind_map(jj)).lt.ind_val(ind_map(jj+1))) jj=jj+1
        endif
        if(qq.lt.ind_val(ind_map(jj))) then
          ind_map(ii)=ind_map(jj)
          ii=jj
          jj=jj+jj
        else
          jj=ir+1
        endif
        goto 20
      endif
      ind_map(ii)=ind
      goto 10
 90 continue

    ! End Routine
    !-----------
    return
  end subroutine ind_sort
  !================================================================


  !================================================================
  logical function chk_Hmap(Hgrid_Out,Hgrid_In,I_HmapOpt,known_Hmap)
    ! 
    ! chk_Hmap:  This is a logical function to test if the input/output
    !            horizontal grids satisfiy the requirements of the 
    !            given Map name for the given known Hmap template.
    !
    !   This is an interface function which calls the Function
    !   specified in the known_Hmap structure to check the output 
    !   model grid and the input dataset grid for the given map 
    !   option name.
    !==============================================================
    !
    ! Passed Variables
    !------------------
    type(packed_hgrid_t) ,intent(in):: Hgrid_Out
    type(dataset_hgrid_t),intent(in):: Hgrid_In
    character(len=*)     ,intent(in):: I_HmapOpt
    type(hmaps_t)        ,intent(in):: known_Hmap
    !
    ! Local Value
    !------------
    character(len=32):: HmapOpt

    HmapOpt = '                                '
    HmapOpt = trim(I_HmapOpt)

    ! DEVNOTE?: Do the compilers used at NCAR support function pointers yet? 
    !            It has been >18 YEARS since the 2003 Standard,  Just sayin..
    !-------------------------------------------------------------------------
    if(trim(known_Hmap%FunctionName).eq.'hmap_rect') then
      ! Maps to read/interpolate rectilinear grid datasets
      !----------------------------------------------------
      chk_Hmap = hmap_rect(Hgrid_Out,Hgrid_In,trim(HmapOpt),known_Hmap)
    elseif(trim(known_Hmap%FunctionName).eq.'hmap_native') then
      ! No Mapping, read data on NATIVE model grids.
      !----------------------------------------------------
      chk_Hmap = hmap_native(Hgrid_Out,Hgrid_In,trim(HmapOpt),known_Hmap)
    elseif(trim(known_Hmap%FunctionName).eq.'hmap_USER') then
      ! USER defined Maps
      !----------------------------------------------------
      chk_Hmap = hmap_USER(Hgrid_Out,Hgrid_In,trim(HmapOpt),known_Hmap)
    else
      ! ERROR: Unknown FunctionName
      !----------------------------------------------------
      chk_Hmap = .false.
    endif

    ! End Routine
    !------------------
    return
  end function chk_Hmap
  !================================================================


  !================================================================
  logical function init_Hmap(Hgrid_Out,Hgrid_In,I_HmapFunction,I_HmapOpt,HmapData)
    ! 
    ! init_Hmap:  This is a logical function to fill-in/initialize the 
    !             HmapData data structure with info for reading and 
    !             mapping data from the input Hgrid to the output
    !             Hgrid for the given mapping option.
    !
    !   This is an interface function which calls the specified Function
    !   with the input/output grids and mapping option. 
    !
    !   DEVNOTE: For now the mapping data is passed thru a specific 
    !            data structure (hmap_rect_data_t). For more general options 
    !            this will be something that needs to be sorted out. 
    !            ?? Can I pass HmapData as a generic class [ class(*) ] so 
    !            that the called HmapFunctions can initialze a datastructure
    !            tailored to the function needs? 
    !==============================================================
    !
    ! Passed Variables
    !------------------
    type(packed_hgrid_t)  ,intent(in   ):: Hgrid_Out
    type(dataset_hgrid_t) ,intent(in   ):: Hgrid_In
    character(len=*)      ,intent(in   ):: I_HmapFunction
    character(len=*)      ,intent(in   ):: I_HmapOpt
    type(hmap_rect_data_t),intent(inout):: HmapData
    !
    ! Local Values
    !---------------
    character(len=32):: HmapFunction
    character(len=32):: HmapOpt

    HmapFunction = trim(I_HmapFunction)
    HmapOpt      = trim(I_HmapOpt)

    ! DEVNOTE?: Do the compilers used at NCAR support function pointers yet? 
    !            It has been >18 YEARS since the 2003 Standard,  Just sayin..
    !-------------------------------------------------------------
    if(trim(HmapFunction).eq.'hmap_rect') then
      ! Maps to read/interpolate rectilinear grid datasets
      !----------------------------------------------------
      init_Hmap = hmap_rect(Hgrid_Out,Hgrid_In,trim(HmapOpt),HmapData)
    elseif(trim(HmapFunction).eq.'hmap_native') then
      ! No Mapping, read data on NATIVE model grids.
      !----------------------------------------------------
      init_Hmap = hmap_native(Hgrid_Out,Hgrid_In,trim(HmapOpt),HmapData)
    elseif(trim(HmapFunction).eq.'hmap_USER') then
      ! USER defined Maps
      !----------------------------------------------------
      init_Hmap = hmap_USER(Hgrid_Out,Hgrid_In,trim(HmapOpt),HmapData)
    else
      ! ERROR: Unknown FunctionName
      !----------------------------------------------------
      init_Hmap = .false.
    endif

    ! End Routine
    !------------------
    return
  end function init_Hmap
  !================================================================


  !================================================================
  logical function apply_Hmap(Hout,DPdesc,I_HmapFunction,I_HmapOpt,HmapData)
    ! 
    ! apply_Hmap:  This is a logical function to use the given Hmap
    !              to Horizontally process/map the data for the given 
    !              variable index and timelevel from the input source to 
    !              the output grid. The data values are all assumed to be real*8.
    !              if (nlev=1) then that data is processed as HORIZ_ONLY. 
    !
    !   This is an interface function which calls the specified Function
    !   with the input/output grids and mapping option. 
    !
    !   DEVNOTE: For now the mapping data is passed thru a specific 
    !            data structure (hmap_rect_data_t). For more general options 
    !            this will be something that needs to be sorted out. 
    !            ?? Can I pass HmapData as a generic class [ class(*) ] so 
    !            that the called HmapFunctions can initialze a datastructure
    !            tailored to the function needs? 
    !==============================================================
    !
    ! Passed Variables
    !------------------
    real(r8)              ,intent(out  ):: Hout(:,:)
    type(PathDesc_t)      ,intent(inout):: DPdesc
    character(len=*)      ,intent(in   ):: I_HmapFunction
    character(len=*)      ,intent(in   ):: I_HmapOpt
    type(hmap_rect_data_t),intent(inout):: HmapData
    !
    ! Local Values
    !--------------
    character(len=32):: HmapFunction
    character(len=32):: HmapOpt

    HmapFunction = trim(I_HmapFunction)
    HmapOpt      = trim(I_HmapOpt)

    ! DEVNOTE?: Do the compilers used at NCAR support function pointers yet? 
    !            It has been >18 YEARS since the 2003 Standard,  Just sayin..
    !-------------------------------------------------------------
    if(trim(HmapFunction).eq.'hmap_rect') then
      ! Maps to read/interpolate rectilinear grid datasets
      !----------------------------------------------------
      apply_Hmap = hmap_rect(Hout,DPdesc,trim(HmapOpt),HmapData)
    elseif(trim(HmapFunction).eq.'hmap_native') then
      ! No Mapping, read data on NATIVE model grids.
      !----------------------------------------------------
      apply_Hmap = hmap_native(Hout,DPdesc,trim(HmapOpt),HmapData)
    elseif(trim(HmapFunction).eq.'hmap_USER') then
      ! USER defined Maps
      !----------------------------------------------------
      apply_Hmap = hmap_USER(Hout,DPdesc,trim(HmapOpt),HmapData)
    else
      ! ERROR: Unknown FunctionName
      !----------------------------------------------------
      apply_Hmap = .false.
    endif

    ! End Routine
    !------------------
    return
  end function apply_Hmap
  !================================================================


  !================================================================
  logical function hmap_rect_chk(Hgrid_m,Hgrid_d,I_HmapOpt,DataMap)
    !
    ! hmap_rect_chk: Check the input Hgrid_d and output Hgrid_m grids
    !                check them against the given HmapOpt for 
    !                datasets with rectilinear grids.
    !==============================================================
    !
    ! Passed variables  
    !-----------------
    type(packed_hgrid_t) ,intent(in):: Hgrid_m
    type(dataset_hgrid_t),intent(in):: Hgrid_d
    character(len=*)     ,intent(in):: I_HmapOpt
    type(hmaps_t)        ,intent(in):: DataMap
    !
    ! Local Values
    !---------------
    character(len=32):: HmapOpt

    HmapOpt = trim(I_HmapOpt)

    ! First check that the in/out grid type match the map
    !------------------------------------------------------
    hmap_rect_chk = .true.

    if((trim(Hgrid_d%Type).ne.trim(DataMap%InType )).or. &
       (trim(Hgrid_m%Type).ne.trim(DataMap%OutType))     ) then
      hmap_rect_chk = .false.
      return
    endif

    ! Now check everything else
    !--------------------------

 
    ! End Routine
    !-----------
    return
  end function hmap_rect_chk
  !================================================================


  !================================================================
  logical function hmap_rect_init(Hgrid_m,Hgrid_d,I_HmapOpt,HmapData)
    ! 
    ! hmap_rect_init: Initialize the HMAP values needed to implement horizontal 
    !                 interpolation of input data for the gridpoints contained 
    !                 in this PE for a rectilinear DataSet grid.
    !
    !       Initialize the (Hmap%) datastructure values needed for the specified
    !       horizontal interpolation method.
    !
    !==============================================================
    !
    ! Passed Variables
    !------------------
    type(packed_hgrid_t)  ,intent(in   ):: Hgrid_m
    type(dataset_hgrid_t) ,intent(in   ):: Hgrid_d
    character(len=*)      ,intent(in   ):: I_HmapOpt
    type(hmap_rect_data_t),intent(inout):: HmapData
    !
    ! Local Values
    !--------------
    real(r8),allocatable:: X_lon    (:)
    real(r8),allocatable:: X_lat    (:)
    integer ,allocatable:: X_lonLo  (:)
    integer ,allocatable:: X_lonHi  (:)
    integer ,allocatable:: X_latLo  (:)
    integer ,allocatable:: X_latHi  (:)
    real(r8),allocatable:: X_lonFrac(:)
    real(r8),allocatable:: X_latFrac(:)
    real(r8),allocatable:: X_lonSize(:)
    real(r8),allocatable:: X_latSize(:)
    real(r8),allocatable:: D_Vparity(:)
    integer ,allocatable:: D_lonInd (:)
    integer ,allocatable:: D_latInd (:)
    integer ,allocatable:: D_List   (:,:)
    real(r8),allocatable:: X_Vparity(:)
    integer ,allocatable:: X_lonInd (:)
    integer ,allocatable:: X_latInd (:)
    integer ,allocatable:: X_List   (:,:)
    integer ,allocatable:: indMap   (:)

    real(r8):: LonMin,LonMax,LatMin,LatMax
    logical :: LatUP,LatPOLE,LonROT
    logical :: newPt
    integer :: Nhalo
    integer :: lonLo,lonHi,latLo,latHi
    integer :: npts,nind
    integer :: nh,nc,nn,kk,kk_old
    integer :: ii,jj

    integer             :: M_ncol
    real(r8),allocatable:: M_lon(:)
    real(r8),allocatable:: M_lat(:)
    integer             :: DS_nlon
    integer             :: DS_nlat
    real(r8),allocatable:: DS_lon(:)
    real(r8),allocatable:: DS_lat(:)

    type(iosystem_desc_t),pointer:: pio_subsystem
    integer :: nl

    character(len=32):: HmapOpt

    HmapOpt = trim(I_HmapOpt)
 
    ! Check that DS_lon(DS_nlon) and DS_lat(DS_nlat) are in degrees
    !---------------------------------------------------------------------------


    ! Local aliases for input Grid information
    !-------------------------------------------------
    allocate( M_lon(Hgrid_m%ncol))
    allocate( M_lat(Hgrid_m%ncol))
    allocate(DS_lon(Hgrid_d%nlon))
    allocate(DS_lat(Hgrid_d%nlat))
    M_ncol    = Hgrid_m%ncol
    M_lon(:)  = Hgrid_m%Lon(:)
    M_lat(:)  = Hgrid_m%Lat(:)
    DS_nlon   = Hgrid_d%nlon
    DS_nlat   = Hgrid_d%nlat
    DS_lon(:) = Hgrid_d%Lon(:)
    DS_lat(:) = Hgrid_d%Lat(:)

    ! We require that the longitudes be in the range [0,360]
    !-------------------------------------------------------
    do nc=1,Hgrid_m%ncol
      if(M_lon(nc).lt.  0._r8) M_lon(nc) = M_lon(nc) + 360._r8
      if(M_lon(nc).gt.360._r8) M_lon(nc) = M_lon(nc) - 360._r8
    end do

    ! Initialiaze HmapData values from given grid and interpolation 
    ! information.
    !--------------------------------------------------------------
    HmapData%ncol     = M_ncol
    HmapData%nlon     = DS_nlon
    HmapData%nlat     = DS_nlat
    HmapData%HmapOpt  = trim(HmapOpt)

    if(allocated(HmapData%lon)) deallocate(HmapData%lon)
    if(allocated(HmapData%lat)) deallocate(HmapData%lat)
    allocate(HmapData%lon(HmapData%nlon))
    allocate(HmapData%lat(HmapData%nlat))
    HmapData%lon(:) = DS_lon(:)
    HmapData%lat(:) = DS_lat(:)

    ! Some preliminaries
    !--------------------
    if((trim(HmapData%HmapOpt).eq.'BILINEAR_INTERP'     ).or. &
       (trim(HmapData%HmapOpt).eq.'VECT_BILINEAR_INTERP')     ) then
      Nhalo = 0
    elseif((trim(HmapData%HmapOpt).eq.'BICUBIC_INTERP'     ).or. & 
           (trim(HmapData%HmapOpt).eq.'VECT_BICUBIC_INTERP')     ) then
      Nhalo = 1
    elseif((trim(HmapData%HmapOpt).eq.'BIQUINTIC_INTERP'     ).or. &
           (trim(HmapData%HmapOpt).eq.'VECT_BIQUINTIC_INTERP')     ) then
      Nhalo = 2
    else
      ! ERROR: Unknown Interpolation value
      !---------------------------------------
      if(masterproc) then
        write(iulog, *)'hmap_rect_init: Unknown HmapOpt = ',trim(HmapData%HmapOpt)
      endif
      call endrun('UNKNOWN HMAP MAP OPTION')
    endif

    ! Set info about Input lat/lon grid
    !-------------------------------------
    if(HmapData%lat(1).gt.HmapData%lat(HmapData%nlat)) then
      LatUP  = .false.
      LatMax = HmapData%lat(1)
      LatMin = HmapData%lat(HmapData%nlat)
    else
      LatUp  = .true.
      LatMax = HmapData%lat(HmapData%nlat)
      LatMin = HmapData%lat(1)
    endif

    if((LatMin.gt.-90._r8).and.(LatMax.gt.90._r8)) then
      LatPOLE = .false.
    else
      LatPOLE = .true.
    endif
 
    LonMin = HmapData%lon(1)
    LonMax = HmapData%lon(HmapData%nlon)
    if(LonMin.lt.0._r8) then
      LonROT = .true.
    else
      LonROT = .false.
    endif
 
    ! Generate augmented lat/lon vectors with wrap-around/halo points
    !------------------------------------------------------------------
    allocate(X_lon((1-Nhalo):(HmapData%nlon+1+Nhalo)))
    X_lon(1:HmapData%nlon) = HmapData%lon(1:HmapData%nlon)
    X_lon(HmapData%nlon+1) = HmapData%lon(1) + 360._r8
    do nh=1,Nhalo
      X_lon(              1-nh) = X_lon(HmapData%nlon+1-nh) - 360._r8
      X_lon(HmapData%nlon+1+nh) = X_lon(              1+nh) + 360._r8
    end do

    if(LatUP) then
      if(LatPOLE) then
        allocate(X_lat((1-Nhalo):(HmapData%nlat+Nhalo)))
        X_lat(1:HmapData%nlat) = HmapData%lat(1:HmapData%nlat)
        do nh=1,Nhalo
          X_lat(            1-nh) = -180._r8 - HmapData%lat(            1+nh)
          X_lat(HmapData%nlat+nh) =  180._r8 - HmapData%lat(HmapData%nlat-nh)
        end do
      else ! .not.LatPOLE
        allocate(X_lat((0-Nhalo):(HmapData%nlat+1+Nhalo)))
        X_lat(1:HmapData%nlat) = HmapData%lat(1:HmapData%nlat)
        X_lat(              0) = -180._r8 - HmapData%lat(            1)
        X_lat(HmapData%nlat+1) =  180._r8 - HmapData%lat(HmapData%nlat)
        do nh=1,Nhalo
          X_lat(              0-nh) = -180._r8 - HmapData%lat(            1+nh)
          X_lat(HmapData%nlat+1+nh) =  180._r8 - HmapData%lat(HmapData%nlat-nh)
        end do
      endif
    else ! .not.LatUP
      if(LatPOLE) then
        allocate(X_lat((1-Nhalo):(HmapData%nlat+Nhalo)))
        X_lat(1:HmapData%nlat) = HmapData%lat(1:HmapData%nlat)
        do nh=1,Nhalo
          X_lat(            1-nh) =  180._r8 - HmapData%lat(            1+nh)
          X_lat(HmapData%nlat+nh) = -180._r8 - HmapData%lat(HmapData%nlat-nh)
        end do
      else ! .not.LatPOLE
        allocate(X_lat((0-Nhalo):(HmapData%nlat+1+Nhalo)))
        X_lat(1:HmapData%nlat) = HmapData%lat(1:HmapData%nlat)
        X_lat(              0) =  180._r8 - HmapData%lat(            1)
        X_lat(HmapData%nlat+1) = -180._r8 - HmapData%lat(HmapData%nlat)
        do nh=1,Nhalo
          X_lat(              0-nh) =  180._r8 - HmapData%lat(            1+nh)
          X_lat(HmapData%nlat+1+nh) = -180._r8 - HmapData%lat(HmapData%nlat-nh)
        end do
      endif
    endif
 
    ! Loop thru Grid columns and locate their postion 
    ! in the augmented recilinear input grid.
    !-------------------------------------------------
    allocate(X_lonLo  (HmapData%ncol))
    allocate(X_lonHi  (HmapData%ncol))
    allocate(X_latLo  (HmapData%ncol))
    allocate(X_latHi  (HmapData%ncol))
    allocate(X_lonFrac(HmapData%ncol))
    allocate(X_latFrac(HmapData%ncol))
    allocate(X_lonSize(HmapData%ncol))
    allocate(X_latSize(HmapData%ncol))

    do nc=1,HmapData%ncol
    do ii=1,HmapData%nlon
      if((M_lon(nc).ge.X_lon(ii)).and.(M_lon(nc).le.X_lon(ii+1))) then
        X_lonLo  (nc) = ii
        X_lonHi  (nc) = ii+1
        X_lonSize(nc) = (X_lon(ii+1)-X_lon(ii))
        X_lonFrac(nc) = (M_lon(nc)-X_lon(ii))/(X_lon(ii+1)-X_lon(ii))
        exit
      endif
    end do
    end do

    if(LatUP) then
      if(LatPOLE) then
        do nc=1, HmapData%ncol
        do jj=1,(HmapData%nlat-1)
          if((M_lat(nc).ge.X_lat(jj)).and.(M_lat(nc).le.X_lat(jj+1))) then
            X_latLo  (nc) = jj
            X_latHi  (nc) = jj+1
            X_latSize(nc) = abs(X_lat(jj+1)-X_lat(jj))
            X_latFrac(nc) = (M_lat(nc)-X_lat(jj))/(X_lat(jj+1)-X_lat(jj))
            exit
          endif
        end do
        end do
      else ! .not.LatPOLE
        do nc=1,HmapData%ncol
        do jj=0,HmapData%nlat
          if((M_lat(nc).ge.X_lat(jj)).and.(M_lat(nc).le.X_lat(jj+1))) then
            X_latLo  (nc) = jj
            X_latHi  (nc) = jj+1
            X_latSize(nc) = abs(X_lat(jj+1)-X_lat(jj))
            X_latFrac(nc) = (M_lat(nc)-X_lat(jj))/(X_lat(jj+1)-X_lat(jj))
            exit
          endif
        end do
        end do
      endif
    else ! .not.LatUP
      if(LatPOLE) then
        do nc=1, HmapData%ncol
        do jj=1,(HmapData%nlat-1)
          if((M_lat(nc).le.X_lat(jj)).and.(M_lat(nc).ge.X_lat(jj+1))) then
            X_latLo  (nc) = jj
            X_latHi  (nc) = jj+1
            X_latSize(nc) = abs(X_lat(jj+1)-X_lat(jj))
            X_latFrac(nc) = (M_lat(nc)-X_lat(jj))/(X_lat(jj+1)-X_lat(jj))
            exit
          endif
        end do
        end do
      else ! .not.LatPOLE
        do nc=1,HmapData%ncol
        do jj=0,HmapData%nlat
          if((M_lat(nc).le.X_lat(jj)).and.(M_lat(nc).ge.X_lat(jj+1))) then
            X_latLo  (nc) = jj
            X_latHi  (nc) = jj+1
            X_latSize(nc) = abs(X_lat(jj+1)-X_lat(jj))
            X_latFrac(nc) = (M_lat(nc)-X_lat(jj))/(X_lat(jj+1)-X_lat(jj))
            exit
          endif
        end do
        end do
      endif
    endif
 
    ! Now for each grid column, create a list of neighboring values for 
    ! the specified number of halo points and accumulate a degenerate 
    ! list of augmented input grid indices along the way.
    ! (degenerate = input grid points can be duplicated)
    !-------------------------------------------------------------------
    HmapData%NbrCount = (2+2*Nhalo)**2
    allocate(D_List   (HmapData%NbrCount,HmapData%ncol))
    allocate(D_Vparity(HmapData%NbrCount*HmapData%ncol))
    allocate(D_lonInd (HmapData%NbrCount*HmapData%ncol))
    allocate(D_latInd (HmapData%NbrCount*HmapData%ncol))

    npts=0
    D_Vparity(:) = 1._r8
    do nc=1,HmapData%ncol
      ! set the first 4 bounding points as a 
      ! minimum for Bi-Linear Interpolation
      !---------------------------------------
      nind=4
      npts=npts+1
      D_List  (1,nc) = npts
      D_lonInd(npts) = X_lonLo(nc)
      D_latInd(npts) = X_latLo(nc)
 
      npts=npts+1
      D_List  (2,nc) = npts
      D_lonInd(npts) = X_lonHi(nc)
      D_latInd(npts) = X_latLo(nc)
 
      npts=npts+1
      D_List  (3,nc) = npts
      D_lonInd(npts) = X_lonHi(nc)
      D_latInd(npts) = X_latHi(nc)
 
      npts=npts+1
      D_List  (4,nc) = npts
      D_lonInd(npts) = X_lonLo(nc)
      D_latInd(npts) = X_latHi(nc)
 
      ! Add the specified number of Halo points in each direction with
      ! index numbers spiraling outward in counter-clockwise direction
      !---------------------------------------------------------------
      do nh=1,Nhalo
        lonLo = X_lonLo(nc)-nh
        lonhi = X_lonHi(nc)+nh
        latLo = X_latLo(nc)-nh
        lathi = X_latHi(nc)+nh
 
        do nn=1,(1+2*nh)
          npts=npts+1
          nind=nind+1
          D_List  (nind,nc) = npts
          D_lonInd(npts)    = lonLo + (nn-1)
          D_latInd(npts)    = latLo
        end do
 
        do nn=1,(1+2*nh)
          npts=npts+1
          nind=nind+1
          D_List  (nind,nc) = npts
          D_lonInd(npts)    = lonHi
          D_latInd(npts)    = latLo + (nn-1)
        end do
 
        do nn=1,(1+2*nh)
          npts=npts+1
          nind=nind+1
          D_List  (nind,nc) = npts
          D_lonInd(npts)    = lonHi - (nn-1)
          D_latInd(npts)    = latHi
        end do
 
        do nn=1,(1+2*nh)
          npts=npts+1
          nind=nind+1
          D_List  (nind,nc) = npts
          D_lonInd(npts)    = lonLo
          D_latInd(npts)    = latHi - (nn-1)
        end do
      end do ! nh=1,Nhalo
    end do ! nc=1,HmapData%ncol

    ! Convert lon wrap-aroud indices to the corresponding input grid indices
    !--------------------------------------------------------------------------
    do nn=1,npts
     if(D_lonInd(nn).lt.1) then
       D_lonInd(nn) = D_lonInd(nn) + HmapData%nlon
     endif
     if(D_lonInd(nn).gt.HmapData%nlon) then
       D_lonInd(nn) = D_lonInd(nn) - HmapData%nlon
     endif
    end do

    ! Convert polar wrap-aroud indices to the corresponding input grid indices. 
    ! For polar wrap-around points, the indices are adjusted to the input grid 
    ! indices on the opposite side of the pole. For vector interpolation across 
    ! the pole, the basis directions are flipped so the set the parity to -1.
    !---------------------------------------------------------------------------
    if(LatPOLE) then
      do nn=1,npts
        if(D_latInd(nn).lt.1) then
          D_Vparity(nn) = -1._r8
          D_latInd (nn) = 2 - D_latInd(nn)
          D_lonInd (nn) = D_lonInd(nn) + (HmapData%nlon/2)
          if(D_lonInd(nn).gt.HmapData%nlon) D_lonInd(nn) = D_lonInd(nn) - HmapData%nlon
        endif
        if(D_latInd(nn).gt.HmapData%nlat) then
          D_Vparity(nn) = -1._r8
          D_latInd (nn) = 2*HmapData%nlat - D_latInd(nn)
          D_lonInd (nn) = D_lonInd(nn) + (HmapData%nlon/2)
          if(D_lonInd(nn).gt.HmapData%nlon) D_lonInd(nn) = D_lonInd(nn) - HmapData%nlon
        endif
      end do
    else ! .not.LatPOLE
      do nn=1,npts
        if(D_latInd(nn).lt.1) then
          D_Vparity(nn) = -1._r8
          D_latInd (nn) = 1 - D_latInd(nn)
          D_lonInd (nn) = D_lonInd(nn) + (HmapData%nlon/2)
          if(D_lonInd(nn).gt.HmapData%nlon) D_lonInd(nn) = D_lonInd(nn) - HmapData%nlon
        endif
        if(D_latInd(nn).gt.HmapData%nlat) then
          D_Vparity(nn) = -1._r8
          D_latInd (nn) = 2*HmapData%nlat + 1 - D_latInd(nn)
          D_lonInd (nn) = D_lonInd(nn) + (HmapData%nlon/2)
          if(D_lonInd(nn).gt.HmapData%nlon) D_lonInd(nn) = D_lonInd(nn) - HmapData%nlon
        endif
      end do
    endif

    ! Now go thru and remove the duplicate input grid indices. 
    !---------------------------------------------------------
    allocate(indMap(npts))
    allocate(X_List   (HmapData%NbrCount,HmapData%ncol))
    allocate(X_Vparity(npts))
    allocate(X_lonInd (npts))
    allocate(X_latInd (npts))

    HmapData%npts = 1
    X_Vparity(1) = D_Vparity(1)
    X_lonInd (1) = D_lonInd (1)
    X_latInd (1) = D_latInd (1)
    indMap   (1) = 1
    do nn=2,npts
      ! See if the current point is in the existing list
      !--------------------------------------------------
      newPt = .true.
      do kk=1,HmapData%npts
        if((X_lonInd(kk).eq.D_lonInd(nn)).and. &
           (X_latInd(kk).eq.D_latInd(nn))      ) then
          newPt  = .false.
          kk_old = kk
          exit
        endif
      end do
      ! add a new unique point and update the index map 
      !--------------------------------------------------
      if(newPt) then
        HmapData%npts = HmapData%npts +1
        X_Vparity(HmapData%npts) = D_Vparity(nn)
        X_lonInd (HmapData%npts) = D_lonInd (nn)
        X_latInd (HmapData%npts) = D_latInd (nn)
        indMap  (nn)   = HmapData%npts
      else
        indMap  (nn) = kk_old
      endif
    end do

    ! Update the neighbor list for changes. 
    ! Free up the degenerate lists
    !---------------------------------------
    do nc=1,HmapData%ncol
    do nn=1,HmapData%NbrCount
      X_List(nn,nc) = indMap(D_List(nn,nc))
    end do
    end do

    deallocate(D_List)
    deallocate(D_Vparity)
    deallocate(D_lonInd )
    deallocate(D_latInd )

    ! Create a idof list of 2D array indices used by 
    ! PIO to create the iodesc values.
    !----------------------------------------------------------
    if(allocated(HmapData%dof2D)) deallocate(HmapData%dof2D)
    allocate(HmapData%dof2D(HmapData%npts))
    do nn=1,HmapData%npts
      HmapData%dof2D(nn) = X_lonInd(nn) + HmapData%nlon*(X_latInd(nn)-1)
    end do

    ! Sort the indices into the order that data is layed 
    ! out on the disk for better reads. Ajust neighbor indices 
    ! as needed.          (Does PIO already do this??)
    !---------------------------------------------------------------
    if(allocated(HmapData%NbrList)) deallocate(HmapData%NbrList)
    if(allocated(HmapData%NbrList)) deallocate(HmapData%NbrList)
    if(allocated(HmapData%Vparity)) deallocate(HmapData%Vparity)
    if(allocated(HmapData%lonInd )) deallocate(HmapData%lonInd )
    if(allocated(HmapData%latInd )) deallocate(HmapData%latInd )
    allocate(HmapData%NbrList(HmapData%NbrCount,HmapData%ncol))
    allocate(HmapData%Vparity(HmapData%npts))
    allocate(HmapData%lonInd (HmapData%npts))
    allocate(HmapData%latInd (HmapData%npts))

    call ind_sort(HmapData%dof2D,indMap,HmapData%npts)

    do nn=1,HmapData%npts
      HmapData%Vparity(indMap(nn)) = X_Vparity(nn)
      HmapData%lonInd (indMap(nn)) = X_lonInd (nn)
      HmapData%latInd (indMap(nn)) = X_latInd (nn)
    end do

    do nc=1,HmapData%ncol
    do nn=1,HmapData%NbrCount
      HmapData%NbrList(nn,nc) = indMap(X_List(nn,nc))
    end do
    end do

    deallocate(indMap)
    deallocate(X_List)
    deallocate(X_Vparity)
    deallocate(X_lonInd )
    deallocate(X_latInd )
 
    ! Reload the sorted idof list of 2D array indices used by 
    ! PIO to create the iodesc values.    (NOT USED ANYMORE)
    !---------------------------------------------------------
    do nn=1,HmapData%npts
      HmapData%dof2D(nn) = HmapData%lonInd(nn) + HmapData%nlon*(HmapData%latInd(nn)-1)
    end do
 
    ! The reads will return an array of input values needed for this PE, 
    ! the neighbor list contains the indices of those values that 
    ! contribute to the interpolation at a given grid column.  
    ! Now, we need the interpolation weights that are applied to 
    ! each contributing point as it is summed up.
    !------------------------------------------------------------
    if(allocated(HmapData%NbrIwgt)) deallocate(HmapData%NbrIwgt)
    allocate(HmapData%NbrIwgt(HmapData%NbrCount,HmapData%ncol))
 
    if((trim(HmapData%HmapOpt).eq.'BILINEAR_INTERP'     ).or. &
       (trim(HmapData%HmapOpt).eq.'VECT_BILINEAR_INTERP')     ) then


      ! Initialize Bi-linear interpolation weights  
      !  (1st order polynomials: NO HALO)
      !---------------------------------------------
      do nc=1,HmapData%ncol
        HmapData%NbrIwgt(1,nc) = (1._r8-X_lonFrac(nc))*(1._r8-X_latFrac(nc))
        HmapData%NbrIwgt(2,nc) =        X_lonFrac(nc) *(1._r8-X_latFrac(nc))
        HmapData%NbrIwgt(3,nc) =        X_lonFrac(nc) *       X_latFrac(nc)
        HmapData%NbrIwgt(4,nc) = (1._r8-X_lonFrac(nc))*       X_latFrac(nc)
      end do
    elseif((trim(HmapData%HmapOpt).eq.'BICUBIC_INTERP'     ).or. & 
           (trim(HmapData%HmapOpt).eq.'VECT_BICUBIC_INTERP')     ) then
      ! Initialize Bi-Cubic interpolation weights  
      !  (3rd order polynomials:  HALO=1)
      !---------------------------------------------
      if(masterproc) then
        write(iulog, *)'hmap_rect_init: BICUBIC interpolation not implemented yet.'
      endif
      call endrun('BICUBIC HMAP INTERPOLATION NOT IMPLEMENTED YET')
    elseif((trim(HmapData%HmapOpt).eq.'BIQUINTIC_INTERP'     ).or. &
           (trim(HmapData%HmapOpt).eq.'VECT_BIQUINTIC_INTERP')     ) then
      ! Initialize Bi-Quintic interpolation weights  
      !  (5th order polynomials:  HALO=2)
      !---------------------------------------------
      if(masterproc) then
        write(iulog, *)'hmap_rect_init: BIQUINTIC interpolation not implemented yet.'
      endif
      call endrun('BIQUINTIC HMAP INTERPOLATION NOT IMPLEMENTED YET')
    else
      ! ERROR: Unknown Interpolation value
      !---------------------------------------
      if(masterproc) then
        write(iulog, *)'hmap_rect_init: Unknown HmapOpt = ',trim(HmapData%HmapOpt)
      endif
      call endrun('UNKNOWN HMAP MAP OPTION')
    endif

    ! My Moma always said: Clean up this mess!
    !-------------------------------------------
    deallocate(X_lon)
    deallocate(X_lat)
    deallocate(X_lonLo  )
    deallocate(X_lonHi  )
    deallocate(X_latLo  )
    deallocate(X_latHi  )
    deallocate(X_lonSize)
    deallocate(X_latSize)
    deallocate(X_lonFrac)
    deallocate(X_latFrac)
    deallocate( M_lon)
    deallocate( M_lat)
    deallocate(DS_lon)
    deallocate(DS_lat)
 
    hmap_rect_init = .true.

    ! End Routine
    !-----------
    return
  end function hmap_rect_init
  !================================================================


  !================================================================
  logical function hmap_rect_apply(OutDat,DPdesc,HmapOpt,HmapData)
    ! 
    ! hmap_rect_apply:
    !
    !==============================================================
 use netcdf                    ! PFC DIAG Temp to test PIO ERROR with nf90
    !
    ! Passed Variables 
    !-------------------
    real(r8)              ,intent(out  ):: OutDat(:,:)
    type(PathDesc_t)      ,intent(inout):: DPdesc
    character(len=*)      ,intent(in   ):: HmapOpt
    type(hmap_rect_data_t),intent(inout):: HmapData
    !
    ! Local Values
    !--------------
    integer             :: nl,nc,nn,nind,nind0
    integer             :: ierr
    real(r8),allocatable:: RawData(:)

!================================================================
!PFC DIAG Temporary bypass to test PIO using nf90 call instead
    real(r8),allocatable:: Vdata  (:)
    integer:: ncid,kk,ii,iret,varid
 IF(.FALSE.) THEN
    ncid = DPdesc%ncid_nf90
    allocate(RawData(HmapData%npts*DPdesc%nlev_d))
    iret = nf90_inq_varid(ncid,trim(DPdesc%varname(1:DPdesc%vlen)),varid)
    if(DPdesc%nlev_d.eq.1) then
      ! HORIZ_ONLY
      !--------------
      nn = 0
      do ii=1,HmapData%npts
        nn = nn + 1
        iret = nf90_get_var(ncid,varid,RawData(nn),                              &
               start=(/HmapData%lonInd(ii),HmapData%latInd(ii),DPdesc%timelevel/))
      end do
    else
      ! 3D ARRAY
      !--------------
      allocate(Vdata(DPdesc%nlev_d))
      do ii=1,HmapData%npts
        iret = nf90_get_var(ncid,varid,Vdata,                                        &
               start=(/HmapData%lonInd(ii),HmapData%latInd(ii),1,DPdesc%timelevel/), &
               count=(/1,1,DPdesc%nlev_d,1/)                                         )
        do kk=1,DPdesc%nlev_d
          nn = ii + (kk-1)*HmapData%npts
          RawData(nn) = Vdata(kk)
        end do
      end do
      deallocate(Vdata)
    endif
 ELSE
    ! skip to the current timelevel
    !--------------------------------
    call pio_setframe(DPdesc%ncid,DPdesc%VarDesc,                       &
                              int(DPdesc%timelevel,kind=pio_offset_kind))
    allocate(RawData(HmapData%npts*DPdesc%nlev_d))

    ! Read in the Raw data for the variable
    !------------------------------------------
    call pio_read_darray(DPdesc%ncid,DPdesc%VarDesc,DPdesc%iodesc,RawData,ierr)
 ENDIF
!PFC DIAG Temporary bypass to test PIO using nf90 call instead
!================================================================

    ! Map values to output grid
    !----------------------------
    OutDat(1:HmapData%ncol,1:DPdesc%nlev_d) = 0._r8
    if(HmapOpt(1:4).eq.'VECT') then
      ! VECTOR MAP
      !--------------
      do nl=1,DPdesc%nlev_d
      do nc=1,HmapData%ncol
      do nn=1,HmapData%NbrCount
        nind0 = HmapData%NbrList(nn,nc)
        nind  = nind0 + HmapData%npts*(nl-1)
        ! If *ANY* missing values are found in the halo, set the 
        ! interp point to the fill_value and move on. 
        !  DEVNOTE: This should be cleaned up to salvage avalaible 
        !           data for interp value.
        !----------------------------------------------------------
        if(RawData(nind).eq.DPdesc%fillvalue) then
          OutDat(nc,nl) = DPdesc%fillvalue
          exit
        else
          OutDat(nc,nl) = OutDat(nc,nl) + HmapData%NbrIwgt(nn,nc)*RawData(nind) &
                                         *HmapData%Vparity(nind0)
        endif
      end do
      end do
      end do
    else
      ! SCALAR MAP
      !--------------
      do nl=1,DPdesc%nlev_d
      do nc=1,HmapData%ncol
      do nn=1,HmapData%NbrCount
        nind0 = HmapData%NbrList(nn,nc)
        nind  = nind0 + HmapData%npts*(nl-1)
        ! If *ANY* missing values are found in the halo, set the 
        ! interp point to the fill_value and move on. 
        !  DEVNOTE: This should be cleaned up to salvage avalaible 
        !           data for interp value.
        !----------------------------------------------------------
        if(RawData(nind).eq.DPdesc%fillvalue) then
          OutDat(nc,nl) = DPdesc%fillvalue
          exit
        else
          OutDat(nc,nl) = OutDat(nc,nl) + HmapData%NbrIwgt(nn,nc)*RawData(nind)
        endif
      end do
      end do
      end do
    endif

    deallocate(RawData)
    hmap_rect_apply = .true.

    ! End Routine
    !-----------
    return
  end function hmap_rect_apply
  !================================================================


  !================================================================
  logical function hmap_native_chk(Hgrid_m,Hgrid_d,I_HmapOpt,DataMap)
    !
    ! hmap_native_chk: Check the input Hgrid_d and output Hgrid_m grids
    !                  check them against the given HmapOpt for 
    !                  datasets with rectilinear grids.
    !==============================================================
    !
    ! Passed variables  
    !-----------------
    type(packed_hgrid_t) ,intent(in):: Hgrid_m
    type(dataset_hgrid_t),intent(in):: Hgrid_d
    character(len=*)     ,intent(in):: I_HmapOpt
    type(hmaps_t)        ,intent(in):: DataMap
    !
    ! Local Values
    !---------------
    character(len=32):: HmapOpt

    HmapOpt = trim(I_HmapOpt)

    ! First check that the in/out grid type match the map
    !------------------------------------------------------


    ! Now check everything else
    !--------------------------

 
    hmap_native_chk = .true.

    ! End Routine
    !-----------
    return
  end function hmap_native_chk
  !================================================================


  !================================================================
  logical function hmap_native_init(Hgrid_m,Hgrid_d,I_HmapOpt,HmapData)
    ! 
    ! hmap_native_init: 
    !
    !==============================================================
    !
    ! Passed Variables
    !------------------
    type(packed_hgrid_t)  ,intent(in   ):: Hgrid_m
    type(dataset_hgrid_t) ,intent(in   ):: Hgrid_d
    character(len=*)      ,intent(in   ):: I_HmapOpt
    type(hmap_rect_data_t),intent(inout):: HmapData
    !
    ! Local Values
    !--------------



    hmap_native_init = .true.

    ! End Routine
    !-----------
    return
  end function hmap_native_init
  !================================================================


  !================================================================
  logical function hmap_native_apply(OutDat,DPdesc,HmapOpt,HmapData)
    ! 
    ! hmap_native_apply:
    !
    !==============================================================
    !
    ! Passed Variables 
    !-------------------
    real(r8)              ,intent(out  ):: OutDat(:,:)
    type(PathDesc_t)      ,intent(inout):: DPdesc
    character(len=*)      ,intent(in   ):: HmapOpt
    type(hmap_rect_data_t),intent(inout):: HmapData
    !
    ! Local Values
    !--------------
    integer             :: nl,nc,nn,nind,nind0
    integer             :: ierr
!    real(r8),allocatable:: RawData(:)
    real(r8),allocatable:: RawData(:,:,:)


    ! Would rather have the native PIO calls here rather than in 
    ! Dataconnector_mod, but the implementation would not be clean 
    ! with the additional arguments that would need to be passed.
    !
    ! FIX this up when all of the other dragons have been slain...
    !-------------------------------------------------------------
    call endrun('hmap_native_apply: Currently, This function should never be called!')



    ! skip to the current timelevel
    !--------------------------------
    call pio_setframe(DPdesc%ncid,DPdesc%VarDesc,                       &
                              int(DPdesc%timelevel,kind=pio_offset_kind))
!    allocate(RawData(HmapData%npts*DPdesc%nlev_d))
!    allocate(RawData(DPdesc%dim1b:DPdesc%dim1e,DPdesc%dim2b:DPdesc%dim2e,DPdesc%dim3b:DPdesc%dim3e))

    ! Read in the Raw data for the variable
    !------------------------------------------
    call pio_read_darray(DPdesc%ncid,DPdesc%VarDesc,DPdesc%iodesc,RawData,ierr)

    ! Map values to output grid
    !----------------------------
    OutDat(1:HmapData%ncol,1:DPdesc%nlev_d) = 0._r8

    deallocate(RawData)
    hmap_native_apply = .true.

    ! End Routine
    !-----------
    return
  end function hmap_native_apply
  !================================================================


  !================================================================
  logical function hmap_USER_chk(Hgrid_m,Hgrid_d,I_HmapOpt,DataMap)
    !
    ! hmap_USER_chk: Check the input Hgrid_d and output Hgrid_m grids
    !                check them against the given HmapOpt and
    !                its requirements.
    !
    !       This is a TEMPLATE/STUB to help users wishing to implement 
    !       their own mapping method.
    !
    !==============================================================
    !
    ! Passed variables  
    !-----------------
    type(packed_hgrid_t) ,intent(in):: Hgrid_m
    type(dataset_hgrid_t),intent(in):: Hgrid_d
    character(len=*)     ,intent(in):: I_HmapOpt
    type(hmaps_t)        ,intent(in):: DataMap
    !
    ! Local Values
    !--------------
    character(len=32):: HmapOpt

    HmapOpt = trim(I_HmapOpt)

    ! First check that the in/out grid type match the map
    !------------------------------------------------------
    hmap_USER_chk = .true.

    if((trim(Hgrid_d%Type).ne.trim(DataMap%InType )).or. &
       (trim(Hgrid_m%Type).ne.trim(DataMap%OutType))     ) then
      hmap_USER_chk = .false.
      return
    endif

    ! Now check everything else
    !--------------------------

 
    ! End Routine
    !-----------
    return
  end function hmap_USER_chk
  !================================================================


  !================================================================
  logical function hmap_USER_init(Hgrid_m,Hgrid_d,I_HmapOpt,HmapData)
    ! 
    ! hmap_USER_init: Initialize the HMAP values needed to implement horizontal 
    !                 interpolation of input data for the gridpoints contained 
    !                 in this PE.
    !
    !       This is a TEMPLATE/STUB to help users wishing to implement 
    !       their own mapping method.
    !
    !       Initialize the (Hmap%) datastructure values needed for the specified
    !       horizontal interpolation method.
    !
    !==============================================================
    !
    ! Passed Variables
    !------------------
    type(packed_hgrid_t)  ,intent(in   ):: Hgrid_m
    type(dataset_hgrid_t) ,intent(in   ):: Hgrid_d
    character(len=*)      ,intent(in   ):: I_HmapOpt
    type(hmap_rect_data_t),intent(inout):: HmapData
    !
    ! Local Values
    !-------------
    character(len=32):: HmapOpt

    HmapOpt = trim(I_HmapOpt)

    hmap_USER_init = .true.

    ! End Routine
    !-----------
    return
  end function hmap_USER_init
  !================================================================


  !================================================================
  logical function hmap_USER_apply(OutDat,DPdesc,HmapOpt,HmapData)
    ! 
    ! hmap_USER_apply: This is an USER based map in which values are 
    !                  interpolated using a user defined mapping.
    !
    !       This is a TEMPLATE/STUB to help users wishing to implement 
    !       their own mapping method.
    !
    !   Constraints: 
    !
    !==============================================================
    !
    ! Passed Variables 
    !-------------------
    real(r8)              ,intent(out  ):: OutDat(:,:)
    type(PathDesc_t)      ,intent(inout):: DPdesc
    character(len=*)      ,intent(in   ):: HmapOpt
    type(hmap_rect_data_t),intent(inout):: HmapData
    !
    ! Local Values
    !--------------


    OutDat(:,:) = 0._r8

    hmap_USER_apply = .true.

    ! End Routine
    !-----------
    return
  end function hmap_USER_apply
  !================================================================

end module Hmaps_mod
