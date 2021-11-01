module Dataconnector_mod
!===========================================================================
!
! MODULE: Dataconnector_mod
!
! DESCRIPTION: Manage the mapping of data values from a given dataset to 
!              model data arrays.
!
!===========================================================================
  ! Useful Modules
  !---------------------
  use pio,           only: pio_offset_kind, pio_double, pio_max_var_dims,       &
                           file_desc_t, var_desc_t, io_desc_t, iosystem_desc_t, &
                           pio_inq_dimid, pio_inq_varid, pio_initdecomp,        &
                           pio_setframe, pio_read_darray,PIO_MAX_NAME
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use spmd_utils,    only: masterproc
  use cam_abortutils,only: endrun
  use cam_logfile,   only: iulog
  use shr_sys_mod,   only: shr_sys_flush
  use Dataset_mod,   only: dataset_t, dataset_var_t, inq_dataset, inq_dataset_var, &
                           inq_dataset_status, inq_dataset_varid, inq_dataset_pio, &
                           inq_dataset_ncid
  use Hmaps_mod,     only: packed_hgrid_t, dataset_hgrid_t, hmap_rect_data_t,    & 
                           known_hgrids_t, known_hmaps_t, init_Hmaps_mod,        &
                           print_known_Hmaps, inq_Hmaps, copy_Hgrid,             &
                           deallocate_Hgrid, chk_Hmap, init_Hmap, apply_Hmap
  use Vmaps_mod,     only: vgrid_t, vmap_t, known_vgrids_t, known_vmaps_t, vmap_data_t,  &
                           init_Vmaps_mod, print_known_Vmaps, inq_Vmaps,         &
                           copy_Vgrid, deallocate_Vgrid, chk_Vmap, init_Vmap,    &
                           apply_Vmap
  use PathDesc_type ,only: PathDesc_t
  use Dataset_mod,   only: inq_dataset_ncid_nf90             !PFC DIAG Temp to test PIO ERROR with nf90

  ! Set all Global values and routines to private by default 
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private
  save

  public :: init_dataconnector_mod
  public :: register_dataPath
  public :: print_dataset_paths
  public :: apply_dataPath
  private:: apply_dataPath_2D
  private:: apply_dataPath_3D
  public :: apply_horiz_dataPath
  private:: apply_horiz_dataPath_2D
  private:: apply_horiz_dataPath_3D
  public :: apply_vert_dataPath
  public :: unpack_dataPath_result
  private:: unpack_dataPath_result_2D
  private:: unpack_dataPath_result_3D
  private:: pack_model_Hgrid_2D
  private:: pack_model_Hgrid_3D
  private:: unpack_model_Hgrid_2D
  private:: unpack_model_Hgrid_3D
  public :: inq_dataPath_dims
!  public :: register_gridPath  !  For model Grid -> Grid maps
!                                  ** Need to add a ListGridPaths_t ---etc...
!  public :: apply_gridPath       !  For model Grid -> Grid maps

  ! Interfaces
  !-----------------
  interface apply_dataPath
    module procedure apply_dataPath_2D
    module procedure apply_dataPath_3D
  end interface 

  interface apply_horiz_dataPath
    module procedure apply_horiz_dataPath_2D
    module procedure apply_horiz_dataPath_3D
  end interface 

  interface unpack_dataPath_result
    module procedure unpack_dataPath_result_2D
    module procedure unpack_dataPath_result_3D
  end interface 

  ! Type Definitions
  !-------------------
  type ListDataPaths_t
    ! Manage List of currently Registered datasets
    !----------------------------------------------
      private
      class(LinkDataPaths_t),pointer:: first => null()
      class(LinkDataPaths_t),pointer:: last  => null()
      class(LinkDataPaths_t),pointer:: curr  => null()
    contains
      procedure,pass:: reset     => ListDataPaths_reset
      procedure,pass:: iterate   => ListDataPaths_iterate
      procedure,pass:: isFirst   => ListDataPaths_isFirst
      procedure,pass:: isLast    => ListDataPaths_isLast
      procedure,pass:: get_ID    => ListDataPaths_getLinkID
      procedure,pass:: set_ID    => ListDataPaths_setLinkID
      procedure,pass:: printList => ListDataPaths_print
!      procedure:: rm_DataSet    ** remove all values associated with given DataSetID
  end type ListDataPaths_t

  type GridID_t
    ! Look up table for registered H/V grids associated with each DataSet grid
    !---------------------------------------------------------------------------
    integer            :: idxHgrid_d
    integer            :: idxHORIZ_ONLY
    integer            :: nvgrids
    integer,allocatable:: idxVgrid_d(:)
  end type GridID_t

  type LinkDataPaths_t
    ! List-Links for currently Registered datasets
    !----------------------------------------------
      private
      integer                       :: DataSetID
      integer                       :: ngrids
      type (GridID_t)   ,allocatable:: GridID(:)
      class(ListHgrid_m_t)  ,pointer:: ListHgrid_m  => null()
      class(ListVgrid_m_t)  ,pointer:: ListVgrid_m  => null()
      class(ListHgrid_d_t)  ,pointer:: ListHgrid_d  => null()
      class(ListVgrid_d_t)  ,pointer:: ListVgrid_d  => null()
      class(ListHmaps_t)    ,pointer:: ListHmaps    => null()
      class(ListVmaps_t)    ,pointer:: ListVmaps    => null()
      class(ListPaths_t)    ,pointer:: ListPaths    => null()
      type (LinkDataPaths_t),pointer:: next         => null()
      type (LinkDataPaths_t),pointer:: prev         => null()
    contains
      procedure,pass:: get_ID    => LinkDataPaths_getID
      procedure,pass:: get_Next  => LinkDataPaths_getNext
      procedure,pass:: set_Next  => LinkDataPaths_setNext
      procedure,pass:: get_Prev  => LinkDataPaths_getPrev
      procedure,pass:: set_Prev  => LinkDataPaths_setPrev
      procedure,pass:: printLink => LinkDataPaths_print
  end type LinkDataPaths_t
  !=============================================================================

  !=============================================================================
  type ListHgrid_m_t
    ! List of currently Registered Model Hgrids
    !-------------------------------------------------
      private
      integer                    :: Counter = 0
      type(LinkHgrid_m_t),pointer:: first   => null()
      type(LinkHgrid_m_t),pointer:: last    => null()
      type(LinkHgrid_m_t),pointer:: curr    => null()
    contains
      procedure,pass:: reset     => ListHgrid_m_reset
      procedure,pass:: iterate   => ListHgrid_m_iterate
      procedure,pass:: isFirst   => ListHgrid_m_isFirst
      procedure,pass:: isLast    => ListHgrid_m_isLast
      procedure,pass:: addLink   => ListHgrid_m_addLink
      procedure,pass:: rmLink    => ListHgrid_m_rmLink
      procedure,pass:: get_ID    => ListHgrid_m_getLinkID
      procedure,pass:: get_Data  => ListHgrid_m_getLinkData
      procedure,pass:: set_ID    => ListHgrid_m_setLinkID
      procedure,pass:: printList => ListHgrid_m_print
  end type ListHgrid_m_t

  type LinkHgrid_m_t
    ! List-Links for currently Registered Model Hgrids
    !-------------------------------------------------
      private
      integer                     :: idxHgrid_m
      type(packed_hgrid_t),pointer:: Hgrid_m      => null()
      type(LinkHgrid_m_t) ,pointer:: next         => null()
      type(LinkHgrid_m_t) ,pointer:: prev         => null()
    contains
      procedure,pass:: get_ID    => LinkHgrid_m_getID
      procedure,pass:: get_Data  => LinkHgrid_m_getData
      procedure,pass:: get_Next  => LinkHgrid_m_getNext
      procedure,pass:: set_Next  => LinkHgrid_m_setNext
      procedure,pass:: get_Prev  => LinkHgrid_m_getPrev
      procedure,pass:: set_Prev  => LinkHgrid_m_setPrev
      procedure,pass:: isEqual   => LinkHgrid_m_isEqual 
      procedure,pass:: printLink => LinkHgrid_m_print
  end type LinkHgrid_m_t
  !=============================================================================

  !=============================================================================
  type ListVgrid_m_t
    ! List of currently Registered Model Vgrids
    !-------------------------------------------------
      private
      integer                    :: Counter = 0
      type(LinkVgrid_m_t),pointer:: first   => null()
      type(LinkVgrid_m_t),pointer:: last    => null()
      type(LinkVgrid_m_t),pointer:: curr    => null()
    contains
      procedure,pass:: reset     => ListVgrid_m_reset
      procedure,pass:: iterate   => ListVgrid_m_iterate
      procedure,pass:: isFirst   => ListVgrid_m_isFirst
      procedure,pass:: isLast    => ListVgrid_m_isLast
      procedure,pass:: addLink   => ListVgrid_m_addLink
      procedure,pass:: rmLink    => ListVgrid_m_rmLink
      procedure,pass:: get_ID    => ListVgrid_m_getLinkID
      procedure,pass:: get_Data  => ListVgrid_m_getLinkData
      procedure,pass:: set_ID    => ListVgrid_m_setLinkID
      procedure,pass:: printList => ListVgrid_m_print
  end type ListVgrid_m_t

  type LinkVgrid_m_t
    ! List-Links for currently Registered Model Vgrids
    !-------------------------------------------------
      private
      integer                     :: idxVgrid_m
      type(vgrid_t)       ,pointer:: Vgrid_m      => null()
      type(LinkVgrid_m_t) ,pointer:: next         => null()
      type(LinkVgrid_m_t) ,pointer:: prev         => null()
    contains
      procedure,pass:: get_ID    => LinkVgrid_m_getID
      procedure,pass:: get_Data  => LinkVgrid_m_getData
      procedure,pass:: get_Next  => LinkVgrid_m_getNext
      procedure,pass:: set_Next  => LinkVgrid_m_setNext
      procedure,pass:: get_Prev  => LinkVgrid_m_getPrev
      procedure,pass:: set_Prev  => LinkVgrid_m_setPrev
      procedure,pass:: isEqual   => LinkVgrid_m_isEqual 
      procedure,pass:: printLink => LinkVgrid_m_print
  end type LinkVgrid_m_t
  !=============================================================================

  !=============================================================================
  type ListHgrid_d_t
    ! List of currently Registered DataSet Hgrids
    !-------------------------------------------------
      private
      integer                    :: Counter = 0
      type(LinkHgrid_d_t),pointer:: first   => null()
      type(LinkHgrid_d_t),pointer:: last    => null()
      type(LinkHgrid_d_t),pointer:: curr    => null()
    contains
      procedure,pass:: reset     => ListHgrid_d_reset
      procedure,pass:: iterate   => ListHgrid_d_iterate
      procedure,pass:: isFirst   => ListHgrid_d_isFirst
      procedure,pass:: isLast    => ListHgrid_d_isLast
      procedure,pass:: addLink   => ListHgrid_d_addLink
      procedure,pass:: rmLink    => ListHgrid_d_rmLink
      procedure,pass:: get_ID    => ListHgrid_d_getLinkID
      procedure,pass:: get_Data  => ListHgrid_d_getLinkData
      procedure,pass:: set_ID    => ListHgrid_d_setLinkID
      procedure,pass:: printList => ListHgrid_d_print
  end type ListHgrid_d_t

  type LinkHgrid_d_t
    ! List-Links for currently Registered DataSet Hgrids
    !-------------------------------------------------
      private
      integer                      :: idxHgrid_d
      type(dataset_hgrid_t),pointer:: Hgrid_d      => null()
      type(LinkHgrid_d_t)  ,pointer:: next         => null()
      type(LinkHgrid_d_t)  ,pointer:: prev         => null()
    contains
      procedure,pass:: get_ID    => LinkHgrid_d_getID
      procedure,pass:: get_Data  => LinkHgrid_d_getData
      procedure,pass:: get_Next  => LinkHgrid_d_getNext
      procedure,pass:: set_Next  => LinkHgrid_d_setNext
      procedure,pass:: get_Prev  => LinkHgrid_d_getPrev
      procedure,pass:: set_Prev  => LinkHgrid_d_setPrev
      procedure,pass:: isEqual   => LinkHgrid_d_isEqual 
      procedure,pass:: printLink => LinkHgrid_d_print
  end type LinkHgrid_d_t
  !=============================================================================

  !=============================================================================
  type ListVgrid_d_t
    ! List of currently Registered DataSet Vgrids
    !-------------------------------------------------
      private
      integer                    :: Counter = 0
      type(LinkVgrid_d_t),pointer:: first   => null()
      type(LinkVgrid_d_t),pointer:: last    => null()
      type(LinkVgrid_d_t),pointer:: curr    => null()
    contains
      procedure,pass:: reset     => ListVgrid_d_reset
      procedure,pass:: iterate   => ListVgrid_d_iterate
      procedure,pass:: isFirst   => ListVgrid_d_isFirst
      procedure,pass:: isLast    => ListVgrid_d_isLast
      procedure,pass:: addLink   => ListVgrid_d_addLink
      procedure,pass:: rmLink    => ListVgrid_d_rmLink
      procedure,pass:: get_ID    => ListVgrid_d_getLinkID
      procedure,pass:: get_Data  => ListVgrid_d_getLinkData
      procedure,pass:: set_ID    => ListVgrid_d_setLinkID
      procedure,pass:: printList => ListVgrid_d_print
  end type ListVgrid_d_t

  type LinkVgrid_d_t
    ! List-Links for currently Registered DataSet Vgrids
    !---------------------------------------------------
      private
      integer                     :: idxVgrid_d
      type(vgrid_t)       ,pointer:: Vgrid_d      => null()
      type(LinkVgrid_d_t) ,pointer:: next         => null()
      type(LinkVgrid_d_t) ,pointer:: prev         => null()
    contains
      procedure,pass:: get_ID    => LinkVgrid_d_getID
      procedure,pass:: get_Data  => LinkVgrid_d_getData
      procedure,pass:: get_Next  => LinkVgrid_d_getNext
      procedure,pass:: set_Next  => LinkVgrid_d_setNext
      procedure,pass:: get_Prev  => LinkVgrid_d_getPrev
      procedure,pass:: set_Prev  => LinkVgrid_d_setPrev
      procedure,pass:: isEqual   => LinkVgrid_d_isEqual 
      procedure,pass:: printLink => LinkVgrid_d_print
  end type LinkVgrid_d_t
  !=============================================================================

  !=============================================================================
  type ListHmaps_t
    ! List of currently Registered Hmaps
    !-------------------------------------------------
      private
      integer                  :: Counter = 0
      type(LinkHmaps_t),pointer:: first   => null()
      type(LinkHmaps_t),pointer:: last    => null()
      type(LinkHmaps_t),pointer:: curr    => null()
    contains
      procedure,pass:: reset     => ListHmaps_reset
      procedure,pass:: iterate   => ListHmaps_iterate
      procedure,pass:: isFirst   => ListHmaps_isFirst
      procedure,pass:: isLast    => ListHmaps_isLast
      procedure,pass:: addLink   => ListHmaps_addLink
      procedure,pass:: rmLink    => ListHmaps_rmLink
      procedure,pass:: get_ID    => ListHmaps_getLinkID
      procedure,pass:: get_Data  => ListHmaps_getLinkData
      procedure,pass:: set_ID    => ListHmaps_setLinkID
      procedure,pass:: printList => ListHmaps_print
  end type ListHmaps_t

  type HdataMap_t
    integer                      :: idxHgrid_d
    integer                      :: idxHgrid_m
    type(dataset_hgrid_t),pointer:: Hgrid_in    => null()
    type(packed_hgrid_t ),pointer:: Hgrid_out   => null()
    character(len=32)            :: HmapOpt
    character(len=32)            :: HmapFunction
    type(hmap_rect_data_t)       :: HmapData
  end type HdataMap_t

  type LinkHmaps_t
    ! List-Links for currently Registered Hmaps
    !-------------------------------------------------
      private
      integer                  :: idxHmap
      type(HdataMap_t) ,pointer:: Hmap
      type(LinkHmaps_t),pointer:: next         => null()
      type(LinkHmaps_t),pointer:: prev         => null()
    contains
      procedure,pass:: get_ID    => LinkHmaps_getID
      procedure,pass:: get_Data  => LinkHmaps_getData
      procedure,pass:: get_Next  => LinkHmaps_getNext
      procedure,pass:: set_Next  => LinkHmaps_setNext
      procedure,pass:: get_Prev  => LinkHmaps_getPrev
      procedure,pass:: set_Prev  => LinkHmaps_setPrev
      procedure,pass:: isEqual   => LinkHmaps_isEqual 
      procedure,pass:: printLink => LinkHmaps_print
  end type LinkHmaps_t
  !=============================================================================

  !=============================================================================
  type ListVmaps_t
    ! List of currently Registered Vmaps
    !-------------------------------------------------
      private
      integer                  :: Counter = 0
      type(LinkVmaps_t),pointer:: first   => null()
      type(LinkVmaps_t),pointer:: last    => null()
      type(LinkVmaps_t),pointer:: curr    => null()
    contains
      procedure,pass:: reset     => ListVmaps_reset
      procedure,pass:: iterate   => ListVmaps_iterate
      procedure,pass:: isFirst   => ListVmaps_isFirst
      procedure,pass:: isLast    => ListVmaps_isLast
      procedure,pass:: addLink   => ListVmaps_addLink
      procedure,pass:: rmLink    => ListVmaps_rmLink
      procedure,pass:: get_ID    => ListVmaps_getLinkID
      procedure,pass:: get_Data  => ListVmaps_getLinkData
      procedure,pass:: set_ID    => ListVmaps_setLinkID
      procedure,pass:: printList => ListVmaps_print
  end type ListVmaps_t

!XXX  type Vmap_t
!XXX    integer              :: idxVgrid_d
!XXX    integer              :: idxVgrid_m
!XXX    type(vgrid_t),pointer:: Vgrid_in    => null()
!XXX    type(vgrid_t),pointer:: Vgrid_out   => null()
!XXX    character(len=32)    :: VmapOpt
!XXX    character(len=32)    :: VmapFunction
!XXX    type(vmap_data_t)    :: VmapData
!XXX  end type Vmap_t

  type LinkVmaps_t
    ! List-Links for currently Registered Vmaps
    !-------------------------------------------------
      private
      integer                  :: idxVmap
      type(vmap_t)     ,pointer:: Vmap
      type(LinkVmaps_t),pointer:: next         => null()
      type(LinkVmaps_t),pointer:: prev         => null()
    contains
      procedure,pass:: get_ID    => LinkVmaps_getID
      procedure,pass:: get_Data  => LinkVmaps_getData
      procedure,pass:: get_Next  => LinkVmaps_getNext
      procedure,pass:: set_Next  => LinkVmaps_setNext
      procedure,pass:: get_Prev  => LinkVmaps_getPrev
      procedure,pass:: set_Prev  => LinkVmaps_setPrev
      procedure,pass:: isEqual   => LinkVmaps_isEqual 
      procedure,pass:: printLink => LinkVmaps_print
  end type LinkVmaps_t
  !=============================================================================

  !=============================================================================
  type ListPaths_t
    ! List of currently Registered Paths
    !-------------------------------------------------
      private
      integer                  :: Counter = 0
      type(LinkPaths_t),pointer:: first   => null()
      type(LinkPaths_t),pointer:: last    => null()
      type(LinkPaths_t),pointer:: curr    => null()
    contains
      procedure,pass:: reset     => ListPaths_reset
      procedure,pass:: iterate   => ListPaths_iterate
      procedure,pass:: isFirst   => ListPaths_isFirst
      procedure,pass:: isLast    => ListPaths_isLast
      procedure,pass:: addLink   => ListPaths_addLink
      procedure,pass:: rmLink    => ListPaths_rmLink
      procedure,pass:: get_ID    => ListPaths_getLinkID
      procedure,pass:: get_Data  => ListPaths_getLinkData
      procedure,pass:: set_ID    => ListPaths_setLinkID
      procedure,pass:: printList => ListPaths_print
  end type ListPaths_t

  type ioPTR_t
    type(io_desc_t),pointer:: iodesc
  end type ioPTR_t
  type Path_t
    integer                  :: ngrid
    integer      ,allocatable:: idxHmap(:)
    integer      ,allocatable:: idxVmap(:)
    type(ioPTR_t),allocatable:: ioPTR  (:)
  end type Path_t

  type LinkPaths_t
    ! List-Links for currently Registered Paths
    !-------------------------------------------------
      private
      integer                  :: idxPath
      type(Path_t)     ,pointer:: Path
      type(LinkPaths_t),pointer:: next         => null()
      type(LinkPaths_t),pointer:: prev         => null()
    contains
      procedure,pass:: get_ID    => LinkPaths_getID
      procedure,pass:: get_Data  => LinkPaths_getData
      procedure,pass:: get_Next  => LinkPaths_getNext
      procedure,pass:: set_Next  => LinkPaths_setNext
      procedure,pass:: get_Prev  => LinkPaths_getPrev
      procedure,pass:: set_Prev  => LinkPaths_setPrev
      procedure,pass:: isEqual   => LinkPaths_isEqual 
      procedure,pass:: printLink => LinkPaths_print
  end type LinkPaths_t
  !=============================================================================

  ! Global Parameters
  !-------------------
  real(r8),parameter:: DIFF_TOLERANCE = 1.0d-6
  real(r8),parameter:: PS0_REF        = 100000._r8

  ! Global Values
  !-----------------
  type(known_hgrids_t ):: known_Hgrids
  type(known_hmaps_t  ):: known_Hmaps
  type(known_vgrids_t ):: known_Vgrids
  type(known_vmaps_t  ):: known_Vmaps
  type(ListDataPaths_t):: My_DataPaths

contains

  !================================================================
  subroutine init_dataconnector_mod
    ! 
    ! init_dataconnector_mod: Get Horizontal and Vertical map information 
    !                         from Hmaps_mod/Vmaps_mod and initialize the 
    !                         tables of known grids and mapping options.
    !==============================================================


    ! Get horizontal grid/map metadata 
    !-------------------------------------
    call init_Hmaps_mod
    call inq_Hmaps(known_Hgrids,known_Hmaps)
    call print_known_Hmaps

    ! Get vertical grid/map metadata 
    !-------------------------------------
    call init_Vmaps_mod
    call inq_Vmaps(known_Vgrids,known_Vmaps)
    call print_known_Vmaps

    ! Initialize a table of known grid connections.
    !-----------------------------------------------

    ! Any other init goes here?
    !---------------------------

    ! End Routine
    !------------------
    return
  end subroutine init_dataconnector_mod
  !================================================================


  !================================================================
  function register_dataPath(DataSetID,Hgrid_m,I_HmapOpt,Vgrid_m,I_VmapOpt)
    ! 
    ! register_dataPath: 
    !==============================================================
    use shr_pio_mod,     only: shr_pio_getiosys
    use cam_instance,    only: atm_id
    use cam_grid_support,only: cam_grid_get_decomp
    !
    ! Passed variables
    !--------------------
    integer             :: DataSetID
    type(packed_hgrid_t):: Hgrid_m
    type(vgrid_t)       :: Vgrid_m
    character(len=*)    :: I_HmapOpt
    character(len=*)    :: I_VmapOpt
    integer             :: register_dataPath
    !
    ! Local values
    !----------------
    type(LinkDataPaths_t),pointer:: DPlink
    type(packed_hgrid_t) ,pointer:: My_Hgrid_m
    type(vgrid_t)        ,pointer:: My_Vgrid_m
    type(dataset_hgrid_t),pointer:: Hgrid_d
    type(vgrid_t)        ,pointer:: Vgrid_d
    type(HdataMap_t)     ,pointer:: new_HdataMap
    type(HdataMap_t)     ,pointer:: Hmap
    type(vmap_t)         ,pointer:: new_Vmap
    type(vmap_t)         ,pointer:: Vmap
    type(Path_t)         ,pointer:: new_Path
    integer                      :: idxHgrid_m
    integer                      :: idxVgrid_m
    integer                      :: idxHgrid_d
    integer                      :: idxVgrid_d
    integer                      :: idxHmap
    integer                      :: idxVmap
    logical                      :: newLink
    logical                      :: Hgrid_match,Vgrid_match
    logical                      :: result
    integer                      :: nn,ii,jj,mm,kk,ll
    character(len=32)            :: HmapOpt
    character(len=32)            :: VmapOpt
    integer                      :: Mlen
    type(iosystem_desc_t),pointer:: pio_subsystem
    integer                      :: pio_iotype
  
    integer                    :: ndof
    integer,allocatable        :: idof(:)
    integer                    :: ndim
    integer                    :: Model_dimSize(3)
    character(len=PIO_MAX_NAME):: Model_dimName(3)
    integer                    ::  Data_dimSize(3)
    character(len=PIO_MAX_NAME)::  Data_dimName(3)
    logical                    :: dim1_found, dim2_found, grids_match

    ! Need to accomotate chararacter(len=*) user interface
    !------------------------------------------------------
    HmapOpt = trim(I_HmapOpt)
    VmapOpt = trim(I_VmapOpt)

    ! Get the DataPath link for the given DataSetID, create and 
    ! initialize  a new link if one is not found for DataSetID
    !----------------------------------------------------------
    DPlink => My_DataPaths%set_ID(DataSetID,newLink)

    ! Get needed PIO info from dataset
    !------------------------------------
    pio_subsystem => shr_pio_getiosys(atm_id)
!???    result = inq_dataset_pio(DataSetID,pio_subsystem,pio_iotype)
!    if(.not.result) then
!      write(iulog,*) ' DataSetID=',DataSetID
!      call endrun('register_dataPath: ERROR get PIO info from DataSetID')
!    endif

    ! The lists need to maintain their own copy of the 
    ! datastructures so first create a local copy for 
    ! each model grid
    !--------------------------------------------------
    allocate(My_Hgrid_m)
    allocate(My_Vgrid_m)
    call copy_Hgrid(My_Hgrid_m,Hgrid_m)
    call copy_Vgrid(My_Vgrid_m,Vgrid_m)

    ! Get Unique Hgrid_m model grid indices
    !---------------------------------------
    idxHgrid_m = DPlink%ListHgrid_m%set_ID(My_Hgrid_m,newLink)
    if(.not.newLink) then
      ! The Grid already exists in the List, so deallocate the copy
      !------------------------------------------------------------
      call deallocate_Hgrid(My_Hgrid_m)
      deallocate(My_Hgrid_m)
    endif

    ! Get Unique Vgrid_m model grid indices
    !---------------------------------------
    idxVgrid_m = DPlink%ListVgrid_m%set_ID(My_Vgrid_m,newLink)
    if(.not.newLink) then
      ! The Grid already exists in the List, so deallocate the copy
      !------------------------------------------------------------
      call deallocate_Vgrid(My_Vgrid_m)
      deallocate(My_Vgrid_m)
    endif

    ! Create a new Path structure
    !-------------------------------
    allocate(new_Path)
    new_Path%ngrid = DPlink%ngrids
    allocate(new_Path%idxHmap(new_Path%ngrid))
    allocate(new_Path%idxVmap(new_Path%ngrid))
    allocate(new_Path%ioPTR  (new_Path%ngrid))
    new_Path%idxHmap(:) = -1
    new_Path%idxVmap(:) = -1
    do nn=1,new_Path%ngrid
      new_Path%ioPTR(nn)%iodesc => null()
    end do

    ! Loop over the number of DataSet grids and add Path 
    ! Map-indicies for each DataSet grid
    !-----------------------------------------------------
    do nn=1,DPlink%ngrids
      idxHmap = -1
      idxVmap = -1

      idxHgrid_d =  DPlink%GridID(nn)%idxHgrid_d
      if(idxHgrid_d.lt.0) cycle

      Hgrid_d    => DPlink%ListHgrid_d%get_Data(idxHgrid_d)

      ! Loop over the the known Horizontal Hmaps and over Horizontal Map Opts
      !-------------------------------------------------------------------------
      do ii=1,known_Hmaps%numDataMaps
      do jj=1,known_Hmaps%DataMap(ii)%numMapOpts

        ! First look for a matching HmapOpt
        !--------------------------------------
        if(trim(HmapOpt).eq.trim(known_Hmaps%DataMap(ii)%Opts(jj))) then

          ! Now, see of the In/Out grids are compatible
          !---------------------------------------------
          Hgrid_match = chk_Hmap(Hgrid_m,Hgrid_d,trim(HmapOpt),known_Hmaps%DataMap(ii))

          if(Hgrid_match) then
            ! Get a unique ID for this Hmap. Note that the HmapData is 
            ! left empty and only filled in when the map is used for the first time.
            !----------------------------------------------------------------------
            allocate(new_HdataMap)
            new_HdataMap%idxHgrid_d       = idxHgrid_d
            new_HdataMap%idxHgrid_m       = idxHgrid_m
            new_HdataMap%Hgrid_in         => DPlink%ListHgrid_d%get_Data(idxHgrid_d)
            new_HdataMap%Hgrid_out        => DPlink%ListHgrid_m%get_Data(idxHgrid_m)
            new_HdataMap%HmapFunction     = trim(known_Hmaps%DataMap(ii)%name)
            new_HdataMap%HmapOpt          = trim(known_Hmaps%DataMap(ii)%Opts(jj))
            new_HdataMap%HmapData%HmapOpt = trim(known_Hmaps%DataMap(ii)%Opts(jj))

            ! Init any mapping data needed (HmapData for the RECTILINEAR mapping).
            !---------------------------------------------------------------------
            result = init_Hmap(Hgrid_m,Hgrid_d,trim(new_HdataMap%HmapFunction), &
                                               trim(new_HdataMap%HmapOpt),      &
                                                    new_HdataMap%HmapData       )

            idxHmap = DPlink%ListHmaps%set_ID(new_HdataMap,newLink)
            if(.not.newLink) then
              ! The Map already exists in the List, so deallocate the copy
              !------------------------------------------------------------
              deallocate(new_HdataMap)
            endif

            ! We have a Hmap index, now lets look for possible Vmaps 
            !----------------------------------------------------------
            if(DPlink%GridID(nn)%nvgrids.eq.0) then
              ! this is a Horizonal Only Grid
              !------------------------------
              allocate(new_Vmap)
              new_Vmap%idxVgrid_d   =  DPlink%GridID(nn)%idxHORIZ_ONLY
              new_Vmap%idxVgrid_m   =  idxVgrid_m
              new_Vmap%Vgrid_in     => DPlink%ListVgrid_d%get_Data(new_Vmap%idxVgrid_d)
              new_Vmap%Vgrid_out    => DPlink%ListVgrid_m%get_Data(new_Vmap%idxVgrid_m)
              new_Vmap%VmapOpt      =  "HORIZ_ONLY"
              new_Vmap%VmapFunction =  "HORIZ_ONLY"

              idxVmap = DPlink%ListVmaps%set_ID(new_Vmap,newLink)
              if(.not.newLink) then
                ! The Map already exists in the List, so deallocate the copy
                !------------------------------------------------------------
                deallocate(new_Vmap)
              endif
            else
              ! Otherwise, loop over possible Vgrids in DataSet
              !-------------------------------------------------
              do mm=1,DPlink%GridID(nn)%nvgrids
                idxVgrid_d = DPlink%GridID(nn)%idxVgrid_d(mm)
                if(idxVgrid_d.lt.0) cycle

                Vgrid_d => DPlink%ListVgrid_d%get_Data(idxVgrid_d)
                  
                ! Scan Vmaps and Options for a match
                !-----------------------------------------------
                do kk=1,known_Vmaps%numVmaps
                do ll=1,known_Vmaps%Vmap(kk)%numMapOpts

                  ! First look for a matching VmapOpt
                  !--------------------------------------
                  if(trim(VmapOpt).eq.trim(known_Vmaps%Vmap(kk)%Opts(ll))) then

                    ! Now, see of the In/Out grids are compatible
                    !---------------------------------------------
                    Mlen=len_trim(VmapOpt)
                    Vgrid_match = chk_Vmap(Vgrid_m,Vgrid_d,VmapOpt(1:Mlen),known_Vmaps%Vmap(kk))

                    if(Vgrid_match) then
                      ! Get a unique ID for this Vmap. 
                      !----------------------------------------
                      allocate(new_Vmap)
                      new_Vmap%idxVgrid_d   =  idxVgrid_d
                      new_Vmap%idxVgrid_m   =  idxVgrid_m
                      new_Vmap%Vgrid_in     => DPlink%ListVgrid_d%get_Data(idxVgrid_d)
                      new_Vmap%Vgrid_out    => DPlink%ListVgrid_m%get_Data(idxVgrid_m)
                      new_Vmap%VmapOpt      =  trim(known_Vmaps%Vmap(kk)%Opts(ll))
                      new_Vmap%VmapFunction =  trim(known_Vmaps%Vmap(kk)%FunctionName)
                      result = init_Vmap(new_Vmap)

                      idxVmap = DPlink%ListVmaps%set_ID(new_Vmap,newLink)
                      if(.not.newLink) then
                        ! The Map already exists in the List, so deallocate the copy
                        !------------------------------------------------------------
                        deallocate(new_Vmap)
                      endif
                    endif
                  endif
                end do ! ll=1,known_Vmaps%Vmap(kk)%numMapOpts
                end do ! kk=1,known_Vmaps%numVmaps

                nullify(Vgrid_d)
              end do ! mm=1,DPlink%GridID(nn)%nvgrids
            endif
          endif
        endif
      end do ! jj=1,known_Hmaps%DataMap(ii)%numMapOpts
      end do ! ii=1,known_Hmaps%numDataMaps
      nullify(Hgrid_d)
      
      ! Update Path values for this DataSet grid
      !------------------------------------------
      new_Path%idxHmap(nn) = idxHmap
      new_Path%idxVmap(nn) = idxVmap

      ! Add a PIO iodesc to the new_Path for this Grid
      !-------------------------------------------------
      if((idxHmap.gt.0).and.(idxVmap.gt.0)) then
        Hmap    => DPlink%ListHmaps%get_Data(idxHmap)
        Vmap    => DPlink%ListVmaps%get_Data(idxVmap)
        Hgrid_d => Hmap%Hgrid_in
        Vgrid_d => Vmap%Vgrid_in

        if(trim(Vgrid_d%Type).eq.'HORIZ_ONLY') then
          ! HORIZ_ONLY decomposition
          !----------------------------
          if(trim(HmapOpt).eq.'NATIVE') then
            ! Check to see if the grids match
            !----------------------------------
            dim1_found =.false.
            dim2_found =.false.
            do ii=1,Hgrid_d%nhdim
              if(trim(Hgrid_m%Hdim1name).eq.trim(Hgrid_d%dimName(ii))) dim1_found=.true.
              if(trim(Hgrid_m%Hdim2name).eq.trim(Hgrid_d%dimName(ii))) dim2_found=.true.
            end do
            if(dim1_found.and.dim2_found) then
              grids_match = .true.
            else
              grids_match = .false.
            endif

            if(.not.grids_match) then
              ! Grids don't match, so unset the H/V map indices
              !-------------------------------------------------
              new_Path%idxHmap(nn) = -1
              new_Path%idxVmap(nn) = -1
            else
              ! Use iodesc from cam_grid_support 
              ! for direct pio reads
              !-----------------------------------
              Model_dimSize(1) =      Hgrid_m%Hdim1len
              Model_dimSize(2) =      Hgrid_m%Hdim2len
              Model_dimName(1) = trim(Hgrid_m%Hdim1name)
              Model_dimName(2) = trim(Hgrid_m%Hdim2name)
  
              ndim = Hgrid_d%nhdim
              do ii=1,ndim
                Data_dimSize (ii) =      Hgrid_d%dimLen (ii)
                Data_dimName (ii) = trim(Hgrid_d%dimName(ii))
              end do
  
              ! These calls should duplicate the original(native) infld usage
              !----------------------------------------------------------------
              if(ndim.ne.2) then
                call cam_grid_get_decomp(Hgrid_m%model_gridid, Model_dimSize(1:2   ), &
                                                                Data_dimSize(1:ndim), &
                                               pio_double, new_Path%ioPTR(nn)%iodesc  )
              else
                call cam_grid_get_decomp(Hgrid_m%model_gridid, Model_dimSize(1:2   ), &
                                                                Data_dimSize(1:ndim), &
                                               pio_double, new_Path%ioPTR(nn)%iodesc, &
                                                  field_dnames=Model_dimName(1:2   )  )
              endif
            endif
          else
            ! Create an iodesc for interpolation data
            !------------------------------------------
            allocate(idof(Hmap%HmapData%npts))
            ndof = 0
            do ii=1,Hmap%HmapData%npts
              ndof = ndof +1
              idof(ndof) = Hmap%HmapData%lonInd(ii) + (Hmap%HmapData%latInd(ii)-1)*Hgrid_d%nlon
            end do
            allocate(new_Path%ioPTR(nn)%iodesc)
            call pio_initdecomp(pio_subsystem,pio_double,(/Hgrid_d%nlon,Hgrid_d%nlat/), &
                                                        idof,new_Path%ioPTR(nn)%iodesc  )
            deallocate(idof)
          endif
        else
          ! 3D decomposition
          !---------------------
          if(trim(HmapOpt).eq.'NATIVE') then
            ! Check to see if the grids match
            !----------------------------------
            dim1_found =.false.
            dim2_found =.false.
            do ii=1,Hgrid_d%nhdim
              if(trim(Hgrid_m%Hdim1name).eq.trim(Hgrid_d%dimName(ii))) dim1_found=.true.
              if(trim(Hgrid_m%Hdim2name).eq.trim(Hgrid_d%dimName(ii))) dim2_found=.true.
            end do
            if(dim1_found.and.dim2_found) then
              grids_match = .true.
            else
              grids_match = .false.
            endif

            if(.not.grids_match) then
              ! Grids don't match, so unset the H/V map indices
              !-------------------------------------------------
              new_Path%idxHmap(nn) = -1
              new_Path%idxVmap(nn) = -1
            else
              ! Use iodesc from cam_grid_support
              ! for direct pio reads
              !-----------------------------------
              if(trim(Hgrid_m%model_gridname).eq.'physgrid') then
                Model_dimSize(1) =      Hgrid_m%Hdim1len
                Model_dimSize(2) =      Vgrid_d%nlev
                Model_dimSize(3) =      Hgrid_m%Hdim2len
                Model_dimName(1) = trim(Hgrid_m%Hdim1name)
                Model_dimName(2) = trim(Vgrid_d%Vdimname )
                Model_dimName(3) = trim(Hgrid_m%Hdim2name)
              else
                Model_dimSize(1) =      Hgrid_m%Hdim1len
                Model_dimSize(2) =      Hgrid_m%Hdim2len
                Model_dimSize(3) =      Vgrid_d%nlev
                Model_dimName(1) = trim(Hgrid_m%Hdim1name)
                Model_dimName(2) = trim(Hgrid_m%Hdim2name)
                Model_dimName(3) = trim(Vgrid_d%Vdimname )
              endif
  
              ndim = Hgrid_d%nhdim
              do ii=1,ndim
                Data_dimSize (ii) =      Hgrid_d%dimLen (ii)
                Data_dimName (ii) = trim(Hgrid_d%dimName(ii))
              end do
              ndim = ndim + 1
              Data_dimSize (ndim) =      Vgrid_d%nlev
              Data_dimName (ndim) = trim(Vgrid_d%Vdimname)

              ! These calls should duplicate the original(native) infld usage
              !----------------------------------------------------------------
              if(ndim.ne.3) then
                call cam_grid_get_decomp(Hgrid_m%model_gridid, Model_dimSize(1:3   ), &
                                                                Data_dimSize(1:ndim), &
                                               pio_double, new_Path%ioPTR(nn)%iodesc, &
                                                  field_dnames=Model_dimName(1:ndim)  )
              else
                call cam_grid_get_decomp(Hgrid_m%model_gridid, Model_dimSize(1:3   ), &
                                                                Data_dimSize(1:ndim), &
                                               pio_double, new_Path%ioPTR(nn)%iodesc, &
                                                  field_dnames=Model_dimName(1:3   ), &
                                                   file_dnames= Data_dimName(1:ndim)  )
              endif
            endif
          else
            ! Create an iodesc for interpolation data
            !------------------------------------------
            allocate(idof(Hmap%HmapData%npts*Vgrid_d%nlev))
            ndof = 0
            do kk=1,Vgrid_d%nlev
            do ii=1,Hmap%HmapData%npts
              ndof = ndof +1
              idof(ndof) = Hmap%HmapData%lonInd(ii) + (Hmap%HmapData%latInd(ii)-1)*Hgrid_d%nlon &
                                                    + (kk-1)*Hgrid_d%nlon*Hgrid_d%nlat
            end do
            end do
            allocate(new_Path%ioPTR(nn)%iodesc)
            call pio_initdecomp(pio_subsystem,pio_double,(/Hgrid_d%nlon,   &
                                                           Hgrid_d%nlat,   &
                                                           Vgrid_d%nlev/), &
                                            idof,new_Path%ioPTR(nn)%iodesc )
            deallocate(idof)
          endif
        endif
        nullify(Hmap   )
        nullify(Vmap   )
        nullify(Hgrid_d)
        nullify(Vgrid_d)
      endif
    end do ! nn=1,DPlink%ngrids

    ! Get unique Data Path index
    !------------------------------
    register_dataPath = DPlink%ListPaths%set_ID(new_Path,newLink)
    if(.not.newLink) then
      ! The Map already exists in the List, so deallocate the copy
      !------------------------------------------------------------
      deallocate(new_Path%idxHmap)
      deallocate(new_Path%idxVmap)
      deallocate(new_Path%ioPTR)
      deallocate(new_Path)
    endif

    nullify(DPlink)
    nullify(My_Hgrid_m)
    nullify(My_Vgrid_m)
    nullify(new_HdataMap)
    nullify(new_Vmap)
    nullify(new_Path)

    ! End Routine
    !------------------
    return
  end function register_dataPath
  !================================================================


  !================================================================
  subroutine print_dataset_paths(DataSetID)
    ! 
    ! print_dataset_paths: 
    !==============================================================
    !
    ! Passed variables
    !--------------------
    integer,intent(in):: DataSetID
    !
    ! Local values
    !----------------
    type(LinkDataPaths_t),pointer:: DPlink
    logical                      :: found,newLink

    ! Get the DataPath link for the given DataSetID, create and 
    ! initialize  a new link if one is not found for DataSetID
    !----------------------------------------------------------
    found = inq_dataset_status(DataSetID)
    if(found) then
      DPlink => My_DataPaths%set_ID(DataSetID,newLink)
      call DPlink%printLink()
      nullify(DPlink)
    else
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    LIST of DataSet Paths: DataSetID=',DataSetID
      write(iulog,*) '          DataSet Not open '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) ' '
    endif

    ! End Routine
    !------------------
    return
  end subroutine print_dataset_paths
  !================================================================


  !================================================================
  subroutine apply_dataPath_2D(DataSetID,idxPath,I_VarName,timelevel,fillvalue, &
                                                  dim1b,dim1e,dim2b,dim2e,Output)
    ! 
    ! apply_dataPath_2D: 
    !==============================================================
    !
    ! Passed variables
    !--------------------
    integer         ,intent(in ):: DataSetID
    integer         ,intent(in ):: idxPath
    character(len=*),intent(in ):: I_VarName
    integer         ,intent(in ):: timelevel
    real(r8)        ,intent(out):: fillvalue
    integer         ,intent(in ):: dim1b,dim1e,dim2b,dim2e
    real(r8)        ,intent(out):: Output(dim1b:dim1e,dim2b:dim2e)
    !
    ! Local values
    !----------------
    logical            :: result
    logical            :: found
    logical            :: newLink
    character(len=80)  :: VarName
    integer            :: idxVar
    type(dataset_var_t):: Var
    integer            :: idxGrid

    real(r8),allocatable:: H_Data(:,:)
    real(r8),allocatable:: V_Data(:,:)

    type(LinkDataPaths_t),pointer:: DPlink
    type(Path_t)         ,pointer:: My_Path
    type(HdataMap_t)     ,pointer:: Hmap
    type(vmap_t)         ,pointer:: Vmap
    type(dataset_hgrid_t),pointer:: Hgrid_d 
    type(packed_hgrid_t) ,pointer:: Hgrid_m
    type(vgrid_t)        ,pointer:: Vgrid_d
    type(PathDesc_t)             :: DPdesc
    integer                      :: ierr

    ! Check that the DataSetID in an open/registered dataset
    !-------------------------------------------------------- 
    found = inq_dataset_status(DataSetID)
    if(.not.found) then
      write(iulog,*) ' DataSetID=',DataSetID
      call endrun('apply_datapath_2D: ERROR Unknown DataSetID')
    endif
    
    ! From the given VarName, get the corrsponding VarID and gridID
    !---------------------------------------------------------------
    VarName = trim(I_VarName)
    found = inq_dataset_varid(DataSetID,VarName,idxVar,Var,idxGrid)
    if(.not.found) then
      write(iulog,*) ' VarName='//trim(VarName)
      call endrun('apply_datapath_2D: ERROR Unknown VarName')
    endif

    ! Get the Path info for the given pathID and gridID
    !-------------------------------------------------------------------
    DPlink => My_DataPaths%set_ID(DataSetID,newLink)
    if(newLink) then
      write(iulog,*) ' newLink=',newLink
      call endrun('apply_datapath_2D: ERROR DPlink is a newLink')
      ! *Remove the new DataSet
      !  call My_DataPaths%rm_DataSet(DataSetID)
      !
      ! ** ERROR STOP??
    endif

    My_Path => DPlink%ListPaths%get_Data(idxPath)
    if((idxGrid.le.0).or.(idxGrid.gt.My_Path%ngrid)) then
      write(iulog,*) ' idxGrid=',idxGrid,'My_Path%ngrid=',My_Path%ngrid
      call endrun('apply_datapath_2D: ERROR idxGrid is out of range')
    endif

    ! Initialize DataSet PathDesc datastructure
    !---------------------------------------------
    Hmap    => DPlink%ListHmaps%get_Data(My_Path%idxHmap(idxGrid))
    Vmap    => DPlink%ListVmaps%get_Data(My_Path%idxVmap(idxGrid))
    Hgrid_d => Hmap%Hgrid_in
    Hgrid_m => Hmap%Hgrid_out
    Vgrid_d => Vmap%Vgrid_in

    DPdesc%ncid      =  inq_dataset_ncid(DataSetID)
    DPdesc%ncid_nf90 =  inq_dataset_ncid_nf90(DataSetID)   ! PFC DIAG Temp to test PIO ERROR with nf90
    DPdesc%varname   =      trim(Var%name)                 ! PFC DIAG Temp to test PIO ERROR with nf90
    DPdesc%vlen      =  len_trim(Var%name)                 ! PFC DIAG Temp to test PIO ERROR with nf90
    DPdesc%timelevel =  timelevel
    DPdesc%ncol      =  Hgrid_m%ncol
    DPdesc%nlev_d    =  Vgrid_d%nlev
    DPdesc%nlev_m    =  -1
    DPdesc%fillvalue =  Var%fill_value
    DPdesc%VarDesc   =  Var%VarDesc
    DPdesc%iodesc    => My_Path%ioPTR(idxGrid)%iodesc

    ! This variable must be on a HORIZ_ONLY grid
    !--------------------------------------------
    if(trim(Vgrid_d%Type).eq.'HORIZ_ONLY') then
      if(trim(Hmap%HmapOpt).eq.'NATIVE') then
        ! For NATIVE grids, read the data in here 
        !---------------------------------------------
        call pio_setframe   (DPdesc%ncid,DPdesc%VarDesc,int(timelevel,kind=pio_offset_kind))
        call pio_read_darray(DPdesc%ncid,DPdesc%VarDesc,DPdesc%iodesc,Output,ierr)
      else
        ! Allocate space and apply Horizontal mapping
        !----------------------------------------------------
        DPdesc%nlev_d=1
        allocate(H_Data(DPdesc%ncol,DPdesc%nlev_d))
  
        result = apply_Hmap(H_Data,DPdesc,trim(Hmap%HmapFunction),             &
                                          trim(Hmap%HmapOpt     ),Hmap%HmapData)

        ! Unpack the data into the PE array
        !-----------------------------------
        call unpack_model_Hgrid_2D(Hgrid_m,H_Data,dim1b,dim1e,dim2b,dim2e,Output)
        deallocate(H_Data)
      endif
    else
     ! This is a Horiz only interface and the variable is not Horiz_only
     !-------------------------------------------------------------------
     call endrun('apply_datapath_2D: ERROR Vgrid_d is not HORIZ_ONLY')
    endif

    ! Return the fill_value from the DataSet
    !---------------------------------------
    fillvalue = DPdesc%fillvalue

    ! Clean up 
    !------------
    nullify(DPlink )
    nullify(My_Path)
    nullify(Hmap   )
    nullify(Vmap   )
    nullify(Hgrid_d)
    nullify(Hgrid_m)
    nullify(Vgrid_d)

    ! End Routine
    !------------------
    return
  end subroutine apply_dataPath_2D
  !================================================================


  !================================================================
  subroutine apply_dataPath_3D(DataSetID,idxPath,I_VarName,timelevel,fillvalue, &
                                      dim1b,dim1e,dim2b,dim2e,dim3b,dim3e,Output)
    ! 
    ! apply_dataPath: 
    !==============================================================
    !
    ! Passed variables
    !--------------------
    integer         ,intent(in ):: DataSetID
    integer         ,intent(in ):: idxPath
    character(len=*),intent(in ):: I_VarName
    integer         ,intent(in ):: timelevel
    real(r8)        ,intent(out):: fillvalue
    integer         ,intent(in ):: dim1b,dim1e,dim2b,dim2e,dim3b,dim3e
    real(r8)        ,intent(out):: Output(dim1b:dim1e,dim2b:dim2e,dim3b:dim3e)
    !
    ! Local values
    !----------------
    logical            :: result
    logical            :: found
    logical            :: newLink
    character(len=80)  :: VarName
    integer            :: idxVar
    type(dataset_var_t):: Var
    integer            :: idxGrid

    real(r8),allocatable:: H_Data(:,:)
    real(r8),allocatable:: V_Data(:,:)

    type(LinkDataPaths_t),pointer:: DPlink
    type(Path_t)         ,pointer:: My_Path
    type(HdataMap_t)     ,pointer:: Hmap
    type(vmap_t)         ,pointer:: Vmap
    type(dataset_hgrid_t),pointer:: Hgrid_d 
    type(packed_hgrid_t) ,pointer:: Hgrid_m
    type(vgrid_t)        ,pointer:: Vgrid_d
    type(vgrid_t)        ,pointer:: Vgrid_m
    type(PathDesc_t)             :: DPdesc
    integer                      :: ierr

!PFC -----------------------------------------------------------------------
  integer :: Ival,ntim0,nlev0,nlat0,nlon0,nn,ll
  real(r8):: Lon0,Lat0,My_diff
!PFC -----------------------------------------------------------------------

    ! Check that the DataSetID in an open/registered dataset
    !-------------------------------------------------------- 
    found = inq_dataset_status(DataSetID)
    if(.not.found) then
      write(iulog,*) ' DataSetID=',DataSetID
      call endrun('apply_datapath_3D: ERROR Unknown DataSetID')
    endif
    
    ! From the given VarName, get the corrsponding VarID and gridID
    !---------------------------------------------------------------
    VarName = trim(I_VarName)
    found = inq_dataset_varid(DataSetID,VarName,idxVar,Var,idxGrid)
    if(.not.found) then
      write(iulog,*) ' VarName='//trim(VarName)
      call endrun('apply_datapath_3D: ERROR Unknown VarName')
    endif

    ! Get the Path info for the given pathID and gridID
    !-------------------------------------------------------------------
    DPlink => My_DataPaths%set_ID(DataSetID,newLink)
    if(newLink) then
      write(iulog,*) ' newLink=',newLink
      call endrun('apply_datapath_3D: ERROR DPlink is a newLink')
      ! *Remove the new DataSet
      !  call My_DataPaths%rm_DataSet(DataSetID)
      !
      ! ** ERROR STOP??
    endif

    My_Path => DPlink%ListPaths%get_Data(idxPath)
    if((idxGrid.le.0).or.(idxGrid.gt.My_Path%ngrid)) then
      write(iulog,*) ' idxGrid=',idxGrid,'My_Path%ngrid=',My_Path%ngrid
      call endrun('apply_datapath_3D: ERROR idxGrid is out of range')
    endif

    ! Initialize DataSet PathDesc datastructure
    !---------------------------------------------
    Hmap    => DPlink%ListHmaps%get_Data(My_Path%idxHmap(idxGrid))
    Vmap    => DPlink%ListVmaps%get_Data(My_Path%idxVmap(idxGrid))
    Hgrid_d => Hmap%Hgrid_in
    Hgrid_m => Hmap%Hgrid_out
    Vgrid_d => Vmap%Vgrid_in
    Vgrid_m => Vmap%Vgrid_out

    DPdesc%ncid      =  inq_dataset_ncid(DataSetID)
    DPdesc%ncid_nf90 =  inq_dataset_ncid_nf90(DataSetID)   ! PFC DIAG Temp to test PIO ERROR with nf90
    DPdesc%varname   =      trim(Var%name)                 ! PFC DIAG Temp to test PIO ERROR with nf90
    DPdesc%vlen      =  len_trim(Var%name)                 ! PFC DIAG Temp to test PIO ERROR with nf90
    DPdesc%timelevel =  timelevel
    DPdesc%ncol      =  Hgrid_m%ncol
    DPdesc%nlev_d    =  Vgrid_d%nlev
    DPdesc%nlev_m    =  Vgrid_m%nlev
    DPdesc%fillvalue =  Var%fill_value
    DPdesc%VarDesc   =  Var%VarDesc
    DPdesc%iodesc    => My_Path%ioPTR(idxGrid)%iodesc

    ! This variable must NOT be on a HORIZ_ONLY grid
    !------------------------------------------------
    if(trim(Vgrid_d%Type).eq.'HORIZ_ONLY') then
      ! This is a 3D only interface and the variable is Horiz_only
      !-------------------------------------------------------------------
      call endrun('apply_datapath_3D: ERROR Vgrid_d is HORIZ_ONLY')
    else
      if(trim(Hmap%HmapOpt).eq.'NATIVE') then
        ! For NATIVE grids, read the data in here     TODO: FIX Output dims not right for general use.....
        !---------------------------------------------
        call pio_setframe   (DPdesc%ncid,DPdesc%VarDesc,int(timelevel,kind=pio_offset_kind))
        call pio_read_darray(DPdesc%ncid,DPdesc%VarDesc,DPdesc%iodesc,Output,ierr)

        if(trim(Vmap%VmapOpt).ne.'INDEX') then
          ! Optionally apply a vertical mapping to the 
          ! NATIVE horizontal grid data
          !-----------------------------------------------
          call endrun('apply_datapath_3D: VmapOpt for NATIVE grid not implemented yet')

          allocate(H_Data(DPdesc%ncol,DPdesc%nlev_d))
          allocate(V_Data(DPdesc%ncol,DPdesc%nlev_m))

          ! Pack data into (Hgrid_m,Vgrid_d) array
          !----------------------------------------
!          call pack_model_Hgrid_3D(Hgrid_m,H_Data,Vgrid_d%nlev,dim1b,dim1e,      &
!                                                               dim2b,dim2e,      &
!                                                               dim3b,dim3e,RawData)

          ! Apply the vertical mapping to (Hgrid_m,Vgrid_m)
          !------------------------------------------------
          result = apply_Vmap(V_Data,H_Data,DPdesc,trim(Vmap%VmapFunction),             &
                                                   trim(Vmap%VmapOpt     ),Vmap%VmapData)

          ! Unpack the data into the PE array
          !-----------------------------------
          call unpack_model_Hgrid_3D(Hgrid_m,V_Data,Vgrid_m%nlev,dim1b,dim1e,      &
                                                                 dim2b,dim2e,      &
                                                                 dim3b,dim3e,Output)
          deallocate(V_Data)
          deallocate(H_Data)
        endif
      else
        ! Apply the Horizontal mapping (Read in data) then apply 
        ! the Vertical mapping to the result
        !--------------------------------------------------------
        allocate(H_Data(DPdesc%ncol,DPdesc%nlev_d))
        allocate(V_Data(DPdesc%ncol,DPdesc%nlev_m))

        result = apply_Hmap(H_Data,DPdesc,trim(Hmap%HmapFunction),             &
                                          trim(Hmap%HmapOpt     ),Hmap%HmapData)

        result = apply_Vmap(V_Data,H_Data,DPdesc,trim(Vmap%VmapFunction),             &
                                                 trim(Vmap%VmapOpt     ),Vmap%VmapData)

        ! Unpack the data into the PE array
        !-----------------------------------
        call unpack_model_Hgrid_3D(Hgrid_m,V_Data,Vgrid_m%nlev,dim1b,dim1e,      &
                                                               dim2b,dim2e,      &
                                                               dim3b,dim3e,Output)
        deallocate(V_Data)
        deallocate(H_Data)
      endif
    endif

    ! Return the fill_value from the DataSet
    !---------------------------------------
    fillvalue = DPdesc%fillvalue

    ! Clean up 
    !------------
    nullify(DPlink )
    nullify(My_Path)
    nullify(Hmap   )
    nullify(Vmap   )
    nullify(Hgrid_d)
    nullify(Hgrid_m)
    nullify(Vgrid_d)
    nullify(Vgrid_m)

    ! End Routine
    !------------------
    return
  end subroutine apply_dataPath_3D
  !================================================================


  !================================================================
  subroutine apply_horiz_dataPath_2D(I_VarName,DataSetID,idxPath,timelevel,  &
                                               dim1b,dim1e,dim2b,dim2e,H_Data)
    ! 
    ! apply_horiz_dataPath_2D: 
    !==============================================================
    !
    ! Passed variables
    !--------------------
    integer         ,intent(in ):: DataSetID
    integer         ,intent(in ):: idxPath
    character(len=*),intent(in ):: I_VarName
    integer         ,intent(in ):: timelevel
    integer         ,intent(in ):: dim1b,dim1e,dim2b,dim2e
    real(r8)        ,intent(out):: H_Data(:,:)
    !
    ! Local values
    !----------------
    logical            :: result
    logical            :: found
    logical            :: newLink
    character(len=80)  :: VarName
    integer            :: idxVar
    type(dataset_var_t):: Var
    integer            :: idxGrid

    real(r8),allocatable:: Output(:,:)

    type(LinkDataPaths_t),pointer:: DPlink
    type(Path_t)         ,pointer:: My_Path
    type(HdataMap_t)     ,pointer:: Hmap
    type(vmap_t)         ,pointer:: Vmap
    type(dataset_hgrid_t),pointer:: Hgrid_d 
    type(packed_hgrid_t) ,pointer:: Hgrid_m
    type(vgrid_t)        ,pointer:: Vgrid_d
    type(PathDesc_t)             :: DPdesc
    integer                      :: ierr

    ! Check that the DataSetID in an open/registered dataset
    !-------------------------------------------------------- 
    found = inq_dataset_status(DataSetID)
    if(.not.found) then
      write(iulog,*) ' DataSetID=',DataSetID
      call endrun('apply_horiz_datapath_2D: ERROR Unknown DataSetID')
    endif
    
    ! From the given VarName, get the corrsponding VarID and gridID
    !---------------------------------------------------------------
    VarName = trim(I_VarName)
    found = inq_dataset_varid(DataSetID,VarName,idxVar,Var,idxGrid)
    if(.not.found) then
      write(iulog,*) ' VarName='//trim(VarName)
      call endrun('apply_horiz_datapath_2D: ERROR Unknown VarName')
    endif

    ! Get the Path info for the given pathID and gridID
    !-------------------------------------------------------------------
    DPlink => My_DataPaths%set_ID(DataSetID,newLink)
    if(newLink) then
      write(iulog,*) ' newLink=',newLink
      call endrun('apply_horiz_datapath_2D: ERROR DPlink is a newLink')
    endif

    My_Path => DPlink%ListPaths%get_Data(idxPath)
    if((idxGrid.le.0).or.(idxGrid.gt.My_Path%ngrid)) then
      write(iulog,*) ' idxGrid=',idxGrid,'My_Path%ngrid=',My_Path%ngrid
      call endrun('apply_horiz_datapath_2D: ERROR idxGrid is out of range')
    endif

    ! Initialize DataSet PathDesc datastructure
    !---------------------------------------------
    Hmap    => DPlink%ListHmaps%get_Data(My_Path%idxHmap(idxGrid))
    Vmap    => DPlink%ListVmaps%get_Data(My_Path%idxVmap(idxGrid))
    Hgrid_d => Hmap%Hgrid_in
    Hgrid_m => Hmap%Hgrid_out
    Vgrid_d => Vmap%Vgrid_in

    DPdesc%ncid      =  inq_dataset_ncid(DataSetID)
    DPdesc%ncid_nf90 =  inq_dataset_ncid_nf90(DataSetID)   ! PFC DIAG Temp to test PIO ERROR with nf90
    DPdesc%varname   =      trim(Var%name)                 ! PFC DIAG Temp to test PIO ERROR with nf90
    DPdesc%vlen      =  len_trim(Var%name)                 ! PFC DIAG Temp to test PIO ERROR with nf90
    DPdesc%timelevel =  timelevel
    DPdesc%ncol      =  Hgrid_m%ncol
    DPdesc%nlev_d    =  Vgrid_d%nlev
    DPdesc%nlev_m    =  -1
    DPdesc%fillvalue =  Var%fill_value
    DPdesc%VarDesc   =  Var%VarDesc
    DPdesc%iodesc    => My_Path%ioPTR(idxGrid)%iodesc

    ! This variable must be on a HORIZ_ONLY grid
    !--------------------------------------------
    if(trim(Vgrid_d%Type).eq.'HORIZ_ONLY') then
      if(trim(Hmap%HmapOpt).eq.'NATIVE') then
        ! For NATIVE grids, read the data in here 
        !---------------------------------------------
        allocate(Output(dim1b:dim1e,dim2b:dim2e))
        call pio_setframe   (DPdesc%ncid,DPdesc%VarDesc,int(timelevel,kind=pio_offset_kind))
        call pio_read_darray(DPdesc%ncid,DPdesc%VarDesc,DPdesc%iodesc,Output,ierr)

        call pack_model_Hgrid_2D(Hgrid_m,H_Data,dim1b,dim1e,dim2b,dim2e,Output)
        deallocate(Output)
      else
        ! Allocate space and apply Horizontal mapping
        !----------------------------------------------------
        DPdesc%nlev_d=1
        result = apply_Hmap(H_Data,DPdesc,trim(Hmap%HmapFunction),             &
                                          trim(Hmap%HmapOpt     ),Hmap%HmapData)
      endif
    else
     ! This is a Horiz only interface and the variable is not Horiz_only
     !-------------------------------------------------------------------
     call endrun('apply_horiz_datapath_2D: ERROR Vgrid_d is not HORIZ_ONLY')
    endif

    ! Clean up 
    !------------
    nullify(DPlink )
    nullify(My_Path)
    nullify(Hmap   )
    nullify(Vmap   )
    nullify(Hgrid_d)
    nullify(Hgrid_m)
    nullify(Vgrid_d)

    ! End Routine
    !------------------
    return
  end subroutine apply_horiz_dataPath_2D
  !================================================================


  !================================================================
  subroutine apply_horiz_dataPath_3D(I_VarName,DataSetID,idxPath,timelevel,    &
                                     dim1b,dim1e,dim2b,dim2e,dim3b,dim3e,H_Data)
    ! 
    ! apply_horiz_dataPath: 
    !==============================================================
    !
    ! Passed variables
    !--------------------
    integer         ,intent(in ):: DataSetID
    integer         ,intent(in ):: idxPath
    character(len=*),intent(in ):: I_VarName
    integer         ,intent(in ):: timelevel
    integer         ,intent(in ):: dim1b,dim1e,dim2b,dim2e,dim3b,dim3e
    real(r8)        ,intent(out):: H_Data(:,:)
    !
    ! Local values
    !----------------
    logical            :: result
    logical            :: found
    logical            :: newLink
    character(len=80)  :: VarName
    integer            :: idxVar
    type(dataset_var_t):: Var
    integer            :: idxGrid

    real(r8),allocatable:: Output(:,:,:)

    type(LinkDataPaths_t),pointer:: DPlink
    type(Path_t)         ,pointer:: My_Path
    type(HdataMap_t)     ,pointer:: Hmap
    type(vmap_t)         ,pointer:: Vmap
    type(dataset_hgrid_t),pointer:: Hgrid_d 
    type(packed_hgrid_t) ,pointer:: Hgrid_m
    type(vgrid_t)        ,pointer:: Vgrid_d
    type(vgrid_t)        ,pointer:: Vgrid_m
    type(PathDesc_t)             :: DPdesc
    integer                      :: ierr

    ! Check that the DataSetID in an open/registered dataset
    !-------------------------------------------------------- 
    found = inq_dataset_status(DataSetID)
    if(.not.found) then
      write(iulog,*) ' DataSetID=',DataSetID
      call endrun('apply_horiz_datapath_3D: ERROR Unknown DataSetID')
    endif
    
    ! From the given VarName, get the corrsponding VarID and gridID
    !---------------------------------------------------------------
    VarName = trim(I_VarName)
    found = inq_dataset_varid(DataSetID,VarName,idxVar,Var,idxGrid)
    if(.not.found) then
      write(iulog,*) ' VarName='//trim(VarName)
      call endrun('apply_horiz_datapath_3D: ERROR Unknown VarName')
    endif

    ! Get the Path info for the given pathID and gridID
    !-------------------------------------------------------------------
    DPlink => My_DataPaths%set_ID(DataSetID,newLink)
    if(newLink) then
      write(iulog,*) ' newLink=',newLink
      call endrun('apply_horiz_datapath_3D: ERROR DPlink is a newLink')
    endif

    My_Path => DPlink%ListPaths%get_Data(idxPath)
    if((idxGrid.le.0).or.(idxGrid.gt.My_Path%ngrid)) then
      write(iulog,*) ' idxGrid=',idxGrid,'My_Path%ngrid=',My_Path%ngrid
      call endrun('apply_horiz_datapath_3D: ERROR idxGrid is out of range')
    endif

    ! Initialize DataSet PathDesc datastructure
    !---------------------------------------------
    Hmap    => DPlink%ListHmaps%get_Data(My_Path%idxHmap(idxGrid))
    Vmap    => DPlink%ListVmaps%get_Data(My_Path%idxVmap(idxGrid))
    Hgrid_d => Hmap%Hgrid_in
    Hgrid_m => Hmap%Hgrid_out
    Vgrid_d => Vmap%Vgrid_in
    Vgrid_m => Vmap%Vgrid_out

    DPdesc%ncid      =  inq_dataset_ncid(DataSetID)
    DPdesc%ncid_nf90 =  inq_dataset_ncid_nf90(DataSetID)   ! PFC DIAG Temp to test PIO ERROR with nf90
    DPdesc%varname   =      trim(Var%name)                 ! PFC DIAG Temp to test PIO ERROR with nf90
    DPdesc%vlen      =  len_trim(Var%name)                 ! PFC DIAG Temp to test PIO ERROR with nf90
    DPdesc%timelevel =  timelevel
    DPdesc%ncol      =  Hgrid_m%ncol
    DPdesc%nlev_d    =  Vgrid_d%nlev
    DPdesc%nlev_m    =  Vgrid_m%nlev
    DPdesc%fillvalue =  Var%fill_value
    DPdesc%VarDesc   =  Var%VarDesc
    DPdesc%iodesc    => My_Path%ioPTR(idxGrid)%iodesc

    ! This variable must NOT be on a HORIZ_ONLY grid
    !------------------------------------------------
    if(trim(Vgrid_d%Type).eq.'HORIZ_ONLY') then
      ! This is a 3D only interface and the variable is Horiz_only
      !-------------------------------------------------------------------
      call endrun('apply_horiz_datapath_3D: ERROR Vgrid_d is HORIZ_ONLY')
    else
      if(trim(Hmap%HmapOpt).eq.'NATIVE') then
        ! For NATIVE grids, read the data in here, then pack it into Hgrid_m format
        !---------------------------------------------------------------------------
        allocate(Output(dim1b:dim1e,dim2b:dim2e,dim3b:dim3e))
        call pio_setframe   (DPdesc%ncid,DPdesc%VarDesc,int(timelevel,kind=pio_offset_kind))
        call pio_read_darray(DPdesc%ncid,DPdesc%VarDesc,DPdesc%iodesc,Output,ierr)

        call pack_model_Hgrid_3D(Hgrid_m,H_Data,DPdesc%nlev_d,             &
                                 dim1b,dim1e,dim2b,dim2e,dim3b,dim3e,Output)
        deallocate(Output)
      else
        ! Apply the Horizontal mapping (Read in data) 
        !--------------------------------------------------------
        result = apply_Hmap(H_Data,DPdesc,trim(Hmap%HmapFunction),             &
                                          trim(Hmap%HmapOpt     ),Hmap%HmapData)
      endif
    endif

    ! Clean up 
    !------------
    nullify(DPlink )
    nullify(My_Path)
    nullify(Hmap   )
    nullify(Vmap   )
    nullify(Hgrid_d)
    nullify(Hgrid_m)
    nullify(Vgrid_d)
    nullify(Vgrid_m)

    ! End Routine
    !------------------
    return
  end subroutine apply_horiz_dataPath_3D
  !================================================================


  !================================================================
  subroutine apply_vert_dataPath(I_VarName,DataSetID,idxPath,timelevel,H_Data,V_Data)
    ! 
    ! apply_vert_dataPath: 
    !==============================================================
    !
    ! Passed variables
    !--------------------
    integer         ,intent(in   ):: DataSetID
    integer         ,intent(in   ):: idxPath
    character(len=*),intent(in   ):: I_VarName
    integer         ,intent(in   ):: timelevel
    real(r8)        ,intent(inout):: H_Data(:,:)
    real(r8)        ,intent(out  ):: V_Data(:,:)
    !
    ! Local values
    !----------------
    logical            :: result
    logical            :: found
    logical            :: newLink
    character(len=80)  :: VarName
    integer            :: idxVar
    type(dataset_var_t):: Var
    integer            :: idxGrid

    type(LinkDataPaths_t),pointer:: DPlink
    type(Path_t)         ,pointer:: My_Path
    type(HdataMap_t)     ,pointer:: Hmap
    type(vmap_t)         ,pointer:: Vmap
    type(dataset_hgrid_t),pointer:: Hgrid_d 
    type(packed_hgrid_t) ,pointer:: Hgrid_m
    type(vgrid_t)        ,pointer:: Vgrid_d
    type(vgrid_t)        ,pointer:: Vgrid_m
    type(PathDesc_t)             :: DPdesc
    integer                      :: ierr

    ! Check that the DataSetID in an open/registered dataset
    !-------------------------------------------------------- 
    found = inq_dataset_status(DataSetID)
    if(.not.found) then
      write(iulog,*) ' DataSetID=',DataSetID
      call endrun('apply_datapath_3D: ERROR Unknown DataSetID')
    endif
    
    ! From the given VarName, get the corrsponding VarID and gridID
    !---------------------------------------------------------------
    VarName = trim(I_VarName)
    found = inq_dataset_varid(DataSetID,VarName,idxVar,Var,idxGrid)
    if(.not.found) then
      write(iulog,*) ' VarName='//trim(VarName)
      call endrun('apply_datapath_3D: ERROR Unknown VarName')
    endif

    ! Get the Path info for the given pathID and gridID
    !-------------------------------------------------------------------
    DPlink => My_DataPaths%set_ID(DataSetID,newLink)
    if(newLink) then
      write(iulog,*) ' newLink=',newLink
      call endrun('apply_datapath_3D: ERROR DPlink is a newLink')
    endif

    My_Path => DPlink%ListPaths%get_Data(idxPath)
    if((idxGrid.le.0).or.(idxGrid.gt.My_Path%ngrid)) then
      write(iulog,*) ' idxGrid=',idxGrid,'My_Path%ngrid=',My_Path%ngrid
      call endrun('apply_datapath_3D: ERROR idxGrid is out of range')
    endif

    ! Initialize DataSet PathDesc datastructure
    !---------------------------------------------
    Hmap    => DPlink%ListHmaps%get_Data(My_Path%idxHmap(idxGrid))
    Vmap    => DPlink%ListVmaps%get_Data(My_Path%idxVmap(idxGrid))
    Hgrid_d => Hmap%Hgrid_in
    Hgrid_m => Hmap%Hgrid_out
    Vgrid_d => Vmap%Vgrid_in
    Vgrid_m => Vmap%Vgrid_out

    DPdesc%ncid      =  inq_dataset_ncid(DataSetID)
    DPdesc%ncid_nf90 =  inq_dataset_ncid_nf90(DataSetID)   ! PFC DIAG Temp to test PIO ERROR with nf90
    DPdesc%varname   =      trim(Var%name)                 ! PFC DIAG Temp to test PIO ERROR with nf90
    DPdesc%vlen      =  len_trim(Var%name)                 ! PFC DIAG Temp to test PIO ERROR with nf90
    DPdesc%timelevel =  timelevel
    DPdesc%ncol      =  Hgrid_m%ncol
    DPdesc%nlev_d    =  Vgrid_d%nlev
    DPdesc%nlev_m    =  Vgrid_m%nlev
    DPdesc%fillvalue =  Var%fill_value
    DPdesc%VarDesc   =  Var%VarDesc
    DPdesc%iodesc    => My_Path%ioPTR(idxGrid)%iodesc

    ! This variable must NOT be on a HORIZ_ONLY grid
    !------------------------------------------------
    if(trim(Vgrid_d%Type).eq.'HORIZ_ONLY') then
      ! This is a 3D only interface and the variable is Horiz_only
      !-------------------------------------------------------------------
      call endrun('apply_vert_dataPath: ERROR Vgrid_d is HORIZ_ONLY')
    else
      ! Apply the Horizontal mapping (Read in data) then apply 
      ! the Vertical mapping to the result
      !--------------------------------------------------------
      result = apply_Vmap(V_Data,H_Data,DPdesc,trim(Vmap%VmapFunction),             &
                                               trim(Vmap%VmapOpt     ),Vmap%VmapData)
    endif

    ! Clean up 
    !------------
    nullify(DPlink )
    nullify(My_Path)
    nullify(Hmap   )
    nullify(Vmap   )
    nullify(Hgrid_d)
    nullify(Hgrid_m)
    nullify(Vgrid_d)
    nullify(Vgrid_m)

    ! End Routine
    !------------------
    return
  end subroutine apply_vert_dataPath
  !================================================================


  !================================================================
  subroutine unpack_dataPath_result_2D(I_VarName,DataSetID,idxPath,timelevel, &
                                       dim1b,dim1e,dim2b,dim2e,DP_Data,Output )
    ! 
    ! unpack_dataPath_result_2D: 
    !==============================================================
    !
    ! Passed variables
    !--------------------
    integer         ,intent(in ):: DataSetID
    integer         ,intent(in ):: idxPath
    character(len=*),intent(in ):: I_VarName
    integer         ,intent(in ):: timelevel
    integer         ,intent(in ):: dim1b,dim1e,dim2b,dim2e
    real(r8)        ,intent(in ):: DP_Data(:,:)
    real(r8)        ,intent(out):: Output(dim1b:dim1e,dim2b:dim2e)
    !
    ! Local values
    !----------------
    logical            :: found
    logical            :: newLink
    character(len=80)  :: VarName
    integer            :: idxVar
    type(dataset_var_t):: Var
    integer            :: idxGrid

    type(LinkDataPaths_t),pointer:: DPlink
    type(Path_t)         ,pointer:: My_Path
    type(HdataMap_t)     ,pointer:: Hmap
    type(packed_hgrid_t) ,pointer:: Hgrid_m

    ! Check that the DataSetID in an open/registered dataset
    !-------------------------------------------------------- 
    found = inq_dataset_status(DataSetID)
    if(.not.found) then
      write(iulog,*) ' DataSetID=',DataSetID
      call endrun('apply_datapath_3D: ERROR Unknown DataSetID')
    endif
    
    ! From the given VarName, get the corrsponding VarID and gridID
    !---------------------------------------------------------------
    VarName = trim(I_VarName)
    found = inq_dataset_varid(DataSetID,VarName,idxVar,Var,idxGrid)
    if(.not.found) then
      write(iulog,*) ' VarName='//trim(VarName)
      call endrun('apply_datapath_3D: ERROR Unknown VarName')
    endif

    ! Get the Path info for the given pathID and gridID
    !-------------------------------------------------------------------
    DPlink => My_DataPaths%set_ID(DataSetID,newLink)
    if(newLink) then
      write(iulog,*) ' newLink=',newLink
      call endrun('apply_datapath_3D: ERROR DPlink is a newLink')
    endif

    My_Path => DPlink%ListPaths%get_Data(idxPath)
    if((idxGrid.le.0).or.(idxGrid.gt.My_Path%ngrid)) then
      write(iulog,*) ' idxGrid=',idxGrid,'My_Path%ngrid=',My_Path%ngrid
      call endrun('apply_datapath_3D: ERROR idxGrid is out of range')
    endif

    ! Initialize DataSet PathDesc datastructure
    !---------------------------------------------
    Hmap    => DPlink%ListHmaps%get_Data(My_Path%idxHmap(idxGrid))
    Hgrid_m => Hmap%Hgrid_out

    ! Unpack the data into the PE array
    !-----------------------------------
    call unpack_model_Hgrid_2D(Hgrid_m,DP_Data,dim1b,dim1e,      &
                                               dim2b,dim2e,Output)
    ! Clean up 
    !------------
    nullify(DPlink )
    nullify(My_Path)
    nullify(Hmap   )
    nullify(Hgrid_m)

    ! End Routine
    !------------------
    return
  end subroutine unpack_dataPath_result_2D
  !================================================================


  !================================================================
  subroutine unpack_dataPath_result_3D(I_VarName,DataSetID,idxPath,timelevel,            &
                                       dim1b,dim1e,dim2b,dim2e,dim3b,dim3e,DP_Data,Output)
    ! 
    ! unpack_dataPath_result_3D: 
    !==============================================================
    !
    ! Passed variables
    !--------------------
    integer         ,intent(in ):: DataSetID
    integer         ,intent(in ):: idxPath
    character(len=*),intent(in ):: I_VarName
    integer         ,intent(in ):: timelevel
    integer         ,intent(in ):: dim1b,dim1e,dim2b,dim2e,dim3b,dim3e
    real(r8)        ,intent(in ):: DP_Data(:,:)
    real(r8)        ,intent(out):: Output(dim1b:dim1e,dim2b:dim2e,dim3b:dim3e)
    !
    ! Local values
    !----------------
    logical            :: found
    logical            :: newLink
    character(len=80)  :: VarName
    integer            :: idxVar
    type(dataset_var_t):: Var
    integer            :: idxGrid

    type(LinkDataPaths_t),pointer:: DPlink
    type(Path_t)         ,pointer:: My_Path
    type(HdataMap_t)     ,pointer:: Hmap
    type(vmap_t)         ,pointer:: Vmap
    type(packed_hgrid_t) ,pointer:: Hgrid_m
    type(vgrid_t)        ,pointer:: Vgrid_m

    ! Check that the DataSetID in an open/registered dataset
    !-------------------------------------------------------- 
    found = inq_dataset_status(DataSetID)
    if(.not.found) then
      write(iulog,*) ' DataSetID=',DataSetID
      call endrun('apply_datapath_3D: ERROR Unknown DataSetID')
    endif
    
    ! From the given VarName, get the corrsponding VarID and gridID
    !---------------------------------------------------------------
    VarName = trim(I_VarName)
    found = inq_dataset_varid(DataSetID,VarName,idxVar,Var,idxGrid)
    if(.not.found) then
      write(iulog,*) ' VarName='//trim(VarName)
      call endrun('apply_datapath_3D: ERROR Unknown VarName')
    endif

    ! Get the Path info for the given pathID and gridID
    !-------------------------------------------------------------------
    DPlink => My_DataPaths%set_ID(DataSetID,newLink)
    if(newLink) then
      write(iulog,*) ' newLink=',newLink
      call endrun('apply_datapath_3D: ERROR DPlink is a newLink')
    endif

    My_Path => DPlink%ListPaths%get_Data(idxPath)
    if((idxGrid.le.0).or.(idxGrid.gt.My_Path%ngrid)) then
      write(iulog,*) ' idxGrid=',idxGrid,'My_Path%ngrid=',My_Path%ngrid
      call endrun('apply_datapath_3D: ERROR idxGrid is out of range')
    endif

    ! Initialize DataSet PathDesc datastructure
    !---------------------------------------------
    Hmap    => DPlink%ListHmaps%get_Data(My_Path%idxHmap(idxGrid))
    Vmap    => DPlink%ListVmaps%get_Data(My_Path%idxVmap(idxGrid))
    Hgrid_m => Hmap%Hgrid_out
    Vgrid_m => Vmap%Vgrid_out

    ! Unpack the data into the PE array
    !-----------------------------------
    call unpack_model_Hgrid_3D(Hgrid_m,DP_Data,Vgrid_m%nlev,dim1b,dim1e,      &
                                                            dim2b,dim2e,      &
                                                            dim3b,dim3e,Output)
    ! Clean up 
    !------------
    nullify(DPlink )
    nullify(My_Path)
    nullify(Hmap   )
    nullify(Vmap   )
    nullify(Hgrid_m)
    nullify(Vgrid_m)

    ! End Routine
    !------------------
    return
  end subroutine unpack_dataPath_result_3D
  !================================================================


  !================================================================
  subroutine pack_model_Hgrid_2D(Hgrid_m,Hdat_Out,dim1b,dim1e,          &
                                                  dim2b,dim2e,Hdat_Model)
    ! 
    ! pack_model_Hgrid: Copy data on Model (PE) grid to output values
    !                   for the given Packed_Model grid Hgrid_m
    !                   
    !==============================================================
    !
    ! Passed Variables
    !--------------------
    type(packed_hgrid_t),intent(in ):: Hgrid_m
    real(r8)            ,intent(out):: Hdat_Out(:,:)
    integer             ,intent(in ):: dim1b,dim1e,dim2b,dim2e
    real(r8)            ,intent(in ):: Hdat_Model(dim1b:dim1e,dim2b:dim2e)
    !
    ! Local Values
    !---------------
    integer:: nn

    ! Zero everything first
    !--------------------------------------------------------------------------
    Hdat_Out(1:Hgrid_m%ncol,1) = 0._r8

    ! For 2D, The maps for each PackType are applied the same 
    !-----------------------------------------------------------
    do nn=1,Hgrid_m%ncol
      Hdat_Out(nn,1) = Hdat_Model(Hgrid_m%Xmap(nn),Hgrid_m%Ymap(nn))
    end do

    ! End Routine
    !------------------
    return
  end subroutine pack_model_Hgrid_2D
  !================================================================


  !================================================================
  subroutine pack_model_Hgrid_3D(Hgrid_m,Hdat_Out,nlev,dim1b,dim1e,          &
                                                       dim2b,dim2e,          &
                                                       dim3b,dim3e,Hdat_Model)
    ! 
    ! pack_model_Hgrid: Copy data on Model (PE) grid to output values
    !                   for the given Packed_Model grid Hgrid_m
    !                   
    !==============================================================
    !
    ! Passed Variables
    !--------------------
    type(packed_hgrid_t),intent(in ):: Hgrid_m
    real(r8)            ,intent(out):: Hdat_Out(:,:)
    integer             ,intent(in ):: nlev
    integer             ,intent(in ):: dim1b,dim1e,dim2b,dim2e,dim3b,dim3e
    real(r8)            ,intent(in ):: Hdat_Model(dim1b:dim1e,dim2b:dim2e,dim3b:dim3e)
    !
    ! Local Values
    !---------------
    integer:: nn,ll

    ! Zero everything first
    !--------------------------------------------------------------------------
    Hdat_Out(1:Hgrid_m%ncol,1:nlev) = 0._r8

    ! The apply maps for each PackType 
    !----------------------------------------------
    if(trim(Hgrid_m%PackType).eq.'UNSTRUCTURED') then
      do ll=1,nlev
      do nn=1,Hgrid_m%ncol
        Hdat_Out(nn,ll) = Hdat_Model(Hgrid_m%Xmap(nn),dim2b+ll-1,Hgrid_m%Ymap(nn))
      end do
      end do
    elseif(trim(Hgrid_m%PackType).eq.'LATLON') then
      do ll=1,nlev
      do nn=1,Hgrid_m%ncol
        Hdat_Out(nn,ll) = Hdat_Model(Hgrid_m%Xmap(nn),Hgrid_m%Ymap(nn),dim3b+ll-1)
      end do
      end do
    else
      call endrun('pack_model_Hgrid_3D: ERROR Unknown PackType = '//trim(Hgrid_m%PackType))
    endif

    ! End Routine
    !------------------
    return
  end subroutine pack_model_Hgrid_3D
  !================================================================


  !================================================================
  subroutine unpack_model_Hgrid_2D(Hgrid_m,Hdat_In,dim1b,dim1e,          &
                                                   dim2b,dim2e,Hdat_Model)
    ! 
    ! unpack_model_Hgrid: Copy input data from the Packed_Model grid to 
    !                     output values for the given Model (PE) grid.
    !
    !==============================================================
    !
    ! Passed Variables
    !--------------------
    type(packed_hgrid_t),intent(in ):: Hgrid_m
    real(r8)            ,intent(in ):: Hdat_In(:,:)
    integer             ,intent(in ):: dim1b,dim1e,dim2b,dim2e
    real(r8)            ,intent(out):: Hdat_Model(dim1b:dim1e,dim2b:dim2e)
    !
    ! Local Values
    !---------------
    integer:: nn

    ! Not gauranteed that all PE gridpoints are used, so zero everything first
    !--------------------------------------------------------------------------
    Hdat_Model(dim1b:dim1e,dim2b:dim2e) = 0._r8

    ! For 2D, The maps for each PackType are applied the same 
    !-----------------------------------------------------------
    do nn=1,Hgrid_m%ncol
      Hdat_Model(Hgrid_m%Xmap(nn),Hgrid_m%Ymap(nn)) = Hdat_In(nn,1)
    end do

    ! End Routine
    !------------------
    return
  end subroutine unpack_model_Hgrid_2D
  !================================================================


  !================================================================
  subroutine unpack_model_Hgrid_3D(Hgrid_m,Hdat_In,nlev,dim1b,dim1e,          &
                                                        dim2b,dim2e,          &
                                                        dim3b,dim3e,Hdat_Model)
    ! 
    ! unpack_model_Hgrid: Copy input data from the Packed_Model grid to 
    !                     output values for the given Model (PE) grid.
    !
    !==============================================================
    !
    ! Passed Variables
    !--------------------
    type(packed_hgrid_t),intent(in ):: Hgrid_m
    real(r8)            ,intent(in ):: Hdat_In(:,:)
    integer             ,intent(in ):: nlev
    integer             ,intent(in ):: dim1b,dim1e,dim2b,dim2e,dim3b,dim3e
    real(r8)            ,intent(out):: Hdat_Model(dim1b:dim1e,dim2b:dim2e,dim3b:dim3e)
    !
    ! Local Values
    !---------------
    integer:: nn,ll

    ! Not gauranteed that all PE gridpoints are used, so zero everything first
    !--------------------------------------------------------------------------
    Hdat_Model(dim1b:dim1e,dim2b:dim2e,dim3b:dim3e) = 0._r8

    ! The apply maps for each PackType 
    !----------------------------------------------
    if(trim(Hgrid_m%PackType).eq.'UNSTRUCTURED') then
      do ll=1,nlev
      do nn=1,Hgrid_m%ncol
        Hdat_Model(Hgrid_m%Xmap(nn),dim2b+ll-1,Hgrid_m%Ymap(nn)) = Hdat_In(nn,ll)
      end do
      end do
    elseif(trim(Hgrid_m%PackType).eq.'LATLON') then
      do ll=1,nlev
      do nn=1,Hgrid_m%ncol
        Hdat_Model(Hgrid_m%Xmap(nn),Hgrid_m%Ymap(nn),dim3b+ll-1) = Hdat_In(nn,ll)
      end do
      end do
    else
      call endrun('unpack_model_Hgrid_3D: ERROR Unknown PackType = '//trim(Hgrid_m%PackType))
    endif

    ! End Routine
    !------------------
    return
  end subroutine unpack_model_Hgrid_3D
  !================================================================


  !================================================================
  logical function inq_dataPath_dims(I_VarName,DataSetID,idxPath,ncol,nlev_m,nlev_d)
    ! 
    ! inq_dataPath_dims: For the given Dataset ID, Variable name, and pathID
    !                    return the dimensions of the associated arrays. 
    !==============================================================
    !
    ! Passed variables
    !--------------------
    character(len=*),intent(in ):: I_VarName
    integer         ,intent(in ):: DataSetID
    integer         ,intent(in ):: idxPath
    integer         ,intent(out):: ncol
    integer         ,intent(out):: nlev_m
    integer         ,intent(out):: nlev_d
    !
    ! Local values
    !----------------
    logical            :: found
    logical            :: newLink
    character(len=80)  :: VarName
    integer            :: idxVar
    type(dataset_var_t):: Var
    integer            :: idxGrid

    type(LinkDataPaths_t),pointer:: DPlink
    type(Path_t)         ,pointer:: My_Path
    type(HdataMap_t)     ,pointer:: Hmap
    type(vmap_t)         ,pointer:: Vmap
    type(packed_hgrid_t) ,pointer:: Hgrid_m
    type(vgrid_t)        ,pointer:: Vgrid_m
    type(vgrid_t)        ,pointer:: Vgrid_d

    ! Check that the DataSetID in an open/registered dataset
    !-------------------------------------------------------- 
    found = inq_dataset_status(DataSetID)
    if(.not.found) then
      write(iulog,*) ' DataSetID=',DataSetID
      call endrun('inq_dataPath_dims: ERROR Unknown DataSetID')
    endif
    
    ! From the given VarName, get the corrsponding VarID and gridID
    !---------------------------------------------------------------
    VarName = trim(I_VarName)
    found = inq_dataset_varid(DataSetID,VarName,idxVar,Var,idxGrid)
    if(.not.found) then
      write(iulog,*) ' VarName='//trim(VarName)
      call endrun('inq_dataPath_dims: ERROR Unknown VarName')
    endif

    ! Get the Path info for the given pathID and gridID
    !-------------------------------------------------------------------
    DPlink => My_DataPaths%set_ID(DataSetID,newLink)
    if(newLink) then
      write(iulog,*) ' newLink=',newLink
      call endrun('inq_dataPath_dims: ERROR DPlink is a newLink')
    endif

    My_Path => DPlink%ListPaths%get_Data(idxPath)
    if((idxGrid.le.0).or.(idxGrid.gt.My_Path%ngrid)) then
      write(iulog,*) ' idxGrid=',idxGrid,'My_Path%ngrid=',My_Path%ngrid
      call endrun('inq_dataPath_dims: ERROR idxGrid is out of range')
    endif

    ! Initialize DataSet PathDesc datastructure
    !---------------------------------------------
    Hmap    => DPlink%ListHmaps%get_Data(My_Path%idxHmap(idxGrid))
    Vmap    => DPlink%ListVmaps%get_Data(My_Path%idxVmap(idxGrid))
    Hgrid_m => Hmap%Hgrid_out
    Vgrid_d => Vmap%Vgrid_in
    if(trim(Vgrid_d%Type).eq.'HORIZ_ONLY') then
      ncol   =  Hgrid_m%ncol
      nlev_m =  1
      nlev_d =  1
    else
      Vgrid_m => Vmap%Vgrid_out
      ncol    =  Hgrid_m%ncol
      nlev_m  =  Vgrid_m%nlev
      nlev_d  =  Vgrid_d%nlev
    endif
    inq_dataPath_dims = .true.

    nullify(DPlink)
    nullify(My_Path)
    nullify(Hmap)
    nullify(Vmap)
    nullify(Hgrid_m)
    nullify(Vgrid_m)
    nullify(Vgrid_d)

    ! End Function
    !-------------
    return
  end function inq_dataPath_dims
  !================================================================


  !========================================================================
  !========================================================================
  ! Procedures For ListDataPaths Structure
  !
  !========================================================================
  !========================================================================
  subroutine ListDataPaths_reset(this)
    ! Reset the Iterator to point at the first DataPath
    !----------------------------------------------------
    class(ListDataPaths_t):: this
    this%curr => this%first
    return
  end subroutine ListDataPaths_reset
  !========================================================================
  subroutine ListDataPaths_iterate(this)
    ! Increment the Iterator to point at the next DataPath
    !----------------------------------------------------
    class(ListDataPaths_t):: this
    this%curr => this%curr%get_Next()
    return
  end subroutine ListDataPaths_iterate
  !========================================================================
  function ListDataPaths_isFirst(this)
    ! Check if the Iterator points to the end of the list.
    !----------------------------------------------------
    class(ListDataPaths_t)        :: this
    logical                       :: ListDataPaths_isFirst
    class(LinkDataPaths_t),pointer:: lastPointer
    if(associated(this%curr)) then
      lastPointer => this%curr%get_Prev()
      ListDataPaths_isFirst = .not.associated(lastPointer)
    else
      ListDataPaths_isFirst = .true.
    endif
    return
  end function ListDataPaths_isFirst
  !========================================================================
  function ListDataPaths_isLast(this)
    ! Check if the Iterator points to the end of the list.
    !----------------------------------------------------
    class(ListDataPaths_t)        :: this
    logical                       :: ListDataPaths_isLast
    class(LinkDataPaths_t),pointer:: nextPointer
    if(associated(this%curr)) then
      nextPointer => this%curr%get_Next()
      ListDataPaths_isLast = .not.associated(nextPointer)
    else
      ListDataPaths_isLast = .true.
    endif
    return
  end function ListDataPaths_isLast
  !========================================================================
  function ListDataPaths_getLinkID(this)
    ! Return the ID DataSet for curr in List
    !---------------------------------------------
    class(ListDataPaths_t):: this
    integer               :: ListDataPaths_getLinkID
    ListDataPaths_getLinkID = this%curr%get_ID()
    return
  end function ListDataPaths_getLinkID
  !========================================================================
  function ListDataPaths_setLinkID(this,DataSetID,newLink)
    ! Scan the List for the Link with the matching DataSetID, if not found 
    ! then create a new Link and add it to the end of the List. 
    ! Return a pointer to the Link for the existing DataSetID or to 
    ! the newly created Link for DataSetID.
    ! Set newLink to indicate if this is a new or existing link
    !--------------------------------------------------------------------
    class(ListDataPaths_t)        :: this
    integer                       :: DataSetID
    class(LinkDataPaths_t),pointer:: ListDataPaths_setLinkID
    logical                       :: newLink
    !
    ! Local Values
    !-------------
    class(LinkDataPaths_t),pointer:: new_DPlink
    type(dataset_t)       ,pointer:: DataSet
    type(dataset_hgrid_t) ,pointer:: new_Hgrid_d
    type(vgrid_t)         ,pointer:: new_Vgrid_d
    logical                       :: link_found
    logical                       :: grid_match, Vgrid_found, match_found
    integer                       :: grid_kk
    integer                       :: nn,kk,dd,vv,jj,ll,ii
    integer,allocatable           :: VgridID(:,:)
    integer                       :: idxHORIZ_ONLY
    integer                       :: idxVdim

    ! Cycle thru the list to find the DataSet with matching ID
    !----------------------------------------------------------
    call this%reset
    link_found = .false.
    do while(associated(this%curr)) 
      if(this%curr%get_ID().eq. DataSetID) then
        link_found = .true.
        exit
      endif
      call this%iterate
    end do

    if(link_found) then
      ! Return the current link for existing DataSetID
      !-------------------------------------------------
      newLink = .false.
      ListDataPaths_setLinkID => this%curr
    else
      ! We need to create a new DataSet Link and all of
      ! the Lists contained in it.
      !--------------------------------------------------
      newLink = .true.
      allocate(new_DPlink)
      new_DPlink%DataSetID = DataSetID
      allocate(new_DPlink%ListHgrid_m)
      allocate(new_DPlink%ListVgrid_m)
      allocate(new_DPlink%ListHgrid_d)
      allocate(new_DPlink%ListVgrid_d)
      allocate(new_DPlink%ListHmaps  )
      allocate(new_DPlink%ListVmaps  )
      allocate(new_DPlink%ListPaths  )
      new_DPlink%ListHgrid_m%Counter = 0
      new_DPlink%ListHgrid_m%Counter = 0
      new_DPlink%ListVgrid_m%Counter = 0
      new_DPlink%ListHgrid_d%Counter = 0
      new_DPlink%ListVgrid_d%Counter = 0
      new_DPlink%ListHmaps%Counter   = 0
      new_DPlink%ListVmaps%Counter   = 0
      new_DPlink%ListPaths%Counter   = 0

      ! Get the meta-data structure for the given DataSetID
      ! and initialize the lists of known dataset grids.
      !-----------------------------------------------------
      DataSet => inq_dataset(DataSetID)

      ! Assume each DataSet dimension can be the vertical coordinate
      ! (wasteful but avoids cumbersome logic). Each dimension will have 
      ! a generic INDEX vertical coordinate and each dimesnsion can have 
      ! multile Vertical grids associated with it.
      !------------------------------------------------------------------
      allocate(VgridID(known_Vgrids%numVgrids,DataSet%ndims))
      VgridID(:,:) = -1

      ! For each dimesion, loop thru known_Vgrids and look for matches
      !---------------------------------------------------------------
      do nn=1,DataSet%ndims
      do kk=1,known_Vgrids%numVgrids
        if(trim(known_Vgrids%Vgrid(kk)%Type).eq.'INDEX') then
          ! Every dimension will have an INDEX coordinate
          !================================================
          allocate(new_Vgrid_d)
          new_Vgrid_d%name        = trim(DataSet%dims(nn)%name)
          new_Vgrid_d%Type        = 'INDEX'
          new_Vgrid_d%Vdimname    = trim(DataSet%dims(nn)%name)
          new_Vgrid_d%nlev        = DataSet%dims(nn)%len
          new_Vgrid_d%ngrid_var   = 0
          new_Vgrid_d%nreq_var    = 0

          ! Add the new Vgrid_d to the List
          !---------------------------------
          VgridID(kk,nn) = new_DPlink%ListVgrid_d%set_ID(new_Vgrid_d,newLink)

          ! If this is not a new link
          ! deallocate the new_Vgrid_d structure
          !-------------------------------------------------------
          if(.not.newLink) then
            ! Clean up new_Vgrid_d  
            !-----------------------
            call deallocate_Vgrid(new_Vgrid_d)
            deallocate(new_Vgrid_d)
          endif
        else
          ! Scan the othe known_Vgrid templates to see if the DataSet supports it.
          !========================================================================
          Vgrid_found = .true.

          ! Loop over grid_var() requirements for the known_Vgrids
          !--------------------------------------------------------
          do jj=1,known_Vgrids%Vgrid(kk)%ngrid_var

            ! Can't rely on DataSets having the PS0 refernce value, 
            ! so we must 'hardwire' in the value. 
            !----------------------------------------------------------
            if(trim(known_Vgrids%Vgrid(kk)%grid_var(jj)).eq.'Pref') then
              ! Set match found for this grid_var
              !------------------------------------
              match_found = .true.
            else
              ! See if grid_var() is present
              !------------------------------------
              match_found = .false.
              do ll=1,DataSet%dims(nn)%ndvar
!PFC                if(trim(known_Vgrids%Vgrid(kk)%grid_var(jj)).eq. &
!PFC                   trim(      DataSet%dims(nn)%var(ll)%name)     ) then
                if((trim(known_Vgrids%Vgrid(kk)%grid_var(jj)).eq.                  &
                    trim(      DataSet%dims(nn)%var(ll)%name)     ).or.            &
                   ((trim(known_Vgrids%Vgrid(kk)%grid_var(jj)).eq.'PRESSURE').and. &
                    (trim(      DataSet%dims(nn)%var(ll)%name).eq.'level'   )     )) then
                  match_found = .true.
                endif
              end do ! ll=1,DataSet%dims(nn)%ndvar
            endif
            if(.not.match_found) Vgrid_found = .false.
          end do ! jj=1,known_Vgrids%Vgrid(kk)%ngrid_var

          ! Now loop over req_var() and see if the required 
          !   values are present in the DataSet
          !----------------------------------------------
!PP          do jj=1,known_Vgrids%Vgrid(kk)%nreq_var
!PP            match_found = inq_dataset_var(DataSet,trim(known_Vgrids%Vgrid(kk)%req_var(jj)))
!PP            if(.not.match_found) Vgrid_found = .false.
!PP          end do 

          if(Vgrid_found) then
            ! A match was found. Create a new Vgrid_d structure and fill in values
            !----------------------------------------------------------------------
            allocate(new_Vgrid_d)
            new_Vgrid_d%name        = trim(known_Vgrids%Vgrid(kk)%name)
            new_Vgrid_d%Type        = trim(known_Vgrids%Vgrid(kk)%Type)
            new_Vgrid_d%Vdimname    = trim(known_Vgrids%Vgrid(kk)%Vdimname)
            new_Vgrid_d%ngrid_var   =      known_Vgrids%Vgrid(kk)%ngrid_var
            if(new_Vgrid_d%ngrid_var.gt.0) then
              allocate(new_Vgrid_d%grid_var(new_Vgrid_d%ngrid_var))
              do ii=1,new_Vgrid_d%ngrid_var
                new_Vgrid_d%grid_var(ii) =trim(known_Vgrids%Vgrid(kk)%grid_var(ii))
              end do
            endif
            new_Vgrid_d%nreq_var    =      known_Vgrids%Vgrid(kk)%nreq_var
            if(new_Vgrid_d%nreq_var.gt.0) then
              allocate(new_Vgrid_d%req_var(new_Vgrid_d%nreq_var))
              do ii=1,new_Vgrid_d%nreq_var
                new_Vgrid_d%req_var(ii) =trim(known_Vgrids%Vgrid(kk)%req_var(ii))
              end do
            endif

            ! Hardwire the refernce Pref value 
            !-----------------------------------
            new_Vgrid_d%Pref = PS0_REF

            ! Fill in Vgrid_d values from the DataSet
            !----------------------------------------
            new_Vgrid_d%nlev = DataSet%dims(nn)%len
            do ll=1,DataSet%dims(nn)%ndvar
              if((trim(DataSet%dims(nn)%var(ll)%name).eq.'hyam').or. &
                 (trim(DataSet%dims(nn)%var(ll)%name).eq.'hyai')     ) then
                allocate(new_Vgrid_d%HYA(new_Vgrid_d%nlev))
                new_Vgrid_d%HYA(:) = DataSet%dims(nn)%var(ll)%val(:)
              endif
              if((trim(DataSet%dims(nn)%var(ll)%name).eq.'hybm').or. &
                 (trim(DataSet%dims(nn)%var(ll)%name).eq.'hybi')     ) then
                allocate(new_Vgrid_d%HYB(new_Vgrid_d%nlev))
                new_Vgrid_d%HYB(:) = DataSet%dims(nn)%var(ll)%val(:)
              endif
              if((trim(DataSet%dims(nn)%var(ll)%name).eq.'pressure').or. &
                 (trim(DataSet%dims(nn)%var(ll)%name).eq.'level')        ) then
                allocate(new_Vgrid_d%PRES(new_Vgrid_d%nlev))
                new_Vgrid_d%PRES(:) = DataSet%dims(nn)%var(ll)%val(:)
              endif
              if(trim(DataSet%dims(nn)%var(ll)%name).eq.'height') then
                allocate(new_Vgrid_d%HGHT(new_Vgrid_d%nlev))
                new_Vgrid_d%HGHT(:) = DataSet%dims(nn)%var(ll)%val(:)
              endif
            end do

            ! Add the new Vgrid_d to the List
            !---------------------------------
            VgridID(kk,nn) = new_DPlink%ListVgrid_d%set_ID(new_Vgrid_d,newLink)
  
            ! If this is not a new link
            ! deallocate the new_Vgrid_d structure
            !-------------------------------------------------------
            if(.not.newLink) then
              ! Clean up new_Vgrid_d  
              !-----------------------
              call deallocate_Vgrid(new_Vgrid_d)
              deallocate(new_Vgrid_d)
            endif
          endif

        endif
      end do ! kk=1,known_Vgrids%numVgrids
      end do ! nn=1,DataSet%ndims

      ! Add a HORIZONTAL ONLY grid to accomodate 2D surface data
      !----------------------------------------------------------
      allocate(new_Vgrid_d)
      new_Vgrid_d%name      = 'HORIZ_ONLY'
      new_Vgrid_d%Type      = 'HORIZ_ONLY'
      new_Vgrid_d%Vdimname  = 'HORIZ_ONLY'
      new_Vgrid_d%nlev      = 1
      new_Vgrid_d%ngrid_var = 0
      new_Vgrid_d%nreq_var  = 0

      idxHORIZ_ONLY = new_DPlink%ListVgrid_d%set_ID(new_Vgrid_d,newLink)
      if(.not.newLink) then
        ! Clean up new_Vgrid_d  
        !-----------------------
        call deallocate_Vgrid(new_Vgrid_d)
        deallocate(new_Vgrid_d)
      endif

      ! Create an array of unique-known horizontal/vertical grid indices associated 
      ! with each dataset gridID. Initialize values before adding grid info.
      !     Note: For unknown grids: idxHgrid_d = -1 
      !                              idxVgrid_d = -1
      !                              idxHORZ_ONLY is a NULL Vgrid for Horiz data arrays
      !
      ! Initialize each grid for HORIZ_ONLY grid values. 
      ! The possible Vertical grid options will be added as needed in the next section. 
      !---------------------------------------------------------------------------------
      new_DPlink%ngrids = DataSet%ngrids
      allocate(new_DPlink%GridID(new_DPlink%ngrids))
      do nn=1,new_DPlink%ngrids
        new_DPlink%GridID(nn)%idxHgrid_d    = -1
        new_DPlink%GridID(nn)%idxHORIZ_ONLY = idxHORIZ_ONLY
        new_DPlink%GridID(nn)%nvgrids       = 0
      end do

      ! Loop over grids found in the DataSet
      !---------------------------------------
      do nn=1,DataSet%ngrids

        ! Use the tables of known Horizontal grids to initialize the
        ! Hgrid_d list of unique dataset grids found.
        !-----------------------------------------------------------
        do kk=1,known_Hgrids%numDataGrids
          if(DataSet%grids(nn)%nsdim.lt.known_Hgrids%DataGrid(kk)%nhdim) then
            ! Not enough spatial dims to qualify
            !-----------------------------------
            grid_match = .false.
          else
            grid_match = .true.
            grid_kk    = kk
            do dd=1,known_Hgrids%DataGrid(kk)%nhdim
!PFC 
              ii = len_trim(known_Hgrids%DataGrid(kk)%dimName(dd))
              if(trim(         known_Hgrids%DataGrid(kk)%dimName(dd)).ne. &
                 DataSet%dims(DataSet%grids(nn)%dimID(dd))%name(1:ii)     ) then
!PFC              if(trim(         known_Hgrids%DataGrid(kk)%dimName(dd)).ne. &
!PFC                 trim(DataSet%dims(DataSet%grids(nn)%dimID(dd))%name)     ) then
                grid_match = .false.
              endif
            end do
          endif

          ! If a match was found....
          !---------------------------
          if(grid_match) then
            ! fill out the contents of a new Hgrid_d from known_Hgrids values.
            !------------------------------------------------------------------
            allocate(new_Hgrid_d)
            new_Hgrid_d%datasetID  = DataSetID
            new_Hgrid_d%name       = trim(known_Hgrids%DataGrid(grid_kk)%name)
            new_Hgrid_d%Type       = trim(known_Hgrids%DataGrid(grid_kk)%Type)
            new_Hgrid_d%nhdim      =      known_Hgrids%DataGrid(grid_kk)%nhdim
            new_Hgrid_d%dimName(1) = trim(known_Hgrids%DataGrid(grid_kk)%dimName(1))
            new_Hgrid_d%dimName(2) = trim(known_Hgrids%DataGrid(grid_kk)%dimName(2))
            new_Hgrid_d%dimID(1:2) =      known_Hgrids%DataGrid(grid_kk)%dimID(1:2)

            ! Add dataset dimID's for each dimension
            !----------------------------------------
            do dd=1,new_Hgrid_d%nhdim
            do vv=1,DataSet%ndims
              if(trim(new_Hgrid_d%dimName(dd)).eq.trim(DataSet%dims(vv)%name)) then
                new_Hgrid_d%dimID (dd) = DataSet%dims(vv)%dimid
                new_Hgrid_d%dimLen(dd) = DataSet%dims(vv)%len
                exit
              endif
            end do
            end do
            
            ! Load Lat/Lon values from DataSet
            !---------------------------------
            do vv=1,DataSet%ndims
            do dd=1,DataSet%dims(vv)%ndvar
              if((trim(DataSet%dims(vv)%var(dd)%name).eq.'lon'      ).or. &
                 (trim(DataSet%dims(vv)%var(dd)%name).eq.'slon'     ).or. &
                 (trim(DataSet%dims(vv)%var(dd)%name).eq.'longitude')     )then
                new_Hgrid_d%nlon = DataSet%dims(vv)%len
                allocate(new_Hgrid_d%Lon(new_Hgrid_d%nlon))
                new_Hgrid_d%Lon(:) = DataSet%dims(vv)%var(dd)%val(:)
              endif
              if((trim(DataSet%dims(vv)%var(dd)%name).eq.'lat'     ).or. &
                 (trim(DataSet%dims(vv)%var(dd)%name).eq.'slat'    ).or. &
                 (trim(DataSet%dims(vv)%var(dd)%name).eq.'latitude')     )then
                new_Hgrid_d%nlat = DataSet%dims(vv)%len
                allocate(new_Hgrid_d%Lat(new_Hgrid_d%nlat))
                new_Hgrid_d%Lat(:) = DataSet%dims(vv)%var(dd)%val(:)
              endif
            end do
            end do

            ! If this grid has a Vertical-Type index, then allocate 
            ! and set the known idxVgrid_d options.
            !-------------------------------------------------------
            if(DataSet%grids(nn)%nsdim.eq.(new_Hgrid_d%nhdim+1)) then
              new_DPlink%GridID(nn)%nvgrids = known_Vgrids%numVgrids
              allocate(new_DPlink%GridID(nn)%idxVgrid_d(new_DPlink%GridID(nn)%nvgrids))
              idxVdim = DataSet%grids(nn)%dimID(DataSet%grids(nn)%nsdim)
              new_DPlink%GridID(nn)%idxVgrid_d(:) = VgridID(:,idxVdim)
            endif

            ! Add unique grids to the Hgrid_d List
            ! there can only be one matching grid, so exit the loop
            !--------------------------------------------------------
            new_DPlink%GridID(nn)%idxHgrid_d = new_DPlink%ListHgrid_d%set_ID(new_Hgrid_d,newLink)

            ! If this is not a new link
            ! deallocate the new_Hgrid_d structure
            !-------------------------------------------------------
            if(.not.newLink) then
              ! Clean up new_Hgrid_d   
              !-----------------------
              call deallocate_Hgrid(new_Hgrid_d)
              deallocate(new_Hgrid_d)
            endif
            exit
          endif ! (grid_match)
        end do ! kk=1,known_Hgrids%numDataGrids
      end do ! nn=1,DataSet%ngrids

      ! Now add the new Link to the end of the List
      !---------------------------------------------
      if(.not.associated(this%last)) then
        ! If the current list is empty, initialize 
        ! the list with this DataSet Link
        !-------------------------------------------
        new_DPlink%prev => null()
        new_DPlink%next => null()
        this%first      => new_DPlink
        this%last       => new_DPlink
        this%curr       => new_DPlink
      else
        ! Append this DataSet Link to the end of the current list.
        !----------------------------------------------------
        new_DPlink%prev => this%last
        new_DPlink%next => null()
        this%last%next  => new_DPlink
        this%last       => new_DPlink
        this%curr       => this%last
      endif

      ! Return the newly created link for DataSetID
      !---------------------------------------------
      ListDataPaths_setLinkID => this%last

      ! Final cleanup
      !---------------
      deallocate(VgridID)
    endif

    return
  end function ListDataPaths_setLinkID
  !========================================================================
  subroutine ListDataPaths_print(this)
    ! Cycle thru the current list and print 
    ! out the meta data for each link 
    !---------------------------------------------
    class(ListDataPaths_t):: this
    !
    ! Local values
    !-------------
    logical:: empty_list

    ! Check to see if the list is empty
    !-----------------------------------
    call this%reset
    empty_list = .false.
    if(this%isFirst().and.this%isLast()) then
      if(.not.associated(this%curr)) then
        empty_list = .true.
      endif
    endif
    
    if(empty_list) then
      ! Print out "Empty List" message here
      !---------------------------------------
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    LIST of DataPath Links:'
      write(iulog,*) '          EMPTY LIST  '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) ' '
    else
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    LIST of DataPath Links:'
      write(iulog,*) ' '

      ! cycle thru the list
      !----------------------
      do while(associated(this%curr))
        ! Issue a print command for each link
        !-------------------------------------
        call this%curr%printLink
        call this%iterate
      end do
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) ' '
    endif

    return
  end subroutine ListDataPaths_print
  !========================================================================



  !========================================================================
  !========================================================================
  ! Procedures For LinkDataPaths Structure
  !
  !========================================================================
  !========================================================================
  function LinkDataPaths_getID(this)
    ! Return the ID for this DataPath
    !----------------------------------------------------
    class(LinkDataPaths_t):: this
    integer               :: LinkDataPaths_getID
    LinkDataPaths_getID = this%DataSetID
    return
  end function LinkDataPaths_getID
  !========================================================================
  function LinkDataPaths_getNext(this)
    ! Return the pointer to the next DataPath
    !----------------------------------------------------
    class(LinkDataPaths_t)        :: this
    class(LinkDataPaths_t),pointer:: LinkDataPaths_getNext
    LinkDataPaths_getNext => this%next
    return
  end function LinkDataPaths_getNext
  !========================================================================
  subroutine LinkDataPaths_setNext(this,next)
    ! Set the pointer for the next DataPath
    !----------------------------------------------------
    class(LinkDataPaths_t)        :: this
    class(LinkDataPaths_t),pointer:: next
    this%next => next
    return
  end subroutine LinkDataPaths_setNext
  !========================================================================
  function LinkDataPaths_getPrev(this)
    ! Return the pointer to the previous DataPath
    !----------------------------------------------------
    class(LinkDataPaths_t)        :: this
    class(LinkDataPaths_t),pointer:: LinkDataPaths_getPrev
    LinkDataPaths_getPrev => this%prev
    return
  end function LinkDataPaths_getPrev
  !========================================================================
  subroutine LinkDataPaths_setPrev(this,prev)
    ! Return the pointer to the previous DataPath
    !----------------------------------------------------
    class(LinkDataPaths_t)        :: this
    class(LinkDataPaths_t),pointer:: prev
    this%prev => prev
    return
  end subroutine LinkDataPaths_setPrev
  !========================================================================
  subroutine LinkDataPaths_print(this)
    ! 
    ! Print out the DataPath metadata for the current link
    !----------------------------------------------------
    class(LinkDataPaths_t):: this
    !
    ! Local Values
    !--------------
    integer:: ii

    ! Print out identifying values
    !--------------------------------
    write(iulog,*) '   ============================================================'
    write(iulog,*) '    DataPath Link METADATA'
    write(iulog,*) '      -------------------------------------------'
    write(iulog,*) '      DataSetID=',this%DataSetID
    write(iulog,*) ' '
    write(iulog,*) '           Table of Unique H/V Grids for each DataSet grid '
    write(iulog,*) '               ngrids= ',this%ngrids
    do ii=1,this%ngrids
      write(iulog,*) '               DataSet Grid ii=',ii,' idxHORIZ_ONLY=',this%GridID(ii)%idxHORIZ_ONLY
      write(iulog,*) '                   idxHgrid_d=',this%GridID(ii)%idxHgrid_d,' nvgrids=',this%GridID(ii)%nvgrids
      if(this%GridID(ii)%nvgrids.eq.0) then
        write(iulog,*) '                   idxVgrid_d= (none)'
      else
        write(iulog,*) '                   idxVgrid_d=',this%GridID(ii)%idxVgrid_d(:)
      endif
    end do
    write(iulog,*) '      -------------------------------------------'
    write(iulog,*) ' '
   
    ! Now print metadata for each of the Lists in this Link
    !-------------------------------------------------------
    call this%ListHgrid_m%printList
    call this%ListVgrid_m%printList
    call this%ListHgrid_d%printList
    call this%ListVgrid_d%printList
    call this%ListHmaps%printList
    call this%ListVmaps%printList
    call this%ListPaths%printList
    write(iulog,*) '   ============================================================'
    write(iulog,*) ' '

    return
  end subroutine LinkDataPaths_print
  !========================================================================



  !========================================================================
  !========================================================================
  ! Procedures For ListHgrid_m Structure
  !
  !========================================================================
  !========================================================================
  subroutine ListHgrid_m_reset(this)
    ! Reset the Iterator to point at the first Hgrid_m Link
    !-------------------------------------------------------
    class(ListHgrid_m_t):: this
    this%curr => this%first
    return
  end subroutine ListHgrid_m_reset
  !========================================================================
  subroutine ListHgrid_m_iterate(this)
    ! Increment the Iterator to point at the next Hgrid_m Link
    !---------------------------------------------------------
    class(ListHgrid_m_t):: this
    this%curr => this%curr%get_Next()
    return
  end subroutine ListHgrid_m_iterate
  !========================================================================
  function ListHgrid_m_isFirst(this)
    ! Check if the Iterator points to the end of the list.
    !----------------------------------------------------
    class(ListHgrid_m_t)        :: this
    logical                     :: ListHgrid_m_isFirst
    class(LinkHgrid_m_t),pointer:: lastPointer
    if(associated(this%curr)) then
      lastPointer => this%curr%get_Prev()
      ListHgrid_m_isFirst = .not.associated(lastPointer)
    else
      ListHgrid_m_isFirst = .true.
    endif
    return
  end function ListHgrid_m_isFirst
  !========================================================================
  function ListHgrid_m_isLast(this)
    ! Check if the Iterator points to the end of the list.
    !----------------------------------------------------
    class(ListHgrid_m_t)        :: this
    logical                     :: ListHgrid_m_isLast
    class(LinkHgrid_m_t),pointer:: nextPointer
    if(associated(this%curr)) then
      nextPointer => this%curr%get_Next()
      ListHgrid_m_isLast = .not.associated(nextPointer)
    else
      ListHgrid_m_isLast = .true.
    endif
    return
  end function ListHgrid_m_isLast
  !========================================================================
  function ListHgrid_m_addLink(this,Hgrid_m)
    ! Add a new Hgrid_m Dataset link to the end of the list.
    !========================================================
    class(ListHgrid_m_t)        :: this
    type(packed_hgrid_t),pointer:: Hgrid_m
    integer                     :: ListHgrid_m_addLink
    !
    ! Local Values
    !--------------
    class(LinkHgrid_m_t),pointer:: newLink

    ! Allocate and a new Hgrid_m Link and initialize 
    ! the pointer with Hgrid_m from the input.
    !----------------------------------------------------
    allocate(newLink)
    newLink%Hgrid_m => Hgrid_m

    ! Increment the Counter to assign a unique index 
    ! to this new Link
    !-----------------------------------------------
    this%Counter       = this%Counter + 1
    newLink%idxHgrid_m = this%Counter

    ! Add it to the list
    !--------------------
    if(.not.associated(this%last)) then
      ! If the current list is empty, initialize 
      ! the list with this Link
      !-------------------------------------------
      newLink%prev => null()
      newLink%next => null()
      this%first   => newLink
      this%last    => newLink
      this%curr    => newLink
    else
      ! Append this Link to the end of the current list.
      !----------------------------------------------------
      newLink%prev   => this%last
      newLink%next   => null()
      this%last%next => newLink
      this%last      => newLink
    endif

    ! Return the DataSet ID to the calling routine
    !-----------------------------------------------
    ListHgrid_m_addLink = newLink%idxHgrid_m

    return
  end function ListHgrid_m_addLink
  !========================================================================
  subroutine ListHgrid_m_rmLink(this,idxHgrid_m)
    ! Remove Hgrid_m Dataset link with given ID
    !========================================================
    class(ListHgrid_m_t)        :: this
    integer                     :: idxHgrid_m
    !
    ! Local Values
    !--------------
    class(LinkHgrid_m_t),pointer:: lastPointer
    class(LinkHgrid_m_t),pointer:: nextPointer
    type(packed_hgrid_t),pointer:: curr_Hgrid_m

    ! cycle thru the list to find the Link with matching ID
    !----------------------------------------------------------
    call this%reset
    do while(this%curr%get_ID() .ne. idxHgrid_m)
      call this%iterate
      if(.not.associated(this%curr)) then
        ! ERROR: idxHgrid_m not found in list
        !-------------------------------------
        write(iulog, *) 'ListHgrid_m ERROR: No ',idxHgrid_m,'index in list'
        call endrun('ERROR: ListHgrid_m Missing Grid Index')
      endif
    end do

    ! Now we want to remove the current 
    ! Link from the list
    !------------------------------------
    if(this%isFirst().and.this%isLast()) then
      ! We are removing the ONLY link in the List
      !--------------------------------------------
      curr_Hgrid_m => this%curr%get_Data()

      ! Deallocate curr Link
      !----------------------
      if(allocated(curr_Hgrid_m%Lon )) deallocate(curr_Hgrid_m%Lon )
      if(allocated(curr_Hgrid_m%Lat )) deallocate(curr_Hgrid_m%Lat )
      if(allocated(curr_Hgrid_m%Xmap)) deallocate(curr_Hgrid_m%Xmap)
      if(allocated(curr_Hgrid_m%Ymap)) deallocate(curr_Hgrid_m%Ymap)
      deallocate(curr_Hgrid_m)
      deallocate(this%curr)
   
      ! Nullify all pointers for an empty list
      !----------------------------------------
      this%first => null()
      this%curr  => null()
      this%last  => null()
    elseif(this%isFirst()) then
      ! We are removing the first Link in the list
      !-----------------------------------------------
      lastPointer => null()
      nextPointer => this%curr%get_Next()
      call nextPointer%set_Prev(lastPointer)

      ! Deallocate curr Link
      !----------------------
      curr_Hgrid_m => this%curr%get_Data()
      if(allocated(curr_Hgrid_m%Lon )) deallocate(curr_Hgrid_m%Lon )
      if(allocated(curr_Hgrid_m%Lat )) deallocate(curr_Hgrid_m%Lat )
      if(allocated(curr_Hgrid_m%Xmap)) deallocate(curr_Hgrid_m%Xmap)
      if(allocated(curr_Hgrid_m%Ymap)) deallocate(curr_Hgrid_m%Ymap)
      deallocate(curr_Hgrid_m)
      deallocate(this%curr)

      ! Reset first and curr pointers
      !-----------------------------------
      this%first => nextPointer
      this%curr  => nextPointer

      nullify(lastPointer)
      nullify(nextPointer)
    elseif(this%isLast()) then
      ! We are removing the last Link in the list
      !-----------------------------------------------
      lastPointer => this%curr%get_Prev()
      nextPointer => null()
      call lastPointer%set_Next(nextPointer)

      ! Deallocate curr
      !--------------------
      curr_Hgrid_m => this%curr%get_Data()
      if(allocated(curr_Hgrid_m%Lon )) deallocate(curr_Hgrid_m%Lon )
      if(allocated(curr_Hgrid_m%Lat )) deallocate(curr_Hgrid_m%Lat )
      if(allocated(curr_Hgrid_m%Xmap)) deallocate(curr_Hgrid_m%Xmap)
      if(allocated(curr_Hgrid_m%Ymap)) deallocate(curr_Hgrid_m%Ymap)
      deallocate(curr_Hgrid_m)
      deallocate(this%curr)

      ! Reset last and curr pointers
      !-----------------------------------
      this%last => lastPointer
      this%curr => lastPointer

      nullify(lastPointer)
      nullify(nextPointer)
    else
      ! We are removing a Link in the middle of the list
      !-----------------------------------------------------
      lastPointer => this%curr%get_Prev()
      nextPointer => this%curr%get_Next()
      call lastPointer%set_Next(nextPointer)
      call nextPointer%set_Prev(lastPointer)

      ! Deallocate curr
      !--------------------
      curr_Hgrid_m => this%curr%get_Data()
      if(allocated(curr_Hgrid_m%Lon )) deallocate(curr_Hgrid_m%Lon )
      if(allocated(curr_Hgrid_m%Lat )) deallocate(curr_Hgrid_m%Lat )
      if(allocated(curr_Hgrid_m%Xmap)) deallocate(curr_Hgrid_m%Xmap)
      if(allocated(curr_Hgrid_m%Ymap)) deallocate(curr_Hgrid_m%Ymap)
      deallocate(curr_Hgrid_m)
      deallocate(this%curr)

      ! Reset curr pointer
      !-----------------------------------
      this%curr => lastPointer

      nullify(lastPointer)
      nullify(nextPointer)
    endif

    return
  end subroutine ListHgrid_m_rmLink
  !========================================================================
  function ListHgrid_m_getLinkID(this)
    ! Return the ID idxHgrid_m for curr in List
    !---------------------------------------------
    class(ListHgrid_m_t):: this
    integer             :: ListHgrid_m_getLinkID
    ListHgrid_m_getLinkID = this%curr%get_ID()
    return
  end function ListHgrid_m_getLinkID
  !========================================================================
  function ListHgrid_m_getLinkData(this,idxHgrid_m)
    ! Given the Link ID, return a pointer to 
    ! the corresponding Hgrid_m data structure.
    !---------------------------------------------
    class(ListHgrid_m_t)        :: this
    integer                     :: idxHgrid_m
    type(packed_hgrid_t),pointer:: ListHgrid_m_getLinkData

    ! cycle thru the list to find the Link with matching ID
    !----------------------------------------------------------
    call this%reset
    do while(this%curr%get_ID() .ne. idxHgrid_m)
      call this%iterate
      if(.not.associated(this%curr)) then
        ! ERROR: idxHgrid_m not found in list
        !-------------------------------------
        write(iulog, *) 'ListHgrid_m ERROR: No ',idxHgrid_m,'index in list'
        call endrun('ERROR: ListHgrid_m Missing Grid Index')
      endif
    end do
    ListHgrid_m_getLinkData => this%curr%get_Data()
    return
  end function ListHgrid_m_getLinkData
  !========================================================================
  function ListHgrid_m_setLinkID(this,new_Hgrid_m,newLink)
    ! Return the ID idxHgrid_m for new_Hgrid_m. 
    ! If the grid exists on the list, then return the corresponding ID, if 
    ! this is a new grid, add it to the end of the list and return the new ID.
    ! Set newLink to indicate if this is a new or existing link
    !-------------------------------------------------------------------------
    class(ListHgrid_m_t)        :: this
    type(packed_hgrid_t),pointer:: new_Hgrid_m
    integer                     :: ListHgrid_m_setLinkID
    logical                     :: newLink
    ! 
    ! Local Values
    !-------------------
    logical:: link_found

    ! Cycle thru the list and see if there is already an  
    ! existing grid the matches the new Hgrid_m
    !----------------------------------------------------
    call this%reset
    link_found = .false.
    do while(associated(this%curr))
      if(this%curr%isEqual(new_Hgrid_m)) then
        link_found = .true.
        exit
      endif
      call this%iterate 
    end do

    if(link_found) then
      ! if the Grid exists, pass back the the corresponding ID. 
      !-----------------------------------------------------------------
      newLink = .false.
      ListHgrid_m_setLinkID = this%curr%get_ID()
    else
      ! Otherwise add the new_Hgrid_m to the end of the 
      ! list and return the new ID
      !--------------------------------------------------
      newLink = .true.
      ListHgrid_m_setLinkID = this%addLink(new_Hgrid_m)
    endif

    return
  end function ListHgrid_m_setLinkID
  !========================================================================
  subroutine ListHgrid_m_print(this)
    ! Cycle thru the current list and print 
    ! out the meta data for each link 
    !---------------------------------------------
    class(ListHgrid_m_t)        :: this
    !
    ! Local values
    !-------------
    logical:: empty_list

    ! Check to see if the list is empty
    !-----------------------------------
    call this%reset
    empty_list = .false.
    if(this%isFirst().and.this%isLast()) then
      if(.not.associated(this%curr)) then
        empty_list = .true.
      endif
    endif
    
    if(empty_list) then
      ! Print out "Empty List" message here
      !---------------------------------------
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    LIST of Registered Model Hgrid_m structures:'
      write(iulog,*) '          EMPTY LIST  '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) ' '
    else
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    LIST of Registered Model Hgrid_m structures:'
      write(iulog,*) ' '

      ! cycle thru the list
      !----------------------
      do while(associated(this%curr))
        ! Issue a print command for each link
        !-------------------------------------
        call this%curr%printLink
        call this%iterate
      end do
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) ' '
    endif

    return
  end subroutine ListHgrid_m_print
  !========================================================================



  !========================================================================
  !========================================================================
  ! Procedures For LinkHgrid_m Structure
  !
  !========================================================================
  !========================================================================
  function LinkHgrid_m_getID(this)
    ! Return the ID for this Hgrid_m Link
    !----------------------------------------------------
    class(LinkHgrid_m_t):: this
    integer             :: LinkHgrid_m_getID
    LinkHgrid_m_getID = this%idxHgrid_m
    return
  end function LinkHgrid_m_getID
  !========================================================================
  function LinkHgrid_m_getData(this)
    ! Return the pointer to the next Hgrid_m Link
    !----------------------------------------------------
    class(LinkHgrid_m_t)        :: this
    type(packed_hgrid_t),pointer:: LinkHgrid_m_getData
    LinkHgrid_m_getData => this%Hgrid_m
    return
  end function LinkHgrid_m_getData
  !========================================================================
  function LinkHgrid_m_getNext(this)
    ! Return the pointer to the next Hgrid_m Link
    !----------------------------------------------------
    class(LinkHgrid_m_t)        :: this
    class(LinkHgrid_m_t),pointer:: LinkHgrid_m_getNext
    LinkHgrid_m_getNext => this%next
    return
  end function LinkHgrid_m_getNext
  !========================================================================
  subroutine LinkHgrid_m_setNext(this,next)
    ! Set the pointer for the next Hgrid_m Link
    !----------------------------------------------------
    class(LinkHgrid_m_t)        :: this
    class(LinkHgrid_m_t),pointer:: next
    this%next => next
    return
  end subroutine LinkHgrid_m_setNext
  !========================================================================
  function LinkHgrid_m_getPrev(this)
    ! Return the pointer to the previous Hgrid_m Link
    !----------------------------------------------------
    class(LinkHgrid_m_t)        :: this
    class(LinkHgrid_m_t),pointer:: LinkHgrid_m_getPrev
    LinkHgrid_m_getPrev => this%prev
    return
  end function LinkHgrid_m_getPrev
  !========================================================================
  subroutine LinkHgrid_m_setPrev(this,prev)
    ! Return the pointer to the previous Hgrid_m Link
    !----------------------------------------------------
    class(LinkHgrid_m_t)        :: this
    class(LinkHgrid_m_t),pointer:: prev
    this%prev => prev
    return
  end subroutine LinkHgrid_m_setPrev
  !========================================================================
  function LinkHgrid_m_isEqual(this,new_Hgrid_m)
    ! Logical function which returns True if the given 
    ! Hgrid_m structure data matches those of the current Link.
    ! Other wise the structures differ and False is returned
    !----------------------------------------------------------
    class(LinkHgrid_m_t)        :: this
    type(packed_hgrid_t),pointer:: new_Hgrid_m
    logical                     :: LinkHgrid_m_isEqual
    !
    ! Local values
    !-----------------
    type(packed_hgrid_t),pointer:: My_Hgrid_m
    real(r8)                    :: diffSumR
    integer                     :: diffSumI
    integer                     :: ii

    ! Get a pointer to the current data structure
    !--------------------------------------------
    My_Hgrid_m => this%get_Data()

    ! Set the restult to true, then scan the structure for differences, 
    ! if/when a difference is found, set the result to false and return 
    ! immediately.
    !-------------------------------------------------------------------
    LinkHgrid_m_isEqual = .true.

    if(trim(My_Hgrid_m%name).ne.trim(new_Hgrid_m%name)) then
      LinkHgrid_m_isEqual = .false.
      return
    endif

    if(trim(My_Hgrid_m%model_gridname).ne.trim(new_Hgrid_m%model_gridname)) then
      LinkHgrid_m_isEqual = .false.
      return
    endif

    if(My_Hgrid_m%model_gridid.ne.new_Hgrid_m%model_gridid) then
      LinkHgrid_m_isEqual = .false.
      return
    endif

    if(My_Hgrid_m%ncol.ne.new_Hgrid_m%ncol) then
      LinkHgrid_m_isEqual = .false.
      return
    else
      diffSumR = 0._r8
      diffSumI = 0
      do ii=1,new_Hgrid_m%ncol
        diffSumR = diffSumR + abs(My_Hgrid_m%Lon(ii)-new_Hgrid_m%Lon(ii)) &
                            + abs(My_Hgrid_m%Lat(ii)-new_Hgrid_m%Lat(ii))
        diffSumI = diffSumI + abs(My_Hgrid_m%Xmap(ii)-new_Hgrid_m%Xmap(ii)) &
                            + abs(My_Hgrid_m%Ymap(ii)-new_Hgrid_m%Ymap(ii))
      end do
      if((diffSumR.gt.DIFF_TOLERANCE).or.(diffSumI.gt.0)) then
        LinkHgrid_m_isEqual = .false.
        return
      endif
    endif
    
    return
  end function LinkHgrid_m_isEqual
  !========================================================================
  subroutine LinkHgrid_m_print(this)
    ! Print out the Hgrid_m metadata for the current link
    !---------------------------------------------
    class(LinkHgrid_m_t):: this

    ! Print out identifying values
    !--------------------------------
    write(iulog,*) '      -------------------------------------------'
    write(iulog,*) '      idxHgrid_m=',this%idxHgrid_m
    write(iulog,*) '                     name= ',trim(this%Hgrid_m%name)
    write(iulog,*) '           model_gridname= ',trim(this%Hgrid_m%model_gridname)
    write(iulog,*) '                     Type= ',trim(this%Hgrid_m%Type)
    write(iulog,*) '                 PackType= ',trim(this%Hgrid_m%PackType)
    write(iulog,*) '                     ncol= ',     this%Hgrid_m%ncol
    write(iulog,*) '      -------------------------------------------'
    write(iulog,*) ' '

    return
  end subroutine LinkHgrid_m_print
  !========================================================================



  !========================================================================
  !========================================================================
  ! Procedures For ListVgrid_m Structure
  !
  !========================================================================
  !========================================================================
  subroutine ListVgrid_m_reset(this)
    ! Reset the Iterator to point at the first DataSet
    !----------------------------------------------------
    class(ListVgrid_m_t):: this
    this%curr => this%first
    return
  end subroutine ListVgrid_m_reset
  !========================================================================
  subroutine ListVgrid_m_iterate(this)
    ! Increment the Iterator to point at the next DataSet
    !----------------------------------------------------
    class(ListVgrid_m_t):: this
    this%curr => this%curr%get_Next()
    return
  end subroutine ListVgrid_m_iterate
  !========================================================================
  function ListVgrid_m_isFirst(this)
    ! Check if the Iterator points to the end of the list.
    !----------------------------------------------------
    class(ListVgrid_m_t)        :: this
    logical                     :: ListVgrid_m_isFirst
    class(LinkVgrid_m_t),pointer:: lastPointer
    if(associated(this%curr)) then
      lastPointer => this%curr%get_Prev()
      ListVgrid_m_isFirst = .not.associated(lastPointer)
    else
      ListVgrid_m_isFirst = .true.
    endif
    return
  end function ListVgrid_m_isFirst
  !========================================================================
  function ListVgrid_m_isLast(this)
    ! Check if the Iterator points to the end of the list.
    !----------------------------------------------------
    class(ListVgrid_m_t)        :: this
    logical                     :: ListVgrid_m_isLast
    class(LinkVgrid_m_t),pointer:: nextPointer
    if(associated(this%curr)) then
      nextPointer => this%curr%get_Next()
      ListVgrid_m_isLast = .not.associated(nextPointer)
    else
      ListVgrid_m_isLast = .true.
    endif
    return
  end function ListVgrid_m_isLast
  !========================================================================
  function ListVgrid_m_addLink(this,Vgrid_m)
    ! Add a new Vgrid_m Dataset link to the end of the list.
    !========================================================
    class(ListVgrid_m_t)        :: this
    type(vgrid_t)       ,pointer:: Vgrid_m
    integer                     :: ListVgrid_m_addLink
    !
    ! Local Values
    !--------------
    class(LinkVgrid_m_t),pointer:: newLink

    ! Allocate and a new Vgrid_m Link and initialize 
    ! the pointer with Vgrid_m from the input.
    !----------------------------------------------------
    allocate(newLink)
    newLink%Vgrid_m => Vgrid_m

    ! Increment the Counter to assign a unique index 
    ! to this new Link
    !-----------------------------------------------
    this%Counter       = this%Counter + 1
    newLink%idxVgrid_m = this%Counter

    ! Add it to the list
    !--------------------
    if(.not.associated(this%last)) then
      ! If the current list is empty, initialize 
      ! the list with this Link
      !-------------------------------------------
      newLink%prev => null()
      newLink%next => null()
      this%first   => newLink
      this%last    => newLink
      this%curr    => newLink
    else
      ! Append this Link to the end of the current list.
      !----------------------------------------------------
      newLink%prev   => this%last
      newLink%next   => null()
      this%last%next => newLink
      this%last      => newLink
    endif

    ! Return the DataSet ID to the calling routine
    !-----------------------------------------------
    ListVgrid_m_addLink = newLink%idxVgrid_m

    return
  end function ListVgrid_m_addLink
  !========================================================================
  subroutine ListVgrid_m_rmLink(this,idxVgrid_m)
    ! Remove Vgrid_m Dataset link with given ID
    !========================================================
    class(ListVgrid_m_t):: this
    integer             :: idxVgrid_m
    !
    ! Local Values
    !--------------
    class(LinkVgrid_m_t),pointer:: lastPointer
    class(LinkVgrid_m_t),pointer:: nextPointer
    type(vgrid_t)       ,pointer:: curr_Vgrid_m

    ! cycle thru the list to find the Link with matching ID
    !----------------------------------------------------------
    call this%reset
    do while(this%curr%get_ID() .ne. idxVgrid_m)
      call this%iterate
      if(.not.associated(this%curr)) then
        ! ERROR: idxVgrid_m not found in list
        !-------------------------------------
        write(iulog, *) 'ListVgrid_m ERROR: No ',idxVgrid_m,'index in list'
        call endrun('ERROR: ListVgrid_m Missing Grid Index')
      endif
    end do

    ! Now we want to remove the current 
    ! Link from the list
    !------------------------------------
    if(this%isFirst().and.this%isLast()) then
      ! We are removing the ONLY link in the List
      !--------------------------------------------
      curr_Vgrid_m => this%curr%get_Data()
   
      ! Deallocate curr Link
      !----------------------
      if(allocated(curr_Vgrid_m%grid_var)) deallocate(curr_Vgrid_m%grid_var)
      if(allocated(curr_Vgrid_m%req_var )) deallocate(curr_Vgrid_m%req_var )
      if(allocated(curr_Vgrid_m%HYA     )) deallocate(curr_Vgrid_m%HYA     )
      if(allocated(curr_Vgrid_m%HYB     )) deallocate(curr_Vgrid_m%HYB     )
      if(allocated(curr_Vgrid_m%PRES    )) deallocate(curr_Vgrid_m%PRES    )
      if(allocated(curr_Vgrid_m%HGHT    )) deallocate(curr_Vgrid_m%HGHT    )
      deallocate(curr_Vgrid_m)
      deallocate(this%curr)

      ! Nullify all pointers for an empty list
      !----------------------------------------
      this%first => null()
      this%curr  => null()
      this%last  => null()
    elseif(this%isFirst()) then
      ! We are removing the first Link in the list
      !-----------------------------------------------
      lastPointer => null()
      nextPointer => this%curr%get_Next()
      call nextPointer%set_Prev(lastPointer)

      ! Deallocate curr Link
      !----------------------
      curr_Vgrid_m => this%curr%get_Data()
      if(allocated(curr_Vgrid_m%grid_var)) deallocate(curr_Vgrid_m%grid_var)
      if(allocated(curr_Vgrid_m%req_var )) deallocate(curr_Vgrid_m%req_var )
      if(allocated(curr_Vgrid_m%HYA     )) deallocate(curr_Vgrid_m%HYA     )
      if(allocated(curr_Vgrid_m%HYB     )) deallocate(curr_Vgrid_m%HYB     )
      if(allocated(curr_Vgrid_m%PRES    )) deallocate(curr_Vgrid_m%PRES    )
      if(allocated(curr_Vgrid_m%HGHT    )) deallocate(curr_Vgrid_m%HGHT    )
      deallocate(curr_Vgrid_m)
      deallocate(this%curr)

      ! Reset fisrt and curr pointers
      !-----------------------------------
      this%first => nextPointer
      this%curr  => nextPointer

      nullify(lastPointer)
      nullify(nextPointer)
    elseif(this%isLast()) then
      ! We are removing the last Link in the list
      !-----------------------------------------------
      lastPointer => this%curr%get_Prev()
      nextPointer => null()
      call lastPointer%set_Next(nextPointer)

      ! Deallocate curr
      !--------------------
      curr_Vgrid_m => this%curr%get_Data()
      if(allocated(curr_Vgrid_m%grid_var)) deallocate(curr_Vgrid_m%grid_var)
      if(allocated(curr_Vgrid_m%req_var )) deallocate(curr_Vgrid_m%req_var )
      if(allocated(curr_Vgrid_m%HYA     )) deallocate(curr_Vgrid_m%HYA     )
      if(allocated(curr_Vgrid_m%HYB     )) deallocate(curr_Vgrid_m%HYB     )
      if(allocated(curr_Vgrid_m%PRES    )) deallocate(curr_Vgrid_m%PRES    )
      if(allocated(curr_Vgrid_m%HGHT    )) deallocate(curr_Vgrid_m%HGHT    )
      deallocate(curr_Vgrid_m)
      deallocate(this%curr)

      ! Reset last and curr pointers
      !-----------------------------------
      this%last => lastPointer
      this%curr => lastPointer

      nullify(lastPointer)
      nullify(nextPointer)
    else
      ! We are removing a Link in the middle of the list
      !-----------------------------------------------------
      lastPointer => this%curr%get_Prev()
      nextPointer => this%curr%get_Next()
      call lastPointer%set_Next(nextPointer)
      call nextPointer%set_Prev(lastPointer)

      ! Deallocate curr
      !--------------------
      curr_Vgrid_m => this%curr%get_Data()
      if(allocated(curr_Vgrid_m%grid_var)) deallocate(curr_Vgrid_m%grid_var)
      if(allocated(curr_Vgrid_m%req_var )) deallocate(curr_Vgrid_m%req_var )
      if(allocated(curr_Vgrid_m%HYA     )) deallocate(curr_Vgrid_m%HYA     )
      if(allocated(curr_Vgrid_m%HYB     )) deallocate(curr_Vgrid_m%HYB     )
      if(allocated(curr_Vgrid_m%PRES    )) deallocate(curr_Vgrid_m%PRES    )
      if(allocated(curr_Vgrid_m%HGHT    )) deallocate(curr_Vgrid_m%HGHT    )
      deallocate(curr_Vgrid_m)
      deallocate(this%curr)

      ! Reset curr pointer
      !-----------------------------------
      this%curr => lastPointer

      nullify(lastPointer)
      nullify(nextPointer)
    endif

    return
  end subroutine ListVgrid_m_rmLink
  !========================================================================
  function ListVgrid_m_getLinkID(this)
    ! Return the ID idxVgrid_m for curr in List
    !---------------------------------------------
    class(ListVgrid_m_t):: this
    integer             :: ListVgrid_m_getLinkID
    ListVgrid_m_getLinkID = this%curr%get_ID()
    return
  end function ListVgrid_m_getLinkID
  !========================================================================
  function ListVgrid_m_getLinkData(this,idxVgrid_m)
    ! Given the Link ID, return a pointer to 
    ! the corresponding Vgrid_m data structure.
    !---------------------------------------------
    class(ListVgrid_m_t) :: this
    integer              :: idxVgrid_m
    type(vgrid_t),pointer:: ListVgrid_m_getLinkData

    ! cycle thru the list to find the Link with matching ID
    !----------------------------------------------------------
    call this%reset
    do while(this%curr%get_ID() .ne. idxVgrid_m)
      call this%iterate
      if(.not.associated(this%curr)) then
        ! ERROR: idxVgrid_m not found in list
        !-------------------------------------
        write(iulog, *) 'ListVgrid_m ERROR: No ',idxVgrid_m,'index in list'
        call endrun('ERROR: ListVgrid_m Missing Grid Index')
      endif
    end do
    ListVgrid_m_getLinkData => this%curr%get_Data()
    return
  end function ListVgrid_m_getLinkData
  !========================================================================
  function ListVgrid_m_setLinkID(this,new_Vgrid_m,newLink)
    ! Return the ID idxVgrid_m for new_Vgrid_m. 
    ! If the grid exists on the list, then return the corresponding ID, if 
    ! this is a new grid, add it to the end of the list and return the new ID.
    ! Set newLink to indicate if this is a new or existing link
    !-------------------------------------------------------------------------
    class(ListVgrid_m_t) :: this
    type(vgrid_t),pointer:: new_Vgrid_m
    integer              :: ListVgrid_m_setLinkID
    logical              :: newLink
    ! 
    ! Local Values
    !-------------------
    logical:: link_found

    ! Cycle thru the list and see if there is already an  
    ! existing grid the matches the new Vgrid_m
    !----------------------------------------------------
    call this%reset
    link_found = .false.
    do while(associated(this%curr))
      if(this%curr%isEqual(new_Vgrid_m)) then
        link_found = .true.
        exit
      endif
      call this%iterate 
    end do

    if(link_found) then
      ! if the Grid exists, pass back the the corresponding ID. 
      !-----------------------------------------------------------------
      newLink = .false.
      ListVgrid_m_setLinkID = this%curr%get_ID()
    else
      ! Otherwise add the new_Vgrid_m to the end of the 
      ! list and return the new ID
      !--------------------------------------------------
      newLink = .true.
      ListVgrid_m_setLinkID = this%addLink(new_Vgrid_m)
    endif

    return
  end function ListVgrid_m_setLinkID
  !========================================================================
  subroutine ListVgrid_m_print(this)
    ! Cycle thru the current list and print 
    ! out the meta data for each link 
    !---------------------------------------------
    class(ListVgrid_m_t)        :: this
    !
    ! Local values
    !-------------
    logical:: empty_list

    ! Check to see if the list is empty
    !-----------------------------------
    call this%reset
    empty_list = .false.
    if(this%isFirst().and.this%isLast()) then
      if(.not.associated(this%curr)) then
        empty_list = .true.
      endif
    endif
    
    if(empty_list) then
      ! Print out "Empty List" message here
      !---------------------------------------
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    LIST of Registered Model Vgrid_m structures:'
      write(iulog,*) '          EMPTY LIST  '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) ' '
    else
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    LIST of Registered Model Vgrid_m structures:'
      write(iulog,*) ' '

      ! cycle thru the list
      !----------------------
      do while(associated(this%curr))
        ! Issue a print command for each link
        !-------------------------------------
        call this%curr%printLink
        call this%iterate
      end do
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) ' '
    endif

    return
  end subroutine ListVgrid_m_print
  !========================================================================



  !========================================================================
  !========================================================================
  ! Procedures For LinkVgrid_m Structure
  !
  !========================================================================
  !========================================================================
  function LinkVgrid_m_getID(this)
    ! Return the ID for this Hgrid_m Link
    !----------------------------------------------------
    class(LinkVgrid_m_t):: this
    integer             :: LinkVgrid_m_getID
    LinkVgrid_m_getID = this%idxVgrid_m
    return
  end function LinkVgrid_m_getID
  !========================================================================
  function LinkVgrid_m_getData(this)
    ! Return the pointer to the next Vgrid_m Link
    !----------------------------------------------------
    class(LinkVgrid_m_t)        :: this
    type(vgrid_t)       ,pointer:: LinkVgrid_m_getData
    LinkVgrid_m_getData => this%Vgrid_m
    return
  end function LinkVgrid_m_getData
  !========================================================================
  function LinkVgrid_m_getNext(this)
    ! Return the pointer to the next Vgrid_m Link
    !----------------------------------------------------
    class(LinkVgrid_m_t)        :: this
    class(LinkVgrid_m_t),pointer:: LinkVgrid_m_getNext
    LinkVgrid_m_getNext => this%next
    return
  end function LinkVgrid_m_getNext
  !========================================================================
  subroutine LinkVgrid_m_setNext(this,next)
    ! Set the pointer for the next Vgrid_m Link
    !----------------------------------------------------
    class(LinkVgrid_m_t)        :: this
    class(LinkVgrid_m_t),pointer:: next
    this%next => next
    return
  end subroutine LinkVgrid_m_setNext
  !========================================================================
  function LinkVgrid_m_getPrev(this)
    ! Return the pointer to the previous Vgrid_m Link
    !----------------------------------------------------
    class(LinkVgrid_m_t)        :: this
    class(LinkVgrid_m_t),pointer:: LinkVgrid_m_getPrev
    LinkVgrid_m_getPrev => this%prev
    return
  end function LinkVgrid_m_getPrev
  !========================================================================
  subroutine LinkVgrid_m_setPrev(this,prev)
    ! Return the pointer to the previous Vgrid_m Link
    !----------------------------------------------------
    class(LinkVgrid_m_t)        :: this
    class(LinkVgrid_m_t),pointer:: prev
    this%prev => prev
    return
  end subroutine LinkVgrid_m_setPrev
  !========================================================================
  function LinkVgrid_m_isEqual(this,new_Vgrid_m)
    ! Logical function which returns True if the given 
    ! Vgrid_m structure data matches those of the current Link.
    ! Other wise the structures differ and False is returned
    !----------------------------------------------------
    class(LinkVgrid_m_t) :: this
    type(vgrid_t),pointer:: new_Vgrid_m
    logical              :: LinkVgrid_m_isEqual
    !
    ! Local values
    !-----------------
    type(vgrid_t),pointer:: My_Vgrid_m
    real(r8)             :: diffSum
    integer              :: ii

    ! Get a pointer to the current data structure
    !--------------------------------------------
    My_Vgrid_m => this%get_Data()

    ! Set the restult to true, then scan the structure for differences, 
    ! if/when a difference is found, set the result to false and return 
    ! immediately.
    !-------------------------------------------------------------------
    LinkVgrid_m_isEqual = .true.

    if(trim(My_Vgrid_m%name).ne.trim(new_Vgrid_m%name)) then
      LinkVgrid_m_isEqual = .false.
      return
    endif

    do ii=1,My_Vgrid_m%ngrid_var
      if(trim(My_Vgrid_m%grid_var(ii)).ne.trim(new_Vgrid_m%grid_var(ii))) then
        LinkVgrid_m_isEqual = .false.
        return
      endif
    end do

    if(My_Vgrid_m%nlev.ne.new_Vgrid_m%nlev) then
      LinkVgrid_m_isEqual = .false.
      return
    endif
    
    return
  end function LinkVgrid_m_isEqual
  !========================================================================
  subroutine LinkVgrid_m_print(this)
    ! Print out the Vgrid_m metadata for the current link
    !---------------------------------------------
    class(LinkVgrid_m_t):: this
    !
    ! Local values
    !---------------
    integer:: ii

    ! Print out identifying values
    !--------------------------------
    write(iulog,*) '      -------------------------------------------'
    write(iulog,*) '      idxVgrid_m=',this%idxVgrid_m
    write(iulog,*) '                     name= ',trim(this%Vgrid_m%name)
    write(iulog,*) '                     type= ',trim(this%Vgrid_m%Type)
    write(iulog,*) '                     nlev= ',this%Vgrid_m%nlev
    write(iulog,*) '                ngrid_var= ',this%Vgrid_m%ngrid_var
    do ii=1,this%Vgrid_m%ngrid_var
      write(iulog,*) '             Grid Variable',ii,                           &
                                         ' name=',trim(this%Vgrid_m%grid_var(ii))
    end do
    write(iulog,*) '                 nreq_var= ',this%Vgrid_m%nreq_var
    do ii=1,this%Vgrid_m%nreq_var
      write(iulog,*) '           Required Variable',ii,                          &
                                           ' name=',trim(this%Vgrid_m%req_var(ii))
    end do
    write(iulog,*) '                      Pref = ',this%Vgrid_m%Pref
    if(allocated(this%Vgrid_m%HYA)) then
      write(iulog,*) '                     HYA = ALLOCATED'
    else
      write(iulog,*) '                     HYA = NOT ALLOCATED'
    endif
    if(allocated(this%Vgrid_m%HYB)) then
      write(iulog,*) '                     HYB = ALLOCATED'
    else
      write(iulog,*) '                     HYB = NOT ALLOCATED'
    endif
    if(allocated(this%Vgrid_m%PRES)) then
      write(iulog,*) '                    PRES = ALLOCATED'
    else
      write(iulog,*) '                    PRES = NOT ALLOCATED'
    endif
    if(allocated(this%Vgrid_m%HGHT)) then
      write(iulog,*) '                    HGHT = ALLOCATED'
    else
      write(iulog,*) '                    HGHT = NOT ALLOCATED'
    endif
    write(iulog,*) '      -------------------------------------------'
    write(iulog,*) ' '

    return
  end subroutine LinkVgrid_m_print
  !========================================================================



  !========================================================================
  !========================================================================
  ! Procedures For ListHgrid_d Structure
  !
  !========================================================================
  !========================================================================
  subroutine ListHgrid_d_reset(this)
    ! Reset the Iterator to point at the first Hgrid_d Link
    !-------------------------------------------------------
    class(ListHgrid_d_t):: this
    this%curr => this%first
    return
  end subroutine ListHgrid_d_reset
  !========================================================================
  subroutine ListHgrid_d_iterate(this)
    ! Increment the Iterator to point at the next Hgrid_d Link
    !----------------------------------------------------------
    class(ListHgrid_d_t):: this
    this%curr => this%curr%get_Next()
    return
  end subroutine ListHgrid_d_iterate
  !========================================================================
  function ListHgrid_d_isFirst(this)
    ! Check if the Iterator points to the begining of the list.
    !----------------------------------------------------------
    class(ListHgrid_d_t)        :: this
    logical                     :: ListHgrid_d_isFirst
    class(LinkHgrid_d_t),pointer:: lastPointer
    if(associated(this%curr)) then
      lastPointer => this%curr%get_Prev()
      ListHgrid_d_isFirst = .not.associated(lastPointer)
    else
      ListHgrid_d_isFirst = .true.
    endif
    return
  end function ListHgrid_d_isFirst
  !========================================================================
  function ListHgrid_d_isLast(this)
    ! Check if the Iterator points to the end of the list.
    !----------------------------------------------------
    class(ListHgrid_d_t)        :: this
    logical                     :: ListHgrid_d_isLast
    class(LinkHgrid_d_t),pointer:: nextPointer
    if(associated(this%curr)) then
      nextPointer => this%curr%get_Next()
      ListHgrid_d_isLast = .not.associated(nextPointer)
    else
      ListHgrid_d_isLast = .true.
    endif
    return
  end function ListHgrid_d_isLast
  !========================================================================
  function ListHgrid_d_addLink(this,Hgrid_d)
    ! Add a new Hgrid_d Dataset link to the end of the list.
    !========================================================
    class(ListHgrid_d_t)         :: this
    type(dataset_hgrid_t),pointer:: Hgrid_d
    integer                      :: ListHgrid_d_addLink
    !
    ! Local Values
    !--------------
    class(LinkHgrid_d_t),pointer:: newLink

    ! Allocate and a new Hgrid_d Link and initialize 
    ! the pointer with Hgrid_d from the input.
    !----------------------------------------------------
    allocate(newLink)
    newLink%Hgrid_d => Hgrid_d

    ! Increment the Counter to assign a unique index 
    ! to this new Link
    !-----------------------------------------------
    this%Counter       = this%Counter + 1
    newLink%idxHgrid_d = this%Counter

    ! Add it to the list
    !--------------------
    if(.not.associated(this%last)) then
      ! If the current list is empty, initialize 
      ! the list with this Link
      !-------------------------------------------
      newLink%prev => null()
      newLink%next => null()
      this%first   => newLink
      this%last    => newLink
      this%curr    => newLink
    else
      ! Append this Link to the end of the current list.
      !----------------------------------------------------
      newLink%prev   => this%last
      newLink%next   => null()
      this%last%next => newLink
      this%last      => newLink
    endif

    ! Return the DataSet ID to the calling routine
    !-----------------------------------------------
    ListHgrid_d_addLink = newLink%idxHgrid_d

    return
  end function ListHgrid_d_addLink
  !========================================================================
  subroutine ListHgrid_d_rmLink(this,idxHgrid_d)
    ! Remove Hgrid_d Dataset link with given ID
    !========================================================
    class(ListHgrid_d_t)        :: this
    integer                     :: idxHgrid_d
    !
    ! Local Values
    !--------------
    class(LinkHgrid_d_t) ,pointer:: lastPointer
    class(LinkHgrid_d_t) ,pointer:: nextPointer
    type(dataset_hgrid_t),pointer:: curr_Hgrid_d

    ! cycle thru the list to find the Link with matching ID
    !----------------------------------------------------------
    call this%reset
    do while(this%curr%get_ID() .ne. idxHgrid_d)
      call this%iterate
      if(.not.associated(this%curr)) then
        ! ERROR: idxHgrid_d not found in list
        !-------------------------------------
        write(iulog, *) 'ListHgrid_d ERROR: No ',idxHgrid_d,'index in list'
        call endrun('ERROR: ListHgrid_d Missing Grid Index')
      endif
    end do

    ! Now we want to remove the current 
    ! Link from the list
    !------------------------------------
    if(this%isFirst().and.this%isLast()) then
      ! We are removing the ONLY link in the List
      !--------------------------------------------
      curr_Hgrid_d => this%curr%get_Data()

      ! Deallocate curr Link
      !----------------------
      if(allocated(curr_Hgrid_d%Lon )) deallocate(curr_Hgrid_d%Lon )
      if(allocated(curr_Hgrid_d%Lat )) deallocate(curr_Hgrid_d%Lat )
      deallocate(curr_Hgrid_d)
      deallocate(this%curr)

      ! Nullify all pointers for an empty list
      !----------------------------------------
      this%first => null()
      this%curr  => null()
      this%last  => null()
    elseif(this%isFirst()) then
      ! We are removing the first Link in the list
      !-----------------------------------------------
      lastPointer => null()
      nextPointer => this%curr%get_Next()
      call nextPointer%set_Prev(lastPointer)

      ! Deallocate curr Link
      !----------------------
      curr_Hgrid_d => this%curr%get_Data()
      if(allocated(curr_Hgrid_d%Lon )) deallocate(curr_Hgrid_d%Lon )
      if(allocated(curr_Hgrid_d%Lat )) deallocate(curr_Hgrid_d%Lat )
      deallocate(curr_Hgrid_d)
      deallocate(this%curr)

      ! Reset fisrt and curr pointers
      !-----------------------------------
      this%first => nextPointer
      this%curr  => nextPointer

      nullify(lastPointer)
      nullify(nextPointer)
    elseif(this%isLast()) then
      ! We are removing the last Link in the list
      !-----------------------------------------------
      lastPointer => this%curr%get_Prev()
      nextPointer => null()
      call lastPointer%set_Next(nextPointer)

      ! Deallocate curr
      !--------------------
      curr_Hgrid_d => this%curr%get_Data()
      if(allocated(curr_Hgrid_d%Lon )) deallocate(curr_Hgrid_d%Lon )
      if(allocated(curr_Hgrid_d%Lat )) deallocate(curr_Hgrid_d%Lat )
      deallocate(curr_Hgrid_d)
      deallocate(this%curr)

      ! Reset last and curr pointers
      !-----------------------------------
      this%last => lastPointer
      this%curr => lastPointer

      nullify(lastPointer)
      nullify(nextPointer)
    else
      ! We are removing a Link in the middle of the list
      !-----------------------------------------------------
      lastPointer => this%curr%get_Prev()
      nextPointer => this%curr%get_Next()
      call lastPointer%set_Next(nextPointer)
      call nextPointer%set_Prev(lastPointer)

      ! Deallocate curr
      !--------------------
      curr_Hgrid_d => this%curr%get_Data()
      if(allocated(curr_Hgrid_d%Lon )) deallocate(curr_Hgrid_d%Lon )
      if(allocated(curr_Hgrid_d%Lat )) deallocate(curr_Hgrid_d%Lat )
      deallocate(curr_Hgrid_d)
      deallocate(this%curr)

      ! Reset curr pointer
      !-----------------------------------
      this%curr => lastPointer

      nullify(lastPointer)
      nullify(nextPointer)
    endif

    return
  end subroutine ListHgrid_d_rmLink
  !========================================================================
  function ListHgrid_d_getLinkID(this)
    ! Return the ID idxHgrid_d for curr in List
    !---------------------------------------------
    class(ListHgrid_d_t):: this
    integer             :: ListHgrid_d_getLinkID
    ListHgrid_d_getLinkID = this%curr%get_ID()
    return
  end function ListHgrid_d_getLinkID
  !========================================================================
  function ListHgrid_d_getLinkData(this,idxHgrid_d)
    ! Given the Link ID, return a pointer to 
    ! the corresponding Hgrid_d data structure.
    !---------------------------------------------
    class(ListHgrid_d_t)         :: this
    integer                      :: idxHgrid_d
    type(dataset_hgrid_t),pointer:: ListHgrid_d_getLinkData

    ! cycle thru the list to find the Link with matching ID
    !----------------------------------------------------------
    call this%reset
    do while(this%curr%get_ID() .ne. idxHgrid_d)
      call this%iterate
      if(.not.associated(this%curr)) then
        ! ERROR: idxHgrid_d not found in list
        !-------------------------------------
        write(iulog, *) 'ListHgrid_d ERROR: No ',idxHgrid_d,'index in list'
        call endrun('ERROR: ListHgrid_d Missing Grid Index')
      endif
    end do
    ListHgrid_d_getLinkData => this%curr%get_Data()
    return
  end function ListHgrid_d_getLinkData
  !========================================================================
  function ListHgrid_d_setLinkID(this,new_Hgrid_d,newLink)
    ! Return the ID idxHgrid_d for new_Hgrid_d. 
    ! If the grid exists on the list, then return the corresponding ID, if 
    ! this is a new grid, add it to the end of the list and return the new ID.
    ! Set newLink to indicate if this is a new or existing link
    !-------------------------------------------------------------------------
    class(ListHgrid_d_t)         :: this
    type(dataset_hgrid_t),pointer:: new_Hgrid_d
    integer                      :: ListHgrid_d_setLinkID
    logical                      :: newLink
    ! 
    ! Local Values
    !-------------------
    logical:: link_found

    ! Cycle thru the list and see if there is already an  
    ! existing grid the matches the new Hgrid_d
    !----------------------------------------------------
    call this%reset
    link_found = .false.
    do while(associated(this%curr))
      if(this%curr%isEqual(new_Hgrid_d)) then
        link_found = .true.
        exit
      endif
      call this%iterate 
    end do

    if(link_found) then
      ! if the Grid exists, pass back the the corresponding ID. 
      !-----------------------------------------------------------------
      newLink = .false.
      ListHgrid_d_setLinkID = this%curr%get_ID()
    else
      ! Otherwise add the new_Hgrid_d to the end of the 
      ! list and return the new ID
      !--------------------------------------------------
      newLink = .true.
      ListHgrid_d_setLinkID = this%addLink(new_Hgrid_d)
    endif

    return
  end function ListHgrid_d_setLinkID
  !========================================================================
  subroutine ListHgrid_d_print(this)
    ! Cycle thru the current list and print 
    ! out the meta data for each link 
    !---------------------------------------------
    class(ListHgrid_d_t)        :: this
    !
    ! Local values
    !-------------
    logical:: empty_list

    ! Check to see if the list is empty
    !-----------------------------------
    call this%reset
    empty_list = .false.
    if(this%isFirst().and.this%isLast()) then
      if(.not.associated(this%curr)) then
        empty_list = .true.
      endif
    endif
    
    if(empty_list) then
      ! Print out "Empty List" message here
      !---------------------------------------
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    LIST of Registered DataSet Hgrid_d structures:'
      write(iulog,*) '          EMPTY LIST  '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) ' '
    else
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    LIST of Registered DataSet Hgrid_d structures:'
      write(iulog,*) ' '

      ! cycle thru the list
      !----------------------
      do while(associated(this%curr))
        ! Issue a print command for each link
        !-------------------------------------
        call this%curr%printLink
        call this%iterate
      end do
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) ' '
    endif

    return
  end subroutine ListHgrid_d_print
  !========================================================================



  !========================================================================
  !========================================================================
  ! Procedures For LinkHgrid_d Structure
  !
  !========================================================================
  !========================================================================
  function LinkHgrid_d_getID(this)
    ! Return the ID for this Hgrid_d Link
    !----------------------------------------------------
    class(LinkHgrid_d_t):: this
    integer             :: LinkHgrid_d_getID
    LinkHgrid_d_getID = this%idxHgrid_d
    return
  end function LinkHgrid_d_getID
  !========================================================================
  function LinkHgrid_d_getData(this)
    ! Return the pointer to the next Hgrid_d Link
    !----------------------------------------------------
    class(LinkHgrid_d_t)         :: this
    type(dataset_hgrid_t),pointer:: LinkHgrid_d_getData
    LinkHgrid_d_getData => this%Hgrid_d
    return
  end function LinkHgrid_d_getData
  !========================================================================
  function LinkHgrid_d_getNext(this)
    ! Return the pointer to the next Hgrid_d Link
    !----------------------------------------------------
    class(LinkHgrid_d_t)        :: this
    class(LinkHgrid_d_t),pointer:: LinkHgrid_d_getNext
    LinkHgrid_d_getNext => this%next
    return
  end function LinkHgrid_d_getNext
  !========================================================================
  subroutine LinkHgrid_d_setNext(this,next)
    ! Set the pointer for the next Hgrid_d Link
    !----------------------------------------------------
    class(LinkHgrid_d_t)        :: this
    class(LinkHgrid_d_t),pointer:: next
    this%next => next
    return
  end subroutine LinkHgrid_d_setNext
  !========================================================================
  function LinkHgrid_d_getPrev(this)
    ! Return the pointer to the previous Hgrid_d Link
    !----------------------------------------------------
    class(LinkHgrid_d_t)        :: this
    class(LinkHgrid_d_t),pointer:: LinkHgrid_d_getPrev
    LinkHgrid_d_getPrev => this%prev
    return
  end function LinkHgrid_d_getPrev
  !========================================================================
  subroutine LinkHgrid_d_setPrev(this,prev)
    ! Return the pointer to the previous Hgrid_d Link
    !----------------------------------------------------
    class(LinkHgrid_d_t)        :: this
    class(LinkHgrid_d_t),pointer:: prev
    this%prev => prev
    return
  end subroutine LinkHgrid_d_setPrev
  !========================================================================
  function LinkHgrid_d_isEqual(this,new_Hgrid_d)
    ! Logical function which returns True if the given 
    ! Hgrid_d structure data matches those of the current Link.
    ! Other wise the structures differ and False is returned
    !----------------------------------------------------
    class(LinkHgrid_d_t):: this
    type(dataset_hgrid_t),pointer:: new_Hgrid_d
    logical             :: LinkHgrid_d_isEqual
    !
    ! Local values
    !-----------------
    type(dataset_hgrid_t),pointer:: My_Hgrid_d
    real(r8)                     :: diffSum
    integer                      :: ii

    ! Get a pointer to the current data structure
    !--------------------------------------------
    My_Hgrid_d => this%get_Data()

    ! Set the restult to true, then scan the structure for differences, 
    ! if/when a difference is found, set the result to false and return 
    ! immediately.
    !-------------------------------------------------------------------
    LinkHgrid_d_isEqual = .true.

    if(trim(My_Hgrid_d%name).ne.trim(new_Hgrid_d%name)) then
      LinkHgrid_d_isEqual = .false.
      return
    endif

    if(My_Hgrid_d%datasetID.ne.new_Hgrid_d%datasetID) then
      LinkHgrid_d_isEqual = .false.
      return
    endif

    do ii=1,My_Hgrid_d%nhdim
      if(trim(My_Hgrid_d%dimName(ii)).ne.trim(new_Hgrid_d%dimName(ii))) then
        LinkHgrid_d_isEqual = .false.
        return
      endif
      if(My_Hgrid_d%dimID(ii).ne.new_Hgrid_d%dimID(ii)) then
        LinkHgrid_d_isEqual = .false.
        return
      endif
    end do

    if(My_Hgrid_d%nlon.ne.new_Hgrid_d%nlon) then
      LinkHgrid_d_isEqual = .false.
      return
    else
      diffSum = 0._r8
      do ii=1,new_Hgrid_d%nlon
        diffSum = diffSum + abs(My_Hgrid_d%Lon(ii)-new_Hgrid_d%Lon(ii))
      end do
      if(diffSum.gt.DIFF_TOLERANCE) then
        LinkHgrid_d_isEqual = .false.
        return
      endif
    endif
    
    if(My_Hgrid_d%nlat.ne.new_Hgrid_d%nlat) then
      LinkHgrid_d_isEqual = .false.
      return
    else
      diffSum = 0._r8
      do ii=1,new_Hgrid_d%nlat
        diffSum = diffSum + abs(My_Hgrid_d%Lat(ii)-new_Hgrid_d%Lat(ii))
      end do
      if(diffSum.gt.DIFF_TOLERANCE) then
        LinkHgrid_d_isEqual = .false.
        return
      endif
    endif
    
    return
  end function LinkHgrid_d_isEqual
  !========================================================================
  subroutine LinkHgrid_d_print(this)
    ! Print out the Hgrid_m metadata for the current link
    !---------------------------------------------
    class(LinkHgrid_d_t):: this
    !
    ! Local Values
    !-------------
    integer:: ii

    ! Print out identifying values
    !--------------------------------
    write(iulog,*) '      -------------------------------------------'
    write(iulog,*) '      idxHgrid_d=',this%idxHgrid_d
    write(iulog,*) '                     name= ',trim(this%Hgrid_d%name)
    write(iulog,*) '                     type= ',trim(this%Hgrid_d%Type)
    write(iulog,*) '                    nhdim= ',this%Hgrid_d%nhdim
    do ii=1,this%Hgrid_d%nhdim
      write(iulog,*) '                     Dim =',ii,' Name= ',trim(this%Hgrid_d%dimName(ii))
    end do
    write(iulog,*) '      -------------------------------------------'
    write(iulog,*) ' '

    return
  end subroutine LinkHgrid_d_print
  !========================================================================



  !========================================================================
  !========================================================================
  ! Procedures For ListVgrid_d Structure
  !
  !========================================================================
  !========================================================================
  subroutine ListVgrid_d_reset(this)
    ! Reset the Iterator to point at the first Vgrid_d Link
    !-------------------------------------------------------
    class(ListVgrid_d_t):: this
    this%curr => this%first
    return
  end subroutine ListVgrid_d_reset
  !========================================================================
  subroutine ListVgrid_d_iterate(this)
    ! Increment the Iterator to point at the next Vgrid_d Link
    !----------------------------------------------------------
    class(ListVgrid_d_t):: this
    this%curr => this%curr%get_Next()
    return
  end subroutine ListVgrid_d_iterate
  !========================================================================
  function ListVgrid_d_isFirst(this)
    ! Check if the Iterator points to the end of the list.
    !----------------------------------------------------
    class(ListVgrid_d_t)        :: this
    logical                     :: ListVgrid_d_isFirst
    class(LinkVgrid_d_t),pointer:: lastPointer
    if(associated(this%curr)) then
      lastPointer => this%curr%get_Prev()
      ListVgrid_d_isFirst = .not.associated(lastPointer)
    else
      ListVgrid_d_isFirst = .true.
    endif
    return
  end function ListVgrid_d_isFirst
  !========================================================================
  function ListVgrid_d_isLast(this)
    ! Check if the Iterator points to the end of the list.
    !----------------------------------------------------
    class(ListVgrid_d_t)        :: this
    logical                     :: ListVgrid_d_isLast
    class(LinkVgrid_d_t),pointer:: nextPointer
    if(associated(this%curr)) then
      nextPointer => this%curr%get_Next()
      ListVgrid_d_isLast = .not.associated(nextPointer)
    else
      ListVgrid_d_isLast = .true.
    endif
    return
  end function ListVgrid_d_isLast
  !========================================================================
  function ListVgrid_d_addLink(this,Vgrid_d)
    ! Add a new Vgrid_d Dataset link to the end of the list.
    !========================================================
    class(ListVgrid_d_t)        :: this
    type(vgrid_t)       ,pointer:: Vgrid_d
    integer                     :: ListVgrid_d_addLink
    !
    ! Local Values
    !--------------
    class(LinkVgrid_d_t),pointer:: newLink

    ! Allocate and a new Vgrid_d Link and initialize 
    ! the pointer with Vgrid_d from the input.
    !----------------------------------------------------
    allocate(newLink)
    newLink%Vgrid_d => Vgrid_d

    ! Increment the Counter to assign a unique index 
    ! to this new Link
    !-----------------------------------------------
    this%Counter       = this%Counter + 1
    newLink%idxVgrid_d = this%Counter

    ! Add it to the list
    !--------------------
    if(.not.associated(this%last)) then
      ! If the current list is empty, initialize 
      ! the list with this Link
      !-------------------------------------------
      newLink%prev => null()
      newLink%next => null()
      this%first   => newLink
      this%last    => newLink
      this%curr    => newLink
    else
      ! Append this Link to the end of the current list.
      !----------------------------------------------------
      newLink%prev   => this%last
      newLink%next   => null()
      this%last%next => newLink
      this%last      => newLink
    endif

    ! Return the DataSet ID to the calling routine
    !-----------------------------------------------
    ListVgrid_d_addLink = newLink%idxVgrid_d

    return
  end function ListVgrid_d_addLink
  !========================================================================
  subroutine ListVgrid_d_rmLink(this,idxVgrid_d)
    ! Remove Vgrid_d Dataset link with given ID
    !========================================================
    class(ListVgrid_d_t):: this
    integer             :: idxVgrid_d
    !
    ! Local Values
    !--------------
    class(LinkVgrid_d_t),pointer:: lastPointer
    class(LinkVgrid_d_t),pointer:: nextPointer
    type(vgrid_t)       ,pointer:: curr_Vgrid_d

    ! cycle thru the list to find the Link with matching ID
    !----------------------------------------------------------
    call this%reset
    do while(this%curr%get_ID() .ne. idxVgrid_d)
      call this%iterate
      if(.not.associated(this%curr)) then
        ! ERROR: idxVgrid_d not found in list
        !-------------------------------------
        write(iulog, *) 'ListVgrid_d ERROR: No ',idxVgrid_d,'index in list'
        call endrun('ERROR: ListVgrid_d Missing Grid Index')
      endif
    end do

    ! Now we want to remove the current 
    ! Link from the list
    !------------------------------------
    if(this%isFirst().and.this%isLast()) then
      ! We are removing the ONLY link in the List
      !--------------------------------------------
      curr_Vgrid_d => this%curr%get_Data()
      
      ! Deallocate curr Link
      !----------------------
      if(allocated(curr_Vgrid_d%grid_var)) deallocate(curr_Vgrid_d%grid_var)
      if(allocated(curr_Vgrid_d%req_var )) deallocate(curr_Vgrid_d%req_var )
      if(allocated(curr_Vgrid_d%HYA     )) deallocate(curr_Vgrid_d%HYA     )
      if(allocated(curr_Vgrid_d%HYB     )) deallocate(curr_Vgrid_d%HYB     )
      if(allocated(curr_Vgrid_d%PRES    )) deallocate(curr_Vgrid_d%PRES    )
      if(allocated(curr_Vgrid_d%HGHT    )) deallocate(curr_Vgrid_d%HGHT    )
      deallocate(curr_Vgrid_d)
      deallocate(this%curr)
   
      ! Nullify all pointers for an empty list
      !----------------------------------------
      this%first => null()
      this%curr  => null()
      this%last  => null()
    elseif(this%isFirst()) then
      ! We are removing the first Link in the list
      !-----------------------------------------------
      lastPointer => null()
      nextPointer => this%curr%get_Next()
      call nextPointer%set_Prev(lastPointer)

      ! Deallocate curr Link
      !----------------------
      curr_Vgrid_d => this%curr%get_Data()
      if(allocated(curr_Vgrid_d%grid_var)) deallocate(curr_Vgrid_d%grid_var)
      if(allocated(curr_Vgrid_d%req_var )) deallocate(curr_Vgrid_d%req_var )
      if(allocated(curr_Vgrid_d%HYA     )) deallocate(curr_Vgrid_d%HYA     )
      if(allocated(curr_Vgrid_d%HYB     )) deallocate(curr_Vgrid_d%HYB     )
      if(allocated(curr_Vgrid_d%PRES    )) deallocate(curr_Vgrid_d%PRES    )
      if(allocated(curr_Vgrid_d%HGHT    )) deallocate(curr_Vgrid_d%HGHT    )
      deallocate(curr_Vgrid_d)
      deallocate(this%curr)

      ! Reset fisrt and curr pointers
      !-----------------------------------
      this%first => nextPointer
      this%curr  => nextPointer

      nullify(lastPointer)
      nullify(nextPointer)
    elseif(this%isLast()) then
      ! We are removing the last Link in the list
      !-----------------------------------------------
      lastPointer => this%curr%get_Prev()
      nextPointer => null()
      call lastPointer%set_Next(nextPointer)

      ! Deallocate curr
      !--------------------
      curr_Vgrid_d => this%curr%get_Data()
      if(allocated(curr_Vgrid_d%grid_var)) deallocate(curr_Vgrid_d%grid_var)
      if(allocated(curr_Vgrid_d%req_var )) deallocate(curr_Vgrid_d%req_var )
      if(allocated(curr_Vgrid_d%HYA     )) deallocate(curr_Vgrid_d%HYA     )
      if(allocated(curr_Vgrid_d%HYB     )) deallocate(curr_Vgrid_d%HYB     )
      if(allocated(curr_Vgrid_d%PRES    )) deallocate(curr_Vgrid_d%PRES    )
      if(allocated(curr_Vgrid_d%HGHT    )) deallocate(curr_Vgrid_d%HGHT    )
      deallocate(curr_Vgrid_d)
      deallocate(this%curr)

      ! Reset last and curr pointers
      !-----------------------------------
      this%last => lastPointer
      this%curr => lastPointer

      nullify(lastPointer)
      nullify(nextPointer)
    else
      ! We are removing a Link in the middle of the list
      !-----------------------------------------------------
      lastPointer => this%curr%get_Prev()
      nextPointer => this%curr%get_Next()
      call lastPointer%set_Next(nextPointer)
      call nextPointer%set_Prev(lastPointer)

      ! Deallocate curr
      !--------------------
      curr_Vgrid_d => this%curr%get_Data()
      if(allocated(curr_Vgrid_d%grid_var)) deallocate(curr_Vgrid_d%grid_var)
      if(allocated(curr_Vgrid_d%req_var )) deallocate(curr_Vgrid_d%req_var )
      if(allocated(curr_Vgrid_d%HYA     )) deallocate(curr_Vgrid_d%HYA     )
      if(allocated(curr_Vgrid_d%HYB     )) deallocate(curr_Vgrid_d%HYB     )
      if(allocated(curr_Vgrid_d%PRES    )) deallocate(curr_Vgrid_d%PRES    )
      if(allocated(curr_Vgrid_d%HGHT    )) deallocate(curr_Vgrid_d%HGHT    )
      deallocate(curr_Vgrid_d)
      deallocate(this%curr)

      ! Reset curr pointer
      !-----------------------------------
      this%curr => lastPointer

      nullify(lastPointer)
      nullify(nextPointer)
    endif

    return
  end subroutine ListVgrid_d_rmLink
  !========================================================================
  function ListVgrid_d_getLinkID(this)
    ! Return the ID idxVgrid_d for curr in List
    !---------------------------------------------
    class(ListVgrid_d_t):: this
    integer             :: ListVgrid_d_getLinkID
    ListVgrid_d_getLinkID = this%curr%get_ID()
    return
  end function ListVgrid_d_getLinkID
  !========================================================================
  function ListVgrid_d_getLinkData(this,idxVgrid_d)
    ! Given the Link ID, return a pointer to 
    ! the corresponding Vgrid_d data structure.
    !---------------------------------------------
    class(ListVgrid_d_t) :: this
    integer              :: idxVgrid_d
    type(vgrid_t),pointer:: ListVgrid_d_getLinkData

    ! cycle thru the list to find the Link with matching ID
    !----------------------------------------------------------
    call this%reset
    do while(this%curr%get_ID() .ne. idxVgrid_d)
      call this%iterate
      if(.not.associated(this%curr)) then
        ! ERROR: idxVgrid_d not found in list
        !-------------------------------------
        write(iulog, *) 'ListVgrid_d ERROR: No ',idxVgrid_d,'index in list'
        call endrun('ERROR: ListVgrid_d Missing Grid Index')
      endif
    end do
    ListVgrid_d_getLinkData => this%curr%get_Data()
    return
  end function ListVgrid_d_getLinkData
  !========================================================================
  function ListVgrid_d_setLinkID(this,new_Vgrid_d,newLink)
    ! Return the ID idxVgrid_d for new_Vgrid_d. 
    ! If the grid exists on the list, then return the corresponding ID, if 
    ! this is a new grid, add it to the end of the list and return the new ID.
    ! Set newLink to indicate if this is a new or existing link
    !-------------------------------------------------------------------------
    class(ListVgrid_d_t) :: this
    type(vgrid_t),pointer:: new_Vgrid_d
    integer              :: ListVgrid_d_setLinkID
    logical              :: newLink
    ! 
    ! Local Values
    !-------------------
    logical:: link_found

    ! Cycle thru the list and see if there is already an  
    ! existing grid the matches the new Vgrid_d
    !----------------------------------------------------
    call this%reset
    link_found = .false.
    do while(associated(this%curr))
      if(this%curr%isEqual(new_Vgrid_d)) then
        link_found = .true.
        exit
      endif
      call this%iterate 
    end do

    if(link_found) then
      ! if the Grid exists, pass back the the corresponding ID. 
      !-----------------------------------------------------------------
      newLink = .false.
      ListVgrid_d_setLinkID = this%curr%get_ID()
    else
      ! Otherwise add the new_Vgrid_d to the end of the 
      ! list and return the new ID
      !--------------------------------------------------
      newLink = .true.
      ListVgrid_d_setLinkID = this%addLink(new_Vgrid_d)
    endif

    return
  end function ListVgrid_d_setLinkID
  !========================================================================
  subroutine ListVgrid_d_print(this)
    ! Cycle thru the current list and print 
    ! out the meta data for each link 
    !---------------------------------------------
    class(ListVgrid_d_t)        :: this
    !
    ! Local values
    !-------------
    logical:: empty_list

    ! Check to see if the list is empty
    !-----------------------------------
    call this%reset
    empty_list = .false.
    if(this%isFirst().and.this%isLast()) then
      if(.not.associated(this%curr)) then
        empty_list = .true.
      endif
    endif
    
    if(empty_list) then
      ! Print out "Empty List" message here
      !---------------------------------------
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    LIST of Registered DataSet Vgrid_d structures:'
      write(iulog,*) '          EMPTY LIST  '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) ' '
    else
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    LIST of Registered DataSet Vgrid_d structures:'
      write(iulog,*) ' '

      ! cycle thru the list
      !----------------------
      do while(associated(this%curr))
        ! Issue a print command for each link
        !-------------------------------------
        call this%curr%printLink
        call this%iterate
      end do
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) ' '
    endif

    return
  end subroutine ListVgrid_d_print
  !========================================================================



  !========================================================================
  !========================================================================
  ! Procedures For LinkVgrid_d Structure
  !
  !========================================================================
  !========================================================================
  function LinkVgrid_d_getID(this)
    ! Return the ID for this Vgrid_d Link
    !----------------------------------------------------
    class(LinkVgrid_d_t):: this
    integer             :: LinkVgrid_d_getID
    LinkVgrid_d_getID = this%idxVgrid_d
    return
  end function LinkVgrid_d_getID
  !========================================================================
  function LinkVgrid_d_getData(this)
    ! Return the pointer to the next Vgrid_d Link
    !----------------------------------------------------
    class(LinkVgrid_d_t)        :: this
    type(vgrid_t)       ,pointer:: LinkVgrid_d_getData
    LinkVgrid_d_getData => this%Vgrid_d
    return
  end function LinkVgrid_d_getData
  !========================================================================
  function LinkVgrid_d_getNext(this)
    ! Return the pointer to the next Vgrid_d Link
    !----------------------------------------------------
    class(LinkVgrid_d_t)        :: this
    class(LinkVgrid_d_t),pointer:: LinkVgrid_d_getNext
    LinkVgrid_d_getNext => this%next
    return
  end function LinkVgrid_d_getNext
  !========================================================================
  subroutine LinkVgrid_d_setNext(this,next)
    ! Set the pointer for the next Vgrid_d Link
    !----------------------------------------------------
    class(LinkVgrid_d_t)        :: this
    class(LinkVgrid_d_t),pointer:: next
    this%next => next
    return
  end subroutine LinkVgrid_d_setNext
  !========================================================================
  function LinkVgrid_d_getPrev(this)
    ! Return the pointer to the previous Vgrid_d Link
    !----------------------------------------------------
    class(LinkVgrid_d_t)        :: this
    class(LinkVgrid_d_t),pointer:: LinkVgrid_d_getPrev
    LinkVgrid_d_getPrev => this%prev
    return
  end function LinkVgrid_d_getPrev
  !========================================================================
  subroutine LinkVgrid_d_setPrev(this,prev)
    ! Return the pointer to the previous Hgrid_d Link
    !----------------------------------------------------
    class(LinkVgrid_d_t)        :: this
    class(LinkVgrid_d_t),pointer:: prev
    this%prev => prev
    return
  end subroutine LinkVgrid_d_setPrev
  !========================================================================
  function LinkVgrid_d_isEqual(this,new_Vgrid_d)
    ! Logical function which returns True if the given 
    ! Vgrid_d structure data matches those of the current Link.
    ! Other wise the structures differ and False is returned
    !----------------------------------------------------
    class(LinkVgrid_d_t) :: this
    type(vgrid_t),pointer:: new_Vgrid_d
    logical              :: LinkVgrid_d_isEqual
    !
    ! Local values
    !-----------------
    type(vgrid_t),pointer:: My_Vgrid_d
    real(r8)             :: diffSum
    integer              :: ii

    ! Get a pointer to the current data structure
    !--------------------------------------------
    My_Vgrid_d => this%get_Data()

    ! Set the restult to true, then scan the structure for differences, 
    ! if/when a difference is found, set the result to false and return 
    ! immediately.
    !-------------------------------------------------------------------
    LinkVgrid_d_isEqual = .true.

    if(trim(My_Vgrid_d%name).ne.trim(new_Vgrid_d%name)) then
      LinkVgrid_d_isEqual = .false.
      return
    endif

    do ii=1,My_Vgrid_d%ngrid_var
      if(trim(My_Vgrid_d%grid_var(ii)).ne.trim(new_Vgrid_d%grid_var(ii))) then
        LinkVgrid_d_isEqual = .false.
        return
      endif
    end do

    if(My_Vgrid_d%nlev.ne.new_Vgrid_d%nlev) then
      LinkVgrid_d_isEqual = .false.
      return
    endif
    
    return
  end function LinkVgrid_d_isEqual
  !========================================================================
  subroutine LinkVgrid_d_print(this)
    ! Print out the Vgrid_d metadata for the current link
    !---------------------------------------------
    class(LinkVgrid_d_t):: this
    !
    ! Local values
    !---------------
    integer:: ii

    ! Print out identifying values
    !--------------------------------
    write(iulog,*) '      -------------------------------------------'
    write(iulog,*) '      idxVgrid_d=',this%idxVgrid_d
    write(iulog,*) '                     name= ',trim(this%Vgrid_d%name)
    write(iulog,*) '                     type= ',trim(this%Vgrid_d%Type)
    write(iulog,*) '                  dimname= ',trim(this%Vgrid_d%Vdimname)
    write(iulog,*) '                     nlev= ',this%Vgrid_d%nlev
    write(iulog,*) '                ngrid_var= ',this%Vgrid_d%ngrid_var
    do ii=1,this%Vgrid_d%ngrid_var
      write(iulog,*) '             Grid Variable',ii,                           &
                                         ' name=',trim(this%Vgrid_d%grid_var(ii))
    end do
    write(iulog,*) '                 nreq_var= ',this%Vgrid_d%nreq_var
    do ii=1,this%Vgrid_d%nreq_var
      write(iulog,*) '           Required Variable',ii,                          &
                                           ' name=',trim(this%Vgrid_d%req_var(ii))
    end do
    write(iulog,*) '                      Pref = ',this%Vgrid_d%Pref
    if(allocated(this%Vgrid_d%HYA)) then
      write(iulog,*) '                     HYA = ALLOCATED'
    else
      write(iulog,*) '                     HYA = NOT ALLOCATED'
    endif
    if(allocated(this%Vgrid_d%HYB)) then
      write(iulog,*) '                     HYB = ALLOCATED'
    else
      write(iulog,*) '                     HYB = NOT ALLOCATED'
    endif
    if(allocated(this%Vgrid_d%PRES)) then
      write(iulog,*) '                    PRES = ALLOCATED'
    else
      write(iulog,*) '                    PRES = NOT ALLOCATED'
    endif
    if(allocated(this%Vgrid_d%HGHT)) then
      write(iulog,*) '                    HGHT = ALLOCATED'
    else
      write(iulog,*) '                    HGHT = NOT ALLOCATED'
    endif
    write(iulog,*) '      -------------------------------------------'
    write(iulog,*) ' '

    return
  end subroutine LinkVgrid_d_print
  !========================================================================



  !========================================================================
  !========================================================================
  ! Procedures For ListHmaps Structure
  !
  !========================================================================
  !========================================================================
  subroutine ListHmaps_reset(this)
    ! Reset the Iterator to point at the first Hmap Link
    !-------------------------------------------------------
    class(ListHmaps_t):: this
    this%curr => this%first
    return
  end subroutine ListHmaps_reset
  !========================================================================
  subroutine ListHmaps_iterate(this)
    ! Increment the Iterator to point at the next Hmap Link
    !----------------------------------------------------------
    class(ListHmaps_t):: this
    this%curr => this%curr%get_Next()
    return
  end subroutine ListHmaps_iterate
  !========================================================================
  function ListHmaps_isFirst(this)
    ! Check if the Iterator points to the end of the list.
    !----------------------------------------------------
    class(ListHmaps_t)        :: this
    logical                   :: ListHmaps_isFirst
    class(LinkHmaps_t),pointer:: lastPointer
    if(associated(this%curr)) then
      lastPointer => this%curr%get_Prev()
      ListHmaps_isFirst = .not.associated(lastPointer)
    else
      ListHmaps_isFirst = .true.
    endif
    return
  end function ListHmaps_isFirst
  !========================================================================
  function ListHmaps_isLast(this)
    ! Check if the Iterator points to the end of the list.
    !----------------------------------------------------
    class(ListHmaps_t)        :: this
    logical                   :: ListHmaps_isLast
    class(LinkHmaps_t),pointer:: nextPointer
    if(associated(this%curr)) then
      nextPointer => this%curr%get_Next()
      ListHmaps_isLast = .not.associated(nextPointer)
    else
      ListHmaps_isLast = .true.
    endif
    return
  end function ListHmaps_isLast
  !========================================================================
  function ListHmaps_addLink(this,Hmap)
    ! Add a new Hmap Dataset link to the end of the list.
    !========================================================
    class(ListHmaps_t)      :: this
    type(HdataMap_t),pointer:: Hmap
    integer                 :: ListHmaps_addLink
    !
    ! Local Values
    !--------------
    class(LinkHmaps_t),pointer:: newLink

    ! Allocate and a new Hmap Link and initialize 
    ! the pointer with Hmap from the input.
    !----------------------------------------------------
    allocate(newLink)
    newLink%Hmap => Hmap

    ! Increment the Counter to assign a unique index 
    ! to this new Link
    !-----------------------------------------------
    this%Counter    = this%Counter + 1
    newLink%idxHmap = this%Counter

    ! Add it to the list
    !--------------------
    if(.not.associated(this%last)) then
      ! If the current list is empty, initialize 
      ! the list with this Link
      !-------------------------------------------
      newLink%prev => null()
      newLink%next => null()
      this%first   => newLink
      this%last    => newLink
      this%curr    => newLink
    else
      ! Append this Link to the end of the current list.
      !----------------------------------------------------
      newLink%prev   => this%last
      newLink%next   => null()
      this%last%next => newLink
      this%last      => newLink
    endif

    ! Return the DataSet ID to the calling routine
    !-----------------------------------------------
    ListHmaps_addLink = newLink%idxHmap

    return
  end function ListHmaps_addLink
  !========================================================================
  subroutine ListHmaps_rmLink(this,idxHmap)
    ! Remove Hmap Dataset link with given ID
    !========================================================
    class(ListHmaps_t):: this
    integer           :: idxHmap
    !
    ! Local Values
    !--------------
    class(LinkHmaps_t),pointer:: lastPointer
    class(LinkHmaps_t),pointer:: nextPointer
    type(HdataMap_t)  ,pointer:: curr_Hmap

    ! cycle thru the list to find the Link with matching ID
    !----------------------------------------------------------
    call this%reset
    do while(this%curr%get_ID() .ne. idxHmap)
      call this%iterate
      if(.not.associated(this%curr)) then
        ! ERROR: idxHmap not found in list
        !-------------------------------------
        write(iulog, *) 'ListHmaps ERROR: No ',idxHmap,'index in list'
        call endrun('ERROR: ListHmaps Missing Grid Index')
      endif
    end do

    ! Now we want to remove the current 
    ! Link from the list
    !------------------------------------
    if(this%isFirst().and.this%isLast()) then
      ! We are removing the ONLY link in the List
      !--------------------------------------------
      curr_Hmap => this%curr%get_Data()
      
      ! Deallocate curr Link
      !----------------------
      if(allocated(curr_Hmap%HmapData%lon    )) deallocate(curr_Hmap%HmapData%lon    )
      if(allocated(curr_Hmap%HmapData%lat    )) deallocate(curr_Hmap%HmapData%lat    )
      if(allocated(curr_Hmap%HmapData%NbrList)) deallocate(curr_Hmap%HmapData%NbrList)
      if(allocated(curr_Hmap%HmapData%Vparity)) deallocate(curr_Hmap%HmapData%Vparity)
      if(allocated(curr_Hmap%HmapData%lonInd )) deallocate(curr_Hmap%HmapData%lonInd )
      if(allocated(curr_Hmap%HmapData%latInd )) deallocate(curr_Hmap%HmapData%latInd )
      if(allocated(curr_Hmap%HmapData%dof2D  )) deallocate(curr_Hmap%HmapData%dof2D  )
      if(allocated(curr_Hmap%HmapData%NbrIwgt)) deallocate(curr_Hmap%HmapData%NbrIwgt)
      deallocate(curr_Hmap)
      deallocate(this%curr)
   
      ! Nullify all pointers for an empty list
      !----------------------------------------
      this%first => null()
      this%curr  => null()
      this%last  => null()
    elseif(this%isFirst()) then
      ! We are removing the first Link in the list
      !-----------------------------------------------
      lastPointer => null()
      nextPointer => this%curr%get_Next()
      call nextPointer%set_Prev(lastPointer)

      ! Deallocate curr Link
      !----------------------
      curr_Hmap => this%curr%get_Data()
      if(allocated(curr_Hmap%HmapData%lon    )) deallocate(curr_Hmap%HmapData%lon    )
      if(allocated(curr_Hmap%HmapData%lat    )) deallocate(curr_Hmap%HmapData%lat    )
      if(allocated(curr_Hmap%HmapData%NbrList)) deallocate(curr_Hmap%HmapData%NbrList)
      if(allocated(curr_Hmap%HmapData%Vparity)) deallocate(curr_Hmap%HmapData%Vparity)
      if(allocated(curr_Hmap%HmapData%lonInd )) deallocate(curr_Hmap%HmapData%lonInd )
      if(allocated(curr_Hmap%HmapData%latInd )) deallocate(curr_Hmap%HmapData%latInd )
      if(allocated(curr_Hmap%HmapData%dof2D  )) deallocate(curr_Hmap%HmapData%dof2D  )
      if(allocated(curr_Hmap%HmapData%NbrIwgt)) deallocate(curr_Hmap%HmapData%NbrIwgt)
      deallocate(curr_Hmap)
      deallocate(this%curr)

      ! Reset fisrt and curr pointers
      !-----------------------------------
      this%first => nextPointer
      this%curr  => nextPointer

      nullify(lastPointer)
      nullify(nextPointer)
    elseif(this%isLast()) then
      ! We are removing the last Link in the list
      !-----------------------------------------------
      lastPointer => this%curr%get_Prev()
      nextPointer => null()
      call lastPointer%set_Next(nextPointer)

      ! Deallocate curr
      !--------------------
      curr_Hmap => this%curr%get_Data()
      if(allocated(curr_Hmap%HmapData%lon    )) deallocate(curr_Hmap%HmapData%lon    )
      if(allocated(curr_Hmap%HmapData%lat    )) deallocate(curr_Hmap%HmapData%lat    )
      if(allocated(curr_Hmap%HmapData%NbrList)) deallocate(curr_Hmap%HmapData%NbrList)
      if(allocated(curr_Hmap%HmapData%Vparity)) deallocate(curr_Hmap%HmapData%Vparity)
      if(allocated(curr_Hmap%HmapData%lonInd )) deallocate(curr_Hmap%HmapData%lonInd )
      if(allocated(curr_Hmap%HmapData%latInd )) deallocate(curr_Hmap%HmapData%latInd )
      if(allocated(curr_Hmap%HmapData%dof2D  )) deallocate(curr_Hmap%HmapData%dof2D  )
      if(allocated(curr_Hmap%HmapData%NbrIwgt)) deallocate(curr_Hmap%HmapData%NbrIwgt)
      deallocate(curr_Hmap)
      deallocate(this%curr)

      ! Reset last and curr pointers
      !-----------------------------------
      this%last => lastPointer
      this%curr => lastPointer

      nullify(lastPointer)
      nullify(nextPointer)
    else
      ! We are removing a Link in the middle of the list
      !-----------------------------------------------------
      lastPointer => this%curr%get_Prev()
      nextPointer => this%curr%get_Next()
      call lastPointer%set_Next(nextPointer)
      call nextPointer%set_Prev(lastPointer)

      ! Deallocate curr
      !--------------------
      curr_Hmap => this%curr%get_Data()
      if(allocated(curr_Hmap%HmapData%lon    )) deallocate(curr_Hmap%HmapData%lon    )
      if(allocated(curr_Hmap%HmapData%lat    )) deallocate(curr_Hmap%HmapData%lat    )
      if(allocated(curr_Hmap%HmapData%NbrList)) deallocate(curr_Hmap%HmapData%NbrList)
      if(allocated(curr_Hmap%HmapData%Vparity)) deallocate(curr_Hmap%HmapData%Vparity)
      if(allocated(curr_Hmap%HmapData%lonInd )) deallocate(curr_Hmap%HmapData%lonInd )
      if(allocated(curr_Hmap%HmapData%latInd )) deallocate(curr_Hmap%HmapData%latInd )
      if(allocated(curr_Hmap%HmapData%dof2D  )) deallocate(curr_Hmap%HmapData%dof2D  )
      if(allocated(curr_Hmap%HmapData%NbrIwgt)) deallocate(curr_Hmap%HmapData%NbrIwgt)
      deallocate(curr_Hmap)
      deallocate(this%curr)

      ! Reset curr pointer
      !-----------------------------------
      this%curr => lastPointer

      nullify(lastPointer)
      nullify(nextPointer)
    endif

    return
  end subroutine ListHmaps_rmLink
  !========================================================================
  function ListHmaps_getLinkID(this)
    ! Return the ID idxHmap for curr in List
    !---------------------------------------------
    class(ListHmaps_t):: this
    integer           :: ListHmaps_getLinkID
    ListHmaps_getLinkID = this%curr%get_ID()
    return
  end function ListHmaps_getLinkID
  !========================================================================
  function ListHmaps_getLinkData(this,idxHmap)
    ! Given the Link ID, return a pointer to 
    ! the corresponding Hmap data structure.
    !---------------------------------------------
    class(ListHmaps_t)      :: this
    integer                 :: idxHmap
    type(HdataMap_t),pointer:: ListHmaps_getLinkData

    ! cycle thru the list to find the Link with matching ID
    !----------------------------------------------------------
    call this%reset
    do while(this%curr%get_ID() .ne. idxHmap)
      call this%iterate
      if(.not.associated(this%curr)) then
        ! ERROR: idxHmap not found in list
        !-------------------------------------
        write(iulog, *) 'ListHmaps ERROR: No ',idxHmap,'index in list'
        call endrun('ERROR: ListHmaps Missing Grid Index')
      endif
    end do
    ListHmaps_getLinkData => this%curr%get_Data()
    return
  end function ListHmaps_getLinkData
  !========================================================================
  function ListHmaps_setLinkID(this,new_Hmap,newLink)
    ! Return the ID idxHmap for new_Hmap. 
    ! If the grid exists on the list, then return the corresponding ID, if 
    ! this is a new grid, add it to the end of the list and return the new ID.
    ! Set newLink to indicate if this is a new or existing link
    !-------------------------------------------------------------------------
    class(ListHmaps_t)      :: this
    type(HdataMap_t),pointer:: new_Hmap
    integer                 :: ListHmaps_setLinkID
    logical                 :: newLink
    ! 
    ! Local Values
    !-------------------
    logical:: link_found

    ! Cycle thru the list and see if there is already an  
    ! existing grid the matches the new Hmap
    !----------------------------------------------------
    call this%reset
    link_found = .false.
    do while(associated(this%curr))
      if(this%curr%isEqual(new_Hmap)) then
        link_found = .true.
        exit
      endif
      call this%iterate 
    end do

    if(link_found) then
      ! if the Grid exists, pass back the the corresponding ID. 
      !-----------------------------------------------------------------
      newLink = .false.
      ListHmaps_setLinkID = this%curr%get_ID()
    else
      ! Otherwise add the new_Hmap to the end of the 
      ! list and return the new ID
      !--------------------------------------------------
      newLink = .true.
      ListHmaps_setLinkID = this%addLink(new_Hmap)
    endif

    return
  end function ListHmaps_setLinkID
  !========================================================================
  subroutine ListHmaps_print(this)
    ! Cycle thru the current list and print 
    ! out the meta data for each link 
    !---------------------------------------------
    class(ListHmaps_t)        :: this
    !
    ! Local values
    !-------------
    logical:: empty_list

    ! Check to see if the list is empty
    !-----------------------------------
    call this%reset
    empty_list = .false.
    if(this%isFirst().and.this%isLast()) then
      if(.not.associated(this%curr)) then
        empty_list = .true.
      endif
    endif
    
    if(empty_list) then
      ! Print out "Empty List" message here
      !---------------------------------------
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    LIST of Registered DataSet Hmap structures:'
      write(iulog,*) '          EMPTY LIST  '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) ' '
    else
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    LIST of Registered DataSet Hmap structures:'
      write(iulog,*) ' '

      ! cycle thru the list
      !----------------------
      do while(associated(this%curr))
        ! Issue a print command for each link
        !-------------------------------------
        call this%curr%printLink
        call this%iterate
      end do
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) ' '
    endif

    return
  end subroutine ListHmaps_print
  !========================================================================



  !========================================================================
  !========================================================================
  ! Procedures For LinkHmaps Structure
  !
  !========================================================================
  !========================================================================
  function LinkHmaps_getID(this)
    ! Return the ID for this Hmap Link
    !----------------------------------------------------
    class(LinkHmaps_t):: this
    integer           :: LinkHmaps_getID
    LinkHmaps_getID = this%idxHmap
    return
  end function LinkHmaps_getID
  !========================================================================
  function LinkHmaps_getData(this)
    ! Return the pointer to the next Hmap Link
    !----------------------------------------------------
    class(LinkHmaps_t)      :: this
    type(HdataMap_t),pointer:: LinkHmaps_getData
    LinkHmaps_getData => this%Hmap
    return
  end function LinkHmaps_getData
  !========================================================================
  function LinkHmaps_getNext(this)
    ! Return the pointer to the next Hmap Link
    !----------------------------------------------------
    class(LinkHmaps_t)        :: this
    class(LinkHmaps_t),pointer:: LinkHmaps_getNext
    LinkHmaps_getNext => this%next
    return
  end function LinkHmaps_getNext
  !========================================================================
  subroutine LinkHmaps_setNext(this,next)
    ! Set the pointer for the next Hmap Link
    !----------------------------------------------------
    class(LinkHmaps_t)        :: this
    class(LinkHmaps_t),pointer:: next
    this%next => next
    return
  end subroutine LinkHmaps_setNext
  !========================================================================
  function LinkHmaps_getPrev(this)
    ! Return the pointer to the previous Hmap Link
    !----------------------------------------------------
    class(LinkHmaps_t)        :: this
    class(LinkHmaps_t),pointer:: LinkHmaps_getPrev
    LinkHmaps_getPrev => this%prev
    return
  end function LinkHmaps_getPrev
  !========================================================================
  subroutine LinkHmaps_setPrev(this,prev)
    ! Return the pointer to the previous Hgrid_d Link
    !----------------------------------------------------
    class(LinkHmaps_t)        :: this
    class(LinkHmaps_t),pointer:: prev
    this%prev => prev
    return
  end subroutine LinkHmaps_setPrev
  !========================================================================
  function LinkHmaps_isEqual(this,new_Hmap)
    ! Logical function which returns True if the given 
    ! Hmap structure data matches those of the current Link.
    ! Other wise the structures differ and False is returned
    !----------------------------------------------------
    class(LinkHmaps_t)      :: this
    type(HdataMap_t),pointer:: new_Hmap
    logical                 :: LinkHmaps_isEqual
    !
    ! Local values
    !-----------------
    type(HdataMap_t),pointer:: My_Hmap
    real(r8)                :: diffSum
    integer                 :: ii

    ! Get a pointer to the current data structure
    !--------------------------------------------
    My_Hmap => this%get_Data()

    ! Set the restult to true, then scan the structure for differences, 
    ! if/when a difference is found, set the result to false and return 
    ! immediately.
    !-------------------------------------------------------------------
    LinkHmaps_isEqual = .true.

    if(trim(My_Hmap%HmapOpt).ne.trim(new_Hmap%HmapOpt)) then
      LinkHmaps_isEqual = .false.
      return
    endif

    if(trim(My_Hmap%HmapFunction).ne.trim(new_Hmap%HmapFunction)) then
      LinkHmaps_isEqual = .false.
      return
    endif

    if(My_Hmap%idxHgrid_d.ne.new_Hmap%idxHgrid_d) then
      LinkHmaps_isEqual = .false.
      return
    endif

    if(My_Hmap%idxHgrid_m.ne.new_Hmap%idxHgrid_m) then
      LinkHmaps_isEqual = .false.
      return
    endif

    return
  end function LinkHmaps_isEqual
  !========================================================================
  subroutine LinkHmaps_print(this)
    ! Print out the Hmap metadata for the current link
    !---------------------------------------------
    class(LinkHmaps_t):: this
    !
    ! Local values
    !---------------
    integer:: ii

    ! Print out identifying values
    !--------------------------------
    write(iulog,*) '      -------------------------------------------'
    write(iulog,*) '      idxHmap=',this%idxHmap
    write(iulog,*) '                 Opt name= ',trim(this%Hmap%HmapOpt)
    write(iulog,*) '                 function= ',trim(this%Hmap%HmapFunction)
    write(iulog,*) '               idxHgrid_d= ',this%Hmap%idxHgrid_d
    write(iulog,*) '               idxHgrid_m= ',this%Hmap%idxHgrid_m
    write(iulog,*) ' '
    write(iulog,*) '            HmapData%nlon= ',this%Hmap%HmapData%nlon
    write(iulog,*) '            HmapData%nlat= ',this%Hmap%HmapData%nlat
    write(iulog,*) '            HmapData%ncol= ',this%Hmap%HmapData%ncol
    write(iulog,*) '         HmapData%HmapOpt= ',trim(this%Hmap%HmapData%HmapOpt)
    write(iulog,*) '      -------------------------------------------'
    write(iulog,*) ' '

    return
  end subroutine LinkHmaps_print
  !========================================================================



  !========================================================================
  !========================================================================
  ! Procedures For ListVmaps Structure
  !
  !========================================================================
  !========================================================================
  subroutine ListVmaps_reset(this)
    ! Reset the Iterator to point at the first Vmap Link
    !-------------------------------------------------------
    class(ListVmaps_t):: this
    this%curr => this%first
    return
  end subroutine ListVmaps_reset
  !========================================================================
  subroutine ListVmaps_iterate(this)
    ! Increment the Iterator to point at the next Vmap Link
    !----------------------------------------------------------
    class(ListVmaps_t):: this
    this%curr => this%curr%get_Next()
    return
  end subroutine ListVmaps_iterate
  !========================================================================
  function ListVmaps_isFirst(this)
    ! Check if the Iterator points to the end of the list.
    !----------------------------------------------------
    class(ListVmaps_t)        :: this
    logical                   :: ListVmaps_isFirst
    class(LinkVmaps_t),pointer:: lastPointer
    if(associated(this%curr)) then
      lastPointer => this%curr%get_Prev()
      ListVmaps_isFirst = .not.associated(lastPointer)
    else
      ListVmaps_isFirst = .true.
    endif
    return
  end function ListVmaps_isFirst
  !========================================================================
  function ListVmaps_isLast(this)
    ! Check if the Iterator points to the end of the list.
    !----------------------------------------------------
    class(ListVmaps_t)        :: this
    logical                   :: ListVmaps_isLast
    class(LinkVmaps_t),pointer:: nextPointer
    if(associated(this%curr)) then
      nextPointer => this%curr%get_Next()
      ListVmaps_isLast = .not.associated(nextPointer)
    else
      ListVmaps_isLast = .true.
    endif
    return
  end function ListVmaps_isLast
  !========================================================================
  function ListVmaps_addLink(this,Vmap)
    ! Add a new Vmap Dataset link to the end of the list.
    !========================================================
    class(ListVmaps_t)  :: this
    type(vmap_t),pointer:: Vmap
    integer             :: ListVmaps_addLink
    !
    ! Local Values
    !--------------
    class(LinkVmaps_t),pointer:: newLink

    ! Allocate and a new Vmap Link and initialize 
    ! the pointer with Vmap from the input.
    !----------------------------------------------------
    allocate(newLink)
    newLink%Vmap => Vmap

    ! Increment the Counter to assign a unique index 
    ! to this new Link
    !-----------------------------------------------
    this%Counter    = this%Counter + 1
    newLink%idxVmap = this%Counter

    ! Add it to the list
    !--------------------
    if(.not.associated(this%last)) then
      ! If the current list is empty, initialize 
      ! the list with this Link
      !-------------------------------------------
      newLink%prev => null()
      newLink%next => null()
      this%first   => newLink
      this%last    => newLink
      this%curr    => newLink
    else
      ! Append this Link to the end of the current list.
      !----------------------------------------------------
      newLink%prev   => this%last
      newLink%next   => null()
      this%last%next => newLink
      this%last      => newLink
    endif

    ! Return the DataSet ID to the calling routine
    !-----------------------------------------------
    ListVmaps_addLink = newLink%idxVmap

    return
  end function ListVmaps_addLink
  !========================================================================
  subroutine ListVmaps_rmLink(this,idxVmap)
    ! Remove Vmap Dataset link with given ID
    !========================================================
    class(ListVmaps_t):: this
    integer           :: idxVmap
    !
    ! Local Values
    !--------------
    class(LinkVmaps_t),pointer:: lastPointer
    class(LinkVmaps_t),pointer:: nextPointer
    type(vmap_t)      ,pointer:: curr_Vmap

    ! cycle thru the list to find the Link with matching ID
    !----------------------------------------------------------
    call this%reset
    do while(this%curr%get_ID() .ne. idxVmap)
      call this%iterate
      if(.not.associated(this%curr)) then
        ! ERROR: idxVmap not found in list
        !-------------------------------------
        write(iulog, *) 'ListVmaps ERROR: No ',idxVmap,'index in list'
        call endrun('ERROR: ListVmaps Missing Grid Index')
      endif
    end do

    ! Now we want to remove the current 
    ! Link from the list
    !------------------------------------
    if(this%isFirst().and.this%isLast()) then
      ! We are removing the ONLY link in the List
      !--------------------------------------------
      curr_Vmap => this%curr%get_Data()
      
      ! Deallocate curr Link
      !----------------------
      deallocate(curr_Vmap)
      deallocate(this%curr)
   
      ! Nullify all pointers for an empty list
      !----------------------------------------
      this%first => null()
      this%curr  => null()
      this%last  => null()
    elseif(this%isFirst()) then
      ! We are removing the first Link in the list
      !-----------------------------------------------
      lastPointer => null()
      nextPointer => this%curr%get_Next()
      call nextPointer%set_Prev(lastPointer)

      ! Deallocate curr Link
      !----------------------
      curr_Vmap => this%curr%get_Data()
      deallocate(curr_Vmap)
      deallocate(this%curr)

      ! Reset fisrt and curr pointers
      !-----------------------------------
      this%first => nextPointer
      this%curr  => nextPointer

      nullify(lastPointer)
      nullify(nextPointer)
    elseif(this%isLast()) then
      ! We are removing the last Link in the list
      !-----------------------------------------------
      lastPointer => this%curr%get_Prev()
      nextPointer => null()
      call lastPointer%set_Next(nextPointer)

      ! Deallocate curr
      !--------------------
      curr_Vmap => this%curr%get_Data()
      deallocate(curr_Vmap)
      deallocate(this%curr)

      ! Reset last and curr pointers
      !-----------------------------------
      this%last => lastPointer
      this%curr => lastPointer

      nullify(lastPointer)
      nullify(nextPointer)
    else
      ! We are removing a Link in the middle of the list
      !-----------------------------------------------------
      lastPointer => this%curr%get_Prev()
      nextPointer => this%curr%get_Next()
      call lastPointer%set_Next(nextPointer)
      call nextPointer%set_Prev(lastPointer)

      ! Deallocate curr
      !--------------------
      curr_Vmap => this%curr%get_Data()
      deallocate(curr_Vmap)
      deallocate(this%curr)

      ! Reset curr pointer
      !-----------------------------------
      this%curr => lastPointer

      nullify(lastPointer)
      nullify(nextPointer)
    endif

    return
  end subroutine ListVmaps_rmLink
  !========================================================================
  function ListVmaps_getLinkID(this)
    ! Return the ID idxVmap for curr in List
    !---------------------------------------------
    class(ListVmaps_t):: this
    integer           :: ListVmaps_getLinkID
    ListVmaps_getLinkID = this%curr%get_ID()
    return
  end function ListVmaps_getLinkID
  !========================================================================
  function ListVmaps_getLinkData(this,idxVmap)
    ! Given the Link ID, return a pointer to 
    ! the corresponding Vmap data structure.
    !---------------------------------------------
    class(ListVmaps_t)  :: this
    integer             :: idxVmap
    type(vmap_t),pointer:: ListVmaps_getLinkData

    ! cycle thru the list to find the Link with matching ID
    !----------------------------------------------------------
    call this%reset
    do while(this%curr%get_ID() .ne. idxVmap)
      call this%iterate
      if(.not.associated(this%curr)) then
        ! ERROR: idxVmap not found in list
        !-------------------------------------
        write(iulog, *) 'ListVmaps ERROR: No ',idxVmap,'index in list'
        call endrun('ERROR: ListVmaps Missing Grid Index')
      endif
    end do
    ListVmaps_getLinkData => this%curr%get_Data()
    return
  end function ListVmaps_getLinkData
  !========================================================================
  function ListVmaps_setLinkID(this,new_Vmap,newLink)
    ! Return the ID idxVmap for new_Vmap. 
    ! If the grid exists on the list, then return the corresponding ID, if 
    ! this is a new grid, add it to the end of the list and return the new ID.
    ! Set newLink to indicate if this is a new or existing link
    !-------------------------------------------------------------------------
    class(ListVmaps_t)  :: this
    type(vmap_t),pointer:: new_Vmap
    integer             :: ListVmaps_setLinkID
    logical             :: newLink
    ! 
    ! Local Values
    !-------------------
    logical:: link_found

    ! Cycle thru the list and see if there is already an  
    ! existing grid the matches the new Vmap
    !----------------------------------------------------
    call this%reset
    link_found = .false.
    do while(associated(this%curr))
      if(this%curr%isEqual(new_Vmap)) then
        link_found = .true.
        exit
      endif
      call this%iterate 
    end do

    if(link_found) then
      ! if the Grid exists, pass back the the corresponding ID. 
      !-----------------------------------------------------------------
      newLink = .false.
      ListVmaps_setLinkID = this%curr%get_ID()
    else
      ! Otherwise add the new_Vmap to the end of the 
      ! list and return the new ID
      !--------------------------------------------------
      newLink = .true.
      ListVmaps_setLinkID = this%addLink(new_Vmap)
    endif

    return
  end function ListVmaps_setLinkID
  !========================================================================
  subroutine ListVmaps_print(this)
    ! Cycle thru the current list and print 
    ! out the meta data for each link 
    !---------------------------------------------
    class(ListVmaps_t):: this
    !
    ! Local values
    !-------------
    logical:: empty_list

    ! Check to see if the list is empty
    !-----------------------------------
    call this%reset
    empty_list = .false.
    if(this%isFirst().and.this%isLast()) then
      if(.not.associated(this%curr)) then
        empty_list = .true.
      endif
    endif
    
    if(empty_list) then
      ! Print out "Empty List" message here
      !---------------------------------------
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    LIST of Registered DataSet Vmap structures:'
      write(iulog,*) '          EMPTY LIST  '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) ' '
    else
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    LIST of Registered DataSet Vmap structures:'
      write(iulog,*) ' '

      ! cycle thru the list
      !----------------------
      do while(associated(this%curr))
        ! Issue a print command for each link
        !-------------------------------------
        call this%curr%printLink
        call this%iterate
      end do
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) ' '
    endif

    return
  end subroutine ListVmaps_print
  !========================================================================



  !========================================================================
  !========================================================================
  ! Procedures For LinkVmaps Structure
  !
  !========================================================================
  !========================================================================
  function LinkVmaps_getID(this)
    ! Return the ID for this Vmap Link
    !----------------------------------------------------
    class(LinkVmaps_t):: this
    integer           :: LinkVmaps_getID
    LinkVmaps_getID = this%idxVmap
    return
  end function LinkVmaps_getID
  !========================================================================
  function LinkVmaps_getData(this)
    ! Return the pointer to the next Vmap Link
    !----------------------------------------------------
    class(LinkVmaps_t)  :: this
    type(vmap_t),pointer:: LinkVmaps_getData
    LinkVmaps_getData => this%Vmap
    return
  end function LinkVmaps_getData
  !========================================================================
  function LinkVmaps_getNext(this)
    ! Return the pointer to the next Vmap Link
    !----------------------------------------------------
    class(LinkVmaps_t)        :: this
    class(LinkVmaps_t),pointer:: LinkVmaps_getNext
    LinkVmaps_getNext => this%next
    return
  end function LinkVmaps_getNext
  !========================================================================
  subroutine LinkVmaps_setNext(this,next)
    ! Set the pointer for the next Vmap Link
    !----------------------------------------------------
    class(LinkVmaps_t)        :: this
    class(LinkVmaps_t),pointer:: next
    this%next => next
    return
  end subroutine LinkVmaps_setNext
  !========================================================================
  function LinkVmaps_getPrev(this)
    ! Return the pointer to the previous Vmap Link
    !----------------------------------------------------
    class(LinkVmaps_t)        :: this
    class(LinkVmaps_t),pointer:: LinkVmaps_getPrev
    LinkVmaps_getPrev => this%prev
    return
  end function LinkVmaps_getPrev
  !========================================================================
  subroutine LinkVmaps_setPrev(this,prev)
    ! Return the pointer to the previous Hgrid_d Link
    !----------------------------------------------------
    class(LinkVmaps_t)        :: this
    class(LinkVmaps_t),pointer:: prev
    this%prev => prev
    return
  end subroutine LinkVmaps_setPrev
  !========================================================================
  function LinkVmaps_isEqual(this,new_Vmap)
    ! Logical function which returns True if the given 
    ! Vmap structure data matches those of the current Link.
    ! Other wise the structures differ and False is returned
    !----------------------------------------------------
    class(LinkVmaps_t)  :: this
    type(vmap_t),pointer:: new_Vmap
    logical             :: LinkVmaps_isEqual
    !
    ! Local values
    !-----------------
    type(vmap_t),pointer:: My_Vmap
    real(r8)            :: diffSum
    integer             :: ii

    ! Get a pointer to the current data structure
    !--------------------------------------------
    My_Vmap => this%get_Data()

    ! Set the restult to true, then scan the structure for differences, 
    ! if/when a difference is found, set the result to false and return 
    ! immediately.
    !-------------------------------------------------------------------
    LinkVmaps_isEqual = .true.

    if(trim(My_Vmap%VmapOpt).ne.trim(new_Vmap%VmapOpt)) then
      LinkVmaps_isEqual = .false.
      return
    endif

    if(trim(My_Vmap%VmapFunction).ne.trim(new_Vmap%VmapFunction)) then
      LinkVmaps_isEqual = .false.
      return
    endif

    if(My_Vmap%idxVgrid_d.ne.new_Vmap%idxVgrid_d) then
      LinkVmaps_isEqual = .false.
      return
    endif

    if(My_Vmap%idxVgrid_m.ne.new_Vmap%idxVgrid_m) then
      LinkVmaps_isEqual = .false.
      return
    endif

    return
  end function LinkVmaps_isEqual
  !========================================================================
  subroutine LinkVmaps_print(this)
    ! Print out the Vmap metadata for the current link
    !---------------------------------------------
    class(LinkVmaps_t):: this
    !
    ! Local values
    !---------------
    integer:: ii

    ! Print out identifying values
    !--------------------------------
    write(iulog,*) '      -------------------------------------------'
    write(iulog,*) '      idxVmap=',this%idxVmap
    write(iulog,*) '                 Opt name= ',trim(this%Vmap%VmapOpt)
    write(iulog,*) '                 function= ',trim(this%Vmap%VmapFunction)
    write(iulog,*) '               idxVgrid_d= ',this%Vmap%idxVgrid_d
    write(iulog,*) '               idxVgrid_m= ',this%Vmap%idxVgrid_m
    write(iulog,*) '      -------------------------------------------'
    write(iulog,*) ' '

    return
  end subroutine LinkVmaps_print
  !========================================================================



  !========================================================================
  !========================================================================
  ! Procedures For ListPaths Structure
  !
  !========================================================================
  !========================================================================
  subroutine ListPaths_reset(this)
    ! Reset the Iterator to point at the first Path Link
    !-------------------------------------------------------
    class(ListPaths_t):: this
    this%curr => this%first
    return
  end subroutine ListPaths_reset
  !========================================================================
  subroutine ListPaths_iterate(this)
    ! Increment the Iterator to point at the next Path Link
    !----------------------------------------------------------
    class(ListPaths_t):: this
    this%curr => this%curr%get_Next()
    return
  end subroutine ListPaths_iterate
  !========================================================================
  function ListPaths_isFirst(this)
    ! Check if the Iterator points to the end of the list.
    !----------------------------------------------------
    class(ListPaths_t)        :: this
    logical                   :: ListPaths_isFirst
    class(LinkPaths_t),pointer:: lastPointer
    if(associated(this%curr)) then
      lastPointer => this%curr%get_Prev()
      ListPaths_isFirst = .not.associated(lastPointer)
    else
      ListPaths_isFirst = .true.
    endif
    return
  end function ListPaths_isFirst
  !========================================================================
  function ListPaths_isLast(this)
    ! Check if the Iterator points to the end of the list.
    !----------------------------------------------------
    class(ListPaths_t)        :: this
    logical                   :: ListPaths_isLast
    class(LinkPaths_t),pointer:: nextPointer
    if(associated(this%curr)) then
      nextPointer => this%curr%get_Next()
      ListPaths_isLast = .not.associated(nextPointer)
    else
      ListPaths_isLast = .true.
    endif
    return
  end function ListPaths_isLast
  !========================================================================
  function ListPaths_addLink(this,Path)
    ! Add a new Path Dataset link to the end of the list.
    !========================================================
    class(ListPaths_t)  :: this
    type(Path_t),pointer:: Path
    integer             :: ListPaths_addLink
    !
    ! Local Values
    !--------------
    class(LinkPaths_t),pointer:: newLink

    ! Allocate and a new Path Link and initialize 
    ! the pointer with Path from the input.
    !----------------------------------------------------
    allocate(newLink)
    newLink%Path => Path

    ! Increment the Counter to assign a unique index 
    ! to this new Link
    !-----------------------------------------------
    this%Counter    = this%Counter + 1
    newLink%idxPath = this%Counter

    ! Add it to the list
    !--------------------
    if(.not.associated(this%last)) then
      ! If the current list is empty, initialize 
      ! the list with this Link
      !-------------------------------------------
      newLink%prev => null()
      newLink%next => null()
      this%first   => newLink
      this%last    => newLink
      this%curr    => newLink
    else
      ! Append this Link to the end of the current list.
      !----------------------------------------------------
      newLink%prev   => this%last
      newLink%next   => null()
      this%last%next => newLink
      this%last      => newLink
    endif

    ! Return the DataSet ID to the calling routine
    !-----------------------------------------------
    ListPaths_addLink = newLink%idxPath

    return
  end function ListPaths_addLink
  !========================================================================
  subroutine ListPaths_rmLink(this,idxPath)
    ! Remove Path Dataset link with given ID
    !========================================================
    class(ListPaths_t):: this
    integer           :: idxPath
    !
    ! Local Values
    !--------------
    class(LinkPaths_t),pointer:: lastPointer
    class(LinkPaths_t),pointer:: nextPointer
    type(Path_t)      ,pointer:: curr_Path

    ! cycle thru the list to find the Link with matching ID
    !----------------------------------------------------------
    call this%reset
    do while(this%curr%get_ID() .ne. idxPath)
      call this%iterate
      if(.not.associated(this%curr)) then
        ! ERROR: idxPath not found in list
        !-------------------------------------
        write(iulog, *) 'ListPaths ERROR: No ',idxPath,'index in list'
        call endrun('ERROR: ListPaths Missing Grid Index')
      endif
    end do

    ! Now we want to remove the current 
    ! Link from the list
    !------------------------------------
    if(this%isFirst().and.this%isLast()) then
      ! We are removing the ONLY link in the List
      !--------------------------------------------
      curr_Path => this%curr%get_Data()
      
      ! Deallocate curr Link
      !----------------------
      if(allocated(curr_Path%idxHmap)) deallocate(curr_Path%idxHmap)
      if(allocated(curr_Path%idxVmap)) deallocate(curr_Path%idxVmap)
      deallocate(curr_Path)
      deallocate(this%curr)
   
      ! Nullify all pointers for an empty list
      !----------------------------------------
      this%first => null()
      this%curr  => null()
      this%last  => null()
    elseif(this%isFirst()) then
      ! We are removing the first Link in the list
      !-----------------------------------------------
      lastPointer => null()
      nextPointer => this%curr%get_Next()
      call nextPointer%set_Prev(lastPointer)

      ! Deallocate curr Link
      !----------------------
      curr_Path => this%curr%get_Data()
      if(allocated(curr_Path%idxHmap)) deallocate(curr_Path%idxHmap)
      if(allocated(curr_Path%idxVmap)) deallocate(curr_Path%idxVmap)
      deallocate(curr_Path)
      deallocate(this%curr)

      ! Reset fisrt and curr pointers
      !-----------------------------------
      this%first => nextPointer
      this%curr  => nextPointer

      nullify(lastPointer)
      nullify(nextPointer)
    elseif(this%isLast()) then
      ! We are removing the last Link in the list
      !-----------------------------------------------
      lastPointer => this%curr%get_Prev()
      nextPointer => null()
      call lastPointer%set_Next(nextPointer)

      ! Deallocate curr
      !--------------------
      curr_Path => this%curr%get_Data()
      if(allocated(curr_Path%idxHmap)) deallocate(curr_Path%idxHmap)
      if(allocated(curr_Path%idxVmap)) deallocate(curr_Path%idxVmap)
      deallocate(curr_Path)
      deallocate(this%curr)

      ! Reset last and curr pointers
      !-----------------------------------
      this%last => lastPointer
      this%curr => lastPointer

      nullify(lastPointer)
      nullify(nextPointer)
    else
      ! We are removing a Link in the middle of the list
      !-----------------------------------------------------
      lastPointer => this%curr%get_Prev()
      nextPointer => this%curr%get_Next()
      call lastPointer%set_Next(nextPointer)
      call nextPointer%set_Prev(lastPointer)

      ! Deallocate curr
      !--------------------
      curr_Path => this%curr%get_Data()
      if(allocated(curr_Path%idxHmap)) deallocate(curr_Path%idxHmap)
      if(allocated(curr_Path%idxVmap)) deallocate(curr_Path%idxVmap)
      deallocate(curr_Path)
      deallocate(this%curr)

      ! Reset curr pointer
      !-----------------------------------
      this%curr => lastPointer

      nullify(lastPointer)
      nullify(nextPointer)
    endif

    return
  end subroutine ListPaths_rmLink
  !========================================================================
  function ListPaths_getLinkID(this)
    ! Return the ID idxPath for curr in List
    !---------------------------------------------
    class(ListPaths_t):: this
    integer           :: ListPaths_getLinkID
    ListPaths_getLinkID = this%curr%get_ID()
    return
  end function ListPaths_getLinkID
  !========================================================================
  function ListPaths_getLinkData(this,idxPath)
    ! Given the Link ID, return a pointer to 
    ! the corresponding Path data structure.
    !---------------------------------------------
    class(ListPaths_t)  :: this
    integer             :: idxPath
    type(Path_t),pointer:: ListPaths_getLinkData

    ! cycle thru the list to find the Link with matching ID
    !----------------------------------------------------------
    call this%reset
    do while(this%curr%get_ID() .ne. idxPath)
      call this%iterate
      if(.not.associated(this%curr)) then
        ! ERROR: idxPath not found in list
        !-------------------------------------
        write(iulog, *) 'ListPaths ERROR: No ',idxPath,'index in list'
        call endrun('ERROR: ListPaths Missing Grid Index')
      endif
    end do
    ListPaths_getLinkData => this%curr%get_Data()
    return
  end function ListPaths_getLinkData
  !========================================================================
  function ListPaths_setLinkID(this,new_Path,newLink)
    ! Return the ID idxPath for new_Path. 
    ! If the grid exists on the list, then return the corresponding ID, if 
    ! this is a new grid, add it to the end of the list and return the new ID.
    ! Set newLink to indicate if this is a new or existing link
    !-------------------------------------------------------------------------
    class(ListPaths_t)  :: this
    type(Path_t),pointer:: new_Path
    integer             :: ListPaths_setLinkID
    logical             :: newLink
    ! 
    ! Local Values
    !-------------------
    logical:: link_found

    ! Cycle thru the list and see if there is already an  
    ! existing grid the matches the new Path
    !----------------------------------------------------
    call this%reset
    link_found = .false.
    do while(associated(this%curr))
      if(this%curr%isEqual(new_Path)) then
        link_found = .true.
        exit
      endif
      call this%iterate 
    end do

    if(link_found) then
      ! if the Grid exists, pass back the the corresponding ID. 
      !-----------------------------------------------------------------
      newLink = .false.
      ListPaths_setLinkID = this%curr%get_ID()
    else
      ! Otherwise add the new_Path to the end of the 
      ! list and return the new ID
      !--------------------------------------------------
      newLink = .true.
      ListPaths_setLinkID = this%addLink(new_Path)
    endif

    return
  end function ListPaths_setLinkID
  !========================================================================
  subroutine ListPaths_print(this)
    ! Cycle thru the current list and print 
    ! out the meta data for each link 
    !---------------------------------------------
    class(ListPaths_t):: this
    !
    ! Local values
    !-------------
    logical:: empty_list

    ! Check to see if the list is empty
    !-----------------------------------
    call this%reset
    empty_list = .false.
    if(this%isFirst().and.this%isLast()) then
      if(.not.associated(this%curr)) then
        empty_list = .true.
      endif
    endif
    
    if(empty_list) then
      ! Print out "Empty List" message here
      !---------------------------------------
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    LIST of Registered DataSet Path structures:'
      write(iulog,*) '          EMPTY LIST  '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) ' '
    else
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) '    LIST of Registered DataSet Path structures:'
      write(iulog,*) ' '

      ! cycle thru the list
      !----------------------
      do while(associated(this%curr))
        ! Issue a print command for each link
        !-------------------------------------
        call this%curr%printLink
        call this%iterate
      end do
      write(iulog,*) ' '
      write(iulog,*) '    ============================================================== '
      write(iulog,*) ' '
    endif

    return
  end subroutine ListPaths_print
  !========================================================================



  !========================================================================
  !========================================================================
  ! Procedures For LinkPaths Structure
  !
  !========================================================================
  !========================================================================
  function LinkPaths_getID(this)
    ! Return the ID for this Path Link
    !----------------------------------------------------
    class(LinkPaths_t):: this
    integer           :: LinkPaths_getID
    LinkPaths_getID = this%idxPath
    return
  end function LinkPaths_getID
  !========================================================================
  function LinkPaths_getData(this)
    ! Return the pointer to the next Path Link
    !----------------------------------------------------
    class(LinkPaths_t)  :: this
    type(Path_t),pointer:: LinkPaths_getData
    LinkPaths_getData => this%Path
    return
  end function LinkPaths_getData
  !========================================================================
  function LinkPaths_getNext(this)
    ! Return the pointer to the next Path Link
    !----------------------------------------------------
    class(LinkPaths_t)        :: this
    class(LinkPaths_t),pointer:: LinkPaths_getNext
    LinkPaths_getNext => this%next
    return
  end function LinkPaths_getNext
  !========================================================================
  subroutine LinkPaths_setNext(this,next)
    ! Set the pointer for the next Path Link
    !----------------------------------------------------
    class(LinkPaths_t)        :: this
    class(LinkPaths_t),pointer:: next
    this%next => next
    return
  end subroutine LinkPaths_setNext
  !========================================================================
  function LinkPaths_getPrev(this)
    ! Return the pointer to the previous Path Link
    !----------------------------------------------------
    class(LinkPaths_t)        :: this
    class(LinkPaths_t),pointer:: LinkPaths_getPrev
    LinkPaths_getPrev => this%prev
    return
  end function LinkPaths_getPrev
  !========================================================================
  subroutine LinkPaths_setPrev(this,prev)
    ! Return the pointer to the previous Hgrid_d Link
    !----------------------------------------------------
    class(LinkPaths_t)        :: this
    class(LinkPaths_t),pointer:: prev
    this%prev => prev
    return
  end subroutine LinkPaths_setPrev
  !========================================================================
  function LinkPaths_isEqual(this,new_Path)
    ! Logical function which returns True if the given 
    ! Path structure data matches those of the current Link.
    ! Other wise the structures differ and False is returned
    !----------------------------------------------------
    class(LinkPaths_t)  :: this
    type(Path_t),pointer:: new_Path
    logical             :: LinkPaths_isEqual
    !
    ! Local values
    !-----------------
    type(Path_t),pointer:: My_Path
    real(r8)            :: diffSum
    integer             :: ii

    ! Get a pointer to the current data structure
    !--------------------------------------------
    My_Path => this%get_Data()

    ! Set the restult to true, then scan the structure for differences, 
    ! if/when a difference is found, set the result to false and return 
    ! immediately.
    !-------------------------------------------------------------------
    LinkPaths_isEqual = .true.

    if(My_Path%ngrid.ne.new_Path%ngrid) then
      LinkPaths_isEqual = .false.
      return
    endif

    do ii=1,My_Path%ngrid
      if(My_Path%idxHmap(ii).ne.new_Path%idxHmap(ii)) then
        LinkPaths_isEqual = .false.
        return
      endif
      if(My_Path%idxVmap(ii).ne.new_Path%idxVmap(ii)) then
        LinkPaths_isEqual = .false.
        return
      endif
    end do

    return
  end function LinkPaths_isEqual
  !========================================================================
  subroutine LinkPaths_print(this)
    ! Print out the Path metadata for the current link
    !---------------------------------------------
    class(LinkPaths_t):: this
    !
    ! Local values
    !---------------
    integer:: ii

    ! Print out identifying values
    !--------------------------------
    write(iulog,*) '      -------------------------------------------'
    write(iulog,*) '      idxPath=',this%idxPath
    write(iulog,*) '                     ngrid= ',this%Path%ngrid
    if(this%Path%ngrid.le.0) then
      write(iulog,*) '                 grid: idxHmap, idxVmap= ** NO PATHS REGISTERED **'
    else
      do ii=1,this%Path%ngrid
        write(iulog,*) '                 grid: idxHmap, idxVmap= ',ii,this%Path%idxHmap(ii), &
                                                                      this%Path%idxVmap(ii)
      end do 
    endif
    write(iulog,*) '      -------------------------------------------'
    write(iulog,*) ' '

    return
  end subroutine LinkPaths_print
  !========================================================================

end module Dataconnector_mod
