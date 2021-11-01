module Dataset_mod
!===========================================================================
!
! MODULE: Dataset_mod
!
! DESCRIPTION: Module to open,close,and manage netCDF datasets. 
!              Internally a linked list of open datasets is maintained. 
!              For each newly opened dataset, a metadata structure describing 
!              the dimensions, variables,and grids in the dataset is initialized
!              is assigend a unique ID index and stored in a linked list of 
!              opened sets. When closed, the metadata structure is deallocated 
!              and the entry is removed from the list.
!              
!
!
! Document the interface:
!
!  public :: open_dataset          
!  public :: change_datafile       
!  public :: check_metadata        
!  public :: close_dataset         
!  public :: inq_dataset         
!  public :: inq_dataset_status
!  public :: inq_dataset_ncid
!  public :: inq_dataset_var
!  public :: inq_dataset_varid
!  public :: inq_dataset_dim
!  public :: inq_dataset_timerecord
!  public :: inq_dataset_step
!  private:: init_dataset_metadata
!  private:: update_dataset_metadata     
!  public :: print_dataset_metadata     
!  private:: set_dataset_time     
!  private:: greg2jday
!        
!
!===========================================================================
  ! Useful Modules
  !---------------------
  use netcdf                    ! PFC DIAG Temp to test PIO ERROR with nf90
  use ESMF
  use pio,           only: file_desc_t, iosystem_desc_t,PIO_NOERR, pio_inquire, pio_double
  use pio,           only: pio_inq_dimid, pio_inq_dimlen, pio_inq_dimname, var_desc_t
  use pio,           only: pio_inq_varndims, pio_inq_varname, pio_inq_vardimid
  use pio,           only: pio_get_var,pio_inq_varid, pio_openfile, pio_closefile
  use pio,           only: pio_inq_att,pio_get_att,PIO_RETURN_ERROR,pio_inq_varnatts,pio_inq_attname
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use cam_abortutils,only: endrun
  use spmd_utils,    only: masterproc    ! NEEDED?
  use cam_logfile,   only: iulog         ! NEEDED?

  ! Set all Global values and routines to private by default 
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private
  save

  public :: open_dataset          
  private:: open_dataset_dsid
  private:: open_dataset_dset
  public :: close_dataset         
  private:: close_dataset_dsid
  private:: close_dataset_dset
  public :: change_datafile       
  private:: change_datafile_dsid
  private:: change_datafile_dset
  public :: check_metadata        
  public :: inq_dataset         
  public :: inq_dataset_status
  public :: inq_dataset_ncid
  public :: inq_dataset_ncid_nf90   ! PFC DIAG Temp function to test PIO ERROR with nf90
  public :: inq_dataset_pio
  public :: inq_dataset_var
  public :: inq_dataset_varid
  public :: inq_dataset_dim
  public :: inq_dataset_timerecord
  public :: inq_dataset_step
  private:: init_dataset_metadata
  private:: update_dataset_metadata     
  public :: print_dataset_metadata
  private:: set_dataset_time     
  private:: greg2jday

  public :: dataset_t
  public :: dataset_dim_t
  public :: dataset_dimvar_t
  public :: dataset_var_t
  public :: dataset_grid_t
  public :: dataset_sclr_t

  ! Interfaces
  !-----------------
  interface open_dataset          
    procedure:: open_dataset_dsid
    procedure:: open_dataset_dset
  end interface

  interface close_dataset          
    procedure:: close_dataset_dsid
    procedure:: close_dataset_dset
  end interface

  interface change_datafile        
    procedure:: change_datafile_dsid
    procedure:: change_datafile_dset
  end interface

  ! Type Definitions
  !-------------------
  type OpenDataSets_t
    ! Linked List of currently Open datasets
    !--------------------------------------
      private
      integer                     :: DataSetID
      type(dataset_t)     ,pointer:: DataSet => null()
      type(OpenDataSets_t),pointer:: next    => null()
      type(OpenDataSets_t),pointer:: prev    => null()
    contains
      procedure:: get_OpenDataSet
      procedure:: set_OpenDataSet
      procedure:: get_DataSetID
      procedure:: set_DataSetID
      procedure:: get_NextDataSet
      procedure:: set_NextDataSet
      procedure:: get_PrevDataSet
      procedure:: set_PrevDataSet
  end type OpenDataSets_t

  type DataSetList_t
    ! Manage List of currently Open datasets
    !-----------------------------------------
      private
      class(OpenDataSets_t),pointer:: firstDS => null()
      class(OpenDataSets_t),pointer:: lastDS  => null()
      class(OpenDataSets_t),pointer:: currDS  => null()
    contains
      procedure:: iterate
      procedure:: reset
      procedure:: isFirst
      procedure:: isLast
      procedure:: add_DataSet
      procedure:: remove_DataSet
      procedure:: get_DataSet
      procedure:: get_DataSet_status
      procedure:: get_DataSet_ncid
      procedure:: get_DataSet_ncid_nf90 ! PFC DIAG Temp function to test PIO ERROR with nf90
      procedure:: get_DataSet_pio
      procedure:: get_DataSet_varid
      procedure:: get_DataSet_timerecord
      procedure:: get_DataSet_timestep
      procedure:: print_DataSet
  end type DataSetList_t

  type dataset_dimvar_t
    ! variables associated with 
    ! a dimension in DataSet 
    !----------------------------
    integer             :: varID   = -1
    type(var_desc_t)    :: VarDesc
    character(len=80)   :: name
    real(r8),allocatable:: val(:)
  end type dataset_dimvar_t

  type dataset_dim_t
    ! dimensions in DataSet 
    !-----------------------------------
    integer                           :: dimid     = -1
    character(len=80)                 :: name
    character(len=256)                :: long_name
    character(len=80)                 :: units
    integer                           :: len       = 0
    logical                           :: hasCoord  = .false.
    integer                           :: idxCoord  = -1
    integer                           :: ndvar     = 0
    type(dataset_dimvar_t),allocatable:: var(:)
  end type dataset_dim_t

  type dataset_var_t
    ! variables in DataSet 
    !--------------------
    integer            :: varID       = -1
    type(var_desc_t)   :: VarDesc
    character(len=80)  :: name
    character(len=256) :: long_name
    character(len=80)  :: units
    real(r8)           :: fill_value 
    integer            :: nvdim      = 0
    integer,allocatable:: dimID(:)
    integer            :: gridID     = -1
  end type dataset_var_t

  type dataset_sclr_t
    ! scalars in DataSet 
    !-------------------
    integer           :: varID      = -1
    type(var_desc_t)  :: VarDesc
    character(len=80) :: name
    character(len=256):: long_name
    real(r8)          :: val
  end type dataset_sclr_t

  type dataset_grid_t
    ! grids in DataSet 
    !--------------------
    integer            :: gridID    = -1
    character(len=80)  :: name
    character(len=256) :: long_name
    logical            :: hasTime   = .false.
    integer            :: nsdim     = 0
    integer,allocatable:: dimID(:)
    integer            :: timeID    = -1
  end type dataset_grid_t

  type dataset_t
    ! Metadata structure for a DataSet
    !----------------------------------
    character(len=256)              :: FilePath
    character(len=80 )              :: FileName
    type(iosystem_desc_t),pointer   :: pio_subsystem => null()
    integer                         :: pio_iotype
    type(file_desc_t)               :: ncid 
    integer                         :: ncid_nf90         !PFC TEMPORARY DIAG duplicate to test PIO error
    integer                         :: ulimID   = -1
    integer                         :: timeID   = -1
    integer                         :: ndims    = 0
    integer                         :: nvars    = 0
    integer                         :: nscalars = 0
    integer                         :: ngrids   = 0
    integer                         :: nsdimMAX = 0
    integer                         :: ntime    = 0
    integer                         :: natts    = 0
    type(dataset_dim_t ),allocatable:: dims   (:)
    type(dataset_var_t ),allocatable:: vars   (:)
    type(dataset_sclr_t),allocatable:: scalars(:)
    type(dataset_grid_t),allocatable:: grids  (:)
    real(r8)            ,allocatable:: time   (:)
    integer             ,allocatable:: date   (:)
    integer             ,allocatable:: datesec(:)
  end type dataset_t

  ! Global Values
  !-----------------
  type(DataSetList_t):: DSlist             ! List to keep track of opened Datasets.
  integer            :: DataSet_Count = 0  ! Counter to assign a unique ID to each dataset

contains

  !================================================================
  subroutine open_dataset_dsid(DataSetID,FilePath,FileName)
    ! 
    ! open_dataset: Given a file name and file path, open it 
    !               and initialize the DataSet structure with 
    !               the contents.
    !==============================================================
    !
    ! Passed Variables
    !------------------
    character(len=*),intent(in ):: FilePath
    character(len=*),intent(in ):: FileName
    integer         ,intent(out):: DataSetID
    !
    ! Local Values
    !-----------------
    type(dataset_t) ,pointer:: DataSet

    allocate(DataSet)
    call open_dataset_dset(FilePath,FileName,DataSet)
    DataSetID = DSlist%add_DataSet(DataSet)

    ! End Routine
    !------------------
    return
  end subroutine open_dataset_dsid
  !=====================================================================


  !================================================================
  subroutine open_dataset_dset(FilePath,FileName,DataSet)
    ! 
    ! open_dataset: Given a file name and file path, open it 
    !               and initialize the DataSet structure with 
    !               the contents.
    !==============================================================
    use shr_pio_mod, only: shr_pio_getiosys, shr_pio_getiotype
    use cam_instance,only: atm_id
    !
    ! Passed Variables
    !------------------
    character(len=*),intent(in   ):: FilePath
    character(len=*),intent(in   ):: FileName
    type(dataset_t) ,intent(inout):: DataSet
    !
    ! Local Values
    !-----------------
!    character(len=*):: fname
    character(len=256):: fname
    integer         :: mode
    integer         :: ret

    ! The DataSet must not contain an opened file. 
    ! Ensure that this is true by closing it first. 
    ! (No effect on unopened DataSet structure)
    !-------------------------------------------------
    call close_dataset_dset(DataSet)

    ! Open the Given netCDF file for reading
    !----------------------------------------
    DataSet%pio_subsystem => shr_pio_getiosys (atm_id)
    DataSet%pio_iotype    =  shr_pio_getiotype(atm_id)
    DataSet%FilePath      =  trim(FilePath)
    DataSet%FileName      =  trim(FileName)


    fname = trim(DataSet%FilePath)//'/'//trim(DataSet%FileName)
    mode     = 0

    ret = pio_openfile(DataSet%pio_subsystem,DataSet%ncid,    &
                       DataSet%pio_iotype   ,trim(fname) ,mode)

    if(ret.ne.PIO_NOERR) then
      call endrun('ERROR: open_dataset() Opening File='//trim(filename))
    endif

!PFC DIAG Temporary to test PIO ERROR with duplicate nf90 calls
    ret = nf90_open(trim(fname),NF90_NOWRITE,DataSet%ncid_nf90)
!PFC DIAG

    ! Initialize metadata values from the file
    !------------------------------------------
    call init_dataset_metadata(DataSet)
    
    if(masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'open_dataset: FilePath=',trim(DataSet%FilePath)
      write(iulog,*) 'open_dataset: FileName=',trim(DataSet%FileName)
      write(iulog,*) ' '
    endif
    
    ! End Routine
    !------------------
    return
  end subroutine open_dataset_dset
  !=====================================================================


  !================================================================
  subroutine close_dataset_dsid(DataSetID)
    ! 
    ! close_dataset: Deallocate the DataSet struct and 
    !                close the file.
    !   Deallocate space before closing the dataset
    !==============================================================
    !
    ! Passed Variables
    !------------------
    integer,intent(inout):: DataSetID

    call DSlist%remove_DataSet(DataSetID)
    DataSetID = -1

    ! End Routine
    !------------------
    return
  end subroutine close_dataset_dsid
  !=====================================================================


  !================================================================
  subroutine close_dataset_dset(DataSet)
    ! 
    ! close_dataset: Deallocate the DataSet struct and 
    !                close the file.
    !   Deallocate space before closing the dataset
    !==============================================================
    !
    ! Passed Variables
    !------------------
    type(dataset_t),intent(inout):: DataSet
    !
    ! Local Values
    !------------------
    integer:: ii,nn

    ! Close it up
    !------------------
    if(DataSet%ndims.gt.0) call pio_closefile(DataSet%ncid)

!PFC DIAG Temporary to test PIO ERROR with duplicate nf90 calls
    if(DataSet%ndims.gt.0) then
      nn = nf90_close(DataSet%ncid_nf90)
    endif
!PFC DIAG

    ! don't need to point at the pio_subsystem anymore
    !---------------------------------------------------
    if(associated(DataSet%pio_subsystem)) nullify(DataSet%pio_subsystem)

    ! Deallocate dims() data structure
    !-----------------------------------
    if(allocated(DataSet%dims)) then
      do nn=1,DataSet%ndims
        if(allocated(DataSet%dims(nn)%var)) then
          do ii=1,DataSet%dims(nn)%ndvar
            if(allocated(DataSet%dims(nn)%var(ii)%val)) then
              deallocate(DataSet%dims(nn)%var(ii)%val)
            endif
          end do
          deallocate(DataSet%dims(nn)%var)
        endif
      end do 
      deallocate(DataSet%dims)
    endif
    DataSet%ndims = 0

    ! Deallocate vars() data structure
    !-----------------------------------
    if(allocated(DataSet%vars)) then
      do ii=1,DataSet%nvars
        if(allocated(DataSet%vars(ii)%dimID)) deallocate(DataSet%vars(ii)%dimID)
      end do
      deallocate(DataSet%vars)
    endif

    ! Deallocate scalars() data structure
    !-----------------------------------
    if(allocated(DataSet%scalars)) deallocate(DataSet%scalars)

    ! Deallocate grids() data structure
    !-----------------------------------
    if(allocated(DataSet%grids)) then
      do ii=1,DataSet%ngrids
        if(allocated(DataSet%grids(ii)%dimID)) deallocate(DataSet%grids(ii)%dimID)
      end do
      deallocate(DataSet%grids)
    endif

    ! Deallocate time arrays
    !------------------------
    if(allocated(DataSet%time   )) deallocate(DataSet%time   )
    if(allocated(DataSet%date   )) deallocate(DataSet%date   )
    if(allocated(DataSet%datesec)) deallocate(DataSet%datesec)

    ! End Routine
    !------------------
    return
  end subroutine close_dataset_dset
  !=====================================================================


  !================================================================
  subroutine change_datafile_dsid(DataSetID,FileName)
    ! 
    ! change_datafile_dsid: Given a file name and a DataSetID, get a 
    !                       pointer to the DataSet structure and then 
    !                       change the current open file and update the 
    !                       file-dependent metadata. It is expected that 
    !                       the variable and grid structures are unchanged.
    !==============================================================
    !
    ! Passed Variables
    !------------------
    integer         ,intent(in):: DataSetID
    character(len=*),intent(in):: FileName
    !
    ! Local Values
    !--------------
    type(dataset_t),pointer:: DataSet

    DataSet => inq_dataset(DataSetID)
    call change_datafile_dset(DataSet,FileName)
    nullify(DataSet)

    ! End Routine
    !------------------
    return
  end subroutine change_datafile_dsid
  !=====================================================================


  !================================================================
  subroutine change_datafile_dset(DataSet,FileName)
    ! 
    ! change_datafile_dset: Assume that the contents in DataSet represent an 
    !                       active/open netCDF dataset. Then close, the current 
    !                       file, open a new one using the given FileName, and 
    !                       then update only the transient DataSet struct values 
    !                       that change from file to file. 
    !==============================================================
    !
    ! Passed Variables
    !------------------
    type(dataset_t) ,intent(inout):: DataSet
    character(len=*),intent(in   ):: FileName
    !
    ! Local Values
    !---------------
!    character(len=*):: fname
    character(len=256):: fname
    integer         :: mode
    integer         :: ret
    integer         :: ndims,nvars,natts,ulimID

    ! Verify that DataSet is initialized for an open netCDF file
    !-------------------------------------------------------------
    ret = pio_inquire(DataSet%ncid, nDimensions=ndims, &
                                     nVariables=nvars, &
                                    nAttributes=natts, &
                                 unlimitedDimID=ulimID )
    if(ret.ne.PIO_NOERR) then
      call endrun('ERROR: change_datafile() Called without an open dataset')
    endif
    if((DataSet%ndims .ne.ndims ).or.(DataSet%nvars .ne.nvars ).or. &
       (DataSet%natts .ne.natts ).or.(DataSet%ulimID.ne.ulimID)     ) then
      call endrun('ERROR: change_datafile() DataSet does not match file contents.')
    endif

    ! Close the current netCDF file
    !-------------------------------
    call pio_closefile(DataSet%ncid)

!PFC DIAG Temporary to test PIO ERROR with duplicate nf90 calls
    ret = nf90_close(DataSet%ncid_nf90)
!PFC DIAG

    ! Open the new file
    !---------------------
    DataSet%FileName      = trim(FileName)
    fname = trim(DataSet%FilePath)//'/'//trim(DataSet%FileName)
    mode     = 0

    ret = pio_openfile(DataSet%pio_subsystem,DataSet%ncid  , &
                       DataSet%pio_iotype   ,trim(fname),mode)

    if(ret.ne.PIO_NOERR) then
      call endrun('ERROR: change_datafile() Opening File='//trim(filename))
    endif
    
!PFC DIAG Temporary to test PIO ERROR with duplicate nf90 calls
    ret = nf90_open(trim(fname),NF90_NOWRITE,DataSet%ncid_nf90)
!PFC DIAG

    ! Do a quick check for a mis-match between DataSet 
    ! and the file contnents.
    !-------------------------------------------------------------
    ret = pio_inquire(DataSet%ncid, nDimensions=ndims, &
                                     nVariables=nvars, &
                                    nAttributes=natts, &
                                 unlimitedDimID=ulimID )
    if(ret.ne.PIO_NOERR) then
      call endrun('ERROR: change_datafile() pio_inquire error')
    endif
    if((DataSet%ndims .ne.ndims ).or.(DataSet%nvars .ne.nvars ).or. &
       (DataSet%natts .ne.natts ).or.(DataSet%ulimID.ne.ulimID)     ) then
      call endrun('ERROR: change_datafile() DataSet does not match file contents.')
    endif

    ! Update metadata values from the file
    !------------------------------------------
    call update_dataset_metadata(DataSet)

    ! End Routine
    !------------------
    return
  end subroutine change_datafile_dset
  !=====================================================================


  !================================================================
  subroutine check_metadata(DataSet)
    ! 
    ! check_metadata: Verify that the metadata in DataSet matches 
    !                       the values in the current file. 
    !==============================================================
    !
    ! Passed Variables
    !------------------
    type(dataset_t),intent(inout):: DataSet

    ! Do Some checks!
    !-------------------

    ! End Routine
    !------------------
    return
  end subroutine check_metadata
  !=====================================================================


  !================================================================
  function inq_dataset(DataSetID,found)
    ! 
    ! inq_dataset: Given the open dataset index, return
    !              a pointer to the DataSet structure
    !==============================================================
    !
    ! Passed Variables
    !------------------
    integer         ,intent(in ):: DataSetID
    logical,optional,intent(out):: found
    type(dataset_t),pointer:: inq_dataset
    logical:: inq_found

    inq_dataset => DSlist%get_DataSet(DataSetID,inq_found)
    if(present(found)) found = inq_found

    ! End Routine
    !------------------
    return
  end function inq_dataset
  !=====================================================================


  !================================================================
  function inq_dataset_status(DataSetID)
    ! 
    ! inq_dataset_status: Check to see if the given dataset index 
    !                     corresponds to an open DataSet.
    !==============================================================
    !
    ! Passed Variables
    !------------------
    integer,intent(in):: DataSetID
    logical            :: inq_dataset_status

    inq_dataset_status = DSlist%get_DataSet_status(DataSetID)

    ! End Routine
    !------------------
    return
  end function inq_dataset_status
  !=====================================================================


  !================================================================
  function inq_dataset_ncid(DataSetID,found)
    ! 
    ! inq_dataset_ncid: Given the open dataset index, return
    !                   the corresponding netcdf ID, ncid
    !==============================================================
    !
    ! Passed Variables
    !------------------
    integer         ,intent(in ):: DataSetID
    logical,optional,intent(out):: found
    type(file_desc_t) :: inq_dataset_ncid
    logical:: inq_found

    inq_dataset_ncid = DSlist%get_DataSet_ncid(DataSetID,inq_found)
    if(present(found)) found = inq_found

    ! End Routine
    !------------------
    return
  end function inq_dataset_ncid
  !=====================================================================

!PFC DIAG Temporary function for testing PIO ERROR with nf90
  !================================================================
  function inq_dataset_ncid_nf90(DataSetID,found)
    ! 
    ! inq_dataset_ncid: Given the open dataset index, return
    !                   the corresponding netcdf ID, ncid
    !==============================================================
    !
    ! Passed Variables
    !------------------
    integer         ,intent(in ):: DataSetID
    logical,optional,intent(out):: found
    integer                     :: inq_dataset_ncid_nf90
    logical:: inq_found

    inq_dataset_ncid_nf90 = DSlist%get_DataSet_ncid_nf90(DataSetID,inq_found)
    if(present(found)) found = inq_found

    ! End Routine
    !------------------
    return
  end function inq_dataset_ncid_nf90
  !=====================================================================
!PFC DIAG Temporary function for testing PIO ERROR with nf90


  !================================================================
  logical function inq_dataset_pio(DataSetID,pio_subsystem,pio_iotype)
    ! 
    ! inq_dataset_pio: Given the open dataset index, return
    !                  the corresponding pio information
    !==============================================================
    !
    ! Passed Variables
    !------------------
    integer                      ,intent(in ):: DataSetID
    type(iosystem_desc_t),pointer,intent(out):: pio_subsystem 
    integer                      ,intent(out):: pio_iotype
    !
    ! Local Values
    !--------------
    logical:: found
    
    pio_subsystem => DSlist%get_DataSet_pio(DataSetID,pio_iotype,found)

    if(found) then
      inq_dataset_pio = .true.
    else
      inq_dataset_pio = .false.
      pio_subsystem   => null()
      pio_iotype      = -1
    endif

    ! End Routine
    !------------------
    return
  end function inq_dataset_pio
  !=====================================================================


  !========================================================================
  logical function inq_dataset_var(I_dataset,var_name,var_index)
    !
    ! inq_dataset_var: Look thru the given dataset structure for a 
    !                  variable with the given name. If it is found, 
    !                  optionally return the variable index in the 
    !                  datastructure.
    !==============================================================
    implicit none
    !
    ! Passed Variables
    !-------------------
    type(dataset_t) ,intent(in) :: I_dataset
    character(len=*),intent(in) :: var_name
    integer,optional,intent(out):: var_index
    !
    ! Local Values
    !---------------
    integer:: ii

    ! Loop thru variable names looking for var_name
    !-----------------------------------------------
    inq_dataset_var = .false.

    if(I_dataset%nvars.gt.0) then
      do ii=1,I_dataset%nvars
        if(trim(I_dataset%vars(ii)%name).eq.trim(var_name)) then
          inq_dataset_var = .true.
          if(present(var_index)) var_index= ii
          exit
        endif
      end do
    endif

    ! End Routine
    !------------------
    return
  end function inq_dataset_var
  !=====================================================================


  !========================================================================
  logical function inq_dataset_varid(dataSetID,var_name,var_index,Var,grid_index)
    !
    ! inq_dataset_varid: Look the dataset structure for dataSetID, then for a 
    !                    variable with the given name. If it is found, 
    !                    return the variable index and the corresponding 
    !                    grid index in the datastructure.
    !==============================================================
    implicit none
    !
    ! Passed Variables
    !-------------------
    integer            ,intent(in) :: dataSetID
    character(len=*)   ,intent(in) :: var_name
    integer            ,intent(out):: var_index
    type(dataset_var_t),intent(out):: Var
    integer            ,intent(out):: grid_index
    !
    ! Local Values
    !---------------
    character(len=80):: My_varName

    My_varName = trim(var_name)
    inq_dataset_varid = DSlist%get_DataSet_varid(dataSetID,My_varName,   &
                                                 var_index,Var,grid_index)

    ! End Routine
    !------------------
    return
  end function inq_dataset_varid
  !=====================================================================


  !========================================================================
  logical function inq_dataset_dim(I_dataset,dim_name,dim_index)
    !
    ! inq_dataset_dim: Look thru the given dataset structure for a 
    !                  dimension with the given name. If it is found, 
    !                  optionally return the dim index in the 
    !                  datastructure.
    !==============================================================
    implicit none
    !
    ! Passed Variables
    !-------------------
    type(dataset_t) ,intent(in) :: I_dataset
    character(len=*),intent(in) :: dim_name
    integer,optional,intent(out):: dim_index
    !
    ! Local Values
    !---------------
    integer:: ii

    ! Loop thru variable names looking for var_name
    !-----------------------------------------------
    inq_dataset_dim = .false.

    if(I_dataset%ndims.gt.0) then
      do ii=1,I_dataset%ndims
        if(trim(I_dataset%dims(ii)%name).eq.trim(dim_name)) then
          inq_dataset_dim = .true.
          if(present(dim_index)) dim_index= ii
          exit
        endif
      end do
    endif

    ! End Routine
    !------------------
    return
  end function inq_dataset_dim
  !=====================================================================


  !========================================================================
  logical function inq_dataset_timerecord(dataSetID,date,datesec,time_index)
    !
    ! inq_dataset_timerecord: Look the dataset structure for dataSetID, then 
    !                         scan date/datesec values to find the given values.
    !                         if the record was found, then return its value
    !                         indicate success. If not, set the return function
    !                         to false and set the time_index to -1.
    !==============================================================
    implicit none
    !
    ! Passed Variables
    !-------------------
    integer,intent(in) :: dataSetID
    integer,intent(in) :: date
    integer,intent(in) :: datesec
    integer,intent(out):: time_index

    inq_dataset_timerecord = DSlist%get_DataSet_timerecord(dataSetID,date,datesec,time_index)

    ! End Routine
    !------------------
    return
  end function inq_dataset_timerecord
  !=====================================================================


  !========================================================================
  logical function inq_dataset_step(dataSetID,OBS_step)
    !
    ! inq_dataset_step: Look the dataset structure for dataSetID, then 
    !                   use date/datesec/time values to determine the 
    !                   time in seconds between time records.
    !                   if the number of time records is lt 2, then 
    !                   set OBS_step = -1 and return false.
    !==============================================================
    implicit none
    !
    ! Passed Variables
    !-------------------
    integer,intent(in) :: dataSetID
    integer,intent(out):: OBS_step

    inq_dataset_step = DSlist%get_DataSet_timestep(dataSetID,OBS_step)

    ! End Routine
    !------------------
    return
  end function inq_dataset_step
  !=====================================================================


  !================================================================
  subroutine init_dataset_metadata(DataSet)
    ! 
    ! init_dataset: Initialize the given datastructure with info about 
    !               the contents of the given netcdf file. 
    !
    !       Given:  It is expected that the DataSet struture will 
    !               have the following values set for a file that
    !               has been opened via open_dataset():
    !
    !                   DataSet%pio_subsystem
    !                   DataSet%pio_iotype
    !                   DataSet%FilePath
    !                   DataSet%FileName
    !                   DataSet%ncid 
    !
    !       This routine reads the contantes of the file and 
    !       initializes the remaining datastructure values needed 
    !       for inquiries about available grids and variables.
    !
    !==============================================================
    integer,parameter:: REFyear  = 0001
    integer,parameter:: REFmonth = 1
    integer,parameter:: REFday   = 1
    !
    ! Passed Variables
    !------------------
    type(dataset_t),intent(inout):: DataSet
    ! 
    ! Local Values
    !---------------
    integer ,allocatable:: date   (:)
    real(r8),allocatable:: datesec(:)
    integer ,allocatable:: varlist(:)
    logical ,allocatable:: Wrk_hasTime(:)
    integer ,allocatable:: Wrk_timeID (:)
    integer ,allocatable:: Wrk_nsdim  (:)
    integer ,allocatable:: Wrk_dimID  (:,:)
    integer ,allocatable:: g_dimID    (:)
    integer:: ret,nn,ii
    logical:: hasDate,hasDateSec
    integer:: tt,year,month,day
    integer:: ndvar,idxCoord
    integer:: dateID, datesecID
    integer:: nscalar
    integer:: Ngrids, Nknown_vars, idxSgrid
    logical:: newGrid,hasTime,g_hasTime
    integer:: nsdimMAX,ngridsMAX,g_nsdim,g_timeID,ndnum,ns,dd,idxGrid
    integer:: Wrk_ngrids
    character(len=80):: AttName
    integer:: natts,na
   
    ! Get the number of dimensions, variables, and ID for unlimited dimension
    !------------------------------------------------------------------------
    ret = pio_inquire(DataSet%ncid, nDimensions=DataSet%ndims, &
                                     nVariables=DataSet%nvars, &
                                    nAttributes=DataSet%natts, &
                                 unlimitedDimID=DataSet%ulimID )

    ! Read in the Dimension information
    !------------------------------------
    if(allocated(DataSet%dims)) deallocate(DataSet%dims)
    allocate(DataSet%dims(DataSet%ndims))
    do nn=1,DataSet%ndims
      DataSet%dims(nn)%dimid = nn
      ret = pio_inq_dimlen (DataSet%ncid,nn,DataSet%dims(nn)%len )
      ret = pio_inq_dimname(DataSet%ncid,nn,DataSet%dims(nn)%name)
    end do

    ! Look for specific dimensions of interest
    !----------------------------------------------------------
    hasTime = (pio_inq_dimid(DataSet%ncid,'time',DataSet%timeID).eq.PIO_NOERR) 
    if(.not.hasTime) DataSet%timeID = -1

    ! If present, Check that the ulimID is the same as timeID.
    !----------------------------------------------------------
    if(hasTime.and.(DataSet%timeID.ne.DataSet%ulimID)) then
      call endrun('ERROR: init_dataset_metadata() timeID not equal to ulimID?')
    endif

    ! Read in the Variable information
    !------------------------------------
    if(allocated(DataSet%vars)) deallocate(DataSet%vars)
    allocate(DataSet%vars(DataSet%nvars))
    nscalar = 0
    do ii=1,DataSet%nvars
      DataSet%vars(ii)%varID = ii

      ! Get the variable name, attributes, 
      ! and number of dimensions
      !-----------------------------------------
      ret = pio_inq_varname (DataSet%ncid,ii,DataSet%vars(ii)%name )
      ret = pio_inq_varndims(DataSet%ncid,ii,DataSet%vars(ii)%nvdim)
      ret = pio_inq_varid   (DataSet%ncid,DataSet%vars(ii)%name,DataSet%vars(ii)%VarDesc)

      ! Create the list of dimension ID's for this variable
      ! Count up scalar variables as they go by.
      !------------------------------------------------------
      if(DataSet%vars(ii)%nvdim.ge.1) then
        if(allocated(DataSet%vars(ii)%dimID)) deallocate(DataSet%vars(ii)%dimID)
        allocate(DataSet%vars(ii)%dimID(DataSet%vars(ii)%nvdim))
        ret = pio_inq_vardimid(DataSet%ncid,ii,                                &
                               DataSet%vars(ii)%dimID(1:DataSet%vars(ii)%nvdim))
      elseif(DataSet%vars(ii)%nvdim.eq.0) then
        nscalar = nscalar + 1
      else
        call endrun('ERROR: init_dataset_metadata() varaible with negative # of dims?')
      endif

      ! Set default variable attributes
      !---------------------------------
      DataSet%vars(ii)%long_name  = trim(DataSet%vars(ii)%name)
      DataSet%vars(ii)%units      = ' '
      DataSet%vars(ii)%fill_value = -9999._r8

      ! Get available variable attributes
      !-----------------------------------
      ret = pio_inq_varnatts(DataSet%ncid,ii,natts)
      do na=1,natts
        ret = pio_inq_attname(DataSet%ncid,ii,na,AttName)
        if(trim(AttName).eq.'long_name') then
          ret = pio_get_att(DataSet%ncid,ii, 'long_name',DataSet%vars(ii)%long_name)
        elseif(trim(AttName).eq.'units') then
          ret = pio_get_att(DataSet%ncid,ii, 'units'    ,DataSet%vars(ii)%units)
        elseif(trim(AttName).eq.'_FillValue') then
          ret = pio_get_att(DataSet%ncid,ii,'_FillValue',DataSet%vars(ii)%fill_value)
        endif
      end do ! na=1,natts
    end do

    ! Create a list of Scalar values
    !----------------------------------
    if(allocated(DataSet%scalars)) deallocate(DataSet%scalars)
    allocate(DataSet%scalars(nscalar))
    DataSet%nscalars = 0
    do ii=1,DataSet%nvars
      if(DataSet%vars(ii)%nvdim.eq.0) then
        DataSet%nscalars = DataSet%nscalars + 1
        DataSet%scalars(DataSet%nscalars)%varID      = ii
        DataSet%scalars(DataSet%nscalars)%name       = DataSet%vars(ii)%name
        DataSet%scalars(DataSet%nscalars)%VarDesc    = DataSet%vars(ii)%VarDesc
        DataSet%scalars(DataSet%nscalars)%long_name  = DataSet%vars(ii)%long_name
        ret = pio_get_var(DataSet%ncid, ii,                    &
                          DataSet%scalars(DataSet%nscalars)%val)
      endif
    end do
    if(nscalar.ne.DataSet%nscalars) then
      call endrun('ERROR: init_dataset_metadata() mis-match in the number of scalar vars?')
    endif

    ! For each dimension, get a list of 1D variables with that dimension. 
    !--------------------------------------------------------------------
    allocate(varlist(DataSet%nvars))
    do nn=1,DataSet%ndims
      DataSet%dims(nn)%hasCoord = .false.
      DataSet%dims(nn)%idxCoord = -1
      DataSet%dims(nn)%ndvar    = 0
      do ii=1,DataSet%nvars
        if((DataSet%vars(ii)%nvdim.eq.1).and.(DataSet%vars(ii)%dimID(1).eq.nn)) then
          DataSet%dims(nn)%ndvar = DataSet%dims(nn)%ndvar + 1
          varlist(DataSet%dims(nn)%ndvar) = DataSet%vars(ii)%varID
          if(trim(DataSet%vars(ii)%name).eq.trim(DataSet%dims(nn)%name)) then
            DataSet%dims(nn)%hasCoord = .true.
            DataSet%dims(nn)%idxCoord = DataSet%vars(ii)%varID
          endif
        endif
      end do

      ! Now go thru and read in 1D variables associated with each dimension.
      !----------------------------------------------------------------------
      if(DataSet%dims(nn)%ndvar.gt.0) then
        if(allocated(DataSet%dims(nn)%var)) deallocate(DataSet%dims(nn)%var)
        allocate(DataSet%dims(nn)%var(DataSet%dims(nn)%ndvar))
        do ii=1,DataSet%dims(nn)%ndvar
          DataSet%dims(nn)%var(ii)%varID   = varlist(ii)
          DataSet%dims(nn)%var(ii)%name    = DataSet%vars(varlist(ii))%name
          DataSet%dims(nn)%var(ii)%VarDesc = DataSet%vars(varlist(ii))%VarDesc
          if(allocated(DataSet%dims(nn)%var(ii)%val)) deallocate(DataSet%dims(nn)%var(ii)%val)
          allocate(DataSet%dims(nn)%var(ii)%val(DataSet%dims(nn)%len))
          ret = pio_get_var(DataSet%ncid,varlist(ii),   &
                            DataSet%dims(nn)%var(ii)%val)
        end do
      endif
    end do ! nn=1,DataSet%ndims
    deallocate(varlist)

    ! Determine the largest spatial dimension of variables, and 
    ! the maximum number of unique grids that can be present
    !----------------------------------------------------------------
    nsdimMAX = 0
    do ii=1,DataSet%nvars
      if(DataSet%vars(ii)%nvdim.ge.2) then
        if(DataSet%vars(ii)%dimID(DataSet%vars(ii)%nvdim).eq.DataSet%timeID) then
          g_nsdim = DataSet%vars(ii)%nvdim - 1
        else
          g_nsdim = DataSet%vars(ii)%nvdim
        endif
        if(g_nsdim.gt.nsdimMAX) nsdimMAX = g_nsdim
      endif
    end do
    if(hasTime) then
      ndnum = DataSet%ndims-1
    else
      ndnum = DataSet%ndims
    endif

    ngridsMAX = ndnum
    do ii=1,(nsdimMAX-1)
      ngridsMAX = ngridsMAX*(ndnum-ii)
    end do

    ! Create works space, then go thru variables and 
    ! tabulate the unique grids present in the file
    !------------------------------------------------------
    allocate(Wrk_hasTime(ngridsMAX         ))
    allocate(Wrk_timeID (ngridsMAX         ))
    allocate(Wrk_nsdim  (ngridsMAX         ))
    allocate(Wrk_dimID  (nsdimMAX,ngridsMAX))
    allocate(  g_dimID  (nsdimMAX          ))

    Wrk_dimID(:,:) = -1
    Wrk_ngrids     = 0
    do ii=1,DataSet%nvars
      if(DataSet%vars(ii)%nvdim.ge.2) then

        ! Get grid values for the current variable
        !-------------------------------------------
        if(DataSet%vars(ii)%dimID( DataSet%vars(ii)%nvdim ).eq. DataSet%timeID ) then
          g_nsdim   = DataSet%vars(ii)%nvdim - 1
          g_hasTime = .true.
          g_timeID  = DataSet%timeID
        else
          g_nsdim = DataSet%vars(ii)%nvdim
          g_hasTime = .false.
          g_timeID  = -1
        endif
        g_dimID(:) = -1
        g_dimID(1:g_nsdim) = DataSet%vars(ii)%dimID(1:g_nsdim)

        ! Check if the grid is in the existing list
        !-------------------------------------------
        newGrid = .true.
        do nn=1,Wrk_ngrids
          if((Wrk_hasTime(nn).eqv.g_hasTime).and. &
             (Wrk_nsdim  (nn).eq. g_nsdim  )      ) then
            ns=0
            do dd=1,g_nsdim
              if(Wrk_dimID(dd,nn).eq.g_dimID(dd)) ns=ns+1
            end do
            if(ns.eq.g_nsdim) then
              newGrid = .false.
              idxGrid = nn
            endif
          endif
          if(.not.newGrid) exit
        end do
    
        ! Add new grids to the list as they occur, 
        ! Assign the unique grid ID for this variable
        !----------------------------------------------
        if(newGrid) then
          Wrk_ngrids = Wrk_ngrids + 1
          Wrk_hasTime(Wrk_ngrids) = g_hasTime
          Wrk_timeID (Wrk_ngrids) = g_timeID
          Wrk_nsdim  (Wrk_ngrids) = g_nsdim
          Wrk_dimID(:,Wrk_ngrids) = g_dimID(:)
          DataSet%vars(ii)%gridID = Wrk_ngrids
        else
          DataSet%vars(ii)%gridID  = idxGrid
        endif
      endif
    end do

    ! Now save the set of unique grids for later use
    !------------------------------------------------
    DataSet%ngrids   = Wrk_ngrids
    DataSet%nsdimMAX = nsdimMAX
    if(allocated(DataSet%grids)) deallocate(DataSet%grids)
    allocate(DataSet%grids(DataSet%ngrids))
    do ii = 1,DataSet%ngrids
      if(allocated(DataSet%grids(ii)%dimID)) deallocate(DataSet%grids(ii)%dimID)
      allocate(DataSet%grids(ii)%dimID(DataSet%nsdimMAX))
      DataSet%grids(ii)%gridID   = ii
      DataSet%grids(ii)%hasTime  = Wrk_hasTime(ii)
      DataSet%grids(ii)%timeID   = Wrk_timeID (ii)
      DataSet%grids(ii)%nsdim    = Wrk_nsdim  (ii)
      DataSet%grids(ii)%dimID(:) = Wrk_dimID  (:,ii)
      DataSet%grids(ii)%name = '('
      do nn=1,DataSet%grids(ii)%nsdim
        DataSet%grids(ii)%name = trim(DataSet%grids(ii)%name)                        // &
                                 trim(DataSet%dims(DataSet%grids(ii)%dimID(nn))%name)
        if(nn.lt.DataSet%grids(ii)%nsdim)                          &
          DataSet%grids(ii)%name = trim(DataSet%grids(ii)%name)//','
      end do
      DataSet%grids(ii)%name = trim(DataSet%grids(ii)%name)//')'
    end do
    deallocate(Wrk_hasTime)
    deallocate(Wrk_timeID)
    deallocate(Wrk_nsdim)
    deallocate(Wrk_dimID)
    deallocate(  g_dimID)

    ! Set the time (unlimited dimension) values for this file
    !------------------------------------------------------------
    call set_dataset_time(DataSet)

    ! End Routine
    !------------------
    return
  end subroutine init_dataset_metadata
  !=====================================================================


  !================================================================
  subroutine update_dataset_metadata(DataSet)
    ! 
    ! update_dataset_metadata: Update transient DataSet struct values 
    !                          from current opened file. 
    !==============================================================
    !
    ! Passed Variables
    !------------------
    type(dataset_t),intent(inout):: DataSet

    ! Set time values from file contents
    !-----------------------------------------
    call set_dataset_time(DataSet)

    ! Update Attributes/other values that change file-to-file
    !----------------------------------------------------------

    ! End Routine
    !------------------
    return
  end subroutine update_dataset_metadata
  !=====================================================================


  !================================================================
  subroutine print_dataset_metadata(DataSetID)
    ! 
    ! print_dataset_metadata: Print the metadata values 
    !                         for the given DataSetID
    ! 
    !==============================================================
    !
    ! Passed Variables
    !------------------
    integer,intent(in):: DataSetID

    call DSlist%print_DataSet(DataSetID)

    ! End Routine
    !------------------
    return
  end subroutine print_dataset_metadata
  !=====================================================================


  !================================================================
  subroutine set_dataset_time(DataSet)
    ! 
    ! set_dataset_time: Update time values for DataSet struct 
    !                   from current opened file. 
    !==============================================================
    integer,parameter:: COORDINATE = 0
    integer,parameter:: DATESECREF = 1
    integer,parameter:: SCALAR     = 2
    integer,parameter:: TIMEINDEX  = 3
    integer,parameter:: NULLTIME   = 4

    integer,parameter:: TimePreference = DATESECREF
 !   integer,parameter:: TimePreference = COORDINATE
    integer,parameter:: REFyear  = 0001
    integer,parameter:: REFmonth = 1
    integer,parameter:: REFday   = 1
    !
    ! Passed Variables
    !------------------
    type(dataset_t),intent(inout):: DataSet
    !
    ! Local Values
    !---------------
    logical            :: hasUTCdate
    integer            :: utcdateID
    integer,allocatable:: UTC_date(:)

    integer:: dateID, datesecID, timeID_dim, timeID_var
    logical:: hasTimeDim, hasTimeVar, hasDate, hasDateSec, hasTimeScalar
    integer:: TimeOption
    integer:: tt, year, month, day, nvdim
    integer:: ret

    ! Check For time dimension and coordinate
    !------------------------------------------------
    ret = pio_inq_dimid(DataSet%ncid,'time',timeID_dim)
    if(ret.ne.PIO_NOERR) then 
      timeID_dim    = -1
      hasTimeDim    = .false.
      DataSet%ntime = 1
      DataSet%timeID = timeID_dim
      DataSet%ulimID = timeID_dim
    else
      hasTimeDim = .true.
      ret = pio_inq_dimlen(DataSet%ncid,timeID_dim, DataSet%ntime)
      DataSet%timeID = timeID_dim
      DataSet%ulimID = timeID_dim
    endif

    hasTimeVar = .false.
    hasDate    = .false.
    hasDateSec = .false.
    hasUTCdate = .false.
    do tt=1,DataSet%nvars
      if(trim(DataSet%vars(tt)%name).eq.'time'    ) hasTimeVar = .true.
      if(trim(DataSet%vars(tt)%name).eq.'date'    ) hasDate    = .true.
      if(trim(DataSet%vars(tt)%name).eq.'datesec' ) hasDateSec = .true.
      if(trim(DataSet%vars(tt)%name).eq.'utc_date') hasUTCdate = .true.
    end do
    if(hasTimeVar) ret = pio_inq_varid(DataSet%ncid,'time'    ,timeID_var)
    if(hasDate   ) ret = pio_inq_varid(DataSet%ncid,'date'    ,dateID    )
    if(hasDateSec) ret = pio_inq_varid(DataSet%ncid,'datesec' ,datesecID )
    if(hasUTCdate) ret = pio_inq_varid(DataSet%ncid,'utc_date',utcdateID )
    if(hasTimeVar) then
      ret = pio_inq_varndims(DataSet%ncid,timeID_var,nvdim)
      if(nvdim.eq.0) then
        hasTimeScalar = .true.
      else
        hasTimeScalar = .false.
      endif
    endif

    ! determine how to fill time values
    !-----------------------------------
    TimeOption = -1
    if(hasTimeDim) then
      if((hasTimeVar).and.(hasDate).and.(hasDateSec)) then
        if(TimePreference.eq.DATESECREF) then
          TimeOption = DATESECREF
        else
          TimeOption = COORDINATE
        endif
      elseif((hasTimeVar).and.(.not.hasDate).and.(.not.hasDateSec)) then
        if(hasUTCdate) then
          TimeOption = DATESECREF
        else
          TimeOption = COORDINATE
        endif
      elseif((.not.hasTimeVar).and.(hasDate).and.(hasDateSec)) then
        TimeOption = DATESECREF
      elseif((.not.hasTimeVar).and.(.not.hasDate).and.(.not.hasDateSec)) then
        if(hasUTCdate) then
          TimeOption = DATESECREF
        else
          TimeOption = TIMEINDEX
        endif
      endif
    else !  if(.not.hasTimeDim) then
      if(hasTimeScalar) then
        TimeOption = SCALAR
      else
        TimeOption = NULLTIME
      endif
    endif

    ! Fill time values for the given option
    !---------------------------------------
    select case(TimeOption)
      case(COORDINATE)
                   ! Read in time() coordinate values from file
                   !-----------------------------------------------
                   ret = pio_inq_dimlen(DataSet%ncid,timeID_dim, DataSet%ntime)
                   if(allocated(DataSet%time)) deallocate(DataSet%time)
                   allocate(DataSet%time(DataSet%ntime))
                   ret = pio_get_var(DataSet%ncid,timeID_var, DataSet%time)
      case(DATESECREF)
                   ! Use date()/datesec() values from file to construct a 
                   ! time() coordinate as the number of days since a specified 
                   ! reference date.
                   !-----------------------------------------------
                   ret = pio_inq_dimlen(DataSet%ncid,timeID_dim, DataSet%ntime)
                   if(allocated(DataSet%time   )) deallocate(DataSet%time   )
                   if(allocated(DataSet%date   )) deallocate(DataSet%date   )
                   if(allocated(DataSet%datesec)) deallocate(DataSet%datesec)
                   allocate(DataSet%time   (DataSet%ntime))
                   allocate(DataSet%date   (DataSet%ntime))
                   allocate(DataSet%datesec(DataSet%ntime))
  
                   if(hasUTCdate) then
                     allocate(UTC_date(DataSet%ntime))
                     ret = pio_get_var(DataSet%ncid,utcdateID,UTC_date)
                     do tt=1,DataSet%ntime
                       DataSet%date   (tt) = UTC_date(tt)/100
                       DataSet%datesec(tt) = 3600*(UTC_date(tt)-(100*DataSet%date(tt)))
                     end do 
                     deallocate(UTC_date)
                   else
                     ret = pio_get_var(DataSet%ncid,dateID   ,DataSet%date   )
                     ret = pio_get_var(DataSet%ncid,datesecID,DataSet%datesec)
                   endif
                   do tt=1,DataSet%ntime
                     year  =     DataSet%date(tt)/10000
                     month = mod(DataSet%date(tt),10000)/100
                     day   = mod(DataSet%date(tt),100  )
                     DataSet%time(tt) = float(DataSet%datesec(tt))/86400._r8 &
                                       +greg2jday(year,month,day         )   &
                                       -greg2jday(REFyear,REFmonth,REFday)
                   end do
      case(SCALAR)
                   ! Use the scalar time value in the file to set time(1)
                   !-----------------------------------------------------
                   DataSet%ntime = 1
                   if(allocated(DataSet%time)) deallocate(DataSet%time)
                   allocate(DataSet%time   (DataSet%ntime))
                   allocate(DataSet%date   (DataSet%ntime))
                   allocate(DataSet%datesec(DataSet%ntime))
                   DataSet%time   (1) = 0._r8
                   DataSet%date   (1) = 0._r8
                   DataSet%datesec(1) = 0._r8
      case(TIMEINDEX)
                   ! Create a time() array and populate it with index 
                   ! values [1,2,....,NTIME]
                   !-----------------------------------------------------
                   if(allocated(DataSet%time)) deallocate(DataSet%time)
                   allocate(DataSet%time   (DataSet%ntime))
                   allocate(DataSet%date   (DataSet%ntime))
                   allocate(DataSet%datesec(DataSet%ntime))
                   do tt=1,DataSet%ntime
                     DataSet%time   (tt) = float(tt)
                     DataSet%date   (tt) = float(tt)
                     DataSet%datesec(tt) = 0
                   end do
      case(NULLTIME)
                   ! Just set a single time value equal to 0.0
                   !-----------------------------------------------------
                   DataSet%ntime = 1
                   if(allocated(DataSet%time)) deallocate(DataSet%time)
                   allocate(DataSet%time   (DataSet%ntime))
                   allocate(DataSet%date   (DataSet%ntime))
                   allocate(DataSet%datesec(DataSet%ntime))
                   DataSet%time   (1) = 0._r8
                   DataSet%date   (1) = 0._r8
                   DataSet%datesec(1) = 0._r8
      case default
                   call endrun('ERROR: set_dataset_time() Unknown TimeOption value')
    end select

    ! End Routine
    !------------------
    return
  end subroutine set_dataset_time
  !=====================================================================


  !========================================================================
  function greg2jday( year, month, day )
    !
    ! greg2jday: Return Julian day number given Gregorian date.
    !
    ! Algorithm from Hatcher,D.A., Simple Formulae for Julian Day Numbers
    ! and Calendar Dates, Q.Jl.R.astr.Soc. (1984) v25, pp 53-55.
    !==============================================================
    implicit none
    !
    ! Passed Variables
    !-------------------
    integer,intent(in):: year, month, day
    integer           :: greg2jday
    !
    ! Local Values
    !---------------
    integer:: ap, mp
    integer:: y, d, n, g

    ! Modify year and month numbers
    !----------------------------------
    ap = year - (12 - month)/10
    mp = MOD( month-3,12 )
    if( mp < 0 ) then
      mp = mp + 12
    endif

    ! Julian day
    !----------------------------------
    y = INT( 365.25_r8 *( ap + 4712 ) )
    d = INT(  30.6_r8  *  mp + 0.5_r8 )
    n = y + d + day  + 59
    g = INT( 0.75_r8*INT( ap/100 + 49 ) ) - 38
    greg2jday = n - g

    ! End Function
    !----------------
    return
  end function greg2jday
  !========================================================================


  !========================================================================
  !========================================================================
  ! Procedures For Data Structure: OpenDataSets_t
  !   
  !========================================================================
  !========================================================================
  function get_OpenDataSet(this)
    ! Return the pointer to the DataSet structure
    !---------------------------------------------
    class(OpenDataSets_t)   :: this
    class(dataset_t),pointer:: get_OpenDataSet
    get_OpenDataSet => this%DataSet
    return
  end function get_OpenDataSet
  !========================================================================
  subroutine set_OpenDataSet(this,DataSet)
    ! Set the pointer to the DataSet structure
    !---------------------------------------------
    class(OpenDataSets_t)  :: this
    type(dataset_t),pointer:: DataSet
    this%DataSet => DataSet
    return
  end subroutine set_OpenDataSet
  !========================================================================
  function get_DataSetID(this)
    ! Return the DataSet ID
    !----------------------------------------
    class(OpenDataSets_t):: this
    integer              :: get_DataSetID
    get_DataSetID = this%DataSetID
    return
  end function get_DataSetID
  !========================================================================
  subroutine set_DataSetID(this,DataSetID)
    ! Set the DataSet ID
    !---------------------------------------
    class(OpenDataSets_t):: this
    integer              :: DataSetID
    this%DataSetID = DataSetID
    return
  end subroutine set_DataSetID
  !========================================================================
  function get_NextDataSet(this)
    ! Return the pointer to the next Dataset
    !----------------------------------------
    class(OpenDataSets_t)        :: this
    class(OpenDataSets_t),pointer:: get_NextDataSet
    get_NextDataSet => this%next
    return
  end function get_NextDataSet
  !========================================================================
  subroutine set_NextDataSet(this,next)
    ! Set the pointer to the next Dataset
    !----------------------------------------
    class(OpenDataSets_t)        :: this
    class(OpenDataSets_t),pointer:: next
    this%next => next
    return
  end subroutine set_NextDataSet
  !========================================================================
  function get_PrevDataSet(this)
    ! Return the pointer to the previous Dataset
    !--------------------------------------------
    class(OpenDataSets_t)        :: this
    class(OpenDataSets_t),pointer:: get_PrevDataSet
    get_PrevDataSet => this%prev
    return
  end function get_PrevDataSet
  !========================================================================
  subroutine set_PrevDataSet(this,prev)
    ! Set the pointer to the previous Dataset
    !--------------------------------------------
    class(OpenDataSets_t)        :: this
    class(OpenDataSets_t),pointer:: prev
    this%prev => prev
    return
  end subroutine set_PrevDataSet
  !========================================================================
  

  !========================================================================
  !========================================================================
  ! Procedures For Data Structure: OpenDataSets_t
  !
  !========================================================================
  !========================================================================
  subroutine iterate(this)
    ! Increment the Iterator to point at the next DataSet
    !----------------------------------------------------
    class(DataSetList_t):: this
    this%currDS => this%currDS%get_NextDataSet()
    return
  end subroutine iterate
  !========================================================================
  subroutine reset(this)
    ! Reset the Iterator to point at the first DataSet
    !----------------------------------------------------
    class(DataSetList_t):: this
    this%currDS => this%firstDS
    return
  end subroutine reset
  !========================================================================
  function isFirst(this)
    ! Check if the Iterator points to the end of the list.
    !----------------------------------------------------
    class(DataSetList_t)         :: this
    logical                      :: isFirst
    class(OpenDataSets_t),pointer:: lastPointer
    lastPointer => this%currDS%get_PrevDataSet()
    isFirst = .not.associated(lastPointer)
    return
  end function isFirst
  !========================================================================
  function isLast(this)
    ! Check if the Iterator points to the end of the list.
    !----------------------------------------------------
    class(DataSetList_t)         :: this
    logical                      :: isLast
    class(OpenDataSets_t),pointer:: nextPointer
    nextPointer => this%currDS%get_NextDataSet()
    isLast = .not.associated(nextPointer)
    return
  end function isLast
  !========================================================================
  function add_DataSet(this,DataSet)
    ! Add a new DataSet to the list.
    !=====================================
    !
    ! Passed Variables
    !------------------
    class(DataSetList_t)        :: this
    type(dataset_t     ),pointer:: DataSet
    integer                     :: add_DataSet
    !
    ! Local Values
    !-------------
    class(OpenDataSets_t),pointer:: new_DataSet

    ! Allocate and a new DataSet node and initialize 
    ! the pointer with the new DataSet from the input.
    !----------------------------------------------------
    allocate(new_DataSet)
    new_DataSet%DataSet => DataSet

    ! Increment the Counter to assign a unique index 
    ! to this new DataSet
    !-----------------------------------------------
    DataSet_Count    = DataSet_Count + 1
    new_DataSet%DataSetID = DataSet_Count

    ! Add it to the list
    !--------------------
    if(.not.associated(this%lastDS)) then
      ! If the current list is empty, initialize 
      ! the list with this DataSet
      !-------------------------------------------
      new_DataSet%prev => null()
      new_DataSet%next => null()
      this%firstDS     => new_DataSet
      this%lastDS      => new_DataSet
      this%currDS      => new_DataSet
    else
      ! Append this DataSet to the end of the current list.
      !----------------------------------------------------
      new_DataSet%prev => this%lastDS
      new_DataSet%next => null()
      this%lastDS%next => new_DataSet
      this%lastDS      => new_DataSet
    endif

    ! Return the DataSet ID to the calling routine
    !-----------------------------------------------
    add_DataSet = new_DataSet%DataSetID

    return
  end function add_DataSet
  !========================================================================
  subroutine remove_DataSet(this,DataSetID)
    class(DataSetList_t):: this
    integer             :: DataSetID
    !
    ! Local Values
    !-----------------
    class(OpenDataSets_t),pointer:: lastPointer
    class(OpenDataSets_t),pointer:: nextPointer
    type(dataset_t)      ,pointer:: curr_DataSet
    logical:: dataset_found

    ! cycle thru the list to find the DataSet with matching ID
    !----------------------------------------------------------
    dataset_found = .false.
    call this%reset
    do while(associated(this%currDS)) 
      if(this%currDS%get_DataSetID() .eq. DataSetID) then
        dataset_found = .true.
        exit
      endif
      call this%iterate
    end do

    if(.not.dataset_found) then
      call endrun('ERROR: remove_dataset() DataSetID NOT Found in List')
    endif

    ! Now we want to remove the current 
    ! DataSet from the list
    !------------------------------------
    if(this%isFirst().and.this%isLast()) then
      ! We are removing the ONLY DataSet in the list
      ! Deallocate currDS 
      !-----------------------------------------------
      curr_DataSet => this%currDS%get_OpenDataSet()
      call close_dataset(curr_DataSet)
      deallocate(curr_DataSet)

      ! Nullify all pointers for an empty list
      !-----------------------------------------------
      this%firstDS => null()
      this%currDS  => null()
      this%lastDS  => null()

    elseif(this%isFirst()) then
      ! We are removing the first DataSet in the list
      !-----------------------------------------------
      lastPointer => null()
      nextPointer => this%currDS%get_NextDataSet()
      call nextPointer%set_PrevDataSet(lastPointer)

      ! Deallocate currDS 
      !--------------------
      curr_DataSet => this%currDS%get_OpenDataSet()
      call close_dataset(curr_DataSet)
      deallocate(curr_DataSet)
      deallocate(this%currDS)

      ! Reset fisrtDS and currDS pointers
      !-----------------------------------
      this%firstDS => nextPointer
      this%currDS  => nextPointer

      nullify(lastPointer)
      nullify(nextPointer)
    elseif(this%isLast()) then
      ! We are removing the last DataSet in the list
      !-----------------------------------------------
      lastPointer => this%currDS%get_PrevDataSet()
      nextPointer => null()
      call lastPointer%set_NextDataSet(nextPointer)

      ! Deallocate currDS 
      !--------------------
      curr_DataSet => this%currDS%get_OpenDataSet()
      call close_dataset(curr_DataSet)
      deallocate(curr_DataSet)
      deallocate(this%currDS)

      ! Reset lastDS and currDS pointers
      !-----------------------------------
      this%lastDS => lastPointer
      this%currDS => lastPointer

      nullify(lastPointer)
      nullify(nextPointer)
    else
      ! We are removing a DataSet in the middle of the list
      !-----------------------------------------------------
      lastPointer => this%currDS%get_PrevDataSet()
      nextPointer => this%currDS%get_NextDataSet()
      call lastPointer%set_NextDataSet(nextPointer)
      call nextPointer%set_PrevDataSet(lastPointer)

      ! Deallocate currDS 
      !--------------------
      curr_DataSet => this%currDS%get_OpenDataSet()
      call close_dataset(curr_DataSet)
      deallocate(curr_DataSet)
      deallocate(this%currDS)
    
      ! Reset currDS pointer
      !-----------------------------------
      this%currDS => lastPointer

      nullify(lastPointer)
      nullify(nextPointer)
    endif

    return
  end subroutine remove_DataSet
  !========================================================================
  function get_DataSet(this,DataSetID,found)
    ! Given the DataSet ID, return a pointer to 
    ! the corresponding DataSet data structure.
    !---------------------------------------------
    class(DataSetList_t):: this
    integer             :: DataSetID
    logical,optional    :: found
    type(dataset_t),pointer:: get_DataSet
    logical:: dataset_found

    ! cycle thru the list to find the DataSet with matching ID
    !----------------------------------------------------------
    dataset_found = .false.
    call this%reset
    do while(associated(this%currDS)) 
      if(this%currDS%get_DataSetID() .eq. DataSetID) then
        dataset_found = .true.
        exit
      endif
      call this%iterate
    end do

    if(dataset_found) then
      get_DataSet => this%currDS%get_OpenDataSet()
    else
      get_DataSet => null()
    endif

    if(present(found)) found = dataset_found

    return
  end function get_DataSet
  !========================================================================
  function get_DataSet_status(this,DataSetID)
    ! Given the DataSet ID, cycle thru the List to see 
    ! if there is a corresponding DataSet data structure.
    !-----------------------------------------------------
    class(DataSetList_t):: this
    integer             :: DataSetID
    logical             :: get_DataSet_status

    ! cycle thru the list to find the DataSet with matching ID
    !----------------------------------------------------------
    get_DataSet_status = .false.
    call this%reset
    do while(associated(this%currDS)) 
      if(this%currDS%get_DataSetID() .eq. DataSetID) then
        get_DataSet_status = .true.
        exit
      endif
      call this%iterate
    end do

    return
  end function get_DataSet_status
  !========================================================================
  function get_DataSet_ncid(this,DataSetID,found)
    ! Given the DataSet ID, return a pointer to 
    ! the corresponding DataSet data structure.
    !---------------------------------------------
    class(DataSetList_t):: this
    integer             :: DataSetID
    logical,optional    :: found
    type(file_desc_t)   :: get_DataSet_ncid
    logical:: dataset_found

    ! cycle thru the list to find the DataSet with matching ID
    !----------------------------------------------------------
    dataset_found = .false.
    call this%reset
    do while(associated(this%currDS)) 
      if(this%currDS%get_DataSetID() .eq. DataSetID) then
        dataset_found = .true.
        exit
      endif
      call this%iterate
    end do

    if(dataset_found) then
      get_DataSet_ncid = this%currDS%DataSet%ncid
    endif

    if(present(found)) found = dataset_found

    return
  end function get_DataSet_ncid
  !========================================================================
!PFC DIAG Temporary function to test PIO ERRO with nf90
  function get_DataSet_ncid_nf90(this,DataSetID,found)
    ! Given the DataSet ID, return a pointer to 
    ! the corresponding DataSet data structure.
    !---------------------------------------------
    class(DataSetList_t):: this
    integer             :: DataSetID
    logical,optional    :: found
    integer             :: get_DataSet_ncid_nf90
    logical:: dataset_found

    ! cycle thru the list to find the DataSet with matching ID
    !----------------------------------------------------------
    dataset_found = .false.
    call this%reset
    do while(associated(this%currDS)) 
      if(this%currDS%get_DataSetID() .eq. DataSetID) then
        dataset_found = .true.
        exit
      endif
      call this%iterate
    end do

    if(dataset_found) then
      get_DataSet_ncid_nf90 = this%currDS%DataSet%ncid_nf90
    endif

    if(present(found)) found = dataset_found

    return
  end function get_DataSet_ncid_nf90
!PFC DIAG Temporary function to test PIO ERRO with nf90
  !========================================================================
  function get_DataSet_pio(this,DataSetID,pio_iotype,found)
    ! Given the DataSet ID, return a pointer to the corresponding 
    ! pio_subsystem for the DataSet data structure, also in the 
    ! argument list, return the pio_iotype used when the file was opened.
    !------------------------------------------------------------------
    class(DataSetList_t) :: this
    integer              :: DataSetID
    integer              :: pio_iotype
    logical,optional     :: found
    type(iosystem_desc_t),pointer:: get_DataSet_pio
    logical:: dataset_found

    ! cycle thru the list to find the DataSet with matching ID
    !----------------------------------------------------------
    dataset_found = .false.
    call this%reset
    do while(associated(this%currDS))
      if(this%currDS%get_DataSetID() .eq. DataSetID) then
        dataset_found = .true.
        exit
      endif
      call this%iterate
    end do

    if(dataset_found) then
      get_DataSet_pio => this%currDS%DataSet%pio_subsystem
      pio_iotype      =  this%currDS%DataSet%pio_iotype
    endif

    if(present(found)) found = dataset_found

    return
  end function get_DataSet_pio
  !========================================================================
  function get_DataSet_varid(this,DataSetID,varName,varID,Var,gridID)
    ! Given the DataSet ID, and a variable name, find the corresponding 
    ! varID and gridID for the variable. return .false. id the dataset
    ! ID is not valid or if the named variable is not found in the dataset.
    !---------------------------------------------------------------------
    class(DataSetList_t):: this
    integer            ,intent(in ):: DataSetID
    character(len=80)  ,intent(in ):: varName
    integer            ,intent(out):: varID
    type(dataset_var_t),intent(out):: Var
    integer            ,intent(out):: gridID
    logical                        :: get_DataSet_varid
    !
    integer:: nn

    ! cycle thru the list to find the DataSet with matching ID
    !----------------------------------------------------------
    get_DataSet_varid = .false.
    call this%reset
    do while(associated(this%currDS)) 
      if(this%currDS%get_DataSetID() .eq. DataSetID) then
        get_DataSet_varid = .true.
        exit
      endif
      call this%iterate
    end do

    if(.not.get_DataSet_varid) then
      varID  = -1
      gridID = -1
      return
    endif

    ! DatsetID is okay, now search for the named varibale
    !-----------------------------------------------------
    get_DataSet_varid = .false.
    varID  = -1
    gridID = -1
    do nn=1,this%currDS%DataSet%nvars
      if(trim(VarName).eq.trim(this%currDS%DataSet%vars(nn)%name)) then
        get_DataSet_varid = .true.
        varID          = this%currDS%DataSet%vars(nn)%varID
        gridID         = this%currDS%DataSet%vars(nn)%gridID
        Var%VarID      =      this%currDS%DataSet%vars(nn)%VarID
        Var%gridID     =      this%currDS%DataSet%vars(nn)%gridID
        Var%VarDesc    =      this%currDS%DataSet%vars(nn)%VarDesc
        Var%name       = trim(this%currDS%DataSet%vars(nn)%name)
        Var%long_name  = trim(this%currDS%DataSet%vars(nn)%long_name)
        Var%units      = trim(this%currDS%DataSet%vars(nn)%units)
        Var%fill_value =      this%currDS%DataSet%vars(nn)%fill_value
        Var%nvdim      =      this%currDS%DataSet%vars(nn)%nvdim
        if(Var%nvdim.gt.0) then
          allocate(Var%dimID(Var%nvdim))
          Var%dimID(:) = this%currDS%DataSet%vars(nn)%dimID(:)
        endif
        return
      endif
    end do
    
    return
  end function get_DataSet_varid
  !========================================================================
  function get_DataSet_timerecord(this,DataSetID,date,datesec,time_index)
    ! Given the DataSet ID, and date/datesec values, find the corresponding 
    ! time record number. return .false. if the dataset ID is not valid or 
    ! if the date/datesec values are not found in the dataset.
    !---------------------------------------------------------------------
    class(DataSetList_t):: this
    integer ,intent(in ):: DataSetID
    integer ,intent(in ):: date
    integer ,intent(in ):: datesec
    integer ,intent(out):: time_index
    logical             :: get_DataSet_timerecord
    !
    integer:: nn

    ! cycle thru the list to find the DataSet with matching ID
    !----------------------------------------------------------
    time_index = -1
    get_DataSet_timerecord = .false.
    call this%reset
    do while(associated(this%currDS)) 
      if(this%currDS%get_DataSetID() .eq. DataSetID) then
        get_DataSet_timerecord = .true.
        exit
      endif
      call this%iterate
    end do

    ! If the DataSetID was not found, or the time 
    ! record does not exist, we are done
    !----------------------------------------------
    if(.not.get_DataSet_timerecord) then
      return
    elseif(this%currDS%DataSet%ntime.lt.1) then
      return
    endif

    ! DatsetID is okay, now search for the named varibale
    !-----------------------------------------------------
    get_DataSet_timerecord = .false.
    time_index = -1
    do nn=1,this%currDS%DataSet%ntime
      if((date   .eq.this%currDS%DataSet%date   (nn)).and. &
         (datesec.eq.this%currDS%DataSet%datesec(nn))      ) then
        get_DataSet_timerecord = .true.
        time_index = nn
        return
      endif
    end do

    return
  end function get_DataSet_timerecord
  !========================================================================
  function get_DataSet_timestep(this,DataSetID,OBS_step)
    ! Given the DataSet ID, and date/datesec values, find the corresponding 
    ! time record number. return .false. if the dataset ID is not valid or 
    ! if the date/datesec values are not found in the dataset.
    !---------------------------------------------------------------------
    class(DataSetList_t):: this
    integer ,intent(in ):: DataSetID
    integer ,intent(out):: OBS_step
    logical             :: get_DataSet_timestep
    !
    ! Local Values
    !----------------
    integer                :: Year1, Month1, Day1, Sec1
    integer                :: Year2, Month2, Day2, Sec2
    type(ESMF_Time)        :: DStime1
    type(ESMF_Time)        :: DStime2
    type(ESMF_TimeInterval):: DSstep
    integer(ESMF_KIND_I8)  :: Step_ESMF
    integer                :: rc

    ! cycle thru the list to find the DataSet with matching ID
    !----------------------------------------------------------
    get_DataSet_timestep = .false.
    call this%reset
    do while(associated(this%currDS)) 
      if(this%currDS%get_DataSetID() .eq. DataSetID) then
        get_DataSet_timestep = .true.
        exit
      endif
      call this%iterate
    end do

    ! If the DataSetID was not found, or the time 
    ! record does not have 2 times, we are done
    !----------------------------------------------
    if(.not.get_DataSet_timestep) then
      return
    elseif(this%currDS%DataSet%ntime.lt.2) then
      return
    endif

    ! DatsetID is okay, now use the first 2 times to determine the 
    ! step size, the size should not require I*8 so change to I*4 
    ! for the return variable.
    !----------------------------------------------------------------
    get_DataSet_timestep = .false.
    Year1  =  this%currDS%DataSet%date(1)/10000
    Month1 = (this%currDS%DataSet%date(1) - Year1*10000)/100
    Day1   =  this%currDS%DataSet%date(1) - Year1*10000 - Month1*100
    Sec1   = this%currDS%DataSet%datesec(1)
    call ESMF_TimeSet(DStime1, yy=Year1, mm=Month1, dd=Day1, s=Sec1,   &
                              calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)

    Year2  =  this%currDS%DataSet%date(2)/10000
    Month2 = (this%currDS%DataSet%date(2) - Year2*10000)/100
    Day2   =  this%currDS%DataSet%date(2) - Year2*10000 - Month2*100
    Sec2   = this%currDS%DataSet%datesec(2)
    call ESMF_TimeSet(DStime2, yy=Year2, mm=Month2, dd=Day2, s=Sec2,   &
                              calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)

    DSstep = DStime2 - DStime1
    call ESMF_TimeIntervalGet(DSstep, s_i8=Step_ESMF, rc=rc )

    OBS_step = Step_ESMF
    get_DataSet_timestep = .true.

    return
  end function get_DataSet_timestep
  !========================================================================
  subroutine print_DataSet(this,DataSetID)
    ! Given the DataSet ID, get a pointer to 
    ! the corresponding DataSet data structure.
    !---------------------------------------------
    class(DataSetList_t):: this
    integer             :: DataSetID
    !
    ! Local Values
    !--------------
    type(dataset_t),pointer:: My_DataSet
    logical                :: dataset_found
    integer                :: nn,ii

    ! Cycle thru the list to find the DataSet with matching ID
    !----------------------------------------------------------
    dataset_found = .false.
    call this%reset
    do while(associated(this%currDS)) 
      if(this%currDS%get_DataSetID() .eq. DataSetID) then
        dataset_found = .true.
        exit
      endif
      call this%iterate
    end do

    if((dataset_found).and.masterproc) then
      ! Get the current Metadata 
      !----------------------------
      My_DataSet => this%currDS%get_OpenDataSet()

      ! Print out contents
      !-------------------
      write(iulog,*) ' '
      write(iulog,*) '============================================================== '
      write(iulog,*) 'METADATA for open_dataset',DataSetID
      write(iulog,*) '  DataSet: FilePath=',trim(My_DataSet%FilePath)
      write(iulog,*) '  DataSet: FileName=',trim(My_DataSet%FileName)
      write(iulog,*) '  DataSet:    ndims=',My_DataSet%ndims   ,'  nvars=',My_DataSet%nvars
      write(iulog,*) '  DataSet: nscalars=',My_DataSet%nscalars,' ngrids=',My_DataSet%ngrids
      write(iulog,*) '  DataSet:    ntime=',My_DataSet%ntime   ,'  natts=',My_DataSet%natts
      write(iulog,*) '  DataSet:   timeID=',My_DataSet%timeID  ,' ulimID=',My_DataSet%ulimID
      write(iulog,*) ' '
      write(iulog,*) '  DIMENSIONS:'
      write(iulog,*) '  ==========='
      do nn=1,My_DataSet%ndims
        write(iulog,*) '   ndim=',nn
        write(iulog,*) '     dimid= ',My_DataSet%dims(nn)%dimid
        write(iulog,*) '      name= ',trim(My_DataSet%dims(nn)%name)
        write(iulog,*) '       len= ',My_DataSet%dims(nn)%len
        write(iulog,*) '     ndvar= ',My_DataSet%dims(nn)%ndvar
        if(My_DataSet%dims(nn)%ndvar.gt.0) then
          do ii=1,My_DataSet%dims(nn)%ndvar
            write(iulog,*) '      ii=',ii,' varID= ',My_DataSet%dims(nn)%var(ii)%varID,  &
                                          ' name= ',trim(My_DataSet%dims(nn)%var(ii)%name)
          end do
        endif
        write(iulog,*) '    -----------------------------------'
      end do
      write(iulog,*) ' '
      write(iulog,*) '  VARIABLES:'
      write(iulog,*) '  =========='
      do nn=1,My_DataSet%nvars
        if(My_DataSet%vars(nn)%nvdim.eq.1) then
          write(iulog,*) '   nvar=',nn
          write(iulog,*) '     varid= ',My_DataSet%vars(nn)%varID
          write(iulog,*) '      name= ',trim(My_DataSet%vars(nn)%name)//'('                      &
                                      //trim(My_DataSet%dims(My_DataSet%vars(nn)%dimID(1))%name)//')'
        elseif(My_DataSet%vars(nn)%nvdim.gt.1) then
          write(iulog,*) '   nvar=',nn
          write(iulog,*) '     varid= ',My_DataSet%vars(nn)%varID
          write(iulog,*) '      name= ',trim(My_DataSet%vars(nn)%name)//                   &
                                        trim(My_DataSet%grids(My_DataSet%vars(nn)%gridID)%name)
        else
          write(iulog,*) '   nvar=',nn
          write(iulog,*) '     varid= ',My_DataSet%vars(nn)%varID
          write(iulog,*) '      name= ',trim(My_DataSet%vars(nn)%name)//'    [SCALAR]'
        endif
        if(My_DataSet%vars(nn)%gridID.gt.0) then
          write(iulog,*) '     hasTime= ',My_DataSet%grids(My_DataSet%vars(nn)%gridID)%hasTime
        endif
        write(iulog,*) '     gridID= ',My_DataSet%vars(nn)%gridID
        write(iulog,*) '    -----------------------------------'
      end do
      write(iulog,*) ' '
      write(iulog,*) '  SCALARS:'
      write(iulog,*) '  ========'
      do nn=1,My_DataSet%nscalars
        write(iulog,*) '   nscalar=',nn
        write(iulog,*) '    varid= ',My_DataSet%scalars(nn)%varID
        write(iulog,*) '     name= ',trim(My_DataSet%scalars(nn)%name)
        write(iulog,*) '      val= ',My_DataSet%scalars(nn)%val
        write(iulog,*) '    -----------------------------------'
      end do
      write(iulog,*) ' '
      write(iulog,*) '  GRIDS:'
      write(iulog,*) '  ======'
      do nn=1,My_DataSet%ngrids
        write(iulog,*) '   ngrid=',nn
        write(iulog,*) '      gridID= ',My_DataSet%grids(nn)%gridID
        write(iulog,*) '        name= ',trim(My_DataSet%grids(nn)%name)
        write(iulog,*) '     hasTime= ',My_DataSet%grids(nn)%hasTime
        write(iulog,*) '       nsdim= ',My_DataSet%grids(nn)%nsdim
        if(My_DataSet%grids(nn)%nsdim.gt.0) then
          write(iulog,*) '     dimID= ',My_DataSet%grids(nn)%dimID(:)
        endif
        write(iulog,*) '    -----------------------------------'
      end do
      write(iulog,*) ' '
      write(iulog,*) '  TIME:'
      write(iulog,*) '  ====='
      write(iulog,*) '   ntime=',My_DataSet%ntime
      if(My_DataSet%ntime.gt.0) then
        write(iulog,*) '   time   (ntime)=',My_DataSet%time(:)
        if(allocated(My_DataSet%date).and.allocated(My_DataSet%datesec)) then
          write(iulog,*) '   date   (ntime)=',My_DataSet%date(:)
          write(iulog,*) '   datesec(ntime)=',My_DataSet%datesec(:)
        endif
      endif
      write(iulog,*) '============================================================== '
      write(iulog,*) ' '

      ! Free up the pointer
      !---------------------
      nullify(My_DataSet)
    elseif(masterproc) then
      ! Print out message that the datasetID does not exist.
      !-----------------------------------------------------
      write(iulog,*) ' '
      write(iulog,*) '============================================================== '
      write(iulog,*) 'METADATA for open_dataset',DataSetID
      write(iulog,*) '  NO OPENED DATASET FOUND'
      write(iulog,*) '============================================================== '
      write(iulog,*) ' '
    endif

    return
  end subroutine print_DataSet
  !========================================================================

end module Dataset_mod
