module OBSdata_mod
!======================================================================
!
! Purpose: 
!
! 
!
!======================================================================
  ! Usefule modules
  !-----------------
  use shr_kind_mod,    only: r8=>SHR_KIND_R8, cs=>SHR_KIND_CS, cl=>SHR_KIND_CL
  use ppgrid          ,only: pver,pcols,begchunk,endchunk
  use cam_grid_support,only: cam_grid_id,cam_grid_get_dim_names
  use cam_grid_support,only: DLEN=>max_hcoordname_len
  use time_manager,    only: timemgr_time_ge, timemgr_time_inc, timemgr_datediff
  use cam_logfile,     only: iulog
  use shr_sys_mod,    only: shr_sys_flush      ! Standardized system subroutines
  use spmd_utils,      only: masterproc
  use cam_abortutils,  only: endrun
  use ncdio_atm,       only: infld, register_infld_grid
  use Dataset_mod,     only: inq_dataset_step, inq_dataset_timerecord, change_datafile
  use Dataset_mod,     only: open_dataset, close_dataset, print_dataset_metadata
  use Dataconnector_mod,only: print_dataset_paths , inq_dataPath_dims
  use Dataconnector_mod,only: apply_horiz_dataPath,apply_vert_dataPath,unpack_dataPath_result

  use netcdf
  use ESMF

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public :: DataFile_Template_t
  public :: OBS_VarConfig_t
  public :: OBS_Data_t
  private:: OBSinput_t
  private:: OBSvar_t

  public :: INV_VarConfig_t
  public :: INV_Data_t
  private:: INVinput_t
  private:: INVvar_t

  private:: ADJvar_t
  private:: ADJ_Data_t

  private:: init_OBS_Data
  private:: set_OBS_time
!xx  private:: get_OBS_data
  private:: get_OBS_data2D
  private:: get_OBS_data3D
  private:: print_OBS_Data
  private:: init_INV_Data
!xx  private:: get_INV_data
  private:: get_INV_data2D
  private:: get_INV_data3D
  private:: print_INV_Data
  private:: lookup_DataType
  private:: replace_string
  private:: days_per_month
  
!xx  interface get_OBS_data
!xx    procedure:: get_OBS_data2D
!xx    procedure:: get_OBS_data3D
!xx  end interface 

!xx  interface get_INV_data
!xx    procedure:: get_INV_data2D
!xx    procedure:: get_INV_data3D
!xx  end interface 

  type DataFile_Template_t
    integer                      :: Num_Templates
    character(len=cl),allocatable:: Path   (:)
    character(len=cs),allocatable:: Name   (:)
    character(len=16),allocatable:: Type   (:)
    character(len=32),allocatable:: TmapOpt(:)
  end type DataFile_Template_t

  type OBS_VarConfig_t
    integer                      :: Num_Vars
    character(len=8 ),allocatable:: Name        (:)
    character(len=32),allocatable:: VmapOpt     (:)
    character(len=32),allocatable:: HmapOpt     (:)
    character(len=32),allocatable:: UseOpt      (:)
    logical          ,allocatable:: HORIZ_ONLY  (:)
    integer          ,allocatable:: SpaceProfile(:)
    real(r8)         ,allocatable:: ScalingCoef (:)
    integer                      :: Num_AdjVar = 0
    character(len=32)            :: AdjOpt     = 'NONE'
    character(len=32),allocatable:: AdjVar      (:)
  end type OBS_VarConfig_t

  type OBSinput_t
    integer          :: dsID
    character(len=cl):: Path_Template
    character(len=cs):: Name_Template
    character(len=32):: TmapOpt
    integer          :: Curr_Year,Curr_Month
    integer          :: Curr_Day ,Curr_Sec
    character(len=cl):: Curr_Path
    character(len=cs):: Curr_Name
    integer          :: timeindex
    integer          :: Last_Year,Last_Month
    integer          :: Last_Day , Last_Sec
    integer          :: Next_Year,Next_Month
    integer          :: Next_Day ,Next_Sec
    integer          :: Step
    real(r8)         :: Tfrac
    integer          :: nobs
  end type OBSinput_t

  type OBSvar_t
    character(len=8 )   :: Name
    character(len=cs)   :: DSname
    character(len=32)   :: HmapOpt
    character(len=32)   :: VmapOpt
    character(len=32)   :: TmapOpt
    character(len=cl)   :: FileString
    character(len=DLEN) :: dim1name
    character(len=DLEN) :: dim2name
    integer             :: idxINPUT
    integer             :: dsID
    integer             :: pathID
    logical             :: HORIZ_ONLY
    logical             :: AdjVar = .false.
    integer             :: idxAdj = -1
    integer             :: nobs
    integer ,allocatable:: idxBUF(:)
    real(r8),allocatable:: buf(:,:,:,:) !(pcols,pver,begchunk:endchunk,nobs)
  end type OBSvar_t

  type ADJvar_t
    character(len=8 )   :: Name
    integer             :: idxVar
    logical ,allocatable:: Present(:)
    real(r8),allocatable:: buf(:,:,:)  !(ncol,nlev_d,nobs)
  end type ADJvar_t

  type ADJ_Data_t
    integer                   :: ncol
    integer                   :: nlev_d
    integer                   :: nlev_m
    integer                   :: nobs
    type(ADJvar_t),allocatable:: var(:)
  end type ADJ_Data_t

  type OBS_Data_t
      private
      integer                     :: NumInput = 0
      integer                     :: NumVars  = 0
      logical         ,allocatable:: OBSupdate   (:)
      logical         ,allocatable:: BufferReset (:)
      integer         ,allocatable:: NstepAdvance(:)
      type(OBSinput_t),allocatable:: input(:)
      type(OBSvar_t  ),allocatable:: var  (:)
      character(len=32)           :: AdjOpt
      integer                     :: NumAdjVar
      type(ADJ_Data_t)            :: ADJ
    contains
      procedure,pass:: init      =>  init_OBS_Data
      procedure,pass:: print     => print_OBS_Data
      procedure,pass:: setTime   =>   set_OBS_time
      procedure,pass:: getData3D =>   get_OBS_data3D
      procedure,pass:: getData2D =>   get_OBS_data2D
  end type OBS_Data_t

  type INV_VarConfig_t
    integer                      :: Num_Vars
    character(len=8 ),allocatable:: Name        (:)
    character(len=32),allocatable:: VmapOpt     (:)
    character(len=32),allocatable:: HmapOpt     (:)
    logical          ,allocatable:: HORIZ_ONLY  (:)
  end type INV_VarConfig_t

  type INVinput_t
    integer          :: dsID
    character(len=cl):: Path_Template
    character(len=cs):: Name_Template
    character(len=cl):: Curr_Path
    character(len=cs):: Curr_Name
    integer          :: timeindex
  end type INVinput_t

  type INVvar_t
    character(len=8 )   :: Name
    character(len=cs)   :: DSname
    character(len=32)   :: HmapOpt
    character(len=32)   :: VmapOpt
    character(len=cl)   :: FileString
    character(len=DLEN) :: dim1name
    character(len=DLEN) :: dim2name
    integer             :: idxINPUT
    integer             :: dsID
    integer             :: pathID
    logical             :: HORIZ_ONLY
    real(r8),allocatable:: val(:,:,:) !(pcols,pver,begchunk:endchunk)
  end type INVvar_t

  type INV_Data_t
      private
      integer                     :: NumInput = 0
      integer                     :: NumVars  = 0
      type(INVinput_t),allocatable:: input(:)
      type(INVvar_t  ),allocatable:: var  (:)
    contains
      procedure,pass:: init      =>  init_INV_Data
      procedure,pass:: print     => print_INV_Data
      procedure,pass:: getData2D =>   get_INV_data2D
      procedure,pass:: getData3D =>   get_INV_data3D
  end type INV_Data_t

  ! Global variables
  !---------------------
  real(r8),allocatable:: Tmp2D(:  ,:)  !(pcols     ,begchunk:endchunk))
  real(r8),allocatable:: Tmp3D(:,:,:)  !(pcols,pver,begchunk:endchunk))

contains
    !==================================================================
    subroutine init_OBS_Data(this,DataFileTemplates,OBS_VarConfig,Year,Month,Day)
      !
      !================================================================
      ! 
      ! Passed variables 
      !------------------
      class(OBS_Data_t)                   :: this
      type(DataFile_Template_t),intent(in):: DataFileTemplates
      type(OBS_VarConfig_t)    ,intent(in):: OBS_VarConfig
      integer                  ,intent(in):: Year, Month, Day
      !
      ! Local values
      !-----------------
      character(len=cl):: Obs_Path, File_Path
      character(len=cs):: Obs_Name, File_Name
      character(len=cl):: FileString
      character(len=16):: DSname
      integer          :: OBS_Step
      character(len=4 ):: Str_Year
      character(len=2 ):: Str_Month, Str_Day, Str_Hour, Str_Mdays
      character(len=5 ):: Str_Sec

      character(len=cl),allocatable:: Var_PathTemplate(:)
      character(len=cs),allocatable:: Var_NameTemplate(:)
      character(len=32),allocatable:: Var_TmapOpt     (:)
      integer          ,allocatable:: Var_OBS_Step    (:)
      integer          ,allocatable:: idxINPUT        (:)

      integer:: grid_id
      logical:: var_found, file_present, match_found
      integer:: Num_Month_Days
      integer:: ncid, ierr
      integer:: nn, tt, match_tt

      integer:: nfound, idxVar
      integer:: Step, ADJ_ncol, ADJ_nlev_m, ADJ_nlev_d
      logical:: found

      ! Allocate space for variables and initialize input values
      !-----------------------------------------------------------------
      if(OBS_VarConfig%Num_Vars.lt.1) then
        call endrun('init_OBS_Data: Num_OBS_vars error ')
      endif

      this%NumVars = OBS_VarConfig%Num_Vars
      if(allocated(this%var)) then
        do nn=1,this%NumVars
          if(allocated(this%var(nn)%buf)) deallocate(this%var(nn)%buf)
        end do
        deallocate(this%var)
      endif
      allocate(this%var(this%NumVars))

      do nn=1,this%NumVars
        this%var(nn)%Name       = trim(OBS_VarConfig%Name      (nn))
        this%var(nn)%HORIZ_ONLY =      OBS_VarConfig%HORIZ_ONLY(nn)
        this%var(nn)%HmapOpt    = trim(OBS_VarConfig%HmapOpt   (nn))
        this%var(nn)%VmapOpt    = trim(OBS_VarConfig%VmapOpt   (nn))
      end do

      ! Initialize OBS Datastructures
      !-----------------------------------
      allocate(Var_PathTemplate(this%NumVars))
      allocate(Var_NameTemplate(this%NumVars))
      allocate(Var_TmapOpt     (this%NumVars))
      allocate(Var_OBS_Step    (this%NumVars))
      do nn=1,this%NumVars

        ! Loop thru the Path/Name templates to find the file that 
        ! contains this OBS variable. For this search the file 
        ! Path/Name's are resolved using the "BEG" begining nudging time.
        !-----------------------------------------------------------------
        var_found = .false.
        do tt=1,DataFileTemplates%Num_Templates

          ! Set known DataSet Type dependent values
          !----------------------------------------
          FileString = 'NONE'
          DSname     = trim(this%var(nn)%Name)
          call lookup_DataType(trim(DataFileTemplates%Type(tt)), &
                               trim(         this%var(nn)%Name), &
                                      FileString,DSname,OBS_Step )

          ! Resolve the Path/Name for the current template
          !-----------------------------------------------
          if(trim(FileString).ne.'NONE') then
            Obs_Path = replace_string(DataFileTemplates%Path(tt),'%v',trim(FileString), &
                                                                  len_trim(FileString)  )
            Obs_Name = replace_string(DataFileTemplates%Name(tt),'%v',trim(FileString), &
                                                                  len_trim(FileString)  )
          else
            Obs_Path = trim(DataFileTemplates%Path(tt))
            Obs_Name = trim(DataFileTemplates%Name(tt))
          endif

          ! Create time strings for the given beginning time, and
          ! resolve the file path/name using the OBS templates
          !--------------------------------------------------
          Num_Month_Days = days_per_month(Year,Month)
          write(Str_Year ,'(i4.4)') Year
          write(Str_Month,'(i2.2)') Month
          write(Str_Day  ,'(i2.2)') Day
          write(Str_Mdays,'(i2.2)') Num_Month_Days
          Str_Hour = '00'
          Str_Sec  = '00000'

          File_Path = replace_string(OBS_Path ,'%y',trim(Str_Year) ,4)
          File_Path = replace_string(File_Path,'%m',trim(Str_Month),2)
          File_Path = replace_string(File_Path,'%d',trim(Str_Day)  ,2)
          File_Path = replace_string(File_Path,'%h',trim(Str_Hour) ,2)
          File_Path = replace_string(File_Path,'%s',trim(Str_Sec)  ,5)
          File_Path = replace_string(File_Path,'%n',trim(Str_Mdays),2)

          File_Name = replace_string(OBS_Name ,'%y',trim(Str_Year) ,4)
          File_Name = replace_string(File_Name,'%m',trim(Str_Month),2)
          File_Name = replace_string(File_Name,'%d',trim(Str_Day)  ,2)
          File_Name = replace_string(File_Name,'%h',trim(Str_Hour) ,2)
          File_Name = replace_string(File_Name,'%s',trim(Str_Sec)  ,5)
          File_Name = replace_string(File_Name,'%n',trim(Str_Mdays),2)

          ! Check if file is present and Inv Variable is in the file
          !-----------------------------------------------------------
          inquire(FILE=trim(File_Path)//trim(File_Name),EXIST=file_present)
          if(file_present) then
            ierr = nf90_open(trim(File_Path)//trim(File_Name),NF90_NOWRITE,ncid)
            ierr = nf90_inq_varid(ncid,trim(DSname),idxVar)
            var_found = (ierr == NF90_NOERR)
            ierr = nf90_close(ncid)
          endif

          ! If the file and variable are found, then initialize the invariant values
          !----------------------------------------------------------------------------
          if(file_present.and.var_found) then
            this%var(nn)%FileString = trim(FileString)
            this%var(nn)%DSname     = trim(DSname    )
            this%var(nn)%TmapOpt    = trim(DataFileTemplates%TmapOpt(tt))
            Var_TmapOpt     (nn)    = trim(DataFileTemplates%TmapOpt(tt))
            Var_PathTemplate(nn)    = trim(Obs_Path)
            Var_NameTemplate(nn)    = trim(Obs_Name)
            Var_OBS_Step    (nn)    =      OBS_Step
            exit
          endif
        end do ! tt=1,DataFileTemplates%Num_Templates

        ! It is an Error if the variable was not found
        !---------------------------------------------
        if(.not.var_found) then
         write(iulog,*) 'init_OBS_Data: OBS Variable Name= '//trim(this%var(nn)%Name)
         call endrun('init_OBS_Data: OBS variable not found in Datasets')
        endif
      end do ! nn=1,this%NumVars

      ! From the files found, Determine the number 
      ! of Unique dataset templates
      !----------------------------------------------------
      allocate(idxINPUT(this%NumVars))
      idxINPUT(1)   = 1
      this%NumInput = 1
      do nn=2,this%NumVars
        match_found = .false.
        do tt=(nn-1),1,-1
          if((trim(Var_PathTemplate(tt)).eq.trim(Var_PathTemplate(nn))).and. &
             (trim(Var_NameTemplate(tt)).eq.trim(Var_NameTemplate(nn)))      ) then
            match_found = .true.
            match_tt = tt
            exit
          endif
        end do
        if(match_found) then
          idxINPUT(nn)  = idxINPUT(match_tt)
        else
          this%NumInput = this%NumInput + 1
          idxINPUT(nn)  = this%NumInput
        endif
      end do ! nn

      ! Allocate and input datastructure to manage input dataset information.
      !---------------------------------------------------------------------
      allocate(this%input(this%NumInput))
      do nn=1,this%NumVars
        this%input(idxINPUT(nn))%Path_Template = trim(Var_PathTemplate(nn))
        this%input(idxINPUT(nn))%Name_Template = trim(Var_NameTemplate(nn))
        this%input(idxINPUT(nn))%TmapOpt       = trim(Var_TmapOpt     (nn))
        this%input(idxINPUT(nn))%Step          =      Var_OBS_Step    (nn)
        this%input(idxINPUT(nn))%timeindex     = 1
        this%input(idxINPUT(nn))%Curr_Year     = Year
        this%input(idxINPUT(nn))%Curr_Month    = Month
        this%input(idxINPUT(nn))%Curr_Day      = Day
        this%input(idxINPUT(nn))%Curr_Sec      = 0
      end do

      ! Resolve the Path/Name values from the templates using the 
      ! begining time for Nudging, then use them to open each unique dataset
      !-----------------------------------------------------------
      do nn=1,this%NumInput

        ! Create time strings for the begining time, and
        ! resolve the file path/name using the OBS templates
        !--------------------------------------------------
        Num_Month_Days = days_per_month(Year,Month)
        write(Str_Year ,'(i4.4)') Year
        write(Str_Month,'(i2.2)') Month
        write(Str_Day  ,'(i2.2)') Day
        write(Str_Mdays,'(i2.2)') Num_Month_Days
        Str_Hour = '00'
        Str_Sec  = '00000'

        File_Path = trim(this%input(nn)%Path_Template)
        File_Path = replace_string(File_Path,'%y',trim(Str_Year) ,4)
        File_Path = replace_string(File_Path,'%m',trim(Str_Month),2)
        File_Path = replace_string(File_Path,'%d',trim(Str_Day)  ,2)
        File_Path = replace_string(File_Path,'%h',trim(Str_Hour) ,2)
        File_Path = replace_string(File_Path,'%s',trim(Str_Sec)  ,5)
        File_Path = replace_string(File_Path,'%n',trim(Str_Mdays),2)

        File_Name = trim(this%input(nn)%Name_Template)
        File_Name = replace_string(File_Name,'%y',trim(Str_Year) ,4)
        File_Name = replace_string(File_Name,'%m',trim(Str_Month),2)
        File_Name = replace_string(File_Name,'%d',trim(Str_Day)  ,2)
        File_Name = replace_string(File_Name,'%h',trim(Str_Hour) ,2)
        File_Name = replace_string(File_Name,'%s',trim(Str_Sec)  ,5)
        File_Name = replace_string(File_Name,'%n',trim(Str_Mdays),2)

   write(iulog,*) ' PFCDIAG: open_dataset='//trim(File_Path),trim(File_Name)
        call open_dataset(this%input(nn)%dsID,trim(File_Path),trim(File_Name))
        this%input(nn)%Curr_Path = trim(File_Path)
        this%input(nn)%Curr_Name = trim(File_Name)

        ! If the OBS_Step is not known from the DataSet type, then try 
        ! to determine it from time values in the dataset
        !------------------------------------------------------------
        if(this%input(nn)%Step.lt.1) then
          if(inq_dataset_step(this%input(nn)%dsID,OBS_Step)) then
            this%input(nn)%Step = OBS_Step
          else
            call endrun('init_OBS_Data: Cannot determine OBS step size')
          endif
        endif
      end do

      do nn=1,this%NumInput
        call print_dataset_metadata(this%input(nn)%dsID)
      end do

      ! Now go back and set dsID and pathID values for var()%
      !------------------------------------------------------------
      do nn=1,this%NumVars
        this%var(nn)%idxINPUT = idxINPUT(nn)
        this%var(nn)%dsID     = this%input(idxINPUT(nn))%dsID
   write(iulog,*) ' PFCDIAG: register:',nn,' dsID=',this%var(nn)%dsID,' Var='//trim(this%var(nn)%Name)
   write(iulog,*) ' PFCDIAG:  HmapOpt='//trim(this%var(nn)%HmapOpt)
   write(iulog,*) ' PFCDIAG:  VmapOpt='//trim(this%var(nn)%VmapOpt)
        if(this%var(nn)%HORIZ_ONLY) then
          call register_infld_grid(this%var(nn)%dsID,this%var(nn)%pathID,       &
                                          'physgrid',trim(this%var(nn)%HmapOpt) )
          this%var(nn)%VmapOpt = 'HORIZ_ONLY'
        else
          call register_infld_grid(this%var(nn)%dsID,this%var(nn)%pathID,        &
                                          'physgrid',trim(this%var(nn)%HmapOpt), &
                                               'lev',trim(this%var(nn)%VmapOpt)  )
        endif
   write(iulog,*) ' PFCDIAG:  PathID=',this%var(nn)%PathID
        grid_id = cam_grid_id('physgrid')
        call cam_grid_get_dim_names(grid_id,this%var(nn)%dim1name, &
                                            this%var(nn)%dim2name  )
      end do

      do nn=1,this%NumInput
        call print_dataset_paths(this%input(nn)%dsID)
      end do

      ! Done with these at this point
      !--------------------------------
      deallocate(Var_PathTemplate)
      deallocate(Var_NameTemplate)
      deallocate(Var_TmapOpt     )
      deallocate(Var_OBS_Step    )
      deallocate(idxINPUT        )

      ! Now allocate space for OBS variable values for later use.
      ! HORIZ_ONLY values are stored in the buffer with a verical dimension of 1
      !------------------------------------------------------------------
      do nn=1,this%NumVars
        if(trim(this%var(nn)%TmapOpt).eq.'INVARIANT') then
          this%var(nn)%nobs = 1
        elseif(trim(this%var(nn)%TmapOpt).eq.'NEXT_OBS') then
          this%var(nn)%nobs = 1
        elseif(trim(this%var(nn)%TmapOpt).eq.'LINEAR') then
          this%var(nn)%nobs = 2
        else
         call endrun('init_OBS_Data_mod: Unknown TmapOpt: '//trim(this%var(nn)%TmapOpt))
        endif
        ! Copy the nobs value back to the input dataset
        !------------------------------------------------
        this%input(this%var(nn)%idxINPUT)%nobs = this%var(nn)%nobs
   write(iulog,*) ' PFCDIAG: Var:',nn,' Var='//trim(this%var(nn)%Name)
   write(iulog,*) ' PFCDIAG: TmapOpt='//trim(this%var(nn)%TmapOpt),' nobs=',this%input(this%var(nn)%idxINPUT)%nobs,this%var(nn)%nobs
   write(iulog,*) ' PFCDIAG: nobs=',this%input(this%var(nn)%idxINPUT)%nobs,this%var(nn)%nobs,this%var(nn)%idxINPUT

        if(this%var(nn)%HORIZ_ONLY) then
          allocate(this%var(nn)%buf(pcols,   1,begchunk:endchunk,this%var(nn)%nobs))
        else
          allocate(this%var(nn)%buf(pcols,pver,begchunk:endchunk,this%var(nn)%nobs))
        endif
        this%var(nn)%buf(:,:,begchunk:endchunk,:) = 0._r8

        allocate(this%var(nn)%idxBUF(this%var(nn)%nobs))
        do tt=1,this%var(nn)%nobs
          this%var(nn)%idxBUF(tt) = tt
        end do
      end do

      ! Alocate some work arrays for future buffer reads
      !--------------------------------------------------
      if(.not.allocated(Tmp2D)) allocate(Tmp2D(pcols     ,begchunk:endchunk))
      if(.not.allocated(Tmp3D)) allocate(Tmp3D(pcols,pver,begchunk:endchunk))

      allocate(this%OBSupdate   (this%NumInput))
      allocate(this%BufferReset (this%NumInput))
      allocate(this%NstepAdvance(this%NumInput))

      ! Set the NEXT/LAST times for each dataset to force 
      ! each buffer to reset for the given initial time.
      !--------------------------------------------------
      do nn=1,this%NumInput
        this%input(nn)%Last_Year = Year -1
        this%input(nn)%Last_Month= 1
        this%input(nn)%Last_Day  = 1
        this%input(nn)%Last_Sec  = 0
        this%input(nn)%Next_Year = Year -1
        this%input(nn)%Next_Month= 2
        this%input(nn)%Next_Day  = 1
        this%input(nn)%Next_Sec  = 0
      end do ! nn=1,this%NumInput

      ! Initialize flags so each OBS variable is independently mapped
      !--------------------------------------------------------------
      do tt=1,this%NumVars
        this%var(tt)%AdjVar = .false.
        this%var(tt)%idxAdj = -1
      end do

      ! Now configure values for mapping options with
      ! inter-dependent variables.
      !---------------------------------------------
      if(trim(OBS_VarConfig%AdjOpt).eq.'NONE') then
        ! There are no dependencies to worry about
        !------------------------------------------
        this%AdjOpt    = trim(OBS_VarConfig%AdjOpt)
        this%NumAdjVar = 0
      else !  (trim(OBS_VarConfig%AdjOpt).NE.'NONE')
        ! Set`config values for the AdjOpt mapping
        !--------------------------------------------
        this%AdjOpt    = trim(OBS_VarConfig%AdjOpt)
        this%NumAdjVar =      OBS_VarConfig%Num_AdjVar
        allocate(this%ADJ%var(this%NumAdjVar))
     
        ! Loop over the list of inter-dependent variables and the OBS variables 
        ! looging for the mathcing variable names. 
        ! Set the variable index for each ADJ variable.
        ! Set the OBS variable flag for each ADJ variable to indicate 
        ! that it is not independently mapped.
        !-----------------------------------------------------------------------
        nfound = 0
        do nn=1,this%NumAdjVar
        do tt=1,this%NumVars
          if(trim(this%var(tt)%Name).eq.trim(OBS_VarConfig%AdjVar(nn))) then
            nfound = nfound + 1
            this%var(tt)%idxAdj     = nn
            this%ADJ%var(nn)%idxVar = tt
            this%var(tt)%AdjVar     = .true.
          endif
        end do
        end do

   write(iulog,*) ' PFCDIAG: AdjOpts= '//trim(this%AdjOpt),' NumAdjVar=',this%NumAdjVar
   do nn=1,this%NumAdjVar
   write(iulog,*) ' PFCDIAG: ADJ%var(nn)%idxVar=',this%ADJ%var(nn)%idxVar
   end do
   do tt=1,this%NumVars
   write(iulog,*) ' PFCDIAG: var(tt)%AdjVar=',this%var(tt)%AdjVar
   end do

        ! We require that all values need for ADJ be in the list of
        ! OBS variables, evenif they are not used directly in the model. 
        !_-------------------------------------------------------------
        if(nfound.ne.this%NumAdjVar) then
          if(masterproc) then
            write(iulog,*) '==============================================================='
            write(iulog,*) 'init_OBS_Data: All of the variables used in the ADJ option list'
            write(iulog,*) '               must be included in the OBS list of variables   '
            write(iulog,*) '               even if they are not used directly for nudging, '
            write(iulog,*) '               etc... Set the Prof and Coef values to Zero when'
            write(iulog,*) '               incuding them in the OBS list to ensure that    '
            write(iulog,*) '               they are not used in the model simulation.      '
            write(iulog,*) '==============================================================='
          endif
          call endrun('init_OBS_Data: ERROR')
        endif

        ! Set values using the first ADJ variable found
        !-----------------------------------------------
        idxVar = this%ADJ%var(1)%idxVar
        this%ADJ%nobs = this%var(idxVar)%nobs
        this%ADJ%var(1)%Name = trim(this%var(idxVar)%Name)
        Step = this%input(this%var(idxVar)%idxINPUT)%Step
        allocate(this%ADJ%var(1)%Present(this%ADJ%nobs))
        this%ADJ%var(1)%Present(:) = .false.

        found = inq_dataPath_dims(trim(this%var(idxVar)%DSname),                    &
                                       this%var(idxVar)%dsID,                       &
                                       this%var(idxVar)%pathID,                     &
                                       this%ADJ%ncol,this%ADJ%nlev_m,this%ADJ%nlev_d)
        if(.not.found) then
          if(masterproc) then
            write(iulog,*) '==============================================================='
            write(iulog,*) 'init_OBS_Data: Error looking up grid dimensions for first AdjVar'
            write(iulog,*) '==============================================================='
          endif
          call endrun('init_OBS_Data: ERROR')
        endif

        ! Loop over the remaining ADV variables
        !---------------------------------------
        do nn=2,this%NumAdjVar
          idxVar = this%ADJ%var(nn)%idxVar

          ! Allocate space for flags to manage buffer indexing
          ! NOTE: the array this%ADJ%var(nn)%buf(ncol,nlev_d,nobs)
          !       is allocated and freed up as needed. 
          !----------------------------------------------------
          allocate(this%ADJ%var(nn)%Present(this%ADJ%nobs))
          this%ADJ%var(nn)%Present(:) = .false.

          ! We require that all of the variables in the 
          ! ADJ list have the same buffer size
          !-----------------------------------------------
          if(this%ADJ%nobs.ne.this%var(idxVar)%nobs) then
            if(masterproc) then
              write(iulog,*) '==============================================================='
              write(iulog,*) 'init_OBS_Data: The "nobs" size of all buffers must be the same '
              write(iulog,*) '               for all AdjVar variables   '
              write(iulog,*) '==============================================================='
            endif
            call endrun('init_OBS_Data: ERROR')
          endif

          ! We require that all of the buffers have the 
          ! same time increment
          !---------------------------------------------
          if(Step.ne.this%input(this%var(idxVar)%idxINPUT)%Step) then
            if(masterproc) then
              write(iulog,*) '==============================================================='
              write(iulog,*) 'init_OBS_Data: The Step size in time between buffer indices    '
              write(iulog,*) '               must be the same for all AdjVar variables.      '
              write(iulog,*) '==============================================================='
            endif
            call endrun('init_OBS_Data: ERROR')
          endif

          ! Get the spatial dims for the ADJ variables 
          ! so we can check if they have the same size.
          !---------------------------------------------
          found = inq_dataPath_dims(trim(this%var(idxVar)%DSname),     &
                                         this%var(idxVar)%dsID,        &
                                         this%var(idxVar)%pathID,      &
                                         ADJ_ncol,ADJ_nlev_m,ADJ_nlev_d)

          ! Monkey business to let HORIZ_ONLY 
          ! variables fall thru the 3D requirement
          !-----------------------------------------------------
          if(this%ADJ%nlev_m.eq.1) this%ADJ%nlev_m = ADJ_nlev_m
          if(this%ADJ%nlev_d.eq.1) this%ADJ%nlev_d = ADJ_nlev_d
          if(ADJ_nlev_m     .eq.1) ADJ_nlev_m      = this%ADJ%nlev_m
          if(ADJ_nlev_d     .eq.1) ADJ_nlev_d      = this%ADJ%nlev_d

          ! We require that the grid dimesnions for all 
          ! the ADJ variables have the same size.
          !---------------------------------------------
          if((.not.found).or.(this%ADJ%ncol  .ne.ADJ_ncol  ).or. &
                             (this%ADJ%nlev_m.ne.ADJ_nlev_m).or. &
                             (this%ADJ%nlev_d.ne.ADJ_nlev_d)     ) then
            if(masterproc) then
              write(iulog,*) '==============================================================='
              write(iulog,*) 'init_OBS_Data: The 3D spatial dimensions must be the same for  '
              write(iulog,*) '               all AdjVar variables'
              write(iulog,*) '==============================================================='
            endif
            call endrun('init_OBS_Data: ERROR')
          endif
        end do ! nn=2,this%NumAdjVar
      endif  ! (trim(OBS_VarConfig%AdjOpt).NE.'NONE')

      ! Initialize the buffer data for the current time
      !--------------------------------------------------
      call this%setTime(Year,Month,Day,0)

      ! End Routine
      !--------------
      return
    end subroutine init_OBS_Data
    !==================================================================


    !==================================================================
    subroutine set_OBS_time(this,Year,Month,Day,Sec) 
      !
      !================================================================
      logical,parameter:: DIAG = .TRUE.
      !
      ! Passed Variables
      !------------------
      class(OBS_Data_t) :: this
      integer,intent(in):: Year
      integer,intent(in):: Month
      integer,intent(in):: Day
      integer,intent(in):: Sec
      !
      ! Local Values
      !--------------
      type(ESMF_Time)        :: Date1,Date2
      type(ESMF_TimeInterval):: DateDiff
      character(len=cl):: File_Path
      character(len=cs):: File_Name

      integer:: YMD, YMD_LAST, YMD_NEXT, YMD_CURR
      integer:: YMDx
      logical:: After_Last, Before_Next, VARflag
      integer:: Curr_Hour, DeltaT, dlt_sec
      integer:: Num_Month_Days, time_record
      real(r8):: dlt_days
      character(len=4 ):: Str_Year
      character(len=2 ):: Str_Month, Str_Day, Str_Hour, Str_Mdays
      character(len=5 ):: Str_Sec
      integer:: nn,vv,tt,kk,idx,rc

      integer             :: idxAdj
      integer             :: idxVar
      integer             :: idxINPUT
      real(r8),allocatable:: H_Tmp  (:,:)
      real(r8),allocatable:: V_Tmp  (:,:)
      real(r8),allocatable:: DP_Data(:,:,:)
      logical             :: buf_present

      ! Loop over OBS inputs and flag any buffer updates that are needed
      !------------------------------------------------------------------
      YMD=(Year*10000) + (Month*100) + Day

      do nn=1,this%NumInput

        ! Determine if the buffer needs to be updated
        !--------------------------------------------
        YMD_LAST=(this%input(nn)%Last_Year *10000) &
                +(this%input(nn)%Last_Month*100  ) &
                + this%input(nn)%Last_Day
        YMD_NEXT=(this%input(nn)%Next_Year *10000) &
                +(this%input(nn)%Next_Month*100  ) &
                + this%input(nn)%Next_Day

        call timemgr_time_ge(YMD_LAST,this%input(nn)%Last_Sec,         &
                             YMD     ,Sec                    ,After_Last)
        call timemgr_time_ge(YMD     ,Sec                    ,           &
                             YMD_NEXT,this%input(nn)%Next_Sec,Before_Next)

        if((After_Last).and.(Before_Next)) then
          ! The buffer brackets the current time, so no update is needed 
          ! just a stright time interpolation of the current buffer values
          !--------------------------------------------------------------
          this%OBSupdate   (nn) = .false.
          this%BufferReset (nn) = .false.
          this%NstepAdvance(nn) = 0
        else
          if(.not.After_Last) then
            ! ERROR: This shoud never happen?? 
            !     time should always be moving forward... right?
            !---------------------------------------------------
            call endrun('set_OBS_time: Buffer sync error ')
          else ! if(.not.Before_Next) then
            ! The conents of the buffer need to be updated
            !-----------------------------------------------
            this%OBSupdate(nn) = .true.

            ! Determine the number of forward steps are 
            ! needed to bracket the given time
            !------------------------------------------------
            call timemgr_datediff(YMD_LAST,this%input(nn)%Last_Sec,        &
                                  YMD     ,Sec                    ,dlt_days)
            this%NstepAdvance(nn) = int(dlt_days*86400._r8)/this%input(nn)%Step

            ! If the Number of steps is large then just do a 
            ! reset of the buffer contents for the current time.
            !----------------------------------------------------
            if(this%NstepAdvance(nn).lt.this%input(nn)%nobs) then
              this%BufferReset(nn) = .false.
            else
              this%BufferReset(nn) = .true.
            endif
          endif
        endif
      end do ! nn=1,this%NumInput

      ! Prepare Work arrays to process ADJ valiables if needed.
      !--------------------------------------------------------
      if(trim(this%AdjOpt).ne.'NONE') then
        do vv=1,this%NumAdjVar
          if(this%var(this%ADJ%var(vv)%idxVar)%HORIZ_ONLY) then
            allocate(this%ADJ%var(vv)%buf(this%ADJ%ncol,1,this%ADJ%nobs))
          else
            allocate(this%ADJ%var(vv)%buf(this%ADJ%ncol,this%ADJ%nlev_d,this%ADJ%nobs))
          endif
          this%ADJ%var(vv)%Present(:) = .false.
        end do
        allocate(H_Tmp  (this%ADJ%ncol,this%ADJ%nlev_d))
        allocate(V_Tmp  (this%ADJ%ncol,this%ADJ%nlev_m))
        allocate(DP_Data(this%ADJ%ncol,this%ADJ%nlev_m,this%NumAdjVar))
      endif

      ! Now go thru and update buffer data as needed
      !-------------------------------------------------
      do nn=1,this%NumInput
        if(.not.this%OBSupdate(nn)) cycle

        if(this%BufferReset(nn)) then
          ! BUFFER RESET:  We must do a full reset of the buffer by starting 
          !                the buffers (NOBS/2) steps before the current time 
          !                and then advancing forward for NOBS steps
          !----------------------------------------------------------------------

          ! Set the current time to the new beginning of the buffer.
          !---------------------------------------------------------
          YMD_LAST=(this%input(nn)%Last_Year *10000) &
                  +(this%input(nn)%Last_Month*100  ) &
                  + this%input(nn)%Last_Day

          dlt_sec = this%input(nn)%Step*(this%NstepAdvance(nn)-(this%input(nn)%nobs/2))
          call timemgr_time_inc(YMD_LAST,this%input(nn)%Last_Sec,           &
                                YMD_CURR,this%input(nn)%Curr_Sec,dlt_sec,0,0)
          
          ! Set Last/Next Time using the Curr value 
          !--------------------------------------------------
          dlt_sec = -this%input(nn)%Step*(this%input(nn)%nobs/2)

          call timemgr_time_inc(YMD_CURR,this%input(nn)%Curr_Sec,           &
                                YMD_LAST,this%input(nn)%Last_Sec,dlt_sec,0,0)

          call timemgr_time_inc(YMD_LAST,this%input(nn)%Last_Sec,                       &
                                YMD_NEXT,this%input(nn)%Next_Sec,this%input(nn)%Step,0,0)

          YMDx = YMD_CURR
          this%input(nn)%Curr_Year =(YMDx/10000)
          YMDx = YMDx-(this%input(nn)%Curr_Year*10000)
          this%input(nn)%Curr_Month=(YMDx/100)
          this%input(nn)%Curr_Day  = YMDx-(this%input(nn)%Curr_Month*100)

          YMDx = YMD_LAST
          this%input(nn)%Last_Year =(YMDx/10000)
          YMDx = YMDx-(this%input(nn)%Last_Year*10000)
          this%input(nn)%Last_Month=(YMDx/100)
          this%input(nn)%Last_Day  = YMDx-(this%input(nn)%Last_Month*100)

          YMDx = YMD_NEXT
          this%input(nn)%Next_Year =(YMDx/10000)
          YMDx = YMDx-(this%input(nn)%Next_Year*10000)
          this%input(nn)%Next_Month=(YMDx/100)
          this%input(nn)%Next_Day  = YMDx-(this%input(nn)%Next_Month*100)

          ! Set the NstepAdvance to the size of the buffer 
          ! and the Curr_Name to force a change_datafile call
          !---------------------------------------------------
          this%NstepAdvance(nn) = this%input(nn)%nobs
          this%input(nn)%Curr_Name = 'RESET'

          if((DIAG).and.(masterproc)) then
            write(iulog,*) ' PFCDIAG: BUFFER RESET new LAST=',this%input(nn)%Last_Year,  &
                                                              this%input(nn)%Last_Month, &
                                                              this%input(nn)%Last_Day,   &
                                                              this%input(nn)%Last_Sec
            write(iulog,*) ' PFCDIAG: BUFFER RESET new NEXT=',this%input(nn)%Next_Year,  &
                                                              this%input(nn)%Next_Month, &
                                                              this%input(nn)%Next_Day,   &
                                                              this%input(nn)%Next_Sec
            write(iulog,*) ' PFCDIAG: BUFFER RESET new CURR=',this%input(nn)%Curr_Year,  &
                                                              this%input(nn)%Curr_Month, &
                                                              this%input(nn)%Curr_Day,   &
                                                              this%input(nn)%Curr_Sec
          endif
        endif ! (this%BufferReset(nn))

        ! Now Loop over the number of steps to advance the buffers
        !-----------------------------------------------------------
        do tt=1,this%NstepAdvance(nn)
          YMD_LAST=(this%input(nn)%Last_Year *10000) &
                  +(this%input(nn)%Last_Month*100  ) &
                  + this%input(nn)%Last_Day
          YMD_CURR=(this%input(nn)%Curr_Year *10000) &
                  +(this%input(nn)%Curr_Month*100  ) &
                  + this%input(nn)%Curr_Day
          YMD_NEXT=(this%input(nn)%Next_Year *10000) &
                  +(this%input(nn)%Next_Month*100  ) &
                  + this%input(nn)%Next_Day

          ! Advance buffer times by one step
          !-------------------------------------
          call timemgr_time_inc(YMD_CURR,this%input(nn)%Curr_Sec,                       &
                                YMD_CURR,this%input(nn)%Curr_Sec,this%input(nn)%Step,0,0)
          call timemgr_time_inc(YMD_LAST,this%input(nn)%Last_Sec,                       &
                                YMD_LAST,this%input(nn)%Last_Sec,this%input(nn)%Step,0,0)
          call timemgr_time_inc(YMD_NEXT,this%input(nn)%Next_Sec,                       &
                                YMD_NEXT,this%input(nn)%Next_Sec,this%input(nn)%Step,0,0)
          
          YMDx = YMD_CURR
          this%input(nn)%Curr_Year =(YMDx/10000)
          YMDx = YMDx-(this%input(nn)%Curr_Year*10000)
          this%input(nn)%Curr_Month=(YMDx/100)
          this%input(nn)%Curr_Day  = YMDx-(this%input(nn)%Curr_Month*100)

          YMDx = YMD_LAST
          this%input(nn)%Last_Year =(YMDx/10000)
          YMDx = YMDx-(this%input(nn)%Last_Year*10000)
          this%input(nn)%Last_Month=(YMDx/100)
          this%input(nn)%Last_Day  = YMDx-(this%input(nn)%Last_Month*100)

          YMDx = YMD_NEXT
          this%input(nn)%Next_Year =(YMDx/10000)
          YMDx = YMDx-(this%input(nn)%Next_Year*10000)
          this%input(nn)%Next_Month=(YMDx/100)
          this%input(nn)%Next_Day  = YMDx-(this%input(nn)%Next_Month*100)

          if((DIAG).and.(masterproc)) then
            write(iulog,*) ' PFCDIAG: ADVANCE 1 STEP: LAST=',this%input(nn)%Last_Year,  &
                                                             this%input(nn)%Last_Month, &
                                                             this%input(nn)%Last_Day,   &
                                                             this%input(nn)%Last_Sec
            write(iulog,*) ' PFCDIAG: ADVANCE 1 STEP: NEXT=',this%input(nn)%Next_Year,  &
                                                             this%input(nn)%Next_Month, &
                                                             this%input(nn)%Next_Day,   &
                                                             this%input(nn)%Next_Sec
            write(iulog,*) ' PFCDIAG: ADVANCE 1 STEP: CURR=',this%input(nn)%Curr_Year,  &
                                                             this%input(nn)%Curr_Month, &
                                                             this%input(nn)%Curr_Day,   &
                                                             this%input(nn)%Curr_Sec
          endif

          ! Resolve File NAME for the CURR time
          !--------------------------------------
          Num_Month_Days = days_per_month(this%input(nn)%Curr_Year, &
                                          this%input(nn)%Curr_Month )
          Curr_Hour = this%input(nn)%Curr_Sec/3600
          write(Str_Mdays,'(i2.2)') Num_Month_Days
          write(Str_Year ,'(i4.4)') this%input(nn)%Curr_Year
          write(Str_Month,'(i2.2)') this%input(nn)%Curr_Month
          write(Str_Day  ,'(i2.2)') this%input(nn)%Curr_Day
          write(Str_Hour ,'(i2.2)') Curr_Hour
          write(Str_Sec  ,'(i5.5)') this%input(nn)%Curr_Sec

          File_name = trim(this%input(nn)%Name_Template)
          File_Name = replace_string(File_Name,'%y',trim(Str_Year) ,4)
          File_Name = replace_string(File_Name,'%m',trim(Str_Month),2)
          File_Name = replace_string(File_Name,'%d',trim(Str_Day)  ,2)
          File_Name = replace_string(File_Name,'%h',trim(Str_Hour) ,2)
          File_Name = replace_string(File_Name,'%s',trim(Str_Sec)  ,5)
          File_Name = replace_string(File_Name,'%n',trim(Str_Mdays),2)

          ! If the file NAME has changed, then change the dataset file
          !----------------------------------------------------------------
          if(trim(File_Name).ne.trim(this%input(nn)%Curr_Name)) then
            call change_datafile(this%input(nn)%dsID,trim(File_Name))
            this%input(nn)%Curr_Name = trim(File_Name)
          endif

          ! Determine the time record for the current time
          !-------------------------------------------------
          if(inq_dataset_timerecord(this%input(nn)%dsID,YMD_CURR,            &
                                         this%input(nn)%Curr_Sec,time_record)) then
            this%input(nn)%timeindex = time_record
          else
            call endrun('set_OBS_time: Error finding time_record ')
          endif

          ! Loop over OBS variables, and update buffer values
          ! if the variable uses this input.
          !--------------------------------------------------
          do vv=1,this%NumVars
            if(this%input(nn)%dsID.ne.this%var(vv)%dsID) cycle

            ! Rotate buffer Indicies for this variable
            !--------------------------------------------
            idx = this%var(vv)%idxBUF(1)
            do kk=1,(this%var(vv)%nobs-1)
              this%var(vv)%idxBUF(kk) = this%var(vv)%idxBUF(kk+1)
            end do
            this%var(vv)%idxBUF(this%var(vv)%nobs) = idx

            ! Process ADJ data differently, otherwise read 
            ! other values into the end of the buffer
            !-------------------------------------------------
            if(this%var(vv)%AdjVar) then
              ! INTER-DEPENDENT(ADJ) VAR: Horizontally map this Var and 
              ! store results for later inter-dependent vertical mapping
              !---------------------------------------------------------
              idxAdj = this%var(vv)%idxAdj

              ! Get Horizontally mapped data from dsID
              ! Store the results for later use
              !----------------------------------------
              if(this%var(vv)%HORIZ_ONLY) then
                call apply_horiz_dataPath(trim(this%var(vv)%DSname),          &
                                               this%var(vv)%dsID,             &
                                               this%var(vv)%pathID,           &
                                               this%input(nn)%timeindex,      &
                                               1,pcols,begchunk,endchunk,H_Tmp)
                this%ADJ%var(idxAdj)%buf(:,1,idx) = H_Tmp(:,1)
              else
                call apply_horiz_dataPath(trim(this%var(vv)%DSname),                 &
                                               this%var(vv)%dsID,                    &
                                               this%var(vv)%pathID,                  &
                                               this%input(nn)%timeindex,             &
                                               1,pcols,1,pver,begchunk,endchunk,H_Tmp)
                this%ADJ%var(idxAdj)%buf(:,:,idx) = H_Tmp(:,:)
              endif

              ! Flag the buffer index for this var as having valid data
              !--------------------------------------------------------
              this%ADJ%var(idxAdj)%Present(idx) = .true.
   
            else ! .not.(this%var(vv)%AdjVar)
              ! INDEPENDENT VAR: Read H/V mapped data into the end 
              ! of the buffer for the given DataSetID and pathID
              !----------------------------------------------------
              if(this%var(vv)%HORIZ_ONLY) then
                Tmp2D(:,begchunk:endchunk) = 0._r8
                call infld(this%var(vv)%DSname,this%var(vv)%dsID,this%var(vv)%pathID, &
                           this%var(vv)%dim1name,this%var(vv)%dim2name,               &
                           1,pcols,begchunk,endchunk,Tmp2D,VARflag,                   &
                           timelevel=this%input(nn)%timeindex                         )
                if(VARflag) then
                  this%var(vv)%buf(:,1,begchunk:endchunk,idx) = Tmp2D(:,begchunk:endchunk)
                else
                  call endrun('set_OBS_time: Variable '//               & 
                               trim(this%var(vv)%DSname)//' is missing ')
                endif
              else
                Tmp3D(:,:,begchunk:endchunk) = 0._r8
                call infld(this%var(vv)%DSname,this%var(vv)%dsID,this%var(vv)%pathID, &
                           this%var(vv)%dim1name,'lev',this%var(vv)%dim2name,         &
                           1,pcols,1,pver,begchunk,endchunk,Tmp3D,VARflag,            &
                           timelevel=this%input(nn)%timeindex                         )
                if(VARflag) then
                  this%var(vv)%buf(:,:,begchunk:endchunk,idx) = Tmp3D(:,:,begchunk:endchunk)
                else
                  call endrun('set_OBS_time: Variable '//               & 
                               trim(this%var(vv)%DSname)//' is missing ')
                endif
              endif
            endif
          end do ! vv=1,this%NumVars
        end do ! tt=1,this%NstepAdvance(nn)

        ! Verify that the buffer brackets the current time
        !--------------------------------------------------
        call timemgr_time_ge(YMD_LAST,this%input(nn)%Last_Sec,         &
                             YMD     ,Sec                    ,After_Last)
        call timemgr_time_ge(YMD     ,Sec                    ,           &
                             YMD_NEXT,this%input(nn)%Next_Sec,Before_Next)
        if((.not.After_Last).or.(.not.Before_Next)) then
          ! if the current time is not bracketed by 
          ! Last/Next, then there is a problem
          !------------------------------------------
          call endrun('set_OBS_time: Buffer Bracketing error ')
        endif
      end do ! nn=1,this%NumInput

      ! Now Process Interdependent ADJ variables,
      ! unpack the data into the ourput PE format, and 
      ! store the results into the OBS variable buffers.
      !--------------------------------------------------
      if(trim(this%AdjOpt).ne.'NONE') then

        ! Loop over buffer indices
        !---------------------------
        do idx = 1,this%ADJ%nobs

          ! We require data to be present in all of the ADJ variable buffers
          !-----------------------------------------------------------------
          buf_present = .true.
          do vv=1,this%NumAdjVar
            if(.not.this%ADJ%var(vv)%Present(idx)) then
              buf_present = .false.
            endif
          end do
          if(.not.buf_present) cycle
       
          ! Apply AdjOpt to data for the current buffer index
          !------------------------------------------------------

          !=======================================================================
          !=======================================================================
          ! DIAG: For Testing, just apply the (INDEX-SHIFT) Vmap to each
          !       ADJ var buffer data, this should recover results of the AdjOpt='NONE'
          !       option if all of the logic os correct.
          !=======================================================================
          do vv=1,this%NumAdjVar
            idxVar   = this%ADJ%var(vv)%idxVar
            idxINPUT = this%var(idxVar)%idxINPUT
            if(this%var(idxVar)%HORIZ_ONLY) then
              ! Just a copy for Horizontal Only data
              !---------------------------------------
              DP_Data(:,1,vv) = this%ADJ%var(vv)%buf(:,1,idx)
            else
              ! Apply Mapping as for AdjOpt = 'NONE'
              !---------------------------------------
              H_Tmp(:,:) = this%ADJ%var(vv)%buf(:,:,idx)
              call apply_vert_dataPath(trim(this%var(idxVar)%DSname),                 &
                                            this%var(idxVar)%dsID,                    &
                                            this%var(idxVar)%pathID,                  &
                                            this%input(idxINPUT)%timeindex,H_Tmp,V_Tmp)
              DP_Data(:,:,vv) = V_Tmp(:,:)
            endif
          end do ! vv=1,this%NumAdjVar
          !=======================================================================
          !=======================================================================


          ! Loop over ADJ variables, unpack values and then store 
          ! the results in the OBS buffer for the varaible.
          !----------------------------------------------------
          do vv=1,this%NumAdjVar
            idxVar     = this%ADJ%var(vv)%idxVar
            idxINPUT   = this%var(idxVar)%idxINPUT
            V_Tmp(:,:) = DP_Data(:,:,vv)
            if(this%var(idxVar)%HORIZ_ONLY) then
              call unpack_dataPath_result(trim(this%var(idxVar)%DSname),            &
                                               this%var(idxVar)%dsID,               &
                                               this%var(idxVar)%pathID,             &
                                               this%input(idxINPUT)%timeindex,      &
                                               1,pcols,begchunk,endchunk,V_Tmp,Tmp2D)
              this%var(idxVar)%buf(:,1,begchunk:endchunk,idx) = Tmp2D(:,begchunk:endchunk)
            else
              call unpack_dataPath_result(trim(this%var(idxVar)%DSname),                   &
                                               this%var(idxVar)%dsID,                      &
                                               this%var(idxVar)%pathID,                    &
                                               this%input(idxINPUT)%timeindex,             &
                                               1,pcols,1,pver,begchunk,endchunk,V_Tmp,Tmp3D)
              this%var(idxVar)%buf(:,:,begchunk:endchunk,idx) = Tmp3D(:,:,begchunk:endchunk)
            endif
          end do ! vv=1,this%NumAdjVar
        end do ! idx = 1,this%ADJ%nobs

        ! Free up ADJ Work arrays 
        !-------------------------
        do vv=1,this%NumAdjVar
          deallocate(this%ADJ%var(vv)%buf)
        end do
        deallocate(H_Tmp)
        deallocate(V_Tmp)
        deallocate(DP_Data)
      endif

      ! Calculate the fractional difference between the 
      ! current time and the bracketing LAST/NEXT values
      !-----------------------------------------------------
      do nn=1,this%NumInput
        call ESMF_TimeSet(Date1,YY=this%input(nn)%Last_Year,  &
                                MM=this%input(nn)%Last_Month, &
                                DD=this%input(nn)%Last_Day,   &
                                 S=this%input(nn)%Last_Sec    )
        call ESMF_TimeSet(Date2,YY=Year,  &
                                MM=Month, &
                                DD=Day,   &
                                 S=Sec    )
        DateDiff =Date2-Date1
        call ESMF_TimeIntervalGet(DateDiff,S=DeltaT,rc=rc)
        this%input(nn)%Tfrac= float(DeltaT)/float(this%input(nn)%Step)
      end do ! nn=1,this%NumInput

      ! End Routine
      !--------------
      return
    end subroutine set_OBS_time
    !==================================================================


    !==================================================================
    subroutine get_OBS_data2D(this,VarName,Output)
      !
      !================================================================
      !
      ! Passed Variables
      !------------------
      class(OBS_Data_t)           :: this
      character(len=*),intent(in ):: VarName
      real(r8)        ,intent(out):: Output(pcols,1,begchunk:endchunk)
      !
      ! Local Values
      !--------------
      integer :: vv,nn
      integer :: idxLAST,idxNEXT
      real(r8):: wgtLAST,wgtNEXT

      ! Loop thue the OBS variables and find the index for the given VarName
      !-----------------------------------------------------------------------
      do vv=1,this%NumVars
        if(trim(VarName).eq.trim(this%var(vv)%Name)) then

          ! Verify the the variable is HORIZ_ONLY
          !---------------------------------------
          if(.not.this%var(vv)%HORIZ_ONLY) then
            call endrun('init_OBSdata_mod: get_OBS_data2D HORIZ_ONLY error ')
          endif

          ! Use the specified TmapOpt option
          !----------------------------------
          if(trim(this%var(vv)%TmapOpt).eq.'NEXT_OBS') then
            ! The NEXT OBS value 
            !---------------------
            idxNEXT = this%var(vv)%idxBUF((this%var(vv)%nobs/2)+1)
            Output(:,1,begchunk:endchunk) = this%var(vv)%buf(:,1,begchunk:endchunk,idxNEXT)
          elseif(trim(this%var(vv)%TmapOpt).eq.'LINEAR') then
            ! Linear interp between LAST/NEXT OBS values
            !-------------------------------------------
            idxLAST = this%var(vv)%idxBUF((this%var(vv)%nobs/2)  )
            idxNEXT = this%var(vv)%idxBUF((this%var(vv)%nobs/2)+1)
            wgtLAST = (1._r8 - this%input(this%var(vv)%idxINPUT)%Tfrac)
            wgtNEXT =          this%input(this%var(vv)%idxINPUT)%Tfrac
            Output(:,1,begchunk:endchunk) =                                              &
                                 wgtLAST*this%var(vv)%buf(:,1,begchunk:endchunk,idxLAST) &
                                +wgtNEXT*this%var(vv)%buf(:,1,begchunk:endchunk,idxNEXT)
          else
            call endrun('init_OBSdata_mod: get_OBS_data2D TmapOpt error ')
          endif
        endif
      end do ! vv=1,this%NumVars

      ! End Routine
      !--------------
      return
    end subroutine get_OBS_data2D
    !==================================================================


    !==================================================================
    subroutine get_OBS_data3D(this,VarName,Output)
      !
      !================================================================
      !
      ! Passed Variables
      !------------------
      class(OBS_Data_t)           :: this
      character(len=*),intent(in ):: VarName
      real(r8)        ,intent(out):: Output(pcols,pver,begchunk:endchunk)
      !
      ! Local Values
      !--------------
      integer :: vv,nn
      integer :: idxLAST,idxNEXT
      real(r8):: wgtLAST,wgtNEXT

      ! Loop thue the OBS variables and find the index for the given VarName
      !-----------------------------------------------------------------------
      do vv=1,this%NumVars
        if(trim(VarName).eq.trim(this%var(vv)%Name)) then

          ! Verify the the variable is HORIZ_ONLY
          !---------------------------------------
          if(this%var(vv)%HORIZ_ONLY) then
            call endrun('init_OBSdata_mod: get_OBS_data3D HORIZ_ONLY error ')
          endif

          ! Use the specified TmapOpt option
          !----------------------------------
          if(trim(this%var(vv)%TmapOpt).eq.'NEXT_OBS') then
            ! The NEXT OBS value 
            !---------------------
            idxNEXT = this%var(vv)%idxBUF((this%var(vv)%nobs/2)+1)
            Output(:,:,begchunk:endchunk) = this%var(vv)%buf(:,:,begchunk:endchunk,idxNEXT)
          elseif(trim(this%var(vv)%TmapOpt).eq.'LINEAR') then
            ! Linear interp between LAST/NEXT OBS values
            !-------------------------------------------
            idxLAST = this%var(vv)%idxBUF((this%var(vv)%nobs/2)  )
            idxNEXT = this%var(vv)%idxBUF((this%var(vv)%nobs/2)+1)
            wgtLAST = (1._r8 - this%input(this%var(vv)%idxINPUT)%Tfrac)
            wgtNEXT =          this%input(this%var(vv)%idxINPUT)%Tfrac
            Output(:,:,begchunk:endchunk) =                                                &
                                   wgtLAST*this%var(vv)%buf(:,:,begchunk:endchunk,idxLAST) &
                                  +wgtNEXT*this%var(vv)%buf(:,:,begchunk:endchunk,idxNEXT)
          else
            call endrun('init_OBSdata_mod: get_OBS_data3D error ')
          endif
        endif
      end do ! vv=1,this%NumVars

      ! End Routine
      !--------------
      return
    end subroutine get_OBS_data3D
    !==================================================================


    !==================================================================
    subroutine print_OBS_Data(this)
      !
      !================================================================
      !
      ! Passed Variables
      !-----------------
      class(OBS_Data_t):: this
      !
      ! Local Values
      !-------------
      integer:: nn

      ! Output the list of OBS datasets, and the number of OBS variables
      !--------------------------------------------------------------------
      if(masterproc) then
        write(iulog,*) 'OBS_Data: '
        write(iulog,*) 'OBS_Data: Input Datasets for OBS variables:'//              &
                                                      ' Num_OBS_input=',this%NumInput
        write(iulog,*) 'OBS_Data:--------------------------------------------------- '
        do nn=1,this%NumInput
          write(iulog,*) 'OBS_Data: Input dataset nn=',nn
          write(iulog,*) 'OBS_Data:    Path_Template=',trim(this%input(nn)%Path_Template)
          write(iulog,*) 'OBS_Data:    Name_Template=',trim(this%input(nn)%Name_Template)
          write(iulog,*) 'OBS_Data:    dsID=',this%input(nn)%dsID, &
                            ' TmapOpt='//trim(this%input(nn)%TmapOpt)
          write(iulog,*) 'OBS_Data:    Step=',this%input(nn)%Step, &
                                     ' nobs=',this%input(nn)%nobs
        end do
        write(iulog,*) 'OBS_Data: '
        write(iulog,*) 'OBS_Data: OBS variables: Num_OBS_Vars=',this%NumVars
        write(iulog,*) 'OBS_Data:--------------------------------------------------- '
        do nn=1,this%NumVars
          write(iulog,*) 'OBS_Data: OBS Variable nn=',nn
          write(iulog,*) 'OBS_Data:     Name='//trim(this%var(nn)%Name  )// &
                                        ' :  '//trim(this%var(nn)%DSname)
          write(iulog,*) 'OBS_Data: idxINPUT =',     this%var(nn)%idxINPUT
          write(iulog,*) 'OBS_Data:    dsID  =',     this%var(nn)%dsID
          write(iulog,*) 'OBS_Data:    pathID=',     this%var(nn)%pathID
          write(iulog,*) 'OBS_Data:   HmapOpt=',trim(this%var(nn)%HmapOpt)
          write(iulog,*) 'OBS_Data:   VmapOpt=',trim(this%var(nn)%VmapOpt)
          write(iulog,*) 'OBS_Data:    AdjVar=',     this%var(nn)%AdjVar
        end do
        write(iulog,*) 'OBS_Data: '
        write(iulog,*) 'OBS_Data: OBS variables: ADJ Option=',trim(this%AdjOpt)
        write(iulog,*) 'OBS_Data: OBS variables: Num_ADJ_Vars=',this%NumAdjVar
        write(iulog,*) 'OBS_Data:--------------------------------------------------- '
        write(iulog,*) 'OBS_Data: '
        endif

        ! End Routine
        !--------------
        return
        end subroutine print_OBS_Data
    !==================================================================


    !==================================================================
    subroutine init_INV_Data(this,DataFileTemplates,INV_VarConfig)
      !
      !================================================================
      ! 
      ! Passed variables 
      !------------------
      class(INV_Data_t)                    :: this
      type (DataFile_Template_t),intent(in):: DataFileTemplates
      type (INV_VarConfig_t)    ,intent(in):: INV_VarConfig
      !
      ! Local values
      !-----------------
      character(len=cl):: Inv_Path, File_Path
      character(len=cs):: Inv_Name, File_Name
      character(len=cl):: FileString
      character(len=16):: DSname
      integer          :: OBS_Step

      character(len=cl),allocatable:: Var_PathTemplate(:)
      character(len=cs),allocatable:: Var_NameTemplate(:)
      integer          ,allocatable:: idxINPUT        (:)

      integer:: grid_id
      logical:: var_found, file_present, match_found
      integer:: ncid, ierr, idxvar
      integer:: nn, tt, match_tt
      logical:: VARflag

      ! Allocate space for variables and initialize input values
      !-----------------------------------------------------------------
      if(INV_VarConfig%Num_Vars.lt.0) then
        call endrun('init_INV_Data: Num_OBS_vars error ')
      elseif(INV_VarConfig%Num_Vars.eq.0) then
        ! There are no variables to initialize, so just return and
        ! hope tha the user knows what they are doing
        !---------------------------------------------------------
        return
      endif

      this%NumVars = INV_VarConfig%Num_Vars
      if(allocated(this%var)) then
        do nn=1,this%NumVars
          if(allocated(this%var(nn)%val)) deallocate(this%var(nn)%val)
        end do
        deallocate(this%var)
      endif
      allocate(this%var(this%NumVars))

      do nn=1,this%NumVars
        this%var(nn)%Name       = trim(INV_VarConfig%Name      (nn))
        this%var(nn)%HORIZ_ONLY =      INV_VarConfig%HORIZ_ONLY(nn)
        this%var(nn)%HmapOpt    = trim(INV_VarConfig%HmapOpt   (nn))
        this%var(nn)%VmapOpt    = trim(INV_VarConfig%VmapOpt   (nn))
      end do

      ! Initialize OBS Datastructures
      !-----------------------------------
      allocate(Var_PathTemplate(this%NumVars))
      allocate(Var_NameTemplate(this%NumVars))
      do nn=1,this%NumVars

        ! Loop thru the Path/Name templates to find the file that 
        ! contains this INV variable. 
        !-----------------------------------------------------------
        var_found = .false.
        do tt=1,DataFileTemplates%Num_Templates

          ! Set known DataSet Type dependent values
          !----------------------------------------
          FileString = 'NONE'
          DSname     = trim(this%var(nn)%Name)
          call lookup_DataType(trim(DataFileTemplates%Type(tt)), &
                               trim(         this%var(nn)%Name), &
                                      FileString,DSname,OBS_Step )

          ! Resolve the Path/Name for the current template
          !-----------------------------------------------
          if(trim(FileString).ne.'NONE') then
            Inv_Path = replace_string(DataFileTemplates%Path(tt),'%v',trim(FileString), &
                                                                  len_trim(FileString)  )
            Inv_Name = replace_string(DataFileTemplates%Name(tt),'%v',trim(FileString), &
                                                                  len_trim(FileString)  )
          else
            Inv_Path = trim(DataFileTemplates%Path(tt))
            Inv_Name = trim(DataFileTemplates%Name(tt))
          endif

          File_Path = trim(Inv_Path)
          File_Name = trim(Inv_Name)

          ! Check if file is present and Inv Variable is in the file
          !-----------------------------------------------------------
          inquire(FILE=trim(File_Path)//trim(File_Name),EXIST=file_present)
          if(file_present) then
            ierr = nf90_open(trim(File_Path)//trim(File_Name),NF90_NOWRITE,ncid)
            ierr = nf90_inq_varid(ncid,trim(DSname),idxvar)
            var_found = (ierr == NF90_NOERR)
            ierr = nf90_close(ncid)
          endif

          ! If the file and variable are found, then initialize the invariant values
          !----------------------------------------------------------------------------
          if(file_present.and.var_found) then
            this%var(nn)%FileString = trim(FileString)
            this%var(nn)%DSname     = trim(DSname    )
            Var_PathTemplate(nn)    = trim(Inv_Path)
            Var_NameTemplate(nn)    = trim(Inv_Name)
            exit
          endif
        end do ! tt=1,DataFileTemplates%Num_Templates

        ! It is an Error if the variable was not found
        !---------------------------------------------
        if(.not.var_found) then
         write(iulog,*) 'init_INV_Data: INV Variable Name= '//trim(this%var(nn)%Name)
         call endrun('init_INV_Data: INV variable not found in Datasets')
        endif
      end do ! nn=1,this%NumVars

      ! From the files found, Determine the number 
      ! of Unique dataset templates
      !----------------------------------------------------
      allocate(idxINPUT(this%NumVars))
      idxINPUT(1)   = 1
      this%NumInput = 1
      do nn=2,this%NumVars
        match_found = .false.
        do tt=(nn-1),1,-1
          if((trim(Var_PathTemplate(tt)).eq.trim(Var_PathTemplate(nn))).and. &
             (trim(Var_NameTemplate(tt)).eq.trim(Var_NameTemplate(nn)))      ) then
            match_found = .true.
            match_tt = tt
            exit
          endif
        end do
        if(match_found) then
          idxINPUT(nn)  = idxINPUT(match_tt)
        else
          this%NumInput = this%NumInput + 1
          idxINPUT(nn)  = this%NumInput
        endif
      end do ! nn

      ! Allocate and init datastructure to manage input dataset information.
      !---------------------------------------------------------------------
      allocate(this%input(this%NumInput))
      do nn=1,this%NumVars
        this%input(idxINPUT(nn))%Path_Template = trim(Var_PathTemplate(nn))
        this%input(idxINPUT(nn))%Name_Template = trim(Var_NameTemplate(nn))
        this%input(idxINPUT(nn))%timeindex     = 1
      end do

      ! Resolve the Path/Name values from the templates using the 
      ! begining time for Nudging, then use them to open each unique dataset
      !-----------------------------------------------------------
      do nn=1,this%NumInput

        ! Create time strings for the begining time, and
        ! resolve the file path/name using the OBS templates
        !--------------------------------------------------
        File_Path = trim(this%input(nn)%Path_Template)
        File_name = trim(this%input(nn)%Name_Template)

        call open_dataset(this%input(nn)%dsID,trim(File_Path),trim(File_Name))
        this%input(nn)%Curr_Path = trim(File_Path)
        this%input(nn)%Curr_Name = trim(File_Name)
      end do

      ! Now go back and set dsID and pathID values for var()%
      !------------------------------------------------------------
      do nn=1,this%NumVars
        this%var(nn)%idxINPUT = idxINPUT(nn)
        this%var(nn)%dsID     = this%input(idxINPUT(nn))%dsID
        call register_infld_grid(this%var(nn)%dsID,this%var(nn)%pathID,        &
                                        'physgrid',trim(this%var(nn)%HmapOpt), &
                                             'lev',trim(this%var(nn)%VmapOpt)  )
        grid_id = cam_grid_id('physgrid')
        call cam_grid_get_dim_names(grid_id,this%var(nn)%dim1name, &
                                            this%var(nn)%dim2name  )
      end do

      do nn=1,this%NumInput
        call print_dataset_metadata(this%input(nn)%dsID)
        call print_dataset_paths   (this%input(nn)%dsID)
      end do

      ! Done with these at this point
      !--------------------------------
      deallocate(Var_PathTemplate)
      deallocate(Var_NameTemplate)
      deallocate(idxINPUT        )

      ! Now allocate space for INV variable values and read them.
      ! HORIZ_ONLY values are stored in the buffer with a verical dimension of 1
      !------------------------------------------------------------------
      if(.not.allocated(Tmp2D)) allocate(Tmp2D(pcols     ,begchunk:endchunk))
      if(.not.allocated(Tmp3D)) allocate(Tmp3D(pcols,pver,begchunk:endchunk))

      do nn=1,this%NumVars
        if(this%var(nn)%HORIZ_ONLY) then
          allocate(this%var(nn)%val(pcols,   1,begchunk:endchunk))
          Tmp2D(:,begchunk:endchunk) = 0._r8
          call infld(this%var(nn)%DSname,this%var(nn)%dsID,this%var(nn)%pathID, &
                     this%var(nn)%dim1name,this%var(nn)%dim2name,               &
                     1,pcols,begchunk,endchunk,Tmp2D,VARflag,                   &
                     timelevel=1                                                )
          if(VARflag) then
            this%var(nn)%val(:,1,begchunk:endchunk) = Tmp2D(:,begchunk:endchunk)
          else
            call endrun('init_INV_Data: Variable '//trim(this%var(nn)%DSname)//' is missing')
          endif
        else
          allocate(this%var(nn)%val(pcols,pver,begchunk:endchunk))
          Tmp2D(:,begchunk:endchunk) = 0._r8
          call infld(this%var(nn)%DSname,this%var(nn)%dsID,this%var(nn)%pathID, &
                     this%var(nn)%dim1name,'lev',this%var(nn)%dim2name,         &
                     1,pcols,1,pver,begchunk,endchunk,Tmp3D,VARflag,            &
                     timelevel=1                                                )
          if(VARflag) then
            this%var(nn)%val(:,:,begchunk:endchunk) = Tmp3D(:,:,begchunk:endchunk)
          else
            call endrun('init_INV_Data: Variable '//trim(this%var(nn)%DSname)//' is missing')
          endif
        endif
      end do

      ! End Routine
      !--------------
      return
    end subroutine init_INV_Data
    !==================================================================


    !==================================================================
    subroutine get_INV_data2D(this,VarName,Output)
      !
      !================================================================
      !
      ! Passed Variables
      !------------------
      class(INV_Data_t)           :: this
      character(len=*),intent(in ):: VarName
      real(r8)        ,intent(out):: Output(1:pcols,begchunk:endchunk)
      !
      ! Local Values
      !--------------
      integer :: vv

      ! Loop thue the OBS variables and find the index for the given VarName
      !-----------------------------------------------------------------------
      do vv=1,this%NumVars
        if(trim(VarName).eq.trim(this%var(vv)%Name)) then

          ! Verify the the variable is HORIZ_ONLY
          !---------------------------------------
          if(.not.this%var(vv)%HORIZ_ONLY) then
            call endrun('get_INV_data2D: HORIZ_ONLY error ')
          endif

          ! copy the values to the output array
          !-------------------------------------
          Output(:,begchunk:endchunk) = this%var(vv)%val(:,1,begchunk:endchunk)
        endif
      end do ! vv=1,this%NumVars

      ! End Routine
      !--------------
      return
    end subroutine get_INV_data2D
    !==================================================================


    !==================================================================
    subroutine get_INV_data3D(this,VarName,Output )
      !
      !================================================================
      !
      ! Passed Variables
      !------------------
      class(INV_Data_t)           :: this
      character(len=*),intent(in ):: VarName
      real(r8)        ,intent(out):: Output(1:pcols,1:pver,begchunk:endchunk)
      !
      ! Local Values
      !--------------
      integer :: vv

      ! Loop thue the OBS variables and find the index for the given VarName
      !-----------------------------------------------------------------------
      do vv=1,this%NumVars
        if(trim(VarName).eq.trim(this%var(vv)%Name)) then

          ! Verify the the variable is HORIZ_ONLY
          !---------------------------------------
          if(this%var(vv)%HORIZ_ONLY) then
            call endrun('get_INV_data3D HORIZ_ONLY error ')
          endif

          ! copy the values to the output array
          !-------------------------------------
          Output(:,:,begchunk:endchunk) = this%var(vv)%val(:,:,begchunk:endchunk)
        endif
      end do ! vv=1,this%NumVars

      ! End Routine
      !--------------
      return
    end subroutine get_INV_data3D
    !==================================================================


    !==================================================================
    subroutine print_INV_Data(this)
      !
      !================================================================
      !
      ! Passed Variables
      !-----------------
      class(INV_Data_t):: this
      !
      ! Local Values
      !-------------
      integer:: nn

      ! Output the list of INV datasets, and the number of INV variables
      !--------------------------------------------------------------------
      if(masterproc) then
        write(iulog,*) 'INV_Data: '
        write(iulog,*) 'INV_Data: Input Datasets for INV variables:'//              &
                                                      ' Num_INV_input=',this%NumInput
        write(iulog,*) 'INV_Data:--------------------------------------------------- '
        do nn=1,this%NumInput
          write(iulog,*) 'INV_Data: Input dataset nn=',nn
          write(iulog,*) 'INV_Data:    Path_Template=',trim(this%input(nn)%Path_Template)
          write(iulog,*) 'INV_Data:    Name_Template=',trim(this%input(nn)%Name_Template)
          write(iulog,*) 'INV_Data:    dsID         =',     this%input(nn)%dsID
        end do
        write(iulog,*) 'INV_Data: '
        write(iulog,*) 'INV_Data: INV variables: Num_INV_Vars=',this%NumVars
        write(iulog,*) 'INV_Data:--------------------------------------------------- '
        do nn=1,this%NumVars
          write(iulog,*) 'INV_Data: INV Variable nn=',nn
          write(iulog,*) 'INV_Data:     Name='//trim(this%var(nn)%Name  )// &
                                        ' :  '//trim(this%var(nn)%DSname)
          write(iulog,*) 'INV_Data: idxINPUT =',     this%var(nn)%idxINPUT
          write(iulog,*) 'INV_Data:    dsID  =',     this%var(nn)%dsID
          write(iulog,*) 'INV_Data:    pathID=',     this%var(nn)%pathID
          write(iulog,*) 'INV_Data:   HmapOpt=',trim(this%var(nn)%HmapOpt)
          write(iulog,*) 'INV_Data:   VmapOpt=',trim(this%var(nn)%VmapOpt)
        end do
        write(iulog,*) 'INV_Data: '
      endif
  
      ! End Routine
      !--------------
      return
    end subroutine print_INV_Data
    !==================================================================


    !================================================================
    subroutine lookup_DataType(FileType,VarName,FileString,DSname,OBSstep)
      !
      !==============================================================
      ! 
      ! Passed variables 
      !-------------------
      character(len=* ),intent(in ):: FileType
      character(len=* ),intent(in ):: VarName
      character(len=cl),intent(out):: FileString
      character(len=16),intent(out):: DSname
      integer          ,intent(out):: OBSstep

      ! Fill in Nudge_Var datastructure for the given FileType
      !-------------------------------------------------------
      if(trim(FileType).eq.'ERA5') then
        OBSstep    = 3600   ! 1 hour
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
        OBSstep    = 6*3600 ! 6 hours
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
        OBSstep = -1
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
      else
        ! For unknown Type, set the FileString to NONE, the DSname equal to 
        ! the VarName, and set the OBSstep=-1 for the init routine to
        ! try and determine OBSstep from the dataset.
        !-------------------------------------------------------------
        OBSstep    = -1
        FileString = 'NONE'
        DSname     = trim(VarName)
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

end module OBSdata_mod
