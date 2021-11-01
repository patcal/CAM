module ncdio_atm

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !MODULE: ncdio_atm
  !
  ! !DESCRIPTION:
  ! Generic interfaces to write fields to PIO files
  !
  ! !USES:

  use pio, only: pio_offset_kind, file_desc_t, var_desc_t, pio_double,        &
       pio_inq_dimid, pio_max_var_dims, io_desc_t, pio_setframe
  use shr_kind_mod,     only: r8 => shr_kind_r8
  use shr_sys_mod,      only: shr_sys_flush      ! Standardized system subroutines
  use shr_scam_mod,     only: shr_scam_getCloseLatLon  ! Standardized system subroutines
  use spmd_utils,       only: masterproc
  use cam_abortutils,   only: endrun
  use scamMod,          only: scmlat,scmlon,single_column
  use cam_logfile,      only: iulog
  use string_utils,     only: to_lower
  use Dataset_mod,      only: dataset_t
  use Dataconnector_mod,only: register_dataPath, apply_dataPath
  use Hmaps_mod,        only: packed_hgrid_t
  use Vmaps_mod,        only: vgrid_t
  !
  ! !PUBLIC TYPES:
  implicit none

  PRIVATE

  save

  logical :: debug = .false.

  !
  !EOP
  !
  interface infld
     module procedure infld_real_1d_2d
     module procedure infld_real_2d_2d
     module procedure infld_real_2d_3d
     module procedure infld_real_3d_3d
     module procedure infld_dataset_3d
     module procedure infld_dataset_horiz
  end interface


  public :: infld
  public :: register_infld_grid
  public :: set_model_Hgrid
  public :: set_model_Vgrid
  private:: infld_dataset_3d
  private:: infld_dataset_horiz

  integer STATUS
  real(r8) surfdat
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: infld_real_1d_2d
  !
  ! !INTERFACE:
  subroutine infld_real_1d_2d(varname, ncid, dimname1,                        &
       dim1b, dim1e, dim2b, dim2e, field, readvar, gridname, timelevel,       &
       fillvalue)
    !
    ! !DESCRIPTION:
    ! Netcdf I/O of initial real field from netCDF file
    ! Read a 1-D field (or slice) into a 2-D variable
    !
    ! !USES
    !

    use pio,              only: pio_get_var, pio_read_darray, pio_setdebuglevel
    use pio,              only: PIO_MAX_NAME, pio_inquire, pio_inq_dimname
    use cam_grid_support, only: cam_grid_check, cam_grid_get_decomp, cam_grid_id, &
                                cam_grid_dimensions
    use cam_pio_utils,    only: cam_pio_check_var, cam_pio_inq_var_fill

    !
    ! !ARGUMENTS:
    implicit none
    character(len=*),  intent(in)     :: varname  ! variable name
    type(file_desc_t), intent(inout)  :: ncid     ! input unit
    character(len=*),  intent(in)     :: dimname1 ! name of 1st array dimensions of field on file (array order)
    integer,           intent(in)     :: dim1b    ! start of first  dimension of array to be returned
    integer,           intent(in)     :: dim1e    ! end   of first  dimension of array to be returned
    integer,           intent(in)     :: dim2b    ! start of second dimension of array to be returned
    integer,           intent(in)     :: dim2e    ! end   of second dimension of array to be returned
    real(r8), target,  intent(out)    :: field(dim1b:dim1e,dim2b:dim2e) ! array to be returned (decomposed or global)
    logical,           intent(out)    :: readvar  ! true => variable is on initial dataset
    character(len=*), optional, intent(in) :: gridname ! Name of variable's grid
    integer, optional, intent(in)     :: timelevel
    real(r8), optional, intent(out)   :: fillvalue
    !
    !EOP
    !
    ! !LOCAL VARIABLES:
    type(io_desc_t), pointer  :: iodesc
    integer                   :: grid_id   ! grid ID for data mapping
    integer                   :: i, j      ! indices
    integer                   :: ierr      ! error status
    type(var_desc_t)          :: varid     ! variable id
    integer                   :: no_fill
    integer                   :: arraydimsize(2) ! field dimension lengths

    integer                   :: ndims ! number of dimensions
    integer                   :: dimids(PIO_MAX_VAR_DIMS) ! file variable dims
    integer                   :: dimlens(PIO_MAX_VAR_DIMS) ! file variable shape
    integer                   :: grid_dimlens(2)

    ! Offsets for reading global variables
    integer                   :: strt(1) = 1 ! start ncol index for netcdf 1-d
    integer                   :: cnt (1) = 1 ! ncol count for netcdf 1-d
    character(len=PIO_MAX_NAME) :: tmpname
    character(len=128)        :: errormsg

    logical                   :: readvar_tmp ! if true, variable is on tape
    character(len=*), parameter :: subname='INFLD_REAL_1D_2D' ! subroutine name

    ! For SCAM
    real(r8)                  :: closelat, closelon
    integer                   :: lonidx, latidx

    nullify(iodesc)

    !
    !-----------------------------------------------------------------------
    !
    !    call pio_setdebuglevel(3)

    !
    ! Error conditions
    !
    if (present(gridname)) then
      grid_id = cam_grid_id(trim(gridname))
    else
      grid_id = cam_grid_id('physgrid')
    end if
    if (.not. cam_grid_check(grid_id)) then
      if(masterproc) then
        if (present(gridname)) then
          write(errormsg, *)': invalid gridname, "',trim(gridname),'", specified for field ',trim(varname)
        else
          write(errormsg, *)': Internal error, no "physgrid" gridname'
        end if
      end if
      call endrun(trim(subname)//errormsg)
    end if

    ! Get the number of columns in the global grid.
    call cam_grid_dimensions(grid_id, grid_dimlens)

    if (debug .and. masterproc) then
      if (present(gridname)) then
        write(iulog, '(5a)') trim(subname),': field = ',trim(varname),', grid = ',trim(gridname)
      else
        write(iulog, '(4a)') trim(subname),': field = ',trim(varname),', grid = physgrid'
      end if
      call shr_sys_flush(iulog)
    end if
    !
    ! Read netCDF file
    !
    !
    ! Check if field is on file; get netCDF variable id
    !
    call cam_pio_check_var(ncid, varname, varid, ndims, dimids, dimlens, readvar_tmp)
    !
    ! If field is on file:
    !
    if (readvar_tmp) then
        if (debug .and. masterproc) then
          write(iulog, '(2a,5(i0,a))') trim(subname),': field(',              &
               dim1b,':',dim1e,',',dim2b,':',dim2e, '), file(',dimlens(1),')'
          call shr_sys_flush(iulog)
        end if
      !
      ! Get array dimension id's and sizes
      !
      arraydimsize(1) = (dim1e - dim1b + 1)
      arraydimsize(2) = (dim2e - dim2b + 1)
      do j = 1, 2
        if (arraydimsize(j) /= size(field, j)) then
          write(errormsg, *) ': Mismatch between array bounds and field size for ', &
               trim(varname), ', dimension', j
          call endrun(trim(subname)//errormsg)
        end if
      end do

      if (ndims > 2) then
        call endrun(trim(subname)//': too many dimensions for '//trim(varname))
      else if (ndims < 1) then
        call endrun(trim(subname)//': too few dimensions for '//trim(varname))
      else
         ! Check that the number of columns in the file matches the number of
         ! columns in the grid object.
         if (dimlens(1) /= grid_dimlens(1)) then
            readvar = .false.
            return
         end if

        ! Check to make sure that the second dimension is time
        if (ndims == 2) then
          ierr = pio_inq_dimname(ncid, dimids(2), tmpname)
          if (trim(tmpname) /= 'time') then
            call endrun(trim(subname)//': dimension mismatch for '//trim(varname))
          end if
        end if
      end if

      if(ndims == 2) then
        if(present(timelevel)) then
          call pio_setframe(ncid, varid, int(timelevel,kind=pio_offset_kind))
        else
          call pio_setframe(ncid, varid, int(1,kind=pio_offset_kind))
        end if
        ndims = ndims - 1
      end if

      ! NB: strt and cnt were initialized to 1
      if (single_column) then
        !!XXgoldyXX: Clearly, this will not work for an unstructured dycore
        call endrun(trim(subname)//': SCAM not supported in this configuration')
      else
        ! All distributed array processing
        call cam_grid_get_decomp(grid_id, arraydimsize, dimlens(1:ndims),    &
             pio_double, iodesc)
        call pio_read_darray(ncid, varid, iodesc, field, ierr)
        if (present(fillvalue)) then
           ierr = cam_pio_inq_var_fill(ncid, varid, fillvalue)
        end if
     end if


      if (masterproc) write(iulog,*) subname//': read field '//trim(varname)

    end if  ! end of readvar_tmp

    readvar = readvar_tmp

    return

  end subroutine infld_real_1d_2d

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: infld_real_2d_2d
  !
  ! !INTERFACE:
  subroutine infld_real_2d_2d(varname, ncid, dimname1, dimname2,              &
       dim1b, dim1e, dim2b, dim2e, field, readvar, gridname, timelevel,       &
       fillvalue)
    !
    ! !DESCRIPTION:
    ! Netcdf I/O of initial real field from netCDF file
    ! Read a 2-D field (or slice) into a 2-D variable
    !
    ! !USES
    !

    use pio,              only: pio_get_var, pio_read_darray, pio_setdebuglevel
    use pio,              only: PIO_MAX_NAME, pio_inquire, pio_inq_dimname
    use cam_grid_support, only: cam_grid_check, cam_grid_get_decomp, cam_grid_id
    use cam_pio_utils,    only: cam_permute_array, calc_permutation
    use cam_pio_utils,    only: cam_pio_check_var, cam_pio_inq_var_fill

    !
    ! !ARGUMENTS:
    implicit none
    character(len=*),  intent(in)     :: varname  ! variable name
    type(file_desc_t), intent(inout)  :: ncid     ! input unit
    character(len=*),  intent(in)     :: dimname1 ! name of 1st array dimensions of field on file (array order)
    character(len=*),  intent(in)     :: dimname2 ! name of 2nd array dimensions of field on file (array order)
    integer,           intent(in)     :: dim1b    ! start of first  dimension of array to be returned
    integer,           intent(in)     :: dim1e    ! end   of first  dimension of array to be returned
    integer,           intent(in)     :: dim2b    ! start of second dimension of array to be returned
    integer,           intent(in)     :: dim2e    ! end   of second dimension of array to be returned
    real(r8), target,  intent(out)    :: field(dim1b:dim1e,dim2b:dim2e) ! array to be returned (decomposed or global)
    logical,           intent(out)    :: readvar  ! true => variable is on initial dataset
    character(len=*), optional, intent(in) :: gridname ! Name of variable's grid
    integer, optional, intent(in)     :: timelevel
    real(r8), optional, intent(out)   :: fillvalue
    !
    !EOP
    !
    ! !LOCAL VARIABLES:
    type(io_desc_t), pointer  :: iodesc
    integer                   :: grid_id   ! grid ID for data mapping
    integer                   :: i, j      ! indices
    integer                   :: ierr      ! error status
    type(var_desc_t)          :: varid     ! variable id

    integer                   :: arraydimsize(2) ! field dimension lengths
    integer                   :: arraydimids(2) ! Dimension IDs
    integer                   :: permutation(2)
    logical                   :: ispermuted

    integer                   :: ndims ! number of dimensions on file
    integer                   :: dimids(PIO_MAX_VAR_DIMS) ! file variable dims
    integer                   :: dimlens(PIO_MAX_VAR_DIMS) ! file variable shape

    ! Offsets for reading global variables
    integer                   :: strt(2) ! start lon, lat indices for netcdf 2-d
    integer                   :: cnt (2) ! lon, lat counts for netcdf 2-d
    character(len=PIO_MAX_NAME) :: tmpname
    character(len=128)        :: errormsg

    real(r8), pointer         :: tmp2d(:,:) ! input data for permutation

    logical                   :: readvar_tmp ! if true, variable is on tape
    character(len=*), parameter :: subname='INFLD_REAL_2D_2D' ! subroutine name
    character(len=PIO_MAX_NAME) :: field_dnames(2)

    ! For SCAM
    real(r8)                  :: closelat, closelon
    integer                   :: lonidx, latidx

    nullify(iodesc)

    !
    !-----------------------------------------------------------------------
    !
    !    call pio_setdebuglevel(3)

    ! Should we be using a different interface?
    if ((trim(dimname1) == trim(dimname2)) .or. (len_trim(dimname2) == 0)) then
      call infld(varname, ncid, dimname1, dim1b, dim1e, dim2b, dim2e,         &
           field, readvar, gridname, timelevel)
    else

      !
      ! Error conditions
      !
      if (present(gridname)) then
        grid_id = cam_grid_id(trim(gridname))
      else
        grid_id = cam_grid_id('physgrid')
      end if
      if (.not. cam_grid_check(grid_id)) then
        if(masterproc) then
          if (present(gridname)) then
            write(errormsg, *)': invalid gridname, "',trim(gridname),'", specified for field ',trim(varname)
          else
            write(errormsg, *)': Internal error, no "physgrid" gridname'
          end if
        end if
        call endrun(trim(subname)//errormsg)
      end if

    if (debug .and. masterproc) then
      if (present(gridname)) then
        write(iulog, '(5a)') trim(subname),': field = ',trim(varname),', grid = ',trim(gridname)
      else
        write(iulog, '(4a)') trim(subname),': field = ',trim(varname),', grid = physgrid'
      end if
      call shr_sys_flush(iulog)
    end if

      !
      ! Read netCDF file
      !
      !
      ! Check if field is on file; get netCDF variable id
      !
      call cam_pio_check_var(ncid, varname, varid, ndims, dimids, dimlens, readvar_tmp)
      !
      ! If field is on file:
      !
      if (readvar_tmp) then
        if (debug .and. masterproc) then
          write(iulog, '(2a,6(i0,a))') trim(subname),': field(',              &
               dim1b,':',dim1e,',',dim2b,':',dim2e,                           &
               '), file(',dimlens(1),',',dimlens(2),')'
          call shr_sys_flush(iulog)
        end if
        !
        ! Get array dimension id's and sizes
        !
        ierr = PIO_inq_dimid(ncid, dimname1, arraydimids(1))
        ierr = PIO_inq_dimid(ncid, dimname2, arraydimids(2))
        arraydimsize(1) = (dim1e - dim1b + 1)
        arraydimsize(2) = (dim2e - dim2b + 1)
        do j = 1, 2
          if (arraydimsize(j) /= size(field, j)) then
            write(errormsg, *) ': Mismatch between array bounds and field size for ', &
                 trim(varname), ', dimension', j
            call endrun(trim(subname)//errormsg)
          end if
        end do

        if (ndims > 3) then
          call endrun(trim(subname)//': too many dimensions for '//trim(varname))
        else if (ndims < 2) then
          call endrun(trim(subname)//': too few dimensions for '//trim(varname))
        else
          ! Check to make sure that the third dimension is time
          if (ndims == 3) then
            ierr = pio_inq_dimname(ncid, dimids(3), tmpname)
            if (trim(tmpname) /= 'time') then
              call endrun(trim(subname)//': dimension mismatch for '//trim(varname))
            end if
          end if
        end if

        if(ndims == 3) then
          if(present(timelevel)) then
            call pio_setframe(ncid, varid, int(timelevel,kind=pio_offset_kind))
          else
            call pio_setframe(ncid, varid, int(1,kind=pio_offset_kind))
          end if
        end if

        field_dnames(1) = dimname1
        field_dnames(2) = dimname2
        if (single_column) then
          ! This could be generalized but for now only handles a single point
          strt(1) = dim1b
          strt(2) = dim2b
          cnt = arraydimsize
          call shr_scam_getCloseLatLon(ncid,scmlat,scmlon,closelat,closelon,latidx,lonidx)
          if (trim(field_dnames(1)) == 'lon') then
            strt(1) = lonidx ! First dim always lon for Eulerian dycore
          else
            call endrun(trim(subname)//': lon should be first dimension for '//trim(varname))
          end if
          if (trim(field_dnames(2)) == 'lat') then
            strt(2) = latidx
          else
            call endrun(trim(subname)//': lat dimension not found for '//trim(varname))
          end if

          ! Check for permuted dimensions ('out of order' array)
          call calc_permutation(dimids, arraydimids, permutation, ispermuted)
          if (ispermuted) then
            call cam_permute_array(strt, permutation)
            call cam_permute_array(cnt, permutation)
            allocate(tmp2d(1:cnt(1), 1:cnt(2)))
            ierr = pio_get_var(ncid, varid, strt, cnt, tmp2d)
            do j = dim2b, dim2e
              do i = dim1b, dim1e
                ! We don't need strt anymore, reuse it
                strt(1) = i - dim1b + 1
                strt(2) = j - dim2b + 1
                call cam_permute_array(strt, permutation)
                field(i,j) = tmp2d(strt(1), strt(2))
              end do
            end do
            deallocate(tmp2d)
          else
            ierr = pio_get_var(ncid, varid, strt, cnt, field)
          end if
        else
          ! All distributed array processing
          call cam_grid_get_decomp(grid_id, arraydimsize, dimlens(1:2),      &
               pio_double, iodesc, field_dnames=field_dnames)
          call pio_read_darray(ncid, varid, iodesc, field, ierr)
          if (present(fillvalue)) then
             ierr = cam_pio_inq_var_fill(ncid, varid, fillvalue)
          end if
        end if

        if (masterproc) write(iulog,*) subname//': read field '//trim(varname)

      end if  ! end of readvar_tmp

      readvar = readvar_tmp

    end if ! end of call infld_real_1d_2d instead

  end subroutine infld_real_2d_2d


  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: infld_real_2d_3d
  !
  ! !INTERFACE:
  subroutine infld_real_2d_3d(varname, ncid, dimname1, dimname2,              &
       dim1b, dim1e, dim2b, dim2e, dim3b, dim3e,                              &
       field, readvar, gridname, timelevel, fillvalue)
    !
    ! !DESCRIPTION:
    ! Netcdf I/O of initial real field from netCDF file
    ! Read a 2-D field (or slice) into a 3-D variable
    !
    ! !USES
    !

    use pio,              only: pio_get_var, pio_read_darray, pio_setdebuglevel
    use pio,              only: PIO_MAX_NAME, pio_inquire, pio_inq_dimname
    use cam_grid_support, only: cam_grid_check, cam_grid_get_decomp, cam_grid_id, &
                                cam_grid_dimensions
    use cam_pio_utils,    only: cam_permute_array, calc_permutation
    use cam_pio_utils,    only: cam_pio_check_var, cam_pio_inq_var_fill

    !
    ! !ARGUMENTS:
    implicit none
    character(len=*),  intent(in)     :: varname  ! variable name
    type(file_desc_t), intent(inout)  :: ncid     ! input unit
    character(len=*),  intent(in)     :: dimname1 ! name of 1st array dimensions of field on file (array order)
    character(len=*),  intent(in)     :: dimname2 ! name of 2nd array dimensions of field on file (array order)
    integer,           intent(in)     :: dim1b    ! start of first  dimension of array to be returned
    integer,           intent(in)     :: dim1e    ! end   of first  dimension of array to be returned
    integer,           intent(in)     :: dim2b    ! start of second dimension of array to be returned
    integer,           intent(in)     :: dim2e    ! end   of second dimension of array to be returned
    integer,           intent(in)     :: dim3b    ! start of third  dimension of array to be returned
    integer,           intent(in)     :: dim3e    ! end   of third  dimension of array to be returned
    real(r8), target,  intent(out)    :: field(dim1b:dim1e,dim2b:dim2e,dim3b:dim3e) ! array to be returned (decomposed or global)
    logical,           intent(out)    :: readvar  ! true => variable is on initial dataset
    character(len=*), optional, intent(in) :: gridname ! Name of variable's grid
    integer, optional, intent(in)     :: timelevel
    real(r8), optional, intent(out)   :: fillvalue
    !
    !EOP
    !
    ! !LOCAL VARIABLES:
    type(io_desc_t), pointer  :: iodesc
    integer                   :: grid_id   ! grid ID for data mapping
    integer                   :: i, j, k   ! indices
    integer                   :: ierr      ! error status
    type(var_desc_t)          :: varid     ! variable id

    integer                   :: arraydimsize(3) ! field dimension lengths
    integer                   :: arraydimids(2) ! Dimension IDs
    integer                   :: permutation(2)
    logical                   :: ispermuted

    integer                   :: ndims ! number of dimensions
    integer                   :: dimids(PIO_MAX_VAR_DIMS) ! file variable dims
    integer                   :: dimlens(PIO_MAX_VAR_DIMS) ! file variable shape
    integer                   :: grid_dimlens(2)

    ! Offsets for reading global variables
    integer                   :: strt(3) = 1 ! start ncol, lev indices for netcdf 2-d
    integer                   :: cnt (3) = 1 ! ncol, lev counts for netcdf 2-d
    character(len=PIO_MAX_NAME) :: tmpname

    real(r8), pointer         :: tmp3d(:,:,:) ! input data for permutation

    logical                   :: readvar_tmp ! if true, variable is on tape
    character(len=*), parameter :: subname='INFLD_REAL_2D_3D' ! subroutine name
    character(len=128)        :: errormsg
    character(len=PIO_MAX_NAME) :: field_dnames(2)
    character(len=PIO_MAX_NAME) :: file_dnames(3)

    ! For SCAM
    real(r8)                  :: closelat, closelon
    integer                   :: lonidx, latidx

    nullify(iodesc)

    !
    !-----------------------------------------------------------------------
    !
    !    call pio_setdebuglevel(3)

    !
    ! Error conditions
    !
    if (present(gridname)) then
      grid_id = cam_grid_id(trim(gridname))
    else
      grid_id = cam_grid_id('physgrid')
    end if
    if (.not. cam_grid_check(grid_id)) then
      if(masterproc) then
        if (present(gridname)) then
          write(errormsg, *)': invalid gridname, "',trim(gridname),'", specified for field ',trim(varname)
        else
          write(errormsg, *)': Internal error, no "physgrid" gridname'
        end if
      end if
      call endrun(trim(subname)//errormsg)
    end if

    ! Get the number of columns in the global grid.
    call cam_grid_dimensions(grid_id, grid_dimlens)

    if (debug .and. masterproc) then
      if (present(gridname)) then
        write(iulog, '(5a)') trim(subname),': field = ',trim(varname),', grid = ',trim(gridname)
      else
        write(iulog, '(4a)') trim(subname),': field = ',trim(varname),', grid = physgrid'
      end if
      call shr_sys_flush(iulog)
    end if

    !
    ! Read netCDF file
    !
    !
    ! Check if field is on file; get netCDF variable id
    !
    call cam_pio_check_var(ncid, varname, varid, ndims, dimids, dimlens, &
       readvar_tmp, dimnames=file_dnames)

    ! If field is on file:
    !
    if (readvar_tmp) then
      if (debug .and. masterproc) then
        write(iulog, '(2a,8(i0,a))') trim(subname),': field(',                &
             dim1b,':',dim1e,',',dim2b,':',dim2e,',',dim3b,':',dim3e,         &
             '), file(',dimlens(1),',',dimlens(2),')'
        call shr_sys_flush(iulog)
      end if
      !
      ! Get array dimension id's and sizes
      !
      arraydimsize(1) = (dim1e - dim1b + 1)
      arraydimsize(2) = (dim2e - dim2b + 1)
      arraydimsize(3) = (dim3e - dim3b + 1)
      do j = 1, 3
        if (arraydimsize(j) /= size(field, j)) then
          write(errormsg, *) ': Mismatch between array bounds and field size for ', &
               trim(varname), ', dimension', j
          call endrun(trim(subname)//errormsg)
        end if
      end do

      if (ndims > 3) then
        call endrun(trim(subname)//': too many dimensions for '//trim(varname))
      else if (ndims < 2) then
        call endrun(trim(subname)//': too few dimensions for '//trim(varname))
      else
         ! Check that the number of columns in the file matches the number of
         ! columns in the grid object.
         if (dimlens(1) /= grid_dimlens(1) .and. dimlens(2) /= grid_dimlens(1)) then
            readvar = .false.
            return
         end if

        ! Check to make sure that the 3rd dimension is time
        if (ndims == 3) then
          ierr = pio_inq_dimname(ncid, dimids(3), tmpname)
          if (to_lower(trim(tmpname)) /= 'time') then
            call endrun(trim(subname)//': dimension mismatch for '//trim(varname))
          end if
        end if
      end if

      if(ndims == 3) then
        if(present(timelevel)) then
          call pio_setframe(ncid, varid, int(timelevel,kind=pio_offset_kind))
        else
          call pio_setframe(ncid, varid, int(1,kind=pio_offset_kind))
        end if
        ndims = ndims - 1
      end if

      field_dnames(1) = dimname1
      field_dnames(2) = dimname2
      ! NB: strt and cnt were initialized to 1
      if (single_column) then
        !!XXgoldyXX: Clearly, this will not work for an unstructured dycore
        ! Check for permuted dimensions ('out of order' array)
!       call calc_permutation(dimids(1:2), arraydimids, permutation, ispermuted)
        call endrun(trim(subname)//': SCAM not supported in this configuration')
      else
        ! All distributed array processing
        call cam_grid_get_decomp(grid_id, arraydimsize, dimlens(1:2),         &
             pio_double, iodesc, field_dnames=field_dnames,                   &
             file_dnames=file_dnames(1:2))
        call pio_read_darray(ncid, varid, iodesc, field, ierr)
          if (present(fillvalue)) then
             ierr = cam_pio_inq_var_fill(ncid, varid, fillvalue)
          end if
      end if

      if (masterproc) write(iulog,*) subname//': read field '//trim(varname)

    end if  ! end of readvar_tmp

    readvar = readvar_tmp

    return

  end subroutine infld_real_2d_3d

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: infld_real_3d_3d
  !
  ! !INTERFACE:
  subroutine infld_real_3d_3d(varname, ncid, dimname1, dimname2, dimname3,    &
       dim1b, dim1e, dim2b, dim2e, dim3b, dim3e,                              &
       field, readvar, gridname, timelevel, fillvalue)
    !
    ! !DESCRIPTION:
    ! Netcdf I/O of initial real field from netCDF file
    ! Read a 3-D field (or slice) into a 3-D variable
    !
    ! !USES
    !

    use pio,              only: pio_get_var, pio_read_darray, pio_setdebuglevel
    use pio,              only: PIO_MAX_NAME, pio_inquire, pio_inq_dimname
    use cam_grid_support, only: cam_grid_check, cam_grid_get_decomp, cam_grid_id
    use cam_pio_utils,    only: cam_permute_array, calc_permutation
    use cam_pio_utils,    only: cam_pio_check_var, cam_pio_inq_var_fill

    !
    ! !ARGUMENTS:
    implicit none
    character(len=*),  intent(in)     :: varname  ! variable name
    type(file_desc_t), intent(inout)  :: ncid     ! input unit
    character(len=*),  intent(in)     :: dimname1 ! name of 1st array dimensions of field on file (array order)
    character(len=*),  intent(in)     :: dimname2 ! name of 2nd array dimensions of field on file (array order)
    character(len=*),  intent(in)     :: dimname3 ! name of 3rd array dimensions of field on file (array order)
    integer,           intent(in)     :: dim1b    ! start of first  dimension of array to be returned
    integer,           intent(in)     :: dim1e    ! end   of first  dimension of array to be returned
    integer,           intent(in)     :: dim2b    ! start of second dimension of array to be returned
    integer,           intent(in)     :: dim2e    ! end   of second dimension of array to be returned
    integer,           intent(in)     :: dim3b    ! start of third  dimension of array to be returned
    integer,           intent(in)     :: dim3e    ! end   of third  dimension of array to be returned
    real(r8), target,  intent(out)    :: field(dim1b:dim1e,dim2b:dim2e,dim3b:dim3e) ! array to be returned (decomposed or global)
    logical,           intent(out)    :: readvar  ! true => variable is on initial dataset
    character(len=*), optional, intent(in) :: gridname ! Name of variable's grid
    integer, optional, intent(in)     :: timelevel
    real(r8), optional, intent(out)   :: fillvalue
    !
    !EOP
    !
    ! !LOCAL VARIABLES:
    type(io_desc_t), pointer  :: iodesc
    integer                   :: grid_id   ! grid ID for data mapping
    integer                   :: i, j, k   ! indices
    integer                   :: ierr      ! error status
    type(var_desc_t)          :: varid     ! variable id

    integer                   :: arraydimsize(3) ! field dimension lengths
    integer                   :: arraydimids(3) ! Dimension IDs
    integer                   :: permutation(3)
    logical                   :: ispermuted

    integer                   :: ndims ! number of dimensions
    integer                   :: pdims ! number of dimensions w/o timeslice
    integer                   :: dimids(PIO_MAX_VAR_DIMS) ! file variable dims
    integer                   :: dimlens(PIO_MAX_VAR_DIMS) ! file variable shape

    ! Offsets for reading global variables
    integer                   :: strt(3)     ! start lon, lev, lat indices for netcdf 3-d
    integer                   :: cnt (3)     ! lon, lat counts for netcdf 3-d
    character(len=PIO_MAX_NAME) :: tmpname

    real(r8), pointer         :: tmp3d(:,:,:) ! input data for permutation

    logical                   :: readvar_tmp ! if true, variable is on tape
    character(len=*), parameter :: subname='INFLD_REAL_3D_3D' ! subroutine name
    character(len=128)        :: errormsg
    character(len=PIO_MAX_NAME) :: field_dnames(3)
    character(len=PIO_MAX_NAME) :: file_dnames(4)

    ! For SCAM
    real(r8)                  :: closelat, closelon
    integer                   :: lonidx, latidx

    nullify(iodesc)

    !
    !-----------------------------------------------------------------------
    !
    !    call pio_setdebuglevel(3)

    ! Should we be using a different interface?
    if ((trim(dimname1) == trim(dimname2)) .or. (len_trim(dimname2) == 0)) then
      call infld(varname, ncid, dimname1, dimname3,                           &
           dim1b, dim1e, dim2b, dim2e, dim3b, dim3e,                          &
           field, readvar, gridname, timelevel)
    else if ((trim(dimname1) == trim(dimname3)) .or. (len_trim(dimname3) == 0)) then
      call infld(varname, ncid, dimname1, dimname2,                           &
           dim1b, dim1e, dim2b, dim2e, dim3b, dim3e,                          &
           field, readvar, gridname, timelevel)
    else

      !
      ! Error conditions
      !
      if (present(gridname)) then
        grid_id = cam_grid_id(trim(gridname))
      else
        grid_id = cam_grid_id('physgrid')
      end if
      if (.not. cam_grid_check(grid_id)) then
        if(masterproc) then
          if (present(gridname)) then
            write(errormsg, *)': invalid gridname, "',trim(gridname),'", specified for field ',trim(varname)
          else
            write(errormsg, *)': Internal error, no "physgrid" gridname'
          end if
        end if
        call endrun(trim(subname)//errormsg)
      end if

      if (debug .and. masterproc) then
        if (present(gridname)) then
          write(iulog, '(5a)') trim(subname),': field = ',trim(varname),', grid = ',trim(gridname)
        else
          write(iulog, '(4a)') trim(subname),': field = ',trim(varname),', grid = physgrid'
        end if
        call shr_sys_flush(iulog)
      end if

      !
      ! Read netCDF file
      !
      !
      ! Check if field is on file; get netCDF variable id
      !
      call cam_pio_check_var(ncid, varname, varid, ndims, dimids, dimlens,    &
           readvar_tmp, dimnames=file_dnames)
      !
      ! If field is on file:
      !
      if (readvar_tmp) then
        if (debug .and. masterproc) then
          write(iulog, '(2a,9(i0,a))') trim(subname),': field(',              &
               dim1b,':',dim1e,',',dim2b,':',dim2e,',',dim3b,':',dim3e,       &
               '), file(',dimlens(1),',',dimlens(2),',',dimlens(3),')'
          call shr_sys_flush(iulog)
        end if
        !
        ! Get array dimension id's and sizes
        !
        ierr = PIO_inq_dimid(ncid, dimname1, arraydimids(1))
        ierr = PIO_inq_dimid(ncid, dimname2, arraydimids(2))
        ierr = PIO_inq_dimid(ncid, dimname3, arraydimids(3))
        arraydimsize(1) = (dim1e - dim1b + 1)
        arraydimsize(2) = (dim2e - dim2b + 1)
        arraydimsize(3) = (dim3e - dim3b + 1)

        do j = 1, 3
          if (arraydimsize(j) /= size(field, j)) then
            write(errormsg, *) ': Mismatch between array bounds and field size for ', &
                 trim(varname), ', dimension', j
            call endrun(trim(subname)//errormsg)
          end if
        end do

        pdims = ndims
        if (ndims > 4) then
          call endrun(trim(subname)//': too many dimensions for '//trim(varname))
        else if (ndims < 3) then
          call endrun(trim(subname)//': too few dimensions for '//trim(varname))
        else
          ! Check to make sure that the fourth dimension is time
          if (ndims == 4) then
            ierr = pio_inq_dimname(ncid, dimids(4), tmpname)
            if (trim(tmpname) /= 'time') then
              call endrun(trim(subname)//': dimension mismatch for '//trim(varname))
            end if
            pdims = 3
          end if
        end if

        if(ndims == 4) then
          if(present(timelevel)) then
            call pio_setframe(ncid, varid, int(timelevel,kind=pio_offset_kind))
          else
            call pio_setframe(ncid, varid, int(1,kind=pio_offset_kind))
          end if
        end if

        field_dnames(1) = dimname1
        field_dnames(2) = dimname2
        field_dnames(3) = dimname3

        if (single_column) then
          ! This could be generalized but for now only handles a single point
          strt(1) = dim1b
          strt(2) = dim2b
          strt(3) = dim3b
          cnt = arraydimsize
          call shr_scam_getCloseLatLon(ncid,scmlat,scmlon,closelat,closelon,latidx,lonidx)
          if (trim(field_dnames(1)) == 'lon') then
            strt(1) = lonidx ! First dim always lon for Eulerian dycore
          else
            call endrun(trim(subname)//': lon should be first dimension for '//trim(varname))
          end if
          if (trim(field_dnames(2)) == 'lat') then
            strt(2) = latidx
          else if (trim(field_dnames(3)) == 'lat') then
            strt(3) = latidx
          else
            call endrun(trim(subname)//': lat dimension not found for '//trim(varname))
          end if

          ! Check for permuted dimensions ('out of order' array)
          call calc_permutation(dimids, arraydimids, permutation, ispermuted)
          if (ispermuted) then
            call cam_permute_array(strt, permutation)
            call cam_permute_array(cnt, permutation)
            allocate(tmp3d(1:cnt(1), 1:cnt(2), 1:cnt(3)))
            ierr = pio_get_var(ncid, varid, strt, cnt, tmp3d)
            do k = dim3b, dim3e
              do j = dim2b, dim2e
                do i = dim1b, dim1e
                  ! We don't need strt anymore, reuse it
                  strt(1) = i - dim1b + 1
                  strt(2) = j - dim2b + 1
                  strt(3) = k - dim3b + 1
                  call cam_permute_array(strt, permutation)
                  field(i,j,k) = tmp3d(strt(1), strt(2), strt(3))
                end do
              end do
            end do
            deallocate(tmp3d)
          else
            ierr = pio_get_var(ncid, varid, strt, cnt, field)
          end if
        else
          ! All distributed array processing
          call cam_grid_get_decomp(grid_id, arraydimsize, dimlens(1:pdims),   &
               pio_double, iodesc, field_dnames=field_dnames,                 &
               file_dnames=file_dnames(1:3))
          call pio_read_darray(ncid, varid, iodesc, field, ierr)
          if (present(fillvalue)) then
             ierr = cam_pio_inq_var_fill(ncid, varid, fillvalue)
          end if
        end if ! end of single column

        if (masterproc) write(iulog,*) subname//': read field '//trim(varname)

      end if  ! end of readvar_tmp

      readvar = readvar_tmp

    end if ! end of call infld_real_2d_3d instead

  end subroutine infld_real_3d_3d


  !================================================================
  subroutine register_infld_grid(DataSetID,PathID,I_HgridName,I_HmapOpt,I_VgridName,I_VmapOpt)
    !
    ! register_infld_grid: Register the given Horizontal/Vertical grid and
    !                      mapping options for the given DataSet ID, and
    !                      return the index for the correponding DataPath.
    !==============================================================
    use cam_grid_support,only: cam_grid_get_dim_names, cam_grid_get_array_bounds, cam_grid_id
    !
    ! Passed variables
    !-------------------
    integer                  ,intent(in ):: DataSetID
    integer                  ,intent(out):: PathID
    character(len=*)         ,intent(in ):: I_HgridName
    character(len=*)         ,intent(in ):: I_HmapOpt
    character(len=*),optional,intent(in ):: I_VgridName
    character(len=*),optional,intent(in ):: I_VmapOpt
    !
    ! Local Values
    !--------------
    type(packed_hgrid_t):: Hgrid_m
    type(vgrid_t       ):: Vgrid_m
    character(len=32)   :: HmapOpt
    character(len=32)   :: VmapOpt
    character(len=32)   :: HgridName
    character(len=32)   :: VgridName

    HgridName = '                                '
    HmapOpt   = '                                '
    HgridName = trim(I_HgridName)
    HmapOpt   = trim(I_HmapOpt)

    if(present(I_VgridName).and.present(I_VmapOpt)) then
      VgridName = '                                '
      VmapOpt   = '                                '
      VgridName = trim(I_VgridName)
      VmapOpt   = trim(I_VmapOpt)
    else
      VgridName = '                                '
      VmapOpt   = '                                '
      VgridName = 'HORIZ_ONLY'
      VmapOpt   = 'HORIZ_ONLY'
    endif

    ! Initialize Hgrid for the given model grid
    !--------------------------------------------
    call set_model_Hgrid(Hgrid_m,trim(HgridName))

    ! If a VgridName and VmapOpt are present, use them,
    ! otherwise set the vertical to HORIZ_ONLY
    !---------------------------------------------------
    if(present(I_VgridName).and.present(I_VmapOpt)) then
      call set_model_Vgrid(Vgrid_m,trim(VgridName))
      PathID = register_dataPath(DataSetID,Hgrid_m,trim(HmapOpt),Vgrid_m,trim(VmapOpt))
    else
      call set_model_Vgrid(Vgrid_m,'HORIZ_ONLY')
      PathID = register_dataPath(DataSetID,Hgrid_m,trim(HmapOpt),Vgrid_m,'HORIZ_ONLY')
    endif

    ! End Routine
    !------------------
    return
  end subroutine register_infld_grid
  !================================================================


  !================================================================
  subroutine set_model_Hgrid(Hgrid_m,hgridname)
    !
    ! set_model_Hgrid: For the given hgridname, determine the packed set
    !                  of unique gridpoints for the current PE and
    !                  mapping indices back to model arrays. Then
    !                  initialize the given model Hgrid datastructure
    !==============================================================
    use cam_grid_support,only: cam_grid_id, cam_grid_check, cam_grid_get_gcid
    use cam_grid_support,only: cam_grid_get_array_bounds,cam_grid_get_local_size
    use cam_grid_support,only: cam_grid_get_latvals, cam_grid_get_lonvals
    use cam_grid_support,only: cam_grid_get_dim_names
    use pio,             only: iMap=>PIO_OFFSET_KIND
    !
    ! Passed variables
    !-------------------
    type(packed_hgrid_t),intent(inout):: Hgrid_m
    character(len=*)    ,intent(in   ):: hgridname
    !
    ! Local Values
    !----------------
    integer(iMap),pointer:: gcid(:)    => null()
    real(r8)     ,pointer:: Lonvals(:) => null()
    real(r8)     ,pointer:: Latvals(:) => null()
    integer:: gridid
    integer:: griddims(2,2)
    integer:: PEsize
    integer:: Lonsize
    integer:: Latsize

    integer             :: G_ncol
    real(r8),allocatable:: G_Lon (:)
    real(r8),allocatable:: G_Lat (:)
    integer ,allocatable:: G_Xmap(:)
    integer ,allocatable:: G_Ymap(:)

    integer:: pcols,begchunk,endchunk,nchunks,ichunk
    integer:: numx,numy,begx,begy,jind,ix,iy,ii
    character(len=8):: dimname1,dimname2

    ! check that the horizontal gridname is a valid grid
    !-----------------------------------------------------
    gridid = cam_grid_id(trim(hgridname))
    if(.not. cam_grid_check(gridid)) then
      write(iulog, *) 'set_model_Hgrid: Internal error, no ',trim(hgridname),' hgridname'
      call endrun('ERROR: No GRIDNAME')
    endif

    ! This is a CAM grid, get the details from the gridid value
    !----------------------------------------------------------
    call cam_grid_get_gcid        (gridid,gcid)
    call cam_grid_get_array_bounds(gridid,griddims)
    call cam_grid_get_dim_names   (gridid,dimname1,dimname2)
    Latvals => cam_grid_get_latvals(gridid)
    Lonvals => cam_grid_get_lonvals(gridid)
    PEsize   = cam_grid_get_local_size(gridid)
    Latsize = size(Latvals)
    Lonsize = size(Lonvals)

    ! CHECK PEsize == size(gcid)
    !---------------------------

    ! Allocate some space
    !----------------------
    allocate(G_Lon (PEsize))
    allocate(G_Lat (PEsize))
    allocate(G_Xmap(PEsize))
    allocate(G_Ymap(PEsize))

    ! Branch to grid-type cases
    !---------------------------
    if((PEsize.eq.Lonsize).and.(PEsize.eq.Latsize)) then
      ! Lat/Lon values are in an indexed map
      !--------------------------------------
      Hgrid_m%PackType = "UNSTRUCTURED"
      pcols    = 1 + griddims(1,2) - griddims(1,1)
      begchunk = griddims(2,1)
      endchunk = griddims(2,2)
      nchunks  = 1 + endchunk - begchunk

      ! create a packed set of active gridpoints,
      ! and a map to the values in the model array.
      !--------------------------------------------
      G_ncol = 0
      do ii = 1, PEsize
        if(gcid(ii).ge.1) then
          ichunk = (ii-1)/pcols
          G_ncol = G_ncol + 1
          G_Xmap(G_ncol) = ii - pcols*ichunk
          G_Ymap(G_ncol) = ichunk + begchunk
          G_Lon (G_ncol) = Lonvals(ii)
          G_Lat (G_ncol) = Latvals(ii)
        endif
      end do
    elseif(PEsize.eq.(Lonsize*Latsize)) then
      ! Lat/Lon values are in a rectilinear grid
      !------------------------------------------
      Hgrid_m%PackType = "LATLON"
      numx = 1 + griddims(1,2) - griddims(1,1)
      numy = 1 + griddims(2,2) - griddims(2,1)
      begx = griddims(1,1)
      begy = griddims(2,1)

      ! create a packed set of active gridpoints,
      ! and a map to the values in the model array.
      !--------------------------------------------
      G_ncol = 0
      do ii = 1, PEsize
        if(gcid(ii).ge.1) then
          jind   = (ii-1)/numx
          ix     = ii - jind*numx
          iy     = 1  + jind
          G_ncol = G_ncol + 1
          G_Xmap(G_ncol) = begx + ix - 1
          G_Ymap(G_ncol) = begy + iy - 1
          G_Lon (G_ncol) = Lonvals(G_Xmap(G_ncol))
          G_Lat (G_ncol) = Latvals(G_Ymap(G_ncol))
        endif
      end do
    else
      ! Unknown Case: Error stop!
      !----------------------------
      call endrun('ERROR: Mis-Match between PEsize and grid dimensions')
    endif

    ! Now Set the Hgrid values
    !-----------------------------
    Hgrid_m%name           = trim(hgridname)
    Hgrid_m%Type           = "PACKED_MODEL"
    Hgrid_m%model_gridname = trim(hgridname)
    Hgrid_m%model_gridid   = gridid
    Hgrid_m%ncol           = G_ncol
    if(allocated(Hgrid_m%Lon )) deallocate(Hgrid_m%Lon )
    if(allocated(Hgrid_m%Lat )) deallocate(Hgrid_m%Lat )
    if(allocated(Hgrid_m%Xmap)) deallocate(Hgrid_m%Xmap)
    if(allocated(Hgrid_m%Ymap)) deallocate(Hgrid_m%Ymap)
    allocate(Hgrid_m%Lon (Hgrid_m%ncol))
    allocate(Hgrid_m%Lat (Hgrid_m%ncol))
    allocate(Hgrid_m%Xmap(Hgrid_m%ncol))
    allocate(Hgrid_m%Ymap(Hgrid_m%ncol))
    Hgrid_m%Xmap(:) = G_Xmap(1:Hgrid_m%ncol)
    Hgrid_m%Ymap(:) = G_Ymap(1:Hgrid_m%ncol)
    Hgrid_m%Lon (:) = G_Lon (1:Hgrid_m%ncol)
    Hgrid_m%Lat (:) = G_Lat (1:Hgrid_m%ncol)

    Hgrid_m%Hdim1name = trim(dimname1)
    Hgrid_m%Hdim2name = trim(dimname2)
    Hgrid_m%Hdim1b    = griddims(1,1)
    Hgrid_m%Hdim1e    = griddims(1,2)
    Hgrid_m%Hdim2b    = griddims(2,1)
    Hgrid_m%Hdim2e    = griddims(2,2)
    Hgrid_m%Hdim1len  = 1 + griddims(1,2) - griddims(1,1)
    Hgrid_m%Hdim2len  = 1 + endchunk - begchunk

    ! Clean up
    !----------------
    nullify(gcid   )
    nullify(Lonvals)
    nullify(Latvals)
    deallocate(G_Lon )
    deallocate(G_Lat )
    deallocate(G_Xmap)
    deallocate(G_Ymap)

    !***************************
    ! (Temporary) Diag output
    !***************************
    if(.FALSE.) then
      write(iulog,*) 'PFCDIAG: Packed gridpoints for ',trim(Hgrid_m%name),' grid'
      write(iulog,*) 'PFCDIAG: Hgrid_m%type          =',trim(Hgrid_m%Type)
      write(iulog,*) 'PFCDIAG: Hgrid_m%PackType      =',trim(Hgrid_m%PackType)
      write(iulog,*) 'PFCDIAG: Hgrid_m%model_gridname=',trim(Hgrid_m%model_gridname)
      write(iulog,*) 'PFCDIAG: Hgrid_m%model_gridid  =',Hgrid_m%model_gridid
      write(iulog,*) 'PFCDIAG: Hgrid_m%ncol          =',Hgrid_m%ncol
      do ii=1,Hgrid_m%ncol
        write(iulog,*) 'PFCDIAG: ii=',ii,' Lon/Lat=',Hgrid_m%Lon (ii),Hgrid_m%Lat (ii), &
                                         ' X,Y Map=',Hgrid_m%Xmap(ii),Hgrid_m%Ymap(ii)
      end do
      write(iulog,*) 'PFCDIAG: '
    endif

    ! End Routine
    !------------------
    return
  end subroutine set_model_Hgrid
  !================================================================


  !================================================================
  subroutine set_model_Vgrid(Vgrid_m,vgridname)
    !
    ! set_model_Vgrid: For the given vgridname (either 'lev' or 'ilev'),
    !                  initialize the given model Vgrid datastructure
    !                  with model levels.
    !==============================================================
    use pmgrid,only: plev, plevp
    use hycoef,only: hyam,hybm,hyai,hybi,ps0
    !
    ! Passed variables
    !-------------------
    type(vgrid_t)   ,intent(inout):: Vgrid_m
    character(len=*),intent(in   ):: vgridname

    ! Initialize Vertical datastructure with model coordinates
    !-----------------------------------------------------------
    if(trim(vgridname).eq.'lev') then
      ! Mid points of Model Layers
      !----------------------------
      Vgrid_m%name      = "HYP_LEV"
      Vgrid_m%Type      = "HYBRIDP"
      Vgrid_m%Vdimname  = 'lev'
      Vgrid_m%nlev      = plev
      Vgrid_m%ngrid_var = 3
      Vgrid_m%nreq_var  = 1
      if(allocated(Vgrid_m%grid_var)) deallocate(Vgrid_m%grid_var)
      if(allocated(Vgrid_m%req_var )) deallocate(Vgrid_m%req_var)
      if(allocated(Vgrid_m%HYA     )) deallocate(Vgrid_m%HYA)
      if(allocated(Vgrid_m%HYB     )) deallocate(Vgrid_m%HYB)
      allocate(Vgrid_m%grid_var(Vgrid_m%ngrid_var))
      allocate(Vgrid_m%req_var (Vgrid_m%nreq_var ))
      allocate(Vgrid_m%HYA     (Vgrid_m%nlev))
      allocate(Vgrid_m%HYB     (Vgrid_m%nlev))
      Vgrid_m%grid_var(1) = "hyam"
      Vgrid_m%grid_var(2) = "hybm"
      Vgrid_m%grid_var(3) = "Pref"
      Vgrid_m%req_var (1) = "PS"
      Vgrid_m%HYA(1:plev) = hyam(1:plev)
      Vgrid_m%HYB(1:plev) = hybm(1:plev)
      Vgrid_m%Pref        = ps0
    elseif(trim(vgridname).eq.'lev-index') then
      ! Mid points of Model Layers
      !----------------------------
      Vgrid_m%name      = "lev"
      Vgrid_m%Type      = "INDEX"
      Vgrid_m%Vdimname  = 'lev'
      Vgrid_m%nlev      = plev
      Vgrid_m%ngrid_var = 1
      Vgrid_m%nreq_var  = 0
      if(allocated(Vgrid_m%grid_var)) deallocate(Vgrid_m%grid_var)
      if(allocated(Vgrid_m%req_var )) deallocate(Vgrid_m%req_var)
      if(allocated(Vgrid_m%HYA     )) deallocate(Vgrid_m%HYA)
      if(allocated(Vgrid_m%HYB     )) deallocate(Vgrid_m%HYB)
      allocate(Vgrid_m%grid_var(Vgrid_m%ngrid_var))
      Vgrid_m%grid_var(1) = "ANYINDEX"
    elseif(trim(vgridname).eq.'ilev') then
      ! Interface points of Model Layers
      !----------------------------------
      Vgrid_m%name      = "HYP_ILEV"
      Vgrid_m%Type      = "HYBRIDP"
      Vgrid_m%Vdimname  = 'ilev'
      Vgrid_m%nlev      = plevp
      Vgrid_m%ngrid_var = 3
      Vgrid_m%nreq_var  = 1
      if(allocated(Vgrid_m%grid_var)) deallocate(Vgrid_m%grid_var)
      if(allocated(Vgrid_m%req_var )) deallocate(Vgrid_m%req_var)
      if(allocated(Vgrid_m%HYA     )) deallocate(Vgrid_m%HYA)
      if(allocated(Vgrid_m%HYB     )) deallocate(Vgrid_m%HYB)
      allocate(Vgrid_m%grid_var(Vgrid_m%ngrid_var))
      allocate(Vgrid_m%req_var (Vgrid_m%nreq_var ))
      allocate(Vgrid_m%HYA     (Vgrid_m%nlev))
      allocate(Vgrid_m%HYB     (Vgrid_m%nlev))
      Vgrid_m%grid_var(1)  = "hyai"
      Vgrid_m%grid_var(2)  = "hybi"
      Vgrid_m%grid_var(3)  = "Pref"
      Vgrid_m%req_var (1)  = "PS"
      Vgrid_m%HYA(1:plevp) = hyai(1:plevp)
      Vgrid_m%HYB(1:plevp) = hybi(1:plevp)
      Vgrid_m%Pref         = ps0
    elseif(trim(vgridname).eq.'ilev-index') then
      ! Interface points of Model Layers
      !----------------------------------
      Vgrid_m%name      = "ilev"
      Vgrid_m%Type      = "INDEX"
      Vgrid_m%Vdimname  = 'ilev'
      Vgrid_m%nlev      = plevp
      Vgrid_m%ngrid_var = 1
      Vgrid_m%nreq_var  = 0
      if(allocated(Vgrid_m%grid_var)) deallocate(Vgrid_m%grid_var)
      if(allocated(Vgrid_m%req_var )) deallocate(Vgrid_m%req_var)
      if(allocated(Vgrid_m%HYA     )) deallocate(Vgrid_m%HYA)
      if(allocated(Vgrid_m%HYB     )) deallocate(Vgrid_m%HYB)
      allocate(Vgrid_m%grid_var(Vgrid_m%ngrid_var))
      Vgrid_m%grid_var(1)  = "ANYINDEX"
    elseif(trim(vgridname).eq.'HORIZ_ONLY') then
      ! Entry for HORIZ_ONLY (2D) arrays
      !----------------------------------
      Vgrid_m%name      = "HORIZ_ONLY"
      Vgrid_m%Type      = "HORIZ_ONLY"
      Vgrid_m%Vdimname  = 'HORIZ_ONLY'
      Vgrid_m%nlev      = 1
      Vgrid_m%ngrid_var = 0
      Vgrid_m%nreq_var  = 0
      if(allocated(Vgrid_m%grid_var)) deallocate(Vgrid_m%grid_var)
      if(allocated(Vgrid_m%req_var )) deallocate(Vgrid_m%req_var)
      if(allocated(Vgrid_m%HYA     )) deallocate(Vgrid_m%HYA)
      if(allocated(Vgrid_m%HYB     )) deallocate(Vgrid_m%HYB)
    else
      ! Unknown Case: Error stop!
      !----------------------------
      call endrun('set_model_Vgrid: ERROR Unknown vgridname = '//trim(vgridname))
    endif

    ! End Routine
    !------------------
    return
  end subroutine set_model_Vgrid
  !================================================================


  !================================================================
  subroutine infld_dataset_3d(VarName,DataSetID,PathId,dimname1,dimname2,dimname3,          &
                                                       dim1b,dim1e,dim2b,dim2e,dim3b,dim3e, &
                                                       field,readvar,timelevel,fillvalue)
    !
    ! infld_dataset_3d: Read a 3D field from dataset into a 3D variable
    !==============================================================
    implicit none
    !
    ! Passed variables
    !------------------
    character(len=*) ,intent(in ):: VarName
    integer          ,intent(in ):: DataSetID
    integer          ,intent(in ):: PathID
    character(len=*) ,intent(in ):: dimname1
    character(len=*) ,intent(in ):: dimname2
    character(len=*) ,intent(in ):: dimname3
    integer          ,intent(in ):: dim1b
    integer          ,intent(in ):: dim1e
    integer          ,intent(in ):: dim2b
    integer          ,intent(in ):: dim2e
    integer          ,intent(in ):: dim3b
    integer          ,intent(in ):: dim3e
    real(r8),target  ,intent(out):: field(dim1b:dim1e,dim2b:dim2e,dim3b:dim3e)
    logical          ,intent(out):: readvar
    integer ,optional,intent(in ):: timelevel
    real(r8),optional,intent(out):: fillvalue
    !
    ! Local Values
    !--------------
    integer :: time
    real(r8):: DS_fillvalue

    if(present(timelevel)) then
      call apply_dataPath(DataSetID,PathID,VarName,timelevel,DS_fillvalue, &
                                  dim1b,dim1e,dim2b,dim2e,dim3b,dim3e,field)
    else
      time = 1
      call apply_dataPath(DataSetID,PathID,VarName,time,DS_fillvalue, &
                             dim1b,dim1e,dim2b,dim2e,dim3b,dim3e,field)
    endif
    if(present(fillvalue)) fillvalue = DS_fillvalue

    readvar = .true.

    ! End Routine
    !------------------
    return
  end subroutine infld_dataset_3d
  !================================================================


  !================================================================
  subroutine infld_dataset_horiz(VarName,DataSetID,PathId,dimname1,dimname2,               &
                                                          dim1b,dim1e,dim2b,dim2e,         &
                                                          field,readvar,timelevel,fillvalue)
    !
    ! infld_dataset_horiz: Read a HORIZ_ONLY field from dataset into a 2D variable
    !==============================================================
    implicit none
    !
    ! Passed variables
    !------------------
    character(len=*) ,intent(in ):: VarName
    integer          ,intent(in ):: DataSetID
    integer          ,intent(in ):: PathID
    character(len=*) ,intent(in ):: dimname1
    character(len=*) ,intent(in ):: dimname2
    integer          ,intent(in ):: dim1b
    integer          ,intent(in ):: dim1e
    integer          ,intent(in ):: dim2b
    integer          ,intent(in ):: dim2e
    real(r8),target  ,intent(out):: field(dim1b:dim1e,dim2b:dim2e)
    logical          ,intent(out):: readvar
    integer ,optional,intent(in ):: timelevel
    real(r8),optional,intent(out):: fillvalue
    !
    ! Local Values
    !--------------
    integer :: time
    real(r8):: DS_fillvalue

    if(present(timelevel)) then
      call apply_dataPath(DataSetID,PathID,VarName,timelevel,DS_fillvalue, &
                                              dim1b,dim1e,dim2b,dim2e,field)
    else
      time = 1
      call apply_dataPath(DataSetID,PathID,VarName,time,DS_fillvalue, &
                                         dim1b,dim1e,dim2b,dim2e,field)
    endif
    if(present(fillvalue)) fillvalue = DS_fillvalue

    readvar = .true.

    ! End Routine
    !------------------
    return
  end subroutine infld_dataset_horiz
  !================================================================

end module ncdio_atm
