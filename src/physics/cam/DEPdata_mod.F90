module DEPdata_mod
!======================================================================
!
! Purpose: Module to maintain sets of input/output data values that 
!          vertical mapping depend on (e.g. PS values). There is not a 
!          gaurantee that the needed values are contained in the same file
!          as the value being interpolated, so this module provides a 
!          common location for these values that all registered mappings
!          access. 
!
!          The user is responsible for updating these values prior to 
!          invoking the mappings that depend on them.
!
!  ** INCOMPLETE** need to finish development to provide PS, PHIS, etc, 
!                  values to vertical interp options that require them.
!
!======================================================================
  ! Usefule modules
  !-----------------
  use shr_kind_mod,   only: r8=>SHR_KIND_R8, cs=>SHR_KIND_CS, cl=>SHR_KIND_CL
  use cam_abortutils, only: endrun



  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public :: init_DEP_Data
!  public :: set_DEP_input
!  public :: set_DEP_output
!  public :: get_DEP_invar
!  public :: get_DEP_outvar


  type DEPvar_t
    character(len=8 )   :: Name
    character(len=32)   :: Shape
    logical             :: CONST
    integer             :: ncol
    integer             :: nlev
    real(r8),allocatable:: val(:,:)
  end type DEPvar_t

  type DEP_Data_t
      private
      integer                   :: num_input
      integer                   :: num_output
      type(DEPvar_t),allocatable:: input (:)
      type(DEPvar_t),allocatable:: output(:)
    contains
      procedure,pass:: init      =>  init_DEP_Data
  end type DEP_Data_t

contains
  !==================================================================
  subroutine init_DEP_Data(this,nlev_in,nlev_out,ncol,                   &
                           ninput ,invar_name ,invar_shape ,invar_const, &
                           noutput,outvar_name,outvar_shape,outvar_const )
    !
    !================================================================
    ! 
    ! Passed variables 
    !------------------
    class(DEP_Data_t) :: this
    integer,intent(in ):: nlev_in,nlev_out,ncol,ninput,noutput
    character(len=8 ),intent(in)::   invar_name(:)
    character(len=8 ),intent(in)::  outvar_name(:)
    character(len=32),intent(in)::  invar_shape(:)
    character(len=32),intent(in):: outvar_shape(:)
    logical          ,intent(in)::  invar_const(:)
    logical          ,intent(in):: outvar_const(:)
    !
    ! Local values
    !------------
    integer:: nn

    ! Allocate storage space for values needed to process input data
    !---------------------------------------------------------------
    if(ninput.gt.0) then
      this%num_input = ninput
      if(allocated(this%input)) deallocate(this%input)
      allocate(this%input(this%num_input))
      do nn=1,this%num_input
        this%input(nn)%Name  = trim(invar_name (nn))
        this%input(nn)%Shape = trim(invar_shape(nn))
        this%input(nn)%CONST =      invar_const(nn)
        if(trim(this%input(nn)%Shape).eq."HORIZ_ONLY") then
          this%input(nn)%ncol = ncol
          this%input(nn)%nlev = 1
        elseif(trim(this%input(nn)%Shape).eq."3D") then
          this%input(nn)%ncol = ncol
          this%input(nn)%nlev = nlev_in
        elseif(trim(this%input(nn)%Shape).eq."VERT_ONLY") then
          this%input(nn)%ncol = 1
          this%input(nn)%nlev = nlev_in
        else
          call endrun()
        endif
        allocate(this%input(nn)%val(this%input(nn)%ncol,this%input(nn)%nlev))
        this%input(nn)%val(:,:) = 0._r8
      end do
    else
      this%num_input = 0
    endif

    ! Allocate storage space for values needed to process output data
    !---------------------------------------------------------------
    if(noutput.gt.0) then
      this%num_output = noutput
      if(allocated(this%output)) deallocate(this%output)
      allocate(this%output(this%num_output))
      do nn=1,this%num_output
        this%output(nn)%Name  = trim(outvar_name (nn))
        this%output(nn)%Shape = trim(outvar_shape(nn))
        this%output(nn)%CONST =      outvar_const(nn)
        if(trim(this%output(nn)%Shape).eq."HORIZ_ONLY") then
          this%output(nn)%ncol = ncol
          this%output(nn)%nlev = 1
        elseif(trim(this%output(nn)%Shape).eq."3D") then
          this%output(nn)%ncol = ncol
          this%output(nn)%nlev = nlev_in
        elseif(trim(this%output(nn)%Shape).eq."VERT_ONLY") then
          this%output(nn)%ncol = 1
          this%output(nn)%nlev = nlev_in
        else
          call endrun()
        endif
        allocate(this%output(nn)%val(this%output(nn)%ncol,this%output(nn)%nlev))
        this%output(nn)%val(:,:) = 0._r8
      end do
    else
      this%num_output = 0
    endif

          
  


    ! End Routine
    !--------------
    return
  end subroutine init_DEP_Data
  !==================================================================




end module DEPdata_mod
