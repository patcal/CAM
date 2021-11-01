module PathDesc_type
!===========================================================================
!
! MODULE: PathDesc_type
!
! DESCRIPTION: A DataPath type created to avoid a circular dependency
!
!===========================================================================
  ! Useful Modules
  !---------------------
  use pio,only: file_desc_t, var_desc_t, io_desc_t
  use shr_kind_mod,  only: r8 => shr_kind_r8

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private
  save

  public :: PathDesc_t

  type PathDesc_t
    type(file_desc_t)       :: ncid
    integer                 :: ncid_nf90  ! PFC DIAG Temp to test PIO ERROR with nf90
    character(len=80)       :: varname    ! PFC DIAG Temp to test PIO ERROR with nf90
    integer                 :: vlen       ! PFC DIAG Temp to test PIO ERROR with nf90
    integer                 :: timelevel
    integer                 :: ncol
    integer                 :: nlev_d
    integer                 :: nlev_m
    type(var_desc_t)        :: VarDesc
    type(io_desc_t) ,pointer:: iodesc
    real(r8)                :: fillvalue
  end type PathDesc_t

end module PathDesc_type
