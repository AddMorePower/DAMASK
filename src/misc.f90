!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @author Philip Eisenlohr, Michigan State University
!> @brief Miscellaneous tools.
!--------------------------------------------------------------------------------------------------
module misc
  use prec
  use constants

  implicit none(type,external)
  private

  interface misc_optional
    module procedure misc_optional_bool
    module procedure misc_optional_pI32
    module procedure misc_optional_pI64
    module procedure misc_optional_real
    module procedure misc_optional_complex
    module procedure misc_optional_str
  end interface misc_optional

  public :: &
    misc_init, &
    misc_selfTest, &
    misc_optional, &
    misc_prefixOptions, &
    misc_ones, &
    misc_zeros

contains


!--------------------------------------------------------------------------------------------------
!> @brief Do self test.
!--------------------------------------------------------------------------------------------------
subroutine misc_init()

  print'(/,1x,a)', '<<<+-  misc init  -+>>>'

  call misc_selfTest()

end subroutine misc_init


!--------------------------------------------------------------------------------------------------
!> @brief Return bool value if given, otherwise default.
!--------------------------------------------------------------------------------------------------
pure elemental function misc_optional_bool(given,default) result(var)

  logical, intent(in), optional :: given
  logical, intent(in)           :: default
  logical                       :: var


  if (present(given)) then
    var = given
  else
    var = default
  end if

end function misc_optional_bool


!--------------------------------------------------------------------------------------------------
!> @brief Return integer(pI32) value if given, otherwise default.
!--------------------------------------------------------------------------------------------------
pure elemental function misc_optional_pI32(given,default) result(var)

  integer(pI32), intent(in), optional :: given
  integer(pI32), intent(in)           :: default
  integer(pI32)                       :: var


  if (present(given)) then
    var = given
  else
    var = default
  end if

end function misc_optional_pI32


!--------------------------------------------------------------------------------------------------
!> @brief Return integer(pI64) value if given, otherwise default.
!--------------------------------------------------------------------------------------------------
pure elemental function misc_optional_pI64(given,default) result(var)

  integer(pI64), intent(in), optional :: given
  integer(pI64), intent(in)           :: default
  integer(pI64)                       :: var


  if (present(given)) then
    var = given
  else
    var = default
  end if

end function misc_optional_pI64


!--------------------------------------------------------------------------------------------------
!> @brief Return real value if given, otherwise default.
!--------------------------------------------------------------------------------------------------
pure elemental function misc_optional_real(given,default) result(var)

  real(pREAL), intent(in), optional :: given
  real(pREAL), intent(in)           :: default
  real(pREAL)                       :: var


  if (present(given)) then
    var = given
  else
    var = default
  end if

end function misc_optional_real


!--------------------------------------------------------------------------------------------------
!> @brief Return complex value if given, otherwise default.
!--------------------------------------------------------------------------------------------------
pure elemental function misc_optional_complex(given,default) result(var)

  complex(pREAL), intent(in), optional :: given
  complex(pREAL), intent(in)           :: default
  complex(pREAL)                       :: var


  if (present(given)) then
    var = given
  else
    var = default
  end if

end function misc_optional_complex

!--------------------------------------------------------------------------------------------------
!> @brief Return string value if given, otherwise default.
!--------------------------------------------------------------------------------------------------
pure function misc_optional_str(given,default) result(var)

  character(len=*), intent(in), optional :: given
  character(len=*), intent(in)           :: default
  character(len=:), allocatable          :: var


  if (present(given)) then
    var = given
  else
    var = default
  end if

end function misc_optional_str

!--------------------------------------------------------------------------------------------------
!> @brief Add prefix to options in string.
!> @detail An option starts with a dash followed by at least one letter.
!--------------------------------------------------------------------------------------------------
pure function misc_prefixOptions(string,prefix) result(prefixed)

  character(len=*), intent(in)  :: string,prefix
  character(len=:), allocatable :: prefixed

  integer :: i,N


  prefixed = ''
  N = len(string)
  do i = 1, N
    prefixed = prefixed//string(i:i)
    if (string(i:i) == '-' .and. verify(string(min(i+1,N):min(i+1,N)),LOWER//UPPER) == 0) &
      prefixed = prefixed//prefix
  end do

end function misc_prefixOptions


!--------------------------------------------------------------------------------------------------
!> @brief 1D array of zeros.
!--------------------------------------------------------------------------------------------------
pure function misc_zeros(N)

  integer, intent(in) :: N                                                                          !< number of zeros
  real(pREAL), dimension(N) :: misc_zeros


  misc_zeros = 0._pREAL

end function misc_zeros


!--------------------------------------------------------------------------------------------------
!> @brief 1D array of ones.
!--------------------------------------------------------------------------------------------------
pure function misc_ones(N)

  integer, intent(in) :: N                                                                          !< number of ones
  real(pREAL), dimension(N) :: misc_ones


  misc_ones = 1._pREAL

end function misc_ones


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of some misc functions.
!--------------------------------------------------------------------------------------------------
subroutine misc_selfTest()

  real(pREAL) :: r
  real(pREAL), dimension(:), allocatable :: rN
  character(len=:),      allocatable :: str,out
  integer :: N


  call random_number(r)
  N = int(r*99._pREAL)
  allocate(rN(N),source=0.0_pREAL)
  call random_number(rN)

  if (test_str('DAMASK') /= 'DAMASK')                        error stop 'optional_str, present'
  if (test_str() /= 'default')                               error stop 'optional_str, not present'
  if (misc_optional(default='default') /= 'default')         error stop 'optional_str, default only'
  if (test_pI32(20191102_pI32) /= 20191102_pI32)             error stop 'optional_int, present'
  if (test_pI32() /= 42_pI32)                                error stop 'optional_int, not present'
  if (misc_optional(default=20191102) /= 20191102)           error stop 'optional_int, default only'
  if (dNeq(test_real(r),r))                                  error stop 'optional_real, present'
  if (dNeq(test_real(),0.0_pREAL))                           error stop 'optional_real, not present'
  if (dNeq(misc_optional(default=r),r))                      error stop 'optional_real, default only'
  if (any(dNeq(misc_optional(default=rN),rN)))               error stop 'optional_real array, default only'
  if (test_bool(r<0.5_pREAL) .neqv. r<0.5_pREAL)             error stop 'optional_bool, present'
  if (.not. test_bool())                                     error stop 'optional_bool, not present'
  if (misc_optional(default=r>0.5_pREAL) .neqv. r>0.5_pREAL) error stop 'optional_bool, default only'

  str='-a -1 -more 123 -flag -'
  out=misc_prefixOptions(str,'p_')
  if (out /= '-p_a -1 -p_more 123 -p_flag -')                error stop 'misc_prefixOptions'

  N = int(r*99._pREAL)
  if (size(misc_zeros(N)) /= N)           error stop 'shape zeros'
  if (size(misc_ones(N)) /= N)            error stop 'shape ones'
  if (any(dNeq(misc_zeros(N),0.0_pREAL))) error stop 'value zeros'
  if (any(dNeq(misc_ones(N),1.0_pREAL)))  error stop 'value ones'

contains

  pure function test_str(str_in) result(str_out)
#ifndef __GFORTRAN__
    import, only: misc_optional_str
#endif
    character(len=:), allocatable          :: str_out
    character(len=*), intent(in), optional :: str_in


    str_out = misc_optional_str(str_in,'default')

  end function test_str


  pure function test_pI32(int_in) result(int_out)
#ifndef __GFORTRAN__
    import, only: pI32, misc_optional_pI32
#endif
    integer(pI32)                       :: int_out
    integer(pI32), intent(in), optional :: int_in


    int_out = misc_optional_pI32(int_in,42_pI32)

  end function test_pI32


  pure function test_real(real_in) result(real_out)
#ifndef __GFORTRAN__
    import, only: pREAL, misc_optional_real
#endif
    real(pREAL)                       :: real_out
    real(pREAL), intent(in), optional :: real_in


    real_out = misc_optional_real(real_in,0.0_pREAL)

  end function test_real


  pure function test_bool(bool_in) result(bool_out)
#ifndef __GFORTRAN__
    import, only: misc_optional_bool
#endif
    logical                       :: bool_out
    logical, intent(in), optional :: bool_in


    bool_out = misc_optional_bool(bool_in,.true.)

  end function test_bool


end subroutine misc_selfTest

end module misc
