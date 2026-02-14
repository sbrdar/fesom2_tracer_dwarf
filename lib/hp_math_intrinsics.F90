!> @file hp_math_intrinsics.F90
!! @brief Half-precision wrappers for intrinsic math functions
!! @details nvfortran supports real(kind=2) for storage but most intrinsic
!!          math functions lack FP16 overloads. This module provides elemental
!!          wrappers that promote to real(4), call the intrinsic, and demote
!!          back. The generic interfaces use the SAME names as the intrinsics,
!!          so existing code just needs 'use hp_math_intrinsics' to get FP16
!!          support transparently.
!!
!!          When USE_HALF_PRECISION is not defined, the module is empty (no-op).

#ifdef USE_HALF_PRECISION

module hp_math_intrinsics
  implicit none
  private

  ! Extend intrinsic generics with real(2) overloads
  interface cos
    module procedure hp_cos
  end interface
  interface sin
    module procedure hp_sin
  end interface
  interface tan
    module procedure hp_tan
  end interface
  interface sqrt
    module procedure hp_sqrt
  end interface
  interface abs
    module procedure hp_abs
  end interface
  interface asin
    module procedure hp_asin
  end interface
  interface atan2
    module procedure hp_atan2
  end interface
  interface sign
    module procedure hp_sign
  end interface
  interface exp
    module procedure hp_exp
  end interface
  interface log
    module procedure hp_log
  end interface

  public :: cos, sin, tan, sqrt, abs, asin, atan2, sign, exp, log

contains

  elemental function hp_cos(x) result(r)
    real(2), intent(in) :: x
    real(2) :: r
    intrinsic :: cos
    r = real(cos(real(x, 4)), 2)
  end function

  elemental function hp_sin(x) result(r)
    real(2), intent(in) :: x
    real(2) :: r
    intrinsic :: sin
    r = real(sin(real(x, 4)), 2)
  end function

  elemental function hp_tan(x) result(r)
    real(2), intent(in) :: x
    real(2) :: r
    intrinsic :: tan
    r = real(tan(real(x, 4)), 2)
  end function

  elemental function hp_sqrt(x) result(r)
    real(2), intent(in) :: x
    real(2) :: r
    intrinsic :: sqrt
    r = real(sqrt(real(x, 4)), 2)
  end function

  elemental function hp_abs(x) result(r)
    real(2), intent(in) :: x
    real(2) :: r
    intrinsic :: abs
    r = real(abs(real(x, 4)), 2)
  end function

  elemental function hp_asin(x) result(r)
    real(2), intent(in) :: x
    real(2) :: r
    intrinsic :: asin
    r = real(asin(real(x, 4)), 2)
  end function

  elemental function hp_atan2(y, x) result(r)
    real(2), intent(in) :: y, x
    real(2) :: r
    intrinsic :: atan2
    r = real(atan2(real(y, 4), real(x, 4)), 2)
  end function

  elemental function hp_sign(a, b) result(r)
    real(2), intent(in) :: a, b
    real(2) :: r
    intrinsic :: sign
    r = real(sign(real(a, 4), real(b, 4)), 2)
  end function

  elemental function hp_exp(x) result(r)
    real(2), intent(in) :: x
    real(2) :: r
    intrinsic :: exp
    r = real(exp(real(x, 4)), 2)
  end function

  elemental function hp_log(x) result(r)
    real(2), intent(in) :: x
    real(2) :: r
    intrinsic :: log
    r = real(log(real(x, 4)), 2)
  end function

end module hp_math_intrinsics

#else

! Empty module when half precision is not enabled
module hp_math_intrinsics
  implicit none
end module hp_math_intrinsics

#endif
