!********************************************************
!*             Equations of degree 2 to 4               *
!* ---------------------------------------------------- *
!* This program calculates the real*8 or complex roots of *
!* algebraic equations of degree 2, 3 and 4.            *
!* ---------------------------------------------------- *
!* Reference: Mathématiques et statitiques - Programmes *
!*            en BASIC. Editions du P.S.I., 1981        *
!*            [BIBLI 13].                               *
!*                                                      *
!*                   F90 version by J-P Moreau, Paris   *
!*                          (www.jpmoreau.fr)           *
!* ---------------------------------------------------- *
module m_poly_roots
  
  use, intrinsic :: iso_fortran_env, only: real64, int64
  implicit none

  real(real64)            :: d_qnan = transfer(9221120237041090560_int64, 1._real64)

  integer, parameter      :: ONLY_REAL_ROOT = 0, ALL_ROOTS = 1
contains

  subroutine polynomial_roots_deg4(coeffs, roots_real, roots_imag)

    implicit none

    real*8, intent(in)        :: coeffs(5)  ! c(1) * x^4 + c(2) * x^3 + ... + c(5)
    real*8, intent(out)       :: roots_real(4), roots_imag(4)

    ! Local variables
    real*8                    :: a(1:4, 0:4), im, ii
    integer                   :: degree, i

    ! Degree must be >= 2 and <= 4
    degree = 4
    do i=1,3
      if (coeffs(i) == 0.0) then
        degree = degree - 1
      else
        exit
      endif
    enddo

    if (degree == 1) then
      roots_real(1) = -coeffs(5) / coeffs(4)
      roots_imag(1) = 0.0
    else
      a = 0.
      do i=0,degree
        a(degree, i) = coeffs(5 - i)
      end do

      ! Calling root search main subroutine
      call Root_4(degree, a, roots_real, im, ii)

      select case(degree)
      case(2)
        roots_imag(1) = im
        roots_imag(2) = -im
      case(3)
        roots_imag(1) = im
        roots_imag(2) = -im
        roots_imag(3) = 0.0
      case(4)
        roots_imag(1) = im
        roots_imag(2) = -im
        roots_imag(3) = ii
        roots_imag(4) = -ii
      end select
    endif

    roots_real(degree+1:4) = d_qnan
    roots_imag(degree+1:4) = d_qnan

  end subroutine


  !*******************************************
  !* This Subroutine calculates the roots of *
  !*    X^2 + A(2,1)*X + A(2,0) = 0          *
  !*    W = Determinant of equation          *
  !*******************************************
  subroutine Root_2(A, r, im)
    implicit none

    real*8, intent(in)      :: A(1:4,0:4)
    real*8, intent(out)   :: r(1:4), im

    call polynomial_roots_deg2([1.0d0, A(2,1), A(2,0)], r(1:2), im)
    if (r(1) < r(2)) r(1:2) = [r(2), r(1)]
  end subroutine


  !********************************************
  !* This subroutine calculates the roots of  *
  !* X^3 + A(3,2)*X^2 + A(3,1)*X + A(3,0) = 0 *
  !********************************************
  subroutine Root_3(A, sw, r, im)
    use m_oned_rootfinding,   only: brent, bisection
    implicit none

    real*8, intent(inout)   :: A(1:4,0:4), r(1:4), im
    integer, intent(in)   :: sw

    ! Local variables
    real*8                  :: aa, am, b, tt, xa, xb, xc
    integer               :: i

    real*8, parameter       :: ERR_TOL = 0.
    integer, parameter    :: MAX_IT = 100

    ! one root equals zero
    if (A(3,0) == 0.0) then
      xc = 0.0
    else
      ! looking for( maximum coefficient in fabsolute magnitude
      am = abs(A(3,0))
      do i=1, 3
        tt = abs(A(3,i))
        if (am < tt) am = tt
      end do
      !Define interval where a real root exists
      ! according to sign of A(3,0)
      if (A(3,0) > 0) then
        aa = -am - 1.0
        b = 0.0
      else
        aa = 0; b = am + 1.0
      end if

      ! Searching for( xc = real root in interval (a,b)

      ! Define segment (xa,xb) including the root
      xa = aa; xb = b

      xc = brent(eval_cubic, xa, xb, ERR_TOL, MAX_IT)

    endif

    ! r(3) is a real root
    r(3) = xc

    if (sw == ALL_ROOTS) then
      ! Calculates the roots of remaining quadratic equation
      ! Define the equation coefficients
      A(2,1) = A(3,2) + xc
      A(2,0) = A(3,1) + A(3,2) * xc + xc * xc
      call Root_2(A,r,im)
    endif

  contains

    pure real*8 function eval_cubic(x) result(Q)
      implicit none

      real*8, intent(in)    :: x

      ! Local variables
      real*8                :: b1, c2

      b1 = x + A(3,2)
      c2 = b1 * x + A(3,1)
      Q = c2 * x + A(3,0)
    end function


  end subroutine

  ! Root search main subroutine
  Subroutine Root_4(degree, A, r, im, ii)
    implicit none

    real*8, intent(inout)   :: A(1:4,0:4), r(1:4), im, ii

    ! Local variables
    integer               :: i, degree
    real*8                  :: aa, b, c, d, k, l, m, q, rr, s

    ! Normalization of coefficients
    do i=0, degree-1
      A(degree,i) = A(degree,i) / A(degree,degree)
    end do

    select case(degree)
    case(2)
      call Root_2(A,r,im)
    case(3)
      call Root_3(A,ALL_ROOTS,r,im)
    case(4)
      aa = A(4,3)
      b = A(4,2)
      c = A(4,1)
      d = A(4,0)
      q = b - (3.0 * aa * aa / 8.0)
      rr = c - (aa * b / 2.0) + (aa * aa * aa / 8.0)
      s = d - (aa * c / 4.0) + (aa * aa * b / 16.0) - (3.0 * aa * aa * aa * aa / 256.0)

      ! Define coefficients of cubic equation
      A(3,2) = q / 2.0
      A(3,1) = (q * q - 4.0 * s) / 16.0
      A(3,0) = -(rr * rr / 64.0)

      ! Calculate a real root of this equation
      if (rr /= 0.0 .or. A(3,1) >= 0.0) then
        ! Calling Root_3() with 2nd parameter sw = -1 to calculate one root only
        ! real root of above equation
        call Root_3(A,ONLY_REAL_ROOT,r,im)
      else
        ! Particular case when this equation is of 2nd order
        A(2,1) = A(3,2); A(2,0) = A(3,1)
        call Root_2(A,r,im)
        ! One takes the positive root
        r(3) = r(1)
      endif

      k = sqrt(r(3))

      ! Calculate L and M if k=0
      if (k == 0.0) then
        rr = sqrt(q * q - 4.0 * s)
      else
        q = q + (4.0 * r(3))
        rr = rr / (2.0 * k)
      end if
!
      l = (q - rr) / 2.0; m = (q + rr) / 2.0

      ! Solving two equations of degree 2
      A(2,1) = 2.0 * k; A(2,0) = l

      ! 1st equation
      call Root_2(A,r,ii)
      ! Transferring solutions in r[3], r[4], ii
      r(3) = r(1) - (A(4,3) / 4.0)
      r(4) = r(2) - (A(4,3) / 4.0)

      ! 2nd equation
      A(2,1) = -2.0 * k; A(2,0) = m
      call Root_2(A,r,im)
      r(2) = r(2) - (A(4,3) / 4.0)
      r(1) = r(1) - (A(4,3) / 4.0)
    end select
  end subroutine

  ! p(x) = c(1) * x^2 + c(2) * x + c(3)
  pure subroutine polynomial_roots_deg2(coeffs, roots, imag)
    implicit none

    real*8, intent(in)      :: coeffs(3)
    real*8, intent(out)     :: roots(2), imag

    ! Local variables
    real*8                  :: D, scaled(2), maxSqrt

    imag = 0.
    if (coeffs(1) == 0) then
      if (coeffs(2) == 0) then
        roots = d_qnan
      else
        roots(1) = -coeffs(3) / coeffs(2)
        roots(2) = roots(1)
      endif
    elseif (coeffs(3) == 0) then
      roots(1) = -coeffs(2) / coeffs(1)
      roots(2) = 0.0
    else
      ! Citarduaq's formula is used to allow for small c(1) in a numerically stable way
      scaled = coeffs(2:3) / coeffs(1)

      maxSqrt = sqrt(huge(scaled(1)))
      if (scaled(1) > maxSqrt .or. scaled(1) < -maxSqrt) then
        ! scaled(1)^2 would overflow, so we let √D ≈ |scaled(1)| instead
        roots(1) = -scaled(2) / scaled(1)
        roots(2) = scaled(2) / roots(1)
      else
        D = scaled(1)**2.0 - 4 * scaled(2)

        if (D <= 0) then
          if (D < 0) imag = sqrt(-D) / 2
          roots = -scaled(1) / 2
        else
          roots(1) = -2 * scaled(2) / (scaled(1) + sign(sqrt(D), scaled(1)))
          roots(2) = scaled(2) / roots(1)
        endif
      endif

    endif
    if (.not. isnan(roots(1)) .and. roots(1) > roots(2)) roots = [roots(2), roots(1)]
  end subroutine

  pure subroutine real_roots_2(roots, coeffs, residual)
    implicit none

    real*8, intent(in)      :: coeffs(3)
    real*8, intent(out)     :: roots(2)
    real*8, intent(out), optional :: residual(2)

    ! Local variables
    real*8                  :: imag

    call polynomial_roots_deg2(coeffs, roots, imag)
    if (imag /= 0.) roots = d_qnan
    if (present(residual)) then
      residual = d_qnan
      where (.not. isnan(roots)) residual = coeffs(1) * roots**2.0 + coeffs(2) * roots + coeffs(3)
    endif

  end

  ! QBC algorithm as described in
  !  To Solve a Real Cubic Equation (Lecture Notes for a Numerical Analysis Course)
  !    W. Kahan
  !    Mathematics Dep’t
  !    University of California
  !    Berkeley CA 94720
  !    Nov. 10, 1986
  pure subroutine polynomial_roots_deg3(coeffs, roots_real, roots_imag, only_real)
    implicit none

    real*8, intent(in)    :: coeffs(4)
    real*8, intent(out)   :: roots_real(3), roots_imag(3)
    logical, intent(in), optional :: only_real

    ! Local variables
    real*8                :: x, x1, x2, y

    call polynomial_roots_deg3_sub(coeffs(1), coeffs(2), coeffs(3), coeffs(4), x, x1, x2, y, &
      merge(only_real, .false., present(only_real)))

    roots_real(1) = x
    roots_imag(1) = 0.
    roots_real(2) = x1
    roots_imag(2) = y
    roots_real(3) = x2
    roots_imag(3) = -y
  end

  pure subroutine polynomial_roots_deg3_sub(A, B, C, D, x, x1, x2, y, only_real)
    implicit none

    real*8, intent(in)    :: A, B, C, D
    real*8, intent(out)   :: x, x1, x2, y
    logical, intent(in), optional :: only_real

    ! Local variables
    real*8                :: Q, Qprime, b1, c2, t, r, s, x0
    real*8, parameter     :: ONE_PLUS_EPS = 1. + 1E-12
    integer, parameter    :: MAX_ITER = 100
    integer               :: count

    x1 = 0
    x2 = 0
    y = 0

    if (A == 0) then
      x = d_qnan
      if (.not. only_real) call qdrtc(B, C, D, x1, x2, y)
    elseif (D == 0) then
      x = 0.
      if (.not. only_real) call qdrtc(A, B, C, x1, x2, y)
    else
      x = -(B / A) / 3
      call eval(x, A, B, C, D, Q, Qprime, b1, c2)
      t = Q / A
      r = abs(t)**(1./3.)
      s = sign(1.0d0, t)

      t = -Qprime / A
      if (t > 0) r = 1.324718 * max(r, sqrt(t))
      x0 = x - s * r
      if (x /= x0) then
        count = 0
        do while (count < MAX_ITER)
          count = count + 1
          x = x0
          call eval(x, A, B, C, D, Q, Qprime, b1, c2)
          if (Qprime == 0) then
            x0 = x
          else
            x0 = x - (Q / Qprime) / ONE_PLUS_EPS
          endif
          if (s * x0 <= s * x) exit
        enddo
        if (x /= 0 .and. abs(A) * x**2 > abs(D / x)) then
          c2 = -D / x
          b1 = (c2 - C) / x
        endif
      endif

      if (.not. only_real) call qdrtc(A, b1, c2, x1, x2, y)
    endif


  contains
    pure subroutine eval(x, A, B, C, D, Q, Qprime, b1, c2)
      implicit none

      real*8, intent(in)    :: x, A, B, C, D
      real*8, intent(out)   :: Q, Qprime, b1, c2

      ! Local variables
      real*8                :: q0

      q0 = A * x
      b1 = q0 + B
      c2 = b1 * x + C
      Qprime = (q0 + b1) * x + c2
      Q = c2 * x + D
    end subroutine

    pure subroutine qdrtc(A, B, C, x1, x2, y)
      implicit none

      real*8, intent(in)    :: A, B, C
      real*8, intent(out)   :: x1, x2, y

      ! Local variables
      real*8                :: roots(2)

      call polynomial_roots_deg2([A, B, C], roots, y)
      x1 = roots(1)
      x2 = roots(2)
    end subroutine
  end

end module
