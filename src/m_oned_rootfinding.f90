module m_oned_rootfinding
  
contains
  real*8 function brent_min(fun, dfun, x0, xTol, maxIt, verbose) result(x_val)
    implicit none

    real*8, external      :: fun, dfun
    real*8, intent(in)    :: x0, xTol
    integer, intent(in)   :: maxIt
    logical, intent(in), optional :: verbose
    
    ! Local variables
    real*8                :: dfun_val, dfun_prev, x_prev
    integer               :: it
    logical               :: verbose_, converged

    verbose_ = merge(verbose, .false., present(verbose))

    x_prev = x0
    x_val = x0
    dfun_prev = dfun(x0)
    dfun_val = dfun_prev

    if (verbose_) then
      write(*,'(A)') '/---------------------------------------------------\'
      write(*,'(A,I4)') '           Starting Brent min search'
      write(*,'(A)') ''
      write(*,'(A)') '      Iter       Sol        DFunVal    Err. est.'
    endif

    converged = .false.
    it = 0
    ! Find bracket using Newton
    do while (dfun_val * dfun_prev > 0 .and. it <= maxIt)
      x_prev = x_val
      dfun_prev = dfun_val

      x_val = x_prev - fun(x_prev) / (2*dfun_prev)
      dfun_val = dfun(x_val)
      
      it = it + 1

      if (verbose_) write(*,'(A,I8,A,1PD10.3,A,1PD10.3,A,1PD9.3)') '  ', it, '   ', &
      x_val, '   ', dfun_val, '  ', abs(x_val - x_prev)

      converged = abs(x_val - x_prev) < xTol
      if (converged) exit
    enddo

    if (.not. converged) then
      if (verbose_) write(*,'(A,I4)') '               Continuing in brent...'

      ! Then apply Brent to derivative
      x_val = brent(dfun, x_val, x_prev, xTol, maxIt - it, dfun_val, dfun_prev, verbose = verbose)
    endif
  end function

  real*8 function brent(fun, xaIn, xbIn, xTol, maxIt, faIn, fbIn, verbose) result(xsol)
  implicit none

  real*8, external        :: fun
  real*8, intent(in)      :: xaIn, xbIn, xTol
  integer, intent(in)   :: maxIt
  real*8, intent(in), optional :: faIn, fbIn
  logical, intent(in), optional :: verbose

  ! Local variables
  real*8                :: xa, xb, fa, fb, fc, xc, xd, xe, xm, p, q, r, s
  real*8                :: toler
  integer               :: iter, flag
  logical               :: verbose_

  verbose_ = merge(verbose, .false., present(verbose))

  xa = xaIn
  xb = xbIn
  if (present(faIn)) then
    fa = faIn
  else
    fa = fun(xa)
  endif
  if (present(fbIn)) then
    fb = fbIn
  else
    fb = fun(xb)
  endif
  fc = fb

  if (fa == 0.0) then
    xsol = xa
    return
  elseif (fb == 0.0) then
    xsol = xb
    return
  elseif (fa * fb > 0.0) then
    xsol = 0.0
    print*, 'Error: brent only possible if function values switch sign'
    return
  endif

  if (verbose_) then
    write(*,'(A)') '/---------------------------------------------------\'
    write(*,'(A,I4)') '               Starting Brent search'
    write(*,'(A)') ''
    write(*,'(A)') '      Iter       Sol         FunVal    Err. est.'
  endif

  iter = 2
  flag = 0
  do while (fb /= 0.0 .and. xa /= xb)

      ! Ensure that b is the best result so far, a is the previous
      ! value of b, and c is on the opposite side of the zero from b.
      if (fb * fc > 0.0) then
          xc = xa;  fc = fa
          xd = xb - xa;  xe = xd
      endif
      if (abs(fc) < abs(fb)) then
          xa = xb;  xb = xc;  xc = xa
          fa = fb;  fb = fc;  fc = fa
      endif

      ! Convergence test and possible exit
      xm = 0.5*(xc - xb)
      toler = 2.0*xTol*max(abs(xb),1.0)
      if ((abs(xm) <= toler) .or. (fb == 0.0)) then
          flag = 0
          exit
      elseif (iter == maxIt) then
          flag = 2
          exit
      endif

      ! Choose bisection or interpolation
      if ((abs(xe) < toler) .or. (abs(fa) <= abs(fb))) then
          ! Bisection
          xd = xm;  xe = xm
      else
          ! Interpolation
          s = fb/fa;
          if (xa == xc) then
            ! Linear interpolation
            p = 2.0*xm*s
            q = 1.0 - s
          else
            ! Inverse quadratic interpolation
            q = fa/fc
            r = fb/fc
            p = s*(2.0*xm*q*(q - r) - (xb - xa)*(r - 1.0))
            q = (q - 1.0)*(r - 1.0)*(s - 1.0)
          endif
          if (p > 0.0) then
            q = -q
          else
            p = -p
          endif
          ! Is interpolated point acceptable
          if ((2.0*p < 3.0*xm*q - abs(toler*q)) .and. (p < abs(0.5*xe*q))) then
            xe = xd;  xd = p/q
          else
            xd = xm;  xe = xm
          endif
      endif ! Interpolation

      ! Next point
      xa = xb
      fa = fb
      if (abs(xd) > toler) then
        xb = xb + xd
      elseif (xb > xc) then
        xb = xb - toler
      else
        xb = xb + toler
      endif
      fb = fun(xb)
      iter = iter + 1

      if (verbose_) write(*,'(A,I8,A,1PD10.3,A,1PD10.3,A,1PD9.3)') '  ', iter, '   ', &
          xb, '   ', fb, '  ', abs(xc - xb)

    end do

    if (verbose_) write(*,'(A)') '\---------------------------------------------------/'

    if (flag == 0) then
      xsol = xb
    else
      ! We guarantee convergence
      xsol = bisection(fun, xb, xc, xTol, maxIt, fb, fc)
    endif

  end function

  real*8 function bisection(fun, xLeft, xRight, xTol, maxIt, fLeft, fRight) result(x)
  implicit none

  real*8, external        :: fun
  real*8, intent(in)      :: xLeft, xRight, xTol
  integer, intent(in)   :: maxIt
  real*8, optional, intent(in) :: fLeft, fRight

  ! Local variables
  real*8                  :: fLeft_, fRight_, xMid, fMid, xRight_, xLeft_
  integer               :: it

  if (present(fLeft)) then
    fLeft_ = fLeft
  else
    fLeft_ = fun(xLeft)
  endif
  if (present(fRight)) then
    fRight_ = fRight
  else
    fRight_ = fun(xRight)
  endif
  if (xRight < xLeft) then
    xMid = xRight
    xRight_ = xLeft
    xLeft_ = xMid
    fMid = fRight
    fRight_ = fLeft
    fLeft_ = fMid
  else
    xLeft_ = xLeft
    xRight_ = xRight
  endif

  if (fLeft_ == 0.0) then
    x = xLeft_
    return
  elseif (fRight_ == 0.0) then
    x = xRight_
    return
  elseif (fLeft_ * fRight_ > 0.0) then
    x = 0.0
    print*, 'Error: bisection only possible if function values switch sign'
    return
  endif

  it = 0
  do while (xRight_ - xLeft_ > xTol)
    it = it + 1
    xMid = (xLeft_ + xRight_) / 2
    fMid = fun(xMid)
    if (fMid == 0.0) then
      x = xMid
      return
    elseif (fMid * fLeft_ < 0.0) then
      xRight_ = xMid
      fRight_ = fMid
    else
      xLeft_ = xMid
      fLeft_ = fMid
    endif

    if (it == maxIt) exit
  enddo

  x = (xLeft_ + xRight_) / 2
  end function

end module