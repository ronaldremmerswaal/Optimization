module m_oned_rootfinding
  abstract interface
    real*8 function der_optim_fun(x, fun)
      real*8, intent(in)  :: x
      real*8, intent(out), optional :: fun
    end function

    real*8 function optim_fun(x)
      real*8, intent(in)  :: x
    end function
  end interface
contains
  real*8 function brent_min(dfun, x0, xTol, maxIt, verbose, maxStep) result(x_val)
    implicit none

    procedure(der_optim_fun) :: dfun
    real*8, intent(in)    :: x0, xTol
    integer, intent(in)   :: maxIt
    logical, intent(in), optional :: verbose
    real*8, intent(in), optional :: maxStep
    
    ! Local variables
    real*8                :: dfun_val, dfun_prev, x_prev, step, maxStep_, fun_val
    integer               :: it
    logical               :: verbose_, converged

    verbose_ = merge(verbose, .false., present(verbose))
    maxStep_ = merge(maxStep, 1D100, present(maxStep))

    x_prev = x0
    x_val = x0
    dfun_val = dfun(x_val, fun_val)
    if (dfun_val==0) return
    dfun_prev = dfun_val

    if (verbose_) then
      write(*,'(A)') '/---------------------------------------------------\'
      write(*,'(A,I4)') '           Starting Brent min search'
      write(*,'(A)') ''
      write(*,'(A)') '      Iter       Sol        DFunVal    Err. est.'
    endif

    converged = .false.
    it = 1
    ! Find bracket using Newton
    do while (dfun_val * dfun_prev > 0 .and. it < maxIt)
      x_prev = x_val
      dfun_prev = dfun_val

      step = - 2 * fun_val / dfun_val
      if (present(maxStep) .and. abs(step) > maxStep_) then
        step = sign(maxStep, step)
      endif
      x_val = x_prev + step
      
      converged = abs(step) < xTol
      if (converged) exit

      dfun_val = dfun(x_val, fun_val)
      
      it = it + 1

      if (verbose_) write(*,'(A,I8,A,1PD10.3,A,1PD10.3,A,1PD10.3)') '  ', it, '   ', &
      x_val, '   ', dfun_val, '  ', abs(step)

    enddo

    if (.not. converged) then
      if (dfun_val * dfun_prev > 0) return
      if (verbose_) write(*,'(A,I4)') '               Continuing in brent...'

      ! Then apply Brent to derivative
      x_val = brent(optim_wrapper, x_val, x_prev, xTol, maxIt - it, dfun_val, dfun_prev, verbose = verbose)
    endif
  contains
    real*8 function optim_wrapper(x) result(ans)
      real*8, intent(in) :: x
      ans = dfun(x)
    end function
  end function
  
  real*8 function newton(dfun, x0, xTol, maxIt, verbose, maxStep) result(x_val)
    implicit none

    procedure(der_optim_fun) :: dfun
    real*8, intent(in)    :: x0, xTol
    integer, intent(in)   :: maxIt
    logical, intent(in), optional :: verbose
    real*8, intent(in), optional :: maxStep
    
    ! Local variables
    real*8                :: dfun_val, x_prev, step, maxStep_, fun_val
    integer               :: it
    logical               :: verbose_, converged

    verbose_ = merge(verbose, .false., present(verbose))
    maxStep_ = merge(maxStep, 1D100, present(maxStep))

    x_prev = x0
    x_val = x0
    dfun_val = dfun(x_val, fun_val)
    if (fun_val==0) return

    if (verbose_) then
      write(*,'(A)') '/---------------------------------------------------\'
      write(*,'(A,I4)') '               Starting Newton'
      write(*,'(A)') ''
      write(*,'(A)') '      Iter       Sol         FunVal    Err. est.'
    endif

    converged = .false.
    ! Find bracket using Newton
    do it=2,maxIt
      x_prev = x_val

      step = - fun_val / dfun_val
      if (present(maxStep) .and. abs(step) > maxStep_) then
        step = sign(maxStep, step)
      endif
      x_val = x_prev + step

      if (verbose_) write(*,'(A,I8,A,1PD10.3,A,1PD10.3,A,1PD10.3)') '  ', it, '   ', &
      x_val, '   ', fun_val, '  ', abs(step)

      converged = abs(step) < xTol .or. fun_val==0
      if (converged) exit
      
      dfun_val = dfun(x_val, fun_val)

    enddo
  end function

  real*8 function brent(fun, xaIn, xbIn, xTol, maxIt, faIn, fbIn, verbose) result(xsol)
    implicit none

    procedure(optim_fun) :: fun
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

    iter = 0

    xa = xaIn
    xb = xbIn
    if (present(faIn)) then
      fa = faIn
    else
      fa = fun(xa)
      iter = iter + 1
    endif
    if (present(fbIn)) then
      fb = fbIn
    else
      fb = fun(xb)
      iter = iter + 1
    endif

    if (fa == 0) then
      xsol = xa
      return
    elseif (fb == 0) then
      xsol = xb
      return
    elseif (fa * fb > 0) then
      xsol = 0
      print*, 'Error: brent only possible if function values switch sign'
      return
    endif

    if (verbose_) then
      write(*,'(A)') '/---------------------------------------------------\'
      write(*,'(A,I4)') '               Starting Brent search'
      write(*,'(A)') ''
      write(*,'(A)') '      Iter       Sol         FunVal    Err. est.'
    endif

    xc = xa;  fc = fa
    xd = xb - xa;  xe = xd

    flag = 0
    do while (fb /= 0 .and. xa /= xb)

      ! Ensure that b is the best result so far, a is the previous
      ! value of b, and c is on the opposite side of the zero from b.
      if (fb * fc > 0) then
          xc = xa;  fc = fa
          xd = xb - xa;  xe = xd
      endif
      if (abs(fc) < abs(fb)) then
          xa = xb;  xb = xc;  xc = xa
          fa = fb;  fb = fc;  fc = fa
      endif

      ! Convergence test and possible exit
      xm = 0.5*(xc - xb)
      if (verbose_) write(*,'(A,I8,A,1PD10.3,A,1PD10.3,A,1PD10.3)') '  ', iter, '   ', &
        xc, '   ', fc, '  ', abs(xm)

      toler = xTol
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
            p = 2*xm*s
            q = 1 - s
          else
            ! Inverse quadratic interpolation
            q = fa/fc
            r = fb/fc
            p = s*(2*xm*q*(q - r) - (xb - xa)*(r - 1))
            q = (q - 1)*(r - 1)*(s - 1)
          endif
          if (p > 0) then
            q = -q
          else
            p = -p
          endif
          ! Is interpolated point acceptable
          if ((2*p < 3*xm*q - abs(toler*q)) .and. (p < abs(0.5*xe*q))) then
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

    end do

    if (verbose_) write(*,'(A)') '\---------------------------------------------------/'

    xsol = xb

  end function

  real*8 function bisection(fun, xLeft, xRight, xTol, maxIt, fLeft, fRight) result(x)
  implicit none

  procedure(optim_fun)  :: fun
  real*8, intent(in)    :: xLeft, xRight, xTol
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