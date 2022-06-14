module m_linesearch

  integer, parameter :: LS_ARMIJO_BACKTRACKING = 0, LS_MORE_THUENTE = 1, LS_NONE = -1
  
  integer, parameter :: MT_BADINPUT = 0
  integer, parameter :: MT_GOOD = 1
  integer, parameter :: MT_INTERVAL = 2
  integer, parameter :: MT_MAX_FUN_EVAL = 3
  integer, parameter :: MT_STEP_TOO_SMALL = 4
  integer, parameter :: MT_STEP_TOO_LARGE = 5
  integer, parameter :: MT_ROUNDING_ERRORS = 6
  integer, parameter :: MT_NO_DESCENT_DIRECTION = 7


  type lsOpts
    integer               :: type               = LS_MORE_THUENTE
    real*8                :: decreaseCondition  = 1E-2
    real*8                :: curvatureCondition = 5E-1
    real*8                :: stepFactor         = 6E-1
    real*8                :: xTol               = 1E-16
    real*8                :: stpMin             = 1E-16
    real*8                :: stpMax             = 1E+16
  end type

contains
  subroutine more_thuente(fun,n,x,f,dginit,s,stp,opts,maxfev,info)
  implicit none

  integer, intent(in) :: n
  real*8, intent(in)   :: s(n), dginit
  real*8, intent(inout):: f,stp
  real*8, intent(inout):: x(n)
  integer, intent(out):: info
  real*8, external     ::  fun
  integer, intent(in) :: maxfev
  type(lsOpts), intent(in) :: opts



!     &*********
!
!     Subroutine cvsrch
!
!     The purpose of cvsrch is to find a step which satisfies
!     a sufficient decrease condition and a curvature condition.
!     The user must provide a subroutine which calculates the
!     function and the gradient.
!
!     At each stage the subroutine updates an interval of
!     uncertainty with endpoints stx and sty. The interval of
!     uncertainty is initially chosen so that it contains a
!     minimizer of the modified function
!
!          f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).
!
!     If a step is obtained for which the modified function
!     has a nonpositive function value and nonnegative derivative,
!     then the interval of uncertainty is chosen so that it
!     contains a minimizer of f(x+stp*s).
!
!     The algorithm is designed to find a step which satisfies
!     the sufficient decrease condition
!
!           f(x+stp*s) .le. f(x) + ftol*stp*(gradf(x)'s),
!
!     and the curvature condition
!
!           abs(gradf(x+stp*s)'s)) .le. gtol*abs(gradf(x)'s).
!
!     If ftol is less than gtol and if, for example, the function
!     is bounded below, then there is always a step which satisfies
!     both conditions. If no step can be found which satisfies both
!     conditions, then the algorithm usually stops when rounding
!     errors prevent further progress. In this case stp only
!     satisfies the sufficient decrease condition.
!
!     The subroutine statement is
!
!        subroutine cvsrch(fcn,n,x,f,g,s,stp,ftol,gtol,xtol,
!                          stpmin,stpmax,maxfev,info,nfev,wa)
!     where
!
!	fcn is the name of the user-supplied subroutine which
!         calculates the function and the gradient.  fcn must
!      	  be declared in an external statement in the user
!         calling program, and should be written as follows.
!
!	  subroutine fcn(n,x,f,g)
!         integer n
!         f
!         x(n),g(n)
!	  ----------
!         Calculate the function at x and
!         return this value in the variable f.
!         Calculate the gradient at x and
!         return this vector in g.
!	  ----------
!	  return
!	  end
!
!       n is a positive integer input variable set to the number
!	  of variables.
!
!	x is an array of length n. On input it must contain the
!	  base point for the line search. On output it contains
!         x + stp*s.
!
!	f is a variable. On input it must contain the value of f
!         at x. On output it contains the value of f at x + stp*s.
!
!	g is an array of length n. On input it must contain the
!         gradient of f at x. On output it contains the gradient
!         of f at x + stp*s.
!
!	s is an input array of length n which specifies the
!         search direction.
!
!	stp is a nonnegative variable. On input stp contains an
!         initial estimate of a satisfactory step. On output
!         stp contains the final estimate.
!
!       ftol and gtol are nonnegative input variables. Termination
!         occurs when the sufficient decrease condition and the
!         directional derivative condition are satisfied.
!
!	xtol is a nonnegative input variable. Termination occurs
!         when the relative width of the interval of uncertainty
!	  is at most xtol.
!
!	stpmin and stpmax are nonnegative input variables which
!	  specify lower and upper bounds for the step.
!
!	maxfev is a positive integer input variable. Termination
!         occurs when the number of calls to fcn is at least
!         maxfev by the end of an iteration.
!
!	info is an integer output variable set as follows:
!
!	  info = 0  Improper input parameters.
!
!	  info = 1  The sufficient decrease condition and the
!                   directional derivative condition hold.
!
!	  info = 2  Relative width of the interval of uncertainty
!		    is at most xtol.
!
!	  info = 3  Number of calls to fcn has reached maxfev.
!
!	  info = 4  The step is at the lower bound stpmin.
!
!	  info = 5  The step is at the upper bound stpmax.
!
!	  info = 6  Rounding errors prevent further progress.
!                   There may not be a step which satisfies the
!                   sufficient decrease and curvature conditions.
!                   Tolerances may be too small.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn.
!
!	wa is a work array of length n.
!
!     Subprograms called
!
!	user-supplied......fcn
!
!	MINPACK-supplied...cstep
!
!	FORTRAN-supplied...abs,max,min
!
!     Argonne National Laboratory. MINPACK Project. June 1983
!     Jorge J. More', David J. Thuente
!
!     &*********
  integer infoc,nfev
  logical brackt,stage1
  real*8 dg,dgm,dgtest,dgx,dgxm,dgy,dgym,   &
  &       finit,ftest1,fm,fx,fxm,fy,fym,p5,p66,stx,sty,       &
  &       stmin,stmax,width,width1,xtrapf
  real*8              :: x0(n)
  real*8                :: ftol,gtol,xtol,stpmin,stpmax

  data p5,p66,xtrapf/0.5,0.66,4.0/
  info = 0
  infoc = 1


  ftol=opts%decreaseCondition
  gtol=opts%curvatureCondition
  xtol=opts%xTol
  stpmin=opts%stpMin
  stpmax=opts%stpMax
  stpmax=opts%stpMax
!
!     Check the input parameters for errors.
!
  if (n .le. 0 .or. stp .le. 0.0 .or. ftol .lt. 0.0 .or.     &
  &    gtol .lt. 0.0 .or. xtol .lt. 0.0 .or. stpmin .lt. 0.0  &
  &    .or. stpmax .lt. stpmin .or. maxfev .le. 0) return


!
!     Compute the initial gradient in the search direction
!     and check that s is a descent direction.
!
    if (dginit >= 0.0) then
      info = 7
      return
    endif


!     Initialize local variables.
!
  brackt = .false.
  stage1 = .true.
  nfev = 0
  finit = f
  dgtest = ftol*dginit
  width = stpmax - stpmin
  width1 = width/0.5d0
  x0 = x
!
!     The variables stx, fx, dgx contain the values of the step,
!     function, and directional derivative at the best step.
!     The variables sty, fy, dgy contain the value of the step,
!     function, and derivative at the other endpoint of
!     the interval of uncertainty.
!     The variables stp, f, dg contain the values of the step,
!     function, and derivative at the current step.
!
  stx = 0.0
  fx = finit
  dgx = dginit
  sty = 0.0
  fy = finit
  dgy = dginit
!
!     Start of iteration.
!
30 continue
!
!        Set the minimum and maximum steps to correspond
!        to the present interval of uncertainty.
!
      if (brackt) then
        stmin = min(stx,sty)
        stmax = max(stx,sty)
      else
        stmin = stx
        stmax = stp + xtrapf*(stp - stx)
      endif
!
!        Force the step to be within the bounds stpmax and stpmin.
!
      stp = max(stp,stpmin)
      stp = min(stp,stpmax)
!
!        If an unusual termination is to occur then let
!        stp be the lowest point obtained so far.
!
      if ((brackt .and. (stp .le. stmin .or. stp .ge. stmax))      &
  &      .or. nfev .ge. maxfev-1 .or. infoc .eq. 0                 &
  &      .or. (brackt .and. stmax-stmin .le. xtol*stmax)) stp = stx
!
!        Evaluate the function and gradient at stp
!        and compute the directional derivative.
!
      x = x0 + stp*s
      f = fun(dg, x, s)

      nfev = nfev + 1
      ftest1 = finit + stp*dgtest
!
!        Test for convergence.
!
      if ((brackt .and. (stp .le. stmin .or. stp .ge. stmax)) .or. infoc .eq. 0) info = 6
      if (stp .eq. stpmax .and. f .le. ftest1 .and. dg .le. dgtest) info = 5
      if (stp .eq. stpmin .and. (f .gt. ftest1 .or. dg .ge. dgtest)) info = 4
      if (nfev .ge. maxfev) info = 3
      if (brackt .and. stmax-stmin .le. xtol*stmax) info = 2
      if (f .le. ftest1 .and. abs(dg) .le. gtol*(-dginit)) info = 1
!
!        Check for termination.
!
      if (info .ne. 0) return
!
!        In the first stage we seek a step for which the modified
!        function has a nonpositive value and nonnegative derivative.
!
      if (stage1 .and. f .le. ftest1 .and. dg .ge. min(ftol,gtol)*dginit) stage1 = .false.
!
!        A modified function is used to predict the step only if
!        we have not obtained a step for which the modified
!        function has a nonpositive function value and nonnegative
!        derivative, and if a lower function value has been
!        obtained but the decrease is not sufficient.
!
      if (stage1 .and. f .le. fx .and. f .gt. ftest1) then
!
!           Define the modified function and derivative values.
!
        fm = f - stp*dgtest
        fxm = fx - stx*dgtest
        fym = fy - sty*dgtest
        dgm = dg - dgtest
        dgxm = dgx - dgtest
        dgym = dgy - dgtest
!
!           Call cstep to update the interval of uncertainty
!           and to compute the new step.
!
        call cstep(stx,fxm,dgxm,sty,fym,dgym,stp,fm,dgm,brackt,stmin,stmax,infoc)
!
!           Reset the function and gradient values for f.
!
        fx = fxm + stx*dgtest
        fy = fym + sty*dgtest
        dgx = dgxm + dgtest
        dgy = dgym + dgtest
      else
!
!           Call cstep to update the interval of uncertainty
!           and to compute the new step.
!
        call cstep(stx,fx,dgx,sty,fy,dgy,stp,f,dg,brackt,stmin,stmax,infoc)
        end if
!
!        Force a sufficient decrease in the size of the
!        interval of uncertainty.
!
      if (brackt) then
        if (abs(sty-stx) .ge. p66*width1) stp = stx + p5*(sty - stx)
        width1 = width
        width = abs(sty-stx)
        end if
!
!        End of iteration.
!
      go to 30
!
!     Last card of subroutine cvsrch.
!
  end
  subroutine cstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,stpmin,stpmax,info)
  integer info
  real*8 stx,fx,dx,sty,fy,dy,stp,fp,dp,stpmin,stpmax
  logical brackt,bound
!     &*********
!
!     Subroutine cstep
!
!     The purpose of cstep is to compute a safeguarded step for
!     a linesearch and to update an interval of uncertainty for
!     a minimizer of the function.
!
!     The parameter stx contains the step with the least function
!     value. The parameter stp contains the current step. It is
!     assumed that the derivative at stx is negative in the
!     direction of the step. If brackt is set true then a
!     minimizer has been bracketed in an interval of uncertainty
!     with endpoints stx and sty.
!
!     The subroutine statement is
!
!       subroutine cstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
!                        stpmin,stpmax,info)
!
!     where
!
!       stx, fx, and dx are variables which specify the step,
!         the function, and the derivative at the best step obtained
!         so far. The derivative must be negative in the direction
!         of the step, that is, dx and stp-stx must have opposite
!         signs. On output these parameters are updated appropriately.
!
!       sty, fy, and dy are variables which specify the step,
!         the function, and the derivative at the other endpoint of
!         the interval of uncertainty. On output these parameters are
!         updated appropriately.
!
!       stp, fp, and dp are variables which specify the step,
!         the function, and the derivative at the current step.
!         If brackt is set true then on input stp must be
!         between stx and sty. On output stp is set to the new step.
!
!       brackt is a logical variable which specifies if a minimizer
!         has been bracketed. If the minimizer has not been bracketed
!         then on input brackt must be set false. If the minimizer
!         is bracketed then on output brackt is set true.
!
!       stpmin and stpmax are input variables which specify lower
!         and upper bounds for the step.
!
!       info is an integer output variable set as follows:
!         If info = 1,2,3,4,5, then the step has been computed
!         according to one of the five cases below. Otherwise
!         info = 0, and this indicates improper input parameters.
!
!     Subprograms called
!
!       FORTRAN-supplied ... abs,max,min,sqrt
!                        ... dble
!
!     Argonne National Laboratory. MINPACK Project. June 1983
!     Jorge J. More', David J. Thuente
!
!     &*********
  real*8 gamma,p,p66,q,r,s,sgnd,stpc,stpf,stpq,theta
  data p66 /0.66/
  info = 0
!
!     Check the input parameters for errors.
!
  if ((brackt .and. (stp .le. min(stx,sty) .or.                   &
  &    stp .ge. max(stx,sty))) .or.                                &
  &    dx*(stp-stx) .ge. 0.0 .or. stpmax .lt. stpmin) return
!
!     Determine if the derivatives have opposite sign.
!
  sgnd = dp*(dx/abs(dx))
!
!     First case. A higher function value.
!     The minimum is bracketed. If the cubic step is closer
!     to stx than the quadratic step, the cubic step is taken,
!     else the average of the cubic and quadratic steps is taken.
!
  if (fp .gt. fx) then
      info = 1
      bound = .true.
      theta = 3*(fx - fp)/(stp - stx) + dx + dp
      s = max(abs(theta),abs(dx),abs(dp))
      gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
      if (stp .lt. stx) gamma = -gamma
      p = (gamma - dx) + theta
      q = ((gamma - dx) + gamma) + dp
      r = p/q
      stpc = stx + r*(stp - stx)
      stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/2)*(stp - stx)
      if (abs(stpc-stx) .lt. abs(stpq-stx)) then
        stpf = stpc
      else
        stpf = stpc + (stpq - stpc)/2
        end if
      brackt = .true.
!
!     Second case. A lower function value and derivatives of
!     opposite sign. The minimum is bracketed. If the cubic
!     step is closer to stx than the quadratic (secant) step,
!     the cubic step is taken, else the quadratic step is taken.
!
  else if (sgnd .lt. 0.0) then
      info = 2
      bound = .false.
      theta = 3*(fx - fp)/(stp - stx) + dx + dp
      s = max(abs(theta),abs(dx),abs(dp))
      gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
      if (stp .gt. stx) gamma = -gamma
      p = (gamma - dp) + theta
      q = ((gamma - dp) + gamma) + dx
      r = p/q
      stpc = stp + r*(stx - stp)
      stpq = stp + (dp/(dp-dx))*(stx - stp)
      if (abs(stpc-stp) .gt. abs(stpq-stp)) then
        stpf = stpc
      else
        stpf = stpq
        end if
      brackt = .true.
!
!     Third case. A lower function value, derivatives of the
!     same sign, and the magnitude of the derivative decreases.
!     The cubic step is only used if the cubic tends to infinity
!     in the direction of the step or if the minimum of the cubic
!     is beyond stp. Otherwise the cubic step is defined to be
!     either stpmin or stpmax. The quadratic (secant) step is also
!     computed and if the minimum is bracketed then the the step
!     closest to stx is taken, else the step farthest away is taken.
!
  else if (abs(dp) .lt. abs(dx)) then
      info = 3
      bound = .true.
      theta = 3*(fx - fp)/(stp - stx) + dx + dp
      s = max(abs(theta),abs(dx),abs(dp))
!
!        The case gamma = 0 only arises if the cubic does not tend
!        to infinity in the direction of the step.
!
      gamma = s*sqrt(max(0.,(theta/s)**2 - (dx/s)*(dp/s)))
      if (stp .gt. stx) gamma = -gamma
      p = (gamma - dp) + theta
      q = (gamma + (dx - dp)) + gamma
      r = p/q
      if (r .lt. 0.0 .and. gamma .ne. 0.0) then
        stpc = stp + r*(stx - stp)
      else if (stp .gt. stx) then
        stpc = stpmax
      else
        stpc = stpmin
        end if
      stpq = stp + (dp/(dp-dx))*(stx - stp)
      if (brackt) then
        if (abs(stp-stpc) .lt. abs(stp-stpq)) then
            stpf = stpc
        else
            stpf = stpq
            end if
      else
        if (abs(stp-stpc) .gt. abs(stp-stpq)) then
            stpf = stpc
        else
            stpf = stpq
            end if
        end if
!
!     Fourth case. A lower function value, derivatives of the
!     same sign, and the magnitude of the derivative does
!     not decrease. If the minimum is not bracketed, the step
!     is either stpmin or stpmax, else the cubic step is taken.
!
  else
      info = 4
      bound = .false.
      if (brackt) then
        theta = 3*(fp - fy)/(sty - stp) + dy + dp
        s = max(abs(theta),abs(dy),abs(dp))
        gamma = s*sqrt((theta/s)**2 - (dy/s)*(dp/s))
        if (stp .gt. sty) gamma = -gamma
        p = (gamma - dp) + theta
        q = ((gamma - dp) + gamma) + dy
        r = p/q
        stpc = stp + r*(sty - stp)
        stpf = stpc
      else if (stp .gt. stx) then
        stpf = stpmax
      else
        stpf = stpmin
        end if
      end if
!
!     Update the interval of uncertainty. This update does not
!     depend on the new step or the case analysis above.
!
  if (fp .gt. fx) then
      sty = stp
      fy = fp
      dy = dp
  else
      if (sgnd .lt. 0.0) then
        sty = stx
        fy = fx
        dy = dx
        end if
      stx = stp
      fx = fp
      dx = dp
      end if
!
!     Compute the new step and safeguard it.
!
  stpf = min(stpmax,stpf)
  stpf = max(stpmin,stpf)
  stp = stpf
  if (brackt .and. bound) then
      if (sty .gt. stx) then
        stp = min(stx+p66*(sty-stx),stp)
      else
        stp = max(stx+p66*(sty-stx),stp)
        end if
      end if
  return
!
!     Last card of subroutine cstep.
!
  end

  subroutine armijo_backtracking(fun, x, f_val, grad_val, dir, step, opts, maxIt, info, nr_fevals)
    implicit none
  
    real*8, external          :: fun
    real*8, intent(inout)     :: x(1:)
    real*8, intent(in)        :: grad_val(1:), dir(1:)
    real*8, intent(inout)     :: f_val, step
    type(lsOpts), intent(in) :: opts
    integer, intent(in)     :: maxIt
    integer, intent(out)    :: info, nr_fevals
  
    ! Local variables
    real*8                    :: f_val0, dir_dot_grad, x_new(size(x, 1)), fHist(3), stepHist(3), numerator, denominator
    integer                 :: it, hdx, nrHist
  
    f_val0 = f_val
    dir_dot_grad = dot_product(dir, grad_val)
    if (dir_dot_grad >= 0.0) then
      info = MT_NO_DESCENT_DIRECTION
      return
    endif
    nr_fevals = 0
  
    nrHist = 1
    hdx = 1
    stepHist(hdx) = 0.0
    fHist(hdx) = f_val
  
    info = MT_MAX_FUN_EVAL
    do it = 1,maxIt
      x_new = x + step * dir
      f_val = fun(x_new)
      nr_fevals = nr_fevals + 1
  
      if (nrHist == 3 .and. minval(fHist) < f_val) then
        nrHist = 0
        f_val = fHist(hdx)
        step = stepHist(hdx)
      endif
  
      nrHist = min(3, nrHist + 1)
      hdx = merge(1, hdx + 1, hdx == 3)
      stepHist(hdx) = step
      fHist(hdx) = f_val
  
      if (f_val < f_val0 + step * opts%decreaseCondition * dir_dot_grad) then
        info = MT_GOOD
        exit
      endif
  
      step = step * opts%stepFactor
      if (nrHist == 3) then
        ! Quadratic fit using previous pairs of (step, f(x + step * dir))
        numerator = fHist(1)*stepHist(2)**2 - fHist(2)*stepHist(1)**2 - fHist(1)*stepHist(3)**2 + &
          fHist(3)*stepHist(1)**2 + fHist(2)*stepHist(3)**2 - fHist(3)*stepHist(2)**2
        denominator = 2*(fHist(1)*stepHist(2) - fHist(2)*stepHist(1) - fHist(1)*stepHist(3) + &
          fHist(3)*stepHist(1) + fHist(2)*stepHist(3) - fHist(3)*stepHist(2))
        if (numerator > 0.0 .and. denominator > 0.0) then
          ! Ensure that we get a local minimum
          step = numerator / denominator
        else
          nrHist = 0
        endif
      endif
    enddo
  
    x = x_new
  
  end subroutine armijo_backtracking
  
end module