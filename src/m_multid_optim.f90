module m_multid_optim
  use m_linesearch

  integer, parameter :: MAX_NR_LBFGS_DIRECTIONS = 5, MAX_NR_BROYDEN_DIRECTIONS = 10

  type optimOpts
    integer               :: maxFEval = 50       ! Maximum nr function evaluations
    real*8                :: fTol     = 0.0      ! Stop if f(x) < fTol (f is signed!)
    real*8                :: gTol     = 1E-12    ! Stop if |g(x)|_2 < gTol
    real*8                :: errTol   = 1E-12    ! Stop if |x - x^*|_2 < errTol (use an error estimate)
    real*8                :: hStep    = 1E-7     ! Step size used for FD approximation of gradient (depends on scaling of f, x)
    integer               :: fdOrder  = 2        ! 1/2
    logical               :: verbose  = .false.  ! Print info
  end type

  type optimInfo
    integer               :: nrFEvals
    integer               :: nrIterations
    integer               :: lineSearchFlag
    real*8                :: fVal
    real*8                :: normGrad
    real*8                :: step
    logical               :: success
  end type


contains
  subroutine optimize(x, fun, opts, info, fun_and_grad, ls_opts)
    implicit none

    real*8, intent(inout) :: x(1:)
    real*8, external      :: fun
    real*8, external, optional :: fun_and_grad
    type(optimOpts), intent(in) :: opts
    type(lsOpts), intent(in), optional :: ls_opts
    type(optimInfo), intent(out):: info

    ! Local variables
    type(lsOpts), parameter :: ls_opts_default = lsOpts()
    type(lsOpts)          :: ls_opts_
    integer               :: it, nr_fevals, nvars, ls_info, nr_directions, last_idx, ls_nr_fevals
    integer               :: actual_max_nr_directions
    real*8                :: fun_val, normGrad, errEst, step, rho(MAX_NR_LBFGS_DIRECTIONS), diag, prev_x(size(x, 1)), prev_fun_val
    real*8                :: search(size(x, 1)), grad(size(x, 1)), prev_grad(size(x, 1))
    real*8                :: Y(size(x, 1), MAX_NR_LBFGS_DIRECTIONS), S(size(x, 1), MAX_NR_LBFGS_DIRECTIONS)
    logical               :: hess_success

    nr_fevals = 0
    nr_directions = 0
    last_idx = 0  ! index of direction which was last added to Y, S
    nvars = size(x, 1)
    ls_info = MT_GOOD

    actual_max_nr_directions = min(MAX_NR_LBFGS_DIRECTIONS, nvars)

    fun_val = evaluate_fun_and_grad(grad, x)

    normGrad = norm2(grad)
    errEst = max(1.0, abs(opts%errTol) * 2.0)

    if (present(ls_opts)) then
      ls_opts_ = ls_opts
    else
      ls_opts_ = ls_opts_default
    endif

    it = 0
    if (opts%verbose) then
      write(*,'(A)') '/------------------------------------------------------------------\'
      write(*,'(A,I4)') '  Starting L-BFGS optimization with N = ', nvars
      write(*,'(A)') ''
      write(*,'(A)') '      Iter   FEvals      FunVal    norm(G)   Err. est.   Stepsize'
      write(*,'(A,I8,A,I6,A,1PD9.3,A,1PD9.3,A,1PD9.3,A,1PD9.3)') '  ', it, '   ', nr_fevals,& 
        '   ', fun_val, '   ', normGrad, '  ', errEst, '  ', 0.0
    endif

    hess_success = .true.
    info%success = .true.
    do while (fun_val > opts%fTol .and. normGrad > opts%gTol .and. errEst > opts%errTol)

      if (nr_fevals >= opts%maxFEval .or. ls_info/=MT_GOOD .or. (it > 0 .and. step==0.0)) then
        info%success = .false.
        if (it > 1 .and. prev_fun_val < fun_val) then
          fun_val = prev_fun_val
          x = prev_x
          grad = prev_grad
        endif
        exit
      endif

      it = it + 1

      if (it > 1) hess_success = update_hessian_inverse(Y, S, rho, diag, nr_directions, &
        last_idx, nvars, step, search, grad, prev_grad)
      if (.not. hess_success) then
        info%success = .false.
        exit
      endif

      ! Find search direction
      !   search = -hessian \ grad
      call apply_hessian_inverse(search, grad, Y, S, rho, diag, nr_directions, last_idx, nvars)
      search = -search

      ! Satisfy sufficient decrease condition (and curvature condition if MT is used)
      step = merge(min(1.0, 1.0 / normGrad), 1.0d0, it==1)
      prev_fun_val = fun_val
      prev_x = x
      prev_grad = grad
      select case(ls_opts_%type)
      case (LS_MORE_THUENTE)
        call more_thuente(evaluate_fun_and_directional_derivative, nvars, x, fun_val, dot_product(grad, search), &
          search, step, ls_opts_, opts%maxFEval - nr_fevals, ls_info)
        if (.not. present(fun_and_grad) .and. nvars > 1) call approx_gradient(grad, fun, fun_val, x)
      case (LS_ARMIJO_BACKTRACKING)
        call armijo_backtracking(fun, x, fun_val, grad, search, step, ls_opts_, opts%maxFEval - nr_fevals, ls_info, ls_nr_fevals)
        if (present(fun_and_grad)) then
          fun_val = evaluate_fun_and_grad(grad, x)
        else
          call approx_gradient(grad, fun, fun_val, x)
        endif
        nr_fevals = nr_fevals + ls_nr_fevals
      end select

      normGrad = norm2(grad)
      errEst = max(normGrad, norm2(search))
      if (opts%verbose) then
        write(*,'(A,I8,A,I6,A,1PD9.3,A,1PD9.3,A,1PD9.3,A,1PD9.3)') '  ', it, '   ', nr_fevals, &
          '   ', fun_val, '   ', normGrad, '  ', errEst, '  ', step
      endif

    enddo

    info%success = info%success .and. .not. isnan(fun_val)

    if (opts%verbose) then
      write(*,'(A)') ''
      if (info%success) then
        write(*,'(A)') '  L-BFGS was successful'
      else
        write(*,'(A,I1)') '  L-BFGS failed with line search flag = ', ls_info
      endif
      write(*,'(A)') '\------------------------------------------------------------------/'
      write(*,'(A)') ''
    endif

    info%step = step
    info%fVal = fun_val
    info%normGrad = normGrad
    info%nrFEvals = nr_fevals
    info%nrIterations = it
    info%lineSearchFlag = ls_info

  contains
    real*8 function evaluate_fun_and_grad(g, x) result(f)
      implicit none

      real*8, intent(out)   :: g(nvars)
      real*8, intent(in)    :: x(nvars)

      if (present(fun_and_grad)) then
        f = fun_and_grad(g, x)
      else
        f = fun(x)
        call approx_gradient(g, fun, f, x)
      endif
      nr_fevals = nr_fevals + 1
    end

    real*8 function evaluate_fun_and_directional_derivative(g, x, dir) result(f)
      implicit none

      real*8, intent(out)   :: g
      real*8, intent(in)    :: x(nvars), dir(nvars)

      ! Local variables
      real*8                :: normDir

      if (present(fun_and_grad)) then
        f = fun_and_grad(grad, x)
        g = dot_product(grad, dir)
      elseif (nvars == 1) then
        f = fun(x)
        ! NB this sets grad from optimize
        call approx_gradient(grad, fun, f, x)
        g = grad(1) * dir(1)
      else
        f = fun(x)
        normDir = 1.0!norm2(dir)
        if (opts%fdOrder == 1) then
          g = normDir * (fun(x + opts%hStep * dir / normDir) - f) / opts%hStep
        else
          g = normDir * (fun(x + opts%hStep * dir / normDir) - fun(x - opts%hStep * dir / normDir)) / (2 * opts%hStep)
        endif
        nr_fevals = nr_fevals + opts%fdOrder
      endif
      nr_fevals = nr_fevals + 1

    end

    subroutine approx_gradient(grad, fun, f0, x)
      implicit none

      real*8, intent(out):: grad(nvars)
      real*8, external   :: fun
      real*8, intent(in) :: x(nvars), f0

      ! Local variables
      integer          :: j
      real*8           :: fdir, fotherdir, xdir(nvars)

      do j=1,nvars
        xdir = x
        xdir(j) = xdir(j) + opts%hStep
        fdir = fun(xdir)
        if (opts%fdOrder == 1) then
          grad(j) = (fdir - f0) / opts%hStep
        else
          xdir(j) = xdir(j) - 2 * opts%hStep
          fotherdir = fun(xdir)
          grad(j) = (fdir - fotherdir) / (2 * opts%hStep)
        endif
      enddo
      nr_fevals = nr_fevals + opts%fdOrder * nvars
    end subroutine

    subroutine apply_hessian_inverse(p, rhs, Y, S, rho, diag, nr_directions, last_idx, nvars)
      implicit none

      real*8, intent(out)     :: p(nvars)
      real*8, intent(in)      :: rhs(nvars), Y(nvars, MAX_NR_LBFGS_DIRECTIONS), S(nvars, MAX_NR_LBFGS_DIRECTIONS)
      real*8, intent(in)      :: rho(MAX_NR_LBFGS_DIRECTIONS), diag
      integer, intent(in)   :: nr_directions, last_idx, nvars

      ! Local variables
      real                  :: alpha(MAX_NR_LBFGS_DIRECTIONS), beta
      integer               :: k, idx

      p = rhs
      if (nr_directions == 0) return

      idx = last_idx
      do k = 1,nr_directions
        if (idx < 1) idx = actual_max_nr_directions
        alpha(idx) = rho(idx) * dot_product(S(:, idx), p)
        p = p - alpha(idx) * Y(:, idx)
        idx = idx - 1
      enddo

      p = p * diag

      do k = 1,nr_directions
        idx = idx + 1
        if (idx > actual_max_nr_directions) idx = 1
        beta = rho(idx) * dot_product(Y(:, idx), p)
        p = p + (alpha(idx) - beta) * S(:, idx)
      enddo
    end subroutine

    logical function update_hessian_inverse(Y, S, rho, diag, nr_directions, last_idx, nvars, step, search,& 
      grad, prev_grad) result(success)
      implicit none

      real*8, intent(inout)   :: Y(nvars, MAX_NR_LBFGS_DIRECTIONS), S(nvars, MAX_NR_LBFGS_DIRECTIONS), &
        rho(MAX_NR_LBFGS_DIRECTIONS), diag
      integer, intent(inout):: nr_directions, last_idx
      integer, intent(in)   :: nvars
      real*8, intent(in)      :: step, search(nvars), grad(nvars), prev_grad(nvars)

      ! Local variables
      real*8                  :: sy_dot, yy_dot

      last_idx = last_idx + 1
      if (last_idx > actual_max_nr_directions) last_idx = 1

      S(:, last_idx) = step * search
      Y(:, last_idx) = grad - prev_grad
      sy_dot = dot_product(S(:, last_idx), Y(:, last_idx))
      yy_dot = dot_product(Y(:, last_idx), Y(:, last_idx))
      success = sy_dot /= 0.0 .and. yy_dot /= 0.0

      if (.not. success) return

      rho(last_idx) = 1.0 / sy_dot
      diag = sy_dot / yy_dot

      nr_directions = min(actual_max_nr_directions, nr_directions + 1)
    end function
  end subroutine
end module