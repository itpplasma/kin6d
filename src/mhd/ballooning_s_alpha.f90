!> Reduced, local, high-n ideal-MHD ballooning calculation in shifted circles.
!> Local magnetic shear s and pressure-gradient drive alpha are dimensionless.
!> A finite angle horizon detects an instability, but absence of a zero
!> is NOT a certificate of continuum stability.
module milca_ballooning_s_alpha
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    integer, parameter, public :: BALLOON_OK = 0, BALLOON_INVALID = 1, &
        BALLOON_NOT_BRACKETED = 2, BALLOON_UNRESOLVED = 3
    public :: ballooning_first_zero, ballooning_critical_alpha, &
        ballooning_beta_limit, milca_ballooning_beta_abi
contains
    !> Integrate the self-adjoint marginal ballooning equation in flux form:
    !> (A F')' + V F = 0, A=1+(s*theta-alpha*sin(theta))**2,
    !> V=alpha*(cos(theta)+(s*theta-alpha*sin(theta))*sin(theta)).
    !> Even solution F(0)=1, F'(0)=0; first zero signals ideal instability.
    pure subroutine ballooning_first_zero(s, alpha, theta_max, nstep, crossed, status)
        real(dp), intent(in) :: s, alpha, theta_max
        integer, intent(in) :: nstep
        logical, intent(out) :: crossed
        integer, intent(out) :: status
        real(dp) :: theta, h, y(2), v(2), k1(2), k2(2), k3(2), k4(2)
        integer :: i
        status = BALLOON_INVALID
        crossed = .false.
        if (.not. all(ieee_is_finite([s, alpha, theta_max]))) return
        if (s <= 0 .or. alpha < 0 .or. theta_max <= 0 .or. nstep < 16) return
        h = theta_max/real(nstep,dp)
        y = [1.0_dp, 0.0_dp]
        theta = 0
        do i = 1, nstep
            call rhs(theta, y, s, alpha, k1)
            call rhs(theta+h/2, y+h*k1/2, s, alpha, k2)
            call rhs(theta+h/2, y+h*k2/2, s, alpha, k3)
            call rhs(theta+h, y+h*k3, s, alpha, k4)
            v = y + (h/6)*(k1+2*k2+2*k3+k4)
            if (.not. all(ieee_is_finite(v))) return
            if (v(1) <= 0) then
                crossed = .true.
                exit
            endif
            y = v
            theta = theta+h
        enddo
        status = BALLOON_OK
    end subroutine

    pure subroutine rhs(theta, y, s, alpha, dy)
        real(dp), intent(in) :: theta, y(2), s, alpha
        real(dp), intent(out) :: dy(2)
        real(dp) :: sn, cs, lambda, a
        sn = sin(theta)
        cs = cos(theta)
        lambda = s*theta-alpha*sn
        a = 1 + lambda*lambda
        dy(1) = y(2)/a
        dy(2) = -alpha*(cs+lambda*sn)*y(1)
    end subroutine

    !> First loss of stability as alpha is raised from zero. Uses Newcomb's
    !> zero-crossing test; later second-stability intervals are not searched.
    !> The result is a finite-domain, finite-step approximation, not a bound.
    pure subroutine ballooning_critical_alpha(s, alpha_max, theta_max, &
            nstep, alpha_crit, status)
        real(dp), intent(in) :: s, alpha_max, theta_max
        integer, intent(in) :: nstep
        real(dp), intent(out) :: alpha_crit
        integer, intent(out) :: status
        real(dp) :: lo, hi, mid, interval
        integer :: i, j, ns, rc
        logical :: crossed
        alpha_crit = 0
        status = BALLOON_INVALID
        if (.not. all(ieee_is_finite([s, alpha_max, theta_max]))) return
        if (s <= 0 .or. alpha_max <= 0 .or. alpha_max > 10 .or. &
            theta_max < 2 .or. nstep < 16) return
        ns = max(80, ceiling(alpha_max/0.025_dp))
        interval = alpha_max/real(ns,dp)
        lo = 0
        do i = 1, ns
            hi = real(i,dp)*interval
            call ballooning_first_zero(s,hi,theta_max,nstep,crossed,rc)
            if (rc /= BALLOON_OK) return
            if (crossed) exit
            lo = hi
        enddo
        if (.not. crossed) then
            status = BALLOON_NOT_BRACKETED
            return
        endif
        do j = 1, 26
            mid = (lo+hi)/2
            call ballooning_first_zero(s,mid,theta_max,nstep,crossed,rc)
            if (rc /= BALLOON_OK) return
            if (crossed) then
                hi = mid
            else
                lo = mid
            endif
        enddo
        alpha_crit = (lo+hi)/2
        status = BALLOON_OK
    end subroutine

    !> Peak and cross-section-averaged toroidal beta for p(r)=p0*(1-r^2/a^2).
    !> alpha=-2*mu0*R*q^2*p'/B_T^2 = 4*R*q^2*beta_avg*r/a^2.
    !> This checks LOCAL high-n ballooning only, NOT global kink/no-wall beta_N.
    pure subroutine ballooning_beta_limit(s, q, r, a, major_r, &
            beta_avg, error_estimate, status)
        real(dp), intent(in) :: s, q, r, a, major_r
        real(dp), intent(out) :: beta_avg, error_estimate
        integer, intent(out) :: status
        real(dp) :: alpha1, alpha2, alpha3, h1, h2, fac
        integer :: rc
        beta_avg = 0; error_estimate = 0
        status = BALLOON_INVALID
        if (.not. all(ieee_is_finite([s,q,r,a,major_r]))) return
        if (s <= 0 .or. q <= 0 .or. a <= 0 .or. major_r <= a .or. &
            r <= 0 .or. r >= a) return
        h1 = 16*acos(-1.0_dp); h2 = 2*h1
        call ballooning_critical_alpha(s, 4.0_dp, h1, 1024, alpha1, rc)
        if (rc /= BALLOON_OK) then
            status = rc
            return
        endif
        call ballooning_critical_alpha(s, 4.0_dp, h1, 2048, alpha2, rc)
        if (rc /= BALLOON_OK) then
            status = rc
            return
        endif
        call ballooning_critical_alpha(s, 4.0_dp, h2, 2048, alpha3, rc)
        if (rc /= BALLOON_OK) then
            status = rc
            return
        endif
        fac = a*a/(4*major_r*q*q*r)
        beta_avg = fac*alpha3
        ! Empirical finite-domain/discretization difference; not rigorous.
        error_estimate = fac*max(abs(alpha3-alpha2),abs(alpha2-alpha1))
        status = BALLOON_OK
        if (error_estimate > 0.02_dp*beta_avg) status = BALLOON_UNRESOLVED
    end subroutine

    !> Simple C ABI for a WASM/JS or embedded teaching frontend.
    subroutine milca_ballooning_beta_abi(s,q,r,a,major_r,beta,err,status) bind(C)
        real(c_double), value, intent(in) :: s,q,r,a,major_r
        real(c_double), intent(out) :: beta,err
        integer(c_int), intent(out) :: status
        integer :: rc
        call ballooning_beta_limit(s,q,r,a,major_r,beta,err,rc)
        status = int(rc,c_int)
    end subroutine
end module
