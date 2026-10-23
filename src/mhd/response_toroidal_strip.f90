module kin6d_response_toroidal_strip
    !! Candidate full-pressure, non-rotating strip action in gauge A_phi = 0.
    !! Per-n convention exp(i n phi); physical components at phi = 0 are
    !! Cartesian (e_R, e_phi, e_Z). The caller supplies a regular mapped
    !! annulus and equilibrium jets, including the Cartesian pressure Hessian.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_polynomial, only: lagrange_deriv_weights
    use fortnum_krylov, only: complex_gmres_operator
    use fortfem_direct_fourier_transform, only: &
        direct_fourier_plan_t, initialize_direct_fourier_plan, &
        direct_fourier_forward, direct_fourier_adjoint
    use fortsparse, only: fortsparse_status_t
    use kin6d_mhd_ampere_n, only: mu0
    use kin6d_response_fluid, only: &
        response_fluid_equilibrium_t, response_fluid_jet_t, &
        response_fluid_result_t, response_fluid_residual, response_fluid_adjoint, &
        response_fluid_stress, response_fluid_stress_adjoint
    use fortnum_krylov, only: complex_linear_operator_t, &
        complex_gmres, complex_fgmres
    implicit none
    private
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)
    type, public :: toroidal_strip_t
        private
        integer :: nr = 0, nt = 0, n = 0
        real(dp), allocatable :: radius(:, :), inverse(:, :, :, :), normal(:, :, :)
        real(dp), allocatable :: dr(:, :)
        complex(dp), allocatable :: dt(:, :)
        type(response_fluid_equilibrium_t), allocatable :: eq(:, :)
    end type toroidal_strip_t
    type, public :: toroidal_strip_work_t
        private
        complex(dp), allocatable :: vector(:, :, :), gradient(:, :, :, :)
        complex(dp), allocatable :: magnetic(:, :, :), current(:, :, :)
        complex(dp), allocatable :: divergence(:, :), gd(:, :, :), gp(:, :, :)
    end type toroidal_strip_work_t
    type, public :: toroidal_strip_status_t
        integer :: info = -1, iterations = 0
        real(dp) :: residual = huge(1.0_dp)
        real(dp) :: div_b = huge(1.0_dp), div_j = huge(1.0_dp)
        logical :: certified = .false.
    end type toroidal_strip_status_t
    public :: toroidal_strip_init, toroidal_strip_action, toroidal_strip_solve
    public :: toroidal_strip_vjp
    public :: toroidal_strip_defects
    public :: toroidal_strip_traces, toroidal_strip_traces_vjp
    public :: toroidal_strip_domains
    type, extends(complex_linear_operator_t) :: toroidal_strip_solve_matvec_t
        complex(dp), pointer :: gamma => null()
        type(toroidal_strip_t), pointer :: s => null()
        real(dp), pointer :: scales(:) => null()
        type(toroidal_strip_work_t), pointer :: work => null()
    contains
        procedure :: apply => toroidal_strip_solve_matvec_apply
    end type
contains
    subroutine toroidal_strip_domains(centers, widths, lo, hi, bounds, owner, info)
        !! Connected unions of rational strips in the caller's radial coordinate.
        !! Clip at the physical domain, merge overlaps/touching intervals and
        !! retain every center, including a rational exactly on either cut.
        !! Topology is discrete: no derivative across a merge/split is asserted.
        real(dp), intent(in) :: centers(:), widths(:), lo, hi
        real(dp), allocatable, intent(out) :: bounds(:, :)
        integer, allocatable, intent(out) :: owner(:)
        integer, intent(out) :: info
        real(dp), allocatable :: intervals(:, :)
        integer, allocatable :: order(:)
        real(dp) :: left, right
        integer :: i, j, k, groups, count
        info = 1
        count = size(centers)
        if (size(widths) /= count) return
        if (.not. ieee_is_finite(lo) .or. .not. ieee_is_finite(hi)) return
        if (lo >= hi) return
        if (any(.not. ieee_is_finite(centers))) return
        if (any(.not. ieee_is_finite(widths))) return
        if (any(widths <= 0.0_dp)) return
        if (any(centers < lo) .or. any(centers > hi)) return
        allocate (intervals(2, count), order(count), owner(count))
        do i = 1, count
            order(i) = i
        end do
        do i = 2, count
            k = order(i); j = i - 1
            do while (j >= 1)
                if (max(lo, centers(order(j)) - widths(order(j))) &
                    <= max(lo, centers(k) - widths(k))) exit
                order(j + 1) = order(j)
                j = j - 1
            end do
            order(j + 1) = k
        end do
        groups = 0
        do i = 1, count
            k = order(i)
            left = max(lo, centers(k) - widths(k))
            right = min(hi, centers(k) + widths(k))
            if (groups == 0) then
                groups = 1
                intervals(1, groups) = left; intervals(2, groups) = right
            else if (left <= intervals(2, groups)) then
                intervals(2, groups) = max(intervals(2, groups), right)
                intervals(1, groups) = min(intervals(1, groups), left)
            else
                groups = groups + 1
                intervals(1, groups) = left; intervals(2, groups) = right
            end if
            owner(k) = groups
        end do
        bounds = intervals(:, 1:groups)
        info = 0
    end subroutine toroidal_strip_domains
    subroutine toroidal_strip_traces(s, work, gamma, x, traces)
        !! Both rings, in one orientation N = grad(u)/|grad(u)|, T=(-N_Z,0,N_R).
        !! traces(:,theta,side) = (A.T,A_phi,B.N,E.T,E_phi,xi.N,traction(3),Ptot1).
        !! E=-gamma A-grad Phi; traction is Eulerian Maxwell-plus-fluid stress.
        !! A matching outer trace uses the same N, so continuity is subtraction;
        !! outward tractions instead add. No sheet or applied traction assumed.
        !! Fixed equilibrium JVP is this linear map. Bounds are unavailable.
        type(toroidal_strip_t), intent(in) :: s
        type(toroidal_strip_work_t), intent(inout) :: work
        complex(dp), intent(in) :: gamma, x(6, s%nr, s%nt)
        complex(dp), intent(out) :: traces(10, s%nt, 2)
        type(response_fluid_jet_t) :: jet
        type(response_fluid_result_t) :: res
        complex(dp) :: electric(3), stress(4)
        real(dp) :: normal(3), tangent(3)
        integer :: side, i, j
        work%vector(1, :, :) = x(1, :, :)
        work%vector(2, :, :) = 0.0_dp
        work%vector(3, :, :) = x(2, :, :)
        call vector_gradient(s, work%vector, work%gradient)
        call curl_from_gradient(s, work%gradient, work%magnetic)
        work%vector = x(3:5, :, :)
        call vector_gradient(s, work%vector, work%gradient)
        call scalar_gradient(s, x(6, :, :), work%gp)
        !$omp parallel do collapse(2) schedule(static) &
        !$omp private(i,jet,res,electric,stress,normal,tangent)
        do side = 1, 2
            do j = 1, s%nt
                i = 1 + (side - 1)*(s%nr - 1)
                normal(1) = s%normal(1, i, j); normal(2) = 0.0_dp
                normal(3) = s%normal(2, i, j)
                tangent(1) = -normal(3); tangent(2) = 0.0_dp
                tangent(3) = normal(1)
                jet = response_fluid_jet_t()
                jet%a(1) = x(1, i, j); jet%a(2) = 0.0_dp
                jet%a(3) = x(2, i, j)
                jet%xi = x(3:5, i, j)
                jet%dxi = work%gradient(:, :, i, j)
                call response_fluid_residual(s%eq(i, j), jet, gamma, res)
                electric = -gamma*jet%a - work%gp(:, i, j)
                call response_fluid_stress(s%eq(i, j)%b0, normal, res%p1, &
                    work%magnetic(:, i, j), stress)
                traces(1, j, side) = sum(tangent*jet%a)
                traces(2, j, side) = jet%a(2)
                traces(3, j, side) = sum(normal*work%magnetic(:, i, j))
                traces(4, j, side) = sum(tangent*electric)
                traces(5, j, side) = electric(2)
                traces(6, j, side) = sum(normal*jet%xi)
                traces(7:10, j, side) = stress
            end do
        end do
        !$omp end parallel do
    end subroutine toroidal_strip_traces
    subroutine toroidal_strip_traces_vjp(s, work, gamma, seed, adj)
        !! Hermitian storage transpose. Physical surface weights multiply seed
        !! before this call; volume Riesz weights are applied by the consumer.
        type(toroidal_strip_t), intent(in) :: s
        type(toroidal_strip_work_t), intent(inout) :: work
        complex(dp), intent(in) :: gamma, seed(10, s%nt, 2)
        complex(dp), intent(out) :: adj(6, s%nr, s%nt)
        type(response_fluid_result_t) :: res
        type(response_fluid_jet_t) :: jet
        complex(dp) :: electric(3), magnetic(3)
        real(dp) :: normal(3), tangent(3)
        integer :: side, i, j
        adj = 0.0_dp; work%gradient = 0.0_dp; work%gp = 0.0_dp
        work%magnetic = 0.0_dp
        !$omp parallel do collapse(2) schedule(static) &
        !$omp private(i,res,jet,electric,magnetic,normal,tangent)
        do side = 1, 2
            do j = 1, s%nt
                i = 1 + (side - 1)*(s%nr - 1)
                normal(1) = s%normal(1, i, j); normal(2) = 0.0_dp
                normal(3) = s%normal(2, i, j)
                tangent(1) = -normal(3); tangent(2) = 0.0_dp
                tangent(3) = normal(1)
                res = response_fluid_result_t()
                call response_fluid_stress_adjoint(s%eq(i, j)%b0, normal, &
                    seed(7:10, j, side), res%p1, magnetic)
                call response_fluid_adjoint(s%eq(i, j), gamma, res, jet)
                electric = tangent*seed(4, j, side)
                electric(2) = seed(5, j, side)
                jet%a = tangent*seed(1, j, side) - conjg(gamma)*electric
                jet%a(2) = jet%a(2) + seed(2, j, side)
                adj(1, i, j) = jet%a(1); adj(2, i, j) = jet%a(3)
                adj(3:5, i, j) = jet%xi + normal*seed(6, j, side)
                work%gradient(:, :, i, j) = jet%dxi
                work%magnetic(:, i, j) = magnetic + normal*seed(3, j, side)
                work%gp(:, i, j) = -electric
            end do
        end do
        !$omp end parallel do
        call scalar_gradient_vjp(s, work%gp, adj(6, :, :))
        call vector_gradient_vjp(s, work%gradient, work%vector)
        adj(3:5, :, :) = adj(3:5, :, :) + work%vector
        call curl_from_gradient_vjp(s, work%magnetic, work%gradient)
        call vector_gradient_vjp(s, work%gradient, work%vector)
        adj(1, :, :) = adj(1, :, :) + work%vector(1, :, :)
        adj(2, :, :) = adj(2, :, :) + work%vector(3, :, :)
    end subroutine toroidal_strip_traces_vjp
    subroutine toroidal_strip_init(s, work, radial, radius, inverse, eq, n, info)
        !! inverse(:, :, i, j) = d(u, theta)/d(R, Z); theta_j = 2 pi(j-1)/nt.
        !! Local five-node radial polynomial derivatives and complete odd-size
        !! Fourier collocation (no omitted Nyquist mode). No rational denominator.
        type(toroidal_strip_t), intent(out) :: s
        type(toroidal_strip_work_t), intent(out) :: work
        real(dp), intent(in) :: radial(:), radius(:, :), inverse(:, :, :, :)
        type(response_fluid_equilibrium_t), intent(in) :: eq(:, :)
        integer, intent(in) :: n
        integer, intent(out) :: info
        type(direct_fourier_plan_t) :: plan
        type(fortsparse_status_t) :: st
        real(dp), allocatable :: angles(:)
        complex(dp), allocatable :: v(:), c(:), dv(:)
        integer, allocatable :: modes(:)
        real(dp) :: weights(5), length
        integer :: nr, nt, i, j, lo
        info = -1
        nr = size(radial); nt = size(radius, 2)
        if (nr < 5 .or. nt < 5 .or. mod(nt, 2) /= 1 .or. n == 0) return
        if (size(radius, 1) /= nr) return
        if (any(shape(inverse) /= [2, 2, nr, nt])) return
        if (any(shape(eq) /= [nr, nt])) return
        if (any(.not. ieee_is_finite(radial))) return
        if (any(radial(2:) <= radial(:nr - 1))) return
        if (any(.not. ieee_is_finite(radius)) .or. any(radius <= 0.0_dp)) return
        if (any(.not. ieee_is_finite(inverse))) return
        do j = 1, nt
            do i = 1, nr
                length = norm2(inverse(1, :, i, j))
                if (length <= tiny(1.0_dp)) return
                if (abs(inverse(1, 1, i, j)*inverse(2, 2, i, j) &
                        - inverse(1, 2, i, j)*inverse(2, 1, i, j)) &
                    <= tiny(1.0_dp)) return
                if (any(.not. ieee_is_finite([eq(i, j)%rho, eq(i, j)%eta, &
                                              eq(i, j)%p, eq(i, j)%adiabatic]))) return
                if (any(.not. ieee_is_finite(eq(i, j)%b0))) return
                if (any(.not. ieee_is_finite(eq(i, j)%j0))) return
                if (any(.not. ieee_is_finite(eq(i, j)%gp))) return
                if (any(.not. ieee_is_finite(eq(i, j)%hp))) return
                if (eq(i, j)%rho <= 0.0_dp .or. eq(i, j)%eta < 0.0_dp &
                    .or. eq(i, j)%p < 0.0_dp .or. eq(i, j)%adiabatic <= 0.0_dp) return
            end do
        end do
        s%nr = nr; s%nt = nt; s%n = n
        s%radius = radius; s%inverse = inverse; s%eq = eq
        allocate (s%normal(2, nr, nt), s%dr(nr, nr), s%dt(nt, nt))
        do j = 1, nt
            do i = 1, nr
                s%normal(:, i, j) = inverse(1, :, i, j)/norm2(inverse(1, :, i, j))
            end do
        end do
        s%dr = 0.0_dp
        do i = 1, nr
            lo = min(max(i - 2, 1), nr - 4)
            call lagrange_deriv_weights(5, radial(i), radial(lo:lo + 4), weights)
            s%dr(i, lo:lo + 4) = weights
        end do
        allocate (angles(nt), modes(nt), v(nt), c(nt), dv(nt))
        do j = 1, nt
            angles(j) = 2*acos(-1.0_dp)*(j - 1)/nt
            modes(j) = j - 1 - (nt - 1)/2
        end do
        call initialize_direct_fourier_plan(plan, angles, modes, nt, status=st)
        if (st%code /= 0) return
        do j = 1, nt
            v = 0.0_dp; v(j) = 1.0_dp
            call direct_fourier_forward(plan, v, c, st)
            if (st%code /= 0) return
            c = iu*modes*c
            call direct_fourier_adjoint(plan, c, dv, st)
            if (st%code /= 0) return
            s%dt(:, j) = dv
        end do
        allocate (work%vector(3, nr, nt), work%gradient(3, 3, nr, nt), &
                  work%magnetic(3, nr, nt), work%current(3, nr, nt), &
                  work%divergence(nr, nt), work%gd(3, nr, nt), work%gp(3, nr, nt))
        info = 0
    end subroutine toroidal_strip_init
    subroutine scalar_gradient(s, v, g)
        type(toroidal_strip_t), intent(in) :: s
        complex(dp), intent(in) :: v(:, :)
        complex(dp), intent(out) :: g(:, :, :)
        complex(dp) :: du, dt
        integer :: i, j, k
        !$omp parallel do collapse(2) private(du,dt,k) schedule(static)
        do j = 1, s%nt
            do i = 1, s%nr
                du = 0.0_dp; dt = 0.0_dp
                do k = max(1, min(i - 2, s%nr - 4)), &
                    min(s%nr, max(i + 2, 5))
                    du = du + s%dr(i, k)*v(k, j)
                end do
                do k = 1, s%nt
                    dt = dt + s%dt(j, k)*v(i, k)
                end do
                g(1, i, j) = s%inverse(1, 1, i, j)*du &
                             + s%inverse(2, 1, i, j)*dt
                g(2, i, j) = iu*s%n*v(i, j)/s%radius(i, j)
                g(3, i, j) = s%inverse(1, 2, i, j)*du &
                             + s%inverse(2, 2, i, j)*dt
            end do
        end do
        !$omp end parallel do
    end subroutine scalar_gradient
    subroutine vector_gradient(s, v, g)
        type(toroidal_strip_t), intent(in) :: s
        complex(dp), intent(in) :: v(3, s%nr, s%nt)
        complex(dp), intent(out) :: g(3, 3, s%nr, s%nt)
        complex(dp) :: du, dt
        integer :: i, j, k, a
        !$omp parallel do collapse(2) private(du,dt,k,a) schedule(static)
        do j = 1, s%nt
            do i = 1, s%nr
                do a = 1, 3
                    du = 0.0_dp; dt = 0.0_dp
                    do k = max(1, min(i - 2, s%nr - 4)), &
                        min(s%nr, max(i + 2, 5))
                        du = du + s%dr(i, k)*v(a, k, j)
                    end do
                    do k = 1, s%nt
                        dt = dt + s%dt(j, k)*v(a, i, k)
                    end do
                    g(a, 1, i, j) = s%inverse(1, 1, i, j)*du &
                                    + s%inverse(2, 1, i, j)*dt
                    g(a, 2, i, j) = iu*s%n*v(a, i, j)/s%radius(i, j)
                    g(a, 3, i, j) = s%inverse(1, 2, i, j)*du &
                                    + s%inverse(2, 2, i, j)*dt
                end do
                ! Cartesian derivative of the rotating cylindrical basis.
                g(1, 2, i, j) = g(1, 2, i, j) - v(2, i, j)/s%radius(i, j)
                g(2, 2, i, j) = g(2, 2, i, j) + v(1, i, j)/s%radius(i, j)
            end do
        end do
        !$omp end parallel do
    end subroutine vector_gradient
    subroutine curl_from_gradient(s, g, v)
        type(toroidal_strip_t), intent(in) :: s
        complex(dp), intent(in) :: g(3, 3, s%nr, s%nt)
        complex(dp), intent(out) :: v(3, s%nr, s%nt)
        integer :: i, j
        !$omp parallel do collapse(2) schedule(static)
        do j = 1, s%nt
            do i = 1, s%nr
                v(1, i, j) = g(3, 2, i, j) - g(2, 3, i, j)
                v(2, i, j) = g(1, 3, i, j) - g(3, 1, i, j)
                v(3, i, j) = g(2, 1, i, j) - g(1, 2, i, j)
            end do
        end do
        !$omp end parallel do
    end subroutine curl_from_gradient
    subroutine toroidal_strip_action(s, work, gamma, x, y, boundary)
        !! x = (A_R,A_Z,xi_R,xi_phi,xi_Z,Phi), y = (Ohm_R,Ohm_Z,
        !! momentum_R,momentum_phi,momentum_Z,Ohm_phi).
        !! boundary replaces only A_t and xi_normal on the two rings;
        !! the other four equations remain active there. These are prescribed
        !! traces for the isolated strip, not an outer-response interface.
        type(toroidal_strip_t), intent(in) :: s
        type(toroidal_strip_work_t), intent(inout) :: work
        complex(dp), intent(in) :: gamma, x(6, s%nr, s%nt)
        complex(dp), intent(out) :: y(6, s%nr, s%nt)
        logical, intent(in) :: boundary
        type(response_fluid_jet_t) :: jet
        type(response_fluid_result_t) :: res
        real(dp) :: nr, nz
        integer :: i, j
        work%vector(1, :, :) = x(1, :, :)
        work%vector(2, :, :) = 0.0_dp
        work%vector(3, :, :) = x(2, :, :)
        call vector_gradient(s, work%vector, work%gradient)
        call curl_from_gradient(s, work%gradient, work%magnetic)
        call vector_gradient(s, work%magnetic, work%gradient)
        call curl_from_gradient(s, work%gradient, work%current)
        work%current = work%current/mu0
        work%vector = x(3:5, :, :)
        call vector_gradient(s, work%vector, work%gradient)
        work%divergence = work%gradient(1, 1, :, :) &
                          + work%gradient(2, 2, :, :) + work%gradient(3, 3, :, :)
        call scalar_gradient(s, work%divergence, work%gd)
        call scalar_gradient(s, x(6, :, :), work%gp)
        !$omp parallel do collapse(2) private(jet,res,nr,nz) schedule(static)
        do j = 1, s%nt
            do i = 1, s%nr
                jet%a = [x(1, i, j), (0.0_dp, 0.0_dp), x(2, i, j)]
                jet%xi = x(3:5, i, j)
                jet%dxi = work%gradient(:, :, i, j)
                jet%gd = work%gd(:, i, j)
                jet%grad_phi = work%gp(:, i, j)
                jet%b = work%magnetic(:, i, j)
                jet%j = work%current(:, i, j)
                call response_fluid_residual(s%eq(i, j), jet, gamma, res)
                y(:, i, j) = [res%ohm(1), res%ohm(3), res%momentum, res%ohm(2)]
                if (boundary .and. (i == 1 .or. i == s%nr)) then
                    nr = s%normal(1, i, j); nz = s%normal(2, i, j)
                    y(1, i, j) = nr*res%ohm(1) + nz*res%ohm(3)
                    y(2, i, j) = -nz*x(1, i, j) + nr*x(2, i, j)
                    y(3, i, j) = nr*x(3, i, j) + nz*x(5, i, j)
                    y(5, i, j) = -nz*res%momentum(1) + nr*res%momentum(3)
                end if
            end do
        end do
        !$omp end parallel do
    end subroutine toroidal_strip_action
    subroutine scalar_gradient_vjp(s, g, v)
        type(toroidal_strip_t), intent(in) :: s
        complex(dp), intent(in) :: g(:, :, :)
        complex(dp), intent(out) :: v(:, :)
        complex(dp) :: du, dt, total
        integer :: i, j, k
        !$omp parallel do collapse(2) private(du,dt,total,k) schedule(static)
        do j = 1, s%nt
            do i = 1, s%nr
                total = -iu*s%n*g(2, i, j)/s%radius(i, j)
                do k = max(1, i - 4), min(s%nr, i + 4)
                    du = s%inverse(1, 1, k, j)*g(1, k, j) &
                         + s%inverse(1, 2, k, j)*g(3, k, j)
                    total = total + s%dr(k, i)*du
                end do
                do k = 1, s%nt
                    dt = s%inverse(2, 1, i, k)*g(1, i, k) &
                         + s%inverse(2, 2, i, k)*g(3, i, k)
                    total = total + conjg(s%dt(k, j))*dt
                end do
                v(i, j) = total
            end do
        end do
        !$omp end parallel do
    end subroutine scalar_gradient_vjp
    subroutine vector_gradient_vjp(s, g, v)
        type(toroidal_strip_t), intent(in) :: s
        complex(dp), intent(in) :: g(3, 3, s%nr, s%nt)
        complex(dp), intent(out) :: v(3, s%nr, s%nt)
        complex(dp) :: du, dt, total
        integer :: i, j, k, a
        !$omp parallel do collapse(2) private(du,dt,total,k,a) schedule(static)
        do j = 1, s%nt
            do i = 1, s%nr
                do a = 1, 3
                    total = -iu*s%n*g(a, 2, i, j)/s%radius(i, j)
                    do k = max(1, i - 4), min(s%nr, i + 4)
                        du = s%inverse(1, 1, k, j)*g(a, 1, k, j) &
                             + s%inverse(1, 2, k, j)*g(a, 3, k, j)
                        total = total + s%dr(k, i)*du
                    end do
                    do k = 1, s%nt
                        dt = s%inverse(2, 1, i, k)*g(a, 1, i, k) &
                             + s%inverse(2, 2, i, k)*g(a, 3, i, k)
                        total = total + conjg(s%dt(k, j))*dt
                    end do
                    v(a, i, j) = total
                end do
                v(1, i, j) = v(1, i, j) + g(2, 2, i, j)/s%radius(i, j)
                v(2, i, j) = v(2, i, j) - g(1, 2, i, j)/s%radius(i, j)
            end do
        end do
        !$omp end parallel do
    end subroutine vector_gradient_vjp
    subroutine curl_from_gradient_vjp(s, v, g)
        type(toroidal_strip_t), intent(in) :: s
        complex(dp), intent(in) :: v(3, s%nr, s%nt)
        complex(dp), intent(out) :: g(3, 3, s%nr, s%nt)
        integer :: i, j
        !$omp parallel do collapse(2) schedule(static)
        do j = 1, s%nt
            do i = 1, s%nr
                g(:, :, i, j) = (0.0_dp, 0.0_dp)
                g(3, 2, i, j) = v(1, i, j); g(2, 3, i, j) = -v(1, i, j)
                g(1, 3, i, j) = v(2, i, j); g(3, 1, i, j) = -v(2, i, j)
                g(2, 1, i, j) = v(3, i, j); g(1, 2, i, j) = -v(3, i, j)
            end do
        end do
        !$omp end parallel do
    end subroutine curl_from_gradient_vjp
    subroutine toroidal_strip_vjp(s, work, gamma, seed, adj, boundary)
        !! Complex cotangent transpose for the stored-grid Hermitian pairing.
        !! A physical weighted adjoint is M^-1 L^H M with the caller's volume
        !! quadrature M; this storage VJP does not assert self-adjointness.
        !! Equilibrium, chart, mesh and gamma are fixed, as in the primal JVP.
        type(toroidal_strip_t), intent(in) :: s
        type(toroidal_strip_work_t), intent(inout) :: work
        complex(dp), intent(in) :: gamma, seed(6, s%nr, s%nt)
        complex(dp), intent(out) :: adj(6, s%nr, s%nt)
        logical, intent(in) :: boundary
        type(response_fluid_result_t) :: res
        type(response_fluid_jet_t) :: jet
        real(dp) :: nr, nz
        integer :: i, j, a
        !$omp parallel do collapse(2) private(jet,res,nr,nz) schedule(static)
        do j = 1, s%nt
            do i = 1, s%nr
                res = response_fluid_result_t()
                res%ohm = [seed(1, i, j), seed(6, i, j), seed(2, i, j)]
                res%momentum = seed(3:5, i, j)
                adj(:, i, j) = 0.0_dp
                if (boundary .and. (i == 1 .or. i == s%nr)) then
                    nr = s%normal(1, i, j); nz = s%normal(2, i, j)
                    res%ohm(1) = nr*seed(1, i, j)
                    res%ohm(3) = nz*seed(1, i, j)
                    res%momentum(1) = -nz*seed(5, i, j)
                    res%momentum(3) = nr*seed(5, i, j)
                    adj(1, i, j) = -nz*seed(2, i, j)
                    adj(2, i, j) = nr*seed(2, i, j)
                    adj(3, i, j) = nr*seed(3, i, j)
                    adj(5, i, j) = nz*seed(3, i, j)
                end if
                call response_fluid_adjoint(s%eq(i, j), gamma, res, jet)
                adj(1, i, j) = adj(1, i, j) + jet%a(1)
                adj(2, i, j) = adj(2, i, j) + jet%a(3)
                adj(3:5, i, j) = adj(3:5, i, j) + jet%xi
                work%gradient(:, :, i, j) = jet%dxi
                work%gd(:, i, j) = jet%gd
                work%gp(:, i, j) = jet%grad_phi
                work%current(:, i, j) = jet%j/mu0
                work%magnetic(:, i, j) = jet%b
            end do
        end do
        !$omp end parallel do
        call scalar_gradient_vjp(s, work%gp, adj(6, :, :))
        call scalar_gradient_vjp(s, work%gd, work%divergence)
        do a = 1, 3
            work%gradient(a, a, :, :) = work%gradient(a, a, :, :) + work%divergence
        end do
        call vector_gradient_vjp(s, work%gradient, work%vector)
        adj(3:5, :, :) = adj(3:5, :, :) + work%vector
        call curl_from_gradient_vjp(s, work%current, work%gradient)
        call vector_gradient_vjp(s, work%gradient, work%vector)
        work%vector = work%vector + work%magnetic
        call curl_from_gradient_vjp(s, work%vector, work%gradient)
        call vector_gradient_vjp(s, work%gradient, work%magnetic)
        adj(1, :, :) = adj(1, :, :) + work%magnetic(1, :, :)
        adj(2, :, :) = adj(2, :, :) + work%magnetic(3, :, :)
    end subroutine toroidal_strip_vjp
    subroutine toroidal_strip_defects(s, work, div_b, div_j)
        !! Maximum nodal div(curl A) [T/m] and div(J) [A/m^3] of the last
        !! action. Collocation is not a commuting FEEC complex; these are
        !! structure-defect indicators, not continuous bounds or certificates.
        type(toroidal_strip_t), intent(in) :: s
        type(toroidal_strip_work_t), intent(inout) :: work
        real(dp), intent(out) :: div_b, div_j
        integer :: i, j
        div_b = 0.0_dp; div_j = 0.0_dp
        call vector_gradient(s, work%magnetic, work%gradient)
        !$omp parallel do collapse(2) reduction(max:div_b) schedule(static)
        do j = 1, s%nt
            do i = 1, s%nr
                div_b = max(div_b, abs(work%gradient(1, 1, i, j) &
                                       + work%gradient(2, 2, i, j) &
                                       + work%gradient(3, 3, i, j)))
            end do
        end do
        !$omp end parallel do
        call vector_gradient(s, work%current, work%gradient)
        !$omp parallel do collapse(2) reduction(max:div_j) schedule(static)
        do j = 1, s%nt
            do i = 1, s%nr
                div_j = max(div_j, abs(work%gradient(1, 1, i, j) &
                                       + work%gradient(2, 2, i, j) &
                                       + work%gradient(3, 3, i, j)))
            end do
        end do
        !$omp end parallel do
    end subroutine toroidal_strip_defects
    subroutine toroidal_strip_solve(s, work, gamma, rhs, x, equation_scale, &
                                    tol, maxit, restart, st)
        !! Same FortNum GMRES as the outer strip-coupling consumer. A tangent
        !! for a supplied RHS perturbation solves this same discrete equation.
        !! Residual is checked again independently of the Krylov recurrence.
        !! No inverse bound or discretisation certificate is available.
        type(toroidal_strip_t), intent(in), target :: s
        type(toroidal_strip_work_t), intent(inout), target :: work
        complex(dp), intent(in), target :: gamma, rhs(6, s%nr, s%nt)
        complex(dp), intent(inout) :: x(6, s%nr, s%nt)
        real(dp), intent(in) :: equation_scale(6, s%nr, s%nt), tol
        integer, intent(in) :: maxit, restart
        type(toroidal_strip_status_t), intent(out) :: st
        complex(dp), allocatable :: b(:), v(:), ax(:)
        real(dp), allocatable , target :: scales(:)
        integer :: nn
        type(toroidal_strip_solve_matvec_t) :: matvec_operator
        if (.not. ieee_is_finite(real(gamma, dp))) return
        if (.not. ieee_is_finite(aimag(gamma))) return
        if (abs(gamma) <= tiny(1.0_dp)) return
        if (any(.not. ieee_is_finite(equation_scale))) return
        if (any(equation_scale <= 0.0_dp)) return
        nn = 6*s%nr*s%nt
        allocate (b(nn), v(nn), ax(nn), scales(nn))
        scales = reshape(equation_scale, [nn])
        b = reshape(rhs, [nn])/scales; v = reshape(x, [nn])
        matvec_operator%gamma => gamma
        matvec_operator%s => s
        nullify(matvec_operator%scales)
        if (allocated(scales)) matvec_operator%scales => scales
        matvec_operator%work => work
        call complex_gmres(matvec_operator, b, v, tol, maxit, restart, &
                                    st%info, st%iterations, st%residual)
        call matvec_operator%apply(v, ax)
        st%residual = norm2(abs(ax - b))/max(norm2(abs(b)), 1.0_dp)
        if (st%info == 0 .and. st%residual > tol) st%info = 1
        call toroidal_strip_defects(s, work, st%div_b, st%div_j)
        x = reshape(v, shape(x))
    end subroutine toroidal_strip_solve
    subroutine toroidal_strip_solve_matvec_apply(self, input, output)
        class(toroidal_strip_solve_matvec_t), intent(inout) :: self
            complex(dp), intent(in) :: input(:)
            complex(dp), intent(out) :: output(:)
            ! Array sequence association preserves the canonical field layout;
            ! FortNum supplies contiguous Krylov columns, without reshape copies.
            call toroidal_strip_action(self%s, self%work, self%gamma, input, output, &
                .true.)
            output = output/self%scales
    end subroutine toroidal_strip_solve_matvec_apply
end module kin6d_response_toroidal_strip
