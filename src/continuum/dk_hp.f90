!> Monoenergetic drift-kinetic coefficients D11, D31, D13, D33 on an
!! axisymmetric surface, the equation of milca_dk_mono, by continuous
!! spectral elements on a mesh of the (theta, xi) cylinder that is fitted to
!! the orbit classes.
!!
!! lambda = (1 - xi^2)/b is constant along orbits, so the layers of the
!! solution lie on level curves of lambda, above all the trapped-passing
!! boundary lambda_c = 1/b_max, whose width shrinks like sqrt(nu). Element
!! rows are bounded by such curves: in the passing region by levels
!! xi = sigma sqrt(1 - Lambda b(theta)), in the trapped region by scaled
!! separatrices xi = zeta xi_c(theta), xi_c = sqrt(1 - b/b_max), so that
!! zeta = +-1 is the separatrix itself. Columns are theta intervals of one
!! period that starts and ends at the maximum theta_x of b, where the
!! separatrix closes in the X point. Trapped elements next to theta_x
!! collapse to that point; its nodes are one unknown, which also makes
!! the discrete functions regular there. Rows are graded geometrically
!! towards the separatrix from both sides and columns towards theta_x.
!! Reflection points need no treatment in (theta, xi): the solution is
!! smooth there, unlike in (theta, lambda).
!!
!! Each element carries tensor Lagrange polynomials of degree p on
!! Gauss-Lobatto nodes in the reference square, mapped to the element by
!! theta linear in the first and xi linear between the bounding curves in
!! the second coordinate. The weak form uses the measure (1/b^2) dtheta dxi/2
!! of the flux-surface and pitch averages, in which streaming is
!! divergence free; there is no inflow boundary, so no upwinding. The
!! density null space is bordered as in milca_dk_mono, and the system is
!! factorised once by a sparse direct solver; adjoint solves reuse the
!! factors through the transpose.
module milca_dk_hp
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre, gauss_lobatto_legendre
    use fortnum_polynomial, only: lagrange_weights, lagrange_deriv_weights
    use fortsparse, only: csc_t, csc_z_t, csc_from_triplet, sparse_solver_t, sparse_factor, &
        sparse_solve_op, sparse_destroy, fortsparse_status_t, status_ok, &
        FORTSPARSE_OP_NORMAL, FORTSPARSE_OP_TRANSPOSE, FORTSPARSE_BACKEND_NATIVE, &
        FORTSPARSE_LU_SUPERNODAL, FORTSPARSE_LU_SCALAR, FORTSPARSE_ORDER_AMD
    use milca_dk_hp_leaves, only: dk_hp_point, dk_hp_stream, dk_hp_ntv
    use milca_dk_mono, only: mono_coefficients_t
    implicit none
    private
    public :: hp_geometry_t, hp_field_t, hp_field, hp_mesh_t, hp_default_mesh, hp_system_t, &
        hp_assemble, hp_factor, hp_solve, hp_coefficients, hp_free, hp_estimate, &
        hp_refine, hp_operator_t, hp_ntv_kernel
    integer, parameter, public :: HP_OK = 0, HP_BAD_INPUT = 1, HP_SINGULAR = 2

    real(dp), parameter :: pi = 4*atan(1.0_dp)


    !> Operator per unit speed: b.grad theta = k b**k_power, nu = nu_D/v; for
    !! one toroidal harmonic n of the NTV problem (milca_dk_ntv) also q, the
    !! speed v, omega_E, the precession prefactor pre and offset, the drift
    !! weight sign wsign and the perturbation c_m exp(i m theta).
    type :: hp_operator_t
        real(dp) :: k = 0, nu = 0, q = 0, v = 1, omega_e = 0, pre = 0, offset = 0, &
            wsign = -1
        integer :: k_power = 0, n = 0
        integer, allocatable :: m(:)
        complex(dp), allocatable :: c(:)
        ! The two further NEO-2-QL sources of milca_dk_ntv (thin fields only):
        ! the parallel-derivative part of the perturbed radial drift and the
        ! perturbed operator acting on the axisymmetric solution f_10.
        logical :: source_dpar = .false., source_dlf10 = .false.
    end type

    !> The orbit geometry of the mesh: on the (theta, xi) cylinder of one
    !> surface, the stream function of the orbits (its bracket is the
    !> streaming term in the measure of the weak form), the row curves
    !> (passing levels of sign kind = +-1 at level val, scaled separatrices
    !> kind = 0 at val = zeta), and the pointwise coefficients of the weak form
    !> per unit dtheta dxi. The separatrix closes at theta_x at the passing
    !> level lambda_c.
    type, abstract :: hp_geometry_t
        real(dp) :: kappa = 0, theta_x = 0, lambda_c = 0
    contains
        procedure(stream_i), deferred :: stream
        procedure(level_i), deferred :: level
        procedure(point_i), deferred :: point
        procedure :: perturbation => op_perturbation
    end type

    abstract interface
        real(dp) function stream_i(g, th, xi)
            import :: dp, hp_geometry_t
            class(hp_geometry_t), intent(in) :: g
            real(dp), intent(in) :: th, xi
        end function stream_i

        subroutine level_i(g, kind, val, th, xi, dxi)
            import :: dp, hp_geometry_t
            class(hp_geometry_t), intent(in) :: g
            integer, intent(in) :: kind
            real(dp), intent(in) :: val, th
            real(dp), intent(out) :: xi, dxi
        end subroutine level_i

        !> collide: coefficient of d_xi w d_xi g; drive_r, drive_e: the radial
        !> and parallel drives of the axisymmetric problem; density: the
        !> measure; phase, precess, drive, torque: the NTV terms of
        !> milca_dk_hp_leaves (dk_hp_ntv).
        subroutine point_i(g, th, xi, op, collide, drive_r, drive_e, density, phase, &
                precess, drive, torque)
            import :: dp, hp_geometry_t, hp_operator_t
            class(hp_geometry_t), intent(in) :: g
            real(dp), intent(in) :: th, xi
            type(hp_operator_t), intent(in) :: op
            real(dp), intent(out) :: collide, drive_r, drive_e, density, phase, precess, &
                drive, torque
        end subroutine point_i
    end interface

    !> b(theta) = a_0 + sum_m a_m cos(m theta) + c_m sin(m theta), the
    !! trigonometric interpolant of the samples, with its single maximum
    !! b_max at theta_x and minimum b_min: the thin-orbit geometry, whose
    !! stream function is -lambda/2, lambda = (1 - xi^2)/b.
    type, extends(hp_geometry_t) :: hp_field_t
        integer :: m = 0
        real(dp) :: b_max = 0, b_min = 0
        real(dp), allocatable :: a(:), c(:)
    contains
        procedure :: stream => field_stream
        procedure :: level => curve
        procedure :: point => field_point
    end type

    !> Column breakpoints theta(1) = theta_x - 2 pi < ... < theta(end) =
    !! theta_x; passing levels lam(1) = 0 < ... < lam(end) = lambda_c;
    !! trapped rows zeta(1) = -1 < ... < zeta(end) = 1; degree p.
    type :: hp_mesh_t
        integer :: p = 0
        real(dp), allocatable :: theta(:), lam(:), zeta(:)
    end type


    !> Assembled system of n pitch unknowns with its element matrices
    !! (columns of mat), loads and density weights for residual estimates;
    !! loc maps local element nodes to unknowns.
    type :: hp_system_t
        integer :: n = 0
        real(dp) :: norm = 0 ! int dtheta dxi/(b b.grad theta)
        real(dp), allocatable :: density(:), mat(:, :), dens(:, :)
        complex(dp), allocatable :: drive(:, :), load(:, :, :)
        integer, allocatable :: loc(:, :)
        type(csc_t) :: a
        type(sparse_solver_t) :: solver
        logical :: factored = .false.
    end type

contains

    !> Field from b on n (odd) uniform points theta_j = 2 pi (j - 1)/n and
    !! kappa = b.grad theta / b; refuses fields with more than one maximum.
    subroutine hp_field(b, kappa, f, status)
        real(dp), intent(in) :: b(:), kappa
        type(hp_field_t), intent(out) :: f
        integer, intent(out) :: status
        integer :: n, k, j, n_fine, n_max, n_min
        real(dp) :: th, d0, d1
        real(dp), allocatable :: theta(:)

        status = HP_BAD_INPUT
        n = size(b)
        if (mod(n, 2) == 0 .or. n < 3 .or. any(b <= 0) .or. &
            abs(kappa) <= tiny(kappa)) return
        f%m = (n - 1)/2
        f%kappa = kappa
        allocate (f%a(0:f%m), f%c(f%m), theta(n))
        theta = [(2*pi*(j - 1)/n, j=1, n)]
        f%a(0) = sum(b)/n
        do k = 1, f%m
            f%a(k) = 2*sum(b*cos(k*theta))/n
            f%c(k) = 2*sum(b*sin(k*theta))/n
        end do
        ! Extrema from sign changes of b' on a fine grid, polished by Newton.
        n_fine = 64*n
        n_max = 0; n_min = 0
        d0 = field_db(f, 0.0_dp)
        do j = 1, n_fine
            th = 2*pi*j/n_fine
            d1 = field_db(f, th)
            if (d0 > 0 .and. d1 <= 0) then
                n_max = n_max + 1
                f%theta_x = newton_extremum(f, th - pi/n_fine)
            else if (d0 < 0 .and. d1 >= 0) then
                n_min = n_min + 1
                f%b_min = field_b(f, newton_extremum(f, th - pi/n_fine))
            end if
            d0 = d1
        end do
        if (n_max /= 1 .or. n_min /= 1) return
        f%b_max = field_b(f, f%theta_x)
        f%lambda_c = 1/f%b_max
        status = HP_OK
    end subroutine hp_field

    !> Geometric mesh: n_bulk passing levels and trapped rows in the bulk,
    !! layers towards the separatrix down to the relative width delta with
    !! ratio q, and n_col columns graded towards theta_x down to sqrt(delta).
    subroutine hp_default_mesh(f, p, n_col, n_bulk, delta, q, mesh)
        class(hp_geometry_t), intent(in) :: f
        integer, intent(in) :: p, n_col, n_bulk
        real(dp), intent(in) :: delta, q
        type(hp_mesh_t), intent(out) :: mesh
        real(dp), allocatable :: layer(:), col(:)
        real(dp) :: lc, h
        integer :: j, n_l

        mesh%p = p
        lc = f%lambda_c
        ! Relative distances to the separatrix: bulk spacing, then q^k.
        n_l = max(0, ceiling(log(delta/0.25_dp)/log(q)))
        layer = [(0.25_dp*q**j, j=0, n_l)]
        mesh%lam = [(lc*(1 - 0.25_dp)*j/n_bulk, j=0, n_bulk), &
            (lc*(1 - layer(j)), j=2, size(layer)), lc]
        mesh%zeta = [-1.0_dp, (-(1 - layer(j)), j=size(layer), 2, -1), &
            (-(1 - 0.25_dp)*(n_bulk - j)/real(n_bulk, dp), j=0, n_bulk), &
            ((1 - 0.25_dp)*j/real(n_bulk, dp), j=1, n_bulk), &
            ((1 - layer(j)), j=2, size(layer)), 1.0_dp]
        ! Columns: uniform interior, geometric towards both ends.
        h = 2*pi/n_col
        n_l = max(0, ceiling(log(0.1_dp*sqrt(delta)/0.5_dp)/log(q)))
        col = [(h*0.5_dp*q**j, j=n_l, 1, -1)]
        mesh%theta = f%theta_x - 2*pi + [0.0_dp, col, &
            (h*(0.5_dp + j), j=0, n_col - 2), 2*pi - h*0.5_dp, &
            (2*pi - col(j), j=size(col), 1, -1), 2*pi]
    end subroutine hp_default_mesh

    !> Element matrices, loads and the map of local nodes to unknowns:
    !! periodic in theta, the trapped nodes on theta_x (the X point) one
    !! unknown. mat holds each element matrix column major.
    subroutine elements(f, mesh, op, sys, mat, status, f10)
        class(hp_geometry_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        type(hp_operator_t), intent(in) :: op
        type(hp_system_t), intent(inout) :: sys
        complex(dp), allocatable, intent(out) :: mat(:, :)
        integer, intent(out) :: status
        ! f_10 per unit speed on the unknowns of the same mesh (source_dlf10).
        real(dp), intent(in), optional :: f10(:)
        integer :: p, n_col, n_row, n_lev, n_th, n_jr, ec, er, e, i, j, a1, a2, &
            n_loc, j_lo, j_hi, n_el
        integer, allocatable :: dof(:, :), curve_kind(:)
        real(dp), allocatable :: curve_val(:)

        status = HP_BAD_INPUT
        p = mesh%p
        if (p < 1 .or. op%nu <= 0 .or. op%k == 0) return
        n_col = size(mesh%theta) - 1
        n_lev = size(mesh%lam)
        ! Curves from xi = -1 to xi = 1: passing levels (sigma = -1), the
        ! trapped rows, passing levels (sigma = +1).
        curve_val = [mesh%lam(1:n_lev - 1), mesh%zeta, mesh%lam(n_lev - 1:1:-1)]
        curve_kind = [(-1, i=1, n_lev - 1), (0, i=1, size(mesh%zeta)), &
            (1, i=1, n_lev - 1)]
        n_row = size(curve_val) - 1
        n_th = n_col*p
        n_jr = n_row*p
        j_lo = (n_lev - 1)*p
        j_hi = (n_lev - 1 + size(mesh%zeta) - 1)*p
        allocate (dof(0:n_th - 1, 0:n_jr))
        sys%n = 0
        do j = 0, n_jr
            do i = 0, n_th - 1
                if (i == 0 .and. j > j_lo .and. j <= j_hi) then
                    dof(i, j) = dof(0, j_lo)
                else
                    sys%n = sys%n + 1
                    dof(i, j) = sys%n
                end if
            end do
        end do

        n_loc = (p + 1)**2
        n_el = n_col*n_row
        allocate (sys%loc(n_loc, n_el), sys%load(n_loc, 4, n_el), &
            sys%dens(n_loc, n_el), mat(n_loc**2, n_el))
        !$omp parallel do default(shared) private(e, ec, er, a1, a2) schedule(dynamic)
        do e = 1, n_el
            ec = 1 + mod(e - 1, n_col)
            er = 1 + (e - 1)/n_col
            do a2 = 0, p
                do a1 = 0, p
                    sys%loc(1 + a1 + (p + 1)*a2, e) = dof(mod((ec - 1)*p + a1, n_th), &
                        (er - 1)*p + a2)
                end do
            end do
            if (present(f10)) then
                call element(f, mesh, op, curve_kind(er:er + 1), curve_val(er:er + 1), &
                    ec, sys%load(:, :, e), sys%dens(:, e), mat(:, e), f10(sys%loc(:, e)))
            else
                call element(f, mesh, op, curve_kind(er:er + 1), curve_val(er:er + 1), &
                    ec, sys%load(:, :, e), sys%dens(:, e), mat(:, e))
            end if
        end do
        !$omp end parallel do
        allocate (sys%drive(sys%n, 4), sys%density(sys%n))
        sys%drive = 0; sys%density = 0
        do e = 1, n_el
            sys%drive(sys%loc(:, e), :) = sys%drive(sys%loc(:, e), :) + sys%load(:, :, e)
            sys%density(sys%loc(:, e)) = sys%density(sys%loc(:, e)) + sys%dens(:, e)
        end do
        sys%norm = sum(sys%density)
        status = HP_OK
    end subroutine elements

    !> Real system of the axisymmetric monoenergetic problem, b.grad theta =
    !! kappa b; drives 1 (radial) and 2 (parallel).
    subroutine hp_assemble(f, mesh, nu, sys, status)
        class(hp_geometry_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        real(dp), intent(in) :: nu
        type(hp_system_t), intent(out) :: sys
        integer, intent(out) :: status
        type(hp_operator_t) :: op
        complex(dp), allocatable :: mat(:, :)
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        integer :: e, i, j, n_loc, nt
        type(fortsparse_status_t) :: st

        op%k = f%kappa
        op%k_power = 1
        op%nu = nu
        call elements(f, mesh, op, sys, mat, status)
        if (status /= HP_OK) return
        status = HP_BAD_INPUT
        sys%mat = real(mat, dp)
        n_loc = (mesh%p + 1)**2
        nt = size(mat)
        allocate (rows(nt + 1), cols(nt + 1))
        do e = 1, size(mat, 2)
            do j = 1, n_loc
                do i = 1, n_loc
                    rows((e - 1)*n_loc**2 + i + n_loc*(j - 1)) = sys%loc(i, e)
                    cols((e - 1)*n_loc**2 + i + n_loc*(j - 1)) = sys%loc(j, e)
                end do
            end do
        end do
        vals = [reshape(sys%mat, [nt]), 0.0_dp]
        ! Constants span the null spaces of A and A^T; the equation of the
        ! first unknown, a combination of the others, pins it instead. A
        ! dense border row would defeat the fill-reducing ordering.
        where (rows(:nt) == 1) vals(:nt) = 0
        rows(nt + 1) = 1; cols(nt + 1) = 1; vals(nt + 1) = maxval(abs(vals(:nt)))
        call csc_from_triplet(sys%n, sys%n, rows, cols, vals, sys%a, st)
        if (.not. status_ok(st)) return
        status = HP_OK
    end subroutine hp_assemble

    !> Galerkin matrix of one element by tensor Gauss quadrature, column
    !> major in s, with its loads (radial and parallel drive, NTV drive and
    !> torque functional) and density weights.
    subroutine element(f, mesh, op, kinds, vals, ec, load, dens, s, f10)
        class(hp_geometry_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        type(hp_operator_t), intent(in) :: op
        real(dp), intent(in) :: vals(2)
        integer, intent(in) :: kinds(2), ec
        complex(dp), intent(out) :: load(:, :), s(:)
        real(dp), intent(out) :: dens(:)
        real(dp), intent(in), optional :: f10(:)
        integer :: p, nq, q1, q2, a1, a2, ia, ic, n_loc, j
        real(dp), allocatable :: node(:), wn(:), xq(:), wq(:), phi(:, :), dphi(:, :), &
            b_ab(:), d_ab(:), da(:), db_(:), psi(:, :)
        real(dp) :: th, th_a, xi_lo, xi_hi, d_lo, d_hi, xi, xi_b, &
            col, dr, de, den, ww, psi_a, psi_b, phase, precess, drive, torque
        complex(dp) :: loc_term, src, s_dpar, s_dlf, s_db
        real(dp) :: bb, bb_t, f10_xi

        p = mesh%p
        n_loc = (p + 1)**2
        nq = 2*p + 3
        allocate (node(p + 1), wn(p + 1), xq(nq), wq(nq), phi(p + 1, nq), &
            dphi(p + 1, nq), b_ab(n_loc), d_ab(n_loc), da(n_loc), db_(n_loc), &
            psi(p + 1, p + 1))
        call gauss_lobatto_legendre(p + 1, node, wn)
        call gauss_legendre(nq, xq, wq)
        do q1 = 1, nq
            call lagrange_weights(p + 1, xq(q1), node, phi(:, q1))
            call lagrange_deriv_weights(p + 1, xq(q1), node, dphi(:, q1))
        end do
        s = 0; load = 0; dens = 0
        th_a = (mesh%theta(ec + 1) - mesh%theta(ec))/2
        ! Streaming as the bracket of the nodal interpolant of the stream
        ! function: the bracket is invariant under the element map and
        ! polynomial, so the quadrature is exact, constants are exact left
        ! and right null vectors and the form is skew to rounding.
        do a1 = 0, p
            th = mesh%theta(ec) + th_a*(1 + node(a1 + 1))
            call f%level(kinds(1), vals(1), th, xi_lo, d_lo)
            call f%level(kinds(2), vals(2), th, xi_hi, d_hi)
            do a2 = 0, p
                psi(a1 + 1, a2 + 1) = f%stream(th, xi_lo + (xi_hi - xi_lo)*(1 + node(a2 + 1))/2)
            end do
        end do
        do q1 = 1, nq
            th = mesh%theta(ec) + th_a*(1 + xq(q1))
            call f%level(kinds(1), vals(1), th, xi_lo, d_lo)
            call f%level(kinds(2), vals(2), th, xi_hi, d_hi)
            xi_b = (xi_hi - xi_lo)/2
            do q2 = 1, nq
                xi = xi_lo + xi_b*(1 + xq(q2))
                src = f%perturbation(th, xi, op)
                ww = wq(q1)*wq(q2)*th_a*xi_b
                call f%point(th, xi, op, col, dr, de, den, phase, precess, drive, torque)
                ! -i n q xi k + i (n/v)(omega_E + pre (1 + wsign xi^2)(cos - offset)/b)
                loc_term = (0, 1)*op%n*(-op%q*phase + (op%omega_e*den + op%pre*precess) &
                    /op%v)
                do a2 = 0, p
                    do a1 = 0, p
                        ia = 1 + a1 + (p + 1)*a2
                        b_ab(ia) = phi(a1 + 1, q1)*phi(a2 + 1, q2)
                        da(ia) = dphi(a1 + 1, q1)*phi(a2 + 1, q2)
                        db_(ia) = phi(a1 + 1, q1)*dphi(a2 + 1, q2)
                        d_ab(ia) = db_(ia)/xi_b
                    end do
                end do
                psi_a = sum(reshape(psi, [n_loc])*da)
                psi_b = sum(reshape(psi, [n_loc])*db_)
                do ic = 1, n_loc
                    s((ic - 1)*n_loc + 1:ic*n_loc) = s((ic - 1)*n_loc + 1:ic*n_loc) &
                        + wq(q1)*wq(q2)*b_ab*(psi_b*da(ic) - psi_a*db_(ic)) &
                        + ww*(col*d_ab*d_ab(ic) + loc_term*b_ab*b_ab(ic))
                end do
                load(:, 1) = load(:, 1) + ww*dr*b_ab
                load(:, 2) = load(:, 2) + ww*de*b_ab
                load(:, 3) = load(:, 3) + ww*drive*src/op%v*b_ab
                if (op%source_dpar .or. present(f10)) then
                    call field_value(f, th, bb, bb_t)
                    call perturbed_sources(op, th, bb_t/bb, s_dpar, s_dlf, s_db)
                    if (op%source_dpar) load(:, 3) = load(:, 3) + ww*drive*s_dpar/op%v*b_ab
                    if (present(f10)) then
                        ! Mirror force of the perturbed parallel gradient on
                        ! f_10, (k/b) s_dlf (1 - xi^2)/2 d_xi f_10, and minus
                        ! the perturbed collision term, (dB/B)(nu/2) L f_10
                        ! with L = -d_xi (1 - xi^2) d_xi, by parts in xi.
                        f10_xi = sum(f10*d_ab)
                        load(:, 3) = load(:, 3) + ww*den*(op%k*bb**op%k_power/bb)*s_dlf &
                            *(1 - xi*xi)/2*f10_xi*b_ab &
                            + ww*den*(s_db/bb)*op%nu/2*(1 - xi*xi)*f10_xi*d_ab
                    end if
                end if
                load(:, 4) = load(:, 4) + ww*torque*conjg(src)*b_ab
                dens = dens + ww*den*b_ab
            end do
        end do
    end subroutine element

    subroutine hp_factor(sys, status)
        type(hp_system_t), intent(inout) :: sys
        integer, intent(out) :: status
        type(fortsparse_status_t) :: st
        status = HP_SINGULAR
        sys%solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        sys%solver%options%ordering = FORTSPARSE_ORDER_AMD
        sys%solver%options%lu_method = FORTSPARSE_LU_SUPERNODAL
        call sparse_factor(sys%solver, sys%a, st)
        if (.not. status_ok(st)) then
            call sparse_destroy(sys%solver)
            sys%solver%backend_id = FORTSPARSE_BACKEND_NATIVE
            sys%solver%options%ordering = FORTSPARSE_ORDER_AMD
            sys%solver%options%lu_method = FORTSPARSE_LU_SCALAR
            call sparse_factor(sys%solver, sys%a, st)
            if (.not. status_ok(st)) return
        end if
        sys%factored = .true.
        status = HP_OK
    end subroutine hp_factor

    !> Solve A x + mu d = rhs, <d, x> = 0, or the same with the transpose of
    !! A (the adjoint); d is the density weight, mu the density source that
    !! makes the right-hand side compatible.
    subroutine hp_solve(sys, rhs, x, mu, status, adjoint)
        type(hp_system_t), intent(inout) :: sys
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: x(:), mu
        integer, intent(out) :: status
        logical, intent(in), optional :: adjoint
        type(fortsparse_status_t) :: st
        integer :: op, it
        real(dp), allocatable :: r(:), dr(:), dx(:)
        op = FORTSPARSE_OP_NORMAL
        if (present(adjoint)) then
            if (adjoint) op = FORTSPARSE_OP_TRANSPOSE
        end if
        status = HP_SINGULAR
        mu = sum(rhs)/sum(sys%density)
        r = rhs - mu*sys%density
        r(1) = 0
        call sparse_solve_op(sys%solver, r, x, op, st)
        if (.not. status_ok(st)) return
        ! The supernodal factors use static pivots; refine on the true residual.
        allocate (dx(size(x)))
        do it = 1, 3
            dr = r - apply(sys%a, x, op == FORTSPARSE_OP_TRANSPOSE)
            call sparse_solve_op(sys%solver, dr, dx, op, st)
            if (.not. status_ok(st)) return
            x = x + dx
        end do
        x = x - dot_product(sys%density, x)/sys%norm
        status = HP_OK
    end subroutine hp_solve

    !> J = <s_j g_i> (D_ji) on the mesh and the dual-weighted residual
    !! estimate eta of J - J_h: the residual of g_h in the space of degree
    !! p + 1 weighted by the adjoint z+ - I z_h of that space, where
    !! A^T z = s_j. eta_el holds element indicators for refinement.
    !! This is an indicator, not a bound.
    subroutine hp_estimate(f, mesh, nu, i, j, value, eta, eta_el, status, n_dof)
        class(hp_geometry_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        real(dp), intent(in) :: nu
        integer, intent(in) :: i, j
        real(dp), intent(out) :: value, eta
        real(dp), allocatable, intent(out) :: eta_el(:)
        integer, intent(out) :: status
        integer, intent(out), optional :: n_dof
        type(hp_system_t) :: sys, fine
        type(hp_mesh_t) :: mesh_fine
        real(dp), allocatable :: g(:), z(:), zf(:), gf(:), izf(:), node(:), &
            node_f(:), w(:), lag(:, :), r(:)
        real(dp) :: mu, mu_z, mu_f
        integer :: p, e, a, n_loc

        call hp_assemble(f, mesh, nu, sys, status)
        if (status /= HP_OK) return
        call hp_factor(sys, status)
        if (status /= HP_OK) return
        if (present(n_dof)) n_dof = sys%n
        allocate (g(sys%n), z(sys%n))
        call hp_solve(sys, real(sys%drive(:, i), dp), g, mu, status)
        if (status == HP_OK) call hp_solve(sys, real(sys%drive(:, j), dp), z, mu_z, status, &
            adjoint=.true.)
        call hp_free(sys)
        if (status /= HP_OK) return
        value = dot_product(real(sys%drive(:, j), dp), g)/sys%norm

        mesh_fine = mesh
        mesh_fine%p = mesh%p + 1
        call hp_assemble(f, mesh_fine, nu, fine, status)
        if (status /= HP_OK) return
        call hp_factor(fine, status)
        if (status /= HP_OK) return
        allocate (zf(fine%n))
        call hp_solve(fine, real(fine%drive(:, j), dp), zf, mu_f, status, adjoint=.true.)
        call hp_free(fine)
        if (status /= HP_OK) return

        ! Interpolation of degree p into degree p + 1 per element.
        p = mesh%p
        allocate (node(p + 1), node_f(p + 2), w(p + 2), lag(p + 2, p + 1))
        call gauss_lobatto_legendre(p + 1, node, w(:p + 1))
        call gauss_lobatto_legendre(p + 2, node_f, w)
        do a = 1, p + 2
            call lagrange_weights(p + 1, node_f(a), node, lag(a, :))
        end do
        n_loc = (p + 2)**2
        allocate (gf(fine%n), izf(fine%n), eta_el(size(fine%loc, 2)), r(n_loc))
        eta = 0
        do e = 1, size(fine%loc, 2)
            gf(fine%loc(:, e)) = reshape(matmul(lag, matmul(reshape(g(sys%loc(:, e)), &
                [p + 1, p + 1]), transpose(lag))), [n_loc])
            izf(fine%loc(:, e)) = reshape(matmul(lag, matmul(reshape(z(sys%loc(:, e)), &
                [p + 1, p + 1]), transpose(lag))), [n_loc])
        end do
        do e = 1, size(fine%loc, 2)
            r = real(fine%load(:, i, e), dp) - mu*fine%dens(:, e) &
                - matmul(reshape(fine%mat(:, e), [n_loc, n_loc]), gf(fine%loc(:, e)))
            eta = eta + dot_product(zf(fine%loc(:, e)) - izf(fine%loc(:, e)), r)/fine%norm
            ! Marking uses node-wise magnitudes; signed element sums cancel.
            eta_el(e) = sum(abs((zf(fine%loc(:, e)) - izf(fine%loc(:, e)))*r))/fine%norm
        end do
    end subroutine hp_estimate

    !> Bisect the rows and columns that carry the largest indicator sums
    !! until they hold the fraction of the total (Doerfler marking). Rows
    !! are bisected in their curve parameter, mirrored in sigma.
    subroutine hp_refine(mesh, eta_el, fraction, out)
        type(hp_mesh_t), intent(in) :: mesh
        real(dp), intent(in) :: eta_el(:), fraction
        type(hp_mesh_t), intent(out) :: out
        integer :: n_col, n_lev, n_z, n_row, e, er, ec, m
        real(dp), allocatable :: s_col(:), s_lam(:), s_z(:), all(:)
        real(dp) :: total, acc, cut

        n_col = size(mesh%theta) - 1
        n_lev = size(mesh%lam); n_z = size(mesh%zeta)
        n_row = 2*(n_lev - 1) + (n_z - 1)
        allocate (s_col(n_col), s_lam(n_lev - 1), s_z(n_z - 1))
        s_col = 0; s_lam = 0; s_z = 0
        do e = 1, size(eta_el)
            ec = 1 + mod(e - 1, n_col)
            er = 1 + (e - 1)/n_col
            s_col(ec) = s_col(ec) + abs(eta_el(e))
            if (er <= n_lev - 1) then
                s_lam(er) = s_lam(er) + abs(eta_el(e))
            else if (er <= n_lev - 1 + n_z - 1) then
                s_z(er - (n_lev - 1)) = s_z(er - (n_lev - 1)) + abs(eta_el(e))
            else
                s_lam(n_row + 1 - er) = s_lam(n_row + 1 - er) + abs(eta_el(e))
            end if
        end do
        ! Trapped rows mirror in xi.
        s_z = s_z + s_z(n_z - 1:1:-1)
        ! Smallest threshold whose larger entries hold the fraction.
        all = [s_col, s_lam, s_z]
        total = sum(all)
        call sort_desc(all)
        acc = 0; cut = all(1)
        do m = 1, size(all)
            acc = acc + all(m)
            cut = all(m)
            if (acc >= fraction*total) exit
        end do
        out%p = mesh%p
        out%theta = bisect(mesh%theta, s_col >= cut)
        out%lam = bisect(mesh%lam, s_lam >= cut)
        out%zeta = bisect(mesh%zeta, s_z >= cut)
    end subroutine hp_refine

    pure function bisect(x, mark) result(y)
        real(dp), intent(in) :: x(:)
        logical, intent(in) :: mark(:)
        real(dp), allocatable :: y(:)
        integer :: k
        y = [x(1)]
        do k = 1, size(mark)
            if (mark(k)) y = [y, (x(k) + x(k + 1))/2]
            y = [y, x(k + 1)]
        end do
    end function bisect

    pure subroutine sort_desc(x)
        real(dp), intent(inout) :: x(:)
        integer :: k, m
        real(dp) :: t
        do k = 2, size(x)
            t = x(k)
            m = k - 1
            do while (m >= 1)
                if (x(m) >= t) exit
                x(m + 1) = x(m)
                m = m - 1
            end do
            x(m + 1) = t
        end do
    end subroutine sort_desc

    !> Flux kernel R = int dtheta dxi (2/b)(1 + wsign xi^2)/b Re(S^* F) of one
    !! toroidal harmonic (the kernel of milca_dk_ntv), op%k the constant
    !! b.grad theta (k_power = 0) or k b (k_power = 1). Without a null space
    !! the complex system is solved directly.
    subroutine hp_ntv_kernel(f, mesh, op, r, status, n_dof)
        class(hp_geometry_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        type(hp_operator_t), intent(in) :: op
        real(dp), intent(out) :: r
        integer, intent(out) :: status
        integer, intent(out), optional :: n_dof
        type(hp_system_t) :: sys
        type(csc_z_t) :: a, scaled
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: st
        complex(dp), allocatable :: mat(:, :), x(:), res(:), dx(:)
        real(dp), allocatable :: rs(:), cs(:)
        integer, allocatable :: rows(:), cols(:)
        integer :: e, i, j, n_loc, nt, it

        status = HP_BAD_INPUT
        if (op%n == 0 .or. .not. allocated(op%m)) return
        if (op%source_dlf10) then
            ! f_10 of the axisymmetric problem on this mesh, in the kernel's
            ! normalisation per unit speed (milca_dk_ntv ntv_f10).
            block
                type(hp_system_t) :: ax
                real(dp), allocatable :: f10(:)
                real(dp) :: mu
                select type (f)
                class is (hp_field_t)
                class default
                    return
                end select
                if (op%q == 0 .or. f%kappa == 0) return
                call hp_assemble(f, mesh, op%nu, ax, status)
                if (status /= HP_OK) return
                call hp_factor(ax, status)
                if (status /= HP_OK) return
                allocate (f10(ax%n))
                call hp_solve(ax, real(ax%drive(:, 1), dp), f10, mu, status)
                call hp_free(ax)
                if (status /= HP_OK) return
                f10 = -2/(op%v*op%q*f%kappa)*f10
                call elements(f, mesh, op, sys, mat, status, f10)
            end block
        else
            call elements(f, mesh, op, sys, mat, status)
        end if
        if (status /= HP_OK) return
        if (present(n_dof)) n_dof = sys%n
        status = HP_BAD_INPUT
        n_loc = (mesh%p + 1)**2
        nt = size(mat)
        allocate (rows(nt), cols(nt))
        do e = 1, size(mat, 2)
            do j = 1, n_loc
                do i = 1, n_loc
                    rows((e - 1)*n_loc**2 + i + n_loc*(j - 1)) = sys%loc(i, e)
                    cols((e - 1)*n_loc**2 + i + n_loc*(j - 1)) = sys%loc(j, e)
                end do
            end do
        end do
        call csc_from_triplet(sys%n, sys%n, rows, cols, reshape(mat, [nt]), a, st)
        if (.not. status_ok(st)) return
        status = HP_SINGULAR
        ! Equilibrated rows and columns keep the static-pivot growth of the
        ! supernodal LU in bounds on graded meshes; without it the factor
        ! falls back to the scalar LU.
        scaled = a
        call equilibrate(scaled, rs, cs)
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%ordering = FORTSPARSE_ORDER_AMD
        solver%options%lu_method = FORTSPARSE_LU_SUPERNODAL
        call sparse_factor(solver, scaled, st)
        if (.not. status_ok(st)) then
            ! Static pivots of the supernodal LU can still break down on wide
            ! orbits at high energy; threshold pivoting does not.
            call sparse_destroy(solver)
            solver%backend_id = FORTSPARSE_BACKEND_NATIVE
            solver%options%ordering = FORTSPARSE_ORDER_AMD
            solver%options%lu_method = FORTSPARSE_LU_SCALAR
            call sparse_factor(solver, scaled, st)
            if (.not. status_ok(st)) return
        end if
        allocate (x(sys%n), dx(sys%n))
        res = rs*sys%drive(:, 3)
        call sparse_solve_op(solver, res, x, FORTSPARSE_OP_NORMAL, st)
        x = cs*x
        do it = 1, 3
            if (.not. status_ok(st)) exit
            res = rs*(sys%drive(:, 3) - apply_z(a, x))
            call sparse_solve_op(solver, res, dx, FORTSPARSE_OP_NORMAL, st)
            x = x + cs*dx
        end do
        call sparse_destroy(solver)
        if (.not. status_ok(st)) return
        r = real(sum(sys%drive(:, 4)*x), dp)
        status = HP_OK
    end subroutine hp_ntv_kernel

    !> Row and column scalings rs, cs with rs A cs of unit row and column
    !> maxima, by three sweeps of Ruiz's iteration; a is replaced.
    subroutine equilibrate(a, rs, cs)
        type(csc_z_t), intent(inout) :: a
        real(dp), allocatable, intent(out) :: rs(:), cs(:)
        real(dp), allocatable :: r(:), c(:)
        integer :: sweep, j, k
        allocate (rs(a%nrow), cs(a%ncol), r(a%nrow), c(a%ncol))
        rs = 1; cs = 1
        do sweep = 1, 3
            r = 0; c = 0
            do j = 1, a%ncol
                do k = a%col_ptr(j), a%col_ptr(j + 1) - 1
                    r(a%row_idx(k)) = max(r(a%row_idx(k)), abs(a%val(k)))
                    c(j) = max(c(j), abs(a%val(k)))
                end do
            end do
            where (r > 0) r = 1/sqrt(r)
            where (c > 0) c = 1/sqrt(c)
            do j = 1, a%ncol
                do k = a%col_ptr(j), a%col_ptr(j + 1) - 1
                    a%val(k) = a%val(k)*r(a%row_idx(k))*c(j)
                end do
            end do
            rs = rs*r
            cs = cs*c
        end do
    end subroutine equilibrate

    pure function apply_z(a, x) result(y)
        type(csc_z_t), intent(in) :: a
        complex(dp), intent(in) :: x(:)
        complex(dp) :: y(a%nrow)
        integer :: j, k
        y = 0
        do j = 1, a%ncol
            do k = a%col_ptr(j), a%col_ptr(j + 1) - 1
                y(a%row_idx(k)) = y(a%row_idx(k)) + a%val(k)*x(j)
            end do
        end do
    end function apply_z

    !> A x or A^T x of a CSC matrix.
    pure function apply(a, x, transpose) result(y)
        type(csc_t), intent(in) :: a
        real(dp), intent(in) :: x(:)
        logical, intent(in) :: transpose
        real(dp) :: y(a%nrow)
        integer :: j, k
        y = 0
        do j = 1, a%ncol
            do k = a%col_ptr(j), a%col_ptr(j + 1) - 1
                if (transpose) then
                    y(j) = y(j) + a%val(k)*x(a%row_idx(k))
                else
                    y(a%row_idx(k)) = y(a%row_idx(k)) + a%val(k)*x(j)
                end if
            end do
        end do
    end function apply

    subroutine hp_free(sys)
        type(hp_system_t), intent(inout) :: sys
        if (sys%factored) call sparse_destroy(sys%solver)
        sys%factored = .false.
    end subroutine hp_free

    !> D_ij = <s_i g_j> on the given mesh; border holds the density-source
    !! multipliers of the two drives.
    subroutine hp_coefficients(f, mesh, nu, c, status)
        class(hp_geometry_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        real(dp), intent(in) :: nu
        type(mono_coefficients_t), intent(out) :: c
        integer, intent(out) :: status
        type(hp_system_t) :: sys
        real(dp), allocatable :: g1(:), g3(:)
        integer :: n

        call hp_assemble(f, mesh, nu, sys, status)
        if (status /= HP_OK) return
        call hp_factor(sys, status)
        if (status /= HP_OK) return
        n = sys%n
        allocate (g1(n), g3(n))
        call hp_solve(sys, real(sys%drive(:, 1), dp), g1, c%border(1), status)
        if (status == HP_OK) call hp_solve(sys, real(sys%drive(:, 2), dp), g3, c%border(2), status)
        call hp_free(sys)
        if (status /= HP_OK) return
        c%d11 = dot_product(real(sys%drive(:, 1), dp), g1)/sys%norm
        c%d31 = dot_product(real(sys%drive(:, 2), dp), g1)/sys%norm
        c%d13 = dot_product(real(sys%drive(:, 1), dp), g3)/sys%norm
        c%d33 = dot_product(real(sys%drive(:, 2), dp), g3)/sys%norm
    end subroutine hp_coefficients

    !> xi and d xi/d theta on a row curve: a passing level of sign kind
    !! (kind = +-1, val = Lambda) or a scaled separatrix (kind = 0, val =
    !! zeta). Differences of b are formed without cancellation near the X
    !! point.
    subroutine curve(g, kind, val, th, xi, dxi)
        class(hp_field_t), intent(in) :: g
        integer, intent(in) :: kind
        real(dp), intent(in) :: val, th
        real(dp), intent(out) :: xi, dxi
        real(dp) :: v, q, db
        v = separatrix2(g, th)
        db = field_db(g, th)
        if (kind /= 0) then
            q = g%b_max*((1/g%b_max - val) + val*v)
            xi = kind*sqrt(q)
            dxi = -kind*val*db/(2*sqrt(q))
        else if (v > 0) then
            xi = val*sqrt(v)
            dxi = -val*db/(2*g%b_max*sqrt(v))
        else
            xi = 0; dxi = 0
        end if
    end subroutine curve

    !> The drive of one toroidal harmonic at (theta, xi), i n sum_m c_m
    !> exp(i m theta) from the operator's amplitudes; a geometry whose orbits
    !> leave the surface overrides it with the amplitudes at the orbit point.
    !> b and db/dtheta of a thin field (zero slope otherwise).
    subroutine field_value(f, th, b, b_t)
        class(hp_geometry_t), intent(in) :: f
        real(dp), intent(in) :: th
        real(dp), intent(out) :: b, b_t
        select type (f)
        class is (hp_field_t)
            b = field_b(f, th)
            b_t = field_db(f, th)
        class default
            b = 1
            b_t = 0
        end select
    end subroutine field_value

    !> At theta with d_theta ln b = slope: the parallel-derivative source
    !> sum_m (i (m - n q) - 2 slope)/q c_m e^{i m theta}, the f_10 mirror
    !> factor sum_m (i (m - n q) - slope) c_m e^{i m theta} and dB/B.
    pure subroutine perturbed_sources(op, th, slope, s_dpar, s_dlf, s_db)
        type(hp_operator_t), intent(in) :: op
        real(dp), intent(in) :: th, slope
        complex(dp), intent(out) :: s_dpar, s_dlf, s_db
        complex(dp) :: e
        integer :: j
        s_dpar = 0; s_dlf = 0; s_db = 0
        do j = 1, size(op%m)
            e = op%c(j)*exp((0, 1)*op%m(j)*th)
            s_dpar = s_dpar + ((0, 1)*(op%m(j) - op%n*op%q) - 2*slope)/op%q*e
            s_dlf = s_dlf + ((0, 1)*(op%m(j) - op%n*op%q) - slope)*e
            s_db = s_db + e
        end do
    end subroutine perturbed_sources

    complex(dp) function op_perturbation(g, th, xi, op) result(src)
        class(hp_geometry_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(hp_operator_t), intent(in) :: op
        integer :: j
        src = 0
        if (.not. allocated(op%m)) return
        do j = 1, size(op%m)
            src = src + (0, 1)*op%n*op%c(j)*exp((0, 1)*op%m(j)*th)
        end do
    end function op_perturbation

    real(dp) function field_stream(g, th, xi) result(psi)
        class(hp_field_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        call dk_hp_stream(field_b(g, th), xi, psi)
    end function field_stream

    subroutine field_point(g, th, xi, op, collide, drive_r, drive_e, density, phase, &
            precess, drive, torque)
        class(hp_field_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(hp_operator_t), intent(in) :: op
        real(dp), intent(out) :: collide, drive_r, drive_e, density, phase, precess, drive, &
            torque
        real(dp) :: bb, k
        bb = field_b(g, th)
        k = op%k*bb**op%k_power
        call dk_hp_point(bb, field_db(g, th), xi, k, op%nu, collide, drive_r, drive_e, density)
        call dk_hp_ntv(bb, xi, k, cos(th), op%offset, op%wsign, phase, precess, drive, torque)
        ! The flux is weighted with the surface measure, which streaming
        ! b.grad theta = k b**k_power divides by b**k_power (Boozer: 2/b^3).
        torque = torque/bb**op%k_power
    end subroutine field_point

    !> xi_c^2 = 1 - b/b_max.
    pure real(dp) function separatrix2(f, th) result(s)
        type(hp_field_t), intent(in) :: f
        real(dp), intent(in) :: th
        real(dp) :: d
        d = modulo(th - f%theta_x + pi, 2*pi) - pi
        s = max(0.0_dp, -field_diff(f, f%theta_x, d)/f%b_max)
    end function separatrix2

    pure real(dp) function field_b(f, th) result(b)
        type(hp_field_t), intent(in) :: f
        real(dp), intent(in) :: th
        integer :: k
        b = f%a(0)
        do k = 1, f%m
            b = b + f%a(k)*cos(k*th) + f%c(k)*sin(k*th)
        end do
    end function field_b

    pure real(dp) function field_db(f, th) result(b)
        type(hp_field_t), intent(in) :: f
        real(dp), intent(in) :: th
        integer :: k
        b = 0
        do k = 1, f%m
            b = b + k*(f%c(k)*cos(k*th) - f%a(k)*sin(k*th))
        end do
    end function field_db

    pure real(dp) function field_d2b(f, th) result(b)
        type(hp_field_t), intent(in) :: f
        real(dp), intent(in) :: th
        integer :: k
        b = 0
        do k = 1, f%m
            b = b - k*k*(f%a(k)*cos(k*th) + f%c(k)*sin(k*th))
        end do
    end function field_d2b

    !> b(y + d) - b(y) from products of sines, exact to rounding for small d.
    pure real(dp) function field_diff(f, y, d) result(s)
        type(hp_field_t), intent(in) :: f
        real(dp), intent(in) :: y, d
        integer :: k
        s = 0
        do k = 1, f%m
            s = s + 2*sin(k*d/2)*(f%c(k)*cos(k*(y + d/2)) - f%a(k)*sin(k*(y + d/2)))
        end do
    end function field_diff

    real(dp) function newton_extremum(f, th0) result(th)
        type(hp_field_t), intent(in) :: f
        real(dp), intent(in) :: th0
        integer :: it
        real(dp) :: step
        th = th0
        do it = 1, 50
            step = field_db(f, th)/field_d2b(f, th)
            th = th - step
            if (abs(step) < 1e-15_dp) exit
        end do
        th = modulo(th, 2*pi)
    end function newton_extremum

end module milca_dk_hp
