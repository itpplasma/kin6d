module kin6d_equilibrium_nonnested
    ! NN1: closed FEEC magnetic states and the quadratic vacuum / N=0
    ! pressure restriction. Generic spaces and sparse solves stay upstream.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_point_basis_t, meridian_point_basis_init, meridian_point_basis, &
        meridian_op_complex, meridian_gs_stiffness, meridian_stiffness, &
        meridian_qpoint, meridian_real_to_complex, meridian_load, meridian_field_i
    use fortsparse, only: csc_t, csc_z_t, csc_matvec, csc_from_triplet, &
        csc_normalize_symmetric_raw, &
        sparse_solver_t, sparse_factor, sparse_solve, sparse_destroy, &
        fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_NATIVE
    use kin6d_gs_direct, only: gs_direct_t, gs_direct_available, gs_direct_factor, &
        gs_direct_solve, gs_direct_free
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_locate
    use kin6d_equilibrium_state, only: equilibrium_field_t, equilibrium_state_t, &
        EQ_BR, EQ_BZ, EQ_BPHI, EQ_P, EQ_PSI, EQ_PSI_R, EQ_PSI_Z, EQ_VALUES
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_nonnested_leaves, only: axisymmetric_pressure, fieldline_velocity
    implicit none
    private
    public :: free_b_field_t, nn_report_t, nn_build, nn_solve, nn_export
    public :: nn_linear_solve, nn_divergence, nn_point, nn_vacuum_mode, nn_riesz
    public :: nn_velocity, nn_fourier_value
    public :: nn_factor_t, nn_factor, nn_factor_solve, nn_factor_solve_real
    public :: nn_factor_free

    type, extends(equilibrium_field_t) :: free_b_field_t
        type(flux_mesh_t) :: mesh
        type(meridian_feec_t) :: fs
        type(meridian_point_basis_t) :: pb
        complex(dp), allocatable :: b(:, :) ! V2, modes 0:N
        real(dp), allocatable :: psi(:)     ! N=0 pressure primitive, internal
        real(dp), allocatable :: p(:)       ! NN2 transported p, V0 Fourier blocks
        real(dp) :: f_vac = 0, pprime = 0, psi_bound = 0, mu0 = 1
        real(dp) :: centre(2) = 0, radial_scale = 1
        procedure(polar_coordinates_i), pointer, nopass :: coordinates => &
            circular_coordinates
    contains
        procedure :: sample => nn_sample
    end type

    abstract interface
        subroutine polar_coordinates_i(centre, scale, x, u, theta)
            import :: dp
            real(dp), intent(in) :: centre(2), scale, x(2)
            real(dp), intent(out) :: u, theta
        end subroutine
    end interface

    type :: nn_factor_t
        ! Reusable factor of a masked symmetric operator (see nn_factor).
        logical :: ready = .false., direct_backend = .false.
        logical, allocatable :: fixed(:)
        type(gs_direct_t) :: direct
        type(sparse_solver_t) :: solver
    end type

    type :: nn_report_t
        real(dp) :: energy_before = 0, energy_after = 0, residual = 0
        real(dp) :: divergence = 0
        integer :: newton_steps = 0
        character(32) :: claim = 'candidate'
    end type
contains
    subroutine nn_build(mesh, order, modes, centre, radial_scale, f_vac, mu0, &
                        field, status)
        type(flux_mesh_t), intent(in) :: mesh
        integer, intent(in) :: order, modes
        real(dp), intent(in) :: centre(2), radial_scale, f_vac, mu0
        type(free_b_field_t), intent(out) :: field
        integer, intent(out) :: status
        status = 1
        if (modes < 0 .or. radial_scale <= 0 .or. mu0 <= 0) return
        if (.not. all(ieee_is_finite([centre, radial_scale, f_vac, mu0]))) return
        field%mesh = mesh
        call meridian_feec_build(mesh%mm, order, field%fs, status)
        if (status /= 0) return
        call meridian_point_basis_init(field%fs, field%pb, status)
        if (status /= 0) return
        allocate (field%b(field%fs%sp(2)%ndof, 0:modes), &
            field%psi(field%fs%n0))
        field%b = 0; field%psi = 0
        field%centre = centre; field%radial_scale = radial_scale
        field%f_vac = f_vac; field%mu0 = mu0
    end subroutine

    subroutine nn_solve(field, n, boundary, seed, report, status, pprime, psi_bound, &
                        solution)
        ! For n=0, boundary/seed are scalar R A_phi; for n>0 they are the
        ! Nedelec A_pol in gauge R A_phi=0. Boundary trace is prescribed.
        ! solution returns that potential (the same unknowns as seed).
        type(free_b_field_t), intent(inout) :: field
        integer, intent(in) :: n
        complex(dp), intent(in) :: boundary(:), seed(:)
        type(nn_report_t), intent(out) :: report
        integer, intent(out) :: status
        real(dp), intent(in), optional :: pprime, psi_bound
        complex(dp), allocatable, intent(out), optional :: solution(:)
        type(csc_t) :: k
        type(csc_z_t) :: zk, curl
        complex(dp), allocatable :: rhs(:), a(:), delta(:), full(:), r(:)
        logical, allocatable :: mask(:)
        real(dp), allocatable :: loads(:, :)
        real(dp) :: pp, en_scale
        integer :: st, t, ndof, failure
        status = 1
        if (.not. allocated(field%b)) return
        if (n < 0 .or. n > ubound(field%b, 2)) return
        pp = field%pprime
        if (present(pprime)) pp = pprime
        if (.not. ieee_is_finite(pp)) return
        if (present(psi_bound)) then
            if (.not. ieee_is_finite(psi_bound)) return
        end if
        if (ubound(field%b, 2) > 0 .and. abs(pp) > 0) return
        if (n == 0) then
            call meridian_gs_stiffness(field%fs, k, st)
            ndof = field%fs%n0
            mask = field%fs%sp(0)%bnd
        else
            call meridian_stiffness(field%fs, 1, n, zk, st)
            if (st /= 0) return
            ndof = field%fs%nn
            call real_leading(zk, ndof, k, st)
            mask = field%fs%sp(1)%bnd(:ndof)
        end if
        if (st /= 0) return
        if (size(seed) /= ndof .or. size(boundary) /= ndof) return
        if (.not. all(ieee_is_finite(real(seed)))) return
        if (.not. all(ieee_is_finite(aimag(seed)))) return
        if (.not. all(ieee_is_finite(real(boundary)))) return
        if (.not. all(ieee_is_finite(aimag(boundary)))) return
        if (n == 0) then
            if (any(abs(aimag(boundary)) > 0) .or. any(abs(aimag(seed)) > 0)) return
        end if
        allocate (rhs(ndof), a(ndof), delta(ndof))
        rhs = 0
        if (n == 0) then
            allocate (loads(field%fs%nl0, field%fs%mm%topology%n_triangles))
            failure = 0
            !$omp parallel do default(shared) private(t, st) schedule(static) &
            !$omp& reduction(max:failure)
            do t = 1, size(loads, 2)
                call pressure_load(t, loads(:, t), st)
                failure = max(failure, st)
            end do
            !$omp end parallel do
            if (failure /= 0) return
            do t = 1, size(loads, 2)
                rhs(field%fs%sp(0)%dof(:, t)) = rhs(field%fs%sp(0)%dof(:, t)) &
                    + loads(:, t)
            end do
        end if
        call meridian_real_to_complex(k, zk)
        a = seed
        where (mask) a = boundary
        r = rhs - csc_matvec(zk, a)
        where (mask) r = 0
        call nn_linear_solve(k, mask, r, delta, st)
        if (st /= 0) return
        report%energy_before = quadratic(a)
        a = a + delta
        report%energy_after = quadratic(a)
        en_scale = max(1.0_dp, abs(report%energy_before))
        if (report%energy_after > report%energy_before + 1.0e-11_dp*en_scale) return
        r = rhs - csc_matvec(zk, a)
        where (mask) r = 0
        report%residual = norm2(abs(r))/max(1.0_dp, norm2(abs(rhs)), &
            norm2(abs(csc_matvec(zk, a))))
        report%newton_steps = 1
        if (report%residual > 1.0e-9_dp) return
        call meridian_op_complex(field%fs%curl, n, curl, st)
        if (st /= 0) return
        allocate (full(field%fs%sp(1)%ndof))
        full = 0
        if (n == 0) then
            full(field%fs%nn + 1:) = a
            field%psi = real(a, dp)
            field%pprime = pp
            if (present(psi_bound)) field%psi_bound = psi_bound
        else
            full(:ndof) = a
        end if
        field%b(:, n) = csc_matvec(curl, full)
        call nn_divergence(field, report%divergence, st)
        if (st /= 0) return
        if (present(solution)) solution = a
        status = 0
    contains
        subroutine pressure_load(t, load, status)
            integer, intent(in) :: t
            real(dp), intent(out) :: load(:)
            integer, intent(out) :: status
            real(dp) :: x(2), weight, phi(3, field%fs%nl0)
            integer :: q
            load = 0; status = 0
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 0, t, q, x, weight, phi, status)
                if (status /= 0) return
                load = load + field%mu0*pp*weight*phi(1, :)
            end do
        end subroutine

        real(dp) function quadratic(v) result(energy)
            complex(dp), intent(in) :: v(:)
            energy = real(dot_product(v, csc_matvec(zk, v))/2 &
                - dot_product(v, rhs), dp)/field%mu0
        end function
    end subroutine

    subroutine nn_vacuum_mode(field, n, potential, report, status, solution)
        type(free_b_field_t), intent(inout) :: field
        integer, intent(in) :: n
        procedure(meridian_field_i) :: potential
        type(nn_report_t), intent(out) :: report
        integer, intent(out) :: status
        complex(dp), allocatable, intent(out), optional :: solution(:)
        type(csc_z_t) :: z
        type(csc_t) :: mass
        complex(dp), allocatable :: rhs(:), a(:), seed(:)
        logical, allocatable :: fixed(:)
        status = 1
        if (n <= 0) return
        call meridian_load(field%fs, 1, potential, rhs, status)
        if (status /= 0) return
        call meridian_real_to_complex(field%fs%mass(1), z)
        call real_leading(z, field%fs%nn, mass, status)
        if (status /= 0) return
        allocate (a(field%fs%nn), seed(field%fs%nn), fixed(field%fs%nn))
        fixed = .false.
        call nn_linear_solve(mass, fixed, rhs(:field%fs%nn), a, status)
        if (status /= 0) return
        seed = 0
        call nn_solve(field, n, a, seed, report, status, solution=solution)
    end subroutine

    subroutine nn_riesz(field, n, rhs, x, status)
        type(free_b_field_t), intent(in) :: field
        integer, intent(in) :: n
        complex(dp), intent(in) :: rhs(:)
        complex(dp), intent(out) :: x(:)
        integer, intent(out) :: status
        type(csc_z_t) :: z
        type(csc_t) :: k
        call meridian_stiffness(field%fs, 0, n, z, status)
        if (status /= 0) return
        call real_leading(z, field%fs%n0, k, status)
        if (status /= 0) return
        call nn_linear_solve(k, field%fs%sp(0)%bnd, rhs, x, status)
    end subroutine

    subroutine nn_divergence(field, error, status)
        type(free_b_field_t), intent(in) :: field
        real(dp), intent(out) :: error
        integer, intent(out) :: status
        type(csc_z_t) :: div
        integer :: n
        error = 0; status = 1
        if (.not. allocated(field%b)) return
        do n = 0, ubound(field%b, 2)
            call meridian_op_complex(field%fs%div, n, div, status)
            if (status /= 0) return
            error = max(error, maxval(abs(csc_matvec(div, field%b(:, n)))) &
                /max(1.0_dp, maxval(abs(field%b(:, n)))))
        end do
        status = 0
    end subroutine

    subroutine nn_export(field, seed, eq)
        type(free_b_field_t), intent(in) :: field
        character(*), intent(in) :: seed
        type(equilibrium_state_t), intent(out) :: eq
        allocate (eq%field, source=field)
        allocate (cylindrical_chart_t :: eq%chart)
        eq%lane = 'free-B'; eq%seed = seed; eq%claim = 'candidate'
        eq%f_vac = field%f_vac; eq%mu0 = field%mu0
        eq%psi_bound = field%psi_bound
    end subroutine

    subroutine nn_point(field, point, value, status)
        type(free_b_field_t), intent(in) :: field
        real(dp), intent(in) :: point(3) ! (R,Z,phi)
        real(dp), intent(out) :: value(EQ_VALUES)
        integer, intent(out) :: status
        real(dp) :: u, theta, xi, eta, x(2), v2(3, field%fs%sp(2)%nloc)
        real(dp) :: v0(3, field%fs%nl0), b(3), psi
        complex(dp) :: bn(3), phase
        integer :: t, n
        value = 0; status = 1
        if (.not. allocated(field%b)) return
        if (.not. all(ieee_is_finite(point))) return
        call field%coordinates(field%centre, field%radial_scale, point(:2), &
            u, theta)
        if (field%mesh%tlast(1) == 0) then
            if (u <= field%mesh%u(1)) return
        end if
        call flux_mesh_locate(field%mesh, point(:2), u, theta, t, xi, eta, status)
        if (status /= 0) return
        call meridian_point_basis(field%fs, field%pb, 2, t, xi, eta, x, v2, status)
        if (status /= 0) return
        b = 0
        do n = 0, ubound(field%b, 2)
            bn = matmul(v2, field%b(field%fs%sp(2)%dof(:, t), n) &
                *field%fs%sp(2)%sgn(:, t))
            phase = exp(cmplx(0.0_dp, n*point(3), dp))
            b = b + merge(1, 2, n == 0)*real(bn*phase, dp)
        end do
        value(EQ_BR) = b(1); value(EQ_BZ) = b(3)
        value(EQ_BPHI) = b(2) + field%f_vac/point(1)
        if (allocated(field%p)) then
            call meridian_point_basis(field%fs, field%pb, 0, t, xi, eta, x, &
                v0, status)
            if (status /= 0) return
            value(EQ_P) = nn_fourier_value(field%p, field%fs%n0, &
                field%fs%sp(0)%dof(:, t), v0(1, :), point(3))
        end if
        if (ubound(field%b, 2) == 0) then
            call meridian_point_basis(field%fs, field%pb, 0, t, xi, eta, x, &
                v0, status)
            if (status /= 0) return
            psi = dot_product(v0(1, :), field%psi(field%fs%sp(0)%dof(:, t)))
            value(EQ_PSI) = psi
            value(EQ_PSI_R) = point(1)*b(3)
            value(EQ_PSI_Z) = -point(1)*b(1)
            if (.not. allocated(field%p)) then
                call axisymmetric_pressure(psi - field%psi_bound, field%pprime, &
                    value(EQ_P))
            end if
        end if
    end subroutine

    pure real(dp) function nn_fourier_value(v, n0, dof, basis, phi) result(f)
        ! Value at phi of a V0 field in FortFEM's real Fourier block layout
        ! (block 0, then cos(n phi), sin(n phi)); size(v)/n0 blocks.
        real(dp), intent(in) :: v(:), basis(:), phi
        integer, intent(in) :: n0, dof(:)
        integer :: n
        f = dot_product(basis, v(dof))
        do n = 1, (size(v)/n0 - 1)/2
            f = f + cos(n*phi)*dot_product(basis, v((2*n - 1)*n0 + dof)) &
                + sin(n*phi)*dot_product(basis, v(2*n*n0 + dof))
        end do
    end function

    subroutine nn_velocity(field, point, velocity, status)
        ! Field-line equations in phi; unavailable where B_phi vanishes.
        type(free_b_field_t), intent(in) :: field
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: velocity(2)
        integer, intent(out) :: status
        real(dp) :: value(EQ_VALUES), background
        velocity = 0
        call nn_point(field, point, value, status)
        if (status /= 0) return
        status = 1
        if (abs(value(EQ_BPHI)) <= tiny(1.0_dp)) return
        background = field%f_vac/point(1)
        call fieldline_velocity(point(1), value(EQ_BR), value(EQ_BZ), &
            value(EQ_BPHI) - background, background, velocity(1), velocity(2))
        status = 0
    end subroutine

    subroutine circular_coordinates(centre, scale, x, u, theta)
        real(dp), intent(in) :: centre(2), scale, x(2)
        real(dp), intent(out) :: u, theta
        u = norm2(x - centre)/scale
        theta = atan2(x(2) - centre(2), x(1) - centre(1))
    end subroutine

    subroutine nn_sample(field, points, values, statuses)
        class(free_b_field_t), intent(in) :: field
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        real(dp) :: point(3)
        integer :: i
        values = 0; statuses = 1
        if (size(points, 1) < 2 .or. size(points, 1) > 3) return
        if (size(values, 1) /= EQ_VALUES) return
        if (size(values, 2) /= size(points, 2)) return
        if (size(statuses) /= size(points, 2)) return
        !$omp parallel do default(shared) private(i, point) schedule(static) &
        !$omp& if(size(points, 2) > 16)
        do i = 1, size(points, 2)
            point = 0
            point(:2) = points(:2, i)
            if (size(points, 1) == 3) point(3) = points(3, i)
            call nn_point(field, point, values(:, i), statuses(i))
        end do
        !$omp end parallel do
    end subroutine

    subroutine real_leading(z, n, a, status)
        type(csc_z_t), intent(in) :: z
        integer, intent(in) :: n
        type(csc_t), intent(out) :: a
        integer, intent(out) :: status
        integer, allocatable :: row(:), col(:)
        real(dp), allocatable :: val(:)
        type(fortsparse_status_t) :: st
        integer :: j, p, count
        allocate (row(z%nnz), col(z%nnz), val(z%nnz))
        count = 0; status = 1
        do j = 1, n
            do p = z%col_ptr(j), z%col_ptr(j + 1) - 1
                if (z%row_idx(p) > n) cycle
                if (abs(aimag(z%val(p))) > 1.0e-12_dp) return
                count = count + 1
                row(count) = z%row_idx(p); col(count) = j
                val(count) = real(z%val(p), dp)
            end do
        end do
        call csc_from_triplet(n, n, row(:count), col(:count), val(:count), a, st)
        if (status_ok(st)) status = 0
    end subroutine

    subroutine nn_linear_solve(k, fixed, rhs, x, status)
        ! Boundary increment is zero. The equation-level tangent/adjoint is
        ! this same symmetric solve, with fixed boundary data differentiated
        ! into the RHS. Never differentiate the Newton iteration.
        type(csc_t), intent(in) :: k
        logical, intent(in) :: fixed(:)
        complex(dp), intent(in) :: rhs(:)
        complex(dp), intent(out) :: x(:)
        integer, intent(out) :: status
        type(nn_factor_t) :: factor
        status = 1
        if (size(fixed) /= size(rhs) .or. size(x) /= size(rhs)) return
        call nn_factor(k, fixed, factor, status)
        if (status /= 0) return
        call nn_factor_solve(factor, rhs, x, status)
        call nn_factor_free(factor)
    end subroutine

    subroutine nn_factor(k, fixed, factor, status)
        ! Symmetric positive definite factor of k with identity rows and
        ! columns at the fixed dofs; reused across right-hand sides.
        type(csc_t), intent(in) :: k
        logical, intent(in) :: fixed(:)
        type(nn_factor_t), intent(out) :: factor
        integer, intent(out) :: status
        type(csc_t) :: a
        type(fortsparse_status_t) :: st
        integer, allocatable :: row(:), col(:)
        real(dp), allocatable :: val(:)
        integer :: n, j, p, count
        real(dp) :: corr
        logical :: ok
        n = size(fixed); status = 1
        if (k%nrow /= n .or. k%ncol /= n) return
        allocate (row(2*(k%nnz + n)), col(2*(k%nnz + n)), val(2*(k%nnz + n)))
        count = 0
        do j = 1, n
            if (fixed(j)) then
                count = count + 1
                row(count) = j; col(count) = j; val(count) = 1
                cycle
            end if
            do p = k%col_ptr(j), k%col_ptr(j + 1) - 1
                if (fixed(k%row_idx(p))) cycle
                count = count + 1
                row(count) = k%row_idx(p); col(count) = j; val(count) = k%val(p)
                ! Explicit zero at the mirror position: assembly can store a
                ! roundoff entry on one side only.
                count = count + 1
                row(count) = j; col(count) = k%row_idx(p); val(count) = 0
            end do
        end do
        call csc_from_triplet(n, n, row(:count), col(:count), val(:count), a, st)
        if (.not. status_ok(st)) return
        ! The stiffness is symmetric only to assembly roundoff (anisotropic
        ! conductivity); the native SPD factor needs bitwise symmetry, so
        ! project explicitly and reject a real asymmetry.
        call csc_normalize_symmetric_raw(n, a%col_ptr, a%row_idx, a%val, &
            1.0e-12_dp, corr, st)
        if (.not. status_ok(st)) return
        factor%fixed = fixed
        factor%direct_backend = gs_direct_available()
        if (factor%direct_backend) then
            call gs_direct_factor(factor%direct, n, a%col_ptr, a%row_idx, a%val, ok)
            if (.not. ok) return
        else
            factor%solver%backend_id = FORTSPARSE_BACKEND_NATIVE
            call sparse_factor(factor%solver, a, st)
            if (.not. status_ok(st)) return
        end if
        factor%ready = .true.
        status = 0
    end subroutine

    subroutine nn_factor_solve(factor, rhs, x, status)
        ! Real and imaginary parts are two solves; fixed entries return zero.
        type(nn_factor_t), intent(inout) :: factor
        complex(dp), intent(in) :: rhs(:)
        complex(dp), intent(out) :: x(:)
        integer, intent(out) :: status
        real(dp), allocatable :: re(:), im(:)
        status = 1
        x = 0
        if (.not. factor%ready) return
        if (size(rhs) /= size(factor%fixed) .or. size(x) /= size(rhs)) return
        allocate (re(size(rhs)), im(size(rhs)))
        call nn_factor_solve_real(factor, real(rhs, dp), re, status)
        if (status /= 0) return
        call nn_factor_solve_real(factor, aimag(rhs), im, status)
        if (status /= 0) return
        x = cmplx(re, im, dp)
    end subroutine

    subroutine nn_factor_solve_real(factor, rhs, x, status)
        type(nn_factor_t), intent(inout) :: factor
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: status
        type(fortsparse_status_t) :: st
        real(dp), allocatable :: b(:)
        logical :: ok
        status = 1
        x = 0
        if (.not. factor%ready) return
        if (size(rhs) /= size(factor%fixed) .or. size(x) /= size(rhs)) return
        b = rhs
        where (factor%fixed) b = 0
        if (factor%direct_backend) then
            call gs_direct_solve(factor%direct, b, x, ok)
        else
            call sparse_solve(factor%solver, b, x, st)
            ok = status_ok(st)
        end if
        if (ok) status = 0
    end subroutine

    subroutine nn_factor_free(factor)
        type(nn_factor_t), intent(inout) :: factor
        if (factor%ready) then
            if (factor%direct_backend) then
                call gs_direct_free(factor%direct)
            else
                call sparse_destroy(factor%solver)
            end if
        end if
        factor%ready = .false.
    end subroutine
end module kin6d_equilibrium_nonnested
