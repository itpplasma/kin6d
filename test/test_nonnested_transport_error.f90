module nt_error_fixture
    ! Fields and loads for the NN2 equilibrated transport bound.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_state, only: equilibrium_field_t, EQ_BR, EQ_BZ, &
        EQ_BPHI, EQ_VALUES
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_source_t
    use kin6d_nonnested_leaves, only: anisotropic_conductivity
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: r0 = 3, ftor = 3

    ! B = grad phi x grad psi + F grad phi about (r0, 0) with psi =
    ! ((R - r0)^2 + Z^2)/2 (1 + 2 eps cos(phi)): axisymmetric for eps = 0.
    type, extends(equilibrium_field_t) :: torus_field_t
        real(dp) :: eps = 0
    contains
        procedure :: sample => torus_sample
    end type

    ! Sovinec et al. (2004) field on the periodic-cylinder chart.
    type, extends(equilibrium_field_t) :: sovinec_field_t
    contains
        procedure :: sample => sovinec_sample
    end type

    ! Strong source of T = cos(pi x) cos(pi y) for chi_perp = 1.
    type, extends(meridian_fourier_source_t) :: sovinec_source_t
    contains
        procedure :: source => sovinec_source
    end type

    ! Weak flux load f = K grad T of T = cos(pi (R - r0)) cos(pi Z).
    type, extends(meridian_fourier_source_t) :: torus_source_t
        type(torus_field_t) :: field
        real(dp) :: chi_par = 1, chi_perp = 1
    contains
        procedure :: source => torus_source
    end type
contains
    subroutine square_map(u, theta, x)
        ! Concentric squares |x|, |y| <= u/2 about the origin.
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = 0.5_dp*u*[cos(theta), sin(theta)] &
            /max(abs(cos(theta)), abs(sin(theta)))
    end subroutine

    subroutine torus_square_map(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        call square_map(u, theta, x)
        x(1) = x(1) + r0
    end subroutine

    pure subroutine exact(x, t, g)
        ! T and its in-plane gradient at x = (x, y) of the centred square.
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: t, g(2)
        t = cos(pi*x(1))*cos(pi*x(2))
        g = -pi*[sin(pi*x(1))*cos(pi*x(2)), cos(pi*x(1))*sin(pi*x(2))]
    end subroutine

    subroutine torus_sample(field, points, values, statuses)
        class(torus_field_t), intent(in) :: field
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        real(dp) :: r, z, phi
        integer :: i
        values = 0; statuses = 1
        if (size(points, 1) /= 3) return
        do i = 1, size(points, 2)
            r = points(1, i); z = points(2, i); phi = points(3, i)
            values(EQ_BR, i) = z*(1 - 2*field%eps*cos(phi))/r
            values(EQ_BZ, i) = -(r - r0)*(1 + 2*field%eps*cos(phi))/r
            values(EQ_BPHI, i) = ftor/r
            statuses(i) = 0
        end do
    end subroutine

    subroutine sovinec_sample(field, points, values, statuses)
        class(sovinec_field_t), intent(in) :: field
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        real(dp) :: x, y
        integer :: i
        values = 0; statuses = 1
        if (size(points, 1) < 2) return
        do i = 1, size(points, 2)
            x = points(1, i); y = points(2, i)
            values(EQ_BR, i) = pi*cos(pi*x)*sin(pi*y)
            values(EQ_BZ, i) = -pi*sin(pi*x)*cos(pi*y)
            statuses(i) = 0
        end do
    end subroutine

    subroutine sovinec_source(self, t, x, f0, f, status)
        class(sovinec_source_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: f0(:), f(:, :)
        integer, intent(out) :: status
        status = 1
        f = 0
        f0 = 2*pi*pi*cos(pi*x(1, :))*cos(pi*x(2, :))
        if (t >= 1) status = 0
    end subroutine

    subroutine torus_source(self, t, x, f0, f, status)
        class(torus_source_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: f0(:), f(:, :)
        integer, intent(out) :: status
        real(dp) :: values(EQ_VALUES, size(x, 2)), k(3, 3), te, g(2)
        integer :: statuses(size(x, 2)), p
        status = 1
        f0 = 0; f = 0
        if (t < 1) return
        call self%field%sample(x, values, statuses)
        if (any(statuses /= 0)) return
        do p = 1, size(x, 2)
            call anisotropic_conductivity(values(EQ_BR, p), values(EQ_BPHI, p), &
                values(EQ_BZ, p), self%chi_par, self%chi_perp, k(1, 1), k(1, 2), &
                k(1, 3), k(2, 2), k(2, 3), k(3, 3))
            k(2, 1) = k(1, 2); k(3, 1) = k(1, 3); k(3, 2) = k(2, 3)
            call exact([x(1, p) - r0, x(2, p)], te, g)
            f(:, p) = matmul(k, [g(1), 0.0_dp, g(2)])
        end do
        status = 0
    end subroutine
end module

program test_nonnested_transport_error
    ! Equilibrated energy bound of the NN2 transport N = 0 block
    ! (nt_solve(..., estimate_error=.true.)) against exact solutions that
    ! vanish on straight-sided square domains, so the discrete Dirichlet
    ! trace is exact and the true energy error
    ! (int R grad e . K_pol grad e dR dZ)^(1/2) is computable:
    !   sovinec: Sovinec et al. (2004) on the periodic-cylinder chart,
    !     chi_par/chi_perp = 1, 1e3, 1e6, 1e9 (P2, nr = 4, 8);
    !   torus: axisymmetric tokamak-like field, weak load, ratio 1e4, modes
    !     = 1 (the N = 0 block decouples), P2 and P3;
    !   helical: phi-dependent K with modes > 0 reports mode coupling.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nt_error_fixture
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_qgrad
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_source_t
    use kin6d_nonnested_transport, only: nt_conductivity_t, nt_report_t, nt_solve
    use kin6d_nonnested_error, only: NN_ERROR_OK, NN_ERROR_MODE_COUPLING
    implicit none
    logical :: ok

    ok = .true.
    call run_sovinec(ok)
    call run_torus(ok)
    call run_helical(ok)
    if (.not. ok) error stop 1
    print '(a)', 'PASS'
contains
    subroutine build(nr, order, cylinder, fs, status)
        integer, intent(in) :: nr, order
        logical, intent(in) :: cylinder
        type(meridian_feec_t), intent(out) :: fs
        integer, intent(out) :: status
        type(flux_mesh_t) :: fm
        real(dp) :: rings(nr)
        integer :: i
        rings = [(real(i, dp)/nr, i=1, nr)]
        ! Straight (P1 geometry) concentric squares: the domain is exact.
        if (cylinder) then
            call flux_mesh_build(rings, 8, 4096, 1, square_map, fm, status)
            if (status /= 0) return
            call meridian_feec_build(fm%mm, order, fs, status, cylinder_r0=1.0_dp)
        else
            call flux_mesh_build(rings, 8, 4096, 1, torus_square_map, fm, status)
            if (status /= 0) return
            call meridian_feec_build(fm%mm, order, fs, status)
        end if
    end subroutine

    real(dp) function energy_error(fs, cond, p, shift) result(err)
        ! (int R grad e . K_pol grad e dR dZ)^(1/2), e = T - p_h (N = 0 block),
        ! K at phi = 0 (K is phi independent in the callers).
        type(meridian_feec_t), intent(in) :: fs
        type(nt_conductivity_t), intent(in) :: cond
        real(dp), intent(in) :: p(:), shift
        real(dp) :: x(2), wq, val(fs%nl0), grad(2, fs%nl0), k(3, 3, 1), te, g(2)
        real(dp) :: e(2), e2
        integer :: t, q, st
        e2 = 0
        do t = 1, fs%mm%topology%n_triangles
            do q = 1, fs%nq
                call meridian_qgrad(fs, t, q, x, wq, val, grad, st)
                call cond%tensor(t, reshape([x, 0.0_dp], [3, 1]), k, st)
                call exact([x(1) - shift, x(2)], te, g)
                e = g - matmul(grad, p(fs%sp(0)%dof(:, t)))
                e2 = e2 + wq*(k(1, 1, 1)*e(1)**2 + 2*k(1, 3, 1)*e(1)*e(2) &
                    + k(3, 3, 1)*e(2)**2)
            end do
        end do
        err = sqrt(e2)
    end function

    subroutine check_case(label, ratio, order, nr, report, error, previous, cap, &
            ok)
        ! A cap of huge(1.0_dp) only checks the guaranteed bound itself.
        character(*), intent(in) :: label
        real(dp), intent(in) :: ratio, error, cap
        integer, intent(in) :: order, nr
        type(nt_report_t), intent(in) :: report
        real(dp), intent(inout) :: previous
        logical, intent(inout) :: ok
        logical :: pass
        pass = report%error_status == NN_ERROR_OK
        pass = pass .and. report%error%bound >= error
        if (cap < huge(1.0_dp)) then
            pass = pass .and. report%error%bound < previous
            pass = pass .and. report%error%bound <= cap*error
        end if
        pass = pass .and. trim(report%error%claim) == 'practical-fixed-dirichlet-trace'
        print '(a, a, i0, a, es8.1, a, i3, 3(a, es10.3), a, f8.3, 2(1x, a))', label, &
            ' P', order, ' ratio', ratio, ' nr', nr, ' error', error, ' bound', &
            report%error%bound, ' flux', report%error%flux_defect, ' index', &
            report%error%bound/error, trim(report%error%scope), &
            merge('ok  ', 'FAIL', pass)
        previous = report%error%bound
        ok = ok .and. pass
    end subroutine

    subroutine run_sovinec(ok)
        logical, intent(inout) :: ok
        type(meridian_feec_t) :: fs
        type(sovinec_field_t), target :: field
        type(nt_conductivity_t) :: cond
        type(sovinec_source_t) :: src
        type(nt_report_t) :: report
        real(dp), allocatable :: p(:), bc(:)
        real(dp) :: previous(4), ratio
        ! Measured index, nr = 4/8/16/32: ratio 1 1.05-1.02, 1e3 3.3/5.0/6.1/
        ! 5.9. The hat-function localisation is not robust in the ratio: at
        ! 1e6 and 1e9, where the P2 solution is polluted (error ~ |||T|||),
        ! the index grows with h (7 -> 160, 7 -> 390); only the guaranteed
        ! bound is checked there.
        real(dp), parameter :: cap(4) = [1.1_dp, 8.0_dp, huge(1.0_dp), &
            huge(1.0_dp)]
        integer :: lev, ir, st
        cond%field => field
        cond%chi_perp = 1
        previous = huge(1.0_dp)
        do lev = 1, 2
            call build(4*2**(lev - 1), 2, .true., fs, st)
            if (st /= 0) error stop 'sovinec mesh'
            allocate (bc(fs%n0))
            bc = 0
            do ir = 1, 4
                ratio = 10.0_dp**(3*(ir - 1))
                cond%chi_par = ratio
                call nt_solve(fs, 0, 1, cond, src, bc, p, report, st, &
                    estimate_error=.true.)
                if (st /= 0) error stop 'sovinec solve'
                call check_case('sovinec', ratio, 2, 4*2**(lev - 1), report, &
                    energy_error(fs, cond, p, 0.0_dp), previous(ir), cap(ir), ok)
            end do
            deallocate (bc)
        end do
    end subroutine

    subroutine run_torus(ok)
        logical, intent(inout) :: ok
        integer, parameter :: modes = 1, nphi = 4
        type(meridian_feec_t) :: fs
        type(torus_field_t), target :: field
        type(nt_conductivity_t) :: cond
        type(torus_source_t) :: src
        type(nt_report_t) :: report
        real(dp), allocatable :: p(:), bc(:)
        real(dp) :: previous
        integer :: lev, st, order
        cond%field => field
        cond%chi_perp = 1
        cond%chi_par = 1.0e4_dp
        src%chi_par = cond%chi_par
        src%chi_perp = cond%chi_perp
        do order = 2, 3
            previous = huge(1.0_dp)
            do lev = 1, 2
                call build(4*2**(lev - 1), order, .false., fs, st)
                if (st /= 0) error stop 'torus mesh'
                allocate (bc(fs%n0*(2*modes + 1)))
                bc = 0
                call nt_solve(fs, modes, nphi, cond, src, bc, p, report, st, &
                    estimate_error=.true.)
                if (st /= 0) error stop 'torus solve'
                ok = ok .and. &
                    trim(report%error%scope) == 'n0-block-of-phi-independent-k'
                ! Measured index, nr = 4/8/16/32: P2 2.4/3.6/3.6/3.1, P3
                ! 1.9/2.7/3.4/3.6.
                call check_case('torus', cond%chi_par, order, 4*2**(lev - 1), &
                    report, energy_error(fs, cond, p, r0), previous, 5.0_dp, ok)
                deallocate (bc)
            end do
        end do
    end subroutine

    subroutine run_helical(ok)
        logical, intent(inout) :: ok
        integer, parameter :: modes = 2, nphi = 8
        type(meridian_feec_t) :: fs
        type(torus_field_t), target :: field
        type(nt_conductivity_t) :: cond
        type(torus_source_t) :: src
        type(nt_report_t) :: report
        real(dp), allocatable :: p(:), bc(:)
        integer :: st
        field%eps = 0.2_dp
        src%field%eps = 0.2_dp
        cond%field => field
        cond%chi_par = 1.0e2_dp
        src%chi_par = cond%chi_par
        call build(4, 2, .false., fs, st)
        if (st /= 0) error stop 'helical mesh'
        allocate (bc(fs%n0*(2*modes + 1)))
        bc = 0
        call nt_solve(fs, modes, nphi, cond, src, bc, p, report, st, &
            estimate_error=.true.)
        if (st /= 0) error stop 'helical solve'
        print '(a, i0, 1x, a)', 'helical status ', report%error_status, &
            trim(report%error%claim)
        ok = ok .and. report%error_status == NN_ERROR_MODE_COUPLING
        ok = ok .and. trim(report%error%claim) == 'unavailable-mode-coupling'
    end subroutine
end program
