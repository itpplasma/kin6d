module nn_transport_fixture
    ! Analytic fields and sources for the NN2 anisotropic transport oracles.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_state, only: equilibrium_field_t, EQ_BR, EQ_BZ, &
        EQ_BPHI, EQ_VALUES
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_source_t
    use kin6d_nonnested_leaves, only: anisotropic_conductivity
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    ! Torus: R0, minor radius a, toroidal F, helical amplitude eps.
    real(dp), parameter :: r0 = 3, minor = 1, ftor = 3, eps = 0.2_dp

    ! B = grad phi x grad psi_h + F grad phi with
    ! psi_h = ((R-R0)^2 + Z^2)/2 + eps ((R-R0)^2 - Z^2) cos(phi):
    ! exactly divergence free, non-axisymmetric (rotating-ellipse type).
    type, extends(equilibrium_field_t) :: helical_field_t
    contains
        procedure :: sample => helical_sample
    end type

    ! Sovinec et al. (2004) field on the periodic-cylinder chart: in-plane
    ! B = z x grad psi, psi = cos(pi x) cos(pi y), slots (R, phi, Z) = (x, z, y).
    type, extends(equilibrium_field_t) :: sovinec_field_t
    contains
        procedure :: sample => sovinec_sample
    end type

    ! Weak flux load f = K grad T_e of the manufactured temperature.
    type, extends(meridian_fourier_source_t) :: manufactured_source_t
        real(dp) :: chi_par = 1, chi_perp = 1
    contains
        procedure :: source => manufactured_source
    end type

    ! Strong source S = 2 pi^2 psi (chi_perp = 1) of the Sovinec test.
    type, extends(meridian_fourier_source_t) :: sovinec_source_t
    contains
        procedure :: source => sovinec_source
    end type
contains
    subroutine circle_map(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [r0 + minor*u*cos(theta), minor*u*sin(theta)]
    end subroutine

    subroutine square_map(u, theta, x)
        ! Concentric squares; the unit-label ring is the square |x|,|y| <= 1/2.
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = 0.5_dp*u*[cos(theta), sin(theta)] &
            /max(abs(cos(theta)), abs(sin(theta)))
    end subroutine

    subroutine helical_sample(field, points, values, statuses)
        class(helical_field_t), intent(in) :: field
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        real(dp) :: r, z, phi, psi_r, psi_z
        integer :: i
        values = 0; statuses = 1
        if (size(points, 1) /= 3) return
        do i = 1, size(points, 2)
            r = points(1, i); z = points(2, i); phi = points(3, i)
            psi_r = (r - r0)*(1 + 2*eps*cos(phi))
            psi_z = z*(1 - 2*eps*cos(phi))
            values(EQ_BR, i) = psi_z/r
            values(EQ_BZ, i) = -psi_r/r
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

    pure subroutine manufactured(x, t, g)
        ! T_e = T0 + T1 cos(phi) + T2 sin(2 phi); g = (d_R, (1/R) d_phi, d_Z).
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: t, g(3)
        real(dp) :: u, v, t0, t1, t2, c1, s1, c2, s2
        u = (x(1) - r0)/minor; v = x(2)/minor
        c1 = cos(x(3)); s1 = sin(x(3)); c2 = cos(2*x(3)); s2 = sin(2*x(3))
        t0 = exp(-u*u - 2*v*v)
        t1 = 0.3_dp*sin(2*u)*v
        t2 = 0.2_dp*cos(pi*v)*u
        t = t0 + t1*c1 + t2*s2
        g(1) = (-2*u*t0 + 0.6_dp*cos(2*u)*v*c1 + 0.2_dp*cos(pi*v)*s2)/minor
        g(3) = (-4*v*t0 + 0.3_dp*sin(2*u)*c1 - 0.2_dp*pi*u*sin(pi*v)*s2)/minor
        g(2) = (-t1*s1 + 2*t2*c2)/x(1)
    end subroutine

    pure real(dp) function manufactured_mode(x, n, part) result(c)
        ! Exact Fourier component of T_e at x = (R, Z).
        real(dp), intent(in) :: x(2)
        integer, intent(in) :: n, part
        real(dp) :: u, v
        u = (x(1) - r0)/minor; v = x(2)/minor
        c = 0
        if (n == 0) c = exp(-u*u - 2*v*v)
        if (n == 1 .and. part == 0) c = 0.3_dp*sin(2*u)*v
        if (n == 2 .and. part == 1) c = 0.2_dp*cos(pi*v)*u
    end function

    subroutine manufactured_source(self, t, x, f0, f, status)
        class(manufactured_source_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: f0(:), f(:, :)
        integer, intent(out) :: status
        type(helical_field_t) :: field
        real(dp) :: values(EQ_VALUES, size(x, 2)), k(3, 3), te, g(3)
        integer :: statuses(size(x, 2)), p
        status = 1
        f0 = 0; f = 0
        if (t < 1) return
        call field%sample(x, values, statuses)
        if (any(statuses /= 0)) return
        do p = 1, size(x, 2)
            call anisotropic_conductivity(values(EQ_BR, p), values(EQ_BPHI, p), &
                values(EQ_BZ, p), self%chi_par, self%chi_perp, k(1, 1), k(1, 2), &
                k(1, 3), k(2, 2), k(2, 3), k(3, 3))
            k(2, 1) = k(1, 2); k(3, 1) = k(1, 3); k(3, 2) = k(2, 3)
            call manufactured(x(:, p), te, g)
            f(:, p) = matmul(k, g)
        end do
        status = 0
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
end module

program test_nonnested_transport
    ! NN2 transport oracles for a given B (see docs/equilibrium-nonnested-nn2.md):
    !   manufactured: 3D helical torus field, weak manufactured load, chi ratio
    !     1 and 1e4, P2/P3 mesh ladders; L2 error and observed order.
    !   sovinec k: Sovinec et al. (2004) perpendicular-pollution test on the
    !     periodic-cylinder chart, chi_par/chi_perp = 1, 1e3, 1e6, 1e9, P_k
    !     mesh ladder; error |T(0,0) - 1| of the exact T = psi.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_transport_fixture
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_lagrange_nodes, meridian_qpoint
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_index, &
        meridian_fourier_value
    use kin6d_nonnested_transport, only: nt_conductivity_t, nt_report_t, nt_solve
    implicit none
    character(32) :: mode, arg
    integer :: k
    logical :: ok
    call get_command_argument(1, mode)
    call get_command_argument(2, arg)
    k = 2
    if (len_trim(arg) > 0) read (arg, *) k
    select case (trim(mode))
    case ('manufactured')
        call run_manufactured(ok)
    case ('sovinec')
        call run_sovinec(k, ok)
    case default
        ok = .false.
    end select
    if (.not. ok) error stop 1
    print '(a)', 'PASS'
contains
    subroutine build(nr, order, square, fs, status)
        integer, intent(in) :: nr, order
        logical, intent(in) :: square
        type(meridian_feec_t), intent(out) :: fs
        integer, intent(out) :: status
        type(flux_mesh_t) :: fm
        real(dp) :: rings(nr)
        integer :: i
        rings = [(real(i, dp)/nr, i=1, nr)]
        if (square) then
            ! The concentric-square rings are exactly straight-sided.
            call flux_mesh_build(rings, 8, 4096, 1, square_map, fm, status)
            if (status /= 0) return
            call meridian_feec_build(fm%mm, order, fs, status, cylinder_r0=1.0_dp)
        else
            call flux_mesh_build(rings, 8, 4096, min(order, 3), circle_map, fm, &
                status)
            if (status /= 0) return
            call meridian_feec_build(fm%mm, order, fs, status)
        end if
    end subroutine

    subroutine run_manufactured(ok)
        logical, intent(out) :: ok
        integer, parameter :: modes = 2, nphi = 8, nlev = 4
        type(meridian_feec_t) :: fs
        type(helical_field_t), target :: field
        type(nt_conductivity_t) :: cond
        type(manufactured_source_t) :: src
        type(nt_report_t) :: report
        real(dp), allocatable :: p(:), bc(:), xy(:, :)
        real(dp) :: err(nlev), interp(nlev), rate, ratio
        integer :: order, lev, st, i, n, part, ir
        ok = .true.
        cond%field => field
        do order = 2, 3
            do ir = 1, 2
                ratio = merge(1.0_dp, 1.0e4_dp, ir == 1)
                cond%chi_par = ratio; cond%chi_perp = 1
                src%chi_par = ratio; src%chi_perp = 1
                do lev = 1, nlev
                    call build(4*2**(lev - 1), order, .false., fs, st)
                    if (st /= 0) then
                        ok = .false.; return
                    end if
                    call meridian_lagrange_nodes(fs, xy)
                    allocate (bc(fs%n0*(2*modes + 1)))
                    do i = 1, fs%n0
                        do n = 0, modes
                            do part = 0, merge(0, 1, n == 0)
                                bc(meridian_fourier_index(fs%n0, n, part, i)) = &
                                    manufactured_mode(xy(:, i), n, part)
                            end do
                        end do
                    end do
                    interp(lev) = l2_error(fs, modes, bc)
                    call nt_solve(fs, modes, nphi, cond, src, bc, p, report, st)
                    deallocate (bc)
                    if (st /= 0 .or. report%residual > 1.0e-13_dp) then
                        print '(a, i0, es10.3)', 'manufactured solve failed ', st, &
                            report%residual
                        ok = .false.; return
                    end if
                    err(lev) = l2_error(fs, modes, p)
                    print '(a, i0, a, es8.1, a, i3, a, i7, a, es10.3, a, es10.3)', &
                        'manufactured P', order, ' ratio', ratio, ' nr', &
                        4*2**(lev - 1), ' ndof', report%ndof, ' L2', err(lev), &
                        ' interpolant', interp(lev)
                end do
                rate = log(err(nlev - 1)/err(nlev))/log(2.0_dp)
                print '(a, i0, a, es8.1, a, f6.2)', 'manufactured P', order, &
                    ' ratio', ratio, ' observed order', rate
                ! The polar axis fan of flux_mesh_t has a non-polynomial
                ! geometry map, limiting global L2 to O(h^3) for P3 (order
                ! 3.98 when the r < 0.3 fan is excluded); quasi-optimality
                ! against the nodal interpolant bounds the pollution.
                ok = ok .and. rate > min(order + 1, 3) - 0.35_dp
                ok = ok .and. err(nlev) <= 5*interp(nlev)
            end do
        end do
    end subroutine

    real(dp) function l2_error(fs, modes, p) result(err)
        ! Relative L2 error over Omega x [0, 2 pi) with a 16-point phi rule.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: modes
        real(dp), intent(in) :: p(:)
        real(dp) :: x(2), wq, phi(3, fs%nl0), th, te, g(3), e2, n2, ph
        integer :: t, q, j, st
        e2 = 0; n2 = 0
        !$omp parallel do default(shared) private(t, q, j, st, x, wq, phi, th, &
        !$omp& te, g, ph) reduction(+:e2, n2) schedule(static)
        do t = 1, fs%mm%topology%n_triangles
            do q = 1, fs%nq
                call meridian_qpoint(fs, 0, t, q, x, wq, phi, st)
                do j = 0, 15
                    ph = 2*pi*j/16
                    th = meridian_fourier_value(fs%n0, modes, p, &
                        fs%sp(0)%dof(:, t), phi(1, :), ph)
                    call manufactured([x, ph], te, g)
                    e2 = e2 + wq*(th - te)**2/16
                    n2 = n2 + wq*te**2/16
                end do
            end do
        end do
        !$omp end parallel do
        err = sqrt(e2/n2)
    end function

    subroutine run_sovinec(order, ok)
        integer, intent(in) :: order
        logical, intent(out) :: ok
        integer, parameter :: nlev = 5, nratio = 4
        type(meridian_feec_t) :: fs
        type(sovinec_field_t), target :: field
        type(nt_conductivity_t) :: cond
        type(sovinec_source_t) :: src
        type(nt_report_t) :: report
        real(dp), allocatable :: p(:), bc(:), xy(:, :)
        real(dp) :: err(nlev, nratio), ratio, worst
        integer :: lev, ir, st, centre(1)
        ok = .true.; worst = 0
        cond%field => field
        cond%chi_perp = 1
        do lev = 1, nlev
            call build(4*2**(lev - 1), order, .true., fs, st)
            if (st /= 0) then
                ok = .false.; return
            end if
            call meridian_lagrange_nodes(fs, xy)
            centre = minloc(norm2(xy, dim=1))
            allocate (bc(fs%n0))
            bc = 0
            do ir = 1, nratio
                ratio = 10.0_dp**(3*(ir - 1))
                cond%chi_par = ratio
                call nt_solve(fs, 0, 1, cond, src, bc, p, report, st)
                if (st /= 0) then
                    ok = .false.; return
                end if
                err(lev, ir) = abs(p(centre(1)) - 1)
                worst = max(worst, report%residual)
                print '(a, i0, a, es8.1, a, i3, a, i7, a, es10.3, a, es9.2)', &
                    'sovinec P', order, ' ratio', ratio, ' nr', 4*2**(lev - 1), &
                    ' ndof', report%ndof, ' |T(0,0)-1|', err(lev, ir), ' res', &
                    report%residual
            end do
            deallocate (bc)
        end do
        do ir = 1, nratio
            print '(a, i0, a, es8.1, a, f6.2, a, es10.3)', 'sovinec P', order, &
                ' ratio', 10.0_dp**(3*(ir - 1)), ' observed order', &
                log(err(nlev - 1, ir)/err(nlev, ir))/log(2.0_dp), ' best', &
                minval(err(:, ir))
            ok = ok .and. minval(err(:, ir)) < err(1, ir)
        end do
        ! Sovinec et al. (2004) regard errors near 1e-2 as acceptable; with
        ! affine meshes not aligned with psi, P2 meets it to 1e6 and P>=3
        ! to 1e9. The ratio-1 gate checks the isotropic discretisation.
        ok = ok .and. minval(err(:, 1)) < 1.0e-7_dp .and. worst <= 1.0e-13_dp
        if (order == 2) ok = ok .and. minval(err(:, 3)) < 1.0e-1_dp
        if (order >= 3) ok = ok .and. minval(err(:, 4)) < 1.0e-2_dp
    end subroutine
end program
