program test_nonnested_beltrami
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_build, nn_solve, &
        nn_report_t, nn_point, nn_linear_solve
    use kin6d_nonnested_beltrami, only: beltrami_system_t, beltrami_report_t, &
        beltrami_prepare, beltrami_solve, &
        beltrami_constraint, BELTRAMI_RESONANT, beltrami_mr_options_t
    use fortfem_meridian_feec_n, only: meridian_load, meridian_lagrange_nodes, &
        meridian_error_norm
    use kin6d_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BPHI, EQ_BZ
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 6.0_dp, mu = 1.1_dp
    type(flux_mesh_t) :: mesh
    type(free_b_field_t) :: field, vacuum, plus, minus
    type(beltrami_system_t) :: system
    type(beltrami_report_t) :: report, ref, other
    type(nn_report_t) :: nnreport
    real(dp), allocatable :: rings(:), data(:), source(:), a(:), da(:), ap(:), am(:)
    real(dp), allocatable :: nodes(:, :)
    complex(dp), allocatable :: load(:), bc(:), seed(:), proj(:)
    logical, allocatable :: fixed(:)
    real(dp) :: lo, hi, mid, qi, qnorm, l2(3), qscale, weight
    real(dp) :: flux, err(3), point(3), value(EQ_VALUES), r, scale, e
    complex(dp) :: bm(3), bn(3)
    integer :: st, k, lev, nr, i, j, n, mode
    logical :: detected, square, cylinder
    character(24) :: arg, kind, solver
    type(beltrami_mr_options_t) :: mr
    call get_command_argument(1, kind)
    call get_command_argument(2, arg); read (arg, *) k
    ! Optional third argument 'mr': the matrix-free MINRES path.
    call get_command_argument(3, solver)
    mr%tolerance = 1.0e-12_dp ! the factored gates below hold to this level
    square = kind == 'manufactured' .or. kind == 'lundquist'
    cylinder = kind == 'lundquist' .or. kind == 'resonance' .or. kind == 'constraints'
    mode = 0
    if (kind == 'manufactured' .or. kind == 'vacuum') mode = 1
    do lev = 1, 3
        nr = 2**lev
        if (allocated(rings)) deallocate (rings, data)
        allocate (rings(nr)); rings = [(real(i, dp)/nr, i=1, nr)]
        ! The square has exactly affine cells and a known physical boundary.
        ! This separates approximation order from curved-density mapping error.
        if (square) then
            call flux_mesh_build(rings, 8*nr, 256, 1, mapping, mesh, st)
        else
            call flux_mesh_build(rings, 8, 256, 2, mapping, mesh, st)
        end if
        if (st /= 0) error stop 'Beltrami mesh'
        if (cylinder) then
            call nn_build(mesh, k, mode, [r0, 0.0_dp], 1.0_dp, 0.0_dp, 1.0_dp, &
                field, st, &
                cylinder_r0=r0)
            flux = 2*pi*bessel_j1(mu)/mu
            if (square) flux = square_flux()
        else
            call nn_build(mesh, k, mode, [r0, 0.0_dp], 1.0_dp, 0.0_dp, 1.0_dp, &
                field, st)
            flux = 0
        end if
        if (st /= 0) error stop 'Beltrami build'
        if (square) field%coordinates => square_coordinates
        n = (2*mode + 1)*field%fs%sp(1)%ndof
        allocate (data(n)); data = 0
        if (kind == 'vacuum' .or. kind == 'lundquist') then
            call meridian_load(field%fs, 1, potential, load, st)
            if (allocated(fixed)) deallocate (fixed, proj)
            allocate (fixed(size(load)), proj(size(load))); fixed = .false.
            call nn_linear_solve(field%fs%mass(1), fixed, load, proj, st)
            if (st /= 0) error stop 'trace projection'
            call meridian_lagrange_nodes(field%fs, nodes)
            if (kind == 'vacuum') then
                data(field%fs%sp(1)%ndof + 1:field%fs%sp(1)%ndof + field%fs%nn) &
                    = 2*real(proj(:field%fs%nn), dp)
                data(field%fs%nn + 1:field%fs%sp(1)%ndof) = nodes(1, :)**2/2
            else
                data = real(proj, dp)
                do i = 1, field%fs%n0
                    r = norm2(nodes(:, i) - [r0, 0.0_dp])
                    data(field%fs%nn + i) = r0*bessel_j0(mu*r)/mu
                end do
            end if
        end if
        if (kind == 'constraints') then
            allocate (field%p(field%fs%n0)); field%p = 1; field%p(1) = 2
            call beltrami_prepare(field, flux, data, 2*mode + 3, system, st)
            if (st == 0) error stop 'accepted nonconstant Taylor pressure'
            field%p = 1
        end if
        if (solver == 'mr') then
            call beltrami_prepare(field, flux, data, 2*mode + 3, system, st, mr=mr)
        else
            call beltrami_prepare(field, flux, data, 2*mode + 3, system, st)
        end if
        if (st /= 0) error stop 'Beltrami prepare'
        if (kind == 'resonance') then
            lo = 0; hi = 8; detected = .false.
            call beltrami_solve(system, lo, field, ref)
            call beltrami_solve(system, hi, field, other)
            if (ref%status /= 0 .or. other%status /= 0) error stop 'resonance endpoints'
            if (ref%negative >= other%negative) error stop 'no inertia crossing'
            do i = 1, 50
                mid = (lo + hi)/2; plus = field
                call beltrami_solve(system, mid, field, report)
                if (report%status == BELTRAMI_RESONANT) then
                    if (maxval(abs(field%b - plus%b)) > 0) error stop &
                        'failed solve mutated B'
                    print '(a,es22.14,a,2es12.4)', 'near-resonant mu ', mid, &
                        ' eigenvalue/residual ', report%nearest_eigenvalue, &
                        report%eigen_residual
                    detected = .true.; exit
                end if
                if (report%status /= 0) error stop 'resonance bisection'
                if (report%negative == ref%negative) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            if (.not. detected) error stop 'missed resonance'
            call beltrami_constraint(system, 'helicity', 0.0_dp, &
                [mid, mid + 0.1_dp], 1.0e-8_dp, field, report)
            if (report%status /= BELTRAMI_RESONANT) &
                error stop 'constraint hid a resonant endpoint'
            if (maxval(abs(field%b - plus%b)) > 0) &
                error stop 'failed constraint mutated B'
            stop
        else if (kind == 'manufactured') then
            if (allocated(source)) deallocate (source)
            allocate (source(n)); source = 0
            do mode = 0, 1
                call meridian_load(field%fs, 1, forcing, load, st)
                if (st /= 0) error stop 'MMS source'
                if (mode == 0) then
                    source(:field%fs%sp(1)%ndof) = real(load, dp)
                else
                    source(field%fs%sp(1)%ndof + 1:2*field%fs%sp(1)%ndof) = &
                        real(load, dp)
                    source(2*field%fs%sp(1)%ndof + 1:) = -aimag(load)
                end if
            end do
            mode = 1
            call beltrami_solve(system, mu, field, report, source=source, &
                potential=a, dmu=da)
        else if (kind == 'vacuum') then
            call beltrami_solve(system, 0.0_dp, field, report)
            vacuum = field; vacuum%b = 0; vacuum%psi = 0
            bc = cmplx(data(field%fs%nn + 1:field%fs%sp(1)%ndof), 0.0_dp, dp)
            seed = 0*bc
            call nn_solve(vacuum, 0, bc, seed, nnreport, st)
            if (st /= 0) error stop 'NN1 vacuum oracle'
            bc = proj(:field%fs%nn); seed = 0*bc
            call nn_solve(vacuum, 1, bc, seed, nnreport, st)
            if (st /= 0) error stop 'NN1 Fourier vacuum oracle'
            e = maxval(abs(vacuum%b - field%b))
            print '(a,es12.4)', 'NN1 coefficient difference ', e
            if (e > 1.0e-9_dp) error stop 'mu=0 NN1 reduction'
        else
            call beltrami_solve(system, mu, field, report, potential=a, dmu=da)
        end if
        if (report%status /= 0) then
            print *, report%status, trim(report%message), report%min_pivot, &
                report%pivot_tolerance
            error stop 'Beltrami solve'
        end if
        if (report%divergence > 2.0e-11_dp) error stop 'div curl identity'
        if (kind == 'vacuum') stop
        if (lev == 1 .and. cylinder) then
            ref = report; plus = field; minus = field
            call beltrami_solve(system, mu + 1.0e-5_dp, plus, other, potential=ap)
            if (other%status /= 0) error stop 'positive tangent probe'
            call beltrami_solve(system, mu - 1.0e-5_dp, minus, other, potential=am)
            if (other%status /= 0) error stop 'negative tangent probe'
            e = norm2((ap - am)/2.0e-5_dp - da)/max(1.0_dp, norm2(da))
            print '(a,es12.4)', 'mu tangent error ', e
            if (e > 1.0e-7_dp) error stop 'equation-level mu tangent'
            if (kind == 'constraints') then
                call beltrami_constraint(system, 'current', ref%current, [0.0_dp, &
                    2.0_dp], &
                    1.0e-10_dp, plus, other)
                if (other%status /= 0 .or. abs(other%mu - mu) > 1.0e-12_dp) &
                    error stop 'current multiplier'
                call beltrami_constraint(system, 'helicity', ref%helicity, [0.0_dp, &
                    2.0_dp], &
                    1.0e-10_dp, plus, other)
                print '(a,es12.4)', 'helicity recovered mu ', other%mu
                if (other%status /= 0 .or. abs(other%mu - mu) > 1.0e-8_dp) &
                    error stop 'helicity multiplier'
                stop
            end if
        end if
        err(lev) = 0; scale = 0
        do i = 1, 16
            do j = 1, 16
                point = [r0 - 0.95_dp + 1.9_dp*(i - 0.37_dp)/16, &
                    -0.95_dp + 1.9_dp*(j - 0.23_dp)/16, 0.37_dp]
                call nn_point(field, point, value, st)
                if (st /= 0) error stop 'sample'
                mode = 0; call reference(point(:2), bm, full=.true.)
                if (kind == 'manufactured') then
                    mode = 1; call reference(point(:2), bn)
                    bm = bm + 2*real(bn*exp(cmplx(0.0_dp, point(3), dp)), dp)
                end if
                err(lev) = err(lev) + sum(([value(EQ_BR), value(EQ_BPHI), &
                    value(EQ_BZ)] &
                    - real(bm, dp))**2)
                scale = scale + sum(abs(bm)**2)
            end do
        end do
        err(lev) = sqrt(err(lev)/scale)
        l2(lev) = 0; qscale = 0
        do mode = 0, ubound(field%b, 2)
            call meridian_error_norm(field%fs, 2, field%b(:, mode), reference_mode, &
                qi, st)
            if (st /= 0) error stop 'field error quadrature'
            call meridian_error_norm(field%fs, 2, 0*field%b(:, mode), &
                reference_mode, qnorm, st)
            weight = merge(1.0_dp, 2.0_dp, mode == 0)
            l2(lev) = l2(lev) + weight*qi**2; qscale = qscale + weight*qnorm**2
        end do
        mode = ubound(field%b, 2); l2(lev) = sqrt(l2(lev)/qscale)
        print '(a,a,i0,a,i0,a,4es12.4)', trim(kind), ' P', k, ' rings=', nr, &
            ' sample/L2/residual/div ', err(lev), l2(lev), report%residual, &
            report%divergence
    end do
    e = log(err(1)/err(3))/log(4.0_dp)
    print '(a,2f9.4)', 'observed sample/L2 orders ', e, log(l2(1)/l2(3))/log(4.0_dp)
    if (e < k - 0.35_dp .or. log(l2(1)/l2(3))/log(4.0_dp) < k - 0.35_dp) &
        error stop 'Beltrami field order'
contains
    subroutine mapping(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        real(dp) :: s
        s = 1
        if (square) s = max(abs(cos(theta)), abs(sin(theta)))
        x = [r0 + u*cos(theta)/s, u*sin(theta)/s]
    end subroutine
    subroutine square_coordinates(centre, radial_scale, x, u, theta)
        real(dp), intent(in) :: centre(2), radial_scale, x(2)
        real(dp), intent(out) :: u, theta
        u = maxval(abs(x - centre))/radial_scale
        theta = atan2(x(2) - centre(2), x(1) - centre(1))
    end subroutine
    real(dp) function square_flux() result(f)
        ! Stokes integral of B/mu around the exact square, composite Simpson.
        integer, parameter :: ns = 4096
        real(dp) :: t, r, w
        integer :: i
        f = 0
        do i = 0, ns
            t = -1 + 2.0_dp*i/ns; r = sqrt(1 + t*t)
            w = 2
            if (mod(i, 2) == 1) w = 4
            if (i == 0 .or. i == ns) w = 1
            f = f + w*bessel_j1(mu*r)/(mu*r)
        end do
        f = 4*f*(2.0_dp/ns)/3
    end function
    subroutine potential(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        if (kind == 'vacuum') then
            v = [cmplx(sin(x(1))*x(2), 0.0_dp, dp), (0.0_dp, 0.0_dp), &
                cmplx(cos(x(2)), 0.0_dp, dp)]
        else
            call reference(x, v, full=.true.)
            v = v/mu
        end if
    end subroutine
    subroutine reference_mode(x, b)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: b(3)
        call reference(x, b)
    end subroutine
    subroutine reference(x, b, full)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: b(3)
        logical, optional, intent(in) :: full
        complex(dp) :: f(3)
        real(dp) :: r, theta
        if (kind == 'lundquist') then
            r = norm2(x - [r0, 0.0_dp]); theta = atan2(x(2), x(1) - r0)
            b = [cmplx(bessel_j1(mu*r)*sin(theta), 0.0_dp, dp), &
                cmplx(bessel_j0(mu*r), 0.0_dp, dp), &
                cmplx(-bessel_j1(mu*r)*cos(theta), 0.0_dp, dp)]
            if (.not. present(full)) b(2) = b(2) - field%f_vac/r0
        else
            call exact(x, b, f)
        end if
    end subroutine
    subroutine forcing(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        complex(dp) :: b(3)
        call exact(x, b, v)
    end subroutine
    subroutine exact(x, b, f)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: b(3), f(3)
        real(dp) :: y, z, p, px, pz, pxx, pzz, pxz, c
        complex(dp) :: dx(3), dz(3), curlb(3)
        y = x(1) - r0; z = x(2)
        p = (1 - y*y)**2*(1 - z*z)**2
        px = -4*y*(1 - y*y)*(1 - z*z)**2; pz = -4*z*(1 - z*z)*(1 - y*y)**2
        pxx = (12*y*y - 4)*(1 - z*z)**2; pzz = (12*z*z - 4)*(1 - y*y)**2
        pxz = 16*y*z*(1 - y*y)*(1 - z*z); c = 0.3_dp
        if (mode == 0) then
            b = [cmplx(-pz/x(1), 0.0_dp, dp), (0.0_dp, 0.0_dp), &
                cmplx(px/x(1), 0.0_dp, dp)]
            curlb = [(0.0_dp, 0.0_dp), &
                cmplx(-(pxx + pzz - px/x(1))/x(1), 0.0_dp, dp), (0.0_dp, 0.0_dp)]
        else
            b = [cmplx(0.0_dp, c*p/x(1), dp), cmplx(pz - c*px, 0.0_dp, dp), &
                cmplx(0.0_dp, -p/x(1), dp)]
            dx = [cmplx(0.0_dp, c*(px/x(1) - p/x(1)**2), dp), &
                cmplx(pxz - c*pxx, 0.0_dp, dp), &
                cmplx(0.0_dp, -px/x(1) + p/x(1)**2, dp)]
            dz = [cmplx(0.0_dp, c*pz/x(1), dp), cmplx(pzz - c*pxz, 0.0_dp, dp), &
                cmplx(0.0_dp, -pz/x(1), dp)]
            curlb = [cmplx(0.0_dp, 1.0_dp, dp)*b(3)/x(1) - dz(2), dz(1) - dx(3), &
                dx(2) + b(2)/x(1) - cmplx(0.0_dp, 1.0_dp, dp)*b(1)/x(1)]
        end if
        f = curlb - mu*b
    end subroutine
end program
