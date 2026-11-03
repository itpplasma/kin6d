program test_vrmhd_harmonics
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_gs_symmetry, only: gs_symmetry_t, GS_ROTATION, GS_TRANSLATION, &
        GS_CHART_RZ, GS_CHART_CARTESIAN
    use kin6d_vrmhd_parameters, only: vrmhd_parameters_t
    use kin6d_vrmhd_point, only: vrmhd_point
    use kin6d_vrmhd_section
    use kin6d_vrmhd_fourier
    use vrmhd_mms_fourier_leaves, only: mms_fourier_rotation, &
        mms_fourier_translation
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    use fortfem_meridian_feec_n, only: meridian_qpoint, meridian_lagrange_nodes
    use fortsparse, only: csc_t, csc_matvec, csc_from_triplet, sparse_solver_t, &
        sparse_factor, sparse_solve, sparse_destroy, FORTSPARSE_SYM_SPD, &
        fortsparse_status_t, status_ok
    use vrmhd_spectrum_support, only: spectrum, block_dofs
    use rm_mhd_oracles, only: uniform_waves, oracle_mu0, oracle_ok
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(24) :: gate
    integer :: kind, degree
    logical :: mixed_family
    call get_command_argument(1, gate)
    mixed_family = index(gate, 'mixed_') == 1
    if (mixed_family) gate = gate(7:)
    if (gate == 'spectra' .or. gate == 'spectra_fine' .or. &
        gate == 'spectra_equal' .or. gate == 'spectra_equal_fine') then
        call spectra(index(gate, 'fine') /= 0, index(gate, 'equal') == 0)
        stop
    end if
    if (gate == 'accumulation') then
        call accumulation()
        stop
    end if
    do kind = GS_ROTATION, GS_TRANSLATION
        do degree = 1, 2
            select case (gate)
            case ('exactness')
                call exactness(kind, degree)
            case ('reduction')
                call reduction(kind, degree, .false.)
                call reduction(kind, degree, .true.)
            case ('decoupling')
                call decoupling(kind, degree)
            case ('manufactured')
                call manufactured(kind, degree)
            case default
                error stop 'unknown RM3 gate'
            end select
        end do
    end do
contains
    subroutine build_fourier(mm, symmetry, degree, harmonics, scale, op, status)
        type(meridian_mesh_t), intent(in) :: mm
        type(gs_symmetry_t), intent(in) :: symmetry
        integer, intent(in) :: degree, harmonics(:)
        real(dp), intent(in) :: scale
        type(vrmhd_fourier_t), intent(out) :: op
        integer, intent(out) :: status
        if (mixed_family) then
            call vrmhd_fourier_build(mm, symmetry, degree + 1, harmonics, scale, &
                op, status, &
                scalar_degree=degree, velocity_degree=degree + 1)
        else
            call vrmhd_fourier_build(mm, symmetry, degree, harmonics, scale, op, status)
        end if
    end subroutine

    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print '(a)', message
            error stop 'RM3 gate failed'
        end if
    end subroutine

    subroutine mesh_for(kind, cells, x0, mm, symmetry, scale)
        integer, intent(in) :: kind, cells
        real(dp), intent(in) :: x0(2)
        type(meridian_mesh_t), intent(out) :: mm
        type(gs_symmetry_t), intent(out) :: symmetry
        real(dp), intent(out) :: scale
        type(mesh_2d_t) :: mesh
        integer :: status
        real(dp) :: shift
        shift = 0; scale = 0.7_dp
        symmetry = gs_symmetry_t(kind, GS_CHART_CARTESIAN, 0.0_dp)
        if (kind == GS_ROTATION) then
            shift = 2; scale = 1; symmetry%chart = GS_CHART_RZ
        end if
        call mesh%create_rectangular(cells + 1, cells + 1, x0(1) + shift, &
            x0(1) + shift + 1, x0(2), x0(2) + 1)
        call meridian_mesh_straight(mesh, 1, mm, status)
        call require(status == 0, 'section mesh')
    end subroutine

    subroutine parameters(p, limited)
        type(vrmhd_parameters_t), intent(out) :: p
        logical, intent(in) :: limited
        p%mi = 1.2_dp; p%mu0 = 0.8_dp; p%eta0 = 0.2_dp; p%nu0 = 0.3_dp
        p%chi_parallel0 = [0.7_dp, 0.4_dp]; p%chi_perpendicular = [0.05_dp, 0.03_dp]
        p%d_parallel = 0.2_dp; p%d_perpendicular = 0.02_dp
        p%nu_parallel = 0.1_dp; p%equilibration_rate = 0.4_dp; p%epsilon_b = 0.1_dp
        p%flux_limited = limited; p%electron_mass = 0.04_dp
        p%gradient_floor = 0.03_dp
    end subroutine

    subroutine state(op, u, ud)
        type(vrmhd_section_t), intent(in) :: op
        real(dp), intent(out) :: u(:), ud(:)
        integer :: i, n0
        n0 = op%scalar_fs%n0
        do i = 1, size(u)
            u(i) = 0.05_dp*sin(real(i, dp)); ud(i) = 0.02_dp*cos(real(i, dp))
        end do
        u(:n0) = u(:n0) + 2
        u(op%offset(5) + 1:op%offset(7)) = &
            u(op%offset(5) + 1:op%offset(7)) + 1.5_dp
    end subroutine

    ! (a + i kappa b) x for a FortFEM per-harmonic operator; abs gives the
    ! roundoff scale |a||x| + |kappa||b||x|.
    function apply(a, b, kappa, x, absolute) result(y)
        type(csc_t), intent(in) :: a, b
        real(dp), intent(in) :: kappa
        complex(dp), intent(in) :: x(:)
        logical, intent(in) :: absolute
        complex(dp) :: y(a%nrow)
        if (absolute) then
            y = cmplx(magnitude(a, abs(x)) + abs(kappa)*magnitude(b, abs(x)), &
                0.0_dp, dp)
        else
            y = cmplx(csc_matvec(a, real(x)), csc_matvec(a, aimag(x)), dp) &
                + cmplx(0.0_dp, kappa, dp)*cmplx(csc_matvec(b, real(x)), &
                csc_matvec(b, aimag(x)), dp)
        end if
    end function

    function magnitude(a, x) result(y)
        type(csc_t), intent(in) :: a
        real(dp), intent(in) :: x(:)
        real(dp) :: y(a%nrow)
        integer :: j, k
        y = 0
        do j = 1, a%ncol
            do k = a%col_ptr(j), a%col_ptr(j + 1) - 1
                y(a%row_idx(k)) = y(a%row_idx(k)) + abs(a%val(k))*x(j)
            end do
        end do
    end function

    ! Harmonic data at the section: cos and sin coefficients of all 36 slots.
    subroutine harmonic_data(basis, orbit, kappa, uc, us, dc, ds)
        real(dp), intent(in) :: basis(:, :), orbit(:, :), kappa, uc(:), us(:)
        real(dp), intent(out) :: dc(36), ds(36)
        dc = matmul(basis, uc) + kappa*matmul(orbit, us)
        ds = matmul(basis, us) - kappa*matmul(orbit, uc)
    end subroutine

    function polynomial(xy, degree) result(f)
        real(dp), intent(in) :: xy(:, :)
        integer, intent(in) :: degree
        complex(dp) :: f(size(xy, 2))
        real(dp) :: d
        d = degree - 1
        f = cmplx(0.3_dp + 0.7_dp*xy(1, :) - 0.4_dp*xy(2, :) &
            + d*(0.5_dp*xy(1, :)**2 - 0.2_dp*xy(1, :)*xy(2, :) + 0.6_dp*xy(2, :)**2), &
            -0.2_dp + 0.1_dp*xy(1, :) + 0.9_dp*xy(2, :) &
            + d*(-0.3_dp*xy(1, :)**2 + 0.8_dp*xy(1, :)*xy(2, :) - 0.1_dp*xy(2, &
            :)**2), dp)
    end function

    subroutine exactness(kind, degree)
        integer, intent(in) :: kind, degree
        type(vrmhd_fourier_t) :: op
        type(meridian_mesh_t) :: mm
        type(gs_symmetry_t) :: symmetry
        complex(dp), allocatable :: f(:), g(:), c(:), a(:), ca(:), dca(:), f0(:), a0(:)
        complex(dp), allocatable :: fp(:), fr(:), gp(:)
        real(dp), allocatable :: uc(:), us(:), basis(:, :), orbit(:, :), psi(:, :)
        real(dp), allocatable :: xa(:, :), xs(:, :)
        real(dp) :: scale, kappa, x(2), weight, dc(36), ds(36), ac(36), as(36)
        real(dp) :: worst(4), size_grad, size_b, sign, size_curl, curl_error
        complex(dp) :: bff(3), field(3)
        integer :: b, i, l, t, q, st, n0, n1, nl, ns, nr
        call mesh_for(kind, 3, [-0.4_dp, -0.3_dp], mm, symmetry, scale)
        call build_fourier(mm, symmetry, degree, [0, 1, 2], scale, op, st)
        call require(st == 0, 'Fourier build')
        associate (fs => op%section%fs)
            n0 = fs%n0; n1 = fs%sp(1)%ndof; nl = op%section%nlocal
            ns = op%section%ndof; nr = op%section%scalar_fs%n0
            ! rho and the edge complex's Lagrange space differ for second-kind
            ! edges; a polynomial of the scalar degree lies in both exactly.
            call meridian_lagrange_nodes(fs, xa)
            call meridian_lagrange_nodes(op%section%scalar_fs, xs)
            fp = polynomial(xa, degree); fr = polynomial(xs, degree)
            allocate (f(n0), a(n1), uc(nl), us(nl), basis(36, nl), orbit(36, nl), &
                psi(3, fs%sp(2)%nloc))
            do i = 1, n0
                f(i) = cmplx(sin(1.3_dp*i), cos(0.7_dp*i), dp)
            end do
            do i = 1, n1
                a(i) = cmplx(cos(1.1_dp*i), sin(0.9_dp*i), dp)
            end do
            f0 = f; a0 = a
            ! FortFEM's straight-cylinder slots (x, z, y) are left-handed.
            sign = 1
            if (kind == GS_TRANSLATION) sign = -1
            do b = 1, size(op%harmonic)
                kappa = vrmhd_fourier_wavenumber(op, b)
                f = f0; a = a0
                ! An invariant field has no sin coefficient.
                if (op%harmonic(b) == 0) then
                    f = real(f); a = real(a)
                end if
                g = apply(fs%grad%a, fs%grad%b, kappa, f, .false.)
                gp = apply(fs%grad%a, fs%grad%b, kappa, fp, .false.)
                if (op%harmonic(b) == 0) then
                    gp = apply(fs%grad%a, fs%grad%b, kappa, &
                        cmplx(real(fp), 0.0_dp, dp), .false.)
                end if
                c = apply(fs%curl%a, fs%curl%b, kappa, g, .false.)
                worst(1) = maxval(abs(c))/maxval(abs(apply(fs%curl%a, fs%curl%b, &
                    kappa, g, .true.)))
                ca = apply(fs%curl%a, fs%curl%b, kappa, a, .false.)
                dca = apply(fs%div%a, fs%div%b, kappa, ca, .false.)
                worst(2) = maxval(abs(dca))/maxval(abs(apply(fs%div%a, fs%div%b, &
                    kappa, ca, .true.)))
                worst(3:4) = 0; size_grad = 0; size_b = 0
                size_curl = 0; curl_error = 0
                do t = 1, fs%mm%topology%n_triangles
                    do q = 1, fs%nq
                        call vrmhd_fourier_data(op, t, q, x, weight, basis, orbit, st)
                        call require(st == 0, 'Fourier data')
                        ! rho = f, A = J = grad_n f; then a second A = a alone.
                        do l = 1, nl
                            i = op%section%map(l, t)
                            if (i <= nr) then
                                uc(l) = real(fr(i)); us(l) = -aimag(fr(i))
                                if (op%harmonic(b) == 0) us(l) = 0
                            else if (i > op%section%offset(7) .and. i <= &
                                    op%section%offset(7) + n1) then
                                uc(l) = real(gp(i - op%section%offset(7))); &
                                    us(l) = -aimag(gp(i - op%section%offset(7)))
                            else if (i > op%section%offset(7) + n1) then
                                uc(l) = real(gp(i - op%section%offset(7) - n1))
                                us(l) = -aimag(gp(i - op%section%offset(7) - n1))
                            else
                                uc(l) = 0; us(l) = 0
                            end if
                        end do
                        call harmonic_data(basis, orbit, kappa, uc, us, dc, ds)
                        size_grad = max(size_grad, maxval(abs(dc(13:15))), &
                            maxval(abs(ds(13:15))))
                        worst(3) = max(worst(3), maxval(abs(dc(5:7) - dc(13:15))), &
                            maxval(abs(ds(5:7) - ds(13:15))), &
                            maxval(abs(dc(8:10) - dc(13:15))), &
                            maxval(abs(ds(8:10) - ds(13:15))))
                        ! Pointwise curl of a gradient is roundoff; second-kind
                        ! basis values carry the conditioning of their moment
                        ! inverse, so it is measured against |curl phi| |coeff|.
                        size_curl = max(size_curl, maxval( &
                            matmul(abs(basis(25:27, :)), abs(uc)) &
                            + abs(kappa)*matmul(abs(orbit(25:27, :)), abs(us))), &
                            maxval(matmul(abs(basis(25:27, :)), abs(us)) &
                            + abs(kappa)*matmul(abs(orbit(25:27, :)), abs(uc))))
                        curl_error = max(curl_error, &
                            maxval(abs(dc(25:27))), maxval(abs(ds(25:27))), &
                            maxval(abs(dc(34:36))), maxval(abs(ds(34:36))))
                        do l = 1, nl
                            i = op%section%map(l, t) - op%section%offset(7)
                            uc(l) = 0; us(l) = 0
                            if (i >= 1 .and. i <= n1) then
                                uc(l) = real(a(i)); us(l) = -aimag(a(i))
                            end if
                        end do
                        call harmonic_data(basis, orbit, kappa, uc, us, ac, as)
                        call meridian_qpoint(fs, 2, t, q, x, weight, psi, st)
                        bff = 0
                        do l = 1, fs%sp(2)%nloc
                            bff = bff + psi(:, l)*real(fs%sp(2)%sgn(l, t), dp) &
                                *ca(fs%sp(2)%dof(l, t))
                        end do
                        field = bff
                        if (kind == GS_TRANSLATION) field = [bff(1), bff(3), bff(2)]
                        field = sign*field
                        size_b = max(size_b, maxval(abs(field)))
                        worst(4) = max(worst(4), &
                            maxval(abs(ac(25:27) - real(field))), &
                            maxval(abs(as(25:27) + aimag(field))))
                    end do
                end do
                worst(3) = max(worst(3)/size_grad, curl_error/size_curl)
                worst(4) = worst(4)/size_b
                print '(a,i0,a,i0,a,i0,a,4es11.3)', 'exactness symmetry ', kind, &
                    ' degree ', degree, ' n ', op%harmonic(b), &
                    ': curl grad, div curl, grad/curl data, curl vs FortFEM ', worst
                call require(all(worst < 1.0e-13_dp), &
                    'discrete sequence identities per harmonic')
            end do
        end associate
    end subroutine

    subroutine reduction(kind, degree, limited)
        integer, intent(in) :: kind, degree
        logical, intent(in) :: limited
        type(vrmhd_fourier_t) :: op, multi
        type(vrmhd_section_t) :: section
        type(vrmhd_parameters_t) :: p
        type(meridian_mesh_t) :: mm
        type(gs_symmetry_t) :: symmetry
        type(csc_t) :: jf, js, jm
        real(dp), allocatable :: u(:), ud(:), du(:), dud(:), r(:), jv(:), rs(:), &
            js_v(:), source(:, :, :), big(:), rm(:), jvm(:), dbig(:)
        real(dp) :: scale, shift, worst(3)
        integer :: st, i, n
        call mesh_for(kind, 2, [-0.4_dp, -0.3_dp], mm, symmetry, scale)
        call build_fourier(mm, symmetry, degree, [0], scale, op, st)
        call require(st == 0, 'n = 0 build')
        call build_fourier(mm, symmetry, degree, [2, 0, 1], scale, multi, st)
        call require(st == 0, 'multi-harmonic build')
        if (mixed_family) then
            call vrmhd_section_build(mm, symmetry, degree + 1, 1, section, st, &
                scalar_degree=degree, velocity_degree=degree + 1)
        else
            call vrmhd_section_build(mm, symmetry, degree, 1, section, st)
        end if
        call require(st == 0, 'section build')
        call parameters(p, limited)
        n = section%ndof
        call require(op%ndof == n .and. multi%ndof == 5*n, 'harmonic layout')
        allocate (u(n), ud(n), du(n), dud(n), r(n), jv(n), rs(n), js_v(n), &
            source(6, section%fs%nq, section%fs%mm%topology%n_triangles), &
            big(5*n), dbig(5*n), rm(5*n), jvm(5*n))
        call state(section, u, ud)
        shift = 0.31_dp
        do i = 1, n
            du(i) = 0.1_dp*cos(real(2*i, dp))
        end do
        dud = shift*du
        source(:3, :, :) = 0.13_dp; source(4:6, :, :) = 0.07_dp
        call vrmhd_section_residual(section, p, u, ud, du, dud, source, rs, js_v, st)
        call require(st == 0, 'section residual')
        call vrmhd_fourier_residual(op, p, u, ud, du, dud, source, r, jv, st)
        call require(st == 0, 'Fourier residual')
        worst(1) = max(maxval(abs(r - rs))/maxval(abs(rs)), &
            maxval(abs(jv - js_v))/maxval(abs(js_v)))
        call vrmhd_section_jacobian(section, p, u, ud, source, shift, js, st)
        call require(st == 0, 'section Jacobian')
        call vrmhd_fourier_jacobian(op, p, u, ud, source, shift, jf, st)
        call require(st == 0, 'Fourier Jacobian')
        rs = csc_matvec(js, du)
        worst(2) = maxval(abs(csc_matvec(jf, du) - rs))/maxval(abs(rs))
        ! The n = 0 block of a multi-harmonic layout is the section operator.
        big = 0; big(2*n + 1:3*n) = du; dbig = shift*big
        call vrmhd_fourier_jacobian(multi, p, u, ud, source, shift, jm, st)
        call require(st == 0, 'multi-harmonic Jacobian')
        call vrmhd_fourier_residual(multi, p, u, ud, big, dbig, source, rm, jvm, st)
        call require(st == 0, 'multi-harmonic residual')
        jvm = csc_matvec(jm, big)
        worst(3) = max(maxval(abs(jvm(2*n + 1:3*n) - rs))/maxval(abs(rs)), &
            maxval(abs(jvm(:2*n)))/maxval(abs(rs)), &
            maxval(abs(jvm(3*n + 1:)))/maxval(abs(rs)), &
            maxval(abs(rm(:2*n))), maxval(abs(rm(3*n + 1:))))
        print '(a,i0,a,i0,l2,a,3es11.3)', 'reduction symmetry ', kind, ' degree ', &
            degree, limited, ': residual, Jacobian, multi-harmonic n=0 ', worst
        ! Equal up to the order of OpenMP atomic accumulation.
        call require(worst(1) < 1.0e-14_dp, 'n = 0 residual is the section residual')
        call require(all(worst(2:3) < 1.0e-12_dp), &
            'n = 0 Jacobian reproduces the section Jacobian')
    end subroutine

    ! Unit-test, unit-trial tangent matrix of the point weak form.
    subroutine point_matrix(p, v, vd, source, shift, kmat)
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: v(33), vd(12), source(6), shift
        real(dp), intent(out) :: kmat(19, 45, 6)
        real(dp) :: du(33), dud(12), test(19), value
        integer :: row, i, j
        do row = 1, 6
            do i = 1, 19
                test = 0; test(i) = 1
                do j = 1, 45
                    du = 0; dud = 0
                    if (j <= 33) du(j) = 1
                    if (j > 33) dud(j - 33) = shift
                    call vrmhd_point(p, row, v, du, vd, dud, test, source, value, &
                        kmat(i, j, row))
                end do
            end do
        end do
    end subroutine

    ! Dense Jacobian from orbit samples of cos/sin products, no averaging
    ! identities: the coupling between harmonics must vanish by quadrature.
    subroutine decoupling(kind, degree)
        integer, intent(in) :: kind, degree
        integer, parameter :: samples = 12
        type(vrmhd_fourier_t) :: op
        type(vrmhd_parameters_t) :: p
        type(meridian_mesh_t) :: mm
        type(gs_symmetry_t) :: symmetry
        type(csc_t) :: jac
        real(dp), allocatable :: u(:), ud(:), source(:, :, :), dense(:, :), e(:)
        real(dp), allocatable :: basis(:, :), orbit(:, :), data(:, :), tests(:, :)
        real(dp), allocatable :: trial(:, :), kd(:, :, :), du(:), dud(:), r(:), &
            jv(:)
        integer, allocatable :: glob(:), owner(:)
        real(dp) :: scale, shift, x(2), weight, v(36), vd(12), kmat(19, 45, 6)
        real(dp) :: tau, kappa, cs, sn, coupled, diagonal, mismatch, apply_error
        integer :: st, nl, ns, nf, t, q, m, b, c, j, i, row, f, nt
        call mesh_for(kind, 2, [-0.4_dp, -0.3_dp], mm, symmetry, scale)
        call build_fourier(mm, symmetry, degree, [0, 1, 2], scale, op, st)
        call require(st == 0, 'Fourier build')
        call parameters(p, .false.)
        nl = op%section%nlocal; ns = op%section%ndof; nf = 5*nl
        nt = op%section%fs%mm%topology%n_triangles
        allocate (u(ns), ud(ns), source(6, op%section%fs%nq, nt), &
            dense(op%ndof, op%ndof), e(op%ndof), basis(36, nl), orbit(36, nl), &
            data(36, nf), tests(19, nf), trial(45, nf), kd(19, nf, 6), &
            glob(nf), owner(op%ndof), du(op%ndof), dud(op%ndof), r(op%ndof), &
            jv(op%ndof))
        call state(op%section, u, ud)
        shift = 0.31_dp
        source(:3, :, :) = 0.13_dp; source(4:6, :, :) = 0.07_dp
        do b = 1, size(op%harmonic)
            owner(op%offset(b) + 1:op%offset(b) + op%copies(b)*ns) = b
        end do
        dense = 0
        do t = 1, nt
            do q = 1, op%section%fs%nq
                call vrmhd_fourier_data(op, t, q, x, weight, basis, orbit, st)
                call require(st == 0, 'Fourier data')
                v = matmul(basis, u(op%section%map(:, t)))
                vd = matmul(basis(:12, :), ud(op%section%map(:, t)))
                call point_matrix(p, v(:33), vd, source(:, q, t), shift, kmat)
                do m = 1, samples
                    tau = 2*pi*op%orbit_scale*(m - 1)/samples
                    f = 0
                    do b = 1, size(op%harmonic)
                        kappa = vrmhd_fourier_wavenumber(op, b)
                        cs = cos(kappa*tau); sn = sin(kappa*tau)
                        do c = 1, op%copies(b)
                            do j = 1, nl
                                f = f + 1
                                glob(f) = op%offset(b) + (c - 1)*ns &
                                    + op%section%map(j, t)
                                if (c == 1) then
                                    data(:, f) = cs*basis(:, j) - kappa*sn*orbit(:, j)
                                else
                                    data(:, f) = sn*basis(:, j) + kappa*cs*orbit(:, j)
                                end if
                                call vrmhd_section_test(op%section%row(j), &
                                    data(:, f), tests(:, f))
                            end do
                        end do
                    end do
                    trial(:33, :) = data(:33, :); trial(34:, :) = data(:12, :)
                    do row = 1, 6
                        kd(:, :, row) = matmul(kmat(:, :, row), trial)
                    end do
                    do f = 1, nf
                        row = op%section%row(mod(f - 1, nl) + 1)
                        dense(glob(f), glob) = dense(glob(f), glob) &
                            + weight/samples*matmul(tests(:, f), kd(:, :, row))
                    end do
                end do
            end do
        end do
        call vrmhd_fourier_jacobian(op, p, u, ud, source, shift, jac, st)
        call require(st == 0, 'Fourier Jacobian')
        coupled = 0; diagonal = maxval(abs(dense)); mismatch = 0
        do j = 1, op%ndof
            e = 0; e(j) = 1
            e = csc_matvec(jac, e)
            mismatch = max(mismatch, maxval(abs(e - dense(:, j))))
            do i = 1, op%ndof
                if (owner(i) /= owner(j)) coupled = max(coupled, abs(dense(i, j)))
            end do
        end do
        do i = 1, op%ndof
            du(i) = 0.1_dp*cos(real(2*i, dp))
        end do
        dud = shift*du
        call vrmhd_fourier_residual(op, p, u, ud, du, dud, source, r, jv, st)
        call require(st == 0, 'matrix-free Fourier residual')
        e = csc_matvec(jac, du)
        apply_error = maxval(abs(jv - e))/maxval(abs(e))
        print '(a,i0,a,i0,a,3es11.3)', 'decoupling symmetry ', kind, ' degree ', &
            degree, ': sampled coupling, sparse vs sampled, sparse vs matrix-free ', &
            coupled/diagonal, mismatch/diagonal, apply_error
        call require(coupled/diagonal < 1.0e-13_dp, &
            'harmonics decouple about an invariant state')
        call require(mismatch/diagonal < 1.0e-12_dp, &
            'block Jacobian equals the orbit-sampled Jacobian')
        call require(apply_error < 1.0e-12_dp, 'sparse versus matrix-free JVP')
    end subroutine

    subroutine exact(kind, x, scale, c0, s0, out)
        integer, intent(in) :: kind
        real(dp), intent(in) :: x(2), scale, c0, s0
        real(dp), intent(out) :: out(84)
        if (kind == GS_ROTATION) then
            call mms_fourier_rotation(x(1), x(2), scale, c0, s0, out)
        else
            call mms_fourier_translation(x(1), x(2), scale, c0, s0, out)
        end if
    end subroutine

    pure function norms(u) result(e)
        real(dp), intent(in) :: u(36)
        real(dp) :: e(8)
        integer :: i
        e(1) = u(1)**2 + sum(u(13:15)**2)
        do i = 2, 4
            e(i) = u(i)**2 + u(i + 14)**2 + u(i + 17)**2 + u(i + 20)**2
        end do
        e(5) = sum(u(5:7)**2) + sum(u(25:27)**2)
        e(6) = sum(u(8:10)**2) + sum(u(34:36)**2)
        e(7) = u(11)**2 + sum(u(28:30)**2)
        e(8) = u(12)**2 + sum(u(31:33)**2)
    end function

    subroutine manufactured(kind, degree)
        integer, intent(in) :: kind, degree
        real(dp) :: errors(8, 3), defect(3), rates(8), rate, target(8)
        integer :: level, cells
        do level = 1, 3
            cells = 2**(level + 3 - degree)
            call manufactured_level(kind, degree, cells, errors(:, level), &
                defect(level))
            print '(a,i0,a,i0,a,i0,a,8es10.3,a,es10.3)', 'manufactured symmetry ', &
                kind, ' degree ', degree, ' cells ', cells, ': errors ', &
                errors(:, level), ' defect ', defect(level)
        end do
        rates = log(errors(:, 2)/errors(:, 3))/log(2.0_dp)
        rate = log(defect(2)/defect(3))/log(2.0_dp)
        print '(a,8f7.3,a,f7.3)', 'finest-pair natural-norm rates ', rates, &
            ' defect rate ', rate
        target = degree
        if (mixed_family) target(2:6) = degree + 1
        call require(all(rates >= target - 0.2_dp), &
            'n = 1 natural-norm design order')
        call require(rate >= degree - 0.2_dp, 'n = 1 linearised weak-defect order')
    end subroutine

    subroutine manufactured_level(kind, degree, cells, error, defect)
        integer, intent(in) :: kind, degree, cells
        real(dp), intent(out) :: error(8), defect
        type(vrmhd_fourier_t) :: op
        type(vrmhd_parameters_t) :: p
        type(meridian_mesh_t) :: mm
        type(gs_symmetry_t) :: symmetry
        type(csc_t) :: mass0, mass1
        type(fortsparse_status_t) :: sparse_status
        real(dp), allocatable :: basis(:, :), orbit(:, :), lift(:, :), u0(:), &
            rhs0(:), du(:), rhs1(:), zero(:), r(:), jv(:), load(:), diagonal(:), &
            source(:, :, :), b0(:, :), b1(:, :), vals0(:), vals1(:)
        integer, allocatable :: rows0(:), cols0(:), rows1(:), cols1(:), glob(:)
        logical, allocatable :: free(:)
        real(dp) :: scale, kappa, x(2), weight, outc(84), outs(84), dc(36), ds(36)
        integer :: st, nl, ns, nt, t, q, i, j, pos0, pos1
        call mesh_for(kind, cells, [0.1_dp, 0.2_dp], mm, symmetry, scale)
        call build_fourier(mm, symmetry, degree, [1], scale, op, st)
        call require(st == 0, 'manufactured build')
        p%mi = 1; p%mu0 = 1; p%eta0 = 0.2_dp; p%nu0 = 0.2_dp
        p%chi_parallel0 = 0.2_dp; p%chi_perpendicular = 0.05_dp
        p%d_parallel = 0.1_dp; p%d_perpendicular = 0.1_dp; p%epsilon_b = 0.1_dp
        kappa = vrmhd_fourier_wavenumber(op, 1)
        nl = op%section%nlocal; ns = op%section%ndof; nt = mm%topology%n_triangles
        allocate (basis(36, nl), orbit(36, nl), lift(72, 2*nl), u0(ns), rhs0(ns), &
            du(2*ns), rhs1(2*ns), zero(2*ns), r(2*ns), jv(2*ns), load(2*ns), &
            diagonal(2*ns), source(6, op%section%fs%nq, nt), b0(nl, nl), &
            b1(2*nl, 2*nl), vals0(nl*nl*nt), vals1(4*nl*nl*nt), &
            rows0(nl*nl*nt), cols0(nl*nl*nt), rows1(4*nl*nl*nt), &
            cols1(4*nl*nl*nt), glob(2*nl), free(2*ns))
        rhs0 = 0; rhs1 = 0; load = 0; diagonal = 0; zero = 0; source = 0
        pos0 = 0; pos1 = 0
        do t = 1, nt
            glob(:nl) = op%section%map(:, t); glob(nl + 1:) = ns + op%section%map(:, t)
            b0 = 0; b1 = 0
            do q = 1, op%section%fs%nq
                call vrmhd_fourier_data(op, t, q, x, weight, basis, orbit, st)
                call require(st == 0, 'manufactured data')
                call exact(kind, x, scale, 1.0_dp, 0.0_dp, outc)
                call exact(kind, x, scale, 0.0_dp, 1.0_dp, outs)
                ! Orbit-averaged natural norm of cos/sin data, factor 1/2.
                lift(:36, :nl) = basis; lift(37:, :nl) = -kappa*orbit
                lift(:36, nl + 1:) = kappa*orbit; lift(37:, nl + 1:) = basis
                call gram_update(basis, weight, b0)
                call gram_update(lift, 0.5_dp*weight, b1)
                rhs0(glob(:nl)) = rhs0(glob(:nl)) + weight*matmul(outc(:36), basis)
                rhs1(glob) = rhs1(glob) + 0.5_dp*weight &
                    *matmul([outc(37:72), outs(37:72)], lift)
                do i = 1, nl
                    load(glob(i)) = load(glob(i)) + 0.5_dp*weight &
                        *dot_product(basis(:12, i), outc(73:84))
                    load(glob(nl + i)) = load(glob(nl + i)) + 0.5_dp*weight &
                        *dot_product(basis(:12, i), outs(73:84))
                    diagonal(glob(i)) = diagonal(glob(i)) + 0.5_dp*weight &
                        *(sum(basis(:, i)**2) + kappa**2*sum(orbit(:, i)**2))
                    diagonal(glob(nl + i)) = diagonal(glob(nl + i)) + 0.5_dp*weight &
                        *(sum(basis(:, i)**2) + kappa**2*sum(orbit(:, i)**2))
                end do
            end do
            do j = 1, nl
                do i = 1, nl
                    pos0 = pos0 + 1; rows0(pos0) = glob(i); cols0(pos0) = glob(j)
                    vals0(pos0) = b0(i, j)
                end do
            end do
            do j = 1, 2*nl
                do i = 1, 2*nl
                    pos1 = pos1 + 1; rows1(pos1) = glob(i); cols1(pos1) = glob(j)
                    vals1(pos1) = b1(i, j)
                end do
            end do
        end do
        call csc_from_triplet(ns, ns, rows0, cols0, vals0, mass0, sparse_status)
        call require(status_ok(sparse_status), 'base projection mass')
        call natural_projection(mass0, rhs0, u0, sparse_status)
        call require(status_ok(sparse_status), 'base projection')
        call csc_from_triplet(2*ns, 2*ns, rows1, cols1, vals1, mass1, sparse_status)
        call require(status_ok(sparse_status), 'harmonic projection mass')
        call natural_projection(mass1, rhs1, du, sparse_status)
        call require(status_ok(sparse_status), 'harmonic projection')
        error = 0
        do t = 1, nt
            do q = 1, op%section%fs%nq
                call vrmhd_fourier_data(op, t, q, x, weight, basis, orbit, st)
                call exact(kind, x, scale, 1.0_dp, 0.0_dp, outc)
                call exact(kind, x, scale, 0.0_dp, 1.0_dp, outs)
                call harmonic_data(basis, orbit, kappa, du(op%section%map(:, t)), &
                    du(ns + op%section%map(:, t)), dc, ds)
                error = error + 0.5_dp*weight*(norms(dc - outc(37:72)) &
                    + norms(ds - outs(37:72)))
            end do
        end do
        error = sqrt(error)
        call vrmhd_fourier_residual(op, p, u0, zero(:ns), du, zero, source, r, jv, &
            st)
        call require(st == 0, 'manufactured linear residual')
        free = .not. op%fixed
        free(op%section%offset(8) + 1:ns) = .not. op%section%fs%sp(1)%bnd
        free(ns + op%section%offset(8) + 1:) = .not. op%section%fs%sp(1)%bnd
        defect = sqrt(sum((jv - load)**2/diagonal, mask=free))
    end subroutine

    subroutine gram_update(basis, weight, gram)
        real(dp), intent(in) :: basis(:,:), weight
        real(dp), intent(inout) :: gram(:,:)
        integer :: rows, cols
        interface
            subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
                import dp
                character, intent(in) :: transa, transb
                integer, intent(in) :: m, n, k, lda, ldb, ldc
                real(dp), intent(in) :: alpha, beta, a(lda,*), b(ldb,*)
                real(dp), intent(inout) :: c(ldc,*)
            end subroutine
        end interface
        rows = size(basis,1); cols = size(basis,2)
        call dgemm('T','N',cols,cols,rows,weight,basis,rows,basis,rows, &
            1.0_dp,gram,cols)
    end subroutine

    subroutine natural_projection(mass, rhs, u, status)
        type(csc_t), intent(in) :: mass
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: u(:)
        type(fortsparse_status_t), intent(out) :: status
        type(sparse_solver_t) :: solver
        real(dp) :: residual
        ! The natural-norm Gram is SPD; general LU wastes fill on its graph.
        solver%options%symmetry = FORTSPARSE_SYM_SPD
        call sparse_factor(solver, mass, status)
        if (status_ok(status)) call sparse_solve(solver, rhs, u, status)
        call sparse_destroy(solver)
        if (.not. status_ok(status)) return
        residual = sqrt(sum((csc_matvec(mass,u)-rhs)**2)/sum(rhs**2))
        print '(a,es20.12)', 'natural-norm projection relative residual: ', residual
        call require(residual < 1.0e-10_dp, 'natural-norm projection solves full Gram')
    end subroutine

    subroutine accumulation()
        !! The small Bloch patch retains the cancellation of a fine physical mesh.
        real(dp) :: b0(3), oracle(3), frequencies(3), optical
        integer :: status, unmatched
        b0 = sqrt(oracle_mu0)*[0.6_dp, 0.48_dp, 0.64_dp]
        call uniform_waves([1.0_dp, 0.0_dp, 0.5_dp], b0, 1.0_dp, 0.6_dp, &
            5.0_dp/3, oracle, status)
        call require(status == oracle_ok, 'fine harmonic uniform-wave oracle')
        oracle = sqrt(oracle)
        call spectrum(96, 2, 2, block_dofs(2, 2, VRMHD_EDGE_SECOND, 1), &
            b0, oracle, .true., frequencies, optical, unmatched, structural=.true., &
            harmonic=1, orbit_scale=2.0_dp, orientation=1, patch=.true.)
        call require(maxval(abs(frequencies/oracle - 1)) < 1.0e-7_dp, &
            'fine harmonic frequencies agree with the uniform-wave oracle')
    end subroutine

    ! RM-B2 per harmonic: uniform slab, B0 oblique, k = (1, 0, kz), kz = +-n/r0,
    ! with the structural energy gate and folded oracle matching of
    ! vrmhd_spectrum_support on every mesh.
    subroutine spectra(fine, mixed)
        logical, intent(in) :: fine, mixed
        real(dp), parameter :: r0 = 2
        integer, allocatable :: cells(:)
        real(dp) :: b0(3), oracle(3), kz, frequencies(3, 4), errors(3, 4), rates(3)
        real(dp) :: optical, condition(4), growth, refinement, full(3), patch(3)
        integer :: harmonic, orientation, mesh, ka, level, unmatched, st, nb, kv
        logical :: ordered
        ! B0 has a y component so k.B0 differs between the in-plane aliases.
        b0 = sqrt(oracle_mu0)*[0.6_dp, 0.48_dp, 0.64_dp]
        ordered = .true.
        kv = merge(2, 1, mixed)
        if (.not. fine) then
            kz = 1/r0
            call uniform_waves([1.0_dp, 0.0_dp, kz], b0, 1.0_dp, 0.6_dp, 5.0_dp/3, &
                oracle, st)
            call require(st == oracle_ok, 'oblique uniform-wave oracle')
            oracle = sqrt(oracle)
            nb = block_dofs(1, kv, VRMHD_EDGE_SECOND, 1, velocity_degree=kv)
            call spectrum(8, kv, 1, nb, b0, oracle, .true., full, optical, unmatched, &
                structural=.true., harmonic=1, orbit_scale=r0, orientation=1, &
                patch=.false., velocity_degree=kv)
            call spectrum(8, kv, 1, nb, b0, oracle, .true., patch, optical, unmatched, &
                structural=.true., harmonic=1, orbit_scale=r0, orientation=1, &
                patch=.true., velocity_degree=kv)
            print '(a,3es20.12)', 'patch / full-period frequency defect: ', &
                abs(patch/full - 1)
            call require(maxval(abs(patch/full - 1)) <= 1.0e-10_dp, &
                'patch and full-period Bloch blocks agree')
        end if
        do harmonic = 1, 2
            ! kz = -n/r0 would make k nearly perpendicular to B0.
            do orientation = 1, 1
                kz = orientation*harmonic/r0
                call uniform_waves([1.0_dp, 0.0_dp, kz], b0, 1.0_dp, 0.6_dp, &
                    5.0_dp/3, oracle, st)
                call require(st == oracle_ok, 'oblique uniform-wave oracle')
                oracle = sqrt(oracle)
                print '(a,i2,f8.4,3es20.12)', &
                    'n, kz, oracle Alfven / fast / slow [1/s]: ', &
                    harmonic, kz, oracle
                do mesh = 1, 2
                    do ka = merge(2, 1, fine), 2
                        if (fine) then
                            cells = [32, 64]
                        else if (ka == 1) then
                            cells = [4, 8, 16, 32]
                        else
                            cells = [8, 16, 32]
                        end if
                        kv = ka + merge(1, 0, mixed)
                        do level = 1, size(cells)
                            call spectrum(cells(level), kv, ka, &
                                block_dofs(ka, kv, VRMHD_EDGE_SECOND, mesh, &
                                velocity_degree=kv), b0, oracle, &
                                .true., frequencies(:, level), optical, unmatched, &
                                perturbed=mesh == 2, structural=.true., &
                                mass_condition=condition(level), harmonic=harmonic, &
                                orbit_scale=r0, orientation=orientation, patch=.true., &
                                velocity_degree=kv)
                            errors(:, level) = abs(frequencies(:, level)/oracle - 1)
                            print '(a,4i4,3es20.12)', &
                                'n, perturbed, degree, cells: frequency errors ', &
                                harmonic*orientation, mesh - 1, ka, cells(level), &
                                errors(:, level)
                            if (level == 1) cycle
                            refinement = real(cells(level), dp)/cells(level - 1)
                            growth = (condition(level)/condition(level - 1)) &
                                **(log(2.0_dp)/log(refinement))
                            call require(growth <= 2.0_dp, &
                                'bounded physical time-mass condition')
                            rates = log(errors(:, level - 1)/errors(:, level)) &
                                /log(refinement)
                            print '(a,3f12.6)', &
                                'Alfven / fast / slow frequency orders: ', &
                                rates
                        end do
                        ! Degree two reaches its design-order pair under spectra_fine.
                        if (ka == 2 .and. .not. fine) cycle
                        ! Every case runs; failures are collected for the report.
                        print '(a,4i4,3f12.6)', &
                            'n, perturbed, degree, cells: finest-pair design orders ', &
                            harmonic, mesh - 1, ka, cells(size(cells)), rates
                        if (minval(rates) <= 2*ka - 0.2_dp) then
                            print '(a)', &
                                'frequency design order gate FAILED for this case'
                            ordered = .false.
                        end if
                    end do
                end do
            end do
        end do
        call require(ordered, 'frequency design order on the finest pair')
    end subroutine
end program test_vrmhd_harmonics
