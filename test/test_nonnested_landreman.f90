module nn_landreman_case
    ! Landreman's exact 3D equilibrium, family 1 (uniform iota = 2; arXiv
    ! 2609.26742, section 2), as an NN2 case. Units mu0 = 1. In Cartesian
    ! (x, y, z) with a = sqrt(1 + eps), b = sqrt(1 - eps), X = x/a, Y = y/b,
    ! s = X^2 + Y^2 and F = sqrt(s(2 - s) - 4 z^2):
    !     B = (a(2zX - FY)/s, b(2zY + FX)/s, 1 - s),
    !     psi = (x^2 + y^2 + 4 z^2 + |B|^2 - 2 + eps^2)/4,  p = p_a - 2 psi,
    !     J = curl B (closed form below), J x B = grad p.
    ! B.grad p = 0, so K grad p = chi_perp grad p and the steady source is
    !     S = -div(K grad p) = -chi_perp div(J x B),
    ! evaluated by complex-step differentiation of the closed-form J x B
    ! (exact to rounding). With Phi = 0 the model current equals curl B iff
    ! J.b = g_par b.E, i.e. g_par = 2 pi R (J.B)/(V B_phi). J.B/B_phi > 0 on
    ! the domain, so g_par = g0 w/w_ref with w = R (J.B)/B_phi and
    ! V = 2 pi w_ref/g0. A vector potential (gauge A_z = 0) is
    !     A_x = b(z^2 Y + G X)/s + y^3/(3 b^2),
    !     A_y = -a(z^2 X - G Y)/s + x - x^3/(3 a^2),
    !     G = z F/2 + (c^2/4) asin(2z/c),  c^2 = s(2 - s),
    ! and the NN1 potential is A - F_vac (Z/R) e_R (curl of that is F_vac/R e_phi).
    ! See docs/equilibrium-nonnested-nn2.md.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), parameter :: eps = 0.02_dp, rc = 0.99_dp, rw = 0.3_dp
    real(dp), parameter :: pa = 0.25_dp, fvac = 1.0_dp, wref = 5.0_dp
    integer, parameter :: nfine = 32 ! phi points of the exact Fourier projections
    real(dp) :: chi_perp = 1, shift = 0
    integer :: block = 1
contains
    subroutine closed_bj(x, b, j)
        complex(dp), intent(in) :: x(3)
        complex(dp), intent(out) :: b(3), j(3)
        complex(dp) :: xx, yy, z, s, f, d
        real(dp) :: ca, cb
        ca = sqrt(1 + eps); cb = sqrt(1 - eps)
        xx = x(1)/ca; yy = x(2)/cb; z = x(3)
        s = xx*xx + yy*yy; f = sqrt(s*(2 - s) - 4*z*z)
        b = [ca*(2*z*xx - f*yy)/s, cb*(2*z*yy + f*xx)/s, 1 - s]
        d = xx*xx - yy*yy
        j(1) = -2*yy*(s + 1 - eps)/(cb*s) + 4*cb*xx*z/(f*s)
        j(2) = 2*xx*(s + 1 + eps)/(ca*s) + 4*ca*yy*z/(f*s)
        j(3) = -2*(s**3 - s*s - eps*s*d + 4*eps*z*z*d - 4*f*xx*yy*eps*z) &
            /(f*ca*cb*s*s)
    end subroutine

    subroutine closed(x, phi, bc, jc, p)
        ! Cylindrical (R, phi, Z) components of B and J, and p, at (R, Z, phi).
        real(dp), intent(in) :: x(2), phi
        real(dp), intent(out) :: bc(3), jc(3), p
        complex(dp) :: b(3), j(3)
        real(dp) :: c, s
        c = cos(phi); s = sin(phi)
        call closed_bj(cmplx([x(1)*c, x(1)*s, x(2)], 0.0_dp, dp), b, j)
        bc = [c*real(b(1)) + s*real(b(2)), -s*real(b(1)) + c*real(b(2)), real(b(3))]
        jc = [c*real(j(1)) + s*real(j(2)), -s*real(j(1)) + c*real(j(2)), real(j(3))]
        p = pa - (x(1)**2 + 4*x(2)**2 + sum(bc*bc) - 2 + eps*eps)/2
    end subroutine

    function potential(x, phi) result(av)
        ! NN1 potential (A_R - F_vac Z/R, A_phi, A_Z) at (R, Z, phi).
        real(dp), intent(in) :: x(2), phi
        real(dp) :: av(3), ca, cb, xc, yc, z, xx, yy, s, c2, f, g, ax, ay
        ca = sqrt(1 + eps); cb = sqrt(1 - eps)
        xc = x(1)*cos(phi); yc = x(1)*sin(phi); z = x(2)
        xx = xc/ca; yy = yc/cb; s = xx*xx + yy*yy
        c2 = s*(2 - s); f = sqrt(c2 - 4*z*z)
        g = z*f/2 + c2/4*asin(2*z/sqrt(c2))
        ax = cb*(z*z*yy + g*xx)/s + yc**3/(3*cb*cb)
        ay = -ca*(z*z*xx - g*yy)/s + xc - xc**3/(3*ca*ca)
        av = [cos(phi)*ax + sin(phi)*ay - fvac*z/x(1), &
            -sin(phi)*ax + cos(phi)*ay, 0.0_dp]
    end function

    function source3(x, phi) result(src)
        ! S = -chi_perp div(J x B) by complex steps in Cartesian.
        real(dp), intent(in) :: x(2), phi
        real(dp) :: src
        real(dp), parameter :: h = 1.0e-30_dp
        complex(dp) :: xc(3), b(3), j(3), f(3)
        integer :: i
        src = 0
        do i = 1, 3
            xc = cmplx([x(1)*cos(phi), x(1)*sin(phi), x(2)], 0.0_dp, dp)
            xc(i) = xc(i) + cmplx(0.0_dp, h, dp)
            call closed_bj(xc, b, j)
            f = [j(2)*b(3) - j(3)*b(2), j(3)*b(1) - j(1)*b(3), j(1)*b(2) - j(2)*b(1)]
            src = src + aimag(f(i))/h
        end do
        src = -chi_perp*src
    end function

    function g_shape3(x, phi) result(w)
        ! g_par/g0 = R (J.B)/(B_phi w_ref): the Ohmic profile for Phi = 0.
        real(dp), intent(in) :: x(2), phi
        real(dp) :: w, bc(3), jc(3), p
        call closed(x, phi, bc, jc, p)
        w = x(1)*dot_product(jc, bc)/(bc(2)*wref)
    end function

    real(dp) function mode_weight(bk, phi) result(c)
        ! Coefficient functional of real block bk (1: mean; 2n: cos; 2n+1: sin).
        integer, intent(in) :: bk
        real(dp), intent(in) :: phi
        if (bk == 1) then
            c = 1
        else if (mod(bk, 2) == 0) then
            c = 2*cos((bk/2)*phi)
        else
            c = 2*sin((bk/2)*phi)
        end if
    end function

    subroutine potential_block(x, v)
        ! Block `block` of the NN1 potential at the point shifted by
        ! shift (1 - u^2) in R (u the disk radius; zero shift on the wall).
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        real(dp) :: phi, acc(3), u
        integer :: j
        u = min(1.0_dp, norm2(x - [rc, 0.0_dp])/rw)
        acc = 0
        do j = 0, nfine - 1
            phi = 2*acos(-1.0_dp)*j/nfine
            acc = acc + mode_weight(block, phi) &
                *potential([x(1) - shift*(1 - u*u), x(2)], phi)
        end do
        v = cmplx(acc/nfine, 0.0_dp, dp)
    end subroutine

    real(dp) function pressure_block(x, bk) result(v)
        real(dp), intent(in) :: x(2)
        integer, intent(in) :: bk
        real(dp) :: phi, bc(3), jc(3), p, u
        integer :: j
        u = min(1.0_dp, norm2(x - [rc, 0.0_dp])/rw)
        v = 0
        do j = 0, nfine - 1
            phi = 2*acos(-1.0_dp)*j/nfine
            call closed([x(1) - shift*(1 - u*u), x(2)], phi, bc, jc, p)
            v = v + mode_weight(bk, phi)*p
        end do
        v = v/nfine
    end function

    subroutine truncated(x, phi, modes, bc, p)
        ! B and p of the exact state projected on Fourier modes 0..modes.
        real(dp), intent(in) :: x(2), phi
        integer, intent(in) :: modes
        real(dp), intent(out) :: bc(3), p
        real(dp) :: b1(3), j1(3), p1, ph, c
        integer :: j, bk
        bc = 0; p = 0
        do j = 0, nfine - 1
            ph = 2*acos(-1.0_dp)*j/nfine
            call closed(x, ph, b1, j1, p1)
            do bk = 1, 2*modes + 1
                c = mode_weight(bk, ph)*merge(1.0_dp, 0.5_dp, bk == 1)
                c = c*mode_weight(bk, phi)
                bc = bc + c*b1/nfine
                p = p + c*p1/nfine
            end do
        end do
    end subroutine

    subroutine mapping(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [rc + rw*u*cos(theta), rw*u*sin(theta)]
    end subroutine
end module

program test_nonnested_landreman
    ! Gate: the coupled NN2 solve reproduces Landreman family 1 (B and p) at
    ! the expected order on a mesh ladder from a seed that is not the target
    ! (interior shifted, V = 0, Phi = 0), with D3 falling. Errors are against
    ! the exact state projected on the retained Fourier modes 0..N; the
    ! untruncated difference is printed for the truncation ledger.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_landreman_case
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes, meridian_qpoint, &
        meridian_load
    use fortsparse, only: csc_t, csc_matvec
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_build, nn_linear_solve
    use kin6d_nonnested_residual, only: nn_residual
    use kin6d_equilibrium_residual, only: equilibrium_residual_t
    use kin6d_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI, EQ_P
    use kin6d_nonnested_coupled, only: nc_problem_t, nc_layout_t, nc_report_t, &
        nc_layout, nc_residual, nc_jacobian, nc_solve, nc_fixed
    implicit none
    type(free_b_field_t) :: field
    type(nc_problem_t) :: problem
    type(nc_layout_t) :: lay
    real(dp), allocatable :: nodes(:, :), x(:)
    integer :: k, modes, nlev, lev, nr0
    real(dp) :: ratio, v_exact, current
    real(dp), allocatable :: eb(:), ep(:), ed(:), ev(:)
    real(dp) :: rate(4)
    character(16) :: arg
    call get_command_argument(1, arg)
    read (arg, *) k
    call get_command_argument(2, arg)
    read (arg, *) ratio
    call get_command_argument(3, arg)
    read (arg, *) nlev
    modes = 2; nr0 = 4
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg)
        read (arg, *) modes
    end if
    if (k < 2 .or. k > 4 .or. ratio < 1 .or. nlev < 2) error stop 'arguments'
    allocate (eb(nlev), ep(nlev), ed(nlev), ev(nlev))
    do lev = 1, nlev
        call run(nr0*2**(lev - 1), lev == 1, eb(lev), ep(lev), ev(lev), ed(lev))
    end do
    rate = log([eb(nlev - 1)/eb(nlev), ep(nlev - 1)/ep(nlev), &
        ed(nlev - 1)/ed(nlev), ev(nlev - 1)/ev(nlev)])/log(2.0_dp)
    print '(a,4f9.3)', 'finest orders B/p/D3/V ', rate
    if (rate(1) < k - 0.3_dp) error stop 'field order'
    ! The pressure order is gated at ratio 1 only: at a large chi_par/chi_perp
    ! the Fourier truncation of b.grad p is amplified (see the NN2 doc).
    if (ratio <= 1.0_dp) then
        if (rate(2) < k + 0.6_dp) error stop 'pressure order'
    end if
    if (rate(3) < k - 0.3_dp) error stop 'D3 order'
    if (any(ed(2:) >= ed(:nlev - 1))) error stop 'D3 not falling'
    if (rate(4) < k - 0.5_dp) error stop 'loop voltage order'
    print '(a)', 'PASS: NN2 coupled Landreman family 1'
contains
    subroutine run(nr, check_jacobian, efield, epres, evolt, ed3)
        integer, intent(in) :: nr
        logical, intent(in) :: check_jacobian
        real(dp), intent(out) :: efield, epres, evolt, ed3
        type(flux_mesh_t) :: fm
        type(nc_report_t) :: report
        type(equilibrium_residual_t) :: stress
        real(dp), allocatable :: rings(:), points(:, :), xe(:)
        real(dp) :: values(EQ_VALUES, 192), bt(3), pt, bf(3), jf(3), pf
        real(dp) :: sum_b, sum_p, etrunc, phi_max, xq(2), wq, basis(3, 64), ph
        integer :: i, j, t, q, st, statuses(192)
        allocate (rings(nr), points(3, 192))
        do i = 1, nr
            rings(i) = real(i, dp)/nr
        end do
        call flux_mesh_build(rings, 8, 1024, k, mapping, fm, st)
        if (st /= 0) error stop 'mesh'
        call nn_build(fm, k, modes, [rc, 0.0_dp], rw, fvac, 1.0_dp, field, st)
        if (st /= 0) error stop 'build'
        call meridian_lagrange_nodes(field%fs, nodes)
        lay = nc_layout(field)
        ! Net toroidal current of the exact state on the discrete domain.
        current = 0
        do t = 1, field%fs%mm%topology%n_triangles
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 0, t, q, xq, wq, basis, st)
                if (st /= 0) error stop 'quadrature'
                do j = 0, nfine - 1
                    ph = 2*acos(-1.0_dp)*j/nfine
                    call closed(xq, ph, bf, jf, pf)
                    current = current + wq*jf(2)/(xq(1)*nfine)
                end do
            end do
        end do
        chi_perp = 1
        problem%chi_perp = chi_perp; problem%chi_par = ratio*chi_perp
        problem%g_par = 1; problem%g_perp = 1/ratio
        problem%source_3d => source3; problem%g_par_shape_3d => g_shape3
        problem%current = current
        v_exact = 2*acos(-1.0_dp)*wref/problem%g_par
        call seed(0.0_dp, xe)
        call seed(0.05_dp, x)
        where (nc_fixed(field)) x = xe
        if (check_jacobian) call jacobian_check(x)
        call nc_solve(field, problem, x, report, st, tol=1.0e-11_dp)
        print '(a,i0,a,i0,a,i0,a,i0,a,15es9.1)', 'P', k, ' N=', modes, ' nr=', nr, &
            ' unknowns ', lay%n, ' Newton ', report%history(:report%newton_steps)
        if (st /= 0) error stop 'coupled Newton'
        if (report%newton_steps > 15) error stop 'Newton steps'
        if (report%divergence > 1.0e-11_dp) error stop 'closed 2-form'
        if (report%multiplier > 1.0e-10_dp*maxval(abs(x(:lay%na)))) &
            error stop 'gauge multiplier'
        do i = 1, 6
            do j = 1, 32
                ph = 2*acos(-1.0_dp)*(j - 0.37_dp)/32
                call mapping(0.12_dp + 0.14_dp*i, ph, points(:2, j + 32*(i - 1)))
                points(3, j + 32*(i - 1)) = 0.43_dp + 1.7_dp*j
            end do
        end do
        call field%sample(points, values, statuses)
        if (any(statuses /= 0)) error stop 'sampling'
        efield = 0; epres = 0; sum_b = 0; sum_p = 0; etrunc = 0
        do i = 1, size(points, 2)
            call truncated(points(:2, i), points(3, i), modes, bt, pt)
            call closed(points(:2, i), points(3, i), bf, jf, pf)
            efield = efield + (values(EQ_BR, i) - bt(1))**2 &
                + (values(EQ_BPHI, i) - bt(2))**2 + (values(EQ_BZ, i) - bt(3))**2
            epres = epres + (values(EQ_P, i) - pt)**2
            etrunc = etrunc + sum((bf - bt)**2)
            sum_b = sum_b + sum(bt*bt); sum_p = sum_p + pt*pt
        end do
        efield = sqrt(efield/sum_b); epres = sqrt(epres/sum_p)
        etrunc = sqrt(etrunc/sum_b)
        evolt = abs(report%voltage - v_exact)/abs(v_exact)
        phi_max = maxval(abs(x(lay%iphi + 1:lay%iphi + lay%nv)))*2*acos(-1.0_dp) &
            /abs(v_exact)
        call nn_residual(field, 0, stress, st)
        if (st /= 0) error stop 'D3'
        ed3 = stress%normalized
        print '(a,i0,a,i0,a,4es11.3,a,es10.2,a,es10.2)', 'P', k, ' nr=', nr, &
            ' B/p/V/D3 ', efield, epres, evolt, ed3, ' Phi ', phi_max, &
            ' truncation ', etrunc
    end subroutine

    subroutine seed(dr, s)
        ! Exact NN1 potential and p, evaluated at R - dr (1 - u^2): the
        ! boundary blocks are the exact wall data (A trace in every Fourier
        ! block, p_b), Phi = 0 and V = 0.
        real(dp), intent(in) :: dr
        real(dp), allocatable, intent(out) :: s(:)
        complex(dp), allocatable :: rhs(:), a(:)
        logical, allocatable :: free(:)
        integer :: bk, n, st, n1
        allocate (s(lay%n))
        s = 0
        shift = dr
        n1 = lay%na1
        allocate (a(n1), free(n1))
        free = .false.
        do bk = 1, lay%nb
            block = bk
            call meridian_load(field%fs, 1, potential_block, rhs, st)
            if (st /= 0) error stop 'potential load'
            call nn_linear_solve(field%fs%mass(1), free, rhs, a, st)
            if (st /= 0) error stop 'potential projection'
            s(lay%ia + (bk - 1)*n1 + 1:lay%ia + bk*n1) = real(a, dp)
            do n = 1, lay%n0
                s(lay%ip + (bk - 1)*lay%n0 + n) = pressure_block(nodes(:, n), bk)
            end do
        end do
        shift = 0
    end subroutine

    subroutine jacobian_check(x0)
        ! Central difference of the residual against the assembled Jacobian
        ! (free directions only) at the 3D seed.
        real(dp), intent(in) :: x0(:)
        type(csc_t) :: jac
        real(dp), allocatable :: d(:), rp(:), rm(:), jd(:), xs(:)
        logical, allocatable :: fixed(:)
        real(dp) :: h, err
        integer :: n, st
        fixed = nc_fixed(field)
        allocate (d(lay%n))
        xs = abs(x0)
        do n = 1, lay%n
            d(n) = sin(1.3_dp*n)
        end do
        d(lay%ia + 1:lay%ia + lay%na) = d(lay%ia + 1:lay%ia + lay%na) &
            *1.0e-2_dp*maxval(xs(lay%ia + 1:lay%ia + lay%na))
        d(lay%il + 1:lay%il + lay%nv) = d(lay%il + 1:lay%il + lay%nv)*1.0e-3_dp
        d(lay%ip + 1:lay%ip + lay%nv) = d(lay%ip + 1:lay%ip + lay%nv) &
            *1.0e-2_dp*maxval(xs(lay%ip + 1:lay%ip + lay%nv))
        d(lay%iphi + 1:lay%iv) = d(lay%iphi + 1:lay%iv)*1.0e-2_dp*v_exact
        where (fixed) d = 0
        call nc_jacobian(field, problem, x0 + d, jac, st)
        if (st /= 0) error stop 'Jacobian'
        jd = csc_matvec(jac, d)
        where (fixed) jd = 0
        h = 1.0e-4_dp
        call nc_residual(field, problem, x0 + d*(1 + h), rp, st)
        if (st /= 0) error stop 'residual +'
        call nc_residual(field, problem, x0 + d*(1 - h), rm, st)
        if (st /= 0) error stop 'residual -'
        err = norm2((rp - rm)/(2*h) - jd)/norm2(jd)
        print '(a,es12.4)', 'Jacobian vs central difference ', err
        if (err > 1.0e-6_dp) error stop 'Jacobian'
    end subroutine
end program
