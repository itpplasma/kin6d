module nn_coupled_case
    ! Solov'ev data for the coupled NN2 solve. With constant F, the exact
    ! state has J = R p' e_phi, Phi = 0 and J_par = p' F B/|B|^2, which is
    ! Ohmic for g_par = g0 (R/R0)^2 and V = 2 pi R0^2 p'/g0. The source
    ! S = -chi_perp p' lap(psi) holds the exact p = p' psi, since
    ! B.grad p = 0.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_solovev_fixture, only: r0, pp, laplacian
    implicit none
    real(dp) :: chi_perp = 1, g0 = 1
contains
    function source(x) result(s)
        real(dp), intent(in) :: x(2)
        real(dp) :: s
        s = -chi_perp*pp*laplacian(x)
    end function
    function g_shape(x) result(w)
        real(dp), intent(in) :: x(2)
        real(dp) :: w
        w = (x(1)/r0)**2
    end function
end module

program test_nonnested_coupled
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_solovev_fixture, only: r0, minor, ftor, pp, mapping, coordinates, exact
    use nn_coupled_case, only: chi_perp, g0, source, g_shape
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes, meridian_qpoint
    use fortsparse, only: csc_t, csc_matvec, csc_transpose, &
        fortsparse_status_t, status_ok
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_build
    use kin6d_nonnested_residual, only: nn_residual
    use kin6d_equilibrium_residual, only: equilibrium_residual_t
    use kin6d_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI, EQ_PSI, &
        EQ_P
    use kin6d_nonnested_coupled, only: nc_problem_t, nc_layout_t, nc_report_t, &
        nc_layout, nc_residual, nc_jacobian, nc_jvp, nc_vjp, nc_solve, nc_fixed
    implicit none
    type(free_b_field_t) :: field
    type(nc_problem_t) :: problem
    type(nc_layout_t) :: lay, lay0
    real(dp), allocatable :: nodes(:, :), x(:), xe(:), x0(:)
    integer :: modes = 0
    real(dp) :: tol_run = 1.0e-11_dp
    real(dp) :: d0(4), dn
    real(dp) :: v_exact
    character(16) :: mode, arg
    integer :: k, lev, nlev
    real(dp) :: ratio, ef(3), ep(3), ev(3), ed(3), eo(3), rate(4)
    call get_command_argument(1, mode)
    call get_command_argument(2, arg)
    read (arg, *) k
    call get_command_argument(3, arg)
    read (arg, *) ratio
    if (k < 2 .or. k > 3 .or. ratio < 1) error stop 'arguments'
    select case (trim(mode))
    case ('solovev')
        nlev = 3
        do lev = 1, nlev
            call run(8*2**(lev - 1), lev == 1, .false., ef(lev), ep(lev), ev(lev), &
                ed(lev), eo(lev))
        end do
        rate = log([ef(1)/ef(3), ep(1)/ep(3), eo(1)/eo(3), ed(1)/ed(3)])/log(4.0_dp)
        print '(a,4f9.3)', 'orders field/p/Phi/D3 ', rate
        if (rate(1) < k - 0.3_dp) error stop 'field order'
        if (rate(2) < k + 0.6_dp) error stop 'pressure order'
        if (rate(4) < k - 0.3_dp) error stop 'D3 order'
        if (ef(3) > merge(5.0e-4_dp, 5.0e-6_dp, k == 2)) error stop 'field error'
        if (maxval(ev) > 1.0e-8_dp) error stop 'loop voltage'
        if (eo(3) >= eo(1) .or. maxval(eo) > 1.0e-5_dp) error stop 'Phi'
        print '(a)', 'PASS: NN2 coupled Solovev nested limit'
    case ('reduction')
        ! An axisymmetric seed solved with modes 0..1 stays axisymmetric and
        ! reproduces the n = 0 solve; the 3D Jacobian is checked at a
        ! non-axisymmetric state.
        tol_run = 1.0e-13_dp
        call run(8, .false., .false., ef(1), ep(1), ev(1), ed(1), eo(1))
        x0 = x; lay0 = lay
        modes = 1
        call run(8, .true., .false., ef(2), ep(2), ev(2), ed(2), eo(2))
        d0 = [blockdiff0(lay%ia, lay0%ia, lay0%na), blockdiff0(lay%ip, lay0%ip, &
            lay0%nv), maxval(abs(x(lay%iphi + 1:lay%iphi + lay0%nv) &
            - x0(lay0%iphi + 1:lay0%iphi + lay0%nv)))*2*acos(-1.0_dp)/abs(v_exact), &
            abs(x(lay%iv) - x0(lay0%iv))/abs(x0(lay0%iv))]
        dn = max(maxval(abs(x(lay%ia + lay0%na + 1:lay%ia + lay%na))) &
            /maxval(abs(x0(:lay0%na))), maxval(abs(x(lay%ip + lay0%nv + 1: &
            lay%ip + lay%nv)))/maxval(abs(x0(lay0%ip + 1:lay0%ip + lay0%nv))))
        print '(a,4es10.2,a,es10.2)', 'N=1 vs N=0 mode-0 A/p/Phi/V ', d0, &
            '  n=1 content ', dn
        print '(a,2es12.4)', 'D3 N=0, N=1 ', ed(1:2)
        if (maxval(d0) > 1.0e-10_dp .or. dn > 1.0e-12_dp) error stop 'reduction'
        if (abs(ed(2) - ed(1)) > 1.0e-8_dp*ed(1)) error stop 'reduction D3'
        print '(a)', 'PASS: NN2 coupled Fourier solve reduces to n = 0'
    case ('matrixfree')
        modes = 1
        if (command_argument_count() >= 4) then
            call get_command_argument(4, arg)
            read (arg, *) modes
        end if
        if (command_argument_count() >= 5) then
            call get_command_argument(5, arg)
            read (arg, *) problem%nfp
        end if
        tol_run = 1.0e-13_dp
        call run(8, .true., .false., ef(1), ep(1), ev(1), ed(1), eo(1))
        print '(a)', 'PASS: dense / FFT-JVP Newton iterates and solution'
    case ('offset')
        call run(16, .false., .true., ef(1), ep(1), ev(1), ed(1), eo(1))
        print '(a,es12.4)', 'offset return difference ', eo(1)
        ! Both solves end at the rounding floor of the backward error; the
        ! forward difference is then bounded by cond * eps, and the
        ! transport condition number grows like chi_par/chi_perp.
        if (eo(1) > max(1.0e-11_dp, 10*ratio*epsilon(1.0_dp))) &
            error stop 'offset return'
        print '(a)', 'PASS: NN2 coupled offset-start return'
    case default
        error stop 'mode'
    end select
contains
    subroutine run(nr, check_jacobian, offset, efield, epres, evolt, ed3, aux)
        integer, intent(in) :: nr
        logical, intent(in) :: check_jacobian, offset
        real(dp), intent(out) :: efield, epres, evolt, ed3
        real(dp), intent(out) :: aux ! max |Phi| in V/(2 pi), or the offset difference
        type(flux_mesh_t) :: fm
        type(nc_report_t) :: report
        type(equilibrium_residual_t) :: stress
        real(dp), allocatable :: rings(:), points(:, :), xd(:), td(:, :), tf(:, :)
        type(nc_report_t) :: dense_report
        real(dp) :: difference
        real(dp) :: values(EQ_VALUES, 128), b(2), psi, area, phi_max
        real(dp) :: sum_b, sum_p, theta, xq(2), wq, basis(3, 64), diffs(4)
        integer :: i, j, t, q, st, statuses(128)
        allocate (rings(nr), points(3, 128))
        do i = 1, nr
            rings(i) = real(i, dp)/nr
        end do
        call flux_mesh_build(rings, 8, 1024, k, mapping, fm, st)
        if (st /= 0) error stop 'mesh'
        call nn_build(fm, k, modes, [r0, 0.0_dp], minor, ftor, 1.0_dp, field, st)
        if (st /= 0) error stop 'build'
        field%coordinates => coordinates
        call meridian_lagrange_nodes(field%fs, nodes)
        lay = nc_layout(field)
        area = 0
        do t = 1, field%fs%mm%topology%n_triangles
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 0, t, q, xq, wq, &
                    basis(:, :field%fs%nl0), st)
                if (st /= 0) error stop 'quadrature'
                area = area + wq
            end do
        end do
        chi_perp = 1; g0 = 1
        problem%chi_perp = chi_perp; problem%chi_par = ratio*chi_perp
        problem%g_par = g0; problem%g_perp = g0/ratio
        problem%source => source; problem%g_par_shape => g_shape
        problem%current = pp*area
        v_exact = 2*acos(-1.0_dp)*r0*r0*pp/g0
        call seed(0.1_dp, 0.0_dp, x)
        if (check_jacobian) call jacobian_check(x)
        if (trim(mode) == 'matrixfree') then
            xd = x
            call nc_solve(field, problem, xd, dense_report, st, &
                tol=tol_run, iterates=td)
            if (st /= 0) error stop 'dense comparison Newton'
            problem%matrix_free = .true.
            call nc_solve(field, problem, x, report, st, tol=tol_run, iterates=tf)
            print '(a,i0,a,40i6)', 'matrix-free status ', st, ' Krylov ', &
                report%krylov_iterations(:report%newton_steps)
            if (st /= 0) error stop 'matrix-free comparison Newton'
            if (report%newton_steps /= dense_report%newton_steps) &
                error stop 'Newton iterate count'
            difference = maxval(abs(td(:, :report%newton_steps) &
                - tf(:, :report%newton_steps)))/maxval(abs(td))
            print '(a,es12.4)', 'relative Newton iterate difference ', difference
            if (difference > 1.0e-10_dp) error stop 'Newton iterate equivalence'
            if (maxval(abs(x - xd))/maxval(abs(xd)) > 1.0e-10_dp) &
                error stop 'solution equivalence'
        end if
        ! The offset comparison converges both solves to the rounding floor.
        if (trim(mode) /= 'matrixfree') call nc_solve(field, problem, x, report, st, &
            tol=merge(1.0e-13_dp, tol_run, offset))
        print '(a,i0,a,i0,a,15es9.1)', 'P', k, ' nr=', nr, ' Newton ', &
            report%history(:report%newton_steps)
        if (st /= 0) error stop 'coupled Newton'
        if (report%newton_steps > 12) error stop 'Newton steps'
        if (report%divergence > 1.0e-12_dp) error stop 'closed 2-form'
        do i = 1, 4
            do j = 1, 32
                theta = 2*acos(-1.0_dp)*(j - 0.37_dp)/32
                call mapping(0.15_dp + 0.18_dp*i, theta, points(:2, j + 32*(i - 1)))
                points(3, j + 32*(i - 1)) = 0.43_dp
            end do
        end do
        call field%sample(points, values, statuses)
        if (any(statuses /= 0)) error stop 'sampling'
        efield = 0; epres = 0; sum_b = 0; sum_p = 0
        do i = 1, 128
            call exact(points(:2, i), psi, b)
            efield = efield + sum(([values(EQ_BR, i), values(EQ_BZ, i)] - b)**2) &
                + (values(EQ_BPHI, i) - ftor/points(1, i))**2
            epres = epres + (values(EQ_P, i) - pp*psi)**2
            sum_b = sum_b + sum(b*b); sum_p = sum_p + (pp*psi)**2
        end do
        efield = sqrt(efield/sum_b); epres = sqrt(epres/sum_p)
        evolt = abs(report%voltage - v_exact)/abs(v_exact)
        phi_max = maxval(abs(x(lay%iphi + 1:lay%iphi + lay%nv)))*2*acos(-1.0_dp) &
            /abs(v_exact)
        call nn_residual(field, 0, stress, st)
        if (st /= 0) error stop 'D3'
        ed3 = stress%normalized
        print '(a,i0,a,i0,a,4es11.3,a,2es10.2)', 'P', k, ' nr=', nr, &
            ' field/p/V/D3 ', efield, epres, evolt, ed3, ' Phi/lambda ', phi_max, &
            report%multiplier
        aux = phi_max
        if (report%multiplier > 1.0e-12_dp*maxval(abs(x(:lay%na)))) &
            error stop 'gauge multiplier'
        if (offset) then
            ! The same state from the exact nodal interpolant and V_exact.
            call seed(0.0_dp, v_exact, xe)
            call nc_solve(field, problem, xe, report, st, tol=1.0e-13_dp)
            if (st /= 0) error stop 'second Newton'
            print '(a,15es9.1)', 'interpolant seed Newton ', &
                report%history(:report%newton_steps)
            ! Phi is compared in units of the loop potential V/(2 pi).
            diffs = [blockdiff(lay%ia, lay%na, maxval(abs(xe(:lay%na)))), &
                blockdiff(lay%ip, lay%nv, maxval(abs(xe(lay%ip + 1:lay%iphi)))), &
                blockdiff(lay%iphi, lay%nv, abs(v_exact)/(2*acos(-1.0_dp))), &
                blockdiff(lay%iv - 1, 1, abs(v_exact))]
            print '(a,4es10.2)', 'offset difference A/p/Phi/V ', diffs
            aux = maxval(diffs)
        end if
    end subroutine

        subroutine seed(shift, voltage, s)
            ! Axis shifted by `shift` (the shift vanishes on the boundary).
            real(dp), intent(in) :: shift, voltage
            real(dp), allocatable, intent(out) :: s(:)
            real(dp) :: u, th, bb(2), ps
            logical, allocatable :: fixed(:)
            integer :: n
            allocate (s(lay%n))
            s = 0
            do n = 1, lay%n0
                call coordinates([r0, 0.0_dp], minor, nodes(:, n), u, th)
                call exact([nodes(1, n) - shift*(1 - u*u), nodes(2, n)], ps, bb)
                s(lay%ia + field%fs%nn + n) = ps
                s(lay%ip + n) = pp*ps
            end do
            s(lay%iv) = voltage
            fixed = nc_fixed(field, problem)
            where (fixed) s = 0
        end subroutine

        real(dp) function blockdiff0(offset1, offset0, n) result(d)
            ! Mode-0 block of the N = 1 state against the N = 0 state.
            integer, intent(in) :: offset1, offset0, n
            d = maxval(abs(x(offset1 + 1:offset1 + n) - x0(offset0 + 1:offset0 + n))) &
                /max(tiny(1.0_dp), maxval(abs(x0(offset0 + 1:offset0 + n))))
        end function

        real(dp) function blockdiff(offset0, n, scale) result(d)
            integer, intent(in) :: offset0, n
            real(dp), intent(in) :: scale
            d = maxval(abs(x(offset0 + 1:offset0 + n) &
                - xe(offset0 + 1:offset0 + n)))/scale
        end function

        subroutine jacobian_check(x0)
            ! Central-difference directional derivative of the residual
            ! against the assembled Jacobian (free directions only).
            real(dp), intent(in) :: x0(:)
            type(csc_t) :: jac, jac_t
            type(fortsparse_status_t) :: sparse_status
            real(dp), allocatable :: d(:), rp(:), rm(:), jd(:), xs(:), adj(:), weight(:)
            logical, allocatable :: fixed(:)
            real(dp) :: h, err, pairing(2)
            integer :: n, st
            fixed = nc_fixed(field, problem)
            allocate (d(lay%n))
            xs = abs(x0)
            do n = 1, lay%n
                d(n) = sin(1.3_dp*n)
            end do
            ! Scale each block by its seed magnitude (V and Phi by V_exact).
            d(lay%ia + 1:lay%ia + lay%na) = d(lay%ia + 1:lay%ia + lay%na) &
                *maxval(xs(lay%ia + 1:lay%ia + lay%na))
            d(lay%il + 1:lay%il + lay%nv) = d(lay%il + 1:lay%il + lay%nv)*1.0e-3_dp
            d(lay%ip + 1:lay%ip + lay%nv) = d(lay%ip + 1:lay%ip + lay%nv) &
                *maxval(xs(lay%ip + 1:lay%ip + lay%nv))
            d(lay%iphi + 1:lay%iv) = d(lay%iphi + 1:lay%iv)*abs(v_exact)
            where (fixed) d = 0
            call nc_jacobian(field, problem, x0 + d, jac, st)
            if (st /= 0) error stop 'Jacobian'
            jd = csc_matvec(jac, d)
            where (fixed) jd = 0
            call nc_jvp(field, problem, x0 + d, d, rp, st)
            if (st /= 0) error stop 'FFT JVP'
            err = norm2(rp - jd)/norm2(jd)
            print '(a,es12.4)', 'FFT JVP vs assembled Jacobian ', err
            if (err > 1.0e-10_dp) error stop 'FFT JVP equivalence'
            weight = [(cos(0.73_dp*n), n=1, lay%n)]
            call nc_vjp(field, problem, x0 + d, weight, adj, st)
            if (st /= 0) error stop 'FFT VJP'
            pairing = [dot_product(weight, rp), dot_product(adj, d)]
            err = abs(pairing(1) - pairing(2))/maxval(abs(pairing))
            print '(a,es12.4)', 'FFT JVP/VJP pairing ', err
            if (err > 1.0e-10_dp) error stop 'FFT adjoint pairing'
            call csc_transpose(jac, jac_t, sparse_status)
            if (.not. status_ok(sparse_status)) error stop 'Jacobian transpose'
            rp = csc_matvec(jac_t, merge(0.0_dp, weight, fixed))
            err = norm2(adj - rp)/norm2(rp)
            print '(a,es12.4)', 'FFT VJP vs assembled transpose ', err
            if (err > 1.0e-10_dp) error stop 'FFT VJP equivalence'
            h = 1.0e-4_dp
            call nc_residual(field, problem, x0 + d*(1 + h), rp, st)
            if (st /= 0) error stop 'residual +'
            call nc_residual(field, problem, x0 + d*(1 - h), rm, st)
            if (st /= 0) error stop 'residual -'
            err = norm2((rp - rm)/(2*h) - jd)/norm2(jd)
            print '(a,es12.4)', 'Jacobian vs central difference ', err
            if (err > 1.0e-7_dp) error stop 'Jacobian'
        end subroutine
end program
