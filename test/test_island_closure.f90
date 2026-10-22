module island_test_data
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_solovev_fixture, only: r0, pp, laplacian
    implicit none
contains
    function source(x) result(s)
        real(dp), intent(in) :: x(2)
        real(dp) :: s
        s = -pp*laplacian(x)
    end function
    function constant_source(x) result(s)
        real(dp), intent(in) :: x(2)
        real(dp) :: s
        s = 3 + 0*x(1)
    end function
    function negative_shape(x) result(w)
        real(dp), intent(in) :: x(2)
        real(dp) :: w
        w = -1 + 0*x(1)
    end function
    function g_shape(x) result(w)
        real(dp), intent(in) :: x(2)
        real(dp) :: w
        w = (x(1)/r0)**2
    end function
end module

program test_island_closure
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use nn_solovev_fixture, only: r0, minor, ftor, pp, mapping, coordinates, exact
    use island_test_data, only: source, constant_source, g_shape, negative_shape
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes, &
        meridian_qpoint, meridian_qgrad
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_build, nn_riesz
    use kin6d_nonnested_coupled, only: nc_problem_t, nc_report_t, nc_layout_t, &
        nc_layout, nc_solve, nc_fixed, nc_validate
    use kin6d_island_case, only: island_case_read
    use kin6d_island_closure, only: ISLAND_INVALID, ISLAND_UNSUPPORTED, &
        ISLAND_DELEGATE_B, ISLAND_DELEGATE_C, ISLAND_NULL_FIELD, ISLAND_BAD_DATA
    use kin6d_island_report, only: island_report_t, island_residual_report
    implicit none
    type(nc_problem_t) :: problem, other
    type(nc_report_t) :: report
    type(free_b_field_t) :: field
    type(flux_mesh_t) :: mesh
    type(nc_layout_t) :: lay
    type(island_report_t) :: diagnostics
    real(dp), allocatable :: nodes(:, :), x(:), xi(:), xv(:), rings(:)
    logical, allocatable :: fixed(:)
    character(160) :: message
    character(24) :: mode
    real(dp) :: ps, b(2), u, th, area, v_exact, error, current_error
    integer :: st, i, lev, nr
    call get_command_argument(1, mode)
    problem%source => source
    select case (trim(mode))
    case ('declaration')
        call declarations()
    case ('drive')
        ! Three meshes, offset starts. Match the *discrete* voltage of the
        ! current-driven state, then independently solve the voltage-driven case.
        do lev = 1, 3
            nr = 2**lev
            call build(nr)
            problem%g_par_shape => g_shape
            problem%chi_par = 10; problem%g_perp = 0.1_dp
            problem%current = pp*area
            v_exact = 2*acos(-1.0_dp)*r0*r0*pp
            x = 0
            do i = 1, lay%n0
                call coordinates([r0, 0.0_dp], minor, nodes(:, i), u, th)
                call exact([nodes(1, i) - 0.03_dp*(1 - u*u), nodes(2, i)], ps, b)
                x(field%fs%nn + i) = ps
                x(lay%ip + i) = pp*ps
            end do
            fixed = nc_fixed(field)
            where (fixed) x = 0
            xi = x
            call nc_solve(field, problem, x, report, st, tol=1.0e-12_dp)
            call require(st == 0, 'current-driven solve')
            call require(maxval(report%residual) <= 1.0e-10_dp, 'backward error')
            call require(abs(x(lay%iv)/v_exact - 1) < 1.0e-6_dp, 'voltage oracle')
            other = problem
            other%closure%drive = 'V'; other%voltage = x(lay%iv)
            xv = xi
            call nc_solve(field, other, xv, report, st, tol=1.0e-12_dp)
            call require(st == 0, 'voltage-driven solve')
            error = maxval(abs(x - xv))/maxval(abs(x))
            call require(error <= 1.0e-7_dp, 'matched V/I states')
           call island_residual_report(field, other, xv, minor, 1.0_dp, diagnostics, st)
            call require(st == 0, 'solved residual report')
            current_error = abs(diagnostics%current/problem%current - 1)
            call require(current_error <= 1.0e-6_dp, 'independent current integral')
            call require(diagnostics%charge_balance <= 1.0e-9_dp, 'charge balance')
            call require(diagnostics%heat_balance <= 1.0e-9_dp, 'heat balance')
            print '(a,i0,4es12.3)', 'V/I nr, state/current/charge/heat ', nr, &
                error, current_error, diagnostics%charge_balance, diagnostics%heat_balance
        end do
    case ('report')
        call build(2)
        call report_oracle()
    case default
        error stop 'mode'
    end select
    print '(a)', 'PASS: I1 '//trim(mode)
contains
    subroutine require(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print '(a,i0,2a)', 'status ', st, ': ', trim(message)
            print '(a)', label
            error stop 'I1 gate'
        end if
    end subroutine

    subroutine declarations()
        other = problem
        call nc_validate(other, st, message)
        call require(st == 0, 'legacy defaults')
        call read_case('&island_case /', 0)
     call read_case("&island_case pressure_closure='B_frozen_fit' /", ISLAND_DELEGATE_B)
        call read_case("&island_case pressure_closure='C', current_closure='C' /", &
            ISLAND_DELEGATE_C)
        call read_case("&island_case pressure_closure='C' /", ISLAND_INVALID)
        call read_case("&island_case current_closure='q_profile' /", ISLAND_INVALID)
        call read_case('&island_case I_phi=1, V=2 /', ISLAND_INVALID)
        call read_case("&island_case drive='both' /", ISLAND_INVALID)
        call read_case("&island_case drive='I_phi', V=2 /", ISLAND_INVALID)
        call read_case("&island_case heat_boundary='neumann' /", ISLAND_UNSUPPORTED)
        call read_case("&island_case target_class='ideal' /", ISLAND_UNSUPPORTED)
        call read_case('&island_case chi_perp=0 /', ISLAND_INVALID)
        call read_case('&island_case unknown=1 /', ISLAND_INVALID)
        call read_case("&island_case pressure_closure(2)='A' /", ISLAND_INVALID)
        call read_case('&island_case V=2 /', 0)
        call require(other%closure%drive == 'V', 'voltage selection')
        other = problem
        other%chi_par = ieee_value(1.0_dp, ieee_quiet_nan)
        call nc_validate(other, st, message)
        call require(st == ISLAND_INVALID, 'nonfinite coefficient')
        other = problem
        nullify (other%source)
        call nc_validate(other, st, message)
        call require(st == ISLAND_INVALID, 'undeclared source')
    end subroutine

    subroutine read_case(text, expected)
        character(*), intent(in) :: text
        integer, intent(in) :: expected
        integer :: unit
        other = problem
        open (newunit=unit, status='scratch', action='readwrite')
        write (unit, '(a)') text
        rewind (unit)
        call island_case_read(unit, other, st, message)
        close (unit)
        call require(st == expected, text)
        if (st /= 0) then
            call require(len_trim(message) > 0, 'clear rejection status')
            call require(other%closure%pressure == 'A', 'transactional rejection')
        end if
    end subroutine

    subroutine build(nr)
        integer, intent(in) :: nr
        real(dp) :: v(3, 64), xq(2), wq
        integer :: t, q
        if (allocated(rings)) deallocate (rings, x)
        rings = [(real(i, dp)/nr, i=1, nr)]
        call flux_mesh_build(rings, 8, 1024, 3, mapping, mesh, st)
        call require(st == 0, 'mesh')
        call nn_build(mesh, 3, 0, [r0, 0.0_dp], minor, ftor, 1.0_dp, field, st)
        call require(st == 0, 'field')
        field%coordinates => coordinates
        call meridian_lagrange_nodes(field%fs, nodes)
        lay = nc_layout(field)
        allocate (x(lay%n))
        area = 0
        do t = 1, field%fs%mm%topology%n_triangles
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 0, t, q, xq, wq, v(:, :field%fs%nl0), st)
                area = area + wq
            end do
        end do
    end subroutine

    subroutine report_oracle()
        ! Independent handwritten oracle: B=F/R e_phi, p=R+2Z,
        ! Phi=0.3R-0.2Z, constant source. Both physical gradients are exact
        ! under the isoparametric mapping. This is force-incompatible data.
        real(dp) :: v1(3, field%fs%sp(1)%nloc), val(field%fs%nl0)
        real(dp) :: grad(2, field%fs%nl0), xq(2), wq, j(3), flux(2), two_pi
        real(dp) :: expected(lay%n), volume, ohmic2, heating, current
        complex(dp) :: rhs(lay%n0), riesz(lay%n0)
        real(dp) :: dual2, means, bn
        integer :: t, q, d1(field%fs%sp(1)%nloc), d0(field%fs%nl0)
        problem%source => constant_source
        problem%current = 1.0e8_dp ! current diagnostics must avoid target cancellation
        problem%g_perp = 0.4_dp; problem%chi_perp = 0.7_dp
        problem%chi_par = 4; problem%g_par = 2
        x = 0; x(lay%iv) = 1.3_dp
        x(lay%ip + 1:lay%iphi) = nodes(1, :) + 2*nodes(2, :)
        x(lay%iphi + 1:lay%iv - 1) = 0.3_dp*nodes(1, :) - 0.2_dp*nodes(2, :)
        call island_residual_report(field, problem, x, 2.0_dp, 3.0_dp, diagnostics, st)
        call require(st == 0, 'report at unconverged state')
        expected = 0; volume = 0; ohmic2 = 0; heating = 0; current = 0
        rhs = 0; means = 0; two_pi = 2*acos(-1.0_dp)
        do t = 1, field%fs%mm%topology%n_triangles
            d0 = field%fs%sp(0)%dof(:, t); d1 = field%fs%sp(1)%dof(:, t)
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 1, t, q, xq, wq, v1, st)
                call meridian_qgrad(field%fs, t, q, xq, wq, val, grad, st)
                bn = ftor/xq(1)
                j = [2/bn - 0.3_dp*problem%g_perp, &
                    problem%g_par*x(lay%iv)/(two_pi*xq(1)), &
                    -1/bn + 0.2_dp*problem%g_perp]
                flux = problem%chi_perp*[1.0_dp, 2.0_dp]
                expected(d1) = expected(d1) - wq*matmul(j, v1)*field%fs%sp(1)%sgn(:, t)
                expected(lay%ip + d0) = expected(lay%ip + d0) &
                    + wq*(matmul(flux, grad) - 3*val)
                expected(lay%iphi + d0) = expected(lay%iphi + d0) &
                    + wq*(j(1)*grad(1, :) + j(3)*grad(2, :))
                current = current + wq*j(2)/xq(1)
                volume = volume + wq; heating = heating + 3*wq
                ohmic2 = ohmic2 + wq*problem%g_perp**2*bn**2*(0.3_dp**2 + 0.2_dp**2)
                rhs(d0) = rhs(d0) + wq*val
                means = means + wq*(xq(1) + 2*xq(2))
            end do
        end do
        call require(abs(diagnostics%outward_current - two_pi*sum( &
            expected(lay%iphi + 1:lay%iv - 1), mask=field%fs%sp(0)%bnd)) &
            < 1.0e-10_dp, 'boundary current reaction oracle')
        error = abs(diagnostics%outward_heat + two_pi*sum( &
            expected(lay%ip + 1:lay%iphi), mask=field%fs%sp(0)%bnd))
        call require(error < 1.0e-10_dp, 'boundary heat reaction oracle')
        fixed = nc_fixed(field, problem)
        where (fixed) expected = 0
        call require(abs(diagnostics%ampere - norm2(expected(:lay%il))) &
            < 1.0e-11_dp, 'Ampere oracle')
        call require(diagnostics%gauge == 0, 'gauge oracle')
        call require(abs(diagnostics%pressure - norm2(expected(lay%ip + 1:lay%iphi))) &
            < 1.0e-11_dp, 'pressure oracle')
        call require(abs(diagnostics%ohm_continuity &
            - norm2(expected(lay%iphi + 1:lay%iv - 1))) &
            < 1.0e-11_dp, 'Ohm/continuity oracle')
        call require(abs(diagnostics%current - current) < 1.0e-11_dp, 'current oracle')
        call require(abs(diagnostics%drive - (current - problem%current)) &
            < 1.0e-11_dp, 'drive oracle')
  call require(abs(diagnostics%heating - two_pi*heating) < 1.0e-10_dp, 'heating oracle')
     call require(diagnostics%epsilon_parallel < 1.0e-12_dp, 'parallel pressure oracle')
     call require(abs(diagnostics%perpendicular_ohmic_force_rms - sqrt(ohmic2/volume)) &
         < 1.0e-11_dp, 'Ohmic force oracle')
        where (field%fs%sp(0)%bnd) rhs = 0
        call nn_riesz(field, 0, rhs, riesz, st)
        dual2 = 5*real(dot_product(rhs, riesz), dp)
        call require(abs(diagnostics%force%dual - sqrt(two_pi*dual2)) &
            < 1.0e-10_dp, 'D3 independent load')
        call require(abs(diagnostics%force%p_mean - means/volume) &
            < 1.0e-11_dp, 'D3 mean pressure')
 call require(diagnostics%algebraic_divergence < 1.0e-12_dp, 'exact-complex divergence')
        problem%g_par_shape => negative_shape
        call nc_solve(field, problem, x, report, st)
        call require(st == ISLAND_BAD_DATA, 'negative conductivity shape')
  call require(index(report%message, 'conductivity') > 0, 'pointwise rejection message')
        nullify (problem%g_par_shape)
        field%f_vac = 0; x = 0
        call nc_solve(field, problem, x, report, st)
        call require(st == ISLAND_NULL_FIELD, 'magnetic null rejection')
        call require(index(report%message, 'Null') > 0, 'magnetic null message')
    end subroutine
end program
