module fitted_island_data
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_state, only: equilibrium_field_t, EQ_BPHI
    use kin6d_nonnested_profile, only: nn_profile_t
    implicit none
    type, extends(equilibrium_field_t) :: seed_t
        real(dp) :: toroidal = 1
    contains
        procedure :: sample => sample_seed
    end type
    type, extends(nn_profile_t) :: target_t
        real(dp) :: amplitude = 0, drive = 0
        logical :: invalid = .false.
    contains
        procedure :: evaluate => evaluate_target
    end type
contains
    subroutine sample_seed(field, points, values, statuses)
        class(seed_t), intent(in) :: field
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        values = 0; statuses = 0
        values(EQ_BPHI, :) = field%toroidal/points(1, :)
    end subroutine
    subroutine evaluate_target(self, x, p, gradp, lambda, gradphi, status)
        class(target_t), intent(in) :: self
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: p(:), gradp(:, :), lambda(:), gradphi(:, :)
        integer, intent(out) :: status
        p = self%amplitude*(x(1, :) + 2*x(2, :))
        gradp(1, :) = self%amplitude; gradp(2, :) = 0
        gradp(3, :) = 2*self%amplitude
        lambda = self%drive; gradphi = 0
        status = 0
        if (self%invalid) status = 1
    end subroutine
    function zero_source(x) result(value)
        real(dp), intent(in) :: x(2)
        real(dp) :: value
        value = 0*x(1)
    end function
end module

program test_island_fitted
    ! Independent public oracles at prescribed, generally force-incompatible
    ! states. A weak fitted-pressure reaction is not physical outward heat.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan
    use fitted_island_data
    use nn_solovev_fixture, only: r0, minor, ftor, mapping, coordinates
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes, meridian_qpoint
    use fortsparse, only: csc_t, csc_matvec
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_build
    use kin6d_nonnested_profile, only: nn_fit_source
    use kin6d_nonnested_coupled, only: nc_problem_t, nc_layout_t, nc_layout, &
        nc_validate, nc_source_data, nc_fixed, nc_reduce, nc_residual, &
        nc_jacobian, nc_jvp, nc_vjp
    use kin6d_nonnested_residual, only: nn_residual, NN_RESIDUAL_ZERO_SCALE
    use kin6d_equilibrium_residual, only: equilibrium_residual_t
    use kin6d_island_case, only: island_case_read
    use kin6d_island_closure, only: ISLAND_INVALID, ISLAND_DELEGATE_B, ISLAND_BAD_DATA
    use kin6d_island_report, only: island_report_t, island_residual_report
    implicit none
    type(nc_problem_t) :: problem, other
    type(flux_mesh_t) :: mesh
    type(free_b_field_t) :: field
    type(nc_layout_t) :: lay
    type(island_report_t) :: report
    real(dp), allocatable :: nodes(:, :), x(:)
    character(24) :: mode
    character(160) :: message = ''
    integer :: st
    call get_command_argument(1, mode)
    call flux_mesh_build([0.5_dp, 1.0_dp], 8, 1024, 3, mapping, mesh, st)
    call require(st == 0, 'mesh')
    call nn_build(mesh, 3, 0, [r0, 0.0_dp], minor, ftor, 1.0_dp, field, st)
    call require(st == 0, 'field')
    field%coordinates => coordinates
    call meridian_lagrange_nodes(field%fs, nodes)
    lay = nc_layout(field); allocate (x(lay%n)); x = 0
    select case (trim(mode))
    case ('readiness')
        call readiness()
    case ('heat')
        call heat()
    case ('current')
        call current()
    case default
        error stop 'mode'
    end select
    print '(a)', 'PASS: I5A '//trim(mode)
contains
    subroutine require(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (ok) return
        print '(a,i0,2a)', 'status ', st, ': ', trim(message)
        print '(a)', label
        error stop 'I5A gate'
    end subroutine

    subroutine fit(p, amplitude, drive, invalid, null_seed)
        type(nc_problem_t), intent(out) :: p
        real(dp), intent(in) :: amplitude, drive
        logical, intent(in), optional :: invalid, null_seed
        type(seed_t) :: seed
        type(target_t) :: target
        seed%toroidal = ftor
        if (present(null_seed)) then
            if (null_seed) seed%toroidal = 0
        end if
        target%amplitude = amplitude; target%drive = drive
        if (present(invalid)) target%invalid = invalid
        p%closure%pressure = 'B_frozen_fit'; allocate (p%fitted)
        call nn_fit_source(seed, target, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
            0.0_dp, p%fitted, st)
        call require(st == 0, 'fit construction')
        ! The fit must retain owned copies, independent of caller mutation.
        target%drive = 999; target%amplitude = 999; seed%toroidal = 0
    end subroutine

    subroutine readiness()
        type(equilibrium_residual_t) :: force
        real(dp) :: scalar, flux(3), drive, before(4)
        integer :: unit
        problem%source => zero_source
        call nc_validate(problem, st, message)
        call require(st == 0, 'scalar A retained')
        problem%closure%pressure = 'B_frozen_fit'
        call nc_validate(problem, st, message)
        call require(st == ISLAND_DELEGATE_B, 'missing fit delegated')
        allocate (problem%fitted)
        call nc_validate(problem, st, message)
        call require(st == ISLAND_DELEGATE_B, 'uninitialised fit delegated')
        call fit(problem, 1.0_dp, -0.37_dp)
        call nc_validate(problem, st, message)
        call require(st == 0, 'ready B admitted')
        call nc_source_data(problem, [r0, 0.0_dp], 0.0_dp, scalar, flux, drive, st)
        call require(st == 0, 'owned fit samples')
        call require(abs(scalar) < tiny(1.0_dp), 'weak fit has no scalar source')
        call require(maxval(abs(flux - [1.0_dp, 0.0_dp, 2.0_dp])) < 1.0e-12_dp &
            .and. abs(drive + 0.37_dp) < 1.0e-12_dp, 'copied fit independent of caller')
        before = [flux, drive]
        other = problem; other%source => zero_source
        call nc_validate(other, st, message)
        call require(st == ISLAND_INVALID, 'conflicting scalar hook rejected')
        other = problem; other%g_par_shape => zero_source
        call nc_validate(other, st, message)
        call require(st == ISLAND_INVALID, 'conflicting shape hook rejected')
        other = problem; other%closure%pressure = 'A'
        call nc_validate(other, st, message)
        call require(st == ISLAND_INVALID, 'fit requires explicit B closure')
        other = problem; other%closure%pressure = 'B_projection'
        call nc_validate(other, st, message)
        call require(st == ISLAND_DELEGATE_B, 'projection remains delegated')
        other = problem; other%chi_par = ieee_value(0.0_dp, ieee_quiet_nan)
        call nc_validate(other, st, message)
        call require(st == ISLAND_INVALID, 'nonfinite live coefficients rejected')
        other = problem
        open (newunit=unit, status='scratch', action='readwrite')
        write (unit, '(a)') "&island_case pressure_closure='A' /"
        rewind (unit)
        call island_case_read(unit, other, st, message); close (unit)
        call require(st == ISLAND_INVALID, 'incompatible case rejected')
        call require(other%closure%pressure == 'B_frozen_fit', 'transactional declaration')
        call nc_source_data(other, [r0, 0.0_dp], 0.0_dp, scalar, flux, drive, st)
        call require(st == 0 .and. maxval(abs([flux, drive] - before)) < 1.0e-12_dp, &
            'transactional fitted data')
        call fit(other, 0.0_dp, 0.0_dp, invalid=.true.)
        call island_residual_report(field, other, x, minor, 1.0_dp, report, st)
        call require(st == ISLAND_BAD_DATA, 'point-data failure is not readiness')
        call fit(other, 0.0_dp, 0.0_dp, null_seed=.true.)
        call nc_source_data(other, [r0, 0.0_dp], 0.0_dp, scalar, flux, drive, st)
        call require(st == ISLAND_BAD_DATA, 'null seed rejected without field floor')
        call fit(other, 0.0_dp, 0.0_dp)
        call island_residual_report(field, other, x, minor, 1.0_dp, report, st, test_modes=-1)
        call require(st /= 0, 'other force-report failures still propagate')
        field%f_vac = ieee_value(0.0_dp, ieee_quiet_nan)
        call nn_residual(field, 0, force, st)
        call require(st /= 0 .and. st /= NN_RESIDUAL_ZERO_SCALE, &
            'nonfinite harmonic is an error, not zero-scale availability')
        call island_residual_report(field, other, x, minor, 1.0_dp, report, st)
        call require(st /= 0 .and. st /= NN_RESIDUAL_ZERO_SCALE, &
            'nonfinite point data still reject full report')
        field%f_vac = ftor
    end subroutine

    subroutine heat()
        real(dp) :: basis(3, field%fs%nl0), xy(2), wq, physical_heating
        integer :: t, q
        call fit(problem, 1.0_dp, 0.0_dp)
        x = 0; x(lay%ip + 1:lay%iphi) = nodes(1, :) + 2*nodes(2, :)
        call island_residual_report(field, problem, x, minor, 1.0_dp, report, st)
        call require(st == 0, 'weak fitted report')
        call require(report%force%sigma > 0 .and. ieee_is_finite(report%force%normalized), &
            'nonzero D3 scale retains finite normalization')
        call require(index(report%force_status, 'finite-test-space') == 1, &
            'nonzero-scale force availability')
        call require(abs(report%weak_heat_reaction) < 1.0e-9_dp &
            .and. report%weak_heat_balance < 1.0e-9_dp, 'matching live/seed weak flux')
        call require(ieee_is_nan(report%heating) .and. ieee_is_nan(report%outward_heat) &
            .and. ieee_is_nan(report%heat_balance), 'physical heat explicitly unavailable')
        call require(index(report%heat_status, 'unavailable') == 1, 'physical heat status')
        ! Here q_seed=(1,0,2), hence cylindrical S=-div(q_seed)=-1/R.
        ! Its nonzero volume integral is not the zero relative weak reaction.
        physical_heating = 0
        do t = 1, field%fs%mm%topology%n_triangles
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 0, t, q, xy, wq, basis, st)
                call require(st == 0, 'heat oracle quadrature')
                physical_heating = physical_heating - 2*acos(-1.0_dp)*wq/xy(1)
            end do
        end do
        call require(abs(physical_heating) > 1.0e-3_dp, 'physical heating differs from reaction')
        print '(a,3es14.5)', 'weak reaction/balance and analytic heating ', &
            report%weak_heat_reaction, report%weak_heat_balance, physical_heating
    end subroutine

    subroutine current()
        type(csc_t) :: jac
        type(equilibrium_residual_t) :: force
        real(dp), allocatable :: direction(:), weight(:), action(:), adjoint(:)
        real(dp), allocatable :: reference(:), plus(:), minus(:)
        logical, allocatable :: fixed(:)
        real(dp) :: basis(3, field%fs%nl0), xy(2), wq, expected, d, error, h
        real(dp) :: physical_b2, volume
        integer :: sign, t, q, i
        d = 0.37_dp; expected = 0; physical_b2 = 0; volume = 0
        do t = 1, field%fs%mm%topology%n_triangles
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 0, t, q, xy, wq, basis, st)
                call require(st == 0, 'current oracle quadrature')
                expected = expected + wq*d*ftor/xy(1)**2
                physical_b2 = physical_b2 + wq*(ftor/xy(1))**2
                volume = volume + wq
            end do
        end do
        call require(physical_b2/volume > 0, 'nonzero physical harmonic magnetic field')
        do sign = -1, 1, 2
            call fit(problem, 0.0_dp, sign*d); x = 0
            call island_residual_report(field, problem, x, minor, 1.0_dp, report, st)
            call require(st == 0, 'signed current report')
            call require(report%force%sigma == 0 .and. report%force%b2_mean == 0 &
                .and. report%force%p_mean == 0 .and. report%force%dual == 0, &
                'retained D3 excludes harmonic background and has zero absolute load')
            call require(ieee_is_nan(report%force%normalized) .and. &
                index(report%force_status, 'normalization unavailable') == 1, &
                'zero-scale force normalization explicitly unavailable')
            call nn_residual(field, 0, force, st)
            call require(st == NN_RESIDUAL_ZERO_SCALE .and. ieee_is_nan(force%normalized), &
                'standalone D3 retains explicit unavailable status')
            call require(abs(report%current - sign*expected) < 1.0e-11_dp, &
                'signed current has meridian measure without 2pi')
            call require(report%charge_balance < 1.0e-11_dp, 'toroidal drive has no outward charge')
        end do
        ! Check nonlinear fit-frozen actions at nonstationary states, including
        ! boundary columns, against assembled and residual finite differences.
        direction = [(0.001_dp*sin(1.3_dp*i), i=1, lay%n)]
        weight = [(cos(0.73_dp*i), i=1, lay%n)]
        x = 0.1_dp*direction; fixed = nc_fixed(field, problem)
        call nc_jacobian(field, problem, x, jac, st)
        call require(st == 0, 'fitted assembled Jacobian')
        reference = csc_matvec(jac, direction)
        where (fixed) reference = 0
        call nc_reduce(field, problem, reference)
        call nc_jvp(field, problem, x, direction, action, st)
        call require(st == 0, 'fitted JVP')
        error = norm2(action - reference)/max(tiny(1.0_dp), norm2(reference))
        call require(error < 1.0e-10_dp, 'fitted JVP/assembled')
        h = 1.0e-4_dp
        call nc_residual(field, problem, x + h*direction, plus, st)
        call require(st == 0, 'positive fitted perturbation')
        call nc_residual(field, problem, x - h*direction, minus, st)
        call require(st == 0, 'negative fitted perturbation')
        error = norm2((plus - minus)/(2*h) - action)/norm2(action)
        call require(error < 1.0e-6_dp, 'fitted JVP/residual finite difference')
        call nc_vjp(field, problem, x, weight, adjoint, st)
        call require(st == 0, 'fitted VJP')
        error = abs(dot_product(weight, action) - dot_product(adjoint, direction)) &
            /max(tiny(1.0_dp), abs(dot_product(weight, action)))
        call require(error < 1.0e-10_dp, 'fitted JVP/VJP pairing')
        print '(a,3es14.5)', 'signed current, pairing error, physical B2 mean ', &
            expected, error, physical_b2/volume
    end subroutine
end program
