program test_vrmhd_evolve
    !! BDF2 evolution gates of docs/vrmhd-section.md#time-evolution.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_vrmhd_parameters, only: vrmhd_parameters_t
    use milca_vrmhd_section, only: vrmhd_section_t, vrmhd_section_basis
    use milca_vrmhd_evolve, only: vrmhd_evolve_t, vrmhd_evolve_build, &
        vrmhd_evolve_expand, vrmhd_evolve_consistent
    use vrmhd_periodic_slab, only: periodic_slab_t, periodic_slab_build, &
        periodic_scalar, periodic_potential
    use fortfem_meridian_feec_n, only: meridian_load
    use fortnum_ode_bdf, only: bdf_options_t, bdf_state_t, bdf_report_t, bdf_init, &
        bdf_integrate_to
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortsparse, only: csc_t, csc_matmul, csc_transpose, sparse_solve_once, &
        csc_from_triplet, &
        fortsparse_status_t, status_ok
    use milca_vrmhd_section, only: vrmhd_section_jacobian, VRMHD_EDGE_SECOND
    use vrmhd_spectrum_support, only: spectrum, block_dofs
    use rm_mhd_oracles, only: uniform_waves, oracle_mu0, oracle_ok
    use milca_gs_symmetry, only: gs_symmetry_t, GS_TRANSLATION, GS_CHART_CARTESIAN
    use milca_vrmhd_section, only: vrmhd_section_build
    use milca_vrmhd_stationary, only: vrmhd_stationary_t, vrmhd_stationary_build, &
        vrmhd_stationary_residual, vrmhd_stationary_jacobian
    use milca_vrmhd_drive, only: vrmhd_drive_build, VRMHD_LOOP_VOLTAGE
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes
    implicit none
    !> Quadrature totals: kinetic, magnetic, internal energy, Ohmic and viscous
    !> dissipation, pressure work, and magnetic energy of B_x, B_y, B_z.
    integer, parameter :: KINETIC = 1, MAGNETIC = 2, INTERNAL = 3, OHMIC = 4, &
        VISCOUS = 5, COMPRESSION = 6, BX = 7, BY = 8, BZ = 9, NLEDGER = 9
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: background(3)
    character(32) :: gate
    call get_command_argument(1, gate)
    select case (trim(gate))
    case ('energy')
        call energy_gate([0.25_dp, 0.125_dp, 0.0625_dp, 0.03125_dp])
    case ('decay')
        call decay_gate()
    case ('steady')
        call steady_gate()
    case ('waves')
        call wave_gate([0.25_dp, 0.125_dp, 0.0625_dp, 0.03125_dp])
    case ('refactor')
        call refactor_gate()
    case default
        error stop 'unknown evolution gate'
    end select
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (ok) return
        print '(a)', message
        error stop 'RM2 evolution gate failed'
    end subroutine

    real(dp) function order(coarse, fine, ratio)
        real(dp), intent(in) :: coarse, fine, ratio
        order = log(abs(coarse)/abs(fine))/log(ratio)
    end function

    subroutine ledger(op, p, u, js, totals)
        !! Physical energies and powers at the production quadrature points,
        !! written from the weak-form definitions, independently of the leaves.
        type(vrmhd_section_t), intent(in) :: op
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: u(:), js(:, :, :)
        real(dp), intent(out) :: totals(NLEDGER)
        real(dp), allocatable :: basis(:, :)
        real(dp) :: x(2), weight, v(36), g(3, 3), s(3, 3), eta, nu, pressure
        integer :: t, q, status, i
        allocate (basis(36, op%nlocal))
        totals = 0
        do t = 1, op%fs%mm%topology%n_triangles
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                call require(status == 0, 'ledger basis')
                v = matmul(basis, u(op%map(:, t)))
                g = reshape(v(16:24), [3, 3])
                s = (g + transpose(g))/2
                do i = 1, 3
                    s(i, i) = s(i, i) - (g(1, 1) + g(2, 2) + g(3, 3))/3
                end do
                eta = p%eta0*(v(11)/p%t0)**(-1.5_dp)
                nu = p%nu0*(v(12)/p%t0)**(-1.5_dp)
                pressure = v(1)*(v(11) + v(12))/p%mi
                totals(KINETIC) = totals(KINETIC) + weight*v(1)*sum(v(2:4)**2)/2
                totals(MAGNETIC) = totals(MAGNETIC) + weight*sum(v(25:27)**2)/(2*p%mu0)
                totals(INTERNAL) = totals(INTERNAL) + weight*1.5_dp*pressure
                totals(OHMIC) = totals(OHMIC) &
                    + weight*eta*dot_product(v(8:10), v(8:10) - js(:, q, t))
                totals(VISCOUS) = totals(VISCOUS) + weight*2*nu*sum(s*s)
                totals(COMPRESSION) = totals(COMPRESSION) &
                    + weight*pressure*(g(1, 1) + g(2, 2) + g(3, 3))
                totals(BX:BZ) = totals(BX:BZ) + weight*v(25:27)**2/(2*p%mu0)
            end do
        end do
    end subroutine

    subroutine uniform_potential(x, a)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: a(3)
        a = cmplx([-background(3)*x(2)/2, background(1)*x(2), &
            background(3)*x(1)/2], 0.0_dp, dp)
    end subroutine

    subroutine background_lift(slab, rho, te, ti, b0, lift)
        !! Uniform state; the uniform-field potential is not periodic, so it
        !! enters through the lift on the full rectangular mesh.
        type(periodic_slab_t), intent(in) :: slab
        real(dp), intent(in) :: rho, te, ti, b0(3)
        real(dp), intent(out) :: lift(:)
        complex(dp), allocatable :: load(:)
        type(fortsparse_status_t) :: st
        integer :: n0, a0, n1, status
        n0 = slab%op%scalar_fs%n0; n1 = slab%op%fs%sp(1)%ndof; a0 = 6*n0
        lift = 0
        lift(:n0) = rho; lift(4*n0 + 1:5*n0) = te; lift(5*n0 + 1:6*n0) = ti
        if (all(b0 == 0)) return
        background = b0
        call meridian_load(slab%op%fs, 1, uniform_potential, load, status)
        call require(status == 0, 'uniform potential load')
        call sparse_solve_once(slab%op%fs%mass(1), real(load, dp), &
            lift(a0 + 1:a0 + n1), st)
        call require(status_ok(st), 'uniform potential projection')
    end subroutine

    subroutine run(slab, p, source, lift, y0, dt, steps, lag, history, final, &
            observable, series, atol, derivative, exact_shift)
        !! Constant-step BDF2 from a consistent start; ledger at every step.
        type(periodic_slab_t), intent(in), target :: slab
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: source(:, :, :), lift(:), y0(:), dt
        integer, intent(in) :: steps, lag
        real(dp), intent(out) :: history(NLEDGER, 0:steps)
        real(dp), intent(out), optional :: final(:)
        !> Row vector w^H B; series(k) = w^H B y(k dt).
        complex(dp), intent(in), optional :: observable(:)
        complex(dp), intent(out), optional :: series(0:)
        real(dp), intent(in), optional :: atol
        !> Consistent initial derivative for a consistent y0, reused across a ladder.
        real(dp), intent(in), optional :: derivative(:)
        logical, intent(in), optional :: exact_shift
        type(vrmhd_evolve_t) :: problem
        type(bdf_options_t) :: options
        type(bdf_state_t) :: state
        type(bdf_report_t) :: report
        type(fortnum_status_t) :: fstatus
        real(dp), allocatable :: y(:), ydot(:), u(:)
        integer :: status, k, newton, krylov
        call vrmhd_evolve_build(problem, slab%op, p, source, lift, slab%reduction, &
            status, lag=lag)
        call require(status == 0, 'evolution problem')
        ! One factor serves both the BDF1 start and the BDF2 shift 1.5/dt,
        ! unless each shift needs its own factor.
        problem%shift_tolerance = 0.5_dp
        if (present(exact_shift)) then
            if (exact_shift) problem%shift_tolerance = 0
        end if
        y = y0
        allocate (ydot(size(y)), u(slab%op%ndof))
        if (present(derivative)) then
            ydot = derivative
        else
            call vrmhd_evolve_consistent(problem, 0.0_dp, y, ydot, status)
            call require(status == 0, 'consistent initial current and derivative')
        end if
        options%adaptive = .false.
        options%initial_step = dt; options%max_step = dt
        options%rtol = 1.0e-8_dp; options%atol = 1.0e-10_dp
        if (present(atol)) options%atol = atol
        call bdf_init(state, 0.0_dp, y, options, fstatus, derivative=ydot)
        call require(fstatus%code == FORTNUM_OK, 'BDF initial state')
        call vrmhd_evolve_expand(problem, state%y, u)
        call ledger(slab%op, p, u, source(:3, :, :), history(:, 0))
        if (present(series)) series(0) = sum(observable*state%y)
        newton = 0; krylov = 0
        do k = 1, steps
            call bdf_integrate_to(problem, state, k*dt, options, report, fstatus)
            if (fstatus%code /= FORTNUM_OK .or. report%accepted /= 1) &
                print '(a,i5,8i7)', '  step report ', k, report%accepted, &
                report%rejected_newton, report%newton_iterations, &
                report%krylov_iterations, report%last_failure, &
                report%callback_status, report%linear_status, fstatus%code
            call require(fstatus%code == FORTNUM_OK, 'BDF step')
            call require(report%accepted == 1, 'one constant step per output')
            newton = newton + report%newton_iterations
            krylov = krylov + report%krylov_iterations
            call vrmhd_evolve_expand(problem, state%y, u)
            call ledger(slab%op, p, u, source(:3, :, :), history(:, k))
            if (present(series)) series(k) = sum(observable*state%y)
        end do
        print '(a,es10.3,a,i5,a,i5,a,i6,a,i4)', '  dt ', dt, ' steps ', steps, &
            ' newton ', newton, ' krylov ', krylov, ' factorizations ', &
            problem%factorizations
        if (present(final)) final = state%y
    end subroutine

    subroutine consistent_start(slab, p, source, lift, y, ydot)
        type(periodic_slab_t), intent(in), target :: slab
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: source(:, :, :), lift(:)
        real(dp), intent(inout) :: y(:)
        real(dp), intent(out) :: ydot(:)
        type(vrmhd_evolve_t) :: problem
        integer :: status
        call vrmhd_evolve_build(problem, slab%op, p, source, lift, slab%reduction, &
            status)
        call require(status == 0, 'start problem')
        call vrmhd_evolve_consistent(problem, 0.0_dp, y, ydot, status)
        call require(status == 0, 'consistent initial current and derivative')
    end subroutine

    pure real(dp) function trapezoid(f, dt) result(total)
        real(dp), intent(in) :: f(0:), dt
        total = dt*(sum(f) - (f(0) + f(ubound(f, 1)))/2)
    end function

    subroutine perturbed_state(slab, amplitude, y)
        type(periodic_slab_t), intent(in) :: slab
        real(dp), intent(in) :: amplitude
        real(dp), intent(out) :: y(:)
        integer :: ns, status
        ns = slab%ns
        y = 0
        call periodic_scalar(slab, rho_wave, y(:ns))
        call periodic_scalar(slab, vx_wave, y(ns + 1:2*ns))
        call periodic_scalar(slab, vy_wave, y(2*ns + 1:3*ns))
        call periodic_scalar(slab, vz_wave, y(3*ns + 1:4*ns))
        call periodic_scalar(slab, te_wave, y(4*ns + 1:5*ns))
        call periodic_scalar(slab, ti_wave, y(5*ns + 1:6*ns))
        call periodic_potential(slab, a_wave, y(6*ns + 1:6*ns + slab%na), status)
        call require(status == 0, 'periodic potential projection')
        y = amplitude*y
    end subroutine

    pure real(dp) function rho_wave(x)
        real(dp), intent(in) :: x(2)
        rho_wave = 0.5_dp*cos(x(1) + x(2))
    end function
    pure real(dp) function vx_wave(x)
        real(dp), intent(in) :: x(2)
        vx_wave = 0.3_dp*sin(x(2))
    end function
    pure real(dp) function vy_wave(x)
        real(dp), intent(in) :: x(2)
        vy_wave = cos(x(1))
    end function
    pure real(dp) function vz_wave(x)
        real(dp), intent(in) :: x(2)
        vz_wave = 0.5_dp*sin(x(1) + x(2))
    end function
    pure real(dp) function te_wave(x)
        real(dp), intent(in) :: x(2)
        te_wave = 0.2_dp*sin(x(1))
    end function
    pure real(dp) function ti_wave(x)
        real(dp), intent(in) :: x(2)
        ti_wave = 0.1_dp*cos(x(2))
    end function
    subroutine a_wave(x, a)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: a(3)
        a = cmplx([0.5_dp*cos(x(2)), 0.3_dp*sin(x(1)), 0.4_dp*cos(x(1) - x(2))], &
            0.0_dp, dp)
    end subroutine

    subroutine refactor_gate()
        !! A lagged factor is kept while FGMRES needs at most krylov_refactor
        !! products per linear solve and replaced at the next linearization
        !! once a solve needs more, even inside the lag and shift windows.
        type(periodic_slab_t), target :: slab
        type(vrmhd_parameters_t) :: p
        type(vrmhd_evolve_t) :: problem
        real(dp), allocatable :: lift(:), y(:), ydot(:), source(:, :, :), v(:), w(:)
        integer :: status, k
        call periodic_slab_build(4, 1, slab, status)
        call require(status == 0, 'periodic slab')
        p%mi = 1; p%mu0 = 1; p%t0 = 0.3_dp; p%eta0 = 0.02_dp; p%nu0 = 0.02_dp
        allocate (lift(slab%op%ndof), y(slab%nred), ydot(slab%nred), &
            source(6, slab%op%fs%nq, slab%op%fs%mm%topology%n_triangles), &
            v(slab%nred), w(slab%nred))
        source = 0
        call background_lift(slab, 1.0_dp, 0.3_dp, 0.3_dp, [0.6_dp, 0.0_dp, 0.8_dp], lift)
        call perturbed_state(slab, 0.1_dp, y)
        call vrmhd_evolve_build(problem, slab%op, p, source, lift, slab%reduction, &
            status, lag=huge(1))
        call require(status == 0, 'evolution problem')
        ydot = 0; v = 1
        call problem%linearize(0.0_dp, y, ydot, 10.0_dp, status)
        call require(status == 0 .and. problem%factorizations == 1, 'first factor')
        do k = 1, 8
            call problem%apply(v, w)
        end do
        call problem%linearize(0.0_dp, y, ydot, 10.0_dp, status)
        call require(status == 0 .and. problem%factorizations == 1, &
            'factor kept after a cheap solve')
        do k = 1, 9
            call problem%apply(v, w)
        end do
        call problem%linearize(0.0_dp, y, ydot, 10.0_dp, status)
        call require(status == 0 .and. problem%factorizations == 2, &
            'factor replaced after an expensive solve')
        print '(a,i0)', 'factorizations ', problem%factorizations
    end subroutine

    subroutine energy_gate(dts)
        !! Ideal drift and the dissipative/driven mechanical ledger over one
        !! constant-step ladder; thresholds: every consecutive order > 1.8.
        real(dp), intent(in) :: dts(:)
        type(periodic_slab_t), target :: slab
        type(vrmhd_parameters_t) :: ideal, resistive
        real(dp), allocatable :: lift(:), y0(:), source(:, :, :), history(:, :)
        real(dp), allocatable :: ydot0(:)
        real(dp) :: drift(size(dts), 2), defect(size(dts)), power(size(dts))
        real(dp), parameter :: final_time = 2.0_dp
        real(dp) :: x(2), weight, scale
        real(dp), allocatable :: basis(:, :)
        integer :: status, level, steps, t, q, case
        ! Semidiscrete conservation holds on any mesh; four cells keep it cheap.
        call periodic_slab_build(4, 1, slab, status)
        call require(status == 0, 'periodic slab')
        ideal%mi = 1; ideal%mu0 = 1; ideal%t0 = 0.3_dp
        resistive = ideal
        resistive%eta0 = 0.02_dp; resistive%nu0 = 0.02_dp
        resistive%chi_perpendicular = [0.01_dp, 0.01_dp]
        allocate (lift(slab%op%ndof), y0(slab%nred), basis(36, slab%op%nlocal), &
            source(6, slab%op%fs%nq, slab%op%fs%mm%topology%n_triangles))
        call background_lift(slab, 1.0_dp, 0.3_dp, 0.3_dp, [0.6_dp, 0.0_dp, 0.8_dp], &
            lift)
        call perturbed_state(slab, 0.1_dp, y0)
        allocate (ydot0(slab%nred))
        do case = 1, 2
            source = 0
            if (case == 2) then
                ! Prescribed current Js for the current-source drive.
                do t = 1, slab%op%fs%mm%topology%n_triangles
                    do q = 1, slab%op%fs%nq
                        call vrmhd_section_basis(slab%op, t, q, x, weight, basis, status)
                        source(:3, q, t) = [0.1_dp*sin(x(1)), 0.0_dp, 0.2_dp*cos(x(2))]
                    end do
                end do
            end if
            if (case == 1) call consistent_start(slab, ideal, source, lift, y0, ydot0)
            if (case == 2) call consistent_start(slab, resistive, source, lift, y0, ydot0)
            do level = 1, size(dts)
                steps = nint(final_time/dts(level))
                if (allocated(history)) deallocate (history)
                allocate (history(NLEDGER, 0:steps))
                if (case == 1) then
                    call run(slab, ideal, source, lift, y0, dts(level), steps, &
                        huge(1), history, derivative=ydot0)
                else
                    call run(slab, resistive, source, lift, y0, dts(level), steps, &
                        huge(1), history, derivative=ydot0)
                end if
                scale = sum(history(KINETIC:INTERNAL, 0))
                drift(level, case) = (sum(history(KINETIC:INTERNAL, steps)) - scale)/scale
                if (case == 2) then
                    power(level) = trapezoid(history(COMPRESSION, :) &
                        - history(OHMIC, :) - history(VISCOUS, :), dts(level))
                    defect(level) = sum(history(KINETIC:MAGNETIC, steps)) &
                        - sum(history(KINETIC:MAGNETIC, 0)) - power(level)
                    print '(a,es10.3,a,es14.6,a,es12.4,a,es12.4)', 'dissipative dt ', &
                        dts(level), ' mechanical power integral ', power(level), &
                        ' ledger defect ', defect(level), ' total drift ', &
                        drift(level, case)
                else
                    print '(a,es10.3,a,es12.4)', 'ideal dt ', dts(level), &
                        ' relative energy drift ', drift(level, case)
                end if
            end do
        end do
        do level = 2, size(dts)
            print '(a,i2,3(a,f7.3))', 'orders level ', level, &
                ' ideal drift ', order(drift(level - 1, 1), drift(level, 1), 2.0_dp), &
                ' dissipative drift ', order(drift(level - 1, 2), drift(level, 2), 2.0_dp), &
                ' ledger defect ', order(defect(level - 1), defect(level), 2.0_dp)
            call require(order(drift(level - 1, 1), drift(level, 1), 2.0_dp) > 1.8_dp, &
                'ideal energy drift order')
            call require(order(drift(level - 1, 2), drift(level, 2), 2.0_dp) > 1.8_dp, &
                'dissipative total energy drift order')
            call require(order(defect(level - 1), defect(level), 2.0_dp) > 1.8_dp, &
                'mechanical energy ledger order')
        end do
    end subroutine

    subroutine a_beltrami(x, a)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: a(3)
        a = cmplx([0.0_dp, sin(x(1)), cos(x(1))], 0.0_dp, dp)
    end subroutine

    subroutine decay_gate()
        !! Force-free Beltrami field B = B1 (0, sin x, cos x) decays at
        !! eta k^2/mu0 with k = 1. The exact BDF2 recurrence for that rate
        !! removes the time error, leaving the spatial error for the gate.
        type(periodic_slab_t), target :: slab
        type(vrmhd_parameters_t) :: p
        real(dp), allocatable :: lift(:), y0(:), source(:, :, :), history(:, :)
        real(dp), parameter :: final_time = 8.0_dp, amplitude = 1.0e-3_dp
        real(dp), parameter :: dt = 0.5_dp
        integer, parameter :: cells(2) = [4, 8]
        real(dp) :: rate(2), error(2, 2), exact, discrete
        integer :: status, level, steps, ns, component
        p%mi = 1; p%mu0 = 1; p%t0 = 1; p%eta0 = 0.05_dp
        exact = p%eta0/p%mu0
        steps = nint(final_time/dt)
        discrete = bdf2_decay_rate(exact, dt, steps)
        print '(a,es14.6,a,es14.6)', 'analytic rate ', exact, &
            ' after the exact BDF2 recurrence ', discrete
        do level = 1, 2
            call periodic_slab_build(cells(level), 2, slab, status)
            call require(status == 0, 'periodic slab')
            if (allocated(lift)) deallocate (lift, y0, source)
            allocate (lift(slab%op%ndof), y0(slab%nred), &
                source(6, slab%op%fs%nq, slab%op%fs%mm%topology%n_triangles))
            source = 0
            call background_lift(slab, 1.0_dp, 1.0_dp, 1.0_dp, [0.0_dp, 0.0_dp, 0.0_dp], &
                lift)
            ns = slab%ns
            y0 = 0
            call periodic_potential(slab, a_beltrami, y0(6*ns + 1:6*ns + slab%na), status)
            call require(status == 0, 'Beltrami projection')
            y0 = amplitude*y0
            if (allocated(history)) deallocate (history)
            allocate (history(NLEDGER, 0:steps))
            call run(slab, p, source, lift, y0, dt, steps, huge(1), history)
            do component = 1, 2
                rate(component) = -log(history(BX + component, steps) &
                    /history(BX + component, 0))/(2*final_time)
                error(component, level) = rate(component)/discrete - 1
            end do
            print '(a,i3,a,es10.3,a,2es14.6,a,2es11.3)', 'cells ', cells(level), &
                ' dt ', dt, ' decay rates B_y, B_z ', rate, &
                ' spatial relative errors ', error(:, level)
        end do
        do component = 1, 2
            print '(a,i2,a,f7.3)', 'spatial order component ', component, ' ', &
                order(error(component, 1), error(component, 2), 2.0_dp)
            call require(order(error(component, 1), error(component, 2), 2.0_dp) &
                > 1.8_dp, 'resistive decay spatial order')
        end do
        call require(maxval(abs(error(:, 2))) <= 1.0e-2_dp, 'resistive decay rate')
    end subroutine

    subroutine steady_gate()
        !! Loop-voltage Ohmic state in a conducting-wall box, translation
        !! symmetry with unit axial period: the stationary Newton solve with
        !! Coulomb multiplier is the reference; BDF2 in Weyl gauge from a
        !! perturbed start must relax to it in every gauge-invariant field.
        type(vrmhd_section_t), target :: op
        type(vrmhd_stationary_t) :: steady
        type(vrmhd_parameters_t) :: p
        type(vrmhd_evolve_t) :: problem
        type(mesh_2d_t) :: mesh
        type(meridian_mesh_t) :: mm
        type(csc_t) :: jac, selection
        type(fortsparse_status_t) :: st
        type(bdf_options_t) :: options
        type(bdf_state_t) :: state
        type(bdf_report_t) :: report
        type(fortnum_status_t) :: fstatus
        real(dp), allocatable :: u(:), boundary(:), load(:), source(:, :, :), js(:, :, :)
        real(dp), allocatable :: harmonic(:), r(:), jv(:), step(:), zero(:), xy(:, :)
        real(dp), allocatable :: lift(:), y(:), ydot(:), evolved(:), stationary(:)
        complex(dp), allocatable :: potential(:)
        integer, allocatable :: rows(:), cols(:)
        real(dp), parameter :: voltage = 2.0e-2_dp, final_time = 16.0_dp
        real(dp) :: weight, difference(3), initial(3)
        integer :: status, n0, n1, nn, a0, j0, i, it, nf, segment
        call mesh%create_rectangular(3, 3, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call meridian_mesh_straight(mesh, 1, mm, status)
        call require(status == 0, 'wall mesh')
        call vrmhd_section_build(mm, gs_symmetry_t(GS_TRANSLATION, &
            GS_CHART_CARTESIAN, 0.0_dp), 1, 1, op, status)
        call require(status == 0, 'wall section')
        call vrmhd_stationary_build(op, steady, status)
        call require(status == 0, 'stationary gauge')
        p%mi = 1; p%mu0 = 1; p%t0 = 1; p%eta0 = 0.2_dp; p%nu0 = 0.2_dp
        p%chi_parallel0 = [0.3_dp, 0.3_dp]; p%chi_perpendicular = [0.3_dp, 0.3_dp]
        p%d_parallel = 0.2_dp; p%d_perpendicular = 0.2_dp
        p%equilibration_rate = 0.5_dp
        n0 = op%scalar_fs%n0; n1 = op%fs%sp(1)%ndof; nn = op%fs%nn
        a0 = 6*n0; j0 = a0 + n1
        allocate (u(steady%ndof), boundary(op%ndof), js(3, op%fs%nq, &
            op%fs%mm%topology%n_triangles), harmonic(n1), load(op%ndof), &
            source(6, op%fs%nq, op%fs%mm%topology%n_triangles), r(steady%ndof), &
            jv(steady%ndof), step(steady%ndof), zero(steady%ndof))
        source = 0; js = 0; zero = 0
        call vrmhd_drive_build(op, VRMHD_LOOP_VOLTAGE, 1.0_dp, voltage, js, &
            harmonic, source, load, status)
        call require(status == 0, 'loop-voltage drive')
        u = 0
        u(:n0) = 1; u(4*n0 + 1:6*n0) = 1
        background = [0.0_dp, 0.0_dp, 1.0_dp]
        call meridian_load(op%fs, 1, uniform_potential, potential, status)
        call require(status == 0, 'guide-field load')
        call sparse_solve_once(op%fs%mass(1), real(potential, dp), u(a0 + 1:j0), st)
        call require(status_ok(st), 'guide-field projection')
        ! Wall Ohm's law: eta J_t = V h_t with h = z/L_z and wall Te = T0.
        do i = 1, op%fs%n0
            if (op%fixed(j0 + nn + i)) u(j0 + nn + i) = voltage/p%eta0
        end do
        boundary = u(:op%ndof)
        do it = 1, 12
            call vrmhd_stationary_residual(op, steady, p, u, zero, boundary, source, &
                load, r, jv, status)
            call require(status == 0, 'stationary residual')
            print '(a,i3,es12.4)', '  stationary Newton ', it, norm2(r)
            if (norm2(r) < 1.0e-11_dp) exit
            call vrmhd_stationary_jacobian(op, steady, p, u, boundary, source, jac, &
                status)
            call require(status == 0, 'stationary Jacobian')
            call sparse_solve_once(jac, -r, step, st)
            call require(status_ok(st), 'stationary solve')
            u = u + step
        end do
        call require(norm2(r) < 1.0e-11_dp, 'stationary Ohmic state')
        stationary = u(:op%ndof)
        print '(a,es12.4)', 'Coulomb multiplier max ', maxval(abs(u(op%ndof + 1:)))
        ! Free coefficients are the reduced state; the wall data are the lift.
        nf = count(.not. op%fixed)
        rows = pack([(i, i=1, op%ndof)], .not. op%fixed)
        cols = [(i, i=1, nf)]
        call csc_from_triplet(op%ndof, nf, rows, cols, [(1.0_dp, i=1, nf)], &
            selection, st)
        call require(status_ok(st), 'free-coefficient selection')
        lift = merge(boundary, 0.0_dp, op%fixed)
        y = stationary(rows)
        call meridian_lagrange_nodes(op%scalar_fs, xy)
        do i = 1, n0
            if (op%fixed(i)) cycle
            weight = 16*xy(1, i)*(1 - xy(1, i))*xy(2, i)*(1 - xy(2, i))
            ! Mass, momentum, both temperatures and the potential are perturbed.
            u(i) = stationary(i) + 0.05_dp*weight
            u(n0 + i) = stationary(n0 + i) + 0.02_dp*weight
            u(3*n0 + i) = stationary(3*n0 + i) - 0.03_dp*weight
            u(4*n0 + i) = stationary(4*n0 + i) + 0.04_dp*weight
            u(5*n0 + i) = stationary(5*n0 + i) - 0.02_dp*weight
        end do
        u(op%ndof + 1:) = 0
        do i = a0 + 1, j0
            if (.not. op%fixed(i)) u(i) = stationary(i) + 1.0e-3_dp*sin(1.7_dp*i)
        end do
        y = u(rows)
        call vrmhd_evolve_build(problem, op, p, source, lift, selection, status, &
            load=load, lag=huge(1))
        call require(status == 0, 'wall evolution problem')
        problem%shift_tolerance = 0.5_dp
        allocate (ydot(nf), evolved(op%ndof))
        call vrmhd_evolve_consistent(problem, 0.0_dp, y, ydot, status)
        call require(status == 0, 'consistent perturbed start')
        call vrmhd_evolve_expand(problem, y, evolved)
        call field_difference(op, evolved, stationary, initial)
        ! Loose error control, but Newton solves each step to 1e-4 of it.
        options%rtol = 1.0e-5_dp; options%atol = 1.0e-8_dp
        options%newton_tolerance = 1.0e-4_dp
        options%initial_step = 1.0e-2_dp
        call bdf_init(state, 0.0_dp, y, options, fstatus, derivative=ydot)
        call require(fstatus%code == FORTNUM_OK, 'steady BDF start')
        ! Segments only report progress; continuation keeps the BDF history.
        do segment = 1, 6
            call bdf_integrate_to(problem, state, final_time*2.0_dp**(segment - 6), &
                options, report, fstatus)
            call require(fstatus%code == FORTNUM_OK, 'steady BDF integration')
            call vrmhd_evolve_expand(problem, state%y, evolved)
            call field_difference(op, evolved, stationary, difference)
            print '(a,f7.3,a,6i6,es10.2,a,3es10.2)', 't ', state%t, &
                ' steps/rejections/Newton/Krylov/factors ', report%accepted, &
                report%rejected_error, report%rejected_newton, &
                report%newton_iterations, report%krylov_iterations, &
                problem%factorizations, state%previous_step, ' differences ', difference
            flush (6)
        end do
        print '(a,3es12.4)', 'initial relative differences fluid / J / B ', initial
        print '(a,3es12.4)', 'final relative differences fluid / J / B   ', difference
        call require(maxval(initial) > 1.0e-3_dp, 'perturbed start')
        call require(maxval(difference) <= 1.0e-7_dp, 'relaxation to stationary state')
    end subroutine

    subroutine field_difference(op, field, stationary, relative)
        !! Fluid coefficients, current coefficients and B at quadrature
        !! points; A itself differs by the Weyl-gauge gradient.
        type(vrmhd_section_t), intent(in) :: op
        real(dp), intent(in) :: field(:), stationary(:)
        real(dp), intent(out) :: relative(3)
        real(dp), allocatable :: basis(:, :)
        real(dp) :: x(2), weight, ve(36), vs(36), scale
        integer :: a0, j0, t, q, status
        a0 = 6*op%scalar_fs%n0; j0 = a0 + op%fs%sp(1)%ndof
        allocate (basis(36, op%nlocal))
        relative(1) = maxval(abs(field(:a0) - stationary(:a0))) &
            /maxval(abs(stationary(:a0)))
        relative(2) = maxval(abs(field(j0 + 1:) - stationary(j0 + 1:))) &
            /maxval(abs(stationary(j0 + 1:)))
        relative(3) = 0; scale = 0
        do t = 1, op%fs%mm%topology%n_triangles
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                call require(status == 0, 'difference basis')
                ve = matmul(basis, field(op%map(:, t)))
                vs = matmul(basis, stationary(op%map(:, t)))
                relative(3) = max(relative(3), maxval(abs(ve(25:27) - vs(25:27))))
                scale = max(scale, maxval(abs(vs(25:27))))
            end do
        end do
        relative(3) = relative(3)/scale
    end subroutine

    subroutine reduced_dense(slab, p, lift, source, shift, matrix)
        !! P^T J P of the production section Jacobian at ydot = 0.
        type(periodic_slab_t), intent(in) :: slab
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: lift(:), source(:, :, :), shift
        real(dp), intent(out) :: matrix(:, :)
        type(csc_t) :: jac, half, restriction, reduced
        type(fortsparse_status_t) :: st
        real(dp), allocatable :: zero(:)
        integer :: status, j, k
        allocate (zero(slab%op%ndof), source=0.0_dp)
        call vrmhd_section_jacobian(slab%op, p, lift, zero, source, shift, jac, status)
        call require(status == 0, 'wave linearization')
        call csc_matmul(jac, slab%reduction, half, st)
        call require(status_ok(st), 'wave reduction')
        call csc_transpose(slab%reduction, restriction, st)
        call require(status_ok(st), 'wave restriction')
        call csc_matmul(restriction, half, reduced, st)
        call require(status_ok(st), 'wave reduced Jacobian')
        matrix = 0
        do j = 1, reduced%ncol
            do k = reduced%col_ptr(j), reduced%col_ptr(j + 1) - 1
                matrix(reduced%row_idx(k), j) = reduced%val(k)
            end do
        end do
    end subroutine

    pure real(dp) function cosine_x(x)
        real(dp), intent(in) :: x(2)
        cosine_x = cos(x(1))
    end function

    subroutine wave_gate(dts)
        !! Small Alfven and fast waves along x on the periodic slab. The
        !! linearized pencil M ydot + K y = 0 has algebraic unknowns (zero rows
        !! and columns of M). Eliminating them leaves ydot_d = A y_d, and a left
        !! eigenvector w of A gives the scalar s = w^H y_d with s' = lambda s
        !! exactly, so the evolution must follow the scalar BDF2 recurrence for
        !! that lambda, and its frequency must converge to Im(lambda) at order
        !! two in dt.
        real(dp), intent(in) :: dts(:)
        type(periodic_slab_t), target :: slab
        type(vrmhd_parameters_t) :: p
        real(dp), allocatable :: lift(:), y0(:), source(:, :, :), history(:, :)
        real(dp), allocatable :: stiffness(:, :), mass(:, :), shifted(:, :)
        real(dp), allocatable :: ode(:, :), coupling(:, :), alphar(:), alphai(:)
        real(dp), allocatable :: left(:, :), work(:), ydot0(:)
        complex(dp), allocatable :: w(:), observable(:), series(:), predicted(:)
        complex(dp) :: lambda
        real(dp), parameter :: final_time = 4.0_dp, amplitude = 1.0e-4_dp
        ! The gate compares with the discrete eigenvalue, so a coarse mesh suffices.
        integer, parameter :: cells = 4
        real(dp) :: b0(3), exact(3), bloch(3), optical, overlap, best, omega_h
        real(dp) :: measured(size(dts)), error(size(dts)), recurrence, dummy(1, 1)
        integer, allocatable :: d(:), a(:), pivot(:)
        integer :: status, n, wave, j, chosen, level, steps, k, info, unmatched, ns
        integer :: nd
        interface
            subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
                import :: dp
                integer, intent(in) :: n, nrhs, lda, ldb
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
                integer, intent(out) :: ipiv(*), info
            end subroutine
            subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, &
                    work, lwork, info)
                import :: dp
                character, intent(in) :: jobvl, jobvr
                integer, intent(in) :: n, lda, ldvl, ldvr, lwork
                real(dp), intent(inout) :: a(lda, *)
                real(dp), intent(out) :: wr(*), wi(*), vl(ldvl, *), vr(ldvr, *), work(*)
                integer, intent(out) :: info
            end subroutine
        end interface
        b0 = sqrt(oracle_mu0)*[0.6_dp, 0.0_dp, 0.8_dp]
        call uniform_waves([1.0_dp, 0.0_dp, 0.0_dp], b0, 1.0_dp, 0.6_dp, &
            5.0_dp/3, exact, status)
        call require(status == oracle_ok, 'uniform-wave oracle')
        exact = sqrt(exact)
        call spectrum(cells, 1, 1, block_dofs(1, 1, VRMHD_EDGE_SECOND, 1), b0, exact, &
            .false., bloch, optical, unmatched, edge_kind=VRMHD_EDGE_SECOND)
        print '(a,3es20.12)', 'Bloch spectrum Alfven / fast / slow: ', bloch
        call periodic_slab_build(cells, 1, slab, status, &
            magnetic_unit=sqrt(oracle_mu0))
        call require(status == 0, 'periodic slab')
        p%mi = 1; p%mu0 = oracle_mu0; p%t0 = 0.3_dp
        n = slab%nred; ns = slab%ns
        allocate (lift(slab%op%ndof), y0(n), stiffness(n, n), mass(n, n), &
            shifted(n, n), observable(n), ydot0(n), &
            source(6, slab%op%fs%nq, slab%op%fs%mm%topology%n_triangles))
        source = 0
        call background_lift(slab, 1.0_dp, 0.3_dp, 0.3_dp, b0, lift)
        call reduced_dense(slab, p, lift, source, 0.0_dp, stiffness)
        call reduced_dense(slab, p, lift, source, 1.0_dp, shifted)
        ! B = -M, so ydot_d = B_dd^{-1} (K_dd - K_da K_aa^{-1} K_ad) y_d.
        mass = -(shifted - stiffness)
        ! The QZ eigenvalues of the singular pencil carry roundoff amplified by
        ! its infinite eigenvalues (|Re lambda| up to 1e-10 |lambda|, varying
        ! with the thread count for identical matrices); the Schur complement
        ! on the differential unknowns leaves a well-conditioned problem.
        d = pack([(j, j=1, n)], any(mass /= 0, dim=2))
        a = pack([(j, j=1, n)], .not. any(mass /= 0, dim=2))
        nd = size(d)
        call require(count(any(mass /= 0, dim=1)) == nd, &
            'algebraic rows and columns of the time mass agree')
        call require(all(d == pack([(j, j=1, n)], any(mass /= 0, dim=1))), &
            'algebraic rows and columns of the time mass agree')
        allocate (ode(nd, nd), alphar(nd), alphai(nd), left(nd, nd), w(nd), &
            work(16*nd), pivot(n))
        shifted(:size(a), :size(a)) = stiffness(a, a)
        coupling = stiffness(a, d)
        call dgesv(size(a), nd, shifted, n, pivot, coupling, size(a), info)
        call require(info == 0, 'algebraic block solve')
        ode = stiffness(d, d) - matmul(stiffness(d, a), coupling)
        shifted(:nd, :nd) = mass(d, d)
        call dgesv(nd, nd, shifted, n, pivot, ode, nd, info)
        call require(info == 0, 'differential mass solve')
        call dgeev('V', 'N', nd, ode, nd, alphar, alphai, left, nd, dummy, 1, work, &
            size(work), info)
        call require(info == 0, 'dense eigensolve')
        do wave = 1, 2
            y0 = 0
            ! The Alfven wave moves v_y; the fast wave along B0's normal moves v_x.
            j = 3 - wave
            call periodic_scalar(slab, cosine_x, y0(j*ns + 1:(j + 1)*ns))
            y0 = amplitude*y0
            call consistent_start(slab, p, source, lift, y0, ydot0)
            best = 0; chosen = 0
            do j = 1, nd
                if (alphai(j) <= 0) cycle
                lambda = cmplx(alphar(j), alphai(j), dp)
                if (abs(aimag(lambda)/exact(wave) - 1) > 0.2_dp) cycle
                w = cmplx(left(:, j), left(:, j + 1), dp)
                overlap = abs(dot_product(w, cmplx(y0(d), 0.0_dp, dp))) &
                    /(norm2(abs(w))*norm2(y0(d)))
                if (overlap > best) then
                    best = overlap; chosen = j
                end if
            end do
            call require(chosen > 0, 'excited discrete wave eigenvalue')
            lambda = cmplx(alphar(chosen), alphai(chosen), dp)
            omega_h = aimag(lambda)
            w = cmplx(left(:, chosen), left(:, chosen + 1), dp)
            ! s = sum(observable*y) = w^H y_d.
            observable = 0
            observable(d) = conjg(w)
            print '(a,i2,a,2es20.12,a,f8.5)', 'wave ', wave, ' lambda ', lambda, &
                ' overlap ', best
            print '(a,es12.4)', '  relative difference to Bloch spectrum ', &
                omega_h/bloch(wave) - 1
            call require(abs(real(lambda))/abs(lambda) <= 1.0e-10_dp, &
                'ideal discrete wave is undamped')
            call require(abs(omega_h/bloch(wave) - 1) <= 1.0e-8_dp, &
                'evolution pencil matches the Bloch spectrum')
            do level = 1, size(dts)
                steps = nint(final_time/dts(level))
                if (allocated(history)) deallocate (history, series, predicted)
                allocate (history(NLEDGER, 0:steps), series(0:steps), &
                    predicted(0:steps))
                call run(slab, p, source, lift, y0, dts(level), steps, huge(1), &
                    history, observable=observable, series=series, atol=1.0e-11_dp, &
                    derivative=ydot0, exact_shift=.true.)
                predicted(0) = series(0)
                predicted(1) = series(0)/(1 - lambda*dts(level))
                do k = 2, steps
                    predicted(k) = (2*predicted(k - 1) - predicted(k - 2)/2) &
                        /(1.5_dp - lambda*dts(level))
                end do
                recurrence = maxval(abs(series - predicted))/maxval(abs(series))
                measured(level) = 0
                do k = 1, steps
                    measured(level) = measured(level) + aimag(log(series(k)/series(k - 1)))
                end do
                measured(level) = measured(level)/final_time
                error(level) = measured(level)/omega_h - 1
                print '(a,es10.3,a,es20.12,a,es11.3,a,es11.3)', '  dt ', dts(level), &
                    ' frequency ', measured(level), ' relative error ', error(level), &
                    ' recurrence defect ', recurrence
                call require(recurrence <= 1.0e-6_dp, 'BDF2 modal recurrence')
                if (level == 1) cycle
                print '(a,f7.3)', '  order ', order(error(level - 1), error(level), &
                    dts(level - 1)/dts(level))
                call require(order(error(level - 1), error(level), &
                    dts(level - 1)/dts(level)) > 1.8_dp, 'wave frequency order in dt')
            end do
        end do
    end subroutine

    real(dp) function bdf2_decay_rate(rate, dt, steps) result(effective)
        !! Backward Euler start, then constant-step BDF2, for a' = -rate a.
        real(dp), intent(in) :: rate, dt
        integer, intent(in) :: steps
        real(dp) :: previous, current, next
        integer :: k
        previous = 1
        current = 1/(1 + rate*dt)
        do k = 2, steps
            next = (2*current - previous/2)/(1.5_dp + rate*dt)
            previous = current; current = next
        end do
        effective = -log(current)/(steps*dt)
    end function
end program test_vrmhd_evolve
