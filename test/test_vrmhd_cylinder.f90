program test_vrmhd_cylinder
    !! Helically symmetric straight-cylinder flux-pumping case
    !! (docs/vrmhd-cylinder.md): profiles and the stationary 2D Ohmic state.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use milca_vrmhd_cylinder
    use milca_vrmhd_parameters, only: vrmhd_parameters_t
    use milca_vrmhd_section, only: vrmhd_section_t, vrmhd_section_build, &
        vrmhd_section_basis
    use milca_vrmhd_stationary, only: vrmhd_stationary_t, vrmhd_stationary_build, &
        vrmhd_stationary_residual, vrmhd_stationary_jacobian
    use milca_vrmhd_drive, only: vrmhd_drive_build, VRMHD_LOOP_VOLTAGE_SOURCE
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use milca_vrmhd_steady, only: vrmhd_steady_solve, vrmhd_steady_options_t, &
        vrmhd_steady_report_t, vrmhd_steady_spectrum
    use milca_vrmhd_evolve, only: vrmhd_evolve_t, vrmhd_evolve_build, &
        vrmhd_evolve_consistent, vrmhd_evolve_expand
    use fortsparse, only: csc_t, csc_from_triplet, fortsparse_status_t, status_ok
    use fortnum_ode_bdf, only: bdf_options_t, bdf_state_t, bdf_report_t, bdf_init, &
        bdf_integrate_to
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none
    character(64) :: gate, arg
    integer :: rings = 16, degree = 1
    real(dp) :: dtau = 0, strength = 0, t_end = 1.0e3_dp, dt_out = 10
    real(dp) :: beta_poloidal = 1, seed = 1.0e-6_dp
    integer :: frozen_eta = 0, k, e, check = 0
    real(dp) :: chi_ratio = 1.0e6_dp
    call get_command_argument(1, gate)
    ! Run options as key=value; the registered tests use the defaults.
    do k = 2, command_argument_count()
        call get_command_argument(k, arg)
        e = index(arg, '=')
        if (e == 0) error stop 'arguments are key=value'
        select case (arg(:e - 1))
        case ('rings'); read (arg(e + 1:), *) rings
        case ('degree'); read (arg(e + 1:), *) degree
        case ('dtau'); read (arg(e + 1:), *) dtau
        case ('cluster'); read (arg(e + 1:), *) strength
        case ('t_end'); read (arg(e + 1:), *) t_end
        case ('dt_out'); read (arg(e + 1:), *) dt_out
        case ('frozen'); read (arg(e + 1:), *) frozen_eta
        case ('beta_p'); read (arg(e + 1:), *) beta_poloidal
        case ('seed'); read (arg(e + 1:), *) seed
        case ('chi_ratio'); read (arg(e + 1:), *) chi_ratio
        case ('check'); read (arg(e + 1:), *) check
        case default; error stop 'unknown argument'
        end select
    end do
    select case (trim(gate))
    case ('profiles')
        call profiles_gate()
    case ('stationary')
        call stationary_gate(rings, .false.)
    case ('linear')
        call stationary_gate(rings, .true.)
    case ('evolve')
        call stationary_gate(rings, .false., .true.)
    case default
        error stop 'choose profiles or stationary'
    end select
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (ok) return
        print '(a)', message
        error stop 'cylinder gate failed'
    end subroutine

    subroutine profiles_gate()
        !! Independent implementation of the JOREK cylinder case generator; the
        !! reference numbers are its output for the same case and grid.
        type(vrmhd_cylinder_case_t) :: c
        type(vrmhd_cylinder_profile_t) :: prof
        real(dp) :: v(CYL_NVALUE), balance, worst
        integer :: status, i
        c%grid = 2001
        call vrmhd_cylinder_profiles(c, prof, status)
        call require(status == 0, 'profiles')
        print '(a,4es24.16)', 'dT, T(0), q_a, q(0.3) [Bz = 1]: ', prof%t(1) &
            /(c%edge_temperature + 1), prof%t(1), 1/(c%r0*prof%btheta_over_r(c%grid)), &
            1/(c%r0*prof%btheta_over_r(601))
        call require(abs(prof%t(1)/(c%edge_temperature + 1) &
            - 0.001138657913975527_dp) < 1.0e-12_dp*0.0011_dp, 'make_case dT')
        call require(abs(1/(c%r0*prof%btheta_over_r(c%grid)) &
            - 3.3519614264670703_dp) < 1.0e-12_dp, 'make_case q_a')
        call require(abs(1/(c%r0*prof%btheta_over_r(601)) &
            - 1.0325368489983022_dp) < 1.0e-12_dp, 'make_case q(0.3)')
        ! Radial force balance and Ohm's law of the interpolated profiles.
        worst = 0
        do i = 2, 99
            call vrmhd_cylinder_value(prof, 0.01_dp*i, v)
            balance = abs(v(CYL_JTHETA_R)*v(CYL_BZ) - v(CYL_DT_R) &
                - v(CYL_BTHETA_R)*v(CYL_JZ))
            worst = max(worst, balance)
        end do
        print '(a,es12.4)', 'force-balance defect: ', worst
        call require(worst < 1.0e-12_dp, 'force balance')
        print '(a,es24.16)', 'axial electric field: ', prof%electric_field
    end subroutine

    subroutine stationary_gate(nr, linear, evolve)
        !! Newton solve of the stationary 2D state from the interpolated
        !! Ohmic profiles; with `linear`, the m = 1 eigenvalue about it with
        !! frozen temperatures (frozen eta), against the reduced-MHD oracle.
        integer, intent(in) :: nr
        logical, intent(in) :: linear
        logical, intent(in), optional :: evolve
        type(vrmhd_cylinder_case_t) :: c
        type(vrmhd_cylinder_profile_t) :: prof
        type(meridian_mesh_t) :: mm
        type(vrmhd_section_t) :: op
        type(vrmhd_stationary_t) :: steady
        type(vrmhd_parameters_t) :: p
        type(vrmhd_steady_options_t) :: options
        type(vrmhd_steady_report_t) :: report
        real(dp), allocatable :: radii(:), u(:), u0(:), residuals(:)
        real(dp), allocatable :: source(:, :, :), load(:), harmonic(:)
        complex(dp), allocatable :: values(:), vectors(:, :)
        logical, allocatable :: frozen(:)
        real(dp) :: r1, change_te, speed
        integer(int64) :: t0, t1, rate
        integer :: status, i, info
        c%grid = 20001
        c%beta_poloidal = beta_poloidal
        if (linear) c%beta_poloidal = 0.01_dp
        c%chi_ratio = chi_ratio
        call vrmhd_cylinder_profiles(c, prof, status)
        call require(status == 0, 'profiles')
        r1 = q_one_radius(c, prof)
        if (strength > 0) then
            radii = vrmhd_cylinder_radii(nr, r1, strength, 0.3_dp)
        else
            radii = [(real(i, dp)/nr, i=1, nr)]
        end if
        call vrmhd_cylinder_mesh(radii, 8, 3.0_dp, mm, status)
        call require(status == 0, 'polar mesh')
        call vrmhd_section_build(mm, vrmhd_cylinder_symmetry(c), degree + 1, 1, op, &
            status, scalar_degree=degree, velocity_degree=degree + 1, &
            coordinate_coframe=.true.)
        call require(status == 0, 'section')
        call vrmhd_stationary_build(op, steady, status)
        call require(status == 0, 'gauge')
        call vrmhd_cylinder_parameters(c, prof, p)
        allocate (source(6, op%fs%nq, op%fs%mm%topology%n_triangles), &
            load(op%ndof), harmonic(op%fs%sp(1)%ndof), u(steady%ndof))
        call vrmhd_cylinder_sources(op, c, prof, source, status)
        call require(status == 0, 'sources')
        call vrmhd_drive_build(op, VRMHD_LOOP_VOLTAGE_SOURCE, prof%axial_period, &
            prof%electric_field*prof%axial_period, source(:3, :, :), harmonic, &
            source, load, status)
        call require(status == 0, 'drive')
        allocate (u0(op%ndof))
        call vrmhd_cylinder_state(op, prof, u0, status)
        call require(status == 0, 'state')
        call field_check(op, prof, u0)
        print '(a,f8.5,a,i0,a,i0,a,i0)', 'q = 1 at r = ', r1, '; triangles ', &
            mm%topology%n_triangles, ' unknowns ', steady%ndof, ' gauge ', steady%ngauge
        u = 0; u(:op%ndof) = u0
        options%dtau = dtau; options%verbose = .true.; options%rtol = 1.0e-12_dp
        options%growth = 4
        call system_clock(t0, rate)
        call vrmhd_steady_solve(op, steady, p, u0, source, load, u, options, report)
        call system_clock(t1)
        call require(report%status == 0, 'steady solve')
        print '(a,f10.2,a,i0,a,i0)', 'steady wall seconds ', real(t1 - t0, dp)/rate, &
            ' iterations ', report%iterations, ' factorizations ', report%factorizations
        change_te = maxval(abs(u(op%offset(5) + 1:op%offset(6)) &
            - u0(op%offset(5) + 1:op%offset(6))))/(prof%temperature_axis/2)
        speed = maxval(abs(u(op%offset(2) + 1:op%offset(5))))
        print '(a,es14.6)', 'max change of Te / Te(0): ', change_te
        print '(a,es14.6)', 'max change of rho: ', maxval(abs(u(op%offset(1) + 1: &
            op%offset(2)) - u0(op%offset(1) + 1:op%offset(2))))
        print '(a,es14.6)', 'max velocity: ', speed
        call report_state(op, c, prof, u0, 'interpolated profiles')
        call report_state(op, c, prof, u(:op%ndof), 'stationary state')
        if (check == 1 .and. .not. linear) then
            ! Isotropic conduction: the discrete state differs from the
            ! interpolated profiles by discretisation error only.
            call require(change_te < 1.0e-3_dp, 'stationary temperature')
            call require(speed < 1.0e-5_dp, 'stationary state at rest')
        end if
        if (.not. (linear .or. present(evolve))) return
        allocate (frozen(op%ndof), source=.false.)
        if (linear .or. frozen_eta == 1) frozen(op%offset(5) + 1:op%offset(7)) = .true.
        allocate (values(4), vectors(steady%ndof, 4), residuals(4))
        call system_clock(t0)
        call vrmhd_steady_spectrum(op, steady, p, u, u0, source, frozen, 1.0e-3_dp, &
            values, vectors, residuals, info)
        call system_clock(t1)
        print '(a,i0,a,f10.2)', 'spectrum info ', info, ' wall seconds ', &
            real(t1 - t0, dp)/rate
        do i = 1, size(values)
            print '(a,2es16.8,es10.2)', '  lambda, residual: ', values(i), residuals(i)
        end do
        call require(info == 0, 'eigenpairs converged')
        if (present(evolve)) then
            ! The unstable eigenvector seeds the evolution without exciting the
            ! fast waves an arbitrary flow would.
            i = maxloc(real(values), 1)
            call evolve_run(op, c, p, source, load, u(:op%ndof), &
                real(vectors(:op%ndof, i)))
            return
        end if
        ! Oracle: Chebyshev eigenvalue solver of visco-resistive reduced MHD of
        ! the same q(r), eta(r) and nu (beta_p = 0.01, N = 300): 7.218665e-4.
        print '(a,f9.5)', 'growth rate / reduced-MHD oracle: ', &
            maxval(real(values))/7.218665e-4_dp
        if (check == 1) call require(abs(maxval(real(values))/7.218665e-4_dp - 1) &
            < 0.05_dp, 'linear m = 1 growth rate')
    end subroutine

    subroutine evolve_run(op, c, p, source, load, u0, mode)
        !! BDF2 evolution from the stationary 2D state seeded with the unstable
        !! m = 1 eigenvector, scaled to a peak velocity `seed`; frozen_eta = 1
        !! keeps both temperatures at their profile.
        type(vrmhd_section_t), intent(in), target :: op
        type(vrmhd_cylinder_case_t), intent(in) :: c
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: source(:, :, :), load(:), u0(:), mode(:)
        type(vrmhd_evolve_t) :: problem
        type(bdf_options_t) :: options
        type(bdf_state_t) :: state
        type(bdf_report_t) :: report
        type(fortnum_status_t) :: st
        type(csc_t) :: reduction
        type(fortsparse_status_t) :: sparse_status
        real(dp), allocatable :: lift(:), y(:), derivative(:), u(:)
        real(dp), parameter :: edges(2) = [0.09_dp, 0.11_dp]
        real(dp) :: q(2, 1), energy(2), t_next
        logical, allocatable :: free(:)
        integer, allocatable :: kept(:)
        integer :: i, status, newton, steps, nfree
        integer(int64) :: t0, t1, rate
        allocate (free(op%ndof))
        free = .not. op%fixed
        if (frozen_eta == 1) free(op%offset(5) + 1:op%offset(7)) = .false.
        kept = pack([(i, i=1, op%ndof)], free)
        nfree = size(kept)
        call csc_from_triplet(op%ndof, nfree, kept, [(i, i=1, nfree)], &
            [(1.0_dp, i=1, nfree)], reduction, sparse_status)
        call require(status_ok(sparse_status), 'reduction')
        lift = u0 + seed*mode/maxval(abs(mode(op%offset(2) + 1:op%offset(5))))
        lift(op%offset(5) + 1:op%offset(7)) = merge(u0(op%offset(5) + 1:op%offset(7)), &
            lift(op%offset(5) + 1:op%offset(7)), frozen_eta == 1)
        y = lift(kept)
        lift(kept) = 0
        call vrmhd_evolve_build(problem, op, p, source, lift, reduction, status, &
            load=load, lag=8)
        call require(status == 0, 'evolution build')
        problem%shift_tolerance = 0.5_dp
        allocate (derivative(nfree))
        call vrmhd_evolve_consistent(problem, 0.0_dp, y, derivative, status)
        call require(status == 0, 'consistent start')
        options%adaptive = .true.; options%initial_step = 1
        options%rtol = 1.0e-4_dp; options%atol = 1.0e-12_dp
        options%max_step = dt_out
        call bdf_init(state, 0.0_dp, y, options, st, derivative=derivative)
        call require(st%code == FORTNUM_OK, 'BDF start')
        print '(a)', 'time, q(0.1) [Bz=1], q(0.1), E_B(m=1), E_K(m=1), '// &
            'steps, newton, factorizations, wall seconds'
        newton = 0; steps = 0
        call system_clock(t0, rate)
        t_next = 0
        allocate (u(op%ndof))
        do
            call vrmhd_evolve_expand(problem, state%y, u)
            call vrmhd_cylinder_diagnostics(op, c%r0, u, edges, q, energy, status)
            call require(status == 0, 'diagnostics')
            call system_clock(t1)
            print '(f12.2,2f11.6,2es13.5,3i8,f10.1)', state%t, q(:, 1), energy, &
                steps, newton, problem%factorizations, real(t1 - t0, dp)/rate
            flush (6)
            if (state%t >= t_end) exit
            t_next = t_next + dt_out
            call bdf_integrate_to(problem, state, t_next, options, report, st)
            steps = steps + report%accepted
            newton = newton + report%newton_iterations
            if (st%code /= FORTNUM_OK) print *, report
            call require(st%code == FORTNUM_OK, 'BDF step')
        end do
    end subroutine

    real(dp) function q_one_radius(c, prof) result(r1)
        type(vrmhd_cylinder_case_t), intent(in) :: c
        type(vrmhd_cylinder_profile_t), intent(in) :: prof
        real(dp) :: lo, hi
        integer :: it
        lo = 0.01_dp; hi = 0.9_dp
        do it = 1, 100
            r1 = (lo + hi)/2
            if (vrmhd_cylinder_q(prof, c, r1) < 1) then
                lo = r1
            else
                hi = r1
            end if
        end do
    end function

    subroutine report_state(op, c, prof, u, label)
        type(vrmhd_section_t), intent(in) :: op
        type(vrmhd_cylinder_case_t), intent(in) :: c
        type(vrmhd_cylinder_profile_t), intent(in) :: prof
        real(dp), intent(in) :: u(:)
        character(*), intent(in) :: label
        real(dp), parameter :: edges(6) = [0.09_dp, 0.11_dp, 0.25_dp, 0.27_dp, &
            0.49_dp, 0.51_dp]
        real(dp) :: q(2, 5), energy(2)
        integer :: status, b
        call vrmhd_cylinder_diagnostics(op, c%r0, u, edges, q, energy, status)
        call require(status == 0, 'diagnostics')
        print '(a)', label
        do b = 1, 5, 2
            print '(a,f6.3,3f12.7)', '  r, q(Bz=1), q, profile q(Bz=1): ', &
                (edges(b) + edges(b + 1))/2, q(:, b), &
                1/(c%r0*vrmhd_cylinder_value_at(prof, (edges(b) + edges(b + 1))/2))
        end do
        print '(a,2es12.4)', '  energy of the m = 1 part (B, v): ', energy
    end subroutine

    real(dp) function vrmhd_cylinder_value_at(prof, r) result(btheta_over_r)
        type(vrmhd_cylinder_profile_t), intent(in) :: prof
        real(dp), intent(in) :: r
        real(dp) :: v(CYL_NVALUE)
        call vrmhd_cylinder_value(prof, r, v)
        btheta_over_r = v(CYL_BTHETA_R)
    end function

    subroutine field_check(op, prof, u)
        type(vrmhd_section_t), intent(in) :: op
        type(vrmhd_cylinder_profile_t), intent(in) :: prof
        real(dp), intent(in) :: u(:)
        real(dp) :: basis(36, op%nlocal), x(2), weight, v(36), w(CYL_NVALUE), e(3)
        real(dp) :: ea, eb, ej
        integer :: t, q, status
        ea = 0; eb = 0; ej = 0
        do t = 1, op%fs%mm%topology%n_triangles
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                v = matmul(basis, u(op%map(:, t)))
                call vrmhd_cylinder_value(prof, norm2(x), w)
                e = [-x(2)*w(CYL_ATHETA_R), x(1)*w(CYL_ATHETA_R), -w(CYL_PSI)]
                ea = max(ea, maxval(abs(v(5:7) - e)))
                e = [-x(2)*w(CYL_BTHETA_R), x(1)*w(CYL_BTHETA_R), w(CYL_BZ)]
                eb = max(eb, maxval(abs(v(25:27) - e)))
                if (t == 1 .and. q == 1) print '(a,6es12.4)', 'B, B exact ', v(25:27), e
                e = [-x(2)*w(CYL_JTHETA_R), x(1)*w(CYL_JTHETA_R), w(CYL_JZ)]
                ej = max(ej, maxval(abs(v(8:10) - e)))
            end do
        end do
        print '(a,3es12.4)', 'interpolant max error A, B, J: ', ea, eb, ej
    end subroutine
end program test_vrmhd_cylinder
