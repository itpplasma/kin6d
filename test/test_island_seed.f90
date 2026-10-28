program test_island_seed
    ! Energy-criterion island seeding (docs/island-seed.md).
    !
    ! Periodic cylinder R0 = 3, radius 1, B_z = 1, prescribed nested
    ! background iota(r) = 0.4 + (0.1/0.5625) r^2 (q = 2 at the mesh ring
    ! r = 3/4) and an n = 1, m = 2 boundary trace of amplitude eps. The
    ! relaxed vacuum mode (nn_solve) has an island at q = 2. The "converged
    ! nested state" is the independent shielded minimiser: least magnetic
    ! energy with the same trace and zero resonant B_r on r = 3/4 (one
    ! Lagrange multiplier; the shielding sheet lies on element edges).
    !   seed   : one resonance; G symmetric positive definite; model exact
    !            against the independent relax_state energy at three
    !            amplitudes; W(shielded) > W(seeded) >= W(vacuum); pendulum
    !            width against field-line tracing of the seeded field for
    !            eps and eps/4 (width ratio exactly 1/2); hand-off to the
    !            NN1 Newton solve (returns the vacuum) and to the ideal
    !            metriplectic relaxation (energy decreases).
    !   none   : iota range without n/m (m <= 4): no seed, state unchanged;
    !            a resonant but axisymmetric state: zero amplitude.
    !   ladder : predicted width on 8/16/32 rings converges.
    !   leaves : generated leaf identities, JVP vs differences, JVP/VJP.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build, flux_mesh_locate
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_build, nn_report_t, &
        nn_solve, nn_linear_solve, nn_velocity
    use kin6d_nonnested_relax, only: relax_system_t, relax_state_t, relax_options_t, &
        relax_report_t, relax_prepare, relax_state, relax_pack_modes, relax_run, &
        relax_free, RELAX_IDEAL
    use kin6d_island_seed, only: seed_options_t, seed_report_t, island_seed, &
        seed_resonant_br, seed_resonance_t, island_resonances, &
        seed_circle_quadrature
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes, meridian_op_complex, &
        meridian_stiffness, meridian_point_basis, meridian_qpoint, &
        meridian_real_to_complex
    use fortsparse, only: csc_t, csc_z_t, csc_matvec, csc_from_triplet, &
        fortsparse_status_t, status_ok
    use kin6d_island_seed_leaves, only: circle_rotation, circle_rotation_jvp, &
        circle_rotation_vjp, seed_potential, seed_potential_jvp, &
        seed_potential_vjp, pendulum_width, pendulum_width_jvp, pendulum_width_vjp
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 3.0_dp, b0 = 1.0_dp
    integer, parameter :: m_res = 2, n_res = 1
    character(16) :: kind
    real(dp) :: iota0, iota1, eps
    ! Seed-case state shared with its helpers (internal procedures cannot nest).
    type(relax_system_t) :: rs
    complex(dp), allocatable :: a0(:), xs(:)
    call get_command_argument(1, kind)
    select case (trim(kind))
    case ('seed')
        call seed_case()
    case ('none')
        call none_case()
    case ('ladder')
        call ladder_case()
    case ('leaves')
        call leaves_case()
    case default
        error stop 'usage: test_island_seed seed|none|ladder|leaves'
    end select
    print '(a)', 'PASS'
contains
    subroutine seed_case()
        type(free_b_field_t) :: field, vacuum
        type(seed_options_t) :: opt
        type(seed_report_t) :: rep, rep4
        type(relax_state_t) :: st_s
        type(relax_options_t) :: ro
        type(relax_report_t) :: rr
        type(nn_report_t) :: nrep
        complex(dp), allocatable :: xv(:), an(:, :), sol(:), an4(:, :)
        complex(dp), allocatable :: delta(:)
        real(dp), allocatable :: a(:)
        real(dp) :: rs_exact, e(0:2), model, defect, dbr_vac, theta, wtrace, wtrace4
        real(dp) :: eseed, ratio
        integer :: st, t
        iota0 = 0.4_dp; iota1 = 0.1_dp/0.5625_dp; eps = 1.0e-4_dp
        rs_exact = sqrt((0.5_dp - iota0)/iota1)
        call fixture(16, field, a0, xv, xs, st)
        if (st /= 0) error stop 'fixture'
        vacuum = field
        call set_mode(vacuum, xv)
        allocate (an(field%fs%nn, 1)); an(:, 1) = xs
        opt%m_max = 4
        call island_seed(field, opt, rep, an)
        print '(a,i0,1x,a)', 'status ', rep%status, trim(rep%message)
        if (rep%status /= 0) error stop 'seed failed'
        if (rep%count /= 1) error stop 'expected exactly one resonance'
        associate (r => rep%res(1))
            print '(a,2i3,4es14.6)', 'm n u iota shear support ', r%m, r%n, r%u, &
                r%iota, r%shear, r%support
            print '(a,2es14.6)', 'circularity, bphi ', r%circularity, r%bphi
            if (r%m /= m_res .or. r%n /= n_res) error stop 'resonance (m, n)'
            if (abs(r%u - rs_exact) > 1.0e-2_dp) error stop 'resonance radius'
            if (abs(r%iota - 0.5_dp) > 1.0e-9_dp) error stop 'resonance iota'
            if (abs(r%shear - 2*iota1*rs_exact) > 2.0e-2_dp*2*iota1*rs_exact) &
                error stop 'shear'
            if (r%circularity > 1.0e-2_dp) error stop 'background not circular'
        end associate
        print '(a,3es12.4)', 'symmetry defect, min eig ratio, gauge defect ', &
            rep%symmetry_defect, rep%min_eigenvalue, rep%gauge_defect
        if (rep%symmetry_defect > 1.0e-13_dp) error stop 'G not symmetric'
        if (rep%min_eigenvalue <= 1.0e-12_dp) error stop 'G not positive definite'
        if (rep%gauge_defect > 1.0e-12_dp) error stop 'gauge transform'
        print '(a,2es14.6)', 'gradient ', rep%gradient
        print '(a,2es14.6)', 'amplitude ', rep%amplitude
        print '(a,2es24.15)', 'energy change model/direct ', rep%energy_model, &
            rep%energy_direct
        if (abs(rep%energy_direct - rep%energy_model) > 1.0e-12_dp*abs(rep%energy_model)) &
            error stop 'internal model energy'
        if (rep%energy_model >= 0) error stop 'seeding must lower the energy'
        ! Independent energy (relax_state: real Fourier blocks, harmonic
        ! included) at t = 0, 1, 2 times the seed: the quadratic model holds.
        call relax_prepare(field, rs, st)
        if (st /= 0) error stop 'relax prepare'
        delta = an(:, 1) - xs
        do t = 0, 2
            call energy_of(xs + t*delta, e(t))
        end do
        model = rep%energy_model
        defect = max(abs(e(1) - e(0) - model), &
            abs(e(2) - e(0) - 2*dot_product(rep%gradient, rep%amplitude) &
            - 2*dot_product(rep%amplitude, matmul(rep%gram, rep%amplitude))))
        print '(a,3es24.15)', 'relax_state W(shielded), W(seeded), W(2 seed) ', e
        print '(a,2es12.4)', 'model defect absolute, relative to W ', defect, &
            defect/abs(e(0))
        if (defect > 1.0e-14_dp*abs(e(0))) error stop 'quadratic model not exact'
        call energy_of(xv, eseed)
        print '(a,es24.15)', 'relax_state W(vacuum) ', eseed
        if (.not. (e(0) > e(1) .and. e(1) >= eseed)) error stop 'energy ordering'
        ! Resonant field: shielded ~ 0, seeded a fraction of the vacuum.
        call seed_resonant_br(vacuum, rep%res(1), opt%ntheta, dbr_vac, theta, st)
        if (st /= 0) error stop 'vacuum dbr'
        ratio = rep%res(1)%dbr_after/dbr_vac
        print '(a,4es14.6)', 'dbr shielded/seeded/vacuum, seeded/vacuum ', &
            rep%res(1)%dbr_before, rep%res(1)%dbr_after, dbr_vac, ratio
        if (rep%res(1)%dbr_before > 5.0e-2_dp*dbr_vac) error stop 'state not shielded'
        if (ratio <= 0 .or. ratio > 1) error stop 'seeded resonant field fraction'
        ! Pendulum width against traced islands of the seeded field.
        wtrace = traced_width(field, rep%res(1)%u, rep%res(1)%theta_o, &
            rep%res(1)%width_after)
        print '(a,3es14.6)', 'width pendulum, traced, ratio ', rep%res(1)%width_after, &
            wtrace, wtrace/rep%res(1)%width_after
        if (abs(wtrace/rep%res(1)%width_after - 1) > 0.05_dp) &
            error stop 'pendulum width vs traced width'
        ! eps/4: the same chain at a quarter of the resonant field.
        call fixture_scaled(0.25_dp, field, an4, st)
        if (st /= 0) error stop 'scaled fixture'
        call island_seed(field, opt, rep4, an4)
        if (rep4%status /= 0 .or. rep4%count /= 1) error stop 'scaled seed'
        print '(a,es24.15)', 'width(eps/4)/width(eps) ', &
            rep4%res(1)%width_after/rep%res(1)%width_after
        if (abs(rep4%res(1)%width_after/rep%res(1)%width_after - 0.5_dp) > 1.0e-9_dp) &
            error stop 'width does not scale as sqrt(dB_r)'
        wtrace4 = traced_width(field, rep4%res(1)%u, rep4%res(1)%theta_o, &
            rep4%res(1)%width_after)
        print '(a,3es14.6)', 'eps/4 width pendulum, traced, ratio ', &
            rep4%res(1)%width_after, wtrace4, wtrace4/rep4%res(1)%width_after
        if (abs(wtrace4/rep4%res(1)%width_after - 1) > 0.05_dp) &
            error stop 'pendulum width vs traced width (eps/4)'
        ! Hand-off 1: the NN1 Newton solve from the seeded potential returns
        ! the vacuum; its energy (no 4 pi, mode only) brackets the seed.
        call nn_solve(field, 1, 0.25_dp*xv, an4(:, 1), nrep, st, solution=sol)
        if (st /= 0) error stop 'Newton hand-off'
        print '(a,2es14.6,es12.4)', 'Newton energy before/after, |x - x_vac| ', &
            nrep%energy_before, nrep%energy_after, &
            maxval(abs(sol - 0.25_dp*xv))/maxval(abs(0.25_dp*xv))
        if (nrep%energy_after >= nrep%energy_before) error stop 'Newton energy'
        if (maxval(abs(sol - 0.25_dp*xv)) > 1.0e-9_dp*maxval(abs(xv))) &
            error stop 'Newton did not return to the vacuum'
        ! Hand-off 2: ideal (topology-preserving) relaxation from the seeded
        ! eps state; the energy must not increase and the step must succeed.
        call relax_pack_modes(rs, a0, an, a)
        call relax_state(rs, RELAX_IDEAL, a, st_s, st)
        if (st /= 0) error stop 'relax state seeded'
        ro%kind = RELAX_IDEAL; ro%dt = 1.0e-3_dp; ro%dt_max = 1.0e-1_dp
        ro%max_steps = 2; ro%handoff = 0
        call relax_run(rs, a, ro, rr)
        print '(a,i0,1x,a)', 'ideal relax status ', rr%status, trim(rr%message)
        print '(a,2es24.15)', 'ideal relax energy initial/final ', &
            rr%initial%energy, rr%final%energy
        if (rr%energy_increase > 1.0e-14_dp) error stop 'ideal relax energy increase'
        if (rr%final%energy > rr%initial%energy) error stop 'ideal relax energy'
        if (abs(rr%initial%energy - st_s%energy) > 1.0e-14_dp*abs(st_s%energy)) &
            error stop 'relax start is not the seeded state'
        call relax_free(rs)
    end subroutine

    subroutine energy_of(x1, w)
        complex(dp), intent(in) :: x1(:)
        real(dp), intent(out) :: w
        complex(dp), allocatable :: m(:, :)
        real(dp), allocatable :: packed(:)
        type(relax_state_t) :: s
        integer :: status
        allocate (m(size(x1), 1)); m(:, 1) = x1
        call relax_pack_modes(rs, a0, m, packed)
        call relax_state(rs, RELAX_IDEAL, packed, s, status)
        if (status /= 0) error stop 'relax_state'
        w = s%energy
    end subroutine

    subroutine fixture_scaled(f, fld, pot, status)
        real(dp), intent(in) :: f
        type(free_b_field_t), intent(inout) :: fld
        complex(dp), allocatable, intent(out) :: pot(:, :)
        integer, intent(out) :: status
        allocate (pot(size(xs), 1)); pot(:, 1) = f*xs
        call set_mode(fld, f*xs)
        status = 0
    end subroutine

    subroutine none_case()
        type(free_b_field_t) :: field, before
        type(seed_options_t) :: opt
        type(seed_report_t) :: rep
        complex(dp), allocatable :: xv(:), an(:, :)
        integer :: st
        iota0 = 0.55_dp; iota1 = 0.2_dp; eps = 1.0e-4_dp
        call fixture(4, field, a0, xv, xs, st, shield=.false.)
        if (st /= 0) error stop 'fixture'
        before = field
        allocate (an(field%fs%nn, 1)); an(:, 1) = xv
        opt%m_max = 4
        call island_seed(field, opt, rep, an)
        print '(a,i0,1x,a,i3)', 'status ', rep%status, trim(rep%message), rep%count
        if (rep%status /= 0 .or. rep%count /= 0) error stop 'unexpected resonance'
        if (any(field%b /= before%b) .or. any(an(:, 1) /= xv)) &
            error stop 'state changed without a resonance'
        ! Resonant background (q = 2 inside) with no n = 1 content: the
        ! gradient vanishes by Fourier orthogonality and nothing is seeded.
        iota0 = 0.4_dp
        call fixture(4, field, a0, xv, xs, st, shield=.false.)
        if (st /= 0) error stop 'fixture'
        field%b(:, 1) = 0
        call island_seed(field, opt, rep)
        print '(a,i3,2es12.4)', 'axisymmetric: count, |g|, |a| ', rep%count, &
            maxval(abs(rep%gradient)), maxval(abs(rep%amplitude))
        if (rep%status /= 0 .or. rep%count /= 1) error stop 'axisymmetric resonance'
        if (any(rep%amplitude /= 0)) error stop 'axisymmetric state was seeded'
        if (any(field%b(:, 1) /= 0)) error stop 'axisymmetric state changed'
    end subroutine

    subroutine ladder_case()
        type(free_b_field_t) :: field
        type(seed_options_t) :: opt
        type(seed_report_t) :: rep
        complex(dp), allocatable :: xv(:), an(:, :)
        real(dp) :: w(3), u(3), dbr(3), dv, th
        type(free_b_field_t) :: vac
        integer :: st, level
        iota0 = 0.4_dp; iota1 = 0.1_dp/0.5625_dp; eps = 1.0e-4_dp
        opt%m_max = 4
        do level = 1, 3
            call fixture(2**(level + 2), field, a0, xv, xs, st)
            if (st /= 0) error stop 'fixture'
            allocate (an(field%fs%nn, 1)); an(:, 1) = xs
            vac = field
            call set_mode(vac, xv)
            call island_seed(field, opt, rep, an)
            if (rep%status /= 0 .or. rep%count /= 1) error stop 'ladder seed'
            call seed_resonant_br(vac, rep%res(1), opt%ntheta, dv, th, st)
            if (st /= 0) error stop 'vacuum dbr'
            w(level) = rep%res(1)%width_after; u(level) = rep%res(1)%u
            dbr(level) = rep%res(1)%dbr_after
            print '(a,i3,6es14.6)', 'rings, u_res, shear, G, dbr, seeded/vacuum, width ', &
                2**(level + 2), u(level), rep%res(1)%shear, rep%gram(1, 1), dbr(level), &
                dbr(level)/dv, w(level)
            deallocate (an)
        end do
        print '(a,2es12.4,f8.3)', 'width differences and observed order ', &
            abs(w(2) - w(1)), abs(w(3) - w(2)), &
            log(abs(w(2) - w(1))/abs(w(3) - w(2)))/log(2.0_dp)
        if (abs(w(3) - w(2)) >= abs(w(2) - w(1))) error stop 'width not converging'
        if (abs(w(3) - w(2)) > 0.01_dp*w(3)) error stop 'finest width change > 1%'
    end subroutine

    subroutine leaves_case()
        ! Generated leaves: identities, JVP against central differences and
        ! the JVP/VJP pairing <y_bar, J x_dot> = <J^T y_bar, x_dot>.
        real(dp) :: x7(7), d7(7), g7(7), y(2), yp(2), ym(2), dy(2), yb(2)
        real(dp) :: x4(4), d4(4), g4(4), z(3), zp(3), zm(3), dz(3), zb(3)
        real(dp) :: x5(5), d5(5), g5(5), w, wp, wm, dw, h, err
        h = 1.0e-6_dp
        ! Pure poloidal field B_theta = 1 on the circle: rotation = R/(r B_phi).
        call circle_rotation(3.0_dp, 0.5_dp, cos(0.3_dp), sin(0.3_dp), -sin(0.3_dp), &
            cos(0.3_dp), 2.0_dp, y(1), y(2))
        if (abs(y(1) - 3.0_dp) > 1.0e-15_dp .or. abs(y(2)) > 1.0e-15_dp) &
            error stop 'circle_rotation identity'
        x7 = [3.1_dp, 0.4_dp, 0.8_dp, 0.6_dp, 0.2_dp, -0.3_dp, 1.3_dp]
        d7 = [0.3_dp, -0.2_dp, 0.1_dp, 0.5_dp, -0.4_dp, 0.7_dp, 0.2_dp]
        yb = [0.7_dp, -1.1_dp]
        call circle_rotation(x7(1) + h*d7(1), x7(2) + h*d7(2), x7(3) + h*d7(3), &
            x7(4) + h*d7(4), x7(5) + h*d7(5), x7(6) + h*d7(6), x7(7) + h*d7(7), &
            yp(1), yp(2))
        call circle_rotation(x7(1) - h*d7(1), x7(2) - h*d7(2), x7(3) - h*d7(3), &
            x7(4) - h*d7(4), x7(5) - h*d7(5), x7(6) - h*d7(6), x7(7) - h*d7(7), &
            ym(1), ym(2))
        call circle_rotation_jvp(x7(1), x7(2), x7(3), x7(4), x7(5), x7(6), x7(7), &
            d7(1), d7(2), d7(3), d7(4), d7(5), d7(6), d7(7), dy(1), dy(2))
        call circle_rotation_vjp(x7(1), x7(2), x7(3), x7(4), x7(5), x7(6), x7(7), &
            yb(1), yb(2), g7(1), g7(2), g7(3), g7(4), g7(5), g7(6), g7(7))
        err = maxval(abs(dy - (yp - ym)/(2*h)))
        print '(a,2es12.4)', 'circle_rotation JVP-FD, pairing ', err, &
            abs(dot_product(yb, dy) - dot_product(g7, d7))
        if (err > 1.0e-8_dp) error stop 'circle_rotation JVP'
        if (abs(dot_product(yb, dy) - dot_product(g7, d7)) > 1.0e-13_dp) &
            error stop 'circle_rotation VJP'
        ! Seed potential: unit direction times (1 - xi^2)^2; zero at |xi| = 1.
        call seed_potential(0.3_dp, 2.0_dp, -0.4_dp, 0.0_dp, z(1), z(2), z(3))
        if (abs(norm2(z) - 1) > 1.0e-15_dp) error stop 'seed_potential direction'
        call seed_potential(0.3_dp, 2.0_dp, -0.4_dp, 1.0_dp, z(1), z(2), z(3))
        if (any(z /= 0)) error stop 'seed_potential support edge'
        x4 = [0.3_dp, 2.0_dp, -0.4_dp, 0.35_dp]; d4 = [0.2_dp, -0.1_dp, 0.3_dp, 0.4_dp]
        zb = [0.5_dp, -0.2_dp, 0.9_dp]
        call seed_potential(x4(1) + h*d4(1), x4(2) + h*d4(2), x4(3) + h*d4(3), &
            x4(4) + h*d4(4), zp(1), zp(2), zp(3))
        call seed_potential(x4(1) - h*d4(1), x4(2) - h*d4(2), x4(3) - h*d4(3), &
            x4(4) - h*d4(4), zm(1), zm(2), zm(3))
        call seed_potential_jvp(x4(1), x4(2), x4(3), x4(4), d4(1), d4(2), d4(3), &
            d4(4), dz(1), dz(2), dz(3))
        call seed_potential_vjp(x4(1), x4(2), x4(3), x4(4), zb(1), zb(2), zb(3), &
            g4(1), g4(2), g4(3), g4(4))
        err = maxval(abs(dz - (zp - zm)/(2*h)))
        print '(a,2es12.4)', 'seed_potential JVP-FD, pairing ', err, &
            abs(dot_product(zb, dz) - dot_product(g4, d4))
        if (err > 1.0e-8_dp) error stop 'seed_potential JVP'
        if (abs(dot_product(zb, dz) - dot_product(g4, d4)) > 1.0e-13_dp) &
            error stop 'seed_potential VJP'
        ! Pendulum width: 4 sqrt(R dbr/(m iota' B_phi)).
        x5 = [3.0_dp, 1.0e-4_dp, 2.0_dp, 0.28_dp, 1.0_dp]
        d5 = [0.1_dp, 2.0e-5_dp, 0.0_dp, -0.03_dp, 0.05_dp]
        call pendulum_width(x5(1), x5(2), x5(3), x5(4), x5(5), w)
        if (abs(w - 4*sqrt(3.0e-4_dp/0.56_dp)) > 1.0e-15_dp) error stop 'pendulum_width'
        call pendulum_width(x5(1) + h*d5(1), x5(2) + h*d5(2), x5(3), x5(4) + h*d5(4), &
            x5(5) + h*d5(5), wp)
        call pendulum_width(x5(1) - h*d5(1), x5(2) - h*d5(2), x5(3), x5(4) - h*d5(4), &
            x5(5) - h*d5(5), wm)
        call pendulum_width_jvp(x5(1), x5(2), x5(3), x5(4), x5(5), d5(1), d5(2), &
            d5(3), d5(4), d5(5), dw)
        call pendulum_width_vjp(x5(1), x5(2), x5(3), x5(4), x5(5), 1.3_dp, &
            g5(1), g5(2), g5(3), g5(4), g5(5))
        err = abs(dw - (wp - wm)/(2*h))
        print '(a,2es12.4)', 'pendulum_width JVP-FD, pairing ', err, &
            abs(1.3_dp*dw - dot_product(g5, d5))
        if (err > 1.0e-8_dp*abs(dw)) error stop 'pendulum_width JVP'
        if (abs(1.3_dp*dw - dot_product(g5, d5)) > 1.0e-14_dp) &
            error stop 'pendulum_width VJP'
    end subroutine

    subroutine fixture(nr, field, a0, xv, xs, status, shield)
        integer, intent(in) :: nr
        type(free_b_field_t), intent(out) :: field
        complex(dp), allocatable, intent(out) :: a0(:), xv(:), xs(:)
        integer, intent(out) :: status
        logical, intent(in), optional :: shield
        type(flux_mesh_t) :: mesh
        type(csc_z_t) :: curl, zk
        type(csc_t) :: k
        real(dp), allocatable :: rings(:), nodes(:, :)
        complex(dp), allocatable :: full(:), ca(:), y(:), rhs(:)
        real(dp) :: r
        integer :: i
        logical :: do_shield
        do_shield = .true.
        if (present(shield)) do_shield = shield
        allocate (rings(nr)); rings = [(real(i, dp)/nr, i=1, nr)]
        call flux_mesh_build(rings, 8, 32*nr, 2, mapping, mesh, status)
        if (status /= 0) return
        call nn_build(mesh, 2, 1, [r0, 0.0_dp], 1.0_dp, r0*b0, 1.0_dp, field, &
            status, cylinder_r0=r0)
        if (status /= 0) return
        ! Background a = R0 A_z with B_theta = a'/R0 and iota = R0 B_theta/(r B_z).
        call meridian_lagrange_nodes(field%fs, nodes)
        allocate (a0(field%fs%n0))
        do i = 1, field%fs%n0
            r = norm2(nodes(:, i) - [r0, 0.0_dp])
            a0(i) = b0*(iota0*r**2/2 + iota1*r**4/4)
        end do
        allocate (full(field%fs%sp(1)%ndof)); full = 0
        full(field%fs%nn + 1:) = a0
        call meridian_op_complex(field%fs%curl, 0, curl, status)
        if (status /= 0) return
        field%b(:, 0) = csc_matvec(curl, full)
        field%psi = real(a0, dp)
        call vacuum_mode(field, xv, status)
        if (status /= 0) return
        xs = xv
        if (.not. do_shield) return
        ! Shielded minimiser: x_s = x_v - y (c.x_v)/(c.y), y = K^{-1} conj(c),
        ! on the analytic q = 2 circle r = 3/4, a mesh ring: the shielding
        ! sheet lies on element edges, where the normal B_r is continuous.
        call resonant_row(field, 0.75_dp, ca, status)
        if (status /= 0) return
        call meridian_stiffness(field%fs, 1, n_res, zk, status)
        if (status /= 0) return
        call leading_real(zk, field%fs%nn, k, status)
        if (status /= 0) return
        rhs = conjg(ca)
        where (field%fs%sp(1)%bnd(:field%fs%nn)) rhs = 0
        allocate (y(field%fs%nn))
        call nn_linear_solve(k, field%fs%sp(1)%bnd(:field%fs%nn), rhs, y, status)
        if (status /= 0) return
        xs = xv - y*sum(ca*xv)/sum(ca*y)
        call set_mode(field, xs)
    end subroutine

    subroutine set_mode(field, x)
        type(free_b_field_t), intent(inout) :: field
        complex(dp), intent(in) :: x(:)
        type(csc_z_t) :: curl
        complex(dp), allocatable :: full(:)
        integer :: status
        call meridian_op_complex(field%fs%curl, n_res, curl, status)
        if (status /= 0) error stop 'curl'
        allocate (full(field%fs%sp(1)%ndof)); full = 0
        full(:field%fs%nn) = x
        field%b(:, n_res) = csc_matvec(curl, full)
    end subroutine

    subroutine vacuum_mode(field, xv, status)
        ! nn_vacuum_mode with the boundary potential of the trace below,
        ! projected by an explicit quadrature loop (no procedure argument).
        type(free_b_field_t), intent(inout) :: field
        complex(dp), allocatable, intent(out) :: xv(:)
        integer, intent(out) :: status
        type(nn_report_t) :: rep
        type(csc_z_t) :: zm
        type(csc_t) :: mass
        complex(dp), allocatable :: load(:), a(:), seed(:)
        real(dp), allocatable :: phi(:, :)
        logical, allocatable :: free(:)
        real(dp) :: x(2), wq
        complex(dp) :: v(3)
        integer :: t, q, i, nn
        nn = field%fs%nn
        allocate (load(field%fs%sp(1)%ndof), phi(3, field%fs%sp(1)%nloc))
        load = 0
        do t = 1, field%fs%mm%topology%n_triangles
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 1, t, q, x, wq, phi, status)
                if (status /= 0) return
                call trace_potential(x, v)
                do i = 1, size(phi, 2)
                    load(field%fs%sp(1)%dof(i, t)) = load(field%fs%sp(1)%dof(i, t)) &
                        + field%fs%sp(1)%sgn(i, t)*wq*sum(v*phi(:, i))
                end do
            end do
        end do
        call meridian_real_to_complex(field%fs%mass(1), zm)
        call leading_real(zm, nn, mass, status)
        if (status /= 0) return
        allocate (a(nn), seed(nn), free(nn)); free = .false.; seed = 0
        call nn_linear_solve(mass, free, load(:nn), a, status)
        if (status /= 0) return
        call nn_solve(field, n_res, a, seed, rep, status, solution=xv)
    end subroutine

    pure subroutine trace_potential(x, v)
        ! Gauge potential (A_R, A_phi, A_Z) = (i R0 B_Z/n, 0, -i R0 B_R/n) of
        ! B = grad psi x e_z, psi = (eps/2) ((x - i y)/a)^m (cos(m theta - phi)).
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        complex(dp) :: z, px, py, br, bz
        z = cmplx(x(1) - r0, -x(2), dp)
        px = eps/2*m_res*z**(m_res - 1)
        py = eps/2*cmplx(0.0_dp, -m_res, dp)*z**(m_res - 1)
        br = py; bz = -px
        v = [cmplx(0.0_dp, r0, dp)*bz/n_res, (0.0_dp, 0.0_dp), &
            -cmplx(0.0_dp, r0, dp)*br/n_res]
    end subroutine

    subroutine resonant_row(field, radius, ca, status)
        ! c(x) = <e^{i m theta} r_hat . B_n(x)>_theta on the circle (the
        ! element-split rule of the module) as a row on the Nedelec
        ! potential: ca = curl_n^T cb.
        type(free_b_field_t), intent(in) :: field
        real(dp), intent(in) :: radius
        complex(dp), allocatable, intent(out) :: ca(:)
        integer, intent(out) :: status
        type(csc_z_t) :: curl
        complex(dp), allocatable :: cb(:)
        real(dp), allocatable :: theta(:), weight(:)
        real(dp) :: x(2), u, th, xi, eta, y(2), v2(3, field%fs%sp(2)%nloc)
        integer :: j, t, l, p
        call seed_circle_quadrature(field, radius, 256, theta, weight, status)
        if (status /= 0) return
        allocate (cb(field%fs%sp(2)%ndof)); cb = 0
        do j = 1, size(theta)
            x = [r0, 0.0_dp] + radius*[cos(theta(j)), sin(theta(j))]
            call field%coordinates(field%centre, field%radial_scale, x, u, th)
            call flux_mesh_locate(field%mesh, x, u, th, t, xi, eta, status)
            if (status /= 0) return
            call meridian_point_basis(field%fs, field%pb, 2, t, xi, eta, y, v2, status)
            if (status /= 0) return
            do l = 1, size(v2, 2)
                cb(field%fs%sp(2)%dof(l, t)) = cb(field%fs%sp(2)%dof(l, t)) &
                    + weight(j)*exp(cmplx(0.0_dp, m_res*theta(j), dp)) &
                    *(cos(theta(j))*v2(1, l) + sin(theta(j))*v2(3, l)) &
                    *field%fs%sp(2)%sgn(l, t)
            end do
        end do
        call meridian_op_complex(field%fs%curl, n_res, curl, status)
        if (status /= 0) return
        allocate (ca(field%fs%nn)); ca = 0
        do j = 1, field%fs%nn
            do p = curl%col_ptr(j), curl%col_ptr(j + 1) - 1
                ca(j) = ca(j) + cb(curl%row_idx(p))*curl%val(p)
            end do
        end do
        status = 0
    end subroutine

    subroutine leading_real(z, n, a, status)
        type(csc_z_t), intent(in) :: z
        integer, intent(in) :: n
        type(csc_t), intent(out) :: a
        integer, intent(out) :: status
        type(fortsparse_status_t) :: fst
        integer, allocatable :: row(:), col(:)
        real(dp), allocatable :: val(:)
        integer :: j, p, c
        allocate (row(z%nnz), col(z%nnz), val(z%nnz))
        c = 0
        do j = 1, n
            do p = z%col_ptr(j), z%col_ptr(j + 1) - 1
                if (z%row_idx(p) > n) cycle
                c = c + 1
                row(c) = z%row_idx(p); col(c) = j; val(c) = real(z%val(p), dp)
            end do
        end do
        call csc_from_triplet(n, n, row(:c), col(:c), val(:c), a, fst)
        status = merge(0, 1, status_ok(fst))
    end subroutine

    real(dp) function traced_width(field, u_res, theta_o, w_guess) result(width)
        ! Radial extent on the phi = 0 O-point ray of the field lines whose
        ! helical angle chi = m theta - n phi librates (excursion < 1.5 pi)
        ! over 160 toroidal transits; edges by bisection.
        type(free_b_field_t), intent(in) :: field
        real(dp), intent(in) :: u_res, theta_o, w_guess
        real(dp) :: lower, upper
        !$omp parallel sections default(shared)
        !$omp section
        lower = trace_edge(field, u_res, theta_o, u_res - 1.5_dp*w_guess)
        !$omp section
        upper = trace_edge(field, u_res, theta_o, u_res + 1.5_dp*w_guess)
        !$omp end parallel sections
        width = upper - lower
    end function

    real(dp) function trace_edge(field, u_res, theta_o, outside) result(r)
        type(free_b_field_t), intent(in) :: field
        real(dp), intent(in) :: u_res, theta_o, outside
        real(dp) :: a, b, mid
        integer :: it
        a = u_res; b = outside
        if (.not. librates(field, theta_o, a)) error stop 'O point does not librate'
        if (librates(field, theta_o, b)) error stop 'bracket outside the island librates'
        do it = 1, 16
            mid = (a + b)/2
            if (librates(field, theta_o, mid)) then
                a = mid
            else
                b = mid
            end if
        end do
        r = (a + b)/2
    end function

    logical function librates(field, theta_o, u)
        type(free_b_field_t), intent(in) :: field
        real(dp), intent(in) :: theta_o, u
        integer, parameter :: turns = 160, steps = 32
        real(dp) :: p(2), k1(2), k2(2), k3(2), k4(2), phi, h, chi0, chi, &
            theta, theta_prev, unwrap
        integer :: i
        h = 2*pi/steps
        p = [r0, 0.0_dp] + u*[cos(theta_o), sin(theta_o)]
        phi = 0
        theta_prev = theta_o; unwrap = theta_o
        chi0 = m_res*theta_o
        librates = .true.
        do i = 1, turns*steps
            call velocity(field, p, phi, k1)
            call velocity(field, p + h/2*k1, phi + h/2, k2)
            call velocity(field, p + h/2*k2, phi + h/2, k3)
            call velocity(field, p + h*k3, phi + h, k4)
            p = p + h/6*(k1 + 2*k2 + 2*k3 + k4)
            phi = phi + h
            theta = atan2(p(2), p(1) - r0)
            unwrap = unwrap + modulo(theta - theta_prev + pi, 2*pi) - pi
            theta_prev = theta
            chi = m_res*unwrap - n_res*phi
            if (abs(chi - chi0) > 1.5_dp*pi) then
                librates = .false.
                return
            end if
        end do
    end function

    subroutine velocity(field, p, phi, v)
        type(free_b_field_t), intent(in) :: field
        real(dp), intent(in) :: p(2), phi
        real(dp), intent(out) :: v(2)
        integer :: status
        call nn_velocity(field, [p, phi], v, status)
        if (status /= 0) error stop 'field line left the domain'
    end subroutine

    subroutine mapping(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [r0 + u*cos(theta), u*sin(theta)]
    end subroutine
end program test_island_seed
