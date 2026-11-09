program test_screening_moments
    ! One executable, one check per CTest name (B05-B07 first-order moment solver).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_screening_moments, only: scr_t, scr_setup, scr_stencil, scr_steady, &
        scr_rate0, scr_transient
    use milca_fp_relaxation, only: fp_relax_t, fp_relax_setup, fp_relax_cells
    use kin6d_fp_relaxation_leaves, only: fp_rel_maxwell, fp_rel_friction, fp_rel_diff_par, &
        fp_rel_diff_perp, fp_rel_cross, fp_rel_face_lo, fp_rel_face_hi, fp_rel_speed
    use kin6d_thermal_screening_leaves, only: scr_friction_dt, scr_diff_par_dt, &
        scr_diff_perp_dt, scr_cross_dt, scr_distortion_k
    use milca_rosenbluth_reference, only: homma_pair_t, homma_radial_quad, &
        homma_force_maxwellian, homma_closure, homma_eps0, homma_vt
    implicit none
    ! B05-B07 parameters (thermal-force-benchmark/cases, conventions.md).
    real(dp), parameter :: amu = 1.66053906660e-27_dp, e_charge = 1.602176634e-19_dp
    real(dp), parameter :: ma = 183.84_dp*amu, mb = 1.67262192369e-27_dp, qb = e_charge
    real(dp), parameter :: nb = 1.0e20_dp, t0 = 50.0_dp*e_charge, lnl = 15.0_dp
    real(dp), parameter :: tail = 1.0e-14_dp, ppu = 4.0_dp, grad100 = 100.0_dp*e_charge
    character(len=32) :: name

    call get_command_argument(1, name)
    select case (trim(name))
    case ('maxwellian')
        call check_maxwellian()
    case ('zero_density')
        call check_zero_density()
    case ('rate0')
        call check_rate0()
    case ('strong_field')
        call check_strong_field()
    case ('flr_ratio')
        call check_flr_ratio()
    case ('convergence_v')
        call check_convergence_v()
    case ('oracle_rosenbluth')
        call check_oracle_rosenbluth()
    case ('oracle_poisson')
        call check_oracle_poisson()
    case ('temperature_derivative')
        call check_temperature_derivative()
    case ('transient')
        call check_transient()
    case ('steady_transient')
        call check_steady_transient()
    case default
        error stop 'usage: test_screening_moments maxwellian|zero_density|rate0|strong_field|'// &
            'flr_ratio|convergence_v|oracle_rosenbluth|oracle_poisson|temperature_derivative|'// &
            'transient|steady_transient'
    end select
    print '(a)', 'PASS '//trim(name)

contains

    subroutine setup(s, z, cpu, omega)
        type(scr_t), intent(out) :: s
        integer, intent(in) :: z, cpu
        real(dp), intent(out) :: omega
        call scr_setup(s, fp_relax_cells(cpu, ma, t0, 0.0_dp, t0, tail), ma, mb, z*e_charge, &
            qb, nb, t0, lnl, tail, ppu)
        omega = z*e_charge/ma
    end subroutine setup

    function pair_of(z) result(p)
        integer, intent(in) :: z
        type(homma_pair_t) :: p
        p%n_b = nb; p%t_b = t0; p%m_b = mb; p%q_b = qb
        p%m_a = ma; p%q_a = z*e_charge; p%lnlambda = lnl
    end function pair_of

    function closure(grad, bmag) result(q)
        real(dp), intent(in) :: grad(3), bmag
        real(dp) :: q(3)
        call homma_closure(pair_of(1), [0.0_dp, 0.0_dp, 1.0_dp], bmag, grad, q)
    end function closure

    ! B04's l = 0 stencil (and its T_b derivative) on the same grid as the solver.
    subroutine l0_stencils(cpu, kl, kd, ku, klt, kdt, kut, f)
        integer, intent(in) :: cpu
        real(dp), allocatable, intent(out) :: kl(:), kd(:), ku(:), klt(:), kdt(:), kut(:)
        type(fp_relax_t), intent(out) :: f
        real(dp), allocatable :: rf(:), at(:), et(:), dpt(:), dperpt(:), pft(:), qft(:)
        integer :: n

        call fp_relax_setup(f, fp_relax_cells(cpu, ma, t0, 0.0_dp, t0, tail), ma, mb, &
            3.0_dp*e_charge, qb, nb, t0, lnl, homma_eps0, t0, t0, 0.0_dp, 1.0_dp, tail)
        n = f%ncell
        call scr_stencil(f%pf, f%qf, f%h, f%dperp, 0.0_dp, kl, kd, ku)
        allocate (rf(n - 1), at(n - 1), et(n - 1), dpt(n), dperpt(n), pft(0:n), qft(0:n))
        rf = f%vface(1:n - 1)
        call scr_friction_dt(rf, t0, mb, nb, f%gam, f%mr, at)
        call scr_cross_dt(rf, t0, mb, nb, f%gam, et)
        call scr_diff_par_dt(f%vc, t0, mb, nb, f%gam, dpt)
        call scr_diff_perp_dt(f%vc, t0, mb, nb, f%gam, dperpt)
        pft = 0.0_dp; qft = 0.0_dp
        call fp_rel_face_lo(rf, at, et, dpt(1:n - 1), f%h, pft(1:n - 1))
        call fp_rel_face_hi(rf, at, et, dpt(2:n), f%h, qft(1:n - 1))
        call scr_stencil(pft, qft, f%h, dperpt, 0.0_dp, klt, kdt, kut)
    end subroutine l0_stencils

    subroutine check_maxwellian()
        ! L f_M = C_T0[f_M] (the Lorentz term annihilates isotropic f, proved in the
        ! generator). The discrete l = 0 residual of f_M at the cell centres, relative
        ! to the operator scale, must fall at second order.
        type(fp_relax_t) :: f
        real(dp), allocatable :: kl(:), kd(:), ku(:), klt(:), kdt(:), kut(:), fm(:), r(:)
        real(dp) :: res(3)
        integer :: k, n, cpu_list(3)

        cpu_list = [16, 32, 64]
        do k = 1, 3
            call l0_stencils(cpu_list(k), kl, kd, ku, klt, kdt, kut, f)
            n = f%ncell
            allocate (fm(n), r(n))
            call fp_rel_maxwell(f%vc, ma, t0, fm)
            r = kd*fm
            r(2:n) = r(2:n) + kl(2:n)*fm(1:n - 1)
            r(1:n - 1) = r(1:n - 1) + ku(1:n - 1)*fm(2:n)
            res(k) = sum(abs(r))/sum(abs(kd*fm))
            print '(a,i4,a,es12.4)', 'L f_M residual cells_per_unit=', cpu_list(k), ' ', res(k)
            deallocate (fm, r)
        end do
        if (res(2) > res(1)/3.0_dp .or. res(3) > res(2)/3.0_dp) error stop 'L f_M not O(h^2)'
        if (res(3) > 1.0e-3_dp) error stop 'L f_M residual above 1e-3 at 64 cells'
    end subroutine check_maxwellian

    subroutine check_zero_density()
        ! n1 is pure l = 1, so its density vanishes by construction; what makes this
        ! hold for the discretised N as well is that K and d_T K conserve density:
        ! every column of the l = 0 stencils sums to zero.
        type(fp_relax_t) :: f
        real(dp), allocatable :: kl(:), kd(:), ku(:), klt(:), kdt(:), kut(:)
        real(dp) :: worst, worst_t
        integer :: j, n

        call l0_stencils(32, kl, kd, ku, klt, kdt, kut, f)
        n = f%ncell
        worst = 0.0_dp; worst_t = 0.0_dp
        do j = 1, n
            worst = max(worst, abs(colsum(kl, kd, ku, j, n))/maxval(abs(kd)))
            worst_t = max(worst_t, abs(colsum(klt, kdt, kut, j, n))/maxval(abs(kdt)))
        end do
        print '(a,2es12.4)', 'largest column sum of K, d_T K (relative):', worst, worst_t
        if (worst > 1.0e-13_dp .or. worst_t > 1.0e-13_dp) error stop 'stencil changes density'
    end subroutine check_zero_density

    real(dp) function colsum(kl, kd, ku, j, n)
        real(dp), intent(in) :: kl(:), kd(:), ku(:)
        integer, intent(in) :: j, n
        colsum = kd(j)
        if (j < n) colsum = colsum + kl(j + 1)
        if (j > 1) colsum = colsum + ku(j - 1)
    end function colsum

    subroutine check_rate0()
        ! d<v>/dt(0) = int v C^q[f_M] = <A^T> over the test Maxwellian: the B03 oracle
        ! average of grad H_T, independent of the l = 1 source P.
        type(scr_t) :: s
        type(homma_pair_t) :: p
        real(dp) :: omega, q(3), rate(3), force(3), ferr, rel(2)
        integer :: k, cpu_list(2)

        q = [1.0e3_dp, -2.0e3_dp, 3.0e3_dp]
        p = pair_of(3)
        p%heat_flux = q
        call homma_force_maxwellian(p, t0, 8.0_dp, 8.0_dp, force, ferr)
        cpu_list = [16, 32]
        do k = 1, 2
            call setup(s, 3, cpu_list(k), omega)
            call scr_rate0(s, omega, [0.0_dp, 0.0_dp, 0.0_dp], q, .false., rate)
            rel(k) = maxval(abs(rate - force/ma))/maxval(abs(force/ma))
            print '(a,i3,a,3es15.7,a,es10.3)', 'cpu=', cpu_list(k), ' rate0 ', rate, ' rel ', rel(k)
        end do
        print '(a,3es15.7)', 'B03 oracle      ', force/ma
        if (rel(2) > 1.0e-3_dp) error stop 't = 0 rate differs from the B03 oracle'
        if (rel(2) > rel(1)/3.0_dp) error stop 't = 0 rate not refining at second order'
    end subroutine check_rate0

    subroutine check_strong_field()
        ! At fixed q and grad T = 0, Omega v_X -> <A^T_y>, the Maxwellian average that
        ! gives gc_thermal = <A^T_y>/Omega. On one grid the limit is the semi-discrete
        ! t = 0 rate, approached as (nu/Omega)^2; the rate itself is the B03 oracle
        ! up to the radial error.
        type(scr_t) :: s
        type(homma_pair_t) :: p
        real(dp) :: omega, q(3), drift(3), rate(3), force(3), ferr, dev(3), scale
        integer :: k

        q = [0.0_dp, 1.0e4_dp, 0.0_dp]
        p = pair_of(3)
        p%heat_flux = q
        call homma_force_maxwellian(p, t0, 8.0_dp, 8.0_dp, force, ferr)
        call setup(s, 3, 32, omega)
        call scr_rate0(s, omega, [0.0_dp, 0.0_dp, 0.0_dp], q, .false., rate)
        do k = 1, 3
            scale = 4.0_dp**(k - 1)
            call scr_steady(s, scale*omega, [0.0_dp, 0.0_dp, 0.0_dp], q, .false., drift)
            dev(k) = drift(1)*scale*omega/rate(2) - 1.0_dp
            print '(a,f5.1,a,es14.6,a,es11.3,a,es11.3)', 'B/B0=', scale, ' v_X=', drift(1), &
                ' Omega v_X/rate_y - 1=', dev(k), ' Omega v_X/<A_y> - 1=', &
                drift(1)*scale*omega/(force(2)/ma) - 1.0_dp
        end do
        if (abs(dev(2)) > abs(dev(1))/8.0_dp .or. abs(dev(3)) > abs(dev(2))/8.0_dp) &
            error stop 'strong-field limit not approached as (nu/Omega)^2'
        if (abs(drift(1)*scale*omega/(force(2)/ma) - 1.0_dp) > 1.0e-3_dp) &
            error stop 'strong-field drift misses gc_thermal'
    end subroutine check_strong_field

    subroutine check_flr_ratio()
        ! With the homma2014 closure, (full - gc)/gc -> (q_b/q_a)(T_a/T_b) = 1/Z as nu/Omega -> 0.
        type(scr_t) :: s
        real(dp) :: omega, g(3), full(3), gcv(3), ratio(2), bmag
        integer :: z, k

        g = [grad100, 0.0_dp, 0.0_dp]
        do z = 1, 3, 2
            call setup(s, z, 32, omega)
            do k = 1, 2
                bmag = 10.0_dp**(k - 1)
                call scr_steady(s, bmag*omega, g, closure(g, bmag), .false., full)
                call scr_steady(s, bmag*omega, g, closure(g, bmag), .true., gcv)
                ratio(k) = (full(1) - gcv(1))/gcv(1)
                print '(a,i0,a,f5.1,a,es14.6,a,es14.6)', 'Z=', z, ' B=', bmag, &
                    ' FLR/gc=', ratio(k), ' Z*FLR/gc-1=', z*ratio(k) - 1.0_dp
            end do
            if (abs(z*ratio(2) - 1.0_dp) > 1.0e-3_dp) error stop 'FLR ratio misses q_b T_a/(q_a T_b)'
            if (abs(z*ratio(2) - 1.0_dp) > abs(z*ratio(1) - 1.0_dp)) &
                error stop 'FLR ratio not approached'
        end do
    end subroutine check_flr_ratio

    subroutine check_convergence_v()
        ! Full-variant B07 drift (W3+, 100 eV/m) under radial refinement: second order.
        type(scr_t) :: s
        real(dp) :: omega, g(3), d(3), vx(3), order
        integer :: k, cpu_list(3)

        cpu_list = [8, 16, 32]
        g = [grad100, 0.0_dp, 0.0_dp]
        do k = 1, 3
            call setup(s, 3, cpu_list(k), omega)
            call scr_steady(s, omega, g, closure(g, 1.0_dp), .false., d)
            vx(k) = d(1)
            print '(a,i3,a,es22.14)', 'cpu=', cpu_list(k), ' v_X=', vx(k)
        end do
        order = log(abs(vx(1) - vx(2))/abs(vx(2) - vx(3)))/log(2.0_dp)
        print '(a,f6.3)', 'observed order ', order
        if (order < 1.8_dp .or. order > 2.6_dp) error stop 'radial convergence is not second order'
    end subroutine check_convergence_v

    subroutine check_oracle_rosenbluth()
        ! lap G^T = 2 H^T for the l = 1 jets: 2 h1 = g1'' + 2 g1'/v - 2 g1/v^2, for
        ! 0.02 <= v/v_t,b <= 1.6 (test ions have v/v_t,b < 0.5).
        type(homma_pair_t) :: p
        real(dp) :: rv(10), lhs, rhs, vt, worst
        integer :: k

        p = pair_of(3)
        vt = homma_vt(p)
        worst = 0.0_dp
        do k = 1, 5
            call homma_radial_quad(p, 0.02_dp*vt*3.0_dp**(k - 1), 8.0_dp, 8.0_dp, rv)
            lhs = 2.0_dp*rv(6)
            rhs = rv(10) + 2.0_dp*rv(9)/(0.02_dp*vt*3.0_dp**(k - 1)) - &
                2.0_dp*rv(8)/(0.02_dp*vt*3.0_dp**(k - 1))**2
            worst = max(worst, abs(lhs - rhs)/abs(lhs))
        end do
        print '(a,es12.4)', 'largest relative Rosenbluth-identity residual ', worst
        if (worst > 1.0e-10_dp) error stop 'oracle violates lap G = 2 H'
    end subroutine check_oracle_rosenbluth

    subroutine check_oracle_poisson()
        ! lap H^T = -4 pi K s: (v^2 h1')' - 2 h1 = -4 pi v^2 K, the derivative by a
        ! central difference of the oracle's v^2 h1', for 0.02 <= v/v_t,b <= 1.6.
        type(homma_pair_t) :: p
        real(dp) :: rp(10), rm(10), r0(10), v, dv, lhs, kq, vt, worst
        integer :: k

        p = pair_of(3)
        vt = homma_vt(p)
        worst = 0.0_dp
        do k = 1, 5
            v = 0.02_dp*vt*3.0_dp**(k - 1)
            dv = 1.0e-4_dp*v
            call homma_radial_quad(p, v + dv, 8.0_dp, 8.0_dp, rp)
            call homma_radial_quad(p, v - dv, 8.0_dp, 8.0_dp, rm)
            call homma_radial_quad(p, v, 8.0_dp, 8.0_dp, r0)
            lhs = ((v + dv)**2*rp(7) - (v - dv)**2*rm(7))/(2.0_dp*dv) - 2.0_dp*r0(6)
            call scr_distortion_k(v, t0, mb, kq)
            worst = max(worst, abs(lhs + 4.0_dp*acos(-1.0_dp)*v**2*kq)/abs(4.0_dp*acos(-1.0_dp)*v**2*kq))
        end do
        print '(a,es12.4)', 'largest relative Poisson residual ', worst
        if (worst > 1.0e-6_dp) error stop 'oracle h1 and K violate lap H = -4 pi f_b'
    end subroutine check_oracle_poisson

    subroutine check_temperature_derivative()
        ! The d_T leaves against a central difference in T_b of the B04 leaves, from
        ! 0.05 v_t,b up (below, cancellation in the B04 leaves swamps the quotient).
        type(homma_pair_t) :: p
        real(dp) :: v, vp, vm, dtb, worst, gam, mr, a(2), b, vt
        integer :: k, j

        p = pair_of(3)
        gam = 3.0_dp**2*e_charge**4*lnl/(4.0_dp*acos(-1.0_dp)*homma_eps0**2*ma**2)
        mr = ma/mb
        dtb = 1.0e-3_dp*t0
        call fp_rel_speed(t0, mb, 2.0_dp, vt)
        worst = 0.0_dp
        do k = 1, 4
            v = 0.05_dp*vt*3.0_dp**(k - 1)
            call fp_rel_speed(t0 + dtb, mb, 2.0_dp, vp)
            call fp_rel_speed(t0 - dtb, mb, 2.0_dp, vm)
            do j = 1, 4
                select case (j)
                case (1)
                    call fp_rel_friction(v, vp, nb, gam, mr, a(1))
                    call fp_rel_friction(v, vm, nb, gam, mr, a(2))
                    call scr_friction_dt(v, t0, mb, nb, gam, mr, b)
                case (2)
                    call fp_rel_diff_par(v, vp, nb, gam, a(1))
                    call fp_rel_diff_par(v, vm, nb, gam, a(2))
                    call scr_diff_par_dt(v, t0, mb, nb, gam, b)
                case (3)
                    call fp_rel_diff_perp(v, vp, nb, gam, a(1))
                    call fp_rel_diff_perp(v, vm, nb, gam, a(2))
                    call scr_diff_perp_dt(v, t0, mb, nb, gam, b)
                case default
                    call fp_rel_cross(v, vp, nb, gam, a(1))
                    call fp_rel_cross(v, vm, nb, gam, a(2))
                    call scr_cross_dt(v, t0, mb, nb, gam, b)
                end select
                worst = max(worst, abs((a(1) - a(2))/(2.0_dp*dtb) - b)/abs(b))
            end do
        end do
        print '(a,es12.4)', 'largest relative difference to the T_b difference quotient ', worst
        if (worst > 1.0e-5_dp) error stop 'T_b derivative leaves disagree with B04 coefficients'
    end subroutine check_temperature_derivative

    subroutine check_transient()
        ! grad T = 0, heat flux from the B05c gradient: over Omega t <= 0.01 the mean
        ! velocity is the collisionless rotation of the initial rate a,
        ! (v_x - i v_y)(t) = a_zeta (exp(i Omega t) - 1)/(i Omega), v_z = a_z t,
        ! up to collisional corrections of relative order nu t.
        type(scr_t) :: s
        real(dp) :: omega, q(3), a(3), vel(3, 0:200), t, worst
        complex(dp) :: w
        integer :: k

        call setup(s, 3, 16, omega)
        q = closure([300.0_dp, 0.0_dp, 5.0_dp]*e_charge, 1.0_dp)
        call scr_rate0(s, omega, [0.0_dp, 0.0_dp, 0.0_dp], q, .false., a)
        call scr_transient(s, omega, [0.0_dp, 0.0_dp, 0.0_dp], q, .false., 0.01_dp/omega, &
            200, vel)
        worst = 0.0_dp
        do k = 50, 200, 50
            t = k*0.01_dp/(omega*200)
            w = cmplx(a(1), -a(2), dp)*(exp(cmplx(0.0_dp, omega*t, dp)) - 1.0_dp)/ &
                cmplx(0.0_dp, omega, dp)
            worst = max(worst, abs(vel(1, k) - real(w))/abs(w), abs(vel(2, k) + aimag(w))/abs(w), &
                abs(vel(3, k) - a(3)*t)/abs(a(3)*t))
            print '(a,f7.4,a,3es15.7)', 'Omega t=', omega*t, ' <v>=', vel(:, k)
            print '(a,7x,a,3es15.7)', '      ', ' rot.=', real(w), -aimag(w), a(3)*t
        end do
        print '(a,es11.3)', 'largest relative difference ', worst
        if (worst > 2.0e-4_dp) error stop 'short-time transient is not the rotated initial rate'
    end subroutine check_transient

    subroutine check_steady_transient()
        ! W10+ (largest nu/Omega): the slope of <X> = <x> + <v_y>/Omega fitted over
        ! 100 <= Omega t <= 700 from TR-BDF2 agrees with the steady solution.
        integer, parameter :: nsteps = 14000
        type(scr_t) :: s
        real(dp) :: omega, g(3), q(3), d(3), dt, sx, st, stt, sxt, slope, rel
        real(dp), allocatable :: vel(:, :), x(:)
        integer :: k, m

        allocate (vel(3, 0:nsteps), x(0:nsteps))
        g = [grad100, 0.0_dp, 0.0_dp]
        q = closure(g, 1.0_dp)
        call setup(s, 10, 16, omega)
        call scr_steady(s, omega, g, q, .false., d)
        call scr_transient(s, omega, g, q, .false., 700.0_dp/omega, nsteps, vel)
        dt = 700.0_dp/omega/nsteps
        x(0) = 0.0_dp
        do k = 1, nsteps
            x(k) = x(k - 1) + 0.5_dp*dt*(vel(1, k) + vel(1, k - 1))
        end do
        x = x + vel(2, :)/omega
        m = 0; st = 0; sx = 0; stt = 0; sxt = 0
        do k = nsteps/7, nsteps
            m = m + 1
            st = st + k*dt; sx = sx + x(k); stt = stt + (k*dt)**2; sxt = sxt + k*dt*x(k)
        end do
        slope = (m*sxt - st*sx)/(m*stt - st**2)
        rel = abs(slope/d(1) - 1.0_dp)
        print '(a,es16.8,a,es16.8,a,es10.3)', 'fit ', slope, ' steady ', d(1), ' rel ', rel
        if (rel > 1.0e-3_dp) error stop 'transient slope differs from the steady drift'
    end subroutine check_steady_transient

end program test_screening_moments
