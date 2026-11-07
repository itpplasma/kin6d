!> Collision operator on the energy grid (#184).
!!   maxwellian_null  test-particle operator at equal temperatures: the
!!                    Maxwellian is a null vector, density is conserved,
!!                    the matrix is symmetric and negative semi-definite,
!!                    on the coarsest and on a fine grid;
!!   rates            the deflection, parallel-diffusion and slowing-down
!!                    rates of the FortSym leaves against Rosenbluth
!!                    potentials of the Maxwellian by quadrature, and the
!!                    momentum and energy exchange of a drifting or hotter
!!                    Maxwellian with its reservoir against the closed forms;
!!   jvp              tangents with respect to temperatures and densities
!!                    against centred differences, test particle and Krook.
program test_dk_collop
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_dk_collop
    use kin6d_dk_collop_leaves, only: collop_tp, collop_tp_series
    use kin6d_rosenbluth_reference, only: homma_pair_t, homma_radial_quad, homma_gamma
    implicit none
    character(64) :: check
    logical :: ok
    real(dp), parameter :: pi = 4*atan(1.0_dp)
    real(dp), parameter :: me = 9.1093837015e-31_dp, mp = 1.67262192369e-27_dp
    real(dp), parameter :: kev = 1e3_dp*collop_qe

    call get_command_argument(1, check)
    select case (trim(check))
    case ('maxwellian_null'); ok = maxwellian_null()
    case ('rates'); ok = rates()
    case ('jvp'); ok = jvp()
    case default
        print '(a)', 'unknown check '//trim(check)
        ok = .false.
    end select
    if (.not. ok) error stop 1

contains

    !> Electrons and deuterons; equal temperatures unless tau_i is given.
    function plasma(t_ratio) result(sp)
        real(dp), intent(in) :: t_ratio
        type(collop_species_t) :: sp(2)
        sp(1) = collop_species_t(me, -1.0_dp, 3e19_dp, 4*kev)
        sp(2) = collop_species_t(2.01410178_dp*1.66053906660e-27_dp, 1.0_dp, 3e19_dp, &
            4*kev*t_ratio)
    end function plasma

    logical function maxwellian_null() result(ok)
        type(collop_grid_t) :: grid
        integer :: k, nelem(2), order(2), status
        nelem = [1, 6]
        order = [3, 6]
        ok = .true.
        do k = 1, 2
            call collop_make_grid(6.0_dp, nelem(k), order(k), grid, status)
            if (status /= COLLOP_OK) error stop 'grid failed'
            ok = null_checks(grid) .and. ok
        end do
    end function maxwellian_null

    logical function null_checks(grid) result(ok)
        type(collop_grid_t), intent(in) :: grid
        type(collop_t) :: op
        integer :: status, n, a, l, info
        real(dp) :: scale, null_defect, density_defect, asym, top
        real(dp), allocatable :: m(:, :), ev(:), work(:)
        real(dp) :: lnl(2, 2)

        lnl = 17
        call collop_build(plasma(1.0_dp), lnl, grid, COLLOP_TEST_PARTICLE, 4, op, status)
        if (status /= COLLOP_OK) error stop 'build failed'
        n = grid%n
        ok = .true.
        do l = 0, 4
            m = op%a(:, :, l)
            scale = maxval(abs(m))
            asym = maxval(abs(m - transpose(m)))/scale
            null_defect = 0
            density_defect = 0
            do a = 1, 2
                if (l == 0) then
                    null_defect = max(null_defect, maxval(abs(sum(m((a - 1)*n + 1:a*n, &
                        (a - 1)*n + 1:a*n), 2))))
                    density_defect = max(density_defect, maxval(abs(sum(m((a - 1)*n + 1:a*n, &
                        :), 1))))
                end if
            end do
            null_defect = null_defect/scale
            density_defect = density_defect/scale
            allocate (ev(2*n), work(10*n))
            call dsyev('N', 'U', 2*n, (m + transpose(m))/2, 2*n, ev, work, size(work), info)
            top = maxval(ev)/scale
            deallocate (ev, work)
            print '(a,i3,a,i2,a,4es10.2)', 'n =', n, ' l =', l, &
                ' null, density, asymmetry, max eigenvalue: ', null_defect, &
                density_defect, asym, top
            ok = ok .and. null_defect < 1e-13_dp .and. density_defect < 1e-13_dp &
                .and. asym < 1e-14_dp .and. top < 1e-14_dp .and. info == 0
        end do
    end function null_checks

    !> Leaf rates against the Rosenbluth potentials G, H of the partner
    !> Maxwellian by quadrature (homma): nu_D = Gamma G'/v^3,
    !> nu_par = Gamma G''/v^2 and the coefficient of h in the radial flux,
    !> -Gamma (G'' v/v_ta^2 + (m_a/m_b) H'), with Gamma = e_a^2 e_b^2 lnL/
    !> (4 pi eps0^2 m_a^2); the leaves carry them per nuhat = n_b Gamma/v_ta^3.
    logical function rates() result(ok)
        type(collop_species_t) :: sp(2)
        type(homma_pair_t) :: p
        real(dp) :: x(6), lor, en, th, rv(10), v, vta, nuhat, err(3), worst(3), rinv, tau
        integer :: a, b, k
        sp = plasma(0.7_dp)
        x = [0.05_dp, 0.3_dp, 0.9_dp, 1.7_dp, 3.0_dp, 5.0_dp]
        worst = 0
        do a = 1, 2
            do b = 1, 2
                p = homma_pair_t(n_b=sp(b)%n, t_b=sp(b)%t, m_b=sp(b)%m, &
                    q_b=sp(b)%z*collop_qe, m_a=sp(a)%m, q_a=sp(a)%z*collop_qe, lnlambda=17.0_dp)
                vta = collop_vth(sp(a))
                rinv = vta/collop_vth(sp(b))
                tau = sp(a)%t/sp(b)%t
                nuhat = sp(b)%n*homma_gamma(p)/vta**3
                do k = 1, size(x)
                    if (x(k)*rinv < 1) then
                        call collop_tp_series(x(k), rinv, tau, lor, en, th)
                    else
                        call collop_tp(x(k), rinv, tau, lor, en, th)
                    end if
                    v = x(k)*vta
                    call homma_radial_quad(p, v, 64.0_dp, 12.0_dp, rv)
                    lor = lor/(x(k)**2*exp(-x(k)**2))
                    en = 2*en/(x(k)**4*exp(-x(k)**2))
                    th = th/(x(k)**2*exp(-x(k)**2))*vta
                    err(1) = abs(lor/(homma_gamma(p)*rv(4)/v**3/nuhat) - 1)
                    err(2) = abs(en/(homma_gamma(p)*rv(5)/v**2/nuhat) - 1)
                    ! Relative to one of its two terms: the sum vanishes for
                    ! equal temperatures.
                    err(3) = abs(th + homma_gamma(p)*(rv(5)*v/vta**2 &
                        + sp(a)%m/sp(b)%m*rv(2))/nuhat)/(homma_gamma(p)*rv(5)*v/vta**2/nuhat)
                    worst = max(worst, err)
                end do
            end do
        end do
        print '(a,3es10.2)', 'leaf rates vs Rosenbluth quadrature (nu_D, nu_par, h-flux): ', &
            worst
        ok = all(worst < 1e-10_dp)
        ok = exchange(sp) .and. ok
    end function rates

    !> Momentum and energy taken by species a from a reservoir b at rest:
    !> R_ab = -n_a n_b e_a^2 e_b^2 lnL u/(3 (2 pi)^(3/2) eps0^2 m_ab S^3),
    !> Q_ab = n_a n_b e_a^2 e_b^2 lnL (T_b - T_a)/((2 pi)^(3/2) eps0^2 m_a m_b S^3),
    !> S^2 = T_a/m_a + T_b/m_b, m_ab = m_a m_b/(m_a + m_b).
    logical function exchange(sp) result(ok)
        type(collop_species_t), intent(in) :: sp(2)
        type(collop_t) :: op
        type(collop_grid_t) :: grid
        type(collop_species_t) :: pair(2)
        real(dp) :: lnl(2, 2), r, q, r_ref, q_ref, s3, c, u, worst(2), vta
        real(dp), allocatable :: h(:), y(:)
        integer :: a, b, status, n
        call collop_make_grid(7.0_dp, 6, 5, grid, status)
        n = grid%n
        lnl = 17
        u = 1e3_dp
        worst = 0
        do a = 1, 2
            do b = 1, 2
                ! Species a with a single partner b, as two entries so that
                ! a = b keeps its own reservoir.
                pair = [sp(a), sp(b)]
                lnl = 17
                call collop_build(pair, lnl, grid, COLLOP_TEST_PARTICLE, 1, op, status)
                if (status /= COLLOP_OK) error stop 'build failed'
                vta = collop_vth(sp(a))
                ! Only the a-on-b block: rows and columns of species 1, minus
                ! its self part, rebuilt for a pure partner b.
                call partner_only(pair, grid, op)
                h = 2*u/vta*grid%x
                y = matmul(op%a(1:n, 1:n, 1), h)
                r = sum(sp(a)%m*vta*grid%x*y)
                h = grid%x**0
                y = matmul(op%a(1:n, 1:n, 0), h)
                q = sum(sp(a)%t*grid%x**2*y)
                s3 = (sp(a)%t/sp(a)%m + sp(b)%t/sp(b)%m)**1.5_dp
                c = sp(a)%n*sp(b)%n*(sp(a)%z*sp(b)%z)**2*collop_qe**4*17 &
                    /((2*pi)**1.5_dp*collop_eps0**2)
                r_ref = -c*u/(3*s3)*(sp(a)%m + sp(b)%m)/(sp(a)%m*sp(b)%m)
                q_ref = c*(sp(b)%t - sp(a)%t)/(sp(a)%m*sp(b)%m*s3)
                if (a == b) then
                    print '(a,2i2,a,es22.14,a,es22.14)', 'pair', a, b, ' R', r, ' ref', r_ref
                    worst(1) = max(worst(1), abs(r/r_ref - 1))
                else
                    print '(a,2i2,a,2es22.14,a,2es22.14)', 'pair', a, b, ' R, Q', r, q, &
                        ' ref', r_ref, q_ref
                    worst(1) = max(worst(1), abs(r/r_ref - 1))
                    worst(2) = max(worst(2), abs(q/q_ref - 1))
                end if
            end do
        end do
        print '(a,2es10.2)', 'friction and heat exchange with a reservoir: ', worst
        ok = all(worst < 1e-12_dp)
    end function exchange

    !> Replace the operator by the test-particle part of species 1 on
    !> species 2 alone.
    subroutine partner_only(pair, grid, op)
        type(collop_species_t), intent(in) :: pair(2)
        type(collop_grid_t), intent(in) :: grid
        type(collop_t), intent(inout) :: op
        type(collop_t) :: self
        real(dp) :: lnl(1, 1)
        integer :: status, n
        n = grid%n
        lnl = 17
        call collop_build(pair(1:1), lnl, grid, COLLOP_TEST_PARTICLE, 1, self, status)
        if (status /= COLLOP_OK) error stop 'build failed'
        op%a(1:n, 1:n, :) = op%a(1:n, 1:n, :) - self%a
    end subroutine partner_only

    logical function jvp() result(ok)
        type(collop_grid_t) :: grid
        type(collop_t) :: op, up, down
        type(collop_species_t) :: sp(2), sp_up(2), sp_down(2)
        real(dp) :: lnl(2, 2), dt(2), dn(2), h, err
        real(dp), allocatable :: da(:, :, :), fd(:, :, :)
        integer :: status, model, k
        integer, parameter :: models(2) = [COLLOP_TEST_PARTICLE, COLLOP_KROOK]
        call collop_make_grid(6.0_dp, 3, 4, grid, status)
        lnl = 17
        sp = plasma(0.7_dp)
        dt = [0.3_dp*sp(1)%t, -0.2_dp*sp(2)%t]
        dn = [0.1_dp*sp(1)%n, 0.25_dp*sp(2)%n]
        h = 1e-5_dp
        ok = .true.
        do k = 1, 2
            model = models(k)
            call collop_build_jvp(sp, lnl, grid, model, 3, dt, dn, op, da, status, &
                nu_krook=1e4_dp)
            if (status /= COLLOP_OK) error stop 'jvp failed'
            sp_up = sp
            sp_down = sp
            sp_up%t = sp%t + h*dt
            sp_up%n = sp%n + h*dn
            sp_down%t = sp%t - h*dt
            sp_down%n = sp%n - h*dn
            call collop_build(sp_up, lnl, grid, model, 3, up, status, nu_krook=1e4_dp)
            call collop_build(sp_down, lnl, grid, model, 3, down, status, nu_krook=1e4_dp)
            fd = (up%a - down%a)/(2*h)
            err = maxval(abs(da - fd))/maxval(abs(fd))
            print '(a,i2,a,es10.2)', 'model', model, ' tangent vs centred difference:', err
            ok = ok .and. err < 1e-8_dp
        end do
    end function jvp

end program test_dk_collop
