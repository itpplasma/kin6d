program test_rosenbluth_reference
    ! Named checks of the Rosenbluth oracle: maxwellian, l1_brute, invariants,
    ! slow_limit. Each stops with a non-zero exit code on failure.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use milca_homma_leaves
    use milca_rosenbluth_reference
    implicit none
    real(dp), parameter :: e_charge = 1.602176634e-19_dp
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: amu = 1.66053906660e-27_dp
    character(len=32) :: name

    if (command_argument_count() < 1) error stop 'usage: test_rosenbluth_reference NAME'
    call get_command_argument(1, name)
    select case (trim(name))
    case ('maxwellian')
        call run_maxwellian()
    case ('l1_brute')
        call run_l1_brute()
    case ('invariants')
        call run_invariants()
    case ('slow_limit')
        call run_slow_limit()
    case default
        error stop 'unknown test name'
    end select
    print '(a)', 'PASS ' // trim(name)

contains

    subroutine fail_if(cond, what)
        logical, intent(in) :: cond
        character(*), intent(in) :: what
        if (cond) then
            print '(a)', 'FAIL ' // what
            error stop 1
        end if
    end subroutine fail_if

    !> Relative difference against a reference value b.
    pure real(dp) function rel(a, b)
        real(dp), intent(in) :: a, b
        rel = abs(a - b)/max(abs(b), tiny(1.0_dp))
    end function rel

    !> Background of benchmark case B01 (H+, n=1e20 m^-3, T=50 eV, lnL=15) with
    !> a W3+ test species.
    function bench_pair(heat_flux) result(p)
        real(dp), intent(in) :: heat_flux(3)
        type(homma_pair_t) :: p
        p%n_b = 1.0e20_dp
        p%t_b = 50.0_dp*e_charge
        p%m_b = 1.67262192369e-27_dp
        p%q_b = e_charge
        p%m_a = 183.84_dp*amu
        p%q_a = 3.0_dp*e_charge
        p%lnlambda = 15.0_dp
        p%heat_flux = heat_flux
    end function bench_pair

    !> Maxwellian background (q=0): radial integrals against closed erf forms.
    subroutine run_maxwellian()
        type(homma_pair_t) :: p
        real(dp) :: rv(10), vt, vk, h0c, h0pc, g0pc, g0ppc
        real(dp) :: speeds(4)
        integer :: k

        p = bench_pair([0.0_dp, 0.0_dp, 0.0_dp])
        vt = homma_vt(p)
        speeds = [0.3_dp, 1.0_dp, 2.5_dp, 4.0_dp]
        do k = 1, 4
            vk = speeds(k)*vt
            call homma_radial_quad(p, vk, 8.0_dp, 7.0_dp, rv)
            call homma_closed_h0(vk, vt, p%n_b, h0c)
            call homma_closed_h0p(vk, vt, p%n_b, h0pc)
            call homma_closed_g0p_over_v(vk, vt, p%n_b, g0pc)
            call homma_closed_g0pp_minus_g0p_over_v(vk, vt, p%n_b, g0ppc)
            call fail_if(rel(rv(1), h0c) > 1.0e-10_dp, 'maxwellian H0 = n erf(x)/v')
            call fail_if(rel(rv(2), h0pc) > 1.0e-10_dp, 'maxwellian dH0/dv')
            call fail_if(rel(rv(4)/vk, g0pc) > 1.0e-10_dp, 'maxwellian G0p/v')
            call fail_if(rel(rv(5) - rv(4)/vk, g0ppc) > 1.0e-10_dp, 'maxwellian G0pp - G0p/v')
        end do
    end subroutine run_maxwellian

    !> Independent 3D quadrature on spheres around v (w = v + r omega, Jacobian r^2,
    !> 1/|v-w| = 1/r). kind: 1 H_M, 2 G_M, 3 H_T, 4 G_T.
    function brute(p, v3, kind) result(total)
        type(homma_pair_t), intent(in) :: p
        real(dp), intent(in) :: v3(3)
        integer, intent(in) :: kind
        real(dp) :: total
        integer, parameter :: np_r = 40, ng = 16, nphi = 32
        real(dp) :: xr(ng), wr(ng), xm(ng), wm(ng), r, wt, a, h, rmax, mu, sn, phi
        real(dp) :: om(3), wv(3), out, vt, c
        integer :: k, j, i, l

        vt = homma_vt(p)
        c = homma_distortion_c(p)
        rmax = sqrt(dot_product(v3, v3)) + 10.0_dp*vt
        call gauss_legendre(ng, xr, wr)
        call gauss_legendre(ng, xm, wm)
        h = rmax/np_r
        total = 0.0_dp
        do k = 1, np_r
            a = (k - 1)*h
            do j = 1, ng
                r = a + (xr(j) + 1.0_dp)*h/2.0_dp
                wt = wr(j)*h/2.0_dp
                do i = 1, ng
                    mu = xm(i)
                    sn = sqrt(1.0_dp - mu**2)
                    do l = 1, nphi
                        phi = 2.0_dp*pi*(l - 1)/nphi
                        om = [sn*cos(phi), sn*sin(phi), mu]
                        wv = v3 + r*om
                        select case (kind)
                        case (1)
                            call homma_bf_hm(r, wv(1), wv(2), wv(3), vt, p%n_b, c, &
                                p%heat_flux(1), p%heat_flux(2), p%heat_flux(3), out)
                        case (2)
                            call homma_bf_gm(r, wv(1), wv(2), wv(3), vt, p%n_b, c, &
                                p%heat_flux(1), p%heat_flux(2), p%heat_flux(3), out)
                        case (3)
                            call homma_bf_ht(r, wv(1), wv(2), wv(3), vt, p%n_b, c, &
                                p%heat_flux(1), p%heat_flux(2), p%heat_flux(3), out)
                        case default
                            call homma_bf_gt(r, wv(1), wv(2), wv(3), vt, p%n_b, c, &
                                p%heat_flux(1), p%heat_flux(2), p%heat_flux(3), out)
                        end select
                        total = total + wt*wm(i)*(2.0_dp*pi/nphi)*out
                    end do
                end do
            end do
        end do
    end function brute

    subroutine run_l1_brute()
        type(homma_pair_t) :: p
        type(homma_coeff_t) :: co
        real(dp) :: v3(3), vt, h, fd, d2, trace
        integer :: i

        p = bench_pair([2.0e3_dp, -1.0e3_dp, 5.0e2_dp])
        vt = homma_vt(p)
        v3 = [0.7_dp, 0.2_dp, -0.4_dp]*vt
        call homma_point(p, v3, 8.0_dp, 7.0_dp, co)

        call fail_if(rel(brute(p, v3, 1), co%h_m) > 1.0e-6_dp, 'brute H_M')
        call fail_if(rel(brute(p, v3, 2), co%g_m) > 1.0e-6_dp, 'brute G_M')
        call fail_if(rel(brute(p, v3, 3), co%h_t) > 1.0e-6_dp, 'brute H_T')
        call fail_if(rel(brute(p, v3, 4), co%g_t) > 1.0e-6_dp, 'brute G_T')

        ! Gradient of H_T by fourth-order differences of the brute-force values.
        h = 1.0e-2_dp*vt
        do i = 1, 3
            fd = (-brute_shift(p, v3, i, 2*h, 3) + 8.0_dp*brute_shift(p, v3, i, h, 3) &
                - 8.0_dp*brute_shift(p, v3, i, -h, 3) + brute_shift(p, v3, i, -2*h, 3))/(12.0_dp*h)
            call fail_if(rel(fd, co%gh_t(i)) > 1.0e-6_dp, 'brute gradient of H_T')
        end do

        ! Diagonal of the Hessian of G_T by fourth-order differences.
        h = 1.0e-2_dp*vt
        trace = 0.0_dp
        do i = 1, 3
            d2 = (-brute_shift(p, v3, i, 2*h, 4) + 16.0_dp*brute_shift(p, v3, i, h, 4) &
                - 30.0_dp*brute(p, v3, 4) + 16.0_dp*brute_shift(p, v3, i, -h, 4) &
                - brute_shift(p, v3, i, -2*h, 4))/(12.0_dp*h**2)
            call fail_if(rel(d2, co%dg_t(i, i)) > 1.0e-5_dp, 'brute diagonal Hessian of G_T')
            trace = trace + co%dg_t(i, i)
        end do
        ! Laplacian of G_T equals 2 H_T.
        call fail_if(rel(trace, 2.0_dp*co%h_t) > 1.0e-8_dp, 'Laplacian G_T = 2 H_T')
    end subroutine run_l1_brute

    !> brute(v3 + shift e_i) for kind.
    function brute_shift(p, v3, i, shift, kind) result(val)
        type(homma_pair_t), intent(in) :: p
        real(dp), intent(in) :: v3(3), shift
        integer, intent(in) :: i, kind
        real(dp) :: val, x(3)
        x = v3
        x(i) = x(i) + shift
        val = brute(p, x, kind)
    end function brute_shift

    subroutine run_invariants()
        type(homma_pair_t) :: p
        type(homma_coeff_t) :: co
        real(dp) :: vt, v3(3), speeds(3), gam, trace_dm, scale, mom, mom_scale, heat
        real(dp) :: xg(16), wg(16), a, hh, x, pt, rv(10), dsym
        integer :: k, i, j, n

        p = bench_pair([2.0e3_dp, -1.0e3_dp, 5.0e2_dp])
        vt = homma_vt(p)
        speeds = [0.2_dp, 1.0_dp, 3.0_dp]
        gam = homma_gamma(p)
        do k = 1, 3
            v3 = [0.6_dp, -0.5_dp, 0.3_dp]*speeds(k)*vt
            call homma_point(p, v3, 8.0_dp, 7.0_dp, co)
            ! Symmetry of the Hessians.
            scale = max(maxval(abs(co%d_m)), tiny(1.0_dp))
            dsym = 0.0_dp
            do j = 1, 3
                do i = 1, 3
                    dsym = max(dsym, abs(co%d_m(i, j) - co%d_m(j, i)))
                    dsym = max(dsym, abs(co%d_t(i, j) - co%d_t(j, i)))
                end do
            end do
            call fail_if(dsym > 1.0e-10_dp*scale, 'symmetry of D')
            ! Positivity of the Maxwellian Hessian: G0'' >= 0 and G0'/v >= 0.
            call homma_radial_quad(p, speeds(k)*vt, 8.0_dp, 7.0_dp, rv)
            call fail_if(rv(5) < 0.0_dp .or. rv(4) < 0.0_dp, 'positivity of G_M curvature')
            ! Trace identity: tr D_M = Gamma H_M.
            trace_dm = co%d_m(1, 1) + co%d_m(2, 2) + co%d_m(3, 3)
            call fail_if(rel(trace_dm, gam*co%h_m) > 1.0e-10_dp, &
                'trace identity tr D_M = Gamma H_M')
        end do

        ! Momentum and heat flux of the distortion, from the generated profile:
        ! int w f_T = (4 pi/3) q int s^3 pt ds, (m/2) int w^2 w f_T = (m/2)(4 pi/3) q int s^5 pt ds.
        call gauss_legendre(16, xg, wg)
        n = 64
        hh = 14.0_dp*vt/n
        mom = 0.0_dp; mom_scale = 0.0_dp; heat = 0.0_dp
        do k = 1, n
            a = (k - 1)*hh
            do j = 1, 16
                x = a + (xg(j) + 1.0_dp)*hh/2.0_dp
                call homma_profile_t(0.0_dp, x, vt, p%n_b, homma_distortion_c(p), pt)
                mom = mom + wg(j)*hh/2.0_dp*x**3*pt
                mom_scale = mom_scale + wg(j)*hh/2.0_dp*x**3*abs(pt)
                heat = heat + wg(j)*hh/2.0_dp*x**5*pt
            end do
        end do
        call fail_if(abs(4.0_dp*pi/3.0_dp*mom) > 1.0e-10_dp*4.0_dp*pi/3.0_dp*mom_scale, &
            'momentum of distortion vanishes')
        call fail_if(rel(0.5_dp*p%m_b*4.0_dp*pi/3.0_dp*heat, 1.0_dp) > 1.0e-10_dp, &
            'heat flux of distortion equals q')
    end subroutine run_invariants

    subroutine run_slow_limit()
        type(homma_pair_t) :: p
        type(homma_coeff_t) :: co, co_small
        real(dp) :: vt, v3(3), cst, pref, expected(3)
        integer :: i

        p = bench_pair([2.0e3_dp, -1.0e3_dp, 5.0e2_dp])
        vt = homma_vt(p)
        ! v = 0 exactly: the limit path, not a small speed.
        v3 = 0.0_dp
        call homma_point(p, v3, 8.0_dp, 7.0_dp, co)
        call fail_if(maxval(abs(co%a_m)) > 0.0_dp .or. maxval(abs(co%d_t)) > 0.0_dp, &
            'A_M(0) and D_T(0) vanish')
        ! Continuity of the isotropic D_M(0) with D_M at a small speed.
        call homma_point(p, [1.0e-3_dp, 0.0_dp, 0.0_dp]*vt, 8.0_dp, 7.0_dp, co_small)
        do i = 1, 3
            call fail_if(rel(co%d_m(i, i), co_small%d_m(i, i)) > 1.0e-5_dp, 'D_M continuous at 0')
        end do
        ! m_a A_T(0) = -C q/(10 pi^{3/2} T^2 vt), C = lnL (q_a q_b/eps0)^2 (1 + m_b/m_a).
        cst = p%lnlambda*(p%q_a*p%q_b/homma_eps0)**2*(1.0_dp + p%m_b/p%m_a)
        pref = -cst/(10.0_dp*pi**1.5_dp*p%t_b**2*vt)
        expected = pref*p%heat_flux
        do i = 1, 3
            call fail_if(rel(p%m_a*co%a_t(i), expected(i)) > 1.0e-7_dp, &
                'slow limit m_a A_T(0)')
        end do
    end subroutine run_slow_limit

end program test_rosenbluth_reference
