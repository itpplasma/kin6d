module milca_rosenbluth_reference
    ! Rosenbluth potentials H=int f/|v-v'|, G=int f|v-v'| of a background that is
    ! a Maxwellian (l=0) plus a heat-flux distortion (l=1, linear in q), and the
    ! Fokker-Planck coefficients A=Gamma(1+m_a/m_b) grad H, D=Gamma/2 hess G.
    ! Radial integrands and their v-derivatives come from kin6d_homma_leaves.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use kin6d_homma_leaves
    implicit none
    private
    public :: homma_pair_t, homma_coeff_t
    public :: homma_vt, homma_distortion_c, homma_gamma, homma_gamma_prime
    public :: homma_radial_quad, homma_point
    public :: homma_force_maxwellian, homma_force_shell, homma_force_beam
    public :: homma_closure, homma_eps0

    real(dp), parameter :: homma_eps0 = 8.8541878128e-12_dp
    integer, parameter :: gl_order = 16
    integer, parameter :: i_h0 = 1, i_dh0 = 2, i_g0 = 3, i_dg0 = 4, i_d2g0 = 5
    integer, parameter :: i_h1 = 6, i_dh1 = 7, i_g1 = 8, i_dg1 = 9, i_d2g1 = 10

    type :: homma_pair_t
        real(dp) :: n_b = 0.0_dp ! background density (m^-3)
        real(dp) :: t_b = 0.0_dp ! background temperature (J)
        real(dp) :: m_b = 0.0_dp ! background mass (kg)
        real(dp) :: q_b = 0.0_dp ! background charge (C)
        real(dp) :: m_a = 0.0_dp ! test mass (kg)
        real(dp) :: q_a = 0.0_dp ! test charge (C)
        real(dp) :: lnlambda = 0.0_dp ! Coulomb logarithm
        real(dp) :: heat_flux(3) = 0.0_dp ! background heat flux q (W m^-2)
    end type homma_pair_t

    type :: homma_coeff_t
        real(dp) :: h_m = 0.0_dp, h_t = 0.0_dp ! H, Maxwellian and distortion parts
        real(dp) :: g_m = 0.0_dp, g_t = 0.0_dp ! G, Maxwellian and distortion parts
        real(dp) :: gh_m(3) = 0.0_dp, gh_t(3) = 0.0_dp
        real(dp) :: dg_m(3, 3) = 0.0_dp, dg_t(3, 3) = 0.0_dp
        real(dp) :: a_m(3) = 0.0_dp, a_t(3) = 0.0_dp
        real(dp) :: d_m(3, 3) = 0.0_dp, d_t(3, 3) = 0.0_dp
        ! Error estimates of a_m, a_t (m s^-2) and d_m, d_t (m^2 s^-3).
        real(dp) :: err_a_m = 0.0_dp, err_a_t = 0.0_dp, err_d_m = 0.0_dp, err_d_t = 0.0_dp
    end type homma_coeff_t

contains

    real(dp) function homma_vt(p)
        type(homma_pair_t), intent(in) :: p
        call homma_pair_vt(p%t_b, p%m_b, homma_vt)
    end function homma_vt

    real(dp) function homma_distortion_c(p)
        type(homma_pair_t), intent(in) :: p
        call homma_pair_c(p%t_b, p%m_b, p%n_b, homma_distortion_c)
    end function homma_distortion_c

    !> Gamma = q_a^2 q_b^2 lnL / (4 pi eps0^2 m_a^2).
    real(dp) function homma_gamma(p)
        type(homma_pair_t), intent(in) :: p
        call homma_pair_gamma(p%q_a, p%q_b, p%lnlambda, homma_eps0, p%m_a, homma_gamma)
    end function homma_gamma

    !> Gamma (1+m_a/m_b): prefactor of the drift A = Gamma(1+m_a/m_b) grad H.
    real(dp) function homma_gamma_prime(p)
        type(homma_pair_t), intent(in) :: p
        call homma_pair_gamma_prime(p%q_a, p%q_b, p%lnlambda, homma_eps0, p%m_a, p%m_b, &
            homma_gamma_prime)
    end function homma_gamma_prime

    !> Ten radial integrals at speed v by Gauss-Legendre panels on [0,v] (s<v) and
    !> [v, v+xmax vt] (s>v). Index order follows i_h0 ... i_d2g1. At v = 0 only the
    !> nonzero limits h0, g0, g0'' and h1' are set.
    subroutine homma_radial_quad(p, v, ppu, xmax, rv)
        type(homma_pair_t), intent(in) :: p
        real(dp), intent(in) :: v, ppu, xmax
        real(dp), intent(out) :: rv(10)
        real(dp) :: vt, c, n, xg(gl_order), wg(gl_order), h, s, wt, f(10), a
        integer :: np, k, j

        if (v < 0.0_dp) error stop 'homma_radial_quad: requires v >= 0'
        vt = homma_vt(p)
        c = homma_distortion_c(p)
        n = p%n_b
        call gauss_legendre(gl_order, xg, wg)
        rv = 0.0_dp
        if (.not. v > 0.0_dp) then
            ! The v -> 0 limits; the remaining entries vanish there.
            np = max(1, ceiling(xmax*ppu))
            h = xmax*vt/np
            do k = 1, np
                a = (k - 1)*h
                do j = 1, gl_order
                    s = a + (xg(j) + 1.0_dp)*h/2.0_dp
                    wt = wg(j)*h/2.0_dp
                    call homma_ih0_zero(v, s, vt, n, c, f(i_h0))
                    call homma_ig0_zero(v, s, vt, n, c, f(i_g0))
                    call homma_ig0_zero_d2(v, s, vt, n, c, f(i_d2g0))
                    call homma_ih1_zero_d1(v, s, vt, n, c, f(i_dh1))
                    rv(i_h0) = rv(i_h0) + wt*f(i_h0)
                    rv(i_g0) = rv(i_g0) + wt*f(i_g0)
                    rv(i_d2g0) = rv(i_d2g0) + wt*f(i_d2g0)
                    rv(i_dh1) = rv(i_dh1) + wt*f(i_dh1)
                end do
            end do
            return
        end if

        np = max(1, ceiling(v/vt*ppu))
        h = v/np
        do k = 1, np
            a = (k - 1)*h
            do j = 1, gl_order
                s = a + (xg(j) + 1.0_dp)*h/2.0_dp
                wt = wg(j)*h/2.0_dp
                call lo_leaves(v, s, vt, n, c, f)
                rv = rv + wt*f
            end do
        end do

        np = max(1, ceiling(xmax*ppu))
        h = xmax*vt/np
        do k = 1, np
            a = v + (k - 1)*h
            do j = 1, gl_order
                s = a + (xg(j) + 1.0_dp)*h/2.0_dp
                wt = wg(j)*h/2.0_dp
                call hi_leaves(v, s, vt, n, c, f)
                rv = rv + wt*f
            end do
        end do
    end subroutine homma_radial_quad

    subroutine lo_leaves(v, s, vt, n, c, f)
        real(dp), intent(in) :: v, s, vt, n, c
        real(dp), intent(out) :: f(10)
        call homma_ih0_lo(v, s, vt, n, c, f(i_h0))
        call homma_ih0_lo_d1(v, s, vt, n, c, f(i_dh0))
        call homma_ig0_lo(v, s, vt, n, c, f(i_g0))
        call homma_ig0_lo_d1(v, s, vt, n, c, f(i_dg0))
        call homma_ig0_lo_d2(v, s, vt, n, c, f(i_d2g0))
        call homma_ih1_lo(v, s, vt, n, c, f(i_h1))
        call homma_ih1_lo_d1(v, s, vt, n, c, f(i_dh1))
        call homma_ig1_lo(v, s, vt, n, c, f(i_g1))
        call homma_ig1_lo_d1(v, s, vt, n, c, f(i_dg1))
        call homma_ig1_lo_d2(v, s, vt, n, c, f(i_d2g1))
    end subroutine lo_leaves

    subroutine hi_leaves(v, s, vt, n, c, f)
        real(dp), intent(in) :: v, s, vt, n, c
        real(dp), intent(out) :: f(10)
        call homma_ih0_hi(v, s, vt, n, c, f(i_h0))
        call homma_ih0_hi_d1(v, s, vt, n, c, f(i_dh0))
        call homma_ig0_hi(v, s, vt, n, c, f(i_g0))
        call homma_ig0_hi_d1(v, s, vt, n, c, f(i_dg0))
        call homma_ig0_hi_d2(v, s, vt, n, c, f(i_d2g0))
        call homma_ih1_hi(v, s, vt, n, c, f(i_h1))
        call homma_ih1_hi_d1(v, s, vt, n, c, f(i_dh1))
        call homma_ig1_hi(v, s, vt, n, c, f(i_g1))
        call homma_ig1_hi_d1(v, s, vt, n, c, f(i_dg1))
        call homma_ig1_hi_d2(v, s, vt, n, c, f(i_d2g1))
    end subroutine hi_leaves

    !> Coefficients at velocity v (vector), v = 0 included as the exact limit.
    !> err_* estimate the error of the reported coefficient in its own units: the
    !> largest component change between (ppu, xmax), (2 ppu, xmax) and (ppu, xmax+1).
    subroutine homma_point(p, v3, ppu, xmax, co)
        type(homma_pair_t), intent(in) :: p
        real(dp), intent(in) :: v3(3), ppu, xmax
        type(homma_coeff_t), intent(out) :: co
        real(dp) :: rv(10), rf(10), rw(10), vm, gam, gdiff
        real(dp) :: ghm(3, 3), dgm(3, 3, 3), ght(3, 3), dgt(3, 3, 3)
        integer :: k

        vm = sqrt(dot_product(v3, v3))
        call homma_radial_quad(p, vm, ppu, xmax, rv)
        call homma_radial_quad(p, vm, 2.0_dp*ppu, xmax, rf)
        call homma_radial_quad(p, vm, ppu, xmax + 1.0_dp, rw)
        call assemble(v3, vm, p%heat_flux, rv, ghm(:, 1), dgm(:, :, 1), ght(:, 1), dgt(:, :, 1))
        call assemble(v3, vm, p%heat_flux, rf, ghm(:, 2), dgm(:, :, 2), ght(:, 2), dgt(:, :, 2))
        call assemble(v3, vm, p%heat_flux, rw, ghm(:, 3), dgm(:, :, 3), ght(:, 3), dgt(:, :, 3))

        co%h_m = rv(i_h0)
        co%g_m = rv(i_g0)
        co%gh_m = ghm(:, 1)
        co%dg_m = dgm(:, :, 1)
        if (vm > 0.0_dp) then
            co%h_t = dot_product(p%heat_flux, v3)*rv(i_h1)/vm
            co%g_t = dot_product(p%heat_flux, v3)*rv(i_g1)/vm
        end if
        co%gh_t = ght(:, 1)
        co%dg_t = dgt(:, :, 1)

        gam = homma_gamma_prime(p)
        gdiff = homma_gamma(p)
        co%a_m = gam*co%gh_m
        co%a_t = gam*co%gh_t
        co%d_m = 0.5_dp*gdiff*co%dg_m
        co%d_t = 0.5_dp*gdiff*co%dg_t
        co%err_a_m = 0.0_dp; co%err_a_t = 0.0_dp; co%err_d_m = 0.0_dp; co%err_d_t = 0.0_dp
        do k = 2, 3
            co%err_a_m = max(co%err_a_m, gam*maxval(abs(ghm(:, k) - ghm(:, 1))))
            co%err_a_t = max(co%err_a_t, gam*maxval(abs(ght(:, k) - ght(:, 1))))
            co%err_d_m = max(co%err_d_m, 0.5_dp*gdiff*maxval(abs(dgm(:, :, k) - dgm(:, :, 1))))
            co%err_d_t = max(co%err_d_t, 0.5_dp*gdiff*maxval(abs(dgt(:, :, k) - dgt(:, :, 1))))
        end do
    end subroutine homma_point

    !> grad H and Hess G, Maxwellian and distortion parts, from the radial integrals.
    subroutine assemble(v3, vm, q, rv, ghm, dgm, ght, dgt)
        real(dp), intent(in) :: v3(3), vm, q(3), rv(10)
        real(dp), intent(out) :: ghm(3), dgm(3, 3), ght(3), dgt(3, 3)
        real(dp) :: o(18)
        integer, parameter :: pi1(6) = [1, 2, 3, 1, 1, 2], pj1(6) = [1, 2, 3, 2, 3, 3]
        integer :: k

        if (vm > 0.0_dp) then
            call homma_assemble(v3(1), v3(2), v3(3), q(1), q(2), q(3), rv(i_dh0), &
                rv(i_dg0), rv(i_d2g0), rv(i_h1), rv(i_dh1), rv(i_g1), rv(i_dg1), &
                rv(i_d2g1), o(1), o(2), o(3), o(4), o(5), o(6), o(7), o(8), o(9), o(10), &
                o(11), o(12), o(13), o(14), o(15), o(16), o(17), o(18))
        else
            call homma_assemble_zero(q(1), q(2), q(3), rv(i_d2g0), rv(i_dh1), o(1), o(2), &
                o(3), o(4), o(5), o(6), o(7), o(8), o(9), o(10), o(11), o(12), o(13), &
                o(14), o(15), o(16), o(17), o(18))
        end if
        ghm = o(1:3)
        ght = o(10:12)
        do k = 1, 6
            dgm(pi1(k), pj1(k)) = o(3 + k)
            dgm(pj1(k), pi1(k)) = o(3 + k)
            dgt(pi1(k), pj1(k)) = o(12 + k)
            dgt(pj1(k), pi1(k)) = o(12 + k)
        end do
    end subroutine assemble

    !> Direction average of grad H_T per unit q at speed v.
    real(dp) function isotropic_ht(p, v, ppu, xmax) result(avg)
        type(homma_pair_t), intent(in) :: p
        real(dp), intent(in) :: v, ppu, xmax
        real(dp) :: rv(10)
        call homma_radial_quad(p, v, ppu, xmax, rv)
        call homma_isotropic_ht(v, rv(i_h1), rv(i_dh1), avg)
    end function isotropic_ht

    !> Maxwellian test species of temperature t_a: m_a Gamma' <grad H_T>. err is the
    !> largest component change under (2 ppu, xmax) and (ppu, xmax+1), in newton.
    subroutine homma_force_maxwellian(p, t_a, ppu, xmax, force, err)
        type(homma_pair_t), intent(in) :: p
        real(dp), intent(in) :: t_a, ppu, xmax
        real(dp), intent(out) :: force(3), err
        real(dp) :: f2(3), f3(3)
        force = maxwellian_mean(ppu, xmax)
        f2 = maxwellian_mean(2.0_dp*ppu, xmax)
        f3 = maxwellian_mean(ppu, xmax + 1.0_dp)
        err = max(maxval(abs(f2 - force)), maxval(abs(f3 - force)))
    contains
        function maxwellian_mean(pp, xm) result(f)
            real(dp), intent(in) :: pp, xm
            real(dp) :: f(3), va, xg(gl_order), wg(gl_order), wt, v, fa, h, a, mean
            integer :: np, k, j
            call homma_pair_vt(t_a, p%m_a, va)
            call gauss_legendre(gl_order, xg, wg)
            mean = 0.0_dp
            np = max(1, ceiling(xm*pp))
            h = xm*va/np
            do k = 1, np
                a = (k - 1)*h
                do j = 1, gl_order
                    v = a + (xg(j) + 1.0_dp)*h/2.0_dp
                    wt = wg(j)*h/2.0_dp
                    call homma_test_maxwell(v, va, fa)
                    mean = mean + wt*fa*isotropic_ht(p, v, pp, xm)
                end do
            end do
            f = p%m_a*homma_gamma_prime(p)*mean*p%heat_flux
        end function maxwellian_mean
    end subroutine homma_force_maxwellian

    !> Isotropic shell of normalised speed x_s.
    subroutine homma_force_shell(p, x_s, ppu, xmax, force, err)
        type(homma_pair_t), intent(in) :: p
        real(dp), intent(in) :: x_s, ppu, xmax
        real(dp), intent(out) :: force(3), err
        real(dp) :: v, f0, f2, f3, scale
        v = x_s*homma_vt(p)
        scale = p%m_a*homma_gamma_prime(p)
        f0 = isotropic_ht(p, v, ppu, xmax)
        f2 = isotropic_ht(p, v, 2.0_dp*ppu, xmax)
        f3 = isotropic_ht(p, v, ppu, xmax + 1.0_dp)
        force = scale*f0*p%heat_flux
        err = scale*max(abs(f2 - f0), abs(f3 - f0))*maxval(abs(p%heat_flux))
    end subroutine homma_force_shell

    !> Beam at normalised speed x_b along unit vector e.
    subroutine homma_force_beam(p, x_b, e, ppu, xmax, force, err)
        type(homma_pair_t), intent(in) :: p
        real(dp), intent(in) :: x_b, e(3), ppu, xmax
        real(dp), intent(out) :: force(3), err
        type(homma_coeff_t) :: co
        call homma_point(p, x_b*homma_vt(p)*e, ppu, xmax, co)
        force = p%m_a*co%a_t
        err = p%m_a*co%err_a_t
    end subroutine homma_force_beam

    !> Homma closure for the background heat flux from grad T (J m^-1).
    subroutine homma_closure(p, b_hat, bmag, grad_t, q_out)
        type(homma_pair_t), intent(in) :: p
        real(dp), intent(in) :: b_hat(3), bmag, grad_t(3)
        real(dp), intent(out) :: q_out(3)
        real(dp) :: lnl
        lnl = p%lnlambda
        call homma_closure_qx(grad_t(1), grad_t(2), grad_t(3), b_hat(1), b_hat(2), b_hat(3), &
            p%n_b, p%t_b, p%m_b, p%q_b, bmag, lnl, homma_eps0, q_out(1))
        call homma_closure_qy(grad_t(1), grad_t(2), grad_t(3), b_hat(1), b_hat(2), b_hat(3), &
            p%n_b, p%t_b, p%m_b, p%q_b, bmag, lnl, homma_eps0, q_out(2))
        call homma_closure_qz(grad_t(1), grad_t(2), grad_t(3), b_hat(1), b_hat(2), b_hat(3), &
            p%n_b, p%t_b, p%m_b, p%q_b, bmag, lnl, homma_eps0, q_out(3))
    end subroutine homma_closure

end module milca_rosenbluth_reference
