! Homogeneous test-particle Fokker-Planck relaxation in the Legendre modes
! l = 0, 1, 2 of the test-ion distribution (Homma 2014 Eq. 7, B04).
! Every pointwise coefficient, flux weight and moment weight comes from the
! FortSym leaves in kin6d_fp_relaxation_leaves; this module holds only the
! cell bookkeeping and the linear algebra.
module kin6d_fp_relaxation
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre, gauss_legendre_ab
    use kin6d_fp_relaxation_leaves
    implicit none
    private
    public :: fp_relax_t, fp_relax_setup, fp_relax_initial, fp_relax_step, &
        fp_relax_run, fp_relax_moments, fp_relax_lmax, fp_relax_nmom, fp_relax_cells

    integer, parameter :: fp_relax_lmax = 2
    integer, parameter :: fp_relax_nmom = 5

    ! Cell i spans [vface(i-1), vface(i)] with uniform width h. The state f(i, l)
    ! is the constant value of the Legendre mode l on that cell, so the cell
    ! integral of v**2 f equals vol(i) * f(i, l). Faces 0 and ncell carry zero flux.
    type :: fp_relax_t
        integer :: ncell = 0
        real(dp) :: ma = 0.0_dp, mb = 0.0_dp, nb = 0.0_dp, tb = 0.0_dp
        real(dp) :: tpar = 0.0_dp, tperp = 0.0_dp, drift = 0.0_dp, dens = 0.0_dp
        real(dp) :: mr = 0.0_dp, gam = 0.0_dp, vt = 0.0_dp, va = 0.0_dp
        real(dp) :: muf = 0.0_dp, tau = 0.0_dp, vmax = 0.0_dp, h = 0.0_dp
        real(dp), allocatable :: vface(:), vc(:), vol(:), w2(:), w3(:), w4(:), w6(:)
        real(dp), allocatable :: pf(:), qf(:), dperp(:)
        real(dp), allocatable :: f(:, :)
    end type fp_relax_t
contains

    subroutine fp_relax_setup(s, ncell, ma, mb, qa, qb, nb, tb, lnl, eps0, &
            tpar, tperp, drift, dens, tail)
        ! The grid spans [0, vmax] with vmax set by the initial tail bound.
        type(fp_relax_t), intent(out) :: s
        integer, intent(in) :: ncell
        real(dp), intent(in) :: ma, mb, qa, qb, nb, tb, lnl, eps0
        real(dp), intent(in) :: tpar, tperp, drift, dens, tail
        real(dp), allocatable :: kf(:), ci(:), a(:), dparc(:), e(:), plo(:), phi(:), rf(:)
        real(dp) :: x
        integer :: k, i

        if (ncell < 4) error stop 'FP relaxation needs at least four cells'
        s%ncell = ncell; s%ma = ma; s%mb = mb; s%nb = nb; s%tb = tb
        s%tpar = tpar; s%tperp = tperp; s%drift = drift; s%dens = dens
        s%mr = ma/mb
        call fp_rel_gamma(qa, qb, lnl, eps0, ma, s%gam)
        call fp_rel_speed(tb, mb, 2.0_dp, s%vt)
        call fp_rel_speed(tb, ma, 3.0_dp, s%va)
        call fp_rel_xs(s%va, s%vt, x)
        call fp_rel_mu(x, s%muf)
        call fp_rel_tau_s(s%va, s%muf, s%gam, nb, s%mr, s%tau)
        call fp_rel_vmax(drift, tpar, ma, tail, s%vmax)
        s%h = s%vmax/real(ncell, dp)

        allocate (kf(0:ncell), s%vface(0:ncell), ci(1:ncell), s%vc(1:ncell))
        allocate (s%vol(1:ncell), s%w2(1:ncell), s%w3(1:ncell), s%w4(1:ncell))
        allocate (s%w6(1:ncell), s%pf(0:ncell), s%qf(0:ncell), s%dperp(1:ncell))
        allocate (s%f(1:ncell, 0:fp_relax_lmax))
        do k = 0, ncell
            kf(k) = real(k, dp)
        end do
        call fp_rel_face(kf, s%h, s%vface)
        do i = 1, ncell
            ci(i) = real(i, dp)
        end do
        call fp_rel_centre(ci, s%h, s%vc)
        call fp_rel_w2(s%vface(0:ncell - 1), s%vface(1:ncell), s%w2)
        s%vol = s%w2
        call fp_rel_w3(s%vface(0:ncell - 1), s%vface(1:ncell), s%w3)
        call fp_rel_w4(s%vface(0:ncell - 1), s%vface(1:ncell), s%w4)
        call fp_rel_w6(s%vface(0:ncell - 1), s%vface(1:ncell), s%w6)
        call fp_rel_diff_perp(s%vc, s%vt, nb, s%gam, s%dperp)

        ! Face flux J = -A f + (D_par f)' + E f: A and E at the face, and the
        ! product (D_par f)' as a difference of D_par f at the two adjacent centres.
        allocate (rf(1:ncell - 1), a(1:ncell - 1), dparc(1:ncell), e(1:ncell - 1))
        allocate (plo(1:ncell - 1), phi(1:ncell - 1))
        call fp_rel_diff_par(s%vc, s%vt, nb, s%gam, dparc)
        rf = s%vface(1:ncell - 1)
        call fp_rel_friction(rf, s%vt, nb, s%gam, s%mr, a)
        call fp_rel_cross(rf, s%vt, nb, s%gam, e)
        call fp_rel_face_lo(rf, a, e, dparc(1:ncell - 1), s%h, plo)
        call fp_rel_face_hi(rf, a, e, dparc(2:ncell), s%h, phi)
        s%pf = 0.0_dp; s%qf = 0.0_dp
        s%pf(1:ncell - 1) = plo
        s%qf(1:ncell - 1) = phi
        s%f = 0.0_dp
    end subroutine fp_relax_setup

    integer function fp_relax_cells(cells_per_unit, ma, tb, drift, tpar, tail) result(n)
        ! Cells across [0, vmax] for cells_per_unit cells per u_a = sqrt(2 T_b/m_a).
        integer, intent(in) :: cells_per_unit
        real(dp), intent(in) :: ma, tb, drift, tpar, tail
        real(dp) :: vmax, unit

        call fp_rel_vmax(drift, tpar, ma, tail, vmax)
        call fp_rel_speed(tb, ma, 2.0_dp, unit)
        n = ceiling(vmax/unit*real(cells_per_unit, dp))
    end function fp_relax_cells

    subroutine fp_relax_initial(s)
        ! Cell values from the drifting bi-Maxwellian: Legendre modes by
        ! Gauss-Legendre quadrature in mu, with the v**2 weight inside the leaves.
        type(fp_relax_t), intent(inout) :: s
        integer, parameter :: nmu = 32, nv = 6
        real(dp) :: xm(nmu), wm(nmu), xv(nv), wv(nv), acc(0:fp_relax_lmax), val
        integer :: i, j, m, l

        call gauss_legendre(nmu, xm, wm)
        do i = 1, s%ncell
            call gauss_legendre_ab(nv, s%vface(i - 1), s%vface(i), xv, wv)
            acc = 0.0_dp
            do j = 1, nv
                do m = 1, nmu
                    call fp_rel_vproj0(xv(j), xm(m), s%ma, s%tpar, s%tperp, s%drift, s%dens, val)
                    acc(0) = acc(0) + wv(j)*wm(m)*val
                    call fp_rel_vproj1(xv(j), xm(m), s%ma, s%tpar, s%tperp, s%drift, s%dens, val)
                    acc(1) = acc(1) + wv(j)*wm(m)*val
                    call fp_rel_vproj2(xv(j), xm(m), s%ma, s%tpar, s%tperp, s%drift, s%dens, val)
                    acc(2) = acc(2) + wv(j)*wm(m)*val
                end do
            end do
            do l = 0, fp_relax_lmax
                s%f(i, l) = acc(l)/s%vol(i)
            end do
        end do
    end subroutine fp_relax_initial

    subroutine fp_relax_step(s, dt)
        ! TR-BDF2 with gamma = 2 - sqrt(2), dt in seconds. The matrix M - cK
        ! has mass-weighted cell volumes on the diagonal; K is the per-mode
        ! operator with zero-flux faces at 0 and vmax.
        type(fp_relax_t), intent(inout) :: s
        real(dp), intent(in) :: dt
        integer :: n, l, i
        real(dp) :: c, a0, a1, a2, ll
        real(dp), allocatable :: kl(:), kd(:), ku(:), sink(:), fpad(:), kf(:)
        real(dp), allocatable :: lo(:), di(:), up(:), rhs(:), y(:), fnew(:)
        real(dp) :: f(1:s%ncell)

        n = s%ncell
        call fp_rel_trap_c(dt, c)
        call fp_rel_bdf2_a0(dt, a0)
        call fp_rel_bdf2_a1(dt, a1)
        call fp_rel_bdf2_a2(dt, a2)
        allocate (kl(n), kd(n), ku(n), sink(n), fpad(0:n + 1), kf(n))
        allocate (lo(n), di(n), up(n), rhs(n), y(n), fnew(n))
        do l = 0, fp_relax_lmax
            ll = real(l*(l + 1), dp)
            call fp_rel_sink(s%h, s%dperp, ll, sink)
            do i = 1, n
                kl(i) = -s%pf(i - 1)
                kd(i) = s%pf(i) - s%qf(i - 1) - sink(i)
                ku(i) = s%qf(i)
            end do
            f = s%f(:, l)
            fpad = 0.0_dp
            fpad(1:n) = f
            do i = 1, n
                kf(i) = kl(i)*fpad(i - 1) + kd(i)*fpad(i) + ku(i)*fpad(i + 1)
            end do
            ! Stage 1: trapezoidal step over gamma*dt.
            lo = -c*kl; up = -c*ku; di = s%vol - c*kd
            rhs = s%vol*f + c*kf
            call thomas(n, lo, di, up, rhs, y)
            ! Stage 2: BDF2 from the three nodes of the TR-BDF2 step.
            lo = -kl; up = -ku; di = a0*s%vol - kd
            rhs = -s%vol*(a1*y + a2*f)
            call thomas(n, lo, di, up, rhs, fnew)
            s%f(:, l) = fnew
        end do
    end subroutine fp_relax_step

    subroutine fp_relax_run(s, times, dt_tau, mom)
        ! Integrates to each output time (in units of tau_s), with steps of at
        ! most dt_tau * tau_s, and returns the moments at those times.
        type(fp_relax_t), intent(inout) :: s
        real(dp), intent(in) :: times(:), dt_tau
        real(dp), intent(out) :: mom(fp_relax_nmom, size(times))
        real(dp) :: t, target_t, step
        integer :: it, nsteps, k

        t = 0.0_dp
        do it = 1, size(times)
            target_t = times(it)
            if (target_t > t) then
                nsteps = max(1, nint((target_t - t)/dt_tau))
                step = (target_t - t)/real(nsteps, dp)
                do k = 1, nsteps
                    call fp_relax_step(s, step*s%tau)
                end do
                t = target_t
            end if
            call fp_relax_moments(s, mom(:, it))
        end do
    end subroutine fp_relax_run

    subroutine fp_relax_moments(s, mom)
        ! mom = (n, u_z, p_par, p_perp, <v^4>) in SI, with n per m^3 of the
        ! normalised initial density.
        type(fp_relax_t), intent(in) :: s
        real(dp), intent(out) :: mom(fp_relax_nmom)
        real(dp) :: a(s%ncell), b(s%ncell), c(s%ncell), d(s%ncell), n, nu, pz, pp, v4

        call fp_rel_mom_n(s%w2, s%f(:, 0), a)
        call fp_rel_mom_u(s%w3, s%f(:, 1), b)
        call fp_rel_mom_pz(s%w4, s%f(:, 0), s%f(:, 2), c)
        call fp_rel_mom_pperp(s%w4, s%f(:, 0), s%f(:, 2), d)
        n = sum(a)
        nu = sum(b)
        pz = sum(c)
        pp = sum(d)
        call fp_rel_mom_v4(s%w6, s%f(:, 0), a)
        v4 = sum(a)
        mom(1) = n
        mom(2) = nu/n
        mom(3) = s%ma*(pz - nu*nu/n)
        mom(4) = 0.5_dp*s%ma*pp
        mom(5) = v4/n
    end subroutine fp_relax_moments

    subroutine thomas(n, lo, di, up, rhs, x)
        ! Tridiagonal solve without pivoting. No pivot failure occurred in the
        ! tests; the convergence checks would show a breakdown as a lost order.
        integer, intent(in) :: n
        real(dp), intent(in) :: lo(n), di(n), up(n), rhs(n)
        real(dp), intent(out) :: x(n)
        real(dp) :: cp(n), dq(n), den
        integer :: i

        cp(1) = up(1)/di(1)
        dq(1) = rhs(1)/di(1)
        do i = 2, n
            den = di(i) - lo(i)*cp(i - 1)
            cp(i) = up(i)/den
            dq(i) = (rhs(i) - lo(i)*dq(i - 1))/den
        end do
        x(n) = dq(n)
        do i = n - 1, 1, -1
            x(i) = dq(i) - cp(i)*x(i + 1)
        end do
    end subroutine thomas

end module kin6d_fp_relaxation
