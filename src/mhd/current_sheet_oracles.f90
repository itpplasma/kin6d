module kin6d_current_sheet_oracles
    ! Linear ideal current-sheet references, mu0 = 1
    ! (docs/current-sheet-equilibrium.md section 5). Oracles only: no solver
    ! uses them as forcing.
    !
    ! Slab (5.1): B = (0, B_y(x), B_z(x)), B_y = b0 (x/a + dplus) for x > 0 and
    ! b0 (x/a - dminus) for x < 0, walls x = +-a + w(+-) cos(ky). The marginal
    ! displacement is xi_x = xi(x) cos(ky), psi = k B_y xi. dplus = dminus = 0 is the
    ! Hahm-Kulsrud case with continuous iota (step in xi at x = 0).
    ! Local formulas are generated (gen/gen_current_sheet.f90).
    !
    ! Screw pinch (5.2): current-free annulus r0 < r < b with B_theta = c/r,
    ! B_z constant on each side of r = sigma, phase mtheta + kz, outer wall
    ! r = b + delta cos(phase). Modified Bessel functions are evaluated by
    ! series and recurrence here (hand-written: FortSym emits no Bessel
    ! runtime); admitted arguments are 0 < |k| r <= 2.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_current_sheet_leaves, only: sheet_slab_side, sheet_slab_interface
    implicit none
    private
    public :: slab_sheet_t, pinch_sheet_t, modified_bessel_ik

    type :: slab_sheet_t
        real(dp) :: a = 1, k = acos(-1.0_dp), b0 = 1
        real(dp) :: dplus = 0, dminus = 0   ! B_y jump: B_y(0+-) = +-b0 d(plus,minus)
        real(dp) :: wp = 0, wm = 0          ! wall displacements at x = +a, -a
        real(dp) :: xi0 = 0                 ! interface displacement (jump case)
        logical :: continuous = .true.      ! dplus = dminus = 0
        integer :: status = -1              ! 0 ready, 1 invalid
    contains
        procedure, non_overridable :: init => slab_init
        procedure, non_overridable :: eval => slab_eval
        procedure, non_overridable :: sheet => slab_sheet_current
        procedure, non_overridable :: by => slab_by
    end type

    type :: pinch_sheet_t
        real(dp) :: r0 = 0.2_dp, sigma = 0.5_dp, b = 1, k = -1.0_dp/3
        integer :: m = 1
        real(dp) :: cm = 0.2_dp, cp = 0.2_dp     ! B_theta r inside/outside
        real(dp) :: bzm = 1, pjump = 0           ! B_z inside; p+ - p-
        real(dp) :: delta = 1e-3_dp              ! outer wall amplitude
        real(dp) :: bzp = 1                      ! from pressure balance
        real(dp) :: coef(4) = 0                  ! A-, B-, A+, B+
        logical :: continuous = .false.
        integer :: status = -1
    contains
        procedure, non_overridable :: init => pinch_init
        procedure, non_overridable :: eval => pinch_eval
        procedure, non_overridable :: sheet => pinch_sheet_current
        procedure, non_overridable :: f => pinch_f
        procedure, non_overridable :: resonant_radius => pinch_resonant_radius
    end type
contains
    ! ---------------------------------------------------------------- slab
    subroutine slab_init(s)
        class(slab_sheet_t), intent(inout) :: s
        real(dp) :: lhs, rhs
        s%status = 1
        if (s%a <= 0 .or. s%k <= 0 .or. s%b0 <= 0 .or. s%dplus < 0 .or. &
            s%dminus < 0) return
        s%continuous = s%dplus == 0 .and. s%dminus == 0
        if (s%continuous) then
            s%xi0 = 0
        else
            call sheet_slab_interface(s%k, s%a, s%b0, s%dplus, s%dminus, s%wp, s%wm, &
                                      lhs, rhs)
            if (lhs == 0) return
            s%xi0 = rhs/lhs
        end if
        s%status = 0
    end subroutine

    pure real(dp) function slab_by(s, x, side) result(by)
        class(slab_sheet_t), intent(in) :: s
        real(dp), intent(in) :: x
        integer, intent(in) :: side   ! +1: x >= 0 branch, -1: x <= 0 branch
        if (side > 0) then
            by = s%b0*(x/s%a + s%dplus)
        else
            by = s%b0*(x/s%a - s%dminus)
        end if
    end function

    !> psi, dpsi/dx and, where F /= 0, xi = psi/F and dxi/dx on one side.
    !> At x = 0 in the continuous case xi is the one-sided limit (the step).
    pure subroutine slab_eval(s, x, side, psi, dpsi, xi, dxi)
        class(slab_sheet_t), intent(in) :: s
        real(dp), intent(in) :: x
        integer, intent(in) :: side
        real(dp), intent(out) :: psi, dpsi, xi, dxi
        real(dp) :: t, psi0, psia, dpsidt, f, fp, d2
        t = side*x
        if (side > 0) then
            psi0 = s%k*s%b0*s%dplus*s%xi0
            psia = s%k*s%b0*(1 + s%dplus)*s%wp
        else
            psi0 = -s%k*s%b0*s%dminus*s%xi0
            psia = -s%k*s%b0*(1 + s%dminus)*s%wm
        end if
        call sheet_slab_side(t, s%k, s%a, psi0, psia, psi, dpsidt)
        dpsi = side*dpsidt
        f = s%k*s%by(x, side)
        fp = s%k*s%b0/s%a
        if (f /= 0) then
            xi = psi/f
            dxi = (dpsi*f - psi*fp)/f**2
        else
            ! psi(0) = 0 with F(0) = 0: l'Hopital, psi'' = k^2 psi = 0 at 0.
            xi = dpsi/fp
            d2 = s%k**2*psi
            dxi = d2/(2*fp)
        end if
    end subroutine

    !> cos(ky) amplitude of the sheet-current perturbation K1_z = [[B_Ly]]
    !> = -[[psi']]/k at x = 0 (Lagrangian; Eulerian [[B1y]] when continuous).
    real(dp) function slab_sheet_current(s) result(k1)
        class(slab_sheet_t), intent(in) :: s
        real(dp) :: psi, dpp, dpm, xi, dxi
        call s%eval(0.0_dp, 1, psi, dpp, xi, dxi)
        call s%eval(0.0_dp, -1, psi, dpm, xi, dxi)
        k1 = -(dpp - dpm)/s%k
    end function

    ! --------------------------------------------------------- screw pinch
    pure real(dp) function pinch_f(s, r, side) result(f)
        class(pinch_sheet_t), intent(in) :: s
        real(dp), intent(in) :: r
        integer, intent(in) :: side
        if (side > 0) then
            f = s%m*s%cp/r**2 + s%k*s%bzp
        else
            f = s%m*s%cm/r**2 + s%k*s%bzm
        end if
    end function

    !> Radius where F = 0 for the inner (minus-side) profile.
    pure real(dp) function pinch_resonant_radius(s) result(r)
        class(pinch_sheet_t), intent(in) :: s
        r = -1
        if (s%k*s%bzm*s%m*s%cm < 0) r = sqrt(-s%m*s%cm/(s%k*s%bzm))
    end function

    subroutine pinch_init(s)
        class(pinch_sheet_t), intent(inout) :: s
        real(dp) :: kap, bz2, mat(4, 4), rhs(4), fm, fpl, cpr(2), cms(2)
        real(dp) :: ip, kp, i0, k0, d1(2), d2(2), v1(2), v2(2)
        integer :: st
        s%status = 1
        kap = abs(s%k)
        if (s%m < 1 .or. kap == 0 .or. s%r0 <= 0 .or. s%sigma <= s%r0 .or. &
            s%b <= s%sigma .or. kap*s%b > 2) return
        bz2 = s%bzm**2 + (s%cm**2 - s%cp**2)/s%sigma**2 - 2*s%pjump
        if (bz2 <= 0) return
        s%bzp = sign(sqrt(bz2), s%bzm)
        s%continuous = s%cp == s%cm .and. s%pjump == 0
        fm = s%f(s%sigma, -1); fpl = s%f(s%sigma, 1)
        if (s%continuous) then
            if (abs(fm) > 1e-12_dp*abs(s%m*s%cm/s%sigma**2)) return
        else
            if (fm == 0 .or. fpl == 0) return
        end if
        ! Columns: chi = A I_m(kap r) + B K_m(kap r); v = (chi, chi') per basis.
        call basis(s%r0, v1, d1)            ! d1 = (I', K') * kap at r0
        call basis(s%b, v2, d2)
        mat = 0; rhs = 0
        if (s%continuous) then
            ! chi- == 0; chi+'(sigma) = 0, chi+'(b) = -F(b) delta.
            call basis(s%sigma, cms, cpr)
            mat(1, 1:2) = cpr
            mat(2, 1:2) = d2
            rhs(2) = -s%f(s%b, 1)*s%delta
            call solve(mat(1:2, 1:2), rhs(1:2), st)
            if (st /= 0) return
            s%coef = [0.0_dp, 0.0_dp, rhs(1), rhs(2)]
        else
            call basis(s%sigma, cms, cpr)
            mat(1, 1:2) = d1
            mat(2, 3:4) = d2
            rhs(2) = -s%f(s%b, 1)*s%delta
            ! chi-'(sigma) F+ - chi+'(sigma) F- = 0
            mat(3, 1:2) = cpr*fpl
            mat(3, 3:4) = -cpr*fm
            ! F+ chi+ - F- chi- - xi_s (c+^2 - c-^2)/sigma^3 = 0,
            ! xi_s = -chi+'(sigma)/F+.
            mat(4, 1:2) = -fm*cms
            mat(4, 3:4) = fpl*cms + cpr/fpl*(s%cp**2 - s%cm**2)/s%sigma**3
            call solve(mat, rhs, st)
            if (st /= 0) return
            s%coef = rhs
        end if
        s%status = 0
    contains
        subroutine basis(r, v, d)
            real(dp), intent(in) :: r
            real(dp), intent(out) :: v(2), d(2)
            call modified_bessel_ik(s%m, kap*r, i0, ip, k0, kp)
            v = [i0, k0]
            d = kap*[ip, kp]
        end subroutine
    end subroutine

    !> chi-hat, chi-hat' and xi-hat = -chi'/F on one side (side = -1: r <=
    !> sigma, +1: r >= sigma). Continuous case at r = sigma: one-sided limit.
    subroutine pinch_eval(s, r, side, chi, dchi, xi)
        class(pinch_sheet_t), intent(in) :: s
        real(dp), intent(in) :: r
        integer, intent(in) :: side
        real(dp), intent(out) :: chi, dchi, xi
        real(dp) :: i0, ip, k0, kp, kap, a, b, f, fp, d2
        kap = abs(s%k)
        if (side > 0) then
            a = s%coef(3); b = s%coef(4)
        else
            a = s%coef(1); b = s%coef(2)
        end if
        call modified_bessel_ik(s%m, kap*r, i0, ip, k0, kp)
        chi = a*i0 + b*k0
        dchi = kap*(a*ip + b*kp)
        f = s%f(r, side)
        if (f /= 0) then
            xi = -dchi/f
        else
            fp = -2*s%m*merge(s%cp, s%cm, side > 0)/r**3
            d2 = -dchi/r + (s%m**2/r**2 + s%k**2)*chi
            xi = -d2/fp
        end if
    end subroutine

    !> cos-phase amplitudes of [[B_L theta]] and [[B_L z]] at r = sigma
    !> (K = r x [[B_L]]); Eulerian [[B1]] in the continuous case.
    subroutine pinch_sheet_current(s, jtheta, jz)
        class(pinch_sheet_t), intent(in) :: s
        real(dp), intent(out) :: jtheta, jz
        real(dp) :: chim, dchim, xim, chip, dchip, xip
        call s%eval(s%sigma, -1, chim, dchim, xim)
        call s%eval(s%sigma, 1, chip, dchip, xip)
        jz = s%k*(chip - chim)
        if (s%continuous) then
            jtheta = s%m*(chip - chim)/s%sigma
        else
            jtheta = s%m*(chip - chim)/s%sigma &
                - xip*(s%cp - s%cm)/s%sigma**2
        end if
    end subroutine

    !> I_m, I_m', K_m, K_m' for integer m >= 0 and 0 < x <= 2: power series
    !> (A&S 9.6.10, 9.6.13, 9.6.11) and upward recurrence for K.
    pure subroutine modified_bessel_ik(m, x, bi, dbi, bk, dbk)
        integer, intent(in) :: m
        real(dp), intent(in) :: x
        real(dp), intent(out) :: bi, dbi, bk, dbk
        real(dp), parameter :: euler = 0.57721566490153286061_dp
        real(dp) :: i0, i1, k0, k1, km, kn, term, h, y, lg, ip1
        integer :: j, n
        bi = iseries(m); ip1 = iseries(m + 1)
        y = x*x/4
        lg = log(x/2) + euler
        ! K0 = -(ln(x/2)+gamma) I0 + sum y^j/(j!)^2 H_j
        i0 = iseries(0); i1 = iseries(1)
        k0 = -lg*i0; term = 1; h = 0
        do j = 1, 60
            term = term*y/real(j, dp)**2
            h = h + 1.0_dp/j
            k0 = k0 + term*h
            if (term*h < 1e-18_dp*abs(k0)) exit
        end do
        ! K1 = 1/x + ln(x/2) I1 - (x/4) sum_j (psi(j+1)+psi(j+2)) y^j/(j!(j+1)!)
        !    with psi(n) = -gamma + H_(n-1).
        k1 = 1/x + log(x/2)*i1
        term = 1; h = 0
        k1 = k1 - x/4*(term*(-2*euler + 1))
        do j = 1, 60
            term = term*y/(real(j, dp)*real(j + 1, dp))
            h = h + 1.0_dp/j
            k1 = k1 - x/4*term*(-2*euler + 2*h + 1.0_dp/(j + 1))
            if (term < 1e-18_dp) exit
        end do
        if (m == 0) then
            bk = k0; dbk = -k1
            dbi = i1
            return
        end if
        km = k0; kn = k1
        do n = 1, m - 1
            term = km + 2*n/x*kn
            km = kn; kn = term
        end do
        bk = kn
        ! K_m' = -K_(m-1) - (m/x) K_m, I_m' = I_(m+1) + (m/x) I_m
        dbk = -km - m/x*kn
        dbi = ip1 + m/x*bi
    contains
        pure real(dp) function iseries(nu) result(v)
            integer, intent(in) :: nu
            real(dp) :: t
            integer :: jj
            t = (x/2)**nu
            do jj = 1, nu
                t = t/jj
            end do
            v = t
            do jj = 1, 80
                t = t*(x*x/4)/(real(jj, dp)*real(jj + nu, dp))
                v = v + t
                if (t < 1e-18_dp*v) exit
            end do
        end function
    end subroutine

    !> Gaussian elimination with partial pivoting; b is overwritten.
    pure subroutine solve(a, b, status)
        real(dp), intent(inout) :: a(:, :), b(:)
        integer, intent(out) :: status
        integer :: n, i, j, p
        real(dp) :: f, row(size(b)), t
        n = size(b); status = 1
        do i = 1, n
            p = maxloc(abs(a(i:n, i)), 1) + i - 1
            if (a(p, i) == 0) return
            if (p /= i) then
                row = a(i, :); a(i, :) = a(p, :); a(p, :) = row
                t = b(i); b(i) = b(p); b(p) = t
            end if
            do j = i + 1, n
                f = a(j, i)/a(i, i)
                a(j, i:n) = a(j, i:n) - f*a(i, i:n)
                b(j) = b(j) - f*b(i)
            end do
        end do
        do i = n, 1, -1
            b(i) = (b(i) - sum(a(i, i + 1:n)*b(i + 1:n)))/a(i, i)
        end do
        status = 0
    end subroutine
end module kin6d_current_sheet_oracles
