module test_geometry_fields
    ! Analytic axisymmetric fields, B = grad phi x grad psi + F/R phi_hat.
    !  solovev: psi = c ((R^2 - R0^2)^2 + 4 R^2 Z^2) - 4 c R0^2 a^2, F constant.
    !  circle:  concentric circles about (R0, 0), psi'(r) = F r/(q_t(r) R0) with
    !           q_t = q0 + q2 r^2, whence q(r) = q_t R0/sqrt(R0^2 - r^2), the
    !           volume V = 2 pi^2 R0 r^2 and |grad psi| constant on a surface.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_state, only: equilibrium_field_t, EQ_BR, EQ_BZ, &
        EQ_BPHI, EQ_P, EQ_PSI, EQ_PSI_R, EQ_PSI_Z
    implicit none
    private
    public :: solovev_t, circle_t, r0, minor, f_tor, c_sol, solovev_b
    public :: circ_r0, circ_f, circ_q0, circ_q2, circ_psi, circ_dpsi, circ_q
    real(dp), parameter :: r0 = 6.2_dp, minor = 0.62_dp, f_tor = 32.86_dp
    real(dp), parameter :: c_sol = f_tor/(12*r0**3)
    real(dp), parameter :: circ_r0 = 1, circ_f = 1, circ_q0 = 1.5_dp, circ_q2 = 3
    type, extends(equilibrium_field_t) :: solovev_t
    contains
        procedure :: sample => solovev_sample
    end type
    type, extends(equilibrium_field_t) :: circle_t
    contains
        procedure :: sample => circle_sample
    end type
contains
    subroutine solovev_sample(field, points, values, statuses)
        class(solovev_t), intent(in) :: field
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        real(dp) :: r, z, psi_r, psi_z
        integer :: i
        values = 0; statuses = 0
        do i = 1, size(points, 2)
            r = points(1, i); z = points(2, i)
            if (r <= 0) then
                statuses(i) = 1
                cycle
            end if
            psi_r = c_sol*(4*r*(r*r - r0*r0) + 8*r*z*z)
            psi_z = c_sol*8*r*r*z
            values(EQ_PSI, i) = c_sol*((r*r - r0*r0)**2 + 4*r*r*z*z &
                - 4*r0*r0*minor*minor)
            values(EQ_PSI_R, i) = psi_r; values(EQ_PSI_Z, i) = psi_z
            values(EQ_BR, i) = -psi_z/r; values(EQ_BZ, i) = psi_r/r
            values(EQ_BPHI, i) = f_tor/r
            values(EQ_P, i) = -16*c_sol*values(EQ_PSI, i)
        end do
    end subroutine

    ! |B| and its (R, Z) gradient in closed form.
    subroutine solovev_b(r, z, bmag, grad)
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: bmag, grad(2)
        real(dp) :: psi_r, psi_z, psi_rr, psi_rz, psi_zz, n, n_r, n_z
        psi_r = c_sol*(4*r*(r*r - r0*r0) + 8*r*z*z); psi_z = c_sol*8*r*r*z
        psi_rr = c_sol*(12*r*r - 4*r0*r0 + 8*z*z)
        psi_rz = 16*c_sol*r*z; psi_zz = 8*c_sol*r*r
        n = psi_r**2 + psi_z**2 + f_tor**2
        n_r = 2*psi_r*psi_rr + 2*psi_z*psi_rz
        n_z = 2*psi_r*psi_rz + 2*psi_z*psi_zz
        bmag = sqrt(n)/r
        grad = [n_r/(2*sqrt(n))/r - sqrt(n)/r**2, n_z/(2*sqrt(n))/r]
    end subroutine

    pure real(dp) function circ_psi(r)
        real(dp), intent(in) :: r
        circ_psi = circ_f/circ_r0*log(circ_q0 + circ_q2*r*r)/(2*circ_q2)
    end function

    pure real(dp) function circ_dpsi(r)
        real(dp), intent(in) :: r
        circ_dpsi = circ_f*r/((circ_q0 + circ_q2*r*r)*circ_r0)
    end function

    pure real(dp) function circ_q(r)
        real(dp), intent(in) :: r
        circ_q = (circ_q0 + circ_q2*r*r)*circ_r0/sqrt(circ_r0**2 - r*r)
    end function

    subroutine circle_sample(field, points, values, statuses)
        class(circle_t), intent(in) :: field
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        real(dp) :: r, z, rho2, qt, psi_r, psi_z
        integer :: i
        values = 0; statuses = 0
        do i = 1, size(points, 2)
            r = points(1, i); z = points(2, i)
            if (r <= 0) then
                statuses(i) = 1
                cycle
            end if
            rho2 = (r - circ_r0)**2 + z*z
            qt = circ_q0 + circ_q2*rho2
            psi_r = circ_f*(r - circ_r0)/(qt*circ_r0)
            psi_z = circ_f*z/(qt*circ_r0)
            values(EQ_PSI, i) = circ_psi(sqrt(rho2))
            values(EQ_PSI_R, i) = psi_r; values(EQ_PSI_Z, i) = psi_z
            values(EQ_BR, i) = -psi_z/r; values(EQ_BZ, i) = psi_r/r
            values(EQ_BPHI, i) = circ_f/r
        end do
    end subroutine
end module test_geometry_fields

program test_equilibrium_geometry
    ! Kinetic geometry view of D7 (src/equilibrium/equilibrium_geometry.f90).
    !   solovev   closed-form B, |B|, grad|B|; q, q', dV/dpsi and flux-surface
    !             averages against independent quadrature on the exact surface
    !             parametrisation; b.grad theta* = B_phi/(q R B); Boozer
    !             Jacobian condition; field-line tangent along B
    !   circle    exact q, q', V', <R>, <1/R^2>, nu, d nu/d psi, PEST = Boozer
    !   drift     bounce-averaged grad-B drift of the geometry view against the
    !             large-aspect eps-expansion (Hinton-Hazeltine form, deeply
    !             trapped limit = omega_gradB_scale of ntv/analytic/regimes.py):
    !             -(q w/(Z e B0 r R0)) [2E/K - 1 + 4 s (E/K + k^2 - 1)], error
    !             O(eps); the minus sign is the convention alpha = phi - nu with
    !             B_phi, q > 0 (b x grad B points up on the outboard side)
    !   state     view of a psi-lane P3 Solov'ev state against the D7 readback
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre_ab
    use kin6d_equilibrium_geometry
    use kin6d_equilibrium_state, only: equilibrium_state_t
    use kin6d_equilibrium_readback, only: equilibrium_readback_t, read_equilibrium
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs
    use kin6d_gs_equilibrium, only: gs_equilibrium
    use test_geometry_fields
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), two_pi = 2*pi
    character(32) :: mode
    call get_command_argument(1, mode)
    select case (trim(mode))
    case ('solovev'); call test_solovev()
    case ('circle'); call test_circle()
    case ('drift'); call test_drift()
    case ('state'); call test_state()
    case default; error stop 'mode: solovev | circle | drift | state'
    end select
contains
    subroutine check(label, value, tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: value, tolerance
        print '(a34,es10.3,a,es9.2)', label, value, '  <=', tolerance
        if (.not. value <= tolerance) error stop 'geometry check failed'
    end subroutine

    ! ---- Solov'ev ---------------------------------------------------------

    subroutine solovev_surface(s, r, z, dl, gp)
        ! Exact surface R^2 = R0^2 + 2 R0 s cos t, R Z = R0 s sin t.
        real(dp), intent(in) :: s
        real(dp), intent(out) :: r(:), z(:), dl(:), gp(:)
        real(dp) :: t, rt, zt, pr, pz
        integer :: i, n
        n = size(r)
        do i = 1, n
            t = two_pi*(i - 1)/n
            r(i) = sqrt(r0*r0 + 2*r0*s*cos(t))
            rt = -r0*s*sin(t)/r(i)
            z(i) = r0*s*sin(t)/r(i)
            zt = r0*s*cos(t)/r(i) - r0*s*sin(t)*rt/r(i)**2
            dl(i) = hypot(rt, zt)
            pr = c_sol*(4*r(i)*(r(i)**2 - r0**2) + 8*r(i)*z(i)**2)
            pz = c_sol*8*r(i)**2*z(i)
            gp(i) = hypot(pr, pz)
        end do
    end subroutine

    subroutine solovev_reference(s, q, vp, avg_r2, avg_b2)
        real(dp), intent(in) :: s
        real(dp), intent(out) :: q, vp, avg_r2, avg_b2
        integer, parameter :: n = 4000
        real(dp) :: r(n), z(n), dl(n), gp(n), w(n)
        call solovev_surface(s, r, z, dl, gp)
        q = f_tor/two_pi*sum(dl/(r*gp))*two_pi/n
        w = r*dl/gp
        vp = two_pi*sum(w)*two_pi/n
        avg_r2 = sum(w*r*r)/sum(w)
        avg_b2 = sum(w*(gp**2 + f_tor**2)/r**2)/sum(w)
    end subroutine

    subroutine test_solovev()
        type(solovev_t) :: field
        type(kinetic_surface_t) :: surf
        type(kinetic_samples_t) :: sm, line
        real(dp), parameter :: s = 0.4_dp, h = 1.0e-3_dp
        real(dp) :: psi0, q, vp, avg_r2, avg_b2, qp, qm, qp2, qm2, dpsi_ds, bm, gb(2)
        real(dp) :: err_b, err_g, err_i, err_j, err_t, err_bz, tgrid(37), dth, tv(3)
        real(dp) :: dq, ratio_err
        integer :: status, i, n
        psi0 = 4*c_sol*r0*r0*(s*s - minor*minor)
        call build_kinetic_surface(field, [r0, 0.0_dp], psi0, surf, status, 128)
        if (status /= 0) error stop 'build'
        call solovev_reference(s, q, vp, avg_r2, avg_b2)
        call check('q', abs(surf%q/q - 1), 1.0e-11_dp)
        call check('dV/dpsi', abs(surf%dv_dpsi/vp - 1), 1.0e-11_dp)
        call solovev_reference(s + h, qp, vp, avg_r2, avg_b2)
        call solovev_reference(s - h, qm, vp, avg_r2, avg_b2)
        call solovev_reference(s + 2*h, qp2, vp, avg_r2, avg_b2)
        call solovev_reference(s - 2*h, qm2, vp, avg_r2, avg_b2)
        dpsi_ds = 8*c_sol*r0*r0*s
        dq = (-qp2 + 8*qp - 8*qm + qm2)/(12*h)/dpsi_ds
        call check('dq/dpsi', abs(surf%dq_dpsi/dq - 1), 1.0e-6_dp)
        call solovev_reference(s, q, vp, avg_r2, avg_b2)
        call check('<1>', abs(surface_average(surf, [(1.0_dp, i=1, 128)]) - 1), &
            1.0e-14_dp)
        call check('<R^2>', abs(surface_average(surf, surf%node%r**2)/avg_r2 - 1), &
            1.0e-11_dp)
        call check('<B^2>', abs(surface_average(surf, surf%node%bmag**2)/avg_b2 - 1), &
            1.0e-11_dp)
        call check('<R^2> (function)', &
            abs(surface_average_function(surf, rsq)/avg_r2 - 1), 1.0e-11_dp)
        call check('<B.grad phi> = 4pi^2 q/V''', abs(surface_average(surf, &
            surf%node%b(3, :)/surf%node%r)/(4*pi*pi*surf%q/surf%dv_dpsi) - 1), &
            1.0e-11_dp)
        ! |B| and grad|B| at arbitrary geometric and PEST angles.
        tgrid = [(-pi + 0.17_dp + two_pi*(i - 1)/37, i=1, 37)]
        err_b = 0; err_g = 0; err_i = 0; err_j = 0; err_bz = 0
        do n = KG_ANGLE_GEOMETRIC, KG_ANGLE_BOOZER
            call sample_kinetic_geometry(field, surf, n, tgrid, sm, status)
            if (status /= 0) error stop 'sample'
            do i = 1, size(tgrid)
                call solovev_b(sm%r(i), sm%z(i), bm, gb)
                err_b = max(err_b, abs(sm%bmag(i)/bm - 1))
                err_g = max(err_g, norm2(sm%grad_b(1:2, i) - gb)/norm2(gb))
                ! b.grad theta* = B_phi/(q R B) from the contour differential.
                if (n == KG_ANGLE_PEST) err_i = max(err_i, abs(sm%b_grad_theta(i) &
                    /(sm%b(3, i)/(surf%q*sm%r(i)*sm%bmag(i))) - 1))
                ! Boozer: d theta_B/d theta* = (h/<h>)/(g/q), h = R B^2 dl/|grad psi|
                ! tested through the independent angle map below.
                err_bz = max(err_bz, abs(sm%grad_psi_mag(i) - hypot(sm%grad_psi(1, i), &
                    sm%grad_psi(2, i))))
            end do
        end do
        call check('|B|', err_b, 1.0e-13_dp)
        call check('grad|B|', err_g, 1.0e-9_dp)
        call check('b.grad theta*', err_i, 1.0e-10_dp)
        ! Boozer condition: the Jacobian 1/(grad psi x grad theta_B . grad zeta_B)
        ! is proportional to 1/B^2, i.e. d theta_B/d theta_g * |grad psi|/(R B^2 dl)
        ! is constant; with the geometry at the nodes, from finite differences of
        ! the returned angles over theta_g (fine grid, second-order).
        block
            integer, parameter :: m = 2001
            real(dp) :: tg(m), dtheta, ref(m), dt_b(m), dt_s(m)
            type(kinetic_samples_t) :: gm
            tg = [(two_pi*(i - 1)/(m - 1), i=1, m)]
            dtheta = two_pi/(m - 1)
            call sample_kinetic_geometry(field, surf, KG_ANGLE_GEOMETRIC, tg, gm, &
                status)
            if (status /= 0) error stop 'boozer sample'
            ! d theta_B/d theta_g and d theta*/d theta_g (central differences).
            do i = 2, m - 1
                dt_b(i) = (gm%theta_boozer(i + 1) - gm%theta_boozer(i - 1))/(2*dtheta)
                dt_s(i) = (gm%theta_pest(i + 1) - gm%theta_pest(i - 1))/(2*dtheta)
                ! |dx/dtheta| from the tangent of the sampled surface.
                ref(i) = hypot(gm%r(i + 1) - gm%r(i - 1), gm%z(i + 1) - gm%z(i - 1)) &
                    /(2*dtheta)
            end do
            ! J B^2 ~ R B^2 dl/(|grad psi| dtheta_B) constant over the surface
            err_j = 0
            ratio_err = 0
            do i = 2, m - 1
                ref(i) = gm%r(i)*gm%bmag(i)**2*ref(i)/gm%grad_psi_mag(i)/dt_b(i)
            end do
            ratio_err = (maxval(ref(2:m - 1)) - minval(ref(2:m - 1))) &
                /sum(ref(2:m - 1))*(m - 2)
            call check('Boozer J B^2 constant', ratio_err, 2.0e-5_dp)
            ! PEST: R dl/(|grad psi| dtheta*) = q-normalised, B_phi dl/(|grad psi|
            ! dtheta*) = q is the straight-field-line condition.
            err_t = 0
            do i = 2, m - 1
                err_t = max(err_t, abs(gm%b(3, i)*hypot(gm%r(i + 1) - gm%r(i - 1), &
                    gm%z(i + 1) - gm%z(i - 1))/(2*dtheta)/gm%grad_psi_mag(i) &
                    /dt_s(i)/surf%q - 1))
            end do
            call check('PEST B_phi dl/|grad psi| = q', err_t, 2.0e-5_dp)
            call check('theta_B(2pi)', abs(gm%theta_boozer(m) - two_pi), 1.0e-12_dp)
        end block
        ! Field line: positions at (R, Z, phi = alpha_0 + q theta*) follow B.
        dth = 1.0e-4_dp
        tgrid(1:3) = [1.3_dp, 1.3_dp + dth, 1.3_dp - dth]
        call sample_kinetic_geometry(field, surf, KG_ANGLE_PEST, tgrid(1:3), line, &
            status)
        if (status /= 0) error stop 'line'
        tv(1) = line%r(2) - line%r(3)                    ! dR
        tv(2) = line%z(2) - line%z(3)                    ! dZ
        tv(3) = line%r(1)*(line%nu(2) - line%nu(3))      ! R dphi
        call check('field-line tangent x B', norm2([ &
            tv(2)*line%b(3, 1) - tv(3)*line%b(2, 1), &
            tv(3)*line%b(1, 1) - tv(1)*line%b(3, 1), &
            tv(1)*line%b(2, 1) - tv(2)*line%b(1, 1)])/(norm2(tv)*line%bmag(1)), &
            1.0e-7_dp)
        ! Unwrapped angles: nu advances 2 pi q per turn, geometry is periodic.
        tgrid(1:2) = [0.7_dp, 0.7_dp + 3*two_pi]
        call sample_kinetic_geometry(field, surf, KG_ANGLE_PEST, tgrid(1:2), line, &
            status)
        call check('nu per turn', abs((line%nu(2) - line%nu(1))/(3*two_pi*surf%q) &
            - 1), 1.0e-13_dp)
        call check('periodic R', abs(line%r(2) - line%r(1)), 1.0e-12_dp)
        print '(a)', 'PASS: kinetic geometry, Solov''ev'
    end subroutine

    real(dp) function rsq(r, z)
        real(dp), intent(in) :: r, z
        rsq = r*r + 0*z
    end function

    ! ---- Large-aspect circle ---------------------------------------------

    real(dp) function nu_exact(a, theta)
        real(dp), intent(in) :: a, theta
        nu_exact = (circ_q0 + circ_q2*a*a)*circ_r0*2/sqrt(circ_r0**2 - a*a) &
            *atan(sqrt((circ_r0 - a)/(circ_r0 + a))*tan(theta/2))
    end function

    subroutine test_circle()
        type(circle_t) :: field
        type(kinetic_surface_t) :: surf
        type(kinetic_samples_t) :: sm
        real(dp), parameter :: a = 0.2_dp, h = 1.0e-3_dp
        real(dp) :: qa, dqdr, tgrid(9), ex, nup, tg_ex, e1, e2, e3, e4
        integer :: status, i
        call build_kinetic_surface(field, [circ_r0, 0.0_dp], circ_psi(a), surf, &
            status, 96)
        if (status /= 0) error stop 'build'
        qa = circ_q(a)
        dqdr = (-circ_q(a + 2*h) + 8*circ_q(a + h) - 8*circ_q(a - h) &
            + circ_q(a - 2*h))/(12*h)
        call check('q', abs(surf%q/qa - 1), 1.0e-12_dp)
        call check('dq/dpsi', abs(surf%dq_dpsi/(dqdr/circ_dpsi(a)) - 1), 1.0e-7_dp)
        call check('dV/dpsi', abs(surf%dv_dpsi/(4*pi*pi*circ_r0*a/circ_dpsi(a)) &
            - 1), 1.0e-12_dp)
        call check('<R>', abs(surface_average(surf, surf%node%r) &
            /((circ_r0**2 + a*a/2)/circ_r0) - 1), 1.0e-13_dp)
        call check('<1/R^2>', abs(surface_average(surf, 1/surf%node%r**2) &
            *circ_r0*sqrt(circ_r0**2 - a*a) - 1), 1.0e-13_dp)
        tgrid = [(0.3_dp*i, i=1, 9)]
        tgrid(9) = pi - 0.2_dp
        call sample_kinetic_geometry(field, surf, KG_ANGLE_GEOMETRIC, tgrid, sm, &
            status)
        if (status /= 0) error stop 'sample'
        e1 = 0; e2 = 0; e3 = 0; e4 = 0
        do i = 1, 9
            ex = nu_exact(a, tgrid(i))
            e1 = max(e1, abs(sm%nu(i)/ex - 1))
            ! d nu/d psi at fixed theta_g by differences of the closed form.
            nup = (-nu_exact(a + 2*h, tgrid(i)) + 8*nu_exact(a + h, tgrid(i)) &
                - 8*nu_exact(a - h, tgrid(i)) + nu_exact(a - 2*h, tgrid(i))) &
                /(12*h)/circ_dpsi(a)
            e2 = max(e2, abs(sm%dalpha_dpsi(i) + nup)/abs(nup))
            ! PEST and Boozer coincide: g and h are both ~ 1/R here.
            e3 = max(e3, abs(sm%theta_boozer(i) - sm%theta_pest(i)))
            e4 = max(e4, abs(sm%theta_pest(i) - 2*atan(sqrt((circ_r0 - a)/ &
                (circ_r0 + a))*tan(tgrid(i)/2))))
        end do
        call check('nu(theta_g)', e1, 1.0e-12_dp)
        call check('d alpha/d psi', e2, 1.0e-7_dp)
        call check('theta_B - theta*', e3, 1.0e-12_dp)
        call check('theta*', e4, 1.0e-12_dp)
        ! Inverse map: PEST angle of a geometric one and back.
        tgrid(1:3) = [0.4_dp, 2.0_dp, 3.0_dp]
        do i = 1, 3
            tg_ex = 2*atan(sqrt((circ_r0 + a)/(circ_r0 - a))*tan(tgrid(i)/2))
            call sample_kinetic_geometry(field, surf, KG_ANGLE_PEST, tgrid(i:i), sm, &
                status)
            e1 = abs(sm%theta_g(1) - tg_ex)
            call check('theta_g(theta*)', e1, 1.0e-12_dp)
        end do
        print '(a)', 'PASS: kinetic geometry, circle'
    end subroutine

    ! ---- Bounce-averaged grad-B drift ---------------------------------------

    subroutine ellip(m, kk, ee)
        real(dp), intent(in) :: m
        real(dp), intent(out) :: kk, ee
        real(dp) :: x(64), w(64)
        call gauss_legendre_ab(64, 0.0_dp, pi/2, x, w)
        kk = sum(w/sqrt(1 - m*sin(x)**2))
        ee = sum(w*sqrt(1 - m*sin(x)**2))
    end subroutine

    real(dp) function bounce_drift(field, surf, kappa2) result(omega)
        ! <(2 - lambda B) drift_alpha>_bounce for lambda from kappa^2 in units of
        ! w/(Z e); theta = theta_b sin(xi) removes the turning-point singularity.
        type(circle_t), intent(in) :: field
        type(kinetic_surface_t), intent(in) :: surf
        real(dp), intent(in) :: kappa2
        type(kinetic_samples_t) :: s
        integer, parameter :: m = 48
        real(dp) :: x(m), w(m), th(m), lam, bmin, bmax, lo, hi, mid, vpar
        real(dp) :: num, den, wt
        integer :: status, i
        call sample_kinetic_geometry(field, surf, KG_ANGLE_GEOMETRIC, &
            [0.0_dp, pi], s, status)
        if (status /= 0) error stop 'bounce'
        bmin = s%bmag(1); bmax = s%bmag(2)
        lam = 1/(bmin + kappa2*(bmax - bmin))
        lo = 0; hi = pi
        do i = 1, 60
            mid = (lo + hi)/2
            call sample_kinetic_geometry(field, surf, KG_ANGLE_GEOMETRIC, [mid], s, &
                status)
            if (1 - lam*s%bmag(1) > 0) then
                lo = mid
            else
                hi = mid
            end if
        end do
        call gauss_legendre_ab(m, -pi/2, pi/2, x, w)
        th = mid*sin(x)
        call sample_kinetic_geometry(field, surf, KG_ANGLE_GEOMETRIC, th, s, status)
        if (status /= 0) error stop 'bounce samples'
        num = 0; den = 0
        do i = 1, m
            vpar = sqrt(max(1 - lam*s%bmag(i), tiny(1.0_dp)))
            wt = w(i)*mid*cos(x(i))/(vpar*s%b_grad_theta(i))
            num = num + wt*(2 - lam*s%bmag(i))*s%drift_alpha(i)
            den = den + wt
        end do
        omega = num/den
    end function

    subroutine test_drift()
        type(circle_t) :: field
        type(kinetic_surface_t) :: surf
        real(dp), parameter :: kappa2(3) = [1.0e-3_dp, 0.3_dp, 0.7_dp]
        real(dp), parameter :: eps(3) = [0.04_dp, 0.02_dp, 0.01_dp]
        real(dp) :: a, kk, ee, ex, b0, shear, err(3, 3), omega, dqdr
        integer :: status, i, j
        do i = 1, 3
            a = eps(i)*circ_r0
            call build_kinetic_surface(field, [circ_r0, 0.0_dp], circ_psi(a), surf, &
                status, 64)
            if (status /= 0) error stop 'build'
            b0 = circ_f/circ_r0
            dqdr = (-circ_q(a + 2e-4_dp) + 8*circ_q(a + 1e-4_dp) &
                - 8*circ_q(a - 1e-4_dp) + circ_q(a - 2e-4_dp))/(12e-4_dp)
            shear = a*dqdr/circ_q(a)
            do j = 1, 3
                call ellip(kappa2(j), kk, ee)
                ex = circ_q(a)/(b0*a*circ_r0)*(2*ee/kk - 1 &
                    + 4*shear*(ee/kk + kappa2(j) - 1))
                omega = bounce_drift(field, surf, kappa2(j))
                err(i, j) = abs(omega/ex + 1)
                print '(a,f6.3,a,f6.3,a,es11.3,a,es11.3,a,es10.3)', ' eps', eps(i), &
                    ' k2', kappa2(j), ' omega', omega, ' analytic', ex, ' rel', &
                    err(i, j)
            end do
        end do
        ! O(eps): the error halves with eps and is small at the smallest eps.
        do j = 1, 3
            call check('error ~ eps', err(3, j)/err(2, j), 0.75_dp)
            call check('error ~ eps (2)', err(2, j)/err(1, j), 0.75_dp)
            call check('error at eps = 0.01', err(3, j), 0.1_dp)
        end do
        print '(a)', 'PASS: bounce-averaged grad-B drift'
    end subroutine

    ! ---- psi-lane state -------------------------------------------------------

    subroutine test_state()
        type(gs_state_t) :: mesh
        type(equilibrium_state_t) :: eq
        type(equilibrium_readback_t) :: rb
        type(kinetic_surface_t) :: surf
        real(dp) :: contour(2, 48), hh
        integer :: j, status, k
        do j = 1, 48
            call solovev_boundary(real(j - 1, dp)/48, contour(:, j))
        end do
        hh = two_pi*minor/48
        call initialize_gs(mesh, contour, hh*hh/2, status, 3, curve=solovev_boundary)
        if (status /= 0) error stop 'mesh'
        call solve_gs(mesh, 1.0_dp, -16*c_sol, 0.0_dp, status)
        if (status /= 0) error stop 'solve'
        call gs_equilibrium(mesh, eq, status, f_tor, 0.0_dp, pprime=-16*c_sol)
        if (status /= 0) error stop 'adapter'
        call read_equilibrium(eq, 65, 0, rb, status)
        if (status /= 0) error stop 'readback'
        do k = 25, 49, 12
            call build_kinetic_surface_from_state(eq, real(k - 1, dp)/64, surf, &
                status, 64)
            if (status /= 0) error stop 'state surface'
            call check('q vs D7 readback', abs(surf%q/rb%q(k) - 1), 1.0e-5_dp)
            call check('<1>', abs(surface_average(surf, [(1.0_dp, j=1, 64)]) - 1), &
                1.0e-14_dp)
            call check('<B.grad phi> = 4pi^2 q/V''', abs(surface_average(surf, &
                surf%node%b(3, :)/surf%node%r)/(4*pi*pi*surf%q/surf%dv_dpsi) - 1), &
                1.0e-8_dp)
        end do
        print '(a)', 'PASS: kinetic geometry, psi-lane state'
    end subroutine

    subroutine solovev_boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: theta
        theta = two_pi*t
        point(1) = sqrt(r0*r0 + 2*r0*minor*cos(theta))
        point(2) = r0*minor*sin(theta)/point(1)
    end subroutine
end program test_equilibrium_geometry
