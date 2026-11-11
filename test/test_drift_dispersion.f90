!> Fluid drift-wave roots against independent oracles: the eigenvalues of the
!> linearised Hasegawa-Wakatani equations in binary128, Eq. (8) of Numata,
!> Ball and Dewar, Phys. Plasmas 14, 102312 (2007), the adiabatic asymptote,
!> the drift and Alfven limits of the drift-Alfven cubic, its
!> Hasegawa-Wakatani limit, and tangents against central differences.
program test_drift_dispersion
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use milca_drift_dispersion
    use milca_drift_fluid_leaves, only: drift_alfven_limits
    implicit none
    integer :: failures

    failures = 0
    call hasegawa_wakatani_roots()
    call hasegawa_wakatani_special_modes()
    call hasegawa_wakatani_tangent()
    call drift_alfven()
    if (failures > 0) error stop 'test_drift_dispersion failed'
    write (*, '(a)') 'test_drift_dispersion passed'
contains
    subroutine hasegawa_wakatani_roots()
        type(hw_parameters_t) :: par
        real(dp), parameter :: alphas(5) = [1e-3_dp, 0.1_dp, 1.0_dp, 30.0_dp, 1e4_dp]
        real(dp), parameter :: kxs(3) = [0.0_dp, 0.15_dp, 1.5_dp]
        real(dp), parameter :: kys(3) = [0.15_dp, 0.6_dp, 1.5_dp]
        complex(dp) :: omega(2), response(2), eq8, asymptote
        complex(qp) :: ref(2), wrong(2)
        real(dp) :: worst, worst_eq8, worst_vec, mutant, k2, b, wstar, kx, ky
        real(dp) :: gmax, gaxis
        integer :: ia, ix, iy, j, status

        worst = 0; worst_eq8 = 0; worst_vec = 0; mutant = 0
        par%kappa = 0.7_dp; par%d = 1e-4_dp; par%order = 2
        do ia = 1, size(alphas)
            par%alpha = alphas(ia)
            do ix = 1, size(kxs)
                do iy = 1, size(kys)
                    kx = kxs(ix); ky = kys(iy)
                    call hw_linear_roots(par, kx, ky, omega, response, status)
                    if (status /= drift_root_ok) failures = failures + 1
                    call matrix_roots(par, kx, ky, .false., ref)
                    call matrix_roots(par, kx, ky, .true., wrong)
                    k2 = kx**2 + ky**2
                    b = par%alpha*(1 + k2)/k2
                    wstar = ky*par%kappa/(1 + k2)
                    do j = 1, 2
                        worst = max(worst, real(abs(omega(j) - ref(j))/abs(ref(j)), dp))
                        eq8 = omega(j)**2 + (0, 1)*omega(j)*(b + 2*par%d*k2**2) &
                              - (0, 1)*b*wstar - par%alpha*par%d*k2*(k2 + 1) &
                              - par%d**2*k2**4
                        worst_eq8 = max(worst_eq8, abs(eq8)/(abs(omega(j))**2 + b*abs(wstar)))
                        worst_vec = max(worst_vec, eigen_defect(par, kx, ky, &
                                                                omega(j), response(j)))
                    end do
                    mutant = max(mutant, real(abs(omega(1) - wrong(1))/abs(ref(1)), dp))
                end do
            end do
        end do
        call check('HW roots vs binary128 matrix eigenvalues', worst, 1e-13_dp)
        call check('HW roots satisfy Numata Eq. (8)', worst_eq8, 1e-13_dp)
        call check('HW eigenvector n/phi', worst_vec, 1e-12_dp)
        ! Mutant oracle: density-gradient drive in the vorticity equation.
        call check_fails('mutant drive in vorticity equation is rejected', mutant, 1e-3_dp)

        ! Adiabatic asymptote for D = 0.
        par%alpha = 1e4_dp; par%d = 0
        kx = 0.3_dp; ky = 0.5_dp
        k2 = kx**2 + ky**2
        wstar = ky*par%kappa/(1 + k2)
        call hw_linear_roots(par, kx, ky, omega, response, status)
        asymptote = cmplx(wstar, wstar**2*k2/(par%alpha*(1 + k2)), dp)
        call check('adiabatic frequency', abs(real(omega(1) - asymptote))/wstar, 1e-8_dp)
        call check('adiabatic growth rate', &
                   abs(aimag(omega(1) - asymptote))/aimag(asymptote), 1e-6_dp)

        ! Numata Sec. II: for alpha = 1, D = 1e-4 the most unstable modes of
        ! the box 2 pi/0.15 lie on kx = 0.
        par%alpha = 1; par%kappa = 1; par%d = 1e-4_dp
        gmax = -huge(1.0_dp); gaxis = -huge(1.0_dp)
        do ix = 0, 20
            do iy = 1, 20
                call hw_linear_roots(par, 0.15_dp*ix, 0.15_dp*iy, omega, response, status)
                gmax = max(gmax, aimag(omega(1)))
                if (ix == 0) gaxis = max(gaxis, aimag(omega(1)))
            end do
        end do
        call check('most unstable mode on kx = 0', gmax - gaxis, 0.0_dp)
        if (gaxis <= 0) failures = failures + 1
        write (*, '(a,es12.5)') 'max growth rate alpha=1 kappa=1 D=1e-4: ', gaxis
    end subroutine

    subroutine hasegawa_wakatani_special_modes()
        type(hw_parameters_t) :: par
        complex(dp) :: omega(2), response(2)
        real(dp) :: nu, b
        integer :: status

        par%alpha = 0.5_dp; par%kappa = 1; par%d = 1e-3_dp; par%order = 3
        nu = par%d*0.4_dp**6
        b = par%alpha*(1 + 0.16_dp)/0.16_dp
        call hw_linear_roots(par, 0.4_dp, 0.0_dp, omega, response, status)
        if (status /= drift_root_ok) failures = failures + 1
        call check('zonal mode, original model', &
                   abs(omega(1) + (0, 1)*nu) + abs(omega(2) + (0, 1)*(b + nu)), 1e-14_dp)
        call check('zonal responses n = phi and n = -k^2 phi', &
                   abs(response(1) - 1) + abs(response(2) + 0.16_dp), 1e-14_dp)
        par%modified = .true.
        call hw_linear_roots(par, 0.4_dp, 0.0_dp, omega, response, status)
        if (status /= drift_root_decoupled) failures = failures + 1
        call check('zonal mode, modified model', maxval(abs(omega + (0, 1)*nu)), 1e-16_dp)
        call hw_linear_roots(par, 0.0_dp, 0.0_dp, omega, response, status)
        if (status /= drift_root_invalid) failures = failures + 1
    end subroutine

    subroutine hasegawa_wakatani_tangent()
        type(hw_parameters_t) :: par, plus, minus
        complex(dp) :: domega(2), wp(2), wm(2), response(2)
        real(dp), parameter :: h = 1e-6_dp, dir(5) = [0.3_dp, -0.7_dp, 2e-4_dp, 0.5_dp, 0.9_dp]
        integer :: status

        par%alpha = 0.8_dp; par%kappa = 1.3_dp; par%d = 1e-3_dp
        plus = par; minus = par
        plus%alpha = par%alpha + h*dir(1); minus%alpha = par%alpha - h*dir(1)
        plus%kappa = par%kappa + h*dir(2); minus%kappa = par%kappa - h*dir(2)
        plus%d = par%d + h*dir(3); minus%d = par%d - h*dir(3)
        call hw_linear_roots_jvp(par, 0.3_dp, 0.7_dp, dir(1), dir(2), dir(3), dir(4), &
                                 dir(5), domega, status)
        call hw_linear_roots(plus, 0.3_dp + h*dir(4), 0.7_dp + h*dir(5), wp, response, status)
        call hw_linear_roots(minus, 0.3_dp - h*dir(4), 0.7_dp - h*dir(5), wm, response, status)
        call check('HW root tangent vs central difference', &
                   maxval(abs(domega - (wp - wm)/(2*h))), 1e-8_dp)
    end subroutine

    subroutine drift_alfven()
        type(drift_alfven_parameters_t) :: par, dpar, plus, minus
        type(hw_parameters_t) :: hw
        complex(dp) :: omega(3), wp(3), wm(3), domega, hw_omega(2), response(2)
        complex(qp) :: z, defect
        real(dp), parameter :: h = 1e-6_dp
        real(dp) :: residual, kx, ky, kpar, k2, om_drift, om_alfven2, worst, scale
        integer :: status, j

        par%beta = 0.065_dp; par%mu = 0.015_dp; par%c = 5.4e-4_dp
        kx = 0.2_dp; ky = 0.4_dp; kpar = 0.3_dp
        k2 = kx**2 + ky**2
        call drift_alfven_roots(par, kpar, kx, ky, omega, residual, status)
        if (status /= drift_root_ok) failures = failures + 1
        call check('cubic residual reported', residual, 1e-14_dp)
        worst = 0
        do j = 1, 3
            z = omega(j)
            defect = z**2*(par%beta*(z - ky) + par%mu*k2*z + (0, 1)*par%c*k2) &
                     - real(kpar, qp)**2*((1 + k2)*z - ky)
            scale = abs(omega(j))**3*(par%beta + par%mu*k2) + kpar**2*ky
            worst = max(worst, real(abs(defect), dp)/scale)
        end do
        call check('roots in the dispersion function (binary128)', worst, 1e-13_dp)
        if (aimag(omega(1)) < aimag(omega(2)) .or. aimag(omega(2)) < aimag(omega(3))) &
            failures = failures + 1

        ! Alfven limit: without drive (ky = 0) and collisions the roots are 0
        ! and +-kpar sqrt((1 + k^2)/(beta + mu k^2)).
        par%c = 0
        call drift_alfven_roots(par, kpar, kx, 0.0_dp, omega, residual, status)
        call drift_alfven_limits(par%beta, par%mu, kpar, kx, 0.0_dp, om_drift, om_alfven2)
        call check('Alfven limit', abs(maxval(real(omega)) - sqrt(om_alfven2)) &
                   + abs(minval(real(omega)) + sqrt(om_alfven2)) &
                   + minval(abs(omega)) + maxval(abs(aimag(omega))), 1e-13_dp)

        ! Drift limit: kpar -> infinity leaves omega = ky/(1 + k^2).
        call drift_alfven_roots(par, 300.0_dp, kx, ky, omega, residual, status)
        call drift_alfven_limits(par%beta, par%mu, 300.0_dp, kx, ky, om_drift, om_alfven2)
        call check('drift-wave limit', minval(abs(omega - om_drift))/om_drift, 1e-5_dp)

        ! Hasegawa-Wakatani limit: beta, mu -> 0 with alpha = kpar^2/c, kappa = 1.
        par%beta = 1e-10_dp; par%mu = 1e-10_dp; par%c = 0.2_dp
        hw%alpha = kpar**2/par%c; hw%kappa = 1; hw%d = 0
        call drift_alfven_roots(par, kpar, kx, ky, omega, residual, status)
        call hw_linear_roots(hw, kx, ky, hw_omega, response, status)
        call check('Hasegawa-Wakatani limit, growing root', &
                   abs(omega(1) - hw_omega(1))/abs(hw_omega(1)), 1e-8_dp)
        call check('Hasegawa-Wakatani limit, damped root', &
                   minval(abs(omega(2:3) - hw_omega(2)))/abs(hw_omega(2)), 1e-7_dp)
        if (aimag(omega(1)) <= 0) failures = failures + 1

        ! Tangent of the most unstable root by the implicit rule.
        par%beta = 0.065_dp; par%mu = 0.015_dp; par%c = 0.05_dp
        dpar%beta = 0.4_dp; dpar%mu = -0.2_dp; dpar%c = 0.3_dp
        plus = par; minus = par
        plus%beta = par%beta + h*dpar%beta; minus%beta = par%beta - h*dpar%beta
        plus%mu = par%mu + h*dpar%mu; minus%mu = par%mu - h*dpar%mu
        plus%c = par%c + h*dpar%c; minus%c = par%c - h*dpar%c
        call drift_alfven_roots(par, kpar, kx, ky, omega, residual, status)
        call drift_alfven_root_tangent(par, kpar, kx, ky, omega(1), dpar, 0.7_dp, &
                                       -0.3_dp, 0.5_dp, domega, status)
        call drift_alfven_roots(plus, kpar + 0.7_dp*h, kx - 0.3_dp*h, ky + 0.5_dp*h, &
                                wp, residual, status)
        call drift_alfven_roots(minus, kpar - 0.7_dp*h, kx + 0.3_dp*h, ky - 0.5_dp*h, &
                                wm, residual, status)
        call check('drift-Alfven root tangent vs central difference', &
                   abs(domega - (wp(1) - wm(1))/(2*h)), 1e-8_dp)
        if (aimag(omega(1)) <= 0) failures = failures + 1
    end subroutine

    !> omega = i lambda for the eigenvalues lambda of d/dt (phi, n) = A (phi, n),
    !> ordered by decreasing growth rate. The mutant moves the drive
    !> -i kappa ky phi from the density to the vorticity equation.
    subroutine matrix_roots(par, kx, ky, mutant, omega)
        type(hw_parameters_t), intent(in) :: par
        real(dp), intent(in) :: kx, ky
        logical, intent(in) :: mutant
        complex(qp), intent(out) :: omega(2)
        complex(qp) :: a(2, 2), trace, det, root, lambda(2)
        real(qp) :: k2, nu

        k2 = real(kx, qp)**2 + real(ky, qp)**2
        nu = par%d*k2**par%order
        a(1, 1) = -par%alpha/k2 - nu
        a(1, 2) = par%alpha/k2
        a(2, 1) = par%alpha
        a(2, 2) = -par%alpha - nu
        if (mutant) then
            a(1, 1) = a(1, 1) + (0, 1)*par%kappa*ky/k2
        else
            a(2, 1) = a(2, 1) - (0, 1)*par%kappa*ky
        end if
        trace = a(1, 1) + a(2, 2)
        det = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)
        root = sqrt(trace**2 - 4*det)
        lambda = [(trace + root)/2, (trace - root)/2]
        ! The root of smaller magnitude from the product avoids cancellation.
        if (abs(lambda(1)) < abs(lambda(2))) lambda = lambda([2, 1])
        lambda(2) = det/lambda(1)
        if (real(lambda(1)) < real(lambda(2))) lambda = lambda([2, 1])
        omega = (0, 1)*lambda
    end subroutine

    real(dp) function eigen_defect(par, kx, ky, omega, response)
        type(hw_parameters_t), intent(in) :: par
        real(dp), intent(in) :: kx, ky
        complex(dp), intent(in) :: omega, response
        complex(dp) :: lambda, row(2)
        real(dp) :: k2, nu

        k2 = kx**2 + ky**2
        nu = par%d*k2**par%order
        lambda = -(0, 1)*omega
        row(1) = (-par%alpha/k2 - nu) + par%alpha/k2*response - lambda
        row(2) = par%alpha - (0, 1)*par%kappa*ky - (par%alpha + nu)*response &
                 - lambda*response
        eigen_defect = maxval(abs(row))/(abs(lambda) + par%alpha*(1 + 1/k2))
    end function

    subroutine check(name, value, tolerance)
        character(*), intent(in) :: name
        real(dp), intent(in) :: value, tolerance
        if (value <= tolerance) then
            write (*, '(a,a,es10.3)') 'ok   ', name//': ', value
        else
            write (*, '(a,a,es10.3,a,es10.3)') 'FAIL ', name//': ', value, ' > ', tolerance
            failures = failures + 1
        end if
    end subroutine

    subroutine check_fails(name, value, tolerance)
        character(*), intent(in) :: name
        real(dp), intent(in) :: value, tolerance
        if (value > tolerance) then
            write (*, '(a,a,es10.3)') 'ok   ', name//': ', value
        else
            write (*, '(a,a,es10.3)') 'FAIL ', name//': ', value
            failures = failures + 1
        end if
    end subroutine
end program test_drift_dispersion
