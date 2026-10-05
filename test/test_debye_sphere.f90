program test_debye_sphere
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortnum_interval, only: interval_t, interval
    use kin6d_debye_sphere_reference, only: debye_sphere_values, &
        debye_sphere_enclosure, debye_sphere_bands
    use kin6d_debye_sphere_leaves, only: sphere_inside, sphere_outside, &
        sphere_inside_jvp, sphere_inside_vjp, sphere_outside_jvp, &
        sphere_outside_vjp, sphere_bare_inside_jvp, sphere_bare_outside_jvp
    implicit none

    call reference_oracles()
    call interface_and_gauss()
    call derivatives()
    call continuous_bands()
    call invalid_inputs()
    write (*, '(a)') 'Matched finite-sphere equilibrium and enclosure oracles pass'
contains
    ! Independent binary128 Green-function radial integration, evaluated via
    ! antiderivatives on each side of s=r. Production uses boundary matching
    ! and an interior series, whereas this integrates the Yukawa convolution.
    function green_potential(r, k) result(phi)
        real(qp), intent(in) :: r, k
        real(qp) :: phi, a, inward, outward, image

        if (r <= 0.0_qp) then
            phi = 3*(1 - (1 + k)*exp(-k))/k**2
            return
        end if
        a = min(r, 1.0_qp)
        inward = exp(-k*r)*(exp(k*a)*(k*a - 1) + 1)/k**2
        outward = exp(k*r)*((k*a + 1)*exp(-k*a) - &
            (k + 1)*exp(-k))/k**2
        image = exp(-k*r)*(1 - (k + 1)*exp(-k))/k**2
        phi = 3*(inward + outward - image)/(2*k*r)
    end function green_potential

    function green_tangent(r, k) result(g)
        real(qp), intent(in) :: r, k
        real(qp) :: g
        real(qp), parameter :: h = 1.0e-11_qp

        g = (green_potential(r, k + h) - green_potential(r, k - h))/(2*h)
    end function green_tangent

    function gauss_potential(r) result(phi)
        real(qp), intent(in) :: r
        real(qp) :: phi

        if (r <= 1.0_qp) then
            phi = 1 + (1 - r*r)/2
        else
            phi = 1/r
        end if
    end function gauss_potential

    subroutine reference_oracles()
        real(dp), parameter :: radii(9) = &
            [0.0_dp, 0.0001_dp, 0.125_dp, 0.75_dp, 1.0_dp, &
             1.0625_dp, 2.0_dp, 3.5_dp, 6.0_dp]
        real(dp), parameter :: kappas(4) = [0.5_dp, 0.75_dp, 1.25_dp, 2.0_dp]
        real(dp) :: phi, bare, tangent
        type(interval_t) :: f, b, g, m, mb, mt
        real(qp) :: exact
        logical :: ok
        integer :: i, j

        do j = 1, size(kappas)
            do i = 1, size(radii)
                call debye_sphere_values(radii(i), kappas(j), phi, bare, tangent, ok)
                call require(ok, 'valid reference evaluation')
                exact = green_potential(real(radii(i), qp), real(kappas(j), qp))
                call close_to(phi, exact, 'independent Yukawa convolution')
                call close_to(bare, gauss_potential(real(radii(i), qp)), &
                    'Gauss potential')
                call debye_sphere_enclosure(interval(radii(i)), interval(kappas(j)), &
                    f, b, g, m, mb, mt, ok)
                call require(ok, 'strict point enclosure')
                call contains(f, exact, 'reference potential interval')
                call contains(b, gauss_potential(real(radii(i), qp)), 'bare interval')
                call contains(g, &
                    green_tangent(real(radii(i), qp), real(kappas(j), qp)), &
                    'kappa sensitivity interval')
                call require(f%lo > 0.0_dp, 'positive matched potential')
                call require(phi < bare, 'screening lowers bare potential')
                call require(g%hi < 0.0_dp, 'increased screening lowers potential')
            end do
        end do
        ! A box crossing the source boundary and a nontrivial parameter range.
        call debye_sphere_enclosure(interval(0.625_dp, 1.5_dp), &
            interval(0.625_dp, 1.75_dp), f, b, g, m, mb, mt, ok)
        call require(ok, 'interface-crossing physical box')
        do j = 0, 12
            do i = 0, 14
                exact = green_potential(0.625_qp + 0.875_qp*i/14, &
                    0.625_qp + 1.125_qp*j/12)
                call contains(f, exact, 'whole-box interior oracle')
            end do
        end do
    end subroutine reference_oracles

    subroutine interface_and_gauss()
        real(dp) :: inside, outside, dr_inside, dr_outside, unused
        real(dp), parameter :: kappas(3) = [0.5_dp, 1.0_dp, 2.0_dp]
        integer :: j

        do j = 1, size(kappas)
            call sphere_inside(1.0_dp, kappas(j), inside)
            call sphere_outside(1.0_dp, kappas(j), outside)
            call close_to(inside, real(outside, qp), 'matched potential at R')
            call sphere_inside_jvp(1.0_dp, kappas(j), 1.0_dp, 0.0_dp, dr_inside)
            call sphere_outside_jvp(1.0_dp, kappas(j), 1.0_dp, 0.0_dp, dr_outside)
            call close_to(dr_inside, real(dr_outside, qp), &
                'no artificial surface charge')
            call sphere_inside_jvp(0.0_dp, kappas(j), 1.0_dp, 0.0_dp, unused)
            call require(abs(unused) < 1.0e-15_dp, 'regular central field')
        end do
        call sphere_bare_inside_jvp(0.5_dp, 1.0_dp, dr_inside)
        call close_to(dr_inside, -0.5_qp, 'Gauss enclosed charge inside')
        call sphere_bare_outside_jvp(2.0_dp, 1.0_dp, dr_outside)
        call close_to(dr_outside, -0.25_qp, 'Gauss total charge outside')
    end subroutine interface_and_gauss

    subroutine derivatives()
        real(dp) :: dr, dk, tangent, gr, gk, r, k
        real(qp) :: oracle, h
        integer :: i

        dr = 0.375_dp; dk = -0.625_dp; k = 1.25_dp
        h = 1.0e-11_qp
        do i = 1, 2
            if (i == 1) then
                r = 0.375_dp
                call sphere_inside_jvp(r, k, dr, dk, tangent)
                call sphere_inside_vjp(r, k, 0.75_dp, gr, gk)
            else
                r = 2.125_dp
                call sphere_outside_jvp(r, k, dr, dk, tangent)
                call sphere_outside_vjp(r, k, 0.75_dp, gr, gk)
            end if
            oracle = (green_potential(real(r, qp) + h*dr, real(k, qp) + h*dk) - &
                green_potential(real(r, qp) - h*dr, real(k, qp) - h*dk))/(2*h)
            call close_to(tangent, oracle, 'independent directional derivative')
            call close_to(0.75_dp*tangent, real(gr, qp)*dr + real(gk, qp)*dk, &
                'Euclidean transpose pairing')
        end do
    end subroutine derivatives

    subroutine continuous_bands()
        real(dp) :: r(7), phi(7), bare(7), tangent(7)
        real(dp) :: lower(3, 7), upper(3, 7), work(3, 6)
        real(qp) :: theta, point, lo, hi, oracle
        logical :: ok
        integer :: cell, j, q

        r = [0.0_dp, 0.125_dp, 0.5_dp, 0.875_dp, 1.2_dp, 2.7_dp, 6.0_dp]
        call debye_sphere_bands(r, 1.75_dp, phi, bare, tangent, lower, upper, work, ok)
        call require(ok, 'nonuniform grid including interface-crossing chord')
        do cell = 1, size(r) - 1
            do j = 0, 64
                theta = real(j, qp)/64
                point = (1 - theta)*real(r(cell), qp) + theta*real(r(cell + 1), qp)
                do q = 1, 3
                    lo = (1 - theta)*real(lower(q, cell), qp) + &
                        theta*real(lower(q, cell + 1), qp)
                    hi = (1 - theta)*real(upper(q, cell), qp) + &
                        theta*real(upper(q, cell + 1), qp)
                    select case (q)
                    case (1)
                        oracle = green_potential(point, 1.75_qp)
                    case (2)
                        oracle = gauss_potential(point)
                    case (3)
                        oracle = green_tangent(point, 1.75_qp)
                    end select
                    call require(oracle >= lo .and. oracle <= hi, &
                        'entire connected radial band contains oracle')
                end do
            end do
        end do
        call require(gauss_potential(0.5_qp) > &
            (gauss_potential(0.0_qp) + gauss_potential(1.0_qp))/2, &
            'endpoint-only bare band would fail without curvature padding')
    end subroutine continuous_bands

    subroutine invalid_inputs()
        real(dp) :: nan, phi, bare, tangent
        real(dp) :: r(3), f(3), b(3), g(3), lo(3, 3), hi(3, 3), work(3, 2)
        type(interval_t) :: fi, bi, gi, m, mb, mt
        logical :: ok

        nan = ieee_value(1.0_dp, ieee_quiet_nan)
        call debye_sphere_values(-0.01_dp, 1.0_dp, phi, bare, tangent, ok)
        call require(.not. ok, 'negative radius rejected')
        call debye_sphere_values(0.0_dp, 0.0_dp, phi, bare, tangent, ok)
        call require(.not. ok, 'unadmitted kappa=0 rejected, bare is separate limit')
        call debye_sphere_values(1.0_dp, nan, phi, bare, tangent, ok)
        call require(.not. ok, 'NaN parameter rejected')
        call debye_sphere_enclosure(interval(nan), interval(1.0_dp), &
            fi, bi, gi, m, mb, mt, ok)
        call require(.not. ok, 'NaN interval rejected')
        call debye_sphere_enclosure(interval(0.0_dp, 6.01_dp), interval(1.0_dp), &
            fi, bi, gi, m, mb, mt, ok)
        call require(.not. ok, 'unsupported radial domain rejected')
        call debye_sphere_enclosure(interval(0.0_dp), interval(0.49_dp, 0.8_dp), &
            fi, bi, gi, m, mb, mt, ok)
        call require(.not. ok, 'unsupported parameter box rejected')
        r = [0.0_dp, 1.0_dp, 2.0_dp]
        call debye_sphere_bands(r, 1.0_dp, f(:2), b, g, lo, hi, work, ok)
        call require(.not. ok, 'shape mismatch rejected before indexing')
        r = [0.0_dp, 1.0_dp, 1.0_dp]
        call debye_sphere_bands(r, 1.0_dp, f, b, g, lo, hi, work, ok)
        call require(.not. ok, 'duplicate knots rejected')
        r = [0.0_dp, nan, 2.0_dp]
        call debye_sphere_bands(r, 1.0_dp, f, b, g, lo, hi, work, ok)
        call require(.not. ok, 'NaN grid rejected')
    end subroutine invalid_inputs

    subroutine contains(box, oracle, message)
        type(interval_t), intent(in) :: box
        real(qp), intent(in) :: oracle
        character(*), intent(in) :: message

        call require(oracle >= real(box%lo, qp) .and. &
            oracle <= real(box%hi, qp), message)
    end subroutine contains

    subroutine close_to(actual, expected, message)
        real(dp), intent(in) :: actual
        real(qp), intent(in) :: expected
        character(*), intent(in) :: message

        call require(abs(real(actual, qp) - expected) < &
            3.0e-13_qp*(1 + abs(expected)), message)
    end subroutine close_to

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (condition) return
        write (*, '(a)') 'FAILED: '//message
        error stop 1
    end subroutine require
end program test_debye_sphere
