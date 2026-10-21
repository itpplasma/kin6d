program test_response_interfaces
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_mhd_ampere_n, only: mu0
    use kin6d_response_fluid, only: response_fluid_equilibrium_t, &
        response_fluid_stress, response_fluid_stress_adjoint, &
        response_cut_balance, response_cut_balance_adjoint
    use kin6d_response_toroidal_strip, only: toroidal_strip_t, toroidal_strip_work_t, &
        toroidal_strip_init, toroidal_strip_traces, toroidal_strip_traces_vjp, &
        toroidal_strip_domains
    implicit none
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)
    character(20) :: command
    call get_command_argument(1, command)
    select case (trim(command))
    case ('traces')
        call trace_oracle()
    case ('cut')
        call cut_oracle()
    case ('domains')
        call domain_oracle()
    case default
        error stop 'expected traces, cut or domains'
    end select
    print '(a)', 'response interfaces: all checks passed'
contains
    subroutine check(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (ok) return
        print '(a)', 'FAIL: '//message
        error stop 1
    end subroutine check

    subroutine trace_oracle()
        integer, parameter :: nr = 9, nt = 17
        complex(dp), parameter :: aa = (0.02_dp, 0.03_dp), az = (-0.03_dp, 0.01_dp)
        complex(dp), parameter :: ac = (0.04_dp, -0.02_dp), xr = (0.1_dp, 0.05_dp)
        complex(dp), parameter :: xp = (0.03_dp, -0.01_dp), xz = (-0.07_dp, 0.04_dp)
        complex(dp), parameter :: ph = (0.02_dp, -0.06_dp), gamma = (0.3_dp, 0.8_dp)
        real(dp), parameter :: bt = 0.4_dp, adiabatic = 5.0_dp/3.0_dp
        type(toroidal_strip_t) :: strip
        type(toroidal_strip_work_t) :: work
        type(response_fluid_equilibrium_t) :: eq(nr, nt)
        real(dp) :: u(nr), radius(nr, nt), inverse(2, 2, nr, nt)
        complex(dp) :: state(6, nr, nt), direction(6, nr, nt), adj(6, nr, nt)
        complex(dp) :: actual(10, nt, 2), exact(10, nt, 2), seed(10, nt, 2)
        complex(dp) :: tangent_trace(10, nt, 2), shifted(10, nt, 2)
        complex(dp) :: avec(3), magnetic(3), electric(3), phi_grad(3), xi(3), p1, pt
        complex(dp) :: lhs, rhs
        real(dp) :: angle, cosine, sine, rr, zz, normal(3), tangent(3), scale, error
        integer :: n, i, j, a, side, info
        u = [(0.3_dp + (i - 1)*0.5_dp/(nr - 1), i=1,nr)]
        do n = -3, 3, 6
            do j = 1, nt
                angle = 2*acos(-1.0_dp)*(j - 1)/nt
                cosine = cos(angle); sine = sin(angle)
                do i = 1, nr
                    rr = 2 + u(i)*cosine; zz = u(i)*sine
                    radius(i, j) = rr
                    inverse(:, 1, i, j) = [cosine, -sine/u(i)]
                    inverse(:, 2, i, j) = [sine, cosine/u(i)]
                    eq(i, j)%b0 = [0.0_dp, bt*rr, 0.0_dp]
                    eq(i, j)%j0 = [0.0_dp, 0.0_dp, 2*bt/mu0]
                    eq(i, j)%p = (4 - bt**2*rr**2)/mu0
                    eq(i, j)%gp = [-2*bt**2*rr/mu0, 0.0_dp, 0.0_dp]
                    eq(i, j)%hp = 0.0_dp
                    eq(i, j)%hp(1, 1) = -2*bt**2/mu0
                    eq(i, j)%hp(2, 2) = -2*bt**2/mu0
                    state(:, i, j) = [aa*rr**2*zz, az*rr*zz**2 + ac*rr**2, &
                        xr*rr, xp*rr**2, xz*zz, ph*rr**2*zz]
                end do
                do side = 1, 2
                    i = 1 + (side - 1)*(nr - 1)
                    rr = radius(i, j); zz = u(i)*sine
                    normal = [cosine, 0.0_dp, sine]
                    tangent = [-sine, 0.0_dp, cosine]
                    avec = [state(1, i, j), (0.0_dp, 0.0_dp), state(2, i, j)]
                    magnetic = [iu*n*(az*zz**2 + ac*rr), &
                        aa*rr**2 - az*zz**2 - 2*ac*rr, -iu*n*aa*rr*zz]
                    phi_grad = [2*ph*rr*zz, iu*n*ph*rr*zz, ph*rr**2]
                    electric = -gamma*avec - phi_grad
                    xi = [xr*rr, xp*rr**2, xz*zz]
                    p1 = 2*bt**2*xr*rr**2/mu0 &
                        - adiabatic*eq(i, j)%p*(2*xr + xz + iu*n*xp*rr)
                    pt = p1 + bt*rr*magnetic(2)/mu0
                    exact(:, j, side) = [sum(tangent*avec), (0.0_dp, 0.0_dp), &
                        sum(normal*magnetic), sum(tangent*electric), electric(2), &
                        sum(normal*xi), -pt*normal, pt]
                    exact(8, j, side) = bt*rr*sum(normal*magnetic)/mu0
                end do
            end do
            call toroidal_strip_init(strip, work, u, radius, inverse, eq, n, info)
            call check(info == 0, 'valid pressure-balanced chart')
            call toroidal_strip_traces(strip, work, gamma, state, actual)
            error = 0.0_dp
            do a = 1, 10
                scale = max(1.0_dp, maxval(abs(exact(a, :, :))))
                error = max(error, maxval(abs(actual(a, :, :) - exact(a, :, :)))/scale)
            end do
            print '(a,i3,a,es12.4)', 'trace n=', n, ' analytic error=', error
            call check(error < 2.0e-10_dp, &
                'both-ring magnetic/electric/mechanical traces')
            do j = 1, nt
                do i = 1, nr
                    do a = 1, 6
                        direction(a, i, j) = cmplx(sin(0.7_dp*a*i + j), &
                            cos(0.6_dp*a*j + i), dp)
                    end do
                end do
                do side = 1, 2
                    do a = 1, 10
                        scale = 1.0_dp
                        if (a >= 7) scale = mu0
                        seed(a, j, side) = scale*cmplx(cos(0.3_dp*a*j - side), &
                            sin(0.2_dp*a*side - j), dp)
                    end do
                end do
            end do
            call toroidal_strip_traces(strip, work, gamma, direction, tangent_trace)
            call toroidal_strip_traces(strip, work, gamma, &
                state + 0.001_dp*direction, shifted)
            call check(maxval(abs((shifted - actual)/0.001_dp - tangent_trace)) &
                /maxval(abs(tangent_trace)) < 2.0e-10_dp, 'trace JVP finite difference')
            call toroidal_strip_traces_vjp(strip, work, gamma, seed, adj)
            lhs = sum(conjg(tangent_trace)*seed); rhs = sum(conjg(direction)*adj)
            error = abs(lhs - rhs)/max(1.0_dp, abs(lhs), abs(rhs))
            print '(a,es12.4)', 'trace storage VJP error=', error
            call check(error < 2.0e-11_dp, 'both-ring trace Hermitian transpose')
        end do
    end subroutine trace_oracle

    subroutine cut_oracle()
        real(dp), parameter :: h = 1.0e-4_dp, edge_p = 10500.0_dp
        real(dp) :: b0(3), normal(3), p0, b0out(3), grad_jump, error
        real(dp) :: plus(3), minus(3), part_p, part_b(3), oracle(3)
        real(dp) :: dpin, dpout, dbin, dbout, distance, xin_part, pin, pout, bin, bout
        complex(dp) :: p1, magnetic(3), stress(4), seed(4), ap, ab(3), lhs, rhs
        complex(dp) :: inner(4), outer(4), xin, balance, pressure_in, pressure_out
        complex(dp) :: ai, ao, ax
        integer :: part
        b0 = [1.1_dp, -0.7_dp, 0.5_dp]; normal = [0.6_dp, 0.0_dp, 0.8_dp]
        p0 = edge_p
        p1 = (1700.0_dp, -800.0_dp)
        magnetic = [(0.13_dp, 0.07_dp), (-0.08_dp, 0.11_dp), (0.04_dp, -0.06_dp)]
        call response_fluid_stress(b0, normal, p1, magnetic, stress)
        do part = 1, 2
            if (part == 1) then
                part_p = real(p1, dp); part_b = real(magnetic, dp)
            else
                part_p = aimag(p1); part_b = aimag(magnetic)
            end if
            plus = nonlinear_stress(p0 + h*part_p, b0 + h*part_b, normal)
            minus = nonlinear_stress(p0 - h*part_p, b0 - h*part_b, normal)
            oracle = (plus - minus)/(2*h)
            if (part == 1) then
                error = maxval(abs(oracle - real(stress(1:3), dp)))
            else
                error = maxval(abs(oracle - aimag(stress(1:3))))
            end if
            call check(error/maxval(abs(oracle)) < 2.0e-8_dp, &
                'full Maxwell tensor derivative, including B0.n tension')
        end do
        seed = [(0.2_dp, 0.3_dp), (-0.7_dp, 0.1_dp), (0.8_dp, -0.5_dp), &
            (0.1_dp, -0.6_dp)]
        call response_fluid_stress_adjoint(b0, normal, seed, ap, ab)
        lhs = sum(conjg(stress)*seed); rhs = conjg(p1)*ap + sum(conjg(magnetic)*ab)
        call check(abs(lhs - rhs)/max(abs(lhs), abs(rhs)) < 2.0e-13_dp, &
            'stress Hermitian transpose')
        normal = [1.0_dp, 0.0_dp, 0.0_dp]
        b0 = [0.0_dp, 2.0_dp, 0.0_dp]; b0out = b0
        ! The retained cut pressure has an equilibrium defect for continuous
        ! B0 with vacuum outside. The nonlinear oracle exposes its exact size.
        call response_fluid_stress(b0, normal, cmplx(edge_p, 0.0_dp, dp), &
            cmplx(b0/2, 0.0_dp, dp), inner)
        call response_fluid_stress(b0out, normal, (0.0_dp, 0.0_dp), &
            cmplx(b0out/2, 0.0_dp, dp), outer)
        call check(abs(inner(4) - outer(4) - edge_p) < 1.0e-8_dp, &
            'nonzero cut pressure is retained as equilibrium traction defect')
        ! A balanced finite-pressure plasma/vacuum interface needs a background
        ! magnetic pressure jump (or specified exterior pressure/traction).
        b0out(2) = sqrt(b0(2)**2 + 2*mu0*edge_p)
        dpin = -3400.0_dp; dpout = 0.0_dp; dbin = 0.2_dp; dbout = -0.1_dp
        grad_jump = dpin - dpout + (b0(2)*dbin - b0out(2)*dbout)/mu0
        xin = (0.08_dp, -0.05_dp)
        pressure_in = -5.0_dp/3.0_dp*edge_p*(0.12_dp, 0.03_dp) - xin*dpin
        pressure_out = (0.0_dp, 0.0_dp)
        call response_fluid_stress(b0, normal, pressure_in, magnetic, inner)
        call response_fluid_stress(b0out, normal, pressure_out, magnetic, outer)
        balance = response_cut_balance(inner(4), outer(4), xin, grad_jump)
        do part = 1, 2
            if (part == 1) then
                xin_part = real(xin, dp); part_b = real(magnetic, dp)
                pin = real(pressure_in, dp); pout = real(pressure_out, dp)
            else
                xin_part = aimag(xin); part_b = aimag(magnetic)
                pin = aimag(pressure_in); pout = aimag(pressure_out)
            end if
            distance = h*xin_part
            bin = b0(2) + dbin*distance + h*part_b(2)
            bout = b0out(2) + dbout*distance + h*part_b(2)
            plus(1) = edge_p + dpin*distance + h*pin - dpout*distance - h*pout &
                + (bin**2 - bout**2)/(2*mu0)
            bin = b0(2) - dbin*distance - h*part_b(2)
            bout = b0out(2) - dbout*distance - h*part_b(2)
            minus(1) = edge_p - dpin*distance - h*pin + dpout*distance + h*pout &
                + (bin**2 - bout**2)/(2*mu0)
            oracle(1) = (plus(1) - minus(1))/(2*h)
            if (part == 1) then
                error = abs(oracle(1) - real(balance, dp))
            else
                error = abs(oracle(1) - aimag(balance))
            end if
            call check(error/max(1.0_dp, abs(oracle(1))) < 2.0e-8_dp, &
                'cut Lagrangian total-pressure balance with actual edge p')
        end do
        call check(abs(balance - response_cut_balance(inner(4), outer(4), &
            xin, 0.0_dp)) > 1000.0_dp, 'gradient-jump traction is consequential')
        call response_cut_balance_adjoint(grad_jump, seed(1), ai, ao, ax)
        lhs = conjg(balance)*seed(1)
        rhs = conjg(inner(4))*ai + conjg(outer(4))*ao + conjg(xin)*ax
        call check(abs(lhs - rhs)/max(1.0_dp, abs(lhs), abs(rhs)) < 2.0e-12_dp, &
            'cut balance Hermitian transpose')
        outer(4) = inner(4) + xin*grad_jump
        call check(abs(response_cut_balance(inner(4), outer(4), xin, grad_jump)) &
            < 1.0e-9_dp, 'balanced boundary with displaced normal pressure')
    end subroutine cut_oracle

    pure function nonlinear_stress(p, b, normal) result(traction)
        real(dp), intent(in) :: p, b(3), normal(3)
        real(dp) :: traction(3), tensor(3, 3)
        integer :: i, j
        do j = 1, 3
            do i = 1, 3
                tensor(i, j) = b(i)*b(j)/mu0
            end do
            tensor(j, j) = tensor(j, j) - p - sum(b**2)/(2*mu0)
        end do
        traction = matmul(tensor, normal)
    end function nonlinear_stress

    subroutine domain_oracle()
        real(dp), allocatable :: bounds(:, :)
        integer, allocatable :: owner(:)
        real(dp) :: centers(3), widths(3)
        integer :: info
        ! Synthetic unit-interval centers: two overlapping edge strips and
        ! one separated interior strip; no equilibrium/profile data.
        centers = [0.998_dp, 0.8_dp, 0.985_dp]
        widths = [0.007_dp, 0.005_dp, 0.01_dp]
        call toroidal_strip_domains(centers, widths, 0.0_dp, 1.0_dp, &
                                    bounds, owner, info)
        call check(info == 0, 'full-cut domains accepted')
        call check(owner(1) == owner(3), 'q=11/3 near-edge strip merges with q=10/3')
        call check(owner(1) /= owner(2), 'q=3 strip remains distinct')
        call check(bounds(2, owner(1)) == 1.0_dp, 'edge strip reaches physical cut')
        call coverage(centers, widths, bounds)
        call toroidal_strip_domains(centers, widths, 0.0_dp, 0.9995_dp, &
            bounds, owner, info)
        call check(info == 0, 'synthetic cut retains q=11/3')
        call check(centers(1) >= bounds(1, owner(1)) .and. &
            centers(1) <= bounds(2, owner(1)), 'q=11/3 covered on full cut')
        call check(bounds(2, owner(1)) == 0.9995_dp, 'synthetic edge extent clipped at cut')
        ! A late, wide interval bridges earlier disjoint domains. Unsorted and
        ! duplicate centers are legitimate; no rational is silently discarded.
        centers = [0.8_dp, 0.1_dp, 0.4_dp]; widths = [0.75_dp, 0.06_dp, 0.01_dp]
        call toroidal_strip_domains(centers, widths, 0.0_dp, 1.0_dp, &
                                    bounds, owner, info)
        call check(info == 0, 'unequal widths accepted')
        call check(all(owner == owner(1)), 'wide interval bridges all neighbors')
        call coverage(centers, widths, bounds)
        centers = [1.0_dp, 0.0_dp, 1.0_dp]; widths = 0.02_dp
        call toroidal_strip_domains(centers, widths, 0.0_dp, 1.0_dp, &
                                    bounds, owner, info)
        call check(info == 0, 'rational exactly at a cut retained')
        call check(owner(1) == owner(3), 'duplicate rational shares domain')
        call coverage(centers, widths, bounds)
        call toroidal_strip_domains([0.25_dp, 0.75_dp], [0.25_dp, 0.25_dp], &
            0.0_dp, 1.0_dp, bounds, owner, info)
        call check(info == 0, 'touching domains accepted')
        call check(all(owner == 1), 'touching domains merge')
        call toroidal_strip_domains([real(dp) ::], [real(dp) ::], &
            0.0_dp, 1.0_dp, bounds, owner, info)
        call check(info == 0, 'empty rational set accepted')
        call toroidal_strip_domains([1.001_dp], [0.01_dp], &
            0.0_dp, 1.0_dp, bounds, owner, info)
        call check(info /= 0, 'outside rational rejected rather than discarded')
        call toroidal_strip_domains([0.998_dp], [-0.01_dp], &
            0.0_dp, 1.0_dp, bounds, owner, info)
        call check(info /= 0, 'negative extent rejected')
        call toroidal_strip_domains([ieee_value(1.0_dp, ieee_quiet_nan)], [0.01_dp], &
            0.0_dp, 1.0_dp, bounds, owner, info)
        call check(info /= 0, 'nonfinite center rejected')
    end subroutine domain_oracle

    subroutine coverage(centers, widths, bounds)
        real(dp), intent(in) :: centers(:), widths(:), bounds(:, :)
        real(dp) :: position
        logical :: reference, actual
        integer :: i, j
        do i = 0, 1000
            position = (i + 0.5_dp)/1001.0_dp
            reference = any(abs(position - centers) <= widths)
            actual = .false.
            do j = 1, size(bounds, 2)
                if (position < bounds(1, j) .or. position > bounds(2, j)) cycle
                actual = .true.
            end do
            call check(reference .eqv. actual, 'domain union against direct membership')
        end do
    end subroutine coverage
end program test_response_interfaces
