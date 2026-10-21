    pure function BILINEAR(b, p, mu, a, c) result(value)
        ! Sesquilinear extension: conjugation on the test jet c.
        real(RK), intent(in) :: b(3), p, mu
        complex(RK), intent(in) :: a(3, 3), c(3, 3)
        complex(RK) :: value
        real(RK) :: ar(9), ai(9), cr(9), ci(9), rr, ii, ir, ri
        ar = reshape(real(a, RK), [9]); ai = reshape(aimag(a), [9])
        cr = reshape(real(c, RK), [9]); ci = reshape(aimag(c), [9])
        call FORM(b(1), b(2), b(3), p, mu, &
            ar(1), ar(2), ar(3), ar(4), ar(5), ar(6), ar(7), ar(8), ar(9), &
            cr(1), cr(2), cr(3), cr(4), cr(5), cr(6), cr(7), cr(8), cr(9), rr)
        call FORM(b(1), b(2), b(3), p, mu, &
            ai(1), ai(2), ai(3), ai(4), ai(5), ai(6), ai(7), ai(8), ai(9), &
            ci(1), ci(2), ci(3), ci(4), ci(5), ci(6), ci(7), ci(8), ci(9), ii)
        call FORM(b(1), b(2), b(3), p, mu, &
            ai(1), ai(2), ai(3), ai(4), ai(5), ai(6), ai(7), ai(8), ai(9), &
            cr(1), cr(2), cr(3), cr(4), cr(5), cr(6), cr(7), cr(8), cr(9), ir)
        call FORM(b(1), b(2), b(3), p, mu, &
            ar(1), ar(2), ar(3), ar(4), ar(5), ar(6), ar(7), ar(8), ar(9), &
            ci(1), ci(2), ci(3), ci(4), ci(5), ci(6), ci(7), ci(8), ci(9), ri)
        value = cmplx(rr + ii, ir - ri, RK)
    end function BILINEAR

    pure subroutine GRADIENT(n, r, xi, dr, dz, a, status)
        ! Components are (R, phi, Z), exp(+i n phi); r > 0.
        integer, intent(in) :: n
        real(RK), intent(in) :: r
        complex(RK), intent(in) :: xi(3), dr(3), dz(3)
        complex(RK), intent(out) :: a(3, 3)
        integer, intent(out) :: status
        real(RK) :: x(3, 2), xr(3, 2), xz(3, 2), g(9, 2)
        status = XI_BAD_INPUT
        a = cmplx(0, 0, RK)
        if (.not. ieee_is_finite(r)) return
        if (r <= 0) return
        x(:, 1) = real(xi, RK); x(:, 2) = aimag(xi)
        xr(:, 1) = real(dr, RK); xr(:, 2) = aimag(dr)
        xz(:, 1) = real(dz, RK); xz(:, 2) = aimag(dz)
        if (.not. all(ieee_is_finite(x))) return
        if (.not. all(ieee_is_finite(xr))) return
        if (.not. all(ieee_is_finite(xz))) return
        call GRAD(x(1,1), x(2,1), x(3,1), x(1,2), x(2,2), x(3,2), &
            xr(1,1), xr(2,1), xr(3,1), xr(1,2), xr(2,2), xr(3,2), &
            xz(1,1), xz(2,1), xz(3,1), xz(1,2), xz(2,2), xz(3,2), r, real(n,RK), &
            g(1,1), g(2,1), g(3,1), g(4,1), g(5,1), g(6,1), g(7,1), g(8,1), g(9,1), &
            g(1,2), g(2,2), g(3,2), g(4,2), g(5,2), g(6,2), g(7,2), g(8,2), g(9,2))
        a = reshape(cmplx(g(:, 1), g(:, 2), RK), [3, 3])
        if (.not. all(ieee_is_finite(real(a,RK)))) return
        if (.not. all(ieee_is_finite(aimag(a)))) return
        status = XI_OK
    end subroutine GRADIENT

    subroutine VOLUME(weight, b, p, mu, a, energy, report)
        ! weight includes the physical volume measure and toroidal period.
        ! Returns half the complex modal Hessian. See docs/ideal-xi.md.
        real(RK), intent(in) :: weight(:), b(3, size(weight)), p(size(weight)), mu
        complex(RK), intent(in) :: a(3, 3, size(weight))
        real(RK), intent(out) :: energy
        type(xi_energy_report_t), intent(out) :: report
        real(RK) :: total
        integer :: q
        report = xi_energy_report_t()
        energy = 0
        if (.not. ieee_is_finite(mu)) return
        if (mu <= 0 .or. size(weight) == 0) return
        if (.not. all(ieee_is_finite(weight))) return
        if (any(weight < 0)) return
        if (.not. all(ieee_is_finite(b))) return
        if (.not. all(ieee_is_finite(p))) return
        if (.not. all(ieee_is_finite(real(a,RK)))) return
        if (.not. all(ieee_is_finite(aimag(a)))) return
        total = 0
        !$omp parallel do default(none) shared(weight,b,p,mu,a) private(q) reduction(+:total)
        do q = 1, size(weight)
            total = total + weight(q)*real(BILINEAR(b(:,q), p(q), mu, &
                a(:,:,q), a(:,:,q)), RK)
        end do
        !$omp end parallel do
        energy = total/2
        if (.not. ieee_is_finite(energy)) return
        report%status = XI_OK
    end subroutine VOLUME
