module kin6d_symbolic_ideal
    ! Authoritative ideal-MHD energy definitions for KIN6D and GLISS.
    use fortsym_expr, only: expr_t, operator(+), operator(-), operator(*), &
        operator(/), operator(**), sqrt, abs
    implicit none
    private
    public :: frozen_energy, coordinate_energy, perpendicular_energy
contains
    function frozen_energy(d, b, p, mu0, sigma, jac) result(density)
        type(expr_t), intent(in) :: d(3, 3), b(3), p, mu0, sigma
        type(expr_t), intent(out) :: jac
        type(expr_t) :: density, bx(3)
        integer :: a
        jac = sigma*(d(1, 1)*(d(2, 2)*d(3, 3) - d(2, 3)*d(3, 2)) &
            - d(1, 2)*(d(2, 1)*d(3, 3) - d(2, 3)*d(3, 1)) &
            + d(1, 3)*(d(2, 1)*d(3, 2) - d(2, 2)*d(3, 1)))
        do a = 1, 3
            bx(a) = d(a, 1)*b(1) + d(a, 2)*b(2) + d(a, 3)*b(3)
        end do
        density = (bx(1)**2 + bx(2)**2 + bx(3)**2)/(2*mu0*jac) - p*jac
    end function frozen_energy
    subroutine coordinate_energy(x, components, two_energy, divergence, three_energy)
        type(expr_t), intent(in) :: x(27)
        type(expr_t), intent(out) :: components(3), two_energy, divergence, three_energy
        type(expr_t) :: field_xi, field_eta, flux_norm, flux_wronskian
        field_xi = (x(2)*x(18) + x(1)*x(19))/x(7)
        field_eta = (x(2)*x(20) + x(1)*x(21))/x(7)
        flux_norm = x(1)**2 + x(2)**2
        flux_wronskian = x(1)*x(4) - x(3)*x(2)
        components(1) = field_xi/sqrt(x(9))
        components(2) = -sqrt(x(9))/(x(8)*x(7)) &
            * (x(7)*field_eta - flux_wronskian*x(16) &
            + x(10)*x(7)*x(16)/x(9) + x(12)*x(8)*x(7)*field_xi/x(9))
        components(3) = (x(6)*x(21) - x(5)*x(20) &
            - (x(1)*x(5) + x(2)*x(6))*x(17) &
            - (x(6)*x(4) + x(5)*x(3))*x(16) - x(11)*x(7)*x(16) &
            + x(13)*x(7)*field_xi)/(x(8)*x(7))
        two_energy = (components(1)**2 + components(2)**2 + components(3)**2 &
            - x(14)*x(16)**2)*abs(x(7))

        ! Divergence of a vector density. The three J-weighted derivatives are
        ! independent inputs so radial/angle product rules stay in the caller.
        divergence = x(22)/x(7) + (x(1)*x(23) - x(2)*x(24) &
            + x(2)*x(25) + x(1)*x(26))/(x(7)*flux_norm)
        three_energy = two_energy/x(27) + x(15)*divergence**2*abs(x(7))
    end subroutine coordinate_energy

    function perpendicular_energy(x) result(energy)
        type(expr_t), intent(in) :: x(7)
        type(expr_t) :: energy, normal, binormal
        normal = x(6)/sqrt(x(3))
        binormal = x(4)*normal + sqrt(x(3))*x(7)/x(2)
        energy = x(5)*abs(x(1))*(normal**2 + binormal**2)/2
    end function perpendicular_energy
end module kin6d_symbolic_ideal
