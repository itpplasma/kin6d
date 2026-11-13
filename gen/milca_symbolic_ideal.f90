module milca_symbolic_ideal
    ! Frozen-in (Kruskal-Kulsrud) ideal-MHD energy density per reference volume,
    ! derived in docs/ideal-map-3d.md: B = F b / J with J = sigma det F, and
    ! w(F) = |F b|^2 / (2 mu J) - p J at frozen pressure p (gamma = 0).
    use fortsym_expr, only: expr_t, operator(+), operator(-), operator(*), &
        operator(/), operator(**)
    implicit none
    private
    public :: frozen_energy, adiabatic_internal_energy
contains
    function frozen_energy(d, b, p, mu, sigma, jac) result(density)
        type(expr_t), intent(in) :: d(3, 3), b(3), p, mu, sigma
        type(expr_t), intent(out) :: jac
        type(expr_t) :: density, fb(3), fb2
        integer :: i
        jac = sigma*(d(1, 1)*(d(2, 2)*d(3, 3) - d(2, 3)*d(3, 2)) &
            - d(1, 2)*(d(2, 1)*d(3, 3) - d(2, 3)*d(3, 1)) &
            + d(1, 3)*(d(2, 1)*d(3, 2) - d(2, 2)*d(3, 1)))
        do i = 1, 3
            fb(i) = d(i, 1)*b(1) + d(i, 2)*b(2) + d(i, 3)*b(3)
        end do
        fb2 = fb(1)**2 + fb(2)**2 + fb(3)**2
        density = fb2/(2*mu*jac) - p*jac
    end function frozen_energy

    function adiabatic_internal_energy(p, gam, jac0, jac) result(density)
        ! Internal energy per reference volume of an adiabatic displacement:
        ! mass conservation and p/rho^gam = const freeze p J^gam along the
        ! map, so the physical pressure is p (jac0/J)^gam when it is p at the
        ! state with Jacobian jac0, and the energy is that pressure times
        ! J/(gam - 1): p jac0 (J/jac0)^(1-gam)/(gam-1).
        type(expr_t), intent(in) :: p, gam, jac0, jac
        type(expr_t) :: density
        density = p*jac0*(jac/jac0)**(1 - gam)/(gam - 1)
    end function adiabatic_internal_energy
end module milca_symbolic_ideal
