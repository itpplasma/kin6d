module milca_symbolic_ideal
    ! Frozen-in (Kruskal-Kulsrud) ideal-MHD energy density per reference volume,
    ! derived in docs/ideal-map-3d.md: B = F b / J with J = sigma det F, and
    ! w(F) = |F b|^2 / (2 mu J) - p J at frozen pressure p (gamma = 0).
    use fortsym_expr, only: expr_t, operator(+), operator(-), operator(*), &
        operator(/), operator(**)
    implicit none
    private
    public :: frozen_energy
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
end module milca_symbolic_ideal
