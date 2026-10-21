module kin6d_symbolic_ideal
    ! One retained gamma=0 frozen-in energy for map and xi consumers.
    use fortsym_expr, only: expr_t, operator(+), operator(-), operator(*), &
        operator(/), operator(**)
    implicit none
    private
    public :: frozen_energy
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
end module kin6d_symbolic_ideal
