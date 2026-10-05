! Shared symbolic physical definitions for KIN6D generators.
! Minimal Q1 version: the foundation (F2) owns and may extend this module.
module kin6d_symbolic_physics
    use fortsym_expr, only: expr_t, operator(+), operator(*), operator(/), operator(**)
    implicit none
    private
    public :: coulomb_potential, radial_effective_potential
contains
    ! Bare two-body Coulomb energy g/r, g = q_i q_j/(4 pi eps0) in declared units.
    function coulomb_potential(g, r) result(v)
        type(expr_t), intent(in) :: g, r
        type(expr_t) :: v
        v = g/r
    end function coulomb_potential

    ! Partial-wave radial potential g/r + hbar^2 l(l+1)/(2 mu r^2).
    function radial_effective_potential(g, r, l, mu, hbar) result(v)
        type(expr_t), intent(in) :: g, r, l, mu, hbar
        type(expr_t) :: v
        v = coulomb_potential(g, r) + hbar**2*l*(l + 1)/(2*mu*r**2)
    end function radial_effective_potential
end module kin6d_symbolic_physics
