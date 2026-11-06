!> Pitch-angle (Lorentz) scattering in lambda = v_perp^2/(v^2 B) at one
!! point of a field line, from its weak form per sign sigma of v_par,
!!   int g C_L f (B/|xi|) d lambda = -2 nu int lambda |xi| g' f' d lambda,
!! with |xi| = sqrt(1 - lambda B). Elements are bounded by lambda levels
!! and carry continuous Lagrange bases of degree p in t = sqrt(1/B - lambda)
!! for each sigma; both sides share the reflection node lambda = 1/B, so f
!! is continuous across the turning point. Elements are uniform in
!! t = sqrt(1/B - lambda), which also carries the quadrature: the measure is
!! smooth in t, so Gauss-Legendre with 2p+2 points per element integrates
!! both forms exactly. Returns the symmetric positive semidefinite
!! stiffness K = -(weak C_L) and the mass M; C_L = -M^{-1} K.
module kin6d_dk_pitch
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_polynomial, only: lagrange_weights, lagrange_deriv_weights
    use kin6d_dk_pitch_leaves, only: dk_pitch_weights
    implicit none
    private
    public :: pitch_operator_t, build_pitch_operator

    type :: pitch_operator_t
        integer :: n = 0 ! unknowns, both sigma
        real(dp), allocatable :: lambda(:), sigma(:)
        real(dp), allocatable :: stiffness(:, :), mass(:, :)
    end type

contains

    subroutine build_pitch_operator(b, nu, n_el, p, op)
        real(dp), intent(in) :: b, nu
        integer, intent(in) :: n_el, p
        type(pitch_operator_t), intent(out) :: op
        integer :: n_side, e, q, a, c, side, nq
        integer, allocatable :: dof(:)
        real(dp) :: t_max, t0, t1, lam, w_mass, w_stiff, s
        real(dp), allocatable :: node(:), tnode(:), tq(:), wq(:), phi(:), dphi(:)

        n_side = n_el*p + 1
        op%n = 2*n_side - 1
        allocate (tnode(p + 1))
        allocate (op%lambda(op%n), op%sigma(op%n), dof(p + 1), node(p + 1), &
            phi(p + 1), dphi(p + 1))
        allocate (op%stiffness(op%n, op%n), op%mass(op%n, op%n))
        op%stiffness = 0; op%mass = 0
        nq = 2*p + 2
        allocate (tq(nq), wq(nq))
        t_max = 1/sqrt(b)

        do side = 1, 2
            s = merge(1.0_dp, -1.0_dp, side == 1)
            do e = 1, n_el
                t0 = t_max*(1 - real(e - 1, dp)/n_el)
                t1 = t_max*(1 - real(e, dp)/n_el)
                do a = 0, p
                    tnode(a + 1) = t0 + (t1 - t0)*a/p
                    node(a + 1) = 1/b - tnode(a + 1)**2
                    dof(a + 1) = global_dof((e - 1)*p + a)
                end do
                ! Element boundaries are lambda levels; the basis is polynomial
                ! in t, d/dlambda = -(1/(2t)) d/dt, so both parities of P_l(xi)
                ! (xi = sigma sqrt(B) t) are represented.
                op%lambda(dof) = node
                op%sigma(dof) = s
                call gauss_legendre_ab(nq, t1, t0, tq, wq)
                do q = 1, nq
                    call dk_pitch_weights(tq(q), b, nu, lam, w_mass, w_stiff)
                    call lagrange_weights(p + 1, tq(q), tnode, phi)
                    call lagrange_deriv_weights(p + 1, tq(q), tnode, dphi)
                    dphi = -dphi/(2*tq(q))
                    do c = 1, p + 1
                        do a = 1, p + 1
                            op%mass(dof(a), dof(c)) = op%mass(dof(a), dof(c)) &
                                + wq(q)*w_mass*phi(a)*phi(c)
                            op%stiffness(dof(a), dof(c)) = &
                                op%stiffness(dof(a), dof(c)) &
                                + wq(q)*w_stiff*dphi(a)*dphi(c)
                        end do
                    end do
                end do
            end do
        end do
        op%sigma(n_side) = 0
    contains
        ! Local index j = 0 (lambda = 0) ... n_side - 1 (reflection node).
        pure integer function global_dof(j)
            integer, intent(in) :: j
            if (side == 1 .or. j == n_side - 1) then
                global_dof = j + 1
            else
                global_dof = n_side + j + 1
            end if
        end function global_dof
    end subroutine build_pitch_operator

end module kin6d_dk_pitch
