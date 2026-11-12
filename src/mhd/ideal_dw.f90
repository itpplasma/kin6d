module milca_ideal_dw
    ! Pointwise delta-W bilinear densities in physical Cartesian components
    ! (docs/stability.md, P1e). g(i,j) = d_j B_i, a(i,j) = d_j xi_i. Values are
    ! twice the delta-W density; complex versions conjugate the test field.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_ideal_dw_leaves, only: dw_bernstein, dw_lagrange
    implicit none
    private
    public :: dw_bernstein_form, dw_lagrange_form
    interface dw_bernstein_form
        module procedure bernstein_real, bernstein_complex
    end interface
    interface dw_lagrange_form
        module procedure lagrange_real, lagrange_complex
    end interface
contains
    pure function bernstein_real(b, g, p, gp, gam, mu, xi, a, eta, c) result(value)
        real(dp), intent(in) :: b(3), g(3, 3), p, gp(3), gam, mu
        real(dp), intent(in) :: xi(3), a(3, 3), eta(3), c(3, 3)
        real(dp) :: value
        call dw_bernstein(b(1), b(2), b(3), g(1,1), g(2,1), g(3,1), g(1,2), &
            g(2,2), g(3,2), g(1,3), g(2,3), g(3,3), p, gp(1), gp(2), gp(3), gam, &
            mu, xi(1), xi(2), xi(3), a(1,1), a(2,1), a(3,1), a(1,2), a(2,2), &
            a(3,2), a(1,3), a(2,3), a(3,3), eta(1), eta(2), eta(3), c(1,1), &
            c(2,1), c(3,1), c(1,2), c(2,2), c(3,2), c(1,3), c(2,3), c(3,3), value)
    end function bernstein_real

    pure function bernstein_complex(b, g, p, gp, gam, mu, xi, a, eta, c) &
            result(value)
        real(dp), intent(in) :: b(3), g(3, 3), p, gp(3), gam, mu
        complex(dp), intent(in) :: xi(3), a(3, 3), eta(3), c(3, 3)
        complex(dp) :: value
        real(dp) :: rr, ii, ir, ri
        rr = bernstein_real(b, g, p, gp, gam, mu, real(xi), real(a), &
            real(eta), real(c))
        ii = bernstein_real(b, g, p, gp, gam, mu, aimag(xi), aimag(a), &
            aimag(eta), aimag(c))
        ir = bernstein_real(b, g, p, gp, gam, mu, aimag(xi), aimag(a), &
            real(eta), real(c))
        ri = bernstein_real(b, g, p, gp, gam, mu, real(xi), real(a), &
            aimag(eta), aimag(c))
        value = cmplx(rr + ii, ir - ri, dp)
    end function bernstein_complex

    pure function lagrange_real(b, p, gam, mu, a, c) result(value)
        real(dp), intent(in) :: b(3), p, gam, mu, a(3, 3), c(3, 3)
        real(dp) :: value
        call dw_lagrange(b(1), b(2), b(3), p, gam, mu, a(1,1), a(2,1), a(3,1), &
            a(1,2), a(2,2), a(3,2), a(1,3), a(2,3), a(3,3), c(1,1), c(2,1), &
            c(3,1), c(1,2), c(2,2), c(3,2), c(1,3), c(2,3), c(3,3), value)
    end function lagrange_real

    pure function lagrange_complex(b, p, gam, mu, a, c) result(value)
        real(dp), intent(in) :: b(3), p, gam, mu
        complex(dp), intent(in) :: a(3, 3), c(3, 3)
        complex(dp) :: value
        real(dp) :: rr, ii, ir, ri
        rr = lagrange_real(b, p, gam, mu, real(a), real(c))
        ii = lagrange_real(b, p, gam, mu, aimag(a), aimag(c))
        ir = lagrange_real(b, p, gam, mu, aimag(a), real(c))
        ri = lagrange_real(b, p, gam, mu, real(a), aimag(c))
        value = cmplx(rr + ii, ir - ri, dp)
    end function lagrange_complex

end module milca_ideal_dw
