!> Hypernetted-chain (HNC) integral equation of the classical OCP (H1a
!> reduced model). Wigner-Seitz units a = 1, n = 3/(4 pi), beta v = Gamma/r.
!>
!> Long-range Coulomb part handled by the Ng (1974) renormalization:
!> beta v = beta v_s + beta v_L with the Ewald split of parameter alpha,
!> c_s = c + beta v_L, gamma_s = h - c_s. Closure (r space):
!> g = exp(-beta v_s + gamma_s); Ornstein-Zernike (k space):
!> c = c_s - beta v_L(k), h = c/(1 - n c), gamma_s = h - c_s. Kernels are the
!> generated leaves of gen/gen_correlation_hierarchy.f90.
!>
!> Radial transform on r_i = i dr, k_j = j pi/R (i, j = 1..m-1, dr = R/m):
!> F(k_j) = (4 pi dr/k_j) sum_i r_i f_i sin(pi i j/m) and the inverse with
!> (dk/(2 pi^2 r_i)) sum_j k_j F_j sin(pi i j/m), dk = pi/R. The sine sums
!> are a DST-I evaluated by one FortNum complex FFT of length 2m on the odd
!> extension. Iteration: damped Picard on gamma_s with step halving when
!> the update norm grows; converged when max |update| <= tolerance.
module kin6d_ocp_hnc
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_fft, only: fortnum_fft_plan_t, fft_c2c_plan_init, fft_c2c
    use kin6d_correlation_hierarchy_leaves, only: ch_hnc_closure, &
        ch_hnc_long_fourier, ch_hnc_short_potential, ch_oz_indirect, &
        ch_excess_energy_density
    implicit none
    private
    public :: hnc_t, hnc_init, hnc_solve, hnc_dst

    real(dp), parameter :: pi = 3.14159265358979323846_dp, density = 3/(4*pi)

    type :: hnc_t
        integer :: m = 0, iterations = 0
        real(dp) :: radius = 0, alpha = 1, gamma = 0, update = 0, energy = 0
        real(dp), allocatable :: r(:), k(:), gs(:), g(:), s(:)
        type(fortnum_fft_plan_t) :: plan
    end type
contains
    subroutine hnc_init(h, m, radius, alpha)
        type(hnc_t), intent(out) :: h
        integer, intent(in) :: m
        real(dp), intent(in) :: radius, alpha
        integer :: i
        h%m = m; h%radius = radius; h%alpha = alpha
        allocate (h%r(m - 1), h%k(m - 1), h%gs(m - 1), h%g(m - 1), h%s(m - 1))
        do i = 1, m - 1
            h%r(i) = i*radius/m; h%k(i) = i*pi/radius
        end do
        h%gs = 0
        call fft_c2c_plan_init(h%plan, 2*m)
    end subroutine

    !> out_j = sum_{i=1}^{m-1} a_i sin(pi i j/m), j = 1..m-1.
    subroutine hnc_dst(h, a, out)
        type(hnc_t), intent(in) :: h
        real(dp), intent(in) :: a(:)
        real(dp), intent(out) :: out(:)
        complex(dp) :: z(0:2*h%m - 1)
        integer :: m
        m = h%m
        z = 0
        z(1:m - 1) = a
        z(2*m - 1:m + 1:-1) = -a
        call fft_c2c(z, -1, h%plan)
        out = -aimag(z(1:m - 1))/2
    end subroutine

    !> Solve at coupling gamma from the current gamma_s (continuation).
    subroutine hnc_solve(h, gamma, tolerance, max_iterations, ok)
        type(hnc_t), intent(inout) :: h
        real(dp), intent(in) :: gamma, tolerance
        integer, intent(in) :: max_iterations
        logical, intent(out) :: ok
        real(dp) :: cs(h%m - 1), ck(h%m - 1), vl(h%m - 1), gk(h%m - 1), new(h%m - 1)
        real(dp) :: lambda, previous, dr, dk, vs(h%m - 1), e(h%m - 1)
        integer :: it
        h%gamma = gamma; ok = .false.
        dr = h%radius/h%m; dk = pi/h%radius
        call ch_hnc_long_fourier(h%k**2, gamma, h%alpha, vl)
        lambda = 0.5_dp; previous = huge(1.0_dp)
        do it = 1, max_iterations
            call ch_hnc_closure(h%r, gamma, h%alpha, h%gs, cs)
            call hnc_dst(h, h%r*cs, ck)
            ck = 4*pi*dr*ck/h%k
            call ch_oz_indirect(ck, vl, density, gk)
            call hnc_dst(h, h%k*gk, new)
            new = dk*new/(2*pi**2*h%r)
            h%update = maxval(abs(new - h%gs))
            if (h%update > previous) lambda = max(1.0e-3_dp, lambda/2)
            previous = h%update
            h%gs = h%gs + lambda*(new - h%gs)
            h%iterations = it
            if (h%update <= tolerance) then
                ok = .true.
                exit
            end if
        end do
        call ch_hnc_short_potential(h%r, gamma, h%alpha, vs)
        h%g = exp(h%gs - vs)
        ! S(k) = 1 + n h(k) from the final closure/OZ state.
        call ch_hnc_closure(h%r, gamma, h%alpha, h%gs, cs)
        call hnc_dst(h, h%r*cs, ck)
        ck = 4*pi*dr*ck/h%k - vl
        h%s = 1 + density*ck/(1 - density*ck)
        call ch_excess_energy_density(h%r, gamma, h%g - 1, density, e)
        h%energy = sum(e)*dr
    end subroutine
end module kin6d_ocp_hnc
