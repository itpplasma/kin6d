program test_mhd_ampere_error
    !! Small independent coupled-energy oracle; no large response iterations.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_stiffness, meridian_load, meridian_real_to_complex
    use fortsparse, only: csc_z_t, csc_matvec, sparse_solve_once, &
        fortsparse_status_t, status_ok
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use kin6d_mhd_ampere_n, only: ampere_n_t, ampere_n_error_t, &
        ampere_n_build, ampere_n_solve, ampere_n_error, &
        ampere_n_equilibrated_error, ampere_n_free
    implicit none
    integer :: mode_n

    call mode_case(1, 1)
    call mode_case(2, -2)
    print '(a)', 'Ampere error and exterior-status checks passed'

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not. condition) then
            print '(a)', message
            error stop 1
        end if
    end subroutine require

    subroutine circle(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [3.0_dp + u*cos(theta), u*sin(theta)]
    end subroutine circle

    subroutine zero_current(x, j)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: j(3)
        j = cmplx(0.0_dp*x(1), 0.0_dp, dp)
    end subroutine zero_current

    subroutine exact_a(x, a)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: a(3)
        a = [(0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)]
        a(1) = cmplx(x(1), 0.0_dp, dp)
    end subroutine exact_a

    subroutine exact_current(x, j)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: j(3)
        j = [(0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)]
        j(1) = cmplx(real(mode_n*mode_n, dp)/x(1), 0.0_dp, dp)
    end subroutine exact_current

    subroutine mode_case(geometry, n)
        integer, intent(in) :: geometry, n
        type(flux_mesh_t) :: mesh
        type(meridian_feec_t) :: fs, fsd
        type(ampere_n_t) :: amp
        type(ampere_n_error_t) :: est, perturbed
        type(csc_z_t) :: stiffness, mass
        type(fortsparse_status_t) :: fst
        complex(dp), allocatable :: exact(:), a(:), f(:), work(:), full(:)
        complex(dp), allocatable :: ca(:), delta(:), exterior(:)
        real(dp) :: expected
        integer :: st, i

        call flux_mesh_build([0.5_dp, 1.0_dp], 4, 8, geometry, circle, mesh, st)
        call require(st == 0, 'mesh build')
        call meridian_feec_build(mesh%mm, 2, fs, st, quad_degree=8)
        call require(st == 0, 'primal FEEC build')
        call meridian_feec_build(mesh%mm, 3, fsd, st, quad_degree=8)
        call require(st == 0, 'dual FEEC build')
        call ampere_n_build(fs, 0, amp, st)
        call require(st /= 0, 'n=0 requires its own gauge and remains rejected')
        call ampere_n_build(fs, n, amp, st)
        call require(st == 0, 'coupled Ampere build')
        allocate (exact(amp%nn), a(amp%nn), full(fs%sp(1)%ndof))
        do i = 1, amp%nn
            exact(i) = cmplx(sin(real(i, dp)), cos(real(2*i, dp)), dp)
        end do
        ! Independent composition of the interior and exterior energy blocks.
        ! The algebraic manufactured load uses this composition, not the
        ! operator retained for production residual evaluation.
        call meridian_stiffness(fs, 1, n, stiffness, st)
        call require(st == 0, 'independent stiffness')
        full = (0.0_dp, 0.0_dp)
        full(:amp%nn) = exact
        work = csc_matvec(stiffness, full)
        f = work(:amp%nn)
        ca = matmul(amp%cb, exact(amp%bdof))
        exterior = matmul(conjg(transpose(amp%cb)), matmul(amp%sinv, ca))
        f(amp%bdof) = f(amp%bdof) + exterior
        call ampere_n_solve(amp, f, a, st)
        call require(st == 0, 'production solve')
        call require(maxval(abs(a - exact)) < 1.0e-9_dp, 'manufactured solution')
        call require(amp%error%discrete_energy_error < 1.0e-9_dp, &
            'production result includes small coupled energy diagnostic')
        call require(amp%error%discrete_claim == 'indicator', &
            'ordinary arithmetic is not certified-discrete')
        call require(amp%error%continuous_claim == 'unavailable-exterior-remainder', &
            'missing exterior remainder remains explicit')
        call require(.not. amp%error%exterior_remainder_available, &
            'unavailable exterior is not zero')
        a = 0.8_dp*exact
        call ampere_n_error(amp, f, a, est, st)
        call require(st == 0, 'perturbed coupled error')
        expected = 0.2_dp*sqrt(real(dot_product(exact, f), dp))
        call require(abs(est%discrete_energy_error - expected) < 1.0e-10_dp*expected, &
            'residual dual norm equals independently computed coupled energy')
        call require(est%discrete_energy_error > amp%error%discrete_energy_error, &
            'perturbed candidate increases coupled energy diagnostic')
        ! The zero-current, zero-candidate continuous problem is exactly known.
        f = (0.0_dp, 0.0_dp)
        call ampere_n_solve(amp, f, a, st, fs, fsd, zero_current)
        call require(st == 0, 'solve with interior reconstruction')
        call require(amp%error%interior%bound < 1.0e-12_dp, &
            'zero field has zero equilibrated diagnostic')
        call require(amp%error%interior_claim == &
            'indicator-fixed-tangential-trace', 'interior claim scope')
        call require(amp%error%derivative_claim == 'unavailable', &
            'primal diagnostic does not certify a derivative')
        a = 0.1_dp*exact
        perturbed = amp%error
        call ampere_n_equilibrated_error(fs, fsd, amp, a, zero_current, perturbed, st)
        call require(st == 0, 'perturbed interior reconstruction')
        call require(perturbed%interior%bound > 1.0e-4_dp, &
            'perturbed candidate increases interior diagnostic')
        call require(perturbed%continuous_claim == 'unavailable-exterior-remainder', &
            'interior reconstruction never upgrades exterior claim')
        if (geometry == 1) then
            ! Independent physical manufactured case: A = R e_R,
            ! B = -i n e_Z, mu0 J = (n^2/R) e_R. The poloidal field
            ! is represented exactly at k=2 on affine cells. A dual scalar
            ! s=0 exactly equilibrates it; its tangential trace is fixed.
            mode_n = n
            call meridian_load(fs, 1, exact_a, work, st)
            call require(st == 0, 'manufactured continuous load')
            call meridian_real_to_complex(fs%mass(1), mass)
            call sparse_solve_once(mass, work, full, fst)
            call require(status_ok(fst), 'independent L2 projection')
            a = full(:amp%nn)
            call ampere_n_equilibrated_error(fs, fsd, amp, a, exact_current, &
                perturbed, st)
            call require(st == 0, 'nonzero manufactured equilibration')
            call require(perturbed%interior%bound < 1.0e-9_dp, &
                'exact physical field recovers zero fixed-trace bound')
        end if
        delta = a(:size(a) - 1)
        call ampere_n_error(amp, f, delta, est, st)
        call require(st /= 0, 'incorrect candidate dimension rejected')
        call require(est%discrete_claim == 'unavailable', 'failed diagnostic status')
        a(1) = cmplx(ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
        call ampere_n_error(amp, f, a, est, st)
        call require(st /= 0, 'nonfinite candidate rejected')
        call require(est%discrete_claim == 'unavailable', 'nonfinite error status')
        call ampere_n_free(amp)
        call require(amp%error%discrete_claim == 'unavailable', 'free clears status')
    end subroutine mode_case
end program test_mhd_ampere_error
