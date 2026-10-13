module kin6d_nonnested_residual
    ! D3 on P_k displacement components times Fourier modes. The Riesz
    ! norm is the sum of scalar cylindrical H1 seminorms, with zero trace.
    ! This is a finite-test-space dual norm, not a continuous majorant.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_meridian_feec_n, only: meridian_qpoint
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use kin6d_mhd_mesh, only: flux_mesh_geometry
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_riesz
    use kin6d_equilibrium_residual, only: equilibrium_residual_t
    use kin6d_equilibrium_core_leaves, only: maxwell_stress_delta
    use kin6d_nonnested_leaves, only: axisymmetric_pressure
    implicit none
    private
    public :: nn_residual
contains
    subroutine nn_residual(field, test_modes, report, status)
        type(free_b_field_t), intent(in) :: field
        integer, intent(in) :: test_modes
        type(equilibrium_residual_t), intent(out) :: report
        integer, intent(out) :: status
        real(dp), allocatable :: xi(:), eta(:), weight(:), means(:, :)
        complex(dp), allocatable :: cells(:, :, :, :), rhs(:, :, :), riesz(:)
        complex(dp), allocatable :: phases(:, :)
        real(dp) :: two_pi, dual2
        integer :: t, n, j, st, failure, nt, nphi, modes
        status = 1
        if (.not. allocated(field%b)) return
        if (test_modes < 0) return
        nt = field%fs%mm%topology%n_triangles
        modes = ubound(field%b, 2)
        nphi = max(1, 2*modes + test_modes + 1)
        two_pi = 2*acos(-1.0_dp)
        call triangle_duffy_quadrature(2*field%fs%k + 2*field%fs%g + 2, &
            xi, eta, weight, st)
        if (st /= 0) return
        if (size(weight) /= field%fs%nq) return
        allocate (phases(0:max(modes, test_modes), nphi))
        do j = 1, nphi
            do n = 0, ubound(phases, 1)
                phases(n, j) = exp(cmplx(0.0_dp, n*two_pi*(j - 1)/nphi, dp))
            end do
        end do
        allocate (cells(field%fs%nl0, 3, 0:test_modes, nt), means(3, nt))
        failure = 0
        !$omp parallel do default(shared) private(t, st) schedule(static) &
        !$omp& reduction(max:failure)
        do t = 1, nt
            call cell_load(t, cells(:, :, :, t), means(:, t), st)
            failure = max(failure, st)
        end do
        !$omp end parallel do
        if (failure /= 0) return
        allocate (rhs(field%fs%n0, 3, 0:test_modes), riesz(field%fs%n0))
        rhs = 0
        do t = 1, nt
            rhs(field%fs%sp(0)%dof(:, t), :, :) = &
                rhs(field%fs%sp(0)%dof(:, t), :, :) + cells(:, :, :, t)
        end do
        dual2 = 0
        do n = 0, test_modes
            do j = 1, 3
                where (field%fs%sp(0)%bnd) rhs(:, j, n) = 0
                call nn_riesz(field, n, rhs(:, j, n), riesz, st)
                if (st /= 0) return
                dual2 = dual2 + merge(1, 2, n == 0) &
                    *real(dot_product(rhs(:, j, n), riesz), dp)
            end do
        end do
        report%volume = two_pi*sum(means(1, :))
        report%b2_mean = two_pi*sum(means(2, :))/report%volume
        report%p_mean = two_pi*sum(means(3, :))/report%volume
        report%sigma = report%b2_mean/field%mu0 + report%p_mean
        report%dual = sqrt(max(0.0_dp, two_pi*dual2))
        if (report%sigma <= 0) return
        report%normalized = report%dual/(report%sigma*sqrt(report%volume))
        status = 0
    contains
        subroutine cell_load(t, load, mean, status)
            integer, intent(in) :: t
            complex(dp), intent(out) :: load(:, :, 0:)
            real(dp), intent(out) :: mean(3)
            integer, intent(out) :: status
            real(dp) :: v2(3, field%fs%sp(2)%nloc), v0(3, field%fs%nl0)
            real(dp) :: grad(2, field%fs%nl0), x(2), jac(2, 2), det, wq
            real(dp) :: b(3), b0, p, psi, tensor(3, 3), tau(9)
            complex(dp) :: bn(3, 0:modes), average(3, 3, 0:test_modes)
            integer :: q, n, j, a, c
            load = 0; mean = 0
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 2, t, q, x, wq, v2, status)
                if (status /= 0) return
                do n = 0, modes
                    bn(:, n) = matmul(v2, field%b(field%fs%sp(2)%dof(:, t), n) &
                        *field%fs%sp(2)%sgn(:, t))
                end do
                call meridian_qpoint(field%fs, 0, t, q, x, wq, v0, status)
                if (status /= 0) return
                psi = dot_product(v0(1, :), field%psi(field%fs%sp(0)%dof(:, t)))
                call axisymmetric_pressure(psi - field%psi_bound, field%pprime, p)
                call flux_mesh_geometry(field%mesh, t, xi(q), eta(q), x, jac)
                det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
                grad(1, :) = (jac(2, 2)*field%fs%tlg(1, :, q) &
                    - jac(2, 1)*field%fs%tlg(2, :, q))/det
                grad(2, :) = (-jac(1, 2)*field%fs%tlg(1, :, q) &
                    + jac(1, 1)*field%fs%tlg(2, :, q))/det
                b0 = field%f_vac/x(1)
                average = 0
                mean(1) = mean(1) + wq
                mean(3) = mean(3) + wq*abs(p)
                do j = 1, nphi
                    b = real(bn(:, 0), dp)
                    do n = 1, modes
                        b = b + 2*real(bn(:, n)*phases(n, j), dp)
                    end do
                    mean(2) = mean(2) + wq*sum(b*b)/nphi
                    call maxwell_stress_delta(b(1), b(2), b(3), b0, p, &
                        field%mu0, tau(1), tau(2), tau(3), tau(4), tau(5), &
                        tau(6), tau(7), tau(8), tau(9))
                    tensor = reshape(tau, [3, 3])
                    do n = 0, test_modes
                        average(:, :, n) = average(:, :, n) &
                            + tensor*conjg(phases(n, j))/nphi
                    end do
                end do
                do n = 0, test_modes
                    do a = 1, field%fs%nl0
                        do c = 1, 3
                            load(a, c, n) = load(a, c, n) + wq &
                                *(grad(1, a)*average(c, 1, n) &
                                + grad(2, a)*average(c, 3, n) &
                                - cmplx(0.0_dp, real(n, dp), dp)/x(1) &
                                *v0(1, a)*average(c, 2, n))
                        end do
                        load(a, 1, n) = load(a, 1, n) &
                            + wq*v0(1, a)*average(2, 2, n)/x(1)
                        load(a, 2, n) = load(a, 2, n) &
                            - wq*v0(1, a)*average(1, 2, n)/x(1)
                    end do
                end do
            end do
            status = 0
        end subroutine
    end subroutine
end module kin6d_nonnested_residual
