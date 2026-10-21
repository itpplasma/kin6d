program benchmark_equilibrium_gpu
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use omp_lib, only: omp_get_max_threads
    use ideal_map_gpu_cpu, only: ideal_hessian_f64, ideal_hvp_f64
    use kin6d_ideal_map_leaves, only: ideal_energy_hessian
    implicit none
    interface
        integer(c_int) function gpu(n, repeats, precision, action, x, y, metrics) &
            bind(c, name='kin6d_gpu_benchmark')
            import c_int, c_double
            integer(c_int), value :: n, repeats, precision, action
            real(c_double), intent(in) :: x(*)
            real(c_double), intent(out) :: y(*), metrics(7)
        end function
    end interface
    real(dp), allocatable :: x(:, :), reference(:, :), candidate(:, :), hessian(:, :)
    real(dp), allocatable :: samples(:)
    real(dp) :: metrics(7), error, tolerance, t, h(9, 9), v(9), check_value
    real(dp) :: step, a(3, 3), plus(3, 3), minus(3, 3), fd, fd_error
    integer :: n, repeats, i, j, r, action, precision, outputs, st
    character(32) :: arg
    character(5), parameter :: names(0:2) = [character(5) :: 'fp64', 'fp32', 'mixed']
    n = 32768; repeats = 21
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read (arg, *) n
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read (arg, *) repeats
    end if
    if (n < 1 .or. n > 1048576 .or. repeats < 1 .or. repeats > 101) &
        error stop 'require 1 <= n <= 1048576, 1 <= repeats <= 101'
    allocate (x(n, 24), reference(n, 45), candidate(n, 45), hessian(n, 45))
    allocate (samples(repeats))
    ! Deterministic well-conditioned deformations, not a solved equilibrium.
    ! SoA (sample, component) matches Fortran storage and CUDA coalescing.
    do j = 1, 9
        do i = 1, n
            x(i, j) = 0.03_dp*sin(real(13*i + 7*j, dp))
            x(i, 15 + j) = sin(real(3*i + 11*j, dp))/3
        end do
    end do
    x(:, 1) = x(:, 1) + 1; x(:, 5) = x(:, 5) + 1; x(:, 9) = x(:, 9) + 1
    x(:, 10) = 0.2_dp; x(:, 11) = 0.1_dp; x(:, 12) = 1
    x(:, 13) = 0.125_dp; x(:, 14) = 1; x(:, 15) = 1
    call cpu(x, hessian, n, 0, .true.)
    ! Independent scalar energy directional curvature, including det(DX).
    fd_error = 0; step = 1e-4_dp
    do i = 1, min(n, 31)
        call unpack(hessian(i, :), h)
        v = x(i, 16:24)
        a = reshape(x(i, 1:9), [3, 3])
        plus = a + step*reshape(v, [3, 3]); minus = a - step*reshape(v, [3, 3])
        fd = (energy(plus, x(i, 10:12), x(i, 13)) &
            - 2*energy(a, x(i, 10:12), x(i, 13)) &
            + energy(minus, x(i, 10:12), x(i, 13)))/step**2
        check_value = dot_product(v, matmul(h, v))
        fd_error = max(fd_error, abs(fd-check_value)/max(1.0_dp, abs(check_value)))
    end do
    if (fd_error > 3e-6_dp) error stop 'independent energy curvature oracle'
    print '(a,i0,a,i0,a,i0,a,es12.4)', '[V] n=', n, ' repeats=', repeats, &
        ' threads=', omp_get_max_threads(), ' energy_fd_error=', fd_error
    print '(a)', '[V] claim=candidate; enclosure=unavailable; no solver/estimator replaced'
    print '(a)', 'action precision cpu_min_s cpu_median_s cpu_max_s gpu_min_s ' // &
        'gpu_median_s gpu_max_s h2d_s d2h_s cold_total_s device_bytes max_scaled_error'
    do action = 0, 1
        outputs = 45
        if (action == 1) outputs = 9
        if (action == 0) reference = hessian
        if (action == 1) then
            ! Independent H*v oracle using the production assembled Hessian.
            do i = 1, n
                call unpack(hessian(i, :), h)
                reference(i, 1:9) = matmul(h, x(i, 16:24))
            end do
        end if
        call cpu(x, candidate, n, action, .false.)
        if (scaled_error(candidate, reference, outputs) > 2e-12_dp) &
            error stop 'generated CPU vs production Hessian'
        do r = 1, repeats
            t = wall()
            call cpu(x, candidate, n, action, .false.)
            samples(r) = wall()-t
        end do
        call sort(samples)
        do precision = 0, 2
            st = gpu(int(n, c_int), int(repeats, c_int), int(precision, c_int), &
                int(action, c_int), x, candidate, metrics)
            if (st /= 0) error stop 'CUDA benchmark failed'
            error = scaled_error(candidate, reference, outputs)
            tolerance = 3e-5_dp
            if (precision == 0) tolerance = 2e-12_dp
            if (.not. all(ieee_is_finite(candidate(:, :outputs)))) &
                error stop 'nonfinite CUDA output'
            if (error > tolerance) then
                print *, 'action, precision, error, tolerance', action, precision, error, tolerance
                error stop 'CUDA vs CPU mismatch'
            end if
            print '(i1,1x,a5,11(1x,es14.6))', action, names(precision), &
                samples(1), samples(1+repeats/2), samples(repeats), &
                metrics(3:5), metrics(2), metrics(6), metrics(1), metrics(7), error
        end do
    end do
contains
    subroutine cpu(x, y, n, action, production)
        integer, intent(in) :: n, action
        real(dp), intent(in) :: x(n, 24)
        real(dp), intent(out) :: y(n, 45)
        logical, intent(in) :: production
        integer :: i
        !$omp parallel do default(none) shared(n, x, y, action, production) private(i)
        do i = 1, n
            if (action == 0) then
                if (production) then
                    call ideal_energy_hessian( &
                        x(i, 1), x(i, 2), x(i, 3), x(i, 4), &
                        x(i, 5), x(i, 6), x(i, 7), x(i, 8), &
                        x(i, 9), x(i, 10), x(i, 11), x(i, 12), &
                        x(i, 13), x(i, 14), x(i, 15), y(i, 1), &
                        y(i, 2), y(i, 3), y(i, 4), y(i, 5), &
                        y(i, 6), y(i, 7), y(i, 8), y(i, 9), &
                        y(i, 10), y(i, 11), y(i, 12), y(i, 13), &
                        y(i, 14), y(i, 15), y(i, 16), y(i, 17), &
                        y(i, 18), y(i, 19), y(i, 20), y(i, 21), &
                        y(i, 22), y(i, 23), y(i, 24), y(i, 25), &
                        y(i, 26), y(i, 27), y(i, 28), y(i, 29), &
                        y(i, 30), y(i, 31), y(i, 32), y(i, 33), &
                        y(i, 34), y(i, 35), y(i, 36), y(i, 37), &
                        y(i, 38), y(i, 39), y(i, 40), y(i, 41), &
                        y(i, 42), y(i, 43), y(i, 44), y(i, 45))
                else
                    call ideal_hessian_f64( &
                        x(i, 1), x(i, 2), x(i, 3), x(i, 4), &
                        x(i, 5), x(i, 6), x(i, 7), x(i, 8), &
                        x(i, 9), x(i, 10), x(i, 11), x(i, 12), &
                        x(i, 13), x(i, 14), x(i, 15), y(i, 1), &
                        y(i, 2), y(i, 3), y(i, 4), y(i, 5), &
                        y(i, 6), y(i, 7), y(i, 8), y(i, 9), &
                        y(i, 10), y(i, 11), y(i, 12), y(i, 13), &
                        y(i, 14), y(i, 15), y(i, 16), y(i, 17), &
                        y(i, 18), y(i, 19), y(i, 20), y(i, 21), &
                        y(i, 22), y(i, 23), y(i, 24), y(i, 25), &
                        y(i, 26), y(i, 27), y(i, 28), y(i, 29), &
                        y(i, 30), y(i, 31), y(i, 32), y(i, 33), &
                        y(i, 34), y(i, 35), y(i, 36), y(i, 37), &
                        y(i, 38), y(i, 39), y(i, 40), y(i, 41), &
                        y(i, 42), y(i, 43), y(i, 44), y(i, 45))
                end if
            else
                call ideal_hvp_f64( &
                    x(i, 1), x(i, 2), x(i, 3), x(i, 4), &
                    x(i, 5), x(i, 6), x(i, 7), x(i, 8), &
                    x(i, 9), x(i, 10), x(i, 11), x(i, 12), &
                    x(i, 13), x(i, 14), x(i, 15), x(i, 16), &
                    x(i, 17), x(i, 18), x(i, 19), x(i, 20), &
                    x(i, 21), x(i, 22), x(i, 23), x(i, 24), &
                    y(i, 1), y(i, 2), y(i, 3), y(i, 4), &
                    y(i, 5), y(i, 6), y(i, 7), y(i, 8), &
                    y(i, 9))
            end if
        end do
        !$omp end parallel do
    end subroutine

    real(dp) function scaled_error(a, b, outputs) result(value)
        real(dp), intent(in) :: a(:, :), b(:, :)
        integer, intent(in) :: outputs
        integer :: i
        value = 0
        do i = 1, size(a, 1)
            value = max(value, maxval(abs(a(i, :outputs)-b(i, :outputs))) &
                /max(1.0_dp, maxval(abs(b(i, :outputs)))))
        end do
    end function

    subroutine unpack(flat, h)
        real(dp), intent(in) :: flat(45)
        real(dp), intent(out) :: h(9, 9)
        integer :: i, j, k
        k = 0
        do j = 1, 9
            do i = 1, j
                k = k + 1
                h(i, j) = flat(k); h(j, i) = flat(k)
            end do
        end do
    end subroutine

    real(dp) function energy(d, b, p) result(w)
        ! Independent oracle only, mu0=sigma=1 in this fixture.
        real(dp), intent(in) :: d(3, 3), b(3), p
        real(dp) :: determinant, field(3)
        determinant = d(1, 1)*(d(2, 2)*d(3, 3)-d(2, 3)*d(3, 2)) &
            - d(1, 2)*(d(2, 1)*d(3, 3)-d(2, 3)*d(3, 1)) &
            + d(1, 3)*(d(2, 1)*d(3, 2)-d(2, 2)*d(3, 1))
        field = matmul(d, b)
        w = dot_product(field, field)/(2*determinant)-p*determinant
    end function

    real(dp) function wall() result(t)
        integer(int64) :: ticks, rate
        call system_clock(ticks, rate)
        t = real(ticks, dp)/real(rate, dp)
    end function

    subroutine sort(x)
        real(dp), intent(inout) :: x(:)
        real(dp) :: t
        integer :: i, j
        do i = 2, size(x)
            t = x(i); j = i-1
            do while (j >= 1)
                if (x(j) <= t) exit
                x(j+1) = x(j); j = j-1
            end do
            x(j+1) = t
        end do
    end subroutine
end program benchmark_equilibrium_gpu
