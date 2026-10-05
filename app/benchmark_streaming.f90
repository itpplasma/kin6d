program benchmark_streaming
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_leaves, only: coordinate, initial_coefficients, stream_step, &
                            stream_factors, stream_cached
    implicit none
    integer, parameter :: n = 131073, steps = 200, repeats = 5
    real(dp), allocatable :: v(:), m(:), a(:), b(:), an(:), bn(:), cv(:), sv(:)
    real(dp) :: start, finish, timings(repeats, 2), checksum(2), scratch
    integer :: mode, trial, j, k
    allocate (v(n), m(n), a(n), b(n), an(n), bn(n), cv(n), sv(n))
    do j = 1, n
        call coordinate(-4.0_dp, 4.0_dp, real(j - 1, dp), real(n - 1, dp), v(j))
    end do
    call stream_factors(v, 0.0625_dp, cv, sv)
    ! Alternate A/B order to reduce drift; include copies in both measurements.
    do trial = 1, repeats
        do k = 1, 2
            mode = k
            if (mod(trial, 2) == 0) mode = 3 - k
            call initial_coefficients(v, 0.5_dp, m, a, b)
            call cpu_time(start)
            do j = 1, steps
                if (mode == 1) then
                    call stream_step(a, b, v, 0.0625_dp, an, bn)
                else
                    call stream_cached(a, b, cv, sv, an, bn)
                end if
                a = an
                b = bn
            end do
            call cpu_time(finish)
            timings(trial, mode) = finish - start
            print '(a,i0,a,i0,a,f10.6)', 'trial=', trial, ',variant=', mode, &
                ',cpu_seconds=', timings(trial, mode)
            checksum(mode) = sum(a) + sum(b)
        end do
    end do
    if (abs(checksum(1) - checksum(2)) > 1.0e-10_dp) &
        error stop 'cached action differs from reference'
    do mode = 1, 2
        do j = 2, repeats
            do k = j, 2, -1
                if (timings(k, mode) >= timings(k - 1, mode)) exit
                scratch = timings(k, mode)
                timings(k, mode) = timings(k - 1, mode)
                timings(k - 1, mode) = scratch
            end do
        end do
    end do
    print '(a)', 'variant,median_cpu_seconds,min_cpu_seconds,max_cpu_seconds'
    print '(a,3(a,f10.6))', 'uncached', ',', timings(3, 1), &
        ',', timings(1, 1), ',', timings(5, 1)
    print '(a,3(a,f10.6))', 'cached', ',', timings(3, 2), &
        ',', timings(1, 2), ',', timings(5, 2)
    print '(a,f8.3)', 'measured_speedup=', timings(3, 1)/timings(3, 2)
end program
