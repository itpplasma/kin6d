program benchmark_nested_oracles
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_nested_oracles, only: nested_oracle_t
    use omp_lib, only: omp_get_wtime, omp_set_num_threads, &
        omp_set_dynamic, omp_get_num_threads
    implicit none
    type(nested_oracle_t) :: o
    real(dp), allocatable :: x(:, :), v(:, :), reference(:, :)
    integer, allocatable :: status(:)
    real(dp) :: t0, elapsed, jac, checksum
    integer :: n, k, threads, j, family, repeats, actual_threads
    character(32) :: arg
    call omp_set_dynamic(.false.)
    n = 1024; repeats = 3
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg); read (arg, *) n
    end if
    if (n < 1) error stop 'benchmark_nested_oracles: sample count'
    allocate (x(3, n), v(20, n), reference(20, n), status(n))
    do family = 1, 2
        o%family = family; o%eps = 0.5_dp; o%tau = 0.4_dp
        do k = 1, n
            call o%surface(0.1_dp, real(k, dp)*0.17_dp, real(k, dp)*0.29_dp, &
                           x(:, k), jac)
        end do
        call omp_set_num_threads(1)
        call o%sample(x, reference, status)
        if (any(status /= 0)) error stop 'benchmark reference admission'
        do j = 0, 3
            threads = 2**j
            call omp_set_num_threads(threads)
            !$omp parallel
            !$omp single
            actual_threads = omp_get_num_threads()
            !$omp end single
            !$omp end parallel
            if(actual_threads/=threads) error stop 'requested thread team unavailable'
            t0 = omp_get_wtime()
            do k = 1, repeats
                call o%sample(x, v, status)
            end do
            elapsed = (omp_get_wtime() - t0)/repeats
            if (any(status /= 0)) error stop 'benchmark sample admission'
          if (maxval(abs(v - reference)) > 1e-13_dp) error stop 'thread result mismatch'
            checksum = sum(v(3:5, :))
            print '(a,3i10,2es22.12)', '[V] family,points,threads,seconds,checksum ', &
                family, n, threads, elapsed, checksum
        end do
    end do
end program
