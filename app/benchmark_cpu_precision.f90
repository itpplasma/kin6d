program benchmark_cpu_precision
    ! Branch-only precision seam: upstream-generated leaves, fixed table reuse.
    ! This measures arithmetic and storage together, excluding one-time casts.
    ! It does not select a production preset or establish a continuous bound.
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32, int64, int16
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortgen_storage, only: pack_fp16, unpack_fp16, pack_bf16, unpack_bf16
    use omp_lib, only: omp_get_max_threads
    use gs_cpu_fp64_batch_mod, only: gs_cpu_fp64_batch
    use gs_cpu_fp32_batch_mod, only: gs_cpu_fp32_batch
    use map_cpu_fp64_batch_mod, only: map_cpu_fp64_batch
    use map_cpu_fp32_batch_mod, only: map_cpu_fp32_batch
    implicit none
    real(dp), allocatable :: x(:, :), y(:, :), oracle(:, :)
    real(sp), allocatable :: xs(:, :), ys(:, :)
    integer(int16), allocatable :: hx(:, :), hy(:, :)
    real(dp) :: best64, best32, t, cast_s, defect, scale, check_value
    integer :: n, repeats, i, j, rep, nx, ny, format
    real(dp) :: half_s, half_cast, half_defect
    character(32) :: mode, arg
    mode = 'gs'; n = 131072; repeats = 20
    call get_command_argument(1, mode)
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg); read (arg, *) n
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg); read (arg, *) repeats
    end if
    if (n < 1 .or. repeats < 1) error stop 'positive size and repeats required'
    select case (trim(mode))
    case ('gs')
        nx = 7; ny = 3
    case ('map')
        nx = 13; ny = 45
    case default
        error stop 'mode must be gs or map'
    end select
    allocate (x(n, nx), xs(n, nx), y(n, ny), ys(n, ny), oracle(n, ny))
    do j = 1, nx
        do i = 1, n
            x(i, j) = 0.03_dp*sin(real(i + 17*j, dp))
        end do
    end do
    if (mode == 'gs') then
        x(:, 1) = 3 + x(:, 1)
        x(:, 2:) = x(:, 2:) - 0.2_dp
        ! Independent analytic source/tangent oracle (not used in timing).
        do j = 1, 3
            oracle(:, j) = x(:, 1)*x(:, 2*j) + x(:, 2*j + 1)/x(:, 1)
        end do
    else
        x(:, 1) = 1 + x(:, 1); x(:, 5) = 1 + x(:, 5)
        x(:, 9) = 1 + x(:, 9)
        x(:, 10:12) = x(:, 10:12) + 0.5_dp
        x(:, 13) = x(:, 13) + 0.1_dp
        ! Independent identity-map oracle with b=(0,0,1), mu0=1.
        ! w = |D b|^2/(2 det D) - p det D gives H11=H77=1, H99=0.
        x(n, :) = 0
        x(n, 1) = 1; x(n, 5) = 1; x(n, 9) = 1
        x(n, 12) = 1; x(n, 13) = 0.1_dp
    end if
    t = wall(); xs = real(x, sp); cast_s = wall() - t
    call run64(); call run32()
    if (.not. all(ieee_is_finite(y))) error stop 'nonfinite fp64 output'
    if (.not. all(ieee_is_finite(ys))) error stop 'nonfinite fp32 output'
    if (mode == 'gs') then
        if (maxval(abs(y - oracle)) > 1e-12_dp) error stop 'GS analytic oracle'
    else
        if (abs(y(n, 1) - 1) > 1e-12_dp .or. &
            abs(y(n, 28) - 1) > 1e-12_dp .or. &
            abs(y(n, 45)) > 1e-12_dp) error stop 'map analytic Hessian oracle'
    end if
    scale = maxval(abs(y))
    defect = maxval(abs(real(ys, dp) - y))/max(scale, tiny(scale))
    ! An empirical comparison on this declared well-conditioned domain.
    ! No rounding majorant, automatic escalation or certificate is claimed.
    if (defect > 1e-2_dp) error stop 'fp32 kernel defect exceeds 1%'
    best64 = huge(1.0_dp); best32 = huge(1.0_dp)
    do rep = 1, repeats
        ! Alternate order to limit systematic cache/clock bias.
        if (mod(rep, 2) == 1) then
            t = wall(); call run64(); best64 = min(best64, wall() - t)
            t = wall(); call run32(); best32 = min(best32, wall() - t)
        else
            t = wall(); call run32(); best32 = min(best32, wall() - t)
            t = wall(); call run64(); best64 = min(best64, wall() - t)
        end if
    end do
    check_value = sum(y)
    print '(a)', '# [V] microkernel comparison; continuous certificate unavailable'
    print '(a)', '# mode,n,threads,fp64_s,candidate_s,speedup,setup_s,'// &
        'scaled_defect,checksum'
    print '(a,a,i0,a,i0,6(a,es23.15))', trim(mode)//'-fp32', ',', n, ',', &
        omp_get_max_threads(), ',', best64, ',', best32, ',', best64/best32, &
        ',', cast_s, ',', defect, ',', check_value
    ! The library codec is intentionally timed inside each tile. This is
    ! portable compressed storage with fp32 arithmetic, not native half math.
    allocate (hx(n, nx), hy(n, ny))
    do format = 1, 2
        t = wall()
        if (format == 1) then
            hx = pack_fp16(x)
        else
            hx = pack_bf16(x)
        end if
        half_cast = wall() - t
        call run16()
        if (format == 1) then
            oracle = unpack_fp16(hy)
        else
            oracle = unpack_bf16(hy)
        end if
        if (.not. all(ieee_is_finite(oracle))) error stop 'nonfinite half output'
        half_defect = maxval(abs(oracle - y))/max(scale, tiny(scale))
        if (half_defect > 3e-2_dp) error stop 'half defect exceeds 3% on fixture'
        half_s = huge(1.0_dp)
        do rep = 1, repeats
            t = wall(); call run16(); half_s = min(half_s, wall() - t)
        end do
        if (format == 1) then
            arg = trim(mode)//'-fp16'
        else
            arg = trim(mode)//'-bf16'
        end if
        print '(a,a,i0,a,i0,6(a,es23.15))', trim(arg), ',', n, ',', &
            omp_get_max_threads(), ',', best64, ',', half_s, ',', best64/half_s, &
            ',', half_cast, ',', half_defect, ',', sum(oracle)
    end do
contains
    function wall() result(t)
        real(dp) :: t
        integer(int64) :: count, rate
        call system_clock(count, rate)
        t = real(count, dp)/real(rate, dp)
    end function
    subroutine run64()
        integer :: first, last, count
        !$omp parallel do private(last, count) schedule(static)
        do first = 1, n, 1024
            last = min(n, first + 1023); count = last - first + 1
            if (mode == 'gs') then
                call gs_cpu_fp64_batch( &
                    count, &
                    x(first:last, 1), &
                    1.0_dp, &
                    x(first:last, 2), &
                    x(first:last, 3), &
                    x(first:last, 4), &
                    x(first:last, 5), &
                    x(first:last, 6), &
                    x(first:last, 7), &
                    y(first:last, 1), &
                    y(first:last, 2), &
                    y(first:last, 3))
            else
                call map_cpu_fp64_batch( &
                    count, &
                    x(first:last, 1), &
                    x(first:last, 2), &
                    x(first:last, 3), &
                    x(first:last, 4), &
                    x(first:last, 5), &
                    x(first:last, 6), &
                    x(first:last, 7), &
                    x(first:last, 8), &
                    x(first:last, 9), &
                    x(first:last, 10), &
                    x(first:last, 11), &
                    x(first:last, 12), &
                    x(first:last, 13), &
                    1.0_dp, &
                    1.0_dp, &
                    y(first:last, 1), &
                    y(first:last, 2), &
                    y(first:last, 3), &
                    y(first:last, 4), &
                    y(first:last, 5), &
                    y(first:last, 6), &
                    y(first:last, 7), &
                    y(first:last, 8), &
                    y(first:last, 9), &
                    y(first:last, 10), &
                    y(first:last, 11), &
                    y(first:last, 12), &
                    y(first:last, 13), &
                    y(first:last, 14), &
                    y(first:last, 15), &
                    y(first:last, 16), &
                    y(first:last, 17), &
                    y(first:last, 18), &
                    y(first:last, 19), &
                    y(first:last, 20), &
                    y(first:last, 21), &
                    y(first:last, 22), &
                    y(first:last, 23), &
                    y(first:last, 24), &
                    y(first:last, 25), &
                    y(first:last, 26), &
                    y(first:last, 27), &
                    y(first:last, 28), &
                    y(first:last, 29), &
                    y(first:last, 30), &
                    y(first:last, 31), &
                    y(first:last, 32), &
                    y(first:last, 33), &
                    y(first:last, 34), &
                    y(first:last, 35), &
                    y(first:last, 36), &
                    y(first:last, 37), &
                    y(first:last, 38), &
                    y(first:last, 39), &
                    y(first:last, 40), &
                    y(first:last, 41), &
                    y(first:last, 42), &
                    y(first:last, 43), &
                    y(first:last, 44), &
                    y(first:last, 45))
            end if
        end do
        !$omp end parallel do
    end subroutine
    subroutine run32()
        integer :: first, last, count
        !$omp parallel do private(last, count) schedule(static)
        do first = 1, n, 1024
            last = min(n, first + 1023); count = last - first + 1
            if (mode == 'gs') then
                call gs_cpu_fp32_batch( &
                    count, &
                    xs(first:last, 1), &
                    1.0_sp, &
                    xs(first:last, 2), &
                    xs(first:last, 3), &
                    xs(first:last, 4), &
                    xs(first:last, 5), &
                    xs(first:last, 6), &
                    xs(first:last, 7), &
                    ys(first:last, 1), &
                    ys(first:last, 2), &
                    ys(first:last, 3))
            else
                call map_cpu_fp32_batch( &
                    count, &
                    xs(first:last, 1), &
                    xs(first:last, 2), &
                    xs(first:last, 3), &
                    xs(first:last, 4), &
                    xs(first:last, 5), &
                    xs(first:last, 6), &
                    xs(first:last, 7), &
                    xs(first:last, 8), &
                    xs(first:last, 9), &
                    xs(first:last, 10), &
                    xs(first:last, 11), &
                    xs(first:last, 12), &
                    xs(first:last, 13), &
                    1.0_sp, &
                    1.0_sp, &
                    ys(first:last, 1), &
                    ys(first:last, 2), &
                    ys(first:last, 3), &
                    ys(first:last, 4), &
                    ys(first:last, 5), &
                    ys(first:last, 6), &
                    ys(first:last, 7), &
                    ys(first:last, 8), &
                    ys(first:last, 9), &
                    ys(first:last, 10), &
                    ys(first:last, 11), &
                    ys(first:last, 12), &
                    ys(first:last, 13), &
                    ys(first:last, 14), &
                    ys(first:last, 15), &
                    ys(first:last, 16), &
                    ys(first:last, 17), &
                    ys(first:last, 18), &
                    ys(first:last, 19), &
                    ys(first:last, 20), &
                    ys(first:last, 21), &
                    ys(first:last, 22), &
                    ys(first:last, 23), &
                    ys(first:last, 24), &
                    ys(first:last, 25), &
                    ys(first:last, 26), &
                    ys(first:last, 27), &
                    ys(first:last, 28), &
                    ys(first:last, 29), &
                    ys(first:last, 30), &
                    ys(first:last, 31), &
                    ys(first:last, 32), &
                    ys(first:last, 33), &
                    ys(first:last, 34), &
                    ys(first:last, 35), &
                    ys(first:last, 36), &
                    ys(first:last, 37), &
                    ys(first:last, 38), &
                    ys(first:last, 39), &
                    ys(first:last, 40), &
                    ys(first:last, 41), &
                    ys(first:last, 42), &
                    ys(first:last, 43), &
                    ys(first:last, 44), &
                    ys(first:last, 45))
            end if
        end do
        !$omp end parallel do
    end subroutine
    subroutine run16()
        ! Fixed caller-stack tile scratch: no allocation in the OpenMP loop.
        real(sp) :: tile_x(1024, 13), tile_y(1024, 45)
        integer :: first, last, count, column
        !$omp parallel do private(last, count, column, tile_x, tile_y) schedule(static)
        do first = 1, n, 1024
            last = min(n, first + 1023); count = last - first + 1
            do column = 1, nx
                if (format == 1) then
                    tile_x(:count, column) = &
                        real(unpack_fp16(hx(first:last, column)), sp)
                else
                    tile_x(:count, column) = &
                        real(unpack_bf16(hx(first:last, column)), sp)
                end if
            end do
            if (mode == 'gs') then
                call gs_cpu_fp32_batch( &
                    count, &
                    tile_x(:count, 1), &
                    1.0_sp, &
                    tile_x(:count, 2), &
                    tile_x(:count, 3), &
                    tile_x(:count, 4), &
                    tile_x(:count, 5), &
                    tile_x(:count, 6), &
                    tile_x(:count, 7), &
                    tile_y(:count, 1), &
                    tile_y(:count, 2), &
                    tile_y(:count, 3))
            else
                call map_cpu_fp32_batch( &
                    count, &
                    tile_x(:count, 1), &
                    tile_x(:count, 2), &
                    tile_x(:count, 3), &
                    tile_x(:count, 4), &
                    tile_x(:count, 5), &
                    tile_x(:count, 6), &
                    tile_x(:count, 7), &
                    tile_x(:count, 8), &
                    tile_x(:count, 9), &
                    tile_x(:count, 10), &
                    tile_x(:count, 11), &
                    tile_x(:count, 12), &
                    tile_x(:count, 13), &
                    1.0_sp, &
                    1.0_sp, &
                    tile_y(:count, 1), &
                    tile_y(:count, 2), &
                    tile_y(:count, 3), &
                    tile_y(:count, 4), &
                    tile_y(:count, 5), &
                    tile_y(:count, 6), &
                    tile_y(:count, 7), &
                    tile_y(:count, 8), &
                    tile_y(:count, 9), &
                    tile_y(:count, 10), &
                    tile_y(:count, 11), &
                    tile_y(:count, 12), &
                    tile_y(:count, 13), &
                    tile_y(:count, 14), &
                    tile_y(:count, 15), &
                    tile_y(:count, 16), &
                    tile_y(:count, 17), &
                    tile_y(:count, 18), &
                    tile_y(:count, 19), &
                    tile_y(:count, 20), &
                    tile_y(:count, 21), &
                    tile_y(:count, 22), &
                    tile_y(:count, 23), &
                    tile_y(:count, 24), &
                    tile_y(:count, 25), &
                    tile_y(:count, 26), &
                    tile_y(:count, 27), &
                    tile_y(:count, 28), &
                    tile_y(:count, 29), &
                    tile_y(:count, 30), &
                    tile_y(:count, 31), &
                    tile_y(:count, 32), &
                    tile_y(:count, 33), &
                    tile_y(:count, 34), &
                    tile_y(:count, 35), &
                    tile_y(:count, 36), &
                    tile_y(:count, 37), &
                    tile_y(:count, 38), &
                    tile_y(:count, 39), &
                    tile_y(:count, 40), &
                    tile_y(:count, 41), &
                    tile_y(:count, 42), &
                    tile_y(:count, 43), &
                    tile_y(:count, 44), &
                    tile_y(:count, 45))
            end if
            do column = 1, ny
                if (format == 1) then
                    hy(first:last, column) = pack_fp16(real(tile_y(:count, column), dp))
                else
                    hy(first:last, column) = pack_bf16(real(tile_y(:count, column), dp))
                end if
            end do
        end do
        !$omp end parallel do
    end subroutine
end program benchmark_cpu_precision
