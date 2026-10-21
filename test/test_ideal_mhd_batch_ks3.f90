program test_ideal_mhd_batch_ks3
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32
    use kin6d_ideal_xi, only: xi_bilinear_batch
    use kin6d_ideal_xi_leaves, only: xi_form_dp, xi_form_jvp_dp
    implicit none
    integer, parameter :: n = 257
    real(dp) :: x(n, 47), y(n, 25), z(n, 25), value, tangent, expected
    real(sp) :: xs(n, 47), ys(n, 25)
    integer :: j, k
    x = 0; x(:, 3) = 1; x(:, 5) = 1; x(:, 47) = 1
    do j = 1, n
        x(j, 12) = 0.25_dp + real(mod(j, 97), dp)/128
        x(j, 21) = -0.125_dp + real(mod(j, 53), dp)/64
        x(j, 35) = 0.125_dp; x(j, 44) = -0.25_dp
    end do
    call xi_bilinear_batch(n, x, y)
    do j = 1, n
        expected = x(j, 12)*x(j, 21)
        if (abs(y(j, 1) - expected) > 2e-12_dp) error stop 'KS3 analytical slab'
        call xi_form_dp(x(j, 1), x(j, 2), x(j, 3), x(j, 4), x(j, 5), x(j, 6), &
                        x(j, 7), x(j, 8), x(j, 9), x(j, 10), x(j, 11), x(j, 12), &
                        x(j, 13), x(j, 14), x(j, 15), x(j, 16), x(j, 17), x(j, 18), &
                        x(j, 19), x(j, 20), x(j, 21), x(j, 22), x(j, 23), value)
        if (abs(y(j, 1) - value) > 2e-12_dp) error stop "KS3 scalar/batch value"
        call xi_form_jvp_dp(x(j, 1), x(j, 2), x(j, 3), x(j, 4), x(j, 5), x(j, 6), &
                            x(j, 7), x(j, 8), x(j, 9), x(j, 10), x(j, 11), x(j, 12), &
                           x(j, 13), x(j, 14), x(j, 15), x(j, 16), x(j, 17), x(j, 18), &
                           x(j, 19), x(j, 20), x(j, 21), x(j, 22), x(j, 23), x(j, 24), &
                           x(j, 25), x(j, 26), x(j, 27), x(j, 28), x(j, 29), x(j, 30), &
                           x(j, 31), x(j, 32), x(j, 33), x(j, 34), x(j, 35), x(j, 36), &
                           x(j, 37), x(j, 38), x(j, 39), x(j, 40), x(j, 41), x(j, 42), &
                            x(j, 43), x(j, 44), x(j, 45), x(j, 46), tangent)
        if (abs(y(j, 2) - tangent) > 2e-12_dp) error stop 'KS3 scalar/batch JVP'
        if (abs(sum(y(j, 3:)*x(j, 24:46)) - tangent) > 2e-12_dp) &
            error stop 'KS3 VJP pairing'
    end do
    xs = real(x, sp)
    call xi_bilinear_batch(n, xs, ys)
    call xi_bilinear_batch(n, xs, z)
    if (maxval(abs(real(ys, dp) - y)) > 512*epsilon(1.0_sp)) error stop 'KS3 fp32'
    if (maxval(abs(z - y)) > 512*epsilon(1.0_sp)) error stop 'KS3 mixed'
    print *, 'PASS KS3 shared batch/scalar/analytical/derivative/kind oracles'
end program test_ideal_mhd_batch_ks3
