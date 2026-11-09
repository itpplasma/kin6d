program test_island_direct_coil_split
    ! Splits the FEEC vacuum-field error at the TIAGO island sample points
    ! into its spatial part (against the phi-truncated Biot-Savart field)
    ! and its Fourier-truncation part (Biot-Savart minus that truncation).
    ! The truncation is an independent DFT of the direct segment field, not
    ! the coil_set_modes boundary data used by the FEEC solve.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_vacuum_fixture, only: pi, field, make_coil, direct_field, solve_vacuum
    use island_a_fixture, only: island_a_sample_points
    use milca_equilibrium_nonnested, only: nn_point
    use milca_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI
    implicit none
    integer, parameter :: nq = 128, nmax = 3
    ! TIAGO's assigned polarization-angle uncertainty in the island benchmark.
    real(dp), parameter :: sigma = 1.0e-4_dp
    real(dp), allocatable :: xyz(:, :), rzp(:, :), exact(:, :), harmonic(:, :, :)
    real(dp) :: spatial, fourier, coarse, fine, pitch, limit, rate
    integer :: i
    call make_coil()
    xyz = island_a_sample_points()
    allocate (rzp(3, size(xyz, 2)), exact(3, size(xyz, 2)), &
        harmonic(3, size(xyz, 2), 0:nmax))
    do i = 1, size(xyz, 2)
        rzp(:, i) = [norm2(xyz(:2, i)), xyz(3, i), atan2(xyz(2, i), xyz(1, i))]
        call direct_field(rzp(:2, i), rzp(3, i), exact(:, i))
        call phi_harmonics(rzp(:, i), harmonic(:, i, :))
    end do

    ! The cross-model TIAGO setting (P2, nr=4, one mode) is dominated by the
    ! spatial error, not by the coil harmonics n >= 2 it omits.
    call solve_vacuum(4, 2, 1)
    call errors(1, spatial, fourier, pitch, limit)
    print '(a,4es12.4)', 'P2 nr=4 N=1 spatial/Fourier [T], pitch/limit [rad] ', &
        spatial, fourier, pitch, limit
    if (spatial < 10*fourier) error stop 'coarse error not spatial'
    if (pitch < 10*sigma) error stop 'coarse pitch error unexpectedly small'
    ! One mode cannot meet sigma at any spatial resolution.
    if (limit < sigma) error stop 'one-mode truncation below sigma'

    ! With three modes the FEEC field converges to the truncated field at P3
    ! order and resolves the direct field below sigma at nr=8.
    call solve_vacuum(4, 3, nmax)
    call errors(nmax, coarse, fourier, pitch, limit)
    call solve_vacuum(8, 3, nmax)
    call errors(nmax, fine, fourier, pitch, limit)
    rate = log(coarse/fine)/log(2.0_dp)
    print '(a,5es12.4)', 'P3 N=3 nr=4/8 spatial, Fourier, pitch, order ', &
        coarse, fine, fourier, pitch, rate
    if (rate < 3) error stop 'P3 spatial order'
    if (pitch > sigma) error stop 'resolved pitch error'
    print '(a)', 'PASS: direct-coil discrepancy is spatial and converges'
contains
    subroutine phi_harmonics(p, h)
        ! Cosine-shifted DFT about p(3): h(:, n) is harmonic n evaluated at p.
        real(dp), intent(in) :: p(3)
        real(dp), intent(out) :: h(:, 0:)
        real(dp) :: b(3), phi
        integer :: j, n
        h = 0
        do j = 0, nq - 1
            phi = 2*pi*j/nq
            call direct_field(p(:2), phi, b)
            do n = 0, ubound(h, 2)
                h(:, n) = h(:, n) + b*cos(n*(p(3) - phi))*merge(1, 2, n == 0)/nq
            end do
        end do
    end subroutine
    subroutine errors(modes, spatial, fourier, pitch, limit)
        integer, intent(in) :: modes
        real(dp), intent(out) :: spatial, fourier, pitch, limit
        real(dp) :: val(EQ_VALUES), fe(3), truncated(3)
        integer :: i, st
        spatial = 0; fourier = 0; pitch = 0; limit = 0
        do i = 1, size(rzp, 2)
            call nn_point(field, rzp(:, i), val, st)
            if (st /= 0) error stop 'sample outside the FEEC annulus'
            fe = [val(EQ_BR), val(EQ_BPHI), val(EQ_BZ)]
            truncated = sum(harmonic(:, i, 0:modes), 2)
            spatial = max(spatial, norm2(fe - truncated))
            fourier = max(fourier, norm2(exact(:, i) - truncated))
            pitch = max(pitch, &
                abs(atan2(fe(3), fe(2)) - atan2(exact(3, i), exact(2, i))), &
                abs(atan2(fe(1), fe(2)) - atan2(exact(1, i), exact(2, i))))
            limit = max(limit, abs(atan2(truncated(3), truncated(2)) &
                - atan2(exact(3, i), exact(2, i))), abs(atan2(truncated(1), &
                truncated(2)) - atan2(exact(1, i), exact(2, i))))
        end do
    end subroutine
end program
