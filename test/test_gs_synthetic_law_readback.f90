program test_gs_synthetic_law_readback
    ! A synthetic zero crossing has a tiny rounded knot mismatch relative to
    ! cancelling million-sized terms. A value-scaled continuity check fails.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_gs_profile_law, only: gs_profile_law_t, read_gs_profile_law, &
        write_gs_profile_law
    implicit none
    type(gs_profile_law_t) :: law, again
    character(*), parameter :: path = 'synthetic-law.dat', copy = 'synthetic-copy.dat'
    integer :: unit, status

    allocate (law%knots(3), law%pp(4, 2), law%ff(4, 2))
    law%normalized = .false.
    law%p_edge = 1; law%f_edge = 2; law%axis = -1
    law%knots = [0.0_dp, 1.0_dp, 2.0_dp]
    law%pp(:, 1) = [1.0e6_dp, -2.0e6_dp, 1.0e6_dp, 0.0_dp]
    law%pp(:, 2) = [4*spacing(1.0e6_dp), 0.0_dp, 0.0_dp, 0.0_dp]
    law%ff = 0
    call write_gs_profile_law(path, law)
    call read_gs_profile_law(path, again, status)
    if (status /= 0) error stop 'synthetic cancellation law rejected'
    call write_gs_profile_law(copy, again)
    call read_gs_profile_law(copy, again, status)
    if (status /= 0) error stop 'rewritten synthetic law rejected'
    if (any(again%knots /= law%knots) .or. any(again%pp /= law%pp) &
        .or. any(again%ff /= law%ff)) error stop 'law changed in write/read'
    if (again%p_edge /= law%p_edge .or. again%f_edge /= law%f_edge &
        .or. again%axis /= law%axis) error stop 'law edge data changed'
    law%pp(1, 2) = 1.0e-3_dp
    call write_gs_profile_law(path, law)
    call read_gs_profile_law(path, again, status)
    if (status == 0) error stop 'discontinuous synthetic law accepted'
    open (newunit=unit, file=path, status='old')
    close (unit, status='delete')
    open (newunit=unit, file=copy, status='old')
    close (unit, status='delete')
    print '(a)', 'synthetic profile law readback: PASS'
end program
