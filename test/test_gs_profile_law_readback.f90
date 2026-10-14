program test_gs_profile_law_readback
    !! Readback of an executed profile law written by gs_phase1 (6c81202, MUMPS)
    !! on the TC24 s_pol 0.997 boundary: the first three intervals of its
    !! profile-law.dat, digits unchanged. At knot 2 p' is 139 while its cubic
    !! terms are 2.2e5, and the written coefficients leave a knot defect of 173
    !! eps relative to the value. A value-scaled continuity check rejected the
    !! law, which stopped every gs_phase1 --evaluate readback of that solve.
    !! The law must read back and survive a write/read cycle unchanged.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_gs_profile_law, only: gs_profile_law_t, read_gs_profile_law, &
                                    write_gs_profile_law
    implicit none
    character(*), parameter :: path = 'test_gs_profile_law_readback.dat'
    character(*), parameter :: copy = 'test_gs_profile_law_readback_copy.dat'
    type(gs_profile_law_t) :: law, again
    integer :: unit, status

    open (newunit=unit, file=path, status='replace', action='write')
    write (unit, '(a)') '4 0'
    write (unit, '(a)') '6.4001140918168367E+003 -3.2854944813396052E+001 ' &
        //'-1.1850683305509241E+001'
    write (unit, '(a)') '0.0 1.9531250000000000E-003 3.9062500000000000E-003 ' &
        //'5.8593750000000000E-003'
    write (unit, '(a)') '1.2271447943921060E+005 2.0501720513996076E+007 ' &
        //'-1.7218473337505180E+010 1.0078545124683417E+012 ' &
        //'-1.7676938357972592E+001 -1.0656206682487825E+003 ' &
        //'8.9163468467457662E+005 -5.2231379323230661E+007'
    write (unit, '(a)') '1.0458274043440996E+005 -3.5223961402051128E+007 ' &
        //'-1.1313075803510994E+010 1.0078545124683456E+012 ' &
        //'-1.6746066433417710E+001 1.8195866191088323E+003 ' &
        //'5.8559144645252207E+005 -5.2231379323230661E+007'
    write (unit, '(a)') '1.3908276268589009E+002 -6.7881684200933501E+007 ' &
        //'-5.4076782695167809E+009 2.7511292762096826E+012 ' &
        //'-1.1347486329751179E+001 3.5093125071615459E+003 ' &
        //'2.7954820823046623E+005 -1.4224831528998795E+008'
    close (unit)

    call read_gs_profile_law(path, law, status)
    if (status /= 0) error stop 'executed 0.997 law rejected by its reader'
    call write_gs_profile_law(copy, law)
    call read_gs_profile_law(copy, again, status)
    if (status /= 0) error stop 'rewritten law rejected by its reader'
    if (any(again%knots /= law%knots) .or. any(again%pp /= law%pp) &
        .or. any(again%ff /= law%ff)) error stop 'law changed in write/read'
    if (again%p_edge /= law%p_edge .or. again%f_edge /= law%f_edge &
        .or. again%axis /= law%axis) error stop 'law edge data changed'
    open (newunit=unit, file=path, status='old')
    close (unit, status='delete')
    open (newunit=unit, file=copy, status='old')
    close (unit, status='delete')
    print '(a)', 'executed profile law readback: PASS'
end program
