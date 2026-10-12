program test_gs_coils
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_axisymmetric_gs_kernel, only: gs_mu0, gs_rect_coil_t, &
        gs_rect_coil_flux, gs_rect_coil_flux_converged
    use kin6d_gs_coils
    implicit none
    integer :: failures

    failures = 0
    call test_rect_flux()
    call test_freegs_parity()
    call test_series_incidence()
    call test_limits()
    call test_table_roundtrip()
    call test_table_errors()
    call test_helmholtz()
    if (failures /= 0) error stop 'gs_coils oracle failures'
    print '(a)', 'gs_coils oracles passed'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            failures = failures + 1
            print '(2a)', 'FAIL: ', what
        end if
    end subroutine require

    subroutine rel_check(a, b, tol, what)
        real(dp), intent(in) :: a, b, tol
        character(len=*), intent(in) :: what
        if (.not. abs(a - b) <= tol*abs(b)) then
            print '(a,2es24.16,es10.2)', what, a, b, abs(a - b)/abs(b)
            call require(.false., what)
        end if
    end subroutine rel_check

    function rect_coil(r, z, dr, dz, turns) result(c)
        real(dp), intent(in) :: r, z, dr, dz, turns
        type(gs_coil_t) :: c
        c%rect = gs_rect_coil_t(r, z, dr, dz, turns)
    end function rect_coil

    subroutine test_rect_flux()
        !! Quadrature coil against fortfem's converged rule and the filament limit.
        type(gs_coil_t) :: coil, fine, coarse
        real(dp), parameter :: pts(2, 4) = reshape([5.0_dp, 0.0_dp, 8.0_dp, -2.0_dp, &
            6.0_dp, 2.5_dp, 4.0_dp, 4.5_dp], [2, 4])
        real(dp) :: ref, e1, e2
        integer :: i, order, status

        coil = rect_coil(6.2_dp, 4.5_dp, 0.7_dp, 0.9_dp, 3.0_dp)
        coil%quad_order = 24
        fine = coil
        call gs_coil_make_freegs_filaments(fine, 0.02_dp)
        coarse = coil
        call gs_coil_make_freegs_filaments(coarse, 0.04_dp)
        do i = 1, 4
            call gs_rect_coil_flux_converged(coil%rect, pts(1, i), pts(2, i), &
                1.0e-13_dp, ref, order, status)
            call require(status == 0, 'fortfem converged status')
            call rel_check(gs_coil_flux(coil, pts(1, i), pts(2, i)), ref, 1.0e-11_dp, &
                'rect flux vs converged quadrature')
            ! The midpoint filament rule is second order: halving the
            ! filament size must cut the error by about four.
            e1 = abs(gs_coil_flux(coarse, pts(1, i), pts(2, i)) - ref)/abs(ref)
            e2 = abs(gs_coil_flux(fine, pts(1, i), pts(2, i)) - ref)/abs(ref)
            call require(e2 < 1.0e-4_dp .and. e1/e2 > 3.0_dp .and. e1/e2 < 5.0_dp, &
                'filament limit converges to rectangle flux')
        end do
    end subroutine test_rect_flux

    subroutine test_freegs_parity()
        !! FreeGS MultiCoil (0.15 m filaments), controlPsi/nfil from
        !! freegs.gradshafranov.Greens with scipy ellipk/ellipe (NumPy
        !! reference script, 5x6 filaments); per ampere-turn.
        type(gs_coil_t) :: coil
        real(dp), parameter :: pts(2, 4) = reshape([5.0_dp, 0.0_dp, 8.0_dp, -2.0_dp, &
            6.0_dp, 2.5_dp, 4.0_dp, 4.5_dp], [2, 4])
        real(dp), parameter :: ref(4) = [5.603334469564452e-07_dp, &
            5.900846227977612e-07_dp, 1.5259770719509673e-06_dp, &
            9.821052962938953e-07_dp]
        real(dp), allocatable :: fr(:), fz(:), fw(:)
        integer :: i

        coil = rect_coil(6.2_dp, 4.5_dp, 0.7_dp, 0.9_dp, 1.0_dp)
        call gs_coil_make_freegs_filaments(coil, 0.15_dp)
        call gs_coil_filament_list(coil, fr, fz, fw)
        call require(size(fr) == 30, 'FreeGS filament count 5 x 6')
        call rel_check(sum(fw), 1.0_dp, 1.0e-14_dp, 'filament turns sum')
        do i = 1, 4
            call rel_check(gs_coil_flux(coil, pts(1, i), pts(2, i)), ref(i), &
                1.0e-12_dp, 'FreeGS MultiCoil parity')
        end do
        coil%rect%turns = 554.0_dp
        call rel_check(gs_coil_flux(coil, pts(1, 1), pts(2, 1)), 554.0_dp*ref(1), &
            1.0e-12_dp, 'flux scales with turns')
    end subroutine test_freegs_parity

    function two_circuit_set() result(set)
        !! CS1U + CS1L in series (circuit 1), one single coil (circuit 2) and
        !! a reversed-winding member.
        type(gs_circuit_set_t) :: set
        allocate (set%coils(3))
        set%coils(1) = rect_coil(1.7_dp, 2.0_dp, 0.3_dp, 1.5_dp, 10.0_dp)
        set%coils(2) = rect_coil(1.7_dp, -2.0_dp, 0.3_dp, 1.5_dp, 10.0_dp)
        set%coils(3) = rect_coil(4.0_dp, 3.0_dp, 0.5_dp, 0.5_dp, 4.0_dp)
        set%coils(1)%name = 'CS1U'
        set%coils(2)%name = 'CS1L'
        set%coils(3)%name = 'PF1'
        allocate (set%circuit_names(2), set%incidence(3, 2), set%i_min(2), set%i_max(2))
        set%circuit_names = ['CS1', 'PF1']
        set%incidence = reshape([1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 2])
        set%i_min = [-4.5e4_dp, -5.2e4_dp]
        set%i_max = [4.5e4_dp, 5.2e4_dp]
    end function two_circuit_set

    subroutine test_series_incidence()
        type(gs_circuit_set_t) :: set
        real(dp) :: r(3), z(3), g(3, 2), gr(3, 2), gz(3, 2), gp(3, 2), gm(3, 2)
        real(dp) :: expect, i_coil(3), h, tot
        integer :: i, k

        set = two_circuit_set()
        r = [3.0_dp, 5.5_dp, 2.9_dp]
        z = [0.5_dp, -1.0_dp, 2.2_dp]
        call gs_circuit_flux_vectors(set, r, z, g)
        do i = 1, 3
            expect = gs_rect_coil_flux(set%coils(1)%rect, r(i), z(i), 12, 12) + &
                gs_rect_coil_flux(set%coils(2)%rect, r(i), z(i), 12, 12)
            call rel_check(g(i, 1), expect, 1.0e-13_dp, 'series flux is the coil sum')
            call rel_check(g(i, 2), gs_rect_coil_flux(set%coils(3)%rect, r(i), z(i), &
                12, 12), 1.0e-13_dp, 'single-coil circuit')
        end do
        ! Total flux from circuit currents equals the sum over coil currents.
        i_coil = gs_coil_currents(set, [2.0e4_dp, -1.0e4_dp])
        call require(all(i_coil == [2.0e4_dp, 2.0e4_dp, -1.0e4_dp]), 'coil currents')
        tot = sum(g(1, :)*[2.0e4_dp, -1.0e4_dp])
        expect = 0.0_dp
        do k = 1, 3
            expect = expect + i_coil(k)*gs_rect_coil_flux(set%coils(k)%rect, r(1), &
                z(1), 12, 12)
        end do
        call rel_check(tot, expect, 1.0e-13_dp, 'psi = sum I_k g_k')
        ! Reversed winding: incidence -1 flips the sign of that member.
        set%incidence(2, 1) = -1.0_dp
        call gs_circuit_flux_vectors(set, r, z, gm)
        expect = gs_rect_coil_flux(set%coils(1)%rect, r(1), z(1), 12, 12) - &
            gs_rect_coil_flux(set%coils(2)%rect, r(1), z(1), 12, 12)
        call rel_check(gm(1, 1), expect, 1.0e-13_dp, 'reversed series member')
        ! Point gradients against central differences.
        set = two_circuit_set()
        call gs_circuit_flux_vectors(set, r, z, g, gr, gz)
        h = 1.0e-5_dp
        call gs_circuit_flux_vectors(set, r + h, z, gp)
        call gs_circuit_flux_vectors(set, r - h, z, gm)
        call require(maxval(abs((gp - gm)/(2*h) - gr)/abs(gr)) < 1.0e-7_dp, &
            'dg/dR against finite difference')
        call gs_circuit_flux_vectors(set, r, z + h, gp)
        call gs_circuit_flux_vectors(set, r, z - h, gm)
        call require(maxval(abs((gp - gm)/(2*h) - gz)/abs(gz)) < 1.0e-7_dp, &
            'dg/dZ against finite difference')
    end subroutine test_series_incidence

    subroutine test_limits()
        type(gs_circuit_set_t) :: set
        real(dp) :: lo(2), hi(2)

        ! ITER CS module: 45 kA x 553 turns = 24.885 MAt.
        call rel_check(gs_amps_to_mat(4.5e4_dp, 553.0_dp), 24.885_dp, 1.0e-14_dp, &
            'A to MA-turns')
        call rel_check(gs_mat_to_amps(24.885_dp, 553.0_dp), 4.5e4_dp, 1.0e-14_dp, &
            'MA-turns to A')
        call rel_check(gs_mat_to_amps(gs_amps_to_mat(-1234.5_dp, 7.0_dp), 7.0_dp), &
            -1234.5_dp, 1.0e-14_dp, 'conversion round trip')
        set = two_circuit_set()
        call require(gs_circuit_turns(set, 1) == 10.0_dp .and. &
            gs_circuit_turns(set, 2) == 4.0_dp, 'circuit turns')
        call gs_circuit_limits_mat(set, lo, hi)
        call rel_check(hi(1), 0.45_dp, 1.0e-14_dp, 'circuit 1 upper limit in MAt')
        call rel_check(lo(2), -0.208_dp, 1.0e-14_dp, 'circuit 2 lower limit in MAt')
        call gs_set_circuit_limits_mat(set, [-0.1_dp, -0.2_dp], [0.3_dp, 0.4_dp])
        call rel_check(set%i_max(1), 3.0e4_dp, 1.0e-14_dp, 'set limit from MAt (max)')
        call rel_check(set%i_min(2), -5.0e4_dp, 1.0e-14_dp, 'set limit from MAt (min)')
        set%coils(1)%rect%turns = -10.0_dp
        call gs_set_circuit_limits_mat(set, [-0.1_dp, -0.2_dp], [0.3_dp, 0.4_dp])
        call require(set%i_min(1) < set%i_max(1), 'negative turns keep i_min < i_max')
    end subroutine test_limits

    function mixed_set() result(set)
        type(gs_circuit_set_t) :: set
        allocate (set%coils(4))
        set%coils(1) = rect_coil(1.7_dp, 2.0_dp, 0.3_dp, 1.5_dp, 553.0_dp)
        set%coils(1)%name = 'CS1U'
        set%coils(2) = rect_coil(6.2_dp, 4.5_dp, 0.7_dp, 0.9_dp, 248.0_dp)
        set%coils(2)%name = 'PF3'
        set%coils(2)%quad_order = 9
        call gs_coil_make_freegs_filaments(set%coils(2), 0.15_dp)
        set%coils(3)%name = 'VS'
        set%coils(3)%kind = gs_coil_filaments
        set%coils(3)%fil_r = [7.0_dp, 7.1_dp]
        set%coils(3)%fil_z = [-3.0_dp, -3.1_dp]
        set%coils(3)%fil_turns = [1.0_dp, 0.5_dp]
        set%coils(4) = rect_coil(1.7_dp, -2.0_dp, 0.3_dp, 1.5_dp, 553.0_dp)
        set%coils(4)%name = 'CS1L'
        allocate (set%circuit_names(2), set%incidence(4, 2), set%i_min(2), set%i_max(2))
        set%circuit_names = ['CS1', 'MIX']
        set%incidence = 0.0_dp
        set%incidence([1, 4], 1) = 1.0_dp
        set%incidence(2, 2) = 1.0_dp
        set%incidence(3, 2) = -0.1_dp + 1.0e-17_dp
        set%i_min = [-4.5e4_dp, -1.0_dp/3.0_dp]
        set%i_max = [4.5e4_dp, 5.2e4_dp]
    end function mixed_set

    subroutine test_table_roundtrip()
        type(gs_circuit_set_t) :: a, b
        character(len=*), parameter :: path = 'gs_coils_roundtrip.tbl'
        character(len=200) :: msg
        integer :: stat, c

        a = mixed_set()
        call gs_circuit_set_check(a, stat, msg)
        call require(stat == 0, 'mixed set is consistent')
        call gs_write_coil_table(path, a, stat)
        call require(stat == 0, 'write coil table')
        call gs_read_coil_table(path, b, stat, msg)
        call require(stat == 0, 'read coil table: '//trim(msg))
        if (stat /= 0) return
        call require(size(b%coils) == 4 .and. size(b%circuit_names) == 2, &
            'table sizes')
        do c = 1, 4
            call require(a%coils(c)%name == b%coils(c)%name .and. &
                a%coils(c)%kind == b%coils(c)%kind .and. &
                a%coils(c)%quad_order == b%coils(c)%quad_order .and. &
                a%coils(c)%fil_size == b%coils(c)%fil_size .and. &
                gs_coil_turns(a%coils(c)) == gs_coil_turns(b%coils(c)), &
                'coil metadata round trip')
            call require(gs_coil_flux(a%coils(c), 3.0_dp, 0.7_dp) == &
                gs_coil_flux(b%coils(c), 3.0_dp, 0.7_dp), 'coil flux round trip')
        end do
        call require(all(a%incidence == b%incidence), 'incidence round trip')
        call require(all(a%i_min == b%i_min) .and. all(a%i_max == b%i_max), &
            'limits round trip')
        call require(all(a%circuit_names == b%circuit_names), 'circuit names')
        call require(all(a%coils(3)%fil_r == b%coils(3)%fil_r) .and. &
            all(a%coils(3)%fil_z == b%coils(3)%fil_z) .and. &
            all(a%coils(3)%fil_turns == b%coils(3)%fil_turns), 'filament round trip')
        call require(a%coils(1)%rect%r == b%coils(1)%rect%r .and. &
            a%coils(1)%rect%dz == b%coils(1)%rect%dz, 'rect geometry')
    end subroutine test_table_roundtrip

    subroutine write_text(path, lines)
        character(len=*), intent(in) :: path, lines(:)
        integer :: u, i
        open (newunit=u, file=path, status='replace', action='write')
        do i = 1, size(lines)
            write (u, '(a)') trim(lines(i))
        end do
        close (u)
    end subroutine write_text

    subroutine test_table_errors()
        type(gs_circuit_set_t) :: set
        character(len=*), parameter :: path = 'gs_coils_bad.tbl'
        character(len=200) :: msg
        integer :: stat

        call write_text(path, [character(len=80) :: 'coil A rect 1 0 1 1 1'])
        call gs_read_coil_table(path, set, stat, msg)
        call require(stat /= 0, 'missing header rejected')
        call write_text(path, [character(len=80) :: 'kin6d-coil-table 1', &
            'coil A rect 1 0 1 1 1', 'circuit C -1 1 B'])
        call gs_read_coil_table(path, set, stat, msg)
        call require(stat /= 0, 'unknown circuit member rejected')
        call write_text(path, [character(len=80) :: 'kin6d-coil-table 1', &
            'coil A rect 1 0 x 1 1'])
        call gs_read_coil_table(path, set, stat, msg)
        call require(stat /= 0, 'bad number rejected')
        call write_text(path, [character(len=80) :: 'kin6d-coil-table 1', &
            '# comment line', 'coil A rect 1 0 1 1 1  # trailing', &
            'circuit C -2 3 A*-1'])
        call gs_read_coil_table(path, set, stat, msg)
        call require(stat == 0, 'comments and factors accepted: '//trim(msg))
        if (stat == 0) call require(set%incidence(1, 1) == -1.0_dp .and. &
            set%i_min(1) == -2.0_dp .and. set%i_max(1) == 3.0_dp, 'factor and limits')
    end subroutine test_table_errors

    subroutine test_helmholtz()
        !! On-axis Bz of a Helmholtz pair (radius a, separation a, series
        !! circuit): Bz = (4/5)^(3/2) mu0 I / a per ampere-turn.  Bz = (1/R)
        !! dpsi/dR is evaluated at R = 0.01 a, where the first correction is
        !! O((R/a)^4).  Br = -(1/R) dpsi/dZ vanishes at Z = 0 by symmetry.
        type(gs_circuit_set_t) :: set
        real(dp) :: r(1), z(1), g(1, 1), gr(1, 1), gz(1, 1), bz, br, expect
        real(dp), parameter :: a = 1.3_dp
        integer :: variant

        expect = 0.8_dp**1.5_dp*gs_mu0/a
        r = 0.01_dp*a
        z = 0.0_dp
        do variant = 1, 2
            allocate (set%coils(2), set%circuit_names(1), set%incidence(2, 1), &
                set%i_min(1), set%i_max(1))
            set%circuit_names = ['HH']
            set%incidence = 1.0_dp
            set%i_min = -1.0_dp
            set%i_max = 1.0_dp
            if (variant == 1) then
                set%coils(1)%kind = gs_coil_filaments
                set%coils(2)%kind = gs_coil_filaments
                set%coils(1)%fil_r = [a]
                set%coils(1)%fil_z = [0.5_dp*a]
                set%coils(1)%fil_turns = [1.0_dp]
                set%coils(2)%fil_r = [a]
                set%coils(2)%fil_z = [-0.5_dp*a]
                set%coils(2)%fil_turns = [1.0_dp]
            else
                set%coils(1) = rect_coil(a, 0.5_dp*a, 0.02_dp, 0.02_dp, 1.0_dp)
                set%coils(2) = rect_coil(a, -0.5_dp*a, 0.02_dp, 0.02_dp, 1.0_dp)
            end if
            call gs_circuit_flux_vectors(set, r, z, g, gr, gz)
            bz = gr(1, 1)/r(1)
            br = -gz(1, 1)/r(1)
            if (variant == 1) then
                call rel_check(bz, expect, 1.0e-6_dp, 'Helmholtz Bz, filament loops')
            else
                call rel_check(bz, expect, 2.0e-4_dp, 'Helmholtz Bz, rectangle coils')
            end if
            call require(abs(br) < 1.0e-8_dp*expect, 'Helmholtz Br vanishes on midplane')
            deallocate (set%coils, set%circuit_names, set%incidence, set%i_min, set%i_max)
        end do
    end subroutine test_helmholtz

end program test_gs_coils
