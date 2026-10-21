module geqdsk_analytic_oracle
    ! Independent polynomial Solov'ev oracle, also the declared continuation.
    ! Same published family as test_equilibrium_core, not a production leaf.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_state, only: equilibrium_field_t, EQ_BR, EQ_BZ, &
                                       EQ_BPHI, EQ_P, EQ_PSI, EQ_PSI_R, EQ_PSI_Z
    implicit none
    real(dp), parameter :: r0 = 2, minor = 0.3_dp, f0 = 4, c = f0/(12*r0**3)
    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 1.2566370614359173e-6_dp
    type, extends(equilibrium_field_t) :: oracle_t
        real(dp) :: pol = 1, tor = 1
    contains
        procedure :: sample => sample_oracle
    end type
contains
    subroutine sample_oracle(field, points, values, statuses)
        class(oracle_t), intent(in) :: field
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        real(dp) :: r, z, flux, pr, pz
        integer :: i
        values = 0; statuses = 0
        do i = 1, size(points, 2)
            r = points(1, i); z = points(2, i)
            flux = c*((r*r - r0*r0)**2 + 4*r*r*z*z - 4*r0*r0*minor*minor)
            pr = c*(4*r*(r*r - r0*r0) + 8*r*z*z); pz = 8*c*r*r*z
            values(EQ_PSI, i) = field%pol*flux
            values(EQ_PSI_R, i) = field%pol*pr
            values(EQ_PSI_Z, i) = field%pol*pz
            values(EQ_BR, i) = -field%pol*pz/r
            values(EQ_BZ, i) = field%pol*pr/r
            values(EQ_BPHI, i) = field%tor*f0/r
            values(EQ_P, i) = -16*c*flux/mu0
        end do
    end subroutine

    real(dp) function q_oracle(s) result(q)
        real(dp), intent(in) :: s
        integer :: j
        real(dp) :: radius
        q = 0
        ! Analytic coarea Jacobian: dl/(R |grad psi|) = dtheta/(8 c R^3).
        do j = 1, 4096
            radius = sqrt(r0*r0 + 2*r0*minor*sqrt(s)*cos(2*pi*j/4096))
            q = q + f0/(8*c*radius**3)/4096
        end do
    end function
end module

program test_geqdsk_adapter
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use geqdsk_analytic_oracle, only: oracle_t, r0, minor, f0, c, pi, mu0, q_oracle
    use kin6d_geqdsk_adapter, only: geqdsk_data_t, geqdsk_import_t, export_geqdsk, &
                                import_geqdsk, geqdsk_grid, resolve_geqdsk, native_cocos
    use kin6d_equilibrium_state, only: equilibrium_state_t, EQ_VALUES, EQ_PSI, EQ_BZ
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs
    use kin6d_gs_equilibrium, only: gs_equilibrium
    use kin6d_gs_ampere_loop, only: gs_boundary_current
    use kin6d_equilibrium_readback, only: equilibrium_readback_t, read_equilibrium
    use kin6d_vacuum_surfaces, only: psi_grid_t, psi_grid_init, eqdsk_read, psi_eval, &
                                     profiles_of_psin, surface_t, surface_build
    implicit none
    character(1024) :: mode, argument
    integer :: polarity, status
    call get_command_argument(1, mode)
    select case (trim(mode))
    case ('roundtrip')
        call get_command_argument(2, argument)
        read (argument, *) polarity
        call roundtrip(polarity, .false.)
    case ('resolve')
        call roundtrip(1, .true.)
    case ('chease')
        call get_command_argument(2, argument)
        call chease(trim(argument))
    case default
        error stop 'roundtrip polarity | resolve | chease path'
    end select
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print *, 'FAIL: ', message
            error stop 1
        end if
    end subroutine

    subroutine roundtrip(polarity, resolving)
        integer, intent(in) :: polarity
        logical, intent(in) :: resolving
        type(gs_state_t) :: mesh, solved
        type(equilibrium_state_t) :: eq, recovered
        type(oracle_t) :: oracle
        type(geqdsk_data_t) :: data
        type(geqdsk_import_t) :: imported, rejected
        type(equilibrium_readback_t) :: rb
        real(dp) :: boundary(2, 96), rectangle(4), ip, t
        real(dp) :: values(EQ_VALUES, 1), point(2, 1), maxpsi, maxq, exact
        real(dp), allocatable :: rho(:)
        integer :: i, j, unit, st(1)
        character(80) :: path
        oracle%pol = 1; oracle%tor = 1
        if (polarity == 2 .or. polarity == 4) oracle%pol = -1
        if (polarity == 3 .or. polarity == 4) oracle%tor = -1
        write (path, '(a,i0,a,l1,a)') 'solovev-', polarity, '-', resolving, '.eqdsk'
        do i = 1, 96
            t = 2*pi*(i - 1)/96
            boundary(1, i) = sqrt(r0*r0 + 2*r0*minor*cos(t))
            boundary(2, i) = r0*minor*sin(t)/boundary(1, i)
        end do
        call initialize_gs(mesh, boundary, 0.0015_dp, status, element_order=2)
        call require(status == 0, 'Solovev mesh')
        call solve_gs(mesh, mu0, -oracle%pol*16*c/mu0, 0.0_dp, status)
        call require(status == 0, 'Solovev native solve')
        call gs_equilibrium(mesh, eq, status, oracle%tor*f0, 0.0_dp, &
                            pprime=-oracle%pol*16*c/mu0)
        call require(status == 0, 'D7 adapter')
        call gs_boundary_current(mesh, mu0, 8, ip, status)
        call require(status == 0, 'physical Ip')
        ! Runtime qualification of the native convention: physical Bz outside
        ! the axis, signed q and Ampere current must have the COCOS 3 signs.
        point(:, 1) = [r0 + minor/2, 0.0_dp]
        call eq%field%sample(point, values, st)
        call require(st(1) == 0, 'sign sample')
        call require(values(EQ_BZ, 1)*oracle%pol > 0, 'native Bz sign')
        call require(ip*oracle%pol < 0, 'native Ampere Ip sign')
        rectangle = [r0 - 1.3_dp*minor, r0 + 1.3_dp*minor, &
                     -1.3_dp*minor, 1.3_dp*minor]
        if (polarity == 1) then
            call export_geqdsk(path, eq, rectangle, 33, 33, boundary, boundary, ip, &
                               r0, native_cocos, data, status)
            call require(status == 1, 'missing exterior continuation rejected')
        end if
        call export_geqdsk(path, eq, rectangle, 33, 33, boundary, boundary, ip, &
                           r0, native_cocos, data, status, continuation=oracle)
        call require(status == 0, 'export with explicit analytic continuation')
        call import_geqdsk(path, native_cocos, imported, status)
        call require(status == 0, 'explicit COCOS import')
        call require(imported%classified_cocos == 3, 'COCOS 3 classification')
        call import_geqdsk(path, 1, rejected, status)
        call require(status == 2, 'wrong COCOS declaration rejected')
        call require(maxval(abs(imported%data%psi - data%psi)) < 1.0e-9_dp, &
                     'SI flux roundtrip')
        call require(maxval(abs(imported%data%f - data%f)) < 1.0e-8_dp, 'SI F')
        call require(maxval(abs(imported%data%p - data%p)) < 1.0e-4_dp, 'SI p')
        call require(maxval(abs(imported%data%q - data%q)) < 1.0e-8_dp, 'signed q')
        call require(abs(imported%data%ip/ip - 1) < 1.0e-8_dp, 'SI Ip')
        call require(imported%data%bcentr*oracle%tor > 0, 'B0 polarity')
        maxpsi = 0; maxq = 0
        do j = 1, 33
            do i = 1, 33
                point(:, 1) = [data%r(i), data%z(j)]
                call oracle%sample(point, values, st)
                maxpsi = max(maxpsi, abs(data%psi(i, j) - values(EQ_PSI, 1)))
            end do
        end do
        do i = 1, 33
            exact = oracle%pol*oracle%tor*q_oracle(real(i - 1, dp)/32)
            maxq = max(maxq, abs(imported%data%q(i) - exact)/abs(exact))
        end do
        print '(a,i0,a,3es14.5)', 'polarity ', polarity, &
            ': max psi error [Wb/rad], relative q error, Ip [A] ', maxpsi, maxq, ip
        call require(maxpsi < 3.0e-4_dp, 'analytic flux discretisation budget')
        call require(maxq < 0.03_dp, 'analytic q discretisation budget')
        call require(maxval(abs(imported%data%p + &
                     oracle%pol*16*c*imported%data%flux/mu0)) < 1.0e-4_dp, 'analytic p')
        if (resolving) then
            call resolve_geqdsk(imported%data, 0.0007_dp, solved, recovered, status)
            call require(status == 0, 'native cold re-solve')
            rho = [(sqrt(real(i - 1, dp)/32), i=1, 33)]
            call read_equilibrium(recovered, 33, 96, rb, status, rho_grid=rho)
            call require(status == 0, 'resolved D7')
          print '(a,3es14.5)', 'native resolve: axis distance [m], flux, q relative ', &
                norm2(rb%axis - data%axis), abs(rb%psi_axis - data%psi_axis), &
                maxval(abs((rb%q - data%q)/data%q))
            call require(norm2(rb%axis - data%axis) < 0.002_dp, 'resolve axis')
            call require(abs(rb%psi_axis - data%psi_axis) < 3.0e-4_dp, 'resolve psi')
            call require(maxval(abs((rb%q - data%q)/data%q)) < 0.03_dp, 'resolve q')
        end if
        open (newunit=unit, file=path, status='old'); close (unit, status='delete')
    end subroutine

    subroutine chease(path)
        character(*), intent(in) :: path
        type(geqdsk_import_t) :: imported
        type(psi_grid_t) :: grid, legacy, half
        type(surface_t) :: surface, half_surface
        type(gs_state_t) :: solved
        type(equilibrium_state_t) :: eq
        type(equilibrium_readback_t) :: rb
        real(dp) :: p1, pr1, pz1, p2, pr2, pz2, f1, f2, pp1, pp2, d1, d2
        real(dp) :: point(2), s, grid_axis_error, grid_q_error, q_error
        real(dp) :: psi_error, axis_error, area, errors(3, 3)
        real(dp), allocatable :: rho(:)
        integer :: i, j, level
        call import_geqdsk(path, 1, imported, status)
        call require(status == 0, 'CHEASE source COCOS 1')
        call require(imported%classified_cocos == 1, 'CHEASE classification')
        call require(maxval(abs(imported%data%ff + imported%original%ff)) &
                     < 1.0e-12_dp, 'FFprime SI/COCOS boundary')
        print '(a,es14.5)', 'explicit native edge pressure offset [Pa]: ', &
            -imported%data%p(65)
        call geqdsk_grid(imported%data, grid, status)
        call require(status == 0, 'canonical SI grid consumer')
        call eqdsk_read(path, legacy, status)
        call require(status == 0, 'legacy fixture reader')
        do j = 2, 64, 7
            do i = 2, 64, 7
                point = [imported%original%r(i), imported%original%z(j)]
                call psi_eval(legacy, point(1), point(2), p1, pr1, pz1)
                call psi_eval(grid, point(1), point(2), p2, pr2, pz2)
                ! Source COCOS 1 -> native 3 changes the flux sign; physical
                ! fields and the pressure/F profiles remain equivalent.
                call require(abs(p1 + p2) < 1.0e-12_dp, 'legacy psi equivalence')
                call require(maxval(abs([pr1 + pr2, pz1 + pz2])) < 1.0e-10_dp, &
                             'legacy physical poloidal field equivalence')
            end do
        end do
        do i = 0, 16
            s = real(i, dp)/16
            call profiles_of_psin(legacy, s, f1, d1, pp1, pr1)
            call profiles_of_psin(grid, s, f2, d2, pp2, pr2)
            call require(maxval(abs([f1 - f2, pp1 - pp2])) < 1.0e-8_dp, &
                         'legacy F and pressure equivalence')
        end do
        call psi_grid_init(half, imported%data%psi(::2, ::2), &
                         imported%data%r(1), imported%data%r(65) - imported%data%r(1), &
                           0.0_dp, imported%data%z(65) - imported%data%z(1), &
                           imported%data%axis(1), imported%data%axis(2), &
                      imported%data%psi_edge, imported%data%ip, imported%data%f, status)
        call require(status == 0, 'coarsened grid error estimate')
        grid_axis_error = norm2([half%r_axis - grid%r_axis, half%z_axis - grid%z_axis])
        grid_q_error = 0
        do i = 1, 15
            s = real(i, dp)/16
            call surface_build(grid, s, 128, surface, status)
            call require(status == 0, 'full grid surface')
            call surface_build(half, s, 128, half_surface, status)
            call require(status == 0, 'half grid surface')
            grid_q_error = max(grid_q_error, abs(surface%q - half_surface%q))
        end do
        print '(a,i0,a,l1)', 'CHEASE source COCOS ', imported%classified_cocos, &
            '; classifier repaired profiles: ', imported%profiles_repaired
        print '(a,3es14.5)', 'grid estimates: axis [m], flux [Wb/rad], max q ', &
            grid_axis_error, abs(half%psi_axis - grid%psi_axis), grid_q_error
        rho = [(sqrt(real(i - 1, dp)/16), i=1, 17)]
        do level = 1, 3
            area = 0.0012_dp/4**(level - 1)
            call resolve_geqdsk(imported%data, area, solved, eq, status, &
                                edge_pressure=0.0_dp)
            call require(status == 0, 'CHEASE native cold solve')
            call read_equilibrium(eq, 17, 96, rb, status, rho_grid=rho)
            call require(status == 0, 'CHEASE native D7')
            psi_error = abs(rb%psi_axis + imported%data%psi_edge &
                            - imported%data%psi_axis)
            axis_error = norm2(rb%axis - imported%data%axis)
            q_error = maxval(abs(rb%q - imported%data%q(::4)))
            errors(:, level) = [axis_error, psi_error, q_error]
            print '(a,i0,a,4es14.5)', 'CHEASE level ', level, &
                ': axis [m], psi_axis absolute [Wb/rad], q_axis, max q error ', &
                axis_error, rb%psi_axis + imported%data%psi_edge, rb%q(1), q_error
        end do
        print '(a,5f12.7)', 'CHEASE reference q at s=0,.25,.5,.75,1: ', &
            imported%data%q(::16)
        print '(a,5f12.7)', 'CHEASE native q at s=0,.25,.5,.75,1:    ', rb%q(::4)
        print '(a,3es14.5)', 'CHEASE fine absolute errors (axis, psi_axis, q): ', &
            errors(:, 3)
        print '(a,3es14.5)', 'CHEASE native refinement changes in these errors: ', &
            abs(errors(:, 3) - errors(:, 2))
        if (errors(1, 3) > grid_axis_error .or. &
            errors(2, 3) > abs(half%psi_axis - grid%psi_axis) .or. &
            errors(3, 3) > grid_q_error) then
            print '(a)', 'FINDING: native mismatch exceeds a grid estimate; '// &
                'candidate comparison, no agreement or enclosure claim.'
        end if
    end subroutine
end program
