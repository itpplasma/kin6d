!> Boozer `.bc` geometry: closed forms on synthetic six- and ten-column
!! (up-down asymmetric) files, refusal outside
!! the tabulated surfaces, JVP/VJP pairing of the metric leaf, and agreement
!! with NEO-2 on a public tokamak file when its path and NEO-2's values are
!! given (args: file s iota diota_ds av_nabla_stor[1/cm] avb2[T^2]).
program test_boozer_file
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_boozer_file
    use kin6d_boozer_file_leaves, only: boozer_metric, boozer_metric_jvp, &
        boozer_metric_vjp
    implicit none
    real(dp), parameter :: r0 = 1.7_dp, a = 0.5_dp, b0 = 2.3_dp, eps = 0.1_dp
    integer :: failures = 0

    call synthetic_file()
    call asymmetric_file()
    call signed_poloidal_modes()
    call high_mode_roundoff()
    call metric_pairing()
    if (command_argument_count() >= 6) call neo2_comparison()
    if (failures > 0) error stop 'FAIL'
    print '(a)', 'All tests passed!'
contains

    subroutine synthetic_file()
        character(*), parameter :: name = 'synthetic.bc'
        type(boozer_file_t) :: bf
        type(boozer_surface_t) :: surf
        type(boozer_point_t) :: p
        real(dp) :: s, theta, grad_s, b2, delta
        integer :: unit, k, status, i
        real(dp), parameter :: ss(9) = [(0.1_dp*k, k=1, 9)]

        ! Circular surfaces r = a sqrt(s): R = r0 + r cos(theta), Z = r sin(theta),
        ! so s = r^2/a^2, |grad s| = 2 sqrt(s)/a and sqrt(g) = -a^2 R/2.
        open (newunit=unit, file=name, status='replace', action='write')
        write (unit, '(a)') 'CC synthetic test file'
        write (unit, '(a)') ' m0b n0b nsurf nper flux/[Tm^2] a/[m] R/[m]'
        write (unit, '(a)') ' 1 0 9 1 -1.0 0.5 1.7'
        do k = 1, 9
            write (unit, '(a)') ' s iota curr_pol/nper curr_tor pprime sqrt g(0,0)'
            write (unit, '(a)') ' [A] [A] dp/ds,[Pa] (dV/ds)/nper'
            write (unit, '(6es24.16)') ss(k), 0.9_dp - 0.3_dp*ss(k), 0.0_dp, &
                0.0_dp, 0.0_dp, 0.0_dp
            write (unit, '(a)') ' m n r/[m] z/[m] (phib-phi)*nper/twopi bmn/[T]'
            write (unit, '(2i5,4es24.16)') 0, 0, r0, 0.0_dp, 0.0_dp, b0
            write (unit, '(2i5,4es24.16)') 1, 0, a*sqrt(ss(k)), a*sqrt(ss(k)), &
                0.0_dp, -eps*b0*sqrt(ss(k))
        end do
        close (unit)
        call load_boozer_file(name, bf, status)
        call check('load status', status == BOOZER_OK)

        do i = 1, 2
            s = merge(0.5_dp, 0.43_dp, i == 1)
            call boozer_surface(bf, s, surf, status)
            call check('surface status', status == BOOZER_OK)
            call close_to('iota', surf%iota, 0.9_dp - 0.3_dp*s, 1e-14_dp)
            call close_to('diota_ds', surf%diota_ds, -0.3_dp, 1e-13_dp)
            theta = 0.7_dp
            p = boozer_point(surf, theta, 0.3_dp)
            call close_to('R', p%r, r0 + a*sqrt(s)*cos(theta), 1e-14_dp)
            call close_to('|grad s|', p%grad_s, 2*sqrt(s)/a, 1e-13_dp)
            call close_to('jacobian', p%jacobian, -a*a*p%r/2, 1e-13_dp)
            call boozer_average(surf, 16, 1, grad_s, b2)
            delta = a*sqrt(s)/r0
            call close_to('<|grad s|>', grad_s, 2*sqrt(s)/a, 1e-13_dp)
            call close_to('<B^2>', b2, b0**2*(1 - eps*sqrt(s)*delta + eps**2*s/2), 1e-13_dp)
        end do

        call boozer_surface(bf, 0.05_dp, surf, status)
        call check('refuse s below the file', status == BOOZER_OUTSIDE)
        call boozer_surface(bf, 0.95_dp, surf, status)
        call check('refuse s above the file', status == BOOZER_OUTSIDE)
        open (newunit=unit, file=name)
        close (unit, status='delete')
    end subroutine synthetic_file

    !> The circular surfaces of `synthetic_file` shifted to Z = z0 by an m = 0
    !> zmnc term, with sine terms in |B| and in the toroidal-angle shift.
    subroutine asymmetric_file()
        character(*), parameter :: name = 'synthetic_asym.bc'
        real(dp), parameter :: z0 = 0.07_dp, beta = 0.03_dp, nu = 0.002_dp
        type(boozer_file_t) :: bf
        type(boozer_surface_t) :: surf
        type(boozer_point_t) :: p
        real(dp) :: s, theta, r
        integer :: unit, k, status
        real(dp), parameter :: ss(9) = [(0.1_dp*k, k=1, 9)]
        real(dp), parameter :: twopi = 8*atan(1.0_dp)

        open (newunit=unit, file=name, status='replace', action='write')
        write (unit, '(a)') 'CC synthetic up-down asymmetric test file'
        write (unit, '(a)') ' m0b n0b nsurf nper flux/[Tm^2] a/[m] R/[m]'
        write (unit, '(a)') ' 1 0 9 1 -1.0 0.5 1.7'
        do k = 1, 9
            write (unit, '(a)') ' s iota curr_pol/nper curr_tor pprime sqrt g(0,0)'
            write (unit, '(a)') ' [A] [A] dp/ds,[Pa] (dV/ds)/nper'
            write (unit, '(6es24.16)') ss(k), 0.9_dp - 0.3_dp*ss(k), 0.0_dp, &
                0.0_dp, 0.0_dp, 0.0_dp
            write (unit, '(a)') ' m n rmnc rmns zmnc zmns vmnc vmns bmnc bmns'
            write (unit, '(2i5,8es24.16)') 0, 0, r0, 0.0_dp, z0, 0.0_dp, nu, &
                0.0_dp, b0, 0.0_dp
            write (unit, '(2i5,8es24.16)') 1, 0, a*sqrt(ss(k)), 0.0_dp, 0.0_dp, &
                a*sqrt(ss(k)), 0.0_dp, 0.0_dp, -eps*b0*sqrt(ss(k)), &
                beta*b0*sqrt(ss(k))
        end do
        close (unit)
        call load_boozer_file(name, bf, status)
        call check('asymmetric load status', status == BOOZER_OK)
        s = 0.43_dp
        call boozer_surface(bf, s, surf, status)
        call check('asymmetric surface status', status == BOOZER_OK)
        theta = 0.7_dp
        r = a*sqrt(s)
        p = boozer_point(surf, theta, 0.3_dp)
        call close_to('asymmetric R', p%r, r0 + r*cos(theta), 1e-14_dp)
        call close_to('asymmetric Z', p%z, z0 + r*sin(theta), 1e-14_dp)
        call close_to('asymmetric |B|', p%bmod, &
            b0*(1 - eps*sqrt(s)*cos(theta) + beta*sqrt(s)*sin(theta)), 1e-14_dp)
        call close_to('asymmetric phi', p%phi, 0.3_dp - twopi*nu, 1e-14_dp)
        call close_to('asymmetric |grad s|', p%grad_s, 2*sqrt(s)/a, 1e-13_dp)
        call close_to('asymmetric jacobian', p%jacobian, -a*a*p%r/2, 1e-13_dp)
        open (newunit=unit, file=name)
        close (unit, status='delete')
    end subroutine asymmetric_file

    subroutine signed_poloidal_modes()
        character(*), parameter :: name = 'synthetic_signed_mode.bc'
        real(dp), parameter :: c = 0.02_dp, sn = -0.03_dp, s = 0.43_dp
        integer, parameter :: modes(3) = [4, 16, 17]
        type(boozer_file_t) :: bf
        type(boozer_surface_t) :: surf
        type(boozer_point_t) :: p
        real(dp) :: ss, phase, theta, phi, expected, power
        integer :: unit, k, status, sign_m, j, m

        theta = 0.7_dp; phi = 0.3_dp
        do j = 1, size(modes)
            m = modes(j); power = 0.5_dp*m
            phase = m*theta - 2*phi
            do sign_m = -1, 1, 2
                open (newunit=unit, file=name, status='replace', action='write')
                write (unit, '(a)') 'CC signed-mode representation of the same regular field'
                write (unit, '(a)') ' m0b n0b nsurf nper flux a R'
                write (unit, '(4i5,3es24.16)') m, 1, 9, 2, -1.0_dp, 0.5_dp, 1.7_dp
                do k = 1, 9
                    ss = 0.1_dp*k
                    write (unit, '(a)') ' s iota Jpol/nper Itor pprime sqrtg'
                    write (unit, '(6es24.16)') ss, 0.9_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp
                    write (unit, '(a)') ' m n rmnc rmns zmnc zmns vmnc vmns bmnc bmns'
                    write (unit, '(2i5,8es24.16)') 0, 0, r0, 0.0_dp, 0.0_dp, 0.0_dp, &
                        0.0_dp, 0.0_dp, b0, 0.0_dp
                    write (unit, '(2i5,8es24.16)') 1, 0, a*sqrt(ss), 0.0_dp, 0.0_dp, &
                        a*sqrt(ss), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp
                    write (unit, '(2i5,8es24.16)') m*sign_m, sign_m, 0.0_dp, 0.0_dp, &
                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, c*ss**power, sign_m*sn*ss**power
                end do
                close (unit)
                call load_boozer_file(name, bf, status)
                call check('signed-mode load status', status == BOOZER_OK)
                call boozer_surface(bf, s, surf, status)
                call check('signed-mode surface status', status == BOOZER_OK)
                call close_to('regular cosine coefficient', surf%bmn(3), c*s**power, 1e-14_dp)
                call close_to('regular sine coefficient', surf%bma(3), sign_m*sn*s**power, 1e-14_dp)
                call close_to('regular cosine derivative', surf%bmn_s(3), power*c*s**(power - 1), 1e-14_dp)
                call close_to('regular sine derivative', surf%bma_s(3), sign_m*power*sn*s**(power - 1), 1e-14_dp)
                p = boozer_point(surf, theta, phi)
                expected = b0 + s**power*(c*cos(phase) + sn*sin(phase))
                call close_to('signed-mode field [T]', p%bmod, expected, 1e-14_dp)
                open (newunit=unit, file=name)
                close (unit, status='delete')
            end do
        end do
    end subroutine signed_poloidal_modes

    subroutine high_mode_roundoff()
        character(*), parameter :: name = 'synthetic_high_mode.bc'
        real(dp), parameter :: ss(3) = [0.0031_dp, 0.056141414141414_dp, 0.6001_dp]
        type(boozer_file_t) :: bf
        type(boozer_surface_t) :: surf
        integer :: unit, k, status
        real(dp) :: s

        open (newunit=unit, file=name, status='replace', action='write')
        write (unit, '(a)') 'CC high-mode coefficient at the numerical-noise floor'
        write (unit, '(a)') ' m0b n0b nsurf nper flux a R'
        write (unit, '(a)') ' 64 1 1000 2 -1.0 0.5 1.7'
        do k = 1, 1000
            s = (k - 0.5_dp)/1000
            write (unit, '(a)') ' s iota Jpol/nper Itor pprime sqrtg'
            write (unit, '(6es24.16)') s, 0.9_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp
            write (unit, '(a)') ' m n rmnc zmns vmns bmnc'
            write (unit, '(2i5,4es24.16)') 0, 0, r0, 0.0_dp, 0.0_dp, b0
            write (unit, '(2i5,4es24.16)') 64, 1, 0.0_dp, 0.0_dp, 0.0_dp, 1e-20_dp
        end do
        close (unit)
        call load_boozer_file(name, bf, status)
        call check('high-mode load status', status == BOOZER_OK)
        do k = 1, size(ss)
            call boozer_surface(bf, ss(k), surf, status)
            call check('high-mode surface status', status == BOOZER_OK)
            call check('high-mode noise remains small [T]', abs(surf%bmn(2)) < 1e-12_dp)
            call check('high-mode derivative remains small [T]', abs(surf%bmn_s(2)) < 1e-9_dp)
        end do
        open (newunit=unit, file=name)
        close (unit, status='delete')
    end subroutine high_mode_roundoff

    subroutine metric_pairing()
        real(dp) :: x(10), v(10), g(10), w(2), jv(2), lhs, rhs
        integer :: k
        call random_seed(put=[(17 + k, k=1, 64)])
        call random_number(x); call random_number(v); call random_number(w)
        x(1) = x(1) + 1
        call boozer_metric_jvp(x(1), x(2), x(3), x(4), x(5), x(6), x(7), &
            x(8), x(9), x(10), v(1), v(2), v(3), v(4), &
            v(5), v(6), v(7), v(8), v(9), v(10), jv(1), jv(2))
        call boozer_metric_vjp(x(1), x(2), x(3), x(4), x(5), x(6), x(7), &
            x(8), x(9), x(10), w(1), w(2), g(1), g(2), &
            g(3), g(4), g(5), g(6), g(7), g(8), g(9), g(10))
        lhs = dot_product(w, jv)
        rhs = dot_product(g, v)
        call close_to('metric JVP/VJP pairing', lhs, rhs, 1e-13_dp)
    end subroutine metric_pairing

    subroutine neo2_comparison()
        character(512) :: file
        real(dp) :: ref(5), grad_s, b2
        type(boozer_file_t) :: bf
        type(boozer_surface_t) :: surf
        integer :: k, status
        character(64) :: arg

        call get_command_argument(1, file)
        do k = 1, 5
            call get_command_argument(k + 1, arg)
            read (arg, *) ref(k)
        end do
        call load_boozer_file(trim(file), bf, status)
        call boozer_surface(bf, ref(1), surf, status)
        call check('NEO-2 surface status', status == BOOZER_OK)
        call boozer_average(surf, 128, 1, grad_s, b2)
        print '(a,4es10.2)', 'relative deviation from NEO-2 (iota, diota/ds, <|grad s|>, <B^2>):', &
            surf%iota/ref(2) - 1, surf%diota_ds/ref(3) - 1, &
            grad_s*boozer_jacobian_ratio(surf, 128)/(100*ref(4)) - 1, b2/ref(5) - 1
        call close_to('NEO-2 iota', surf%iota, ref(2), 1e-6_dp)
        call close_to('NEO-2 diota/ds', surf%diota_ds, ref(3), 1e-4_dp)
        ! NEO-2 divides by the flux-function Jacobian; this file's geometric
        ! Jacobian is smaller by about 1.5e-4, so compare through the ratio.
        call close_to('NEO-2 <|grad s|>', grad_s*boozer_jacobian_ratio(surf, 128)/100, &
            ref(4), 1e-6_dp)
        call close_to('NEO-2 <B^2>', b2, ref(5), 1e-4_dp)
    end subroutine neo2_comparison

    subroutine close_to(label, got, want, rtol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, want, rtol
        if (abs(got - want) > rtol*abs(want)) then
            print '(a,2es24.16)', 'FAIL '//label//': ', got, want
            failures = failures + 1
        end if
    end subroutine close_to

    subroutine check(label, ok)
        character(*), intent(in) :: label
        logical, intent(in) :: ok
        if (.not. ok) then
            print '(a)', 'FAIL '//label
            failures = failures + 1
        end if
    end subroutine check
end program test_boozer_file
