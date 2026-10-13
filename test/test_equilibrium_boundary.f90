program test_equilibrium_boundary
    ! D6 boundary layer (docs/equilibrium-core.md, M4).
    !   exterior  per-n Steklov-Poincare map on a fixed axisymmetric Gamma_c
    !             against an exact decaying field: a tilted current loop
    !             inside Gamma_c has, outside it, B = grad Phi with Phi_n
    !             solving L_n Phi_n = 0 and decaying, so the core traces
    !             u = Phi_n and h = int R B.n N_j ds satisfy S_n u + h = 0
    !             up to discretisation, for n = 0..3 (n = 0 after fitting the
    !             free constant of its line-integral trace), with convergence
    !             under edge refinement;
    !   coils     the two coil sets behind the one interface agree: the n = 0
    !             field of an axisymmetric filament (exact flux) equals that
    !             of the same circle as a 3D polygon (Biot-Savart + DFT) to
    !             the polygon error, and the polygon has no n /= 0 part.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_boundary, only: gamma_c_t, gamma_c_build, &
        gamma_c_exterior, gamma_c_coil_traces, axisymmetric_coils_t, &
        polyline_coils_t, GAMMA_C_SCALAR, GAMMA_C_OK
    use kin6d_coils_3d, only: coil_set_from_polylines
    use kin6d_gs_coils, only: gs_coil_filaments
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(16) :: argument
    call get_command_argument(1, argument)
    select case (trim(argument))
    case ('exterior')
        call exterior_oracle()
    case ('coils')
        call coil_oracle()
    case default
        error stop 'usage: test_equilibrium_boundary exterior|coils'
    end select
contains
    subroutine ellipse(n_edges, gc)
        ! P2 trace of the ellipse R = 1 + 0.45 cos t, Z = 0.6 sin t.
        integer, intent(in) :: n_edges
        type(gamma_c_t), intent(out) :: gc
        real(dp) :: points(2, 2*n_edges), t
        integer :: j, status
        do j = 1, 2*n_edges
            t = 2*pi*real(j - 1, dp)/(2*n_edges)
            points(:, j) = [1.0_dp + 0.45_dp*cos(t), 0.6_dp*sin(t)]
        end do
        call gamma_c_build(points, 2, 1, gc, status)
        if (status /= GAMMA_C_OK) error stop 'ellipse curve'
    end subroutine

    subroutine tilted_loop(coils)
        ! Circle of radius 0.12 about (x, y, z) = (1.05, 0, 0.08) in a plane
        ! tilted against every axis, 512 segments.
        type(polyline_coils_t), intent(out) :: coils
        integer, parameter :: np = 512
        real(dp) :: pts(3, np), normal(3), e1(3), e2(3), s
        integer :: j, status
        normal = [sin(0.5_dp), 0.3_dp, cos(0.5_dp)]
        normal = normal/norm2(normal)
        e1 = [normal(3), 0.0_dp, -normal(1)]
        e1 = e1/norm2(e1)
        e2 = [normal(2)*e1(3) - normal(3)*e1(2), normal(3)*e1(1) - normal(1)*e1(3), &
              normal(1)*e1(2) - normal(2)*e1(1)]
        do j = 1, np
            s = 2*pi*real(j - 1, dp)/np
            pts(:, j) = [1.05_dp, 0.0_dp, 0.08_dp] + 0.12_dp*(cos(s)*e1 + sin(s)*e2)
        end do
        call coil_set_from_polylines(coils%set, [np], pts, status=status)
        if (status /= 0) error stop 'tilted loop'
        coils%nphi = 128
    end subroutine

    subroutine exterior_oracle()
        integer, parameter :: modes(4) = [0, 1, 2, 3]
        integer, parameter :: levels(3) = [24, 48, 96]
        type(gamma_c_t) :: gc
        type(polyline_coils_t) :: coils
        complex(dp), allocatable :: u(:, :, :), h(:, :, :)
        real(dp), allocatable :: s(:, :), one(:)
        complex(dp), allocatable :: r(:)
        real(dp) :: error(size(modes), size(levels)), closure, rate
        integer :: l, m, j, status
        logical :: ok
        call tilted_loop(coils)
        do l = 1, size(levels)
            call ellipse(levels(l), gc)
            call gamma_c_coil_traces(gc, coils, modes, 12, u, h, status, closure)
            if (status /= GAMMA_C_OK) error stop 'coil traces'
            if (closure > 1.0e-10_dp) error stop 'n = 0 potential is not single valued'
            do m = 1, size(modes)
                call gamma_c_exterior(gc, GAMMA_C_SCALAR, modes(m), 2, s, status)
                if (status /= GAMMA_C_OK) error stop 'exterior map'
                ! The n = 0 trace is a line integral, fixed up to a constant
                ! that S_0 (the 3D capacitance) does not annihilate: fit it.
                r = matmul(s, u(:, m, 1)) + h(:, m, 1)
                if (modes(m) == 0) then
                    one = matmul(s, [(1.0_dp, j=1, size(r))])
                    r = r - one*dot_product(one, r)/dot_product(one, one)
                end if
                error(m, l) = norm2(abs(r))/norm2(abs(h(:, m, 1)))
            end do
        end do
        ok = .true.
        do m = 1, size(modes)
            rate = log(error(m, 2)/error(m, 3))/log(2.0_dp)
            print '(a,i2,a,3es11.3,a,f6.2)', 'n =', modes(m), &
                '  |S u + h|/|h| at 24/48/96 edges', error(m, :), '  rate', rate
            ok = ok .and. error(m, 3) < 1.0e-4_dp .and. rate > 2.5_dp
        end do
        if (.not. ok) error stop 'per-n exterior map against the exact field'
        print '(a)', 'PASS: per-n Steklov-Poincare map and coil traces on Gamma_c'
    end subroutine

    subroutine coil_oracle()
        integer, parameter :: np = 2048
        type(axisymmetric_coils_t) :: ring
        type(polyline_coils_t) :: polygon
        complex(dp) :: ba(3, 3, 4, 1), bp(3, 3, 4, 1)
        real(dp) :: pts(3, np), points(2, 4), s, scale, error, leak
        integer :: j, status
        ! Axisymmetric filament, 1 A, at R = 1.6, Z = 0.4.
        allocate (ring%set%coils(1), ring%set%incidence(1, 1))
        ring%set%coils(1)%kind = gs_coil_filaments
        ring%set%coils(1)%fil_r = [1.6_dp]
        ring%set%coils(1)%fil_z = [0.4_dp]
        ring%set%coils(1)%fil_turns = [1.0_dp]
        ring%set%incidence = 1
        ! The same circle as a polygon, current in +phi.
        do j = 1, np
            s = 2*pi*real(j - 1, dp)/np
            pts(:, j) = [1.6_dp*cos(s), 1.6_dp*sin(s), 0.4_dp]
        end do
        call coil_set_from_polylines(polygon%set, [np], pts, status=status)
        if (status /= 0) error stop 'polygon'
        polygon%nphi = 64
        points = reshape([1.0_dp, 0.0_dp, 1.3_dp, -0.3_dp, 0.8_dp, 0.5_dp, &
                          2.2_dp, 0.1_dp], [2, 4])
        call ring%field_modes(points, [0, 1, 2], ba, status)
        if (status /= GAMMA_C_OK) error stop 'axisymmetric modes'
        call polygon%field_modes(points, [0, 1, 2], bp, status)
        if (status /= GAMMA_C_OK) error stop 'polygon modes'
        scale = maxval(abs(ba(:, 1, :, 1)))
        error = maxval(abs(bp(:, 1, :, 1) - ba(:, 1, :, 1)))/scale
        leak = maxval(abs(bp(:, 2:3, :, 1)))/scale
        print '(a,es11.3,a,es11.3)', 'n = 0 polygon vs exact ring', error, &
            '  n /= 0 leak', leak
        if (error > 1.0e-5_dp) error stop 'coil sets disagree at n = 0'
        if (leak > 1.0e-10_dp .or. maxval(abs(ba(:, 2:3, :, 1))) > 0) &
            error stop 'axisymmetric coils have n /= 0 harmonics'
        print '(a)', 'PASS: axisymmetric and polyline coils through one interface'
    end subroutine
end program test_equilibrium_boundary
