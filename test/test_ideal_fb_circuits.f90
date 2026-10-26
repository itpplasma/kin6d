module ideal_fb_circuits_test_field
    ! Smooth nested frozen data (elliptic labels about (3.05, 0), F = 6 T m),
    ! the same simple data as the N4 derivative checks.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_ideal_map, only: ideal_frozen_t
    implicit none
    private
    type, extends(ideal_frozen_t), public :: elliptic_field_t
    contains
        procedure :: evaluate => elliptic_evaluate
    end type
contains
    subroutine elliptic_evaluate(field, point, chi, chi_gradient, b3, pressure)
        class(elliptic_field_t), intent(in) :: field
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: chi, chi_gradient(2), b3, pressure
        chi = 0.5_dp*((point(1) - 3.05_dp)**2 + (point(2)/1.5_dp)**2)
        chi_gradient = [point(1) - 3.05_dp, point(2)/1.5_dp**2]
        b3 = 6.0_dp/point(1)
        pressure = 1.0e4_dp*max(0.0_dp, 1 - chi/0.18_dp)
    end subroutine
end module ideal_fb_circuits_test_field

program test_ideal_fb_circuits
    ! Circuit currents as parameters of the nested-map free-boundary residual
    ! (kin6d_ideal_fb_circuits), synthetic public coils only, nfp = 3, N = 1:
    !   derivatives  JVP in the currents against central differences of the
    !                residual (exact for its quadratic current dependence),
    !                adjoint identity <v, J w> = <J^T v, w>, superposition of
    !                the vacuum traces and of the F_vac-free residual part,
    !                Ampere oracle for F_vac of the toroidal-field circuit.
    !   exterior3d   the same JVP/VJP checks with a frozen full-surface R2
    !                exterior attached (kin6d_equilibrium_boundary_3d form).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_equilibrium_boundary_3d, only: gamma_c_3d_t
    use kin6d_ideal_free_boundary, only: ideal_fb_t, ideal_fb_init, ideal_fb_gradient, &
        ideal_fb_step, ideal_fb_attach_exterior_3d, ideal_fb_full_turn_traces, &
        IDEAL_FB_OK
    use kin6d_ideal_fb_circuits, only: coil_circuit_basis_t, coil_circuit_basis_init, &
        ideal_fb_circuits_t, ideal_fb_circuits_build, ideal_fb_set_currents, &
        ideal_fb_currents_jvp, ideal_fb_currents_vjp
    use kin6d_coils_3d, only: coil_set_t, coil_set_from_polylines, coil_mu0
    use fortfem_axisymmetric_gs_kernel, only: gs_green_gradient
    use ideal_fb_circuits_test_field, only: elliptic_field_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = coil_mu0
    real(dp), parameter :: rc = 3.0_dp, ac = 1.5_dp, bc = 2.0_dp
    real(dp), parameter :: fil(2) = [3.05_dp, 0.05_dp], circulation = -mu0*1.0e6_dp
    integer, parameter :: ntf = 12, nc = 4
    ! Circuits: 1 PF pair (Z = +-2), 2 inner ring (2 turns), 3 window coils
    ! (n = 3), 4 toroidal field (12 planar coils, 3-periodic).
    real(dp), parameter :: p0(nc) = [4.0e5_dp, -1.5e5_dp, 1.0e5_dp, 2.5e6_dp]
    character(16) :: mode
    mode = 'derivatives'
    if (command_argument_count() > 0) call get_command_argument(1, mode)
    select case (trim(mode))
    case ('derivatives')
        call checks(.false.)
    case ('exterior3d')
        call checks(.true.)
    case default
        error stop 'unknown mode'
    end select
contains
    subroutine gamma_c(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [rc + ac*cos(2*pi*t), bc*sin(2*pi*t)]
    end subroutine

    subroutine filament(r, z, k, bcos, bsin)
        ! Fixed part: unit-circulation filament inside Gamma_p times c.
        real(dp), intent(in) :: r, z
        integer, intent(in) :: k
        real(dp), intent(out) :: bcos(3), bsin(3)
        real(dp) :: g, gg(4)
        bcos = 0; bsin = 0
        if (k /= 0) return
        call gs_green_gradient(r, z, fil(1), fil(2), g, gg)
        bcos([1, 3]) = circulation*[gg(2), -gg(1)]/r
    end subroutine

    subroutine synthetic_coils(basis)
        type(coil_circuit_basis_t), intent(out) :: basis
        type(coil_set_t) :: coils
        integer, parameter :: nring = 64
        real(dp), allocatable :: pts(:, :), turns(:)
        integer, allocatable :: npts(:), circuit(:)
        real(dp) :: corner(2, 4), phi0
        integer :: c, j, n, st, ncoil
        ncoil = 3 + 3 + ntf
        allocate (pts(3, 3*nring + 3*4 + ntf*4), npts(ncoil), circuit(ncoil), &
            turns(ncoil))
        n = 0
        do c = 1, 3
            do j = 1, nring
                phi0 = 2*pi*(j - 1)/nring
                select case (c)
                case (1)
                    pts(:, n + j) = [5.2_dp*cos(phi0), 5.2_dp*sin(phi0), 2.0_dp]
                case (2)
                    pts(:, n + j) = [5.2_dp*cos(phi0), 5.2_dp*sin(phi0), -2.0_dp]
                case (3)
                    pts(:, n + j) = [0.8_dp*cos(phi0), 0.8_dp*sin(phi0), 0.0_dp]
                end select
            end do
            n = n + nring; npts(c) = nring
        end do
        circuit(1:3) = [1, 1, 2]; turns(1:3) = [1.0_dp, 1.0_dp, 2.0_dp]
        corner = reshape([-0.35_dp, -0.8_dp, 0.35_dp, -0.8_dp, 0.35_dp, 0.8_dp, &
            -0.35_dp, 0.8_dp], [2, 4])
        do c = 1, 3
            phi0 = 2*pi*(c - 1)/3
            do j = 1, 4
                pts(:, n + j) = [4.9_dp*cos(phi0 + corner(1, j)/4.9_dp), &
                    4.9_dp*sin(phi0 + corner(1, j)/4.9_dp), corner(2, j) + 0.3_dp]
            end do
            n = n + 4; npts(3 + c) = 4
            circuit(3 + c) = 3; turns(3 + c) = 1
        end do
        ! Planar toroidal-field coils: current up the inner leg (R = 0.5),
        ! out along the top, down the outer leg (R = 6), so B_phi > 0 inside.
        corner = reshape([0.5_dp, -3.0_dp, 0.5_dp, 3.0_dp, 6.0_dp, 3.0_dp, &
            6.0_dp, -3.0_dp], [2, 4])
        do c = 1, ntf
            phi0 = 2*pi*(c - 0.5_dp)/ntf
            do j = 1, 4
                pts(:, n + j) = [corner(1, j)*cos(phi0), corner(1, j)*sin(phi0), &
                    corner(2, j)]
            end do
            n = n + 4; npts(6 + c) = 4
            circuit(6 + c) = 4; turns(6 + c) = 1
        end do
        call coil_set_from_polylines(coils, npts, pts, status=st)
        if (st /= 0) error stop 'synthetic coils'
        call coil_circuit_basis_init(basis, coils, circuit, turns, nc, 96, st)
        if (st /= 0) error stop 'circuit basis'
    end subroutine

    subroutine build(fb, circuits)
        type(ideal_fb_t), intent(out) :: fb
        type(ideal_fb_circuits_t), intent(out) :: circuits
        type(gs_state_t) :: mesh
        type(cylindrical_chart_t) :: chart
        type(elliptic_field_t) :: field
        type(coil_circuit_basis_t) :: basis
        real(dp), allocatable :: contour(:, :), chain(:, :)
        real(dp) :: h
        integer :: j, st
        allocate (contour(2, 16), chain(2, 12))
        do j = 1, 16
            call gamma_c(real(j - 1, dp)/16, contour(:, j))
        end do
        do j = 1, 12
            chain(:, j) = [3.05_dp + 0.6_dp*cos(2*pi*(j - 1)/12), &
                0.9_dp*sin(2*pi*(j - 1)/12)]
        end do
        h = 2*pi*ac/16
        call initialize_gs(mesh, contour, 0.5_dp*h*h, st, 2, curve=gamma_c, &
            interior_chain=chain)
        if (st /= 0) error stop 'mesh'
        call ideal_fb_init(fb, mesh, chain, chart, field, 3, 1, 6.0_dp, circulation, &
            st, mu0=mu0)
        if (st /= IDEAL_FB_OK) error stop 'fb init'
        call synthetic_coils(basis)
        call ideal_fb_circuits_build(fb, basis, circuits, st, background=filament)
        if (st /= IDEAL_FB_OK) error stop 'circuit traces'
        print '(a,4es11.3)', 'circuit closure (T m, T) ', circuits%closure
        print '(a,4es11.3)', 'circuit F_k (T m/A)      ', circuits%f_vac
        print '(a,4es11.3)', 'R B_phi^0 spread (T m/A)  ', circuits%f_spread
        ! Unlinked axisymmetric-type circuits: the n = 0 poloidal circulation
        ! on Gamma_c vanishes (Ampere). The window entry also holds the n = 3
        ! midpoint curl-free indicator, a P2 trace discretisation measure.
        if (any(circuits%closure([1, 2, 4]) > 1e-15_dp)) &
            error stop 'circuit trace closure'
    end subroutine

    subroutine residual(fb, circuits, p, g)
        type(ideal_fb_t), intent(inout) :: fb
        type(ideal_fb_circuits_t), intent(in) :: circuits
        real(dp), intent(in) :: p(:)
        real(dp), intent(out) :: g(:)
        integer :: st
        call ideal_fb_set_currents(fb, circuits, p, st)
        if (st /= IDEAL_FB_OK) error stop 'set currents'
        call ideal_fb_gradient(fb, g, st)
        if (st /= IDEAL_FB_OK) error stop 'residual'
    end subroutine

    subroutine attach_frozen_exterior(fb)
        ! Independent SPD full-surface operator with cross-mode coupling
        ! (a frozen R2 exterior form, no BEM assembly).
        type(ideal_fb_t), intent(inout) :: fb
        type(gamma_c_3d_t) :: boundary
        real(dp), allocatable :: v(:), d3(:), l3(:)
        integer :: n, i, st
        n = fb%nt*fb%nmode
        allocate (boundary%s(n, n), v(n), d3(n), l3(n))
        boundary%genus = 1
        v = [(sin(.37_dp*i), i=1, n)]/sqrt(real(n, dp))
        boundary%s = spread(v, 2, n)*spread(v, 1, n)
        do i = 1, n
            boundary%s(i, i) = boundary%s(i, i) + 1 + 0.1_dp*cos(real(i, dp))
        end do
        call ideal_fb_full_turn_traces(fb, fb%drive, fb%load, d3, l3)
        call ideal_fb_attach_exterior_3d(fb, boundary, d3, l3, st)
        if (st /= IDEAL_FB_OK) error stop 'frozen exterior attachment'
    end subroutine

    subroutine checks(exterior)
        logical, intent(in) :: exterior
        type(ideal_fb_t) :: fb
        type(ideal_fb_circuits_t) :: circuits
        real(dp), allocatable :: g0(:), gp(:), gm(:), jv(:), v(:), du(:), g1(:), g2(:)
        real(dp), allocatable :: g12(:), d1(:, :), d2(:, :), d12(:, :), drv0(:, :)
        real(dp), allocatable :: ld0(:, :)
        real(dp) :: dp_(nc), w(nc), jtv(nc), h, err, scale, lhs, rhs, f_tf, f1, f2, f12
        real(dp) :: p1(nc), p2(nc)
        integer :: i, k, st
        call build(fb, circuits)
        ! Ampere: F_k = mu0 (linked current)/(2 pi); TF circuit links 12 A/A.
        f_tf = mu0*ntf/(2*pi)
        print '(a,es12.4)', 'TF circuit F relative to Ampere ', &
            abs(circuits%f_vac(4) - f_tf)/f_tf
        if (abs(circuits%f_vac(4) - f_tf) > 1e-10_dp*f_tf) error stop 'Ampere F_vac'
        if (maxval(abs(circuits%f_vac(1:3))) > 1e-12_dp*f_tf) &
            error stop 'unlinked circuit F_vac'
        print '(a,es12.4)', 'F_vac(p0) (T m) ', dot_product(p0, circuits%f_vac)
        allocate (g0(fb%n), gp(fb%n), gm(fb%n), jv(fb%n), v(fb%n), du(fb%n), &
            g1(fb%n), g2(fb%n), g12(fb%n))
        call ideal_fb_set_currents(fb, circuits, p0, st)
        if (st /= IDEAL_FB_OK) error stop 'set currents p0'
        ! Generic state: perturbed map, multipliers and potential.
        do i = 1, fb%n
            du(i) = 1.0e-3_dp*sin(1.3_dp*i)
        end do
        du(fb%np + 2*fb%nx + 1:) = 1.0e2_dp*du(fb%np + 2*fb%nx + 1:)
        du(fb%np + fb%nx + 1:fb%np + 2*fb%nx) = &
            1.0e5_dp*du(fb%np + fb%nx + 1:fb%np + 2*fb%nx)
        call ideal_fb_step(fb, du, 1.0_dp)
        if (exterior) call attach_frozen_exterior(fb)
        call residual(fb, circuits, p0, g0)
        scale = maxval(abs(g0))
        ! (ii) JVP against central differences per circuit and along a mix.
        ! R is quadratic in p (F_vac enters the shell energy quadratically),
        ! so the central difference is exact up to roundoff.
        do k = 0, nc
            if (k == 0) then
                dp_ = p0*[0.3_dp, -0.7_dp, 1.1_dp, 0.2_dp]
            else
                dp_ = 0; dp_(k) = p0(k)
            end if
            h = 1.0e-2_dp
            call residual(fb, circuits, p0 + h*dp_, gp)
            call residual(fb, circuits, p0 - h*dp_, gm)
            call ideal_fb_set_currents(fb, circuits, p0, st)
            call ideal_fb_currents_jvp(fb, circuits, dp_, jv, st)
            if (st /= IDEAL_FB_OK) error stop 'jvp'
            err = maxval(abs((gp - gm)/(2*h) - jv))/max(maxval(abs(jv)), tiny(1.0_dp))
            print '(a,i0,a,es11.3,a,es11.3)', 'circuit ', k, ' |JVP - FD|/|JVP| ', &
                err, '  |JVP| ', maxval(abs(jv))
            if (err > 1e-8_dp) error stop 'JVP against central differences'
        end do
        ! Adjoint identity <v, J w> = <J^T v, w>.
        v = [(cos(0.71_dp*i), i=1, fb%n)]
        w = [0.4_dp, -1.3_dp, 0.9_dp, 0.05_dp]*p0
        call ideal_fb_currents_jvp(fb, circuits, w, jv, st)
        call ideal_fb_currents_vjp(fb, circuits, v, jtv, st)
        if (st /= IDEAL_FB_OK) error stop 'vjp'
        lhs = dot_product(v, jv); rhs = dot_product(jtv, w)
        err = abs(lhs - rhs)/max(abs(lhs), abs(rhs))
        print '(a,2es24.15,a,es10.2)', '<v,Jw> <J^T v,w> ', lhs, rhs, ' rel ', err
        if (err > 1e-12_dp) error stop 'adjoint identity'
        ! (iii) Superposition. Vacuum data are linear in p.
        p1 = [1.0_dp, 0.0_dp, -2.0_dp, 0.3_dp]*p0
        p2 = [-0.5_dp, 1.0_dp, 0.7_dp, 0.4_dp]*p0
        call ideal_fb_set_currents(fb, circuits, 0*p0, st); drv0 = fb%drive
        ld0 = fb%load
        call ideal_fb_set_currents(fb, circuits, p1, st); d1 = fb%drive; f1 = fb%f_vac
        g1(1:size(fb%load)) = reshape(fb%load, [size(fb%load)])
        call ideal_fb_set_currents(fb, circuits, p2, st); d2 = fb%drive; f2 = fb%f_vac
        g2(1:size(fb%load)) = reshape(fb%load, [size(fb%load)])
        call ideal_fb_set_currents(fb, circuits, p1 + p2, st); d12 = fb%drive
        f12 = fb%f_vac
        g12(1:size(fb%load)) = reshape(fb%load, [size(fb%load)])
        err = maxval(abs(d12 - drv0 - (d1 - drv0) - (d2 - drv0))) &
            /maxval(abs(d12 - drv0))
        k = size(fb%load)
        g0(:k) = reshape(ld0, [k])
        err = max(err, maxval(abs(g12(:k) - g1(:k) - g2(:k) + g0(:k))) &
            /maxval(abs(g12(:k) - g0(:k))), abs(f12 - f1 - f2)/abs(f12))
        print '(a,es11.3)', 'trace superposition defect ', err
        if (err > 1e-13_dp) error stop 'vacuum trace superposition'
        ! Residual: affine in the currents of circuits with F_k = 0, so the
        ! second difference along PF/inner/window directions vanishes.
        p1 = [0.5_dp, 0.0_dp, -1.0_dp, 0.0_dp]*p0
        p2 = [0.0_dp, 0.8_dp, 0.6_dp, 0.0_dp]*p0
        call residual(fb, circuits, p0, g0)
        call residual(fb, circuits, p0 + p1, g1)
        call residual(fb, circuits, p0 + p2, g2)
        call residual(fb, circuits, p0 + p1 + p2, g12)
        err = maxval(abs(g12 - g1 - g2 + g0))/maxval(abs(g12 - g0))
        print '(a,es11.3)', 'residual superposition defect (F_k = 0 circuits) ', err
        if (err > 1e-11_dp) error stop 'residual superposition'
        ! F_vac direction: R is exactly quadratic, third difference vanishes.
        dp_ = 0; dp_(4) = 0.3_dp*p0(4)
        call residual(fb, circuits, p0 + dp_, g1)
        call residual(fb, circuits, p0 + 2*dp_, g2)
        call residual(fb, circuits, p0 + 3*dp_, g12)
        err = maxval(abs(g12 - 3*g2 + 3*g1 - g0))/maxval(abs(g12 - g0))
        print '(a,es11.3)', 'residual third difference along TF current ', err
        if (err > 1e-10_dp) error stop 'quadratic F_vac dependence'
        print '(a)', 'circuit parameter derivatives ok'
    end subroutine
end program test_ideal_fb_circuits
