module ideal_fb_traces3d_test_support
    ! Frozen data of test_ideal_fb_circuits, in (R, Z) for the cylindrical
    ! chart and pulled back to the normal plane (q1, q2) = (Z, R - 3) of an
    ! nfp = 3 helical-axis curve chart, so that the coupling surface
    ! H(gamma x [0, 2 pi)) is not axisymmetric. The curve chart is the real
    ! plasma chart of the nested-map lane (kin6d_equilibrium_chart), not a
    ! test extension of the cylindrical chart.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_ideal_map, only: ideal_frozen_t
    use kin6d_equilibrium_chart, only: curve_chart_t
    implicit none
    private
    public :: helical_chart
    type, extends(ideal_frozen_t), public :: elliptic_field_t
    contains
        procedure :: evaluate => elliptic_evaluate
    end type
    type, extends(ideal_frozen_t), public :: elliptic_curve_field_t
    contains
        procedure :: evaluate => elliptic_curve_evaluate
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

    subroutine elliptic_curve_evaluate(field, point, chi, chi_gradient, b3, pressure)
        ! Orientation-changing pullback (R, Z) = (3 + q2, q1), as
        ! relabelled_landreman_t of test_curve_chart_solver.
        class(elliptic_curve_field_t), intent(in) :: field
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: chi, chi_gradient(2), b3, pressure
        type(elliptic_field_t) :: original
        real(dp) :: gradient(2)
        call original%evaluate([3 + point(2), point(1)], chi, gradient, b3, pressure)
        chi = -chi; chi_gradient = -gradient([2, 1])
    end subroutine

    subroutine helical_chart(chart, eps)
        ! Axis R = 3 + eps cos(3 phi), Z = eps sin(3 phi) as full-turn
        ! Cartesian Fourier modes, default (z) guide.
        type(curve_chart_t), intent(out) :: chart
        real(dp), intent(in) :: eps
        real(dp) :: c(3, 0:4), s(3, 0:4)
        integer :: st
        c = 0; s = 0
        c(1, 1) = 3; s(2, 1) = 3
        c(1, 2) = eps/2; s(2, 2) = -eps/2
        c(1, 4) = eps/2; s(2, 4) = eps/2
        s(3, 3) = eps
        call chart%init(c, s, 3, st)
        if (st /= 0) error stop 'helical curve chart'
    end subroutine
end module ideal_fb_traces3d_test_support

program test_ideal_fb_traces3d
    ! Vacuum traces on a general coupling surface Gamma_c through the
    ! Cartesian basis contract (ideal_fb_drive_traces_3d), synthetic public
    ! coils only (the circuits of test_ideal_fb_circuits), nfp = 3, N = 1:
    !   axisymmetric  cylindrical chart: per-circuit and background full-turn
    !                 traces, F_k and closures equal those of the meridional
    !                 path (ideal_fb_drive_traces) to 1e-12.
    !   helical       nfp = 3 helical-axis curve chart (non-axisymmetric
    !                 Gamma_c, ideal_fb_init on the real chart): Gauss
    !                 (net flux), Ampere (F_k), superposition of circuits,
    !                 thread invariance, and JVP/adjoint identities of the
    !                 residual with a frozen R2 exterior.
    !   makegrid      the same coils written as a MAKEGRID file with groups
    !                 and read through libneo (neo_makegrid_coils): circuits,
    !                 turns, identical traces, and kin6d Biot-Savart of each
    !                 group against libneo's neo_biotsavart.
    !   scaling       (not a CTest) trace-build strong scaling over threads.
    !   exterior      R2 exterior of the helical Gamma_c: the field-period
    !                 (block-circulant) build against the exact full-surface
    !                 DtN of the same circulant operators (roundoff), against
    !                 the dense all-pairs build (pair quadrature tolerance),
    !                 symmetry and semidefiniteness, and thread invariance.
    !   exterior-scaling (not a CTest) cost and strong scaling of both builds.
    !   curve-hessian gradient and assembled KKT matrix of the residual on
    !                 the helical curve chart (N = 1, 2) against central
    !                 differences, with a frozen R2 exterior and coil traces.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    !$ use omp_lib, only: omp_set_num_threads, omp_get_max_threads, omp_get_wtime
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs
    use kin6d_equilibrium_chart, only: cylindrical_chart_t, curve_chart_t, &
        ideal_chart_t
    use kin6d_equilibrium_boundary_3d, only: gamma_c_3d_t
    use kin6d_ideal_free_boundary_bem3d, only: ideal_fb_use_bem3d, ideal_fb_surface, &
        coupling_surface_t
    use fortfem_parametric_bem_high_order, only: bem_pair_options_t, &
        bem_pair_report_t, assemble_parametric_bem_p2_circulant
    use fortfem_parametric_laplace_bem_3d, only: bem_operators_t, &
        bem_exterior_dtn_t, build_exterior_dtn
    use fortfem_bem_circulant, only: bem_circulant_t
    use kin6d_ideal_free_boundary, only: ideal_fb_t, ideal_fb_init, ideal_fb_gradient, &
        ideal_fb_hessian, &
        ideal_fb_step, ideal_fb_attach_exterior_3d, ideal_fb_full_turn_traces, &
        IDEAL_FB_OK
    use kin6d_ideal_fb_circuits, only: coil_circuit_basis_t, coil_circuit_basis_init, &
        ideal_fb_circuits_t, ideal_fb_circuits_build, ideal_fb_set_currents, &
        ideal_fb_currents_jvp, ideal_fb_currents_vjp
    use kin6d_coils_3d, only: coil_set_t, coil_set_from_polylines, coil_mu0, segment_b
    use kin6d_makegrid_circuits, only: makegrid_circuits_read
    use neo_biotsavart, only: coils_t, compute_magnetic_field, clight
    use neo_makegrid_coils, only: makegrid_coils_t, read_makegrid_coils, &
        makegrid_to_biotsavart
    use fortfem_axisymmetric_gs_kernel, only: gs_green_gradient
    use ideal_fb_traces3d_test_support, only: elliptic_field_t, &
        elliptic_curve_field_t, helical_chart
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = coil_mu0
    real(dp), parameter :: rc = 3.0_dp, ac = 1.5_dp, bc = 2.0_dp
    real(dp), parameter :: fil(2) = [3.05_dp, 0.05_dp], circulation = -mu0*1.0e6_dp
    integer, parameter :: ntf = 12, nc = 4, nphi = 96
    real(dp), parameter :: p0(nc) = [4.0e5_dp, -1.5e5_dp, 1.0e5_dp, 2.5e6_dp]
    character(len=4), parameter :: names(nc) = ['PF  ', 'RING', 'WIN ', 'TF  ']
    character(16) :: mode
    mode = 'axisymmetric'
    if (command_argument_count() > 0) call get_command_argument(1, mode)
    select case (trim(mode))
    case ('axisymmetric')
        call axisymmetric
    case ('helical')
        call helical
    case ('makegrid')
        call makegrid
    case ('scaling')
        call scaling
    case ('curve-hessian')
        call curve_derivatives
    case ('exterior')
        call exterior(6, 9, .true.)
    case ('exterior-scaling')
        call exterior(12, 18, .false.)
    case default
        error stop 'unknown mode'
    end select
contains
    subroutine gamma_c(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [rc + ac*cos(2*pi*t), bc*sin(2*pi*t)]
    end subroutine

    subroutine gamma_c_curve(t, point)
        ! The ellipse of gamma_c in (q1, q2) = (Z, R - 3), counter-clockwise.
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [bc*cos(2*pi*t), ac*sin(2*pi*t)]
    end subroutine

    subroutine gamma_c_shifted(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        call gamma_c_curve(t, point)
        point(1) = point(1) + 3
    end subroutine

    subroutine filament(r, z, k, bcos, bsin)
        real(dp), intent(in) :: r, z
        integer, intent(in) :: k
        real(dp), intent(out) :: bcos(3), bsin(3)
        real(dp) :: g, gg(4)
        bcos = 0; bsin = 0
        if (k /= 0) return
        call gs_green_gradient(r, z, fil(1), fil(2), g, gg)
        bcos([1, 3]) = circulation*[gg(2), -gg(1)]/r
    end subroutine

    subroutine coil_polylines(pts, npts, circuit, turns, nring)
        ! test_ideal_fb_circuits coils: 1 PF pair (Z = +-2), 2 inner ring
        ! (2 turns), 3 window coils (n = 3), 4 toroidal field (12 planar).
        real(dp), allocatable, intent(out) :: pts(:, :), turns(:)
        integer, allocatable, intent(out) :: npts(:), circuit(:)
        integer, intent(in) :: nring
        real(dp) :: corner(2, 4), phi0, rr, zz
        integer :: c, j, n, ncoil
        ncoil = 3 + 3 + ntf
        allocate (pts(3, 3*nring + 3*4 + ntf*4), npts(ncoil), circuit(ncoil), &
            turns(ncoil))
        n = 0
        do c = 1, 3
            rr = merge(0.8_dp, 5.2_dp, c == 3)
            zz = merge(2.0_dp, -2.0_dp, c == 1); if (c == 3) zz = 0
            do j = 1, nring
                phi0 = 2*pi*(j - 1)/nring
                pts(:, n + j) = [rr*cos(phi0), rr*sin(phi0), zz]
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
    end subroutine

    subroutine synthetic_basis(basis, nring, weights)
        ! weights present: one circuit carrying sum_k weights(k) b_k.
        type(coil_circuit_basis_t), intent(out) :: basis
        integer, intent(in) :: nring
        real(dp), intent(in), optional :: weights(nc)
        type(coil_set_t) :: coils
        real(dp), allocatable :: pts(:, :), turns(:)
        integer, allocatable :: npts(:), circuit(:)
        integer :: st
        call coil_polylines(pts, npts, circuit, turns, nring)
        call coil_set_from_polylines(coils, npts, pts, status=st)
        if (st /= 0) error stop 'synthetic coils'
        if (present(weights)) then
            turns = turns*weights(circuit); circuit = 1
            call coil_circuit_basis_init(basis, coils, circuit, turns, 1, nphi, st)
        else
            call coil_circuit_basis_init(basis, coils, circuit, turns, nc, nphi, st)
        end if
        if (st /= 0) error stop 'circuit basis'
    end subroutine

    subroutine build_fb(fb, chart, ncontour, ntor)
        ! Cylindrical chart: Gamma_c and Gamma_p in (R, Z). Curve chart: the
        ! same ellipses in the normal plane (q1, q2) = (Z, R - 3), counter-
        ! clockwise, with the pulled-back frozen field.
        type(ideal_fb_t), intent(out) :: fb
        class(ideal_chart_t), intent(in) :: chart
        integer, intent(in), optional :: ncontour, ntor
        type(gs_state_t) :: mesh
        type(elliptic_field_t) :: field
        type(elliptic_curve_field_t) :: curve_field
        real(dp), allocatable :: contour(:, :), chain(:, :)
        real(dp) :: h, t
        integer :: j, st, nb, nt
        logical :: curve
        nb = 16; nt = 1
        if (present(ntor)) nt = ntor
        if (present(ncontour)) nb = ncontour
        curve = .false.
        select type (chart)
        type is (curve_chart_t)
            curve = .true.
        end select
        allocate (contour(2, nb), chain(2, 12))
        do j = 1, nb
            if (curve) then
                call gamma_c_curve(real(j - 1, dp)/nb, contour(:, j))
            else
                call gamma_c(real(j - 1, dp)/nb, contour(:, j))
            end if
        end do
        do j = 1, 12
            t = 2*pi*(j - 1)/12
            chain(:, j) = [3.05_dp + 0.6_dp*cos(t), 0.9_dp*sin(t)]
            if (curve) chain(:, j) = [0.9_dp*cos(t), 0.05_dp + 0.6_dp*sin(t)]
        end do
        h = 2*pi*ac/nb
        if (curve) then
            ! initialize_gs meshes R > 0 planes: mesh the normal plane shifted
            ! by q1 + 3, then translate the nodes back.
            contour(1, :) = contour(1, :) + 3; chain(1, :) = chain(1, :) + 3
            call initialize_gs(mesh, contour, 0.5_dp*h*h, st, 2, &
                curve=gamma_c_shifted, interior_chain=chain)
            if (st /= 0) error stop 'mesh'
            mesh%points(1, :) = mesh%points(1, :) - 3; chain(1, :) = chain(1, :) - 3
            call ideal_fb_init(fb, mesh, chain, chart, curve_field, 3, nt, 6.0_dp, &
                circulation, st, mu0=mu0, nphi=8)
        else
            call initialize_gs(mesh, contour, 0.5_dp*h*h, st, 2, curve=gamma_c, &
                interior_chain=chain)
            if (st /= 0) error stop 'mesh'
            call ideal_fb_init(fb, mesh, chain, chart, field, 3, nt, 6.0_dp, &
                circulation, st, mu0=mu0)
        end if
        if (st /= IDEAL_FB_OK) error stop 'fb init'
    end subroutine

    real(dp) function full_turn_defect(fb, d1, l1, d2, l2) result(err)
        ! Relative difference of two trace pairs in the R2 full-turn pairing.
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(in) :: d1(:, :), l1(:, :), d2(:, :), l2(:, :)
        real(dp), allocatable :: a1(:), b1(:), a2(:), b2(:)
        integer :: n
        n = fb%nt*fb%nmode
        allocate (a1(n), b1(n), a2(n), b2(n))
        call ideal_fb_full_turn_traces(fb, d1, l1, a1, b1)
        call ideal_fb_full_turn_traces(fb, d2, l2, a2, b2)
        err = max(norm2(a1 - a2)/max(norm2(a1), tiny(1.0_dp)), &
            norm2(b1 - b2)/max(norm2(b1), tiny(1.0_dp)))
    end function

    subroutine axisymmetric
        type(ideal_fb_t) :: fb
        type(cylindrical_chart_t) :: chart
        type(coil_circuit_basis_t) :: basis
        type(ideal_fb_circuits_t) :: c2, c3
        real(dp) :: err, worst
        integer :: k, st
        call build_fb(fb, chart)
        call synthetic_basis(basis, 64)
        call ideal_fb_circuits_build(fb, basis, c2, st, background=filament)
        if (st /= IDEAL_FB_OK) error stop 'meridional traces'
        call ideal_fb_circuits_build(fb, basis, c3, st, background=filament, &
            nphi_surface=nphi)
        if (st /= IDEAL_FB_OK) error stop 'surface traces'
        worst = 0
        do k = 1, nc
            err = full_turn_defect(fb, c2%drive(:, :, k), c2%load(:, :, k), &
                c3%drive(:, :, k), c3%load(:, :, k))
            ! The TF traces vanish on an axisymmetric Gamma_c (|.| ~ 1e-22,
            ! roundoff), so defects are relative to the largest circuit.
            err = err*max(norm2(c2%drive(:, :, k)), norm2(c2%load(:, :, k))) &
                /max(maxval(norm2(reshape(c2%drive, [fb%nt*fb%nmode, nc]), 1)), &
                maxval(norm2(reshape(c2%load, [fb%nt*fb%nmode, nc]), 1)))
            print '(a,i0,a,es10.2)', 'circuit ', k, &
                ' full-turn traces 3D vs meridional ', err
            worst = max(worst, err)
        end do
        err = full_turn_defect(fb, c2%drive0, c2%load0, c3%drive0, c3%load0)
        print '(a,es10.2)', 'background full-turn traces 3D vs meridional ', err
        worst = max(worst, err)
        err = maxval(abs(c2%f_vac - c3%f_vac))/maxval(abs(c2%f_vac))
        print '(a,es10.2)', 'F_k 3D vs meridional ', err
        worst = max(worst, err)
        print '(a,4es10.2,a,4es10.2)', 'closure meridional ', c2%closure, ' 3D ', &
            c3%closure
        if (any(abs(c2%closure - c3%closure) &
                > 1e-12_dp*max(abs(c2%closure), 1e-3_dp))) &
            error stop 'closure 3D vs meridional'
        if (worst > 1e-12_dp) error stop 'axisymmetric 3D traces vs meridional path'
        call derivative_checks(fb, c3)
        print '(a)', 'axisymmetric surface traces ok'
    end subroutine

    subroutine helical
        type(ideal_fb_t) :: fb, fb_axi
        type(curve_chart_t) :: chart
        type(cylindrical_chart_t) :: cyl
        type(coil_circuit_basis_t) :: basis, sum_basis
        type(ideal_fb_circuits_t) :: c3, c3n, csum, caxi
        real(dp) :: f_tf, err, w(nc), dsum(1), flux
        integer :: k, st, nthreads
        call helical_chart(chart, 0.15_dp)
        call build_fb(fb, chart)
        if (allocated(fb%s_trace)) error stop 'curve chart got the per-n exterior'
        block
            real(dp), allocatable :: g(:)
            allocate (g(fb%n))
            call ideal_fb_gradient(fb, g, st)
            if (st == IDEAL_FB_OK) error stop 'residual without an exterior'
        end block
        call build_fb(fb_axi, cyl)
        call synthetic_basis(basis, 64)
        !$ call omp_set_num_threads(1)
        call ideal_fb_circuits_build(fb, basis, c3, st, background=filament, &
            nphi_surface=nphi)
        if (st /= IDEAL_FB_OK) error stop 'helical surface traces'
        nthreads = 1
        !$ nthreads = 8
        !$ call omp_set_num_threads(nthreads)
        call ideal_fb_circuits_build(fb, basis, c3n, st, background=filament, &
            nphi_surface=nphi)
        if (st /= IDEAL_FB_OK) error stop 'helical surface traces (threads)'
        if (any(c3n%drive /= c3%drive) .or. any(c3n%load /= c3%load) .or. &
            any(c3n%f_vac /= c3%f_vac) .or. any(c3n%drive0 /= c3%drive0)) &
            error stop 'thread count changed the traces'
        print '(a,i0,a)', 'traces at 1 and ', nthreads, ' threads bitwise equal'
        ! The surface is not axisymmetric: traces differ from those of the
        ! cylindrical chart on the same reference boundary.
        call ideal_fb_circuits_build(fb_axi, basis, caxi, st, nphi_surface=nphi)
        err = full_turn_defect(fb, c3%drive(:, :, 4), c3%load(:, :, 4), &
            caxi%drive(:, :, 4), caxi%load(:, :, 4))
        print '(a,es10.2)', 'TF traces helical vs axisymmetric Gamma_c ', err
        if (err < 1e-3_dp) error stop 'helical surface indistinguishable'
        ! Gauss: the n = 0 load summed over the nodes is the net flux / 2 pi.
        ! The window coils (0.25 m from the helical Gamma_c, 0.14 rad wide)
        ! are under-resolved by the nphi = 96 trapezoid rule; their flux
        ! defect has to fall with nphi instead (quadrature, not geometry).
        do k = 1, nc
            flux = sum(c3%load(:, 1, k))/sum(abs(c3%load(:, 1, k)))
            print '(a,i0,a,es10.2)', 'circuit ', k, ' net flux / sum |load_0| ', flux
            if (k /= 3 .and. abs(flux) > 1e-12_dp) &
                error stop 'Gauss (net flux) on helical Gamma_c'
        end do
        ! Ampere: F_k = mu0 (linked current)/(2 pi), whatever the surface.
        ! The curve-chart angle q3 is not the cylindrical angle, so the
        ! nphi-point trapezoid rule over the TF ripple is exact only to its
        ! quadrature error (1e-11 at nphi = 96, roundoff at 192): Ampere is
        ! checked on the converged nphi = 384 build of the window check.
        f_tf = mu0*ntf/(2*pi)
        print '(a,4es11.3)', 'F_k (T m/A)        ', c3%f_vac
        print '(a,4es11.3)', 'F spread (T m/A)   ', c3%f_spread
        print '(a,4es11.3)', 'closure (T m, T)   ', c3%closure
        block
            type(ideal_fb_circuits_t) :: cn
            real(dp) :: f96
            f96 = abs(sum(c3%load(:, 1, 3)))/sum(abs(c3%load(:, 1, 3)))
            call ideal_fb_circuits_build(fb, basis, cn, st, nphi_surface=4*nphi)
            flux = abs(sum(cn%load(:, 1, 3)))/sum(abs(cn%load(:, 1, 3)))
            print '(a,i0,a,es10.2)', 'window net flux at nphi ', 4*nphi, ' ', flux
            if (flux > 1e-4_dp*f96) error stop 'window flux does not converge in nphi'
            print '(a,i0,a,es11.3)', 'TF F relative to Ampere (nphi ', 4*nphi, ') ', &
                abs(cn%f_vac(4) - f_tf)/f_tf
            print '(a,4es11.3)', 'F spread (nphi 384) ', cn%f_spread
            if (abs(cn%f_vac(4) - f_tf) > 1e-10_dp*f_tf) error stop 'Ampere F_vac'
            if (maxval(abs(cn%f_vac(1:2))) > 1e-12_dp*f_tf) error stop 'unlinked F_vac'
            if (maxval(cn%f_spread([1, 2, 4])) > 1e-10_dp*f_tf) error stop 'F_vac spread'
        end block
        ! closure = max(n = 0 Ampere closure, n = 3 curl-free indicator).
        ! On the helical surface every circuit has n = 3 content in B.x_t, so
        ! the indicator is the P2 trace discretisation error: it must fall
        ! with refinement, while the n = 0 closure alone (N = 0) vanishes.
        block
            type(ideal_fb_t) :: other
            type(ideal_fb_circuits_t) :: cf
            call build_fb(other, chart, 32)
            call ideal_fb_circuits_build(other, basis, cf, st, nphi_surface=nphi)
            print '(a,4es11.3)', 'closure, 32 edges   ', cf%closure
            if (any(cf%closure > 0.6_dp*c3%closure)) &
                error stop 'curl-free indicator does not fall with refinement'
            call build_fb(other, chart, 16, 0)
            call ideal_fb_circuits_build(other, basis, cf, st, nphi_surface=nphi)
            print '(a,4es11.3)', 'closure, n = 0 only ', cf%closure
            if (any(cf%closure([1, 2, 4]) > 1e-15_dp)) error stop 'n = 0 trace closure'
        end block
        ! Superposition: one circuit carrying sum_k w_k b_k against sum_k w_k traces_k.
        w = [1.0_dp, -2.0_dp, 0.7_dp, 0.3_dp]
        call synthetic_basis(sum_basis, 64, w)
        call ideal_fb_circuits_build(fb, sum_basis, csum, st, nphi_surface=nphi)
        if (st /= IDEAL_FB_OK) error stop 'sum traces'
        err = 0
        dsum = 0
        block
            real(dp), allocatable :: d(:, :), l(:, :)
            d = 0*csum%drive(:, :, 1); l = d
            do k = 1, nc
                d = d + w(k)*c3%drive(:, :, k); l = l + w(k)*c3%load(:, :, k)
            end do
            err = full_turn_defect(fb, csum%drive(:, :, 1), csum%load(:, :, 1), d, l)
            dsum = abs(csum%f_vac(1) - dot_product(w, c3%f_vac))/abs(csum%f_vac(1))
        end block
        err = max(err, dsum(1))
        print '(a,es10.2)', 'circuit superposition defect (helical Gamma_c) ', err
        if (err > 1e-14_dp) error stop 'superposition of circuit traces'
        call derivative_checks(fb, c3)
        print '(a)', 'helical surface traces ok'
    end subroutine

    subroutine attach_frozen_exterior(fb)
        ! Independent SPD full-surface operator (frozen R2 exterior form).
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

    subroutine derivative_checks(fb, circuits)
        ! JVP in the currents against central differences of the residual
        ! (exact for its quadratic current dependence) and the adjoint
        ! identity, with a frozen R2 exterior, as in test_ideal_fb_circuits.
        type(ideal_fb_t), intent(inout) :: fb
        type(ideal_fb_circuits_t), intent(in) :: circuits
        real(dp), allocatable :: gp(:), gm(:), jv(:), v(:), du(:)
        real(dp) :: dp_(nc), w(nc), jtv(nc), h, err, lhs, rhs
        integer :: i, k, st
        allocate (gp(fb%n), gm(fb%n), jv(fb%n), v(fb%n), du(fb%n))
        call ideal_fb_set_currents(fb, circuits, p0, st)
        if (st /= IDEAL_FB_OK) error stop 'set currents p0'
        do i = 1, fb%n
            du(i) = 1.0e-3_dp*sin(1.3_dp*i)
        end do
        du(fb%np + 2*fb%nx + 1:) = 1.0e2_dp*du(fb%np + 2*fb%nx + 1:)
        du(fb%np + fb%nx + 1:fb%np + 2*fb%nx) = &
            1.0e5_dp*du(fb%np + fb%nx + 1:fb%np + 2*fb%nx)
        call ideal_fb_step(fb, du, 1.0_dp)
        call attach_frozen_exterior(fb)
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
            print '(a,i0,a,es11.3)', 'circuit ', k, ' |JVP - FD|/|JVP| ', err
            if (err > 1e-8_dp) error stop 'JVP against central differences'
        end do
        v = [(cos(0.71_dp*i), i=1, fb%n)]
        w = [0.4_dp, -1.3_dp, 0.9_dp, 0.05_dp]*p0
        call ideal_fb_currents_jvp(fb, circuits, w, jv, st)
        call ideal_fb_currents_vjp(fb, circuits, v, jtv, st)
        if (st /= IDEAL_FB_OK) error stop 'vjp'
        lhs = dot_product(v, jv); rhs = dot_product(jtv, w)
        err = abs(lhs - rhs)/max(abs(lhs), abs(rhs))
        print '(a,es10.2)', '<v,Jw> vs <J^T v,w> rel ', err
        if (err > 1e-12_dp) error stop 'adjoint identity'
    end subroutine

    subroutine write_makegrid(path, nring)
        ! W7-X style: fourth column = winding count, group on the terminator.
        character(len=*), intent(in) :: path
        integer, intent(in) :: nring
        real(dp), allocatable :: pts(:, :), turns(:)
        integer, allocatable :: npts(:), circuit(:)
        integer :: unit, c, j, n
        call coil_polylines(pts, npts, circuit, turns, nring)
        open (newunit=unit, file=path, status='replace', action='write')
        write (unit, '(a)') 'periods 1'
        write (unit, '(a)') 'begin filament'
        write (unit, '(a)') 'mirror NIL'
        n = 0
        do c = 1, size(npts)
            do j = 1, npts(c)
                write (unit, '(4es25.16e3)') pts(:, n + j), turns(c)
            end do
            write (unit, '(4es25.16e3,1x,i0,1x,a)') pts(:, n + 1), 0.0_dp, &
                circuit(c), trim(names(circuit(c)))
            n = n + npts(c)
        end do
        write (unit, '(a)') 'end'
        close (unit)
    end subroutine

    subroutine makegrid
        character(len=*), parameter :: path = 'ideal_fb_traces3d_coils.synthetic'
        type(ideal_fb_t) :: fb
        type(curve_chart_t) :: chart
        type(coil_set_t) :: coils
        type(coil_circuit_basis_t) :: direct, read_basis
        type(ideal_fb_circuits_t) :: cd, cr
        type(makegrid_coils_t) :: mc
        type(coils_t) :: group_polygon
        real(dp), allocatable :: turns(:), pts(:, :), turns0(:)
        integer, allocatable :: circuit(:), npts(:), circuit0(:)
        character(len=16), allocatable :: group_name(:)
        real(dp) :: x(3), b_kin6d(3), b_libneo(3), err, worst
        integer :: ncircuit, st, k, c, s, i
        call write_makegrid(path, 64)
        call makegrid_circuits_read(path, coils, circuit, turns, ncircuit, st, &
            group_name=group_name)
        if (st /= 0) error stop 'makegrid read'
        call coil_polylines(pts, npts, circuit0, turns0, 64)
        if (ncircuit /= nc .or. coils%ncoil /= size(npts)) error stop 'circuit count'
        if (any(circuit /= circuit0) .or. any(turns /= turns0)) &
            error stop 'circuits or winding counts'
        do k = 1, nc
            if (trim(group_name(k)) /= trim(names(k))) error stop 'group names'
        end do
        print '(a,i0,a,i0,a)', 'MAKEGRID: ', coils%ncoil, ' coils in ', ncircuit, &
            ' circuits, winding counts and names ok'
        ! Traces from the file equal those from the direct construction.
        call helical_chart(chart, 0.15_dp)
        call build_fb(fb, chart)
        call synthetic_basis(direct, 64)
        call coil_circuit_basis_init(read_basis, coils, circuit, turns, ncircuit, &
            nphi, st)
        if (st /= 0) error stop 'basis from MAKEGRID'
        call ideal_fb_circuits_build(fb, direct, cd, st, nphi_surface=nphi)
        call ideal_fb_circuits_build(fb, read_basis, cr, st, nphi_surface=nphi)
        if (st /= IDEAL_FB_OK) error stop 'traces from MAKEGRID'
        err = max(maxval(abs(cd%drive - cr%drive))/maxval(abs(cd%drive)), &
            maxval(abs(cd%load - cr%load))/maxval(abs(cd%load)), &
            maxval(abs(cd%f_vac - cr%f_vac))/maxval(abs(cd%f_vac)))
        print '(a,es10.2)', 'traces MAKEGRID vs direct coils ', err
        if (err > 1e-14_dp) error stop 'MAKEGRID traces'
        ! Independent Biot-Savart: libneo (Gaussian units, current/c) of each
        ! group polygon against kin6d segments with winding counts.
        call read_makegrid_coils(path, mc)
        worst = 0
        do k = 1, nc
            call makegrid_to_biotsavart(mc, group_polygon, group=k)
            do i = 1, 7
                x = [3.0_dp + 1.2_dp*cos(1.1_dp*i), 0.7_dp*sin(0.9_dp*i), &
                    1.5_dp*sin(1.7_dp*i)]
                b_kin6d = 0
                do c = 1, coils%ncoil
                    if (circuit(c) /= k) cycle
                    do s = coils%first(c), coils%first(c + 1) - 1
                        b_kin6d = b_kin6d + turns(c)*segment_b(coils%a(:, s), &
                            coils%l(:, s), x)
                    end do
                end do
                b_libneo = 1.0e-7_dp*clight*compute_magnetic_field(group_polygon, x)
                worst = max(worst, norm2(b_kin6d - b_libneo)/norm2(b_kin6d))
            end do
        end do
        print '(a,es10.2)', 'group field kin6d vs libneo Biot-Savart ', worst
        if (worst > 1e-12_dp) error stop 'Biot-Savart kin6d vs libneo'
        print '(a)', 'MAKEGRID circuits ok'
    end subroutine

    subroutine scaling
        ! Strong scaling of the per-circuit trace build on the helical
        ! surface (finer coils, more toroidal points); results must be
        ! bitwise equal for every thread count.
        type(ideal_fb_t) :: fb
        type(curve_chart_t) :: chart
        type(coil_circuit_basis_t) :: basis
        type(ideal_fb_circuits_t) :: ref, run
        integer, parameter :: counts(5) = [1, 2, 4, 8, 16]
        real(dp) :: t0, t1, t_ref
        integer :: i, st, rep
        call helical_chart(chart, 0.15_dp)
        call build_fb(fb, chart)
        call synthetic_basis(basis, 2048)
        t_ref = 0
        do i = 1, size(counts)
            !$ call omp_set_num_threads(counts(i))
            t1 = huge(1.0_dp)
            do rep = 1, 3
                t0 = 0
                !$ t0 = omp_get_wtime()
                call ideal_fb_circuits_build(fb, basis, run, st, nphi_surface=384)
                if (st /= IDEAL_FB_OK) error stop 'scaling traces'
                !$ t1 = min(t1, omp_get_wtime() - t0)
            end do
            if (i == 1) then
                ref = run; t_ref = t1
            end if
            print '(a,i3,a,f9.4,a,f6.2,a,l1)', 'threads ', counts(i), '  time ', t1, &
                ' s  speedup ', t_ref/t1, '  bitwise ', &
                all(run%drive == ref%drive) .and. all(run%load == ref%load) &
                .and. all(run%f_vac == ref%f_vac)
        end do
    end subroutine
    subroutine curve_derivatives
        type(ideal_fb_t) :: fb
        type(curve_chart_t) :: chart
        type(coil_circuit_basis_t) :: basis
        type(ideal_fb_circuits_t) :: circuits
        integer, allocatable :: rowptr(:), colidx(:)
        real(dp), allocatable :: val(:), g(:), gp(:), gm(:), v(:), hv(:), du(:)
        real(dp) :: lp, lm, l0, h, fd, an, err, scale, worst(4)
        integer :: ntor, i, st, p, b0, b(5)
        call helical_chart(chart, 0.15_dp)
        call synthetic_basis(basis, 64)
        do ntor = 1, 2
            call build_fb(fb, chart, 8, ntor)
            call ideal_fb_circuits_build(fb, basis, circuits, st, background=filament, &
                nphi_surface=nphi)
            if (st /= IDEAL_FB_OK) error stop 'curve traces'
            call ideal_fb_set_currents(fb, circuits, p0, st)
            if (st /= IDEAL_FB_OK) error stop 'curve currents'
            allocate (g(fb%n), gp(fb%n), gm(fb%n), v(fb%n), hv(fb%n), du(fb%n))
            do i = 1, fb%n
                du(i) = 1.0e-3_dp*sin(1.3_dp*i)
            end do
            du(fb%np + 2*fb%nx + 1:) = 1.0e2_dp*du(fb%np + 2*fb%nx + 1:)
            du(fb%np + fb%nx + 1:fb%np + 2*fb%nx) = 1.0e5_dp &
                *du(fb%np + fb%nx + 1:fb%np + 2*fb%nx)
            call ideal_fb_step(fb, du, 1.0_dp)
            call attach_frozen_exterior(fb)
            do i = 1, fb%n
                v(i) = cos(0.71_dp*i)
            end do
            v(fb%np + 2*fb%nx + 1:) = 1.0e2_dp*v(fb%np + 2*fb%nx + 1:)
            call ideal_fb_gradient(fb, g, st, l0)
            if (st /= IDEAL_FB_OK) error stop 'curve base state'
            h = 1.0e-6_dp
            call ideal_fb_step(fb, v, h)
            call ideal_fb_gradient(fb, gp, st, lp)
            call ideal_fb_step(fb, v, -2*h)
            call ideal_fb_gradient(fb, gm, st, lm)
            call ideal_fb_step(fb, v, h)
            fd = (lp - lm)/(2*h); an = dot_product(g, v)
            err = abs(fd - an)/max(abs(an), 1.0e-30_dp)
            print '(a,i2,a,2es16.8,a,es9.2)', 'curve N =', ntor, ' dL.v fd/an', fd, an, &
                ' rel', err
            if (err > 1.0e-6_dp) error stop 'curve-chart gradient mismatch'
            call ideal_fb_hessian(fb, rowptr, colidx, val)
            hv = 0
            do i = 1, fb%n
                do p = rowptr(i), rowptr(i + 1) - 1
                    hv(i) = hv(i) + val(p)*v(colidx(p))
                end do
            end do
            gp = (gp - gm)/(2*h)
            scale = maxval(abs(hv))
            b = [0, fb%np, fb%np + fb%nx, fb%np + 2*fb%nx, fb%n]
            do b0 = 1, 4
                worst(b0) = maxval(abs(gp(b(b0) + 1:b(b0 + 1)) - hv(b(b0) + 1:b(b0 + 1)))) &
                    /max(maxval(abs(hv(b(b0) + 1:b(b0 + 1)))), tiny(1.0_dp))
            end do
            err = maxval(abs(gp - hv))/scale
            print '(a,i2,a,es9.2,a,4es9.2)', 'curve N =', ntor, ' |Hv - fd|/|Hv|', err, &
                '  by block (plasma, shell map, multiplier, potential)', worst
            if (err > 1.0e-5_dp) error stop 'curve-chart Hessian mismatch'
            deallocate (g, gp, gm, v, hv, du)
        end do
        print '(a)', 'curve-chart derivatives ok'
    end subroutine

    subroutine exterior(ncontour, nphi_s, gate)
        type(ideal_fb_t) :: fbd, fbp, fbt
        type(curve_chart_t) :: chart
        type(bem_pair_options_t) :: options
        type(bem_pair_report_t) :: rd, rp
        integer, parameter :: counts(5) = [1, 2, 4, 8, 16]
        integer, intent(in) :: ncontour, nphi_s
        logical, intent(in) :: gate
        real(dp), allocatable :: s1(:, :)
        real(dp) :: t0, t_dense, t_per, t_ref, err, scale, emin
        integer :: st, i
        call helical_chart(chart, 0.15_dp)
        call build_fb(fbd, chart, ncontour, 1)
        call build_fb(fbp, chart, ncontour, 1)
        print '(a,i0,a,i0,a,i0,a)', 'Gamma_c: ', fbd%nedge, ' P2 edges x ', nphi_s, &
            ' slices, ', 2*fbd%nedge*nphi_s, ' panels, nfp = 3'
        t0 = wall()
        call ideal_fb_use_bem3d(fbp, nphi_s, options, rp, st, periodic=.true.)
        t_per = wall() - t0
        if (st /= IDEAL_FB_OK) error stop 'periodic R2 build'
        print '(a,f9.3,a,i0,a,i0,a,es9.2)', 'periodic build ', t_per, ' s, pairs ', &
            rp%accepted_pairs, ', max degree ', rp%maximum_degree_used, &
            ', error ratio ', rp%error_ratio
        call exact_full_surface(fbp, nphi_s, options, s1)
        scale = maxval(abs(fbp%exterior_3d%s))
        err = maxval(abs(s1 - fbp%exterior_3d%s))/scale
        print '(a,es10.2)', 'periodic vs full-surface DtN of the same operators ', err
        if (gate .and. err > 1e-11_dp) error stop 'zero-harmonic reduction'
        t0 = wall()
        call ideal_fb_use_bem3d(fbd, nphi_s, options, rd, st)
        t_dense = wall() - t0
        if (st /= IDEAL_FB_OK) error stop 'dense R2 build'
        print '(a,f9.3,a,i0,a,i0)', 'dense build    ', t_dense, ' s, pairs ', &
            rd%accepted_pairs, ', max degree ', rd%maximum_degree_used
        err = maxval(abs(fbd%exterior_3d%s - fbp%exterior_3d%s))/scale
        print '(a,es10.2,a,f6.2)', 'periodic vs dense all-pairs build ', err, &
            '  speedup ', t_dense/t_per
        if (gate .and. err > 10*options%relative_tolerance) &
            error stop 'periodic and dense R2 differ beyond the pair tolerance'
        err = maxval(abs(fbp%exterior_3d%s - transpose(fbp%exterior_3d%s)))/scale
        emin = smallest_eigenvalue(fbp%exterior_3d%s)/scale
        print '(a,es10.2,a,es10.2)', 'asymmetry ', err, '  min eigenvalue / max|S| ', &
            emin
        if (gate .and. (err > 1e-13_dp .or. emin < -1e-10_dp)) &
            error stop 'R2 exterior not symmetric positive semidefinite'
        t_ref = 0
        do i = 1, size(counts)
            if (gate .and. counts(i) /= 2) cycle
            !$ call omp_set_num_threads(counts(i))
            call build_fb(fbt, chart, ncontour, 1)
            t0 = wall()
            call ideal_fb_use_bem3d(fbt, nphi_s, options, rp, st, periodic=.true.)
            t_per = wall() - t0
            if (st /= IDEAL_FB_OK) error stop 'periodic R2 build (threads)'
            if (t_ref == 0) t_ref = t_per
            err = maxval(abs(fbt%exterior_3d%s - fbp%exterior_3d%s))/scale
            print '(a,i3,a,f9.3,a,f6.2,a,es9.2)', 'threads ', counts(i), '  periodic ', &
                t_per, ' s  speedup ', t_ref/t_per, '  max rel diff ', err
            if (gate .and. err > 1e-13_dp) error stop 'thread count changed R2'
        end do
        print '(a)', 'periodic R2 exterior ok'
    end subroutine

    subroutine exact_full_surface(fb, nphi_s, options, s)
        ! R^T S R with S the full-surface DtN of the expanded circulant
        ! operators: the exact-algebra oracle of the zero-harmonic reduction.
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: nphi_s
        type(bem_pair_options_t), intent(in) :: options
        real(dp), allocatable, intent(out) :: s(:, :)
        type(coupling_surface_t) :: surface
        type(bem_circulant_t) :: cop
        type(bem_operators_t) :: op
        type(bem_exterior_dtn_t) :: dtn
        type(bem_pair_report_t) :: report
        integer, allocatable :: faces(:, :), td(:, :), fd(:, :)
        real(dp), allocatable :: restriction(:, :), full(:, :)
        integer :: st
        call ideal_fb_surface(fb, nphi_s, surface, faces, restriction, st)
        if (st /= 0) error stop 'surface'
        call assemble_parametric_bem_p2_circulant(surface, faces, fb%nfp, options, &
            cop, td, fd, report, st)
        if (st /= 0) error stop 'circulant assembly'
        call cop%expand(op, st)
        if (st /= 0) error stop 'expand'
        call build_exterior_dtn(op, dtn, st)
        if (st /= 0) error stop 'full DtN'
        call dtn%matrix(full, st)
        s = matmul(transpose(restriction), matmul(full, restriction))
    end subroutine

    real(dp) function smallest_eigenvalue(a) result(lambda)
        real(dp), intent(in) :: a(:, :)
        real(dp), allocatable :: b(:, :), w(:), work(:)
        integer :: n, info, lwork
        external :: dsyev
        n = size(a, 1); b = 0.5_dp*(a + transpose(a))
        lwork = max(1, 3*n)
        allocate (w(n), work(lwork))
        call dsyev('N', 'U', n, b, n, w, work, lwork, info)
        if (info /= 0) error stop 'dsyev'
        lambda = w(1)
    end function

    real(dp) function wall() result(t)
        integer(8) :: c, r
        call system_clock(c, r)
        t = real(c, dp)/real(r, dp)
    end function
end program test_ideal_fb_traces3d
