module ideal_fb_test_field
    ! Smooth nested frozen data for the derivative checks: elliptic labels
    ! about (3.05, 0), F = 6 T m, pressure linear in the label.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_ideal_map, only: ideal_frozen_t
    implicit none
    private
    type, extends(ideal_frozen_t), public :: screw_field_t
        real(dp) :: major = 20, pitch = 0.5_dp
    contains
        procedure :: evaluate => screw_evaluate
    end type
    type, extends(ideal_frozen_t), public :: simple_field_t
        logical :: solovev = .false.
    contains
        procedure :: evaluate => simple_evaluate
    end type
contains
    subroutine screw_evaluate(field, point, chi, chi_gradient, b3, pressure)
        class(screw_field_t), intent(in) :: field
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: chi, chi_gradient(2), b3, pressure
        real(dp) :: x, y
        x = point(1) - field%major; y = point(2)
        ! Unit-radius uniform-current screw pinch, mu0=1. The toroidal
        ! embedding approaches the exact cylinder as major -> infinity.
        chi = field%major*(x*x + y*y - 1)/2
        chi_gradient = field%major*[x, y]
        b3 = 10*field%pitch*field%major/point(1)
        pressure = 1 - x*x - y*y
    end subroutine

    subroutine simple_evaluate(field, point, chi, chi_gradient, b3, pressure)
        class(simple_field_t), intent(in) :: field
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: chi, chi_gradient(2), b3, pressure
        if (field%solovev) then
            ! Solov'ev polynomial: Delta* chi = 16 C R^2, F constant.
            ! Its axis is deliberately displaced by +37 mm in Z.
            chi = 0.04_dp*((point(1)**2 - 3.05_dp**2)**2 &
                           + 4*point(1)**2*(point(2) - 0.037_dp)**2 - 10.0_dp)
            chi_gradient = 0.04_dp*[4*point(1)*(point(1)**2 - 3.05_dp**2 &
                                                + 2*(point(2) - 0.037_dp)**2), &
                                    8*point(1)**2*(point(2) - 0.037_dp)]
            b3 = 6.0_dp/point(1)
            pressure = -16*0.04_dp*min(chi, 0.0_dp)/(4.0e-7_dp*acos(-1.0_dp))
            return
        end if
        chi = 0.5_dp*((point(1) - 3.05_dp)**2 + (point(2)/1.5_dp)**2)
        chi_gradient = [point(1) - 3.05_dp, point(2)/1.5_dp**2]
        b3 = 6.0_dp/point(1)
        pressure = 1.0e4_dp*max(0.0_dp, 1 - chi/0.18_dp)
    end subroutine
end module ideal_fb_test_field

program test_ideal_free_boundary
    ! N4 building blocks of the 3D free boundary (kin6d_ideal_free_boundary).
    !   derivatives      gradient of L against central differences of L, KKT
    !                    Hessian against differences of the gradient (N = 0
    !                    and N = 1, nfp = 3, perturbed state, all blocks)
    !   vacuum <k>       exact n = 0 vacuum oracle: PF filaments outside
    !                    Gamma_c plus a current filament inside Gamma_p, Gamma_p
    !                    a flux surface of that field. The shell solution must
    !                    be the filament + coil field itself (zero exterior
    !                    reaction); L2 error over the shell on three meshes
    !                    (48, 96, 192 Gamma_c edges), observed order >= k - 0.6
    !   vacuum3d <k>     n = 3 vacuum drive from a 3-fold window-coil set:
    !                    the shell field must not depend on the choice of
    !                    Gamma_c (two coupling curves, common probe points),
    !                    and the n /= 0 coil traces must be curl free
    !   slide            n = 0 relabelling gauge: the slide moments of every
    !                    curl-free displacement (translations, elongation,
    !                    radial stretch) vanish to discretisation accuracy and
    !                    converge under refinement, while a relabelling slide
    !                    has O(1) moments (the gauge fixes labels, not shape)
    !   vertical         slide-gauged Newton matrix against differences of the
    !                    gradient (N = 0, 1), and quadratic return of the gauged
    !                    Newton step to a manufactured displaced Solov'ev state
    !   hypotheses       ideal_fb_hypotheses: (V), (R) and In(S_red); the
    !                    FB-D2 bordering relation against a dense eigen-solve of
    !                    the slide-bordered KKT matrix (N = 0, 1)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_ideal_map, only: ideal_frozen_t
    use kin6d_equilibrium_boundary_3d, only: gamma_c_3d_t
    use kin6d_ideal_free_boundary, only: ideal_fb_t, ideal_fb_init, &
        ideal_fb_set_drive, ideal_fb_gradient, ideal_fb_hessian, ideal_fb_step, &
        ideal_fb_solve_vacuum, ideal_fb_shell_field, IDEAL_FB_OK, &
        ideal_fb_slide_system, ideal_fb_hypotheses_t, ideal_fb_hypotheses, &
        ideal_fb_hypotheses_print, ideal_fb_mode_inertia, ideal_fb_mode_lanczos, &
        ideal_fb_newton_step, ideal_fb_solve_shell_map, ideal_fb_solve_multipliers, &
        ideal_fb_solve, ideal_fb_report_t, ideal_fb_mode_energy, &
        ideal_fb_mode_stability, ideal_fb_stability_t, IDEAL_FB_BAD_INPUT, &
        IDEAL_FB_NOT_CONVERGED, ideal_fb_attach_exterior_3d, ideal_fb_mode_weight
    use kin6d_curved_triangle, only: triangle_locate
    use kin6d_coils_3d, only: coil_set_t, coil_set_from_polylines, coil_set_modes
    use fortfem_axisymmetric_gs_kernel, only: gs_green_gradient, gs_green
    use fortnum_symmetric_eigen, only: dense_oracle => symmetric_eigen
    use fortnum_quadrature, only: gauss_legendre
    use fortfem_gauss_quadrature_2d, only: gauss_quadrature_triangle_t
    use kin6d_triangle_basis, only: triangle_basis_jet
    use ideal_fb_kink_leaves, only: cylinder_kink_energy, &
        cylinder_kink_energy_jvp, cylinder_kink_energy_vjp
    use kin6d_s0_spectrum, only: s0_pencil_t, s0_report_t, s0_prepare, &
        s0_solve, s0_destroy
    use fortnum_eigen, only: eigen_options_t, EIGEN_OK, EIGEN_NEAREST
    use fortsparse, only: csc_z_t, csc_from_triplet, fortsparse_status_t, FORTSPARSE_OK
    use ideal_fb_test_field, only: simple_field_t, screw_field_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 4.0e-7_dp*pi
    ! Gamma_c ellipse, PF filaments (R, Z, A) and the inner current filament.
    real(dp), parameter :: rc = 3.0_dp, zc = 0.0_dp, ac = 1.5_dp, bc = 2.0_dp
    real(dp), parameter :: pf(3, 3) = reshape([5.2_dp, 2.0_dp, 4.0e5_dp, &
        5.2_dp, -2.0_dp, 4.0e5_dp, 1.2_dp, 0.0_dp, -3.0e5_dp], [3, 3])
    real(dp), parameter :: fil(2) = [3.05_dp, 0.05_dp], f_vac = 6.0_dp
    real(dp) :: circulation = -mu0*1.0e6_dp
    real(dp) :: gamma_scale = 1, area_factor = 0.5_dp, window_current = 1.0e5_dp
    type(coil_set_t) :: windows
    logical :: windows_on = .false.
    character(16) :: mode, argument
    integer :: degree
    real(dp) :: kink_major = 20
    mode = 'derivatives'; degree = 2
    if (command_argument_count() > 0) call get_command_argument(1, mode)
    if (command_argument_count() > 1) then
        call get_command_argument(2, argument)
        read (argument, *) degree
    end if
    select case (trim(mode))
    case ('derivatives')
        call derivative_checks()
    case ('frozen-exterior')
        call derivative_checks(.true.)
    case ('vertical')
        call vertical_checks()
    case ('vacuum')
        call vacuum_ladder()
    case ('vacuum3d')
        call vacuum_3d()
    case ('slide')
        call slide_checks()
    case ('kink')
        call kink_checks()
    case ('spectrum')
        call spectrum_checks()
    case ('kink_growth')
        call kink_checks(.true.)
    case ('hypotheses')
        call hypothesis_checks()
    case default
        error stop 'unknown mode'
    end select
contains
    subroutine spectrum_checks()
        ! Independent dense Schur solve/eigenspectrum, without the production
        ! mode extraction or sparse inertia. The full Hessian is the shared
        ! physical definition; indexing and elimination here are an oracle.
        type(ideal_fb_t) :: fb
        type(ideal_fb_stability_t) :: report
        integer, allocatable :: rp(:), ci(:), keep(:), piv(:)
        real(dp), allocatable :: val(:), h(:, :), c(:, :), b(:, :), s(:, :)
        real(dp), allocatable :: mass(:), w(:), vectors(:, :)
        real(dp) :: shifts(3), lambda, allowance, norm_s
        integer :: i, j, p, phase, n, np, nr, st, count, k, expected(3)
        integer :: counts(3, 3)
        interface
            subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
                import dp
                integer, intent(in) :: n, nrhs, lda, ldb
                integer, intent(out) :: ipiv(*), info
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            end subroutine
        end interface
        call build(8, 6, 2, .false., 3, 1, fb, solovev=.true.)
        call ideal_fb_set_drive(fb, sampler, st)
        if (st /= IDEAL_FB_OK) error stop 'spectrum drive'
        call ideal_fb_solve_vacuum(fb, st)
        if (st /= IDEAL_FB_OK) error stop 'spectrum vacuum'
        call ideal_fb_solve_multipliers(fb, st)
        if (st /= IDEAL_FB_OK) error stop 'spectrum multipliers'
        call ideal_fb_hessian(fb, rp, ci, val)
        allocate (keep(fb%n)); keep = 0; count = 0
        do i = 1, fb%n
            if (i <= fb%np) then
                phase = mod((i - 1)/2, fb%nmode) + 1
            else if (i <= fb%np + 2*fb%nx) then
                phase = mod((i - fb%np - 1)/3, fb%nmode) + 1
            else
                phase = mod(i - fb%np - 2*fb%nx - 1, fb%nmode) + 1
            end if
            if (phase == 1) cycle
            count = count + 1; keep(i) = count
        end do
        n = count; np = 4*fb%plasma%nodes; nr = n - np
        allocate (h(n, n), c(nr, nr), b(nr, np), piv(nr), s(np, np), &
            mass(np), w(np), vectors(np, np))
        h = 0
        do i = 1, fb%n
            if (keep(i) == 0) cycle
            do p = rp(i), rp(i + 1) - 1
                j = ci(p)
                if (keep(j) /= 0) h(keep(i), keep(j)) = val(p)
            end do
        end do
        h = 0.5_dp*(h + transpose(h))
        do i = 1, np
            mass(i) = abs(h(i, i))
        end do
        c = h(np + 1:, np + 1:); b = h(np + 1:, :np)
        call dgesv(nr, np, c, nr, piv, b, nr, st)
        if (st /= 0) error stop 'dense eliminated block singular'
        s = h(:np, :np) - matmul(h(:np, np + 1:), b)
        do j = 1, np
            do i = 1, np
                s(i, j) = s(i, j)/sqrt(mass(i)*mass(j))
            end do
        end do
        norm_s = maxval(sum(abs(s), dim=2))
        call dense_oracle(s, w, vectors, st)
        if (st /= 0) error stop 'dense spectrum failed'
        shifts = [0.0_dp, 1.e-4_dp, 1.e-2_dp]
        call ideal_fb_mode_inertia(fb, 2, shifts, counts, st)
        if (st /= IDEAL_FB_OK) error stop 'sparse inertia unavailable'
        do k = 1, 3
            expected = [count_negative(w - shifts(k)), 0, count_positive(w - shifts(k))]
            print '(a,es12.4,6i7)', 'inertia dense/sparse ', shifts(k), &
                expected, counts(:, k)
            if (any(counts(:, k) /= expected)) error stop 'Schur inertia mismatch'
        end do
        call ideal_fb_mode_lanczos(fb, 2, lambda, report, st)
        print '(a,i0,4es20.10)', 'Lanczos status/value/radius/backward/dense ', &
            st, lambda, report%residual_radius, report%backward_error, w(minloc(abs(w), 1))
        if (st /= IDEAL_FB_OK) error stop 'Schur Lanczos failed'
        allowance = report%residual_radius + 50*np*epsilon(1.0_dp)*norm_s
        if (minval(abs(w - lambda)) > allowance) error stop 'eigenpair outside dense oracle'
        if (report%backward_error > 1.e-10_dp) error stop 'backward gate failed'
        fb%plasma%x(1, 2, 1) = 1.e-3_dp
        call ideal_fb_mode_inertia(fb, 2, shifts, counts, st)
        if (st /= IDEAL_FB_BAD_INPUT) error stop 'nonaxisymmetric inertia admitted'
    end subroutine

    integer function count_negative(x) result(n)
        real(dp), intent(in) :: x(:)
        n = count(x < 0)
    end function

    integer function count_positive(x) result(n)
        real(dp), intent(in) :: x(:)
        n = count(x > 0)
    end function

    subroutine kink_curve(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [kink_major + 2*cos(2*pi*t), 2*sin(2*pi*t)]
    end subroutine

    subroutine kink_sampler(r, z, k, bc, bs)
        ! Curl-free toroidal ring field of unit circulation 2*pi. Its
        ! large-aspect-ratio limit is B_theta=1/r. No conducting wall.
        real(dp), intent(in) :: r, z
        integer, intent(in) :: k
        real(dp), intent(out) :: bc(3), bs(3)
        real(dp) :: g, gg(4)
        bc = 0; bs = 0
        if (k /= 0) return
        call gs_green_gradient(r, z, kink_major, 0.0_dp, g, gg)
        bc([1, 3]) = 2*pi*[gg(2), -gg(1)]/r
    end subroutine

    subroutine kink_checks(growth)
        logical, intent(in), optional :: growth
        ! Independent cylinder oracle: m=1 rigid helical displacement, ka=0.1.
        ! delta W / (pi*a^2*L/(2*mu0)) = q^2-1 + f(ka)*(q-1)^2,
        ! f=K1/(K1+ka*K0). Its ka -> 0 limit is 2*q*(q-1).
        ! The two real phases are necessary; keeping only cosine replaces
        ! the helical drive by a standing perturbation.
        type(ideal_fb_t) :: fb
        type(gs_state_t) :: mesh
        type(cylindrical_chart_t) :: chart
        type(screw_field_t) :: field
        type(ideal_fb_stability_t) :: ep, em, single, probe
        real(dp), allocatable :: chain(:, :), contour(:, :), u(:)
        real(dp) :: q, scale, exact, err, previous(2), lambda, wp, wv, ws
        real(dp) :: jv(4), gq, gf, factor, k0, k1, t, gx(96), gw(96)
        real(dp) :: growth_errors(3), physical_lambda, physical_mass, density, r
        integer :: last_level, alevel, cell, point
        logical :: rates
        integer :: level, iq, j, a, st
        rates = .false.
        if (present(growth)) rates = growth
        last_level = 2
        if (rates) last_level = 3
        density = 1.3_dp
        call gauss_legendre(96, gx, gw)
        k0 = 0; k1 = 0
        do j = 1, 96
            t = 6*(gx(j) + 1)
            k0 = k0 + 6*gw(j)*exp(-0.1_dp*cosh(t))
            k1 = k1 + 6*gw(j)*exp(-0.1_dp*cosh(t))*cosh(t)
        end do
        factor = k1/(k1 + 0.1_dp*k0)
        previous = huge(1.0_dp)
        do level = 1, last_level
            kink_major = 40.0_dp*2**level
            field%major = kink_major
            allocate (contour(2, 16), chain(2, 16))
            do j = 1, 16
                call kink_curve(real(j - 1, dp)/16, contour(:, j))
                chain(:, j) = [kink_major + cos(2*pi*(j - 1)/16), &
                    sin(2*pi*(j - 1)/16)]
            end do
            do iq = 1, 2
                q = merge(0.5_dp, 1.5_dp, iq == 1)
                field%pitch = q
                call initialize_gs(mesh, contour, 0.12_dp, st, 3, &
                    curve=kink_curve, interior_chain=chain)
                if (st /= 0) error stop 'kink mesh'
                call ideal_fb_init(fb, mesh, chain, chart, field, &
                    nint(kink_major/10), 1, 10*q*kink_major, 2*pi, st, mu0=1.0_dp)
                if (st /= IDEAL_FB_OK) error stop 'kink init'
                call ideal_fb_set_drive(fb, kink_sampler, st)
                if (st /= IDEAL_FB_OK) error stop 'kink drive'
                call ideal_fb_solve_vacuum(fb, st)
                if (st /= IDEAL_FB_OK) error stop 'kink vacuum'
                allocate (u(4*fb%plasma%nodes))
                u = 0
                do a = 1, fb%plasma%nodes
                    u(4*a - 3) = 1; u(4*a) = 1
                end do
                call ideal_fb_mode_energy(fb, 2, u, ep, st)
                if (st /= IDEAL_FB_OK) then
                    print *, 'energy status/amplitudes', st, &
                        maxval(abs(fb%plasma%x(:, 2:, :))), &
                        maxval(abs(fb%xs(:, 2:, :))), maxval(abs(fb%phi(2:, :)))
                    error stop 'kink energy plus'
                end if
                u(4::4) = -1
                call ideal_fb_mode_energy(fb, 3, u, em, st)
                if (st /= IDEAL_FB_OK) error stop 'kink energy minus'
                if (ep%total < em%total) em = ep
                scale = 10*pi*pi
                call cylinder_kink_energy(q, factor, wp, wv, exact, ws)
                call cylinder_kink_energy_jvp(q, factor, 0.7_dp, 0.0_dp, &
                    jv(1), jv(2), jv(3), jv(4))
                call cylinder_kink_energy_vjp(q, factor, 1.0_dp, -2.0_dp, &
                    3.0_dp, -4.0_dp, gq, gf)
                if (abs(dot_product(jv, [1.0_dp, -2.0_dp, 3.0_dp, -4.0_dp]) &
                    - 0.7_dp*gq) > 1.0e-12_dp) error stop 'kink JVP/VJP'
                err = abs(em%total/scale - exact)
                print '(a,f6.1,a,f4.1,a,4es13.5)', 'kink A=', kink_major, &
                    ' q=', q, ' plasma,vacuum,total,limit ', em%plasma/scale, &
                    em%vacuum/scale, em%total/scale, exact
                if (em%total*exact <= 0) error stop 'external kink threshold sign'
                if (abs(em%constraint)/scale > 1.0e-8_dp) &
                    error stop 'kink constraint work'
                if (err > 0.03_dp*abs(exact)) error stop 'cylinder energy limit'
                if (err >= previous(iq)) error stop 'cylinder aspect convergence'
                previous(iq) = err
                if (abs(em%plasma/scale - wp) > 0.03_dp*abs(wp)) &
                    error stop 'cylinder plasma energy'
                if (abs(em%vacuum/scale - wv) > 0.03_dp*abs(wv)) &
                    error stop 'cylinder vacuum energy'
                if (rates .and. iq == 1) then
                    ! Rigid circular polarization has |xi_real|^2=1 at
                    ! every toroidal angle. Integrate physical rho*dV on
                    ! the actual curved plasma mesh, not coefficient metric.
                    physical_mass = 0
                    do cell = 1, fb%plasma%cells
                        do point = 1, fb%plasma%nq
                            r = dot_product(fb%plasma%reference(1, &
                                fb%plasma%conn(:,cell)),fb%plasma%shape(:,point))
                            physical_mass = physical_mass+density*r &
                                *fb%plasma%weight(point)*fb%plasma%det(point,cell) &
                                *2*pi/fb%plasma%nfp
                        end do
                    end do
                    ! Two real phases are an exact Fourier degeneracy.
                    ! This is the projected rigid subspace, not the minimum
                    ! over the full displacement space.
                    call rigid_spectrum(em%total,physical_mass,physical_lambda)
                    alevel = level
                    growth_errors(alevel) = abs(physical_lambda-exact/density) &
                        /abs(exact/density)
                    print '(a,f7.1,a,4es18.9)', '[V] rigid kink R=',kink_major, &
                        ' gamma,S0-lambda,oracle-lambda,relative-error=', &
                        sqrt(-physical_lambda),physical_lambda,exact/density, &
                        growth_errors(alevel)
                    if (growth_errors(alevel)>0.03_dp) error stop 'physical kink rate'
                end if
                if (iq == 1 .and. .not. rates) then
                    u(4::4) = 0
                    call ideal_fb_mode_energy(fb, 2, u, single, st)
                    if (st /= IDEAL_FB_OK) error stop 'standing energy'
                    print '(a,es13.5)', 'cosine-only energy/scale ', single%total/scale
                    if (single%total <= 0) error stop 'standing-wave control'
                    ! A shift near a negative Rayleigh witness must return a
                    ! negative eigenvalue from the full coupled harmonic.
                    call ideal_fb_mode_stability(fb, 2, lambda, st, &
                        shift=2*em%total/em%metric_norm, report=probe)
                    print '(a,es13.5,a,i0,a,es13.5)', 'paired eigenvalue ', lambda, &
                        ' status ', st, ' residual ', probe%relative_residual
                    if (st /= IDEAL_FB_OK) error stop 'kink eigen convergence'
                    if (lambda >= 0) error stop 'missing helical instability'
                    call ideal_fb_mode_stability(fb, 2, lambda, st, iterations=1, &
                        shift=-0.1_dp)
                    if (st /= IDEAL_FB_NOT_CONVERGED) &
                        error stop 'unconverged probe accepted'
                    fb%plasma%x(1, 2, 1) = 1.0e-3_dp
                    call ideal_fb_mode_stability(fb, 2, lambda, st)
                    if (st /= IDEAL_FB_BAD_INPUT) &
                        error stop 'nonaxisymmetric probe accepted'
                end if
                deallocate (u)
            end do
            deallocate (chain, contour)
        end do
        if (rates) then
            print '(a,2f12.6)', '[V] rigid kink aspect-error orders=', &
                log(growth_errors(:2)/growth_errors(2:))/log(2.0_dp)
        end if
    end subroutine

    subroutine rigid_spectrum(energy,physical_mass,lambda)
        real(dp), intent(in) :: energy,physical_mass
        real(dp), intent(out) :: lambda
        type(csc_z_t) :: stiffness,mass
        type(fortsparse_status_t) :: status
        type(s0_pencil_t) :: pencil
        type(s0_report_t) :: report
        type(eigen_options_t) :: options
        complex(dp) :: kv(4),mv(4),values(1),vectors(2,1),start(2)
        real(dp) :: residual(1),allowance
        kv = 0; mv = 0
        kv([1,4]) = 2*energy; mv([1,4]) = physical_mass
        call csc_from_triplet(2,2,[1,2,1,2],[1,1,2,2],kv,stiffness,status)
        if(status%code/=FORTSPARSE_OK) error stop 'rigid stiffness'
        call csc_from_triplet(2,2,[1,2,1,2],[1,1,2,2],mv,mass,status)
        if(status%code/=FORTSPARSE_OK) error stop 'rigid physical mass'
        call s0_prepare(pencil,stiffness,mass, &
            cmplx(2.2_dp*energy/physical_mass,0.0_dp,dp),report)
        if(report%status/=EIGEN_OK) error stop 'rigid shift factor'
        options%max_basis=2; options%which=EIGEN_NEAREST
        start=[cmplx(1.0_dp,0.0_dp,dp),cmplx(0.0_dp,1.0_dp,dp)]
        call s0_solve(pencil,start,1.0e-10_dp,options,values,vectors,residual,report)
        if(report%status/=EIGEN_OK) error stop 'rigid S0 spectrum'
        lambda=real(values(1),dp)
        allowance=residual(1)/sqrt(physical_mass) &
            +128*epsilon(lambda)*abs(2*energy/physical_mass)
        if(abs(values(1)-2*energy/physical_mass)>allowance) &
            error stop 'rigid original-equation tolerance'
        call s0_destroy(pencil)
    end subroutine rigid_spectrum

    subroutine field_ext(r, z, b, psi)
        ! Poloidal field and flux of the PF filaments plus the unit-circulation
        ! filament scaled by c: B_unit = (G_Z, -G_R)/R, psi_unit = -G.
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: b(2), psi
        real(dp) :: g, gg(4)
        integer :: k
        b = 0; psi = 0
        do k = 1, 3
            call gs_green_gradient(r, z, pf(1, k), pf(2, k), g, gg)
            psi = psi + mu0*pf(3, k)*g
            b = b + mu0*pf(3, k)*[-gg(2), gg(1)]/r
        end do
        call gs_green_gradient(r, z, fil(1), fil(2), g, gg)
        psi = psi - circulation*g
        b = b + circulation*[gg(2), -gg(1)]/r
    end subroutine

    subroutine sampler(r, z, k, bcos, bsin)
        real(dp), intent(in) :: r, z
        integer, intent(in) :: k
        real(dp), intent(out) :: bcos(3), bsin(3)
        real(dp) :: b(2), psi
        complex(dp) :: bn(3, 1, 1)
        integer :: st
        bcos = 0; bsin = 0
        if (k == 0) then
            call field_ext(r, z, b, psi)
            bcos = [b(1), 0.0_dp, b(2)]
        end if
        if (.not. windows_on) return
        call coil_set_modes(windows, window_current*[1.0_dp, 1.0_dp, 1.0_dp], &
            reshape([r, z], [2, 1]), 64, [3*k], bn, st)
        if (st /= 0) error stop 'coil modes'
        if (k == 0) then
            bcos = bcos + real(bn(:, 1, 1))
        else
            bcos = 2*real(bn(:, 1, 1)); bsin = -2*aimag(bn(:, 1, 1))
        end if
    end subroutine

    subroutine gamma_c(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [rc + gamma_scale*ac*cos(2*pi*t), zc + gamma_scale*bc*sin(2*pi*t)]
    end subroutine

    subroutine solovev_chain(n, chain)
        ! Gamma_p on the zero-pressure level chi = 0 of the Solov'ev data,
        ! traced along rays from (3.05, 0.037): the frozen field is then
        ! tangent to Gamma_p, as the slide gauge assumes (b.n = 0, FB-C2).
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: chain(:, :)
        type(simple_field_t) :: field
        real(dp) :: lo, hi, mid, chi, grad(2), b3, p, t
        integer :: j, it
        field%solovev = .true.
        allocate (chain(2, n))
        do j = 1, n
            t = 2*pi*(j - 1)/n
            lo = 0.0_dp; hi = 1.0_dp
            do it = 1, 200
                mid = 0.5_dp*(lo + hi)
                call field%evaluate([3.05_dp + mid*cos(t), 0.037_dp + mid*sin(t)], &
                    chi, grad, b3, p)
                if (chi < 0) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            chain(:, j) = [3.05_dp, 0.037_dp] + mid*[cos(t), sin(t)]
        end do
    end subroutine

    subroutine flux_surface_chain(n, chain)
        ! Level psi_ext(fil + 0.55 e_R) traced along rays from the filament.
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: chain(:, :)
        real(dp) :: level, lo, hi, mid, b(2), v, t
        integer :: j, it
        allocate (chain(2, n))
        call field_ext(fil(1) + 0.55_dp, fil(2), b, level)
        do j = 1, n
            t = 2*pi*(j - 1)/n
            lo = 0.05_dp; hi = 1.2_dp
            do it = 1, 200
                mid = 0.5_dp*(lo + hi)
                call field_ext(fil(1) + mid*cos(t), fil(2) + mid*sin(t), b, v)
                ! -c G dominates: psi_ext falls away from the filament (c < 0).
                if ((v - level)*sign(1.0_dp, -circulation) > 0) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            chain(:, j) = fil + mid*[cos(t), sin(t)]
        end do
    end subroutine

    subroutine build(nc, nchain, order, flux_chain, nfp, ntor, fb, solovev)
        integer, intent(in) :: nc, nchain, order, nfp, ntor
        logical, intent(in) :: flux_chain
        type(ideal_fb_t), intent(out) :: fb
        logical, intent(in), optional :: solovev
        type(gs_state_t) :: mesh
        type(cylindrical_chart_t) :: chart
        type(simple_field_t) :: field
        real(dp), allocatable :: contour(:, :), chain(:, :)
        real(dp) :: h
        integer :: j, st
        allocate (contour(2, nc))
        do j = 1, nc
            call gamma_c(real(j - 1, dp)/nc, contour(:, j))
        end do
        if (present(solovev)) then
            if (solovev) call solovev_chain(nchain, chain)
        end if
        if (allocated(chain)) then
            continue
        else if (flux_chain) then
            call flux_surface_chain(nchain, chain)
        else
            allocate (chain(2, nchain))
            do j = 1, nchain
                chain(:, j) = [3.05_dp + 0.6_dp*cos(2*pi*(j - 1)/nchain), &
                    0.9_dp*sin(2*pi*(j - 1)/nchain)]
            end do
        end if
        if (present(solovev)) field%solovev = solovev
        h = 2*pi*ac/nc
        call initialize_gs(mesh, contour, area_factor*h*h, st, order, curve=gamma_c, &
            interior_chain=chain)
        if (st /= 0) error stop 'mesh'
        call ideal_fb_init(fb, mesh, chain, chart, field, nfp, ntor, f_vac, &
            circulation, st, mu0=mu0)
        if (st /= IDEAL_FB_OK) then
            print *, 'init status', st
            error stop 'fb init'
        end if
    end subroutine

    subroutine slide_checks()
        type(ideal_fb_t) :: fb
        real(dp), allocatable :: u(:, :)
        real(dp) :: ratio(2), slide, curl_free
        integer :: level, k, a
        do level = 1, 2
            call build(24*2**level, 16*2**level, 3, .false., 1, 0, fb)
            allocate (u(2, fb%plasma%nodes))
            ! Relabelling slide: unit tangent of the reference surfaces.
            do a = 1, fb%plasma%nodes
                u(:, a) = [-fb%plasma%direction(2, a), fb%plasma%direction(1, a)]
            end do
            slide = slide_moments(fb, u)
            curl_free = 0
            do k = 1, 4
                select case (k)
                case (1)
                    u(1, :) = 1; u(2, :) = 0
                case (2)
                    u(1, :) = 0; u(2, :) = 1
                case (3)
                    u(1, :) = 0; u(2, :) = fb%plasma%reference(2, :)
                case (4)
                    u(1, :) = fb%plasma%reference(1, :) - 3.05_dp; u(2, :) = 0
                end select
                curl_free = max(curl_free, slide_moments(fb, u))
            end do
            ratio(level) = curl_free/slide
            print '(a,i0,a,i0,a,es10.3,a,es10.3)', 'slide gauge level ', level, &
                ' layers ', size(fb%slide, 3), ' slide moment ', slide, &
                ' curl-free / slide ', ratio(level)
            deallocate (u)
        end do
        if (.not. ratio(2) < 1.0e-3_dp) error stop 'slide moments see a curl-free motion'
        if (.not. ratio(2) < 0.5_dp*ratio(1)) error stop 'slide moments do not converge'
        if (.not. slide > 0.1_dp) error stop 'slide moments miss the relabelling slide'
    end subroutine

    real(dp) function slide_moments(fb, v) result(largest)
        ! Largest |M_l(v)| of the n = 0 slide-gauge rows for a nodal field v.
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(in) :: v(:, :)
        largest = maxval(abs(sum(sum(fb%slide*spread(v, 3, size(fb%slide, 3)), &
            dim=1), dim=1)))
    end function

    subroutine attach_frozen_oracle(fb)
        !! Independent SPD rank-one form, including cross-mode coupling.
        !! Energy difference checks the full-turn/nfp convention separately.
        type(ideal_fb_t), intent(inout) :: fb
        type(gamma_c_3d_t) :: boundary
        real(dp), allocatable :: trace(:), drive(:), load(:), v(:), g(:), u(:)
        real(dp) :: before, after, old_energy, new_energy, allowance
        integer :: n, i, m, k, st
        n = fb%nt*fb%nmode
        allocate (boundary%s(n, n), v(n), load(n), g(fb%n))
        boundary%genus = 1
        v = [(sin(.37_dp*i), i=1, n)]/sqrt(real(n, dp))
        boundary%s = spread(v, 2, n)*spread(v, 1, n)
        do i = 1, n
            boundary%s(i, i) = boundary%s(i, i) + 1
        end do
        trace = reshape(transpose(fb%phi(:, fb%trace)), [n])
        drive = reshape(fb%drive, [n]); u = trace - drive
        old_energy = 0
        do m = 1, fb%nmode
            k = m/2
            old_energy = old_energy + ideal_fb_mode_weight(fb, m)/fb%mu0* &
                (.5_dp*dot_product(u((m - 1)*fb%nt + 1:m*fb%nt), &
                matmul(fb%s_trace(:, :, k), u((m - 1)*fb%nt + 1:m*fb%nt))) - &
                dot_product(fb%load(:, m), fb%phi(m, fb%trace)))
            load((m - 1)*fb%nt + 1:m*fb%nt) = &
                fb%nfp*ideal_fb_mode_weight(fb, m)*fb%load(:, m)
        end do
        new_energy = (.5_dp*(sum(u*u) + dot_product(v, u)**2) - &
                      dot_product(load, trace))/(fb%mu0*fb%nfp)
        call ideal_fb_gradient(fb, g, st, before)
        if (st /= 0) error stop 'frozen exterior baseline'
        call ideal_fb_attach_exterior_3d(fb, boundary, drive, load, st)
        if (st /= 0) error stop 'frozen exterior attachment'
        call ideal_fb_gradient(fb, g, st, after)
        if (st /= 0) error stop 'frozen exterior attached state'
        allowance = 1e-11_dp*max(1.0_dp, abs(before), abs(after), abs(old_energy), abs(new_energy))
        print '(a,es15.7)', 'frozen exterior once-only nfp energy defect ', &
            abs(after - before - old_energy + new_energy)/max(1.0_dp, abs(new_energy))
        if (abs(after - before - old_energy + new_energy) > allowance) &
            error stop 'frozen exterior full-turn energy normalization'
    end subroutine

    subroutine derivative_checks(frozen_exterior)
        logical, intent(in), optional :: frozen_exterior
        type(ideal_fb_t) :: fb
        integer, allocatable :: rowptr(:), colidx(:)
        real(dp), allocatable :: val(:), g(:), gp(:), gm(:), v(:), hv(:), du(:)
        real(dp) :: lp, lm, l0, h, fd, an, err, scale
        integer :: ntor, i, st, p
        do ntor = 0, 1
            windows_on = ntor > 0
            if (windows_on) call make_windows(1.0e5_dp)
            if (present(frozen_exterior)) then
                call build(8, 6, 2, .false., 3, ntor, fb)
            else
                call build(24, 16, 2, .false., 3, ntor, fb)
            end if
            call ideal_fb_set_drive(fb, sampler, st)
            allocate (g(fb%n), gp(fb%n), gm(fb%n), v(fb%n), hv(fb%n), du(fb%n))
            ! Perturbed base state: map, multipliers, potential.
            do i = 1, fb%n
                du(i) = 1.0e-3_dp*sin(1.3_dp*i)
            end do
            du(fb%np + 2*fb%nx + 1:) = 1.0e2_dp*du(fb%np + 2*fb%nx + 1:)
            du(fb%np + fb%nx + 1:fb%np + 2*fb%nx) = 1.0e5_dp &
                *du(fb%np + fb%nx + 1:fb%np + 2*fb%nx)
            call ideal_fb_step(fb, du, 1.0_dp)
            if (present(frozen_exterior)) call attach_frozen_oracle(fb)
            do i = 1, fb%n
                v(i) = cos(0.71_dp*i)
            end do
            v(fb%np + 2*fb%nx + 1:) = 1.0e2_dp*v(fb%np + 2*fb%nx + 1:)
            call ideal_fb_gradient(fb, g, st, l0)
            if (st /= IDEAL_FB_OK) error stop 'base state folded'
            h = 1.0e-6_dp
            call ideal_fb_step(fb, v, h)
            call ideal_fb_gradient(fb, gp, st, lp)
            call ideal_fb_step(fb, v, -2*h)
            call ideal_fb_gradient(fb, gm, st, lm)
            call ideal_fb_step(fb, v, h)
            fd = (lp - lm)/(2*h); an = dot_product(g, v)
            err = abs(fd - an)/max(abs(an), 1.0e-30_dp)
            print '(a,i2,a,2es16.8,a,es9.2)', 'N =', ntor, ' dL.v fd/an', fd, an, &
                ' rel', err
            if (err > 1.0e-6_dp) error stop 'gradient mismatch'
            call ideal_fb_hessian(fb, rowptr, colidx, val)
            hv = 0
            do i = 1, fb%n
                do p = rowptr(i), rowptr(i + 1) - 1
                    hv(i) = hv(i) + val(p)*v(colidx(p))
                end do
            end do
            gp = (gp - gm)/(2*h)
            scale = maxval(abs(hv))
            err = maxval(abs(gp - hv))/scale
            print '(a,i2,a,es9.2,a,es9.2)', 'N =', ntor, ' |Hv - fd|/|Hv|', err, &
                '  |Hv|', scale
            call block_report('plasma', gp(1:fb%np), hv(1:fb%np))
            call block_report('shell map', gp(fb%np + 1:fb%np + fb%nx), &
                hv(fb%np + 1:fb%np + fb%nx))
            call block_report('multiplier', gp(fb%np + fb%nx + 1:fb%np + 2*fb%nx), &
                hv(fb%np + fb%nx + 1:fb%np + 2*fb%nx))
            call block_report('potential', gp(fb%np + 2*fb%nx + 1:), &
                hv(fb%np + 2*fb%nx + 1:))
            if (err > 1.0e-5_dp) error stop 'Hessian mismatch'
            ! Symmetry of the assembled matrix.
            call check_symmetry(rowptr, colidx, val)
            deallocate (g, gp, gm, v, hv, du)
        end do
        print '(a)', 'derivatives ok'
    end subroutine

    subroutine vertical_checks()
        ! Manufactured discrete load at a displaced Solov'ev state, retained
        ! unchanged throughout each Newton solve. This tests the nonlinear
        ! equations of the slide-gauged system, not a claim of a
        ! coil-supported analytic free boundary.
        type(ideal_fb_t) :: fb
        type(ideal_fb_report_t) :: report
        type(ideal_fb_report_t) :: checked_report
        integer, allocatable :: rp(:), ci(:), srp(:), sci(:)
        real(dp), allocatable :: val(:), sv(:), g(:), gp(:), gm(:), v(:), du(:)
        real(dp), allocatable :: hv(:), hh(:), fd(:), x0(:, :, :), d(:, :), q(:, :, :)
        real(dp) :: h, error(0:3), residual(0:3), scale, fd_error, dz, fz(2)
        integer :: ntor, i, p, st, a, it, sign_z, nl, n, k, j
        do ntor = 0, 1
            call build(16, 12, 2, .false., 3, ntor, fb, solovev=.true.)
            call ideal_fb_set_drive(fb, sampler, st)
            if (st /= IDEAL_FB_OK) error stop 'vertical drive'
            n = fb%n; nl = size(fb%slide, 3)
            allocate (g(n), gp(n), gm(n), du(n), v(n + nl), hv(n + nl), hh(n), &
                      fd(n + nl))
            ! Nonaxisymmetric perturbations exercise the off-diagonal blocks.
            do i = 1, n
                du(i) = 2.0e-4_dp*sin(0.71_dp*i)
            end do
            call ideal_fb_step(fb, du, 1.0_dp)
            call ideal_fb_gradient(fb, g, st)
            call ideal_fb_hessian(fb, rp, ci, val)
            call ideal_fb_slide_system(fb, rp, ci, val, 1.0_dp, srp, sci, sv)
            call check_symmetry(srp, sci, sv)
            ! Bordered residual [g(U) + G^T mu; G (U - U_start)] along (du, nu).
            v = [(cos(0.31_dp*i), i=1, n + nl)]
            do a = 1, fb%ns
                i = fb%np + fb%nx + 2 + 3*fb%nmode*(a - 1)
                v(i) = 1.0e3_dp
            end do
            du = v(:n)
            h = 1.0e-6_dp
            call ideal_fb_step(fb, du, h)
            call ideal_fb_gradient(fb, gp, st)
            call ideal_fb_step(fb, du, -2*h)
            call ideal_fb_gradient(fb, gm, st)
            call ideal_fb_step(fb, du, h)
            fd = 0
            fd(:n) = (gp - gm)/(2*h)
            do k = 1, nl
                do a = 1, fb%plasma%nodes
                    do j = 1, 2
                        i = j + 2*fb%nmode*(a - 1)
                        fd(i) = fd(i) + fb%slide(j, a, k)*v(n + k)
                        fd(n + k) = fd(n + k) + fb%slide(j, a, k)*du(i)
                    end do
                end do
            end do
            hv = 0
            do i = 1, n + nl
                do p = srp(i), srp(i + 1) - 1
                    hv(i) = hv(i) + sv(p)*v(sci(p))
                end do
            end do
            scale = maxval(abs(hv))
            fd_error = maxval(abs(fd - hv))/scale
            print '(a,i0,a,es12.4)', 'gauged N=', ntor, ' Jacobian FD error ', fd_error
            if (fd_error > 1.0e-6_dp) error stop 'gauged Jacobian mismatch'
            ! Net vertical force row (n = 0 plasma Z rows) on its own scale.
            hh = 0
            do i = 1, n
                do p = rp(i), rp(i + 1) - 1
                    hh(i) = hh(i) + val(p)*du(ci(p))
                end do
            end do
            fz = 0
            do a = 1, fb%plasma%nodes
                i = 2 + 2*fb%nmode*(a - 1)
                fz = fz + [fd(i) - dot_product(v(n + 1:n + nl), fb%slide(2, a, :)), &
                    hh(i)]
            end do
            fd_error = abs(fz(1) - fz(2))/max(abs(fz(2)), 1.0_dp)
            print '(a,es12.4)', 'vertical force row FD error ', fd_error
            if (fd_error > 1.0e-5_dp) error stop 'vertical force Jacobian mismatch'
            ! The affine manufactured load makes this exact discrete target.
            call ideal_fb_gradient(fb, g, st)
            fb%shift = g
            call ideal_fb_solve(fb, 1.0e-10_dp, 0, report, st, &
                residual_tolerance=1.0e-12_dp)
            if (st /= IDEAL_FB_OK .or. report%iterations /= 0) &
                error stop 'zero-step exact target not recognized'
            if (report%residual /= 0 .or. report%step /= 0) &
                error stop 'zero-step exact target report'
            call ideal_fb_solve(fb, 1.0e-10_dp, 2, report, st, &
                residual_tolerance=1.0e-12_dp)
            if (st /= IDEAL_FB_OK .or. report%iterations /= 0) &
                error stop 'exact manufactured target not recognized'
            x0 = fb%plasma%x
            ! Orthonormal slide rows: the displaced start keeps the target in
            ! the gauge (zero slide moments of target - start).
            q = fb%slide
            do k = 1, nl
                do j = 1, k - 1
                    q(:, :, k) = q(:, :, k) - sum(q(:, :, j)*q(:, :, k))*q(:, :, j)
                end do
                q(:, :, k) = q(:, :, k)/norm2(q(:, :, k))
            end do
            do sign_z = -1, 1, 2
                ! 2 mm: from 4 mm the slide-gauged iterates pass 3e-10 at step 2,
                ! where the quadratic test below would ask for 1e-17 at step 3,
                ! below the 2e-13 roundoff floor of the state.
                dz = sign_z*0.002_dp
                allocate (d(2, fb%plasma%nodes))
                d(1, :) = 0; d(2, :) = dz
                do k = 1, nl
                    d = d - sum(q(:, :, k)*d)*q(:, :, k)
                end do
                fb%plasma%x(1:2, 1, :) = x0(1:2, 1, :) + d
                print '(a,es12.4,a,es12.4)', 'gauged start: mean dZ ', &
                    sum(d(2, :))/size(d, 2), ' max |d - dZ| ', &
                    maxval(abs(d - spread([0.0_dp, dz], 2, size(d, 2))))
                deallocate (d)
                ! Keep copies of the other target variables: the loop restores
                ! all state via an inverse accumulated step below.
                block
                    real(dp), allocatable :: xs0(:, :, :), phi0(:, :), lam0(:, :, :)
                    xs0 = fb%xs; phi0 = fb%phi; lam0 = fb%lam
                    call ideal_fb_solve_shell_map(fb, st)
                    if (st /= IDEAL_FB_OK) error stop 'displaced shell map'
                    call ideal_fb_solve_vacuum(fb, st)
                    if (st /= IDEAL_FB_OK) error stop 'displaced vacuum'
                    call ideal_fb_solve_multipliers(fb, st)
                    if (st /= IDEAL_FB_OK) error stop 'displaced multipliers'
                    block
                        real(dp), allocatable :: x1(:, :, :), xs1(:, :, :)
                        real(dp), allocatable :: phi1(:, :), lam1(:, :, :)
                        x1 = fb%plasma%x; xs1 = fb%xs
                        phi1 = fb%phi; lam1 = fb%lam
                        ! A zero-step check must evaluate a displaced state,
                        ! retaining failure and leaving the plasma map alone.
                        call ideal_fb_solve(fb, 1.0e-10_dp, 0, report, st, &
                            residual_tolerance=1.0e-20_dp)
                        if (st /= IDEAL_FB_NOT_CONVERGED) &
                            error stop 'zero-step displaced target status'
                        if (.not. ieee_is_finite(report%residual)) &
                            error stop 'zero-step displaced target not finite'
                        if (report%residual <= 1.0e-20_dp .or. &
                            report%residual >= 1.0_dp) &
                            error stop 'zero-step displaced target residual'
                        if (report%iterations /= 0 .or. report%step /= 0) &
                            error stop 'zero-step displaced target step'
                        if (any(fb%plasma%x /= x1)) &
                            error stop 'zero-step check moved the plasma'
                        ! Independently recheck the state after one capped
                        ! step: the solve must report that returned state.
                        call ideal_fb_solve(fb, 1.0e-10_dp, 1, report, st, &
                            residual_tolerance=1.0e-20_dp)
                        if (st /= IDEAL_FB_NOT_CONVERGED .or. &
                            report%iterations /= 1) &
                            error stop 'one-step cap status'
                        call ideal_fb_solve(fb, 1.0e-10_dp, 0, checked_report, st, &
                            residual_tolerance=1.0e-20_dp)
                        if (abs(checked_report%residual - report%residual) &
                            > 1.0e-10_dp*report%residual) &
                            error stop 'capped solve reported stale residual'
                        ! The same final state must be admitted if it meets
                        ! the requested absolute target on that last step.
                        fb%plasma%x = x1; fb%xs = xs1
                        fb%phi = phi1; fb%lam = lam1
                        call ideal_fb_solve(fb, 1.0e-10_dp, 1, checked_report, st, &
                            residual_tolerance=2*report%residual)
                        if (st /= IDEAL_FB_OK .or. checked_report%iterations /= 1) &
                            error stop 'last-step convergence not recognized'
                        fb%plasma%x = x1; fb%xs = xs1
                        fb%phi = phi1; fb%lam = lam1
                        ! A loose relative/step tolerance must not bypass an
                        ! explicit residual target. The manufactured load is
                        ! an independent state oracle for the damped solve.
                        fb%tangential_damping = 1.0e-3_dp
                        call ideal_fb_solve(fb, 1.0_dp, 20, report, st, &
                                    step_tolerance=1.0_dp, residual_tolerance=1.0e-8_dp)
                        if (st /= IDEAL_FB_OK .or. report%residual > 1.0e-8_dp) &
                            error stop 'explicit vertical residual target'
                        if (maxval(abs(fb%plasma%x - x0)) > 1.0e-8_dp) &
                            error stop 'damped manufactured target error'
                        fb%plasma%x = x1; fb%xs = xs1
                        fb%phi = phi1; fb%lam = lam1
                        call ideal_fb_solve(fb, 1.0_dp, 20, report, st, &
                                   step_tolerance=1.0_dp, residual_tolerance=1.0e-20_dp)
                        if (st == IDEAL_FB_OK) &
                            error stop 'unmet residual target reported success'
                        fb%tangential_damping = 0
                        fb%plasma%x = x1; fb%xs = xs1
                        fb%phi = phi1; fb%lam = lam1
                    end block
                    do it = 0, 3
                        error(it) = max(maxval(abs(fb%plasma%x - x0)), &
                            maxval(abs(fb%xs - xs0)), &
                            maxval(abs(fb%phi - phi0))/f_vac, &
                            maxval(abs(fb%lam - lam0))/1.0e6_dp)
                        call ideal_fb_gradient(fb, g, st)
                        residual(it) = maxval(abs(g))
                        print '(a,3i3,2es14.5)', 'vertical Newton ', ntor, sign_z, it, &
                            error(it), residual(it)
                        if (it == 3) exit
                        call ideal_fb_newton_step(fb, st)
                        if (st /= IDEAL_FB_OK) error stop 'vertical Newton solve'
                    end do
                    if (error(1) > 0.1_dp*error(0)) error stop 'vertical return'
                    do it = 1, 2
                        if (error(it) > 1.0e-10_dp) then
                            if (error(it + 1) > 100*error(it)**2) &
                                error stop 'vertical Newton not quadratic'
                        end if
                    end do
                    if (error(3) > 1.0e-9_dp) error stop 'vertical target error'
                    fb%xs = xs0; fb%phi = phi0; fb%lam = lam0
                    fb%plasma%x = x0
                end block
            end do
            deallocate (g, gp, gm, du, v, hv, hh, fd)
        end do
        print '(a)', 'slide-gauged Jacobian and vertical quadratic return ok'
    end subroutine

    subroutine hypothesis_checks()
        ! FB-D2 bordering relation as an independent oracle of
        ! ideal_fb_hypotheses: the Sylvester inertia of the slide-bordered KKT
        ! matrix from a dense Bunch-Kaufman LDL^T (after a symmetric diagonal
        ! equilibration, a congruence) must equal In(S_red) + (r, r + dim
        ! Phi~, 0), with In(S_red) from the plasma Schur complement and its QR
        ! gauge basis. The hypotheses (V) and (R) of that relation are checked.
        type(ideal_fb_t) :: fb
        type(ideal_fb_hypotheses_t) :: hyp
        integer, allocatable :: rp(:), ci(:), srp(:), sci(:), ipiv(:)
        real(dp), allocatable :: val(:), sv(:), k(:, :), work(:), scale(:)
        real(dp) :: sscale, query(1), det
        integer :: ntor, st, i, p, m, info, inertia(3)
        interface
            subroutine dsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
                import :: dp
                character, intent(in) :: uplo
                integer, intent(in) :: n, lda, lwork
                real(dp), intent(inout) :: a(lda, *)
                integer, intent(out) :: ipiv(*), info
                real(dp), intent(out) :: work(*)
            end subroutine
        end interface
        do ntor = 0, 1
            call build(16, 12, 2, .false., 3, ntor, fb, solovev=.true.)
            call ideal_fb_set_drive(fb, sampler, st)
            if (st /= IDEAL_FB_OK) error stop 'hypotheses drive'
            call ideal_fb_solve_vacuum(fb, st)
            if (st /= IDEAL_FB_OK) error stop 'hypotheses vacuum'
            call ideal_fb_solve_multipliers(fb, st)
            if (st /= IDEAL_FB_OK) error stop 'hypotheses multipliers'
            call ideal_fb_hypotheses(fb, hyp, st)
            call ideal_fb_hypotheses_print(hyp)
            if (st /= IDEAL_FB_OK) error stop 'hypotheses status'
            if (.not. hyp%s_min > 0) error stop 'S_h not positive'
            if (.not. hyp%vacuum_min > 0) error stop 'A_h + S_h not positive'
            if (.not. hyp%full_rank) error stop '[grad C; M] rank deficient'
            if (hyp%n_zero /= 0) error stop 'S_red singular'
            call ideal_fb_hessian(fb, rp, ci, val)
            sscale = 0
            do i = 1, fb%np
                do p = rp(i), rp(i + 1) - 1
                    if (ci(p) == i) sscale = sscale + abs(val(p))
                end do
            end do
            call ideal_fb_slide_system(fb, rp, ci, val, sscale/fb%np, srp, sci, sv)
            m = size(srp) - 1
            allocate (k(m, m), scale(m), ipiv(m))
            k = 0; scale = 1
            do i = 1, m
                do p = srp(i), srp(i + 1) - 1
                    k(i, sci(p)) = sv(p)
                    if (sci(p) == i .and. abs(sv(p)) > 0) scale(i) = 1/sqrt(abs(sv(p)))
                end do
            end do
            do i = 1, m
                k(:, i) = scale*k(:, i)*scale(i)
            end do
            call dsytrf('L', m, k, m, ipiv, query, -1, info)
            allocate (work(max(1, nint(query(1)))))
            call dsytrf('L', m, k, m, ipiv, work, size(work), info)
            if (info < 0) error stop 'dense KKT factorisation'
            inertia = 0
            i = 1
            do while (i <= m)
                if (ipiv(i) > 0) then
                    inertia = inertia + sign_count(k(i, i), 1)
                    i = i + 1
                else
                    det = k(i, i)*k(i + 1, i + 1) - k(i + 1, i)**2
                    if (det < 0) then
                        inertia(1:2) = inertia(1:2) + 1
                    else
                        inertia = inertia + sign_count(k(i, i) + k(i + 1, i + 1), 2)
                    end if
                    i = i + 2
                end if
            end do
            print '(a,i0,a,i0,a,3(1x,i0),a,3(1x,i0))', 'bordering N=', ntor, &
                ' size ', m, ' dense LDL^T In(K)', inertia, ' predicted', &
                hyp%kkt_inertia
            if (any(inertia /= hyp%kkt_inertia)) error stop 'bordering inertia mismatch'
            deallocate (k, scale, ipiv, work)
        end do
        print '(a)', 'hypotheses and bordering inertia ok'
    end subroutine

    pure function sign_count(d, times) result(add)
        ! (positive, negative, zero) increment of an inertia count.
        real(dp), intent(in) :: d
        integer, intent(in) :: times
        integer :: add(3)
        add = 0
        if (d > 0) then
            add(1) = times
        else if (d < 0) then
            add(2) = times
        else
            add(3) = times
        end if
    end function

    subroutine block_report(name, fd, an)
        character(*), intent(in) :: name
        real(dp), intent(in) :: fd(:), an(:)
        print '(2x,a12,es10.2)', name, maxval(abs(fd - an)) &
            /max(maxval(abs(an)), tiny(1.0_dp))
    end subroutine

    subroutine check_symmetry(rowptr, colidx, val)
        integer, intent(in) :: rowptr(:), colidx(:)
        real(dp), intent(in) :: val(:)
        real(dp) :: worst, big, other
        integer :: i, p, q, j
        worst = 0; big = maxval(abs(val))
        do i = 1, size(rowptr) - 1
            do p = rowptr(i), rowptr(i + 1) - 1
                j = colidx(p); other = 0
                do q = rowptr(j), rowptr(j + 1) - 1
                    if (colidx(q) == i) other = val(q)
                end do
                worst = max(worst, abs(val(p) - other))
            end do
        end do
        print '(a,es9.2)', '  asymmetry/max ', worst/big
        if (worst > 1.0e-12_dp*big) error stop 'KKT matrix not symmetric'
    end subroutine

    subroutine vacuum_ladder()
        type(ideal_fb_t) :: fb
        real(dp) :: errors(3), closure
        integer :: level, st, nc
        windows_on = .false.
        do level = 1, 3
            nc = 48*2**(level - 1)
            call build(nc, nc, degree, .true., 1, 0, fb)
            call ideal_fb_set_drive(fb, sampler, st, closure)
            call ideal_fb_solve_vacuum(fb, st)
            if (st /= IDEAL_FB_OK) error stop 'vacuum solve'
            errors(level) = shell_error(fb)
            print '(a,i4,a,es10.3,a,es10.3)', 'Gamma_c edges', nc, &
                '  rel L2 error', errors(level), '  closure', closure
        end do
        print '(a,2f7.2)', 'observed orders', log(errors(1:2)/errors(2:3))/log(2.0_dp)
        ! The filament field varies on the scale of its distance to Gamma_p;
        ! the observed orders approach k from below on this ladder.
        if (log(errors(2)/errors(3))/log(2.0_dp) < degree - 0.6_dp) &
            error stop 'vacuum order below k - 0.6'
        if (errors(3) > 1.0e-3_dp) error stop 'vacuum field error'
    end subroutine

    real(dp) function shell_error(fb) result(rel)
        ! L2 error of B over the shell against the exact field.
        type(ideal_fb_t), intent(in) :: fb
        type(gauss_quadrature_triangle_t) :: rule
        real(dp) :: x(3), b(3), be(2), psi, err2, ref2, w, shape(10), grad(2, 10)
        real(dp) :: dg(2, 2)
        integer :: e, g, st
        call rule%init(7)
        err2 = 0; ref2 = 0
        do e = 1, fb%ncell
            do g = 1, rule%n_points
                call ideal_fb_shell_field(fb, e, [rule%xi(g), rule%eta(g)], 0.0_dp, &
                    x, b, st)
                if (st /= 0) error stop 'shell field'
                call triangle_basis_jet([rule%xi(g), rule%eta(g)], shape(:fb%nloc), &
                    grad(:, :fb%nloc))
                dg = matmul(fb%reference(:, fb%conn(:, e)), &
                    transpose(grad(:, :fb%nloc)))
                w = rule%weights(g)*abs(dg(1, 1)*dg(2, 2) - dg(1, 2)*dg(2, 1))*x(1)
                call field_ext(x(1), x(2), be, psi)
                err2 = err2 + w*((b(1) - be(1))**2 + (b(3) - be(2))**2 &
                    + (b(2) - f_vac/x(1))**2)
                ref2 = ref2 + w*sum(be**2)
            end do
        end do
        rel = sqrt(err2/ref2)
        block
            real(dp) :: part(3)
            integer :: kind
            part = 0
            do e = 1, fb%ncell
                kind = 3
                if (any(fb%role(fb%conn(:, e)) == 1)) kind = 1
                if (any(fb%role(fb%conn(:, e)) == 3)) kind = 2
                do g = 1, rule%n_points
                    call ideal_fb_shell_field(fb, e, [rule%xi(g), rule%eta(g)], &
                        0.0_dp, x, b, st)
                    call triangle_basis_jet([rule%xi(g), rule%eta(g)], &
                        shape(:fb%nloc), grad(:, :fb%nloc))
                    dg = matmul(fb%reference(:, fb%conn(:, e)), &
                        transpose(grad(:, :fb%nloc)))
                    w = rule%weights(g)*abs(dg(1, 1)*dg(2, 2) - dg(1, 2)*dg(2, 1)) &
                        *x(1)
                    call field_ext(x(1), x(2), be, psi)
                    part(kind) = part(kind) + w*((b(1) - be(1))**2 &
                        + (b(3) - be(2))**2 + (b(2) - f_vac/x(1))**2)
                end do
            end do
            print '(a,3es10.2)', '  error near Gamma_p, Gamma_c, rest', &
                sqrt(part/ref2)
        end block
    end function

    subroutine make_windows(current)
        ! Three rectangular window coils at phi = 0, 120, 240 deg, R = 4.9 m.
        real(dp), intent(in) :: current
        real(dp), allocatable :: pts(:, :)
        integer :: npts(3), c, j, st
        real(dp) :: phi0, corner(2, 5)
        character(8) :: names(3)
        corner = reshape([-0.35_dp, -0.8_dp, 0.35_dp, -0.8_dp, 0.35_dp, 0.8_dp, &
            -0.35_dp, 0.8_dp, -0.35_dp, -0.8_dp], [2, 5])
        allocate (pts(3, 15))
        do c = 1, 3
            phi0 = 2*pi*(c - 1)/3
            do j = 1, 5
                pts(:, 5*(c - 1) + j) = [4.9_dp*cos(phi0 + corner(1, j)/4.9_dp), &
                    4.9_dp*sin(phi0 + corner(1, j)/4.9_dp), corner(2, j) + 0.3_dp]
            end do
            npts(c) = 5
            write (names(c), '(a,i0)') 'w', c
        end do
        call coil_set_from_polylines(windows, npts, pts, names, st)
        if (st /= 0) error stop 'window coils'
        gamma_scale = 1
        if (current < 0) continue
    end subroutine

    subroutine vacuum_3d()
        ! Same shell problem with Gamma_c scaled by 1 and 1.12: the field at
        ! probe points inside the smaller shell must agree.
        type(ideal_fb_t) :: fa, fbig
        real(dp) :: probes(2, 6), ba(3), bb(3), b2(3), worst, scale, closure
        integer :: k, st
        windows_on = .true.
        call make_windows(1.0e5_dp)
        probes = reshape([4.35_dp, 0.0_dp, 4.1_dp, 1.0_dp, 3.0_dp, 1.75_dp, &
            1.8_dp, 0.4_dp, 3.0_dp, -1.75_dp, 4.1_dp, -1.0_dp], [2, 6])
        gamma_scale = 1
        call build(24*degree, 24, degree, .false., 3, 1, fa)
        call ideal_fb_set_drive(fa, sampler, st, closure)
        print '(a,es10.3)', 'trace closure / curl-free mismatch (T, T m)', closure
        call ideal_fb_solve_vacuum(fa, st)
        if (st /= IDEAL_FB_OK) error stop 'vacuum solve a'
        gamma_scale = 1.12_dp
        call build(28*degree, 24, degree, .false., 3, 1, fbig)
        call ideal_fb_set_drive(fbig, sampler, st)
        call ideal_fb_solve_vacuum(fbig, st)
        if (st /= IDEAL_FB_OK) error stop 'vacuum solve b'
        worst = 0; scale = 0
        ! n = 3 part: (B(phi) - B(phi + pi/3))/2 removes n = 0 (and n = 6).
        do k = 1, size(probes, 2)
            call probe(fa, probes(:, k), 0.3_dp, ba)
            call probe(fa, probes(:, k), 0.3_dp + pi/3, b2)
            ba = 0.5_dp*(ba - b2)
            call probe(fbig, probes(:, k), 0.3_dp, bb)
            call probe(fbig, probes(:, k), 0.3_dp + pi/3, b2)
            bb = 0.5_dp*(bb - b2)
            print '(a,2f6.2,a,3es11.3,a,3es11.3)', 'probe', probes(:, k), ' a', ba, &
                ' b', bb
            worst = max(worst, maxval(abs(ba - bb))); scale = max(scale, maxval(abs(ba)))
        end do
        print '(a,es10.3)', 'Gamma_c independence of the n = 3 field |dB|/|B|', &
            worst/scale
        if (worst > 2.0e-2_dp*scale) error stop 'Gamma_c dependence'
    end subroutine

    subroutine probe(fb, y, phi_hat, b)
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(in) :: y(2), phi_hat
        real(dp), intent(out) :: b(3)
        real(dp) :: xi(2), x(3)
        integer :: e, st
        do e = 1, fb%ncell
            call triangle_locate(fb%reference(:, fb%conn(:, e)), y, xi, st)
            if (st /= 0) cycle
            if (minval([xi, 1 - sum(xi)]) < -1.0e-10_dp) cycle
            call ideal_fb_shell_field(fb, e, xi, phi_hat, x, b, st)
            return
        end do
        error stop 'probe not in shell'
    end subroutine
end program test_ideal_free_boundary
