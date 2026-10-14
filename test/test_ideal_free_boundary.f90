module ideal_fb_test_field
    ! Smooth nested frozen data for the derivative checks: elliptic labels
    ! about (3.05, 0), F = 6 T m, pressure linear in the label.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_ideal_map, only: ideal_frozen_t
    implicit none
    private
    type, extends(ideal_frozen_t), public :: simple_field_t
        logical :: solovev = .false.
    contains
        procedure :: evaluate => simple_evaluate
    end type
contains
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
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_ideal_map, only: ideal_frozen_t
    use kin6d_ideal_free_boundary, only: ideal_fb_t, ideal_fb_init, &
        ideal_fb_set_drive, ideal_fb_gradient, ideal_fb_hessian, ideal_fb_step, &
        ideal_fb_solve_vacuum, ideal_fb_shell_field, IDEAL_FB_OK, &
        ideal_fb_vertical_system, ideal_fb_vertical_gradient, ideal_fb_vertical_expand, &
        ideal_fb_newton_step, ideal_fb_solve_shell_map, ideal_fb_solve_multipliers, &
        ideal_fb_solve, ideal_fb_report_t
    use kin6d_curved_triangle, only: triangle_locate
    use kin6d_coils_3d, only: coil_set_t, coil_set_from_polylines, coil_set_modes
    use fortfem_axisymmetric_gs_kernel, only: gs_green_gradient, gs_green
    use fortfem_gauss_quadrature_2d, only: gauss_quadrature_triangle_t
    use kin6d_triangle_basis, only: triangle_basis_jet
    use ideal_fb_test_field, only: simple_field_t
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
    mode = 'derivatives'; degree = 2
    if (command_argument_count() > 0) call get_command_argument(1, mode)
    if (command_argument_count() > 1) then
        call get_command_argument(2, argument)
        read (argument, *) degree
    end if
    select case (trim(mode))
    case ('derivatives')
        call derivative_checks()
    case ('vertical')
        call vertical_checks()
    case ('vacuum')
        call vacuum_ladder()
    case ('vacuum3d')
        call vacuum_3d()
    case default
        error stop 'unknown mode'
    end select
contains
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
        if (flux_chain) then
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

    subroutine derivative_checks()
        type(ideal_fb_t) :: fb
        integer, allocatable :: rowptr(:), colidx(:)
        real(dp), allocatable :: val(:), g(:), gp(:), gm(:), v(:), hv(:), du(:)
        real(dp) :: lp, lm, l0, h, fd, an, err, scale
        integer :: ntor, i, st, p
        do ntor = 0, 1
            windows_on = ntor > 0
            if (windows_on) call make_windows(1.0e5_dp)
            call build(24, 16, 2, .false., 3, ntor, fb)
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
        ! equations, not a claim of a coil-supported analytic free boundary.
        type(ideal_fb_t) :: fb
        type(ideal_fb_report_t) :: report
        integer, allocatable :: rp(:), ci(:), vrp(:), vci(:)
        real(dp), allocatable :: val(:), vv(:), g(:), gp(:), gm(:), v(:), du(:)
        real(dp), allocatable :: reduced(:), plus(:), minus(:), hv(:), x0(:, :, :)
        real(dp) :: h, error(0:3), residual(0:3), scale, fd_error, dz
        integer :: ntor, i, p, st, a, z, it, sign_z
        do ntor = 0, 1
            call build(16, 12, 2, .false., 3, ntor, fb, solovev=.true.)
            call ideal_fb_set_drive(fb, sampler, st)
            if (st /= IDEAL_FB_OK) error stop 'vertical drive'
            fb%rigid_vertical = .true.; fb%tangential_gauge = 100
            allocate (g(fb%n), gp(fb%n), gm(fb%n), du(fb%n), v(fb%n + 1), &
                      reduced(fb%n + 1), plus(fb%n + 1), minus(fb%n + 1), hv(fb%n + 1))
            ! Nonaxisymmetric perturbations exercise the off-diagonal blocks.
            do i = 1, fb%n
                du(i) = 2.0e-4_dp*sin(0.71_dp*i)
            end do
            call ideal_fb_step(fb, du, 1.0_dp)
            call ideal_fb_gradient(fb, g, st)
            call ideal_fb_hessian(fb, rp, ci, val)
            call ideal_fb_vertical_system(fb, rp, ci, val, val, vrp, vci, vv)
            call check_symmetry(vrp, vci, vv)
            v = [(cos(0.31_dp*i), i=1, fb%n + 1)]
            do a = 1, fb%plasma%nodes
                z = 2 + 2*fb%nmode*(a - 1)
                v(z) = 0
            end do
            do a = 1, fb%ns
                i = fb%np + fb%nx + 2 + 3*fb%nmode*(a - 1)
                v(i) = 1.0e3_dp
            end do
            v(fb%n + 1) = 0.73_dp
            call ideal_fb_vertical_expand(fb, v, du)
            h = 1.0e-6_dp
            call ideal_fb_step(fb, du, h)
            call ideal_fb_gradient(fb, gp, st)
            call ideal_fb_vertical_gradient(fb, gp, plus)
            call ideal_fb_step(fb, du, -2*h)
            call ideal_fb_gradient(fb, gm, st)
            call ideal_fb_vertical_gradient(fb, gm, minus)
            call ideal_fb_step(fb, du, h)
            hv = 0
            do i = 1, fb%n + 1
                do p = vrp(i), vrp(i + 1) - 1
                    hv(i) = hv(i) + vv(p)*v(vci(p))
                end do
            end do
            scale = maxval(abs(hv))
            fd_error = maxval(abs((plus - minus)/(2*h) - hv))/scale
           print '(a,i0,a,es12.4)', 'vertical N=', ntor, ' Jacobian FD error ', fd_error
            if (fd_error > 1.0e-6_dp) error stop 'vertical Jacobian mismatch'
            ! Explicitly check the force row on its own scale.
            fd_error = abs((plus(fb%n + 1) - minus(fb%n + 1))/(2*h) &
                           - hv(fb%n + 1))/max(abs(hv(fb%n + 1)), 1.0_dp)
            print '(a,es12.4)', 'vertical force row FD error ', fd_error
            if (fd_error > 1.0e-5_dp) error stop 'vertical force Jacobian mismatch'
            ! The affine manufactured load makes this exact discrete target.
            call ideal_fb_gradient(fb, g, st)
            fb%shift = g
            call ideal_fb_solve(fb, 1.0e-10_dp, 2, report, st, &
                residual_tolerance=1.0e-12_dp)
            if (st /= IDEAL_FB_OK .or. report%iterations /= 0) &
                error stop 'exact manufactured target not recognized'
            x0 = fb%plasma%x
            do sign_z = -1, 1, 2
                dz = sign_z*0.004_dp
                fb%plasma%x(2, 1, :) = x0(2, 1, :) + dz
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
                    do it = 0, 3
                        error(it) = max(maxval(abs(fb%plasma%x - x0)), &
                            maxval(abs(fb%xs - xs0)), &
                            maxval(abs(fb%phi - phi0))/f_vac, &
                            maxval(abs(fb%lam - lam0))/1.0e6_dp)
                        call ideal_fb_gradient(fb, g, st)
                        call ideal_fb_vertical_gradient(fb, g, reduced)
                        residual(it) = maxval(abs(reduced))
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
            deallocate (g, gp, gm, du, v, reduced, plus, minus, hv)
        end do
        print '(a)', 'vertical coupled Jacobian and quadratic return ok'
    end subroutine

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
