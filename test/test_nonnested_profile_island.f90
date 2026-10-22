module profile_island_fixture
    ! Known-perturbation oracle, not a nonlinear saturation prediction.
    ! Canonical plane u=(R^2-R0^2)/(2 R0), v=Z, action
    ! x=(u^2+v^2-rs^2)/(2 rs), alpha=theta-phi. B_phi=R makes
    ! phi the flow time. The helical-frame Hamiltonian is
    ! H=x^2/2+eps cosh(x) cos(alpha), and dt= dphi/rs.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_state, only: equilibrium_field_t, EQ_VALUES, &
                                       EQ_BR, EQ_BZ, EQ_BPHI, EQ_P
    use kin6d_nonnested_profile, only: nn_profile_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), rs = 0.7_dp
    real(dp) :: major = 6
    type, extends(equilibrium_field_t) :: island_field_t
        real(dp) :: r0 = 6, eps = 1.0e-5_dp
    contains
        procedure :: sample => island_sample
    end type
    type, extends(nn_profile_t) :: island_profile_t
        real(dp) :: r0 = 6, eps = 1.0e-5_dp
    contains
        procedure :: evaluate => island_profile
    end type
contains
    subroutine mapping(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [sqrt(major**2 + 2*major*u*cos(theta)), u*sin(theta)]
    end subroutine
    subroutine kernel(r0, eps, point, b, p, gp)
        real(dp), intent(in) :: r0, eps, point(3)
        real(dp), intent(out) :: b(3), p, gp(3)
        real(dp) :: u, v, radius2, action, alpha, hu, hv, hx, ha
        u = (point(1)**2 - r0*r0)/(2*r0); v = point(2)
        radius2 = u*u + v*v
        action = (radius2 - rs*rs)/(2*rs)
        alpha = atan2(v, u) - point(3)
        p = action**2/2 + eps*cosh(action)*cos(alpha)
        hx = action + eps*sinh(action)*cos(alpha)
        ha = -eps*cosh(action)*sin(alpha)
        hu = hx*u/rs - ha*v/radius2
        hv = hx*v/rs + ha*u/radius2
        gp = [hu*point(1)/r0, -ha/point(1), hv]
        ! Full Hamiltonian adds rs*x to cancel the rotating-frame term.
        hu = hu + u; hv = hv + v
        b = [-r0*hv/point(1), point(1), hu]
    end subroutine
    subroutine island_sample(field, points, values, statuses)
        class(island_field_t), intent(in) :: field
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        real(dp) :: b(3), p, gp(3)
        integer :: i
        values = 0; statuses = 0
        do i = 1, size(points, 2)
            call kernel(field%r0, field%eps, points(:, i), b, p, gp)
            values(EQ_BR, i) = b(1); values(EQ_BPHI, i) = b(2)
            values(EQ_BZ, i) = b(3); values(EQ_P, i) = p
        end do
    end subroutine
    subroutine island_profile(self, x, p, gradp, lambda, gradphi, status)
        class(island_profile_t), intent(in) :: self
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: p(:), gradp(:, :), lambda(:), gradphi(:, :)
        integer, intent(out) :: status
        real(dp) :: b(3)
        integer :: i
        lambda = 0; gradphi = 0; status = 0
        do i = 1, size(x, 2)
            call kernel(self%r0, self%eps, x(:, i), b, p(i), gradp(:, i))
        end do
    end subroutine
end module

program test_nonnested_profile_island
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use profile_island_fixture
    use fortnum_ode_cash_karp, only: cash_karp_step
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes, meridian_qpoint
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_value
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_build
    use kin6d_nonnested_profile, only: nn_fitted_source_t, nn_fit_source
    use kin6d_nonnested_transport, only: nt_conductivity_t, nt_report_t, nt_solve
    implicit none
    real(dp) :: eps, width(4), exact_width, lo, hi, x, slope, physical, local_width
    real(dp) :: rminus, rplus, phase, phase_x, d, err, olderr, outer(3)
    integer :: i, j, lev
    ! Independent outer-region finite-difference BVP U''-U=0, wall
    ! U(+-0.15)=cosh(0.15). It fixes the resonant perturbation U(0)=1.
    do lev = 1, 3
        call outer_solution(16*2**(lev - 1), outer(lev))
    end do
    if (abs(outer(3) - 1) > 1.0e-6_dp) error stop 'independent outer solution'
    if (abs(outer(3) - 1) >= abs(outer(1) - 1)) error stop 'outer convergence'
    print '(a,3es13.4)', 'outer resonant amplitude errors ', abs(outer - 1)
    do j = 1, 3
        ! First canonical slab/cylinder, then circular toroidal geometry at
        ! aspect ratio 6 and 60. The physical radial conversion is explicit.
        major = merge(6.0_dp, 60.0_dp, j == 2)
        do i = 1, 4
            eps = 1.0e-6_dp*4**(i - 1)
            lo = 0; hi = 0.1_dp
            do lev = 1, 50
                x = (lo + hi)/2
                if (x*x/2 - eps*cosh(x) < eps) then
                    lo = x
                else
                    hi = x
                end if
            end do
            x = (lo + hi)/2
            width(i) = 2*x
            exact_width = 4*sqrt(eps*outer(3))
            physical = width(i); local_width = exact_width
            if (j > 1) then
                rminus = sqrt(rs*rs - 2*rs*x)
                rplus = sqrt(rs*rs + 2*rs*x)
                physical = sqrt(major**2 - 2*major*rminus) &
                           - sqrt(major**2 - 2*major*rplus)
                local_width = exact_width*major/sqrt(major**2 - 2*major*rs)
                width(i) = physical
            end if
            if (abs(physical - local_width) > max(0.01_dp*local_width, 1.0e-4_dp)) &
                error stop 'G4 width'
            ! Independent time integration straddles the contour-derived
            ! separatrix. Use a helical frame and scaled time to avoid a
            ! trace duration proportional to 1/sqrt(eps).
            call orbit(0.995_dp*x/sqrt(eps), eps, 2048, phase)
            call orbit(1.005_dp*x/sqrt(eps), eps, 2048, phase_x)
            if (phase > pi .or. phase_x < pi) error stop 'G4 separatrix trace'
            call orbit(0.995_dp*x/sqrt(eps), eps, 4096, d)
            if (abs(d - phase) > 1.0e-5_dp) error stop 'G4 trace refinement'
            ! The O/X points of the generated flow are alpha=pi and 0.
            ! Check their flow and Hessian classification, not a coil phase copy.
            call fixed_points(eps, phase, phase_x)
        if (abs(phase - pi) > 0.02_dp .or. abs(phase_x) > 0.02_dp) error stop 'G4 phase'
            print '(a,2i3,4es13.4)', 'G4 geometry/amplitude eps/width/oracle/phase ', &
                j, i, eps, physical, local_width, phase
        end do
        slope = log(width(4)/width(1))/log(64.0_dp)
        if (abs(slope - 0.5_dp) > 0.05_dp) error stop 'G4 sqrt amplitude'
        print '(a,f9.5)', 'G4 log slope ', slope
        if (j == 1) cycle
        olderr = huge(1.0_dp)
        do lev = 1, 2
            call fitted_transport(4*2**lev, eps, err)
            if (err >= olderr) error stop 'island fitted target convergence'
            olderr = err
        end do
        if (err > 0.01_dp) error stop 'island fitted pressure target'
    end do
    print *, 'PASS: G4 known-forced-field geometry and frozen pressure fit'
contains
    subroutine outer_solution(n, center)
        integer, intent(in) :: n
        real(dp), intent(out) :: center
        real(dp) :: h, dl(2*n - 2), du(2*n - 2), diag(2*n - 1), rhs(2*n - 1, 1)
        integer :: info
        interface
            subroutine dgtsv(n, nrhs, dl, d, du, b, ldb, info)
                import :: dp
                integer, intent(in) :: n, nrhs, ldb
                real(dp), intent(inout) :: dl(*), d(*), du(*), b(ldb, *)
                integer, intent(out) :: info
            end subroutine
        end interface
        h = 0.15_dp/n
        dl = -1; du = -1; diag = 2 + h*h; rhs = 0
        rhs(1, 1) = cosh(0.15_dp); rhs(2*n - 1, 1) = cosh(0.15_dp)
        call dgtsv(2*n - 1, 1, dl, diag, du, rhs, 2*n - 1, info)
        if (info /= 0) error stop 'outer BVP solve'
        center = rhs(n, 1)
    end subroutine
    subroutine rhs(t, y, dy, ctx)
        real(dp), intent(in) :: t, y(:)
        real(dp), intent(out) :: dy(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: eps, x
        eps = 0
        if (present(ctx)) then
            select type (ctx)
            type is (real(dp))
                eps = ctx
            end select
        end if
        x = sqrt(eps)*y(1)
        dy = [cosh(x)*sin(y(2)), y(1) + sqrt(eps)*sinh(x)*cos(y(2))]
    end subroutine
    subroutine orbit(start, eps, steps, span)
        real(dp), intent(in) :: start, eps
        integer, intent(in) :: steps
        real(dp), intent(out) :: span
        real(dp) :: y(2), stage(2, 6), tmp(2), next(2), error(2), h
        integer :: i, nfev
        y = [start, pi]; span = 0; nfev = 0; h = 40.0_dp/steps
        do i = 1, steps
          call cash_karp_step(rhs, (i - 1)*h, y, h, .false., stage(:, 1), stage(:, 2), &
    stage(:, 3), stage(:, 4), stage(:, 5), stage(:, 6), tmp, next, error, nfev, ctx=eps)
            y = next
            span = max(span, abs(y(2) - pi))
        end do
    end subroutine
    subroutine fixed_points(eps, ophase, xphase)
        real(dp), intent(in) :: eps
        real(dp), intent(out) :: ophase, xphase
        real(dp) :: left, right, mid, f
        integer :: i
        ! Roots of alpha force, bracketed independently around each point.
        left = 2; right = 4
        do i = 1, 45
            mid = (left + right)/2; f = sin(mid)
            if (f > 0) then
                left = mid
            else
                right = mid
            end if
        end do
        ophase = (left + right)/2
        left = -1; right = 1
        do i = 1, 45
            mid = (left + right)/2
            if (sin(mid) < 0) then
                left = mid
            else
                right = mid
            end if
        end do
        xphase = (left + right)/2
        if ((1 - eps)*eps*cos(ophase) >= 0) error stop 'O stability'
        if ((1 + eps)*eps*cos(xphase) <= 0) error stop 'X stability'
    end subroutine
    subroutine fitted_transport(nr, eps, error)
        integer, intent(in) :: nr
        real(dp), intent(in) :: eps
        real(dp), intent(out) :: error
        type(flux_mesh_t) :: mesh
        type(mesh_2d_t) :: topo
        type(free_b_field_t) :: space
        type(island_field_t), target :: seed
        type(island_profile_t) :: target
        type(nn_fitted_source_t) :: fit
        type(nt_conductivity_t) :: tensor
        type(nt_report_t) :: report
        real(dp), allocatable :: nodes(:, :), bc(:), p(:)
        real(dp) :: rings(nr + 1), u, v, action, r, xq(2), wq, basis(3, 10)
        real(dp) :: ph, point(3), pe, gp(3), b(3), e2, p2, par
        integer :: i, fan, st, t, q, j
        rings = [(0.4_dp + 0.6_dp*(i - 1)/nr, i=1, nr + 1)]
        call flux_mesh_build(rings, 8*nr, 8*nr, 3, mapping, mesh, st)
        if (st /= 0) error stop 'island mesh'
        fan = mesh%tlast(1)
        topo%n_triangles = mesh%mm%topology%n_triangles - fan
        topo%n_vertices = mesh%mm%topology%n_vertices - 1; topo%has_triangles = .true.
        topo%vertices = mesh%mm%topology%vertices(:, 2:)
        topo%triangles = mesh%mm%topology%triangles(:, fan + 1:) - 1
        mesh%mm%topology = topo
        mesh%mm%nodes = mesh%mm%nodes(:, :, fan + 1:)
        mesh%tlast = mesh%tlast - fan; mesh%vfirst = mesh%vfirst - 1
        mesh%ring_tri = mesh%ring_tri - fan
        call nn_build(mesh, 3, 1, [major, 0.0_dp], 1.0_dp, 0.0_dp, 1.0_dp, space, st)
        if (st /= 0) error stop 'island space'
        seed%r0 = major; seed%eps = eps
        target%r0 = major; target%eps = eps
     call nn_fit_source(seed, target, 100.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, fit, st)
        if (st /= 0) error stop 'island fit'
        call meridian_lagrange_nodes(space%fs, nodes)
        allocate (bc(3*space%fs%n0)); bc = 0
        do i = 1, space%fs%n0
            u = (nodes(1, i)**2 - major**2)/(2*major); v = nodes(2, i)
            r = sqrt(u*u + v*v); action = (r*r - rs*rs)/(2*rs)
            bc(i) = action*action/2
            bc(space%fs%n0 + i) = eps*cosh(action)*u/r
            bc(2*space%fs%n0 + i) = eps*cosh(action)*v/r
        end do
        tensor%field => seed; tensor%chi_par = 100; tensor%chi_perp = 1
        call nt_solve(space%fs, 1, 8, tensor, fit, bc, p, report, st)
        if (st /= 0) error stop 'island frozen transport solve'
        e2 = 0; p2 = 0; par = 0
        do t = 1, space%fs%mm%topology%n_triangles
            do q = 1, space%fs%nq
                call meridian_qpoint(space%fs, 0, t, q, xq, wq, basis, st)
                if (st /= 0) error stop 'island quadrature'
                do j = 0, 15
                    ph = 2*pi*j/16; point = [xq, ph]
                    call kernel(major, eps, point, b, pe, gp)
                    par = max(par, abs(dot_product(b, gp)))
                    r = meridian_fourier_value(space%fs%n0, 1, p, &
                                              space%fs%sp(0)%dof(:, t), basis(1, :), ph)
                    e2 = e2 + wq*(r - pe)**2/16; p2 = p2 + wq*pe**2/16
                end do
            end do
        end do
        error = sqrt(e2/p2)
     print '(a,f6.1,i5,3es13.4)', 'G4 frozen p target R0/nr/error/parallel/residual ', &
            major, nr, error, par, report%residual
        if (par > 1.0e-10_dp) error stop 'island invariant pressure'
    end subroutine
end program
