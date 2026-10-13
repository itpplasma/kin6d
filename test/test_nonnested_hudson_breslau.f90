module hb_fixture
    ! Hudson & Breslau (PRL 100, 095001, 2008) model field and heat-flux
    ! problem on the periodic-cylinder chart of the NN2 transport operator.
    !
    ! Field-line Hamiltonian chi = psi^2/2 + sum_j chi_j cos(m_j theta - n_j z)
    ! with the iota = psi profile, A = psi grad theta - chi grad z.  On the
    ! straight cylinder (x, y, z) with the area-preserving flux label
    ! psi = psi_c + (r^2 - 1)/2 (so r dr dtheta = dpsi dtheta) and B_z = 1 the
    ! field B = grad psi x grad theta - grad chi x grad z is exactly
    ! divergence free and has the field lines
    !     d theta/dz = d chi/d psi = psi,  d psi/dz = - d chi/d theta,
    ! i.e. the paper's lines.  Only the metric differs from the paper's
    ! Cartesian (psi, theta) metric (g^{psi psi} = r^2, g^{theta theta} =
    ! 1/r^2, both within 20 percent of 1 on the annulus psi in [0.5, 0.68]).
    ! The paper's perturbation is 2 chi_21 = 3 chi_32 = k, k = 2.1e-3, just
    ! above the critical value 2.039e-3 of the iota = 0.5607 KAM surface.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_state, only: equilibrium_field_t, EQ_BR, EQ_BZ, &
        EQ_BPHI, EQ_VALUES
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_source_t
    use kin6d_nonnested_leaves, only: anisotropic_conductivity
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: psi_c = 0.59_dp, psi_l = 0.5_dp, psi_u = 0.68_dp
    real(dp), parameter :: r_l = sqrt(1 + 2*(psi_l - psi_c))
    real(dp), parameter :: r_u = sqrt(1 + 2*(psi_u - psi_c))
    real(dp), parameter :: big_l = log(r_u/r_l)
    ! Heat flux of the integrable problem (T = ln(r_u/r)/big_l, chi_perp = 1)
    ! per unit axial length (FortFEM forms average over phi): 2 pi/big_l.
    real(dp), parameter :: q_integrable = 2*pi/big_l

    type, extends(equilibrium_field_t) :: hb_field_t
        integer :: nm = 0
        integer :: mpol(4) = 0, ntor(4) = 0
        real(dp) :: amp(4) = 0
    contains
        procedure :: sample => hb_sample
        procedure :: set_hb
    end type

    ! Weak flux load f = K grad w of the lift w = ln(r_u/r)/big_l; with the
    ! discrete solution p, load . p is the heat flux through the annulus.
    type, extends(meridian_fourier_source_t) :: flux_lift_t
        class(hb_field_t), pointer :: field => null()
        real(dp) :: chi_par = 1, chi_perp = 1
    contains
        procedure :: source => flux_lift_source
    end type
contains
    subroutine set_hb(self, k, scheme)
        class(hb_field_t), intent(inout) :: self
        real(dp), intent(in) :: k
        character(*), intent(in) :: scheme
        self%nm = 0
        select case (scheme)
        case ('hb')       ! 2 chi_21 = 3 chi_32 = k
            self%nm = 2
            self%mpol(1:2) = [2, 3]; self%ntor(1:2) = [1, 2]
            self%amp(1:2) = [k/2, k/3]
        case ('island53') ! one isolated (5,3) island at psi = 0.6
            self%nm = 1
            self%mpol(1) = 5; self%ntor(1) = 3; self%amp(1) = k
        end select
    end subroutine

    subroutine hb_sample(field, points, values, statuses)
        class(hb_field_t), intent(in) :: field
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        real(dp) :: x, y, z, r, th, psi, br, bt
        integer :: i, j
        values = 0; statuses = 1
        if (size(points, 1) /= 3) return
        do i = 1, size(points, 2)
            x = points(1, i); y = points(2, i); z = points(3, i)
            r = hypot(x, y); th = atan2(y, x)
            psi = psi_c + (r*r - 1)/2
            br = 0
            do j = 1, field%nm
                br = br + field%mpol(j)*field%amp(j) &
                     *sin(field%mpol(j)*th - field%ntor(j)*z)/r
            end do
            bt = r*psi
            values(EQ_BR, i) = (br*x - bt*y)/r
            values(EQ_BZ, i) = (br*y + bt*x)/r
            values(EQ_BPHI, i) = 1
            statuses(i) = 0
        end do
    end subroutine

    subroutine flux_lift_source(self, t, x, f0, f, status)
        class(flux_lift_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: f0(:), f(:, :)
        integer, intent(out) :: status
        real(dp) :: values(EQ_VALUES, size(x, 2)), k(3, 3), g(3), r2
        integer :: statuses(size(x, 2)), p
        status = 1
        f0 = 0; f = 0
        if (t < 1 .or. .not. associated(self%field)) return
        call self%field%sample(x, values, statuses)
        if (any(statuses /= 0)) return
        do p = 1, size(x, 2)
            call anisotropic_conductivity(values(EQ_BR, p), values(EQ_BPHI, p), &
                values(EQ_BZ, p), self%chi_par, self%chi_perp, k(1, 1), k(1, 2), &
                k(1, 3), k(2, 2), k(2, 3), k(3, 3))
            k(2, 1) = k(1, 2); k(3, 1) = k(1, 3); k(3, 2) = k(2, 3)
            r2 = x(1, p)**2 + x(2, p)**2
            ! (d_R, (1/R) d_phi, d_Z) = (d_x, d_z, d_y) of w(r)
            g = [-x(1, p), 0.0_dp, -x(2, p)]/(r2*big_l)
            f(:, p) = matmul(k, g)
        end do
        status = 0
    end subroutine

    subroutine annulus_map(u, theta, x)
        ! Ring label u in [1/2, 1] -> radius r_l .. r_u; the axis label 0 is
        ! only used for the discarded fan vertex.
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        real(dp) :: r
        r = r_l + (r_u - r_l)*(u - 0.5_dp)/0.5_dp
        x = r*[cos(theta), sin(theta)]
    end subroutine

    subroutine build_annulus(nr, nth, order, fs, status)
        ! Annular mesh r_l <= r <= r_u of nr rings (nr - 1 radial cell layers)
        ! and nth angular cells: the flux mesh without its axis fan.
        integer, intent(in) :: nr, nth, order
        type(meridian_feec_t), intent(out) :: fs
        integer, intent(out) :: status
        type(flux_mesh_t) :: fm
        real(dp), allocatable :: rings(:), vert(:, :), nodes(:, :, :)
        integer, allocatable :: tri(:, :)
        integer :: i, nfan, nt
        status = 1
        if (nr < 2 .or. nth < 3) return
        allocate (rings(nr))
        rings = [(0.5_dp + 0.5_dp*real(i - 1, dp)/real(nr - 1, dp), i=1, nr)]
        call flux_mesh_build(rings, nth, nth, min(order, 3), annulus_map, fm, &
            status)
        if (status /= 0) return
        status = 1
        nfan = fm%nth(1)
        nt = fm%mm%topology%n_triangles - nfan
        tri = fm%mm%topology%triangles(:, nfan + 1:) - 1
        vert = fm%mm%topology%vertices(:, 2:)
        nodes = fm%mm%nodes(:, :, nfan + 1:)
        fm%mm%topology%n_vertices = size(vert, 2)
        fm%mm%topology%n_triangles = nt
        call move_alloc(tri, fm%mm%topology%triangles)
        call move_alloc(vert, fm%mm%topology%vertices)
        call move_alloc(nodes, fm%mm%nodes)
        call meridian_feec_build(fm%mm, order, fs, status, cylinder_r0=1.0_dp)
    end subroutine

    pure real(dp) function t_exact(r) result(t)
        real(dp), intent(in) :: r
        t = log(r_u/r)/big_l
    end function

    pure real(dp) function psi_of_r(r) result(psi)
        real(dp), intent(in) :: r
        psi = psi_c + (r*r - 1)/2
    end function
end module

program test_nonnested_hudson_breslau
    ! Hudson & Breslau (2008) anisotropic heat transport in a chaotic field as
    ! an oracle of the NN2 operator K = chi_par b b + chi_perp (I - b b)
    ! (docs/equilibrium-nonnested-nn2.md, docs/benchmarks.md).
    !   integrable : k = 0.  T = ln(r_u/r)/ln(r_u/r_l) is exact for every
    !     chi_par/chi_perp (grad r is perpendicular to B); L2 error, observed
    !     order and the flux against 4 pi^2/ln(r_u/r_l).
    !   chaos k    : the paper's field at k (default 2.1e-3).  Rigorous
    !     (Dirichlet principle): Q/Q_int >= 1 and nondecreasing in the ratio.
    !     Paper: Q/Q_int about 2 for chi_perp/chi_par = 1e-10.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use hb_fixture
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_lagrange_nodes, &
        meridian_qpoint
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_index, &
        meridian_fourier_value, meridian_fourier_load
    use kin6d_nonnested_transport, only: nt_conductivity_t, nt_report_t, nt_solve
    implicit none
    character(32) :: mode, arg
    real(dp) :: kpert
    integer :: order, nr, nth, modes
    logical :: ok
    call get_command_argument(1, mode)
    order = 3; nr = 7; nth = 48; modes = 3; kpert = 2.1e-3_dp
    call get_command_argument(2, arg)
    if (len_trim(arg) > 0) read (arg, *) order
    call get_command_argument(3, arg)
    if (len_trim(arg) > 0) read (arg, *) nr
    call get_command_argument(4, arg)
    if (len_trim(arg) > 0) read (arg, *) nth
    call get_command_argument(5, arg)
    if (len_trim(arg) > 0) read (arg, *) modes
    call get_command_argument(6, arg)
    if (len_trim(arg) > 0) read (arg, *) kpert
    select case (trim(mode))
    case ('integrable')
        call run_integrable(order, ok)
    case ('chaos')
        call run_chaos(order, nr, nth, modes, kpert, ok)
    case default
        ok = .false.
    end select
    if (.not. ok) error stop 1
    print '(a)', 'PASS'
contains
    subroutine solve(fs, modes, field, ratio, p, report, qflux, st)
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: modes
        type(hb_field_t), intent(in), target :: field
        real(dp), intent(in) :: ratio
        real(dp), allocatable, intent(out) :: p(:)
        type(nt_report_t), intent(out) :: report
        real(dp), intent(out) :: qflux
        integer, intent(out) :: st
        type(nt_conductivity_t) :: cond
        type(flux_lift_t) :: lift
        real(dp), allocatable :: bc(:), xy(:, :), load(:)
        integer :: i, s2
        cond%field => field
        cond%chi_par = ratio; cond%chi_perp = 1
        lift%field => field
        lift%chi_par = ratio; lift%chi_perp = 1
        call meridian_lagrange_nodes(fs, xy)
        allocate (bc(fs%n0*(2*modes + 1)))
        bc = 0
        do i = 1, fs%n0
            if (fs%sp(0)%bnd(i) .and. norm2(xy(:, i)) < 0.5_dp*(r_l + r_u)) &
                bc(meridian_fourier_index(fs%n0, 0, 0, i)) = 1
        end do
        call nt_solve(fs, modes, 4*modes + 4, cond, lift, bc, p, report, st)
        qflux = 0
        if (st /= 0) return
        call meridian_fourier_load(fs, modes, 4*modes + 4, lift, load, s2)
        if (s2 /= 0) then
            st = 1; return
        end if
        qflux = dot_product(load, p)
    end subroutine

    subroutine profile(fs, modes, p, nbin, tbar, tvar, psi)
        ! Volume-binned mean and variance of T in nbin equal-psi bins.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: modes, nbin
        real(dp), intent(in) :: p(:)
        real(dp), intent(out) :: tbar(nbin), tvar(nbin), psi(nbin)
        real(dp) :: x(2), wq, phi(3, fs%nl0), th, w(nbin), s1(nbin), s2(nbin)
        real(dp) :: ps
        integer :: t, q, j, st, ib
        w = 0; s1 = 0; s2 = 0
        do t = 1, fs%mm%topology%n_triangles
            do q = 1, fs%nq
                call meridian_qpoint(fs, 0, t, q, x, wq, phi, st)
                ps = psi_of_r(norm2(x))
                ib = min(nbin, max(1, 1 + int((ps - psi_l)/(psi_u - psi_l)*nbin)))
                do j = 0, 15
                    th = meridian_fourier_value(fs%n0, modes, p, &
                        fs%sp(0)%dof(:, t), phi(1, :), 2*pi*j/16)
                    w(ib) = w(ib) + wq
                    s1(ib) = s1(ib) + wq*th
                    s2(ib) = s2(ib) + wq*th*th
                end do
            end do
        end do
        do ib = 1, nbin
            psi(ib) = psi_l + (ib - 0.5_dp)*(psi_u - psi_l)/nbin
            tbar(ib) = s1(ib)/w(ib)
            tvar(ib) = max(0.0_dp, s2(ib)/w(ib) - tbar(ib)**2)
        end do
    end subroutine

    real(dp) function l2_error_exact(fs, modes, p) result(err)
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: modes
        real(dp), intent(in) :: p(:)
        real(dp) :: x(2), wq, phi(3, fs%nl0), th, te, e2, n2
        integer :: t, q, j, st
        e2 = 0; n2 = 0
        !$omp parallel do default(shared) private(t, q, j, st, x, wq, phi, th, &
        !$omp& te) reduction(+:e2, n2) schedule(static)
        do t = 1, fs%mm%topology%n_triangles
            do q = 1, fs%nq
                call meridian_qpoint(fs, 0, t, q, x, wq, phi, st)
                te = t_exact(norm2(x))
                do j = 0, 3
                    th = meridian_fourier_value(fs%n0, modes, p, &
                        fs%sp(0)%dof(:, t), phi(1, :), 2*pi*j/4)
                    e2 = e2 + wq*(th - te)**2/4
                    n2 = n2 + wq*te**2/4
                end do
            end do
        end do
        !$omp end parallel do
        err = sqrt(e2/n2)
    end function

    subroutine run_integrable(order, ok)
        integer, intent(in) :: order
        logical, intent(out) :: ok
        integer, parameter :: nlev = 3, nratio = 3
        type(meridian_feec_t) :: fs
        type(hb_field_t) :: field
        type(nt_report_t) :: report
        real(dp), allocatable :: p(:)
        real(dp) :: err(nlev, nratio), qerr(nlev, nratio), ratio, qflux
        integer :: lev, ir, st
        ok = .true.
        call field%set_hb(0.0_dp, 'hb')
        do lev = 1, nlev
            call build_annulus(3*2**(lev - 1) + 1, 24*2**(lev - 1), order, fs, st)
            if (st /= 0) then
                ok = .false.; return
            end if
            do ir = 1, nratio
                ratio = 10.0_dp**(4*(ir - 1))
                call solve(fs, 1, field, ratio, p, report, qflux, st)
                if (st /= 0) then
                    ok = .false.; return
                end if
                err(lev, ir) = l2_error_exact(fs, 1, p)
                qerr(lev, ir) = abs(qflux/q_integrable - 1)
                print '(a, i0, a, es8.1, a, i3, a, i7, a, es10.3, a, es10.3, a, &
                    &es9.2)', 'integrable P', order, ' ratio', ratio, ' nr', &
                    3*2**(lev - 1) + 1, ' ndof', report%ndof, ' L2', err(lev, ir), &
                    ' |Q/Q_int-1|', qerr(lev, ir), ' res', report%residual
                ok = ok .and. report%residual <= 1.0e-12_dp
            end do
        end do
        do ir = 1, nratio
            print '(a, i0, a, es8.1, a, f6.2, a, f6.2)', 'integrable P', order, &
                ' ratio', 10.0_dp**(4*(ir - 1)), ' L2 order', &
                log(err(nlev - 1, ir)/err(nlev, ir))/log(2.0_dp), ' Q order', &
                log(qerr(nlev - 1, ir)/qerr(nlev, ir))/log(2.0_dp)
        end do
        ! Exact oracle: P3 order 4 at ratio 1 and 1e4 (measured 3.99/4.01),
        ! error below 1e-4 at the finest mesh up to ratio 1e8 (measured 2e-6).
        ok = ok .and. err(nlev, 1) < 1.0e-6_dp .and. err(nlev, nratio) < 1.0e-4_dp
        ok = ok .and. log(err(nlev - 1, 1)/err(nlev, 1))/log(2.0_dp) > 3.5_dp
    end subroutine

    subroutine run_chaos(order, nr, nth, modes, kpert, ok)
        integer, intent(in) :: order, nr, nth, modes
        real(dp), intent(in) :: kpert
        logical, intent(out) :: ok
        integer, parameter :: nratio = 5, nbin = 18
        type(meridian_feec_t) :: fs
        type(hb_field_t) :: field
        type(nt_report_t) :: report
        real(dp), allocatable :: p(:)
        real(dp) :: ratio, qflux, qprev, tbar(nbin), tvar(nbin), psi(nbin)
        integer :: ir, st, ib
        ok = .true.
        call field%set_hb(kpert, 'hb')
        call build_annulus(nr, nth, order, fs, st)
        if (st /= 0) then
            ok = .false.; return
        end if
        qprev = 0
        do ir = 1, nratio
            ratio = 10.0_dp**(2*(ir - 1))
            call solve(fs, modes, field, ratio, p, report, qflux, st)
            if (st /= 0) then
                ok = .false.; return
            end if
            print '(a, i0, a, es8.1, a, i7, a, f8.4, a, es9.2)', 'chaos P', order, &
                ' ratio', ratio, ' ndof', report%ndof, ' Q/Q_int', &
                qflux/q_integrable, ' res', report%residual
            ok = ok .and. qflux/q_integrable >= 1 - 1.0e-6_dp &
                 .and. qflux >= qprev*(1 - 1.0e-6_dp)
            qprev = qflux
        end do
        call profile(fs, modes, p, nbin, tbar, tvar, psi)
        do ib = 1, nbin
            print '(a, f7.4, a, f8.4, a, es10.3)', 'psi', psi(ib), ' Tbar', &
                tbar(ib), ' var', tvar(ib)
        end do
    end subroutine
end program
