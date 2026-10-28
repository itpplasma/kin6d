program test_ideal_map_parallel
    ! Parallel-only Levenberg-Marquardt term of the ideal map Newton
    ! (docs/relaxation-comparators.md item 2; docs/ideal-map-3d.md).
    !   unit [n]         : M_par is symmetric, 0 <= M_par <= M_total, nearly
    !                      vanishes on transversal (nu) increments of the
    !                      identity map (b . grad chi = 0 at the nodes) and is
    !                      about <b_phi^2> M_total on omega increments.
    !   compare [n] [k]  : nested A3 seed, P_k, one mesh of n boundary nodes:
    !                      Newton steps, energy, exact-flux error and the
    !                      generalized spectrum of the converged Hessian
    !                      (H v = lambda M_total v): lowest curvatures, their
    !                      parallel fraction v.M_par v and the condition
    !                      numbers of the trial matrices with diagonal and
    !                      with parallel damping.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_ideal_map, only: ideal_map_t, ideal_map_init, ideal_map_freeze, &
        ideal_map_hessian, ideal_map_solve, ideal_map_solve_report_t, &
        ideal_map_parallel_mass, ideal_map_rdof, ideal_map_point
    use kin6d_ideal_seed, only: nested_seed_t, ideal_seed_init
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_gauss_quadrature_2d, only: gauss_quadrature_triangle_t
    use ideal_map_cases, only: setup_case, boundary_point, phi_edge, &
        solovev_profile_t, pi, a, psi_exact, psi_axis
    use kin6d_triangle_basis, only: triangle_basis_jet, TRIANGLE_MAX_NODES
    implicit none
    character(16) :: mode, arg
    integer :: n = 32, degree = 2
    type(nested_seed_t) :: seed
    type(solovev_profile_t) :: profile
    integer :: status
    call get_command_argument(1, mode)
    if (command_argument_count() > 1) then
        call get_command_argument(2, arg); read (arg, *) n
    end if
    if (command_argument_count() > 2) then
        call get_command_argument(3, arg); read (arg, *) degree
    end if
    call setup_case(.false.)
    call ideal_seed_init(seed, boundary_point, phi_edge(), profile, status)
    if (status /= 0) error stop 'seed'
    select case (trim(mode))
    case ('unit')
        call unit_checks()
    case ('compare')
        call compare()
    case default
        error stop 'mode'
    end select
    print '(a)', 'PASS'
contains
    subroutine make_map(map)
        type(ideal_map_t), intent(out) :: map
        type(gs_state_t) :: gs
        type(cylindrical_chart_t) :: chart
        real(dp), allocatable :: contour(:, :)
        real(dp) :: h
        integer :: j
        allocate (contour(2, n))
        do j = 1, n
            call boundary_point(real(j - 1, dp)/n, contour(:, j))
        end do
        h = 2*pi*a/n
        call initialize_gs(gs, contour, h*h/2, status, degree, curve=boundary_point, &
            geometry_only=.true.)
        if (status /= 0) error stop 'mesh'
        call ideal_map_init(map, gs, chart, 1, 0, status, quadrature_order=7)
        if (status /= 0) error stop 'map init'
        if (degree >= 3) map%gauge_weight = 1.0e-2_dp
        call ideal_map_freeze(map, seed, status)
        if (status /= 0) error stop 'freeze'
    end subroutine

    function quad(m, v) result(s)
        type(csr_matrix_t), intent(in) :: m
        real(dp), intent(in) :: v(:)
        real(dp) :: s
        integer :: i, k
        s = 0
        do i = 1, m%n
            do k = m%rowptr(i), m%rowptr(i + 1) - 1
                s = s + v(i)*m%val(k)*v(m%colidx(k))
            end do
        end do
    end function

    subroutine unit_checks()
        type(ideal_map_t) :: map
        type(csr_matrix_t) :: mp, mt
        real(dp), allocatable :: v(:), w(:)
        real(dp) :: ratio, asym
        integer :: trial, i, k, node, kk
        call make_map(map)
        call ideal_map_parallel_mass(map, map%x, mp)
        call ideal_map_parallel_mass(map, map%x, mt, total=.true.)
        if (size(mp%val) /= size(mt%val)) error stop 'pattern'
        asym = 0
        do i = 1, mp%n
            do k = mp%rowptr(i), mp%rowptr(i + 1) - 1
                do kk = mp%rowptr(mp%colidx(k)), mp%rowptr(mp%colidx(k) + 1) - 1
                    if (mp%colidx(kk) == i) asym = max(asym, abs(mp%val(k) - mp%val(kk)))
                end do
            end do
        end do
        print '(a,es12.4)', 'M_par asymmetry ', asym/maxval(abs(mp%val))
        if (asym > 1.0e-14_dp*maxval(abs(mp%val))) error stop 'M_par symmetry'
        allocate (v(map%nred), w(map%nred))
        do trial = 1, 5
            call random_number(v); v = v - 0.5_dp
            ratio = quad(mp, v)/quad(mt, v)
            print '(a,es12.4)', 'random v.M_par v / v.M_total v ', ratio
            if (ratio < -1.0e-14_dp .or. ratio > 1 + 1.0e-12_dp) error stop 'bounds'
            ! Transversal increments of the identity map are perpendicular to
            ! B at the nodes (d_a || grad chi); between nodes the interpolated
            ! direction leaves an O(h^2) parallel part.
            w = 0
            do node = 1, map%nodes
                w(ideal_map_rdof(map, 1, 1, node)) = v(ideal_map_rdof(map, 1, 1, node))
            end do
            ratio = quad(mp, w)/quad(mt, w)
            print '(a,es12.4)', 'nu-only ratio (b . grad chi = 0) ', ratio
            if (ratio < 0 .or. ratio > 1.0e-3_dp) error stop 'transversal not perpendicular'
            w = 0
            do node = 1, map%nodes
                w(ideal_map_rdof(map, 2, 1, node)) = v(ideal_map_rdof(map, 2, 1, node))
            end do
            ratio = quad(mp, w)/quad(mt, w)
            print '(a,es12.4)', 'omega-only ratio (b_phi^2) ', ratio
            if (ratio < 0.9_dp .or. ratio > 1) error stop 'omega parallel fraction'
        end do
    end subroutine

    subroutine solve(map, kappa, report, error)
        type(ideal_map_t), intent(inout) :: map
        real(dp), intent(in) :: kappa
        type(ideal_map_solve_report_t), intent(out) :: report
        real(dp), intent(out) :: error
        integer(int64) :: t0, t1, rate
        call make_map(map)
        map%parallel_damping = kappa
        call system_clock(t0, rate)
        call ideal_map_solve(map, 1.0e-11_dp, 600, report, status, .false.)
        call system_clock(t1)
        if (status /= 0) error stop 'solve'
        error = flux_error(map)
        print '(a,es9.2,a,i4,a,es24.16,a,es12.4,a,f8.2)', 'kappa ', kappa, &
            ' newton ', report%iterations, ' W ', report%energy, ' d3D_exact ', &
            error, ' wall_s ', real(t1 - t0, dp)/rate
    end subroutine

    real(dp) function flux_error(map)
        ! ||psi_exact(X) - chi|| / (|psi_axis| sqrt(A)) over the reference
        ! (errors(1) of test_ideal_map_axisymmetric's defects).
        type(ideal_map_t), intent(in) :: map
        type(gauss_quadrature_triangle_t) :: rule
        real(dp) :: q(3), y(2), chi, ex, g(2), area, w, b3, p, s, det, dg(2, 2)
        real(dp) :: shape(TRIANGLE_MAX_NODES), dshape(2, TRIANGLE_MAX_NODES)
        integer :: e, k
        call rule%init(7)
        s = 0; area = 0
        do e = 1, map%cells
            do k = 1, rule%n_points
                call ideal_map_point(map, e, [rule%xi(k), rule%eta(k)], 0.0_dp, q, y)
                call triangle_basis_jet([rule%xi(k), rule%eta(k)], &
                    shape(:map%nloc), dshape(:, :map%nloc))
                dg = matmul(map%reference(:, map%conn(:, e)), &
                    transpose(dshape(:, :map%nloc)))
                det = dg(1, 1)*dg(2, 2) - dg(1, 2)*dg(2, 1)
                w = rule%weights(k)*det
                call seed%evaluate(y, chi, g, b3, p)
                call psi_exact(q(1), q(2), ex, g)
                s = s + w*(ex - chi)**2; area = area + w
            end do
        end do
        flux_error = sqrt(s/area)/abs(psi_axis)
    end function

    subroutine dense_free(m, free, d)
        type(csr_matrix_t), intent(in) :: m
        logical, intent(in) :: free(:)
        real(dp), allocatable, intent(out) :: d(:, :)
        integer, allocatable :: idx(:)
        integer :: i, k, nf
        nf = count(free)
        allocate (idx(m%n)); idx = 0
        idx(pack([(i, i=1, m%n)], free)) = [(i, i=1, nf)]
        allocate (d(nf, nf)); d = 0
        do i = 1, m%n
            if (.not. free(i)) cycle
            do k = m%rowptr(i), m%rowptr(i + 1) - 1
                if (free(m%colidx(k))) d(idx(i), idx(m%colidx(k))) = m%val(k)
            end do
        end do
        d = (d + transpose(d))/2
    end subroutine

    subroutine generalized(a, b, w, vectors)
        real(dp), intent(in) :: a(:, :), b(:, :)
        real(dp), allocatable, intent(out) :: w(:)
        real(dp), allocatable, intent(out), optional :: vectors(:, :)
        interface
            subroutine dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, &
                    lwork, info)
                import :: dp
                integer, intent(in) :: itype, n, lda, ldb, lwork
                character, intent(in) :: jobz, uplo
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
                real(dp), intent(out) :: w(*), work(*)
                integer, intent(out) :: info
            end subroutine
        end interface
        real(dp), allocatable :: aa(:, :), bb(:, :), work(:)
        integer :: nf, info
        character :: jobz
        nf = size(a, 1)
        aa = a; bb = b
        allocate (w(nf), work(max(1, 66*nf)))
        jobz = merge('V', 'N', present(vectors))
        call dsygv(1, jobz, 'U', nf, aa, nf, bb, nf, w, work, size(work), info)
        if (info /= 0) error stop 'dsygv'
        if (present(vectors)) vectors = aa
    end subroutine

    subroutine compare()
        type(ideal_map_t) :: map
        type(ideal_map_solve_report_t) :: report
        type(csr_matrix_t) :: h, mp, mt
        real(dp), allocatable :: dh(:, :), dp_(:, :), dt(:, :), w(:), v(:, :)
        real(dp), allocatable :: diag(:, :), wd(:)
        real(dp) :: kappas(3) = [0.0_dp, 1.0e-2_dp, 1.0e-1_dp], errors(3), energies(3)
        real(dp) :: scale, pscale, mu
        integer :: i, j, iterations(3)
        logical, allocatable :: free(:)
        do i = 3, 1, -1
            call solve(map, kappas(i), report, errors(i))
            iterations(i) = report%iterations; energies(i) = report%energy
            ! Lowest curvatures at this end point (positive: a local minimum).
            call ideal_map_hessian(map, map%x, h)
            call ideal_map_parallel_mass(map, map%x, mt, total=.true.)
            free = .not. map%fixed
            call dense_free(h, free, dh)
            call dense_free(mt, free, dt)
            if (size(dh, 1) <= 6000) then
                call generalized(dh, dt, w)
                print '(a,3es12.4)', '  end point lowest curvatures ', w(1:3)
                if (w(1) <= 0) error stop 'end point is not a local minimum'
            end if
        end do
        ! map holds the baseline (kappa = 0) solution.
        call ideal_map_hessian(map, map%x, h)
        call ideal_map_parallel_mass(map, map%x, mp)
        call ideal_map_parallel_mass(map, map%x, mt, total=.true.)
        free = .not. map%fixed
        call dense_free(h, free, dh)
        call dense_free(mp, free, dp_)
        call dense_free(mt, free, dt)
        if (size(dh, 1) > 6000) return
        call generalized(dh, dt, w, v)
        print '(a,i0,a,es12.4)', 'free dofs ', size(w), '  largest curvature ', w(size(w))
        print '(a)', ' i  curvature   /largest   parallel fraction v.Mpar v (v.Mt v = 1)'
        do i = 1, 8
            print '(i3,3es12.4)', i, w(i), w(i)/w(size(w)), &
                dot_product(v(:, i), matmul(dp_, v(:, i)))
        end do
        ! Trial-matrix conditioning in the M_total metric.
        allocate (diag(size(dh, 1), size(dh, 1))); diag = 0
        do j = 1, size(dh, 1)
            diag(j, j) = abs(dh(j, j))
        end do
        scale = maxval(abs([(dh(j, j), j=1, size(dh, 1))]))
        pscale = maxval([(dp_(j, j), j=1, size(dh, 1))])
        do i = 1, 4
            mu = 10.0_dp**(-2*i - 2)
            call generalized(dh + mu*diag, dt, wd)
            print '(a,es9.2,a,es12.4,a,es12.4)', 'diagonal LM mu ', mu, &
                '  lambda_min ', wd(1), '  cond ', wd(size(wd))/wd(1)
        end do
        do i = 1, 4
            mu = 10.0_dp**(-i - 1)
            call generalized(dh + (mu*scale/pscale)*dp_, dt, wd)
            print '(a,es9.2,a,es12.4,a,es12.4)', 'parallel LM kappa ', mu, &
                '  lambda_min ', wd(1), '  cond ', wd(size(wd))/wd(1)
        end do
        print '(a,3i5)', 'newton steps kappa 0, 1e-2, 1e-1 ', iterations
        print '(a,3es24.16)', 'energies ', energies
        print '(a,3es12.4)', 'd3D_exact ', errors
        ! Gates: the damped solves converge (solve stops otherwise) to a
        ! stationary point whose oracle error agrees with the baseline at
        ! the discretisation level.
        if (any(abs(errors(2:3) - errors(1)) > 0.05_dp*errors(1))) &
            error stop 'parallel damping end point'
    end subroutine
end program
