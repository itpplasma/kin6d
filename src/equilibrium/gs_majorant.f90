module kin6d_gs_majorant
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, gs_cell_data, apply_gs, &
                                  gs_linear_response, gs_recovery_t
    use kin6d_grad_shafranov_leaves, only: gs_source, gs_source_jvp
    use kin6d_triangle_basis, only: triangle_bernstein_controls
    use fortfem_equilibrated_flux_2d, only: equilibrated_flux_estimate, &
        equilibrated_flux_estimate_t, EQUILIBRATION_OK
    use fortfem_triangle_rt_arbitrary_order, only: triangle_rt_basis_t, &
        initialize_triangle_raviart_thomas, evaluate_triangle_raviart_thomas
    use kin6d_gs_error_leaves, only: gs_error_coefficient, gs_error_flux
    use fortfem_generated_equilibrated_flux_normal_component, only: &
        generated_equilibrated_flux_normal_component
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_sources
    use kin6d_gs_nonlinear, only: gs_forcing_t
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_krylov, only: real_conjugate_gradient_operator, KRYLOV_OK, &
        KRYLOV_MAX_ITERATIONS
    use kin6d_gs_csr, only: gs_csr_build, gs_csr_apply_full, csr_matrix_t
    use kin6d_gs_options, only: gs_use_direct
    use kin6d_gs_direct, only: gs_direct_t, gs_direct_available, gs_direct_factor, &
        gs_direct_solve, gs_direct_free
    use kin6d_gs_timers, only: gs_phase_tick
    use kin6d_gs_coupling, only: gs_coupling_t, gs_coupling_bem_nodes, &
        gs_coupled_apply, gs_coupling_steklov
    use kin6d_gs_coils, only: gs_circuit_set_t, gs_circuit_flux_vectors
    use kin6d_gs_free_boundary, only: gs_fb_system_t, gs_fb_eval_t, gs_fb_factor_t, &
        gs_fb_evaluate, gs_fb_factor, gs_fb_solve, gs_fb_release, GS_FB_OK
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gs_bem_curve_t, &
        gs_bem_panel_point, gs_bem_trace_node_index
    use fortfem_sparse_direct, only: sparse_direct_factor_t, &
        sparse_direct_factor_csc, sparse_direct_solve_factored, sparse_direct_free
    use fortnum_quadrature, only: gauss_legendre
    use kin6d_sort_keys, only: sort_key_order
    use, intrinsic :: iso_fortran_env, only: i64 => int64
    implicit none
    private
    public :: gs_functional_majorant
    ! Free-boundary FEM-BEM error estimates (DESIGN section 4, step P10).
    public :: gs_coupled_estimator_t, gs_coupled_majorant_t, gs_coupled_goal_t
    public :: gs_fb_majorant_t, gs_point_source_fn, gs_goal_weight_fn
    public :: gs_coupled_estimator_build, gs_coupled_estimator_free
    public :: gs_coupled_majorant, gs_coupled_goal, gs_coupled_ip_moments
    public :: gs_coupled_energy_apply, gs_coupled_reference_solve
    public :: gs_fb_majorant, gs_fb_inf_sup, gs_fb_goal, gs_coil_current_bounds

    integer, parameter, public :: GS_CM_OK = 0
    integer, parameter, public :: GS_CM_BAD_INPUT = 1
    integer, parameter, public :: GS_CM_BEM_FAILED = 2
    integer, parameter, public :: GS_CM_CELL_FAILED = 3
    integer, parameter, public :: GS_CM_SOLVE_FAILED = 4
    integer, parameter, public :: GS_CM_FB_FAILED = 5

    abstract interface
        !> Frozen source f of -div(grad psi / R) = f at a point of cell `cell`
        !! where the represented flux is psi.
        function gs_point_source_fn(cell, r, z, psi) result(f)
            import :: dp
            integer, intent(in) :: cell
            real(dp), intent(in) :: r, z, psi
            real(dp) :: f
        end function gs_point_source_fn
        !> Interior weight q of a goal Q(psi) = int q psi dA + <q_b, psi>.
        function gs_goal_weight_fn(r, z) result(q)
            import :: dp
            real(dp), intent(in) :: r, z
            real(dp) :: q
        end function gs_goal_weight_fn
    end interface

    !> Mesh data of the coupled estimates: a reference BEM on Gamma_c with
    !! `reference_panels` P2 panels per mesh edge (the stand-in for the exact
    !! Steklov-Poincare map S), the boundary moments of the recovered flux,
    !! the coupled Friedrichs constant and reference energy factor.
    type :: gs_coupled_estimator_t
        integer :: n = 0, nt = 0, nx = 0, nk = 0, ndof = 0, nel = 0
        integer :: reference_panels = 0
        type(gs_coupling_t) :: coupling
        type(gs_bem_curve_t) :: curve
        real(dp), allocatable :: s_ref(:, :)        !< (nx, nx) reference S
        real(dp), allocatable :: s_chol(:, :)       !< its Cholesky factor
        real(dp), allocatable :: restriction(:, :)  !< (nx, nt) FEM -> BEM trace
        real(dp), allocatable :: s_trace(:, :)      !< R^T S_ref R (nt, nt)
        real(dp), allocatable :: coil_moments(:, :) !< (nx, nk) <(S+d_nu)G_k, phi_j>
        type(gs_direct_t) :: energy_direct
        !> Sparse LU of the reference energy operator K + R^T S_ref R.
        type(sparse_direct_factor_t) :: energy_factor
        real(dp) :: cf2 = 0.0_dp                    !< coupled Friedrichs C_F^2
        integer :: cf_iterations = 0
        logical :: ready = .false.
    end type gs_coupled_estimator_t

    type :: gs_coupled_majorant_t
        real(dp) :: majorant = huge(1.0_dp)  !< D + oscillation + balance + N
        real(dp) :: flux_defect = 0.0_dp     !< D
        !> Source oscillation plus C_F ||div sigma + Pi_k f||.
        real(dp) :: balance_defect = 0.0_dp
        real(dp) :: neumann_defect = 0.0_dp  !< N
        ! A finite inverse norm has no control of untested boundary directions.
        character(32) :: exterior_status = 'unavailable-bem-remainder'
        integer :: iterations = 0
        logical :: converged = .false.
    end type gs_coupled_majorant_t

    type :: gs_coupled_goal_t
        real(dp) :: value = 0.0_dp       !< Q(psi_h)
        real(dp) :: correction = 0.0_dp  !< rho(psi_h; z_h)
        real(dp) :: bound = 0.0_dp       !< M_psi M_z (0: no bound available)
        real(dp) :: adjoint_majorant = 0.0_dp
        real(dp) :: lower = 0.0_dp, upper = 0.0_dp
    end type gs_coupled_goal_t

    type :: gs_fb_majorant_t
        type(gs_coupled_majorant_t) :: frozen  !< M0: frozen-source majorant
        real(dp) :: beta = 0.0_dp              !< inf-sup estimate of J
        real(dp) :: bound = huge(1.0_dp)       !< M0 / beta
        integer :: beta_iterations = 0
    end type gs_fb_majorant_t

    interface
        subroutine dpotrf(uplo, n, a, lda, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
        end subroutine dpotrf
        subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dpotrs
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface
contains
    ! ------------------------------------------------------------------
    ! Coupled FEM-BEM majorant (DESIGN section 4). Energy norm
    !   |||v|||^2 = int |grad v|^2 / R dA + <S v, v>_Gamma_c.
    ! For every sigma in H(div) and psi the solution of
    !   a(psi, v) + <S psi, v> = (f, v) + <g, v>,
    ! the error identity is
    !   a(e, v) + <S e, v> = (sigma - grad psi_h/R, grad v)
    !       + (f + div sigma, v) + <g - S psi_h - sigma.n, v>.
    ! FortFEM supplies a patchwise RT flux and the first two terms, with
    ! coupled C_F for the measured divergence defect and curved-cell source
    ! oscillation (local Payne--Weinberger only for affine cells). The last
    ! term is measured in the reference Steklov inverse norm below. Boundary
    ! RT moments are free recovery unknowns, never Dirichlet data for psi.
    ! Practical, not rigorous: S and its dual norm are those of a reference
    ! BEM with more panels, and C_F is a discrete eigenvalue estimate.
    ! ------------------------------------------------------------------

    subroutine gs_coupled_estimator_build(state, coupling, est, status, set, &
                                          reference_panels, cf2)
        !! Reference BEM, coil moments, C_F^2 and the reference energy factor.
        !! cf2, when given, replaces the eigenvalue estimate of C_F^2.
        type(gs_state_t), intent(in) :: state
        type(gs_coupling_t), intent(in) :: coupling
        type(gs_coupled_estimator_t), intent(inout) :: est
        integer, intent(out) :: status
        type(gs_circuit_set_t), intent(in), optional :: set
        integer, intent(in), optional :: reference_panels
        real(dp), intent(in), optional :: cf2
        real(dp), allocatable :: nodes(:, :), gx(:), gw(:), gamma(:, :)
        real(dp), allocatable :: xq(:, :), nq(:, :), wq(:), bq(:, :), gq(:, :)
        real(dp), allocatable :: gr(:, :), gz(:, :)
        real(dp), allocatable :: x(:), y(:), w(:)
        real(dp) :: point(2), tau(2), speed, basis(3)
        integer :: m, np, p, kk, j, l, k, info, stat, ng, point_id
        integer, allocatable :: bq_index(:, :)
        real(dp) :: t0
        logical :: direct_ok

        t0 = -1
        call gs_phase_tick('', t0)
        call gs_coupled_estimator_free(est)
        status = GS_CM_BAD_INPUT
        if (.not. state%ready .or. coupling%n_trace == 0) return
        if (state%element_order /= 2 .and. state%element_order /= 3) return
        m = 3
        if (present(reference_panels)) m = reference_panels
        if (m < 1) return
        est%coupling = coupling
        est%reference_panels = m
        est%n = size(state%points, 2)
        est%nt = coupling%n_trace
        est%ndof = size(state%triangles, 1)
        est%nel = size(state%triangles, 2)
        np = coupling%n_edges*m
        est%nx = 2*np

        status = GS_CM_BEM_FAILED
        call gs_coupling_bem_nodes(state, coupling, m, nodes, est%restriction)
        call gs_coupling_steklov(nodes, coupling%flux_degree, est%curve, est%s_ref, stat)
        if (stat /= 0) return
        call gs_phase_tick('est.steklov', t0)
        est%s_chol = est%s_ref
        call dpotrf('U', est%nx, est%s_chol, est%nx, info)
        if (info /= 0) return
        est%s_trace = matmul(transpose(est%restriction), &
                             matmul(est%s_ref, est%restriction))
        est%s_trace = 0.5_dp*(est%s_trace + transpose(est%s_trace))
        call gs_phase_tick('est.chol_strace', t0)

        ! Boundary quadrature on the reference panels for coil moments.
        ng = 8
        allocate (gx(ng), gw(ng))
        call gauss_legendre(ng, gx, gw)
        gx = 0.5_dp*(gx + 1.0_dp)
        gw = 0.5_dp*gw
        allocate (xq(2, np*ng), nq(2, np*ng), wq(np*ng), bq(3, np*ng), &
                  bq_index(3, np*ng))
        point_id = 0
        do p = 1, np
            do l = 1, ng
                point_id = point_id + 1
                call gs_bem_panel_point(est%curve, p, gx(l), point, tau, speed)
                basis = [(1.0_dp - gx(l))*(1.0_dp - 2.0_dp*gx(l)), &
                         4.0_dp*gx(l)*(1.0_dp - gx(l)), &
                         gx(l)*(2.0_dp*gx(l) - 1.0_dp)]
                xq(:, point_id) = point
                nq(:, point_id) = [tau(2), -tau(1)]
                wq(point_id) = gw(l)*speed
                bq(:, point_id) = basis
                do kk = 1, 3
                    bq_index(kk, point_id) = gs_bem_trace_node_index(est%curve, p, kk)
                end do
            end do
        end do
        est%nk = 0
        if (present(set)) then
            if (.not. allocated(set%incidence)) return
            est%nk = size(set%incidence, 2)
            allocate (gamma(est%nx, est%nk), gq(np*ng, est%nk), gr(np*ng, est%nk), &
                      gz(np*ng, est%nk))
            call gs_circuit_flux_vectors(set, nodes(1, :), nodes(2, :), gamma)
            call gs_circuit_flux_vectors(set, xq(1, :), xq(2, :), gq, gr, gz)
            est%coil_moments = matmul(est%s_ref, gamma)
            do point_id = 1, np*ng
                do k = 1, est%nk
                    do kk = 1, 3
                        j = bq_index(kk, point_id)
                        est%coil_moments(j, k) = est%coil_moments(j, k) &
                            + wq(point_id)*bq(kk, point_id) &
                            *(nq(1, point_id)*gr(point_id, k) &
                              + nq(2, point_id)*gz(point_id, k))/xq(1, point_id)
                    end do
                end do
            end do
        else
            allocate (est%coil_moments(est%nx, 0))
        end if

        call gs_phase_tick('est.boundary_maps', t0)
        status = GS_CM_SOLVE_FAILED
        call factor_energy(stat)
        if (stat /= 0) return
        call gs_phase_tick('est.factor_energy', t0)
        status = GS_CM_CELL_FAILED
        call triangle_duffy_quadrature(8, x, y, w, stat)
        if (stat /= 0) return
        if (present(cf2)) then
            est%cf2 = cf2
        else
            call friedrichs_estimate(stat)
            if (stat /= 0) return
        end if
        call gs_phase_tick('est.friedrichs', t0)
        est%ready = .true.
        status = GS_CM_OK
    contains
        subroutine factor_energy(fstat)
            !! CSC of K + scatter(S_trace) (symmetric: CSR of K is its CSC).
            integer, intent(out) :: fstat
            integer(i64), allocatable :: keys(:)
            real(dp), allocatable :: raw(:), values(:)
            integer, allocatable :: order(:), colptr(:), rowind(:)
            integer :: nk_, nnz, a_, b_, i_, p_, row, col
            integer(i64) :: previous
            nk_ = size(est%coupling%val) + est%nt**2
            allocate (keys(nk_), raw(nk_), order(nk_))
            nnz = 0
            do i_ = 1, est%n
                do p_ = est%coupling%row_ptr(i_), est%coupling%row_ptr(i_ + 1) - 1
                    nnz = nnz + 1
                    keys(nnz) = int(est%coupling%col(p_) - 1, i64)*est%n + i_ - 1
                    raw(nnz) = est%coupling%val(p_)
                end do
            end do
            do b_ = 1, est%nt
                do a_ = 1, est%nt
                    nnz = nnz + 1
                    keys(nnz) = int(est%coupling%trace_node(b_) - 1, i64)*est%n &
                                + est%coupling%trace_node(a_) - 1
                    raw(nnz) = est%s_trace(a_, b_)
                end do
            end do
            call sort_key_order(keys, order)
            allocate (rowind(nk_), values(nk_), colptr(est%n + 1))
            colptr = 0
            previous = -1
            nnz = 0
            do i_ = 1, nk_
                if (keys(order(i_)) /= previous) then
                    nnz = nnz + 1
                    previous = keys(order(i_))
                    col = int(previous/est%n) + 1
                    row = int(mod(previous, int(est%n, i64))) + 1
                    rowind(nnz) = row
                    values(nnz) = 0.0_dp
                    colptr(col + 1) = colptr(col + 1) + 1
                end if
                values(nnz) = values(nnz) + raw(order(i_))
            end do
            colptr(1) = 1
            do i_ = 1, est%n
                colptr(i_ + 1) = colptr(i_ + 1) + colptr(i_)
            end do
            if (gs_direct_available()) then
                call gs_direct_factor(est%energy_direct, est%n, colptr, &
                                      rowind(:nnz), values(:nnz), direct_ok)
                if (direct_ok) then
                    fstat = 0
                    return
                end if
                call gs_direct_free(est%energy_direct)
            end if
            call sparse_direct_factor_csc(est%energy_factor, est%n, colptr, &
                                          rowind(:nnz), values(:nnz), fstat)
        end subroutine factor_energy

        subroutine friedrichs_estimate(fstat)
            !! C_F^2 = max ||v||^2 / |||v|||^2: power iteration on A^-1 M with
            !! A = K + S_ref, M the L2 mass. The Rayleigh quotient approaches
            !! the discrete maximum from below (practical estimate).
            integer, intent(out) :: fstat
            type(csr_matrix_t) :: mass
            real(dp), allocatable :: mblocks(:, :, :), v(:), mv(:), av(:), z(:)
            real(dp) :: theta, previous, shape(10), grads(2, 10), hess(2, 2, 10)
            real(dp) :: pt(2), det
            integer :: e, qq, a, b, it, cg_stat, cstat

            fstat = 1
            allocate (mblocks(est%ndof, est%ndof, est%nel))
            mblocks = 0.0_dp
            fstat = 0
            !$omp parallel do default(shared) &
            !$omp private(e, qq, a, b, shape, grads, hess, pt, det, cstat) &
            !$omp schedule(static)
            do e = 1, est%nel
                do qq = 1, size(w)
                    call gs_cell_data(state, e, [x(qq), y(qq)], pt, det, &
                                      shape(:est%ndof), grads(:, :est%ndof), &
                                      hess(:, :, :est%ndof), cstat)
                    if (cstat /= 0) then
                        !$omp atomic write
                        fstat = cstat
                        exit
                    end if
                    do b = 1, est%ndof
                        do a = 1, est%ndof
                            mblocks(a, b, e) = mblocks(a, b, e) &
                                + w(qq)*det*shape(a)*shape(b)
                        end do
                    end do
                end do
            end do
            !$omp end parallel do
            if (fstat /= 0) return
            call gs_csr_build(mass, state%triangles, mblocks, est%n)
            allocate (v(est%n), mv(est%n), av(est%n), z(est%n))
            v = 1.0_dp
            theta = 0.0_dp
            fstat = 1
            do it = 1, 60
                call gs_csr_apply_full(mass, v, mv)
                call gs_coupled_reference_solve(est, mv, z, cg_stat)
                if (cg_stat /= 0) return
                call gs_coupled_energy_apply(est, z, av)
                previous = theta
                call gs_csr_apply_full(mass, z, mv)
                theta = dot_product(z, mv)/dot_product(z, av)
                v = z/norm2(z)
                est%cf_iterations = it
                if (abs(theta - previous) <= 1.0e-7_dp*theta) exit
            end do
            if (.not. (theta > 0.0_dp)) return
            est%cf2 = theta
            fstat = 0
        end subroutine friedrichs_estimate
    end subroutine gs_coupled_estimator_build

    subroutine gs_coupled_estimator_free(est)
        type(gs_coupled_estimator_t), intent(inout) :: est

        if (est%ready) then
            call sparse_direct_free(est%energy_factor)
            call gs_direct_free(est%energy_direct)
        end if
        est%ready = .false.
        est%cf2 = 0.0_dp
        if (allocated(est%s_ref)) deallocate (est%s_ref)
        if (allocated(est%s_chol)) deallocate (est%s_chol)
        if (allocated(est%restriction)) deallocate (est%restriction)
        if (allocated(est%s_trace)) deallocate (est%s_trace)
        if (allocated(est%coil_moments)) deallocate (est%coil_moments)
    end subroutine gs_coupled_estimator_free

    subroutine gs_coupled_energy_apply(est, input, output)
        !! output = (K + scatter(R^T S_ref R gather)) input: the reference
        !! coupled energy operator.
        type(gs_coupled_estimator_t), intent(in) :: est
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        integer :: i, p
        real(dp) :: acc

        !$omp parallel do default(shared) private(i, p, acc) schedule(static)
        do i = 1, est%n
            acc = 0.0_dp
            do p = est%coupling%row_ptr(i), est%coupling%row_ptr(i + 1) - 1
                acc = acc + est%coupling%val(p)*input(est%coupling%col(p))
            end do
            output(i) = acc
        end do
        !$omp end parallel do
        output(est%coupling%trace_node) = output(est%coupling%trace_node) &
            + matmul(est%s_trace, input(est%coupling%trace_node))
    end subroutine gs_coupled_energy_apply

    subroutine gs_coupled_reference_solve(est, rhs, solution, status)
        !! (K + R^T S_ref R) solution = rhs: sparse LU and two refinement steps.
        type(gs_coupled_estimator_t), intent(inout) :: est
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: solution(:)
        integer, intent(out) :: status
        real(dp), allocatable :: residual(:), correction(:)
        integer :: sweep

        allocate (residual(est%n), correction(est%n))
        call energy_solve(rhs, solution, status)
        do sweep = 1, 2
            if (status /= 0) return
            call gs_coupled_energy_apply(est, solution, residual)
            residual = rhs - residual
            call energy_solve(residual, correction, status)
            solution = solution + correction
        end do

    contains
        subroutine energy_solve(b, x, sstat)
            real(dp), intent(in) :: b(:)
            real(dp), intent(out) :: x(:)
            integer, intent(out) :: sstat
            logical :: ok
            if (est%energy_direct%factored) then
                call gs_direct_solve(est%energy_direct, b, x, ok)
                sstat = merge(0, 1, ok)
            else
                call sparse_direct_solve_factored(est%energy_factor, b, x, sstat)
            end if
        end subroutine energy_solve
    end subroutine gs_coupled_reference_solve

    function gs_coupled_ip_moments(est, mu0) result(moments)
        !! Boundary moments of Ip(psi) = (<S psi, 1> - <g, 1>)/mu0 (Ampere on
        !! Gamma_c, the conormal flux being g - S psi): S_ref 1 / mu0. The
        !! constant -<g, 1>/mu0 is sum(boundary load)/mu0 with a minus sign.
        type(gs_coupled_estimator_t), intent(in) :: est
        real(dp), intent(in) :: mu0
        real(dp) :: moments(est%nx)

        moments = sum(est%s_ref, dim=2)/mu0
    end function gs_coupled_ip_moments

    subroutine gs_coupled_majorant(state, est, psi, source, boundary_load, res, &
                                   status)
        !! Patchwise H(div) flux with the physical FEM--BEM interface residual.
        !! Boundary patch moments are free in the RT recovery. This does not
        !! impose a Dirichlet condition on psi: its arbitrary trace enters the
        !! separately measured ||sigma.n + S_ref psi - g||_{S_ref^-1} term.
        !! Both source balance and curved-cell oscillation use the COUPLED
        !! Friedrichs constant. Reference BEM, that constant and quadrature
        !! remain approximate, so this is a practical energy estimate.
        type(gs_state_t), intent(in) :: state
        type(gs_coupled_estimator_t), intent(inout) :: est
        real(dp), intent(in) :: psi(:), boundary_load(:)
        procedure(gs_point_source_fn) :: source
        type(gs_coupled_majorant_t), intent(out) :: res
        integer, intent(out) :: status
        type(equilibrated_flux_estimate_t) :: recovered
        type(triangle_rt_basis_t) :: rt
        real(dp), allocatable :: x(:), y(:), w(:), coefficient(:, :), f(:, :)
        real(dp), allocatable :: flux(:, :, :), lower(:), geometry(:, :, :)
        real(dp), allocatable :: mu(:), smu(:), gx(:), gw(:), rt_value(:, :)
        real(dp), allocatable :: rt_div(:)
        real(dp) :: point(2), det, shape(10), grads(2, 10), hess(2, 2, 10)
        real(dp) :: vals(10), controls(10), grad(2), xi(2), normal(2), sigma(2)
        real(dp) :: t, basis(3), normal_flux
        integer :: e, q, nd, stat, failed, p, i, panel, side, j, node, ng

        status = GS_CM_BAD_INPUT
        if (.not. est%ready) return
        if (size(psi) /= est%n .or. size(boundary_load) /= est%nx) return
        if (.not. all(ieee_is_finite(psi))) return
        if (.not. all(ieee_is_finite(boundary_load))) return
        if (.not. (est%cf2 > 0.0_dp)) return
        nd = est%ndof
        status = GS_CM_CELL_FAILED
        call triangle_duffy_quadrature(2*state%element_order + 6, x, y, w, stat)
        if (stat /= 0) return
        allocate (coefficient(size(w), est%nel), f(size(w), est%nel), &
            flux(2, size(w), est%nel), lower(est%nel), geometry(2, nd, est%nel))
        failed = 0
        !$omp parallel do schedule(static) reduction(max:failed) &
        !$omp private(q, point, det, shape, grads, hess, vals, controls, grad, stat)
        do e = 1, est%nel
            geometry(:, :, e) = state%points(:, state%triangles(:, e))
            call triangle_bernstein_controls(geometry(1, :, e), controls(:nd))
            if (minval(controls(:nd)) <= 0.0_dp) then
                failed = 1
                cycle
            end if
            call gs_error_coefficient(maxval(controls(:nd)), lower(e))
            vals(:nd) = psi(state%triangles(:, e))
            do q = 1, size(w)
                call gs_cell_data(state, e, [x(q), y(q)], point, det, &
                    shape(:nd), grads(:, :nd), hess(:, :, :nd), stat)
                if (stat /= 0) then
                    failed = 1
                    exit
                end if
                call gs_error_coefficient(point(1), coefficient(q, e))
                f(q, e) = source(e, point(1), point(2), &
                    dot_product(shape(:nd), vals(:nd)))
                grad = matmul(grads(:, :nd), vals(:nd) - vals(1))
                call gs_error_flux(point(1), grad(1), grad(2), &
                    flux(1, q, e), flux(2, q, e))
            end do
        end do
        !$omp end parallel do
        if (failed /= 0) return
        status = GS_CM_SOLVE_FAILED
        call equilibrated_flux_estimate(state%points, state%triangles(1:3, :), &
            state%element_order, x, y, w, coefficient, f, flux, lower, &
            sqrt(est%cf2), recovered, stat, geometry_nodes=geometry)
        if (stat /= EQUILIBRATION_OK) return

        ! Piola preserves sigma.n ds = sigma_hat.n_hat dt exactly, including
        ! curved edges. Integrate these RT normal moments against the SAME
        ! reference BEM trace basis used by S_ref and boundary_load.
        call initialize_triangle_raviart_thomas(state%element_order, rt, stat)
        if (stat /= 0) return
        allocate (rt_value(2, size(recovered%flux_coefficients, 1)), &
            rt_div(size(recovered%flux_coefficients, 1)))
        mu = matmul(est%s_ref, matmul(est%restriction, &
            psi(est%coupling%trace_node))) - boundary_load
        ng = 8
        allocate (gx(ng), gw(ng))
        call gauss_legendre(ng, gx, gw)
        gx = (gx + 1.0_dp)/2
        gw = gw/2
        do p = 1, est%coupling%n_edges*est%reference_panels
            i = (p - 1)/est%reference_panels + 1
            panel = p - (i - 1)*est%reference_panels
            e = est%coupling%edge_element(i)
            side = est%coupling%edge_side(i)
            do q = 1, ng
                t = (real(panel - 1, dp) + gx(q))/est%reference_panels
                select case (side)
                case (1)
                    xi = [t, 0.0_dp]
                    normal = [0.0_dp, -1.0_dp]
                case (2)
                    xi = [1.0_dp - t, t]
                    normal = [1.0_dp, 1.0_dp]
                case (3)
                    xi = [0.0_dp, 1.0_dp - t]
                    normal = [-1.0_dp, 0.0_dp]
                case default
                    return
                end select
                call evaluate_triangle_raviart_thomas(rt, xi(1), xi(2), &
                    rt_value, rt_div, stat)
                if (stat /= 0) return
                sigma = matmul(rt_value, recovered%flux_coefficients(:, e))
                call generated_equilibrated_flux_normal_component(sigma, normal, &
                    normal_flux)
                basis = [(1.0_dp - gx(q))*(1.0_dp - 2.0_dp*gx(q)), &
                    4.0_dp*gx(q)*(1.0_dp - gx(q)), gx(q)*(2.0_dp*gx(q) - 1.0_dp)]
                do j = 1, 3
                    node = gs_bem_trace_node_index(est%curve, p, j)
                    mu(node) = mu(node) + gw(q)/est%reference_panels &
                        *basis(j)*normal_flux
                end do
            end do
        end do
        smu = mu
        call dpotrs('U', est%nx, 1, est%s_chol, est%nx, smu, est%nx, stat)
        if (stat /= 0) return
        res%flux_defect = recovered%flux_defect
        res%balance_defect = recovered%oscillation + recovered%divergence_defect
        res%neumann_defect = sqrt(max(0.0_dp, dot_product(mu, smu)))
        res%majorant = res%flux_defect + res%balance_defect + res%neumann_defect
        res%iterations = 0
        res%converged = .true.
        if (.not. ieee_is_finite(res%majorant)) return
        status = GS_CM_OK
    end subroutine gs_coupled_majorant

    subroutine goal_loads(state, est, psi, source, weight, loads, value, status)
        !! loads(:, 1) = int f phi_i, loads(:, 2) = int q phi_i, value = int q psi.
        type(gs_state_t), intent(in) :: state
        type(gs_coupled_estimator_t), intent(in) :: est
        real(dp), intent(in) :: psi(:)
        procedure(gs_point_source_fn) :: source
        procedure(gs_goal_weight_fn) :: weight
        real(dp), allocatable, intent(out) :: loads(:, :)
        real(dp), intent(out) :: value
        integer, intent(out) :: status
        real(dp), allocatable :: x(:), y(:), w(:)
        real(dp) :: shape(10), grads(2, 10), hess(2, 2, 10), pt(2), det, vals(10)
        real(dp) :: f, q, ph
        integer :: e, qq, a, nd

        call triangle_duffy_quadrature(14, x, y, w, status)
        if (status /= 0) return
        nd = est%ndof
        allocate (loads(est%n, 2))
        loads = 0.0_dp
        value = 0.0_dp
        do e = 1, est%nel
            vals(:nd) = psi(state%triangles(:, e))
            do qq = 1, size(w)
                call gs_cell_data(state, e, [x(qq), y(qq)], pt, det, shape(:nd), &
                                  grads(:, :nd), hess(:, :, :nd), status)
                if (status /= 0) return
                ph = dot_product(shape(:nd), vals(:nd))
                f = source(e, pt(1), pt(2), ph)
                q = weight(pt(1), pt(2))
                value = value + w(qq)*det*q*ph
                do a = 1, nd
                    loads(state%triangles(a, e), 1) = loads(state%triangles(a, e), 1) &
                        + w(qq)*det*f*shape(a)
                    loads(state%triangles(a, e), 2) = loads(state%triangles(a, e), 2) &
                        + w(qq)*det*q*shape(a)
                end do
            end do
        end do
    end subroutine goal_loads

    subroutine gs_coupled_goal(state, est, psi, source, boundary_load, weight, &
                               goal_boundary, primal_majorant, res, status)
        !! Goal Q(psi) = int q psi dA + <q_b, psi> (q_b with reference moments
        !! goal_boundary) of the frozen coupled problem. With z_h the discrete
        !! adjoint (K + S_ref) z_h = Q', Q(psi) - Q(psi_h) = r(z_h) + r(z - z_h)
        !! and |r(z - z_h)| <= |||psi - psi_h||| |||z - z_h||| <= M_psi M_z, so
        !! Q(psi) lies in Q(psi_h) + rho +- M_psi M_z (rho = r(z_h)).
        type(gs_state_t), intent(in) :: state
        type(gs_coupled_estimator_t), intent(inout) :: est
        real(dp), intent(in) :: psi(:), boundary_load(:), goal_boundary(:)
        procedure(gs_point_source_fn) :: source
        procedure(gs_goal_weight_fn) :: weight
        real(dp), intent(in) :: primal_majorant
        type(gs_coupled_goal_t), intent(out) :: res
        integer, intent(out) :: status
        type(gs_coupled_majorant_t) :: adjoint
        real(dp), allocatable :: loads(:, :), rhs(:), z(:), apsi(:), bnd(:)
        real(dp) :: value
        integer :: stat

        status = GS_CM_BAD_INPUT
        if (.not. est%ready .or. size(goal_boundary) /= est%nx) return
        if (size(boundary_load) /= est%nx .or. size(psi) /= est%n) return
        status = GS_CM_CELL_FAILED
        call goal_loads(state, est, psi, source, weight, loads, value, stat)
        if (stat /= 0) return
        allocate (bnd(est%n), rhs(est%n), z(est%n), apsi(est%n))
        bnd = 0.0_dp
        bnd(est%coupling%trace_node) = matmul(transpose(est%restriction), &
                                              goal_boundary)
        rhs = loads(:, 2) + bnd
        status = GS_CM_SOLVE_FAILED
        call gs_coupled_reference_solve(est, rhs, z, stat)
        if (stat /= 0) return
        res%value = value + dot_product(goal_boundary, &
            matmul(est%restriction, psi(est%coupling%trace_node)))
        bnd = 0.0_dp
        bnd(est%coupling%trace_node) = matmul(transpose(est%restriction), &
                                              boundary_load)
        call gs_coupled_energy_apply(est, psi, apsi)
        res%correction = dot_product(z, loads(:, 1) + bnd - apsi)
        call gs_coupled_majorant(state, est, z, adjoint_source, goal_boundary, &
                                 adjoint, stat)
        if (stat /= GS_CM_OK) then
            status = stat
            return
        end if
        res%adjoint_majorant = adjoint%majorant
        res%bound = primal_majorant*adjoint%majorant
        res%lower = res%value + res%correction - res%bound
        res%upper = res%value + res%correction + res%bound
        status = GS_CM_OK
    contains
        function adjoint_source(cell, r, zz, p) result(f)
            integer, intent(in) :: cell
            real(dp), intent(in) :: r, zz, p
            real(dp) :: f
            f = weight(r, zz) + 0.0_dp*(cell + p)
        end function adjoint_source
    end subroutine gs_coupled_goal

    ! ------------------------------------------------------------------
    ! Nonlinear free boundary (DESIGN section 4, "Nonlinear"). M0 is the
    ! frozen-source majorant at psi_h; with fixed coil currents J can be
    ! indefinite, so the coercivity of the fixed-boundary path is replaced
    ! by the inf-sup constant beta of J in the energy pairing, and the
    ! practical bound is |||e||| <~ M0 / beta (first order in the error).
    ! ------------------------------------------------------------------

    subroutine gs_fb_majorant(state, sys, est, currents, res, status, forcing, &
                              max_iterations)
        !! Frozen-source majorant of the free-boundary solution state%psi and
        !! the inf-sup estimate. The frozen source is the law inside the
        !! plasma (plasma cells, s < 1) plus forcing(R, Z), which must be the
        !! pointwise form of sys%extra_load (zero when that load is zero).
        type(gs_state_t), intent(inout) :: state
        type(gs_fb_system_t), intent(in) :: sys
        type(gs_coupled_estimator_t), intent(inout) :: est
        real(dp), intent(in) :: currents(:)
        type(gs_fb_majorant_t), intent(out) :: res
        integer, intent(out) :: status
        procedure(gs_forcing_t), optional :: forcing
        integer, intent(in), optional :: max_iterations
        type(gs_fb_eval_t) :: ev
        type(gs_profile_law_t) :: law
        real(dp), allocatable :: psi(:)
        real(dp) :: sigma
        integer :: stat

        status = GS_CM_BAD_INPUT
        if (.not. est%ready .or. size(currents) /= est%nk) return
        psi = state%psi
        status = GS_CM_FB_FAILED
        call gs_fb_evaluate(state, sys, psi, currents, ev, stat, jacobian=.true.)
        if (stat /= GS_FB_OK) return
        law = sys%law
        law%axis = ev%psi_axis - ev%psi_bound
        sigma = real(sys%cfg%axis_sign, dp)
        call gs_coupled_majorant(state, est, psi, frozen_source, &
                                 matmul(est%coil_moments, currents), res%frozen, stat)
        if (stat /= GS_CM_OK) then
            status = stat
            return
        end if
        call gs_fb_inf_sup(state, sys, ev, res%beta, res%beta_iterations, stat, &
                           max_iterations)
        if (stat /= GS_CM_OK) then
            status = stat
            return
        end if
        res%bound = res%frozen%majorant/res%beta
        status = GS_CM_OK
    contains
        function frozen_source(cell, r, z, p) result(f)
            integer, intent(in) :: cell
            real(dp), intent(in) :: r, z, p
            real(dp) :: f, pp, ff
            f = 0.0_dp
            if (ev%plasma_cell(cell)) then
                if (sigma*(p - ev%psi_bound) > 0.0_dp) then
                    call gs_law_sources(law, p - ev%psi_bound, pp, ff)
                    call gs_source(r, sys%mu0, pp, ff, f)
                end if
            end if
            if (present(forcing)) f = f + forcing(r, z)
        end function frozen_source
    end subroutine gs_fb_majorant

    subroutine gs_fb_inf_sup(state, sys, ev, beta, iterations, status, &
                             max_iterations)
        !! beta = min_v sup_w <J v, w> / (|||v||| |||w|||), the smallest
        !! singular value of J in the energy pairing of A = K + S_h. Power
        !! iteration on T = J^-1 A J^-T A, self-adjoint and positive in the A
        !! inner product, whose largest eigenvalue is 1/beta^2. The Rayleigh
        !! quotient approaches it from below, so beta is a practical
        !! (slightly high) estimate, not a bound.
        type(gs_state_t), intent(in) :: state
        type(gs_fb_system_t), intent(in) :: sys
        type(gs_fb_eval_t), intent(in) :: ev
        real(dp), intent(out) :: beta
        integer, intent(out) :: iterations, status
        integer, intent(in), optional :: max_iterations
        type(gs_fb_factor_t) :: fac
        real(dp), allocatable :: v(:), av(:), t(:), u(:), wv(:)
        real(dp) :: theta, previous, energy
        integer :: it, stat, n, cap

        status = GS_CM_FB_FAILED
        beta = 0.0_dp
        iterations = 0
        cap = 200
        if (present(max_iterations)) cap = max_iterations
        n = sys%n
        call gs_fb_factor(state, sys, ev, fac, stat)
        if (stat /= GS_FB_OK) return
        allocate (v(n), av(n), t(n), u(n), wv(n))
        v = 1.0_dp + state%points(1, :) - 0.3_dp*state%points(2, :)
        call gs_coupled_apply(sys%coupling, v, av)
        v = v/sqrt(dot_product(v, av))
        theta = 0.0_dp
        do it = 1, cap
            call gs_coupled_apply(sys%coupling, v, av)
            call gs_fb_solve(sys, ev, fac, av, t, stat, transpose=.true.)
            if (stat /= GS_FB_OK) exit
            call gs_coupled_apply(sys%coupling, t, u)
            call gs_fb_solve(sys, ev, fac, u, wv, stat)
            if (stat /= GS_FB_OK) exit
            previous = theta
            theta = dot_product(av, wv)
            call gs_coupled_apply(sys%coupling, wv, u)
            energy = dot_product(wv, u)
            if (.not. (energy > 0.0_dp)) exit
            v = wv/sqrt(energy)
            iterations = it
            if (abs(theta - previous) <= 1.0e-6_dp*theta) exit
        end do
        call gs_fb_release(fac)
        if (stat /= GS_FB_OK .or. .not. (theta > 0.0_dp)) return
        beta = 1.0_dp/sqrt(theta)
        status = GS_CM_OK
    end subroutine gs_fb_inf_sup

    subroutine gs_fb_goal(state, sys, est, currents, weight, goal_boundary, &
                          residual_majorant, res, status)
        !! Goal of the nonlinear free-boundary solution. With the discrete
        !! adjoint J^T z_h = Q', to first order Q(psi) - Q(psi_h) = <z, r>, r
        !! the residual functional, |<z, r>| <= |||z||| ||r||_* <= |||z_h||| M0
        !! (M0 = residual_majorant, the frozen majorant, bounds ||r||_*):
        !! res%bound = |||z_h||| M0, an estimate (z ~ z_h), not a guaranteed
        !! bound. rho = -z_h . F_ref(psi_h) (reference S) is the computable part,
        !! small by Galerkin orthogonality; no remainder bound is claimed for it
        !! because the adjoint residual carries the interface and point-
        !! evaluation (axis, boundary flux) terms, which are not in H^-1.
        type(gs_state_t), intent(inout) :: state
        type(gs_fb_system_t), intent(in) :: sys
        type(gs_coupled_estimator_t), intent(in) :: est
        real(dp), intent(in) :: currents(:), goal_boundary(:), residual_majorant
        procedure(gs_goal_weight_fn) :: weight
        type(gs_coupled_goal_t), intent(out) :: res
        integer, intent(out) :: status
        type(gs_fb_eval_t) :: ev
        type(gs_fb_factor_t) :: fac
        real(dp), allocatable :: loads(:, :), psi(:), rhs(:), z(:), apsi(:), bnd(:)
        real(dp), allocatable :: plasma_load(:)
        real(dp) :: value
        integer :: stat

        status = GS_CM_BAD_INPUT
        if (.not. est%ready .or. size(goal_boundary) /= est%nx) return
        if (size(currents) /= est%nk) return
        psi = state%psi
        status = GS_CM_FB_FAILED
        call fb_evaluate_loads(state, sys, psi, currents, ev, plasma_load, stat)
        if (stat /= GS_FB_OK) return
        call gs_fb_factor(state, sys, ev, fac, stat)
        if (stat /= GS_FB_OK) return
        status = GS_CM_CELL_FAILED
        call goal_loads(state, est, psi, zero_source, weight, loads, value, stat)
        if (stat /= 0) return
        allocate (bnd(est%n), rhs(est%n), z(est%n), apsi(est%n))
        bnd = 0.0_dp
        bnd(est%coupling%trace_node) = matmul(transpose(est%restriction), &
                                              goal_boundary)
        rhs = loads(:, 2) + bnd
        status = GS_CM_FB_FAILED
        call gs_fb_solve(sys, ev, fac, rhs, z, stat, transpose=.true.)
        call gs_fb_release(fac)
        if (stat /= GS_FB_OK) return
        res%value = value + dot_product(goal_boundary, &
            matmul(est%restriction, psi(est%coupling%trace_node)))
        bnd = 0.0_dp
        bnd(est%coupling%trace_node) = matmul(transpose(est%restriction), &
                                              matmul(est%coil_moments, currents))
        call gs_coupled_energy_apply(est, psi, apsi)
        res%correction = dot_product(z, plasma_load + bnd - apsi)
        call gs_coupled_apply(sys%coupling, z, apsi)
        res%adjoint_majorant = sqrt(max(0.0_dp, dot_product(z, apsi)))
        res%bound = res%adjoint_majorant*residual_majorant
        res%lower = res%value - res%bound
        res%upper = res%value + res%bound
        status = GS_CM_OK
    contains
        function zero_source(cell, r, zz, p) result(f)
            integer, intent(in) :: cell
            real(dp), intent(in) :: r, zz, p
            real(dp) :: f
            f = 0.0_dp*(cell + r + zz + p)
        end function zero_source
    end subroutine gs_fb_goal

    subroutine fb_evaluate_loads(state, sys, psi, currents, ev, load, status)
        !! Evaluation with Jacobian and load = extra + b(psi) = A_h psi - F(psi)
        !! - sum_k I_k load_k.
        type(gs_state_t), intent(inout) :: state
        type(gs_fb_system_t), intent(in) :: sys
        real(dp), intent(in) :: psi(:), currents(:)
        type(gs_fb_eval_t), intent(inout) :: ev
        real(dp), allocatable, intent(out) :: load(:)
        integer, intent(out) :: status

        call gs_fb_evaluate(state, sys, psi, currents, ev, status, jacobian=.true.)
        if (status /= GS_FB_OK) return
        allocate (load(size(psi)))
        call gs_coupled_apply(sys%coupling, psi, load)
        load = load - ev%residual - matmul(sys%coil_load, currents)
    end subroutine fb_evaluate_loads

    subroutine gs_coil_current_bounds(d, lambda, scale, free, row_bounds, bounds, &
                                      status, sensitivity)
        !! Coil-current goal (DESIGN 4): at a Gauss-Newton optimum with scaled
        !! rows (W_r = I), dI = -H^-1 D^T dr on the free set, H = D^T D +
        !! lambda^2 diag(1/scale^2). |dI_k| <= sum_i |(H^-1 D^T)_ki| b_i for
        !! |dr_i| <= b_i, attained for dr_i = b_i sign((H^-1 D^T)_ki). Fixed
        !! circuits (free = .false.) get 0.
        real(dp), intent(in) :: d(:, :), lambda, scale(:), row_bounds(:)
        logical, intent(in) :: free(:)
        real(dp), intent(out) :: bounds(:)
        integer, intent(out) :: status
        real(dp), intent(out), optional :: sensitivity(:, :)
        real(dp), allocatable :: h(:, :), rhs(:, :), df(:, :)
        integer, allocatable :: idx(:), ipiv(:)
        integer :: nf, info, k

        status = GS_CM_BAD_INPUT
        bounds = 0.0_dp
        if (size(d, 2) /= size(scale) .or. size(free) /= size(scale)) return
        if (size(row_bounds) /= size(d, 1) .or. size(bounds) /= size(scale)) return
        if (any(row_bounds < 0.0_dp) .or. any(scale <= 0.0_dp)) return
        idx = pack([(k, k=1, size(scale))], free)
        nf = size(idx)
        if (present(sensitivity)) sensitivity = 0.0_dp
        status = GS_CM_OK
        if (nf == 0) return
        df = d(:, idx)
        h = matmul(transpose(df), df)
        do k = 1, nf
            h(k, k) = h(k, k) + (lambda/scale(idx(k)))**2
        end do
        rhs = transpose(df)
        allocate (ipiv(nf))
        call dgesv(nf, size(d, 1), h, nf, ipiv, rhs, nf, info)
        if (info /= 0) then
            status = GS_CM_SOLVE_FAILED
            return
        end if
        bounds(idx) = matmul(abs(rhs), row_bounds)
        if (present(sensitivity)) sensitivity(idx, :) = rhs
    end subroutine gs_coil_current_bounds

    subroutine gs_functional_majorant(state, mu0, pprime, ffprime, majorant, &
                                     flux_defect, balance_defect, status, &
                                     recovery_converged, recovery_iterations, law, &
                                     stability_lower_bound, forcing, &
                                     estimate_stability, stability_estimate)
        ! Repin functional majorant for -div(grad(psi)/R)=source, homogeneous
        ! Dirichlet trace, on the represented domain. Recover y in continuous
        ! vector Pk, k=2 or 3 (a subspace of H(div)) by minimizing D^2+C_F^2*B^2.
        ! Returns D+C_F*B, with ordinary quadrature; no rigorous enclosure.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: mu0, pprime, ffprime
        real(dp), intent(out) :: majorant, flux_defect, balance_defect
        integer, intent(out) :: status
        logical, intent(out), optional :: recovery_converged
        integer, intent(out), optional :: recovery_iterations
        type(gs_profile_law_t), intent(in), optional :: law
        real(dp), intent(in), optional :: stability_lower_bound
        procedure(gs_forcing_t), optional :: forcing
        logical, intent(in), optional :: estimate_stability
        real(dp), intent(out), optional :: stability_estimate
        real(dp), allocatable :: x(:), y(:), w(:)
        real(dp), allocatable :: rhs(:), recovered(:), rhs_local(:, :)
        real(dp), allocatable :: source_mass(:, :, :)
        type(gs_recovery_t), pointer :: rec
        type(gs_recovery_t), target :: local_rec
        real(dp) :: cf2, lo(2), hi(2), residual, alpha, signature
        real(dp) :: flux_sum, balance_sum
        logical :: practical, solved, build, cached
        integer :: e, j, n, iterations, ndof, nel, failed_element, failed_status
        majorant = huge(1.0_dp); flux_defect = 0; balance_defect = 0
        status = 1
        practical = .false.
        if (present(estimate_stability)) practical = estimate_stability
        if (.not. present(law)) practical = .false.
        if (present(law)) then
            if (law%normalized) practical = .false.
        end if
        if (present(stability_lower_bound)) practical = .false.
        if (present(stability_estimate)) stability_estimate = 0
        if (present(recovery_converged)) recovery_converged = .false.
        if (present(recovery_iterations)) recovery_iterations = 0
        if (.not. state%ready) return
        if (state%element_order /= 2 .and. state%element_order /= 3) return
        ndof = size(state%triangles, 1)
        if (practical) then
            allocate (source_mass(ndof, ndof, size(state%triangles, 2)))
            source_mass = 0
        end if
        nel = size(state%triangles, 2)
        if (any(abs(pack(state%psi, state%boundary)) > tiny(1.0_dp))) return
        if (mu0 <= 0) return
        if (.not. all(ieee_is_finite([mu0, pprime, ffprime]))) return
        cached = associated(state%recovery)
        if (cached) then
            rec => state%recovery
        else
            rec => local_rec
        end if
        signature = sum(state%points) + sum(state%area)
        build = .true.
        if (rec%ready) build = rec%nodes /= size(state%psi) .or. rec%cells /= nel &
            .or. rec%order /= state%element_order &
            .or. abs(rec%signature - signature) > 0
        if (build) then
            ! Mesh-only part: the Poincare constant of the enclosing box.
            lo = huge(1.0_dp); hi = -huge(1.0_dp)
            !$omp parallel default(shared)
            block
                real(dp) :: bezier(2, ndof), block_lo(2), block_hi(2)
                integer :: ie, jj
                block_lo = huge(1.0_dp); block_hi = -huge(1.0_dp)
                !$omp do schedule(static)
                do ie = 1, nel
                    bezier = state%points(:, state%triangles(:, ie))
                    do jj = 1, 2
                        call triangle_bernstein_controls( &
                            state%points(jj, state%triangles(:, ie)), bezier(jj, :))
                    end do
                    block_lo = min(block_lo, minval(bezier, dim=2))
                    block_hi = max(block_hi, maxval(bezier, dim=2))
                end do
                !$omp end do
                !$omp critical(gs_majorant_box)
                lo = min(lo, block_lo)
                hi = max(hi, block_hi)
                !$omp end critical(gs_majorant_box)
            end block
            !$omp end parallel
            if (lo(1) <= 0 .or. any(hi <= lo)) return
            ! Zero extension to the enclosing rectangle gives this Poincare bound.
            cf2 = hi(1)/(acos(-1.0_dp)**2*sum(1/(hi - lo)**2))
        else
            cf2 = rec%cf2
        end if
        n = 2*size(state%psi)
        if (build) then
            rec%ready = .false.
            call gs_direct_free(rec%factor)
            if (allocated(rec%matrix)) deallocate (rec%matrix)
            if (allocated(rec%indices)) deallocate (rec%indices)
            if (allocated(rec%diagonal)) deallocate (rec%diagonal)
            allocate (rec%matrix(2*ndof, 2*ndof, nel), rec%indices(2*ndof, nel), &
                      rec%diagonal(n))
            rec%matrix = 0; rec%diagonal = 0
        end if
        allocate (rhs(n), recovered(n), rhs_local(2*ndof, nel))
        rhs = 0; recovered = 0; rhs_local = 0
        call triangle_duffy_quadrature(8, x, y, w, status)
        if (status /= 0) return
        failed_element = huge(1); failed_status = 1
        !$omp parallel default(shared)
        block
            real(dp) :: shape(ndof), gradients(2, ndof), hessians(2, 2, ndof)
            real(dp) :: d(ndof), point(2), det, grad(2), source, weight, pp, ff, psi
            real(dp) :: dpp, dff, derivative
            integer :: ie, q, jj, kk, aa, bb, row, col, stat, snapshot
            !$omp do schedule(static)
            do ie = 1, nel
                !$omp atomic read
                snapshot = failed_element
                if (ie > snapshot) cycle
                if (build) then
                    do jj = 1, ndof
                        rec%indices(2*jj - 1:2*jj, ie) = &
                            [2*state%triangles(jj, ie) - 1, 2*state%triangles(jj, ie)]
                    end do
                end if
                d = state%psi(state%triangles(:, ie)) &
                    - state%psi(state%triangles(1, ie))
                do q = 1, size(w)
                    call gs_cell_data(state, ie, [x(q), y(q)], point, det, shape, &
                        gradients, hessians, stat)
                    if (stat /= 0) then
                        call record_failure(ie, stat)
                        exit
                    end if
                    weight = w(q)*det
                    grad = matmul(gradients, d)
                    pp = pprime; ff = ffprime
                    if (present(law)) then
                        psi = state%psi(state%triangles(1, ie)) + dot_product(shape, d)
                        if (practical) then
                            call gs_law_sources(law, psi, pp, ff, dpp, dff)
                            call gs_source_jvp(point(1), mu0, pp, ff, &
                                0.0_dp, 0.0_dp, dpp, dff, derivative)
                            do kk = 1, ndof
                                do jj = 1, ndof
                                    source_mass(jj, kk, ie) = source_mass(jj, kk, ie) &
                                        + weight*max(0.0_dp, derivative) &
                                        *shape(jj)*shape(kk)
                                end do
                            end do
                        else
                            call gs_law_sources(law, psi, pp, ff)
                        end if
                    end if
                    call gs_source(point(1), mu0, pp, ff, source)
                    if (present(forcing)) source = source + forcing(point(1), point(2))
                    do jj = 1, ndof
                        do aa = 1, 2
                            row = 2*(jj - 1) + aa
                            rhs_local(row, ie) = rhs_local(row, ie) + weight &
                                *(shape(jj)*grad(aa) - cf2*source*gradients(aa, jj))
                            if (.not. build) cycle
                            do kk = 1, ndof
                                do bb = 1, 2
                                    col = 2*(kk - 1) + bb
                                    rec%matrix(row, col, ie) = &
                                        rec%matrix(row, col, ie) &
                                        + weight*cf2*gradients(aa, jj)*gradients(bb, kk)
                                    if (aa == bb) rec%matrix(row, col, ie) = &
                                        rec%matrix(row, col, ie) &
                                        + weight*point(1)*shape(jj)*shape(kk)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$omp end do
        end block
        !$omp end parallel
        if (failed_element /= huge(1)) then
            status = failed_status
            return
        end if
        ! Ordered gather: the sums do not depend on the thread count.
        do e = 1, nel
            do j = 1, 2*ndof
                rhs(rec%indices(j, e)) = rhs(rec%indices(j, e)) + rhs_local(j, e)
            end do
            if (build) then
                do j = 1, 2*ndof
                    rec%diagonal(rec%indices(j, e)) = rec%diagonal(rec%indices(j, e)) &
                        + rec%matrix(j, j, e)
                end do
            end if
        end do
        if (build) then
            call gs_csr_build(rec%system, rec%indices, rec%matrix, n)
            rec%nodes = size(state%psi); rec%cells = nel
            rec%order = state%element_order
            rec%signature = signature; rec%cf2 = cf2
            rec%ready = .true.
        end if
        if (present(recovery_converged)) recovery_converged = .true.
        if (norm2(rhs) > tiny(1.0_dp)) then
            solved = .false.
            if (gs_use_direct() .and. gs_direct_available()) &
                call direct_recovery(solved)
            if (solved) then
                status = KRYLOV_OK
            else
                recovered = 0
                call real_conjugate_gradient_operator(action, rhs, recovered, &
                    1.0e-11_dp*min(1.0_dp, norm2(rhs)), 6000, status, iterations, &
                    residual, precondition)
            end if
            if (present(recovery_converged)) recovery_converged = status == KRYLOV_OK
            if (present(recovery_iterations)) recovery_iterations = iterations
            ! The functional estimate is valid for every continuous recovered
            ! vector field, not only the minimizer. An iteration cap can make
            ! it less sharp, but must not discard an otherwise valid estimate.
            ! Algebraic breakdown and invalid arguments remain failures.
            if (status /= KRYLOV_OK .and. status /= KRYLOV_MAX_ITERATIONS) return
            if (.not. all(ieee_is_finite(recovered))) then
                status = 1
                return
            end if
        end if
        flux_sum = 0; balance_sum = 0
        failed_element = huge(1); failed_status = 1
        !$omp parallel default(shared)
        block
            real(dp) :: shape(ndof), gradients(2, ndof), hessians(2, 2, ndof)
            real(dp) :: d(ndof), local(2, ndof), point(2), det, grad(2), source
            real(dp) :: weight, flux(2), divergence, pp, ff, psi
            integer :: ie, q, stat, snapshot
            !$omp do schedule(static) reduction(+: flux_sum, balance_sum)
            do ie = 1, nel
                !$omp atomic read
                snapshot = failed_element
                if (ie > snapshot) cycle
                local = reshape(recovered(rec%indices(:, ie)), [2, ndof])
                d = state%psi(state%triangles(:, ie)) &
                    - state%psi(state%triangles(1, ie))
                do q = 1, size(w)
                    call gs_cell_data(state, ie, [x(q), y(q)], point, det, shape, &
                        gradients, hessians, stat)
                    if (stat /= 0) then
                        call record_failure(ie, stat)
                        exit
                    end if
                    weight = w(q)*det
                    grad = matmul(gradients, d)
                    flux = matmul(local, shape)
                    divergence = sum(local*gradients)
                    pp = pprime; ff = ffprime
                    if (present(law)) then
                        psi = state%psi(state%triangles(1, ie)) + dot_product(shape, d)
                        call gs_law_sources(law, psi, pp, ff)
                    end if
                    call gs_source(point(1), mu0, pp, ff, source)
                    if (present(forcing)) source = source + forcing(point(1), point(2))
                    flux_sum = flux_sum + weight*point(1) &
                        *sum((flux - grad/point(1))**2)
                    balance_sum = balance_sum + weight*(source + divergence)**2
                end do
            end do
            !$omp end do
        end block
        !$omp end parallel
        if (failed_element /= huge(1)) then
            status = failed_status
            return
        end if
        ! A capped recovery CG leaves KRYLOV_MAX_ITERATIONS in status; the
        ! estimate is valid (see recovery_converged), so the pass resets it.
        status = 0
        flux_defect = flux_sum; balance_defect = balance_sum
        flux_defect = sqrt(flux_defect)
        balance_defect = sqrt(cf2*balance_defect)
        majorant = flux_defect + balance_defect
        if (.not. ieee_is_finite(majorant)) then
            status = 1
            return
        end if
        if (present(law)) then
            ! Without branch stability this is the frozen-source defect, not
            ! a nonlinear error bound. A caller-supplied alpha is conditional
            ! on the documented strong-monotonicity / normalized-axis assumptions.
            status = 4
            if (practical) then
                call estimate_coercivity(alpha)
                if (present(stability_estimate)) stability_estimate = alpha
                if (alpha > 0) then
                    majorant = majorant/alpha
                    status = 5
                end if
            end if
            if (present(stability_lower_bound)) then
                if (.not. ieee_is_finite(stability_lower_bound)) return
                if (stability_lower_bound <= 0) return
                majorant = majorant/stability_lower_bound
                status = 0
            end if
        end if
    contains
        subroutine record_failure(element, code)
            ! The lowest failing element decides the status, as in a serial pass.
            integer, intent(in) :: element, code
            !$omp critical(gs_majorant_failure)
            if (element < failed_element) then
                failed_element = element
                failed_status = code
            end if
            !$omp end critical(gs_majorant_failure)
        end subroutine

        subroutine direct_recovery(converged)
            ! The recovery system is symmetric positive definite: one Cholesky
            ! factor (cached with the mesh), one back-substitution and iterative
            ! refinement. When the factor or the refinement fails, the caller
            ! runs the CG path.
            logical, intent(out) :: converged
            real(dp), allocatable :: residual_vector(:), correction(:)
            real(dp) :: rnorm, previous, target
            logical :: ok
            integer :: sweep
            converged = .false.
            rnorm = huge(1.0_dp)
            if (rec%factor%failed) return
            if (.not. rec%factor%factored) then
                call gs_direct_factor(rec%factor, n, rec%system%rowptr, &
                                      rec%system%colidx, rec%system%val, ok)
                if (.not. ok) return
            end if
            allocate (residual_vector(n), correction(n))
            call gs_direct_solve(rec%factor, rhs, recovered, ok)
            iterations = 1
            previous = huge(1.0_dp)
            target = 1.0e-11_dp*min(1.0_dp, norm2(rhs))
            do sweep = 0, 4
                if (.not. ok) exit
                call gs_csr_apply_full(rec%system, recovered, residual_vector)
                residual_vector = rhs - residual_vector
                rnorm = norm2(residual_vector)
                if (rnorm <= 1.0e-2_dp*target .or. sweep == 4) exit
                if (rnorm > 0.5_dp*previous) exit
                previous = rnorm
                call gs_direct_solve(rec%factor, residual_vector, correction, ok)
                recovered = recovered + correction
                iterations = iterations + 1
            end do
            converged = ok .and. rnorm <= target
            ! An uncached (hand-built) state keeps no factor.
            if (.not. cached) call gs_direct_free(rec%factor)
        end subroutine

        subroutine estimate_coercivity(alpha)
            ! Numerical positive-part comparison: J >= A - M_+.
            ! Power iteration is self-adjoint in the A pairing. Its result is
            ! a practical discrete stability estimate, not a continuous bound.
            real(dp), intent(out) :: alpha
            real(dp), allocatable :: v(:), av(:), mv(:), z(:), defect(:), ad(:)
            real(dp) :: theta, radius, energy
            real(dp) :: element(ndof)
            integer :: iteration, cell, node, i, j, stat, count
            count = size(state%psi)
            allocate (v(count), av(count), mv(count), defect(count), ad(count))
            v = state%psi
            if (norm2(v) <= tiny(1.0_dp)) v = 1
            where (state%boundary) v = 0
            alpha = 0
            do iteration = 1, 60
                call apply_gs(state, v, av)
                energy = dot_product(v, av)
                if (energy <= 0) return
                v = v/sqrt(energy); av = av/sqrt(energy)
                mv = 0
                do cell = 1, size(state%triangles, 2)
                    element = 0
                    do j = 1, ndof
                        do i = 1, ndof
                            element(i) = element(i) + source_mass(i, j, cell)* &
                                v(state%triangles(j, cell))
                        end do
                    end do
                    do j = 1, ndof
                        node = state%triangles(j, cell)
                        if (state%boundary(node)) cycle
                        mv(node) = mv(node) + element(j)
                    end do
                end do
                theta = dot_product(v, mv)
                if (norm2(mv) <= tiny(1.0_dp)) then
                    alpha = 1
                    return
                end if
                call gs_linear_response(state, mv, z, stat, &
                    relative_tolerance=1.0e-9_dp)
                if (stat /= 0) return
                defect = z - theta*v
                call apply_gs(state, defect, ad)
                radius = sqrt(max(0.0_dp, dot_product(defect, ad)))
                if (radius <= 1.0e-5_dp*max(1.0_dp, abs(theta))) then
                    alpha = 1 - theta - radius
                    return
                end if
                v = z
            end do
        end subroutine

        subroutine action(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            real(dp) :: value(2*ndof)
            integer :: cell, i, j
            output = 0
            do cell = 1, size(rec%indices, 2)
                value = 0
                do j = 1, 2*ndof
                    do i = 1, 2*ndof
                        value(i) = value(i) &
                            + rec%matrix(i, j, cell)*input(rec%indices(j, cell))
                    end do
                end do
                do i = 1, 2*ndof
                    output(rec%indices(i, cell)) = &
                        output(rec%indices(i, cell)) + value(i)
                end do
            end do
        end subroutine
        subroutine precondition(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            output = input/rec%diagonal
        end subroutine
    end subroutine
end module
