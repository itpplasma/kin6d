!> Energy-coupled NTV kernel on the spectral-element meshes of milca_dk_hp:
!! the operator of ntv_energy_coefficients (milca_dk_ntv) with the angular
!! part of hp_ntv_kernel, each speed node on its own mesh and geometry.
!!
!! Species a, speed nodes k = 2..n of the collop grid (h = 0 at x = 0 as in
!! the Legendre kernel), v_ak = v_a x_k. Node (a, k) carries h_ak on the mesh
!! mid(k, a) with the geometry geom(k, a). Divided by the streaming
!! coefficient (ratio of ntv_collision_rows), the equation of node (a, k)
!! for the drive A_j of species b reads
!!   A_ak h_ak - sum_p sum_k' E_a^(p)(k, k') P_{p<-k}^T M_p P_{p<-k'} h_ak'
!!     = delta_ab tz_a x_k^(2 j) load3_ak,
!! A_ak the node matrix of hp_node_matrix at nu = nu(k, a), v = v_ak and the
!! omega_e, pre of ops(k, a); E^(p) the energy-scattering rows with the speed
!! basis function of grid node p inside the speed integral (energy_w), M_p
!! the mass of the mesh and measure of grid node p (node 2 for p = 1, x = 0),
!! P = hp_transfer between the meshes; load3_ak the drive per unit speed.
!! On one mesh and measure the energy part is -sum_k' E_a(k, k') M h_ak'.
!! Responses
!!   K(a, b, j) = tz_a/(2 sqrt(pi) J_0) sum_k d4_k Re(load4_ak . h_ak),
!! q with d6_k, as ntv_energy_coefficients.
!!
!! Collision models: Lorentz (pitch-angle scattering only), test particle
!! (with energy scattering) and full linearised operator, which adds to
!! node (a, k)
!!   - sum_p P_{p<-k}^T sum_l D^l_p sum_b sum_k' F_ab^{l,(p)}(k, k') G^l_bk' h_bk',
!! l = 0..l_field: G^l_bk' the pitch moments (2 l + 1)/2 int P_l(xi) h dxi of
!! node (b, k') on theta lines common to all meshes, D^l_p their load
!! int measure P_l(xi) phi_i dxi on the mesh and measure of grid node p, and
!! F^(p) the field-particle rows with phi_p inside the speed integral
!! (field_w). Pitch-angle scattering is lumped at the speed nodes; energy
!! scattering and the field part conserve density up to the dropped x = 0
!! node.
!!
!! The coupled system is solved by GMRES, right-preconditioned by a coarse
!! Galerkin correction (per node the pitch-theta functions P_l(xi) cos/sin(m
!! theta), l < 3, m <= 4 and m within 2 of each perturbation harmonic, and
!! Chebyshev polynomials of the stream function on the counter-passing,
!! trapped and co-passing classes) followed by three symmetric block
!! Gauss-Seidel passes over the speed nodes with the factors of the diagonal
!! blocks A_ak - sum_p E_a^(p)(k, k) P_{p<-k}^T M_p P_{p<-k}. The field
!! part enters the coarse matrix only: in the sweeps it tripled their cost
!! and saved few iterations. GMRES minimises the residual weighted by the
!! speed mass v_a d3_k.
module milca_dk_hp_energy
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_krylov, only: complex_gmres_operator, KRYLOV_OK
    use fortnum_special_legendre, only: legendre_p
    use fortsparse, only: csc_t, csc_z_t, csc_from_triplet, csc_matmul, csc_transpose, &
        fortsparse_status_t, status_ok, sparse_factor_stats, sparse_factor_stats_t
    use milca_dk_hp, only: hp_geometry_t, hp_operator_t, hp_mesh_t, hp_terms_t, hp_ntv_terms, &
        hp_node_matrix, hp_mass, hp_transfer, hp_ntv_factor_t, hp_ntv_factor, hp_ntv_solve, &
        hp_ntv_free, hp_lines_t, hp_lines, hp_line_ops, HP_OK, HP_BAD_INPUT, HP_SINGULAR, &
        HP_NO_CONVERGENCE
    use milca_dk_ntv, only: ntv_collision_rows, NTV_OK
    use milca_dk_collop, only: collop_species_t, collop_grid_t, COLLOP_LORENTZ, &
        COLLOP_TEST_PARTICLE, COLLOP_FULL
    implicit none
    private
    public :: hp_energy_system_t, hp_energy_node_t, hp_energy_setup, hp_energy_collide, &
        hp_energy_coefficients, hp_energy_free
    !> The node factors exceed the memory budget.
    integer, parameter, public :: HP_NO_MEMORY = 4

    real(dp), parameter :: pi = 4*atan(1.0_dp)
    real(dp), parameter :: gmres_tolerance = 1e-8_dp
    integer, parameter :: gmres_restart = 40, default_iterations = 400
    ! Coarse space per node: P_l(xi) cos/sin(m theta), l < n_pitch, m <=
    ! n_mode and within n_side of each perturbation harmonic; Chebyshev
    ! degree < n_cheb of the stream function per class.
    integer, parameter :: n_pitch = 3, n_mode = 4, n_side = 2, n_cheb = 6
    ! Symmetric Gauss-Seidel passes per application: energy scattering
    ! couples the speed nodes stiffly when collisional.
    integer, parameter :: n_sgs = 3
    ! Highest Legendre degree of the field part unless l_field is given.
    integer, parameter :: default_l_field = 3

    !> One speed node: its terms on its mesh, the node matrix a, the place of
    !> its unknowns (off + 1 .. off + n) in the coupled vector, its factor and
    !> coarse functions z (orthonormal columns) from coarse column c_off + 1;
    !> with the field part its line moments g(l) and loads d(l).
    type :: hp_energy_node_t
        type(hp_terms_t) :: terms
        type(csc_t), allocatable :: g(:), d(:)
        type(csc_z_t) :: a
        integer :: off = 0, n = 0, mesh = 0, fac = 0, c_off = 0
        real(dp) :: nu = 0, v = 0
        real(dp), allocatable :: z(:, :)
    end type

    !> P_{p<-k}: unknowns of the mesh of grid node p from those of node k.
    type :: transfer_t
        logical :: identity = .true.
        type(csc_t) :: t
    end type

    type :: hp_energy_system_t
        integer :: ns = 0, nx = 0, n = 0, model = 0, n_dof = 0, n_coarse = 0, n_factors = 0, &
            l_field = 0, n_t = 0
        ! factor_peak: the most factor bytes held at once during setup.
        integer(int64) :: factor_bytes = 0, factor_peak = 0
        type(hp_energy_node_t), allocatable :: node(:, :)
        type(transfer_t), allocatable :: tr(:, :, :)
        type(hp_ntv_factor_t), allocatable :: fac(:)
        ! x of the grid nodes 1..n; per species vt, tz; per node d3, d4, d6
        ! and nu(k, a); energy_w(p, k, k', a) and field_w(p, k, k', a, b, l)
        ! as from ntv_collision_rows.
        real(dp), allocatable :: x(:), vt(:), tz(:), d3(:), d4(:), d6(:), nu(:, :), &
            energy(:, :, :), energy_w(:, :, :, :), field_w(:, :, :, :, :, :)
        complex(dp), allocatable :: coarse(:, :)
        integer, allocatable :: piv(:)
    end type

contains

    !> The coupled system of one surface: node terms, matrices, transfers,
    !> node factors and coarse matrix. geom, ops and mid are indexed (k, a),
    !> k = 1..n - 1 the speed nodes 2..n of grid; nu and v of ops are set
    !> here. Every node holds its own factor; if these exceed memory_budget
    !> bytes, HP_NO_MEMORY with sys%n_dof and sys%n_factors (the node count)
    !> set and no factor kept. The decision is taken from one factor per
    !> mesh, before the other factors exist. Limitation: there is no
    !> fallback for a smaller budget: factors shared by neighbouring speed
    !> nodes as block preconditioner left GMRES unconverged. l_field (default
    !> 3, at least 1 so that density and momentum return) is the highest
    !> Legendre degree of the field part of the full operator.
    subroutine hp_energy_setup(geom, ops, meshes, mid, grid, species, lnlambda, model, &
            memory_budget, sys, status, l_field)
        class(hp_geometry_t), intent(in) :: geom(:, :)
        type(hp_operator_t), intent(in) :: ops(:, :)
        type(hp_mesh_t), intent(in) :: meshes(:)
        integer, intent(in) :: mid(:, :), model
        type(collop_grid_t), intent(in) :: grid
        type(collop_species_t), intent(in) :: species(:)
        real(dp), intent(in) :: lnlambda(:, :)
        integer(int64), intent(in) :: memory_budget
        type(hp_energy_system_t), intent(out) :: sys
        integer, intent(out) :: status
        integer, intent(in), optional :: l_field
        real(dp), allocatable :: rows(:, :, :, :, :), ratio(:, :), lorentz(:), &
            full(:, :, :, :, :), r1(:, :), r2(:), r3(:), r4(:), r5(:), r6(:)
        real(dp) :: scale
        type(hp_lines_t) :: lines
        integer :: ns, nx, a, k, l, p, j, st, off

        status = HP_BAD_INPUT
        ns = size(species)
        nx = grid%n - 1
        if (all(model /= [COLLOP_LORENTZ, COLLOP_TEST_PARTICLE, COLLOP_FULL])) return
        sys%l_field = default_l_field
        if (present(l_field)) sys%l_field = l_field
        if (sys%l_field < 1) return
        if (ns < 1 .or. nx < 1) return
        if (any(shape(geom) /= [nx, ns]) .or. any(shape(ops) /= [nx, ns]) .or. &
            any(shape(mid) /= [nx, ns])) return
        if (any(mid < 1) .or. any(mid > size(meshes))) return
        if (any(ops%source_dpar) .or. any(ops%source_dlf10)) return
        sys%ns = ns
        sys%nx = nx
        sys%n = grid%n
        sys%model = model
        ! The rows of the test-particle part: with the field part their
        ! l-dependence would not be the Lorentz part alone.
        call ntv_collision_rows(species, lnlambda, grid, min(model, COLLOP_TEST_PARTICLE), 3, &
            rows, ratio, sys%d3, sys%d4, sys%d6, sys%vt, sys%tz, st, energy=sys%energy, &
            nu=sys%nu, energy_w=sys%energy_w)
        if (st /= NTV_OK) return
        if (model == COLLOP_FULL) then
            call ntv_collision_rows(species, lnlambda, grid, model, sys%l_field + 1, full, &
                r1, r2, r3, r4, r5, r6, st, field_w=sys%field_w)
            if (st /= NTV_OK) return
        end if
        ! Pitch-angle scattering enters the node matrices pointwise as nu C;
        ! that needs the l-dependence of the rows to be exactly the lumped
        ! Lorentz part -l (l + 1)/2 nu.
        do a = 1, ns
            do k = 1, nx
                scale = maxval(abs(rows(k, :, a, a, :)))
                do l = 1, 2
                    lorentz = (rows(k, :, a, a, l) - rows(k, :, a, a, 0))/ratio(k, a)
                    lorentz(k) = lorentz(k) + l*(l + 1)/2.0_dp*sys%nu(k, a)
                    if (maxval(abs(lorentz))*ratio(k, a) > 1e-10_dp*scale) return
                end do
            end do
        end do
        sys%x = grid%x

        allocate (sys%node(nx, ns))
        off = 0
        do a = 1, ns
            do k = 1, nx
                associate (nd => sys%node(k, a))
                    nd%mesh = mid(k, a)
                    call hp_ntv_terms(geom(k, a), meshes(nd%mesh), ops(k, a), nd%terms, status)
                    if (status /= HP_OK) return
                    status = HP_BAD_INPUT
                    nd%nu = sys%nu(k, a)
                    nd%v = sys%vt(a)*grid%x(k + 1)
                    call hp_node_matrix(nd%terms, nd%nu, nd%v, ops(k, a)%omega_e, &
                        ops(k, a)%pre, ops(k, a)%n, ops(k, a)%q, nd%a)
                    if (nd%a%nrow /= nd%terms%n) return
                    nd%n = nd%terms%n
                    nd%off = off
                    off = off + nd%n
                end associate
            end do
        end do
        sys%n_dof = off

        ! Transfers P_{p<-k} wherever grid node p measures node k.
        allocate (sys%tr(grid%n, nx, ns))
        do a = 1, ns
            do k = 1, nx
                do p = 1, grid%n
                    if (all(sys%energy_w(p, k, :, a) == 0) .and. &
                        all(sys%energy_w(p, :, k, a) == 0) .and. .not. field_from(sys, p, k, a)) &
                        cycle
                    j = max(p - 1, 1)
                    if (mid(j, a) == mid(k, a)) cycle
                    sys%tr(p, k, a)%identity = .false.
                    call hp_transfer(geom(k, a), meshes(mid(j, a)), meshes(mid(k, a)), &
                        sys%node(j, a)%terms%loc, sys%node(k, a)%terms%loc, &
                        sys%tr(p, k, a)%t, status)
                    if (status /= HP_OK) return
                end do
            end do
        end do

        ! Pitch moments of every node and their loads on the theta lines
        ! common to all meshes.
        if (model == COLLOP_FULL) then
            call hp_lines(meshes, lines)
            sys%n_t = size(lines%theta)
            do a = 1, ns
                do k = 1, nx
                    call hp_line_ops(geom(k, a), meshes(mid(k, a)), ops(k, a), lines, &
                        sys%l_field, sys%node(k, a)%g, sys%node(k, a)%d, status)
                    if (status /= HP_OK) return
                end do
            end do
        end if

        call node_factors(sys, memory_budget, status)
        if (status /= HP_OK) return
        do a = 1, ns
            do k = 1, nx
                call coarse_space(geom(k, a), meshes(mid(k, a)), sys%node(k, a)%terms, &
                    ops(k, a)%m, sys%node(k, a)%z)
            end do
        end do
        call coarse_matrix(sys, status)
    end subroutine hp_energy_setup

    !> Factors of the diagonal blocks, one per node. Whether all of them fit
    !> the budget is estimated from one factor per mesh before the others are
    !> built, so that a refusal costs only those; if the estimate is low the
    !> build stops once the bytes held pass the budget. Either way all factors
    !> are freed and HP_NO_MEMORY returned.
    subroutine node_factors(sys, budget, status)
        type(hp_energy_system_t), intent(inout) :: sys
        integer(int64), intent(in) :: budget
        integer, intent(out) :: status
        integer, allocatable :: rep(:)
        integer(int64), allocatable :: bytes(:), rep_bytes(:)
        integer(int64) :: need
        logical, allocatable :: built(:)
        integer :: a, k, m, i

        ! Factor slot of node (k, a): (a - 1) nx + k.
        allocate (sys%fac(sys%nx*sys%ns), bytes(sys%nx*sys%ns), built(sys%nx*sys%ns))
        bytes = 0
        built = .false.
        sys%factor_bytes = 0
        sys%factor_peak = 0
        sys%n_factors = sys%nx*sys%ns
        do a = 1, sys%ns
            do k = 1, sys%nx
                sys%node(k, a)%fac = (a - 1)*sys%nx + k
            end do
        end do
        allocate (rep(maxval(sys%node%mesh)), rep_bytes(maxval(sys%node%mesh)))
        rep = 0
        rep_bytes = 0
        reps: do a = 1, sys%ns
            do k = 1, sys%nx
                m = sys%node(k, a)%mesh
                if (rep(m) /= 0) cycle
                i = (a - 1)*sys%nx + k
                rep(m) = i
                call factor_node(sys, k, a, sys%fac(i), bytes(i), status)
                if (status /= HP_OK) exit reps
                built(i) = .true.
                rep_bytes(m) = bytes(i)
                sys%factor_bytes = sys%factor_bytes + bytes(i)
                sys%factor_peak = max(sys%factor_peak, sys%factor_bytes)
                if (sys%factor_bytes > budget) exit reps
            end do
        end do reps
        if (status == HP_OK) then
            need = 0
            do a = 1, sys%ns
                do k = 1, sys%nx
                    need = need + rep_bytes(sys%node(k, a)%mesh)
                end do
            end do
            if (sys%factor_bytes > budget .or. need > budget) then
                status = HP_NO_MEMORY
            else
                call build_factors(sys, budget, built, bytes, status)
                if (status == HP_OK) return
            end if
        end if
        do i = 1, size(sys%fac)
            if (.not. built(i)) cycle
            call hp_ntv_free(sys%fac(i))
            sys%factor_bytes = sys%factor_bytes - bytes(i)
        end do
        deallocate (sys%fac)
    end subroutine node_factors

    !> The factors not yet built, in parallel; once the bytes held pass the
    !> budget no further factor is started, HP_NO_MEMORY.
    subroutine build_factors(sys, budget, built, bytes, status)
        type(hp_energy_system_t), intent(inout) :: sys
        integer(int64), intent(in) :: budget
        logical, intent(inout) :: built(:)
        integer(int64), intent(inout) :: bytes(:)
        integer, intent(out) :: status
        integer, allocatable :: fails(:)
        integer(int64) :: b
        integer :: a, k, i, st
        logical :: over, stop_now

        allocate (fails(sys%nx*sys%ns))
        fails = HP_OK
        over = sys%factor_bytes > budget
        !$omp parallel do default(shared) schedule(dynamic) private(i, a, k, st, b, stop_now)
        do i = 1, sys%nx*sys%ns
            k = mod(i - 1, sys%nx) + 1
            a = (i - 1)/sys%nx + 1
            if (built(i)) cycle
            !$omp atomic read
            stop_now = over
            if (stop_now) cycle
            call factor_node(sys, k, a, sys%fac(i), b, st)
            fails(i) = st
            if (st /= HP_OK) cycle
            !$omp critical (energy_factor_bytes)
            built(i) = .true.
            bytes(i) = b
            sys%factor_bytes = sys%factor_bytes + b
            sys%factor_peak = max(sys%factor_peak, sys%factor_bytes)
            if (sys%factor_bytes > budget) then
                !$omp atomic write
                over = .true.
            end if
            !$omp end critical (energy_factor_bytes)
        end do
        !$omp end parallel do
        status = maxval(fails)
        if (status == HP_OK .and. over) status = HP_NO_MEMORY
    end subroutine build_factors

    !> Factor of the diagonal block of node (k, a) and its bytes.
    subroutine factor_node(sys, k, a, fac, bytes, status)
        type(hp_energy_system_t), intent(in) :: sys
        integer, intent(in) :: k, a
        type(hp_ntv_factor_t), intent(inout) :: fac
        integer(int64), intent(out) :: bytes
        integer, intent(out) :: status
        type(csc_z_t) :: shifted
        type(sparse_factor_stats_t) :: stats
        type(fortsparse_status_t) :: fst

        bytes = 0
        call diagonal_block(sys, k, a, shifted, status)
        if (status == HP_OK) call hp_ntv_factor(shifted, fac, status)
        if (status /= HP_OK) return
        call sparse_factor_stats(fac%solver, stats, fst)
        ! Without a native byte count, 16 bytes per stored entry.
        if (status_ok(fst) .and. stats%factor_storage_bytes >= 0) then
            bytes = stats%factor_storage_bytes
        else if (status_ok(fst) .and. stats%factor_nnz >= 0) then
            bytes = 16*stats%factor_nnz
        else
            bytes = 16_int64*sys%node(k, a)%a%nnz
        end if
    end subroutine factor_node

    !> The diagonal block of node (k, a) of the coupled operator, A_ak -
    !> sum_p E_a^(p)(k, k) P_{p<-k}^T M_p P_{p<-k}: on other meshes than its
    !> own the mass seen through the transfer is rank deficient, and a factor
    !> of A_ak - E_a(k, k) M_ak would miss that where energy scattering is stiff.
    subroutine diagonal_block(sys, k, a, d, status)
        type(hp_energy_system_t), intent(in) :: sys
        integer, intent(in) :: k, a
        type(csc_z_t), intent(out) :: d
        integer, intent(out) :: status
        type(fortsparse_status_t) :: st
        type(csc_t) :: m, mp, pt, e
        integer, allocatable :: rows(:), cols(:)
        complex(dp), allocatable :: vals(:)
        real(dp) :: w
        integer :: p, j, i

        status = HP_BAD_INPUT
        associate (nd => sys%node(k, a))
            rows = [(nd%a%row_idx(i), i=1, nd%a%nnz)]
            cols = [((j, i=nd%a%col_ptr(j), nd%a%col_ptr(j + 1) - 1), j=1, nd%a%ncol)]
            vals = nd%a%val(:nd%a%nnz)
        end associate
        if (sys%model /= COLLOP_LORENTZ) then
            do p = 1, sys%n
                w = sys%energy_w(p, k, k, a)
                if (w == 0) cycle
                call hp_mass(sys%node(max(p - 1, 1), a)%terms, m)
                if (sys%tr(p, k, a)%identity) then
                    e = m
                else
                    call csc_matmul(m, sys%tr(p, k, a)%t, mp, st)
                    if (.not. status_ok(st)) return
                    call csc_transpose(sys%tr(p, k, a)%t, pt, st)
                    if (.not. status_ok(st)) return
                    call csc_matmul(pt, mp, e, st)
                    if (.not. status_ok(st)) return
                end if
                rows = [rows, e%row_idx(:e%nnz)]
                cols = [cols, ((j, i=e%col_ptr(j), e%col_ptr(j + 1) - 1), j=1, e%ncol)]
                vals = [vals, cmplx(-w*e%val(:e%nnz), 0.0_dp, dp)]
            end do
        end if
        call csc_from_triplet(sys%node(k, a)%n, sys%node(k, a)%n, rows, cols, vals, d, st)
        if (status_ok(st)) status = HP_OK
    end subroutine diagonal_block

    !> Coarse functions of one node at its unknowns, orthonormalised; pert
    !> the poloidal harmonics of the perturbation.
    subroutine coarse_space(g, mesh, terms, pert, z)
        class(hp_geometry_t), intent(in) :: g
        type(hp_mesh_t), intent(in) :: mesh
        type(hp_terms_t), intent(in) :: terms
        integer, intent(in) :: pert(:)
        real(dp), allocatable, intent(out) :: z(:, :)
        real(dp), allocatable :: s(:), t(:), w(:, :)
        integer, allocatable :: cls(:), modes(:)
        integer :: p, n_col, j_lo, j_hi, e, a1, a2, i, j, l, m, c, d, nc
        real(dp) :: lo, hi

        p = mesh%p
        n_col = size(mesh%theta) - 1
        j_lo = (size(mesh%lam) - 1)*p
        j_hi = j_lo + (size(mesh%zeta) - 1)*p
        allocate (cls(terms%n), s(terms%n), t(terms%n))
        ! Row index of each unknown: below the trapped rows counter-passing,
        ! above co-passing.
        do e = 1, size(terms%loc, 2)
            do a2 = 0, p
                do a1 = 0, p
                    i = terms%loc(1 + a1 + (p + 1)*a2, e)
                    j = ((e - 1)/n_col)*p + a2
                    cls(i) = merge(1, merge(3, 2, j > j_hi), j < j_lo)
                end do
            end do
        end do
        do i = 1, terms%n
            s(i) = g%stream(terms%theta(i), terms%xi(i))
        end do
        ! Collisional solutions follow the drive in theta: its harmonics and
        ! their sidebands from the variation of b join the lowest modes.
        modes = [(m, m=0, n_mode)]
        do i = 1, size(pert)
            do m = abs(pert(i)) - n_side, abs(pert(i)) + n_side
                if (m >= 0 .and. all(modes /= m)) modes = [modes, m]
            end do
        end do
        allocate (w(terms%n, n_pitch*(2*size(modes) - 1) + 3*n_cheb))
        nc = 0
        do l = 0, n_pitch - 1
            do j = 1, size(modes)
                m = modes(j)
                nc = nc + 1
                w(:, nc) = [(legendre_p(l, 0, terms%xi(i)), i=1, terms%n)]*cos(m*terms%theta)
                if (m == 0) cycle
                nc = nc + 1
                w(:, nc) = [(legendre_p(l, 0, terms%xi(i)), i=1, terms%n)]*sin(m*terms%theta)
            end do
        end do
        do c = 1, 3
            if (.not. any(cls == c)) cycle
            lo = minval(s, cls == c)
            hi = maxval(s, cls == c)
            t = 0
            if (hi > lo) t = max(-1.0_dp, min(1.0_dp, 2*(s - lo)/(hi - lo) - 1))
            do d = 0, n_cheb - 1
                nc = nc + 1
                w(:, nc) = merge(cos(d*acos(t)), 0.0_dp, cls == c)
            end do
        end do
        call orthonormal(w(:, :nc), z)
    end subroutine coarse_space

    !> Orthonormal basis of the span of the columns of w (modified
    !> Gram-Schmidt, twice), dependent columns dropped.
    subroutine orthonormal(w, z)
        real(dp), intent(in) :: w(:, :)
        real(dp), allocatable, intent(out) :: z(:, :)
        real(dp) :: v(size(w, 1)), tmp(size(w, 1), size(w, 2)), n0
        integer :: c, i, m, pass
        m = 0
        do c = 1, size(w, 2)
            v = w(:, c)
            n0 = norm2(v)
            if (n0 == 0) cycle
            do pass = 1, 2
                do i = 1, m
                    v = v - dot_product(tmp(:, i), v)*tmp(:, i)
                end do
            end do
            if (norm2(v) <= 1e-8_dp*n0) cycle
            m = m + 1
            tmp(:, m) = v/norm2(v)
        end do
        z = tmp(:, :m)
    end subroutine orthonormal

    !> Coarse matrix Z^T N Z, factorised.
    subroutine coarse_matrix(sys, status)
        type(hp_energy_system_t), intent(inout) :: sys
        integer, intent(out) :: status
        complex(dp), allocatable :: h(:), y(:)
        integer :: a, k, j, c, i, nc, info
        external :: zgetrf

        nc = 0
        do a = 1, sys%ns
            do k = 1, sys%nx
                sys%node(k, a)%c_off = nc
                nc = nc + size(sys%node(k, a)%z, 2)
            end do
        end do
        sys%n_coarse = nc
        allocate (sys%coarse(nc, nc), sys%piv(nc))
        sys%coarse = 0
        !$omp parallel do default(shared) schedule(dynamic) private(i, a, j, k, c, h, y)
        do i = 1, sys%nx*sys%ns
            j = mod(i - 1, sys%nx) + 1
            a = (i - 1)/sys%nx + 1
            allocate (h(sys%n_dof))
            h = 0
            associate (nj => sys%node(j, a))
                do c = 1, size(nj%z, 2)
                    h(nj%off + 1:nj%off + nj%n) = nj%z(:, c)
                    do k = 1, sys%nx
                        if (.not. coupled(sys, k, j, a)) cycle
                        associate (nk => sys%node(k, a))
                            allocate (y(nk%n))
                            call row_apply(sys, a, k, h, y, j)
                            sys%coarse(nk%c_off + 1:nk%c_off + size(nk%z, 2), nj%c_off + c) = &
                                matmul(transpose(nk%z), y)
                            deallocate (y)
                        end associate
                    end do
                end do
            end associate
            deallocate (h)
        end do
        !$omp end parallel do
        call coarse_field(sys)
        status = HP_SINGULAR
        call zgetrf(nc, nc, sys%coarse, nc, sys%piv, info)
        if (info /= 0) return
        status = HP_OK
    end subroutine coarse_matrix

    !> Field part of the coarse matrix, -sum_p sum_l F_ab^{l,(p)}(k, j)
    !> (Z_ak^T P_{p<-k}^T D^l_p)(G^l_bj Z_bj) for row node (k, a) and column
    !> node (j, b), from the products with the line moments: the coarse
    !> functions never leave their mesh.
    subroutine coarse_field(sys)
        type(hp_energy_system_t), intent(inout) :: sys
        type :: lines_t
            real(dp), allocatable :: v(:, :, :)
        end type
        type(lines_t), allocatable :: gz(:, :)
        real(dp), allocatable :: zd(:, :, :, :), s(:, :), pz(:, :), blk(:, :)
        integer, allocatable :: ps(:)
        real(dp) :: w
        integer :: i, a, k, b, j, p, l, q, pk, nck, ncj

        if (sys%model /= COLLOP_FULL) return
        allocate (gz(sys%nx, sys%ns))
        !$omp parallel do default(shared) schedule(dynamic) private(i, a, k, l)
        do i = 1, sys%nx*sys%ns
            k = mod(i - 1, sys%nx) + 1
            a = (i - 1)/sys%nx + 1
            associate (nd => sys%node(k, a))
                allocate (gz(k, a)%v(sys%n_t, size(nd%z, 2), 0:sys%l_field))
                do l = 0, sys%l_field
                    gz(k, a)%v(:, :, l) = sparse_times(nd%g(l), nd%z, .false.)
                end do
            end associate
        end do
        !$omp end parallel do
        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp private(i, a, k, b, j, p, l, q, pk, nck, ncj, w, ps, zd, s, pz, blk)
        do i = 1, sys%nx*sys%ns
            k = mod(i - 1, sys%nx) + 1
            a = (i - 1)/sys%nx + 1
            ps = pack([(p, p=1, sys%n)], [(field_from(sys, p, k, a), p=1, sys%n)])
            associate (nk => sys%node(k, a))
                nck = size(nk%z, 2)
                allocate (zd(sys%n_t, nck, 0:sys%l_field, size(ps)), s(sys%n_t, nck))
                do q = 1, size(ps)
                    pk = max(ps(q) - 1, 1)
                    if (sys%tr(ps(q), k, a)%identity) then
                        pz = nk%z
                    else
                        pz = sparse_times(sys%tr(ps(q), k, a)%t, nk%z, .false.)
                    end if
                    do l = 0, sys%l_field
                        zd(:, :, l, q) = sparse_times(sys%node(pk, a)%d(l), pz, .true.)
                    end do
                end do
                do b = 1, sys%ns
                    do j = 1, sys%nx
                        associate (nj => sys%node(j, b))
                            ncj = size(nj%z, 2)
                            allocate (blk(nck, ncj))
                            blk = 0
                            do l = 0, sys%l_field
                                s = 0
                                do q = 1, size(ps)
                                    w = sys%field_w(ps(q), k, j, a, b, l)
                                    if (w /= 0) s = s + w*zd(:, :, l, q)
                                end do
                                blk = blk + matmul(transpose(s), gz(j, b)%v(:, :, l))
                            end do
                            sys%coarse(nk%c_off + 1:nk%c_off + nck, nj%c_off + 1:nj%c_off &
                                + ncj) = sys%coarse(nk%c_off + 1:nk%c_off + nck, &
                                nj%c_off + 1:nj%c_off + ncj) - blk
                            deallocate (blk)
                        end associate
                    end do
                end do
                deallocate (zd, s)
            end associate
        end do
        !$omp end parallel do
    end subroutine coarse_field

    !> A x, or A^T x with transpose, for a real sparse A and dense x.
    function sparse_times(m, x, transpose) result(y)
        type(csc_t), intent(in) :: m
        real(dp), intent(in) :: x(:, :)
        logical, intent(in) :: transpose
        real(dp), allocatable :: y(:, :)
        integer :: j, i
        if (transpose) then
            allocate (y(m%ncol, size(x, 2)))
            y = 0
            do j = 1, m%ncol
                do i = m%col_ptr(j), m%col_ptr(j + 1) - 1
                    y(j, :) = y(j, :) + m%val(i)*x(m%row_idx(i), :)
                end do
            end do
        else
            allocate (y(m%nrow, size(x, 2)))
            y = 0
            do j = 1, m%ncol
                do i = m%col_ptr(j), m%col_ptr(j + 1) - 1
                    y(m%row_idx(i), :) = y(m%row_idx(i), :) + m%val(i)*x(j, :)
                end do
            end do
        end if
    end function sparse_times

    !> Line moments gm(:, l, k, a) = G^l_ak h_ak of the coupled vector h.
    subroutine field_moments(sys, h, gm)
        type(hp_energy_system_t), intent(in) :: sys
        complex(dp), intent(in) :: h(:)
        complex(dp), allocatable, intent(out) :: gm(:, :, :, :)
        integer :: a, k
        allocate (gm(sys%n_t, 0:sys%l_field, sys%nx, sys%ns))
        do a = 1, sys%ns
            do k = 1, sys%nx
                associate (nd => sys%node(k, a))
                    call node_moments(sys, k, a, h(nd%off + 1:nd%off + nd%n), gm(:, :, k, a))
                end associate
            end do
        end do
    end subroutine field_moments

    !> g(:, l) = G^l x for x on the mesh of node (k, a).
    subroutine node_moments(sys, k, a, x, g)
        type(hp_energy_system_t), intent(in) :: sys
        integer, intent(in) :: k, a
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: g(:, 0:)
        integer :: l
        do l = 0, sys%l_field
            g(:, l) = csc_r_apply(sys%node(k, a)%g(l), x)
        end do
    end subroutine node_moments

    !> y += sign sum_p P_{p<-k}^T sum_l D^l_p sum_b sum_k' F_ab^{l,(p)}(k, k')
    !> gm(:, l, k', b), gm the line moments of the coupled vector.
    subroutine field_row(sys, a, k, gm, y, sign)
        type(hp_energy_system_t), intent(in) :: sys
        integer, intent(in) :: a, k
        complex(dp), intent(in) :: gm(:, 0:, :, :)
        complex(dp), intent(inout) :: y(:)
        real(dp), intent(in) :: sign
        complex(dp), allocatable :: t(:), s(:)
        real(dp) :: w
        integer :: p, pk, l, b, j

        allocate (s(sys%n_t))
        do p = 1, sys%n
            if (.not. field_from(sys, p, k, a)) cycle
            pk = max(p - 1, 1)
            allocate (t(sys%node(pk, a)%n))
            t = 0
            do l = 0, sys%l_field
                s = 0
                do b = 1, sys%ns
                    do j = 1, sys%nx
                        w = sys%field_w(p, k, j, a, b, l)
                        if (w /= 0) s = s + w*gm(:, l, j, b)
                    end do
                end do
                t = t + csc_r_apply(sys%node(pk, a)%d(l), s)
            end do
            y = y + sign*moved(sys%tr(p, k, a), t, .true.)
            deallocate (t)
        end do
    end subroutine field_row

    pure function csc_r_apply(a, x) result(y)
        type(csc_t), intent(in) :: a
        complex(dp), intent(in) :: x(:)
        complex(dp) :: y(a%nrow)
        integer :: j, i
        y = 0
        do j = 1, a%ncol
            do i = a%col_ptr(j), a%col_ptr(j + 1) - 1
                y(a%row_idx(i)) = y(a%row_idx(i)) + a%val(i)*x(j)
            end do
        end do
    end function csc_r_apply

    !> Whether grid node p measures the field part of row node (k, a).
    logical function field_from(sys, p, k, a)
        type(hp_energy_system_t), intent(in) :: sys
        integer, intent(in) :: p, k, a
        field_from = .false.
        if (sys%model == COLLOP_FULL) field_from = any(sys%field_w(p, k, :, a, :, :) /= 0)
    end function field_from

    logical function coupled(sys, k, j, a)
        type(hp_energy_system_t), intent(in) :: sys
        integer, intent(in) :: k, j, a
        coupled = k == j
        if (sys%model /= COLLOP_LORENTZ) coupled = coupled .or. any(sys%energy_w(:, k, j, a) /= 0)
    end function coupled

    !> y = (N h) of node (k, a), h the coupled vector; with only, the part
    !> of the unknowns of node (only, a) alone; the field part only with the
    !> line moments gm of h.
    subroutine row_apply(sys, a, k, h, y, only, gm)
        type(hp_energy_system_t), intent(in) :: sys
        integer, intent(in) :: a, k
        complex(dp), intent(in) :: h(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(in), optional :: only
        complex(dp), intent(in), optional :: gm(:, 0:, :, :)
        y = 0
        associate (nd => sys%node(k, a))
            if (present(only)) then
                if (only == k) y = csc_z_apply(nd%a, h(nd%off + 1:nd%off + nd%n))
            else
                y = csc_z_apply(nd%a, h(nd%off + 1:nd%off + nd%n))
            end if
        end associate
        call energy_row(sys, a, k, h, y, -1.0_dp, only)
        if (present(gm)) call field_row(sys, a, k, gm, y, -1.0_dp)
    end subroutine row_apply

    !> y += sign sum_p sum_k' E_a^(p)(k, k') P_{p<-k}^T M_p P_{p<-k'} h_ak'.
    subroutine energy_row(sys, a, k, h, y, sign, only)
        type(hp_energy_system_t), intent(in) :: sys
        integer, intent(in) :: a, k
        complex(dp), intent(in) :: h(:)
        complex(dp), intent(inout) :: y(:)
        real(dp), intent(in) :: sign
        integer, intent(in), optional :: only
        complex(dp), allocatable :: t(:)
        real(dp) :: w
        integer :: p, j, pk
        logical :: any_j

        if (sys%model == COLLOP_LORENTZ) return
        do p = 1, sys%n
            if (all(sys%energy_w(p, k, :, a) == 0)) cycle
            pk = max(p - 1, 1)
            allocate (t(sys%node(pk, a)%n))
            t = 0
            any_j = .false.
            do j = 1, sys%nx
                if (present(only)) then
                    if (j /= only) cycle
                end if
                w = sys%energy_w(p, k, j, a)
                if (w == 0) cycle
                any_j = .true.
                associate (nj => sys%node(j, a))
                    t = t + w*moved(sys%tr(p, j, a), h(nj%off + 1:nj%off + nj%n), .false.)
                end associate
            end do
            if (any_j) then
                t = element_apply(sys%node(pk, a)%terms%loc, sys%node(pk, a)%terms%m, t)
                y = y + sign*moved(sys%tr(p, k, a), t, .true.)
            end if
            deallocate (t)
        end do
    end subroutine energy_row

    !> P x, or P^T x with transpose, of a transfer.
    function moved(tr, x, transpose) result(y)
        type(transfer_t), intent(in) :: tr
        complex(dp), intent(in) :: x(:)
        logical, intent(in) :: transpose
        complex(dp), allocatable :: y(:)
        integer :: j, i
        if (tr%identity) then
            y = x
            return
        end if
        if (transpose) then
            allocate (y(tr%t%ncol))
            y = 0
            do j = 1, tr%t%ncol
                do i = tr%t%col_ptr(j), tr%t%col_ptr(j + 1) - 1
                    y(j) = y(j) + tr%t%val(i)*x(tr%t%row_idx(i))
                end do
            end do
        else
            allocate (y(tr%t%nrow))
            y = 0
            do j = 1, tr%t%ncol
                do i = tr%t%col_ptr(j), tr%t%col_ptr(j + 1) - 1
                    y(tr%t%row_idx(i)) = y(tr%t%row_idx(i)) + tr%t%val(i)*x(j)
                end do
            end do
        end if
    end function moved

    pure function csc_z_apply(a, x) result(y)
        type(csc_z_t), intent(in) :: a
        complex(dp), intent(in) :: x(:)
        complex(dp) :: y(a%nrow)
        integer :: j, i
        y = 0
        do j = 1, a%ncol
            do i = a%col_ptr(j), a%col_ptr(j + 1) - 1
                y(a%row_idx(i)) = y(a%row_idx(i)) + a%val(i)*x(j)
            end do
        end do
    end function csc_z_apply

    !> Sum of the element matrices mats (column major) applied to x.
    pure function element_apply(loc, mats, x) result(y)
        integer, intent(in) :: loc(:, :)
        real(dp), intent(in) :: mats(:, :)
        complex(dp), intent(in) :: x(:)
        complex(dp) :: y(size(x))
        complex(dp) :: t(size(loc, 1))
        integer :: e, nl, i
        nl = size(loc, 1)
        y = 0
        do e = 1, size(loc, 2)
            t = matmul(reshape(mats(:, e), [nl, nl]), x(loc(:, e)))
            ! Element by element: the nodes collapsed on the X point repeat
            ! an unknown in loc, which a vector subscript would not sum.
            do i = 1, nl
                y(loc(i, e)) = y(loc(i, e)) + t(i)
            end do
        end do
    end function element_apply

    !> The collision operator in the weak form of the node equations, y_ak =
    !> -nu_ak C_ak h_ak + sum_p sum_k' E_a^(p)(k, k') P^T M_p P h_ak' and the
    !> field part (the negative of its part of the coupled operator), h and
    !> y coupled vectors.
    subroutine hp_energy_collide(sys, h, y)
        type(hp_energy_system_t), intent(in) :: sys
        complex(dp), intent(in) :: h(:)
        complex(dp), intent(out) :: y(:)
        complex(dp), allocatable :: gm(:, :, :, :)
        integer :: i, a, k

        if (sys%model == COLLOP_FULL) call field_moments(sys, h, gm)
        !$omp parallel do default(shared) schedule(dynamic) private(i, a, k)
        do i = 1, sys%nx*sys%ns
            k = mod(i - 1, sys%nx) + 1
            a = (i - 1)/sys%nx + 1
            associate (nd => sys%node(k, a))
                y(nd%off + 1:nd%off + nd%n) = -nd%nu*element_apply(nd%terms%loc, &
                    nd%terms%c, h(nd%off + 1:nd%off + nd%n))
                call energy_row(sys, a, k, h, y(nd%off + 1:nd%off + nd%n), 1.0_dp)
                if (allocated(gm)) call field_row(sys, a, k, gm, y(nd%off + 1:nd%off &
                    + nd%n), 1.0_dp)
            end associate
        end do
        !$omp end parallel do
    end subroutine hp_energy_collide

    subroutine full_apply(sys, h, y)
        type(hp_energy_system_t), intent(in) :: sys
        complex(dp), intent(in) :: h(:)
        complex(dp), intent(out) :: y(:)
        complex(dp), allocatable :: gm(:, :, :, :)
        integer :: i, a, k

        if (sys%model == COLLOP_FULL) call field_moments(sys, h, gm)
        !$omp parallel do default(shared) schedule(dynamic) private(i, a, k)
        do i = 1, sys%nx*sys%ns
            k = mod(i - 1, sys%nx) + 1
            a = (i - 1)/sys%nx + 1
            associate (nd => sys%node(k, a))
                if (allocated(gm)) then
                    call row_apply(sys, a, k, h, y(nd%off + 1:nd%off + nd%n), gm=gm)
                else
                    call row_apply(sys, a, k, h, y(nd%off + 1:nd%off + nd%n))
                end if
            end associate
        end do
        !$omp end parallel do
    end subroutine full_apply

    !> u = Z Ac^-1 Z^T r, then n_sgs symmetric block Gauss-Seidel passes on r - N u.
    subroutine precondition(sys, r, u, status)
        type(hp_energy_system_t), intent(inout) :: sys
        complex(dp), intent(in) :: r(:)
        complex(dp), intent(out) :: u(:)
        integer, intent(out) :: status
        complex(dp), allocatable :: c(:), res(:), d(:), y(:), dx(:)
        integer :: a, k, sweep, it, info
        external :: zgetrs

        allocate (c(sys%n_coarse), res(sys%n_dof), d(sys%n_dof))
        do a = 1, sys%ns
            do k = 1, sys%nx
                associate (nd => sys%node(k, a))
                    c(nd%c_off + 1:nd%c_off + size(nd%z, 2)) = &
                        matmul(transpose(nd%z), r(nd%off + 1:nd%off + nd%n))
                end associate
            end do
        end do
        call zgetrs('N', sys%n_coarse, 1, sys%coarse, sys%n_coarse, sys%piv, c, &
            sys%n_coarse, info)
        do a = 1, sys%ns
            do k = 1, sys%nx
                associate (nd => sys%node(k, a))
                    u(nd%off + 1:nd%off + nd%n) = matmul(nd%z, c(nd%c_off + 1:nd%c_off &
                        + size(nd%z, 2)))
                end associate
            end do
        end do
        call full_apply(sys, u, res)
        res = r - res
        d = 0
        status = HP_OK
        do a = 1, sys%ns
            do sweep = 1, 2*n_sgs
                do it = 1, sys%nx
                    k = merge(it, sys%nx + 1 - it, mod(sweep, 2) == 1)
                    associate (nd => sys%node(k, a))
                        allocate (y(nd%n), dx(nd%n))
                        call row_apply(sys, a, k, d, y)
                        call hp_ntv_solve(sys%fac(nd%fac), res(nd%off + 1:nd%off + nd%n) - y, &
                            dx, .false., status, 0)
                        if (status /= HP_OK) return
                        d(nd%off + 1:nd%off + nd%n) = d(nd%off + 1:nd%off + nd%n) + dx
                        deallocate (y, dx)
                    end associate
                end do
            end do
        end do
        u = u + d
    end subroutine precondition

    !> K(a, b, j) and q(a, b, j) of the drives A_j of species b (module
    !> header), node_r(k, a, b, j) = Re(load4_ak . h_ak); iterations the
    !> largest GMRES count of the 2 ns solves (at most max_iterations each,
    !> default 400, else HP_NO_CONVERGENCE); solutions(:, b, j) the coupled
    !> vectors h on request.
    subroutine hp_energy_coefficients(sys, j0, k, q, node_r, iterations, status, &
            solutions, max_iterations)
        type(hp_energy_system_t), intent(inout) :: sys
        real(dp), intent(in) :: j0
        real(dp), intent(out) :: k(:, :, :), q(:, :, :)
        real(dp), allocatable, intent(out) :: node_r(:, :, :, :)
        integer, intent(out) :: iterations, status
        complex(dp), allocatable, intent(out), optional :: solutions(:, :, :)
        integer, intent(in), optional :: max_iterations
        complex(dp), allocatable :: b(:), u(:), h(:)
        real(dp), allocatable :: wt(:)
        real(dp) :: scale, residual
        integer :: sb, jf, a, kk, its, info, max_its, st

        status = HP_BAD_INPUT
        iterations = 0
        if (.not. allocated(sys%node)) return
        if (any(shape(k) /= [sys%ns, sys%ns, 2]) .or. any(shape(q) /= shape(k)) .or. &
            j0 <= 0) return
        max_its = default_iterations
        if (present(max_iterations)) max_its = max_iterations
        if (max_its < 1) return
        allocate (node_r(sys%nx, sys%ns, sys%ns, 2), b(sys%n_dof), u(sys%n_dof), &
            h(sys%n_dof), wt(sys%n_dof))
        if (present(solutions)) allocate (solutions(sys%n_dof, sys%ns, 2))
        ! GMRES minimises the residual of the Galerkin form (rows times their
        ! speed mass v_a d3_k): the lumped rows of nodes far in the Maxwellian
        ! tail are huge and would dominate the plain residual.
        do a = 1, sys%ns
            do kk = 1, sys%nx
                associate (nd => sys%node(kk, a))
                    wt(nd%off + 1:nd%off + nd%n) = sys%vt(a)*sys%d3(kk)
                end associate
            end do
        end do
        do sb = 1, sys%ns
            do jf = 1, 2
                b = 0
                do kk = 1, sys%nx
                    associate (nd => sys%node(kk, sb))
                        b(nd%off + 1:nd%off + nd%n) = sys%tz(sb)*sys%x(kk + 1)**(2*jf) &
                            *nd%terms%load3/nd%v
                    end associate
                end do
                b = wt*b
                ! Unit right-hand side: the GMRES target is absolute below one.
                scale = norm2(abs(b))
                u = 0
                h = 0
                if (scale > 0) then
                    b = b/scale
                    st = HP_OK
                    call complex_gmres_operator(operator, b, u, gmres_tolerance, max_its, &
                        min(gmres_restart, max_its), info, its, residual)
                    status = st
                    if (status /= HP_OK) return
                    iterations = max(iterations, its)
                    status = HP_NO_CONVERGENCE
                    if (info /= KRYLOV_OK) return
                    call precondition(sys, u/wt, h, status)
                    if (status /= HP_OK) return
                    h = scale*h
                end if
                if (present(solutions)) solutions(:, sb, jf) = h
                do a = 1, sys%ns
                    do kk = 1, sys%nx
                        associate (nd => sys%node(kk, a))
                            node_r(kk, a, sb, jf) = real(sum(nd%terms%load4 &
                                *h(nd%off + 1:nd%off + nd%n)), dp)
                        end associate
                    end do
                    k(a, sb, jf) = sys%tz(a)/(2*sqrt(pi)*j0)*sum(sys%d4*node_r(:, a, sb, jf))
                    q(a, sb, jf) = sys%tz(a)/(2*sqrt(pi)*j0)*sum(sys%d6*node_r(:, a, sb, jf))
                end do
            end do
        end do
        status = HP_OK
    contains
        subroutine operator(input, output)
            complex(dp), intent(in) :: input(:)
            complex(dp), intent(out) :: output(:)
            complex(dp), allocatable :: t(:)
            integer :: s
            allocate (t(size(input)))
            call precondition(sys, input/wt, t, s)
            if (s /= HP_OK) st = s
            call full_apply(sys, t, output)
            output = wt*output
        end subroutine operator
    end subroutine hp_energy_coefficients

    subroutine hp_energy_free(sys)
        type(hp_energy_system_t), intent(inout) :: sys
        integer :: i
        if (allocated(sys%fac)) then
            do i = 1, size(sys%fac)
                call hp_ntv_free(sys%fac(i))
            end do
            deallocate (sys%fac)
        end if
    end subroutine hp_energy_free

end module milca_dk_hp_energy
