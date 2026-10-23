module kin6d_nonnested_beltrami
    ! Taylor closure: <curl v,B> - mu <v,B> = <v,f>, B=curl A+h.
    ! Boundary tangential A fixes normal B. Only the n=0 gradient kernel
    ! needs a Coulomb multiplier; n>0 uses the exact R A_phi=0 gauge.
    use, intrinsic :: iso_fortran_env, only: dp => real64, i8 => int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_meridian_feec_n, only: meridian_load
    use fortfem_meridian_fourier_curl, only: meridian_fourier_curl_curl, &
        meridian_fourier_curl, &
        meridian_fourier_pairing, &
        meridian_fourier_diagonal_mass
    use fortsparse, only: csc_t, csc_matmul, csc_matvec, csc_transpose, &
        csc_from_triplet, csc_normalize_symmetric_raw, &
        sparse_solver_t, &
        sparse_factor, sparse_solve, &
        sparse_destroy, sparse_inertia, sparse_factor_stats, &
        sparse_factor_stats_t, &
        fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_NATIVE, &
        FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_INERTIA_RELIABLE, &
        FORTSPARSE_ORDER_AMD, native_ordering_available, &
        FORTSPARSE_SINGULAR
    use fortnum_eigen, only: lanczos, eigen_options_t, eigen_stats_t, &
        EIGEN_NEAREST, EIGEN_OK
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_divergence
    use kin6d_nonnested_mapped, only: nn_map_t, map_metric_t, &
        METRIC_IDENTITY, nn_harmonic_forms
    use fortnum_eigen, only: eigen_real_operator_dp_t
    use fortnum_eigen, only: eigen_real_solve_dp_t
    use fortfem_meridian_feec_n, only: meridian_field_t
    implicit none
    private
    public :: beltrami_system_t, beltrami_report_t, beltrami_prepare, beltrami_solve
    public :: beltrami_constraint, BELTRAMI_OK, BELTRAMI_INVALID, BELTRAMI_RESONANT
    integer, parameter :: BELTRAMI_OK = 0, BELTRAMI_INVALID = 1, BELTRAMI_RESONANT = 2
    real(dp), parameter :: pi = acos(-1.0_dp)
    type :: beltrami_system_t
        private
        type(csc_t) :: energy, helicity, curl, gauge
        real(dp), allocatable :: energy_h(:), helicity_h(:), lift(:)
        integer, allocatable :: index(:), gauge_index(:)
        integer :: nf = 0, ng = 0, modes = 0, n1 = 0, n2 = 0, nfp = 1
        real(dp) :: harmonic = 0, flux = 0, h_energy = 0
        logical :: ready = .false., closed = .false.
    end type
    type :: beltrami_report_t
        integer :: status = BELTRAMI_INVALID, sparse_status = 0
        integer(i8) :: negative = -1, zero = -1, positive = -1
        integer :: inertia_status = -1, nfp = 1
        real(dp) :: residual = 0, divergence = 0, helicity = 0, current = 0, energy = 0
        real(dp) :: min_pivot = -1, pivot_tolerance = -1, mu = 0
        real(dp) :: nearest_eigenvalue = huge(1.0_dp), eigen_residual = huge(1.0_dp)
        character(24) :: claim = 'candidate', bound_status = 'unavailable'
        character(160) :: message = 'invalid or unprepared Beltrami problem'
    end type
    type, extends(eigen_real_operator_dp_t) :: beltrami_solve_apply_operator_t
        type(csc_t), pointer :: a => null()
    contains
        procedure :: apply => beltrami_solve_apply_operator_evaluate
    end type
    type, extends(eigen_real_solve_dp_t) :: beltrami_solve_inverse_operator_t
        type(sparse_solver_t), pointer :: solver => null()
    contains
        procedure :: solve => beltrami_solve_inverse_operator_evaluate
    end type
    type, extends(meridian_field_t) :: beltrami_prepare_harmonic_unit_t
        type(free_b_field_t), pointer :: field => null()
    contains
        procedure :: evaluate => beltrami_prepare_harmonic_unit_evaluate
    end type
contains
    subroutine beltrami_prepare(field, flux, data, nphi, system, status, map, nfp)
        type(free_b_field_t), intent(in), target :: field
        real(dp), intent(in) :: flux, data(:)
        integer, intent(in) :: nphi
        type(beltrami_system_t), intent(out) :: system
        integer, intent(out) :: status
        type(nn_map_t), optional, intent(in), target :: map
        integer, optional, intent(in) :: nfp
        type(nn_map_t), target :: identity
        type(map_metric_t) :: metric
        type(csc_t) :: pair
        type(fortsparse_status_t) :: fst
        logical, allocatable :: fixed(:)
        real(dp), allocatable :: g(:)
        complex(dp), allocatable :: load(:)
        integer :: n, b, st, j, harmonic
        real(dp) :: harmonic_energy, harmonic_volume
        type(beltrami_prepare_harmonic_unit_t) :: harmonic_unit_context
        status = BELTRAMI_INVALID
        if (.not. allocated(field%b)) return
        if (.not. ieee_is_finite(flux)) return
        if (.not. all(ieee_is_finite(data))) return
        if (.not. ieee_is_finite(field%pprime)) return
        if (abs(field%pprime) > 0) return
        if (allocated(field%p)) then
            if (size(field%p) /= field%fs%n0*(2*ubound(field%b, 2) + 1)) return
            if (.not. all(ieee_is_finite(field%p))) return
            if (any(abs(field%p(:field%fs%n0) - field%p(1)) > 0)) return
            if (size(field%p) > field%fs%n0) then
                if (any(abs(field%p(field%fs%n0 + 1:)) > 0)) return
            end if
        end if
        ! One solid torus: extra meridional holes need additional periods.
        if (field%mesh%tlast(1) == 0) return
        if (present(nfp)) system%nfp = nfp
        if (system%nfp < 1) return
        system%modes = ubound(field%b, 2)
        system%n1 = field%fs%sp(1)%ndof; system%n2 = field%fs%sp(2)%ndof
        n = (2*system%modes + 1)*system%n1
        if (size(data) /= n .or. nphi < 2*system%modes + 1) return
        metric%fs => field%fs
        if (present(map)) then
            if (.not. allocated(map%d)) return
            if (size(map%d, 1) /= field%fs%n0 .or. size(map%d, 2) /= 2) return
            if (lbound(map%d, 3) /= 0 .or. ubound(map%d, 3) /= 2*map%modes) return
            if (field%fs%cyl_r0 > 0) return
            do harmonic = 1, map%modes
                if (mod(harmonic, system%nfp) == 0) cycle
                if (maxval(abs(map%d(:, :, 2*harmonic - 1:2*harmonic))) &
                    > 1.0e-13_dp*max(1.0_dp, maxval(abs(map%d)))) return
            end do
            metric%map => map
        else
            identity%modes = 0
            allocate (identity%d(field%fs%n0, 2, 0:0)); identity%d = 0
            metric%map => identity
            metric%kind = METRIC_IDENTITY
        end if
        call meridian_fourier_curl(field%fs, system%modes, system%curl, st)
        if (st /= 0) return
        if (present(map)) then
            call meridian_fourier_curl_curl(field%fs, system%modes, nphi, metric, &
                system%energy, st)
        else
            call meridian_fourier_diagonal_mass(field%fs, 2, system%modes, pair, st)
            if (st /= 0) return
            call csc_matmul(pair, system%curl, system%energy, fst)
            if (.not. status_ok(fst)) return
            call csc_transpose(system%curl, pair, fst)
            if (.not. status_ok(fst)) return
            block
                type(csc_t) :: product
                call csc_matmul(pair, system%energy, product, fst)
                if (.not. status_ok(fst)) return
                system%energy = product
            end block
            st = 0
        end if
        if (st /= 0) return
        call meridian_fourier_pairing(field%fs, system%modes, pair, st)
        if (st /= 0) return
        call csc_matmul(pair, system%curl, system%helicity, fst)
        if (.not. status_ok(fst)) return
        call csc_matmul(field%fs%mass(1), field%fs%grad%a, system%gauge, fst)
        if (.not. status_ok(fst)) return
        allocate (fixed(n), system%index(n), system%lift(n), &
            system%gauge_index(field%fs%n0))
        do b = 0, 2*system%modes
            fixed(b*system%n1 + 1:(b + 1)*system%n1) = field%fs%sp(1)%bnd
            if (b > 0) then
                if (mod((b + 1)/2, system%nfp) /= 0) then
                    if (any(abs(data(b*system%n1 + 1:(b + 1)*system%n1)) > 0)) return
                    fixed(b*system%n1 + 1:(b + 1)*system%n1) = .true.
                end if
            end if
            if (b > 0) then
                if (any(abs(data(b*system%n1 + field%fs%nn + 1: &
                    (b + 1)*system%n1)) > 0)) return
                fixed(b*system%n1 + field%fs%nn + 1:(b + 1)*system%n1) = .true.
            end if
        end do
        system%closed = .not. any(abs(pack(data, fixed)) > 0)
        system%lift = merge(data, 0.0_dp, fixed)
        system%nf = 0; system%index = 0
        do j = 1, n
            if (fixed(j)) cycle
            system%nf = system%nf + 1; system%index(j) = system%nf
        end do
        system%ng = 0; system%gauge_index = 0
        do j = 1, field%fs%n0
            if (field%fs%sp(0)%bnd(j)) cycle
            system%ng = system%ng + 1; system%gauge_index(j) = system%nf + system%ng
        end do
        call nn_harmonic_forms(field, metric, nphi, g, harmonic_energy, &
            harmonic_volume, st)
        if (st /= 0) return
        system%flux = flux
        harmonic_unit_context%field => field
        call meridian_load(field%fs, 2, harmonic_unit_context, load, st)
        if (st /= 0) return
        block
            real(dp), allocatable :: boundary_b(:)
            boundary_b = csc_matvec(system%curl, system%lift)
            system%harmonic = (flux - dot_product(real(load, dp), &
                boundary_b(:system%n2)))/harmonic_volume
        end block
        if (.not. ieee_is_finite(system%harmonic)) return
        system%h_energy = system%harmonic**2*harmonic_energy
        call csc_transpose(system%curl, pair, fst)
        if (.not. status_ok(fst)) return
        system%energy_h = system%harmonic*csc_matvec(pair, g)
        harmonic_unit_context%field => field
        call meridian_load(field%fs, 1, harmonic_unit_context, load, st)
        if (st /= 0) return
        allocate (system%helicity_h(n)); system%helicity_h = 0
        system%helicity_h(:system%n1) = system%harmonic*real(load, dp)
        system%ready = .true.; status = BELTRAMI_OK
    end subroutine
    subroutine beltrami_solve(system, mu, field, report, source, potential, dmu)
        type(beltrami_system_t), intent(in) :: system
        real(dp), intent(in) :: mu
        type(free_b_field_t), intent(inout) :: field
        type(beltrami_report_t), intent(out) :: report
        ! source is the weak <v,f> load in the full real Fourier V1 layout.
        real(dp), optional, intent(in) :: source(:)
        real(dp), allocatable, optional, intent(out) :: potential(:), dmu(:)
        type(csc_t), target :: a
        type(sparse_solver_t), target :: solver
        type(fortsparse_status_t) :: fst
        type(sparse_factor_stats_t) :: stats
        type(eigen_options_t) :: eigen_options
        type(eigen_stats_t) :: eigen_stats
        real(dp), allocatable :: eigenvector(:, :)
        real(dp) :: eigenvalue(1), eigen_error(1), operator_scale
        integer :: eigen_status
        integer, allocatable :: row(:), col(:)
        real(dp), allocatable :: val(:), rhs(:), x(:), full(:), r(:), bh(:)
        real(dp), allocatable :: load(:), dx(:)
        real(dp) :: correction
        integer :: cap, k, j, p, i, ii, jj, n, st
        type(beltrami_solve_apply_operator_t) :: apply_operator_context
        type(beltrami_solve_inverse_operator_t) :: inverse_operator_context
        report%mu = mu; report%nfp = system%nfp
        if (.not. system%ready .or. .not. ieee_is_finite(mu)) return
        if (.not. allocated(field%b)) return
        if (size(field%b, 1) /= system%n2) return
        if (ubound(field%b, 2) /= system%modes) return
        if (lbound(field%b, 2) /= 0 .or. field%fs%sp(1)%ndof /= system%n1) return
        n = system%nf + system%ng
        allocate (load(size(system%lift)))
        load = -system%energy_h + mu*system%helicity_h
        if (present(source)) then
            if (size(source) /= size(load)) return
            if (.not. all(ieee_is_finite(source))) return
            do j = 1, system%modes
                if (mod(j, system%nfp) == 0) cycle
                if (any(abs(source((2*j - 1)*system%n1 + 1: &
                    (2*j + 1)*system%n1)) > 0)) return
            end do
            load = load + source
        end if
        load = load - csc_matvec(system%energy, system%lift) &
            + mu*csc_matvec(system%helicity, system%lift)
        cap = 2*(system%energy%nnz + system%helicity%nnz) + 2*system%gauge%nnz
        allocate (row(cap), col(cap), val(cap), rhs(n), x(n)); k = 0; rhs = 0
        do j = 1, size(load)
            jj = system%index(j)
            if (jj > 0) rhs(jj) = load(j)
        end do
        call put(system%energy, 1.0_dp)
        call put(system%helicity, -mu)
        do j = 1, system%gauge%ncol
            jj = system%gauge_index(j)
            if (jj == 0) cycle
            do p = system%gauge%col_ptr(j), system%gauge%col_ptr(j + 1) - 1
                i = system%gauge%row_idx(p); ii = system%index(i)
                rhs(jj) = rhs(jj) - system%gauge%val(p)*system%lift(i)
                if (ii == 0) cycle
                k = k + 1; row(k) = ii; col(k) = jj; val(k) = system%gauge%val(p)
                k = k + 1; row(k) = jj; col(k) = ii; val(k) = system%gauge%val(p)
            end do
        end do
        call csc_from_triplet(n, n, row(:k), col(:k), val(:k), a, fst)
        if (.not. status_ok(fst)) return
        call csc_normalize_symmetric_raw(n, a%col_ptr, a%row_idx, a%val, &
            1.0e-12_dp, correction, fst)
        if (.not. status_ok(fst)) then
            report%message = fst%msg
            return
        end if
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%refine = .true.
        solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        solver%options%ldlt_blocked = .true.
        if (native_ordering_available(FORTSPARSE_ORDER_AMD)) &
            solver%options%ordering = FORTSPARSE_ORDER_AMD
        call sparse_factor(solver, a, fst)
        report%sparse_status = fst%code; report%message = fst%msg
        if (.not. status_ok(fst)) then
            if (fst%code == FORTSPARSE_SINGULAR) report%status = BELTRAMI_RESONANT
            call sparse_destroy(solver)
            return
        end if
        call sparse_inertia(solver, report%negative, report%zero, report%positive, &
            report%inertia_status, fst)
        call sparse_factor_stats(solver, stats, fst)
        report%min_pivot = stats%inertia_min_pivot
        report%pivot_tolerance = stats%inertia_tolerance
        ! Pivot diagnostics are numerical reliability checks, not a spectral
        ! gap certificate. Reject uncertain/near-zero pivots without shifting.
        if (report%inertia_status /= FORTSPARSE_INERTIA_RELIABLE .or. &
            report%zero > 0 .or. report%min_pivot <= 1000*report%pivot_tolerance) then
            report%status = BELTRAMI_RESONANT
            report%message = 'singular or numerically near-resonant Beltrami operator'
            call sparse_destroy(solver)
            return
        end if
        call sparse_solve(solver, rhs, x, fst)
        if (.not. status_ok(fst)) then
            report%message = fst%msg
            call sparse_destroy(solver)
            return
        end if
        r = csc_matvec(a, x) - rhs
        report%residual = norm2(r)/max(1.0_dp, norm2(rhs))
        if (report%residual > 1.0e-9_dp) then
            ! A residual failure near resonance needs a spectral diagnostic,
            ! not an invented field or a shifted operator. FortNum applies
            ! shift-invert Lanczos to this same unshifted LDLT factor.
            operator_scale = 0
            do j = 1, a%ncol
                operator_scale = max(operator_scale, &
                    sum(abs(a%val(a%col_ptr(j):a%col_ptr(j + 1) - 1))))
            end do
            eigen_options%which = EIGEN_NEAREST
            eigen_options%max_basis = min(24, n); eigen_options%max_restarts = 4
            eigen_options%operator_scale = operator_scale
            allocate (eigenvector(n, 1))
            eigenvalue = huge(1.0_dp); eigen_error = huge(1.0_dp)
        apply_operator_context%a => a
        inverse_operator_context%solver => solver
            call lanczos(apply_operator_context, x, eigenvalue, eigenvector, &
                eigen_error, &
                1.0e-10_dp, eigen_options, eigen_status, eigen_stats, &
                solve_shift=inverse_operator_context, sigma=0.0_dp)
            report%nearest_eigenvalue = eigenvalue(1)
            report%eigen_residual = eigen_error(1)
            report%message = &
                'equation residual failed; spectral diagnostic inconclusive'
            if (eigen_status == EIGEN_OK) then
                if (abs(eigenvalue(1)) + eigen_error(1) &
                    <= sqrt(epsilon(1.0_dp))*operator_scale) then
                    report%status = BELTRAMI_RESONANT
                    report%message = 'near-singular Beltrami operator detected by '// &
                        'FortNum eigen solve'
                end if
            end if
            call sparse_destroy(solver)
            return
        end if
        full = system%lift
        do j = 1, size(full)
            if (system%index(j) > 0) full(j) = full(j) + x(system%index(j))
        end do
        if (present(dmu)) then
            load = csc_matvec(system%helicity, full) + system%helicity_h
            rhs = 0
            do j = 1, size(full)
                if (system%index(j) > 0) rhs(system%index(j)) = load(j)
            end do
            allocate (dx(n))
            call sparse_solve(solver, rhs, dx, fst)
            if (.not. status_ok(fst)) then
                call sparse_destroy(solver)
                return
            end if
            allocate (dmu(size(full))); dmu = 0
            do j = 1, size(full)
                if (system%index(j) > 0) dmu(j) = dx(system%index(j))
            end do
        end if
        call sparse_destroy(solver)
        if (report%residual > 1.0e-9_dp) return
        report%helicity = 2*pi*(dot_product(full, csc_matvec(system%helicity, full)) &
            + 2*dot_product(full, system%helicity_h))
        report%energy = pi*(dot_product(full, csc_matvec(system%energy, full)) &
            + 2*dot_product(full, system%energy_h) + system%h_energy)/field%mu0
        report%current = mu*system%flux/field%mu0
        if (present(source)) report%current = report%current &
            + sum(source(system%n1 - size(system%gauge_index) + 1: &
            system%n1))/field%mu0
        bh = csc_matvec(system%curl, full)
        field%b(:, 0) = cmplx(bh(:system%n2), 0.0_dp, dp)
        do j = 1, system%modes
            field%b(:, j) = cmplx(bh((2*j - 1)*system%n2 + 1:2*j*system%n2), &
                -bh(2*j*system%n2 + 1:(2*j + 1)*system%n2), dp)/2
        end do
        field%psi = full(field%fs%nn + 1:system%n1)
        field%f_vac = system%harmonic
        call nn_divergence(field, report%divergence, st)
        if (st /= 0) return
        if (present(potential)) potential = full
        report%status = BELTRAMI_OK
        report%message = 'Beltrami candidate; bounds unavailable'
    contains
        subroutine put(m, scale)
            type(csc_t), intent(in) :: m
            real(dp), intent(in) :: scale
            integer :: j, p, ii, jj
            do j = 1, m%ncol
                jj = system%index(j)
                if (jj == 0) cycle
                do p = m%col_ptr(j), m%col_ptr(j + 1) - 1
                    ii = system%index(m%row_idx(p))
                    if (ii == 0) cycle
                    k = k + 1; row(k) = ii; col(k) = jj; val(k) = scale*m%val(p)
                    k = k + 1; row(k) = jj; col(k) = ii; val(k) = 0
                end do
            end do
        end subroutine
    end subroutine
    subroutine beltrami_constraint(system, kind, target, bracket, tolerance, field, &
        report)
        type(beltrami_system_t), intent(in) :: system
        character(*), intent(in) :: kind
        real(dp), intent(in) :: target, bracket(2), tolerance
        type(free_b_field_t), intent(inout) :: field
        type(beltrami_report_t), intent(out) :: report
        type(free_b_field_t) :: candidate
        type(beltrami_report_t) :: left_report, right_report
        real(dp) :: lo, hi, mu, fl, fh, fm
        integer :: j
        if (.not. system%ready) return
        if (.not. system%closed) return ! Relative helicity/current: closed wall.
        if (.not. all(ieee_is_finite([target, bracket, tolerance]))) return
        if (tolerance <= 0 .or. bracket(1) >= bracket(2)) return
        if (kind == 'current') then
            if (abs(system%flux) <= tiny(1.0_dp)) return
            mu = target*field%mu0/system%flux
            if (mu < bracket(1) .or. mu > bracket(2)) return
            call beltrami_solve(system, mu, field, report)
            return
        end if
        if (kind /= 'helicity') return
        ! Restrict a helicity bracket to one side of zero. On each side
        ! the inertia of the positive-energy/helicity pencil is monotone,
        ! so equal endpoint inertia excludes an intervening resonance.
        if (bracket(1) < 0 .and. bracket(2) > 0) return
        candidate = field; lo = bracket(1); hi = bracket(2)
        call beltrami_solve(system, lo, candidate, left_report)
        if (left_report%status /= BELTRAMI_OK) then
            report = left_report
            return
        end if
        call beltrami_solve(system, hi, candidate, right_report)
        if (right_report%status /= BELTRAMI_OK) then
            report = right_report
            return
        end if
        if (left_report%negative /= right_report%negative) then
            report%status = BELTRAMI_RESONANT
            report%message = 'constraint bracket crosses a Beltrami resonance'
            return
        end if
        fl = left_report%helicity - target; fh = right_report%helicity - target
        if (fl*fh > 0) return
        do j = 1, 80
            mu = (lo + hi)/2
            call beltrami_solve(system, mu, candidate, report)
            if (report%status /= 0) return
            if (report%negative /= left_report%negative) then
                report%status = BELTRAMI_RESONANT
                report%message = 'constraint iteration changed Beltrami branch'
                return
            end if
            fm = report%helicity - target
            if (abs(fm) <= tolerance*max(1.0_dp, abs(target))) then
                field = candidate
                return
            end if
            if (fl*fm <= 0) then
                hi = mu
            else
                lo = mu; fl = fm
            end if
        end do
        report%status = BELTRAMI_INVALID
        report%message = 'helicity constraint did not converge'
    end subroutine
    subroutine beltrami_solve_apply_operator_evaluate(self, x, y)
        class(beltrami_solve_apply_operator_t), intent(inout) :: self
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: y(:)
            y = csc_matvec(self%a, x)
    end subroutine beltrami_solve_apply_operator_evaluate
    subroutine beltrami_solve_inverse_operator_evaluate(self, x, y, info)
        class(beltrami_solve_inverse_operator_t), intent(inout) :: self
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: y(:)
            integer, intent(out) :: info
            type(fortsparse_status_t) :: fst
            call sparse_solve(self%solver, x, y, fst)
            info = fst%code
    end subroutine beltrami_solve_inverse_operator_evaluate
    subroutine beltrami_prepare_harmonic_unit_evaluate(self, x, v)
        class(beltrami_prepare_harmonic_unit_t), intent(in) :: self
            real(dp), intent(in) :: x(2)
            complex(dp), intent(out) :: v(3)
            real(dp) :: radius
            radius = x(1)
            if (self%field%fs%cyl_r0 > 0) radius = self%field%fs%cyl_r0
            v = 0; v(2) = 1/radius
    end subroutine beltrami_prepare_harmonic_unit_evaluate
end module
