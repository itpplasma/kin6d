program test_gll_block_preconditioner
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use omp_lib, only: omp_set_num_threads
    use kin6d_ideal_gll_spike
    use kin6d_landreman_benchmark, only: landreman_frozen_t
    implicit none
    real(dp), parameter :: shift = 0.01_dp, delta = 1.0_dp/64
    integer, parameter :: counts(3) = [1, 3, 16]
    type(landreman_frozen_t) :: field
    type(gll_map_t) :: map, negative
    type(gll_block_jacobi_t) :: pc, scalar, rejected, rational, threaded
    real(dp), allocatable :: x(:), gradient(:), matrix(:, :), v(:), w(:), &
        y(:), z(:), sum_y(:)
    real(dp) :: energy, r, zcoord, s, block_error, inverse_error, symmetry, value, best
    real(dp) :: lower(3, 3), local_rhs(3), local_solution(3), expected(3)
    integer :: status, nd, i, j, node, indices(3), e, a, i1, i2, q1, q2
    integer :: chosen_e, chosen_q, n, m, t

    field%delta = delta
    call build_gll_map(2, 1, delta, field, map, x, status)
    call require(status == 0, 'p2 structural fixture geometry')
    nd = size(x)
    allocate(gradient(nd), matrix(nd, nd), v(nd), w(nd), y(nd), z(nd), sum_y(nd))
    do i = 1, map%nodes
        r = x(i); zcoord = x(map%nodes + i)
        s = ((r*r - 1)**2/4 + zcoord**2)/delta
        if (.not. map%fixed(i)) x(i) = x(i) &
            + 1.e-4_dp*(1 - s)*cos(17*zcoord)
        if (.not. map%fixed(map%nodes + i)) x(map%nodes + i) = &
            x(map%nodes + i) + 1.e-4_dp*(1 - s)*sin(11*r)
        if (.not. map%fixed(2*map%nodes + i)) x(2*map%nodes + i) = &
            1.e-3_dp*(1 - s)*sin(13*r + 9*zcoord)
    end do
    call evaluate_gll_map(map, x, energy, gradient, status, .true.)
    call require(status == 0, 'structural fixture linearization')
    call build_gll_block_jacobi(map, shift, pc, status)
    call require(status == 0 .and. pc%ready, 'physical nodal block admission')
    call build_gll_block_jacobi(map, shift, scalar, status, diagonal_only=.true.)
    call require(status == 0 .and. scalar%ready, 'physical diagonal admission')
    ! Colour-parallel accumulation visits every node's cells in colour order for
    ! any thread count, so blocks and factors do not depend on the threads.
    do t = 1, 3
        call omp_set_num_threads(counts(t))
        call build_gll_block_jacobi(map, shift, threaded, status)
        call require(status == 0 .and. all(threaded%blocks == pc%blocks) &
            .and. all(threaded%factor == pc%factor), 'thread-independent blocks')
    end do
    ! An independent gather/transpose path constructs only this tiny oracle.
    do j = 1, nd
        v = 0; v(j) = 1
        call apply_gll_hessian(map, v, matrix(:, j))
        where (.not. map%fixed) matrix(:, j) = matrix(:, j) + shift*v
    end do
    block_error = 0
    do node = 1, map%nodes
        indices = [node, map%nodes + node, 2*map%nodes + node]
        block_error = max(block_error, maxval(abs(pc%blocks(:, :, node) &
            - matrix(indices, indices)))/max(1.0_dp, &
            maxval(abs(matrix(indices, indices)))))
        do j = 1, 3
            call require(abs(scalar%blocks(j, j, node) - matrix(indices(j), &
                indices(j))) <= 5.e-12_dp*max(1.0_dp, &
                abs(matrix(indices(j), indices(j)))), 'scalar diagonal oracle')
        end do
    end do
    call require(block_error <= 5.e-12_dp, 'HVP-column principal block oracle')
    do i = 1, nd
        v(i) = sin(0.37_dp*i); w(i) = cos(0.11_dp*i)
    end do
    call apply_gll_block_jacobi(pc, v, y, status)
    call require(status == 0, 'coupled inverse application')
    inverse_error = 0
    do node = 1, map%nodes
        indices = [node, map%nodes + node, 2*map%nodes + node]
        inverse_error = max(inverse_error, &
            norm2(matmul(matrix(indices, indices), y(indices)) - v(indices)))
    end do
    inverse_error = inverse_error/max(1.0_dp, norm2(v))
    call require(inverse_error <= 5.e-12_dp, 'independent block inverse residual')
    call apply_gll_block_jacobi(scalar, v, sum_y, status)
    call require(status == 0, 'scalar inverse application')
    do i = 1, nd
        call require(abs(matrix(i, i)*sum_y(i) - v(i)) <= 5.e-13_dp, &
            'independent scalar inverse')
    end do
    call require(maxval(abs(pack(y - v, map%fixed))) <= 1.e-13_dp, &
        'fixed row identity')
    call apply_gll_block_jacobi(pc, w, z, status)
    call require(status == 0, 'second coupled application')
    symmetry = abs(dot_product(v, z) - dot_product(w, y)) &
        /max(1.0_dp, norm2(v)*norm2(z), norm2(w)*norm2(y))
    call require(symmetry <= 5.e-13_dp, 'coupled inverse pairing symmetry')
    call require(dot_product(v, y) > 0, 'positive preconditioner action')
    call apply_gll_block_jacobi(pc, v + w, sum_y, status)
    call require(status == 0, 'coupled sum application')
    call require(norm2(sum_y - y - z) <= 5.e-13_dp*max(1.0_dp, norm2(sum_y)), &
        'coupled inverse linearity')
    v = 0
    where (map%fixed) v = 1
    call apply_gll_block_jacobi(pc, v, y, status)
    call require(status == 0 .and. maxval(abs(y - v)) <= 1.e-13_dp, &
        'fixed columns do not enter free response')

    ! Exact-rational coupled matrix with independently constructed Cholesky.
    rational%nodes = 1; rational%ready = .true.
    allocate(rational%factor(3, 3, 1))
    lower = 0
    lower(1, 1) = 2; lower(2, 1) = 0.5_dp; lower(3, 1) = 0.25_dp
    lower(2, 2) = sqrt(2.75_dp)
    lower(3, 2) = -0.375_dp/lower(2, 2)
    lower(3, 3) = sqrt(2 - lower(3, 1)**2 - lower(3, 2)**2)
    rational%factor(:, :, 1) = lower
    local_rhs = [0.875_dp, -1.4375_dp, 1.75_dp]
    expected = [0.25_dp, -0.5_dp, 0.75_dp]
    call apply_gll_block_jacobi(rational, local_rhs, local_solution, status)
    call require(status == 0 .and. &
        maxval(abs(local_solution - expected)) <= 1.e-13_dp, &
        'independent rational coupled inverse')

    ! Deliberately nonphysical cache: a free RR diagonal equals -1 at the
    ! unchanged shift.  Admission must fail without altering that operator.
    negative = map
    negative%hessian = 0
    node = findloc(map%fixed(:map%nodes), .false., dim=1)
    n = map%reference%ni(1); m = map%reference%no(1)
    best = 0; chosen_e = 0; chosen_q = 0
    do e = 1, map%reference%mesh%nelem
        do a = 1, n*n
            if (map%reference%mesh%dof(a, e) /= node) cycle
            i1 = mod(a - 1, n) + 1; i2 = (a - 1)/n + 1
            do q2 = 1, m
                do q1 = 1, m
                    value = map%reference%line(1)%b(q1, i1, 1) &
                        *map%reference%line(2)%b(q2, i2, 1)
                    if (abs(value) <= best) cycle
                    best = abs(value); chosen_e = e; chosen_q = q1 + m*(q2 - 1)
                end do
            end do
        end do
    end do
    call require(best > 0.1_dp, 'well-sized negative fixture basis')
    negative%hessian(1, 1, chosen_q, chosen_e) = -(shift + 1)/best**2
    call build_gll_block_jacobi(negative, shift, rejected, status)
    call require(status == 4 .and. .not. rejected%ready, &
        'negative block rejected without damping change')
    call apply_gll_block_jacobi(rejected, v, y, status)
    call require(status /= 0, 'rejected factor cannot be applied')
    negative%hessian(1, 1, chosen_q, chosen_e) = &
        ieee_value(0.0_dp, ieee_quiet_nan)
    call build_gll_block_jacobi(negative, shift, rejected, status)
    call require(status == 2 .and. .not. rejected%ready, 'nonfinite cache rejected')
    print '(a,3es17.9)', 'block_oracle_scaled_error,inverse_residual,symmetry=', &
        block_error, inverse_error, symmetry
    print '(a)', 'PASS: exact coupled blocks, constraints and failure semantics'
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (ok) return
        print '(a)', message
        error stop 'GLL coupled block oracle'
    end subroutine
end program
