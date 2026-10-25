module kin6d_ideal_gll_spike
    ! Opt-in composition of FortFEM tensor operations and KIN6D physics leaves.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_tensor_mesh, only: tensor_mesh_t, ogrid_tensor_mesh, build_tensor_mesh
    use fortfem_tensor_operator, only: tensor_operator_t, build_tensor_operator
    use fortfem_tensor_contractions, only: tensor_contract
    use kin6d_ideal_map, only: ideal_frozen_t
    use kin6d_ideal_map_leaves, only: ideal_reference_form
    use kin6d_ideal_gll_leaves, only: gll_energy, gll_gradient, gll_hessian, &
        gll_landreman_map
    implicit none
    private
    public :: gll_map_t, build_gll_map, evaluate_gll_map, apply_gll_hessian
    public :: diffusion, regularized_diffusion
    public :: gll_block_jacobi_t, build_gll_block_jacobi, apply_gll_block_jacobi
    integer, parameter :: component(7) = [1, 1, 2, 3, 1, 2, 3]
    integer, parameter :: derivative(7) = [0, 1, 1, 1, 2, 2, 2]
    type :: gll_map_t
        type(tensor_operator_t) :: reference
        real(dp), allocatable :: flux(:, :, :), pressure(:, :), hessian(:, :, :, :)
        logical, allocatable :: fixed(:)
        integer :: nodes = 0, pin = 0
    end type
    type :: gll_block_jacobi_t
        real(dp), allocatable :: blocks(:, :, :), factor(:, :, :)
        integer :: nodes = 0
        logical :: ready = .false.
        logical :: diagonal_only = .false.
    end type
    interface
        subroutine dpotrf(uplo, n, a, lda, info)
            import dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
        end subroutine
        subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
            import dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine
    end interface
contains
    subroutine diffusion(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        integer :: i
        k = 0; c = 0
        do i = 1, size(x)
            k(i, i) = 1
        end do
    end subroutine

    subroutine regularized_diffusion(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        call diffusion(x, k, c)
        c = 1
    end subroutine

    subroutine build_gll_map(p, n, delta, field, map, x, status)
        integer, intent(in) :: p, n
        real(dp), intent(in) :: delta
        class(ideal_frozen_t), intent(in) :: field
        type(gll_map_t), intent(out) :: map
        real(dp), allocatable, intent(out) :: x(:)
        integer, intent(out) :: status
        type(tensor_mesh_t) :: disk, mesh
        real(dp), allocatable :: geometry(:, :, :)
        real(dp) :: chi, grad(2), b3, pressure, dg(2, 2)
        integer :: a, e, q, nd
        call ogrid_tensor_mesh(2, n, p, disk, status)
        if (status /= 0) return
        geometry = disk%geometry
        do e = 1, disk%nelem
            do a = 1, size(geometry, 2)
                call gll_landreman_map(disk%geometry(1, a, e), &
                    disk%geometry(2, a, e), delta, geometry(1, a, e), geometry(2, a, e))
            end do
        end do
        call build_tensor_mesh(2, p, p, disk%vertices, geometry, mesh, status)
        if (status /= 0) return
        mesh%boundary = disk%boundary
        ! Gauss overintegration: nq=2p+3 per direction; not GLL collocation.
        call build_tensor_operator(mesh, diffusion, map%reference, status, nq=2*p + 3)
        if (status /= 0) return
        nd = mesh%ndof
        map%nodes = nd
        allocate(x(3*nd), map%fixed(3*nd))
        x(1:nd) = mesh%coordinates(1, :)
        x(nd + 1:2*nd) = mesh%coordinates(2, :)
        x(2*nd + 1:) = 0
        map%fixed = .false.
        map%fixed(:nd) = mesh%boundary
        map%fixed(nd + 1:2*nd) = mesh%boundary
        map%pin = minloc(sum((mesh%coordinates - &
            spread([1.0_dp, 0.0_dp], 2, nd))**2, dim=1), dim=1)
        map%fixed(2*nd + map%pin) = .true.
        q = product(map%reference%no)
        allocate(map%flux(3, q, mesh%nelem), map%pressure(q, mesh%nelem), &
            map%hessian(7, 7, q, mesh%nelem))
        !$omp parallel do default(shared) private(e,q,chi,grad,b3,pressure,dg)
        do e = 1, mesh%nelem
            do q = 1, product(map%reference%no)
                call field%evaluate(map%reference%points(:, q, e), chi, grad, &
                    b3, pressure)
                dg = map%reference%jac(:2, :2, q, e)
                call ideal_reference_form(dg(1, 1), dg(2, 1), dg(1, 2), dg(2, 2), &
                    grad(1), grad(2), b3, map%flux(1, q, e), map%flux(2, q, e), &
                    map%flux(3, q, e))
                map%pressure(q, e) = pressure
            end do
        end do
        !$omp end parallel do
    end subroutine

    subroutine evaluate_gll_map(map, x, energy, gradient, status, linearize)
        type(gll_map_t), intent(inout) :: map
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: energy, gradient(:)
        integer, intent(out) :: status
        logical, intent(in) :: linearize
        real(dp) :: cell_energy(map%reference%mesh%nelem), &
            cell_minjac(map%reference%mesh%nelem)
        integer :: colour, ce, e, q, c, k, i, j, n, m, nl, nq, nd, ws
        real(dp), parameter :: twopi = 2*acos(-1.0_dp)
        n = map%reference%ni(1); m = map%reference%no(1)
        nl = n*n; nq = m*m; nd = map%nodes; ws = max(n, m)**3
        gradient = 0; cell_energy = 0; cell_minjac = huge(1.0_dp)
        !$omp parallel default(shared) private(colour,ce,e,q,c,k,i,j)
        block
            real(dp), allocatable :: z(:, :), g(:, :), u(:), v(:), t1(:), t2(:)
            real(dp) :: val, jac, gg(7), hh(28), b(3), weight
            allocate(z(nq, 7), g(nq, 7), u(nl), v(nl), t1(ws), t2(ws))
            do colour = 1, map%reference%ncolours
                !$omp do schedule(static)
                do ce = map%reference%colour_ptr(colour), &
                    map%reference%colour_ptr(colour + 1) - 1
                    e = map%reference%colour_cells(ce)
                    do c = 1, 7
                        u = x((component(c) - 1)*nd + map%reference%mesh%dof(:, e))
                        call channel(map, c, u, z(:, c), t1, t2, .false.)
                    end do
                    do q = 1, nq
                        b = map%flux(:, q, e)
                        call gll_energy(z(q, 1), z(q, 2), z(q, 3), z(q, 4), &
                            z(q, 5), z(q, 6), z(q, 7), b(1), b(2), b(3), &
                            map%pressure(q, e), 1.0_dp, val, jac)
                        cell_minjac(e) = min(cell_minjac(e), jac)
                        weight = twopi*map%reference%weight(q, e)
                        cell_energy(e) = cell_energy(e) + weight*val
                        call gll_gradient(z(q, 1), z(q, 2), z(q, 3), z(q, 4), &
                            z(q, 5), z(q, 6), z(q, 7), b(1), b(2), b(3), &
                            map%pressure(q, e), 1.0_dp, &
                            gg(1), gg(2), gg(3), gg(4), gg(5), gg(6), gg(7))
                        g(q, :) = weight*gg
                        if (.not. linearize) cycle
                        call gll_hessian(z(q, 1), z(q, 2), z(q, 3), z(q, 4), &
                            z(q, 5), z(q, 6), z(q, 7), b(1), b(2), b(3), &
                            map%pressure(q, e), 1.0_dp, &
                            hh(1), hh(2), hh(3), hh(4), hh(5), hh(6), hh(7), &
                            hh(8), hh(9), hh(10), hh(11), hh(12), hh(13), hh(14), &
                            hh(15), hh(16), hh(17), hh(18), hh(19), hh(20), hh(21), &
                            hh(22), hh(23), hh(24), hh(25), hh(26), hh(27), hh(28))
                        k = 0
                        do j = 1, 7
                            do i = 1, j
                                k = k + 1
                                map%hessian(i, j, q, e) = weight*hh(k)
                                map%hessian(j, i, q, e) = weight*hh(k)
                            end do
                        end do
                    end do
                    do c = 1, 7
                        call channel(map, c, g(:, c), v, t1, t2, .true.)
                        do i = 1, nl
                            j = (component(c) - 1)*nd + map%reference%mesh%dof(i, e)
                            gradient(j) = gradient(j) + v(i)
                        end do
                    end do
                end do
                !$omp end do
            end do
        end block
        !$omp end parallel
        energy = sum(cell_energy)
        status = 0
        if (minval(cell_minjac) <= 0 .or. .not. ieee_is_finite(energy)) status = 1
        where (map%fixed) gradient = 0
    end subroutine

    subroutine apply_gll_hessian(map, x, y)
        type(gll_map_t), intent(in) :: map
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer :: colour, ce, e, c, i, j, q, n, m, nl, nq, nd, ws
        n = map%reference%ni(1); m = map%reference%no(1)
        nl = n*n; nq = m*m; nd = map%nodes; ws = max(n, m)**3
        y = 0
        !$omp parallel default(shared) private(colour,ce,e,c,i,j,q)
        block
            real(dp), allocatable :: z(:, :), f(:, :), u(:), v(:), t1(:), t2(:)
            allocate(z(nq, 7), f(nq, 7), u(nl), v(nl), t1(ws), t2(ws))
            do colour = 1, map%reference%ncolours
                !$omp do schedule(static)
                do ce = map%reference%colour_ptr(colour), &
                    map%reference%colour_ptr(colour + 1) - 1
                    e = map%reference%colour_cells(ce)
                    do c = 1, 7
                        do i = 1, nl
                            j = (component(c) - 1)*nd + map%reference%mesh%dof(i, e)
                            u(i) = x(j)
                            if (map%fixed(j)) u(i) = 0
                        end do
                        call channel(map, c, u, z(:, c), t1, t2, .false.)
                    end do
                    do q = 1, nq
                        f(q, :) = matmul(map%hessian(:, :, q, e), z(q, :))
                    end do
                    do c = 1, 7
                        call channel(map, c, f(:, c), v, t1, t2, .true.)
                        do i = 1, nl
                            j = (component(c) - 1)*nd + map%reference%mesh%dof(i, e)
                            y(j) = y(j) + v(i)
                        end do
                    end do
                end do
                !$omp end do
            end do
        end block
        !$omp end parallel
        where (map%fixed) y = x
    end subroutine

    subroutine build_gll_block_jacobi(map, shift, pc, status, diagonal_only)
        ! Exact node-principal blocks of P H P + shift*P + Q.  This is a
        ! local coupled preconditioner, not the production gauge Schur solve.
        type(gll_map_t), intent(in) :: map
        real(dp), intent(in) :: shift
        type(gll_block_jacobi_t), intent(out) :: pc
        integer, intent(out) :: status
        logical, intent(in), optional :: diagonal_only
        real(dp) :: jet(7, 3), work(7, 3), value, ds, dt
        integer :: n, m, e, a, i1, i2, q1, q2, q, node, i, j, info
        logical :: fixed(3)
        status = 1
        if (.not. ieee_is_finite(shift) .or. shift <= 0) return
        if (.not. allocated(map%hessian) .or. .not. allocated(map%fixed)) return
        if (map%nodes < 1 .or. size(map%fixed) /= 3*map%nodes) return
        if (.not. all(ieee_is_finite(map%hessian))) then
            status = 2
            return
        end if
        n = map%reference%ni(1); m = map%reference%no(1)
        pc%nodes = map%nodes
        if (present(diagonal_only)) pc%diagonal_only = diagonal_only
        allocate(pc%blocks(3, 3, pc%nodes))
        if (.not. pc%diagonal_only) allocate(pc%factor(3, 3, pc%nodes))
        pc%blocks = 0
        do e = 1, map%reference%mesh%nelem
            do q = 1, m*m
                if (maxval(abs(map%hessian(:, :, q, e) - &
                    transpose(map%hessian(:, :, q, e)))) > &
                    128*epsilon(1.0_dp)*max(1.0_dp, &
                    maxval(abs(map%hessian(:, :, q, e))))) then
                    status = 3
                    return
                end if
            end do
            do i2 = 1, n
                do i1 = 1, n
                    a = i1 + n*(i2 - 1)
                    node = map%reference%mesh%dof(a, e)
                    do q2 = 1, m
                        do q1 = 1, m
                            q = q1 + m*(q2 - 1)
                            value = map%reference%line(1)%b(q1, i1, 1) &
                                *map%reference%line(2)%b(q2, i2, 1)
                            ds = map%reference%line(1)%d(q1, i1, 1) &
                                *map%reference%line(2)%b(q2, i2, 1)
                            dt = map%reference%line(1)%b(q1, i1, 1) &
                                *map%reference%line(2)%d(q2, i2, 1)
                            jet = 0
                            jet(:, 1) = [value, ds, 0.0_dp, 0.0_dp, &
                                dt, 0.0_dp, 0.0_dp]
                            jet(:, 2) = [0.0_dp, 0.0_dp, ds, 0.0_dp, &
                                0.0_dp, dt, 0.0_dp]
                            jet(:, 3) = [0.0_dp, 0.0_dp, 0.0_dp, ds, &
                                0.0_dp, 0.0_dp, dt]
                            work = matmul(map%hessian(:, :, q, e), jet)
                            do j = 1, 3
                                do i = 1, j
                                    if (pc%diagonal_only .and. i /= j) cycle
                                    pc%blocks(i, j, node) = pc%blocks(i, j, node) &
                                        + dot_product(jet(:, i), work(:, j))
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
        do node = 1, pc%nodes
            do j = 1, 3
                do i = j + 1, 3
                    pc%blocks(i, j, node) = pc%blocks(j, i, node)
                end do
                pc%blocks(j, j, node) = pc%blocks(j, j, node) + shift
                fixed(j) = map%fixed((j - 1)*pc%nodes + node)
            end do
            do j = 1, 3
                if (.not. fixed(j)) cycle
                pc%blocks(j, :, node) = 0; pc%blocks(:, j, node) = 0
                pc%blocks(j, j, node) = 1
            end do
            if (.not. all(ieee_is_finite(pc%blocks(:, :, node)))) then
                status = 2
                return
            end if
            if (pc%diagonal_only) then
                do j = 1, 3
                    if (pc%blocks(j, j, node) > 0) cycle
                    status = 4
                    return
                end do
            else
                pc%factor(:, :, node) = pc%blocks(:, :, node)
                call dpotrf('L', 3, pc%factor(:, :, node), 3, info)
                if (info /= 0) then
                    status = 4
                    return
                end if
            end if
        end do
        pc%ready = .true.
        status = 0
    end subroutine

    subroutine apply_gll_block_jacobi(pc, x, y, status)
        type(gll_block_jacobi_t), intent(in) :: pc
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer, intent(out) :: status
        real(dp) :: local(3, 1)
        integer :: node, j, info
        status = 1
        y = 0
        if (.not. pc%ready) return
        if (size(x) /= 3*pc%nodes .or. size(y) /= size(x)) return
        status = 2
        if (.not. all(ieee_is_finite(x))) return
        do node = 1, pc%nodes
            do j = 1, 3
                local(j, 1) = x((j - 1)*pc%nodes + node)
            end do
            if (pc%diagonal_only) then
                do j = 1, 3
                    local(j, 1) = local(j, 1)/pc%blocks(j, j, node)
                end do
            else
                call dpotrs('L', 3, 1, pc%factor(:, :, node), 3, local, 3, info)
                if (info /= 0) then
                    status = 4
                    y = 0
                    return
                end if
            end if
            do j = 1, 3
                y((j - 1)*pc%nodes + node) = local(j, 1)
            end do
        end do
        if (.not. all(ieee_is_finite(y))) then
            y = 0
            return
        end if
        status = 0
    end subroutine

    subroutine channel(map, c, u, v, t1, t2, transposed)
        type(gll_map_t), intent(in) :: map
        integer, intent(in) :: c
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: v(:)
        real(dp), intent(inout) :: t1(:), t2(:)
        logical, intent(in) :: transposed
        real(dp) :: b1(map%reference%no(1), map%reference%ni(1)), &
            b2(map%reference%no(1), map%reference%ni(1)), identity(1, 1)
        b1 = map%reference%line(1)%b(:, :, 1)
        b2 = map%reference%line(2)%b(:, :, 1); identity = 1
        if (derivative(c) == 1) b1 = map%reference%line(1)%d(:, :, 1)
        if (derivative(c) == 2) b2 = map%reference%line(2)%d(:, :, 1)
        if (transposed) then
            call tensor_contract(map%reference%no, map%reference%ni, &
                transpose(b1), transpose(b2), identity, u, v, t1, t2)
        else
            call tensor_contract(map%reference%ni, map%reference%no, &
                b1, b2, identity, u, v, t1, t2)
        end if
    end subroutine
end module
