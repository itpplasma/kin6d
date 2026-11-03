module vrmhd_periodic_slab
    !! Doubly periodic translation slab [0, 2 pi]^2 as a coefficient
    !! identification P on the rectangular section mesh: u = lift + P y.
    !! Edge signs follow the Bloch transform of vrmhd_spectrum_support.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_gs_symmetry, only: gs_symmetry_t, GS_TRANSLATION, GS_CHART_CARTESIAN
    use kin6d_vrmhd_section, only: vrmhd_section_t, vrmhd_section_build
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes, meridian_load
    use fortsparse, only: csc_t, csc_from_triplet, csc_transpose, csc_matmul, &
        csc_matvec, sparse_solve_once, fortsparse_status_t, status_ok
    implicit none
    private
    type, public :: periodic_slab_t
        type(vrmhd_section_t) :: op
        integer :: n = 0, ns = 0, na = 0, nred = 0
        real(dp) :: length = 0
        !> Full-to-periodic identification and its vector-potential block.
        type(csc_t) :: reduction, potential
        integer, allocatable :: scalar_class(:)
        real(dp), allocatable :: scalar_nodes(:, :)
    end type
    abstract interface
        subroutine vector_field_t(x, a)
            import :: dp
            real(dp), intent(in) :: x(2)
            complex(dp), intent(out) :: a(3)
        end subroutine
    end interface
    public :: periodic_slab_build, periodic_scalar, periodic_potential
    public :: vector_field_t
contains
    subroutine periodic_slab_build(n, degree, slab, status, magnetic_unit)
        integer, intent(in) :: n, degree
        type(periodic_slab_t), intent(out) :: slab
        integer, intent(out) :: status
        !> Physical A per reduced A coefficient; J uses its reciprocal.
        real(dp), intent(in), optional :: magnetic_unit
        type(mesh_2d_t) :: mesh
        type(meridian_mesh_t) :: mm
        type(fortsparse_status_t) :: st
        real(dp), allocatable :: xy(:, :), vals(:), avals(:)
        integer, allocatable :: rows(:), cols(:), arows(:), acols(:)
        integer, allocatable :: edge_class(:), key_class(:)
        real(dp), allocatable :: edge_sign(:)
        real(dp) :: h, start(2), tangent(2), midpoint(2), direction, sgn
        integer :: n0, n1, nn, ke, k0, ks, t, i, j, m, idx, g, key, nclass, pos
        integer :: field, apos, a0
        status = 1
        if (n < 1) return
        if (present(magnetic_unit)) then
            if (magnetic_unit <= 0) return
        end if
        slab%n = n
        slab%length = 2*acos(-1.0_dp)
        h = slab%length/n
        call mesh%create_rectangular(n + 1, n + 1, 0.0_dp, slab%length, 0.0_dp, &
            slab%length)
        call meridian_mesh_straight(mesh, 1, mm, status)
        if (status /= 0) return
        call vrmhd_section_build(mm, gs_symmetry_t(GS_TRANSLATION, &
            GS_CHART_CARTESIAN, 0.0_dp), degree, 1, slab%op, status)
        if (status /= 0) return
        status = 1
        associate (op => slab%op)
            n0 = op%scalar_fs%n0; n1 = op%fs%sp(1)%ndof; nn = op%fs%nn
            ke = op%fs%ke; k0 = op%fs%k0; ks = op%scalar_fs%k0
            call meridian_lagrange_nodes(op%scalar_fs, xy)
            slab%ns = (n*ks)**2
            allocate (slab%scalar_class(n0), slab%scalar_nodes(2, slab%ns))
            do i = 1, n0
                slab%scalar_class(i) = nodal_class(xy(:, i), ks)
                slab%scalar_nodes(:, slab%scalar_class(i)) = xy(:, i)
            end do
            ! Edge moments are keyed by their midpoint on the half-cell lattice.
            allocate (edge_class(n1), edge_sign(n1), key_class(4*n*n*ke))
            edge_class = 0; edge_sign = 0; key_class = 0; nclass = 0
            do t = 1, op%fs%mm%topology%n_triangles
                do i = 1, 3
                    j = mod(i, 3) + 1
                    start = op%fs%mm%topology%vertices(:, op%fs%mm%topology%triangles(i, t))
                    tangent = op%fs%mm%topology%vertices(:, &
                        op%fs%mm%topology%triangles(j, t)) - start
                    midpoint = start + tangent/2
                    direction = sign(1.0_dp, tangent(1))
                    if (abs(tangent(1)) < epsilon(h)) direction = sign(1.0_dp, tangent(2))
                    key = modulo(nint(2*midpoint(1)/h), 2*n) &
                        + 2*n*modulo(nint(2*midpoint(2)/h), 2*n)
                    do m = 1, ke
                        idx = (i - 1)*ke + m
                        g = op%fs%sp(1)%dof(idx, t)
                        sgn = real(op%fs%sp(1)%sgn(idx, t), dp)*direction**m
                        if (key_class(key*ke + m) == 0) then
                            nclass = nclass + 1
                            key_class(key*ke + m) = nclass
                        end if
                        if (edge_class(g) /= 0) then
                            if (edge_class(g) /= key_class(key*ke + m)) return
                            if (edge_sign(g) /= sgn) return
                        end if
                        edge_class(g) = key_class(key*ke + m); edge_sign(g) = sgn
                    end do
                end do
                do idx = 3*ke + 1, op%fs%nln
                    g = op%fs%sp(1)%dof(idx, t)
                    if (edge_class(g) /= 0) return
                    nclass = nclass + 1
                    edge_class(g) = nclass; edge_sign(g) = 1
                end do
            end do
            deallocate (xy)
            call meridian_lagrange_nodes(op%fs, xy)
            do i = 1, op%fs%n0
                edge_class(nn + i) = nclass + nodal_class(xy(:, i), k0)
                edge_sign(nn + i) = 1
            end do
            if (any(edge_class == 0)) return
            slab%na = nclass + (n*k0)**2
            slab%nred = 6*slab%ns + 2*slab%na
            allocate (rows(6*n0 + 2*n1), cols(6*n0 + 2*n1), vals(6*n0 + 2*n1), &
                arows(n1), acols(n1), avals(n1))
            pos = 0
            do field = 1, 6
                do i = 1, n0
                    pos = pos + 1
                    rows(pos) = (field - 1)*n0 + i
                    cols(pos) = (field - 1)*slab%ns + slab%scalar_class(i)
                    vals(pos) = 1
                end do
            end do
            apos = 0
            do field = 1, 2
                a0 = 6*n0 + (field - 1)*n1
                do g = 1, n1
                    pos = pos + 1
                    rows(pos) = a0 + g
                    cols(pos) = 6*slab%ns + (field - 1)*slab%na + edge_class(g)
                    vals(pos) = edge_sign(g)
                    if (field == 1) then
                        apos = apos + 1
                        arows(apos) = g; acols(apos) = edge_class(g)
                        avals(apos) = edge_sign(g)
                    end if
                    if (present(magnetic_unit)) then
                        if (field == 1) then
                            vals(pos) = vals(pos)*magnetic_unit
                            avals(apos) = avals(apos)*magnetic_unit
                        else
                            vals(pos) = vals(pos)/magnetic_unit
                        end if
                    end if
                end do
            end do
            call csc_from_triplet(op%ndof, slab%nred, rows, cols, vals, &
                slab%reduction, st)
            if (.not. status_ok(st)) return
            call csc_from_triplet(n1, slab%na, arows, acols, avals, slab%potential, st)
            if (.not. status_ok(st)) return
        end associate
        status = 0
    contains
        integer function nodal_class(x, degree) result(label)
            real(dp), intent(in) :: x(2)
            integer, intent(in) :: degree
            label = 1 + modulo(nint(degree*x(1)/h), n*degree) &
                + n*degree*modulo(nint(degree*x(2)/h), n*degree)
        end function
    end subroutine

    subroutine periodic_scalar(slab, f, y)
        !! Nodal interpolation of a periodic scalar into one scalar block.
        type(periodic_slab_t), intent(in) :: slab
        interface
            pure real(dp) function f(x)
                import :: dp
                real(dp), intent(in) :: x(2)
            end function
        end interface
        real(dp), intent(out) :: y(slab%ns)
        integer :: i
        do i = 1, slab%ns
            y(i) = f(slab%scalar_nodes(:, i))
        end do
    end subroutine

    subroutine periodic_potential(slab, a, y, status)
        !! L2 projection of a periodic vector potential onto the periodic A block.
        type(periodic_slab_t), intent(in) :: slab
        procedure(vector_field_t) :: a
        real(dp), intent(out) :: y(slab%na)
        integer, intent(out) :: status
        complex(dp), allocatable :: load(:)
        type(csc_t) :: transpose, half, mass
        type(fortsparse_status_t) :: st
        y = 0
        call meridian_load(slab%op%fs, 1, a, load, status)
        if (status /= 0) return
        status = 1
        call csc_transpose(slab%potential, transpose, st)
        if (.not. status_ok(st)) return
        call csc_matmul(slab%op%fs%mass(1), slab%potential, half, st)
        if (.not. status_ok(st)) return
        call csc_matmul(transpose, half, mass, st)
        if (.not. status_ok(st)) return
        call sparse_solve_once(mass, csc_matvec(transpose, real(load, dp)), y, st)
        if (status_ok(st)) status = 0
    end subroutine
end module vrmhd_periodic_slab
