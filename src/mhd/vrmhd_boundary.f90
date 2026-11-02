module kin6d_vrmhd_boundary
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_gs_symmetry, only: GS_ROTATION
    use kin6d_vrmhd_section, only: vrmhd_section_t, vrmhd_section_trace
    use kin6d_vrmhd_section_leaves, only: vrmhd_ampere_surface, &
        vrmhd_grad_rotation, vrmhd_grad_screw
    use fortfem_meridian_feec_n, only: meridian_point_basis_t, &
        meridian_point_basis_init
    use fortfem_meridian_curved_mesh, only: meridian_geometry_point
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private
    abstract interface
        subroutine vrmhd_trace_i(x, b)
            import :: dp
            real(dp), intent(in) :: x(2)
            real(dp), intent(out) :: b(3)
        end subroutine
    end interface
    public :: vrmhd_trace_i, vrmhd_ampere_boundary
contains
    subroutine vrmhd_ampere_boundary(op, magnetic_trace, residual, status)
        type(vrmhd_section_t), intent(in) :: op
        procedure(vrmhd_trace_i) :: magnetic_trace
        real(dp), intent(out) :: residual(op%ndof)
        integer, intent(out) :: status
        type(meridian_point_basis_t) :: pb
        integer, allocatable :: adjacent(:)
        real(dp), parameter :: vertex(2, 3) = reshape([0.0_dp, 0.0_dp, &
            1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        real(dp) :: nodes(op%fs%k + 5), weights(op%fs%k + 5)
        real(dp) :: trace(3, op%fs%sp(1)%nloc), ref(2), direction(2), tangent(2)
        real(dp) :: x(2), jac(2, 2), normal(3), b(3), value, weight
        integer :: edge, e, t, a, c, i, q, row, n0, n1
        residual = 0
        call meridian_point_basis_init(op%fs, pb, status)
        if (status /= 0) return
        call gauss_legendre_ab(size(nodes), 0.0_dp, 1.0_dp, nodes, weights)
        n0 = op%fs%n0; n1 = op%fs%sp(1)%ndof
        do edge = 1, op%fs%mm%topology%n_boundary_edges
            e = op%fs%mm%topology%boundary_edges(edge)
            call op%fs%mm%topology%get_edge_triangles(e, adjacent)
            status = 1
            if (size(adjacent) /= 1) return
            t = adjacent(1)
            do a = 1, 3
                c = mod(a, 3) + 1
                if (.not. any(op%fs%mm%topology%edges(:, e) == &
                    op%fs%mm%topology%triangles(a, t))) cycle
                if (.not. any(op%fs%mm%topology%edges(:, e) == &
                    op%fs%mm%topology%triangles(c, t))) cycle
                exit
            end do
            if (a > 3) return
            direction = vertex(:, c) - vertex(:, a)
            do q = 1, size(nodes)
                ref = vertex(:, a) + nodes(q)*direction
                call meridian_geometry_point(op%fs%mm, t, ref(1), ref(2), x, jac)
                tangent = matmul(jac, direction)
                if (op%symmetry%kind == GS_ROTATION) then
                    call vrmhd_grad_rotation(x(1), x(2), op%symmetry%k, &
                        tangent(2), -tangent(1), normal(1), normal(2), normal(3))
                    weight = weights(q)*x(1)
                else
                    call vrmhd_grad_screw(x(1), x(2), op%symmetry%k, &
                        tangent(2), -tangent(1), normal(1), normal(2), normal(3))
                    weight = weights(q)
                end if
                call vrmhd_section_trace(op, pb, t, ref, x, trace, status)
                if (status /= 0) return
                call magnetic_trace(x, b)
                do i = 1, op%fs%sp(1)%nloc
                    call vrmhd_ampere_surface(b(1), b(2), b(3), trace(1, i), &
                        trace(2, i), trace(3, i), normal(1), normal(2), &
                        normal(3), value)
                    row = 6*n0 + n1 + op%fs%sp(1)%dof(i, t)
                    residual(row) = residual(row) + weight*value
                end do
            end do
        end do
        status = 0
    end subroutine
end module kin6d_vrmhd_boundary
