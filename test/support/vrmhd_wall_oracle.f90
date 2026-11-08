module vrmhd_wall_oracle
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_gs_symmetry, only: GS_ROTATION
    use milca_vrmhd_section, only: vrmhd_section_t
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes
    implicit none
    private
    public :: current_wall_dofs
    abstract interface
        subroutine current_datum(x, data, current)
            import :: dp
            real(dp), intent(in) :: x(2), data(:)
            real(dp), intent(out) :: current(3)
        end subroutine
    end interface
contains
    subroutine current_wall_dofs(op, datum, data, cochain)
        type(vrmhd_section_t), intent(in) :: op
        procedure(current_datum) :: datum
        real(dp), intent(in) :: data(:)
        real(dp), intent(out) :: cochain(op%fs%sp(1)%ndof)
        real(dp) :: nodes(12), weights(12), x(2), current(3), horizontal(2)
        real(dp) :: start(2), tangent(2), a, pitch
        real(dp), allocatable :: xy(:, :)
        integer :: t, i, j, q, m, edge, vertex, idx
        ! Edge moments against shifted Legendre P_(m-1) on the straight edge
        ! and orbit values at the Lagrange nodes are the FortFEM dofs.
        if (op%fs%g /= 1) error stop 'straight trace oracle'
        call gauss_legendre_ab(size(nodes), 0.0_dp, 1.0_dp, nodes, weights)
        cochain = 0
        pitch = op%symmetry%k
        call meridian_lagrange_nodes(op%fs, xy)
        do vertex = 1, op%fs%n0
            if (.not. op%fs%sp(0)%bnd(vertex)) cycle
            x = xy(:, vertex)
            call datum(x, data, current)
            cochain(op%fs%nn + vertex) = orbit_component(x, current)
        end do
        do t = 1, op%fs%mm%topology%n_triangles
            do i = 1, 3
                if (.not. op%fs%sp(1)%bnd(op%fs%sp(1)%dof((i - 1)*op%fs%ke + 1, t))) &
                    cycle
                j = mod(i, 3) + 1
                start = op%fs%mm%topology%vertices(:, op%fs%mm%topology%triangles(i, t))
                tangent = op%fs%mm%topology%vertices(:, &
                    op%fs%mm%topology%triangles(j, t)) - start
                do m = 1, op%fs%ke
                    idx = (i - 1)*op%fs%ke + m
                    edge = op%fs%sp(1)%dof(idx, t)
                    do q = 1, size(nodes)
                        x = start + nodes(q)*tangent
                        call datum(x, data, current)
                        a = orbit_component(x, current)
                        if (op%symmetry%kind == GS_ROTATION) then
                            horizontal = current([1, 3])
                        else
                            horizontal = current(:2) &
                                + pitch*a*[x(2), -x(1)]/(1 + pitch**2*sum(x*x))
                        end if
                        cochain(edge) = cochain(edge) + weights(q) &
                            *legendre(m - 1, 2*nodes(q) - 1) &
                            *dot_product(horizontal, tangent)*op%fs%sp(1)%sgn(idx, t)
                    end do
                end do
            end do
        end do
    contains
        pure real(dp) function legendre(degree, s) result(value)
            integer, intent(in) :: degree
            real(dp), intent(in) :: s
            real(dp) :: previous, older
            integer :: l
            value = 1
            if (degree == 0) return
            older = 1; value = s
            do l = 1, degree - 1
                previous = value
                value = ((2*l + 1)*s*previous - l*older)/(l + 1)
                older = previous
            end do
        end function

        function orbit_component(x, current) result(a)
            real(dp), intent(in) :: x(2), current(3)
            real(dp) :: a
            if (op%symmetry%kind == GS_ROTATION) then
                a = x(1)*current(2)
            else
                a = dot_product([-pitch*x(2), pitch*x(1), 1.0_dp], current)
            end if
        end function
    end subroutine
end module vrmhd_wall_oracle
