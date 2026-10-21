module kin6d_equilibrium_boundary_3d
    !! R2 of free-boundary-theory.md. Full-surface Cartesian Laplace BEM;
    !! R maps the consumer trace into continuous P2 BEM nodal values.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_parametric_laplace_bem_3d, only: bem_surface_t, bem_operators_t, &
                                                 bem_exterior_dtn_t, build_exterior_dtn
    use fortfem_parametric_bem_high_order, only: bem_pair_options_t, &
                                           bem_pair_report_t, assemble_parametric_bem_p2
    use fortfem_exterior_harmonic_fields, only: assemble_exterior_harmonics, &
                                                prescribe_exterior_periods
    implicit none
    private
    public :: gamma_c_3d_t, gamma_c_3d_build, gamma_c_3d_energy, assignment(=)
    public :: gamma_c_3d_harmonics
    type :: gamma_c_3d_t
        real(dp), allocatable :: restriction(:, :), s(:, :)
        real(dp), allocatable :: cross(:, :), gram(:, :), coefficients(:)
        integer :: genus = -1
    end type
    interface assignment(=)
        module procedure assign_boundary
    end interface
contains
    subroutine assign_boundary(lhs, rhs)
        type(gamma_c_3d_t), intent(out) :: lhs
        type(gamma_c_3d_t), intent(in) :: rhs
        if (allocated(rhs%restriction)) lhs%restriction = rhs%restriction
        if (allocated(rhs%s)) lhs%s = rhs%s
        if (allocated(rhs%cross)) lhs%cross = rhs%cross
        if (allocated(rhs%gram)) lhs%gram = rhs%gram
        if (allocated(rhs%coefficients)) lhs%coefficients = rhs%coefficients
        lhs%genus = rhs%genus
    end subroutine
    subroutine gamma_c_3d_build(surface,faces,restriction,options,boundary, &
            report,status)
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: faces(:, :)
        real(dp), intent(in) :: restriction(:, :)
        type(bem_pair_options_t), intent(in) :: options
        type(gamma_c_3d_t), intent(out) :: boundary
        type(bem_pair_report_t), intent(out) :: report
        integer, intent(out) :: status
        type(bem_operators_t) :: op
        type(bem_exterior_dtn_t) :: exterior
        real(dp), allocatable :: s(:, :)
        integer, allocatable :: td(:, :), fd(:, :)
        status = 1; report = bem_pair_report_t()
        if (size(restriction, 2) < 1) return
        if (.not. all(ieee_is_finite(restriction))) return
        call assemble_parametric_bem_p2(surface,faces,options,op,td,fd, &
            report,status)
        if (status /= 0) return
        status = 1
        if (size(restriction, 1) /= size(op%w, 1)) return
        call build_exterior_dtn(op, exterior, status)
        if (status /= 0) return
        call exterior%matrix(s, status)
        if (status /= 0) return
        boundary%restriction = restriction
        boundary%s = matmul(transpose(restriction), matmul(s, restriction))
        boundary%genus = op%genus
    end subroutine
    subroutine gamma_c_3d_harmonics(boundary, surface, faces, trace_dofs, degree, &
                 vertices, offsets, periods, prescribed, independent, tolerance, status)
        !! Periods must be measured on the same oriented wires and cycles.
        type(gamma_c_3d_t), intent(inout) :: boundary
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: faces(:, :), trace_dofs(:, :), degree, offsets(:)
        integer, intent(in) :: independent(:)
        real(dp), intent(in) :: vertices(:, :), periods(:, :), prescribed(:), tolerance
        integer, intent(out) :: status
        real(dp), allocatable :: normal(:, :), gram(:, :), coefficients(:)
        real(dp) :: asymmetry
        integer :: g
        status = 1
        if (.not. allocated(boundary%restriction)) return
        g = boundary%genus
        if (g < 1 .or. size(periods, 2) /= g) return
        if (size(offsets) /= g + 1) return
        allocate (coefficients(g))
        call prescribe_exterior_periods(periods, prescribed, independent, tolerance, &
                                        coefficients, status)
        if (status /= 0) return
        call assemble_exterior_harmonics(surface, faces, trace_dofs, degree, &
                                     vertices, offsets, normal, gram, asymmetry, status)
        if (status /= 0) return
        status = 1
        if (size(normal, 1) /= size(boundary%restriction, 1)) return
        boundary%cross = -matmul(transpose(boundary%restriction), normal)
        boundary%gram = gram; boundary%coefficients = coefficients
        status = 0
    end subroutine
    subroutine gamma_c_3d_energy(boundary,trace,drive,normal_load,mu0,energy, &
            load,status)
        !! Positive scalar vacuum energy and derivative. Coil drive is the
        !! SINGLE-VALUED remainder; harmonic periods are separate data.
        type(gamma_c_3d_t), intent(in) :: boundary
        real(dp), intent(in) :: trace(:), drive(:), normal_load(:), mu0
        real(dp), intent(out) :: energy, load(:)
        integer, intent(out) :: status
        real(dp), allocatable :: u(:), su(:), cross_load(:)
        integer :: n
        status = 1; energy = 0; load = 0
        if (.not. allocated(boundary%s)) return
        n = size(boundary%s, 1)
        if (size(trace) /= n .or. size(drive) /= n .or. size(normal_load) /= n) return
        if (size(load) /= n) return
        if (.not. ieee_is_finite(mu0)) return
        if (mu0 <= 0) return
        if (.not. all(ieee_is_finite(trace))) return
        if (.not. all(ieee_is_finite(drive))) return
        if (.not. all(ieee_is_finite(normal_load))) return
        u = trace - drive; su = matmul(boundary%s, u)
        energy = (0.5_dp*dot_product(u, su) - dot_product(normal_load, trace))/mu0
        load = (su - normal_load)/mu0
        if (allocated(boundary%coefficients)) then
            if (.not. allocated(boundary%cross)) return
            if (.not. allocated(boundary%gram)) return
            cross_load = matmul(boundary%cross, boundary%coefficients)
            energy = energy + (dot_product(u, cross_load) + &
                               0.5_dp*dot_product(boundary%coefficients, &
                                      matmul(boundary%gram, boundary%coefficients)))/mu0
            load = load + cross_load/mu0
        end if
        status = 0
    end subroutine
end module
