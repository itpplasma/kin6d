module kin6d_nonnested_exterior
    !! NN2 on a coil-free toroidal interior with a full 3D exterior response.
    !! B = B_coil - grad u outside Gamma_c. The response has zero periods,
    !! so this entry point admits net toroidal plasma current I_phi = 0.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_parametric_laplace_bem_3d, only: bem_surface_t, bem_operators_t, &
                                   bem_exterior_dtn_t, build_exterior_dtn, assignment(=)
    use fortfem_parametric_bem_high_order, only: bem_pair_options_t, &
                                           bem_pair_report_t, assemble_parametric_bem_p2
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_generated_bem_frame, only: bem_frame
    use fortfem_meridian_feec_n, only: meridian_qpoint, meridian_point_basis
    use fortfem_meridian_fourier_curl, only: meridian_fourier_curl
    use fortnum_linalg, only: dense_solve
    use fortsparse, only: csc_t, csc_matvec, csc_transpose, csc_from_triplet, &
                          fortsparse_status_t, status_ok
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_triangle_basis, only: triangle_basis_jet
    use kin6d_mhd_mesh, only: flux_mesh_ring_edge, flux_mesh_geometry
    use kin6d_equilibrium_nonnested, only: free_b_field_t
    use kin6d_nonnested_coupled, only: nc_problem_t, nc_layout_t, nc_layout, nc_trig, &
        nc_validate
    use kin6d_coils_3d, only: coil_set_t, coil_set_b_cart
    implicit none
    private
    public :: nn_exterior_t, nn_attach_exterior, nn_exterior_trace
    public :: nn_attach_exterior_operators
    real(dp), parameter :: pi = acos(-1.0_dp)
    type :: nn_exterior_t
        type(bem_exterior_dtn_t) :: dtn
        real(dp), allocatable :: normal_map(:, :), coil_normal(:), potential_map(:, :)
        real(dp), allocatable :: flux_pairing(:, :)
        integer, allocatable :: magnetic_dofs(:)
        integer :: panels = 0, trace_dofs = 0, flux_dofs = 0
        integer(int64) :: operator_bytes = 0
        real(dp) :: assembly_seconds = 0, dtn_seconds = 0, coupling_seconds = 0
    end type
    type, extends(bem_surface_t) :: coupling_surface_t
        real(dp), allocatable :: edge(:, :, :), parameters(:, :, :)
        integer, allocatable :: panel_edge(:), triangles(:)
        integer :: order = 0
    contains
        procedure :: evaluate => surface_evaluate
    end type
contains
    pure subroutine edge_basis(order, t, w, dw)
        integer, intent(in) :: order
        real(dp), intent(in) :: t
        real(dp), intent(out) :: w(0:order), dw(0:order)
        real(dp) :: values(10), gradient(2, 10)
        integer :: index(4), nloc
        nloc = 6; index = [1, 4, 2, 0]
        if (order == 3) then
            nloc = 10; index = [1, 4, 5, 2]
        end if
        call triangle_basis_jet([t, 0.0_dp], values(:nloc), gradient(:, :nloc))
        w = values(index(:order + 1)); dw = gradient(1, index(:order + 1))
    end subroutine

    pure subroutine surface_evaluate(self, panel, uv, x, a, b, status)
        class(coupling_surface_t), intent(in) :: self
        integer, intent(in) :: panel
        real(dp), intent(in) :: uv(2)
        real(dp), intent(out) :: x(3), a(3), b(3)
        integer, intent(out) :: status
        type(cylindrical_chart_t) :: chart
        real(dp) :: p(2), deriv(2, 2), q(3), qd(3, 2), w(0:3), dw(0:3)
        real(dp) :: d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        integer :: e, k
        e = self%panel_edge(panel)
        deriv(:, 1) = self%parameters(:, 2, panel) - self%parameters(:, 1, panel)
        deriv(:, 2) = self%parameters(:, 3, panel) - self%parameters(:, 1, panel)
        p = self%parameters(:, 1, panel) + matmul(deriv, uv)
        call edge_basis(self%order, p(1), w(:self%order), dw(:self%order))
        q(:2) = matmul(self%edge(:, :, e), w(:self%order)); q(3) = p(2)
        do k = 1, 2
            qd(:2, k) = matmul(self%edge(:, :, e), dw(:self%order))*deriv(1, k)
            qd(3, k) = deriv(2, k)
        end do
        call chart%jet(q, x, d1, d2, d3)
        a = matmul(d1, qd(:, 1)); b = matmul(d1, qd(:, 2)); status = 0
    end subroutine

    subroutine make_surface(field, nphi, surface, faces, status)
        type(free_b_field_t), intent(in) :: field
        integer, intent(in) :: nphi
        type(coupling_surface_t), intent(out) :: surface
        integer, allocatable, intent(out) :: faces(:, :)
        integer, intent(out) :: status
        real(dp) :: p0, p1, th0, th1, x(2), jac(2, 2)
        integer :: ne, nf, e, j, p, i, t, after, next, ids(4), order
        status = 1
        order = field%mesh%mm%geo_order
        if (order < 2 .or. order > 3 .or. nphi < 4) return
        if (field%mesh%tlast(1) == 0) return ! no internal conductor boundary
        ne = field%mesh%nth(field%mesh%nring); nf = 2*ne*nphi
        surface%order = order
        allocate (surface%edge(2, 0:order, ne), surface%parameters(2, 3, nf), &
                  surface%panel_edge(nf), surface%triangles(ne), faces(3, nf))
        do e = 1, ne
            call flux_mesh_ring_edge(field%mesh, field%mesh%nring, e, t, th0, th1)
            surface%triangles(e) = t
            do i = 0, order
                call flux_mesh_geometry(field%mesh, t, real(i, dp)/order, &
                    0.0_dp, x, jac)
                surface%edge(:, i, e) = x
            end do
        end do
        p = 0
        do j = 1, nphi
            after = mod(j, nphi) + 1
            p0 = 2*pi*(j - 1)/nphi; p1 = 2*pi*j/nphi
            do e = 1, ne
                next = mod(e, ne) + 1
                ids = [e + (j - 1)*ne, e + (after - 1)*ne, &
                       next + (j - 1)*ne, next + (after - 1)*ne]
                p = p + 1; faces(:, p) = ids([1, 2, 3]); surface%panel_edge(p) = e
                surface%parameters(:, :, p) = &
                    reshape([0.0_dp, p0, 0.0_dp, p1, 1.0_dp, p0], [2, 3])
                p = p + 1; faces(:, p) = ids([2, 4, 3]); surface%panel_edge(p) = e
                surface%parameters(:, :, p) = &
                    reshape([0.0_dp, p1, 1.0_dp, p1, 1.0_dp, p0], [2, 3])
            end do
        end do
        status = 0
    end subroutine

    subroutine nn_attach_exterior(field, problem, coils, currents, nphi, degree, &
                                  options, exterior, report, status)
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(inout) :: problem
        type(coil_set_t), intent(in) :: coils
        real(dp), intent(in) :: currents(:)
        integer, intent(in) :: nphi, degree
        type(bem_pair_options_t), intent(in) :: options
        type(nn_exterior_t), intent(out) :: exterior
        type(bem_pair_report_t), intent(out) :: report
        integer, intent(out) :: status
        type(coupling_surface_t) :: surface
        type(bem_operators_t) :: op
        type(csc_t) :: curl, ct
        type(fortsparse_status_t) :: sparse_status
        type(nc_layout_t) :: lay
        real(dp), allocatable :: moment(:, :), normal(:, :), steklov(:, :), load(:)
        real(dp), allocatable :: coil_normal(:)
        integer, allocatable :: faces(:, :), td(:, :), fd(:, :)
        integer :: j, p, n
        integer(int64) :: tick0, tick1, rate
        status = 1; report = bem_pair_report_t()
        if (.not. allocated(field%b)) return
        if (allocated(problem%coil_drive)) return
        if (allocated(problem%floating)) return
        if (problem%closure%drive /= 'I_phi') return
        if (abs(problem%current) > 0 .or. nphi <= 2*ubound(field%b, 2)) return
        if (size(currents) /= coils%ncoil .or. degree < 2) return
        if (.not. all(ieee_is_finite(currents))) return
        call make_surface(field, nphi, surface, faces, status)
        if (status /= 0) return
        call system_clock(tick0, rate)
        call assemble_parametric_bem_p2(surface, faces, options, op, td, fd, &
                                        report, status)
        call system_clock(tick1)
        exterior%assembly_seconds = real(tick1 - tick0, dp)/rate
        if (status /= 0) return
        call system_clock(tick0)
        call build_exterior_dtn(op, exterior%dtn, status)
        if (status /= 0) return
        call exterior%dtn%matrix(steklov, status)
        if (status /= 0) return
        call system_clock(tick1)
        exterior%dtn_seconds = real(tick1 - tick0, dp)/rate
        exterior%flux_pairing = transpose(op%mass)
        exterior%panels = size(faces, 2); exterior%trace_dofs = size(op%w, 1)
        exterior%flux_dofs = size(op%v, 1)
        exterior%operator_bytes = 8_int64*(size(op%v, kind=int64) &
                               + size(op%k, kind=int64) + size(op%adjoint, kind=int64) &
                                   + size(op%w, kind=int64) + size(op%mass, kind=int64))
        call system_clock(tick0)
        lay = nc_layout(field)
        call trace_moments(field, surface, td, coils, currents, degree, &
                           moment, coil_normal, status)
        if (status /= 0) return
        call meridian_fourier_curl(field%fs, lay%modes, curl, status)
        if (status /= 0) return
        n = size(steklov, 1)
        allocate (normal(n, lay%na)); normal = 0
        !$omp parallel do default(shared) private(j, p) schedule(static)
        do j = 1, lay%na
            do p = curl%col_ptr(j), curl%col_ptr(j + 1) - 1
                normal(:, j) = normal(:, j) + curl%val(p)*moment(:, curl%row_idx(p))
            end do
        end do
        !$omp end parallel do
        deallocate (moment)
        call volume_coil_load(field, coils, currents, nphi, load, status)
        if (status /= 0) return
        call csc_transpose(curl, ct, sparse_status)
        status = 1
        if (.not. status_ok(sparse_status)) return
        load = csc_matvec(ct, load)
        call nn_attach_exterior_operators(field, problem, steklov, normal, &
            coil_normal, load, exterior, status)
        if (status /= 0) return
        call system_clock(tick1)
        exterior%coupling_seconds = real(tick1 - tick0, dp)/rate
        status = 0
    end subroutine

    subroutine nn_attach_exterior_operators(field, problem, steklov, normal, &
            coil_normal, volume_drive, exterior, status)
        ! Existing PR47 Schur composition, factored from nn_attach_exterior.
        ! S is the positive exterior energy map, normal=T in the same dual
        ! scalar trace basis and physical full-surface measure. volume_drive
        ! is C^T<V2,B_coil-B_h> in NN2's toroidal-mean volume measure.
        ! This algebraic entry does not construct/qualify a BEM dtn object.
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(inout) :: problem
        real(dp), intent(in) :: steklov(:, :), normal(:, :), coil_normal(:)
        real(dp), intent(in) :: volume_drive(:)
        type(nn_exterior_t), intent(inout) :: exterior
        integer, intent(out) :: status
        type(nc_layout_t) :: lay
        type(fortsparse_status_t) :: sparse_status
        real(dp), allocatable :: rhs(:, :), solution(:, :), h(:, :), drive(:)
        integer, allocatable :: ids(:), rows(:), cols(:)
        logical, allocatable :: mask(:)
        character(160) :: message
        integer :: n, count, i, j, k
        call nc_validate(problem, status, message)
        if (status /= 0) return
        status = 1
        if (.not. allocated(field%b)) return
        if (field%fs%cyl_r0 > 0) return
        if (.not. ieee_is_finite(field%mu0)) return
        if (field%mu0 <= 0) return
        if (allocated(problem%coil_drive)) return
        if (allocated(problem%floating)) return
        if (problem%closure%drive /= 'I_phi') return
        if (abs(problem%current) > 0) return
        lay = nc_layout(field); n = size(steklov, 1)
        if (n < 1 .or. size(steklov, 2) /= n) return
        if (any(shape(normal) /= [n, lay%na])) return
        if (size(coil_normal) /= n .or. size(volume_drive) /= lay%na) return
        if (.not. all(ieee_is_finite(steklov))) return
        if (.not. all(ieee_is_finite(normal))) return
        if (.not. all(ieee_is_finite(coil_normal))) return
        if (.not. all(ieee_is_finite(volume_drive))) return
        if (maxval(abs(steklov - transpose(steklov))) &
            > 1.0e-11_dp*max(1.0_dp, maxval(abs(steklov)))) return
        ! curl A has a normal trace determined solely by tangential A.
        allocate (mask(lay%na))
        do j = 0, lay%nb - 1
            mask(j*lay%na1 + 1:(j + 1)*lay%na1) = field%fs%sp(1)%bnd
        end do
        do j = 1, lay%na
            if (mask(j)) cycle
            if (maxval(abs(normal(:, j))) > 1.e-10_dp) return
        end do
        ids = pack([(j, j=1, lay%na)], mask); count = size(ids)
        if (count < 1) return
        allocate (rhs(n, count + 1), solution(n, count + 1))
        rhs(:, :count) = normal(:, ids); rhs(:, count + 1) = coil_normal
        call dense_solve(steklov, rhs, solution, status)
        if (status /= 0) return
        exterior%normal_map = normal(:, ids); exterior%magnetic_dofs = ids
        exterior%potential_map = solution; exterior%coil_normal = coil_normal
        exterior%trace_dofs = n
        ! Exactly one full-turn/phi-mean conversion; not a fitted multiplier.
        h = matmul(transpose(exterior%normal_map), solution(:, :count))/(2*pi)
        drive = matmul(transpose(exterior%normal_map), solution(:, count + 1))/(2*pi)
        allocate (rows(count*count), cols(count*count)); k = 0
        do j = 1, count
            do i = 1, count
                k = k + 1; rows(k) = ids(i); cols(k) = ids(j)
            end do
        end do
        call csc_from_triplet(lay%na, lay%na, rows, cols, reshape(h, [count*count]), &
            problem%exterior, sparse_status)
        status = 1
        if (.not. status_ok(sparse_status)) return
        problem%coil_drive = volume_drive
        problem%coil_drive(ids) = problem%coil_drive(ids) + drive
        status = 0
    end subroutine

    subroutine frame(rz, phi, x, rotation)
        real(dp), intent(in) :: rz(2), phi
        real(dp), intent(out) :: x(3), rotation(3, 3)
        type(cylindrical_chart_t) :: chart
        real(dp) :: d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        call chart%jet([rz, phi], x, d1, d2, d3)
        rotation(:, 1) = d1(:, 1)
        rotation(:, 2) = d1(:, 3)/rz(1)
        rotation(:, 3) = d1(:, 2)
    end subroutine

    subroutine real_modes(modes, angle, tr)
        integer, intent(in) :: modes
        real(dp), intent(in) :: angle
        real(dp), intent(out) :: tr(0:2*modes)
        real(dp) :: derivative(2*modes + 1)
        call nc_trig(modes, angle, tr, derivative)
    end subroutine

    subroutine trace_moments(field, surface, td, coils, currents, degree, &
                             moment, coil_normal, status)
        type(free_b_field_t), intent(in) :: field
        type(coupling_surface_t), intent(in) :: surface
        integer, intent(in) :: td(:, :), degree
        type(coil_set_t), intent(in) :: coils
        real(dp), intent(in) :: currents(:)
        real(dp), allocatable, intent(out) :: moment(:, :), coil_normal(:)
        integer, intent(out) :: status
        real(dp), allocatable :: xi(:), eta(:), weights(:)
        real(dp) :: local(6, field%fs%sp(2)%nloc, 0:2*ubound(field%b, 2)), drive(6)
        integer :: p, a, j, bk, row, col, t, ierr, failure, nb2, nl
        call triangle_duffy_quadrature(degree, xi, eta, weights, status)
        if (status /= 0) return
        nb2 = field%fs%sp(2)%ndof; nl = field%fs%sp(2)%nloc
        allocate (moment(maxval(td), nb2*(2*ubound(field%b, 2) + 1)), &
                  coil_normal(maxval(td)))
        moment = 0; coil_normal = 0; failure = 0
        !$omp parallel do default(shared) private(p, a, j, bk, row, col, t, ierr, &
        !$omp& local, drive) schedule(static) reduction(max:failure)
        do p = 1, size(td, 2)
            call panel_moments(p, local, drive, ierr)
            failure = max(failure, ierr)
            if (ierr /= 0) cycle
            t = surface%triangles(surface%panel_edge(p))
            do a = 1, 6
                row = td(a, p)
                !$omp atomic update
                coil_normal(row) = coil_normal(row) + drive(a)
                do bk = 0, 2*ubound(field%b, 2)
                    do j = 1, nl
                        col = bk*nb2 + field%fs%sp(2)%dof(j, t)
                        !$omp atomic update
                        moment(row, col) = moment(row, col) + local(a, j, bk)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        status = failure
    contains
        subroutine panel_moments(p, local, drive, status)
            integer, intent(in) :: p
            real(dp), intent(out) :: local(:, :, 0:), drive(:)
            integer, intent(out) :: status
            real(dp) :: x(3), a(3), b(3), normal(3), ncy(3), jacobian, bc(3)
            real(dp) :: par(2), rz(2), basis(3, nl), v(6), grad(2, 6)
            real(dp) :: rot(3, 3), unused(3), tr(0:2*ubound(field%b, 2)), w
            integer :: q, t, i, j, bk
            t = surface%triangles(surface%panel_edge(p)); local = 0; drive = 0
            status = 0
            do q = 1, size(weights)
                call surface%evaluate(p, [xi(q), eta(q)], x, a, b, status)
                if (status /= 0) return
                call bem_frame(a(1), a(2), a(3), b(1), b(2), b(3), &
                               jacobian, normal(1), normal(2), normal(3))
                par = surface%parameters(:, 1, p) &
                   + xi(q)*(surface%parameters(:, 2, p) - surface%parameters(:, 1, p)) &
                    + eta(q)*(surface%parameters(:, 3, p) - surface%parameters(:, 1, p))
                call meridian_point_basis(field%fs, field%pb, 2, t, par(1), &
                                          0.0_dp, rz, basis, status)
                if (status /= 0) return
                call frame(rz, par(2), unused, rot)
                ncy = matmul(transpose(rot), normal)
                call triangle_basis_jet([xi(q), eta(q)], v, grad)
                call real_modes(ubound(field%b, 2), par(2), tr)
                call coil_set_b_cart(coils, currents, x, bc)
                w = weights(q)*jacobian
                drive = drive + w*v*dot_product(normal, bc)
                do bk = 0, ubound(tr, 1)
                    do j = 1, nl
                        do i = 1, 6
                            local(i, j, bk) = local(i, j, bk) + w*v(i)*tr(bk) &
                                 *dot_product(ncy, basis(:, j))*field%fs%sp(2)%sgn(j, t)
                        end do
                    end do
                end do
            end do
        end subroutine
    end subroutine

    subroutine volume_coil_load(field, coils, currents, nphi, load, status)
        type(free_b_field_t), intent(in) :: field
        type(coil_set_t), intent(in) :: coils
        real(dp), intent(in) :: currents(:)
        integer, intent(in) :: nphi
        integer, intent(out) :: status
        real(dp), allocatable, intent(out) :: load(:)
        real(dp), allocatable :: local(:, :, :)
        real(dp) :: x(3), rz(2), bc(3), b(3), rotation(3, 3), w, phi
        real(dp) :: basis(3, field%fs%sp(2)%nloc), tr(0:2*ubound(field%b, 2))
        integer :: t, q, j, bk, i, nl, nb2, nt, failure, ierr
        nl = field%fs%sp(2)%nloc; nb2 = field%fs%sp(2)%ndof
        nt = field%fs%mm%topology%n_triangles
        allocate (local(nl, 0:2*ubound(field%b, 2), nt))
        local = 0; failure = 0
        !$omp parallel do default(shared) private(t, q, j, bk, i, x, rz, bc, b, &
        !$omp& rotation, w, phi, basis, tr, ierr) schedule(static) &
        !$omp& reduction(max:failure)
        do t = 1, nt
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 2, t, q, rz, w, basis, ierr)
                failure = max(failure, ierr)
                if (ierr /= 0) cycle
                do j = 1, nphi
                    phi = 2*pi*(j - 1)/nphi
                    call frame(rz, phi, x, rotation)
                    call coil_set_b_cart(coils, currents, x, bc)
                    b = matmul(transpose(rotation), bc)
                    b(2) = b(2) - field%f_vac/rz(1)
                    call real_modes(ubound(field%b, 2), phi, tr)
                    do bk = 0, ubound(tr, 1)
                        do i = 1, nl
                            local(i, bk, t) = local(i, bk, t) + w/nphi*tr(bk) &
                                   *dot_product(basis(:, i), b)*field%fs%sp(2)%sgn(i, t)
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        status = failure
        if (status /= 0) return
        allocate (load(nb2*(2*ubound(field%b, 2) + 1))); load = 0
        do t = 1, nt
            do bk = 0, 2*ubound(field%b, 2)
                load(bk*nb2 + field%fs%sp(2)%dof(:, t)) = &
                    load(bk*nb2 + field%fs%sp(2)%dof(:, t)) + local(:, bk, t)
            end do
        end do
    end subroutine

    subroutine nn_exterior_trace(exterior, magnetic, potential, flux, load, &
                                 interface_defect, status, physical_defect)
        !! The algebraic weak normal defect uses S u and is distinct from a
        !! continuous field/trace error or a virtual-casing comparison.
        type(nn_exterior_t), intent(in) :: exterior
        real(dp), intent(in) :: magnetic(:)
        real(dp), allocatable, intent(out) :: potential(:), flux(:), load(:)
        real(dp), intent(out), optional :: physical_defect
        real(dp), intent(out) :: interface_defect
        integer, intent(out) :: status
        real(dp), allocatable :: normal(:)
        status = 1; interface_defect = huge(1.0_dp)
        if (.not. allocated(exterior%potential_map)) return
        if (size(magnetic) < maxval(exterior%magnetic_dofs)) return
        if (.not. all(ieee_is_finite(magnetic))) return
        potential = matmul(exterior%potential_map(:, :size(exterior%magnetic_dofs)), &
                           magnetic(exterior%magnetic_dofs)) &
                    - exterior%potential_map(:, size(exterior%magnetic_dofs) + 1)
        allocate (flux(exterior%flux_dofs), load(exterior%trace_dofs))
        call exterior%dtn%apply(potential, flux, load, status)
        if (status /= 0) return
        normal = matmul(exterior%normal_map, magnetic(exterior%magnetic_dofs)) &
                 - exterior%coil_normal
        interface_defect = norm2(load - normal)/max(norm2(normal), tiny(1.0_dp))
        if (present(physical_defect)) physical_defect = &
            norm2(normal + matmul(exterior%flux_pairing, flux)) &
            /max(norm2(exterior%coil_normal), tiny(1.0_dp))
    end subroutine
end module
