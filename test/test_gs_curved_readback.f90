program test_gs_curved_readback
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
        write_gs_mesh, read_gs_mesh, sample_gs
    use kin6d_gs_field_readback, only: sample_gs_many
    use kin6d_gs_readback, only: gs_readback_t, read_gs_equilibrium
    implicit none
    real(dp), parameter :: r0 = 6.2_dp, a = 0.62_dp, f = 32.86_dp
    real(dp), parameter :: c = f/(12*r0**3)
    type(gs_state_t) :: state, restored, split, flat
    type(gs_readback_t) :: result, reversed
    real(dp), allocatable :: contour(:, :)
    real(dp) :: h, s, angle, radius, exact_q, exact_phi, exact_volume, flux
    real(dp) :: q_error(3), volume_error(3), phi_error(3), axis_error(3)
    real(dp) :: label_error(3), rates(4)
    real(dp) :: values(4), restored_values(4), point(2)
    real(dp) :: query_points(2, 20), batch(4, 20)
    real(dp) :: split_query(2, 1), split_batch(4, 1)
    real(dp) :: near_query(2, 1), near_batch(4, 1)
    integer :: batch_status(20)
    integer :: split_status(1), near_status(1)
    integer :: n, mesh, j, i, status, degree, nnode, ncell
    character(8) :: argument
    ! The scalar evaluator accepts a small barycentric tolerance outside a
    ! straight edge; bulk readback must preserve the same field there.
    flat%element_order = 2
    flat%mapped = .false.
    allocate (flat%points(2, 6), flat%psi(6), flat%boundary(6), &
        flat%triangles(6, 1), flat%gradient(2, 3, 1), flat%area(1), &
        flat%curved_cell(1))
    flat%points = reshape([2.0_dp, 0.0_dp, 3.0_dp, 0.0_dp, &
        2.0_dp, 1.0_dp, 2.5_dp, 0.0_dp, 2.5_dp, 0.5_dp, &
        2.0_dp, 0.5_dp], [2, 6])
    flat%psi = sum(flat%points, dim=1)
    flat%boundary = .true.
    flat%triangles(:, 1) = [1, 2, 3, 4, 5, 6]
    flat%gradient(:, 1, 1) = [-1.0_dp, -1.0_dp]
    flat%gradient(:, 2, 1) = [1.0_dp, 0.0_dp]
    flat%gradient(:, 3, 1) = [0.0_dp, 1.0_dp]
    flat%area = 0.5_dp
    flat%curved_cell = .false.
    near_query(:, 1) = [2.5_dp, -5.0e-11_dp]
    call sample_gs(flat, near_query(1, 1), near_query(2, 1), f, &
        values(1), values(2), values(3), values(4), status)
    if (status /= 0) error stop 'scalar edge tolerance'
    call sample_gs_many(flat, near_query, f, near_batch, near_status)
    if (near_status(1) /= status) error stop 'batch edge tolerance'
    if (maxval(abs(values - near_batch(:, 1))) > 1.0e-12_dp) &
        error stop 'batch edge field'
    degree = 2
    if (command_argument_count() > 0) then
        call get_command_argument(1, argument)
        read (argument, *) degree
    end if
    print '(a)', 'n axis_error volume_rel Phi_rel q_rel_max s_tor_abs_max'
    do mesh = 1, 3
        n = 48*2**(mesh - 1)
        allocate (contour(2, n))
        do j = 1, n
            call boundary(real(j - 1, dp)/n, contour(:, j))
        end do
        h = 2*acos(-1.0_dp)*a/n
        call initialize_gs(state, contour, h*h/2, status, degree, curve=boundary)
        if (status /= 0) error stop 'readback mesh'
        call solve_gs(state, 1.0_dp, -16*c, 0.0_dp, status)
        if (status /= 0) error stop 'readback solve'
        call read_gs_equilibrium(state, f, 0.0_dp, 257, 128, result, status)
        if (status /= 0) then
            print *, 'readback status', status
            error stop 'curved readback'
        end if
        exact_volume = 0; exact_phi = 0
        do i = 1, 4096
            angle = 2*acos(-1.0_dp)*(i - 0.5_dp)/4096
            radius = sqrt(r0*r0 + 2*r0*a*cos(angle))
            exact_phi = exact_phi + sin(angle)**2/radius**3
            ! Green's theorem: volume = pi * integral R^2 dZ.
            exact_volume = exact_volume + radius*r0*a*cos(angle) &
                + r0*r0*a*a*sin(angle)**2/radius
        end do
        exact_phi = exact_phi/4096*(2*acos(-1.0_dp)*f*r0*r0*a*a)
        exact_volume = exact_volume/4096*2*acos(-1.0_dp)**2
        axis_error(mesh) = norm2(result%axis - [r0, 0.0_dp])
        volume_error(mesh) = abs(result%volume/exact_volume - 1)
        phi_error(mesh) = abs(result%phi_edge/exact_phi - 1)
        q_error(mesh) = 0; label_error(mesh) = 0
        do j = 2, size(result%q)
            ! Compare at the actual returned physical psi, not a fitted axis label.
            s = 1 + result%psi(j)/(4*c*r0*r0*a*a)
            exact_q = 0; flux = 0
            do i = 1, 4096
                angle = 2*acos(-1.0_dp)*(i - 0.5_dp)/4096
                radius = sqrt(r0*r0 + 2*r0*a*sqrt(s)*cos(angle))
                exact_q = exact_q + 1/radius**3
                flux = flux + sin(angle)**2/radius**3
            end do
            exact_q = exact_q/4096*1.5_dp*r0**3
            flux = flux/4096*(2*acos(-1.0_dp)*f*r0*r0*a*a*s)
            label_error(mesh) = max(label_error(mesh), abs(result%s_tor(j) &
                - flux/exact_phi))
            if (s < 0.05_dp .or. s > 0.98_dp) cycle
            q_error(mesh) = max(q_error(mesh), abs(result%q(j)/exact_q - 1))
        end do
        print '(i4,5es14.5)', n, axis_error(mesh), volume_error(mesh), &
            phi_error(mesh), q_error(mesh), label_error(mesh)
        if (mesh == 1) then
            call write_gs_mesh(state, 'curved_readback_mesh.dat', status)
            if (status /= 0) error stop 'native export'
            call read_gs_mesh(restored, 'curved_readback_mesh.dat', status)
            if (status /= 0) error stop 'native readback'
            do i = 1, 17
                query_points(:, i) = [r0 + 0.53_dp*cos(real(i, dp)), &
                    0.51_dp*sin(real(i, dp))]
            end do
            query_points(:, 18) = state%points(:, state%triangles(1, 1))
            query_points(:, 19) = [r0 + 2*a, 0.0_dp]
            ! Inside the curved edge, beyond its affine chord.
            call boundary(0.5_dp/n, query_points(:, 20))
            query_points(1, 20) = query_points(1, 20) - 2.0e-4_dp
            call sample_gs_many(restored, query_points, f, batch, batch_status)
            if (batch_status(20) /= 0) error stop 'curved boundary bulge'
            do i = 1, 20
                point = query_points(:, i)
                call sample_gs(restored, point(1), point(2), f, values(1), values(2), &
                    values(3), values(4), status)
                if (batch_status(i) /= status) error stop 'batch location status'
                if (status == 0) then
                    if (maxval(abs(values - batch(:, i))) > 1.0e-12_dp) &
                        error stop 'batch physical field'
                end if
            end do
            ! The walking locator may exit the first component even when a
            ! point belongs to a separate affine component.
            nnode = size(restored%points, 2)
            ncell = size(restored%triangles, 2)
            split%element_order = restored%element_order
            split%mapped = restored%mapped
            allocate (split%points(2, 2*nnode), split%psi(2*nnode), &
                split%boundary(2*nnode), split%triangles(size(restored%triangles, 1), &
                2*ncell), split%gradient(2, 3, 2*ncell), split%area(2*ncell), &
                split%curved_cell(2*ncell))
            split%points(:, :nnode) = restored%points
            split%points(:, nnode + 1:) = restored%points
            split%points(1, nnode + 1:) = split%points(1, nnode + 1:) + 3*a
            split%psi = [restored%psi, restored%psi]
            split%boundary = [restored%boundary, restored%boundary]
            split%triangles(:, :ncell) = restored%triangles
            split%triangles(:, ncell + 1:) = restored%triangles + nnode
            split%gradient(:, :, :ncell) = restored%gradient
            split%gradient(:, :, ncell + 1:) = restored%gradient
            split%area = [restored%area, restored%area]
            split%curved_cell = [restored%curved_cell, restored%curved_cell]
            split_query(:, 1) = query_points(:, 1) + [3*a, 0.0_dp]
            call sample_gs_many(split, split_query, f, split_batch, split_status)
            call sample_gs(split, split_query(1, 1), split_query(2, 1), f, &
                values(1), values(2), values(3), values(4), status)
            if (status /= 0 .or. split_status(1) /= 0) error stop 'split mesh location'
            if (maxval(abs(values - split_batch(:, 1))) > 1.0e-12_dp) &
                error stop 'split mesh field'
            do i = 1, 17
                point = [r0 + 0.53_dp*cos(real(i, dp)), 0.51_dp*sin(real(i, dp))]
                call sample_gs(state, point(1), point(2), f, values(1), values(2), &
                    values(3), values(4), status)
                if (status /= 0) error stop 'native sample'
                call sample_gs(restored, point(1), point(2), f, restored_values(1), &
                    restored_values(2), restored_values(3), restored_values(4), status)
                if (status /= 0) error stop 'native reloaded sample'
                if (maxval(abs(values - restored_values)) > 1.0e-12_dp) &
                    error stop 'native export readback round trip'
            end do
            call read_gs_equilibrium(state, -f, 0.0_dp, 257, 128, reversed, status)
            if (status /= 0) error stop 'signed readback'
            if (maxval(abs(result%q + reversed%q)) > 1.0e-13_dp) &
                error stop 'signed q'
            if (abs(result%phi_edge + reversed%phi_edge) > 1.0e-13_dp) &
                error stop 'signed toroidal flux'
        end if
        deallocate (contour)
    end do
    do mesh = 2, 3
        rates = log([axis_error(mesh - 1)/axis_error(mesh), &
            volume_error(mesh - 1)/volume_error(mesh), &
            phi_error(mesh - 1)/phi_error(mesh), &
            q_error(mesh - 1)/q_error(mesh)])/log(2.0_dp)
        print '(a,4f10.4)', 'axis, volume, Phi, q rates: ', rates
        if (rates(2) < 3.7_dp .or. rates(3) < 3.7_dp) error stop 'geometry measures'
        if (q_error(mesh) >= q_error(mesh - 1)) error stop 'q refinement'
    end do
    if (q_error(3) > 3.0e-4_dp) error stop 'q accuracy'
    if (label_error(3) > 1.0e-6_dp) error stop 'toroidal flux label'
contains
    subroutine boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: theta
        theta = 2*acos(-1.0_dp)*t
        point(1) = sqrt(r0*r0 + 2*r0*a*cos(theta))
        point(2) = r0*a*sin(theta)/point(1)
    end subroutine
end program
