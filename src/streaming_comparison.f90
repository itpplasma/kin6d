module kin6d_streaming_comparison
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, hull, intersect, is_empty, &
        linear_interpolation_remainder
    use kin6d_streaming_leaves, only: streaming_rotation, streaming_rotation_interval, &
        streaming_cell_curvature_interval, streaming_add_interval, &
        streaming_density_cell, streaming_density_cell_interval, &
        streaming_density_radius_interval, streaming_enlarge_interval, &
        streaming_interpolate, streaming_interpolate_interval, &
        streaming_mode_accumulate, streaming_mode_accumulate_interval, &
        streaming_x_curvature_interval, streaming_particle_cell_midpoint, &
        streaming_particle_cell_midpoint_interval, streaming_particle_weight, &
        streaming_particle_weight_interval
    implicit none
    private
    public :: streaming_state_t, init_streaming, advance_streaming
    public :: evaluate_streaming, evaluate_streaming_interval, streaming_radius
    public :: evaluate_streaming_candidate_interval
    public :: evaluate_streaming_initial_interval, streaming_shape
    public :: get_streaming_initial
    public :: density_streaming, density_streaming_chord, density_coefficients_streaming
    public :: streaming_particles_t, init_streaming_particles, get_streaming_particles
    public :: streaming_particle_info, streaming_particles_source_matches

    ! Supplied initial data defines the model exactly, without an omitted-mode tail.
    ! At absolute t, the spectral nodal ODE is solved by a generated linear map.
    ! Each physical Fourier harmonic is C2 inside each original velocity cell;
    ! the P1 chord remainder does not assume smoothness across initial knots.
    ! Generated derivative leaves use Euclidean pairings; no continuous tangent
    ! or particle pointwise distribution certificate is claimed here.
    type :: streaming_state_t
        private
        real(dp), allocatable :: v(:), c0(:, :), s0(:, :), c(:, :), s(:, :)
        type(interval_t), allocatable :: ci(:, :), si(:, :), cell_radius(:)
        real(dp), allocatable :: dc(:), ds(:)
        type(interval_t), allocatable :: dci(:), dsi(:)
        type(interval_t) :: density_radius
        real(dp) :: t = 0.0_dp
        logical :: initialized = .false., prepared = .false.
    end type streaming_state_t

    type :: streaming_particles_t
        private
        integer :: nx = 0, nv = 0, cells = 0
        logical :: global_velocity = .false., prepared = .false.
        real(dp), allocatable :: grid(:), c0(:, :), s0(:, :)
        real(dp), allocatable :: x(:), v(:), w(:)
        type(interval_t), allocatable :: xi(:), vi(:), wi(:)
    end type streaming_particles_t

contains

    subroutine init_streaming(state, velocity, cosine0, sine0, ok)
        type(streaming_state_t), intent(out) :: state
        real(dp), intent(in) :: velocity(:), cosine0(0:, :), sine0(0:, :)
        logical, intent(out) :: ok
        integer :: j, modes, nodes

        ok = .false.
        nodes = size(velocity)
        modes = size(cosine0, 1) - 1
        if (nodes < 2 .or. modes < 0) return
        if (any(shape(cosine0) /= shape(sine0))) return
        if (size(cosine0, 2) /= nodes) return
        if (.not. all(ieee_is_finite(velocity))) return
        if (.not. all(ieee_is_finite(cosine0))) return
        if (.not. all(ieee_is_finite(sine0))) return
        if (velocity(1) < -1.0_dp .or. velocity(1) > -1.0_dp) return
        if (velocity(nodes) < 1.0_dp .or. velocity(nodes) > 1.0_dp) return
        if (any(sine0(0, :) < 0.0_dp .or. sine0(0, :) > 0.0_dp)) return
        do j = 1, nodes - 1
            if (velocity(j + 1) <= velocity(j)) return
        end do
        allocate (state%v(nodes), state%c0(0:modes, nodes), &
                  state%s0(0:modes, nodes), state%c(0:modes, nodes), &
                  state%s(0:modes, nodes), state%ci(0:modes, nodes), &
                  state%si(0:modes, nodes), state%cell_radius(nodes - 1), &
                  state%dc(0:modes), state%ds(0:modes), &
                  state%dci(0:modes), state%dsi(0:modes))
        state%v = velocity
        state%c0 = cosine0
        state%s0 = sine0
        state%initialized = .true.
        call advance_streaming(state, 0.0_dp, ok)
    end subroutine init_streaming

    subroutine advance_streaming(state, t, ok, certify)
        type(streaming_state_t), intent(inout) :: state
        real(dp), intent(in) :: t
        logical, intent(out) :: ok
        logical, intent(in), optional :: certify
        type(interval_t) :: a, b, bound, remainder, acc, next
        integer :: j, k

        ok = .false.
        state%prepared = .false.
        if (.not. state%initialized) return
        if (.not. ieee_is_finite(t)) return
        state%t = t
        do j = 1, size(state%v)
            do k = 0, size(state%c, 1) - 1
                call streaming_rotation(t, state%v(j), real(k, dp), &
                    state%c0(k, j), state%s0(k, j), state%c(k, j), state%s(k, j))
            end do
        end do
        if (.not. all(ieee_is_finite(state%c))) return
        if (.not. all(ieee_is_finite(state%s))) return
        ok = .true.
        if (present(certify)) then
            if (.not. certify) return
        end if
        ok = .false.
        do j = 1, size(state%v)
            do k = 0, size(state%c, 1) - 1
                call streaming_rotation_interval(interval(t), interval(state%v(j)), &
                    interval(real(k, dp)), interval(state%c0(k, j)), &
                    interval(state%s0(k, j)), a, b)
                if (.not. finite_interval(a) .or. .not. finite_interval(b)) return
                ! Include the untrusted candidate; this also covers its nodal
                ! discrepancy throughout a P1 cell by convex interpolation.
                state%ci(k, j) = hull(a, interval(state%c(k, j)))
                state%si(k, j) = hull(b, interval(state%s(k, j)))
            end do
        end do
        do j = 1, size(state%v) - 1
            acc = interval(0.0_dp)
            do k = 0, size(state%c, 1) - 1
                call streaming_cell_curvature_interval(interval(t), &
                    interval(real(k, dp)), interval(state%v(j)), &
                    interval(state%v(j + 1)), interval(state%c0(k, j)), &
                    interval(state%c0(k, j + 1)), interval(state%s0(k, j)), &
                    interval(state%s0(k, j + 1)), bound)
                if (.not. finite_interval(bound)) return
                ! Enclosure rounding can make the lower endpoint slightly negative.
                bound%lo = max(0.0_dp, bound%lo)
                remainder = linear_interpolation_remainder(interval(state%v(j)), &
                    interval(state%v(j + 1)), bound)
                if (.not. finite_interval(remainder)) return
                call streaming_add_interval(acc, remainder, next)
                acc = next
            end do
            if (.not. finite_interval(acc)) return
            state%cell_radius(j) = interval(0.0_dp, acc%hi)
        end do
        call prepare_density(state, ok)
        state%prepared = ok
    end subroutine advance_streaming

    subroutine prepare_density(state, ok)
        type(streaming_state_t), intent(inout) :: state
        logical, intent(out) :: ok
        type(interval_t) :: acc, next
        real(dp) :: fnext
        integer :: j, k

        ok = .false.
        do k = 0, size(state%c, 1) - 1
            state%dc(k) = 0.0_dp
            state%ds(k) = 0.0_dp
            state%dci(k) = interval(0.0_dp)
            state%dsi(k) = interval(0.0_dp)
            do j = 1, size(state%v) - 1
                call streaming_density_cell(state%v(j), state%v(j + 1), &
                    state%c(k, j), state%c(k, j + 1), state%dc(k), fnext)
                state%dc(k) = fnext
                call streaming_density_cell(state%v(j), state%v(j + 1), &
                    state%s(k, j), state%s(k, j + 1), state%ds(k), fnext)
                state%ds(k) = fnext
                call streaming_density_cell_interval(interval(state%v(j)), &
                    interval(state%v(j + 1)), state%ci(k, j), state%ci(k, j + 1), &
                    state%dci(k), next)
                state%dci(k) = next
                call streaming_density_cell_interval(interval(state%v(j)), &
                    interval(state%v(j + 1)), state%si(k, j), state%si(k, j + 1), &
                    state%dsi(k), next)
                state%dsi(k) = next
            end do
            if (.not. finite_interval(state%dci(k))) return
            if (.not. finite_interval(state%dsi(k))) return
        end do
        acc = interval(0.0_dp)
        do j = 1, size(state%v) - 1
            call streaming_density_radius_interval(interval(state%v(j)), &
                interval(state%v(j + 1)), state%cell_radius(j), acc, next)
            acc = next
        end do
        if (.not. finite_interval(acc)) return
        state%density_radius = interval(0.0_dp, acc%hi)
        ok = .true.
    end subroutine prepare_density

    subroutine evaluate_streaming(state, x, v, value, lower, upper, ok)
        type(streaming_state_t), intent(in) :: state
        real(dp), intent(in) :: x, v
        real(dp), intent(out) :: value, lower, upper
        logical, intent(out) :: ok
        type(interval_t) :: enclosure

        value = 0.0_dp
        lower = 0.0_dp
        upper = 0.0_dp
        call evaluate_streaming_interval(state, interval(x), interval(v), enclosure, ok)
        if (.not. ok) return
        call evaluate_center(state, x, v, .false., value)
        if (.not. ieee_is_finite(value)) then
            ok = .false.
            return
        end if
        enclosure = hull(enclosure, interval(value))
        lower = enclosure%lo
        upper = enclosure%hi
    end subroutine evaluate_streaming

    subroutine evaluate_streaming_interval(state, x, v, value, ok)
        type(streaming_state_t), intent(in) :: state
        type(interval_t), intent(in) :: x, v
        type(interval_t), intent(out) :: value
        logical, intent(out) :: ok

        ok = .false.
        value = interval(0.0_dp)
        if (.not. state%prepared) return
        call evaluate_packet(state, x, v, .false., value, ok)
    end subroutine evaluate_streaming_interval

    ! Strict nodal/P1 packet without the continuous velocity chord margin.
    ! It is intended for uniform-band construction with streaming_radius;
    ! the packet alone does not enclose the off-node continuous solution.
    subroutine evaluate_streaming_candidate_interval(state, x, v, value, ok)
        type(streaming_state_t), intent(in) :: state
        type(interval_t), intent(in) :: x, v
        type(interval_t), intent(out) :: value
        logical, intent(out) :: ok

        ok = .false.
        value = interval(0.0_dp)
        if (.not. state%prepared) return
        call evaluate_packet(state, x, v, .false., value, ok, include_radius=.false.)
    end subroutine evaluate_streaming_candidate_interval

    subroutine evaluate_streaming_initial_interval(state, x, v, value, ok)
        type(streaming_state_t), intent(in) :: state
        type(interval_t), intent(in) :: x, v
        type(interval_t), intent(out) :: value
        logical, intent(out) :: ok

        ok = .false.
        value = interval(0.0_dp)
        if (.not. state%initialized) return
        call evaluate_packet(state, x, v, .true., value, ok)
    end subroutine evaluate_streaming_initial_interval

    subroutine evaluate_packet(state, x, v, initial, value, ok, include_radius)
        type(streaming_state_t), intent(in) :: state
        type(interval_t), intent(in) :: x, v
        logical, intent(in) :: initial
        type(interval_t), intent(out) :: value
        logical, intent(out) :: ok
        logical, intent(in), optional :: include_radius
        type(interval_t) :: vc, a, b, al, ar, bl, br, acc, next, lower, upper
        integer :: j, k, first, last
        logical :: widen

        ok = .false.
        value = interval(0.0_dp)
        if (.not. finite_interval(x) .or. .not. finite_interval(v)) return
        if (v%lo < -1.0_dp .or. v%hi > 1.0_dp) return
        widen = .not. initial
        if (present(include_radius)) then
            if (.not. include_radius) widen = .false.
        end if
        first = locate(state%v, v%lo)
        last = locate(state%v, v%hi)
        do j = first, last
            vc = intersect(v, interval(state%v(j), state%v(j + 1)))
            if (is_empty(vc)) cycle
            acc = interval(0.0_dp)
            do k = 0, size(state%c0, 1) - 1
                if (initial) then
                    al = interval(state%c0(k, j)); ar = interval(state%c0(k, j + 1))
                    bl = interval(state%s0(k, j)); br = interval(state%s0(k, j + 1))
                else
                    al = state%ci(k, j); ar = state%ci(k, j + 1)
                    bl = state%si(k, j); br = state%si(k, j + 1)
                end if
                call streaming_interpolate_interval(vc, interval(state%v(j)), &
                    interval(state%v(j + 1)), al, ar, a)
                call streaming_interpolate_interval(vc, interval(state%v(j)), &
                    interval(state%v(j + 1)), bl, br, b)
                call streaming_mode_accumulate_interval(x, interval(real(k, dp)), &
                    a, b, acc, next)
                acc = next
            end do
            if (widen) then
                call streaming_enlarge_interval(acc, state%cell_radius(j), lower, upper)
                acc = hull(lower, upper)
            end if
            if (.not. finite_interval(acc)) return
            if (j == first) then
                value = acc
            else
                value = hull(value, acc)
            end if
        end do
        ok = finite_interval(value)
    end subroutine evaluate_packet

    subroutine evaluate_center(state, x, v, initial, value)
        type(streaming_state_t), intent(in) :: state
        real(dp), intent(in) :: x, v
        logical, intent(in) :: initial
        real(dp), intent(out) :: value
        real(dp) :: al, ar, bl, br, a, b, next
        integer :: j, k

        j = locate(state%v, v)
        value = 0.0_dp
        do k = 0, size(state%c0, 1) - 1
            if (initial) then
                al = state%c0(k, j); ar = state%c0(k, j + 1)
                bl = state%s0(k, j); br = state%s0(k, j + 1)
            else
                al = state%c(k, j); ar = state%c(k, j + 1)
                bl = state%s(k, j); br = state%s(k, j + 1)
            end if
            call streaming_interpolate(v, state%v(j), state%v(j + 1), al, ar, a)
            call streaming_interpolate(v, state%v(j), state%v(j + 1), bl, br, b)
            call streaming_mode_accumulate(x, real(k, dp), a, b, value, next)
            value = next
        end do
    end subroutine evaluate_center

    subroutine density_streaming(state, x, value, lower, upper, ok)
        type(streaming_state_t), intent(in) :: state
        real(dp), intent(in) :: x
        real(dp), intent(out) :: value, lower, upper
        logical, intent(out) :: ok
        real(dp) :: fnext
        type(interval_t) :: acc, next, lo, hi
        integer :: k

        ok = .false.
        value = 0.0_dp
        lower = 0.0_dp
        upper = 0.0_dp
        if (.not. state%prepared) return
        if (.not. ieee_is_finite(x)) return
        acc = interval(0.0_dp)
        do k = 0, size(state%dc) - 1
            call streaming_mode_accumulate(x, real(k, dp), state%dc(k), state%ds(k), &
                value, fnext)
            value = fnext
            call streaming_mode_accumulate_interval(interval(x), &
                interval(real(k, dp)), state%dci(k), state%dsi(k), acc, next)
            acc = next
        end do
        if (.not. ieee_is_finite(value)) return
        call streaming_enlarge_interval(acc, state%density_radius, lo, hi)
        acc = hull(hull(lo, hi), interval(value))
        if (.not. finite_interval(acc)) return
        lower = acc%lo
        upper = acc%hi
        ok = .true.
    end subroutine density_streaming

    subroutine density_streaming_chord(state, xl, xr, padding, ok)
        type(streaming_state_t), intent(in) :: state
        real(dp), intent(in) :: xl, xr
        real(dp), intent(out) :: padding
        logical, intent(out) :: ok
        type(interval_t) :: acc, next
        integer :: k

        padding = 0.0_dp
        ok = .false.
        if (.not. state%prepared) return
        if (.not. ieee_is_finite(xl) .or. .not. ieee_is_finite(xr)) return
        if (xr < xl) return
        acc = interval(0.0_dp)
        do k = 0, size(state%dc) - 1
            call streaming_x_curvature_interval(interval(real(k, dp)), &
                state%dci(k), state%dsi(k), acc, next)
            acc = next
        end do
        if (.not. finite_interval(acc)) return
        acc%lo = max(0.0_dp, acc%lo)
        next = linear_interpolation_remainder(interval(xl), interval(xr), acc)
        if (.not. finite_interval(next)) return
        padding = next%hi
        ok = .true.
    end subroutine density_streaming_chord

    subroutine density_coefficients_streaming(state, c, s, clo, chi, slo, shi, &
                                              radius, ok)
        type(streaming_state_t), intent(in) :: state
        real(dp), intent(out) :: c(0:), s(0:), clo(0:), chi(0:), slo(0:), shi(0:)
        real(dp), intent(out) :: radius
        logical, intent(out) :: ok
        integer :: n

        ok = .false.
        radius = 0.0_dp
        if (.not. state%prepared) return
        n = size(state%dc)
        if (size(c) /= n .or. size(s) /= n .or. size(clo) /= n .or. &
            size(chi) /= n .or. size(slo) /= n .or. size(shi) /= n) return
        c = state%dc; s = state%ds
        clo = state%dci%lo; chi = state%dci%hi
        slo = state%dsi%lo; shi = state%dsi%hi
        radius = state%density_radius%hi
        ok = .true.
    end subroutine density_coefficients_streaming

    subroutine streaming_radius(state, radius, ok)
        type(streaming_state_t), intent(in) :: state
        real(dp), intent(out) :: radius
        logical, intent(out) :: ok
        type(interval_t) :: bound
        integer :: j

        ! Velocity interpolation remainder only. Nodal rotation and evaluation
        ! arithmetic are separately enclosed by the strict sampling packet.
        ok = state%prepared
        radius = 0.0_dp
        if (.not. ok) return
        bound = state%cell_radius(1)
        do j = 2, size(state%cell_radius)
            bound = hull(bound, state%cell_radius(j))
        end do
        radius = bound%hi
    end subroutine streaming_radius

    subroutine streaming_shape(state, modes, nodes, ok)
        type(streaming_state_t), intent(in) :: state
        integer, intent(out) :: modes, nodes
        logical, intent(out) :: ok

        modes = 0
        nodes = 0
        ok = state%initialized
        if (.not. ok) return
        modes = size(state%c0, 1) - 1
        nodes = size(state%v)
    end subroutine streaming_shape

    subroutine get_streaming_initial(state, velocity, cosine0, sine0, ok)
        type(streaming_state_t), intent(in) :: state
        real(dp), intent(out) :: velocity(:), cosine0(0:, :), sine0(0:, :)
        logical, intent(out) :: ok

        ok = .false.
        if (.not. state%initialized) return
        if (size(velocity) /= size(state%v)) return
        if (any(shape(cosine0) /= shape(state%c0))) return
        if (any(shape(sine0) /= shape(state%s0))) return
        velocity = state%v
        cosine0 = state%c0
        sine0 = state%s0
        ok = .true.
    end subroutine get_streaming_initial

    subroutine init_streaming_particles(initial, nx, nv, particles, ok, global_velocity)
        type(streaming_state_t), intent(in) :: initial
        integer, intent(in) :: nx, nv
        type(streaming_particles_t), intent(out) :: particles
        logical, intent(out) :: ok
        logical, intent(in), optional :: global_velocity
        type(interval_t) :: measure_i, value_i
        real(dp) :: vl, vr, measure, value
        integer :: cells, n, i, j, cell, p, k
        logical :: global, valid

        ok = .false.
        if (.not. initial%initialized) return
        if (nx < 1 .or. nv < 1) return
        global = .false.
        if (present(global_velocity)) global = global_velocity
        cells = size(initial%v) - 1
        if (global) then
            ! A global midpoint rule may cross P1 knots only in this exactly
            ! v-independent admitted case; arbitrary P1 sources use cell loading.
            do j = 2, size(initial%v)
                do k = 0, size(initial%c0, 1) - 1
                    if (initial%c0(k, j) < initial%c0(k, 1) .or. &
                        initial%c0(k, j) > initial%c0(k, 1)) return
                    if (initial%s0(k, j) < initial%s0(k, 1) .or. &
                        initial%s0(k, j) > initial%s0(k, 1)) return
                end do
            end do
            cells = 1
        end if
        if (nx > huge(n)/nv) return
        n = nx*nv
        if (n > huge(n)/cells) return
        n = n*cells
        allocate (particles%grid(size(initial%v)), &
                  particles%c0(0:size(initial%c0, 1) - 1, size(initial%v)), &
                  particles%s0(0:size(initial%s0, 1) - 1, size(initial%v)), &
                  particles%x(n), particles%v(n), particles%w(n), &
                  particles%xi(n), particles%vi(n), particles%wi(n))
        particles%grid = initial%v
        particles%c0 = initial%c0
        particles%s0 = initial%s0
        particles%nx = nx
        particles%nv = nv
        particles%cells = cells
        particles%global_velocity = global
        p = 0
        do cell = 1, cells
            if (global) then
                vl = -1.0_dp; vr = 1.0_dp
            else
                vl = initial%v(cell); vr = initial%v(cell + 1)
            end if
            do j = 0, nv - 1
                do i = 0, nx - 1
                    p = p + 1
                    call streaming_particle_cell_midpoint(real(i, dp), real(nx, dp), &
                        real(j, dp), real(nv, dp), vl, vr, particles%x(p), &
                        particles%v(p), measure)
                    call streaming_particle_cell_midpoint_interval( &
                        interval(real(i, dp)), interval(real(nx, dp)), &
                        interval(real(j, dp)), &
                        interval(real(nv, dp)), interval(vl), interval(vr), &
                        particles%xi(p), particles%vi(p), measure_i)
                    particles%xi(p) = hull(particles%xi(p), interval(particles%x(p)))
                    particles%vi(p) = intersect(hull(particles%vi(p), &
                        interval(particles%v(p))), interval(-1.0_dp, 1.0_dp))
                    call evaluate_center(initial, particles%x(p), particles%v(p), &
                        .true., value)
                    call streaming_particle_weight(measure, value, particles%w(p))
                    call evaluate_streaming_initial_interval(initial, particles%xi(p), &
                        particles%vi(p), value_i, valid)
                    if (.not. valid) return
                    call streaming_particle_weight_interval(measure_i, value_i, &
                        particles%wi(p))
                    particles%wi(p) = hull(particles%wi(p), interval(particles%w(p)))
                    if (.not. finite_interval(particles%xi(p))) return
                    if (.not. finite_interval(particles%vi(p))) return
                    if (.not. finite_interval(particles%wi(p))) return
                    if (.not. ieee_is_finite(particles%w(p))) return
                end do
            end do
        end do
        particles%prepared = .true.
        ok = .true.
    end subroutine init_streaming_particles

    subroutine streaming_particle_info(particles, nx, nv, cells, global_velocity, &
                                       ok, modes)
        type(streaming_particles_t), intent(in) :: particles
        integer, intent(out) :: nx, nv, cells
        logical, intent(out) :: global_velocity, ok
        integer, intent(out), optional :: modes

        nx = 0; nv = 0; cells = 0
        global_velocity = .false.
        if (present(modes)) modes = 0
        ok = particles%prepared
        if (.not. ok) return
        nx = particles%nx; nv = particles%nv; cells = particles%cells
        global_velocity = particles%global_velocity
        if (present(modes)) modes = size(particles%c0, 1) - 1
    end subroutine streaming_particle_info

    subroutine get_streaming_particles(particles, x0, velocity, weight, ok, &
                                      x0_interval, velocity_interval, weight_interval)
        type(streaming_particles_t), intent(in) :: particles
        real(dp), intent(out) :: x0(:), velocity(:), weight(:)
        logical, intent(out) :: ok
        type(interval_t), intent(out), optional :: x0_interval(:), velocity_interval(:)
        type(interval_t), intent(out), optional :: weight_interval(:)
        integer :: n

        ok = .false.
        if (.not. particles%prepared) return
        n = size(particles%x)
        if (size(x0) /= n .or. size(velocity) /= n .or. size(weight) /= n) return
        if (present(x0_interval)) then
            if (size(x0_interval) /= n) return
        end if
        if (present(velocity_interval)) then
            if (size(velocity_interval) /= n) return
        end if
        if (present(weight_interval)) then
            if (size(weight_interval) /= n) return
        end if
        x0 = particles%x; velocity = particles%v; weight = particles%w
        if (present(x0_interval)) x0_interval = particles%xi
        if (present(velocity_interval)) velocity_interval = particles%vi
        if (present(weight_interval)) weight_interval = particles%wi
        ok = .true.
    end subroutine get_streaming_particles

    logical function streaming_particles_source_matches(initial, particles) result(ok)
        type(streaming_state_t), intent(in) :: initial
        type(streaming_particles_t), intent(in) :: particles

        ok = .false.
        if (.not. initial%initialized .or. .not. particles%prepared) return
        if (size(initial%v) /= size(particles%grid)) return
        if (any(shape(initial%c0) /= shape(particles%c0))) return
        if (any(shape(initial%s0) /= shape(particles%s0))) return
        if (any(initial%v < particles%grid .or. initial%v > particles%grid)) return
        if (any(initial%c0 < particles%c0 .or. initial%c0 > particles%c0)) return
        if (any(initial%s0 < particles%s0 .or. initial%s0 > particles%s0)) return
        ok = .true.
    end function streaming_particles_source_matches

    pure logical function finite_interval(a) result(ok)
        type(interval_t), intent(in) :: a

        ok = .false.
        if (is_empty(a)) return
        if (.not. ieee_is_finite(a%lo) .or. .not. ieee_is_finite(a%hi)) return
        ok = .true.
    end function finite_interval

    integer function locate(grid, value) result(cell)
        real(dp), intent(in) :: grid(:), value
        integer :: left, right, middle

        left = 1
        right = size(grid)
        do while (right - left > 1)
            middle = (left + right)/2
            if (grid(middle) <= value) then
                left = middle
            else
                right = middle
            end if
        end do
        cell = min(left, size(grid) - 1)
    end function locate

end module kin6d_streaming_comparison
