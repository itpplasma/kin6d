program test_streaming_comparison
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
    use fortnum_interval, only: interval_t, interval
    use kin6d_streaming_comparison, only: &
        streaming_state_t, init_streaming, advance_streaming, evaluate_streaming, &
        evaluate_streaming_interval, evaluate_streaming_candidate_interval, &
        streaming_radius, streaming_particles_t, &
        init_streaming_particles, get_streaming_particles, streaming_particle_info, &
        streaming_particles_source_matches, density_streaming, &
        density_coefficients_streaming, density_streaming_chord
    use kin6d_streaming_leaves, only: &
        streaming_rotation, streaming_rotation_jvp, streaming_rotation_vjp, &
        streaming_rotation_interval, streaming_characteristic, &
        streaming_characteristic_jvp, streaming_characteristic_vjp, &
        streaming_characteristic_interval
    implicit none

    call arbitrary_initial_data()
    call velocity_slice_chords()
    call refinement()
    call mean_and_zero_time()
    call derivative_oracles()
    call particle_oracles()
    call particle_quadrature()
    call density_oracles()
    call invalid_inputs()
    write (*, '(a)') &
        'Streaming characteristic, reconstruction and derivative oracles pass'

contains

    subroutine fixture(grid, a0, b0)
        real(dp), intent(out) :: grid(5), a0(0:3, 5), b0(0:3, 5)

        grid = [-1.0_dp, -0.625_dp, -0.125_dp, 0.5_dp, 1.0_dp]
        a0(0, :) = [0.75_dp, -0.25_dp, 1.5_dp, 0.375_dp, -0.5_dp]
        a0(1, :) = [0.25_dp, -0.5_dp, 0.75_dp, 0.5_dp, -0.25_dp]
        a0(2, :) = [-0.125_dp, 0.375_dp, 0.25_dp, -0.5_dp, 0.125_dp]
        a0(3, :) = [0.375_dp, 0.125_dp, -0.25_dp, 0.5_dp, 0.25_dp]
        b0(0, :) = 0.0_dp
        b0(1, :) = [-0.25_dp, 0.5_dp, 0.125_dp, -0.75_dp, 0.375_dp]
        b0(2, :) = [0.5_dp, -0.25_dp, 0.75_dp, 0.125_dp, -0.5_dp]
        b0(3, :) = [-0.5_dp, 0.25_dp, 0.125_dp, -0.125_dp, 0.375_dp]
    end subroutine fixture

    subroutine arbitrary_initial_data()
        type(streaming_state_t) :: state
        real(dp) :: grid(5), a0(0:3, 5), b0(0:3, 5)
        real(dp), parameter :: times(4) = [-1.25_dp, 0.0_dp, 0.875_dp, 3.125_dp]
        real(dp) :: x, velocity, value, lo, hi, radius
        real(qp) :: exact, candidate
        type(interval_t) :: box, xcell, vcell
        logical :: ok
        integer :: observation, j, r

        call fixture(grid, a0, b0)
        call init_streaming(state, grid, a0, b0, ok)
        call require(ok, 'arbitrary nonuniform initial packet')
        do observation = 1, size(times)
            ! Each call selects absolute time from the initial packet.
            call advance_streaming(state, times(observation), ok)
            call require(ok, 'absolute signed streaming time')
            call streaming_radius(state, radius, ok)
            call require(ok, 'streaming radius prepared')
            call require(radius >= 0.0_dp, 'negative streaming radius')
            do j = 1, size(grid) - 1
                do r = 0, 9
                    velocity = grid(j) + (grid(j + 1) - grid(j))*real(r, dp)/9
                    x = -1.375_dp + real(3*j + r, dp)/7
                    exact = characteristic_value(grid, a0, b0, x, velocity, &
                                                 times(observation))
                    candidate = nodal_reconstruction(grid, a0, b0, x, velocity, &
                                                     times(observation))
                    call evaluate_streaming(state, x, velocity, value, lo, hi, ok)
                    call require(ok, 'off-node streaming evaluation')
                    call close_to(value, candidate, 'P1 rotated-node reconstruction')
                    call contains(interval(lo, hi), exact, 'continuous characteristic')
                    call require(abs(exact - candidate) <= real(radius, qp) + &
                                 256*epsilon(1.0_qp)*(1 + abs(exact)), &
                                 'whole-cell interpolation error radius')
                end do
            end do
            xcell = interval(-0.875_dp, 0.75_dp)
            vcell = interval(-0.8125_dp, 0.8125_dp)
            call evaluate_streaming_interval(state, xcell, vcell, box, ok)
            call require(ok, 'whole-box spanning nonuniform velocity cells')
            do j = 0, 12
                x = xcell%lo + (xcell%hi - xcell%lo)*real(j, dp)/12
                do r = 0, 12
                    velocity = vcell%lo + (vcell%hi - vcell%lo)*real(r, dp)/12
                    exact = characteristic_value(grid, a0, b0, x, velocity, &
                                                 times(observation))
                    call contains(box, exact, 'whole-box characteristic containment')
                end do
            end do
        end do
        ! A spatial box wider than one period must contain interior extrema.
        call evaluate_streaming_interval(state, interval(-8.0_dp, 9.0_dp), &
                                         interval(-1.0_dp, 1.0_dp), box, ok)
        call require(ok, 'multi-period spatial box')
        do r = 0, 16
            x = -7.75_dp + real(r, dp)
            velocity = -1.0_dp + real(r, dp)/8
            exact = characteristic_value(grid, a0, b0, x, velocity, times(4))
            call contains(box, exact, 'multi-period characteristic containment')
        end do
    end subroutine arbitrary_initial_data

    subroutine velocity_slice_chords()
        type(streaming_state_t) :: state
        real(dp) :: grid(5), a0(0:3, 5), b0(0:3, 5), radius, velocity, x, time
        type(interval_t) :: vertices(5)
        real(qp) :: fraction, lower, upper, exact, candidate
        logical :: ok
        integer :: observation, j, r

        call fixture(grid, a0, b0)
        ! A wide cell with tiny neighbours exercises the shared-knot margin.
        grid = [-1.0_dp, -0.96875_dp, 0.75_dp, 0.78125_dp, 1.0_dp]
        call init_streaming(state, grid, a0, b0, ok)
        call require(ok, 'coarse velocity chord initial packet')
        x = 0.6875_dp
        do observation = 1, 3
            time = real(observation - 2, dp)*0.875_dp
            call advance_streaming(state, time, ok)
            call require(ok, 'coarse velocity chord evolution')
            call streaming_radius(state, radius, ok)
            call require(ok, 'uniform whole-velocity interpolation margin')
            do j = 1, size(grid)
                call evaluate_streaming_candidate_interval(state, interval(x), &
                                                           interval(grid(j)), &
                                                           vertices(j), ok)
                call require(ok, 'strict candidate vertex interval')
                candidate = nodal_reconstruction(grid, a0, b0, x, grid(j), time)
                call contains(vertices(j), candidate, 'candidate vertex reconstruction')
            end do
            do j = 1, size(grid) - 1
                do r = 0, 20
                    velocity = grid(j) + (grid(j + 1) - grid(j))*real(r, dp)/20
                    ! Match the exact real value of the sampled binary64 velocity.
                    fraction = (real(velocity, qp) - real(grid(j), qp))/ &
                               (real(grid(j + 1), qp) - real(grid(j), qp))
                    lower = (1 - fraction)*real(vertices(j)%lo, qp) + &
                            fraction*real(vertices(j + 1)%lo, qp) - real(radius, qp)
                    upper = (1 - fraction)*real(vertices(j)%hi, qp) + &
                            fraction*real(vertices(j + 1)%hi, qp) + real(radius, qp)
                    candidate = nodal_reconstruction(grid, a0, b0, x, velocity, time)
                    exact = characteristic_value(grid, a0, b0, x, velocity, time)
                    call require(candidate >= lower .and. candidate <= upper, &
                                 'connected candidate velocity chord')
                    call require(exact >= lower .and. exact <= upper, &
                                 'connected continuous velocity band including knots')
                end do
            end do
        end do
    end subroutine velocity_slice_chords

    subroutine refinement()
        type(streaming_state_t) :: state
        real(dp), allocatable :: grid(:), a0(:, :), b0(:, :)
        real(dp) :: radii(3), x, value, lo, hi
        real(qp) :: exact
        logical :: ok
        integer :: level, cells, j

        do level = 1, 3
            cells = 4*2**level
            allocate (grid(cells + 1), a0(0:3, cells + 1), b0(0:3, cells + 1))
            do j = 1, cells + 1
                grid(j) = -1.0_dp + 2.0_dp*real(j - 1, dp)/real(cells, dp)
            end do
            ! Multiple signed harmonics, unrelated to the demonstration profile.
            a0(0, :) = -0.375_dp
            a0(1, :) = 0.625_dp
            a0(2, :) = -0.25_dp
            a0(3, :) = 0.125_dp
            b0(0, :) = 0.0_dp
            b0(1, :) = -0.125_dp
            b0(2, :) = 0.375_dp
            b0(3, :) = -0.5_dp
            call init_streaming(state, grid, a0, b0, ok)
            call require(ok, 'refinement initial packet')
            call advance_streaming(state, -0.875_dp, ok)
            call require(ok, 'refinement evolution')
            call streaming_radius(state, radii(level), ok)
            call require(ok, 'refinement certificate')
            call require(radii(level) > 0.0_dp, 'nontrivial curvature radius')
            x = 0.6875_dp
            call evaluate_streaming(state, x, 0.0_dp, value, lo, hi, ok)
            call require(ok, 'refinement node query')
            exact = characteristic_value(grid, a0, b0, x, 0.0_dp, -0.875_dp)
            call contains(interval(lo, hi), exact, 'refinement continuous oracle')
            deallocate (grid, a0, b0)
        end do
        call require(radii(2) <= 0.26_dp*radii(1), 'first h-squared refinement')
        call require(radii(3) <= 0.26_dp*radii(2), 'second h-squared refinement')
    end subroutine refinement

    subroutine mean_and_zero_time()
        type(streaming_state_t) :: state
        real(dp) :: grid(3), a0(0:0, 3), b0(0:0, 3), radius, value, lo, hi
        real(qp) :: exact
        logical :: ok
        integer :: j

        grid = [-1.0_dp, -0.25_dp, 1.0_dp]
        a0(0, :) = [-0.5_dp, 0.875_dp, 0.25_dp]
        b0 = 0.0_dp
        call init_streaming(state, grid, a0, b0, ok)
        call require(ok, 'mean-only initialization')
        do j = 1, 3
            call advance_streaming(state, real(j - 2, dp)*1.0e4_dp, ok)
            call require(ok, 'mean-only signed time')
            call streaming_radius(state, radius, ok)
            call require(ok, 'mean-only certificate')
            call require(radius < 1.0e-12_dp, 'mean must have zero curvature')
            call evaluate_streaming(state, 2.0_dp, 0.125_dp, value, lo, hi, ok)
            call require(ok, 'mean-only evaluation')
            exact = characteristic_value(grid, a0, b0, 2.0_dp, 0.125_dp, 0.0_dp)
            call contains(interval(lo, hi), exact, 'invariant mean-only interpolation')
        end do
    end subroutine mean_and_zero_time

    subroutine derivative_oracles()
        real(dp), parameter :: v = -0.625_dp, k = 3.0_dp
        real(dp), parameter :: a0 = -0.375_dp, b0 = 0.875_dp
        real(dp), parameter :: dt = 0.25_dp, dv = -0.125_dp
        real(dp), parameter :: da0 = 0.0625_dp, db0 = -0.25_dp
        real(dp), parameter :: ba = -0.5_dp, bb = 0.375_dp
        real(dp) :: t, a, b, da, db, gt, gv, gk, ga0, gb0, back_a, back_b
        real(qp) :: theta, c, s, aq, bq, dtheta, ja, jb, phase_dual, lhs, rhs
        type(interval_t) :: ai, bi
        integer :: j

        do j = 1, 3
            t = real(j - 2, dp)*1.375_dp
            theta = real(k, qp)*real(v, qp)*real(t, qp)
            c = cos(theta); s = sin(theta)
            aq = real(a0, qp)*c - real(b0, qp)*s
            bq = real(a0, qp)*s + real(b0, qp)*c
            call streaming_rotation(t, v, k, a0, b0, a, b)
            call close_to(a, aq, 'independent signed rotation cosine')
            call close_to(b, bq, 'independent signed rotation sine')
            call require(abs(a*a + b*b - a0*a0 - b0*b0) < 1.0e-13_dp, &
                         'mode norm conserved by streaming')
            call streaming_rotation(-t, v, k, a, b, back_a, back_b)
            call close_to(back_a, real(a0, qp), 'rotation time reversal cosine')
            call close_to(back_b, real(b0, qp), 'rotation time reversal sine')
            call streaming_rotation_interval(interval(t), interval(v), interval(k), &
                                             interval(a0), interval(b0), ai, bi)
            call contains(ai, aq, 'rigorous rotated cosine')
            call contains(bi, bq, 'rigorous rotated sine')

            ! Integer mode choice is fixed: dk=0 throughout the physical direction.
            dtheta = real(k, qp)*(real(v, qp)*real(dt, qp) + &
                                  real(t, qp)*real(dv, qp))
            ja = real(da0, qp)*c - real(db0, qp)*s - bq*dtheta
            jb = real(da0, qp)*s + real(db0, qp)*c + aq*dtheta
            call streaming_rotation_jvp(t, v, k, a0, b0, dt, dv, 0.0_dp, &
                                        da0, db0, da, db)
            call close_to(da, ja, 'analytical streaming JVP cosine')
            call close_to(db, jb, 'analytical streaming JVP sine')
            call streaming_rotation_vjp(t, v, k, a0, b0, ba, bb, &
                                        gt, gv, gk, ga0, gb0)
            phase_dual = -real(ba, qp)*bq + real(bb, qp)*aq
            call close_to(gt, real(k*v, qp)*phase_dual, 'rotation time VJP')
            call close_to(gv, real(k*t, qp)*phase_dual, 'rotation velocity VJP')
            call close_to(ga0, real(ba, qp)*c + real(bb, qp)*s, 'rotation input VJP')
            call close_to(gb0, -real(ba, qp)*s + real(bb, qp)*c, 'rotation sine VJP')
            lhs = real(ba, qp)*real(da, qp) + real(bb, qp)*real(db, qp)
            rhs = real(gt, qp)*real(dt, qp) + real(gv, qp)*real(dv, qp) + &
                  real(ga0, qp)*real(da0, qp) + real(gb0, qp)*real(db0, qp)
            call require(abs(lhs - rhs) < 1.0e-13_qp, 'raw leaf JVP/VJP pairing')
        end do
    end subroutine derivative_oracles

    subroutine particle_oracles()
        real(dp) :: x, dx, gx0, gv, gt, t, v, x0, bx
        real(qp) :: exact, tangent, lhs, rhs
        type(interval_t) :: enclosure
        integer :: j

        ! These inputs cross both periodic seams. The mathematical map is unwrapped.
        x0 = 7.0_dp; v = -0.625_dp; bx = -0.375_dp
        do j = 1, 3
            t = real(j - 2, dp)*15.0_dp
            exact = real(x0, qp) + real(v, qp)*real(t, qp)
            call streaming_characteristic(x0, v, t, x)
            call close_to(x, exact, 'unwrapped particle trajectory')
            call streaming_characteristic_interval(interval(x0), interval(v), &
                                                   interval(t), enclosure)
            call contains(enclosure, exact, 'particle trajectory interval')
            tangent = 0.25_qp + real(t, qp)*0.125_qp + real(v, qp)*(-0.5_qp)
            call streaming_characteristic_jvp(x0, v, t, 0.25_dp, 0.125_dp, &
                                              -0.5_dp, dx)
            call close_to(dx, tangent, 'analytical particle JVP')
            call streaming_characteristic_vjp(x0, v, t, bx, gx0, gv, gt)
            call close_to(gx0, real(bx, qp), 'particle initial-position VJP')
            call close_to(gv, real(bx, qp)*real(t, qp), 'particle velocity VJP')
            call close_to(gt, real(bx, qp)*real(v, qp), 'particle time VJP')
            lhs = real(bx, qp)*real(dx, qp)
            rhs = real(gx0, qp)*0.25_qp + real(gv, qp)*0.125_qp - real(gt, qp)*0.5_qp
            call require(abs(lhs - rhs) < 1.0e-13_qp, 'particle JVP/VJP pairing')
        end do
        call streaming_characteristic_interval(interval(-1.0_dp, 0.25_dp), &
                                               interval(-0.75_dp, 1.0_dp), &
                                               interval(-2.0_dp, 1.0_dp), enclosure)
        do j = 0, 12
            exact = (-1.0_qp + 1.25_qp*real(j, qp)/12) + &
                    (-0.75_qp + 1.75_qp*real(j, qp)/12)* &
                    (-2.0_qp + 3.0_qp*real(12 - j, qp)/12)
            call contains(enclosure, exact, 'uncertain particle parameter box')
        end do
    end subroutine particle_oracles

    subroutine particle_quadrature()
        integer, parameter :: nx = 8, nv = 3, cells = 4, count = nx*nv*cells
        type(streaming_state_t) :: state, other
        type(streaming_particles_t) :: particles
        real(dp) :: grid(5), a0(0:3, 5), b0(0:3, 5)
        real(dp) :: x0(count), velocity(count), weight(count), changed(0:3, 5)
        type(interval_t) :: xi(count), vi(count), wi(count)
        real(qp) :: piq, xq, vq, phase_volume, exact_weight
        integer :: reported_nx, reported_nv, reported_cells, p, j, ix, iv
        integer :: seen(nx, nv, cells)
        logical :: ok, global_velocity

        call fixture(grid, a0, b0)
        call init_streaming(state, grid, a0, b0, ok)
        call require(ok, 'particle initial data')
        call init_streaming_particles(state, nx, nv, particles, ok)
        call require(ok, 'cell-aligned deterministic particle initialization')
        call streaming_particle_info(particles, reported_nx, reported_nv, &
                                     reported_cells, global_velocity, ok)
        call require(ok, 'particle packet layout information')
        call require(reported_nx == nx, 'particle x count')
        call require(reported_nv == nv, 'particle cell-velocity count')
        call require(reported_cells == cells, 'particle cell count')
        call require(.not. global_velocity, 'nonuniform input uses cell-aligned rule')
        call get_streaming_particles(particles, x0, velocity, weight, ok, xi, vi, wi)
        call require(ok, 'particle packet and rigorous weight extraction')
        piq = acos(-1.0_qp)
        seen = 0
        do p = 1, count
            j = oracle_cell(grid, velocity(p))
            ix = nint(real(x0(p), qp)*nx/(2*piq) + 0.5_qp)
            iv = nint((real(velocity(p), qp) - real(grid(j), qp))*nv/ &
                      (real(grid(j + 1), qp) - real(grid(j), qp)) + 0.5_qp)
            call require(ix >= 1 .and. ix <= nx, 'particle spatial midpoint index')
            call require(iv >= 1 .and. iv <= nv, 'particle velocity midpoint index')
            seen(ix, iv, j) = seen(ix, iv, j) + 1
            xq = (real(ix, qp) - 0.5_qp)*2*piq/nx
            vq = real(grid(j), qp) + (real(iv, qp) - 0.5_qp)* &
                 (real(grid(j + 1), qp) - real(grid(j), qp))/nv
            phase_volume = (2*piq/nx)* &
                           (real(grid(j + 1), qp) - real(grid(j), qp))/nv
            exact_weight = initial_value_qp(grid, a0, b0, xq, vq)*phase_volume
            call close_to(x0(p), xq, 'particle spatial midpoint oracle')
            call close_to(velocity(p), vq, 'particle velocity midpoint oracle')
            call close_to(weight(p), exact_weight, 'particle signed quadrature weight')
            call contains(xi(p), xq, 'particle midpoint position interval')
            call contains(vi(p), vq, 'particle midpoint velocity interval')
            call contains(wi(p), exact_weight, &
                          'particle initial-function weight interval')
        end do
        call require(all(seen == 1), 'complete nonduplicated midpoint particle lattice')
        call require(any(weight < 0.0_dp), &
                     'signed initial data retain negative weights')
        call require(streaming_particles_source_matches(state, particles), &
                     'particle source association')
        call advance_streaming(state, -0.875_dp, ok)
        call require(ok, 'particle source evolution')
        call require(streaming_particles_source_matches(state, particles), &
                     'particle packet remains tied to immutable initial data')
        changed = a0; changed(2, 3) = changed(2, 3) + 0.125_dp
        call init_streaming(other, grid, changed, b0, ok)
        call require(ok, 'changed particle source fixture')
        call require(.not. streaming_particles_source_matches(other, particles), &
                     'different initial source rejected')
        call get_streaming_particles(particles, x0(:count - 1), velocity, weight, ok)
        call require(.not. ok, 'wrong particle output extent rejected')
        call init_streaming_particles(other, nx, nv, particles, ok, &
                                      global_velocity=.true.)
        call require(.not. ok, 'global velocity rule rejected for nonconstant source')
        call init_streaming_particles(state, 0, nv, particles, ok)
        call require(.not. ok, 'zero spatial particle count rejected')
        call init_streaming_particles(state, nx, -1, particles, ok)
        call require(.not. ok, 'negative velocity particle count rejected')
        do j = 2, size(grid)
            a0(:, j) = a0(:, 1)
            b0(:, j) = b0(:, 1)
        end do
        call init_streaming(other, grid, a0, b0, ok)
        call require(ok, 'velocity-independent global particle source')
        call init_streaming_particles(other, nx, nv, particles, ok, &
                                      global_velocity=.true.)
        call require(ok, 'admitted global particle velocity rule')
        call streaming_particle_info(particles, reported_nx, reported_nv, &
                                     reported_cells, global_velocity, ok)
        call require(ok, 'global particle packet information')
        call require(global_velocity, 'global particle rule identified')
        call require(reported_cells == 1, 'global rule has one quadrature domain')
        call get_streaming_particles(particles, x0(:nx*nv), velocity(:nx*nv), &
                                     weight(:nx*nv), ok, xi(:nx*nv), &
                                     vi(:nx*nv), wi(:nx*nv))
        call require(ok, 'global particle packet extraction')
        seen = 0
        do p = 1, nx*nv
            ix = nint(real(x0(p), qp)*nx/(2*piq) + 0.5_qp)
            iv = nint((real(velocity(p), qp) + 1)*nv/2 + 0.5_qp)
            call require(ix >= 1 .and. ix <= nx, 'global spatial midpoint index')
            call require(iv >= 1 .and. iv <= nv, 'global velocity midpoint index')
            seen(ix, iv, 1) = seen(ix, iv, 1) + 1
            xq = (real(ix, qp) - 0.5_qp)*2*piq/nx
            vq = -1 + (real(iv, qp) - 0.5_qp)*2/nv
            exact_weight = initial_value_qp(grid, a0, b0, xq, vq)*(2*piq/nx)*(2.0_qp/nv)
            call close_to(weight(p), exact_weight, 'global particle quadrature weight')
            call contains(xi(p), xq, 'global spatial midpoint interval')
            call contains(vi(p), vq, 'global velocity midpoint interval')
            call contains(wi(p), exact_weight, 'global particle weight interval')
        end do
        call require(all(seen(:, :, 1) == 1), 'complete global particle lattice')
    end subroutine particle_quadrature

    subroutine density_oracles()
        type(streaming_state_t) :: state
        real(dp) :: grid(5), a0(0:3, 5), b0(0:3, 5), time, x, value, lo, hi
        real(dp) :: c(0:3), s(0:3), clo(0:3), chi(0:3), slo(0:3), shi(0:3), radius
        real(dp) :: endpoint_lo(2), endpoint_hi(2), padding, fraction
        complex(qp) :: exact_coefficient, candidate_coefficient, phase
        real(qp) :: exact
        logical :: ok
        integer :: observation, j, k

        call fixture(grid, a0, b0)
        call init_streaming(state, grid, a0, b0, ok)
        call require(ok, 'density source initialization')
        do observation = 1, 3
            time = real(observation - 2, dp)*0.875_dp
            call advance_streaming(state, time, ok)
            call require(ok, 'density signed evolution')
            call density_coefficients_streaming(state, c, s, clo, chi, slo, shi, &
                                                radius, ok)
            call require(ok, 'density coefficient extraction')
            do k = 0, 3
                exact_coefficient = density_coefficient(grid, a0, b0, time, k, .false.)
                candidate_coefficient = &
                    density_coefficient(grid, a0, b0, time, k, .true.)
                call close_to(c(k), real(candidate_coefficient, qp), &
                              'candidate density cosine integral')
                call close_to(s(k), -aimag(candidate_coefficient), &
                              'candidate density sine integral')
                call contains(interval(clo(k), chi(k)), &
                              real(candidate_coefficient, qp), &
                              'density cosine arithmetic')
                call contains(interval(slo(k), shi(k)), &
                              -aimag(candidate_coefficient), 'density sine arithmetic')
                call contains(interval(clo(k) - radius, chi(k) + radius), &
                              real(exact_coefficient, qp), 'continuous density cosine')
                call contains(interval(slo(k) - radius, shi(k) + radius), &
                              -aimag(exact_coefficient), 'continuous density sine')
            end do
            do j = 0, 16
                x = -0.875_dp + real(j, dp)/4
                exact = 0.0_qp
                do k = 0, 3
                    exact_coefficient = &
                        density_coefficient(grid, a0, b0, time, k, .false.)
                    phase = exp(cmplx(0.0_qp, real(k, qp)*real(x, qp), qp))
                    exact = exact + real(exact_coefficient*phase, qp)
                end do
                call density_streaming(state, x, value, lo, hi, ok)
                call require(ok, 'normalized physical density query')
                call contains(interval(lo, hi), exact, &
                              'continuous density integral oracle')
            end do
        end do
        ! A coarse spatial chord must cover the interior cosine maximum at x=0.
        a0 = 0.0_dp; b0 = 0.0_dp; a0(1, :) = 1.0_dp
        call init_streaming(state, grid, a0, b0, ok)
        call require(ok, 'density interior-maximum fixture')
        call advance_streaming(state, 0.0_dp, ok)
        call require(ok, 'density initial-time preparation')
        call density_streaming(state, -1.0_dp, value, &
                               endpoint_lo(1), endpoint_hi(1), ok)
        call require(ok, 'left density chord endpoint')
        call density_streaming(state, 1.0_dp, value, endpoint_lo(2), endpoint_hi(2), ok)
        call require(ok, 'right density chord endpoint')
        call density_streaming_chord(state, -1.0_dp, 1.0_dp, padding, ok)
        call require(ok, 'whole-spatial-cell density chord padding')
        call require(padding > 0.4_dp, 'density chord cannot use endpoint bands alone')
        do j = 0, 20
            fraction = real(j, dp)/20
            x = -1.0_dp + 2*fraction
            lo = (1 - fraction)*endpoint_lo(1) + fraction*endpoint_lo(2) - padding
            hi = (1 - fraction)*endpoint_hi(1) + fraction*endpoint_hi(2) + padding
            call contains(interval(lo, hi), cos(real(x, qp)), &
                          'density chord contains interior continuous curve')
        end do
    end subroutine density_oracles

    function density_coefficient(grid, a0, b0, time, mode, candidate) result(value)
        ! TEST ONLY: exact P1 characteristic integral, using its complex antiderivative.
        ! For candidate=true, use independently rotated-node trapezoidal integration.
        real(dp), intent(in) :: grid(:), a0(0:, :), b0(0:, :), time
        integer, intent(in) :: mode
        logical, intent(in) :: candidate
        complex(qp) :: value, zl, zr, slope, p, left, right
        real(qp) :: vl, vr
        integer :: j

        p = cmplx(0.0_qp, -real(mode, qp)*real(time, qp), qp)
        value = cmplx(0.0_qp, 0.0_qp, qp)
        do j = 1, size(grid) - 1
            vl = real(grid(j), qp); vr = real(grid(j + 1), qp)
            zl = cmplx(real(a0(mode, j), qp), -real(b0(mode, j), qp), qp)
            zr = cmplx(real(a0(mode, j + 1), qp), -real(b0(mode, j + 1), qp), qp)
            if (candidate) then
                value = value + (vr - vl)*(zl*exp(p*vl) + zr*exp(p*vr))/4
            else if (mode == 0 .or. abs(time) < tiny(time)) then
                value = value + (vr - vl)*(zl + zr)/4
            else
                slope = (zr - zl)/(vr - vl)
                left = exp(p*vl)*(zl/p - slope/p**2)
                right = exp(p*vr)*(zr/p - slope/p**2)
                value = value + (right - left)/2
            end if
        end do
    end function density_coefficient

    subroutine invalid_inputs()
        type(streaming_state_t) :: state, empty
        real(dp) :: grid(5), a0(0:3, 5), b0(0:3, 5), bad_grid(5), bad_a(0:3, 5)
        real(dp) :: nan, value, lo, hi, radius
        type(interval_t) :: box, reversed
        logical :: ok

        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        call fixture(grid, a0, b0)
        call advance_streaming(empty, 0.0_dp, ok)
        call require(.not. ok, 'uninitialized evolution rejected')
        call evaluate_streaming(empty, 0.0_dp, 0.0_dp, value, lo, hi, ok)
        call require(.not. ok, 'uninitialized evaluation rejected')
        call streaming_radius(empty, radius, ok)
        call require(.not. ok, 'uninitialized radius rejected')
        bad_grid = grid; bad_grid(3) = bad_grid(2)
        call init_streaming(state, bad_grid, a0, b0, ok)
        call require(.not. ok, 'collapsed velocity cell rejected')
        bad_grid = grid; bad_grid(1) = -2.0_dp
        call init_streaming(state, bad_grid, a0, b0, ok)
        call require(.not. ok, 'wrong velocity endpoint rejected')
        bad_grid = grid; bad_grid(4) = nan
        call init_streaming(state, bad_grid, a0, b0, ok)
        call require(.not. ok, 'NaN velocity rejected')
        bad_a = a0; bad_a(2, 3) = nan
        call init_streaming(state, grid, bad_a, b0, ok)
        call require(.not. ok, 'NaN initial coefficient rejected')
        b0(0, 2) = 0.25_dp
        call init_streaming(state, grid, a0, b0, ok)
        call require(.not. ok, 'nonzero mean-sine rejected')
        b0(0, 2) = 0.0_dp
        call init_streaming(state, grid, a0(:, :4), b0(:, :4), ok)
        call require(.not. ok, 'wrong coefficient count rejected')
        call init_streaming(state, grid, a0, b0(:2, :), ok)
        call require(.not. ok, 'inconsistent mode counts rejected')
        call init_streaming(state, grid(:1), a0(:, :1), b0(:, :1), ok)
        call require(.not. ok, 'one velocity node rejected')
        call init_streaming(state, grid, a0, b0, ok)
        call require(ok, 'valid packet after failures')
        call advance_streaming(state, nan, ok)
        call require(.not. ok, 'NaN observation time rejected')
        call advance_streaming(state, 0.5_dp, ok)
        call require(ok, 'valid time after rejected time')
        call advance_streaming(state, 0.875_dp, ok, certify=.false.)
        call require(ok, 'candidate-only evolution succeeds')
        call streaming_radius(state, radius, ok)
        call require(.not. ok, 'candidate-only radius remains unavailable')
        call evaluate_streaming(state, 0.0_dp, 0.0_dp, value, lo, hi, ok)
        call require(.not. ok, 'candidate-only bounded query rejected')
        call density_streaming(state, 0.0_dp, value, lo, hi, ok)
        call require(.not. ok, 'candidate-only density bound rejected')
        call advance_streaming(state, 0.5_dp, ok)
        call require(ok, 'certified preparation restored after candidate-only step')
        ! Finite x can overflow k*x for the active k=2 and k=3 harmonics.
        call evaluate_streaming(state, 0.75_dp*huge(0.0_dp), 0.25_dp, &
                                value, lo, hi, ok)
        call require(.not. ok, 'nonfinite point candidate from finite phase rejected')
        call density_streaming(state, 0.75_dp*huge(0.0_dp), value, lo, hi, ok)
        call require(.not. ok, 'nonfinite density candidate from finite phase rejected')
        call evaluate_streaming(state, nan, 0.0_dp, value, lo, hi, ok)
        call require(.not. ok, 'NaN position rejected')
        call evaluate_streaming(state, 0.0_dp, nan, value, lo, hi, ok)
        call require(.not. ok, 'NaN sampling velocity rejected')
        call evaluate_streaming(state, 0.0_dp, 1.125_dp, value, lo, hi, ok)
        call require(.not. ok, 'sampling beyond velocity domain rejected')
        call evaluate_streaming_interval(state, interval(0.0_dp), &
                                         interval(-1.125_dp, 0.5_dp), box, ok)
        call require(.not. ok, 'interval beyond velocity domain rejected')
        reversed%lo = 1.0_dp; reversed%hi = -1.0_dp
        call evaluate_streaming_interval(state, reversed, interval(0.0_dp), box, ok)
        call require(.not. ok, 'reversed spatial interval rejected')
    end subroutine invalid_inputs

    function characteristic_value(grid, a0, b0, x, velocity, time) result(value)
        ! TEST ONLY: advect the supplied piecewise-linear initial function.
        ! The binary64 input bytes denote exact reals; no analytic-profile shortcut.
        real(dp), intent(in) :: grid(:), a0(0:, :), b0(0:, :), x, velocity, time
        real(qp) :: value, shifted_x

        shifted_x = real(x, qp) - real(velocity, qp)*real(time, qp)
        value = initial_value_qp(grid, a0, b0, shifted_x, real(velocity, qp))
    end function characteristic_value

    function initial_value_qp(grid, a0, b0, xq, vq) result(value)
        real(dp), intent(in) :: grid(:), a0(0:, :), b0(0:, :)
        real(qp), intent(in) :: xq, vq
        real(qp) :: value, fraction, aq, bq
        complex(qp) :: coefficient, exponential
        integer :: j, k

        j = oracle_cell(grid, real(vq, dp))
        fraction = (vq - real(grid(j), qp))/ &
                   (real(grid(j + 1), qp) - real(grid(j), qp))
        value = 0.0_qp
        do k = 0, size(a0, 1) - 1
            aq = (1 - fraction)*real(a0(k, j), qp) + fraction*real(a0(k, j + 1), qp)
            bq = (1 - fraction)*real(b0(k, j), qp) + fraction*real(b0(k, j + 1), qp)
            coefficient = cmplx(aq, -bq, qp)
            exponential = exp(cmplx(0.0_qp, real(k, qp)*xq, qp))
            value = value + real(coefficient*exponential, qp)
        end do
    end function initial_value_qp

    function nodal_reconstruction(grid, a0, b0, x, velocity, time) result(value)
        ! TEST ONLY: independently reconstruct the rotated nodal candidate.
        real(dp), intent(in) :: grid(:), a0(0:, :), b0(0:, :), x, velocity, time
        real(qp) :: value, fraction
        complex(qp) :: zl, zr, coefficient, exponential
        integer :: j, k

        j = oracle_cell(grid, velocity)
        fraction = (real(velocity, qp) - real(grid(j), qp))/ &
                   (real(grid(j + 1), qp) - real(grid(j), qp))
        value = 0.0_qp
        do k = 0, size(a0, 1) - 1
            zl = cmplx(real(a0(k, j), qp), -real(b0(k, j), qp), qp)* &
                 exp(cmplx(0.0_qp, -real(k, qp)*real(grid(j), qp)*real(time, qp), qp))
            zr = cmplx(real(a0(k, j + 1), qp), -real(b0(k, j + 1), qp), qp)* &
                 exp(cmplx(0.0_qp, &
                           -real(k, qp)*real(grid(j + 1), qp)*real(time, qp), qp))
            coefficient = (1 - fraction)*zl + fraction*zr
            exponential = exp(cmplx(0.0_qp, real(k, qp)*real(x, qp), qp))
            value = value + real(coefficient*exponential, qp)
        end do
    end function nodal_reconstruction

    integer function oracle_cell(grid, velocity) result(cell)
        real(dp), intent(in) :: grid(:), velocity
        integer :: j

        cell = size(grid) - 1
        do j = 1, size(grid) - 1
            if (velocity <= grid(j + 1)) then
                cell = j
                return
            end if
        end do
    end function oracle_cell

    subroutine close_to(value, exact, description)
        real(dp), intent(in) :: value
        real(qp), intent(in) :: exact
        character(*), intent(in) :: description

        call require(ieee_is_finite(value), description//' finite value')
        call require(abs(real(value, qp) - exact) <= &
                     128*real(epsilon(1.0_dp), qp)*(1 + abs(exact)), description)
    end subroutine close_to

    subroutine contains(enclosure, exact, description)
        type(interval_t), intent(in) :: enclosure
        real(qp), intent(in) :: exact
        character(*), intent(in) :: description

        call require(ieee_is_finite(enclosure%lo), description//' finite lower')
        call require(ieee_is_finite(enclosure%hi), description//' finite upper')
        call require(enclosure%lo <= enclosure%hi, description//' ordered interval')
        call require(real(enclosure%lo, qp) <= exact, description//' lower endpoint')
        call require(real(enclosure%hi, qp) >= exact, description//' upper endpoint')
    end subroutine contains

    subroutine require(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        if (.not. condition) then
            write (*, '(a)') 'FAILED: '//description
            error stop 1
        end if
    end subroutine require
end program test_streaming_comparison
