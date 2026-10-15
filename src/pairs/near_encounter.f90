!> Encounter-resolved near/far splitting of the periodic Coulomb Hamiltonian.
!>
!> H = T + U_far + U_near with U_near = sum_pairs S(r) g erfc(alpha r)/r, S a
!> C2 switch (1 below r1, 0 above r2 = 2 r1; leaf near_switch_pair). The sum
!> U_far + U_near is exactly the bare Ewald/PME interaction: nothing is
!> softened. One step dt is the Strang composition
!>   kick_far(dt/2) o flow_near(dt) o kick_far(dt/2),
!> where flow_near integrates T + U_near: particles without a candidate near
!> pair drift exactly; each cluster of candidate pairs (straight-line closest
!> approach over [0,dt] below 1.5 r2) takes M equal KDK substeps resolving
!> its two-body periapsis Kepler, crossing and switch-crossing times (eta).
!> Pair radii r1 = max(speed tau, (2|g| tau^2)^(1/3)), tau = dt/eta_outer,
!> keep the far force smooth on the outer step. Exactness of the cluster
!> decomposition is checked after the flow: every particle's deviation from
!> its straight line must keep all non-candidate pairs outside r2, otherwise
!> the pair is added and the flow repeated (fail closed after 4 rounds).
!> Fixed particles (an impurity) are external sources inside clusters and
!> never merge clusters. Derivative: tangent of the same flow with the
!> cluster structure and substep counts held at the primal values.
module kin6d_near_encounter
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR, &
                              status_set
    use kin6d_coulomb_leaves, only: ewald_real_pair, ewald_real_pair_jvp
    use kin6d_ocp_leaves, only: near_switch_pair, near_switch_pair_jvp
    implicit none
    private
    public :: near_split_t, near_t, near_radii, near_pair, near_pair_jvp
    public :: near_setup, near_flow

    type :: near_split_t
        logical :: active = .false.
        real(dp) :: base = 0.0_dp, tau2 = 0.0_dp
    end type

    !> Scratch of the near flow, sized once (pair list grows by doubling).
    type :: near_t
        type(near_split_t) :: split
        real(dp) :: dt = 0.0_dp, eta = 0.05_dp
        integer :: n = 0, pairs = 0, max_substeps = 100000
        integer, allocatable :: pair(:, :), head(:), next(:), parent(:)
        integer, allocatable :: order(:), first(:), pair_order(:), pair_first(:)
        integer, allocatable :: cursor(:)
        real(dp), allocatable :: x0(:), v0(:), dx0(:), dv0(:), dev(:), f(:), df(:)
        logical, allocatable :: active(:)
        !> Statistics of the last flow: clusters, active particles, the
        !> largest substep count and verification rounds.
        integer :: clusters = 0, members = 0, deepest = 0, rounds = 0
    end type
contains
    !> Switch radii of a pair of coupling g.
    pure subroutine near_radii(split, g, r1, r2)
        type(near_split_t), intent(in) :: split
        real(dp), intent(in) :: g
        real(dp), intent(out) :: r1, r2
        r1 = max(split%base, (2*abs(g)*split%tau2)**(1.0_dp/3))
        r2 = 2*r1
    end subroutine

    !> Near part of the real-space Ewald pair at separation d (force on the
    !> particle at d, potential); zero beyond r2.
    pure subroutine near_pair(split, d, alpha, g, f, u)
        type(near_split_t), intent(in) :: split
        real(dp), intent(in) :: d(3), alpha, g
        real(dp), intent(out) :: f(3), u
        real(dp) :: r, r1, r2
        f = 0; u = 0
        call near_radii(split, g, r1, r2)
        r = sqrt(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
        if (r >= r2) return
        if (r <= r1) then
            call ewald_real_pair(d(1), d(2), d(3), alpha, g, f(1), f(2), f(3), u)
        else
            call near_switch_pair(d(1), d(2), d(3), alpha, g, r1, r2, f(1), f(2), &
                                  f(3), u)
        end if
    end subroutine

    pure subroutine near_pair_jvp(split, d, alpha, g, dd, df)
        type(near_split_t), intent(in) :: split
        real(dp), intent(in) :: d(3), alpha, g, dd(3)
        real(dp), intent(out) :: df(3)
        real(dp) :: r, r1, r2, du
        df = 0
        call near_radii(split, g, r1, r2)
        r = sqrt(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
        if (r >= r2) return
        if (r <= r1) then
            call ewald_real_pair_jvp(d(1), d(2), d(3), alpha, g, dd(1), dd(2), dd(3), &
                                     0.0_dp, 0.0_dp, df(1), df(2), df(3), du)
        else
            call near_switch_pair_jvp(d(1), d(2), d(3), alpha, g, r1, r2, dd(1), &
                dd(2), dd(3), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, df(1), df(2), df(3), du)
        end if
    end subroutine

    !> Outer step dt, reference relative speed and outer/inner step factors.
    subroutine near_setup(near, n, dt, speed, eta_outer, eta_inner, ok)
        type(near_t), intent(out) :: near
        integer, intent(in) :: n
        real(dp), intent(in) :: dt, speed, eta_outer, eta_inner
        logical, intent(out) :: ok
        real(dp) :: tau
        ok = n >= 1 .and. dt > 0 .and. speed > 0 .and. eta_outer > 0 .and. &
             eta_inner > 0
        if (.not. ok) return
        tau = dt/eta_outer
        near%split = near_split_t(.true., speed*tau, tau*tau)
        near%dt = dt; near%eta = eta_inner; near%n = n
        allocate (near%pair(2, max(16, n)), near%next(n), near%parent(n))
        allocate (near%head(1), near%order(n), near%first(n + 1))
        allocate (near%pair_order(max(16, n)), near%pair_first(n + 1), near%cursor(n + 1))
        allocate (near%x0(3*n), near%v0(3*n), near%dx0(3*n), near%dv0(3*n))
        allocate (near%dev(n), near%f(3*n), near%df(3*n), near%active(n))
    end subroutine

    !> Flow of T + U_near over dt for positions x(3n) and velocities v(3n);
    !> optional tangents (dx, dv) follow the same substeps. xs/vs are caller
    !> scratch (3n) holding the start state for verification and repeats.
    subroutine near_flow(near, x, v, xs, vs, mass, charge, fixed, coulomb, alpha, &
                         length, dt, status, dx, dv)
        type(near_t), intent(inout) :: near
        real(dp), intent(inout) :: x(:), v(:), xs(:), vs(:)
        real(dp), intent(in) :: mass(:), charge(:), coulomb, alpha, length, dt
        logical, intent(in) :: fixed(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(inout), optional :: dx(:), dv(:)
        integer :: n, i, round, cells, violator(2)
        real(dp) :: rs, vmax, gmax, r1, r2
        logical :: tangent, ok
        n = near%n; tangent = present(dx)
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid near-flow input')
        if (size(x) /= 3*n .or. size(v) /= 3*n) return
        xs = x; vs = v
        if (tangent) then
            near%dx0 = dx; near%dv0 = dv
        end if
        vmax = 0; gmax = 0
        do i = 1, n
            if (fixed(i)) cycle
            vmax = max(vmax, norm2(v(3*i - 2:3*i)))
            gmax = max(gmax, abs(coulomb)*charge(i)**2)
        end do
        call near_radii(near%split, gmax, r1, r2)
        rs = 1.5_dp*r2 + 2*vmax*dt
        cells = int(length/rs)
        if (cells < 3) cells = 0
        call build_cells(cells)
        near%pairs = 0
        call search()
        do round = 1, 4
            near%rounds = round
            call integrate(ok)
            if (.not. ok) return
            call verify(violator)
            if (violator(1) == 0) then
                call status_set(status, FORTNUM_OK, '')
                return
            end if
            call add_pair(violator(1), violator(2))
            x = xs; v = vs
            if (tangent) then
                dx = near%dx0; dv = near%dv0
            end if
        end do
        call status_set(status, FORTNUM_DOMAIN_ERROR, &
            'near-encounter cluster decomposition not verified after 4 rounds')
    contains
        subroutine build_cells(nc)
            integer, intent(in) :: nc
            integer :: c
            if (nc == 0) return
            if (size(near%head) /= nc**3) then
                deallocate (near%head); allocate (near%head(nc**3))
            end if
            near%head = 0
            do i = 1, n
                if (fixed(i)) cycle
                c = cell_of(xs(3*i - 2:3*i), nc)
                near%next(i) = near%head(c); near%head(c) = i
            end do
        end subroutine

        integer function cell_of(p, nc) result(c)
            real(dp), intent(in) :: p(3)
            integer, intent(in) :: nc
            integer :: q(3), a
            do a = 1, 3
                q(a) = min(nc - 1, int((p(a)/length - floor(p(a)/length))*nc))
            end do
            c = 1 + q(1) + nc*(q(2) + nc*q(3))
        end function

        !> Visit every mobile pair in neighbouring cells (each unordered pair
        !> once) and every (mobile, fixed) pair: mode 1 candidate search,
        !> mode 2 verification (stops at the first violation).
        subroutine for_neighbours(mode)
            integer, intent(in) :: mode
            integer :: a, b, c, oa, ob, oc, c1, c2, p, q, s
            if (cells == 0) then
                do p = 1, n
                    if (fixed(p)) cycle
                    do q = p + 1, n
                        if (.not. fixed(q)) call visit(mode, p, q)
                    end do
                end do
            else
                do c = 0, cells - 1
                    do b = 0, cells - 1
                        do a = 0, cells - 1
                            c1 = 1 + a + cells*(b + cells*c)
                            do oc = -1, 1
                                do ob = -1, 1
                                    do oa = -1, 1
                                        c2 = 1 + modulo(a + oa, cells) + cells*( &
                                             modulo(b + ob, cells) + cells*modulo(c + oc, cells))
                                        if (c2 < c1) cycle
                                        p = near%head(c1)
                                        do while (p > 0)
                                            q = near%head(c2)
                                            if (c2 == c1) q = near%next(p)
                                            do while (q > 0)
                                                call visit(mode, p, q)
                                                if (violator(1) > 0) return
                                                q = near%next(q)
                                            end do
                                            p = near%next(p)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end if
            do s = 1, n
                if (.not. fixed(s)) cycle
                do p = 1, n
                    if (.not. fixed(p)) call visit(mode, p, s)
                end do
            end do
        end subroutine

        subroutine visit(mode, i, j)
            integer, intent(in) :: mode, i, j
            real(dp) :: g, a1, a2, d0(3), w(3), dmin, devsum
            if (violator(1) > 0) return
            g = coulomb*charge(i)*charge(j)
            call near_radii(near%split, g, a1, a2)
            d0 = xs(3*i - 2:3*i) - xs(3*j - 2:3*j)
            d0 = d0 - length*anint(d0/length)
            w = vs(3*i - 2:3*i) - vs(3*j - 2:3*j)
            dmin = straight_min(d0, w, dt)
            if (mode == 1) then
                if (dmin < 1.5_dp*a2) call add_pair(i, j)
                return
            end if
            devsum = near%dev(i) + near%dev(j)
            if (devsum <= 0.5_dp*a2 .or. dmin - devsum >= a2) return
            if (is_candidate(i, j)) return
            violator = [i, j]
        end subroutine

        subroutine search()
            violator = 0
            call for_neighbours(1)
        end subroutine

        subroutine add_pair(i, j)
            integer, intent(in) :: i, j
            integer, allocatable :: grown(:, :)
            if (near%pairs == size(near%pair, 2)) then
                allocate (grown(2, 2*size(near%pair, 2)))
                grown(:, :near%pairs) = near%pair(:, :near%pairs)
                call move_alloc(grown, near%pair)
                deallocate (near%pair_order); allocate (near%pair_order(size(near%pair, 2)))
            end if
            near%pairs = near%pairs + 1
            ! Mobile endpoint first; the second may be fixed.
            if (fixed(i)) then
                near%pair(:, near%pairs) = [j, i]
            else
                near%pair(:, near%pairs) = [i, j]
            end if
        end subroutine

        logical function is_candidate(i, j) result(found)
            integer, intent(in) :: i, j
            integer :: r, k, p
            found = .false.
            r = 0
            if (near%active(i)) r = find(i)
            if (r == 0 .and. near%active(j)) r = find(j)
            if (r == 0) return
            do k = near%pair_first(r), near%pair_first(r + 1) - 1
                p = near%pair_order(k)
                if ((near%pair(1, p) == i .and. near%pair(2, p) == j) .or. &
                    (near%pair(1, p) == j .and. near%pair(2, p) == i)) then
                    found = .true.; return
                end if
            end do
        end function

        integer function find(i) result(r)
            integer, intent(in) :: i
            r = i
            do while (near%parent(r) /= r)
                near%parent(r) = near%parent(near%parent(r))
                r = near%parent(r)
            end do
        end function

        subroutine verify(bad)
            integer, intent(out) :: bad(2)
            integer :: p, q
            real(dp) :: big
            violator = 0
            call for_neighbours(2)
            ! Particles that moved far from their line: check all partners
            ! (pairs beyond the search radius are covered only for dev <= r2/4).
            big = 0.25_dp*r2
            do p = 1, n
                if (violator(1) > 0) exit
                if (near%dev(p) <= big) cycle
                do q = 1, n
                    if (q /= p .and. .not. (fixed(p) .and. fixed(q))) call visit(2, p, q)
                end do
            end do
            bad = violator
        end subroutine

        subroutine integrate(ok)
            logical, intent(out) :: ok
            integer :: p, a, b, r, k
            ok = .true.
            do i = 1, n
                near%parent(i) = i
            end do
            near%active = .false.; near%dev = 0
            do p = 1, near%pairs
                a = near%pair(1, p); b = near%pair(2, p)
                near%active(a) = .true.
                if (.not. fixed(b)) then
                    near%active(b) = .true.
                    if (find(a) /= find(b)) near%parent(find(a)) = find(b)
                end if
            end do
            ! Counting sort of members and pairs by cluster root.
            near%first = 0; near%pair_first = 0
            do i = 1, n
                if (near%active(i)) then
                    r = find(i); near%first(r + 1) = near%first(r + 1) + 1
                end if
            end do
            do p = 1, near%pairs
                r = find(near%pair(1, p)); near%pair_first(r + 1) = near%pair_first(r + 1) + 1
            end do
            near%first(1) = 1; near%pair_first(1) = 1
            do k = 2, n + 1
                near%first(k) = near%first(k) + near%first(k - 1)
                near%pair_first(k) = near%pair_first(k) + near%pair_first(k - 1)
            end do
            near%cursor(:n) = near%first(:n)
            do i = 1, n
                if (.not. near%active(i)) cycle
                r = find(i); near%order(near%cursor(r)) = i
                near%cursor(r) = near%cursor(r) + 1
            end do
            near%cursor(:n) = near%pair_first(:n)
            do p = 1, near%pairs
                r = find(near%pair(1, p)); near%pair_order(near%cursor(r)) = p
                near%cursor(r) = near%cursor(r) + 1
            end do
            near%clusters = 0; near%members = 0; near%deepest = 0
            do r = 1, n
                if (near%first(r + 1) == near%first(r)) cycle
                near%clusters = near%clusters + 1
                near%members = near%members + near%first(r + 1) - near%first(r)
                call cluster(r, ok)
                if (.not. ok) return
            end do
            do i = 1, n
                if (near%active(i) .or. fixed(i)) cycle
                x(3*i - 2:3*i) = x(3*i - 2:3*i) + dt*v(3*i - 2:3*i)
                if (tangent) dx(3*i - 2:3*i) = dx(3*i - 2:3*i) + dt*dv(3*i - 2:3*i)
            end do
        end subroutine

        !> M equal KDK substeps of the cluster rooted at r.
        subroutine cluster(r, ok)
            integer, intent(in) :: r
            logical, intent(out) :: ok
            integer :: k, m, s, i1, i2
            real(dp) :: h, shortest
            i1 = near%first(r); i2 = near%first(r + 1) - 1
            shortest = huge(1.0_dp)
            do k = near%pair_first(r), near%pair_first(r + 1) - 1
                shortest = min(shortest, encounter_time(near%pair_order(k)))
            end do
            ok = shortest > 0 .and. ieee_is_finite(shortest)
            if (.not. ok) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    'coincident charges in a near encounter; no softening applied')
                return
            end if
            m = max(1, ceiling(dt/(near%eta*shortest)))
            if (dt/(near%eta*shortest) > near%max_substeps) then
                ok = .false.
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    'near encounter needs more than max_substeps substeps')
                return
            end if
            near%deepest = max(near%deepest, m)
            h = dt/m
            call cluster_force(r, ok)
            if (.not. ok) return
            call kick(i1, i2, h/2)
            do s = 1, m
                do k = i1, i2
                    i = near%order(k)
                    x(3*i - 2:3*i) = x(3*i - 2:3*i) + h*v(3*i - 2:3*i)
                    if (tangent) dx(3*i - 2:3*i) = dx(3*i - 2:3*i) + h*dv(3*i - 2:3*i)
                    near%dev(i) = max(near%dev(i), norm2(x(3*i - 2:3*i) - &
                        xs(3*i - 2:3*i) - (s*h)*vs(3*i - 2:3*i)))
                end do
                call cluster_force(r, ok)
                if (.not. ok) return
                if (s < m) then
                    call kick(i1, i2, h)
                else
                    call kick(i1, i2, h/2)
                end if
            end do
        end subroutine

        subroutine kick(i1, i2, tau)
            integer, intent(in) :: i1, i2
            real(dp), intent(in) :: tau
            integer :: kk, j
            do kk = i1, i2
                j = near%order(kk)
                v(3*j - 2:3*j) = v(3*j - 2:3*j) + tau*near%f(3*j - 2:3*j)/mass(j)
                if (tangent) dv(3*j - 2:3*j) = dv(3*j - 2:3*j) + &
                    tau*near%df(3*j - 2:3*j)/mass(j)
            end do
        end subroutine

        subroutine cluster_force(r, ok)
            integer, intent(in) :: r
            logical, intent(out) :: ok
            integer :: k, p, a, b
            real(dp) :: d(3), f(3), u, dd(3), df(3), g
            do k = near%first(r), near%first(r + 1) - 1
                a = near%order(k); near%f(3*a - 2:3*a) = 0
                if (tangent) near%df(3*a - 2:3*a) = 0
            end do
            ok = .true.
            do k = near%pair_first(r), near%pair_first(r + 1) - 1
                p = near%pair_order(k); a = near%pair(1, p); b = near%pair(2, p)
                d = x(3*a - 2:3*a) - x(3*b - 2:3*b)
                d = d - length*anint(d/length)
                if (maxval(abs(d)) <= 0) then
                    ok = .false.
                    call status_set(status, FORTNUM_DOMAIN_ERROR, &
                        'coincident charges in a near encounter; no softening applied')
                    return
                end if
                g = coulomb*charge(a)*charge(b)
                call near_pair(near%split, d, alpha, g, f, u)
                near%f(3*a - 2:3*a) = near%f(3*a - 2:3*a) + f
                if (.not. fixed(b)) near%f(3*b - 2:3*b) = near%f(3*b - 2:3*b) - f
                if (tangent) then
                    dd = dx(3*a - 2:3*a)
                    if (.not. fixed(b)) dd = dd - dx(3*b - 2:3*b)
                    call near_pair_jvp(near%split, d, alpha, g, dd, df)
                    near%df(3*a - 2:3*a) = near%df(3*a - 2:3*a) + df
                    if (.not. fixed(b)) near%df(3*b - 2:3*b) = near%df(3*b - 2:3*b) - df
                end if
            end do
        end subroutine

        !> Shortest of the two-body periapsis Kepler time sqrt(mu rp^3/|g|),
        !> the periapsis crossing time rp/v_max and the switch crossing time
        !> (r2 - r1)/v_max of pair p from the start state (repulsive g: rp is
        !> the exact isolated-pair closest approach; attractive: the straight
        !> line closest approach).
        real(dp) function encounter_time(p) result(t)
            integer, intent(in) :: p
            integer :: a, b
            real(dp) :: d0(3), w(3), g, mu, e, l2, rnow, rp, vfast, a1, a2
            a = near%pair(1, p); b = near%pair(2, p)
            d0 = xs(3*a - 2:3*a) - xs(3*b - 2:3*b)
            d0 = d0 - length*anint(d0/length)
            w = vs(3*a - 2:3*a) - vs(3*b - 2:3*b)
            g = coulomb*charge(a)*charge(b)
            mu = mass(a)
            if (.not. fixed(b)) mu = mass(a)*mass(b)/(mass(a) + mass(b))
            rnow = norm2(d0)
            if (g > 0) then
                e = mu*sum(w**2)/2 + g/rnow
                l2 = mu**2*sum(cross(d0, w)**2)
                rp = (g + sqrt(g*g + 2*e*l2/mu))/(2*e)
                if (dot_product(d0, w) >= 0) rp = rnow
                vfast = sqrt(2*e/mu)
            else
                rp = straight_min(d0, w, dt)
                vfast = norm2(w) + sqrt(2*abs(g)/(mu*max(rp, tiny(1.0_dp))))
            end if
            call near_radii(near%split, g, a1, a2)
            t = sqrt(mu*rp**3/max(abs(g), tiny(1.0_dp)))
            if (vfast > 0) t = min(t, rp/vfast, (a2 - a1)/vfast)
        end function
    end subroutine

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function

    !> Closest distance of d0 + w t over t in [0, dt].
    pure real(dp) function straight_min(d0, w, dt) result(dmin)
        real(dp), intent(in) :: d0(3), w(3), dt
        real(dp) :: t, w2
        w2 = sum(w**2)
        t = 0
        if (w2 > 0) t = min(dt, max(0.0_dp, -dot_product(d0, w)/w2))
        dmin = norm2(d0 + t*w)
    end function
end module kin6d_near_encounter
