!> Finite-orbit geometry of the orbit-class spectral elements
!! (milca_ntv_orbit_geometry) in its thin-orbit limit. On a synthetic Boozer
!! file with B = 2.3 (1 - 0.1 sqrt(s) cos(theta)), constant G, I = 0 and no
!! magnetic shear (with shear the orbit also precesses through -q' <theta
!! psi_dot>, which the local thin field lacks) the
!! orbit has theta_dot = v_par iota B/G + H_psi and zeta_dot = v_par B/G, so in
!! the frame exp(-i n zeta) it precesses at n q H_psi with
!!   H_psi = (m v^2/(2 Z e B)) (1 + xi^2) dB/dpsi.
!! The thin field (hp_field_t) with Boozer streaming, the full drift weight
!! and the precession prefactor
!!   pre = q m v^2 dB/ds(cos)/(2 Z e B_ref psi'),
!! (offset 0) solves the same problem to first order in the orbit width, so
!! the NTV kernels of the two geometries agree when the orbit is thin.
!!
!! The cases dwell_area, deposit_pairing and deposit_converges test the
!! radial deposit of the kernel (hp_orbit_bands, hp_deposit) on the same
!! file with the deuteron's own mass, where a 10 keV orbit is about 0.07
!! wide in s.

!> The orbit geometry with the drive removed on one side of xi = 0
!! (one_sided_t) or on one side of an orbit level (level_cut_t).
module ntv_one_sided_geometry
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_dk_hp, only: hp_operator_t
    use milca_ntv_orbit_geometry, only: hp_orbit_geometry_t
    implicit none
    type, extends(hp_orbit_geometry_t) :: one_sided_t
        real(dp) :: side = 1
    contains
        procedure :: point_source => one_sided_source
    end type
    type, extends(hp_orbit_geometry_t) :: level_cut_t
        real(dp) :: cut = 0, keep = 1
    contains
        procedure :: point_source => level_cut_source
    end type
contains
    subroutine one_sided_source(g, th, xi, op, src, collide, drive_r, drive_e, density, &
            phase, precess, drive, torque, level, valid)
        class(one_sided_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(hp_operator_t), intent(in) :: op
        complex(dp), intent(out) :: src
        real(dp), intent(out) :: collide, drive_r, drive_e, density, phase, precess, drive, &
            torque
        real(dp), intent(out), optional :: level
        logical, intent(out), optional :: valid
        call g%hp_orbit_geometry_t%point_source(th, xi, op, src, collide, drive_r, drive_e, &
            density, phase, precess, drive, torque, level, valid)
        if (g%side*xi < 0) src = 0
    end subroutine one_sided_source

    !> Drive only at levels above cut (keep = 1) or below it (keep = -1).
    subroutine level_cut_source(g, th, xi, op, src, collide, drive_r, drive_e, density, &
            phase, precess, drive, torque, level, valid)
        class(level_cut_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(hp_operator_t), intent(in) :: op
        complex(dp), intent(out) :: src
        real(dp), intent(out) :: collide, drive_r, drive_e, density, phase, precess, drive, &
            torque
        real(dp), intent(out), optional :: level
        logical, intent(out), optional :: valid
        real(dp) :: lv
        logical :: has
        call g%hp_orbit_geometry_t%point_source(th, xi, op, src, collide, drive_r, drive_e, &
            density, phase, precess, drive, torque, lv, has)
        if (present(level)) level = lv
        if (present(valid)) valid = has
        if (g%keep*(lv - g%cut) < 0) src = 0
    end subroutine level_cut_source
end module ntv_one_sided_geometry

program test_ntv_orbit_geometry
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_boozer_file, only: boozer_file_t, load_boozer_file, BOOZER_OK
    use milca_boozer_orbit, only: orbit_particle_t, orbit_state, orbit_point_t, orbit_point, &
        orbit_frequencies, ORBIT_OK
    use milca_dk_hp, only: hp_field_t, hp_field, hp_mesh_t, hp_default_mesh, hp_operator_t, &
        hp_ntv_kernel, hp_bands_t, hp_orbit_bands, HP_OK
    use ntv_one_sided_geometry, only: one_sided_t, level_cut_t
    use milca_ntv_orbit_geometry, only: hp_orbit_geometry_t, hp_orbit_geometry, &
        hp_orbit_perturbation
    implicit none
    real(dp), parameter :: pi = 4*atan(1.0_dp), e = 1.602176634e-19_dp
    real(dp), parameter :: md = 3.3435837768e-27_dp, s_star = 0.45_dp, b_ref = 2.3_dp
    character(*), parameter :: file = 'ntv_orbit_geometry.bc'
    type(boozer_file_t), target :: bf, pert
    type(orbit_particle_t) :: p
    type(hp_orbit_geometry_t) :: g
    type(hp_field_t) :: f
    type(hp_mesh_t) :: mesh_g, mesh_f
    type(hp_operator_t) :: op_g, op_f
    real(dp) :: b, b_s, gg, i, g_s, i_s, iota, psi_prime, v, r_g, r_f, samples(33), &
        omega_d, factor
    integer :: status, j
    integer(8) :: t0, t1, rate
    character(32) :: arg

    call write_file(file)
    call load_boozer_file(file, bf, status)
    if (status /= BOOZER_OK) error stop 'cannot read the test file'
    factor = 1e-3_dp
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg)
        select case (trim(arg))
        case ('dwell_area')
            if (.not. dwell_area()) error stop 1
            stop
        case ('deposit_pairing')
            if (.not. deposit_pairing()) error stop 1
            stop
        case ('deposit_converges')
            if (.not. deposit_converges()) error stop 1
            stop
        case default
            read (arg, *) factor
        end select
    end if
    p%mass = md*factor
    p%charge = e
    v = sqrt(2*10e3_dp*e/md)
    p%energy = p%mass*v**2/2
    p%s_ref = s_star
    call hp_orbit_geometry(bf, p, s_star, b_ref, g, status)
    if (status /= 0) error stop 'orbit geometry failed'
    call orbit_state(bf, s_star, 0.0_dp, b, b_s, gg, i, g_s, i_s, iota, psi_prime, status)
    do j = 1, 33
        call orbit_state(bf, s_star, 2*pi*(j - 1)/33, b, b_s, gg, i, g_s, i_s, iota, &
            psi_prime, status)
        samples(j) = b/b_ref
    end do
    call hp_field(samples, iota*b_ref/gg, f, status)
    if (status /= HP_OK) error stop 'hp_field failed'

    ! dB/ds = -0.23/(2 sqrt(s)) cos(theta) for this file.
    op_f%pre = (1/iota)*p%mass*v**2*(-0.23_dp/(2*sqrt(s_star)))/(2*e*b_ref*psi_prime)
    omega_d = abs(op_f%pre)/v
    op_f%k = iota*b_ref/gg
    op_f%k_power = 1
    op_f%q = 1/iota
    op_f%n = 3
    op_f%m = [3]
    op_f%c = [(1e-3_dp, 0.0_dp)]
    op_f%v = v
    op_f%wsign = 1
    op_f%offset = 0
    op_f%nu = 0.1_dp*op_f%n*omega_d
    op_g = op_f
    op_g%pre = 0

    call hp_default_mesh(f, 4, 16, 6, 1e-3_dp, 0.3_dp, mesh_f)
    call hp_ntv_kernel(f, mesh_f, op_f, r_f, status)
    if (status /= HP_OK) error stop 'thin kernel failed'
    call system_clock(t0, rate)
    call hp_default_mesh(g, 4, 16, 6, 1e-3_dp, 0.3_dp, mesh_g)
    call hp_ntv_kernel(g, mesh_g, op_g, r_g, status)
    call system_clock(t1)
    if (status /= HP_OK) error stop 'orbit kernel failed'
    print '(a,2es20.11,a,es9.2,a,f8.2)', 'R thin field, orbit geometry ', r_f, r_g, ' rel ', &
        abs(r_g/r_f - 1), ' orbit s ', real(t1 - t0, dp)/rate
    if (abs(r_g/r_f - 1) > 1e-3_dp) error stop 1

    ! The drive taken along the orbit from a perturbation file whose m = 3
    ! amplitude at s* equals the operator's: the same kernel for thin orbits.
    call write_pert('ntv_orbit_pert.bc')
    call load_boozer_file('ntv_orbit_pert.bc', pert, status)
    if (status /= BOOZER_OK) error stop 'cannot read the perturbation file'
    call hp_orbit_perturbation(g, pert, 3, b_ref)
    op_g%c = [cmplx(1e-3_dp*b_ref*s_star**1.5_dp, 0.0_dp, dp)/b_ref]
    op_g%m = [-3]
    call hp_ntv_kernel(g, mesh_g, op_g, r_f, status)
    if (status /= HP_OK) error stop 'orbit kernel with perturbation file failed'
    g%pert => null()
    call hp_ntv_kernel(g, mesh_g, op_g, r_g, status)
    print '(a,2es20.11,a,es9.2)', 'R drive along orbit, drive at s* ', r_f, r_g, ' rel ', &
        abs(r_f/r_g - 1)
    if (abs(r_f/r_g - 1) > 1e-4_dp) error stop 1

contains

    !> Geometry and operator of the deuteron at scale times 10 keV on s_star,
    !> with the collision frequency nu_rel times the bounce frequency of a
    !> trapped orbit halfway between the separatrix and the bottom of the well.
    subroutine finite_width(scale, nu_rel, gf, op)
        real(dp), intent(in) :: scale, nu_rel
        type(hp_orbit_geometry_t), intent(out) :: gf
        type(hp_operator_t), intent(out) :: op
        type(orbit_particle_t) :: q
        real(dp) :: omega_b, omega_zeta
        logical :: trapped
        integer :: st
        q%mass = md
        q%charge = e
        q%energy = 10e3_dp*e*scale
        q%s_ref = s_star
        call hp_orbit_geometry(bf, q, s_star, b_ref, gf, st)
        if (st /= 0) error stop 'orbit geometry failed'
        call orbit_state(bf, s_star, 0.0_dp, b, b_s, gg, i, g_s, i_s, iota, psi_prime, st)
        op%k = iota*b_ref/gg
        op%k_power = 1
        op%q = 1/iota
        op%n = 3
        op%m = [3]
        op%c = [(1e-3_dp, 0.0_dp)]
        op%v = gf%v_star
        op%wsign = 1
        q%mu = 0.5_dp*(gf%lambda_c + b_ref/b)*gf%e_star/b_ref
        call orbit_frequencies(bf, q, s_star, 1.0_dp, 64, omega_b, omega_zeta, trapped, st)
        if (st /= ORBIT_OK .or. .not. trapped) error stop 'no trapped reference orbit'
        op%nu = nu_rel*omega_b/gf%v_star
    end subroutine finite_width

    !> Hat function of a trapped level node, constant beyond the end nodes.
    pure real(dp) function hat(l, lo, peak, hi) result(h)
        real(dp), intent(in) :: l, lo, peak, hi
        h = 1
        if (l < peak .and. lo < peak) h = max(0.0_dp, (l - lo)/(peak - lo))
        if (l > peak .and. hi > peak) h = max(0.0_dp, (hi - l)/(hi - peak))
    end function hat

    !> The time tables of hp_orbit_bands are phase-space measures (Liouville):
    !> time(k, b) = int dtheta dxi |measure| w_b(mu_hat) [s in box k], here by
    !> the midpoint rule in (theta, xi), which uses neither orbit integrals
    !> nor turning points: a passing band between its two level curves, the
    !> trapped bands with their hat functions between the separatrix curves.
    !> Errors are relative to the band's measure and split in two. The excess
    !> is the error of the orbit model itself: the level integral of the
    !> geometry's dwell with 12 Gauss levels per interval and every column
    !> bisected three times, against the phase-space measure. It stays below
    !> 1 %, and 5 % for the node on the separatrix, whose period diverges,
    !> from 10 keV (orbit about 0.07 wide in s) down to 1/16 of it, with the
    !> boxes scaled with the width. Without the caps beyond the turning points
    !> it is first order in the width: 11 % in the deepest node and 6 % on
    !> the separatrix at 10 keV. The floor is the quadrature of the tables
    !> (two and three levels per band) against the fine integral: up to 2 %
    !> where a box edge or the separatrix cuts a band.
    logical function dwell_area() result(ok)
        use fortnum_quadrature, only: gauss_legendre
        integer, parameter :: n_th = 3200, n_xi = 192, n_g = 12, n_run = 3
        type(hp_orbit_geometry_t) :: gf
        type(hp_operator_t) :: op
        type(hp_mesh_t) :: mesh
        type(hp_bands_t) :: bands
        type(orbit_point_t) :: pt
        real(dp) :: edges(6), tau(5), cut(-1:1), th, xi, xa, xb, d, l, w, xg(n_g), wg(n_g), &
            a, b, scale, floor, excess, tol
        real(dp), allocatable :: panels(:), ref(:, :), fine(:, :)
        integer :: nb, bd, jt, jx, part, k, st, n_pass, run, g, side, b1, b2, nx
        ok = .true.
        call gauss_legendre(n_g, xg, wg)
        do run = 1, n_run
            scale = 0.25_dp**(run - 1)
            edges = [0.1_dp, s_star + sqrt(scale)*[-0.03_dp, -0.01_dp, 0.01_dp, 0.03_dp], 0.9_dp]
            call finite_width(scale, 1e-2_dp, gf, op)
            call hp_default_mesh(gf, 4, 8, 3, 1e-2_dp, 0.3_dp, mesh)
            call hp_orbit_bands(gf, mesh, edges, bands, st)
            if (st /= HP_OK) error stop 'hp_orbit_bands failed'
            panels = halved(halved(halved(mesh%theta)))
            nb = size(bands%kind)
            n_pass = size(mesh%lam) - 1
            ok = ok .and. bands%unplaced == 0
            if (allocated(ref)) deallocate (ref, fine)
            allocate (ref(5, nb), fine(5, nb))
            fine = 0
            do bd = 1, nb
                do side = 1, 2
                    a = merge(bands%lo(bd), bands%level(bd), side == 1)
                    b = merge(bands%level(bd), bands%hi(bd), side == 1)
                    if (b <= a) cycle
                    do g = 1, n_g
                        call gf%dwell(bands%kind(bd), a + (b - a)*(1 + xg(g))/2, panels, edges, &
                            tau, st)
                        w = (b - a)/2*wg(g)
                        if (bands%kind(bd) == 0) w = w*merge(1 + xg(g), 1 - xg(g), side == 1)/2
                        fine(:, bd) = fine(:, bd) + w*tau
                    end do
                end do
            end do
            ref = 0
            !$omp parallel do default(shared) reduction(+:ref) &
            !$omp private(jt, th, bd, b1, b2, nx, cut, part, xa, xb, jx, xi, l, pt, st, w, k, d)
            do jt = 1, n_th
                th = gf%theta_x - 2*pi*(jt - 0.5_dp)/n_th
                do bd = 1, nb
                    if (bands%kind(bd) /= 0) then
                        call gf%level(bands%kind(bd), bands%lo(bd), th, cut(-1), d)
                        call gf%level(bands%kind(bd), bands%hi(bd), th, cut(1), d)
                        cut(0) = (cut(-1) + cut(1))/2
                        b1 = bd
                        b2 = bd
                        nx = n_xi/8
                    else if (bd == n_pass + 1) then
                        call gf%level(0, -1.0_dp, th, cut(-1), d)
                        call gf%level(0, 1.0_dp, th, cut(1), d)
                        cut(0) = 0
                        b1 = n_pass + 1
                        b2 = nb - n_pass
                        nx = n_xi
                    else
                        cycle
                    end if
                    do part = 0, 1
                        xa = cut(part - 1)
                        xb = cut(part)
                        do jx = 1, nx
                            xi = xa + (xb - xa)*(jx - 0.5_dp)/nx
                            l = gf%mu_hat(th, xi, pt, st)
                            if (st /= ORBIT_OK) cycle
                            w = abs(gf%measure(th, xi, op)*(xb - xa))/nx*2*pi/n_th
                            k = count(edges <= pt%s)
                            if (k < 1 .or. k >= size(edges)) cycle
                            if (b1 == b2) then
                                ref(k, b1) = ref(k, b1) + w
                            else
                                ref(k, b1:b2) = ref(k, b1:b2) + w*[(hat(l, bands%lo(g), &
                                    bands%level(g), bands%hi(g)), g=b1, b2)]
                            end if
                        end do
                    end do
                end do
            end do
            !$omp end parallel do
            do bd = 1, nb
                floor = maxval(abs(bands%time(:, bd) - fine(:, bd)))/sum(ref(:, bd))
                excess = maxval(abs(fine(:, bd) - ref(:, bd)))/sum(ref(:, bd))
                tol = merge(0.05_dp, 0.01_dp, bd == n_pass + 1)
                print '(a,es8.1,a,i3,i3,2f9.5,a,5f7.4,a,es9.2,a,es9.2,a,f5.2,l2)', 'energy', &
                    scale, ' band', bd, bands%kind(bd), bands%lo(bd), bands%hi(bd), ' share', &
                    bands%time(:, bd)/sum(ref(:, bd)), ' floor', floor, ' excess', excess, &
                    ' tol', tol, excess < tol .and. floor < 0.025_dp
                ok = ok .and. excess < tol .and. floor < 0.025_dp
            end do
        end do
    end function dwell_area

    !> Mean and variance of s over one orbit of level and sign sg (0: both
    !> legs of a trapped orbit), time weighted, by the midpoint rule on
    !> orbit_point. The turning points of the test field are analytic:
    !> B(s*, theta) = 2.3 - 0.23 sqrt(s*) cos(theta) = b_ref/level.
    subroutine orbit_moments(gf, level, sg, mean, var)
        type(hp_orbit_geometry_t), intent(in) :: gf
        real(dp), intent(in) :: level
        integer, intent(in) :: sg
        real(dp), intent(out) :: mean, var
        integer, parameter :: n = 4000
        type(orbit_particle_t) :: q
        type(orbit_point_t) :: pt
        real(dp) :: th_t, u, th, w, m0, m1, m2
        integer :: jt, leg, st
        q = gf%p
        q%mu = level*gf%e_star/b_ref
        m0 = 0; m1 = 0; m2 = 0
        if (sg /= 0) then
            th_t = pi
        else
            th_t = acos((2.3_dp - b_ref/level)/(0.23_dp*sqrt(s_star)))
        end if
        do leg = 1, merge(2, 1, sg == 0)
            do jt = 1, n
                u = pi*((jt - 0.5_dp)/n - 0.5_dp)
                th = th_t*sin(u)
                call orbit_point(bf, q, s_star, th, merge(3.0_dp - 2*leg, real(sg, dp), sg == 0), &
                    pt, st)
                if (st /= ORBIT_OK) error stop 'orbit_point failed'
                w = cos(u)/abs(pt%theta_dot)
                m0 = m0 + w
                m1 = m1 + w*(pt%s - s_star)
                m2 = m2 + w*(pt%s - s_star)**2
            end do
        end do
        mean = s_star + m1/m0
        var = m2/m0 - (m1/m0)**2
    end subroutine orbit_moments

    !> Bands, time tables and torques belong together. (i) Per band, the
    !> mean s of its time table against the orbit at the band's level:
    !> passing rows lie on their own side of s*, 0.02 to 0.13 away, and
    !> neighbouring rows differ by more than the tolerance; trapped orbits
    !> are centred on s*, so their rms width is compared as well, loosely,
    !> since the reference has the legs without the caps at s* (not for the
    !> nodes on the separatrix and at the bottom of the well, which have no
    !> orbit, nor the one next to the bottom, whose hat reaches zero width).
    !> (ii) The first radial moment of the deposit is that of the bands.
    !> (iii) A drive on one side of xi = 0 leaves the passing rows of the
    !> other side without torque, and the passing deposit on one side of s*.
    !> (iv) A drive on one side of a trapped level node leaves the trapped
    !> bands wholly on the other side without torque, exactly, and the
    !> others with torque: a Gauss point feeds the two nodes around its
    !> own level.
    logical function deposit_pairing() result(ok)
        integer, parameter :: n_box = 64
        type(hp_orbit_geometry_t) :: gf
        type(one_sided_t) :: one
        type(level_cut_t) :: cut
        type(hp_operator_t) :: op
        type(hp_mesh_t) :: mesh
        type(hp_bands_t) :: bands
        type(orbit_point_t) :: pt
        type(orbit_particle_t) :: q
        real(dp) :: edges(n_box + 1), mid(n_box), dep(n_box), dep2(2), pass(2), r, mean, rms, &
            m_ref, v_ref, total, scale, moment, e_mean, e_rms
        real(dp), allocatable :: m_band(:)
        real(dp), allocatable :: g_far(:), g_near(:)
        integer :: st, bd, nb, n_pass, k, side, box
        logical :: good
        logical, allocatable :: far(:), trap(:)
        call finite_width(1.0_dp, 1e-2_dp, gf, op)
        call hp_default_mesh(gf, 4, 16, 6, 1e-3_dp, 0.3_dp, mesh)
        edges = [(0.25_dp + 0.4_dp*k/n_box, k=0, n_box)]
        mid = (edges(:n_box) + edges(2:))/2
        call hp_ntv_kernel(gf, mesh, op, r, st, edges=edges, dep=dep, bands=bands)
        if (st /= HP_OK) error stop 'orbit kernel failed'
        nb = size(bands%kind)
        n_pass = size(mesh%lam) - 1
        allocate (m_band(nb))
        e_mean = 0
        e_rms = 0
        do bd = 1, nb
            total = sum(bands%time(:, bd))
            mean = sum(bands%time(:, bd)*mid)/total
            m_band(bd) = mean
            rms = sqrt(sum(bands%time(:, bd)*(mid - mean)**2)/total)
            if (bd == n_pass + 1 .or. bd == nb - n_pass) cycle
            call orbit_moments(gf, bands%level(bd), bands%kind(bd), m_ref, v_ref)
            e_mean = max(e_mean, abs(mean - m_ref))
            if (bands%kind(bd) == 0 .and. bd < nb - n_pass - 1) &
                e_rms = max(e_rms, abs(rms/sqrt(v_ref) - 1))
            print '(a,i3,i3,f9.5,a,2f9.5,a,2es10.3)', 'band', bd, bands%kind(bd), &
                bands%level(bd), ' mean s', mean, m_ref, ' rms', rms, sqrt(v_ref)
        end do
        ok = e_mean < 2.5e-3_dp .and. e_rms < 0.12_dp
        print '(a,es9.2,a,es9.2,l2)', 'largest error of mean s ', e_mean, ', of trapped rms ', &
            e_rms, ok

        scale = sum(abs(bands%torque))
        moment = sum(bands%torque*m_band) + bands%stray*mid(bands%home)
        good = abs(sum(dep*mid) - moment) < 1e-10_dp*scale .and. abs(sum(dep) - r) < 1e-10_dp*scale
        print '(a,2es22.14,a,es10.3,l2)', 'first moment of deposit, of bands ', sum(dep*mid), &
            moment, ' sum |G| ', scale, good
        ok = ok .and. good

        cut%hp_orbit_geometry_t = gf
        cut%cut = bands%level((nb + 1)/2)
        trap = bands%kind == 0
        do side = 1, -1, -2
            cut%keep = side
            call hp_ntv_kernel(cut, mesh, op, r, st, edges=edges, dep=dep, bands=bands)
            if (st /= HP_OK) error stop 'level-cut kernel failed'
            far = trap .and. merge(bands%hi <= cut%cut, bands%lo >= cut%cut, side == 1)
            g_far = pack(bands%torque, far)
            g_near = pack(bands%torque, trap .and. .not. far)
            scale = sum(abs(bands%torque))
            good = size(g_far) > 1 .and. size(g_near) > 1 .and. all(g_far == 0) &
                .and. all(abs(g_near) > 1e-3_dp*scale/size(g_near)) &
                .and. sum(abs(g_near)) > 1e-3_dp*scale
            print '(a,i3,a,f9.5,a,i3,es10.3,a,i3,2es10.3,l2)', 'drive on side', side, &
                ' of level ', cut%cut, ' trapped bands beyond: n, max |G|', size(g_far), &
                maxval(abs(g_far)), ' others: n, min and sum |G|/sum |G| ', size(g_near), &
                minval(abs(g_near))/scale, sum(abs(g_near))/scale, good
            ok = ok .and. good
        end do

        ! The side of s* of a co-passing orbit, from the orbit itself.
        q = gf%p
        q%mu = 0.5_dp*gf%lambda_c*gf%e_star/b_ref
        call orbit_point(bf, q, s_star, 0.0_dp, 1.0_dp, pt, st)
        if (st /= ORBIT_OK) error stop 'orbit_point failed'
        one%hp_orbit_geometry_t = gf
        do side = 1, -1, -2
            one%side = side
            box = merge(2, 1, (pt%s > s_star) .eqv. (side == 1))
            call hp_ntv_kernel(one, mesh, op, r, st, edges=[0.1_dp, s_star, 0.9_dp], dep=dep2, &
                bands=bands)
            if (st /= HP_OK) error stop 'one-sided kernel failed'
            pass = 0
            do bd = 1, nb
                if (bands%kind(bd) /= 0) pass = pass + bands%torque(bd)*bands%time(:, bd) &
                    /sum(bands%time(:, bd))
            end do
            good = all(pack(bands%torque, bands%kind == -side) == 0) &
                .and. pass(3 - box) == 0 .and. abs(pass(box)) > 1e-3_dp*sum(abs(bands%torque)) &
                .and. abs(sum(dep2) - r) < 1e-10_dp*sum(abs(bands%torque))
            print '(a,i3,a,2es11.3,a,2es11.3,l2)', 'drive on side', side, ' passing deposit ', &
                pass, ' deposit ', dep2, good
            ok = ok .and. good
        end do
    end function deposit_pairing

    !> Every interval of the mesh bisected.
    subroutine double_mesh(mesh)
        type(hp_mesh_t), intent(inout) :: mesh
        mesh%theta = halved(mesh%theta)
        mesh%lam = halved(mesh%lam)
        mesh%zeta = halved(mesh%zeta)
    end subroutine double_mesh

    pure function halved(x) result(y)
        real(dp), intent(in) :: x(:)
        real(dp) :: y(2*size(x) - 1)
        y(1::2) = x
        y(2::2) = (x(:size(x) - 1) + x(2:))/2
    end function halved

    !> At a collision frequency of 1e-3 of the bounce frequency the deposit
    !> over fixed boxes converges under bisection of all rows and columns:
    !> the largest change of a box, on the scale of the band torques, falls
    !> by more than 3 per bisection (measured at degree 2: 3.9e-2, 7.4e-3,
    !> 1.5e-3). At 1e-4 of the energy the orbit is a hundredth as wide and
    !> the whole kernel stays in the box of s*.
    logical function deposit_converges() result(ok)
        integer, parameter :: n_box = 16
        type(hp_orbit_geometry_t) :: gf
        type(hp_operator_t) :: op
        type(hp_mesh_t) :: mesh
        type(hp_bands_t) :: bands
        real(dp) :: edges(n_box + 1), dep(n_box, 0:3), r, diff(3), scale
        integer :: st, k, lvl, n, home
        edges = [(0.3_dp + 0.02_dp*k, k=0, n_box)]
        call finite_width(1.0_dp, 1e-3_dp, gf, op)
        call hp_default_mesh(gf, 2, 8, 3, 1e-2_dp, 0.3_dp, mesh)
        do lvl = 0, 3
            call hp_ntv_kernel(gf, mesh, op, r, st, n, edges=edges, dep=dep(:, lvl), bands=bands)
            if (st /= HP_OK) error stop 'orbit kernel failed'
            scale = sum(abs(bands%torque))
            if (lvl > 0) diff(lvl) = maxval(abs(dep(:, lvl) - dep(:, lvl - 1)))/scale
            print '(a,i2,a,i7,a,es16.8,a,es10.3,a,es10.3)', 'bisections', lvl, ' unknowns', n, &
                ' R', r, ' sum |G|', scale, ' change', merge(diff(max(lvl, 1)), 0.0_dp, lvl > 0)
            call double_mesh(mesh)
        end do
        ok = diff(1) > 3*diff(2) .and. diff(2) > 3*diff(3)

        call finite_width(1e-4_dp, 1e-3_dp, gf, op)
        call hp_default_mesh(gf, 2, 8, 3, 1e-2_dp, 0.3_dp, mesh)
        call hp_ntv_kernel(gf, mesh, op, r, st, edges=edges, dep=dep(:, 0), bands=bands)
        if (st /= HP_OK) error stop 'thin orbit kernel failed'
        home = count(edges <= s_star)
        scale = sum(abs(bands%torque))
        print '(a,2es16.8,a,es9.2)', 'thin orbit: R, deposit at s* ', r, dep(home, 0), &
            ' elsewhere ', sum(abs(dep(:, 0))) - abs(dep(home, 0))
        ok = ok .and. abs(dep(home, 0) - r) < 1e-10_dp*scale .and. bands%home == home &
            .and. count(dep(:, 0) /= 0) == 1 .and. abs(r) > 0
    end function deposit_converges

    !> One harmonic m = 3, n = 3 with amplitude 1e-3 b0 s^1.5.
    subroutine write_pert(name)
        character(*), intent(in) :: name
        integer :: unit, k
        real(dp) :: s
        open (newunit=unit, file=name, status='replace', action='write')
        write (unit, '(a)') 'CC synthetic test file'
        write (unit, '(a)') ' m0b n0b nsurf nper flux/[Tm^2] a/[m] R/[m]'
        write (unit, '(a)') ' 3 3 19 1 -1.3 0.5 1.7'
        do k = 1, 19
            s = 0.05_dp*k
            write (unit, '(a)') ' s iota curr_pol/nper curr_tor pprime sqrt g(0,0)'
            write (unit, '(a)') ' [A] [A] dp/ds,[Pa] (dV/ds)/nper'
            write (unit, '(6es24.16)') s, 0.765_dp, -1.2e7_dp, 0.0_dp, 0.0_dp, 0.0_dp
            write (unit, '(a)') ' m n r/[m] z/[m] (phib-phi)*nper/twopi bmn/[T]'
            write (unit, '(2i5,4es24.16)') 3, 3, 0.0_dp, 0.0_dp, 0.0_dp, 1e-3_dp*b_ref*s**1.5_dp
        end do
        close (unit)
    end subroutine write_pert

    subroutine write_file(name)
        character(*), intent(in) :: name
        integer :: unit, k
        real(dp) :: s
        open (newunit=unit, file=name, status='replace', action='write')
        write (unit, '(a)') 'CC synthetic test file'
        write (unit, '(a)') ' m0b n0b nsurf nper flux/[Tm^2] a/[m] R/[m]'
        write (unit, '(a)') ' 1 0 19 1 -1.3 0.5 1.7'
        do k = 1, 19
            s = 0.05_dp*k
            write (unit, '(a)') ' s iota curr_pol/nper curr_tor pprime sqrt g(0,0)'
            write (unit, '(a)') ' [A] [A] dp/ds,[Pa] (dV/ds)/nper'
            write (unit, '(6es24.16)') s, 0.765_dp, -1.2e7_dp, 0.0_dp, 0.0_dp, 0.0_dp
            write (unit, '(a)') ' m n r/[m] z/[m] (phib-phi)*nper/twopi bmn/[T]'
            write (unit, '(2i5,4es24.16)') 0, 0, 1.7_dp, 0.0_dp, 0.0_dp, 2.3_dp
            write (unit, '(2i5,4es24.16)') 1, 0, 0.5_dp*sqrt(s), 0.5_dp*sqrt(s), 0.0_dp, &
                -0.23_dp*sqrt(s)
        end do
        close (unit)
    end subroutine write_file
end program test_ntv_orbit_geometry
