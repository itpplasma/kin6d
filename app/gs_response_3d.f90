program gs_response_3d
    !! L4: per-n linear response of an axisymmetric G-EQDSK equilibrium to 3D
    !! coils.  The coil field (L3, closed-form Biot-Savart) drives
    !!  - model = 'shield': resonant sheet currents on the rational surfaces
    !!    that cancel the total resonant flux (ideal shielding of a vacuum
    !!    plasma, no volume current);
    !!  - model = 'ideal': in addition the volume current of the linear ideal
    !!    perturbed equilibrium (kin6d_mhd_ideal_n, MEPHIT form) with the
    !!    field-aligned boundary sheet (boundary_sheet).
    !! closure selects the treatment of the rational surfaces in the ideal
    !! model: 'sheet' (L4: one sheet per rational ring, zero resonant flux
    !! there) or 'frobenius' (L6, kin6d_mhd_closure_n: current-free layer
    !! |psin - psin_s| < layer_delta with edge sheets fixed by vanishing
    !! large-solution coefficients; rings graded geometrically with
    !! layer_ratio away from each layer, layer_nfit fit rings per side).
    !! Fields of all plasma currents come from the per-n FEEC Ampere operator
    !! on a flux-aligned mesh of the plasma with the exact exterior (BEM).
    !!
    !! Usage: gs_response_3d input.nml   (namelist /l4/, see below)
    !!
    !! Output prefix_res.dat, one line per (phasing, rational):
    !!   q psin m area Phi_vac b_vac lambda b_shield
    !! with Phi_vac the resonant flux harmonic of the coil field (T m^2 per
    !! rad^2, peak amplitude, kin6d_vacuum_surfaces phi_mn), area =
    !! (1/2 pi) int R ds (m^2 per rad), b_vac = Phi_vac/area (T), lambda the
    !! sheet stream-function amplitude (A, peak), and b_shield = -M_ss lambda/
    !! area the resonant flux of the shielding sheet alone per area (the
    !! quantity GPEC reports as singflx).  m is the resonant m = n q of the
    !! signed q of kin6d_vacuum_surfaces (GPEC labels it -m).
    !! With closure = 'frobenius' lambda and b_shield are the equivalent
    !! sheet amplitude lambda_eq of the jump of d Phi/d psin across
    !! psin_s -+ spot/(n |dq/dpsin|) (GPEC sing_spot) and its self flux per
    !! area: the net current of the layer (sum of the edge sheets) carried
    !! to spot with the small-solution power (spot/delta)^(mu - 1/2).
    !! prefix_layer.dat gives per (phasing, rational) D_I, mu, delta, the
    !! jump, the large-solution residuals |a1/Phi| on both sides, the two
    !! edge-sheet amplitudes and the jump of the fitted outer slopes.
    !! closure = 'resistive' (L5, kin6d_mhd_layer_n) uses the same layers and
    !! edge sheets with the resistive matching conditions: continuous
    !! reconnected flux A and net layer current = Delta_layer A, the
    !! resistive-inertial(-viscous) layer response of a plasma rotating with
    !! Omega(psin) to the drive of frequency omega_drive (rad/s), from the
    !! table `profiles` (columns psin, eta [Ohm m], rho [kg/m^3], Omega
    !! [rad/s]; '#' comments), with rho scaled by 1 + 2 q^2 if inertia_q2
    !! and viscosity layer_nu (m^2/s).  prefix_resist.dat gives per
    !! (phasing, rational) Delta_layer (per psin), the layer scale L (psin),
    !! Q = eps^(-2/3), the reconnected flux A (T m^2/rad^2) and A/area (T),
    !! the layer jump of d Phi/d psin at the layer edges, the same carried to
    !! spot, and |Phi_vac|.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_point_basis_t, meridian_point_basis_init
    use kin6d_coils_3d
    use kin6d_vacuum_surfaces
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build, ring_labels_merge, &
        ring_labels_smooth
    use kin6d_mhd_ampere_n, only: ampere_n_t, ampere_n_build, ampere_n_solve, &
        ampere_n_curl
    use kin6d_mhd_response_n, only: ring_angles_t, ring_angles_set, &
        ring_points, ring_flux_harmonic, ring_flux_from_values, &
        ring_sheet_load, ring_area, resonant_surfaces_missing
    use kin6d_mhd_ideal_n, only: sample_ring_t, sample_ring_build, &
        ideal_options_t, ideal_response, sample_field
    use kin6d_mhd_perturbation_n, only: perturbation_t, perturbation_on_surface, &
        perturbation_add_spectra, perturbation_export
    use kin6d_mhd_response_n, only: sheet_closure_t, closure_shield_init, &
        closure_fluxes
    use kin6d_mhd_closure_n, only: rational_layer_t, mercier_di, layer_labels, layer_xmax, &
        layer_init, layer_closure, layer_coefficients, layer_slope, &
        layer_regular_slope, layer_jump, layer_kink, layer_large, layer_inside
    use kin6d_mhd_layer_n, only: layer_params_t, layer_scales, layer_delta_x => &
        layer_delta, layer_metrics, resistive_closure
    implicit none
    integer, parameter :: maxr = 32, maxp = 4
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: mu0_l = 4.0e-7_dp*pi
    character(len=512) :: eqdsk = '', coils = '', currents(maxp) = '', prefix = 'l4'
    character(len=512) :: nml
    integer :: n = 3, nrad = 64, n_core = 8, nth_max = 512, order = 2
    integer :: nphi = 64, nfine = 1024, angle = 2
    integer :: rat_num(maxr) = 0, rat_den(maxr) = 0
    real(dp) :: psin_rat_min = 0.05_dp, psin_rat_max = 1.0_dp, min_gap = 0.3_dp
    character(len=16) :: model = 'shield'
    integer :: np_sample = 512, gmres_maxit = 300
    real(dp) :: gmres_tol = 1.0e-9_dp
    logical :: boundary_sheet = .true.
    real(dp) :: qlow = 0.0_dp, rat_pad = 0.0_dp
    integer :: rat_grade = 1
    logical :: edge_moved_current = .false.
    logical :: rat_auto = .true.
    character(len=16) :: closure = 'sheet'
    real(dp) :: layer_delta = 1.0e-3_dp, layer_ratio = 2.0_dp, spot = 5.0e-4_dp
    integer :: layer_nfit = 3
    !! L7 export: psin of the surfaces, points per surface, spectral window
    real(dp) :: export_psin(maxr) = 0.0_dp
    integer :: export_np = 256, export_mmax = 40
    character(len=512) :: profiles = ''
    real(dp) :: omega_drive = 0.0_dp, layer_nu = 0.0_dp, eta_scale = 1.0_dp
    logical :: inertia_q2 = .false.
    logical :: zero_pressure = .false.
    namelist /l4/ eqdsk, coils, currents, prefix, n, nrad, n_core, nth_max, &
        order, nphi, nfine, angle, rat_num, rat_den, psin_rat_min, &
        psin_rat_max, min_gap, model, np_sample, gmres_maxit, gmres_tol, &
        boundary_sheet, qlow, rat_pad, rat_grade, edge_moved_current, closure, &
        layer_delta, layer_ratio, layer_nfit, spot, rat_auto, export_psin, &
        export_np, export_mmax, profiles, omega_drive, layer_nu, eta_scale, &
        inertia_q2, zero_pressure

    type(psi_grid_t) :: g
    type(coil_set_t) :: set
    type(surface_t) :: sref, surf
    type(flux_mesh_t) :: fm
    type(meridian_feec_t) :: fs
    type(meridian_point_basis_t) :: pb
    type(ampere_n_t) :: amp
    type(ring_angles_t), allocatable :: ra(:)
    real(dp), allocatable :: w(:), ul(:), u_rat(:), psin_rat(:), q_rat(:)
    real(dp), allocatable :: xg(:, :), nds(:, :), thg(:), tau(:), wq(:), area(:)
    integer, allocatable :: tg(:), m_rat(:), ring_of(:), ipiv(:)
    complex(dp), allocatable :: f(:), a(:), b2(:, :), mm(:, :), lu(:, :)
    complex(dp), allocatable :: phivac(:), lam(:), bc(:, :, :), bn(:)
    complex(dp), allocatable :: phil3(:), bs(:, :, :), bb(:, :)
    complex(dp) :: bmn1(1), phimn1(1)
    real(dp), allocatable :: xs(:, :)
    type(sample_ring_t), allocatable :: rings(:)
    type(ideal_options_t) :: opt
    type(sheet_closure_t) :: cl
    type(rational_layer_t), allocatable :: lays(:)
    real(dp), allocatable :: di_rat(:), dq_rat(:), jump_unit(:)
    complex(dp), allocatable :: b2cl(:, :), phicl(:), phitot(:), lam_cl(:)
    logical :: frob, resist
    real(dp), allocatable :: prof_tab(:, :), lscale(:), qlay(:), gpsi(:)
    complex(dp), allocatable :: dlayer(:), lcen(:)
    real(dp), allocatable :: xr(:, :)
    complex(dp), allocatable :: bcr(:, :, :), avol(:)
    real(dp) :: qref, target, root, t0, t1, res
    integer :: ulay, ures, unit, ur, st, i, j, k, nr, nrat, ip, nphase, info, kk, its, nsr
    logical :: found
    ! Discrete residual-energy indicator of the app's own unit-sheet Ampere
    ! solves (algebraic residual of the discrete system, not a discretization
    ! bound): worst value, solve count and number of solves without indicator.
    real(dp) :: ampere_indicator = 0.0_dp
    integer :: ampere_solves = 0, ampere_missing = 0

    if (command_argument_count() /= 1) error stop 'usage: gs_response_3d input.nml'
    call get_command_argument(1, nml)
    open (newunit=unit, file=trim(nml), status='old', action='read')
    read (unit, nml=l4)
    close (unit)
    nphase = count(len_trim(currents) > 0)
    if (n == 0 .or. nphase == 0) error stop 'need n /= 0 and current tables'

    call eqdsk_read(trim(eqdsk), g, st)
    if (st /= 0) error stop 'cannot read EQDSK'
    ! discriminating study: p0' = 0 (J0 and the Mercier index D_I = -1/4, mu = 1/2)
    if (zero_pressure .and. allocated(g%pcoef)) g%pcoef = 0.0_dp
    call coil_set_read(set, trim(coils), st)
    if (st /= 0) error stop 'cannot read coil table'

    ! rational surfaces of this n
    call surface_build(g, 0.5_dp, nfine, sref, st)
    if (st /= 0) error stop 'reference surface'
    qref = sref%q
    nrat = count(rat_den > 0)
    allocate (psin_rat(nrat), q_rat(nrat), m_rat(nrat))
    k = 0
    do i = 1, nrat
        target = sign(real(rat_num(i), dp)/rat_den(i), qref)
        found = find_surface_q(g, target, psin_rat_min, psin_rat_max, 256, &
            1.0e-13_dp, root, st)
        if (.not. found) then
            write (*, '(a,i0,a,i0,a)') 'rational q = ', rat_num(i), '/', &
                rat_den(i), ' not in range, skipped'
            cycle
        end if
        if (mod(n*rat_num(i), rat_den(i)) /= 0) error stop 'q not resonant for n'
        k = k + 1
        psin_rat(k) = root
        q_rat(k) = target
        m_rat(k) = nint(n*target)
        write (*, '(a,f10.6,a,f14.10,a,i0)') 'rational q = ', target, &
            ' at psin = ', root, ', m = ', m_rat(k)
    end do
    nrat = k
    resist = trim(closure) == 'resistive' .and. trim(model) == 'ideal'
    frob = (trim(closure) == 'frobenius' .and. trim(model) == 'ideal') .or. resist
    if (trim(closure) /= 'sheet' .and. trim(closure) /= 'frobenius' .and. &
        trim(closure) /= 'resistive') &
        error stop 'closure must be sheet, frobenius or resistive'
    if (resist) then
        call read_profiles(trim(profiles), prof_tab)
        write (*, '(a,i0,a)') 'resistive closure: ', size(prof_tab, 2), &
            ' profile points'
    end if

    ! Every other resonant surface m/n of the plasma (the sidebands m +- 1 of
    ! the toroidal coupling, e.g. q = 3 at the edge for n q = 2) carries the
    ! same 1/(u - u_s) pole of the ideal volume current.  Left unshielded the
    ! pole makes the response depend on where the nearest half ring falls
    ! (|b_shield| of the q = 2 surface not mesh converged at q0 1.55, A = 10),
    ! so the sheet closure shields them all; rat_auto = .false. keeps only the
    ! listed surfaces.  (The Frobenius closure takes its list as given.)
    if (rat_auto .and. .not. frob .and. trim(model) == 'ideal') &
        call add_resonant_surfaces()

    ! rings: uniform in u = sqrt(psin) plus the rational surfaces
    allocate (u_rat(nrat))
    u_rat = sqrt(psin_rat(1:nrat))
    if (frob) then
        ! Mercier index and shear of each rational, rings of the layers
        allocate (di_rat(nrat), dq_rat(nrat))
        block
            real(dp), allocatable :: extra(:), lab(:)
            real(dp) :: qa
            extra = u_rat
            do i = 1, nrat
                call mercier_di(g, psin_rat(i), nfine, di_rat(i), qa, dq_rat(i), st)
                if (st /= 0) error stop 'Mercier index'
                write (*, '(a,f10.6,a,f10.6,a,f10.6,a,es11.3)') 'rational q = ', &
                    q_rat(i), ': D_I = ', di_rat(i), ', mu = ', sqrt(-di_rat(i)), &
                    ', dq/dpsin = ', dq_rat(i)
                ! grade up to the uniform spacing, at most 0.45 of the gap to the
                ! neighbouring rationals
                qa = layer_xmax(psin_rat(i), &
                    [psin_rat(:i - 1), psin_rat(i + 1:nrat)], layer_delta, &
                    layer_ratio, layer_nfit, nrad)
                call layer_labels(psin_rat(i), layer_delta, layer_ratio, qa, &
                    (0.5_dp/nrad)**2, 1 - 2.0e-4_dp, lab)
                extra = [extra, sqrt(lab)]
            end do
            call ring_labels_merge(nrad, extra, min_gap, ul)
            call ring_labels_smooth(ul, 2*layer_ratio)
        end block
    else if (rat_pad > 0.0_dp) then
        ! rings at u_s +- rat_pad 2^-j/nrad (j < rat_grade): the mesh is graded
        ! geometrically toward every rational surface, symmetrically about it
        block
            real(dp), allocatable :: extra(:)
            integer :: jg
            extra = u_rat
            do jg = 0, rat_grade - 1
                extra = [extra, u_rat - rat_pad*0.5_dp**jg/nrad, &
                    min(u_rat + rat_pad*0.5_dp**jg/nrad, &
                    1.0_dp - 0.5_dp*rat_pad*0.5_dp**jg/nrad)]
            end do
            call ring_labels_merge(nrad, extra, min_gap, ul)
        end block
    else
        call ring_labels_merge(nrad, u_rat, min_gap, ul)
    end if
    nr = size(ul)
    allocate (ring_of(nrat))
    do i = 1, nrat
        ring_of(i) = minloc(abs(ul - u_rat(i)), 1)
    end do
    write (*, '(a,i0,a,i0)') 'rings: ', nr, ', rational rings: ', nrat

    call cpu_time(t0)
    call flux_mesh_build(ul, n_core, nth_max, min(order, 2), polar_map, fm, st)
    if (st /= 0) error stop 'mesh'
    call meridian_feec_build(fm%mm, order, fs, st)
    if (st /= 0) error stop 'FEEC'
    call ampere_n_build(fs, n, amp, st)
    if (st /= 0) error stop 'Ampere operator'
    call meridian_point_basis_init(fs, pb, st)
    call cpu_time(t1)
    write (*, '(a,i0,a,i0,a,i0,a,f10.2,a)') 'triangles ', &
        fs%mm%topology%n_triangles, ', unknowns ', amp%nn, ', boundary panels ', &
        amp%bd%n_panels, ', setup cpu ', t1 - t0, ' s'

    if (frob) then
        allocate (lays(nrat))
        do i = 1, nrat
            call layer_init(ul, psin_rat(i), abs(q_rat(i)), dq_rat(i), di_rat(i), &
                m_rat(i), layer_delta, layer_ratio, layer_nfit, lays(i), st)
            if (st /= 0) error stop 'layer fit rings'
        end do
    end if

    if (trim(model) == 'ideal') then
        opt%boundary_sheet = boundary_sheet
        opt%qlow = qlow
        opt%moved_current = edge_moved_current
        if (frob) then
            opt%res_m = m_rat(1:nrat)
            opt%res_eps = lays(1:nrat)%eps
        end if
        nsr = nr
        allocate (rings(nsr))
        do j = 2, nr
            call sample_ring_build(g, opt, fm, fs, pb, 0.5_dp*(ul(j - 1) + ul(j)), &
                j, np_sample, rings(j - 1), st)
            if (st /= 0) error stop 'half ring'
            if (frob) rings(j - 1)%active = .not. layer_inside(lays, ul, rings(j - 1)%u)
        end do
        call sample_ring_build(g, opt, fm, fs, pb, ul(nr), nr, np_sample, &
            rings(nsr), st)
        if (st /= 0) error stop 'boundary sample ring'
        allocate (avol(amp%nn))
        avol = (0.0_dp, 0.0_dp)
        write (*, '(a,i0,a,i0,a)') 'ideal model: ', nsr, ' sample rings of ', &
            np_sample, ' points'
    else if (trim(model) /= 'shield') then
        error stop 'model must be shield or ideal'
    end if

    ! ring angles and unit-sheet fields
    allocate (ra(nrat), area(nrat), b2(fs%sp(2)%ndof, nrat), f(amp%nn), a(amp%nn))
    do i = 1, nrat
        call surface_build(g, ul(ring_of(i))**2, nfine, surf, st)
        if (st /= 0) error stop 'rational surface'
        call ring_angles_set(surf%q, surf%theta_s(:, angle), surf%nu, ra(i))
        area(i) = ring_area(fm, ring_of(i))
        f = (0.0_dp, 0.0_dp)
        call ring_sheet_load(fs, pb, fm, ring_of(i), ra(i), n, m_rat(i), &
            (1.0_dp, 0.0_dp), f)
        call ampere_n_solve(amp, f, a, st)
        if (st /= 0) error stop 'sheet solve'
        call note_ampere_indicator()
        block
            complex(dp), allocatable :: bb(:)
            call ampere_n_curl(amp, a, bb)
            b2(:, i) = bb
        end block
    end do
    allocate (mm(nrat, nrat), lu(nrat, nrat), ipiv(nrat), phivac(nrat), lam(nrat), &
        phil3(nrat))
    do j = 1, nrat
        do i = 1, nrat
            mm(i, j) = ring_flux_harmonic(fs, pb, fm, ring_of(i), ra(i), n, &
                b2(:, j), m_rat(i))
        end do
    end do
    write (*, '(a,es10.2)') 'reciprocity |M - M^H|/|M| ', &
        maxval(abs(mm - conjg(transpose(mm))))/maxval(abs(mm))

    ! closure of the ideal model: sheets, their vacuum fields and constraints
    if (frob) then
        call layer_closure(g, ul, lays, nfine, angle, cl, st)
        if (st /= 0) error stop 'layer closure'
        allocate (b2cl(fs%sp(2)%ndof, cl%ns), jump_unit(nrat))
        do i = 1, cl%ns
            f = (0.0_dp, 0.0_dp)
            call ring_sheet_load(fs, pb, fm, cl%sheet_ring(i), cl%sheet_ra(i), n, &
                cl%sheet_m(i), (1.0_dp, 0.0_dp), f)
            call ampere_n_solve(amp, f, a, st)
            if (st /= 0) error stop 'edge sheet solve'
            call note_ampere_indicator()
            block
                complex(dp), allocatable :: bb(:)
                call ampere_n_curl(amp, a, bb)
                b2cl(:, i) = bb
            end block
        end do
        ! jump of d Phi/d psin across the rational of the unit rational sheet
        do i = 1, nrat
            phicl = closure_fluxes(fs, pb, fm, cl, n, b2(:, i))
            jump_unit(i) = layer_jump(lays(i), layer_block(phicl, i))
        end do
        if (resist) then
            call resistive_layers()
            call resistive_closure(lays, dlayer, cmplx(jump_unit, 0.0_dp, dp), cl)
        end if
    else if (trim(model) == 'ideal') then
        call closure_shield_init(ring_of, ra, m_rat(1:nrat), cl)
        b2cl = b2
    end if

    write (*, '(a,es10.3,a,i0,a,i0,a)') 'Ampere discrete residual-energy '// &
        'indicator (max over sheet solves) ', ampere_indicator, ', solves ', &
        ampere_solves, ', without indicator ', ampere_missing, &
        '; claim indicator, continuous error unavailable-exterior-remainder'
    open (newunit=ur, file=trim(prefix)//'_res.dat', status='replace')
    write (ur, '(a,es14.6,a,i0,a,i0,a)') '# ampere_discrete_energy_error_max ', &
        ampere_indicator, ' sheet_solves ', ampere_solves, ' without_indicator ', &
        ampere_missing, ' claim indicator (algebraic residual energy of the '// &
        'discrete Ampere solves) continuous_claim unavailable-exterior-remainder'
    write (ur, '(a)') '# phasing q psin m area[m^2/rad] Re,Im Phi_vac[T m^2/rad^2] '// &
        'Re,Im b_vac[T] Re,Im lambda[A] Re,Im b_shield[T] '// &
        'Re,Im Phi_vac(L3 surface path) (peak amplitudes)'
    if (resist) then
        open (newunit=ures, file=trim(prefix)//'_resist.dat', status='replace')
        write (ures, '(a)') '# phasing q psin m Re,Im Delta_layer[1/psin] '// &
            'L[psin] Q Re,Im A[T m^2/rad^2] Re,Im A/area[T] '// &
            'Re,Im jump[dPhi/dpsin] Re,Im jump(spot) |Phi_vac| |a1_in-a1_out|/|A| '// &
            'Re,Im b(0)/b0 Re,Im A b(0)/b0/area[T]'
    end if
    if (frob) then
        open (newunit=ulay, file=trim(prefix)//'_layer.dat', status='replace')
        write (ulay, '(a)') '# phasing q psin m D_I mu delta[psin] '// &
            'spot_x[psin] Re,Im jump[dPhi/dpsin] |a1/Phi|_in |a1/Phi|_out '// &
            'Re,Im lambda_in[A] Re,Im lambda_out[A] Re,Im lambda_eq[A] '// &
            'Re,Im jump_unit Re,Im jump_fit(spot)'
    end if
    allocate (w(set%ncoil))
    do ip = 1, nphase
        call coil_read_currents(trim(currents(ip)), set%ncoil, w, st)
        if (st /= 0) error stop 'cannot read currents'
        call cpu_time(t0)
        do i = 1, nrat
            call ring_points(fm, ring_of(i), xg, nds, thg, tg, tau, wq)
            if (allocated(bc)) deallocate (bc)
            allocate (bc(3, 1, size(thg)))
            call coil_set_modes(set, w, xg, nphi, [n], bc, st)
            if (st /= 0) error stop 'coil projection'
            allocate (bn(size(thg)))
            do kk = 1, size(thg)
                bn(kk) = (bc(1, 1, kk)*nds(1, kk) + bc(3, 1, kk)*nds(2, kk))*xg(1, kk)
            end do
            phivac(i) = ring_flux_from_values(ra(i), n, m_rat(i), thg, bn, wq)
            deallocate (bn)
            ! independent L3 path: coil field on the surface_build points and
            ! the phi_mn of kin6d_vacuum_surfaces
            call surface_build(g, ul(ring_of(i))**2, nfine, surf, st)
            allocate (xs(2, nfine), bs(3, 1, nfine), bb(2, nfine))
            xs(1, :) = surf%r
            xs(2, :) = surf%z
            call coil_set_modes(set, w, xs, nphi, [n], bs, st)
            bb(1, :) = bs(1, 1, :)
            bb(2, :) = bs(3, 1, :)
            call surface_harmonics(surf, n, bb, angle, m_rat(i), m_rat(i), &
                bmn1, phimn1)
            phil3(i) = phimn1(1)
            deallocate (xs, bs, bb)
        end do
        if (trim(model) == 'ideal') then
            ! coil field on the sample rings
            allocate (xr(2, nsr*np_sample), bcr(3, 1, nsr*np_sample))
            do j = 1, nsr
                xr(:, (j - 1)*np_sample + 1:j*np_sample) = rings(j)%x
            end do
            call coil_set_modes(set, w, xr, nphi, [n], bcr, st)
            if (st /= 0) error stop 'coil projection on sample rings'
            do j = 1, nsr
                rings(j)%bc = bcr(:, 1, (j - 1)*np_sample + 1:j*np_sample)
            end do
            deallocate (xr, bcr)
        end if
        if (trim(model) == 'ideal') then
            if (frob) then
                ! coil flux on the closure flux rings
                if (allocated(phicl)) deallocate (phicl)
                allocate (phicl(cl%nf))
                do k = 1, cl%nf
                    call ring_points(fm, cl%flux_ring(k), xg, nds, thg, tg, tau, wq)
                    if (allocated(bc)) deallocate (bc)
                    allocate (bc(3, 1, size(thg)), bn(size(thg)))
                    call coil_set_modes(set, w, xg, nphi, [n], bc, st)
                    if (st /= 0) error stop 'coil projection'
                    do kk = 1, size(thg)
                        bn(kk) = (bc(1, 1, kk)*nds(1, kk) + bc(3, 1, kk)*nds(2, kk)) &
                            *xg(1, kk)
                    end do
                    phicl(k) = ring_flux_from_values(cl%flux_ra(k), n, cl%flux_m(k), &
                        thg, bn, wq)
                    deallocate (bn)
                end do
            else
                phicl = phivac
            end if
        end if
        call cpu_time(t1)
        lu = mm
        lam = -phivac
        call zgesv(nrat, 1, lu, nrat, ipiv, lam, nrat, info)
        if (info /= 0) error stop 'shielding system singular'
        if (trim(model) == 'ideal') then
            call cpu_time(t0)
            if (allocated(lam_cl)) deallocate (lam_cl, phitot)
            allocate (lam_cl(cl%ns), phitot(cl%nf))
            call ideal_response(fs, pb, fm, g, amp, opt, n, rings, cl, b2cl, phicl, &
                avol, lam_cl, phitot, gmres_tol, gmres_maxit, info, its, res)
            call cpu_time(t1)
            write (*, '(a,i0,a,i0,a,es10.2,a,f10.1,a)') 'ideal response: info ', &
                info, ', GMRES iterations ', its, ', residual ', res, ', cpu ', &
                t1 - t0, ' s'
            if (frob) then
                call layer_report(ip)
            else
                lam = lam_cl
            end if
        end if
        call export_l7(ip)
        do i = 1, nrat
            write (ur, '(i3,f10.6,f14.10,i5,es16.8,10es16.8)') ip, q_rat(i), &
                ul(ring_of(i))**2, m_rat(i), area(i), phivac(i), &
                phivac(i)/area(i), lam(i), -mm(i, i)*lam(i)/area(i), phil3(i)
            write (*, '(a,i0,a,f8.4,3(a,es11.3),a,es9.2)') 'phasing ', ip, &
                ' q ', q_rat(i), ' |b_vac| ', abs(phivac(i))/area(i), &
                ' |lambda| ', abs(lam(i)), ' |b_shield| ', &
                abs(mm(i, i)*lam(i))/area(i), ' L3-path rel. diff ', &
                abs(phivac(i) - phil3(i))/abs(phil3(i))
        end do
        write (*, '(a,f10.2,a)') 'coil drive cpu ', t1 - t0, ' s'
    end do
    close (ur)
    if (frob) close (ulay)
    if (resist) close (ures)

contains

    subroutine note_ampere_indicator()
        ampere_solves = ampere_solves + 1
        if (amp%error%discrete_status /= 0) then
            ampere_missing = ampere_missing + 1
        else
            ampere_indicator = max(ampere_indicator, amp%error%discrete_energy_error)
        end if
    end subroutine note_ampere_indicator

    function layer_block(phi, l) result(pb2)
        !! Fluxes of layer l on its fit rings, (nfit, side).
        complex(dp), intent(in) :: phi(:)
        integer, intent(in) :: l
        complex(dp) :: pb2(lays(l)%nfit, 2)
        integer :: base, ll
        base = 0
        do ll = 1, l - 1
            base = base + 2*lays(ll)%nfit
        end do
        pb2(:, 1) = phi(base + 1:base + lays(l)%nfit)
        pb2(:, 2) = phi(base + lays(l)%nfit + 1:base + 2*lays(l)%nfit)
    end function layer_block

    subroutine layer_report(iph)
        !! Shielding current of every rational from the fitted outer solution:
        !! jump of d Phi/d psin at psin_s -+ spot/(n |q'|), the equivalent
        !! rational-ring sheet amplitude (stored in lam) and the layer data.
        integer, intent(in) :: iph
        complex(dp) :: cf(3, 2), jmp, jfit, ph(layer_nfit, 2), a1c(2)
        real(dp) :: xs, res1(2)
        integer :: l, side
        do l = 1, nrat
            ph = layer_block(phitot, l)
            cf = layer_coefficients(lays(l), ph)
            xs = spot/(n*lays(l)%dq)
            ! fitted outer slopes at spot (diagnostic)
            jfit = layer_slope(lays(l), cf(:, 1), xs) &
                + layer_slope(lays(l), cf(:, 2), xs)
            ! net layer current: jump of d Phi/d psin across the layer, carried
            ! to spot with the small-solution power |x|^(mu - 1/2)
            jmp = layer_kink(lays(l), ph, resist)*(xs/lays(l)%delta)**lays(l)%eps
            ! large-solution residual relative to the flux at the layer edge
            a1c = layer_large(lays(l), ph)
            do side = 1, 2
                res1(side) = abs(a1c(side))/max(abs(ph(1, side)), tiny(1.0_dp))
            end do
            lam(l) = jmp/jump_unit(l)
            write (ulay, &
                '(i3,f10.6,f14.10,i5,3es14.6,es12.4,2es16.8,2es11.3,10es16.8)') &
                iph, q_rat(l), psin_rat(l), m_rat(l), lays(l)%di, &
                sqrt(-lays(l)%di), lays(l)%delta, xs, jmp, res1, &
                lam_cl(2*l - 1), lam_cl(2*l), lam(l), jump_unit(l), jfit
            write (*, '(a,f8.4,a,es11.3,a,2es9.2)') '  layer q ', q_rat(l), &
                ' |lambda_eq| ', abs(lam(l)), ' large residual in/out ', res1
            if (resist) then
                block
                    complex(dp) :: arec, jedge
                    arec = 0.5_dp*(a1c(1) + a1c(2))
                    jedge = (lam_cl(2*l - 1) + lam_cl(2*l))*jump_unit(l)
                    write (ures, '(i3,f10.6,f14.10,i5,2es14.6,2es12.4,14es14.6)') &
                        iph, q_rat(l), psin_rat(l), m_rat(l), dlayer(l), &
                        lscale(l), qlay(l), arec, arec/area(l), jedge, jmp, &
                        abs(phivac(l)), abs(a1c(1) - a1c(2))/max(abs(arec), &
                        tiny(1.0_dp)), lcen(l), arec*lcen(l)/area(l)
                    write (*, '(a,f8.4,a,es11.3,a,es11.3,a,es11.3)') '  resistive q ', &
                        q_rat(l), ' |A/area| ', abs(arec)/area(l), ' |jump| ', &
                        abs(jedge), ' |A/Phi_vac| ', abs(arec/phivac(l))
                end block
            end if
        end do
    end subroutine layer_report

    subroutine export_l7(iph)
        !! L7: displacement, Eulerian/Lagrangian field and curvature coupling
        !! of the total response (coil + plasma + sheets) on the surfaces
        !! export_psin, written to prefix_l7_ph<iph>.dat with a .sha256 sidecar.
        integer, intent(in) :: iph
        type(perturbation_t), allocatable :: pr(:)
        type(surface_t) :: sx
        type(sample_ring_t) :: sr
        complex(dp), allocatable :: b2tot(:), bt(:), bf(:, :), bcx(:, :, :), dbe(:, :)
        character(len=512) :: fname
        integer :: ne, ie, js, s2, st2, ngood
        real(dp) :: ue

        ne = count(export_psin > 0.0_dp)
        if (ne == 0) return
        allocate (b2tot(fs%sp(2)%ndof), pr(ne))
        b2tot = (0.0_dp, 0.0_dp)
        if (trim(model) == 'ideal') then
            call ampere_n_curl(amp, avol, bt)
            b2tot = bt
            do s2 = 1, cl%ns
                b2tot = b2tot + lam_cl(s2)*b2cl(:, s2)
            end do
        else
            do s2 = 1, nrat
                b2tot = b2tot + lam(s2)*b2(:, s2)
            end do
        end if
        ngood = 0
        do ie = 1, ne
            ue = sqrt(export_psin(ie))
            js = nr
            do s2 = 2, nr
                if (ul(s2) >= ue) then
                    js = s2
                    exit
                end if
            end do
            call surface_build(g, export_psin(ie), export_np, sx, st2)
            if (st2 /= 0) error stop 'export surface'
            call sample_ring_build(g, opt, fm, fs, pb, ue, js, export_np, sr, st2)
            if (st2 /= 0) error stop 'export sample ring'
            allocate (bf(3, export_np), bcx(3, 1, export_np), dbe(3, export_np))
            call sample_field(fs, sr, b2tot, bf)
            call coil_set_modes(set, w, sr%x, nphi, [n], bcx, st2)
            if (st2 /= 0) error stop 'export coil projection'
            dbe = bf + bcx(:, 1, :)
            call perturbation_on_surface(g, opt, sx, n, dbe, pr(ngood + 1), st2)
            if (st2 == 0) then
                call perturbation_add_spectra(sx, pr(ngood + 1), -export_mmax, &
                    export_mmax)
                ngood = ngood + 1
            else
                write (*, '(a,f10.6,a,i0)') 'L7 export: psin ', export_psin(ie), &
                    ' skipped, status ', st2
            end if
            deallocate (bf, bcx, dbe)
        end do
        if (ngood == 0) return
        write (fname, '(a,a,i0,a)') trim(prefix), '_l7_ph', iph, '.dat'
        call perturbation_export(trim(fname), pr(1:ngood), 'gs_response_3d '// &
            trim(nml)//' model='//trim(model), 3, g%s_pol, sign(1.0_dp, g%ip), st2)
        write (*, '(a,i0,a,a,a,i0)') 'L7 export: ', ngood, ' surfaces to ', &
            trim(fname), ', status ', st2
    end subroutine export_l7

    subroutine read_profiles(path, tab)
        !! Columns psin, eta [Ohm m], rho [kg/m^3], Omega [rad/s].
        character(len=*), intent(in) :: path
        real(dp), allocatable, intent(out) :: tab(:, :)
        character(len=512) :: line
        real(dp) :: v(4), buf(4, 20000)
        integer :: u, ios, c
        open (newunit=u, file=path, status='old', action='read')
        c = 0
        do
            read (u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            read (line, *, iostat=ios) v
            if (ios /= 0) error stop 'profiles: bad line'
            c = c + 1
            buf(:, c) = v
        end do
        close (u)
        if (c < 2) error stop 'profiles: need at least two rows'
        tab = buf(:, 1:c)
    end subroutine read_profiles

    function prof_at(x) result(v)
        !! eta, rho, Omega at psin x (linear interpolation, clamped).
        real(dp), intent(in) :: x
        real(dp) :: v(3), t
        integer :: k
        k = 1
        do while (k < size(prof_tab, 2) - 1 .and. prof_tab(1, k + 1) < x)
            k = k + 1
        end do
        t = (x - prof_tab(1, k))/(prof_tab(1, k + 1) - prof_tab(1, k))
        t = min(1.0_dp, max(0.0_dp, t))
        v = (1 - t)*prof_tab(2:4, k) + t*prof_tab(2:4, k + 1)
    end function prof_at

    subroutine resistive_layers()
        !! Layer response Delta_layer (per psin) of every rational.
        type(layer_params_t) :: par
        real(dp) :: v(3), vp, bth, eps, pv, lx, rho
        complex(dp) :: gu
        integer :: l
        allocate (dlayer(nrat), lcen(nrat), lscale(nrat), qlay(nrat), gpsi(nrat))
        do l = 1, nrat
            v = prof_at(psin_rat(l))
            call layer_metrics(g, psin_rat(l), nfine, vp, bth, gpsi(l), st)
            if (st /= 0) error stop 'layer metrics'
            rho = v(2)
            if (inertia_q2) rho = rho*(1 + 2*q_rat(l)**2)
            ! static drive omega_drive in a plasma rotating with Omega:
            ! gamma = -i (omega - n Omega)
            par = layer_params_t(cmplx(0.0_dp, -(omega_drive - n*v(3)), dp), &
                eta_scale*v(1)/mu0_l, rho, layer_nu, abs(n*dq_rat(l))*bth*gpsi(l))
            call layer_scales(par, lx, eps, pv, gu)
            dlayer(l) = layer_delta_x(par, st, lcen(l))/gpsi(l)
            if (st /= 0) error stop 'layer solve'
            lscale(l) = lx*gpsi(l)
            qlay(l) = eps**(-2.0_dp/3)
            write (*, '(a,f8.4,a,es10.3,a,es10.3,a,es10.3,a,2es11.3,a,es10.3,a,f8.3)') &
                'layer q ', q_rat(l), ': eta ', v(1), ' rho ', rho, ' Omega ', v(3), &
                ' Delta[1/psin] ', dlayer(l), ' L[psin] ', lscale(l), ' Q ', qlay(l)
            write (*, '(a,2es11.3,a,es10.3)') '    b(0)/b0 ', lcen(l), &
                ' inertial width x_ii[psin] ', abs(par%gamma)*sqrt(mu0_l*rho) &
                /par%fprime*gpsi(l)
        end do
    end subroutine resistive_layers

    subroutine polar_map(u, theta, x)
        !! Flux-surface map: psin = u^2 on the ray of geometric angle theta.
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        real(dp) :: rho
        logical :: ok
        if (u <= 0.0_dp) then
            x = [g%r_axis, g%z_axis]
            return
        end if
        call flux_ray_point(g, u*u, theta, rho, ok)
        if (.not. ok) error stop 'flux surface ray'
        x = [g%r_axis + rho*cos(theta), g%z_axis + rho*sin(theta)]
    end subroutine polar_map



    subroutine add_resonant_surfaces()
        real(dp), allocatable :: px(:), qx(:)
        integer, allocatable :: mx(:)
        integer :: ia

        call resonant_surfaces_missing(g, n, psin_rat_min, psin_rat_max, nfine, &
            m_rat(1:nrat), px, qx, mx)
        do ia = 1, size(px)
            write (*, '(a,f10.6,a,f14.10,a,i0,a)') 'rational q = ', qx(ia), &
                ' at psin = ', px(ia), ', m = ', mx(ia), ' (added, rat_auto)'
        end do
        psin_rat = [psin_rat(1:nrat), px]
        q_rat = [q_rat(1:nrat), qx]
        m_rat = [m_rat(1:nrat), mx]
        nrat = nrat + size(px)
    end subroutine add_resonant_surfaces

end program gs_response_3d
