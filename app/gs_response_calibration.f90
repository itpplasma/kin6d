program gs_response_calibration
    !! Independent TC24 unit-sheet magnetic solve/readback. This executable
    !! does not solve a pressure-retaining coupled plasma response.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_result_io, only: run_t, result_file_t
    use kin6d_vacuum_surfaces, only: psi_grid_t, eqdsk_read, flux_ray_point, &
        surface_t, surface_build, find_surface_q, profiles_of_psin, angle_hamada
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build, ring_labels_merge
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_point_basis_t, meridian_point_basis_init
    use kin6d_mhd_ampere_n, only: ampere_n_t, ampere_n_build, ampere_n_solve, &
        ampere_n_curl
    use kin6d_mhd_response_n, only: ring_angles_t, ring_angles_set, &
        ring_flux_harmonic, ring_sheet_load, ring_area, resonant_surfaces_missing
    use kin6d_response_calibration, only: spot_weights, spot_jump, spot_distance
    implicit none
    character(512) :: eqdsk = '', argument, message
    integer :: n = 3, nrad = 64, nth_max = 512, order = 2, nfine = 1024
    real(dp) :: qabs = 3.0_dp, spot = 5.0e-4_dp, stencil_fraction = 0.1_dp
    namelist /spot_calibration/ eqdsk, n, nrad, nth_max, order, nfine, qabs, &
        spot, stencil_fraction
    type(psi_grid_t) :: g
    type(surface_t) :: surf, plus, minus
    type(flux_mesh_t) :: fm
    type(meridian_feec_t) :: fs
    type(meridian_point_basis_t) :: pb
    type(ampere_n_t) :: amp
    type(ring_angles_t) :: ra
    type(run_t) :: run
    type(result_file_t) :: output
    real(dp) :: root, qtarget, dq, xs, h, nodes(5, 2), weights(5, 2), area
    real(dp) :: fedge, dfedge, pedge, dpedge
    real(dp), allocatable :: labels(:), extra(:), rats(:), qs(:)
    integer, allocatable :: modes(:)
    complex(dp), allocatable :: f(:), a(:), field(:)
    complex(dp) :: flux(5, 2), self, jump
    integer :: u, ios, st, i, side, ring, js, m
    logical :: found, validate
    call run%start('gs_response_calibration')
    validate = .false.
    if (run%argument_count() > 1) error stop 'only --validate is supported'
    if (run%argument_count() == 1) then
        call run%argument(1, argument)
        if (trim(argument) /= '--validate') error stop 'unknown argument'
        validate = .true.
    end if
    if (.not. run%has_config()) error stop 'need --config calibration.nml'
    call run%open_config(u)
    read (u, nml=spot_calibration, iostat=ios, iomsg=message)
    call run%close_config(u, 'spot_calibration', ios, message)
    if (ios /= 0) error stop 'spot_calibration group required'
    if (n == 0 .or. nrad < 8 .or. nth_max < 16 .or. nfine < 64 &
        .or. order < 1 .or. order > 2) error stop 'invalid resolution/mode'
    if (.not. all(ieee_is_finite([qabs, spot, stencil_fraction]))) &
        error stop 'nonfinite control'
    if (qabs <= 0 .or. spot <= 0 .or. stencil_fraction <= 0 &
        .or. stencil_fraction >= 0.4_dp) error stop 'invalid spot controls'
    if (abs(n*qabs - nint(n*qabs)) > 1.0e-10_dp) error stop 'nonresonant q'
    call eqdsk_read(trim(eqdsk), g, st)
    if (st /= 0) error stop 'EQDSK read'
    call profiles_of_psin(g, 1.0_dp, fedge, dfedge, pedge, dpedge)
    if (.not. ieee_is_finite(pedge)) error stop 'nonfinite edge pressure'
    call surface_build(g, 0.5_dp, nfine, surf, st)
    if (st /= 0) error stop 'reference surface'
    qtarget = sign(qabs, surf%q)
    found = find_surface_q(g, qtarget, 0.05_dp, 0.999999_dp, 256, &
        1.0e-13_dp, root, st)
    if (.not. found .or. st /= 0) error stop 'target rational not found'
    m = nint(n*qtarget)
    h = 1.0e-5_dp
    call surface_build(g, root + h, nfine, plus, st)
    if (st /= 0) error stop 'shear plus'
    call surface_build(g, root - h, nfine, minus, st)
    if (st /= 0) error stop 'shear minus'
    dq = (plus%q - minus%q)/(2*h)
    if (.not. ieee_is_finite(dq)) error stop 'nonfinite shear'
    if (abs(dq) <= tiny(1.0_dp)) error stop 'zero shear'
    call spot_distance(n, dq, spot, xs, st)
    if (st /= 0) error stop 'spot distance'
    do side = 1, 2
        do i = 1, 5
            nodes(i, side) = root + merge(-xs, xs, side == 1) &
                + (i - 3)*stencil_fraction*xs
        end do
    end do
    if (minval(nodes) <= 0 .or. maxval(nodes) >= 1) error stop 'spot outside cut'
    call spot_weights(nodes, root, n, dq, spot, weights, st)
    if (st /= 0) error stop 'spot stencil'
    call resonant_surfaces_missing(g, n, 0.05_dp, 0.999999_dp, 256, &
        [integer ::], rats, qs, modes)
    extra = [sqrt(rats), sqrt(root), sqrt(nodes(:, 1)), sqrt(nodes(:, 2))]
    call ring_labels_merge(nrad, extra, 0.0_dp, labels)
    call run%open_resolved(u)
    write (u, nml=spot_calibration)
    call run%close_resolved(u)
    print '(a,f12.8,a,es14.6,a,i0,a,es14.6)', 'VALIDATED psin=', root, &
        ' spot_distance=', xs, ' rationals=', size(rats), ' p_edge[Pa]=', pedge
    if (validate) stop
    call flux_mesh_build(labels, 8, nth_max, order, polar_map, fm, st)
    if (st /= 0) error stop 'mesh'
    call meridian_feec_build(fm%mm, order, fs, st)
    if (st /= 0) error stop 'FEEC'
    call ampere_n_build(fs, n, amp, st)
    if (st /= 0) error stop 'Ampere'
    call meridian_point_basis_init(fs, pb, st)
    if (st /= 0) error stop 'point basis'
    allocate (f(amp%nn), a(amp%nn))
    f = 0
    js = minloc(abs(labels**2 - root), 1)
    call surface_build(g, labels(js)**2, nfine, surf, st)
    if (st /= 0) error stop 'sheet surface'
    call ring_angles_set(surf%q, surf%theta_s(:, angle_hamada), surf%nu, ra)
    call ring_sheet_load(fs, pb, fm, js, ra, n, m, (1.0_dp, 0.0_dp), f)
    call ampere_n_solve(amp, f, a, st)
    if (st /= 0) error stop 'unit sheet solve'
    call ampere_n_curl(amp, a, field)
    self = ring_flux_harmonic(fs, pb, fm, js, ra, n, field, m)
    area = ring_area(fm, js)
    do side = 1, 2
        do i = 1, 5
            ring = minloc(abs(labels**2 - nodes(i, side)), 1)
            call surface_build(g, labels(ring)**2, nfine, surf, st)
            if (st /= 0) error stop 'sample surface'
            call ring_angles_set(surf%q, surf%theta_s(:, angle_hamada), surf%nu, ra)
            flux(i, side) = ring_flux_harmonic(fs, pb, fm, ring, ra, n, field, m)
        end do
    end do
    jump = spot_jump(weights, flux)
    print '(a,i0,a,es14.6,a,4es18.9)', 'unknowns=', amp%nn, ' area=', area, &
        ' unit_flux,unit_jump=', self, jump
    call execute_command_line('mkdir -p output', exitstat=ios)
    if (ios /= 0) error stop 'output directory'
    call output%create(run)
    call output%attribute('claim', &
        'unit-sheet magnetic candidate; coupled solve unavailable')
    call output%attribute('derivatives', &
        'fixed samples linear JVP/VJP; map/shear unavailable')
    call output%write('psin', nodes, '1', 'spot sampling surfaces', &
        'sample,side', 'candidate')
    call output%write('flux_real', real(flux), 'T m2/rad2', 'unit-sheet flux real', &
        'sample,side', 'candidate')
    call output%write('flux_imag', aimag(flux), 'T m2/rad2', &
        'unit-sheet flux imaginary', &
        'sample,side', 'candidate')
    call output%write('unit_flux', [real(self), aimag(self)], 'T m2/rad2', &
        'unit sheet self flux', 'complex', 'candidate')
    call output%write('unit_jump', [real(jump), aimag(jump)], 'T m2/rad2', &
        'direct unit sheet spot derivative jump', 'complex', 'candidate')
    call output%write('area', area, 'm2/rad', 'ring area', 'candidate')
    call output%write('dq', dq, '1', 'signed q shear', 'candidate')
    call output%write('edge_pressure', pedge, 'Pa', 'retained imported edge pressure', &
        'candidate')
    call output%close()
contains
    subroutine polar_map(label, theta, point)
        real(dp), intent(in) :: label, theta
        real(dp), intent(out) :: point(2)
        real(dp) :: radius
        logical :: ok
        if (label <= 0) then
            point = [g%r_axis, g%z_axis]
            return
        end if
        call flux_ray_point(g, label*label, theta, radius, ok)
        if (.not. ok) error stop 'flux surface ray'
        point = [g%r_axis + radius*cos(theta), g%z_axis + radius*sin(theta)]
    end subroutine polar_map
end program gs_response_calibration
