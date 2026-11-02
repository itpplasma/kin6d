module fb3d_synthetic_field
    ! Frozen data of the synthetic nfp = 3 free-boundary case in the normal
    ! plane (q1, q2) of the helical-axis curve chart: uniform toroidal flux
    ! density b3 = phiedge/(pi a^2) on the reference disc of radius a, zero
    ! pressure, and the poloidal label chi with field-line rotation
    ! d theta_c/d q3 = -iota_VMEC(s), s = |q|^2/a^2 (theta_c counter-clockwise
    ! in (q1, q2), so the VMEC angle is pi/2 - theta_c):
    !   chi = sgn b3 a^2/2 int_0^s iota ds',  grad chi = sgn b3 iota(s) q.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_ideal_map, only: ideal_frozen_t
    implicit none
    private
    type, extends(ideal_frozen_t), public :: iota_field_t
        real(dp), allocatable :: s(:), iota(:), integral(:)
        real(dp) :: a2 = 1, b3 = 1, sgn = -1, centre(2) = 0
    contains
        procedure :: evaluate => iota_evaluate
    end type
contains
    subroutine iota_evaluate(field, point, chi, chi_gradient, b3, pressure)
        class(iota_field_t), intent(in) :: field
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: chi, chi_gradient(2), b3, pressure
        real(dp) :: s, w, iota, integral, y(2)
        integer :: k, n
        n = size(field%s)
        y = point - field%centre
        s = sum(y**2)/field%a2
        if (s >= field%s(n)) then
            iota = field%iota(n)
            integral = field%integral(n) + iota*(s - field%s(n))
        else
            k = 1
            do while (field%s(k + 1) < s)
                k = k + 1
            end do
            w = (s - field%s(k))/(field%s(k + 1) - field%s(k))
            iota = (1 - w)*field%iota(k) + w*field%iota(k + 1)
            integral = field%integral(k) + (s - field%s(k))*0.5_dp*(field%iota(k) + iota)
        end if
        chi = field%sgn*field%b3*field%a2/2*integral
        chi_gradient = field%sgn*field%b3*iota*y
        b3 = field%b3
        pressure = 0
    end subroutine
end module fb3d_synthetic_field

program fb3d_synthetic
    ! KIN6D 3D free boundary on the helical-axis curve chart for the public
    ! synthetic nfp = 3 modular-coil case (benchmarks/fb3d-synthetic): coils
    ! from the MAKEGRID file (libneo reader, one circuit, EXTCUR = 1), traces
    ! on the non-axisymmetric Gamma_c, field-period R2 exterior, Newton on the
    ! one residual of kin6d_ideal_free_boundary. Writes <out>.lcfs (boundary
    ! points x y z |B| on full-turn phi-hat samples) and <out>.json.
    ! Usage: fb3d_synthetic coils iota.dat phiedge out [nedge nchain ntor
    !        nphi_map nphi_bem nphi_surface sign max_iterations q1c q2c
    !        tangential_damping globalisation reference_lcfs kick_q1 mirror]
    ! (q1c, q2c): centre of the reference disc in the normal plane. sign = -1
    ! (default) is the coil field's rotation sense in the chart frame
    ! (benchmarks/fb3d-synthetic/iota_sign.py); tangential_damping defaults to
    ! 0 (off). globalisation: -1 (default) the ideal_fb_init default (2 for
    ! ntor >= 2 on this chart, 0 otherwise), 0 backtracking and Levenberg
    ! damping, 1 dogleg trust region on the Gauss-Newton model of the merit
    ! (ideal_fb_t%trust_region), 2 the trust region with the negative-curvature
    ! escape on stalls and saddles (ideal_fb_t%negative_curvature); the gate
    ! is the same. 0 or 1 is the opt-out for a wanted unstable equilibrium,
    ! which globalisation 2 never returns.
    ! reference_lcfs ('-' none): a previous <out>.lcfs of the same resolution
    ! whose pointwise distance is reported; kick_q1 (m, default 0): a rigid
    ! q1 shift of the start (q1 is the stellarator-symmetry-odd direction of
    ! the chart at phi-hat = 0), to probe symmetry-broken equilibria.
    ! mirror (default 0): 1 meshes the normal plane exactly symmetric under
    ! q1 -> -q1 (initialize_gs mirror_x), so that the discrete problem is
    ! stellarator symmetric at every resolution; 0 the direct triangulation.
    ! The JSON "state" block is the acceptance check beyond the residual:
    ! In_-(S_red) and its lowest eigenvalue (0 negative: a strict local
    ! minimum of the reduced energy, i.e. a stable discrete equilibrium), the
    ! reduced energy L, the stellarator-symmetry defect of the plasma map
    ! (ideal_fb_symmetry_defect, component 1), the axis Z at phi-hat = 0 (0
    ! for a symmetric state) and the distance to reference_lcfs. The vacuum
    ! oracle is benchmarks/fb3d-synthetic/vacuum_reference.py.
    use, intrinsic :: iso_fortran_env, only: dp => real64, output_unit
    !$ use omp_lib, only: omp_get_max_threads
    use fb3d_synthetic_field, only: iota_field_t
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs
    use kin6d_equilibrium_chart, only: curve_chart_t
    use kin6d_ideal_free_boundary, only: ideal_fb_t, ideal_fb_init, ideal_fb_solve, &
        ideal_fb_report_t, ideal_fb_shell_field, IDEAL_FB_OK, ideal_fb_gradient, &
        ideal_fb_reduced_spectrum, ideal_fb_symmetry_defect, ideal_fb_solve_shell_map
    use kin6d_ideal_free_boundary_bem3d, only: ideal_fb_use_bem3d
    use kin6d_ideal_fb_circuits, only: coil_circuit_basis_t, coil_circuit_basis_init, &
        ideal_fb_circuits_t, ideal_fb_circuits_build, ideal_fb_set_currents
    use kin6d_coils_3d, only: coil_set_t, coil_mu0
    use kin6d_makegrid_circuits, only: makegrid_circuits_read
    use fortfem_parametric_bem_high_order, only: bem_pair_options_t, bem_pair_report_t
    !$ use omp_lib, only: omp_get_max_threads
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 1.0_dp, eps = 0.15_dp
    real(dp), parameter :: rc = 0.2_dp, shift = 0.5_dp
    integer, parameter :: nfp = 3
    character(256) :: coils_path, iota_path, out, arg, reference_path
    type(curve_chart_t) :: chart
    type(iota_field_t) :: field
    type(gs_state_t) :: mesh
    type(ideal_fb_t) :: fb
    type(ideal_fb_report_t) :: report
    type(coil_set_t) :: coils
    type(coil_circuit_basis_t) :: basis
    type(ideal_fb_circuits_t) :: circuits
    type(bem_pair_options_t) :: options
    type(bem_pair_report_t) :: bem_report
    real(dp), allocatable :: contour(:, :), chain(:, :), turns(:)
    integer, allocatable :: circuit(:)
    real(dp) :: phiedge, aref, h, t_start, t_setup, t_traces, t_bem, t_solve, sgn
    real(dp) :: centre(2), tangential_damping, kick, energy, lambda_min
    real(dp) :: defect, boundary_defect, ref_rms, ref_max
    real(dp), allocatable :: gtmp(:)
    integer :: n_negative, mirrored, spectrum_status
    integer :: nedge, nchain, ntor, nphi_map, nphi_bem, nphi_surface, ncircuit
    integer :: st, j, threads, solve_status, max_iterations, globalisation, mirror
    nedge = 12; nchain = 16; ntor = 1; nphi_map = 8; nphi_bem = 18
    nphi_surface = 384; sgn = -1; max_iterations = 200
    ! Optional step regularisation of the Newton matrix along the reference
    ! tangent of the plasma nodes (ideal_fb_t%tangential_damping, faded by 100
    ! after every full step); it never enters the residual or the gate.
    tangential_damping = 0
    globalisation = -1; mirror = 0
    call get_command_argument(1, coils_path)
    call get_command_argument(2, iota_path)
    call get_command_argument(3, arg); read (arg, *) phiedge
    call get_command_argument(4, out)
    if (command_argument_count() >= 5) call int_arg(5, nedge)
    if (command_argument_count() >= 6) call int_arg(6, nchain)
    if (command_argument_count() >= 7) call int_arg(7, ntor)
    if (command_argument_count() >= 8) call int_arg(8, nphi_map)
    if (command_argument_count() >= 9) call int_arg(9, nphi_bem)
    if (command_argument_count() >= 10) call int_arg(10, nphi_surface)
    if (command_argument_count() >= 11) then
        call get_command_argument(11, arg); read (arg, *) sgn
    end if
    if (command_argument_count() >= 12) call int_arg(12, max_iterations)
    centre = 0
    if (command_argument_count() >= 14) then
        call get_command_argument(13, arg); read (arg, *) centre(1)
        call get_command_argument(14, arg); read (arg, *) centre(2)
    end if
    if (command_argument_count() >= 15) then
        call get_command_argument(15, arg); read (arg, *) tangential_damping
    end if
    if (command_argument_count() >= 16) call int_arg(16, globalisation)
    reference_path = '-'
    if (command_argument_count() >= 17) call get_command_argument(17, reference_path)
    kick = 0
    if (command_argument_count() >= 18) then
        call get_command_argument(18, arg); read (arg, *) kick
    end if
    if (command_argument_count() >= 19) call int_arg(19, mirror)
    threads = 1
    !$ threads = omp_get_max_threads()
    !$ threads = omp_get_max_threads()
    t_start = wall()
    call helical(chart)
    call read_iota(field)
    ! Reference disc of the toroidal flux phiedge at b3 = 1 T-equivalent
    ! density of the coil field on the axis (|B| ~ 1 T at R0 = 1 m).
    aref = sqrt(phiedge/pi)
    field%a2 = aref**2; field%b3 = phiedge/(pi*field%a2); field%sgn = sgn
    field%centre = centre
    allocate (contour(2, nedge), chain(2, nchain))
    do j = 1, nedge
        call gamma_c_shifted(real(j - 1, dp)/nedge, contour(:, j))
    end do
    do j = 1, nchain
        chain(:, j) = centre + aref*[cos(2*pi*(j - 1)/nchain), &
            sin(2*pi*(j - 1)/nchain)]
        chain(1, j) = chain(1, j) + shift
    end do
    h = 2*pi*rc/nedge
    if (mirror == 1) then
        call initialize_gs(mesh, contour, 0.5_dp*h*h, st, 2, curve=gamma_c_shifted, &
            interior_chain=chain, mirror_x=shift)
    else
        call initialize_gs(mesh, contour, 0.5_dp*h*h, st, 2, curve=gamma_c_shifted, &
            interior_chain=chain)
    end if
    if (st /= 0) error stop 'mesh'
    mesh%points(1, :) = mesh%points(1, :) - shift; chain(1, :) = chain(1, :) - shift
    call ideal_fb_init(fb, mesh, chain, chart, field, nfp, ntor, 0.0_dp, 0.0_dp, st, &
        mu0=coil_mu0, nphi=nphi_map)
    if (st /= IDEAL_FB_OK) error stop 'ideal_fb_init on the curve chart'
    fb%tangential_damping = tangential_damping
    if (globalisation < 0) then
        globalisation = merge(2, merge(1, 0, fb%trust_region), fb%negative_curvature)
    else
        fb%trust_region = globalisation == 1 .or. globalisation == 2
        fb%negative_curvature = globalisation == 2
    end if
    t_setup = wall() - t_start
    call makegrid_circuits_read(trim(coils_path), coils, circuit, turns, ncircuit, st)
    if (st /= 0) error stop 'MAKEGRID coils'
    call coil_circuit_basis_init(basis, coils, circuit, turns, ncircuit, 64, st)
    if (st /= 0) error stop 'circuit basis'
    call ideal_fb_circuits_build(fb, basis, circuits, st, nphi_surface=nphi_surface)
    if (st /= IDEAL_FB_OK) error stop 'circuit traces'
    call ideal_fb_set_currents(fb, circuits, [(1.0_dp, j=1, ncircuit)], st)
    if (st /= IDEAL_FB_OK) error stop 'currents'
    t_traces = wall() - t_start - t_setup
    print '(a,es12.4,a,es10.2)', 'F_vac (T m) ', fb%f_vac, '  closure ', &
        maxval(circuits%closure)
    call ideal_fb_use_bem3d(fb, nphi_bem, options, bem_report, st, periodic=.true.)
    if (st /= IDEAL_FB_OK) error stop 'R2 exterior'
    t_bem = wall() - t_start - t_setup - t_traces
    print '(a,i0,a,f9.2,a,i0)', 'R2 exterior: ', 2*fb%nedge*nphi_bem, ' panels, ', &
        t_bem, ' s, pairs ', bem_report%accepted_pairs
    print '(a,i0,a,i0,a,i0)', 'unknowns ', fb%n, '  plasma nodes ', fb%plasma%nodes, &
        '  shell nodes ', fb%ns
    if (kick /= 0) then
        fb%plasma%x(1, 1, :) = fb%plasma%x(1, 1, :) + kick
        call ideal_fb_solve_shell_map(fb, st)
        if (st /= IDEAL_FB_OK) error stop 'kicked shell map'
        print '(a,es10.2)', 'start kicked by q1 ', kick
    end if
    flush (output_unit)
    call ideal_fb_solve(fb, 1e-9_dp, max_iterations, report, solve_status, &
        verbose=.true., residual_tolerance=1e-8_dp)
    t_solve = wall() - t_start - t_setup - t_traces - t_bem
    print '(a,i0,a,i0,a,es12.4)', 'solve status ', solve_status, ' iterations ', &
        report%iterations, ' residual ', report%residual
    call write_lcfs()
    call state_check()
    call write_json()
contains
    subroutine int_arg(i, v)
        integer, intent(in) :: i
        integer, intent(out) :: v
        character(32) :: s
        call get_command_argument(i, s); read (s, *) v
    end subroutine

    subroutine helical(c)
        type(curve_chart_t), intent(out) :: c
        real(dp) :: ac(3, 0:4), as(3, 0:4)
        ac = 0; as = 0
        ac(1, 1) = r0; as(2, 1) = r0
        ac(1, 2) = eps/2; as(2, 2) = -eps/2
        ac(1, 4) = eps/2; as(2, 4) = eps/2
        as(3, 3) = eps
        call c%init(ac, as, nfp, st)
        if (st /= 0) error stop 'helical chart'
    end subroutine

    subroutine gamma_c_shifted(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = rc*[cos(2*pi*t), sin(2*pi*t)]
        point(1) = point(1) + shift
    end subroutine

    subroutine read_iota(f)
        ! VMEC iota.dat (s, iota on the full grid) and its running integral.
        type(iota_field_t), intent(inout) :: f
        real(dp) :: row(2)
        integer :: u, n, ios, k
        character(256) :: line
        open (newunit=u, file=trim(iota_path), status='old', action='read')
        n = 0
        do
            read (u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (line(1:1) == '#') cycle
            n = n + 1
        end do
        allocate (f%s(n), f%iota(n), f%integral(n))
        rewind (u)
        k = 0
        do
            read (u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (line(1:1) == '#') cycle
            read (line, *) row
            k = k + 1; f%s(k) = row(1); f%iota(k) = row(2)
        end do
        close (u)
        f%integral(1) = 0
        do k = 2, n
            f%integral(k) = f%integral(k - 1) + 0.5_dp*(f%iota(k) + f%iota(k - 1)) &
                *(f%s(k) - f%s(k - 1))
        end do
    end subroutine

    subroutine write_lcfs()
        ! Gamma_p points X = H(q(phi-hat), phi-hat) and |B| of the vacuum
        ! side at the Gamma_p nodes of the shell cells, full-turn samples.
        integer, parameter :: nout = 72
        real(dp), parameter :: corner(2, 10) = reshape([0.0_dp, 0.0_dp, 1.0_dp, &
            0.0_dp, 0.0_dp, 1.0_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 10])
        real(dp) :: phi, x(3), b(3), q(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        real(dp) :: tr
        logical, allocatable :: done(:)
        integer :: u, e, a, k, m, node
        open (newunit=u, file=trim(out)//'.lcfs', status='replace', action='write')
        write (u, '(a)') '# phi_hat x y z |B| (Gamma_p nodes, vacuum side)'
        allocate (done(fb%ns))
        do k = 0, nout - 1
            phi = 2*pi*k/nout
            done = .false.
            do e = 1, fb%ncell
                do a = 1, min(fb%nloc, 6)
                    node = fb%conn(a, e)
                    if (fb%role(node) /= 1 .or. done(node)) cycle
                    done(node) = .true.
                    q = 0
                    do m = 1, fb%nmode
                        tr = 1
                        if (m > 1) then
                            tr = cos((m/2)*nfp*phi)
                            if (mod(m, 2) == 1) tr = sin((m/2)*nfp*phi)
                        end if
                        q = q + fb%xs(:, m, node)*tr
                    end do
                    q(3) = q(3) + phi
                    call fb%plasma%chart%jet(q, x, d1, d2, d3)
                    call ideal_fb_shell_field(fb, e, corner(:, a), phi, q, b, st)
                    if (st /= 0) error stop 'shell field'
                    if (abs(q(1) - hypot(x(1), x(2))) + abs(q(2) - x(3)) > 1e-9_dp) &
                        error stop 'shell-field node does not match the boundary node'
                    write (u, '(5es24.15)') phi, x, norm2(b)
                end do
            end do
        end do
        close (u)
    end subroutine

    subroutine state_check()
        ! Acceptance quantities beyond the residual (see the header).
        integer :: nr, k
        real(dp), allocatable :: a(:, :), b(:, :)
        allocate (gtmp(fb%n))
        call ideal_fb_gradient(fb, gtmp, st, lagrangian=energy)
        call ideal_fb_reduced_spectrum(fb, n_negative, lambda_min, spectrum_status)
        call ideal_fb_symmetry_defect(fb, 1, defect, boundary_defect, mirrored, st)
        print '(a,i0,a,es11.3,a,es10.2,a,es10.2,a,i0)', 'state: In_-(S_red) ', &
            n_negative, ' min ', lambda_min, ' symmetry defect ', defect, &
            ' Gamma_p ', boundary_defect, ' mirrored nodes ', mirrored
        ref_rms = -1; ref_max = -1
        if (trim(reference_path) == '-') return
        call read_lcfs(trim(out)//'.lcfs', a, nr)
        call read_lcfs(trim(reference_path), b, k)
        if (nr /= k .or. nr == 0) return
        ref_max = sqrt(maxval(sum((a(2:4, :) - b(2:4, :))**2, dim=1)))
        ref_rms = sqrt(sum((a(2:4, :) - b(2:4, :))**2)/nr)
        print '(a,2es11.3)', 'distance to reference rms max (m)', ref_rms, ref_max
    end subroutine

    subroutine read_lcfs(path, rows, n)
        character(*), intent(in) :: path
        real(dp), allocatable, intent(out) :: rows(:, :)
        integer, intent(out) :: n
        character(512) :: line
        integer :: u, ios, k
        n = 0
        open (newunit=u, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        do
            read (u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (line(1:1) /= '#') n = n + 1
        end do
        allocate (rows(5, n))
        rewind (u)
        k = 0
        do
            read (u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (line(1:1) == '#') cycle
            k = k + 1
            read (line, *) rows(:, k)
        end do
        close (u)
    end subroutine

    subroutine write_json()
        ! axis_RZ_phi0: (R, Z) of the image of the plasma axis node (the
        ! extremum of the frozen label chi, where b_pol = 0) at phi-hat = 0.
        integer :: u, m
        real(dp) :: axis_rz(2), q(3), x(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        q = 0
        do m = 1, fb%nmode
            ! Mode m = 1 is n = 0, 2k is cos and 2k + 1 sin of k nfp phi-hat,
            ! all equal to 1 or 0 at phi-hat = 0.
            if (mod(m, 2) == 1 .and. m > 1) cycle
            q(1:2) = q(1:2) + fb%plasma%x(1:2, m, fb%plasma%axis_node)
        end do
        call fb%plasma%chart%jet(q, x, d1, d2, d3)
        axis_rz = [hypot(x(1), x(2)), x(3)]
        open (newunit=u, file=trim(out)//'.json', status='replace', action='write')
        write (u, '(a)') '{'
        write (u, '(a,i0,a)') ' "status": ', solve_status, ','
        write (u, '(a,i0,a)') ' "iterations": ', report%iterations, ','
        write (u, '(a,es14.6,a)') ' "residual": ', report%residual, ','
        write (u, '(a,es14.6,a)') ' "merit": ', report%merit, ','
        write (u, '(a,i0,a)') ' "unknowns": ', fb%n, ','
        write (u, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)') ' "nedge": ', nedge, &
            ', "nchain": ', nchain, ', "ntor": ', ntor, ', "nphi_map": ', nphi_map, &
            ', "nphi_bem": ', nphi_bem, ', "nphi_surface": ', nphi_surface, ','
        write (u, '(a,i0,a)') ' "panels": ', 2*fb%nedge*nphi_bem, ','
        write (u, '(a,i0,a)') ' "threads": ', threads, ','
        write (u, '(a,f12.4,a)') ' "setup_s": ', t_setup, ','
        write (u, '(a,f12.4,a)') ' "traces_s": ', t_traces, ','
        write (u, '(a,f12.4,a)') ' "exterior_s": ', t_bem, ','
        write (u, '(a,f12.4,a)') ' "solve_s": ', t_solve, ','
        write (u, '(a,f12.4,a)') ' "total_s": ', wall() - t_start, ','
        write (u, '(a,es14.6,a)') ' "f_vac": ', fb%f_vac, ','
        write (u, '(a,es14.6,a)') ' "tangential_damping": ', tangential_damping, ','
        write (u, '(a,i0,a)') ' "globalisation": ', globalisation, ','
        write (u, '(a,i0,a)') ' "mirror_mesh": ', mirror, ','
        write (u, '(a,i0,a)') ' "rejected": ', report%rejected, ','
        write (u, '(a,i0,a)') ' "factorisations": ', report%factorisations, ','
        write (u, '(a,i0,a)') ' "escapes": ', report%escapes, ','
        write (u, '(a,i0,a)') ' "curvature_steps": ', report%curvature_steps, ','
        write (u, '(a,es14.6,a)') ' "kick_q1": ', kick, ','
        write (u, '(a,es14.6,a)') ' "sign": ', sgn, ','
        write (u, '(a,2(es22.14,a))') ' "axis_RZ_phi0": [', axis_rz(1), ',', &
            axis_rz(2), '],'
        write (u, '(a)') ' "state": {'
        write (u, '(a,i0,a)') '  "spectrum_status": ', spectrum_status, ','
        write (u, '(a,i0,a)') '  "reduced_negative": ', n_negative, ','
        write (u, '(a,es22.14,a)') '  "reduced_min": ', lambda_min, ','
        write (u, '(a,es22.14,a)') '  "reduced_energy": ', energy, ','
        write (u, '(a,es14.6,a)') '  "symmetry_defect": ', defect, ','
        write (u, '(a,es14.6,a)') '  "symmetry_defect_gamma_p": ', boundary_defect, ','
        write (u, '(a,i0,a,i0,a)') '  "mirrored_nodes": ', mirrored, ', "nodes": ', &
            fb%plasma%nodes, ','
        write (u, '(a,es14.6,a)') '  "axis_Z_phi0": ', axis_rz(2), ','
        write (u, '(a,es14.6,a,es14.6)') '  "reference_rms": ', ref_rms, &
            ', "reference_max": ', ref_max
        write (u, '(a)') ' }'
        write (u, '(a)') '}'
        close (u)
    end subroutine

    real(dp) function wall() result(t)
        integer(8) :: c, r
        call system_clock(c, r)
        t = real(c, dp)/real(r, dp)
    end function
end program fb3d_synthetic
