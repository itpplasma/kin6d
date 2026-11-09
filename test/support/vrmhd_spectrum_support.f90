module vrmhd_spectrum_support
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use milca_gs_symmetry, only: gs_symmetry_t, GS_TRANSLATION, GS_CHART_CARTESIAN
    use milca_vrmhd_parameters, only: vrmhd_parameters_t
    use milca_vrmhd_section, only: vrmhd_section_t, vrmhd_section_build, &
        vrmhd_section_basis, vrmhd_section_residual, vrmhd_section_jacobian, &
        vrmhd_section_mass, VRMHD_EDGE_FIRST, VRMHD_EDGE_SECOND
    use milca_vrmhd_fourier, only: vrmhd_fourier_t, vrmhd_fourier_build, &
        vrmhd_fourier_data, vrmhd_fourier_jacobian, vrmhd_fourier_wavenumber
    use rm_mhd_oracles, only: uniform_waves, oracle_mu0, oracle_ok
    use vrmhd_structural_spectrum, only: structural_spectrum, bloch_projection, &
        bloch_rows, bloch_sandwich, bloch_adjoint, check_bloch_action, &
        harmonic_block
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes, meridian_load
    use fortsparse, only: csc_t, csc_matvec, sparse_solve_once, &
        fortsparse_status_t, status_ok
    implicit none
    private
    public :: spectrum, require, spectrum_result_t, edge_cell_dofs, block_dofs
    type :: spectrum_result_t
        integer :: finite = 0, zeros = 0, positive = 0, unmatched = 0
        integer :: magnetic = 0, velocity = 0, mixed = 0, shear_nulls = 0
        integer :: svd_zeros = 0, energetic_nulls = 0
        integer :: static_magnetic = 0, static_velocity = 0, static_mixed = 0
        integer :: signed_counts(3, 2) = 0
        real(dp) :: reality = 0, backward_error = 0, equilibrium = 0
        real(dp) :: analytic_null_defect = 0, shear_ratio = 0, seconds = 0
        logical :: real_ok = .false., count_ok = .false., null_split_ok = .false.
    end type
    real(dp) :: background_field(3)
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (ok) return
        print '(a)', message
        error stop 'RM2 linear spectrum gate'
    end subroutine

    pure integer function edge_cell_dofs(ka, edge_kind) result(count)
        !! A coefficients per periodic cell: three edges, two triangles and
        !! the orbit Lagrange nodes of the magnetic complex.
        integer, intent(in) :: ka, edge_kind
        if (edge_kind == VRMHD_EDGE_FIRST) then
            count = 3*ka + 2*ka*(ka - 1) + ka**2
        else
            count = 3*(ka + 1) + 2*(ka**2 - 1) + (ka + 1)**2
        end if
    end function

    pure integer function block_dofs(ks, ka, edge_kind, repeat, velocity_degree) &
            result(count)
        !! Bloch block: three thermodynamic fields, three velocities, A and J.
        integer, intent(in) :: ks, ka, edge_kind, repeat
        integer, intent(in), optional :: velocity_degree
        integer :: kv
        kv = ks
        if (present(velocity_degree)) kv = velocity_degree
        count = repeat**2*(3*ks**2 + 3*kv**2 + 2*edge_cell_dofs(ka, edge_kind))
    end function

    subroutine spectrum(n, ka, ks, nb, b0, exact, audit, frequencies, &
            optical_fraction, unmatched, result, perturbed, structural, mass_condition, &
            edge_kind, harmonic, orbit_scale, orientation, patch, velocity_degree)
        integer, intent(in) :: n, ka, ks, nb
        real(dp), intent(in) :: b0(3), exact(3)
        logical, intent(in) :: audit
        real(dp), intent(out) :: frequencies(3), optical_fraction
        real(dp), intent(out), optional :: mass_condition
        integer, intent(out) :: unmatched
        type(spectrum_result_t), intent(out), optional :: result
        logical, intent(in), optional :: perturbed
        logical, intent(in), optional :: structural
        integer, intent(in), optional :: edge_kind, velocity_degree
        !! Harmonic n of the section x Fourier operator with kz = orientation n/orbit_scale
        !! along the invariant direction (structural path only).
        integer, intent(in), optional :: harmonic, orientation
        real(dp), intent(in), optional :: orbit_scale
        !! Assemble only 2 repeat cells per direction: the Bloch block of a
        !! translation-invariant operator does not depend on the period count.
        logical, intent(in), optional :: patch
        type(spectrum_result_t) :: measured
        type(vrmhd_section_t) :: reference_op, assembly
        type(vrmhd_fourier_t) :: fop
        integer :: repeat, cell, ix, iy, tick0, tick1, tick_rate, cells
        real(dp), allocatable :: reference_vertices(:, :), orbit(:, :), nodes(:, :)
        real(dp) :: h, extent, kz, frame(3, 3), local(3)
        real(dp), parameter :: mass_shift = 2.0_dp**30
        logical :: harmonic_mode, small
        real(dp) :: displacement(2, 0:1, 0:1)
        logical :: study, distort
        integer :: ns, na, nv, kv, nf, te, ti, off(9), order
        type(mesh_2d_t) :: mesh
        type(meridian_mesh_t) :: mm
        type(vrmhd_section_t) :: op
        type(vrmhd_parameters_t) :: p
        type(csc_t) :: jac, time_mass
        integer, allocatable :: column(:)
        complex(dp), allocatable :: row_phase(:)
        complex(qp), allocatable :: column_norm(:)
        logical :: reduced
        real(dp), allocatable :: u(:), zero(:), source(:, :, :), r(:), unused(:)
        real(dp), allocatable :: assembly_source(:, :, :)
        real(dp), allocatable :: mr(:), mi(:), basis(:, :), xy(:, :)
        complex(dp), allocatable :: load(:)
        type(fortsparse_status_t) :: st
        complex(dp), allocatable :: transform(:, :), coefficients(:), periodic(:, :)
        real(dp), allocatable :: phase_x(:)
        complex(dp) :: stiffness(nb, nb), mass(nb, nb), alpha(nb), beta(nb)
        complex(dp) :: vectors(nb, nb), dummy(1, 1), dummy_right(1, 1), work(16*nb), lambda(nb)
        complex(dp) :: original_k(nb, nb), original_m(nb, nb), residual(nb)
        complex(dp) :: diagnostic_alpha(nb), diagnostic_beta(nb)
        complex(dp) :: optical_vector(nb), projected(nb), nulls(nb, 3)
        complex(dp) :: phase, field(3), mean_field(3)
        real(dp) :: rwork(8*nb), length, x(2), start(2), tangent(2), midpoint(2)
        real(dp) :: weight, direction, volume, energy, max_energy, fraction, scale
        real(dp) :: positive(nb), tolerance, reality, equilibrium_defect
        real(dp) :: backward_error, diagnostic_reality, column_scale(nb), row_scale(nb)
        real(dp) :: magnetic_energy(3), optical_energy(3)
        integer :: n0, nn, n1, a0, j0, t, i, j, q, edge, group
        integer :: status, info, finite_count, zero_count, wave_count, m, idx, group0
        integer :: kind, ke, ninterior, k0
        interface
            subroutine zggev(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, &
                    vl, ldvl, vr, ldvr, work, lwork, rwork, info)
                import :: dp
                character, intent(in) :: jobvl, jobvr
                integer, intent(in) :: n, lda, ldb, ldvl, ldvr, lwork
                complex(dp), intent(inout) :: a(lda, *), b(ldb, *)
                complex(dp), intent(out) :: alpha(*), beta(*), vl(ldvl, *)
                complex(dp), intent(out) :: vr(ldvr, *), work(*)
                real(dp), intent(out) :: rwork(*)
                integer, intent(out) :: info
            end subroutine
        end interface
        call system_clock(tick0, tick_rate)
        background_field = b0
        study = present(result); distort = .false.
        if (present(perturbed)) distort = perturbed
        repeat = 1
        if (distort) repeat = 2
        kind = VRMHD_EDGE_SECOND
        if (present(edge_kind)) kind = edge_kind
        kv = ks
        if (present(velocity_degree)) kv = velocity_degree
        ns = (repeat*ks)**2; nv = (repeat*kv)**2
        na = repeat**2*edge_cell_dofs(ka, kind)
        nf = 3*ns + 3*nv; te = ns + 3*nv; ti = te + ns
        off = [0, ns, ns + nv, ns + 2*nv, te, ti, nf, nf + na, nf + 2*na]
        call require(nb == nf + 2*na, 'complete study block size')
        length = 2*acos(-1.0_dp)
        unmatched = 0
        harmonic_mode = present(harmonic)
        call require(harmonic_mode .eqv. present(orbit_scale), &
            'harmonic and orbit scale')
        small = .false.
        if (present(patch)) small = patch
        h = length/n; cells = n; extent = length
        if (small) then
            cells = 2*repeat; extent = cells*h
        end if
        call mesh%create_rectangular(cells + 1, cells + 1, 0.0_dp, extent, &
            0.0_dp, extent)
        reference_vertices = mesh%vertices
        if (distort) then
            ! The fixed pattern repeats across both periodic boundaries.
            displacement(:,0,0) = [0.0_dp, 0.0_dp]
            displacement(:,1,0) = [0.13_dp, -0.11_dp]
            displacement(:,0,1) = [-0.09_dp, 0.16_dp]
            displacement(:,1,1) = [0.18_dp, 0.07_dp]
            do i = 1, mesh%n_vertices
                ix = nint(reference_vertices(1,i)/h)
                iy = nint(reference_vertices(2,i)/h)
                mesh%vertices(:,i) = mesh%vertices(:,i) &
                    + h*displacement(:,modulo(ix,2),modulo(iy,2))
            end do
        end if
        call meridian_mesh_straight(mesh, 1, mm, status)
        call require(status == 0, 'periodic wave mesh')
        if (harmonic_mode) then
            call require(kind == VRMHD_EDGE_SECOND .and. .not. study &
                .and. present(structural), &
                'harmonic block uses the admitted section interface')
            call vrmhd_fourier_build(mm, &
                gs_symmetry_t(GS_TRANSLATION, GS_CHART_CARTESIAN, 0.0_dp), &
                ka, [harmonic], orbit_scale, fop, status, &
                scalar_degree=ks, velocity_degree=kv)
            call require(status == 0 .and. harmonic > 0, 'periodic wave harmonic')
            op = fop%section
            kz = vrmhd_fourier_wavenumber(fop, 1)
            if (present(orientation)) kz = sign(kz, real(orientation, dp))
        else
            call vrmhd_section_build(mm, &
                gs_symmetry_t(GS_TRANSLATION, GS_CHART_CARTESIAN, 0.0_dp), &
                ka, 1, op, status, scalar_degree=ks, study_pairing=study, &
                edge_kind=kind, velocity_degree=kv)
        end if
        call require(status == 0, 'periodic wave section')
        reduced = .false.
        if (present(structural)) reduced = structural
        if (reduced .and. .not. harmonic_mode) then
            ! At the uniform state on affine cells every Jacobian and time-mass
            ! integrand is a product of two basis polynomials with constant
            ! coefficients, so degree 2 max(ks, k0) is integrated exactly.
            call require(mm%geo_order == 1, 'affine cells for the reduced quadrature')
            call vrmhd_section_build(mm, &
                gs_symmetry_t(GS_TRANSLATION, GS_CHART_CARTESIAN, 0.0_dp), &
                ka, 1, assembly, status, scalar_degree=ks, study_pairing=study, &
                edge_kind=kind, quad_degree=2*max(ks, kv, op%fs%k0), &
                velocity_degree=kv)
            call require(status == 0, 'periodic wave section, exact bilinear quadrature')
            call require(assembly%ndof == op%ndof .and. all(assembly%map == op%map), &
                'reduced-quadrature section shares the coefficient numbering')
        else
            assembly = op
        end if
        ke = op%fs%ke; ninterior = op%fs%nln - 3*ke; k0 = op%fs%k0
        if (distort) then
            mesh%vertices = reference_vertices
            call meridian_mesh_straight(mesh, 1, mm, status)
            call require(status == 0, 'reference mesh for periodic identification')
            call vrmhd_section_build(mm, &
                gs_symmetry_t(GS_TRANSLATION, GS_CHART_CARTESIAN, 0.0_dp), &
                ka, 1, reference_op, status, scalar_degree=ks, study_pairing=.true., &
                edge_kind=kind, velocity_degree=kv)
            call require(status == 0, 'reference coefficient topology')
            mm = op%fs%mm
        end if
        n0 = op%scalar_fs%n0; nn = op%fs%nn; n1 = op%fs%sp(1)%ndof
        a0 = op%offset(7); j0 = op%offset(8)
        p%mi = 1; p%mu0 = oracle_mu0
        allocate (u(op%ndof), zero(op%ndof), r(op%ndof), unused(op%ndof), &
            source(6, op%fs%nq, mm%topology%n_triangles), mr(op%ndof), mi(op%ndof), &
            transform(op%ndof, nb), coefficients(op%ndof), phase_x(op%ndof), &
            basis(36, op%nlocal))
        allocate (assembly_source(6, assembly%fs%nq, mm%topology%n_triangles))
        u = 0; zero = 0; source = 0; transform = 0; assembly_source = 0
        u(:n0) = 1; u(op%offset(5) + 1:op%offset(7)) = 0.3_dp
        do group = 1, 6
            order = ks
            if (group >= 2 .and. group <= 4) order = kv
            if (group >= 2 .and. group <= 4 .and. allocated(op%velocity_fs)) then
                if (distort) then
                    call meridian_lagrange_nodes(reference_op%velocity_fs, xy)
                else
                    call meridian_lagrange_nodes(op%velocity_fs, xy)
                end if
            else if (distort) then
                call meridian_lagrange_nodes(reference_op%scalar_fs, xy)
            else
                call meridian_lagrange_nodes(op%scalar_fs, xy)
            end if
            do i = 1, size(xy, 2)
                group0 = nodal_group(xy(:, i), order, length/n, repeat)
                phase = exp(cmplx(0.0_dp, xy(1, i), dp))
                transform(op%offset(group) + i, off(group) + group0) = phase
                phase_x(op%offset(group) + i) = xy(1, i)
            end do
            deallocate (xy)
        end do
        if (distort) then
            call meridian_lagrange_nodes(reference_op%fs, xy)
        else
            call meridian_lagrange_nodes(op%fs, xy)
        end if
        do i = 1, op%fs%n0
            group0 = repeat**2*(3*ke + 2*ninterior) &
                + nodal_group(xy(:, i), k0, length/n, repeat)
            phase = exp(cmplx(0.0_dp, xy(1, i), dp))
            transform(a0 + nn + i, nf + group0) = phase
            transform(j0 + nn + i, nf + na + group0) = phase
            phase_x([a0 + nn + i, j0 + nn + i]) = xy(1, i)
        end do
        do t = 1, mm%topology%n_triangles
            do i = 1, 3
                j = mod(i, 3) + 1
                start = reference_vertices(:, mm%topology%triangles(i, t))
                tangent = reference_vertices(:, mm%topology%triangles(j, t)) - start
                midpoint = start + tangent/2
                group = 3
                if (abs(tangent(2)) < epsilon(length)) group = 1
                if (abs(tangent(1)) < epsilon(length)) group = 2
                direction = sign(1.0_dp, tangent(1))
                if (group == 2) direction = sign(1.0_dp, tangent(2))
                do m = 1, ke
                    idx = (i - 1)*ke + m
                    edge = op%fs%sp(1)%dof(idx, t)
                    phase = real(op%fs%sp(1)%sgn(idx, t), dp)*direction**m &
                        *norm2(tangent)*exp(cmplx(0.0_dp, midpoint(1), dp))
                    cell = 1 + modulo(nint(midpoint(1)/h - 0.5_dp), repeat) &
                        + repeat*modulo(nint(midpoint(2)/h - 0.5_dp), repeat)
                    if (group == 1) cell = 1 + modulo(nint(midpoint(1)/h - 0.5_dp), repeat) &
                        + repeat*modulo(nint(midpoint(2)/h), repeat)
                    if (group == 2) cell = 1 + modulo(nint(midpoint(1)/h), repeat) &
                        + repeat*modulo(nint(midpoint(2)/h - 0.5_dp), repeat)
                    group0 = (cell - 1)*3*ke + (group - 1)*ke + m
                    transform(a0 + edge, nf + group0) = phase
                    transform(j0 + edge, nf + na + group0) = phase
                    phase_x([a0 + edge, j0 + edge]) = midpoint(1)
                end do
            end do
            midpoint = sum(reference_vertices(:, mm%topology%triangles(:,t)), dim=2)/3
            do i = 3*ke + 1, op%fs%nln
                edge = op%fs%sp(1)%dof(i, t)
                ix = modulo((t-1)/2, cells); iy = (t-1)/(2*cells)
                cell = 1 + modulo(ix,repeat) + repeat*modulo(iy,repeat)
                group0 = repeat**2*3*ke + ((cell-1)*2 + mod(t-1,2)) &
                    *ninterior + i - 3*ke
                phase = (length/n)*exp(cmplx(0.0_dp, midpoint(1), dp))
                transform(a0 + edge, nf + group0) = phase
                transform(j0 + edge, nf + na + group0) = phase
                phase_x([a0 + edge, j0 + edge]) = midpoint(1)
            end do
        end do
        call require(all(count(abs(transform) > 0, dim=2) == 1), 'every Bloch coefficient included once')
        call require(all(sum(abs(transform), dim=1) > 0), 'every Bloch column represented')
        call bloch_rows(transform, column, row_phase, column_norm)
        call meridian_load(op%fs, 1, equilibrium_potential, load, status)
        call require(status == 0, 'analytic uniform-field potential load')
        call sparse_solve_once(op%fs%mass(1), real(load, dp), u(a0 + 1:j0), st)
        call require(status_ok(st), 'uniform-field potential mass projection')
        ! Entropy is pressure balanced; the gauge is the discrete gradient of exp(i x).
        nulls = 0
        nulls(1:ns, 1) = 1; nulls(te + 1:ti, 1) = -0.6_dp
        nulls(te + 1:ti, 2) = 1; nulls(ti + 1:nf, 2) = -1
        if (distort) then
            deallocate(xy)
            call meridian_lagrange_nodes(op%fs, xy)
        end if
        coefficients = 0
        coefficients(a0 + 1:j0) = cmplx( &
            csc_matvec(op%fs%grad%a, cos(xy(1, :))), &
            csc_matvec(op%fs%grad%a, sin(xy(1, :))), dp)
        nulls(:, 3) = cmplx(bloch_adjoint(column, row_phase, &
            cmplx(coefficients,kind=qp), nb)/column_norm, kind=dp)
        ! Translation invariance makes each Bloch block an invariant periodic subspace.
        call vrmhd_section_residual(op, p, u, zero, zero, zero, source, r, unused, status)
        call require(status == 0, 'uniform equilibrium residual')
        if (small) then
            ! Over a patch the Bloch phases do not cancel; the equilibrium is
            ! the periodic (k = 0) projection of the assembled residual.
            allocate(periodic, mold=transform)
            do i = 1, op%ndof
                periodic(i, :) = transform(i, :)*exp(cmplx(0.0_dp, -phase_x(i), dp))
            end do
            equilibrium_defect = maxval(abs(bloch_projection(periodic,cmplx(r, &
                kind=qp))))
        else
            equilibrium_defect = maxval(abs(bloch_projection(transform,cmplx(r, &
                kind=qp))))
        end if
        print '(a,i4,es20.12)', 'periodic equilibrium residual: ', n, equilibrium_defect
        call require(equilibrium_defect < 1.0e-12_dp, 'uniform periodic equilibrium')
        if (harmonic_mode) then
            call harmonic_frame(frame)
            print '(a,9f8.4)', 'B0 / grad x / grad y / orbit wave direction: ', frame
            call require(maxval(abs(frame(:, 1) - b0)) <= 1.0e-12_dp*norm2(b0), &
                'equilibrium field is the oracle field')
            call require(maxval(abs(frame(:, 2:3) - reshape([1, 0, 0, 0, 1, 0], &
                [3, 2]))) <= 1.0e-12_dp, 'section axes are the oracle x and y')
            call require(maxval(abs(local - [0, 0, 1])) <= 1.0e-12_dp, &
                'invariant direction is the oracle z')
            call vrmhd_fourier_jacobian(fop, p, u, zero, source, 0.0_dp, jac, status)
            call require(status == 0, 'assembled harmonic linearization')
            original_k = harmonic_block(jac, transform, kz)
            ! A power-of-two shift isolates the time mass without cancellation.
            call vrmhd_fourier_jacobian(fop, p, u, zero, source, mass_shift, jac, &
                status)
            call require(status == 0, 'assembled harmonic time mass')
            original_m = (harmonic_block(jac, transform, kz) - original_k)/mass_shift
            call structural_spectrum(op, assembly, p, u, source, transform, column, &
                row_phase, original_k, original_m, ns, na, n, repeat, ks, b0, &
                exact, frequencies, &
                mass_condition, fourier=fop, kz=kz)
            optical_fraction = 0
            return
        end if
        call vrmhd_section_jacobian(assembly, p, u, zero, assembly_source, 0.0_dp, jac, &
            status)
        call require(status == 0, 'assembled production linearization')
        stiffness = bloch_sandwich(jac, column, row_phase, nb)
        call vrmhd_section_jacobian(assembly, p, u, zero, assembly_source, 1.0_dp, &
            time_mass, status, mass_only=.true.)
        call require(status == 0, 'assembled production time mass')
        mass = bloch_sandwich(time_mass, column, row_phase, nb)
        ! The production quadrature confirms the exact reduced-degree assembly.
        call check_bloch_action(op, p, u, column, row_phase, stiffness, .false., &
            'Jacobian')
        call check_bloch_action(op, p, u, column, row_phase, mass, .true., 'time mass')
        original_k = stiffness; original_m = mass
        if (present(structural)) then
            if (structural) then
                call structural_spectrum(op, assembly, p, u, assembly_source, &
                    transform, column, row_phase, original_k, original_m, ns, na, n, &
                    repeat, ks, b0, exact, frequencies, mass_condition)
                optical_fraction = 0
                return
            end if
        end if
        call require(kv == ks, 'legacy pairing diagnostics use a shared fluid space')
        ! Balance magnetic units before QZ; this is an invertible pencil equivalence.
        column_scale = 1; row_scale = 1
        column_scale(6*ns + 1:6*ns + na) = sqrt(p%mu0)
        column_scale(6*ns + na + 1:) = 1/sqrt(p%mu0)
        row_scale(6*ns + 1:) = 1/sqrt(p%mu0)
        do j = 1, nb
            stiffness(:, j) = row_scale*original_k(:, j)*column_scale(j)
            mass(:, j) = row_scale*original_m(:, j)*column_scale(j)
        end do
        call zggev('N', 'V', nb, stiffness, nb, mass, nb, alpha, beta, dummy, 1, &
            vectors, nb, work, size(work), rwork, info)
        call require(info == 0, 'assembled generalized eigenproblem')
        do j = 1, nb
            vectors(:, j) = column_scale*vectors(:, j)
        end do
        scale = 1; tolerance = 1024*epsilon(scale)
        finite_count = 0; zero_count = 0; wave_count = 0; reality = 0
        lambda = 0; backward_error = 0
        do j = 1, nb
            if (abs(beta(j)) < tolerance*(abs(alpha(j)) + abs(beta(j)))) cycle
            finite_count = finite_count + 1
            lambda(j) = alpha(j)/beta(j)
            residual = matmul(original_k, vectors(:, j)) &
                - lambda(j)*matmul(original_m, vectors(:, j))
            backward_error = max(backward_error, sqrt(sum(abs(residual)**2)) &
                /((sqrt(sum(abs(original_k)**2)) &
                + abs(lambda(j))*sqrt(sum(abs(original_m)**2))) &
                *sqrt(sum(abs(vectors(:, j))**2))))
            print '(a,i3,4es22.14)', 'lambda and omega^2, real / imaginary: ', &
                j, real(lambda(j), dp), aimag(lambda(j)), &
                real(-lambda(j)**2, dp), aimag(-lambda(j)**2)
            reality = max(reality, abs(aimag(-lambda(j)**2)))
            if (abs(lambda(j)) < tolerance) zero_count = zero_count + 1
            if (aimag(lambda(j)) <= tolerance) cycle
            wave_count = wave_count + 1
            positive(wave_count) = aimag(lambda(j))
        end do
        print '(a,i4,3i4,es20.12)', 'finite / zero / positive modes, Im(omega^2): ', &
            n, finite_count, zero_count, wave_count, reality
        print '(a,es20.12)', 'maximum generalized eigenpair backward error: ', backward_error
        stiffness = original_k; mass = original_m
        call zggev('N', 'N', nb, stiffness, nb, mass, nb, diagnostic_alpha, &
            diagnostic_beta, dummy, 1, dummy_right, 1, work, size(work), rwork, info)
        call require(info == 0, 'diagnostic SI-unit pencil')
        diagnostic_reality = 0
        do j = 1, nb
            if (abs(diagnostic_beta(j)) < tolerance &
                *(abs(diagnostic_alpha(j)) + abs(diagnostic_beta(j)))) cycle
            diagnostic_reality = max(diagnostic_reality, &
                abs(aimag(-(diagnostic_alpha(j)/diagnostic_beta(j))**2)))
        end do
        print '(a,es20.12)', 'diagnostic SI-unit omega^2 imaginary defect: ', &
            diagnostic_reality
        if (wave_count > 0) print '(a,*(es20.12))', &
            'positive angular frequencies [1/s], unsorted: ', positive(:wave_count)
        measured%reality = reality; measured%backward_error = backward_error
        measured%finite = finite_count; measured%zeros = zero_count
        measured%positive = wave_count; measured%equilibrium = equilibrium_defect
        measured%real_ok = reality < tolerance
        measured%count_ok = finite_count == zero_count + 2*wave_count
        if (.not. study) then
            call require(measured%real_ok, 'ideal squared frequencies real to roundoff')
            call require(measured%count_ok, 'periodic finite spectrum count')
        end if
        do i = 1, 3
            idx = minloc(abs(positive(:wave_count) - exact(i)), dim=1)
            frequencies(i) = positive(idx)
        end do
        optical_fraction = 0
        if (.not. study) then
            max_energy = 0; optical_vector = 0; optical_energy = 0
            do j = 1, nb
                if (abs(beta(j)) < tolerance*(abs(alpha(j)) + abs(beta(j)))) cycle
                if (abs(lambda(j)) >= tolerance) cycle
                scale = sqrt(sum(abs(vectors(6*ns + 1:6*ns + na, j))**2))
                if (scale < tiny(scale)) cycle
                coefficients = matmul(transform, vectors(:, j))
                ! LAPACK normalizes coefficients, whose magnetic units vary with the mesh.
                coefficients = coefficients/scale
                energy = 0; mean_field = 0; volume = 0; magnetic_energy = 0
                do t = 1, mm%topology%n_triangles
                    do q = 1, op%fs%nq
                        call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                        call require(status == 0, 'null-mode physical magnetic readback')
                        field = matmul(basis(25:27, :), coefficients(op%map(:, t)))
                        energy = energy + weight*sum(abs(field)**2)
                        magnetic_energy = magnetic_energy + weight*abs(field)**2
                        mean_field = mean_field + weight*exp(cmplx(0.0_dp, -x(1), dp))*field
                        volume = volume + weight
                    end do
                end do
                if (energy <= max_energy) cycle
                max_energy = energy
                fraction = sum(abs(mean_field)**2)/(volume*energy)
                optical_fraction = fraction
                optical_vector = vectors(:, j)/scale
                optical_energy = magnetic_energy/energy
            end do
            print '(a,i4,2es20.12)', 'magnetic null energy / smooth Fourier fraction: ', &
                n, max_energy, optical_fraction
            if (audit .and. zero_count > 3 .and. max_energy > 1.0e-12_dp) then
                print '(a,3es20.12)', 'magnetic null Bx / By / Bz energy fractions: ', &
                    optical_energy
                coefficients = matmul(transform, optical_vector)
                call vrmhd_section_residual(op, p, u, zero, real(coefficients, dp), &
                    zero, source, r, mr, status)
                call require(status == 0, 'null-mode real matrix-free Jacobian')
                call vrmhd_section_residual(op, p, u, zero, aimag(coefficients), &
                    zero, source, r, mi, status)
                call require(status == 0, 'null-mode imaginary matrix-free Jacobian')
                projected = bloch_projection(transform,cmplx(mr,mi,qp))
                scale = sqrt(sum(abs(original_k)**2))*sqrt(sum(abs(optical_vector)**2))
                print '(a,es20.12)', 'magnetic null matrix-free relative residual: ', &
                    sqrt(sum(abs(projected)**2))/scale
                call require(sqrt(sum(abs(projected)**2))/scale < tolerance, &
                    'extra null also satisfies production matrix-free Jacobian')
            end if
        end if
        if (audit) call shear_null(op, p, u, source, transform, original_k, ns, measured)
        if (audit) call audit_block(original_k, alpha, beta, lambda, &
            exact, tolerance, zero_count, nulls, unmatched, study, measured)
        if (study) then
            ! The three permitted nulls are removed by subspace projection below.
            call classify_unmatched(op, p, transform, original_k, original_m, &
                vectors, alpha, beta, lambda, exact, tolerance, nulls, measured)
            measured%unmatched = unmatched
            call system_clock(tick1)
            measured%seconds = real(tick1-tick0,dp)/tick_rate
            result = measured
        end if
    contains
        subroutine harmonic_frame(axes)
            !! Point-data B0, grad x and grad y of the first quadrature point;
            !! local returns the orbit gradient of the constant 1.
            real(dp), intent(out) :: axes(3, 3)
            real(dp) :: weight0, x0(2), values(op%nlocal, 2), ones(op%nlocal)
            integer :: l, st
            allocate(orbit(36, op%nlocal))
            call meridian_lagrange_nodes(op%scalar_fs, nodes)
            call vrmhd_fourier_data(fop, 1, 1, x0, weight0, basis, orbit, st)
            call require(st == 0, 'harmonic frame data')
            values = 0; ones = 0
            do l = 1, op%nlocal
                if (op%map(l, 1) > n0) cycle
                values(l, :) = nodes(:, op%map(l, 1)); ones(l) = 1
            end do
            axes(:, 1) = matmul(basis(25:27, :), u(op%map(:, 1)))
            axes(:, 2) = matmul(basis(13:15, :), values(:, 1))
            axes(:, 3) = matmul(basis(13:15, :), values(:, 2))
            local = matmul(orbit(13:15, :), ones)
        end subroutine
    end subroutine

    subroutine shear_null(op, p, u, source, transform, stiffness, ns, measured)
        type(vrmhd_section_t), intent(in) :: op
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: u(:), source(:, :, :)
        complex(dp), intent(in) :: transform(:, :), stiffness(:, :)
        integer, intent(in) :: ns
        type(spectrum_result_t), intent(inout) :: measured
        complex(dp) :: shear(size(stiffness, 1), ns), vh(ns, ns), left(1, 1)
        complex(dp) :: work(16*size(stiffness, 1)), mode(size(stiffness, 1))
        complex(dp) :: coefficients(op%ndof), point(36), projected(size(mode))
        real(dp) :: singular(ns), rwork(5*ns), zero(op%ndof), r(op%ndof)
        real(dp) :: real_jv(op%ndof), imag_jv(op%ndof), basis(36, op%nlocal)
        real(dp) :: x(2), weight, velocity_energy, gradient_energy, scale
        integer :: info, t, q, status
        interface
            subroutine zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, &
                    work, lwork, rwork, info)
                import :: dp
                character, intent(in) :: jobu, jobvt
                integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
                complex(dp), intent(inout) :: a(lda, *)
                real(dp), intent(out) :: s(*), rwork(*)
                complex(dp), intent(out) :: u(ldu, *), vt(ldvt, *), work(*)
                integer, intent(out) :: info
            end subroutine
        end interface
        shear = stiffness(:, 3*ns + 1:4*ns)
        call zgesvd('N', 'A', size(mode), ns, shear, size(mode), singular, left, 1, &
            vh, ns, work, size(work), rwork, info)
        call require(info == 0, 'velocity-only null diagnostic SVD')
        print '(a,*(es20.12))', 'vz-only coupling singular values: ', singular
        measured%shear_ratio = singular(ns)/singular(1)
        measured%shear_nulls = count(singular < 1024*epsilon(1.0_dp)*singular(1))
        mode = 0; mode(3*ns + 1:4*ns) = conjg(vh(ns, :))
        coefficients = matmul(transform, mode); zero = 0
        call vrmhd_section_residual(op, p, u, zero, real(coefficients, dp), zero, &
            source, r, real_jv, status)
        call require(status == 0, 'shear real matrix-free Jacobian')
        call vrmhd_section_residual(op, p, u, zero, aimag(coefficients), zero, &
            source, r, imag_jv, status)
        call require(status == 0, 'shear imaginary matrix-free Jacobian')
        projected = bloch_projection(transform,cmplx(real_jv,imag_jv,qp))
        scale = sqrt(sum(abs(stiffness)**2))*sqrt(sum(abs(mode)**2))
        print '(a,es20.12)', 'vz null production Jacobian relative residual: ', &
            sqrt(sum(abs(projected)**2))/scale
        print '(a,es20.12)', 'vz null relative residual in shear coupling: ', &
            sqrt(sum(abs(projected)**2)/sum(abs(stiffness(:, 3*ns + 1:4*ns))**2))
        print '(a,es20.12)', 'vz null assembled / matrix-free relative difference: ', &
            sqrt(sum(abs(projected - matmul(stiffness, mode))**2))/scale
        velocity_energy = 0; gradient_energy = 0
        do t = 1, op%fs%mm%topology%n_triangles
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                call require(status == 0, 'shear continuum diagnostic basis')
                point = matmul(basis, coefficients(op%map(:, t)))
                velocity_energy = velocity_energy + weight*abs(point(4))**2
                gradient_energy = gradient_energy + weight*abs(point(18))**2
            end do
        end do
        print '(a,es20.12)', 'vz null continuum Alfven derivative / velocity [1/s]: ', &
            abs(background_field(1))*sqrt(gradient_energy/(p%mu0*velocity_energy))
    end subroutine

    integer function nodal_group(x, degree, h, repeat) result(group)
        real(dp), intent(in) :: x(2), h
        integer, intent(in) :: degree, repeat
        group = 1 + modulo(nint(degree*x(1)/h), repeat*degree) &
            + repeat*degree*modulo(nint(degree*x(2)/h), repeat*degree)
    end function

    subroutine equilibrium_potential(x, a)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: a(3)
        ! FortFEM straight-cylinder slots are (x, z, y).
        a = cmplx([-background_field(3)*x(2)/2, &
            background_field(1)*x(2) - background_field(2)*x(1), &
            background_field(3)*x(1)/2], 0.0_dp, dp)
    end subroutine

    subroutine audit_block(stiffness, alpha, beta, lambda, exact, tolerance, &
            zero_count, nulls, unmatched, study, measured)
        integer, intent(in) :: zero_count
        logical, intent(in) :: study
        type(spectrum_result_t), intent(inout) :: measured
        complex(dp), intent(in) :: stiffness(:, :), alpha(:), beta(:), lambda(:)
        complex(dp), intent(in) :: nulls(:, :)
        real(dp), intent(in) :: exact(3), tolerance
        integer, intent(out) :: unmatched
        character(6), parameter :: branches(3) = ['Alfven', 'fast  ', 'slow  ']
        complex(dp) :: defect(size(lambda))
        real(dp) :: distances(3), separation, error
        integer :: j, branch, sign_index, counts(3, 2), zeros
        ! At nonzero k there are two pressure-balanced entropy modes and one gauge.
        do j = 1, 3
            defect = matmul(stiffness, nulls(:, j))
            error = sqrt(sum(abs(defect)**2)) &
                /(sqrt(sum(abs(stiffness)**2))*sqrt(sum(abs(nulls(:, j))**2)))
            print '(a,i2,es20.12)', 'analytic entropy / entropy / gauge null defect: ', &
                j, error
            measured%analytic_null_defect = max(measured%analytic_null_defect, error)
            call require(error < tolerance, 'expected continuum null modes')
        end do
        print '(a)', 'complete computed Bloch spectrum: k=(1,0,0) [1/m]'
        print '(a,es20.12)', 'resolved oracle frequency cutoff [1/s]: ', maxval(exact)
        print '(a)', 'expected nulls: entropy=2, gauge=1, translation=0 (k is nonzero)'
        counts = 0; zeros = 0; unmatched = max(0, zero_count - 3)
        do j = 1, size(lambda)
            if (abs(beta(j)) < tolerance*(abs(alpha(j)) + abs(beta(j)))) then
                print '(a,i3,a)', 'eigenvalue ', j, ': infinite, algebraic current constraint'
                cycle
            end if
            if (abs(lambda(j)) < tolerance) then
                zeros = zeros + 1
                print '(a,i3,a,2es22.14)', 'eigenvalue ', j, &
                    ': zero subspace, omega real / imaginary [1/s] ', &
                    aimag(lambda(j)), -real(lambda(j), dp)
                cycle
            end if
            distances = abs(abs(aimag(lambda(j))) - exact)
            branch = minloc(distances, dim=1)
            separation = minval(abs(exact(branch) - pack(exact, &
                [(sign_index /= branch, sign_index=1, 3)])))/2
            sign_index = 1
            if (aimag(lambda(j)) > 0) sign_index = 2
            print '(a,i3,3a,3es22.14)', 'eigenvalue ', j, ': ', branches(branch), &
                ', omega real / imaginary, oracle [1/s] ', aimag(lambda(j)), &
                -real(lambda(j), dp), sign(1.0_dp, aimag(lambda(j)))*exact(branch)
            if (distances(branch) < separation) then
                counts(branch, sign_index) = counts(branch, sign_index) + 1
            else if (abs(lambda(j)) < maxval(exact)) then
                unmatched = unmatched + 1
            end if
        end do
        print '(a,3i4)', 'computed / expected / unmatched zero multiplicity: ', &
            zeros, 3, max(0, zero_count - 3)
        print '(a,i4)', 'unmatched eigenvalues below resolved cutoff: ', unmatched
        measured%signed_counts = counts
        if (.not. study) call require(all(counts == 1), 'each signed oracle wave occurs once')
    end subroutine
    subroutine classify_unmatched(op, p, transform, stiffness, mass, vectors, &
            alpha, beta, lambda, exact, tolerance, nulls, measured)
        type(vrmhd_section_t), intent(in) :: op
        type(vrmhd_parameters_t), intent(in) :: p
        complex(dp), intent(in) :: transform(:, :), stiffness(:, :), mass(:, :)
        complex(dp), intent(in) :: vectors(:, :), alpha(:), beta(:), lambda(:)
        complex(dp), intent(in) :: nulls(:, :)
        real(dp), intent(in) :: exact(3), tolerance
        type(spectrum_result_t), intent(inout) :: measured
        type(spectrum_result_t) :: static
        complex(dp), allocatable :: a(:, :), vh(:, :), modes(:, :), eq(:, :)
        complex(dp), allocatable :: ek(:, :), eb(:, :), total(:, :), energetic(:, :)
        complex(dp), allocatable :: whitened(:, :), work(:), coefficients(:)
        complex(dp), allocatable :: curls(:, :), curl_mass(:, :)
        complex(dp) :: known(size(lambda), 3), left(1, 1), v(size(lambda))
        real(dp) :: singular(size(lambda)), rw(8*size(lambda)), scales(size(lambda))
        real(dp) :: energies(size(lambda)), ratio(size(lambda)), distance(3), separation
        real(dp) :: magnetic, kinetic, physical, norm, energy_floor, reference_energy
        integer :: nb, ns, na, j, i, k, d, energetic_count, info, zeros
        interface
            subroutine zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, &
                    work, lwork, rwork, info)
                import :: dp
                character, intent(in) :: jobu, jobvt
                integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
                complex(dp), intent(inout) :: a(lda, *)
                real(dp), intent(out) :: s(*), rwork(*)
                complex(dp), intent(out) :: u(ldu, *), vt(ldvt, *), work(*)
                integer, intent(out) :: info
            end subroutine
            subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
                import :: dp
                character, intent(in) :: jobz, uplo
                integer, intent(in) :: n, lda, lwork
                complex(dp), intent(inout) :: a(lda, *)
                real(dp), intent(out) :: w(*), rwork(*)
                complex(dp), intent(out) :: work(*)
                integer, intent(out) :: info
            end subroutine
        end interface
        nb = size(lambda)
        ns = 0
        do j = 1, nb
            if (any(abs(transform(1:op%scalar_fs%n0, j)) > 0)) ns = ns + 1
        end do
        na = (nb - 6*ns)/2
        allocate(a(nb,nb), vh(nb,nb), work(32*nb))
        scales = 1
        scales(6*ns+1:6*ns+na) = sqrt(p%mu0)
        scales(6*ns+na+1:) = 1/sqrt(p%mu0)
        do j = 1, nb
            a(:,j) = stiffness(:,j)*scales(j)
        end do
        do i = 1, nb
            norm = sqrt(sum(abs(a(i,:))**2))
            if (norm > tiny(norm)) a(i,:) = a(i,:)/norm
        end do
        call zgesvd('N', 'A', nb, nb, a, nb, singular, left, 1, vh, nb, &
            work, size(work), rw, info)
        call require(info == 0, 'static coupling nullspace SVD')
        d = count(singular < tolerance*singular(1))
        measured%svd_zeros = d
        measured%null_split_ok = d == measured%zeros
        call require(d >= 3, 'three analytic nulls in static coupling')
        known = nulls
        do j = 1, 3
            known(:,j) = known(:,j)/scales
            do i = 1, j-1
                known(:,j) = known(:,j)-known(:,i)*dot_product(known(:,i), known(:,j))
            end do
            known(:,j) = known(:,j)/sqrt(sum(abs(known(:,j))**2))
        end do
        allocate(modes(nb,d), eq(nb,d))
        modes = conjg(transpose(vh(nb-d+1:nb,:)))
        do j = 1, d
            do i = 1, 3
                modes(:,j) = modes(:,j)-known(:,i)*dot_product(known(:,i), modes(:,j))
            end do
        end do
        ! Rank after quotienting the allowed entropy/entropy/gauge subspace.
        call zgesvd('S', 'N', nb, d, modes, nb, singular, eq, nb, left, 1, &
            work, size(work), rw, info)
        call require(info == 0, 'allowed-null quotient SVD')
        ! The three verified analytic nulls fix the quotient dimension exactly.
        d = measured%svd_zeros-3
        print '(a,es20.12)', 'allowed-null projection subspace defect: ', &
            maxval(singular(d+1:measured%svd_zeros))
        call require(maxval(singular(d+1:measured%svd_zeros)) < sqrt(tolerance), &
            'analytic nulls align with the static SVD subspace')
        if (d > 0) call require(singular(d) > 0.5_dp, 'resolved static null quotient')
        zeros = max(0, measured%zeros-3)
        if (d > 0) then
            modes = eq
            do j = 1, size(modes,2)
                modes(:,j) = modes(:,j)*scales
            end do
            allocate(ek(d,d), eb(d,d), total(d,d))
            call energy_grams(modes(:,:d), ek, eb)
            total = ek+eb
            call zheev('V','U',d,total,d,energies,work,size(work),rw,info)
            call require(info == 0, 'null kinetic plus magnetic energy metric')
            reference_energy = maxval(abs(mass(ns+1:4*ns,ns+1:4*ns)))
            energy_floor = 1.0e-10_dp*max(maxval(energies(:d)),reference_energy)
            energetic_count = count(energies(:d) > energy_floor)
            measured%energetic_nulls = energetic_count
            if (energetic_count > 0) then
                allocate(energetic(d,energetic_count), whitened(energetic_count,energetic_count))
                k = 0
                do j = 1, d
                    if (energies(j) <= energy_floor) cycle
                    k = k+1
                    energetic(:,k) = total(:,j)/sqrt(energies(j))
                end do
                whitened = matmul(transpose(conjg(energetic)),matmul(eb,energetic))
                call zheev('N','U',energetic_count,whitened,energetic_count, &
                    ratio,work,size(work),rw,info)
                call require(info == 0, 'magnetic energy fractions on null quotient')
                do j = 1, energetic_count
                    call classify(ratio(j), static)
                end do
            end if
            static%mixed = static%mixed+max(0,d-energetic_count)
        end if
        measured%static_magnetic = static%magnetic
        measured%static_velocity = static%velocity
        measured%static_mixed = static%mixed
        if (measured%null_split_ok) then
            call require(static%magnetic+static%velocity+static%mixed == zeros, &
                'static quotient dimension agrees with QZ null count')
            measured%magnetic = static%magnetic
            measured%velocity = static%velocity
            measured%mixed = static%mixed
        else
            ! Different numerical null counts cannot support a common energy split.
            measured%mixed = zeros
        end if
        do j = 1, nb
            if (abs(beta(j)) < tolerance*(abs(alpha(j))+abs(beta(j)))) cycle
            if (abs(lambda(j)) < tolerance) cycle
            distance = abs(abs(aimag(lambda(j)))-exact)
            i = minloc(distance,dim=1)
            separation = minval(abs(exact(i)-pack(exact,[(k /= i,k=1,3)])))/2
            if (distance(i) < separation .or. abs(lambda(j)) >= maxval(exact)) cycle
            v = vectors(:,j)
            kinetic = real(dot_product(v(ns+1:4*ns), &
                matmul(mass(ns+1:4*ns,ns+1:4*ns),v(ns+1:4*ns))),dp)
            coefficients = matmul(transform(:,6*ns+1:6*ns+na),v(6*ns+1:6*ns+na))
            block
                complex(dp) :: curl(op%fs%sp(2)%ndof), image(op%fs%sp(2)%ndof)
                curl = cmplx(csc_matvec(op%fs%curl%a,real(coefficients(6*op%scalar_fs%n0+1: &
                    6*op%scalar_fs%n0+op%fs%sp(1)%ndof),dp)), &
                    csc_matvec(op%fs%curl%a,aimag(coefficients(6*op%scalar_fs%n0+1: &
                    6*op%scalar_fs%n0+op%fs%sp(1)%ndof))),dp)
                image = cmplx(csc_matvec(op%fs%mass(2),real(curl,dp)), &
                    csc_matvec(op%fs%mass(2),aimag(curl)),dp)
                magnetic = real(dot_product(curl,image),dp)/p%mu0
            end block
            physical = kinetic+magnetic
            energy_floor = 1.0e-10_dp*maxval(abs(mass(ns+1:4*ns,ns+1:4*ns))) &
                *sum(abs(v/scales)**2)
            if (physical > energy_floor) then
                call classify(magnetic/physical,measured)
            else
                measured%mixed = measured%mixed+1
            end if
        end do
        print '(a,5i5)', 'SVD nulls / energetic / magnetic / velocity / mixed: ', &
            measured%svd_zeros, measured%energetic_nulls, measured%magnetic, &
            measured%velocity, measured%mixed
    contains
        subroutine energy_grams(m, kinetic_gram, magnetic_gram)
            complex(dp), intent(in) :: m(:,:)
            complex(dp), intent(out) :: kinetic_gram(:,:), magnetic_gram(:,:)
            integer :: j, a0, j0
            complex(dp) :: full(op%ndof), curl(op%fs%sp(2)%ndof)
            complex(dp) :: image(op%fs%sp(2)%ndof)
            a0 = 6*op%scalar_fs%n0; j0 = a0+op%fs%sp(1)%ndof
            allocate(curls(size(curl),size(m,2)),curl_mass(size(curl),size(m,2)))
            kinetic_gram = matmul(transpose(conjg(m(ns+1:4*ns,:))), &
                matmul(mass(ns+1:4*ns,ns+1:4*ns),m(ns+1:4*ns,:)))
            do j = 1, size(m,2)
                full = matmul(transform,m(:,j))
                curl = cmplx(csc_matvec(op%fs%curl%a,real(full(a0+1:j0),dp)), &
                    csc_matvec(op%fs%curl%a,aimag(full(a0+1:j0))),dp)
                image = cmplx(csc_matvec(op%fs%mass(2),real(curl,dp)), &
                    csc_matvec(op%fs%mass(2),aimag(curl)),dp)
                curls(:,j) = curl; curl_mass(:,j) = image
            end do
            magnetic_gram = matmul(transpose(conjg(curls)),curl_mass)/p%mu0
        end subroutine
    end subroutine

    subroutine classify(fraction, measured)
        real(dp), intent(in) :: fraction
        type(spectrum_result_t), intent(inout) :: measured
        if (fraction >= 0.9_dp) then
            measured%magnetic = measured%magnetic+1
        else if (fraction <= 0.1_dp) then
            measured%velocity = measured%velocity+1
        else
            measured%mixed = measured%mixed+1
        end if
    end subroutine
end module vrmhd_spectrum_support
