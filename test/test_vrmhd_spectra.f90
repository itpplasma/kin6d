program test_vrmhd_spectra
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_gs_symmetry, only: gs_symmetry_t, GS_TRANSLATION, GS_CHART_CARTESIAN
    use kin6d_vrmhd_parameters, only: vrmhd_parameters_t
    use kin6d_vrmhd_section, only: vrmhd_section_t, vrmhd_section_build, &
        vrmhd_section_basis, vrmhd_section_residual, vrmhd_section_jacobian, &
        vrmhd_section_mass
    use rm_mhd_oracles, only: uniform_waves, oracle_mu0, oracle_ok
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    use fortsparse, only: csc_t
    implicit none
    integer, parameter :: levels = 5
    real(dp) :: frequencies(3, levels), errors(3, levels), rates(3), exact(3), b0(3)
    real(dp) :: optical(levels)
    integer :: level, status, unmatched
    b0 = sqrt(oracle_mu0)*[0.6_dp, 0.0_dp, 0.8_dp]
    call uniform_waves([1.0_dp, 0.0_dp, 0.0_dp], b0, 1.0_dp, 0.6_dp, &
        5.0_dp/3, exact, status)
    call require(status == oracle_ok, 'existing uniform-wave oracle')
    exact = sqrt(exact)
    print '(a,3es20.12)', 'oracle Alfven / fast / slow angular frequencies [1/s]: ', exact
    do level = 1, levels
        call spectrum(2**(level + 1), b0, exact, level == levels, &
            frequencies(:, level), optical(level), unmatched)
        errors(:, level) = abs(frequencies(:, level)/exact - 1)
        print '(a,i4,3es20.12)', 'wave frequencies [1/s], cells per period: ', &
            2**(level + 1), frequencies(:, level)
        print '(a,i4,3es20.12)', 'relative frequency errors: ', &
            2**(level + 1), errors(:, level)
        if (level == 1) cycle
        rates = log(errors(:, level - 1)/errors(:, level))/log(2.0_dp)
        print '(a,3f12.6)', 'Alfven / fast / slow frequency orders: ', rates
        if (level == levels) &
            call require(minval(rates) > 1.8_dp, 'P1 frequency design order')
    end do
    call require(unmatched == 0, 'no unmatched eigenvalue below resolved cutoff')
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (ok) return
        print '(a)', message
        error stop 'RM2 linear spectrum gate'
    end subroutine

    subroutine spectrum(n, b0, exact, audit, frequencies, optical_fraction, unmatched)
        integer, intent(in) :: n
        real(dp), intent(in) :: b0(3), exact(3)
        logical, intent(in) :: audit
        real(dp), intent(out) :: frequencies(3), optical_fraction
        integer, intent(out) :: unmatched
        integer, parameter :: nb = 14
        type(mesh_2d_t) :: mesh
        type(meridian_mesh_t) :: mm
        type(vrmhd_section_t) :: op
        type(vrmhd_parameters_t) :: p
        type(csc_t) :: jac
        real(dp), allocatable :: u(:), zero(:), source(:, :, :), r(:), unused(:)
        real(dp), allocatable :: mr(:), mi(:), basis(:, :)
        complex(dp), allocatable :: transform(:, :), image(:), coefficients(:)
        complex(dp) :: stiffness(nb, nb), mass(nb, nb), alpha(nb), beta(nb)
        complex(dp) :: vectors(nb, nb), dummy(1, 1), work(16*nb), lambda(nb)
        complex(dp) :: original_k(nb, nb), original_m(nb, nb), residual(nb)
        complex(dp) :: diagnostic_alpha(nb), diagnostic_beta(nb)
        complex(dp) :: optical_vector(nb), projected(nb)
        complex(dp) :: phase, field(3), mean_field(3)
        real(dp) :: rwork(8*nb), length, x(2), start(2), tangent(2), midpoint(2)
        real(dp) :: weight, direction, volume, energy, max_energy, fraction, scale
        real(dp) :: positive(3), tolerance, reality, equilibrium_defect
        real(dp) :: backward_error, diagnostic_reality, column_scale(nb), row_scale(nb)
        real(dp) :: magnetic_energy(3), optical_energy(3)
        integer :: n0, nn, n1, a0, j0, t, i, j, q, row, entry, edge, group
        integer :: status, info, finite_count, zero_count, wave_count
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
        length = 2*acos(-1.0_dp)
        unmatched = 0
        call mesh%create_rectangular(n + 1, n + 1, 0.0_dp, length, 0.0_dp, length)
        call meridian_mesh_straight(mesh, 1, mm, status)
        call require(status == 0, 'periodic wave mesh')
        call vrmhd_section_build(mm, &
            gs_symmetry_t(GS_TRANSLATION, GS_CHART_CARTESIAN, 0.0_dp), 1, 1, op, status)
        call require(status == 0, 'periodic wave section')
        n0 = op%fs%n0; nn = op%fs%nn; n1 = op%fs%sp(1)%ndof
        a0 = 6*n0; j0 = a0 + n1
        p%mi = 1; p%mu0 = oracle_mu0
        allocate (u(op%ndof), zero(op%ndof), r(op%ndof), unused(op%ndof), &
            source(6, op%fs%nq, mm%topology%n_triangles), mr(op%ndof), mi(op%ndof), &
            transform(op%ndof, nb), image(op%ndof), coefficients(op%ndof), &
            basis(36, op%nlocal))
        u = 0; zero = 0; source = 0; transform = 0
        u(:n0) = 1; u(4*n0 + 1:6*n0) = 0.3_dp
        do i = 1, n0
            x = mm%topology%vertices(:, i)
            phase = exp(cmplx(0.0_dp, x(1), dp))
            do group = 1, 6
                transform((group - 1)*n0 + i, group) = phase
            end do
            u(a0 + nn + i) = b0(1)*x(2)
            transform(a0 + nn + i, 10) = phase
            transform(j0 + nn + i, 14) = phase
        end do
        do t = 1, mm%topology%n_triangles
            do i = 1, 3
                j = mod(i, 3) + 1
                start = mm%topology%vertices(:, mm%topology%triangles(i, t))
                tangent = mm%topology%vertices(:, mm%topology%triangles(j, t)) - start
                midpoint = start + tangent/2
                edge = op%fs%sp(1)%dof(i, t)
                direction = real(op%fs%sp(1)%sgn(i, t), dp)
                u(a0 + edge) = direction*b0(3) &
                    *(-midpoint(2)*tangent(1) + midpoint(1)*tangent(2))/2
                group = 3
                if (abs(tangent(2)) < epsilon(length)) group = 1
                if (abs(tangent(1)) < epsilon(length)) group = 2
                if (group == 2) then
                    direction = direction*sign(1.0_dp, tangent(2))
                else
                    direction = direction*sign(1.0_dp, tangent(1))
                end if
                phase = direction*norm2(tangent)*exp(cmplx(0.0_dp, midpoint(1), dp))
                transform(a0 + edge, 6 + group) = phase
                transform(j0 + edge, 10 + group) = phase
            end do
        end do
        ! Translation invariance makes each Bloch block an invariant periodic subspace.
        call vrmhd_section_residual(op, p, u, zero, zero, zero, source, r, unused, status)
        call require(status == 0, 'uniform equilibrium residual')
        equilibrium_defect = maxval(abs(matmul(transpose(conjg(transform)), r)))
        print '(a,i4,es20.12)', 'periodic equilibrium residual: ', n, equilibrium_defect
        call require(equilibrium_defect < 1.0e-12_dp, 'uniform periodic equilibrium')
        call vrmhd_section_jacobian(op, p, u, zero, source, 0.0_dp, jac, status)
        call require(status == 0, 'assembled production linearization')
        do j = 1, nb
            image = 0
            do i = 1, op%ndof
                do entry = jac%col_ptr(i), jac%col_ptr(i + 1) - 1
                    row = jac%row_idx(entry)
                    image(row) = image(row) + jac%val(entry)*transform(i, j)
                end do
            end do
            stiffness(:, j) = matmul(transpose(conjg(transform)), image)
            call vrmhd_section_mass(op, p, u, real(transform(:, j), dp), mr, status)
            call require(status == 0, 'production real mass action')
            call vrmhd_section_mass(op, p, u, aimag(transform(:, j)), mi, status)
            call require(status == 0, 'production imaginary mass action')
            image = cmplx(mr, mi, dp)
            mass(:, j) = matmul(transpose(conjg(transform)), image)
        end do
        original_k = stiffness; original_m = mass
        ! Balance magnetic units before QZ; this is an invertible pencil equivalence.
        column_scale = 1; row_scale = 1
        column_scale(7:10) = sqrt(p%mu0)
        column_scale(11:14) = 1/sqrt(p%mu0)
        row_scale(7:14) = 1/sqrt(p%mu0)
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
            call require(wave_count <= 3, 'no extra oscillatory branch')
            positive(wave_count) = aimag(lambda(j))
        end do
        print '(a,i4,3i4,es20.12)', 'finite / zero / positive modes, Im(omega^2): ', &
            n, finite_count, zero_count, wave_count, reality
        print '(a,es20.12)', 'maximum generalized eigenpair backward error: ', backward_error
        stiffness = original_k; mass = original_m
        call zggev('N', 'N', nb, stiffness, nb, mass, nb, diagnostic_alpha, &
            diagnostic_beta, dummy, 1, dummy, 1, work, size(work), rwork, info)
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
        if (wave_count == 3) print '(a,3es20.12)', &
            'positive angular frequencies [1/s], unsorted: ', positive
        call require(reality < tolerance, 'ideal squared frequencies real to roundoff')
        call require(finite_count == zero_count + 2*wave_count .and. wave_count == 3, &
            'periodic finite spectrum count')
        do i = 1, 2
            do j = i + 1, 3
                if (positive(i) <= positive(j)) cycle
                weight = positive(i); positive(i) = positive(j); positive(j) = weight
            end do
        end do
        frequencies = positive([2, 3, 1])
        max_energy = 0; optical_fraction = 0; optical_vector = 0; optical_energy = 0
        do j = 1, nb
            if (abs(beta(j)) < tolerance*(abs(alpha(j)) + abs(beta(j)))) cycle
            if (abs(lambda(j)) >= tolerance) cycle
            scale = sqrt(sum(abs(vectors(7:10, j))**2))
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
        if (audit .and. zero_count > 3) then
            call require(max_energy > 1.0e-12_dp, 'identify extra magnetic null branch')
            print '(a,3es20.12)', 'magnetic null Bx / By / Bz energy fractions: ', &
                optical_energy
            coefficients = matmul(transform, optical_vector)
            call vrmhd_section_residual(op, p, u, zero, real(coefficients, dp), &
                zero, source, r, mr, status)
            call require(status == 0, 'null-mode real matrix-free Jacobian')
            call vrmhd_section_residual(op, p, u, zero, aimag(coefficients), &
                zero, source, r, mi, status)
            call require(status == 0, 'null-mode imaginary matrix-free Jacobian')
            projected = matmul(transpose(conjg(transform)), cmplx(mr, mi, dp))
            scale = sqrt(sum(abs(original_k)**2))*sqrt(sum(abs(optical_vector)**2))
            print '(a,es20.12)', 'magnetic null matrix-free relative residual: ', &
                sqrt(sum(abs(projected)**2))/scale
            call require(sqrt(sum(abs(projected)**2))/scale < tolerance, &
                'extra null also satisfies production matrix-free Jacobian')
        end if
        if (audit) call audit_block(n, original_k, alpha, beta, lambda, &
            exact, tolerance, zero_count, unmatched)
    end subroutine

    subroutine audit_block(n, stiffness, alpha, beta, lambda, exact, tolerance, &
            zero_count, unmatched)
        integer, intent(in) :: n, zero_count
        complex(dp), intent(in) :: stiffness(14, 14), alpha(14), beta(14), lambda(14)
        real(dp), intent(in) :: exact(3), tolerance
        integer, intent(out) :: unmatched
        character(6), parameter :: branches(3) = ['Alfven', 'fast  ', 'slow  ']
        complex(dp) :: nulls(14, 3), defect(14)
        real(dp) :: distances(3), separation, error, h
        integer :: j, branch, sign_index, counts(3, 2), zeros
        ! At nonzero k there are two pressure-balanced entropy modes and one gauge.
        nulls = 0
        nulls(1, 1) = 1; nulls(5, 1) = -0.6_dp
        nulls(5, 2) = 1; nulls(6, 2) = -1
        h = 2*acos(-1.0_dp)/n
        nulls(7, 3) = cmplx(0.0_dp, 2*sin(h/2)/h, dp)
        nulls(9, 3) = nulls(7, 3)/sqrt(2.0_dp)
        do j = 1, 3
            defect = matmul(stiffness, nulls(:, j))
            error = sqrt(sum(abs(defect)**2)) &
                /(sqrt(sum(abs(stiffness)**2))*sqrt(sum(abs(nulls(:, j))**2)))
            print '(a,i2,es20.12)', 'analytic entropy / entropy / gauge null defect: ', &
                j, error
            call require(error < tolerance, 'expected continuum null modes')
        end do
        print '(a)', 'complete computed Bloch spectrum: k=(1,0,0) [1/m]'
        print '(a,es20.12)', 'resolved oracle frequency cutoff [1/s]: ', maxval(exact)
        print '(a)', 'expected nulls: entropy=2, gauge=1, translation=0 (k is nonzero)'
        counts = 0; zeros = 0; unmatched = max(0, zero_count - 3)
        do j = 1, 14
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
        call require(all(counts == 1), 'each signed oracle wave occurs once')
    end subroutine
end program test_vrmhd_spectra
