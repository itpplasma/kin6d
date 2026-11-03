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
    real(dp) :: frequencies(3, 3), errors(3, 3), rates(3), exact(3), b0(3)
    real(dp) :: optical(3)
    integer :: level, status
    b0 = sqrt(oracle_mu0)*[0.6_dp, 0.0_dp, 0.8_dp]
    call uniform_waves([1.0_dp, 0.0_dp, 0.0_dp], b0, 1.0_dp, 0.6_dp, &
        5.0_dp/3, exact, status)
    call require(status == oracle_ok, 'existing uniform-wave oracle')
    exact = sqrt(exact)
    print '(a,3es20.12)', 'oracle Alfven / fast / slow angular frequencies [1/s]: ', exact
    do level = 1, 3
        call spectrum(2**(level + 1), b0, frequencies(:, level), optical(level))
        errors(:, level) = abs(frequencies(:, level)/exact - 1)
        print '(a,i4,3es20.12)', 'wave frequencies [1/s], cells per period: ', &
            2**(level + 1), frequencies(:, level)
        print '(a,i4,3es20.12)', 'relative frequency errors: ', &
            2**(level + 1), errors(:, level)
        if (level == 1) cycle
        rates = log(errors(:, level - 1)/errors(:, level))/log(2.0_dp)
        print '(a,3f12.6)', 'Alfven / fast / slow frequency orders: ', rates
        call require(minval(rates) > 1.8_dp, 'P1 frequency design order')
        call require(optical(level) < optical(level - 1), &
            'cell-scale null mode must lose smooth Fourier content')
    end do
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (ok) return
        print '(a)', message
        error stop 'RM2 linear spectrum gate'
    end subroutine

    subroutine spectrum(n, b0, frequencies, optical_fraction)
        integer, intent(in) :: n
        real(dp), intent(in) :: b0(3)
        real(dp), intent(out) :: frequencies(3), optical_fraction
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
        complex(dp) :: phase, field(3), mean_field(3)
        real(dp) :: rwork(8*nb), length, x(2), start(2), tangent(2), midpoint(2)
        real(dp) :: weight, direction, volume, energy, max_energy, fraction, scale
        real(dp) :: positive(3), tolerance, reality, equilibrium_defect
        real(dp) :: backward_error, diagnostic_reality, column_scale(nb), row_scale(nb)
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
        call zggev('N', 'V', nb, stiffness, nb, mass, nb, alpha, beta, dummy, 1, &
            vectors, nb, work, size(work), rwork, info)
        call require(info == 0, 'assembled generalized eigenproblem')
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
        column_scale = 1; row_scale = 1
        column_scale(7:10) = sqrt(p%mu0)
        column_scale(11:14) = 1/sqrt(p%mu0)
        row_scale(7:14) = 1/sqrt(p%mu0)
        do j = 1, nb
            stiffness(:, j) = row_scale*original_k(:, j)*column_scale(j)
            mass(:, j) = row_scale*original_m(:, j)*column_scale(j)
        end do
        call zggev('N', 'N', nb, stiffness, nb, mass, nb, diagnostic_alpha, &
            diagnostic_beta, dummy, 1, dummy, 1, work, size(work), rwork, info)
        call require(info == 0, 'diagnostic change of magnetic units')
        diagnostic_reality = 0
        do j = 1, nb
            if (abs(diagnostic_beta(j)) < tolerance &
                *(abs(diagnostic_alpha(j)) + abs(diagnostic_beta(j)))) cycle
            diagnostic_reality = max(diagnostic_reality, &
                abs(aimag(-(diagnostic_alpha(j)/diagnostic_beta(j))**2)))
        end do
        print '(a,es20.12)', 'diagnostic magnetic-unit-scaled omega^2 imaginary defect: ', &
            diagnostic_reality
        if (wave_count == 3) print '(a,3es20.12)', &
            'positive angular frequencies [1/s], unsorted: ', positive
        call require(reality < tolerance, 'ideal squared frequencies real to roundoff')
        call require(finite_count == 10 .and. zero_count == 4 .and. wave_count == 3, &
            'periodic finite spectrum count')
        do i = 1, 2
            do j = i + 1, 3
                if (positive(i) <= positive(j)) cycle
                weight = positive(i); positive(i) = positive(j); positive(j) = weight
            end do
        end do
        frequencies = positive([2, 3, 1])
        max_energy = 0; optical_fraction = 0
        do j = 1, nb
            if (abs(beta(j)) < tolerance*(abs(alpha(j)) + abs(beta(j)))) cycle
            if (abs(lambda(j)) >= tolerance) cycle
            coefficients = matmul(transform, vectors(:, j))
            energy = 0; mean_field = 0; volume = 0
            do t = 1, mm%topology%n_triangles
                do q = 1, op%fs%nq
                    call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                    call require(status == 0, 'null-mode physical magnetic readback')
                    field = matmul(basis(25:27, :), coefficients(op%map(:, t)))
                    energy = energy + weight*sum(abs(field)**2)
                    mean_field = mean_field + weight*exp(cmplx(0.0_dp, -x(1), dp))*field
                    volume = volume + weight
                end do
            end do
            if (energy <= max_energy) cycle
            max_energy = energy
            fraction = sum(abs(mean_field)**2)/(volume*energy)
            optical_fraction = fraction
        end do
        print '(a,i4,2es20.12)', 'magnetic null energy / smooth Fourier fraction: ', &
            n, max_energy, optical_fraction
        call require(max_energy > 1.0e-12_dp, 'identify magnetic optical null branch')
        ! A resolved null branch would retain most energy in its smooth Fourier field.
        call require(optical_fraction < 0.5_dp, 'no resolved spurious magnetic null mode')
    end subroutine
end program test_vrmhd_spectra
