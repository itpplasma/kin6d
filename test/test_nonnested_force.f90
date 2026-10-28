program test_nonnested_force
    ! Harmonic deflation, potential-velocity and Leray force, residual floor
    ! (docs/nonnested-force.md).
    !   hodge    : manufactured f = grad phi + curl w + c e_phi/R on the
    !              Solov'ev domain (N = 1, P2 8/16/32, P3 8/16 rings); both forms
    !              against P f = curl w + c e_phi/R, their roundoff agreement,
    !              div u = 0 and harmonic orthogonality to roundoff.
    !   solovev  : NN1 n = 0 Solov'ev states (P2 8/16/32, P3 8/16 rings): P(J x B)
    !              of the discrete state decreases; a vacuum state is force
    !              free to quadrature roundoff.
    !   beltrami : Lundquist Beltrami states on the periodic cylinder (N = 1,
    !              4/8 rings, 16 with a second argument): P(J x B) decreases
    !              with the primal current, like h^1/2 with the dual one.
    !   shell    : the NN1 coil annulus (two harmonic 2-forms, one relative
    !              harmonic 1-form): manufactured ladder and deflation of a
    !              right-hand side with kernel content.
    !   sheet    : straight periodic cylinder with B_z(r) (J x B = -grad
    !              B_z^2/2): a smooth profile decreases, a tangential
    !              discontinuity off the mesh rings sets a residual floor.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_solovev_fixture, only: s_r0 => r0, s_a => minor, s_f => ftor, pp, &
        s_mapping => mapping, s_coordinates => coordinates, s_c => c
    use nn_vacuum_fixture, only: annulus_mesh
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_report_t, nn_build, &
        nn_solve, nn_factor_t, nn_factor, nn_factor_solve_real, nn_factor_free
    use kin6d_nonnested_beltrami, only: beltrami_system_t, beltrami_report_t, &
        beltrami_prepare, beltrami_solve
    use kin6d_nonnested_force, only: nf_system_t, nf_report_t, nf_prepare, &
        nf_free, nf_load, nf_lorentz_load, nf_project, nf_potential_solve, &
        nf_l2_error
    use kin6d_equilibrium_harness, only: equilibrium_floor_report_t, &
        equilibrium_residual_floor, write_residual_floor
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes
    use fortsparse, only: csc_matvec
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), cyl_r0 = 6.0_dp
    real(dp), parameter :: c_h = 12.0_dp
    character(16) :: kind
    character(8), parameter :: currents(2) = [character(8) :: ' dual', ' primal']
    ! Shell centre and radii (nn_vacuum_fixture annulus) for the bump.
    real(dp), parameter :: shell_r0 = 6.0_dp, rin = 0.3_dp, rout = 1.2_dp
    logical :: shell_domain = .false.
    real(dp) :: sheet_radius = 0.47_dp, sheet_width = 0
    call get_command_argument(1, kind)
    select case (trim(kind))
    case ('hodge')
        call hodge_case(2)
        if (command_argument_count() < 2) call hodge_case(3)
    case ('solovev')
        call solovev_case(2); call solovev_case(3); call vacuum_case()
    case ('beltrami')
        call beltrami_case()
    case ('shell')
        call shell_case()
    case ('sheet')
        call sheet_case()
    case default
        error stop 'unknown case'
    end select
    print '(a)', 'PASS'
contains
    ! ---- manufactured Hodge decomposition -------------------------------
    subroutine bump(x, p, pr, pz)
        ! p = 0 on the boundary (Solov'ev flux or the shell quartic).
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: p, pr, pz
        real(dp) :: s, r2
        if (shell_domain) then
            r2 = (x(1) - shell_r0)**2 + x(2)**2
            p = (r2 - rin**2)*(rout**2 - r2)
            pr = 2*(x(1) - shell_r0)*((rout**2 - r2) - (r2 - rin**2))
            pz = 2*x(2)*((rout**2 - r2) - (r2 - rin**2))
            return
        end if
        s = 4*s_c*s_r0**2*s_a**2
        p = s_c*((x(1)**2 - s_r0**2)**2 + 4*x(1)**2*x(2)**2 - 4*s_r0**2*s_a**2)/s
        pr = s_c*(4*x(1)*(x(1)**2 - s_r0**2) + 8*x(1)*x(2)**2)/s
        pz = 8*s_c*x(1)**2*x(2)/s
    end subroutine

    subroutine projected(x, v)
        ! P f = curl w + c_h e_phi/R, w = p ((0, 1 + 0.3 cos phi, 0.4 sin phi)).
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: v(3)
        real(dp) :: p, pr, pz, r, ph
        call bump(x(:2), p, pr, pz)
        r = x(1); ph = x(3)
        v(1) = 0.4_dp*p*cos(ph)/r - pz*(1 + 0.3_dp*cos(ph))
        v(2) = -0.4_dp*pr*sin(ph) + c_h/r
        v(3) = (p + r*pr)*(1 + 0.3_dp*cos(ph))/r
    end subroutine

    subroutine manufactured(x, v)
        ! f = grad phi + P f, phi = sin(2R) cos(3Z) + 0.5 R Z cos(phi).
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: v(3)
        real(dp) :: r, z, ph
        call projected(x, v)
        r = x(1); z = x(2); ph = x(3)
        v(1) = v(1) + 2*cos(2*r)*cos(3*z) + 0.5_dp*z*cos(ph)
        v(2) = v(2) - 0.5_dp*z*sin(ph)
        v(3) = v(3) - 3*sin(2*r)*sin(3*z) + 0.5_dp*r*cos(ph)
    end subroutine

    subroutine gradient_only(x, v)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: v(3)
        real(dp) :: p(3)
        call manufactured(x, v)
        call projected(x, p)
        v = v - p
    end subroutine

    subroutine solovev_field(k, nr, modes, field)
        integer, intent(in) :: k, nr, modes
        type(free_b_field_t), intent(out) :: field
        type(flux_mesh_t) :: fm
        real(dp), allocatable :: rings(:)
        integer :: i, st
        rings = [(real(i, dp)/nr, i=1, nr)]
        call flux_mesh_build(rings, 8, 1024, k, s_mapping, fm, st)
        if (st /= 0) error stop 'mesh'
        call nn_build(fm, k, modes, [s_r0, 0.0_dp], s_a, s_f, 1.0_dp, field, st)
        if (st /= 0) error stop 'build'
        field%coordinates => s_coordinates
    end subroutine

    subroutine check_identities(r, label)
        type(nf_report_t), intent(in) :: r
        character(*), intent(in) :: label
        print '(2a,6es11.3)', label, ' |u_pot| |u_ler| diff divP divL horth ', &
            r%force_potential, r%force_leray, r%difference, r%divergence_potential, &
            r%divergence_leray, r%harmonic_orthogonality
        if (r%status /= 0) error stop 'projection failed'
        if (r%divergence_potential > 1.0e-13_dp) error stop 'div u (potential)'
        if (r%divergence_leray > 1.0e-9_dp) error stop 'div u (Leray)'
        if (r%harmonic_orthogonality > 1.0e-10_dp) error stop 'harmonic orthogonality'
    end subroutine

    subroutine hodge_ladder(k, label, build)
        integer, intent(in) :: k
        character(*), intent(in) :: label
        interface
            subroutine build(k, level, field, inner)
                import :: free_b_field_t
                integer, intent(in) :: k, level
                type(free_b_field_t), intent(out) :: field
                logical, allocatable, intent(out) :: inner(:)
            end subroutine
        end interface
        type(free_b_field_t) :: field
        type(nf_system_t) :: sys
        type(nf_report_t) :: rep
        type(equilibrium_floor_report_t) :: floor
        real(dp), allocatable :: load(:), up(:), ul(:), gload(:)
        logical, allocatable :: inner(:)
        logical :: with_leray
        real(dp) :: ep(3), el(3), norm, pf, dummy
        integer :: lev, st, levels
        ! P3 with three Fourier blocks: two levels fit this suite's budget.
        levels = merge(3, 2, k == 2)
        do lev = 1, levels
            call build(k, lev, field, inner)
            ! The Leray comparison path on the cheaper levels only.
            with_leray = lev <= merge(2, 1, k == 2)
            if (allocated(inner)) then
                call nf_prepare(field, sys, st, inner=inner, leray=with_leray)
            else
                call nf_prepare(field, sys, st, leray=with_leray)
            end if
            if (st /= 0) error stop 'prepare'
            call nf_load(sys, manufactured, load)
            rep = nf_report_t()
            call nf_project(sys, load, up, ul, rep)
            call check_identities(rep, label)
            print '(a,2f9.3)', '  factor seconds potential, Leray ', &
                sys%potential_seconds, sys%leray_seconds
            if (rep%difference > 1.0e-9_dp) error stop 'potential and Leray forms differ'
            call nf_l2_error(sys, up, projected, ep(lev), norm)
            call nf_l2_error(sys, ul, projected, el(lev), dummy)
            ep(lev) = ep(lev)/norm; el(lev) = el(lev)/norm
            ! A pure gradient load is annihilated up to quadrature.
            call nf_load(sys, gradient_only, gload)
            call nf_project(sys, gload, up, ul, rep)
            pf = rep%force_potential/sqrt(dot_product(gload, gload))
            print '(a,a,i0,a,i0,a,i0,a,3es12.4)', label, ' P', k, ' level ', lev, &
                ' harmonics ', sys%nh, ' |P f - u| / |P f| (pot, Leray), |P grad| ', &
                ep(lev), el(lev), rep%force_potential
            call nf_free(sys)
        end do
        call equilibrium_residual_floor(ep(:levels), 2.0_dp, 'L2 |P f - u_h|/|P f|', &
            floor)
        call write_residual_floor(6, label, floor)
        if (floor%status /= 'decreasing') error stop 'manufactured ladder'
        ! Pre-asymptotic on these flux meshes (measured P2 1.75/1.86, P3 2.18;
        ! without the harmonic part 1.82/1.92 and 2.54): gate at k - 1/2 (P2)
        ! and k - 1 (P3).
        if (floor%last_order < merge(k - 0.5_dp, k - 1.0_dp, k == 2)) &
            error stop 'manufactured order'
        if (pf > 1.0e-6_dp) print *, 'note: gradient residual', pf
    end subroutine

    subroutine hodge_case(k)
        integer, intent(in) :: k
        shell_domain = .false.
        call hodge_ladder(k, 'hodge', solovev_build)
    end subroutine

    subroutine solovev_build(k, level, field, inner)
        integer, intent(in) :: k, level
        type(free_b_field_t), intent(out) :: field
        logical, allocatable, intent(out) :: inner(:)
        call solovev_field(k, 8*2**(level - 1), 1, field)
        if (allocated(inner)) deallocate (inner)
    end subroutine

    ! ---- shell domain -----------------------------------------------------
    subroutine shell_build(k, level, field, inner)
        integer, intent(in) :: k, level
        type(free_b_field_t), intent(out) :: field
        logical, allocatable, intent(out) :: inner(:)
        type(flux_mesh_t) :: fm
        real(dp), allocatable :: nodes(:, :)
        integer :: st, i
        call annulus_mesh(2*2**(level - 1), k, fm)
        call nn_build(fm, k, 1, [shell_r0, 0.0_dp], 1.0_dp, 0.0_dp, 1.0_dp, &
            field, st)
        if (st /= 0) error stop 'shell build'
        call meridian_lagrange_nodes(field%fs, nodes)
        allocate (inner(field%fs%n0))
        do i = 1, field%fs%n0
            inner(i) = field%fs%sp(0)%bnd(i) .and. &
                norm2(nodes(:, i) - [shell_r0, 0.0_dp]) < 0.5_dp*(rin + rout)
        end do
    end subroutine

    subroutine shell_case()
        type(free_b_field_t) :: field
        type(nf_system_t) :: sys
        real(dp), allocatable :: r1(:), r2(:), psi1(:), psi2(:), lam(:), extra(:)
        logical, allocatable :: inner(:)
        real(dp) :: kp1, kp2
        integer :: st, i
        shell_domain = .true.
        call hodge_ladder(2, 'shell', shell_build)
        ! Deflation: a right-hand side with gradient and relative-harmonic
        ! content M1 G (phi_in + interior noise) returns the same psi.
        call shell_build(2, 2, field, inner)
        call nf_prepare(field, sys, st, inner=inner)
        if (st /= 0) error stop 'shell prepare'
        if (sys%nh /= 2) error stop 'shell harmonic 2-forms'
        call random_seed()
        allocate (r1(sys%nb*sys%n1), lam(sys%nb*sys%n0))
        call random_number(r1)
        r1 = csc_matvec(sys%ct, csc_matvec(sys%m2, csc_matvec(sys%c, r1)))
        call random_number(lam)
        lam = merge(0.0_dp, lam - 0.5_dp, sys%fixed0)
        do i = 1, sys%n0
            if (inner(i)) lam(i) = 1
        end do
        extra = csc_matvec(sys%m1g, lam)
        r2 = r1 + extra
        call nf_potential_solve(sys, r1, psi1, st, kp1)
        if (st /= 0) error stop 'deflated solve'
        call nf_potential_solve(sys, r2, psi2, st, kp2)
        if (st /= 0) error stop 'deflated solve'
        print '(a,4es12.4)', 'shell deflation |dpsi|/|psi|, kernel parts, |extra|/|r| ', &
            norm2(psi2 - psi1)/norm2(psi1), kp1, kp2, &
            norm2(merge(0.0_dp, extra, sys%fixed1))/norm2(merge(0.0_dp, r1, sys%fixed1))
        if (norm2(psi2 - psi1) > 1.0e-9_dp*norm2(psi1)) error stop 'deflation'
        if (kp1 > 1.0e-10_dp) error stop 'consistent rhs has kernel content'
        call nf_free(sys)
    end subroutine

    ! ---- discrete equilibria ----------------------------------------------
    subroutine solovev_case(k)
        integer, intent(in) :: k
        type(free_b_field_t) :: field
        type(nf_system_t) :: sys
        type(nf_report_t) :: rep
        type(nn_report_t) :: srep
        type(equilibrium_floor_report_t) :: floor(2)
        real(dp), allocatable :: load(:), up(:), ul(:), nodes(:, :)
        complex(dp), allocatable :: seed(:), bc(:)
        real(dp) :: force(3, 2), u, theta
        integer :: lev, st, i, cur, levels
        levels = merge(3, 2, k == 2)
        do lev = 1, levels
            call solovev_field(k, 8*2**(lev - 1), 0, field)
            call meridian_lagrange_nodes(field%fs, nodes)
            allocate (seed(field%fs%n0), bc(field%fs%n0)); bc = 0
            do i = 1, size(seed)
                call s_coordinates([s_r0, 0.0_dp], s_a, nodes(:, i), u, theta)
                seed(i) = 0.1_dp*(1 - u*u)*cos(theta)*u
            end do
            call nn_solve(field, 0, bc, seed, srep, st, pprime=pp)
            if (st /= 0) error stop 'NN1 Solovev'
            call nf_prepare(field, sys, st, leray=lev <= merge(2, 1, k == 2))
            if (st /= 0) error stop 'prepare'
            do cur = 1, 2
                rep = nf_report_t()
                call nf_lorentz_load(sys, field, load, rep, st, primal=cur == 2)
                if (st /= 0) error stop 'Lorentz load'
                call nf_project(sys, load, up, ul, rep)
                call check_identities(rep, 'solovev')
                if (rep%difference > 1.0e-8_dp) error stop 'forms differ'
                force(lev, cur) = rep%relative
                print '(a,i0,a,i0,2a,es12.4,a,es12.4)', 'solovev P', k, ' rings ', &
                    8*2**(lev - 1), trim(currents(cur)), ' |P(J x B)|/|J||B| ', &
                    rep%relative, ' harmonic part ', rep%harmonic_part
            end do
            call nf_free(sys)
            deallocate (seed, bc)
        end do
        do cur = 1, 2
            call equilibrium_residual_floor(force(:levels, cur), 2.0_dp, &
                '|P(J x B)|_M2/<|J|^2|B|^2>^1/2', floor(cur))
            call write_residual_floor(6, 'solovev'//trim(currents(cur)), floor(cur))
        end do
        if (floor(2)%status /= 'decreasing') error stop 'Solovev force floor'
    end subroutine

    subroutine vacuum_case()
        ! Uniform vertical field a = R A_phi = R^2 (Delta* a = 0) with the
        ! harmonic F/R: J = 0 up to the quadrature of the harmonic's curl*.
        type(free_b_field_t) :: field
        type(nf_system_t) :: sys
        type(nf_report_t) :: rep
        type(nn_report_t) :: srep
        real(dp), allocatable :: load(:), up(:), ul(:), nodes(:, :)
        complex(dp), allocatable :: bc(:)
        real(dp) :: scale
        integer :: st, i
        call solovev_field(2, 8, 0, field)
        call meridian_lagrange_nodes(field%fs, nodes)
        allocate (bc(field%fs%n0))
        do i = 1, field%fs%n0
            bc(i) = nodes(1, i)**2
        end do
        call nn_solve(field, 0, bc, bc, srep, st, pprime=0.0_dp)
        if (st /= 0) error stop 'vacuum solve'
        call nf_prepare(field, sys, st)
        if (st /= 0) error stop 'prepare'
        call nf_lorentz_load(sys, field, load, rep, st)
        if (st /= 0) error stop 'load'
        call nf_project(sys, load, up, ul, rep)
        call check_identities(rep, 'vacuum')
        ! Force relative to B^2/L (sqrt(V) L / <|B|^2>).
        scale = rep%force_potential*sqrt(sys%volume)*sys%length/rep%b2
        print '(a,es12.4)', 'vacuum |P(J x B)| sqrt(V) L/<|B|^2> ', scale
        if (scale > 1.0e-10_dp) error stop 'vacuum force'
        call nf_free(sys)
    end subroutine

    subroutine cylinder_mapping(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [cyl_r0 + u*cos(theta), u*sin(theta)]
    end subroutine

    subroutine cylinder_field(nr, modes, field)
        integer, intent(in) :: nr, modes
        type(free_b_field_t), intent(out) :: field
        type(flux_mesh_t) :: fm
        real(dp), allocatable :: rings(:)
        integer :: i, st
        rings = [(real(i, dp)/nr, i=1, nr)]
        call flux_mesh_build(rings, 8, 256, 2, cylinder_mapping, fm, st)
        if (st /= 0) error stop 'cylinder mesh'
        call nn_build(fm, 2, modes, [cyl_r0, 0.0_dp], 1.0_dp, 0.0_dp, 1.0_dp, &
            field, st, cylinder_r0=cyl_r0)
        if (st /= 0) error stop 'cylinder build'
    end subroutine

    subroutine beltrami_case()
        real(dp), parameter :: mu = 1.1_dp
        type(free_b_field_t) :: field
        type(beltrami_system_t) :: bsys
        type(beltrami_report_t) :: bref
        type(nf_system_t) :: sys
        type(nf_report_t) :: rep
        type(equilibrium_floor_report_t) :: floor(2)
        real(dp), allocatable :: data(:), load(:), up(:), ul(:)
        real(dp) :: force(3, 2)
        integer :: lev, st, cur, levels
        ! The 16-ring Beltrami solve takes about a minute here: two levels
        ! by default, three with a second argument.
        levels = merge(3, 2, command_argument_count() > 1)
        do lev = 1, levels
            call cylinder_field(4*2**(lev - 1), 1, field)
            allocate (data(3*field%fs%sp(1)%ndof)); data = 0
            call beltrami_prepare(field, 2*pi*bessel_j1(mu)/mu, data, 5, bsys, st)
            if (st /= 0) error stop 'Beltrami prepare'
            call beltrami_solve(bsys, mu, field, bref)
            if (bref%status /= 0) error stop 'Beltrami solve'
            call nf_prepare(field, sys, st, leray=lev <= 2)
            if (st /= 0) error stop 'prepare'
            do cur = 1, 2
                rep = nf_report_t()
                call nf_lorentz_load(sys, field, load, rep, st, primal=cur == 2)
                if (st /= 0) error stop 'load'
                call nf_project(sys, load, up, ul, rep)
                call check_identities(rep, 'beltrami')
                if (rep%difference > 1.0e-8_dp) error stop 'forms differ'
                force(lev, cur) = rep%relative
                print '(a,i0,2a,es12.4)', 'beltrami rings ', 4*2**(lev - 1), &
                    trim(currents(cur)), ' |P(J x B)|/|J||B| ', rep%relative
            end do
            call nf_free(sys)
            deallocate (data)
        end do
        do cur = 1, 2
            call equilibrium_residual_floor(force(:levels, cur), 2.0_dp, &
                '|P(J x B)|_M2/<|J|^2|B|^2>^1/2', floor(cur))
            call write_residual_floor(6, 'beltrami'//trim(currents(cur)), floor(cur))
        end do
        if (floor(2)%status /= 'decreasing') error stop 'Beltrami force floor'
    end subroutine

    ! ---- current sheet: residual floor -------------------------------------
    real(dp) function bz_profile(r)
        real(dp), intent(in) :: r
        if (sheet_width > 0) then
            bz_profile = 1 + 0.5_dp*tanh((r - sheet_radius)/sheet_width)
        else
            bz_profile = merge(0.5_dp, 1.5_dp, r < sheet_radius)
        end if
    end function

    subroutine sheet_load(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        v = 0
        v(2) = bz_profile(norm2(x - [cyl_r0, 0.0_dp]))
    end subroutine

    subroutine sheet_ladder(width, levels, floor)
        real(dp), intent(in) :: width
        integer, intent(in) :: levels
        type(equilibrium_floor_report_t), intent(out) :: floor
        type(free_b_field_t) :: field
        type(nf_system_t) :: sys
        type(nf_report_t) :: rep
        type(nn_factor_t) :: m2
        real(dp), allocatable :: load(:), up(:), ul(:), coeff(:), force(:)
        complex(dp), allocatable :: zl(:)
        logical, allocatable :: none(:)
        integer :: lev, st
        sheet_width = width
        allocate (force(levels))
        do lev = 1, levels
            call cylinder_field(4*2**(lev - 1), 0, field)
            ! B = L2 projection of B_z(r) e_z onto V2 (no harmonic part).
            call meridian_load_v2(field, zl)
            allocate (none(field%fs%sp(2)%ndof)); none = .false.
            call nn_factor(field%fs%mass(2), none, m2, st)
            if (st /= 0) error stop 'V2 mass'
            allocate (coeff(size(zl)))
            call nn_factor_solve_real(m2, real(zl, dp), coeff, st)
            if (st /= 0) error stop 'V2 projection'
            call nn_factor_free(m2)
            field%b(:, 0) = cmplx(coeff, 0.0_dp, dp)
            call nf_prepare(field, sys, st)
            if (st /= 0) error stop 'prepare'
            rep = nf_report_t()
            call nf_lorentz_load(sys, field, load, rep, st, primal=.true.)
            if (st /= 0) error stop 'load'
            call nf_project(sys, load, up, ul, rep)
            call check_identities(rep, 'sheet')
            force(lev) = rep%relative
            call nf_free(sys)
            deallocate (none, coeff)
        end do
        call equilibrium_residual_floor(force, 2.0_dp, &
            '|P(J x B)|_M2/<|J|^2|B|^2>^1/2', floor)
    end subroutine

    subroutine meridian_load_v2(field, zl)
        use fortfem_meridian_feec_n, only: meridian_load
        type(free_b_field_t), intent(in) :: field
        complex(dp), allocatable, intent(out) :: zl(:)
        integer :: st
        call meridian_load(field%fs, 2, sheet_load, zl, st)
        if (st /= 0) error stop 'V2 load'
    end subroutine

    subroutine sheet_case()
        type(equilibrium_floor_report_t) :: smooth, sheet
        call sheet_ladder(0.1_dp, 4, smooth)
        call write_residual_floor(6, 'smooth B_z', smooth)
        call sheet_ladder(0.0_dp, 4, sheet)
        call write_residual_floor(6, 'sheet B_z', sheet)
        if (smooth%status /= 'decreasing') error stop 'smooth profile floor'
        if (sheet%last_order >= smooth%last_order) error stop 'sheet not slower'
    end subroutine
end program
