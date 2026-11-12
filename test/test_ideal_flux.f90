program test_ideal_flux
    ! P1e.3 checks of the flux-coordinate leaves and assembly. The leaves are
    ! compared with an independent dual-basis evaluation of B and with central
    ! differences on a non-axisymmetric analytic chart; the assembly on a
    ! straight-cylinder chart must reproduce the cylinder results of P1e.2.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortsparse, only: csc_z_t, fortsparse_status_t, FORTSPARSE_OK
    use milca_ideal_flux, only: flux_space_t, flux_space_init, flux_assemble, &
        flux_point_t, flux_equilibrium
    use milca_ideal_flux_leaves, only: flux_jets_array
    implicit none
    interface
        subroutine zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
                rwork, info)
            import :: dp
            integer, intent(in) :: itype, n, lda, ldb, lwork
            character, intent(in) :: jobz, uplo
            complex(dp), intent(inout) :: a(lda, *), b(ldb, *)
            real(dp), intent(out) :: w(*), rwork(*)
            complex(dp), intent(inout) :: work(*)
            integer, intent(out) :: info
        end subroutine zhegv
    end interface
    real(dp), parameter :: pi = acos(-1.0_dp)
    ! Analytic chart parameters.
    real(dp), parameter :: r0 = 3.0_dp, aa = 1.0_dp, dd = 0.15_dp, kap = 0.2_dp, &
        phi_edge = 1.3_dp, iota0 = 0.4_dp, iota2 = 0.3_dp
    ! Cylinder chart (S2.1-3 of docs/s0-spectrum.md): a = 2, R0 = 2.5, n = -1.
    real(dp), parameter :: kink_reference = -0.1949752661_dp
    real(dp) :: cyl_a, cyl_r0
    integer :: cyl_case
    character(32) :: mode
    mode = 'all'
    if (command_argument_count() > 0) call get_command_argument(1, mode)
    if (mode == 'all' .or. mode == 'leaves') call leaves()
    if (mode == 'all' .or. mode == 'kink') call cylinder_kink()
    if (mode == 'all' .or. mode == 'uniform') call uniform_harmonics()
    print '(a)', 'PASS ideal flux ('//trim(mode)//')'
contains
    subroutine leaves()
        real(dp) :: u(3), xu(3, 3), xuu(3, 3, 3), lam(7), jac, b(3), g(3, 3), &
            q(3), t(12, 12), bref(3), dual(3, 3), bp(3), bm(3), du(3, 3), gfd(3, 3), &
            h, psi1, psi2, iota, iota1, xp(3), xm(3), afd(3, 3), jets(12), d(3)
        integer :: k, s
        h = 1e-5_dp
        do s = 1, 4
            u = [0.3_dp + 0.15_dp*s, 0.7_dp*s, 0.4_dp*s - 0.2_dp]
            call chart(u, xu, xuu, psi1, psi2, iota, iota1, lam)
            call flux_jets_array(xu, xuu, psi1, psi2, iota, iota1, 0.0_dp, lam, &
                jac, b, g, q, t)
            ! Independent B from the dual basis of the numerically inverted chart.
            dual = transpose(inverse(xu))
            bref = psi1*cross(dual(:, 1), dual(:, 2) + lam(3)*dual(:, 1) &
                + lam(1)*dual(:, 2) + lam(2)*dual(:, 3)) &
                + iota*psi1*cross(dual(:, 3), dual(:, 1))
            if (maxval(abs(b - bref)) > 1e-13_dp*maxval(abs(bref))) &
                error stop 'B differs from the dual-basis formula'
            if (abs(g(1,1) + g(2,2) + g(3,3)) > 1e-12_dp*maxval(abs(g))) &
                error stop 'div B /= 0'
            ! d_u B by central differences, then g = sum_u d_u B grad u.
            do k = 1, 3
                bp = field(u + h*unit(k)); bm = field(u - h*unit(k))
                du(:, k) = (bp - bm)/(2*h)
            end do
            gfd = matmul(du, transpose(dual))
            if (maxval(abs(g - gfd)) > 1e-7_dp*maxval(abs(g))) &
                error stop 'grad B differs from central differences'
            ! xi from densities D = d0 + rho d1 (componentwise) and its gradient.
            jets = 0
            d = density(u)
            jets(1:3) = d
            do k = 1, 3
                jets(3*k + 1:3*k + 3) = (density(u + h*unit(k)) &
                    - density(u - h*unit(k)))/(2*h)
            end do
            do k = 1, 3
                xp = displacement(u + h*unit(k)); xm = displacement(u - h*unit(k))
                du(:, k) = (xp - xm)/(2*h)
            end do
            afd = matmul(du, transpose(dual))
            if (maxval(abs(matmul(t(1:3, :), jets) - displacement(u))) > 1e-13_dp) &
                error stop 'xi from densities'
            if (maxval(abs(reshape(matmul(t(4:12, :), jets), [3, 3]) - afd)) &
                > 1e-7_dp*maxval(abs(afd))) error stop 'grad xi from densities'
            if (abs(jac*(t(4, 4) + t(8, 4) + t(12, 4)) - 1) > 1e-12_dp) &
                error stop 'sqrt(g) div xi /= d_rho D^rho'
            if (abs(jac*(t(4, 8) + t(8, 8) + t(12, 8)) - 1) > 1e-12_dp) &
                error stop 'sqrt(g) div xi /= d_theta D^theta'
            if (abs(jac*(t(4, 12) + t(8, 12) + t(12, 12)) - 1) > 1e-12_dp) &
                error stop 'sqrt(g) div xi /= d_zeta D^zeta'
        end do
        print '(a)', 'leaves: B, grad B, xi and grad xi agree; div B = 0'
    end subroutine leaves

    function field(u) result(bv)
        real(dp), intent(in) :: u(3)
        real(dp) :: bv(3), xu(3, 3), xuu(3, 3, 3), lam(7), jac, g(3, 3), q(3), &
            t(12, 12), psi1, psi2, iota, iota1
        call chart(u, xu, xuu, psi1, psi2, iota, iota1, lam)
        call flux_jets_array(xu, xuu, psi1, psi2, iota, iota1, 0.0_dp, lam, jac, &
            bv, g, q, t)
    end function field

    function density(u) result(d)
        real(dp), intent(in) :: u(3)
        real(dp) :: d(3)
        d = [0.2_dp + u(1)*cos(u(2)), -0.4_dp*u(1)**2 + sin(u(3)), &
            0.3_dp*u(1)*sin(u(2) - u(3))]
    end function density

    function displacement(u) result(xi)
        real(dp), intent(in) :: u(3)
        real(dp) :: xi(3), xu(3, 3), xuu(3, 3, 3), lam(7), psi1, psi2, iota, &
            iota1, jac
        call chart(u, xu, xuu, psi1, psi2, iota, iota1, lam)
        jac = dot_product(xu(:, 1), cross(xu(:, 2), xu(:, 3)))
        xi = matmul(xu, density(u))/jac
    end function displacement

    subroutine chart(u, xu, xuu, psi1, psi2, iota, iota1, lam)
        ! R = r0 + a rho cos t + d rho cos(t - 2z), Z = a rho sin t
        ! - d rho sin(t - 2z), phi = z; lambda = kap rho^2 sin(t - z);
        ! psi = phi_edge rho^2 / 2, iota = iota0 + iota2 rho^2.
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: xu(3, 3), xuu(3, 3, 3), psi1, psi2, iota, &
            iota1, lam(7)
        real(dp) :: r, t, z, rr(3), rz(3, 3), zz(3), zzz(3, 3), cr(3), c, s, c2, s2
        integer :: i, j
        r = u(1); t = u(2); z = u(3)
        c = cos(t); s = sin(t); c2 = cos(t - 2*z); s2 = sin(t - 2*z)
        ! R and Z with first and second derivatives in (rho, theta, zeta).
        rr = [aa*c + dd*c2, -aa*r*s - dd*r*s2, 2*dd*r*s2]
        rz = 0
        rz(1, 2) = -aa*s - dd*s2; rz(1, 3) = 2*dd*s2
        rz(2, 2) = -aa*r*c - dd*r*c2; rz(2, 3) = 2*dd*r*c2; rz(3, 3) = -4*dd*r*c2
        zz = [aa*s - dd*s2, aa*r*c - dd*r*c2, 2*dd*r*c2]
        zzz = 0
        zzz(1, 2) = aa*c - dd*c2; zzz(1, 3) = 2*dd*c2
        zzz(2, 2) = -aa*r*s + dd*r*s2; zzz(2, 3) = -2*dd*r*s2; zzz(3, 3) = 4*dd*r*s2
        do j = 1, 3
            do i = 1, j - 1
                rz(j, i) = rz(i, j); zzz(j, i) = zzz(i, j)
            end do
        end do
        cr = [cos(z), sin(z), 0.0_dp]
        block
            real(dp) :: big, e(3), ep(3)
            big = r0 + aa*r*c + dd*r*c2
            e = cr; ep = [-sin(z), cos(z), 0.0_dp]
            do i = 1, 3
                xu(:, i) = rr(i)*e + [0.0_dp, 0.0_dp, zz(i)]
            end do
            xu(:, 3) = xu(:, 3) + big*ep
            do j = 1, 3
                do i = 1, 3
                    xuu(:, i, j) = rz(i, j)*e + [0.0_dp, 0.0_dp, zzz(i, j)]
                    if (i == 3) xuu(:, i, j) = xuu(:, i, j) + rr(j)*ep
                    if (j == 3) xuu(:, i, j) = xuu(:, i, j) + rr(i)*ep
                    if (i == 3 .and. j == 3) xuu(:, i, j) = xuu(:, i, j) - big*e
                end do
            end do
        end block
        psi1 = phi_edge*r; psi2 = phi_edge
        iota = iota0 + iota2*r*r; iota1 = 2*iota2*r
        ! lambda_theta, lambda_zeta, d(rho,theta,zeta) lambda_theta,
        ! d(rho,zeta) lambda_zeta.
        lam = [kap*r*r*cos(t - z), -kap*r*r*cos(t - z), &
            2*kap*r*cos(t - z), -kap*r*r*sin(t - z), kap*r*r*sin(t - z), &
            -2*kap*r*cos(t - z), -kap*r*r*sin(t - z)]
    end subroutine chart

    subroutine cylinder_kink()
        real(dp), allocatable :: lam(:)
        cyl_a = 2; cyl_r0 = 2.5_dp; cyl_case = 1
        call spectrum([1], [-1], 32, 3, 8, 0.0_dp, lam)
        print '(a,es20.11,es10.2)', 'kink, m = 1     ', lam(1), &
            abs(lam(1)/kink_reference - 1)
        if (abs(lam(1)/kink_reference - 1) > 1e-6_dp) error stop 'kink, m = 1'
        call spectrum([0, 1, 2], [-1, -1, -1], 32, 3, 16, 0.0_dp, lam)
        print '(a,es20.11,es10.2)', 'kink, m = 0..2  ', lam(1), &
            abs(lam(1)/kink_reference - 1)
        if (abs(lam(1)/kink_reference - 1) > 1e-6_dp) error stop 'kink, m = 0..2'
    end subroutine cylinder_kink

    subroutine uniform_harmonics()
        ! The straight cylinder decouples harmonics: the spectrum with m = 1, 2, 3
        ! is the union of the single-harmonic spectra.
        real(dp), allocatable :: lam(:), single(:), union(:)
        integer :: m, k
        cyl_a = 1; cyl_r0 = 1; cyl_case = 2
        allocate(union(0))
        do m = 1, 3
            call spectrum([m], [1], 12, 3, 8, 5.0_dp/3, single)
            union = [union, single]
        end do
        call spectrum([1, 2, 3], [1, 1, 1], 12, 3, 16, 5.0_dp/3, lam)
        call sort(union)
        if (size(lam) /= size(union)) error stop 'dof count'
        do k = 1, size(lam)
            if (abs(lam(k) - union(k)) > 1e-9_dp*max(1.0_dp, abs(union(k)))) &
                error stop 'harmonics couple in a straight cylinder'
        end do
        print '(a,es10.2)', 'uniform: union of single-harmonic spectra, max dev', &
            maxval(abs(lam - union)/max(1.0_dp, abs(union)))
    end subroutine uniform_harmonics

    subroutine cylinder(rho, theta, zeta, point)
        real(dp), intent(in) :: rho, theta, zeta
        type(flux_point_t), intent(out) :: point
        real(dp) :: r, bz, bzp, p, pp, c, s
        r = cyl_a*rho; c = cos(theta); s = sin(theta)
        if (cyl_case == 1) then
            bz = sqrt(0.25_dp + r*r); bzp = r/bz
            p = 1.5_dp*(cyl_a**2 - r*r); pp = -3*r
        else
            bz = 1; bzp = 0; p = 0.3_dp; pp = 0
        end if
        point%xu(:, 1) = [cyl_a*c, cyl_a*s, 0.0_dp]
        point%xu(:, 2) = [-r*s, r*c, 0.0_dp]
        point%xu(:, 3) = [0.0_dp, 0.0_dp, cyl_r0]
        point%xuu = 0
        point%xuu(:, 1, 2) = [-cyl_a*s, cyl_a*c, 0.0_dp]
        point%xuu(:, 2, 1) = point%xuu(:, 1, 2)
        point%xuu(:, 2, 2) = [-r*c, -r*s, 0.0_dp]
        ! psi' = a r Bz, iota = R0 Btheta/(r Bz) in rho; zeta = z/R0.
        point%psi1 = cyl_a*r*bz
        point%psi2 = cyl_a**2*(bz + r*bzp)
        if (cyl_case == 1) then
            point%iota = cyl_r0/bz
            point%iota1 = -cyl_a*cyl_r0*bzp/bz**2
        else
            point%iota = 0; point%iota1 = 0
        end if
        point%lam = 0
        point%p = p; point%pp = cyl_a*pp; point%density = 1
    end subroutine cylinder

    subroutine spectrum(m, n, ne, deg, ntheta, gam, lam)
        integer, intent(in) :: m(:), n(:), ne, deg, ntheta
        real(dp), intent(in) :: gam
        real(dp), allocatable, intent(out) :: lam(:)
        type(flux_space_t) :: space
        type(csc_z_t) :: stiffness, mass
        type(fortsparse_status_t) :: status
        complex(dp), allocatable :: ka(:, :), ma(:, :), work(:)
        real(dp), allocatable :: rwork(:), edges(:)
        integer :: ok, j, k, info
        edges = [(real(j, dp)/ne, j = 0, ne)]
        call flux_space_init(space, edges, deg, m, n, 1, ntheta, 1, ok)
        if (ok /= 0) error stop 'space'
        call flux_assemble(space, cylinder, gam, 1.0_dp, stiffness, mass, status)
        if (status%code /= FORTSPARSE_OK) error stop 'assembly'
        allocate(ka(space%ndof, space%ndof), ma(space%ndof, space%ndof))
        ka = 0; ma = 0
        do j = 1, space%ndof
            do k = stiffness%col_ptr(j), stiffness%col_ptr(j + 1) - 1
                ka(stiffness%row_idx(k), j) = stiffness%val(k)
                ma(mass%row_idx(k), j) = mass%val(k)
            end do
        end do
        allocate(lam(space%ndof), work(4*space%ndof), rwork(3*space%ndof))
        call zhegv(1, 'N', 'U', space%ndof, ka, space%ndof, ma, space%ndof, lam, &
            work, size(work), rwork, info)
        if (info /= 0) error stop 'dense generalized eigensolve'
    end subroutine spectrum

    pure function unit(k) result(e)
        integer, intent(in) :: k
        real(dp) :: e(3)
        e = 0; e(k) = 1
    end function unit

    pure function cross(p, q) result(w)
        real(dp), intent(in) :: p(3), q(3)
        real(dp) :: w(3)
        w = [p(2)*q(3) - p(3)*q(2), p(3)*q(1) - p(1)*q(3), p(1)*q(2) - p(2)*q(1)]
    end function cross

    function inverse(a) result(b)
        ! Gauss-Jordan with partial pivoting, independent of the cross-product
        ! dual basis used by the leaves.
        real(dp), intent(in) :: a(3, 3)
        real(dp) :: b(3, 3), w(3, 6), row(6)
        integer :: i, k, piv
        w(:, 1:3) = a; w(:, 4:6) = 0
        do i = 1, 3
            w(i, 3 + i) = 1
        end do
        do k = 1, 3
            piv = k - 1 + maxloc(abs(w(k:3, k)), 1)
            row = w(k, :); w(k, :) = w(piv, :); w(piv, :) = row
            w(k, :) = w(k, :)/w(k, k)
            do i = 1, 3
                if (i /= k) w(i, :) = w(i, :) - w(i, k)*w(k, :)
            end do
        end do
        b = w(:, 4:6)
    end function inverse

    subroutine sort(x)
        real(dp), intent(inout) :: x(:)
        real(dp) :: v
        integer :: i, j
        do i = 2, size(x)
            v = x(i); j = i - 1
            do while (j >= 1)
                if (x(j) <= v) exit
                x(j + 1) = x(j); j = j - 1
            end do
            x(j + 1) = v
        end do
    end subroutine sort
end program test_ideal_flux
