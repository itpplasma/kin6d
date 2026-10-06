! Radial Lagrange finite-element discrete-variable representation (FEM-DVR).
!
! Elements on [0,R] carry Gauss-Lobatto-Legendre (GLL) Lagrange polynomials;
! bridge functions join neighbouring elements continuously. The node r=0 is
! removed, so every basis function vanishes at the origin (u_l(0)=0 regular
! partial-wave condition); the node r=R is removed (hard wall u(R)=0).
!
! In the first element every retained basis function is a polynomial with a
! zero at r=0. Hence chi_i chi_j/r and chi_i chi_j/r^2 are polynomials of
! degree <= 2p-1 and the p+1 point GLL rule integrates the Coulomb and
! centrifugal matrix elements there exactly: the singularity is handled by the
! basis, not softened. The kinetic matrix int chi_i' chi_j' dr is exact on
! every element. The remaining approximation is the usual diagonal (lumped)
! GLL overlap/potential quadrature on smooth elements, which converges
! spectrally with the order and is measured by the oracles in the tests.
!
! The GLL rule and the normalized associated Legendre recurrence used by the
! encounter are generic numerics pending an upstream FortNum home; they are
! kept private and minimal here (see docs/coulomb-encounter.md, open issues).
module kin6d_radial_dvr
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_coulomb_encounter_leaves, only: ce_radial_potential, ce_kinetic_scale
    use fortnum_quadrature, only: gauss_lobatto_legendre
    use fortnum_polynomial, only: barycentric_weights, lagrange_differentiation_matrix
    implicit none
    private
    public :: radial_mesh_t, build_radial_mesh, mesh_spacing
    public :: radial_band_hamiltonian, radial_band_matvec, radial_eigen_window
    public :: radial_locate, radial_dimension

    type :: radial_mesh_t
        integer :: order = 0          ! GLL points per element (polynomial degree + 1)
        integer :: nelem = 0
        integer :: n = 0              ! retained DVR functions
        real(dp) :: rmax = 0
        real(dp), allocatable :: edges(:)        ! (0:nelem)
        real(dp), allocatable :: r(:), w(:)      ! retained nodes and summed weights
        real(dp), allocatable :: tband(:, :)     ! lower band of -d^2/dr^2 (kd+1,n)
        real(dp), allocatable :: xi(:), wi(:), bary(:)  ! reference GLL rule
    end type radial_mesh_t

    interface
        subroutine dsbevd(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, iwork, &
                          liwork, info)
            import :: dp
            character, intent(in) :: jobz, uplo
            integer, intent(in) :: n, kd, ldab, ldz, lwork, liwork
            real(dp), intent(inout) :: ab(ldab, *)
            real(dp), intent(out) :: w(*), z(ldz, *), work(*)
            integer, intent(out) :: iwork(*), info
        end subroutine dsbevd
        subroutine dsbevx(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, &
                          abstol, m, w, z, ldz, work, iwork, ifail, info)
            import :: dp
            character, intent(in) :: jobz, range, uplo
            integer, intent(in) :: n, kd, ldab, ldq, il, iu, ldz
            real(dp), intent(inout) :: ab(ldab, *)
            real(dp), intent(out) :: q(ldq, *), w(*), z(ldz, *), work(*)
            real(dp), intent(in) :: vl, vu, abstol
            integer, intent(out) :: m, iwork(*), ifail(*), info
        end subroutine dsbevx
    end interface
contains
    pure integer function radial_dimension(mesh)
        type(radial_mesh_t), intent(in) :: mesh
        radial_dimension = mesh%n
    end function radial_dimension

    ! Element size resolving wavenumber kmax with a declared number of GLL
    ! nodes per local wavelength.
    pure real(dp) function mesh_spacing(kmax, order, nodes_per_wavelength)
        real(dp), intent(in) :: kmax, nodes_per_wavelength
        integer, intent(in) :: order
        mesh_spacing = real(order - 1, dp)*2*acos(-1.0_dp)/(nodes_per_wavelength*kmax)
    end function mesh_spacing

    ! Graded mesh: first element h0, geometric growth up to hmax, ending at rmax.
    subroutine build_radial_mesh(mesh, rmax, order, hmax, h0, growth, ok)
        type(radial_mesh_t), intent(out) :: mesh
        real(dp), intent(in) :: rmax, hmax, h0, growth
        integer, intent(in) :: order
        logical, intent(out) :: ok
        real(dp), allocatable :: edges(:)
        real(dp) :: h, rr, dref(order, order)
        integer :: ne, e, a, b, q, deg, kd, i, j
        ok = .false.
        if (order < 3 .or. order > 40) return
        if (.not. (rmax > 0 .and. hmax > 0 .and. h0 > 0 .and. growth >= 1)) return
        if (.not. ieee_is_finite(rmax/min(h0, hmax))) return
        if (rmax/min(h0, hmax) > 1.0e6_dp) return
        allocate (edges(0:ceiling(rmax/min(h0, hmax)) + 2))
        edges(0) = 0; rr = 0; h = min(h0, hmax); ne = 0
        do while (rr < rmax*(1 - 1.0e-12_dp))
            ne = ne + 1
            rr = min(rmax, rr + h)
            edges(ne) = rr
            h = min(hmax, h*growth)
        end do
        ! Merge a short last element into its predecessor.
        if (ne >= 2) then
            if (edges(ne) - edges(ne - 1) < 0.3_dp*(edges(ne - 1) - edges(ne - 2))) then
                edges(ne - 1) = edges(ne); ne = ne - 1
            end if
        end if
        deg = order - 1; kd = deg
        mesh%order = order; mesh%nelem = ne; mesh%rmax = rmax
        mesh%n = ne*deg - 1
        allocate (mesh%edges(0:ne), mesh%xi(order), mesh%wi(order), mesh%bary(order))
        mesh%edges = edges(0:ne)
        call gauss_lobatto_legendre(order, mesh%xi, mesh%wi)
        call barycentric_weights(order, mesh%xi, mesh%bary)
        call lagrange_differentiation_matrix(order, mesh%xi, dref)  ! dref(q,a) = L_a'(xi_q)
        allocate (mesh%r(mesh%n), mesh%w(mesh%n), mesh%tband(kd + 1, mesh%n))
        mesh%w = 0; mesh%tband = 0
        do e = 1, ne
            h = mesh%edges(e) - mesh%edges(e - 1)
            do a = 1, order
                i = (e - 1)*deg + a - 1
                if (i < 1 .or. i > mesh%n) cycle
                mesh%r(i) = mesh%edges(e - 1) + h*(mesh%xi(a) + 1)/2
                mesh%w(i) = mesh%w(i) + h*mesh%wi(a)/2
            end do
        end do
        do e = 1, ne
            h = mesh%edges(e) - mesh%edges(e - 1)
            do b = 1, order
                j = (e - 1)*deg + b - 1
                if (j < 1 .or. j > mesh%n) cycle
                do a = b, order
                    i = (e - 1)*deg + a - 1
                    if (i < 1 .or. i > mesh%n) cycle
                    mesh%tband(1 + i - j, j) = mesh%tband(1 + i - j, j) + &
                        2/h*sum(mesh%wi*dref(:, a)*dref(:, b))
                end do
            end do
        end do
        do j = 1, mesh%n
            do i = j, min(mesh%n, j + kd)
                mesh%tband(1 + i - j, j) = mesh%tband(1 + i - j, j)/ &
                    sqrt(mesh%w(i)*mesh%w(j))
            end do
        end do
        ok = all(mesh%w > 0) .and. all(ieee_is_finite(mesh%tband))
    end subroutine build_radial_mesh

    ! Lower band of H_l = (hbar^2/2mu) T + diag(V_l(r_i)) in the normalized DVR.
    subroutine radial_band_hamiltonian(mesh, l, g, mu, hbar, hb)
        type(radial_mesh_t), intent(in) :: mesh
        integer, intent(in) :: l
        real(dp), intent(in) :: g, mu, hbar
        real(dp), intent(out) :: hb(:, :)
        real(dp) :: c, v
        integer :: i
        call ce_kinetic_scale(mu, hbar, c)
        hb = c*mesh%tband
        do i = 1, mesh%n
            call ce_radial_potential(mesh%r(i), real(l, dp), g, mu, hbar, v)
            hb(1, i) = hb(1, i) + v
        end do
    end subroutine radial_band_hamiltonian

    ! y = A x for a symmetric matrix in lower band storage, several columns.
    pure subroutine radial_band_matvec(hb, x, y)
        real(dp), intent(in) :: hb(:, :), x(:, :)
        real(dp), intent(out) :: y(:, :)
        integer :: i, j, n, kd
        n = size(hb, 2); kd = size(hb, 1) - 1
        y = 0
        do j = 1, n
            y(j, :) = y(j, :) + hb(1, j)*x(j, :)
            do i = j + 1, min(n, j + kd)
                y(i, :) = y(i, :) + hb(1 + i - j, j)*x(j, :)
                y(j, :) = y(j, :) + hb(1 + i - j, j)*x(i, :)
            end do
        end do
    end subroutine radial_band_matvec

    ! Eigenpairs of H_l with energy in (emin, emax] (all if the window covers the
    ! Gershgorin interval), with column residual norms ||H v_j - e_j v_j||_2 and the
    ! orthogonality defect max|V^T V - I| of the returned columns.
    subroutine radial_eigen_window(mesh, l, g, mu, hbar, emax, e, v, resid, &
                                   orthogonality, ok, emin, estimator_seconds)
        type(radial_mesh_t), intent(in) :: mesh
        integer, intent(in) :: l
        real(dp), intent(in) :: g, mu, hbar, emax
        real(dp), intent(in), optional :: emin
        real(dp), allocatable, intent(out) :: e(:), v(:, :), resid(:)
        real(dp), intent(out) :: orthogonality
        logical, intent(out) :: ok
        real(dp), intent(out), optional :: estimator_seconds
        real(dp), allocatable :: hb(:, :), ab(:, :), q(:, :), z(:, :), work(:), w(:)
        real(dp), allocatable :: hv(:, :)
        integer, allocatable :: iwork(:), ifail(:)
        real(dp) :: vl, vu, glo, ghi, off, estimator_start, estimator_end
        integer :: n, kd, m, info, j, i, i0
        ok = .false.; orthogonality = huge(1.0_dp)
        if (present(estimator_seconds)) estimator_seconds = 0
        n = mesh%n; kd = mesh%order - 1
        allocate (hb(kd + 1, n), ab(kd + 1, n))
        call radial_band_hamiltonian(mesh, l, g, mu, hbar, hb)
        if (.not. all(ieee_is_finite(hb))) return
        glo = huge(1.0_dp); ghi = -huge(1.0_dp)
        do i = 1, n
            off = 0
            do j = max(1, i - kd), min(n, i + kd)
                if (j < i) off = off + abs(hb(1 + i - j, j))
                if (j > i) off = off + abs(hb(1 + j - i, i))
            end do
            glo = min(glo, hb(1, i) - off); ghi = max(ghi, hb(1, i) + off)
        end do
        vl = glo - 1 - abs(glo)*1.0e-6_dp
        if (present(emin)) vl = max(vl, emin)
        vu = min(emax, ghi + 1 + abs(ghi)*1.0e-6_dp)
        if (.not. (vu > vl)) return
        ! Count the window from eigenvalues only. Small windows use bisection plus
        ! inverse iteration (high accuracy for long-time phases); large windows
        ! use divide-and-conquer on the full band (measured faster there).
        ab = hb
        allocate (w(n), work(max(7*n, 1 + 5*n + 2*n*n)), iwork(5*n + 3), ifail(n))
        allocate (z(1, 1))
        call dsbevd('N', 'L', n, kd, ab, kd + 1, w, z, 1, work, size(work), iwork, &
                    size(iwork), info)
        if (info /= 0) return
        m = count(w > vl .and. w <= vu)
        deallocate (z)
        ab = hb
        if (5*m <= n) then
            allocate (q(n, n), z(n, max(1, m)))
            call dsbevx('V', 'V', 'L', n, kd, ab, kd + 1, q, n, vl, vu, 0, 0, &
                        2*tiny(1.0_dp), m, w, z, n, work, iwork, ifail, info)
            deallocate (q)
            if (info /= 0) return
            i0 = 1
        else
            allocate (z(n, n))
            call dsbevd('V', 'L', n, kd, ab, kd + 1, w, z, n, work, size(work), iwork, &
                        size(iwork), info)
            if (info /= 0) return
            i0 = n + 1; m = 0
            do j = 1, n
                if (w(j) > vl .and. w(j) <= vu) then
                    i0 = min(i0, j); m = m + 1
                end if
            end do
        end if
        deallocate (work)
        allocate (e(m), v(n, m), resid(m), hv(n, m))
        e = w(i0:i0 + m - 1); v = z(:, i0:i0 + m - 1)
        deallocate (z)
        call cpu_time(estimator_start)
        call radial_band_matvec(hb, v, hv)
        do j = 1, m
            resid(j) = norm2(hv(:, j) - e(j)*v(:, j))
        end do
        orthogonality = 0
        do j = 1, m
            do i = 1, j
                off = dot_product(v(:, i), v(:, j))
                if (i == j) off = off - 1
                orthogonality = max(orthogonality, abs(off))
            end do
        end do
        call cpu_time(estimator_end)
        if (present(estimator_seconds)) &
            estimator_seconds = estimator_end-estimator_start
        ok = all(ieee_is_finite(e)) .and. all(ieee_is_finite(resid))
    end subroutine radial_eigen_window

    ! Values of the retained basis functions at radius rq: u(rq) = sum val*a(idx),
    ! with a the normalized DVR coefficients. Exact interpolation of the represented
    ! polynomial; idx entries of removed boundary nodes carry val=0.
    pure subroutine radial_locate(mesh, rq, idx, val)
        type(radial_mesh_t), intent(in) :: mesh
        real(dp), intent(in) :: rq
        integer, intent(out) :: idx(:)
        real(dp), intent(out) :: val(:)
        real(dp) :: x, h, den, t
        integer :: lo, hi, mid, e, a, gnode, deg
        logical :: at_node
        deg = mesh%order - 1
        idx = 1; val = 0
        if (.not. (rq >= 0 .and. rq <= mesh%rmax)) return
        lo = 1; hi = mesh%nelem
        do while (lo < hi)
            mid = (lo + hi)/2
            if (rq > mesh%edges(mid)) then
                lo = mid + 1
            else
                hi = mid
            end if
        end do
        e = lo
        h = mesh%edges(e) - mesh%edges(e - 1)
        x = 2*(rq - mesh%edges(e - 1))/h - 1
        den = 0; at_node = .false.
        do a = 1, mesh%order
            if (.not. (abs(x - mesh%xi(a)) > 0)) then
                val = 0; val(a) = 1; at_node = .true.; exit
            end if
            t = mesh%bary(a)/(x - mesh%xi(a))
            val(a) = t; den = den + t
        end do
        if (.not. at_node) val = val/den
        do a = 1, mesh%order
            gnode = (e - 1)*deg + a - 1
            if (gnode < 1 .or. gnode > mesh%n) then
                val(a) = 0
            else
                idx(a) = gnode
                val(a) = val(a)/sqrt(mesh%w(gnode))
            end if
        end do
    end subroutine radial_locate
end module kin6d_radial_dvr
