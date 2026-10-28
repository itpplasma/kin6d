module kin6d_ideal_fb_circuits
    !! Circuit currents as parameters of the nested-map free-boundary residual
    !! (docs/free-boundary-nested-map-3d.md). The vacuum source is linear in
    !! the circuit currents p = (I_1..I_K):
    !!     B_vac = B_0 + sum_k I_k b_k,
    !! b_k the field of circuit k at unit current (coils with winding counts),
    !! B_0 a fixed part (the unit-circulation filament times c). The residual
    !! R(x; p) of kin6d_ideal_free_boundary sees B_vac only through the Gamma_c
    !! traces (drive a, load h) and the toroidal period F_vac, all linear in
    !! B_vac, so they are precomputed once per circuit:
    !!     a(p) = a_0 + sum_k I_k a_k,  h(p) = h_0 + sum_k I_k h_k,
    !!     F_vac(p) = F_0 + sum_k I_k F_k.
    !! dR/dp is exact: the trace part is linear, the F_vac part is the
    !! generated shell Hessian leaf (ideal_fb_traces_jvp/vjp).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_coils_3d, only: coil_set_t, segment_b
    use kin6d_ideal_free_boundary, only: ideal_fb_t, ideal_fb_field_basis_t, &
        ideal_fb_sampler, ideal_fb_drive_traces, ideal_fb_set_traces, &
        ideal_fb_traces_jvp, ideal_fb_traces_vjp, ideal_fb_set_drive, &
        ideal_fb_drive_traces_3d, ideal_fb_set_drive_3d, &
        IDEAL_FB_OK, IDEAL_FB_BAD_INPUT
    implicit none
    private
    public :: coil_circuit_basis_t, coil_circuit_basis_init
    public :: ideal_fb_circuits_t, ideal_fb_circuits_build, ideal_fb_set_currents
    public :: ideal_fb_currents_jvp, ideal_fb_currents_vjp
    real(dp), parameter :: pi = acos(-1.0_dp)

    !> Unit-current circuit fields of a coils_3d filament set. Coil c belongs
    !> to circuit circuit(c) with turns(c) windings (0 = not connected), so
    !> b_k = sum_{circuit(c) = k} turns(c) B_c[1 A]. field_scale converts
    !> the SI Biot-Savart field to the consumer's mu0 (mu0/coil_mu0).
    !> Toroidal harmonics come from an nphi-point rule over one full turn
    !> (aliasing error exp(-nphi d/R), d the coil distance).
    type, extends(ideal_fb_field_basis_t) :: coil_circuit_basis_t
        type(coil_set_t) :: coils
        integer :: ncircuit = 0, nphi = 0
        integer, allocatable :: circuit(:)
        real(dp), allocatable :: turns(:)
        real(dp) :: field_scale = 1
    contains
        procedure :: nfield => coil_nfield
        procedure :: sample => coil_sample
        procedure :: sample_points => coil_sample_points
    end type

    !> Precomputed per-circuit traces on Gamma_c of one ideal_fb_t:
    !> drive(:, :, k), load(:, :, k), f_vac(k) for unit current in circuit k,
    !> and the current-independent part (index 0 arrays).
    type :: ideal_fb_circuits_t
        integer :: ncircuit = 0
        real(dp), allocatable :: drive(:, :, :), load(:, :, :), f_vac(:)
        real(dp), allocatable :: drive0(:, :), load0(:, :)
        real(dp) :: f_vac0 = 0
        ! Diagnostics: trace closure / curl-free mismatch and spread of
        ! R B_phi^(0) over the trace, per circuit.
        real(dp), allocatable :: closure(:), f_spread(:)
    end type

contains

    subroutine coil_circuit_basis_init(basis, coils, circuit, turns, ncircuit, nphi, &
                                       status, field_scale)
        type(coil_circuit_basis_t), intent(out) :: basis
        type(coil_set_t), intent(in) :: coils
        integer, intent(in) :: circuit(:), ncircuit, nphi
        real(dp), intent(in) :: turns(:)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: field_scale
        status = IDEAL_FB_BAD_INPUT
        if (coils%ncoil < 1 .or. ncircuit < 1 .or. nphi < 2) return
        if (size(circuit) /= coils%ncoil .or. size(turns) /= coils%ncoil) return
        if (any(circuit < 0) .or. any(circuit > ncircuit)) return
        if (.not. all(ieee_is_finite(turns))) return
        basis%coils = coils
        basis%circuit = circuit; basis%turns = turns
        basis%ncircuit = ncircuit; basis%nphi = nphi
        if (present(field_scale)) basis%field_scale = field_scale
        status = IDEAL_FB_OK
    end subroutine

    integer function coil_nfield(self)
        class(coil_circuit_basis_t), intent(in) :: self
        coil_nfield = self%ncircuit
    end function

    subroutine coil_sample(self, rz, nfp, ntor, bcos, bsin, status)
        !! Real cos/sin harmonics n = k nfp of every circuit field at each
        !! meridional point. Points are independent (OpenMP over points), so
        !! the result does not depend on the thread count.
        class(coil_circuit_basis_t), intent(in) :: self
        real(dp), intent(in) :: rz(:, :)
        integer, intent(in) :: nfp, ntor
        real(dp), intent(out) :: bcos(:, 0:, :, :), bsin(:, 0:, :, :)
        integer, intent(out) :: status
        real(dp), allocatable :: cp(:), sp(:), ck(:, :), sk(:, :)
        real(dp) :: p(3), bc(3), bcyl(3), bk(3, self%ncircuit)
        integer :: q, j, c, s, k, nk, np
        status = IDEAL_FB_BAD_INPUT
        np = self%nphi
        if (size(rz, 1) /= 2 .or. nfp < 1 .or. ntor < 0 .or. np < 2) return
        if (size(bcos, 1) /= 3 .or. ubound(bcos, 2) /= ntor) return
        if (size(bcos, 3) /= self%ncircuit .or. size(bcos, 4) /= size(rz, 2)) return
        if (any(shape(bsin) /= shape(bcos))) return
        allocate (cp(np), sp(np), ck(np, 0:ntor), sk(np, 0:ntor))
        do j = 1, np
            cp(j) = cos(2*pi*(j - 1)/np); sp(j) = sin(2*pi*(j - 1)/np)
            do k = 0, ntor
                nk = mod(k*nfp*(j - 1), np)
                ck(j, k) = merge(1.0_dp, 2.0_dp, k == 0)*cos(2*pi*nk/np)/np
                sk(j, k) = 2*sin(2*pi*nk/np)/np
            end do
        end do
        !$omp parallel do default(none) schedule(dynamic, 1) &
        !$omp shared(self, rz, ntor, np, cp, sp, ck, sk, bcos, bsin) &
        !$omp private(q, j, c, s, k, p, bc, bcyl, bk)
        do q = 1, size(rz, 2)
            bcos(:, :, :, q) = 0; bsin(:, :, :, q) = 0
            do j = 1, np
                p = [rz(1, q)*cp(j), rz(1, q)*sp(j), rz(2, q)]
                bk = 0
                do c = 1, self%coils%ncoil
                    if (self%circuit(c) == 0) cycle
                    bc = 0
                    do s = self%coils%first(c), self%coils%first(c + 1) - 1
                        bc = bc + segment_b(self%coils%a(:, s), self%coils%l(:, s), p)
                    end do
                    bk(:, self%circuit(c)) = bk(:, self%circuit(c)) + self%turns(c)*bc
                end do
                do c = 1, self%ncircuit
                    bcyl = self%field_scale*[bk(1, c)*cp(j) + bk(2, c)*sp(j), &
                        -bk(1, c)*sp(j) + bk(2, c)*cp(j), bk(3, c)]
                    do k = 0, ntor
                        bcos(:, k, c, q) = bcos(:, k, c, q) + ck(j, k)*bcyl
                        bsin(:, k, c, q) = bsin(:, k, c, q) + sk(j, k)*bcyl
                    end do
                end do
            end do
            bsin(:, 0, :, q) = 0
        end do
        !$omp end parallel do
        status = IDEAL_FB_OK
    end subroutine

    subroutine coil_sample_points(self, x, nfp, ntor, b, status)
        !! Unit-current circuit fields at Cartesian points (Cartesian
        !! components), no toroidal sampling. Points are independent (OpenMP
        !! over points), so the result does not depend on the thread count.
        class(coil_circuit_basis_t), intent(in) :: self
        real(dp), intent(in) :: x(:, :)
        integer, intent(in) :: nfp, ntor
        real(dp), intent(out) :: b(:, :, :)
        integer, intent(out) :: status
        real(dp) :: bc(3)
        integer :: q, c, s
        status = IDEAL_FB_BAD_INPUT
        if (size(x, 1) /= 3 .or. nfp < 1 .or. ntor < 0) return
        if (any(shape(b) /= [3, self%ncircuit, size(x, 2)])) return
        !$omp parallel do default(none) schedule(dynamic, 16) &
        !$omp shared(self, x, b) private(q, c, s, bc)
        do q = 1, size(x, 2)
            b(:, :, q) = 0
            do c = 1, self%coils%ncoil
                if (self%circuit(c) == 0) cycle
                bc = 0
                do s = self%coils%first(c), self%coils%first(c + 1) - 1
                    bc = bc + segment_b(self%coils%a(:, s), self%coils%l(:, s), x(:, q))
                end do
                b(:, self%circuit(c), q) = b(:, self%circuit(c), q) + self%turns(c)*bc
            end do
            b(:, :, q) = self%field_scale*b(:, :, q)
        end do
        !$omp end parallel do
        status = IDEAL_FB_OK
    end subroutine

    subroutine ideal_fb_circuits_build(fb, basis, circuits, status, background, &
                                       nphi_surface)
        !! Per-circuit traces (zero circulation: circuits lie outside Gamma_c)
        !! and the fixed part from the optional background sampler, which
        !! carries the unit-circulation filament times fb%circulation as in
        !! ideal_fb_set_drive. Without a background the fixed part is zero.
        !! With nphi_surface the traces are taken on the coupling surface of
        !! the plasma chart through the Cartesian basis contract
        !! (ideal_fb_drive_traces_3d, nphi_surface toroidal points), which
        !! also holds for a non-axisymmetric Gamma_c.
        type(ideal_fb_t), intent(inout) :: fb
        class(ideal_fb_field_basis_t), intent(in) :: basis
        type(ideal_fb_circuits_t), intent(out) :: circuits
        integer, intent(out) :: status
        procedure(ideal_fb_sampler), optional :: background
        integer, intent(in), optional :: nphi_surface
        real(dp), allocatable :: zeros(:)
        integer :: nc
        status = IDEAL_FB_BAD_INPUT
        nc = basis%nfield()
        if (nc < 1) return
        circuits%ncircuit = nc
        allocate (circuits%drive(fb%nt, fb%nmode, nc), &
            circuits%load(fb%nt, fb%nmode, nc), &
            circuits%f_vac(nc), circuits%closure(nc), circuits%f_spread(nc), zeros(nc))
        zeros = 0
        if (present(nphi_surface)) then
            call ideal_fb_drive_traces_3d(fb, basis, zeros, nphi_surface, &
                circuits%drive, circuits%load, circuits%closure, circuits%f_vac, &
                circuits%f_spread, status)
        else
            call ideal_fb_drive_traces(fb, basis, zeros, circuits%drive, &
                circuits%load, circuits%closure, circuits%f_vac, circuits%f_spread, &
                status)
        end if
        if (status /= IDEAL_FB_OK) return
        allocate (circuits%drive0(fb%nt, fb%nmode), circuits%load0(fb%nt, fb%nmode))
        circuits%drive0 = 0; circuits%load0 = 0; circuits%f_vac0 = 0
        if (present(background)) then
            ! The fixed part uses ideal_fb_set_drive itself, then restores fb.
            block
                real(dp), allocatable :: d(:, :), l(:, :)
                integer :: st
                d = fb%drive; l = fb%load
                if (present(nphi_surface)) then
                    call ideal_fb_set_drive_3d(fb, background, nphi_surface, status)
                else
                    call ideal_fb_set_drive(fb, background, status)
                end if
                if (status /= IDEAL_FB_OK) return
                circuits%drive0 = fb%drive; circuits%load0 = fb%load
                call ideal_fb_set_traces(fb, d, l, fb%f_vac, st)
            end block
        end if
        status = IDEAL_FB_OK
    end subroutine

    subroutine ideal_fb_set_currents(fb, circuits, currents, status)
        !! Install a(p), h(p), F_vac(p) for the circuit currents p (A, or the
        !! consumer's current unit consistent with the basis).
        type(ideal_fb_t), intent(inout) :: fb
        type(ideal_fb_circuits_t), intent(in) :: circuits
        real(dp), intent(in) :: currents(:)
        integer, intent(out) :: status
        real(dp), allocatable :: d(:, :), l(:, :)
        real(dp) :: f
        integer :: k
        status = IDEAL_FB_BAD_INPUT
        if (size(currents) /= circuits%ncircuit) return
        if (.not. all(ieee_is_finite(currents))) return
        d = circuits%drive0; l = circuits%load0; f = circuits%f_vac0
        do k = 1, circuits%ncircuit
            d = d + currents(k)*circuits%drive(:, :, k)
            l = l + currents(k)*circuits%load(:, :, k)
            f = f + currents(k)*circuits%f_vac(k)
        end do
        call ideal_fb_set_traces(fb, d, l, f, status)
    end subroutine

    subroutine ideal_fb_currents_jvp(fb, circuits, dcurrents, jv, status)
        !! jv = (dR/dp) dp at the current state of fb.
        type(ideal_fb_t), intent(in) :: fb
        type(ideal_fb_circuits_t), intent(in) :: circuits
        real(dp), intent(in) :: dcurrents(:)
        real(dp), intent(out) :: jv(:)
        integer, intent(out) :: status
        real(dp), allocatable :: d(:, :), l(:, :)
        real(dp) :: f
        integer :: k
        status = IDEAL_FB_BAD_INPUT
        if (size(dcurrents) /= circuits%ncircuit) return
        allocate (d(fb%nt, fb%nmode), l(fb%nt, fb%nmode))
        d = 0; l = 0; f = 0
        do k = 1, circuits%ncircuit
            d = d + dcurrents(k)*circuits%drive(:, :, k)
            l = l + dcurrents(k)*circuits%load(:, :, k)
            f = f + dcurrents(k)*circuits%f_vac(k)
        end do
        call ideal_fb_traces_jvp(fb, d, l, f, jv, status)
    end subroutine

    subroutine ideal_fb_currents_vjp(fb, circuits, v, jtv, status)
        !! jtv = (dR/dp)^T v at the current state of fb.
        type(ideal_fb_t), intent(in) :: fb
        type(ideal_fb_circuits_t), intent(in) :: circuits
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: jtv(:)
        integer, intent(out) :: status
        real(dp), allocatable :: d(:, :), l(:, :)
        real(dp) :: f
        integer :: k
        status = IDEAL_FB_BAD_INPUT
        if (size(jtv) /= circuits%ncircuit) return
        allocate (d(fb%nt, fb%nmode), l(fb%nt, fb%nmode))
        call ideal_fb_traces_vjp(fb, v, d, l, f, status)
        if (status /= IDEAL_FB_OK) return
        !$omp parallel do default(shared) private(k) schedule(static)
        do k = 1, circuits%ncircuit
            jtv(k) = sum(d*circuits%drive(:, :, k)) + sum(l*circuits%load(:, :, k)) &
                + f*circuits%f_vac(k)
        end do
        !$omp end parallel do
    end subroutine
end module kin6d_ideal_fb_circuits
