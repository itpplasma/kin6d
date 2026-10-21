module kin6d_map_3d_equilibrium
    ! Cylindrical Fourier maps -> the common physical-field evaluator.
    ! Retain the frozen form; never Fourier-fit B from a sampled candidate.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_equilibrium_state, only: equilibrium_field_t, equilibrium_state_t, &
                      EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI, EQ_P, EQ_PSI, EQ_PSI_R, EQ_PSI_Z
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_ideal_map, only: ideal_map_t, ideal_frozen_t
    use kin6d_triangle_basis, only: triangle_basis_jet
    use kin6d_ideal_map_leaves, only: ideal_reference_form
    use kin6d_map_3d_readback_leaves, only: map_3d_field, map_3d_inverse_step, &
                                            map_3d_determinant
    implicit none
    private
    public :: map_3d_field_t, map_3d_equilibrium, map_3d_coordinates, map_3d_sample
    type, extends(equilibrium_field_t) :: map_3d_field_t
        private
        type(ideal_map_t) :: map
        class(ideal_frozen_t), allocatable :: frozen
    contains
        procedure :: sample => sample_many
    end type
contains
    subroutine map_3d_equilibrium(map, frozen, equilibrium, status, toroidal_flux)
        type(ideal_map_t), intent(in) :: map
        class(ideal_frozen_t), intent(in) :: frozen
        type(equilibrium_state_t), intent(out) :: equilibrium
        integer, intent(out) :: status
        real(dp), intent(in), optional :: toroidal_flux
        type(map_3d_field_t), allocatable :: field
        real(dp) :: y(2), dg(2, 2), chi, cg(2), b3, p, b(3), mismatch, scale
        integer :: e, g
        status = 1
        if (.not. allocated(map%chart)) return
        select type (chart => map%chart)
        type is (cylindrical_chart_t)
        class default
            return
        end select
        if (.not. allocated(map%x)) return
        if (.not. allocated(map%b)) return
        ! A caller must supply the very form frozen into the candidate.
        mismatch = 0; scale = max(tiny(1.0_dp), maxval(abs(map%b)), maxval(abs(map%p)))
        !$omp parallel do schedule(static) private(e,g,y,dg,chi,cg,b3,p,b) &
        !$omp reduction(max:mismatch)
        do e = 1, map%cells
            do g = 1, map%nq
                y = matmul(map%reference(:, map%conn(:, e)), map%shape(:, g))
                dg = matmul(map%reference(:, map%conn(:, e)), &
                    transpose(map%dshape(:, :, g)))
                call frozen%evaluate(y, chi, cg, b3, p)
                call ideal_reference_form(dg(1, 1), dg(2, 1), dg(1, 2), dg(2, 2), &
                                          cg(1), cg(2), b3, b(1), b(2), b(3))
                if (.not. all(ieee_is_finite([b,p]))) mismatch=huge(1.0_dp)
                mismatch = max(mismatch, maxval(abs(b - map%b(:, g, e))), &
                    abs(p - map%p(g, e)))
            end do
        end do
        !$omp end parallel do
        if (mismatch > 128*epsilon(1.0_dp)*scale) return
        allocate (field)
        field%map = map
        allocate (field%frozen, source=frozen)
        allocate (equilibrium%chart, source=map%chart)
        call move_alloc(field, equilibrium%field)
        equilibrium%mu0 = map%mu0
        equilibrium%lane = 'map'; equilibrium%claim = 'candidate'
        equilibrium%seed = 'represented Fourier map and frozen reference form'
        if (present(toroidal_flux)) equilibrium%phi_tor = toroidal_flux
        status = 0
    end subroutine

    subroutine map_3d_coordinates(map, cell, xi, phi, q, derivative, reference, dg)
        type(ideal_map_t), intent(in) :: map
        integer, intent(in) :: cell
        real(dp), intent(in) :: xi(2), phi
        real(dp), intent(out) :: q(3), derivative(3, 3), reference(2), dg(2, 2)
        real(dp) :: shape(10), gradient(2, 10), mode, slope, coefficients(3, 10)
        real(dp) :: phase
        integer :: k, n
        n = map%nloc
        call triangle_basis_jet(xi, shape(:n), gradient(:, :n))
        q = 0; derivative = 0
        do k = 1, map%nmode
            mode = 1; slope = 0
            if (k > 1) then
                phase = (k/2)*map%nfp*phi
                if (mod(k, 2) == 0) then
                    mode = cos(phase); slope = -(k/2)*map%nfp*sin(phase)
                else
                    mode = sin(phase); slope = (k/2)*map%nfp*cos(phase)
                end if
            end if
            coefficients(:, :n) = map%x(:, k, map%conn(:, cell))
            q = q + mode*matmul(coefficients(:, :n), shape(:n))
            derivative(:, :2) = derivative(:, :2) &
                          + mode*matmul(coefficients(:, :n), transpose(gradient(:, :n)))
            derivative(:, 3) = derivative(:, 3) &
                               + slope*matmul(coefficients(:, :n), shape(:n))
        end do
        q(3) = q(3) + phi; derivative(3, 3) = derivative(3, 3) + 1
        reference = matmul(map%reference(:, map%conn(:, cell)), shape(:n))
        dg = matmul(map%reference(:, map%conn(:, cell)), transpose(gradient(:, :n)))
    end subroutine

    subroutine map_3d_sample(map, frozen, cell, xi, phi, point, value, status)
        type(ideal_map_t), intent(in) :: map
        class(ideal_frozen_t), intent(in) :: frozen
        integer, intent(in) :: cell
        real(dp), intent(in) :: xi(2), phi
        real(dp), intent(out) :: point(3), value(EQ_VALUES)
        integer, intent(out) :: status
        real(dp) :: q(3), d(3, 3), reference(2), g(2, 2), chi, cg(2), b3, p
        real(dp) :: jac, chi_phi, det
        call map_3d_coordinates(map, cell, xi, phi, q, d, reference, g)
        point = q
        value = 0; status = 1
        call map_3d_determinant(d(1, 1), d(2, 1), d(3, 1), d(1, 2), d(2, 2), &
                                d(3, 2), d(1, 3), d(2, 3), d(3, 3), det)
        if (.not. -map%sigma*q(1)*det > 0) return
        call frozen%evaluate(reference, chi, cg, b3, p)
        call map_3d_field(d(1, 1), d(2, 1), d(3, 1), d(1, 2), d(2, 2), &
                          d(3, 2), d(1, 3), d(2, 3), d(3, 3), g(1, 1), g(2, 1), &
                          g(1, 2), g(2, 2), cg(1), cg(2), b3, q(1), map%sigma, &
                          value(EQ_BR), value(EQ_BZ), value(EQ_BPHI), &
                          value(EQ_PSI_R), value(EQ_PSI_Z), chi_phi, jac)
        ! In (R,Z,phi), the flux gauge BR=-psi_Z/R has psi=-sigma*chi.
        value(EQ_PSI_R:EQ_PSI_Z) = -map%sigma*value(EQ_PSI_R:EQ_PSI_Z)
        value(EQ_PSI) = -map%sigma*chi
        value(EQ_P) = p
        status = 1
        if (jac <= 0 .or. q(1) <= 0) return
        if (.not. all(ieee_is_finite(value))) return
        status = 0
    end subroutine

    subroutine sample_many(field, points, values, statuses)
        class(map_3d_field_t), intent(in) :: field
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        real(dp) :: target(3), xi(2), phi, point(3)
        integer :: i, cell, seed, local
        values = 0; statuses = 1
        if (size(points, 1) < 2 .or. size(points, 1) > 3) return
        if (size(values, 1) /= EQ_VALUES) return
        if (size(points, 2) /= size(values, 2)) return
        if (size(statuses) /= size(points, 2)) return
        seed = 1
        !$omp parallel do schedule(static) firstprivate(seed) &
        !$omp private(i,target,xi,phi,point,cell,local)
        do i = 1, size(points, 2)
            target = 0; target(:size(points, 1)) = points(:, i)
            if (.not. all(ieee_is_finite(target))) cycle
            if (target(1) <= 0) cycle
            call locate(field%map, target, seed, cell, xi, phi, local)
            if (local /= 0) cycle
            seed = cell
            call map_3d_sample(field%map, field%frozen, cell, xi, phi, point, &
                               values(:, i), statuses(i))
        end do
        !$omp end parallel do
    end subroutine

    subroutine locate(map, target, seed, cell, xi, phi, status)
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(in) :: target(3)
        integer, intent(in) :: seed
        integer, intent(out) :: cell, status
        real(dp), intent(out) :: xi(2), phi
        real(dp) :: q(3), d(3, 3), reference(2), g(2, 2), residual(3), step(3)
        real(dp) :: damping, scale, pi, det
        integer :: attempt, iteration
        status = 1; cell = 0; pi = acos(-1.0_dp)
        scale = max(1.0_dp, maxval(abs(target(:2))))
        do attempt = 0, map%cells - 1
            cell = 1 + modulo(seed - 1 + merge(attempt/2, &
                          map%cells - (attempt + 1)/2, mod(attempt, 2) == 0), map%cells)
            xi = 1.0_dp/3; phi = target(3)
            do iteration = 1, 20
                call map_3d_coordinates(map, cell, xi, phi, q, d, reference, g)
                residual = q - target
                residual(3) = modulo(residual(3) + pi, 2*pi) - pi
                if (maxval(abs(residual)) < 5e-13_dp*scale) then
                    if (minval([xi, 1 - sum(xi)]) >= -5e-10_dp) then
                        status = 0; return
                    end if
                    exit
                end if
                call map_3d_determinant(d(1, 1), d(2, 1), d(3, 1), d(1, 2), d(2, 2), &
                                        d(3, 2), d(1, 3), d(2, 3), d(3, 3), det)
                if (.not. abs(det) > tiny(det)) exit
                call map_3d_inverse_step(d(1, 1), d(2, 1), d(3, 1), &
                                 d(1, 2), d(2, 2), d(3, 2), d(1, 3), d(2, 3), d(3, 3), &
                       residual(1), residual(2), residual(3), step(1), step(2), step(3))
                if (.not. all(ieee_is_finite(step))) exit
                damping = max(1.0_dp, maxval(abs(step)))
                xi = xi - step(:2)/damping; phi = phi - step(3)/damping
                if (minval(xi) < -2 .or. maxval(xi) > 3) exit
            end do
        end do
        cell = 0
    end subroutine
end module
