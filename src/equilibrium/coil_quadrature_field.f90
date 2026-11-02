module kin6d_coil_quadrature_field
    ! Exact polygon source, evaluated by FortNum Gauss-Legendre quadrature.
    ! Quadrature order is input data and must be checked by refinement.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_coils_3d, only: coil_set_t, coil_mu0
    use kin6d_equilibrium_state, only: equilibrium_field_t, EQ_BR, EQ_BZ, EQ_BPHI
    use kin6d_island_leaves, only: biot_savart_integrand
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    private
    public :: coil_quadrature_field_t, coil_quadrature_init
    type, extends(equilibrium_field_t) :: coil_quadrature_field_t
        real(dp), allocatable :: node(:, :), dl(:, :)
        real(dp) :: ftor = 0
    contains
        procedure :: sample => quadrature_sample
    end type
contains
    subroutine coil_quadrature_init(field, coils, currents, order, ftor, status)
        type(coil_quadrature_field_t), intent(out) :: field
        type(coil_set_t), intent(in) :: coils
        real(dp), intent(in) :: currents(:), ftor
        integer, intent(in) :: order
        integer, intent(out) :: status
        real(dp), allocatable :: nodes(:), weights(:)
        real(dp) :: scale
        integer :: k, s, q, j
        status = 1
        if (order < 1 .or. size(currents) /= coils%ncoil) return
        allocate (nodes(order), weights(order))
        call gauss_legendre(order, nodes, weights)
        allocate (field%node(3, order*coils%nseg), field%dl(3, order*coils%nseg))
        do k = 1, coils%ncoil
            scale = coil_mu0*currents(k)/(8*acos(-1.0_dp))
            do s = coils%first(k), coils%first(k + 1) - 1
                do q = 1, order
                    j = (s - 1)*order + q
                    field%node(:, j) = coils%a(:, s) + &
                                       (nodes(q) + 1)*coils%l(:, s)/2
                    field%dl(:, j) = scale*weights(q)*coils%l(:, s)
                end do
            end do
        end do
        field%ftor = ftor; status = 0
    end subroutine

    subroutine quadrature_sample(field, points, values, statuses)
        class(coil_quadrature_field_t), intent(in) :: field
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        real(dp) :: phi, p(3), d(3), cart(3), bx, by, bz
        integer :: i, j
        values = 0; statuses = 1
        if (.not. allocated(field%node)) return
        do i = 1, size(points, 2)
            if (.not. all(ieee_is_finite(points(:, i)))) cycle
            if (points(1, i) <= 0) cycle
            phi = 0
            if (size(points, 1) >= 3) phi = points(3, i)
            p = [points(1, i)*cos(phi), points(1, i)*sin(phi), points(2, i)]
            cart = 0
            do j = 1, size(field%node, 2)
                d = p - field%node(:, j)
                if (norm2(d) <= tiny(1.0_dp)) exit
                call biot_savart_integrand(d(1), d(2), d(3), &
                             field%dl(1, j), field%dl(2, j), field%dl(3, j), bx, by, bz)
                cart = cart + [bx, by, bz]
            end do
            if (j <= size(field%node, 2)) cycle
            values(EQ_BR, i) = cart(1)*cos(phi) + cart(2)*sin(phi)
            values(EQ_BPHI, i) = -cart(1)*sin(phi) + cart(2)*cos(phi) &
                                 + field%ftor/points(1, i)
            values(EQ_BZ, i) = cart(3)
            statuses(i) = 0
        end do
    end subroutine
end module
