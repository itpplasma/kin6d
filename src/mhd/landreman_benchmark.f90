module kin6d_landreman_benchmark
    ! Numerical benchmark data, mu0=1. No validated continuation/error bound.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_landreman_leaves, only: landreman_axisymmetric, landreman_surface, &
                                      landreman_current_ray, landreman_toroidal_external
    use kin6d_ideal_map, only: ideal_frozen_t
    use fortnum_quadrature, only: gauss_legendre
    use fortfem_axisymmetric_gs_kernel, only: gs_green_gradient
    implicit none
    private
    public :: landreman_frozen_t, landreman_drive_t, landreman_volume_boundary
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter, public :: LANDREMAN_DELTA = 0.01_dp
    type, extends(ideal_frozen_t) :: landreman_frozen_t
        real(dp) :: delta = LANDREMAN_DELTA
    contains
        procedure :: evaluate => frozen_evaluate
    end type
    type :: landreman_drive_t
        real(dp) :: delta = LANDREMAN_DELTA, circulation = 0, fit_error = huge(1.0_dp)
        real(dp), allocatable :: source(:, :), current(:)
    contains
        procedure :: fit => drive_fit
        procedure :: external => drive_external
        procedure :: sample => drive_sample
    end type
contains
    subroutine frozen_evaluate(field, point, chi, chi_gradient, b3, pressure)
        class(landreman_frozen_t), intent(in) :: field
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: chi, chi_gradient(2), b3, pressure
        real(dp) :: b(3), j(3), psi, gp(2)
        call landreman_axisymmetric(point(1), point(2), field%delta, psi, &
                             pressure, b(1), b(2), b(3), j(1), j(2), j(3), gp(1), gp(2))
        ! Map convention b_pol=(-chi_Z,chi_R); physical B_pol=b_pol/R.
        chi = -psi; chi_gradient = -gp; b3 = b(2)
    end subroutine

    subroutine landreman_volume_boundary(theta, delta, nq, bj, ip)
        ! Toroidal J volume integral, evaluated on S0. Boundary-centred
        ! Duffy rays remove the 1/distance singularity of the ring field.
        ! Poloidal J produces BJ_phi=(F(psi)-F(delta))/R inside and zero
        ! outside: its continuous boundary value is exactly zero.
        real(dp), intent(in) :: theta, delta
        integer, intent(in) :: nq
        real(dp), intent(out) :: bj(3), ip
        real(dp) :: gx(nq), gw(nq), r, z, dr, dz, rp, zp, wt, g, gg(4)
        real(dp) :: angle, rho, ww
        integer :: a, k
        call gauss_legendre(nq, gx, gw)
        call landreman_surface(theta, delta, r, z, dr, dz)
        bj = 0; ip = 0
        do a = 1, nq
            angle = (gx(a) + 1)*pi/2
            do k = 1, nq
                rho = (gx(k) + 1)/2
                call landreman_current_ray(theta, delta, rho, angle, rp, zp, wt)
                ww = gw(a)*gw(k)*pi*wt/4
                call gs_green_gradient(r, z, rp, zp, g, gg)
                bj(1) = bj(1) - ww*gg(2)/r
                bj(3) = bj(3) + ww*gg(1)/r
                ip = ip + ww
            end do
        end do
    end subroutine

    subroutine drive_fit(drive, ns, nf, nq, status)
        class(landreman_drive_t), intent(inout) :: drive
        integer, intent(in) :: ns, nf, nq
        integer, intent(out) :: status
        real(dp), allocatable :: a(:, :), rhs(:, :), singular(:), work(:)
        real(dp) :: theta, r, z, dr, dz, psi, p, b(3), j(3), gp(2), bj(3)
        real(dp) :: ip, g, gg(4), query(1)
        integer :: i, k, rank, lwork
        interface
            subroutine dgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, &
                              rank, work, lwork, info)
                import dp
                integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
                real(dp), intent(out) :: s(*), work(*)
                real(dp), intent(in) :: rcond
                integer, intent(out) :: rank, info
            end subroutine
        end interface
        status = 1
        if (ns < 4 .or. nf < ns .or. nq < 4) return
        if (.not. (drive%delta > 0 .and. 9*drive%delta < 0.25_dp)) return
        if (allocated(drive%source)) deallocate (drive%source)
        if (allocated(drive%current)) deallocate (drive%current)
        drive%fit_error = huge(1.0_dp)
        allocate (drive%source(2, ns), drive%current(ns), a(2*nf, ns), &
                  rhs(2*nf, 1), singular(ns))
        ! Sources on u^2+z^2=9 delta, outside Gamma_c (2.56 delta).
        do k = 1, ns
            call landreman_surface(2*pi*(k - 1)/ns, 9*drive%delta, &
                                   drive%source(1, k), drive%source(2, k), dr, dz)
        end do
        !$omp parallel do private(theta,r,z,dr,dz,psi,p,b,j,gp,bj,ip,g,gg,k)
        do i = 1, nf
            theta = 2*pi*(i - 1)/nf
            call landreman_surface(theta, drive%delta, r, z, dr, dz)
            call landreman_axisymmetric(r, z, drive%delta, psi, p, &
                                       b(1), b(2), b(3), j(1), j(2), j(3), gp(1), gp(2))
            call landreman_volume_boundary(theta, drive%delta, nq, bj, ip)
            rhs(2*i - 1, 1) = b(1) - bj(1); rhs(2*i, 1) = b(3) - bj(3)
            do k = 1, ns
                call gs_green_gradient(r, z, drive%source(1, k), &
                                       drive%source(2, k), g, gg)
                a(2*i - 1, k) = -gg(2)/r; a(2*i, k) = gg(1)/r
            end do
        end do
        !$omp end parallel do
        call landreman_volume_boundary(0.0_dp, drive%delta, nq, bj, ip)
        drive%circulation = -ip
        call dgelss(2*nf, ns, 1, a, 2*nf, rhs, 2*nf, singular, 1.0e-12_dp, &
                    rank, query, -1, status)
        if (status /= 0) return
        lwork = ceiling(query(1)); allocate (work(lwork))
        call dgelss(2*nf, ns, 1, a, 2*nf, rhs, 2*nf, singular, 1.0e-12_dp, &
                    rank, work, lwork, status)
        if (status /= 0) return
        drive%current = rhs(:ns, 1)
        drive%fit_error = 0
        do i = 1, nf
            theta = 2*pi*(i - 0.5_dp)/nf
            call landreman_surface(theta, drive%delta, r, z, dr, dz)
            call landreman_axisymmetric(r, z, drive%delta, psi, p, &
                                       b(1), b(2), b(3), j(1), j(2), j(3), gp(1), gp(2))
            call landreman_volume_boundary(theta, drive%delta, nq, bj, ip)
            call drive%external(r, z, j)
            drive%fit_error = max(drive%fit_error, norm2(bj + j - b))
        end do
    end subroutine

    subroutine drive_external(drive, r, z, b)
        class(landreman_drive_t), intent(in) :: drive
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: b(3)
        real(dp) :: g, gg(4)
        integer :: k
        b = 0
        call landreman_toroidal_external(r, drive%delta, b(2))
        do k = 1, size(drive%current)
            call gs_green_gradient(r, z, drive%source(1, k), drive%source(2, k), g, gg)
            b(1) = b(1) - drive%current(k)*gg(2)/r
            b(3) = b(3) + drive%current(k)*gg(1)/r
        end do
    end subroutine

    subroutine drive_sample(drive, r, z, k, bc, bs)
        ! Existing N4 sampler needs the external field PLUS its internal
        ! circulation filament; BJ's remaining exterior multipoles are solved
        ! by the BEM. Adding BJ itself here would double count plasma response.
        class(landreman_drive_t), intent(in) :: drive
        real(dp), intent(in) :: r, z
        integer, intent(in) :: k
        real(dp), intent(out) :: bc(3), bs(3)
        real(dp) :: g, gg(4)
        bc = 0; bs = 0
        if (k /= 0) return
        call drive%external(r, z, bc)
        call gs_green_gradient(r, z, 1.0_dp, 0.0_dp, g, gg)
        bc(1) = bc(1) + drive%circulation*gg(2)/r
        bc(3) = bc(3) - drive%circulation*gg(1)/r
    end subroutine
end module kin6d_landreman_benchmark
