!> Thermal neoclassical transport of electrons and one ion species on a 3D
!! surface from monoenergetic coefficients with the momentum correction of
!! Sugama and Nishimura, Phys. Plasmas 9, 4637 (2002) (SN): a model separate
!! from the full linearised operator of ntv_surface_coefficients, cheaper
!! because every speed node is an independent monoenergetic solve.
!!
!! Per species a and speed node x_k (v = v_a x_k), mono_coefficients_3d at
!! nu = nu_D^a(x_k)/v gives d11, d13, d33; the leaf dk_mono_thermal_point
!! turns them into SN's monoenergetic M, N, L (Eq. 43), the energy integrals
!! M_aj, N_aj, L_aj follow from Eq. (36) on the nodes of the speed grid, and
!! Eqs. (C8)-(C13) give the transport matrix (C5) with the 13-moment friction
!! of electrons and ions to lowest order in (m_e/m_i)^(1/2). SI units, radial
!! label the poloidal flux per radian, the DKES source of SN Eq. (40) taken as
!! c_s v^2 times the radial drift of dk_mono with c_s = 2 G m_a/(Z_a e B0)
!! (G the covariant toroidal Boozer component), the scale for which SN's
!! symmetric-case identity (D4) holds; forces X_a1 = -p_a'/n_a - Z_a e Phi',
!! X_a2 = -T_a', X_E = <B E_par>/<B^2>^(1/2).
module milca_dk_mono_thermal
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_dk_mono, only: mono_coefficients_t, mono_coefficients_3d, &
        mono_surface_operators, MONO_OK
    use milca_dk_mono_leaves, only: dk_mono_thermal_point
    use milca_dk_collop, only: collop_species_t, collop_grid_t, collop_nuhat, &
        collop_vth, collop_qe
    implicit none
    private
    public :: mono_thermal_t, mono_thermal_coefficients
    integer, parameter, public :: THERMAL_OK = 0, THERMAL_BAD_INPUT = 1, &
        THERMAL_SINGULAR = 2

    !> l: rows (Gamma_e, q_e/T_e, Gamma_i, q_i/T_i, J_E^BS), columns
    !> (X_e1, X_e2, X_i1, X_i2, X_E), Eq. (C5); sigma_s the 13-moment Spitzer
    !> conductivity (C7), so <j_par B>/<B E_par> = sigma_s + l(5, 5).
    !> mnl(j, k, a): M_aj, N_aj, L_aj (k = 1, 2, 3) of Eq. (36), species
    !> a = 1 electrons, 2 ions; pas(j, a) the energy integral of D11/T alone
    !> (pitch-angle scattering without the momentum correction), j = 1..3.
    type :: mono_thermal_t
        real(dp) :: l(5, 5) = 0, sigma_s = 0, b2 = 0
        real(dp) :: mnl(3, 3, 2) = 0, pas(3, 2) = 0
    end type

    real(dp), parameter :: pi = 4*atan(1.0_dp)

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface

contains

    !> species(1) electrons (z = -1), species(2) the ion species; b, nper, q,
    !> i_over_g, kappa as in mono_coefficients_3d, b0 the field unit [T],
    !> lnlambda(a, b) the Coulomb logarithms, n_leg Legendre polynomials.
    subroutine mono_thermal_coefficients(b, nper, q, i_over_g, kappa, b0, species, &
            lnlambda, grid, n_leg, res, status)
        real(dp), intent(in) :: b(:, :), q, i_over_g, kappa, b0
        integer, intent(in) :: nper, n_leg
        type(collop_species_t), intent(in) :: species(:)
        real(dp), intent(in) :: lnlambda(:, :)
        type(collop_grid_t), intent(in) :: grid
        type(mono_thermal_t), intent(out) :: res
        integer, intent(out) :: status
        real(dp) :: tau(2), ups2, g_cov, w(size(b))
        integer :: a

        status = THERMAL_BAD_INPUT
        if (size(species) /= 2 .or. any(shape(lnlambda) /= [2, 2])) return
        if (species(1)%z >= 0 .or. species(2)%z <= 0 .or. b0 <= 0) return
        w = reshape(1/b**2, [size(b)])
        res%b2 = sum(w*reshape(b**2, [size(b)]))/sum(w)
        call pfirsch_schlueter(b, nper, q, i_over_g, ups2, status)
        if (status /= THERMAL_OK) return
        ! kappa = iota B0/(G + iota I).
        g_cov = b0/(kappa*q*(1 + i_over_g/q))
        do a = 1, 2
            call energy_integrals(a, res%mnl(:, :, a), res%pas(:, a), status)
            if (status /= THERMAL_OK) return
            tau(a) = 3*sqrt(pi)/(4*collop_nuhat(species(a), species(a), lnlambda(a, a)))
        end do
        call assemble(res, species, tau, b0, status)
    contains
        !> M_aj, N_aj, L_aj of species a on the nodes 2..n of the speed grid.
        subroutine energy_integrals(a, mnl, pas, status)
            integer, intent(in) :: a
            real(dp), intent(out) :: mnl(3, 3), pas(3)
            integer, intent(out) :: status
            type(mono_coefficients_t) :: c
            real(dp) :: x, v, nud, nu, cs, vm, vn, vl, weight, y, gch
            integer :: k, j, bsp, st
            mnl = 0
            pas = 0
            status = THERMAL_SINGULAR
            ! Drift scale of SN's DKES source sigma_1 in units of the radial
            ! drift of dk_mono, fixed by SN's symmetric-case identity N/M =
            ! B_zeta/(e chi' <B^2>) (Eqs. (46), (D4)); test dk_mono_thermal_
            ! symmetric checks it on an axisymmetric and a quasi-helical field.
            cs = 2*g_cov*species(a)%m/(species(a)%z*collop_qe*b0)
            do k = 2, grid%n
                x = grid%x(k)
                v = collop_vth(species(a))*x
                nud = 0
                do bsp = 1, 2
                    y = x*collop_vth(species(a))/collop_vth(species(bsp))
                    gch = (erf(y) - 2*y*exp(-y*y)/sqrt(pi))/(2*y*y)
                    nud = nud + collop_nuhat(species(a), species(bsp), lnlambda(a, bsp)) &
                        *(erf(y) - gch)/x**3
                end do
                nu = nud/v
                call mono_coefficients_3d(b, nper, q, i_over_g, kappa, nu, n_leg, c, st)
                if (st /= MONO_OK) return
                call dk_mono_thermal_point(species(a)%m, species(a)%t, nud, v, nu, b0, &
                    res%b2, cs, c%d11, c%d13, c%d33, ups2, vm, vn, vl)
                ! n (2/sqrt(pi)) int dK sqrt(K) exp(-K) (K - 5/2)^(j - 1), dK sqrt(K)
                ! = 2 x^2 dx.
                do j = 1, 3
                    weight = species(a)%n*2/sqrt(pi)*grid%w(k)*2*x*x*exp(-x*x) &
                        *(x*x - 2.5_dp)**(j - 1)
                    mnl(j, 1) = mnl(j, 1) + weight*vm
                    mnl(j, 2) = mnl(j, 2) + weight*vn
                    mnl(j, 3) = mnl(j, 3) + weight*vl
                    pas(j) = pas(j) + weight*cs*cs*v**3*c%d11/(4*species(a)%t)
                end do
            end do
            status = THERMAL_OK
        end subroutine energy_integrals
    end subroutine mono_thermal_coefficients

    !> <u^2> with b^2 d_par(u/b) = -b_drift/b and <b u> = 0 (SN Eq. (12) in
    !> the units of dk_mono), the Pfirsch-Schlueter function of L.
    subroutine pfirsch_schlueter(b, nper, q, i_over_g, ups2, status)
        real(dp), intent(in) :: b(:, :), q, i_over_g
        integer, intent(in) :: nper
        real(dp), intent(out) :: ups2
        integer, intent(out) :: status
        real(dp), allocatable :: d(:, :), b_par(:), b_drift(:), bv(:), sys(:, :), &
            rhs(:, :), wgt(:), u(:)
        integer, allocatable :: piv(:)
        integer :: n, st
        status = THERMAL_BAD_INPUT
        call mono_surface_operators(b, nper, q, i_over_g, d, b_par, b_drift, st)
        if (st /= MONO_OK) return
        n = size(b)
        bv = reshape(b, [n])
        ! Bordered by a constant source and sum(w/b) = 0, w = u/b.
        allocate (sys(n + 1, n + 1), rhs(n + 1, 1), piv(n + 1))
        sys = 0
        sys(1:n, 1:n) = d
        sys(1:n, n + 1) = 1
        sys(n + 1, 1:n) = 1
        rhs(1:n, 1) = -b_drift/bv**3
        rhs(n + 1, 1) = 0
        status = THERMAL_SINGULAR
        call dgesv(n + 1, 1, sys, n + 1, piv, rhs, n + 1, st)
        if (st /= 0) return
        wgt = 1/bv**2
        wgt = wgt/sum(wgt)
        u = bv*rhs(1:n, 1)
        u = u - bv*sum(wgt*bv*u)/sum(wgt*bv*bv)
        ups2 = sum(wgt*u*u)
        status = 0
    end subroutine pfirsch_schlueter

    !> Eqs. (C7)-(C13) for electrons (1) and ions (2).
    subroutine assemble(res, species, tau, b0, status)
        type(mono_thermal_t), intent(inout) :: res
        type(collop_species_t), intent(in) :: species(:)
        real(dp), intent(in) :: tau(2), b0
        integer, intent(out) :: status
        real(dp) :: cl(2, 2, 2), cm(2, 2, 2), cn(2, 2, 2), lam(2, 2, 2)
        real(dp), dimension(2, 2) :: e11, ae, ai, aei, ki, block
        real(dp) :: bsq, z, ne, me, e, l11, l12, l22, row(2)
        integer :: a
        status = THERMAL_SINGULAR
        bsq = b0*b0*res%b2
        z = species(2)%z
        ne = species(1)%n
        me = species(1)%m
        e = collop_qe
        l11 = z
        l12 = 1.5_dp*z
        l22 = sqrt(2.0_dp) + 3.25_dp*z
        lam(:, :, 1) = reshape([l11, -l12, -l12, l22], [2, 2])
        lam(:, :, 2) = reshape([0.0_dp, 0.0_dp, 0.0_dp, sqrt(2.0_dp)], [2, 2])
        e11 = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2])
        do a = 1, 2
            associate (mnl => res%mnl(:, :, a), nm => species(a)%n*species(a)%m)
                cl(:, :, a) = reshape([mnl(1, 3), mnl(2, 3), mnl(2, 3), mnl(3, 3)], [2, 2])
                cm(:, :, a) = tau(a)/(nm*bsq)*reshape([mnl(1, 1), mnl(2, 1), mnl(2, 1), &
                    mnl(3, 1)], [2, 2])
                cn(:, :, a) = tau(a)/nm*reshape([mnl(1, 2), mnl(2, 2), mnl(2, 2), mnl(3, 2)], &
                    [2, 2])
            end associate
        end do
        if (.not. inverse(cm(:, :, 1) + lam(:, :, 1), ae)) return
        if (.not. inverse(cm(:, :, 2) + lam(:, :, 2), ai)) return
        ! (M_e^-1 + Lambda_e^-1)^-1 = M_e (M_e + Lambda_e)^-1 Lambda_e.
        aei = matmul(cm(:, :, 1), matmul(ae, lam(:, :, 1)))
        res%l = 0
        block = cl(:, :, 1) - species(1)%n*species(1)%m/(tau(1)*bsq) &
            *matmul(cn(:, :, 1), matmul(ae, cn(:, :, 1)))
        res%l(1:2, 1:2) = block
        ki = matmul(cn(:, :, 2), ai)
        block = cl(:, :, 2) - species(2)%n*species(2)%m/(tau(2)*bsq)*matmul(ki, cn(:, :, 2)) &
            + ne*me/(tau(1)*bsq)*matmul(ki, matmul(e11, matmul(aei, matmul(e11, &
            matmul(ai, cn(:, :, 2))))))
        res%l(3:4, 3:4) = block
        block = -ne*me/(tau(1)*bsq)*matmul(cn(:, :, 1), matmul(ae, matmul(lam(:, :, 1), matmul(e11, &
            matmul(ai, cn(:, :, 2))))))
        res%l(1:2, 3:4) = block
        res%l(3:4, 1:2) = transpose(block)
        row = ne*e/sqrt(bsq)*matmul([1.0_dp, 0.0_dp], matmul(ae, cn(:, :, 1)))
        res%l(5, 1:2) = row
        res%l(1:2, 5) = -row
        row = -ne*e/sqrt(bsq)*matmul([1.0_dp, 0.0_dp], matmul(ae, matmul(cm(:, :, 1), &
            matmul(e11, matmul(ai, cn(:, :, 2))))))
        res%l(5, 3:4) = row
        res%l(3:4, 5) = -row
        ! L_EE = -(n_e e^2 tau_ee/m_e) [Lambda_e^-1 - (M_e + Lambda_e)^-1]_11.
        block = ne*e*e*tau(1)/me*(ae - inverse_of(lam(:, :, 1)))
        res%l(5, 5) = block(1, 1)
        res%sigma_s = ne*e*e*tau(1)/me*l22/(l11*l22 - l12*l12)
        status = THERMAL_OK
    end subroutine assemble

    logical function inverse(m, inv) result(ok)
        real(dp), intent(in) :: m(2, 2)
        real(dp), intent(out) :: inv(2, 2)
        real(dp) :: det
        det = m(1, 1)*m(2, 2) - m(1, 2)*m(2, 1)
        ok = abs(det) > tiny(det)
        inv = 0
        if (.not. ok) return
        inv = reshape([m(2, 2), -m(2, 1), -m(1, 2), m(1, 1)], [2, 2])/det
    end function inverse

    function inverse_of(m) result(inv)
        real(dp), intent(in) :: m(2, 2)
        real(dp) :: inv(2, 2)
        logical :: ok
        ok = inverse(m, inv)
    end function inverse_of
end module milca_dk_mono_thermal
