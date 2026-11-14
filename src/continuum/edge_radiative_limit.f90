!> A uniform edge-shell electron-energy balance, not a universal Greenwald law.
!> A fold of the thermally stable hot branch is a model-dependent density limit.
module milca_edge_radiative_limit
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    real(dp), parameter :: elementary_charge = 1.602176634e-19_dp
    integer, parameter, public :: EDGE_OK = 0, EDGE_INVALID = 1, &
        EDGE_NO_HOT_BRANCH = 2, EDGE_NO_FOLD = 3
    type, public :: edge_model_t
        real(dp) :: j_parallel = 0.0_dp ! A / m^2 (local current density)
        real(dp) :: p_external = 0.0_dp ! W / m^3 deposited in shell
        real(dp) :: zeff = 1.0_dp
        real(dp) :: coulomb_log = 15.0_dp
        real(dp) :: chi_e = 1.0_dp ! m^2 / s, prescribed heat diffusivity
        real(dp) :: shell_width = 0.1_dp ! m, transport length
        real(dp) :: t_wall = 2.0_dp ! eV; reference sink temperature
        real(dp) :: impurity_fraction = 0.005_dp ! n_Z/n_e
        real(dp) :: cooling_peak = 1e-31_dp ! W m^3; USER FIT, not atomic data
        real(dp) :: cooling_tpeak = 15.0_dp ! eV
        real(dp) :: cooling_low_power = 2.0_dp
        real(dp) :: cooling_high_power = 1.0_dp
    end type
    public :: edge_equilibrium_density, edge_density_limit, &
        greenwald_reference_density, milca_edge_density_abi
contains
    !> Trace-impurity radiation: smooth two-power surrogate with its maximum
    !> L(Tp)=Lpeak; coefficients MUST be fitted to declared atomic data.
    pure subroutine cooling(model, t, l, dl)
        type(edge_model_t), intent(in) :: model
        real(dp), intent(in) :: t
        real(dp), intent(out) :: l, dl
        real(dp) :: x, a, b, u
        a = model%cooling_low_power
        b = model%cooling_high_power
        x = t/model%cooling_tpeak
        u = x**(a+b)
        l = model%cooling_peak*(a+b)*x**a/(b+a*u)
        dl = l*(a-a*(a+b)*u/(b+a*u))/t
    end subroutine

    pure logical function valid(model,t) result(ok)
        type(edge_model_t), intent(in) :: model
        real(dp), intent(in) :: t
        ok = all(ieee_is_finite([t,model%j_parallel,model%p_external, &
            model%zeff,model%coulomb_log,model%chi_e,model%shell_width, &
            model%t_wall,model%impurity_fraction,model%cooling_peak, &
            model%cooling_tpeak,model%cooling_low_power,model%cooling_high_power]))
        if (.not. ok) return
        ok = t > model%t_wall .and. t > 0 .and. model%t_wall >= 0 .and. &
            model%p_external >= 0 .and. model%zeff >= 1 .and. &
            model%coulomb_log > 0 .and. model%chi_e > 0 .and. &
            model%shell_width > 0 .and. model%impurity_fraction > 0 .and. &
            model%cooling_peak > 0 .and. model%cooling_tpeak > 0 .and. &
            model%cooling_low_power > 0 .and. model%cooling_high_power > 0
    end function

    !> Edge energy equation (at fixed n, fixed j):
    !> F(T,n)=P_ext + eta_Spitzer(T)*j^2 - n*D(T) - n^2*L_eff(T).
    !> D(T)=(3/2)*e*chi_e*(T-T_wall)/width^2, for electron heat only.
    !> L_eff=fz*Lz(T)+5.35e-37*sqrt(T_eV/1000) [W m^3].
    !> n_eq(T) solves F=0. Stability requires F_T(T,n_eq)<0.
    pure subroutine edge_equilibrium_density(model,t,n_eq,net_dT,status)
        type(edge_model_t), intent(in) :: model
        real(dp), intent(in) :: t
        real(dp), intent(out) :: n_eq, net_dT
        integer, intent(out) :: status
        real(dp) :: eta, qheat, heat, dheat, lz, dlz, lff, ll, dll
        real(dp) :: d, dd, discriminant
        n_eq = 0; net_dT = 0
        status = EDGE_INVALID
        if (.not. valid(model,t)) return
        eta = 5.2e-5_dp*model%zeff*model%coulomb_log/t**1.5_dp
        heat = eta*model%j_parallel**2
        dheat = -1.5_dp*heat/t
        qheat = model%p_external + heat
        if (qheat <= 0) return
        call cooling(model,t,lz,dlz)
        lff = 5.35e-37_dp*sqrt(t/1000.0_dp)
        ll = model%impurity_fraction*lz+lff
        dll = model%impurity_fraction*dlz+lff/(2*t)
        dd = 1.5_dp*elementary_charge*model%chi_e/model%shell_width**2
        d = dd*(t-model%t_wall)
        discriminant = sqrt(d*d+4*ll*qheat)
        n_eq = 2*qheat/(d+discriminant)
        net_dT = dheat - n_eq*dd - n_eq*n_eq*dll
        if (.not. all(ieee_is_finite([n_eq,net_dT]))) return
        if (n_eq <= 0) return
        status = EDGE_OK
    end subroutine

    !> First thermal fold reached by lowering T from the hot stable branch.
    !> Output is edge LOCAL density, NOT line-averaged density. A flat or
    !> always-stable scan has no detected fold, not an infinite Greenwald limit.
    pure subroutine edge_density_limit(model,t_min,t_max,n_limit,t_fold,status)
        type(edge_model_t), intent(in) :: model
        real(dp), intent(in) :: t_min,t_max
        real(dp), intent(out) :: n_limit,t_fold
        integer, intent(out) :: status
        real(dp) :: high,low,t,n,ft,base,root
        integer :: i,j,rc
        n_limit = 0; t_fold = 0
        status = EDGE_INVALID
        if (.not. all(ieee_is_finite([t_min,t_max]))) return
        if (t_min <= model%t_wall .or. t_min >= t_max) return
        call edge_equilibrium_density(model,t_max,n,ft,rc)
        if (rc /= EDGE_OK) return
        if (ft >= 0) then
            status = EDGE_NO_HOT_BRANCH
            return
        endif
        high = t_max
        base = log(t_max/t_min)/256
        do i = 1,256
            low = t_max*exp(-i*base)
            call edge_equilibrium_density(model,low,n,ft,rc)
            if (rc /= EDGE_OK) return
            if (ft >= 0) exit
            high = low
        enddo
        if (ft < 0) then
            status = EDGE_NO_FOLD
            return
        endif
        do j = 1,45
            t = 0.5_dp*(high+low)
            call edge_equilibrium_density(model,t,n,ft,rc)
            if (rc /= EDGE_OK) return
            if (ft >= 0) then
                low = t
            else
                high = t
            endif
        enddo
        root = 0.5_dp*(high+low)
        call edge_equilibrium_density(model,root,n_limit,ft,rc)
        if (rc /= EDGE_OK) return
        t_fold = root
        status = EDGE_OK
    end subroutine

    !> Empirical Greenwald reference only, input I_p in amperes and a in m.
    !> Result is line-averaged n_e in m^-3, unlike edge_density_limit.
    pure subroutine greenwald_reference_density(ip,minor_a,n_gw,status)
        real(dp), intent(in) :: ip,minor_a
        real(dp), intent(out) :: n_gw
        integer, intent(out) :: status
        n_gw=0; status=EDGE_INVALID
        if (.not. all(ieee_is_finite([ip,minor_a]))) return
        if (ip <= 0 .or. minor_a <= 0) return
        n_gw = (ip/1e6_dp)/(acos(-1.0_dp)*minor_a**2)*1e20_dp
        status=EDGE_OK
    end subroutine

    !> Browser-friendly fixed ABI: 14 inputs, 3 physical outputs, integer status.
    !> input = [j, Pext, Zeff, lnL, chi, width, Tw, fZ, Lpeak, Tpeak,
    !>          low_power, high_power, Tmin, Tmax].
    subroutine milca_edge_density_abi(input,n_limit,t_fold,status) bind(C)
        real(c_double), intent(in) :: input(14)
        real(c_double), intent(out) :: n_limit,t_fold
        integer(c_int), intent(out) :: status
        type(edge_model_t) :: model
        integer :: rc
        model%j_parallel=input(1); model%p_external=input(2)
        model%zeff=input(3); model%coulomb_log=input(4)
        model%chi_e=input(5); model%shell_width=input(6)
        model%t_wall=input(7); model%impurity_fraction=input(8)
        model%cooling_peak=input(9); model%cooling_tpeak=input(10)
        model%cooling_low_power=input(11)
        model%cooling_high_power=input(12)
        call edge_density_limit(model,input(13),input(14),n_limit,t_fold,rc)
        status=int(rc,c_int)
    end subroutine
end module
