module milca_ideal_xi
    ! D5 first increment: plasma-volume Hessian in physical displacement jets.
    use, intrinsic :: iso_fortran_env, only: real32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use milca_ideal_xi_leaves, only: xi_form_dp, xi_form_sp, &
        xi_gradient_dp, xi_gradient_sp, xi_form_jvp_dp, xi_form_jvp_sp, &
        xi_form_vjp_dp, xi_form_vjp_sp
    implicit none
    private
    public :: xi_bilinear, xi_gradient_n, xi_volume_energy, xi_energy_report_t
    integer, parameter, public :: XI_OK = 0, XI_BAD_INPUT = 1
    type :: xi_energy_report_t
        integer :: status = XI_BAD_INPUT
        character(16) :: claim = 'candidate'
        character(16) :: bound_status = 'unavailable'
    end type xi_energy_report_t
    interface xi_bilinear
        module procedure bilinear_dp, bilinear_sp
    end interface
    interface xi_gradient_n
        module procedure gradient_dp, gradient_sp
    end interface
    interface xi_volume_energy
        module procedure volume_dp, volume_sp
    end interface
    public :: xi_bilinear_batch
    !> SoA columns: inputs B(3), p, mu, A(3,3), C(3,3), 23 tangents and one
    !> cotangent; outputs value, JVP and VJP(23).
    interface xi_bilinear_batch
        module procedure batch_dp, batch_sp, batch_mixed
    end interface

contains
#define RK real64
#define BILINEAR bilinear_dp
#define GRADIENT gradient_dp
#define VOLUME volume_dp
#define FORM xi_form_dp
#define GRAD xi_gradient_dp
#include "ideal_xi_impl.inc"
#undef RK
#undef BILINEAR
#undef GRADIENT
#undef VOLUME
#undef FORM
#undef GRAD
#define RK real32
#define BILINEAR bilinear_sp
#define GRADIENT gradient_sp
#define VOLUME volume_sp
#define FORM xi_form_sp
#define GRAD xi_gradient_sp
#include "ideal_xi_impl.inc"
#undef RK
#undef BILINEAR
#undef GRADIENT
#undef VOLUME
#undef FORM
#undef GRAD

    subroutine batch_dp(n, x, y)
        integer, intent(in) :: n
        real(real64), intent(in) :: x(n, 47)
        real(real64), intent(out) :: y(n, 25)
#define RK real64
#define FORM xi_form_dp
#define JVP xi_form_jvp_dp
#define VJP xi_form_vjp_dp
#include "ideal_xi_batch.inc"
#undef RK
#undef FORM
#undef JVP
#undef VJP
    end subroutine batch_dp

    subroutine batch_sp(n, x, y)
        integer, intent(in) :: n
        real(real32), intent(in) :: x(n, 47)
        real(real32), intent(out) :: y(n, 25)
#define RK real32
#define FORM xi_form_sp
#define JVP xi_form_jvp_sp
#define VJP xi_form_vjp_sp
#include "ideal_xi_batch.inc"
#undef RK
#undef FORM
#undef JVP
#undef VJP
    end subroutine batch_sp

    subroutine batch_mixed(n, x, y)
        ! FP32 storage, FP64 evaluation and output.
        integer, intent(in) :: n
        real(real32), intent(in) :: x(n, 47)
        real(real64), intent(out) :: y(n, 25)
        call batch_dp(n, real(x, real64), y)
    end subroutine batch_mixed
end module milca_ideal_xi
