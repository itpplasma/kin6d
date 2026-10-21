module kin6d_ideal_xi
    ! D5 first increment: plasma-volume Hessian in physical displacement jets.
    use, intrinsic :: iso_fortran_env, only: real32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_ideal_xi_leaves, only: xi_form_dp, xi_form_sp, &
        xi_gradient_dp, xi_gradient_sp
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
end module kin6d_ideal_xi
