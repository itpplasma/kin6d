module milca_gs_symmetry
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_gs_symmetry_leaves, only: gs_screw_profile_source, &
        gs_translation_source, gs_translation_source_jvp, &
        gs_screw_field, gs_translation_field
    use kin6d_grad_shafranov_leaves, only: gs_source, gs_source_jvp, gs_field
    implicit none
    private
    integer, parameter, public :: GS_ROTATION = 1, GS_TRANSLATION = 2, GS_SCREW = 3
    integer, parameter, public :: GS_CHART_RZ = 1, GS_CHART_CARTESIAN = 2
    integer, parameter, public :: GS_UNSUPPORTED_SYMMETRY = 23
    type, public :: gs_symmetry_t
        integer :: kind = GS_ROTATION
        integer :: chart = GS_CHART_RZ
        real(dp) :: k = 0
    contains
        procedure :: valid => symmetry_valid
    end type
    public :: gs_reduced_source, gs_reduced_field
contains
    logical function symmetry_valid(symmetry) result(valid)
        class(gs_symmetry_t), intent(in) :: symmetry
        valid = .false.
        if (.not. ieee_is_finite(symmetry%k)) return
        select case (symmetry%kind)
        case (GS_ROTATION)
            if (symmetry%chart /= GS_CHART_RZ) return
            if (abs(symmetry%k) > 0) return
        case (GS_TRANSLATION)
            if (symmetry%chart /= GS_CHART_CARTESIAN) return
            if (abs(symmetry%k) > 0) return
        case (GS_SCREW)
            if (symmetry%chart /= GS_CHART_CARTESIAN) return
        case default
            return
        end select
        valid = .true.
    end function

    pure subroutine gs_reduced_field(symmetry, x, y, px, py, f, b1, b2, b3)
        type(gs_symmetry_t), intent(in) :: symmetry
        real(dp), intent(in) :: x, y, px, py, f
        real(dp), intent(out) :: b1, b2, b3
        select case (symmetry%kind)
        case (GS_ROTATION)
            call gs_field(x, px, py, f, b1, b2, b3)
        case (GS_TRANSLATION)
            call gs_translation_field(px, py, f, b1, b2, b3)
        case (GS_SCREW)
            call gs_screw_field(x, y, symmetry%k, px, py, f, b1, b2, b3)
        end select
    end subroutine

    subroutine gs_reduced_source(symmetry, x, y, mu0, pp, ff, f, dpp, dff, source, tangent)
        type(gs_symmetry_t), intent(in) :: symmetry
        real(dp), intent(in) :: x, y, mu0, pp, ff, f, dpp, dff
        real(dp), intent(out) :: source, tangent
        select case (symmetry%kind)
        case (GS_ROTATION)
            call gs_source(x, mu0, pp, ff, source)
            call gs_source_jvp(x, mu0, pp, ff, 0.0_dp, 0.0_dp, dpp, dff, tangent)
        case (GS_TRANSLATION)
            call gs_translation_source(mu0, pp, ff, source)
            call gs_translation_source_jvp(mu0, pp, ff, 0.0_dp, dpp, dff, tangent)
        case (GS_SCREW)
            call gs_screw_profile_source(x, y, symmetry%k, mu0, pp, ff, f, &
                dpp, dff, source, tangent)
        end select
    end subroutine
end module milca_gs_symmetry
