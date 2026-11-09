module kin6d_vrmhd_drive
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_gs_symmetry, only: GS_ROTATION
    use kin6d_vrmhd_section, only: vrmhd_section_t, vrmhd_section_basis
    use kin6d_vrmhd_section_leaves, only: vrmhd_harmonic_rotation, &
        vrmhd_harmonic_screw, vrmhd_pairing
    implicit none
    private
    integer, parameter, public :: VRMHD_CURRENT_SOURCE = 1, VRMHD_LOOP_VOLTAGE = 2
    public :: vrmhd_drive_build
contains
    subroutine vrmhd_drive_build(op, convention, axial_period, voltage, js, &
            harmonic, source, load, status)
        type(vrmhd_section_t), intent(in) :: op
        integer, intent(in) :: convention
        real(dp), intent(in) :: axial_period, voltage
        real(dp), intent(in) :: js(3, op%fs%nq, op%fs%mm%topology%n_triangles)
        real(dp), intent(out) :: harmonic(op%fs%sp(1)%ndof)
        real(dp), intent(inout) :: source(6, op%fs%nq, op%fs%mm%topology%n_triangles)
        real(dp), intent(out) :: load(op%ndof)
        integer, intent(out) :: status
        real(dp) :: basis(36, op%nlocal), x(2), weight, h(3), period, value, a
        integer :: t, q, i
        status = 1
        harmonic = 0; load = 0
        if (axial_period <= 0) return
        period = axial_period
        if (op%symmetry%kind == GS_ROTATION) period = 2*acos(-1.0_dp)
        select case (convention)
        case (VRMHD_CURRENT_SOURCE)
            if (abs(voltage) > 0) return
            source(:3, :, :) = js
        case (VRMHD_LOOP_VOLTAGE)
            source(:3, :, :) = 0
        case default
            return
        end select
        ! Section edges have zero coordinate increment; each periodic orbit has one.
        harmonic(op%fs%nn + 1:) = 1
        if (convention == VRMHD_CURRENT_SOURCE) then
            status = 0
            return
        end if
        do t = 1, op%fs%mm%topology%n_triangles
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                if (status /= 0) return
                a = dot_product(op%fs%tl(:, q), &
                    harmonic(op%fs%nn + op%fs%sp(0)%dof(:, t)))/period
                if (op%symmetry%kind == GS_ROTATION) then
                    call vrmhd_harmonic_rotation(x(1), a, h(1), h(2), h(3))
                else
                    call vrmhd_harmonic_screw(x(1), a, h(1), h(2), h(3))
                end if
                do i = 1, op%nlocal
                    if (op%row(i) /= 3) cycle
                    call vrmhd_pairing(h(1), h(2), h(3), &
                        basis(5, i), basis(6, i), basis(7, i), value)
                    load(op%map(i, t)) = load(op%map(i, t)) + voltage*weight*value
                end do
            end do
        end do
        status = 0
    end subroutine
end module kin6d_vrmhd_drive
