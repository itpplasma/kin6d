module milca_vrmhd_boundary
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_vrmhd_section, only: vrmhd_section_t
    implicit none
    private
    public :: vrmhd_current_lift
contains
    subroutine vrmhd_current_lift(op, trace, wall)
        type(vrmhd_section_t), intent(in) :: op
        real(dp), intent(in) :: trace(op%fs%sp(1)%ndof)
        real(dp), intent(out) :: wall(op%fs%sp(1)%ndof)
        wall = 0
        where (op%fs%sp(1)%bnd) wall = trace
    end subroutine
end module milca_vrmhd_boundary
