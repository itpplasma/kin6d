module kin6d_gs_readback
    ! psi-lane entry to the core readback (docs/equilibrium-core.md, D7): the
    ! state is wrapped by the psi adapter and read by read_equilibrium.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_equilibrium_state, only: equilibrium_state_t
    use kin6d_equilibrium_readback, only: gs_readback_t => equilibrium_readback_t, &
        read_equilibrium, gs_readback_q => readback_q
    use kin6d_gs_equilibrium, only: gs_equilibrium, gs_profile_f, &
        gs_readback_resolution
    implicit none
    private
    public :: gs_readback_t, read_gs_equilibrium, gs_profile_f, gs_readback_q
    public :: gs_readback_resolution
contains
    subroutine read_gs_equilibrium(state, f_edge, ffprime, radial_points, &
                                   angular_points, result, status, law, &
                                   psi_bound, plasma_cell, rho_grid)
        !! Flux-surface readback: q, normalised toroidal flux, volume and edge
        !! toroidal flux on rho = sqrt((psi - psi_a)/(psi_b - psi_a)) levels.
        !!
        !! Fixed boundary (psi_bound absent): psi = 0 on the mesh boundary, the
        !! whole mesh is plasma; volume and phi_edge are whole-mesh sums and
        !! flux_consistency compares phi_edge with the contour integral.
        !!
        !! Free boundary (psi_bound present, psi = 0 on the symmetry axis): the
        !! plasma boundary is the level psi = psi_bound and the plasma lies in
        !! the cells flagged in plasma_cell (default all), which must contain
        !! every flux surface inside it. The law is in the flux measured from
        !! psi_bound (law%axis = psi_a - psi_b). On a plasma sub-mesh phi_edge
        !! is the contour integral itself, volume and flux_consistency are NaN
        !! (the cut cells need the level-cut quadrature of the solver).
        !! result%psi and psi_axis are absolute. rho_grid (optional) replaces
        !! the uniform rho levels (strictly increasing from 0 to 1).
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: f_edge, ffprime
        integer, intent(in) :: radial_points, angular_points
        type(gs_readback_t), intent(out) :: result
        integer, intent(out) :: status
        type(gs_profile_law_t), intent(in), optional :: law
        real(dp), intent(in), optional :: psi_bound
        logical, intent(in), optional :: plasma_cell(:)
        real(dp), intent(in), optional :: rho_grid(:)
        type(equilibrium_state_t) :: eq
        status = 1
        if (present(plasma_cell) .and. .not. present(psi_bound)) return
        if (.not. all(ieee_is_finite([f_edge, ffprime]))) return
        if (abs(f_edge) <= tiny(1.0_dp)) return
        call gs_equilibrium(state, eq, status, f_edge, ffprime, law=law, &
            psi_bound=psi_bound, plasma_cell=plasma_cell)
        if (status /= 0) then
            status = 1
            return
        end if
        call read_equilibrium(eq, radial_points, angular_points, result, status, &
                              rho_grid)
    end subroutine
end module kin6d_gs_readback
