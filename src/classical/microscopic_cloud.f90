!> C2: isolated electron cloud in a fixed neutralizing uniform sphere, with an
!> optional fixed point impurity at the origin. Electrons, sphere and impurity
!> are declared through kin6d_species (atomic units: q_e=-1, m_e=1,
!> 1/(4 pi eps0)=1); dynamics use the shared kin6d_nbody KDK path, any N.
module kin6d_microscopic_cloud
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_species, only: units_t, species_t, atomic_units, electron
    use kin6d_nbody, only: nbody_system_t, nbody_isolated, nbody_set_sphere, &
                           nbody_add_source
    use kin6d_cloud_leaves, only: cloud_initial_position, cloud_initial_velocity
    implicit none
    private
    public :: cloud_model_t, cloud_system, cloud_initialize

    type :: cloud_model_t
        !> Total background charge and radius of the sphere (units of e, a0).
        real(dp) :: sphere_charge = 32.0_dp
        real(dp) :: sphere_radius = 3.0_dp
        !> Fixed point impurity at the origin (units of e); 0 means none.
        real(dp) :: impurity_charge = 0.0_dp
    end type
contains
    subroutine cloud_system(model, particles, system, ok)
        type(cloud_model_t), intent(in) :: model
        integer, intent(in) :: particles
        type(nbody_system_t), intent(out) :: system
        logical, intent(out) :: ok
        type(species_t), allocatable :: electrons(:)
        type(units_t) :: units
        ok = .false.
        if (particles < 1 .or. model%sphere_charge < 0) return
        if (model%impurity_charge > 0) return
        units = atomic_units()
        allocate (electrons(particles), source=electron())
        call nbody_isolated(system, units, electrons, ok)
        if (.not. ok) return
        call nbody_set_sphere(system, model%sphere_charge, model%sphere_radius, ok)
        if (.not. ok) return
        if (abs(model%impurity_charge) > 0) call nbody_add_source(system, &
            [0.0_dp, 0.0_dp, 0.0_dp], model%impurity_charge, ok)
    end subroutine

    !> Deterministic Fibonacci-sphere start for n = size(state)/6 electrons.
    subroutine cloud_initialize(state)
        real(dp), intent(out) :: state(:)
        integer :: p, j, n
        n = size(state)/6
        do p = 1, n
            j = 3*p - 2
            call cloud_initial_position(real(p, dp), real(n, dp), state(j), &
                                        state(j + 1), state(j + 2))
            call cloud_initial_velocity(real(p, dp), state(3*n + j), &
                                        state(3*n + j + 1), state(3*n + j + 2))
        end do
    end subroutine
end module
