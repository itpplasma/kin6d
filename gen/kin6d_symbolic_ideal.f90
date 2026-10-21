module kin6d_symbolic_ideal
    ! Stable KS3/map import path; physics lives in the pinned consumer package.
    use ideal_mhd_symbolic, only: frozen_energy
    implicit none
    private
    public :: frozen_energy
end module kin6d_symbolic_ideal
