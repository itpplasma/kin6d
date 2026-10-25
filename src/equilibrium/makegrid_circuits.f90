module kin6d_makegrid_circuits
    !! Circuits of a MAKEGRID coils file, read with the libneo reader
    !! (neo_makegrid_coils, https://github.com/itpplasma/libneo). Coil c becomes
    !! a closed coils_3d polygon in circuit(c) = its MAKEGRID group (the EXTCUR
    !! index) with turns(c) = its fourth-column value, which must be constant
    !! along the coil (the winding count for W7-X style files). The circuit
    !! field b_k = sum_{circuit(c) = k} turns(c) B_c[1 A] is then the field per
    !! unit EXTCUR. Coils without a group (group 0) get circuit 0, which
    !! coil_circuit_basis_t treats as not connected.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_makegrid_coils, only: makegrid_coils_t, read_makegrid_coils
    use kin6d_coils_3d, only: coil_set_t, coil_set_from_polylines, coil_name_len
    implicit none
    private
    public :: makegrid_circuits_read

contains

    subroutine makegrid_circuits_read(path, coils, circuit, turns, ncircuit, status, &
                                      group_name, periods)
        !! status: 0 ok, 1 unreadable or malformed file, 2 a coil whose
        !! fourth column varies, 3 a coil with fewer than three points.
        character(len=*), intent(in) :: path
        type(coil_set_t), intent(out) :: coils
        integer, allocatable, intent(out) :: circuit(:)
        real(dp), allocatable, intent(out) :: turns(:)
        integer, intent(out) :: ncircuit, status
        character(len=*), allocatable, intent(out), optional :: group_name(:)
        integer, intent(out), optional :: periods
        type(makegrid_coils_t) :: mc
        integer, allocatable :: npts(:)
        real(dp), allocatable :: pts(:, :)
        character(len=coil_name_len), allocatable :: names(:)
        integer :: c, i0, i1, n, ierr

        ncircuit = 0
        call read_makegrid_coils(path, mc, ierr)
        status = 1
        if (ierr /= 0) return
        status = 3
        if (mc%n_coils < 1) return
        allocate (npts(mc%n_coils), circuit(mc%n_coils), turns(mc%n_coils), &
            names(mc%n_coils))
        npts = mc%last - mc%first + 1
        if (any(npts < 3)) return
        allocate (pts(3, sum(npts)))
        n = 0
        do c = 1, mc%n_coils
            i0 = mc%first(c); i1 = mc%last(c)
            status = 2
            if (maxval(abs(mc%current(i0:i1) - mc%current(i0))) > 0) return
            pts(:, n + 1:n + npts(c)) = transpose(reshape([mc%x(i0:i1), mc%y(i0:i1), &
                mc%z(i0:i1)], [npts(c), 3]))
            n = n + npts(c)
            circuit(c) = max(0, mc%group(c))
            turns(c) = mc%current(i0)
            names(c) = ''
            if (mc%group(c) >= 1) names(c) = mc%group_name(mc%group(c))(1:coil_name_len)
        end do
        call coil_set_from_polylines(coils, npts, pts, names, status)
        if (status /= 0) then
            status = 3
            return
        end if
        ncircuit = mc%n_groups
        if (present(group_name)) then
            allocate (group_name(mc%n_groups))
            group_name = mc%group_name
        end if
        if (present(periods)) periods = mc%periods
        status = 0
    end subroutine

end module kin6d_makegrid_circuits
