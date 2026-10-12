module kin6d_gs_coils
    !! Coils and circuits for the free-boundary Grad-Shafranov problem.
    !!
    !! Flux convention.  Every coil flux is the poloidal flux per radian
    !! psi = R A_phi (Wb/rad) at a field point (R, Z) per ampere of the
    !! circuit current I_k flowing in the conductor,
    !!
    !!     psi_c(R, Z) = mu0 * turns_c * avg_{cross-section} G(R, Z; R', Z'),
    !!
    !! with the Green's function G of fortfem_axisymmetric_gs_kernel (no mu0).
    !! The flux vector of circuit k is g_k = sum_c A(c,k) psi_c, where A is
    !! the coil-by-circuit incidence matrix: coil c carries the current
    !! sum_k A(c,k) I_k.  A series pair (CS1U, CS1L) is one column with two
    !! entries +1; a reversed winding has the entry -1.  The total flux is
    !! psi = sum_k I_k g_k.
    !!
    !! Coil representations (gs_coil_t%kind)
    !!   gs_coil_rect            rectangle, tensor Gauss-Legendre average of
    !!                           quad_order x quad_order points (fortfem rule)
    !!   gs_coil_rect_filaments  rectangle replaced by a grid of equal-current
    !!                           filaments at the cell centres, nr = max(1,
    !!                           ceiling(dr/fil_size)), nz likewise.  This is
    !!                           the FreeGS MultiCoil discretisation of the
    !!                           ITER parity runs (fil_size = 0.15 m); the
    !!                           flux is mu0 * turns / (nr nz) * sum G.
    !!   gs_coil_filaments       explicit filaments (fil_r, fil_z, fil_turns);
    !!                           the flux is mu0 * sum fil_turns * G.
    !! Filament positions must not coincide with a field point.
    !!
    !! Current limits.  The limits i_min, i_max of a circuit are in A (of the
    !! conductor current I_k).  The MA-turn figure of a circuit is
    !! I_k * turns_k / 1e6 with turns_k the turns of its first member coil.
    !!
    !! Text coil table (neutral format; units m, A, turns).  Lines; '#'
    !! starts a comment; tokens are blank separated; names contain no blank.
    !!
    !!     kin6d-coil-table 1
    !!     coil NAME rect R Z DR DZ TURNS [quad=N] [fil=SIZE]
    !!     coil NAME filament R Z TURNS
    !!     circuit NAME IMIN IMAX MEMBER [MEMBER ...]
    !!
    !! 'rect' is a rectangle with centre (R, Z) and full sizes (DR, DZ);
    !! 'quad=N' sets the Gauss order (default 12), 'fil=SIZE' selects the
    !! FreeGS filament grid.  'filament' lines with the same NAME accumulate
    !! the filaments of one coil.  A MEMBER is a coil name or NAME*FACTOR
    !! (incidence entry FACTOR, default 1).  A coil that belongs to no
    !! circuit is uncontrolled (zero incidence row).  Circuits refer to coils
    !! declared earlier in the file.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_axisymmetric_gs_kernel, only: gs_mu0, gs_rect_coil_t, &
        gs_rect_coil_flux, gs_rect_coil_flux_gradient, gs_green, gs_green_gradient
    implicit none
    private

    integer, parameter, public :: gs_coil_rect = 1
    integer, parameter, public :: gs_coil_rect_filaments = 2
    integer, parameter, public :: gs_coil_filaments = 3
    integer, parameter, public :: gs_name_len = 32

    type, public :: gs_coil_t
        character(len=gs_name_len) :: name = ''
        integer :: kind = gs_coil_rect
        type(gs_rect_coil_t) :: rect
        integer :: quad_order = 12
        real(dp) :: fil_size = 0.0_dp
        real(dp), allocatable :: fil_r(:), fil_z(:), fil_turns(:)
    end type gs_coil_t

    type, public :: gs_circuit_set_t
        type(gs_coil_t), allocatable :: coils(:)
        character(len=gs_name_len), allocatable :: circuit_names(:)
        real(dp), allocatable :: incidence(:, :) !< (coil, circuit)
        real(dp), allocatable :: i_min(:), i_max(:) !< A
    end type gs_circuit_set_t

    public :: gs_coil_turns, gs_coil_make_freegs_filaments, gs_coil_filament_list
    public :: gs_coil_flux, gs_coil_flux_gradient
    public :: gs_circuit_turns, gs_circuit_flux_vectors
    public :: gs_coil_currents
    public :: gs_amps_to_mat, gs_mat_to_amps
    public :: gs_circuit_limits_mat, gs_set_circuit_limits_mat
    public :: gs_circuit_set_check
    public :: gs_read_coil_table, gs_write_coil_table

contains

    pure function gs_amps_to_mat(amps, turns) result(mat)
        !! Conductor current (A) to MA-turns.
        real(dp), intent(in) :: amps, turns
        real(dp) :: mat
        mat = amps*turns*1.0e-6_dp
    end function gs_amps_to_mat

    pure function gs_mat_to_amps(mat, turns) result(amps)
        !! MA-turns to conductor current (A); turns must be nonzero.
        real(dp), intent(in) :: mat, turns
        real(dp) :: amps
        amps = mat*1.0e6_dp/turns
    end function gs_mat_to_amps

    pure function gs_coil_turns(coil) result(turns)
        !! Total number of turns of a coil.
        type(gs_coil_t), intent(in) :: coil
        real(dp) :: turns
        if (coil%kind == gs_coil_filaments) then
            turns = sum(coil%fil_turns)
        else
            turns = coil%rect%turns
        end if
    end function gs_coil_turns

    subroutine gs_coil_make_freegs_filaments(coil, fil_size)
        !! Switch a rectangle to the FreeGS MultiCoil filament grid.
        type(gs_coil_t), intent(inout) :: coil
        real(dp), intent(in) :: fil_size
        coil%kind = gs_coil_rect_filaments
        coil%fil_size = fil_size
    end subroutine gs_coil_make_freegs_filaments

    subroutine gs_coil_filament_list(coil, r, z, w)
        !! Positions and turns weights of the filaments of a coil of kind
        !! gs_coil_rect_filaments or gs_coil_filaments.
        type(gs_coil_t), intent(in) :: coil
        real(dp), allocatable, intent(out) :: r(:), z(:), w(:)
        integer :: nr, nz, i, j, k

        if (coil%kind == gs_coil_filaments) then
            r = coil%fil_r
            z = coil%fil_z
            w = coil%fil_turns
            return
        end if
        nr = max(1, ceiling(abs(coil%rect%dr)/coil%fil_size))
        nz = max(1, ceiling(abs(coil%rect%dz)/coil%fil_size))
        allocate (r(nr*nz), z(nr*nz), w(nr*nz))
        k = 0
        do j = 1, nz
            do i = 1, nr
                k = k + 1
                r(k) = coil%rect%r + coil%rect%dr*((i - 0.5_dp)/nr - 0.5_dp)
                z(k) = coil%rect%z + coil%rect%dz*((j - 0.5_dp)/nz - 0.5_dp)
                w(k) = coil%rect%turns/real(nr*nz, dp)
            end do
        end do
    end subroutine gs_coil_filament_list

    function gs_coil_flux(coil, r, z) result(psi)
        !! Flux per radian (Wb/rad) at (r, z) per ampere of coil current.
        type(gs_coil_t), intent(in) :: coil
        real(dp), intent(in) :: r, z
        real(dp) :: psi
        real(dp), allocatable :: fr(:), fz(:), fw(:)
        integer :: k

        if (coil%kind == gs_coil_rect) then
            psi = gs_rect_coil_flux(coil%rect, r, z, coil%quad_order, coil%quad_order)
        else
            call gs_coil_filament_list(coil, fr, fz, fw)
            psi = 0.0_dp
            do k = 1, size(fr)
                psi = psi + fw(k)*gs_green(r, z, fr(k), fz(k))
            end do
            psi = gs_mu0*psi
        end if
    end function gs_coil_flux

    subroutine gs_coil_flux_gradient(coil, r, z, psi, dpsi_dr, dpsi_dz)
        !! Coil flux and its derivatives with respect to the field point.
        !! B_R = -dpsi_dz / R and B_Z = dpsi_dr / R per ampere.
        type(gs_coil_t), intent(in) :: coil
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: psi, dpsi_dr, dpsi_dz
        real(dp), allocatable :: fr(:), fz(:), fw(:)
        real(dp) :: g, gg(4), grad6(6)
        integer :: k

        if (coil%kind == gs_coil_rect) then
            call gs_rect_coil_flux_gradient(coil%rect, r, z, coil%quad_order, &
                coil%quad_order, psi, grad6)
            dpsi_dr = grad6(1)
            dpsi_dz = grad6(2)
        else
            call gs_coil_filament_list(coil, fr, fz, fw)
            psi = 0.0_dp
            dpsi_dr = 0.0_dp
            dpsi_dz = 0.0_dp
            do k = 1, size(fr)
                call gs_green_gradient(r, z, fr(k), fz(k), g, gg)
                psi = psi + fw(k)*g
                dpsi_dr = dpsi_dr + fw(k)*gg(1)
                dpsi_dz = dpsi_dz + fw(k)*gg(2)
            end do
            psi = gs_mu0*psi
            dpsi_dr = gs_mu0*dpsi_dr
            dpsi_dz = gs_mu0*dpsi_dz
        end if
    end subroutine gs_coil_flux_gradient

    subroutine gs_circuit_flux_vectors(set, r, z, g, dg_dr, dg_dz)
        !! Circuit flux vectors g(i, k) = sum_c A(c,k) psi_c(r(i), z(i)) on a
        !! set of points (boundary-curve nodes or quadrature points) and,
        !! when both dg_dr and dg_dz are given, their point gradients.
        type(gs_circuit_set_t), intent(in) :: set
        real(dp), intent(in) :: r(:), z(:)
        real(dp), intent(out) :: g(:, :)
        real(dp), intent(out), optional :: dg_dr(:, :), dg_dz(:, :)
        real(dp) :: psi, pr, pz
        integer :: c, i, k
        logical :: gradient

        gradient = present(dg_dr) .and. present(dg_dz)
        g = 0.0_dp
        if (present(dg_dr)) dg_dr = 0.0_dp
        if (present(dg_dz)) dg_dz = 0.0_dp
        ! Points are independent; per point the coils are summed in their
        ! order, so the result does not depend on the thread count.
        !$omp parallel do default(shared) private(i, c, k, psi, pr, pz) &
        !$omp schedule(dynamic, 16) if (size(r) > 8)
        do i = 1, size(r)
            do c = 1, size(set%coils)
                if (all(set%incidence(c, :) == 0.0_dp)) cycle
                if (gradient) then
                    call gs_coil_flux_gradient(set%coils(c), r(i), z(i), psi, pr, pz)
                else
                    psi = gs_coil_flux(set%coils(c), r(i), z(i))
                    pr = 0.0_dp
                    pz = 0.0_dp
                end if
                do k = 1, size(set%incidence, 2)
                    g(i, k) = g(i, k) + set%incidence(c, k)*psi
                    if (gradient) then
                        dg_dr(i, k) = dg_dr(i, k) + set%incidence(c, k)*pr
                        dg_dz(i, k) = dg_dz(i, k) + set%incidence(c, k)*pz
                    end if
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine gs_circuit_flux_vectors

    pure function gs_coil_currents(set, i_circuit) result(i_coil)
        !! Coil currents (A) from the circuit currents (A).
        type(gs_circuit_set_t), intent(in) :: set
        real(dp), intent(in) :: i_circuit(:)
        real(dp) :: i_coil(size(set%coils))
        i_coil = matmul(set%incidence, i_circuit)
    end function gs_coil_currents

    pure function gs_circuit_turns(set, k) result(turns)
        !! Turns of circuit k: turns of its first member coil (0 if none).
        type(gs_circuit_set_t), intent(in) :: set
        integer, intent(in) :: k
        real(dp) :: turns
        integer :: c

        turns = 0.0_dp
        do c = 1, size(set%coils)
            if (set%incidence(c, k) /= 0.0_dp) then
                turns = gs_coil_turns(set%coils(c))
                return
            end if
        end do
    end function gs_circuit_turns

    subroutine gs_circuit_limits_mat(set, lo, hi)
        !! Circuit current limits in MA-turns, lo <= hi even for negative turns.
        type(gs_circuit_set_t), intent(in) :: set
        real(dp), intent(out) :: lo(:), hi(:)
        real(dp) :: a, b
        integer :: k

        do k = 1, size(set%i_min)
            a = gs_amps_to_mat(set%i_min(k), gs_circuit_turns(set, k))
            b = gs_amps_to_mat(set%i_max(k), gs_circuit_turns(set, k))
            lo(k) = min(a, b)
            hi(k) = max(a, b)
        end do
    end subroutine gs_circuit_limits_mat

    subroutine gs_set_circuit_limits_mat(set, lo, hi)
        !! Set the circuit limits from MA-turn bounds.
        type(gs_circuit_set_t), intent(inout) :: set
        real(dp), intent(in) :: lo(:), hi(:)
        real(dp) :: a, b
        integer :: k

        do k = 1, size(lo)
            a = gs_mat_to_amps(lo(k), gs_circuit_turns(set, k))
            b = gs_mat_to_amps(hi(k), gs_circuit_turns(set, k))
            set%i_min(k) = min(a, b)
            set%i_max(k) = max(a, b)
        end do
    end subroutine gs_set_circuit_limits_mat

    subroutine gs_circuit_set_check(set, stat, message)
        !! Consistency of array sizes, filament data and limits.
        type(gs_circuit_set_t), intent(in) :: set
        integer, intent(out) :: stat
        character(len=*), intent(out) :: message
        integer :: c, nc, nk

        stat = 0
        message = ''
        if (.not. allocated(set%coils) .or. .not. allocated(set%incidence) .or. &
            .not. allocated(set%i_min) .or. .not. allocated(set%i_max)) then
            stat = 1
            message = 'circuit set is not fully allocated'
            return
        end if
        nc = size(set%coils)
        nk = size(set%incidence, 2)
        if (size(set%incidence, 1) /= nc .or. size(set%i_min) /= nk .or. &
            size(set%i_max) /= nk) then
            stat = 2
            message = 'incidence matrix and limit sizes disagree'
            return
        end if
        do c = 1, nc
            if (set%coils(c)%kind == gs_coil_filaments) then
                if (.not. allocated(set%coils(c)%fil_r)) stat = 3
            else if (set%coils(c)%quad_order < 1) then
                stat = 3
            else if (set%coils(c)%kind == gs_coil_rect_filaments .and. &
                    .not. set%coils(c)%fil_size > 0.0_dp) then
                stat = 3
            end if
            if (stat /= 0) then
                message = 'invalid coil data: '//trim(set%coils(c)%name)
                return
            end if
        end do
        if (any(set%i_min > set%i_max)) then
            stat = 4
            message = 'a circuit has i_min > i_max'
        end if
    end subroutine gs_circuit_set_check

    ! ------------------------------------------------------------------
    ! Coil table
    ! ------------------------------------------------------------------

    subroutine gs_write_coil_table(path, set, stat)
        !! Write the coil table; reals use 17 significant digits (exact).
        character(len=*), intent(in) :: path
        type(gs_circuit_set_t), intent(in) :: set
        integer, intent(out) :: stat
        character(len=512) :: line
        character(len=64) :: member
        real(dp), allocatable :: fr(:), fz(:), fw(:)
        integer :: u, c, k, f

        open (newunit=u, file=path, status='replace', action='write', iostat=stat)
        if (stat /= 0) return
        write (u, '(a)') 'kin6d-coil-table 1'
        write (u, '(a)') '# units: m, A, turns; flux per ampere includes turns'
        do c = 1, size(set%coils)
            associate (coil => set%coils(c))
                if (coil%kind == gs_coil_filaments) then
                    call gs_coil_filament_list(coil, fr, fz, fw)
                    do f = 1, size(fr)
                        write (u, '(4a,3(1x,es24.16e3))') 'coil ', &
                            trim(coil%name), ' filament', '', fr(f), fz(f), fw(f)
                    end do
                else
                    write (line, '(3a,5(1x,es24.16e3),a,i0)') 'coil ', trim(coil%name), &
                        ' rect', coil%rect%r, coil%rect%z, coil%rect%dr, coil%rect%dz, &
                        coil%rect%turns, ' quad=', coil%quad_order
                    if (coil%kind == gs_coil_rect_filaments) then
                        line = trim(line)//' fil='//real_text(coil%fil_size)
                    end if
                    write (u, '(a)') trim(line)
                end if
            end associate
        end do
        do k = 1, size(set%circuit_names)
            write (line, '(3a,2(1x,es24.16e3))') 'circuit ', &
                trim(set%circuit_names(k)), '', set%i_min(k), set%i_max(k)
            do c = 1, size(set%coils)
                if (set%incidence(c, k) == 0.0_dp) cycle
                if (set%incidence(c, k) == 1.0_dp) then
                    member = trim(set%coils(c)%name)
                else
                    member = trim(set%coils(c)%name)//'*'// &
                        real_text(set%incidence(c, k))
                end if
                line = trim(line)//' '//trim(member)
            end do
            write (u, '(a)') trim(line)
        end do
        close (u)
    end subroutine gs_write_coil_table

    subroutine gs_read_coil_table(path, set, stat, message)
        !! Read a coil table (format in the module header).  stat = 0 on
        !! success; otherwise message names the offending line.
        character(len=*), intent(in) :: path
        type(gs_circuit_set_t), intent(out) :: set
        integer, intent(out) :: stat
        character(len=*), intent(out) :: message
        type(gs_coil_t), allocatable :: coils(:)
        character(len=gs_name_len), allocatable :: cnames(:)
        real(dp), allocatable :: inc(:, :), imin(:), imax(:)
        character(len=512) :: line
        character(len=128), allocatable :: tok(:)
        real(dp) :: x(5), fac
        integer :: u, ios, nline, ntok, nc, nk, c, k, t, star, pos
        logical :: ok, header

        stat = 0
        message = ''
        nc = 0
        nk = 0
        nline = 0
        header = .false.
        allocate (coils(0), cnames(0), imin(0), imax(0), inc(0, 0))
        open (newunit=u, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            stat = 10
            message = 'cannot open coil table '//trim(path)
            return
        end if
        do
            read (u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            nline = nline + 1
            pos = index(line, '#')
            if (pos > 0) line(pos:) = ''
            call split_tokens(line, tok, ntok)
            if (ntok == 0) cycle
            if (.not. header) then
                if (ntok /= 2 .or. tok(1) /= 'kin6d-coil-table' .or. tok(2) /= '1') then
                    call fail(11, 'missing header "kin6d-coil-table 1"')
                    exit
                end if
                header = .true.
                cycle
            end if
            select case (trim(tok(1)))
            case ('coil')
                if (ntok < 4) then
                    call fail(12, 'short coil line'); exit
                end if
                c = find_name(coils, nc, tok(2))
                select case (trim(tok(3)))
                case ('rect')
                    if (c /= 0 .or. ntok < 8) then
                        call fail(13, 'bad or duplicate rect coil'); exit
                    end if
                    call parse_reals(tok(4:8), x, 5, ok)
                    if (.not. ok) then
                        call fail(14, 'bad number in rect coil'); exit
                    end if
                    call append_coil(coils, nc, tok(2))
                    c = nc
                    coils(c)%rect = gs_rect_coil_t(x(1), x(2), x(3), x(4), x(5))
                    do t = 9, ntok
                        if (tok(t)(1:5) == 'quad=') then
                            read (tok(t)(6:), *, iostat=ios) coils(c)%quad_order
                            ok = ios == 0
                        else if (tok(t)(1:4) == 'fil=') then
                            read (tok(t)(5:), *, iostat=ios) coils(c)%fil_size
                            coils(c)%kind = gs_coil_rect_filaments
                            ok = ios == 0
                        else
                            ok = .false.
                        end if
                        if (.not. ok) exit
                    end do
                    if (.not. ok) then
                        call fail(15, 'bad rect option'); exit
                    end if
                case ('filament')
                    if (ntok /= 6) then
                        call fail(16, 'filament line needs R Z TURNS'); exit
                    end if
                    call parse_reals(tok(4:6), x, 3, ok)
                    if (.not. ok) then
                        call fail(14, 'bad number in filament'); exit
                    end if
                    if (c == 0) then
                        call append_coil(coils, nc, tok(2))
                        c = nc
                        coils(c)%kind = gs_coil_filaments
                        allocate (coils(c)%fil_r(0), coils(c)%fil_z(0), &
                            coils(c)%fil_turns(0))
                    else if (coils(c)%kind /= gs_coil_filaments) then
                        call fail(13, 'filament added to a rect coil'); exit
                    end if
                    coils(c)%fil_r = [coils(c)%fil_r, x(1)]
                    coils(c)%fil_z = [coils(c)%fil_z, x(2)]
                    coils(c)%fil_turns = [coils(c)%fil_turns, x(3)]
                case default
                    call fail(12, 'unknown coil type'); exit
                end select
            case ('circuit')
                if (ntok < 5) then
                    call fail(17, 'circuit needs name, limits and members'); exit
                end if
                call parse_reals(tok(3:4), x, 2, ok)
                if (.not. ok .or. any(cnames(1:nk) == tok(2))) then
                    call fail(18, 'bad limits or duplicate circuit'); exit
                end if
                nk = nk + 1
                cnames = [cnames, tok(2)(1:gs_name_len)]
                imin = [imin, x(1)]
                imax = [imax, x(2)]
                call grow_incidence(inc, nc, nk)
                do t = 5, ntok
                    fac = 1.0_dp
                    star = index(tok(t), '*')
                    if (star > 0) then
                        read (tok(t)(star + 1:), *, iostat=ios) fac
                        if (ios /= 0) then
                            call fail(19, 'bad member factor'); exit
                        end if
                        tok(t)(star:) = ''
                    end if
                    c = find_name(coils, nc, tok(t))
                    if (c == 0) then
                        call fail(20, 'unknown coil '//trim(tok(t))); exit
                    end if
                    inc(c, nk) = fac
                end do
                if (stat /= 0) exit
            case default
                call fail(12, 'unknown record '//trim(tok(1))); exit
            end select
        end do
        close (u)
        if (stat == 0 .and. .not. header) then
            stat = 11
            message = 'empty coil table'
        end if
        if (stat /= 0) return
        ! Coils declared after the last circuit add zero rows.
        call grow_incidence(inc, nc, nk)
        allocate (set%coils(nc), set%circuit_names(nk), set%incidence(nc, nk), &
            set%i_min(nk), set%i_max(nk))
        set%coils = coils(1:nc)
        set%circuit_names = cnames(1:nk)
        set%incidence = inc(1:nc, 1:nk)
        set%i_min = imin(1:nk)
        set%i_max = imax(1:nk)
        call gs_circuit_set_check(set, stat, message)
        if (stat /= 0) stat = stat + 30

    contains

        subroutine fail(code, text)
            integer, intent(in) :: code
            character(len=*), intent(in) :: text
            character(len=16) :: num
            stat = code
            write (num, '(i0)') nline
            message = 'line '//trim(num)//': '//text
        end subroutine fail

    end subroutine gs_read_coil_table

    function real_text(x) result(text)
        !! Exactly round-tripping decimal text without blanks.
        real(dp), intent(in) :: x
        character(len=24) :: text
        write (text, '(es24.16e3)') x
        text = adjustl(text)
    end function real_text

    subroutine split_tokens(line, tok, ntok)
        character(len=*), intent(in) :: line
        character(len=128), allocatable, intent(out) :: tok(:)
        integer, intent(out) :: ntok
        integer :: i, start, n

        n = len_trim(line)
        allocate (tok(n/2 + 1))
        ntok = 0
        i = 1
        do while (i <= n)
            if (line(i:i) == ' ' .or. line(i:i) == char(9)) then
                i = i + 1
                cycle
            end if
            start = i
            do while (i <= n)
                if (line(i:i) == ' ' .or. line(i:i) == char(9)) exit
                i = i + 1
            end do
            ntok = ntok + 1
            tok(ntok) = line(start:i - 1)
        end do
    end subroutine split_tokens

    subroutine parse_reals(tok, x, n, ok)
        character(len=*), intent(in) :: tok(:)
        real(dp), intent(out) :: x(:)
        integer, intent(in) :: n
        logical, intent(out) :: ok
        integer :: i, ios

        ok = .true.
        x = 0.0_dp
        do i = 1, n
            read (tok(i), *, iostat=ios) x(i)
            if (ios /= 0) ok = .false.
        end do
    end subroutine parse_reals

    pure function find_name(coils, nc, name) result(c)
        type(gs_coil_t), intent(in) :: coils(:)
        integer, intent(in) :: nc
        character(len=*), intent(in) :: name
        integer :: c, i
        c = 0
        do i = 1, nc
            if (trim(coils(i)%name) == trim(name)) then
                c = i
                return
            end if
        end do
    end function find_name

    subroutine append_coil(coils, nc, name)
        type(gs_coil_t), allocatable, intent(inout) :: coils(:)
        integer, intent(inout) :: nc
        character(len=*), intent(in) :: name
        type(gs_coil_t), allocatable :: tmp(:)
        integer :: i

        allocate (tmp(nc + 1))
        do i = 1, nc
            call move_alloc(coils(i)%fil_r, tmp(i)%fil_r)
            call move_alloc(coils(i)%fil_z, tmp(i)%fil_z)
            call move_alloc(coils(i)%fil_turns, tmp(i)%fil_turns)
            tmp(i)%name = coils(i)%name
            tmp(i)%kind = coils(i)%kind
            tmp(i)%rect = coils(i)%rect
            tmp(i)%quad_order = coils(i)%quad_order
            tmp(i)%fil_size = coils(i)%fil_size
        end do
        nc = nc + 1
        tmp(nc)%name = name(1:min(len(name), gs_name_len))
        call move_alloc(tmp, coils)
    end subroutine append_coil

    subroutine grow_incidence(inc, nc, nk)
        !! Resize the incidence matrix to nc x nk, keeping existing entries.
        real(dp), allocatable, intent(inout) :: inc(:, :)
        integer, intent(in) :: nc, nk
        real(dp), allocatable :: tmp(:, :)
        integer :: a, b

        allocate (tmp(nc, nk))
        tmp = 0.0_dp
        a = min(nc, size(inc, 1))
        b = min(nk, size(inc, 2))
        tmp(1:a, 1:b) = inc(1:a, 1:b)
        call move_alloc(tmp, inc)
    end subroutine grow_incidence

end module kin6d_gs_coils
