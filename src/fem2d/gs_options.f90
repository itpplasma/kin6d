module kin6d_gs_options
    ! Run-time switches for the Grad--Shafranov linear algebra. The defaults are
    ! the fast path (CSR operators, cached Cholesky factor). The old element-wise
    ! operator and the Jacobi-CG solvers stay selectable for comparison:
    !   KIN6D_GS_OPERATOR = csr (default) | element
    !   KIN6D_GS_LINEAR   = direct (default when MUMPS is built) | iterative
    ! The variables are read once, on first use, outside parallel regions.
    implicit none
    private
    public :: gs_use_csr, gs_use_direct, gs_set_options, gs_options_text
    logical :: initialised = .false.
    logical :: csr_on = .true., direct_on = .true.
contains
    subroutine read_environment()
        character(32) :: value
        integer :: length
        if (initialised) return
        initialised = .true.
        call get_environment_variable('KIN6D_GS_OPERATOR', value, length)
        if (length > 0) then
            select case (trim(value))
            case ('csr'); csr_on = .true.
            case ('element'); csr_on = .false.
            case default
                error stop 'KIN6D_GS_OPERATOR must be csr (default) or element'
            end select
        end if
        call get_environment_variable('KIN6D_GS_LINEAR', value, length)
        if (length > 0) then
            select case (trim(value))
            case ('direct'); direct_on = .true.
            case ('iterative'); direct_on = .false.
            case default
                error stop 'KIN6D_GS_LINEAR must be direct (default) or iterative'
            end select
        end if
    end subroutine

    logical function gs_use_csr() result(on)
        call read_environment()
        on = csr_on
    end function

    logical function gs_use_direct() result(on)
        call read_environment()
        on = direct_on
    end function

    subroutine gs_set_options(csr, direct)
        ! Programmatic override (tests); the environment is not read afterwards.
        logical, intent(in) :: csr, direct
        initialised = .true.
        csr_on = csr
        direct_on = direct
    end subroutine

    function gs_options_text() result(text)
        character(:), allocatable :: text
        call read_environment()
        if (csr_on) then
            text = 'operator=csr'
        else
            text = 'operator=element'
        end if
        if (direct_on) then
            text = text // ' linear=direct'
        else
            text = text // ' linear=iterative'
        end if
    end function
end module
