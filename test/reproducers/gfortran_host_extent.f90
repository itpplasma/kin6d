module gfortran_host_extent
    implicit none
contains
    subroutine matcher(n)
        integer, intent(in) :: n
        real :: oracle(n)
        logical :: candidates(size(oracle))
        oracle = 1
        candidates = .false.
        print *, 'expected extents: 3 3; actual:', size(oracle), size(candidates)
        candidates(1) = .true.
        if (.not. assigned()) error stop 'host association'
    contains
        logical function assigned()
            assigned = candidates(1)
        end function
    end subroutine
end module

program reproduce_host_extent
    use gfortran_host_extent, only: matcher
    implicit none
    call matcher(3)
end program
