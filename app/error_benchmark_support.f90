module error_benchmark_support
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: row, header, require
contains
    subroutine header()
        print '(a)', 'case,tier,trial,solve_seconds,estimator_seconds,'// &
            'certificate_seconds,total_seconds,radius,scale,radius_over_scale,'// &
            'budget,status,inputs,timing_method'
    end subroutine

    subroutine row(name, tier, trial, solve, estimate, certificate, radius, &
                   scale, status, inputs, method)
        character(*), intent(in) :: name, tier, status, inputs, method
        integer, intent(in) :: trial
        real(dp), intent(in) :: solve, estimate, certificate, radius, scale
        print '(*(g0,:,","))', name, tier, trial, solve, estimate, certificate, &
            solve + estimate + certificate, radius, scale, radius/scale, &
            1.e-12_dp + .01_dp*scale, status, inputs, method
    end subroutine

    subroutine require(ok)
        logical, intent(in) :: ok
        if (.not. ok) error stop 'error-control benchmark premise failed'
    end subroutine
end module
