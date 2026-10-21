module kin6d_s0_spectrum
    ! S0 consumer of the existing MHD energy and physical inertia operators.
    use, intrinsic :: iso_fortran_env, only: dp => real64, i8 => int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_eigen, only: lanczos, krylov_schur, eigen_options_t, &
        eigen_stats_t, EIGEN_OK, EIGEN_INVALID
    use fortsparse, only: csc_z_t, csc_is_valid, csc_conjugate_transpose, &
        csc_matvec, sparse_solver_t, sparse_factor, sparse_solve, sparse_destroy, &
        sparse_inertia, sparse_factor_stats, sparse_factor_stats_t, &
        fortsparse_status_t, FORTSPARSE_OK, FORTSPARSE_BACKEND_NATIVE, &
        FORTSPARSE_SYM_SPD, FORTSPARSE_SYM_HERMITIAN, FORTSPARSE_SYM_GENERAL, &
        FORTSPARSE_LU_SUPERNODAL, FORTSPARSE_ORDER_AMD, native_ordering_available, &
        FORTSPARSE_INERTIA_UNAVAILABLE
    implicit none
    private
    public :: s0_pencil_t, s0_report_t, s0_prepare, s0_solve, s0_destroy
    type :: s0_pencil_t
        private
        type(csc_z_t) :: stiffness, mass
        type(sparse_solver_t) :: mass_factor, shift_factor
        complex(dp) :: shift = 0
        logical :: ready = .false.
        integer(i8) :: below = -1, at = -1, above = -1
        integer :: inertia_status = FORTSPARSE_INERTIA_UNAVAILABLE
        real(dp) :: stiffness_scale = 0, mass_scale = 0
    end type s0_pencil_t
    type :: s0_report_t
        integer :: status = EIGEN_INVALID, sparse_status = FORTSPARSE_OK
        integer(i8) :: below = -1, at = -1, above = -1
        integer :: inertia_status = FORTSPARSE_INERTIA_UNAVAILABLE
        type(eigen_stats_t) :: eigen
        type(sparse_factor_stats_t) :: factor
        character(16) :: claim = 'candidate', bound_status = 'unavailable'
        character(120) :: message = 'invalid S0 input or unprepared pencil'
    end type s0_report_t
contains
    subroutine s0_prepare(pencil, stiffness, mass, shift, report)
        type(s0_pencil_t), intent(inout) :: pencil
        type(csc_z_t), intent(in) :: stiffness, mass
        complex(dp), intent(in) :: shift
        type(s0_report_t), intent(out) :: report
        type(csc_z_t) :: shifted
        type(fortsparse_status_t) :: status
        integer :: j
        call s0_destroy(pencil)
        report = s0_report_t()
        if (.not. ieee_is_finite(real(shift))) return
        if (.not. ieee_is_finite(aimag(shift))) return
        if (.not. hermitian(stiffness)) return
        if (.not. hermitian(mass)) return
        if (stiffness%nrow /= mass%nrow .or. stiffness%nnz /= mass%nnz) return
        if (any(stiffness%col_ptr /= mass%col_ptr)) return
        if (any(stiffness%row_idx /= mass%row_idx)) return
        ! Same full CSC pattern, including structural zeros, is the assembly
        ! contract. No dense logical pairwise matrix or Jacobian is formed.
        pencil%stiffness = stiffness; pencil%mass = mass; pencil%shift = shift
        pencil%mass_factor%backend_id = FORTSPARSE_BACKEND_NATIVE
        pencil%shift_factor%backend_id = FORTSPARSE_BACKEND_NATIVE
        pencil%mass_factor%options%symmetry = FORTSPARSE_SYM_SPD
        call sparse_factor(pencil%mass_factor,mass,status)
        report%sparse_status = status%code
        report%message = status%msg
        if (status%code /= FORTSPARSE_OK) return
        shifted = stiffness
        shifted%val = stiffness%val-shift*mass%val
        if (abs(aimag(shift)) <= 0) then
            pencil%shift_factor%options%symmetry = FORTSPARSE_SYM_HERMITIAN
        else
            pencil%shift_factor%options%symmetry = FORTSPARSE_SYM_GENERAL
            pencil%shift_factor%options%lu_method = FORTSPARSE_LU_SUPERNODAL
            if (native_ordering_available(FORTSPARSE_ORDER_AMD)) &
                pencil%shift_factor%options%ordering = FORTSPARSE_ORDER_AMD
        end if
        call sparse_factor(pencil%shift_factor,shifted,status)
        report%sparse_status = status%code
        report%message = status%msg
        call sparse_factor_stats(pencil%shift_factor,report%factor,status)
        if (abs(aimag(shift)) <= 0) then
            call sparse_inertia(pencil%shift_factor,pencil%below,pencil%at, &
                pencil%above,pencil%inertia_status,status)
            report%below = pencil%below; report%at = pencil%at
            report%above = pencil%above; report%inertia_status = pencil%inertia_status
        end if
        if (report%sparse_status /= FORTSPARSE_OK) return
        ! ||A||_2 <= ||A||_1 for Hermitian A: a supplied, independently
        ! computable scale for LE's normwise stopping criterion.
        do j = 1, mass%ncol
            pencil%stiffness_scale = max(pencil%stiffness_scale, &
                sum(abs(stiffness%val(stiffness%col_ptr(j): &
                    stiffness%col_ptr(j+1)-1))))
            pencil%mass_scale = max(pencil%mass_scale, &
                sum(abs(mass%val(mass%col_ptr(j):mass%col_ptr(j+1)-1))))
        end do
        pencil%ready = .true.
        report%status = EIGEN_OK
    end subroutine s0_prepare

    subroutine s0_solve(pencil, start, tolerance, options, omega2, vectors, &
            residuals, report)
        type(s0_pencil_t), intent(inout) :: pencil
        complex(dp), intent(in) :: start(:)
        real(dp), intent(in) :: tolerance
        type(eigen_options_t), intent(in) :: options
        complex(dp), intent(out) :: omega2(:), vectors(:,:)
        real(dp), intent(out) :: residuals(:)
        type(s0_report_t), intent(out) :: report
        type(eigen_options_t) :: selected
        type(fortsparse_status_t) :: status
        real(dp) :: real_values(size(omega2))
        report = s0_report_t()
        omega2 = 0; vectors = 0; residuals = huge(1.0_dp)
        if (.not. pencil%ready) return
        if (size(start) /= pencil%mass%nrow) return
        report%message = ''
        report%below = pencil%below; report%at = pencil%at
        report%above = pencil%above; report%inertia_status = pencil%inertia_status
        call sparse_factor_stats(pencil%shift_factor,report%factor,status)
        selected = options
        selected%operator_scale = pencil%stiffness_scale
        selected%mass_scale = pencil%mass_scale
        if (abs(aimag(pencil%shift)) <= 0) then
            call lanczos(apply_k,start,real_values,vectors,residuals,tolerance, &
                selected,report%status,report%eigen,mass=apply_m, &
                solve_mass=solve_m,solve_shift=solve_shift, &
                sigma=real(pencil%shift,dp))
            omega2 = cmplx(real_values,0.0_dp,dp)
        else
            call krylov_schur(apply_k,start,omega2,vectors,residuals,tolerance, &
                selected,report%status,report%eigen,mass=apply_m, &
                solve_mass=solve_m,solve_shift=solve_shift,sigma=pencil%shift)
        end if
    contains
        subroutine apply_k(x,y)
            complex(dp), intent(in) :: x(:)
            complex(dp), intent(out) :: y(:)
            y = csc_matvec(pencil%stiffness,x)
        end subroutine apply_k
        subroutine apply_m(x,y)
            complex(dp), intent(in) :: x(:)
            complex(dp), intent(out) :: y(:)
            y = csc_matvec(pencil%mass,x)
        end subroutine apply_m
        subroutine solve_m(x,y,info)
            complex(dp), intent(in) :: x(:)
            complex(dp), intent(out) :: y(:)
            integer, intent(out) :: info
            type(fortsparse_status_t) :: status
            call sparse_solve(pencil%mass_factor,x,y,status)
            info = status%code
        end subroutine solve_m
        subroutine solve_shift(x,y,info)
            complex(dp), intent(in) :: x(:)
            complex(dp), intent(out) :: y(:)
            integer, intent(out) :: info
            type(fortsparse_status_t) :: status
            call sparse_solve(pencil%shift_factor,x,y,status)
            info = status%code
        end subroutine solve_shift
    end subroutine s0_solve

    logical function hermitian(a) result(ok)
        type(csc_z_t), intent(in) :: a
        type(csc_z_t) :: adjoint
        type(fortsparse_status_t) :: status
        ok = .false.
        if (.not. csc_is_valid(a)) return
        if (a%nrow /= a%ncol .or. a%nrow < 1) return
        if (.not. all(ieee_is_finite(real(a%val)))) return
        if (.not. all(ieee_is_finite(aimag(a%val)))) return
        call csc_conjugate_transpose(a,adjoint,status)
        if (status%code /= FORTSPARSE_OK) return
        if (any(a%col_ptr /= adjoint%col_ptr)) return
        if (any(a%row_idx /= adjoint%row_idx)) return
        ! A^H is already conjugated by the library.
        ok = maxval(abs(a%val-adjoint%val)) <= &
            64*epsilon(1.0_dp)*max(tiny(1.0_dp),maxval(abs(a%val)))
    end function hermitian

    subroutine s0_destroy(pencil)
        type(s0_pencil_t), intent(inout) :: pencil
        call sparse_destroy(pencil%mass_factor)
        call sparse_destroy(pencil%shift_factor)
        pencil%ready = .false.
        pencil%below = -1; pencil%at = -1; pencil%above = -1
        pencil%inertia_status = FORTSPARSE_INERTIA_UNAVAILABLE
        pencil%stiffness_scale = 0; pencil%mass_scale = 0
    end subroutine s0_destroy
end module kin6d_s0_spectrum
