module kin6d_s0_volume
    ! Element-local Galerkin assembly; topology and physical jets belong to caller.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_ideal_xi, only: xi_bilinear
    use kin6d_s0_leaves, only: s0_mass_dp
    use fortsparse, only: csc_z_t, csc_from_triplet, fortsparse_status_t, &
        FORTSPARSE_OK, status_set, FORTSPARSE_INVALID_MATRIX
    implicit none
    private
    public :: s0_assemble_volume
contains
    subroutine s0_assemble_volume(ndof, dof, weight, b, p, rho, mu, value, &
            gradient, stiffness, mass, status)
        integer, intent(in) :: ndof, dof(:,:)
        real(dp), intent(in) :: weight(:,:), b(:,:,:), p(:,:), rho(:,:), mu
        complex(dp), intent(in) :: value(:,:,:,:), gradient(:,:,:,:,:)
        type(csc_z_t), intent(out) :: stiffness, mass
        type(fortsparse_status_t), intent(out) :: status
        integer, allocatable :: rows(:), cols(:)
        complex(dp), allocatable :: kval(:), mval(:)
        complex(dp) :: ke, me
        real(dp) :: rr, ii, ir, ri
        integer :: nb, nq, ne, e, i, j, q, at, mirror, ntrip
        call status_set(status,FORTSPARSE_INVALID_MATRIX,'invalid S0 volume input')
        nb = size(dof,1); ne = size(dof,2); nq = size(weight,1)
        if (min(ndof,nb,ne,nq) < 1) return
        if (any(dof < 0) .or. any(dof > ndof)) return
        if (any(shape(weight) /= [nq,ne])) return
        if (any(shape(b) /= [3,nq,ne])) return
        if (any(shape(p) /= [nq,ne]) .or. any(shape(rho) /= [nq,ne])) return
        if (any(shape(value) /= [3,nb,nq,ne])) return
        if (any(shape(gradient) /= [3,3,nb,nq,ne])) return
        if (.not. ieee_is_finite(mu)) return
        if (mu <= 0) return
        if (.not. all(ieee_is_finite(weight))) return
        if (.not. all(ieee_is_finite(b))) return
        if (.not. all(ieee_is_finite(p))) return
        if (.not. all(ieee_is_finite(rho))) return
        if (any(weight <= 0) .or. any(rho <= 0)) return
        if (.not. all(ieee_is_finite(real(value,dp)))) return
        if (.not. all(ieee_is_finite(aimag(value)))) return
        if (.not. all(ieee_is_finite(real(gradient,dp)))) return
        if (.not. all(ieee_is_finite(aimag(gradient)))) return
        if (real(nb,dp)*nb*ne > real(huge(ntrip),dp)) return
        ntrip = nb*nb*ne
        allocate(rows(ntrip),cols(ntrip),kval(ntrip),mval(ntrip))
        rows = 1; cols = 1; kval = 0; mval = 0
        ! One triangle evaluated and reflected before global summation. No
        ! allocation or mutable shared scratch in the element loop.
        !$omp parallel do default(none) shared(nb,ne,nq,dof,weight,b,p,rho,mu, &
        !$omp value,gradient,rows,cols,kval,mval) &
        !$omp private(e,i,j,q,at,mirror,ke,me,rr,ii,ir,ri)
        do e = 1, ne
            do j = 1, nb
                if (dof(j,e) == 0) cycle
                do i = 1, j
                    if (dof(i,e) == 0) cycle
                    ke = 0; me = 0
                    do q = 1, nq
                        ke = ke+weight(q,e)*xi_bilinear(b(:,q,e),p(q,e),mu, &
                            gradient(:,:,j,q,e),gradient(:,:,i,q,e))
                        call s0_mass_dp(rho(q,e), &
                            real(value(1,j,q,e)),real(value(2,j,q,e)), &
                            real(value(3,j,q,e)),real(value(1,i,q,e)), &
                            real(value(2,i,q,e)),real(value(3,i,q,e)),rr)
                        call s0_mass_dp(rho(q,e), &
                            aimag(value(1,j,q,e)),aimag(value(2,j,q,e)), &
                            aimag(value(3,j,q,e)),aimag(value(1,i,q,e)), &
                            aimag(value(2,i,q,e)),aimag(value(3,i,q,e)),ii)
                        call s0_mass_dp(rho(q,e), &
                            aimag(value(1,j,q,e)),aimag(value(2,j,q,e)), &
                            aimag(value(3,j,q,e)),real(value(1,i,q,e)), &
                            real(value(2,i,q,e)),real(value(3,i,q,e)),ir)
                        call s0_mass_dp(rho(q,e), &
                            real(value(1,j,q,e)),real(value(2,j,q,e)), &
                            real(value(3,j,q,e)),aimag(value(1,i,q,e)), &
                            aimag(value(2,i,q,e)),aimag(value(3,i,q,e)),ri)
                        me = me+weight(q,e)*cmplx(rr+ii,ir-ri,dp)
                    end do
                    at = (e-1)*nb*nb+(j-1)*nb+i
                    mirror = (e-1)*nb*nb+(i-1)*nb+j
                    rows(at) = dof(i,e); cols(at) = dof(j,e)
                    rows(mirror) = dof(j,e); cols(mirror) = dof(i,e)
                    kval(at) = ke; mval(at) = me
                    kval(mirror) = conjg(ke); mval(mirror) = conjg(me)
                end do
            end do
        end do
        !$omp end parallel do
        if (.not. all(ieee_is_finite(real(kval,dp)))) return
        if (.not. all(ieee_is_finite(aimag(kval)))) return
        if (.not. all(ieee_is_finite(real(mval,dp)))) return
        if (.not. all(ieee_is_finite(aimag(mval)))) return
        call csc_from_triplet(ndof,ndof,rows,cols,kval,stiffness,status)
        if (status%code /= FORTSPARSE_OK) return
        call csc_from_triplet(ndof,ndof,rows,cols,mval,mass,status)
    end subroutine s0_assemble_volume
end module kin6d_s0_volume
