program test_vrmhd_fourier
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_gs_symmetry, only: gs_symmetry_t, GS_TRANSLATION, GS_CHART_CARTESIAN
    use kin6d_vrmhd_section, only: vrmhd_section_t, vrmhd_section_build
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes
    use vrmhd_structural_spectrum, only: fourier_content
    use vrmhd_spectrum_support, only: require
    implicit none
    type(mesh_2d_t) :: mesh
    type(meridian_mesh_t) :: mm
    type(vrmhd_section_t) :: op
    real(dp), allocatable :: xy(:,:), xa(:,:), shares(:,:), captured(:)
    complex(dp), allocatable :: transform(:,:)
    complex(dp) :: modes(2,2)
    real(dp) :: length, expected(4), capture
    integer :: status, n0, a0, nn
    length = 2*acos(-1.0_dp)
    call mesh%create_rectangular(5,5,0.0_dp,length,0.0_dp,length)
    call meridian_mesh_straight(mesh,1,mm,status)
    call require(status == 0, 'Fourier oracle mesh')
    call vrmhd_section_build(mm, &
        gs_symmetry_t(GS_TRANSLATION,GS_CHART_CARTESIAN,0.0_dp),2,1,op,status)
    call require(status == 0, 'Fourier oracle section')
    call meridian_lagrange_nodes(op%scalar_fs,xy)
    call meridian_lagrange_nodes(op%fs,xa)
    n0 = op%scalar_fs%n0; a0 = 6*n0; nn = op%fs%nn
    allocate(transform(op%ndof,2))
    transform = 0; modes = 0
    modes(1,1) = 1; modes(2,2) = 1
    transform(4*n0+1:5*n0,1) = xy(1,:)
    transform(a0+nn+1:a0+nn+op%fs%n0,2) = xa(1,:)**2/2
    ! p=x and A_z=x^2/2 give p=x and B_y=-x exactly in these spaces.
    ! Integral x exp(-i m x) / (2 pi) = i/m for nonzero integer m.
    expected = [0.0_dp,0.0_dp,0.1_dp,0.9_dp]
    capture = (1.0_dp+1.0_dp/9)/(length**2/3)
    call fourier_content(op,transform,modes,4,1,2,shares,captured)
    print '(a,es20.12)', 'analytic Fourier share maximum error: ', &
        maxval(abs(shares-spread(expected,2,2)))
    print '(a,es20.12)', 'analytic captured L2 fraction maximum error: ', &
        maxval(abs(captured-capture))
    call require(maxval(abs(shares-spread(expected,2,2))) < 1.0e-8_dp, &
        'pressure and magnetic Fourier powers agree with analytic integrals')
    call require(maxval(abs(captured-capture)) < 1.0e-8_dp, &
        'captured power uses the full physical L2 norm')
end program
