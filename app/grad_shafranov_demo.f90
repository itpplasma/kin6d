program grad_shafranov_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
                                   solve_gs_profiles, gs_profile_integral, write_gs_mesh
   use kin6d_result_io, only: run_t, result_file_t
   use kin6d_grad_shafranov_leaves, only: gs_field
   use fortplot_system_runtime, only: create_directory_runtime
   implicit none
   type(gs_state_t) :: state
   type(run_t) :: run
   type(result_file_t) :: result
   integer :: boundary_nodes = 128, unit, status, i
   real(dp) :: r0 = 6.2_dp, a = 0.62_dp, kappa = 1.0_dp, delta = 0.0_dp
   real(dp) :: b0 = 5.3_dp, q_seed = 1.5_dp, max_area = 0.002_dp
   real(dp) :: mu0 = 1.2566370614359173e-6_dp, pprime, ffprime, theta
   real(dp), allocatable :: contour(:, :)
   real(dp), allocatable :: flux(:), pp(:), ff(:), pressure(:), fpol(:)
   real(dp), allocatable :: nodal_p(:), nodal_f(:), element_r(:), element_z(:)
   real(dp), allocatable :: element_br(:), element_bz(:), element_bt(:)
   real(dp) :: element_psi, element_f, gradient(2)
   character(16) :: pressure_mode = 'zero_beta'
   character(512) :: mesh_file = 'output/grad_shafranov_mesh.dat'
   character(512) :: boundary_file = ''
   character(512) :: profile_file = ''
   integer :: profile_count
   logical :: directory_ok
   character(256) :: message
   namelist /grad_shafranov/ r0, a, kappa, delta, b0, q_seed, max_area, &
      boundary_nodes, pressure_mode, mesh_file, boundary_file, profile_file
   call run%start('grad_shafranov')
   status = 0; message = ''
   call run%open_config(unit)
   if (run%has_config()) read (unit, nml=grad_shafranov, &
                               iostat=status, iomsg=message)
   call run%close_config(unit, 'grad_shafranov', status, message)
   if (boundary_nodes < 8 .or. r0 <= a .or. a <= 0.0_dp) &
      error stop 'invalid GS fixed boundary'
    if (q_seed == 0.0_dp) error stop 'GS source q_seed cannot vanish'
    if (abs(delta) > 1 .or. kappa <= 0) error stop 'invalid GS shaping'
   if (len_trim(boundary_file) > 0) then
      open (newunit=unit, file=trim(boundary_file), status='old')
      read (unit, *) boundary_nodes
      allocate (contour(2, boundary_nodes))
      do i = 1, boundary_nodes
         read (unit, *) contour(:, i)
      end do
      close (unit)
   else
      allocate (contour(2, boundary_nodes))
      do i = 1, boundary_nodes
         theta = 2*acos(-1.0_dp)*real(i - 1, dp)/boundary_nodes
         contour(:, i) = [r0 + a*cos(theta + asin(delta)*sin(theta)), &
                          kappa*a*sin(theta)]
      end do
   end if
   select case (trim(pressure_mode))
   case ('zero_beta')
      pprime = 0.0_dp; ffprime = -2*b0/q_seed
   case ('finite_beta')
      pprime = -2*b0/(q_seed*mu0*r0**2); ffprime = 0.0_dp
   case default
      error stop 'unsupported GS source model'
   end select
   call run%open_resolved(unit)
   write (unit, nml=grad_shafranov)
   call run%close_resolved(unit)
   call initialize_gs(state, contour, max_area, status)
   if (status /= 0) error stop 'GS meshing failed'
   if (len_trim(profile_file) > 0) then
      open (newunit=unit, file=trim(profile_file), status='old')
      read (unit, *) profile_count
      if (profile_count < 2) error stop 'GS profile requires >=2 flux points'
      allocate (flux(profile_count), pp(profile_count), ff(profile_count), &
                pressure(profile_count), fpol(profile_count))
      do i = 1, profile_count
         read (unit, *) flux(i), pp(i), ff(i), pressure(i), fpol(i)
      end do
        close (unit)
        if (any(.not. ieee_is_finite(pressure)) .or. &
            any(.not. ieee_is_finite(fpol))) error stop 'nonfinite GS p/F reference'
        ! Nontrivial initialization uses the axis source, not the known answer.
      ! A zero initial state can otherwise select a trivial edge-source root.
      pprime = pp(1); ffprime = ff(1)
      call solve_gs(state, mu0, pp(1), ff(1), status)
      if (status == 0) call solve_gs_profiles(state, mu0, flux, pp, ff, status)
   else
      call solve_gs(state, mu0, pprime, ffprime, status)
   end if
   if (status /= 0) then
      print *, 'GS failure status/iterations/residual:', status, &
         state%iterations, state%residual
      error stop 'GS candidate solve failed'
   end if
   allocate (nodal_p(size(state%psi)), nodal_f(size(state%psi)), &
             element_r(size(state%triangles, 2)), element_z(size(state%triangles, 2)), &
           element_br(size(state%triangles, 2)), element_bz(size(state%triangles, 2)), &
             element_bt(size(state%triangles, 2)))
   do i = 1, size(state%psi)
      if (len_trim(profile_file) > 0) then
         nodal_p(i) = pressure(profile_count) + &
                      gs_profile_integral(flux, pp, state%psi(i))
         nodal_f(i) = reconstructed_f(state%psi(i))
      else
         nodal_p(i) = pprime*state%psi(i)
         if ((r0*b0)**2 + 2*ffprime*state%psi(i) <= 0) &
            error stop 'GS toroidal flux function square is nonpositive'
         nodal_f(i) = sign(sqrt((r0*b0)**2 + 2*ffprime*state%psi(i)), b0)
      end if
   end do
   do i = 1, size(element_r)
      element_r(i) = sum(state%points(1, state%triangles(:, i)))/3
      element_z(i) = sum(state%points(2, state%triangles(:, i)))/3
      element_psi = sum(state%psi(state%triangles(:, i)))/3
      if (len_trim(profile_file) > 0) then
         element_f = reconstructed_f(element_psi)
      else
         element_f = sign(sqrt((r0*b0)**2 + 2*ffprime*element_psi), b0)
      end if
      gradient = matmul(state%gradient(:, :, i), &
                        state%psi(state%triangles(:, i)))
      call gs_field(element_r(i), gradient(1), gradient(2), element_f, &
                    element_br(i), element_bz(i), element_bt(i))
   end do
   call create_directory_runtime('output', directory_ok)
   if (.not. directory_ok) error stop 'GS output directory unavailable'
   call write_gs_mesh(state, trim(mesh_file), status)
   if (status /= 0) error stop 'GS mesh export failed'
   call result%create(run)
   call result%attribute('claim', 'candidate')
   call result%attribute('equations', 'Delta*psi=-mu0*R^2*pprime-FFprime')
   call result%attribute('conventions', 'SI; phi CCW; psi flux per radian')
   call result%attribute('geometry_error', 'polygonal LCFS; bound unavailable')
   call result%attribute('continuous_error', 'bound unavailable; refine independently')
   call result%attribute('profile_model', &
          'constant source or piecewise-linear flux table; constant endpoint extension')
    call result%attribute('fieldlaw', &
        'derivative-controlled reconstruction; edge p/F; signed F branch')
    call result%attribute('fortfem_revision', &
        '2d7baa19fe1bef02ce279fa407830c8731bdc663')
   call result%write('R', state%points(1, :), 'm', 'cylindrical radius')
   call result%write('Z', state%points(2, :), 'm', 'vertical position')
   call result%write('psi', state%psi, 'Wb/rad', 'poloidal flux per radian')
   call result%write('pressure', nodal_p, 'Pa', 'pressure at nodes')
   call result%write('F', nodal_f, 'T m', 'toroidal flux function at nodes')
   call result%write('element_R', element_r, 'm', 'triangle centroid R')
   call result%write('element_Z', element_z, 'm', 'triangle centroid Z')
   call result%write('element_B_R', element_br, 'T', 'radial B at centroids')
   call result%write('element_B_Z', element_bz, 'T', 'vertical B at centroids')
   call result%write('element_B_phi', element_bt, 'T', 'CCW toroidal B at centroids')
   call result%write('triangles', real(state%triangles, dp), '1', &
                     'one-based triangle node indices')
   call result%write('pprime', pprime, 'Pa/(Wb/rad)', 'constant pressure derivative')
   call result%write('ffprime', ffprime, 'T', 'constant F*dF/dpsi')
   if (len_trim(profile_file) > 0) then
      call result%write('F_edge', fpol(size(fpol)), 'T m', &
                        'toroidal flux function at edge')
   else
      call result%write('F_edge', r0*b0, 'T m', 'toroidal flux function at edge')
   end if
   if (len_trim(profile_file) > 0) then
      call result%write('profile_flux', flux, 'Wb/rad', 'prescribed profile flux')
      call result%write('profile_pprime', pp, 'Pa/(Wb/rad)', 'pressure derivative')
      call result%write('profile_ffprime', ff, 'T', 'F*dF/dpsi')
      call result%write('profile_pressure', pressure, 'Pa', 'pressure data')
      call result%write('profile_F', fpol, 'T m', 'toroidal flux function data')
      open (newunit=unit, file=trim(mesh_file)//'.profiles.dat', &
            status='replace')
      write (unit, *) profile_count
      do i = 1, profile_count
         write (unit, '(5ES25.16)') flux(i), pp(i), ff(i), &
            pressure(profile_count) + gs_profile_integral(flux, pp, flux(i)), &
            reconstructed_f(flux(i))
      end do
      close (unit)
   end if
    call result%write('residual', state%residual, '1', 'relative weak equation residual')
   call result%close()
   print *, 'GS candidate nodes/triangles/residual:', size(state%psi), &
      size(state%triangles, 2), state%residual
   print *, 'q_seed is a source parameter, not a prescribed q profile.'
contains
   real(dp) function reconstructed_f(psi) result(value)
      real(dp), intent(in) :: psi
      real(dp) :: square
      square = fpol(profile_count)**2 + 2*gs_profile_integral(flux, ff, psi)
      if (square <= 0 .or. fpol(profile_count) == 0) &
         error stop 'GS reconstructed F branch is undefined'
      value = sign(sqrt(square), fpol(profile_count))
   end function reconstructed_f
end program grad_shafranov_demo
