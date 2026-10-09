program grad_shafranov_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
      solve_gs_profiles, gs_profile_integral, write_gs_mesh, gs_element_jet, &
      gs_classical_current, gs_current_jump, gs_axis_candidate, gs_axis_safety_factor
   use kin6d_gs_inverse, only: gs_inverse_result_t, solve_gs_constant_q
   use kin6d_result_io, only: run_t, result_file_t
   use kin6d_grad_shafranov_leaves, only: gs_field
   use fortplot_system_runtime, only: create_directory_runtime
   implicit none
   type(gs_state_t) :: state
   type(gs_inverse_result_t)::inverse
   type(run_t) :: run
   type(result_file_t) :: result
   integer :: boundary_nodes = 128, element_order = 1, unit, status, i
   real(dp) :: r0 = 6.2_dp, a = 0.62_dp, kappa = 1.0_dp, delta = 0.0_dp
   real(dp) :: b0 = 5.3_dp, q_seed = 1.5_dp, max_area = 0.002_dp
   real(dp) :: nonlinear_tolerance = 1.0e-8_dp
   character(32)::solve_mode='forward'
   real(dp)::q_target=1.5_dp,toroidal_flux_target=0,q_tolerance=1.0e-3_dp,flux_tolerance=1.0e-6_dp
   integer::inverse_flux_points=9,inverse_angular_points=256
   logical::has_profiles, prescribed_q
   real(dp) :: mu0 = 1.2566370614359173e-6_dp, pprime, ffprime, theta
   real(dp), allocatable :: contour(:, :)
   real(dp), allocatable :: flux(:), pp(:), ff(:), pressure(:), fpol(:)
   real(dp), allocatable :: nodal_p(:), nodal_f(:), element_r(:), element_z(:)
   real(dp), allocatable :: element_br(:), element_bz(:), element_bt(:)
   real(dp), allocatable :: element_current(:), sheet_r(:), sheet_z(:), sheet_current(:)
   real(dp) :: element_psi, element_f, gradient(2), hessian(2, 2)
   real(dp) :: axis(2), axis_flux, axis_q, normal(2), edge_point(2), toward(2)
   integer :: edge, neighbor, face, face_count, axis_element, axis_status
   character(16) :: pressure_mode = 'zero_beta'
   character(512) :: mesh_file = 'output/grad_shafranov_mesh.dat'
   character(512) :: boundary_file = ''
   character(512) :: profile_file = ''
   integer :: profile_count
   logical :: directory_ok
   character(256) :: message
   namelist /grad_shafranov/ r0, a, kappa, delta, b0, q_seed, max_area, &
      boundary_nodes, element_order, pressure_mode, mesh_file, boundary_file, profile_file, &
      nonlinear_tolerance,solve_mode,q_target,toroidal_flux_target,q_tolerance,flux_tolerance, &
      inverse_flux_points,inverse_angular_points
   call run%start('grad_shafranov')
   status = 0; message = ''
   call run%open_config(unit)
   if (run%has_config()) read (unit, nml=grad_shafranov, &
                               iostat=status, iomsg=message)
   call run%close_config(unit, 'grad_shafranov', status, message)
   if (boundary_nodes < 8 .or. r0 <= a .or. a <= 0.0_dp) &
      error stop 'invalid GS fixed boundary'
    if (q_seed == 0.0_dp.and.trim(solve_mode)=='forward') &
        error stop 'GS source q_seed cannot vanish'
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
   prescribed_q=trim(solve_mode)=='prescribed_constant_q'
   if(trim(solve_mode)/='forward'.and..not.prescribed_q)error stop 'invalid GS solve mode'
   if(prescribed_q.and.(trim(pressure_mode)/='zero_beta'.or.len_trim(profile_file)>0)) &
       error stop 'prescribed-q mode requires zero pressure and no forward profile input'
   has_profiles=len_trim(profile_file)>0.or.prescribed_q
   if(prescribed_q)then
      pprime=0; ffprime=0
   else
   select case (trim(pressure_mode))
   case ('zero_beta')
      pprime = 0.0_dp; ffprime = -2*b0/q_seed
   case ('finite_beta')
      pprime = -2*b0/(q_seed*mu0*r0**2); ffprime = 0.0_dp
   case default
      error stop 'unsupported GS source model'
   end select
   end if
   call run%open_resolved(unit)
   write (unit, nml=grad_shafranov)
   call run%close_resolved(unit)
   call initialize_gs(state, contour, max_area, status, element_order)
   if (status /= 0) error stop 'GS meshing failed'
   if(prescribed_q)then
      call solve_gs_constant_q(state,q_target,toroidal_flux_target,inverse,status, &
          inverse_flux_points,inverse_angular_points,nonlinear_tolerance,q_tolerance,flux_tolerance)
      if(status==0)then
         profile_count=size(inverse%flux)
         flux=inverse%flux; ff=inverse%ffprime
         allocate(pp(profile_count),pressure(profile_count),fpol(profile_count))
         pp=0; pressure=0
         fpol=sign(sqrt(inverse%f_squared),q_target)
         pprime=0; ffprime=0
      end if
   else if (len_trim(profile_file) > 0) then
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
      if (status == 0) call solve_gs_profiles(state, mu0, flux, pp, ff, status, &
                                             nonlinear_tolerance)
   else
      call solve_gs(state, mu0, pprime, ffprime, status)
   end if
   if (status /= 0) then
      if(prescribed_q)print *, 'inverse failure discrete indicators:',inverse%iterations, &
          inverse%coupled_residual,inverse%knot_indicator,inverse%profile_indicator, &
          inverse%angular_indicator,inverse%flux_indicator
      print *, 'GS failure status/iterations/residual:', status, &
         state%iterations, state%residual
      error stop 'GS candidate solve failed'
   end if
   allocate (nodal_p(size(state%psi)), nodal_f(size(state%psi)), &
             element_r(size(state%triangles, 2)), element_z(size(state%triangles, 2)), &
           element_br(size(state%triangles, 2)), element_bz(size(state%triangles, 2)), &
             element_bt(size(state%triangles, 2)), &
             element_current(size(state%triangles, 2)))
   do i = 1, size(state%psi)
      if (has_profiles) then
         nodal_p(i) = pressure(profile_count) + &
                      gs_profile_integral(flux, pp, state%psi(i))
         nodal_f(i) = reconstructed_f(state%psi(i))
      else
         nodal_p(i) = pprime*state%psi(i)
         nodal_f(i) = constant_source_f(state%psi(i))
      end if
   end do
   do i = 1, size(element_r)
      element_r(i) = sum(state%points(1, state%triangles(1:3, i)))/3
      element_z(i) = sum(state%points(2, state%triangles(1:3, i)))/3
      call gs_element_jet(state, i, element_r(i), element_z(i), &
                          element_psi, gradient, hessian, status)
      if (status /= 0) error stop 'GS element field evaluation failed'
      if (has_profiles) then
         element_f = reconstructed_f(element_psi)
      else
         element_f = constant_source_f(element_psi)
      end if
      call gs_field(element_r(i), gradient(1), gradient(2), element_f, &
                    element_br(i), element_bz(i), element_bt(i))
      call gs_classical_current(state, i, element_r(i), element_z(i), &
                                mu0, element_current(i), status)
      if (status /= 0) error stop 'GS classical current evaluation failed'
   end do
   face_count = count(state%neighbors > 0)/2
   allocate (sheet_r(face_count), sheet_z(face_count), sheet_current(face_count))
   face = 0
   do i = 1, size(state%triangles, 2)
      do edge = 1, 3
         neighbor = state%neighbors(edge, i)
         if (neighbor <= i) cycle
         face = face+1
         edge_point = state%points(:, state%triangles(edge, i))
         normal = state%points(:, state%triangles(mod(edge, 3)+1, i))-edge_point
         edge_point = edge_point+0.5_dp*normal
         normal = [normal(2), -normal(1)]/norm2(normal)
         toward = sum(state%points(:, state%triangles(1:3, neighbor)), dim=2) &
                - sum(state%points(:, state%triangles(1:3, i)), dim=2)
         if (dot_product(normal, toward) < 0) normal = -normal
         sheet_r(face) = edge_point(1); sheet_z(face) = edge_point(2)
         call gs_current_jump(state, i, neighbor, edge_point(1), edge_point(2), &
                              normal, mu0, sheet_current(face), status)
         if (status /= 0) error stop 'GS current sheet evaluation failed'
      end do
   end do
   call gs_axis_candidate(state, axis, axis_flux, hessian, axis_element, axis_status)
   if (axis_status == 0) then
      if (has_profiles) then
         element_f = reconstructed_f(axis_flux)
      else
         element_f = constant_source_f(axis_flux)
      end if
      call gs_axis_safety_factor(axis(1), element_f, hessian, axis_q, axis_status)
   end if
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
   call result%attribute('finite_element', &
       merge('P2 continuous quadratic', 'P1 continuous affine   ', element_order == 2))
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
   call result%write('element_J_phi', element_current, 'A m-2', &
                     'classical CCW toroidal current inside each element')
   if (face_count > 0) then
      call result%write('sheet_R', sheet_r, 'm', 'interior face midpoint R')
      call result%write('sheet_Z', sheet_z, 'm', 'interior face midpoint Z')
      call result%write('sheet_K_phi', sheet_current, 'A m-1', &
                        'distributional CCW current sheet at face midpoints')
   end if
   call result%write('smooth_axis_status', real(axis_status, dp), '1', &
                     '0 regular candidate; 3 P1; 4 unqualified minimum; 5 nonunique')
   if (axis_status == 0) then
      call result%write('axis_RZ', axis, 'm', 'stationary interior P2 axis candidate')
      call result%write('axis_psi', axis_flux, 'Wb/rad', 'P2 flux at axis candidate')
      call result%write('axis_q', axis_q, '1', 'signed smooth-axis candidate')
   end if
   call result%write('triangles', real(state%triangles, dp), '1', &
                     'one-based triangle node indices')
   call result%write('pprime', pprime, 'Pa/(Wb/rad)', 'constant pressure derivative')
   if(.not.has_profiles)call result%write('ffprime', ffprime, 'T', 'constant F*dF/dpsi')
   if (has_profiles) then
      call result%write('F_edge', fpol(size(fpol)), 'T m', &
                        'toroidal flux function at edge')
   else
      call result%write('F_edge', r0*b0, 'T m', 'toroidal flux function at edge')
   end if
   if (has_profiles) then
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
   if(prescribed_q)then
      call result%write('q_target',q_target,'1','canonical signed prescribed safety factor')
      call result%write('Phi_target',toroidal_flux_target,'Wb','signed full physical toroidal flux')
      call result%write('Phi_native',inverse%phi_native,'Wb','toroidal flux before positive homogeneity')
      call result%write('positive_normalization',inverse%positive_normalization,'1','zero-pressure homogeneity')
      call result%write('q_knot_error',inverse%knot_indicator,'1','native q sampled-knot relative error')
      call result%write('q_half_knot_error',inverse%half_indicator,'1','native q sampled half-knot relative error')
      call result%write('q_validation_error',inverse%validation_indicator,'1','native q sampled uniform rho_pol error')
      call result%write('coarea_angular_indicator',inverse%angular_indicator,'1','arc Gauss order N/16 versus2N/16 indicator')
      call result%write('Phi_quadrature_indicator',inverse%flux_indicator,'1','degree8 versus12 area indicator')
   end if
    call result%write('residual', state%residual, '1', 'relative weak equation residual')
   call result%close()
   print *, 'GS candidate nodes/triangles/residual:', size(state%psi), &
      size(state%triangles, 2), state%residual
   if(prescribed_q)then
      print *, 'prescribed q/Phi discrete indicators:',inverse%knot_indicator, &
          inverse%profile_indicator,inverse%angular_indicator,inverse%flux_indicator
   else
      print *, 'q_seed is a source parameter, not a prescribed q profile.'
   end if
contains
   real(dp) function constant_source_f(psi) result(value)
      real(dp), intent(in) :: psi
      real(dp) :: square
      square = (r0*b0)**2+2*ffprime*psi
      if (.not. ieee_is_finite(square)) error stop 'GS nonfinite F square'
      if (square <= 0) error stop 'GS toroidal flux function square is nonpositive'
      value = sign(sqrt(square), b0)
   end function constant_source_f
   real(dp) function reconstructed_f(psi) result(value)
      real(dp), intent(in) :: psi
      real(dp) :: square
      square = fpol(profile_count)**2 + 2*gs_profile_integral(flux, ff, psi)
      if (square <= 0 .or. fpol(profile_count) == 0) &
         error stop 'GS reconstructed F branch is undefined'
      value = sign(sqrt(square), fpol(profile_count))
   end function reconstructed_f
end program grad_shafranov_demo
