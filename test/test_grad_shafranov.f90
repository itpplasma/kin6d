program test_grad_shafranov
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
                   gs_linear_response, sample_gs, solve_gs_profiles, gs_profile_integral
   use kin6d_grad_shafranov_leaves, only: gs_solovev, gs_field_jvp, gs_field_vjp, &
      gs_gold_hoyle, gs_gold_hoyle_jvp
   implicit none
   type(gs_state_t) :: coarse, fine
   real(dp) :: boundary(2, 4), ec, ef, exact, pr, pz, br, bz, bt
   real(dp) :: lhs, rhs, dr, dpr, dpz, df, ar, apr, apz, af
   real(dp) :: radius, k, cylinder_bz, cylinder_bt, cylinder_psi
   real(dp), allocatable :: tangent(:), adjoint(:), source(:), goal(:)
   real(dp), parameter :: c = 0.002_dp, r0 = 6.2_dp
   integer :: status, i
   k = 1/(1.5_dp*r0)
   do i = 1, 20
      radius = 0.62_dp*i/20
      call gs_gold_hoyle(radius, r0, 5.3_dp, k, &
         cylinder_bz, cylinder_bt, cylinder_psi)
      call gs_gold_hoyle_jvp(radius, r0, 5.3_dp, k, &
         1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, dr, dpr, dpz)
      ! Cylindrical curl crossed with B vanishes; q is independently prescribed.
      if (abs(-dr*cylinder_bz - (cylinder_bt/radius + dpr)* &
         cylinder_bt) > 1.0e-12_dp) error stop 'Gold-Hoyle physical force balance'
      if (abs(radius*cylinder_bz/(r0*cylinder_bt) - 1.5_dp) > &
         1.0e-12_dp) error stop 'Gold-Hoyle independent constant q'
      if (abs(dpz - r0*cylinder_bt) > 1.0e-12_dp) &
         error stop 'Gold-Hoyle cylindrical flux derivative'
   end do
   ! Integral of the independent linear derivative 2*psi+3, edge=2.
   if (abs(gs_profile_integral([0.0_dp, 1.0_dp, 2.0_dp], &
                            [3.0_dp, 5.0_dp, 7.0_dp], 0.5_dp) + 8.25_dp) > 1.0e-13_dp) &
      error stop 'GS derivative-controlled profile integral'
   if (abs(gs_profile_integral([0.0_dp, 1.0_dp, 2.0_dp], &
      [3.0_dp, 5.0_dp, 7.0_dp], -1.0_dp) + 13.0_dp) > 1.0e-13_dp) &
      error stop 'GS derivative extension below axis'
   if (abs(gs_profile_integral([0.0_dp, 1.0_dp, 2.0_dp], &
      [3.0_dp, 5.0_dp, 7.0_dp], 3.0_dp) - 7.0_dp) > 1.0e-13_dp) &
      error stop 'GS derivative extension beyond edge'
   boundary = reshape([5.8_dp, -0.4_dp, 6.6_dp, -0.4_dp, &
                       6.6_dp, 0.4_dp, 5.8_dp, 0.4_dp], [2, 4])
   call oracle_case(coarse, 0.0032_dp, ec)
   call oracle_case(fine, 0.0008_dp, ef)
   if (ef >= 0.65_dp*ec) error stop 'GS independent oracle refinement'
   if (ef > 0.004_dp) error stop 'GS independent oracle error'
   call sample_gs(fine, 6.25_dp, 0.075_dp, 32.86_dp, exact, br, bz, bt, status)
   if (status /= 0) error stop 'GS physical sample unavailable'
   if (abs(bt - 32.86_dp/6.25_dp) > 1.0e-12_dp) error stop 'GS field units'
   allocate (source(size(fine%psi)), goal(size(fine%psi)))
   do i = 1, size(source)
      source(i) = sin(real(i, dp))
      goal(i) = cos(real(i, dp))
   end do
   where (fine%boundary) source = 0.0_dp
   where (fine%boundary) goal = 0.0_dp
   call gs_linear_response(fine, source, tangent, status)
   if (status /= 0) error stop 'GS equation tangent failed'
   call gs_linear_response(fine, goal, adjoint, status)
   if (status /= 0) error stop 'GS equation adjoint failed'
   lhs = dot_product(goal, tangent); rhs = dot_product(source, adjoint)
   if (abs(lhs - rhs) > 1.0e-8_dp*max(abs(lhs), abs(rhs))) &
      error stop 'GS implicit tangent/adjoint pairing'
   call gs_field_jvp(6.0_dp, 2.0_dp, 3.0_dp, 30.0_dp, &
                     0.5_dp, 0.25_dp, -0.75_dp, 2.0_dp, dr, dpr, dpz)
   call gs_field_vjp(6.0_dp, 2.0_dp, 3.0_dp, 30.0_dp, &
                     2.0_dp, -1.0_dp, 0.5_dp, ar, apr, apz, af)
   lhs = 2*dr - dpr + 0.5_dp*dpz
   rhs = 0.5_dp*ar + 0.25_dp*apr - 0.75_dp*apz + 2*af
   if (abs(lhs - rhs) > 1.0e-13_dp) error stop 'GS generated field JVP/VJP'
   ! Independent flux-dependent oracle: psi=R^2 has Delta*psi=0.
   ! pprime=-a/mu0 and FFprime=a*psi cancel exactly in the named PDE.
   ! Initial interior values are zero; only the exact boundary is supplied.
   fine%psi = 0.0_dp
   where (fine%boundary) fine%psi = fine%points(1, :)**2
   call solve_gs_profiles(fine, 1.0_dp, [0.0_dp, 50.0_dp], &
                          [-0.01_dp, -0.01_dp], [0.0_dp, 0.5_dp], status)
   if (status /= 0) error stop 'GS flux-dependent equation failed'
   if (maxval(abs(fine%psi - fine%points(1, :)**2)) > 0.002_dp) &
      error stop 'GS flux-dependent independent oracle'
   print *, 'PASS: GS analytic refinement and equation derivatives', ec, ef
contains
   subroutine oracle_case(state, area, error)
      type(gs_state_t), intent(out) :: state
      real(dp), intent(in) :: area
      real(dp), intent(out) :: error
      real(dp), parameter :: mu0 = 1.2566370614359173e-6_dp
      real(dp) :: value, qr, qz
      integer :: j, stat
      call initialize_gs(state, boundary, area, stat)
      if (stat /= 0) error stop 'GS oracle mesh rejected'
      do j = 1, size(state%psi)
         call gs_solovev(state%points(1, j), state%points(2, j), &
                         r0, c, value, qr, qz)
         if (state%boundary(j)) state%psi(j) = value
      end do
      call solve_gs(state, mu0, -16*c/mu0, 0.0_dp, stat)
      if (stat /= 0) error stop 'GS oracle solve failed'
      error = 0.0_dp
      do j = 1, size(state%psi)
         call gs_solovev(state%points(1, j), state%points(2, j), &
                         r0, c, value, qr, qz)
         error = max(error, abs(state%psi(j) - value))
      end do
   end subroutine oracle_case
end program test_grad_shafranov
