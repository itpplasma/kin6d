program test_gs_symbolic
   use fortsym_arena, only: arena_t
   use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
                           operator(*), operator(/), operator(**)
   use fortsym_parse, only: parse_expr
   use fortsym_diff, only: diff
   use fortsym_engine_native, only: native_engine_t, make_native_engine
   use fortsym_check, only: suite_t, suite_begin, check_identity, suite_end
   implicit none
   type(arena_t), target :: arena
   type(native_engine_t) :: eng
   type(suite_t) :: suite
   type(expr_t) :: r, z, p, bz, bt, k, br, bs, alpha
   call arena%init()
   eng = make_native_engine(arena)
   r = sym(arena, 'R'); z = sym(arena, 'Z'); k = sym(arena, 'k')
   call suite_begin(suite, 'axisymmetric equilibrium exact identities')
   p = parsed('C*((R**2-R0**2)**2+4*R**2*Z**2)')
   call check_identity(suite, eng, 'Solovev Delta_star', delta(p) - parsed('16*C*R**2'))
   br = -diff(p, z)/r; bs = diff(p, r)/r
   call check_identity(suite, eng, 'axisymmetric divergence', diff(r*br, r)/r + diff(bs, z))
   call check_identity(suite, eng, 'toroidal curl equals -Delta_star/R', &
                       diff(br, z) - diff(bs, r) + delta(p)/r)
   p = parsed('D*(R**2*log(R/R0)-(R**2-R0**2)/2+Z**2)')
   call check_identity(suite, eng, 'zero-beta logarithmic Delta_star', &
                       delta(p) - parsed('4*D'))
   bz = parsed('B0/(1+k**2*R**2)'); bt = k*r*bz
   alpha = parsed('2*k/(1+k**2*R**2)')
   call check_identity(suite, eng, 'Gold-Hoyle curl theta', -diff(bz, r) - alpha*bt)
   call check_identity(suite, eng, 'Gold-Hoyle curl z', diff(r*bt, r)/r - alpha*bz)
   call check_identity(suite, eng, 'Gold-Hoyle force balance', &
                       -diff(bz, r)*bz - diff(r*bt, r)/r*bt)
   call check_identity(suite, eng, 'Gold-Hoyle constant q', &
                       r*bz - parsed('R0*(1/(k*R0))')*bt)
   call suite_end(suite, 'exact-identities.json')
contains
   function parsed(text) result(e)
      character(*), intent(in)::text
      type(expr_t)::e
      character(:), allocatable::message
      logical::good
      e = parse_expr(arena, text, good, message)
      if (.not. good) then
         print *, message
         error stop 1
      end if
   end function
   function delta(e) result(v)
      type(expr_t), intent(in)::e
      type(expr_t)::v
      v = diff(diff(e, r), r) - diff(e, r)/r + diff(diff(e, z), z)
   end function
end program
