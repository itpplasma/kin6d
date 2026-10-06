! Single symbolic source of the microscopic Coulomb physics and its scales.
!
! Conventions (all consumers):
! - g = q_i q_j/(4 pi eps0) is the pair coupling in declared reduced units;
!   g > 0 is repulsive. U(r) = g/r. The force on the particle at relative
!   position x (from the partner/source at the origin) is F = g x/|x|^3.
! - Uniform sphere of total source charge Q, radius R: g = q Q/(4 pi eps0).
! - Periodic cube of volume V with a neutralizing uniform background; Ewald
!   splitting parameter alpha (1/length); wave vectors k = 2 pi n/L, n /= 0.
!   The pair potential is phi_E = real(x) + sum_k recip(k) cos(k.x)
!   + background; self = lim_{r->0}(phi_E - g/r) is the image/self term.
! - Energies use the same reduced units as g/length; hbar is reduced too.
! Every function returns FortSym expressions in the arena of its arguments.
module kin6d_symbolic_physics
    use, intrinsic :: iso_fortran_env, only: int64
    use fortsym_expr, only: expr_t, pi_expr, rat, sqrt, exp, erf, erfc, cos, sin, &
                            operator(+), operator(-), operator(*), operator(/), &
                            operator(**)
    implicit none
    private
    public :: coulomb_potential, coulomb_force
    public :: sphere_potential_inside, sphere_force_inside
    public :: gaussian_coulomb_potential, gaussian_coulomb_series
    public :: ewald_real_potential, ewald_real_force, ewald_reciprocal_weight
    public :: ewald_self_potential, ewald_background_potential
    public :: jellium_coefficient, partial_wave_potential
    public :: yukawa_potential, centrifugal_potential, screened_partial_wave_potential
    public :: reduced_mass, bohr_radius, rydberg_energy, sommerfeld_parameter
    public :: wigner_seitz_radius, coupling_parameter, fermi_energy
    public :: degeneracy_parameter, density_parameter, debye_length
    public :: plasma_frequency, modal_poisson
    public :: phase_rotation, free_streaming_rotation, krook_node
    public :: krook_node_residual, p1_weight, p1_blend, p1_interpolate
    public :: yukawa_zero_mode, debye_structure_factor, bohm_gross_frequency
contains
    !> Bare Coulomb pair potential g/r.
    function coulomb_potential(g, r) result(u)
        type(expr_t), intent(in) :: g, r
        type(expr_t) :: u
        u = g/r
    end function

    !> Force g x/|x|^3 = -grad(g/|x|) on the particle at x.
    function coulomb_force(g, x) result(f)
        type(expr_t), intent(in) :: g, x(3)
        type(expr_t) :: f(3), r2
        integer :: j
        r2 = x(1)**2 + x(2)**2 + x(3)**2
        do j = 1, 3
            f(j) = g*x(j)/(r2*sqrt(r2))
        end do
    end function

    !> Potential energy inside a uniform sphere (radius R), as a function of
    !> the squared distance r2 from its centre: g (3 R^2 - r2)/(2 R^3).
    !> Outside the sphere the potential is coulomb_potential(g, r).
    function sphere_potential_inside(g, radius, r2) result(u)
        type(expr_t), intent(in) :: g, radius, r2
        type(expr_t) :: u
        u = g*(3*radius**2 - r2)/(2*radius**3)
    end function

    !> Force g x/R^3 inside a uniform sphere (Gauss law).
    function sphere_force_inside(g, radius, x) result(f)
        type(expr_t), intent(in) :: g, radius, x(3)
        type(expr_t) :: f(3)
        integer :: j
        do j = 1, 3
            f(j) = g*x(j)/radius**3
        end do
    end function

    !> Bare Coulomb g/r averaged over two isotropic Gaussian densities whose
    !> separation has variance s per axis and mean distance q:
    !> g erf(q/sqrt(2 s))/q (exact; finite at q = 0).
    function gaussian_coulomb_potential(g, s, q) result(u)
        type(expr_t), intent(in) :: g, s, q
        type(expr_t) :: u
        u = g*erf(q/sqrt(2*s))/q
    end function

    !> The same in T = q^2/(2 s): g sqrt(2/(pi s)) sum_j (-T)^j/(j! (2j+1)),
    !> j = 0..order (removable q = 0 form; truncation < T^(order+1)/(order+1)!).
    function gaussian_coulomb_series(g, s, t, order) result(u)
        type(expr_t), intent(in) :: g, s, t
        integer, intent(in) :: order
        type(expr_t) :: u, term
        integer :: j
        integer(int64) :: fact
        fact = 1
        u = rat(g%a, 1_int64, 1_int64)
        do j = 1, order
            fact = fact*j
            term = rat(g%a, (-1_int64)**j, fact*(2*j + 1))*t**j
            u = u + term
        end do
        u = g*sqrt(2/(pi_expr(g%a)*s))*u
    end function

    !> Ewald real-space (short-range) pair potential g erfc(alpha r)/r.
    function ewald_real_potential(g, alpha, r) result(u)
        type(expr_t), intent(in) :: g, alpha, r
        type(expr_t) :: u
        u = g*erfc(alpha*r)/r
    end function

    !> Ewald real-space force -grad(g erfc(alpha r)/r).
    function ewald_real_force(g, alpha, x) result(f)
        type(expr_t), intent(in) :: g, alpha, x(3)
        type(expr_t) :: f(3), r, a
        integer :: j
        r = sqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))
        a = g*(erfc(alpha*r)/r**3 + &
               2*alpha*exp(-alpha*alpha*r*r)/(sqrt(pi_expr(g%a))*r*r))
        do j = 1, 3
            f(j) = a*x(j)
        end do
    end function

    !> Bare Coulomb Fourier coefficient on a periodic volume:
    !> g/r = sum_k [4 pi g/(V k^2)] exp(i k.x) (k=0 removed by neutrality).
    function jellium_coefficient(g, volume, k2) result(c)
        type(expr_t), intent(in) :: g, volume, k2
        type(expr_t) :: c
        c = 4*pi_expr(g%a)*g/(volume*k2)
    end function

    !> Ewald reciprocal coefficient of cos(k.x), summed over all k /= 0.
    function ewald_reciprocal_weight(g, alpha, volume, k2) result(w)
        type(expr_t), intent(in) :: g, alpha, volume, k2
        type(expr_t) :: w
        w = jellium_coefficient(g, volume, k2)*exp(-k2/(4*alpha*alpha))
    end function

    !> Self term of the Gaussian screening charge: -g alpha/sqrt(pi).
    function ewald_self_potential(g, alpha) result(u)
        type(expr_t), intent(in) :: g, alpha
        type(expr_t) :: u
        u = -g*alpha/sqrt(pi_expr(g%a))
    end function

    !> Neutralizing-background (k=0) term of the pair potential:
    !> -pi g/(V alpha^2).
    function ewald_background_potential(g, alpha, volume) result(u)
        type(expr_t), intent(in) :: g, alpha, volume
        type(expr_t) :: u
        u = -pi_expr(g%a)*g/(volume*alpha*alpha)
    end function

    !> Radial effective potential of partial wave l in the relative
    !> coordinate: g/r + hbar^2 l(l+1)/(2 mu r^2).
    function partial_wave_potential(g, hbar, mu, l, r) result(v)
        type(expr_t), intent(in) :: g, hbar, mu, l, r
        type(expr_t) :: v
        v = coulomb_potential(g, r) + centrifugal_potential(hbar, mu, l, r)
    end function

    !> Centrifugal term hbar^2 l(l+1)/(2 mu r^2) of partial wave l.
    function centrifugal_potential(hbar, mu, l, r) result(v)
        type(expr_t), intent(in) :: hbar, mu, l, r
        type(expr_t) :: v
        v = hbar**2*l*(l + 1)/(2*mu*r**2)
    end function

    !> Statically screened (Debye/Yukawa) pair potential g exp(-r/lambda)/r:
    !> the bare Coulomb potential times the screening factor; lambda -> oo
    !> recovers coulomb_potential.
    function yukawa_potential(g, lambda, r) result(u)
        type(expr_t), intent(in) :: g, lambda, r
        type(expr_t) :: u
        u = coulomb_potential(g, r)*exp(-r/lambda)
    end function

    !> Radial effective potential of partial wave l for the screened pair:
    !> g exp(-r/lambda)/r + hbar^2 l(l+1)/(2 mu r^2).
    function screened_partial_wave_potential(g, hbar, mu, l, r, lambda) result(v)
        type(expr_t), intent(in) :: g, hbar, mu, l, r, lambda
        type(expr_t) :: v
        v = yukawa_potential(g, lambda, r) + centrifugal_potential(hbar, mu, l, r)
    end function

    function reduced_mass(m1, m2) result(mu)
        type(expr_t), intent(in) :: m1, m2
        type(expr_t) :: mu
        mu = m1*m2/(m1 + m2)
    end function

    !> Bohr radius of the pair, hbar^2/(mu |g|).
    function bohr_radius(hbar, mu, g) result(a)
        type(expr_t), intent(in) :: hbar, mu, g
        type(expr_t) :: a
        a = hbar**2/(mu*magnitude(g))
    end function

    !> Rydberg energy mu g^2/(2 hbar^2) (binding energy of the 1s state).
    function rydberg_energy(hbar, mu, g) result(e)
        type(expr_t), intent(in) :: hbar, mu, g
        type(expr_t) :: e
        e = mu*g**2/(2*hbar**2)
    end function

    !> Sommerfeld parameter eta = g/(hbar v) at relative speed v.
    function sommerfeld_parameter(g, hbar, v) result(eta)
        type(expr_t), intent(in) :: g, hbar, v
        type(expr_t) :: eta
        eta = g/(hbar*v)
    end function

    !> Wigner-Seitz radius a = (3/(4 pi n))^(1/3) at number density n.
    function wigner_seitz_radius(n) result(a)
        type(expr_t), intent(in) :: n
        type(expr_t) :: a
        a = (3/(4*pi_expr(n%a)*n))**rat(n%a, 1_int64, 3_int64)
    end function

    !> Classical coupling Gamma = |g|/(a T), T in energy units.
    function coupling_parameter(g, a, temperature) result(gamma)
        type(expr_t), intent(in) :: g, a, temperature
        type(expr_t) :: gamma
        gamma = magnitude(g)/(a*temperature)
    end function

    !> Fermi energy of spin-1/2 fermions, hbar^2 (3 pi^2 n)^(2/3)/(2 m).
    function fermi_energy(hbar, m, n) result(e)
        type(expr_t), intent(in) :: hbar, m, n
        type(expr_t) :: e
        e = hbar**2*(3*pi_expr(n%a)**2*n)**rat(n%a, 2_int64, 3_int64)/(2*m)
    end function

    !> Degeneracy Theta = T/E_F.
    function degeneracy_parameter(temperature, fermi) result(theta)
        type(expr_t), intent(in) :: temperature, fermi
        type(expr_t) :: theta
        theta = temperature/fermi
    end function

    !> Density parameter r_s = a/a_B.
    function density_parameter(a, bohr) result(rs)
        type(expr_t), intent(in) :: a, bohr
        type(expr_t) :: rs
        rs = a/bohr
    end function

    !> Debye length sqrt(T/(4 pi n g_self)), g_self = q^2/(4 pi eps0).
    function debye_length(temperature, n, g_self) result(lambda)
        type(expr_t), intent(in) :: temperature, n, g_self
        type(expr_t) :: lambda
        lambda = sqrt(temperature/(4*pi_expr(n%a)*n*g_self))
    end function

    !> Plasma frequency sqrt(4 pi n g_self/m).
    !> |g| written as sqrt(g^2): equal for real g and its derivative g/|g|
    !> is representable in generated JVP/VJP leaves (undefined at g=0).
    function magnitude(g) result(a)
        type(expr_t), intent(in) :: g
        type(expr_t) :: a
        a = sqrt(g**2)
    end function

    function plasma_frequency(n, g_self, m) result(omega)
        type(expr_t), intent(in) :: n, g_self, m
        type(expr_t) :: omega
        omega = sqrt(4*pi_expr(n%a)*n*g_self/m)
    end function

    !> Periodic 1D modal Poisson law -eps phi'' = rho, E = -phi', for one
    !> mode k /= 0 in the real basis rho = rc cos(kx) + rs sin(kx).
    !> Returns [phi_c, phi_s, E_c, E_s] in the same basis. A complex
    !> coefficient rho_k of exp(ikx) maps to rc = Re rho_k, rs = -Im rho_k
    !> (then phi_k = phi_c - i phi_s, E_k = E_c - i E_s).
    function modal_poisson(eps, k, rc, rs) result(mode)
        type(expr_t), intent(in) :: eps, k, rc, rs
        type(expr_t) :: mode(4)
        mode = [rc/(eps*k**2), rs/(eps*k**2), -rs/(eps*k), rc/(eps*k)]
    end function

    !> Rotation of the Fourier pair (a, b) of a cos + b sin by the angle with
    !> cosine c and sine s: [a c - b s, a s + b c].
    function phase_rotation(a, b, c, s) result(rotated)
        type(expr_t), intent(in) :: a, b, c, s
        type(expr_t) :: rotated(2)
        rotated = [a*c - b*s, a*s + b*c]
    end function

    !> Exact free streaming f_t + v f_x = 0 of one mode a cos(kx) + b sin(kx)
    !> over time t: f(x - v t) is the rotation by phase = k v t.
    function free_streaming_rotation(a, b, phase) result(rotated)
        type(expr_t), intent(in) :: a, b, phase
        type(expr_t) :: rotated(2)
        rotated = phase_rotation(a, b, cos(phase), sin(phase))
    end function

    !> Stationary Krook node v f_x + nu f = s for one mode with kv = k v:
    !> [nu, kv; -kv, nu][a; b] = [sc; ss], determinant nu^2 + kv^2 > 0.
    function krook_node(nu, kv, sc, ss) result(node)
        type(expr_t), intent(in) :: nu, kv, sc, ss
        type(expr_t) :: node(2)
        node = [(nu*sc - kv*ss)/(nu**2 + kv**2), (kv*sc + nu*ss)/(nu**2 + kv**2)]
    end function

    !> Residual of the krook_node equations at the candidate (a, b).
    function krook_node_residual(nu, kv, a, b, sc, ss) result(residual)
        type(expr_t), intent(in) :: nu, kv, a, b, sc, ss
        type(expr_t) :: residual(2)
        residual = [nu*a + kv*b - sc, nu*b - kv*a - ss]
    end function

    !> Affine (P1) cell weight of v on [vl, vr].
    function p1_weight(v, vl, vr) result(w)
        type(expr_t), intent(in) :: v, vl, vr
        type(expr_t) :: w
        w = (v - vl)/(vr - vl)
    end function

    !> P1 blend of the endpoint values al, ar at cell weight w.
    function p1_blend(al, ar, w) result(value)
        type(expr_t), intent(in) :: al, ar, w
        type(expr_t) :: value
        value = al + (ar - al)*w
    end function

    !> P1 interpolant of nodal values al, ar on [vl, vr] at v.
    function p1_interpolate(v, vl, vr, al, ar) result(value)
        type(expr_t), intent(in) :: v, vl, vr, al, ar
        type(expr_t) :: value
        value = p1_blend(al, ar, p1_weight(v, vl, vr))
    end function

    !> k=0 Fourier coefficient 4 pi g/(V kappa^2) of the periodic image sum
    !> of yukawa_potential (lambda = 1/kappa). The neutral periodic screened-Poisson Green
    !> function is the image sum minus this constant (k /= 0 modes only).
    function yukawa_zero_mode(g, kappa, volume) result(c)
        type(expr_t), intent(in) :: g, kappa, volume
        type(expr_t) :: c
        c = 4*pi_expr(g%a)*g/(volume*kappa**2)
    end function

    !> Debye-Hueckel (weak-coupling OCP) static structure factor
    !> S(k) = k^2/(k^2 + kappa^2), kappa = 1/lambda_D of the responsive species.
    function debye_structure_factor(k2, kappa2) result(s)
        type(expr_t), intent(in) :: k2, kappa2
        type(expr_t) :: s
        s = k2/(k2 + kappa2)
    end function

    !> Warm-fluid (Bohm-Gross) Langmuir frequency omega_p sqrt(1 + 3 k^2 lambda^2).
    function bohm_gross_frequency(omega_p, k, lambda) result(omega)
        type(expr_t), intent(in) :: omega_p, k, lambda
        type(expr_t) :: omega
        omega = omega_p*sqrt(1 + 3*k**2*lambda**2)
    end function
end module kin6d_symbolic_physics
