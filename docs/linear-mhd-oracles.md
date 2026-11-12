# Linear cylinder and island references

These test-side references supply candidate comparison values for the linear
wave and instability benchmark RM-B2 and the Rutherford island benchmark RM-B5
in [the nonlinear resistive MHD design](nonlinear-resistive-mhd.md).
They do not establish either gate for a nonlinear solver. Here MHD means
magnetohydrodynamics, and all dimensional interfaces use SI units. The same
library holds the toroidal references of [the torus case](vrmhd-torus.md):
the Bussac internal-kink energy and threshold, the large-aspect-ratio
Shafranov shift, the geodesic-acoustic frequency and the Suydam index, in
the nondimensional units of that page.

Build `test_rm_oracles` with the native CPU preset, then run:

```sh
ctest --test-dir build -L rm-oracle -j 8 --output-on-failure
ctest --test-dir build -R '^rm_oracle_table$' -j 8 -V
```

The second command prints the reference table. With the pinned FortSym build
configured as described in [tooling.md](tooling.md), also build
`gen_rm_oracles`; the label then includes a regeneration comparison.
The mathematical source is `gen/gen_rm_oracles.f90`, the emitted leaves are
`src/generated/milca_rm_oracle_leaves.f90`, and the shooting and root-finding
composition is `test/support/rm_mhd_oracles.f90`. This library is linked only
by the oracle test executable. FortNum supplies adaptive Cash--Karp integration
and Brent root finding; the reference has no dependency on a nonlinear MHD
residual.

## Geometry and normalization

The periodic cylinder has radius `a` and length `2 pi R0`. A perturbation has
phase `m theta - n z/R0`, where `m` and `n` are positive integers, `theta` is
the azimuthal angle and `z` is axial distance. Thus the signed axial wave
number is `kz = -n/R0`, and

```text
F(r) = m Btheta(r)/r + kz Bz(r),      F(rs) = 0,
q(r) = r Bz(r)/(R0 Btheta(r)).
```

The resonant radius is `rs`; primes below mean radial derivatives. `F` has
units T/m and `F'` has units T/m^2. The magnetic permeability is
`mu0 = 4 pi 10^-7 H/m`, mass density is `rho` in kg/m^3, electrical
resistivity is `eta` in Ohm m, and magnetic diffusivity is `D = eta/mu0` in
m^2/s. Growth rates and angular frequencies are in 1/s.

The Wesson-type table uses the reduced, zero-beta cylinder equation with
constant axial field. It is a strong-guide-field approximation. A constant
`Bz` and a nonzero axial current do not constitute an exact zero-pressure
force-balanced equilibrium of full cylindrical MHD. Full zero-pressure
Newcomb comparisons instead require a force-balanced supplied field, such
as the force-free Bessel field used in the independent test below.

For the table, `jz = j0 (1-r^2/a^2)^nu`, with `nu = q(a)/q0 - 1`, where
`q0 = q(0)`. Ampere's law gives

```text
1/q = [1-(1-r^2/a^2)^(nu+1)] / [(nu+1) q0 (r^2/a^2)],
J = (mu0 R0/Bz) jz = (2/q0) (1-r^2/a^2)^nu,
Fn = m/q - n,                      F = (Bz/R0) Fn.
```

`J` and `Fn` are dimensionless. The generated axis expansion avoids
cancellation in `1/q` for `r/a < 0.02`. The reference accepts `q(a)>2q0`,
so that `J'` is finite at the wall. It requests `J''` only inside the plasma;
this derivative diverges at the wall when `1<nu<2`.

## Waves and resistive rates

For a uniform field vector `B`, arbitrary wave vector `k`, pressure `p` and
adiabatic index `Gamma`, define `vA^2 = |B|^2/(mu0 rho)` and
`cs^2 = Gamma p/rho`. `uniform_waves` returns squared angular frequencies
in the order Alfven, fast, slow:

```text
omegaA^2 = (k.B)^2/(mu0 rho),
omegaf,s^2 = [k^2(vA^2+cs^2)
              +/- sqrt(k^4(vA^2+cs^2)^2-4 cs^2 k^2 omegaA^2)]/2.
```

The implementation uses a nonnegative equivalent discriminant involving
`|k cross B|^2`, and obtains the slow root from the product of the roots to
avoid subtractive cancellation. Its independent oracle diagonalizes the
Cartesian linearized momentum/induction force matrix with LAPACK. Parallel,
perpendicular, cold, unmagnetized, zero-wave-vector and coincident-speed
limits are checked separately. Invalid density or pressure is rejected.

The inviscid, constant-psi Furth--Killeen--Rosenbluth (FKR) tearing rate is

```text
C = 2 pi GammaEuler(3/4)/GammaEuler(1/4),
gammaFKR = C^(-4/5) Delta'^(4/5) D^(3/5) [F'(rs)^2/(mu0 rho)]^(1/5).
```

`GammaEuler` is Euler's gamma function; `C^(-4/5)` is approximately 0.55.
`Delta'`, defined below, has units 1/m and must be positive for this branch.
This is the dimensional form of [FKR, Phys. Fluids 6, 459 (1963)](
https://doi.org/10.1063/1.1706761). The coefficient is checked independently
against the small-growth limit of the complete resistive-inertial dispersion
relation of [Ara et al., Ann. Phys. 112, 443 (1978)](
https://doi.org/10.1016/S0003-4916(78)80007-4), already implemented by
`ri_dispersion_gamma`:

```text
d = [mu0 rho D^2/F'^2]^(1/6),       Q = gamma d^2/D,
Delta' d = -(pi/8) Q^(5/4)
           GammaEuler[(Q^(3/2)-1)/4]/GammaEuler[(Q^(3/2)+5)/4].
```

Here `d` is a layer length. Its `Q -> 0` limit is
`Delta' d = C Q^(5/4)`. The test uses `D=10^-10 m^2/s` and allows a relative
difference of `3e-5`, including the finite-`Q` correction, rather than
asserting equality with an asymptotic expression at finite resistivity.

The marginal-ideal resistive internal-kink reference is

```text
gammaK = [D F'(rs)^2/(mu0 rho)]^(1/3)
       = tauH^(-2/3) tauR^(-1/3),
tauH = sqrt(mu0 rho)/(rs |F'|),     tauR = rs^2/D.
```

This is the zero-ideal-energy limit of Coppi et al., Sov. J. Plasma Phys. 2,
533 (1976), and the `Q=1` pole of Ara's dispersion relation. The physical
scales and Newcomb equation are also set out in the public
[PPPL internal-kink tutorial, equations (2)--(6) and (65)--(67)](
https://w3.pppl.gov/cemm/Project/Giant_Sawteeth_Tutorial_01.pdf).
It is not a general finite-ideal-energy or viscous kink growth rate.
The table labels it as an asymptote; it must not be mistaken for the
eigenvalue of the supplied finite-aspect-ratio equilibrium.

## Outer shooting and the resonant logarithm

`shoot_outer` integrates `psi'' + P psi' - U psi = 0` separately from axis
and wall to `rs +/- cut`. The default axis condition is the regular
`psi proportional to r^m`, started at `10^-5 a`; the conducting wall has
`psi(a)=0`. Optional logarithmic derivatives allow exact current-jump and
no-wall matching in independent tests. The integrator has relative tolerance
`tolerance`, absolute tolerance `0.01 tolerance` and maximum step `0.02 a`.

For the reduced cylinder, `P=1/r` and
`U=m^2/r^2 + m J'/(r Fn)`. For full zero-pressure Newcomb, the callback
`newcomb_coefficients` transforms the radial displacement `xi_r` to
`psi = r F xi_r/m`. Define `K2=m^2/r^2+kz^2` and `v=K2'/K2`. Then

```text
P = -1/r-v,
U = K2-2/r^2-v/r + N/F,
N = F''-(v-1/r)F' + 2 kz^2 (kz Bz-m Btheta/r)/(r^2 K2).
```

The caller supplies the field and its `F` derivatives consistently. At a
simple resonance the pole residue is `h=N(rs)/F'(rs)`; the reduced equation
has `h=m J'(rs)/(rs Fn'(rs))`. At exactly `F=0` the helper returns the
regular part of `U` and the residue, not a finite value of the singular
coefficient. Only isolated simple resonances are supported.

With `x=r-rs`, write `P=P0+O(x)` and `U=h/x+U0+O(x)`. The explicit
Frobenius basis used on both sides is

```text
L = 1+h x log(|x|/a)+c x^2 log(|x|/a)+e x^2,
S = x+(h-P0)x^2/2,
c = h(h-P0)/2,                   e = (U0-P0 h-3c)/2,
psi = psi_s [L+b_minus,plus S].
```

The two coefficients are fitted from both `psi` and `psi'` at the cut;
`Delta'=b_plus-b_minus` removes the common logarithmic divergence. The
returned slopes are these finite expansion coefficients, not unregularized
logarithmic derivatives. Their sum, denoted `Sigma`, uses the same scale `a`
inside the logarithm. `P0` and `U0` are estimated by symmetric evaluation
around the resonance at a distance at most `10^-4 a`. An outer branch with
vanishing `psi_s` returns `oracle_singular`, rather than a finite large
tearing index. In particular, the reduced `1/1` internal-kink branch has no
finite `Delta'` in this formulation.

Three independent solutions check the shooting. A manufactured equation has
exact branches `1+h x log|x|` and `x`, testing the logarithmic finite part.
A full force-free field `Bz=J0(lambda r)`, `Btheta=J1(lambda r)` uses Bessel
functions `J0,J1`, `lambda=1.8 m^-1` and `kz=-1 m^-1`. With
`kappa=sqrt(lambda^2-kz^2)`, the independent flux branches are
`H_Z=lambda Z1(kappa r)+kz r kappa Z1'(kappa r)`, for `Z=J,Y`.
The inner branch is `H_J`, and the outer branch is
`H_J(r) H_Y(a)-H_Y(r) H_J(a)`, so its logarithmic derivative jump is analytic.

Finally, a top-hat axial current has radius `r1`, constant `q0` inside,
`q=q0(r/r1)^2` outside and `rs=r1 sqrt(2/q0)` for `m/n=2/1`.
Integrating the delta-function current gradient gives
`[psi']=-g psi/r1`, `g=2/(1-q0/2)`. Match the axis branch `r^2` through this
jump and the exterior powers `r^2,r^-2`. Continuation to infinity gives

```text
rs Delta'_infinity = 16(q0-1)/(q0-2)^2,
Delta'_wall = Delta'_infinity - 4 rs^3/(a^4-rs^4).
```

Thus this family has the no-wall marginal point `q0=1`, also stated for a
top-hat profile by [Aydemir et al., IAEA 1988 proceedings, vol. 2, p. 141](
https://web2.ph.utexas.edu/~morrison/89IAEA_morrison.pdf).
The test brackets that point and cross-checks the existing response-strip
normalization. Absolute tolerances `3e-5 /m` for the power and Bessel
solutions cover the finite resonance cut; they are looser than the observed
integration error. Separate ladders refine the ODE tolerance and cut, requiring
agreement within `10^-5 /m` and improved differences. These are numerical
convergence checks, not rigorous error bounds.

## Rutherford growth and saturation

`W` always denotes the full radial island width. The Rutherford reference is
`dW/dt=1.22 D Delta'(W)`; a negative index therefore gives decay. The width
convention and coefficient are stated explicitly in
[Hastie et al., IAEA 2004 PD/1-1, introduction](
https://www.osti.gov/etdeweb/servlets/purl/20644171), citing Rutherford,
Phys. Fluids 16, 1903 (1973), and White et al., Phys. Fluids 20, 800 (1977).

For a symmetric slab, the Militello--Porcelli / Escande--Ottaviani result is
`Wsat=-Delta'/(0.41 beta)`, where `beta=J''/J<0` at the resonance.
The independent Harris-sheet case has `B_y proportional to tanh(x/a)`,
`beta=-2/a^2`, and checks this width with specified positive `Delta'`.
See [Militello and Porcelli, Phys. Plasmas 11, L13 (2004)](
https://doi.org/10.1063/1.1677089) and [Escande and Ottaviani,
Phys. Lett. A 323, 278 (2004)](https://doi.org/10.1016/j.physleta.2004.02.010).

For the cylindrical example, `cylinder_saturated_width` solves the uniform
resistivity, fixed equilibrium current-source result of
[Arcis, Escande and Ottaviani, Phys. Plasmas 13, 052305 (2006),
equation (58), model A](https://arxiv.org/abs/physics/0603204):

```text
0 = 1.22 Delta' + W [alpha/2 (alpha log(W/a)+Sigma/2)
                     -2.21 alpha^2+0.40 alpha/rs+beta/2],
Jstar = J-2n/m,      alpha = J'/Jstar,      beta = J''/Jstar.
```

This source defines `Sigma` by subtracting `2 alpha [1+log(|x|/a)]`
from the sum of normalized flux derivatives, exactly the slope sum used
above. The formula is invariant under a change of normalization length.
Its symmetric limit has coefficient `1/(2*1.22)`, agreeing with the published
rounded `0.41` to relative `4e-4`; the test allows `5e-4`, including roundoff.
Other resistivity/temperature-evolution models have different saturation
laws; this reference chooses model A explicitly. It omits higher orders in
`W/a`, assumes a single small island with small `Delta' a`, and excludes
the `m=1` internal kink.

The cylinder example uses the published `nu=2` family. Its conducting-wall
`2/1` marginal point is approximately `rs/a=0.88` in Hastie et al., section 5;
the test checks opposite signs at `0.875` and `0.885`. A near-marginal case
at `rs/a=0.8788` supplies a positive, small saturated width, also checked
under uniform rescaling of every length. The requested `q0=0.9,1.1`,
`q(a)=3` profiles are printed with `a=1 m`, `Bz=1 T`, `R0=3 m`,
`rho=10^-6 kg/m^3`, and `eta=10^-8 Ohm m`. Unavailable tearing indices are
marked by status: `-2` means no resonance and `-3` means singular outer
flux; zero-valued fields on those rows are placeholders. These references
remain candidate values without certified bounds or a nonlinear validation.
