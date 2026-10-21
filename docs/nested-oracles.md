# Exact nested three-dimensional equilibrium references

[V] `kin6d_nested_oracles` evaluates both families in
[Landreman, arXiv:2609.26742v1](https://arxiv.org/html/2609.26742v1) and their
extensions in [Issan et al., arXiv:2610.07304v1](https://arxiv.org/html/2610.07304v1).
The common parameter `tau=0` selects Landreman; nonzero `tau` selects Issan.
The retained source is `gen/gen_nested_oracles.f90`, with checked-in output in
`src/generated/kin6d_nested_oracle_leaves.f90`. All vectors and spatial
derivatives are Cartesian. Units have mu0=1. These are ideal-MHD references,
not kinetic equilibria or a new equilibrium solver.

[V] The four inputs in `examples/nested-oracles/` select genuinely 3D cases:
eps=0.5 for the integer-transform family and eps=1.2 for the sheared family.
Both have two field periods. The Issan inputs use tau=0.4. The eps=0
axisymmetric limit remains available for reduction checks.

## Domain and mathematical definition

[V] For family 1, let a=sqrt(1+eps), b=sqrt(1-eps),
Omega²=1-tau², s=x²/a²+y²/b² and
F²=s(2-s)-4z²-4 tau s z. The domain requires 0<=eps<1, |tau|<1, F>0,
and 0<delta<(1-Omega² eps)²/4. The fields are

\[
 B=\left(\frac{(2z+\tau s)x-(a/b)Fy}{s},
         \frac{(2z+\tau s)y+(b/a)Fx}{s},1-s-2\tau z\right),
\]
\[
 Q=\frac{\Omega^2}{2}(x^2+y^2+4z^2)+2\tau z,\qquad
 \psi=\frac{\Omega^2}{2}\left(Q+\frac{|B|^2}{2}-1
       +\frac{\tau^2}{2\Omega^2}+\frac{\Omega^2\epsilon^2}{2}\right),
 \qquad p=p_a-2\psi/\Omega^2.
\]

[V] The generator derives J=curl B, grad p and grad B. Its surface map is
Issan (10), with u=-Omega² eps/2+rho cos(theta), v=rho sin(theta).
Thus psi=rho², and the positive volume density is
rho ab/Omega³. The supplied pressure is p_a=p_boundary+2 delta/Omega².
The exact averages are

\[
 V=2\pi^2ab\delta/\Omega^3,\quad
 \langle p\rangle=p_a-\delta/\Omega^2,\quad
 \langle B^2\rangle=1-\Omega^2\epsilon^2/2+\delta/\Omega^2,
 \quad \beta=2\langle p\rangle/\langle B^2\rangle,\quad\iota=2.
\]

[V] Family 2 requires eps>=0, lambda>0, delta<1/2 and
S>asin(sqrt(2 delta)). Its physical domain is the connected torus around
sigma=S, |lambda z|<pi/2. It excludes the focal segment x=0,
|y|<=sqrt(eps). Define omega=x+i y,
K=conj(omega) sqrt(1+eps/conj(omega)²), with the square root having positive
real part, and Xi=omega K+pi/2-S+i tau. Then

\[
 W=\frac{i\sin\Xi}{2K},\quad B_x+iB_y=e^{-i\lambda z}W,\quad
 B_z=\lambda^{-1}\Re(e^{-i\lambda z}\cos\Xi),\quad
 \psi=\tfrac12(\sin^2(\lambda z)+\lambda^2B_z^2),\quad
 p=p_a-\psi/\lambda^2.
\]

[V] The emitted real implementation evaluates the same branch continuously
around the torus, including x=0 outside the focal segment. It uses a real
square root of 1+eps/conj(omega)² before multiplying by conj(omega), avoiding
an inconsistent principal square root of K². The independent test evaluates
the complex formula directly.

[V] The inverse elliptic map uses X=-rho cos(theta), Y=rho sin(theta),
nu=eps sin(2t)/2+tau, and

\[
 \sigma=S+\arctan\frac{\sinh\nu\,Y}{\cosh\nu\sqrt{1-Y^2}}
            -\arcsin\frac{X}{\sqrt{\cosh^2\nu-Y^2}},\qquad
 h=\sqrt{4\sigma^2+\epsilon^2},\quad G=(h+\epsilon\cos2t)/2.
\]

[V] The position is (sqrt((h-eps)/2) cos t,
sqrt((h+eps)/2) sin t, -asin(Y)/lambda), psi=rho²/2, and the positive
volume density is rho G/[lambda h sqrt((cosh²nu-rho²)(1-Y²))].
Volume averages use radial Gauss and periodic angular quadrature of generated
fields and this density. They are E1 evaluations, with doubled-resolution
differences reported separately from any rigorous error bound.

## Exact symbolic proof and independent checks

[V] The family-1 proof is factored through the cylindrical field before its
constant horizontal stretch. FortSym proves the cylindrical divergence and
all three tension components exactly. They give
(B0·grad)B0=(-Omega² R,0,-4 Omega² z-2 tau). The constant stretch preserves
divergence and commutes with the horizontal part of this affine vector field.
Consequently (B·grad)B=-grad Q. Since p is a constant minus Q+|B|²/2,
J cross B=grad p follows from the vector identity for magnetic tension.

[V] For family 2, use the Wirtinger derivatives K_omega=0 and
K_conj(omega)=conj(omega)/K. FortSym checks the following four identities,
with the integration constant allowed to be complex:

\[
 2W_\omega=i\cos\Xi,\qquad W_{\omega\omega}=-K^2W,\qquad
 W_{\omega\bar\omega}=-|\omega|^2W,\qquad
 WW_{\omega\omega}-(W_\omega)^2=1/4.
\]

[V] The first identity cancels the horizontal divergence against partial_z Bz.
Using the remaining identities in B·grad gives horizontal tension
(grad_x+i grad_y)(|W|²/2), and vertical tension
-sin(lambda z) cos(lambda z)/lambda. Thus tension=-grad Q for
Q=-|W|²/2+sin²(lambda z)/(2 lambda²), and
Q+|B|²/2=psi/lambda² establishes force balance. The proof requires no
stellarator symmetry and therefore admits the Issan complex phase.

[V] These are exact, factored symbolic proofs; the generator refuses UNKNOWN
and accepts no randomized probe. They are not a claim that FortSym reduced
the entire expanded Cartesian residual. The retained independent test checks
the emitted real fields against complex formulas, checks their currents and
pressure gradients with five-point spatial differences, checks the map
Jacobian with spatial differences, and integrates volumes by a distinct
height-integral formula. Tests include both signs of coordinate quadrants,
both families, symmetry and broken symmetry. Generated parameter JVP/VJP
pairings, independent directional differences and interval containment are
also checked.

## Error and derivative contracts

[V] Floating field values are accompanied by strict outward interval leaves
for psi, p, B, J, grad p and grad B. `point_enclosure` encloses evaluation at
the exact stored Cartesian coordinates and parameters, including outward
formation of the axis pressure. It does not enclose the surface map's rounded
coordinates, an approximate Fourier boundary, a numerical equilibrium, or
parameter uncertainty. The selected component units differ; `field_rows` in
the result identifies every row.

[V] Smooth field and surface leaves expose generated JVP/VJP. The numerical
rotation-number driver is an oracle diagnostic with unavailable validated
ODE/sensitivity bounds. A general continuous equilibrium existence or error
certificate is unavailable here; no numerical equilibrium is promoted by an
exact-reference comparison alone.

[V] The off-axis transform of family 2 is a rotation number. Its generated
ODE is Landreman (54)/Issan (21), integrated by FortNum. The driver reports
axis/edge values and the change when doubling the transit count and tightening
ODE tolerance. Off-axis iota is E1, not a closed-form E0 scalar. Differences
are indicators, not certificates, and finite-transit bias remains separate
from local integration error. The independent on-axis check uses the
published 2.28690 value at eps=2, S=1.

## Boundary interchange and execution

[V] Build `nested_oracle_export`, then run, for example:

```sh
mkdir -p /tmp/issan-sheared
cd /tmp/issan-sheared
OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 \
  /path/to/build/nested_oracle_export \
  --config /path/to/examples/nested-oracles/issan-sheared.nml
```

[V] Results retain literal/resolved namelist, revisions, compiler information,
field evaluation intervals, volume/iota refinement indicators, boundary
coefficients and a sampled boundary error in `output/nested_oracle_export.nc`.
`output/input.nested_oracle` is a geometry-only VMEC namelist with
RBC/RBS/ZBC/ZBS(n,m), NFP=2, and LASYM selected by tau. The convention is
cos(m theta-n NFP phi), positive geometric Cartesian phi and clockwise
poloidal theta. The material chart parameter t is converted to geometric phi
by a bracketed solve before projection. Truncation is checked on a denser,
staggered grid against a requested 1% of the declared minor scale; this is
explicitly a sampled indicator.

[V] The namelist is the boundary interchange for VMEC++ and DESC's VMEC-input
reader; GVEC provides `pygvec convert-params --vmec input.nested_oracle`.
The latter converter also handles its differing angle conventions. Preserve
all four coefficient sets for tau!=0. The files supply geometry: a matched
external equilibrium run must additionally supply the toroidal flux and
pressure/iota profiles in that solver's flux coordinate. The exporter does
not silently substitute dummy profiles. Time-to-accuracy comparisons require
those matched constraints and a solved external state.

[V] `equilibrium_nested_boundary_input_1/2` reconstruct the exported namelists
with a separate Fortran consumer and compare physical boundary positions.
`equilibrium_nested_oracle_batch` checks thread invariance. The same executable
with a larger sample count reports wall time for 1/2/4/8 threads, without
allocation inside the point loop. `equilibrium_nested_oracle_regeneration`
reruns the exact proofs and compares regenerated leaves.

[V] The current map solver freezes an axisymmetric two-coordinate reference
(`ideal_frozen_t%evaluate(point(2))`, `src/mhd/ideal_map.f90`). This oracle
addition leaves that solver contract intact. Existing supported fixed-boundary
map cases remain the solver checks; a claim of a solved eps>0 equilibrium
requires a compatible reference/chart adapter and the D3/D4 error report.

## Admission of solved-case regression tests

[V] The paper-to-generator check uses Issan (4)--(7) for the integer family
and (15)--(21) for the sheared family. In the integer field, the generated
factors `(1+eps)/sqrt(1-eps**2)` and `(1-eps)/sqrt(1-eps**2)` equal a/b and
b/a on the admitted domain. In the sheared field, replacing the symbolic
phase by Re(omega K)-S absorbs the pi/2 in Xi: sin(Xi)=cos(phase+i nu) and
cos(Xi)=-sin(phase+i nu). Expanding these two identities gives the retained
real field. Setting tau=0 gives the corresponding Landreman definitions.
The surface inversion solves Issan (19) on the principal branch, and the
generated `rhs=2 G sqrt((1-Y**2)/(cosh(nu)**2-rho**2))` is (21), since
`2 G=h+eps cos(2t)` and `rho**2=2 psi`. No equation discrepancy was found in
these implemented definitions. The eight exact FortSym lemmas and the
independent complex-formula test check different parts of this argument.

[V] Exact-oracle tests are not fixed-boundary equilibrium solves. Three
current interfaces prevent registering the requested non-axisymmetric
Landreman/Issan solved-case gates as already supported:

- `ideal_frozen_t%evaluate` accepts two coordinates, with a reference field
  and pressure independent of the toroidal coordinate. Neither a 3D oracle
  initializer nor its frozen-form pullback is supplied here.
- `map_equilibrium` in `src/mhd/map_equilibrium.f90` returns
  `EQ_UNAVAILABLE` for `ntor/=0` or a non-cylindrical chart.
- `ideal_map_residual` in `src/mhd/ideal_map.f90` rejects a non-cylindrical
  chart; the D7 labelled-field contract and its surface readback are
  axisymmetric (`src/equilibrium/equilibrium_state.f90`).

[H] A solved-case gate needs an admitted reference/chart pullback, an
independent 3D field/surface/iota readback and a compatible force diagnostic,
followed by a small native convergence ladder. Extending those scientific
contracts is solver work. A boundary export, oracle interpolation, or an
axisymmetric reduction must not be counted as that solve. The sheared family
also needs pressure and iota tabulated in the receiving solver's toroidal
flux coordinate; the existing geometry-only export does not supply them.


[V] Optional installed-reader verification is retained in
`test/external/check_nested_boundary_consumers.py`. Run it with `vmecpp`,
`desc` or `gvec` and the exported input path, using that code's environment.
It records installed versions, checks nonsymmetric coefficient preservation,
and for DESC compares physical positions against the Fortran-produced NetCDF
points. These optional external readers are not KIN6D runtime dependencies.

Chris&AI
