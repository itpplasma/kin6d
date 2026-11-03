# Invariant visco-resistive MHD volume assembly

The invariant section assembler evaluates the mixed-current volume residual,
its matrix-free state Jacobian action, and the state-dependent mass terms of
[the Cartesian weak form](vrmhd-weak-form.md). It accepts the rotation,
translation and screw cases of `gs_symmetry_t`. This is a candidate operator;
the stationary wrapper adds fixed boundary values and a mixed Coulomb gauge.
Continuous error certificates are not supplied.

## Representation and measure

Each harmonic block stores six scalar H1 coefficient vectors, in the order
`rho, vx, vy, vz, Te, Ti`, followed by the H(curl) vectors `A, J`. Each
H(curl) vector stores the section Nedelec coefficients followed by the scalar
H1 coefficients of `a = xi.A` or `xi.J`. Here H1 controls values and first
derivatives, H(curl) controls values and curls, and xi is the case's Killing
field. The harmonic count is explicit in the layout and assembly loops;
counts other than one are rejected until nonzero harmonics are implemented.

The section components specify the pullback of the horizontal one-form,
orthogonal to xi. On the screw chart, with `s=1+k^2(x^2+y^2)`,

```
A = (u1, u2, k(y u1 - x u2)) + a (-k y, k x, 1)/s.
```

For rotation the physical section at orbit coordinate zero is `(R,0,Z)`;
`A=(u1,a/R,u2)`. Physical vectors and all gradients supplied to the weak-form
leaves are Cartesian. Velocity components on the section are transported by
the isometry; they are not independently constant Cartesian components along
an orbit. Scalar gradients and vector gradients/curls include that transport.
FortSym differentiates the Killing-field ansatz, including screw twist. Its
identities check that only the section curl of the Nedelec field is needed,
that a scalar gradient has zero curl, and that
`curl(xi/s)=2 k xi/s^2` for the Cartesian screw chart.

FortFEM supplies the section topology, oriented Nedelec and Lagrange bases,
quadrature, Piola maps and discrete curl/gradient. The scalar and vector
fields use the same element and harmonic traversal. The measure is
`R dR dZ` for rotation and `dx dy` for the Cartesian translation/screw chart:
as in [the GS reduction](grad-shafranov.md), quotient area is
`dA_q=dx dy/|xi|`, so `|xi| dA_q=dx dy`. No additional screw norm factor is
applied. Integration is per unit orbit parameter; no orbit period is included.

## Polynomial degrees and quadrature

`vrmhd_section_build` takes the magnetic edge order `degree = k_A` and
optional `scalar_degree = k_s`, which defaults to `k_A`. Ordinary construction
admits equal degrees only. Unequal pairs require `study_pairing=.true.` and
`k_s >= k_A-1`; they remain available for comparison, not as the default.
The [pairing study](#pairing-study-for-the-ideal-uniform-equilibrium) found
`k` energy-bearing static modes per cell for equal degree `k`, fewer than for
any tested unequal pair. Static modes are permitted only with the structural
and dissipative checks below; their presence alone does not establish spectral
pollution or stability of the full system.

[FortFEM](https://github.com/lazy-fortran/fortfem)'s
`fortfem_meridian_feec_n` uses continuous Lagrange `P_k` for its scalar
space, first-kind Nedelec order `k` for section one-forms, and discontinuous
`P_{k-1}` for their section curl. The edge space contains `[P_{k-1}]^2`
and the homogeneous degree-`k` rotational enrichment; it has `k(k+2)` local
coefficients. The face space contains the rotated edge space (Raviart–Thomas
degree `k-1` in the standalone basis API) and discontinuous `P_{k-1}` for
the orbit component. In particular, FortFEM's complex order `1` means the
lowest-order Nedelec/Raviart–Thomas pair, rather than a section curl of degree
one.

One such complex cannot provide both requested orders. `op%fs` is the order
`k_A` magnetic complex for A and J, including their orbit Lagrange `P_{k_A}`
component; `op%scalar_fs` supplies continuous `P_{k_s}` for density, each
velocity component, both temperatures and the Coulomb multiplier. Both
complexes use the same mesh and quadrature rule. Basis values and gradients
of the magnetic orbit component always come from the magnetic complex.
FortFEM supports complex orders through eight, so this construction admits
`1 <= k_A = k_s <= 8` ordinarily, with the stated study exception.

For `p=max(k_s,k_A)` and geometry degree `g`, both complexes request
triangle quadrature degree `4 p + 2 g + 2`. On affine translation elements,
the largest polynomial product is the skew momentum advection term
`rho (grad v) v.w`, of degree `4 k_s-1`; the other polynomial mass and
magnetic products have no higher degree for the admitted pairs. The extra
quadrature degree oversamples geometry. Rotation/screw metric factors and
state-dependent transport closures are nonpolynomial, so no finite polynomial
rule makes those integrals exact. Energy exchange still uses identical
quadrature points and weights on both sides of each pairing. The generated
pointwise leaves take physical values, gradients and curls; their interface
does not depend on finite-element degree or coefficient layout.

For smooth solutions, the target approximation orders are `k_s` for density,
each velocity component and both temperatures in H1, `k_A` for A in H(curl),
and `k_A` for recovered J in L2. H1 and H(curl) include the field value and
its first gradient or curl, respectively. The spectral frequency targets are
`2 k_A`, hence two and four for equal degrees one and two. These remain
targets until the complete spectral ladder passes. The zero continuous
Coulomb multiplier has no nonzero reference error from which to infer a rate.

## Operator contract

`vrmhd_section_build` constructs the spaces. `vrmhd_section_residual` returns
`M(U) U_dot + R(U)` and its directional derivative with respect to both
coefficient arguments. Parameters, geometry and prescribed sources are fixed.
`vrmhd_section_mass` returns `M(U) U_dot` through the same generated mass
leaves. `vrmhd_section_jacobian` assembles `R_U + (M U_dot)_U + shift M`
from those JVP leaves, for small reference problems. `vrmhd_section_vjp`
applies its sparse transpose. The duality is the Euclidean coefficient/test
pairing; physical Riesz maps belong to the consumer.

Sources have six channels at each quadrature point: `Js(1:3), Srho, Qe, Qi`,
with units and signs from the weak-form document. Returned rows include every
boundary test. The `fixed` mask identifies scalar boundary values and the
tangential traces of both A and J; it is available for consumers that impose
Dirichlet data, and is not applied by the volume operator. Ideal-wall assembly
supplies both trace constraints. Current is `J_wall + J0`, with `J0` in
H0(curl); its wall coefficients are data and its interior coefficients satisfy
weak Ampere. Constrained consumers use only zero-tangential tests.

## Boundary and stationary gauge

For positive resistivity, `v.n=B.n=0` implies `(v cross B)_t=0`.
Ohm's law with `n cross (E-V h)=0` prescribes `J_t=Js_t+(V h)_t/eta`.
For ideal studies the equilibrium supplies `J_t`. `vrmhd_current_lift` extracts
prescribed wall coefficients and sets every interior coefficient to zero.
The same total current enters every magnetic and thermal term. Both induction
and Ampere use H0(curl) tests, so Ampere has no surface residual. The explicit
power of the prescribed affine current constraint and checked literature
locations are in [the weak-form ledger](vrmhd-weak-form.md#conducting-wall-current-and-its-power).

`vrmhd_stationary_build` uses an H1 scalar multiplier with zero boundary
trace and the matrix `G_ij = integral W_i.grad(lambda_j) dx`. The stationary
residual adds `G lambda` to induction and `G^T A` as the gauge equation.
This is the weak Coulomb condition `div A = 0`. FortFEM provides
lowest-order graph tree–cotree metadata; the mixed form also covers the
higher-order scalar gradient kernel without a separate graph reduction.
The multiplier of a compatible continuous solution is zero. Its computed
value and Jacobian nonsingularity require separate checks.

The stationary singular-value gate removes fixed trace coefficients and
eliminates the algebraic current with its L2 mass matrix before measuring the
primal Jacobian. The physical Riesz matrices are H1 for the scalar, velocity
and multiplier fields and H(curl) for `A`, with unit reference scales on both
trial and test sides. If `H=L L^T` is this Gram matrix and `S` the reduced
Jacobian defined in [the weak form](vrmhd-weak-form.md), the gate measures
`L^(-1) S L^(-T)`. It requires a minimum singular value above `1e-4` and
rejects a decrease greater than 25 percent per halving of cell size. It also
reports the maximum singular value and minimum mode's field composition.
The current mass matrix is checked separately after diagonal normalization:
a condition-number increase greater than `1/0.75` per refinement fails.
The raw coefficient condition number is diagnostic only because section edge
and orbit nodal bases have different mesh scaling. These are numerical gates,
not proofs of uniform continuous bounds.
The Jacobian is evaluated after Newton convergence of the manufactured
stationary problem, including its lifted boundary values and prescribed current trace.

The wrapper applies `fixed` to scalar Dirichlet values, all velocity
components, and the tangential traces of A and J. The homogeneous velocity
case is a no-slip ideal wall; nonzero manufactured traces describe a prescribed
boundary-value problem. Only interior current tests remain free.
Boundary values are lifted before volume evaluation; fixed rows are identity
rows and their columns are removed from free equations. The load is an
assembled volume right-hand side, with no Ampere surface contribution.
Prescribed loads have zero state derivative. The sparse Jacobian and
matrix-free derivative use the same lifting and mixed gauge. Time-dependent
assembly keeps Weyl gauge and its vector-potential mass row.

## Driven-state conventions

A current source uses `eta (J - Js)` at zero loop voltage; this is the
JOREK comparison convention. A loop-voltage drive instead uses a harmonic
one-form `h`, normalized by `integral_loop h = 1`, and voltage `V` in volts.
Rotation uses `h = grad(phi)/(2 pi)`. Translation and screw symmetry use
`h = grad(z)/L_z` on a domain with axial period `L_z`. The latter is
invariant under screw motion and curl-free. The screw field `xi/|xi|^2`
has nonzero curl and is unsuitable for this drive.

In Weyl gauge write `A = Atilde - t V h`. Then `B = curl Atilde` and
`Atilde_dot = V h + v cross B - eta (J - Js)`. The induction residual
therefore gains `-V h.W`, and the magnetic power ledger gains `V h.J`.
A stationary Ohmic state satisfies `eta (J - Js) = V h`. The two limiting
conventions are `Js=J, V=0` and `Js=0, eta J=V h`.
For GS consistency, the loop-voltage convention prescribes the electron heat
sink `Qe=-eta |J_GS|^2` from the analytic GS current at quadrature points to
balance Ohmic heating; the ideal and prescribed-current conventions keep `Qe=0`.

Discrete harmonic coefficients must be edge integrals of the multivalued
coordinate, divided by its period, with `curl_h h=0` to roundoff. A
projection with a small nonzero curl does not meet this contract.
`vrmhd_drive_build` selects `VRMHD_CURRENT_SOURCE` or `VRMHD_LOOP_VOLTAGE`:
the former retains the prescribed `Js` and requires zero voltage; the latter
sets `Js=0` and returns the positive voltage load subtracted by the stationary
wrapper. Other source channels are retained.

The returned harmonic cochain stores zero increments on section edges and a
unit integral on every periodic orbit edge. It uses the full coordinate
one-form, not the horizontal `A/J` decomposition above. Its metric-free
FortFEM curl is zero, and its prescribed Cartesian lift divides by the orbit
period (`2 pi` or the supplied axial period). This exact affine harmonic
component of `Atilde-t V h` is evaluated separately from the polynomial
horizontal state space: projecting it into that space would destroy exact
closedness for screw symmetry.

## Verification limits

The equal-degree default passes the energy ledger and the earlier RM2 gates
described below at degree one. The spectral gate remains incomplete: its
first failure is folded-count equality at 16 cells per period on the
structured degree-two mesh. Fourier attribution does not support treating its
three positive frequencies beyond the base-wavenumber branches as one
folded-alias family. The complete structured degree-one ladder passes; degree
two and the perturbed meshes remain unqualified, as detailed below.

The `rm2` tests compare the complete assembled volume residual, including mass,
with independently written numerical weak forms; check second-order centered
finite differences, sparse versus matrix-free action, transpose duality and
one/four-thread agreement; and compare the invariant geometry with independent
Cartesian isometry flows. A separate FortSym oracle differentiates explicit
three-dimensional nonpolynomial fields directly in Cartesian coordinates to
obtain manufactured strong sources. Three section refinements check H1 norms
for each scalar/velocity component, H(curl) norms for A/J, and decreasing weak
residual defects after natural-norm projection. These are representation and
residual consistency checks, not convergence of a solved nonlinear boundary
problem.
The reported weak defect uses diagonal natural-norm scaling, not a certified
dual norm. Separate stationary tests check weak Ampere on zero-tangential
tests against the independent analytic current and magnetic field, and the lifted/gauged
sparse Jacobian against its matrix-free action and centered differences.
The solved manufactured test uses H1 for scalar and velocity fields,
H(curl) for `A` and L2 for the weakly recovered `J`, as explained in the
weak-form document. It passes over three refinements in all three symmetries,
with quadratic Newton convergence. The manufactured problem retains its
analytic fields and supplies their current trace as prescribed data; it is a
lifted boundary-value problem, not a homogeneous conducting-wall solution.
Its diagnostic interior `curl J` rate remains about 0.4–0.5; only the L2
current error is a convergence gate.
The harmonic gate checks closedness and an independent Cartesian voltage load.
The reduced stability and diagonally normalized current-mass gates pass on
meshes with 4, 8 and 16 cells per side for rotation, translation and screw.
With prescribed current trace the smallest reduced singular value exceeds
0.088; its worst refinement ratio exceeds 0.935, while the largest normalized
mass-condition growth remains below 1.317.

The GS consistency gate solves the RM0 flux equation with analytic pressure
profiles, projects its flux lift into the section spaces and recovers current
by weak Ampere's law with the prescribed equilibrium current trace. Rotation
uses `psi=0.1 R^2(log R-1/2)`, `Bz=0.2 log R`, `Jphi=-0.2/R` and
`p=6-0.02(log R)^2`. Screw symmetry has pitch `k=0.7/m`,
`psi=0.25(x^2+y^2)`, `B=(-0.2 y,0.2 x,1)` and `J=(0,0,0.4)`;
its pressure is `p=6-0.04(x^2+y^2)`. Coordinates and field values use SI
units and unit reference scales; both species have constant temperature
`Te=Ti=1 J` and `mi=mu0=1` in their respective SI units, with `rho=p/2`.
The pressure is a function of the monotone radial flux in each case.
The full gauged residual, measured in the dual of the constrained physical
Riesz norm, converges at least at first order over three refinements in the
ideal limit, with `Js=J_GS`, and with loop voltage and the prescribed heat sink
for both symmetries.

The voltage drive retains the specified Ohmic heating term and preserves the
separately prescribed heat sink. The skew momentum split closes the
hydrodynamic energy ledger, including its explicit mass-source/diffusion
coupling. The full magnetic ledger includes the explicit power of the
prescribed current lift, as explained in the weak-form document.

## Structural spectral gate

The linear-spectrum test uses the production Jacobian and mass in a periodic
Bloch block of the uniform translation slab, with the equilibrium and
`uniform_waves` oracle specified in the pairing study below. A Bloch field
repeats between blocks with the phase of its wave vector. The gate retains
every coefficient; it does not discard cell-scale modes.

Write the balanced spatial operator as `L=-R_U`, where `R_U` is the residual
Jacobian. An invertible change of units scales A by `sqrt(mu0)` and J by
`1/sqrt(mu0)`. The algebraic current is eliminated through its mass matrix,
giving the physical time mass `M_r` and the spatial Schur complement `L_r`.
The evolution generator is `F=M_r^-1 L_r`. For this equilibrium the conserved
quadratic wave energy is

```
Q = (1/2) integral [rho0 |v'|^2 + |curl A'|^2/mu0 + |p'|^2/(gamma p0)] dx,
p' = 0.6 rho' + Te' + Ti',  rho0=1, gamma p0=1.
```

Primes denote perturbations; density is in kg/m^3, pressure in Pa,
temperatures in joules, and `mi=1 kg` in the reference normalization.
The Hermitian Gram matrix `E` represents `2 Q` in the balanced coefficients.
The first gate checks inverse-free block identities at relative Frobenius-norm
tolerance `1e-14` on every mesh. Here `*` means conjugate transpose; transpose
is equivalent in real sine/cosine coordinates. In the balanced `L`, the
nonzero ideal dynamical couplings are continuity/velocity `L_rho,v`,
momentum/density and temperatures `L_v,rho`, `L_v,Te`, `L_v,Ti`,
each temperature/velocity `L_Te,v`, `L_Ti,v`, Lorentz `L_vJ` and induction
`L_Av`. The algebraic Ampere row contains `L_JA` and `L_JJ`. The background
velocity and current vanish, so momentum advection and magnetic self-couplings
vanish; the test checks that no remaining block is unpaired.

1. Lorentz and induction require
   `||L_vJ + L_Av^*||_F / ||L_vJ||_F <= 1e-14`.
   Magnetic energy tests induction with current, before eliminating current.
   For pressure/compression define `C_v=(2/3)(L_Te,v+L_Ti,v)`:
   adding the two internal-energy equations gives the pressure rate directly,
   without a spatial mass inverse. Each scalar field `f=(rho,Te,Ti)` requires
   `L_v,f + beta_f C_v^*=0`, with `beta=(0.6,1,1)`, at the same relative
   tolerance. The separate species advection/compression blocks must each
   equal `0.75 L_rho,v`, and `C_v=L_rho,v` in this normalization.
2. The independently assembled energy Grams must agree with the production
   masses and Ampere stiffness. The scalar Gram `S` comes from the finite-element
   L2 mass, independently of the production time-mass leaves. The kinetic
   blocks equal `S`; pressure energy is `beta beta^T` times `S`.
   The primitive species mass diagonal is `1.5 S`, with density-rate block
   `0.45 S`. Induction time mass equals the Ampere current mass `-L_JJ`,
   and magnetic energy equals the Ampere curl-curl stiffness `L_JA`.
   These comparisons use relative tolerance `1e-14`.
3. The energy Gram must be Hermitian at relative tolerance `1e-14`. Its kernel
   must equal the analytic entropy and curl-free gauge space `Z` described below:
   each column `Z_j` requires `||E Z_j||/(||E||_F ||Z_j||) <= 1e-14`, and
   Cholesky must succeed on the complementary energy Gram. The smallest
   complementary eigenvalue relative to the largest is reported.
   The primitive time mass `M_r` must admit a lower-upper (LU) factorization;
   its estimated 1-norm condition number must grow by at most a factor of two
   per halving of cell width, separately for each degree and mesh family.

The primitive mass is legitimately nonsymmetric. The generated species mass
is `3/2 (rho Ts_dot + Ts rho_dot)/mi`, where `Ts` is either species temperature,
so at this equilibrium `M_Ts,rho=0.45 S` while `M_rho,Ts=0`.
Its invertibility, the paired-block identities and Gram consistency give:

```
Inverse-free adjoint pairs and matching Grams give E F + F^* E = 0, F=M_r^-1 L_r.
On the E-quotient, F_q is skew-adjoint in the positive inner product x^* E_q y.
Hence F_q has imaginary spectrum and i F_q has real angular frequencies.
```

Here `F_q` and `E_q` are the induced evolution generator and energy Gram after
quotienting static energy-free directions. In particular,
`E_q^(1/2) i F_q E_q^(-1/2)` is Hermitian. This argument requires neither
symmetry of the primitive time mass nor positivity of the full energy Gram.

For a block with `r` cells in each periodic direction and equal degree `k`,
the scalar dimension is `s=(r k)^2`. At each scalar coefficient the entropy
directions are `(rho',Te',Ti')=(1,-0.6,0)` and `(0,1,-1)`; both have zero
pressure perturbation. The `s` gauge directions are discrete gradients of
the scalar Bloch basis, projected into the same vector-potential block.
The Bloch phase is nontrivial on every tested mesh, so the scalar gradient is
injective and no constant or harmonic gauge direction is added. Thus the
expected energy nullity is `3 s`: structured blocks (`r=1`) have 3 and 12
nulls at degrees one and two; perturbed blocks (`r=2`) have 12 and 48.
These counts are independent of the number of cells per period in the ladder.
A QR factorization into unitary and upper triangular matrices checks the
analytic basis rank and provides its Euclidean complement as quotient
representatives. Since `E Z=0`, every vector is E-orthogonal to `Z`;
E-orthogonality alone does not select a complement.
Positivity on the chosen complement proves that there are no additional
energy nulls. All analytic nulls must also be static under the full spatial
operator, with algebraic current restored.

The formed-product defect `||E F + (E F)^*||_F / ||E F||_F` and
`100 epsilon kappa_1(M_r)` are diagnostics only. Machine epsilon `epsilon`
refers to double precision; `kappa_1` is LAPACK's `zgecon` 1-norm condition
estimate from the LU factor used for `F`. Forward products involving
curl-curl entries growing as inverse cell width squared can amplify cancellation
roundoff even when this mass condition estimate remains bounded.
After the structural gate passes,
a Hermitian eigensolve uses the roundoff Hermitian representative on the
positive-energy quotient. The nonsymmetric solve's
`max |Im(omega^2)| / max |omega^2|` remains diagnostic.

The second gate places the cutoff in a folded-oracle spectral gap. The largest
oracle frequency to be counted is the largest in the retained band no greater
than the base fast-wave frequency `sqrt(1.8) rad/s`. The cutoff is its midpoint
with the next larger folded-oracle frequency; both endpoints and the cutoff
are reported for every block. The one-sample structured degree-one band has
no larger frequency, so adjacent aliases with integer indices `jx,jy=-1,0,1`
provide the upper endpoint without adding branches to the retained band.
The computed signed nonzero count below that cutoff must equal the folded
oracle count on resolved meshes: those with at least 16 cells per period in
the 4, 8, 16, 32 ladder. One-to-one branch matching, including repeated
frequencies, applies to the same resolved meshes. Coarse counts are diagnostics.
For a block of `r` cells in each direction, its wave vectors are
`(1 + n jx/r, n jy/r, 0) m^-1`, where `n` is the number of cells per period.
This is `q0 + 2 pi j/(r h)`, with base wave vector `q0=(1,0,0) m^-1`,
integer vector `j=(jx,jy,0)` and cell width `h=2 pi/n m`.
The resolved test band has `r k` samples in each direction, with integer
indices `-floor(r k/2)` through `floor((r k-1)/2)`. The structured mesh uses
`r=1`; the periodically perturbed mesh uses `r=2`. Matching retains the prior
isolation window: the frequency error must be smaller than half the distance
to the nearest distinct oracle frequency. A bipartite assignment preserves
the folded multiplicity instead of independently accepting nearest branches.
The ladder is 4, 8, 16, 32 cells per period at `k=1,2`, on both meshes;
the base-wavenumber frequency order must exceed `2 k-0.2` on the finest pair.
The two coarse meshes retain nearest positive frequency estimates in the ladder
without asserting count equality or branch identity. The isolation window is
unchanged.

The third gate counts energetic static modes per cell against the study's
`k`, separately from entropy and gauge. Every static direction, with its
algebraic current restored, must satisfy
`||L x|| / (||L||_F ||x||) <= 1e-12` in the same balanced units.
For magnetic diffusivity `eta/mu0=1e-3 m^2/s` and kinematic viscosity
`nu/rho0=1e-3 m^2/s`, the full dissipative generator must have no positive
real eigenvalue beyond `1e-12 ||F_d||_F`, where `F_d` denotes that generator.
All energetic modes, including the ideal static kernel, must acquire decay;
entropy and curl-free gauge remain static. The slowest energetic decay is
reported in s^-1 and relative to `(eta/mu0) k_cell^2`, with
`k_cell=2 pi/h`. The slowest rate is taken over the whole energetic quotient.

The structured `k=1`, four-cell case exceeds the preserved branch-isolation
window. Its Alfven frequency is approximately
`0.6994805151 rad/s`, against `0.6 rad/s`; their difference exceeds half the
Alfven/slow-oracle separation. The number of signed waves below the cutoff
is correct, but only four of six have an isolated oracle match. This is a
coarse-mesh dispersion observation, not evidence by itself of growth or a
spurious branch; isolation is required only on the resolved meshes.
Test-side Bloch projections and structural identity reductions accumulate in
quadruple precision (`real128`), converting the double-precision inputs exactly.
The projected matrices and LAPACK solves remain double precision. This reduces
summation roundoff from the comparisons without changing their `1e-14`
tolerances or the production operators.
The structured degree-one ladder passes through 32 cells per period, including
the finest-pair order check. Structured degree two passes both coarse meshes,
whose counts are diagnostic. At 16 cells per period its structural identities,
analytic nulls, quotient positivity and time-mass inversion pass, but twelve
computed signed waves lie below the gap cutoff, against ten analytic waves.
The next larger folded frequency belongs to the slow branch at
`(-15,-16,0) m^-1`, so the gap endpoints are approximately `1.341640786500`
and `6.508898217884 rad/s`, with midpoint `3.925269502192 rad/s`.
This is a resolved-mesh count failure away from an oracle boundary.

A diagnostic projects each below-cutoff mode's physical wave-energy fields
`(v', curl A'/sqrt(mu0), p'/sqrt(gamma p0))` onto
`exp(i (q0+G).x)`, where `G` is a retained reciprocal-block offset,
using the existing element quadrature. It sums the squared
L2 Fourier coefficients over those seven components, then normalizes over the
retained alias set to obtain each alias's share. It separately reports the
fraction of the full physical L2 norm captured by that set. Thus a large
retained share does not imply that the retained set captures the entire mode.
These fields avoid arbitrary entropy and vector-potential gauge content;
the balanced A coefficients already include the `sqrt(mu0)` scaling.
An independent test uses exactly representable `p'=x` and `A_z'=x^2/2` on
`[0,2 pi]^2`. Their x Fourier coefficients are `i/m` up to the magnetic
component's sign, for nonzero integer x wavenumber `m`; coefficients with
nonzero y wavenumber vanish. The test checks both shares and captured power.

The hypothesis that the three positive modes beyond the base-wavenumber
branches all belong to `G=(0,-16) m^-1` is refuted by this diagnostic. At 16 cells per period,
that alias dominates the `0.614683369776 rad/s` mode, but the
`1.131816854126 rad/s` mode splits its retained power between that alias
and `G=(-16,0) m^-1`; the latter alias dominates the
`2.624802298648 rad/s` mode. The diagnostic reports every retained share
and the captured fraction for both frequency signs.
The conditional replacement of count equality by dominant-alias matching
is therefore not applied. In particular, no Fourier-share threshold or
matching gate is relaxed. The origin of these optical or mixed modes remains
unresolved; dominant Fourier power alone does not prove their branch identity.
The ladder stops at the unchanged count gate; no complete spectral
qualification is claimed for degree two or either perturbed ladder.
No full-system inf-sup bound or continuous spectral non-pollution result is
claimed. The independent folding test uses analytic frequencies at
`(1,0,0)`, `(-1,0,0)`, `(1,-2,0)` and `(-1,-2,0) m^-1`.

## Pairing study for the ideal uniform equilibrium

Equal degree is the default because it minimizes the measured energetic
static kernel among the tested Lagrange pairings. This choice does not waive
the structural spectral gate above.
Here pairs are written **`(k_s,k_A)`**, with continuous Lagrange degree `k_s`
for the fluid and multiplier, and first-kind Nedelec order `k_A` for the planar
vector potential and current. Their orbit components remain continuous degree
`k_A`. The [complete table](data/vrmhd-pairings.csv) covers degrees one through
four with `k_s >= k_A-1`, on both meshes and both resolutions.
On the structured mesh, every pairing retains energy-bearing static defects
at both resolutions. Equal degrees one through four give respectively one,
two, three and four magnetic nulls in the static quotient. Taking
`k_A=k_s+1` increases this count to eight, thirteen and eighteen; taking
`k_s>k_A` instead gives four to twenty-eight kinetic nulls. Refinement can
reduce the physical-frequency errors without removing these defects. The
magnetic and kinetic counts per cell are also unchanged by the prescribed
mesh perturbation at either resolution.

The equilibrium has mass density `1 kg/m^3`, pressure `0.6 Pa`, adiabatic
exponent `5/3`, and Alfven velocity
`B0/sqrt(mu0*rho0)=(0.6,0,0.8) m/s`. The wave vector is
`k=(1,0,0) m^-1` on a square period of side `2 pi m`. The analytic Alfven,
fast and slow angular frequencies are respectively `0.6`, `sqrt(1.8)` and
`sqrt(0.2) rad/s`. Each direction has four or eight cells per period; one cell
contains two triangles. Degrees of freedom per cell include the algebraic
current and are `6 k_s^2 + 2(3 k_A^2+k_A)`.

The second mesh repeats a two-by-two displacement pattern. In units of the
cell width, the vertex displacements for parity `(0,0)`, `(1,0)`, `(0,1)` and
`(1,1)` are `(0,0)`, `(0.13,-0.11)`, `(-0.09,0.16)` and `(0.18,0.07)`.
Periodic boundary copies receive the same displacement. Its Bloch block retains
four cells and all their coefficients; the structured block retains one.
Unequal degrees remain an explicit study option. Ordinary construction
defaults to `k_s=k_A` and rejects unequal degrees.

The table preserves the historical generalized eigensolver, zero tolerance
`1024 epsilon`, cutoff `sqrt(1.8) rad/s`, nearest-branch matching and allowance
of three null modes. It records failures of those historical reality,
finite-count and signed-wave diagnostics instead of stopping the measurement.
These raw flags are not the restated structural spectral gate. `unmatched`
is the historical audit's count, including
zero multiplicity beyond three. It does not count extra eigenvalues that the
existing nearest-branch rule accepts. Additional cell patterns in entropy and
gauge fields, and folded Fourier branches in a supercell, prevent interpreting
every extra eigenvalue as spectral pollution.

`alfven_error`, `fast_error` and `slow_error` are dimensionless relative errors
of the nearest positive frequency to each analytic frequency. They measure
frequency approximation, without establishing branch identity or absence of
extra modes. The diagnostic `omega2_imag_s_2` is the largest absolute imaginary part of
squared angular frequency in `s^-2`. Backward error and analytic-null defects
are normalized; equilibrium defect is the largest projected residual entry in
the solver reference units.
`seconds` measures one complete case, including assembly and diagnostics.

The energy split uses kinetic energy `E_K` and magnetic energy `E_B`, obtained
from the production mass and discrete curl operators. A singular-value
decomposition (SVD) of the physically scaled static coupling identifies its
nullspace. After quotienting the two analytic entropy modes and one gauge
mode, the positive part of `E_K+E_B` is whitened, then `E_B` is diagonalized.
This avoids choosing an arbitrary basis of a degenerate zero eigenspace.
`static_magnetic` and `static_velocity` count fractions
`E_B/(E_K+E_B) >= 0.9` and `<= 0.1`, respectively; `static_mixed` includes the
remaining quotient directions and directions with negligible physical energy.
The energy floor is `1e-10` times the larger of the largest null energy and the
largest velocity-mass entry, in the scaled coefficient coordinates.

`magnetic`, `velocity` and `mixed` sum to `unmatched`. They apply the same split
to unmatched finite modes and to the static quotient when the SVD and generalized
zero counts agree (`null_split_ok`). Otherwise the zero part is unresolved and
included in `mixed`; the separate static counts remain available. Thus `mixed`
is a conservative catch-all, including energy-free entropy/gauge patterns, not
solely modes with comparable kinetic and magnetic energy. `vz_blind` separately
counts singular values below `1024 epsilon` times the largest singular value
in the axial-velocity coupling. This is independent of the generalized zero
classification. The smallest shear direction is also checked against the
production matrix-free Jacobian.

The full study is registered as `vrmhd_pairing_benchmark`, with the `benchmark`
and `performance` labels and without `rm2`; the small admission/measurement
smoke test remains separate. Reproduce the table with `BLA_VENDOR=OpenBLAS`,
eight OpenMP threads and one OpenBLAS thread:

```sh
cmake --build build --target vrmhd_pairing_study -j 8
OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 build/vrmhd_pairing_study docs/data/vrmhd-pairings.csv
python3 docs/plot_vrmhd_pairings.py --output /tmp/vrmhd-pairings.png
```

The optional [plotting script](plot_vrmhd_pairings.py) requires NumPy and Matplotlib. It writes a raw
count/frequency comparison and a companion `-energy` plot of the static split;
generated images are not tracked.

### Literature and scope of the compatibility conditions

`H(curl)` and `H(div)` mean square-integrable fields with square-integrable curl
and divergence, respectively; `H1` requires square-integrable first derivatives.
Raviart--Thomas (RT) and Brezzi--Douglas--Marini (BDM) elements have continuous
normal traces; Nedelec elements have continuous tangential traces. `CG` and
`DG` denote continuous and discontinuous polynomial spaces. Locations below
were checked in accessible text; the two original-source access gaps remain.

- [Rappaz (1977), Numerische Mathematik 28, 15--24](https://doi.org/10.1007/BF01403854):
  the publisher abstract states a sufficient non-pollution condition, without
  specifying the pairing. The full theorem was unavailable.
  [Ida, Todoroki and Sanuki, NIFS-PROC-040, p. 108 after (20)](https://www.nifs.ac.jp/report/NIFS-PROC-040.pdf#page=116)
  attribute continuous degree `p` derivative-bearing displacement and
  discontinuous degree `p-1` algebraic displacement, `p>=1`, to Rappaz.
  This is a secondary attribution for a one-dimensional formulation.
- [Gruber and Rappaz (1985), Section 4.1.1, p. 80, (4.8), publisher preview](https://page-one.springer.com/pdf/preview/10.1007/978-3-642-86708-8_5):
  continuous linear triangles for all displacement components give a deficient
  discrete Alfven constraint. Section 4.2's full hybrid pairing and condition
  were inaccessible. [Gruber (1978), abstract](https://doi.org/10.1016/0021-9991(78)90076-1)
  describes independent displacement derivatives, auxiliary constraints, and
  cellwise-constant Lagrangian terms at lowest hybrid order; it gives no
  velocity/edge-potential degree rule.
- [Gawlik and Gay-Balmaz, incompressible MHD (2021), Section 3, (20)--(23), Proposition 3.1](https://arxiv.org/html/2012.04122):
  velocity and magnetic field use `RT_s`, mean-zero pressure `DG_s`, density
  `DG_m`, and auxiliary fields `NED_s`, with `s,m>=0`. The paper numbers the
  lowest Nedelec space by zero. Curl inclusion and the discrete divergence
  range yield pointwise divergence constraints and conservation, not a
  spectral non-pollution theorem.
- [Gawlik and Gay-Balmaz, compressible MHD (2021), Section 3](https://arxiv.org/html/2105.00785):
  velocity uses `CG_(r+1)` in `H1_0`, magnetic field `RT_r`, density `DG_s`,
  auxiliary fields `NED_r`, with `r,s>=0`. The discrete de Rham inclusion and
  balance identities concern this formulation; velocity is not in `H(div)`
  alone, and spectral correctness for vector potential is not established.
- [Claes, De Jonghe and Keppens, *Legolas*, Appendix A.1, (A1)--(A2)](https://arxiv.org/html/2010.14148):
  quadratic density, tangential velocity, temperature and normal potential;
  cubic Hermite normal velocity and tangential potentials. The stated condition
  is exact representability of incompressible eigenfunctions in a one-dimensional
  equilibrium with two Fourier directions, not a triangle degree rule.
- [Sovinec, Patil and Cauilan (2026 preprint), Appendix A, pp. 25--28, Figure 14](https://arxiv.org/pdf/2606.28613#page=26):
  cubic continuous velocity/pressure, discontinuous cubic normal potential,
  continuous quartic tangential potentials and quartic damping-gauge scalar.
  Total pressure reproduces oscillatory frequencies, but five solenoidal
  mesh-scale potential modes remain uncoupled. Equal continuous potential and
  velocity degrees add oscillatory pollution. None of the tested formulations
  is strictly stable near marginality; this is not a spectral-correctness theorem.

### What an H(div) velocity study would require

[FortFEM](https://github.com/lazy-fortran/fortfem) already supplies triangular RT degree `r>=0` and BDM degree `r>=1`
in `triangle_rt_arbitrary_order` and `triangle_bdm_arbitrary_order`; their local
sizes are `(r+1)(r+3)` and `(r+1)(r+2)`. The meridian complex's face space is
rotated first-kind Nedelec, hence RT degree `k_A-1` in the plane, with a
`DG_(k_A-1)` orbit component. BDM assembly and orientation transforms are
available through `assembly_full_vector_arbitrary_order_2d`. Existing triangle tests
exercise RT degrees zero through three and BDM degrees one through four;
this inspection does not certify every order.

Changing the velocity requires new coefficient maps, mass and divergence
couplings, boundary traces, and consistent induction/Lorentz pairings. Normal
continuity alone does not make the velocity conforming for symmetric viscous
stress: tangential jumps need consistent facet terms, such as an interior
penalty or a mixed stress formulation. Strong normal no-slip data would need
weak tangential enforcement, for example Nitsche terms or a hybrid trace.
FortFEM has scalar interior-penalty and Nitsche interface primitives, but the complete
vector stress, tangential boundary treatment and matching viscous-heating
ledger are not supplied by those primitives. Generic facet assembly belongs
in FortFEM; MHD composition, generated pointwise stress/heating derivatives,
and the energy and no-slip tests belong here. No such change is implemented.

Equal-degree admission follows the measured static counts; exactly static
kernels are acceptable only subject to the structural and dissipative gates.
An H(div) velocity formulation remains a separate possible design study, not
an implemented change. Conservation or a discrete de Rham sequence alone does
not prove spectral non-pollution of the coupled MHD system.
