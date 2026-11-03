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
vector-potential tangential trace; it is available for consumers that impose
Dirichlet data, and is not applied by the volume operator. Ideal-wall or
periodic assembly must supply its trace constraints and any surface terms.
In particular, unrestricted auxiliary-current tests need the Ampere boundary
term if the magnetic tangential trace is nonzero.

## Boundary and stationary gauge

`vrmhd_ampere_boundary` integrates the residual contribution
`-integral_boundary (n cross B).K dS`, including every auxiliary-current test.
It accepts a prescribed invariant magnetic trace. The volume and surface
terms together represent `mu0 J - curl B`. Surface integration uses the
outward conormal and measure of the actual swept boundary; in screw
coordinates the section normal alone is insufficient.

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
stationary problem, including its lifted boundary values and Ampere trace.

The wrapper applies `fixed` to scalar Dirichlet values, all velocity
components, and the vector-potential tangential trace. The homogeneous
velocity case is a no-slip ideal wall; nonzero manufactured traces describe
a prescribed boundary-value problem. Auxiliary-current tests remain free.
Boundary values are lifted before volume evaluation; fixed rows are identity
rows and their columns are removed from free equations. The load is an
assembled right-hand side: subtract the Ampere surface residual from that
load. Prescribed loads have zero state derivative. The sparse Jacobian and
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
dual norm. Separate stationary tests check the Ampere surface sign against
the independent analytic current and magnetic field, and the lifted/gauged
sparse Jacobian against its matrix-free action and centered differences.
The solved manufactured test uses H1 for scalar and velocity fields,
H(curl) for `A` and L2 for the weakly recovered `J`, as explained in the
weak-form document. It passes over three refinements in all three symmetries,
with quadratic Newton convergence. The diagnostic interior `curl J` rate is
about 0.4–0.5, like the whole-mesh rate, so it is not a boundary-layer effect.
The harmonic gate checks closedness and an independent Cartesian voltage load.
The reduced stability and diagonally normalized current-mass gates pass on
meshes with 4, 8 and 16 cells per side for rotation, translation and screw.
The smallest reduced singular value is above 0.089; its worst refinement
ratio is 0.937, while the largest normalized mass-condition growth is 1.065.

The GS consistency gate solves the RM0 flux equation with analytic pressure
profiles, projects its flux lift into the section spaces and recovers current
by weak Ampere's law with the prescribed magnetic boundary trace. Rotation
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
coupling. The full magnetic ledger still fails because unrestricted recovered
current is not an admissible induction test at the fixed wall trace, as
explained in the weak-form document. Linear spectra remain unverified.
