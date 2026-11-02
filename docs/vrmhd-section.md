# Invariant visco-resistive MHD volume assembly

The invariant section assembler evaluates the mixed-current volume residual,
its matrix-free state Jacobian action, and the conservative mass terms of
[the Cartesian weak form](vrmhd-weak-form.md). It accepts the rotation,
translation and screw cases of `gs_symmetry_t`. This is a candidate operator;
boundary elimination, a stationary gauge, nonlinear solves and continuous
error certificates are not supplied.

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

The induction equation uses Weyl gauge. A stationary solve must separately
remove its gradient nullspace. A nonzero prescribed Ohmic field
`E_loop xi/|xi|^2` is not a zero of the present stationary induction residual.
In rotation it can represent a changing potential with stationary magnetic
field; in screw symmetry its nonzero curl also changes the magnetic field by
Faraday's law. A driven stationary gate therefore needs an explicit compatible
electric-field/gauge and boundary convention before it can be asserted.

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
dual norm. Manufactured nonlinear solves, GS stationary consistency, stationary
gauge nonsingularity, the semidiscrete energy ledger and linear spectra remain
unverified. Conservative momentum alone does not enforce the kinetic
continuity product discussed in the weak-form document.
