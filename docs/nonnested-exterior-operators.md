# Frozen exterior operator composition

The nonnested equilibrium model [NN2](equilibrium-nonnested-nn2.md) can
compose a frozen exterior magnetic response with its interior equations.
The implementation below checks this linear composition. Full boundary
element assembly, general interface geometry and a finite-pressure free
plasma boundary remain unqualified under [gate G8](island-equilibrium.md#verification-ladder-before-machines).

`nn_attach_exterior` retains the existing swept-meridian interface and coil
construction using [FortFEM](https://github.com/lazy-fortran/fortfem).
`nn_attach_exterior_operators` factors only its frozen algebra. Its caller
supplies a qualified positive symmetric exterior energy matrix S, the full
physical normal-moment map T from the interior magnetic-potential
coefficient vector a, the coil normal moments g, and the interior coil
volume load `volume_drive`. Finiteness, dimensions, symmetry and the restriction to
boundary columns are checked. These checks do not certify positive
definiteness, quadrature accuracy or the physical construction of T.

## Orientation and measure

The existing consumer uses B_response = -grad(u), with

    S u = T a - g.

The exterior plus-gradient potential is therefore Phi_ext = -u. Neither
potential is the ohmic potential in NN2's electrical closure. The normal
points out of the computational interior. S is the positive exterior
energy map, whose sign differs from an exterior-domain outward derivative.
There is no sign multiplier to adjust when comparing curves.

The scalar trace basis is the basis shared by S, T and g. With dimensionless
trace functions, u has units T m, T a and g have units T m^2, and S has units
m. With vacuum permeability mu0, the physical exterior magnetic energy is

    (T a - g)^T S^(-1) (T a - g) / (2 mu0).

NN2 uses the toroidal mean: its energy is the physical energy above divided
by 2 pi. Its gradient and Hessian contain that same factor. In particular,

    H = T^T S^(-1) T / (2 pi),
    d = T^T S^(-1) g / (2 pi).

The residual gains `(H a - d - volume_drive)/mu0`. Its Jacobian and
Jacobian-vector product (JVP) gain `H/mu0`, and the vector-Jacobian product
(VJP) gains its actual transpose. Fixed-row masking and
the reverse of row reduction follow the existing NN2 coefficient pairing.
No new constitutive expression is introduced by this linear composition.

The response is a single-valued scalar gradient with zero periods. The
consumer therefore retains the zero-net-current I_phi restriction and rejects
floating internal conductors. A nonzero plasma current requires an explicit
topological harmonic response. The known coil field is kept separate.

## Independent oracle and its limits

`nonnested_exterior_frozen_sphere_oracle` supplies synthetic boundary-column
maps in the trace basis n_x, n_y, n_z of a unit sphere. The exact exterior
potential is u = (c dot x)/r^3, where c is a constant vector, r = |x|
and n = x/r. Its field is

    B_response = [3 (c dot n) n - c] / r^3.

Thus S = 2 M, where M = (4 pi/3) I is the trace mass matrix. A nonzero coil trace and a nonzero state
produce a nonzero response. Six-axis angular integration is exact for the
degree-two energy integrands. The independent radial integral is
`integral_1^infinity r^(-4) dr = 1/3`. The test compares this physical energy
with the Schur energy, checks its derivative, and checks normal and tangential
signs, zero net flux and a closed-loop gradient period. The supplied T is
synthetic and does not test a toroidal finite element exterior calculus
(FEEC) trace construction.

The polynomial-order-two and -three (P2/P3) `nonnested_actions_exterior` cases compare the attached operator's JVP
with the assembled Jacobian and central differences, and its VJP with the
assembled transpose and Euclidean pairing. These are operator checks at
nonstationary states. They are not equilibrium solutions, continuous error
estimates, rigorous error bounds or GPU benchmarks.

## Fitted current and the vacuum shortcut

An allocated frozen fit excludes the pure-vacuum elimination of Ohmic
potential and loop voltage, even when its scalar and weak pressure loads
vanish. The harmonic regression uses B_phi=F/R, with cylindrical radius R
and constant toroidal-field parameter F. For the frozen parallel-current
drive d_seed=0.37 and parallel conductivity g_par, zero net toroidal current
requires loop voltage V=-2*pi*F*d_seed/g_par. This checks cancellation of the
frozen drive by the Ohmic current. The [force residual D3](equilibrium-core.md)
is normalized using the declared perturbation field; that scale is zero at
the exactly harmonic state, so its normalized force remains unavailable.

The pressure row is R_p(B,p)=A_p(B)*p-F_p, where F_p includes the complete
scalar and frozen weak-flux drive on free rows. These frozen loads depend
only on the declared data and fixed geometry, not on the current B. If that
assembled drive and all pressure seed/boundary coefficients are exactly zero,
p=0 is an invariant for every admitted nonnull B. Its derivative with respect to the other fields
is zero there for the pressure residual at fixed p=0. The reduced Newton
system can therefore hold every pressure coefficient at zero, avoiding
roundoff leakage from mixed LU. This
uses exact zero tests and preserves all original residual and merit gates.
With a frozen fit, electric potential, loop voltage and signed-current
equations remain live. Without a fit, zero pressure and zero electrical drive
also preserve the exact zero-potential/zero-voltage vacuum subspace. Nonzero frozen pressure drive is
covered by a separate forced control and excludes this reduction. This is an
algebraic subspace of the declared model, not a vacuum-collar support law.
