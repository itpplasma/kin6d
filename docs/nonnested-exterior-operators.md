# Frozen exterior operator composition

This slice recovers `nonnested_exterior.f90` from draft PR #47 at
8029fc50660c777191364cdbac2cf6708d0a4fc2 onto the current NN2 closure and
matrix-free actions. It retains the current dependency pins. It does not
qualify full BEM assembly, general interface geometry, a finite-beta free
plasma boundary, or G8. The full free-boundary work remains draft until G8
passes on a qualified dependency revision.

`nn_attach_exterior` retains the existing swept-meridian interface and coil
construction using [FortFEM](https://github.com/lazy-fortran/fortfem).
`nn_attach_exterior_operators` factors only its frozen algebra. Its caller
supplies a qualified positive symmetric exterior energy matrix S, the full
physical normal-moment map T, the coil normal moments g, and the interior
coil volume load. Finiteness, dimensions, symmetry and the restriction to
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
m. The physical exterior magnetic energy is

    (T a - g)^T S^(-1) (T a - g) / (2 mu0).

NN2 uses the toroidal mean, so the entire quadratic energy and both of its
derivatives contain exactly one factor 1/(2 pi). In particular,

    H = T^T S^(-1) T / (2 pi),
    d = T^T S^(-1) g / (2 pi).

The residual gains `(H a - d - volume_drive)/mu0`. Its Jacobian and JVP
gain `H/mu0`, and the VJP gains its actual transpose. Fixed-row masking and
the reverse of row reduction follow the existing NN2 coefficient pairing.
No new constitutive expression is introduced by this linear composition.

The response is a single-valued scalar gradient with zero periods. The
consumer therefore retains the zero-net-current I_phi restriction and rejects
floating internal conductors. A nonzero plasma current requires an explicit
topological harmonic response. The known coil field is kept separate.

## Independent oracle and its limits

`nonnested_exterior_frozen_sphere_oracle` supplies synthetic boundary-column
maps in the trace basis n_x, n_y, n_z of a unit sphere. The exact exterior
potential is u = (c dot x)/r^3 and its field is

    B_response = [3 (c dot n) n - c] / r^3.

Thus S = 2 M with M = (4 pi/3) I. A nonzero coil trace and a nonzero state
produce a nonzero response. Six-axis angular integration is exact for the
degree-two energy integrands. The independent radial integral is
`integral_1^infinity r^(-4) dr = 1/3`. The test compares this physical energy
with the Schur energy, checks its derivative, and checks normal and tangential
signs, zero net flux and a closed-loop gradient period. The supplied T is
synthetic and does not test a toroidal FEEC trace construction.

P2/P3 `nonnested_actions_exterior` cases compare the attached operator's JVP
with the assembled Jacobian and central differences, and its VJP with the
assembled transpose and Euclidean pairing. These are operator checks at
nonstationary states. They are not equilibrium solutions, continuous error
estimates, rigorous error bounds or GPU benchmarks.

## Fitted current and the vacuum shortcut

An allocated frozen fit excludes the pure-vacuum elimination of pressure,
Ohmic potential and loop voltage, even when its scalar and weak pressure loads
vanish. The public harmonic regression uses B_phi=F/R and d_seed=0.37: zero
net toroidal current requires V=-2*pi*F*d_seed/g_par. This checks cancellation
of the frozen drive by the Ohmic current, rather than reproduction of the
original target lambda. At the exactly harmonic state the declared perturbation
D3 scale is zero, so its normalized force remains explicitly unavailable.

The pressure row is R_p(B,p)=A_p(B)*p-F_p, where F_p includes the complete
scalar and frozen weak-flux drive on free rows. If that assembled drive and
all pressure seed/boundary coefficients are exactly zero, p=0 is an invariant
for every admitted nonnull B. Its derivative with respect to the other fields
is zero there. The existing reduced Newton system can therefore keep only
those pressure increments zero, avoiding roundoff leakage from mixed LU. This
uses exact zero tests and preserves all original residual and merit gates.
Electric potential, loop voltage and signed-current equations remain live,
including when a frozen fit is allocated. Nonzero frozen pressure drive is
covered by a separate forced control and excludes this reduction. This is an
algebraic subspace of the declared model, not a vacuum-collar support law.
