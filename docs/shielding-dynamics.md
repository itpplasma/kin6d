# Planar Debye shielding relaxation

The Chapter 2 transient is a **linear, overdamped continuum** experiment.
It does not represent microscopic Coulomb particles, collisionless Vlasov
response or spherical screening. It complements the separate matched
three-dimensional charged-sphere reference.

## Named model

Lengths are normalized by the electron Debye length, time by
`lambda_D^2/D`, potential by `k_B T/e`, and field by `k_B T/(e lambda_D)`.
The periodic domain is `[-pi,pi]`. A fixed external charge perturbation has
exactly supplied, finite cosine coefficients `rho_k`; its zero mode is
compensated. The electron excess is `eta=delta n_e/n_0`, initially zero.

Linearized electron drift-diffusion with the Einstein relation and unscreened
Poisson coupling gives

\[
\eta_t=\eta_{xx}-\eta+\rho_{\rm ext},\qquad
-\phi_{xx}=\rho_{\rm ext}-\eta,\qquad E=-\phi_x.
\]

The potential has zero mean. The fixed source and thermal reservoir make this
an open, overdamped model; energy conservation of a closed particle system is
not asserted. Constant temperature, constant diffusion/mobility, small density
and potential perturbations, and instantaneous electrostatics are declared
physical assumptions. No Yukawa interaction or prescribed screening cloud
enters the time evolution. At equilibrium the coupled model produces the
periodic linear Debye/Helmholtz response.

For each nonzero cosine mode,

\[
\dot a_k=\rho_k-\mu_k a_k,\qquad \mu_k=k^2+1,
\qquad \phi_k=(\rho_k-a_k)/k^2,\quad E_k=(\rho_k-a_k)/k.
\]

The application uses a generated explicit Euler candidate, with maximum
normalized timestep `0.0025`. All mode, derivative, residual, enclosure and
coordinate leaves are retained in `gen/gen_shielding_dynamics.f90` and emitted
into `src/generated/kin6d_shielding_dynamics_leaves.f90`.

## Continuous numerical certificate

The checker reconstructs the actual floating-point endpoint coefficients by
an affine segment on the exact represented endpoint times. Its residual is

\[
r_k=\dot a_{h,k}+\mu_k a_{h,k}-\rho_k.
\]

This residual is affine in time, so strict endpoint intervals bound its
supremum on the complete segment. For an incoming error radius `e_k`, variation
of constants and positive damping imply

\[
e_k^{\rm new}\le e^{-\mu_k h}e_k+h\sup|r_k|.
\]

All arithmetic in this recurrence and residual evaluation is outward interval
arithmetic. Candidate arithmetic is untrusted, and its roundoff appears in the
reconstructed residual. The checker never receives the exact time solution.
The initial radius is zero because the specified initial state is exactly zero.
A failed partial advance invalidates the certificate packet; it cannot be
sampled as a certified result.

Summing the modal radii gives global pointwise bounds for electron excess,
potential and electric field, with the potential/field radii weighted by
`1/k^2` and `1/k`. Strict modal sampling encloses reconstruction roundoff.
Numerical bands cover this finite-Fourier continuous model; there is no omitted
source tail because the supplied finite coefficient vector defines the source.
A different external source would need its own reconstruction bound.

Straight plot segments receive the `M Delta x^2/8` remainder from the strict
candidate curvature, in addition to the pointwise numerical radii. Thus the
visible spatial bands enclose the whole curve between plot nodes. The potential
history receives an analogous temporal chord margin. For zero initial state,
`|phi_k''| <= (k^2+1)|rho_k|/k^2`, since the derivative itself obeys a contractive
modal equation. This bound does not require evaluating the exact time profile.

## Derivatives, verification and evidence

Smooth modal production leaves expose generated primal/JVP/VJP using Euclidean
input/output pairings. A derivative consumer would solve the equation-level
linearized evolution; no continuous sensitivity enclosure is presently
claimed. Model-reduction and uncertain physical-input bounds are also
**unavailable**, separately from the supported continuous numerical bounds.

The registered `shielding_dynamics_oracles` test checks a real128 exact ODE
oracle for independent signed Fourier inputs, electron/potential/field interval
containment, midpoint display chords, Euler refinement, equilibrium Helmholtz
balance, finite-difference derivatives, JVP/VJP pairing and failure semantics.
The analytical time solution and equilibrium formula are independent benchmark
oracles, not inputs required by the general transient checker.

`shielding_dynamics_demo` writes:

- `output/debye-shielding.png` and `output/debye-shielding.mp4`;
- `output/shielding-dynamics.csv`, including full spatial numerical bands;
- `output/shielding-performance.csv` with per-frame solve-only and certified
  step costs, extra seconds and multipliers;
- matched setup-plus-solve times for candidate-only and certified evolution,
  extra seconds, their ratio, and separate checked solve/sampling cost.

The timing receipt is a single local CPU trial and excludes plotting/encoding
from both solve baselines. It is not a speed ranking. `--snapshot [frame]`
retains the same numerical evidence and renders a selected frame without a
movie encoder. The 121 frames span normalized times 0 through 4 at 24 fps.
Generated images, movies and CSV evidence remain outside source history.
