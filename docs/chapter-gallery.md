# Chapter 1–3 animation gallery

The gallery accompanies Chapters 1–3 of the plasma script in
`/home/ert/proj/plasma-skript`. Its purpose is to let students see individual
particles, continuous distributions and measurable fields side by side.
The examples use declared reduced models; each result retains its own
certification status.

## Reproduce the gallery

```sh
cmake --preset cpu
cmake --build --preset cpu
ctest --preset cpu
./scripts/render-gallery.sh
```

The applications and scientific postprocessing are Fortran. FortSym emits
the mathematical leaves, FortNum provides rigorous arithmetic and the
oscillator time step, and FortPlot draws the frames. FFmpeg encodes MP4 and
looping GIF files. All generated media and numerical receipts stay under
the ignored `output/` directory.

## Models shown

| Chapter | Animation | Mathematical meaning |
|---|---|---|
| 1 | `collective-response` | Contact interaction compared with an unscreened repulsive Coulomb encounter; three-dimensional particle mechanics with a planar orbit imposed by symmetric initial data |
| 1 | `particles-to-moments` | Weighted characteristic particles and a Fourier/P1 Eulerian distribution, followed by density, current and bulk velocity |
| 2 | `debye-potential-reduction` | Bare uniformly charged three-dimensional sphere compared with its matched linear Debye–Hückel equilibrium; a parameter sweep, not a transient |
| 2 | `debye-shielding` | Computed periodic planar linear drift-diffusion–Poisson relaxation toward screening equilibrium |
| 3 | `plasma-oscillation` | A numerically integrated cold electron slab against fixed ions, with visible particles, Eulerian slab profiles and a harmonic-oscillator oracle |
| Comparison | `particles-vlasov` | Prescribed-zero-field free streaming represented by weighted particles and an Eulerian Fourier/P1 state |
| Comparison | `free-streaming` | The earlier Maxwellian-envelope free-streaming laboratory |
| 1, quantum comparison | `quantum-correlations` | Interacting opposite-spin two-electron wavefunction in a periodic 3D Coulomb box; numerical finite-basis evolution with continuous bounds and a computed noninteracting control |
| 1, classical comparison | `classical-coulomb-pair` | Independent actual two-electron trajectories in the same 3D box and preparation; particle projections and matching relative-position slab bins |

The bounded-velocity distribution in the moments example is specified by its
finite Fourier/P1 coefficients. It differs from the script's illustrative
unbounded Maxwellian-plus-beam profile. No omitted velocity tail is silently
included in its certificate.

The [matched pair experiment](quantum-teaching-slice.md) is a microscopic
cold-quench example, not a hot bulk-degeneracy or Debye-screening simulation.
Its quantum bands cover the declared finite basis; its classical intervals
describe sampling uncertainty. Render these two additions alone with
`./scripts/render-gallery.sh --only quantum_reference_demo,classical_pair_demo`.
The 125/343-mode comparison reports a basis-convergence indicator separately.

The shielding relaxation solves, in normalized units,

\[
\partial_t\eta=\partial_x^2\eta-\eta+\rho_{\rm ext},
\qquad -\partial_x^2\phi=\rho_{\rm ext}-\eta.
\]

Here \(\eta\) is the electron-density excess. The periodic source has zero
mean and a finite supplied Fourier expansion. This is an overdamped,
linearized mobility model with electrostatic feedback; it is distinct from
collisionless Vlasov dynamics and from microscopic Coulomb trajectories.
Its equilibrium solves the corresponding screened-Poisson equation.

The oscillating slab uses \(\tau=\omega_{pe}t\) and displacement normalized
by its initial amplitude. Its restoring field follows from planar charge
separation. A uniform shift of particles in a fully periodic homogeneous
three-dimensional plasma has different boundary conditions and is not this
example. The slab derivation is given in
[Fitzpatrick's plasma-frequency lecture](https://farside.ph.utexas.edu/teaching/plasma/lectures1/node6.html).

## Error and derivative interpretation

Connected bands include the declared reconstruction error and a margin over
the entire displayed segment. They are deterministic numerical enclosures
for the named reduced model. They do not establish a bound on a different
physical model, statistical fluctuations or a model-reduction remainder.

The particle Dirac measure has no pointwise distribution band. Particle
moments are weak observables with their own quadrature status. Density and
bulk velocity obtained from the Eulerian distribution have separate bounds;
division requires a strictly positive enclosed density.

The analytic sphere reference is evaluated rigorously. It is an independent
limiting-case reference, not evidence of a general radial boundary-value
solver or of time-dependent finite-particle shielding.

Smooth mathematical leaves retain generated primal/JVP/VJP families.
Discrete loading choices, contact events and display wrapping have explicit
derivative domains. A generated derivative is not automatically a certified
continuous sensitivity. Each example's numerical receipt and documentation
state the available scope.

## Microscopic comparison

The supplementary isolated Coulomb cloud evolves individual electrons in a
fixed positive spherical background, with and without a weak negative central
impurity. It shows actual interactions and an impurity response. Its
trajectories and empirical radial profiles are exploratory candidates;
energy diagnostics do not replace a trajectory certificate. It is not a
converged Debye-shielding or finite-\(N\) fluctuation benchmark.

The matched two-electron pair now adds periodic Coulomb/Ewald mechanics;
many-particle ensembles and the \(\Lambda\)-scan remain the route to
quantitative microscopic shielding. See
[the microscopic many-body contract](microscopic-many-body.md).

## Cost receipts

Measure numerical solve-only time and certificate time separately, retaining
absolute extra time and their ratio with identical inputs and reuse. Plotting
and encoding are additional stages. An unmeasured baseline stays unavailable;
a tiny analytical or scalar example does not establish a general performance
ranking.

## Integrated verification

The quantum/classical extension passes all 21 registered native CPU cases in
Release and checked Debug. All seventeen retained generators reproduce their
leaves byte for byte at FortSym `a6142712af998b727c245fb6cfbb3745b296ec24`.
FortNum `015bb8b1490d99d55493b3936e0dd6adf77bc293` supplies the verified
midpoint integrator; FortPlot `87ff7a43809e868e0153f400a21cdb31b9a92d84`
repairs subplot error-bar placement. Independent audits cover the Hamiltonian,
evolution, Fermi EOS, correlation helper and small cutoff oracle. The final
gallery composition also passes independent review, including raw-norm and
physical-Hamiltonian perturbation bridges, complete plot-cell allowances,
and malformed pair-basis rejection. Both 61-frame additions pass encoded
first/middle/last and grayscale inspection, with all reported bands and
displayed particle projections inside their axes. The render manifest retains
their exact executable, source and artifact identities.

The earlier 2026-10-05 gallery integration passed all 14 registered native CTest cases
in Release and in checked Debug (`-O0 -fcheck=all -fbacktrace`). The six new
registered cases also pass through the verified public Fo driver recorded in
[tooling](tooling.md). The existing resident lane still has its documented
CMake-only input-capture blocker; these are native and focused CLI results,
not a resident green status or formal Lean proof.

Independent audits cover the moment quotient, spherical equilibrium series and
chords, shielding residual/contractivity, oscillator Hermite/Bernstein residual
and Coulomb endpoint enclosure. Review inputs were frozen at base
`564d7edf9d3ab265d0c9aed199983ef370fe50cc`; candidate-only cloud dynamics retain
that status. Each source retains its independent behavioral oracle.

That integration pinned FortSym at `44da68f44c1ef3bc7617a307409c812a6ea4d515`, including the
rigorous large-rational and exact-exponent emission repairs. It pinned FortPlot
at `9949d818b2ba5f43d582d6aff58c8032cdeec3de`, including requested-movie failure
status, retained-raster clearing and procedural subplot-quiver repairs.
Their native upstream regression
suites pass. Final movie acceptance includes encoded first/middle/last frames,
grayscale readability, clipping, legends and full connected band geometry.
That build retained an executable-stack warning for GNU Fortran's internal
animation callbacks. The current build selects heap trampolines when supported
by the compiler, preserving the project profiles.

The render manifest retains executable and source hashes, run resources,
application wall time and GIF conversion time. An application's wall time is
repeated on its copied artifact records; count the corresponding `runs` entry
once rather than summing those copies. For certification overhead use the
numerical performance CSVs, excluding plotting, I/O and encoding.

The quiver repair covers the procedural plotting API used by the gallery.
Direct object-method subplot routing remains adjacent upstream work; it is not
claimed by this consumer check.
