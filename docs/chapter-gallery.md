# Animation and figure gallery

The gallery accompanies the plasma lecture script in
`/home/ert/proj/plasma-skript`, starting with Chapters 1–3. It lets students
see individual particles, continuous distributions, quantum probability
densities and measurable fields side by side. Every demo runs a declared
model and retains its own claim; the status labels (certified,
certified-discrete, reference comparison, indicator, statistical, candidate)
are defined in
[quantum-microscopic-validation.md](quantum-microscopic-validation.md#status-vocabulary).

## Reproduce the gallery

```sh
cmake --preset cpu
cmake --build --preset cpu
ctest --preset cpu
./scripts/render-gallery.sh            # all demos
./scripts/render-gallery.sh --list     # targets and primary movies
./scripts/render-gallery.sh --only coulomb_encounter_demo,packet_models_demo
```

The applications and scientific postprocessing are Fortran. FortSym emits
the mathematical leaves, FortNum provides rigorous arithmetic and
integrators, FortPlot draws the frames and FFmpeg encodes MP4 and looping GIF
files. Every demo writes `output/<target>.nc` with its literal and resolved
configuration and provenance ([interface-reproducibility.md](interface-reproducibility.md));
the renderer records its sha256 beside every image, movie and CSV. The
style rules are in [figure-style.md](figure-style.md). All generated media
and receipts stay under the ignored `output/` directory. The classical
correlation map is the slowest demo (about 12 minutes); most others take
seconds to a minute.

## Continuum and reduced-model demos

| Target | Movie / figures | Teaching message |
|---|---|---|
| `kin6d_demo` | `free-streaming.mp4`, `final-snapshot.png` | Free streaming shears f(x,v) into ever finer velocity stripes; the certified sup-error bound grows with them and the slice error stays inside |
| `streaming_comparison_demo` | `particles-vlasov.mp4` | Weighted particles and an Eulerian Fourier/P1 grid carry the same free-streaming f; their densities agree within both certified bands |
| `streaming_moments_demo` | `particles-to-moments.mp4` | Velocity moments compress f(x,v) into density and bulk velocity with whole-cell bounds; they follow the flow but discard the velocity structure |
| `debye_sphere_demo` | `debye-potential-reduction.mp4` | As R/lambda_D sweeps from 0.5 to 2, linear electron screening confines the sphere potential within a Debye length (a sequence of equilibria, not a transient) |
| `shielding_dynamics_demo` | `debye-shielding.mp4` | Overdamped electrons pile up on a fixed source and relax exponentially to the linear Debye equilibrium |
| `slab_oscillation_demo` | `plasma-oscillation.mp4` | A displaced cold electron slab oscillates at omega_pe; the certified trajectory is compared with cos(tau) |
| `steady_laboratory_demo` | `fp-steady.png`, `poisson-fp-steady.png`, `steady-observables.png` | Halving the cell width shrinks the certified OU band about fourfold; a driven Poisson--Fokker--Planck plasma settles to a self-consistent steady state |
| `transport_demo` | `general-transport.png`, `general-transport-krook-limit.png` | Refining from 16 to 64 cells tightens the certified bands of f and of its collision-rate derivative; the field-free limit reproduces the Krook value 3/4 |

## Classical microscopic demos

| Target | Movie / figures | Teaching message |
|---|---|---|
| `interaction_demo` | `collective-response.mp4` | Neutral contact acts only on touch, while unscreened Coulomb repulsion deflects the pair smoothly from a distance |
| `classical_pair_demo` | `classical-coulomb-pair.mp4` | In a periodic Coulomb box, repulsion empties small separations of a classical two-electron ensemble, so g(x) dips at x = 0 (simultaneous statistical bars) |
| `microscopic_cloud_demo` | `microscopic-cloud.mp4` | In an exploratory 32-electron cloud, a fixed negative impurity pushes nearby electrons outward |
| `ocp_oscillation_demo` | `ocp-oscillation-landau.mp4`, `-response.png`, `-structure.png`, `-particles.png` | Tens of thousands of real electrons oscillate collectively and phase-mix without a Vlasov equation; the paired response is compared with linear Vlasov and the Landau root |
| `ocp_shielding_demo` | `ocp-shielding.mp4`, `-profiles.png`, `-scan.png`, `-particles.png` | Debye shielding emerges from particles: the averaged shell density around a weak impurity follows the periodic screened-Poisson Green function |
| `correlation_map_classical_demo` | `correlation-map-g2.png`, `-triplet.png`, `-size.png`, `-ybg.png`, `-energy.png` | As Gamma grows the OCP develops a correlation hole and shell structure; three-body correlation is mostly Kirkwood superposition, yet the superposition closure fails above Gamma ~ 10 while HNC stays within 1% |

The `ocp_*` figure names carry the configured `label` as prefix; the
defaults are `ocp-oscillation-landau` and `ocp-shielding`.

## Quantum demos

| Target | Movie / figures | Teaching message |
|---|---|---|
| `quantum_reference_demo` | `quantum-correlations.mp4`, `quantum-correlations.png` | A Coulomb quench of an opposite-spin pair depletes the k = 0 pair by up to about 10% into the first momentum shells and opens a correlation hole in g(x), with certified finite-basis bounds |
| `coulomb_encounter_demo` | `coulomb-encounter-headon.mp4`, `coulomb-encounter-offaxis.mp4` and their `-frame-*.png`, `coulomb-turning-cut.png`, `coulomb-mott-eta-scan.png`, `coulomb-b-scan.png`, `coulomb-kepler.png` | Identical particles interfere: singlet and triplet show exchange fringes and a triplet nodal plane that the distinguishable pair and the classical ensemble lack; outgoing bins follow the Mott oscillation; a Rydberg packet breathes, collapses and revives |
| `packet_models_demo` | `packet-density.mp4`, `packet-density-mid.png`, `packet-density-final.png`, `packet-error-vs-eta.png`, `packet-angular-eta1-eta10.png` | One Gaussian stays a blob while the exact density forms a fringed crescent; Gaussian models lose accuracy as eta grows while the classical ensemble improves |
| `fermion_ensemble_demo` | `fermion-occupations.png`, `fermion-exchange-hole.png`, `fermion-eos-scan.png` | Free same-spin fermions avoid each other within k_F r < 3; canonical box occupations lie on Fermi--Dirac and the finite-N EOS converges to the thermodynamic limit |
| `relativistic_fermi_demo` | `relativistic-fermi-eos.png` | Electron pressure at 1e7 K crosses from n k_B T through the nonrelativistic Fermi law to the softer ultrarelativistic law above x_F ~ 1 |
| `jellium_cell_demo` | `jellium-kick.mp4`, `jellium-pair-correlation.png`, `jellium-response.png`, `jellium-cumulants.png`, `jellium-mode-history.png` | The opposite-spin Coulomb hole deepens with r_s while the exchange hole barely changes; RPA overscreens; a kicked finite cell rings without Landau damping |
| `mean_field_cell_demo` | `mean-field-correlation.png`, `mean-field-pair-correlation.png`, `mean-field-response.png`, `mean-field-estimator.png` | HF has the exchange hole but no Coulomb hole; TDHF tracks the exact N = 2 response far better than RPA; open-shell UHF breaks translation symmetry |
| `screening_demo` | `screening-hierarchy.png` | Static screening goes from Debye (hot) to Thomas--Fermi (degenerate); Lindhard adds Friedel oscillations of period pi/k_F that heating smears out |
| `wigner_demo` | `wigner-two-stream.mp4`, `wigner-two-stream.png`, `wigner-dispersion.png` | Wigner--Poisson reproduces quantum Langmuir and two-stream roots; hbar stabilizes short wavelengths, f_W converges to Vlasov as hbar^2 and goes negative where the beams overlap |
| `hybrid_encounter_demo` | `hybrid-beam.mp4`, `hybrid-beam-frame*.png`, `hybrid-coulomb-log.png`, `hybrid-relaxation-rate.png` | The quantum Coulomb logarithm saturates while the classical one keeps growing; promoting close encounters to quantum kernels slows a beam less, matching the all-quantum rate within 4e-3 |

## Model notes

The bounded-velocity distribution in the moments example is specified by its
finite Fourier/P1 coefficients. It differs from the script's illustrative
unbounded Maxwellian-plus-beam profile, and no omitted velocity tail is
silently included in its certificate.

The shielding relaxation solves, in normalized units,

\[
\partial_t\eta=\partial_x^2\eta-\eta+\rho_{\rm ext},
\qquad -\partial_x^2\phi=\rho_{\rm ext}-\eta.
\]

Here \(\eta\) is the electron-density excess. The periodic source has zero
mean and a finite supplied Fourier expansion. This is an overdamped,
linearized mobility model with electrostatic feedback; it is distinct from
collisionless Vlasov dynamics and from the microscopic shielding demo, whose
screening emerges from Coulomb trajectories.

The oscillating slab uses \(\tau=\omega_{pe}t\) and displacement normalized
by its initial amplitude. Its restoring field follows from planar charge
separation. A uniform shift of particles in a fully periodic homogeneous
three-dimensional plasma has different boundary conditions and is not this
example; the periodic OCP demo uses a finite-k displacement instead. The
slab derivation is given in
[Fitzpatrick's plasma-frequency lecture](https://farside.ph.utexas.edu/teaching/plasma/lectures1/node6.html).

The pair quench ([quantum-reference-gallery.md](quantum-reference-gallery.md),
[classical-pair.md](classical-pair.md)) is a microscopic cold-quench example,
not a hot bulk-degeneracy or Debye-screening simulation. Its quantum bands
cover the declared finite basis; its classical intervals describe sampling
uncertainty.

## Error and derivative interpretation

Connected bands include the declared reconstruction error and a margin over
the entire displayed segment. They are deterministic numerical enclosures
for the named reduced model. They do not establish a bound on a different
physical model, statistical fluctuations or a model-reduction remainder.
Statistical bands in the microscopic demos are styled differently and state
their sampling law.

The particle Dirac measure has no pointwise distribution band. Particle
moments are weak observables with their own quadrature status. Density and
bulk velocity obtained from the Eulerian distribution have separate bounds;
division requires a strictly positive enclosed density.

Smooth mathematical leaves retain generated primal/JVP/VJP families.
Discrete loading choices, contact events, spin sectors and display wrapping
have explicit derivative domains. A generated derivative is not
automatically a certified continuous sensitivity. Each demo's document
states the available scope.

## Cost receipts

Measure numerical solve-only time and certificate time separately, retaining
absolute extra time and their ratio with identical inputs and reuse. Plotting
and encoding are additional stages. An unmeasured baseline stays unavailable;
a small analytical or scalar example does not establish a general
performance ranking.

The render manifest retains executable and source hashes, run resources,
application wall time and GIF conversion time. An application's wall time is
repeated on its copied artifact records; count the corresponding `runs` entry
once rather than summing those copies. For certification overhead use the
numerical performance CSVs, excluding plotting, I/O and encoding.

## Integrated verification

The current integration (hierarchy campaign, 2026-10-06) passes 38/38
registered CTest cases in Release and Debug, and all retained generators
reproduce their committed leaves at FortSym
`a6142712af998b727c245fb6cfbb3745b296ec24`. It pins FortNum
`e97a493bfe67a404daa2d7c2090ff3eef9c32b75`, FortPlot
`cad62c675002029ab7a0d875d80baeb155fca8d1` and FortIO
`38add03dc691a7cd4715cd3e611ad42422e46091`. Each demo's own oracles are listed in
its document; earlier integration receipts (pins, reviews, frame
inspections) are in the git history of this file and in
[hierarchy-campaign.md](hierarchy-campaign.md). The renderer runs no tests:
pass a controller receipt with `--verification-receipt`.

Known FortPlot limitations visible in the current figures: lines are not
clipped at the axes, `text()` annotations are not rendered, log ranges are
clamped to 1e±25, and some TeX glyphs (\langle, \hbar, arrows) are missing.
