# kin6d documentation

Pages are grouped by topic. Pages that describe implemented code name the
modules, applications and tests involved; design and scope pages say so in
their opening status line. Start with the [README](../README.md) for an
overview and [scope and design principles](scope.md) for the long-term aims.

## Project-level documents

- [Scope and design principles](scope.md)
- [DESIGN.md](../DESIGN.md): canonical architecture
- [ROADMAP.md](../ROADMAP.md): staged delivery plan
- [Software engineering policy](software-engineering.md),
  [tooling and development workflow](tooling.md),
  [figure style](figure-style.md), [literature policy](literature-policy.md)
- [Configuration, I/O and reproducibility](interface-reproducibility.md)

## Kinetic laboratory (1D–1V)

- [First 1X1V laboratory](first-example.md)
- [Steady Poisson and Fokker–Planck laboratory](steady-laboratory.md)
- [General damped transport](general-transport.md) and its
  [certificate argument](transport-certificate.md)
- [Particles and Eulerian free streaming](streaming-comparison.md),
  [particles, distribution and moments](streaming-moments.md)
- [Finite spherical source and Debye equilibrium](debye-sphere.md),
  [planar shielding relaxation](shielding-dynamics.md),
  [cold electron slab](slab-oscillation.md),
  [electromagnetic wave at a plasma edge](em-slab-reflection.md)
- [Student laboratory](student-laboratory.md),
  [runnable demonstrations](demonstrations.md),
  [animation and figure gallery](chapter-gallery.md)

## Axisymmetric equilibrium

- [Fixed-boundary Grad–Shafranov](grad-shafranov.md)
- [Prescribed-q inverse equilibrium](gs-inverse.md)
- [Free-boundary Grad–Shafranov](gs-free-boundary.md)
- [Native P2/P3 coarea evaluation](gs-coarea.md)

## Three-dimensional equilibrium

- [Common equilibrium core](equilibrium-core.md) (design) and
  [verification harness](equilibrium-harness.md)
- [Ideal frozen-in map equilibrium](ideal-map-3d.md),
  [boundary and profile initializer](ideal-seed-3d.md),
  [GLL map experiment](ideal-map-gll.md)
- Free boundary: [formulation and proof obligations](free-boundary-theory.md),
  [map-lane free boundary](ideal-free-boundary.md),
  [circuit currents as parameters](free-boundary-nested-map-3d.md),
  [chart-backed coupling surface](free-boundary-chart-surface.md)
- [Exact nested 3D references (Landreman, Issan et al.)](nested-oracles.md)

## Non-nested equilibrium and islands

- [Equilibria with islands](island-equilibrium.md) (design)
- [NN1: vacuum and Solov'ev on the FEEC complex](equilibrium-nonnested.md),
  [mapped 3D domains](nn1-mapped-3d.md),
  [quasi-axisymmetric vacuum benchmark](nn1-qa-benchmark.md)
- [NN2: resistive relaxed lane](equilibrium-nonnested-nn2.md),
  [frozen fitted profiles](nonnested-profiles.md)
- [Taylor (Beltrami) closure](equilibrium-beltrami.md),
  [metriplectic relaxation](metriplectic-relaxation.md)
- [Relaxation comparators and adopted algorithms](relaxation-comparators.md) (design)
- [G6 sample readback](g6-readback.md),
  [external sampled force diagnostic](g6-sampled-force.md)

## Stability and response

- [Stability, response and bifurcation](stability.md) (specification)
- [Ideal displacement energy](ideal-xi.md),
  [S0 physical-mass spectrum](s0-spectrum.md)
- [Linear per-n MHD response](mhd-linear-response.md),
  [unified response operator family](response-unified.md) (design)
- [Strip jets and preconditioning](response-strip-jets.md),
  [response calibration and qualification](response-qualification.md)
- [Drift kinetics: neoclassical transport and NTV](drift-kinetics.md) (design)

## Microscopic classical and quantum references

- [Microscopic many-body reference](microscopic-many-body.md) and
  [hierarchy campaign](hierarchy-campaign.md)
- Classical: [two-body interaction laboratory](interaction-laboratory.md),
  [matched classical pair quench](classical-pair.md),
  [isolated Coulomb cloud](microscopic-cloud.md),
  [periodic OCP oscillation and shielding](ocp-microscopic.md),
  [OCP correlation hierarchy](correlation-hierarchy.md)
- Quantum: [validation plan](quantum-microscopic-validation.md),
  [periodic two-electron quench](quantum-reference-gallery.md),
  [two-body Coulomb encounter](coulomb-encounter.md),
  [wave-packet models](packet-models.md),
  [finite-basis evolution](quantum-evolution.md),
  [fermionic Fock reference](quantum-fock.md),
  [fermionic correlations](quantum-correlations.md),
  [quantum versus classical correlations](quantum-classical-correlations.md),
  [ideal fermion ensembles](fermion-ensemble.md),
  [ideal Fermi gas](fermi-gas.md),
  [jellium cell](jellium-cell.md),
  [mean-field models](mean-field-cell.md),
  [screening and Wigner to Vlasov](screening-wigner.md),
  [hybrid quantum encounters](hybrid-encounter.md)
- Proton–proton fusion: [two-body potential](pp-fusion.md),
  [pionless EFT](pp-fusion-eft.md)

## Numerics, derivatives and error control

- [Numerical foundation](numerics.md),
  [numerical-method catalogue](numerical-methods.md),
  [geometry and charts](geometry.md)
- [Differentiability](differentiability.md),
  [adjoints and structural certificates](adjoints-structure.md)
- [Certification](certification.md),
  [error control and performance](error-control-strategy.md),
  [performance engineering](performance.md)
- [Analytical benchmark atlas](benchmarks.md),
  [theory map](theory-map.md), [references](references.md)

## Long-range scope (design only)

- [Capability catalogue](capability-catalogue.md),
  [model family](model-family.md), [species models](species-models.md),
  [moment models](moment-models.md), [representations](representations.md),
  [collisions](collisions.md), [coupling of background and response](coupling.md)
- [Transport and multiscale coupling](transport-multiscale.md),
  [adaptive multiscale and uncertainty](adaptive-multiscale-uq.md),
  [time-dependent extension](transients.md)
- [Equilibria, imported candidates and inference](equilibria-inference.md),
  [real-device applicability](devices.md)
- [Edge physics](edge-physics.md),
  [reacting and partially ionized plasmas](reacting-partially-ionized.md),
  [radiation](radiation.md), [quantum plasma](quantum-plasma.md),
  [space plasmas](space-plasmas.md), [extensions](extensions.md)
- [Programme estimates](programme-estimates.md)
