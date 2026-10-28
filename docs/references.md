# References and validation scope

Checked on 2026-10-05. These sources motivate individual components and
candidate comparisons. They do not establish novelty of the complete project
or a certification theorem for its final nonlinear kinetic model.

## Numerical and proof tools

- [PETSc Landau collision operator](https://petsc.org/release/manualpages/LANDAU/)
  documents its velocity-space discretization and residual/Jacobian support.
  The [multispecies tutorial](https://petsc.org/release/src/ts/utils/dmplexlandau/tutorials/ex1.c.html)
  includes species/grid configuration and a relativistic-correction option.
  Pin a version and identify the exact equations before using an option as
  an independent physical benchmark. Different discretizations strengthen
  numerical comparisons; they do not supply a rigorous continuous enclosure.
- [PETSc SNES](https://petsc.org/release/manual/snes/) and
  [KSP](https://petsc.org/release/manual/ksp/) are the proposed numerical
  nonlinear and Krylov solvers. Their stopping criteria are candidate-solve
  diagnostics, not proof of a continuous solution.
- [RadiiPolynomial.jl](https://olivierhnt.github.io/RadiiPolynomial.jl/stable/)
  illustrates validated numerical constructions for computer-assisted proofs.
  An infinite-dimensional application still needs model-specific space,
  inverse, complement, and derivative estimates.
- [LeanCert](https://github.com/alerad/leancert) provides checked interval,
  integration, and root-bound tools in Lean. Evaluate the smallest required
  operations for R0 and pin the version. No general Vlasov–Maxwell–collision
  theorem is inferred from this tool's capabilities.

## Kinetic representations and collision physics

- [OSHUN documentation](https://oshun.readthedocs.io/en/latest/) describes
  spherical-harmonic momentum expansion, relativistic Vlasov treatment, and
  a nonrelativistic linearized Fokker–Planck treatment in that documented
  version. It is a representation and limiting-case reference.
- [VlaPy](https://github.com/joglekara/VlaPy) implements 1X1V
  Vlasov–Poisson–Fokker–Planck evolution and analytical regression cases.
  Its periodic transport and particular collision model must overlap a
  kin6d case before comparison. It is not a certified steady-state oracle.
- Braams and Karney,
  [Conductivity of a Relativistic Plasma](https://arxiv.org/abs/plasm-ph/9502001),
  reformulate relativistic collisions in spherical harmonics with five scalar
  potentials and one-dimensional integral formulas. This is the main basis
  motivation for early Beliaev–Budker work. Its conductivity application does
  not establish a final arbitrary multispecies implementation.
- Daniel, Taitano, and Chacón,
  [A fully implicit, scalable, conservative nonlinear relativistic Fokker–Planck 0D-2P solver](https://arxiv.org/abs/1902.10241),
  demonstrate nonlinear relativistic potential solvers and conservation/
  positivity methods. Their electron–ion treatment uses a Lorentz operator;
  this is not an arbitrary nonlinear multispecies BB implementation.
- Duerinckx and Winter,
  [Well-posedness of the Lenard–Balescu equation with smooth interactions](https://arxiv.org/abs/2111.13320),
  analyze collective screening under specific smooth-interaction assumptions.
  This theory is not transferred to singular Coulomb, relativistic,
  magnetized, bounded tokamak problems without new analysis.
- Gould and DeWitt,
  [Convergent Kinetic Equation for a Classical Plasma](https://doi.org/10.1103/PhysRev.155.68),
  give a weak-coupling construction combining complementary contributions.
  A production matched operator must derive its actual screening and overlap
  terms; the shorthand `Boltzmann + LB - overlap` is not that derivation.
- Baalrud and Daligault,
  [Mean Force Kinetic Theory](https://arxiv.org/abs/1904.09208),
  motivate a correlation-based closure spanning different coupling regimes.
  The required equilibrium pair information and nonequilibrium validity
  remain part of any proposed benchmark's assumptions.

## Structure and transient methods

- Kraus, Kormann, Morrison, and Sonnendrücker,
  [GEMPIC](https://arxiv.org/abs/1609.03053),
  develop geometric Vlasov–Maxwell particle-in-cell discretization. The
  compatible field complex is relevant; its particle representation differs
  from kin6d's deterministic distribution coefficients.
- Kraus and Hirvijoki,
  [Metriplectic Integrators for the Landau Collision Operator](https://arxiv.org/abs/1707.01801),
  derive discrete conservation and entropy properties for a collision
  discretization. These properties motivate future time integration and
  invariants; they are distinct from a quantitative perturbation-stability
  or long-time trajectory enclosure theorem.

## Fortran ecosystem snapshot

This snapshot records earlier source inspection. Current integration pins and
executed gates are recorded in [the steady laboratory](steady-laboratory.md).
GitHub main revisions were checked through the repository API. fortnum,
fortsym, and fo source were also read locally at matching revisions. This
is source inspection, not a claim that kin6d integration tests have run.

| Repository | Checked main revision |
|---|---|
| fortplot | [5daec546](https://github.com/lazy-fortran/fortplot/tree/5daec54649e305fd09bc1ef8345511d5692a1889) |
| fortnum | [1a52e13c](https://github.com/lazy-fortran/fortnum/tree/1a52e13c5e28b09d083947a6ec7e8a1dcf9d502e) |
| fortio | [38add03d](https://github.com/lazy-fortran/fortio/tree/38add03dc691a7cd4715cd3e611ad42422e46091) |
| fortad | [c7e01cc9](https://github.com/lazy-fortran/fortad/tree/c7e01cc9f00efd288502ba003270bcf00b42f6f0) |
| fortsym | [a66f3756](https://github.com/lazy-fortran/fortsym/tree/a66f3756585752fb3e3aa9cf4d27646919656734) |
| fortarray | [e121e314](https://github.com/lazy-fortran/fortarray/tree/e121e31460624b30bb1eeba7fa0c0a0d5bbd5188) |
| fo | [1c45a559](https://github.com/lazy-fortran/fo/tree/1c45a5595fd1e9c0a05b2d4a71650826efd857e9) |

The [fortsym rigorous-runtime contract](https://github.com/lazy-fortran/fortsym/blob/a66f3756585752fb3e3aa9cf4d27646919656734/doc/rigorous-runtime.md)
defines supported expressions, exact literals, runtime requirements, and
refusals. It distinguishes paired Fortran float/enclosure emission from
the CUDA numerical target. Unsupported transcendental enclosure emission and
AD-to-DAG integration must be treated as generic upstream work when needed.

The [fortnum verified sources](https://github.com/lazy-fortran/fortnum/tree/1a52e13c5e28b09d083947a6ec7e8a1dcf9d502e/src/verified)
provide interval/ball arithmetic, rigorous radix-2 transforms, and verified
linear algebra with documented arithmetic assumptions. Lean soundness bridges
are optional additional assurance; the production containment argument follows
the Fortran verified-arithmetic contract in [certification](certification.md).

The [fo CMake plan](https://github.com/lazy-fortran/fo/blob/1c45a5595fd1e9c0a05b2d4a71650826efd857e9/doc/CMAKE_GREMLIN_PLAN.md)
distinguishes present delegation to CMake/CTest from future standalone native
compatibility. The actual public help was read for command semantics.
[Tooling policy](tooling.md) defines their use and the upstream-main workflow.

## Benchmark selection

GENE-X's [collision module](https://phoenix-public.pages.mpcdf.de/genex/module/op_coll_m.html)
documents collision-base and concrete BGK/LBD/Lorentz/FPL types. Ulbl,
Michels, and Jenko's [LBD paper](https://onlinelibrary.wiley.com/doi/full/10.1002/ctpp.202100180)
supports the discrete mixing-quantity conservation approach. The
[improved Sugama model](https://arxiv.org/abs/1906.07427) supplies a reduced
linearized comparison reference. [Collision specification](collisions.md)
uses these within their assumptions rather than claiming a final full-orbit
certification theorem from them.

fortplot's [fill-between example](https://github.com/lazy-fortran/fortplot/tree/5daec54649e305fd09bc1ef8345511d5692a1889/example/fortran/fill_between_demo)
and error-bar example provide plotting facilities for interval displays.
Their numerical enclosure and interpolation semantics must still be checked
by the kin6d postprocessor contract; a rendering routine supplies no error
theorem.

NEO-2, CQL3D, COGENT, and KAMEL are proposed limiting-case references from the
project brief. Each comparison needs a written mapping of geometry, collision
model, normalization, boundary data, asymptotic ordering, and observables.
Do not interpret disagreement outside that intersection as a solver defect.
Select and archive a particular executable case before treating it as a gate.

The project brief's broader novelty claims, cost multipliers for collision
models, general tokamak closure remainders, and a complete Lean formalization
of the final PDE remain unestablished. Prior-art and analytical work should
answer those questions when their milestones require them.

## Species reductions and device missions

[Brizard and Tronci's guiding-center Vlasov–Maxwell formulation](https://arxiv.org/abs/1602.05030)
derives polarization/magnetization and conservation from a variational
formulation. It motivates consistent reduced field coupling; arbitrary
mixed FK/GK/DK collision and field combinations remain project obligations.
[ORB5's model description](https://arxiv.org/abs/1905.01906) provides an
existing gyrokinetic/kinetic-electron model precedent, not the full kin6d
per-species feature envelope or continuous certificates.

[The AP semi-Lagrangian Vlasov–Maxwell paper](https://www.sciencedirect.com/science/article/pii/S0021999125001238)
is a numerical precedent for the quasineutral limit. Its stability and
conservation results do not establish kin6d's continuous-error bounds or
resolve arbitrary fast solutions on coarse time/space representations.

[Device envelope](devices.md#device-envelope) links official experimental
descriptions and primary reactor studies next to their parameter claims.
The SPARC overview source is DOI 10.1017/S0022377820001257. ARC uses
arXiv:1409.3540v2; Helios uses arXiv:2512.08027v3. The cited 2025 Stellaris
study specifies four periods and a 12.7 m major radius. These are named
design points, not timeless specifications or measured reactor operation.
No equilibrium input or production benchmark has yet been selected.

## Analytical and bulk-model references

- [Landreman 2026, arXiv:2609.26742](https://arxiv.org/abs/2609.26742) gives
  analytic, non-axisymmetric toroidal MHD equilibria (steady Euler flows)
  with exact nested surfaces. They are counterexamples to Grad's conjecture
  and exact 3D benchmarks. HTML v2 reviewed 2026-10-12; its ε = 0 equations (7)–(8)
  are now re-derived and checked symbolically by `gen/gen_landreman.f90`
  (2026-10-20). The [free-boundary benchmark](ideal-free-boundary.md#landreman-axisymmetric-benchmark)
  supplies numerical evidence conditional on G4; the full 3D families are not
  re-derived by that generator.
- [Gómez-Serrano, Liehr and Taylor 2026, arXiv:2609.24739](https://arxiv.org/abs/2609.24739)
  give a Lean-certified family of smooth counterexamples to Grad's conjecture
  with cyclic symmetry. Abstract only; not yet read in full.
- [Issan, Citrin, Landreman and Tracey 2026, arXiv:2610.07304](https://arxiv.org/abs/2610.07304)
  (v1, 5 Oct 2026) generalise both Landreman families with a symmetry-breaking
  parameter (no stellarator symmetry). Abstract only.

The [Analytical Benchmark Atlas](benchmarks.md) identifies closed forms,
semi-analytic evaluation, asymptotes, identities, and numerical comparisons.
Guazzotto–Freidberg exact Grad–Shafranov geometry uses
[DOI 10.1017/S002237782100009X](https://doi.org/10.1017/S002237782100009X).
The [Landau time-domain construction](https://arxiv.org/abs/2503.22501),
[Bernstein time-domain construction](https://arxiv.org/abs/2608.12994), and
[cylindrical pinch constructions](https://arxiv.org/abs/2608.06780) are specific
semi-analytic candidates requiring rigorous evaluation before certification.
Their existence does not establish all atlas cases or their error theorems.

[Model-family references](model-family.md#independent-references-and-certification)
define the checked local ARES scope and documented GPEC/JOREK intersections.
The [VMEC++ documentation](https://proximafusion.github.io/vmecpp/) describes
an external ideal-MHD equilibrium candidate source, not a rigorous continuous
equilibrium enclosure. External solvers retain their own model/topology scope.

The user's expanded first-ring literature list and curated second ring belong
to the external corpus under the [literature policy](literature-policy.md).
Metadata resolution, lawful OA provenance, attachment hashes, review queues
and Zotero imports are outside the solver repository. DOI metadata and a PDF
download do not substitute for reading and checking a benchmark's assumptions.

## Three-dimensional equilibrium comparators

Checked 2026-10-12; "abstract" means the arXiv/Crossref abstract or metadata
only.

- [Blickhan, Stratton and Kaptanoglu, arXiv:2510.26986](https://arxiv.org/abs/2510.26986),
  "MRX: A differentiable 3D MHD equilibrium solver without nested flux
  surfaces" (v1 30 Oct 2025, v3 28 Sep 2026). Full text of v3 read
  2026-10-28. Open source: [github.com/ToBlick/mrx](https://github.com/ToBlick/mrx)
  (MIT, JAX). Fixed boundary only (free boundary is future work). B-spline
  FEEC on a fixed map from VMEC/GVEC with C¹ polar splines; div B = 0 to
  machine precision; harmonic forms deflated from every solve. Ideal
  relaxation δB = curl(v × B) with divergence-free v; pressure is a Lagrange
  multiplier (Leray projection or the potential-velocity form
  curl (curl curl)⁻¹ curl(J × B)), so the final profile depends on the initial
  condition. Inexact Newton-MR with MINRES on the second variation, a
  Levenberg–Marquardt term on the parallel velocity only (κ = 3), matrix-free
  sum-factorized operators, and preconditioners from metric lumping (Kronecker
  products/sums of 1-D matrices plus a dense axis block) and geometric
  multigrid on the div-div stiffness (in the code; the paper defers the
  description). Island seeding at rational surfaces with amplitudes from a
  quadratic energy criterion; resistive phase guided by a prescribed current
  carries states of different topology to one equilibrium. NCSX li383
  (β ≈ 4.25%): normalized squared force residual 1e-8 in 5.4 min at 1.8e5
  degrees of freedom on one H100; floor 4e-11 to 1.5e-10 that does not
  decrease under h or p refinement (current sheets at rational surfaces).
  No rotating-ellipse case. VMEC, GVEC and DESC are stated to assume nested
  flux surfaces.
- [Zhang, Zhu and Hegna, arXiv:2412.18922](https://arxiv.org/abs/2412.18922)
  (25 Dec 2024), published as Computer Physics Communications 320, 109961
  (March 2026), [DOI 10.1016/j.cpc.2025.109961](https://doi.org/10.1016/j.cpc.2025.109961).
  Full text read 2026-10-28. Code NTEC (C++, Eigen, Intel MKL PARDISO); no
  code or data availability is stated, so it enters only through published
  results. Quintic Bell (C¹) triangles in the poloidal plane times toroidal
  Fourier modes on a pseudo-flux map from VMEC; div B kept by a cleaning
  term. Pressure relaxes towards B·∇p = 0 by a hyperbolic (damped) parallel
  diffusion equation; B by resistive induction with a viscous flow; linear
  terms implicit, nonlinear explicit, thousands of iterations. Reports
  islands and stochastic regions: fixed-boundary classical stellarator
  against SIESTA (boundary given analytically, eq. 3-1), free-boundary CFQS
  against published HINT results at ⟨β⟩ ≈ 1.4–1.6%, helical-core tokamak and
  RFP quasi-single-helicity states.
- [Balkovic et al., arXiv:2607.27135](https://arxiv.org/abs/2607.27135),
  "SPECTRE: A robust solver for 3D equilibria with arbitrary topology"
  (29 Jul 2026). Abstract and section list read: MRxMHD, fixed and free
  boundary, QA vacuum, free-boundary W7-X vacuum, classical-stellarator
  finite-β free boundary, rotating-ellipse free boundary, finite-β QI with a
  core island compared with HINT.
- [Hudson et al., Phys. Plasmas 19, 112502 (2012)](https://doi.org/10.1063/1.4765691),
  "Computation of multi-region relaxed magnetohydrodynamic equilibria"
  (stepped-pressure SPEC; abstract).
- [Cooper et al., Comput. Phys. Commun. 180, 1524 (2009)](https://doi.org/10.1016/j.cpc.2009.04.006),
  "Three-dimensional anisotropic pressure free boundary equilibria" (ANIMEC;
  Crossref metadata only).
- HINT/HINT2 and PIES: Hayashi et al., Contrib. Plasma Phys. 42, 309 (2002),
  "Island formation for finite pressure equilibria in non-axisymmetric tori
  obtained by the HINT code" (title); Suzuki et al., Nucl. Fusion 46, L19
  (2006), HINT2 (abstract); "Pressure-induced breaking of equilibrium flux surfaces in the W7AS
  stellarator", Nucl. Fusion (2007), chaotic regions (abstract; the code is
  not named there); Salas et al. 2004, TJ-II with PIES (title). Island and stochastic
  capability is taken from these records only.

## Cartesian maps, edge data and broader scope

[Plunk et al., Figure-8](https://arxiv.org/abs/2411.16411),
[Hindenlang–Plunk–Maj, flexible frame](https://arxiv.org/abs/2410.17595), and
[Babin et al., invertible mappings](https://arxiv.org/abs/2411.04683) motivate
the selectable Cartesian embedding. They supply candidate constructions;
kin6d still needs rigorous map and transformed-operator bounds.

[OPEN-ADAS ADF11](https://open.adas.ac.uk/adf11) documents derived coefficient
classes, while [ADF04](https://open.adas.ac.uk/adf04) documents resolved
level/rate data whose coverage must be inspected per file. These do not
provide arbitrary differential product cross sections by default.
[EIRENE source/manuals](https://github.com/iterorganization/EIRENE) provide
a neutral backend/reference; sampling uncertainty is statistical.
[GENE-X](https://genex.ipp.mpg.de/about) and
[XGC](https://xgc.pppl.gov/html/general_info.html) supply kinetic edge
precedents, not certificates for unrestricted FK plasma–neutral–wall equations.

The broader application and actuator/diagnostic tables are planned scientific
requirements, not assertions of code parity, prior-art completeness or
existing proof coverage. Pin and review the exact primary formula/data/code
source before activating each benchmark.

## FortFEM source snapshot

[Spatial numerics](numerics.md) pins FortFEM main
`2d7baa19fe1bef02ce279fa407830c8731bdc663`, its scope and spline source files.
The inspection confirms mathematical ingredients and an assembly-oriented
spline path; it is not a build/performance audit or a continuous certificate.
Sum factorization, partial assembly and batching are client-driven production
goals to add upstream when the next workload requires them. FortFEM retains
generic spatial ownership; plasma closures remain outside it.

## Performance infrastructure

[FortNum benchmarks](https://github.com/lazy-fortran/fortnum/blob/main/benchmark/README.md)
document native measurements, controlled reference evidence and optional
Python gates. [FortAD vector-JVP](https://github.com/lazy-fortran/fortad/blob/main/docs/design/vector-jvp-benchmark.md)
and [FortIO performance](https://github.com/lazy-fortran/fortio/blob/main/docs/performance.md)
are upstream starting points; pin and verify the selected client workflow.
No mandatory Python gate is imported into kin6d and no kin6d frontier
measurement is claimed. [Performance policy](performance.md) defines that
future evidence and promotion contract.

Chris&AI
