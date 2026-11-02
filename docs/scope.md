# Scope and design principles

This page summarizes the long-term scope of kin6d and the principles that
shape its design. It describes intent, not implemented capability; for what
runs today see the [README](../README.md) and the topic pages listed in the
[documentation index](README.md). The canonical architecture is
[DESIGN.md](../DESIGN.md) and the staged plan is [ROADMAP.md](../ROADMAP.md).

## Aim

kin6d is meant to become one plasma code in which equilibrium, stability,
fluid, kinetic and microscopic models are configurations of shared operators
rather than separate programs. Given a question and a required accuracy, the
user (or eventually an adaptive controller) should be able to choose the
cheapest model and discretization that answers it, with an explicit estimate
of the remaining error.

The deepest classical reference is the relativistic Maxwell–Lorentz
many-particle system; relativistic Vlasov–Maxwell/Fokker–Planck kinetics is
the highest continuum reference. A quantum branch covers degeneracy, exchange
and quantum kinetic models for dense plasmas, but stops short of QED. Fusion
plasmas are the main application; space, low-temperature, beam and laser
plasmas are in the intended envelope where a classical or nonideal model
applies. The full list is in the
[capability catalogue](capability-catalogue.md).

## Principles

[DESIGN.md](../DESIGN.md#non-negotiable-scientific-principles) states these
in full. In short:

- One mathematical definition per operator. Pointwise kernels are written
  once symbolically and generated as primal, JVP and VJP code with FortSym.
- Differentiability from the start. Linear, nonlinear and eigenvalue solves
  are differentiated implicitly at the equation level, not by differentiating
  through iterations. FortAD and finite differences serve as independent
  derivative checks. See [differentiability](differentiability.md).
- Explicit error. A result reports numerical, model-reduction, input and (for
  Monte Carlo) statistical errors separately, or states that a bound is
  unavailable. The fast solver is untrusted; where possible a separate
  interval-arithmetic checker encloses the answer. Lean is optional
  assurance for selected theorems. See [certification](certification.md).
- Structure preservation. Reversible dynamics is written in Hamiltonian or
  Poisson form and dissipation in metric or metriplectic form, with the
  corresponding discrete conservation and positivity properties.
- Stationary states first. Equilibria are solutions of a residual
  `R(U; θ) = 0` of a named model; time evolution uses the same residual.
  Imported equilibria (EQDSK, VMEC and others) are treated as candidates to be
  evaluated or re-solved, and nested flux surfaces are an optional structure,
  not an assumption. See [equilibria and inference](equilibria-inference.md).
- Independent verification and reproducibility. Each case is to be tested
  against analytic, manufactured or external references, and each result
  file records its configuration and provenance. See
  [interface and reproducibility](interface-reproducibility.md).
- Performance at matched accuracy. Speed is compared at the same model,
  observable, accuracy and claim strength. See [performance](performance.md).

## Model and representation axes

The model hierarchy for charged species runs from MHD through extended and
multi-fluid models, gyrofluid and gyromoment hierarchies, to drift, gyro and
full kinetics. Species may use different models in one calculation. The field
model is chosen independently (prescribed, electrostatic, quasineutral
electromagnetic, Darwin, full Maxwell). See [model family](model-family.md),
[species models](species-models.md) and [moment models](moment-models.md).

The numerical representation of a distribution is a separate choice: full-f
or delta-f; grid, semi-Lagrangian, particle or hybrid carrier; PIC, PIF or
particle-in-basis field coupling; deterministic quadrature, quasi-Monte Carlo
or Monte Carlo sampling. See [representations](representations.md).

Analysis type (spectra, forced response, continuation, bifurcation) and
macroscopic transport are further axes; see [stability](stability.md) and
[transport and multiscale coupling](transport-multiscale.md). Adaptive
refinement in space, phase space, model and correlation order with a common
error budget is described in
[adaptive multiscale and uncertainty](adaptive-multiscale-uq.md).

## Planned progression

The kinetic line is steady state first: 1X1V steady problems, then 0X3V
collisions, 1X3V transport, drift kinetics for real devices, gyromoments and
gyrokinetics, and global full kinetics. Small transient and turbulence
examples (Hasegawa–Mima, Hasegawa–Wakatani, interchange) are to be built
early on the same operators as test cases. Equilibrium and linear stability
form a parallel line that currently receives most of the development effort.
Details and ordering are in [ROADMAP.md](../ROADMAP.md).

## Geometry and numerics

Physical vectors are Cartesian, and computational coordinates are general
smooth maps; no core model may require a global toroidal angle or nested
flux surfaces. FortFEM provides the generic spatial discretizations (FEEC,
B-splines, DG/HDG, boundary operators); kin6d provides the plasma physics.
See [geometry](geometry.md), [numerics](numerics.md) and the
[numerical-method catalogue](numerical-methods.md).

Collision models are planned from Krook and Lorentz operators up to
nonlinear and relativistic Landau operators, with the correct adjoint,
nullspace and entropy structure; see [collisions](collisions.md).

## Boundaries

QED processes, general relativity, reactor neutronics, detailed material
mechanics and line-by-line non-LTE radiation transfer are outside the scope
and are left to external codes. Reconstruction and diagnostic inference
belong to the separate TIAGO project, which is to consume kin6d models.

## Software stack

Scientific code and checkers are Fortran; no Python is needed to build, run
or inspect results. The supporting libraries are FortSym, FortNum, FortFEM,
FortSparse, FortIO, FortPlot, FortAD and FortArray, with external solver and
accelerator libraries used behind replaceable interfaces. See
[tooling](tooling.md) and [software engineering](software-engineering.md).
