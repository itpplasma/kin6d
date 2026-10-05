# Moment models: fluid, gyrofluid, Landau-fluid and gyromoments

Status: design specification. No production gyrofluid/gyromoment branch is yet
implemented in KIN6D.

## Role in the hierarchy

Moment models are not separate lightweight codes. They are controlled reduced
states sitting between bulk fluid descriptions and kinetic distributions.

The intended hierarchy is

\[
{\rm MHD/fluid}
\rightarrow
{\rm anisotropic\ fluid}
\rightarrow
{\rm gyrofluid/Landau\ fluid}
\rightarrow
{\rm gyromoments}
\rightarrow
{\rm GK}
\rightarrow
{\rm FK}.
\]

The first useful turbulence demonstrations are low-dimensional, but their
operators and free-energy structure must be compatible with the later hierarchy.

## Parent structure

Where available, derive reduced moments from a variational gyrokinetic or
kinetic parent rather than guessing closures term by term. A typical
gyrokinetic expansion is

\[
F_s(R,v_\parallel,\mu)
=
\sum_{p,j} N_s^{pj}(R)\,
H_p(v_\parallel)\,L_j(\mu B),
\]

or an explicitly declared equivalent basis.

The kinetic equation generates an infinite hierarchy for \(N_s^{pj}\).
Selecting a finite index set plus a closure gives a moment model.

Very small retained sets give gyrofluid equations. Larger sets produce a
gyromoment representation. A closure may represent unresolved phase mixing,
collisions or both, but these physical effects must not be conflated.

## Fluid and anisotropic moment branches

Ordinary fluid options include density, momentum, scalar pressure and selected
heat fluxes. Space-plasma and strongly magnetized applications also need
anisotropic pressure states such as

\[
p_\parallel,\qquad p_\perp
\]

and CGL/six/ten-moment closures.

A ten-moment pressure tensor does not close itself; the third moment or a
nonlocal closure remains a physical model choice.

Braginskii, CGL and Landau-fluid closures are not a monotone accuracy ladder.
Their assumptions differ.

## Gyrofluid

Gyrofluid models should retain a declared subset of gyrocenter moments together
with the consistent polarization, magnetization, FLR and field couplings of
their ordering.

The defining structural requirement is a correct free-energy/energy-transfer
ledger. Scott's GEM work is an important reference for deriving coupled terms
from free-energy consistency. Kendl's TIFF/T3FF are useful full-f/full-k
gyrofluid comparison cases when their assumptions match.

"Full-k" or "full-f" in a reference model does not mean full kinetic; it names
particular retained nonlinear/FLR amplitudes within a gyrofluid ordering.

## Landau-fluid closure

Hammett--Perkins-type closures approximate collisionless phase mixing and
Landau response through rational/nonlocal moment closures. They are especially
valuable because they bridge cheap fluid dynamics and kinetic linear response.

A Landau-fluid closure must declare:

- the response function/rational approximation;
- parallel geometry and boundary realization;
- retained moments;
- collision treatment;
- free-energy/damping interpretation.

Landau damping represented by closure is not Coulomb entropy production.
On open SOL field lines, a periodic Fourier closure is not automatically valid.

## Gyromoment hierarchy

The gyromoment branch should make the retained index set explicit and
configurable. Increasing moments is useful only if the numerics also control:

- streaming recurrence/reflection from the truncation;
- collision coupling;
- phase-space measure;
- nonlinear products;
- moment-tail reconstruction.

The hierarchy is valuable because the same calculation can be repeated at
increasing moment count while holding geometry and physical inputs fixed.

A successful finite hierarchy is a certified solution of that finite model.
Calling it GK requires a separate tail/closure remainder argument or convergence
evidence with the correct status.

## Hamiltonian, metric and open-system structure

For a closed reversible moment system, seek

\[
\dot U=J(U)\nabla H(U),\qquad J^\dagger=-J.
\]

For collisional/dissipative pieces,

\[
\dot U=M(U)\nabla S(U)
\]

with the correct symmetry, positivity and degeneracies for the selected model.

Sources, profile maintenance, sheath losses, curvature drive and other open
terms are ports. Their energy/free-energy contribution is shown explicitly.

A structure claim requires more than a visually conserved scalar. If a Poisson
bracket is claimed, Jacobi or the appropriate reduced structure must be
established for the semi-discrete form.

## Numerical machinery

The lowest-dimensional slab models may use Fourier actions directly. General
geometry should use the common FortFEM/FCI infrastructure.

Reusable primitives include:

- perpendicular Laplacian/Helmholtz solves;
- Poisson brackets;
- parallel derivative/FCI actions;
- gyroaverage/FLR operators;
- moment coupling matrices;
- closures;
- electromagnetic field blocks.

All continuous project-owned local leaves are FortSym-generated with
primal/JVP/VJP.

The same matrix-free/tensor performance rules as the kinetic code apply. Small
moment count should produce compact dense/tensor blocks rather than forcing a
general sparse representation.

## Early turbulence ladder

The first demonstration sequence is:

1. Hasegawa--Mima;
2. Hasegawa--Wakatani;
3. interchange/blob;
4. minimal ITG gyrofluid;
5. electromagnetic drift-Alfvén;
6. a Scott/Kendl-compatible gyrofluid case.

Each begins with a linear dispersion/growth benchmark and an exact
semi-discrete invariant/source/sink ledger before nonlinear turbulence.

These cases are deliberately cheap enough for student use and aggressive
profiling. They do not displace the direct steady kinetic research spine.

## Link to GK

The scientifically valuable bridge is not simply that "gyrofluid resembles
GK". For common physical observables and inputs, retain a sequence

\[
N_{\rm moments}=N_1<N_2<\cdots
\]

and compare to a gyrokinetic reference. Record numerical error separately from
moment truncation/model discrepancy.

Eventually the adjoint of the higher-fidelity problem can weight the unresolved
moment defect and guide which moments matter to the requested observable.

## Validation

Required evidence includes:

- analytic wave/instability limits;
- exact low-moment invariant identities;
- Scott/GEM or Kendl TIFF/T3FF overlap when exact equations match;
- gyromoment/GK cross-comparison;
- collision/entropy limits;
- adjoint pairing;
- performance scaling with moment count.

A named external code is an oracle only for the overlapping model, not for the
entire gyrofluid category.

Chris&AI
