# Theory map: verified model reductions

Status: plan, 2026-10-06. Planned items are not support claims.

KIN6D certifies numerical solutions of a named model. This document plans
the matching layer for the models themselves: every model the code solves is
connected to its parent by an explicit derivation, and every derivation step
carries a checked status. The goal is a verified map of plasma theory from
the microscopic parents to the reduced models, with assumptions, orderings,
remainders and errata recorded, so that a model error can be traced to the
step that introduced it.

The numerical certificates ([certification.md](certification.md)) answer
"is this the solution of the model?". The theory map answers "is this the
model we meant, and what did each reduction drop?". The measured model
defects of the microscopic campaign
([hierarchy-campaign.md](hierarchy-campaign.md)) are the empirical side of
the same question.

## Objects

- **Node:** a model with its state variables, equations, invariants and
  structure (Hamiltonian/Poisson bracket, conservation laws, H-theorem,
  Onsager symmetry), its domain and the physical constants it keeps.
- **Edge:** a reduction from a parent node with
  - small parameters and their ordering (e.g. epsilon = rho/L,
    1/Lambda, hbar/(p L), Knudsen number);
  - the derivation steps, each an identity or an expansion at a stated order;
  - what is dropped, the remainder order, and the validity conditions
    (uniformity, time scales, secular terms);
  - which invariants and structure survive exactly, approximately, or not;
  - the source (textbook/paper, page and equation) it is checked against;
  - the numerical companion experiment in KIN6D that measures its defect.

## Skeleton

Classical branch:

1. Maxwell--Lorentz N-body (relativistic) -> Darwin -> electrostatic N-body.
2. N-body -> Klimontovich -> BBGKY hierarchy (exact).
3. BBGKY -> Vlasov (mean field, 1/Lambda -> 0); -> Lenard--Balescu (weak
   coupling, pair cumulant, Bogoliubov time-scale separation); -> Landau
   (static screening cutoff); -> Boltzmann (short-range); Landau ->
   Rosenbluth Fokker--Planck.
4. Equilibrium BBGKY -> YBG equations -> closures (Kirkwood, HNC,
   Debye--Hückel).
5. Magnetized branch: full orbit -> guiding center (Lie-transform
   perturbation theory) -> drift kinetics -> gyrokinetics (gyrokinetic
   ordering, delta-f and full-f) -> gyrofluid/gyromoments.
6. Moment branch: kinetic -> moment hierarchy -> closures (Chapman--Enskog,
   Braginskii, Grad, Landau-fluid) -> two-fluid -> resistive and ideal MHD.
7. Linear response: Vlasov -> dielectric function -> Landau roots,
   Bohm--Gross and cold limits.

Quantum branch:

1. Nonrelativistic electrostatic many-electron Schrödinger/von Neumann
   equation (the QED boundary is external).
2. -> quantum BBGKY for reduced density matrices -> cumulant truncations,
   Hartree/TDHF, quantum Lenard--Balescu, Uehling--Uhlenbeck, and the G1--G2
   nonequilibrium Green-function closures (second Born, GW, T-matrix).
3. -> Wigner representation -> Vlasov as hbar -> 0 (Moyal expansion).
4. Two-body: Schrödinger -> partial waves -> Coulomb scattering kernel ->
   quantum Coulomb logarithm (Bloch correction ln(eta) - Re psi(1 + i eta)).
5. Equilibrium: ideal Fermi gas -> Thomas--Fermi -> Lindhard/RPA ->
   correlated references; Wigner--Kirkwood and Kelbg as semiclassical limits.

Cross edges link the branches (Wigner -> Vlasov, quantum -> classical
collision kernels, Kelbg pair statistics). The implemented KIN6D cases
attach to nodes as numerical companions: H1a to step 4, Q1/QG6 to quantum
step 4, Q7 to quantum step 3, C3 to classical step 7, H1b to the
quantum-classical cross edges.

## Status ladder for derivations

Each step and each edge carries one status, mirroring the numerical states:

1. **stated:** written with assumptions and source;
2. **algebra-checked:** every identity verified symbolically under declared
   assumptions, with agreement of independent engines;
3. **limit-checked:** known special cases and limits reproduced;
4. **defect-measured:** the reduction's error measured by a KIN6D companion
   experiment across its regime;
5. **formally proved:** a Lean theorem covering the step's mathematical
   content, with its hypotheses;
6. **erratum:** a discrepancy with the cited source, recorded with the
   corrected statement and its evidence.

An algebra check does not prove an analytic step. Interchanging limits,
uniform validity of an asymptotic series, existence of solutions, secular
growth and remainder estimates need either a formal proof or an explicit
"unproved assumption" entry.

## Tools: FortSym first, Lean selectively

**FortSym is the primary engine.** It already drives several algebra engines
behind one representation, compares their zero tests, carries an assumption
context, computes series, and generates the numerical leaves of the code.
Deriving the model and generating its implementation from the same symbolic
definitions is the DRY principle extended to theory: a reduced equation that
is checked symbolically is the one that is compiled. Derivation scripts are
Fortran generators like the existing ones and run in CI.

Required FortSym extensions, added upstream when a step needs them:

- noncommutative operator algebra: commutators, creation/annihilation
  operators, truncated Moyal star products;
- multi-parameter asymptotic series with ordering bookkeeping;
- Poisson brackets, symplectic forms and Lie transforms (guiding center,
  gyrokinetics);
- velocity-space Gaussian integrals and moments; gyroaverages (Bessel
  functions); the plasma dispersion function;
- tensor/index calculus; distributions (delta, principal value) for
  Klimontovich and Landau;
- export of a machine-readable step record (assumptions, order, residual,
  engine agreement) and, where possible, a proof trace.

**Lean is selective.** It proves what algebra cannot: the logical structure
and the analytic steps of selected reductions, and the theorems that
numerical certificates rely on. Mathlib's PDE and asymptotic-analysis
coverage is thin, so a full formalization of, say, gyrokinetics is not a
realistic target. Good Lean targets:

- structure of collision operators: conservation of number, momentum and
  energy, and the H-theorem of the Landau form;
- the unitary Duhamel bound and the Hermite/Bernstein residual theorems used
  by the quantum and ODE certificates;
- finite-order Moyal/Wigner identities and the hbar^2 Wigner -> Vlasov
  remainder in a simple setting;
- series identities such as Re psi(1 + i y) = -gamma + y^2 sum 1/(n(n^2+y^2))
  behind the Bloch correction;
- existing certificate theorems (contraction, SPD and adjoint bounds).

Where both apply, FortSym does the computation and emits the identity, and
Lean checks it (ring and field normalization, or a certificate of the
rewriting) rather than re-deriving it.

## First milestones

Each milestone pairs a derivation with an existing KIN6D experiment:

- T1 Vlasov linear response -> dielectric function -> Landau root,
  Bohm--Gross and cold limits (C3, Q7).
- T2 Wigner -> Vlasov at O(hbar^2) via the Moyal expansion (Q7, observed
  rate 1.99).
- T3 partial-wave momentum-transfer sum -> quantum Coulomb logarithm and the
  Bloch identity (Q1, QG6); NRL, Gericke--Murillo--Schlanges and
  Brown--Preston--Singleton as literature checks.
- T4 first equilibrium BBGKY equation -> YBG -> Kirkwood/HNC closures (H1a).
- T5 Lorentz force -> guiding-center drifts and the magnetic-moment invariant
  to first order, then the gyrokinetic ordering.
- T6 Lenard--Balescu -> Landau -> Rosenbluth; conservation and H-theorem,
  with a Lean target for the structural identities.

Errata found against cited textbooks or papers are recorded in a ledger
next to the derivation, with the exact reference and the checked correction.
