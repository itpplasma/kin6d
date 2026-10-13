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

## Scope

The map covers all plasma physics needed for the research and teaching it
serves, from first principles without gaps: every model used by KIN6D, by
its consumers (e.g. NEO-2, SIMPLE, MEPHIT, TIAGO) and by the lectures is a
node connected to the microscopic parents by checked edges. A typical
end-to-end claim is a chain such as

  many-electron Schrödinger (nonrelativistic electrostatic limit of QED)
  -> two-body Coulomb scattering -> quantum Coulomb logarithm
  + classical N-body -> BBGKY -> Lenard--Balescu -> Landau
  -> Rosenbluth Fokker--Planck -> linearized operator with
  momentum-conserving terms -> drift-kinetic (guiding-center) ordering
  -> NEO-2-type linearized collision operator in its Legendre/Laguerre
  representation,

where every edge is derived and checked here and every numerical stage has
its certificate, so that the operator in a production neoclassical code is
traceable both theoretically and numerically to quantum mechanics. A gap in
such a chain is an explicit open item, never an implicit assumption.

## Reproduce, do not rely

Literature results are intake material, not trust anchors. Every proof or
derivation used by the map is reproduced here:

1. preserve the source and its provenance (reference, edition, equation and
   page numbers, checksum where available);
2. state a scoped claim with explicit hypotheses, solution concept and
   non-claims;
3. reproduce it: machine-checked algebra (FortSym), executable oracles and
   KIN6D defect experiments, and a Lean replay where feasible; otherwise a
   complete written proof in our notation;
4. obtain an independent review (`proof-audit` protocol);
5. promote only when these gates pass; disagreements with the source become
   errata records.

A cited theorem is therefore a status on the way to reproduction, never a
terminal status. The operating structure follows the established
frontier-mathematics projects (e.g. the Stafford project): claim records with
status axes (leaf, claim, parent, review verdict, evidence-gate verdict),
proof and checker paths and non-claims; a generated claim DAG; an append-only
frontier ledger and decision log of rejected routes; raw literature intake
kept separate from promoted claims. The `math-frontier`, `proof-audit` and
`proof-map` skills supply the generic protocols. Proposed layout:

- `theory/claims/` one file per claim (node, edge or step);
- `theory/intake/` preserved literature sources and scoped intake notes;
- `theory/reviews/` independent audits;
- `theory/history/` claim DAG, frontier ledger, decision log, errata;
- `theory/derivations/` FortSym derivation programs run in CI;
- `proof/` the existing Lean package, extended per claim.

## Shared map with the courses

One map serves research and teaching. Its node and edge registry (ids, names,
key equations, orderings, course coverage and public evidence) lives in the
public plasma Skript (`map/plasma-models.yaml`) and is rendered there as an
interactive page (https://krystophny.github.io/plasma-skript/map/index.html,
`?course=plasma|fusion|kinetic`, deep links `#node-<id>` and
`#edge-<from>--<to>`) and a static overview figure. The Plasma Physics, Fusion
Physics and Kinetic Theory courses tag the nodes and edges they cover; their
derivations are checked with SymPy, which is the course-level tool.

KIN6D does not copy the registry. Its theory claims, FortSym derivations,
Lean proofs and defect experiments attach to the same node and edge ids.
Because KIN6D is private, the public map shows research status as text only
and never links into this repository. Where a course SymPy derivation and a
KIN6D FortSym derivation cover the same edge, each is an independent oracle
for the other.

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
- **Solution concept per node:** function spaces and the notion of solution
  (strong/classical, weak, renormalized in the DiPerna--Lions sense,
  H-solutions for Landau with Coulomb interaction, measure-valued for N-body
  empirical measures). A step proved for smooth solutions is recorded as such;
  its use where only weak or renormalized solutions are known is a separate
  claim. Weak formulations are also the link to the FEEC discretizations.
- **Formulation level per edge:** whether the reduction is performed on the
  action or bracket (Low Lagrangian, gyrokinetic field theory,
  metriplectic/GENERIC structure for collisions) or on the equations. Action-
  and bracket-level reductions carry exact energy-momentum conservation and
  the H-theorem into the reduced model; equation-level truncations must state
  which they lose.

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
5. **literature intake:** a cited theorem with its hypotheses and
   solution concept, preserved for reproduction (not yet a trust anchor) (e.g. mean-field limits of Dobrushin, Braun--Hepp and
   Lazarovici--Pickl; Hartree limits of Erdős--Yau, Rodnianski--Schlein and
   Pickl; the Lions--Paul semiclassical Wigner limit; Lanford's short-time
   Boltzmann theorem; partial gyrokinetic limits);
   **reproduced:** the intake proof replayed here (checked algebra and a
   complete written proof, or Lean) and independently reviewed;
6. **formally proved:** a Lean theorem covering the step's mathematical
   content, with its hypotheses;
7. **open problem:** no rigorous justification is known (for example the
   derivation of the Landau equation from many-body Coulomb dynamics, or
   uniform-in-time validity of many asymptotic models); the step stays
   formal and says so;
8. **erratum:** a discrepancy with the cited source, recorded with the
   corrected statement and its evidence.

Asymptotic steps record their method (multiple scales, matched asymptotics
for boundary layers such as sheaths, Lie-transform perturbation theory), the
solvability conditions at each order, checks for secular terms, and the time
scale on which the expansion is claimed. The algebra of each order is
machine-checked; uniform validity is an analytic statement with its own
status.

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
- T7 drift kinetics -> radially local DK in (theta, alpha, lambda, x) with
  the finite-shear term -> bounce-averaged precession (NEO-RT Om_tB, shear
  on/off) and the quasilinear switch set; trapped branch limit-checked
  ([DK0](../theory/derivations/drift-kinetics-dk0.md),
  [drift-kinetics.md](drift-kinetics.md)).

Errata found against cited textbooks or papers are recorded in a ledger
next to the derivation, with the exact reference and the checked correction.
