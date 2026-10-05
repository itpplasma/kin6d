# Microscopic hierarchy campaign

Status: active controller state, started 2026-10-05. This file is the
authoritative task list for consolidating the teaching cases into one
microscopic foundation and building the classical/quantum references that
feed the BBGKY/correlation hierarchy. Planned items are not support claims.

## Goal

One physical foundation, many validated cases. Every microscopic case, classical
or quantum, uses the same species/unit record, the same Coulomb definitions,
the same splitting integrators, the same configuration and result/provenance
contract. On that foundation build references that quantify, as functions of
plasma parameters (Gamma, Theta, r_s, eta):

1. the size and influence of quantum effects (diffraction, exchange,
   degeneracy) relative to the classical finite-N parent;
2. the size and influence of correlations beyond the pair level
   (three-body cumulants, Kirkwood superposition defect) in classical and
   quantum systems;
3. the error of each reduced model (Gaussian packets, Hartree/HF, Wigner,
   Thomas-Fermi/Lindhard, Debye-Hueckel, Landau/Lenard-Balescu closures)
   against those references.

These references are the building blocks for a later BBGKY/cumulant hierarchy
and for adaptive classical/quantum and correlation-order switching
(docs/adaptive-multiscale-uq.md).

## Foundation decisions (F)

- F1 `kin6d_species`/units: one record with charge, mass, hbar, Coulomb
  coupling g = q_i q_j/(4 pi eps0) in declared reduced units and an SI map.
  Physical scales (Bohr radius, Rydberg, Sommerfeld eta, Gamma, Theta, r_s,
  lambda_D, omega_p) are generated from this record, never retyped.
- F2 one Coulomb generator `gen/gen_coulomb.f90`: bare pair potential/force,
  uniform-sphere background, periodic Ewald real/reciprocal/self/background,
  jellium Fourier coefficient, radial partial-wave potential. Primal, JVP/VJP
  and strict interval leaves. All consumers use it; per-case copies are removed.
- F3 one Hamiltonian splitting path (FortNum geometric Strang/KDK) for every
  classical particle case, general N, tiled direct sum and Ewald periodic sum,
  explicit encounter-resolved step control (no softening).
- F4 shared generator support (one emit helper, one FortSym revision source).
- F5 namelist configuration and FortIO NetCDF-4/HDF5 results with embedded
  literal/resolved configuration and provenance for every demo.
- F6 generic Hermite/Bernstein ODE residual certificate moved upstream to
  FortNum; kin6d slab/quantum consumers use it.
- F7 one plot style and robust MP4/GIF pipeline (fortplot owns fixes).

## Case ladder

Classical (C) and quantum (Q) references, each with oracle, derivative and
error status as in AGENTS.md.

| ID | Case | Oracle | State |
|---|---|---|---|
| C0 | Two-body Rutherford encounter (relative coordinate) | Rutherford orbit, energy/angular momentum | exists, migrate to F |
| C1 | Periodic two-electron quench ensemble | Ewald agreement, energy | exists, migrate to F |
| C2 | Isolated electron cloud in neutralizing sphere | Gauss law, harmonic orbit | exists, generalize N |
| C3 | Periodic 3D OCP plasma oscillation (M3) | finite-k dielectric response | planned |
| C4 | Emergent Debye shielding, Lambda scan (M4) | periodic screened Poisson Green function | planned |
| C5 | g(r), S(k), triplet g3, Kirkwood defect vs Gamma (M5) | Debye-Hueckel limit, HNC/literature OCP | planned |
| Q0 | Periodic opposite-spin pair quench | finite-basis certificate | exists, migrate to F |
| Q1 | Two-electron Coulomb encounter, partial waves in relative coordinate, singlet/triplet/classical, eta scan | Coulomb wave functions, Mott formula, free Gaussian | planned |
| Q1b | Attractive electron-proton Rydberg packet (oscillation, revivals) | hydrogen spectrum, Kepler period | planned |
| Q2 | Gaussian/coherent-sum packet models vs Q1 | Q1 reference | planned |
| Q3 | Finite-N ideal fermions: exchange hole, canonical sums | ideal-gas identities | planned (EOS exists) |
| Q4 | Interacting N=2..4 jellium cell, matrix-free Fock evolution, weak perturbation | tiny-basis diagonalization | planned |
| Q5 | Hartree/HF vs Q4 | Q4 | planned |
| Q7 | Debye/finite-T Fermi/Thomas-Fermi/Lindhard screening; Wigner -> Vlasov as hbar->0 | analytic limits | planned |
| Q8 | Relativistic ideal Fermi EOS | classical/NR/UR limits | planned |
| H1 | Correlation map: pair and three-body cumulant size, quantum-classical difference vs (Gamma, Theta) | C5, Q4 references | planned |

## Execution

Controller owns this file, integration and promotion. Workers run in
isolated worktrees, return evidence (tests, numbers, figures) and never
self-promote. Waves:

1. F1-F4 foundation; F6 upstream FortNum; F7 plotting/MP4.
2. F5 I/O; Q1/Q1b; C3-C4; Q3-Q4.
3. Q2, Q5, Q7, Q8, C5, H1.
4. Plan/doc consolidation, gallery render, review.

Build: `FC=gfortran cmake --preset cpu`. Regeneration:
`-DKIN6D_FORTSYM_BUILD=<fortsym build>` at the pinned revision; generators
run from the repository root and must reproduce committed leaves exactly.

## Log

- 2026-10-05: baseline 21/21 CTest green with gfortran 15; regenerating all
  17 leaf families with FortSym a614271 reproduces the committed leaves
  bit-identically. Fixed generator link to use the platform shared-library
  suffix.
