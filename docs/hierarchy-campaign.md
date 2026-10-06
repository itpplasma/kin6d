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
| Q4 | Interacting N=2..4 jellium cell, matrix-free Fock evolution, weak perturbation | tiny-basis diagonalization | exists (docs/jellium-cell.md) |
| Q5 | Hartree/HF vs Q4 | Q4 | exists (worker, unpromoted; docs/mean-field-cell.md) |
| Q7 | Debye/finite-T Fermi/Thomas-Fermi/Lindhard screening; Wigner -> Vlasov as hbar->0 | analytic limits | exists (worker, unpromoted) |
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
- 2026-10-05 (worker, F1-F4): `gen/kin6d_gen_support.f90` is the one emit
  helper for all generators; the FortSym revision comes only from
  `KIN6D_FORTSYM_REVISION` (configured include). `gen/kin6d_symbolic_physics.f90`
  holds the single symbolic Coulomb/sphere/Ewald/jellium/partial-wave and
  scale definitions; `gen_coulomb` and `gen_species` emit them.
  `src/species.f90` (units with SI map, species, pair/plasma scales),
  `src/ewald.f90` (alpha/cutoffs from a tolerance with a shell-sum truncation
  bound) and `src/nbody.f90` (one FortNum Strang KDK path, any N, tiled
  direct and structure-factor Ewald sums, force/step tangents,
  encounter-resolved step, fail-closed coincidence) are shared by C0, C1, C2;
  Q0 and the Debye bare sphere use the shared definitions. New
  `test_nbody` (Wigner and rock-salt Madelung constants incl. Evjen direct
  sum, Ewald bound, Newton/momentum/energy, FD derivatives, CODATA/hand
  scales). 22/22 CTest green in Release and Debug; regeneration of all 19
  leaf families is deterministic. C1 output changes: Ewald 0.3825/2/4 chosen
  from tolerance 1e-13 (bound 1.7e-15 |g|/L), encounter step re-evaluated
  each step (max energy defect 6.2e-5 vs 2.6e-5, one histogram count moved).
  Open: modal Poisson still duplicated (gen_poisson, gen_shielding_dynamics);
  fp/laboratory/poisson keep their own emit convention; no VJP of the KDK step.
- 2026-10-05 (worker Q3/Q8): Q3 `src/fermion_ensemble.f90` +
  `gen_fermion_ensemble` (docs/fermion-ensemble.md). It provides canonical
  fixed-N box/trap sums through a subtraction-free recursion, interval
  enclosures with a rigorous omitted-orbital remainder, exchange hole,
  density matrix, Boltzmann and ground preparations, and the box EOS with
  T/density derivatives. Oracles: brute-force determinants, cold
  thermodynamic-limit hole (deviation 8.6e-5 at N=1021), grand-canonical
  EOS (e_N-e_TL ~ 1/N), classical limits and finite differences. The
  alternating power-sum recursion cancels up to 5e36 cold; it is used only
  as a cross-check. Q8: relativistic option in `kin6d_fermi_gas`
  (`gen_relativistic_fermi`, convexity tail, graded Fermi-window mesh;
  nonrelativistic path bit-identical). Certified to x_F=10 at 1e7 K;
  NR/UR/Chandrasekhar/classical limits checked by an independent binary128
  oracle. Demos `fermion_ensemble_demo` and `relativistic_fermi_demo` write
  to output/. Open: binary64 range limit for cold large-N states;
  the nonrelativistic certificate fails below theta of about 1e-3 (graded
  mesh not applied); fortplot subplot axis scales leak between panels.
- 2026-10-05 (controller): FortNum main 7f6a0ba adds the generic
  `fortnum_ode_residual_certificate` (F6; 7.8x faster than the kin6d copy
  at equal radius on the 125-mode case) and BSD-sed-portable tests. fortplot
  branch `mp4-quality` (local): BT.709-tagged yuv420p H.264, faststart,
  1 s keyframes, CRF 18, even-size padding, exclusive frame directory,
  edge tick labels, `fill_between` alpha; follow-up for subplot colorbars,
  glyphs, fallback status, log-scale leakage and large-value ticks.
  docs/figure-style.md added. Merged F1-F4 and Q3/Q8; 24/24 CTest green.
  Controller follow-ups: apply the graded Fermi-window mesh to the
  nonrelativistic certificate; per-degree scaling for cold large-N
  canonical sums; FortNum interval layer compiled without
  `-ffp-contract=off` in FortNum's own build (audit directed rounding).
- 2026-10-05 (worker, Q1/Q1b): `src/coulomb_encounter.f90`,
  `src/radial_dvr.f90`, `src/coulomb_kepler.f90`, leaves from
  `gen/gen_coulomb_encounter.f90`, docs/coulomb-encounter.md. Partial waves
  on a GLL FEM-DVR (exact singular integrals in the origin element), banded
  LAPACK eigen windows, exact discrete evolution; product/singlet/triplet
  sectors; classical Wigner ensemble through the nbody interaction step;
  JVP/VJP of bin probabilities in (b, sigma, g). Oracles: hydrogen levels
  9e-13, Coulomb phases 1e-9, F_l shape 1e-8, free Gaussian 5e-7, Mott ratio
  0.018 abs (60-120 deg), Rutherford fluence 1% (120-180 deg), eta=10
  quantum vs classical within 4 sd, derivatives vs FD 1.6e-8, Kepler
  autocorrelation vs exact levels 1e-10 to 1.1 T_rev. 23/23 CTest green in
  Release and Debug. Open: GLL and normalized Legendre belong in FortNum; no
  continuum discretization bound.
- 2026-10-05 (F5 worker): `kin6d_io` (`src/result_io.f90`) plus FortIO
  pinned at 38add03 via FetchContent. All 22 applications read `--config FILE`
  namelists (defaults unchanged, one group per app) and write
  `output/<app>.nc` (NetCDF-4/HDF5) with units, per-variable `kin6d_status`,
  literal/resolved configuration, build-time and runtime git commit/tree
  state, compiler version/options, FortSym/FortNum/fortplot/FortIO pins, UTC
  time and command line. quantum_reference_demo and classical_pair_demo share
  `&periodic_pair` (`src/periodic_pair_config.f90`). The classical demo
  rejects a mass/coupling it cannot realize; the old
  `&quantum_reference length=` spelling is rejected. With defaults, every
  migrated app's numerical CSV/text output matched the pre-migration
  binaries byte for byte. Only timing columns and a few
  comment headers that now print the configured values differ. CTest
  `result_io_contract` checks FortIO-reader, ncdump and h5dump round trips,
  resolve idempotence and rejection cases; the gallery records `.nc` sha256.
  Open: in-file groups, generated schema, external-input hashes, hardware
  records; new apps (coulomb_encounter, ocp_*, jellium_cell, screening,
  wigner) must adopt the module and join the `kin6d_io` list in CMake.
- 2026-10-05 (worker, Q7): `src/screening_hierarchy.f90`,
  `src/wigner_poisson.f90`, `src/wigner_dispersion.f90` (independent Z-ODE
  dispersion oracle), leaves from `gen/gen_screening.f90` and
  `gen/gen_wigner.f90` (deterministic regeneration); FortNum fft compiled as
  `kin6d_fft`; docs/screening-wigner.md. Static ideal-gas RPA (no LFC):
  T=0 momentum route vs closed form 4.7e-15, momentum vs Maldague routes
  8.7e-15, Sommerfeld 2e-8, Dawson limit 3e-5, Friedel period 3.1382/kF,
  rigorous truncation tails, quadrature indicators. Wigner-Poisson split
  step (shared streaming_rotation leaf, exact eta-space kick, forward
  tangent in hbar): Landau/quantum-Langmuir roots within 1e-4, two-stream
  growth within 1.7e-4 for 26 (k, hbar), O(hbar^2) Wigner->Vlasov (rate
  1.99/2.00), tangent = FD to 1e-6, implicit root derivative, unclipped
  negativity 0.25 at t=20 (Vlasov 1e-13 while resolved). Both demos use
  result_io. 28/28 CTest Release and Debug. Finding: fermi_gas candidate k_s
  off by 1.8e-3 at theta~0.025 (inside its certified 1% enclosure). Open:
  Z function and graded Gauss rules belong in FortNum; modal Poisson copy;
  no rigorous Lindhard/Wigner enclosure; screening figure uses log10 data
  workaround for the old fortplot pin.
- 2026-10-06 (controller, pause): merged F5 (namelist + FortIO results),
  Q1/Q1b (partial-wave encounter, Rydberg packet), Q3/Q8, Q7 (screening,
  Wigner-Poisson); kin6d uses FortNum 0f643d0 GLL/Legendre and fortplot
  67834eb; 28/28 CTest green. Unmerged worker branches (resume/integrate):
  worktree-agent-a555f8713dd7a8d74 (Q4 jellium + F6 migration),
  worktree-agent-a6bde7f3f0e416042 (C3/C4 OCP), worktree-agent-a13ae492e23305233
  (Q2 packet models). Remaining: C5/H1 correlation map (g2/g3 cumulants vs
  Gamma, Theta), Q5 Hartree/HF, QG6 hybrid encounters, modal Poisson
  dedup, fermi_gas k_s candidate (1.8e-3 near Theta 0.025) and cold
  certificate mesh, coulomb_encounter_demo result_io adoption, figure pass
  (Q1 movie: CM relative density singlet/triplet/classical, log scale,
  near-closest-approach frames; fine Mott bins), plan/doc consolidation,
  gallery render and slopbox upload.
- 2026-10-06 (controller): Q4 worker finished (branch
  worktree-agent-a555f8713dd7a8d74, 27/27 tests). Before merging, verify
  and push FortNum branch `symmetric-lanczos-eigen` (80b7340, worktree
  /Users/ert/code/lazy-fortran/fortnum-symmetric-eigen) to FortNum main;
  the Q4 branch pins it. Q4 open: Theta=0.5 N=3/4 thermal rows unconverged
  (indicators only); response only at q=(100).
- 2026-10-06 (worker, F6 + Q4): F6: `quantum_evolution` and
  `slab_oscillation` use FortNum `fortnum_ode_residual_certificate`; the
  generated qe residual/Bernstein/radius-interval and slab residual/radius
  leaves are gone. Quantum radii unchanged (two-level 1.0976e-7, tangent
  1.9888e-7, n=125 8.719e-6; ~15 ms/cell, dominated by the kin6d dense
  interval action). The slab moved to the componentwise first-order
  reconstruction: O(h^3) instead of O(h^2) residual, demo radius 6.61e-3 ->
  3.30e-4. Q4: `src/jellium_cell.f90` on new engine pieces (sector bases,
  sparse coupled action with fail-closed leakage, factorized 1/2/3-RDMs,
  `add_fock_one_body`) and new FortNum `fortnum_symmetric_eigen` (dense
  eigen + block Lanczos; branch `symmetric-lanczos-eigen` 80b7340, NOT
  pushed: the pin KIN6D_FORTNUM_REVISION=80b7340 needs the controller to push
  it, local builds use -DFETCHCONTENT_SOURCE_DIR_FORTNUM). `cutoff=2`
  (M=38) default: N=2 E0(rs=2)=-1.668e-2 Ha (cutoff 1/2/3: -1.426/-1.668/
  -1.709e-2); N=4 sector dims up to 73815 states. chi/chi0_Lindhard at
  q=(100): N=2 0.97 (rs 0.5) -> 0.69 (rs 8); RPA overscreens (0.49 at rs 8).
  H1: |lambda3|/|Gamma3| N=3 ground multiplet 0.060 (rs 0.5) -> 0.236
  (rs 8); N=4 0.29 -> 0.41 (multiplet mixture dominated). dE/dg HF vs FD
  7e-10 rel, dg_opp/dg 3e-10 rel, VJP=JVP. Kick band 2e-6 at t=8 (800 RK
  steps, certificate 0.16 s). Theta=0.5 ensembles for N=3,4 are NOT
  converged with nev=8 per sector (truncation estimate 0.05-0.2 / 2-4);
  treat those rows as indicators. 27/27 CTest Release and Debug.
- 2026-10-06 (worker, Q5): `src/mean_field_cell.f90` + `gen_mean_field_cell`
  (docs/mean-field-cell.md): Hartree, RHF, UHF (DIIS + trust-region
  second-order SCF, instability following, multi-start, aufbau-repair
  attempt), restricted Mermin HF and static TDHF/TDH response in the Q4
  basis; Slater states expanded in the engine basis for the full-parent
  residual ||(H - E_HF) Phi|| and a Temple bound (interval-evaluated,
  Lanczos gap input = candidate). New test `mean_field_cell_oracles`
  (functional = <Phi|H|Phi>, brute-force two-mode UHF incl. SDW minimum,
  Koopmans, g->0 Lindhard, TDH = RPA to 1e-12, variational/Temple, HF dE/dg
  vs FD, saddle following, Slater cumulants, Mermin limits/dF/dT = -S).
  `jellium_cumulants` now supports N = 2. cutoff 2: E_c/N(UHF) N=2 -9.3 ->
  -5.9 mHa (r_s 0.5 -> 8; E_HF(N=2) = 0), N=4 -29.8 -> -7.9 mHa; N=2 Temple
  bound/true error 1.25 -> 2.2; N=3, 4 UHF minima break translation symmetry
  (plane-wave determinant is a saddle) and are non-aufbau, so Temple does
  not apply there. chi/chi0 at r_s 8, N=2: exact 0.69, TDHF 0.66, RPA 0.49.
  Mermin vs exact canonical thermal energy is dominated by the
  grand-canonical orbital ensemble (factor 2.7, N=2, r_s 0.5).
- 2026-10-06 (controller): merged Q4 (FortNum 80b7340 eigen/Lanczos pushed),
  cleanup (one modal Poisson / streaming / Krook / P1 definition; fermi_gas
  graded Fermi-surface mesh, candidate k_s now 3e-10 vs Sommerfeld, cold
  certificate to theta 1e-4; encounter demo uses result_io) and Q5 (Hartree,
  RHF/UHF, Mermin, static TDHF vs exact cell; HF symmetry breaking for open
  shells at r_s 0.5). 31/31 CTest, deterministic regeneration. Figure-pass
  items: fortplot log-axis tick labels on sub-decade ranges ("2x10^0");
  Q5 N=4 response panel clipping/legend and symmetry-broken TDHF note;
  Q1 movie redesign; Mott fine bins; Wigner run into saturation.
