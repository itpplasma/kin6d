# Interacting jellium cell (Q4)

Status: finite-basis reference, 2026-10-05. Hierarchy-campaign item Q4 and
the first H1 cumulant data. Planned items below are not support claims.

## Model

N electrons (N = 2, 3, 4) in a periodic cube of side L with a uniform
neutralizing background and bare Coulomb interaction. Spin orbitals are
plane waves `V^(-1/2) exp(i 2 pi n.x/L)` with spin; the mode set is the
spherical cutoff `|n|^2 <= c` (c = 1, 2, 3 gives 7, 19, 27 modes, M = 14, 38,
54 spin orbitals) or the cubic cutoff `max|n_i| <= c` (`jellium_modes`; for
c = 1 the cubic set equals the spherical c = 3 set). The Hamiltonian is the
Fock engine's jellium Galerkin projection (`prepare_jellium_hamiltonian`,
shared jellium coefficient `4 pi g/(V q^2)`, q = 0 removed). Atomic units
come from `kin6d_species`; L follows from the shared r_s definition and
`Theta = T/T_F` uses the shared Fermi energy.

Total momentum K and S_z are exact symmetries. Every solve uses a sector
basis (`init_fock_sector_basis`) with the matrix-free coupled action
`fock_matvec_between`, which iterates only nonzero coefficients and fails
closed if one leaves the target basis. No dense many-body matrix is formed.

## Quantities

- Low-lying spectrum of every (S_z, K) sector by FortNum
  `block_lanczos_lowest` (block 8, full reorthogonalization, recomputed
  residuals <= 1e-10 relative). A residual r places an eigenvalue within r
  of the Ritz value; that the kept values are the lowest of their sector is
  a candidate claim.
- Canonical ensembles at Theta from the kept spectrum, with the omitted
  Boltzmann weight estimate `sum_s (d_s - kept_s) exp(-beta(top_s - E0))/Z`
  (valid under the candidate ordering claim). Theta = 0 means the equal
  mixture of the ground multiplet (degenerate open shells for N = 3, 4).
- Shell-averaged pair correlations `g_same`, `g_opp` from the 2-RDM: the
  operator `sum u_ijkl a_i^+ a_j^+ a_l a_k` with the radial-shell average of
  `exp(i Q.r)` (shells inside `|r| < L/2`), normalized by `2 (N/2)^2`.
- Static isothermal response `chi(q) = sum_n w_n chi_n`,
  `chi_n = -(4/V) <A psi_n|(H - E_n)^(-1)|A psi_n>` for
  `A = integral cos(q.x) n(x)`, by GMRES on the K +- q union basis. Only q
  with odd `2 k.q + q^2` (q L/2pi = (100), (111)) are used, so ideal
  denominators never vanish. Oracles: the canonical ideal determinant sum of
  the same cell and orbital set (`jellium_ideal_response`, exact g -> 0
  limit), the grand-canonical cell Lindhard function and RPA
  `chi0/(1 - c S0)`.
- Kicked density mode: `delta<A>(t)/kappa = 2 Im(exp(i E0 t) C(t))`,
  `C(t) = <A psi0|exp(-iHt)|A psi0>` (first order in kappa). Candidate:
  FortNum Cash-Karp RK5(4) on the coupled action; certificate: FortNum
  `certify_linear_trajectory` with the strict interval action
  `fock_matvec_between_interval` (skew-adjoint, mu = 0). The band
  `2 ||A psi0|| e(t)` is certified-discrete for the supplied Ritz pair; the
  Ritz error of (psi0, E0) and the O(kappa^2) kick are separate. A finite
  spectrum shows recurrences, not damping.
- Connected cumulants `lambda2 = Gamma2 - gamma^gamma` and
  `lambda3 = Gamma3 - det3(gamma) - sum_ab (-1)^(a+b) gamma lambda2`
  (`quantum_cumulant_norms`, generated Wick/product leaves) from the
  factorized 1/2/3-RDMs (`fock_transition_rdms`, `fock_rdm3_entries`).
  Fixed-N states have nonzero `lambda3` even when `Gamma3 = 0` (N = 2);
  cumulants are additive over independent subsystems (tested).
- Coupling derivatives for a nondegenerate state: Hellmann-Feynman
  `dE/dg = <psi|V1|psi>`; `psi' = -(H - E)^+ (V1 - dE/dg) psi` (projected
  GMRES); `dg_b/dg = 2 <psi'|O_b|psi>` from the transition 2-RDM; the VJP
  solves the adjoint system with the matrix-free bin operator.

## Verification (`jellium_cell_oracles`)

Full parent-basis legacy assembly plus FortNum dense eigen equals the union
of all sector Lanczos spectra (N = 2, 3, M = 14, 1e-9 relative); the N = 2,
K = 0, S_z = 0 sector equals the Q0 pair-sector constructor (1e-10); the
closed-shell T = 0 ideal response equals `-4/(V e_q)`; at g = 1e-8 the
interacting canonical response equals the ideal determinant sum (1e-5
relative); ideal singlet `g_opp = 1`, `g_same = 0`; `g_same(0) = 0`
(Pauli) and a Coulomb hole for interacting N = 3; RDM traces; Slater
cumulants vanish; cumulant additivity; dE/dg and dg/dg against central
differences (1e-7, 1e-6) and VJP/JVP pairing (1e-8); the certified kick band
contains a dense spectral oracle. Generic engine pieces have their own
tests in `quantum_fock_oracles` (unchanged expectations) and FortNum.

## Results and evidence

`jellium_cell_demo [--config FILE]` (namelist `&jellium_cell`) writes
`output/jellium-{cutoff-scan,pair-correlation,response,cumulants,
mode-history,derivative,cost-ledger}.csv`, `output/jellium_cell_demo.nc`
and `output/jellium-{pair-correlation,response,cumulants,mode-history}.png`.
Key numbers are in the hierarchy-campaign log entry for Q4.

## Open

Spin-flip and point-group symmetries are not exploited; the 62-bit backend
limits M; no momentum-cutoff remainder or thermodynamic-limit bound;
ensemble cumulants include mixing (multiplet) contributions; the
eigenvector derivative needs a nondegenerate state; no equilibrium QMC/PIMC
comparison yet; classical finite-N comparison (QG4 gate) pending C3-C5.
