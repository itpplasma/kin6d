# Mean-field models of the jellium cell (Q5)

Status: finite-basis model comparison, 2026-10-06. Hierarchy-campaign item
Q5 (QG5/QG7 step "Hartree/HF versus the exact cell"). Planned items below are
not support claims.

## Model and assumptions

The reference is the Q4 cell (docs/jellium-cell.md): N electrons, periodic
cube, neutralizing background, bare Coulomb interaction, plane-wave spin
orbitals in the spherical cutoff. Every mean-field model uses the SAME basis,
N, spin filling, r_s and cutoff, and the same generated definitions: kinetic
energy `qf_kinetic`, interaction coefficient `qf_coulomb` (q = 0 removed) and
the engine's density operator `jellium_density_model`. Each row changes one
assumption relative to the exact sector solution:

| Model | Assumption | Code |
|---|---|---|
| Hartree | product of orbitals, direct term only (self-interaction kept) | `mf_ground(exchange=.false.)` |
| RHF | one Slater determinant, C_up = C_dn (N even) | `mf_ground(restricted=.true.)` |
| UHF | one determinant, separate real orbitals per spin, collinear | `mf_ground` |
| Mermin HF | grand-canonical orbital ensemble at Theta, paramagnetic | `mf_mermin` |
| TDHF (static) | linear response of the UHF state, exchange in kernel and eps (RPAx) | `mf_static_response` |
| TDH = RPA | same with exchange switched off (Hartree reference) | `mf_static_response` |

Orbitals are real combinations of the cell modes; complex (spiral) and
spin-mixing (GHF) determinants are not searched. In the homogeneous cell the
Hartree potential vanishes (q = 0 removed), so Hartree orbitals are plane
waves and Hartree energy equals the ideal kinetic energy; Hartree matters
through the static response (TDH = RPA) and for broken-symmetry states.

## SCF and multiple solutions

T = 0: Roothaan iterations with fixed aufbau filling per spin, Pulay DIIS
(8 vectors) on the commutator `F P - P F` and a constant virtual level
shift; if DIIS oscillates, retries without DIIS and with 4x/16x shifts.
Converged orbitals are canonicalized inside the occupied and virtual spaces.
`mf_ground` scans the spin fillings nup = ceil(N/2)..N (RHF: N/2) and
`nstarts` starts: the plane-wave determinant and symmetry-broken starts
(eigenvectors of T + E_F R, R a fixed pseudo-random symmetric matrix). The
lowest converged energy is kept; the number of distinct converged energies
and of converged runs is reported. The lowest solution is a candidate global
minimum. The real orbital Hessian (A+B) of the TDHF step gives the lowest
stability eigenvalue; zero modes (degenerate open-shell families) are
reported. Mermin: damped density mixing with mu by bisection on 2 sum f = N.

## Error estimator for HF as a model of the full parent

The UHF determinant Phi is expanded exactly in the engine's determinant basis
(`mf_slater_vector`, amplitudes det C_occ[S,:]); `mf_slater_sector` restricts
to the (S_z, K) sectors in which Phi has weight. With the matrix-free Fock
action, rho = <Phi|H|Phi> (= E_HF, tested) and sigma = ||(H - rho) Phi||.
Variationally E0 <= rho. Temple's inequality gives
E0 >= rho - sigma^2/(beta - rho) for rho < beta <= E1 of that space; beta is
the Lanczos E1 minus its residual and the Temple formula is evaluated in
outward interval arithmetic. The window [Temple, rho] is a rigorous
finite-basis statement only under the candidate claim that no eigenvalue of
the space lies between E0 and the Lanczos E1 (ordering not proved); it is
labelled "Temple bound, candidate gap input". Kato/Weinstein
[rho - sigma, rho + sigma] contains some eigenvalue and is only an
indicator for the ground state. Neither bounds the cutoff error.

## Derivatives

dE_HF/dg = <Phi|V1|Phi> (Hellmann-Feynman at the SCF stationary point;
orbital relaxation does not enter): `mf_coupling_derivative`, checked by
central differences of warm-started SCF energies. The exact N = 2 dE/dg
comes from Q4 `jellium_coupling_derivative`, giving dE_c/dg.

## Verification (`mean_field_cell_oracles`)

Mean-field energy functional = <Phi|H|Phi> of the exact engine for converged
and arbitrary symmetry-broken determinants (N = 2, 3; 1e-12); N = 2
opposite-spin UHF in a two-mode cell = brute-force angle minimization of
<Phi|H|Phi> (grid plus golden section, g = 1 and 30; the strong-coupling
minimum is a spin-density wave with zero Hartree energy); Koopmans
(frozen-orbital N-1 energies, 1e-12) and E = sum_occ (eps + t)/2;
g -> 0: HF = ideal ground energy, TDHF = canonical same-cell Lindhard;
exchange off: TDH = same-cell RPA (1e-12); homogeneous Hartree energy = 0;
E_HF >= exact E0 and Temple lower <= sector E0 <= rho; UHF <= RHF;
Hellmann-Feynman = FD (1e-7); all Slater cumulants vanish (N = 2, 3);
Mermin at g = 0 = ideal grand-canonical energy, Theta -> 0 = RHF,
dF/dT = -S. `jellium_cumulants` now treats N = 2 (Gamma3 = 0) instead of
failing.

## Results and evidence

`mean_field_cell_demo [--config FILE]` (namelist `&mean_field_cell`)
writes `output/mean-field-{energy,estimator,response,thermal,cumulants,
derivative,pair-correlation,cost}.csv`, `output/mean_field_cell_demo.nc` and
`output/mean-field-{correlation,pair-correlation,response,estimator}.png`.
Key numbers are in the hierarchy-campaign log entry for Q5.

## Open

Open-shell exact rows (N = 3, 4) are ground-multiplet averages, a single
determinant is not; compare g(r)/chi there with care. No GHF/complex
orbitals, no thermal UHF; Mermin is grand-canonical, the exact thermal rows
are canonical (and unconverged for N = 3, 4 with nev = 8, see Q4). The
Temple gap input is a Lanczos candidate; no cutoff/thermodynamic-limit
bound. TDHF response is static only (no dynamical RPAx spectrum).
