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

T = 0: Roothaan iterations with fixed filling per spin, Pulay DIIS
(8 vectors) on the commutator `F P - P F` and a constant virtual level
shift. DIIS stalls on the flat landscapes of degenerate open shells, so
`mf_newton` finishes with trust-region second-order steps
x = -(A+B)^(-1) F_ov in the eigenbasis of the real orbital Hessian.
`mf_ground` scans the spin fillings nup = ceil(N/2)..N (RHF: N/2) and
`nstarts` starts: the plane-wave determinant and symmetry-broken starts
(eigenvectors of T + E_F R, R a fixed pseudo-random symmetric matrix). Every
run ends with instability following (steps along negative Hessian modes
until (A+B) >= 0) and an aufbau-repair attempt. The lowest converged energy
is kept; the number of distinct converged energies and of converged runs is
reported. The result is a stable local minimum in the real collinear
determinant manifold and a candidate global minimum.

Two findings matter for later packet/orbital models. (1) With 19 modes the
plane-wave UHF determinant of the open shells N = 3, 4 is a saddle (lowest
Hessian eigenvalue about -2.6e-3 Ha); the stable minimum breaks translation
symmetry (Phi spans 37-147 total-momentum sectors) already at r_s = 0.5.
(2) The stable minima are not aufbau (an empty orbital lies below an
occupied one, e.g. 3.448 vs 3.512 Ha at N = 3, r_s = 1) and refilling
relaxes back: the no-unfilled-shell theorem needs a pointwise nonnegative
interaction, and the neutralized periodic Coulomb kernel is not.
Mermin: damped density mixing with mu by bisection on 2 sum f = N.

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
indicator for the ground state. Neither bounds the cutoff error. When Phi
breaks translation symmetry its space is a union of K sectors whose ground
states are nearly degenerate, rho >= beta and Temple does not apply; then
only sigma (indicator) and the variational upper bound remain.

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
Hellmann-Feynman = FD (1e-7); 19-mode N = 3 plane-wave saddle followed to
a lower stable minimum (= <Phi|H|Phi>); all Slater cumulants vanish (N = 2, 3);
Mermin at g = 0 = ideal grand-canonical energy, Theta -> 0 = RHF,
dF/dT = -S. `jellium_cumulants` now treats N = 2 (Gamma3 = 0) instead of
failing.

## Results and evidence

`mean_field_cell_demo [--config FILE]` (namelist `&mean_field_cell`)
writes `output/mean-field-{energy,estimator,response,thermal,cumulants,
derivative,pair-correlation,cost}.csv`, `output/mean_field_cell_demo.nc` and
`output/mean-field-{correlation,pair-correlation,response,estimator}.png`.
Default `cutoff = 2` (19 modes, M = 38), r_s = 0.5, 1, 2, 4, 8, about 3 min:

| quantity | N = 2 | N = 3 | N = 4 |
|---|---|---|---|
| E_c/N = (E_exact - E_UHF)/N, r_s 0.5 -> 8 [mHa] | -9.3 -> -5.9 | -9.8 -> -6.0 | -29.8 -> -7.9 |
| Temple bound / true HF error | 1.25 (r_s 0.5) -> 2.2 (r_s 8) | n/a (broken K) | n/a |
| chi/chi0: exact, TDHF, TDH (= RPA) at r_s 8 | 0.69, 0.66, 0.49 | 0.54, 0.55, 0.32 | 0.34, 0.42, 0.20 |
| exact |lambda2|/|Gamma2|, r_s 0.5 -> 8 (HF: 0) | 0.03 -> 0.31 | 0.06 -> 0.28 | 0.48 -> 0.53 |

N = 2: E_HF = 0 exactly (both electrons at k = 0, no exchange between
opposite spins), so E_c = E_exact; HF g_opp = 1 versus an exact Coulomb hole
g_opp(0.02 L) = 0.86 (r_s 0.5) -> 0.11 (r_s 8); dE_c/dg = -0.036 -> -0.019 Ha.
TDHF is within 0.1% (r_s 0.5) to 5% (r_s 8) of the exact N = 2 response,
RPA is 3-30% low. UHF Hellmann-Feynman versus FD: 1e-10 relative. Mermin
HF versus the exact canonical thermal energy at Theta = 0.5 differs by a
factor 2.7 for N = 2 at r_s 0.5: the grand-canonical orbital ensemble, not
exchange, dominates in so small a cell (U_ideal,GC = 6.39 vs U_exact 2.40
Ha). The exact Theta = 0.5 rows of N = 3, 4 are unconverged (Q4).

## Open

Open-shell exact rows (N = 3, 4) are ground-multiplet averages, a single
determinant is not; compare g(r)/chi there with care. No GHF/complex
orbitals, no thermal UHF; Mermin is grand-canonical, the exact thermal rows
are canonical (and unconverged for N = 3, 4 with nev = 8, see Q4). The
Temple gap input is a Lanczos candidate; no cutoff/thermodynamic-limit
bound. TDHF response is static only (no dynamical RPAx spectrum).
