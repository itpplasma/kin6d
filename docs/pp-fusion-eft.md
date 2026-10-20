# Proton–proton fusion in pionless EFT

Status: **verified implementation of a candidate EFT calculation**. Pions are
integrated out; the degrees of freedom are nonrelativistic nucleons, point
Coulomb photons, and an external weak current. No nuclear potential or fusion
rate is fitted. This supplements the [Reid68 calculation](pp-fusion.md).
The EFT and lattice uncertainties below are estimates, not certificates.

## Primary sources and equation audit

The arXiv PDFs were downloaded and their equations read on 2026-10-20.
Equation numbers below refer to these versions; Butler–Chen was also checked
against the published journal PDF.

| Source | Version; equations used | DOI |
|---|---|---|
| Kong–Ravndal, LO | [nucl-th/9902064v2](https://arxiv.org/pdf/nucl-th/9902064v2), (21)–(23), (28) | [10.1016/S0375-9474(99)00314-0](https://doi.org/10.1016/S0375-9474(99)00314-0) |
| Kong–Ravndal, NLO | [nucl-th/0004038v1](https://arxiv.org/pdf/nucl-th/0004038v1), (30)–(33), (54)–(56), (60) | [10.1103/PhysRevC.64.044002](https://doi.org/10.1103/PhysRevC.64.044002) |
| Butler–Chen (BC) | [nucl-th/0101017v1](https://arxiv.org/pdf/nucl-th/0101017v1), (4)–(10), (13)–(14) | [10.1016/S0370-2693(01)01152-2](https://doi.org/10.1016/S0370-2693(01)01152-2) |
| Ando et al. | [0801.4330v1](https://arxiv.org/pdf/0801.4330v1), (22)–(30), (31) | [10.1016/j.physletb.2008.08.040](https://doi.org/10.1016/j.physletb.2008.08.040) |
| Savage et al., NPLQCD | [1610.04545v2](https://arxiv.org/pdf/1610.04545v2), (7)–(9), (13)–(14), supplement (19) | [10.1103/PhysRevLett.119.062002](https://doi.org/10.1103/PhysRevLett.119.062002) |
| Chen–Liu–Yu (CLY) | [1209.2552v2](https://arxiv.org/pdf/1209.2552v2), (3), (7a–d), (8)–(13), (19)–(21) | [10.1016/j.physletb.2013.02.019](https://doi.org/10.1016/j.physletb.2013.02.019) |
| De-Leon–Gazit | [2207.10176v2](https://arxiv.org/pdf/2207.10176v2), (1)–(2), uncertainty discussion | [10.1016/j.physletb.2023.138093](https://doi.org/10.1016/j.physletb.2023.138093) |

The PLB **450**, 320 (1999) Kong–Ravndal paper concerns scattering lengths;
it is not the LO fusion paper. The later De-Leon–Gazit work uses a tritium
calibration and a different expansion arrangement. Its result is context,
not an additional input to this calculation.

## Inputs and conventions

Masses denote rest energies; length is fm, energy MeV, time s, and
1 barn = 100 fm². Set m = m_p/(ℏc), μ = μ_MeV/(ℏc), and
γ = sqrt(2 m_p m_n B_d/(m_p+m_n))/(ℏc). The reduced-mass convention is
also stated immediately after CLY (3). The shared weak/physical constants
remain in `src/quantum/pp_fusion.f90`; no competing constants module is added.

| Input | Default | Origin/type |
|---|---:|---|
| B_d | 2.224575 ± 0.000009 MeV | measured; Wiringa et al., Table X, [10.1103/PhysRevC.51.38](https://doi.org/10.1103/PhysRevC.51.38), [nucl-th/9408016](https://arxiv.org/pdf/nucl-th/9408016) |
| a_pp | −7.8063 ± 0.0026 fm | Coulomb-modified fit; Bergervoet et al., [10.1103/PhysRevC.38.15](https://doi.org/10.1103/PhysRevC.38.15), reproduced in Wiringa Table VIII and Ando §5 |
| r_pp | 2.794 ± 0.014 fm | same scattering analysis and table |
| ρ_d | 1.764 ± 0.002 fm | deuteron **pole** effective range; BC (7), NPLQCD supplement input list |
| L_1,A(μ) | 3.9(0.2)(1.0)(0.4)(0.9) fm³ | lattice/EFT matching; Savage **v2**, (14) |
| μ_MeV | 140 MeV | PDS convention μ ≈ physical m_π; explicit numerical convention in CLY (11) |
| g_A, ft | 1.2695, 3071.4 s | frozen Solar Fusion II convention, [10.1103/RevModPhys.83.195](https://doi.org/10.1103/RevModPhys.83.195), (20)–(22) |
| G_V | 1.14974681 × 10⁻¹¹ MeV⁻² | derived from that superallowed ft, including the Fermi isospin factor 2 |
| α, ℏc, ℏ, m_p, m_n, m_e | existing CODATA constants | [10.1103/RevModPhys.97.025002](https://doi.org/10.1103/RevModPhys.97.025002); values in pp-fusion.md |
| f_pp^R(0) | approximately 0.1444092 | computed from B_d and the existing beta-decay phase integral, including its frozen outer factor 1.0162 |

ρ_d is the pole-residue parameter, not the zero-energy triplet effective
range 1.760 fm in Wiringa Table VIII. In the S-wave pole approximation
A_S² = 2γ/(1−γρ_d). No potential-dependent D-state probability is required.
The quoted pp parameters retain the literature electromagnetic convention;
higher electromagnetic corrections are not independently recomputed here.

**Hypothesis:** the solely two-body axial quantity changes mildly from the
lattice calculation at m_π ≈ 806 MeV to physical quark masses. NPLQCD uses
this assumption and a 50% additive uncertainty on its short-distance quantity.
Its L_1,A is not a direct physical-mass lattice measurement. The four v2
uncertainties are statistical, fitting/analysis, mass extrapolation, and EFT
matching/higher orders. The v1 numbers 0.1 and 0.3 have been superseded by
0.2 and 0.4. μ ≈ 140 MeV is a subtraction/matching convention, not a pion
propagating in this model. The paper does not specify a more precise μ.

## Threshold derivation and counterterm conversion

With χ = αm/γ, BC (5) fixes the reduced-amplitude normalization. The Coulomb
ladder sum in Kong–Ravndal (22), (28), or NLO (30), (31), gives

\[
 I(\chi)=\int_0^\infty\frac{2x}{e^x-1}
 \frac{\exp[(x/\pi)\arctan(\pi\chi/x)]}{x^2+\pi^2\chi^2}\,dx
 =\frac1\chi-e^\chi E_1(\chi),\qquad J=\chi I.
\]

Define A = exp(χ) − γ a_pp J and z = γρ_d. The LO result is Λ_LO=A.
At zero Coulomb, χ→0 and J→1, so Λ_LO→1−γa_pp. This limit holds with
that scattering length fixed; physically switching electromagnetism changes
which scattering length should be supplied.

BC (9), (10) prescribe the PDS matching, with γ_E the Euler constant:

\[
 D=\frac1a-\mu+\alpha m\left[\ln\frac{\mu\sqrt\pi}{\alpha m}
               +1-\frac32\gamma_E\right],\quad
 C_0^{pp}=\frac{4\pi}{mD},\quad
 C_2^{pp}=\frac{m r_{pp}}{8\pi}(C_0^{pp})^2,\quad
 C_2^d=\frac{2\pi\rho_d}{m(\mu-\gamma)^2},
\]
\[
 \overline L_{1,A}=-\frac{\mu-\gamma}{m C_0^{pp}}
 \left[\frac{L_{1,A}(\mu)}{g_A}-\frac m2(C_2^{pp}+C_2^d)\right],
 \qquad B=-\gamma^2a\overline L_{1,A}.
\]

C_0 has units fm², C_2 fm⁴, L_1,A fm³, and Lbar fm; A,B,z,Λ are
dimensionless. The g_A division belongs to the conversion, not to an extra
factor applied to the lattice input. Keep Lbar fixed when running L_1,A
between scales. Setting L_1,A(μ)=0 leaves range terms in Lbar and defines
a scale-dependent diagnostic, not a unique physical one-body observable.

Expanding BC (7) through ε², with ε subsequently set to one, gives

\[
 c_0=A,\qquad c_1=\tfrac12Az+B,\qquad
 c_2=\tfrac38Az^2+\tfrac12Bz,\qquad
 \Lambda^{[n]}=\sum_{j=0}^n c_j.
\]

CLY (12) requires expanding the **rate** to the same order:

\[
 R_0=A^2,\quad R_1=A^2+2Ac_1,\quad
 R_2=A^2+2Ac_1+c_1^2+2Ac_2.
\]

These R_n are the exported Λ². Squaring the truncated amplitude instead
would give 7.03519 at default NNLO, adding uncomputed higher-order terms.
The separately exported range-resummed diagnostic is
Λ_res=(A+B)/sqrt(1−z), R_res=Λ_res²; it is **not** another EFT order.
BC's K_1,A contribution starts at ε³ in this expansion and involves an
additional pp shape parameter through (9), (10). Without these inputs,
NNLO is the last complete order of this scheme. No N3LO/N4LO claim is made.

The rate uses the existing Solar Fusion II (20) leaf without refitting:

\[
 S(0)=\frac{6\pi^2m_p\alpha\ln2\,\hbar}{100\hbar c}
       \frac{g_A^2 R\,f_{pp}^R(0)}{\gamma^3 ft}\quad[\mathrm{MeV\,b}].
\]

The weak transition conserves baryon number 2→2, charge 2→1+1 and lepton
number 0→−1+1. Positron/neutrino energies share m_p−m_n+B_d+E.
This is a stationary reaction amplitude, not a time-dependent nuclear model.

## Energy slope derived from the same amplitudes

CLY (7a–d), (8)–(11) give the energy-dependent amplitudes. Put t=p² with
p in fm⁻¹ and define F(k)=C_η(k)² exp[(αm/k) atan(k/γ)]. Their real
Coulomb Green integral, after removing its principal-value pole, gives

\[
 J(t)=\frac{2\gamma}{\pi}\operatorname{PV}\int_0^\infty
 \frac{k^2F(k)}{(k^2+\gamma^2)(k^2-t)}\,dk,\quad J(0)=J,
 \quad J_p=\frac{2\gamma}{\pi}\int_0^\infty
 \frac{F(k)}{k^2(k^2+\gamma^2)}\,dk.
\]

For repulsive Coulomb the origin is exponentially suppressed, permitting the
threshold derivative. Penetration-dependent imaginary terms are flat at
threshold and do not change these Taylor coefficients. Expanding CLY (9),
(10) uses αm h(η)=t/(3αm)+O(t²), not h=0. Consequently

\[
 A_p=-\frac{e^\chi(1+\chi/3)}{\gamma^2}-\gamma aJ_p
           +\frac{\gamma a^2J}{3\alpha m},\quad
 B_p=\frac{\gamma^2a^2\overline L}{3\alpha m}
           -\frac12\gamma a^2r_{pp}J,\quad
 C_p=-\frac12\gamma^2a^2r_{pp}\overline L,
\]
\[
 (c_{0,p},c_{1,p},c_{2,p})=
 (A_p,\;zA_p/2+B_p,\;3z^2A_p/8+zB_p/2+C_p),\qquad
 \frac{S'}S=\frac{m_p}{(\hbar c)^2}\frac{R_p}{R}+\frac{f'}f.
\]

R_p is the exact derivative of the order-truncated R above. C_p is required
at NNLO; dropping it changes the slope. The resummed diagnostic uses
(A_p+B_p+C_p)/sqrt(1−z). Its derivative is also calculated, but the same
missing-order qualifications apply. The neutral energy-slope limit is
nonuniform and is not exposed by this repulsive-Coulomb evaluator.

## Results and distinct uncertainty components

S and its estimates below use units 10⁻²⁵ MeV b. Errors are displayed on S;
the JSON separately records amplitude and Λ² estimates.

| Order | Λ², L_1,A=0 | S, L_1,A=0 | Λ², lattice L_1,A | S, lattice | S truncation estimate, lattice | S′/S, lattice (MeV⁻¹) |
|---|---:|---:|---:|---:|---:|---:|
| LO | 6.292551 | 3.597888 | 6.292551 | 3.597888 | 2.732435 | 11.69508 |
| NLO | 6.500133 | 3.716577 | 6.673984 | 3.815980 | 0.830925 | 11.27150 |
| NNLO | 6.806831 | 3.891938 | 7.020265 | 4.013973 | 0.269040 | 11.27989 |
| Range-resummed diagnostic | 7.018690 | 4.013072 | 7.259466 | 4.150740 | 0.273226 | 11.26190 |
| Reid68, one-body | 6.907589 | 3.950291 | — | — | unavailable | 11.2322 |
| Adelberger 2011 | — | 4.01 ± 0.04 | — | — | literature recommendation | 11.2 ± 0.1 |

**Hypothesis/estimate:** with breakdown scale Λ_b=140 MeV, Q=γℏc/Λ_b≈0.32645,
use δΛ_n=Q^(n+1) max_{j≤n}(|c_j|/Q^j). The symmetric displayed rate
estimate is 2|Λ^[n]|δΛ_n+δΛ_n², mapped through the weak prefactor. The
resummed diagnostic keeps the NNLO power-counting estimate. The breakdown
scale stays fixed if the PDS subtraction scale changes. This conservative
rule is not a probabilistic model or a proof of the omitted remainder.

At NNLO the propagated lattice S components, in the same units, are
0.00629 (statistical), 0.03147 (analysis), 0.01259 (mass), and 0.02832
(EFT matching). Other default input contributions are below 0.00351 each.
They are absolute linear sensitivity estimates, not jointly sampled errors.
NPLQCD's systematic already includes scattering-input variations; its matching
error can overlap the EFT truncation estimate. No quadrature sum or combined
confidence interval is asserted. Custom configurations retain reference input
errors from the default sources; those are not refitted uncertainties.
Slope uncertainty, full high-level slope input sensitivities, rigorous numerical
bounds, and higher electromagnetic/model-error bounds are **unavailable**.

## Source differences that must not be hidden

Ando's LO resums ranges and gives Λ≈2.641 (its first input set), whereas
Kong–Ravndal's zero-range LO gives approximately 2.51. This is a counting
choice, not an inconsistency. Ando's l_1A in (29) is not BC's dimensional
L_1,A. Its tritium-based calibration is not substituted for the lattice input.

KR NLO (55), evaluated with γℏc=45.7 MeV, a=−7.82 fm, ρ_d=1.76 fm,
r_pp=2.79 fm, m_N=938.9 MeV and α=1/137, gives Λ²≈6.47852;
the paper's rounded prose quotes 6.45. BC (7),(9),(10), with its stated
inputs and μ=140 MeV, gives c_1/c_0≈0.037899 and c_2/c_0≈0.028593;
(13) prints 0.039 and 0.029. These rounded statements are not used as
high-precision calibration data. The tests retain tight independent equation
oracles and comparisons at the precision justified by each printed number.

More substantially, Savage v2's text near (14) prints the resummed relation
Λ≈2.62+0.0105 L_1,A and quotes Λ=2.659(2)(9)(5). Direct BC PDS matching
and range resummation with the defaults here give 2.649281+0.0115535 L_1,A,
hence 2.694340 for L_1,A=3.9. The NNLO amplitude before resummation is
approximately 2.652393. These are different approximations; a reconciliation
of the resummed numerical conversion is **unresolved**, not attributed to a
specific typo. Savage supplement (19) also prints the triplet C_2 twice;
BC (9),(10) and CLY (11) explicitly distinguish pp and triplet terms and are
the implemented conversion authority. The lattice counterterm is used as
quoted; neither it nor a range is adjusted to obtain 2.659 or S=4.01.

## Implementation, independent checks and reproduction

`gen/gen_pp_fusion_eft.f90` retains all pointwise formulas and generated
primal/JVP/VJP leaves. The public rate Jacobian covers (B_d,a,r,ρ,L_1,A),
with Euclidean parameter pairing; J′(χ)=(1+1/χ)J−1/χ comes from the integral,
not differentiated quadrature iterations. The high-level Jacobian covers Λ²
and S; generated slope leaves also retain products, but their full input chain
is explicitly unavailable. No new upstream special function is duplicated.

With x=πχ tanθ, the production Coulomb integral is
J=(2/π)∫₀^(π/2) C₀²(χ tanθ/2) exp[χ tanθ(π/2−θ)]dθ.
Multiplication by tan²θ gives γ²J_p. FortNum owns C₀² and Gauss–Legendre
quadrature. Orders 128,256,512,1024 are tried until successive integrals
agree; this is an indicator, not an enclosure. The supported domain is
χ=0 for the neutral overlap alone or 0.01≤χ≤1; nonpositive pole residue,
nonfinite inputs and failed quadrature are rejected before export.
The schedule uses bounded O(N_q) scratch, O(N_q²) rule setup and O(N_q)
OpenMP quadrature work. No global matrix or per-node allocation is introduced.

CTest `pp_fusion_eft_oracles` checks independent E1 values, neutral and unit
limits, KR/BC/Ando published comparisons, an optional 70-digit mpmath oracle
using unreduced C_0,C_2 and adaptive phase space, finite-energy principal-value
integrals for the slope, all five rate-input derivatives, transpose pairing,
and agreement with a freshly solved Reid68 result inside the EFT estimate.
`pp_fusion_eft_export_contract` checks JSON/native output and invalid input.
The previous pp tests remain unchanged. No full ctest suite is run.

```sh
cmake --build build --target kin6d_check_quantum -j 8
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg \
  ctest --test-dir build -j 32 -L quantum -LE slow --output-on-failure
# Focused development gates:
fo test pp_fusion_eft_oracles pp_fusion_eft_export_contract
# Optional independent oracle regeneration:
python3 test/reference/pp_fusion_eft_oracle.py
# Symbolic regeneration, using the exact pinned FortSym build:
cmake --build build --target gen_pp_fusion_eft -j 8
build/gen_pp_fusion_eft
# Export only the requested interchange file to fusion-skript:
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 build/pp_fusion_eft_export \
  --json docs/data/pp_fusion_eft.json
cp docs/data/pp_fusion_eft.json \
  ~/proj/fusion-skript/derivations/chapters/data/pp_fusion_eft.json
```

`--config FILE` accepts `&pp_fusion_eft binding,a_pp,r_pp,rho_d,l1a,mu_mev /`.
Defaults live in the Fortran configuration type. JSON records all orders,
separate error components, input sources, comparison, build/dependency/runtime
provenance and literal/resolved namelists. `output/pp_fusion_eft_export.nc`
retains the native FortIO result. Generated PDFs/images are not committed.
Gremlin capture fails at Fo `src/run/fo_cmake_generation.f90:834` on the
external MUMPS CMake input; native CMake/CTest and ordinary Fo are used.
The installed driver SHA256 is
`e011c0d03e8a3e3a73a53fc01e191658335fcb068faf4b471ab9f0645989d5fc`.

Validation on GNU Fortran 16.2.1 Release: 22/22 fast quantum-class tests;
4/4 pp tests plus the layering gate; 2/2 EFT tests in runtime-checked Debug;
2/2 focused ordinary Fo tests. FortSym regeneration is byte-identical.
The changed files have no short-circuit lint findings; repository-wide lint
retains pre-existing warnings and findings outside this change.

A single local export measurement including the fresh Reid68 comparison and
I/O took 0.74 s with one OpenMP thread and 0.74 s with four (OpenBLAS one
thread), with peak RSS 36716/37096 KiB.
The maximum relative change in exported overlaps, rates and slopes was
1.2e-16. These observations establish no speedup or performance ranking.
