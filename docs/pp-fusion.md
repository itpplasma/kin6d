# Proton–proton fusion from a two-body potential

This is a **candidate, one-body Gamow–Teller calculation** for teaching the
reaction p+p → d+e⁺+νₑ. It solves the nuclear bound and Coulomb-distorted
scattering problems; it does not insert a published overlap or S factor.
The default result is approximately

| Quantity | Calculated | Independent comparison |
|---|---:|---|
| B_d | 2.22460208 MeV | Reid68: 2.2246 MeV |
| P_D | 0.06469611 | Reid68: about 6.47% |
| a_pp, Coulomb modified | −7.748503 fm | low-energy pp scale about −7.8 fm |
| Λ²(0), one body | 6.907589 | historical one-body value 6.92; modern 6.975 |
| S_pp(0), one body | 3.950291 × 10⁻²⁵ MeV b | Solar Fusion II: 4.01 × 10⁻²⁵ MeV b |
| S′(0)/S(0) | 11.2322 MeV⁻¹ | Bahcall–May: 11.2 ± 0.1 MeV⁻¹ |
| f_pp^R(0) | 0.14444456 | Solar Fusion II rounded convention: 0.144 |

The comparison with the full recommended rate is at the few-percent level.
It does not establish its quoted one-percent uncertainty for this older,
one-body model. P_D is potential dependent, not a directly measurable fraction.
No two-body axial current, vacuum polarization, magnetic-moment interaction,
finite-size nuclear Coulomb correction or plasma screening is included.
Numerical, model and physical-input uncertainty bounds are **unavailable**;
mesh/domain differences and residuals are explicitly indicators.

## Sources and fixed conventions

* R. V. Reid, *Local phenomenological nucleon–nucleon potentials*, Ann. Phys.
  **50**, 411 (1968), [doi:10.1016/0003-4916(68)90126-7](https://doi.org/10.1016/0003-4916(68)90126-7).
  The singlet and coupled triplet soft-core coefficients below are retained
  without refitting. The triplet formulas are also printed in Rakityansky and
  Sofianos, *Jost function for coupled partial waves*,
  [arXiv:nucl-th/9710001](https://arxiv.org/abs/nucl-th/9710001).
* R. Schiavilla et al., *Weak capture of protons by protons*, Phys. Rev. C
  **58**, 1263 (1998), [doi:10.1103/PhysRevC.58.1263](https://doi.org/10.1103/PhysRevC.58.1263),
  especially the unit-amplitude scattering convention and Eq. (59).
* E. G. Adelberger et al., *Solar fusion cross sections II*, Rev. Mod. Phys.
  **83**, 195 (2011), [doi:10.1103/RevModPhys.83.195](https://doi.org/10.1103/RevModPhys.83.195),
  Eqs. (19)–(22). These supply the frozen weak-input convention, radiative
  factor, one-body/full comparisons and slope comparison.
* CODATA 2022, [doi:10.1103/RevModPhys.97.025002](https://doi.org/10.1103/RevModPhys.97.025002),
  supplies α=0.0072973525643, ℏc=197.3269804 MeV fm,
  ℏ=6.582119569×10⁻²² MeV s, m_pc²=938.27208943 MeV,
  m_nc²=939.56542194 MeV and m_ec²=0.51099895069 MeV.

Mass symbols below denote rest energies and c=1. Length is fm, energy MeV,
time s; 1 barn=100 fm². The deuteron uses Reid's fitted kinetic coefficient
K_d=41.47 MeV fm². The pp kinetic coefficient uses the physical proton mass,
K_pp=(ℏc)²/m_p. This explicit historical/model convention is not a fit of the
fusion observable. In particular, the frozen older singlet potential with
modern pp kinematics and point Coulomb need not reproduce a modern fitted
pp scattering length to its experimental uncertainty. The test accepts its
literature scale within 0.1 fm and separately checks its actual model value
against an independent ODE solution.

The nuclear-strength controls s_s,s_t default to one. Values near one are
provided to evaluate local model sensitivities, not to calibrate the answer.

## Radial Hamiltonian and tensor reduction

Write Ψ_d=(u_S(r)Y_{01}^{1M}+u_D(r)Y_{21}^{1M})/r with the spin-isospin
factors implicit. Normalize ∫₀∞(u_S²+u_D²)dr=1, so u has units fm⁻¹/².
With x=0.7 r/fm and Y_n=e^(−n x)/x, the singlet strong potential is

    V_s = −10.463 Y_1 −1650.6 Y_4 +6484.2 Y_7                  [MeV].

For the triplet, V=V_c+V_t S_12+V_ls L·S, where

    V_c  = −10.463 Y_1 +105.468 Y_2 −3187.8 Y_4 +9924.3 Y_6,
    V_t  = −10.463 [(1+3/x+3/x²)Y_1 −(12/x+3/x²)Y_4]
           +351.77 Y_4 −1673.5 Y_6,
    V_ls = 708.91 Y_4 −2713.1 Y_6.                            [MeV]

In the ordered (³S₁,³D₁) basis,

    S_12 = [0, sqrt(8); sqrt(8), −2],     L·S = diag(0,−3),

because L·S=[J(J+1)−L(L+1)−S(S+1)]/2. Hence

    H_d = −K_d d²/dr² I + s_t [V_c, sqrt(8)V_t;
                                sqrt(8)V_t, V_c−2V_t−3V_ls]
                            + diag(0,6K_d/r²),
    H_d u = −B_d u,                 P_D = ∫u_D² dr.

The singlet equation is

    [−K_pp d²/dr² + αℏc/r + s_s V_s] χ = E χ.

The origin condition is u(0)=χ(0)=0. There is no softened Coulomb core.
The pointwise expressions and their radial/kinetic JVP/VJP products are
FortSym definitions in `gen/gen_pp_fusion.f90`, with checked-in generated
output. The potential singularities are never evaluated at the removed
origin node. FEM-DVR quadrature of the short-range potentials is an
approximation checked by h/p refinement; no exact-quadrature claim is made.

`kin6d_radial_dvr` supplies the existing normalized GLL basis with coefficients
c_i=√w_i u(r_i). A scalar nodal potential adds to the Hamiltonian diagonal.
The coupled matrix interleaves the two channels at every radius: its lower
half-bandwidth is 2(order−1), not the radial dimension. It stores O(n p)
entries. LAPACK's existing band eigensolver retains a single bound eigenpair;
its reduction workspace is O(n²), explicitly accounted for below. The
Coulomb-only wrapper preserves its original assembly arithmetic and tests.

## Coulomb normalization and scattering length

Define k=√(E/K_pp), η=αℏc/(2K_pp k), ρ=kr and

    C_η² = 2πη/[exp(2πη)−1],
    F₀′G₀−F₀G₀′ = 1,                    (prime in ρ)
    χ(r) → F₀(η,kr) cos δ + G₀(η,kr) sin δ.

FortNum owns the bounded S-wave Coulomb implementation; its Frobenius/log
series and domain are documented in upstream `docs/coulomb.md`. Committed
mpmath data check F,G and both derivatives, including η=0 and attractive
η. The pp application remains inside its supported domain.

For scattering the radial mesh retains R, while the bound mesh removes it.
Set the temporary boundary value y(R)=1 and solve the interior band system
(H−E)c=0. The weak boundary residual determines the derivative rather than a
finite difference:

    y′(R) = √w_R [(H−E)c]_R / K_pp.

This follows by integration by parts: the interior differential residual
vanishes and the boundary term is K_pp φ_i(R)y′(R). Match y=A F+B G,
using the Wronskian, to obtain

    A = G y′/k − G′,       B = F′ − F y′/k,
    χ = y / sqrt(A²+B²),   δ = atan2(B,A).

These formulas set the scattering amplitude, not an arbitrary origin slope.
They are generated leaves with matching derivatives. The independently
computed δ at E=1 keV is about 1.91956×10⁻¹⁰ radians. Its smallness is the
Coulomb barrier suppression, not evidence of a negligible nuclear overlap.

The Coulomb-modified effective-range function is

    K_C(k) = C_η² k cot δ + 2kη h(η),
    h(η) = Re ψ(1+iη) − ln η,
    K_C(k) = −1/a_pp + (r_pp/2)k² + O(k⁴).

The code extrapolates K_C as a quadratic in E to E=0 and returns −1/K_C(0).
It does not report the neutral scattering length or k cot δ as a_pp.

## Overlap and weak rate

In the allowed, one-body approximation only the deuteron S component enters
this initial ¹S₀ matrix element. The D component affects its normalization
and dynamics through the coupled bound-state problem. Schiavilla's convention is

    Λ(E) = sqrt(γ³/2) exp(iδ)/(C_η k) ∫₀∞ u_S(r) χ(r;E) dr,
    γ = sqrt(B_d/K_d).

The phase drops out of Λ². The exported integrand is
sqrt(γ³/2) u_S χ/(C_η k), in fm⁻¹. It integrates to Λ at the representative
energy, not to its square. `u_pp` is χ with unit asymptotic Coulomb amplitude;
`u_pp_over_Ck_fm` is the useful barrier-divided wave, in fm. Neither is
normalized to unit probability in the finite radial box.

At zero energy, Adelberger's Eq. (20) gives in the units used here

    S_pp(0) = [6π² m_p α ln2 ℏ/(100 ℏc)]
              × g_A² Λ²(0) f_pp^R(0) / [γ³ ft].              [MeV b]

The phase-space factor MULTIPLIES this expression. It is not its denominator.
This is equivalent to a rate proportional to G_F² |V_ud|² g_A² Λ² f_pp^R.
The corrected superallowed ft=3071.4 s and g_A=1.2695 are deliberately frozen
to Solar Fusion II. Rather than combining incompatible radiative conventions
for G_F and V_ud, the effective vector coupling is inferred from

    K = 2π³ ln2 ℏ/m_e⁵,    |M_F|²=2,    ft = K/(2 G_V²).

`pp_vector_coupling()` returns G_V=1.14974681×10⁻¹¹ MeV⁻². Thus the weak normalization is independently supplied by beta decay,
not inferred from the target S factor. Solar Fusion II relates its convention
to |V_ud|=0.97418; the inner electroweak correction is already in G_V.

Let W₀=(m_p−m_n+B_d+E)/m_e, W=E_e/m_e, p=√(W²−1), β=p/W.
Neglecting recoil and neutrino mass, energy conservation gives

    f_pp^R(E) = 1.0162 ∫₁^W₀ W p (W₀−W)² F_+(W) dW,
    F_+(W) = [2πα/β]/[exp(2πα/β)−1].

F_+ accounts for repulsion of the positron by the charge-one daughter. The
outer radiative factor 1.0162 is external literature input, not a computed
QED correction. The integral is evaluated with 128-point Gauss–Legendre
quadrature; the independent adaptive quadrature differs by about 1.3×10⁻⁸
in f. Finite-size, recoil and more complete radiative effects retain
unavailable model-error status.

At E=h,2h,4h with h=0.0001 MeV, a quadratic extrapolation gives
Q(0)=(8Q(h)−6Q(2h)+Q(4h))/3. Comparing with 2Q(h)−Q(2h) supplies an
energy-extrapolation indicator. Applying the derivative weights
(−2,5/2,−1/2)/h to S(E) gives S′(0). At these energies the finite-E factor
1/[1−exp(−2πη)] differs from unity by less than 2×10⁻¹¹; neglecting it does
not affect the reported slope precision. This slope concerns the declared
one-body, recoil-free model.

## Derivatives, structure and error accounting

The two public strength products return derivatives of
(B_d,P_D,a_pp,Λ²(0),S_pp(0)). All constants, energy nodes and meshes are frozen.
The coefficient pairing is Euclidean after the DVR mass normalization; the
physical radial pairing is ∫u v dr. Real symmetric potentials and symmetric
kinetic assembly preserve self-adjointness. No dissipative, open-system or
stochastic term is introduced into either stationary Hamiltonian.

For a normalized isolated eigenpair Hc=εc,

    dε = cᵀ(dH)c,
    (H−εI)dc = −(dH−dεI)c,        cᵀdc=0.

`radial_bound_response` fixes a nonzero coefficient for the singular solve,
then projects onto c-perp to enforce the normalization tangent. The sign
convention is u_S(3 fm)>0. For scattering, y(R)=1 is fixed, so

    (H_II−E) dc_I = −(dH c)_I,       dc_R=0.

Generated matching/overlap/phase-space products finish the chain. The VJP
contracts the transpose of the resulting five-by-two parameter response;
no state Jacobian or differentiated LAPACK iteration is stored. Parameter
products are tested against Richardson-extrapolated perturbed solves and a
JVP/VJP dot-product identity. Full derivatives in physical constants and
pointwise certified sensitivity bands are unavailable.

The source ledger for the reaction conserves baryon number (2→2), electric
charge (2→1+1) and lepton number (0→−1+1); the neutrino carries the remaining
energy m_p−m_n+B_d+E−E_e. The two-body Schrödinger calculation is a stationary
nuclear input to that weak transition, not a simulation of time-resolved weak
conversion or a model of the solar medium.

Default numerical settings are order 8, graded elements starting at 0.05 fm,
growing by 1.15 to 0.8 fm, R=100 fm. The measured Λ² change from R=80 to
100 fm is about 7.2×10⁻⁶ (1.1 ppm); increasing order and reducing element
size provides a separate discretization check. The eigen residual is about
4×10⁻¹⁰ MeV, the scattering residual about 3×10⁻¹² MeV. These residuals
and comparisons are not deterministic certificates. The numerical bound,
model-reduction bound and input uncertainty are kept separate and unavailable;
there is no stochastic component. No norm radius is exported as a pointwise band.

The default target is a few-percent teaching observable; tests require ppm
level mesh/domain stability for the declared candidate. The banded solves
use O(np²) factorization, O(np) band storage and bounded RHS scratch. The
inherited bound eigensolver uses O(n²) reduction workspace. OpenMP distributes
nodal assembly and the three independent scattering energies; BLAS threading
is set separately. No performance leadership or rigorous-accuracy timing is
claimed. A single local Release measurement (GNU 16.2.1, Ryzen 9 5950X,
OpenBLAS one thread) of the export including parameter derivatives and I/O
was 1.22 s with one OpenMP thread and 1.00 s with four; peak RSS was about
37 MiB. Every scalar and radial value was identical between thread counts.
These single samples establish only this small-workload comparison.
The focused test is required to stay in the fast tier (<60 s).

## Reproduction and export

No Python or CAS runtime is needed to build, run or test the checked-in
scientific code. The optional independent oracle is
`test/reference/pp_fusion_oracle.py`: SciPy BVP for the coupled bound state,
DOP853 for scattering, mpmath Coulomb matching and adaptive overlap/phase
quadrature. Its committed numeric output is consumed by the Fortran test.
It deliberately shares neither DVR assembly nor generated potential code.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target test_pp_fusion test_coulomb_encounter pp_fusion_export -j 16
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 ctest --test-dir build -j 4 \
  -R '^(pp_fusion_oracles|pp_fusion_export_contract|coulomb_encounter_oracles|kin6d_layering)$' \
  --output-on-failure
# Equivalent focused ordinary Fo entry point:
fo test pp_fusion_oracles
mkdir -p docs/data
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 build/pp_fusion_export \
  --json docs/data/pp_fusion_kin6d.json
mkdir -p ~/proj/fusion-skript/derivations/chapters/data
cp docs/data/pp_fusion_kin6d.json \
  ~/proj/fusion-skript/derivations/chapters/data/pp_fusion_kin6d.json
```

The exporter also writes `output/pp_fusion_export.nc` through FortIO, embedding
literal/resolved namelists and the native provenance contract. The JSON is a
slide interchange copy with node/quadrature arrays, the four channel-potential
entries (DD excludes the centrifugal term), observables, statuses, compiler,
source/dependency revisions and the generating command. It excludes r=0,
where the unsoftened potential is singular. Its radial data have no pointwise
enclosure. The two JSON copies are byte-identical.

Optional configuration is `--config FILE`, group `&pp_fusion`, with
`rmax,hmax,h0,order,energy,strength`. Defaults live only in the Fortran type;
the exporter materializes them in the resolved namelist. Unsupported inputs
fail before any JSON is written. Retained generators are rebuilt by setting
`KIN6D_FORTSYM_BUILD` to a build of the exact CMake-pinned FortSym revision:

```sh
cmake -S . -B build -DKIN6D_FORTSYM_BUILD=/path/to/pinned/fortsym-build
cmake --build build --target gen_pp_fusion -j 8
build/gen_pp_fusion
```

FortNum's bounded Coulomb functionality is tested upstream separately. The
Gremlin lane was attempted with the pinned installed Fo driver, but capture
rejected an external MUMPS CMake package under `/mnt/storage`; ordinary Fo
and native CMake/CTest provide the reported behavioral evidence. No resident
Gremlin success is claimed.
