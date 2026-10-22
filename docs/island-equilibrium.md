# General three-dimensional equilibrium with islands

[H] Design accepted with controller decisions, 2026-10-22. This document proposes the next
non-nested increments of [KIN6D](https://github.com/itpplasma/kin6d), rather
than claiming an implemented machine equilibrium solver. It admits islands,
separatrices and stochastic regions without requiring a global flux label or
piecewise-constant pressure. Static isotropic MHD is the first target. Flows,
anisotropic pressure, bootstrap current and kinetic closures require separately
specified equations and are deferred.

[V] [The common core](equilibrium-core.md) owns the physical state D1,
parametrisation D2, Maxwell-stress residual D3, certificate D4, operator chain
D5, boundary D6, evaluator D7 and harness D8. [NN1](equilibrium-nonnested.md)
already defines the closed magnetic two-form and prescribed-trace vacuum solve.
[NN2](equilibrium-nonnested-nn2.md) defines coupled pressure transport,
current continuity and loop voltage. Their equations, generators and tests
remain authoritative. This document specifies closure selection, compatibility,
scaling and the additional gates, without creating a second equilibrium core.

[V] Evidence notation throughout: **[V]** denotes a checked source, repository
contract or mathematical deduction. **[H]** denotes a proposed implementation,
acceptance threshold or resource estimate that remains to be demonstrated.
A source comparison does not establish that KIN6D meets the proposed gate.

## Decisions (controller review 2026-10-22)

[H] Default A (anisotropic pressure transport with declared sources and heat
boundary data) + D (Ohm's law with V or I_phi) is the “relaxed equilibrium
as it is”. An ideal statement requires a measured limit sequence.

[H] Implement B (“tell it the profile”) first by fitting S and conductivity/
current drive for λ on the reference label to reproduce p_target/λ_target
on the seed, freezing, solving A+D and reporting target/fit defects.
This retains a smooth differentiable branch. Π_B is diagnostic first (G1b).
Topology detection inside Newton is non-smooth, so a projected solve is later and optional (I5c).
C is regional and opt-in, and SPEC single volume (G3) verifies only C.

[H] Chris's decision of 2026-10-22 supersedes the staged I7a/I7b plan and
axisymmetric Γ_c fit check: I7 uses full non-axisymmetric 3D BEM from the start,
coupling the KIN6D interior through FortFEM's exterior Dirichlet-to-Neumann
(DtN) API. Reuse the adapter in [draft PR #47](https://github.com/itpplasma/kin6d/pull/47)
(lane `bem3d5`), subject to its consumer gates. A prescribed coil trace on a
large Ω may remain an independent verification oracle.

[H] Prioritise I3/I4 after I1: matrix-free generated JVPs, mode-block
preconditioning and FFT tiles gate all machine work. P2/P3 remain the verified
reference path, with P3 the minimum for the existing high-anisotropy gates.
Start machines with W7-X vacuum 5/5 (EIM), then 5/4, 5/6 and
edge-island candidates, then finite beta. Use Davies et al.'s public normalised
currents with synthetic/private filaments outside the public repo. LHD n=1 RMP
follows only with an LHD-lane calibrated RMP coil set and I11 gates.

[H] Chris's design decision of 2026-10-22 is to retain general FEM,
non-nested fields and explicit error control as KIN6D's purpose. Slower solves
than VMEC are acceptable, but the nested path must pursue VMEC-comparable
DOF counts at matched accuracy. His 19:30 direction makes **high order the
default path** within general isoparametric FEM: p=4–8 and higher
isoparametric/spectral elements and/or periodic B-spline (IGA) angular bases,
with GVEC-like radial splines times high-continuity angular bases. Fourier
is reserved for exact symmetry directions on this new path. Existing NN2
Fourier coupling remains the reference implementation. Nested Newton uses
the exact generated energy Hessian, including a matrix-free Hessian action.

[H] Pursue IH1 → IH2 → IH3 below: high-order nested efficiency first, then
field-aligned/curve-following charts constructed from a VMEC++ or KIN6D nested
solution, then estimator-driven hp-adaptive non-nested solves. These local
charts supersede NN2's unaligned-mesh-only design restriction without imposing
nested topology on the solved field. Use high p in smooth regions and local
h-refinement near islands, separatrices and current sheets, keeping DOFs away
from those regions close to the nested case. High order must demonstrate
control of parallel-to-perpendicular pollution through chi_parallel/chi_perp
of 1e6, rather than assuming alignment alone suffices.

## Literature and the model distinctions

[V] Zotero was consulted first through local metadata/attachment packets on
2026-10-22. The table uses the paper bodies identified in the source ledger
below, rather than treating abstracts or code names as model specifications.
An item key means metadata exists, not that its attachment has been read.
“Not found” means the recorded title query found no item, not library absence.
No Zotero records were modified. Public bibliographic identifiers are retained
here, while downloaded papers and extraction scratch stay outside Git.

[V] The HINT 1989 citation is **Harafuji, Hayashi and Sato**, not Amano.
The original's body was unavailable in this survey, so the HINT-family
algorithm below is checked against Suzuki et al. 2006, §§2–3. Later
implementations must not be inferred to have every restriction of the original.
Likewise, SIESTA uses nested *background coordinates* but permits a non-nested
magnetic field. BMW/EXTENDER reconstruction can reveal islands without solving
a new equilibrium for the reconstructed field.

[H] In the final column, “KIN6D” means the proposed formulation, subject to the
gates below. Avoidance of a representation restriction is not an existence
or uniqueness theorem, and it does not remove static MHD's compatibility laws.

| code | formulation | what it imposes | limitation for us | how KIN6D avoids it |
|---|---|---|---|---|
| HINT / HINT2 [H89], [H06], [H17] | Alternating pressure and resistive/viscous magnetic relaxation on a real-space mesh | HINT2 averages pressure along a finite field-line length and sets it to zero for sufficiently short connections to the exterior | Profile and edge depend on averaging length and relaxation prescription, which are additional model inputs | Declare transport coefficients, sources, open-line boundary conditions and their limits explicitly, with separate pressure and force residuals |
| SIESTA [S11] | Preconditioned energy relaxation with interlaced resistive steps | VMEC-derived coordinate mesh, fixed conducting outer boundary in the 2011 formulation, adiabatic pressure evolution and enhanced numerical resistivity | A nested seed is useful but cannot define the general computational domain or arbitrary final profiles | Physical-space FEEC mesh independent of field topology, explicit closure targets and Γ_c boundary coupling |
| PIES [P86], [P89], [E05] | Alternates magnetic differential equations for pressure-driven current and Ampere inversion | Prescribed profiles with field-line/near-magnetic-coordinate treatment and compatibility in nonintegrable regions | Small denominators and pressure treatment still require a regularisation decision, and convergence is not guaranteed by iteration | Retain its current-continuity physics with finite anisotropic regularisation or invariant-projected profiles, solved as a coupled residual |
| SPEC / MRxMHD [M06], [M12] | Regional constrained energy extrema | Constant pressure and Beltrami field per relaxed volume, helicity/flux constraints, ideal interfaces and total-pressure jump balance | Global stepped-pressure model and selected barriers are not the requested general pressure closure | Continuous-pressure options are the default, while Taylor regions are explicit local choices with interface checks |
| SPEC islands and singular currents [L15], [L16], [Q21], [T23] | MRxMHD island, shielding and saturation studies | Relaxation volumes and surviving interfaces select which islands/current sheets are admissible | Barrier breakup, finite pressure jumps and current flattening affect existence and saturated widths | Diagnose unresolved barriers and current layers, compare only matched closures and report topology and closure sensitivity |
| Free-boundary SPEC [F20] | Vacuum region inside a fixed computational boundary and iterative virtual casing | Exterior coil consistency plus the same MRxMHD plasma constraints | Free boundary removes a geometric restriction but does not change the pressure model | Reuse D6 shell/exterior coupling with the selected pressure/current closure |
| SPECTRE [SP26] | Revised MRxMHD force representation and trust-region least-squares solve | The same constant-pressure subvolumes, Beltrami equations and ideal interfaces as SPEC | Improved numerical robustness does not make pressure continuous within a relaxed volume | Use its island cases as matched regional/vacuum comparators, without importing a stepped-pressure default |
| BMW [B17] | Volume-current and coil Biot–Savart vector potential, then B = curl A | Input currents from VMEC and their numerical reconstruction | Reconstructed islands are not evidence of self-consistent pressure/current force balance | Reconstruct a candidate/seed and re-solve the selected equations, retaining independent D3 checks |
| EXTENDER / EXTENDER_P [E05], [B17] | Generalised virtual-casing field extension from VMEC or PIES | Input equilibrium and casing/interpolation assumptions | Extension and differentiation errors can change islands or divergence, and the interior model is inherited | Use only as a field adapter or exterior check, with commuting magnetic projection and separate trace error |
| JOREK [J17], M3D-C1 [C09], NIMROD [N04] | Time-dependent resistive or extended MHD can approach saturated states | Initial state, resistivity, viscosity, transport, drive, boundary and possibly nonzero flow | Saturation can be time dependent or flow balanced and is not automatically static force balance | Direct stationary solve with an explicit model, while checking matched dynamic end states through D3 and closure residuals |
| MRX [R25] | Admissible-variation magnetic relaxation with FEEC and energy/helicity structure | Initial magnetic data, boundary conditions and choice of variation metric; pressure arises from incompressible projection | Initial invariants select the state, and arbitrary final pressure/current targets are not native independent inputs | Offer explicit pressure/current data and record their compatibility rather than identifying incompressible multiplier pressure with transport pressure |
| Compatible magneto-friction [A26] | FEEC magnetic relaxation with discrete energy dissipation and helicity conservation | Force-free relaxation, specified boundary conditions and a discrete topological barrier | Preserving helicity alone neither preserves every field-line invariant nor supplies finite-pressure closure | Reuse structural identities as tests and keep pressure/current closure separate from the relaxation metric |
| Fractal-pressure ideal theory [K17] | Removes resonant pressure gradients in intervals around dense rationals | A pressure regularity/Diophantine prescription | Shows why an arbitrary smooth profile cannot simply survive all resonances | State the finite regularisation scale and limit path, rejecting unsupported smooth ideal claims |
| VMEC++ / GVEC / DESC [D23] | Nested-surface equilibrium representations | A global foliation and profile laws on it | Islands and stochastic volumes lie outside that representation | Use only for the nested limit with matched pressure, current/transform, flux and boundary data |

[V] These references describe different physical constraint sets. In particular,
MRX's energy-dissipating admissible variations, He et al.'s helicity-preserving
magneto-friction, and B. F. Kraus–Hudson's fractal-pressure construction are not
one algorithm. The latter is not the unrelated M. Kraus kinetic discretisation
work. No unidentified “Chen scheme” is used as an authority for a solver claim.
Energy dissipation is useful structural evidence but cannot certify pressure
compatibility or convergence to the required stationary solution.

### Source ledger: DOI, Zotero and body checks

[V] DOI links below identify published works, with arXiv DOIs explicitly marked
as preprints. Section/page locators identify the body material actually checked.
For paired references, the formulation claim uses the checked body rather than
an unexamined historical item. Missing original bodies are stated explicitly.

- **H89:** Harafuji, Hayashi & Sato, JCP 81, 169–192 (1989),
  [10.1016/0021-9991(89)90069-7](https://doi.org/10.1016/0021-9991(89)90069-7).
  Zotero BFWKJRQ7, 5QTTVZRA, YPLPQFYK exist. Original body unavailable.
  The attribution is checked against publisher metadata and HINT-family sources.
- **H06:** Suzuki et al., NF 46, L19–L24 (2006),
  [10.1088/0029-5515/46/11/L01](https://doi.org/10.1088/0029-5515/46/11/L01).
  Zotero 5KIJFWCX, CQMQLVLK, YJTINDIU exist. The 5KIJFWCX PDF, §§2–3,
  gives alternating steps, the dl/B average and the connection-length cutoff.
- **H17:** Suzuki, PPCF 59, 054008 (2017),
  [10.1088/1361-6587/aa5adc](https://doi.org/10.1088/1361-6587/aa5adc).
  Not found by title query. Bibliographic continuation for tokamak RMP work,
  not a body-verified basis for additional algorithmic claims here.
- **S11:** Hirshman, Sanchez & Cook, PoP 18, 062504 (2011),
  [10.1063/1.3597155](https://doi.org/10.1063/1.3597155).
  Zotero AY6VN2C3, V4AGY6L2 exist. Saved full HTML, §§II–III, VII and IX,
  gives adiabatic pressure, background coordinates, wall conditions and resistive steps.
- **P86:** Reiman & Greenside, CPC 43, 157–167 (1986),
  [10.1016/0010-4655(86)90059-7](https://doi.org/10.1016/0010-4655(86)90059-7).
  Zotero 3K9778E8 PDF, §§2–3, gives magnetic differential equations,
  field-line treatment and Ampere iteration. JMGYPYJM and CM724Z8Q also exist.
- **P89:** Greenside et al., JCP 81, 102–136 (1989),
  [10.1016/0021-9991(89)90066-1](https://doi.org/10.1016/0021-9991(89)90066-1).
  Zotero 2Q4LVYQ6 PDF, §2, details current/profile iteration and near-magnetic coordinates.
- **M06:** Hole, Hudson & Dewar, JPP 72, 1167–1171 (2006),
  [10.1017/S0022377806005861](https://doi.org/10.1017/S0022377806005861).
  Not found by title query. Historical regional-relaxation citation, with a
  corrigendum linked by the publisher. The model equations here use M12's body.
- **M12:** Hudson et al., PoP 19, 112502 (2012),
  [10.1063/1.4765691](https://doi.org/10.1063/1.4765691).
  Zotero IPPUZV64 and ZQVUJQJE PDFs exist. IPPUZV64, §§II–III,
  gives ergodic pressure compatibility, Beltrami regions and interface constraints.
- **L15:** Loizu et al., PoP 22, 022501 (2015),
  [10.1063/1.4906888](https://doi.org/10.1063/1.4906888).
  Zotero ZZZUBN4J PDF, §§II–IV: island/shielding solutions and distinct
  delta-function shielding versus Pfirsch–Schlüter 1/x currents. JE3PHSZZ also exists.
- **L16:** Loizu, Hudson & Nührenberg, PoP 23, 112505 (2016),
  [10.1063/1.4967709](https://doi.org/10.1063/1.4967709).
  Zotero 6IQ8MT4Z PDF, §§II–IV: vacuum stellarator/island verification against
  analytic/Biot–Savart fields and partially relaxed equilibria.
- **Q21:** Qu, Hudson, Dewar, Loizu & Hole, PPCF (2021),
  [10.1088/1361-6587/ac2afc](https://doi.org/10.1088/1361-6587/ac2afc).
  Not found by title query. [Author body](https://arxiv.org/html/2108.02056v2),
  §§3–4: HKT islands/chaos, interface breakup and a maximum admissible pressure jump.
- **F20:** Hudson et al., PPCF 62, 084002 (2020),
  [10.1088/1361-6587/ab9a61](https://doi.org/10.1088/1361-6587/ab9a61).
  Not found by title query. [Author manuscript](https://pure.mpg.de/rest/items/item_3245698_3/component/file_3246059/content),
  §II: computational boundary between plasma and coils, vacuum and virtual casing.
- **SP26:** Balkovic et al., SPECTRE (2026 preprint),
  [10.48550/arXiv.2607.27135](https://doi.org/10.48550/arXiv.2607.27135).
  Zotero WB4EHRK9 exists. [Body](https://arxiv.org/html/2607.27135v1), §§1.2, 2,
  3.3–3.5: unchanged MRxMHD equations, revised force, trust region and island examples.
- **E05:** Drevlak, Monticello & Reiman, NF 45, 731–740 (2005),
  [10.1088/0029-5515/45/7/022](https://doi.org/10.1088/0029-5515/45/7/022).
  Not found by title query. Original body unavailable.
  [Maintainer documentation](https://princetonuniversity.github.io/STELLOPT/EXTENDER.html)
  and B17's body check the virtual-casing reconstruction role.
- **B17:** Lore et al., scraper-element modelling, IEEE TPS 46 (2018),
  [10.1109/TPS.2017.2780624](https://doi.org/10.1109/TPS.2017.2780624).
  Not found by title query. [Author manuscript](https://scipub.euro-fusion.org/wp-content/uploads/eurofusion/WPS1CPR17_17542_submitted-4.pdf),
  magnetic-field modelling section and Figs. 4–5: BMW volume-current A and
  EXTENDER alternatives. Zotero TAG8AT3M is a related Geiger internal-island
  paper, not this item, and is not substituted for it.
- **J17:** Orain et al., NF 57, 022013 (2017),
  [10.1088/0029-5515/57/2/022013](https://doi.org/10.1088/0029-5515/57/2/022013).
  Zotero V7HS7DC6 PDF, model section: nonlinear RMP response with transport,
  resistivity and flows. These equations exceed static isotropic force balance.
- **C09:** Breslau, Ferraro & Jardin, PoP 16, 092503 (2009),
  [10.1063/1.3224035](https://doi.org/10.1063/1.3224035).
  Zotero UVZK9RB3 and J786WHEP exist. UVZK9RB3 PDF, §II,
  specifies the time-dependent M3D-C1 equation set and its representation.
- **N04:** Sovinec et al., JCP 195, 355–386 (2004),
  [10.1016/j.jcp.2003.10.004](https://doi.org/10.1016/j.jcp.2003.10.004).
  Zotero EDERZD9H/CV22MNMK exist. CV22MNMK PDF, §§2–3 and pp. 374–375,
  checks resistive MHD, Fourier/high-order elements and perpendicular pollution.
- **R25:** Blickhan, Stratton & Kaptanoglu, MRX (2025 preprint),
  [10.48550/arXiv.2510.26986](https://doi.org/10.48550/arXiv.2510.26986).
  Zotero 65IGZEMA exists. [Body](https://arxiv.org/html/2510.26986v1), §§III–IV
  and V.3: admissible variations, helicity, Leray pressure and profile selection.
- **A26:** He, Farrell, Hu & Andrews, SISC 48, B165–B183 (2026),
  [10.1137/25M1727540](https://doi.org/10.1137/25M1727540).
  Not found by title query. [Author body](https://arxiv.org/html/2501.11654v2),
  §§2–4: discrete dissipation, helicity, harmonic forms and tangential discontinuities.
- **K17:** B. F. Kraus & Hudson, PoP 24, 092519 (2017),
  [10.1063/1.4986493](https://doi.org/10.1063/1.4986493).
  Not found by title query. [Author body](https://arxiv.org/html/1709.10446v1),
  §§II–III: small divisors and the pressure construction around dense rationals.
- **H22:** Huang et al., PoP 29, 032513 (2022),
  [10.1063/5.0067898](https://doi.org/10.1063/5.0067898).
  Zotero ZXKZTTVM/JAUPBI2M PDFs, §§II–III: HKT current-sheet convergence
  using a nested GS solver and SPEC. This is not a JPP saturated-tearing paper.
- **T23:** Loizu & Bonfiglio, JPP 89, 905890507 (2023),
  [10.1017/S0022377823000934](https://doi.org/10.1017/S0022377823000934).
  Not found by title query. [Body](https://arxiv.org/html/2305.16785v2),
  §§2–4: cylinder, outer-region Delta-prime, SPEC and resistive-MHD saturation.
  Large islands expose the flat-current approximation's limitation.
- **D23:** Panici et al., DESC Part 1, JPP 89 (2023),
  [10.1017/S0022377823000272](https://doi.org/10.1017/S0022377823000272).
  Part 1 not found by title query (Part 2 exists as IAZ3ECDE). Publisher body,
  §§1–2, and R25, §I.3, check the nested-limit role.
  No VMEC++/GVEC performance or non-nested capability is inferred from it.

## Continuous problem and admissible state

[V] On a fixed computational domain Ω, with cylindrical chart (R,Z,φ), seek
B, p and declared closure variables. Coordinates cover physical space, not
magnetic surfaces. The first domain excludes R=0. Physical field components
are Cartesian at the D7 interface, with cylindrical components used internally.
The static smooth target is

\[
 \nabla\cdot B=0,\qquad \mu_0 J=\nabla\times B,\qquad
 J\times B=\nabla p,\qquad B\cdot\nabla p=0.
\]

[V] The last equation follows from the preceding force balance. It is not an
optional profile convention. In the weak/current-sheet class, the D3 stress
and interface traction define equilibrium. A pointwise J requires regularity
or a declared reconstruction, since H(div) alone does not give curl B in L².
Do not sample an unqualified strong current across a tangential discontinuity.

[H] Use NN1's lane-internal A and harmonic cohomology coefficients c_j:

\[
 B=\nabla\times A+\sum_j c_j H_j,\qquad
 A\in H(\operatorname{curl}),\quad B\in H(\operatorname{div}).
\]

[V] The exact complex gives div curl = 0, and each H_j is closed. Retain the
NN1 gauge, boundary trace and harmonic F_vac/R case rather than exporting A
as a new core state. Additional holes/cavities require an explicit harmonic
basis with independent flux/circulation normalisations. Harmonic coefficients
must be constrained by the coil/flux contract, without double-counting coil
flux. The state and its provenance still pass through D1/D2.

[V] Writing J = B×∇p/|B|² + λB, for a compatible pressure and nonzero B,
current continuity requires

\[
 B\cdot\nabla\lambda=-\nabla\cdot(B\times\nabla p/|B|^2).
\]

[V] A closed line requires the right-hand side divided by |B| to have zero
line integral. On nested surfaces, Fourier inversion has denominators
mι−n. A generic resonant pressure drive produces Pfirsch–Schlüter singular
currents, distinct from the delta-function shielding current of ideal
island suppression [L15]. Rationals are dense when ι has shear. Removing a
few resonances does not establish a globally smooth ideal solution [K17].

[V] A continuous p is constant on the closure of each dense field-line orbit.
Thus a connected ergodic volume has constant p, while an island may support
p(ψ_island) on its internal surfaces. B·∇p=0 alone does **not** require
constant pressure throughout an integrable island. Full island flattening
needs an additional transport, mixing or relaxation assumption. In a smooth
ergodic volume, current continuity also makes λ constant, so a local
Beltrami field follows from the equations without imposing a global partition.

## Closure menu and compatibility

[H] A case declares one pressure closure, one compatible current closure,
regularisation parameters, boundary data and the target solution class.
Region-wise combinations require interface conditions, not overlapping source
assignments. Masks must represent physical regions and remain fixed during a
smooth Newton step. Topology detection may suggest a later mask update, but
cannot silently change the problem during differentiation or convergence.
A general static equilibrium is not uniquely selected by force balance alone.

### Case declaration (user view)

[H] Proposed namelist interface: each region selects A/B + D or opt-in C with
regional invariants. This example selects B's frozen fit, then solves A+D.

```fortran
&island_case
 domain='large_coil_free_torus', chart='cylindrical_R_Z_phi', nfp=5
 boundary='full_3d_bem_exterior_dtn', coils='synthetic_filaments'
 region(1)='plasma', pressure_closure(1)='B_frozen_fit', current_closure(1)='D'
 reference_label='seed_label', p_target='p_profile', lambda_target='lambda_profile'
 fit_on_seed=.true., freeze_fit=.true.
 source='fit_S', conductivity_drive='fit_lambda', heat_boundary='declared_heat_data'
 ! A instead uses declared sources; C instead uses a marked Taylor region.
 chi_perp_over_parallel=1.e-4, g_perp_over_parallel=1.e-4, I_phi=0.0 ! or V
 limit_path='refine_then_reduce_both_ratios'
 h_levels=1,2,4, p_levels=3,4,5, K_levels=8,16,24
 observables='B,p,D3,width,phase,flattening,target_defect,fit_defect,trace_defect'
/
```

### A. Pressure transport with sources

[V] Reuse NN2's scalar pressure model, p = n0 T at uniform prescribed density:

\[
 -\nabla\cdot(K\nabla p)=S,\qquad
 K=\chi_\parallel bb+\chi_\perp(I-bb),\quad b=B/|B|.
\]

[H] In the active transport domain, admissibility requires positive finite
coefficients, specified S and Dirichlet/Robin heat data, or compatible Neumann
data plus a pressure mean.
For pure Neumann data, total source must equal outward heat flux. A varying
density or equation of state is a new closure, not a reinterpretation of p.
A null in |B| returns unsupported until a null-compatible transport model is
supplied. No small denominator is hidden by replacing |B| with a floor.

[V] At finite anisotropy, this elliptic pressure equation is regularised but
can have B·∇p ≠ 0. In particular, the NN2 Ampere current ansatz then solves
a regularised stationary model rather than exact static ideal MHD. Its
smooth force defect contains −b(b·∇p) and the perpendicular potential-current
force. D3 must measure both, even when all NN2 algebraic rows converge.
Exact force balance can instead be imposed as an additional compatibility
condition only if sources and boundary data permit it. Reject incompatible
data rather than treating an overdetermined system as an equilibrium.

[H] The ideal claim is a measured sequence in chi_perp/chi_parallel and
current regularisation, with S, heat boundary data, total current and the
path held explicit. Refine h, polynomial degree and modes sufficiently at
each ratio before extrapolation. Raising chi_parallel at fixed chi_perp and
lowering chi_perp at fixed chi_parallel give different source balances and
must not be conflated. Quote only the path actually computed.

### B. Prescribed profiles: frozen fitted sources, then projection

[H] The user supplies a reference label s_ref(x), pressure target f(x) =
p_target(s_ref(x)) and current target λ_target(s_ref(x)), with units and
normalisations. s_ref may come from a vacuum/nested seed or a declared
field-line-averaged label. It is reference data, not an assertion that the
final field has global surfaces. Prescribing ι where it does not exist is
rejected. Use λ = J_parallel/|B| or integral current moments instead.

[H] I5a implements B first by fitting S and conductivity/current drive on
s_ref to reproduce p_target and λ_target on the seed. Freeze these inputs,
solve A+D and report seed-fit and solved-target defects separately, including
incompatible fits. Do not re-fit inside Newton. This keeps the smooth branch
and its equation-level derivatives without topology detection in the solve.
I5b adds Π_B as a diagnostic gated by G1b on the operator itself. I5c may later
use it in an optional solve with topology-change and unavailable-derivative status.

[H] Define V_B = {q in L²(Ω): B·∇q = 0 weakly}, with compatible boundary
data and means, and define Π_B f as the orthogonal projection onto V_B in
a declared positive weighted L² inner product. With volume weight one,
this is conditional averaging over invariant components. On a closed line,
it is the invariant-flow average with time measure dl/|B|, not uniform
arclength. On an ergodic volume it is a volume average. On a surviving
irrational torus it is its invariant surface measure average.

[H] For the projection diagnostic or optional solve, record either p = Π_B[p_target(s_ref)] or
p = p_target(Π_B s_ref). They differ for nonlinear targets. The default is
the first, because it preserves the component mean of the supplied pressure.
On surfaces for which s_ref is already constant, the requested profile is
enforced exactly. Otherwise “enforced on surfaces” means the explicitly
specified invariant average, and its mismatch from the pointwise reference
is reported. A finite trace-time average is an approximation to Π_B, not a
proof that a chaotic orbit has explored its invariant component.

[H] To request flat pressure over a whole island chain, merge its internal
invariant-surface sets into a declared mixing component and take one weighted
mean there. Without that explicit choice, pressure may vary between internal
island surfaces. Chaotic components flatten only to the extent supported by
the stated projection/trace tolerance. Open field lines inherit compatible
heat/pressure boundary data rather than an arbitrary closed-orbit average.

[H] A later alternative is a penalised projection: with κ>0, solve
−div(K∇p)+κp=κf, taking chi_parallel large and chi_perp positive.
At fixed chi_perp its limiting constrained minimiser includes perpendicular
smoothing and is generally **not** Π_B f. Recovering the unsmoothed projection
requires a separately resolved chi_perp→0 limit and boundary compatibility.

[V] A pressure projection does not solve the resonant-current problem. First
compute a particular λ_PS from current continuity. Its component mean is
undetermined, so only the invariant homogeneous part λ0 in
λ=λ_PS+λ0 may be matched to Π_B λ_target. If the magnetic differential
equation is incompatible, use declared finite current regularisation (D)
or admit a weak current sheet. Do not overwrite λ_PS with a target profile.
Dense-rational singularities remain a possible failure of a smooth ideal
limit even when the projected pressure is well defined.

### C. Optional Taylor/Beltrami regions

[H] Inside a user-marked region T, request curl B = μ_T B and constant p_T.
A regional helicity and fluxes select μ_T, or μ_T and a field amplitude/flux
may be supplied when the boundary-value problem is solvable. Gauge-invariant
relative helicity is required when normal flux crosses a reference boundary.
A closed relaxed region normally has B·n=0 and fixed helicity/flux data.
A geometric mask cut through field lines is not automatically such a region.

[V] Across a magnetic interface require continuous normal B and the weak
Maxwell traction, including [[p+|B|²/(2μ0)]]=0 when B·n=0. Tangential jumps
carry declared sheet currents. If B·n is nonzero, use the full traction and
compatible pressure/current traces, rather than imposing the magnetic-surface
jump formula. Internal pressure laws outside T remain continuous-pressure
closures A or B. Selecting C globally would recover a Taylor model, not the
default general equilibrium. SPEC single-volume comparisons use exactly this
restricted problem and cannot verify A or B.

### D. Current closure and drive

[V] NN2 defines anisotropic Ohm/current continuity with electric potential Φ,
loop voltage V and finite perpendicular conductivity regularisation. Reuse its
actual current law and generated derivatives. The inductive field is
V∇φ/(2π) for a full-torus loop voltage, with Φ single valued. Prescribe either
V or net toroidal plasma current I_phi and solve for the other. I_phi is the
meridian integral of J_phi, averaged over toroidal angle. Do not prescribe
both independently unless a compatibility test supplies the redundant datum.

[H] Specify potential gauge and electrical boundary conditions, conductivity
profiles and heat balance separately. Finite perpendicular conductivity
regularises closed-line current solvability but introduces a perpendicular
force defect, which must decrease on an ideal-limit sequence. A prescribed
net current supplies one global constraint, not a complete parallel-current
profile. Profile constraints from B require the homogeneous-current freedom
and must be rejected if they conflict with the chosen Ohm law or regional
helicity. Bootstrap is a later explicitly coupled current provider, with its
own residual and uncertainty, never a hidden correction to an ideal solve.

## Free boundary and external drive

[H] I7 uses **D6 Γ_c FEM–BEM coupling with full non-axisymmetric 3D BEM**
from the start through FortFEM's exterior DtN API. The fixed boundary Γ_c of Ω
lies in vacuum, inside every external coil and outside plasma/current support,
without an axisymmetry requirement. Ω contains a plasma region and a vacuum
collar where p=0 and J=0. No globally nested plasma boundary is required.
The plasma extent needs an explicit support/edge law as part of the closure,
with matching to the vacuum region. For B, open-line pressure may be tied to
zero exterior pressure. For A, specify the heat-loss boundary and plasma
support problem, solving a sharp edge with its weak traction conditions when
that model is chosen. A source-free collar is not automatically vacuum at
finite perpendicular transport. Verify vanishing pressure and current there,
or report the exterior-model defect and enlarge the domain. Never zero the
collar after solving. No free-boundary claim is made without this support law.

[V] Coils supply B_coil through Biot–Savart. Existing `coils_3d` polylines and
D6's axisymmetric coils are the first adapters. Filaments must stay outside
Ω or inside excluded cavities with separately specified topology. Finite-build
coils integrate a declared current density over their cross-section/volume,
with polygon, quadrature and source-geometry errors reported independently.
No physical coil lies in a cell declared current-free vacuum.

[H] Decompose exterior response from imposed coil field, enforce the normal
and tangential field matching through Γ_c, and carry noncontractible fluxes
explicitly. Reuse M4's D6 contract and the exterior DtN adapter in
[KIN6D draft PR #47](https://github.com/itpplasma/kin6d/pull/47) (`bem3d5`),
subject to consumer verification. Its unmerged status does not establish a
verified free-boundary capability of this design's baseline.

[V] [FortFEM's general parametric 3D Laplace BEM](https://github.com/lazy-fortran/fortfem/blob/main/doc/PARAMETRIC_LAPLACE_BEM_3D.md)
provides P1 potential/P0 flux traces on geometry curved to any order and a
cached exterior DtN API ([#84](https://github.com/lazy-fortran/fortfem/pull/84)).
Curved P2 potential/P1 flux traces are available through
[#93](https://github.com/lazy-fortran/fortfem/pull/93).
The refined-solve flux error estimate in
[#88](https://github.com/lazy-fortran/fortfem/pull/88) applies to P0 fluxes and
is conditional on saturation and quadrature convergence. It is not a
certificate and does not apply to P1 fluxes. Rotational block-circulant P2 assembly
([#105](https://github.com/lazy-fortran/fortfem/pull/105), queued) reduces
repeated-sector work, but the KIN6D gate still stalls on coincident pairs at
quadrature degree 96.

[H] Close generic coincident-pair quadrature and any further high-order trace
gaps in FortFEM, without duplicating BEM in KIN6D. Add a compressed operator
(H-matrix/FMM/ACA) there only if measured dense memory and solve/assembly cost
at W7-X surface sizes require it. Consumer accuracy gates remain unchanged.
A prescribed coil normal trace on a large coil-free Ω is an optional
verification oracle, exact for vacuum with the harmonic/flux data retained.
At finite beta, estimate its plasma-response trace defect against the coupled
DtN solve and by Ω enlargement, or mark it unavailable until measured.

[H] Virtual casing is retained as an independent exterior oracle and optional
adapter for a sufficiently smooth enclosing source surface. It is not the
primary nonlinear coupling: Γ_c already supplies the shared D6 residual,
coil traces and derivative contract without repeatedly reconstructing a
plasma casing. Virtual casing still needs consistent surface fields,
near-singular quadrature, orientation and circulation treatment [F20], [E05].
Agreement with a casing field alone cannot verify interior force balance.

[H] Field-period symmetry uses φ in [0,2π/nfp) and modes n=k nfp. Full torus
uses nfp=1, including nonsymmetric sine/cosine components. nfp specifies the
symmetry retained by **all** coils, masks, sources and targets, not merely
by the unperturbed device. Any n=1 LHD RMP requires nfp=1 even though the
main helix has ten toroidal periods. Voltage remains a full-torus voltage,
and period-integrated energies/currents must be converted with declared factors.

## Discretisation and scalable nonlinear solve

[V] NN2 stores dense local derivative blocks scaling as
O(E d² M²), where E is the meridian element count, M=2K+1 is the real
Fourier-block count and d is the local coupled row count per block.
`nonnested_coupled.f90`'s `element` constructs coupled magnetic, pressure and
potential blocks, with the gauge/multiplier handled in the global system.
Compression afterwards does not remove the peak local quadratic storage.
This is unsuitable for machine-scale Fourier resolution.

[H] Use matrix-free Newton–Krylov, often called JFNK, with **analytic generated
JVP actions**, not subtractive finite-difference Jacobians. Gather a vector,
apply discrete curl, synthesize Fourier fields, apply quadrature-point JVP
leaves, transform/test and scatter. Generate primal/JVP/VJP leaves from the
same KIN6D FortSym definitions. Gauge and harmonic constraints, current rows,
boundary response and closure derivatives participate in the same action.
Keep the small assembled NN2 solve as an independent discrete comparator.

[H] Use flexible GMRES with a bounded restart and a block-triangular physics
preconditioner: mode-wise curl-curl/gauge solves, anisotropic scalar solves
for p and Φ, and a bordered Schur solve for V/current constraints. Assemble
only these sparse meridian blocks in FortSparse. Start with coefficients
averaged over φ, including the full axisymmetric field. D5's axisymmetric
force operator may augment the magnetic block where its constraints match.
It is not the whole NN2 Jacobian. Strong mode coupling requires measured
adjacent-mode blocks or a coarse correction, without restoring dense M² storage.

[H] Use FFT synthesis/testing on tiled quadrature data, with O(M log M)
toroidal work rather than direct O(M²) transforms. The rational nonlinear
coefficients 1/|B| and b b have unbounded harmonic content, so the usual
3/2 padding is only a starting rule. Increase toroidal quadrature at fixed
K until residuals and observables stabilise independently of the mode ladder.
Partial assembly stores fixed geometry/basis data and bounded tile scratch,
not E×quadrature×all-modes×all-variables global copies.

[H] At fixed degree and Krylov-vector count, coefficient storage is O(D M),
with D the globally assembled spatial unknown count per real block. Operator
work is O(E Q [d M + M log M]) for fixed quadrature count Q. A Krylov solve
costs its measured number of JVP and preconditioner actions, plus reductions.
No mode-independent iteration count or wall time is assumed. Mode-block
preconditioner factor fill can dominate memory and must be measured separately.

### High-order efficiency gates (IH1–IH3)

[H] IH1 closes the proposed nested DOF gate: at the same G6 accuracy and
matched model/input data, KIN6D uses **at most 2× VMEC++'s independent real
solution DOFs**. This factor is a hypothesis, not measured performance.
[The nested benchmark](../benchmarks/nested/README.md) (lane N1) reports DOFs
and their counting conventions, achieved accuracy, setup/solve wall time,
cores and peak RSS per code for VMEC++, GVEC, DESC and KIN6D, including its
high-order variant once available. Report constraint/auxiliary unknowns
separately so different state parametrisations remain interpretable. A failed
2× gate leaves IH1 open for basis/order refinement, with no efficiency claim.

[H] Upstream prerequisites are high-order/IGA compatible bases and generic
p-multigrid or low-order-refined preconditioning in
[FortFEM](https://github.com/lazy-fortran/fortfem), and sum-factorisation kernels
generated through [FortSym](https://github.com/lazy-fortran/fortsym), with
missing generic generation support added there. I3/I4 (matrix-free lane I34)
consume these operators rather than implementing them privately in KIN6D.
Target O(p^4) work per 3D tensor-product element at quadrature order O(p),
and measure operator cost per element versus p=2,3,4,6,8 and higher before
claiming that scaling. High-order states/observables must agree with the
converged P2/P3 reference within combined estimated errors on existing ladders,
with all stricter G0–G8 gates retained.

[H] IH3 freezes each variant's DOFs across chi_parallel/chi_perp=1e2,1e3,
1e4,1e5,1e6 and compares low and high order at matched total DOFs. Use an
independent refined oracle at each ratio, fixed sources/boundary data and
declared pressure/flux scales with absolute floors. Require flattening error
|G−G_oracle| ≤0.01 and perpendicular heat-flux error ≤1% of the declared
scale throughout, with high order reducing both errors relative to low order
at 1e6. Keep the actual physical ratio dependence of G separate from pollution.
Apply I1's closure-specific flattening interpretation, including no imposed
whole-island flattening for B without declared mixing. Failure keeps the
high-anisotropy efficiency claim open.

### Machine resolution and cost envelope

[H] The following are sizing examples, **not verified sufficient meshes or
machine runtimes**. Let D=10^6 real unknowns per Fourier block, E=10^5,
d=50 effective local rows and Q=25 meridian quadrature points. These are
rounded planning assumptions, not exact P3 local-space dimensions. One double
vector costs 8 D M bytes. Reserve 70 vectors for restart-30 flexible GMRES:
31 Krylov vectors, 30 preconditioned vectors and nine state/residual/work
vectors. Count the
actual allocations and preconditioner scratch before requesting a node.

| target | retained Fourier ladder K; physical n_max | M at middle K | one vector | 70-vector floor | one dense E(dM)² double-block estimate |
|---|---|---:|---:|---:|---:|
| W7-X, nfp=5 | K=8,16,24; n_max=40,80,120 | 33 | 0.264 GB | 18.48 GB | 2.18 TB |
| LHD RMP, nfp=1 | K=24,32,48; n_max=24,32,48 | 65 | 0.520 GB | 36.40 GB | 8.45 TB |

[H] Values use decimal GB/TB and omit sparse factors, geometry, BEM and
quadrature buffers. They are memory floors for the matrix-free route and
illustrations of the dense-block problem, not job requests. At Q=25 the
middle grids have 82.5/162.5 million element-quadrature-mode samples per
residual/JVP before padding. With 3/2 padding these are roughly 125/245
million samples. Measure bytes/sample and operations/sample on the verified
cylinder before converting to seconds. Reuse factors only while coefficient
changes preserve the declared preconditioning quality.

[H] Machine order is W7-X vacuum 5/5 (EIM), then 5/4, 5/6 and edge-island
candidates, then finite beta, including internal-island cases. Use public
normalised currents with synthetic or privately held filaments, and run
private geometry outside the public repository. Use the convention n/m
in these machine names, so ι=n/m, while the helical phase below is mθ−nφ.
The meridian mesh must resolve m=4,5,6 chains and local widths independently
of the Fourier ladder. Fivefold periodic runs exclude symmetry-breaking
responses by construction, so a full-torus check is needed if the input or
solution under investigation breaks that symmetry.

[V] Davies et al., [arXiv:2607.14722](https://arxiv.org/html/2607.14722v1),
§4.2/Table 2, gives normalised filament currents, including EIM's unit NPC
and zero planar currents. Body checked after a local Zotero identifier query
found no match. Currents do not supply geometry.

[H] LHD uses the l=2, toroidal-period-10 main helix with an n=1 RMP. Retain
physical toroidal harmonics at least through 2×10 and their sidebands:
K=24 includes 20±1, while 32 and 48 test additional nonlinear couplings.
Do not keep only multiples of ten or equate the poloidal helical number l
with a toroidal mode. These cutoffs are initial ladders and must be extended
if width, phase, pressure or force residual remains sensitive to the tail.
Begin only after the LHD lane supplies a calibrated RMP coil set and I11 gates pass.

[V] NN2's recorded high-anisotropy closed-line tests pass at P3, while P2 can
show mesh-dependent Φ pollution and fail the voltage-order gate even with
a converged nonlinear solve. At fixed unaligned polynomial space,
chi_parallel→∞ can lock the pressure into the wrong discrete near-kernel.
P3 is a starting requirement, not a proof of uniform accuracy.

[H] IH1–IH3 use high order, h/p and anisotropy ladders, with the same B in
all magnetic/transport operators. IH2's field-aligned local charts must preserve
commuting magnetic operators, conservation and consistent primal/JVP/VJP
semantics across separatrices without demanding global nested coordinates.
An FCI parallel operator remains an alternative with those same gates. Test
an asymptotic-preserving transport correction in its owning library against
the existing closed-line oracle, rather than changing a pressure gauge ad hoc.

[H] Element/quadrature tiles, Fourier batches, tracing seeds and estimator
patches use OpenMP on CPU. Hot pointwise kernels are FortSym-generated with
checked-in output and no per-call allocations. Native CUDA is optional behind
a thin C ABI, with CPU comparisons at matched observable accuracy. Add
missing generic code-emission targets in FortSym. Report 1/2/4/8-thread
scaling and GPU transfer/solve cost. A serial production hot path is a defect,
while a small serial oracle is legitimate. Large cases go to Slurm, and
resources follow measured peak memory, never these rounded floors alone.

## Estimators, observables and certificate status

[V] Reuse D3's Maxwell-stress residual with its background-field subtraction
and σ normalisation. NN1's finite-test-space dual norm is not automatically
a continuum majorant. Report the test space, quadrature ladder, normalisation
and boundary terms. D4 needs a majorant and a justified stability/inf-sup bound
on the gauge complement before a residual becomes a state-error bound.
At bifurcations or separatrix changes this bound may be unavailable.

[H] Every result includes separate Ampere/current, pressure-transport or
projection, Ohm/Φ, total-current, interface and Γ_c residuals. Measure
algebraic div B from the exact complex and physical sampled divergence with
geometry/trace error separately. Use ε_parallel = L RMS(b·∇p)/p_scale and
ε_force from D3, with explicit positive pressure/length scales. Neither a
small Newton step nor a small algebraic residual can waive these checks.

[H] D7 topology output records the section, seed set, field-line integrator,
tolerance, trace length and chart conventions. Find elliptic/hyperbolic fixed
points of the appropriate q-turn map and track the same island through each
ladder. Report O/X positions, residues, island radial full width W on a
specified chord and helical phase α=mθ−nφ relative to the coil phase.
Phases are compared modulo 2π and become unavailable as W approaches zero.
Separatrix width and finite-time libration width are different diagnostics.
Chaotic layers may prevent a unique separatrix width, which is reported as
unavailable or as an explicitly defined finite-time interval.

[V] NN2's flattening diagnostic G is the absolute central pressure slope at
the O point divided by the absolute slope across a larger chord. Its current
fixture uses offsets ±0.03 m and ±0.12 m respectively in
`test_nonnested_island.f90`. G≈0 indicates local flattening, while G≈1
indicates an unchanged local slope. A zero denominator makes G undefined.

[H] Generalise those offsets as declared fractions of W, frozen across a
comparison. Also report pressure spread along traced lines and across each
mixing/chaotic component, with independent line-length refinement. Do not
force G→0 for closure B without whole-island mixing. Use separate h, p, K,
quadrature, coil discretisation and trace ladders for W, α, G and pressure.
Error bars are conservative refinement envelopes plus tracing/input budgets,
with assumed asymptotic rates stated. They are practical estimates, not
certified continuous bounds or statistical confidence intervals.

[H] Claim labels remain D4's candidate, residual verified, practical estimate
or certified with named discrete/continuous scope. Current-sheet strength,
closure-limit extrapolation and topology classification each have their own
available/failed/unavailable status. Keep numerical, input and model errors
separate. Equation-level tangents/adjoints apply on a fixed smooth closure
branch, while topology/mask changes and an unbounded inverse return
unavailable derivatives rather than differentiated solver iterations.

## Verification ladder before machines

[H] These are proposed gates for new I-lane cases, not replacements for any
stricter existing NN1/NN2 test. Use dimensionless scales stated in each fixture.
All gates need nonstationary starts, source/input provenance, full residual
reports and an independent oracle. Tolerances are frozen before running.
Unless stated otherwise, require blockwise nonlinear backward error ≤1e−10,
scaled algebraic divergence ≤1e−12, and at least three mesh/mode levels.
For exact smooth Pk cases require B order ≥k−0.3 and falling D3 with order
≥k−0.3. A finite-test D3 threshold is only a fixture gate, not a certificate.

| gate | problem and independent oracle | acceptance before promotion |
|---|---|---|
| G0 | Existing NN1 vacuum/Solov'ev and NN2 manufactured, transport, coupled, island and closed-line cases | Existing thresholds unchanged, including P3 high-ratio gates; assembled and matrix-free states agree to 1e−8 relative |
| G1a | Closure A manufactured anisotropic transport with analytic B,p, finite-ratio sources and heat boundaries | L² p order ≥k+0.6 at ratio one, existing high-ratio oracle gates, source/heat balance ≤1e−9; force-incompatible manufactured data are labelled as such |
| G1b | Π_B diagnostic operator itself on integrable field, closed rational lines, an analytic helical island and a declared mixing component | Projection idempotence and mean preservation ≤1e−9 in the stated norms; target/invariant mismatch ≤1e−6 after a trace ladder; island-internal profile and whole-island-mixed cases distinguished; no projected-solve claim |
| G1c | Closure C manufactured forced curl problem for operator checks, then unforced Taylor eigenfields | Forced source residual ≤1e−10; physical unforced case independently meets Beltrami and D3 gates, with no manufactured force admitted as physical equilibrium |
| G1d | Closure D manufactured Φ/V/current-continuity system and equivalent voltage/current-driven solutions | Charge balance ≤1e−9, I_phi error ≤1e−6, matched V/I-driven fields agree to 1e−7, incompatible overconstraint rejected |
| G2 | Force-free Chandrasekhar–Kendall periodic cylinder, B built independently from Bessel/Helmholtz functions | Relative B error ≤1e−4 and expected order, relative curl B−μB ≤1e−6, D3 ≤1e−6; phase/amplitude fixed by declared symmetry/flux normalisations |
| G3 | Single-volume SPEC at the same wall, μ or helicity, flux and pressure constant, verifying only C | B/energy agreement ≤1e−3 after independent ladders, each solver meets its own residual gate; exclude unmatched regional invariants |
| G4 | Weakly forced resonant island in slab/cylinder, then large-aspect-ratio tokamak; analytic Hamiltonian width and independent outer-region solution | W within 1% or 1e−4 minor-radius units, whichever larger, phase within 0.02 rad; four amplitudes recover W proportional to sqrt(perturbation) with log slope 0.5±0.05 |
| G5 | Ideal shielding/current sheets from L15 and H22 versus reconnecting/island branch; separate closures and topology | Correct branch and current jump/outer field within 1% after convergence; do not demand a smooth pointwise current at the sheet |
| G6 | Nested limit versus VMEC++/GVEC/DESC through D7 and common-core residuals | Matched B,p RMS ≤1e−3, flux ≤1e−4, ι on surviving surfaces ≤1e−3; zero unintended modes ≤1e−10 for axisymmetric data |
| G7 | Island W, phase and p for A, B and C, with h/p/K/ratio and trace ladders | Finest two widths differ ≤1%, phases ≤0.02 rad, independent trace contribution ≤0.2% of W; normalised pressure error ≤1% where an oracle exists |
| G8 | Synthetic free boundary with non-axisymmetric Γ_c and full 3D exterior DtN coupling versus direct coil/vacuum and virtual-casing oracles | Field/trace error ≤1e−3, interface/Γ_c defect ≤1e−6, displaced seed returns within 1e−3 minor radii; independent geometry/trace/quadrature ladders and changing Γ_c in the vacuum collar keep observables within estimated errors |
| G9 | Matched nonlinear saturated tearing cylinder from T23, not an inferred “Huang JPP 2022” case | Freeze current/helicity/flux and closure comparison before running; W agreement ≤2% in the small-island regime, report large-island closure discrepancy rather than retuning |
| G10 | W7-X 5/5, 5/4, 5/6 and LHD n=1 RMP, after G0–G9 | All model/residual statuses retained, numerical W and pressure uncertainty ≤1%, phase ≤0.02 rad where defined, extended mode/mesh ladder required if a gate fails |

[V] For G4 an independent local field-line Hamiltonian
H(x,α)=a x²/2+ε cos α has separatrix full width W=4 sqrt(|ε/a|),
with canonical x,α and the physical radial conversion specified. This is a
known-perturbation oracle, not a prediction of nonlinear saturation.
Delta-prime is [∂r ψ_tilde]/ψ_tilde at the resonant surface from independent
outer solutions, with signs, wall and matching normalisations fixed.
Delta-prime alone does not determine a forced or saturated island width.
G9 supplies the additional current-relaxation/saturation assumptions.

[H] Pressure checks must reflect the closure: A tests finite-width transport
and a resolved anisotropy limit, B tests frozen fits and target defects, then
diagnostic invariant averages/mixing, and C tests constant regional p and traction.
The limiting ideal claim requires ε_parallel and the perpendicular Ohmic
force contribution each ≤1e−6, with decreasing D3 on spatial/mode ladders.
If these conditions cannot be reached, return a regularised candidate rather
than declaring an ideal equilibrium. Machine comparisons never replace G0–G9.

## Library boundaries and bounded milestones

[V] Generic numerical machinery stays in its owning library, with physics
composition in KIN6D. Reuse the pinned interfaces and check upstream before
introducing a workaround. Implementation PRs include generated leaves and
focused affected tests, while this design PR changes documentation only.

| owner | responsibility added only for an admitted milestone |
|---|---|
| [FortFEM](https://github.com/lazy-fortran/fortfem) | Commuting FEEC complexes, Piola/trace operators, partial assembly and Fourier testing, regional interface quadrature, generic anisotropic estimators and Γ_c/3D BEM actions |
| [FortSparse](https://github.com/lazy-fortran/fortsparse) | Sparse meridian block setup, factor reuse, gauge-aware block preconditioner and bordered Schur solves, memory/scaling instrumentation |
| [FortNum](https://github.com/lazy-fortran/fortnum) | FFT, Krylov orchestration where owned there, field-line ODE/root routines, Bessel functions, quadrature and verified arithmetic, without island physics policy |
| [FortSym](https://github.com/lazy-fortran/fortsym) | Generic symbolic differentiation and CPU/CUDA emission, primal/JVP/VJP consistency tools; no KIN6D-specific equilibrium model |
| KIN6D `gen/`, equilibrium and MHD layers | Pressure/current/Beltrami derivations and emitted leaves, closure admissibility and masks, physical targets/coil composition, shared D3–D8 adapters and topology observables |
| libneo adapter | Mature equilibrium interchange and conventions, upstream fixes for exercised defects, with no imported solver ontology imposed on native states |

[H] Each row below is one bounded PR with its own gate. Upstream generic
changes land in the owner first and the consumer pin is updated only after
rechecking the actual consumer. I-names are independent of equilibrium-core
D1–D8. All dependencies are explicit, and the controller chooses promotion.

| milestone | smallest complete change | dependency and gate |
|---|---|---|
| I1 | Closure declaration/admissibility and residual/report extension using existing A/D equations | Design approved, then G0/G1a/G1d; no new solver family |
| I2 | Independent CK and matched single-volume Taylor oracles and closure C | I1, G1c/G2/G3; regional interfaces beyond the single volume remain unavailable |
| I3 | Matrix-free generated JVP/VJP action over existing coupled unknowns | I1, assembled comparison on G0/G1, JVP finite-difference check ≤1e−7 and adjoint pairing ≤1e−10 |
| I4 | Mode-block preconditioner and FFT/tiled action | I3, G0–G3 unchanged; memory grows ≤2.5× when M doubles at fixed mesh, with timed 1/2/4/8-thread ladder |
| IH1 | High-order nested path with exact-Hessian Newton and upstream sum-factorised actions | Upstream prerequisites above and I3/I4 (I34); existing-ladder agreement with P2/P3 within estimated error, G6, N1 matched-accuracy report and ≤2× VMEC++ DOF gate closed; operator cost per element versus p reported |
| IH2 | Construct local field-aligned/curve-following charts from a VMEC++ or KIN6D nested solution | IH1, I1 and I3/I4; positive chart Jacobians, commuting/Piola and conservative transfer checks, G0/G6 unchanged within estimated error; charts remain regular through the G4 separatrix while the solved field may change topology |
| IH3 | hp-adaptive non-nested path with a refinement indicator from I6 island diagnostics/error estimator, high p in smooth regions and local h-refinement at islands/separatrices/current sheets | IH2, I1, I3/I4 and I6; G4/G5/G7 and fixed-DOF pollution gate above; report DOFs inside/outside the marked region, with unmarked-region DOFs ≤1.25× the nested seed at matched accuracy [H] |
| I5a | B via frozen fitted S and conductivity/current drive on the reference label, solving A+D | I1/I3, G1a/G1d and seed-target reproduction on G4; expose seed-fit and solved-target defects |
| I5b | Π_B diagnostic and explicit whole-island mixing diagnostic | I5a, G1b on the projection operator itself and G4; expose trace, target and current-solvability defects |
| I5c | Optional projected solve mode | I5b, G1b/G4 plus topology-change and derivative-status checks; not a prerequisite of machine work |
| I6 | Width/phase/G estimator integration and resonant/current-sheet cases | I2/I4/I5a/I5b, G4/G5/G7 and independent nested-limit G6; practical error labels retained |
| I7 | Free-B Γ_c coupling through D6 with full non-axisymmetric 3D BEM and exterior DtN from the start, reusing PR #47 (`bem3d5`) | I4/I6, M4 and accepted FortFEM/consumer adapter with coincident-pair quadrature resolved, G8 plus coil and casing checks; optional large-Ω prescribed-trace oracle, no duplicate exterior solver |
| I8 | Marked regional Taylor interface coupled to A/B exterior pressure closure | I2/I5a/I7, interface MMS and traction/jump gates of G5/G8, helicity/flux compatibility checked |
| I9 | Small-island nonlinear saturation comparison and machine-size synthetic memory trial | I6–I8, G9 plus measured action/preconditioner memory and cost; no production machine claim |
| I10 | W7-X vacuum 5/5 EIM, then 5/4, 5/6, edge-island candidates, then finite beta | I9 and prior required gates including I7 full 3D BEM/DtN, G10 per case; public currents, synthetic/private filaments outside public repo |
| I11 | Full-torus LHD helix plus n=1 RMP and sidebands | I9 and prior required gates plus calibrated RMP coil set owned by LHD lane, G10 with K≥24,32,48 extended as needed; no symmetry reduction to nfp=10 |

[H] Open mathematical risks are the smooth ideal limit at dense resonances,
closed-line/chaotic averaging without a uniform finite trace bound, and loss
of invertibility near bifurcations. Open computational risks are anisotropy-
uniform preconditioning, Fourier aliasing, factor memory and scalable BEM.
They are reported as failed or unavailable gates, not removed by changing
tolerances. No W7-X/LHD result, continuous certificate or speed advantage is
established by this document. Controller design acceptance does not promote
an implementation or machine campaign without its required gates.

[H89]: https://doi.org/10.1016/0021-9991(89)90069-7
[H06]: https://doi.org/10.1088/0029-5515/46/11/L01
[H17]: https://doi.org/10.1088/1361-6587/aa5adc
[S11]: https://doi.org/10.1063/1.3597155
[P86]: https://doi.org/10.1016/0010-4655(86)90059-7
[P89]: https://doi.org/10.1016/0021-9991(89)90066-1
[M06]: https://doi.org/10.1017/S0022377806005861
[M12]: https://doi.org/10.1063/1.4765691
[L15]: https://doi.org/10.1063/1.4906888
[L16]: https://doi.org/10.1063/1.4967709
[Q21]: https://doi.org/10.1088/1361-6587/ac2afc
[F20]: https://doi.org/10.1088/1361-6587/ab9a61
[SP26]: https://doi.org/10.48550/arXiv.2607.27135
[E05]: https://doi.org/10.1088/0029-5515/45/7/022
[B17]: https://doi.org/10.1109/TPS.2017.2780624
[J17]: https://doi.org/10.1088/0029-5515/57/2/022013
[C09]: https://doi.org/10.1063/1.3224035
[N04]: https://doi.org/10.1016/j.jcp.2003.10.004
[R25]: https://doi.org/10.48550/arXiv.2510.26986
[A26]: https://doi.org/10.1137/25M1727540
[K17]: https://doi.org/10.1063/1.4986493
[H22]: https://doi.org/10.1063/5.0067898
[T23]: https://doi.org/10.1017/S0022377823000934
[D23]: https://doi.org/10.1017/S0022377823000272
