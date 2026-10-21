# Common equilibrium core

Status: decision proposal for review, 2026-10-13; no code change. Refines
[R12](../ROADMAP.md) and [equilibria-inference.md](equilibria-inference.md).
Sources outside the repo: 3D TC24 design (§2, §7), free-boundary GS design
(P1–P13), GVEC FB-cut diagnosis (2026-10-13).

## Problem

KIN6D has 2D ψ Grad–Shafranov (fixed, prescribed-q, free boundary with coils
and inverse; CHEASE/FreeGS/CEDRES++ class), the frozen-in map lane N1
(`ideal_map`, branch `gs-3d-nonlinear-n1`; VMEC/GVEC/DESC class), the planned
per-n linear response lane L3–L8 (MARS-F/GPEC class), and a planned
non-nested Eulerian branch (MRX/HINT class). If each keeps its own state,
residual, boundary coupling, readback and tests, KIN6D becomes four separate
codes. Decision: they are constrained parametrisations of one problem; each
lane keeps its unknowns and seed, and shares the pieces D1–D8.

## D1. One physical state

The state is 𝒮 = (B, p, Λ, 𝒞, π):

- B: closed 2-form in physical space, held as parametric H(div) data on a
  (possibly deformed) mesh (Piola-mapped RT ⊕ DG per toroidal mode n), or as
  map plus reference form. Evaluation exposes Cartesian components
  ([geometry.md](geometry.md)).
- p: pressure in L², optionally p(label).
- Λ: label field(s) where they exist, profile laws (`gs_profile_law_t`),
  invariants Φ_tor, ψ_pol, I_p, F_vac. May be empty in the non-nested lane.
- 𝒞: chart X = H(q) with jet H, dH, d²H, d³H, i.e. `ideal_chart_t` lifted
  unchanged out of `ideal_map`. Cylindrical first; the curve frame (G-frame,
  Figure-8) is a later chart and changes nothing else.
- π: provenance (lane, seed, continuation history, claim label).

The core stores B, never a vector potential; its harmonic content (toroidal
flux, poloidal circulation) is data in Λ. Lane gauges stay inside lanes (D2).
The absolute poloidal flux is defined from B, ψ(R,Z) = (1/2π)∫_disk B·dS =
R A_φ, with ψ = 0 on the symmetry axis (already the `gs_free` gauge). This
removes the "ψ = 0 on the mesh boundary" assumption of `gs_readback` and
`gs_coarea`.

## D2. Parametrisations are reduced spaces

Each lane is a map P: U_lane → 𝒮 with a declared constraint set.

| Lane | Unknown U | B = P(U) | Constraint set | Gauge | Seed |
|---|---|---|---|---|---|
| ψ (2D) | nodal ψ, ψ_b, ψ_a | ∇ψ×∇φ + F(ψ)∇φ | axisymmetry; p, F as functions of ψ (or q prescribed) | ψ = 0 on the symmetry axis (free boundary) or on Γ (fixed boundary) | analytic or previous ψ |
| map (nested) | X: Ω̂ → R³, P_k × Fourier | DX b̂ / det DX | topology, fluxes and p per label frozen | parallel relabelling (fixed transversal + axis expansion), ω_0 pin | nested seed from boundary + profiles, never the target |
| free B (non-nested) | B ∈ H(div) with div B = 0 (or A ∈ H(curl) internally) | identity | only div B = 0, B·∇p = 0 by a declared p-model | A-gauge internal, harmonic part from Λ | 𝒮 from ψ or map lane |

The constraint sets differ. Only states satisfying all of them (the
axisymmetric, nested, smooth ideal state) are fixed points of every lane.
Elsewhere the lanes give different equilibria, stored separately
([equilibria-inference.md](equilibria-inference.md)). A lane is admitted with
a P and a declared constraint set, never with its own residual.

## D3. One residual: weak Maxwell-stress form

\[
\mathcal R(B,p)[\xi] = \int_\Omega \Big( \tfrac{1}{\mu_0}B\cdot\nabla\xi\cdot B
 - \big(\tfrac{|B|^2}{2\mu_0}+p\big)\nabla\cdot\xi \Big)\,dV
 + \text{boundary terms},
\]

equal to −⟨J×B − ∇p, ξ⟩ for smooth fields.
- It needs no derivatives of B or p, so it is defined for the ψ lane (P_k ψ
  gives tangential B jumps and face currents), the map lane and class (b)
  states with current sheets. A strong J×B − ∇p is not.
- For the map lane it is exactly −δW[ξ] (TC24 §7.1): Newton gradient and
  residual are one object. For a smooth ψ-lane state it is the GS residual
  tested with v = ξ·∇ψ/(μ0R²).
- Scale-aware normalisation (GVEC FB-cut lesson): with B_0 the prescribed
  curl-free background (e.g. TF F_vac/R ∇φ), drop ∇·T(B_0) = 0 analytically
  and evaluate only cross and quadratic terms in b = B − B_0. Normalise by
  σ = ⟨|b|²⟩/μ0 + ⟨p⟩ (B_p²-scale in a tokamak), never by ⟨B²⟩. Report ‖𝓡‖ in
  the dual norm of the displacement space with σ. Gate edge F consistency
  |F_edge − F_vac|/F_vac ≤ 1e-6 separately.
- Plasma-boundary position is never a [[B²]] relaxation; it is set through
  the Γ_c coupling (D6).

## D4. One certificate interface

Generalise `gs_majorant` into the interface ‖U − U_h‖ ≲ 𝓜(U_h)/β_h. 𝓜 is a
computable majorant of 𝓡 in the dual norm (ψ lane: the existing Repin
D + C_F B plus the P10 coupling term N). β_h is an inf-sup estimate of the
lane Jacobian on the gauge complement, replacing coercivity α where the
problem is indefinite (free boundary, vertical mode). Lanes implement 𝓜 and
β_h. The core owns the claim labels (candidate,
residual verified, practical estimate, certified); a lane without β_h reports
"residual verified" at most. ψ-lane statuses 4 and 5 are unchanged.

## D5. One operator chain

The linearisation of 𝓡 (second variation of W) at an axisymmetric state, per
toroidal mode n, is the ideal δW_n. It is assembled once (L4) from the
generated K-E1 leaves and used as (1) the per-n block preconditioner of the
map-lane Newton, exact at the reference; (2) the regular-space linear ideal
response operator; (3) the base of the closures L5 (resistive) and L6 (ideal,
broken + Frobenius) for MARS-F/GPEC-class response.

[V] The first local D5 increment is described in [ideal-xi.md](ideal-xi.md):
the gamma=0 plasma-volume xi form, exact FortSym Hessian identities and
per-n cylindrical displacement jets. Global FE assembly, gauge reduction and
boundary coupling remain next increments.

Required code tests: Hess W |_{X=id} on mode n equals L4_n to roundoff
(TC24 §7.4). At N = 0 the map lane is prescribed-q GS (the generator proves
the pointwise reduction), so at a solution its reduced Hessian equals the
linearised prescribed-q GS operator pulled back by δψ = −ξ·∇ψ on the gauge
complement. Forward GS is not the partner: p, F fixed as functions of the
moving flux add profile-response terms.

## D6. One boundary layer

- A fixed axisymmetric surface Γ_c in vacuum encloses the plasma, so the
  exterior separates in n. It carries the per-n Steklov–Poincaré map S_n of
  the decaying exterior field (F-n2/F-n3: kernel Q_{|n|−1/2}, symmetric
  Costabel coupling) acting on the tangential trace of the interior field, and
  the coil drive g_{k,n}.
- One coil layer: `gs_coils` rectangles and `coils_3d` polylines share an
  interface returning B_coil and its per-n traces on Γ_c.
- The ψ coupling is the scalar n = 0 instance; by the F-n2 correction
  ψ = R A_φ is the n = 1 vector mode of the kernel. `gs_coupling` becomes an
  adapter.
- Every lane exposes the field in a vacuum shell between Γ_p and Γ_c. The
  free-boundary Newton acts on the Γ_c trace with exact coil Green's
  functions; the plasma position follows from flux matching (FreeGS and ψ-lane
  conditioning), not from [[B²]] targets.
- Γ_p conditions are lane-specific: implicit in the ψ lane (weak continuity,
  cut-cell interface term); natural traction [[p + B²/2μ0]] = 0 with B·n = 0
  on a material surface in the map lane; absent in the free-B lane, where
  vacuum is the region with p = 0 and J = 0.
- The coil inverse acts on functionals of 𝒮 (disk flux, point B), not nodal
  ψ; its Gauss–Newton, BVLS and GSVD parts are lane-independent.

## D7. One readback and export layer

Every consumer reads 𝒮 through one evaluator: batched point values of B, p,
J, ∇B (as `sample_gs_many`); the magnetic axis as a fixed point of the
poloidal field-line map; label level sets where Λ has a label; q, Φ_tor and
volume from flux integrals of B over them; SFL angles (PEST, Hamada, Boozer);
exports (G-EQDSK for axisymmetric states with ψ from D1, Boozer `.bc` through
libneo, L7 spectra).

Kinetic geometry view (`src/equilibrium/equilibrium_geometry.f90`, axisymmetric states):
on a surface psi = psi_0 of the evaluator it returns B, |B|, grad|B| (five-point
stencil of the evaluator), grad psi, b.grad theta, (b x grad psi).grad theta,
the field-line label alpha = phi - nu with grad alpha and d alpha/d psi (shear),
q, q', dV/dpsi, flux-surface averages (`surface_average`) and the grad-B drift
factors (b x grad|B|).grad psi/|B|^2 and (b x grad|B|).grad alpha/|B|^2, at
geometric, PEST or Boozer angles (unwrapped for field lines). Boozer is
theta_B with Jacobian ~ 1/B^2 and zeta_B = phi + q (theta_B - theta*); the
covariant components F, I are not exported. Assumes surfaces star-shaped about
the axis; curvature drift at finite beta (needs p') is not included.
Tests: `equilibrium_core_d7_geometry_{solovev,circle,drift,state}`.

In the map lane q per label is an input (frozen fluxes); readback recomputes
it from B in physical space as the cross-lane check and never copies the
reference. `vacuum_surfaces` (ψ grid from G-EQDSK) becomes an import adapter
that builds a core state, not a second readback.

## D8. One verification harness

One case registry, `run(case, lane, resolution) -> 𝒮 + report`, reports per
lane:

1. Oracle error through D7 (Solovev A3, shaped Solovev, Cerfon–Freidberg,
   cylinder, Landreman E0, exact vacuum coil fields, FreeGS FB-cut parity) and
   the σ-normalised D3 residual.
2. Reductions: map lane at N = 0 against prescribed-q GS (existing N1
   ladders); no leakage for N ≥ 1 with axisymmetric data; free-B lane seeded
   from an axisymmetric state stays there.
3. Offset-start return: displaced-axis, -boundary or -shape seeds (15 mm,
   2 % elongation) must return; reference stationarity alone is not
   acceptance (GVEC lesson). N1's interior-shift seed is the first case.
4. Cross-lane: 𝒮 of lane A is evaluated by D3/D7 in lane B's discretisation
   and seeds lane B, which must stay stationary to discretisation accuracy.
5. The D5 identities.

Registered first-stage free-boundary case: `ideal_fb_landreman_axisymmetric`
([data and evidence](ideal-free-boundary.md#landreman-axisymmetric-benchmark)).
It checks Landreman E0 through the N = 0 map solve, source-fit vacuum drive,
S0 distance and field/pressure refinement, and solves the same ring-coil
drive through `gs_free`. Offset recovery and a shared registry runner remain open; the closed-line family supplies
stationarity evidence only.

Solvers are Newton(–Krylov) with exact derivatives and shared globalisation
(Armijo or trust region, quadrature tightening, det DX > 0 guard); heuristic
gradient descent is not an acceptance path.

M6's [D8 harness and case registry](equilibrium-harness.md) routes the retained
map and manufactured psi ladders through one status/report contract. It adds
psi return ladders for +/-15 mm axis offsets and 2% elongation. Reports name
their residual norms and unavailable state/derivative bounds explicitly;
this does not establish all cross-lane and free-B gates listed above.

## Module classification

The line counts are from `origin/main` (and N1) on 2026-10-13. The
shared/lane split is a rough estimate of the code after the lift, ±30 %.

| Module | Lines | Class | Shared after lift | Lane-specific |
|---|---:|---|---:|---|
| `grad_shafranov` | 2325 | split: mesh, curved cells, jets, axis search, sampling → core; Δ* operator, profile loads, response cache → ψ | ~1000 | ~1300 ψ |
| `gs_nonlinear` | 1166 | ψ lane; Newton globalisation lifted to the shared driver; level-cut quadrature stays generic in FortFEM | ~200 | ~950 ψ |
| `gs_profile_law` | 242 | core as is (profile data of Λ) | 240 | – |
| `gs_readback` | 221 | lift to B-based (D7); drop the ψ = 0 boundary assumption | ~200 | – |
| `gs_coarea` | 863 | lift: label level-set integration of any scalar label; ψ is one label | ~650 | ~200 ψ |
| `gs_field_readback` | 186 | lift: batched evaluator of 𝒮 | 186 | – |
| `gs_axis_recovery` | 98 | lift into D7 | 98 | – |
| `gs_majorant` | 460 | ψ-lane implementation of the D4 interface | ~60 | ~400 ψ |
| `gs_free_boundary` | 2030 | ψ lane (LCFS modes, cut at ψ_b, interface mass, ψ-level refit); Newton/factor driver → shared | ~400 | ~1600 ψ |
| `gs_coupling` | 575 | lift to per-n Γ_c layer (D6); ψ adapter remains | ~400 | ~150 ψ |
| `gs_coils` | 651 | core as is (n = 0 coil set); merged interface with `coils_3d` | ~600 | – |
| `gs_coil_inverse` | 979 | lift: rows as functionals of 𝒮; GN/BVLS/GSVD generic | ~650 | ~300 ψ |
| `gs_prescribed_q` | 318 | ψ lane | – | 318 ψ |
| `gs_inverse` | 447 | ψ lane (constant-q entry point) | – | 447 ψ |
| `gs_csr`, `gs_direct`, `gs_options`, `gs_timers` | 462 | core infrastructure as is (rename later) | 462 | – |
| `coils_3d` | 318 | core as is (boundary layer) | 318 | – |
| `vacuum_surfaces` | 446 | lift: import adapter + SFL harmonics on D7 | ~350 | ~100 import |
| `ideal_map` (N1) | 973 | map lane; `ideal_chart_t` → core chart; Hessian assembly → D5 | ~150 | ~820 map |

Sum: 12.8k lines today; about 6.0k shared, 5.7k ψ-specific, 0.8k
map-specific. The core is mostly existing `gs_*` code behind ψ-specific
interfaces, not new code.

Library boundaries are unchanged: FortFEM owns FEEC spaces, Piola maps
(F-n1), commuting projections, level-cut quadrature and per-n BEM (F-n2/F-n3);
FortNum Q_{n−1/2}, elliptic integrals, BVLS, GSVD and Krylov; FortSparse the
factorisations; libneo G-EQDSK/Boozer I/O. Energy and residual leaves are
physics and are generated by fortsym into KIN6D `gen/`.

## Migration order (K6 stays reproducible)

The public CTest `gs_k6_mini_regression` now uses an independent synthetic
Solov'ev fixture, documented in [its recipe](../test/fixtures/gs_synthetic/README.md).
The historical migration gates below describe the original private K6 campaign;
its machine-derived inputs and reference outputs stay in the private TC24
repository. Public regression results do not establish agreement with K6.

K6, the accepted TC24 phase-4b state, runs `gs_phase1` → cubic-profile
Newton (`gs_nonlinear`, P3, direct factor) → `gs_readback`, `gs_coarea`,
`gs_majorant`. No step changes the numerics of `grad_shafranov` or
`gs_nonlinear`; only interfaces and callers move. M0 adds a CTest K6-mini case
(same law and path, coarse mesh) that must reproduce `gs_phase1` outputs to
≤ 1e-12 relative at every step. The full K6 rerun (≈ 100 s, 48 threads) is
repeated at M2 and M4 and must match the accepted Bpol, q and ψ figures to the
last reported digit.

| Step | Content | Gate |
|---|---|---|
| M0 | `equilibrium_state_t`, chart module (`ideal_chart_t` moved), evaluator interface; K6-mini regression test | builds; K6-mini bitwise |
| M1 | Adapters ψ → 𝒮 and map → 𝒮 (read-only mesh view of `gs_state_t`, no type move); D3 residual evaluator with σ normalisation | D3 of exact Solovev → 0 at the expected order in both lanes; K6-mini |
| M2 | D7 readback from 𝒮 (axis, label sets, q, Φ_tor, volume, ψ from disk flux); `gs_readback` becomes a thin caller | equals existing readback on fixed-boundary meshes to ≤ 1e-12; `gs_free` gains q; full K6 rerun |
| M3 | D4 interface; `gs_majorant` registers as ψ implementation; claim labels in core | status 4/5 unchanged on existing tests |
| M4 | D6 per-n Γ_c layer; `gs_coupling` → adapter; unified coil set | P6–P9 tests and `gs_free_app_contract` unchanged; FreeGS parity numbers unchanged; full K6 rerun |
| M5 | D5: L4 assembled in core from K-E1 leaves; identity tests | Hess W = L4_n to roundoff; N = 0 map Hessian = prescribed-q Jacobian |
| M6 | D8 harness and case registry; offset-start return tests for the ψ free-boundary and map lanes | existing ladders reproduced through the harness |
| M7 | Physical move of the mesh container out of `gs_state_t` (composition) | K6-mini bitwise; full suite |

Implemented: M0 `equilibrium_state` (𝒮, evaluator `equilibrium_field_t`),
`equilibrium_chart` (`ideal_chart_t`, cylindrical; N1 still has its own
copy), CTest `gs_k6_mini_regression`. M1 `gs_equilibrium` (ψ → 𝒮 over a
geometry view), `equilibrium_residual` (D3, fortsym-proved stress leaf);
`map_equilibrium` now adapts N = 0 cylindrical maps, retaining the reference
form and represented geometry without solver factors. It recomputes physical B
with the generated Piola leaf (including omega shear), checks the supplied
reference against the frozen quadrature data, and reports a candidate. Curved
charts and nonzero Fourier maps return unavailable; their D3 displacement
spaces and D7 contracts still need extension. P2/P3 deformed Solovev tests
check D3 convergence, independent physical fields, generated JVP/VJP and
OpenMP parity/scaling. Continuous field bounds remain unavailable. M2 `equilibrium_readback` on
`labelled_field_t`; `read_gs_equilibrium` is a thin caller and takes
`psi_bound` and `plasma_cell` for free-boundary states. The N=0 map adapter
provides label primitives on the Pk projection of its frozen label, physical
volume/flux quadrature, and contour-integrated toroidal flux density. The latter
uses the physical B rather than assuming R B_phi is constant on a candidate's
label. Label projection and quadrature errors are not enclosed; readback reports
`candidate` and `bound_status = unavailable`, retaining the source claim
separately. The current contour calculation requires labels star-shaped about
the axis and rejects unresolved intersections. M3
`equilibrium_certificate` (`equilibrium_certifier_t`: a lane returns 𝓜 and
β_h with its kind; `certify` forms 𝓜/β_h and the claim label, written to π);
`gs_certificate` wraps `gs_functional_majorant` (status 4 → residual
verified, 5 → practical estimate, 0 → certified) and `gs_fb_majorant`
(practical estimate) without changing them. M4 `equilibrium_boundary`
(D6): Γ_c chain from a mesh boundary, P2 panels and trace restriction,
per-(kind, n) cached Steklov–Poincaré map (`GAMMA_C_SCALAR`: L_n of the
scalar potential, any n; `GAMMA_C_FLUX`: ψ, n = 0) on the panel and FEM
trace, quadrature and trace loads, one coil interface
(`axisymmetric_coils_t` over `gs_coils`, `polyline_coils_t` over
`coils_3d`) and per-n coil traces on Γ_c (Φ_n and ∫ R B·n N_j ds).
`gs_coupling` is the n = 0 adapter (same chain, panels, S and coil load
arithmetic). CTest `equilibrium_core_d6_*`: S_n u + h = 0 for a tilted
loop inside Γ_c, n = 0..3, third order; ring vs polygon coils at n = 0.

### Which upcoming steps wait for what

| Step | May proceed now | Must wait for |
|---|---|---|
| P12 (ITER FB producer, iter_tc24) | yes; pins a kin6d commit, ψ lane frozen | M2 only for the q output |
| N2 (NL-1 axisymmetric ladder, 3D mode) | solver work | M1 + M2: it must report through the core residual and readback, not ψ-style readback |
| N3 (NL-2, Landreman/VMEC/DESC/GVEC) | – | N2, M6 (offset-start return in the harness) |
| N4 (3D free boundary) | – | M4; no map-lane-specific vacuum; D3 normalisation and return tests (M6) |
| L4 | – | M0 (chart, state); it *is* M5 |
| L5, L6 | – | L4 only |
| L7 (spectra, Boozer, exports) | – | M2; built on D7, not on `vacuum_surfaces` ψ grids |
| L8 | – | L7 |
| Non-nested branch | [NN1](equilibrium-nonnested.md): prescribed-trace vacuum and the N=0 constant-p' restriction | General finite-pressure/free-boundary work still needs M3, M4 and M6; seeds come from 𝒮 and the p-model is declared first |

L3 (vacuum response) needs no core piece and continues. Its Γ_c/BEM_n and
coil code is written as the D6 layer, which M4 then reuses for `gs_coupling`.

## What is and is not novel

Not novel: GS FEM-BEM coupling on a fixed curve (CEDRES++, FEEQS), the
frozen-in energy over maps (Kruskal–Kulsrud, VMEC, GVEC, DESC), per-n toroidal
Green's functions (VACUUM), δW (Bernstein et al.), FEEC magnetic relaxation in
physical space (MRX), the Maxwell-stress weak form, Newton–Krylov.

New in combination, as far as we know (literature check required before any
claim; no performance claim against MRX or VMEC++/GVEC/DESC without
time-to-accuracy measurements):

1. One residual and certificate interface across ψ, nested-map and Eulerian
   lanes, with D8 cross-lane checks; the families above certify, if at all,
   only within their own representation.
2. A nested-map lane on an unstructured curved FE reference mesh with an
   abstract chart and no flux coordinates in the solver (VMEC: flux surfaces;
   GVEC: radial B-spline; DESC: Zernike). Benefit near separatrix-adjacent
   boundaries is unproven.
3. The discrete identity between nonlinear map Hessian and per-n linear ideal
   operator, used as preconditioner and consistency test.
4. 3D free boundary through a fixed Γ_c FEM shell with per-n exact exterior
   instead of virtual casing on the moving plasma boundary (NESTOR, BIEST).
   The conditioning advantage is a hypothesis until NL-3a passes
   offset-start return tests.

## Open mathematical questions

1. **Push-forward exactness.** B = DX b̂/det DX is exactly divergence-free
   cellwise and has continuous normal flux. With b̂ = rot χ̂ + b̂_3 and X
   polynomial per cell, it lies in parametric (Piola-mapped) RT on the
   deformed mesh. Open: (a) transfer to a fixed Eulerian mesh (free-B seed,
   comparison) is exact only with a conservative 2-form remap (supermesh) or a
   commuting projection on the deformed geometry, cost and accuracy
   unmeasured; (b) Fourier truncation of B itself breaks exactness, so export
   carries the map and b̂, not truncated B modes.
2. **Gauges.** Map lane: parallel relabelling; ω_0 = f(χ) is a continuum null
   direction at n = 0 that N1 leaves to damping. An unbiased gauge keeping P_k
   order is not settled (candidate: flux-aligned reference mesh K-M1). Free-B
   lane: B with a div multiplier versus A ∈ H(curl) with tree–cotree gauge is
   undecided; the harmonic data in Λ is common.
3. **Certificate for the map lane.** W is non-convex; δW has continuous
   spectrum and, with sheared ι, no uniform inf-sup in H¹ at rationals, so β_h
   may tend to 0 with h. Expected: "residual verified", with error bounds only
   away from rationals. Which displacement norm makes the D3 dual norm stable
   is open.
4. **ψ lane against ideal constraints.** Forward GS fixes p(ψ), F(ψ); the
   map lane fixes mass and flux per label, hence D5 pairs it with prescribed-q
   GS. Whether a current-constrained f(χ̂) map lane reproduces the forward
   ψ-lane Jacobian exactly is open.
5. **Pressure in the free-B lane.** B·∇p = 0 makes the Eulerian Jacobian
   singular on rational and ergodic sets. Anisotropic diffusion (χ_∥/χ_⊥),
   MRX-type relaxation and resistive closure are model choices (R12); the core
   does not decide.
6. **3D free-boundary conditioning.** That shell + Γ_c keeps the ψ lane's
   position conditioning when Γ_p is a material surface of the map is a
   hypothesis; NL-3a with offset starts decides it.
   Gauge, saddle and convergence obligations: [free-boundary-theory.md](free-boundary-theory.md).

Chris&AI

The [production scaling ledger](benchmarks.md#production-estimator-scaling)
times D3 on a solved GS field. Profile conversion, stress-load conversion
and estimator geometry hulls use OpenMP. The estimator cost is the paired
difference between calls with and without `error_estimate`; it includes both
Riesz-component reconstructions and additional field sampling.
