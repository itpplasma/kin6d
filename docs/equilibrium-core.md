# Common equilibrium core

Status: architecture decision proposal, 2026-10-13, for review. It changes no
code. It refines R12 ([ROADMAP](../ROADMAP.md#r12--equilibrium-mhdkinetic-mhd-vacuum-and-circuits))
and [equilibria-inference.md](equilibria-inference.md) for the static ideal
equilibrium family. The source designs are the 3D TC24 design (linear lane §2,
nonlinear lane §7), the free-boundary GS design (P1–P13) and the GVEC
free-boundary diagnosis of 2026-10-13. All three sit outside the repository.

## Problem

KIN6D now has three ways to represent an equilibrium, and a fourth is planned:

- 2D ψ Grad–Shafranov: fixed boundary, prescribed q, free boundary with coils
  and inverse. This is the CHEASE/FreeGS/CEDRES++ class.
- The frozen-in map lane N1 (`ideal_map`, branch `gs-3d-nonlinear-n1`). This
  is the VMEC/GVEC/DESC class.
- The per-n linear response lane (L3–L8). This is the MARS-F/GPEC class.
- A non-nested Eulerian branch, still to come. This is the MRX/HINT class.

If each one keeps its own state, residual, boundary coupling, readback and
tests, KIN6D becomes four separate codes. This note makes them constrained
parametrisations of one problem. Each lane keeps its own unknowns and seed.
They share the physical state, the residual, the certificate interface, the
linearised operator, the vacuum/coil boundary layer, readback/export and the
verification harness.

## D1. One physical state

\[
\mathcal S = (B, p, \Lambda, \mathcal C, \pi).
\]

- B is a closed 2-form in physical space. For export it is held as parametric
  H(div) data on a (possibly deformed) mesh: Piola-mapped RT ⊕ DG per toroidal
  mode n, or the map together with its reference form. Evaluation exposes
  Cartesian components, as [geometry.md](geometry.md) requires.
- p is the pressure in L², optionally as p(label).
- Λ is the label and profile data: label field(s) where they exist, profile
  laws (`gs_profile_law_t`), and the invariants Φ_tor, ψ_pol, I_p and F_vac.
  Λ may be empty in the non-nested lane.
- 𝒞 is the chart, X = H(q), with the jet H, dH, d²H, d³H. It is
  `ideal_chart_t`, lifted unchanged out of `ideal_map`. Cylindrical is the
  first chart. The curve frame (G-frame, Figure-8) is a later chart and
  changes nothing else.
- π is the provenance: lane, seed, continuation history and claim label.

The core stores B, never a vector potential. Its harmonic content (toroidal
flux and poloidal circulation for a toroidal shell) is data in Λ. That is the
only gauge the core knows. Each lane's gauge stays inside the lane (D2).

The absolute poloidal flux is defined from B as
ψ(R,Z) = (1/2π) ∫_disk B·dS = R A_φ, with ψ = 0 on the symmetry axis. This
removes the "ψ = 0 on the mesh boundary" assumption of `gs_readback` and
`gs_coarea`. It is already the gauge of `gs_free`.

## D2. Parametrisations are reduced spaces

Each lane is a map P: U_lane → 𝒮 together with a constraint set.

| Lane | Unknown U | B = P(U) | Constraint set | Gauge | Seed |
|---|---|---|---|---|---|
| ψ (2D) | nodal ψ, ψ_b, ψ_a | ∇ψ×∇φ + F(ψ)∇φ | axisymmetry; p, F as functions of ψ (or q prescribed) | ψ = 0 on axis of symmetry (free) or on Γ (fixed) | analytic or previous ψ |
| map (nested) | X: Ω̂ → R³, P_k × Fourier | DX b̂ / det DX | topology, fluxes and p per label frozen | parallel relabelling (fixed transversal + axis expansion), ω_0 pin | nested seed from boundary + profiles, never the target |
| free B (non-nested) | B ∈ H(div) with div B = 0 (or A ∈ H(curl) internally) | identity | only div B = 0, B·∇p = 0 by a declared p-model | A-gauge internal, harmonic part from Λ | 𝒮 from ψ or map lane |

The constraint sets differ. Only states that satisfy every lane's constraints
are fixed points of all lanes; the axisymmetric, nested, smooth ideal state is
the common one. Away from it the lanes give different equilibria (D7), and
these are stored separately as [equilibria-inference.md](equilibria-inference.md)
requires. A lane is admitted only with a P and a declared constraint set. The
lane does not get its own residual.

## D3. One residual: weak Maxwell-stress form

\[
\mathcal R(B,p)[\xi] = \int_\Omega \Big( \tfrac{1}{\mu_0}B\cdot\nabla\xi\cdot B
 - \big(\tfrac{|B|^2}{2\mu_0}+p\big)\nabla\cdot\xi \Big)\,dV
 + \text{boundary terms},
\]

which equals −⟨J×B − ∇p, ξ⟩ for smooth fields.

- It needs only B, p ∈ L². It is therefore defined for the ψ lane, where P_k ψ
  has tangential B jumps and face currents, for the map lane, and for class (b)
  states with current sheets. A strong J×B − ∇p is not.
- For the map lane it is exactly −δW[ξ] (TC24 §7.1), so the Newton gradient
  and the residual are the same object. For the ψ lane, restricted to
  axisymmetric ξ, it is the weak GS residual.
- Scale-aware normalisation (lesson of the GVEC FB-cut diagnosis).
  - Let B_0 be the prescribed curl-free background, for example the TF field
    F_vac/R ∇φ. Remove ∇·T(B_0) = 0 analytically and evaluate only the cross
    and quadratic terms in b = B − B_0.
  - Normalise by σ = ⟨|b|²⟩/μ0 + ⟨p⟩, which is B_p²-scale in a tokamak, and
    never by ⟨B²⟩.
  - Report ‖𝓡‖ in the dual norm of the displacement space together with σ.
    Separately, gate edge F consistency |F_edge − F_vac|/F_vac ≤ 1e-6.
- Position residuals on the plasma boundary are never [[B²]] relaxations. The
  position is set through the Γ_c coupling (D6).

## D4. One certificate interface

Generalise `gs_majorant` into the interface

\[
\|U-U_h\| \lesssim \mathcal M(U_h)/\beta_h ,
\]

- 𝓜 is a computable majorant of 𝓡 in the dual norm. For the ψ lane it is the
  existing Repin D + C_F B, extended by the coupling term N of P10.
- β_h is an inf-sup estimate of the lane Jacobian on the gauge complement,
  replacing the coercivity α where the problem is indefinite (free boundary,
  vertical mode).

Lanes implement 𝓜 and β_h. The core owns the claim labels, from candidate to
residual verified to practical estimate to certified, and the rule that a lane
without β_h reports "residual verified" at most. The ψ lane keeps statuses 4
and 5 unchanged.

## D5. One operator chain

The second variation of 𝓡 at an axisymmetric state, per toroidal mode n, is
the ideal δW_n. It is assembled once (L4) from the generated K-E1 leaves and
used for three things:

1. the Newton preconditioner of the map lane (per-n block-diagonal, exact at
   the reference);
2. the linear ideal response operator, regular space without closure;
3. through the closures L5 (resistive) and L6 (ideal, broken + Frobenius), the
   MARS-F/GPEC-class response.

Identities required as code tests:

- Hess W |_{X=id}, restricted to mode n, equals L4_n to roundoff (TC24 §7.4).
- At N = 0 the map lane is the prescribed-q GS problem. The generator already
  proves the pointwise reduction. Its reduced Hessian must therefore equal the
  Jacobian of the prescribed-q ψ lane on the gauge complement. This test is
  sharper than a comparison with the forward ψ lane, whose Jacobian belongs to
  a different constraint set (p, F fixed as functions of the moving flux) and
  differs from δW_0 by the profile-response terms.

## D6. One boundary layer

- A fixed artificial surface Γ_c lies in vacuum and encloses the plasma. Γ_c
  is axisymmetric, so the exterior separates in n.
- It carries the per-n Steklov–Poincaré map S_n of the exterior decaying field
  (F-n2/F-n3: kernel Q_{|n|−1/2}, symmetric Costabel coupling), coupled to the
  tangential trace of the interior field.
- The coil drive is g_{k,n}. The coils form one layer: axisymmetric
  rectangles (`gs_coils`) and 3D polylines (`coils_3d`) share one coil-set
  interface that returns B_coil and its per-n traces on Γ_c. Circuit incidence
  and limits stay as they are.
- The current ψ coupling is the n = 0 instance in its scalar form. By the
  F-n2 correction, ψ = R A_φ is the n = 1 vector mode of the kernel, so
  `gs_coupling` becomes an adapter onto the general layer.
- Every lane exposes the same thing to the layer: the field in a vacuum shell
  between Γ_p and Γ_c. The free-boundary Newton acts on the Γ_c trace with
  exact coil Green's functions, and the plasma-boundary position follows from
  flux matching. This is the conditioning of FreeGS and of the ψ lane, not of
  [[B²]] targets.
- Γ_p conditions are lane-specific:
  - ψ lane: implicit (weak continuity, cut-cell interface term);
  - map lane: natural traction [[p + B²/2μ0]] = 0 with B·n = 0 on a material
    surface;
  - free-B lane: absent, since vacuum is the region where p = 0 and J = 0.
- The inverse coil problem (`gs_coil_inverse`) acts on functionals of 𝒮
  (flux through a disk, B at a point), not on nodal ψ. Its Gauss–Newton,
  BVLS and GSVD parts are lane-independent.

## D7. One readback and export layer

Every consumer reads 𝒮 through one evaluator:

- point values of B, p, J and ∇B (batched, as `sample_gs_many`);
- the magnetic axis as a fixed point of the poloidal field-line map;
- label level sets where Λ provides a label;
- q, Φ_tor and volume from flux integrals of B over those sets;
- SFL angles (PEST, Hamada, Boozer);
- exports: G-EQDSK (axisymmetric states only, ψ from D1), Boozer `.bc`
  through libneo, and the L7 spectra.

In the map lane, q per label is an input (frozen fluxes). Readback must
recompute it from B in physical space, because that is the cross-lane check;
it must not copy the reference. `vacuum_surfaces` (ψ-grid from G-EQDSK) is
kept as an import-side adapter that builds a core state. It is not a second
readback.

## D8. One verification harness

There is one case registry: `run(case, lane, resolution) -> 𝒮 + report`. It
reports the following for each lane.

1. Oracle error through D7. The cases are Solovev A3, shaped Solovev,
   Cerfon–Freidberg, the cylinder, the Landreman E0 families, vacuum coil
   fields (curl-free, exact), and FreeGS FB-cut parity.
2. The D3 residual, normalised by σ.
3. Reductions:
   - the map lane at N = 0 against prescribed-q GS (the existing N1 ladders);
   - the map lane with N ≥ 1 and axisymmetric data must show no leakage;
   - the free-B lane seeded from an axisymmetric state must stay there.
4. Offset-start return tests. Seeds with a displaced axis, boundary or shape
   (for example 15 mm, or 2 % elongation) must return to the reference.
   Stationarity of the reference alone is not acceptance (GVEC lesson).
   The N1 interior-shift seed is the first instance.
5. Cross-lane checks. The 𝒮 of lane A is evaluated by D3 and D7 in lane B's
   discretisation, and it seeds lane B, which must stay stationary to
   discretisation accuracy.
6. The D5 identities.

Solvers are Newton or Newton–Krylov with exact derivatives and shared
globalisation (Armijo or trust region, quadrature tightening, det DX > 0
guard). Gradient descent with heuristic step control is not an acceptance
path.

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

Sum: about 12.8k lines today. About 6.0k would be shared, 5.7k ψ-specific
and 0.8k map-specific. The core is not new code. Most of it already exists in
`gs_*` behind ψ-specific interfaces.

Library boundaries (fortran rules) are unchanged.

- FortFEM: generic FEEC spaces and Piola maps (F-n1 twisted per-n complex),
  commuting projections, level-cut quadrature, and per-n BEM operators
  (F-n2/F-n3).
- FortNum: Q_{n−1/2}, elliptic integrals, BVLS, GSVD, Krylov.
- FortSparse: factorisations.
- fortsym: the generator. The energy and residual leaves are physics, so they
  live in KIN6D `gen/`.
- libneo: G-EQDSK and Boozer I/O.

## Migration order (K6 stays reproducible)

K6 is the accepted TC24 phase-4b state. Its path is `gs_phase1` → forward
cubic-profile Newton (`gs_nonlinear`) on P3 with a direct factor, then
`gs_readback`, `gs_coarea` and `gs_majorant`. The rules are:

- No migration step changes the numerics of `grad_shafranov` or `gs_nonlinear`.
  Only interfaces and callers move.
- A CTest K6-mini case (same law, same path, coarse mesh) is added in M0. It
  must reproduce `gs_phase1` outputs to ≤ 1e-12 relative at every step.
- The full K6 rerun (≈ 100 s on 48 threads) is repeated at M2 and M4. It must
  match the accepted Bpol, q and ψ figures to the last reported digit.

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
| Non-nested branch | – | M0–M4 and M6; its seed comes from 𝒮 of the ψ or map lane; its p-model is declared first |

L3 (vacuum response) depends on no core piece beyond M4's coil interface and
may continue. Its Γ_c and BEM_n code is written as the D6 layer from the start.

## What is and is not novel

These parts are not novel:

- the GS FEM-BEM coupling on a fixed curve (CEDRES++, FEEQS);
- the frozen-in energy over maps (Kruskal–Kulsrud, VMEC, GVEC, DESC);
- per-n toroidal Green's functions (Chance's VACUUM);
- the ideal energy principle and δW (Bernstein et al.);
- FEEC magnetic relaxation in physical space (MRX);
- the Maxwell-stress weak form;
- Newton–Krylov.

These parts we believe are new in combination, and they still need a
literature check before any claim:

1. One residual and certificate interface across 2D ψ, nested-map and
   Eulerian lanes, with the cross-lane checks of D8. Each code family above
   certifies, if at all, only within its own representation.
2. A nested-map lane on an unstructured curved FE reference mesh with an
   abstract chart and no flux coordinates in the solver. VMEC uses flux
   surfaces, GVEC a radial B-spline in a flux-like coordinate, DESC Zernike.
   Whether this helps near separatrix-adjacent boundaries is unproven.
3. The discrete identity between the nonlinear map Hessian and the per-n
   linear ideal operator, used both as preconditioner and as consistency test.
4. 3D free boundary through a fixed Γ_c FEM shell plus a per-n exact
   exterior, instead of virtual casing on the moving plasma boundary
   (NESTOR, BIEST). The motivation is conditioning (GVEC diagnosis). The
   advantage is a hypothesis until NL-3a passes offset-start return tests.

Not claimed: better performance than MRX in the non-nested regime, or than
VMEC++/GVEC/DESC in the smooth nested regime. Those claims need time-to-accuracy
measurements.

## Open mathematical questions

1. **Push-forward exactness.** B = DX b̂/det DX is exactly divergence-free
   cellwise and has continuous normal flux. With b̂ = rot χ̂ + b̂_3 and X
   polynomial per cell, it lies in parametric (Piola-mapped) RT on the
   deformed mesh. The open parts:
   - (a) Transfer to a fixed Eulerian mesh, for the free-B lane or for
     comparison, is exact only with a conservative 2-form remap (supermesh) or
     a commuting projection evaluated on the deformed geometry. Its cost and
     accuracy are unmeasured.
   - (b) Fourier truncation of B itself breaks exactness. Export must
     therefore carry the map and b̂, not truncated B modes.
2. **Gauges.**
   - Map lane: the parallel relabelling null space. ω_0 = f(χ) is a continuum
     null direction at n = 0, which N1 currently leaves to damping. A gauge
     that is unbiased and keeps P_k order is not settled (flux-aligned
     reference mesh K-M1 is the candidate).
   - Free-B lane: B with a div constraint (multiplier) against A ∈ H(curl)
     with tree–cotree gauge. Undecided; the harmonic data in Λ is common.
3. **Certificate for the map lane.** W is non-convex. δW has continuous
   spectrum, and with sheared ι it has no uniform inf-sup in H¹ at rationals.
   β_h may then tend to 0 with h. The expected outcome is "residual verified"
   for the map lane with error bounds only away from rationals. Which
   displacement-space norm makes the D3 dual norm stable is open.
4. **ψ lane against ideal constraints.** The forward ψ lane fixes p(ψ) and
   F(ψ); the map lane fixes mass and flux per label. Their Jacobians differ by
   profile-response terms. D5 therefore ties the map lane to prescribed-q GS,
   not to forward GS. Whether a current-constrained f(χ̂) map lane reproduces
   the forward ψ lane's Jacobian exactly is open.
5. **Pressure in the free-B lane.** B·∇p = 0 gives a singular Eulerian
   Jacobian on rational and ergodic sets. Possible p-models are anisotropic
   diffusion (χ_∥/χ_⊥), MRX-type relaxation and resistive closure. These are
   model choices (R12), and the core does not decide between them.
6. **Free-boundary conditioning in 3D.** That the shell + Γ_c formulation
   keeps the ψ lane's position conditioning when the plasma boundary is a
   material surface of the map is a hypothesis. NL-3a with offset starts
   decides it.

Chris&AI
