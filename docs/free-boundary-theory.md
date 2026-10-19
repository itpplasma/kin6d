# Free-boundary equilibrium: native formulation and proof obligations

Status: theory specification, stage 1 (2026-10-19). No code change.

This document owns the mathematical formulation of KIN6D's free-boundary map-lane equilibrium and the numbered proof-obligation map. Other ownership is unchanged:

- implementation and measurements: [ideal-free-boundary.md](ideal-free-boundary.md) (N4);
- boundary layer: [equilibrium-core.md](equilibrium-core.md) D6;
- plasma energy: [ideal-map-3d.md](ideal-map-3d.md) (N1);
- certificate semantics: [certification.md](certification.md);
- figure-8 gate: [geometry.md](geometry.md#figure-8-and-geometry-gates).

Derivations live in `theory/derivations/free-boundary-*.md`:

- [first variation](../theory/derivations/free-boundary-first-variation.md);
- [gauge, bordered system and vertical saddle](../theory/derivations/free-boundary-gauge-vertical.md);
- [Newton, Kantorovich, descent and Picard](../theory/derivations/free-boundary-iterations.md);
- [exact benchmark from Landreman's equilibria](../theory/derivations/free-boundary-landreman-benchmark.md).

The formulation is native and chart-abstract. It is not VMEC/NESTOR (virtual casing on the moving boundary, descent), DESC (collocation) or GVEC (frozen-B² boundary push).

## 1. Formulation

### Geometry, chart-abstract

- Reference domain Ω̂ = Ω̂_p ∪ Ω̂_s: a curved P_k reference mesh of the region inside Γ̂_c, with the plasma boundary Γ̂_p as an internal chain, times a reference angle φ̂ (Fourier, n ∈ nfp ℤ).
- Map X = H∘q. The chart H is an `ideal_chart_t` and enters only through its jet (H, dH, d²H, d³H). Two charts are needed:
  - cylindrical (implemented);
  - curve frame H = Γ(q3) + q1 N(q3) + q2 B(q3) (planned; required for figure-8, [Plunk et al. 2025](https://doi.org/10.1088/1361-6587/adb64b), [Hindenlang, Plunk and Maj, arXiv:2410.17595](https://arxiv.org/abs/2410.17595)).
- Every statement below uses only X, DX, det DX, dH and the chart coordinate q3. No statement needs a global cylindrical angle.
- Chart hypothesis (per case): det dH > 0 on the image (within the frame's injectivity radius for a curve frame), and B·∇q3 > 0 on the plasma (field lines advance in q3).

### State and constraints

- Plasma frozen data:
  - closed reference 2-form b̂ with b̂·n̂ = 0 on Γ̂_p, label χ̂, p(χ̂), γ = 0;
  - B = DX b̂/det DX, p∘X = p̂.
  
  Topology, toroidal flux and the flux and pressure of every label are fixed (Kruskal–Kulsrud).
- Shell potential Φ = F τ₁ + c τ₂ + Φ̃, where τ₁, τ₂ are fixed closed 1-forms on the shell with unit periods. The code uses τ₁ = dq3 and τ₂ = dθ̂/2π (`src/mhd/ideal_free_boundary.f90:707`); harmonicity is not required. F is the toroidal constant and c = −μ0 I_p the poloidal circulation. Both are fixed.
- Shell map: the harmonic extension of the Γ_p displacement, imposed by multipliers. It is a continuum gauge.
- Exterior: an abstract symmetric Steklov–Poincaré operator S on the fixed closed surface Γ_c, ⟨Su, v⟩ = ∫_ext ∇w_u·∇w_v. Coil drive: traces (a, h) of the field of all sources outside Γ_c plus a unit-circulation field B_f. Realisations:
  - **R1**, per-n toroidal-mode BEM. Implemented. Verified cylindrical: L_n = −div(R∇·) + n²/R on an axisymmetric Γ_c (FortFEM `toroidal_mode_boundary_operators_2d`, header; `src/equilibrium/equilibrium_boundary.f90:2–28`; `src/mhd/ideal_free_boundary.f90:435–479`). The coil sampler takes (R, Z, k) and returns cos/sin(nφ) harmonics (`ideal_free_boundary.f90:87–96`).
  - **R2**, chart-free 3D Laplace BEM on a general Γ_c, symmetric coupling S = W + (½ − K′)V⁻¹(½ − K). Not implemented. FortFEM `2ed9d024` has 3D Laplace Galerkin/DtN operators only for an exactly parametrised circular torus (`laplace_torus_curved_bem_3d.f90`, `major_radius, minor_radius` arguments). A general curved Γ_c is missing upstream.

### Functional

Fixed coil currents (default):

\[
L[X,\tilde\Phi]=W_p[X]-\frac1{\mu_0}\Big(\int_V\frac{|\nabla\Phi|^2}{2}dV+\frac12\langle S u,u\rangle_{\Gamma_c}-\langle h,\tilde\Phi\rangle_{\Gamma_c}\Big)+\lambda\cdot C(X),\qquad u=\tilde\Phi|_{\Gamma_c}-a ,
\]

with W_p the N1 energy and C the shell-map constraints. L is concave in Φ̃ (vacuum co-energy) and a saddle overall.

- Variant: fixed fluxes (ideal wall at Γ_c, or fixed coil fluxes) with the true vacuum energy. Same equilibria, different Hessian (obligation FB-A4).

### Gauge and vertical coordinate

- Position gauge: δq3 = 0 on plasma nodes. It removes the parallel relabelling λB.
- Slide gauge: bordered moments M_l(ξ) = ∫ Λ_l(χ) ∇q3·(∇×ξ) dV = ∫ φ_l(χ)(∇χ×∇q3)·ξ dV, with φ_l = Λ_l′ a C¹ partition in ρ and Λ_l(χ_b) = 0.
  - They remove the remaining exact null ξ_f = f(χ)∇χ×∇q3/(B·∇q3).
  - They vanish on every curl-free displacement.
  - In the cylindrical axisymmetric case this is the implemented `slide_directions`. For a 3D base state or a curve frame, the moments act on all reference Fourier modes.
- The vertical mode is not gauged. It is an ordinary negative direction of the reduced Hessian at fixed currents. The bordered Newton matrix is nonsingular iff the reduced Hessian on the gauge complement has no zero eigenvalue (Lemma B1).
- The rigid-Z coordinate is a diagnostic (F_Z, K_Z) and an optional coordinate within the slide-gauged space. Combined with the normal-only restriction (`tangential_gauge`) it is not a gauge (derivation, section 7).

### Algorithm

1. Solve Φ̃ and λ at the seed.
2. Run Newton on the slide-gauged KKT system: exact generated Hessian, symmetric indefinite LDLᵀ/LU, residual merit ½‖F‖²_D in a scaling D frozen per iteration, trust radius. Never energy descent.
3. Continue in parameters (current, β, coil amplitude). Monitor the inertia of the reduced Hessian, so that a fold or a marginal crossing is detected rather than stepped over.
4. Accept with an absolute residual target in σV units. Report F_Z/K_Z.
5. Per state, attempt a discrete Kantorovich/radii-polynomial certificate (FB-E2). The claim is "candidate" until it passes.

## 2. Proof obligations

Status codes:

- **L**: proved in Lean;
- **P**: proved on paper here;
- **C**: literature, with a verified primary reference (not reproduced);
- **N**: numerical evidence only;
- **O**: open;
- **H**: hypothesis, per case.

| Id | Statement (continuous unless marked discrete) | Status | Checked by | Remains |
|---|---|---|---|---|
| FB-A1 | First variation of L (chart- and S-abstract): Φ̃ ⇒ curl/div-free vacuum, B·n = 0 on Γ_p natural, B continuous on Γ_c; X ⇒ J×B = ∇p and [[p + B²/2μ0]] = 0; shell map is a gauge | P (C²; [first variation](../theory/derivations/free-boundary-first-variation.md) §3–6); C [KK 1958](https://doi.org/10.1063/1.1705884), [BFKK 1958](https://doi.org/10.1098/rspa.1958.0023) | discrete: `ideal_fb_derivatives` (N); D3 stress leaf (FortSym) | proof-audit; Sobolev version |
| FB-A2 | S symmetric and positive definite on H^{1/2}(Γ_c) (3D decaying exterior); discrete S_h likewise | P (continuous); discrete O | code symmetrises S_h (`ideal_free_boundary.f90:475`) | report the asymmetry defect before symmetrisation and λ_min(S_h) > 0; R2 not built |
| FB-A3 | Reduced Hessian at equilibrium = δW_F + δW_S + δW_V at fixed currents; discrete Hessian = discrete δW_n | C [BFKK 1958](https://doi.org/10.1098/rspa.1958.0023), [Chance 1997](https://doi.org/10.1063/1.872380) | none (D5/M5 planned) | reproduce; FB form of the D5 identity |
| FB-A4 | Fixed-flux and fixed-current formulations have the same equilibria; Hessians differ by the response term; wall stabilises n = 0 vertical | H (sketch, first variation §7) | none | derivation |
| FB-B1 | Discrete residual = exact gradient of the discrete Lagrangian; Hessian = its exact derivative | P for leaves (FortSym `diff` of one density, `gen/gen_vacuum_shell.f90`, `gen/gen_ideal_map.f90`); assembly N | `ideal_fb_derivatives`, `ideal_fb_vertical_newton` (FD) | structural proof of the assembly (FortAD oracle); third-derivative leaves for FB-E2 |
| FB-B2 | Identity-map reduction of the vacuum leaf to R\|∇Φ\|²/2μ0 | P (FortSym native zero test, `gen/gen_vacuum_shell.f90:52`) | generator | analogous curve-frame reduction once that chart exists |
| FB-C1 | In the position gauge, relabelling null = {f(χ)∇χ×∇q3/(B·∇q3)} for any chart | P under H0 ([gauge](../theory/derivations/free-boundary-gauge-vertical.md) §1) | none | — |
| FB-C2 | Slide moments vanish on every curl-free displacement | P (§2 i–ii) | `ideal_fb_slide_gauge_curl_free` (N: 1.9e-3 → 5.4e-4) | discrete defect order |
| FB-C3 | Moment Gram matrix SPD; full family G_∞ gives a direct sum with the null | P for smooth ξ; axis regularity H (§2 iii–iv) | none | Sobolev version |
| FB-C4 | Discrete: L + 1 layer moments remove every discrete near-null | N (six fixture starts converge) | `gs_free_3d_nl3a_fixture_return` | smallest gauged eigenvalue under refinement |
| FB-C5 | At a gauged solution the multiplier vanishes, so it is an equilibrium | P continuum (§3); discrete μ = O(h) H | none | report μ |
| FB-C6 | 3D/figure-8 slide gauge on all reference modes | P (formula); not implemented | none | implementation and test |
| FB-D1 | Bordering lemma; inertia In(K) = In(ZᵀHZ) + (r, r, 0) | P (§4); Lean target | none | `proof/Kin6d/Bordered.lean` |
| FB-D2 | In(K_full) = In(H_ΦΦ) + In(S), H_ΦΦ ≺ 0; fixed-current vertical is a negative direction of S \| ker M | P (inertia); C [Haynsworth 1968](https://doi.org/10.1016/0024-3795(68)90050-5); negativity N (K_Z < 0) | `ideal_fb_rigid_stiffness`, `ideal_fb_mode_stability` | Morse index = 1 per state (H) |
| FB-D3 | Newton direction is a descent direction for the residual merit in a fixed scaling | P (§6) | none | check that `ideal_fb_solve` compares merits in one scaling |
| FB-D4 | Normal-only + rigid-Z restriction is not equivalent to ∇L = 0 | P (§7), counterexample N | rp/tangential_gauge run | remove from production path |
| FB-E1 | Local quadratic convergence of gauged Newton if S\|ker M is nonsingular | P ([iterations](../theory/derivations/free-boundary-iterations.md) §1–2); C [Deuflhard 2011](https://doi.org/10.1007/978-3-642-23899-4); nonsingularity N | `ideal_fb_vertical_newton` (quadratic return, manufactured) | per-state nonsingularity |
| FB-E2 | Discrete Kantorovich/radii-polynomial certificate: existence, uniqueness in B(x₀, ρ₋), rate; norm ‖D·‖₂ | C (theorem); application O | none | outward Y, third-derivative/Lipschitz bound, β_h enclosure |
| FB-E3 | Continuous certificate | O (two-norm gap; no uniform inf-sup at rationals, core OQ3) | none | research problem |
| FB-F1 | Descent converges to a saddle only from a positive-codimension set | P (linear); C [Lee et al. 2019](https://doi.org/10.1007/s10107-019-01374-3) | none | — |
| FB-F2 | Descent rate (κ−1)/(κ+1), κ ≥ K_R/\|K_Z\|; heavy ball (√κ−1)/(√κ+1) | P; C [Polyak 1964](https://doi.org/10.1016/0041-5553(64)90137-5) | none | — |
| FB-F3 | Picard/frozen-B² updates: linear rate ρ(DΦ), oscillation from non-symmetric DΦ | C [Ortega–Rheinboldt 2000](https://doi.org/10.1137/1.9780898719468); application to GVEC H | none | FD spectrum of DΦ at a converged state |
| FB-G1 | Inside virtual casing B_in χ_Ω = BS[μ0Jχ_Ω] + BS[(B_in×n)δ] | P ([benchmark](../theory/derivations/free-boundary-landreman-benchmark.md) §1); C [Shafranov–Zakharov 1972](https://doi.org/10.1088/0029-5515/12/5/009), [Lazerson 2012](https://doi.org/10.1088/0741-3335/54/12/122002) | none | — |
| FB-G2 | Exterior B_J + B̃_x gives an exact surface-current-free free-boundary equilibrium with p̃ = p − p(S0); necessity and uniqueness | P conditional on G4 (§2) | none | — |
| FB-G3 | Circulations: c = −μ0 I_tor from B_J; toroidal constant from B̃_x only | P (§3) | none | — |
| FB-G4 | Harmonic continuation of B_x reaches the chosen Γ_c | P only in a thin unquantified shell (Cauchy–Kovalevskaya, [John 1982](https://doi.org/10.1007/978-1-4684-9333-7)); O for a given Γ_c | none | validated continuation or numerical evidence (ε = 0 first) |
| FB-G5 | Benchmark data satisfy every KIN6D constraint; target is stationary; nondegenerate | P (stationarity); nondegeneracy H (at risk for ι = 2, closed lines) | none | use sheared family/Issan for return tests |
| FB-G6 | Approximate drive: computable flux and traction defect on S0 | P (residual statement) | none | error bound needs FB-E2 β_h |
| FB-H1 | Axisymmetric well-posedness and solver literature | C (see §3) | — | read full texts; no 3D free-boundary existence theorem adopted |

No obligation is promoted on numerical evidence. Lean coverage of this map is currently empty. FB-D1 is the first realistic Lean target.

## 3. Axisymmetric literature (FB-H1)

Bibliographic metadata verified by Crossref (2026-10-19). The content claims are limited to what was read.

- [Temam 1975, Arch. Rational Mech. Anal. 60, 51–73](https://doi.org/10.1007/BF00281469).
  - Read: the abstract (Springer page).
  - Content: reduces the confined-plasma free-boundary problem to a non-classical nonlinear eigenvalue problem, and proves existence of solutions as critical points of functionals.
  - Gives: existence for the 2D (ψ-lane) problem. It does not give uniqueness or 3D results.
- [Temam 1977, Comm. PDE 2, 563–585](https://doi.org/10.1080/03605307708820039).
  - Read: title and metadata only. Content not checked.
- [Heumann, Blum, Boulbe, Faugeras et al. 2015, J. Plasma Phys. 81](https://doi.org/10.1017/S0022377814001251).
  - Read: the abstract.
  - Content: a survey of CEDRES++ free-boundary methods with P1 FEM and FEM–BEM coupling, with coil currents or circuit voltages as data.
  - Gives: the 2D solver precedent for the ψ lane, and KIN6D's ψ-lane comparator.
- [Faugeras and Heumann 2017, J. Comput. Phys. 343, 201–216](https://doi.org/10.1016/j.jcp.2017.04.047).
  - Read: title and metadata.
  - Content (as titled): FEM–BEM coupling for axisymmetric free-boundary equilibrium in unbounded domains. This is the 2D precedent for R1.
- Not adopted: 3D existence near axisymmetry with current sheets (Bruno and Laurence 1996, Commun. Pure Appl. Math. 49, 717–764, as cited by Landreman; not verified here).
- Smooth nested 3D existence is settled by exact families ([Landreman](https://arxiv.org/abs/2609.26742); [Issan et al.](https://arxiv.org/abs/2610.07304); [Gómez-Serrano, Liehr and Taylor](https://arxiv.org/abs/2609.24739)). None addresses the free-boundary problem with given coils. The practical open question is the convergence of the free-boundary solve, which FB-D to FB-F address.

## 4. Figure-8 and non-cylindrical geometry

| Piece | Chart dependence today | Required |
|---|---|---|
| plasma energy, gradients, Hessian | chart jets only (`ideal_map`, `ideal_fb_init` takes `class(ideal_chart_t)`) | curve-frame chart in `equilibrium_chart` |
| vacuum shell leaf | DX only (`gen_vacuum_shell`) | none |
| shell multivalued part | F dq3 + c dθ̂/2π, chart coordinate (`:707`) | none (harmonicity not needed) |
| exterior S | R1 per-n cylindrical, axisymmetric Γ_c (verified) | R2 chart-free 3D BEM on a general Γ_c (FortFEM upstream) |
| coil drive | sampler (R, Z, k) with cos/sin(nφ) (`:87–96`, `:481–580`) | surface traces (a, h) on Γ_c points; polyline coils already Cartesian |
| unit-circulation field | n = 0 filament in the cylindrical construction | closed filament along the axis curve, Biot–Savart |
| slide gauge | n = 0 rows, 2D reference flux (`:1345–1393`) | full 3D moments with ∇q3 (FB-C6) |
| shell field output | cylindrical components (`:2451–2500`) | Cartesian (DESIGN rule) |
| Γ_c | scaled 2D boundary | closed surface inside the frame tube and outside the plasma |

## 5. Implementation gap list (not implemented in this stage)

1. Curve-frame chart H = Γ(q3) + q1 N(q3) + q2 B(q3), with jets to third order, in `equilibrium_chart`. Includes an injectivity-radius check and a FortSym identity-reduction test (gate `g_frame_figure8_equilibrium`, [ROADMAP](../ROADMAP.md) R12).
2. Chart-free exterior operator R2: 3D Laplace Calderón operators on general curved Γ_c panels in FortFEM (upstream), symmetric coupling, plain surface traces of coil fields, and a closed-filament B_f.
3. 3D slide gauge on all reference modes, with the moment formula ∫φ_l(∇χ×∇q3)·ξ (FB-C6). Test: curl-free defect and Gram positivity on a 3D frozen state.
4. Remove `tangential_gauge` from the production path. Keep F_Z and K_Z as reported diagnostics (FB-D4).
5. Merit acceptance in one frozen scaling per iteration (FB-D3); verify in `ideal_fb_solve`.
6. Report μ (gauge multiplier), the asymmetry of S_h before symmetrisation, λ_min(S_h), and the inertia of the reduced Hessian (FB-A2, FB-C5, FB-D2).
7. Third-derivative (or Hessian-Lipschitz) generated leaves, interval evaluation, and outward residual for a discrete Kantorovich certificate. Register it as a D4 `equilibrium_certifier_t` (FB-E2).
8. 3D frozen reference input for N4: N1/N3 frozen data in place of a pushed-forward 2D gs_free state, needed for Landreman ε > 0 and every stellarator.
9. Landreman benchmark harness:
   - Biot–Savart B_J of the closed-form J by high-order quadrature;
   - continuation fit of B̃_x and its S0 defect (FB-G4, FB-G6);
   - ε = 0 Solov'ev case first, cross-checked with `gs_free`;
   - case registered in D8.
10. Fixed-flux variant (FB-A4) only if continuation through vertical marginality needs it.
11. Lean: `Bordered.lean` (FB-D1).

Chris&AI
