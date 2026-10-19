# Free-boundary equilibrium: native formulation and proof obligations

Status: theory specification, stage 1 (2026-10-19), repaired after the independent proof audit of 2026-10-19. No code change. Implementation: `src/mhd/ideal_free_boundary.f90`.

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
- Shell potential Φ = F τ₁ + c τ₂ + Φ̃, where τ₁, τ₂ are fixed closed 1-forms on the shell with periods (2π, 0) and (0, 1) on the toroidal and poloidal cycles (verified, repair 2026-10-19; τ₁ = dq3 has period 2π). The code uses τ₁ = dq3 and τ₂ = dθ̂/2π (`src/mhd/ideal_free_boundary.f90:707`); harmonicity is not required. F is the toroidal constant (2πF = toroidal circulation) and c = −μ0 I_p the poloidal circulation. Both are fixed.
- Period hypothesis H_per (hypothesis, per case): the periods of the exterior coil plus filament field B_e + c B_f on the cycles of Γ_c equal (2πF, c), and a is defined relative to F τ₁ + c τ₂. Otherwise a sheet current remains on Γ_c. τ₁ = dq3 is a fixed physical form, not transported; its change under the shell map is the exact form d(θ·∇q3), absorbed by Φ̃ (verified).
- Shell map: the harmonic extension of the Γ_p displacement, imposed by multipliers (rows ∇C). It is a continuum gauge of the reduced functional (L at Φ̃ = Φ̃*(X)).
- Exterior: an abstract symmetric Steklov–Poincaré operator S on the fixed closed surface Γ_c, ⟨Su, v⟩ = ∫_ext ∇w_u·∇w_v. Coil drive: traces (a, h) of the field of all sources outside Γ_c plus a unit-circulation field B_f. Realisations:
  - **R1**, per-n toroidal-mode BEM. Implemented. Verified cylindrical: L_n = −div(R∇·) + n²/R on an axisymmetric Γ_c (FortFEM `toroidal_mode_boundary_operators_2d`, header; `src/equilibrium/equilibrium_boundary.f90:2–28`; `src/mhd/ideal_free_boundary.f90:435–479`). The coil sampler takes (R, Z, k) and returns cos/sin(nφ) harmonics (`src/mhd/ideal_free_boundary.f90:87–96`).
  - **R2**, chart-free 3D Laplace BEM on a general Γ_c, symmetric coupling S = W + (½ − K′)V⁻¹(½ − K). Not implemented. FortFEM `2ed9d024` has 3D Laplace Galerkin/DtN operators only for an exactly parametrised circular torus (`laplace_torus_curved_bem_3d.f90`, `major_radius, minor_radius` arguments). A general curved Γ_c is missing upstream.

### Functional

Fixed coil currents (default):

\[
L[X,\tilde\Phi]=W_p[X]-\frac1{\mu_0}\Big(\int_V\frac{|\nabla\Phi|^2}{2}dV+\frac12\langle S u,u\rangle_{\Gamma_c}-\langle h,\tilde\Phi\rangle_{\Gamma_c}\Big)+\lambda\cdot C(X),\qquad u=\tilde\Phi|_{\Gamma_c}-a ,
\]

with W_p the N1 energy and C the shell-map constraints. L is strictly concave in Φ̃ in the continuum (vacuum co-energy) and a saddle overall. Discretely, H_ΦΦ = −(A_h + S_h)/μ0 ≺ 0 requires A_h + S_h ≻ 0, a per-state hypothesis (FB-A2, open).

- Variant: fixed fluxes (ideal wall at Γ_c, or fixed coil fluxes) with the true vacuum energy. Same equilibria, different Hessian (obligation FB-A4).

### Gauge and vertical coordinate

- Position gauge: δq3 = 0 on plasma nodes. It removes the parallel relabelling λB. Stationarity in the position gauge implies full stationarity, because ⟨∇L, λB⟩ = 0 at Φ̃ = Φ̃*(X) (verified, repair 2026-10-19).
- Slide gauge: bordered moments M_l(ξ) = ∫ Λ_l(χ) ∇q3·(∇×ξ) dV = ∫ φ_l(χ)(∇χ×∇q3)·ξ dV, with φ_l = Λ_l′ a C¹ partition in ρ and Λ_l(χ_b) = 0.
  - They remove the remaining exact null ξ_f = f(χ)∇χ×∇q3/(B·∇q3).
  - They vanish on every curl-free displacement.
  - In the cylindrical axisymmetric case this is the implemented `slide_directions`. For a 3D base state or a curve frame, the moments act on all reference Fourier modes.
- The vertical mode is not gauged. It is an ordinary negative direction of the reduced Hessian S_red = ZᵀS_mapZ at fixed currents, Z a basis of ker[∇C; M] on the map unknowns. The bordered Newton matrix is nonsingular iff [∇C; M] has full row rank and the gauged KKT block is nonsingular (Lemma B1); if A_h + S_h ≻ 0, the latter is equivalent to S_red having no zero eigenvalue.
- The rigid-Z coordinate is a diagnostic (F_Z, K_Z) and an optional coordinate within the slide-gauged space. Combined with the normal-only restriction (`tangential_gauge`) it is not a gauge (derivation, section 7).

### Algorithm

1. Solve Φ̃ and λ at the seed.
2. Run Newton on the slide-gauged KKT system: exact generated Hessian, symmetric indefinite LDLᵀ/LU, residual merit ½‖F‖²_D in a scaling D frozen per iteration, trust radius. Never energy descent.
3. Continue in parameters (current, β, coil amplitude). Monitor the inertia of the reduced Hessian, so that a fold, a branch point (e.g. the up-down-symmetric pitchfork at K_Z = 0) or a marginal crossing is detected rather than stepped over.
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

The column "Audit 2026-10-19" records the verdict of the independent proof audit of 2026-10-19 after the repairs of the same date: **PASS** (checked, no repair), **valid with repair** (repaired here), **conditional** (holds under the named hypotheses), **open**, or **not audited**. The audit checked continuum, smooth statements; no Sobolev, discrete-assembly or continuum-certificate statement was audited.

| Id | Statement (continuous unless marked discrete) | Status | Audit 2026-10-19 | Checked by | Remains |
|---|---|---|---|---|---|
| FB-A1 | First variation of L (chart- and S-abstract): Φ̃ ⇒ curl/div-free vacuum, B·n = 0 on Γ_p natural, B continuous on Γ_c (tangential part under H_per); X ⇒ J×B = ∇p and [[p + B²/2μ0]] = 0; shell map is a gauge of the reduced functional | P (C²; [first variation](../theory/derivations/free-boundary-first-variation.md) §3–6); C [KK 1958](https://doi.org/10.1063/1.1705884), [BFKK 1958](https://doi.org/10.1098/rspa.1958.0023) | valid with repair (τ₁ period 2π; H_per; non-transported dq3 absorbed by Φ̃; gauge wording) | discrete: `ideal_fb_derivatives` (N); D3 stress leaf (FortSym) | Sobolev version |
| FB-A2 | S symmetric and positive on H^{1/2}(Γ_c) (3D decaying exterior); discrete A_h + S_h ≻ 0 (so H_ΦΦ ≺ 0) | P (continuous, positivity; coercivity not shown); discrete O, used as hypothesis by FB-D2, FB-E1 | PASS (continuous); discrete open | code symmetrises S_h (`src/mhd/ideal_free_boundary.f90:475`) | report the asymmetry defect before symmetrisation, λ_min(S_h) and λ_min(A_h + S_h); R2 not built |
| FB-A3 | Reduced Hessian at equilibrium = δW_F + δW_S + δW_V at fixed currents (BFKK, wall at infinity); discrete Hessian = discrete δW_n | C [BFKK 1958](https://doi.org/10.1098/rspa.1958.0023), [Chance 1997](https://doi.org/10.1063/1.872380) | not audited | none (D5/M5 planned) | reproduce; FB form of the D5 identity |
| FB-A4 | Fixed-flux and fixed-current formulations have the same equilibria; Hessians differ by the response term; wall stabilises n = 0 vertical | H (sketch, first variation §7) | not audited | none | derivation |
| FB-B1 | Discrete residual = exact gradient of the discrete Lagrangian; Hessian = its exact derivative | P for leaves (FortSym `diff` of one density, `gen/gen_vacuum_shell.f90`, `gen/gen_ideal_map.f90`); assembly N | not audited | `ideal_fb_derivatives`, `ideal_fb_vertical_newton` (FD) | structural proof of the assembly (FortAD oracle); third-derivative leaves for FB-E2 |
| FB-B2 | Identity-map reduction of the vacuum leaf to R\|∇Φ\|²/2μ0 | P (FortSym native zero test, `gen/gen_vacuum_shell.f90:52`) | not audited | generator | analogous curve-frame reduction once that chart exists |
| FB-C1 | In the position gauge, relabelling null = {f(χ)∇χ×∇q3/(B·∇q3)} for any chart; invariance on the position-gauged set; position-gauged stationarity ⇒ full stationarity | P under H0 ([gauge](../theory/derivations/free-boundary-gauge-vertical.md) §1) | valid with repair ("only if" dropped; invariance restricted to position-gauged X; stationarity bridge added) | none | — |
| FB-C2 | Slide moments vanish on every curl-free displacement | P (§2 i–ii) | PASS | `ideal_fb_slide_gauge_curl_free` (N: 1.9e-3 → 5.4e-4) | discrete defect order |
| FB-C3 | Moment Gram matrix SPD (at X₀); full family G_∞ gives a direct sum with the null | P for smooth ξ; axis regularity H (§2 iii–iv) | PASS (smooth, continuum) | none | Sobolev version |
| FB-C4 | Discrete: L + 1 layer moments remove every discrete near-null | N (six fixture starts converge) | not audited | `gs_free_3d_nl3a_fixture_return` | smallest gauged eigenvalue under refinement |
| FB-C5 | At a gauged solution near X₀ the multiplier vanishes, so it is an equilibrium | P continuum, local (radius unquantified), via the FB-C1 bridge (§3); discrete μ = O(h) H | valid with repair (local; depends on FB-C1 bridge) | none | report μ; quantify the neighbourhood |
| FB-C6 | 3D/figure-8 slide gauge on all reference modes | P (formula); not implemented | not audited (note: curl-free motions free only modulo ξ_f in 3D) | none | implementation and test |
| FB-D1 | Bordering lemma; inertia In(K) = In(ZᵀHZ) + (r, r, 0) on ker M, rank M = r | P (§4); Lean target | PASS | none | `proof/Kin6d/Bordered.lean` |
| FB-D2 | In(K_full) = (0, dim Φ̃, 0) + In(S_red) + (r, r, 0), S_red = ZᵀS_mapZ, Z a basis of ker[∇C; M] on map unknowns; physical Morse index neg(S_red); fixed-current vertical is a negative direction of S_red | P conditional on A_h + S_h ≻ 0 (FB-A2 discrete) and full row rank of [∇C; M]; C [Haynsworth 1968](https://doi.org/10.1016/0024-3795(68)90050-5); negativity N (K_Z < 0) | first invalid bridge (index off by rank C), repaired; conditional | `ideal_fb_rigid_stiffness`, `ideal_fb_mode_stability` | check both hypotheses per state; neg(S_red) = 1 per state (H) |
| FB-D3 | Newton direction is a descent direction for the residual merit in a fixed scaling | P (§6) | PASS | none | check that `ideal_fb_solve` compares merits in one scaling |
| FB-D4 | Normal-only + rigid-Z restriction is not equivalent to ∇L = 0 | P (§7), counterexample N | not audited | rp/tangential_gauge run | remove from production path |
| FB-E1 | Local quadratic convergence of pure gauged Newton if the gauged KKT block is nonsingular (⇔ S_red nonsingular when A_h + S_h ≻ 0, [∇C; M] full rank) | P ([iterations](../theory/derivations/free-boundary-iterations.md) §1–2); C [Deuflhard 2011](https://doi.org/10.1007/978-3-642-23899-4); nonsingularity N; full-step acceptance for the globalised solver O | theorem PASS; applicability valid with repair, conditional | `ideal_fb_vertical_newton` (quadratic return, manufactured) | per-state nonsingularity; full-step-acceptance lemma for `ideal_fb_solve` |
| FB-E2 | Discrete Kantorovich/radii-polynomial certificate: existence, uniqueness in B(x₀, ρ₋), rate; norm ‖D·‖₂ | C (theorem; Deuflhard theorem number, content not verified); application O | statement matches classical form; content not read | none | outward Y, third-derivative/Lipschitz bound, β_h enclosure |
| FB-E3 | Continuous certificate | O (two-norm gap; no uniform inf-sup at rationals, core OQ3) | open | none | research problem |
| FB-F1 | Descent converges to a saddle, locally, only from a local C¹ manifold of codimension ≥ neg(S) (non-hyperbolic: ξ_f nulls) | P (local, centre-stable manifold); global measure-zero claim O ([Lee et al. 2019](https://doi.org/10.1007/s10107-019-01374-3) needs a globally Lipschitz gradient on R^d) | valid with repair (local only) | none | global statement on U |
| FB-F2 | Descent rate (κ−1)/(κ+1); κ ≥ K_R/\|K_Z\| if the rigid R and Z shifts have equal P-norm and the vertical mode is stable; heavy ball (√κ−1)/(√κ+1) | P; C [Polyak 1964](https://doi.org/10.1016/0041-5553(64)90137-5); fixture K_R/\|K_Z\| ≈ 556 ⇒ ≈ 6396 plain / 271 heavy-ball iterations for d = 10 (P = I) | valid with repair (equal-P-norm hypothesis) | none | — |
| FB-F3 | Picard/frozen-B² updates: linear rate ρ(DΦ), oscillation from non-symmetric DΦ | C [Ortega–Rheinboldt 2000](https://doi.org/10.1137/1.9780898719468) (Ostrowski clause, content not verified); application to GVEC H | cited/hypothesis as labelled | none | FD spectrum of DΦ at a converged state |
| FB-G1 | Inside virtual casing B_in χ_Ω = BS[μ0Jχ_Ω] + BS[(B_in×n)δ]; B_J continuous on R³ | P ([benchmark](../theory/derivations/free-boundary-landreman-benchmark.md) §1); C [Shafranov–Zakharov 1972](https://doi.org/10.1088/0029-5515/12/5/009), [Lazerson 2012](https://doi.org/10.1088/0741-3335/54/12/122002) | PASS (regularity statement corrected: C¹ Newton potential, B_J continuous) | none | — |
| FB-G2 | Exterior B_J + B̃_x gives an exact surface-current-free free-boundary equilibrium with p̃ = p − p(S0); necessity and uniqueness | P conditional on G4 (§2) | PASS, conditional on G4 | none | — |
| FB-G3 | Circulations: c = −μ0 I_tor from B_J; 2πF = ∮_γ B_in·dl on the reference toroidal cycle; from B̃_x only iff lk(γ, axis) = 0 (else + k μ0 I_tor) | P (§3); linking number H per case | valid with repair (linking condition) | none | compute lk of the reference toroidal cycle |
| FB-G4 | Harmonic continuation of B_x to the open D_c, C¹ up to Γ_c | P only in a thin unquantified shell (Cauchy–Kovalevskaya; Holmgren, [John 1982](https://doi.org/10.1007/978-1-4684-9333-7) ch. 3, content not verified); O for a given Γ_c | valid as stated (restated to match the converse) | none | validated continuation or numerical evidence (ε = 0 first) |
| FB-G5 | Benchmark data satisfy every KIN6D constraint; target is stationary; nondegenerate | P stationarity conditional on G4; nondegeneracy H (at risk for ι = 2, closed lines) | conditional on G4 | none | use sheared family/Issan for return tests |
| FB-G6 | Approximate drive: computable flux and traction defect on S0 | P (residual statement) | verified (traction expansion) | none | error bound needs FB-E2 β_h |
| FB-H1 | Axisymmetric well-posedness and solver literature | C (see §3) | metadata verified; content not audited | — | read full texts; no 3D free-boundary existence theorem adopted |

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
| shell multivalued part | F dq3 + c dθ̂/2π, chart coordinate (`src/mhd/ideal_free_boundary.f90:707`) | none (harmonicity not needed) |
| exterior S | R1 per-n cylindrical, axisymmetric Γ_c (verified) | R2 chart-free 3D BEM on a general Γ_c (FortFEM upstream) |
| coil drive | sampler (R, Z, k) with cos/sin(nφ) (`src/mhd/ideal_free_boundary.f90:87–96`, `:481–580`) | surface traces (a, h) on Γ_c points; polyline coils already Cartesian |
| unit-circulation field | n = 0 filament in the cylindrical construction | closed filament along the axis curve, Biot–Savart |
| slide gauge | n = 0 rows, 2D reference flux (`src/mhd/ideal_free_boundary.f90:1345–1393`) | full 3D moments with ∇q3 (FB-C6) |
| shell field output | cylindrical components (`src/mhd/ideal_free_boundary.f90:2451–2500`) | Cartesian (DESIGN rule) |
| Γ_c | scaled 2D boundary | closed surface inside the frame tube and outside the plasma |

## 5. Implementation gap list (not implemented in this stage)

1. Curve-frame chart H = Γ(q3) + q1 N(q3) + q2 B(q3), with jets to third order, in `equilibrium_chart`. Includes an injectivity-radius check and a FortSym identity-reduction test (gate `g_frame_figure8_equilibrium`, [ROADMAP](../ROADMAP.md) R12).
2. Chart-free exterior operator R2: 3D Laplace Calderón operators on general curved Γ_c panels in FortFEM (upstream), symmetric coupling, plain surface traces of coil fields, and a closed-filament B_f.
3. 3D slide gauge on all reference modes, with the moment formula ∫φ_l(∇χ×∇q3)·ξ (FB-C6). Test: curl-free defect and Gram positivity on a 3D frozen state.
4. Remove `tangential_gauge` from the production path. Keep F_Z and K_Z as reported diagnostics (FB-D4).
5. Merit acceptance in one frozen scaling per iteration (FB-D3); verify in `ideal_fb_solve`. State and prove the full-step-acceptance lemma for the globalised solver (FB-E1).
6. Report μ (gauge multiplier), the asymmetry of S_h before symmetrisation, λ_min(S_h), λ_min(A_h + S_h), the row rank of [∇C; M], and the inertia of S_red (FB-A2, FB-C5, FB-D2).
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
