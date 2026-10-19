# Free boundary: Newton, Newton–Kantorovich, and why descent and Picard fail at a saddle

- Scope: the discrete gauged free-boundary system F(x) = 0 of [free-boundary-gauge-vertical.md](free-boundary-gauge-vertical.md) (bordered slide gauge, vacuum eliminated or kept in the KKT matrix), and the fixed-point iterations used by other codes.
- Obligations FB-E1 to FB-E3 and FB-F1 to FB-F3 of [free-boundary-theory.md](../../docs/free-boundary-theory.md). Status labels: **proved**, **cited**, **hypothesis**, **open**. Repairs from the independent proof audit of 2026-10-19 are marked **verified** (re-derived in the audit and the repair) or **hypothesis**; "content not verified" marks a citation whose cited statement was not read.
- Notation: x* is a solution, J = F′, ‖·‖ a fixed norm on Rᴺ. The recommended norm is ‖x‖_D = ‖Dx‖₂ with D = diag(|H_ii|)^{1/2} frozen at the start state (energy units, so mixed map, potential and multiplier blocks are comparable).

## 1. Smoothness of the discrete system (proved)

- Every term of L_h is a finite sum, over Gauss points, of rational functions of the map coefficients with denominators det DX > 0. The chart jet H(q) enters through C^∞ functions (cos, sin for the cylindrical chart; the curve and frame for a curve frame).
- The exterior terms ½⟨S_h u, u⟩ − ⟨h, Φ̃⟩ are quadratic and linear with fixed matrices.
- Hence L_h ∈ C^∞ on the open set U = {all quadrature Jacobians > 0, chart regular}. F and J are Lipschitz on every compact convex subset of U.

## 2. Local quadratic convergence of the gauged Newton method (proved; nonsingularity is a hypothesis)

**Theorem.** Let F ∈ C²(U), F(x*) = 0, J(x*) nonsingular, and ‖J(x*)⁻¹(J(x) − J(y))‖ ≤ ω‖x − y‖ on B(x*, ρ) ⊂ U. If ‖x₀ − x*‖ < min(ρ, 2/(3ω)), the Newton iterates stay in the ball and satisfy ‖x_{k+1} − x*‖ ≤ ω‖x_k − x*‖²/(2(1 − ω‖x_k − x*‖)).

*Proof.* Let e_k = x_k − x*. By the Banach lemma, ‖J(x_k)⁻¹J(x*)‖ ≤ 1/(1 − ω‖e_k‖). Then

\[
e_{k+1}=J(x_k)^{-1}\!\int_0^1\!\big(J(x_k)-J(x^*+te_k)\big)e_k\,dt ,
\]

and the affine-covariant Lipschitz bound ∫₀¹ ω(1 − t)‖e_k‖²dt = ω‖e_k‖²/2 gives the estimate. With ω‖e₀‖ < 2/3, the contraction factor ω‖e₀‖/(2(1 − ω‖e₀‖)) < 1, and induction keeps every iterate in the ball. ∎

- Standard; see [Deuflhard 2011](https://doi.org/10.1007/978-3-642-23899-4) and [Ortega and Rheinboldt 2000](https://doi.org/10.1137/1.9780898719468) (**cited**; the proof above is reproduced and does not rely on them).
- Applicability:
  - F ∈ C² holds by section 1.
  - J(x*) = K_full is nonsingular iff the constraint rows R = [∇C; M] (shell constraints and slide moments) have full row rank and the gauged block diag(Z, I)ᵀ H diag(Z, I) is nonsingular (Lemma B1 on the whole KKT matrix). If in addition H_ΦΦ is nonsingular (sufficient: A_h + S_h ≻ 0, FB-A2 discrete, **open**), this holds iff the reduced Hessian S_red = ZᵀS_map Z is nonsingular (gauge file, section 5; **verified**, repair 2026-10-19). Full row rank of R and A_h + S_h ≻ 0 are per-state **hypotheses**. Nonsingularity itself, for a given state, is a **hypothesis** with numerical evidence: "Both signs of a 4 mm shift return quadratically" ([ideal-free-boundary.md](../../docs/ideal-free-boundary.md), line 338), and the 2–5 iteration returns of the slide-gauged fixture starts (lines 285–291).
- Globalised solver (**open**, repair 2026-10-19). The theorem is for pure Newton. `ideal_fb_solve` adds a trust radius and Armijo backtracking on the residual merit. Transferring the local quadratic rate needs the standard full-step-acceptance lemma: near x*, φ(x + d) = O(φ(x)²), so α = 1 passes Armijo with σ < ½, and the trust radius is eventually inactive. This lemma is not reproduced here, and its hypotheses (one frozen scaling per iteration, FB-D3) are not yet checked in the code.
- Without the slide gauge, J(x*) is singular, or near-singular on the discrete near-null ξ_f. The theorem then does not apply, and the drift observed along the null (lines 262–276) is expected.

## 3. A posteriori Newton–Kantorovich certificate per equilibrium (theorem cited; implementation open)

**Theorem (Kantorovich, affine-covariant form, cited).** Let x₀ ∈ U, J(x₀) nonsingular,

\[
\|J(x_0)^{-1}F(x_0)\|\le\eta,\qquad \|J(x_0)^{-1}(J(x)-J(y))\|\le\omega\|x-y\|\ \ (x,y\in B(x_0,\rho)),
\]

with h = ωη ≤ ½ and ρ ≥ ρ₋ = (1 − √(1 − 2h))/ω. Then:
- F has a zero x* in the closed ball B̄(x₀, ρ₋);
- the zero is unique in B(x₀, ρ₊) ∩ B(x₀, ρ), with ρ₊ = (1 + √(1 − 2h))/ω;
- Newton from x₀ converges to x*, quadratically when h < ½.

[Deuflhard 2011](https://doi.org/10.1007/978-3-642-23899-4), chapter 2 (content not verified: the theorem number and exact statement in the 2011 edition have not been read). The proof is not yet reproduced here: FB-E2 remains "cited".

- Equivalent radii-polynomial form ([certification.md](../../docs/certification.md), "Nonlinear steady states"): with A ≈ J(x₀)⁻¹, Y ≥ ‖AF(x₀)‖ and Z(r) ≥ sup_{B(x₀,r)}‖I − AJ(x)‖, any r with Y + Z(r)r < r proves existence and uniqueness in B̄(x₀, r).
- What a KIN6D certificate needs (all **open**):
  1. Y: an outward-rounded residual F(x₀) with interval quadrature of the generated leaves, and an enclosure of AF(x₀). The solve with the MUMPS factor is untrusted; an iterative-refinement residual bound is acceptable.
  2. ω or Z(r): a bound on the third derivative of L_h over the ball. FortSym currently emits value, gradient and Hessian leaves (`gen/gen_vacuum_shell.f90`, `gen/gen_ideal_map.f90`); a third-derivative-contraction or HVP-of-Hessian leaf with interval evaluation, or a Lipschitz bound per cell, is missing.
  3. ‖J(x₀)⁻¹‖: an enclosure of the smallest singular value on the gauged space, i.e. the β_h of the D4 interface ([equilibrium-core.md](../../docs/equilibrium-core.md), D4). Since σ_min ≤ |λ_v| (λ_v the vertical eigenvalue), β_h is at most of the |λ_v| scale, estimated by |K_Z| (numerical identification); the certificate radius is then limited by roughly 1/|K_Z| or worse (repair 2026-10-19, see gauge file section 5, conditioning).
- Claim level: this is "certified-discrete" for the discrete gauged system at fixed resolution. It says nothing about the continuum solution.

## 4. Continuous certificate (open)

- L is C² in a map space only with control of DX pointwise (W^{1,∞}). δW is a quadratic form on H¹ displacements (two-norm discrepancy, as in nonlinear elasticity).
- With sheared ι, δW has continuous spectrum accumulating at 0 near rational surfaces, so there is no uniform inf-sup in H¹ ([equilibrium-core.md](../../docs/equilibrium-core.md), open question 3).
- A continuous Newton–Kantorovich certificate therefore needs: a pair of norms where both the residual bound and the Lipschitz bound hold; a complement estimate for unresolved Fourier modes and elements; and a stability constant that does not collapse at rationals. None is available. FB-E3 is **open**.

## 5. Descent cannot converge at a saddle (proved for the linearisation; nonlinear statement cited)

- Preconditioned descent on the reduced functional (vacuum solved exactly every step):

\[
x_{k+1}=G(x_k)=x_k-\tau P^{-1}\nabla L_{\rm red}(x_k),\qquad P\succ0 .
\]

- The linearisation at x* is DG = I − τP⁻¹S, with S the Hessian of L_red on the map unknowns (shell map eliminated by its constraints, Φ̃ eliminated; no slide gauge). Its negative index is at least neg(S_red) of the gauge file, section 5. P⁻¹S is similar to P^{-1/2}SP^{-1/2}, so its eigenvalues λ_i are real, and by Sylvester's law their signs are the inertia of S.
- At a saddle (neg(S) ≥ 1, e.g. the fixed-current vertical mode), DG has the eigenvalue 1 − τλ_i > 1 for every τ > 0. x* is a linearly repelling fixed point along that eigenvector.
- Local statement (proved from the cited centre-stable manifold theorem; repair 2026-10-19). x* is not hyperbolic: the relabelling nulls ξ_f give λ = 0, i.e. the eigenvalue 1 of DG. If τ is such that DG(x*) is invertible (τλ_i ≠ 1), the centre-stable manifold theorem for C¹ local diffeomorphisms gives: orbits that remain in a small neighbourhood of x* lie on a local C¹ manifold of codimension #{i : |1 − τλ_i| > 1} ≥ neg(S) > 0. The codimension can exceed neg(S) when τλ_i > 2.
- No global statement is established here. [Lee, Panageas, Piliouras, Simchowitz, Jordan and Recht 2019](https://doi.org/10.1007/s10107-019-01374-3) (**cited**) give "the set of starts converging to a strict saddle has measure zero" for a C¹ map G with det DG ≠ 0 everywhere and a globally Lipschitz gradient on the whole of R^d. L_red lives on the open set U (it blows up as det DX → 0), and G(U) ⊂ U is not shown. So the global null-set claim is **open** for this problem.
- In floating point, rounding perturbs iterates off the local centre-stable manifold, so descent is not expected to settle at a vertically unstable fixed-current equilibrium; monotone energy decrease drives it along the vertical mode (heuristic, consistent with N4 observations).
- Exact symmetry nulls (ξ_f) give λ = 0, i.e. neutral eigenvalue 1. ∇L ⊥ ξ_f exactly on the position-gauged set (gauge file, section 1), so the descent step has no ξ_f component at the current iterate. That descent therefore does not drift along ξ_f is **heuristic**: the direction is linearly neutral, not damped, and only near-null discretely. This is why descent codes can work without a relabelling gauge while Newton cannot.
- Simultaneous descent in (x, Φ̃) on the full Lagrangian would also diverge in Φ̃: the co-energy is maximised there. Any scheme must eliminate or maximise the vacuum.

## 6. Linear rate of descent at a minimum (proved)

- If S ≻ 0 on the gauged space (S_red ≻ 0), with generalised eigenvalues 0 < λ₁ ≤ … ≤ λ_N for (S, P), then the error contracts asymptotically by ρ(τ) = max_i|1 − τλ_i|. This is minimised at τ* = 2/(λ₁ + λ_N), where ρ* = (κ − 1)/(κ + 1), κ = λ_N/λ₁.
  - Proof: |1 − τλ| is maximised at the spectrum ends, and the two ends balance at τ*.
- Reducing the error by 10^{-d} needs k ≥ d ln 10 / ln((κ + 1)/(κ − 1)) ≈ (d ln 10/2)κ iterations.
- Heavy-ball / second-order Richardson (VMEC's scheme is of this type) gives (√κ − 1)/(√κ + 1) with optimal parameters ([Polyak 1964](https://doi.org/10.1016/0041-5553(64)90137-5), **cited**).
- Bound on κ (repair 2026-10-19). Let v_R, v_Z be the unit rigid R and Z shifts (with shell, Φ̃ and λ re-solved; they lie in the gauged space by gauge file section 2 i). The generalised Rayleigh quotients vᵀSv/vᵀPv lie in [λ₁, λ_N]. The bound κ ≥ K_R/|K_Z| therefore needs two hypotheses: v_RᵀPv_R = v_ZᵀPv_Z (true for P = I on nodal coordinates, false for a general preconditioner P), and the vertical mode stable (hypothetically, at the same softness). In general κ ≥ (K_R/|K_Z|)·(v_ZᵀPv_Z / v_RᵀPv_R) (**verified**).
  - Illustration (numerical, fixture with the vertical mode hypothetically stable at the same softness, P = I): K_R/|K_Z| = 1.5e7/2.7e4 ≈ 556 gives ≈ 6396 plain-descent and ≈ 271 heavy-ball iterations for d = 10 (**verified**, recomputed 2026-10-19 from [ideal-free-boundary.md](../../docs/ideal-free-boundary.md), lines 82–84).
  - Even a perfect preconditioner of all other modes cannot remove this ratio, unless P itself carries the vertical stiffness.

## 7. Picard / frozen-field boundary updates (theorem cited; application hypothesis)

- A boundary update x_{k+1} = Φ(x_k) converges locally linearly, with asymptotic rate ρ(DΦ(x*)), if ρ < 1. It is a point of repulsion if ρ > 1 ([Ortega and Rheinboldt 2000](https://doi.org/10.1137/1.9780898719468), Ostrowski's theorem, §10.1; **cited**, content not verified: the point-of-repulsion clause and section number were not read).
- If Φ is a relaxation of the true gradient, Φ = x − τP⁻¹∇L_red, then section 5 applies: no convergence at a saddle.
- If Φ freezes the external B² while pushing the boundary toward pressure balance, then DΦ = I − τP⁻¹(S + E). E is the omitted response of the vacuum field to the boundary displacement.
  - S + E is not symmetric in general, so complex eigenvalues (oscillating iterates) and ρ ≥ 1 are both possible.
  - The fixed points coincide with equilibria only where the frozen field is updated consistently.
- The controller reports that upstream GVEC's free-boundary update is of this frozen, non-variational type and oscillates. This was not re-verified here: the local GVEC checkout (`/home/ert/code/gvec`, `bd52d45e`) contains no free-boundary source. The claim "its DΦ has eigenvalues near or outside the unit circle" is a **hypothesis**. The test: assemble DΦ by finite differences at a converged state and compute its spectrum.

## 8. Summary for the algorithm

- Pure Newton on the slide-gauged KKT system is locally quadratic at saddles and minima alike, provided R has full row rank and the gauged KKT block is nonsingular (equivalently S_red nonsingular when A_h + S_h ≻ 0). For the globalised solver this also needs the full-step-acceptance lemma (open).
- Descent and frozen Picard updates converge at best linearly, at a rate set by the stiffness ratio, and, locally, only from a positive-codimension set at the fixed-current vertical saddle.
- A per-state Kantorovich certificate is the planned discrete claim. A continuous certificate remains open.
