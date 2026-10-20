# Free boundary: relabelling gauge, bordered Newton system and the vertical saddle

- Scope: the map-lane free-boundary Lagrangian L of [free-boundary-first-variation.md](free-boundary-first-variation.md), γ = 0, nested plasma region Ω with label χ (χ_a on the axis, χ_b on Γ_p), chart-abstract (q3 is the chart's third coordinate: the physical φ for the cylindrical chart, the curve parameter for a curve frame).
- Obligations FB-C1 to FB-C6 and FB-D1 to FB-D4 of [free-boundary-theory.md](../../docs/free-boundary-theory.md). Status labels: **proved**, **cited**, **hypothesis**, **open**. Repairs from the independent proof audit of 2026-10-19 are marked **verified** (re-derived in the audit and the repair) or **hypothesis**.
- Standing hypotheses (H0): B, χ smooth on Ω̄; ∇χ ≠ 0 off the axis; each level set of χ is a connected torus; the axis is nondegenerate (χ − χ_a ≍ r², r the distance to the axis); B·∇q3 > 0 on Ω (field lines advance monotonically in q3). H0 is checked pointwise per case; for a curve frame, B·∇q3 > 0 is a **hypothesis** to be checked on the frozen state.

## 1. The relabelling symmetry (proved)

- Let Ψ_t be the flow of a reference vector field v tangent to Γ̂_p. If Ψ_t preserves the frozen data, Ψ_t^* b̂ = b̂ and χ̂∘Ψ_t = χ̂ (p̂ = p(χ̂)), then W_p[X∘Ψ_t] = W_p[X] for every X. (The converse is neither proved nor used; repair 2026-10-19.)
  - Reason: B = X_* b̂ and p∘X = p̂ are then unchanged, and W_p is the integral of an X-independent function of (B, p) over X(Ω̂).
- The vacuum part depends on X only through Γ_p = X(Γ̂_p) and the shell map. Γ_p is invariant, because v is tangent to Γ̂_p. The shell map is a continuum gauge (first-variation file, section 5). So L is invariant in the continuum.
- Infinitesimally, using db̂ = 0 and Cartan's formula: L_v b̂ = d(i_v b̂) = 0 and v·∇χ̂ = 0. By naturality of the push-forward, the physical displacement ξ = X_* v = DX v satisfies the same at X:

\[
\xi\cdot\nabla\chi=0,\qquad \nabla\times(\xi\times B)=0 ,
\]

  that is δB = 0 and δp = −ξ·∇p = 0. The set 𝔤_X of such ξ is the relabelling algebra.
- For γ > 0, the frozen mass 3-form adds div(ρξ) = 0. This is not needed in N4, which uses γ = 0.

### Position gauge leaves exactly one function of χ (proved under H0)

- N4 moves nodes only in the chart's (q1, q2) directions: ξ = dH(δq1, δq2, 0), i.e. ξ·∇q3 = 0 with ∇q3 the dual covector of ∂_{q3}H.
- Let ξ ∈ 𝔤_X with ξ·∇q3 = 0. Since ξ ⊥ ∇χ and ξ ⊥ ∇q3, and these are independent off the axis, ξ = α ∇χ×∇q3.
- Then (∇χ×∇q3)×B = ∇q3(B·∇χ) − ∇χ(B·∇q3) = −(B·∇q3)∇χ, so ξ×B = −g∇χ with g = α B·∇q3.
- The condition ∇×(ξ×B) = −∇g×∇χ = 0 holds iff dg∧dχ = 0, i.e. g is constant on each connected level set: g = f(χ). Hence

\[
\mathfrak g_X\cap\{\xi\cdot\nabla q_3=0\}=\Big\{\xi_f=f(\chi)\,\frac{\nabla\chi\times\nabla q_3}{B\cdot\nabla q_3}\Big\}.
\]

- Axisymmetric cylindrical case: ∇χ×∇φ = ∇^⊥χ/R and B·∇φ = B_φ/R, so ξ_f = f∇^⊥χ/B_φ with b3 = B_φ (the toroidal-flux density in dR dZ). This is the null of `slide_directions` (`src/mhd/ideal_free_boundary.f90:1345–1360`) and of [ideal-free-boundary.md](../../docs/ideal-free-boundary.md), "Relabelling (slide) gauge".
- Consequences:
  - The general 3D null (stellarator, figure-8 chart) is still one function of χ, but ξ_f has every reference Fourier mode when B and χ do. The present n = 0-only gauge is the axisymmetric special case (gap FB-C6).
  - The parallel relabelling ξ = λB is in 𝔤_X, since i_B(i_B dV) = 0. The position gauge removes it because λB·∇q3 ≠ 0 under H0.

### Invariance identities (proved)

- Let P_g = {X : q3∘X = q̂3 on plasma nodes} be the position-gauged set. The reference generator v_f = (DX)⁻¹ξ_f of ξ_f is X-independent only on P_g, because ∇q3 pulls back to the fixed dq̂3 exactly when q3∘X = q̂3. Hence L[X∘Ψ_t] = L[X] for all X ∈ P_g and t (restriction **verified**, repair 2026-10-19; for X ∉ P_g the flow Ψ_t still preserves L, but its push-forward is not ξ_f). Differentiating at t = 0:

\[
\langle\nabla L(X),\xi_f\rangle=0\quad\forall f,\ \forall X\in P_g .
\]

  The gradient is orthogonal to the orbit on all of P_g, not only at equilibria.
- Differentiating again in X at a critical point: H(X*)ξ_f = 0. The relabelling null is an exact Hessian kernel of the continuum problem. Discretely it is only near-null, because ξ_f is not an exact symmetry of the quadrature and P_k map.

### Position-gauged stationarity is full stationarity (proved; repair 2026-10-19)

- Every smooth ξ on Ω splits as ξ = ξ′ + λB with λ = ξ·∇q3/(B·∇q3) (smooth, B·∇q3 > 0 under H0) and ξ′·∇q3 = 0. In the shell the map is not position-gauged, so no splitting is needed there.
- For every X, at Φ̃ = Φ̃*(X): ⟨∇L, λB⟩ = 0. By the first variation (first-variation file, section 6), δL[λB] = −∫λB·(J×B − ∇p) dV − ∮λ(B·n)[…] dS. Here B·(J×B) = 0, B·∇p = p′(χ)B·∇χ = 0 (B = X_*b̂ and b̂(∇χ̂) = 0 for every X), and B·n = 0 on Γ_p. The vacuum part depends only on the normal trace ξ·n on Γ_p, which is 0 (**verified**).
- Hence, if (X, Φ̃) is stationary for all single-valued δΦ̃ and all position-gauged ξ′, it is stationary for all ξ, and FB-A1 applies. This is the step used by section 3 and by the first-variation result.

## 2. The slide gauge functional (proved)

For Λ ∈ C¹([χ_a, χ_b]) with Λ(χ_b) = 0, define

\[
M_\Lambda(\xi)=\int_\Omega\Lambda(\chi)\,\nabla q_3\cdot(\nabla\times\xi)\,dV .
\]

∇q3 is single valued although q3 is multivalued.

(i) **Curl-free displacements are free.** M_Λ(ξ) = 0 for every ξ with ∇×ξ = 0 in Ω, trivially. This covers translations, gradients of potentials (elongation (a x, −a z)-type strains) and the toroidal harmonic field. A rigid poloidal rotation (curl = 2ω) is not free, which is correct: it is not a symmetry of the problem with coils. In 3D, a curl-free motion with ξ·∇q3 ≠ 0 is represented in the position gauge by ξ − λB (section 1), which is not curl-free in general. The gauge stays complete by (iv), but "translations need no relabelling" then holds only modulo an added ξ_f (**verified** note, audit 2026-10-19; relevant to FB-C6).

(ii) **Derivative-free form.** With ∇·(ξ×A) = A·∇×ξ − ξ·∇×A, A = Λ∇q3 and ∇×A = Λ′∇χ×∇q3:

\[
M_\Lambda(\xi)=\int_\Omega\Lambda'(\chi)\,(\nabla\chi\times\nabla q_3)\cdot\xi\,dV+\oint_{\Gamma_p}\Lambda(\chi_b)\,(\xi\times\nabla q_3)\cdot n\,dS ,
\]

and the boundary term vanishes. In the cylindrical n = 0 case, (∇χ×e_φ) = ∇^⊥χ and dV = R dR dZ dφ give M_Λ = 2π∫ Λ′ ∇^⊥χ·u dR dZ. This is the code's M_l with φ_l = Λ_l′ (`src/mhd/ideal_free_boundary.f90:1350–1356`), up to sign and 2π. The code's claim "M_l vanishes for every curl-free u" is therefore **proved** in the continuum.

(iii) **Transversality.** For ξ_f, by (ii) and the coarea formula,

\[
M_\Lambda(\xi_f)=\int_{\chi_a}^{\chi_b}\Lambda'(\chi)f(\chi)\,w(\chi)\,d\chi,\qquad
w(\chi)=\oint_{\chi}\frac{|\nabla\chi\times\nabla q_3|^2}{(B\cdot\nabla q_3)\,|\nabla\chi|}\,dS>0 .
\]

For the code's partition φ_l ≥ 0 (C¹, one per element layer, Λ_l(χ) = −∫_χ^{χ_b}φ_l), the gauge matrix A_lk = M_{Λ_l}(ξ_{φ_k}) = ∫φ_lφ_k w dχ is a Gram matrix. It is symmetric positive definite whenever the φ_l are linearly independent. The code normalises its directions and uses a Cholesky factor of the analogous matrix (`project_slide`, `src/mhd/ideal_free_boundary.f90:1430–1470`).

(iv) **Exact continuum decomposition.** Let G_∞ = {ξ : M_Λ(ξ) = 0 ∀Λ}. Define m_ξ(χ) = ∫_{χ′<χ}∇q3·∇×ξ dV. Up to sign, this equals the circulation-type quantity ∮_χ(∇χ×∇q3)·ξ/|∇χ| dS, by (ii) applied to the step Λ = 1_{χ′<χ} (limit of C¹ functions). Then:
- G_∞ = {ξ : m_ξ ≡ 0}.
- Every smooth ξ splits uniquely as ξ = ξ_f + g with f = m_ξ/w and g ∈ G_∞. Uniqueness: ξ_f ∈ G_∞ ⇒ fw = 0 ⇒ f = 0.
- At the axis, m_ξ = O(χ − χ_a) (volume of a thin tube times a bounded integrand) and w ≍ (χ − χ_a) under H0. So f stays bounded. Smoothness of f across the axis and the Sobolev version are a **hypothesis**.
- Hence G_∞ removes exactly the relabelling null and leaves every curl-free displacement free (**proved**, continuum, smooth ξ).

(v) **Finite gauge.** The code imposes only L + 1 moments. They fix the w-weighted projection of f onto span{φ_l}; f orthogonal to it stays null in the continuum. Discretely those f oscillate within one element layer. The claim that they are not near-null of the discrete Hessian ("stiff") is a **hypothesis** with numerical evidence only: all six fixture starts converge to ≤ 4.4e-11 σV ([ideal-free-boundary.md](../../docs/ideal-free-boundary.md), lines 285–291). The discrete curl-free defect ratio 1.9e-3 → 5.4e-4 under refinement (line 284) is consistent with (i) and is not a proof.

## 3. The gauged stationary point is an equilibrium (proved, continuum)

The bordered system solves ∇L(X) + Σ_l μ_l ∇M_l = 0 on position-gauged variations, and M_l(X − X₀) = 0, with X ∈ P_g. Pair the first equation with ξ_{φ_k}:
- ⟨∇L, ξ_{φ_k}⟩ = 0 by the invariance identity (X ∈ P_g);
- ⟨∇M_l, ξ_{φ_k}⟩ = A_lk(X). A is a Gram matrix only at X₀, where the moments are fixed functionals built from X₀ (section 2 iii). For X ≠ X₀ it is nonsingular by continuity only, in a neighbourhood of X₀ whose radius is not quantified (**proved**, local; repair 2026-10-19).

Hence, for X in that neighbourhood, μ = 0 and ∇L(X) = 0 on position-gauged variations. By "Position-gauged stationarity is full stationarity" (section 1), ∇L(X) = 0 for all ξ. The result is local near X₀ and depends on that step (**verified**, audit 2026-10-19). The same argument applies to the code's projected form (`project_slide`): if the projected gradient Π g vanishes, then g = Σ μ_l s_l, so μA = 0 and μ = 0. The gauge multiplier vanishes at a solution, which justifies taking the merit with the multiplier removed ([ideal-free-boundary.md](../../docs/ideal-free-boundary.md), "Relabelling (slide) gauge"). Discretely, μ = O(discretisation) is a **hypothesis**.

## 4. Bordered matrix lemma (proved)

Let H ∈ R^{n×n} be symmetric, M ∈ R^{r×n} and K = [H Mᵀ; M 0].

**Lemma B1.** K is nonsingular iff rank M = r and ZᵀHZ is nonsingular, where Z is a basis of ker M.

*Proof.* (⇐) Let K(x, y) = 0. Then Mx = 0, so x = Zz, and ZᵀHZz = −ZᵀMᵀy = 0, so z = 0. Then Mᵀy = 0, and y = 0 by full row rank. (⇒) If rank M < r, there is y ≠ 0 with Mᵀy = 0, and K(0, y) = 0. If ZᵀHZz = 0 with z ≠ 0, then HZz ∈ (ker M)^⊥ = range Mᵀ, say HZz = −Mᵀy. Then K(Zz, y) = 0. ∎

**Lemma B2 (inertia).** If rank M = r, then In(K) = In(ZᵀHZ) + (r, r, 0).

*Proof.* Complete Z to a basis [Z Y] of Rⁿ, with C = MY ∈ R^{r×r} nonsingular. The congruence diag([Z Y], I) gives

\[
\begin{pmatrix}Z^THZ&Z^THY&0\\Y^THZ&Y^THY&C^T\\0&C&0\end{pmatrix}.
\]

The entries ZᵀHY and YᵀHY are removed by congruences with the C and Cᵀ blocks:
- add −ZᵀHY C⁻¹ times the third block row to the first block row, and symmetrically for the columns;
- remove YᵀHY with −½ YᵀHY C⁻¹ on the third block row, and symmetrically.

The result is diag(ZᵀHZ, [0 Cᵀ; C 0]). The second block has eigenvalues ±σ_i(C): r positive and r negative. Sylvester's law of inertia completes the proof. ∎

- Lean: B1 and B2 are checked in `proof/Kin6d/Bordered.lean` (`lemma_B1`, `lemma_B2`; real matrices, finite index types, inertia counted variationally and shown equal to eigenvalue counts for symmetric matrices).

## 5. Vacuum elimination and the vertical saddle

Repaired 2026-10-19 (audit: first invalid bridge). The earlier text took x = (plasma map and multipliers), so its Schur complement was still bordered by the shell constraints, and neg(S|_{ker M}) exceeded the physical index by rank C.

- Unknowns: map coefficients q (plasma and shell nodes), shell-constraint multipliers λ (constraint rows ∇C), slide multipliers μ (rows M), and Φ̃. The constraint rows do not involve Φ̃, so

\[
K_{\rm full}=\begin{pmatrix}H_{qq}&H_{q\Phi}&\nabla C^T&M^T\\H_{\Phi q}&H_{\Phi\Phi}&0&0\\\nabla C&0&0&0\\M&0&0&0\end{pmatrix},
\]

  with H_qq the Hessian of the Lagrangian (including λ·∇²C).
- **Hypotheses (open, per state):** (V) H_ΦΦ ≺ 0, i.e. A_h + S_h ≻ 0 (first-variation file, section 3; FB-A2 discrete); (R) the stacked rows R = [∇C; M] have full row rank r = r_C + r_M on the map unknowns.
- Define the map Schur complement S_map = H_qq − H_qΦ H_ΦΦ⁻¹ H_Φq and Z a basis of ker R on the map unknowns. Under (V), Haynsworth inertia additivity ([Haynsworth 1968](https://doi.org/10.1016/0024-3795(68)90050-5), **cited**; it follows from the same congruence argument with a nonsingular pivot block) eliminates Φ̃. Under (R), Lemma B2 applied to the bordered Schur complement gives

\[
\mathrm{In}(K_{\rm full})=\mathrm{In}(H_{\Phi\Phi})+\mathrm{In}(Z^TS_{\rm map}Z)+(r,r,0),\qquad \mathrm{In}(H_{\Phi\Phi})=(0,\dim\tilde\Phi,0).
\]

  (**verified**, repair 2026-10-19: the bordered Schur complement in (q, λ, μ) has exactly the form of Lemma B2 with H = S_map and M = R.)
- S_red := ZᵀS_map Z is the reduced plasma Hessian at fixed coil currents on the gauged space (shell map eliminated by its constraints, Φ̃ eliminated, slide gauge imposed): the discrete δW in the sense of FB-A3. The physical Morse index is neg(S_red). If R has full rank this equals neg(S_b|_{ker M}) − r_C, where S_b is the Schur complement still bordered by ∇C.
- Indefiniteness of K_full says nothing about physical stability: it always has dim Φ̃ + r negative eigenvalues besides neg(S_red).
- Vertical saddle (numerical evidence, not proof):
  - fixture: K_Z = −2.7e4 against K_R = +1.5e7 J/m² per period;
  - ITER: K_Z = −5.4e6 against K_R = +2.5e8;
  - n = 0 shifted inverse iteration finds the vertical negative eigenvalue (−3.9e-7, ITER r0).
  
  So neg(S_red) ≥ 1 at fixed currents ([ideal-free-boundary.md](../../docs/ideal-free-boundary.md), lines 80–86 and NL-3a). K_R and K_Z re-solve shell, Φ̃ and λ for a rigid shift, so they are Rayleigh numerators of S_red (rigid shifts are curl-free and satisfy M = 0 by section 2 i). That neg(S_red) = 1 at the n = 0 fixture states (no other n = 0 negative mode) is a **hypothesis**: the probe found one.
- Nonsingularity condition (proved, conditional on (V)): K_full is nonsingular iff R has full row rank and S_red is nonsingular. Proof: Lemma B1 on the whole KKT matrix with constraint rows [R 0] gives "full rank and diag(Z, I)ᵀ H diag(Z, I) nonsingular"; under (V) the latter is nonsingular iff its Schur complement S_red is. Without (V), only the Lemma B1 form on the whole KKT matrix holds. The negative vertical eigenvalue is harmless.
  - Singularity occurs exactly when S_red has a zero eigenvalue. Along a continuation path this need not be a sign crossing. It includes marginal vertical stability (K_Z → 0), folds of the equilibrium branch, branch points (e.g. the up-down-symmetric pitchfork at K_Z = 0) and any other marginal mode (**verified** correction, repair 2026-10-19).
  - The eigenvalues depend continuously on the parameters, so the Morse index is constant along a path on which S_red stays nonsingular.
- Conditioning (proved, linear algebra): for the Newton correction, |δZ| ≤ ‖δx‖ ≤ ‖K⁻¹‖‖F‖ with ‖K⁻¹‖₂ = 1/σ_min(K). Since σ_min ≤ |λ_v| for the vertical eigenvalue λ_v, ‖F‖/|λ_v| is a **lower** estimate of this bound. It is the bound only if the vertical mode is the smallest singular direction (**hypothesis**, per state). The observed "F_Z/K_Z is the vertical uncertainty" uses |K_Z| as an estimate of |λ_v| (numerical; repair 2026-10-19).

## 6. Globalisation at a saddle (proved)

- An energy merit is invalid at a saddle: L decreases away from x* along the vertical mode.
- The residual merit φ = ½‖F‖²_D, in a fixed scaling D, is valid. The Newton direction d = −J⁻¹F satisfies ⟨∇φ, d⟩ = ⟨JᵀD²F, −J⁻¹F⟩ = −‖F‖²_D < 0 whenever J is nonsingular, so it is a descent direction for φ at every non-solution. Armijo backtracking on φ is therefore well defined.
- Caveat: the N4 merit Σ g²/|H_ii| uses a diagonal scaling D² = diag(1/|H_ii|). The descent property holds for any fixed SPD D. If D is updated with H every iteration, the merit changes between iterations. Acceptance must compare φ_k(x_k + αd) with φ_k(x_k) in the same D_k, which is a **hypothesis** to verify in `ideal_fb_solve`.

## 7. The vertical-coordinate restriction is not a gauge (proved)

- The "Coupled vertical Newton coordinate" solves Pᵀ∇L = 0 with P = [T v]: T holds the normal n = 0 node directions, v the rigid Z shift ([ideal-free-boundary.md](../../docs/ideal-free-boundary.md), lines 303–330).
- This is the stationarity of L restricted to range P. It implies ∇L = 0 only if ∇L has no component orthogonal to range P. The tangential components (including the net R force where the surface tangent is R) are discarded.
  - Counterexample (numerical, from the N4 doc): `tangential_gauge` stops the rp start with Γ_p moved −0.18 mm instead of −15 mm and a net rigid R force of 1.4e3 J/m.
- The slide gauge removes only ξ_f (section 2), and its multiplier vanishes at solutions (section 3).
- Recommendation: the production algorithm is slide-gauged Newton without `tangential_gauge`. The vertical row is then not needed for nonsingularity (B1). It can stay as a reported diagnostic (F_Z, K_Z), or as an explicit coordinate within the slide-gauged space. Within that space the system remains a true gauge plus a change of basis, if v is not in span(T) and the restriction is to the full slide-gauged space rather than to normals only.
