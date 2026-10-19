# Exact 3D free-boundary benchmark from an exact interior equilibrium

- Question: given an exact smooth 3D equilibrium (B, p) on a nested domain, can a flux surface S0 be the free boundary of an exact free-boundary equilibrium, driven by an external field with sources outside a fixed coupling surface Γ_c?
- Answer: **yes, conditionally.**
  - The construction below is proved, except for one analytic hypothesis (G4): the harmonic continuation of the virtual-casing field must reach Γ_c.
  - G4 is proved only for a thin unquantified shell around S0. Whether it reaches a given axisymmetric (realisation R1) or general (R2) Γ_c is open, and must be checked per case.
- Obligations FB-G1 to FB-G6 of [free-boundary-theory.md](../../docs/free-boundary-theory.md). Status labels: **proved**, **cited**, **hypothesis**. Units μ0 = 1 in the Landreman formulas; restore μ0 by B → B/√μ0.

## 0. The source equilibrium (verified against the paper)

[Landreman, arXiv:2609.26742v2](https://arxiv.org/abs/2609.26742), section 2 (HTML v2 read 2026-10-19), family 1. Choose 0 < ε < 1, a = √(1 + ε), b = √(1 − ε). Then

\[
B=\Big(\frac{2zx-(a/b)Fy}{s},\ \frac{2zy+(b/a)Fx}{s},\ 1-s\Big),\quad s=\frac{x^2}{a^2}+\frac{y^2}{b^2},\quad F=\sqrt{1-(1-s)^2-4z^2},
\]

on U_ε = {F² > 0} (eqs. 2–3).
- Flux label ψ = (x² + y² + 4z² + |B|² − 2 + ε²)/4, and p = p_a − 2ψ (eq. 5).
- (∇×B)×B = ∇p exactly (eq. 6). In field-line labels, ψ = (u + ε/2)² + v² ≥ 0 (eq. 16).
- Nested solid tori Ω_δ = {ψ ≤ δ} for 0 < δ < (1 − ε)²/4 (eqs. 18–19).
- ι = 2 on every surface (closed field lines, eq. 26). Two field periods; stellarator symmetric. ε → 0 gives the Solov'ev field Ψ0 = (R² − 1)²/4 + Z², I = √(1 − 4Ψ0) (eqs. 7–8).
- KIN6D's closed form in `test/test_nonnested_landreman.f90:2–20` writes the same field with X = x/a, Y = y/b. That the two are algebraically identical is not re-checked here.
- The sheared family (section 3) and the non-symmetric extension ([Issan et al., arXiv:2610.07304](https://arxiv.org/abs/2610.07304)) satisfy the same hypotheses (smooth, analytic, nested). The construction applies verbatim; their parameter ranges are not re-derived here.

Pick S0 = {ψ = δ}, Ω = Ω_δ, and set

\[
\tilde p=p-p(\delta)=2(\delta-\psi)\ \in[0,2\delta]\ \text{on }\Omega,\qquad \tilde p=0\ \text{outside}.
\]

∇p̃ = ∇p, so force balance is unchanged, and p̃ ≥ 0 (**proved**, from ψ ∈ [0, δ]). On S0, p̃ = 0 but ∇p̃ = −2∇ψ ≠ 0, and J ≠ 0 (finite edge pressure gradient and edge current). The map lane admits this. J·n = 0 on S0 because J·∇p = 0 and p′ ≠ 0.

## 1. Inside virtual casing (FB-G1, proved)

Let B_in be C¹ on Ω̄ with ∇·B_in = 0 and B_in·n = 0 on S0 (n outward), and J = ∇×B_in/μ0. Let G = B_in χ_Ω (zero outside). With ∇χ_Ω = −n δ_{S0} (distributions):

\[
\nabla\cdot G=-(B_{in}\cdot n)\,\delta_{S_0}=0,\qquad
\nabla\times G=\mu_0J\chi_\Omega+(B_{in}\times n)\,\delta_{S_0}.
\]

- G ∈ L²(R³) has compact support. A divergence-free, curl-free L² field on R³ vanishes (Fourier transform). So G is the Biot–Savart field of its curl:

\[
B_{in}\chi_\Omega=B_J+B_K,\qquad B_J={\rm BS}[\mu_0J\chi_\Omega],\quad B_K={\rm BS}[\mu_0K\delta_{S_0}],\quad \mu_0K=B_{in}\times n .
\]

- Hence B_K = B_in − B_J inside Ω and B_K = −B_J outside.
- B_x := B_K|_Ω is curl- and divergence-free in Ω: its sources are on S0.
- This is the virtual-casing principle in its interior form: [Shafranov and Zakharov 1972](https://doi.org/10.1088/0029-5515/12/5/009), 3D form in [Lazerson 2012](https://doi.org/10.1088/0741-3335/54/12/122002) (**cited**; the proof above is self-contained).
- Regularity: J is bounded and piecewise smooth, so B_J is C^{0,α}(R³), continuous across S0 ([Gilbarg and Trudinger 2001](https://doi.org/10.1007/978-3-642-61798-0), Lemma 4.1, Newton potential of a bounded density).

## 2. Exterior construction and exact pressure balance (FB-G2, proved conditional on G4)

**Hypothesis G4.** Let D_c be the solid torus bounded by Γ_c, with Ω̄ ⊂ D_c. B_x extends to a curl- and divergence-free field B̃_x on an open set containing D̄_c.

**Theorem G2.** Under G4, let B_out = B_J + B̃_x on the shell V = D_c ∖ Ω̄. Then:
1. ∇×B_out = 0 and ∇·B_out = 0 in V (B_J is curl-free outside Ω̄).
2. B_out = B_in on S0, as full vectors: B_J is continuous across S0, and B̃_x is continuous with B_x.
3. Hence B_out·n = 0 on S0, there is no surface current, and |B_out|² = |B_in|² on S0. So

\[
[\![\tilde p+|B|^2/2\mu_0]\!]_{S_0}=\tilde p|_{S_0}+\frac{|B_{in}|^2-|B_{out}|^2}{2\mu_0}=0 .
\]

4. With coil field B_e := B̃_x|_{D_c} (all sources outside Γ_c) and internal field B_i := B_J outside Γ_c (decaying, dipole-like), (B_in, p̃ on Ω; B_out on V) satisfies conditions 1–4 of the first variation ([free-boundary-first-variation.md](free-boundary-first-variation.md), section 6). It is therefore a stationary point of the continuous KIN6D Lagrangian with this drive.

**Converse and uniqueness (proved).**
- If any vacuum field B_out on V equals B_in on S0 as a vector, then B_out − B_J on V and B_x on Ω glue across S0 without jumps. The glued field is curl- and divergence-free in the sense of distributions, hence harmonic componentwise (Weyl's lemma). So G4 is necessary for a surface-current-free exterior.
- B_out is unique: two such fields differ by a gradient ∇Φ with zero Cauchy data on S0, which vanishes on the connected V by unique continuation for harmonic functions.
- Allowing a surface current with |B_out| = |B_in| but different tangential direction gives other, current-sheet solutions. They are excluded here.

## 3. Circulations and topology (FB-G3, proved)

- The shell V has H₁ ≅ Z² (toroidal and poloidal cycles). D_c, a solid torus, has H₁ ≅ Z (toroidal).
- Poloidal circulation of B_out = μ0 × (flux of J through a meridian section of Ω) = μ0 I_tor. It comes entirely from B_J. B̃_x is curl-free on D_c, which contains a spanning disk of every poloidal loop, so it contributes 0. This fixes c (with the sign convention of N4: c = −μ0 I_p for current along +φ).
- Toroidal circulation of B_J on toroidal loops in V is 0.
  - The flux of the divergence-free J (tangent to S0) through Σ ∩ Ω, with Σ a disk spanning the loop, is a relative cycle with zero intersection with the core circle.
  - Axisymmetric check: B_{J,φ} = 0 outside Ω.
- Hence 2πF := toroidal circulation of B̃_x = ∮_γ B_in·dl, for γ ⊂ S0 in the class (1 toroidal, 0 poloidal). At ε = 0, F = I(δ) = √(1 − 4δ). F is the TF-coil content of the drive, the current through the hole of D_c, which must avoid the z-axis.
- The single-valued remainder of B̃_x on D_c is ∇χ_x. Its traces on Γ_c are the coil drive (a, h) of the first-variation file. For R1 they are per-n Fourier traces in the physical φ: the Landreman field has two field periods, so n ∈ 2Z; with nfp = 2 the drive has k = n/2. For R2 they are plain surface traces.

## 4. The continuation hypothesis G4 (partly proved, mostly open)

- Proved: S0 is a real-analytic surface. ψ is a real-analytic function on U_ε (elementary functions with F > 0, s > 0), and ∇ψ ≠ 0 on S0 (paper: ∇ψ vanishes only on the axis). The Cauchy data of the exterior potential (∂_nΦ = 0, tangential gradient = B_in,t) are real-analytic.
- By Cauchy–Kovalevskaya, a harmonic Φ_out with these data exists in a neighbourhood of S0. By Holmgren's theorem it is unique ([John 1982](https://doi.org/10.1007/978-1-4684-9333-7), ch. 3; **cited**). So G4 holds for every Γ_c in some shell {0 < dist(x, S0) < r₀}, with r₀ > 0 not quantified. This does not guarantee an axisymmetric Γ_c (R1) enclosing the non-axisymmetric S0 inside that shell; for small ε it is plausible (**hypothesis**).
- The singularities of the continuation are the "coils" of the benchmark. If they lie between S0 and the chosen Γ_c, the case is inadmissible. This is the same coil-inside-Γ_c violation the N4 RMP stage met ([ideal-free-boundary.md](../../docs/ideal-free-boundary.md), "RMP stage").
- Proposed evidence route (numerical, not proof): represent B̃_x in D_c by sources on a surface outside Γ_c (method of fundamental solutions, or toroidal harmonics). Fit it to B_x = B_in − B_J at interior points of Ω, with B_J by high-order Biot–Savart quadrature. Monitor coefficient decay and the fit residual. Report the drive defect of section 5.
- The ε = 0 (Solov'ev) case reduces to a 2D Cauchy problem for Δ*ψ_out = 0, the cheapest first check.
- A quantitative r₀ (Cauchy–Kovalevskaya majorants or a validated harmonic-continuation enclosure) would make G4 a theorem per case. **Open.**

## 5. Approximate drive and computable defect (FB-G6, proved as a residual statement)

- For any curl-free B̃_x^h on D_c (e.g. the fit of section 4), let d = B_J + B̃_x^h − B_in on S0.
- The Landreman state then satisfies the KIN6D stationarity conditions up to:
  - a flux leak (B_out^h·n = d·n on S0);
  - a traction defect (|B_in|² − |B_out^h|²)/2μ0 = −B_in·d_t/μ0 + O(|d|²).
- Both are computable from closed forms and quadrature. Converted into an error bound on the discrete solution, the defect gives ≲ ‖F(x_Landreman)‖/β_h. This requires the β_h of FB-E2 (**open**). Without it the benchmark is a residual test, not an error oracle.

## 6. KIN6D data for the benchmark (FB-G5)

| KIN6D datum | Value | Status |
|---|---|---|
| frozen b̂, χ̂ on the reference plasma domain | pull-back of B_in, ψ through the chosen reference map | proved (exact data); needs a 3D frozen reference (N1/N3), see gap list |
| pressure per label | p̃ = 2(δ − ψ) | proved |
| poloidal circulation c | −μ0 I_tor (flux of J through a section) | proved |
| toroidal constant F | (1/2π)∮_γ B_in·dl, γ ⊂ S0 toroidal | proved; ε = 0: √(1 − 4δ) |
| drive (a, h) on Γ_c | traces of ∇χ_x = B̃_x − F τ₁ | conditional on G4 |
| exact target | X* with X*_* b̂ = B_in, X*(Γ̂_p) = S0 | proved stationary (G2) |

- Nondegeneracy (**hypothesis, at risk** for family 1). ι = 2 is constant and the field lines are closed. Perturbations constant along field lines make interchange-type terms of δW degenerate, so S|_{ker M} may have additional (near-)zero modes. The bordered Newton system would then be singular there (Lemma B1), and the Kantorovich certificate unavailable.
  - Use family 1 as a stationarity and residual test (start at the target, expect zero residual to discretisation order).
  - Use the sheared family or Issan et al. for offset-start return and certificates.
- Physical stability is not required. The fixed-"coil" state may be vertically or otherwise unstable; Newton converges to it if it is nondegenerate (gauge file, section 5).
- The ε = 0 case is an axisymmetric Solov'ev free boundary with an exactly known interior. It also tests the ψ lane (`gs_free`) against the same drive: a cross-lane check (D8, item 4).

## 7. What is not claimed

- No claim that a physically realisable coil set produces B̃_x.
- No claim for a given Γ_c until G4 is checked.
- No claim that the discrete N4 solution equals the Landreman state beyond discretisation and drive-defect order.
- No uniqueness of the free-boundary solution for the given drive: other equilibria with the same frozen data and drive may exist.
