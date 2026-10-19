# Free boundary: first variation of the plasma–shell–exterior Lagrangian

- Scope: static ideal MHD, scalar pressure, nested plasma region Ω (solid torus) bounded by the material surface Γ_p, vacuum shell V between Γ_p and a fixed closed surface Γ_c, coils outside Γ_c. SI units. Smooth (C²) data up to Γ_p; the Sobolev version is not claimed.
- Geometry is chart-abstract: the map is X = H∘q with H an `ideal_chart_t` jet (cylindrical, or the curve frame H = Γ(q3) + q1 N(q3) + q2 B(q3)). Only X, DX, det DX and dH appear. The exterior is an abstract Steklov–Poincaré operator S on Γ_c; per-n toroidal modes are one realisation.
- Owning specification: [free-boundary-theory.md](../../docs/free-boundary-theory.md), obligations FB-A1, FB-A2, FB-B1. Implementation: [ideal-free-boundary.md](../../docs/ideal-free-boundary.md).
- Status labels: **proved** (complete argument here), **cited** (verified primary reference, not reproduced), **hypothesis**.

## 1. Data and unknowns

- Reference domain Ω̂ = Ω̂_p ∪ Ω̂_s with internal surface Γ̂_p and outer surface Γ̂_c. X: Ω̂ → R³ is continuous, cellwise smooth, J = σ det DX > 0. X is fixed on Γ̂_c.
- Frozen plasma data on Ω̂_p: a closed 2-form b̂ with b̂·n̂ = 0 on Γ̂_p, a label χ̂ with b̂(∇χ̂) = 0, pressure p(χ̂) (γ = 0, as in N1/N4). B = DX b̂ / J (push-forward), p∘X = p̂.
- Shell potential Φ = F τ₁ + c τ₂ + Φ̃ on V. τ₁, τ₂ are fixed closed 1-forms on the shell with periods (1, 0) on a toroidal and (0, 1) on a poloidal cycle; Φ̃ is single valued. The code uses τ₁ = dq3 (chart coordinate) and τ₂ = dθ̂/2π (reference angle), `src/mhd/ideal_free_boundary.f90:707`. Neither needs to be harmonic: Φ̃ absorbs the difference. F (2π F = toroidal circulation) and c (poloidal circulation) are fixed.
- Exterior (outside Γ_c): the field is B_e + B_i. B_e is the coil field (sources outside Γ_c, harmonic in the solid region D_c inside Γ_c). B_i = c B_f + ∇w has its sources inside Γ_c: B_f is a fixed unit-circulation field (a filament inside Γ_c), w is harmonic outside Γ_c and decays.
- Abstract exterior operator. For u ∈ H^{1/2}(Γ_c) let w_u be the decaying exterior harmonic extension and S u = −∂_n w_u (n outward from D_c). Then ⟨S u, v⟩ = ∫_{R³∖D_c} ∇w_u·∇w_v dV, so S is symmetric. In 3D, ⟨S u, u⟩ > 0 for u ≠ 0, because a constant trace has a nonzero capacity potential (**proved** by Green's identity on a large ball, with decay O(1/r) of w and O(1/r²) of ∇w).
  - Realisation R1 (implemented): axisymmetric Γ_c, physical cylindrical φ, per-mode operator L_n u = −div(R∇u) + n²u/R, symmetric Costabel coupling (`fortfem_toroidal_mode_boundary_operators_2d`, header lines 1–15; `src/equilibrium/equilibrium_boundary.f90:2–28`; `src/mhd/ideal_free_boundary.f90:435–479`). This realisation requires a global cylindrical angle on Γ_c and an axisymmetric Γ_c (**verified** against source).
  - Realisation R2 (not implemented): 3D Laplace BEM on a general closed Γ_c, S = W + (½ − K′)V⁻¹(½ − K), with V, K, K′, W the single-layer, double-layer, adjoint double-layer and hypersingular operators. It needs no chart.
- Coil drive on Γ_c: a = the single-valued trace of the coil potential plus the filament potential, h = (B_e + c B_f)·n. Both are fixed data.

## 2. The Lagrangian at fixed coil currents

\[
L[X,\tilde\Phi]=W_p[X]-\frac1{\mu_0}E_v[X,\tilde\Phi],\qquad
W_p=\int_{\hat\Omega_p}\Big(\frac{|DX\,\hat b|^2}{2\mu_0 J}-\hat p\,J\Big)ds,
\]
\[
E_v=\int_V\frac{|\nabla\Phi|^2}2\,dV+\frac12\langle S u,u\rangle_{\Gamma_c}-\langle h,\tilde\Phi\rangle_{\Gamma_c},\qquad u=\tilde\Phi|_{\Gamma_c}-a .
\]

This is the N4 functional of [ideal-free-boundary.md](../../docs/ideal-free-boundary.md) ("Formulation"), written for an abstract S. Both shell terms are pulled back to Ω̂_s through X. The pulled-back density |cof(DX) g|²/(2J), g = DΦ, is the generated leaf `gen/gen_vacuum_shell.f90`.

## 3. Variation in Φ̃: vacuum magnetostatics (proved)

For δΦ̃ single valued:

\[
\delta_{\tilde\Phi}E_v=\int_V\nabla\Phi\cdot\nabla\delta\tilde\Phi\,dV+\langle S u-h,\delta\tilde\Phi\rangle_{\Gamma_c}
=-\int_V\Delta\Phi\,\delta\tilde\Phi\,dV+\oint_{\Gamma_p}\partial_{n_V}\Phi\,\delta\tilde\Phi\,dS+\oint_{\Gamma_c}(\partial_n\Phi+Su-h)\,\delta\tilde\Phi\,dS .
\]

- Stationarity for all δΦ̃ gives ΔΦ = 0 in V, so B = ∇Φ is curl- and divergence-free.
- On Γ_p it gives B·n = 0 as a **natural** condition: the flux does not cross the material boundary, and no constraint is imposed.
- On Γ_c it gives ∂_nΦ = h − S u = (B_e + cB_f)·n + ∂_n w_u, the continuity of B·n with the exterior field B_e + cB_f + ∇w_u.
- Tangential continuity holds because the traces agree: Φ|_{Γ_c} = a + u, the coil and filament potentials plus w_u.
- E_v is strictly convex in Φ̃: ∫_V|∇Φ̃|² has the constants as its kernel, and ⟨S·,·⟩ is positive on constants. So Φ̃* is unique, and L is strictly concave in Φ̃ (a saddle, maximised in Φ̃). The discrete vacuum block is negative definite, as stated in the N4 doc.
- Chart independence: only V, ∇Φ and S enter. τ₁ need not be harmonic.

## 4. Variation of the plasma map: force balance and boundary traction (proved)

Let X_t = (I + tξ)∘X with ξ a smooth displacement field and ξ = 0 on Γ_c. Frozen-in transport gives the Lagrangian variations ΔB = (Dξ)B − (∇·ξ)B, Δp = 0 (γ = 0) and Δ(dV) = (∇·ξ)dV. Then

\[
\delta\!\int_\Omega\frac{|B|^2}{2\mu_0}dV=\frac1{\mu_0}\int_\Omega\Big(B\cdot(D\xi)B-\frac{|B|^2}2\nabla\cdot\xi\Big)dV,\qquad
\delta\!\int_\Omega(-p)\,dV=-\int_\Omega p\,\nabla\cdot\xi\,dV,
\]

so δW_p = ∫_Ω T:Dξ dV with T = (B⊗B − ½|B|²I)/μ0 − pI. Integrating by parts, with ∇·(B⊗B − ½|B|²I) = (∇×B)×B for ∇·B = 0:

\[
\delta W_p[\xi]=-\int_\Omega\xi\cdot(J\times B-\nabla p)\,dV-\oint_{\Gamma_p}\Big(p+\frac{|B|^2}{2\mu_0}\Big)\xi\cdot n\,dS ,
\]

where B·n = 0 makes Tn = −(p + |B|²/2μ0)n on Γ_p. This is the D3 Maxwell-stress residual ([equilibrium-core.md](../../docs/equilibrium-core.md) D3). For γ > 0 with p ∝ J^{-γ} and energy p/(γ−1), Δp = −γp∇·ξ gives the same −∫p∇·ξ ([grad-shafranov-variational-bridge.md](grad-shafranov-variational-bridge.md)).

## 5. Shape derivative of the vacuum co-energy (proved)

Let e(X) = E_v[X, Φ̃*(X)]. Φ̃* is stationary, so (envelope argument) de = ∂_X E_v at a fixed pulled-back potential. Transport the shell by a velocity θ (the shell map) with θ = ξ on Γ_p and θ = 0 on Γ_c. The exterior terms depend only on the trace on the fixed Γ_c and do not change. The Dirichlet energy at a transported potential has the classical derivative

\[
\frac{d}{dt}\int_{V_t}\frac{|\nabla\Phi_t|^2}2dV=-\int_V T_v:D\theta\,dV,\qquad T_v=\nabla\Phi\otimes\nabla\Phi-\tfrac12|\nabla\Phi|^2I .
\]

Here ∇·T_v = (ΔΦ)∇Φ = 0, so

\[
de[\theta]=-\oint_{\partial V}\theta\cdot T_v n_V\,dS=\oint_{\Gamma_p}\frac{|B_v|^2}2\,\theta\cdot n_V\,dS=-\oint_{\Gamma_p}\frac{|B_v|^2}2\,\xi\cdot n\,dS .
\]

The second equality uses ∂_{n_V}Φ = 0 on Γ_p, and n_V = −n is the shell normal there. The closed forms τ₁ and τ₂ are transported as forms, so their periods F and c do not change.

- Only ξ·n on Γ_p remains. The interior values of the shell map are a continuum gauge: the energy does not depend on them. The discrete shell energy does, which is why N4 imposes the harmonic-extension constraint and re-solves its multipliers ([ideal-free-boundary.md](../../docs/ideal-free-boundary.md), "Solver revision", item 2).

## 6. Result (proved)

\[
\delta L[\xi]=-\int_\Omega\xi\cdot(J\times B-\nabla p)\,dV-\oint_{\Gamma_p}\xi\cdot n\Big[p+\frac{|B|^2}{2\mu_0}-\frac{|B_v|^2}{2\mu_0}\Big]dS .
\]

(X, Φ̃) is a stationary point of L, for all ξ vanishing on Γ_c and all δΦ̃, if and only if

1. J×B = ∇p in Ω (B·n = 0 on Γ_p holds by construction of the push-forward);
2. ∇×B_v = 0 and ∇·B_v = 0 in V, with B_v·n = 0 on Γ_p;
3. B_v·n and B_v tangential are continuous on Γ_c with the field of the fixed coils plus a decaying field of internal sources;
4. [[p + |B|²/2μ0]] = 0 on Γ_p.

- The map unknowns are q, with ξ = dH(q)δq. Stationarity in q is equivalent to stationarity in ξ wherever dH is invertible: on the whole image for the cylindrical chart away from R = 0, and within the injectivity radius for a curve frame (**hypothesis** per case, checked by det dH > 0).
- The tangential field may jump on Γ_p (a surface current) while 4 holds. The vacuum circulation c is fixed data, and the interior circulation follows from X, so a mismatch is an admitted surface current rather than an inconsistency.
- Literature: frozen-in energy principle, [Kruskal and Kulsrud 1958](https://doi.org/10.1063/1.1705884) (**cited**); free-boundary energy with vacuum, [Bernstein, Frieman, Kruskal and Kulsrud 1958](https://doi.org/10.1098/rspa.1958.0023) (**cited**). The argument above is reproduced in KIN6D notation and does not rely on them.

## 7. Fixed flux versus fixed current (hypothesis)

- An ideal wall at Γ_c, or fixed coil fluxes, give the true magnetic energy (minimised in the vacuum field) in place of the co-energy (maximised in Φ̃).
- The first variation, and hence the equilibria with matching data, coincide (Legendre duality between the two vacuum problems). The second variations differ by the response term. With fixed flux the n = 0 vertical mode can be wall-stabilised, so the Hessian can be positive where the fixed-current Hessian is a saddle.
- Not derived here: FB-A4.

## 8. Second variation (cited, not reproduced)

- At an equilibrium, the Hessian of the reduced functional X ↦ L[X, Φ̃*(X)] on displacements is the ideal δW = δW_F + δW_S + δW_V of [Bernstein et al. 1958](https://doi.org/10.1098/rspa.1958.0023), with δW_V the vacuum energy at fixed coil currents.
- The per-n vacuum Green's-function evaluation of δW_V is [Chance 1997](https://doi.org/10.1063/1.872380) (VACUUM).
- Identifying the discrete Hessian with the discrete δW_n is FB-A3 (open; D5/M5 of the equilibrium core).

## 9. Discrete counterpart and checks

- Discrete L_h: the same densities at Gauss points, cellwise P_k × Fourier map, S_h symmetrised, `src/mhd/ideal_free_boundary.f90:475–476`. Its gradient is the generated derivative of the same leaf, assembled by the chain rule through the chart jet; see FB-B1.
- `ideal_fb_derivatives` checks gradient against central differences and Hessian against gradient differences (numerical evidence, not proof).
- `gen/gen_vacuum_shell.f90:52` proves with the FortSym native zero test that the leaf reduces to R|∇Φ|²/2μ0 at the identity cylindrical map.
- Not checked: positivity of S_h (FB-A2). Symmetrisation hides any assembly asymmetry and needs a separate defect report.
