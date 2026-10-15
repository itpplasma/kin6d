# Grad–Shafranov, weak form and nested ideal-MHD variation

- Scope: static, scalar-pressure, axisymmetric equilibrium; smooth nested surfaces, fixed physical boundary, R>0 and SI units. Flow, anisotropy, islands and free-boundary vacuum work require separate derivations.
- Canonical field and signed flux: follow `grad-shafranov.md`; physical phi is CCW and psi is poloidal flux per radian.
- This file separates hand derivations, exact local algebra and unresolved solver correspondence. No identity count establishes all code modes or continuous convergence.

## Strong equation and physical force

Write B=(-psi_Z/R,F/R,psi_R/R) in the RH cylindrical basis (R,phi,Z), with F=F(psi). Cylindrical curl gives

\[
\mu_0 j=(-F'\psi_Z/R,-\Delta^*\psi/R,F'\psi_R/R),\quad
\Delta^*\psi=\psi_{RR}-\psi_R/R+\psi_{ZZ}.
\]

Therefore

\[
 j\times B-\nabla p=-\frac{\Delta^*\psi+FF'+\mu_0R^2p'}{\mu_0R^2}\nabla\psi.
\]

- The toroidal force vanishes identically. Arbitrary symbolic derivative jets check both remaining components independently of an exact equilibrium fixture.
- Away from critical points, force balance implies GS. At the magnetic axis, force vanishes with grad(psi), so dividing by that gradient is invalid; the scalar equation follows by smooth continuation from surrounding points.
- At a separatrix saddle, the physical GS operator remains regular for R>0. A nested flux-coordinate inverse loses rank; that is a different failure.

## Weak equation and boundary terms

\[
\frac{\Delta^*\psi}{R}=\nabla_{R,Z}\cdot\frac{\nabla_{R,Z}\psi}{R},\quad
\int_D\frac{\nabla\psi\cdot\nabla v}{R}\,dA
=\int_D\left(\mu_0Rp'+\frac{FF'}R\right)v\,dA.
\]

- Use prescribed psi trace and v in H_0^1(D). The unrestricted Green boundary term is integral_boundary v*(normal dot grad(psi))/R dl; it vanishes for zero test trace.
- For 0<R_min<=R<=R_max, the weighted principal bilinear form is positive/coercive on a bounded Dirichlet domain by Poincare. Boundary regularity, source regularity and nonlinear solvability remain separate hypotheses.
- Local product-rule algebra can be CAS checked. Integration by parts, trace theorems and admissible function spaces are analytic arguments, not CAS consequences.
- KIN P1/Gauss assembly, CHEASE Hermite/Gauss assembly and FreeGS centered differences discretize different representations of this equation. Matching the continuous equation does not identify their discrete residuals.

## Scalar stationary functional

\[
\mathcal G[\psi]=\int_D\left[\frac{|\nabla\psi|^2}{2R}-\mu_0Rp(\psi)-\frac{F(\psi)^2}{2R}\right]dA,
\]
\[
\delta\mathcal G[v]=\int_D\left[\frac{\nabla\psi\cdot\nabla v}{R}-\left(\mu_0Rp'+\frac{FF'}R\right)v\right]dA.
\]

- Fix D, the Dirichlet trace and the functions p(psi), F(psi). This variation gives the weak equation.
- The quadratic gradient term is positive; the source potential need not be convex. A stationary point is not automatically a minimum, unique branch or stable dynamical equilibrium.
- Second variation: integral_D [|grad(v)|^2/R -(mu0 R p''+((F')^2+FF'')/R)*v^2] dA. Constant p' and FF' give a positive Hessian for fixed Dirichlet trace; nonlinear laws require a separate bound.
- The negative F^2 term makes this a GS source functional, **not** magnetic energy with B substituted. It must not be equated directly with the nested-surface MHD energy used by VMEC/GVEC.

## Nested ideal-MHD energy to force balance

For a material displacement xi, preserve magnetic flux with delta B=curl(xi cross B). A fixed-boundary displacement makes the magnetic surface-work term vanish. The pointwise magnetic identity is

\[
\frac{B\cdot\delta B}{\mu_0}=\nabla\cdot\frac{(\xi\times B)\times B}{\mu_0}-\xi\cdot(j\times B).
\]

For advected entropy/mass with material pressure variation delta_L p=-gamma*p*div(xi), the internal-energy variation obeys

\[
\delta\int\frac p{\gamma-1}dV=-\int p\,\nabla\cdot\xi\,dV
=\int\xi\cdot\nabla p\,dV-\int_{\partial V}p\xi\cdot n\,dS.
\]

Thus, under the stated admissible variations and vanished boundary work,

\[
\delta W=-\int\xi\cdot(j\times B-\nabla p)\,dV.
\]

- Gamma=0 is regular: advected prescribed pressure has energy density -p. Gamma=1 needs its own logarithmic/isothermal functional; do not evaluate p/(gamma-1) there.
- For shell-volume constrained pressure p=M/(V')^gamma, the shell energy M*(V')^(1-gamma)/(gamma-1) has derivative -p with respect to V'. This establishes the pressure force in VMEC's constrained profile model.
- Flux-preserving geometry and lambda variations in a nested representation are restrictions of the ideal variations. Their stationarity supplies projected force balance; completeness of the variation space is an additional condition.
- Axisymmetric force reduction then gives GS. This bridge is through force balance, not an identification of W with the scalar functional G.
- Finite Fourier/radial trial spaces can hide unresolved force harmonics. Native projected stationarity is not an error bound on continuous physical force.
- Prescribed iota/current constraints, numerical de-aliasing constraints, regularity projections and free-boundary forces each need an explicit derivation before claiming equivalence for those branches.

## Signed profiles and flux constraints

- Let Phi be full signed toroidal flux, s=Phi/Phi_edge and psi signed poloidal flux per radian. For consistently oriented nested surfaces,

\[
\frac{d\Phi}{d\psi}=2\pi q,\quad \iota=1/q,\quad
\frac{d\psi}{ds}=\frac{\Phi_{edge}}{2\pi q(s)},\quad
p'(\psi)=\frac{2\pi q}{\Phi_{edge}}\frac{dp}{ds}.
\]

- Establish the flux relation by oriented field-line winding/coarea on each regular surface. Positive arc length with |Bpol| loses the poloidal-circulation sign; it cannot independently establish signed q.
- In rho=sqrt(s), dp/drho=2rho*dp/ds and d2p/drho2=2dp/ds+4rho^2*d2p/ds2. A cubic in s becomes a piecewise even sextic in rho.
- Constant-q input is an inverse constraint. Prescribing boundary, p, total toroidal flux and q requires releasing F/current; it is not generally the same problem as prescribed p'(psi), FF'(psi) and edge F.
- NCURR=0 prescribed iota and NCURR=1 prescribed enclosed-current profiles are different constraints. A q seed, total-current scaling and a point q target do not prescribe the full q function.
- Derivative normalization includes signed flux span: if x=(psi-psi_axis)/(psi_edge-psi_axis), then dp/dpsi=(dp/dx)/(psi_edge-psi_axis). Knots, units, offsets, branch sign and extrapolation are parts of the source law.
- Equal knot values do not imply equal between-knot derivatives, endpoints or continuous force. First establish literal functions and chain rules, then assess interpolation/derivative error separately.

## Code correspondence and closure obligations

| Lane | Executed structure | Established / still open |
|---|---|---|
| KIN6D | `src/equilibrium/grad_shafranov.f90` weighted P1 assembly; generated source/field/primitive leaves | Strong/source algebra and assembly inspected. Weak boundary theorem, quadrature/continuous error and nonlinear implicit derivatives are separate. |
| FreeGS | `gradshafranov.py::GSsparse`, `equilibrium.py::solve`, `jtor.py::ProfilesPprimeFfprime` | Native second-order operator/current law inspected. Curved-LCFS mask is external; signed q and boundary refinement remain open. |
| Public CHEASE | `setupa/setupb/curent/matrix/surface` | Hermite weak/Picard paths inspected with native normalization. Every source/constraint selector and full SI current correspondence are not proved. |
| MARS CHEASE | Monolithic `SETUPA/SETUPB/CURENT/SURFACE` | Shared algorithm ancestry with public CHEASE. Native exports/endpoint derivatives differ; no separate independence vote. |
| Original VMEC | `bcovar.f`, `profil1d.f`, `add_fluxes.f90`, force projections | Continuous nested field and pressure constraint inspected. Full discrete force variation/current-constrained branches remain open. |
| VMEC++ | `radial_profiles.cc`, `ideal_mhd_model.cc` and local force kernels | Same variational ancestry; exact output-map/profile functions must match native pin. Complete discrete constraints and all modes remain open. |
| DESC | `_field.py`, `_equil.py`, `ForceBalance`, profile classes | Force collocation differs from energy objective. Existing ten exact identities are local algebra, not full collocation/GS equivalence. |
| GVEC | `mhd3d_evalfunc.F90`, hmap and B-spline profiles | Nested MHD variation; discretization independent of VMEC. Documentation force-sign discrepancies are scoped; full discrete force/regularity derivation remains open. |

- Source pins/file anchors and literature differences live in writeups' individual correspondence ledgers. Their local tests do not certify other source commits or modes.
- Physical comparison obligations: common points, native pressure derivatives/current, signed flux and q, boundary refinement, off-grid force and derivative budgets. Solver iteration success is a different gate.
- This first slice excludes perturbation and NTV; its chart/constraint contracts are prerequisites for those later comparisons.

## Exact algebra evidence and proof limits

- `test/test_gs_variational_identities.f90` checks fourteen exact identities: three arbitrary-jet forces, divergence form, Green product rule, nonlinear polynomial functional variation, ideal magnetic/material pressure variation, two toroidal-flux/rho chain rules and four direct generated-KIN source/field expressions.
- Fresh standalone native CMake/CTest: 3/3 gates passed, including intentional wrong-field-sign and wrong-weak-source-weight rejection. Fourteen exact native decisions; zero randomized probes. The three gates are registered in the optional pinned FortSym configuration.
- Source reader scopes `gs_source` and `gs_field` from the actual generated file using `logical_lines/find_assignment`; no expression is copied from the producer. Statement/control-flow and floating rounding semantics are outside this expression check.
- Generated segment primitive readback is currently blocked in pinned FortSym a614271 by `real(2,kind=kind(0.0_dp))` parsing. The generator's existing exact primitive check is separate; do not claim source-readback closure for this leaf.
- Native exact zero decisions only; UNKNOWN/probe does not count as an identity. Require R!=0, mu0!=0, gamma!=1 and the declared chart/regularity assumptions.
- Polynomial differential fixtures test nontrivial local algebra with exact coefficients; they do not prove identities for all admissible functions without the accompanying calculus argument.
- No claim of Lean proof, integration by parts proof, solver existence/uniqueness, compiler correctness or complete third-party discrete equivalence.

## Source and literature anchors

- KIN6D inspected base: `5d75dfc6d9a43b0c3af1f61b44163a08f8c0a2a4`; generated source and exact FortSym pin are hashed in the TC24 closure receipt.
- [GS source ledger](https://github.com/itpplasma/writeups/blob/main/theory/gs_solver_source_correspondence.md): KIN/FreeGS/public and MARS CHEASE pins, routine paths and weak/source normalizations.
- [Original VMEC ledger](https://github.com/itpplasma/writeups/blob/main/theory/vmec2000_equilibrium_source_correspondence.md): native magnetic field, shell pressure energy, profile families and flux/current constraints.
- [DESC ledger](https://github.com/itpplasma/writeups/blob/main/theory/desc_equilibrium_source_correspondence.md): collocation residual, native signed flux, axisymmetric reduction and separate Energy objective.
- [GVEC ledger](https://github.com/itpplasma/writeups/blob/main/theory/gvec_equilibrium_source_correspondence.md): fixed-pressure nested energy, native field/force and documented variation-sign differences.
- [Hirshman and Whitson (1983)](https://doi.org/10.1063/1.864116): variational equilibrium genealogy. Possessing/citing this reference does not establish every implemented discrete variation.
- [CHEASE paper (1996)](https://doi.org/10.1016/0010-4655(96)00046-X): Hermite/Gauss weak architecture; per-selector/source normalization remains pinned-code work.
- [DESC Part I (2023)](https://doi.org/10.1017/S0022377823000272): force collocation/basis; the implementation's least-squares stationarity and physical force zero are distinct.
