# 3D free boundary of the nested map with circuit currents as parameters

Status: first layer, branch `lane/fb3d-nested-map-20261025`. Module
`kin6d_ideal_fb_circuits` (`src/mhd/ideal_fb_circuits.f90`) on top of the N4
residual in `kin6d_ideal_free_boundary` ([ideal-free-boundary.md](ideal-free-boundary.md)),
its formulation theory ([free-boundary-theory.md](free-boundary-theory.md)) and the
R2 coupling surface ([free-boundary-chart-surface.md](free-boundary-chart-surface.md)).
This note fixes the residual that the forward W7-X case and the coil-current
inverse share. There is one state and one residual R(x; p), with the same
Γ_c exterior, readback and certificates. The coil currents p enter only as
parameters of R, and the inverse uses the JVP/VJP in p of this R rather than a
separate solver.

## 1. State, residual and exterior

The state is the N4 KKT vector

x = (X, x_s, λ, Φ̃),

with X the plasma map coefficients in the position gauge, x_s the harmonic
shell map between Γ_p and the fixed coupling surface Γ_c, λ the shell-map
multipliers and Φ̃ the single-valued shell potential (real Fourier modes in
nfp·φ̂). The vacuum field in the shell is

B = ∇Φ,  Φ = F_vac q₃ + c θ/2π + Φ̃,

where F_vac is the toroidal period divided by 2π and c is the poloidal
circulation, equal to μ0 I_p. The residual is the gradient of the Lagrangian

L(x; a, h, F_vac) = W_p[X] − E_v[x; F_vac]/μ0 + λ·C(x),
R(x; p) = ∇_x L(x; a(p), h(p), F_vac(p)).

The exterior part of E_v couples the Γ_c trace u = Φ̃|Γc − a to the decaying
exterior through one of three interchangeable symmetric operators. All three
are existing code.

| exterior | operator | E_v exterior term | where |
|---|---|---|---|
| per-n (axisymmetric Γ_c) | Galerkin Steklov–Poincaré S_n per harmonic | Σ_m c_m [½⟨S_n u_m, u_m⟩ − ⟨h_m, Φ̃_m⟩] | N4, FortFEM toroidal-mode BEM |
| R2 (general 3D Γ_c) | full-surface DtN S on P2 panels | [½⟨S u, u⟩ − ⟨h³, Φ̃⟩ + harmonic periods]/nfp | `kin6d_equilibrium_boundary_3d`, FortFEM parametric BEM |
| frozen (#88) | any fixed SPD S of R2 shape | as R2 | `ideal_fb_attach_exterior_3d` |

Here c_m = 2π/nfp for n = 0 and π/nfp otherwise. The R2 full-turn pairing
is a³ = a and h³_m = nfp c_m h_m (`ideal_fb_full_turn_traces`; both the R2
adapter and the circuit layer use this single definition). On the Φ̃ rows at
the Γ_c trace nodes, the exterior contributes

−∂E_v/∂Φ̃|Γc = −(c_m/μ0)(S_n(Φ̃_m − a_m) − h_m)            (per n),
−∂E_v/∂Φ̃|Γc = −(S(Φ̃ − a³) − h³ + K κ)/(μ0 nfp)           (R2),

with K κ the fixed harmonic-period cross term of R2, which is independent of
p. Every other block of R is free of a and h.

## 2. Vacuum source linear in the circuit currents

Circuit k is a group of filaments c with winding counts N_c (a MAKEGRID
group with its EXTCUR, or a FreeGS-style circuit with multipliers). At unit
circuit current, its field is

b_k = Σ_{c ∈ k} N_c B_c[1 A],

where B_c is the closed-form straight-segment Biot–Savart field of
`kin6d_coils_3d`. The external field is

B_vac(p) = B_0 + Σ_k I_k b_k,  p = (I_1, …, I_K),

with B_0 fixed. On the present path B_0 is the unit-circulation filament
times c, which makes the trace single valued.

R depends on B_vac only through three functionals on Γ_c, all linear in
B_vac. They are those of `ideal_fb_set_drive`, now implemented once as
`ideal_fb_drive_traces` for any number of fields.

- Normal load: h_{i,m} = ∫_Γc (B_vac·n)_m N_i R ds.
- Single-valued trace: for n = 0, a = ∫ (B_vac − c ∇θ/2π)·t ds from trace node
  1, whose closure after one turn is the poloidal circulation of B_vac minus c.
  It vanishes for circuits outside Γ_c by Ampère's law. For n ≠ 0, a_n follows
  from B_φ = (1/R)∂Φ/∂φ, so that a = R B_φ,sin/n and −R B_φ,cos/n.
- Toroidal period: F_vac = R B_φ^(0) on Γ_c. For a field that is curl free
  across the shell, ∂(R B_φ^(0)) vanishes along the meridian, so F_vac is
  constant and equals μ0 I_link/2π, where I_link is the circuit current that
  links the torus. The layer reports the spread over the trace nodes.

Superposition then gives

a(p) = a_0 + Σ I_k a_k,  h(p) = h_0 + Σ I_k h_k,  F_vac(p) = F_0 + Σ I_k F_k.

`ideal_fb_circuits_build` evaluates all K circuits in one pass over the Γ_c
sample points. Per point, the coil basis takes an nphi-point toroidal DFT
over a full turn. The aliasing error is exp(−nphi d/R), with d the distance
from the coils to Γ_c. The result is (a_k, h_k, F_k) for every k.
`ideal_fb_set_currents` installs a(p), h(p) and F_vac(p) through
`ideal_fb_set_traces`. With an R2 exterior attached, it also updates the
full-turn drive and load.

## 3. Derivatives in p

Write e = (δa, δh, δF) = Σ_k δI_k (a_k, h_k, F_k). Then

J_p δp = ∂R/∂a δa + ∂R/∂h δh + ∂R/∂F δF.

The trace part is linear and touches only the Φ̃ rows on Γ_c:

δR_Φ̃(m, Γc) = (c_m/μ0)(S_n δa_m + δh_m)                     (per n),
δR_Φ̃(Γc) = (S δa³ + δh³)/(μ0 nfp)                             (R2/frozen).

F_vac enters only the shell density |cof(DX)·g|²/(2μ0 J) through
g = ∇̂Φ̃ + F_vac Q_{3,:} + c ∇̂θ/2π and through the chain row
dv(9+i, 6+3(i−1)) = F_vac of the shell map. With gd the generated leaf
gradient, hd its generated Hessian and v the local variables,

∂(gd·dv)/∂F = (hd_{:,g} Q_{3,:})·dv + gd·∂dv/∂F,  with ∂dv(9+i, 6+3(i−1))/∂F = 1.

`shell_point(..., fdot)` evaluates this from the existing generated
`vacuum_energy_hessian`, so no new handwritten pointwise derivative is
introduced. `ideal_fb_traces_jvp` assembles −(∂g_x/∂F, ∂g_Φ/∂F) on the x_s
and Φ̃ rows. Since g is linear in F, the density is quadratic in F, so R is
exactly quadratic in p and affine along circuits with F_k = 0.

The VJP is the transpose, ⟨v, J_p e⟩ = ⟨ā, δa⟩ + ⟨h̄, δh⟩ + F̄ δF, with

ā_m = (c_m/μ0) S_nᵀ v_Γc,m,  h̄_m = (c_m/μ0) v_Γc,m,  F̄ = −⟨v_x, ∂g_x/∂F⟩ − ⟨v_Φ, ∂g_Φ/∂F⟩,

or the R2 form of the same expressions. Contracting with the circuit columns
gives (J_pᵀ v)_k = ⟨ā, a_k⟩ + ⟨h̄, h_k⟩ + F̄ F_k (`ideal_fb_currents_vjp`,
OpenMP over k).

For the inverse (c), let O(x) be the observables (probe fields or boundary
points), with ∂O/∂x read back from the same state. At an equilibrium
R(x*; p) = 0 the implicit-function theorem gives

dO/dp = −(∂O/∂x) (∂R/∂x)⁻¹ J_p.

∂R/∂x is the symmetric KKT Hessian of L, which `ideal_fb_hessian` assembles
and the Newton solve factors as LDLᵀ/LU. The adjoint form solves
(∂R/∂x) μ = (∂O/∂x)ᵀ with the same factor and returns
dO/dp = −J_pᵀ μ = −`ideal_fb_currents_vjp`(μ). The cost is one solve per
observable block instead of one per circuit. On the slide-gauged system the
bordered matrix of `ideal_fb_slide_system` takes the place of ∂R/∂x.

## 4. Readback and certificates

The layer does not change the definition of R. It changes only how
(a, h, F_vac) are produced, so every existing readback and certificate
applies unchanged at a given p. These include the residual norm in σV units,
the solved-vacuum field and traction defects, the hypotheses (V), (R) and
In(S_red) with slide bordering, the harmonic-pair stability audit and the
frozen-exterior energy normalisation. The layer adds three per-circuit
diagnostics: the n = 0 poloidal closure (zero for unlinked circuits), the
n ≠ 0 curl-free indicator of the P2 trace, which is a discretisation measure
and not roundoff, and the spread of R B_φ^(0) together with an Ampère check
of F_k.

## 5. Reused pieces and new code

| piece | status |
|---|---|
| N4 residual, Hessian, Newton/KKT, slide gauge (`kin6d_ideal_free_boundary`) | reused; `ideal_fb_set_drive` now wraps the batched `ideal_fb_drive_traces` |
| per-n exterior S_n (FortFEM toroidal-mode BEM) | reused |
| R2 general 3D Γ_c (`kin6d_equilibrium_boundary_3d`, FortFEM parametric BEM) and its adapter | reused; full-turn mapping factored into `ideal_fb_full_turn_traces` |
| frozen exterior (#88) | reused; traces now replaceable after attachment (`ideal_fb_set_traces`) |
| Biot–Savart segments (`kin6d_coils_3d`) | reused |
| field-basis contract `ideal_fb_field_basis_t`, coil circuit basis, circuit traces, currents → traces, JVP/VJP | new |
| exact ∂/∂F_vac of the shell gradient (existing generated Hessian leaf) | new |

Generic pieces stay in their owning libraries. Linear algebra goes to fortnum,
BEM to FortFEM and coil file readers to libneo. libneo origin/main 7db0217 has
`neo_biotsavart::load_coils_from_file`, which reads a point list with one
current column (x y z I, no group or winding columns), and a Python
STELLOPT/MAKEGRID parser `python/libneo/coils.py`. It has no Fortran MAKEGRID
group reader. Any such reader belongs in libneo and should output the neutral
polyline table of `coils_3d` with circuit and turns columns.

## 6. Parallelism and generated kernels

- The coil basis samples all circuits at once with OpenMP over sample points.
  Each point is independent, so the result does not depend on the thread
  count. Circuits share the per-point segment sweep.
- The ∂/∂F tangent runs inside the existing OpenMP cell loop of the shell
  energy. The circuit VJP contraction runs with OpenMP over circuits.
- Strong scaling of the trace build has not been measured yet. The present
  test cases are small, so performance claims are open.
- The segment Biot–Savart kernel is the existing handwritten closed form in
  `coils_3d`. Emitting it with FortSym, including its derivative with respect
  to the filament points for later coil-shape parameters, is a follow-up. The
  F tangent uses only existing generated leaves.

## 7. Ideas from other codes (re-derived, no code copied)

- VMEC/NESTOR and VMEC++ `mgrid` treat the vacuum field as a precomputed
  per-group field that is scaled by EXTCUR and summed (Hirshman, van Rij and
  Merkel, Comput. Phys. Commun. 43 (1986) 143; VMEC++ `mgrid_provider.cc`
  AccumulateCircuit). KIN6D uses the same linearity in a different form. It
  precomputes the per-circuit Γ_c traces (a_k, h_k, F_k) rather than a
  volume grid, so no grid interpolation enters the residual.
- FreeGS circuits are groups of coils that share one current with per-coil
  multipliers (FreeGS `Circuit`). Here this becomes the winding count N_c in
  b_k.
- DESC treats coil currents as differentiable parameters of the external
  field in free-boundary optimisation (DESC repository, PlasmaControl/DESC).
  This is a hypothesis about its exact interface. KIN6D derives the
  parameter derivative J_p analytically (section 3) instead of using
  automatic differentiation, and it reuses the symmetric KKT factor for the
  adjoint.
- GVEC has no free-boundary path that is used here.

## 8. Verification (synthetic public coils only)

| CTest | checks |
|---|---|
| `ideal_fb_landreman_axisymmetric_circuits` | Landreman fitted drive written as 24 ring circuits plus a TF circuit (unit 1/R), currents = fit coefficients and F(δ). The traces equal the sampler traces to roundoff (7.7e-15 relative), and the unchanged gates and refinement checks of `ideal_fb_landreman_axisymmetric` pass. |
| `ideal_fb_circuits_derivatives` | nfp = 3, N = 1, four circuits (PF pair, 2-turn inner ring, 3 window coils, 12 planar TF coils). The JVP matches central differences of R to roundoff (≤ 2e-12 relative; exact because R is quadratic in p). The adjoint identity holds to 3e-15. Trace superposition holds to 6e-16, the residual second difference along F_k = 0 circuits is 3e-15, and the third difference along the TF current is 2e-15. F_TF matches Ampère μ0·12/2π to 8e-14. |
| `ideal_fb_circuits_exterior3d` | the same checks with a frozen full-surface R2 operator attached |
| `ideal_fb_landreman_helical_circuits` | nfp = 3, N = 1, Landreman drive plus six alternating-sign window coils (n = 3(2j+1)) with \|B_3\| = 1e-4 on Γ_c. Newton converges in 5 steps to an absolute residual of 1.4e-12, the axisymmetric budgets still pass, and the n = 3 boundary amplitude is linear in the current (ratio 2.00000 for I/(I/2)). |

The last row is numerical evidence for a small nonaxisymmetric response, not
a validated 3D equilibrium. The Landreman family has closed field lines, and
the per-n exterior assumes an axisymmetric Γ_c. Same-sign windows with
\|B\| = 1e-3 add an n = 0 radial field. With them, Newton first shifts the plasma
by 9 mm, about 0.09 a, and then stalls at merit 3.7e-2. This is recorded as an
observation. That the vertical mode under fixed currents causes the stall is a
hypothesis that has not been traced.

## 9. Gaps before the W7-X forward case (b) and the inverse (c)

1. General 3D Γ_c traces (done, step b start). `ideal_fb_drive_traces_3d`
   evaluates the same three functionals on Γ_c = H(γ × [0, 2π)) of the plasma
   chart H, the geometry of `ideal_fb_surface`, at Cartesian points
   X = H(γ(t), φ) with x_t = dH(γ′, 0) and x_φ = dH e_3:
   nfp c_m h_{i,m} = ∫ B·(x_φ × x_t) N_i trig_m dt dφ, the n = 0 drive from
   ⟨B·x_t⟩ along γ, the n ≠ 0 drive from B·x_φ = ∂Φ/∂φ at the trace nodes, and
   F_vac = ⟨B·x_φ⟩. ⟨·⟩ and the harmonics use the nphi-point trapezoid rule
   over the full turn. The field basis has a Cartesian contract
   (`sample_points`); coils evaluate Biot–Savart directly (OpenMP over
   points), other bases fall back to their harmonics.
   `ideal_fb_circuits_build(..., nphi_surface=)` and `ideal_fb_set_drive_3d`
   use it. With the cylindrical chart the result equals the meridional path
   to roundoff (test `ideal_fb_traces3d_axisymmetric`).
2. R2 cost (field-period reduction done). `ideal_fb_use_bem3d(...,
   periodic=.true.)` assembles FortFEM's block-circulant operators from the
   target pairs of one field period only
   (`assemble_parametric_bem_p2_circulant`). The traces repeat in every
   period, so only the zero sector harmonic acts:
   R^T S R = nfp R_0^T S_0 R_0 with S_0 = W_0 + J_0^T V_0^{-1} J_0,
   J_0 = M_0/2 - K_0 and X_0 the sum of the circulant blocks
   (`gamma_c_3d_build_periodic`). This is exact algebra (test
   `ideal_fb_exterior3d_periodic`: roundoff against the full-surface DtN of the
   expanded operators and against the dense all-pairs build). The pair work
   falls by nfp and the dense DtN by nfp^3. The pair quadrature itself
   (FortFEM `converge_pair`) remains the cost driver.
3. Curve-chart free boundary (done). `ideal_fb_init` takes any
   `ideal_chart_t`. The per-n exterior is built only for the cylindrical chart
   itself; any other chart (curve chart, nested-map chart) gets no exterior
   until the R2 one is attached, and the residual, Newton step and traces
   return bad input until then. The test-only cylindrical subclass with an
   overridden jet is gone; `ideal_fb_traces3d_helical` uses the helical-axis
   curve chart.
   On the synthetic nfp = 3 case (`fb3d_synthetic`) the frozen iota needs the
   sign of the coil field's rotation in the chart frame. A multi-turn trace of
   the vacuum field in chart coordinates
   (`benchmarks/fb3d-synthetic/iota_sign.py`) turns clockwise in (q1, q2) at
   0.190-0.191 per unit phi-hat, while VMEC++ reports iota = +0.190-0.201
   with theta counter-clockwise in (R, Z). The field therefore takes
   sign = -1, the default of `fb3d_synthetic`. Newton from the disc start is
   not yet robust. Its steps are dominated by the motion of single plasma
   nodes (a Gamma_p node, sometimes the axis node), along which the merit
   grows like alpha^4, until the Levenberg ladder runs out. These directions
   are not a gauge freedom. At the converged N = 2 state the reduced Hessian
   S_red of `ideal_fb_hypotheses` (slide gauge kept) is positive definite,
   and its smallest eigenvalues (about 21, 22, 91 and 97 J/m^2 per field
   period at 12 edges and chain 16) belong to smooth n = 0 motions: the
   rigid Z and R shifts and a shift of the core about the axis. No localised
   axis or boundary-tangential mode is soft there, so fixing those nodes
   would constrain the physical labelling and change the discrete
   equilibrium. The stalled states lie about 2 cm (rms over Gamma_p) from
   the solution, and their S_red has three or four negative eigenvalues. The
   one or two largest in magnitude, -2e5 to -4e5, are concentrated on a
   single Gamma_p node and spread over all Fourier modes. In the stall where
   the curvature was split, the plasma block alone is positive along them,
   so the negative curvature comes from the eliminated shell, vacuum and
   exterior block. A stall is therefore a nonzero local minimum of the
   residual merit, which is a globalisation problem and not a missing gauge.
   Converged runs at 1, 10, 15 and 16 threads agree to 1e-6 mm on Gamma_p,
   and the image of the axis node lies at R = 1.1410 m, Z = 0 at
   phi-hat = 0 (`axis_RZ_phi0` in the JSON), 0.6 mm from the traced vacuum
   axis (R = 1.1404 m) and 2.2 cm outboard of the VMEC++ axis (1.119 m). The
   R2 exterior is assembled with atomic accumulation of the pair blocks
   (FortFEM `assemble_parametric_bem_p2_circulant`), so its roundoff differs
   from run to run even at a fixed thread count, and whether a run converges
   or stalls depends on that roundoff.
   With `ideal_fb_t%negative_curvature` (`fb3d_synthetic` globalisation 2,
   opt-in) a stalled merit hands over to `ideal_fb_negative_curvature`,
   which takes exact More-Sorensen trust-region steps on the second-order
   model of the reduced energy E_red(q) (shell map, Phi~ and lambda solved
   exactly, slide gauge) while S_red has a negative eigenvalue, and a state
   that passes the residual gate with negative curvature is left in the same
   way. In sCluster job 3370099 all 15 glob-2 runs at R1 converged (1 to 16
   threads, starts kicked by 8 mm along q1 included) and 2 of 3 at R2. The
   R2 failure parked the map at a fold, which the Jacobian guard of the
   escape now prevents, and with it job 3370145 converged in 8 of 8 glob-2
   runs (R1, R2 and R3, kicks up to 15 mm). All R1 runs reached the same state within
   4e-7 mm, in 34 to 87 Newton steps against 93 to 152 for the trust region
   alone.
   The case is stellarator symmetric: the coils (make_coils.py), the chart
   frame, which maps (q1, q2, phi-hat) to (-q1, q2, -phi-hat) under
   (x, y, z) -> (x, -y, -z), the frozen labels and, at R1 only, the reference
   mesh (49 of 49 nodes mirrored; 40 of 77 at R2 and 48 of 97 at R3). The
   gate also admits symmetry-broken equilibria. A start kicked by 8 mm
   without the escape converged (residual 2e-13) to a state with boundary
   symmetry defect 11.6 mm, one negative eigenvalue of S_red and an E_red
   0.067 J above the symmetric state, which is a saddle of the reduced
   energy, and the extra state of job 3369792 at 8 threads has the same
   defect and the same 12.98 mm distance from the symmetric state. Starts
   kicked by 8 mm stall under the trust region alone at R1 and R2, at
   states 24 to 31 mm from symmetry with 5 to 51 negative eigenvalues. The
   acceptance check is therefore In_-(S_red) = 0, the boundary symmetry
   defect, the vacuum oracle and, where available, the distance to a
   reference (`fb3d_synthetic` "state" block and
   `benchmarks/fb3d-synthetic/state_check.py`). The stable state has a
   boundary symmetry defect below 1e-7 mm at R1 and 0.9 mm at R2, where the
   mesh is not mirror symmetric. The single R3 run converged to a stable
   state with a 5.4 mm defect, and whether that is the mesh asymmetry or a
   second stable branch is open.
4. Coil input. The MAKEGRID group reader is `neo_makegrid_coils` in
   [libneo](https://github.com/itpplasma/libneo) (PR #469). KIN6D compiles
   its two dependency-free modules from the pinned `KIN6D_LIBNEO_REVISION`
   (CMakeLists.txt) and `makegrid_circuits_read` turns a file into coils_3d
   polygons with circuit = MAKEGRID group and turns = fourth column (per unit
   EXTCUR). The W7-X sign and EXTCUR convention (A·turn versus A) is still
   settled by T1 of the scoping note, and W7-X data stay in the private data
   repository.
5. Accuracy at W7-X distances. The aliasing error exp(−nphi d/R) and the
   Γ_c clearance need checks against VMEC++ at matched accuracy. Field-period
   symmetric coil evaluation would cut the sampling cost by nfp.
6. Inverse API. This needs observable operators O(x) with ∂O/∂x and an
   adjoint solve through the existing KKT factor, followed by the
   Gauss–Newton or GSVD loop in fortnum on dO/dp.
