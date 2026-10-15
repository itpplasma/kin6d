# 3D free boundary of the map lane (N4)

Status: first implementation on branch `gs-3d-free-boundary` (built on the N1
map lane and the L4 per-n BEM). Design: kin6d_3d_tc24 DESIGN §7.3, row N4;
[equilibrium-core.md](equilibrium-core.md) D3, D6, D8. Module
`kin6d_ideal_free_boundary` (`src/mhd/ideal_free_boundary.f90`), vacuum leaves
`gen/gen_vacuum_shell.f90`, driver `app/gs_free_3d.f90`.

## Formulation

One curved P2/P3 reference mesh of the region inside a fixed axisymmetric
coupling curve Γ_c, with the plasma boundary Γ_p as an internal chain
(`initialize_gs(..., interior_chain=)`).

- Plasma cells: N1's frozen-in energy W_p (unchanged, same generated
  leaves) in the position gauge. Every node, Γ_p included, carries its full
  poloidal position (q1, q2) per Fourier mode, and the toroidal displacement ω
  is zero. Relabelling along B needs a toroidal component, so this gauge has
  no relabelling null, and rigid motions are represented exactly. The Hessian
  is assembled in these coefficients (`plasma_hessian`); the gradient is N1's
  full-coefficient gradient.
- Shell cells (Γ_p to Γ_c): curl-free field B = ∇Φ,
  Φ = F_vac q3 + c θ/2π + Φ̃, with q3 the physical toroidal angle of the shell
  map, c the poloidal circulation (−μ0 I_p for current along +φ), θ the
  reference poloidal angle about the reference axis and Φ̃ single valued
  (P_k × real Fourier in φ̂, n ∈ nfp ℤ). The shell map is the harmonic
  extension of the Γ_p displacement (mesh motion only), imposed by
  multipliers.
- Outside Γ_c: per n the decaying exterior enters through the exact
  Steklov–Poincaré map S_n (FortFEM `tm_bem_assemble`,
  `tm_bem_steklov_symmetric`, curved P2 panels), and the coils (PF/CS at
  n = 0 from `gs_coils`, RMP at n ∈ nfp ℤ from `coils_3d`) through their traces
  on Γ_c. At n = 0 the multivalued part is carried by a unit-circulation
  filament inside Γ_c, so the exterior unknown is single valued.

Equilibrium is a stationary point of

L = W_p[X] − (1/μ0) E_v[X, Φ̃] + λ·C(X),
E_v = ∫_shell |∇Φ|²/2 + Σ_m c_m [⟨S_n u_m, u_m⟩/2 − ⟨h_m, Φ̃_m⟩_Γc],

with W_p the plasma energy at frozen fluxes, −E_v/μ0 the vacuum co-energy at
fixed coil currents and fixed c, u_m = Φ̃_m − a_m on Γ_c (a_m the external
trace, h_m the external normal field) and C the shell-map constraints. The Φ̃
equations are vacuum magnetostatics with B·n = 0 on Γ_p and continuity of B·n
and B_t on Γ_c. The Γ_p map equations give [[p + B²/2μ0]] = 0 as the natural
traction. The plasma position is not a [[B²]] relaxation: it follows from
the coupled Newton system with the coil traces fixed on Γ_c. Everything is
per field period (c_m = 2π/nfp at n = 0, π/nfp otherwise). The vacuum
density |cof(DX) g|²/(2μ0 J) is a generated leaf; FortSym proves its
reduction to R|∇Φ|²/2μ0 for the identity cylindrical map.

The KKT matrix is symmetric and indefinite: the vacuum block is negative
definite, and with fixed coil currents the vertical mode is a negative
direction of the reduced plasma Hessian. Newton therefore uses an LDLᵀ
(MUMPS) or LU (SuperLU) factor of the exact matrix, never an energy descent.
Globalisation:

- merit Σ g²/|H_ii| (primal rows) and Σ g²|H_xx| (constraint rows);
- Levenberg damping on the plasma map rows only;
- trust radii of 2 % of the Γ_c extent for displacement and 0.5 rad for ω.

Φ̃ and λ are first solved at fixed geometry. fortsparse reads CSC, so blocks
that are not symmetric are passed transposed.

Why not N1's gauge (transversal displacement plus ω, axis expansion).
Material can slide along surfaces there only through ω₀, whose ω₀ = f(χ)
part is a continuum null. Measured on the gs_free fixture (+15 mm R start):

- With N1's 1e-4 gauge term, Newton (and N1's own solver, even at fixed
  boundary) drifts along that null and drags the surfaces 1–5 cm.
- With ω₀ pinned, the map is over-constrained: it stalls 7 mm short, although
  the exact target is a discrete equilibrium (started there it stays within
  16 µm).
- With ω₀ free, Newton crawls (−1.6 mm after 15 iterations).

The position gauge returns in about 3 Newton steps.

Vertical mode. `ideal_fb_rigid_stiffness` reports the generalised force and
stiffness of rigid R and Z shifts, with shell map, Φ̃ and λ re-solved. On the
fixture: K_R = +1.5e7 and K_Z = −2.7e4 J/m² per field period (decay index
−0.44). The vertical mode is unstable with fixed currents and about 600×
softer than the radial one. Newton converges to the saddle. Near the roundoff
floor of the merit the vertical position creeps by about 3e-5 m per step, and
the remaining vertical force over K_Z (0.1–0.3 mm) is the vertical
uncertainty. Stopping rule: relative merit, absolute residual in σV units, or
a stagnating merit with sub-step_tolerance increments.

## Inputs and the offset start

`gs_free_3d --config in.nml` (`&gs_free_3d`) reads a gs_free result (`mesh.dat`,
`profile-law.dat`, `summary.dat`, coil table, Γ_c curve). The frozen data are
the 2D state pushed forward by an offset map T: a translation (offset_r,
offset_z) or an elongation about the axis. So X = T⁻¹ is the exact target and
the seed X = id is displaced from it. Profiles are clamped to the plasma side
and p is measured from p_edge, since the map lane sees absolute pressure at
Γ_p and GS sees only p′. Outputs:

- `summary-3d.dat`: Newton report, boundary distance to the 2D boundary
  (ψ/|∇ψ| at every φ plane), interior surface defect, Γ_p normal amplitude per
  n, timings;
- `boundary-3d.dat`.

## Verification

- `ideal_fb_derivatives`: gradient of L against central differences and the
  KKT Hessian against gradient differences, all blocks, N = 0 and N = 1
  (nfp = 3); symmetry to roundoff.
- `ideal_fb_vacuum_exact_n0_p{2,3}`: exact n = 0 oracle (PF filaments plus a
  current filament, Γ_p one of its flux surfaces; the shell solution must be
  that field). Γ_c with 48/96/192 edges: P2 3.0e-3, 1.4e-3, 3.8e-4; P3 3.8e-4,
  9.9e-5, 1.6e-5 (relative L2 of B). The orders (1.9 and 2.6 on the finest
  pair) are still rising: the filament field varies on its distance to Γ_p.
- `ideal_fb_vacuum_gamma_c_independence_n3`: n = 3 field of a 3-fold window
  coil set at shell probes, Γ_c scaled by 1 and 1.12: 2.7e-4 relative (P3),
  1.0e-2 (P2).
- `gs_free_3d_nl3a_fixture_return`: NL-3a on the gs_free fixture (see below).

## NL-3a results (2026-10-13, sCluster, MUMPS, unregistered dev runs in
`~/N4work/runs`)

Distances are the gauge-invariant boundary distance of the mapped Γ_p nodes to
the 2D boundary (max / rms) and the interior surface defect. "Return"
compares each offset start with the unshifted start r0.

| case | Γ_p distance max / rms (mm) | surface defect (mm) |
|---|---|---|
| fixture P3 (limiter, 340 kA), r0 | 0.007 / 0.005 | 0.007 |
| fixture ±15 mm R, ±15 mm Z | ≤ 0.020 / 0.014 | ≤ 0.037 |
| fixture +2 % elongation (150 it) | 0.020 / 0.014 | 0.022 |
| ITER FB-cut 0.997 v2 (gs_free P3 L2, 15.52 MA), r0 | 0.94 / 0.085 | 0.61 |
| ITER ±15 mm R, ±15 mm Z | 0.95–1.17 / 0.085–0.094 | 0.42–0.60 |
| ITER +2 % elongation (150 it) | 8.4 / 1.1 | 11.9 (not converged) |

ITER runs: 219k unknowns, about 100 s for 40 iterations on 16 threads.

- The ITER r0 distance is the map-lane/ψ-lane discretisation difference,
  largest near the X-point side of the 0.997 cut. All four R/Z offsets return
  to the r0 state within about 0.2 mm.
- Rigid stiffness on ITER: K_R = +2.5e8 and K_Z = −5.4e6 J/m² per period,
  vertically unstable with fixed currents.
- The circulation from the 2D field matches −μ0 I_p to 1e-6, and the vacuum
  shell field matches the 2D field to 7e-4 (relative L2).

Mode stability (`stability=.true.`, shifted inverse iteration on the Schur
complement, metric |diag H|) at the ITER r0 state:

- n = 0: negative eigenvalues down to −3.9e-7 (the vertical mode);
- n = 3: none negative for shifts 0 … −0.1, smallest +3.7e-5.

DCON reports this FB-cut state n = 3 unstable without a wall (−0.097). The
discrepancy is open: the probe may lack resolution of the edge-current drive
(128 chain nodes, P3, one shell layer), or the two differ in model.

RMP stage. The ELM coils (bl2, R ≈ 7.7–8.7 m) lie inside the P12 coupling
ellipse, which violates the formulation (coils must be outside Γ_c). The first
scan used that ellipse and is invalid. `gamma_c_scale = 1.08` (the 2D boundary
scaled about the axis, 0.31 m coil clearance) is the RMP setting.

FB-taper(0.997, w = 0.0228), the adopted stable family. Reference: KIN6D 2D
`kin6d_fbtaper0997_w023n_fwdref_20261013_01/cases/fwd_p3_l3`.

- Resolution limits accuracy and convergence. With 128 chain and 192 Γ_c
  nodes, r0 sits 1.77 mm (max) from the 2D boundary and Newton creeps at
  2e-5 σV. With 256 and 384 nodes (863k unknowns, 8 min on 32 threads) it sits
  0.066 mm max / 6 µm rms, with surface defect 0.03 mm.
- At 128/192, the ±15 mm R and +15 mm Z starts return to the r0 state (≈ 1.9
  mm max, 0.17 mm rms). −15 mm Z and +2 % elongation are not converged.
- n = 3 stability probe: smallest eigenvalue +3.7e-5 (stable; DCON +0.206).
- The RMP ladder (0.25 … 0.005 × 90 kAt, Γ_c = 1.08 × LCFS, coil clearance
  0.31 m) at 128/192 is not converged (merit ratio 0.05–0.9 after 60
  iterations). Response/amplitude is 8.2e-3 m at 0.25 and a noisy
  2.9–5.1e-3 m below. The small-amplitude limit must be rerun at 256/384
  (about 2.6M unknowns at N = 1).

## Solver revision (2026-10-13, FB-taper at 256/384 nodes)

Where the iterations went, and the fixes:

1. **Tangential continuum null.** The position gauge still has one null family
   at n = 0: tangential flow on each surface weighted by 1/B_φ, which keeps
   the surfaces, the poloidal flux and the toroidal-flux 2-form. Newton
   directions of 2–7 cm along it were cut back to about 1 mm.
   - Fix: `tangential_damping` adds μ_t |H| t tᵀ to the plasma node blocks of
     the Newton matrix only (the gradient is untouched, so the equilibrium is
     unchanged).
   - μ_t = 1e-7 lies below the vertical eigenvalue (about 4e-7 in this metric);
     the ±15 mm Z starts still return.
   - Effect on r0: residual 1.0e-6 → 3.1e-9 σV in 40 iterations; Γ_p distance
     6.9e-5 m either way.
2. **Toroidal-field-scale mesh-motion rows.** The shell-map rows carry
   toroidal-field-scale forces on shell nodes that the multipliers absorb
   exactly. Their O(ε²) change after an RMP step outweighed the whole linear
   residual, so every full step was rejected.
   - Fix: Cᵀ is factorised once (it depends only on the reference); after every
     trial step λ is re-solved, those rows are zero, and the Γ_p plasma rows
     carry the corrected boundary force.
3. **One MUMPS handle with `sparse_refactor`** (fortsparse 7300414): the
   analysis is kept on 39 of 40 factorisations. Per iteration at 863k
   unknowns, N = 0: assembly 2.2 s, factor+solve 2.6 s; the line search is
   negligible.
4. **Restart and the ε → 0 limit.** `save_state`/`load_state` restart from a
   converged axisymmetric state. `linear_response` makes one exact Newton step
   restricted to n ≠ 0; at an axisymmetric state that is the exact ε → 0 limit
   of this discrete system.
5. **Γ_c for in-vessel coils.** `gamma_c_scale` and `gamma_c_modes` give a
   scaled 2D boundary with low-pass Fourier smoothing; the driver checks coil
   clearance and plasma margin. The production setting is 1.10 with M = 12:
   coil clearance 0.27 m, plasma margin 0.18 m. Pass `max_area` explicitly so
   the plasma mesh does not follow the Γ_c spacing.

FB-taper 0.997 results (fwd_p3_l3 reference, 256/384, μ_t = 1e-7):

- N = 0: r0 6.9e-5 m (max) / 5.9e-6 m (rms) from the 2D boundary.
- ±15 mm Z starts return to 7.1e-5 / 6.9e-5 m.
- RMP n = 3 ladder (bl2, phasing (0,0,0), restarted from the converged base):
  Γ_p normal amplitude per unit scale (90 kAt = 1) is 1.449e-2, 1.464e-2,
  1.466e-2, 1.472e-2, 1.468e-2 and 1.484e-2 m at 0.25, 0.1, 0.05, 0.02, 0.01
  and 0.005. The ε → 0 limit (linearised step) is 1.511e-2 m. Residuals are
  4–7e-8 σV after 40 iterations.
- The 0.25 point is 1.4 % below the small-amplitude values (finite-amplitude
  physics, reported separately). The remaining 2 % gap at the smallest scales
  sits at the N = 0 floor, where the vertical position still creeps by about
  2e-5 m per step.
- `rigid_vertical = .true.` applies dZ = -F_Z/K_Z (rigid n = 0 Z shift,
  K_Z from `ideal_fb_rigid_stiffness`, shell map, vacuum and multipliers
  re-solved) after every accepted step, kept only if the merit decreases. On
  the base state F_Z = -9 J/m and K_Z = -5.2e6 J/m^2 give dZ = -1.7e-6 m, which
  the merit rejects; the Gamma_p mean Z still creeps by +9.5e-6 m per step
  against that sign. The creep is therefore not the rigid vertical mode but a
  soft shape/gauge direction; the rigid correction is a no-op there.
- **Poloidal-slide gauge (`tangential_gauge`).** Over the 40-iteration ladder
  the Γ_p shape was fixed (normal change 0.26 µm rms) while the nodes slid
  along it (tangential rms 0.5 mm, mean 0.7 mm), and the n = 3 tangential
  content was still converging; the slow decay (0.6 % per step) is the
  λ/(λ + μ_t) relaxation of that null. With `tangential_gauge = 100` the n = 0
  tangential component of every plasma node is penalised in the Newton matrix
  and removed from the gradient of the merit. Result (rigid_vertical, restart
  from the converged base): s = 0.01 gives 1.51144e-2 m per unit scale against
  the linear limit 1.51140e-2 m (3e-5 relative, 54 iterations, residual
  9.7e-13 σV); s = 0.05 gives 1.50206e-2 m (-0.6 %, finite amplitude, 50
  iterations, 9.1e-13 σV). The old 2–4 % gap (1.4685e-2 at 40 iterations)
  was the gauge slide.
- The merit with the gauge projected out is nearly blind to the near-null
  vertical mode (a 15 mm Z start converged at iteration 1 with the base
  stopping rule). With `rigid_vertical` and the gauge the solve is therefore
  not converged while |F_Z/K_Z| > 1e-6 m, F_Z taken before the projection. The
  fixture-return test case `gzm` uses it.

## Relabelling (slide) gauge (2026-10-15)

Without `tangential_gauge`, the n = 0 relabelling null is now fixed exactly.
A map of the reference onto itself that keeps χ and the toroidal-flux 2-form
b3 dR dZ moves points with v = f(χ) ∇^⊥χ / b3 for any f. The fields, pressure
and surfaces do not change, so this is an exact continuum symmetry. It is
near-null once discretised.

- Diagnosis on the gs_free fixture (r0, SuperLU). The Newton direction was
  1.5 mm tangential and 2e-5 m normal. Its plasma block matched central
  differences of the gradient to 9e-4. Along it the merit followed (1 - α)²
  up to α = 1/32 and grew 1.9e4-fold at α = 1: the straight tangent line
  leaves the curved surfaces at second order. Adding the second-order path
  term (H a = -g'') gave merit 0.27 at α = 1/2. However, the iterate drifted
  along the null and the Newton target kept moving, so the solve crept:
  -0.3 % merit per step, step 1.1e-5 m, just above `step_tolerance`. The
  discrete minimiser along the null carries no physics.
- Gauge. The b3-weighted slide moments of the displacement from the start,
  M_l(u) = ∫ φ_l(ρ) b3 v · u dR dZ, are held at zero. φ_l is the C¹ partition
  of unity in ρ = (1 - T)^1/2 with one layer per element layer (the nested
  ω gauge of the map lane). The moments enter the Newton matrix as bordered
  Lagrange rows, and the merit is taken with their multiplier removed.
  Integration by parts shows that M_l vanishes for every curl-free u
  (translations, elongation). Offset starts therefore need no relabelling,
  and a smooth f restores the moments for any continuum equilibrium.
  `ideal_fb_slide_gauge_curl_free` checks this: curl-free/slide moment ratio
  1.9e-3 → 5.4e-4 under refinement.
- Result (fixture, max_iterations 150). All six ungauged cases converge to
  residual ≤ 4.4e-11 σV: r0 in 2 iterations, rp 5, rm 4, zp 4, zm 4, el 18.
  They reach the same state: Γ_p distance 1.97–1.98e-5 m and surface defect
  2.02–2.04e-5 m. The previous r0 value of 2.7e-6 m was the start state
  before convergence.
- A first version that held the Euclidean nodal projection onto
  P_k(2T - 1) v at zero did not leave translations free. It returned rp to
  1.15e-4 m, which shows that the gauge metric matters.
- `tangential_gauge` is not a gauge. It removes every n = 0 tangential
  component of the plasma gradient. Where the surface tangent is R, these
  components include the net force on the plasma, so an offset start stops
  short of the equilibrium. With rp, the gauge converged (3e-13 σV,
  projected) with Γ_p moved only -0.18 mm instead of -15 mm, a net rigid R
  force of 1.4e3 J/m against ≤ 1.3e-4 J/m at the slide-gauged states, and
  a distance of 17 mm. The `gr0` and `gzm` cases still use it and pass: gr0
  starts at the equilibrium, and in gzm the joint vertical coordinate carries
  the net Z force explicitly.

## Coupled vertical Newton coordinate

With `rigid_vertical = .true.` and `tangential_gauge > 0`, the Newton
unknowns are the normal n = 0 plasma-node displacements, unchanged shell,
vacuum and nonaxisymmetric coordinates, and one global plasma Z shift.
Writing the physical increment as `dU = T dy + v dz`, with fixed reference
normals in T and unit Z entries on the plasma nodes in v, the equations are
`T^T g = 0` and `F_Z = v^T g = 0`. Their symmetric Jacobian is

```
[ T^T H T   T^T H v ]
[ v^T H T   v^T H v ]
```

Here H is the existing generated Lagrangian Hessian, including shell ties,
magnetostatics and multipliers. Unused tangential slots are identity rows with
zero right-hand sides; they do not enter the physical increment. This is an
exact normal-coordinate gauge, rather than a finite penalty. Levenberg and
tangential damping do not act on the global shift or its force row. The full
coupled increment goes through the existing Jacobian/trust/backtracking checks;
there is no separate post-step rigid correction in this mode.

The merit includes the unprojected net vertical force, weighted by
`cap_x**2 / energy_scale`, as well as the normal-projected residual. All stopping
paths require the last full Newton shift estimate to be below 1 micrometre,
and retain the pre-existing `|F_Z/K_Z| <= 1 micrometre` guard when the initial
rigid-stiffness diagnostic is available. That finite-difference diagnostic is
used only for stopping, never for the augmented equation or its Jacobian.
Backtracking cannot satisfy the shift guard simply by reducing the step length.
Without a tangential gauge, `rigid_vertical` retains its earlier behavior.
Restart files still contain physical map coefficients, not solver coordinates.

`ideal_fb_vertical_newton` checks the bordered Jacobian against finite
differences, its force row independently and symmetry, for N = 0 and N = 1
(nfp = 3). A retained constant manufactured load makes a discrete state with
vertically displaced Solov'ev frozen data an exact target. Both signs of a
4 mm shift return quadratically after shell/vacuum equilibration. This is a
manufactured discrete solver oracle, not an analytic coil-supported free
boundary or a continuous equilibrium certificate.
