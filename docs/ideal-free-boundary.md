# 3D free boundary of the map lane (N4)

Status: first implementation on branch `gs-3d-free-boundary` (built on the N1
map lane and the L4 per-n BEM). Design: kin6d_3d_tc24 DESIGN §7.3, row N4;
[equilibrium-core.md](equilibrium-core.md) D3, D6, D8. Module
`kin6d_ideal_free_boundary` (`src/ideal_free_boundary.f90`), vacuum leaves
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
