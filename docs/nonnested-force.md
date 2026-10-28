# Projected force, harmonic deflation and the residual floor

Status: implemented (2026-10-28) in `src/mhd/nonnested_force.f90`
(`kin6d_nonnested_force`) and `equilibrium_residual_floor` in
`src/equilibrium/equilibrium_harness.f90`; oracles in
`test/test_nonnested_force.f90`. Items 3 and 5 of
[relaxation comparators](relaxation-comparators.md#algorithms-to-adopt).
All numbers are candidate diagnostics; nothing here is an enclosure.

## Force

The force of a state on the NN1 real-Fourier complex is the projection
P(J x B) of J x B onto the velocities u in V2_0 (RT x DG with zero normal
trace, all Fourier blocks) with div u = 0. P removes every gradient, so
P(J x B - grad p) = P(J x B): the pressure is the multiplier and need not be
supplied. It is the force measure of the MRX relaxation (arXiv:2510.26986,
Leray projection with an incompressible multiplier pressure) applied to a
kin6d state. Two discrete forms of the same projection are provided:

- Leray (mixed): [M2 D^T 0; D 0 l; 0 l^T 0] (u, phi, mu) = (F, 0, 0), with
  F = <V2, J x B>, D the metric-free coefficient divergence V2 -> V3 and l
  the n = 0 integral functional that removes the constant from ker D^T on
  V2_0. kin6d had no Leray projection in any lane; this is the comparison
  path.
- Potential velocity: u = C psi + sum_j <F, h_j> h_j with K psi = C^T F on
  V1_0, K = C^T M2 C, i.e. curl (curl curl)^-1 curl (J x B) + <J x B, h> h.
  K is singular (gradients, and in a shell the relative harmonic 1-form);
  the solve is deflated by the gauge multiplier lambda in V0_0,
  [K M1 G; G^T M1 0] (psi, lambda) = (C^T F, 0), so psi is M1-orthogonal to
  ker K and a right-hand side with kernel content returns the same psi
  (its kernel part is reported).

On the exact complex ker D in V2_0 = C V1_0 + H, so both forms are the same
M2-orthogonal projection: they agree to solver roundoff, and D u = 0 holds as
a matrix identity for the potential form.

Harmonic 2-forms are computed per geometry (`nf_prepare`): the toroidal
class (an n = 0 DG field with toroidal flux) always, and in a shell
(`inner` wall nodes) also the poloidal class C a_in with a = R A_phi = 1 on
the inner wall. Each candidate g becomes h = g - C psi(g) (the discrete
Hodge projection) and is M2-orthonormalised. The shell's relative harmonic
1-form G phi_in is deflated by tying the n = 0 multipliers of the inner wall
to one unknown, the NN2 floating-conductor tie. `nf_deflate` removes the
harmonic part of any V2 vector; `nf_potential_solve` is the deflated
singular curl-curl solve.

Current. `nf_lorentz_load` has two discrete currents. The default (dual)
J = M1^-1 (C^T M2 B + f_vac C^T <V2, e_phi/R>) in V1_0 is the FEEC curl* used
by the metriplectic relaxation. Its zero tangential trace forces J to zero
at the wall, so where n x B /= 0 there it leaves an O(1) layer of width h and
P(J x B) decreases only like h^1/2 (measured below). `primal = .true.` uses
J = C B1 in V2 with B1 the full-V1 L2 projection of B (no boundary
condition), which has no such layer and converges at order k - 1.

## Verification (this VM, 2 threads, `ctest -R nonnested_force`)

Solov'ev domain, P2 rings 8/16/32 and P3 rings 8/16, N = 1, manufactured
f = grad phi + curl w + 12 e_phi/R with w = p(R, Z)(0, 1 + 0.3 cos phi,
0.4 sin phi), p = 0 on the wall, so P f = curl w + 12 e_phi/R:

| quantity | P2 | P3 |
|---|---|---|
| L2 error of u / \|P f\| | 3.44e-3, 1.03e-3, 2.83e-4 | 3.90e-4, 8.63e-5 |
| observed orders | 1.75, 1.86 | 2.18 |
| \|u_potential - u_Leray\| / \|u\| | <= 3.4e-13 | 2.4e-13 |
| max \|D u\| / (\|D\| max \|u\|), potential / Leray | <= 1.7e-15 / 2e-16 | 7e-16 / 7e-17 |
| harmonic orthogonality \|<h, C psi>\| / \|C psi\| | <= 7e-29 | <= 7e-28 |
| \|P grad phi\| (gradient load) | <= 7e-14 | <= 1.3e-13 |

The orders are pre-asymptotic on these flux meshes (without the harmonic
term 1.82/1.92 and 2.54); the tests gate at k - 1/2 (P2) and k - 1 (P3).
Shell (NN1 coil annulus, P2, 2/4/8 rings): two harmonic 2-forms, errors
6.79e-2, 1.75e-2, 4.42e-3 (orders 1.95, 1.99); a right-hand side with
M1 G (phi_in + interior noise) added (1.2e-3 of its norm) returns the same
psi to 1.4e-12, and that kernel part is reported (1.25e-3) while a
consistent right-hand side has 2.5e-16.

Factorisation cost (P2, 8/16/32 rings, N = 1): potential 0.11/0.47/2.3 s,
Leray 0.17/2.0/47 s (native supernodal LU, AMD; no METIS on this VM). The
P3 32-ring potential system with three Fourier blocks needs more than 6 GB
and is not run here. The Leray path is therefore factored only on the
cheaper levels of the tests (`leray = .false.` skips it).

Discrete equilibria, |P(J x B)|_M2 / <|J|^2 |B|^2>^1/2:

| case | dual current | primal current |
|---|---|---|
| Solov'ev NN1 P2, 8/16/32 | 8.83e-4, 6.46e-4, 1.53e-4 (0.45, 2.08) floor | 1.09e-3, 5.26e-4, 2.45e-4 (1.05, 1.10) |
| Solov'ev NN1 P3, 8/16/32 | 6.83e-4, 4.94e-4, 1.83e-4 (0.47, 1.43) floor | 3.51e-5, 1.09e-5, 2.92e-6 (1.69, 1.90) |
| Lundquist Beltrami, cylinder P2, 4/8/16 | 1.78e-2, 1.33e-2, 9.59e-3 (0.43, 0.47) floor | 2.06e-2, 1.07e-2, 5.26e-3 (0.94, 1.03) |
| vacuum a = R^2 with F/R, P2 8 rings | 1.1e-13 (relative to B^2/L) | |

Both force forms agree to <= 1.1e-12 on every state; the harmonic part of
J x B is <= 1e-12 (no net toroidal force). The registered tests run the P3
Solov'ev ladder on 8/16 rings (the 32-ring P3 level took 8 min here) and the
Beltrami ladder on 4/8 rings (`test_nonnested_force beltrami 3` adds 16).

## Residual floor

`equilibrium_residual_floor(residuals, ratio, norm, report)` records the
observed orders log(r_(i-1)/r_i)/log(ratio) of a refinement ladder and
reports `decreasing`, or `floor` with the finest residual when some
refinement reduces it by less than ratio^-min_order (min_order 1/2 by
default); `write_residual_floor` prints it in the D8 `[V]` style. Gates
report the floor instead of chasing it. On the straight periodic cylinder
(P2, 4/8/16/32 rings, primal current) a smooth B_z(r) (tanh, width 0.1)
decreases with orders 2.19, 2.42, 2.48, while a tangential discontinuity of
B_z at r = 0.47 (off the mesh rings; J x B = -grad B_z^2/2 in the continuum,
a current sheet balanced by a pressure jump) gives 3.97e-2, 1.49e-2,
1.66e-2, 2.65e-2: a floor. The dual current's h^1/2 wall layer is also
flagged (table above): it is a boundary artefact of J in V1_0, not physics,
and the primal current removes it. The ideal-map ladders of
`test_ideal_map_axisymmetric` report their D3 residual through the same
routine; the A3 P2 seed ladder (32/64/128 nodes) gives 3.24e-2, 1.29e-2,
3.46e-3 (orders 1.32, 1.90, decreasing).

## Open

- The Leray saddle LU dominates at 32 rings; a better ordering (METIS) or
  an iterative Schur solve would make it a practical production path. The
  potential form is the intended one.
- FortFEM exposes no Fourier-block divergence; `divergence_blocks` applies
  FortFEM's documented block rule to `fs%div` and is a candidate for
  upstream.
- The force is a diagnostic; no lane uses it as a Newton residual yet.
