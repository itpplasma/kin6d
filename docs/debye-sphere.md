# Finite spherical source and Debye equilibrium

A finite charge source and a plasma response change the potential in different
ways. The source radius makes the bare potential finite at the center. The
linear electron response then suppresses the exterior Coulomb tail. The
`debye_sphere_demo` teaching example separates these effects using a uniformly
charged, permeable three-dimensional sphere and its matched Debye equilibrium.

This is an **analytical reference evaluation**. Its animation sweeps an
equilibrium parameter; it does not compute a time-dependent shielding process
or an emergent particle simulation. The separate slab shielding example owns
the computed transient. The dots and plus signs in the sphere cross-section
illustrate mobile electrons and a fixed continuous source, respectively.
They do not generate the radial field.

## Model and normalization

The source has positive total charge Q, radius R, and uniform volume density
3Q/(4πR³). Electrons occupy the source interior and the exterior. Ions supply a
uniform fixed background. The radial potential is regular at the origin,
continuous with its radial derivative at R, and zero at infinity. There is
no surface charge at the source boundary.

Set s=r/R, κ=R/λD, and Φ=φ/[Q/(4πε₀R)]. The bare potential is

\[
\Phi_C(s)=\frac{3-s^2}{2}\quad(s\leq1),\qquad
\Phi_C(s)=\frac1s\quad(s\geq1).
\]

The screened potential solves

\[
\Delta_s\Phi_D-\kappa^2\Phi_D=-3\quad(s<1),\qquad
\Delta_s\Phi_D-\kappa^2\Phi_D=0\quad(s>1).
\]

Matching the regular interior and decaying exterior solutions gives

\[
\Phi_D(s)=\frac3{\kappa^2}
\left[1-(1+\kappa)e^{-\kappa}
\frac{\sinh(\kappa s)}{\kappa s}\right]\quad(s\leq1),
\]

\[
\Phi_D(s)=\frac3{2\kappa^3s}
\left[(\kappa-1)e^{-\kappa(s-1)}
+(\kappa+1)e^{-\kappa(s+1)}\right]\quad(s\geq1).
\]

The interior quotient takes its regular value 1 at s=0. The executable admits
s∈[0,6] and κ∈[0.5,2]. The separate bare curve represents the unscreened
limit; κ=0 is rejected by the screened evaluator.

The electron response view reports nₑ/n₀=1+χΦD with χ=0.02. Here
χ=eQ/(4πε₀RkBTₑ), so the maximum bare thermal potential is 3χ/2=0.03.
This states the small-source ordering needed for linearization. The κ sweep
can be read as changing the background density at fixed source and electron
temperature. Nonlinear Boltzmann response, finite boundaries, kinetic dynamics,
physical-data uncertainty, and model-reduction errors are outside this reference.

## Evaluation and continuous bands

[The retained generator](../gen/gen_debye_sphere.f90) is the mathematical
source for floating-point primal/JVP/VJP and outward interval leaves. The
interior uses the regular polynomial

\[
S_{12}(z)=\sum_{j=0}^{12}\frac{z^{2j}}{(2j+1)!},\qquad z=\kappa s,
\]

with an explicit remainder. For z≥0 and derivative order m=0,…,3,

\[
|[S-S_{12}]^{(m)}(z)|\leq
T_m(z)=\frac{e^z z^{26-m}}{(27-m)!}.
\]

Termwise differentiation gives coefficient
(2j)ₘ/(2j+1)!≤1/(2j+1−m)!. Extending the positive parity series to all
powers and comparing factorials with the exponential series gives the stated
bound, including the origin by continuity.

Writing A=3(1+κ)e⁻κ/κ² and B=|∂κA|, the omitted interior terms contribute
at most A T₀ to ΦD and B T₀+A s T₁ to ∂κΦD. The curvature margins use
Aκ²T₂ and Bκ²T₂+A(2κT₂+κ²sT₃), respectively. The generated polynomial
curvatures and these remainder terms bound the exact reference curvature.
The exterior formula is evaluated directly with outward interval arithmetic.

[The reference wrapper](../src/debye_sphere_reference.f90) validates the
admitted finite domain and splits radial intervals crossing s=1. Each plotted
node first receives its strict reference interval. Every radial cell then
adds the chord margin M h²/8, with M bounding the second radial derivative;
each vertex is widened by the maximum adjacent margin. Linear interpolation
of the resulting bounds encloses the reference throughout every cell.
The exact potential and κ tangent are C1 at s=1, with bounded piecewise
second derivatives, so the same chord theorem covers cells crossing R.

The bands cover analytical evaluation, series truncation, arithmetic, and
radial display reconstruction. They do not certify a numerical Poisson/BVP
solve. κ derivatives have interval bounds; emitted JVP/VJP use a Euclidean
scalar pairing on each approximation branch. Source-radius, charge, uncertain
input, and transient sensitivities are not reported. No formal-proof status
is claimed.

## Running and outputs

Build with the project's CPU preset, then run:

```sh
build/debye_sphere_demo
build/debye_sphere_demo --snapshot 97
```

The default sweeps 193 equilibria from κ=0.5 to 2 at 24 frames per second
(8 s) and holds the last one for 1 s. `--snapshot [frame]` evaluates and
renders only the requested frame, with 97 as the default. One 2x2 view uses
fixed axes: the screened potential on the z=0 plane (equal aspect, with r=R
and r=λ_D marked), radial bare and screened potentials, the electron density
response and the certified band half-widths, which are below a pixel in the
curve panels. Before 2026-10-06 the sweep had 121 equilibria and the density
was a separate `debye-sphere-response` movie.

| Output | Content |
|---|---|
| `output/debye-potential-reduction.mp4/.png` | Potential map, radial potentials, electron density, band half-widths |
| `output/debye-potential-reduction.csv` | Plotted nodes and full potential/density/κ-tangent bands |
| `output/debye-sphere-evidence.csv` | Per-frame candidate, derivative and complete-band costs and widths |

The gallery renderer produces GIFs from these movies. Generated images,
movies, GIFs and run data remain ignored artifacts. Normalized axes use unit
[1], and every animation title displays the changing κ rather than a physical
time. The source diagram is a static cross-section illustration.

CPU_TIME records three repeated floating-point primal and primal-plus-tangent
evaluations per frame, and one complete reference/band construction. The ledger
retains absolute extra cost and multipliers against both baselines. Setup,
rendering and encoding are excluded from those scientific timers. This is
instrumentation for a reference example, not a controlled performance ranking.
All scientific buffers are preallocated; no trajectory is stored.

## Independent verification

[The focused test](../test/test_debye_sphere.f90) evaluates the radial Yukawa
Green-function convolution in binary128 using separate antiderivatives split
at the observation radius. This oracle uses volume integration rather than
the production boundary-matching/series derivation. It checks strict reference
intervals at the center, small positive radius, boundary and exterior, Gauss's
law, matching of potential and field, directional derivatives, transpose
pairing, whole radial/parameter boxes, and connected bands on nonuniform cells
crossing the source boundary. Invalid domains, NaNs, duplicate knots and
incompatible buffer shapes are rejected.

The isolated native focused oracle passes. Registered CPU/Debug, upstream
arithmetic, independent proof review and rendered-gallery evidence belong to
the controller's integration packet; they are not inferred from that local pass.
