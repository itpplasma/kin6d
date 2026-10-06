# Matched classical two-electron Coulomb quench

This teaching candidate evolves genuine three-dimensional classical particle
pairs, not PIC macro-particles. It is the classical companion to the
opposite-spin quantum preparation in [quantum-reference-gallery.md](quantum-reference-gallery.md).
It is not a thermodynamic Debye-screening experiment.

Each independent preparation contains two equal-mass like charges, with
`L=2*pi`, `m=1`, Coulomb coupling `g=2`, and zero initial momenta. The two
positions are independent uniform torus positions. A uniform neutralizing
background removes the zero Fourier mode. Spin is a label in the classical
experiment and does not change its force. The quantum companion starts with
one spin-up and one spin-down particle in the zero-momentum orbitals; its
initial positive Wigner preparation has the same classical position and
momentum law. Classical markers represent particle trajectories; quantum
plots represent probabilities.

## Periodic interaction and energy ledger

The pair Green function is the periodic Ewald Coulomb sum:

\[
G(r)=\sum_{n\in\mathbb Z^3}\frac{\operatorname{erfc}(\alpha|r+nL|)}{|r+nL|}
 +\frac{4\pi}{L^3}\sum_{k\ne0}\frac{e^{-|k|^2/(4\alpha^2)}}{|k|^2}
 \cos(k\cdot r)-\frac{\pi}{\alpha^2 L^3},\quad k=2\pi n/L,\quad n\in\mathbb Z^3.
\]

The force is `-g grad G`. The last constant is the uniform-background
contribution. Position-independent self-image constants do not affect this
quench's force or Hamiltonian-drift diagnostic. There is no Yukawa potential,
minimum-image force replacement, artificial wall or softened encounter.
Wrapping selects a chart for the periodic sum and the displayed particles;
it does not restrict the particles to a line or sheet.

The candidate truncates real lattice components at two and reciprocal
components at four, using `alpha=.4`. One member of every reciprocal `+/-k`
pair is evaluated with doubled weight. Generated phasor recurrences avoid
one trigonometric call per reciprocal vector. The Ewald decomposition is an
algorithm, not an added physical screening length. Independent splitting
checks support convergence, but no rigorous truncation remainder is currently
provided. The certified status remains **unavailable**.

Translation invariance permits exact reduction to center of mass and relative
coordinate `r=x1-x2`: the center is fixed and `r''=2 F(r)`. The monitored
Hamiltonian is `|r'|^2/4+g G(r)`. Generated velocity-Verlet leaves advance the
relative state; the frame interval is divided into steps no larger than
`.02` and `sqrt(r_initial^3/4)/50`. The schedule resolves close initial pairs
more finely without rejecting, replacing or softening them. A failed trajectory
writes its sample/frame to `classical-pair-failed.txt` and stops the gallery;
a failed sample is never silently discarded. Numerical timestep and
interaction-truncation errors have **unavailable** deterministic enclosures.

## Sampling uncertainty and common observable

The sample size is fixed at 4096 before examining results. The deterministic
Park--Miller pseudo-random stream starts at integer seed 104729; its integer
arithmetic is portable. It approximates independent uniform sampling; the
confidence interpretation assumes the declared iid uniform preparation law,
not a theorem that a finite pseudo-random sequence is genuinely independent.

Sixteen bins partition relative `x/L` in `[-1/2,1/2)`, and each bin integrates
over the full transverse torus. These are rectangular **three-dimensional
slab averages**, with widths/L `(1/16,1,1)`, not pointwise line densities. The
reported correlation is `g_bin=16 count/4096`. The initial law has `g=1`, but
the initial finite histogram retains its shot noise. Quantum Fourier voxel
averaging uses exactly these same cells and normalization.

For every one of the 16 bins at all 61 fixed frames, the plotted sampling
radius is `16 sqrt(11/(2*4096))`. Hoeffding's inequality and a union bound give
failure probability at most `2*16*61*exp(-11)<.05` for the numerical-map bin
probabilities under iid sampling. The same draws across frames are dependent;
the union bound does not require independence between frames. These are
simultaneous statistical confidence intervals, **not deterministic numerical
error bands**. They do not cover timestep error, Ewald truncation, model error
or pseudo-random discrepancy. No extrapolation to exact physical trajectories
or arbitrarily selected frames is certified.

## Derivatives and verification

Retained `gen/gen_classical_pair.f90` definitions emit primal/JVP/VJP leaves.
The periodic force exposes its position JVP on a fixed wrapping chart.
Production force and its tangent share the same Ewald definition. Histogram
binning is discontinuous at cell faces and its ordinary derivative is
unavailable; no differentiability claim is made for the histogram or adaptive
schedule. A smooth statistical estimator and equation-level trajectory
sensitivities are subsequent gates, rather than fabricated finite differences
reported as production AD.

`test_classical_pair` checks force oddness, potential evenness, agreement of
different Ewald splittings, an independent energy finite-difference gradient,
force JVP finite differences, primal-leaf JVP/VJP duality, and a short-run
Hamiltonian diagnostic. These are behavioral tests, not certificate proofs.

## Reproduction and performance

`classical_pair_demo [--snapshot [frame]]` accepts frame indices 1 through 61.
The movie shows each of the 61 shared frames three times at 24 fps (8
physical frames per second, no interpolation) and holds the last for 1 s; the
frame count stays tied to the shared `&periodic_pair` configuration and the
Hoeffding coverage.
It outputs `classical-coulomb-pair.mp4`, or the snapshot PNG, alongside bin,
particle and timing CSV files under `output/`. The left panel shows the 3D
particle projection `(x+z/4,y+z/6)/L` for twelve independent pairs. The right
panel shows all 4096 pairs' slab correlation with statistical interval bars.
Frames carry physical time from zero to sixteen; they are not an equilibrium
parameter sweep. All arrays are prepared before encoding and callbacks clear
the canvas each frame.

The timing CSV separates trajectory solve CPU time from histogram/energy
extra CPU time. Rendering and movie encoding are additional, separate costs.
Ewald setup is performed once and reciprocal weights are reused. Generated
sources retain the FortSym provenance; large movies and sample CSV files stay
outside source history. No performance-leadership claim follows from this
small reference.
