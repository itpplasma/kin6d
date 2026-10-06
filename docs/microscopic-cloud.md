# Exploratory isolated Coulomb cloud

This animation is a **candidate illustration**, not an implemented periodic
Debye-shielding benchmark or a promoted microscopic production solver. It shows
32 real, mutually interacting particles moving in three positions and three
velocities. Every electron has normalized mass 1 and charge -1; the Coulomb
prefactor is 1. Particle weights and screened pair potentials are absent.

A fixed, uniformly charged positive sphere has radius 3 and total charge 32.
Its force follows Gauss' law inside and outside the sphere; its surface is not
a reflecting wall. Two runs start from the same deterministic positions and
velocities (a nonthermal Fibonacci shell at radius 2.5, with deterministic
velocity components of amplitude 0.4). The second additionally contains a fixed negative point impurity
with electron--impurity coupling 2. The isolated system need not be neutral
with that impurity. There is no periodic minimum-image approximation, cutoff,
thermostat, effective collision operator, radiation reaction or particle self
force. Mutual electrostatic forces are instantaneous: this is a declared
nonrelativistic Coulomb limit, not full Maxwell--Lorentz dynamics.

The two projected position panels show actual three-dimensional dynamics.
The dotted circle shows the sphere's equator at z=0, not a confinement boundary
or the full projected silhouette. The radial empirical CDF uses actual 3D
radii. It is neither a continuum pointwise density nor a statistical confidence
statement. This nonthermal, small-N confined cloud does not establish an
equilibrium screening length, a Vlasov limit or tokamak collisionality. Those
claims require the periodic, weak-coupling, ensemble/size/box gates in
[the microscopic plan](microscopic-many-body.md).

FortSym emits the shared Coulomb force/potential and its primal/JVP/VJP, and
the cloud's background, subflows and observable leaves. The Gauss field is
continuous at the sphere surface; its spatial derivative changes there.
The norm and point-force derivative require a nonzero separation. Unresolved
pair or impurity encounters fail explicitly rather than using softening.
Force-leaf derivatives use Euclidean pairings. Trajectory derivatives and
their enclosures are unavailable in this exploratory executable.

FortNum owns Strang kick/drift composition. The default step is 0.001, with
100 steps per animation frame, over normalized time 0 to 12. The energy curve
is a Hamiltonian conservation **diagnostic**, never a trajectory enclosure.
Numerical trajectory, model-reduction and statistical bounds are explicitly
unavailable. No interval band is manufactured from the energy defect or a
step-refinement difference. A future certifier must close force, time and
stability budgets before this candidate can become a certified result.

Registered tests retain independent Gauss-law/potential oracles, Coulomb
action--reaction, an exact single-electron harmonic orbit, an independently
derived background directional derivative and the JVP/VJP pairing identity.
Those tests verify the named cases; they do not certify a 32-electron trajectory.

Run `./build/microscopic_cloud_demo` after the CPU build, or add
`--snapshot [frame]` for one of frames 1 through 241 without an encoder.
Since 2026-10-06 the default records 241 frames of 50 substeps instead of
121 frames of 100 substeps: the same 12000 steps (identical trajectories and
energy defect), sampled twice as often for a 10 s movie at 24 fps plus a 1 s
hold. The cloud panels use equal aspect, so the sphere equator is a circle.
The Fortran postprocessor writes `output/microscopic-cloud.mp4`, a PNG,
particle states and CSV diagnostics. These exploratory CSVs are retained
diagnostic evidence; they do not establish the planned canonical FortIO
scientific-result interchange contract. Generated media stay outside Git.
The diagnostic CSV records per-frame solve-only CPU time and extra Hamiltonian
diagnostic time with its multiplier. Certificate cost and derivative-solve
cost are unavailable. No performance ranking follows from these samples.
