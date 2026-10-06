# Q0: periodic quantum two-electron quench

The runnable reference is an interacting quantum two-electron Coulomb quench
in a periodic **three-dimensional** box with a uniform neutralizing background.
It evolves a many-body wavefunction, not particle trajectories or an imposed
exchange/correlation hole. The default box, coupling, initial preparation and
sixteen relative-position bins match [the classical pair experiment](classical-pair.md).
Both runs start from an uncorrelated opposite-spin pair, and the same bare
Coulomb coupling is switched on at time zero. Two opposite spins isolate
interacting diffraction and correlation from same-spin exchange; localized
exchange scattering is the separate case [Q1](coulomb-encounter.md). The quench
is a microscopic correlation reference, not incoming wave-packet scattering,
thermodynamic Debye shielding or a finite-temperature equilibrium plasma.

Build and test through the repository's native CMake/CTest path, then run
`quantum_reference_demo --snapshot 31` for a middle-frame PNG or
`quantum_reference_demo` for `output/quantum-correlations.mp4`. The gallery
renderer produces its GIF. `--config input.nml` changes the declared parameters
and supplied complex initial amplitudes. Generated media are not committed.

## Model and preparation

The shared Fock engine admits arbitrary supported spin-orbital lists and
fixed particle count subject to explicit state and memory budgets. This gallery
selects the useful N=2 benchmark with one up-spin and one down-spin electron;
that restriction belongs to this experiment, not the engine. Spatial orbitals
are

\[
 \varphi_{\mathbf n}(\mathbf x)=L^{-3/2}
 e^{2\pi i\mathbf n\cdot\mathbf x/L},\qquad
 |n_\alpha|\leq K,\quad\alpha=x,y,z.
\]

The default cube half-width K=2 retains 125 spatial momenta. K=1 and
K=3 retain 27 and 343. The earlier spherical 7/19/27 cases remain independent
small full-Fock comparison fixtures. The nonrelativistic Hamiltonian has kinetic energy
\(\hbar^2|\mathbf k|^2/(2m)\) and the **unscreened** periodic Coulomb
coefficients \(4\pi g/(L^3|\mathbf q|^2)\) for nonzero transfer. The removed
zero mode declares the neutral background and potential gauge; it is not
Debye screening. The shared Fock generator owns the coefficients and CAR signs.

Total momentum and both spin counts are preserved. The initial zero-momentum
Slater determinant therefore lies in the exact invariant sector
\(|\mathbf n\uparrow,-\mathbf n\downarrow\rangle\), of dimension 125 in
the default basis. A physical ordered-sector constructor uses the same retained
kinetic/Coulomb leaves as the full Fock engine and avoids a large bitmask basis
or pair-integral tensor. Small full-Fock selected-block tests check its
coefficients, creation-order signs and every outgoing term independently.
Reducing this exact sector does not approximate the admitted finite-basis model.

Default parameters are L=2*pi, m=hbar=1, g=2 and t=0..16, in declared reduced
units. A reference length ell0 and mass m0 give the time unit
`t0=m0*ell0**2/hbar0` and coupling unit `hbar0**2/(m0*ell0)`.
This is not a 10 keV plasma example or a thermodynamic screening calculation.

The default state is the k=0 opposite-spin Slater determinant. It has uniform,
independent initial positions and zero momentum, matching the positive initial
Wigner preparation of the classical experiment. General complex amplitudes
may instead be supplied by their integer momentum triplets in the admitted
invariant sector. The app rejects momenta outside the retained cube and
duplicate triplets. It preserves
raw supplied amplitudes; physical expectation values divide by their exact
initial squared norm. There is no silent normalization or known-solution input.

The app amplitudes refer directly to ordered states
`a_dagger(k,up)*a_dagger(-k,down)|0>`. Full-Fock oracle amplitudes instead use
increasing orbital order, requiring the creation-order sign when comparing
those representations. The physical constructor includes no extra canonical
orientation sign in the app; applying that sign twice changes the interference.

## Observables and matching bins

Write the physical sector coefficients as c_n and their initial squared norm
as Z. The unlike-spin relative pair density, normalized to the uniform law,
is

\[
 g_{\uparrow\downarrow}(\mathbf r)
 =Z^{-1}\left|\sum_{\mathbf n}c_{\mathbf n}
 e^{2\pi i\mathbf n\cdot\mathbf r/L}\right|^2.
\]

Both one-body densities are exactly V*n_sigma=1 by normalization, translation
symmetry and fixed spin count. Each paired spin orbital has occupation
|c_n|^2/Z. Same-spin pair density is zero because this experiment has one
particle of each spin. These are structural consequences of the named model,
not analytical solutions used by the time checker.

The plotted pair quantity averages over a relative-x window of width L/16 and
the full transverse torus. Its value is the integral of the same three-dimensional
pair law as the classical histogram. A Fourier transfer d acquires the exact
factor product over axes of `sinc(pi*d_alpha*w_alpha)`, where
w=(1/16,1,1). Complete transverse periods eliminate nonzero transverse transfers
exactly. Bin centers are `x/L=-1/2+(j-1/2)/16`, j=1..16.
No particle is confined to a line or sheet.

The unlike-spin cumulant is C_updown=g_updown-1. The same-spin Wick cumulant
is -1+|sum_n p_n exp(2*pi*i*n.r/L)|^2, with p_n=|c_n|^2/Z, averaged over the
same window. It records fixed-number correlation beyond the Wick product;
it is not a same-spin pair-density hole produced by two up-spin particles.
Wick products are observables here, not an enforced dynamics closure.
Fixed-N pure states and fractional grand-canonical occupations have different
contraction identities, so a fixed-N identity is never applied to a
finite-temperature ideal-gas Wick state.
The noninteracting comparison is computed by the same Fock and evolution
engine at coupling zero with identical amplitudes and bins.

## Deterministic numerical bands

The shared Hermitian evolution checker supplies a continuous-time Hilbert
radius r around the reconstructed candidate for its supplied binary64
Hamiltonian. Strict assembly additionally encloses the mathematical finite-basis
Coulomb matrix. If delta bounds the Frobenius norm of the difference between
these matrices, unitary perturbation gives the additional radius
`t*delta*sqrt(Z)/hbar`. Exact represented step times, including multiplication
rounding in the displayed time label, determine this budget.

The checker now defaults to deterministic aggregate-rounding actions (Tier D)
and retains dense outward actions (Tier S) as an audit. This preserves the
finite-state and assembly budgets above. The matched timing fixtures at t=1
are distinct from the longer gallery defaults; see
[the strategy and measurement ledger](error-control-strategy.md).

For a candidate squared norm z and a bounded observable A,

\[
 |\langle\psi,A\psi\rangle-\langle\psi_h,A\psi_h\rangle|
 \leq\|A\|(2\sqrt z\,r+r^2).
\]

Division by Z gives the normalized expectation budget. Occupation-projector
norm is 1. For a point pair-density observable in this sector the norm is at
most the number of momenta. Full-transverse averaging makes its matrix block
diagonal in (n_y,n_z); the norm is at most the largest number of n_x values
in a block, **5 for the default 125-mode cube**. This bound follows from the
matrix row sums and |sinc|<=1 and holds for arbitrary supplied sector amplitudes.
It does not assume their future solution. The exchange-amplitude observable
has norm 1, so its squared-modulus error is bounded separately by the
expectation budget times `(1+z/Z)`.

FortSym generates strict complex products, normalization, Fourier phases,
observable budgets and curvature sums from the retained source. Strict
candidate arithmetic and the state-to-observable budgets are both included.
The general full-RDM path independently contracts gamma/Gamma using the
shared Wick source; it does not replace Gamma by a closure.
It requires distinct momentum/spin orbitals and a complete, unique list of
ordered pair labels. Repeated labels invalidate the packet explicitly.

Connected bands enclose the whole displayed center-coordinate cells. For
Fourier coefficients a_d, the curvature bound is
`sum_d (2*pi*d_x)**2*(abs(Re(a_d))+abs(Im(a_d)))`.
Add M*gap^2/8 using the maximum actual represented plotting gap. This concerns
the smoothly translated slab functional; it does not reinterpret a coarse
histogram as a point density. Discrete orbital occupations use error bars on
a log axis, restricted to the orbitals that reach occupation 1e-4; a fourth
panel follows the depletion 1 - n_0 and one orbital per shell |k|^2 = 1-3
with their bands (display floor 1e-6). The uniform one-body density is not
plotted. The movie runs at 24 fps and shows each certified frame for three
video frames (no interpolation), then holds the last frame for 1 s.

The default requested radius is 1% of scale 1 plus an absolute floor 1e-12.
A Hamiltonian row-norm scale chooses a starting work estimate; the app
doubles steps until all displayed
bands and occupations pass that target. Valid wider certificates and failed
accuracy targets remain distinct. Finite-basis-to-continuum error, physical
model inadequacy and input uncertainty are **unavailable**. Basis differences
are refinement indicators, not certified remainders. Classical sampling
intervals retain their separate probabilistic meaning.

## Derivatives, evidence and cost

Generated smooth observable leaves expose primal/JVP/VJP with the real
Euclidean pairing. Shared evolution supports equation-level Hermitian tangent
and initial-state adjoint paths. This app does not report observable sensitivity
bands or a certified Hamiltonian-parameter gradient; their status is unavailable.

The independent observable test uses a complex two-determinant state and
real128 plane-wave amplitudes, analytic slab integrals, a perturbed-state norm
bridge, complete-cell chord containment, transverse-mode elimination, a
Hamiltonian norm-defect oracle and finite-difference/JVP/VJP checks. It compares
the exact-sector and full-RDM contraction without using a time solution in the
checker. Shared Fock tests independently verify CAR signs and invariant blocks;
shared evolution tests independently verify its time certificate.

Outputs retain resolved namelists, full observable and occupation intervals,
and a separate sixteen-bin CSV for matched comparisons. `quantum-performance.csv`
records candidate-only initialization/solve, certified initialization/solve,
extra time and multiplier, all observable postprocessing, common model setup,
point-only and strict matrix assembly, and adaptive-search cost. Warm totals
reuse supplied matrix packets; cold totals include the corresponding constructor
and common setup costs. The total including all bands is reported against
solve alone. CSV/movie encoding is excluded explicitly.
These are single local CPU trials, not a speed ranking. The gallery manifest
retains executable/source/compiler/dependency provenance. The certified
finite-basis reference currently favors clarity and correctness; interval
residual and postprocessing costs remain measured optimization targets.

The reviewed local Release run records 1.155 s for the matched quantum
candidate initialization and evolution and 104.967 s including certification
and all plotted bands, 90.89 times the solve alone; the adaptive first result
costs 157.35 s before rendering, with maximum display radius 0.002757 on
scale 1 at 125 modes and 16 steps per frame. The classical companion costs
33.508 s plus 2.318 s for sampling and diagnostics (factor 1.069), which is
not a deterministic trajectory certificate. These single CPU trials are the
optimization baseline for the dense interval residual action, not a
workhorse performance claim.

Retained cutoff fixtures compare the same physical quench and sixteen bins
at t=0,1,4,16. At t=16 the earlier spherical differences are about 0.3883
for 7 versus 19 momenta and 0.06767 for 19 versus 27. The larger cube comparison
gives about 0.05358 for 27 versus 125 and 0.00517 for 125 versus 343. The
125-mode strict numerical bin radius is about 0.00225 at 1024 steps; the
343-mode run has an explicitly unavailable certificate and a 1024/2048-step
difference about 2.4e-8. These indicators improve the resolution evidence
without certifying a continuum remainder. All deterministic bands still cover
the admitted finite-basis model only. Small fixtures independently evolve the
physical signed Hamiltonian in real128 and check analytical slab integrals.
