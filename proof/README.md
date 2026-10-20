# Continuous reference proofs

This optional package is checked with the explicit `formal` CMake preset:
`cmake --preset formal`, `cmake --build --preset formal`, `ctest --preset formal`.
Ordinary CPU/debug builds retain numerical certification without requiring
Lean/Lake. The theorem scopes below do not extend to the general solver.

Lean checks a pointwise error enclosure for the bounded-velocity prescribed-field
Krook reference. For `v ∈ [-1,1]`, define

```
M(v) = 3/4 * (1-v²)
f_a(x,v) = M(v) * (1+a*cos(x))
forcing_a(x,v) = -v*M(v)*a*sin(x) + f_a(x,v)
```

The equation is `v*d_x f + f = forcing_a`, with spatial period `2*pi`,
`nu=1` and prescribed `E=0`. The total forcing includes the fixed reservoir
term. The proofs establish existence through the explicit profile, uniqueness
among differentiable periodic solutions, and nonnegativity for `|a|≤1`.
For the reconstructed candidate `f_b`, they bound the continuous residual by
`3/2*|b-a|` and the error at every spatial and velocity point by the requested
radius. The residual estimate covers the entire domain rather than sampled
nodes.

`PeriodicTransport.lean` proves the underlying maximum principle for arbitrary
positive `nu`, arbitrary real velocity and positive spatial period. A periodic
continuous function attains a global maximum; its spatial derivative vanishes
there, so damping bounds its amplitude. Applying this argument to the negative
function gives the absolute-value inequality, and applying it to the difference
between a candidate and exact solution gives the residual-to-error theorem.
No finite-dimensional inverse estimate or assumed continuous coercivity bound
is used.

`Certificate.lean` interprets all four certificate fields as exact rational
numbers. Its validity predicate checks a sufficient continuous residual bound,
the radius, and amplitudes within the nonnegative-profile range. The soundness
theorem embeds the rational inequalities into the real numbers and applies the
continuous theorem. `tests/Accepted.lean` contains kernel-checked applications
with candidate amplitudes `3/8` and `7/16`, exact amplitude `1/2`, and radii
`3/16` and `3/32`. These fixtures are explicit approximate reconstructions,
not assertions about a numerical solver's convergence.

`PoissonMode.lean` additionally proves the normalized unit Fourier-mode
component `phi=a*cos(x)`, `E=a*sin(x)`, `rho=a*cos(x)`, `epsilon0=k=1`.
The continuous derivative signs, Poisson equation, exact net-charge
compatibility, zero-mean potential gauge, and integrated field energy
`integral E² = pi*a²` are checked. This component theorem does not supply a
nonlinear Vlasov–Poisson error enclosure.

`Bordered.lean` proves the discrete bordering lemmas of the free-boundary
Newton system (obligation FB-D1, `theory/derivations/free-boundary-gauge-vertical.md`
§4) for real matrices over arbitrary finite index types. With
`K = [H Mᵀ; M 0]`, `M : m × n` and `Z : n × k` whose columns form a basis of
`ker M` (`IsKernelBasis`):

- `lemma_B1`: `K` is nonsingular iff `rank M = card m` and `Zᵀ H Z` is
  nonsingular. `H` need not be symmetric.
- `lemma_B2`: for symmetric `H` and `rank M = card m = r`,
  `In(K) = In(Zᵀ H Z) + (r, r, 0)`.

Inertia is defined variationally: `posIndex` is the largest dimension of a
subspace on which `xᵀ A x` is positive definite, `negIndex A = posIndex (-A)`
and `zeroIndex A = card - rank A`. `posIndex_eq_card_pos`,
`negIndex_eq_card_neg` and `zeroIndex_eq_card_zero` show that for symmetric
matrices these are the counts of positive, negative and zero eigenvalues of
mathlib's spectral theorem, so `lemma_B2` is the eigenvalue statement used in
the derivation. `tests/Bordered.lean` applies both lemmas to explicit 3×3
examples, with an explicit inverse and an explicit kernel vector as independent
oracles for the nonsingular and singular cases.

These are finite-dimensional matrix lemmas only. They do not cover the PDE,
the discretisation, the vacuum Schur-complement step (FB-D2, Haynsworth with
a nonsingular pivot), or any check that a computed KKT matrix satisfies the
hypotheses (full row rank of `M`, a valid kernel basis, symmetric `H`).

## Reproduce

The package pins Lean `4.33.0`, mathlib commit
`db584cd6d46c92f209a44c0f1c829460d327499d`, and its transitive dependency
revisions in `lake-manifest.json`. From this directory, run:

```sh
lake update
lake exe cache get
./check.sh
```

The public mathlib cache is optional; building the pinned sources produces
the same kernel obligations. `check.sh` builds the theorem library, checks
both accepted numerical applications and the bordering-lemma applications,
prints their axioms, and requires Lean
to reject independently corrupted residual and radius certificates. The
accepted module also proves that the actual continuous residual and solution
error exceed the corrupted values at `(x,v)=(0,0)`. Those concrete mathematical
counterexamples establish the rejection tests' behavioral meaning.

The Fortran executable `formal_krook_certificate` solves the scalar amplitude
equation for exact source amplitude `1/2`, checks IEEE binary64 representation,
and writes the candidate's actual coefficient bits as an exact dyadic rational.
The emitted `output/ComputedKrook.lean` proves `Certificate.Valid` using
`norm_num` and applies `Certificate.sound` to the declared smooth polynomial
reconstruction. Run `./check.sh ../build/formal_krook_certificate` after the
numerical build to check that computed certificate as well as the fixtures.
The default candidate is exact, so its checked error radius is zero. The generic
Gaussian/sampled-state certifier, interval arithmetic implementation, FP model,
self-consistent electric acceleration and omitted velocity tails are outside
these theorem applications.

The printed axioms are `propext`, `Classical.choice` and `Quot.sound`, the
standard Lean/mathlib foundations. There are no project axioms, `sorry`
terms or native-decision shortcuts in the accepted proofs. Implicit undeclared
parameters are disabled so a misspelled constant cannot become a hidden
assumption. The trust base also includes Lean's kernel, the modeled equation,
input interpretation and the pinned libraries; source hashes alone prove no
scientific inequality.

Chris&AI
