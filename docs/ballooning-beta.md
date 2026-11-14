# Real-time local ideal-ballooning beta estimate

**Scope:** a genuine reduced ideal-MHD stability calculation, not an empirical
Troyon coefficient or a global no-wall kink calculation. This is a local,
large-aspect-ratio, shifted-circular, high-toroidal-mode-number `s-alpha`
model. It **does not replace** MILCA's 3D and cylindrical energy-principle
solvers; it is a fast lower-fidelity diagnostic and teaching model.

## Derivation and model

The ideal-MHD energy principle, in the standard shifted-circle local
ballooning limit, reduces at marginality to

```math
\frac{d}{d\theta}\left[A(\theta)\frac{dF}{d\theta}\right]
+V(\theta)F=0,
\quad A=1+\Lambda^2,\quad
\Lambda=s\theta-\alpha\sin\theta,
\quad V=\alpha(\cos\theta+\Lambda\sin\theta).
```

The dimensionless inputs are the local shear
`s=(r/q)(dq/dr)` and pressure drive
`alpha=-2 mu0 R q^2 (dp/dr)/B_T^2`.
Equivalently, the variational energy is
`delta W ~ integral [A(F')^2-V F^2] dtheta`.

We integrate the even marginal equation using `G=A F'` and fourth-order
Runge–Kutta:

```text
F' = G/A,   G' = -V F,   F(0)=1,   G(0)=0.
```

A zero of the even solution on positive theta is the Newcomb
zero-crossing indication of instability. A scan from alpha=0 detects the
**first** crossing (first-stability boundary); it deliberately does not
identify later second-stability branches.

For a **declared** parabolic pressure profile
`p(r)=p0(1-r^2/a^2)`, circular area averaging gives
`beta_avg=mu0*p0/B_T^2=beta_axis/2`. Eliminating `p0` gives

```math
\boxed{\langle\beta_T\rangle_{\rm crit}
= \frac{a^2\alpha_{\rm crit}(s)}{4Rrq^2}.}
```

The interface `ballooning_beta_limit(s,q,r,a,R,...)` returns the **volume
(or circular cross-section) averaged toroidal beta of that specific profile**,
not a universal beta limit. In particular this model tests *high-n ideal
ballooning*, not external/internal kink, resistive modes, edge-localized
modes, equilibrium feasibility, nor the empirical Troyon `beta_N`.

## Accuracy and status

The fast integrator is deterministic, allocation-free, has no I/O, and
offers a `bind(C)` scalar entry point for a WebAssembly host. The convenience
call checks step refinement (1024 to 2048 steps over 16π) and angle-horizon
sensitivity (16π to 32π). The returned difference is a **numerical estimate,
not a rigorous error enclosure or a physical-model uncertainty**.

Return codes: `BALLOON_OK=0` means a finite-domain first crossing was located
and the two sensitivity checks agree to 2%; `INVALID=1`, `NOT_BRACKETED=2`,
and `UNRESOLVED=3` are not physical limits. A search without a zero cannot
be interpreted as a proof of stability, and zero-crossing sensitivity must
not be conflated with a mathematically certified continuum bound.

For reproducible benchmarking record `(s,q,r,a,R)`, `alpha_max`, angular
domain, steps, compiler revision and timing. Do not claim real-time WASM
performance until a target-toolchain benchmark is run. No full-WASM build is
delivered here; all downstream dependencies should be eligible to run there.

## Independent checks

- At alpha=0 the exact even solution is F=1, with no zero.
- For s=1, the well-known first-stability boundary is near alpha=0.6;
  numerical outputs at 16π and 32π are compared, not fitted to this number.
- Doubling q at fixed local shear/profile reduces this profile-based beta
  estimate by a factor of four.
- The P1e `milca_ideal_cylinder` and global energy-principle work provide
  distinct future model-level comparators, not the same stability prediction.

## Sources and provenance

- J. W. Connor, R. J. Hastie, J. B. Taylor, classical high-n ideal-ballooning
  theory, *Physical Review Letters* **40**, 396 (1978).
- J. P. Freidberg, *Ideal MHD* (Cambridge University Press, 2014),
  high-n ballooning energy principle.
- [Classic Problems in MHD, lecture 40](https://magnetohydrodynamics.physics.wisc.edu/lecture40.html),
  including the shifted-circle reduction and s-alpha definitions.
- [GPEC.jl ballooning documentation](https://openfusiontoolkit.github.io/GPEC/dev/ballooning/),
  the general ideal-ballooning variational equation.
- [VAFT formula reference](https://vest-tokamak.github.io/vaft/reference/formula/stability/),
  independently documented marginal s-alpha ODE and zero-crossing convention.

Implementation derived here from public equations, not copied from another
implementation. No non-public code of others was used or in context.
