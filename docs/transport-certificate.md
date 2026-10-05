# Arbitrary-source damped transport certificate

This increment separates a general continuous residual certificate from the
laboratory's exact-profile comparisons. The production checker consumes the
supplied source and candidate coefficient bytes; it never receives an exact
solution. Mathematical claims below concern prescribed E=0 transport on a
finite velocity domain. Full VPFP stability, velocity tails and physical
nonnegativity have separate obligations.

## Continuous problem and input interpretation

Solve `L f = s`, `L=nu+v*d_x`, with constant nu>0, spatial period `2*pi` and
v in `[-1,1]`. Every supplied binary64 parameter, velocity node and coefficient
is interpreted as its exact real value. The strictly increasing grid has exact
endpoints -1 and 1. Source and candidate coefficients interpolate linearly
between velocity nodes and use a finite spatial Fourier expansion:

```
s(x,v) = s_0(v) + sum[k=1..K] (s_c,k(v)*cos(k*x)+s_s,k(v)*sin(k*x)),
fbar(x,v) = a_0(v) + sum[k=1..K] (a_k(v)*cos(k*x)+b_k(v)*sin(k*x)).
```

The mode-zero sine coefficients are exactly zero. Missing spatial modes are
zero by the declared source definition, rather than ignored source tails.
An analytic expression used to construct nodal values is only a constructor;
the certificate describes the resulting P1/Fourier source. Certifying an
external analytic source additionally requires a validated difference bound
and its propagation through the inverse. That interface is not supplied here.
No velocity boundary condition is imposed because this operator has no
velocity derivative.

The source is continuous in velocity and smooth in x. The characteristic
resolvent

```
f(x,v) = integral[t=0..infinity] exp(-nu*t)*s(x-v*t,v) dt
```

exists absolutely, is periodic, and satisfies L f=s by integration by parts.
Spatial derivatives may be passed under the integral because the finite
Fourier source and its derivatives are uniformly bounded. This includes v=0,
where the equation is nu*f=s. The resolvent has continuous sup-norm bound
`||L^-1|| <= 1/nu`. A periodic differentiable homogeneous solution has zero
maximum and minimum by the damped maximum principle, giving uniqueness.
The production certifier does not evaluate this solution formula.

## Whole-cell residual proof

Use residual `r=s-L*fbar`. For mode k its cosine and sine coefficients are

```
r_c = s_c - nu*a - k*v*b,
r_s = s_s - nu*b + k*v*a.
```

They are quadratic polynomials in each velocity cell. Write
`t=(v-v0)/(v1-v0)` and express each polynomial in the nonnegative Bernstein
basis `(1-t)^2, 2*t*(1-t), t^2`, whose sum is one. For r_c the three exact
coefficients are

```
B0 = s_c0 - nu*a0 - k*v0*b0,
B1 = (s_c0+s_c1-nu*(a0+a1)-k*(v0*b1+v1*b0))/2,
B2 = s_c1 - nu*a1 - k*v1*b1.
```

For r_s replace a with b in the damping term and replace the last product
term by `+k*v*a`. Each coefficient is rigorously enclosed by a FortSym-emitted
interval leaf. The maximum of their magnitudes bounds that harmonic over the
entire velocity cell. Summing the cosine magnitudes, the nonzero-mode sine
magnitudes and the mean with outward arithmetic gives a bound R_j valid for
all x in the cell. Taking `R=max_j R_j` proves

```
||f-fbar||_infinity <= R/nu.
```

This controls every spatial and velocity point of the declared problem,
including cell interiors and the continuous complement. Nodal inverse
accuracy, a small nodal residual or mesh agreement does not replace the
whole-cell argument. The source and candidate can be independently changed;
a claimed radius is checked against freshly evaluated raw data.

## Workhorse and strict audit

The workhorse prepares compact rigorous Bernstein coefficients, with fixed
scratch and one cell/mode traversal. It performs no transcendental evaluation,
point sampling, exact-state comparison or global matrix construction.
Both current paths use FortNum outward interval arithmetic. This is an
initial practical rigorous path; a generic center-plus-aggregate
roundoff replacement is not asserted to exist in this operator yet.

`strict=.true.` selects a distinct direct interval expression evaluation over
each whole velocity cell. It encloses the P1 coefficient ranges and evaluates
`s_c-nu*a-k*v*b` and `s_s-nu*b+k*v*a` as interval expressions. Losing their
shared-velocity correlations can make this audit radius wider than the
Bernstein radius. Both paths guarantee containment; radius equality or strict
mode being tighter is not an acceptance requirement. Bernstein coefficients
may cost more per cell; measured cost at a matched observable accuracy target
determines the practical tradeoff. The flag names the independent box audit,
not a switch from unvalidated floating arithmetic. Nonfinite inputs,
nonpositive nu, invalid dimensions/grid, nonzero mode-zero sine data and
nonfinite certificate arithmetic fail explicitly.

## Tangents and observables

For a declared source direction ds and damping direction dnu on the fixed
domain, the continuous tangent obeys `L df=ds-dnu*f`. Its numerical candidate
has the same P1/Fourier representation. The checker bounds the whole-cell
residual `r_t=ds-dnu*fbar-L*dfbar` and includes primal uncertainty:

```
||df-dfbar||_infinity <= (||r_t||_infinity+abs(dnu)*r_f)/nu,
r_f = certified primal radius.
```

Source derivative coefficients also denote exact supplied P1/Fourier data.
This certifies the solution's declared direction, rather than differentiating
iterations. Parameter positivity is required at the evaluation point; integer
mode choices and domain changes are not derivatives supported by this gate.
Generated nodal JVP/VJP checks remain distinct from a continuous physical
adjoint certificate. The continuous adjoint differential expression under
`integral dx dv` and periodic boundaries is `nu-v*d_x`.

For a bounded linear observable `Q(f)=integral q*f`, rigorously enclose
Q(fbar) and provide a verified upper enclosure of `||q||_L1`. The observable
helper widens the candidate interval by this norm times r_f. It validates
finite arithmetic, not the scientific correctness of a caller-supplied norm
or candidate integral. A useful nontrivial example is

```
Q(f) = (1/(2*pi))*integral[-1,1] integral[0,2*pi] cos(x)*f(x,v) dx dv,
Q(fbar) = (1/2)*integral[-1,1] a_1(v) dv.
```

A norm bound of 2 follows from `abs(cos)<=1`; the generated P1 cell integral
and outward sum enclose the candidate observable. The same construction with
tangent coefficients gives its directional derivative interval. Report state,
field/source, cell quantities, this observable and requested tangent outputs
with their respective supported scopes and statuses. An L2-to-pointwise
assumption is never used.

## Prepared packets and verification

Certificate preparation establishes an enclosure for the supplied source/
candidate/parameter packet. The caller must keep that packet and the returned
certificate immutable. The sampling wrappers check layout and touched data but do
not rescan or hash every coefficient per plotted cell. Callers must freeze
those bytes until postprocessing finishes and retain their provenance; any
change requires new preparation. `transport_certificate_accepts` deliberately
recomputes the global certificate, so its cost must be reported rather than
hidden in a sampling benchmark.

Whole-x-cell sample bands enclose the reconstructed candidate with generated
interval interpolation and Fourier accumulation, then widen by the prepared
uniform radius. Plot these as full-cell step envelopes. Source input ranges
have their own exact P1/Fourier interval evaluation and no solution-error
radius. The tangent packet, including its supplied coefficients, must also
remain immutable after preparation.

The registered native test `general_transport_and_tangent_oracles` passes in
CPU Release and Debug. It checks independently inverted rational node/JVP/VJP
fixtures, an arbitrary unsolved quadratic residual's exact Bernstein
coefficients and midpoint values, binary128 continuous mode solutions between
velocity nodes, both arithmetic paths, and full-cell state/tangent bands.
The weighted observable has independent identities
`Q=atan(1/nu)/4` and `dQ/dnu=-1/(4*(nu^2+1))` for its test-only source.
Invalid grids/domains, nonzero mode-zero sine data, NaNs and corrupted
candidate/source/radius packets are rejected or require wider recomputed
bounds. Known solutions belong only to these independent benchmarks.
A claimed analytic-source certificate needs an additional oracle for the
source reconstruction difference.

Integrated [measurements](general-transport.md) record solver, primal/tangent
certificate preparation, acceptance rechecks, batched postprocessing, achieved
observable/gradient widths and attributed buffers separately. Apply the
default 1% per-quantity accuracy and amplified arithmetic
budgets in [certification](certification.md); do not replace a proved rounding
bound by that percentage. Lean is optional assurance and is not required to
promote a valid Fortran continuous certificate.

An independent proof audit passed against base
`f30129e0bb59b9a213a1f30a6001e2117d3e6fd3` and the six-file frozen source
manifest digest
`10a4d1871e515f7d73f0cf5dd45dfee368953db03579a57e84ffcfbdbbd6589a`.
The audit covers the continuous inverse, Bernstein and box residuals, tangent
coupling, sampling and observable argument under the stated arithmetic and
immutable-packet premises. This evidence paragraph updates reporting after
the frozen mathematical text; no proof or numerical definition changed.
Integrated Release/Debug gates pass 4/4 and optional formal gates pass 5/5.

Chris&AI
