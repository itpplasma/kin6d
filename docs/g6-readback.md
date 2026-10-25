# G6 physical sample readback

The `g6` mode of `test_nonnested_coupled` restores the public Solovev
analytic regression and optional physical sample comparison from d35db4e.
For example, after building its native CMake target:

```sh
./build/test_nonnested_coupled g6 3 1e6 16 0 samples.txt
./build/test_nonnested_coupled g6 3 1e6 16 0 samples.txt 1e-8
```

The last argument, when supplied, is a positive finite **absolute RMS
mu0*p tolerance in T^2** for an identically zero reference pressure. It has
no default. Zero-pressure reference rows are rejected without this explicit
tolerance, and a missing, nonfinite, zero or negative tolerance is rejected.
For a nonzero reference pressure the unchanged relative RMS criterion is
sqrt(sum((p_KIN-p_ref)^2)/sum(p_ref^2)) <= 1e-3, regardless of the optional
absolute tolerance. The magnetic-field gate is unchanged: total Cartesian
vector relative RMS <= 1e-3. Output names the pressure error kind, value,
units and applied tolerance. These are equally weighted sample diagnostics,
not volume errors, error bounds or certificates.

Every row must contain exactly 14 finite values:

    s x y z Bx By Bz mu0*p mu0*Jx mu0*Jy mu0*Jz mu0*grad(p)x mu0*grad(p)y mu0*grad(p)z

At least three rows and a nonzero magnetic reference norm are required.
Malformed rows, nonfinite values, failed physical-domain samples and invalid
norm arithmetic are rejected. The sample diagnostics use scaled sums of
squares and guard overflowing field/pressure differences; invalid arithmetic
cannot satisfy a gate through a NaN comparison.

The current callback compares B and p only. Parsed J and grad(p) do not
establish a common D3 residual. The restored analytic flux/iota check uses
KIN u=0.30,0.55,0.80; a future matched reference must evaluate iota at
s=exact_flux(u)/exact_flux(1), not compare different default surfaces.
No VMEC/GVEC provider chain, imported common state, mapped QA solve or QA
mode completeness follows from this source-reader correction. The optional
zero-pressure path is a metric contract and does not convert the Solovev
caller into a QA solver.
