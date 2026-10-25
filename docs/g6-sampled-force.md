# External sampled physical force diagnostic

The optional G6 physical reader reports an additional external-row diagnostic:

    f_i = (mu0 J_i) cross B_i - mu0 grad(p)_i
    force_rms = sqrt(sum_i |f_i|^2 / number_of_rows)

All vectors are Cartesian at the same physical position. B is in T, mu0*J
in T/m, and both mu0*grad(p) and the force residual in T^2/m. The input
already includes mu0, so no second factor is applied. The divisor counts
physical rows, not Cartesian components. Samples have equal weight; the
result is not a volume norm or a continuous bound.

The recovered case contract provides no physical-force normalization or
acceptance tolerance. This metric is therefore absolute and diagnostic only.
It does not add or relax an admission gate, and a nonzero value does not fail
the existing B/p comparison. Unavailable force arithmetic is named explicitly
and its value is NaN, while the independent B/p admission retains its own
status. Unsafe products/subtractions, nonfinite inputs, empty accumulation
and an unrepresentable RMS do not yield a finite partial result.

The parsed columns describe the external/reference state. The current KIN
sampler supplies only B and p, without native J and grad(p) on those points.
This diagnostic therefore measures external sample consistency, not a common
KIN/reference force difference, a full imported state, NN2 D3 or a validated
force-balance certificate. Matched physical points, profiles, units, wall,
flux and angle conventions remain prerequisites for cross-code interpretation.
No VMEC/GVEC provider chain or new reference run is qualified by this change.

The public tests independently differentiate the Solovev flux and pressure,
which give curl(B)=(0,p'R,0) and zero exact force. A genuinely 3D control uses
B=(y,2z,3x), mu0*p=10+x+2y-3z in SI coordinates. Its analytic curl is
(-2,-3,-1), and the force vectors at e_x,e_y,e_z are (-10,4,3),(-1,-3,6),
(1,-2,-1), giving RMS sqrt(59) T^2/m. Component signs, a proper rotation and
physical rescaling are checked. For example B-scale2 and length-scale2/3
imply current-scale3, pressure-scale4 and force-scale6. Roundoff allowances
inside these analytic tests are not production acceptance tolerances.

This independent comparison diagnostic lives in test support and does not
change a production solver kernel, generated operator, derivative contract,
NN2 normalization, current or pressure/sign convention.
