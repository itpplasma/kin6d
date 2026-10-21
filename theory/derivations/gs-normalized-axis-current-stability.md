# Normalized profiles, moving axis and prescribed current

[V] Status: the identities and conditional theorem below are proved. The
continuous inverse bound, regularity and profile premises for production K6
remain **open**. This document does not promote its sampled residual to a
solution-error bound. The implementation retains certificate status 4.

## Model, trace and gauge

Let the physical meridian domain Omega satisfy 0<R_min<=R<=R_max. Fix zero
Dirichlet flux on its entire boundary; this fixes the additive flux gauge.
Write L u=-div(grad u/R), with energy norm
||v||_E^2=integral_Omega |grad v|^2/R. For fixed profile shapes P and H, set

    s=1-u/a,  f(R,s)=mu0 R P(s)+H(s)/R.

The unknown source amplitude lambda scales both derivative profiles. The
unknown a is the flux at a selected interior magnetic axis x. The coupled
equations are

    F1 = L u-lambda f(R,1-u/a) = 0,
    F2 = a-u(x) = 0,
    F3 = grad u(x) = 0,                         (two components)
    F4 = lambda integral_Omega f(R,1-u/a)/mu0-I_target = 0.

[V] These signs and normalization are those of `gs_normalized_flux` and
`gs_normalize_current` in `src/equilibrium/gs_nonlinear.f90`. A branch also
specifies the selected axis, its interior neighborhood, a separated nonzero
axis flux, and continuation provenance. A distinct stationary solution is a
distinct branch; this system does not select a branch by itself.

## Exact first variation and amplitude direction

For a perturbation (v,b,ell,h) of (u,a,lambda,x),

    ds = -v/a+u b/a^2,
    DF1 = L v-ell f-lambda f_s ds,
    DF2 = b-v(x)-grad u(x).h,
    DF3 = grad v(x)+Hess u(x) h,
    DF4 = ell integral f/mu0+lambda integral f_s ds/mu0.

[V] The variation of u(x) includes the moving evaluation point. At a stationary
axis its first-order motion term vanishes in DF2; it remains in DF3. Eliminating
h requires an invertible axis Hessian. That elimination never removes v(x).
Freezing a omits the u b/a^2 contribution and changes the derivative.

[V] Under (u,a,lambda)->(c u,c a,c lambda), c>0, s stays exactly fixed.
F1--F3 vanish along this family when they vanish at c=1; the current scales
by c. Consequently the direction (u,a,lambda,0) is a null direction of the
first three equations, while DF4 on it equals I_target. A nonzero prescribed
current removes this amplitude null direction. This proves neither a lower
singular-value bound nor absence of another kernel/complement instability.

## Why the energy-only argument fails

[V] In a disk B_r(x) inside Omega, choose a radial function that is one for
rho<=epsilon, log(r/rho)/log(r/epsilon) for epsilon<rho<r, and zero outside
B_r. It belongs to H1_0 and has unweighted gradient energy
2 pi/log(r/epsilon), which tends to zero, while its value at x stays one.
The weighted energy is at most 1/R_min times this expression. Smooth
approximations retain the conclusion. Point evaluation is therefore unbounded
on the two-dimensional energy space, even when the axis point is fixed.

[V] A Jacobian coercivity constant obtained for an absolute-flux law cannot
control the normalized moving-axis system in this space without an additional
argument. A sampled FE axis gradient, a finite-dimensional eigenvalue, or a
frozen-source H(div) residual does not supply that argument.

## Conditional local existence and error theorem

[V] Consider a Banach space X embedded continuously in C2(Omega), with the
declared homogeneous trace, and a data space Y for which the coupled map
F:X x R^4 -> Y x R^4 is continuously differentiable on a closed ball about U0.
The norms include explicit scales for a, lambda and x. Require a!=0 throughout
the ball, a fixed axis neighborhood and all profile/regularity premises.
Suppose A=DF(U0) has a bounded two-sided inverse and proved constants satisfy

    ||A^-1 F(U0)|| <= eta,
    ||A^-1 (DF(U)-DF(V))|| <= gamma ||U-V||,
    2 gamma eta < 1.

For gamma>0 set r=2 eta/(1+sqrt(1-2 gamma eta)); for gamma=0 set r=eta.
Require the entire ball of radius r to lie in the declared domain.
Then T(U)=U-A^-1 F(U) obeys

    ||T(U)-U0|| <= eta+gamma r^2/2 = r,
    ||DT(U)|| <= gamma r < 1.

The integral form of Taylor's theorem proves the first inequality and the
derivative bound proves contraction. The Banach fixed-point theorem gives a
unique zero in this ball and ||U*-U0||<=r. A proved embedding constant
||v||_E<=C_EX ||v||_X gives the energy bound C_EX r. Uniqueness outside this
ball is not claimed. These are analytic premises, not values supplied by a
converged Newton solve.

## Premises still open for K6

[H] A suitable continuous space and elliptic regularity estimate on the actual
piecewise curved domain must be established. Schauder estimates on a smooth
domain cannot simply be applied across unexamined boundary corners. A
piecewise C0 finite-element candidate also needs a proved reconstruction into
X and a reconstruction-error budget before this theorem can use it.

[H] The full constrained inverse, including its infinite-dimensional
complement and the axis Hessian, needs an outward bound. The finite Newton
matrix and removal of the amplitude direction establish neither.

[H] The source builder must enclose every profile piece intersected by the
whole-cell normalized flux range, including knot crossings and infinite
extensions. `law_valid` accepts small floating value jumps; it proves neither
exact continuity nor the differentiability/Lipschitz premises above. Source
remainder and derivative bounds must describe the actual represented law.

[H] The mapped residual must use the actual cell maps, verify positive R and
Jacobian throughout cells, and reconcile affine/mapped shared-edge traces.
`gs_cell_data` currently selects affine maps on unflagged cells, whereas
`gs_equilibrated_majorant` supplies all high-order geometry nodes to recovery.
Replacing the represented geometry without an error budget changes the model.

[V] `gs_error_control.f90` refuses normalized laws and `gs_certificate.f90`
restricts the outward residual upgrade to the constant-profile path. Those
guards remain necessary. Numerical, model, input/reconstruction and stochastic
uncertainties retain separate ledgers. No production nonlinear certificate is
implemented by this conditional theorem.
