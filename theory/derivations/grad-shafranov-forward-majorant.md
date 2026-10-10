# Weighted forward Grad–Shafranov energy majorant

This derivation proves a continuous forward error bound under explicit premises.
The accompanying native FortSym test checks local algebra and negative premise
witnesses. It does not implement an estimator, validate an arbitrary mesh, or
certify an unknown equilibrium or inverse-q solution.

## Declared operator and domain

Let `Omega` be a bounded Lipschitz domain in physical `(R,Z)`, contained in a
rectangle of side lengths `L_R,L_Z`, with `0<r0<=R<=r1`. The candidate and exact
solution must have the same specified Dirichlet trace on the entire boundary.
Set `V=H0^1(Omega)`, `w=1/R`, and

\[
 a(u,\phi)=\int_\Omega w\nabla u\cdot\nabla\phi,
 \qquad \|\phi\|_a^2=a(\phi,\phi),
 \qquad \|\phi\|_w^2=\int_\Omega w\phi^2.
\]

The executed forward weak equation is

\[
 a(u,\phi)=\int_\Omega wS(x,u)\phi,\qquad
 S(x,t)=\mu_0R^2p'(t)+H(t),\quad H=FF',\quad\phi\in V.       \tag{1}
\]

Its strong sign is `-Delta_star u=mu0*R^2*pprime+FFprime`. The code's load
density is `w*S=mu0*R*pprime+FFprime/R`. Physical profile knots are fixed inputs
in this theorem. The same argument allows a spatial manufactured source `S`.

## Explicit stability constant

Zero extension to the containing rectangle and its sine eigenfunctions give

\[
 C_P^2=\frac1{\pi^2(L_R^{-2}+L_Z^{-2})},\quad
 \|\phi\|_{L^2}\le C_P\|\nabla\phi\|_{L^2},\quad
 \|\phi\|_w\le C_w\|\phi\|_a,\quad C_w^2=\frac{r_1}{r_0}C_P^2. \tag{2}
\]

These are conservative upper constants, not fitted polygon eigenvalues.
Outward bounds on the containing rectangle and radii are required for a
numerical guarantee. Using the rational lower bound `pi>3` remains valid.

Suppose `S` is globally Lipschitz in its scalar argument and has upper
one-sided slope `m`:

\[
 [S(x,s)-S(x,t)](s-t)\le m(s-t)^2,\qquad
 m_+=\max(0,m),\qquad \theta=1-m_+C_w^2>0.              \tag{3}
\]

For the continuous piecewise-linear profiles in `gs_profile_value`, a global
choice is the maximum of zero and every segment slope
`mu0*R^2*diff(pprime)/diff(flux)+diff(FFprime)/diff(flux)`, maximized at the
appropriate endpoint of `[r0,r1]`. Zero covers the constant endpoint extensions.
These are the retained linear-profile rules, not the current cubic-law rules.
Linear profiles use common flux knots. If separate knots are admitted later,
take their sorted union and both slopes on every resulting segment. Large
negative slopes do not harm this one-sided condition.

For the operator `A(u)=a(u,.)-(w*S(u),.)`, (2)–(3) imply

\[
 \langle A(u)-A(v),u-v\rangle\ge\theta\|u-v\|_a^2.         \tag{4}
\]

Indeed the lower estimate is `||e||_a^2-m*||e||_w^2`, which is at least
`||e||_a^2-m_+*||e||_w^2 >= theta*||e||_a^2`. A negative source slope
cannot improve this energy coercivity constant using an upper Poincare bound.

Continuity, linear growth and this strong monotonicity give a unique weak
solution with the declared trace. A fixed source has `m=0,theta=1`.
Failure of the sufficient condition `theta>0` is inconclusive: a sharper
certified spectral bound may still establish stability. An iteration's
convergence does not establish it.

## Functional majorant and proof

Let `v` be a conforming `H1` candidate with the specified trace, and let
`y` be any globally conforming vector field in `H(div,Omega)`. Set

\[
 d=\|w^{-1/2}(y-w\nabla v)\|_{L^2},\qquad
 b=\|w^{-1/2}(wS(x,v)+\operatorname{div}y)\|_{L^2}.
\]

Under (2)–(3), the guaranteed bound is

\[
 \boxed{\|u-v\|_a\le\frac{d+C_wb}{\theta}}
 \le\frac{\sqrt{(1+\beta)d^2+(1+1/\beta)C_w^2b^2}}{\theta},
 \qquad\beta>0.                                        \tag{5}
\]

Write `rho_v(phi)=(w*S(v),phi)-a(v,phi)` and `e=u-v`.
Subtracting (1) and using (4) gives
`theta*||e||_a^2 <= rho_v(e)`. Global flux conformity and zero test trace give

\[
 \rho_v(\phi)=(wS(v)+\operatorname{div}y,\phi)
                  +(y-w\nabla v,\nabla\phi).
\]

Weighted Cauchy–Schwarz and (2) bound this by
`(d+C_w*b)*||phi||_a`. Division and Young's inequality prove (5).
No Galerkin orthogonality or mesh convergence rate is used. Exact flux
equilibration is unnecessary; even `y=0` is valid, although often loose.
An RT/BDM reconstruction can improve sharpness, but a broken flux with omitted
normal jumps is not in `H(div)` and cannot be substituted into (5).

For fixed source and matching trace, `rho_v=a(u-v,.)`; hence
`||u-v||_a=||rho_v||_(a*)` by Riesz representation. Restriction to an FE test
space gives a lower bound on this continuous dual norm. The native relative
Euclidean coefficient residual is neither that norm nor a continuous upper
error bound.

## Residual jumps and explicit fallback

On a valid nonoverlapping affine triangulation, define

\[
 r_T=wS(v)+\operatorname{div}(w\nabla v)|_T,\qquad
 j_E=(w\nabla v|_{T_+})\cdot n_++(w\nabla v|_{T_-})\cdot n_-.
\]

Elementwise Green identities give
`rho_v(phi)=sum_T(r_T,phi)-sum_E(j_E,phi)` for zero-boundary tests.
P1 has a `-v_R/R^2` volume divergence and possible interface sheets even
though its classical Hessian is zero. P2 also retains interface flux jumps.

A coarse bound can retain all jumps with explicit constants. For each interior
edge select one adjacent triangle, its opposite-vertex altitude `h_E`, and
the largest distance `D_E` from that vertex to any triangle point.
Applying the divergence theorem to `(x-x_opposite)*phi^2` yields

\[
 \|\phi\|_E^2\le\frac2{h_E}
 (\|\phi\|_T^2+D_E\|\phi\|_T\|\nabla\phi\|_T).
\]

Global norm bounds from (2) then give

\[
 C_E^2=\frac{2r_1}{h_E}(C_P^2+D_E C_P),\qquad
 \|u-v\|_a\le\theta^{-1}\left[
  \sqrt{r_1}C_P\|r\|_{L^2}+\sum_E C_E\|j_E\|_{L^2(E)}\right]. \tag{6}
\]

This fallback is proved but can be pessimistic. Standard adaptive indicators
with local mesh-size factors require their own reliability/interpolation
constants before becoming guaranteed upper bounds.

## Data, integration and geometry premises

- For an approximate source `S_tilde`, add
  `C_w*||sqrt(w)*(S(v)-S_tilde(v))||` to (5)'s numerator. The target source's
  stability condition must still hold.
- All norm integrals, geometry predicates, input conversions and reported
  endpoints need verified arithmetic. Paired quadrature differences are
  estimates; they are not guaranteed remainder bounds. Direct enclosure of
  the continuous residual of stored `v` already includes its algebraic and
  quadrature-induced state error; do not add it twice.
- If traces differ, construct an explicit `H1` lifting `ell`, apply (5) to
  `v+ell`, and add `||ell||_a`. Maximum boundary flux error alone does not
  bound this lifting energy.
- A curved domain requires a validated geometry map. Pullback diffusion is
  `abs(det(DF))*w(F)*DF^(-1)*DF^(-T)`. Coefficient/source discrepancy,
  transformed trace and map conditioning are separate consistency bounds.
  A polygon norm does not cover missing physical volume.
- Exact Galerkin integration and absolute global source Lipschitz bound `L`
  give quasi-optimality factor `(1+L*C_w^2)/theta`. Here `L` is distinct from
  the one-sided `m`. A P2 energy `h^2` rate also needs appropriate regularity
  and mesh constants; no pointwise Hessian/current rate follows automatically.

## Executed curved-P2/P3 functional and cubic laws

`gs_functional_majorant` minimizes and evaluates a slightly different valid
functional: its balance defect uses the **unweighted** physical-area norm.
For `g=mu0*R*pprime+FFprime/R`, set

\[
 C_F=\sqrt{r_1}C_P,\quad
 D=\|\sqrt R(y-\nabla v/R)\|_{L^2},\quad
 B=C_F\|g(v)+\operatorname{div}y\|_{L^2}.
\]

Since `||phi||_L2 <= C_P*sqrt(r1)*||phi||_a`, the same Green identity proves
`||u-v||_a <= (D+B)/theta`. The executed `cf2=r1*C_P^2` is therefore correct;
it is not `C_w^2`, and its unweighted balance norm must not be identified with
the weighted `b` in (5). The recovery objective is `D^2+B^2`; the returned
numerator is `D+B`. Any conforming recovered flux is admissible, so an
incomplete recovery solve affects sharpness rather than the continuous
inequality. Ordinary order-eight Duffy quadrature makes the executed value an
estimate of this functional, not a guaranteed quadrature upper bound.

The same proof applies on the represented curved domain when the P2/P3 maps
are globally one-to-one, orientation preserving, nonoverlapping, and share
their edge traces. The continuous vector Lagrange recovery then belongs to
`H1(Omega)^2`, hence `H(div,Omega)`. Physical gradients and divergence must use
the same inverse map as the field. Bernstein geometry controls enclose every
mapped point and supply the containing rectangle used by the code. Positive
Jacobians at quadrature nodes alone do not prove global map admissibility.
Comparison with a different physical LCFS additionally needs the geometry and
trace consistency terms above; it cannot inherit a represented-domain bound.

For the executed cubic law, on a coordinate interval write
`P(t)=a0+a1*t+a2*t^2+a3*t^3`, where `t=x-knot_i`. Its exact primitive has
derivative `P`, and its slope is `a1+2*a2*t+3*a3*t^2`. For fixed physical knots,
the one-sided source slope at each radius is the quadratic

\[
 Q_R(t)=\mu_0R^2 P_p'(t)+P_H'(t).
\]

On a **declared finite solution-value interval**, its maximum is computable
by evaluating interval endpoints and the stationary point of this quadratic
when that point lies inside the interval. First maximize at `R=r0,r1`:
`Q_R(t)` is affine in `R^2`, so these endpoints suffice even when the pressure
slope changes sign. Split at all knots and include end-segment extrapolation.
Taking the maximum of zero and the resulting upper slope supplies (3) for
pairs of states whose values stay in that interval; an absolute maximum
similarly supplies the Lipschitz `L`. The piecewise cubic must be exactly
continuous at every knot: integrating its segment derivatives then bounds all
secants, including those crossing knots. Reader acceptance within a floating
continuity tolerance does not imply this premise. For a law with small jumps,
use an exactly continuous target law and the approximate-source consistency
term above rather than assigning a finite secant bound from segment slopes.
This conditional estimate alone does not establish existence or that an unknown
solution stays in the declared interval. The cubic reader extrapolates its
first/last polynomial; it does **not** use constant endpoint extensions, and a
cubic cannot generally be assigned a finite global slope from its knot values.

For normalized laws with a **fixed nonzero** axis `a`, `x=1-psi/a`. Multiply
the source derivative by `-1/a` before taking its maximum; the law API also
admits positive axes, for which multiplying an already computed upper slope
by this negative factor would reverse the inequality. The physical primitives
carry the factor
`-a`, which cancels this coordinate derivative: `dp/dpsi=pprime(x)` and
`d(F^2/2)/dpsi=FFprime(x)`. `gen_grad_shafranov` proves these cubic derivative
and normalized-primitive identities with FortSym before emitting the kernels.
If the axis is recomputed with the state, the additional term is
`delta x=-delta psi/a+psi*delta a/a^2`. That nonlocal term is absent from the
scalar fixed-axis monotonicity argument.

The API reports status 4 for a frozen-source diagnostic without supplied
stability. Status 5 uses a numerical discrete positive-part power iteration;
neither an iterate's Rayleigh quotient nor its small eigen-residual bounds
the largest continuous eigenvalue. A caller-supplied positive stability factor
is a conditional premise, not something the API proves. Measured effectivity
against an exact solution tests sharpness; a finite CHEASE comparison also
contains reference error. Neither supplies the missing nonlinear premise.

## Current, field and inverse limits

An energy bound `E` implies a distributional toroidal-current error bound
`||j_phi(v)-j_phi(u)||_(a*) <= E/mu0` for zero-boundary scalar tests.
It includes interface sheets. It is not a pointwise current or physical
normal-trace bound. Poloidal field error in toroidal physical volume satisfies
`||Bpol(v)-Bpol(u)||_(L2(V))^2=2*pi*||u-v||_a^2 <= 2*pi*E^2`.

Toroidal-field/pressure, Maxwell stress and q/coarea need additional profile,
product, admissibility and nonvanishing-gradient bounds. Inverse q/Phi and
moving normalized-axis profile knots require a separate coupled continuous
inverse/inf-sup bound. Axis point evaluation is not bounded on general 2D `H1`;
a discrete shape vector does not supply that bound.

In the executed prescribed-q path the exact coarea identities remain
`q=F*K/(2*pi)` and `Phi_edge=int_Omega F/R dR dZ`. In a regular star-shaped
chart `(rho,theta)`, the delta-function change of variables gives
`K=int rho/(R*partial_rho psi) dtheta`. For a curved physical map `X` of an
affine chart `A`, the area-density ratio is
`det(DX)/det(DA)*R_chart/R_physical`, precisely the native cell-arc factor.
This identity requires the same regular contour in both maps and a strictly
positive outward radial derivative (otherwise the denominator is its absolute
value); it gives no contour-error bound near a critical
point. The independent q probes and area-versus-contour Phi check are
consistency/refinement evidence, not bounds on the coupled inverse problem.

A coupled bound must control the weak GS residual **and** the q-moment/Phi
residuals with a branch-specific inverse bound for their joint derivative.
It must include source coefficients, solved F_edge, moving axis, and contour
events in a norm where those operations are continuous. The fixed-source
energy norm alone does not do this. In particular, zero moment residuals do
not control the off-moment q profile; its independent defect remains necessary.
No such continuous coupled inverse constant is supplied by the current code
or the E1/E2/TC24 finite-reference effectivity measurements.

## Native exact controls

`test_gs_forward_majorant_identities` proves 19 identities through the pinned
native FortSym engine with no probe fallback. Three native nonzero decisions
and three expected-failure CTest fixtures reject omitted jump, boundary and
quadrature premises. They do not exercise a production estimator API.

The fixed-source oracle on `[1,2]x[0,1]` is
`u=-R^4+5R^2-4`, `mu0=1,pprime=8,FFprime=0`, with its specified nonzero trace.
Set `v=u+alpha*(R-1)*(2-R)*Z*(1-Z)` and `y=grad(u)/R`.
Then `div(y)=-8R`, and the majorant in (5) equals the energy error:

\[
 \|u-v\|_a^2=\alpha^2\left(\frac{49}{30}\log2-\frac{67}{60}\right).
\]

Exact polynomial antiderivatives verify this value. The candidate is conforming
quartic, not a claim of native P2 representability. A separate manufactured
spatial datum plus `lambda*psi` verifies the nonlinear residual and second
variation sign, without claiming it is a physical KIN profile.

Negative witnesses are a tent's nonzero normal jump `4/3`, harmonic `u=Z`
versus zero candidate with mismatched trace, and the positive polynomial
`(t-1/2)^2*(t^2-t+1/6)^2`. One- and two-node Gauss rules both give zero for
that polynomial, while its exact integral is `11/15120`.

FortSym establishes algebra, not the Sobolev inequalities or mesh/trace
premises. The analytic proof supplies that bridge. A numerical unknown-state
certificate remains unavailable until those premises and integrals are checked.

## Inverse-produced Solovev error relative to the analytic target

This exact-case theorem bounds the error to the declared analytic Solovev state,
not to an unidentified solution of the interpolated prescribed-q input. It does
not establish coupled inverse reliability for other cases. Let `U` be the
polynomial target restricted to the valid represented domain `Omega_h`, `v` its
conforming zero-trace candidate, and `y` any conforming `H(div)` recovery. For
`U=psi_scale*((R/R0)^4-2*(R/R0)^2+4*(R*Z/R0^2)^2+c)`, its exact load is
`g_*=-16*psi_scale*R/R0^4`. This definition also accounts for rounded physical
source inputs. With the executed unweighted balance convention, set

\[
 D=\|\sqrt R(y-\nabla v/R)\|_2,\quad
 B=C_F\|g_v+\operatorname{div}y\|_2,\quad
 C=C_F\|g_*-g_v\|_2.                                      \tag{12}
\]

The target generally has nonzero trace on the approximate curved boundary.
For any explicit continuous lifting `ell` of that trace,

\[
 \boxed{\|U-v\|_a\le\sqrt{(D+B+C)^2+\|\ell\|_a^2}}.       \tag{13}
\]

To prove this, let `h` be the harmonic extension of `U`'s boundary trace.
Then `e0=U-v-h` has zero trace, and its residual is the fixed-source residual
for `g_*`, because `a(h,phi)=0`. The functional identity gives
`||e0||_a<=D+B+C`; orthogonality gives
`||U-v||_a^2=||e0||_a^2+||h||_a^2`. The Dirichlet minimum principle implies
`||h||_a<=||ell||_a`. No axis point estimate or unknown stability constant is
used. Replacing `B+C` by `C_F||g_*+div y||_2` retains source cancellation.

An explicit lifting uses only the actual P3 geometry and analytic polynomial.
Start with continuous reference vertex hats carrying `U` at boundary vertices
and zero at other vertices. Each boundary-edge trace `b(t)=U(F_edge(t))` has
degree at most 12. Subtract its endpoint linear interpolant to obtain `c(t)`.
On the adjacent reference triangle extend it as
`(lambda_i+lambda_j)^12*c(lambda_j/(lambda_i+lambda_j))` and by zero elsewhere.
It is polynomial, has the required edge trace, and vanishes on the other two
edges. Summing these edge terms and the vertex hats gives the exact continuous
lifting, including rounded nonzero vertex traces.

Its mapped energy integrand is rational, not polynomial:
`|adj(J)^T grad(ell_hat)|^2/(R*det J)`. For each cell let `r_T,d_T` be positive
lower Bernstein controls for `R,det J`, `M_ab` upper absolute controls for the
Jacobian components, and `g_x,g_y` upper absolute controls for reference lifting
derivatives. Its energy contribution is bounded by

\[
 \frac{(M_{22}g_x+M_{21}g_y)^2+(M_{12}g_x+M_{11}g_y)^2}
      {2r_Td_T}.                                           \tag{14}
\]

The determinant is degree 4, the lifting degree 12, and its derivatives degree 11.
For constant candidate pprime, with `D_T,R_T` upper determinant/radial controls
and `f_T` the maximum absolute recovered FFprime cubic over the entire flux
range enclosed by candidate controls, a source correction upper bound is

\[
 C^2\le C_F^2\sum_T\frac{D_T}{2}
 \left(\left|-16\,psi_{scale}/R_0^4-\mu_0p'_v\right|R_T
                         +f_T/r_T\right)^2.                \tag{15}
\]

Cubic endpoint and stationary values determine `f_T`, including the executed
polynomial extrapolation if the candidate range exceeds its knots. Equations
(13)--(15) retain represented-domain geometry and source mismatch separately.
They cover every represented cell, not additional exact-domain volume outside
`Omega_h`. Sound numerical upper evaluation of the norms and controls remains
necessary to call a computed value guaranteed; ordinary paired quadrature and
floating control evaluations are estimates, not arithmetic enclosures. The
unmodified status-4 estimator remains a general frozen-source diagnostic.

## Primary references

- [Repin, A posteriori analysis for PDEs](https://math.aalto.fi/opetus/lme/ApostPDERepin.pdf): functional majorants and monotonicity.
- [Valdman (2009)](https://onlinelibrary.wiley.com/doi/10.1155/2009/164519): conforming approximation and `H(div)` majorant minimization.
- [Ern–Vohralík, unified framework](https://who.rocq.inria.fr/Martin.Vohralik/Files/Pub_proc_FVCA_11_3.pdf): flux conformity, explicit trace inequalities and guaranteed reconstruction bounds.
