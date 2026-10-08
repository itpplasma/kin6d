# Axisymmetric Grad–Shafranov derivation and exact oracles

- Scope: static scalar-pressure MHD equilibrium; no flow, pressure anisotropy or 3D equilibrium solver.
- Assume smooth nested surfaces, R>0, axisymmetry, and F=F(psi), p=p(psi). Axis and boundary limits require separate regularity checks.
- SI units; physical phi increases counter-clockwise viewed from +Z. Cylindrical basis (e_R,e_phi,e_Z) is right-handed.
- psi is signed poloidal flux per radian. Full poloidal flux is 2*pi*psi; full toroidal flux is the surface integral of B_phi dR dZ.
- Convention: B=grad(psi) cross grad(phi)+F grad(phi). Therefore B_R=-psi_Z/R, B_Z=psi_R/R, B_phi=F/R.
- Canonical theta is counter-clockwise in (R-R0,Z); (r,theta,phi) is left-handed. This is COCOS 3; native solver conventions need explicit reversible maps.

## Force balance

From cylindrical curl and axisymmetry,

\[
\mu_0 j_R=-F'\psi_Z/R=F'B_R,\qquad
\mu_0 j_Z=F'\psi_R/R=F'B_Z,
\]
\[
\mu_0 j_\phi=-\Delta^*\psi/R,\qquad
\Delta^*\psi=\psi_{RR}-\psi_R/R+\psi_{ZZ}.
\]

Consequently,

\[
\mathbf j\times\mathbf B
=-\frac{\Delta^*\psi+FF'}{\mu_0R^2}\nabla\psi.
\]

With grad(p)=p' grad(psi), equilibrium requires

\[
\boxed{\Delta^*\psi+\mu_0R^2p'(\psi)+F(\psi)F'(\psi)=0.}
\]

- div(B)=(-psi_ZR+psi_RZ)/R=0 by commuting mixed derivatives.
- The equation holds away from critical points of psi; smooth continuation supplies the axis value. Dividing by grad(psi) at the axis is invalid.
- Forward inputs: boundary flux trace, p'(psi), FF'(psi), one F value and pressure offset. Current and q are outputs.
- Integrate F^2(psi)=F_edge^2+2*integral_edge^psi FF'(u)du; reject F^2<=0. Integrate pressure with an explicit pressure offset.
- Prescribing q as well generally changes the unknown F/current profile. A source amplitude named q_scale is not an exact constant-q constraint.

For continuous P2 flux, the classical element current is
`j_phi=-(psi_RR+psi_ZZ-psi_R/R)/(mu0*R)`. Across an interior face with unit
normal pointing left to right, the additional distributional sheet is
`K_phi=-normal·(grad(psi)_right-grad(psi)_left)/(mu0*R)`. Polynomial degree
does not remove these face terms or turn a vertex cusp into a smooth axis.

At a regular stationary minimum with positive Hessian H, local elliptic coarea
gives `q_axis=F_axis/(R_axis*sqrt(det(H)))` on the positive-outward-flux branch.
This formula requires a smooth quadratic neighbourhood. The native P2 axis
routine qualifies only strict element-interior minima; it does not assign this
limit to vertex/edge minima or claim a continuum error bound. A direct contour
integration on an independently parameterized ellipse checks the axis formula
and both the physical field and F signs in `test_gs_p2`.

## Polynomial GS oracle

Define C>0 and

\[
\Psi=C[(R^2-R_0^2)^2+4R^2Z^2],\qquad
\Delta^*\Psi=16CR^2.
\]

- Set p'=-16C/mu0 and F constant. With p_edge=0, p=p'*(Psi-Psi_edge) is nonnegative inside Psi_edge.
- An arbitrary rectangular box uses the exact nonconstant Dirichlet trace. It tests the PDE, not an arbitrary physical plasma/vacuum boundary.
- Near the axis, Psi=4CR0^2*((R-R0)^2+Z^2)+higher orders; q_axis=F/(8CR0^3).
- Choose C=F/(8*q_axis*R0^3). The full q profile remains a derived output.

The exact closed surface Psi_edge=4CR0^2*a^2 has parameterization

\[
R(t)=\sqrt{R_0^2+2R_0a\cos t},\qquad
Z(t)=R_0a\sin t/R(t).
\]

- Require 2a<R0. This finite-aspect-ratio LCFS is not an exact circle.
- Interior surfaces replace a by a*sqrt(s_pol), where s_pol=Psi/Psi_edge.
- For general surface radius r, |partial(R,Z)/partial(r,t)|=R0^2*r/R^2. Its orientation is positive in (R,Z).
- Field-line dphi/dt=F/(8C*R^3), since B_pol points along the positive t tangent.
- With angular mean over t, q(s_pol)=F/(8C)*mean(R^-3).
- The full toroidal flux is Phi(s_pol)=2*pi*F*R0^2*a^2*s_pol*mean(sin(t)^2/R^3).
- Obtain the flux formula by Green's theorem: Phi=F*integral(log(R)*dZ)=F*integral(-Z*R_t/R*dt), then substitute the surface map.
- Equivalently dPhi/dPsi=2*pi*q. This identity independently checks flux quadrature and profile maps.
- VMEC++ input uses p(s_tor), iota(s_tor)=1/q and full Phi_edge, with s_tor=Phi/Phi_edge. Refine boundary Fourier and radial interpolation separately.

## Logarithmic zero-beta GS oracle

For D>0 and R>0,

\[
\Psi=D\left[R^2\log(R/R_0)-\frac{R^2-R_0^2}{2}+Z^2\right],
\qquad \Delta^*\Psi=4D.
\]

- Set p=0 and FF'=-4D. F^2=F_axis^2-8D*Psi; restrict the domain to F^2>0.
- Near the axis Psi=D*((R-R0)^2+Z^2)+higher orders, hence q_axis=F_axis/(2DR0).
- D=F_axis/(2*q_axis*R0) sets this axis limit, not the full q profile.
- Exact rectangular traces provide a second source/field oracle with a nonconstant F profile.

## Strict cylindrical zero-beta limit

For **right-handed** straight cylindrical coordinates (r,theta_RH,z), define

\[
B_r=0,\quad B_z=\frac{B_0}{1+k^2r^2},\quad B_\theta=krB_z.
\]

Then

\[
(\nabla\times B)_\theta=-\partial_rB_z,
\quad(\nabla\times B)_z=r^{-1}\partial_r(rB_\theta),
\quad\nabla\times B=\frac{2k}{1+k^2r^2}B.
\]

- This Gold–Hoyle field is force-free: j cross B=0 and p=0. A uniform axial field plus B_theta proportional to r is not force-free.
- Declare axial period 2*pi*R0. Native q_RH=r*B_z/(R0*B_theta_RH)=1/(k*R0), with its regular axis limit understood. The native C0 value +1.5 is not canonical toroidal q.
- Native cylindrical flux analogue per radian: psi_RH=R0*B0/(2k)*log(1+k^2*r^2); dpsi_RH/dr=R0*B_theta_RH.
- This is an exact straight-cylinder oracle. Embedding it in a finite-radius torus does not preserve exact force balance.
- Fixed circular toroidal LCFS cases must solve GS; their interior surfaces need not be concentric circles.

### Orientation when comparing C0 with the canonical toroidal limit

- Align the straight axial direction with physical +phi; do not reverse the physical field. At fixed phi, e_r=cos(theta_c)e_R+sin(theta_c)e_Z and e_theta_c=-sin(theta_c)e_R+cos(theta_c)e_Z. Their cross product is -e_phi, so (e_r,e_theta_c,e_phi) is left-handed.
- A right-handed cylinder with e_z=e_phi therefore has theta_c=-theta_RH and e_theta_RH=-e_theta_c. Consequently B_theta_c=-B_theta_RH, B_phi=B_z and F=R0*B_z in the strict cylindrical limit.
- The canonical field formula gives B_theta_c=(1/R0)*dpsi_c/dr. Hence psi_c=-psi_RH, up to its declared gauge, and q_c=r*B_z/(R0*B_theta_c)=-q_RH. Native +1.5 maps to canonical -1.5 for this axial alignment.
- This is a coordinate/flux-analogue map in the straight limit, not a finite-R toroidal force-balanced construction. Toroidal curvature and R dependence must be solved separately.
- The physical helicity coefficient alpha=2*k/(1+k^2*r^2) in curl(B)=alpha*B is unchanged. Using right-handed cylindrical curl formulas with the left-handed canonical components would introduce an artificial sign error.
- Evidence: the basis/cross-product and flux-orientation argument above is a hand derivation. `test/check_cylinder_orientation.py` independently checks Cartesian bases, field pullback, both q/k and B0 polarities, and a Cartesian finite-difference curl; it is a numerical diagnostic, not an additional native CAS identity. The eight registered identities retain their original RH cylinder scope.

## Checked scope and evidence

- [Variational bridge and profile constraints](grad-shafranov-variational-bridge.md) derives the strong/weak/source-functional and nested ideal-MHD relationships. Its separate native gate checks fourteen local/source identities and two rejection controls; it does not prove boundary theorems or complete third-party discretizations.

- Native FortSym `check_identity` verifies eight exact identities: polynomial and logarithmic Delta_star, axisymmetric divergence/toroidal curl, Gold–Hoyle theta/z curl, force balance and constant q.
- These use symbolic differentiation and exact native zero decisions; no randomized probes. They are algebra checks under R>0, R0>0 and nonzero denominators, not Lean proofs or compiler verification.
- Persistent source: `test/test_gs_symbolic.f90`; registered CTest `grad_shafranov_symbolic_identities` passes 8/8 native FortSym identities. Configure the optional pinned FortSym generator build to register this gate.
- Retained generator and registered symbolic gate execute FortSym `a6142712af998b727c245fb6cfbb3745b296ec24`, KIN6D's pinned source. The initial scratch algebra used `5be9575b707f7d2dff998c238046b223cb5f3a19`; it is superseded by the retained eight-identity gate. Fo driver SHA256 `3510671d89c610075d759861f010e70952783eb19daacc25e294924240ccdd3e`. CMake/CTest are build/test authorities.
- Weak-form and energy identities are hand derivations; the eight CAS checks do not prove integration by parts, boundary conditions or minimization.
- Surface Jacobian, Green-theorem flux map and axis limits above are explicit hand derivations, not among the eight CAS identities. A separate midpoint-quadrature/finite-difference smoke check at s_pol=.001,.1,.5,1 verifies the surface trace and Jacobian and dPhi/dPsi=2*pi*q (maximum relative q discrepancy 2.37e-11). Retained source: `test/check_gs_surface_map.py`; run with Python/NumPy. This smoke check is not a registered native gate or proof; quadrature/limit convergence remains an admission gate.
- Numerical validation additionally needs mesh, boundary and interpolation convergence; common off-grid B/p/j/force-balance observables; and independent total-current/full-flux/q reconstruction.

- [Distributional P1 current](grad-shafranov-p1-current.md): independent cell/jump/boundary integrals and fifteen exact antiderivatives; native Release5/5 and checked Debug4/4 focused gates pass. Pointwise smooth current, force and axis-q remain unadmitted.

## Primary references

- [Cerfon and Freidberg, analytical Grad–Shafranov equilibria (2010)](https://doi.org/10.1063/1.3328818): Solov'ev-type analytical families; the explicit oracle above is checked directly.
- [Pataki et al., high-order Grad–Shafranov solver (2013)](https://arxiv.org/abs/1210.2113): independent analytical-solution comparison and numerical convergence.
- [Vandas and Romashets, uniform-twist flux ropes (2017)](https://doi.org/10.1051/0004-6361/201731412): cylindrical Gold–Hoyle field and its toroidal limitations.
