# Two-body interaction laboratory

This Chapter 1 demonstration compares genuine, unscreened repulsive Coulomb
scattering with an ideal neutral hard-sphere contact. It is a microscopic
verification case, not a plasma shielding model. The animation uses the x-y
projection of a Cartesian 3D trajectory whose initial data happen to give the
invariant plane z=0. The tests also use nonplanar Cartesian initial data.

## Named model and ledger

Two particles have equal masses m1=m2=1, center of mass at rest, and positive
charge product C=1 in normalized nonrelativistic electrostatic units. With
relative position r=x1-x2 and relative velocity u, reduced mass mu=1/2 gives

    r' = u, u' = alpha r/|r|^3, alpha=C/mu=2.

The physical positions are x1=r/2 and x2=-r/2. The forces are equal and opposite,
F1=C r/|r|^3 and F2=-F1. The physical energy is
mu*(|u|^2/2+alpha/|r|). The plotted specific energy is the expression in
parentheses. There are no sources, reservoirs, imposed confining forces,
softening, screening lengths or imposed Yukawa interactions. Radiation,
relativistic corrections and physical-data uncertainty have unavailable bounds.

The candidate uses fixed-step velocity Verlet from retained FortSym definitions.
The neutral comparison follows free motion until the first approaching,
nongrazing contact at |r|=1, then reflects u elastically in the contact normal.
It is a separate ideal contact model, not a neutral plasma collision operator.
Neutral trajectory certification and contact-event sensitivities are explicitly
unavailable. Neutral outlines show physical radius 1/2; filled center markers
identify particle positions. Coulomb markers depict point particles.

## Continuous Coulomb endpoint certificate

Binary64 initial r, u, alpha and each positive step h are interpreted as exact
real parameters. The checker does not use the Rutherford orbit, energy drift,
or the candidate as its correctness premise. It starts with outward initial
point intervals and propagates boxes for the continuous ODE.

For an incoming enclosure B and a proposed step [0,h], it constructs a box Y
and verifies strict componentwise inclusion

    B + [0,h] F(Y) is contained in the interior of Y,

where F(r,u)=(u,a(r)). It also verifies the outward lower bound for |r|^2 over
Y is positive. Thus F is continuously differentiable on a neighborhood of Y.
Local existence and uniqueness hold there. The integral equation and the
strict inward margin prevent the flow starting anywhere in B from exiting Y
before h; compactness and the nonsingularity guard extend it through h.
This establishes a full step tube independently of the candidate.

For each coordinate the endpoint enclosures are then

    r(h) in R0 + h U0 + h^2 a(Yr)/2,
    u(h) in U0 + h a(R0) + h^2 Da(Yr) Yu/2.

These are the second-order integral Taylor remainders with nonnegative
integration weights. FortSym generates the force, its directional derivative,
and both interval Taylor expressions from the same symbolic source. FortNum's
declared outward-arithmetic contract is the checker arithmetic premise.

Nonfinite boxes, a singularity guard failure, or unsuccessful self inclusion
set the trajectory status to unavailable for that and every subsequent step.
The finite Verlet candidate can continue. A caller must inspect the returned
`bounded` flag before using boxes; old boxes after failure carry no current
claim. A candidate-only advance also permanently clears the flag. Reinitializing
starts a new independent initial-value problem.

The app shows boxes only at observation endpoints. Colored trajectory segments
and the energy curve are candidate polylines, with no continuous-between-frame
band claim. Force arrows are candidate forces at candidate positions; their
true-flow force bands and energy observable bands are unavailable. Each
particle-coordinate error radius is evaluated outward as
max_i sup |Ri-r_candidate_i|/2. The requested radius is 0.01 of coordinate scale
1, with absolute floor 1e-12; enclosure availability and budget satisfaction are
separate CSV statuses. Exact observation time is step count times exact h; the
CSV supplies an outward interval for it alongside the displayed float time.
The neutral candidate is evaluated at the displayed represented time.

## Derivatives, independent evidence and cost

Primal, JVP and VJP leaves share the retained FortSym source. Their pairing is
ordinary Euclidean Cartesian input/output pairing. They differentiate the
smooth force and fixed-step numerical map, including h and alpha; they do not
claim a certified continuous-flow derivative or a physical trajectory adjoint.
The nongrazing contact leaf derivatives have their stated smooth branch domain;
the application offers no complete contact-event sensitivity.

`test_interaction` compares arbitrary 3D numerical trajectories with an
independent real128 repulsive Kepler hyperbola, checks the Rutherford
impact/deflection identity, and verifies second-order refinement. It checks true
positions and velocities inside propagated boxes throughout the six-unit run,
force antisymmetry and analytical force values, a directional finite-difference
JVP, Euclidean VJP pairing, elastic contact, and explicit failure cases. The
hyperbola is test-only and is never passed to the generator or checker.

The state has fixed-size coordinate boxes and no hot-call allocation. The app
reports separately measured CPU time for candidate advancement and for the
same advancement with checker attempts, excluding sampling and rendering.
It reports their extra cost and factor; a factor of -1 means an unresolved
baseline timer. If the checker fails, subsequent attempts cease, as stated in
the timing output. These are single-case measurements, not performance leadership
claims. The final GIF is produced from the native FortPlot MP4 outside source
history.
