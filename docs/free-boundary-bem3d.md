# Full-surface exterior coupling

[V] `kin6d_equilibrium_boundary_3d` supplies the Cartesian R2 boundary adapter.
It assembles the curved continuous-P2/discontinuous-P1 scalar Laplace BEM in
[FortFEM](https://github.com/lazy-fortran/fortfem), factors its flux block once,
and constructs the consumer trace matrix `R^T S R`. The restriction maps
consumer coefficients to BEM nodes; the load uses its transpose.

[V] `gamma_c_3d_energy` evaluates, with `u = trace - drive`,

    E = (u^T S u/2 - h^T trace)/mu0,
    dE/dtrace = (S u - h)/mu0.

[V] `gamma_c_3d_harmonics` additionally constructs decaying exterior fields
from interior current cycles, checks every prescribed boundary period, and
adds `u^T C c + c^T G c/2` to the numerator. The scalar load gains `C c/mu0`.
The coordinates are physical circulations. There are genus-g independent
exterior periods, not 2g independent boundary periods. A nonextendable period
pair is rejected. Total cut-flux constraints need the scalar contribution
as well; this entry point prescribes periods only.

[V] The scalar coil drive is the single-valued remainder. It must not include
a second copy of the prescribed harmonic field. The chosen harmonic
representatives and the shell's closed reference forms can differ by an exact
form: include that change of primitive in the harmonic trace bookkeeping.
Neither a multivalued angular potential nor a fitted rescaling of a field is
an admissible substitute for these data.

[V] `ideal_fb_attach_exterior_3d` replaces the exterior quadratic form in the
existing ideal-map free-boundary KKT functional. Its gradient and Hessian use
the same full matrix, including couplings between real Fourier modes. Both
the potential block solve and the nonlinear solver consume those assembled
rows. The plasma and shell shape variations retain their existing generated
energy derivatives. Full-turn surface integrals are divided by `nfp` once,
matching the plasma energy per field period.

[V] `ideal_fb_use_bem3d` builds this restriction on the fixed coupling surface
of an existing `ideal_fb_t`. Curved P2/P3 meridional edges are swept through
its chart. A caller chooses azimuthal resolution and BEM pair tolerances;
the swept geometry itself is evaluated through the chart jet. Existing
Fourier drives must be set before attaching, and are converted with their
physical full-turn pairing. Other callers can supply a general parametric
surface and restriction directly to `gamma_c_3d_build`.

[V] The shell formulation places the DtN surface on fixed Gamma_c, enclosing
the free plasma Gamma_p. The generic boundary adapter can act on Gamma_p,
but the present solver integration does not differentiate a BEM assembled
on moving Gamma_p. `ideal_fb_init` still requires the cylindrical shell
initializer. Arbitrary chart shell initialization and assembled BEM geometry
sensitivities are unavailable. Per-harmonic S0 reduction and the old
per-harmonic hypothesis report reject an attached full matrix.

[V] Native gates are `equilibrium_boundary_3d_d6` (sphere energy and trace
variation), `equilibrium_boundary_3d_harmonic_data`, and
`ideal_fb_bem3d_derivatives` (the complete coupled functional and KKT action).
The last includes an axisymmetric smooth-trace energy comparison with the
existing per-n boundary operator. Quadrature is controlled numerically;
continuous field and equilibrium error certificates remain unavailable.

[H] `ideal_fb_landreman_bem3d_axisymmetric` runs the existing epsilon=0
Landreman drive and free-boundary/GS comparison through the full-surface BEM.
Its registered accuracy gates are unchanged. This is an analytic limit of
the 3D path, not evidence for epsilon>0.

[H] The non-axisymmetric Landreman benchmark still needs a checked extension
of the virtual-casing remainder to the chosen Gamma_c (G4), its two cycle
periods/linking convention, and traces of that single-valued remainder.
The current adapter does not construct or certify that continuation. See
[the benchmark derivation](../theory/derivations/free-boundary-landreman-benchmark.md).
