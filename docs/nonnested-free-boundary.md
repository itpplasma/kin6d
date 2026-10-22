# Coil-driven NN2 boundary coupling

[V] `coil_set_read_makegrid` reads explicit MAKEGRID filaments in metres with
signed effective filament currents in A-turns. The closing zero-current row
supplies the group ID and name. The reader checks closure, constant current,
finite values and the final `end` record, and rejects implicit mirror rules.
The `periods` header is metadata and does not replicate geometry. Only synthetic
filaments are used in the tests.

[V] `coil_w7x_normalized_currents` supplies EIM, DMM, FTM, SVM, WNJ and YFK
ratios in the order NPC1–5, PCA, PCB from Davies et al.,
[Table 2, arXiv:2607.14722v1](https://arxiv.org/html/2607.14722v1#S4.T2).
These are filament-current ratios. No winding factor is applied. A caller
must supply the group-to-circuit mapping, orientation and dimensional scale
for its geometry. These data alone do not define an operating configuration.

[V] `nn_attach_exterior` uses [FortFEM](https://github.com/lazy-fortran/fortfem)'s
full three-dimensional continuous-P2/discontinuous-P1 scalar Laplace BEM on
fixed Gamma_c. It sweeps the NN2 mesh's outer polynomial edges through a full
turn, so FEM and BEM use the same physical interface. The kernel and all panel
quadrature remain in FortFEM. There is no axisymmetric DtN reduction.

## Variational coupling

[V] Let C map NN2's real Fourier vector-potential coefficients a to magnetic
two-form coefficients. Let T a contain the moments of n.curl A against the
BEM trace basis, and let g contain the same moments of the direct coil field.
The harmonic reference F_vac/R e_phi has zero normal trace on this swept
surface. For the exterior scalar potential u, B_response = -grad u and
S u = T a - g, where S is FortFEM's positive exterior energy map. Thus the
normal-response energy is

    E_ext = (T a - g)^T S^(-1) (T a - g) / (2 mu0).

[V] The interior vacuum functional is the squared difference from the coil
field, int_Omega |B - B_coil|^2/(2 mu0). Its volume derivative is Ampere's
law because curl B_coil = 0 in the coil-free interior. Its boundary term,
combined with E_ext, imposes tangential continuity of the response. The
assembled NN2 magnetic residual gains

    [T^T S^(-1) T a - T^T S^(-1) g - C^T <V2,B_coil-B_h>] / mu0.

[V] Full-turn BEM integrals are divided by 2*pi once to match NN2's phi-mean
volume pairing. `nc_jacobian` includes the same quadratic boundary block.
The magnetic trace and gauge tests are released, with the constant gauge
multiplier and the curl-free toroidal period of A pinned. The pressure and
electric-potential boundary data keep the NN2 closure contract. In the exact discrete vacuum subspace (zero source,
zero pressure, electric potential and loop voltage), those homogeneous blocks
are eliminated to preserve zero current under mixed sparse LU. Only boundary
columns enter the dense exterior Schur complement, which is stored in the existing sparse matrix container.

[V] This entry point requires I1's `I_phi` drive with zero net toroidal
plasma current and rejects an internal conductor. A scalar exterior response has zero periods,
whereas a nonzero plasma current needs a separately prescribed topological
harmonic response. The known coil field is kept separate and can carry
nonzero periods. Coil clearance is a caller responsibility.

## Evidence and limits

[H] `nonnested_exterior_3d_vacuum_pilot` exercises a synthetic toroidal coil
set and a remote non-axisymmetric coil through NN2, with a Jacobian check,
a direct Biot–Savart field comparison and a perturbed coefficient seed.
`nn_exterior_trace` reports the algebraic defect of S u = T a - g. It returns
FortFEM's physical outward derivative separately from the positive energy load.
An algebraic defect is not a continuous interface-error estimate.

[H] G8 additionally requires field and trace errors <= 1e-3, an independently
measured interface defect <= 1e-6, a virtual-casing comparison, and return of
a geometrically displaced seed within 1e-3 minor radii. The pilot does not
establish those gates. In particular, coefficient perturbations are not a
measurement of geometric displacement, and the current module has no exterior
point-field evaluator or continuous error certificate. Generic quadrature and
representation gaps belong in FortFEM. This work must remain draft until G8
passes on a qualified dependency revision.
