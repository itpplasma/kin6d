# Physical space and computational charts

Status: geometry contract and staged options. Cartesian Euclidean physical
space is fundamental. The initial toroidal Fourier chart is an optimization,
not a requirement that every domain have a cylindrical angle or flux surface.

## Physical contract

Physics consumes Cartesian physical position, momentum, fields, force and
boundary normals on a domain Omega in R3. Computational coordinates xi map
through x=X(xi), with a checked Jacobian, inverse where required, orientation,
metric, measures and boundary map. Derive operators and transforms from that
map; hand-coded coordinate-specific numerical leaves are forbidden.

Record det(DX), its sign, and the positive volume measure |det(DX)| dxi.
Local invertibility alone does not prove global injectivity or nonoverlap.
Certification includes element validity, map/inverse errors, adjacency and
boundary representation. An imported CAD/STL wall is a geometry candidate;
its approximation error is included or explicitly outside the model claim.

The cylindrical chart X(R,Z,phi)=(R cos phi,R sin phi,Z) is useful for
ordinary toroidal machines. Coordinate ordering matters for orientation:
(R,phi,Z) and (R,Z,phi) have opposite signed Jacobians, while both have
positive volume measure R. Keep operator orientation and curl conventions
consistent with the declared ordering. Cylindrical-axis/polar chart
singularities need patches or regularity treatment; Cartesian physical space
has no magnetic-axis or X-point coordinate singularity.

Reserve a small field/domain evaluator first, with affine Cartesian and
cylindrical analytic cases. Add mapping operations only when the next case
requires them. Do not implement arbitrary topology, CAD import, a mesh
registry, or every coordinate family before the first executable slice.

## Chart choices

| Computational structure | Trigger |
|---|---|
| Cartesian affine elements | Slab verification, current sheets, space-plasma boxes |
| Cylindrical chart plus Fourier | First toroidal production case; default fusion optimization |
| Periodic cylinder | Helical/Bessel benchmark and cylinder-axis treatment |
| Spherical chart or patched spherical shell | A solar-wind/planetary case requiring it |
| Curve-following or G-frame map | Nonplanar/Figure-8 geometry benchmark |
| General curvilinear/unstructured physical mesh | Next strongly shaped boundary, wall, or non-toroidal domain |
| Flux coordinates | Optional nested-surface operator/diagnostic acceleration |
| Field-line/FCI map | Optional parallel-operator acceleration, independent of chart |

Periodic Fourier expansion is available in a valid computational periodic
coordinate, including a generalized curve coordinate. Its metric and field
coefficients can couple modes. Such a coordinate is not automatically the
cylindrical angle, and its zero mode does not generally mean axisymmetry.
Local elements in the periodic direction remain an option when Fourier
representation or nonperiodic boundaries are inappropriate.

Prefer mapped tensor patches where practical; use fitted/cut/unstructured
pieces when required rather than imposing tetrahedra everywhere. FortFEM owns
generic maps/spaces/traces; kin6d supplies the physical field and boundary law.
[Spatial numerics](numerics.md) defines FEEC/IGA/DG and batching choices.

## Magnetic topology and equilibrium representations

Flux coordinates are permitted only as an optional representation when the
state actually supports them. A VMEC/DESC/GVEC-like nested-surface state may be
used for import, initialization, compression, preconditioning or a cheaper
equilibrium model, then converted to the common physical-space field contract.
No downstream MHD/kinetic component may infer that nested surfaces exist merely
because the source file used flux coordinates.

A general three-dimensional equilibrium mesh must remain meaningful through
islands, stochastic/ergodic regions, separatrices and wall intersections. This
is a representational requirement, not permission to violate force balance.
For static smooth ideal MHD, B·grad p = 0; pressure must therefore be compatible
with the field-line topology, including flattening on ergodic components where
appropriate. Relaxed-region or resistive equilibria are different physical
models, not numerical workarounds for this constraint.

Different stationary solutions on the same external data are represented as
different states with their own provenance. Geometry stores maps/topology
diagnostics needed by a state; it does not own an inference mixture over those
states.

## Momentum frames and vector components

Physical momentum and the physical-space chart are independent choices.
Collision algebra acts on Euclidean momentum distributions, with any
computational momentum frame explicitly transformed. Cartesian components
are the common field/force API. Discretizations may store Cartesian,
covariant/contravariant, or co-rotating components when beneficial, provided
their transformations and geometric derivatives are defined and generated.

The initial toroidal FK expansion retains the co-rotating cylindrical
momentum frame from [the design](../DESIGN.md#cylindrical-momentum-frame),
so n=0 correctly represents axisymmetric scalar distributions. Frame
connection terms remain necessary. A globally Cartesian momentum harmonic
basis does not make axisymmetry n=0 at fixed Cartesian momentum.

Likewise an axisymmetric magnetic vector can have phi-dependent Cartesian
components. Axisymmetry is rotational equivariance, not constant Cartesian
components. The optimized zero-mode field coefficients may use cylindrical
components and transform to Cartesian at the physics boundary. Do not
truncate Cartesian vector Fourier modes blindly. Compare Cartesian and
cylindrical representations of the same physical orbit/field independently.

A locally B-aligned angular frame remains a later optimization with explicit
connection terms. A Cartesian physical API removes derivatives of the
physical Cartesian basis, but not derivatives of a chosen moving
computational frame or coordinate map.

## Field-line maps

An optional map follows dx/ds=B/|B| in physical Cartesian space. A
cross-sectional map at fixed toroidal increment is valid only where the
selected field component permits it; it can fail near a zero component.
Use arc length, alternate sections or another method rather than assuming
every field line crosses each toroidal plane. B=0 also invalidates b itself
and reduced magnetized models; physical coordinates do not remove that
singularity.

Maps can terminate at a wall, with certified intersection location, normal,
incidence, interpolation and boundary action. Islands, stochastic/open lines,
separatrices and reconnection need no global nested-surface assumption.
Preserve map validity only over its actual field/time state; rebuilding a
map for evolving fields has a consistency defect that enters the certificate.

FCI is one parallel-operator implementation for DK/GK, anisotropic fluid and
MHD terms. FK follows full physical characteristics in arbitrary directions
and does not require field-line coordinates. An FCI stencil also needs
measures, transverse interpolation, boundary handling and adjoint/conservation
properties; nearest traced values alone do not define a certified operator.
[XGC meshing](https://xgc.pppl.gov/html/meshing_tutorial.html) illustrates
physical cylindrical meshing through regions where flux coordinates fail.

## Figure-8 and geometry gates

Figure-8 is a mandatory planned geometry stress test alongside conventional
tokamaks, W7-X and LHD. No core physics component may require a global
cylindrical toroidal angle. An analytic curve-following case comes before a
frozen public Figure-8 equilibrium, followed by the same geometry/small-case/
production gate ladder as applicable. Do not certify an arbitrary schematic
curve map merely because it looks like a Figure-8.

[The Figure-8 study](https://arxiv.org/abs/2411.16411),
[GVEC's general frame](https://arxiv.org/abs/2410.17595), and
[boundary-conforming map construction](https://arxiv.org/abs/2411.04683)
are proposed source families for selecting and reviewing the map. A
curve-plus-transverse-frame expression is a candidate embedding; injectivity,
periodicity and nondegenerate metric require independent evidence.

Geometry tests cover affine maps, cylindrical transformations, circular and
helical fields, large-aspect-ratio limits, period maps, shaped boundaries,
and later Figure-8 and open-field wall maps. Their bands enclose field/map
observables, not an uncomputed plasma solution. Unsupported charts stay
planned rather than forcing exotic domains into the toroidal storage layout.

Chris&AI
