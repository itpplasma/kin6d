# Distributional current of a continuous P1 flux

- Scope: axisymmetric physical `(R,Z)`, `R>0`; CCW physical phi and signed flux per radian. Synthetic validation only: no equilibrium solve, smooth Hessian, continuous-force error estimate or axis-q admission.
- Native field correspondence: `src/generated/kin6d_grad_shafranov_leaves.f90::gs_field` gives `BR=-psi_Z/R`, `BZ=psi_R/R`, `Bphi=F/R`. The test calls this native leaf and checks its signs/units against independent affine fields.
- Numerical test: `test/test_gs_p1_weak_current.f90`. Optional native FortSym exact check: `test/test_gs_p1_current_identities.f90`; fifteen antiderivative identities, no probes.

## Cell, edge and boundary identity

Let continuous piecewise-affine `psi_h` have constant gradient `g_K` in each triangle K. Put `w=grad(psi_h)/R`. The poloidal field has continuous normal trace across an edge, but its tangential trace may jump.

1. Cylindrical curl gives `mu0*j_phi = -div(w)`. Within a cell, `div(w)=-g_K,R/R^2`, so `mu0*j_phi,cell=g_K,R/R^2`. A zero affine Hessian does **not** imply zero current.
2. On an interior edge choose normal n from K to L. The distributional divergence contributes `[(g_L-g_K)·n]/R` times the edge delta measure. Current contributes its negative.
3. Integrate by parts in every cell, then sum. For continuous test v:

```text
mu0 <j_phi,h,v> =
    sum_K integral_K (g_K,R/R^2) v dR dZ
  - sum_interior_edges integral_edge ((g_L-g_K)·n/R) v ds
  = integral_domain grad(psi_h)·grad(v)/R dR dZ
  - integral_boundary v grad(psi_h)·n/R ds.
```

- Changing edge orientation also exchanges its cell labels; the jump measure is invariant. Triangle winding and edge parametrization do not change the identity.
- For zero-boundary v the last term vanishes. For constant v the weak volume term vanishes; discarding the boundary term would give the wrong total current.
- This is a current **distribution**, including current sheets. It is not a pointwise smooth `j×B` assertion at an edge: B itself has a trace jump there. It does not provide a derivative-error bound for an equilibrium candidate.

## Independent exact fixture

- Domain `[1,2]×[0,1]`; vertices A=(1,0), B=(2,0), C=(2,1), D=(1,1); triangles ABC and ACD, shared edge AC.
- Affine fixture: nodal `psi=R+2Z`, both cell gradients `(1,2)`, no jump.
- Non-affine fixture: nodal `psi=RZ`, P1 interpolants `psi_ABC=2Z`, `psi_ACD=R+Z-1`; gradients `(0,2)` and `(1,1)`.
- Normal ABC→ACD is `(-1,1)/sqrt(2)`. Parameterize AC by `(R,Z)=(1+t,t)`, `0<=t<=1`. Its normalized current measure is `2 dt/(1+t)`.
- Test weights: `1`, `R`, `Z`, `RZ`, and `b=(R-1)(2-R)Z(1-Z)`. The last vanishes on every exterior edge.
- Write `L=log(2)`. C and E below denote mu0 times the cell and interior-edge current integrals; W is weak volume and D outward boundary flux; **C+E=W-D**.

| Affine weight | C | E | W | D |
|---|---|---|---|---|
| 1 | 1/2 | 0 | 0 | -1/2 |
| R | L | 0 | L | 0 |
| Z | 1/4 | 0 | 2L | 2L-1/4 |
| RZ | L/2 | 0 | 2+L/2 | 2 |
| b | (3L-2)/6 | 0 | (3L-2)/6 | 0 |

| Non-affine weight | C | E | W | D |
|---|---|---|---|---|
| 1 | 1-L | 2L | 0 | -1-L |
| R | 2L-1 | 2 | 2L-1 | -2 |
| Z | L-1/2 | 2-2L | 1 | L-1/2 |
| RZ | 1/4 | 1 | 7/4 | 1/2 |
| b | (25-36L)/6 | 8L-11/2 | 2L-4/3 | 0 |

The tables follow independently by integrating each polynomial weight in Z over `[0,R-1]` or `[R-1,1]`, then in R, and by evaluating the four exterior edge traces. For example:

```text
C_nonaffine[1] = integral_1^2 (2-R)/R^2 dR = 1-L
E_nonaffine[b] = 2 integral_0^1 t^2(1-t)^2/(1+t) dt = 8L-11/2
C_nonaffine[b] = integral_1^2
    (-2R^3+15R^2-43R+58-36/R+8/R^2)/6 dR = (25-36L)/6.
```

- Native FortSym differentiates all five cell antiderivatives for each fixture and all five edge antiderivatives. It requires exact zero decisions; log arguments are positive on the stated integration intervals.
- The numerical oracle separately integrates cell and edge/boundary terms using FortNum 32-point Gauss–Legendre rules, rather than evaluating these antiderivatives. Synthetic test formulas are deliberate independent oracles, not duplicated production physics.
- Both triangle winding reversal and cell-owner permutation are tested. Reversed line parametrization retains outward geometric normals.
- Registered negative controls omit edge current, omit the cylindrical cell current, or reverse current sign. These must fail the named identity. Their nonzero exits are expected CTest outcomes; the positive registered test must also pass.

## Admission limits

- This gate closes local signed-field/current-identity checks for these synthetic P1 fixtures. It does not certify the native FEM assembly, nonlinear branch, mesh/source quadrature error, current recovery, smooth force or an actual device current.
- Raw P1 has no classical axis Hessian. A smooth recovered axis-q limit needs a declared recovery and separate mesh/derivative/axis error budget; polynomial reproduction alone is insufficient.
- Existing native weak residuals, exact Solovev field/flux refinement and recovered Ampere controls remain separate evidence. No P1 equilibrium state or evaluator is changed by this test.
