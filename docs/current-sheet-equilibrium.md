# Ideal equilibria with current sheets at rational surfaces (class b)

Status: design, stage 0 (2026-10-28). Nothing on this page is implemented
unless a module or test is named. This page owns the formulation of ROADMAP
R12 item 6 (weak ideal MHD, class (b) of the
[solution-space contract](equilibria-inference.md#solution-space-contract))
in the frozen-in map lane. The plasma energy and its discretization are owned
by [ideal-map-3d.md](ideal-map-3d.md); the plasma–vacuum interface by
[free-boundary-theory.md](free-boundary-theory.md); the linear per-n response
with shielding sheets by [mhd-linear-response.md](mhd-linear-response.md).

## 1. Problem

Three-dimensional ideal equilibria with nested surfaces and sheared rotational
transform generally have no smooth solution. At a rational surface
ι = n/m the resonant normal field must vanish (ideal shielding). This needs a
surface current K = n × [[B]], so B is tangentially discontinuous there. Force
balance across the surface is [[p + B²/2]] = 0 (μ0 = 1 on this page). The
Pfirsch–Schlüter current ∝ p'/(mι − n) diverges unless p' = 0 at the rational
or the surface carries a pressure jump. A smooth discretization cannot
represent such a state: under refinement its residual or error stalls at a
floor, or decreases at a fractional rate set by the singularity rather than
by the element order.

Goal: a formulation whose discrete solutions converge at the design order of
the element to a well-defined weak (class (b)) solution, for tokamaks with
RMPs and for stellarators, plus a quantified relation to the classical ideal
shielding limit.

## 2. The frozen-in map and its interfaces

The map lane minimizes

    W[X] = ∫ ( |DX b|²/(2J) − p J ) ds,    B = DX b/J,   J = det DX,

over maps X from reference coordinates s into space, with a frozen closed
2-form b and frozen pressure p(s) (γ = 0, Kruskal–Kulsrud). Facts used below,
all exact consequences of this definition:

1. **Labels are material.** A reference surface S = {s1 = s1*} is mapped to
   X(S). The toroidal and poloidal fluxes enclosed by every label surface equal
   their reference values for every admissible map. Hence ι(s1) is fixed by b,
   and so is the label s1* of every rational surface. Placing an interface at
   a rational label is therefore preserved exactly by ideal motion; no
   tracking or remeshing is needed. (This is the owner's key observation and
   it is correct.)
2. **B·n = 0 on interfaces is built in** whenever b is tangent to the label
   surfaces (b¹ = 0). It holds on both sides, for every map.
3. **B ∈ H(div) for C⁰ maps.** If X is continuous across S and b·n is
   continuous (here zero), B is divergence free in the sense of distributions:
   face fluxes are reference fluxes.
4. **Tangential discontinuity.** If X is only C⁰ across S (broken normal
   derivative), DX b/J differs on the two sides, so B can jump tangentially.
   The surface current K = n × [[B]] is then represented without smoothing.
5. **Force balance on S is the natural interface condition.** Varying X by
   δX supported near S and integrating by parts on each side gives
   δW = ∫ (−∇p + J×B)·δX dV + ∮_S [[p + B²/2]] (n·δX) dA, using B·n = 0 on
   S. Stationarity enforces [[p + B²/2]] = 0 weakly on S, exactly as the
   plasma boundary map equations do at the plasma–vacuum interface
   ([free-boundary-theory.md](free-boundary-theory.md)). Nothing has to be
   added to the energy; the pressure jump of a sheet is allowed and balanced.

The interface-aware current is J = curl B on each side plus K δ_S, which is
the class (b) "J = curl B in an interface-aware sense".

## 3. Evaluation of the key idea: continuous ι does not work

The idea "interface at the rational label, C⁰ map, smooth elsewhere" fails
when the frozen ι is continuous through the rational. The linear ideal
response then has a **step** in the normal displacement at the rational
(Hahm & Kulsrud 1985 and §5.1 below): ξ_n → ±C on the two sides. A step is
not a C⁰ map, and the C⁰ constraint is not closed in the energy: near the
resonance the bending energy is ∫ F² ξ'² with F = k·B ≈ F' x. A layer of width
ε that joins the two values costs ∫₀^ε F'²x² (C/ε)² dx = F'² C² ε/3 → 0. The
infimum over C⁰ maps equals the infimum over step functions and is not
attained. Discrete C⁰ solutions form a layer of width ~h, and their error
decreases like h^(1/2) independently of the element order and of interface
alignment.

Measured (scratch prototype of the linear slab problem of §5.1, P2 elements,
L2 relative errors of ξ and of the field (ψ, ψ'/k) against the oracle):

| ι at rational | mesh | 8 | 32 | 128 | 256 | rate (field) |
|---|---|---|---|---|---|---|
| continuous | node at the rational | 1.1e-2 | 4.1e-3 | 2.0e-3 | 1.4e-3 | 0.5 |
| continuous | rational inside an element | 1.1e-2 | 3.3e-3 | 1.6e-3 | 1.2e-3 | 0.5 |
| jump Δ+ = 0.2, Δ− = 0.05 | node at the interface | 1.0e-2 | 1.6e-3 | 1.4e-4 | 3.7e-5 | 1.9 (ξ: 2.9) |
| same | interface inside an element | 1.4e-2 | 7.1e-3 | 3.4e-3 | 2.3e-3 | 0.5 (ξ: 1.1) |

A broken (double-valued) normal displacement at the rational would represent
the linear step exactly, but it opens a gap or an overlap of plasma at finite
amplitude; it is admissible for the linear response only (that is what the
existing per-n response does with its shielding sheets).

## 4. Formulation adopted: interfaces with a frozen ι jump

Loizu, Hudson, Bhattacharjee, Lazerson and Helander
([Phys. Plasmas 22, 090704, 2015](https://doi.org/10.1063/1.4931094)) showed
that ideal equilibria with continuously nested surfaces exist, perturbatively
and nonlinearly (SPEC), when ι is **discontinuous** across the resonant
surface. In the map lane this is a property of the frozen data:

- The interface label s1 = σ carries a jump of the frozen ι,
  ι(σ−) < n/m < ι(σ+) (or reversed). No surface is exactly resonant. The jump
  is an equilibrium sheet current K₀ = n × [[B₀]] already in the unperturbed
  state, balanced by [[p + B₀²/2]] = 0 (a pressure jump is allowed).
- Away from interfaces b, p are smooth; at the interface b may jump in its
  tangential components (b¹ = 0 on both sides, so div b = 0 holds).
- X is C⁰ across σ and smooth elsewhere. With mesh faces on σ, the exact
  solution is piecewise smooth on the element partition and Lagrange elements
  of order k converge at their design order (k + 1 for the map in L2, k for B).
  The prototype rows above confirm this.

The jump size Δι is a model input, like the pressure profile. Two limits make
it physically interpretable:

- **Classical shielding limit.** As Δι → 0 the solution outside a layer of
  width O(Δι) converges to the singular continuous-ι ideal solution, and the
  sheet current converges to the classical ideal shielding current (§5;
  exactly independent of Δ in the symmetric slab, O(Δ) otherwise).
- **Existence margin.** Nonlinearly the layer has |∂ξ_n/∂s1| ≈ C/(Δ a), so an
  invertible map needs C ≲ Δ a: the admissible perturbation amplitude scales
  linearly with the jump (estimate from §5.1, to be measured in stage 2).
  Beyond it the map folds and no nested equilibrium with that jump exists.

**Pressure.** With Δι > 0, p may be smooth through the interface: the
Pfirsch–Schlüter current is finite, of order p'/(mΔι), because mι − n is
bounded away from zero on both sides. To approach the Δι → 0 limit with
bounded currents use p' = 0 in a neighbourhood of the interface, or put the
pressure gradient into a pressure jump at the interface balanced by [[B²]].
Both are frozen-data choices; the solver is unchanged.

**What is carried by the map lane today.** The existing P_k continuous
Lagrange map is already C⁰ with broken derivatives across every element face.
Required changes are small and local:

1. reference meshes whose faces lie on the interface labels (for the 2D
   triangle × Fourier lane: an internal chain on the curve χ = χ_σ, as the
   free-boundary lane already does for Γ̂_p);
2. frozen data b, p evaluated per region (cell tag), never by thresholding a
   label value, so that the jump is sharp under curved-geometry approximation;
3. a sheet-current diagnostic: traces of B from both sides of interface
   faces, K = n × [[B]], resolved in the Fourier harmonics of the rational;
4. gauges: the transversal direction ∇χ/|∇χ| is continuous across χ_σ; the
   ω (toroidal relabelling) moments are taken per region.

## 5. Analytic references (stage 1)

### 5.1 Slab (Hahm & Kulsrud 1985 geometry, with jump)

Plasma −a < x < a between walls, z invariant, B₀ = (0, B_y(x), B_z(x)),
B_y = B0 (x/a + Δ± ), Δ+ ≥ 0 for x > 0 and −Δ− ≤ 0 for x < 0 (Δ± = 0:
Hahm–Kulsrud, iota continuous with the resonance k·B = 0 at x = 0). Any B_z(x),
p(x) with p + (B_y² + B_z²)/2 piecewise constant and [[p + B²/2]] = 0 is an
equilibrium. Walls x = ±a + δ± cos ky. With γ = 0 (or incompressible) the
marginal ideal displacement ξ_x = ξ(x) cos ky satisfies Newcomb's slab
equation (F²ξ'/k²)' = F²ξ, F = k B_y, independent of B_z and p. With
ψ = Fξ and F'' = 0 on each side, ψ'' = k²ψ. Interface conditions at x = 0:
ξ continuous (material interface) and continuity of the Lagrangian total
pressure perturbation, which here is [[F²ξ']] = 0. Then

    ψ+(x) = [F+ ξ0 sinh k(a − x) + ψa+ sinh kx]/sinh ka,
    ψ−(x) = [F− ξ0 sinh k(a + x) − ψa− sinh kx]/sinh ka,

with F± = ±k B0 Δ±, ψa± = F(±a) δ±, and ξ0 from

    F+(ψ+'(0) − F'ξ0) = F−(ψ−'(0) − F'ξ0),   F' = k B0/a,

a scalar linear equation. The Lagrangian field perturbation is
(−k B_y ξ sin ky, −B_y ξ' cos ky, 0) (it is incompressible, B1z = 0), and the
sheet-current perturbation is K1_z = [[B_Ly]] = −[[B_y ξ']] cos ky.

- Hahm–Kulsrud (Δ± = 0): ψ(0) = 0, ψ = ψa sinh(k|x|)/sinh(ka) for symmetric
  walls δ+ = −δ− = −δ, ξ steps by ±a ψa'(0)/(k B0) at 0, and
  K1 = 2 k B0 δ/sinh(ka) (cos ky amplitude).
- Symmetric jump Δ+ = Δ− = Δ, symmetric walls: ξ0 = 0, ψ(0) = 0 on both sides,
  K1 = 2 k B0 (1 + Δ) δ/sinh(ka) — the Hahm–Kulsrud sheet current (with the
  wall value of B_y), independent of Δ; ξ is continuous and joins the step
  over a layer of width Δ a.
- General walls: ξ0 → k a (δ+ + δ−)/(2 sinh ka) as Δ → 0, ψ(0±) = O(Δ), so
  the shielding condition ψ(0) = 0 is recovered at rate O(Δ).

### 5.2 Screw pinch with a resonant wall perturbation

Annulus r0 < r < b, periodic length 2πR0, a current-free plasma
B_θ = c±/r, B_z = Bz± constant on each side of the interface r = σ, p
constant on each side, [[c²/σ² + Bz² + 2p]] = 0. Perturbation phase
α = mθ + kz, F = m B_θ/r + k B_z, q = r B_z/(R0 B_θ). Inner wall fixed,
outer wall r = b + δ cos α. Static ideal response with γ = 0 in a
current-free, constant-pressure region is a vacuum field B1 = ∇χ,
χ = χ̂(r) sin α, χ̂ = A I_m(|k| r) + B K_m(|k| r), and ξ_r = ξ̂ cos α with
χ̂' = −F ξ̂. Four conditions fix (A±, B±):

    χ̂−'(r0) = 0,   χ̂+'(b) = −F+(b) δ,
    χ̂−'(σ)/F−(σ) = χ̂+'(σ)/F+(σ) = −ξ̂σ,
    [[F χ̂ − ξ̂σ B_θ²/r]](σ) = 0    (Lagrangian total pressure).

Sheet current: K = r̂ × [[B_L]] with [[B_Lθ]] = [[(m/r) χ̂ − ξ̂ B_θ/r]],
[[B_Lz]] = [[k χ̂]] (cos α amplitudes). Continuous ι (F(σ) = 0, ideal
shielding χ̂'(σ) = 0): the inner region is fully shielded (χ̂− ≡ 0), the outer
χ̂+ has χ̂+'(σ) = 0, χ̂+'(b) = −F(b)δ, the Eulerian sheet current is
((m/σ) χ̂+(σ), k χ̂+(σ)) and ξ̂ steps from 0 to −χ̂+''(σ)/F'(σ). With a jump,
χ̂±'(σ) = O(Δ) and the sheet current converges to the shielding value at
rate O(Δ). I_m, K_m are evaluated by series and recurrence (a documented
hand-written exception: FortSym has no Bessel functions), checked by the
Wronskian I_m K_m' − I_m' K_m = −1/x.

### 5.3 Published nonlinear results

Loizu et al. ([Phys. Plasmas 22, 022501, 2015](https://doi.org/10.1063/1.4906888):
slab Hahm–Kulsrud with SPEC, ideal limit, δ-function and 1/x currents;
[090704](https://doi.org/10.1063/1.4931094): nonlinear existence with
discontinuous ι) are published literature. Under the
[source policy](literature-policy.md#benchmark-sources) they are context for a
matched-model comparison, not an oracle. SPEC solves MRxMHD with a finite
number of Taylor regions; only its ideal limit (many interfaces, μ → ideal)
matches this model. Their tables were not accessible from the development
environment (no network route to the publishers or arXiv), so no published
number is transcribed here yet; the stage-1 tests use only §5.1–5.2.

## 6. Composition

- **Fixed and free boundary.** Interior interfaces and the plasma boundary
  are the same kind of object: a material label surface with B·n = 0 and the
  natural traction [[p + B²/2]] = 0. The free-boundary coupling
  (shell potential, Steklov map, coils) is untouched; RMP coils enter as the
  n ≠ 0 exterior drive. A tokamak RMP case puts interfaces at the q = m/n
  labels with frozen jumps.
- **Linear per-n response.** In the small-amplitude limit the map-lane sheet
  current must agree with the shielding sheets of
  [mhd-linear-response.md](mhd-linear-response.md) as Δι → 0 (rate O(Δι),
  §5). This is the stage-3 cross-check.
- **FEEC non-nested lane.** Its Eulerian A ∈ H(curl) on a fixed mesh cannot
  hold a moving sheet. It does not represent class (b) states; a map-lane
  class (b) state can be its initial guess or comparator, but NN results are
  class (c) or regularized states and are labelled so.
- **Island seeding (PR #125) and resistive lanes** change topology; they are
  a different model class and are not compared with class (b) except in
  limits stated per case.

## 7. Proven, assumed, open

Proven on this page (exact or elementary): items 1–5 of §2; the C⁰
non-closure argument of §3; the slab and screw-pinch linear solutions and
their Δ → 0 limits (§5.1, §5.2; derivations to be checked by the stage-1
tests).

Assumed: existence of nonlinear minimizers with Δι > 0 below the existence
margin (supported by Loizu et al., not proven here); the interface condition
for curved interfaces in 3D (it follows from item 5, but its discrete
consistency is to be measured).

Open, with the known obstruction:

- **Many rationals in 3D.** A stellarator boundary drives every rational in
  the ι range. Only rationals resolved by the Fourier set (|m| ≤ M) can be
  forced, but their number grows like M². Interfaces are placed at the
  low-order rationals with significant resonant forcing; higher ones need
  p' = 0 locally, or the h-convergence holds at fixed M only. Convergence
  in M remains open.
- **Close rationals.** Two interfaces closer than their layers Δι-width
  interact; the element size between them bounds the resolution.
- **Choice of Δι.** It is a model parameter. Reporting results at several Δι
  and the extrapolation Δι → 0 with the derived rate is the honest output.
- **Islands.** The map preserves topology by construction; island formation
  is outside class (b) (no assumption is made, it is a model property).

## 8. Verification plan

| Stage | Case | Evidence |
|---|---|---|
| 1 | slab §5.1, screw pinch §5.2 | oracle module, FortSym pointwise leaves, interface conditions and limits as tests |
| 2 | slab, nonlinear frozen-in map, interface node at σ | design-order convergence of ξ, B and K1 against §5.1 (linear limit by δ-differences), self-convergence at finite δ, h^(1/2) or fold with continuous ι, h^(1/2) with the jump inside an element |
| 3 | tokamak + RMP | sheet currents at q = m/n vs the per-n response as Δι → 0 |
| 4 | stellarator with one low-order rational | small-resolution convergence; SPEC/SPECTRE ideal-limit comparison only where assumptions match |
| 5 | Newton-MR on the symmetric map Hessian | iteration counts and identical end points |

Chris&AI
