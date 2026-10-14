# Direct spot calibration and held TC24 qualification

The response calibration path reads a derivative jump directly from sampled
resonant fluxes on both sides of the GPEC spot. It uses no Frobenius coefficient,
Mercier-power transport, Delta-prime fit, limiting intercept or code-agreement
factor. All scientific outputs are candidates; continuous numerical, model and
input enclosures remain unavailable.

`kin6d_response_calibration` constructs separate derivative stencils at
`psi_s +/- spot/(abs(n)*abs(dq/dpsi_N))`. Each stencil must bracket its spot and
stay on its side of the rational. FortNum supplies the local polynomial
**derivative weights**. The jump is outer derivative minus inner derivative.
At fixed samples its JVP is the same linear action; `spot_jump_vjp` is its
Hermitian storage transpose. Sample-coordinate, equilibrium and shear
derivatives are unavailable through this wrapper.

For a unit outward-normal sheet `N x grad Lambda`, measured with the same
magnetic operator, angles and spot stencil as the response, the equivalent
shielding field is

```
b = -unit_self_flux * response_jump / (area * unit_jump).
```

The minus sign, peak-amplitude convention and area are those of
`kin6d_mhd_response_n`. Units are T for b, A for equivalent Lambda,
T m2/rad2 for flux, and m2/rad for ring area. The normalization primal and
real JVP/VJP come from `gen/gen_response_calibration.f90` and its committed
FortSym output. Zero/nonfinite unit jumps, zero unit self flux and nonpositive
area are refused. This unit-sheet normalization is not a calibration to GPEC.

The independent cylinder oracle has continuous normal B and a unit jump of
magnetic scalar potential. For `k=abs(n)/R0`, `z=k*r_s`, the self flux is
`2*mu0*R0*z^2*I_m'(z)*K_m'(z)`. Its derivative at either spot follows from the
modified-Bessel equation. The test evaluates I from its series and K from its
integral, independently of FEEC/BEM. Another test reads an actual resolved
single-ring R-B1 strip and checks its spot jump and integrated current against
the closed ideal step-profile solution, retaining a complex drive.

## TC24 magnetic calibration executable

`gs_response_calibration --config FILE [--validate]` reads `/spot_calibration/`:
`eqdsk`, `n`, `qabs`, `nrad`, `nth_max`, `order`, `nfine`, `spot`,
`stencil_fraction`. `--validate` checks the retained equilibrium, target
rational, shear and spot coverage before a matrix is assembled. The solve
places all discovered rationals and both five-point stencils on the mesh,
solves one unit sheet and writes its direct spot jump, self flux, area, samples
and retained imported edge pressure to `output/gs_response_calibration.nc`.

This is a magnetic calibration executable. It does not solve a full-pressure
coupled response, enforce interface traction or supply a cut-support model.
Its time/RSS measurements must therefore retain that narrower scope. Refine
sample spacing, radial and poloidal controls before using its unit-sheet
normalization in a qualified coupled response.

## Frozen qualification harness

`response_qualify --config qualification.nml` is the Fortran authority for the
predeclared TC24 gates in the held root's `RULE.md`. It reads `/qualification/`
with up to 288 records: three cuts, eight axes, three levels, two phasings and
two target rationals. Integer identities are cut (c0995,c0997,full), axis
(radial,core,grading,harmonic,poloidal,extent,resistivity,frequency), level,
phase (000,010), and target (3,10/3). Missing records remain unavailable;
duplicate/out-of-range identities are errors.

Each record carries complex b, vacuum b, total resonant flux and layer current;
true residual and balance defect; actual inner/outer extents; and explicit
availability, coupled, calibrated, full-pressure, all-rational, independent
interface-test and whole-edge-strip status. A single-ring result may be
compared diagnostically but cannot pass coupled qualification. Status flags
are evidence supplied by the solver, not certificates established by this
postprocessor. Preserve the diagnostic evidence that supports them.

All ladder gates use complex differences, including the vacuum-scale absolute
floor. Five mesh axes require a final change <=5% and non-increasing changes,
except that both changes <=0.5% suffice. Extent requires pairwise spread <=5%;
identical capped distances require independent passed interface evidence or a
whole-edge-strip solve for each repeated pair. S and gamma require a final
change <=1% and non-increasing changes. Every member needs true residual <=1e-8
and balance <=1e-6, full pressure, calibration and all rational sidebands.

GPEC amplitude ratio and phase remain diagnostic on failed rows. A code
verdict requires all gates and independently admitted MARS normalization;
both/neither code within 10% gives no verdict. If more than six of twelve rows
fail, all verdicts are suppressed and the lane reports `STOP_NO_TIEBREAK`.
The NetCDF output preserves input, provenance, gate status and ladder changes.

The optional `scripts/prepare_response_qualification.py ROOT --output FILE`
serializes the existing frozen `ladder.json` and six raw GPEC references to
this namelist. It executes no solver or scheduler. Legacy
`runs/<case>/<case>_res.dat` gives diagnostic single-ring b and vacuum b;
`runs/<case>/qualification.json` can supply coupled records after integration.
That JSON must declare unchanged `controls`, `cut`, `spot`, and records with
`phase`, `q` and the fields above; complex values are `[real,imaginary]`.
The adapter retains MARS estimates as unadmitted, never silently normalizes
or fits them. Direct namelist use needs no Python.

`--calibration-decks DIR` emits unit-sheet validation decks for the frozen
outer mesh controls and both targets. These are not coupled solver decks.
The held registry proposal and production root must remain blocked until the
coupled executable/decks, TC24 resources and every frozen physics ladder are
actually available. The controller alone registers and launches the root.
