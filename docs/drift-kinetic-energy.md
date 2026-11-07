# Axisymmetric energy-dependent drift kinetics

`kin6d_dk_energy` solves the single-species, axisymmetric drift-kinetic equation
with Lorentz, test-particle or full linearised self-collisions. Speed uses the
collision operator's consistent Galerkin grid, including its polynomial tail;
pitch uses Legendre polynomials and the periodic poloidal angle uses Fourier
collocation. The bordered system uses pivoted band LU and iterative refinement.

`energy_axi_coefficients` returns the parallel-flow coefficients `D31` and
`D32` in the normalization defined in the module header, and the rotation
coefficient `k = 5/2 - D32/D31`. The two drives have speed weights 1 and `x^2`,
where `x` is speed divided by the species thermal speed. Implicit tangents use
`ds = dn/n - 2 dT/T`, where `n` is density and `T` is temperature.
The Coulomb logarithm, geometry, species mass and charge, and discretization
are held fixed. Thus the log-density tangent
equals the `s` tangent and the log-temperature tangent is minus twice it.
The returned residual measures the componentwise backward error of the discrete
linear solve; it does not bound the error in the continuous kinetic solution.
The certification status is `unavailable`.

`energy_axi_batch` shares the speed grid and collision model across surfaces.
It assembles one reference collision operator and rescales it for each species
state. Each surface returns its own status, so an invalid surface does not stop
the remaining profile. A supplied cached operator must have the identical grid,
including nodes, element edges, basis nodes and differentiation matrix.

`kin6d_dk_force_balance` exposes `energy_force_balance_cgs` for the single-ion
force balance used by [KAMEL](https://github.com/itpplasma/KAMEL). Its inputs use
Gaussian centimetre-gram-second units: temperature in erg, density in cm^-3,
magnetic field in gauss, charge in statcoulomb, lengths in centimetres and
velocity in cm/s. `nr` and `tr` are the density and temperature derivatives
with respect to the same physical radius. The safety factor `q` and toroidal
field component `bz` retain their signs. Radial gradients follow increasing
effective minor radius; the signed poloidal velocity follows KAMEL's formula
and coordinate convention. The outputs are radial electric field `er`
in statvolt/cm and poloidal velocity `vth` in cm/s. The optional `direction`
argument has the same fields as the input and may vary any of its twelve
components; the returned `der` and `dvth` are directional derivatives.
No smoothing or coordinate-sign
conversion is applied.

`kin6d_dk_energy_profile_io` provides `energy_profile_write_cgs` for an
in-process caller. It accepts the batch responses, their solver statuses and
the Gaussian-unit force-balance inputs above; each input's `k` is replaced by
the corresponding solver response. `p%r` contains finite, nonnegative, strictly
increasing effective radii in cm.
The caller creates `output/` and supplies a `run_t` with its captured literal
and resolved configuration, including geometry, collision model and grid.
The exporter writes `output/<run name>.nc` through the common result contract,
retaining the inputs, flow coefficients, density/temperature scaling tangents,
backward errors, force balance and separate per-surface statuses.

If every surface succeeds, `output/Er.dat` and `output/Vth.dat` contain
KAMEL's two headerless columns: effective radius in cm, then electric field
in statV/cm or poloidal velocity in cm/s. These legacy paths are shared across
run names; concurrent exports use separate working directories. A failed solve
or force balance
returns `PROFILE_INCOMPLETE`, removes previous legacy files and retains the
successful rows in NetCDF; unavailable values are NaN. Invalid array sizes or
radius ordering return `PROFILE_BAD_INPUT` before changing files. Legacy-file
I/O errors return `PROFILE_IO_ERROR`; the canonical result records the export
status. `PROFILE_OK` records export success and is not a numerical certificate.

The local rotation-ratio and force-balance values, tangents and adjoints come from the
retained definition in `gen/gen_dk_energy.f90`. The force-balance value test
uses a synthetic input evaluated independently with
[KAMEL's profile equations](https://github.com/itpplasma/KAMEL/blob/266210c5f030dae287c7d1a7cdcd6a23af599b4d/python/profile_processor/profile_processor.py).
Other checks cover the momentum-conserving rigid shift, finite-difference
tangents, adjoint pairing, cache reuse and invalid inputs.

The high-collisionality regression compares the flow with the full ion-ion
collision limit in [Landreman's refined calculation](https://terpconnect.umd.edu/~mattland/assets/notes/20131019-03-PfirschSchluterRegimeIonNeoclassicalFlow.pdf),
equations (2) and (5), which neglects ion-electron collisions as this
single-ion path does. The limit is `k = k0 + k1 R`, with `k0 = -1.97213`,
`k1 = -0.09361` and the geometry factor below. Here `B` is magnetic-field
strength and `grad_parallel` means differentiation along the magnetic field:
`R = <B^2> <(grad_parallel ln B)^2> / <(grad_parallel B)^2>`.
Angle brackets denote the flux-surface average. The regression uses the
Boozer parallel derivative and surface measure of the solver when evaluating
this geometry factor, and extrapolates two finite-collisionality solutions in
inverse collisionality squared. The often quoted `-2.1` is an approximate
limit; the refined coefficients give about `-2.06574` when `R` is one.
This checks asymptotic consistency without supplying a continuous error bound.

This path does not yet provide a continuous error certificate, a qualified
resolution indicator, or the surface-profile time and memory comparison with
NEO-2 required by [#216](https://github.com/itpplasma/kin6d/issues/216).
It does not establish full three-dimensional neoclassical toroidal viscosity
or accuracy on the AUG 30835 and ITER TC24 cases.
