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

The local rotation-ratio and force-balance values, tangents and adjoints come from the
retained definition in `gen/gen_dk_energy.f90`. The force-balance value test
uses a synthetic input evaluated independently with
[KAMEL's profile equations](https://github.com/itpplasma/KAMEL/blob/266210c5f030dae287c7d1a7cdcd6a23af599b4d/python/profile_processor/profile_processor.py).
Other checks cover the momentum-conserving rigid shift, finite-difference
tangents, adjoint pairing, cache reuse and invalid inputs.

This path does not yet provide a continuous error certificate, a qualified
resolution indicator, or the surface-profile time and memory comparison with
NEO-2 required by [#216](https://github.com/itpplasma/kin6d/issues/216).
It does not establish full three-dimensional neoclassical toroidal viscosity
or accuracy on the AUG 30835 and ITER TC24 cases.
