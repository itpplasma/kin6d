# Island reconstruction with TIAGO

This optional benchmark reconstructs the perturbation amplitude of a synthetic
helical-coil annulus from point polarization angles. KIN6D supplies its existing
fixed-wall nonnested stationary field and implicit derivatives. TIAGO supplies
its existing motional Stark effect (MSE) polarization operator and maximum a
posteriori (MAP) optimizer. Neither implementation is copied into the adapter.

Configure the normal KIN6D build with `-DKIN6D_TIAGO_SOURCE_DIR=/path/to/tiago`,
using TIAGO with the [shared core targets](https://github.com/itpplasma/tiago/pull/76)
that provide `cmake/TiagoInferenceCore.cmake`. The
adapter links the existing `fortnum` target selected by KIN6D. Run:

```sh
cmake --build build --target test_tiago_island
ctest --test-dir build -R '^tiago_island_(vacuum|pressure)$' --output-on-failure
```

The slow `tiago_island_direct` case takes its observations from the segment
Biot–Savart field of the coil with all its toroidal harmonics instead of the
fitted model, and fits them with the P3, nr=8, three-mode provider. At the
default P2, nr=4, one-mode resolution the same fit gives amplitude 1.020 and
a held-out error of 69 sigma. `island_direct_coil_split` shows that the
coarse field alone has a pitch error of 54 sigma at those points. The direct
case uses tolerances of 1e-2 on the amplitude and 1 sigma on held-out angles,
because its observations come from a different model than the fit.

The two cases use zero pressure source or the regularized finite-pressure
stationary model. Each checks zero and nonzero perturbations. Sixteen synthetic
channels are fitted; eight at a different toroidal plane are held out. Checks
include independently initialized finite differences of the composed field and
polarization map, tangent/adjoint pairing, parameter recovery and held-out
prediction. The single inferred parameter is the coil perturbation amplitude;
pressure source, current drive and parallel diffusion are fixed. The assigned
angle uncertainty is 0.0001 radians. The data are noiseless and generated in the
same model family, so these checks do not establish posterior coverage or
cross-model topology discrimination.

The geometry, coils and observations are synthetic. This is not a machine
reconstruction or a moving-boundary island equilibrium. In particular, the
finite-pressure stationary closure is not automatically ideal MHD or
multi-region relaxed MHD. Physical force-balance defects and numerical
refinement remain separate acceptance conditions; see the
[nonnested model](../../docs/equilibrium-nonnested-nn2.md).
