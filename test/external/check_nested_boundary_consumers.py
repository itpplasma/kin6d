"""Optional reader checks; scientific definitions and exports remain in Fortran.

Run in an environment containing the selected external code. No external
solver convergence or matched-profile claim is made by this reader fixture.
"""
import argparse
import importlib.metadata
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("consumer", choices=("vmecpp", "desc", "gvec"))
parser.add_argument("input", type=Path)
args = parser.parse_args()
path = args.input.resolve()
record = {"consumer": args.consumer, "input": str(path), "claim": "boundary reader"}
if args.consumer == "vmecpp":
    import numpy as np
    import vmecpp

    value = vmecpp.VmecInput.from_file(path)
    assert value.nfp == 2
    assert np.all(np.isfinite(value.rbc))
    assert np.all(np.isfinite(value.zbs))
    if value.lasym:
        assert np.max(np.abs(value.rbs)) + np.max(np.abs(value.zbc)) > 1e-5
    record.update(version=importlib.metadata.version("vmecpp"), lasym=value.lasym)
elif args.consumer == "desc":
    import numpy as np
    from desc.geometry import FourierRZToroidalSurface
    from desc.grid import Grid
    from netCDF4 import Dataset

    value = FourierRZToroidalSurface.from_input_file(str(path))
    assert value.NFP == 2
    with Dataset(path.parent / "nested_oracle_export.nc") as nc:
        points = np.asarray(nc["boundary_points"][:])
        error_indicator = float(nc["boundary_sampled_error"][:])
    count = round(np.sqrt(len(points)))
    assert count * count == len(points)
    indices = np.arange(0, len(points), max(1, len(points) // 157))
    theta = 2 * np.pi * (indices % count) / count
    phi = 2 * np.pi * (indices // count) / count
    grid = Grid(np.column_stack((np.ones(len(indices)), theta, phi)), sort=False)
    values = value.compute(["R", "Z"], grid=grid)
    expected = points[indices]
    discrepancy = np.hypot(
        np.asarray(values["R"]) - np.hypot(expected[:, 0], expected[:, 1]),
        np.asarray(values["Z"]) - expected[:, 2],
    ).max()
    assert discrepancy < max(1e-10, 2 * error_indicator)
    record.update(version=importlib.metadata.version("desc-opt"),
                  sampled_position_error=float(discrepancy), sym=value.sym)
else:
    import f90nml
    import gvec

    value = gvec.util.parameters_from_vmec(f90nml.read(path)["indata"], path.name)
    flipped = gvec.util.flip_parameters_zeta(value)
    assert flipped["nfp"] == 2
    if f90nml.read(path)["indata"].get("lasym", False):
        sine_r = max(map(abs, flipped["X1_b_sin"].values()))
        cosine_z = max(map(abs, flipped["X2_b_cos"].values()))
        assert sine_r + cosine_z > 1e-5
        record.update(sine_R_norm=sine_r, cosine_Z_norm=cosine_z)
    output = path.parent / "gvec-boundary.yaml"
    gvec.util.write_parameters(flipped, output)
    loaded = gvec.util.read_parameters(output)
    assert loaded["nfp"] == 2
    assert loaded["X1_b_cos"] == flipped["X1_b_cos"]
    assert loaded["X2_b_sin"] == flipped["X2_b_sin"]
    record.update(version=importlib.metadata.version("gvec"), output=str(output))
print("[V] " + json.dumps(record, sort_keys=True))
