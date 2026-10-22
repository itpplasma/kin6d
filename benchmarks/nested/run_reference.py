"""Optional external-code adapter; inputs and exact oracles remain in Fortran.

Run each point in a separate Slurm task. Equal-weight interior samples are
diagnostics, not volume norms or continuous bounds. VMEC sampling uses SIMSOPT's
native radial splines, not a DESC fit or an additional equilibrium solve.
"""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import resource
import time

import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("consumer", choices=("vmecpp", "desc", "gvec"))
parser.add_argument("input", type=Path)
parser.add_argument("resolution", type=int)
parser.add_argument("output", type=Path)
parser.add_argument("--clockwise", action="store_true",
                    help="VMEC boundary uses clockwise poloidal angle (integer oracles)")
parser.add_argument("--axisymmetric", action="store_true",
                    help="retain only n=0 in VMEC's discretization")
parser.add_argument("--ns", type=int, help="VMEC radial points (default: 4*resolution+1)")
parser.add_argument("--ftol", type=float, default=1e-8,
                    help="VMEC squared-force stopping tolerance")
args = parser.parse_args()
if args.consumer != "vmecpp" and (args.axisymmetric or args.ns is not None
                                  or args.ftol != 1e-8):
    parser.error("--axisymmetric, --ns and --ftol apply only to VMEC++")
if args.ns is not None and args.ns < 5:
    parser.error("--ns must be at least 5")
if not np.isfinite(args.ftol) or args.ftol <= 0:
    parser.error("--ftol must be finite and positive")
args.input = args.input.resolve()
args.output = args.output.resolve()
args.output.mkdir(parents=True, exist_ok=True)
os.chdir(args.output)
m = args.resolution
mu0 = 4 * np.pi * 1e-7  # VMEC/GVEC convention, also used by the Fortran export.
s = np.array([0.25, 0.5, 0.75])
angles = 2 * np.pi * (np.arange(24) + 0.37) / 24
record = dict(consumer=args.consumer, input=str(args.input), resolution=m,
              cores=int(os.environ.get("OMP_NUM_THREADS", "1")),
              claim="candidate; sampled diagnostics; bounds unavailable")
record["input_sha256"] = hashlib.sha256(args.input.read_bytes()).hexdigest()
Path("original.input").write_bytes(args.input.read_bytes())


def cylindrical_to_cartesian(values, phi):
    """Coordinate adapter only; Cartesian physical vectors are canonical."""
    values = np.asarray(values)
    return np.column_stack((values[:, 0] * np.cos(phi) - values[:, 1] * np.sin(phi),
                            values[:, 0] * np.sin(phi) + values[:, 1] * np.cos(phi),
                            values[:, 2]))


started = time.perf_counter()
if args.consumer == "vmecpp":
    import vmecpp

    value = vmecpp.VmecInput.from_file(args.input)
    value.mpol = m + 1
    value.ntor = 0 if args.axisymmetric else m
    value.ns_array = np.array([args.ns if args.ns is not None else 4 * m + 1])
    value.ftol_array = np.array([args.ftol])
    value.niter_array = np.array([5000])
    value.save("resolved.json")
    result = vmecpp.run(value, max_threads=record["cores"], verbose=True)
    record["solve_wall_s"] = time.perf_counter() - started
    result.wout.save("wout_reference.nc")
    record.update(version=importlib.metadata.version("vmecpp"),
                  phiedge=float(result.wout.phi[-1]),
                  fsqr=float(result.wout.fsqr), fsqz=float(result.wout.fsqz),
                  fsql=float(result.wout.fsql), ier_flag=int(result.wout.ier_flag))
    from simsopt.mhd import Vmec
    from simsopt.mhd.vmec_diagnostics import vmec_splines, vmec_compute_geometry

    value = Vmec("wout_reference.nc")
    splines = vmec_splines(value)
    data = vmec_compute_geometry(splines, s, angles, angles / value.wout.nfp)
    shape = data.X.shape
    def vector(*names):
        return np.column_stack([np.asarray(getattr(data, key)).ravel() for key in names])
    xyz = vector("X", "Y", "Z")
    field = vector("B_X", "B_Y", "B_Z")
    # Curl from SIMSOPT's Cartesian field derivatives; no force-balance assumption.
    current = np.column_stack((
        (data.grad_B__ZY - data.grad_B__YZ).ravel(),
        (data.grad_B__XZ - data.grad_B__ZX).ravel(),
        (data.grad_B__YX - data.grad_B__XY).ravel()))
    pressure = np.broadcast_to(mu0 * splines.pressure(s)[:, None, None], shape).ravel()
    gradp = vector("grad_s_X", "grad_s_Y", "grad_s_Z")
    gradp *= np.broadcast_to(mu0 * splines.d_pressure_d_s(s)[:, None, None], shape).ravel()[:, None]
    labels = np.broadcast_to(s[:, None, None], shape).ravel()
    record["iota_sample"] = splines.iota(s).tolist()
    # VMEC enforces its own negative-Jacobian poloidal orientation, flipping
    # a clockwise input's theta and iota together (boundaries.cc:flipTheta).
    record["iota_input_convention"] = ((-1 if args.clockwise else 1) * splines.iota(s)).tolist()
elif args.consumer == "desc":
    import desc
    import f90nml
    desc.set_device("cpu")
    from desc.equilibrium import Equilibrium
    from desc.grid import LinearGrid

    namelist = f90nml.read(args.input)["indata"]
    axis = [namelist[key] for key in ("raxis_cc", "raxis_cs", "zaxis_cc", "zaxis_cs")
            if key in namelist]
    automatic_axis = bool(axis) and all(np.all(np.asarray(x, dtype=float) == 0) for x in axis)
    # A zero VMEC axis asks VMEC to construct a guess; in DESC it is a
    # literal curve through R=0. Use DESC's own boundary-centre guess.
    options = {"axis": None} if automatic_axis else {}
    value = Equilibrium.from_input_file(str(args.input), **options)
    record["automatic_axis"] = automatic_axis
    value.change_resolution(L=2*m, M=m, N=m)
    value, result = value.solve(ftol=1e-6, xtol=1e-6, gtol=1e-6,
                                maxiter=100, verbose=3)
    record["solve_wall_s"] = time.perf_counter() - started
    record["optimizer_success"] = bool(result["success"])
    value.save("equilibrium.h5")
    grid = LinearGrid(rho=np.sqrt(s), theta=angles, zeta=angles/value.NFP,
                      NFP=value.NFP, sym=False)
    data = value.compute(["R", "phi", "Z", "B", "J", "p", "grad(p)", "iota"], grid=grid)
    phi = np.asarray(data["phi"])
    xyz = np.column_stack((data["R"]*np.cos(phi), data["R"]*np.sin(phi), data["Z"]))
    field = cylindrical_to_cartesian(data["B"], phi)
    current = mu0 * cylindrical_to_cartesian(data["J"], phi)
    gradp = mu0 * cylindrical_to_cartesian(data["grad(p)"], phi)
    pressure = mu0 * np.asarray(data["p"])
    labels = grid.nodes[:, 0]**2
    record.update(version=importlib.metadata.version("desc-opt"),
                  iota_sample=np.asarray(data["iota"])[grid.unique_rho_idx].tolist())
else:
    import f90nml
    import gvec

    parameters = gvec.util.parameters_from_vmec(f90nml.read(args.input)["indata"], args.input.name)
    parameters = gvec.util.flip_parameters_zeta(parameters)
    if args.clockwise:
        # GVEC's RZ chart embeds phi=-zeta and requires positive Jacobians.
        # Reverse theta as well for a clockwise VMEC boundary, transforming
        # iota with both flips. Physical Cartesian fields are unchanged.
        parameters = gvec.util.flip_parameters_theta(parameters)
    parameters.update(X1_mn_max=[m, m], X2_mn_max=[m, m], LA_mn_max=[m, m],
                      sgrid=dict(grid_type=0, nElems=max(3, m)),
                      totalIter=20000, minimize_tol=1e-6, logIter=100)
    gvec.util.write_parameters(parameters, "resolved.ini")
    result = gvec.run(parameters, runpath=Path.cwd()/"solve", quiet=False)
    record["solve_wall_s"] = time.perf_counter() - started
    record.update(success=bool(result.max_force <= 1e-6),
                  native_force=float(result.max_force), iterations=int(result.GVEC_iter_used))
    value = result.state
    data = gvec.evaluate(value, "pos", "B", "J", "p", "dp_dr", "grad_rho", "iota",
                         rho=np.sqrt(s), theta=angles, zeta=angles/value.nfp)
    def vector(name):
        return np.asarray(data[name].transpose("rad", "pol", "tor", "xyz")).reshape(-1, 3)
    xyz = vector("pos")
    field = vector("B")
    current = mu0 * vector("J")
    template = data.pos.isel(xyz=0, drop=True)
    radial_derivative = data.dp_dr.broadcast_like(template).transpose("rad", "pol", "tor")
    gradp = mu0 * vector("grad_rho") * np.asarray(radial_derivative).ravel()[:, None]
    pressure = mu0 * np.asarray(data.p.broadcast_like(template).transpose("rad", "pol", "tor")).ravel()
    labels = np.broadcast_to(s[:, None, None], template.shape).ravel()
    record.update(version=importlib.metadata.version("gvec"),
                  iota_sample=np.asarray(data.iota).tolist())

np.savetxt("samples.txt", np.column_stack((labels, xyz, field, pressure, current, gradp)))
record["total_wall_s"] = time.perf_counter() - started
record["process_peak_rss_kib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
Path("metrics.json").write_text(json.dumps(record, indent=2) + "\n")
print("[V] " + json.dumps(record, sort_keys=True))
