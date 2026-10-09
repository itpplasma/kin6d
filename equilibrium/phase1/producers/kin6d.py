"""Phase 1 producer: native Fortran solve, contours and physical field readback.

Sign ledger: equilibrium/CASE_CONTRACT.md, "Signs, units and selectors":
COCOS 3, BR=-psi_Z/R, BZ=psi_R/R, Bphi=F/R, q=dPhi/dpsi/(2*pi).
The native solver already uses these signs and SI; the producer applies no flips.
Set KIN6D_BINARY to the built gs_phase1 executable when installed elsewhere.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import io
import os
from pathlib import Path
import subprocess

import numpy as np


@dataclass
class Result:
    axis: tuple[float, float]
    volume: float
    phi_edge: float
    dof: int
    wall_s: float
    native_converged: bool
    meta: dict
    _binary: Path
    _workdir: Path
    _f_edge: float
    _ffprime: float
    _profile: np.ndarray

    def evaluate(self, R, Z) -> dict:
        r, z = np.broadcast_arrays(np.asarray(R, float), np.asarray(Z, float))
        if not r.size:
            return {key: np.empty(r.shape) for key in ("psi", "BR", "BZ", "Bphi")}
        request = io.StringIO()
        request.write(f"{r.size}\n")
        np.savetxt(request, np.column_stack((r.ravel(), z.ravel())), fmt="%.17e")
        environment = _environment()
        environment["TMPDIR"] = str(self._workdir / "tmp")
        process = subprocess.run(
            [str(self._binary), "--evaluate", "mesh.dat", repr(self._f_edge),
             repr(self._ffprime)], cwd=self._workdir, input=request.getvalue(),
            text=True, capture_output=True, check=True, timeout=300, env=environment,
        )
        values = np.loadtxt(io.StringIO(process.stdout), ndmin=2)
        if values.shape != (r.size, 4):
            raise RuntimeError("Incomplete native field readback")
        return {key: values[:, i].reshape(r.shape)
                for i, key in enumerate(("psi", "BR", "BZ", "Bphi"))}

    def q(self, s_tor) -> np.ndarray:
        return np.interp(np.asarray(s_tor, float), self._profile[:, 1],
                         self._profile[:, 2], left=np.nan, right=np.nan)


def _environment() -> dict:
    return {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1", "BLIS_NUM_THREADS": "1"}


def run(case: dict, resolution: dict, workdir: Path) -> Result:
    """Solve constant pprime/ffprime data with a periodic spline boundary.

    source contains pprime and ffprime in SI (directly or in coefficients).
    Resolution keys: boundary_nodes, max_area, radial_points, angular_points,
    curved. Defaults are resolved by Fortran and recorded in output/gs_phase1.nc.
    Dense boundary samples must follow a uniform periodic parameter, omit the
    duplicate endpoint, and run counter-clockwise in the R,Z plane.
    """
    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    binary = Path(os.environ.get(
        "KIN6D_BINARY", Path(__file__).resolve().parents[3] / "build" / "gs_phase1"
    )).resolve()
    if not binary.is_file():
        raise FileNotFoundError(f"Build the gs_phase1 CMake target or set KIN6D_BINARY: {binary}")
    if case["source"]["kind"] != "solovev":
        raise ValueError("Phase 1 producer supports constant-source Solovev laws")
    if float(case["psi_edge"]) != 0:
        raise ValueError("Phase 1 requires psi_edge=0")
    source = case["source"].get("coefficients", case["source"])
    pprime, ffprime = float(source["pprime"]), float(source["ffprime"])
    f_edge = float(case["F_edge"])
    boundary = np.column_stack((case["boundary"]["R"], case["boundary"]["Z"]))
    if len(boundary) > 4 and np.array_equal(boundary[0], boundary[-1]):
        boundary = boundary[:-1]
    if boundary.ndim != 2 or boundary.shape[1] != 2 or len(boundary) < 4:
        raise ValueError("A closed R,Z curve requires at least four samples")
    if not np.isfinite(boundary).all():
        raise ValueError("Boundary samples must be finite")
    with (workdir / "boundary.dat").open("w") as stream:
        stream.write(f"{len(boundary)}\n")
        np.savetxt(stream, boundary, fmt="%.17e")
    allowed = {"boundary_nodes", "max_area", "radial_points", "angular_points", "curved"}
    if set(resolution) - allowed:
        raise ValueError(f"Unsupported KIN6D settings: {set(resolution) - allowed}")
    lines = ["&gs_phase1", " boundary_file='boundary.dat', mesh_file='mesh.dat',",
             f" pprime={pprime:.17e}, ffprime={ffprime:.17e}, f_edge={f_edge:.17e},"]
    for key, value in resolution.items():
        if key == "curved":
            literal = ".true." if value else ".false."
        elif key == "max_area":
            literal = f"{float(value):.17e}"
        else:
            if int(value) != value:
                raise ValueError(f"{key} must be an integer")
            literal = str(int(value))
        lines.append(f" {key}={literal},")
    lines.append("/\n")
    namelist = "\n".join(lines)
    (workdir / "input.nml").write_text(namelist)
    environment = _environment()
    (workdir / "tmp").mkdir(exist_ok=True)
    environment["TMPDIR"] = str(workdir / "tmp")
    process = subprocess.run([str(binary), "--config", "input.nml"], cwd=workdir,
                             text=True, capture_output=True, timeout=300, env=environment)
    (workdir / "native.log").write_text(process.stdout + process.stderr)
    process.check_returncode()
    summary = dict(line.split(maxsplit=1)
                   for line in (workdir / "summary.dat").read_text().splitlines())
    profile = np.loadtxt(workdir / "profiles.dat", ndmin=2)
    if np.any(np.diff(profile[:, 1]) <= 0) or not np.isfinite(profile).all():
        raise RuntimeError("Nonmonotone or nonfinite native toroidal-flux profile")
    residual = float(summary["residual"])
    meta = {
        "code_commit": summary["code_commit"], "code_tree": summary["code_tree"],
        "binary_path": str(binary), "settings": dict(resolution),
        "resolved_input": "output/gs_phase1.nc:kin6d_config_resolved",
        "input_hash": hashlib.sha256(
            namelist.encode() + (workdir / "boundary.dat").read_bytes()).hexdigest(),
        "case_id": case["id"], "convention": "COCOS 3, SI, psi Wb/rad, psi_edge=0",
        "sign_ledger": "equilibrium/CASE_CONTRACT.md#signs-units-and-selectors",
        "residual": residual, "flux_consistency": float(summary["flux_consistency"]),
        "threads": 1, "iterations": int(summary["iterations"]),
    }
    return Result((float(summary["axis_R"]), float(summary["axis_Z"])),
                  float(summary["volume"]), float(summary["phi_edge"]),
                  int(summary["dof"]), float(summary["wall_s"]),
                  bool(np.isfinite(residual) and residual < 1e-8), meta,
                  binary, workdir, f_edge, ffprime, profile)
