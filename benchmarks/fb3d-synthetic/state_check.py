"""Acceptance of a converged fb3d_synthetic state beyond the residual.

The residual gate admits several discrete equilibria (states 13-16 mm apart
were seen), so a run is accepted only if, in addition to the residual,
- the state is a strict local minimum of the reduced energy in the slide gauge
  (In_-(S_red) = 0, "reduced_negative" in <prefix>.json);
- it is stellarator symmetric, as the coils, the chart and the frozen labels
  are: the boundary <prefix>.lcfs at phi-hat and its image under
  (x, y, z) -> (x, -y, -z) at -phi-hat lie within --symmetry (Hausdorff
  distance, largest over phi-hat). The direct triangulation is mirror
  symmetric at R1 (12 edges, chain 16) but not at R2 or R3, so the discrete
  problem carries an O(h) asymmetry there unless the mesh is mirrored
  (fb3d_synthetic argument 19); the plasma-map defect and the axis Z of the
  JSON are exact checks only when "mirrored_nodes" equals "nodes";
- its phi = 0 section agrees with the coil-field flux surface of the vacuum
  oracle (vacuum_reference.py, <prefix>.oracle.json) within --oracle-rms and
  --oracle-max;
- if a reference boundary was given to fb3d_synthetic, the largest pointwise
  distance to it is below --reference.
Prints one JSON line per prefix; the exit status is the number of rejected runs.
"""
import argparse
import json
import sys

import numpy as np


def boundary_symmetry_defect(path):
    k = np.loadtxt(path)
    m = len(np.unique(k[:, 0]))
    p = k[:, 1:4].reshape(m, -1, 3)
    d = 0.0
    for i in range(m):
        dist = np.linalg.norm((p[i] * [1, -1, -1])[:, None] - p[(-i) % m][None], axis=2)
        d = max(d, dist.min(1).max(), dist.min(0).max())
    return float(d)


ap = argparse.ArgumentParser()
ap.add_argument("prefix", nargs="+")
ap.add_argument("--residual", type=float, default=1e-8)
ap.add_argument("--symmetry", type=float, default=2e-3, help="m")
ap.add_argument("--oracle-rms", type=float, default=3.5e-3, help="m")
ap.add_argument("--oracle-max", type=float, default=7.0e-3, help="m")
ap.add_argument("--reference", type=float, default=1e-6, help="m")
a = ap.parse_args()

rejected = 0
for prefix in a.prefix:
    r = json.load(open(prefix + ".json"))
    s = r.get("state", {})
    checks = {
        "residual": r["status"] == 0 and r["residual"] <= a.residual,
        "stable": s.get("spectrum_status") == 0 and s.get("reduced_negative") == 0,
    }
    sym = boundary_symmetry_defect(prefix + ".lcfs")
    checks["symmetric"] = bool(sym <= a.symmetry)
    if s.get("mirrored_nodes") == s.get("nodes"):
        checks["symmetric"] = checks["symmetric"] and s["symmetry_defect"] <= a.symmetry \
            and abs(s["axis_Z_phi0"]) <= a.symmetry
    if s.get("reference_max", -1) >= 0:
        checks["reference"] = s["reference_max"] <= a.reference
    try:
        o = json.load(open(prefix + ".oracle.json"))[prefix]
        checks["oracle"] = o["rms_m"] <= a.oracle_rms and o["max_m"] <= a.oracle_max
    except (OSError, KeyError, ValueError):
        o = None
    ok = all(checks.values())
    rejected += not ok
    print(json.dumps(dict(run=prefix, accepted=ok, checks=checks,
                          reduced_energy=s.get("reduced_energy"),
                          boundary_symmetry_defect=sym,
                          symmetry_defect=s.get("symmetry_defect"),
                          axis_Z_phi0=s.get("axis_Z_phi0"),
                          reference_max=s.get("reference_max"), oracle=o)))
sys.exit(min(rejected, 125))
