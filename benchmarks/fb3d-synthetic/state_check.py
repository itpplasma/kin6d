"""Acceptance of a converged fb3d_synthetic state beyond the residual.

The residual gate admits several discrete equilibria (states 13-16 mm apart
were seen), so a run is accepted only if, in addition to the residual,
- the state is a strict local minimum of the reduced energy in the slide gauge
  (In_-(S_red) = 0, "reduced_negative" in <prefix>.json);
- it is stellarator symmetric, as the coils, the chart, the reference mesh and
  the frozen labels are: the plasma-map symmetry defect and the axis Z at
  phi-hat = 0 (a symmetry plane) are below --symmetry;
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

ap = argparse.ArgumentParser()
ap.add_argument("prefix", nargs="+")
ap.add_argument("--residual", type=float, default=1e-8)
ap.add_argument("--symmetry", type=float, default=1e-4, help="m")
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
        "symmetric": s.get("mirrored_nodes") == s.get("nodes")
        and s.get("symmetry_defect", 1.0) <= a.symmetry
        and abs(s.get("axis_Z_phi0", 1.0)) <= a.symmetry,
    }
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
                          symmetry_defect=s.get("symmetry_defect"),
                          axis_Z_phi0=s.get("axis_Z_phi0"),
                          reference_max=s.get("reference_max"), oracle=o)))
sys.exit(min(rejected, 125))
