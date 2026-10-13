"""Smoke check of nn1_spec_case output: KIN6D NN1 field vs the direct Biot-Savart sum of the same
filaments (columns 7-9 of the output), relative to the RMS non-toroidal-field deviation.
usage: check_nn1_spec_case.py out.B.txt [tolerance]"""
import math
import sys

rows = [list(map(float, ln.split())) for ln in open(sys.argv[1]) if not ln.startswith("#")]
tol = float(sys.argv[2]) if len(sys.argv) > 2 else 0.15
assert len(rows) >= 100, f"too few points: {len(rows)}"
num = den = 0.0
for r in rows:
    R, bk, bd = r[0], r[3:6], r[6:9]
    num += sum((a - b) ** 2 for a, b in zip(bk, bd))
    den += bd[0] ** 2 + bd[1] ** 2 + (bd[2] - 10.0 / R) ** 2
rel = math.sqrt(num / den)
print(f"relative RMS deviation from direct Biot-Savart {rel:.4e} (tolerance {tol})")
if not rel < tol:
    sys.exit(1)
