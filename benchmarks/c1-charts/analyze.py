#!/usr/bin/env python3
"""Summarise c1_charts sweep CSVs: error versus DOFs per (chart, p, chi), and
the DOFs, PCG iterations and solver time (solve minus coefficient callback)
needed to reach a target relative L2 error, interpolated log-log between
adjacent converged rungs. Rows that did not converge are flagged, not used."""
import argparse
import csv
import math
from collections import defaultdict


def load(paths):
    rows = []
    for path in paths:
        with open(path) as handle:
            for row in csv.DictReader(handle):
                if row.get("status", "").strip() != "0" and row.get("status", "").strip() != "1":
                    continue
                rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", nargs="+")
    parser.add_argument("--target", type=float, action="append")
    args = parser.parse_args()
    targets = args.target or [1e-2, 1e-3]
    groups = defaultdict(list)
    for row in load(args.csv):
        key = (row["chart"], int(row["p"]), float(row["chi_ratio"]))
        solver = float(row["t_solve_s"]) - float(row["t_coef_in_solve_s"])
        groups[key].append((int(row["ndof"]), float(row["l2_rel"]), int(row["iterations"]),
                            solver, int(row["status"]), int(row["n"])))
    print("chart,p,chi,ladder(n:ndof:l2_rel:iters:solver_s[:NC])")
    for key in sorted(groups, key=lambda k: (k[2], k[1], k[0])):
        ladder = sorted(groups[key])
        text = " ".join(f"{n}:{d}:{e:.2e}:{i}:{t:.2f}" + (":NC" if s else "")
                        for d, e, i, t, s, n in ladder)
        print(f"{key[0]},{key[1]},{key[2]:.0e},{text}")
    print()
    print("target,chi,p,chart,ndof_at_target,iters_at_target,solver_s_at_target")
    for target in targets:
        for key in sorted(groups, key=lambda k: (k[2], k[1], k[0])):
            ladder = [r for r in sorted(groups[key]) if r[4] == 0]
            hit = None
            for a, b in zip(ladder, ladder[1:]):
                if a[1] > target >= b[1] and b[1] > 0 and a[1] > 0:
                    w = math.log(target / a[1]) / math.log(b[1] / a[1])
                    interp = lambda x, y: math.exp(math.log(x) + w * (math.log(y) - math.log(x)))
                    hit = (interp(a[0], b[0]), interp(max(a[2], 1), max(b[2], 1)),
                           interp(max(a[3], 1e-3), max(b[3], 1e-3)))
                    break
            if ladder and ladder[0][1] <= target:
                hit = ("<=" + str(ladder[0][0]), ladder[0][2], ladder[0][3])
            if hit is None:
                print(f"{target:.0e},{key[2]:.0e},{key[1]},{key[0]},not reached,,")
            elif isinstance(hit[0], str):
                print(f"{target:.0e},{key[2]:.0e},{key[1]},{key[0]},{hit[0]},{hit[1]},{hit[2]:.2f}")
            else:
                print(f"{target:.0e},{key[2]:.0e},{key[1]},{key[0]},{hit[0]:.0f},{hit[1]:.0f},{hit[2]:.2f}")


if __name__ == "__main__":
    main()
