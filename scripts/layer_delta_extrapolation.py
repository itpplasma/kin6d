#!/usr/bin/env python3
"""Delta asymptotics of the ideal layer jump (closure = 'frobenius').

S(delta) = layer_kink * (xs/delta)^eps is the small-solution jump extracted
from the two-sided fit (columns 1, sigma t, t^(1+eps), t^2, sigma t^(2+eps),
t = |x|/delta at t_k = ratio^(k-1), k = 1..nfit).  An even term a |x|^p of the
outer flux that is not in this span enters S as

    S - J_inf = 2 (1 + eps) a pi_p xs^eps delta^(p - 1 - eps),

pi_p the weight of t^p on the t^(1+eps) column (ratio, nfit); odd terms drop
out because the rings are symmetric.  Mode 'model' checks this on the
manufactured flux of test_mhd_closure layer_kink_finite_mu at three deltas
and prints the fitted exponent.  Mode 'fit' reads prefix_layer.dat of TC24
runs at several deltas and extrapolates S(delta) = J_inf + C delta^p for a
given p (default 1, the mu = 1/2 form) per (phasing, rational, ratio).

usage: layer_delta_extrapolation.py model
       layer_delta_extrapolation.py fit [--p P] [--pfree] file_d<delta>_ra<r>_layer.dat ...
"""
import argparse
import re
import sys

import numpy as np


def kink_weights(eps, ratio, nfit):
    t = ratio ** np.arange(nfit)
    rows = []
    for sg in (-1.0, 1.0):
        for tk in t:
            rows.append([1.0, sg * tk, tk ** (1 + eps), tk**2, sg * tk ** (2 + eps)])
    return t, np.linalg.pinv(np.array(rows))[2]


def kink(phi_of_x, eps, delta, ratio, nfit, xs):
    t, w = kink_weights(eps, ratio, nfit)
    ph = np.concatenate([phi_of_x(-delta * t), phi_of_x(delta * t)])
    return 2 * (1 + eps) * (w @ ph) / delta * (xs / delta) ** eps


def model():
    bb, cc, b1, b2, c2, c4 = 0.8, -2.3, 3.0, 200.0, 40.0, 50.0
    xs, deltas = 5e-4, np.array([1e-3, 5e-4, 2.5e-4])
    print("term                     mu    ratio  predicted p-1-eps  fitted")
    for mu in (0.5, 0.553, 0.747):
        eps = mu - 0.5
        jt = 2 * (1 + eps) * bb * xs**eps
        cases = {
            "test flux (|x|^(3+eps))": (
                lambda x: cc * x + bb * abs(x) ** (1 + eps) * (1 + b1 * x + b2 * x * x)
                + c2 * x * x + c4 * x**3, 2.0),
            "large |x|^-eps, fixed A": (
                lambda x: bb * abs(x) ** (1 + eps) + 1e-6 * abs(x) ** (-eps),
                -1 - 2 * eps),
            "large tail |x|^(2-eps)": (
                lambda x: bb * abs(x) ** (1 + eps) + 1e-3 * abs(x) ** (2 - eps),
                1 - 2 * eps),
        }
        for name, (f, pred) in cases.items():
            for ratio in (1.5, 2.0):
                s = np.array([kink(f, eps, d, ratio, 5, xs) for d in deltas]) - jt
                if mu == 0.5 and name.startswith("large"):
                    continue
                if np.all(np.abs(s) < 1e-13 * jt):
                    print(f"{name:24s} {mu:5.3f} {ratio:4.1f}   {pred:7.3f}   exact")
                    continue
                pf = np.log(s[0] / s[1]) / np.log(2.0)
                pf2 = np.log(s[1] / s[2]) / np.log(2.0)
                print(f"{name:24s} {mu:5.3f} {ratio:4.1f}   {pred:7.3f}   "
                      f"{pf:7.3f} {pf2:7.3f}")


def read_layer(files):
    data = {}
    for f in files:
        m = re.search(r"_d([\d.e-]+)_ra([\d.]+)_layer\.dat$", f)
        if not m:
            sys.exit(f"cannot parse delta/ratio from {f}")
        d, r = float(m[1]), m[2]
        for line in open(f):
            if line.startswith("#") or not line.strip():
                continue
            v = [float(x) for x in line.split()]
            key = (int(v[0]), round(abs(v[1]), 3), r)
            data.setdefault(key, {})[d] = (complex(v[8], v[9]), v[5])
    return data


def fit(files, p, pfree):
    rows = {1: "010", 2: "U", 3: "M", 4: "L"}
    data = read_layer(files)
    out = {}
    for (ph, q, r), sd in sorted(data.items()):
        d = np.array(sorted(sd))
        s = np.array([sd[x][0] for x in d])
        mu = sd[d[0]][1]
        pp = p
        if pfree and len(d) >= 3:
            ps = np.linspace(-0.5, 3.0, 701)
            res = []
            for pt in ps:
                a = np.stack([np.ones_like(d), d**pt], 1).astype(complex)
                c = np.linalg.lstsq(a, s, rcond=None)[0]
                res.append(np.linalg.norm(a @ c - s))
            pp = ps[int(np.argmin(res))]
        a = np.stack([np.ones_like(d), d**pp], 1).astype(complex)
        c = np.linalg.lstsq(a, s, rcond=None)[0]
        out[(ph, q, r)] = (c[0], mu, pp, s)
    print("row  q     mu    ratio  p    |S| (delta ascending)          |J_inf|")
    for (ph, q, r), (j, mu, pp, s) in out.items():
        print(f"{rows.get(ph, ph):4s} {q:5.3f} {mu:5.3f} {r:4s} {pp:5.2f} "
              + " ".join(f"{abs(x):.4f}" for x in s) + f"   {abs(j):.4f}")
    print("ratio spread |J(1.5) - J(2.0)|/mean|J| (%):")
    for (ph, q, r), (j, *_rest) in out.items():
        if r != "1.5" or (ph, q, "2.0") not in out:
            continue
        j2 = out[(ph, q, "2.0")][0]
        print(f"  {rows.get(ph, ph):4s} q {q:5.3f}: "
              f"{100 * abs(j - j2) / (0.5 * (abs(j) + abs(j2))):6.1f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["model", "fit"])
    ap.add_argument("files", nargs="*")
    ap.add_argument("--p", type=float, default=1.0)
    ap.add_argument("--pfree", action="store_true")
    a = ap.parse_args()
    if a.mode == "model":
        model()
    else:
        fit(a.files, a.p, a.pfree)


if __name__ == "__main__":
    main()
