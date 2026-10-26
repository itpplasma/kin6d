#!/usr/bin/env python3
"""Figures for the closure A+D island study (test/study_island_closure_a.f90).

Usage:
  plot_island_closure_a.py poincare OUT.png LABEL=DIR [LABEL=DIR ...]
  (a LABEL starting with "zoom" restricts that panel to the island region)
  plot_island_closure_a.py flattening OUT.png DIR [DIR ...]
  plot_island_closure_a.py convergence OUT.png CONNECT_DIR LABEL=DIR [...]
  plot_island_closure_a.py fitzpatrick OUT.png LEVEL=DIR,DIR,... [...]

poincare: p(R, Z, phi=0)/max p (computed, viridis) with the Poincare section
of the same converged field (one point per toroidal turn).
convergence: mesh/order study at fixed chi_par/chi_perp (scan directories
labelled e.g. "P3 nr16", each with eps = 0 and eps = 1 runs): min p/p_max and
the p < 0.01 p_max area against the ratio, G against the ratio, and the
eps = 0 diameter profile at 1e6 with the wall-connected lines of CONNECT_DIR.
fitzpatrick: the half-flattening ratio against W for each mesh level with
the fitted exponent and its standard error.
flattening: chord profiles p/p_max on Z = 0 through the O point, the O-point
slope ratio G against chi_par/chi_perp, G against W (chi_par/chi_perp)^(1/4)
and the half-flattening ratio against W with the Fitzpatrick (1995) w_c
proportional to (chi_perp/chi_par)^(1/4), i.e. slope -4 in log ratio vs log W.
"""
import os
import sys

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Okabe-Ito, colour-vision safe; markers/line styles add redundancy.
OI = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]
MARK = ["o", "s", "^", "D", "v", "P", "X"]
LS = ["-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1)), "-"]
plt.rcParams.update({"font.size": 10, "axes.labelsize": 10, "legend.fontsize": 8})


def poincare(out, items):
    n = len(items)
    nc = min(n, 3)
    nrow = -(-n // nc)
    fig, axes = plt.subplots(nrow, nc, figsize=(3.6 * nc, 3.9 * nrow),
                             constrained_layout=True, squeeze=False)
    axes = axes.ravel()
    for ax in axes[n:]:
        ax.set_visible(False)
    axes = axes[:n]
    for ax, item in zip(axes, items):
        label, d = item.rsplit("=", 1)
        p = np.loadtxt(os.path.join(d, "pressure.dat"))
        ng = int(round(np.sqrt(len(p))))
        R = p[:, 0].reshape(ng, ng)
        Z = p[:, 1].reshape(ng, ng)
        P = p[:, 2].reshape(ng, ng)
        P = P / np.nanmax(P)
        pmin = np.nanmin(P)
        cs = ax.contourf(R, Z, P, levels=np.linspace(0, 1, 21), cmap="viridis",
                         extend="min")
        if pmin < 0:
            ax.contour(R, Z, P, levels=[0], colors=OI[3], linewidths=0.8)
            label += "\n" + rf"min $p/p_{{\max}}={pmin:.3f}$, pink: $p=0$"
        ax.contour(R, Z, P, levels=np.linspace(0.05, 0.95, 10), colors="w",
                   linewidths=0.4)
        pc = np.loadtxt(os.path.join(d, "poincare.dat"))
        ax.plot(pc[:, 1], pc[:, 2], ".", ms=0.7, mew=0, color="k", rasterized=True)
        ax.set_aspect("equal")
        if label.startswith("zoom"):
            # Inner side of the island (O point at Z = 0, R ~ 5.5 m).
            ax.set_xlim(5.15, 5.85)
            ax.set_ylim(-0.4, 0.4)
            ax.lines[-1].set_markersize(1.5)
        ax.set_xlabel("R [m]")
        ax.set_title(label, fontsize=10)
    for ax in axes[::nc]:
        ax.set_ylabel("Z [m]")
    fig.colorbar(cs, ax=axes, shrink=0.8, label=r"$p/p_{\max}$ at $\varphi=0$")
    fig.savefig(out, dpi=200)


def load_scan(d):
    s = np.atleast_2d(np.loadtxt(os.path.join(d, "scan.dat")))
    c = np.loadtxt(os.path.join(d, "chords.dat"))
    return s, c


def half_ratio(ratio, G):
    """log-linear interpolation of the first crossing of G = 1/2."""
    for i in range(1, len(G)):
        if G[i - 1] > 0.5 >= G[i]:
            t = (G[i - 1] - 0.5) / (G[i - 1] - G[i])
            return 10 ** ((1 - t) * np.log10(ratio[i - 1]) + t * np.log10(ratio[i]))
    return np.nan


def flattening(out, dirs):
    scans = [load_scan(d) for d in dirs]
    fig, ax = plt.subplots(2, 2, figsize=(8.6, 7.0), constrained_layout=True)
    # (a) chord profiles for the eps = 1 run (or the last one given).
    pick = [i for i, (s, _) in enumerate(scans) if np.isclose(s[0, 0], 1.0)]
    s, c = scans[pick[0] if pick else -1]
    ratios = np.unique(c[:, 1])
    show = [r for r in ratios if np.isclose(np.log10(r) % 2, 0)]
    for j, r in enumerate(show):
        m = np.isclose(c[:, 1], r)
        ax[0, 0].plot(c[m, 2], c[m, 3], ls=LS[j % 7], color=OI[j % 7],
                      label=rf"$\chi_\parallel/\chi_\perp=10^{{{np.log10(r):.0f}}}$")
    row = s[np.isclose(s[:, 1], show[-1])][0]
    rO = 6.0 - row[3]
    for x in (rO - row[2] / 2, rO + row[2] / 2):
        ax[0, 0].axvline(x, color="0.5", lw=0.8, ls=":")
    ax[0, 0].axvline(rO, color="0.3", lw=0.8)
    ax[0, 0].set_xlabel(r"$R$ [m] ($Z=0$, $\varphi=0$, through the O point)")
    ax[0, 0].set_ylabel(r"$p/p_{\max}$")
    ax[0, 0].set_title(rf"(a) chord profiles, $\epsilon={s[0, 0]:g}$, "
                       r"island $|R-R_O|<W/2$ dotted", fontsize=9)
    ax[0, 0].legend(frameon=False)
    # (b) G against the ratio; (c) collapse; (d) half-flattening ratio.
    Ws, Rh = [], []
    for j, (s, _) in enumerate(scans):
        eps, ratio, W, G = s[0, 0], s[:, 1], s[:, 2], s[:, 4]
        Wm = np.median(W)
        kw = dict(color=OI[j % 7], marker=MARK[j % 7], ls=LS[j % 7], ms=4,
                  label=rf"$\epsilon={eps:g}$, $W={Wm:.3f}$ m")
        ax[0, 1].semilogx(ratio, G, **kw)
        ax[1, 0].semilogx(W * ratio ** 0.25, G, **kw)
        Ws.append(Wm)
        Rh.append(half_ratio(ratio, G))
    for a in (ax[0, 1], ax[1, 0]):
        a.axhline(0.5, color="0.5", lw=0.8, ls=":")
        a.set_ylabel(r"$G$ (O-point slope / chord slope)")
        a.legend(frameon=False)
    ax[0, 1].set_xlabel(r"$\chi_\parallel/\chi_\perp$")
    ax[0, 1].set_title("(b) flattening against anisotropy", fontsize=9)
    ax[1, 0].set_xlabel(r"$W\,(\chi_\parallel/\chi_\perp)^{1/4}$ [m]"
                        r"  $\propto W/w_c$")
    ax[1, 0].set_title("(c) Fitzpatrick scaling collapse", fontsize=9)
    ax[1, 0].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    Ws, Rh = np.array(Ws), np.array(Rh)
    # Fit only scans that reach full flattening (min G < 0.25): a declared
    # resolution criterion, since an under-resolved island saturates.
    full = np.array([s[:, 4].min() < 0.25 for s, _ in scans])
    ax[1, 1].loglog(Ws[np.isfinite(Rh) & ~full], Rh[np.isfinite(Rh) & ~full], "s",
                    mfc="none", color=OI[1], label="no full flattening (excluded)")
    ok = np.isfinite(Rh) & full
    ax[1, 1].loglog(Ws[ok], Rh[ok], "o", color=OI[0], label="computed ($G=1/2$)")
    if ok.sum() >= 2:
        slope, icpt = np.polyfit(np.log(Ws[ok]), np.log(Rh[ok]), 1)
        w = np.geomspace(Ws[ok].min(), Ws[ok].max(), 20)
        ax[1, 1].loglog(w, np.exp(icpt) * w ** slope, "-", color=OI[0],
                        label=f"fit slope {slope:.2f}")
        k = np.argmin(np.abs(np.log(Ws[ok]) - np.log(Ws[ok]).mean()))
        ax[1, 1].loglog(w, Rh[ok][k] * (w / Ws[ok][k]) ** -4, "--", color=OI[1],
                        label=r"$w_c\propto(\chi_\perp/\chi_\parallel)^{1/4}$: slope $-4$")
        print(f"half-flattening fit: slope {slope:.3f} over W = {Ws[ok]}")
    log_width_axis(ax[1, 1])
    ax[1, 1].set_xlabel(r"island width $W$ [m]")
    ax[1, 1].set_ylabel(r"$(\chi_\parallel/\chi_\perp)$ at $G=1/2$")
    ax[1, 1].set_title("(d) critical anisotropy against width", fontsize=9)
    ax[1, 1].legend(frameon=False)
    for (s, _), r in zip(scans, Rh):
        print(f"eps {s[0, 0]:g} W {np.median(s[:, 2]):.4f} ratio_half {r:.4g}")
    fig.savefig(out, dpi=200)


def log_width_axis(ax):
    """Plain decimal ticks on a log width axis without crowded minor labels."""
    lo, hi = ax.get_xlim()
    ticks = [t for t in (0.03, 0.05, 0.07, 0.1, 0.12, 0.15, 0.2, 0.3) if lo <= t <= hi]
    ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(ticks))
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())


def load_plane(d):
    return np.atleast_2d(np.loadtxt(os.path.join(d, "plane.dat")))


def convergence(out, connect_dir, items):
    runs = []
    for item in items:
        label, d = item.split("=", 1)
        runs.append((label, d))
    fig, ax = plt.subplots(2, 2, figsize=(8.6, 7.2), constrained_layout=True)
    for j, (label, d) in enumerate(runs):
        kw = dict(color=OI[j % 7], marker=MARK[j % 7], ms=3.5)
        for eps, ls in ((0.0, "--"), (1.0, "-")):
            dd = d.replace("eps1.0", f"eps{eps:.1f}").replace("eps0.0", f"eps{eps:.1f}")
            if not os.path.exists(os.path.join(dd, "plane.dat")):
                continue
            pl = load_plane(dd)
            lab = f"{label}" if eps == 1.0 else None
            ax[0, 0].semilogx(pl[:, 1], pl[:, 2], ls=ls, label=lab, **kw)
            ax[1, 1].semilogx(pl[:, 1], pl[:, 4], ls=ls, label=lab, **kw)
            if eps == 1.0:
                s = np.atleast_2d(np.loadtxt(os.path.join(dd, "scan.dat")))
                ax[0, 1].semilogx(s[:, 1], s[:, 4], ls=ls, label=label, **kw)
            else:
                c = np.loadtxt(os.path.join(dd, "diameter.dat"))
                m = np.isclose(c[:, 1], 1e6)
                if m.any():
                    R, P = c[m, 2], c[m, 3].copy()
                    P[np.r_[np.diff(R) > 0.1, False]] = np.nan  # the coil hole
                    ax[1, 0].plot(R, P, color=OI[j % 7], ls=LS[j % 7], lw=1,
                                  label=label)
    ax[0, 0].axhline(0, color="0.5", lw=0.8, ls=":")
    ax[0, 0].set_ylabel(r"min $p/p_{\max}$ on $\varphi=0$")
    ax[0, 0].set_title(r"(a) undershoot (solid $\epsilon=1$, dashed $\epsilon=0$)",
                       fontsize=9)
    ax[0, 1].set_ylabel(r"$G$ (O-point slope / chord slope)")
    ax[0, 1].set_title(r"(b) flattening, $\epsilon=1$", fontsize=9)
    ax[1, 1].set_ylabel(r"area fraction with $p<0.01\,p_{\max}$")
    ax[1, 1].set_title(r"(d) drained area (solid $\epsilon=1$, dashed $\epsilon=0$)",
                       fontsize=9)
    for a in (ax[0, 0], ax[0, 1], ax[1, 1]):
        a.set_xlabel(r"$\chi_\parallel/\chi_\perp$")
    ax[0, 0].legend(frameon=False, fontsize=7)
    cn = np.loadtxt(os.path.join(connect_dir, "connect.dat"))
    nturn = 100
    openl = (cn[:, 2] < nturn) | (cn[:, 3] < nturn)
    order = np.argsort(cn[:, 0])
    R, op = cn[order, 0], openl[order]
    for side in (R < 6.0, R > 6.0):
        r, o = R[side], op[side]
        edges = np.flatnonzero(np.diff(np.r_[0, o.astype(int), 0]))
        for a, b in zip(edges[::2], edges[1::2]):
            ax[1, 0].axvspan(r[a], r[b - 1], color="0.88", lw=0, zorder=0)
    ax[1, 0].axvspan(5.7, 6.3, color="0.6", lw=0, zorder=0)
    ax[1, 0].set_xlim(4.8, 7.2)
    ax[1, 0].set_xlabel(r"$R$ [m] ($Z=0$, $\varphi=0$)")
    ax[1, 0].set_ylabel(r"$p/p_{\max}$")
    ax[1, 0].set_title(r"(c) $\epsilon=0$, $\chi_\parallel/\chi_\perp=10^6$; "
                       "grey: lines reaching a wall", fontsize=9)
    ax[1, 0].legend(frameon=False, fontsize=7)
    fig.savefig(out, dpi=200)


def fit_exponent(W, Rh):
    """Least-squares slope of log Rh against log W and its standard error."""
    x, y = np.log(W), np.log(Rh)
    A = np.vstack([x, np.ones_like(x)]).T
    coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    n = len(x)
    if n < 3:
        return coef[0], np.nan, coef[1]
    sigma2 = np.sum((y - A @ coef) ** 2) / (n - 2)
    cov = sigma2 * np.linalg.inv(A.T @ A)
    return coef[0], np.sqrt(cov[0, 0]), coef[1]


def fitzpatrick(out, items):
    fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.6), constrained_layout=True)
    for j, item in enumerate(items):
        label, ds = item.split("=", 1)
        W, Rh = [], []
        for d in ds.split(","):
            if not os.path.exists(os.path.join(d, "scan.dat")):
                continue
            s = np.atleast_2d(np.loadtxt(os.path.join(d, "scan.dat")))
            if len(s) < 2:
                continue
            W.append(s[0, 2])
            Rh.append(half_ratio(s[:, 1], s[:, 4]))
            if j == len(items) - 1:
                ax[0].semilogx(s[:, 1], s[:, 4], lw=1,
                               label=rf"$W={s[0, 2]:.3f}$ m")
        W, Rh = np.array(W), np.array(Rh)
        ok = np.isfinite(Rh)
        kw = dict(color=OI[j % 7], marker=MARK[j % 7], ls="none")
        ax[1].loglog(W[ok], Rh[ok], **kw)
        if ok.sum() >= 2:
            slope, err, icpt = fit_exponent(W[ok], Rh[ok])
            w = np.geomspace(W[ok].min(), W[ok].max(), 20)
            ax[1].loglog(w, np.exp(icpt) * w ** slope, color=OI[j % 7], ls=LS[j % 7],
                         label=f"{label}: slope {slope:.2f} $\\pm$ {err:.2f}")
            print(f"{label}: slope {slope:.4f} +- {err:.4f} (n={ok.sum()}, "
                  f"W {W[ok].min():.4f}..{W[ok].max():.4f})")
        for wi, ri in zip(W, Rh):
            print(f"  {label} W {wi:.4f} ratio_half {ri:.4g}")
    ax[0].axhline(0.5, color="0.5", lw=0.8, ls=":")
    ax[0].set_xlabel(r"$\chi_\parallel/\chi_\perp$")
    ax[0].set_ylabel(r"$G$")
    ax[0].set_title(f"(a) G, finest level ({items[-1].split('=')[0]})", fontsize=9)
    ax[0].legend(frameon=False, fontsize=7, ncol=2)
    log_width_axis(ax[1])
    ax[1].set_xlabel(r"island width $W$ [m]")
    ax[1].set_ylabel(r"$\chi_\parallel/\chi_\perp$ at $G=1/2$")
    ax[1].set_title(r"(b) Fitzpatrick: slope $-4$ expected", fontsize=9)
    ax[1].legend(frameon=False, fontsize=7)
    fig.savefig(out, dpi=200)


if __name__ == "__main__":
    if sys.argv[1] == "poincare":
        poincare(sys.argv[2], sys.argv[3:])
    elif sys.argv[1] == "convergence":
        convergence(sys.argv[2], sys.argv[3], sys.argv[4:])
    elif sys.argv[1] == "fitzpatrick":
        fitzpatrick(sys.argv[2], sys.argv[3:])
    else:
        flattening(sys.argv[2], sys.argv[3:])
