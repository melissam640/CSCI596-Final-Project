#!/usr/bin/env python3
"""Transverse-stress maps of compressed finite pillars (numpy + matplotlib only).

Visualizes the boundary-induced transverse tension of Khodadad et al. (PRL 2026): with clamped
platens the pillar interior develops sigma_xx > 0 (a tensile "hotspot") that drives axial splitting.

Usage:
  python stress_map.py OUTDIR LABEL=stress_compress.lammpstrj [LABEL=...]
         [--strains 2,5,8,10,12,15] [--strain-per-frame 0.005] [--bin 4] [--slab 6]
         [--vmax 0.05] [--edge 8] [--grip 8] [--no-kscan]

Input: the per-atom stress dump of inputs/in.compress_finite (columns id element x y z
f_sav[1..6] = time-averaged stress x volume, atm*A^3, LAMMPS sign: positive = tension).
Frame 0 is the unloaded reference; frame k is at k x strain-per-frame engineering strain.

Per atom, the change from frame 0 is taken (delta sigma, removes the residual surface stress),
then binned in a slab: sigma_bin = n_bulk * mean(delta s_i), with n_bulk the bulk-cell number
density (each atom owns 1/n_bulk of volume). Sums and counts are 3x3 box-smoothed before
dividing. sigma_c (applied compressive stress) = -mean delta sigma_zz over the mid-height
cross-section. Maps are normalized by sigma_c, as in the paper's Fig. 2a.

Writes to OUTDIR:
  sxx_xz.png       delta sigma_xx / sigma_c in the slab |y - y_c| < slab (x across, z up),
                   black line = zero contour, hatched = opened bins (density < 50 % of frame 0)
  syy_yz.png       delta sigma_yy / sigma_c in the slab |x - x_c| < slab (y across, z up)
  sxx_xsection.png delta sigma_xx / sigma_c in the mid-height xy slab
  k_vs_strain.png  k vs strain per run. Primary (robust): k_core = max over z of the core-averaged
                   (r < R/2) delta sigma_xx / sigma_c, in z bins of 2*bin, grips and ends excluded.
                   Secondary (noisier): max over the smoothed slab bins.
  k_<label>.dat    strain sigma_c kcore_xx z_at kcore_yy z_at k_xx x_at z_at k_yy y_at z_at
                   open_bins nu_x nu_y
  stress_<label>_<strain>.xyz  extended XYZ with dsxx dsyy dszz (GPa, per-atom) for OVITO
Also prints the apparent Poisson ratios (mid-height lateral rms strain / axial strain, < 2 %).
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, MUTED, GRID = "#1f1f1e", "#6b6b68", "#e4e4e0"
COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]
N_BULK = 1408 / 16924.65          # atoms / A^3, bulk tobermorite-11A cell (info/tobe.xyz)
ATM2GPA = 1.01325e-4
OPEN_RATIO = 0.5


def read_frames(path, wanted=None):
    """Yield (k, box[3,2], elem[N], xyz[N,3], s[N,6]) for frames k in wanted (all if None)."""
    k = -1
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                return
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            k += 1
            f.readline()
            f.readline()
            n = int(f.readline())
            f.readline()
            box = np.array([[float(v) for v in f.readline().split()[:2]] for _ in range(3)])
            cols = f.readline().split()[2:]
            if wanted is not None and k not in wanted:
                for _ in range(n):
                    f.readline()
                if k >= max(wanted):
                    return
                continue
            rows = [f.readline().split() for _ in range(n)]
            c = {name: i for i, name in enumerate(cols)}
            ids = np.array([int(r[c["id"]]) for r in rows])
            assert np.all(np.diff(ids) > 0), "dump must be sorted by id"
            elem = np.array([r[c["element"]] for r in rows])
            num = [c[a] for a in "xyz"] + [c[f"f_sav[{i}]"] for i in range(1, 7)]
            a = np.array([[r[i] for i in num] for r in rows], float)
            yield k, box, elem, a[:, :3], a[:, 3:]


def smooth3(a):
    """3x3 box sum (zero padded)."""
    p = np.pad(a, 1)
    return sum(p[i:i + a.shape[0], j:j + a.shape[1]] for i in range(3) for j in range(3))


def bin2d(u, v, w, ue, ve):
    cnt, _, _ = np.histogram2d(u, v, bins=[ue, ve])
    tot = np.zeros((len(w), len(ue) - 1, len(ve) - 1))
    for i, wi in enumerate(w):
        tot[i], _, _ = np.histogram2d(u, v, bins=[ue, ve], weights=wi)
    return cnt, tot


class Run:
    def __init__(self, label, path, spf):
        self.label, self.path, self.spf = label, path, spf
        _, box, elem, x0, s0 = next(read_frames(path, {0}))
        self.elem, self.x0, self.s0 = elem, x0, s0
        self.xc, self.yc = x0[:, 0].mean(), x0[:, 1].mean()
        self.L0 = x0[:, 2].max() - x0[:, 2].min()
        r = np.hypot(x0[:, 0] - self.xc, x0[:, 1] - self.yc)
        self.R = np.percentile(r, 99.5)
        mid = np.abs(x0[:, 2] - self.zmid(x0)) < 10
        self.rms0 = np.sqrt(np.mean((x0[mid, :2] - [self.xc, self.yc]) ** 2, axis=0))
        self.mid0 = mid

    @staticmethod
    def zmid(x):
        return 0.5 * (x[:, 2].min() + x[:, 2].max())

    def fields(self, x, s):
        """Per-atom delta stress (GPa, as if each atom had volume 1/n_bulk) and sigma_c."""
        ds = (s - self.s0) * N_BULK * ATM2GPA
        zm = self.zmid(x)
        r = np.hypot(x[:, 0] - self.xc, x[:, 1] - self.yc)
        sec = (np.abs(x[:, 2] - zm) < 5) & (r < self.R)
        sig_c = -ds[sec, 2].mean()
        return ds, zm, sig_c

    def nu(self, x, strain):
        rms = np.sqrt(np.mean((x[self.mid0, :2] - [self.xc, self.yc]) ** 2, axis=0))
        return (rms / self.rms0 - 1) / max(strain, 1e-9)


def slab_map(run, x, x0, ds, zm, mode, a):
    """Binned, smoothed field and density ratio for one slab. Returns (ue, ve, field[3], ratio)."""
    half = run.R + 6
    ue = np.arange(-half, half + 1e-9, a.bin)
    if mode == "z":
        ve = ue
    else:
        hz = 0.5 * run.L0 + 8
        ve = np.arange(-hz, hz + 1e-9, a.bin)
    zm0 = run.zmid(x0)

    def coords(xx, zmid):
        dx, dy, dz = xx[:, 0] - run.xc, xx[:, 1] - run.yc, xx[:, 2] - zmid
        if mode == "y":      # slab normal y: (x, z)
            return np.abs(dy) < a.slab, dx, dz
        if mode == "x":      # slab normal x: (y, z)
            return np.abs(dx) < a.slab, dy, dz
        return np.abs(dz) < a.slab, dx, dy

    m, u, v = coords(x, zm)
    cnt, tot = bin2d(u[m], v[m], ds[m].T, ue, ve)
    # reference density mapped affinely onto the current (shortened) pillar length
    x0a = x0.copy()
    x0a[:, 2] = zm + (x0[:, 2] - zm0) * (np.ptp(x[:, 2]) / np.ptp(x0[:, 2]))
    m0, u0, v0 = coords(x0a, zm)
    cnt0, _ = bin2d(u0[m0], v0[m0], [], ue, ve)
    sc = smooth3(cnt)
    with np.errstate(invalid="ignore", divide="ignore"):
        f = np.array([smooth3(t) / sc for t in tot])
        f[:, sc < 9 * a.mincount] = np.nan
        ratio = smooth3(cnt) / smooth3(cnt0)
        ratio[smooth3(cnt0) < 9 * a.mincount] = np.nan
    return ue, ve, f, ratio


def interior(ue, ve, mode, run, a, L):
    uc, vc = 0.5 * (ue[1:] + ue[:-1]), 0.5 * (ve[1:] + ve[:-1])
    U, V = np.meshgrid(uc, vc, indexing="ij")
    if mode == "z":
        return np.hypot(U, V) < run.R - a.edge
    return (np.abs(U) < run.R - a.edge) & (np.abs(V) < 0.5 * L - a.grip - a.edge)


def write_xyz(fn, elem, x, ds, strain, box):
    with open(fn, "w") as f:
        L = box[:, 1] - box[:, 0]
        f.write(f"{len(elem)}\n")
        f.write(f'Lattice="{L[0]} 0 0 0 {L[1]} 0 0 0 {L[2]}" '
                f'Origin="{box[0, 0]} {box[1, 0]} {box[2, 0]}" pbc="F F F" '
                f'Properties=species:S:1:pos:R:3:dsxx:R:1:dsyy:R:1:dszz:R:1 strain={strain:.4f}\n')
        for e, p, d in zip(elem, x, ds):
            f.write(f"{e} {p[0]:.4f} {p[1]:.4f} {p[2]:.4f} {d[0]:.5g} {d[1]:.5g} {d[2]:.5g}\n")


def style(ax):
    ax.tick_params(colors=MUTED, labelsize=7)
    for sp in ax.spines.values():
        sp.set_color(GRID)
    ax.set_aspect("equal")


def map_figure(fn, panels, labels, strains, mode, comp, a):
    nr, nc = len(labels), len(strains)
    tall = mode != "z"
    fig, axes = plt.subplots(nr, nc, figsize=(1.9 * nc + 1.1, (3.6 if tall else 2.2) * nr + 0.4),
                             squeeze=False)
    im = None
    for row, lab in enumerate(labels):
        for col, eps in enumerate(strains):
            ax = axes[row, col]
            ue, ve, f, ratio, sig_c = panels[lab][eps][mode]
            g = f[comp] / sig_c if sig_c > 0 else f[comp] * np.nan
            im = ax.pcolormesh(ue, ve, g.T, cmap="RdBu_r", vmin=-a.vmax, vmax=a.vmax,
                               shading="flat", rasterized=True)
            uc, vc = 0.5 * (ue[1:] + ue[:-1]), 0.5 * (ve[1:] + ve[:-1])
            if np.isfinite(g).sum() > 4:
                ax.contour(uc, vc, np.ma.masked_invalid(g).T, levels=[0.0], colors=INK, linewidths=0.6)
            op = np.ma.masked_invalid(ratio).T < OPEN_RATIO
            if op.any():
                ax.contourf(uc, vc, op.astype(float), levels=[0.5, 1.5], colors="none", hatches=["////"])
            style(ax)
            if row == 0:
                ax.set_title(f"{eps:g} %  (σ_c {sig_c:.2f} GPa)", color=INK, fontsize=8)
            if col == 0:
                yl = "z − z_mid (Å)" if tall else "y − y_c (Å)"
                ax.set_ylabel(f"{lab}\n{yl}", color=MUTED, fontsize=8)
            if row == nr - 1:
                ax.set_xlabel({"y": "x − x_c (Å)", "x": "y − y_c (Å)", "z": "x − x_c (Å)"}[mode],
                              color=MUTED, fontsize=8)
    cb = fig.colorbar(im, ax=axes, shrink=0.8, pad=0.02)
    name = "xx" if comp == 0 else "yy"
    cb.set_label(f"Δσ_{name} / σ_c  (red = tension)", color=MUTED, fontsize=8)
    cb.ax.tick_params(colors=MUTED, labelsize=7)
    where = {"y": f"slab |y − y_c| < {a.slab:g} Å", "x": f"slab |x − x_c| < {a.slab:g} Å",
             "z": f"mid-height slab |z − z_mid| < {a.slab:g} Å"}[mode]
    fig.suptitle(f"Transverse stress Δσ_{name}/σ_c, {where}\nblack line = zero contour; "
                 f"hatched = density < {OPEN_RATIO:.0%} of the (affinely mapped) frame 0",
                 color=INK, fontsize=8, x=0.01, ha="left")
    fig.savefig(fn, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("runs", nargs="+", help="LABEL=stress_compress.lammpstrj")
    ap.add_argument("--strains", default="2,5,8,10,12,15", help="percent, for the maps")
    ap.add_argument("--strain-per-frame", type=float, default=0.005)
    ap.add_argument("--bin", type=float, default=4.0, help="bin size (A)")
    ap.add_argument("--slab", type=float, default=6.0, help="slab half-thickness (A)")
    ap.add_argument("--mincount", type=float, default=3.0, help="min mean atoms per bin")
    ap.add_argument("--vmax", type=float, default=0.05, help="color limit for delta sigma / sigma_c")
    ap.add_argument("--edge", type=float, default=8.0, help="exclude bins within this of the surface (A)")
    ap.add_argument("--grip", type=float, default=8.0, help="grip thickness excluded at each end (A)")
    ap.add_argument("--kmin-sc", type=float, default=1.0,
                    help="k_vs_strain.png: hide k where sigma_c < this (GPa); k diverges as sigma_c -> 0")
    ap.add_argument("--no-kscan", action="store_true", help="only the map strains, not every frame")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    strains = [float(s) for s in a.strains.split(",")]
    map_frames = {int(round(s / 100 / a.strain_per_frame)): s for s in strains}

    labels, panels, kdata = [], {}, {}
    for spec in a.runs:
        lab, path = spec.split("=", 1)
        labels.append(lab)
        run = Run(lab, path, a.strain_per_frame)
        print(f"{lab}: N={len(run.elem)} L0={run.L0:.1f} A R(99.5%)={run.R:.1f} A")
        panels[lab], rows = {}, []
        wanted = None if not a.no_kscan else set(map_frames) | {0}
        for k, box, elem, x, s in read_frames(path, wanted):
            if k == 0:
                continue
            eps = k * a.strain_per_frame
            ds, zm, sig_c = run.fields(x, s)
            ue, ve, f, ratio = slab_map(run, x, run.x0, ds, zm, "y", a)
            inn = interior(ue, ve, "y", run, a, x[:, 2].max() - x[:, 2].min())
            gx = np.where(inn, f[0], np.nan)
            ix = np.unravel_index(np.nanargmax(gx), gx.shape)
            ue2, ve2, f2, ratio2 = slab_map(run, x, run.x0, ds, zm, "x", a)
            gy = np.where(inn, f2[1], np.nan)
            iy = np.unravel_index(np.nanargmax(gy), gy.shape)
            nopen = int(np.nansum((ratio < OPEN_RATIO) & inn) + np.nansum((ratio2 < OPEN_RATIO) & inn))
            uc = 0.5 * (ue[1:] + ue[:-1])
            vc = 0.5 * (ve[1:] + ve[:-1])
            nux, nuy = run.nu(x, eps)
            # core-averaged profile along z (many atoms per bin -> robust k)
            Lc = x[:, 2].max() - x[:, 2].min()
            rr = np.hypot(x[:, 0] - run.xc, x[:, 1] - run.yc)
            core = rr < 0.5 * run.R
            zb = np.arange(-0.5 * Lc, 0.5 * Lc + 1e-9, 2 * a.bin)
            zz = x[core, 2] - zm
            cn, _ = np.histogram(zz, zb)
            zcen = 0.5 * (zb[1:] + zb[:-1])
            okz = (cn > 20) & (np.abs(zcen) < 0.5 * Lc - a.grip - a.edge)
            prof = []
            for comp in (0, 1):
                t, _ = np.histogram(zz, zb, weights=ds[core, comp])
                pr = np.where(okz, t / np.maximum(cn, 1), np.nan)
                j = np.nanargmax(pr)
                prof += [pr[j] / sig_c, zcen[j]]
            rows.append((eps, sig_c, *prof, gx[ix] / sig_c, uc[ix[0]], vc[ix[1]],
                         gy[iy] / sig_c, uc[iy[0]], vc[iy[1]], nopen, nux, nuy))
            if k in map_frames:
                pe = map_frames[k]
                fz = slab_map(run, x, run.x0, ds, zm, "z", a)
                panels[lab][pe] = {"y": (ue, ve, f, ratio, sig_c), "x": (ue2, ve2, f2, ratio2, sig_c),
                                   "z": (*fz, sig_c)}
                write_xyz(os.path.join(a.outdir, f"stress_{lab}_{pe:g}.xyz"), elem, x, ds, eps, box)
        missing = set(map_frames.values()) - set(panels[lab])
        if missing:
            raise SystemExit(f"{path}: strains {sorted(missing)} % not in the dump")
        rows = np.array(rows)
        kdata[lab] = rows
        np.savetxt(os.path.join(a.outdir, f"k_{lab}.dat"), rows, fmt="%.5g",
                   header="strain sigma_c_GPa kcore_xx z_at kcore_yy z_at k_xx x_at z_at "
                          "k_yy y_at z_at open_bins nu_x nu_y")
        lo = rows[:, 0] <= 0.02
        if lo.any():
            print(f"{lab}: apparent Poisson ratio (mid-height, strain <= 2 %): "
                  f"nu_x {np.mean(rows[lo, 13]):.3f}  nu_y {np.mean(rows[lo, 14]):.3f}")
        for r in rows:
            if np.any(np.isclose(100 * r[0], strains)):
                print(f"  {100 * r[0]:5.1f} %  sigma_c {r[1]:6.2f} GPa  kcore_xx {100 * r[2]:6.2f} % "
                      f"(z {r[3]:+.0f})  kcore_yy {100 * r[4]:6.2f} %  | slab-bin max k_xx "
                      f"{100 * r[6]:6.2f} % at (x {r[7]:+.0f}, z {r[8]:+.0f})  open bins {int(r[12])}")

    shown = sorted(set(map_frames.values()))
    map_figure(os.path.join(a.outdir, "sxx_xz.png"), panels, labels, shown, "y", 0, a)
    map_figure(os.path.join(a.outdir, "syy_yz.png"), panels, labels, shown, "x", 1, a)
    map_figure(os.path.join(a.outdir, "sxx_xsection.png"), panels, labels, shown, "z", 0, a)

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(6.4, 6.0), sharex=True)
    for lab, col in zip(labels, COLORS):
        r = kdata[lab]
        kx, ky = np.where(r[:, 1] >= a.kmin_sc, r[:, 2], np.nan), np.where(r[:, 1] >= a.kmin_sc, r[:, 4], np.nan)
        ax.plot(100 * r[:, 0], 100 * kx, "-", color=col, lw=1.6, label=f"{lab} core k_xx (normal to layers)")
        ax.plot(100 * r[:, 0], 100 * ky, ":", color=col, lw=1.2, label=f"{lab} core k_yy")
        ax2.plot(100 * r[:, 0], r[:, 1], "-", color=col, lw=1.6, label=lab)
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.set_ylabel(f"k = max core Δσ / σ_c (%), σ_c ≥ {a.kmin_sc:g} GPa", color=MUTED)
    ax.set_title("Boundary-induced transverse tension", color=INK, loc="left")
    ax2.set_ylabel("σ_c, mid-height (GPa)", color=MUTED)
    ax2.set_xlabel("engineering strain (%)", color=MUTED)
    for axx in (ax, ax2):
        axx.grid(color=GRID, lw=0.8)
        axx.tick_params(colors=MUTED)
        for sp in axx.spines.values():
            sp.set_color(GRID)
        axx.legend(frameon=False, labelcolor=INK, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "k_vs_strain.png"), dpi=150)


if __name__ == "__main__":
    main()
