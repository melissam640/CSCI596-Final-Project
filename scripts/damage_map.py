#!/usr/bin/env python3
"""Shear-damage maps of compressed pillars from LAMMPS custom dumps (numpy only).

Usage:
  python damage_map.py OUTDIR LABEL=traj_compress.lammpstrj [LABEL=...]
         [--strains 5,8,12,14,20,30] [--rc 4.0] [--crack LABEL]
         [--crack-a 10] [--crack-tilt 45]

Reference state is frame 0 of each dump (post-NPT). Frame k is at k x 0.1 % strain
(dump every 400 steps x 0.25 fs, erate 1e-5/fs). Per atom and per requested strain:
  eta    von Mises atomic shear strain (Shimizu, Ogata & Li 2007; = OVITO Atomic Strain)
  d2min  non-affine residual (Falk & Langer 1998), A^2
  broken 1 = Si that lost a frame-0 Si-O bond (d0 < 2.0 A, now > 2.3 A),
         2 = Ca whose Ca-O coordination (< 2.9 A) changed by >= 2
Writes to OUTDIR: shear_slices.png (edge-on y-slab, x across, z up), shear_xsection.png
(mid-height xy-slab), damage_vs_strain.png, damage_summary.dat, and one extended-XYZ file
per run and strain (damage_<label>_<strain>.xyz, for OVITO).
Stress curves for damage_vs_strain.png are read from results/stress_strain_<label>.dat.
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, MUTED, GRID = "#1f1f1e", "#6b6b68", "#e4e4e0"
COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]   # same categorical order as plot_compare.py
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRAIN_PER_FRAME = 1.0e-3
ETA_DAMAGED = 0.2
SIO_CUT, SIO_BROKEN, CAO_CUT = 2.0, 2.3, 2.9   # Si-O cutoff as in build_pillar.py


def read_frames(path, wanted):
    """Return {frame_index: (box[3,2], elem[N], xyz[N,3])} for the wanted frame indices."""
    out, k = {}, -1
    with open(path) as f:
        while len(out) < len(wanted):
            line = f.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            k += 1
            f.readline()
            f.readline()
            n = int(f.readline())
            f.readline()
            box = np.array([[float(v) for v in f.readline().split()[:2]] for _ in range(3)])
            cols = f.readline().split()[2:]
            if k not in wanted:
                for _ in range(n):
                    f.readline()
                continue
            rows = [f.readline().split() for _ in range(n)]
            c = {name: i for i, name in enumerate(cols)}
            ids = np.array([int(r[c["id"]]) for r in rows])
            assert np.all(np.diff(ids) > 0), "dump must be sorted by id"
            elem = np.array([r[c["element"]] for r in rows])
            xyz = np.array([[float(r[c[a]]) for a in "xyz"] for r in rows])
            out[k] = (box, elem, xyz)
    missing = sorted(set(wanted) - set(out))
    if missing:
        raise SystemExit(f"{path}: frames {missing} not found (only {k + 1} frames)")
    return out


def min_image_z(d, lz):
    d[..., 2] -= lz * np.round(d[..., 2] / lz)
    return d


def pairs_within(xa, xb, lz, rc, same=False, chunk=256):
    """All (i, j) with |xa[i] - xb[j]| < rc, z periodic. Brute force in chunks."""
    I, J = [], []
    for s in range(0, len(xa), chunk):
        d = min_image_z(xb[None, :, :] - xa[s:s + chunk, None, :], lz)
        r2 = np.einsum("ijk,ijk->ij", d, d)
        m = r2 < rc * rc
        if same:
            m[np.arange(m.shape[0]), np.arange(s, s + m.shape[0])] = False
        i, j = np.nonzero(m)
        I.append(i + s)
        J.append(j)
    return np.concatenate(I), np.concatenate(J)


def atomic_strain(x0, lz0, x1, lz1, I, J):
    """Per-atom eta and d2min, frame 0 -> frame 1, over the frame-0 neighbor pairs (I, J)."""
    n = len(x0)
    d0 = min_image_z(x0[J] - x0[I], lz0)
    d1 = min_image_z(x1[J] - x1[I], lz1)
    V = np.zeros((n, 3, 3))
    W = np.zeros((n, 3, 3))
    np.add.at(V, I, d0[:, :, None] * d0[:, None, :])
    np.add.at(W, I, d0[:, :, None] * d1[:, None, :])
    nnb = np.bincount(I, minlength=n)
    ok = (nnb >= 4) & (np.linalg.cond(V) < 1e4)
    X = np.full((n, 3, 3), np.nan)              # row convention: d1 ~ d0 @ X, F = X^T
    X[ok] = np.linalg.solve(V[ok], W[ok])
    FtF = np.einsum("nij,nkj->nik", X, X)       # (X X^T) = F^T F in column convention
    E = 0.5 * (FtF - np.eye(3))
    eta = np.sqrt(E[:, 0, 1] ** 2 + E[:, 0, 2] ** 2 + E[:, 1, 2] ** 2 +
                  ((E[:, 0, 0] - E[:, 1, 1]) ** 2 + (E[:, 1, 1] - E[:, 2, 2]) ** 2 +
                   (E[:, 2, 2] - E[:, 0, 0]) ** 2) / 6.0)
    res = d1 - np.einsum("pi,pij->pj", d0, np.nan_to_num(X[I]))
    d2min = np.bincount(I, weights=np.einsum("pi,pi->p", res, res), minlength=n)
    d2min[~ok] = np.nan
    return eta, d2min, X[:, 2, 2]


def coordination(x, elem, a, b, lz, rc):
    ia, ib = np.nonzero(elem == a)[0], np.nonzero(elem == b)[0]
    i, j = pairs_within(x[ia], x[ib], lz, rc)
    return ia, ia[i], ib[j]


def analyse(path, strains, rc):
    wanted = [int(round(s * 10)) for s in strains]
    data = read_frames(path, {0} | set(wanted))
    box0, elem, x0 = data[0]
    lz0 = box0[2, 1] - box0[2, 0]
    I, J = pairs_within(x0, x0, lz0, rc, same=True)
    _, si0, o0 = coordination(x0, elem, "Si", "O", lz0, SIO_CUT)
    ca_idx, ca_i, _ = coordination(x0, elem, "Ca", "O", lz0, CAO_CUT)
    cn_ca0 = np.bincount(ca_i, minlength=len(x0))[ca_idx]
    res = []
    for k in wanted:
        box, _, x = data[k]
        lz = box[2, 1] - box[2, 0]
        eta, d2min, fzz = atomic_strain(x0, lz0, x, lz, I, J)
        dsio = np.linalg.norm(min_image_z(x[o0] - x[si0], lz), axis=1)
        broken = np.zeros(len(x), int)
        broken[np.unique(si0[dsio > SIO_BROKEN])] = 1
        _, ca_i1, _ = coordination(x, elem, "Ca", "O", lz, CAO_CUT)
        cn_ca = np.bincount(ca_i1, minlength=len(x))[ca_idx]
        broken[ca_idx[np.abs(cn_ca - cn_ca0) >= 2]] = 2
        res.append(dict(frame=k, strain=k * STRAIN_PER_FRAME, box=box, x=x, eta=eta,
                        d2min=d2min, fzz=fzz, broken=broken,
                        n_sio=int((dsio > SIO_BROKEN).sum())))
    return elem, box0, x0, res


def write_xyz(fn, elem, r):
    box = r["box"]
    L = box[:, 1] - box[:, 0]
    with open(fn, "w") as f:
        f.write(f"{len(elem)}\n")
        f.write(f'Lattice="{L[0]} 0 0 0 {L[1]} 0 0 0 {L[2]}" '
                f'Origin="{box[0, 0]} {box[1, 0]} {box[2, 0]}" pbc="F F T" '
                f'Properties=species:S:1:pos:R:3:eta:R:1:d2min:R:1:broken:I:1 '
                f'strain={r["strain"]:.4f}\n')
        for e, p, a, d, b in zip(elem, r["x"], np.nan_to_num(r["eta"]),
                                 np.nan_to_num(r["d2min"]), r["broken"]):
            f.write(f"{e} {p[0]:.4f} {p[1]:.4f} {p[2]:.4f} {a:.5f} {d:.4f} {b}\n")


def style(ax):
    ax.tick_params(colors=MUTED, labelsize=7)
    for sp in ax.spines.values():
        sp.set_color(GRID)
    ax.set_aspect("equal")


def slab_figure(fn, runs, labels, strains, mode, crack, crack_a, crack_tilt):
    nr, nc = len(runs), len(strains)
    fig, axes = plt.subplots(nr, nc, figsize=(2.1 * nc + 0.9, 3.4 * nr if mode == "y" else 2.4 * nr),
                             squeeze=False)
    sc = None
    for row, (lab, (elem, box0, x0, res)) in enumerate(zip(labels, runs)):
        for col, r in enumerate(res):
            ax = axes[row, col]
            x, box = r["x"], r["box"]
            xc, yc = x[:, 0].mean(), x[:, 1].mean()
            lz = box[2, 1] - box[2, 0]
            zc = box[2, 0] + 0.5 * lz
            if mode == "y":
                m = np.abs(x[:, 1] - yc) < 3.0
                u, v = x[m, 0] - xc, x[m, 2]
            else:
                dz = min_image_z(x - [0, 0, zc], lz)[:, 2]
                m = np.abs(dz) < 3.0
                u, v = x[m, 0] - xc, x[m, 1] - yc
            order = np.argsort(np.nan_to_num(r["eta"][m]))
            sc = ax.scatter(u[order], v[order], c=r["eta"][m][order], s=2.0, cmap="magma_r",
                            vmin=0, vmax=0.3, linewidths=0, rasterized=True)
            b = r["broken"][m] > 0
            ax.scatter(u[b], v[b], s=9, facecolors="none", edgecolors="#2a78d6", linewidths=0.5)
            if mode == "y" and lab == crack:
                # frame-0 crack trace mapped affinely with the box in z
                lz0 = box0[2, 1] - box0[2, 0]
                t = np.radians(crack_tilt)
                s = np.array([-crack_a, crack_a])
                zs0 = 0.5 * lz0 - s * np.sin(t)           # relative to zlo
                ax.plot(s * np.cos(t), box[2, 0] + zs0 * lz / lz0, "--", color="#1baf7a", lw=0.9)
            style(ax)
            if row == 0:
                ax.set_title(f"{100 * r['strain']:.0f} % strain", color=INK, fontsize=9)
            if col == 0:
                ax.set_ylabel(f"{lab}\n" + ("z (Å)" if mode == "y" else "y − y_c (Å)"),
                              color=MUTED, fontsize=8)
            if row == nr - 1:
                ax.set_xlabel("x − x_c (Å)", color=MUTED, fontsize=8)
    cb = fig.colorbar(sc, ax=axes, shrink=0.8, pad=0.02)
    cb.set_label("von Mises shear strain η", color=MUTED, fontsize=8)
    cb.ax.tick_params(colors=MUTED, labelsize=7)
    what = "edge-on slab |y − y_c| < 3 Å" if mode == "y" else "mid-height slab |z − z_c| < 3 Å"
    fig.suptitle(f"Atomic shear strain vs frame 0, {what}; blue rings = broken Si–O / Ca–O change",
                 color=INK, fontsize=9, x=0.01, ha="left")
    fig.savefig(fn, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("runs", nargs="+", help="LABEL=traj_compress.lammpstrj")
    ap.add_argument("--strains", default="5,8,12,14,20,30", help="percent, comma separated")
    ap.add_argument("--rc", type=float, default=4.0, help="neighbor cutoff for atomic strain (Å)")
    ap.add_argument("--crack", default=None, help="label of the run with the penny crack")
    ap.add_argument("--crack-a", type=float, default=10.0)
    ap.add_argument("--crack-tilt", type=float, default=45.0)
    a = ap.parse_args()
    strains = [float(s) for s in a.strains.split(",")]
    os.makedirs(a.outdir, exist_ok=True)
    labels, runs = [], []
    for spec in a.runs:
        lab, path = spec.split("=", 1)
        labels.append(lab)
        runs.append(analyse(path, strains, a.rc))

    lines = ["# label strain_pct median_eta core_mean_eta core_p99_eta shell_mean_eta "
             "frac_eta_gt_0.2 core_median_Fzz broken_SiO_bonds Si_broken Ca_cn_changed"]
    print(f"{'run':8s} {'ε%':>5s} {'η med':>6s} {'core η':>7s} {'core p99':>8s} {'shell η':>7s} "
          f"{'f(η>.2)':>7s} {'Fzz':>6s} {'SiO brk':>7s} {'Ca Δcn':>6s}")
    for lab, (elem, box0, x0, res) in zip(labels, runs):
        for r in res:
            x = r["x"]
            rr = np.hypot(x[:, 0] - x[:, 0].mean(), x[:, 1] - x[:, 1].mean())
            core, eta = rr < 10.0, r["eta"]
            row = (np.nanmedian(eta), np.nanmean(eta[core]), np.nanpercentile(eta[core], 99),
                   np.nanmean(eta[~core]), np.nanmean(eta > ETA_DAMAGED),
                   np.nanmedian(r["fzz"][core]), r["n_sio"],
                   int((r["broken"] == 1).sum()), int((r["broken"] == 2).sum()))
            r["frac"] = row[4]
            print(f"{lab:8s} {100 * r['strain']:5.1f} {row[0]:6.3f} {row[1]:7.3f} {row[2]:8.3f} "
                  f"{row[3]:7.3f} {row[4]:7.3f} {row[5]:6.3f} {row[6]:7d} {row[8]:6d}")
            lines.append(f"{lab} {100 * r['strain']:.1f} " + " ".join(f"{v:.5g}" for v in row))
            write_xyz(os.path.join(a.outdir, f"damage_{lab}_{100 * r['strain']:.0f}.xyz"), elem, r)
    open(os.path.join(a.outdir, "damage_summary.dat"), "w").write("\n".join(lines) + "\n")

    shown = [r["strain"] for r in runs[0][3]]
    slab_figure(os.path.join(a.outdir, "shear_slices.png"), runs, labels, shown, "y",
                a.crack, a.crack_a, a.crack_tilt)
    slab_figure(os.path.join(a.outdir, "shear_xsection.png"), runs, labels, shown, "z",
                a.crack, a.crack_a, a.crack_tilt)

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax2 = ax.twinx()
    for lab, (_, _, _, res), col in zip(labels, runs, COLORS):
        e = [100 * r["strain"] for r in res]
        ax.plot(e, [100 * r["frac"] for r in res], "o-", color=col, lw=1.6, label=f"{lab} η > {ETA_DAMAGED}")
        ax.plot(e, [100 * r["n_sio"] / len(r["x"]) for r in res], "s:", color=col, lw=1.0, ms=4,
                label=f"{lab} broken Si–O (per 100 atoms)")
        ss = os.path.join(REPO, "results", f"stress_strain_{lab}.dat")
        if os.path.exists(ss):
            d = np.loadtxt(ss)
            ax2.plot(100 * d[:, 1], d[:, 2], color=col, lw=0.8, alpha=0.45)
    ax.set_xlabel("engineering strain (%)", color=MUTED)
    ax.set_ylabel("damaged atoms (%)", color=MUTED)
    ax2.set_ylabel("compressive stress (GPa), thin lines", color=MUTED)
    ax.set_title("Shear damage vs strain", color=INK, loc="left")
    ax.grid(color=GRID, lw=0.8)
    for axx in (ax, ax2):
        axx.tick_params(colors=MUTED)
        for sp in axx.spines.values():
            sp.set_color(GRID)
    ax.legend(frameon=False, labelcolor=INK, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "damage_vs_strain.png"), dpi=150)


if __name__ == "__main__":
    main()
