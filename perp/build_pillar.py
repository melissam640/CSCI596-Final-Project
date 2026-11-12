#!/usr/bin/env python3
"""Carve a cylindrical tobermorite-11A pillar with crystal c along x and the pillar axis along z.

Usage: python build_pillar.py [--R 20.0] [--nc 4] [--vac 40] [--out DIR] [--dry]
       [--crack-a 10 --crack-h 4 [--crack-tilt 90 --crack-az 0 --crack-zc 0.5]]

Reads the rotated bulk cell written by perp/rotate_cell.py (rows c' = x, b' in yz, a' = +z).
The pillar axis is crystal a (periodic length nc*|a|), the layers are perpendicular to x.
The axis is placed on the Ca-sheet midplane nearest the cell centre (in x), so a penny crack
with its normal along x (tilt 90, az 0, the defaults) lies in that Ca sheet.

Same unit-keeping and charge-balance rules as scripts/build_pillar.py: whole SiO4 tetrahedra
and whole H2O/OH units are kept; the formal charge is balanced by protonating the outermost
non-bridging Si-O (if negative) or deleting the outermost Ca (if positive). Whole units whose
Si/Ca/water-O centre lies in the crack disk are removed before charge balance.
Writes pillar.data (LAMMPS atom_style charge) and pillar0.xyz, and prints R0 for -var R0.
"""
import argparse
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TYPES = ["Ca", "Si", "O", "H"]
MASS = {"Ca": 40.078, "Si": 28.0855, "O": 15.999, "H": 1.008}
QFORMAL = {"Ca": 2, "Si": 4, "O": -2, "H": 1}
R_SIO, R_OH = 2.0, 1.2


def read_frame1(path):
    with open(path) as f:
        n = int(f.readline())
        hdr = f.readline()
        lat = np.array(hdr.split('Lattice="')[1].split('"')[0].split(), float).reshape(3, 3)
        sp, xyz = [], []
        for _ in range(n):
            w = f.readline().split()
            sp.append(w[0])
            xyz.append([float(v) for v in w[1:4]])
    return np.array(sp), np.array(xyz), lat  # rows c' (x), b' (yz), a (+z)


def ca_sheet_x(sp0, x0, px, target, gap=2.5):
    """x of the Ca-sheet midplane nearest target (sheets are planes x = const, period px)."""
    v = np.sort(x0[sp0 == "Ca", 0] % px)
    groups = np.split(v, np.where(np.diff(v) > gap)[0] + 1)
    if len(groups) > 1 and v[0] + px - v[-1] <= gap:   # cluster wraps through x = 0
        groups[0] = np.concatenate([groups[-1] - px, groups[0]])
        groups.pop()
    mids = np.sort(np.array([g.mean() for g in groups]) % px)
    cand = (mids[None, :] + px * np.arange(-1, 2)[:, None]).ravel()
    xs = cand[np.argmin(np.abs(cand - target))]
    dsi = (x0[sp0 == "Si", 0] - xs + px / 2) % px - px / 2
    print(f"Ca-sheet midplanes in cell x (A): {' '.join(f'{m:.2f}' for m in mids)}; "
          f"pillar axis at x = {xs:.2f} A (cell centre {target:.2f}); "
          f"nearest Si planes {dsi[dsi < 0].max():+.2f} / {dsi[dsi > 0].min():+.2f} A")
    return xs


def build(args):
    sp0, x0, lat = read_frame1(args.xyz)
    u, v, w = lat                       # u = c' (x), v = b' (yz), w = a (+z, pillar axis)
    assert abs(w[0]) < 1e-6 and abs(w[1]) < 1e-6 and w[2] > 0, "third lattice row must be +z"
    assert abs(u[1]) < 1e-6 and abs(u[2]) < 1e-6, "first lattice row must be along x"
    frac = np.linalg.solve(lat.T, x0.T).T % 1.0
    x0 = frac @ lat
    Lz = args.nc * w[2]

    # supercell: i,j in [-m,m] along u,v; k in [0,nc) along w; then wrap z into [0, Lz).
    # v has a z component, but each lattice site mod nc*w still appears exactly once.
    m = int(np.ceil((args.R + 15) / min(u[0], abs(v[1])))) + 1
    shifts = np.array([i * u + j * v + k * w for i in range(-m, m + 1)
                       for j in range(-m, m + 1) for k in range(args.nc)])
    X = (x0[None, :, :] + shifts[:, None, :]).reshape(-1, 3)
    S = np.tile(sp0, len(shifts))
    X[:, 2] %= Lz
    # axis: x on the Ca sheet nearest the cell centre, y at the cell centre
    ctr = 0.5 * (u + v)
    X[:, 0] -= ca_sheet_x(sp0, x0, u[0], ctr[0])
    X[:, 1] -= ctr[1]
    r = np.hypot(X[:, 0], X[:, 1])
    crop = r < args.R + 6.0
    X, S, r = X[crop], S[crop], r[crop]

    # optional penny-shaped crack: disk of radius a, opening h, normal at tilt from z, az in xy
    # (tilt 90, az 0 = normal along x = layer normal: the disk lies in the axis Ca sheet)
    incrack = np.zeros(len(S), bool)
    if args.crack_a > 0:
        th, ph = np.radians(args.crack_tilt), np.radians(args.crack_az)
        n = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
        assert args.crack_a < args.R - 5, "crack must stay >= 5 A inside the pillar surface"
        assert args.crack_a * np.sin(th) + args.crack_h / 2 < Lz / 2, "crack overlaps its z image"
        d = X - np.array([0.0, 0.0, args.crack_zc * Lz])
        d[:, 2] -= Lz * np.round(d[:, 2] / Lz)
        dn = d @ n
        dp = np.linalg.norm(d - dn[:, None] * n, axis=1)
        incrack = (np.abs(dn) < args.crack_h / 2) & (dp < args.crack_a)

    def pair_d(i_idx, j_idx):
        d = X[j_idx][None, :, :] - X[i_idx][:, None, :]
        d[..., 2] -= Lz * np.round(d[..., 2] / Lz)
        return np.linalg.norm(d, axis=2)

    iSi, iO, iH, iCa = (np.where(S == e)[0] for e in ("Si", "O", "H", "Ca"))
    dSiO = pair_d(iSi, iO) < R_SIO           # [nSi, nO]
    dOH = pair_d(iO, iH)                      # [nO, nH]
    h_owner = iO[np.argmin(dOH, axis=0)]      # O index of each H
    h_ok = dOH.min(axis=0) < R_OH               # crop-edge H may lose their O
    assert np.all(h_ok | (r[iH] > args.R + 2)), "H without O partner inside pillar"
    o_nsi = dSiO.sum(axis=0)                  # Si neighbors per O

    def select(excl):
        """Whole-unit selection inside R, skipping Si/Ca/free-O centers flagged in excl."""
        keep = np.zeros(len(S), bool)
        si_keep = (r[iSi] < args.R) & ~excl[iSi]
        keep[iSi[si_keep]] = True
        keep[iCa[(r[iCa] < args.R) & ~excl[iCa]]] = True
        o_from_si = dSiO[si_keep].any(axis=0)
        o_free = (o_nsi == 0) & (r[iO] < args.R) & ~excl[iO]   # water / Ca-OH oxygens
        keep[iO[o_from_si | o_free]] = True
        keep[iH] = keep[h_owner] & h_ok
        return keep, si_keep

    keep, si_keep = select(incrack)
    gone = select(np.zeros(len(S), bool))[0] & ~keep
    ncrack = {e: int((gone & (S == e)).sum()) for e in TYPES}

    # charge balance
    q = sum(QFORMAL[e] for e in S[keep])
    newH = []
    if q < 0:
        # non-bridging O: bonded to exactly one kept Si and no H
        nkeptSi = dSiO[si_keep].sum(axis=0)
        nH = np.bincount(np.searchsorted(iO, h_owner[h_ok]), minlength=len(iO))
        cand = np.where(keep[iO] & (nkeptSi == 1) & (nH == 0))[0]
        cand = cand[np.argsort(-r[iO[cand]])]
        if len(cand) < -q:
            raise SystemExit(f"not enough non-bridging O to protonate (q={q})")
        sik = iSi[si_keep]
        for oc in cand[:-q]:
            o = iO[oc]
            s = sik[np.where(dSiO[si_keep][:, oc])[0][0]]
            vv = X[o] - X[s]
            vv[2] -= Lz * np.round(vv[2] / Lz)
            vv /= np.linalg.norm(vv)
            rad = np.array([X[o, 0], X[o, 1], 0.0]) / max(r[o], 1e-6)
            uu = vv + rad
            uu /= np.linalg.norm(uu)
            newH.append(X[o] + 0.97 * uu)
    elif q > 0:
        if q % 2:
            raise SystemExit(f"odd positive charge q={q}; adjust R")
        ca = np.where(keep & (S == "Ca"))[0]
        ca = ca[np.argsort(-r[ca])][: q // 2]
        keep[ca] = False

    Xk, Sk = X[keep], S[keep]
    if newH:
        Xk = np.vstack([Xk, np.array(newH)])
        Sk = np.concatenate([Sk, np.array(["H"] * len(newH))])
    Xk[:, 2] %= Lz
    ndens = len(sp0) / abs(np.linalg.det(lat))
    return Xk, Sk, Lz, q, len(newH), ncrack, ndens


def min_dist(X, Lz, cut=1.5):
    best = np.inf
    for i in range(0, len(X), 500):
        d = X[None, :, :] - X[i:i + 500, None, :]
        d[..., 2] -= Lz * np.round(d[..., 2] / Lz)
        dd = np.linalg.norm(d, axis=2)
        idx = np.arange(i, min(i + 500, len(X)))
        dd[np.arange(len(idx)), idx] = np.inf
        best = min(best, dd.min())
    return best


def crack_label(args):
    if args.crack_a <= 0:
        return ""
    return (f" crack a={args.crack_a} h={args.crack_h} tilt={args.crack_tilt}"
            f" az={args.crack_az} zc={args.crack_zc:.4g}")


def write(args, X, S, Lz):
    half = args.R + args.vac
    order = np.argsort([TYPES.index(e) for e in S], kind="stable")
    X, S = X[order], S[order]
    # box [0, 2*half] in x,y with the pillar axis at its center (viewers assume a 0-origin cell)
    X = X + np.array([half, half, 0.0])
    with open(f"{args.out}/pillar.data", "w") as f:
        f.write(f"tobermorite-11A pillar R={args.R} nc={args.nc} (axis = a = z, c = x){crack_label(args)}\n\n")
        f.write(f"{len(S)} atoms\n{len(TYPES)} atom types\n\n")
        f.write(f"0.0 {2*half:.6f} xlo xhi\n0.0 {2*half:.6f} ylo yhi\n")
        f.write(f"0.0 {Lz:.6f} zlo zhi\n\nMasses\n\n")
        for t, e in enumerate(TYPES, 1):
            f.write(f"{t} {MASS[e]}  # {e}\n")
        f.write("\nAtoms # charge\n\n")
        for i, (e, x) in enumerate(zip(S, X), 1):
            f.write(f"{i} {TYPES.index(e) + 1} 0.0 {x[0]:.6f} {x[1]:.6f} {x[2]:.6f}\n")
    with open(f"{args.out}/pillar0.xyz", "w") as f:
        f.write(f"{len(S)}\nLattice=\"{2*half} 0 0 0 {2*half} 0 0 0 {Lz}\" "
                f"Properties=species:S:1:pos:R:3\n")
        for e, x in zip(S, X):
            f.write(f"{e} {x[0]:.6f} {x[1]:.6f} {x[2]:.6f}\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--xyz", default=os.path.join(HERE, "tobe_rot.xyz"))
    p.add_argument("--R", type=float, default=20.0)
    p.add_argument("--nc", type=int, default=4, help="unit cells along the pillar axis (crystal a)")
    p.add_argument("--vac", type=float, default=40.0)
    p.add_argument("--out", default=".")
    p.add_argument("--dry", action="store_true", help="report only, no files")
    p.add_argument("--crack-a", type=float, default=0.0, help="penny crack radius (A); 0 = none")
    p.add_argument("--crack-h", type=float, default=4.0, help="crack opening thickness (A)")
    p.add_argument("--crack-tilt", type=float, default=90.0,
                   help="crack normal vs z (deg); 90 = normal along x, in the Ca sheet")
    p.add_argument("--crack-az", type=float, default=0.0, help="azimuth of the normal in xy (deg)")
    p.add_argument("--crack-zc", type=float, default=0.5, help="crack center, fraction of Lz")
    args = p.parse_args()
    X, S, Lz, q, nadd, ncrack, ndens = build(args)
    els, cnt = np.unique(S, return_counts=True)
    qf = sum(QFORMAL[e] * n for e, n in zip(els, cnt))
    rr = np.hypot(X[:, 0], X[:, 1])
    R0 = np.sqrt(len(S) / (ndens * np.pi * Lz))
    print(f"R={args.R} Lz={Lz:.4f} N={len(S)} " + " ".join(f"{e}{n}" for e, n in zip(els, cnt)))
    print(f"formal charge before balance={q}, H added={nadd}, final={qf}, r_max={rr.max():.2f}")
    print(f"R0 (bulk-density equivalent radius, pass as -var R0) = {R0:.2f} A")
    if args.crack_a > 0:
        print(f"penny{crack_label(args)}: {sum(ncrack.values())} atoms removed "
              f"({' '.join(f'{e}{n}' for e, n in ncrack.items())})")
    if not args.dry:
        print(f"min interatomic distance={min_dist(X, Lz):.3f} A")
        write(args, X, S, Lz)


if __name__ == "__main__":
    main()
