#!/usr/bin/env python3
"""Carve a cylindrical tobermorite-11A pillar (axis = crystal c = z) from info/tobe.xyz.

Usage: python build_pillar.py [--R 20.2] [--nc 4] [--vac 40] [--target 10000] [--out DIR]
       [--crack-a 10 --crack-h 4 --crack-tilt 45 --crack-az 0 --crack-zc 0.5]

Keeps whole SiO4 tetrahedra and whole H2O/OH units, then charge-balances the
formal charges (Ca+2 Si+4 O-2 H+1) by protonating outermost non-bridging Si-O
(if negative) or deleting outermost Ca (if positive).
Optional embedded penny-shaped crack (--crack-a > 0): whole units whose Si/Ca/
water-O center lies in the disk are removed before charge balance.
Writes pillar.data (LAMMPS atom_style charge) and pillar0.xyz.
"""
import argparse
import numpy as np

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
    return np.array(sp), np.array(xyz), lat  # lat rows = a, b, c


def build(args):
    sp0, x0, lat = read_frame1(args.xyz)
    a, b, c = lat
    assert abs(c[0]) < 1e-6 and abs(c[1]) < 1e-6 and abs(a[2]) < 1e-6 and abs(b[2]) < 1e-6
    frac = np.linalg.solve(lat.T, x0.T).T % 1.0
    x0 = frac @ lat
    Lz = args.nc * c[2]

    # supercell: i,j in [-m,m] along a,b; k in [0,nc) along c
    m = int(np.ceil((args.R + 15) / min(np.linalg.norm(a), b[1]))) + 1
    shifts = np.array([i * a + j * b + k * c for i in range(-m, m + 1)
                       for j in range(-m, m + 1) for k in range(args.nc)])
    X = (x0[None, :, :] + shifts[:, None, :]).reshape(-1, 3)
    S = np.tile(sp0, len(shifts))
    # center the cylinder axis on the centroid of the central cell
    ctr = (0.5 * (a + b))[:2]
    X[:, :2] -= ctr
    r = np.hypot(X[:, 0], X[:, 1])
    crop = r < args.R + 6.0
    X, S, r = X[crop], S[crop], r[crop]

    # optional penny-shaped crack: disk of radius a, opening h, normal tilted from z
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
    ncrack = select(np.zeros(len(S), bool))[0].sum() - keep.sum()

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
            v = X[o] - X[s]
            v[2] -= Lz * np.round(v[2] / Lz)
            v /= np.linalg.norm(v)
            rad = np.array([X[o, 0], X[o, 1], 0.0]) / max(r[o], 1e-6)
            u = v + rad
            u /= np.linalg.norm(u)
            newH.append(X[o] + 0.97 * u)
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
    return Xk, Sk, Lz, q, len(newH), ncrack


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
            f" az={args.crack_az} zc={args.crack_zc}")


def write(args, X, S, Lz):
    half = args.R + args.vac
    order = np.argsort([TYPES.index(e) for e in S], kind="stable")
    X, S = X[order], S[order]
    # box [0, 2*half] in x,y with the pillar axis at its center (viewers assume a 0-origin cell)
    X = X + np.array([half, half, 0.0])
    with open(f"{args.out}/pillar.data", "w") as f:
        f.write(f"tobermorite-11A pillar R={args.R} nc={args.nc} (axis = c = z){crack_label(args)}\n\n")
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
    p.add_argument("--xyz", default="/home1/mm18756/tobermolite/info/tobe.xyz")
    p.add_argument("--R", type=float, default=20.2)
    p.add_argument("--nc", type=int, default=4)
    p.add_argument("--vac", type=float, default=40.0)
    p.add_argument("--out", default=".")
    p.add_argument("--dry", action="store_true", help="report only, no files")
    p.add_argument("--crack-a", type=float, default=0.0, help="penny crack radius (A); 0 = none")
    p.add_argument("--crack-h", type=float, default=4.0, help="crack opening thickness (A)")
    p.add_argument("--crack-tilt", type=float, default=45.0, help="crack normal vs z (deg)")
    p.add_argument("--crack-az", type=float, default=0.0, help="azimuth of the tilt in xy (deg)")
    p.add_argument("--crack-zc", type=float, default=0.5, help="crack center, fraction of Lz")
    args = p.parse_args()
    X, S, Lz, q, nadd, ncrack = build(args)
    els, cnt = np.unique(S, return_counts=True)
    qf = sum(QFORMAL[e] * n for e, n in zip(els, cnt))
    rr = np.hypot(X[:, 0], X[:, 1])
    print(f"R={args.R} Lz={Lz:.4f} N={len(S)} " + " ".join(f"{e}{n}" for e, n in zip(els, cnt)))
    print(f"formal charge before balance={q}, H added={nadd}, final={qf}, r_max={rr.max():.2f}")
    if args.crack_a > 0:
        print(f"penny{crack_label(args)}: {ncrack} atoms removed")
    if not args.dry:
        print(f"min interatomic distance={min_dist(X, Lz):.3f} A")
        write(args, X, S, Lz)


if __name__ == "__main__":
    main()
