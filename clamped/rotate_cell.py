#!/usr/bin/env python3
"""Rotate the bulk tobermorite-11A cell by +90 deg about y so that crystal c lies along x.

Usage: python rotate_cell.py [--xyz ../info/tobe.xyz] [--out tobe_rot.xyz]

R_y(+90): (x, y, z) -> (z, y, -x).  Lattice rows a, b, c become
  c' = (|c|, 0, 0), b' = (0, b_y, -b_x), a' = (0, 0, -|a|).
The rotated cell is written with rows [c', b', -a'], so the third (periodic, pillar-axis)
vector is +z and of length |a|.  The layers (normal c) are now perpendicular to x.
Writes extended XYZ (species + positions) readable by clamped/build_pillar.py and OVITO.
"""
import argparse
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROT = np.array([[0.0, 0.0, 1.0],     # x' = z
                [0.0, 1.0, 0.0],     # y' = y
                [-1.0, 0.0, 0.0]])   # z' = -x


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


def min_dist_pbc(x, lat):
    """Minimum interatomic distance in a periodic cell (images -1..1)."""
    sh = np.array([i * lat[0] + j * lat[1] + k * lat[2]
                   for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
    best = np.inf
    for s in sh:
        d = np.linalg.norm(x[None, :, :] + s - x[:, None, :], axis=2)
        if not s.any():
            np.fill_diagonal(d, np.inf)
        best = min(best, d.min())
    return best


def sheets(v, period, gap=2.5):
    """Cluster 1D coordinates (mod period) into sheets; return sorted midplanes."""
    z = np.sort(v % period)
    groups = np.split(z, np.where(np.diff(z) > gap)[0] + 1)
    if len(groups) > 1 and z[0] + period - z[-1] <= gap:
        groups[0] = np.concatenate([groups[-1] - period, groups[0]])
        groups.pop()
    return np.sort(np.array([g.mean() for g in groups]) % period)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--xyz", default=os.path.join(HERE, "..", "info", "tobe.xyz"))
    p.add_argument("--out", default=os.path.join(HERE, "tobe_rot.xyz"))
    args = p.parse_args()

    sp, x0, lat = read_frame1(args.xyz)
    a, b, c = lat
    assert abs(c[0]) < 1e-6 and abs(c[1]) < 1e-6 and abs(a[2]) < 1e-6 and abs(b[2]) < 1e-6
    frac = np.linalg.solve(lat.T, x0.T).T % 1.0
    x0 = frac @ lat

    xr = x0 @ ROT.T
    ar, br, cr = lat @ ROT.T
    latr = np.array([cr, br, -ar])           # rows: c' (x), b' (yz), -a' (+z)
    # wrap into the new cell; the new cell is the same lattice, so fractions just permute
    fr = np.linalg.solve(latr.T, xr.T).T % 1.0
    xr = fr @ latr

    # ---- checks ----
    assert np.allclose(cr, [np.linalg.norm(c), 0, 0], atol=1e-6), cr
    assert abs(latr[2, 0]) < 1e-6 and abs(latr[2, 1]) < 1e-6 and latr[2, 2] > 0, latr[2]
    v0, v1 = abs(np.linalg.det(lat)), abs(np.linalg.det(latr))
    assert abs(v0 - v1) < 1e-6 * v0, (v0, v1)
    fexp = np.column_stack([frac[:, 2], frac[:, 1], (-frac[:, 0]) % 1.0])
    df = (fr - fexp + 0.5) % 1.0 - 0.5
    assert np.abs(df).max() < 1e-8, np.abs(df).max()
    d0, d1 = min_dist_pbc(x0, lat), min_dist_pbc(xr, latr)
    assert abs(d0 - d1) < 1e-6, (d0, d1)
    ca_x = sheets(xr[sp == "Ca", 0], latr[0, 0])
    ca_z = sheets(x0[sp == "Ca", 2], c[2])
    assert len(ca_x) == 2 and np.allclose(ca_x, ca_z, atol=1e-6), (ca_x, ca_z)

    print(f"cell volume {v1:.3f} A^3 (unchanged), min distance {d1:.4f} A (unchanged)")
    print("rotated lattice rows (A):")
    for name, row in zip(("c' ", "b' ", "-a'"), latr):
        print(f"  {name} {row[0]:9.4f} {row[1]:9.4f} {row[2]:9.4f}")
    print(f"Ca-sheet midplanes x (A, mod {latr[0, 0]:.2f}): {' '.join(f'{m:.2f}' for m in ca_x)}")

    with open(args.out, "w") as f:
        f.write(f"{len(sp)}\n")
        f.write('Lattice="' + " ".join(f"{v:.6f}" for v in latr.ravel()) + '" pbc="T T T" '
                "Properties=species:S:1:pos:R:3 "
                'comment="tobermorite-11A rotated +90 deg about y: c->x, a->-z; rows c\', b\', -a\'"\n')
        for e, x in zip(sp, xr):
            f.write(f"{e} {x[0]:.8f} {x[1]:.8f} {x[2]:.8f}\n")
    print(f"wrote {args.out} ({len(sp)} atoms)")


if __name__ == "__main__":
    main()
