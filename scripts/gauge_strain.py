#!/usr/bin/env python3
"""Gauge strain and mid-height bulge of a compressed finite pillar (numpy only).

Usage: python gauge_strain.py traj_compress.lammpstrj OUT.dat [--strain-per-frame 0.001] [--half 2]

Gauge strain = -(relative change of the z distance between the centroids of the atoms that start
(frame 0) in two thin slabs |z - (z_mid +- L/4)| < half). It measures the pillar's own strain and
leaves out the rigid grips. Bulge = relative change of the mean radius of the atoms that start in
the mid-height slab |z - z_mid| < half, about their own xy centroid.
Writes OUT.dat (platen_strain gauge_strain bulge) and prints every 10th frame.
"""
import argparse
import numpy as np


def frames(path):
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                return
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            f.readline()
            f.readline()
            n = int(f.readline())
            for _ in range(4):
                f.readline()
            cols = f.readline().split()[2:]
            c = [cols.index(a) for a in "xyz"]
            a = np.loadtxt([f.readline() for _ in range(n)], usecols=c)
            yield a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("traj")
    ap.add_argument("out")
    ap.add_argument("--strain-per-frame", type=float, default=1.0e-3)
    ap.add_argument("--half", type=float, default=2.0, help="slab half-thickness (A)")
    a = ap.parse_args()
    rows = []
    for k, x in enumerate(frames(a.traj)):
        if k == 0:
            z = x[:, 2]
            zmid, L = 0.5 * (z.min() + z.max()), z.max() - z.min()
            lo = np.abs(z - (zmid - L / 4)) < a.half
            hi = np.abs(z - (zmid + L / 4)) < a.half
            mid = np.abs(z - zmid) < a.half
        g = x[hi, 2].mean() - x[lo, 2].mean()
        m = x[mid, :2]
        r = np.hypot(*(m - m.mean(0)).T).mean()
        if k == 0:
            g0, r0 = g, r
        rows.append((k * a.strain_per_frame, 1 - g / g0, r / r0 - 1))
        if k % 10 == 0:
            print(f"{100 * rows[-1][0]:5.1f} %  gauge {100 * rows[-1][1]:6.2f} %  "
                  f"bulge {100 * rows[-1][2]:6.2f} %", flush=True)
    np.savetxt(a.out, rows, header="platen_strain gauge_strain bulge", fmt="%.6g")


if __name__ == "__main__":
    main()
