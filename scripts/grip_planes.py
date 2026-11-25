#!/usr/bin/env python3
"""Equal-solid grip planes for a finite pillar (numpy only).

Usage: python grip_planes.py FRAME.lammpstrj [--grip 8]

Reads the first frame of a dump (columns with x y z) and returns the z planes at which each end's
grip holds the atoms of a solid slab `grip` thick, i.e. N_target = (N / L) * grip, with L the
atom z extent. Measuring the grip from the outermost atom instead (the in.compress_finite default)
under-fills it when an end face is rough or domed. Prints the -var gtop / gbot arguments for
in.compress_finite and the counts and depths of both definitions.
"""
import argparse
import numpy as np


def first_frame_z(path):
    with open(path) as f:
        for _ in range(3):
            f.readline()
        n = int(f.readline())
        for _ in range(4):
            f.readline()
        cols = f.readline().split()[2:]
        iz = cols.index("z")
        return np.array([float(f.readline().split()[iz]) for _ in range(n)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("--grip", type=float, default=8.0, help="equivalent solid thickness (A)")
    a = ap.parse_args()
    z = np.sort(first_frame_z(a.dump))
    lo, hi = z[0], z[-1]
    target = int(round(len(z) / (hi - lo) * a.grip))
    # plane halfway between the target-th atom and the next one from each end
    zbot = 0.5 * (z[target - 1] + z[target])
    ztop = 0.5 * (z[-target] + z[-target - 1])
    print(f"N = {len(z)}, z {lo:.3f} .. {hi:.3f}, target per grip = {target} atoms")
    print(f"outermost-atom grips ({a.grip} A): bottom {(z < lo + a.grip).sum()}, top {(z > hi - a.grip).sum()}")
    print(f"equal-solid grips: bottom z < {zbot:.3f} ({zbot - lo:.2f} A deep), "
          f"top z > {ztop:.3f} ({hi - ztop:.2f} A deep)")
    print(f"-var gbot {zbot:.4f} -var gtop {ztop:.4f}")


if __name__ == "__main__":
    main()
