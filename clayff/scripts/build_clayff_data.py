#!/usr/bin/env python3
"""Convert the ReaxFF pillar (atom_style charge, types Ca Si O H) into a ClayFF data file.

Usage: python build_clayff_data.py [--src ../run1/pillar.data] [--out DIR]

Uses the identical atoms/coordinates as the ReaxFF run, typed by connectivity:
  Ca touching a water O (< RCAW)       -> ca2  (aqueous Ca2+, interlayer)
  other Ca                             -> cao  (octahedral Ca, intralayer)
  Si                                   -> st
  O with 2 Si                          -> ob   (bridging)
  O with 1 Si, 0 H                     -> nbo  (non-bridging; charge solved for neutrality)
  O with 1 Si, 1 H / H on it           -> oh / ho  (silanol)
  O with 0 Si, 2 H / H on it           -> o* / h*  (water, flexible SPC)
Charges and LJ: Cygan, Liang, Kalinichev, J. Phys. Chem. B 108, 1255 (2004).
Writes pillar_clayff.data (atom_style full) with O-H bonds and H-O-H angles.
"""
import argparse
from collections import Counter
import numpy as np

# name, mass, charge (nbo filled in below)
TYPES = [("cao", 40.078, 1.36), ("ca2", 40.078, 2.0), ("st", 28.0855, 2.1),
         ("ob", 15.999, -1.05), ("nbo", 15.999, None), ("oh", 15.999, -0.95),
         ("ho", 1.008, 0.425), ("ow", 15.999, -0.82), ("hw", 1.008, 0.41)]
TID = {n: i + 1 for i, (n, _, _) in enumerate(TYPES)}
R_SIO, R_OH, RCAW = 2.0, 1.2, 2.9


def read_charge_data(path):
    L = open(path).read().split("\n")
    box = {}
    for l in L:
        w = l.split()
        if len(w) == 4 and w[2][1:] == "lo" and w[3][1:] == "hi":
            box[w[2][0]] = (float(w[0]), float(w[1]))
    i0 = next(k for k, l in enumerate(L) if l.startswith("Atoms")) + 2
    A = np.array([l.split()[:6] for l in L[i0:] if l.strip()], float)
    A = A[np.argsort(A[:, 0])]
    S = np.array(["Ca", "Si", "O", "H"])[A[:, 1].astype(int) - 1]
    return S, A[:, 3:6], box


def typify(S, X, Lz):
    def d(i, j):
        v = X[j][None] - X[i][:, None]
        v[..., 2] -= Lz * np.round(v[..., 2] / Lz)
        return np.linalg.norm(v, axis=2)

    iSi, iO, iH, iCa = (np.where(S == e)[0] for e in ("Si", "O", "H", "Ca"))
    DOH = d(iO, iH)
    assert np.all((DOH < R_OH).sum(0) == 1), "every H must have exactly one O"
    h_o = iO[np.argmin(DOH, axis=0)]
    nSi = (d(iO, iSi) < R_SIO).sum(1)
    nH = (DOH < R_OH).sum(1)

    T = np.empty(len(S), object)
    T[iSi] = "st"
    for o, ns, nh in zip(iO, nSi, nH):
        T[o] = {(2, 0): "ob", (1, 0): "nbo", (1, 1): "oh", (0, 2): "ow"}[(int(ns), int(nh))]
    T[iH] = np.where(T[h_o] == "ow", "hw", "ho")
    wat = iO[T[iO] == "ow"]
    T[iCa] = np.where((d(iCa, wat) < RCAW).any(1), "ca2", "cao")

    bonds = [(o, h) for h, o in zip(iH, h_o)]
    angles = []
    for o in wat:
        hs = iH[h_o == o]
        angles.append((hs[0], o, hs[1]))
    return T, bonds, angles


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--src", default="/scratch1/anakano/tobermolite/run1/pillar.data")
    p.add_argument("--out", default=".")
    args = p.parse_args()

    S, X, box = read_charge_data(args.src)
    Lz = box["z"][1] - box["z"][0]
    T, bonds, angles = typify(S, X, Lz)
    cnt = Counter(T)
    qfix = sum(q * cnt[n] for n, _, q in TYPES if q is not None)
    qnbo = -qfix / cnt["nbo"]
    Q = {n: (qnbo if q is None else q) for n, _, q in TYPES}
    qtot = sum(Q[t] for t in T)
    print(f"N={len(S)} " + " ".join(f"{n}={cnt[n]}" for n, _, _ in TYPES))
    print(f"q(nbo)={qnbo:.6f}  total charge={qtot:.2e}")

    # one molecule per water; everything else (silicate/Ca framework) is molecule 1
    mol = np.ones(len(S), int)
    for k, (_, o, _) in enumerate(angles):
        mol[o] = k + 2
    # z image flag on H so each O-H bond is intact when unwrapped (O keeps image 0)
    iz = np.zeros(len(S), int)
    for o, h in bonds:
        mol[h] = mol[o]
        iz[h] = int(np.round((X[o, 2] - X[h, 2]) / Lz))

    with open(f"{args.out}/pillar_clayff.data", "w") as f:
        f.write(f"tobermorite-11A pillar, ClayFF typing of {args.src}  q(nbo)={qnbo:.6f}\n\n")
        f.write(f"{len(S)} atoms\n{len(bonds)} bonds\n{len(angles)} angles\n\n")
        f.write(f"{len(TYPES)} atom types\n1 bond types\n1 angle types\n\n")
        for ax in "xyz":
            f.write(f"{box[ax][0]:.6f} {box[ax][1]:.6f} {ax}lo {ax}hi\n")
        f.write("\nMasses\n\n")
        for n, m, _ in TYPES:
            f.write(f"{TID[n]} {m}  # {n}\n")
        f.write("\nAtoms # full\n\n")
        for i in range(len(S)):
            f.write(f"{i+1} {mol[i]} {TID[T[i]]} {Q[T[i]]:.10f} "
                    f"{X[i,0]:.6f} {X[i,1]:.6f} {X[i,2]:.6f} 0 0 {iz[i]}\n")
        f.write("\nBonds\n\n")
        for k, (o, h) in enumerate(bonds, 1):
            f.write(f"{k} 1 {o+1} {h+1}\n")
        f.write("\nAngles\n\n")
        for k, (h1, o, h2) in enumerate(angles, 1):
            f.write(f"{k} 1 {h1+1} {o+1} {h2+1}\n")


if __name__ == "__main__":
    main()
