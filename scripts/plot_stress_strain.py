#!/usr/bin/env python3
"""Parse the compression stage of log.lammps and plot engineering stress vs strain.

Usage: python plot_stress_strain.py [log.lammps] [out.png]
Writes stress_strain.dat (time_ps strain scomp_GPa temp pe) and a PNG.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def compression_block(path):
    rows, cols, inblock = [], None, False
    for line in open(path):
        w = line.split()
        if w[:2] == ["Step", "Time"] and "v_strain" in w:
            cols, inblock, rows = w, True, []
            continue
        if inblock:
            if w and w[0] == "Loop":
                inblock = False
                continue
            try:
                rows.append([float(v) for v in w])
            except ValueError:
                pass  # warnings interleaved with thermo output
    if cols is None:
        raise SystemExit("no compression thermo block (v_strain) in log")
    return cols, np.array([r for r in rows if len(r) == len(cols)])


def main():
    log = sys.argv[1] if len(sys.argv) > 1 else "log.lammps"
    png = sys.argv[2] if len(sys.argv) > 2 else "stress_strain.png"
    cols, d = compression_block(log)
    c = {k: i for i, k in enumerate(cols)}
    t, e, s = d[:, c["Time"]] / 1000.0, d[:, c["v_strain"]], d[:, c["v_scomp"]]
    np.savetxt("stress_strain.dat", np.c_[t, e, s, d[:, c["Temp"]], d[:, c["PotEng"]]],
               header="time_ps strain scomp_GPa temp_K pe_kcalmol", fmt="%.6g")
    ipk = np.argmax(s)
    lin = e < 0.02
    E = np.polyfit(e[lin], s[lin], 1)[0] if lin.sum() > 5 else np.nan
    print(f"peak stress {s[ipk]:.3f} GPa at strain {e[ipk]:.4f}; initial modulus {E:.1f} GPa")

    fig, ax = plt.subplots(figsize=(6, 4.2))
    ax.plot(100 * e, s, lw=1.2)
    ax.plot(100 * e[ipk], s[ipk], "o", ms=5)
    ax.set_xlabel("engineering strain (%)")
    ax.set_ylabel("compressive stress (GPa)")
    ax.set_title("Tobermorite 11 Å pillar, ReaxFF, 10 K, 1%/ps")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(png, dpi=150)


if __name__ == "__main__":
    main()
