#!/usr/bin/env python3
"""Overlay stress-strain curves from several compression logs (e.g. ReaxFF vs ClayFF).

Usage: python plot_compare.py OUT.png LABEL=log.compress [LABEL=log.compress ...]
Writes OUT.png and, per run, stress_strain_<label>.dat next to it; prints peak stress
and initial (strain < 2%) modulus for each run.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_stress_strain import compression_block

COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]   # fixed categorical order
INK, MUTED, GRID = "#1f1f1e", "#6b6b68", "#e4e4e0"


def main():
    png, runs = sys.argv[1], [a.split("=", 1) for a in sys.argv[2:]]
    outdir = os.path.dirname(os.path.abspath(png))
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for (label, log), col in zip(runs, COLORS):
        cols, d = compression_block(log)
        c = {k: i for i, k in enumerate(cols)}
        t, e, s = d[:, c["Time"]] / 1000.0, d[:, c["v_strain"]], d[:, c["v_scomp"]]
        np.savetxt(os.path.join(outdir, f"stress_strain_{label}.dat"),
                   np.c_[t, e, s, d[:, c["Temp"]], d[:, c["PotEng"]]],
                   header="time_ps strain scomp_GPa temp_K pe_kcalmol", fmt="%.6g")
        ipk = np.argmax(s)
        lin = e < 0.02
        E = np.polyfit(e[lin], s[lin], 1)[0] if lin.sum() > 5 else np.nan
        print(f"{label:8s} peak {s[ipk]:.3f} GPa at strain {e[ipk]:.4f}; "
              f"initial modulus {E:.1f} GPa; final strain {e[-1]:.4f}")
        ax.plot(100 * e, s, lw=1.6, color=col, label=f"{label} (E₀ ≈ {E:.0f} GPa)")
        ax.plot(100 * e[ipk], s[ipk], "o", ms=6, color=col, mec="white", mew=1.5)
        ax.annotate(f"{s[ipk]:.1f} GPa", (100 * e[ipk], s[ipk]), xytext=(6, 6),
                    textcoords="offset points", color=INK, fontsize=9)
    ax.set_xlabel("engineering strain (%)", color=MUTED)
    ax.set_ylabel("compressive stress (GPa)", color=MUTED)
    ax.set_title("Tobermorite 11 Å pillar, 10 K, 1%/ps", color=INK, loc="left")
    ax.grid(color=GRID, lw=0.8)
    ax.tick_params(colors=MUTED)
    for sp in ax.spines.values():
        sp.set_color(GRID)
    ax.legend(frameon=False, labelcolor=INK)
    fig.tight_layout()
    fig.savefig(png, dpi=150)


if __name__ == "__main__":
    main()
