# Tobermorite-11 Å nanopillar compression: ReaxFF vs ClayFF (LAMMPS, USC CARC Discovery)

Workflow for carving a cylindrical tobermorite-11 Å pillar from a bulk cell, equilibrating it
at 10 K, compressing it uniaxially along the pillar axis (crystal **c** = z) at 1 %/ps, and
comparing stress–strain curves from two force fields. Both reference runs were completed on
Discovery (October 2026). The workflow is meant to be rerun with **other pillar sizes and other
starting configurations**.

This file is for a human reader and for Claude Code. Section 8 has the Claude-specific rules.

---

## 1. What's in the tarball

```
tobermolite/
├── README.md                    this file
├── info/
│   ├── tobe.xyz                 bulk tobermorite-11 Å, extended XYZ (frame 1 used; 1408 atoms/cell)
│   └── ffield-CaSiOH-lammps.reax.txt   ReaxFF Ca/Si/O/H (Manzano et al., JACS 2012)
├── scripts/
│   ├── build_pillar.py          carve + charge-balance pillar -> pillar.data (atom_style charge)
│   ├── plot_stress_strain.py    single log -> stress_strain.dat + PNG (also provides the log parser)
│   └── plot_compare.py          overlay several logs (e.g. ReaxFF vs ClayFF), prints peak + modulus
├── inputs/                      ReaxFF LAMMPS inputs
│   ├── in.equil                 CG minimize -> Berendsen 5 ps -> NPT(z) 10 ps -> NVT 5 ps
│   └── in.compress              fix deform z, erate -1e-5 /fs, 30 ps -> 30 % strain
├── slurm/run_pillar.slurm       ReaxFF job script (stage = equil | compress)
├── lammps/                      LAMMPS 30Sep2026 stable, packages REAXFF OPENMP (CPU, MPI)
├── clayff/
│   ├── scripts/build_clayff_data.py  ReaxFF pillar.data -> ClayFF pillar_clayff.data (atom_style full)
│   ├── inputs/ff.clayff         ClayFF parameters (Cygan, Liang, Kalinichev, JPCB 2004)
│   ├── inputs/in.equil, in.compress  same protocol, 0.5 fs timestep, PPPM
│   ├── slurm/run_pillar.slurm   ClayFF job script
│   └── lammps/                  LAMMPS 30Sep2026 stable, packages KSPACE MOLECULE OPENMP
└── reference/                   outputs of the reference runs, for checking a rerun
    ├── reaxff/  pillar.data log.equil log.compress species.out stress_strain_ReaxFF.dat
    ├── clayff/  pillar_clayff.data log.equil log.compress stress_strain_ClayFF.dat
    └── stress_strain_reaxff_vs_clayff.png
```

Not included: the literature PDFs (copyright), trajectories (~250–400 MB per run), restart
files and `relaxed.data`/`compressed.data`. All of these can be regenerated.

## 2. Reference result (R = 20 Å, 4 unit cells along c, 10 140 atoms, 10 K, 1 %/ps)

| Force field | Peak stress | Strain at peak | Initial modulus (fit, strain < 2 %) | Plateau 20–30 % |
|---|---|---|---|---|
| ReaxFF | 11.4 GPa | 11.9 % | 124 GPa | ~7–7.5 GPa |
| ClayFF | 8.8 GPa  | 11.4 % | 108 GPa | ~6–7 GPa |

The two curves almost overlap up to ~6 % strain. ClayFF then softens and peaks lower. Both
fail by stepwise load drops at about the same strain. In ReaxFF, about 7 of 870 water
molecules dissociate during compression (`species.out`).

Cost on 1 node, 32 MPI ranks (`main` partition):

| | equil | compress |
|---|---|---|
| ReaxFF (0.25 fs, QEq every step) | 2.0 h | 1.8 h (0.40 ns/day) |
| ClayFF (0.5 fs, PPPM 1e-5) | 10 min | 15 min |

Cost scales roughly linearly with atom count. The ReaxFF compress stage of a pillar with 4×
the atoms will not fit in the 8 h limit on one node. Request more nodes, or use `oneweek`.

## 3. Setup on Discovery (once)

The scripts contain absolute paths of the original user (`/home1/anakano/...`,
`/scratch1/anakano/...`) and the account `anakano_81`. Extract to `~/tobermolite` and rewrite
the paths:

```bash
cd ~ && tar xzf tobermolite-compression.tar.gz          # creates ~/tobermolite
cd ~/tobermolite
grep -rlI --exclude-dir=lammps --exclude-dir=reference --exclude=README.md -e /home1/anakano -e /scratch1/anakano . \
  | xargs sed -i -e "s#/home1/anakano#/home1/$USER#g" -e "s#/scratch1/anakano#/scratch1/$USER#g"
sacctmgr -n show assoc user=$USER format=account%30      # your SLURM account(s)
sed -i "s/--account=anakano_81/--account=<YOUR_ACCOUNT>/" slurm/run_pillar.slurm clayff/slurm/run_pillar.slurm
```

The bundled `lmp` binaries are dynamically linked against the cluster's
`gcc/13.3.0` and `openmpi/5.0.5` modules, and the SLURM scripts load those modules. Check the
binaries with:

```bash
module purge; module load gcc/13.3.0 openmpi/5.0.5
~/tobermolite/lammps/bin/lmp -h | grep -A2 "Installed packages"          # OPENMP REAXFF
~/tobermolite/clayff/lammps/bin/lmp -h | grep -A2 "Installed packages"   # KSPACE MOLECULE OPENMP
```

The Python scripts need only numpy and matplotlib. The default `python/3.11.9` module from the
`usc` collection has both.

## 4. Running one size (the reference case)

Run in `/scratch1/$USER`. The trajectories are large. `/scratch1` is **not backed up**,
so each SLURM script copies logs and the small data files back to `~/tobermolite/.../results/run1`.

### ReaxFF

```bash
mkdir -p /scratch1/$USER/tobermolite/run1 && cd $_
python ~/tobermolite/scripts/build_pillar.py --R 20.0 --nc 4 --vac 40 --out .
#   -> pillar.data, pillar0.xyz; prints composition, charge balance, min distance
jid=$(sbatch --parsable --job-name=tobe-equil ~/tobermolite/slurm/run_pillar.slurm equil)
sbatch --dependency=afterok:$jid --job-name=tobe-compress ~/tobermolite/slurm/run_pillar.slurm compress
```

### ClayFF (uses the *same atoms* as the ReaxFF pillar)

```bash
mkdir -p /scratch1/$USER/tobermolite/clayff/run1 && cd $_
python ~/tobermolite/clayff/scripts/build_clayff_data.py --src ../../run1/pillar.data --out .
jid=$(sbatch --parsable --job-name=clayff-equil ~/tobermolite/clayff/slurm/run_pillar.slurm equil)
sbatch --dependency=afterok:$jid --job-name=clayff-compress ~/tobermolite/clayff/slurm/run_pillar.slurm compress
```

### Plot

```bash
cd ~/tobermolite && mkdir -p results/compare
python scripts/plot_compare.py results/compare/stress_strain.png \
  ReaxFF=results/run1/log.compress ClayFF=clayff/results/run1/log.compress
```

This prints peak stress, strain at peak and initial modulus for each run, and writes
`stress_strain_<label>.dat` (time_ps, strain, stress_GPa, T, PE) next to the PNG. The plot
title is hard-coded in `plot_compare.py`. Edit it for other conditions.

## 5. Changing the system size

`build_pillar.py` options:

| Option | Meaning | Reference |
|---|---|---|
| `--R` | pillar radius (Å); whole SiO₄ tetrahedra and H₂O/OH units are kept | 20.0 |
| `--nc` | number of unit cells along c (periodic pillar length = nc × 23.42 Å) | 4 |
| `--vac` | vacuum beyond R on each side in x, y (box = 2(R+vac)) | 40 |
| `--xyz` | source bulk structure | info/tobe.xyz |
| `--dry` | report composition and charge only, write nothing | |

Use `--dry` to scan R values first. The script charge-balances by protonating the outermost
non-bridging O or removing the outermost Ca. It exits with *"odd positive charge"* for some
R values. Change R by ~0.1 Å and retry.

**Things that must be changed together with the size:**

1. **`R0` in all four `in.*` files** (`variable R0 index 20.35`). This is the radius of
   a bulk-density cylinder holding the same number of atoms, and it is used for the stress
   cross-section `σ = pzz·lx·ly / (π R0²)`. Recompute it for every new pillar:
   ```bash
   python -c "import math; N=10140; nc=4; n=1408/(28.01*25.80*23.42); print(math.sqrt(N/(n*math.pi*nc*23.42)))"
   # -> 20.35 ; replace N (atom count printed by build_pillar.py) and nc
   ```
   You can pass it on the command line (`-var R0 23.1`) instead of editing. The SLURM
   scripts don't do this, so either edit the inputs or add `-var R0 ...` to the `srun` line.
   For a different bulk structure, replace 1408 and the cell volume.
2. **Box size for ClayFF.** ClayFF runs fully periodic (`p p p`, PPPM). The vacuum must stay
   large enough that periodic images of the pillar don't interact. Keep `--vac` ≥ 40 Å. The
   ReaxFF runs use shrink-wrapped `s s p` boundaries, so the vacuum only matters there for
   the stress normalization (which uses R0, not the box).
3. **MPI layout.** `clayff/inputs/in.*` hard-code `processors 2 2 8`, which requires
   exactly 32 ranks. Change it if you change `--ntasks`, and keep 2×2 in x, y so that the cuts
   go through the pillar axis. ReaxFF uses `balance … shift xy` and a `fix balance` instead.
4. **SLURM resources.** Scale `--time` (and nodes or ranks) with atom count, using the
   costs in §2.
5. **Run directory and copy-back destination.** Each SLURM script copies results to a fixed
   `results/run1`. For a size sweep, use one scratch dir per case (e.g. `R25_nc6`) and change
   `dest=` in the SLURM script (for example to `.../results/$(basename $PWD)`). Otherwise every
   case overwrites the same backup.

Strain rate and total strain are controlled by `erate` and `ncomp` in `in.compress`
(`ncomp × timestep × |erate|` = final strain; ReaxFF timestep 0.25 fs, ClayFF 0.5 fs).
The temperature is set by `T`. If you raise T, read the note about `fix reaxff/species`
in §7.

## 6. Using a different atomic configuration

- `build_pillar.py` reads **frame 1** of an extended-XYZ file with a `Lattice="..."` header,
  elements **Ca, Si, O, H** only. It requires **c ∥ z** and a, b in the xy plane (the script
  asserts this). For other orientations, rotate the cell first.
- Bonding is detected by distance: Si–O < 2.0 Å, O–H < 1.2 Å. Check these cutoffs for strongly
  distorted or amorphous structures.
- `build_clayff_data.py` assigns ClayFF types from coordination. It only knows these
  O environments: (n_Si, n_H) = (2,0) ob, (1,0) nbo, (1,1) oh, (0,2) water. A structure with
  Ca–OH hydroxyls (0 Si, 1 H), three-coordinated O, etc. fails with a `KeyError`. Add the
  missing types to `TYPES`, the lookup table, and `ff.clayff`. The nbo charge is solved for
  exact neutrality and printed. Check it is still physically reasonable (reference:
  −1.366 e). Ca next to water O (< 2.9 Å) is typed `ca2` (interlayer), and all other Ca `cao`.
- To compare force fields on identical atoms, always build the ClayFF data from the
  ReaxFF `pillar.data`, as in §4.
- **Finite (non-periodic) pillars between platens**: `clamped/build_pillar.py --R .. --L ..`
  (c ∥ x orientation), then `inputs/in.equil_finite` and `inputs/in.compress_finite`
  (`-var ends clamped|wall -var R0 ..`, `boundary s s s`), analysed with `scripts/stress_map.py`.
  Motivation, settings and results are in `SPLIT_PLAN.md`.

## 7. Known warnings and caveats

- ReaxFF prints at startup *"Van der Waals parameters for element SI indicate inner
  wall+shielding…"* and *"Changed valency_val to valency_boc for X"*. These come from the
  force-field file. The reference runs were stable (no lost atoms, energy smooth).
- `fix reaxff/species` in the ReaxFF `in.compress` forces neighbor rebuilds to
  `every 400 check no`, so the log says *"Dangerous builds not checked"*. With a 2 Å skin
  at 10 K this was fine. At higher temperature or strain rate, drop the species fix or reduce
  its interval (`fix reaxff/species 100 1 100 …`).
- The stress is the full virial `pzz` (kinetic part included, negligible at 10 K) divided by
  a **fixed** cross-section π R0². It is an engineering stress, and the modulus is the slope
  below 2 % strain.
- Disk use: each run leaves ~250 MB (ReaxFF) to ~400 MB (ClayFF) of trajectories on
  `/scratch1`. Copy anything worth keeping to `/project2/<group>/` (use `cp -r`, never
  `cp -a` or `mv`, into `/project2`).

## 8. Notes for Claude Code

- These are CARC Discovery rules. The organization `CLAUDE.md` at `/etc/claude-code/CLAUDE.md`
  is authoritative.
- Never run `lmp`/`srun` on a login node. Every LAMMPS run goes through `sbatch` with the
  scripts above. Show the SLURM script (and any edits) to the user and get approval for each
  `sbatch`. Use `--partition=debug` (≤ 1 h) for a quick "does it start" test of new inputs,
  e.g. with `ncomp`/`nnpt` reduced via `-var`.
- `build_pillar.py`, `build_clayff_data.py` and the plot scripts are light (seconds) and can
  run on the login node after the user approves the prompt.
- Before submitting a new size, check:
  `R0` recomputed (§5.1); `processors` matches `--ntasks` (ClayFF); `--time` scaled;
  `dest=` points at a per-case directory; `--account` is the user's own (`sacctmgr show assoc user=$USER`).
- Monitor with `squeue -u $USER` and `tail log.<stage>` at long intervals, not tight loops.
  The API rate limit is shared across the user's sessions.
- Check a finished run with `grep "Total wall" log.*` and the `ERROR`/`WARNING` lines, then
  `plot_compare.py`. Compare against `reference/` when rerunning the reference case. Peak
  stress should match to within ~0.1 GPa.

## 9. References

- Force field (ReaxFF Ca/Si/O/H): H. Manzano *et al.*, *J. Am. Chem. Soc.* **134**, 2208 (2012).
- ClayFF: R. T. Cygan, J.-J. Liang, A. G. Kalinichev, *J. Phys. Chem. B* **108**, 1255 (2004).
- Tobermorite background: Jiang *et al.*, *Phys. Chem. Chem. Phys.* (2018). This is the
  `Jiang-tobermolite-PCCP18.pdf` in the original project's `papers/`.
