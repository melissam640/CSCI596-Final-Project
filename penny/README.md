# rot_penny_x10 build files (penny_R43)

Everything needed to rebuild the rot_penny_x10 pillar (R = 43 Å, 9 unit cells along the axis,
Ca-sheet penny a = 21.5 Å), built 2026-11-12. See `CRACK_PLAN.md`, section "rot_penny_x10".

| File | What it is |
|---|---|
| `build_pillar.py` | the builder **as of commit e59bd93**, the version used for this build (`clamped/build_pillar.py` has changed since, for finite pillars) |
| `tobe_rot.xyz` | rotated bulk cell (c ∥ x, a ∥ z); the builder reads it from its own directory by default |
| `rotate_cell.py`, `tobe.xyz` | how `tobe_rot.xyz` was made: `python rotate_cell.py --xyz tobe.xyz --out tobe_rot.xyz` (copied from `clamped/` and `info/`) |
| `build.log` | output of the original build (731 s) |
| `pillar.data`, `pillar0.xyz` | the built pillar (LAMMPS data, atom_style charge; XYZ for viewing) |

Rebuild (writes `pillar.data` and `pillar0.xyz` into `--out`; ~12–14 min, ~1 GB peak memory, use a compute node):

    python build_pillar.py --R 43.0 --nc 9 --vac 40 --crack-a 21.5 --crack-h 4 --out <dir>

Expected (`build.log`): N = 125917 (Ca 10811, Si 16848, O 62424, H 35834), 2354 H added,
final charge 0, 133 Ca removed by the penny, min distance 0.970 Å, **R0 = 43.72 Å**.

Rebuild check (2026-11-24, job 14764470, debug partition, numpy 2.4.6): `build.log` identical, and all
123 563 crystal atoms identical in `pillar.data`. The 2354 charge-balancing H go on the same O columns
(same x, y and count per column, same z modulo half a cell), but 2221 of them sit on a different, symmetry-
equivalent O along z. The builder ranks candidate O by radius with an unstable `np.argsort` (line 149),
and O copies along the axis tie, so which copy gets the H depends on the machine/numpy sort. The rebuilt
pillar is therefore equivalent, not byte-identical. **`pillar.data` here is the exact file that was run.**

Next steps used `inputs/in.equil` and `inputs/in.compress` (unchanged since) via
`slurm/run_pillar.slurm`, on 2 × 64 ranks with `-var R0 43.72` (jobs 14059080, 14059081).
