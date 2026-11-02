# Penny-crack rerun of the ReaxFF tobermorite-11A pillar workflow

Status as of 2026-10-30. Baseline for the changes below is the original
`tobermolite-compression.tar.gz`.

## Changes made

### `scripts/build_pillar.py`: optional embedded penny-shaped crack

New options (no crack options = same intact pillar as before):

| Flag | Default | Meaning |
|---|---|---|
| `--crack-a` | `0` (no crack) | Crack radius (Å) |
| `--crack-h` | `4.0` | Crack opening thickness (Å) |
| `--crack-tilt` | `45` | Angle between crack normal and z (loading) axis (deg) |
| `--crack-az` | `0` | Azimuth of the tilt in the xy plane (deg) |
| `--crack-zc` | `0.5` | Crack center height, fraction of Lz |

How it works:

- Atoms within `h/2` of the tilted crack plane and within `a` of the pillar
  axis are flagged. The calculation handles periodicity in z.
- Atom selection is now a function, `select(excl)`. It removes whole units
  (SiO4 tetrahedron, Ca, or H2O/OH) whose Si, Ca or water-O center is
  flagged, so no broken tetrahedra or orphan H are left behind.
- Charge balancing runs after crack removal, so the structure stays neutral.
- Safety checks:
  - `a < R - 5`: the crack stays at least 5 Å inside the pillar surface.
  - `a*sin(tilt) + h/2 < Lz/2`: the crack does not overlap its periodic z image.
- Prints the number of atoms removed by the crack and writes the crack
  parameters into the `pillar.data` header line.
- Default `--xyz` path changed from `/home1/anakano/...` to `/home1/mm18756/...`.

### `slurm/run_pillar.slurm`: supports multiple runs

- Account: `anakano_81` → `anakano_429`.
- Arguments after the stage name are passed to LAMMPS (`"${@:2}"`), e.g.
  `sbatch run_pillar.slurm equil -var nnpt 400`.
- Results are copied to `results/<run directory name>` instead of always
  `results/run1`, so a crack run does not overwrite the intact run.
- Paths in the usage comments and the `srun -in` path changed from `anakano` to `mm18756`.

### `inputs/in.equil`, `inputs/in.compress`

- Only the `ffield` path changed (`anakano` → `mm18756`). Simulation
  settings are unchanged.

## Current state

- `/scratch1/mm18756/tobermolite/run1` is the intact pillar (R = 20.0, nc = 4,
  N = 10140, Lz = 93.68 Å). Its results are backed up in `results/run1`.
  run1 compress took 1 h 50 m on 32 tasks.
- Account check done (2026-10-30): `sacctmgr` lists `anakano_429` (QOS normal).
- Crack parameters chosen for a single first run, `crack1`: a = 10 Å, h = 4 Å,
  tilt = 45° (shear-maximizing), az = 0, zc = 0.5. Sweeps over tilt or radius
  are deferred until crack1 is compared with run1.
- Crack location in the built `pillar.data` (pillar axis at x = y = 60 Å, Lz = 93.68 Å):
  - Center (60, 60, 46.84) Å, which is on the pillar axis at Lz/2.
  - Normal (0.707, 0, 0.707), tilted 45° from z toward +x. The disk plane contains y
    and follows z − 46.84 ≈ −(x − 60).
  - Extent: y 50–70 Å. In (x, z) it runs from about (52.9, 53.9) to (67.1, 39.8) Å.
  - The edge is 10 Å inside the side surface. The crack does not reach the z
    boundaries. Look along y to see it edge-on.
  - After NPT, Lz = 93.05 Å (z 0.31–93.37), and the crack stays at mid-height
    (z ≈ 46.8 Å).
- Viewing the crack in OVITO 3.x:
  - **Which file.** Use `/scratch1/mm18756/tobermolite/crack1/traj_compress.lammpstrj`.
    It has the cell and an `element` column. Frame k = k × 0.1% strain, so frame 0
    is the post-NPT state, frame 80 is the kink at 8% and frame 142 is the peak.
    For the as-built crack use `pillar0.xyz`. Avoid `traj_*.xyz` because it has no
    cell. When importing, check that `element` is mapped to Particle Type.
  - **Edge-on view (best for watching opening/sliding).** Use the *Front* viewport,
    which looks along +y with x across and z up. Add *Slice*: normal (0, 1, 0),
    distance 60, slab width 6 Å. The crack shows as an empty band through the
    pillar center, running from upper-left (52.9, 53.9) to lower-right (67.1, 39.8).
    This slice stays valid for every frame, because compression only changes z.
  - **Face-on view (shows the disk outline).** Add *Slice*: normal (0.7071, 0, 0.7071),
    distance 75.54 (= n · center), slab width 4 Å. Then use *Adjust View* to set
    the view direction to (−0.7071, 0, −0.7071). The crack appears as a hole about
    20 Å across in the cross-section. This slice is only correct near frame 0,
    because the crack center moves down as Lz shrinks.
  - **Highlight the crack-face atoms in any frame.** Use *Expression selection*,
    then *Assign color*. Let dx = X − 60, dy = Y − 60 and dz = Z − (zlo + Lz/2).
    For a normal tilted 45° in xz, the distance from the crack plane is
    dn = (dx + dz)/√2, and the distance from the disk axis within the plane is
    ρ² = (dx − dz)²/2 + dy². Expressions:
    - Crack faces (|dn| < 4, ρ < 12; 230 atoms at frame 0):
      `abs((Position.X-60) + (Position.Z-CellOrigin.Z-0.5*CellSize.Z)) < 5.66 && 0.5*((Position.X-60) - (Position.Z-CellOrigin.Z-0.5*CellSize.Z))^2 + (Position.Y-60)^2 < 144`
    - Original crack region (|dn| < 2, ρ < 10). This should be nearly empty: 27 atoms
      at frame 0, compared with about 115 for the same region in run1:
      `abs((Position.X-60) + (Position.Z-CellOrigin.Z-0.5*CellSize.Z)) < 2.83 && 0.5*((Position.X-60) - (Position.Z-CellOrigin.Z-0.5*CellSize.Z))^2 + (Position.Y-60)^2 < 100`
    - The z center follows the shrinking box. The 45° tilt is fixed, though, so at
      large strain the crack, which flattens as it is compressed, drifts out of the
      selection. Widen the 5.66 if the faces drop out. Turning on *Make semi-transparent* for the
    unselected atoms (*Compute property* → Transparency = 0.8 for
    `Selection==0`) makes the faces visible inside the pillar.
  - The pillar's xy centroid stays within 0.2 Å of (60, 60): (60.17, 59.96) at
    frame 0 and (60.12, 59.94) at frame 300, so x = y = 60 works for all 301 frames.
- Dry run (`--dry`, no files written): N = 10047, 106 atoms removed by the
  crack, 237 H added, final charge 0.
- Built 2026-10-30 in `/scratch1/mm18756/tobermolite/crack1` (`pillar.data`,
  `pillar0.xyz`). The output matches the dry run, the minimum interatomic distance is
  0.970 Å (the added O–H) and the crack label is in the `pillar.data` header.
- Submitted 2026-10-30 from `crack1/`: equil job 13616286, and compress job
  13616287 (`--dependency=afterok:13616286`). Results will be copied to `results/crack1/`.

- Both jobs COMPLETED with exit code 0 (equil 1 h 43 m, compress 1 h 49 m). Results
  were copied to `results/crack1/`.

## Results (2026-11-02)

### 1. Crack still open after equilibration: yes

Note: `relaxed.data` is written right after the CG minimization (timestep 778),
*before* the Berendsen/NPT/NVT stages. The post-NPT state is `thermal.restart`, and it
is the same as frame 0 of `traj_compress.lammpstrj`, so that frame was checked too.

I counted atoms in the original crack region (centered on the pillar xy centroid at
`zlo + 0.5 Lz`, normal tilted 45° in xz) and in the same region of run1:

| State | run1 slab (\|dn\|<2, ρ<10) | crack1 slab | run1 core (\|dn\|<1, ρ<8) | crack1 core |
|---|---|---|---|---|
| as built (`pillar.data`) | 109 | 17 | 33 | 0 |
| minimized (`relaxed.data`) | 114 | 20 | 37 | 0 |
| after NPT+NVT (lammpstrj frame 0) | 115 | 30 | 37 | 0 |

- The crack core stays completely empty. The 4 Å slab holds about 26% of the
  intact atom count (30 vs 115), up from 16% as built. That is some closure at the
  crack faces/rim, mostly O (20) and Si (5), with no Ca.
- Lz after NPT: crack1 93.05 Å, run1 93.17 Å (93.68 Å as built).
- Script: `/scratch1/mm18756/tmp/crack_open.py` (not backed up).

### 2. Stress–strain comparison

Plot: `results/compare_crack1.png`. Curves: `results/stress_strain_{intact,crack}.dat`.

| Run | Peak stress (GPa) | Strain at peak | Initial modulus, ε<2% (GPa) | Mean stress, ε=20–30% (GPa) |
|---|---|---|---|---|
| intact (run1) | 11.44 | 0.120 | 124.2 | 7.53 |
| crack (crack1) | 10.09 | 0.142 | 113.5 | 8.98 |
| crack / intact | −11.8% | +0.022 | −8.6% | +19% |

- The two curves are nearly the same up to about 5% strain. Above that the crack1
  curve softens, with a kink near 8% strain (about 7.5 GPa), and peaks later and lower.
- After the peak, the intact pillar drops to about 7.5 GPa, but crack1 levels off
  higher, at about 9 GPa.

### 3. Damage maps (2026-11-02)

Script: `scripts/damage_map.py` (numpy only). It computes per-atom von Mises shear
strain η (Shimizu, Ogata & Li 2007; the same as OVITO Atomic Strain, cutoff 4 Å) and D²min
relative to frame 0 of `traj_compress.lammpstrj`. It also flags Si that lost a frame-0
Si–O bond (d > 2.3 Å) and Ca whose Ca–O coordination (< 2.9 Å) changed by 2 or more.

```
python scripts/damage_map.py results/damage \
  intact=/scratch1/mm18756/tobermolite/run1/traj_compress.lammpstrj \
  crack=/scratch1/mm18756/tobermolite/crack1/traj_compress.lammpstrj \
  --strains 5,8,12,14,20,30 --crack crack
```

It takes about 17 s on the login node. Outputs go to `results/damage/` (not tracked by git):
`shear_slices.png` (edge-on y-slab), `shear_xsection.png` (mid-height xy-slab),
`damage_vs_strain.png`, `damage_summary.dat`, and `damage_<run>_<strain>.xyz` (extended
XYZ with `eta d2min broken` columns, for OVITO).

Checks: frame 0 against itself gives η = 0 and F_zz = 1.000. At 5 % strain the core median
F_zz is 0.951.

| Run | ε | core p99 η | fraction η > 0.2 | core median F_zz | broken Si–O | Ca Δcn ≥ 2 |
|---|---|---|---|---|---|---|
| intact | 8 % | 0.13 | 4.3 % | 0.920 | 0 | 22 |
| crack | 8 % | 0.27 | 6.9 % | 0.942 | 0 | 23 |
| intact | 14 % | 0.45 | 29 % | 0.871 | 0 | 66 |
| crack | 14 % | 0.57 | 24 % | 0.899 | 0 | 60 |
| intact | 30 % | 1.52 | 64 % | 0.642 | 1 | 147 |
| crack | 30 % | 1.99 | 62 % | 0.860 | 4 | 117 |

Findings:

- **Intact (run1).** Shear stays spread out until the peak (12 %). At 14 % a horizontal
  band forms across the full width at z ≈ 55–70 Å. By 20–30 % there are two wide bands.
  The layers are ⊥ c = z, so this is layer-parallel (interlayer) collapse, not a 45° shear plane.
- **Crack (crack1).** At 8 % (the kink) shear concentrates at the crack rim. Core p99 η is
  twice the intact value. By 12 % the damage zone runs **horizontally** across the pillar at
  the crack height (z ≈ 40–50 Å). It does not follow the 45° crack plane. By 30 % the band
  has widened and the pillar bulges (x up to about ±28 Å).
- **Higher post-peak stress in crack1.** At 30 % the crack1 core is compressed much less
  than run1 (core F_zz 0.86 vs 0.64). Shortening is concentrated in the single crack band
  while the rest of the pillar stays nearly elastic. In run1, two bands collapse.
- **Bonds.** At most 4 Si–O bonds break. The silicate chains survive, so the damage is
  Ca–O and interlayer rearrangement.

## Next steps

- Done (§3): the kink is rim damage that turns into a horizontal band at the crack height.
- Then decide on the deferred tilt/radius sweeps.
