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

## Next steps

- Look at the crack1 trajectory near 8% strain and at the peak (14%) to see whether
  the kink is crack-face sliding or rim growth, and why its post-peak stress stays higher.
- Then decide on the deferred tilt/radius sweeps.
