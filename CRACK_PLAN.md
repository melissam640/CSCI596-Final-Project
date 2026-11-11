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

## crack2: horizontal penny crack through a Ca sheet

Motivation: in crack1 the failure became a horizontal (layer-parallel) band (§3). crack2
puts the crack flat (tilt 0°, normal ∥ z ∥ c), on the midplane of a Ca sheet.

- New `build_pillar.py` flag `--crack-snap-ca`. It clusters the bulk-cell Ca z values into
  sheets (midplanes at 5.91 and 17.63 Å in the 23.42 Å cell) and moves the crack centre to the
  sheet nearest `--crack-zc`. Mid-height (zc 0.5) is an interlayer, so snapping is required.
  The builder now also prints the removed atoms per element.
- Regression check: crack1 parameters without snapping still give N = 10047, 106 removed
  (Ca9 Si15 O49 H33) and 237 H added.
- Build (2026-11-02), `/scratch1/mm18756/tobermolite/crack2`:
  `build_pillar.py --R 20.0 --nc 4 --vac 40 --crack-a 10 --crack-h 4 --crack-tilt 0 --crack-zc 0.56 --crack-snap-ca`
  - The crack centre snapped to z = 52.75 Å (zc = 0.5631). The nearest Si planes are at
    −2.16 and +2.14 Å, just outside the ±2 Å opening, so no SiO4 is cut.
  - 27 atoms were removed, all Ca. N = 10167, 278 H added (outer surface, same rule as
    crack1), final charge 0, minimum interatomic distance 0.970 Å.
  - Check against run1, r < 10 Å and |z − 52.75| < 2: Ca 27 → 0, Si 0 → 0. Si total at
    r < 10 Å is 320 in both pillars.
- Submitted 2026-11-02 from `crack2/`: equil job 13678139, and compress job 13678140
  (`--dependency=afterok:13678139`). Results will be copied to `results/crack2/`.
- Viewing in OVITO: this is a horizontal disk, so look along y (*Front* viewport) with a
  y-slab at 60. To highlight it, use
  `abs(Position.Z - (CellOrigin.Z + 0.5631*CellSize.Z)) < 2 && (Position.X-60)^2 + (Position.Y-60)^2 < 100`.
  This assumes the crack scales with the box in z.
- Both jobs COMPLETED with exit code 0 (equil 1 h 48 m, compress 1 h 56 m). The only
  warnings are the usual ReaxFF force-field ones. Results were copied to `results/crack2/`.
  R0 stays at 20.35 Å, the same as for crack1, so the stresses are directly comparable.

### crack2 results (2026-11-03)

**1. Crack after equilibration.** Scripts: `/scratch1/mm18756/tmp/crack_open2.py`
(crack_open.py with tilt 0 and zc as an argument) and `crack2_gap.py` (not backed up).
Region: flat slab |z − zc| < 2 Å, r < 10 Å, centred at `zlo + 0.5631 Lz`.

| State | run1 slab | crack2 slab | z-gap at r < 8 Å, run1 | z-gap, crack2 |
|---|---|---|---|---|
| as built | 84 (Ca28 O56) | 56 (Ca1 O55) | 2.08 Å | 2.42 Å |
| minimized | 91 (Ca30 O61) | 53 (Ca2 O51) | | |
| after NPT+NVT (frame 0) | 93 (Ca30 Si2 O61) | 58 (Ca3 O55) | 1.95 Å | 1.84 Å |

- The crack removed only Ca. The sheet O stay, so crack2 is a disc of Ca vacancies, not an
  open slit like crack1. The vacancies survive equilibration: 3 Ca return, compared with 30
  in run1.
- The z-gap is the largest empty z interval through the sheet midplane. The two O layers move
  inward during equilibration, from 2.42 to 1.84 Å, which is slightly below the intact value
  of 1.95 Å. So there is no geometric opening left, only the missing Ca.
- Lz after NPT: 93.08 Å (run1 93.17, crack1 93.05).

**2. Stress–strain.** Plot: `results/compare_crack2.png` (intact, crack1, crack2). Curves:
`results/stress_strain_{intact,crack1,crack2}.dat`. `stress_strain_crack1.dat` holds the same
data as the older `stress_strain_crack.dat`.

| Run | Peak (GPa) | Strain at peak | E₀, ε<2 % (GPa) | σ at 8 % / 10 % | Mean σ, ε=20–30 % |
|---|---|---|---|---|---|
| intact | 11.44 | 0.120 | 124.2 | 8.25 / 10.29 | 7.53 |
| crack1 | 10.09 | 0.142 | 113.5 | 7.07 / 7.73 | 8.98 |
| crack2 | 9.67 | 0.195 | 111.8 | 7.08 / 8.04 | 9.02 |

- crack2 follows crack1 up to about 9 %: the same modulus loss (−10 %) and the same kink near
  8 % at about 7–8 GPa. So in both cracks the kink comes from a weakened Ca sheet / interlayer,
  not from the tilted slit geometry.
- After the kink, crack2 has **no sharp peak**. The stress rises slowly, from about 8.0 GPa at
  10 % to a broad maximum of 9.7 GPa at 19.5 %, and then falls gently to about 8.4 GPa at 30 %.
  crack1 has a small peak at 14 % and then drops. The plateau averages are the same (about 9 GPa).

**3. Damage maps.** Output in `results/damage_crack2/` (not tracked by git; the crack1 maps in
`results/damage/` were left as they were):

```
python scripts/damage_map.py results/damage_crack2 \
  intact=/scratch1/mm18756/tobermolite/run1/traj_compress.lammpstrj \
  crack2=/scratch1/mm18756/tobermolite/crack2/traj_compress.lammpstrj \
  --strains 5,8,12,14,20,30 --crack crack2 --crack-tilt 0 --crack-zc 0.5631
```

| Run | ε | core p99 η | fraction η > 0.2 | core median F_zz | broken Si–O | Ca Δcn ≥ 2 |
|---|---|---|---|---|---|---|
| crack2 | 5 % | 0.15 | 1.6 % | 0.961 | 0 | 8 |
| crack2 | 8 % | 0.26 | 6.4 % | 0.940 | 0 | 24 |
| crack2 | 12 % | 0.50 | 21 % | 0.922 | 0 | 50 |
| crack2 | 14 % | 0.70 | 26 % | 0.916 | 1 | 48 |
| crack2 | 20 % | 1.09 | 38 % | 0.883 | 3 | 80 |
| crack2 | 30 % | 1.71 | 68 % | 0.728 | 6 | 129 |

(The intact rows are the same as in §3.)

- At 5–8 % shear is concentrated in the Ca-vacancy disc. Core p99 η is already 1.7× the
  intact value at 5 % and 2× at 8 %, about the same as crack1 at the kink.
- At 12 % a horizontal band runs across the full width at the crack height (z ≈ 45–60 Å).
  This is the same layer-parallel failure as crack1, but here it starts exactly on the
  weakened sheet.
- The band stays single and widens up to 20 %. By 30 % a second band forms lower down
  (z ≈ 20–40 Å) and the pillar bulges. Core F_zz at 30 % is 0.73, between crack1 (0.86,
  one band) and intact (0.64, two bands).
- At most 6 Si–O bonds break. As in crack1, the silicate chains survive and the damage is
  Ca–O and interlayer rearrangement.
- Caveat: the `shear_xsection.png` slab is at mid-height (z ≈ 46.5 Å), about 6 Å below the
  crack plane. It still cuts the band, but it is not a section through the crack itself.

**Interpretation.** Removing Ca from a single sheet is enough to reproduce crack1's
softening and kink. A slit is not needed. The flat Ca-vacancy disc delays and spreads the
post-yield response (a broad maximum at 19.5 % instead of a peak at 14 %), probably because
the band starts on a plane that is already perpendicular to the load and keeps hardening as
the O layers close up.

## Next steps

- Done (§3): the kink is rim damage that turns into a horizontal band at the crack height.
- Done (crack2): a flat Ca-sheet vacancy disc gives the same kink as the tilted slit.
- Then decide on the deferred tilt/radius sweeps. A crack radius sweep at tilt 0 on a Ca
  sheet would isolate how the defect size controls the yield point.

## rot_penny: rotated pillar (c ∥ x, axis ∥ z) with a Ca-sheet penny

Motivation: so far every pillar had c ∥ z, which puts the layers perpendicular to the load.
In this pillar the crystal is rotated +90° about y, so c ∥ x and the load axis z lies in the
layer plane. The crack2-style penny (a disc of Ca vacancies in one Ca sheet) rotates with the
layers, so its normal is now x. The new code lives in `prep/`. `scripts/build_pillar.py` is unchanged.

- `prep/rotate_cell.py` maps (x, y, z) → (z, y, −x) and writes `prep/tobe_rot.xyz` with rows
  c' = (23.42, 0, 0), b' = (0, 25.80, 16.54) and −a' = (0, 0, 28.01). The pillar axis is crystal a.
  The script checks that the cell volume (16 924.65 Å³), the fractional coordinates and the minimum
  distance (1.0014 Å) are unchanged, and that the Ca sheets are now the planes x = 5.91 and 17.63 Å.
- `prep/build_pillar.py` is `scripts/build_pillar.py` adapted to this cell:
  - The periodic vector is the third row (+z), and b' may have a z component.
  - The axis goes on the Ca-sheet midplane nearest the cell centre (x = 5.91 Å), which replaces
    `--crack-snap-ca`.
  - The crack defaults are tilt 90 and az 0, so the normal is x.
  - It prints R0 for `-var R0`.
- Build (2026-11-10), `/scratch1/mm18756/tobermorite/rot_penny`:
  `prep/build_pillar.py --R 20.0 --nc 4 --vac 40 --crack-a 10 --crack-h 4`
  - Lz = 4 × 28.01 = 112.04 Å, N = 12541 (Ca995 Si1600 O6112 H3834), 442 H added, final charge 0,
    minimum distance 0.970 Å. R0 = 20.69 Å (intact rotated pillar: N = 12512, R0 = 20.67).
  - The penny removed 29 atoms, all Ca. The nearest Si planes are −2.16 / +2.14 Å from the sheet,
    so no SiO4 is cut.
  - Check in `pillar0.xyz`, slab |x − 60| < 2, (y − 60)² + (z − 56.02)² < 100: Ca 0, O 55. The same disc
    shifted 30 Å in z gives Ca 27 and O 56, so as in crack2 it is a Ca-vacancy disc, not an open slit.
    The Ca x-profile shows sheets at x − 60 = 0, ±11.7 and ±23.4 (each one a ±1.5 Å Ca double plane).
- Inputs and SLURM are unchanged. Pass `-var R0 20.69` to both stages.
- Viewing in OVITO: the disc lies in the yz plane through the axis. Look along x, or slice
  with normal (1, 0, 0) at distance 60. Highlight it with
  `abs(Position.X-60) < 2 && (Position.Y-60)^2 + (Position.Z-CellOrigin.Z-0.5*CellSize.Z)^2 < 100`.
- Damage maps: `damage_map.py --crack-tilt 90 --crack-zc 0.5` draws the trace as a vertical
  line at x = 0 in the y-slab view.
- Submitted 2026-11-10 from `/scratch1/mm18756/tobermorite/rot_penny`: equil job 13974972, and
  compress job 13974973 (`--dependency=afterok:13974972`), both with `-var R0 20.69`. Results
  will be copied to `results/rot_penny/`. Before submitting, the stale `/home1/.../tobermolite`
  paths in `slurm/run_pillar.slurm` and `inputs/in.*` were changed to `tobermorite`.
- Both jobs COMPLETED with exit code 0 (equil 2 h 09 m, compress 2 h 18 m). The only warnings
  are the usual ReaxFF force-field ones. Results were copied to `results/rot_penny/`. `sacct`
  SubmitLine shows `-var R0 20.69` on both jobs, and the stress back-computed from
  Pzz·Lx·Ly in `log.compress` gives R0 = 20.690. (The log still echoes `index 20.35`, because
  a command-line `-var` overrides an `index` variable silently.)

### rot_penny results (2026-11-11)

**1. Crack after equilibration.** Script: `/scratch1/mm18756/tmp/rot_open.py` (not backed up).
Region: |x − x_c| < 2 Å, ρ_yz < 10 Å, centred at `zlo + 0.5 Lz`. There is no intact rotated run,
so the reference is the same disc shifted +30 Å in z in the same pillar, which lies on the same Ca sheet.

| State | Lz (Å) | disc slab | reference slab (+30 Å) | x-gap at ρ < 8 Å, disc / reference |
|---|---|---|---|---|
| as built | 112.04 | 53 (O53) | 82 (Ca26 O56) | 2.42 / 2.08 Å |
| minimized | 112.04 | 49 (O49) | 91 (Ca27 Si5 O59) | 2.62 / 2.01 Å |
| after NPT+NVT (frame 0) | 105.39 | 46 (Ca1 O45) | 97 (Ca30 Si4 O63) | 1.92 / 2.03 Å |

- The Ca vacancies survive. One Ca moves back in, compared with 3 in crack2. As in crack2, the two O
  layers close in to about the intact gap, so the result is a Ca-vacancy disc with no geometric slit.
- **Lz shrinks 5.9 % during NPT** (112.04 → 105.39 Å). The c ∥ z pillars shrank only 0.5–0.7 %.
  At the end of the Berendsen stage (fixed box) the axial stress was −2.2 GPa (tension) along a.
  It was −0.74 GPa along c in run1 and −0.80 GPa in crack2. So the layer-parallel direction has a much
  larger residual surface or force-field stress, which NPT relaxes. Compression strain is
  measured from the post-NPT length.

**2. Stress–strain.** Plot: `results/compare_rot_penny.png` (intact c ∥ z, crack2, rot_penny).
Curve: `results/stress_strain_rot_penny.dat`.

| Run | Peak (GPa) | Strain at peak | E₀, ε<2 % | σ at 5 / 8 / 10 / 12 % | Mean σ, ε=20–30 % |
|---|---|---|---|---|---|
| intact (c ∥ z) | 11.44 | 0.120 | 124.2 | 5.36 / 8.25 / 10.29 / 11.12 | 7.53 |
| crack2 (c ∥ z) | 9.67 | 0.195 | 111.8 | 5.08 / 7.08 / 8.04 / 8.26 | 9.02 |
| rot_penny (c ∥ x) | 16.62 | 0.174 | 52.6 | 3.04 / 5.79 / 7.90 / 10.41 | 12.54 |

- The response stiffens strongly with strain (convex). The tangent modulus is 53 GPa at 0–2 %,
  65 at 2–5 %, 90 at 5–8 %, 117 at 8–12 % and 126 GPa at 12–16 %. So E₀ alone understates the
  stiffness. There is no early kink like the 8 % kink in crack1/crack2.
- The rotated pillar peaks at 16.6 GPa (+45 % over intact c ∥ z) at 17.4 %. It then drops sharply,
  from 16.4 at 18 % to 13.3 GPa at 20 %, and levels off at about 12.5 GPa. That plateau is the
  highest of all the runs.

**3. Damage maps.** Output in `results/damage_rot_penny/` (not tracked by git), with an added
17 % frame (near the peak):

```
python scripts/damage_map.py results/damage_rot_penny \
  intact=/scratch1/mm18756/tobermolite/run1/traj_compress.lammpstrj \
  rot_penny=/scratch1/mm18756/tobermorite/rot_penny/traj_compress.lammpstrj \
  --strains 5,8,12,14,17,20,30 --crack rot_penny --crack-tilt 90 --crack-zc 0.5
```

(Takes about 34 s. The intact rows match §3; intact at 17 % is p99 η 0.70, f 38 %, F_zz 0.856, Ca Δcn 75.)

| Run | ε | core p99 η | fraction η > 0.2 | core median F_zz | broken Si–O | Ca Δcn ≥ 2 |
|---|---|---|---|---|---|---|
| rot_penny | 5 % | 0.12 | 1.0 % | 0.951 | 0 | 0 |
| rot_penny | 8 % | 0.19 | 2.6 % | 0.922 | 0 | 0 |
| rot_penny | 12 % | 0.26 | 7.6 % | 0.884 | 0 | 5 |
| rot_penny | 14 % | 0.29 | 13 % | 0.867 | 0 | 12 |
| rot_penny | 17 % | 0.36 | 23 % | 0.840 | 0 | 33 |
| rot_penny | 20 % | 0.69 | 40 % | 0.850 | 4 | 58 |
| rot_penny | 30 % | 1.64 | 57 % | 0.810 | 45 | 75 |

- At 5–8 % shear is concentrated around the penny (the yz disc at x = 0, z ≈ 47–65 Å). Core p99 η is
  1.3–1.4× the intact value. Unlike crack1/crack2, this does not produce a stress kink.
- At 12–17 % the strain stays spread through the pillar, in columns that follow the now-vertical layers.
  Up to the peak, damage is low: Ca Δcn is 5–33 vs 39–75 intact, and f(η>0.2) is 23 % vs 38 % at 17 %.
- After the peak (20 %), a **horizontal band forms across the full width at the penny height**
  (z ≈ 45–70 Å). Because the layers are vertical here, the band cuts *across* the layers. By 30 % the
  pillar bulges at the band (x up to about ±28 Å) and the rest stays nearly elastic (core F_zz 0.81).
- **Silicate chains break.** 45 Si–O bonds are broken at 30 %, compared with ≤ 6 in every c ∥ z run.
  Loading in the layer plane makes the band cut the chains, while in the c ∥ z runs it opened the
  interlayer/Ca–O.

**Interpretation.** With the load in the layer plane, the pillar is stiffer at large strain and about
45 % stronger, and the Ca-vacancy disc no longer causes the 8 % kink. It still seeds the post-peak
band, which forms at the penny height. Failure moves from interlayer collapse (c ∥ z) to a band that
cuts across the layers and breaks Si–O. Caveat: there is no intact rotated run. So the effect of
orientation and the effect of the penny cannot be separated here, and "+45 %" compares against the
intact c ∥ z pillar.

### Next steps (rotated pillar)

- Run the intact rotated pillar (`prep/build_pillar.py --R 20.0 --nc 4 --vac 40`, N = 12512,
  `-var R0 20.67`) so the penny effect can be measured within the c ∥ x orientation.
- Check whether the soft toe region (E₀ ≈ 53 GPa) comes from the 5.9 % NPT contraction along a.
  One test is a longer NPT, or a strain referenced to the as-built Lz.
