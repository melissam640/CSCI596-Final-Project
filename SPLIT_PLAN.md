# Axial splitting of a finite, clamped tobermorite pillar (no penny crack)

Status as of 2026-11-17. Follows `CRACK_PLAN.md` (penny-crack series).

## Motivation

None of the penny-crack runs (crack1, crack2, rot_penny, rot_penny_x10) showed brittle **axial
splitting**, the dominant compressive failure mode of brittle solids. Each failed by a horizontal
band at the defect. The papers in `/scratch1/mm18756/papers` point to three changes:

- **Khodadad, Barthelat, Clayton, Gazonas, Dayal, PRL (2026), "Compressive splitting in brittle
  solids: the inverse of wrinkling in sheets".** Axial splitting is a deterministic, flaw-free
  effect of the **loading boundaries**. Clamped platens stop the Poisson expansion at the ends, which
  creates a central tensile "hotspot" of transverse stress σ_xx > 0 (their Fig. 2a). A Mode-I crack
  parallel to the load starts there. The hotspot is k = σ_xx,max/σ_c ≈ 0.5–3 % for ν = 0.2–0.3.
  For low ν it is largest near L/B ≈ 1.7, and slender specimens with L/B ≥ 3 converge. The tension
  disappears for ν = 0, for frictionless ends and under lateral confinement (σ_conf/σ_c > k).
  Splitting strength is C = σ_t/k(L/B, ν).
- **Yuan, Nomura, Nakano, APL (2012)** (ZnO nanowires) and **Östlund et al., AFM (2009)** (Si pillars)
  found that brittle cleavage needs a **large diameter and low temperature**. Small or hot pillars go
  through structural transformation or ductile flow instead. Östlund's > 400 nm pillars split along
  the axis, starting in the interior.
- **Aditya et al., JPCL (2023)** (MoSe₂ wrinkles) is the tension–compression mirror image of the PRL
  mechanism: constrained lateral contraction under tension leads to wrinkling.

All earlier pillars here were **periodic along z** and compressed with `fix deform`. They are
infinitely long with no ends, so this boundary-induced tension could not appear in them. The penny was
also a flaw-based trigger, which the PRL argues is not needed.

## Design (decided 2026-11-17)

| Choice | Value | Reason |
|---|---|---|
| Penny crack | none | test the flaw-free, boundary-driven trigger |
| Orientation | c ∥ x (layers vertical), axis = crystal a = z, `perp/` builder | the weak interlayer is parallel to the load, so σ_xx (normal to the layers) can open an axial split along it |
| Geometry | finite cylinder, L/D = 2 (R = 60 Å, L = 240 Å, about 225k atoms) | standard specimen aspect ratio, near the k peak; large diameter for brittle cleavage |
| Ends | rigid **clamped** grips (8 Å), plus a **frictionless-wall** control | clamped gives the hotspot; frictionless should remove it |
| T, rate | 10 K, 1 %/ps (as before) | low T for brittleness; rate kept comparable with the earlier runs |
| Force field | ReaxFF (Manzano 2012) | bond breaking, comparable with the earlier runs |

## Code changes

- `perp/build_pillar.py --L <Å>`: finite (non-periodic) pillar. Whole units are kept with r < R and
  |z| < L/2. The charge balance protonates the outermost non-bridging O by **depth to the nearest free
  surface** (side or end face), and end-face H point along ±z. `pillar.data` has a z box of
  [0, L + 2·vac] and is meant for `boundary s s s`. `min_dist` now uses a numpy cell list (scipy is
  not installed), and the periodic mode gives the same atoms as before (`--R 20 --nc 4`: N = 12512,
  identical atom set).
- `inputs/in.equil_finite`: `boundary s s s`. CG, then Berendsen 5 ps, then NVT 10 ps at 10 K. No NPT
  (there is no periodic direction; the free ends relax the axial surface stress). Writes `thermal.restart`.
- `inputs/in.compress_finite`, with `-var ends clamped|wall` and **`-var R0` required**:
  - clamped: atoms within `grip` = 8 Å of each end form rigid grips. The bottom grip is frozen
    (`fix setforce` reports its load). The top grip moves with `fix move linear 0 0 v`, v = erate·L0,
    with no lateral motion. NVT applies to the middle only.
  - wall: two frictionless `fix indent … plane z` walls. The top wall moves at v. NVT on all atoms,
    with lateral momentum and rotation removed.
  - Load = ½(F_top − F_bot) (forces on the platens; signed, compression > 0, tension < 0),
    σ = F/(πR0²) in GPa, strain = platen displacement / L0. (Was ½(|F_top| + |F_bot|) until the
    debug run showed this hides tension; fixed 2026-11-17.)
    Thermo columns: `v_strain v_ftop v_fbot v_scomp`.
  - Stage A: an `nhold` = 2000-step hold with no load, which gives the reference stress. Stage B:
    `ncomp` = 100000 steps (25 %).
  - Dumps: `traj_compress.lammpstrj` every 400 steps (0.1 %), for η. `stress_compress.lammpstrj`
    every 2000 steps (0.5 %) with `f_sav[1..6]` = `stress/atom` averaged over the last 20×10 steps.
    The first stress frame is the unloaded reference.
- `slurm/run_pillar.slurm`: the copy-back matches `equil*` / `compress*`, so the new stage names are
  backed up too.
- `scripts/stress_map.py` (new): Δσ = σ(ε) − σ(0) per atom, binned in slabs (σ_bin = n_bulk·⟨Δs_i⟩,
  3×3 smoothing), normalized by σ_c (= −⟨Δσ_zz⟩ over the mid-height cross-section). Writes
  `sxx_xz.png` (paper Fig. 2a analog: y-slab, zero contour, hatched opened bins), `syy_yz.png`,
  `sxx_xsection.png`, `k_vs_strain.png`, `k_<label>.dat` and `stress_<label>_<strain>.xyz` (OVITO).
  - The robust k is `kcore` = the peak along z of Δσ_xx averaged over the core r < R/2. The
    slab-bin maximum is also reported but is noisier.
  - Also prints apparent Poisson ratios ν_x, ν_y (mid-height rms lateral strain over axial strain,
    ≤ 2 %), to compare with the paper's k(L/B = 2, ν).
  - Tested on synthetic dumps: an imposed 2 % σ_xx hotspot is recovered as kcore = 1.6–2.5 %, with a
    noise floor (k_yy with no signal) of 0.1–0.8 % at 0.1 σ_c per-atom noise. With very heavy noise
    (0.5 σ_c per atom) on the 9k-atom pillar the floor is about 3–5 %. So **always compare the
    clamped k with the wall-control k**, not with zero.
- `scripts/damage_map.py`: detects `BOX BOUNDS ss ss ss` and then applies no z minimum image. Adds
  `--strain-per-frame`. Output on rot_penny is unchanged (checked at 8 and 20 %).

## Run sequence

1. **Debug** (`split_debug/`, `debug.slurm`, debug partition): R = 20, L = 80, N = 8969, R0 = 20.71 Å.
   Shortened equil, then 0.5 % of compression in each end mode.
2. **Pilot** (R = 20, L = 80, full 25 %, clamped + wall, 32 ranks). Check that the hotspot appears in
   the clamped run only.
3. **Production** (`split_R60/`, `split_R60_wall/`): `perp/build_pillar.py --R 60 --L 240 --vac 40`
   (build on a compute node). One equil (2 nodes × 64 ranks), then `damp_finite` (`-var tdamp 800
   -var ndamp 48000`), then two compress jobs from the same `thermal.restart`. Account `anakano_429`.

## Viewing in OVITO

- `stress_<label>_<strain>.xyz`: color by `dsxx` (GPa). Add *Spatial binning* (2D, x–z, slab
  |y − 60| < 6) to reduce per-atom noise.
- Split: in `traj_compress.lammpstrj` look along y (Front viewport) with a y-slab at the axis. An
  axial split along the vertical layers appears as a vertical gap through the centre, parallel to z.

## Results

### Debug (`split_debug/`, jobs 14317303 + 14319281, 2026-11-17)

- 14317303 ran the equilibration and then failed at compress: `reset_timestep` after `fix ave/atom`. Fixed by moving it
  before the averaging fix. 14319281 ran both 0.5 % compressions without errors (about 23 steps/s on
  32 ranks, 8969 atoms). It wrote 11 trajectory frames + 11 stress frames per mode, and `stress_map.py` ran
  on both (needs `OPENBLAS_NUM_THREADS=1` on the login node).
- The shortened equil (188 CG steps + 0.2 + 0.2 ps) left the pillar still relaxing: T 40 → 18 K and PE
  falling. The free (wall) pillar shrank about 1 Å while the wall moved 0.4 Å, so the wall barely touched it. The clamped
  pillar went into about 2 GPa tension within 50 fs of the hold. The stress maps are relaxation drift
  (k 100–600 %, ν_x ≈ 3.8, ν_y ≈ −4.7, the same in both modes), not response to the load.
- That tension exposed the `abs()` load bug (fixed above).
- For the pilot: use the full equil and check that PE, T and the pillar length plateau before compressing.

### Pilot (`split_pilot/`, R = 20, L = 80, 8969 atoms, R0 = 20.71 Å)

- Equil (job 14320574, 54 min, 32 ranks): CG converged on energy tolerance. NVT held T = 10.0 K, and PE
  drift fell to about 7 kcal/mol per 2.5 ps. The mean length plateaued at 79.4–79.7 Å (82.3 Å at the start of MD).
  The pillar was still **ringing in its fundamental axial mode**: period about 3.3 ps, sd 0.38 Å
  (about ±0.7 % strain) at 15 ps, decaying slowly. Berendsen and NVT do not damp this coherent mode.
  The restart was taken at Lz = 79.30 Å (−0.5 % from the mean).
- New stage `inputs/in.damp_finite` (job 14321933, 20 min): 4 ps Langevin (τ = 250 fs, near-critical
  for this mode), then 2.5 ps NVT. It overwrites `thermal.restart`; the undamped one is kept as
  `thermal_nvt.restart`. Result: the length sd per ps is 0.02–0.04 Å, the 1-ps means over the last 3.5 ps
  stay within 79.50–79.61 Å, T = 9.97 K and PE is flat. The coherent mode is gone. **For R60, use
  `-var tdamp 800 -var ndamp 48000`** (P ≈ 10 ps).
- Compression: (to be filled in)
