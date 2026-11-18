# Axial splitting of a finite, clamped tobermorite pillar (no penny crack)

Status as of 2026-11-18 (pilot analysed; production not started). Follows `CRACK_PLAN.md` (penny-crack series).

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
| Ends | rigid **clamped** grips (8 Å) only | clamped gives the hotspot. A frictionless-wall control was run in the pilot but **dropped for production (2026-11-18)**, see below |
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
3. **Production** (`split_R60/`): `perp/build_pillar.py --R 60 --L 240 --vac 40` (build on a compute
   node). One equil (2 nodes × 64 ranks), then `damp_finite` (`-var tdamp 800 -var ndamp 48000`), then
   one **clamped** compress job. Account `anakano_429`. (No wall run; `-var ends wall` is still in
   `in.compress_finite` if a control is wanted later.)

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
- Compression (jobs 14323290 clamped, 14323291 wall, about 1 h 26 min each on 32 ranks, both
  finished cleanly at 25 % platen strain). Plots: `results/compare_split_pilot.png`,
  `results/k_vs_strain_split_pilot.png`. Curves: `results/stress_strain_split_pilot_{clamped,wall}.dat`.
  Maps (untracked): `results/stress_split_pilot/`, `results/damage_split_pilot/`.

| Run | Peak (GPa) | Platen strain at peak | Gauge strain at peak | After peak |
|---|---|---|---|---|
| clamped | 19.6 | 14.6 % | ≈ 18 % | drops to 12.0 GPa by 20 %, then 13.7 GPa at 25 % |
| wall | no peak (14.2 at 25 %) | — | 14 % at the end | still rising |

  "Gauge strain" = relative change of the z distance between the centroids of the atoms that start
  in two thin slabs at z_mid ± L/4. It measures the pillar's own strain and leaves out the ends.

**1. The wall control was not loaded like the clamped run.** Platen strain is not pillar strain in the
wall run. The `fix indent` wall (K = 100 kcal/mol/Å³, F ∝ K·δ²) starts at the outermost atom, which
is an end-face H. The first ~3 % of platen travel only flattens the H-terminated end faces: the gauge
strain is 0.3 % at 3 %, 3.6 % at 10 % and 14.1 % at 25 %, and the damage maps show η > 0.3 only in
the top and bottom 2–3 Å. In the clamped run the gauge strain is ≈ 1.27× the platen strain, because
the rigid grips take no strain. At the same gauge strain the two runs carry about the same stress
(clamped 12.2 GPa at gauge 12.7 %, wall 13.7 GPa at gauge 13.6 %). So the bulk responses match, and the
apparent "1.5 GPa" wall modulus is an artifact of contact. The wall pillar never got to the clamped
failure strain (gauge ≈ 18 %). The ν values that `stress_map.py` prints for the wall (ν_x 0.04,
ν_y −0.06 at ≤ 2 %) are meaningless for the same reason. Clamped: ν_x 0.19, ν_y 0.49.

**2. No detectable transverse-tension hotspot.** k values below are the mean ± sd of kcore over frames
with σ_c ≥ 1 GPa:

| Run | platen strain | σ_c (GPa) | kcore_xx (%) | kcore_yy (%) |
|---|---|---|---|---|
| clamped | 3–8 % | 2.3–7.9 | 4.2 ± 3.8 | 7.9 ± 2.3 |
| clamped | 8–14.5 % | 7.9–15.6 | −1.4 ± 1.8 | 3.1 ± 1.4 |
| wall | 15–25 % | 5.0–12.0 | −0.2 ± 2.3 | 3.8 ± 4.0 |

At comparable σ_c the clamped k_xx (normal to the layers) is ≈ 0 and is no larger than the wall
control. The paper's predicted k(L/B = 2, ν ≈ 0.2) ≈ 1–2 % sits inside the ±2–4 % frame-to-frame scatter,
so the pilot could not resolve it, which `stress_map.py` testing already suggested for a 9k-atom pillar.
The `sxx_xz.png` maps show no coherent central tensile region. The red/blue pattern is bin noise, plus
strong vertical stripes along the layers (the stiff silicate chains and the interlayer alternate).
The slightly positive k_yy in both modes is about the same size, so it is not boundary-induced.
`stress_map.py`: `k_vs_strain.png` now hides k where σ_c < `--kmin-sc` (default 1 GPa), because k
diverges as σ_c → 0 (the wall run reached 5800 % at σ_c ≈ 0).

**3. Clamped failure mode: central shear/barreling, not axial splitting.** η rises in a broad central
zone that grows from ~10 % platen strain (η median 0.10 → 0.17 at 16 %), and the grips stay
undeformed. At 14–16 % the high-η region (η > 0.3) is concentrated in the lower-left half of the
y-slab (x < 0, z ≈ 50–90 Å). After the peak the
middle bulges (mid-height radius +7 % at 15 %, +23 % at 25 %; the wall run reaches only +4 %). The
`stress_map.py` density check finds no opened bins (0 everywhere). Bond damage comes late: no Si–O
breaks up to 16 %, 1 at 18 % and 18 at 25 %. Ca–O coordination changes rise from 9 (10 %) to 75 (25 %). So the
drop at the peak is shear/rearrangement in the Ca–interlayer region and not cleavage. That is the
same family as the penny runs' band failure. There is no vertical gap along the layers.

**Takeaways for production.**
- With R = 20 Å, the pilot can't confirm or rule out the hotspot (noise ≳ predicted k). R60 has
  9× the cross-section and 3× the length, so it averages over about 25× more atoms per core bin. That
  should bring the noise below 1 %, but only if the wall control is fixed first.
- **Decision (2026-11-18): drop the wall control and run only the clamped case at R60.** The pilot wall
  run was contact-dominated (see 1). Fixing it would need a stiffer wall or a pre-contact stage.
  Without the control, a hotspot in the clamped run must be judged against that run's own noise. Call
  it real only if all of these hold:
  - kcore_xx is positive and roughly constant over the pre-peak frames with σ_c ≥ 1 GPa (mean
    clearly larger than the frame-to-frame sd), not swinging in sign as in the pilot;
  - z_at, the height of the k_xx maximum, stays near mid-height, and the `sxx_xz.png` maps show one
    coherent central tensile region;
  - the size is above the `stress_map.py` synthetic noise floor for that atom count, and is
    compared with the paper's k(L/B = 2, ν) using the clamped ν_x, ν_y.
  If the result is ambiguous, rerun a wall control with a stiffer wall (`-var Kw 1000`+) placed at
  the mean end-face height.
- Also report gauge strain as well as platen strain (the grips make platen strain ≈ 0.8× gauge).
- At R = 20 the shear-band mode wins in the clamped run before any splitting. Whether a larger
  diameter changes that is the actual question for R60 (Östlund / Yuan size effect).
