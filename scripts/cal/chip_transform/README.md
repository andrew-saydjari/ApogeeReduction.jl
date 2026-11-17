# Comb-measured chip transform: rebuilding the calibration before a release

## What it calibrates

Per fiber, a wavelength solution is one quartic in the green-chip coordinate shared by the
three chips; chips R and B are mapped onto it by `x' = a_c + b_c * x`
(`transform_x_chips`, `x = (pixel - 1024) / 2048`).

- The **sky-line fit** (`get_and_save_sky_wavecal`, cporder 0) frees the offsets `a_c` and
  forces the chip scales `b_c = 1`.
- The **FPI comb fit** (cporder 1) measures `b_c` from ~500 uniform peaks per chip. The scales
  are not 1 (APO `b_R - 1 ~ -2.6e-4`, `b_B - 1 ~ +1.1e-4`; LCO `-0.9e-4`, `-0.4e-4`). They are
  stable to ~1e-5 for hundreds of nights and have a reproducible per-fiber pattern.
- When `b = 1` is forced, the shared quartic absorbs the scale error as 65-95 m/s of per-chip
  curvature.

The calibration file holds the comb's chip scales per fiber, by default in one all-time block
per telescope. The sky fit (`pipeline_2d_1d.jl --chip_transform_cal <file>`) fixes `b_c` to the
file's values instead of 1:

- the offsets stay free;
- the number of fitted parameters is unchanged;
- the stored `nlParams` take the cporder-1 layout `(a_R, b_R, a_B, b_B)`.

**It is ON by default.** The pinned default is
`/mnt/ceph/users/sdssv/work/asaydjari/cal_ref/chip_transform/chiptransform_pass1dr21_v2_allblock.h5`,
set in `pipeline_2d_1d.jl`'s arg table and in `airflow/dags/ar_common.py` (`CHIP_TRANSFORM_CAL`).

- **Off switch:** `--chip_transform_cal ""`, `CHIP_TRANSFORM_CAL = ""`, or
  `AR_CHIP_TRANSFORM_CAL=""` for the SLURM scripts. This gives chip scale 1, exactly the
  behaviour from before the calibration existed.
- **Hard errors, never a silent fallback to scale 1:**
  - a configured file that is missing or unreadable;
  - a file with no block for the exposure's telescope and MJD.
- **Checkpointing:** each `waveCalSkyLine` product records `chip_transform_id`
  (`<version>:<FNV-1a fingerprint of the file>`, or `none`). A stored product is reused only if
  it was made with the same file. Changing the file (or turning the feature on or off)
  recomputes the sky wavecal even under `commit_same` / `commit_exists`.
  - Products downstream of it (`waveCalNightskyDither`, `waveCalFPI`, `waveCalNightfpiDither`
    and the reinterpolated `ar1D`) are checkpointed on the code commit only. After changing the
    calibration file at the same code commit, re-run `pipeline_2d_1d.jl` with
    `--checkpoint_mode clobber` so they are rebuilt from the new sky solutions.

Nothing learned is written into the source. Every number comes from the comb fits of a
reduction, through `build_chip_transform.jl`, into a versioned HDF5 file that carries its own
provenance:

- code commit and tree state;
- input root and MJD range;
- block rule and all its parameters;
- blocks and counts;
- coverage choice;
- creation date.

## Stages (`build_chip_transform.jl --stage ...`; each is resumable)

| stage | what | output in `--workdir` |
|---|---|---|
| `collect` (1) | Per night with a nightly wavelength file: the comb fit's per-fiber `(a_R, b_R, a_B, b_B)`, comb peaks used per fiber/chip, and residual RMS. Nights the FPI gate refused are kept and marked. With `--run_missing true`, nights that have FPI exposures and no gate verdict are re-fitted with the production comb fit (`FPI_COMB_FIT_KWARGS`, same gate). The fit writes only into `--workdir`; the reduction root is never modified. | `stage1/<tele>/<mjd>.h5`, `stage1_<tele>_nights.csv` |
| `derive` (2) | Robust per-fiber chip transform, one all-time block per telescope by default. See the rules below. | `stage2_<tele>.h5`, `stage2_<tele>_{blocks,night_series,rejected,fiber_model_cv,detected_change_points}.csv`, `qa/qa_<tele>_*.png` |
| `write` (3) | The calibration file with provenance. `--cover_from` / `--cover_to` extend the first and last blocks, e.g. back over the plate era. | `<version>.h5` |
| `check` (4) | Re-evaluates sky (and, with `--arc_lines_dir`, arc) solutions on a sample of nights under four hypotheses: `scale1` (production sky), `cal_scale` (this change), `cal_scale_off` (scale and offset from the file), and `free_scale` (cporder-1 freedom, production arcs). Metrics are on held-out lines (exact leave-one-line-out), zero-point bias, and the shape against the comb where one exists. | `stage4/<tag>/`, `stage4_<tag>_{per_night,summary_vs_scale1,summary_vs_free_scale,flat_rv}.csv` |

### Stage 2 rules

All parameters are command-line flags. They are recorded in the file.

- **Fiber cut on a night.** A fiber is dropped on a night if any chip has fewer than
  `--fiber_min_peaks` comb peaks used. It is also dropped if its scale is more than
  `--clip_nsigma` robust sigma off that night's quadratic-in-fiber trend.
- **Night cut.** A night is dropped if fewer than `--night_min_fibers` fibers remain. A night
  is also dropped if its median scale is a `--clip_nsigma` outlier within its block.
- **Blocks.** The default is `--block_mode none`: one all-time block per telescope.
  - Change points of the scale are detected and reported in every mode
    (`stage2_<tele>_detected_change_points.csv`, dotted lines in the QA plot, and the
    file's provenance). They are kept for matching against a shared time-block definition
    (traces and other calibrations).
  - `--block_mode detect` uses them as block boundaries. Detection is binary segmentation on
    the nightly fiber-median `b_R` and `b_B`:
    - noise comes from successive differences;
    - outlier nights are screened against a running median;
    - a split is accepted if both sides have at least `--min_block_nights` nights, the step
      is significant at more than `--block_nsigma`, and its size exceeds `--block_min_step`.
  - Boundaries can also be supplied (`--block_mode supplied` or `supplied+detect`,
    `--boundaries`), e.g. from a shared time-block definition.
  - Step changes of the offsets are reported as hardware events (QA only).
- **Per block.** For each fiber, the median over nights after `--clip_nsigma` clipping, with
  its standard error.
- **Fiber model** (`--fiber_model auto`). A split-half cross-validation over the block's nights
  estimates the error of the per-fiber medians and of a smooth polynomial in fiber
  (`--fiber_poly_degree`). The model with the smaller error on both chips is used. Fibers
  without usable comb data get the smooth value and are flagged `fiber_filled`.

## Running a full rebuild before a release

```bash
cd <ApogeeReduction.jl checkout at the release commit>
# local (one node): N single-threaded workers
bash scripts/cal/chip_transform/run_chip_transform.sh <reduction root> <workdir> chiptransform_<release> 6 \
     --cover_from <release start MJD> --cover_to 99999 \
     --cover_note "one all-time block extended over the plate era; plate-era check: <path to stage-4 tables>"
# Slurm (every task = one worker; the script holds the #SBATCH header)
sbatch scripts/cal/chip_transform/run_chip_transform.sh <reduction root> <workdir> chiptransform_<release> 0 --cover_from ...
```

Then:

1. Read `qa/qa_<tele>_night_series.png` (the detected change points against hardware events),
   `stage2_<tele>_detected_change_points.csv` and `stage2_<tele>_rejected.csv`.
2. Read `stage4_check_summary_vs_scale1.csv`. Held-out line RMS must not get worse, and the
   `cal_scale` zero-point shifts are the expected bias.
3. Install `<version>.h5` in `/mnt/ceph/users/sdssv/work/asaydjari/cal_ref/chip_transform/`.
   Then point the defaults at it:
   - `--chip_transform_cal` in `pipeline_2d_1d.jl`'s arg table;
   - `CHIP_TRANSFORM_CAL` in `airflow/dags/ar_common.py`.

   Never overwrite an installed version; install a new file name. Products record the file,
   version and fingerprint in their `metadata`.

A plate-era validation set (nights without a comb, plus FPS controls) is passed as
`--check_nights <file>` with lines `<tele> <mjd> <category>`. Arc lines are passed as
`--arc_lines_dir` holding `arclines_<tele>_<mjd>.h5` (datasets `fiber, chip, group, pix, wave`
of identified lamp lines).

## Measured cost (pass-1 DR21, ccalin051, 3 single-threaded workers)

Measured on the full pass-1 DR21 archive:

- 3733 nights: APO 2084, LCO 1649.
- 2111 comb nights used: APO 1046, LCO 1065.
- 196 nights were refused by the FPI gate.

| stage | wall time | peak memory |
|---|---|---|
| 1 `collect`, cold | ~9 min (APO 4.0 + LCO 4.6) | 3.2 GB main + ~1.5 GB per worker |
| 1 `collect`, warm (all nights cached) | 2.6 min, 1 worker | 1.9 GB |
| 2 `derive`, both telescopes | 2.1 min, serial | 1.9 GB |
| 3 `write` | 11 s | |
| 4 `check`, 320 nights (up to 4 sky exposures and all arc sequences per night) | 9.2 min (~5 s per night per worker) | 1.9 GB main |
| comb re-fit of one missing night (`--run_missing`) | 116-138 s per night, single-threaded | |

- The comb re-fit reproduces the stored pass-1 chip scales to <= 7e-12.
- Stages 1 and 4 are embarrassingly parallel over nights and scale with `--workers` or the
  number of Slurm tasks. Stages 2 and 3 are serial and cheap.
- A full rebuild from stored comb fits is therefore ~25 min on one node with 3 workers.
