# Coherent detector-defect region finder

An **offline QA tool**.  It reads calibration and 2D products from a completed
reduction and reports coherent defect *regions*.  It does not import
`ApogeeReduction`, does not change any reduction behaviour, does not add bits to
`bad_pix_bits`, and writes nothing into a reduction directory.

## Why it exists

`bad_pix_bits` is built entirely from per-pixel thresholds: dark rate too high or
negative, flat response below `flat_frac_cut`, SUTR chi2, saturation.  Those cuts are
blind to a defect where no individual pixel is far out of range but a contiguous block
of them is consistently, mildly off.

The worked example is APO chip G, columns ≈512–531, rows ≈1362–1377.  The block is
raised ~4 % in `flat_im`, slightly *negative* in `dark_rate`, and has its counts
destroyed in dome and quartz `dimage`.  `pix_bitmask` there reads 8
(`pix_not_dark_corr_bits`), which is not in `bad_pix_bits`; 1.3 % of the surrounding
box carries any bad bit and none of it is inside the block.  The traces of fibres 201
and 202 cross it, the notch survives into `ar1Dcal` and `ar1Dunical` with
`mask_1d = GOOD` and negative flux from a flat lamp, and 84 % / 78 % of those two
fibres' entire chi2 comes from the corresponding wavelength band.

## Method

For each (telescope, chip) and each product:

1. **Robust detrend along the dispersion direction.**  Per row, a piecewise-linear
   interpolation of 64-column block medians models the smooth level; the noise scale
   comes from first differences within 256-column blocks, which is immune both to the
   smooth trend and to a flat-valued defect inflating its own sigma.
2. **Winsorise at z = ±4.**  Without this a single hot pixel (z up to 2.7e4 in a real
   `dark_rate` map) dominates every box that contains it.  Isolated extremes are the
   existing per-pixel cuts' business; this tool is about coherence.
3. **Subtract the 129×129 local mean.**  Removes edge roll-off and broad gradients,
   which are real but are not defect regions.  A block the filter bank can resolve
   loses only a few percent of its amplitude.
4. **Rectangular matched filter**, `Σz/√N`, at 13 scales from 3×3 to 65×5 — square
   scales for blobs, thin ones for bad columns and rows.
5. **Standardise each scale's output against the rest of its own row.**  The analytic
   unit-variance normalisation is not usable on real detector data: measured on an APO
   chip-G dome flat it exceeds 6 over 59 % of the array.
6. **Threshold |Sn| ≥ 6**, 8-connected components, minimum area 6 px.
7. The 2D-flat channel repeats this per epoch-frame and keeps only pixels firing in
   **≥ 2 independent frames**, so a cosmic ray cannot survive.

Regions are then graded by how many independent products confirm them
(`dark_rate`, `flat_im`, the 2D flats) and screened for compactness, and each is
reported with its MJD extent, the fibre traces that cross it, and the fraction of it
that AR already flags.

## Coordinates

`dimage` / `ivarimage` / `pix_bitmask` / `dark_rate` are `(2560, 2048)`; columns
2049:2560 are the reference array.  `flat_im` is `(2040, 2040)` and corresponds to
`dimage[5:2044, 5:2044]` (`src/ar3D.jl`).  1D extraction indexes `dimage[xpix, ypix]`
directly, so the 1D column index is the detector column and `extract_trace_centers` is
in detector rows.  `extract_trace_coords[:,:,1]` in `ar1Dunical` is the *stacked* x,
with `detector column = 2049 - x_stack`, and `[:,:,3]` is chip 1=R, 2=G, 3=B.

Everything is carried internally in science coordinates `1:2040` and **reported in
detector coordinates** (`detector = science + 4`).

## Running

```
AR_APRED=<outdir>/apred AR_QA_OUT=<workdir> julia --project=. run_sweep.jl
AR_APRED=<outdir>/apred AR_QA_OUT=<workdir> julia --project=. null_test.jl
AR_APRED=<outdir>/apred AR_QA_OUT=<workdir> julia --project=. injection_test.jl
```

`NMJD` sets how many MJDs the epoch channel samples (default 14).  The sweep is
single-threaded, needs ~1 GB, and takes about a minute per telescope-chip at
`NMJD=14`; it is comfortably a login-node job.

Outputs: `regions.jls` (full region records), `masks_<tele>_<chip>.h5` (detection
masks and evidence maps), plus `NULL_TEST.txt` and `INJECTION_TEST.txt`.

## Validation

* `null_test.jl` re-runs the whole chain on the same residual map with each row
  shuffled — marginal distribution preserved exactly, all spatial coherence destroyed.
  Anything reported there is a false positive by construction.
* `injection_test.jl` adds synthetic blocks of known size and amplitude at locations
  the finder calls clean and measures the recovered fraction, giving a completeness
  surface in (size, per-pixel amplitude).

## Status

The masks this produces are a **proposal**, not a pipeline change.  Nothing here is
wired into a reduction.
