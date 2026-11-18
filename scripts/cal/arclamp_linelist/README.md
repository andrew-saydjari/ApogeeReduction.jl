# Arc-lamp line list (stage-1 arc wavelength calibration)

`data/arclamp_linelist_v1.csv` is the frozen ThAr + UNe line list read by
`load_arclamp_linelist` (`src/arclamp_wavecal.jl`). It is built by `build_arclamp_linelist.jl`
in this folder. The file is versioned: its FNV-1a fingerprint and the selection used are
recorded in every `arclampWaveCalNight_*` product. Changing the list or the selection moves the
arc's absolute zero point by tens of m/s per chip, so a change is a release-level decision.

## What is in the file

One row per recurring arc feature (lamp × chip). Columns:

- `assoc_wave`: the wavelength used to associate peaks with the feature. This is the lab
  wavelength for candidate fit lines and the empirical wavelength for everything else.
- `med_pix`: the feature's median pixel. It is used only for the seed polynomial.
- `lab_id`, `lab_wave`, `lab_species`, `lab_class`, `lab_status`: the matched laboratory line.
  Wavelengths are vacuum Å from the clean tables (NIST SRD 161 / Redman 2014 Th, Redman 2011
  U/Ar, NIST ASD Ne), never the apogee_drp copies.
- Quality and selection flags:
  - `store`: lab-matched. These rows are kept in the line table.
  - `shape_clean`, `shape_blend`, `strict_clean`: the data-derived line-shape and stability
    classes.
  - `prelim_classB`: the September lab class-B list.
  - `drp_rederived`: one of the 8 lines whose wavelength apogee_drp re-derived from data in 2021.
    These are blends.
  - `neon`.
- Fit selections:
  - `use_default`: shape-clean Th/U/Ar, without the blended lines, without neon. 93 lines:
    R 20 / G 21 / B 52.
  - `use_with_blended`: 96 lines.
  - `use_with_neon`: 127 lines.

  Switching selection is one argument:
  - `load_arclamp_linelist(; selection = "use_with_neon")`;
  - or `pipeline_2d_1d.jl --arclamp_linelist_selection use_with_neon`.

  A measured per-lamp neon offset, from a calibration file still to be written, enters as
  `load_arclamp_linelist(; selection = "use_with_neon", neon_offset_A = Dict("ThAr" => ..., "UNe" => ...))`.

## Rebuild

```
julia --project=. scripts/cal/arclamp_linelist/build_arclamp_linelist.jl \
    --empirical <empirical_linelist.csv> --labfield <labfield.csv> --lineset <lineset.csv> \
    --fitsets <fitsets.csv> --status <line_status_pooled.csv> \
    --out data/arclamp_linelist_v1.csv --version v1
```

The v1 inputs are listed in the file header, with their sha256. They are:

| input | produced by (offline study, 2026-11-16/17) |
|---|---|
| empirical feature catalog (`v2_N8_clean`) | local-noise finder at 8σ over 44 archive nights; features recurring across nights, fibers and both telescopes |
| `labfield.csv`, `lineset.csv` | the clean laboratory tables, restricted to the band and the two lamps |
| `fitsets.csv`, `line_status_pooled.csv` | line-shape and stability classification of every lab line on 40 nights (3.74 M line × fiber measurements): width against the local LSF, stacked-profile residual and asymmetry, night-to-night stability, offset |

**Known gap.** The empirical catalog and the line classification are data-derived. The scripts
that re-derive them from the full dataset are still the offline study scripts. They are not yet
in this repository. Porting them here is stage 2 of the design (`ARCLAMP_PRODUCTION_DESIGN.md`
§3), and it must merge before the arc solution is switched on. Until then this builder is
deterministic given its recorded inputs, but it does not start from raw data.

## Without the file

The arc fit needs a line list. The file is in the repository, so it is always present. The arc
product is off by default (`--arclamp_wavecal false`), and nothing reads it yet.
