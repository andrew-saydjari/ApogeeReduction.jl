# Code Development Plan

## Introduction
- Purpose of the code
- Overview of the code base

## Big Picture Questions
- When do we split from all fibers in exposure to single object files
- What tuple do we want to use as an identifier (file naming)? Currently using (tele,mjd,chip,expid).
- What foldering do we want? apMADGICS folders per fiber, because that is how the priors should be loaded at that stage. Nightly processing lends itself to MJD foldering.
- What information gets carried where? Do we carry telescope header information to each stage? What meta information is in a lookup table and what needs to follow spectra?
- In each subheading, maybe we should list the files we intend to have as inputs/outputs?

## Almanac Level
- Tracking all raw files on disk
- Fiber mapping/target information
- Bad exposure listing (prefilters, updated based on outlier detection from the pipeline)
- Final QA pass to ensure all files exist/are not corrupted
- Make flowchat showing number of exposures/spectra dropped and why
- Identify changes in file structure, flag to datamodel (datamodel automation)

## 3D -> 2D
- Better zeropointing [initial complete, AKS]
- 1/f correction [chip a APO complete, AKS]
- reference array based masking/correction
- nonlinearity corrections
- correcting the first/second read

## 2D Cal
- darkRate subtraction [in progress, AKS]
- dark bad pix mask [initial complete, AKS]
- flat fielding (no clear route, internal flats are not flat)
- flat bad pix mask

## 2D -> 1D
- find traces [in progress, KM]
- tweak traces to exposures
- extract flux
- PSF model
- LSF model

## Wavelength Cal
- fit FPI and arclamp peaks
- fit FPI cavity [initial complete, AKS]
- fit wavelength solution tweaks

## Fluxing
- relative fluxing so "sky" has same flux key
- fiber fluxing might have a multiplicative component for each connector stage
- use transfer function from apMADGICS to handle color dependence
- learn single flux-ing number per night to convert to physical units
- consider fiber flux versus PSF fluxing for sky fibers
- consider tellurics-based calibration (AMO motivated from CKJ)

## Tellurics
- apMADGICS needs to decouple its prior from the current DRP telluric model outputs

## QA
- wavelength calibration stability
- fluxing stability
- darks clean
- flats clean
- bad pixel mask evolution
- sky flux variation over field of view