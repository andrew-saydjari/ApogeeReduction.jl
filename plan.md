# Code Development Plan

## Big Picture Questions
- I am leaning towards a slight refactoring, where we run everything per chip. It would really help with overheads. It would just mean we launch 6 jobs nightly. [agreed]
- Related to the first+next question. When do we want the chips combined?
- When do we split from all fibers in exposure to single object files [at dither combination]
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
- speed up (this is currently bottlenecking my development on the darks)
- compute MJD-mid at this stage based on firstind selection

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
- should we be evolving the LSF?
- bad pixel mask evolution
- sky flux variation over field of view
- fix SlackThreads.jl for new API [completed, AKS]
