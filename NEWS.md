# rastersample 0.1.2

# rastersample 0.1.1

## Documentation

* Updated README.Rmd to suppress build messages and warnings in output
* Updated installation instructions in README
* Minor corrections to tests/README.md

# rastersample 0.1.0

## New Features

### Core Functions

* Added `spatial_sample()` function for flexible spatial sampling with multiple methods:
  - Random sampling via `dplyr::slice_sample()`
  - Biased sampling with threshold filtering
  - Stratified random sampling
  - Conditioned Latin Hypercube (CLH) sampling via `clhs::clhs()`
  - Spatially balanced sampling via `MBHdesign::quasiSamp()`
  - Stratified spatially balanced sampling
  - Transect/line sampling (random and balanced)

* Added `stratify()` function for creating strata from numeric vectors or `SpatRaster` objects:
  - Equal-split quantile-based stratification
  - Custom probability-based stratification
  - Custom value-based stratification

### Input/Output Options

* `spatial_sample()` supports both `data.frame` and `SpatRaster` inputs
* Added `as_raster` parameter to return results as `SpatRaster` objects
* Added `drop_na` parameter for handling missing values
* Added `type` parameter to switch between point and line sampling

### Internal Functions

* `spatialsample_random()`: Simple random sampling
* `spatialsample_biased()`: Biased sampling with threshold
* `spatialsample_stratified()`: Stratified random sampling
* `spatialsample_clh()`: Conditioned Latin Hypercube sampling
* `spatialsample_balanced()`: Spatially balanced sampling
* `spatialsample_transect()`: Transect-based line sampling

## Testing

* Comprehensive test suite using testthat (65+ test assertions)
* Tests for all sampling methods
* Tests for stratification approaches
* Tests for error handling and edge cases
* Tests for NA value handling

## Documentation

* Full roxygen2 documentation for all exported functions
* Examples for common use cases
* Vignette-ready documentation structure

## Bug Fixes

* Fixed issue in `stratify()` where numeric vector inputs were not properly handled
* Added proper handling of `.data`, `cell`, `x`, and `y` variables for NSE evaluation

## Infrastructure

* MIT License
* GitHub repository setup
* testthat framework configured
* R CMD check passes with no errors, warnings, or notes
* Minimum R version: 4.1.0 (uses native pipe `|>`)

---

# rastersample 0.0.0.9000

* Initial development version
* Package skeleton created
