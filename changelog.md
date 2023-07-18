# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.8]

### Fixed
- missing data_range specification in structural_similarity call

## [0.3.7]

### Removed
- Info about `chess filter` in the cli help.

### Changed
- `chess pairs` now produces 0-based ranges.

## [0.3.6]

### Fixed
- All NaN results in `chess sim` runs on cooler files, fixed through fanc 0.9.9.
- Cooler file loading time via fanc, significant speed up in fanc 0.9.8

### Changed
- Switch to matplotlib 'pdf' backend

## [0.3.5]

### Added
- Raise error if too small (<20 bins) regions are used in `chess sim` before running comparisons.

## [0.3.4]

### Fixed
- "ValueError: Image must contain only positive values" in `chess extract`

## [0.3.3]

### Added
- Explicit error for zero valid pairs in the pairs bed input for `chess sim`

## [0.3.2]

### Added
- Do output directory checks in chess extract and chess crosscorrelate

### Removed
- `chess filter`

## [0.3.1]

### Added
- Check whether output directories exist before running any computations.
- Allow paths to output files without basenames.

## [0.3.0]

### Added
- `--window-size`, `--sigma-spatial`, `--size-medianfilter` and `--closing-square` optional parameters for more control in the feature extraction.
- Accept data in FANC, Juicer and Cooler formats (.hic, .cool, .mcool)
- Raise a warning message if the pairs file contains chromosomes that are not specified in the provided contact data.

### Fixed
- `chess pairs`: Fix error on empty lines in chromosome sizes file.
- `chess pairs`: Catch OSError raised by pybedtools when the provided path is not recognized as a UCSC genome id.

## [0.2.0]

### Added
- `chess extract`: Extract specific regions that are significantly different
- `chess crosscorrelate`: Get structural clusters from the extracted submatrices