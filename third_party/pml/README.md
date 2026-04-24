# PML - FLINT Extras (bundled for DRsolve)

This directory contains the bundled PML library (FLINT extras) used by DRsolve. The original upstream project can be found at: <https://github.com/vneiger/pml>

DRsolve utilizes the PML library to accelerate univariate polynomial matrix determinant computation over prime fields.

## About PML

PML is a supplementary library for FLINT (Fast Library for Number Theory) that provides additional functionality including:
- Extended CRT (Chinese Remainder Theorem) operations
- Enhanced matrix and polynomial operations
- FFT-based operations
- Additional modular arithmetic utilities

## DRsolve-specific modifications

The main functional change in this bundled copy is in `src/nmod_poly_mat_extra/nmod_poly_mat_det.c`:
- added `nmod_poly_mat_det_generic`
- added `nmod_poly_mat_det_generic_knowing_degree`
- kept determinant computation usable from DRsolve for generic polynomial matrices over prime fields

There are also a few small compatibility-oriented header/config adjustments so the bundled copy builds cleanly inside DRsolve, but those are secondary to the determinant changes above.

## Building

This bundled copy is intended to be built as part of the DRsolve project. For standalone compilation, please refer to the original PML documentation.

## License

See the [COPYING](COPYING) and [COPYING_FLINT](COPYING_FLINT) files for licensing information.

## Original Authors

See the [AUTHORS](AUTHORS) file for the original PML contributors.
