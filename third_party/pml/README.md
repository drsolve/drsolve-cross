# PML - FLINT Extras (Windows Compatibility Fork)

This is a modified version of the PML library (FLINT extras) with minor adjustments to support Windows compilation. DixonRes uses the PML library to accelerate univariate polynomial matrix determinant computation over prime fields.

## About PML

PML is a supplementary library for FLINT (Fast Library for Number Theory) that provides additional functionality including:
- Extended CRT (Chinese Remainder Theorem) operations
- Enhanced matrix and polynomial operations
- FFT-based operations
- Additional modular arithmetic utilities

## Modifications for Windows Compatibility

The following changes have been made to enable compilation on Windows:

### 1. Machine Vectors Disabled (machine_vectors.h)
Completely disabled machine vector operations (AVX2, AVX512) to avoid compilation issues on Windows:
- All vector type definitions and operations are commented out
- `PML_HAVE_MACHINE_VECTORS`, `PML_HAVE_AVX2`, and `PML_HAVE_AVX512` are undefined

### 2. Macro Redefinition Protections (config.h)
Added `#ifndef` guards to diagnostic and optimization macros to prevent conflicts with FLINT's header files:
- `DIAGNOSTIC_PUSH`, `DIAGNOSTIC_POP`
- `PUSH_OPTIONS`, `POP_OPTIONS`
- `OPTIMIZE_O2`, `OPTIMIZE_OSIZE`, `OPTIMIZE_UNROLL_LOOPS`

### 3. Macro Guards in pml.h
Added `#ifndef` guards for feature detection macros to avoid conflicts with FLINT:
- `PML_HAVE_MACHINE_VECTORS`
- `PML_HAVE_AVX2`
- `PML_HAVE_AVX512`

## Building

This modified version is designed to be built via CMake as part of the DixonRes project. For standalone compilation, please refer to the original PML documentation.

## License

See the [COPYING](COPYING) and [COPYING_FLINT](COPYING_FLINT) files for licensing information.

## Original Authors

See the [AUTHORS](AUTHORS) file for the original PML contributors.
