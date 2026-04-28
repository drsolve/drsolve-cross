# DRsolve: Dixon Resultant & Polynomial System Solver

A cross-platform C library and command-line tool for computing Dixon resultants and solving polynomial systems over finite fields and the rationals ℚ, based on the FLINT and PML libraries.

## Library & Platform Support

DRsolve is distributed as both a **shared library** (`libdrsolve.so` / `libdrsolve.dylib` / `libdrsolve-1.dll`) and a **static library** (`libdrsolve.a` / `libdrsolve-1.a`), alongside the `drsolve` command-line executable.

| Platform | Shared lib | Static lib | CLI |
|---|---|---|---|
| Linux (x86-64, ARM64) | `libdrsolve.so` | `libdrsolve.a` | `drsolve` |
| macOS (x86-64, Apple Silicon) | `libdrsolve.dylib` | `libdrsolve.a` | `drsolve` |
| Windows (x86-64) | `libdrsolve-1.dll` | `libdrsolve-1.a` | `drsolve.exe` + `drsolve_win_gui.exe` |

---

## Features

- Dixon resultant computation for variable elimination
- Polynomial system solver for n×n systems
- Dixon with triangular ideal reduction
- Finite fields:
  - Prime fields F_p (any size): Implemented with FLINT modular arithmetic, optionally accelerated by PML.
  - Extension fields F_{p^k}: Further optimized for binary fields F_{2^n} with n in {8, 16, 32, 64, 128}.
- Rational field ℚ: Rational reconstruction via multi-prime CRT. Set `field_size = 0` to enable.
- Complexity analysis — estimates Dixon matrix size, Bezout degree bound, and operation count before computing.
- Command-line input or file input. Automatic output to solution files.

---

## Dependencies

- **FLINT** (recommended version: 3.5.0)  
  <https://github.com/flintlib/flint>

Optional:
- **PML** (used automatically if available)  
  <https://github.com/vneiger/pml>

---

## Build

DRsolve uses **CMake** (≥ 3.16) as its primary build system.

### Linux / macOS

```bash
git clone https://github.com/drsolve/drsolve-cross.git && cd drsolve-cross
cmake -B build
cmake --build build -j$(nproc)
ctest --test-dir build          # optional: run tests
sudo cmake --install build      # optional: install to /usr/local
```

FLINT must be installed (e.g. via your package manager) before configuring.  
If FLINT is in a non-standard location, pass `-DFLINT_ROOT=/path/to/flint`.

```bash
# Example: FLINT installed under $HOME/.local
cmake -B build -DFLINT_ROOT=$HOME/.local
cmake --build build -j$(nproc)
```

### Windows

**Option A — GUI installer (recommended)**  
Download from [drsolve/drsolve-win](https://github.com/drsolve/drsolve-win).

**Option B — MSYS2/UCRT64**

```bash
pacman -S mingw-w64-ucrt-x86_64-cmake \
          mingw-w64-ucrt-x86_64-gcc \
          mingw-w64-ucrt-x86_64-flint
cmake -B build -G "MinGW Makefiles"
cmake --build build -j$(nproc)
```

**Option C — cross-compile from Linux/macOS (MinGW-w64)**  
Bundled `third_party/` libraries are used automatically:

```bash
cmake -B build-win \
      -DCMAKE_TOOLCHAIN_FILE="$(pwd)/cmake/toolchain-mingw64.cmake"
cmake --build build-win -j$(nproc)
```

For full build options and advanced configurations, see [BUILDING.md](BUILDING.md).

---

## Usage

### Dixon Resultant (Basic)
```bash
./drsolve "polynomials" "eliminate_vars" field_size
```
Examples:
```bash
./drsolve "x+y+z, x*y+y*z+z*x, x*y*z+1" "x,y" 257
./drsolve "x^2+y^2+z^2-1, x^2+y^2-2*z^2, x+y+z" "x,y" 0
```

---

### Polynomial System Solver (n equations in n variables)
```bash
./drsolve --solve "polynomials" field_size
```
Example:
```bash
./drsolve --solve "x^2 + y^2 + z^2 - 6, x + y + z - 4, x*y*z - x - 1" 257
```

---

### Complexity Analysis

Estimates the difficulty of a Dixon resultant computation **without** performing it.
Reports equation count, variable count, degree sequence, Dixon matrix size
(via Hessenberg recurrence), Bezout degree bound, and complexity in bits.

```bash
./drsolve --comp "polynomials" "eliminate_vars" field_size
./drsolve -c     "polynomials" "eliminate_vars" field_size
```

Examples:
```bash
./drsolve --comp "x^3+y^3+z^3, x^2*y+y^2*z+z^2*x, x+y+z-1" "x,y" 257
```

**Custom omega** — set the matrix-multiplication exponent used in the complexity formula
(default: 2.81):
```bash
./drsolve --comp --omega 2.81 "polynomials" "eliminate_vars" field_size
./drsolve -c -w 2.0            "polynomials" "eliminate_vars" field_size
```

File input uses the same format as Dixon mode:
```bash
./drsolve --comp example.dr          # output: example_comp.dr
```

---

### Extension Fields
```bash
./drsolve "x + y^2 + t, x*y + t*y + 1" "x" 2^8
```
The default settings use `t` as the extension field generator and FLINT's built-in field polynomial.
```bash
./drsolve --solve "x^2 + t*y, x*y + t^2" "2^8: t^8 + t^4 + t^3 + t + 1"
```
(with AES custom polynomial for F_256)

---

### Dixon with Ideal Reduction
```bash
./drsolve --ideal "ideal_generators" "polynomials" "eliminate_vars" field_size
```
Example:
```bash
./drsolve --ideal "a2^3=2*a1+1, a3^3=a1*a2+3" "a1^2+a2^2+a3^2-10, a3^3-a1*a2-3" "a3" 257
```

### Field-equation reduction mode
After each multiplication, reduces x^q → x for every variable.
```bash
./drsolve --field-equation "polynomials" "eliminate_vars" field_size
./drsolve --field-equation -r "[d1,d2,...,dn]" field_size
```
Example:
```bash
./drsolve --field-equation "x0*x2+x1, x0*x1*x2+x2+1, x1*x2+x0+1" "x0,x1" 2
./drsolve --field-equation -r [3]*5 2
```

---

### Silent Mode
```bash
./drsolve --silent [--solve|--comp|-c] <arguments>
```
No console output is produced; the solution/report file is still generated.

---

### Method Selection
```bash
./drsolve --method <num> <args>
./drsolve --step1 <num> --step4 <num> <args>
```

Available determinant methods:
- `0` = Recursive
- `1` = Kronecker+HNF
- `2` = Interpolation
- `3` = Sparse interpolation

`--method` sets both step 1 and step 4 for backward compatibility.

---

## Random Mode

Generate random polynomial systems with specified degrees for testing and benchmarking.

### Basic Usage
```bash
./drsolve --random "[d1,d2,...,dn]" field_size
./drsolve -r       "[d]*n"          field_size
```

- `[d1,d2,...,dn]`: degree list (comma-separated) for n polynomials
- `[d]*n`: all n polynomials have same degree d
- `field_size`: field size (prime or extension); use `0` for Q

### Combine with Compute Flags
```bash
# Random + Dixon elimination
./drsolve -r --solve "[d1,...,dn]" field_size

# Random + complexity analysis
./drsolve -r --comp  "[d]*n" field_size
./drsolve -r -c --omega 2.81 "[4]*5" 257   # custom omega

# Random + Dixon with ideal reduction
./drsolve -r "[d1,d2,d3]" "ideal_generators" field_size
```

### Examples
```bash
# 3 polynomials (deg 3,3,2) in GF(257)
./drsolve --random "[3,3,2]" 257

# 3 polynomials (deg 3,3,2) over Q
./drsolve --random "[3,3,2]" 0

# Solve 3 quadratic system in GF(257)
./drsolve -r --solve "[2]*3" 257

# Complexity analysis of 4 quartic polynomials
./drsolve -r --comp --omega 2.81 "[4]*4" 257

# GF(2^8) with degrees 3 and 2
./drsolve -r "[3,2]" 2^8
```

## File Input Format

### Dixon Mode / Complexity Mode (multiline)
```
Line 1 : field size (prime or p^k; use 0 for Q; generator defaults to 't')
Line 2+: polynomials (comma-separated, may span multiple lines)
Last   : variables to ELIMINATE (comma-separated)
         (#eliminate = #equations - 1)
```
Example:
```bash
./drsolve       example.dr
./drsolve --comp example.dr
```

### Polynomial Solver Mode (multiline)
```
Line 1 : field size
Line 2+: polynomials
         (n equations in n variables)
```

---

## Output

| Mode | Command-line input | File input `example.dr` |
|---|---|---|
| Dixon / Solver | `solution_YYYYMMDD_HHMMSS.dr` | `example_solution.dr` |
| Complexity | `comp_YYYYMMDD_HHMMSS.dr` | `example_comp.dr` |

Each output file contains field information, input polynomials, computation time,
and the resultant, solutions, or complexity report.

### Complexity report contents
- Equation count, variable list, elimination variable list, remaining variables
- Degree sequence of input polynomials
- Bezout bound (product of degrees)
- Dixon matrix size (Hessenberg recurrence)
- Resultant degree estimate
- Complexity in log₂ bits (with the omega value used)

---

## Notes
- All computation modes generate a solution/report file by default
- Extension fields are slower than prime fields due to polynomial arithmetic
- The optional PML library only accelerates well-determined systems over prime fields
- Complexity analysis does not run any polynomial arithmetic; it parses only
- Over Q (`field_size=0`), `--solve`, `--ideal`, and `--field-equation` are not yet supported

---

## Feature Support by Field

| Feature | F_p (p<2^63) | F_p (p>2^63) | F_{p^k} (p<2^63) | Q |
|---|---|---|---|---|
| Dixon resultant | ✅ | ✅ | ✅ | ✅ |
| Complexity analysis (`--comp`) | ✅ | ✅ | ✅ | ✅ |
| Random mode (`-r`) | ✅ | ✅ | ✅ | ✅ |
| Polynomial solver (`--solve`) | ✅ | ❌ | ✅ | ❌ |
| Ideal reduction (`--ideal`) | ✅ | ❌ | ✅ | ❌ |
| Field-equation reduction | ✅ | ❌ | ✅ | ❌ |
| PML acceleration | ✅ | ❌ | ❌ | ❌ |

---

## License
DRsolve is distributed under the GNU General Public License version 2.0 (GPL-2.0-or-later). See the file COPYING.
