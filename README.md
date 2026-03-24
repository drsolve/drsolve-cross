# DixonRes: Dixon Resultant & Polynomial System Solver

A cross-platform C library and command-line tool for computing Dixon resultants and solving polynomial systems over finite fields and the rationals ℚ, based on the FLINT and PML libraries.

Website: <https://dixonres.github.io>

## Library & Platform Support

DixonRes is distributed as both a **shared library** (`libdixon.so` / `libdixon.dylib` / `libdixon-1.dll`) and a **static library** (`libdixon.a` / `libdixon-1.a`), alongside the `dixon` command-line executable.

| Platform | Shared lib | Static lib | CLI |
|---|---|---|---|
| Linux (x86-64, ARM64) | `libdixon.so` | `libdixon.a` | `dixon` |
| macOS (x86-64, Apple Silicon) | `libdixon.dylib` | `libdixon.a` | `dixon` |
| Windows (x86-64) | `libdixon-1.dll` | `libdixon-1.a` | `dixon.exe` + `dixon_win_gui.exe` |

Pre-built static binaries for Linux and macOS are available on the [Releases](https://github.com/DixonRes/DixonRes/releases) page.  
A Windows GUI installer is available at [DixonRes/DixonRes-Windows](https://github.com/DixonRes/DixonRes-Windows).

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

- **FLINT** (recommended version: 3.4.0)  
  <https://github.com/flintlib/flint>

Optional:
- **PML** (used automatically if available)  
  <https://github.com/vneiger/pml>

---

## Build

DixonRes uses **CMake** (≥ 3.16) as its primary build system.

### Linux / macOS

```bash
git clone https://github.com/DixonRes/DixonRes.git && cd DixonRes
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
Download from [DixonRes/DixonRes-Windows](https://github.com/DixonRes/DixonRes-Windows).

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
./dixon "polynomials" "eliminate_vars" field_size
```
Examples:
```bash
./dixon "x+y+z, x*y+y*z+z*x, x*y*z+1" "x,y" 257
./dixon "x^2+y^2+z^2-1, x^2+y^2-2*z^2, x+y+z" "x,y" 0
```

---

### Polynomial System Solver (n equations in n variables)
```bash
./dixon --solve "polynomials" field_size
```
Example:
```bash
./dixon --solve "x^2 + y^2 + z^2 - 6, x + y + z - 4, x*y*z - x - 1" 257
```

---

### Complexity Analysis

Estimates the difficulty of a Dixon resultant computation **without** performing it.
Reports equation count, variable count, degree sequence, Dixon matrix size
(via Hessenberg recurrence), Bezout degree bound, and complexity in bits.

```bash
./dixon --comp "polynomials" "eliminate_vars" field_size
./dixon -c     "polynomials" "eliminate_vars" field_size
```

Examples:
```bash
./dixon --comp "x^3+y^3+z^3, x^2*y+y^2*z+z^2*x, x+y+z-1" "x,y" 257
```

**Custom omega** — set the matrix-multiplication exponent used in the complexity formula
(default: 2.3):
```bash
./dixon --comp --omega 2.373 "polynomials" "eliminate_vars" field_size
./dixon -c -w 2.0            "polynomials" "eliminate_vars" field_size
```

File input uses the same format as Dixon mode:
```bash
./dixon --comp example.dat          # output: example_comp.dat
```

---

### Extension Fields
```bash
./dixon "x + y^2 + t, x*y + t*y + 1" "x" 2^8
```
The default settings use `t` as the extension field generator and FLINT's built-in field polynomial.
```bash
./dixon --solve "x^2 + t*y, x*y + t^2" "2^8: t^8 + t^4 + t^3 + t + 1"
```
(with AES custom polynomial for F_256)

---

### Dixon with Ideal Reduction
```bash
./dixon --ideal "ideal_generators" "polynomials" "eliminate_vars" field_size
```
Example:
```bash
./dixon --ideal "a2^3=2*a1+1, a3^3=a1*a2+3" "a1^2+a2^2+a3^2-10, a3^3-a1*a2-3" "a3" 257
```

### Field-equation reduction mode
After each multiplication, reduces x^q → x for every variable.
```bash
./dixon --field-equation "polynomials" "eliminate_vars" field_size
./dixon --field-equation -r "[d1,d2,...,dn]" field_size
```
Example:
```bash
./dixon --field-equation "x0*x2+x1, x0*x1*x2+x2+1, x1*x2+x0+1" "x0,x1" 2
./dixon --field-equation -r [3]*5 2
```

---

### Silent Mode
```bash
./dixon --silent [--solve|--comp|-c] <arguments>
```
No console output is produced; the solution/report file is still generated.

---

## Random Mode

Generate random polynomial systems with specified degrees for testing and benchmarking.

### Basic Usage
```bash
./dixon --random "[d1,d2,...,dn]" field_size
./dixon -r       "[d]*n"          field_size
```

- `[d1,d2,...,dn]`: degree list (comma-separated) for n polynomials
- `[d]*n`: all n polynomials have same degree d
- `field_size`: field size (prime or extension); use `0` for Q

### Combine with Compute Flags
```bash
# Random + Dixon elimination
./dixon -r --solve "[d1,...,dn]" field_size

# Random + complexity analysis
./dixon -r --comp  "[d]*n" field_size
./dixon -r -c --omega 2.373 "[4]*5" 257   # custom omega

# Random + Dixon with ideal reduction
./dixon -r "[d1,d2,d3]" "ideal_generators" field_size
```

### Examples
```bash
# 3 polynomials (deg 3,3,2) in GF(257)
./dixon --random "[3,3,2]" 257

# 3 polynomials (deg 3,3,2) over Q
./dixon --random "[3,3,2]" 0

# Solve 3 quadratic system in GF(257)
./dixon -r --solve "[2]*3" 257

# Complexity analysis of 4 quartic polynomials
./dixon -r --comp --omega 2.373 "[4]*4" 257

# GF(2^8) with degrees 3 and 2
./dixon -r "[3,2]" 2^8
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
./dixon       example.dat
./dixon --comp example.dat
```

### Polynomial Solver Mode (multiline)
```
Line 1 : field size
Line 2+: polynomials
         (n equations in n variables)
```

---

## Output

| Mode | Command-line input | File input `example.dat` |
|---|---|---|
| Dixon / Solver | `solution_YYYYMMDD_HHMMSS.dat` | `example_solution.dat` |
| Complexity | `comp_YYYYMMDD_HHMMSS.dat` | `example_comp.dat` |

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
DixonRes is distributed under the GNU General Public License version 2.0 (GPL-2.0-or-later). See the file COPYING.