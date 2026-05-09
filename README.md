# DRSolve: Dixon Resultant & Polynomial System Solver
A C implementation for computing Dixon resultants and solving polynomial systems over finite fields and the rationals ℚ, based on the FLINT and PML libraries.

Website: <https://drsolve.github.io>

## Features
- Dixon resultant computation for variable elimination
- Polynomial system solver for n×n systems
- Dixon with triangular ideal reduction
- Finite fields:
  - Prime fields F_p (any size): Implemented with FLINT modular arithmetic, optionally accelerated by PML.
  - Extension fields F_{p^k}: Further optimized for binary fields F_{2^n} with n in {8, 16, 32, 64, 128}.
- Rational field ℚ: Rational reconstruction via multi-prime CRT. Set field_size = 0 to enable.
- Complexity analysis — estimates Dixon matrix size, Bezout degree bound, and operation count before computing
- Command line input or file input. Automatic output to solution files

---

## Dependencies
- **FLINT** (recommended version: 3.5.0)  
  <https://github.com/flintlib/flint>

```bash
git clone https://github.com/flintlib/flint.git && cd flint
./bootstrap.sh
./configure 
make
make install
```
  
- **PML** (built in)  
  <https://github.com/vneiger/pml>

---

## Build
```bash
git clone https://github.com/drsolve/drsolve.git && cd drsolve
./configure
make
make check                         # optional
make install                       # optional
```
For more options, run `./configure --help` or `make help`.
We also provide a Windows GUI at [drsolve-win](https://github.com/drsolve/drsolve-win) or [drsolve-cross](https://github.com/drsolve/drsolve-cross).

---

## Usage

### Dixon Resultant (Basic)
```bash
./drsolve "polynomials" "eliminate_vars" field_size
./drsolve -o output.dr "polynomials" "eliminate_vars" field_size
```
Examples:
```bash
./drsolve "x+y+z, x*y+y*z+z*x, x*y*z+1" "x,y" 257
./drsolve -o output.dr "x^2+y^2+z^2-1, x^2+y^2-2*z^2, x+y+z" "x,y" 0
```

---

### Polynomial System Solver (n equations in n variables)
```bash
./drsolve "polynomials" field_size
./drsolve -s "polynomials" field_size
```
`-s` is the short form of solver mode. `--solve` still works for compatibility, and is optional when only `"polynomials" field_size` are provided.
Example:
```bash
./drsolve -s "x^2 + y^2 + z^2 - 6, x + y + z - 4, x*y*z - x - 1" 257
```

## File Input Format

### Auto-detected file input
Without flags, drsolve inspects the first non-space character of line 1:
- digit => solver mode
- otherwise => elimination mode

### Elimination / complexity / ideal mode (multiline)
```
Line 1 : variables to ELIMINATE (comma-separated)
Line 2 : field size (prime or p^k; use 0 for Q; generator defaults to 't')
Line 3+: polynomials (comma-separated, may span multiple lines)
         (#eliminate = #equations - 1)
```
Example:
```bash
# example.dr
x0,x1
257
x0^3+x1^3+x2^3, x0*x1+x1*x2+x2*x1, x1*x2*x0+1
```
Run:
```bash
./drsolve example.dr
./drsolve -f example.dr -o my_result.dr
```

### Polynomial solver mode (multiline)
```
Line 1 : field size
Line 2+: polynomials
         (n equations in n variables)
```
Example:
```bash
# example_solve.dr
257
x^2+y^2+z^2-6, x+y+z-4, x*y*z-x-1
```
Run:
```bash
./drsolve example_solve.dr
./drsolve -s -f example_solve.dr -o my_solutions.dr
```
---

### Verbosity
```bash
./drsolve -v 0 <arguments>
./drsolve -v 1 <arguments>
./drsolve -v 2 <arguments>
./drsolve -v 3 <arguments>
```
`-v 0` prints nothing but still writes the output file. `-v 1` is the default. `-v 2` restores the debug-level console output and timing. `-v 3` additionally dumps small intermediate matrices.

For `--comp`: `-v 0` prints only the final overall complexity, `-v 1` prints the per-step / per-method summary, `-v 2` adds formulas and parameter values, and `-v 3` adds extra detail.

Example:
```bash
./drsolve -v 2 -f in.dr -o out.dr
```
---

### Complexity Analysis
Estimates the difficulty of a Dixon resultant computation **without** performing it.
Reports equation count, variable count, degree sequence, Dixon matrix size
(via Hessenberg recurrence), Bezout degree bound, and complexity in bits.

```bash
./drsolve --comp "polynomials" "eliminate_vars" field_size
./drsolve -c     "polynomials" "eliminate_vars" field_size
./drsolve --comp -f example.dr -o report.dr
```

Examples:
```bash
./drsolve --comp "x^3+y^3+z^3, x^2*y+y^2*z+z^2*x, x+y+z-1" "x,y" 257
```

**Custom omega** — set the matrix-multiplication exponent used in the complexity formula
(default: 2.3):
```bash
./drsolve --comp --omega 2.373 "polynomials" "eliminate_vars" field_size
./drsolve -c -w 2.0            "polynomials" "eliminate_vars" field_size
```

File input uses the same elimination-file format shown above:
```bash
./drsolve --comp example.dr           # default output: out/example_comp.dr
./drsolve --comp -f example.dr -o report.dr
```

---

### Extension Fields
```bash
./drsolve "x + y^2 + t, x*y + t*y + 1" "y" 2^8
```
The default settings use `t` as the extension field generator and FLINT's built-in field polynomial.
```bash
./drsolve -s "x^2 + t*y, x*y + t^2" "2^8: t^8 + t^4 + t^3 + t + 1"
```
(with AES custom polynomial for F_256)

---

### Dixon with Ideal Reduction
```bash
./drsolve --ideal "ideal_generators" "polynomials" "eliminate_vars" field_size
./drsolve --ideal -f input.dr -o output.dr
```
Example:
```bash
./drsolve --ideal "a2^3=2*a1+1, a3^3=a1*a2+3" "a1^2+a2^2+a3^2-10, a3^3-a1*a2-3" "a3" 257
```

### Field-equation reduction mode 
After each multiplication, reduces x^q -> x for every variable.
```bash
./drsolve --field-equation "polynomials" "eliminate_vars" field_size
./drsolve --field-eqution -r "[d1,d2,...,dn]" field_size
```
Example:
```bash
./drsolve --field-eqution "x0*x2+x1, x0*x1*x2+x2+1, x1*x2+x0+1" "x0,x1" 2
./drsolve --field-eqution -r [3]*5 2
```

### Method Selection
```bash
./drsolve --method <num> --threads <num> <args>
./drsolve --dixon <args>
./drsolve --macaulay <args>
./drsolve --subres <args>
```
Available methods: 0. Recursive; 1. Kronecker+HNF; 2. Interpolation; 3. sparse interpolation; 5. fast recursive Dixon construction.

For convenience, 2 equations + 1 elimination variable auto-enable `--subres`, and 3/4 equations with standard Dixon shape auto-enable the fast recursive Dixon construction unless a method is explicitly selected.

**Note:** Only the `Interpolation` method supports multi-threading. The default method HNF or sparse interpolation does not support parallel acceleration.

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
./drsolve -r -s "[d1,...,dn]" field_size

# Random + complexity analysis
./drsolve -r --comp  "[d]*n" field_size
./drsolve -r -c --omega 2.373 "[4]*5" 257   # custom omega

# Random + Dixon with ideal reduction
./drsolve -r "[d1,d2,d3]" "ideal_generators" field_size
```

### Examples
```bash
# 3 polynomials (deg 3,3,2) in GF(257)
./drsolve --random "[3,3,2]" 257

# Solve 3 quadratic system in GF(257)
./drsolve -r -s "[2]*3" 257

# Complexity analysis of 4 quartic polynomials
./drsolve -r --comp --omega 2.373 "[4]*4" 257
```

---

## SageMath Interface

`drsolve_sage_interface.sage` is a thin wrapper that lets you call drsolve directly from a SageMath session.

### Quick start

```python
load("drsolve_sage_interface.sage")
set_dixon_path("./drsolve")   # set once per session

R.<x, y, z> = GF(257)[]
F = [x + y + z - 3, x*y + y*z + z*x - 3, x*y*z - 1]

res  = DixonResultant(F, [x, y])   # Dixon resultant, eliminating x and y
sols = DixonSolve(F)               # enumerate all solutions
info = DixonComplexity(F, [x, y])  # complexity estimate (no arithmetic)

# iterative elimination: output is a plain string, feed into the next call
res1 = DixonResultant([x+y+z, x*y+y*z+z*x+1], [x])
res2 = DixonResultant([res1, y*z-1], [y])
```

### API reference

| Function | Description | Returns |
|---|---|---|
| `DixonResultant(F, elim_vars, ...)` | Dixon resultant, eliminating the specified variables. | String or `None` |
| `DixonSolve(F, ...)` | Solve an n×n system, enumerate all solutions. | List of `{var: val}` dicts; `[]`; or `"infinite"` |
| `DixonComplexity(F, elim_vars, ...)` | Estimate complexity without any polynomial arithmetic. | Dict with `complexity_log2`, `bezout_bound`, `matrix_size`, … |
| `DixonIdeal(F, ideal_gens, elim_vars, ...)` | Dixon resultant with triangular ideal reduction. `ideal_gens`: list of strings like `"a^3=2*b+1"`. | String or `None` |
| `set_dixon_path(p)` / `get_dixon_path()` | Set / get the default path to the `drsolve` binary. | — |
| `ToDixon(...)` / `ToDixonSolver(...)` | Write an input file without running the binary. | — |

`field_size` accepts an integer, `"p^k"` string, `GF(...)` object, or `0` for ℚ; inferred from the polynomial ring if omitted. All main functions also accept `debug=True` and `timeout` (seconds).

---

## Output

| Mode | Command-line input | File input `example.dr` |
|---|---|---|
| Dixon / Solver | `out/solution_YYYYMMDD_HHMMSS.dr` | `out/example_solution.dr` |
| Complexity | `out/comp_YYYYMMDD_HHMMSS.dr` | `out/example_comp.dr` |

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
- All computation modes generate a solution/report file by default under `out/`
- Use `-o output.dr` to override the default output filename
- Extension fields are slower than prime fields due to polynomial arithmetic
- The optional PML library only accelerates well-determined systems over prime fields
- Complexity analysis does not run any polynomial arithmetic; it parses only
- Over Q (field_size=0), `--ideal` and `--field-equation` are not supported

---

## Feature Support by Field

| Feature | F_p (p<2^63) | F_p (p>2^63) | F_{p^k} (p<2^63) | Q |
|---|---|---|---|---|
| Dixon resultant | ✅ | ✅ | ✅ | ✅ |
| Complexity analysis (`--comp`) | ✅ | ✅ | ✅ | ✅ |
| Random mode (`-r`) | ✅ | ✅ | ✅ | ✅ |
| Polynomial solver (`-s` / `--solve`) | ✅ | ✅ | ✅ | ✅ |
| Ideal reduction (`--ideal`) | ✅ | ❌ | ✅ | ❌ |
| Field-equation reduction | ✅ | ❌ | ✅ | ❌ |
| PML acceleration | ✅ | ✅ | ❌ | ✅ |

---

## License
DRSolve is distributed under the GNU General Public License version 2.0 (GPL-2.0-or-later). See the file COPYING.
