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
### BASIC USAGE

#### Elimination / resultant mode
```bash
./drsolve "polynomials" "eliminate_vars" field_size
```
Example:
```bash
./drsolve "x+y+z, x*y+y*z+z*x, x*y*z+1" "x,y" 257
```
- Default output file: `out/solution_YYYYMMDD_HHMMSS.dr`

#### Polynomial system solver
```bash
./drsolve "polynomials" field_size
```
Example:
```bash
./drsolve "x^2+y^2+z^2-6, x+y+z-4, x*y*z-x-1" 0
```
- Writes all solutions to `out/solution_YYYYMMDD_HHMMSS.dr`

### FILE FORMAT
#### File input/output
```bash
./drsolve input_file
./drsolve -f input_file -o output.dr
```
#### Dixon resultant elimination (multiline)
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
- If line 1 lists `n` variables for `n` equations, compatibility mode uses the first `n-1` variables

#### Polynomial solver mode (multiline)
```
Line 1 : field size
Line 2+: polynomials (comma-separated, may span multiple lines)
         (n equations in n variables)
```
Example:
```bash
# example_solve.dr
0
x^2+y^2+z^2-6, x+y+z-4, x*y*z-x-1
```
Run:
```bash
./drsolve example_solve.dr
./drsolve -f example_solve.dr -o my_solutions.dr
```

### OTHER MODES

#### Extension fields
```bash
./drsolve "x + y^2 + t, x*y + t*y + 1" "y" 2^8
```
The default settings use `t` as the extension field generator and FLINT's built-in field polynomial.
```bash
./drsolve -s "x^2 + t*y, x*y + t^2" "2^8: t^8 + t^4 + t^3 + t + 1"
```
- Example: AES polynomial for `GF(2^8)`
- In `Q` and prime fields, `t` is treated as an ordinary variable; only extension fields reserve it as the generator

#### Complexity analysis
Estimates the difficulty of a Dixon resultant computation without performing it.
Reports equation count, variable count, degree sequence, Dixon matrix size,
Bezout degree bound, and complexity in bits.

```bash
./drsolve -c "polynomials" "eliminate_vars" field_size
./drsolve -c -f input.dr
```
- Prints complexity information
- Default output file: `out/comp_YYYYMMDD_HHMMSS.dr`
- Add `--omega <value>` or `-w <value>` to set the matrix-multiplication exponent

Example:
```bash
./drsolve -c "x^3+y^3+z^3, x^2*y+y^2*z+z^2*x, x+y+z-1" "x,y" 257
```

#### Random mode
Generate random polynomial systems for testing and benchmarking.

```bash
./drsolve -r       "[d1,d2,...,dn]" field_size
./drsolve -r       "[d]*n"          field_size
./drsolve -r -n 4 --density 0.5 "[d]*3" field_size
./drsolve -r -s    "[d1,...,dn]"    field_size
./drsolve -r -c    "[d]*n"         field_size
```
- Add `-n <num_vars>` to set the variable count
- Add `--density <ratio>` with `0 <= ratio <= 1`
- Add `--seed <num>` for reproducible random systems
- Mixed degree specifications such as `"[2]*5+[3]*6"` are supported

Examples:
```bash
./drsolve -r "[3,3,2]" 257
./drsolve -r "[3]*3" 0
./drsolve -r -n 4 --density 0.5 "[3]*3" 257
./drsolve -r --seed 12345 "[3]*3" 257
./drsolve -r "[2]*4+[3]*2" 257
./drsolve -r -s "[2]*3" 257
./drsolve -r --comp --omega 2.373 "[4]*4" 257
```

#### Dixon with ideal reduction
```bash
./drsolve --ideal "ideal_generators" "polynomials" "eliminate_vars" field_size
./drsolve --ideal -f input.dr -o output.dr
```
- `ideal_generators` is a comma-separated list of relations with `=`
- In file mode, lines after the first two lines containing `=` are treated as ideal generators

Example:
```bash
./drsolve --ideal "a2^3=2*a1+1, a3^3=a1*a2+3" "a1^2+a2^2+a3^2-10, a3^3-a1*a2-3" "a3" 257
```

#### Field-equation reduction mode
After each multiplication, reduces `x^q -> x` for every variable.

```bash
./drsolve --field-equation "polynomials" "eliminate_vars" field_size
./drsolve --field-equation -r "[d1,d2,...,dn]" field_size
```

Example:
```bash
./drsolve --field-equation "x0*x2+x1, x0*x1*x2+x2+1, x1*x2+x0+1" "x0,x1" 2
./drsolve --field-equation -r "[3]*5" 2
```

### OPTIONS

#### Method selection
```bash
./drsolve --method <num> <args>
./drsolve --step1 <num> --step4 <num> <args>
```
- Available methods: `0.Recursive`, `1.Kronecker+HNF`, `2.Interpolation`, `3.Sparse interpolation`, `4.Bareiss`, `5.Recursive Dixon construction`
- `--method` sets both Step 1 and Step 4 for backward compatibility

#### Resultant construction
```bash
./drsolve --dixon <args>
./drsolve --macaulay <args>
./drsolve --subres <args>
```
- `--dixon`, `--macaulay`, and `--subres` are direct method selectors
- `--subres` is for exactly 2 polynomials and 1 elimination variable

#### Verbosity
```bash
./drsolve -v 0 <arguments>
./drsolve -v 1 <arguments>
./drsolve -v 2 <arguments>
./drsolve -v 3 <arguments>
```
`-v 0` prints nothing but still writes the output file. `-v 1` is the default. `-v 2` restores the debug-level console output and timing. `-v 3` additionally dumps small intermediate matrices.

Example:
```bash
./drsolve -v 2 -f in.dr -o out.dr
```

#### Diagnostics
```bash
./drsolve --time <args>
./drsolve -v 2 <args>
```
- `--time` prints per-step timing
- Compatibility flags `--silent`, `--debug`, `--solve-verbose`, and `--solve` are still accepted


#### Parallelism
```bash
./drsolve --threads <num> <args>
```
- Sets the number of threads for parallel computation

### EXAMPLES

```bash
./drsolve "x+y+z, x*y+y*z+z*x, x*y*z+1" "x,y" 257
./drsolve "x^2+y^2+z^2-1, x^2+y^2-2*z^2, x+y+z" "x,y" 0
./drsolve "x^2+y^2+z^2-6, x+y+z-4, x*y*z-x-1" 257
./drsolve --comp "x^2+y^2+1, x*y+z, x+y+z^2" "x,y" 257
./drsolve --random "[3,3,2]" 257
./drsolve --ideal "a2^3=2*a1+1, a3^3=a1*a2+3" "a1^2+a2^2+a3^2-10, a3^3-a1*a2-3" "a3" 257
./drsolve --field-equation "x0*x2+x1, x0*x1*x2+x2+1, x1*x2+x0+1" "x0,x1" 2
./drsolve -v 0 "x+y^2+t, x*y+t*y+1" "y" 2^8
./drsolve "x^2 + t*y, x*y + t^2" "2^8: t^8 + t^4 + t^3 + t + 1"
./drsolve example.dr
./drsolve -v 2 -f in.dr -o out.dr
./drsolve example_solve.dr
```

---

## SageMath Interface

`drsolve_sage_interface.sage` lets you call DRsolve directly from SageMath with Sage polynomial objects.

- Load the interface with `load("drsolve_sage_interface.sage")`, then set the binary path once with `set_dixon_path("./drsolve")`.
- Main entry points:
  - `DixonRes(F, elim_vars, ...)` / `DixonResultant(...)`
  - `DixonSolve(F, ...)`
  - `DixonComplexity(F, elim_vars, ...)`
  - `DixonIdeal(F, ideal_gens, elim_vars, ...)`
- Common options include `field_size`, `verbosity`, `time`, `threads`, `debug`, `live_output`, and `timeout`.
- `field_size` may be an integer prime, a string such as `"2^8"` or `"2^8: t^8+t^4+t^3+t+1"`, a Sage `GF(...)` object, or `0` for ℚ. If omitted, it is inferred from the Sage polynomial ring when possible.
- Resultants are returned as strings, so iterative elimination works naturally by feeding one `DixonRes(...)` output into the next call.
- For a fuller Sage reference with examples and options, see `index.txt` or the top docstring in `drsolve_sage_interface.sage`.

---

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
