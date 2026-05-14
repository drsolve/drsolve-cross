"""
SageMath interface for drsolve
================================

Usage
-----

Functions
---------

DixonRes(F, elim_vars, ...)
    Compute a resultant / elimination polynomial.

DixonSolve(F, ...)
    Solve a polynomial system.

DixonComplexity(F, elim_vars, ...)
    Run complexity analysis.

DixonIdeal(F, ideal_gens, elim_vars, ...)
    Compute a resultant with ideal reduction.
    
set_dixon_path(path)
    Set the default `drsolve` executable path.

get_dixon_path()
    Return the current default executable path.

Common options
--------------

debug=True
    Print wrapper-level diagnostics.

live_output=True
    Stream `drsolve` stdout/stderr directly.

verbosity=0..3
    Pass `-v` to `drsolve`. If omitted, keep the default `-v 1`.

time=True
    Pass `--time`.

threads=<n>
    Pass `--threads <n>`.

method=0..5
    Pass `--method <n>`.

step1=0..4, step4=0..4
    Pass `--step1` / `--step4`.

resultant_method="dixon" | "macaulay" | "subres"
    Select the resultant construction method.

field_equation=True
    Enable field-equation reduction.

fast_ksy=True/False, fast_ksy_col=<idx>
    Control the method-5 KSY precondition.

rational_root_scan="auto" | "off" | "force"
    Solver option for rational root scanning.

rational_only=True
    Keep only exact rational solutions; requires `field_size=0`.

Examples
--------

Basic
-----

load("drsolve_sage_interface.sage")
set_dixon_path("./drsolve")

R.<x, y, z> = GF(257)[]
F = [x + y + z - 3, x*y + y*z + z*x - 3, x*y*z - 1]
res  = DixonRes(F, [x, y])
sols = DixonSolve(F)
info = DixonComplexity(F, [x, y])
print(res, "\n", sols, "\n", info, "\n")

Iterative elimination
---------------------

load("drsolve_sage_interface.sage")
set_dixon_path("./drsolve")
R.<x, y, z> = GF(17)[]
f1 = x + y + z
f2 = x*y + y*z + z*x + 1
f3 = y*z - 1
f4 = z - 2
res1 = DixonRes([f1, f2], [x])
res2 = DixonRes([res1, f3], [y])
res3 = DixonRes([res2, f4], [z])
print("res1 =", res1, "\nres2 =", res2, "\nres3 =", res3)

Method selection
----------------

load("drsolve_sage_interface.sage")
set_dixon_path("./drsolve")

R.<x, y, z> = GF(257)[]
F = [x + y + z - 3, x*y + y*z + z*x - 3, x*y*z - 1]

res_macaulay = DixonRes(F, [x, y], resultant_method="macaulay", verbosity=2)
res_fast = DixonRes(F, [x, y], method=5, fast_ksy=True, fast_ksy_col=0, time=True)
info = DixonComplexity(F, [x, y], omega=2.81, threads=4, verbosity=2)

print("macaulay =", res_macaulay)
print("fast =", res_fast)
print("complexity =", info)

Ideal reduction
---------------

load("drsolve_sage_interface.sage")
set_dixon_path("./drsolve")

R.<x0, x1, x2> = GF(257)[]
F = [x0^2 + x1^2 + x2^2 - 10, x2^3 - x0*x1 - 3]
I = [[x1^3, 2*x0 + 1], [x2^3, x0*x1 + 3]]

res = DixonIdeal(F, I, [x2], verbosity=2, time=True)
print(res)

Extension field
---------------

load("drsolve_sage_interface.sage")
set_dixon_path("./drsolve")

K.<z8> = GF(2^8, modulus=x^8 + x^4 + x^3 + x + 1)
R.<x, y> = PolynomialRing(K, 2)
F = [x^2 + z8*y + 1, x*y + z8^2]

res = DixonRes(F, [y], field_size=K, verbosity=2)
sols = DixonSolve(F, field_size=K, verbosity=2)
print(res)
print(sols)

Rational solving
----------------

load("drsolve_sage_interface.sage")
set_dixon_path("./drsolve")

R.<x, y> = QQ[]
F = [x^2 - 1, y - x]

sols_exact = DixonSolve(F, field_size="0", rational_only=True, verbosity=2)
sols_all = DixonSolve(F, field_size="0", rational_root_scan="force", verbosity=2)
print("exact =", sols_exact)
print("all =", sols_all)

"""

import os
import re
import glob
import tempfile
import subprocess
from itertools import product as _iproduct


# ---------------------------------------------------------------------------
# Global default path — set once with set_dixon_path()
# ---------------------------------------------------------------------------

_default_dixon_path = "drsolve"

def set_dixon_path(path):
    """
    Set the default path to the drsolve binary for all subsequent calls.

    Example
    -------
        set_dixon_path("./drsolve")
        set_dixon_path("/usr/local/bin/drsolve")
    """
    global _default_dixon_path
    _default_dixon_path = path

def get_dixon_path():
    """Return the currently configured default drsolve binary path."""
    return _default_dixon_path

def _resolve_dixon_path(dixon_path):
    """
    Return *dixon_path* if explicitly supplied (not None), otherwise fall back
    to the module-level default set by set_dixon_path().
    """
    if dixon_path is not None:
        return dixon_path
    return _default_dixon_path


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _field_size_str(field_size):
    """
    Convert various field-size representations to the string drsolve expects.

    0 / "0"       -> "0"    (rational / Q mode)
    257           -> "257"
    (2, 8)        -> "2^8"
    GF(257)       -> "257"
    GF(2^8)       -> "2^8"
    "2^8"         -> "2^8"  (passed through)
    """
    def _replace_symbol(expr_str, old_symbol, new_symbol):
        if not old_symbol or old_symbol == new_symbol:
            return expr_str
        pattern = r"(?<![A-Za-z0-9_])%s(?![A-Za-z0-9_])" % re.escape(old_symbol)
        return re.sub(pattern, new_symbol, expr_str)

    def _finite_field_to_str(field):
        p = int(field.characteristic())
        k = int(field.degree())
        if k <= 1:
            return str(p)

        field_str = "%d^%d" % (p, k)
        try:
            modulus = field.modulus()
            gen_name = field.variable_name()
            modulus_var = str(modulus.parent().gen())
            modulus_str = _replace_symbol(str(modulus), modulus_var, gen_name)
            return "%s: %s" % (field_str, modulus_str)
        except Exception:
            return field_str

    try:
        from sage.rings.finite_rings.finite_field_base import FiniteField
        if isinstance(field_size, FiniteField):
            return _finite_field_to_str(field_size)
    except ImportError:
        pass

    if isinstance(field_size, tuple) and len(field_size) == 2:
        p, k = int(field_size[0]), int(field_size[1])
        return "%d^%d" % (p, k) if k > 1 else str(p)

    if isinstance(field_size, tuple) and len(field_size) in (3, 4):
        p, k = int(field_size[0]), int(field_size[1])
        gen_name = field_size[3] if len(field_size) == 4 else None
        modulus_str = str(field_size[2])
        if gen_name is not None:
            gen_name = str(gen_name)
            modulus_str = re.sub(r"\b[a-zA-Z]\b", gen_name, modulus_str, count=1)
        return "%d^%d: %s" % (p, k, modulus_str)

    if isinstance(field_size, str):
        return field_size

    return str(int(field_size))


def _relation_pair(item):
    if isinstance(item, (list, tuple)) and len(item) == 2:
        return item[0], item[1]

    lhs = getattr(item, "lhs", None)
    rhs = getattr(item, "rhs", None)
    if callable(lhs) and callable(rhs):
        try:
            return lhs(), rhs()
        except TypeError:
            pass

    return None


def _iter_polynomial_parts(item):
    if isinstance(item, str):
        return

    pair = _relation_pair(item)
    if pair is not None:
        lhs, rhs = pair
        for part in (lhs, rhs):
            if not isinstance(part, str):
                yield part
        return

    yield item


def _poly_to_str(f):
    """Sage polynomial, equation-like pair, or plain string -> drsolve string."""
    if isinstance(f, str):
        return f.strip()

    pair = _relation_pair(f)
    if pair is not None:
        lhs, rhs = pair
        return str(lhs - rhs)

    return str(f)


def _ideal_gen_to_str(gen):
    """Normalize an ideal generator to the `lhs = rhs` form drsolve expects."""
    if isinstance(gen, str):
        text = gen.strip()
        if "=" in text:
            return text
        return "%s = 0" % text

    pair = _relation_pair(gen)
    if pair is not None:
        lhs, rhs = pair
        return "%s = %s" % (lhs, rhs)

    return "%s = 0" % gen


def _elim_vars_to_str(elim_vars):
    return ", ".join(str(v) for v in elim_vars)

def _infer_field_size(F, field_size, extra_items=None):
    """
    If field_size is None or 0 and F contains at least one Sage polynomial,
    extract the base ring from it and use that.  Otherwise return field_size
    unchanged so the caller can pass an explicit value.
    """
    if field_size is not None and field_size != 0:
        return field_size
    for item in list(F) + list(extra_items or []):
        for part in _iter_polynomial_parts(item):
            try:
                return part.base_ring()
            except AttributeError:
                pass
    # Fallback: 0 means rational / Q mode
    return 0


def _run(cmd, timeout, debug, live_output=False):
    if not isinstance(cmd, (list, tuple)):
        raise ValueError("_run: cmd must be a list, not a string (shell injection risk)")
    if debug:
        print("[debug] command: %s" % " ".join(cmd))
    try:
        if live_output:
            proc = subprocess.run(
                cmd,
                shell=False,
                text=True,
                timeout=timeout,
            )
        else:
            proc = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        if debug:
            print("[debug] return code: %d" % proc.returncode)
            if live_output:
                print("[debug] live_output=True; drsolve stdout/stderr were streamed directly")
            elif proc.stdout.strip():
                print("[debug] stdout:\n%s" % proc.stdout.rstrip())
            if proc.stderr.strip():
                print("[debug] stderr:\n%s" % proc.stderr.rstrip())
        return proc
    except FileNotFoundError:
        raise RuntimeError(
            "drsolve binary not found: '%s'.\n"
            "Set the correct path with set_dixon_path(...) or pass dixon_path=..." % cmd[0]
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError("drsolve timed out after %d s." % timeout)


def _find_output_file_tagged(finput, tag):
    """
    Reproduce generate_tagged_filename() from the C source:
    insert *tag* before the last dot-extension.

      /tmp/drsolve_in.dat  +  "_solution"  ->  /tmp/drsolve_in_solution.dat
      /tmp/drsolve_in      +  "_solution"  ->  /tmp/drsolve_in_solution
    """
    dot = finput.rfind(".")
    if dot != -1:
        return finput[:dot] + tag + finput[dot:]
    return finput + tag


def _find_output_file_tagged_default_out(finput, tag):
    """
    Match the current default drsolve behavior:
    auto-generated outputs are written under ./out/ using the BASENAME only.

      /tmp/drsolve_in.dat  + "_solution" -> out/drsolve_in_solution.dat
    """
    return os.path.join("out", os.path.basename(_find_output_file_tagged(finput, tag)))


def _locate_output_file(finput, tag, debug):
    """
    Try several candidate paths for the output file drsolve wrote.

    Priority:
    1. The tagged path in ./out/ (current file-input default mode):
         out/basename(generate_tagged_filename(finput, tag))
    2. The legacy tagged path next to the input file:
         generate_tagged_filename(finput, tag)
    3. The most recently modified <tag[1:]>_*.dr/.dat in ./out/
       (CLI mode uses generate_timestamped_filename under ./out/)
    4. The legacy most recently modified <tag[1:]>_*.dr/.dat in the CWD
    5. Any <tag[1:]>*.dr/.dat in the same directory as finput

    Returns the path string if found, else None.
    """
    # 1. tagged path in ./out/
    tagged_out = _find_output_file_tagged_default_out(finput, tag)
    if debug:
        print("[debug] looking for out-tagged output file: %s" % tagged_out)
    if os.path.isfile(tagged_out):
        if debug:
            print("[debug] found out-tagged output file: %s" % tagged_out)
        return tagged_out

    # 2. legacy tagged path next to the input file
    tagged = _find_output_file_tagged(finput, tag)
    if debug:
        print("[debug] looking for legacy tagged output file: %s" % tagged)
    if os.path.isfile(tagged):
        if debug:
            print("[debug] found legacy tagged output file: %s" % tagged)
        return tagged

    # 3. timestamped file in ./out/
    stem = tag.lstrip("_")
    candidates = []
    for ext in ("dr", "dat"):
        candidates.extend(glob.glob(os.path.join("out", "%s_*.%s" % (stem, ext))))
    candidates = sorted(candidates, key=os.path.getmtime)
    if debug:
        print("[debug] out glob 'out/%s_*.{dr,dat}': %s" % (stem, candidates))
    if candidates:
        if debug:
            print("[debug] using most recent out file: %s" % candidates[-1])
        return candidates[-1]

    # 4. legacy timestamped file in CWD
    candidates = []
    for ext in ("dr", "dat"):
        candidates.extend(glob.glob("%s_*.%s" % (stem, ext)))
    candidates = sorted(candidates, key=os.path.getmtime)
    if debug:
        print("[debug] legacy CWD glob '%s_*.{dr,dat}': %s" % (stem, candidates))
    if candidates:
        if debug:
            print("[debug] using most recent legacy CWD file: %s" % candidates[-1])
        return candidates[-1]

    # 5. same dir as finput
    d = os.path.dirname(finput) or "."
    candidates = []
    for ext in ("dr", "dat"):
        candidates.extend(glob.glob(os.path.join(d, "%s_*.%s" % (stem, ext))))
    candidates = sorted(candidates, key=os.path.getmtime)
    if debug:
        print("[debug] dir glob '%s/%s_*.{dr,dat}': %s" % (d, stem, candidates))
    if candidates:
        return candidates[-1]

    if debug:
        print("[debug] output file NOT found for tag '%s'" % tag)
    return None


# ---------------------------------------------------------------------------
# File writers  (now accept mixed Sage-polynomial / string lists)
# ---------------------------------------------------------------------------

def _check_ring_consistency(F, extra_items=None):
    """
    Verify that all *Sage polynomial* elements in F share the same parent.
    Plain strings are skipped.  Raises AssertionError on mismatch.
    """
    ring = None
    for item in list(F) + list(extra_items or []):
        for part in _iter_polynomial_parts(item):
            try:
                r = part.parent()
            except AttributeError:
                continue
            if ring is None:
                ring = r
            else:
                assert r == ring, (
                    "Polynomial ring mismatch: %s vs %s" % (ring, r)
                )


def _coerce_str_list(items):
    if items is None:
        return []
    if isinstance(items, str):
        return [part.strip() for part in items.split(",") if part.strip()]
    return [str(item) for item in items]


def _make_temp_path(prefix, suffix):
    fd, path = tempfile.mkstemp(prefix=prefix, suffix=suffix)
    os.close(fd)
    try:
        os.remove(path)
    except OSError:
        pass
    return path


def _cleanup_file(path):
    if path and os.path.exists(path):
        try:
            os.remove(path)
        except OSError:
            pass


def ToDixon(F, elim_vars, field_size=257, finput="/tmp/drsolve_in.dat", debug=False):
    """
    Write a drsolve input file (resultant / complexity mode).

    Each element of F may be a Sage polynomial **or** a plain string
    (e.g. the output of a previous DixonRes call).

    Format:
      Line 1      : comma-separated variables to ELIMINATE
      Line 2      : field size
      Lines 3..n  : one polynomial per line
    """
    _check_ring_consistency(F)

    with open(finput, "w") as fd:
        fd.write(_elim_vars_to_str(elim_vars) + "\n")
        fd.write(_field_size_str(field_size) + "\n")
        fd.write(", ".join(_poly_to_str(f) for f in F) + "\n")

    if debug:
        print("[debug] wrote input file: %s" % finput)
        with open(finput) as fh:
            print("[debug] --- input file content ---")
            print(fh.read().rstrip())
            print("[debug] --- end ---")

    return finput

def ToDixonSolver(F, field_size=257, finput="/tmp/drsolve_solve_in.dat", debug=False):
    """
    Write a drsolve solver input file (no elimination-variable line).

    Each element of F may be a Sage polynomial **or** a plain string.

    Format:
      Line 1     : field size
      Lines 2..n : one polynomial per line
    """
    _check_ring_consistency(F)

    with open(finput, "w") as fd:
        fd.write(_field_size_str(field_size) + "\n")
        fd.write(", ".join(_poly_to_str(f) for f in F) + "\n")

    if debug:
        print("[debug] wrote solver input file: %s" % finput)
        with open(finput) as fh:
            print("[debug] --- input file content ---")
            print(fh.read().rstrip())
            print("[debug] --- end ---")

    return finput


def ToDixonIdeal(
    F,
    ideal_gens,
    elim_vars,
    field_size=257,
    finput="/tmp/drsolve_ideal_in.dat",
    debug=False,
):
    """
    Write a drsolve ideal-reduction input file.

    Format:
      Line 1     : variables to eliminate
      Line 2     : field size
      Lines 3..n : ideal generators and polynomials, one per line
    """
    _check_ring_consistency(F, extra_items=ideal_gens)

    with open(finput, "w") as fd:
        fd.write(_elim_vars_to_str(elim_vars) + "\n")
        fd.write(_field_size_str(field_size) + "\n")
        for gen in ideal_gens or []:
            fd.write(_ideal_gen_to_str(gen) + "\n")
        for poly in F:
            fd.write(_poly_to_str(poly) + "\n")

    if debug:
        print("[debug] wrote ideal input file: %s" % finput)
        with open(finput) as fh:
            print("[debug] --- input file content ---")
            print(fh.read().rstrip())
            print("[debug] --- end ---")

    return finput


def _append_option(cmd, flag, value):
    if value is None:
        return
    cmd.extend([flag, str(value)])


def _build_common_cli_flags(
    verbosity=None,
    time=False,
    threads=None,
    resultant_method=None,
    method=None,
    step1=None,
    step4=None,
    field_equation=False,
    fast_ksy=None,
    fast_ksy_col=None,
    rational_root_scan=None,
):
    cmd = []

    if verbosity is not None:
        _append_option(cmd, "-v", int(verbosity))

    if time:
        cmd.append("--time")
    if threads is not None:
        _append_option(cmd, "--threads", int(threads))

    if resultant_method is not None:
        method_name = str(resultant_method).lower()
        if method_name not in ("dixon", "macaulay", "subres"):
            raise ValueError("resultant_method must be one of: dixon, macaulay, subres")
        cmd.append("--" + method_name)

    if method is not None:
        _append_option(cmd, "--method", int(method))
    if step1 is not None:
        _append_option(cmd, "--step1", int(step1))
    if step4 is not None:
        _append_option(cmd, "--step4", int(step4))

    if field_equation:
        cmd.append("--field-equation")

    if fast_ksy is True:
        cmd.append("--fast-ksy")
    elif fast_ksy is False:
        cmd.append("--no-fast-ksy")

    if fast_ksy_col is not None:
        _append_option(cmd, "--fast-ksy-col", int(fast_ksy_col))

    if rational_root_scan is not None:
        mode = str(rational_root_scan).lower()
        if mode not in ("auto", "off", "force"):
            raise ValueError("rational_root_scan must be one of: auto, off, force")
        cmd.extend(["--rational-root-scan", mode])

    return cmd


def _run_file_mode(
    api_name,
    writer,
    parser,
    dixon_path,
    cmd_prefix,
    finput,
    foutput,
    timeout,
    debug,
    live_output,
):
    created_input = finput is None
    created_output = foutput is None

    if finput is None:
        finput = _make_temp_path("drsolve_", ".in")
    if foutput is None:
        foutput = _make_temp_path("drsolve_", ".out")

    writer(finput)
    cmd = [dixon_path] + list(cmd_prefix) + ["-f", finput, "-o", foutput]
    proc = _run(cmd, timeout, debug, live_output=live_output)

    if proc.returncode != 0:
        print("[%s] drsolve exited with code %d" % (api_name, proc.returncode))
        if proc.stderr:
            print(proc.stderr)
        if created_input and not debug:
            _cleanup_file(finput)
        if created_output and not debug:
            _cleanup_file(foutput)
        return None

    result = parser(foutput, debug)

    if created_input and not debug:
        _cleanup_file(finput)
    if created_output and not debug:
        _cleanup_file(foutput)

    return result


# ---------------------------------------------------------------------------
# Output parsers
# ---------------------------------------------------------------------------

def _parse_resultant_file(foutput, debug):
    """
    Parse a drsolve resultant output file.
    Returns the resultant as a raw string, or None if not found.
    """
    try:
        with open(foutput) as f:
            content = f.read()
    except FileNotFoundError:
        if debug:
            print("[debug] output file does not exist: %s" % foutput)
        return None

    if debug:
        print("[debug] --- raw output file content ---")
        print(content.rstrip())
        print("[debug] --- end ---")

    m = re.search(r"Resultant:\n(.*)", content, re.DOTALL)
    if not m:
        if debug:
            print("[debug] regex 'Resultant:\\n...' did NOT match")
        return None

    raw = m.group(1).strip()
    if debug:
        print("[debug] parsed resultant: %r" % raw[:120])
    return raw


def _parse_solutions_file(foutput, debug):
    """
    Parse a drsolve solver output file.

    Returns
    -------
    list of dict  {var_name: value_str}   one dict per solution
    []                                    no solutions
    "infinite"                            positive-dimensional system
    None                                  parse failure
    """
    try:
        with open(foutput) as f:
            content = f.read()
    except FileNotFoundError:
        if debug:
            print("[debug] solver output file does not exist: %s" % foutput)
        return None

    if debug:
        print("[debug] --- raw solver output file ---")
        print(content.rstrip())
        print("[debug] --- end ---")

    if "has no solutions" in content:
        return []
    if "positive dimension" in content or "positive-dimensional" in content:
        return "infinite"

    solutions = []

    compat = re.search(
        r"=== Compatibility View ===(.*?)=== Solution Complete ===",
        content, re.DOTALL,
    )
    if compat:
        var_vals = {}
        for line in compat.group(1).strip().splitlines():
            m = re.match(r"(\w+)\s*=\s*\{(.*?)\}", line.strip())
            if m:
                var_name = m.group(1)
                raw = m.group(2).strip()
                var_vals[var_name] = [v.strip() for v in raw.split(",")] if raw else []
        if debug:
            print("[debug] compatibility view vars: %s" % list(var_vals.keys()))
        if var_vals:
            keys = list(var_vals.keys())
            for combo in _iproduct(*[var_vals[k] for k in keys]):
                solutions.append(dict(zip(keys, combo)))
            return solutions

    blocks = re.findall(
        r"Solution set \d+:\n(.*?)(?=\nSolution set|\n===|$)",
        content, re.DOTALL,
    )
    for block in blocks:
        sol = {}
        for line in block.strip().splitlines():
            m = re.match(r"\s*(\w+)\s*=\s*(.+)", line)
            if m:
                sol[m.group(1)] = m.group(2).strip()
        if sol:
            solutions.append(sol)

    return solutions if solutions else None


def _parse_complexity_file(foutput, debug):
    """Parse a drsolve complexity output file."""
    try:
        with open(foutput) as f:
            content = f.read()
    except FileNotFoundError:
        if debug:
            print("[debug] complexity output file does not exist: %s" % foutput)
        return None

    if debug:
        print("[debug] --- raw complexity output file ---")
        print(content.rstrip())
        print("[debug] --- end ---")

    result = {}

    m = re.search(r"Complexity \(log2, omega=([\d.]+)\):\s*([\d.]+|inf)", content)
    if m:
        result["omega"]           = float(m.group(1))
        result["complexity_log2"] = float(m.group(2))

    m = re.search(r"Bezout bound.*?:\s*(\d+)", content)
    if m:
        result["bezout_bound"] = int(m.group(1))

    m = re.search(r"Dixon matrix size:\s*(\d+)", content)
    if m:
        result["matrix_size"] = int(m.group(1))

    m = re.search(r"Degree sequence:\s*\[(.*?)\]", content)
    if m:
        result["degrees"] = [int(d) for d in m.group(1).split(",")]

    return result or None


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def DixonRes(
    F,
    elim_vars,
    field_size=None,
    dixon_path=None,
    finput=None,
    foutput=None,
    debug=False,
    live_output=False,
    timeout=None,
    verbosity=None,
    time=False,
    threads=None,
    resultant_method=None,
    method=None,
    step1=None,
    step4=None,
    field_equation=False,
    fast_ksy=None,
    fast_ksy_col=None,
):
    """
    Compute the Dixon resultant of a polynomial system.

    Parameters
    ----------
    F          : list of Sage polynomials **and/or plain strings**.
                 Strings are accepted so that the output of a previous
                 DixonRes call can be fed in directly, enabling
                 iterative / cascaded elimination.
    elim_vars  : variables to eliminate (Sage vars or strings)
    field_size : prime, prime power, (p,k), GF(...), or 0 for Q.
                 If None (default), inferred from the first Sage polynomial
                 in F; must be supplied explicitly when F is all strings.
    dixon_path : path to the drsolve executable.
                 If None (default), uses the global default set by
                 set_dixon_path() (initially "./drsolve").
    finput     : temporary input file path; autogenerated when None
    foutput    : output file path; autogenerated when None
    debug      : print wrapper-level diagnostics; does not change drsolve verbosity
    live_output : stream drsolve stdout/stderr directly during execution
    timeout    : seconds before aborting
    verbosity  : drsolve verbosity level (0..3); omitted means drsolve default (-v 1)
    time       : pass --time
    threads    : pass --threads
    resultant_method : one of "dixon", "macaulay", "subres"
    method     : pass --method <0..5>
    step1      : pass --step1 <0..4>
    step4      : pass --step4 <0..4>
    field_equation : pass --field-equation
    fast_ksy   : True -> --fast-ksy, False -> --no-fast-ksy
    fast_ksy_col : pass --fast-ksy-col

    Returns
    -------
    str   raw resultant polynomial string
    None  on failure

    Example – iterative elimination
    --------------------------------
        set_dixon_path("./drsolve")
        R.<x, y, z> = GF(257)[]
        F = [x + y + z, x*y + y*z + z*x, x*y*z + 1]

        res  = DixonRes(F, [x, y])
        res2 = DixonRes([res, str(z)], ["z"], field_size=257)
    """
    dixon_path = _resolve_dixon_path(dixon_path)
    field_size = _infer_field_size(F, field_size)
    cmd_prefix = _build_common_cli_flags(
        verbosity=verbosity,
        time=time,
        threads=threads,
        resultant_method=resultant_method,
        method=method,
        step1=step1,
        step4=step4,
        field_equation=field_equation,
        fast_ksy=fast_ksy,
        fast_ksy_col=fast_ksy_col,
    )

    return _run_file_mode(
        "DixonRes",
        lambda path: ToDixon(F, elim_vars, field_size, path, debug=debug),
        _parse_resultant_file,
        dixon_path,
        cmd_prefix,
        finput,
        foutput,
        timeout,
        debug,
        live_output,
    )


def DixonSolve(
    F,
    field_size=None,
    dixon_path=None,
    finput=None,
    foutput=None,
    debug=False,
    live_output=False,
    timeout=None,
    verbosity=None,
    time=False,
    threads=None,
    method=None,
    step1=None,
    step4=None,
    rational_root_scan=None,
    rational_only=False,
    solve_verbose=False,
    field_equation=False,
    fast_ksy=None,
    fast_ksy_col=None,
):
    """
    Solve a polynomial system with drsolve.

    Parameters
    ----------
    F          : list of Sage polynomials and/or plain strings
    field_size : prime or prime power; inferred from F[0].base_ring() if None
    dixon_path : path to the drsolve executable.
                 If None (default), uses the global default set by
                 set_dixon_path() (initially "./drsolve").
    finput     : temporary input file; autogenerated when None
    foutput    : output file path; autogenerated when None
    debug      : print wrapper-level diagnostics; does not change drsolve verbosity
    live_output : stream drsolve stdout/stderr directly during execution
    timeout    : seconds before aborting
    verbosity  : drsolve verbosity level (0..3); omitted means drsolve default (-v 1)
    time       : pass --time
    threads    : pass --threads
    method     : pass --method <0..5>
    step1      : pass --step1 <0..4>
    step4      : pass --step4 <0..4>
    rational_root_scan : one of "auto", "off", "force"
    rational_only : use --solve-rational-only (requires field_size=0)
    solve_verbose : pass --solve-verbose
    field_equation : pass --field-equation
    fast_ksy   : True -> --fast-ksy, False -> --no-fast-ksy
    fast_ksy_col : pass --fast-ksy-col

    Returns
    -------
    list of dict  {var_name: value_str}  one dict per solution
    []                                   no solutions
    "infinite"                           positive-dimensional system
    None                                 failure
    """
    dixon_path = _resolve_dixon_path(dixon_path)
    field_size = _infer_field_size(F, field_size)

    cmd_prefix = _build_common_cli_flags(
        verbosity=verbosity,
        time=time,
        threads=threads,
        method=method,
        step1=step1,
        step4=step4,
        field_equation=field_equation,
        fast_ksy=fast_ksy,
        fast_ksy_col=fast_ksy_col,
        rational_root_scan=rational_root_scan,
    )
    if rational_only:
        cmd_prefix.append("--solve-rational-only")
    elif solve_verbose:
        cmd_prefix.append("--solve-verbose")
    else:
        cmd_prefix.append("--solve")

    return _run_file_mode(
        "DixonSolve",
        lambda path: ToDixonSolver(F, field_size, path, debug=debug),
        _parse_solutions_file,
        dixon_path,
        cmd_prefix,
        finput,
        foutput,
        timeout,
        debug,
        live_output,
    )


def DixonComplexity(
    F,
    elim_vars,
    field_size=None,
    omega=None,
    dixon_path=None,
    finput=None,
    foutput=None,
    debug=False,
    live_output=False,
    timeout=60,
    verbosity=None,
    time=False,
    threads=None,
    resultant_method=None,
    method=None,
    step1=None,
    step4=None,
    field_equation=False,
    fast_ksy=None,
    fast_ksy_col=None,
):
    """
    Run drsolve complexity analysis.

    Parameters
    ----------
    F          : list of Sage polynomials and/or plain strings
    elim_vars  : variables to eliminate
    field_size : prime or prime power; inferred from F if None
    omega      : matrix-multiplication exponent (default: drsolve built-in)
    dixon_path : path to the drsolve executable.
                 If None (default), uses the global default set by
                 set_dixon_path() (initially "./drsolve").
    finput     : temporary input file; autogenerated when None
    foutput    : output file path; autogenerated when None
    debug      : print wrapper-level diagnostics; does not change drsolve verbosity
    live_output : stream drsolve stdout/stderr directly during execution
    timeout    : seconds before aborting
    verbosity  : drsolve verbosity level (0..3); omitted means drsolve default (-v 1)
    time       : pass --time
    threads    : pass --threads
    resultant_method : one of "dixon", "macaulay", "subres"
    method     : pass --method <0..5>
    step1      : pass --step1 <0..4>
    step4      : pass --step4 <0..4>
    field_equation : pass --field-equation
    fast_ksy   : True -> --fast-ksy, False -> --no-fast-ksy
    fast_ksy_col : pass --fast-ksy-col

    Returns
    -------
    dict with keys: complexity_log2, omega, bezout_bound, matrix_size, degrees
    None on failure
    """
    dixon_path = _resolve_dixon_path(dixon_path)
    field_size = _infer_field_size(F, field_size)
    if field_size == 0 or field_size is None:
        field_size = 257   # complexity mode needs a finite field

    cmd_prefix = _build_common_cli_flags(
        verbosity=verbosity,
        time=time,
        threads=threads,
        resultant_method=resultant_method,
        method=method,
        step1=step1,
        step4=step4,
        field_equation=field_equation,
        fast_ksy=fast_ksy,
        fast_ksy_col=fast_ksy_col,
    )
    cmd_prefix.append("--comp")
    if omega is not None:
        cmd_prefix.extend(["--omega", str(omega)])

    return _run_file_mode(
        "DixonComplexity",
        lambda path: ToDixon(F, elim_vars, field_size, path, debug=debug),
        _parse_complexity_file,
        dixon_path,
        cmd_prefix,
        finput,
        foutput,
        timeout,
        debug,
        live_output,
    )


def DixonIdeal(
    F,
    ideal_gens,
    elim_vars,
    field_size=None,
    dixon_path=None,
    finput=None,
    foutput=None,
    debug=False,
    live_output=False,
    timeout=None,
    verbosity=None,
    time=False,
    threads=None,
    resultant_method=None,
    method=None,
    step1=None,
    step4=None,
    field_equation=False,
    fast_ksy=None,
    fast_ksy_col=None,
):
    """
    Dixon resultant with ideal reduction.

    Parameters
    ----------
    F          : list of Sage polynomials and/or plain strings
    ideal_gens : list whose elements may be:
                 - strings like "a^3=2*b+1"
                 - Sage polynomials / expressions, interpreted as `gen = 0`
                 - pairs like [lhs, rhs] or (lhs, rhs)
                 - Sage equality/relational objects when available
    elim_vars  : variables to eliminate
    field_size : prime or prime power; inferred from F if None
    dixon_path : path to the drsolve executable.
                 If None (default), uses the global default set by
                 set_dixon_path() (initially "./drsolve").
    finput     : temporary input file; autogenerated when None
    foutput    : output file path; autogenerated when None
    debug      : print wrapper-level diagnostics; does not change drsolve verbosity
    live_output : stream drsolve stdout/stderr directly during execution
    timeout    : seconds
    verbosity  : drsolve verbosity level (0..3); omitted means drsolve default (-v 1)
    time       : pass --time
    threads    : pass --threads
    resultant_method : one of "dixon", "macaulay", "subres"
    method     : pass --method <0..5>
    step1      : pass --step1 <0..4>
    step4      : pass --step4 <0..4>
    field_equation : pass --field-equation
    fast_ksy   : True -> --fast-ksy, False -> --no-fast-ksy
    fast_ksy_col : pass --fast-ksy-col

    Returns
    -------
    str resultant string, or None on failure
    """
    dixon_path = _resolve_dixon_path(dixon_path)
    field_size = _infer_field_size(F, field_size, extra_items=ideal_gens)
    if field_size == 0 or field_size is None:
        raise ValueError("DixonIdeal requires a finite field; drsolve does not support --ideal over Q.")

    cmd_prefix = _build_common_cli_flags(
        verbosity=verbosity,
        time=time,
        threads=threads,
        resultant_method=resultant_method,
        method=method,
        step1=step1,
        step4=step4,
        field_equation=field_equation,
        fast_ksy=fast_ksy,
        fast_ksy_col=fast_ksy_col,
    )
    cmd_prefix.append("--ideal")

    return _run_file_mode(
        "DixonIdeal",
        lambda path: ToDixonIdeal(F, ideal_gens, elim_vars, field_size, path, debug=debug),
        _parse_resultant_file,
        dixon_path,
        cmd_prefix,
        finput,
        foutput,
        timeout,
        debug,
        live_output,
    )


# Backward-compatible alias matching the README/API name.
DixonResultant = DixonRes
