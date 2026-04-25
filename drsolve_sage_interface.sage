"""
SageMath interface for drsolve
================================

Usage
-----

# Basic
load("drsolve_sage_interface.sage")
set_dixon_path("./drsolve")

R.<x, y, z> = GF(257)[]
F = [x + y + z - 3, x*y + y*z + z*x - 3, x*y*z - 1]
res  = drsolveultant(F, [x, y])
sols = DixonSolve(F)
info = DixonComplexity(F, [x, y])
print(res, "\n", sols, "\n", info, "\n")

# Iterative elimination
load("dixon_sage_interface.sage")
set_dixon_path("./drsolve")
R.<x, y, z> = GF(17)[]
f1 = x + y + z
f2 = x*y + y*z + z*x + 1
f3 = y*z - 1
f4 = z - 2
res1 = drsolveultant([f1, f2], [x])
res2 = drsolveultant([res1, f3], [y])
res3 = drsolveultant([res2, f4], [z])
print("res1 =", res1, "\nres2 =", res2, "\nres3 =", res3)

# Pass debug=True to any function for verbose diagnostics.

"""

import os
import re
import glob
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
    try:
        from sage.rings.finite_rings.finite_field_base import FiniteField
        if isinstance(field_size, FiniteField):
            p = int(field_size.characteristic())
            k = int(field_size.degree())
            return "%d^%d" % (p, k) if k > 1 else str(p)
    except ImportError:
        pass

    if isinstance(field_size, tuple) and len(field_size) == 2:
        p, k = int(field_size[0]), int(field_size[1])
        return "%d^%d" % (p, k) if k > 1 else str(p)

    if isinstance(field_size, str):
        return field_size

    return str(int(field_size))


def _poly_to_str(f):
    """Sage polynomial or plain string -> string drsolve can parse."""
    return str(f)


def _elim_vars_to_str(elim_vars):
    return ", ".join(str(v) for v in elim_vars)


def _infer_field_size(F, field_size):
    """
    If field_size is None or 0 and F contains at least one Sage polynomial,
    extract the base ring from it and use that.  Otherwise return field_size
    unchanged so the caller can pass an explicit value.
    """
    if field_size is not None and field_size != 0:
        return field_size
    for f in F:
        if not isinstance(f, str):
            try:
                return f.base_ring()
            except AttributeError:
                pass
    # Fallback: 0 means rational / Q mode
    return 0


def _run(cmd, timeout, debug):
    if debug:
        print("[debug] command: %s" % " ".join(cmd))
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        if debug:
            print("[debug] return code: %d" % proc.returncode)
            if proc.stdout.strip():
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


def _locate_output_file(finput, tag, debug):
    """
    Try several candidate paths for the output file drsolve wrote.

    Priority:
    1. The tagged path (file-input mode):
         generate_tagged_filename(finput, tag)
    2. The most recently modified  <tag[1:]>_*.dat  in the CWD
       (CLI mode uses generate_timestamped_filename which writes to CWD)
    3. Any <tag[1:]>*.dat in the same directory as finput

    Returns the path string if found, else None.
    """
    # 1. tagged path next to the input file
    tagged = _find_output_file_tagged(finput, tag)
    if debug:
        print("[debug] looking for tagged output file: %s" % tagged)
    if os.path.isfile(tagged):
        if debug:
            print("[debug] found tagged output file: %s" % tagged)
        return tagged

    # 2. timestamped file in CWD  (solution_*.dr/.dat / comp_*.dr/.dat)
    stem = tag.lstrip("_")
    candidates = []
    for ext in ("dr", "dat"):
        candidates.extend(glob.glob("%s_*.%s" % (stem, ext)))
    candidates = sorted(candidates, key=os.path.getmtime)
    if debug:
        print("[debug] CWD glob '%s_*.{{dr,dat}}': %s" % (stem, candidates))
    if candidates:
        if debug:
            print("[debug] using most recent CWD file: %s" % candidates[-1])
        return candidates[-1]

    # 3. same dir as finput
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

def _check_ring_consistency(F):
    """
    Verify that all *Sage polynomial* elements in F share the same parent.
    Plain strings are skipped.  Raises AssertionError on mismatch.
    """
    ring = None
    for f in F:
        if isinstance(f, str):
            continue
        try:
            r = f.parent()
        except AttributeError:
            continue
        if ring is None:
            ring = r
        else:
            assert r == ring, (
                "Polynomial ring mismatch: %s vs %s" % (ring, r)
            )


def ToDixon(F, elim_vars, field_size=257, finput="/tmp/drsolve_in.dat", debug=False):
    """
    Write a drsolve input file (resultant / complexity mode).

    Each element of F may be a Sage polynomial **or** a plain string
    (e.g. the output of a previous drsolveultant call).

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

def drsolveultant(
    F,
    elim_vars,
    field_size=None,
    dixon_path=None,
    finput="/tmp/drsolve_in.dat",
    debug=False,
    timeout=600,
):
    """
    Compute the Dixon resultant of a polynomial system.

    Parameters
    ----------
    F          : list of Sage polynomials **and/or plain strings**.
                 Strings are accepted so that the output of a previous
                 drsolveultant call can be fed in directly, enabling
                 iterative / cascaded elimination.
    elim_vars  : variables to eliminate (Sage vars or strings)
    field_size : prime, prime power, (p,k), GF(...), or 0 for Q.
                 If None (default), inferred from the first Sage polynomial
                 in F; must be supplied explicitly when F is all strings.
    dixon_path : path to the drsolve executable.
                 If None (default), uses the global default set by
                 set_dixon_path() (initially "./drsolve").
    finput     : temporary input file path
    debug      : print detailed diagnostics
    timeout    : seconds before aborting

    Returns
    -------
    str   raw resultant polynomial string
    None  on failure

    Example – iterative elimination
    --------------------------------
        set_dixon_path("./drsolve")
        R.<x, y, z> = GF(257)[]
        F = [x + y + z, x*y + y*z + z*x, x*y*z + 1]

        res  = drsolveultant(F, [x, y])
        res2 = drsolveultant([res, str(z)], ["z"], field_size=257)
    """
    dixon_path = _resolve_dixon_path(dixon_path)
    field_size = _infer_field_size(F, field_size)

    ToDixon(F, elim_vars, field_size, finput, debug=debug)

    cmd  = [dixon_path, finput]
    proc = _run(cmd, timeout, debug)

    if proc.returncode != 0:
        print("[drsolveultant] drsolve exited with code %d" % proc.returncode)
        if proc.stderr:
            print(proc.stderr)
        return None

    foutput = _locate_output_file(finput, "_solution", debug)
    if foutput is None:
        print("[drsolveultant] could not locate output file")
        return None

    result = _parse_resultant_file(foutput, debug)

    if not debug:
        try:
            os.remove(foutput)
        except OSError:
            pass

    return result


def DixonSolve(
    F,
    field_size=None,
    dixon_path=None,
    finput="/tmp/drsolve_solve_in.dat",
    debug=False,
    timeout=600,
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
    finput     : temporary input file
    debug      : print detailed diagnostics
    timeout    : seconds before aborting

    Returns
    -------
    list of dict  {var_name: value_str}  one dict per solution
    []                                   no solutions
    "infinite"                           positive-dimensional system
    None                                 failure
    """
    dixon_path = _resolve_dixon_path(dixon_path)
    field_size = _infer_field_size(F, field_size)

    ToDixonSolver(F, field_size, finput, debug=debug)

    cmd  = [dixon_path, finput]
    proc = _run(cmd, timeout, debug)

    if proc.returncode != 0:
        print("[DixonSolve] drsolve exited with code %d" % proc.returncode)
        if proc.stderr:
            print(proc.stderr)
        return None

    foutput   = _locate_output_file(finput, "_solution", debug)
    if foutput is None:
        print("[DixonSolve] could not locate output file")
        return None

    solutions = _parse_solutions_file(foutput, debug)

    if not debug:
        try:
            os.remove(foutput)
        except OSError:
            pass

    return solutions


def DixonComplexity(
    F,
    elim_vars,
    field_size=None,
    omega=None,
    dixon_path=None,
    finput="/tmp/drsolve_comp_in.dat",
    debug=False,
    timeout=120,
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
    finput     : temporary input file
    debug      : print detailed diagnostics
    timeout    : seconds before aborting

    Returns
    -------
    dict with keys: complexity_log2, omega, bezout_bound, matrix_size, degrees
    None on failure
    """
    dixon_path = _resolve_dixon_path(dixon_path)
    field_size = _infer_field_size(F, field_size)
    if field_size == 0 or field_size is None:
        field_size = 257   # complexity mode needs a finite field

    ToDixon(F, elim_vars, field_size, finput, debug=debug)

    cmd = [dixon_path, "--comp"]
    if omega is not None:
        cmd += ["--omega", str(omega)]
    cmd.append(finput)

    proc = _run(cmd, timeout, debug)

    if proc.returncode != 0:
        print("[DixonComplexity] drsolve exited with code %d" % proc.returncode)
        return None

    foutput = _locate_output_file(finput, "_comp", debug)
    if foutput is None:
        print("[DixonComplexity] could not locate output file")
        return None

    result = _parse_complexity_file(foutput, debug)

    if not debug:
        try:
            os.remove(foutput)
        except OSError:
            pass

    return result


def DixonIdeal(
    F,
    ideal_gens,
    elim_vars,
    field_size=None,
    dixon_path=None,
    debug=False,
    timeout=600,
):
    """
    Dixon resultant with ideal reduction.

    Parameters
    ----------
    F          : list of Sage polynomials and/or plain strings
    ideal_gens : list of strings like "a^3=2*b+1", or Sage expressions
    elim_vars  : variables to eliminate
    field_size : prime or prime power; inferred from F if None
    dixon_path : path to the drsolve executable.
                 If None (default), uses the global default set by
                 set_dixon_path() (initially "./drsolve").
    debug      : print detailed diagnostics
    timeout    : seconds

    Returns
    -------
    str resultant string, or None on failure
    """
    dixon_path = _resolve_dixon_path(dixon_path)
    field_size = _infer_field_size(F, field_size)
    if field_size == 0 or field_size is None:
        field_size = 257

    ideal_str = ", ".join(str(g) for g in ideal_gens)
    polys_str = ", ".join(_poly_to_str(f) for f in F)
    elim_str  = _elim_vars_to_str(elim_vars)
    field_str = _field_size_str(field_size)

    if debug:
        print("[debug] ideal_str : %s" % ideal_str)
        print("[debug] polys_str : %s" % polys_str)
        print("[debug] elim_str  : %s" % elim_str)
        print("[debug] field_str : %s" % field_str)

    cmd = [dixon_path, "--ideal", ideal_str, polys_str, elim_str, field_str]
    proc = _run(cmd, timeout, debug)

    if proc.returncode != 0:
        print("[DixonIdeal] drsolve exited with code %d" % proc.returncode)
        return None

    candidates = []
    for ext in ("dr", "dat"):
        candidates.extend(glob.glob("solution_*.%s" % ext))
    candidates = sorted(candidates, key=os.path.getmtime)
    if debug:
        print("[debug] solution_*.{dr,dat} in CWD: %s" % candidates)
    if not candidates:
        print("[DixonIdeal] output file not found")
        return None

    foutput = candidates[-1]
    result  = _parse_resultant_file(foutput, debug)

    if not debug:
        try:
            os.remove(foutput)
        except OSError:
            pass

    return result


# Backward-compatible alias matching the README/API name.
DixonResultant = drsolveultant
