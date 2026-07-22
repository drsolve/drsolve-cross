#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include <flint/flint.h>

#include "dixon_complexity.h"
#include "drsolve_cli.h"
#include "fmpq_acb_roots.h"

#define PROGRAM_VERSION "0.4.3"
#define DEFAULT_OUTPUT_DIR "out"

/* =========================================================================
 * Print usage
 * ========================================================================= */

void drsolve_cli_print_version(void)
{
    printf("===============================================\n");
    printf("DRSolve v%s\n", PROGRAM_VERSION);
    printf("FLINT version: %s (Recommended: 3.6.0)\n", FLINT_VERSION);
#ifdef HAVE_PML
    printf("PML support: ENABLED\n");
#else
    printf("PML support: DISABLED\n");
#endif
    printf("===============================================\n");
}

static void print_short_usage(const char *prog_name)
{
    printf("USAGE:\n");
    printf("  %s \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("  %s \"polynomials\" field_size\n", prog_name);
    printf("  %s input_file -o output_file\n", prog_name);
    printf("FILE FORMAT:\n");
    printf("  Dixon resultant elimination:\n");
    printf("    Line 1 : variables TO ELIMINATE (comma-separated)\n");
    printf("    Line 2 : field size (prime or p^k, 0 means rational)\n");
    printf("    Line 3+: polynomials (comma-separated, may span multiple lines)\n");
    printf("  Solver mode:\n");
    printf("    Line 1 : field size (prime or p^k, 0 means rational)\n");
    printf("    Line 2+: polynomials (comma-separated, may span multiple lines)\n");
    printf("OPTIONS:\n");
    printf("  -r \"[d1,d2,...,dn]\" random polynomial generation\n");
    printf("  -s  solving mode (auto-enables when no vars given)\n");
    printf("  -c, --comp, --complexity  complexity analysis mode\n");
    printf("EXAMPLES:\n");
    printf("  Dixon resultant elimination:\n");
    printf("    %s \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    %s \"x^2+y^2+z^2-1, x^2+y^2-2*z^2, x+y+z\" \"x,y\" 0\n", prog_name);
    printf("  Polynomial system solving:\n");
    printf("    %s \"x^3+y^2+z-8, x+y+z-6, x*y*z-6\" 0\n", prog_name);
    printf("    %s \"x^2 + t*y, x*y + t^2\" \"2^8: t^8 + t^4 + t^3 + t + 1\"\n", prog_name);
    printf("  Random input:\n");
    printf("    %s -r \"[3]*4\" 257\n", prog_name);
    printf("    %s -r -s \"[2,3,4]\" 2^8 --seed 1234\n", prog_name);
    printf("  Complexity analysis:\n");
    printf("    %s -c \"x^2+y^2+1, x*y+z, x+y+z^2\" \"x,y\" 257\n", prog_name);
    printf("    %s -c -r \"[10]*10\" 257\n", prog_name);
    printf("  File input:\n");
    printf("    %s example.dr\n", prog_name);
    printf("    %s example_solve.dr -o solution.dr\n", prog_name);
    printf("OTHER OPTIONS:\n");
    printf("  --method <n>      Determinant method selection (0:Expansion, 1:HNF, 2:Interpolation, 3:Sparse, 4:Bareiss, 5:Fdixon)\n");
    printf("  --step1, --step4  Override method <n> for specific algorithm steps\n");
    printf("  --cache <num>     Determinant memoization cache entry limit (default: 1024)\n");
    printf("  --threads <num>   Set number of threads for parallel computation\n");
    printf("  --max-primes <n>  Maximum primes for rational reconstruction (Q default: 64; large-prime fallback: 256)\n");
    printf("  --dixon           Use Dixon resultant (default)\n");
    printf("  --macaulay        Use Macaulay resultant\n");
    printf("  --subres          Use Subresultant (2 polys)\n");
    printf("  --field-equation  After each multiplication, reduces x^q -> x for every variable\n");
    printf("  --ideal <args>    After each multiplication, reduces using the given substitution\n");
    printf("  --complex         Output complex solutions (2x2 solver or complex roots over Q)\n");
    printf("  --test <n>        Run built-in tests (1: Dixon matrix size, 2: Bezout bound, 3: solver correctness, 4: performance)\n");
    printf("  --time            Print per-step timing information\n");
    printf("  -v, --verbose <n> Verbosity level (0:silent, 1:default, 2:detailed, 3:debug)\n");
    printf("  -h, --help        Show full detailed help information\n");
    printf("  -V, --version     Print version and build information\n");
}

void drsolve_cli_print_usage(const char *prog_name)
{
    printf("USAGE:\n");
    printf("  Elimination / resultant mode:\n");
    printf("    %s \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    %s -o output.dr \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    Example: %s \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    Example: %s \"x^2+y^2+z^2-1, x^2+y^2-2*z^2, x+y+z\" \"x,y\" 0\n", prog_name);
    printf("    Example: %s --resultant-only \"x^2+y^2-1, x-y\" \"x\" 0\n", prog_name);
    printf("    Example: %s --approx-roots --root-precision 256 \"x+y, x-y^2+2\" \"x\" 0\n", prog_name);
    printf("    Example: %s -o out/elimination.dr \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    -> Default output file: %s/solution_YYYYMMDD_HHMMSS.dr\n", DEFAULT_OUTPUT_DIR);
    printf("\n");

    printf("  Polynomial system solver:\n");
    printf("    %s \"polynomials\" field_size\n", prog_name);
    printf("    %s -s \"polynomials\" field_size\n", prog_name);
    printf("    %s --solve-rational-only \"polynomials\" 0\n", prog_name);
    printf("    %s -v 2 -s \"polynomials\" field_size\n", prog_name);
    printf("    %s -v 3 \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    Example: %s \"x^2-1, y-x\" 257\n", prog_name);
    printf("    Example: %s -s \"x^2+y^2-1, x-y\" 0\n", prog_name);
    printf("    Example: %s --solve-rational-only \"x^2-1, y-x\" 0\n", prog_name);
    printf("    Example: %s -v 2 -s \"x^2-1, y-x\" 257\n", prog_name);
    printf("    Example: %s -v 3 \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    Example: %s \"x^2 + t*y, x*y + t^2\" \"2^8: t^8 + t^4 + t^3 + t + 1\"\n", prog_name);
    printf("    Example: %s -s --specialize-vars \"z=0,w=1\" \"x+z, y+w\" 257\n", prog_name);
    printf("    -> Writes all solutions to %s/solution_YYYYMMDD_HHMMSS.dr\n", DEFAULT_OUTPUT_DIR);
    printf("    -> `-s` / `--solve` is optional here; `--solve-rational-only` keeps only exact rational solutions\n");
    printf("    -> `--specialize-vars a=0,b=1` fixes variables before solving; in extension fields integers use the generator expansion (F_2^k: 2=t, 3=1+t, 4=t^2)\n");
    printf("    -> `-v 2` matches the old debug / verbose solver output\n");
    printf("    -> `-v 3` also dumps small Step 1/2/3 matrices (<= 10 x 10)\n");
    printf("    -> In extension fields, 't' is the field generator; in Q and prime fields it is an ordinary variable\n");
    printf("\n");

    printf("  File input:\n");
    printf("    %s input_file\n", prog_name);
    printf("    %s -f input_file\n", prog_name);
    printf("    %s -s input_file\n", prog_name);
    printf("    %s -s -f input_file -o output.dr\n", prog_name);
    printf("    Example: %s example.dr\n", prog_name);
    printf("    Example: %s -f example.dr\n", prog_name);
    printf("    Example: %s -s example_solve.dr\n", prog_name);
    printf("    Example: %s -s -f example_solve.dr -o out/example_solve_copy.dr\n", prog_name);
    printf("    -> Without flags, auto-detects solver mode when line 1 starts with a digit; otherwise uses elimination mode\n");
    printf("    -> If elimination vars count equals equation count n, auto-adjusts to eliminate the first n-1 variables\n");
    printf("    -> Input file may be given directly or with `-f`; output file must be given with `-o` when overriding the default\n");
    printf("    -> Default file-mode outputs are input_file_solution.dr or input_file_comp.dr\n");
    printf("\n");

    printf("FILE FORMAT (auto-detected for input_file):\n");
    printf("  Solver mode (line 1 starts with a digit):\n");
    printf("    Line 1 : field size\n");
    printf("    Line 2+: polynomials (one per line or comma-separated)\n");
    printf("  Elimination / complexity / ideal mode (otherwise):\n");
    printf("    Line 1 : variables TO ELIMINATE (comma-separated)\n");
    printf("    Line 2 : field size (prime or p^k; generator defaults to 't')\n");
    printf("    Line 3+: polynomials (comma-separated, may span multiple lines)\n");
    printf("    -> If line 1 lists n vars for n equations, compatibility mode uses the first n-1 variables\n");
    printf("\n");

    printf("MODES:\n");
    printf("  Complexity analysis:\n");
    printf("    %s --comp|--complexity \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    %s -c    \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    %s --comp -f input.dr\n", prog_name);
    printf("    Example: %s --comp \"x^2+y^2+1, x*y+z, x+y+z^2\" \"x,y\" 257\n", prog_name);
    printf("    Example: %s -c \"x^2+y^2+1, x*y+z, x+y+z^2\" \"x,y\" 257\n", prog_name);
    printf("    Example: %s --comp -f example.dr\n", prog_name);
    printf("    Example: %s -r --comp --omega 2.81 \"[2]*3\" 257\n", prog_name);
    printf("    -> Prints complexity info; saves to %s/comp_YYYYMMDD_HHMMSS.dr by default\n",
           DEFAULT_OUTPUT_DIR);
    printf("    Add --omega <value> (or -w <value>) to set omega (default: %.4g)\n",
           DIXON_OMEGA);
    printf("    Add --time to print per-step timing; use -v 2 for the old debug-level diagnostics\n");
    printf("\n");

    printf("  Dixon with ideal reduction:\n");
    printf("    %s --ideal \"ideal_generators\" \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    %s --ideal -f input.dr\n", prog_name);
    printf("    Example: %s --ideal \"a2^3=2*a1+1, a3^3=a1*a2+3\" \"a1^2+a2^2+a3^2-10, a3^3-a1*a2-3\" \"a3\" 257\n", prog_name);
    printf("    Example: %s --ideal -f example.dr\n", prog_name);
    printf("    -> ideal_generators: comma-separated relations with '=' (e.g. \"a2^3=2*a1+1, a3^3=a1*a2+3\")\n");
    printf("    -> In file mode, lines after the first two lines containing '=' are ideal generators; others are polynomials\n");
    printf("\n");

    printf("  Field-equation reduction mode (combine with any compute flag):\n");
    printf("    %s --field-equation \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    %s --field-equation input_file\n", prog_name);
    printf("    %s --field-equation -r \"[d1,d2,...,dn]\" field_size\n", prog_name);
    printf("    Example: %s --field-equation \"x0*x2+x1, x0*x1*x2+x2+1, x1*x2+x0+1\" \"x0,x1\" 2\n", prog_name);
    printf("    Example: %s --field-equation example.dr\n", prog_name);
    printf("    Example: %s --field-equation -r \"[2]*3\" 257\n", prog_name);
    printf("    -> After each multiplication, reduces x^q -> x for every variable\n");
    printf("    %s --field-equation-s \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    %s --field-equation-s input_file\n", prog_name);
    printf("    %s --field-equation-s -r \"[d1,d2,...,dn]\" field_size\n", prog_name);
    printf("    Example: %s --field-equation-s \"x0*x2+x1, x0*x1*x2+x2+1, x1*x2+x0+1\" \"x0,x1\" 2\n", prog_name);
    printf("    Example: %s --field-equation-s example.dr\n", prog_name);
    printf("    Example: %s --field-equation-s -r \"[2]*3\" 257\n", prog_name);
    printf("    -> Reduce only the final resultant/output via x^q -> x\n");
    printf("\n");

    printf("  Random mode (combine with any compute flag):\n");
    printf("    %s --random \"[d1,d2,...,dn]\" field_size\n", prog_name);
    printf("    %s -r       \"[d]*n\"          field_size\n", prog_name);
    printf("    %s -r -n 4 --density 0.5 \"[d]*3\" field_size\n", prog_name);
    printf("    %s -r -s    \"[d1,...,dn]\" field_size\n", prog_name);
    printf("    %s -r --comp  \"[d]*n\"        field_size\n", prog_name);
    printf("    Example: %s --random \"[3,3,2]\" 257\n", prog_name);
    printf("    Example: %s -r \"[3]*3\" 257\n", prog_name);
    printf("    Example: %s -r -n 4 --density 0.5 \"[3]*3\" 257\n", prog_name);
    printf("    Example: %s -r --seed 12345 \"[3]*3\" 257\n", prog_name);
    printf("    Example: %s -r \"[2]*4+[3]*2\" 257\n", prog_name);
    printf("    Example: %s -r -s \"[2]*3\" 257\n", prog_name);
    printf("    Example: %s -r --comp --omega 2.81 \"[2]*3\" 257\n", prog_name);
    printf("    -> Add -n <num_vars> to set the total variable count (must satisfy num_vars >= #equations-1)\n");
    printf("    -> Add --density <ratio> with 0 <= ratio <= 1 to choose the fraction of all monomials used (default: 1)\n");
    printf("    -> Add --seed <num> to generate the same random system reproducibly across runs\n");
    printf("    -> Mixed degree specs such as \"[2]*5+[3]*6\" are supported\n");
    printf("\n");

    printf("OPTIONS:\n");
    printf("  Verbosity:\n");
    printf("    %s -v 0 <args>\n", prog_name);
    printf("    %s -v 1 <args>\n", prog_name);
    printf("    %s -v 2 <args>\n", prog_name);
    printf("    %s -v 3 <args>\n", prog_name);
    printf("    Example: %s -v 0 \"x+y^2+t, x*y+t*y+1\" \"y\" 2^8\n", prog_name);
    printf("    Example: %s -v 1 \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    Example: %s -v 2 -s \"x^2-1, y-x\" 257\n", prog_name);
    printf("    Example: %s -v 3 \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    -> `-v 0` matches `--silent` and prints nothing\n");
    printf("    -> `-v 1` is the default output level\n");
    printf("    -> `-v 2` matches the old `--debug` output and also enables per-step timing\n");
    printf("    -> `-v 2` also prints detailed profiling for recursive Dixon construction (block counts, tuple counts, per-phase timings)\n");
    printf("    -> `-v 3` additionally prints the cancellation matrix, Dixon matrix, maximal-rank submatrix when each is <= 10 x 10, plus recursive fast-Dixon trace lines\n");
    printf("\n");

    printf("  Diagnostics:\n");
    printf("    %s --test <n>\n", prog_name);
    printf("    %s --complex \"f1, f2\" 0\n", prog_name);
    printf("    %s --time <args>\n", prog_name);
    printf("    %s -v 2 <args>\n", prog_name);
    printf("    Example: %s --test 0\n", prog_name);
    printf("    Example: %s --complex \"x^2+y^2-1, x-y\" 0\n", prog_name);
    printf("    Example: %s --time \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    Example: %s -v 2 \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    -> --test values: 0 help, 1 matrix size, 2 Bezout bound, 3 solver correctness, 4 solver performance, 5 XHash\n");
    printf("    -> --complex supports Q-only 2x2 solver mode; in elimination mode it prints complex roots when the resultant is univariate\n");
    printf("    -> --time prints per-step timing; interpolation steps also show CPU/Wall/Threads\n");
    printf("    -> `--silent`, `--debug`, `--solve-verbose` and `--solve` remain accepted for compatibility\n");
    printf("\n");

    printf("  Root search controls:\n");
    printf("    %s --rational-root-scan <auto|off|force> <args>\n", prog_name);
    printf("    %s --no-rational-root-scan <args>\n", prog_name);
    printf("    %s --force-rational-root-scan <args>\n", prog_name);
    printf("    Example: %s --rational-root-scan off \"x^2-1, y-x\" 0\n", prog_name);
    printf("    Example: %s --no-rational-root-scan \"x^2-1, y-x\" 0\n", prog_name);
    printf("    Example: %s --force-rational-root-scan \"x^2-1, y-x\" 0\n", prog_name);
    printf("    -> controls the exhaustive Rational Root Theorem scan used before approximate real-root finding\n");
    printf("    -> default `auto` skips the scan when candidate count exceeds %d; `off` disables it; `force` runs it up to the hard cap %d\n",
           FMPQ_ROOT_SEARCH_AUTO_MAX_CANDIDATES,
           FMPQ_ROOT_SEARCH_HARD_MAX_CANDIDATES);
    printf("\n");

    printf("  Method selection:\n");
    printf("    %s --method <num> <args>\n", prog_name);
    printf("    %s --fq-det-method <auto|hnf|iter> <args>\n", prog_name);
    printf("    %s --step1 <num> --step4 <num> <args>\n", prog_name);
    printf("    Example: %s --method 4 \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    Example: %s --fq-det-method hnf \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    Example: %s --step1 4 --step4 4 \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    Example: %s --cache 256 \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    Example: %s --array-limit-k 16 \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    Example: %s --fast-ksy --fast-ksy-col 0 --method 5 \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    Example: %s --step3-verify-second \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    -> Available methods: 0.Minor expansion; 1.HNF; 2.Interpolation; 3.Sparse interpolation; 4.Bareiss; 5.Recursive Dixon construction; 6.Balanced split Laplace (experimental)\n");
    printf("    -> --method sets both step 1 and step 4 for backward compatibility\n");
    printf("    -> --fq-det-method (auto|hnf|iter) controls the prime-field univariate polynomial-matrix determinant backend used in fq_poly_mat_det\n");
    printf("    -> --cache sets the determinant memoization cache entry cap (method 0 / unified expansion path)\n");
    printf("    -> --array-limit-k <k> caps optimized extension-field array multiplication tables at 2^k entries (0-62)\n");
    printf("    -> --fast-ksy enables a KSY precondition check for method 5 submatrix extraction; --no-fast-ksy disables it\n");
    printf("    -> --fast-ksy-col <idx> selects which fast-Dixon column is treated as the constant column for the KSY check (default: 0)\n");
    printf("    -> --step3-verify-second enables the second Step 3 verification pass; default is off\n");
    printf("    -> --no-step3-verify-second disables the second Step 3 verification pass\n");
    printf("\n");

    printf("  Resultant construction:\n");
    printf("    %s --dixon <args>\n", prog_name);
    printf("    %s --macaulay <args>\n", prog_name);
    printf("    %s --subres <args>\n", prog_name);
    printf("    Example: %s --dixon \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    Example: %s --macaulay \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    Example: %s --subres \"x^2+y, x+y+1\" \"x\" 257\n", prog_name);
    printf("    -> --dixon / --macaulay / --subres are direct method selectors\n");
    printf("    -> --subres is for exactly 2 polynomials and 1 elimination variable\n");
    printf("\n");

    printf("  Parallelism:\n");
    printf("    %s --threads <num> <args>\n", prog_name);
    printf("    Example: %s --threads 2 \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    -> Set number of threads for parallel computation\n");
}


static const char *display_prog_name(const char *argv0)
{
    const char *env_name = getenv("DIXON_DISPLAY_NAME");
    if (env_name && env_name[0] != '\0') return env_name;
    return argv0;
}

static int validate_cli_options(int argc, char *argv[])
{
    enum { OPT_FLAG = 1000 };
    static const struct option options[] = {
        {"silent", no_argument, NULL, OPT_FLAG},
        {"solve-verbose", no_argument, NULL, OPT_FLAG},
        {"solve-rational-only", no_argument, NULL, OPT_FLAG},
        {"complex", no_argument, NULL, OPT_FLAG},
        {"resultant-only", no_argument, NULL, OPT_FLAG},
        {"no-roots", no_argument, NULL, OPT_FLAG},
        {"approx-roots", no_argument, NULL, OPT_FLAG},
        {"root-precision", required_argument, NULL, OPT_FLAG},
        {"solve", no_argument, NULL, 's'},
        {"comp", no_argument, NULL, 'c'},
        {"complexity", no_argument, NULL, 'c'},
        {"random", no_argument, NULL, 'r'},
        {"ideal", no_argument, NULL, OPT_FLAG},
        {"field-equation", no_argument, NULL, OPT_FLAG},
        {"field-equation-s", no_argument, NULL, OPT_FLAG},
        {"time", no_argument, NULL, OPT_FLAG},
        {"no-rational-root-scan", no_argument, NULL, OPT_FLAG},
        {"force-rational-root-scan", no_argument, NULL, OPT_FLAG},
        {"rational-root-scan", required_argument, NULL, OPT_FLAG},
        {"debug", no_argument, NULL, OPT_FLAG},
        {"verbose", required_argument, NULL, 'v'},
        {"dixon", no_argument, NULL, OPT_FLAG},
        {"macaulay", no_argument, NULL, OPT_FLAG},
        {"subres", no_argument, NULL, OPT_FLAG},
        {"resultant", required_argument, NULL, OPT_FLAG},
        {"resultant-method", required_argument, NULL, OPT_FLAG},
        {"omega", required_argument, NULL, 'w'},
        {"method", required_argument, NULL, OPT_FLAG},
        {"fq-det-method", required_argument, NULL, OPT_FLAG},
        {"fast-ksy", no_argument, NULL, OPT_FLAG},
        {"ksy-precondition", no_argument, NULL, OPT_FLAG},
        {"fast-ksy-col", required_argument, NULL, OPT_FLAG},
        {"no-fast-ksy", no_argument, NULL, OPT_FLAG},
        {"step3-verify-second", no_argument, NULL, OPT_FLAG},
        {"no-step3-verify-second", no_argument, NULL, OPT_FLAG},
        {"step1", required_argument, NULL, OPT_FLAG},
        {"step4", required_argument, NULL, OPT_FLAG},
        {"threads", required_argument, NULL, OPT_FLAG},
        {"max-primes", required_argument, NULL, OPT_FLAG},
        {"cache", required_argument, NULL, OPT_FLAG},
        {"array-limit-k", required_argument, NULL, OPT_FLAG},
        {"nvars", required_argument, NULL, 'n'},
        {"num-vars", required_argument, NULL, 'n'},
        {"density", required_argument, NULL, OPT_FLAG},
        {"seed", required_argument, NULL, OPT_FLAG},
        {"specialize-vars", required_argument, NULL, OPT_FLAG},
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'V'},
        {"test", no_argument, NULL, OPT_FLAG},
        {"test-solver", no_argument, NULL, OPT_FLAG},
        {NULL, 0, NULL, 0}
    };
    char *args[argc + 1];
    int option_index = 0;
    int parsed;

    memcpy(args, argv, (size_t) argc * sizeof(*args));
    args[argc] = NULL;
    optind = 1;
    opterr = 0;

    while ((parsed = getopt_long(argc, args, "-scrv:w:n:f:o:hV", options,
                                  &option_index)) != -1) {
        if (parsed == '?') {
            const char *bad = (optind > 0 && optind <= argc) ? args[optind - 1] : "";
            fprintf(stderr, "Error: Unknown option or missing argument '%s'. Use --help to list supported options.\n",
                    bad);
            optind = 1;
            return 0;
        }
    }

    optind = 1;
    return 1;
}

int main(int argc, char *argv[])
{
    const char *prog_name = display_prog_name(argv[0]);

    if (argc == 1) {
        drsolve_cli_print_version();
        print_short_usage(prog_name);
        return 0;
    }

    if (argc == 2 &&
        (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)) {
        drsolve_cli_print_version();
        drsolve_cli_print_usage(prog_name);
        return 0;
    }

    if (!validate_cli_options(argc, argv)) {
        return 1;
    }

    if (argc == 2 &&
        (strcmp(argv[1], "--version") == 0 || strcmp(argv[1], "-V") == 0)) {
        drsolve_cli_print_version();
        return 0;
    }

    return drsolve_cli_main(argc, argv, prog_name);
}
