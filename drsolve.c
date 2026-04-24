#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include <flint/fq_nmod.h>
#include <ctype.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "dixon_flint.h"
#include "dixon_interface_flint.h"
#include "fq_mvpoly.h"
#include "fq_unified_interface.h"
#include "fq_multivariate_interpolation.h"
#include "unified_mpoly_resultant.h"
#include "dixon_with_ideal_reduction.h"
#include "dixon_complexity.h"
#include "polynomial_system_solver.h"
#include "rational_system_solver.h"
#include "dixon_test.h"

#define PROGRAM_VERSION "0.2.0"

#ifdef _WIN32
#define DIXON_NULL_DEVICE "NUL"
#else
#define DIXON_NULL_DEVICE "/dev/null"
#endif

static char *dixon_arb_to_string(const arb_t value, slong digits)
{
    char *buffer = NULL;
#ifdef _WIN32
    FILE *mem = tmpfile();
    long length;
    size_t read_size;

    if (!mem) return NULL;
    arb_fprintd(mem, value, digits);
    fflush(mem);
    if (fseek(mem, 0, SEEK_END) != 0) {
        fclose(mem);
        return NULL;
    }
    length = ftell(mem);
    if (length < 0) {
        fclose(mem);
        return NULL;
    }
    if (fseek(mem, 0, SEEK_SET) != 0) {
        fclose(mem);
        return NULL;
    }

    buffer = (char *) malloc((size_t) length + 1);
    if (!buffer) {
        fclose(mem);
        return NULL;
    }

    read_size = fread(buffer, 1, (size_t) length, mem);
    buffer[read_size] = '\0';
    fclose(mem);
#else
    size_t size = 0;
    FILE *mem = open_memstream(&buffer, &size);
    if (!mem) return NULL;
    arb_fprintd(mem, value, digits);
    fclose(mem);
#endif
    return buffer;
}

static int dixon_string_is_zeroish(const char *s)
{
    if (!s) return 0;
    while (*s && isspace((unsigned char)*s)) s++;
    if (*s == '+') s++;
    while (*s && isspace((unsigned char)*s)) s++;
    int saw_digit = 0;
    while (*s) {
        if (isdigit((unsigned char)*s)) {
            saw_digit = 1;
            if (*s != '0') return 0;
        } else if (*s == '.' || isspace((unsigned char)*s)) {
        } else {
            return 0;
        }
        s++;
    }
    return saw_digit;
}

static void dixon_fprint_arb_pretty(FILE *fp, const arb_t value, slong digits)
{
    char *raw = dixon_arb_to_string(value, digits);
    if (!raw) {
        arb_fprintd(fp, value, digits);
        return;
    }

    char *pm = strstr(raw, " +/- ");
    if (pm && (arb_is_exact(value) || dixon_string_is_zeroish(pm + 5))) {
        size_t len = (size_t)(pm - raw);
        while (len > 0 && isspace((unsigned char)raw[len - 1])) len--;
        fprintf(fp, "%.*s", (int)len, raw);
    } else {
        fputs(raw, fp);
    }
    free(raw);
}

/* =========================================================================
 * Print usage
 * ========================================================================= */

static void print_version()
{
    printf("===============================================\n");
    printf("drsolve v%s\n", PROGRAM_VERSION);
    printf("FLINT version: %s (Recommended: 3.4.0)\n", FLINT_VERSION);
#ifdef HAVE_PML
    printf("PML support: ENABLED\n");
#else
    printf("PML support: DISABLED\n");
#endif
    printf("===============================================\n");
}

static void print_usage(const char *prog_name)
{
    printf("USAGE:\n");
    printf("  Basic Dixon:\n");
    printf("    %s \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    -> Default output file: solution+timestamp.dat\n");

    printf("  Polynomial system solver:\n");
    printf("    %s \"polynomials\" field_size\n", prog_name);
    printf("    %s --solve \"polynomials\" field_size\n", prog_name);
    printf("    %s --solve-verbose \"polynomials\" field_size\n", prog_name);
    printf("    -> Writes all solutions to solution+timestamp.dat\n");
    printf("    -> `--solve` is optional here; `--solve-verbose` keeps full solver logs\n");

    printf("  Complexity analysis:\n");
    printf("    %s --comp \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    %s -c    \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    %s --comp input_file\n", prog_name);
    printf("    -> Prints complexity info; saves to comp+timestamp.dat\n");
    printf("    Add --omega <value> (or -w <value>) to set omega (default: %.4g)\n",
           DIXON_OMEGA);

    printf("  Dixon with ideal reduction:\n");
    printf("    %s --ideal \"ideal_generators\" \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    %s --ideal input_file\n", prog_name);
    printf("    -> ideal_generators: comma-separated relations with '=' (e.g. \"a2^3=2*a1+1, a3^3=a1*a2+3\")\n");
    printf("    -> In file mode, lines after the first two lines containing '=' are ideal generators; others are polynomials\n");
    
    printf("  Field-equation reduction mode (combine with any compute flag):\n");
    printf("    %s --field-equation \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    %s --field-eqution input_file\n", prog_name);
    printf("    %s --field-eqution -r \"[d1,d2,...,dn]\" field_size\n", prog_name);
    printf("    -> After each multiplication, reduces x^q -> x for every variable\n");

    printf("  Random mode (combine with any compute flag):\n");
    printf("    %s --random \"[d1,d2,...,dn]\" field_size\n", prog_name);
    printf("    %s -r       \"[d]*n\"          field_size\n", prog_name);
    printf("    %s -r --solve \"[d1,...,dn]\" field_size\n", prog_name);
    printf("    %s -r --comp  \"[d]*n\"        field_size\n", prog_name);

    printf("  File input:\n");
    printf("    %s input_file\n", prog_name);
    printf("    %s --solve input_file\n", prog_name);
    printf("    -> Without flags, auto-detects solver mode when line 1 starts with a digit; otherwise uses elimination mode\n");
    printf("    -> Output saved as input_file_solution+timestamp.dat\n");

    printf("  Silent mode:\n");
    printf("    %s --silent [--solve|--comp|-c] <args>\n", prog_name);
    printf("    -> No console output; solution file is still generated\n");

    printf("  Method selection:\n");
    printf("    %s --method <num> <args>\n", prog_name);
    printf("    %s --step1 <num> --step4 <num> <args>\n", prog_name);
    printf("    -> Available methods: 0.Recursive; 1.Kronecker+HNF; 2.Interpolation; 3.Sparse interpolation\n");
    printf("    -> --method sets both step 1 and step 4 for backward compatibility\n");

    printf("  Process count:\n");
    printf("    %s --threads <num> <args>\n", prog_name);
    printf("    -> Set number of threads for parallel computation\n");

    printf("FILE FORMAT (auto-detected for input_file):\n");
    printf("  Solver mode (line 1 starts with a digit):\n");
    printf("    Line 1 : field size\n");
    printf("    Line 2+: polynomials (one per line or comma-separated)\n");
    printf("  Elimination / complexity / ideal mode (otherwise):\n");
    printf("    Line 1 : variables TO ELIMINATE (comma-separated)\n");
    printf("    Line 2 : field size (prime or p^k; generator defaults to 't')\n");
    printf("    Line 3+: polynomials (comma-separated, may span multiple lines)\n");

    printf("EXAMPLES:\n");
    printf("  %s \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("  %s \"x^2+y^2+z^2-1, x^2+y^2-2*z^2, x+y+z\" \"x,y\" 0\n", prog_name);
    printf("  %s \"x^2+y^2+z^2-6, x+y+z-4, x*y*z-x-1\" 257\n", prog_name);
    printf("  %s --comp \"x^2+y^2+1, x*y+z, x+y+z^2\" \"x,y\" 257\n", prog_name);
    printf("  %s --random \"[3,3,2]\" 257\n", prog_name);
    printf("  %s -r \"[3]*3\" 0\n", prog_name);
    printf("  %s -r --solve \"[2]*3\" 257\n", prog_name);
    printf("  %s -r --comp --omega 2.373 \"[4]*4\" 257\n", prog_name);
    printf("  %s --ideal \"a2^3=2*a1+1, a3^3=a1*a2+3\" \"a1^2+a2^2+a3^2-10, a3^3-a1*a2-3\" \"a3\" 257\n", prog_name);
    printf("  %s --field-eqution \"x0*x2+x1, x0*x1*x2+x2+1, x1*x2+x0+1\" \"x0,x1\" 2\n", prog_name);
    printf("  %s --silent \"x+y^2+t, x*y+t*y+1\" \"x\" 2^8\n", prog_name);
    printf("  %s \"x^2 + t*y, x*y + t^2\" \"2^8: t^8 + t^4 + t^3 + t + 1\"\n", prog_name);
    printf("  (AES polynomial for GF(2^8), 't' is the field extension generator)\n");
    printf("  %s example.dr\n", prog_name);
    printf("  %s example_solve.dr\n", prog_name);
}

/* =========================================================================
 * Utility helpers
 * ========================================================================= */
static char *trim(char *str)
{
    char *end;
    while (isspace((unsigned char)*str)) str++;
    if (*str == 0) return str;
    end = str + strlen(str) - 1;
    while (end > str && isspace((unsigned char)*end)) end--;
    end[1] = '\0';
    return str;
}

static const char *display_prog_name(const char *argv0)
{
    const char *env_name = getenv("DIXON_DISPLAY_NAME");
    if (env_name && env_name[0] != '\0') return env_name;
    return argv0;
}

static const char *det_method_name_cli(int method)
{
    switch (method) {
        case 0: return "Recursive expansion";
        case 1: return "Kronecker+HNF";
        case 2: return "Interpolation";
        case 3: return "sparse interpolation";
        default: return "Default";
    }
}

static int check_prime_power(const fmpz_t n, fmpz_t prime, ulong *power)
{
    if (fmpz_cmp_ui(n, 1) <= 0) return 0;
    if (fmpz_is_probabprime(n)) {
        fmpz_set(prime, n);
        *power = 1;
        return 1;
    }
    fmpz_factor_t factors;
    fmpz_factor_init(factors);
    fmpz_factor(factors, n);
    if (factors->num == 1) {
        fmpz_set(prime, factors->p + 0);
        *power = factors->exp[0];
        fmpz_factor_clear(factors);
        return 1;
    }
    fmpz_factor_clear(factors);
    return 0;
}

static int parse_field_polynomial(nmod_poly_t modulus, const char *poly_str,
                                  mp_limb_t prime, const char *var_name)
{
    if (!poly_str || !var_name) return 0;
    nmod_poly_zero(modulus);
    char *work_str = strdup(poly_str);

    char *dst = work_str;
    for (char *src = work_str; *src; src++)
        if (!isspace(*src)) *dst++ = *src;
    *dst = '\0';

    char *token = work_str;
    while (*token) {
        if (*token == '+' || *token == '-') { token++; continue; }

        char *term_end = token;
        while (*term_end && *term_end != '+' && *term_end != '-') term_end++;

        char term_char = *term_end;
        *term_end = '\0';

        mp_limb_t coeff  = 1;
        ulong      degree = 0;
        char      *var_pos = strstr(token, var_name);

        if (!var_pos) {
            if (strlen(token) > 0) coeff = strtoul(token, NULL, 10);
            degree = 0;
        } else {
            if (var_pos > token) {
                size_t coeff_len = var_pos - token;
                char  *coeff_str = malloc(coeff_len + 1);
                strncpy(coeff_str, token, coeff_len);
                coeff_str[coeff_len] = '\0';
                size_t len = strlen(coeff_str);
                if (len > 0 && coeff_str[len - 1] == '*') coeff_str[len - 1] = '\0';
                if (strlen(coeff_str) > 0) coeff = strtoul(coeff_str, NULL, 10);
                else                        coeff = 1;
                free(coeff_str);
            }
            char *degree_pos = var_pos + strlen(var_name);
            if (*degree_pos == '^') degree = strtoul(degree_pos + 1, NULL, 10);
            else                    degree = 1;
        }

        mp_limb_t existing = nmod_poly_get_coeff_ui(modulus, degree);
        nmod_poly_set_coeff_ui(modulus, degree, (existing + coeff) % prime);

        *term_end = term_char;
        token = term_end;
    }
    free(work_str);
    return 1;
}

static int parse_field_size(const char *field_str, fmpz_t prime, ulong *power,
                            char **field_poly, char **gen_var)
{
    if (!field_str || strlen(field_str) == 0) return 0;
    if (field_poly) *field_poly = NULL;
    if (gen_var)    *gen_var    = NULL;

    if (strcmp(field_str, "0") == 0) {
        fmpz_zero(prime);
        *power = 1;
        return 1;
    }

    const char *colon = strchr(field_str, ':');
    if (colon) {
        size_t  size_len   = colon - field_str;
        char   *size_part  = malloc(size_len + 1);
        strncpy(size_part, field_str, size_len);
        size_part[size_len] = '\0';
        char *trimmed_size = trim(size_part);

        if (field_poly) {
            const char *poly_start = colon + 1;
            while (isspace(*poly_start)) poly_start++;
            *field_poly = strdup(poly_start);
            if (gen_var && *field_poly) {
                const char *p = *field_poly;
                while (*p && !isalpha(*p)) p++;
                if (*p && isalpha(*p)) {
                    char var_name[2] = { *p, '\0' };
                    *gen_var = strdup(var_name);
                }
            }
        }
        int result = parse_field_size(trimmed_size, prime, power, NULL, NULL);
        free(size_part);
        return result;
    }

    const char *caret = strchr(field_str, '^');
    if (caret) {
        char *prime_str = malloc(caret - field_str + 1);
        strncpy(prime_str, field_str, caret - field_str);
        prime_str[caret - field_str] = '\0';

        fmpz_t p;
        fmpz_init(p);
        if (fmpz_set_str(p, prime_str, 10) != 0) {
            fmpz_clear(p); free(prime_str); return 0;
        }
        char *endptr;
        ulong k = strtoul(caret + 1, &endptr, 10);
        if (*endptr != '\0' || k == 0) {
            fmpz_clear(p); free(prime_str); return 0;
        }
        if (!fmpz_is_probabprime(p)) {
            fmpz_clear(p); free(prime_str); return 0;
        }
        fmpz_set(prime, p);
        *power = k;
        fmpz_clear(p);
        free(prime_str);
        return 1;
    } else {
        fmpz_t field_size;
        fmpz_init(field_size);
        int success = fmpz_set_str(field_size, field_str, 10);
        if (success != 0) { fmpz_clear(field_size); return 0; }
        int result = check_prime_power(field_size, prime, power);
        fmpz_clear(field_size);
        return result;
    }
}

static char *generate_tagged_filename(const char *input_filename, const char *tag)
{
    if (!input_filename || !tag) return NULL;

    const char *dot = strrchr(input_filename, '.');
    size_t input_len = strlen(input_filename);
    size_t tag_len = strlen(tag);

    if (dot) {
        size_t base_len = (size_t) (dot - input_filename);
        size_t ext_len = input_len - base_len;
        char *output = malloc(base_len + tag_len + ext_len + 1);
        if (!output) return NULL;

        memcpy(output, input_filename, base_len);
        memcpy(output + base_len, tag, tag_len);
        memcpy(output + base_len + tag_len, dot, ext_len + 1);
        return output;
    }

    char *output = malloc(input_len + tag_len + 1);
    if (!output) return NULL;

    memcpy(output, input_filename, input_len);
    memcpy(output + input_len, tag, tag_len + 1);
    return output;
}

static char *generate_timestamped_filename(const char *prefix)
{
    char format[128];
    char buffer[128];
    time_t now = time(NULL);
    struct tm *t = localtime(&now);

    if (t) {
        snprintf(format, sizeof(format), "%s_%%Y%%m%%d_%%H%%M%%S.dr", prefix);
        strftime(buffer, sizeof(buffer), format, t);
    } else {
        snprintf(buffer, sizeof(buffer), "%s.dr", prefix);
    }

    return strdup(buffer);
}

static void print_field_label(FILE *out, const fmpz_t prime, ulong power)
{
    if (fmpz_is_zero(prime)) {
        fprintf(out, "Q");
        return;
    }

    fprintf(out, "F_");
    fmpz_fprint(out, prime);
    if (power > 1) {
        fprintf(out, "^%lu", power);
        if (fmpz_abs_fits_ui(prime)) {
            mp_limb_t prime_ui = fmpz_get_ui(prime);
            mp_limb_t field_size = 1;
            int overflow = 0;

            for (ulong i = 0; i < power; i++) {
                if (prime_ui != 0 && field_size > UWORD_MAX / prime_ui) {
                    overflow = 1;
                    break;
                }
                field_size *= prime_ui;
            }

            if (!overflow) {
                fprintf(out, " (size %lu)", field_size);
            }
        }
    }
}

/* =========================================================================
 * Complexity analysis: save to file
 * ========================================================================= */
static void save_comp_result_to_file(
        const char   *filename,
        const char   *polys_str,
        const char   *vars_str,
        const fmpz_t  prime,
        ulong         power,
        slong         num_polys,
        slong         num_all_vars,
        char        **all_vars,
        char        **elim_var_list,
        slong         num_elim,
        const long   *degrees,
        const fmpz_t  matrix_size,
        long          bezout_bound,
        double        complexity,
        double        omega,
        double        comp_time)
{
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Warning: Cannot create output file '%s'\n", filename);
        return;
    }

    fprintf(fp, "Dixon Complexity Analysis\n");
    fprintf(fp, "=========================\n");
    fprintf(fp, "Field: ");
    print_field_label(fp, prime, power);
    fprintf(fp, "\n");
    fprintf(fp, "Polynomials: %s\n", polys_str);
    fprintf(fp, "Eliminate:   %s\n", vars_str);
    fprintf(fp, "Computation time: %.3f seconds\n\n", comp_time);

    fprintf(fp, "--- Input Summary ---\n");
    fprintf(fp, "Equations : %ld\n", num_polys);
    fprintf(fp, "Variables (%ld): ", num_all_vars);
    for (slong i = 0; i < num_all_vars; i++) {
        if (i > 0) fprintf(fp, ", ");
        fprintf(fp, "%s", all_vars[i]);
    }
    fprintf(fp, "\n");
    fprintf(fp, "Elim vars (%ld): ", num_elim);
    for (slong i = 0; i < num_elim; i++) {
        if (i > 0) fprintf(fp, ", ");
        fprintf(fp, "%s", elim_var_list[i]);
    }
    fprintf(fp, "\n");

    /* Remaining vars */
    fprintf(fp, "Remaining vars: ");
    int first = 1;
    for (slong i = 0; i < num_all_vars; i++) {
        int is_elim = 0;
        for (slong j = 0; j < num_elim; j++)
            if (strcmp(all_vars[i], elim_var_list[j]) == 0) { is_elim = 1; break; }
        if (!is_elim) {
            if (!first) fprintf(fp, ", ");
            fprintf(fp, "%s", all_vars[i]);
            first = 0;
        }
    }
    if (first) fprintf(fp, "(none)");
    fprintf(fp, "\n\n");

    fprintf(fp, "--- Degree & Size ---\n");
    fprintf(fp, "Degree sequence: [");
    for (slong i = 0; i < num_polys; i++) {
        if (i > 0) fprintf(fp, ", ");
        fprintf(fp, "%ld", degrees[i]);
    }
    fprintf(fp, "]\n");
    fprintf(fp, "Bezout bound (degree product): %ld\n", bezout_bound);
    fprintf(fp, "Dixon matrix size: ");
    fmpz_fprint(fp, matrix_size);
    fprintf(fp, "\n");
    fprintf(fp, "Resultant degree estimate (Bezout): %ld\n", bezout_bound);
    fprintf(fp, "\n--- Complexity ---\n");
    fprintf(fp, "Complexity (log2, omega=%.4g): %.6f\n", omega, complexity);

    fclose(fp);
}

/* =========================================================================
 * Complexity analysis: run and display
 * ========================================================================= */
static void run_complexity_analysis(
        const char      *polys_str,
        const char      *vars_str,
        const fmpz_t     prime,
        ulong            power,
        const fq_nmod_ctx_t ctx,
        const char      *output_filename,
        int              silent_mode,
        double           comp_time,
        double           omega)
{
    /* ---- split polynomials and elimination variables ---- */
    slong   num_polys;
    char  **poly_arr = split_string(polys_str, &num_polys);

    slong   num_elim;
    char  **elim_arr = split_string(vars_str, &num_elim);

    /* ---- get generator name ---- */
    char *gen_name = (ctx == NULL) ? NULL : get_generator_name(ctx);

    /* ---- collect all variables ---- */
    char  **all_vars;
    slong   num_all_vars;
    collect_variables((const char **)poly_arr, num_polys,
                      gen_name, &all_vars, &num_all_vars);

    /* ---- compute degree of each polynomial ---- */
    if (num_polys <= 0) {
        if (!silent_mode) fprintf(stderr, "Error: no polynomials to analyze\n");
        free_split_strings(poly_arr, num_polys);
        free_split_strings(elim_arr, num_elim);
        if (gen_name) free(gen_name);
        for (slong i = 0; i < num_all_vars; i++) free(all_vars[i]);
        free(all_vars);
        return;
    }
    long *degrees = calloc((size_t)num_polys, sizeof(long));
    for (slong i = 0; i < num_polys; i++)
        degrees[i] = get_poly_total_degree(poly_arr[i], gen_name);

    /* ---- Bezout bound = product of degrees ---- */
    long bezout = 1;
    for (slong i = 0; i < num_polys; i++) bezout *= degrees[i];

    /* ---- Dixon matrix size via Hessenberg recurrence ---- */
    fmpz_t matrix_size;
    fmpz_init(matrix_size);
    /* suppress internal prints from dixon_size */
    {
        fflush(stdout);
        int orig_stdout = dup(STDOUT_FILENO);
        int devnull     = open(DIXON_NULL_DEVICE, O_WRONLY);
        if (devnull != -1) { dup2(devnull, STDOUT_FILENO); close(devnull); }
        dixon_size(matrix_size, degrees, (int)num_polys, 0);
        fflush(stdout);
        dup2(orig_stdout, STDOUT_FILENO);
        close(orig_stdout);
    }

    /* ---- Dixon complexity ---- */
    double complexity = dixon_complexity(degrees, (int)num_polys,
                                        (int)num_all_vars, omega);

    /* ---- console output (concise) ---- */
    if (!silent_mode) {
        printf("\n=== Complexity Analysis ===\n");
        printf("Equations: %ld  |  Variables: %ld  |  Eliminate: %ld\n",
               num_polys, num_all_vars, num_elim);

        printf("All vars : ");
        for (slong i = 0; i < num_all_vars; i++) {
            if (i > 0) printf(", ");
            printf("%s", all_vars[i]);
        }
        printf("\n");

        printf("Degrees  : [");
        for (slong i = 0; i < num_polys; i++) {
            if (i > 0) printf(", ");
            printf("%ld", degrees[i]);
        }
        printf("]\n");

        printf("Bezout bound      : %ld\n", bezout);

        printf("Dixon matrix size : ");
        fmpz_print(matrix_size);
        printf("\n");

        printf("Resultant deg est : %ld  (Bezout bound)\n", bezout);

        if (isfinite(complexity))
            printf("Complexity (log2) : %.4f  (omega=%.4g)\n",
                   complexity, omega);
        else
            printf("Complexity (log2) : inf / undefined\n");

        if (output_filename)
            printf("Report saved to   : %s\n", output_filename);
        printf("===========================\n");
    }

    /* ---- save to file ---- */
    if (output_filename) {
        save_comp_result_to_file(
            output_filename, polys_str, vars_str,
            prime, power,
            num_polys, num_all_vars, all_vars,
            elim_arr, num_elim,
            degrees, matrix_size, bezout,
            complexity, omega, comp_time);
    }

    /* ---- cleanup ---- */
    fmpz_clear(matrix_size);
    free(degrees);
    for (slong i = 0; i < num_all_vars; i++) free(all_vars[i]);
    free(all_vars);
    if (gen_name) free(gen_name);
    free_split_strings(poly_arr, num_polys);
    free_split_strings(elim_arr, num_elim);
}


/* =========================================================================
 * Random polynomial generation helpers
 * ========================================================================= */

/*
 * Parse degree list in any of:
 *   "3,3,2"          plain comma-separated
 *   "[3,3,2]"        with brackets
 *   "[3]*5"          repeat shorthand → [3,3,3,3,3]
 *   "[3,2]*4"        repeat a sub-list → [3,2,3,2,3,2,3,2]
 *
 * Returns malloc'd array; caller frees.  *count_out = number of entries.
 */
static long *parse_degree_list(const char *deg_str, slong *count_out)
{
    slong cap = 32;
    long *degs = malloc((size_t)cap * sizeof(long));
    slong count = 0;

    /* Helper: append one value */
#define PUSH(v) do { \
    if (count >= cap) { cap *= 2; degs = realloc(degs, (size_t)cap * sizeof(long)); } \
    degs[count++] = (v); } while(0)

    /* Strip optional outer brackets and scan for "[...]*N" repeat syntax */
    const char *p = deg_str;
    while (isspace((unsigned char)*p)) p++;

    /* Check for [...]* repeat */
    if (*p == '[') {
        const char *bracket_start = p + 1;
        /* find matching ']' */
        const char *bracket_end = strchr(bracket_start, ']');
        if (bracket_end) {
            /* collect base list between brackets */
            slong base_cap = 16;
            long *base = malloc((size_t)base_cap * sizeof(long));
            slong base_n = 0;

            const char *q = bracket_start;
            while (q < bracket_end) {
                while (q < bracket_end &&
                       (isspace((unsigned char)*q) || *q == ',')) q++;
                if (q >= bracket_end) break;
                char *end;
                long d = strtol(q, &end, 10);
                if (end == q) break;
                if (d > 0) {
                    if (base_n >= base_cap) { base_cap *= 2; base = realloc(base, (size_t)base_cap * sizeof(long)); }
                    base[base_n++] = d;
                }
                q = end;
            }

            /* check for '*N' after ']' */
            const char *after = bracket_end + 1;
            while (isspace((unsigned char)*after)) after++;
            long repeat = 1;
            if (*after == '*') {
                char *end;
                repeat = strtol(after + 1, &end, 10);
                if (repeat <= 0) repeat = 1;
            } else if (base_n == 1) {
                /* "[d]*N" only makes sense with an explicit repeat;
                   without '*', treat as a regular bracketed list        */
                repeat = 1;
            }

            for (long r = 0; r < repeat; r++)
                for (slong j = 0; j < base_n; j++)
                    PUSH(base[j]);

            free(base);
            *count_out = count;
            return degs;
        }
        /* malformed bracket — fall through to plain parsing */
        p++;  /* skip the '[' */
    }

    /* Plain comma-separated parsing (also handles "[3,3,2]" without '*') */
    while (*p) {
        /* skip non-digit (brackets, spaces, commas) */
        while (*p && !isdigit((unsigned char)*p) && *p != '-') p++;
        if (!*p) break;
        char *end;
        long d = strtol(p, &end, 10);
        if (end == p) { p++; continue; }
        if (d > 0) PUSH(d);
        else fprintf(stderr, "Warning: degree %ld ignored (must be > 0)\n", d);
        p = end;
    }
#undef PUSH

    *count_out = count;
    return degs;
}


/*
 * Generate a random polynomial system and return it as strings.
 *
 * solve_mode=0 (Dixon/comp): n polys in n vars; elim_vars = x0..x_{n-2}
 * solve_mode=1 (solver):     n polys in n vars; elim_vars = all (solver ignores it)
 *
 * All three output strings are malloc'd; caller must free them.
 * Returns 1 on success, 0 on failure.
 */
static int generate_random_poly_strings(
        const long *degrees, slong npolys,
        const fq_nmod_ctx_t ctx,
        int is_solve_mode,
        int silent_mode,
        char **polys_str_out,
        char **elim_vars_str_out,
        char **all_vars_str_out)
{
    slong nvars = npolys;
    slong npars = 0;

    /* ---- x0, x1, ..., x_{n-1} ---- */
    char **var_names = malloc((size_t)nvars * sizeof(char *));
    for (slong i = 0; i < nvars; i++) {
        var_names[i] = malloc(16);
        snprintf(var_names[i], 16, "x%ld", i);
    }

    /* ---- all_vars ---- */
    size_t av_cap = (size_t)nvars * 8 + 2;
    char  *all_vars = malloc(av_cap);
    all_vars[0] = '\0';
    for (slong i = 0; i < nvars; i++) {
        if (i > 0) strcat(all_vars, ",");
        strcat(all_vars, var_names[i]);
    }

    /* ---- elim_vars ---- */
    char *elim_vars;
    if (!is_solve_mode && nvars >= 2) {
        size_t ev_cap = (size_t)(nvars - 1) * 8 + 2;
        elim_vars = malloc(ev_cap);
        elim_vars[0] = '\0';
        for (slong i = 0; i < nvars - 1; i++) {
            if (i > 0) strcat(elim_vars, ",");
            strcat(elim_vars, var_names[i]);
        }
    } else {
        elim_vars = strdup(all_vars);
    }

    slong *slong_deg = malloc((size_t)npolys * sizeof(slong));
    for (slong i = 0; i < npolys; i++) slong_deg[i] = (slong)degrees[i];

    flint_rand_t rstate;
    flint_rand_init(rstate);
    flint_rand_set_seed(rstate, (ulong)time(NULL),
                        (ulong)time(NULL) ^ (ulong)clock());

    if (!silent_mode) printf("Generating random polynomial system...\n");

    fflush(stdout);
    int orig_stdout = dup(STDOUT_FILENO);
    int devnull = open(DIXON_NULL_DEVICE, O_WRONLY);
    if (devnull != -1) { dup2(devnull, STDOUT_FILENO); close(devnull); }

    fq_mvpoly_t *polys = NULL;
    generate_polynomial_system(&polys, nvars, npolys, npars,
                               slong_deg, 0.8, ctx, rstate);

    fflush(stdout);
    dup2(orig_stdout, STDOUT_FILENO);
    close(orig_stdout);

    char *gen_name = get_generator_name(ctx);

    size_t total_len = 4;
    char **poly_strs = malloc((size_t)npolys * sizeof(char *));
    for (slong i = 0; i < npolys; i++) {
        char *s = fq_mvpoly_to_string(&polys[i], NULL, gen_name);
        poly_strs[i] = (s && strlen(s) > 0) ? s : (free(s), strdup("0"));
        total_len += strlen(poly_strs[i]) + 3;
    }

    char *polys_str = malloc(total_len);
    polys_str[0] = '\0';
    for (slong i = 0; i < npolys; i++) {
        if (i > 0) strcat(polys_str, ", ");
        strcat(polys_str, poly_strs[i]);
        free(poly_strs[i]);
    }
    free(poly_strs);

    if (!silent_mode) {
        printf("System: %ld equations, %ld variables\n", npolys, nvars);
        printf("Degrees: [");
        for (slong i = 0; i < npolys; i++) {
            if (i > 0) printf(", ");
            printf("%ld", slong_deg[i]);
        }
        printf("]\n");
        if (!is_solve_mode)
            printf("Eliminate: %s  |  Remaining: %s\n",
                   elim_vars, var_names[nvars - 1]);
    }

    /* ---- clear ---- */
    for (slong i = 0; i < npolys; i++) fq_mvpoly_clear(&polys[i]);
    free(polys);
    free(slong_deg);
    if (gen_name) free(gen_name);
    for (slong i = 0; i < nvars; i++) free(var_names[i]);
    free(var_names);
    flint_rand_clear(rstate);

    *polys_str_out      = polys_str;
    *elim_vars_str_out  = elim_vars;
    *all_vars_str_out   = all_vars;
    return 1;
}

/* =========================================================================
 * File reading helpers
 * ========================================================================= */
static char *read_entire_line(FILE *fp)
{
    if (!fp) return NULL;
    size_t capacity = 4096, length = 0;
    char  *line = malloc(capacity);
    if (!line) return NULL;
    int c;
    while ((c = fgetc(fp)) != EOF && c != '\n' && c != '\r') {
        if (length + 1 >= capacity) {
            capacity *= 2;
            char *nl = realloc(line, capacity);
            if (!nl) { free(line); return NULL; }
            line = nl;
        }
        line[length++] = (char)c;
    }
    if (c == '\r') {
        int next = fgetc(fp);
        if (next != '\n' && next != EOF) ungetc(next, fp);
    }
    if (length == 0 && c == EOF) { free(line); return NULL; }
    line[length] = '\0';
    return line;
}

static void free_input_lines(char **lines, int line_count)
{
    if (!lines) return;
    for (int i = 0; i < line_count; i++) free(lines[i]);
    free(lines);
}

static int collect_input_lines(FILE *fp, char ***lines_out, int *line_count_out)
{
    char **lines = NULL;
    int line_count = 0, line_capacity = 10;

    lines = malloc((size_t) line_capacity * sizeof(char *));
    if (!lines) return 0;

    char *line;
    while ((line = read_entire_line(fp)) != NULL) {
        char *trimmed = trim(line);
        if (strlen(trimmed) == 0 || trimmed[0] == '#') {
            free(line);
            continue;
        }
        if (line_count >= line_capacity) {
            line_capacity *= 2;
            char **nl = realloc(lines, (size_t) line_capacity * sizeof(char *));
            if (!nl) {
                free_input_lines(lines, line_count);
                free(line);
                return 0;
            }
            lines = nl;
        }
        lines[line_count++] = strdup(trimmed);
        free(line);
    }

    *lines_out = lines;
    *line_count_out = line_count;
    return 1;
}

static char *join_polynomial_lines(char **lines, int start, int end)
{
    size_t total_len = 0;
    for (int i = start; i < end; i++)
        total_len += strlen(lines[i]) + 3;

    char *poly_buffer = malloc(total_len + 1);
    if (!poly_buffer) return NULL;

    poly_buffer[0] = '\0';
    for (int i = start; i < end; i++) {
        if (i > start) {
            size_t prev_len = strlen(poly_buffer);
            int prev_comma = (prev_len > 0 &&
                (poly_buffer[prev_len - 1] == ',' ||
                 (prev_len > 1 && poly_buffer[prev_len - 2] == ',' &&
                  poly_buffer[prev_len - 1] == ' ')));
            int curr_comma = (lines[i][0] == ',');
            if (!prev_comma && !curr_comma)      strcat(poly_buffer, ", ");
            else if (!prev_comma && curr_comma)  strcat(poly_buffer, " ");
        }
        strcat(poly_buffer, lines[i]);
    }

    return poly_buffer;
}

static int line_starts_with_digit(const char *line)
{
    return (line && line[0] != '\0' && isdigit((unsigned char) line[0]));
}

static int read_multiline_file(FILE *fp, char **field_str, char **polys_str,
                               char **vars_str, char **ideal_str,
                               char **allvars_str)
{
    char **lines = NULL;
    int line_count = 0;

    if (!collect_input_lines(fp, &lines, &line_count)) return 0;

    if (line_count < 3) {
        fprintf(stderr, "Error: Elimination file must contain at least 3 non-empty lines\n");
        fprintf(stderr, "  Line 1 : variables to ELIMINATE\n");
        fprintf(stderr, "  Line 2 : field size\n");
        fprintf(stderr, "  Line 3+: polynomials\n");
        free_input_lines(lines, line_count);
        return 0;
    }

    *vars_str  = lines[0];
    *field_str = lines[1];
    *polys_str = join_polynomial_lines(lines, 2, line_count);
    if (!*polys_str) {
        free_input_lines(lines, line_count);
        return 0;
    }

    for (int i = 2; i < line_count; i++) free(lines[i]);
    free(lines);

    *ideal_str    = NULL;
    *allvars_str  = NULL;
    return 1;
}

static int read_solver_file(FILE *fp, char **field_str, char **polys_str)
{
    char **lines = NULL;
    int line_count = 0;

    if (!collect_input_lines(fp, &lines, &line_count)) return 0;

    if (line_count < 2) {
        fprintf(stderr, "Error: Solver file must contain at least 2 non-empty lines\n");
        fprintf(stderr, "  Line 1 : field size\n");
        fprintf(stderr, "  Line 2+: polynomials\n");
        free_input_lines(lines, line_count);
        return 0;
    }

    *field_str = lines[0];
    *polys_str = join_polynomial_lines(lines, 1, line_count);
    if (!*polys_str) {
        free_input_lines(lines, line_count);
        return 0;
    }

    for (int i = 1; i < line_count; i++) free(lines[i]);
    free(lines);
    return 1;
}

/* =========================================================================
 * Read file for --ideal mode:
 *   Line 1 : variables to ELIMINATE
 *   Line 2 : field size
 *   Line 3+: polys (no '=') or ideal generators (has '='), mixed
 * ========================================================================= */
static int read_ideal_file(FILE *fp, char **field_str, char **polys_str,
                           char **vars_str, char **ideal_str)
{
    char **lines = NULL;
    int line_count = 0;

    if (!collect_input_lines(fp, &lines, &line_count)) return 0;

    if (line_count < 3) {
        fprintf(stderr, "Error: --ideal file needs at least 3 non-empty lines\n");
        fprintf(stderr, "  Line 1 : variables to ELIMINATE\n");
        fprintf(stderr, "  Line 2 : field size\n");
        fprintf(stderr, "  Line 3+: polynomials and/or ideal generators (lines with '=' are ideal)\n");
        free_input_lines(lines, line_count);
        return 0;
    }

    *vars_str  = lines[0];
    *field_str = lines[1];

    size_t poly_len  = 0, ideal_len = 0;
    for (int i = 2; i < line_count; i++) {
        if (strchr(lines[i], '='))  ideal_len += strlen(lines[i]) + 3;
        else                         poly_len  += strlen(lines[i]) + 3;
    }

    char *poly_buf  = malloc(poly_len  + 4);
    char *ideal_buf = malloc(ideal_len + 4);
    if (!poly_buf || !ideal_buf) {
        free(poly_buf); free(ideal_buf);
        for (int i = 0; i < line_count; i++) free(lines[i]);
        free(lines);
        return 0;
    }
    poly_buf[0]  = '\0';
    ideal_buf[0] = '\0';

    int poly_first = 1, ideal_first = 1;
    for (int i = 2; i < line_count; i++) {
        if (strchr(lines[i], '=')) {
            if (!ideal_first) strcat(ideal_buf, ", ");
            strcat(ideal_buf, lines[i]);
            ideal_first = 0;
        } else {
            if (!poly_first) strcat(poly_buf, ", ");
            strcat(poly_buf, lines[i]);
            poly_first = 0;
        }
        free(lines[i]);
    }

    *polys_str = poly_buf;
    *ideal_str = (ideal_len > 0) ? ideal_buf : (free(ideal_buf), NULL);

    free(lines);
    return 1;
}

static int read_auto_input_file(FILE *fp, int *solve_mode_out,
                                char **field_str, char **polys_str,
                                char **vars_str, char **ideal_str,
                                char **allvars_str)
{
    char **lines = NULL;
    int line_count = 0;

    if (!collect_input_lines(fp, &lines, &line_count)) return 0;
    if (line_count < 2) {
        fprintf(stderr, "Error: Input file must contain at least 2 non-empty lines\n");
        fprintf(stderr, "  Solver mode      : line 1 = field size, line 2+ = polynomials\n");
        fprintf(stderr, "  Elimination mode : line 1 = elimination vars, line 2 = field size, line 3+ = polynomials\n");
        free_input_lines(lines, line_count);
        return 0;
    }

    *ideal_str = NULL;
    *allvars_str = NULL;

    if (line_starts_with_digit(lines[0])) {
        *solve_mode_out = 1;
        *vars_str = NULL;
        *field_str = lines[0];
        *polys_str = join_polynomial_lines(lines, 1, line_count);
        if (!*polys_str) {
            free_input_lines(lines, line_count);
            return 0;
        }
        for (int i = 1; i < line_count; i++) free(lines[i]);
        free(lines);
        return 1;
    }

    if (line_count < 3) {
        fprintf(stderr, "Error: Elimination file must contain at least 3 non-empty lines\n");
        fprintf(stderr, "  Line 1 : variables to ELIMINATE\n");
        fprintf(stderr, "  Line 2 : field size\n");
        fprintf(stderr, "  Line 3+: polynomials\n");
        free_input_lines(lines, line_count);
        return 0;
    }

    *solve_mode_out = 0;
    *vars_str = lines[0];
    *field_str = lines[1];
    *polys_str = join_polynomial_lines(lines, 2, line_count);
    if (!*polys_str) {
        free_input_lines(lines, line_count);
        return 0;
    }
    for (int i = 2; i < line_count; i++) free(lines[i]);
    free(lines);
    return 1;
}

/* =========================================================================
 * Result saving helpers (unchanged from original)
 * ========================================================================= */
static void save_solver_result_to_file(const char *filename,
                                       const char *polys_str,
                                       const fmpz_t prime, ulong power,
                                       const polynomial_solutions_t *sols,
                                       double cpu_time, double wall_time, int threads_num)
{
    FILE *out_fp = fopen(filename, "w");
    if (!out_fp) {
        fprintf(stderr, "Warning: Could not create output file '%s'\n", filename);
        return;
    }

    fprintf(out_fp, "Polynomial System Solver\n");
    fprintf(out_fp, "========================\n");
    fprintf(out_fp, "Field: ");
    print_field_label(out_fp, prime, power);
    if (!fmpz_is_zero(prime) && power > 1) {
        fprintf(out_fp, "\nField extension generator: t");
    }
    fprintf(out_fp, "\n");
    fprintf(out_fp, "Polynomials: %s\n", polys_str);
    fprintf(out_fp, "CPU time: %.3f seconds | Wall time: %.3f seconds | Threads: %d\n", cpu_time, wall_time, threads_num);
    fprintf(out_fp, "\nSolutions:\n==========\n");

    if (!sols) { fprintf(out_fp, "Solution structure is null\n"); fclose(out_fp); return; }

    fprintf(out_fp, "\n=== Polynomial System Solutions ===\n");
    if (!sols->is_valid) {
        fprintf(out_fp, "Solving failed");
        if (sols->error_message) fprintf(out_fp, ": %s", sols->error_message);
        fprintf(out_fp, "\n");
        fclose(out_fp); return;
    }
    if (sols->has_no_solutions == -1) {
        fprintf(out_fp, "System has positive dimension; finite solution listing skipped\n");
        fclose(out_fp); return;
    }
    if (sols->has_no_solutions == 1) {
        fprintf(out_fp, "System has no solutions over the finite field\n");
        fclose(out_fp); return;
    }
    if (sols->num_variables == 0) {
        fprintf(out_fp, "No variables\n"); fclose(out_fp); return;
    }
    if (sols->num_solution_sets == 0) {
        fprintf(out_fp, "No solutions found\n"); fclose(out_fp); return;
    }

    fprintf(out_fp, "Found %ld complete solution set(s):\n", sols->num_solution_sets);
    for (slong set = 0; set < sols->num_solution_sets; set++) {
        fprintf(out_fp, "\nSolution set %ld:\n", set + 1);
        for (slong var = 0; var < sols->num_variables; var++) {
            fprintf(out_fp, "  %s = ", sols->variable_names[var]);
            slong num_sols = sols->solutions_per_var[set * sols->num_variables + var];
            if (num_sols == 0) {
                fprintf(out_fp, "no solution");
            } else if (num_sols == 1) {
                char *sol_str = fq_nmod_get_str_pretty(sols->solution_sets[set][var][0], sols->ctx);
                fprintf(out_fp, "%s", sol_str); free(sol_str);
            } else {
                fprintf(out_fp, "{");
                for (slong sol = 0; sol < num_sols; sol++) {
                    if (sol > 0) fprintf(out_fp, ", ");
                    char *sol_str = fq_nmod_get_str_pretty(sols->solution_sets[set][var][sol], sols->ctx);
                    fprintf(out_fp, "%s", sol_str); free(sol_str);
                }
                fprintf(out_fp, "}");
            }
            fprintf(out_fp, "\n");
        }
    }

    fprintf(out_fp, "\n=== Compatibility View ===\n");
    for (slong var = 0; var < sols->num_variables; var++) {
        fprintf(out_fp, "%s = {", sols->variable_names[var]);
        slong total_printed = 0;
        for (slong set = 0; set < sols->num_solution_sets; set++) {
            slong num_sols = sols->solutions_per_var[set * sols->num_variables + var];
            for (slong sol = 0; sol < num_sols; sol++) {
                if (total_printed > 0) fprintf(out_fp, ", ");
                char *sol_str = fq_nmod_get_str_pretty(sols->solution_sets[set][var][sol], sols->ctx);
                fprintf(out_fp, "%s", sol_str); free(sol_str);
                total_printed++;
            }
        }
        fprintf(out_fp, "}");
        if (total_printed > 1)      fprintf(out_fp, " (%ld solutions)", total_printed);
        else if (total_printed == 0) fprintf(out_fp, " (no solutions)");
        fprintf(out_fp, "\n");
    }
    fprintf(out_fp, "=== Solution Complete ===\n\n");
    fclose(out_fp);
}

static void save_rational_solver_result_to_file(const char *filename,
                                                const char *polys_str,
                                                const rational_solutions_t *sols,
                                                double cpu_time, double wall_time, int threads_num)
{
    FILE *out_fp = fopen(filename, "w");
    if (!out_fp) {
        fprintf(stderr, "Warning: Could not create output file '%s'\n", filename);
        return;
    }

    fprintf(out_fp, "Rational Polynomial System Solver\n");
    fprintf(out_fp, "==================================\n");
    fprintf(out_fp, "Field: Q\n");
    fprintf(out_fp, "Polynomials: %s\n", polys_str);
    fprintf(out_fp, "CPU time: %.3f seconds | Wall time: %.3f seconds | Threads: %d\n", cpu_time, wall_time, threads_num);
    fprintf(out_fp, "\nSolutions:\n==========\n");

    if (!sols) { fprintf(out_fp, "Solution structure is null\n"); fclose(out_fp); return; }

    fprintf(out_fp, "\n=== Rational Polynomial System Solutions ===\n");
    if (!sols->is_valid) {
        fprintf(out_fp, "Solving failed");
        if (sols->error_message) fprintf(out_fp, ": %s", sols->error_message);
        fprintf(out_fp, "\n");
        fclose(out_fp); return;
    }
    if (sols->has_no_solutions == -1) {
        fprintf(out_fp, "System has positive dimension; finite solution listing skipped\n");
        fclose(out_fp); return;
    }
    if (sols->has_no_solutions == 1) {
        if (sols->real_root_summary && sols->num_solution_sets == 0 && sols->num_real_solution_sets == 0) {
            fprintf(out_fp, "No exact rational solution sets were assembled\n");
            fprintf(out_fp, "%s\n", sols->real_root_summary);
        } else {
            fprintf(out_fp, "System has no solutions over the rational numbers\n");
        }
        fclose(out_fp); return;
    }
    if (sols->num_variables == 0) {
        fprintf(out_fp, "No variables\n"); fclose(out_fp); return;
    }
    if (sols->num_solution_sets > 0) {
        fprintf(out_fp, "Found %ld exact rational solution set(s):\n", sols->num_solution_sets);
        for (slong set = 0; set < sols->num_solution_sets; set++) {
            fprintf(out_fp, "\nSolution set %ld:\n", set + 1);
            for (slong var = 0; var < sols->num_variables; var++) {
                fprintf(out_fp, "  %s = ", sols->variable_names[var]);
                slong num_sols = sols->solutions_per_var[set * sols->num_variables + var];
                if (num_sols == 0) {
                    fprintf(out_fp, "no solution");
                } else if (num_sols == 1) {
                    char *sol_str = fmpq_get_str(NULL, 10, sols->solution_sets[set][var][0]);
                    fprintf(out_fp, "%s", sol_str);
                    flint_free(sol_str);
                } else {
                    fprintf(out_fp, "{");
                    for (slong sol = 0; sol < num_sols; sol++) {
                        if (sol > 0) fprintf(out_fp, ", ");
                        char *sol_str = fmpq_get_str(NULL, 10, sols->solution_sets[set][var][sol]);
                        fprintf(out_fp, "%s", sol_str);
                        flint_free(sol_str);
                    }
                    fprintf(out_fp, "}");
                }
                fprintf(out_fp, "\n");
            }
            if (sols->solution_residuals && sols->solution_residuals[set]) {
                fprintf(out_fp, "  residuals: [");
                for (slong eq = 0; eq < sols->num_equations; eq++) {
                    if (eq > 0) fprintf(out_fp, ", ");
                    fprintf(out_fp, "eq%ld=", eq + 1);
                    dixon_fprint_arb_pretty(out_fp, sols->solution_residuals[set][eq], 8);
                }
                fprintf(out_fp, "]\n");
            }
        }
    }

    if (sols->num_real_solution_sets > 0) {
        fprintf(out_fp, "Found %ld approximate real solution set(s):\n", sols->num_real_solution_sets);
        for (slong set = 0; set < sols->num_real_solution_sets; set++) {
            fprintf(out_fp, "\nSolution set %ld:\n", set + 1);
            for (slong var = 0; var < sols->num_variables; var++) {
                fprintf(out_fp, "  %s = ", sols->variable_names[var]);
                slong num_sols = sols->real_solutions_per_var[set * sols->num_variables + var];
                if (num_sols == 0) {
                    fprintf(out_fp, "no solution");
                } else if (num_sols == 1) {
                    dixon_fprint_arb_pretty(out_fp, sols->real_solution_sets[set][var][0], 20);
                } else {
                    fprintf(out_fp, "{");
                    for (slong sol = 0; sol < num_sols; sol++) {
                        if (sol > 0) fprintf(out_fp, ", ");
                        dixon_fprint_arb_pretty(out_fp, sols->real_solution_sets[set][var][sol], 20);
                    }
                    fprintf(out_fp, "}");
                }
                fprintf(out_fp, "\n");
            }
            if (sols->real_solution_residuals && sols->real_solution_residuals[set]) {
                fprintf(out_fp, "  residuals: [");
                for (slong eq = 0; eq < sols->num_equations; eq++) {
                    if (eq > 0) fprintf(out_fp, ", ");
                    fprintf(out_fp, "eq%ld=", eq + 1);
                    dixon_fprint_arb_pretty(out_fp, sols->real_solution_residuals[set][eq], 8);
                }
                fprintf(out_fp, "]\n");
            }
        }
        if (sols->real_root_summary) {
            fprintf(out_fp, "\n%s\n", sols->real_root_summary);
        }
    }
    if (sols->num_solution_sets == 0 && sols->num_real_solution_sets == 0) {
        if (sols->real_root_summary) {
            fprintf(out_fp, "%s\n", sols->real_root_summary);
        } else {
            fprintf(out_fp, "No solutions found\n");
        }
        fclose(out_fp); return;
    }

    fprintf(out_fp, "\n=== Projection View (all solution sets) ===\n");
    for (slong var = 0; var < sols->num_variables; var++) {
        fprintf(out_fp, "%s = {", sols->variable_names[var]);
        slong total_printed = 0;
        for (slong set = 0; set < sols->num_solution_sets; set++) {
            slong num_sols = sols->solutions_per_var[set * sols->num_variables + var];
            for (slong sol = 0; sol < num_sols; sol++) {
                if (total_printed > 0) fprintf(out_fp, ", ");
                char *sol_str = fmpq_get_str(NULL, 10, sols->solution_sets[set][var][sol]);
                fprintf(out_fp, "%s", sol_str);
                flint_free(sol_str);
                total_printed++;
            }
        }
        for (slong set = 0; set < sols->num_real_solution_sets; set++) {
            slong num_sols = sols->real_solutions_per_var[set * sols->num_variables + var];
            for (slong sol = 0; sol < num_sols; sol++) {
                if (total_printed > 0) fprintf(out_fp, ", ");
                dixon_fprint_arb_pretty(out_fp, sols->real_solution_sets[set][var][sol], 20);
                total_printed++;
            }
        }
        fprintf(out_fp, "}");
        if (total_printed > 1)      fprintf(out_fp, " (%ld solutions)", total_printed);
        else if (total_printed == 0) fprintf(out_fp, " (no solutions)");
        fprintf(out_fp, "\n");
    }
    if (sols->real_root_summary) {
        fprintf(out_fp, "\n%s\n", sols->real_root_summary);
    }
    fprintf(out_fp, "=== Solution Complete ===\n\n");
    fclose(out_fp);
}

static int redirect_fd_to_devnull(int target_fd, int *orig_fd)
{
    int devnull;

    if (!orig_fd) {
        return 0;
    }

    fflush(stdout);
    fflush(stderr);

    *orig_fd = dup(target_fd);
    if (*orig_fd == -1) {
        return 0;
    }

    devnull = open(DIXON_NULL_DEVICE, O_WRONLY);
    if (devnull == -1) {
        close(*orig_fd);
        *orig_fd = -1;
        return 0;
    }

    if (dup2(devnull, target_fd) == -1) {
        close(devnull);
        close(*orig_fd);
        *orig_fd = -1;
        return 0;
    }

    close(devnull);
    return 1;
}

static void restore_fd(int target_fd, int orig_fd);

static int redirect_stdio_to_devnull(int *orig_stdout, int *orig_stderr)
{
    if (!orig_stdout || !orig_stderr) {
        return 0;
    }

    *orig_stdout = -1;
    *orig_stderr = -1;

    if (!redirect_fd_to_devnull(STDOUT_FILENO, orig_stdout)) {
        return 0;
    }

    if (!redirect_fd_to_devnull(STDERR_FILENO, orig_stderr)) {
        restore_fd(STDOUT_FILENO, *orig_stdout);
        *orig_stdout = -1;
        return 0;
    }

    return 1;
}

static void restore_fd(int target_fd, int orig_fd)
{
    fflush(stdout);
    fflush(stderr);

    if (orig_fd != -1) {
        dup2(orig_fd, target_fd);
        close(orig_fd);
    }
}

static void restore_stdio(int orig_stdout, int orig_stderr)
{
    restore_fd(STDOUT_FILENO, orig_stdout);
    restore_fd(STDERR_FILENO, orig_stderr);
}

static void print_solution_set_brief(FILE *out_fp,
                                     const polynomial_solutions_t *sols,
                                     slong set_idx)
{
    for (slong var = 0; var < sols->num_variables; var++) {
        slong num_sols = sols->solutions_per_var[set_idx * sols->num_variables + var];

        if (var > 0) {
            fprintf(out_fp, ", ");
        }

        fprintf(out_fp, "%s=", sols->variable_names[var]);
        if (num_sols == 0) {
            fprintf(out_fp, "no solution");
        } else if (num_sols == 1) {
            char *sol_str = fq_nmod_get_str_pretty(sols->solution_sets[set_idx][var][0], sols->ctx);
            if (sol_str) {
                fprintf(out_fp, "%s", sol_str);
                free(sol_str);
            }
        } else {
            fprintf(out_fp, "{");
            for (slong sol = 0; sol < num_sols; sol++) {
                char *sol_str = fq_nmod_get_str_pretty(sols->solution_sets[set_idx][var][sol], sols->ctx);
                if (sol > 0) {
                    fprintf(out_fp, ", ");
                }
                if (sol_str) {
                    fprintf(out_fp, "%s", sol_str);
                    free(sol_str);
                }
            }
            fprintf(out_fp, "}");
        }
    }
}

static void print_polynomial_solutions_brief(const polynomial_solutions_t *sols)
{
    if (!sols) {
        printf("Status: solve failed (no solution data available)\n");
        return;
    }

    if (!sols->is_valid) {
        printf("Status: solve failed");
        if (sols->error_message) {
            printf(" (%s)", sols->error_message);
        }
        printf("\n");
        return;
    }

    if (sols->has_no_solutions == -1) {
        printf("Status: positive-dimensional system\n");
        if (sols->elimination_summary) {
            printf("Elimination: %s\n", sols->elimination_summary);
        }
        if (sols->total_combinations > 0) {
            printf("Resultants: %ld/%ld non-zero combination(s)\n",
                   sols->successful_combinations, sols->total_combinations);
        }
        return;
    }

    if (sols->has_no_solutions == 1) {
        printf("Status: no finite-field solutions\n");
        if (sols->elimination_summary) {
            printf("Elimination: %s\n", sols->elimination_summary);
        }
        if (sols->total_combinations > 0) {
            printf("Resultants: %ld/%ld non-zero combination(s)\n",
                   sols->successful_combinations, sols->total_combinations);
        }
        if (sols->num_base_solutions > 0) {
            printf("Base roots: %ld\n", sols->num_base_solutions);
        }
        return;
    }

    printf("Status: found %ld solution set(s)\n", sols->num_solution_sets);
    if (sols->variable_order) {
        printf("Variables: %s\n", sols->variable_order);
    }
    if (sols->elimination_summary) {
        printf("Elimination: %s\n", sols->elimination_summary);
    }
    if (sols->total_combinations > 0) {
        printf("Resultants: %ld/%ld non-zero combination(s)\n",
               sols->successful_combinations, sols->total_combinations);
    }
    if (sols->num_base_solutions > 0) {
        printf("Base roots: %ld\n", sols->num_base_solutions);
    }
    if (sols->checked_solution_sets >= 0 && sols->verified_solution_sets >= 0) {
        printf("Verification: %ld/%ld candidate set(s) passed\n",
               sols->verified_solution_sets, sols->checked_solution_sets);
    }

    printf("Solutions:\n");
    for (slong set = 0; set < sols->num_solution_sets; set++) {
        printf("  [%ld] ", set + 1);
        print_solution_set_brief(stdout, sols, set);
        printf("\n");
    }
}

static void save_result_to_file(const char *filename,
                                const char *polys_str,
                                const char *vars_str,
                                const char *ideal_str,
                                const char *allvars_str,
                                const fmpz_t prime, ulong power,
                                const char *result,
                                double cpu_time, double wall_time, int threads_num)
{
    FILE *out_fp = fopen(filename, "w");
    if (!out_fp) {
        fprintf(stderr, "Warning: Could not create output file '%s'\n", filename);
        return;
    }

    fprintf(out_fp, "Dixon Resultant Computation\n");
    fprintf(out_fp, "==========================\n");
    fprintf(out_fp, "Field: ");
    print_field_label(out_fp, prime, power);
    if (!fmpz_is_zero(prime) && power > 1) {
        fprintf(out_fp, "\nField extension generator: t");
    }
    fprintf(out_fp, "\n");

    if (ideal_str && allvars_str) {
        fprintf(out_fp, "Mode: Dixon with ideal reduction\n");
        fprintf(out_fp, "Ideal generators: %s\n", ideal_str);
        fprintf(out_fp, "All variables: %s\n", allvars_str);
    } else {
        fprintf(out_fp, "Mode: Basic Dixon resultant\n");
    }
    fprintf(out_fp, "Variables eliminated: %s\n", vars_str);
    fprintf(out_fp, "Polynomials: %s\n", polys_str);
    fprintf(out_fp, "CPU time: %.3f seconds | Wall time: %.3f seconds | Threads: %d\n", cpu_time, wall_time, threads_num);
    fprintf(out_fp, "\nResultant:\n%s\n", result);
    fclose(out_fp);
}

static int count_comma_separated_items(const char *str)
{
    if (!str || strlen(str) == 0) return 0;
    int count = 1;
    for (const char *p = str; *p; p++)
        if (*p == ',') count++;
    return count;
}

static int parse_x_index_alias(const char *name, slong *idx_out)
{
    char *endptr;
    long idx;

    if (!name || name[0] != 'x' || !isdigit((unsigned char) name[1]))
        return 0;

    idx = strtol(name + 1, &endptr, 10);
    if (*endptr != '\0' || idx < 0)
        return 0;

    *idx_out = (slong) idx;
    return 1;
}

static char *normalize_elimination_vars(const char *polys_str,
                                        const char *vars_str,
                                        int silent_mode)
{
    slong num_polys, num_vars, num_all_vars;
    char **poly_arr = NULL, **vars_arr = NULL, **all_vars = NULL;
    char **mapped = NULL;
    char *result = NULL;
    int changed = 0;

    if (!polys_str || !vars_str) return NULL;

    poly_arr = split_string(polys_str, &num_polys);
    vars_arr = split_string(vars_str, &num_vars);
    collect_variables((const char **) poly_arr, num_polys, NULL, &all_vars, &num_all_vars);

    mapped = (char **) malloc((size_t) num_vars * sizeof(char *));
    for (slong i = 0; i < num_vars; i++) {
        slong alias_idx;
        int found = 0;

        for (slong j = 0; j < num_all_vars; j++) {
            if (strcmp(vars_arr[i], all_vars[j]) == 0) {
                mapped[i] = strdup(vars_arr[i]);
                found = 1;
                break;
            }
        }

        if (!found && parse_x_index_alias(vars_arr[i], &alias_idx) &&
            alias_idx >= 0 && alias_idx < num_all_vars) {
            mapped[i] = strdup(all_vars[alias_idx]);
            if (strcmp(mapped[i], vars_arr[i]) != 0) changed = 1;
            found = 1;
        }

        if (!found) {
            mapped[i] = strdup(vars_arr[i]);
        }
    }

    if (changed) {
        size_t total_len = 1;
        for (slong i = 0; i < num_vars; i++) total_len += strlen(mapped[i]) + 2;
        result = (char *) malloc(total_len);
        result[0] = '\0';
        for (slong i = 0; i < num_vars; i++) {
            if (i > 0) strcat(result, ",");
            strcat(result, mapped[i]);
        }
        if (!silent_mode)
            printf("Normalized elimination variables: %s -> %s\n", vars_str, result);
    }

    for (slong i = 0; i < num_polys; i++) free(poly_arr[i]);
    free(poly_arr);
    for (slong i = 0; i < num_vars; i++) free(vars_arr[i]);
    free(vars_arr);
    for (slong i = 0; i < num_all_vars; i++) free(all_vars[i]);
    free(all_vars);
    for (slong i = 0; i < num_vars; i++) free(mapped[i]);
    free(mapped);

    return result;
}

static int generate_random_poly_strings_rational(
        const long *degrees, slong npolys,
        int is_solve_mode,
        int silent_mode,
        char **polys_str_out,
        char **elim_vars_str_out,
        char **all_vars_str_out)
{
    slong nvars = npolys;
    char **var_names = malloc((size_t) nvars * sizeof(char *));
    char coeff_chars[] = { -2, -1, 0, 1, 2 };
    char *all_vars = NULL;
    char *elim_vars = NULL;
    char *polys_str = NULL;
    size_t all_cap, elim_cap, polys_cap = 256;

    srand((unsigned) (time(NULL) ^ clock()));

    for (slong i = 0; i < nvars; i++) {
        var_names[i] = malloc(16);
        snprintf(var_names[i], 16, "x%ld", i);
    }

    all_cap = (size_t) nvars * 8 + 2;
    all_vars = malloc(all_cap);
    all_vars[0] = '\0';
    for (slong i = 0; i < nvars; i++) {
        if (i > 0) strcat(all_vars, ",");
        strcat(all_vars, var_names[i]);
    }

    if (!is_solve_mode && nvars >= 2) {
        elim_cap = (size_t) (nvars - 1) * 8 + 2;
        elim_vars = malloc(elim_cap);
        elim_vars[0] = '\0';
        for (slong i = 0; i < nvars - 1; i++) {
            if (i > 0) strcat(elim_vars, ",");
            strcat(elim_vars, var_names[i]);
        }
    } else {
        elim_vars = strdup(all_vars);
    }

    polys_str = malloc(polys_cap);
    polys_str[0] = '\0';

    if (!silent_mode) printf("Generating random polynomial system over Q...\n");

    fflush(stdout);
    for (slong i = 0; i < npolys; i++) {
        char poly_buf[8192];
        int first_term = 1;
        slong degree = degrees[i] > 0 ? degrees[i] : 1;
        slong term_count = FLINT_MAX(3, degree + 2);

        poly_buf[0] = '\0';
        for (slong t = 0; t < term_count; t++) {
            int coeff = 0;
            slong *exp = calloc((size_t) nvars, sizeof(slong));
            slong total_deg;
            char term_buf[512];
            int has_var = 0;

            while (coeff == 0) {
                coeff = coeff_chars[rand() % 5];
            }

            if (t == 0) {
                exp[i % nvars] = degree;
            } else {
                total_deg = rand() % (degree + 1);
                for (slong left = total_deg; left > 0; left--) {
                    exp[rand() % nvars]++;
                }
            }

            term_buf[0] = '\0';
            if (first_term) {
                if (coeff < 0) strcat(term_buf, "-");
            } else {
                strcat(term_buf, coeff < 0 ? " - " : " + ");
            }

            if (abs(coeff) != 1) {
                char coeff_buf[32];
                snprintf(coeff_buf, sizeof(coeff_buf), "%d", abs(coeff));
                strcat(term_buf, coeff_buf);
            }

            for (slong j = 0; j < nvars; j++) {
                if (exp[j] == 0) continue;
                if (abs(coeff) != 1 || has_var) strcat(term_buf, "*");
                strcat(term_buf, var_names[j]);
                if (exp[j] != 1) {
                    char exp_buf[32];
                    snprintf(exp_buf, sizeof(exp_buf), "^%ld", exp[j]);
                    strcat(term_buf, exp_buf);
                }
                has_var = 1;
            }

            if (!has_var && abs(coeff) == 1) strcat(term_buf, "1");
            strcat(poly_buf, term_buf);
            first_term = 0;
            free(exp);
        }

        if (strlen(polys_str) + strlen(poly_buf) + 4 > polys_cap) {
            while (strlen(polys_str) + strlen(poly_buf) + 4 > polys_cap)
                polys_cap *= 2;
            polys_str = realloc(polys_str, polys_cap);
        }
        if (i > 0) strcat(polys_str, ", ");
        strcat(polys_str, poly_buf);
    }

    if (!silent_mode) {
        printf("System: %ld equations, %ld variables\n", npolys, nvars);
        printf("Degrees: [");
        for (slong i = 0; i < npolys; i++) {
            if (i > 0) printf(", ");
            printf("%ld", degrees[i]);
        }
        printf("]\n");
        if (!is_solve_mode)
            printf("Eliminate: %s  |  Remaining: %s\n",
                   elim_vars, var_names[nvars - 1]);
    }

    for (slong i = 0; i < nvars; i++) free(var_names[i]);
    free(var_names);

    *polys_str_out = polys_str;
    *elim_vars_str_out = elim_vars;
    *all_vars_str_out = all_vars;
    return 1;
}

/* =========================================================================
 * main
 * ========================================================================= */
int main(int argc, char *argv[])
{
    clock_t start_time = clock();
    const char *prog_name = display_prog_name(argv[0]);
    
    // Initialize wall clock time at the very beginning
    struct timeval program_start;
    gettimeofday(&program_start, NULL);

    if (argc == 1) { print_version(); print_usage(prog_name); return 0; }

    /* ---- parse leading flags ---- */
    int    silent_mode = 0;
    int    solve_mode  = 0;
    int    solve_verbose_mode = 0;
    int    comp_mode   = 0;
    int    rand_mode   = 0;   /* --random / -r */
    int    ideal_mode  = 0;   /*  --ideal flag */
    int    field_eq_mode = 0; /* --field-equation */
    int    arg_offset  = 0;
    double omega       = DIXON_OMEGA;   /* default, overridden by --omega */
    int    det_method_step1 = -1;  /* determinant method override for step 1 */
    int    det_method_step4 = -1;  /* determinant method override for step 4 */
    int    num_threads = -1;  /* number of threads, -1 means use default */

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--silent") == 0) {
            silent_mode = 1; arg_offset++;
        } else if (strcmp(argv[i], "--solve-verbose") == 0) {
            solve_mode = 1; solve_verbose_mode = 1; arg_offset++;
        } else if (strcmp(argv[i], "--solve") == 0) {
            solve_mode = 1; arg_offset++;
        } else if (strcmp(argv[i], "--comp") == 0 ||
                   strcmp(argv[i], "-c")     == 0) {
            comp_mode = 1; arg_offset++;
        } else if (strcmp(argv[i], "--random") == 0 ||
                   strcmp(argv[i], "-r")       == 0) {
            rand_mode = 1; arg_offset++;
        } else if (strcmp(argv[i], "--ideal") == 0) {  /* <<< NEW >>> */
            ideal_mode = 1; arg_offset++;
        } else if (strcmp(argv[i], "--field-equation") == 0 ||
                   strcmp(argv[i], "--field-eqution")  == 0) {
            field_eq_mode = 1; arg_offset++;
        } else if ((strcmp(argv[i], "--omega") == 0 ||
                    strcmp(argv[i], "-w")      == 0) && i + 1 < argc) {
            char *endptr = NULL;
            double val = strtod(argv[i + 1], &endptr);
            if (endptr && *endptr == '\0' && val > 0.0) {
                omega = val;
            } else {
                fprintf(stderr, "Warning: invalid --omega value '%s', "
                                "using default %.4g\n", argv[i + 1], omega);
            }
            arg_offset += 2;
            i++;          /* skip the value token */
        } else if ((strcmp(argv[i], "--method") == 0) && i + 1 < argc) {
            char *endptr = NULL;
            long val = strtol(argv[i + 1], &endptr, 10);
            if (endptr && *endptr == '\0' && val >= 0 && val <= 3) {
                det_method_step1 = (int)val;
                det_method_step4 = (int)val;
            } else {
                fprintf(stderr, "Warning: invalid --method value '%s', "
                                "must be 0-3. Using default.\n", argv[i + 1]);
            }
            arg_offset += 2;
            i++;          /* skip the value token */
        } else if ((strcmp(argv[i], "--step1") == 0) && i + 1 < argc) {
            char *endptr = NULL;
            long val = strtol(argv[i + 1], &endptr, 10);
            if (endptr && *endptr == '\0' && val >= 0 && val <= 3) {
                det_method_step1 = (int)val;
            } else {
                fprintf(stderr, "Warning: invalid --step1 value '%s', "
                                "must be 0-3. Using default.\n", argv[i + 1]);
            }
            arg_offset += 2;
            i++;          /* skip the value token */
        } else if ((strcmp(argv[i], "--step4") == 0) && i + 1 < argc) {
            char *endptr = NULL;
            long val = strtol(argv[i + 1], &endptr, 10);
            if (endptr && *endptr == '\0' && val >= 0 && val <= 3) {
                det_method_step4 = (int)val;
            } else {
                fprintf(stderr, "Warning: invalid --step4 value '%s', "
                                "must be 0-3. Using default.\n", argv[i + 1]);
            }
            arg_offset += 2;
            i++;          /* skip the value token */
        } else if ((strcmp(argv[i], "--threads") == 0) && i + 1 < argc) {
            char *endptr = NULL;
            long val = strtol(argv[i + 1], &endptr, 10);
            if (endptr && *endptr == '\0' && val > 0) {
                num_threads = (int)val;
            } else {
                fprintf(stderr, "Warning: invalid --threads value '%s', "
                                "must be positive integer. Using default.\n", argv[i + 1]);
            }
            arg_offset += 2;
            i++;          /* skip the value token */
        } else {
            break;
        }
    }

    int    effective_argc = argc - arg_offset;
    char **effective_argv = argv + arg_offset;


    
    /* ---- help ---- */
    if (effective_argc >= 2 &&
        (strcmp(effective_argv[1], "--help") == 0 ||
         strcmp(effective_argv[1], "-h")     == 0)) {
        if (!silent_mode)  { print_version(); print_usage(prog_name); }
        return 0;
    }

    /* ---- version banner ---- */
    if (!silent_mode) {
        print_version();
    }

    /* ---- test modes ---- */
    if (argc >= 2 && strcmp(argv[1], "--test") == 0) {
        if (!silent_mode) {
            if (argc >= 3) test_dixon(atoi(argv[2]));
            else           test_dixon(0);
        }
        return 0;
    }
    if (effective_argc >= 2 &&
        strcmp(effective_argv[1], "--test-solver") == 0) {
        if (!silent_mode) test_polynomial_solver();
        return 0;
    }

    /* ---- input variables ---- */
    char *polys_str     = NULL;
    char *vars_str      = NULL;
    char *ideal_str     = NULL;
    char *allvars_str   = NULL;
    char *field_str     = NULL;
    char *vars_str_override = NULL;
    int   need_free     = 0;
    int   rand_generated = 0;  /* 1 = polys_str/vars_str were malloc'd by random gen */
    char *deg_str       = NULL; /* degree list string when rand_mode */
    char *input_filename  = NULL;
    char *output_filename = NULL;
    mp_limb_t prime = 0;
    ulong     power = 0;

    /* ---- determine input mode ---- */
    if (rand_mode) {
        /* --random / -r: positional args are "d1,d2,...,dn" field_size */
        if (effective_argc != 3) {
            if (!silent_mode) {
                fprintf(stderr, "Error: Random mode requires exactly:\n");
                fprintf(stderr, "  %s [--solve|--comp] --random \"d1,d2,...,dn\" field_size\n",
                        prog_name);
            }
            return 1;
        }
        deg_str   = effective_argv[1];
        field_str = effective_argv[2];

        output_filename = generate_timestamped_filename(comp_mode ? "comp" : "solution");
    /*  --ideal mode */
    } else if (ideal_mode) {
        if (effective_argc == 2) {
            /* File input: auto-detect ideal generators by '=' */
            FILE *fp = fopen(effective_argv[1], "r");
            if (!fp) {
                if (!silent_mode)
                    fprintf(stderr, "Error: Cannot open file '%s'\n", effective_argv[1]);
                return 1;
            }
            if (!silent_mode)
                printf("Reading from file (--ideal mode): %s\n", effective_argv[1]);

            input_filename  = strdup(effective_argv[1]);
            output_filename = generate_tagged_filename(input_filename, "_solution");

            if (!read_ideal_file(fp, &field_str, &polys_str,
                                 &vars_str, &ideal_str)) {
                fclose(fp); return 1;
            }
            fclose(fp);
            need_free = 1;

            if (!ideal_str) {
                if (!silent_mode)
                    fprintf(stderr, "Warning: No ideal generators found in file "
                                    "(no lines with '='). Running basic Dixon.\n");
            }

        } else if (effective_argc == 5) {
            /* CLI: --ideal "generators" "polys" "elim_vars" field_size */
            ideal_str = effective_argv[1];
            polys_str = effective_argv[2];
            vars_str  = effective_argv[3];
            field_str = effective_argv[4];

            output_filename = generate_timestamped_filename("solution");

        } else {
            if (!silent_mode) {
                fprintf(stderr, "Error: --ideal mode requires either:\n");
                fprintf(stderr, "  %s --ideal \"ideal_generators\" \"polynomials\" \"eliminate_vars\" field_size\n",
                        prog_name);
                fprintf(stderr, "  %s --ideal input_file\n", prog_name);
            }
            return 1;
        }
    } else if (comp_mode) {
        /* Complexity analysis mode: same argument format as basic Dixon */
        if (effective_argc == 2) {
            /* file input */
            FILE *fp = fopen(effective_argv[1], "r");
            if (!fp) {
                if (!silent_mode)
                    fprintf(stderr, "Error: Cannot open file '%s'\n",
                            effective_argv[1]);
                return 1;
            }
            if (!silent_mode)
                printf("Reading from file: %s\n", effective_argv[1]);

            input_filename  = strdup(effective_argv[1]);
            output_filename = generate_tagged_filename(input_filename, "_comp");

            if (!read_multiline_file(fp, &field_str, &polys_str,
                                     &vars_str, &ideal_str, &allvars_str)) {
                fclose(fp); return 1;
            }
            fclose(fp);
            need_free = 1;

        } else if (effective_argc == 4) {
            /* command line: polys vars field_size */
            polys_str = effective_argv[1];
            vars_str  = effective_argv[2];
            field_str = effective_argv[3];

            output_filename = generate_timestamped_filename("comp");

        } else {
            if (!silent_mode) {
                fprintf(stderr, "Error: Complexity mode requires:\n");
                fprintf(stderr, "  %s --comp \"polynomials\" \"eliminate_vars\" field_size\n",
                        prog_name);
                fprintf(stderr, "  %s --comp input_file\n", prog_name);
            }
            return 1;
        }

    } else if (solve_mode) {
        if (effective_argc == 2) {
            FILE *fp = fopen(effective_argv[1], "r");
            if (!fp) {
                if (!silent_mode)
                    fprintf(stderr, "Error: Cannot open file '%s'\n",
                            effective_argv[1]);
                return 1;
            }
            if (!silent_mode)
                printf("Reading polynomial system from file: %s\n",
                       effective_argv[1]);

            input_filename  = strdup(effective_argv[1]);
            output_filename = generate_tagged_filename(input_filename, "_solution");

            if (!read_solver_file(fp, &field_str, &polys_str)) {
                fclose(fp); return 1;
            }
            fclose(fp);
            need_free = 1;

        } else if (effective_argc == 3) {
            polys_str = effective_argv[1];
            field_str = effective_argv[2];
            output_filename = generate_timestamped_filename("solution");

        } else {
            if (!silent_mode) {
                fprintf(stderr, "Error: Solver mode requires either:\n");
                fprintf(stderr, "  %s --solve \"polynomials\" field_size\n", prog_name);
                fprintf(stderr, "  %s --solve input_file\n", prog_name);
            }
            return 1;
        }

    } else {
        /* Basic Dixon, or auto-detect solve/file mode when flags are omitted */
        if (effective_argc == 3) {
            /* Auto-detect solve mode: only polynomials and field size */
            solve_mode = 1;
            polys_str = effective_argv[1];
            field_str = effective_argv[2];
            output_filename = generate_timestamped_filename("solution");
        } else if (effective_argc == 2) {
            FILE *fp = fopen(effective_argv[1], "r");
            int detected_solve_mode = 0;
            if (!fp) {
                if (!silent_mode)
                    fprintf(stderr, "Error: Cannot open file '%s'\n",
                            effective_argv[1]);
                return 1;
            }
            if (!silent_mode)
                printf("Reading from file: %s\n", effective_argv[1]);

            input_filename  = strdup(effective_argv[1]);
            output_filename = generate_tagged_filename(input_filename, "_solution");

            if (!read_auto_input_file(fp, &detected_solve_mode, &field_str, &polys_str,
                                      &vars_str, &ideal_str, &allvars_str)) {
                fclose(fp); return 1;
            }
            fclose(fp);
            solve_mode = detected_solve_mode;
            need_free = 1;

            if (!silent_mode)
                printf("Detected file mode: %s\n", solve_mode ? "solver" : "elimination");

        } else if (effective_argc == 4 || effective_argc == 5) {
            polys_str = effective_argv[1];
            vars_str  = effective_argv[2];
            if (effective_argc == 5) {
                ideal_str = effective_argv[3];
                field_str = effective_argv[4];
            } else {
                field_str = effective_argv[3];
            }
            output_filename = generate_timestamped_filename("solution");

        } else {
            if (!silent_mode) print_usage(prog_name);
            return 1;
        }
    }

    /* ---- parse field size ---- */
    fmpz_t p_fmpz;
    fmpz_init(p_fmpz);
    char *field_poly_str = NULL;
    char *gen_var_name   = NULL;

    if (!parse_field_size(field_str, p_fmpz, &power,
                          &field_poly_str, &gen_var_name)) {
        if (!silent_mode) {
            fprintf(stderr, "Error: Invalid field size '%s'\n", field_str);
            fprintf(stderr, "Field size must be 0 (for Q), a prime, prime power (e.g. 256), or p^k (e.g. 2^8)\n");
        }
        if (need_free) {
            free(field_str); free(polys_str);
            if (vars_str)    free(vars_str);
            if (ideal_str)   free(ideal_str);
        }
        if (input_filename)  free(input_filename);
        if (output_filename) free(output_filename);
        fmpz_clear(p_fmpz);
        if (field_poly_str) free(field_poly_str);
        if (gen_var_name)   free(gen_var_name);
        return 1;
    }

    int rational_mode = fmpz_is_zero(p_fmpz);
    int large_prime_mode = (!rational_mode && power == 1 && !fmpz_abs_fits_ui(p_fmpz));
    int ctx_initialized = 0;

    if (!rational_mode && !large_prime_mode) {
        prime = fmpz_get_ui(p_fmpz);
    }

    if (!silent_mode) {
        if (!comp_mode && !solve_mode) {
            printf("=== Dixon Resultant Computation ===\n");
            printf("Field: ");
            print_field_label(stdout, p_fmpz, power);
            printf("\n");
            if (!rational_mode && power > 1) {
                printf("Field extension generator: t\n");
            }
        } else if (comp_mode) {
            printf("Mode: Complexity analysis  |  Field: ");
            print_field_label(stdout, p_fmpz, power);
            printf("\n");
        } else {
            printf("Mode: Polynomial system solver  |  Field: ");
            print_field_label(stdout, p_fmpz, power);
            printf("\n");
        }
    }

    if (rational_mode) {
        if (ideal_str) {
            fprintf(stderr, "Error: field_size=0 currently does not support --ideal.\n");
            goto cleanup_fail;
        }
        if (field_eq_mode) {
            fprintf(stderr, "Error: field_size=0 currently does not support --field-equation.\n");
            goto cleanup_fail;
        }
    }
    if (large_prime_mode) {
        if (solve_mode) {
            fprintf(stderr, "Error: large prime fields beyond the nmod limit currently support Dixon resultant only; solver mode is not implemented.\n");
            goto cleanup_fail;
        }
        if (ideal_str) {
            fprintf(stderr, "Error: large prime fallback currently does not support --ideal.\n");
            goto cleanup_fail;
        }
        if (field_eq_mode) {
            fprintf(stderr, "Error: large prime fallback currently does not support --field-equation.\n");
            goto cleanup_fail;
        }
        if (rand_mode) {
            fprintf(stderr, "Error: large prime fallback currently does not support --random.\n");
            goto cleanup_fail;
        }
        if (comp_mode) {
            fprintf(stderr, "Error: large prime fallback currently does not support --comp.\n");
            goto cleanup_fail;
        }
        if (!silent_mode) {
            printf("Prime exceeds nmod limit; using Q reconstruction fallback before reducing modulo the target prime.\n");
        }
    }
    if (!rational_mode && power > 1 && !fmpz_abs_fits_ui(p_fmpz)) {
        fprintf(stderr, "Error: extension fields with characteristic beyond the nmod limit are not supported.\n");
        goto cleanup_fail;
    }

    /* ---- initialize finite field ---- */
    fq_nmod_ctx_t ctx;

    if (!rational_mode && !large_prime_mode && power > 1 && field_poly_str) {
        const char *var_name = gen_var_name ? gen_var_name : "t";
        if (!silent_mode) {
            printf("Using custom field polynomial: %s\n", field_poly_str);
            printf("Generator variable: %s\n", var_name);
        }
        nmod_poly_t modulus;
        nmod_poly_init(modulus, prime);
        if (!parse_field_polynomial(modulus, field_poly_str, prime, var_name)) {
            fprintf(stderr, "Error: Failed to parse field polynomial\n");
            nmod_poly_clear(modulus); fmpz_clear(p_fmpz);
            if (field_poly_str) free(field_poly_str);
            if (gen_var_name)   free(gen_var_name);
            return 1;
        }
        if (nmod_poly_degree(modulus) != (slong)power) {
            fprintf(stderr, "Error: Polynomial degree %ld doesn't match power %lu\n",
                    nmod_poly_degree(modulus), power);
            nmod_poly_clear(modulus); fmpz_clear(p_fmpz);
            if (field_poly_str) free(field_poly_str);
            if (gen_var_name)   free(gen_var_name);
            return 1;
        }
        fq_nmod_ctx_init_modulus(ctx, modulus, var_name);
        ctx_initialized = 1;
        nmod_poly_clear(modulus);

    } else if (!rational_mode && !large_prime_mode) {
        fq_nmod_ctx_init(ctx, p_fmpz, power, "t");
        ctx_initialized = 1;

        if (!silent_mode && power > 1) {
            printf("Using FLINT's default irreducible polynomial:\n  ");
            const nmod_poly_struct *modulus = ctx->modulus;
            int first_term = 1;
            for (slong i = nmod_poly_degree(modulus); i >= 0; i--) {
                mp_limb_t coeff = nmod_poly_get_coeff_ui(modulus, i);
                if (coeff != 0) {
                    if (!first_term) printf(" + ");
                    if      (i == 0) printf("%lu", coeff);
                    else if (i == 1) { if (coeff == 1) printf("t"); else printf("%lu*t", coeff); }
                    else             { if (coeff == 1) printf("t^%ld", i); else printf("%lu*t^%ld", coeff, i); }
                    first_term = 0;
                }
            }
            printf("\n");
        }
    }

    /* ---- activate field-equation reduction mode ---- */
    if (!rational_mode && field_eq_mode) {
        fq_mvpoly_set_field_equation_reduction(1);
        if (!silent_mode)
            printf("Reduction: field equations enabled\n");
    }

    /* ======================================================
     * If --random/-r, generate polynomial strings now that ctx is ready
     * ====================================================== */
    if (rand_mode) {
        slong  npolys_rand;
        long  *degrees_rand = parse_degree_list(deg_str, &npolys_rand);

        if (npolys_rand < 2) {
            if (!silent_mode)
                fprintf(stderr, "Error: --random requires at least 2 degrees "
                                "(e.g. \"3,3,2\")\n");
            free(degrees_rand);
            fq_nmod_ctx_clear(ctx);
            fmpz_clear(p_fmpz);
            if (field_poly_str) free(field_poly_str);
            if (gen_var_name)   free(gen_var_name);
            if (output_filename) free(output_filename);
            return 1;
        }

        char *gen_polys = NULL, *gen_elim = NULL, *gen_allvars = NULL;
        int generated_ok;
        if (rational_mode) {
            generated_ok = generate_random_poly_strings_rational(
                               degrees_rand, npolys_rand,
                               solve_mode, silent_mode,
                               &gen_polys, &gen_elim, &gen_allvars);
        } else {
            generated_ok = generate_random_poly_strings(
                               degrees_rand, npolys_rand,
                               ctx,
                               solve_mode,
                               silent_mode,
                               &gen_polys, &gen_elim, &gen_allvars);
        }
        if (!generated_ok) {
            if (!silent_mode)
                fprintf(stderr, "Error: random polynomial generation failed\n");
            free(degrees_rand);
            if (!rational_mode) fq_nmod_ctx_clear(ctx);
            fmpz_clear(p_fmpz);
            if (field_poly_str) free(field_poly_str);
            if (gen_var_name)   free(gen_var_name);
            if (output_filename) free(output_filename);
            return 1;
        }

        free(degrees_rand);
        polys_str    = gen_polys;
        vars_str     = gen_elim;
        allvars_str  = gen_allvars;
        rand_generated = 1;
    }

    if (polys_str && vars_str) {
        char *normalized_vars = normalize_elimination_vars(polys_str, vars_str, silent_mode);
        if (normalized_vars) {
            if (need_free || rand_generated) {
                free(vars_str);
            } else {
                vars_str_override = normalized_vars;
            }
            vars_str = normalized_vars;
        }
    }

    /* ======================================================
     * Set determinant method and thread count
     * ====================================================== */
    dixon_global_method_step1 = -1;
    dixon_global_method_step4 = -1;
    dixon_global_method = -1;
    if (det_method_step1 != -1) {
        dixon_global_method_step1 = (det_method_t)det_method_step1;
        dixon_global_method = dixon_global_method_step1;
        if (!silent_mode) {
            printf("Step 1 determinant method: %s\n",
                   det_method_name_cli(det_method_step1));
        }
    }
    if (det_method_step4 != -1) {
        dixon_global_method_step4 = (det_method_t)det_method_step4;
        dixon_global_method = dixon_global_method_step4;
        if (!silent_mode) {
            printf("Step 4 determinant method: %s\n",
                   det_method_name_cli(det_method_step4));
        }
    }
    if (num_threads != -1) {
        fq_interpolation_set_threads(num_threads);
        // fq_nmod_poly_mat_det_set_threads(num_threads);
        if (!silent_mode) {
            printf("Using %d threads\n", num_threads);
        }
    }

    /* ======================================================
     * EXECUTE requested mode
     * ====================================================== */

    char *result = NULL;
    polynomial_solutions_t *solutions = NULL;
    rational_solutions_t *rational_solutions = NULL;

    if (comp_mode) {
        /* ---- Complexity analysis ---- */
        if (!silent_mode) {
            int var_count  = count_comma_separated_items(vars_str);
            printf("Eliminate   (%d): %s\n", var_count,  vars_str);
            printf("Omega            : %.4g\n", omega);
            printf("--------------------------------\n");
        }

        clock_t end_time  = clock();
        double comp_time  = (double)(end_time - start_time) / CLOCKS_PER_SEC;

        run_complexity_analysis(polys_str, vars_str,
                                p_fmpz, power, ctx_initialized ? ctx : NULL,
                                output_filename, silent_mode,
                                comp_time, omega);

    } else if (solve_mode) {
        /* ---- Polynomial system solver ---- */
        if (!silent_mode) {
            int poly_count = count_comma_separated_items(polys_str);
            printf("\nEquations: %d\n", poly_count);
            printf("--------------------------------\n");
        }

        int orig_stdout = -1, orig_stderr = -1;
        int suppress_solver_stdout = !silent_mode && !solve_verbose_mode;
        int suppress_solver_trace = silent_mode;

        if (rational_mode) {
            rational_solver_set_realtime_progress(0);

            if (suppress_solver_trace) {
                redirect_stdio_to_devnull(&orig_stdout, &orig_stderr);
            } else if (suppress_solver_stdout) {
                redirect_fd_to_devnull(STDOUT_FILENO, &orig_stdout);
            }

            rational_solutions = solve_rational_polynomial_system_string(polys_str);

            if (suppress_solver_trace) {
                restore_stdio(orig_stdout, orig_stderr);
            } else if (suppress_solver_stdout) {
                restore_fd(STDOUT_FILENO, orig_stdout);
            }
        } else {
            polynomial_solver_set_realtime_progress(0);

            if (suppress_solver_trace) {
                redirect_stdio_to_devnull(&orig_stdout, &orig_stderr);
            } else if (suppress_solver_stdout) {
                redirect_fd_to_devnull(STDOUT_FILENO, &orig_stdout);
            }

            solutions = solve_polynomial_system_string(polys_str, ctx);

            if (suppress_solver_trace) {
                restore_stdio(orig_stdout, orig_stderr);
            } else if (suppress_solver_stdout) {
                restore_fd(STDOUT_FILENO, orig_stdout);
            }
        }

    } else if (ideal_str) {
        /* ---- Dixon with ideal reduction ---- */
        if (!silent_mode) {
            printf("Task: Dixon + ideal reduction  |  Eliminate: %s\n", vars_str);
            printf("Ideal: %s\n", ideal_str);
            printf("--------------------------------\n");
        }
        result = dixon_with_ideal_reduction_str(polys_str, vars_str, ideal_str, ctx);

    } else {
        /* ---- Basic Dixon resultant ---- */
        if (!silent_mode) {
            int poly_count = count_comma_separated_items(polys_str);
            int var_count  = count_comma_separated_items(vars_str);
            printf("Task: Dixon resultant  |  Equations: %d  |  Eliminate: %s\n",
                   poly_count, vars_str);
            if (var_count != poly_count - 1)
                printf("WARNING: Dixon method requires eliminating exactly %d variables "
                       "for %d equations!\n", poly_count - 1, poly_count);
            printf("--------------------------------\n");
        }

        int orig_stdout = -1, orig_stderr = -1, devnull = -1;
        if (silent_mode) {
            fflush(stdout); fflush(stderr);
            orig_stdout = dup(STDOUT_FILENO);
            orig_stderr = dup(STDERR_FILENO);
            devnull = open(DIXON_NULL_DEVICE, O_WRONLY);
            if (devnull != -1) {
                dup2(devnull, STDOUT_FILENO);
                dup2(devnull, STDERR_FILENO);
                close(devnull);
            }
        }

        if (rational_mode) {
            if (output_filename) {
                result = dixon_str_rational_with_file(polys_str, vars_str, output_filename);
            } else {
                result = dixon_str_rational(polys_str, vars_str);
            }
        } else if (large_prime_mode) {
            result = dixon_str_large_prime(polys_str, vars_str, p_fmpz);
        } else {
            result = dixon_str(polys_str, vars_str, ctx);
        }

        if (silent_mode && orig_stdout != -1) {
            fflush(stdout); fflush(stderr);
            dup2(orig_stdout, STDOUT_FILENO);
            dup2(orig_stderr, STDERR_FILENO);
            close(orig_stdout); close(orig_stderr);
        }
    }

    /* ---- compute total time ---- */
    clock_t cpu_end_time   = clock();
    double  cpu_time       = (double)(cpu_end_time - start_time) / CLOCKS_PER_SEC;
    
    // Compute wall time
    struct timeval program_end;
    gettimeofday(&program_end, NULL);
    double wall_time = (program_end.tv_sec - program_start.tv_sec) + 
                      (program_end.tv_usec - program_start.tv_usec) / 1000000.0;

    /* ---- Get actual used threads ---- */
    int total_threads = 1;
    #ifdef _OPENMP
    if (g_interpolation_threads > 0) {
        total_threads = g_interpolation_threads;
    } else {
        total_threads = omp_get_max_threads();
        // Default to half threads if not set
        int half_threads = (total_threads + 1) / 2;
        total_threads = half_threads;
    }
    #endif

    /* ---- output results ---- */
    if (comp_mode) {
        /* already printed inside run_complexity_analysis */
        ;
    } else if (solve_mode) {
        if (rational_mode) {
            if (rational_solutions) {
                if (!silent_mode) {
                    print_rational_solutions(rational_solutions);
                }

                if (output_filename) {
                    save_rational_solver_result_to_file(output_filename, polys_str,
                                                       rational_solutions,
                                                       cpu_time, wall_time, total_threads);
                    if (!silent_mode)
                        printf("Result saved to: %s\n", output_filename);
                }
                rational_solutions_clear(rational_solutions);
                free(rational_solutions);
            } else {
                if (!silent_mode)
                    fprintf(stderr, "\nError: Rational polynomial system solving failed\n");
            }
        } else {
            if (solutions) {
                if (!silent_mode) {
                    print_polynomial_solutions(solutions);
                }

                if (output_filename) {
                    save_solver_result_to_file(output_filename, polys_str,
                                               p_fmpz, power, solutions,
                                               cpu_time, wall_time, total_threads);
                    if (!silent_mode)
                        printf("Result saved to: %s\n", output_filename);
                }
                polynomial_solutions_clear(solutions);
                free(solutions);
            } else {
                if (!silent_mode)
                    fprintf(stderr, "\nError: Polynomial system solving failed\n");
            }
        }
    } else {
        if (result) {
            if (output_filename && !rational_mode) {
                save_result_to_file(output_filename, polys_str, vars_str,
                                    ideal_str, allvars_str, p_fmpz, power,
                                    result, cpu_time, wall_time, total_threads);
                if (!silent_mode)
                    printf("\nResult saved to: %s\n", output_filename);
                
                FILE *fp_append = fopen(output_filename, "a");
                if (fp_append) {
                    append_roots_to_file_from_result(result, polys_str, vars_str, ctx, fp_append);
                    fclose(fp_append);
                }
            } else if (output_filename && rational_mode) {
                if (!silent_mode)
                    printf("\nResult saved to: %s\n", output_filename);
            }
            free(result);
        } else {
            if (!silent_mode)
                fprintf(stderr, "\nError: Computation failed\n");
        }
    }

    printf("Total - CPU time: %.3f seconds | Wall time: %.3f seconds | Threads: %d\n", 
           cpu_time, wall_time, total_threads);

    /* ---- cleanup ---- */
    if (ctx_initialized) fq_nmod_ctx_clear(ctx);
    if (field_poly_str) free(field_poly_str);
    if (gen_var_name)   free(gen_var_name);
    fmpz_clear(p_fmpz);

    if (need_free) {
        free(field_str);
        free(polys_str);
        if (vars_str)   free(vars_str);
        if (ideal_str)  free(ideal_str);
        if (allvars_str) free(allvars_str);
    }
    if (rand_generated) {
        free(polys_str);
        free(vars_str);
        free(allvars_str);
    }
    if (!need_free && !rand_generated && vars_str_override)
        free(vars_str_override);
    if (input_filename)  free(input_filename);
    if (output_filename) free(output_filename);

    cleanup_unified_workspace();
    flint_cleanup();
    return 0;

cleanup_fail:
    if (field_poly_str) free(field_poly_str);
    if (gen_var_name)   free(gen_var_name);
    fmpz_clear(p_fmpz);

    if (need_free) {
        free(field_str);
        free(polys_str);
        if (vars_str)    free(vars_str);
        if (ideal_str)   free(ideal_str);
        if (allvars_str) free(allvars_str);
    }
    if (!need_free && !rand_generated && vars_str_override)
        free(vars_str_override);
    if (input_filename)  free(input_filename);
    if (output_filename) free(output_filename);
    flint_cleanup();
    return 1;
}
