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
#include "large_prime_system_solver.h"
#include "polynomial_system_solver.h"
#include "rational_system_solver.h"
#include "dixon_test.h"

#define PROGRAM_VERSION "0.2.3"

#ifdef _WIN32
#define DIXON_NULL_DEVICE "NUL"
#else
#define DIXON_NULL_DEVICE "/dev/null"
#endif

/* =========================================================================
 * Print usage
 * ========================================================================= */

static void print_version()
{
    printf("===============================================\n");
    printf("DRSolve v%s\n", PROGRAM_VERSION);
    printf("FLINT version: %s (Recommended: 3.5.0)\n", FLINT_VERSION);
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
    printf("    %s -o output.dr \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    -> Default output file: solution_YYYYMMDD_HHMMSS.dr\n");

    printf("  Polynomial system solver:\n");
    printf("    %s \"polynomials\" field_size\n", prog_name);
    printf("    %s -s \"polynomials\" field_size\n", prog_name);
    printf("    %s --solve-rational-only \"polynomials\" 0\n", prog_name);
    printf("    %s -v 2 -s \"polynomials\" field_size\n", prog_name);
    printf("    -> Writes all solutions to solution_YYYYMMDD_HHMMSS.dr\n");
    printf("    -> `-s` / `--solve` is optional here; `--solve-rational-only` keeps only exact rational solutions\n");
    printf("    -> `-v 2` matches the old debug / verbose solver output\n");

    printf("  Complexity analysis:\n");
    printf("    %s --comp \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    %s -c    \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    %s --comp -f input.dr\n", prog_name);
    printf("    -> Prints complexity info; saves to comp_YYYYMMDD_HHMMSS.dr by default\n");
    printf("    Add --omega <value> (or -w <value>) to set omega (default: %.4g)\n",
           DIXON_OMEGA);
    printf("    Add --time to print per-step timing; use -v 2 for the old debug-level diagnostics\n");

    printf("  Dixon with ideal reduction:\n");
    printf("    %s --ideal \"ideal_generators\" \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    %s --ideal -f input.dr\n", prog_name);
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
    printf("    %s -r -n 4 --density 0.5 \"[d]*3\" field_size\n", prog_name);
    printf("    %s -r -s    \"[d1,...,dn]\" field_size\n", prog_name);
    printf("    %s -r --comp  \"[d]*n\"        field_size\n", prog_name);
    printf("    -> Add -n <num_vars> to set the total variable count (must satisfy num_vars >= #equations-1)\n");
    printf("    -> Add --density <ratio> with 0 <= ratio <= 1 to choose the fraction of all monomials used (default: 1)\n");
    printf("    -> Mixed degree specs such as \"[2]*5+[3]*6\" are supported\n");

    printf("  File input:\n");
    printf("    %s input_file\n", prog_name);
    printf("    %s -f input_file\n", prog_name);
    printf("    %s -s input_file\n", prog_name);
    printf("    %s -s -f input_file -o output.dr\n", prog_name);
    printf("    -> Without flags, auto-detects solver mode when line 1 starts with a digit; otherwise uses elimination mode\n");
    printf("    -> If elimination vars count equals equation count n, auto-adjusts to eliminate the first n-1 variables\n");
    printf("    -> Input file may be given directly or with `-f`; output file must be given with `-o` when overriding the default\n");
    printf("    -> Default file-mode outputs are input_file_solution.dr or input_file_comp.dr\n");

    printf("  Verbosity:\n");
    printf("    %s -v 0 <args>\n", prog_name);
    printf("    %s -v 1 <args>\n", prog_name);
    printf("    %s -v 2 <args>\n", prog_name);
    printf("    -> `-v 0` matches `--silent` and prints nothing\n");
    printf("    -> `-v 1` is the default output level\n");
    printf("    -> `-v 2` matches the old `--debug` output and also enables per-step timing\n");

    printf("  Diagnostics:\n");
    printf("    %s --time <args>\n", prog_name);
    printf("    %s -v 2 <args>\n", prog_name);
    printf("    -> --time prints per-step timing; interpolation steps also show CPU/Wall/Threads\n");
    printf("    -> `--silent`, `--debug`, `--solve-verbose` and `--solve` remain accepted for compatibility\n");

    printf("  Method selection:\n");
    printf("    %s --method <num> <args>\n", prog_name);
    printf("    %s --step1 <num> --step4 <num> <args>\n", prog_name);
    printf("    -> Available methods: 0.Recursive; 1.Kronecker+HNF; 2.Interpolation; 3.Sparse interpolation\n");
    printf("    -> --method sets both step 1 and step 4 for backward compatibility\n");
    printf("  Resultant construction:\n");
    printf("    %s --resultant dixon|macaulay|subres <args>\n", prog_name);
    printf("    %s --macaulay <args>\n", prog_name);
    printf("    %s --subres <args>\n", prog_name);
    printf("    -> --macaulay is shorthand for --resultant macaulay\n");
    printf("    -> --subres is for exactly 2 polynomials and 1 elimination variable\n");

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
    printf("    -> If line 1 lists n vars for n equations, compatibility mode uses the first n-1 variables\n");

    printf("EXAMPLES:\n");
    printf("  %s \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("  %s \"x^2+y^2+z^2-1, x^2+y^2-2*z^2, x+y+z\" \"x,y\" 0\n", prog_name);
    printf("  %s -s \"x^2+y^2+z^2-6, x+y+z-4, x*y*z-x-1\" 257\n", prog_name);
    printf("  %s --comp \"x^2+y^2+1, x*y+z, x+y+z^2\" \"x,y\" 257\n", prog_name);
    printf("  %s --random \"[3,3,2]\" 257\n", prog_name);
    printf("  %s -r \"[3]*3\" 0\n", prog_name);
    printf("  %s -r -n 4 --density 0.5 \"[3]*3\" 257\n", prog_name);
    printf("  %s -r \"[2]*4+[3]*2\" 257\n", prog_name);
    printf("  %s -r -s \"[2]*3\" 257\n", prog_name);
    printf("  %s -r --comp --omega 2.81 \"[4]*4\" 257\n", prog_name);
    printf("  %s --ideal \"a2^3=2*a1+1, a3^3=a1*a2+3\" \"a1^2+a2^2+a3^2-10, a3^3-a1*a2-3\" \"a3\" 257\n", prog_name);
    printf("  %s --field-eqution \"x0*x2+x1, x0*x1*x2+x2+1, x1*x2+x0+1\" \"x0,x1\" 2\n", prog_name);
    printf("  %s -v 0 \"x+y^2+t, x*y+t*y+1\" \"y\" 2^8\n", prog_name);
    printf("  %s \"x^2 + t*y, x*y + t^2\" \"2^8: t^8 + t^4 + t^3 + t + 1\"\n", prog_name);
    printf("  (AES polynomial for GF(2^8), 't' is the field extension generator)\n");
    printf("  %s example.dr\n", prog_name);
    printf("  %s -v 2 -f in.dr -o out.dr\n", prog_name);
    printf("  %s example_solve.dr\n", prog_name);
}

/* =========================================================================
 * Utility helpers
 * ========================================================================= */

static char *dixon_arb_to_string(const arb_t value, slong digits)
{
    char *buffer = NULL;
#ifdef _WIN32
    char *raw = arb_get_str(value, digits, 0);
    size_t start = 0;
    size_t end;
    size_t len;

    if (!raw) return NULL;

    end = strlen(raw);
    if (end >= 2 && raw[0] == '[' && raw[end - 1] == ']') {
        start = 1;
        end--;
    }

    len = end - start;
    buffer = (char *) malloc(len + 1);
    if (!buffer) {
        flint_free(raw);
        return NULL;
    }

    memcpy(buffer, raw + start, len);
    buffer[len] = '\0';
    flint_free(raw);
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

static const char *resultant_method_heading(resultant_method_t method)
{
    switch (method) {
        case RESULTANT_METHOD_MACAULAY:
            return "Macaulay Resultant Computation";
        case RESULTANT_METHOD_SUBRES:
            return "Subresultant Resultant Computation";
        case RESULTANT_METHOD_DIXON:
        default:
            return "Dixon Resultant Computation";
    }
}

static const char *resultant_method_task_label(resultant_method_t method)
{
    switch (method) {
        case RESULTANT_METHOD_MACAULAY:
            return "Macaulay resultant";
        case RESULTANT_METHOD_SUBRES:
            return "Subresultant resultant";
        case RESULTANT_METHOD_DIXON:
        default:
            return "Dixon resultant";
    }
}

static int validate_subres_input(const char *polys_str,
                                 const char *vars_str,
                                 int silent_mode)
{
    slong num_polys = 0, num_vars = 0;
    char **poly_array = split_string(polys_str, &num_polys);
    char **vars_array = split_string(vars_str, &num_vars);
    int ok = 1;

    if (num_polys != 2) {
        if (!silent_mode) {
            fprintf(stderr,
                    "Error: --subres supports exactly 2 polynomials, but got %ld.\n",
                    num_polys);
        }
        ok = 0;
    }

    if (ok && num_vars != 1) {
        if (!silent_mode) {
            fprintf(stderr,
                    "Error: --subres requires exactly 1 elimination variable, but got %ld.\n",
                    num_vars);
        }
        ok = 0;
    }

    free_split_strings(poly_array, num_polys);
    free_split_strings(vars_array, num_vars);
    return ok;
}

static char *compute_subres_resultant_str(const char *polys_str,
                                          const char *vars_str,
                                          const fq_nmod_ctx_t ctx)
{
    slong num_polys = 0, num_vars = 0;
    char **poly_array = split_string(polys_str, &num_polys);
    char **vars_array = split_string(vars_str, &num_vars);
    char *result = NULL;

    if (num_polys == 2 && num_vars == 1) {
        result = bivariate_resultant(poly_array[0], poly_array[1], vars_array[0], ctx);
    }

    free_split_strings(poly_array, num_polys);
    free_split_strings(vars_array, num_vars);
    return result;
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

static char *choose_output_filename(const char *requested_output,
                                    const char *input_filename,
                                    const char *timestamp_prefix,
                                    const char *tag)
{
    if (requested_output) {
        return strdup(requested_output);
    }
    if (input_filename && tag) {
        return generate_tagged_filename(input_filename, tag);
    }
    if (timestamp_prefix) {
        return generate_timestamped_filename(timestamp_prefix);
    }
    return NULL;
}

static int parse_verbose_level(const char *value, int *verbose_level)
{
    char *endptr = NULL;
    long parsed = strtol(value, &endptr, 10);

    if (!value || !verbose_level) return 0;
    if (!endptr || *endptr != '\0' || parsed < 0 || parsed > 2) {
        return 0;
    }

    *verbose_level = (int) parsed;
    return 1;
}

static int parse_positive_slong_option(const char *value, slong *out)
{
    char *endptr = NULL;
    long parsed = strtol(value, &endptr, 10);

    if (!value || !out) return 0;
    if (!endptr || *endptr != '\0' || parsed <= 0) return 0;

    *out = (slong) parsed;
    return 1;
}

static int parse_density_option(const char *value, double *density_out)
{
    char *endptr = NULL;
    double density = strtod(value, &endptr);

    if (!value || !density_out) return 0;
    if (!endptr || *endptr != '\0' || density < 0.0 || density > 1.0) return 0;

    *density_out = density;
    return 1;
}

static int append_text(char **buffer, size_t *capacity, size_t *length,
                       const char *text)
{
    size_t text_len;
    char *grown;

    if (!buffer || !capacity || !length || !text) return 0;

    text_len = strlen(text);
    if (*capacity == 0) {
        *capacity = text_len + 32;
        *buffer = (char *) malloc(*capacity);
        if (!*buffer) return 0;
        (*buffer)[0] = '\0';
        *length = 0;
    }

    if (*length + text_len + 1 > *capacity) {
        while (*length + text_len + 1 > *capacity) {
            *capacity *= 2;
        }
        grown = (char *) realloc(*buffer, *capacity);
        if (!grown) return 0;
        *buffer = grown;
    }

    memcpy(*buffer + *length, text, text_len + 1);
    *length += text_len;
    return 1;
}

static char *join_x_var_range(slong start, slong end_exclusive)
{
    char *result = NULL;
    size_t capacity = 0;
    size_t length = 0;

    if (end_exclusive < start) end_exclusive = start;

    for (slong i = start; i < end_exclusive; i++) {
        char piece[64];

        if (i > start && !append_text(&result, &capacity, &length, ",")) {
            free(result);
            return NULL;
        }

        snprintf(piece, sizeof(piece), "x%ld", i);
        if (!append_text(&result, &capacity, &length, piece)) {
            free(result);
            return NULL;
        }
    }

    if (!result) {
        result = strdup("");
    }

    return result;
}

static int build_random_system_strings(slong nvars,
                                       slong num_elim_vars,
                                       char **elim_vars_out,
                                       char **all_vars_out,
                                       char **remaining_vars_out)
{
    char *all_vars = NULL;
    char *elim_vars = NULL;
    char *remaining_vars = NULL;

    if (!elim_vars_out || !all_vars_out || !remaining_vars_out) return 0;
    if (nvars <= 0 || num_elim_vars < 0 || num_elim_vars > nvars) return 0;

    all_vars = join_x_var_range(0, nvars);
    elim_vars = join_x_var_range(0, num_elim_vars);
    remaining_vars = join_x_var_range(num_elim_vars, nvars);

    if (!all_vars || !elim_vars || !remaining_vars) {
        free(all_vars);
        free(elim_vars);
        free(remaining_vars);
        return 0;
    }

    *elim_vars_out = elim_vars;
    *all_vars_out = all_vars;
    *remaining_vars_out = remaining_vars;
    return 1;
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
 * Random polynomial generation helpers
 * ========================================================================= */

/*
 * Parse degree list in any of:
 *   "3,3,2"              plain comma-separated
 *   "[3,3,2]"            with brackets
 *   "[3]*5"              repeat shorthand
 *   "[3,2]*4"            repeat a sub-list
 *   "[2]*5+[3]*6"        sum of repeated blocks
 *
 * Returns malloc'd array; caller frees.  *count_out = number of entries.
 */
static long *parse_degree_list(const char *deg_str, slong *count_out)
{
    slong cap = 32;
    slong count = 0;
    long *degs = (long *) malloc((size_t) cap * sizeof(long));
    const char *p = deg_str;

    if (count_out) *count_out = 0;
    if (!deg_str || !degs) return NULL;

#define PUSH(v) do {     if (count >= cap) {         long *grown;         cap *= 2;         grown = (long *) realloc(degs, (size_t) cap * sizeof(long));         if (!grown) { free(degs); return NULL; }         degs = grown;     }     degs[count++] = (v); } while (0)

    while (*p) {
        while (*p && (isspace((unsigned char) *p) || *p == ',' || *p == '+')) p++;
        if (!*p) break;

        if (*p == '[') {
            const char *bracket_end = strchr(p + 1, ']');
            slong base_cap = 8;
            slong base_count = 0;
            long *base = NULL;

            if (!bracket_end) {
                fprintf(stderr, "Warning: malformed degree block '%s'\n", p);
                break;
            }

            base = (long *) malloc((size_t) base_cap * sizeof(long));
            if (!base) {
                free(degs);
                return NULL;
            }

            {
                const char *q = p + 1;
                while (q < bracket_end) {
                    char *endptr = NULL;
                    long d;

                    while (q < bracket_end &&
                           (isspace((unsigned char) *q) || *q == ',')) q++;
                    if (q >= bracket_end) break;

                    d = strtol(q, &endptr, 10);
                    if (endptr == q) {
                        q++;
                        continue;
                    }

                    if (d > 0) {
                        if (base_count >= base_cap) {
                            long *grown;
                            base_cap *= 2;
                            grown = (long *) realloc(base, (size_t) base_cap * sizeof(long));
                            if (!grown) {
                                free(base);
                                free(degs);
                                return NULL;
                            }
                            base = grown;
                        }
                        base[base_count++] = d;
                    } else {
                        fprintf(stderr, "Warning: degree %ld ignored (must be > 0)\n", d);
                    }
                    q = endptr;
                }
            }

            if (base_count > 0) {
                long repeat = 1;
                const char *after = bracket_end + 1;

                while (isspace((unsigned char) *after)) after++;
                if (*after == '*') {
                    char *endptr = NULL;
                    repeat = strtol(after + 1, &endptr, 10);
                    if (repeat <= 0) repeat = 1;
                    p = (endptr && endptr != after + 1) ? endptr : after + 1;
                } else {
                    p = after;
                }

                for (long r = 0; r < repeat; r++) {
                    for (slong j = 0; j < base_count; j++) {
                        PUSH(base[j]);
                    }
                }
            } else {
                p = bracket_end + 1;
            }

            free(base);
            continue;
        }

        {
            char *endptr = NULL;
            long d = strtol(p, &endptr, 10);

            if (endptr == p) {
                p++;
                continue;
            }

            if (d > 0) {
                PUSH(d);
            } else {
                fprintf(stderr, "Warning: degree %ld ignored (must be > 0)\n", d);
            }
            p = endptr;
        }
    }
#undef PUSH

    if (count == 0) {
        free(degs);
        return NULL;
    }

    if (count_out) *count_out = count;
    return degs;
}


/*
 * Generate a random polynomial system and return it as strings.
 *
 * All three output strings are malloc'd; caller must free them.
 * Returns 1 on success, 0 on failure.
 */
static void free_enumerated_monomials(monomial_t *monomials, slong count)
{
    if (!monomials) return;
    for (slong i = 0; i < count; i++) {
        free(monomials[i].exponents);
    }
    free(monomials);
}

static int append_signed_monomial_text(char **buffer,
                                       size_t *capacity,
                                       size_t *length,
                                       int coeff,
                                       const slong *exp,
                                       slong nvars,
                                       int *first_term)
{
    int abs_coeff;
    int has_var = 0;

    if (!buffer || !capacity || !length || !first_term) return 0;
    if (coeff == 0) return 0;

    abs_coeff = coeff < 0 ? -coeff : coeff;
    for (slong j = 0; j < nvars; j++) {
        if (exp[j] != 0) {
            has_var = 1;
            break;
        }
    }

    if (*first_term) {
        if (coeff < 0 && !append_text(buffer, capacity, length, "-")) return 0;
    } else {
        if (!append_text(buffer, capacity, length, coeff < 0 ? " - " : " + ")) return 0;
    }

    if (abs_coeff != 1 || !has_var) {
        char coeff_buf[32];
        snprintf(coeff_buf, sizeof(coeff_buf), "%d", abs_coeff);
        if (!append_text(buffer, capacity, length, coeff_buf)) return 0;
    }

    if (has_var) {
        int wrote_any_var = 0;
        for (slong j = 0; j < nvars; j++) {
            char piece[64];
            if (exp[j] == 0) continue;

            if ((abs_coeff != 1) || wrote_any_var) {
                if (!append_text(buffer, capacity, length, "*")) return 0;
            }

            snprintf(piece, sizeof(piece), "x%ld", j);
            if (!append_text(buffer, capacity, length, piece)) return 0;

            if (exp[j] != 1) {
                snprintf(piece, sizeof(piece), "^%ld", exp[j]);
                if (!append_text(buffer, capacity, length, piece)) return 0;
            }
            wrote_any_var = 1;
        }
    }

    *first_term = 0;
    return 1;
}

static char *build_random_system_spec(const char *deg_spec,
                                      slong nvars,
                                      double density_ratio)
{
    size_t spec_len;
    char *spec;

    if (!deg_spec) return strdup("random system");

    spec_len = strlen(deg_spec) + 128;
    spec = (char *) malloc(spec_len);
    if (!spec) return NULL;

    snprintf(spec, spec_len,
             "random system spec: degrees=%s, variables=%ld, density=%.6g",
             deg_spec, nvars, density_ratio);
    return spec;
}

static int generate_random_poly_strings(
        const long *degrees, slong npolys,
        slong nvars,
        slong num_elim_vars,
        double density_ratio,
        const fq_nmod_ctx_t ctx,
        int silent_mode,
        char **polys_str_out,
        char **elim_vars_str_out,
        char **all_vars_str_out)
{
    slong npars = 0;
    char *all_vars = NULL;
    char *elim_vars = NULL;
    char *remaining_vars = NULL;
    slong *slong_deg = NULL;
    fq_mvpoly_t *polys = NULL;
    char *gen_name = NULL;
    char **poly_strs = NULL;
    char *polys_str = NULL;
    size_t total_len = 4;
    flint_rand_t rstate;
    int orig_stdout = -1;
    int devnull = -1;
    int rand_initialized = 0;

    if (!degrees || npolys <= 0 || nvars <= 0 || num_elim_vars < 0 || num_elim_vars > nvars)
        return 0;

    if (!build_random_system_strings(nvars, num_elim_vars,
                                     &elim_vars, &all_vars, &remaining_vars)) {
        return 0;
    }

    slong_deg = (slong *) malloc((size_t) npolys * sizeof(slong));
    if (!slong_deg) goto fail;
    for (slong i = 0; i < npolys; i++) slong_deg[i] = (slong) degrees[i];

    flint_rand_init(rstate);
    rand_initialized = 1;
    flint_rand_set_seed(rstate, (ulong) time(NULL),
                        (ulong) time(NULL) ^ (ulong) clock());

    if (!silent_mode) printf("Generating random polynomial system...\n");

    fflush(stdout);
    orig_stdout = dup(STDOUT_FILENO);
    devnull = open(DIXON_NULL_DEVICE, O_WRONLY);
    if (devnull != -1) {
        dup2(devnull, STDOUT_FILENO);
        close(devnull);
    }

    generate_polynomial_system(&polys, nvars, npolys, npars,
                               slong_deg, density_ratio, ctx, rstate);

    fflush(stdout);
    if (orig_stdout != -1) {
        dup2(orig_stdout, STDOUT_FILENO);
        close(orig_stdout);
        orig_stdout = -1;
    }

    gen_name = get_generator_name(ctx);
    poly_strs = (char **) malloc((size_t) npolys * sizeof(char *));
    if (!poly_strs) goto fail;

    for (slong i = 0; i < npolys; i++) {
        char *s = fq_mvpoly_to_string(&polys[i], NULL, gen_name);
        poly_strs[i] = (s && strlen(s) > 0) ? s : (free(s), strdup("0"));
        if (!poly_strs[i]) goto fail;
        total_len += strlen(poly_strs[i]) + 3;
    }

    polys_str = (char *) malloc(total_len);
    if (!polys_str) goto fail;
    polys_str[0] = '\0';
    for (slong i = 0; i < npolys; i++) {
        if (i > 0) strcat(polys_str, ", ");
        strcat(polys_str, poly_strs[i]);
    }

    if (!silent_mode) {
        printf("System: %ld equations, %ld variables\n", npolys, nvars);
        printf("Degrees: [");
        for (slong i = 0; i < npolys; i++) {
            if (i > 0) printf(", ");
            printf("%ld", slong_deg[i]);
        }
        printf("]\n");
        printf("Density: %.2f%% of all monomials up to each polynomial degree\n",
               density_ratio * 100.0);
        if (num_elim_vars > 0) {
            printf("Eliminate: %s\n", elim_vars);
        }
        if (strlen(remaining_vars) > 0) {
            printf("Remaining: %s\n", remaining_vars);
        }
    }

    for (slong i = 0; i < npolys; i++) free(poly_strs[i]);
    free(poly_strs);
    for (slong i = 0; i < npolys; i++) fq_mvpoly_clear(&polys[i]);
    free(polys);
    free(slong_deg);
    if (gen_name) free(gen_name);
    free(remaining_vars);
    if (rand_initialized) flint_rand_clear(rstate);

    *polys_str_out = polys_str;
    *elim_vars_str_out = elim_vars;
    *all_vars_str_out = all_vars;
    return 1;

fail:
    if (orig_stdout != -1) {
        fflush(stdout);
        dup2(orig_stdout, STDOUT_FILENO);
        close(orig_stdout);
    }
    if (poly_strs) {
        for (slong i = 0; i < npolys; i++) free(poly_strs[i]);
        free(poly_strs);
    }
    if (polys) {
        for (slong i = 0; i < npolys; i++) fq_mvpoly_clear(&polys[i]);
        free(polys);
    }
    free(polys_str);
    free(slong_deg);
    free(all_vars);
    free(elim_vars);
    free(remaining_vars);
    if (gen_name) free(gen_name);
    if (rand_initialized) flint_rand_clear(rstate);
    return 0;
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
        *field_str = lines[0];
        *polys_str = join_polynomial_lines(lines, 1, line_count);
        if (!*polys_str) {
            free_input_lines(lines, line_count);
            return 0;
        }
        *solve_mode_out = 1;
        *vars_str = NULL;
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

    if (sols->num_candidate_solution_lines > 0) {
        fprintf(out_fp, "Candidate sets: %ld\n", sols->num_candidate_solution_lines);
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
        if (sols->num_candidate_solution_lines > 0) {
            fprintf(out_fp, "Candidate sets: %ld\n", sols->num_candidate_solution_lines);
        }
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
        if (sols->num_solution_sets == 0 && sols->num_candidate_solution_lines > 0) {
            fprintf(out_fp, "Candidate sets: %ld\n", sols->num_candidate_solution_lines);
        }
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

    fprintf(out_fp, "%s\n", resultant_method_heading(g_resultant_method));
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

static int auto_adjust_elimination_vars_for_msolve_compat(const char *polys_str,
                                                          const char *field_str,
                                                          const char *input_vars_str,
                                                          char **vars_str_out,
                                                          int silent_mode)
{
    slong num_polys = 0;
    char **poly_arr = NULL;
    char **input_vars = NULL;
    char *vars_str = NULL;
    int ok = 0;

    slong num_input_vars = 0;

    if (vars_str_out) *vars_str_out = NULL;

    if (!polys_str || !input_vars_str || !vars_str_out) return 0;

    poly_arr = split_string(polys_str, &num_polys);
    if (!poly_arr || num_polys <= 0) goto cleanup;

    input_vars = split_string(input_vars_str, &num_input_vars);
    if (!input_vars) goto cleanup;

    (void) field_str;

    if (num_input_vars == num_polys && num_polys > 1) {
        size_t total_len = 1;
        for (slong i = 0; i < num_polys - 1; i++)
            total_len += strlen(input_vars[i]) + 2;

        vars_str = (char *) malloc(total_len);
        if (!vars_str) goto cleanup;

        vars_str[0] = '\0';
        for (slong i = 0; i < num_polys - 1; i++) {
            if (i > 0) strcat(vars_str, ",");
            strcat(vars_str, input_vars[i]);
        }

    }

    *vars_str_out = vars_str;
    vars_str = NULL;
    ok = 1;

cleanup:
    if (poly_arr) free_split_strings(poly_arr, num_polys);
    if (input_vars) free_split_strings(input_vars, num_input_vars);
    free(vars_str);
    return ok;
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
        slong nvars,
        slong num_elim_vars,
        double density_ratio,
        int silent_mode,
        char **polys_str_out,
        char **elim_vars_str_out,
        char **all_vars_str_out)
{
    const int coeff_choices[] = { -2, -1, 1, 2 };
    char *all_vars = NULL;
    char *elim_vars = NULL;
    char *remaining_vars = NULL;
    char *polys_str = NULL;
    size_t polys_cap = 0;
    size_t polys_len = 0;

    if (!degrees || npolys <= 0 || nvars <= 0 || num_elim_vars < 0 || num_elim_vars > nvars)
        return 0;

    if (!build_random_system_strings(nvars, num_elim_vars,
                                     &elim_vars, &all_vars, &remaining_vars)) {
        return 0;
    }

    srand((unsigned) (time(NULL) ^ clock()));

    if (!silent_mode) printf("Generating random polynomial system over Q...\n");

    for (slong i = 0; i < npolys; i++) {
        monomial_t *monomials = NULL;
        slong monomial_count = 0;
        slong degree = degrees[i] > 0 ? degrees[i] : 1;
        slong target_terms;
        slong *indices = NULL;
        slong leading_idx = -1;
        char *poly_buf = NULL;
        size_t poly_cap = 0;
        size_t poly_len = 0;
        int first_term = 1;
        slong selected = 0;

        enumerate_all_monomials(&monomials, &monomial_count, nvars, degree);
        if (!monomials || monomial_count <= 0) {
            free(indices);
            free(poly_buf);
            free_enumerated_monomials(monomials, monomial_count);
            free(polys_str);
            free(all_vars);
            free(elim_vars);
            free(remaining_vars);
            return 0;
        }

        target_terms = (slong) (density_ratio * monomial_count);
        if (target_terms < 1) target_terms = 1;
        if (target_terms > monomial_count) target_terms = monomial_count;

        indices = (slong *) malloc((size_t) monomial_count * sizeof(slong));
        if (!indices) {
            free_enumerated_monomials(monomials, monomial_count);
            free(polys_str);
            free(all_vars);
            free(elim_vars);
            free(remaining_vars);
            return 0;
        }

        for (slong j = 0; j < monomial_count; j++) indices[j] = j;
        for (slong j = monomial_count - 1; j > 0; j--) {
            slong k = (slong) (rand() % (j + 1));
            slong tmp = indices[j];
            indices[j] = indices[k];
            indices[k] = tmp;
        }

        for (slong j = 0; j < monomial_count; j++) {
            if (monomials[indices[j]].total_degree == degree) {
                leading_idx = indices[j];
                break;
            }
        }
        if (leading_idx < 0) leading_idx = indices[0];

        if (!append_signed_monomial_text(&poly_buf, &poly_cap, &poly_len,
                                         coeff_choices[rand() % 4],
                                         monomials[leading_idx].exponents,
                                         nvars, &first_term)) {
            free(indices);
            free(poly_buf);
            free_enumerated_monomials(monomials, monomial_count);
            free(polys_str);
            free(all_vars);
            free(elim_vars);
            free(remaining_vars);
            return 0;
        }
        selected++;

        for (slong j = 0; j < monomial_count && selected < target_terms; j++) {
            slong idx = indices[j];
            if (idx == leading_idx) continue;
            if (!append_signed_monomial_text(&poly_buf, &poly_cap, &poly_len,
                                             coeff_choices[rand() % 4],
                                             monomials[idx].exponents,
                                             nvars, &first_term)) {
                free(indices);
                free(poly_buf);
                free_enumerated_monomials(monomials, monomial_count);
                free(polys_str);
                free(all_vars);
                free(elim_vars);
                free(remaining_vars);
                return 0;
            }
            selected++;
        }

        if (i > 0 && !append_text(&polys_str, &polys_cap, &polys_len, ", ")) {
            free(indices);
            free(poly_buf);
            free_enumerated_monomials(monomials, monomial_count);
            free(polys_str);
            free(all_vars);
            free(elim_vars);
            free(remaining_vars);
            return 0;
        }
        if (!append_text(&polys_str, &polys_cap, &polys_len, poly_buf ? poly_buf : "0")) {
            free(indices);
            free(poly_buf);
            free_enumerated_monomials(monomials, monomial_count);
            free(polys_str);
            free(all_vars);
            free(elim_vars);
            free(remaining_vars);
            return 0;
        }

        free(indices);
        free(poly_buf);
        free_enumerated_monomials(monomials, monomial_count);
    }

    if (!silent_mode) {
        printf("System: %ld equations, %ld variables\n", npolys, nvars);
        printf("Degrees: [");
        for (slong i = 0; i < npolys; i++) {
            if (i > 0) printf(", ");
            printf("%ld", degrees[i]);
        }
        printf("]\n");
        printf("Density: %.2f%% of all monomials up to each polynomial degree\n",
               density_ratio * 100.0);
        if (num_elim_vars > 0) {
            printf("Eliminate: %s\n", elim_vars);
        }
        if (strlen(remaining_vars) > 0) {
            printf("Remaining: %s\n", remaining_vars);
        }
    }

    free(remaining_vars);
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

    /* ---- parse flags ---- */
    int    verbose_level = 1;
    int    solve_mode  = 0;
    int    solve_rational_only_mode = 0;
    int    comp_mode   = 0;
    int    rand_mode   = 0;   /* --random / -r */
    int    ideal_mode  = 0;   /*  --ideal flag */
    int    field_eq_mode = 0; /* --field-equation */
    int    time_mode   = 0;   /* --time */
    double omega       = DIXON_OMEGA;   /* default, overridden by --omega */
    int    det_method_step1 = -1;  /* determinant method override for step 1 */
    int    det_method_step4 = -1;  /* determinant method override for step 4 */
    int    num_threads = -1;  /* number of threads, -1 means use default */
    resultant_method_t resultant_method = RESULTANT_METHOD_DIXON;
    const char *cli_input_filename = NULL;
    const char *cli_output_filename = NULL;
    char *positional_args[argc];
    int positional_count = 0;
    slong random_nvars = 0;
    int random_nvars_given = 0;
    double random_density = 1.0;
    int random_density_given = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--silent") == 0) {
            verbose_level = 0;
        } else if (strcmp(argv[i], "--solve-verbose") == 0) {
            solve_mode = 1;
            verbose_level = 2;
        } else if (strcmp(argv[i], "--solve-rational-only") == 0) {
            solve_mode = 1;
            solve_rational_only_mode = 1;
        } else if (strcmp(argv[i], "--solve") == 0 ||
                   strcmp(argv[i], "-s")      == 0) {
            solve_mode = 1;
        } else if (strcmp(argv[i], "--comp") == 0 ||
                   strcmp(argv[i], "-c")     == 0) {
            comp_mode = 1;
        } else if (strcmp(argv[i], "--random") == 0 ||
                   strcmp(argv[i], "-r")       == 0) {
            rand_mode = 1;
        } else if (strcmp(argv[i], "--ideal") == 0) {  /* <<< NEW >>> */
            ideal_mode = 1;
        } else if (strcmp(argv[i], "--field-equation") == 0 ||
                   strcmp(argv[i], "--field-eqution")  == 0) {
            field_eq_mode = 1;
        } else if (strcmp(argv[i], "--time") == 0) {
            time_mode = 1;
        } else if (strcmp(argv[i], "--debug") == 0) {
            verbose_level = 2;
        } else if ((strcmp(argv[i], "--verbose") == 0 ||
                    strcmp(argv[i], "-v")        == 0) && i + 1 < argc) {
            if (!parse_verbose_level(argv[i + 1], &verbose_level)) {
                fprintf(stderr, "Error: invalid verbose level '%s'; expected 0, 1, or 2.\n",
                        argv[i + 1]);
                return 1;
            }
            i++;
        } else if (strcmp(argv[i], "--verbose") == 0 ||
                   strcmp(argv[i], "-v")        == 0) {
            fprintf(stderr, "Error: %s requires an integer argument 0, 1, or 2.\n", argv[i]);
            return 1;
        } else if (strcmp(argv[i], "--macaulay") == 0) {
            resultant_method = RESULTANT_METHOD_MACAULAY;
        } else if (strcmp(argv[i], "--subres") == 0) {
            resultant_method = RESULTANT_METHOD_SUBRES;
        } else if ((strcmp(argv[i], "--resultant") == 0 ||
                    strcmp(argv[i], "--resultant-method") == 0) && i + 1 < argc) {
            if (strcmp(argv[i + 1], "macaulay") == 0) {
                resultant_method = RESULTANT_METHOD_MACAULAY;
            } else if (strcmp(argv[i + 1], "subres") == 0) {
                resultant_method = RESULTANT_METHOD_SUBRES;
            } else if (strcmp(argv[i + 1], "dixon") == 0) {
                resultant_method = RESULTANT_METHOD_DIXON;
            } else {
                fprintf(stderr, "Warning: invalid resultant method '%s', using dixon.\n",
                                argv[i + 1]);
            }
            i++;
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
            i++;          /* skip the value token */
        } else if ((strcmp(argv[i], "-n") == 0 ||
                    strcmp(argv[i], "--nvars") == 0 ||
                    strcmp(argv[i], "--num-vars") == 0) && i + 1 < argc) {
            if (!parse_positive_slong_option(argv[i + 1], &random_nvars)) {
                fprintf(stderr, "Error: invalid variable count '%s'; expected a positive integer.\n",
                        argv[i + 1]);
                return 1;
            }
            random_nvars_given = 1;
            i++;
        } else if (strcmp(argv[i], "-n") == 0 ||
                   strcmp(argv[i], "--nvars") == 0 ||
                   strcmp(argv[i], "--num-vars") == 0) {
            fprintf(stderr, "Error: %s requires a positive integer argument.\n", argv[i]);
            return 1;
        } else if (strcmp(argv[i], "--density") == 0 && i + 1 < argc) {
            if (!parse_density_option(argv[i + 1], &random_density)) {
                fprintf(stderr, "Error: invalid density '%s'; expected 0 <= density <= 1.\n",
                        argv[i + 1]);
                return 1;
            }
            random_density_given = 1;
            i++;
        } else if (strcmp(argv[i], "--density") == 0) {
            fprintf(stderr, "Error: --density requires a numeric argument between 0 and 1.\n");
            return 1;
        } else if (strcmp(argv[i], "-f") == 0 && i + 1 < argc) {
            cli_input_filename = argv[i + 1];
            i++;
        } else if (strcmp(argv[i], "-f") == 0) {
            fprintf(stderr, "Error: -f requires an input filename.\n");
            return 1;
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            cli_output_filename = argv[i + 1];
            i++;
        } else if (strcmp(argv[i], "-o") == 0) {
            fprintf(stderr, "Error: -o requires an output filename.\n");
            return 1;
        } else {
            positional_args[positional_count++] = argv[i];
        }
    }

    int silent_mode = (verbose_level == 0);
    int solve_verbose_mode = (verbose_level >= 2);
    int debug_mode = (verbose_level >= 2);

    if (!rand_mode && (random_nvars_given || random_density_given)) {
        fprintf(stderr, "Error: -n/--nvars and --density may only be used together with --random.\n");
        return 1;
    }


    
    /* ---- help ---- */
    if (positional_count >= 1 &&
        (strcmp(positional_args[0], "--help") == 0 ||
         strcmp(positional_args[0], "-h")     == 0)) {
        if (!silent_mode)  { print_version(); print_usage(prog_name); }
        return 0;
    }

    /* ---- version banner ---- */
    if (!silent_mode) {
        print_version();
    }

    /* ---- test modes ---- */
    if (positional_count >= 1 && strcmp(positional_args[0], "--test") == 0) {
        if (!silent_mode) {
            if (positional_count >= 2) test_dixon(atoi(positional_args[1]));
            else           test_dixon(0);
        }
        return 0;
    }
    if (positional_count >= 1 &&
        strcmp(positional_args[0], "--test-solver") == 0) {
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
    int   vars_str_generated = 0;
    int   need_free     = 0;
    int   rand_generated = 0;  /* 1 = polys_str/vars_str were malloc'd by random gen */
    int   rand_comp_direct = 0;
    char *deg_str       = NULL; /* degree list string when rand_mode */
    char *rand_comp_spec = NULL;
    char *input_filename  = NULL;
    char *output_filename = NULL;
    mp_limb_t prime = 0;
    ulong     power = 0;
    const char *file_input_arg = cli_input_filename;
    long *rand_comp_degrees = NULL;
    slong rand_comp_npolys = 0;
    slong rand_comp_nvars = 0;

    /* ---- determine input mode ---- */
    if (rand_mode) {
        /* --random / -r: positional args are "d1,d2,...,dn" field_size */
        if (file_input_arg || positional_count != 2) {
            if (!silent_mode) {
                fprintf(stderr, "Error: Random mode requires exactly:\n");
                fprintf(stderr, "  %s [-s|--comp] --random \"d1,d2,...,dn\" field_size\n",
                        prog_name);
            }
            return 1;
        }
        deg_str   = positional_args[0];
        field_str = positional_args[1];

        output_filename = choose_output_filename(cli_output_filename, NULL,
                                                 comp_mode ? "comp" : "solution", NULL);
    /*  --ideal mode */
    } else if (ideal_mode) {
        if (file_input_arg && positional_count > 0) {
            if (!silent_mode)
                fprintf(stderr, "Error: use either -f input_file or a bare input_file, not both.\n");
            return 1;
        }

        if (file_input_arg || positional_count == 1) {
            /* File input: auto-detect ideal generators by '=' */
            const char *input_arg = file_input_arg ? file_input_arg : positional_args[0];
            FILE *fp = fopen(input_arg, "r");
            if (!fp) {
                if (!silent_mode)
                    fprintf(stderr, "Error: Cannot open file '%s'\n", input_arg);
                return 1;
            }
            if (!silent_mode)
                printf("Reading from file (--ideal mode): %s\n", input_arg);

            input_filename  = strdup(input_arg);
            output_filename = choose_output_filename(cli_output_filename, input_filename,
                                                     NULL, "_solution");

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

        } else if (positional_count == 4) {
            /* CLI: --ideal "generators" "polys" "elim_vars" field_size */
            ideal_str = positional_args[0];
            polys_str = positional_args[1];
            vars_str  = positional_args[2];
            field_str = positional_args[3];

            output_filename = choose_output_filename(cli_output_filename, NULL,
                                                     "solution", NULL);

        } else {
            if (!silent_mode) {
                fprintf(stderr, "Error: --ideal mode requires either:\n");
                fprintf(stderr, "  %s --ideal \"ideal_generators\" \"polynomials\" \"eliminate_vars\" field_size\n",
                        prog_name);
                fprintf(stderr, "  %s --ideal [-f] input_file\n", prog_name);
            }
            return 1;
        }
    } else if (comp_mode) {
        /* Complexity analysis mode: same argument format as basic Dixon */
        if (file_input_arg && positional_count > 0) {
            if (!silent_mode)
                fprintf(stderr, "Error: use either -f input_file or a bare input_file, not both.\n");
            return 1;
        }

        if (file_input_arg || positional_count == 1) {
            /* file input */
            const char *input_arg = file_input_arg ? file_input_arg : positional_args[0];
            FILE *fp = fopen(input_arg, "r");
            if (!fp) {
                if (!silent_mode)
                    fprintf(stderr, "Error: Cannot open file '%s'\n",
                            input_arg);
                return 1;
            }
            if (!silent_mode)
                printf("Reading from file: %s\n", input_arg);

            input_filename  = strdup(input_arg);
            output_filename = choose_output_filename(cli_output_filename, input_filename,
                                                     NULL, "_comp");

            if (!read_multiline_file(fp, &field_str, &polys_str,
                                     &vars_str, &ideal_str, &allvars_str)) {
                fclose(fp); return 1;
            }
            fclose(fp);
            need_free = 1;

        } else if (positional_count == 3) {
            /* command line: polys vars field_size */
            polys_str = positional_args[0];
            vars_str  = positional_args[1];
            field_str = positional_args[2];

            output_filename = choose_output_filename(cli_output_filename, NULL,
                                                     "comp", NULL);

        } else {
            if (!silent_mode) {
                fprintf(stderr, "Error: Complexity mode requires:\n");
                fprintf(stderr, "  %s --comp \"polynomials\" \"eliminate_vars\" field_size\n",
                        prog_name);
                fprintf(stderr, "  %s --comp [-f] input_file\n", prog_name);
            }
            return 1;
        }

    } else if (solve_mode) {
        if (file_input_arg && positional_count > 0) {
            if (!silent_mode)
                fprintf(stderr, "Error: use either -f input_file or a bare input_file, not both.\n");
            return 1;
        }

        if (file_input_arg || positional_count == 1) {
            const char *input_arg = file_input_arg ? file_input_arg : positional_args[0];
            FILE *fp = fopen(input_arg, "r");
            int detected_solve_mode = 0;
            if (!fp) {
                if (!silent_mode)
                    fprintf(stderr, "Error: Cannot open file '%s'\n",
                            input_arg);
                return 1;
            }
            if (!silent_mode)
                printf("Reading polynomial system from file: %s\n",
                       input_arg);

            input_filename  = strdup(input_arg);
            output_filename = choose_output_filename(cli_output_filename, input_filename,
                                                     NULL, "_solution");

            if (!read_auto_input_file(fp, &detected_solve_mode, &field_str, &polys_str,
                                      &vars_str, &ideal_str, &allvars_str)) {
                fclose(fp); return 1;
            }
            fclose(fp);
            need_free = 1;

            if (vars_str) {
                free(vars_str);
                vars_str = NULL;
            }
            if (ideal_str) {
                free(ideal_str);
                ideal_str = NULL;
            }
            if (allvars_str) {
                free(allvars_str);
                allvars_str = NULL;
            }

        } else if (positional_count == 2) {
            polys_str = positional_args[0];
            field_str = positional_args[1];
            output_filename = choose_output_filename(cli_output_filename, NULL,
                                                     "solution", NULL);

        } else {
            if (!silent_mode) {
                fprintf(stderr, "Error: Solver mode requires either:\n");
                fprintf(stderr, "  %s -s \"polynomials\" field_size\n", prog_name);
                fprintf(stderr, "  %s -s [-f] input_file\n", prog_name);
            }
            return 1;
        }

    } else {
        /* Basic Dixon, or auto-detect solve/file mode when flags are omitted */
        if (file_input_arg && positional_count > 0) {
            if (!silent_mode)
                fprintf(stderr, "Error: use either -f input_file or a bare input_file, not both.\n");
            return 1;
        }

        if (positional_count == 2) {
            /* Auto-detect solve mode: only polynomials and field size */
            solve_mode = 1;
            polys_str = positional_args[0];
            field_str = positional_args[1];
            output_filename = choose_output_filename(cli_output_filename, NULL,
                                                     "solution", NULL);
        } else if (file_input_arg || positional_count == 1) {
            const char *input_arg = file_input_arg ? file_input_arg : positional_args[0];
            FILE *fp = fopen(input_arg, "r");
            int detected_solve_mode = 0;
            if (!fp) {
                if (!silent_mode)
                    fprintf(stderr, "Error: Cannot open file '%s'\n",
                            input_arg);
                return 1;
            }
            if (!silent_mode)
                printf("Reading from file: %s\n", input_arg);

            input_filename  = strdup(input_arg);
            output_filename = choose_output_filename(cli_output_filename, input_filename,
                                                     NULL, "_solution");

            if (!read_auto_input_file(fp, &detected_solve_mode, &field_str, &polys_str,
                                      &vars_str, &ideal_str, &allvars_str)) {
                fclose(fp); return 1;
            }
            fclose(fp);
            solve_mode = detected_solve_mode;
            need_free = 1;

            if (!silent_mode)
                printf("Detected file mode: %s\n", solve_mode ? "solver" : "elimination");

        } else if (positional_count == 3 || positional_count == 4) {
            polys_str = positional_args[0];
            vars_str  = positional_args[1];
            if (positional_count == 4) {
                ideal_str = positional_args[2];
                field_str = positional_args[3];
            } else {
                field_str = positional_args[2];
            }
            output_filename = choose_output_filename(cli_output_filename, NULL,
                                                     "solution", NULL);

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
            printf("=== %s ===\n", resultant_method_heading(resultant_method));
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
    } else if (solve_rational_only_mode) {
        fprintf(stderr, "Error: --solve-rational-only is only supported for field_size=0.\n");
        goto cleanup_fail;
    }
    if (large_prime_mode) {
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
        slong npolys_rand = 0;
        long *degrees_rand = parse_degree_list(deg_str, &npolys_rand);
        slong nvars_rand = random_nvars_given ? random_nvars : npolys_rand;
        slong num_elim_rand = comp_mode ? (npolys_rand - 1)
                                        : (solve_mode ? nvars_rand : (npolys_rand - 1));

        if (!degrees_rand || npolys_rand < 2) {
            if (!silent_mode)
                fprintf(stderr, "Error: --random requires at least 2 degrees "
                                "(e.g. \"3,3,2\")\n");
            free(degrees_rand);
            goto cleanup_fail;
        }

        if (nvars_rand < npolys_rand - 1) {
            if (!silent_mode)
                fprintf(stderr, "Error: random mode requires -n/--nvars >= #equations-1 = %ld (got %ld).\n",
                        npolys_rand - 1, nvars_rand);
            free(degrees_rand);
            goto cleanup_fail;
        }

        if (num_elim_rand < 0 || num_elim_rand > nvars_rand) {
            if (!silent_mode)
                fprintf(stderr, "Error: invalid random elimination-variable count %ld for %ld variables.\n",
                        num_elim_rand, nvars_rand);
            free(degrees_rand);
            goto cleanup_fail;
        }

        if (comp_mode) {
            char *gen_elim = NULL, *gen_allvars = NULL, *gen_remaining = NULL;

            rand_comp_spec = build_random_system_spec(deg_str, nvars_rand, random_density);
            if (!rand_comp_spec ||
                !build_random_system_strings(nvars_rand, num_elim_rand,
                                             &gen_elim, &gen_allvars, &gen_remaining)) {
                free(rand_comp_spec);
                rand_comp_spec = NULL;
                free(gen_elim);
                free(gen_allvars);
                free(gen_remaining);
                free(degrees_rand);
                if (!silent_mode)
                    fprintf(stderr, "Error: random complexity setup failed\n");
                goto cleanup_fail;
            }

            if (!silent_mode) {
                printf("Random system spec: %ld equations, %ld variables\n",
                       npolys_rand, nvars_rand);
                printf("Degrees: [");
                for (slong i = 0; i < npolys_rand; i++) {
                    if (i > 0) printf(", ");
                    printf("%ld", degrees_rand[i]);
                }
                printf("]\n");
                printf("Density: %.2f%% of all monomials up to each polynomial degree\n",
                       random_density * 100.0);
                if (strlen(gen_elim) > 0) {
                    printf("Eliminate: %s\n", gen_elim);
                }
                if (strlen(gen_remaining) > 0) {
                    printf("Remaining: %s\n", gen_remaining);
                }
            }

            free(gen_remaining);
            vars_str = gen_elim;
            allvars_str = gen_allvars;
            rand_generated = 1;
            rand_comp_direct = 1;
            rand_comp_degrees = degrees_rand;
            rand_comp_npolys = npolys_rand;
            rand_comp_nvars = nvars_rand;

        } else {
            char *gen_polys = NULL, *gen_elim = NULL, *gen_allvars = NULL;
            int generated_ok;

            if (rational_mode) {
                generated_ok = generate_random_poly_strings_rational(
                                   degrees_rand, npolys_rand,
                                   nvars_rand, num_elim_rand,
                                   random_density,
                                   silent_mode,
                                   &gen_polys, &gen_elim, &gen_allvars);
            } else {
                generated_ok = generate_random_poly_strings(
                                   degrees_rand, npolys_rand,
                                   nvars_rand, num_elim_rand,
                                   random_density,
                                   ctx,
                                   silent_mode,
                                   &gen_polys, &gen_elim, &gen_allvars);
            }
            if (!generated_ok) {
                if (!silent_mode)
                    fprintf(stderr, "Error: random polynomial generation failed\n");
                free(degrees_rand);
                goto cleanup_fail;
            }

            free(degrees_rand);
            polys_str = gen_polys;
            vars_str = gen_elim;
            allvars_str = gen_allvars;
            rand_generated = 1;
        }
    }

    if (!solve_mode && polys_str && vars_str) {
        char *compat_vars = NULL;
        if (!auto_adjust_elimination_vars_for_msolve_compat(polys_str, field_str, vars_str,
                                                            &compat_vars, silent_mode)) {
            if (!silent_mode)
                fprintf(stderr, "Error: Failed to analyze elimination variables for msolve compatibility\n");
            goto cleanup_fail;
        }
        if (compat_vars) {
            if (need_free || rand_generated || vars_str_generated) {
                free(vars_str);
            }
            vars_str_generated = 1;
            vars_str = compat_vars;
        }
    }

    if (polys_str && vars_str) {
        char *normalized_vars = normalize_elimination_vars(polys_str, vars_str, silent_mode);
        if (normalized_vars) {
            if (need_free || rand_generated || vars_str_generated) {
                free(vars_str);
                vars_str_generated = 1;
            } else {
                vars_str_override = normalized_vars;
            }
            vars_str = normalized_vars;
        }
    }

    if (resultant_method == RESULTANT_METHOD_SUBRES) {
        if (comp_mode) {
            fprintf(stderr, "Error: --subres does not support --comp.\n");
            goto cleanup_fail;
        }
        if (solve_mode) {
            fprintf(stderr, "Error: --subres does not support --solve.\n");
            goto cleanup_fail;
        }
        if (ideal_str) {
            fprintf(stderr, "Error: --subres does not support --ideal.\n");
            goto cleanup_fail;
        }
        if (rational_mode) {
            fprintf(stderr, "Error: --subres currently requires a finite field.\n");
            goto cleanup_fail;
        }
        if (large_prime_mode) {
            fprintf(stderr, "Error: --subres currently does not support large-prime fallback.\n");
            goto cleanup_fail;
        }
        if (!validate_subres_input(polys_str, vars_str, silent_mode)) {
            goto cleanup_fail;
        }
    }

    /* ======================================================
     * Set determinant method and thread count
     * ====================================================== */
    dixon_global_method_step1 = -1;
    dixon_global_method_step4 = -1;
    dixon_global_method = -1;
    g_resultant_method = resultant_method;
    g_dixon_verbose_level = verbose_level;
    g_dixon_show_step_timing = (!silent_mode) && (time_mode || debug_mode);
    g_dixon_debug_mode = debug_mode;
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
    large_prime_solutions_t *large_prime_solutions = NULL;

    dixon_clear_last_root_report();

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

        if (rand_comp_direct) {
            run_complexity_analysis_from_degrees(rand_comp_degrees,
                                                 rand_comp_npolys,
                                                 rand_comp_nvars,
                                                 count_comma_separated_items(vars_str),
                                                 p_fmpz, power,
                                                 ctx_initialized ? ctx : NULL,
                                                 output_filename, silent_mode,
                                                 comp_time, omega,
                                                 rand_comp_spec);
        } else {
            run_complexity_analysis(polys_str, vars_str,
                                    p_fmpz, power, ctx_initialized ? ctx : NULL,
                                    output_filename, silent_mode,
                                    comp_time, omega);
        }

    } else if (solve_mode) {
        /* ---- Polynomial system solver ---- */
        if (!silent_mode) {
            int poly_count = count_comma_separated_items(polys_str);
            printf("\nEquations: %d\n", poly_count);
            printf("--------------------------------\n");
        }

        int orig_stdout = -1, orig_stderr = -1;
        int suppress_solver_internal_output = !solve_verbose_mode;

        if (rational_mode) {
            rational_solver_set_realtime_progress(solve_verbose_mode);
            rational_solver_set_internal_trace(solve_verbose_mode);
            rational_solver_set_exact_only(solve_rational_only_mode);

            if (suppress_solver_internal_output) {
                redirect_stdio_to_devnull(&orig_stdout, &orig_stderr);
            }

            rational_solutions = solve_rational_polynomial_system_string(polys_str);

            if (suppress_solver_internal_output) {
                restore_stdio(orig_stdout, orig_stderr);
            }
        } else if (large_prime_mode) {
            large_prime_solver_set_realtime_progress(solve_verbose_mode);

            if (suppress_solver_internal_output) {
                redirect_stdio_to_devnull(&orig_stdout, &orig_stderr);
            }

            large_prime_solutions = solve_large_prime_polynomial_system_string(polys_str, p_fmpz);

            if (suppress_solver_internal_output) {
                restore_stdio(orig_stdout, orig_stderr);
            }
        } else {
            polynomial_solver_set_realtime_progress(solve_verbose_mode);
            polynomial_solver_set_internal_trace(solve_verbose_mode);

            if (suppress_solver_internal_output) {
                redirect_stdio_to_devnull(&orig_stdout, &orig_stderr);
            }

            solutions = solve_polynomial_system_string(polys_str, ctx);

            if (suppress_solver_internal_output) {
                restore_stdio(orig_stdout, orig_stderr);
            }
        }

    } else if (ideal_str) {
        /* ---- Dixon with ideal reduction ---- */
        if (!silent_mode) {
            printf("Task: Dixon + ideal reduction  |  Eliminate: %s\n", vars_str);
            printf("Ideal: %s\n", ideal_str);
            printf("--------------------------------\n");
        }

        int orig_stdout = -1, orig_stderr = -1;
        if (silent_mode) {
            redirect_stdio_to_devnull(&orig_stdout, &orig_stderr);
        }

        result = dixon_with_ideal_reduction_str(polys_str, vars_str, ideal_str, ctx);

        if (silent_mode) {
            restore_stdio(orig_stdout, orig_stderr);
        }

    } else {
        /* ---- Basic resultant ---- */
        if (!silent_mode) {
            int poly_count = count_comma_separated_items(polys_str);
            int var_count  = count_comma_separated_items(vars_str);
            printf("Task: %s  |  Equations: %d  |  Eliminate: %s\n",
                   resultant_method_task_label(resultant_method),
                   poly_count, vars_str);
            if (resultant_method != RESULTANT_METHOD_SUBRES &&
                var_count != poly_count - 1)
                printf("WARNING: resultant mode requires eliminating exactly %d variables "
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
            if (resultant_method == RESULTANT_METHOD_SUBRES) {
                result = compute_subres_resultant_str(polys_str, vars_str, ctx);
            } else {
                result = dixon_str(polys_str, vars_str, ctx);
            }
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
        } else if (large_prime_mode) {
            if (large_prime_solutions) {
                if (!silent_mode) {
                    print_large_prime_solutions(large_prime_solutions);
                }

                if (output_filename) {
                    save_large_prime_solver_result_to_file(output_filename, polys_str,
                                                           p_fmpz, large_prime_solutions,
                                                           cpu_time, wall_time, total_threads);
                    if (!silent_mode)
                        printf("Result saved to: %s\n", output_filename);
                }
                large_prime_solutions_clear(large_prime_solutions);
                free(large_prime_solutions);
            } else {
                if (!silent_mode)
                    fprintf(stderr, "\nError: Large-prime polynomial system solving failed\n");
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
                    if (large_prime_mode) {
                        large_prime_print_roots_from_resultant_string(result, polys_str, vars_str,
                                                                      p_fmpz, fp_append, !silent_mode);
                    } else {
                        const char *cached_root_report = dixon_get_last_root_report();
                        if (cached_root_report && cached_root_report[0] != '\0') {
                            fputs(cached_root_report, fp_append);
                        } else {
                            append_roots_to_file_from_result(result, polys_str, vars_str, ctx, fp_append);
                        }
                    }
                    fclose(fp_append);
                } else if (large_prime_mode && !silent_mode) {
                    large_prime_print_roots_from_resultant_string(result, polys_str, vars_str,
                                                                  p_fmpz, NULL, 1);
                }
            } else if (large_prime_mode) {
                large_prime_print_roots_from_resultant_string(result, polys_str, vars_str,
                                                              p_fmpz, NULL, !silent_mode);
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

    if (!silent_mode) {
        printf("Total - CPU time: %.3f seconds | Wall time: %.3f seconds | Threads: %d\n",
               cpu_time, wall_time, total_threads);
    }

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
    if (!need_free && !rand_generated && vars_str_generated && vars_str)
        free(vars_str);
    free(rand_comp_degrees);
    free(rand_comp_spec);
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
    if (rand_generated) {
        free(polys_str);
        free(vars_str);
        free(allvars_str);
    }
    if (!need_free && !rand_generated && vars_str_override)
        free(vars_str_override);
    if (!need_free && !rand_generated && vars_str_generated && vars_str)
        free(vars_str);
    free(rand_comp_degrees);
    free(rand_comp_spec);
    if (input_filename)  free(input_filename);
    if (output_filename) free(output_filename);
    flint_cleanup();
    return 1;
}
