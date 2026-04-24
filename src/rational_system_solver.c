#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "rational_system_solver.h"
#include <ctype.h>
#include <stdarg.h>

static int rational_solver_realtime_progress_enabled = 0;

#define RATIONAL_REAL_ROOT_PREC 256
#define RATIONAL_REAL_ROOT_DIGITS 20
#define RATIONAL_NUMERIC_NEWTON_TOL 1e-10
#define RATIONAL_NUMERIC_VERIFY_TOL 1e-7
#define RATIONAL_NUMERIC_NEWTON_MAX_ITER 60

typedef struct {
    const char *cursor;
    rational_variable_info_t *vars;
    slong num_vars;
    const double *values;
    int ok;
} rational_numeric_parser_t;

static int rational_verify_numeric_full_solution(char **poly_strings, slong num_polys,
                                                 rational_variable_info_t *vars, slong num_vars,
                                                 const double *values, double tol,
                                                 double *max_residual_out);

void rational_solver_set_realtime_progress(int enabled) {
    rational_solver_realtime_progress_enabled = enabled;
}

static void rational_solver_progress(const char *fmt, ...) {
    if (!rational_solver_realtime_progress_enabled) {
        return;
    }

    va_list args;
    va_start(args, fmt);
    fprintf(stderr, "Progress: ");
    vfprintf(stderr, fmt, args);
    fprintf(stderr, "\n");
    fflush(stderr);
    va_end(args);
}

static char *rational_arb_to_string(const arb_t value, slong digits) {
    char *buffer = NULL;
#ifdef _WIN32
    FILE *mem = tmpfile();
    long length;
    size_t read_size;

    if (!mem) {
        return NULL;
    }

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
    if (!mem) {
        return NULL;
    }

    arb_fprintd(mem, value, digits);
    fclose(mem);
#endif
    return buffer;
}

static int rational_string_is_zeroish(const char *s) {
    if (!s) return 0;
    while (*s && isspace((unsigned char) *s)) s++;
    if (*s == '+') s++;
    while (*s && isspace((unsigned char) *s)) s++;
    int saw_digit = 0;
    while (*s) {
        if (isdigit((unsigned char) *s)) {
            saw_digit = 1;
            if (*s != '0') return 0;
        } else if (*s == '.' || isspace((unsigned char) *s)) {
            /* ok */
        } else {
            return 0;
        }
        s++;
    }
    return saw_digit;
}

static char *rational_arb_to_pretty_string(const arb_t value, slong digits) {
    char *raw = rational_arb_to_string(value, digits);
    if (!raw) {
        return NULL;
    }

    char *pm = strstr(raw, " +/- ");
    if (!pm) {
        return raw;
    }

    if (!arb_is_exact(value) && !rational_string_is_zeroish(pm + 5)) {
        return raw;
    }

    size_t len = (size_t) (pm - raw);
    while (len > 0 && isspace((unsigned char) raw[len - 1])) {
        len--;
    }

    char *pretty = (char *) malloc(len + 1);
    if (!pretty) {
        return raw;
    }
    memcpy(pretty, raw, len);
    pretty[len] = '\0';
    free(raw);
    return pretty;
}

static void rational_print_arb_pretty(const arb_t value, slong digits) {
    char *pretty = rational_arb_to_pretty_string(value, digits);
    if (pretty) {
        printf("%s", pretty);
        free(pretty);
    } else {
        arb_printd(value, digits);
    }
}

static void rational_fprint_arb_pretty(FILE *fp, const arb_t value, slong digits) {
    char *pretty = rational_arb_to_pretty_string(value, digits);
    if (pretty) {
        fprintf(fp, "%s", pretty);
        free(pretty);
    } else {
        arb_fprintd(fp, value, digits);
    }
}

static int rational_append_unique_real_root(arb_t **roots, slong *count, slong *cap, const arb_t value) {
    if (!roots || !count || !cap) {
        return 0;
    }

    for (slong i = 0; i < *count; i++) {
        if (arb_overlaps((*roots)[i], value)) {
            return 1;
        }
    }

    if (*count >= *cap) {
        slong new_cap = (*cap == 0) ? 4 : (*cap * 2);
        arb_t *new_roots = (arb_t *) realloc(*roots, new_cap * sizeof(arb_t));
        if (!new_roots) {
            return 0;
        }

        *roots = new_roots;
        *cap = new_cap;
    }

    arb_init((*roots)[*count]);
    arb_set((*roots)[*count], value);
    (*count)++;
    return 1;
}

static void rational_clear_real_root_array(arb_t *roots, slong count) {
    if (!roots) {
        return;
    }

    for (slong i = 0; i < count; i++) {
        arb_clear(roots[i]);
    }
    free(roots);
}

static void rational_skip_ws(rational_numeric_parser_t *parser) {
    while (parser->cursor && isspace((unsigned char) *parser->cursor)) {
        parser->cursor++;
    }
}

static double rational_parse_numeric_expression(rational_numeric_parser_t *parser);

static double rational_lookup_numeric_variable(rational_numeric_parser_t *parser, const char *name, size_t len) {
    for (slong i = 0; i < parser->num_vars; i++) {
        if (strlen(parser->vars[i].name) == len &&
            strncmp(parser->vars[i].name, name, len) == 0) {
            return parser->values[i];
        }
    }
    parser->ok = 0;
    return 0.0;
}

static double rational_parse_numeric_primary(rational_numeric_parser_t *parser) {
    rational_skip_ws(parser);

    if (*parser->cursor == '\0') {
        parser->ok = 0;
        return 0.0;
    }

    if (*parser->cursor == '(') {
        parser->cursor++;
        double value = rational_parse_numeric_expression(parser);
        rational_skip_ws(parser);
        if (*parser->cursor == ')') {
            parser->cursor++;
        } else {
            parser->ok = 0;
        }
        return value;
    }

    if (*parser->cursor == '+' || *parser->cursor == '-') {
        int sign = (*parser->cursor == '-') ? -1 : 1;
        parser->cursor++;
        return sign * rational_parse_numeric_primary(parser);
    }

    if (isalpha((unsigned char) *parser->cursor) || *parser->cursor == '_') {
        const char *start = parser->cursor;
        while (isalnum((unsigned char) *parser->cursor) || *parser->cursor == '_') {
            parser->cursor++;
        }
        return rational_lookup_numeric_variable(parser, start, (size_t) (parser->cursor - start));
    }

    if (isdigit((unsigned char) *parser->cursor) || *parser->cursor == '.') {
        char *endptr = NULL;
        double numerator = strtod(parser->cursor, &endptr);
        if (endptr == parser->cursor) {
            parser->ok = 0;
            return 0.0;
        }
        parser->cursor = endptr;
        rational_skip_ws(parser);
        if (*parser->cursor == '/') {
            parser->cursor++;
            rational_skip_ws(parser);
            double denominator = strtod(parser->cursor, &endptr);
            if (endptr == parser->cursor || fabs(denominator) < 1e-30) {
                parser->ok = 0;
                return 0.0;
            }
            parser->cursor = endptr;
            numerator /= denominator;
        }
        return numerator;
    }

    parser->ok = 0;
    return 0.0;
}

static double rational_parse_numeric_power(rational_numeric_parser_t *parser) {
    double base = rational_parse_numeric_primary(parser);
    rational_skip_ws(parser);

    while (*parser->cursor == '^') {
        parser->cursor++;
        rational_skip_ws(parser);

        char *endptr = NULL;
        long exponent = strtol(parser->cursor, &endptr, 10);
        if (endptr == parser->cursor) {
            parser->ok = 0;
            return 0.0;
        }
        parser->cursor = endptr;
        base = pow(base, (double) exponent);
        rational_skip_ws(parser);
    }

    return base;
}

static double rational_parse_numeric_term(rational_numeric_parser_t *parser) {
    double value = rational_parse_numeric_power(parser);
    rational_skip_ws(parser);

    while (*parser->cursor == '*' || *parser->cursor == '/') {
        char op = *parser->cursor;
        parser->cursor++;
        double rhs = rational_parse_numeric_power(parser);
        if (op == '*') {
            value *= rhs;
        } else {
            if (fabs(rhs) < 1e-30) {
                parser->ok = 0;
                return 0.0;
            }
            value /= rhs;
        }
        rational_skip_ws(parser);
    }

    return value;
}

static double rational_parse_numeric_expression(rational_numeric_parser_t *parser) {
    double value = rational_parse_numeric_term(parser);
    rational_skip_ws(parser);

    while (*parser->cursor == '+' || *parser->cursor == '-') {
        char op = *parser->cursor;
        parser->cursor++;
        double rhs = rational_parse_numeric_term(parser);
        if (op == '+') {
            value += rhs;
        } else {
            value -= rhs;
        }
        rational_skip_ws(parser);
    }

    return value;
}

static double rational_evaluate_polynomial_double(const char *poly_str,
                                                  rational_variable_info_t *vars,
                                                  slong num_vars,
                                                  const double *values,
                                                  int *ok_out) {
    rational_numeric_parser_t parser;
    parser.cursor = poly_str;
    parser.vars = vars;
    parser.num_vars = num_vars;
    parser.values = values;
    parser.ok = 1;

    double value = rational_parse_numeric_expression(&parser);
    rational_skip_ws(&parser);

    if (!parser.ok || *parser.cursor != '\0') {
        if (ok_out) *ok_out = 0;
        return 0.0;
    }

    if (ok_out) *ok_out = 1;
    return value;
}

static int rational_contains_any_variable(const char *poly_str,
                                          rational_variable_info_t *vars,
                                          slong num_vars) {
    for (slong i = 0; i < num_vars; i++) {
        if (rational_contains_identifier_improved(poly_str, vars[i].name)) {
            return 1;
        }
    }
    return 0;
}

static int rational_real_root_matches_rational(const arb_t root,
                                               const fmpq_t *rational_roots,
                                               slong num_rational_roots,
                                               slong prec) {
    for (slong i = 0; i < num_rational_roots; i++) {
        arb_t rational_as_real;
        arb_init(rational_as_real);
        arb_set_fmpq(rational_as_real, rational_roots[i], prec);
        int match = arb_overlaps(root, rational_as_real);
        arb_clear(rational_as_real);
        if (match) {
            return 1;
        }
    }
    return 0;
}

static double rational_fmpq_to_double(const fmpq_t value) {
    char *buffer = fmpq_get_str(NULL, 10, value);
    double result = buffer ? strtod(buffer, NULL) : 0.0;
    if (buffer) {
        flint_free(buffer);
    }
    return result;
}

static double rational_arb_to_double(const arb_t value) {
    char *buffer = rational_arb_to_string(value, 40);
    double result = buffer ? strtod(buffer, NULL) : 0.0;
    if (buffer) {
        free(buffer);
    }
    return result;
}

static int rational_solution_vectors_close(const double *lhs, const double *rhs, slong n, double tol) {
    for (slong i = 0; i < n; i++) {
        if (fabs(lhs[i] - rhs[i]) > tol) {
            return 0;
        }
    }
    return 1;
}

static int rational_solve_linear_system_double(double *A, double *b, double *x, slong n) {
    for (slong i = 0; i < n; i++) {
        slong pivot = i;
        double pivot_abs = fabs(A[i * n + i]);
        for (slong r = i + 1; r < n; r++) {
            double cand = fabs(A[r * n + i]);
            if (cand > pivot_abs) {
                pivot = r;
                pivot_abs = cand;
            }
        }

        if (pivot_abs < 1e-14) {
            return 0;
        }

        if (pivot != i) {
            for (slong c = i; c < n; c++) {
                double tmp = A[i * n + c];
                A[i * n + c] = A[pivot * n + c];
                A[pivot * n + c] = tmp;
            }
            double tmpb = b[i];
            b[i] = b[pivot];
            b[pivot] = tmpb;
        }

        double diag = A[i * n + i];
        for (slong r = i + 1; r < n; r++) {
            double factor = A[r * n + i] / diag;
            A[r * n + i] = 0.0;
            for (slong c = i + 1; c < n; c++) {
                A[r * n + c] -= factor * A[i * n + c];
            }
            b[r] -= factor * b[i];
        }
    }

    for (slong i = n - 1; i >= 0; i--) {
        double sum = b[i];
        for (slong c = i + 1; c < n; c++) {
            sum -= A[i * n + c] * x[c];
        }
        if (fabs(A[i * n + i]) < 1e-14) {
            return 0;
        }
        x[i] = sum / A[i * n + i];
    }
    return 1;
}

static void rational_set_real_root_summary(rational_solutions_t *sols,
                                           const char *var_name,
                                           const arb_t *roots,
                                           slong num_roots,
                                           slong num_rational_roots_used) {
    if (!sols || !var_name || !roots || num_roots <= 0) {
        return;
    }

    if (sols->real_root_summary) {
        free(sols->real_root_summary);
        sols->real_root_summary = NULL;
    }

    size_t total_len = strlen(var_name) + 192;
    for (slong i = 0; i < num_roots; i++) {
        char *root_str = rational_arb_to_string(roots[i], RATIONAL_REAL_ROOT_DIGITS);
        if (root_str) {
            total_len += strlen(root_str) + 4;
            free(root_str);
        }
    }

    char *summary = (char *) malloc(total_len);
    if (!summary) {
        return;
    }

    if (num_rational_roots_used > 0) {
        snprintf(summary, total_len,
                 "Detected %ld real root(s) for resultant variable %s. Expanded %ld rational branch(es); additional irrational real branches were handled by numerical back-substitution when possible.",
                 num_roots, var_name, num_rational_roots_used);
    } else {
        snprintf(summary, total_len,
                 "Detected %ld real root(s) for resultant variable %s and attempted numerical back-substitution on those branches.",
                 num_roots, var_name);
    }

    strcat(summary, " Real roots ≈ {");
    for (slong i = 0; i < num_roots; i++) {
        char *root_str = rational_arb_to_string(roots[i], RATIONAL_REAL_ROOT_DIGITS);
        if (i > 0) {
            strcat(summary, ", ");
        }
        if (root_str) {
            strcat(summary, root_str);
            free(root_str);
        } else {
            strcat(summary, "?");
        }
    }
    strcat(summary, "}");

    sols->real_root_summary = summary;
}

static char *rational_join_variable_names(rational_variable_info_t *vars, slong num_vars) {
    if (!vars || num_vars <= 0) {
        return NULL;
    }

    size_t total_len = 1;
    for (slong i = 0; i < num_vars; i++) {
        total_len += strlen(vars[i].name) + 2;
    }

    char *result = (char *) malloc(total_len);
    if (!result) {
        return NULL;
    }

    result[0] = '\0';
    for (slong i = 0; i < num_vars; i++) {
        if (i > 0) {
            strcat(result, ", ");
        }
        strcat(result, vars[i].name);
    }

    return result;
}

static char *rational_format_equation_index_list(const slong *indices, slong count) {
    if (!indices || count <= 0) {
        return NULL;
    }

    size_t total_len = 1;
    for (slong i = 0; i < count; i++) {
        total_len += 24;
    }

    char *result = (char *) malloc(total_len);
    if (!result) {
        return NULL;
    }

    result[0] = '\0';
    for (slong i = 0; i < count; i++) {
        char piece[32];
        if (i > 0) {
            strcat(result, ",");
        }
        snprintf(piece, sizeof(piece), "%ld", indices[i] + 1);
        strcat(result, piece);
    }

    return result;
}

static void rational_append_string_entry(char ***items, slong *count, slong *cap, char *entry) {
    if (!entry) {
        return;
    }

    if (*count >= *cap) {
        slong new_cap = (*cap == 0) ? 8 : (*cap * 2);
        char **new_items = (char **) realloc(*items, new_cap * sizeof(char *));
        if (!new_items) {
            free(entry);
            return;
        }
        *items = new_items;
        *cap = new_cap;
    }

    (*items)[*count] = entry;
    (*count)++;
}

static void rational_store_resultant_step(rational_solutions_t *sols, char *buffer) {
    if (!sols || !buffer) {
        free(buffer);
        return;
    }

    rational_append_string_entry(&sols->resultant_steps,
                                &sols->num_resultant_steps,
                                &sols->resultant_steps_cap,
                                buffer);
}

static int rational_append_real_solution_values(rational_solutions_t *sols,
                                                const double *values,
                                                slong num_values) {
    if (!sols || !values || num_values != sols->num_variables) {
        return 0;
    }

    slong new_count = sols->num_real_solution_sets + 1;
    arb_t ***new_sets = (arb_t ***) realloc(sols->real_solution_sets, new_count * sizeof(arb_t **));
    slong *new_counts = (slong *) realloc(sols->real_solutions_per_var, new_count * sols->num_variables * sizeof(slong));
    arb_t **new_residuals = (arb_t **) realloc(sols->real_solution_residuals, new_count * sizeof(arb_t *));
    if (!new_sets || !new_counts || !new_residuals) {
        if (new_sets) sols->real_solution_sets = new_sets;
        if (new_counts) sols->real_solutions_per_var = new_counts;
        if (new_residuals) sols->real_solution_residuals = new_residuals;
        return 0;
    }

    sols->real_solution_sets = new_sets;
    sols->real_solutions_per_var = new_counts;
    sols->real_solution_residuals = new_residuals;

    slong set_idx = sols->num_real_solution_sets;
    sols->real_solution_sets[set_idx] = (arb_t **) calloc(sols->num_variables, sizeof(arb_t *));
    sols->real_solution_residuals[set_idx] = NULL;
    if (!sols->real_solution_sets[set_idx]) {
        return 0;
    }

    for (slong var = 0; var < sols->num_variables; var++) {
        sols->real_solution_sets[set_idx][var] = (arb_t *) malloc(sizeof(arb_t));
        if (!sols->real_solution_sets[set_idx][var]) {
            return 0;
        }
        arb_init(sols->real_solution_sets[set_idx][var][0]);
        arb_set_d(sols->real_solution_sets[set_idx][var][0], values[var]);
        sols->real_solutions_per_var[set_idx * sols->num_variables + var] = 1;
    }

    sols->num_real_solution_sets = new_count;
    return 1;
}

static void rational_clear_residual_array(arb_t **residuals, slong num_sets, slong num_equations) {
    if (!residuals) {
        return;
    }
    for (slong set = 0; set < num_sets; set++) {
        if (residuals[set]) {
            for (slong eq = 0; eq < num_equations; eq++) {
                arb_clear(residuals[set][eq]);
            }
            free(residuals[set]);
        }
    }
    free(residuals);
}

static void rational_add_resultant_step(rational_solutions_t *sols, const char *fmt, ...) {
    if (!sols) {
        return;
    }

    va_list args;
    va_list args_copy;
    va_start(args, fmt);
    va_copy(args_copy, args);
    int needed = vsnprintf(NULL, 0, fmt, args_copy);
    va_end(args_copy);

    if (needed < 0) {
        va_end(args);
        return;
    }

    char *buffer = (char *) malloc((size_t) needed + 1);
    if (!buffer) {
        va_end(args);
        return;
    }

    vsnprintf(buffer, (size_t) needed + 1, fmt, args);
    va_end(args);

    rational_store_resultant_step(sols, buffer);
    rational_solver_progress("%s", buffer);
}

static char *rational_format_solution_set_line_from_arrays(char **variable_names, slong num_variables,
                                                           fmpq_t **solution_set, slong *solutions_per_var) {
    if (!variable_names || !solution_set || !solutions_per_var || num_variables <= 0) {
        return NULL;
    }

    size_t total_len = 1;
    for (slong var = 0; var < num_variables; var++) {
        total_len += strlen(variable_names[var]) + 16;
        slong num_sols = solutions_per_var[var];
        if (num_sols == 0) {
            total_len += strlen("no solution");
        } else {
            for (slong sol = 0; sol < num_sols; sol++) {
                char *sol_str = fmpq_get_str(NULL, 10, solution_set[var][sol]);
                if (sol_str) {
                    total_len += strlen(sol_str) + 2;
                    flint_free(sol_str);
                }
            }
        }
    }

    char *line = (char *) malloc(total_len);
    if (!line) {
        return NULL;
    }

    line[0] = '\0';
    for (slong var = 0; var < num_variables; var++) {
        slong num_sols = solutions_per_var[var];

        if (var > 0) {
            strcat(line, ", ");
        }
        strcat(line, variable_names[var]);
        strcat(line, " = ");

        if (num_sols == 0) {
            strcat(line, "no solution");
        } else if (num_sols == 1) {
            char *sol_str = fmpq_get_str(NULL, 10, solution_set[var][0]);
            if (sol_str) {
                strcat(line, sol_str);
                flint_free(sol_str);
            }
        } else {
            strcat(line, "{");
            for (slong sol = 0; sol < num_sols; sol++) {
                char *sol_str = fmpq_get_str(NULL, 10, solution_set[var][sol]);
                if (sol > 0) {
                    strcat(line, ", ");
                }
                if (sol_str) {
                    strcat(line, sol_str);
                    flint_free(sol_str);
                }
            }
            strcat(line, "}");
        }
    }

    return line;
}

static void rational_add_candidate_solution(rational_solutions_t *sols, char *line, int pass) {
    if (!sols || !line) {
        return;
    }

    if (sols->num_candidate_solution_lines >= sols->candidate_solution_lines_cap) {
        slong new_cap = (sols->candidate_solution_lines_cap == 0)
                        ? 8 : (sols->candidate_solution_lines_cap * 2);
        char **new_lines = (char **) realloc(sols->candidate_solution_lines,
                                             new_cap * sizeof(char *));
        int *new_pass = (int *) realloc(sols->candidate_solution_pass,
                                        new_cap * sizeof(int));
        if (!new_lines || !new_pass) {
            free(line);
            if (new_lines) sols->candidate_solution_lines = new_lines;
            if (new_pass) sols->candidate_solution_pass = new_pass;
            return;
        }
        sols->candidate_solution_lines = new_lines;
        sols->candidate_solution_pass = new_pass;
        sols->candidate_solution_lines_cap = new_cap;
    }

    sols->candidate_solution_lines[sols->num_candidate_solution_lines] = line;
    sols->candidate_solution_pass[sols->num_candidate_solution_lines] = pass;
    sols->num_candidate_solution_lines++;
}

static void rational_merge_solver_logs(rational_solutions_t *dst, const rational_solutions_t *src) {
    if (!dst || !src) {
        return;
    }

    for (slong i = 0; i < src->num_resultant_steps; i++) {
        rational_store_resultant_step(dst, strdup(src->resultant_steps[i]));
    }
}

int rational_contains_variable(const char *poly_str, const char *var_name) {
    const char *ptr = poly_str;
    size_t var_len = strlen(var_name);
    
    while ((ptr = strstr(ptr, var_name)) != NULL) {
        if ((ptr == poly_str || !isalnum(*(ptr-1))) && 
            !isalnum(*(ptr + var_len))) {
            return 1;
        }
        ptr++;
    }
    return 0;
}

int rational_contains_identifier_improved(const char *poly_str, const char *identifier) {
    const char *ptr = poly_str;
    size_t id_len = strlen(identifier);
    
    while ((ptr = strstr(ptr, identifier)) != NULL) {
        if ((ptr == poly_str || !isalnum(*(ptr-1))) && 
            !isalnum(*(ptr + id_len))) {
            return 1;
        }
        ptr++;
    }
    return 0;
}

char** rational_extract_variables_improved(char **poly_strings, slong num_polys, slong *num_vars_out) {
    char **temp_vars = (char**) malloc(100 * sizeof(char*));
    slong temp_count = 0;
    
    printf("=== Improved Rational Variable Extraction ===\n");
    
    for (slong poly_idx = 0; poly_idx < num_polys; poly_idx++) {
        const char *poly = poly_strings[poly_idx];
        size_t len = strlen(poly);
        
        for (size_t i = 0; i < len; i++) {
            if (isalpha(poly[i]) || poly[i] == '_') {
                size_t start = i;
                while (i < len && (isalnum(poly[i]) || poly[i] == '_')) {
                    i++;
                }
                
                size_t id_len = i - start;
                char *identifier = (char*) malloc(id_len + 1);
                strncpy(identifier, poly + start, id_len);
                identifier[id_len] = '\0';
                
                if (isdigit(identifier[0])) {
                    free(identifier);
                    continue;
                }
                
                int already_exists = 0;
                for (slong j = 0; j < temp_count; j++) {
                    if (strcmp(temp_vars[j], identifier) == 0) {
                        already_exists = 1;
                        break;
                    }
                }
                
                if (!already_exists) {
                    temp_vars[temp_count] = strdup(identifier);
                    temp_count++;
                    printf("Found variable: %s\n", identifier);
                }
                
                free(identifier);
                i--;
            }
        }
    }
    
    char **final_vars = NULL;
    if (temp_count > 0) {
        final_vars = (char**) malloc(temp_count * sizeof(char*));
        for (slong i = 0; i < temp_count; i++) {
            final_vars[i] = temp_vars[i];
        }
    }
    
    free(temp_vars);
    *num_vars_out = temp_count;
    
    printf("Total variables found: %ld\n", temp_count);
    for (slong i = 0; i < temp_count; i++) {
        printf("  Variable %ld: %s\n", i, final_vars[i]);
    }
    
    return final_vars;
}

slong rational_get_variable_max_degree_in_polynomial(const char *poly_str, const char *var_name) {
    const char *ptr = poly_str;
    slong max_degree = 0;
    size_t var_len = strlen(var_name);
    
    while ((ptr = strstr(ptr, var_name)) != NULL) {
        if ((ptr == poly_str || !isalnum(*(ptr-1))) && 
            !isalnum(*(ptr + var_len))) {
            
            const char *exp_ptr = ptr + var_len;
            while (isspace(*exp_ptr)) exp_ptr++;
            
            slong degree = 1;
            if (*exp_ptr == '^') {
                exp_ptr++;
                while (isspace(*exp_ptr)) exp_ptr++;
                if (isdigit(*exp_ptr)) {
                    degree = strtol(exp_ptr, NULL, 10);
                }
            }
            
            if (degree > max_degree) {
                max_degree = degree;
            }
        }
        ptr++;
    }
    
    return max_degree;
}

static int rational_compare_variables_by_degree(const void *a, const void *b) {
    const rational_variable_info_t *va = (const rational_variable_info_t*)a;
    const rational_variable_info_t *vb = (const rational_variable_info_t*)b;
    
    if (va->max_degree < vb->max_degree) return -1;
    if (va->max_degree > vb->max_degree) return 1;
    return strcmp(va->name, vb->name);
}

void rational_generate_equation_combinations(slong num_equations, slong target_equations, 
                                             rational_equation_combination_t **combinations, slong *num_combinations) {
    slong total_combinations = 1;
    for (slong i = 0; i < target_equations; i++) {
        total_combinations = total_combinations * (num_equations - i) / (i + 1);
    }
    
    *num_combinations = total_combinations;
    *combinations = (rational_equation_combination_t*) malloc(total_combinations * sizeof(rational_equation_combination_t));
    
    slong comb_idx = 0;
    slong *current_combination = (slong*) malloc(target_equations * sizeof(slong));
    
    for (slong i = 0; i < target_equations; i++) {
        current_combination[i] = i;
    }
    
    while (comb_idx < total_combinations) {
        (*combinations)[comb_idx].equation_indices = (slong*) malloc(target_equations * sizeof(slong));
        (*combinations)[comb_idx].num_equations = target_equations;
        (*combinations)[comb_idx].tried = 0;
        (*combinations)[comb_idx].success = 0;
        
        for (slong i = 0; i < target_equations; i++) {
            (*combinations)[comb_idx].equation_indices[i] = current_combination[i];
        }
        comb_idx++;
        
        slong k = target_equations - 1;
        while (k >= 0 && current_combination[k] == num_equations - target_equations + k) {
            k--;
        }
        
        if (k < 0) break;
        
        current_combination[k]++;
        for (slong i = k + 1; i < target_equations; i++) {
            current_combination[i] = current_combination[i-1] + 1;
        }
    }
    
    free(current_combination);
}

void rational_free_equation_combinations(rational_equation_combination_t *combinations, slong num_combinations) {
    if (combinations) {
        for (slong i = 0; i < num_combinations; i++) {
            if (combinations[i].equation_indices) {
                free(combinations[i].equation_indices);
            }
        }
        free(combinations);
    }
}

void rational_solutions_init(rational_solutions_t *sols, slong num_vars) {
    sols->num_variables = num_vars;
    sols->is_valid = 0;
    sols->has_no_solutions = 0;
    sols->error_message = NULL;
    sols->num_solution_sets = 0;
    sols->solution_residuals = NULL;
    sols->real_solution_sets = NULL;
    sols->real_solutions_per_var = NULL;
    sols->num_real_solution_sets = 0;
    sols->real_solution_residuals = NULL;
    sols->real_solution_precision = RATIONAL_REAL_ROOT_PREC;
    sols->real_root_summary = NULL;
    sols->num_equations = 0;
    sols->variable_order = NULL;
    sols->elimination_summary = NULL;
    sols->total_combinations = 0;
    sols->successful_combinations = 0;
    sols->num_base_solutions = 0;
    sols->checked_solution_sets = -1;
    sols->verified_solution_sets = -1;
    sols->resultant_steps = NULL;
    sols->num_resultant_steps = 0;
    sols->resultant_steps_cap = 0;
    sols->candidate_solution_lines = NULL;
    sols->candidate_solution_pass = NULL;
    sols->num_candidate_solution_lines = 0;
    sols->candidate_solution_lines_cap = 0;
    
    if (num_vars > 0) {
        sols->variable_names = (char**) calloc(num_vars, sizeof(char*));
        sols->solution_sets = NULL;
        sols->solutions_per_var = NULL;
    } else {
        sols->variable_names = NULL;
        sols->solution_sets = NULL;
        sols->solutions_per_var = NULL;
    }
}

void rational_solutions_clear(rational_solutions_t *sols) {
    if (sols->variable_names) {
        for (slong i = 0; i < sols->num_variables; i++) {
            if (sols->variable_names[i]) {
                free(sols->variable_names[i]);
            }
        }
        free(sols->variable_names);
    }
    
    if (sols->solution_sets) {
        for (slong set = 0; set < sols->num_solution_sets; set++) {
            if (sols->solution_sets[set]) {
                for (slong var = 0; var < sols->num_variables; var++) {
                    if (sols->solution_sets[set][var]) {
                        for (slong sol = 0; sol < sols->solutions_per_var[set * sols->num_variables + var]; sol++) {
                            fmpq_clear(sols->solution_sets[set][var][sol]);
                        }
                        free(sols->solution_sets[set][var]);
                    }
                }
                free(sols->solution_sets[set]);
            }
        }
        free(sols->solution_sets);
    }
    
    if (sols->solutions_per_var) {
        free(sols->solutions_per_var);
    }

    rational_clear_residual_array(sols->solution_residuals, sols->num_solution_sets, sols->num_equations);

    if (sols->real_solution_sets) {
        for (slong set = 0; set < sols->num_real_solution_sets; set++) {
            if (sols->real_solution_sets[set]) {
                for (slong var = 0; var < sols->num_variables; var++) {
                    if (sols->real_solution_sets[set][var]) {
                        for (slong sol = 0; sol < sols->real_solutions_per_var[set * sols->num_variables + var]; sol++) {
                            arb_clear(sols->real_solution_sets[set][var][sol]);
                        }
                        free(sols->real_solution_sets[set][var]);
                    }
                }
                free(sols->real_solution_sets[set]);
            }
        }
        free(sols->real_solution_sets);
    }

    if (sols->real_solutions_per_var) {
        free(sols->real_solutions_per_var);
    }

    rational_clear_residual_array(sols->real_solution_residuals, sols->num_real_solution_sets, sols->num_equations);
    
    if (sols->error_message) {
        free(sols->error_message);
    }

    if (sols->variable_order) {
        free(sols->variable_order);
    }

    if (sols->elimination_summary) {
        free(sols->elimination_summary);
    }

    if (sols->real_root_summary) {
        free(sols->real_root_summary);
    }

    if (sols->resultant_steps) {
        for (slong i = 0; i < sols->num_resultant_steps; i++) {
            free(sols->resultant_steps[i]);
        }
        free(sols->resultant_steps);
    }

    if (sols->candidate_solution_lines) {
        for (slong i = 0; i < sols->num_candidate_solution_lines; i++) {
            free(sols->candidate_solution_lines[i]);
        }
        free(sols->candidate_solution_lines);
    }

    if (sols->candidate_solution_pass) {
        free(sols->candidate_solution_pass);
    }
    
    memset(sols, 0, sizeof(rational_solutions_t));
}

static void rational_compute_exact_solution_residuals(rational_solutions_t *sols,
                                                      char **poly_strings,
                                                      rational_variable_info_t *vars) {
    if (!sols || !poly_strings || !vars || sols->num_solution_sets <= 0 || sols->num_equations <= 0) {
        return;
    }

    rational_clear_residual_array(sols->solution_residuals, sols->num_solution_sets, sols->num_equations);
    sols->solution_residuals = (arb_t **) calloc(sols->num_solution_sets, sizeof(arb_t *));
    if (!sols->solution_residuals) {
        return;
    }

    double *values = (double *) calloc(sols->num_variables, sizeof(double));
    if (!values) {
        free(sols->solution_residuals);
        sols->solution_residuals = NULL;
        return;
    }

    for (slong set = 0; set < sols->num_solution_sets; set++) {
        int complete = 1;
        for (slong var = 0; var < sols->num_variables; var++) {
            if (sols->solutions_per_var[set * sols->num_variables + var] <= 0) {
                complete = 0;
                break;
            }
            values[var] = rational_fmpq_to_double(sols->solution_sets[set][var][0]);
        }
        if (!complete) {
            continue;
        }

        sols->solution_residuals[set] = (arb_t *) malloc(sols->num_equations * sizeof(arb_t));
        if (!sols->solution_residuals[set]) {
            continue;
        }

        for (slong eq = 0; eq < sols->num_equations; eq++) {
            arb_init(sols->solution_residuals[set][eq]);
            int ok = 0;
            double residual = rational_evaluate_polynomial_double(poly_strings[eq], vars, sols->num_variables, values, &ok);
            arb_set_d(sols->solution_residuals[set][eq], ok ? fabs(residual) : INFINITY);
        }
    }

    free(values);
}

static void rational_compute_real_solution_residuals(rational_solutions_t *sols,
                                                     char **poly_strings,
                                                     rational_variable_info_t *vars) {
    if (!sols || !poly_strings || !vars || sols->num_real_solution_sets <= 0 || sols->num_equations <= 0) {
        return;
    }

    rational_clear_residual_array(sols->real_solution_residuals, sols->num_real_solution_sets, sols->num_equations);
    sols->real_solution_residuals = (arb_t **) calloc(sols->num_real_solution_sets, sizeof(arb_t *));
    if (!sols->real_solution_residuals) {
        return;
    }

    double *values = (double *) calloc(sols->num_variables, sizeof(double));
    if (!values) {
        free(sols->real_solution_residuals);
        sols->real_solution_residuals = NULL;
        return;
    }

    for (slong set = 0; set < sols->num_real_solution_sets; set++) {
        int complete = 1;
        for (slong var = 0; var < sols->num_variables; var++) {
            if (sols->real_solutions_per_var[set * sols->num_variables + var] <= 0) {
                complete = 0;
                break;
            }
            values[var] = rational_arb_to_double(sols->real_solution_sets[set][var][0]);
        }
        if (!complete) {
            continue;
        }

        sols->real_solution_residuals[set] = (arb_t *) malloc(sols->num_equations * sizeof(arb_t));
        if (!sols->real_solution_residuals[set]) {
            continue;
        }

        for (slong eq = 0; eq < sols->num_equations; eq++) {
            arb_init(sols->real_solution_residuals[set][eq]);
            int ok = 0;
            double residual = rational_evaluate_polynomial_double(poly_strings[eq], vars, sols->num_variables, values, &ok);
            arb_set_d(sols->real_solution_residuals[set][eq], ok ? fabs(residual) : INFINITY);
        }
    }

    free(values);
}

static void rational_print_residuals(const arb_t *residuals, slong num_equations) {
    if (!residuals || num_equations <= 0) {
        return;
    }

        printf("    residuals: [");
    for (slong eq = 0; eq < num_equations; eq++) {
        if (eq > 0) {
            printf(", ");
        }
        printf("eq%ld=", eq + 1);
        rational_print_arb_pretty(residuals[eq], 8);
    }
    printf("]\n");
}

static void rational_filter_real_solutions_by_verification(rational_solutions_t *sols,
                                                           char **poly_strings,
                                                           rational_variable_info_t *sorted_vars) {
    if (!sols || !poly_strings || !sorted_vars || sols->num_real_solution_sets == 0) {
        return;
    }

    int *valid = (int *) calloc(sols->num_real_solution_sets, sizeof(int));
    if (!valid) {
        return;
    }

    double *values = (double *) calloc(sols->num_variables, sizeof(double));
    if (!values) {
        free(valid);
        return;
    }

    slong keep_count = 0;
    for (slong set = 0; set < sols->num_real_solution_sets; set++) {
        int complete = 1;
        for (slong var = 0; var < sols->num_variables; var++) {
            if (sols->real_solutions_per_var[set * sols->num_variables + var] <= 0) {
                complete = 0;
                break;
            }
            values[var] = rational_arb_to_double(sols->real_solution_sets[set][var][0]);
        }

        double max_residual = 0.0;
        if (complete &&
            rational_verify_numeric_full_solution(poly_strings, sols->num_equations,
                                                  sorted_vars, sols->num_variables,
                                                  values, RATIONAL_NUMERIC_VERIFY_TOL,
                                                  &max_residual)) {
            valid[set] = 1;
            keep_count++;
        }
    }

    if (keep_count == sols->num_real_solution_sets) {
        free(valid);
        free(values);
        return;
    }

    if (keep_count == 0) {
        if (sols->real_solution_sets) {
            for (slong set = 0; set < sols->num_real_solution_sets; set++) {
                if (sols->real_solution_sets[set]) {
                    for (slong var = 0; var < sols->num_variables; var++) {
                        if (sols->real_solution_sets[set][var]) {
                            for (slong sol = 0; sol < sols->real_solutions_per_var[set * sols->num_variables + var]; sol++) {
                                arb_clear(sols->real_solution_sets[set][var][sol]);
                            }
                            free(sols->real_solution_sets[set][var]);
                        }
                    }
                    free(sols->real_solution_sets[set]);
                }
            }
            free(sols->real_solution_sets);
            sols->real_solution_sets = NULL;
        }
        free(sols->real_solutions_per_var);
        sols->real_solutions_per_var = NULL;
        rational_clear_residual_array(sols->real_solution_residuals, sols->num_real_solution_sets, sols->num_equations);
        sols->real_solution_residuals = NULL;
        sols->num_real_solution_sets = 0;
        free(valid);
        free(values);
        return;
    }

    arb_t ***new_sets = (arb_t ***) calloc(keep_count, sizeof(arb_t **));
    slong *new_counts = (slong *) calloc(keep_count * sols->num_variables, sizeof(slong));
    if (!new_sets || !new_counts) {
        free(new_sets);
        free(new_counts);
        free(valid);
        free(values);
        return;
    }

    slong out = 0;
    for (slong set = 0; set < sols->num_real_solution_sets; set++) {
        if (valid[set]) {
            new_sets[out] = sols->real_solution_sets[set];
            for (slong var = 0; var < sols->num_variables; var++) {
                new_counts[out * sols->num_variables + var] =
                    sols->real_solutions_per_var[set * sols->num_variables + var];
            }
            sols->real_solution_sets[set] = NULL;
            out++;
        } else if (sols->real_solution_sets[set]) {
            for (slong var = 0; var < sols->num_variables; var++) {
                if (sols->real_solution_sets[set][var]) {
                    for (slong sol = 0; sol < sols->real_solutions_per_var[set * sols->num_variables + var]; sol++) {
                        arb_clear(sols->real_solution_sets[set][var][sol]);
                    }
                    free(sols->real_solution_sets[set][var]);
                }
            }
            free(sols->real_solution_sets[set]);
        }
    }

    free(sols->real_solution_sets);
    free(sols->real_solutions_per_var);
    rational_clear_residual_array(sols->real_solution_residuals, sols->num_real_solution_sets, sols->num_equations);
    sols->real_solution_residuals = NULL;

    sols->real_solution_sets = new_sets;
    sols->real_solutions_per_var = new_counts;
    sols->num_real_solution_sets = keep_count;

    free(valid);
    free(values);
}

void print_rational_solutions(const rational_solutions_t *sols) {
    if (!sols) {
        printf("No solution data available.\n");
        return;
    }

    if (sols->num_resultant_steps > 0) {
        for (slong i = 0; i < sols->num_resultant_steps; i++) {
            printf("Progress: %s\n", sols->resultant_steps[i]);
        }
        printf("\n");
    }

    if (!sols->is_valid) {
        printf("Solve failed");
        if (sols->error_message) {
            printf(": %s", sols->error_message);
        }
        printf("\n");
        return;
    }

    if (sols->has_no_solutions == -1) {
        printf("System has positive dimension.\n");
        return;
    }

    if (sols->num_variables == 0) {
        printf("No variables.\n");
        return;
    }

    if (sols->has_no_solutions == 1 &&
        sols->num_solution_sets == 0 &&
        sols->num_real_solution_sets == 0) {
        if (sols->real_root_summary) {
            printf("No exact rational solution sets were assembled.\n");
            printf("%s\n", sols->real_root_summary);
        } else {
            printf("No solutions over the rational numbers.\n");
        }
        return;
    }

    if (sols->num_candidate_solution_lines > 0 && sols->num_solution_sets > 0) {
        printf("Candidate sets:\n");
        for (slong i = 0; i < sols->num_candidate_solution_lines; i++) {
            printf("  [%ld] %s [%s]\n",
                   i + 1,
                   sols->candidate_solution_lines[i],
                   sols->candidate_solution_pass[i] ? "PASS" : "FAIL");
        }
    }

    if (sols->num_solution_sets > 0) {
        printf("Found %ld exact rational solution set(s):\n", sols->num_solution_sets);
        for (slong set = 0; set < sols->num_solution_sets; set++) {
            printf("  [%ld] ", set + 1);
            for (slong var = 0; var < sols->num_variables; var++) {
                slong num_sols = sols->solutions_per_var[set * sols->num_variables + var];

                if (var > 0) {
                    printf(", ");
                }

                printf("%s = ", sols->variable_names[var]);
                if (num_sols == 0) {
                    printf("no solution");
                } else if (num_sols == 1) {
                    fmpq_print(sols->solution_sets[set][var][0]);
                } else {
                    printf("{");
                    for (slong sol = 0; sol < num_sols; sol++) {
                        if (sol > 0) printf(", ");
                        fmpq_print(sols->solution_sets[set][var][sol]);
                    }
                    printf("}");
                }
            }
            printf("\n");
            rational_print_residuals(sols->solution_residuals ? sols->solution_residuals[set] : NULL,
                                     sols->num_equations);
        }
    }

    if (sols->num_real_solution_sets > 0) {
        printf("Found %ld approximate real solution set(s):\n", sols->num_real_solution_sets);
        for (slong set = 0; set < sols->num_real_solution_sets; set++) {
            printf("  [%ld] ", set + 1);
            for (slong var = 0; var < sols->num_variables; var++) {
                slong num_sols = sols->real_solutions_per_var[set * sols->num_variables + var];

                if (var > 0) {
                    printf(", ");
                }

                printf("%s = ", sols->variable_names[var]);
                if (num_sols == 0) {
                    printf("no solution");
                } else if (num_sols == 1) {
                    rational_print_arb_pretty(sols->real_solution_sets[set][var][0], RATIONAL_REAL_ROOT_DIGITS);
                } else {
                    printf("{");
                    for (slong sol = 0; sol < num_sols; sol++) {
                        if (sol > 0) printf(", ");
                        rational_print_arb_pretty(sols->real_solution_sets[set][var][sol], RATIONAL_REAL_ROOT_DIGITS);
                    }
                    printf("}");
                }
            }
            printf("\n");
            rational_print_residuals(sols->real_solution_residuals ? sols->real_solution_residuals[set] : NULL,
                                     sols->num_equations);
        }
    }

    if (sols->real_root_summary) {
        printf("%s\n", sols->real_root_summary);
    }
}

int rational_extract_and_sort_variables(char **poly_strings, slong num_polys,
                                         rational_variable_info_t **sorted_vars_out, slong *num_vars_out) {
    
    printf("=== Analyzing Rational Polynomial Variables (improved) ===\n");
    
    char **all_vars = rational_extract_variables_improved(poly_strings, num_polys, num_vars_out);
    
    if (*num_vars_out == 0) {
        *sorted_vars_out = NULL;
        return 1;
    }
    
    for (slong i = 0; i < num_polys; i++) {
        slong total_degree = 0;
        
        for (slong j = 0; j < *num_vars_out; j++) {
            slong var_degree = rational_get_variable_max_degree_in_polynomial(poly_strings[i], all_vars[j]);
            if (var_degree > total_degree) {
                total_degree = var_degree;
            }
        }
        
        printf("Polynomial %ld: estimated degree = %ld\n", i + 1, total_degree);
    }
    
    printf("Discovered variables (%ld): ", *num_vars_out);
    for (slong i = 0; i < *num_vars_out; i++) {
        if (i > 0) printf(", ");
        printf("%s", all_vars[i]);
    }
    printf("\n");
    
    rational_variable_info_t *var_info = (rational_variable_info_t*) malloc(*num_vars_out * sizeof(rational_variable_info_t));
    
    for (slong i = 0; i < *num_vars_out; i++) {
        var_info[i].name = strdup(all_vars[i]);
        var_info[i].index = i;
        var_info[i].max_degree = 0;
        
        for (slong j = 0; j < num_polys; j++) {
            slong degree = rational_get_variable_max_degree_in_polynomial(poly_strings[j], all_vars[i]);
            if (degree > var_info[i].max_degree) {
                var_info[i].max_degree = degree;
            }
        }
        
        free(all_vars[i]);
    }
    
    free(all_vars);
    
    qsort(var_info, *num_vars_out, sizeof(rational_variable_info_t), rational_compare_variables_by_degree);
    
    printf("\nVariables sorted by maximum degree:\n");
    for (slong i = 0; i < *num_vars_out; i++) {
        printf("  %s: max degree %ld\n", var_info[i].name, var_info[i].max_degree);
    }
    
    *sorted_vars_out = var_info;
    return 1;
}

static int rational_parse_string_to_fmpq_poly(const char *poly_str, const char *var_name,
                                               fmpq_poly_t poly) {
    fmpq_poly_init(poly);
    fmpq_poly_zero(poly);
    
    if (strcmp(poly_str, "0") == 0 || strcmp(poly_str, "") == 0) {
        return 1;
    }
    
    const char *ptr = poly_str;
    fmpq_t coeff, exp_val;
    fmpq_init(coeff);
    fmpq_init(exp_val);
    
    while (*ptr) {
        while (isspace(*ptr)) ptr++;
        
        if (*ptr == '\0') break;
        
        // Skip any parentheses
        if (*ptr == '(' || *ptr == ')') {
            ptr++;
            continue;
        }
        
        int sign = 1;
        if (*ptr == '+') {
            ptr++;
            while (isspace(*ptr)) ptr++;
        } else if (*ptr == '-') {
            sign = -1;
            ptr++;
            while (isspace(*ptr)) ptr++;
        }
        
        fmpq_t term_coeff;
        fmpq_init(term_coeff);
        fmpq_set_si(term_coeff, 1, 1);
        
        slong term_degree = 0;
        
        if (isdigit(*ptr) || *ptr == '.') {
            const char *num_start = ptr;
            while (isdigit(*ptr) || *ptr == '.' || *ptr == '/') ptr++;
            size_t num_len = ptr - num_start;
            char *num_str = (char*) malloc(num_len + 1);
            strncpy(num_str, num_start, num_len);
            num_str[num_len] = '\0';
            
            char *slash = strchr(num_str, '/');
            if (slash) {
                *slash = '\0';
                fmpz_t num, den;
                fmpz_init(num);
                fmpz_init(den);
                fmpz_set_str(num, num_str, 10);
                fmpz_set_str(den, slash + 1, 10);
                fmpq_set_fmpz_frac(term_coeff, num, den);
                fmpz_clear(num);
                fmpz_clear(den);
            } else {
                fmpz_t num;
                fmpz_init(num);
                fmpz_set_str(num, num_str, 10);
                fmpq_set_fmpz(term_coeff, num);
                fmpz_clear(num);
            }
            free(num_str);
        }
        
        while (isspace(*ptr)) ptr++;
        
        if (*ptr == '*') {
            ptr++;
            while (isspace(*ptr)) ptr++;
        }
        
        if (strncmp(ptr, var_name, strlen(var_name)) == 0) {
            ptr += strlen(var_name);
            term_degree = 1;
            
            while (isspace(*ptr)) ptr++;
            if (*ptr == '^') {
                ptr++;
                while (isspace(*ptr)) ptr++;
                if (isdigit(*ptr)) {
                    term_degree = strtol(ptr, (char**)&ptr, 10);
                }
            }
        } else if (*ptr == '^') {
            // Skip constant exponentiation: const^exp is still a const
            ptr++;
            while (isspace(*ptr)) ptr++;
            if (isdigit(*ptr)) {
                while (isdigit(*ptr)) ptr++;
            }
        }
        
        if (sign < 0) {
            fmpq_neg(term_coeff, term_coeff);
        }
        
        // Add to the existing coefficient instead of replacing
        fmpq_t existing_coeff;
        fmpq_init(existing_coeff);
        fmpq_poly_get_coeff_fmpq(existing_coeff, poly, term_degree);
        fmpq_add(existing_coeff, existing_coeff, term_coeff);
        fmpq_poly_set_coeff_fmpq(poly, term_degree, existing_coeff);
        fmpq_clear(existing_coeff);
        
        fmpq_clear(term_coeff);
    }
    
    fmpq_clear(coeff);
    fmpq_clear(exp_val);
    
    return 1;
}

fmpq_t* rational_solve_univariate_equation_all_roots(const char *poly_str, const char *var_name,
                                                      slong *num_roots_out) {
    printf("Solving univariate equation: %s (variable: %s)\n", poly_str, var_name);
    
    fmpq_poly_t fmpq_poly;
    rational_parse_string_to_fmpq_poly(poly_str, var_name, fmpq_poly);
    
    fmpq_roots_t roots;
    fmpq_roots_init(&roots);
    slong num_roots = fmpq_poly_roots(&roots, fmpq_poly, 0);
    
    fmpq_t *result_roots = NULL;
    if (num_roots > 0) {
        result_roots = (fmpq_t*) malloc(roots.num_roots * sizeof(fmpq_t));
        for (slong i = 0; i < roots.num_roots; i++) {
            fmpq_init(result_roots[i]);
            fmpq_set(result_roots[i], roots.roots[i]);
        }
        *num_roots_out = roots.num_roots;
        
        printf("Found %ld rational roots:\n", roots.num_roots);
        for (slong i = 0; i < roots.num_roots; i++) {
            printf("  Root %ld: ", i + 1);
            fmpq_print(roots.roots[i]);
            printf("\n");
        }
    } else {
        *num_roots_out = 0;
        printf("No rational roots found.\n");
    }
    
    fmpq_roots_clear(&roots);
    fmpq_poly_clear(fmpq_poly);
    
    return result_roots;
}

arb_t* rational_solve_univariate_equation_all_real_roots(const char *poly_str, const char *var_name,
                                                          slong *num_roots_out, slong prec) {
    printf("Solving univariate equation over R: %s (variable: %s)\n", poly_str, var_name);

    fmpq_poly_t fmpq_poly;
    rational_parse_string_to_fmpq_poly(poly_str, var_name, fmpq_poly);

    fmpq_roots_t rational_roots;
    acb_roots_t complex_roots;
    arb_roots_t real_roots;
    fmpq_roots_init(&rational_roots);
    acb_roots_init(&complex_roots);
    arb_roots_init(&real_roots);

    fmpq_poly_roots(&rational_roots, fmpq_poly, 0);
    fmpq_poly_acb_roots(&complex_roots, fmpq_poly, prec);
    acb_roots_to_real(&real_roots, &complex_roots, prec);

    arb_t *result_roots = NULL;
    slong result_count = 0;
    slong result_cap = 0;

    for (slong i = 0; i < rational_roots.num_roots; i++) {
        arb_t as_real;
        arb_init(as_real);
        arb_set_fmpq(as_real, rational_roots.roots[i], prec);
        rational_append_unique_real_root(&result_roots, &result_count, &result_cap, as_real);
        arb_clear(as_real);
    }

    for (slong i = 0; i < real_roots.num_roots; i++) {
        rational_append_unique_real_root(&result_roots, &result_count, &result_cap, real_roots.roots[i]);
    }

    *num_roots_out = result_count;

    if (result_count > 0) {
        printf("Found %ld real roots:\n", result_count);
        for (slong i = 0; i < result_count; i++) {
            printf("  Root %ld: ", i + 1);
            arb_printd(result_roots[i], RATIONAL_REAL_ROOT_DIGITS);
            printf("\n");
        }
    } else {
        printf("No real roots found.\n");
    }

    fmpq_roots_clear(&rational_roots);
    acb_roots_clear(&complex_roots);
    arb_roots_clear(&real_roots);
    fmpq_poly_clear(fmpq_poly);

    return result_roots;
}

char* rational_substitute_variable_in_polynomial(const char *poly_str, const char *var_name,
                                                  const fmpq_t value) {
    char *val_str = fmpq_get_str(NULL, 10, value);
    size_t val_len = val_str ? strlen(val_str) : 1;
    size_t estimated = strlen(poly_str) + (val_len + 2) * 8 + 32;
    char *result = (char*) malloc(estimated);
    result[0] = '\0';
    const char *ptr = poly_str;
    
    while (*ptr) {
        size_t var_len = strlen(var_name);
        if (strncmp(ptr, var_name, var_len) == 0 &&
            (ptr == poly_str || !isalnum(*(ptr-1))) &&
            !isalnum(*(ptr + var_len))) {
            const char *after_var = ptr + var_len;
            while (isspace((unsigned char) *after_var)) after_var++;

            if (*after_var == '^') {
                after_var++;
                while (isspace((unsigned char) *after_var)) after_var++;
                if (isdigit((unsigned char) *after_var)) {
                    char *endptr = NULL;
                    long exponent = strtol(after_var, &endptr, 10);
                    fmpq_t pow_val;
                    fmpq_init(pow_val);
                    fmpq_one(pow_val);
                    for (long i = 0; i < exponent; i++) {
                        fmpq_mul(pow_val, pow_val, value);
                    }
                    char *pow_str = fmpq_get_str(NULL, 10, pow_val);
                    strcat(result, pow_str ? pow_str : "0");
                    if (pow_str) flint_free(pow_str);
                    fmpq_clear(pow_val);
                    ptr = endptr;
                    continue;
                }
            }

            strcat(result, "(");
            strcat(result, val_str ? val_str : "0");
            strcat(result, ")");
            ptr += var_len;
        } else {
            size_t len = strlen(result);
            result[len] = *ptr;
            result[len + 1] = '\0';
            ptr++;
        }
    }

    if (val_str) {
        flint_free(val_str);
    }
    
    return result;
}

static int rational_verify_numeric_full_solution(char **poly_strings, slong num_polys,
                                                 rational_variable_info_t *vars, slong num_vars,
                                                 const double *values, double tol,
                                                 double *max_residual_out) {
    double max_residual = 0.0;
    for (slong eq = 0; eq < num_polys; eq++) {
        int ok = 0;
        double residual = rational_evaluate_polynomial_double(poly_strings[eq], vars, num_vars, values, &ok);
        if (!ok) {
            return 0;
        }
        residual = fabs(residual);
        if (residual > max_residual) {
            max_residual = residual;
        }
    }

    if (max_residual_out) {
        *max_residual_out = max_residual;
    }
    return max_residual <= tol;
}

static int rational_newton_solve_reduced_system(char **poly_strings,
                                                const slong *eq_indices,
                                                slong num_eqs,
                                                rational_variable_info_t *vars,
                                                slong num_vars,
                                                const double *fixed_values,
                                                slong fixed_count,
                                                const double *initial_guess,
                                                double *full_solution_out) {
    slong num_unknowns = num_vars - fixed_count;
    if (num_unknowns <= 0 || num_eqs != num_unknowns) {
        return 0;
    }

    double *x = (double *) malloc(num_unknowns * sizeof(double));
    double *F = (double *) malloc(num_unknowns * sizeof(double));
    double *J = (double *) malloc(num_unknowns * num_unknowns * sizeof(double));
    double *rhs = (double *) malloc(num_unknowns * sizeof(double));
    double *delta = (double *) calloc(num_unknowns, sizeof(double));
    double *full_values = (double *) calloc(num_vars, sizeof(double));
    if (!x || !F || !J || !rhs || !delta || !full_values) {
        free(x); free(F); free(J); free(rhs); free(delta); free(full_values);
        return 0;
    }

    for (slong i = 0; i < fixed_count; i++) {
        full_values[i] = fixed_values[i];
    }
    for (slong i = 0; i < num_unknowns; i++) {
        x[i] = initial_guess[i];
        full_values[fixed_count + i] = x[i];
    }

    int converged = 0;
    for (slong iter = 0; iter < RATIONAL_NUMERIC_NEWTON_MAX_ITER; iter++) {
        double max_F = 0.0;
        for (slong eq = 0; eq < num_unknowns; eq++) {
            int ok = 0;
            F[eq] = rational_evaluate_polynomial_double(poly_strings[eq_indices[eq]], vars, num_vars, full_values, &ok);
            if (!ok || !isfinite(F[eq])) {
                converged = 0;
                goto cleanup;
            }
            if (fabs(F[eq]) > max_F) {
                max_F = fabs(F[eq]);
            }
        }

        if (max_F < RATIONAL_NUMERIC_NEWTON_TOL) {
            converged = 1;
            break;
        }

        for (slong row = 0; row < num_unknowns; row++) {
            rhs[row] = -F[row];
            for (slong col = 0; col < num_unknowns; col++) {
                double saved = x[col];
                double h = 1e-6 * fmax(1.0, fabs(saved));

                x[col] = saved + h;
                full_values[fixed_count + col] = x[col];
                int ok_plus = 0;
                double f_plus = rational_evaluate_polynomial_double(poly_strings[eq_indices[row]], vars, num_vars, full_values, &ok_plus);

                x[col] = saved - h;
                full_values[fixed_count + col] = x[col];
                int ok_minus = 0;
                double f_minus = rational_evaluate_polynomial_double(poly_strings[eq_indices[row]], vars, num_vars, full_values, &ok_minus);

                x[col] = saved;
                full_values[fixed_count + col] = saved;

                if (!ok_plus || !ok_minus || !isfinite(f_plus) || !isfinite(f_minus)) {
                    converged = 0;
                    goto cleanup;
                }
                J[row * num_unknowns + col] = (f_plus - f_minus) / (2.0 * h);
            }
        }

        memset(delta, 0, num_unknowns * sizeof(double));
        if (!rational_solve_linear_system_double(J, rhs, delta, num_unknowns)) {
            converged = 0;
            goto cleanup;
        }

        double max_delta = 0.0;
        for (slong i = 0; i < num_unknowns; i++) {
            x[i] += delta[i];
            full_values[fixed_count + i] = x[i];
            if (!isfinite(x[i])) {
                converged = 0;
                goto cleanup;
            }
            if (fabs(delta[i]) > max_delta) {
                max_delta = fabs(delta[i]);
            }
        }

        if (max_delta < RATIONAL_NUMERIC_NEWTON_TOL) {
            converged = 1;
            break;
        }
    }

    if (converged) {
        double max_residual = 0.0;
        if (rational_verify_numeric_full_solution(poly_strings, num_eqs, vars, num_vars, full_values,
                                                  RATIONAL_NUMERIC_VERIFY_TOL, &max_residual)) {
            for (slong i = 0; i < num_vars; i++) {
                full_solution_out[i] = full_values[i];
            }
        } else {
            converged = 0;
        }
    }

cleanup:
    free(x);
    free(F);
    free(J);
    free(rhs);
    free(delta);
    free(full_values);
    return converged;
}

static int rational_full_solution_already_present(const rational_solutions_t *sols,
                                                  const double *values,
                                                  slong num_vars,
                                                  double tol) {
    if (!sols || !values) {
        return 0;
    }

    for (slong set = 0; set < sols->num_real_solution_sets; set++) {
        double *existing = (double *) malloc(num_vars * sizeof(double));
        if (!existing) {
            return 0;
        }
        int complete = 1;
        for (slong var = 0; var < num_vars; var++) {
            if (sols->real_solutions_per_var[set * num_vars + var] <= 0) {
                complete = 0;
                break;
            }
            existing[var] = rational_arb_to_double(sols->real_solution_sets[set][var][0]);
        }
        int match = complete && rational_solution_vectors_close(existing, values, num_vars, tol);
        free(existing);
        if (match) {
            return 1;
        }
    }

    return 0;
}

static void rational_try_numeric_backsolve_from_real_root(char **poly_strings, slong num_polys,
                                                          rational_variable_info_t *sorted_vars, slong num_vars,
                                                          const arb_t base_root,
                                                          rational_solutions_t *sols) {
    if (!poly_strings || !sorted_vars || !sols || num_vars <= 1) {
        return;
    }

    slong num_unknowns = num_vars - 1;
    double base_value = rational_arb_to_double(base_root);
    double *fixed_values = (double *) malloc(sizeof(double));
    double *candidate = (double *) calloc(num_vars, sizeof(double));
    slong *active_eq_indices = (slong *) malloc(num_polys * sizeof(slong));
    if (!fixed_values || !candidate || !active_eq_indices) {
        free(fixed_values); free(candidate); free(active_eq_indices);
        return;
    }

    fixed_values[0] = base_value;
    candidate[0] = base_value;

    slong num_active = 0;
    rational_variable_info_t *remaining_vars = sorted_vars + 1;
    for (slong eq = 0; eq < num_polys; eq++) {
        if (!rational_contains_any_variable(poly_strings[eq], remaining_vars, num_unknowns)) {
            int ok = 0;
            double constant_val = rational_evaluate_polynomial_double(poly_strings[eq], sorted_vars, num_vars, candidate, &ok);
            if (!ok || fabs(constant_val) > RATIONAL_NUMERIC_VERIFY_TOL) {
                free(fixed_values); free(candidate); free(active_eq_indices);
                return;
            }
            continue;
        }
        active_eq_indices[num_active++] = eq;
    }

    if (num_active < num_unknowns) {
        free(fixed_values); free(candidate); free(active_eq_indices);
        return;
    }

    rational_equation_combination_t *combinations = NULL;
    slong num_combinations = 0;
    if (num_active == num_unknowns) {
        num_combinations = 1;
        combinations = (rational_equation_combination_t *) malloc(sizeof(rational_equation_combination_t));
        combinations[0].equation_indices = (slong *) malloc(num_unknowns * sizeof(slong));
        combinations[0].num_equations = num_unknowns;
        combinations[0].tried = 0;
        combinations[0].success = 0;
        for (slong i = 0; i < num_unknowns; i++) {
            combinations[0].equation_indices[i] = i;
        }
    } else {
        rational_generate_equation_combinations(num_active, num_unknowns, &combinations, &num_combinations);
    }

    double scale = fmax(2.0, fabs(base_value) + 1.0);
    double grid5[] = {-scale, -0.5 * scale, 0.0, 0.5 * scale, scale};
    double grid3[] = {-scale, 0.0, scale};
    double *grid = (num_unknowns <= 2) ? grid5 : grid3;
    slong grid_count = (num_unknowns <= 2) ? 5 : 3;
    slong total_seed_count = 1;
    for (slong i = 0; i < num_unknowns; i++) {
        total_seed_count *= grid_count;
    }

    double *initial_guess = (double *) calloc(num_unknowns, sizeof(double));
    double *full_solution = (double *) calloc(num_vars, sizeof(double));
    slong *eq_indices = (slong *) calloc(num_unknowns, sizeof(slong));
    if (!initial_guess || !full_solution || !eq_indices) {
        free(initial_guess); free(full_solution); free(eq_indices);
        rational_free_equation_combinations(combinations, num_combinations);
        free(fixed_values); free(candidate); free(active_eq_indices);
        return;
    }

    for (slong comb_idx = 0; comb_idx < num_combinations; comb_idx++) {
        for (slong i = 0; i < num_unknowns; i++) {
            eq_indices[i] = active_eq_indices[combinations[comb_idx].equation_indices[i]];
        }

        for (slong seed_idx = 0; seed_idx < total_seed_count; seed_idx++) {
            slong code = seed_idx;
            for (slong j = 0; j < num_unknowns; j++) {
                initial_guess[j] = grid[code % grid_count];
                code /= grid_count;
            }

            if (rational_newton_solve_reduced_system(poly_strings, eq_indices, num_unknowns,
                                                     sorted_vars, num_vars,
                                                     fixed_values, 1,
                                                     initial_guess, full_solution)) {
                double max_residual = 0.0;
                if (rational_verify_numeric_full_solution(poly_strings, num_polys,
                                                          sorted_vars, num_vars,
                                                          full_solution, RATIONAL_NUMERIC_VERIFY_TOL,
                                                          &max_residual) &&
                    !rational_full_solution_already_present(sols, full_solution, num_vars, 1e-6)) {
                    rational_append_real_solution_values(sols, full_solution, num_vars);
                }
            }
        }
    }

    free(initial_guess);
    free(full_solution);
    free(eq_indices);
    rational_free_equation_combinations(combinations, num_combinations);
    free(fixed_values);
    free(candidate);
    free(active_eq_indices);
}

char* rational_eliminate_variable_dixon_with_selection(char **poly_strings, slong num_polys, 
                                                        const char **elim_vars, slong num_elim_vars,
                                                        rational_equation_combination_t *combination) {
    printf("  Entering rational_eliminate_variable_dixon_with_selection...\n");
    printf("    combination->num_equations = %ld, num_elim_vars = %ld\n", combination->num_equations, num_elim_vars);
    fflush(stdout);
    
    if (combination->num_equations == num_elim_vars + 1) {
        printf("    Using rational dixon_str_rational for elimination\n");
        fflush(stdout);
        
        size_t total_poly_len = 0;
        for (slong i = 0; i < combination->num_equations; i++) {
            total_poly_len += strlen(poly_strings[combination->equation_indices[i]]) + 2;
        }
        char *combined_polys = (char*) malloc(total_poly_len + 1);
        combined_polys[0] = '\0';
        for (slong i = 0; i < combination->num_equations; i++) {
            if (i > 0) {
                strcat(combined_polys, ", ");
            }
            strcat(combined_polys, poly_strings[combination->equation_indices[i]]);
        }
        
        size_t total_vars_len = 0;
        for (slong i = 0; i < num_elim_vars; i++) {
            total_vars_len += strlen(elim_vars[i]) + 2;
        }
        char *combined_vars = (char*) malloc(total_vars_len + 1);
        combined_vars[0] = '\0';
        for (slong i = 0; i < num_elim_vars; i++) {
            if (i > 0) {
                strcat(combined_vars, ", ");
            }
            strcat(combined_vars, elim_vars[i]);
        }
        
        printf("    Calling dixon_str_rational_with_file...\n");
        printf("    Polys: %s\n", combined_polys);
        printf("    Elim:  %s\n", combined_vars);
        fflush(stdout);
        
        char *result = dixon_str_rational_with_file(combined_polys, combined_vars, NULL);
        
        printf("    dixon_str_rational_with_file returned: %s\n", result ? result : "NULL");
        fflush(stdout);
        
        free(combined_polys);
        free(combined_vars);
        return result;
    } else {
        printf("    Number of equations does not match, returning 0\n");
        fflush(stdout);
        return strdup("0");
    }
}

int rational_verify_solution_set(char **original_polys, slong num_polys,
                                 rational_variable_info_t *sorted_vars, slong num_vars,
                                 fmpq_t **solution_values) {
    if (num_vars <= 0) {
        return 0;
    }

    printf("    Verifying rational solution: ");
    for (slong v = 0; v < num_vars; v++) {
        if (v > 0) printf(", ");
        printf("%s=", sorted_vars[v].name);
        fmpq_print(solution_values[v][0]);
    }
    printf("\n");

    double *values = (double *) calloc(num_vars, sizeof(double));
    if (!values) {
        return 0;
    }

    for (slong var_idx = 0; var_idx < num_vars; var_idx++) {
        values[var_idx] = rational_fmpq_to_double(solution_values[var_idx][0]);
    }

    double max_residual = 0.0;
    int ok = rational_verify_numeric_full_solution(original_polys, num_polys,
                                                   sorted_vars, num_vars,
                                                   values, 1e-10, &max_residual);
    free(values);

    if (!ok) {
        printf("    Verification FAILED (max residual %.3e)\n", max_residual);
        return 0;
    }

    printf("    Verification PASSED\n");
    return 1;
}

void rational_filter_solutions_by_verification(rational_solutions_t *sols,
                                                char **original_polys, slong num_polys,
                                                rational_variable_info_t *sorted_vars) {
    if (!sols->is_valid || sols->has_no_solutions || sols->num_solution_sets == 0) {
        return;
    }
    
    printf("\n=== Verifying Rational Solution Sets ===\n");
    
    slong verified_count = 0;
    int *valid_sets = (int*) calloc(sols->num_solution_sets, sizeof(int));
    
    for (slong set = 0; set < sols->num_solution_sets; set++) {
        printf("  Checking rational solution set %ld...\n", set + 1);
        int pass = 0;
        char *candidate_line = rational_format_solution_set_line_from_arrays(
            sols->variable_names,
            sols->num_variables,
            sols->solution_sets[set],
            &sols->solutions_per_var[set * sols->num_variables]
        );
        
        int is_complete = 1;
        for (slong var = 0; var < sols->num_variables; var++) {
            if (sols->solutions_per_var[set * sols->num_variables + var] == 0) {
                printf("    Solution set %ld incomplete (missing value for %s)\n", 
                       set + 1, sols->variable_names[var]);
                is_complete = 0;
                break;
            }
        }
        
        if (!is_complete) {
            rational_add_candidate_solution(sols, candidate_line, 0);
            continue;
        }
        
        if (rational_verify_solution_set(original_polys, num_polys, sorted_vars, sols->num_variables,
                                        sols->solution_sets[set])) {
            valid_sets[set] = 1;
            verified_count++;
            pass = 1;
        }
        rational_add_candidate_solution(sols, candidate_line, pass);
    }
    
    printf("  Verification complete: %ld out of %ld solution sets are valid\n", 
           verified_count, sols->num_solution_sets);
    sols->checked_solution_sets = sols->num_solution_sets;
    sols->verified_solution_sets = verified_count;
    
    if (verified_count == 0) {
        sols->has_no_solutions = 1;
        
        if (sols->solution_sets) {
            for (slong set = 0; set < sols->num_solution_sets; set++) {
                if (sols->solution_sets[set]) {
                    for (slong var = 0; var < sols->num_variables; var++) {
                        if (sols->solution_sets[set][var]) {
                            for (slong sol = 0; sol < sols->solutions_per_var[set * sols->num_variables + var]; sol++) {
                                fmpq_clear(sols->solution_sets[set][var][sol]);
                            }
                            free(sols->solution_sets[set][var]);
                        }
                    }
                    free(sols->solution_sets[set]);
                }
            }
            free(sols->solution_sets);
            sols->solution_sets = NULL;
        }
        
        if (sols->solutions_per_var) {
            free(sols->solutions_per_var);
            sols->solutions_per_var = NULL;
        }
        
        sols->num_solution_sets = 0;
        
    } else if (verified_count < sols->num_solution_sets) {
        fmpq_t ***new_solution_sets = (fmpq_t***) malloc(verified_count * sizeof(fmpq_t**));
        slong *new_solutions_per_var = (slong*) calloc(verified_count * sols->num_variables, sizeof(slong));
        
        slong new_idx = 0;
        for (slong set = 0; set < sols->num_solution_sets; set++) {
            if (valid_sets[set]) {
                new_solution_sets[new_idx] = sols->solution_sets[set];
                for (slong var = 0; var < sols->num_variables; var++) {
                    new_solutions_per_var[new_idx * sols->num_variables + var] = 
                        sols->solutions_per_var[set * sols->num_variables + var];
                }
                sols->solution_sets[set] = NULL;
                new_idx++;
            } else {
                if (sols->solution_sets[set]) {
                    for (slong var = 0; var < sols->num_variables; var++) {
                        if (sols->solution_sets[set][var]) {
                            for (slong sol = 0; sol < sols->solutions_per_var[set * sols->num_variables + var]; sol++) {
                                fmpq_clear(sols->solution_sets[set][var][sol]);
                            }
                            free(sols->solution_sets[set][var]);
                        }
                    }
                    free(sols->solution_sets[set]);
                }
            }
        }
        
        free(sols->solution_sets);
        free(sols->solutions_per_var);
        sols->solution_sets = new_solution_sets;
        sols->solutions_per_var = new_solutions_per_var;
        sols->num_solution_sets = verified_count;
    }
    
    free(valid_sets);
    printf("=== Verification Complete ===\n\n");
}

int rational_solve_by_back_substitution_recursive_enhanced(char **original_polys, slong num_polys,
                                                             rational_variable_info_t *sorted_vars, slong num_vars,
                                                             fmpq_t *base_solutions, slong num_base_solutions,
                                                             rational_solutions_t *sols) {
    if (sols->has_no_solutions == -1) {
        printf("Back substitution SKIPPED: system already has dimension > 0 status\n");
        return 1;
    }
    
    if (num_vars <= 1) {
        return 1; 
    }
    
    printf("Starting rational back substitution process...\n");
    printf("Base variable %s has %ld solutions\n", sorted_vars[0].name, num_base_solutions);
    
    fmpq_t ***all_solution_sets = NULL;
    slong *all_solutions_per_var = NULL;
    slong total_solution_sets = 0;
    
    for (slong base_idx = 0; base_idx < num_base_solutions; base_idx++) {
        printf("\n--- Processing rational base solution %ld: ", base_idx + 1);
        fmpq_print(base_solutions[base_idx]);
        printf(" ---\n");
        
        char **reduced_polys = (char**) malloc(num_polys * sizeof(char*));
        slong num_nonzero_polys = 0;
        
        for (slong poly_idx = 0; poly_idx < num_polys; poly_idx++) {
            char *substituted = rational_substitute_variable_in_polynomial(original_polys[poly_idx],
                                                                          sorted_vars[0].name,
                                                                          base_solutions[base_idx]);
            printf("Equation %ld after substituting %s: %s\n", 
                   poly_idx, sorted_vars[0].name, substituted);
            
            if (strcmp(substituted, "0") != 0) {
                reduced_polys[num_nonzero_polys] = substituted;
                num_nonzero_polys++;
            } else {
                printf("Equation %ld became zero (satisfied)\n", poly_idx);
                free(substituted);
            }
        }
        
        if (num_nonzero_polys == 0) {
            printf("All equations satisfied by base solution\n");
            if (total_solution_sets == 0) {
                all_solution_sets = (fmpq_t***) malloc(sizeof(fmpq_t**));
                all_solutions_per_var = (slong*) calloc(num_vars, sizeof(slong));
            } else {
                all_solution_sets = (fmpq_t***) realloc(all_solution_sets, 
                                                          (total_solution_sets + 1) * sizeof(fmpq_t**));
                all_solutions_per_var = (slong*) realloc(all_solutions_per_var,
                                                         (total_solution_sets + 1) * num_vars * sizeof(slong));
            }
            
            all_solution_sets[total_solution_sets] = (fmpq_t**) calloc(num_vars, sizeof(fmpq_t*));
            all_solution_sets[total_solution_sets][0] = (fmpq_t*) malloc(sizeof(fmpq_t));
            fmpq_init(all_solution_sets[total_solution_sets][0][0]);
            fmpq_set(all_solution_sets[total_solution_sets][0][0], base_solutions[base_idx]);
            all_solutions_per_var[total_solution_sets * num_vars + 0] = 1;
            
            for (slong v = 1; v < num_vars; v++) {
                all_solutions_per_var[total_solution_sets * num_vars + v] = 0;
            }
            
            total_solution_sets++;
            free(reduced_polys);
            continue;
        }
        
        slong target_equations = num_vars - 1;
        
        if (num_nonzero_polys >= target_equations) {
            rational_equation_combination_t *combinations = NULL;
            slong num_combinations = 0;
            
            if (num_nonzero_polys == target_equations) {
                num_combinations = 1;
                combinations = (rational_equation_combination_t*) malloc(sizeof(rational_equation_combination_t));
                combinations[0].equation_indices = (slong*) malloc(target_equations * sizeof(slong));
                combinations[0].num_equations = target_equations;
                for (slong i = 0; i < target_equations; i++) {
                    combinations[0].equation_indices[i] = i;
                }
            } else {
                rational_generate_equation_combinations(num_nonzero_polys, target_equations, 
                                                         &combinations, &num_combinations);
            }
            
            printf("Generated %ld equation combinations to try\n", num_combinations);
            
            int found_finite_solution = 0;
            
            for (slong comb_idx = 0; comb_idx < num_combinations; comb_idx++) {
                printf("  Trying combination %ld: equations ", comb_idx + 1);
                for (slong i = 0; i < target_equations; i++) {
                    if (i > 0) printf(", ");
                    printf("%ld", combinations[comb_idx].equation_indices[i]);
                }
                printf("\n");
                
                char **final_polys = (char**) malloc(target_equations * sizeof(char*));
                for (slong i = 0; i < target_equations; i++) {
                    final_polys[i] = strdup(reduced_polys[combinations[comb_idx].equation_indices[i]]);
                }

                {
                    char *base_str = fmpq_get_str(NULL, 10, base_solutions[base_idx]);
                    char *eq_list = rational_format_equation_index_list(combinations[comb_idx].equation_indices,
                                                               target_equations);
                    if (target_equations == 1) {
                        rational_add_resultant_step(sols,
                                           "Base %s = %s -> solve reduced eq(%s)",
                                           sorted_vars[0].name,
                                           base_str ? base_str : "?",
                                           eq_list ? eq_list : "?");
                    } else {
                        rational_add_resultant_step(sols,
                                           "Base %s = %s -> compute reduced resultants from eq(%s)",
                                           sorted_vars[0].name,
                                           base_str ? base_str : "?",
                                           eq_list ? eq_list : "?");
                    }
                    if (eq_list) free(eq_list);
                    if (base_str) flint_free(base_str);
                }
                
                rational_variable_info_t *remaining_vars = (rational_variable_info_t*) malloc((num_vars - 1) * sizeof(rational_variable_info_t));
                for (slong i = 1; i < num_vars; i++) {
                    remaining_vars[i - 1] = sorted_vars[i];
                }
                
                rational_solutions_t *test_sols = solve_rational_polynomial_system_array_with_vars(final_polys, 
                                                                                           target_equations,
                                                                                           remaining_vars,
                                                                                           num_vars - 1);
                
                if (test_sols && test_sols->is_valid) {
                    rational_merge_solver_logs(sols, test_sols);
                    if (test_sols->has_no_solutions == -1) {
                        printf("  Combination %ld: dimension > 0\n", comb_idx + 1);
                    } else if (test_sols->has_no_solutions == 1) {
                        printf("  Combination %ld: no solutions\n", comb_idx + 1);
                    } else {
                        printf("  ✓ Combination %ld succeeded with finite solutions!\n", comb_idx + 1);
                        found_finite_solution++;
                        
                        if (test_sols->num_solution_sets > 0) {
                            slong old_total = total_solution_sets;
                            total_solution_sets += test_sols->num_solution_sets;
                            
                            if (old_total == 0) {
                                all_solution_sets = (fmpq_t***) malloc(total_solution_sets * sizeof(fmpq_t**));
                                all_solutions_per_var = (slong*) calloc(total_solution_sets * num_vars, sizeof(slong));
                            } else {
                                all_solution_sets = (fmpq_t***) realloc(all_solution_sets, 
                                                                          total_solution_sets * sizeof(fmpq_t**));
                                all_solutions_per_var = (slong*) realloc(all_solutions_per_var,
                                                                         total_solution_sets * num_vars * sizeof(slong));
                            }
                            
                            for (slong ts = 0; ts < test_sols->num_solution_sets; ts++) {
                                slong new_set_idx = old_total + ts;
                                all_solution_sets[new_set_idx] = (fmpq_t**) calloc(num_vars, sizeof(fmpq_t*));
                                
                                all_solution_sets[new_set_idx][0] = (fmpq_t*) malloc(sizeof(fmpq_t));
                                fmpq_init(all_solution_sets[new_set_idx][0][0]);
                                fmpq_set(all_solution_sets[new_set_idx][0][0], base_solutions[base_idx]);
                                all_solutions_per_var[new_set_idx * num_vars + 0] = 1;
                                
                                for (slong v = 1; v < num_vars; v++) {
                                    slong src_var_idx = v - 1;
                                    slong num_sols = test_sols->solutions_per_var[ts * (num_vars - 1) + src_var_idx];
                                    all_solutions_per_var[new_set_idx * num_vars + v] = num_sols;
                                    if (num_sols > 0) {
                                        all_solution_sets[new_set_idx][v] = (fmpq_t*) malloc(num_sols * sizeof(fmpq_t));
                                        for (slong s = 0; s < num_sols; s++) {
                                            fmpq_init(all_solution_sets[new_set_idx][v][s]);
                                            fmpq_set(all_solution_sets[new_set_idx][v][s], 
                                                     test_sols->solution_sets[ts][src_var_idx][s]);
                                        }
                                    } else {
                                        all_solution_sets[new_set_idx][v] = NULL;
                                    }
                                }
                            }
                        }

                        if (test_sols->num_real_solution_sets > 0) {
                            double *combined_values = (double *) calloc(num_vars, sizeof(double));
                            if (combined_values) {
                                combined_values[0] = rational_fmpq_to_double(base_solutions[base_idx]);
                                for (slong ts = 0; ts < test_sols->num_real_solution_sets; ts++) {
                                    int complete_real = 1;
                                    for (slong v = 1; v < num_vars; v++) {
                                        slong src_var_idx = v - 1;
                                        if (test_sols->real_solutions_per_var[ts * (num_vars - 1) + src_var_idx] <= 0) {
                                            complete_real = 0;
                                            break;
                                        }
                                        combined_values[v] = rational_arb_to_double(
                                            test_sols->real_solution_sets[ts][src_var_idx][0]
                                        );
                                    }
                                    if (complete_real &&
                                        !rational_full_solution_already_present(sols, combined_values, num_vars, 1e-6)) {
                                        rational_append_real_solution_values(sols, combined_values, num_vars);
                                    }
                                }
                                free(combined_values);
                            }
                        }
                    }
                }
                
                rational_solutions_clear(test_sols);
                free(test_sols);
                
                for (slong i = 0; i < target_equations; i++) {
                    free(final_polys[i]);
                }
                free(final_polys);
                free(remaining_vars);
                
                if (found_finite_solution > 0) {
                    break;
                }
            }
            
            rational_free_equation_combinations(combinations, num_combinations);
        }
        
        for (slong i = 0; i < num_nonzero_polys; i++) {
            free(reduced_polys[i]);
        }
        free(reduced_polys);
    }
    
    if (total_solution_sets > 0) {
        sols->solution_sets = all_solution_sets;
        sols->solutions_per_var = all_solutions_per_var;
        sols->num_solution_sets = total_solution_sets;
    }
    
    return 1;
}

int rational_solve_by_elimination_enhanced(char **poly_strings, slong num_polys,
                                             rational_variable_info_t *sorted_vars, slong num_vars,
                                             rational_solutions_t *sols) {
    printf("\n=== Rational Elimination Solver ===\n");
    printf("Number of variables: %ld\n", num_vars);
    printf("Number of equations: %ld\n", num_polys);
    
    printf("Sorted variables order:\n");
    for (slong i = 0; i < num_vars; i++) {
        printf("  [%ld] %s (max_degree=%ld)\n", i, sorted_vars[i].name, sorted_vars[i].max_degree);
    }
    fflush(stdout);
    
    if (num_vars == 0) {
        printf("No variables to solve for.\n");
        return 1;
    }
    
    if (num_vars == 1) {
        printf("Single variable case: solving directly\n");
        slong num_roots = 0;
        fmpq_t *roots = rational_solve_univariate_equation_all_roots(poly_strings[0], sorted_vars[0].name, &num_roots);
        slong num_real_roots = 0;
        arb_t *real_roots = rational_solve_univariate_equation_all_real_roots(
            poly_strings[0], sorted_vars[0].name, &num_real_roots, sols->real_solution_precision
        );
        sols->num_base_solutions = num_real_roots;
        
        if (num_roots > 0) {
            sols->solution_sets = (fmpq_t***) calloc(num_roots, sizeof(fmpq_t**));
            sols->solutions_per_var = (slong*) calloc(num_roots, sizeof(slong));
            for (slong i = 0; i < num_roots; i++) {
                sols->solution_sets[i] = (fmpq_t**) calloc(1, sizeof(fmpq_t*));
                sols->solution_sets[i][0] = (fmpq_t*) malloc(sizeof(fmpq_t));
                fmpq_init(sols->solution_sets[i][0][0]);
                fmpq_set(sols->solution_sets[i][0][0], roots[i]);
                sols->solutions_per_var[i] = 1;
            }
            sols->num_solution_sets = num_roots;
        }

        if (num_real_roots > 0) {
            for (slong i = 0; i < num_real_roots; i++) {
                if (!rational_real_root_matches_rational(real_roots[i], roots, num_roots, sols->real_solution_precision)) {
                    double value = rational_arb_to_double(real_roots[i]);
                    rational_append_real_solution_values(sols, &value, 1);
                }
            }
        }

        if (roots) {
            for (slong i = 0; i < num_roots; i++) {
                fmpq_clear(roots[i]);
            }
            free(roots);
        }
        rational_clear_real_root_array(real_roots, num_real_roots);

        if (num_roots == 0 && num_real_roots == 0) {
            sols->has_no_solutions = 1;
        }
        
        return 1;
    }
    
    slong num_elim_vars = num_vars - 1;
    const char **elim_vars = (const char**) malloc(num_elim_vars * sizeof(char*));
    for (slong i = 1; i < num_vars; i++) {
        elim_vars[i - 1] = sorted_vars[i].name;
    }
    char *keep_var = sorted_vars[0].name;
    
    printf("Eliminating variables: ");
    for (slong i = 0; i < num_elim_vars; i++) {
        if (i > 0) printf(", ");
        printf("%s", elim_vars[i]);
    }
    printf("\nKeeping variable: %s\n", keep_var);
    fflush(stdout);
    
    slong target_equations = num_vars;
    
    if (num_polys < target_equations) {
        fprintf(stderr, "Error: Need at least %ld equations for %ld variables, got %ld\n",
                target_equations, num_vars, num_polys);
        free(elim_vars);
        sols->has_no_solutions = -1;
        return 0;
    }
    
    rational_equation_combination_t *combinations = NULL;
    slong num_combinations = 0;
    
    if (num_polys == target_equations) {
        num_combinations = 1;
        combinations = (rational_equation_combination_t*) malloc(sizeof(rational_equation_combination_t));
        combinations[0].equation_indices = (slong*) malloc(target_equations * sizeof(slong));
        combinations[0].num_equations = target_equations;
        for (slong i = 0; i < target_equations; i++) {
            combinations[0].equation_indices[i] = i;
        }
    } else {
        rational_generate_equation_combinations(num_polys, target_equations, 
                                                 &combinations, &num_combinations);
    }
    
    printf("Generated %ld equation combinations to try\n", num_combinations);
    sols->total_combinations = num_combinations;
    
    int success = 0;
    char *working_resultant = NULL;
    slong successful_combinations = 0;
    
    for (slong comb_idx = 0; comb_idx < num_combinations; comb_idx++) {
        printf("\n=== Trying combination %ld ===\n", comb_idx + 1);
        printf("Equations: ");
        for (slong i = 0; i < combinations[comb_idx].num_equations; i++) {
            if (i > 0) printf(", ");
            printf("%ld", combinations[comb_idx].equation_indices[i] + 1);
        }
        printf("\n");
        
        combinations[comb_idx].tried = 1;
        
        char *resultant = rational_eliminate_variable_dixon_with_selection(
            poly_strings, num_polys,
            elim_vars, num_elim_vars,
            &combinations[comb_idx]
        );
        
        printf("Resultant polynomial: %s\n", resultant);
        
        if (resultant && strcmp(resultant, "0") != 0) {
            combinations[comb_idx].success = 1;
            successful_combinations++;
            
            if (!working_resultant) {
                working_resultant = resultant;
            } else {
                free(resultant);
            }
        } else {
            if (resultant) free(resultant);
        }
    }
    
    free(elim_vars);
    
    if (successful_combinations == 0) {
        printf("All %ld equation combinations resulted in zero resultant\n", num_combinations);
        sols->successful_combinations = 0;
        sols->has_no_solutions = -1;
        rational_free_equation_combinations(combinations, num_combinations);
        return 1;
    }
    
    sols->successful_combinations = successful_combinations;
    printf("Found %ld non-zero resultant(s) out of %ld combinations, proceeding with solving...\n", 
           successful_combinations, num_combinations);
    
    if (working_resultant) {
        slong num_roots = 0;
        fmpq_t *base_roots = rational_solve_univariate_equation_all_roots(working_resultant, keep_var, &num_roots);
        slong num_real_roots = 0;
        arb_t *real_base_roots = rational_solve_univariate_equation_all_real_roots(
            working_resultant, keep_var, &num_real_roots, sols->real_solution_precision
        );
        sols->num_base_solutions = (num_real_roots > 0) ? num_real_roots : num_roots;
        
        if (num_roots > 0) {
            rational_solve_by_back_substitution_recursive_enhanced(
                poly_strings, num_polys, sorted_vars, num_vars,
                base_roots, num_roots, sols
            );
            
            for (slong i = 0; i < num_roots; i++) {
                fmpq_clear(base_roots[i]);
            }
            free(base_roots);
            success = 1;
        } else {
            printf("No rational roots found for the resultant.\n");
        }

        if (num_real_roots > num_roots) {
            rational_set_real_root_summary(sols, keep_var, real_base_roots, num_real_roots, num_roots);
            for (slong i = 0; i < num_real_roots; i++) {
                if (!rational_real_root_matches_rational(real_base_roots[i], base_roots, num_roots,
                                                         sols->real_solution_precision)) {
                    rational_try_numeric_backsolve_from_real_root(poly_strings, num_polys,
                                                                  sorted_vars, num_vars,
                                                                  real_base_roots[i], sols);
                }
            }
        }

        rational_clear_real_root_array(real_base_roots, num_real_roots);
        
        free(working_resultant);
    }
    
    rational_free_equation_combinations(combinations, num_combinations);
    
    if (!success && sols->num_real_solution_sets == 0) {
        sols->has_no_solutions = 1;
    } else if (sols->num_real_solution_sets > 0) {
        sols->has_no_solutions = 0;
    }
    
    return success || (sols->num_real_solution_sets > 0);
}

rational_solutions_t* solve_rational_polynomial_system_array_with_vars(char **poly_strings, slong num_polys,
                                                                       rational_variable_info_t *original_vars, slong num_original_vars) {
    rational_solutions_t *sols = (rational_solutions_t*) malloc(sizeof(rational_solutions_t));
    rational_solutions_init(sols, num_original_vars);
    
    for (slong i = 0; i < num_original_vars; i++) {
        sols->variable_names[i] = strdup(original_vars[i].name);
    }
    sols->num_equations = num_polys;
    
    rational_solve_by_elimination_enhanced(poly_strings, num_polys, original_vars, num_original_vars, sols);
    sols->is_valid = 1;
    
    if (sols->num_solution_sets > 0) {
        rational_filter_solutions_by_verification(sols, poly_strings, num_polys, original_vars);
    }
    if (sols->num_real_solution_sets > 0) {
        rational_filter_real_solutions_by_verification(sols, poly_strings, original_vars);
    }
    if (sols->num_solution_sets > 0) {
        rational_compute_exact_solution_residuals(sols, poly_strings, original_vars);
    }
    if (sols->num_real_solution_sets > 0) {
        rational_compute_real_solution_residuals(sols, poly_strings, original_vars);
    }
    return sols;
}

rational_solutions_t* solve_rational_polynomial_system_array(char **poly_strings, slong num_polys) {
    rational_variable_info_t *sorted_vars = NULL;
    slong num_vars = 0;
    
    if (!rational_extract_and_sort_variables(poly_strings, num_polys, &sorted_vars, &num_vars)) {
        return NULL;
    }
    
    rational_solutions_t *result = solve_rational_polynomial_system_array_with_vars(
        poly_strings, num_polys, sorted_vars, num_vars
    );
    
    for (slong i = 0; i < num_vars; i++) {
        free(sorted_vars[i].name);
    }
    free(sorted_vars);
    
    return result;
}

rational_solutions_t* solve_rational_polynomial_system_string(const char *poly_string) {
    slong num_polys = 0;
    char **poly_strings = split_string((char*)poly_string, &num_polys);
    
    if (!poly_strings || num_polys == 0) {
        return NULL;
    }
    
    rational_solutions_t *result = solve_rational_polynomial_system_array(poly_strings, num_polys);
    
    free_split_strings(poly_strings, num_polys);
    
    return result;
}

void test_rational_polynomial_solver(void) {
    printf("\n=== Testing Rational Polynomial System Solver ===\n");
    
    char *test1[] = {
        "x + y - 3",
        "x - y + 1"
    };
    
    printf("\nTest 1: Linear system\n");
    rational_solutions_t *sols1 = solve_rational_polynomial_system_array(test1, 2);
    if (sols1) {
        print_rational_solutions(sols1);
        rational_solutions_clear(sols1);
        free(sols1);
    }
    
    char *test2[] = {
        "x^2 - 2",
        "x - 1"
    };
    
    printf("\nTest 2: Quadratic with no rational solution\n");
    rational_solutions_t *sols2 = solve_rational_polynomial_system_array(test2, 2);
    if (sols2) {
        print_rational_solutions(sols2);
        rational_solutions_clear(sols2);
        free(sols2);
    }

    char *test3[] = {
        "x^2 - 2"
    };

    printf("\nTest 3: Univariate irrational real roots\n");
    rational_solutions_t *sols3 = solve_rational_polynomial_system_array(test3, 1);
    if (sols3) {
        print_rational_solutions(sols3);
        rational_solutions_clear(sols3);
        free(sols3);
    }
    
    printf("\n=== Test Complete ===\n");
}
