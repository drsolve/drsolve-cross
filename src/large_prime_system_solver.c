#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "large_prime_system_solver.h"

#include "dixon_interface_flint.h"
#include "rational_system_solver.h"

#include <ctype.h>
#include <stdarg.h>

#include <flint/flint.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_mpoly.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_poly_factor.h>

#define LP_SOLVE_ERROR   (-1)
#define LP_SOLVE_NONE    0
#define LP_SOLVE_FOUND   1
#define LP_SOLVE_DIMPOS  2

typedef struct {
    fmpz_t *values;
    slong *mults;
    slong count;
    slong alloc;
} lp_root_list_t;

typedef struct {
    large_prime_solutions_t *sols;
    char **original_polys;
    slong num_original_polys;
    rational_variable_info_t *original_vars;
    slong num_original_vars;
} lp_solver_state_t;

typedef struct {
    const char *cursor;
    rational_variable_info_t *vars;
    slong num_vars;
    fmpz_t *values;
    const fmpz_mod_ctx_struct *ctx;
    int ok;
} lp_eval_parser_t;

typedef struct {
    const char *cursor;
    const char *var_name;
    const fmpz_mod_ctx_struct *ctx;
    int ok;
} lp_univar_parser_t;

static int large_prime_solver_realtime_progress_enabled = 0;

static void lp_progress(const char *fmt, ...)
{
    if (!large_prime_solver_realtime_progress_enabled)
        return;

    va_list args;
    va_start(args, fmt);
    fprintf(stderr, "Progress: ");
    vfprintf(stderr, fmt, args);
    fprintf(stderr, "\n");
    fflush(stderr);
    va_end(args);
}

void large_prime_solver_set_realtime_progress(int enabled)
{
    large_prime_solver_realtime_progress_enabled = enabled;
}

static void lp_set_error(large_prime_solutions_t *sols, const char *message)
{
    if (!sols)
        return;

    if (sols->error_message) {
        free(sols->error_message);
        sols->error_message = NULL;
    }

    if (message)
        sols->error_message = strdup(message);

    sols->is_valid = 0;
}

static void lp_solutions_init(large_prime_solutions_t *sols, slong num_vars, const fmpz_t prime)
{
    sols->num_variables = num_vars;
    sols->num_solution_sets = 0;
    sols->num_equations = 0;
    sols->variable_names = (char **) calloc((size_t) (num_vars > 0 ? num_vars : 1), sizeof(char *));
    sols->solution_values = NULL;
    sols->is_valid = 0;
    sols->has_no_solutions = 0;
    sols->error_message = NULL;
    fmpz_init(sols->prime);
    fmpz_set(sols->prime, prime);
    fmpz_mod_ctx_init(sols->ctx, prime);
}

static void lp_free_string_array(char **strings, slong count)
{
    if (!strings)
        return;

    for (slong i = 0; i < count; i++)
        free(strings[i]);
    free(strings);
}

static void lp_free_rational_vars(rational_variable_info_t *vars, slong num_vars)
{
    if (!vars)
        return;

    for (slong i = 0; i < num_vars; i++)
        free(vars[i].name);
    free(vars);
}

static int lp_is_name_in_list(const char *name, char **list, slong count)
{
    for (slong i = 0; i < count; i++) {
        if (strcmp(name, list[i]) == 0)
            return 1;
    }
    return 0;
}

static void lp_skip_ws(const char **cursor)
{
    while (**cursor && isspace((unsigned char) **cursor))
        (*cursor)++;
}

static int lp_is_identifier_char(char c)
{
    return isalnum((unsigned char) c) || c == '_';
}

static void lp_set_small_sample_value(fmpz_t out, ulong value, const fmpz_mod_ctx_t ctx)
{
    fmpz_t tmp;
    fmpz_init(tmp);
    fmpz_set_ui(tmp, value);
    fmpz_mod_set_fmpz(out, tmp, ctx);
    fmpz_clear(tmp);
}

static int lp_parse_integer_token(const char **cursor, fmpz_t value)
{
    const char *start = *cursor;
    char *buf;
    size_t len;

    if (!isdigit((unsigned char) **cursor))
        return 0;

    while (isdigit((unsigned char) **cursor))
        (*cursor)++;

    len = (size_t) (*cursor - start);
    buf = (char *) malloc(len + 1);
    if (!buf)
        return 0;

    memcpy(buf, start, len);
    buf[len] = '\0';
    if (fmpz_set_str(value, buf, 10) != 0) {
        free(buf);
        return 0;
    }
    free(buf);
    return 1;
}

static int lp_lookup_eval_variable(lp_eval_parser_t *parser, const char *name, size_t len, fmpz_t result)
{
    for (slong i = 0; i < parser->num_vars; i++) {
        if (strlen(parser->vars[i].name) == len &&
            strncmp(parser->vars[i].name, name, len) == 0) {
            fmpz_set(result, parser->values[parser->vars[i].index]);
            return 1;
        }
    }
    parser->ok = 0;
    return 0;
}

static void lp_eval_expression(lp_eval_parser_t *parser, fmpz_t result);

static void lp_fmpz_mod_pow_si(fmpz_t result, const fmpz_t base, long exponent, const fmpz_mod_ctx_t ctx)
{
    fmpz_t acc, pow_base, inv;
    ulong e;

    fmpz_init(acc);
    fmpz_init(pow_base);
    fmpz_set_ui(acc, 1);
    fmpz_set(pow_base, base);

    if (exponent < 0) {
        fmpz_init(inv);
        if (!fmpz_is_zero(base) && fmpz_mod_is_invertible(base, ctx)) {
            fmpz_mod_inv(inv, base, ctx);
            fmpz_set(pow_base, inv);
        } else {
            fmpz_zero(acc);
            exponent = 0;
        }
        fmpz_clear(inv);
        e = (ulong) (-exponent);
    } else {
        e = (ulong) exponent;
    }

    while (e > 0) {
        if (e & 1)
            fmpz_mod_mul(acc, acc, pow_base, ctx);
        e >>= 1;
        if (e)
            fmpz_mod_mul(pow_base, pow_base, pow_base, ctx);
    }

    fmpz_set(result, acc);
    fmpz_clear(pow_base);
    fmpz_clear(acc);
}

static void lp_eval_primary(lp_eval_parser_t *parser, fmpz_t result)
{
    lp_skip_ws(&parser->cursor);

    if (!parser->ok) {
        fmpz_zero(result);
        return;
    }

    if (*parser->cursor == '\0') {
        parser->ok = 0;
        fmpz_zero(result);
        return;
    }

    if (*parser->cursor == '(') {
        parser->cursor++;
        lp_eval_expression(parser, result);
        lp_skip_ws(&parser->cursor);
        if (*parser->cursor == ')') {
            parser->cursor++;
        } else {
            parser->ok = 0;
        }
        return;
    }

    if (*parser->cursor == '+' || *parser->cursor == '-') {
        char op = *parser->cursor++;
        lp_eval_primary(parser, result);
        if (op == '-')
            fmpz_mod_neg(result, result, parser->ctx);
        return;
    }

    if (isalpha((unsigned char) *parser->cursor) || *parser->cursor == '_') {
        const char *start = parser->cursor;
        while (lp_is_identifier_char(*parser->cursor))
            parser->cursor++;
        if (!lp_lookup_eval_variable(parser, start, (size_t) (parser->cursor - start), result))
            fmpz_zero(result);
        return;
    }

    if (isdigit((unsigned char) *parser->cursor)) {
        fmpz_t tmp;
        fmpz_init(tmp);
        if (!lp_parse_integer_token(&parser->cursor, tmp)) {
            parser->ok = 0;
            fmpz_zero(result);
        } else {
            fmpz_mod_set_fmpz(result, tmp, parser->ctx);
        }
        fmpz_clear(tmp);
        return;
    }

    parser->ok = 0;
    fmpz_zero(result);
}

static void lp_eval_power(lp_eval_parser_t *parser, fmpz_t result)
{
    fmpz_t base;

    fmpz_init(base);
    lp_eval_primary(parser, base);
    lp_skip_ws(&parser->cursor);

    while (parser->ok && *parser->cursor == '^') {
        char *endptr = NULL;
        long exponent;

        parser->cursor++;
        lp_skip_ws(&parser->cursor);
        exponent = strtol(parser->cursor, &endptr, 10);
        if (endptr == parser->cursor) {
            parser->ok = 0;
            break;
        }
        parser->cursor = endptr;
        lp_fmpz_mod_pow_si(base, base, exponent, parser->ctx);
        lp_skip_ws(&parser->cursor);
    }

    fmpz_set(result, base);
    fmpz_clear(base);
}

static void lp_eval_term(lp_eval_parser_t *parser, fmpz_t result)
{
    fmpz_t rhs;

    fmpz_init(rhs);
    lp_eval_power(parser, result);
    lp_skip_ws(&parser->cursor);

    while (parser->ok && (*parser->cursor == '*' || *parser->cursor == '/')) {
        char op = *parser->cursor++;
        lp_eval_power(parser, rhs);
        if (!parser->ok)
            break;
        if (op == '*') {
            fmpz_mod_mul(result, result, rhs, parser->ctx);
        } else {
            if (fmpz_is_zero(rhs) || !fmpz_mod_is_invertible(rhs, parser->ctx)) {
                parser->ok = 0;
                break;
            }
            fmpz_mod_inv(rhs, rhs, parser->ctx);
            fmpz_mod_mul(result, result, rhs, parser->ctx);
        }
        lp_skip_ws(&parser->cursor);
    }

    fmpz_clear(rhs);
}

static void lp_eval_expression(lp_eval_parser_t *parser, fmpz_t result)
{
    fmpz_t rhs;

    fmpz_init(rhs);
    lp_eval_term(parser, result);
    lp_skip_ws(&parser->cursor);

    while (parser->ok && (*parser->cursor == '+' || *parser->cursor == '-')) {
        char op = *parser->cursor++;
        lp_eval_term(parser, rhs);
        if (!parser->ok)
            break;
        if (op == '+')
            fmpz_mod_add(result, result, rhs, parser->ctx);
        else
            fmpz_mod_sub(result, result, rhs, parser->ctx);
        lp_skip_ws(&parser->cursor);
    }

    fmpz_clear(rhs);
}

static int lp_evaluate_polynomial_mod(const char *poly_str,
                                      rational_variable_info_t *vars, slong num_vars,
                                      fmpz_t *values,
                                      const fmpz_mod_ctx_t ctx,
                                      fmpz_t result)
{
    lp_eval_parser_t parser;

    parser.cursor = poly_str;
    parser.vars = vars;
    parser.num_vars = num_vars;
    parser.values = values;
    parser.ctx = ctx;
    parser.ok = 1;

    lp_eval_expression(&parser, result);
    lp_skip_ws(&parser.cursor);
    if (*parser.cursor != '\0')
        parser.ok = 0;

    if (!parser.ok) {
        fmpz_zero(result);
        return 0;
    }
    return 1;
}

static int lp_polynomial_is_probably_zero(const char *poly_str,
                                          lp_solver_state_t *state,
                                          rational_variable_info_t *remaining_vars,
                                          slong num_remaining,
                                          fmpz_t *current_values)
{
    fmpz_t *trial_values;
    fmpz_t eval;
    int all_zero = 1;

    if (!poly_str)
        return -1;

    if (strcmp(poly_str, "0") == 0)
        return 1;

    trial_values = (fmpz_t *) malloc((size_t) (state->num_original_vars > 0 ? state->num_original_vars : 1) * sizeof(fmpz_t));
    if (!trial_values)
        return -1;

    for (slong i = 0; i < state->num_original_vars; i++) {
        fmpz_init(trial_values[i]);
        fmpz_set(trial_values[i], current_values[i]);
    }

    fmpz_init(eval);
    for (slong trial = 0; trial < (num_remaining > 0 ? 3 : 1); trial++) {
        for (slong i = 0; i < num_remaining; i++) {
            ulong sample = (ulong) ((trial + 2) * (i + 3) + 1);
            lp_set_small_sample_value(trial_values[remaining_vars[i].index], sample, state->sols->ctx);
        }

        if (!lp_evaluate_polynomial_mod(poly_str, state->original_vars, state->num_original_vars,
                                        trial_values, state->sols->ctx, eval)) {
            all_zero = -1;
            break;
        }
        if (!fmpz_is_zero(eval)) {
            all_zero = 0;
            break;
        }
    }

    fmpz_clear(eval);
    for (slong i = 0; i < state->num_original_vars; i++)
        fmpz_clear(trial_values[i]);
    free(trial_values);

    return all_zero;
}

static char *lp_substitute_variable_in_polynomial(const char *poly_str, const char *var_name, const fmpz_t value)
{
    char *value_str = fmpz_get_str(NULL, 10, value);
    size_t value_len = value_str ? strlen(value_str) : 1;
    size_t estimate = strlen(poly_str) + (value_len + 2) * 8 + 32;
    size_t var_len = strlen(var_name);
    char *result = (char *) malloc(estimate);
    const char *ptr = poly_str;

    if (!result) {
        if (value_str) flint_free(value_str);
        return NULL;
    }

    result[0] = '\0';
    while (*ptr) {
        if (strncmp(ptr, var_name, var_len) == 0 &&
            (ptr == poly_str || !lp_is_identifier_char(*(ptr - 1))) &&
            !lp_is_identifier_char(ptr[var_len])) {
            strcat(result, "(");
            strcat(result, value_str ? value_str : "0");
            strcat(result, ")");
            ptr += var_len;
        } else {
            size_t len = strlen(result);
            result[len] = *ptr;
            result[len + 1] = '\0';
            ptr++;
        }
    }

    if (value_str) flint_free(value_str);
    return result;
}

static int lp_duplicate_nonzero_equations(lp_solver_state_t *state,
                                          char **poly_strings, slong num_polys,
                                          rational_variable_info_t *remaining_vars, slong num_remaining,
                                          fmpz_t *current_values,
                                          char ***active_polys_out,
                                          slong *active_count_out)
{
    char **active = (char **) calloc((size_t) (num_polys > 0 ? num_polys : 1), sizeof(char *));
    slong active_count = 0;

    *active_polys_out = NULL;
    *active_count_out = 0;

    if (!active)
        return 0;

    for (slong i = 0; i < num_polys; i++) {
        int zero_status = lp_polynomial_is_probably_zero(poly_strings[i], state,
                                                         remaining_vars, num_remaining,
                                                         current_values);
        if (zero_status < 0) {
            lp_free_string_array(active, active_count);
            return 0;
        }
        if (!zero_status) {
            active[active_count] = strdup(poly_strings[i]);
            if (!active[active_count]) {
                lp_free_string_array(active, active_count);
                return 0;
            }
            active_count++;
        }
    }

    *active_polys_out = active;
    *active_count_out = active_count;
    return 1;
}

static int lp_reduce_system_by_root(lp_solver_state_t *state,
                                    char **poly_strings, slong num_polys,
                                    rational_variable_info_t *vars, slong num_vars,
                                    const fmpz_t root,
                                    fmpz_t *current_values,
                                    char ***reduced_polys_out, slong *reduced_count_out)
{
    char **reduced = (char **) calloc((size_t) (num_polys > 0 ? num_polys : 1), sizeof(char *));
    slong reduced_count = 0;

    *reduced_polys_out = NULL;
    *reduced_count_out = 0;

    if (!reduced)
        return 0;

    for (slong i = 0; i < num_polys; i++) {
        char *substituted = lp_substitute_variable_in_polynomial(poly_strings[i], vars[0].name, root);
        int zero_status;

        if (!substituted) {
            lp_free_string_array(reduced, reduced_count);
            return 0;
        }

        zero_status = lp_polynomial_is_probably_zero(substituted, state, vars + 1, num_vars - 1, current_values);
        if (zero_status < 0) {
            free(substituted);
            lp_free_string_array(reduced, reduced_count);
            return 0;
        }
        if (zero_status) {
            free(substituted);
        } else {
            reduced[reduced_count++] = substituted;
        }
    }

    *reduced_polys_out = reduced;
    *reduced_count_out = reduced_count;
    return 1;
}

static int lp_join_strings_by_indices(char **strings, const slong *indices, slong count, char **joined_out)
{
    size_t total = 1;
    char *joined;

    *joined_out = NULL;

    for (slong i = 0; i < count; i++) {
        total += strlen(strings[indices[i]]) + 2;
    }

    joined = (char *) malloc(total);
    if (!joined)
        return 0;

    joined[0] = '\0';
    for (slong i = 0; i < count; i++) {
        if (i > 0)
            strcat(joined, ",");
        strcat(joined, strings[indices[i]]);
    }

    *joined_out = joined;
    return 1;
}

static int lp_join_elimination_vars(rational_variable_info_t *vars, slong num_vars, char **joined_out)
{
    size_t total = 1;
    char *joined;

    *joined_out = NULL;

    for (slong i = 1; i < num_vars; i++)
        total += strlen(vars[i].name) + 2;

    joined = (char *) malloc(total);
    if (!joined)
        return 0;

    joined[0] = '\0';
    for (slong i = 1; i < num_vars; i++) {
        if (i > 1)
            strcat(joined, ",");
        strcat(joined, vars[i].name);
    }

    *joined_out = joined;
    return 1;
}

static void lp_root_list_init(lp_root_list_t *roots)
{
    roots->values = NULL;
    roots->mults = NULL;
    roots->count = 0;
    roots->alloc = 0;
}

static void lp_root_list_clear(lp_root_list_t *roots)
{
    if (!roots)
        return;

    for (slong i = 0; i < roots->count; i++)
        fmpz_clear(roots->values[i]);
    free(roots->values);
    free(roots->mults);
    roots->values = NULL;
    roots->mults = NULL;
    roots->count = 0;
    roots->alloc = 0;
}

static int lp_root_list_append(lp_root_list_t *roots, const fmpz_t value, slong mult)
{
    for (slong i = 0; i < roots->count; i++) {
        if (fmpz_equal(roots->values[i], value)) {
            roots->mults[i] += mult;
            return 1;
        }
    }

    if (roots->count >= roots->alloc) {
        slong new_alloc = (roots->alloc == 0) ? 4 : 2 * roots->alloc;
        fmpz_t *new_values = (fmpz_t *) realloc(roots->values, (size_t) new_alloc * sizeof(fmpz_t));
        slong *new_mults = (slong *) realloc(roots->mults, (size_t) new_alloc * sizeof(slong));
        if (!new_values || !new_mults) {
            if (new_values && new_values != roots->values)
                roots->values = new_values;
            if (new_mults && new_mults != roots->mults)
                roots->mults = new_mults;
            return 0;
        }
        roots->values = new_values;
        roots->mults = new_mults;
        roots->alloc = new_alloc;
    }

    fmpz_init(roots->values[roots->count]);
    fmpz_set(roots->values[roots->count], value);
    roots->mults[roots->count] = mult;
    roots->count++;
    return 1;
}

static void lp_univar_expression(lp_univar_parser_t *parser, fmpz_mod_poly_t result);

static void lp_univar_primary(lp_univar_parser_t *parser, fmpz_mod_poly_t result)
{
    lp_skip_ws(&parser->cursor);
    fmpz_mod_poly_zero(result, parser->ctx);

    if (!parser->ok || *parser->cursor == '\0') {
        parser->ok = 0;
        return;
    }

    if (*parser->cursor == '(') {
        parser->cursor++;
        lp_univar_expression(parser, result);
        lp_skip_ws(&parser->cursor);
        if (*parser->cursor == ')')
            parser->cursor++;
        else
            parser->ok = 0;
        return;
    }

    if (*parser->cursor == '+' || *parser->cursor == '-') {
        char op = *parser->cursor++;
        lp_univar_primary(parser, result);
        if (op == '-')
            fmpz_mod_poly_neg(result, result, parser->ctx);
        return;
    }

    if (isdigit((unsigned char) *parser->cursor)) {
        fmpz_t tmp;
        fmpz_init(tmp);
        if (!lp_parse_integer_token(&parser->cursor, tmp)) {
            parser->ok = 0;
        } else {
            fmpz_mod_poly_set_coeff_fmpz(result, 0, tmp, parser->ctx);
        }
        fmpz_clear(tmp);
        return;
    }

    if (isalpha((unsigned char) *parser->cursor) || *parser->cursor == '_') {
        const char *start = parser->cursor;
        size_t len;
        while (lp_is_identifier_char(*parser->cursor))
            parser->cursor++;
        len = (size_t) (parser->cursor - start);
        if (strlen(parser->var_name) == len && strncmp(start, parser->var_name, len) == 0) {
            fmpz_t one;
            fmpz_init(one);
            fmpz_set_ui(one, 1);
            fmpz_mod_poly_set_coeff_fmpz(result, 1, one, parser->ctx);
            fmpz_clear(one);
        } else {
            parser->ok = 0;
        }
        return;
    }

    parser->ok = 0;
}

static void lp_univar_power(lp_univar_parser_t *parser, fmpz_mod_poly_t result)
{
    fmpz_mod_poly_t base, tmp;

    fmpz_mod_poly_init(base, parser->ctx);
    fmpz_mod_poly_init(tmp, parser->ctx);
    lp_univar_primary(parser, base);
    lp_skip_ws(&parser->cursor);

    while (parser->ok && *parser->cursor == '^') {
        char *endptr = NULL;
        long exponent;
        parser->cursor++;
        lp_skip_ws(&parser->cursor);
        exponent = strtol(parser->cursor, &endptr, 10);
        if (endptr == parser->cursor || exponent < 0) {
            parser->ok = 0;
            break;
        }
        parser->cursor = endptr;
        fmpz_mod_poly_pow(tmp, base, (ulong) exponent, parser->ctx);
        fmpz_mod_poly_swap(base, tmp, parser->ctx);
        lp_skip_ws(&parser->cursor);
    }

    fmpz_mod_poly_set(result, base, parser->ctx);
    fmpz_mod_poly_clear(tmp, parser->ctx);
    fmpz_mod_poly_clear(base, parser->ctx);
}

static void lp_univar_term(lp_univar_parser_t *parser, fmpz_mod_poly_t result)
{
    fmpz_mod_poly_t rhs, tmp, scalar;

    fmpz_mod_poly_init(rhs, parser->ctx);
    fmpz_mod_poly_init(tmp, parser->ctx);
    fmpz_mod_poly_init(scalar, parser->ctx);
    lp_univar_power(parser, result);
    lp_skip_ws(&parser->cursor);

    while (parser->ok && (*parser->cursor == '*' || *parser->cursor == '/')) {
        char op = *parser->cursor++;
        lp_univar_power(parser, rhs);
        if (!parser->ok)
            break;
        if (op == '*') {
            fmpz_mod_poly_mul(tmp, result, rhs, parser->ctx);
            fmpz_mod_poly_swap(result, tmp, parser->ctx);
        } else {
            fmpz_t c0, inv;
            if (fmpz_mod_poly_degree(rhs, parser->ctx) != 0) {
                parser->ok = 0;
                break;
            }
            fmpz_init(c0);
            fmpz_init(inv);
            fmpz_mod_poly_get_coeff_fmpz(c0, rhs, 0, parser->ctx);
            if (fmpz_is_zero(c0) || !fmpz_mod_is_invertible(c0, parser->ctx)) {
                parser->ok = 0;
                fmpz_clear(inv);
                fmpz_clear(c0);
                break;
            }
            fmpz_mod_inv(inv, c0, parser->ctx);
            fmpz_mod_poly_zero(scalar, parser->ctx);
            fmpz_mod_poly_set_coeff_fmpz(scalar, 0, inv, parser->ctx);
            fmpz_mod_poly_mul(tmp, result, scalar, parser->ctx);
            fmpz_mod_poly_swap(result, tmp, parser->ctx);
            fmpz_clear(inv);
            fmpz_clear(c0);
        }
        lp_skip_ws(&parser->cursor);
    }

    fmpz_mod_poly_clear(scalar, parser->ctx);
    fmpz_mod_poly_clear(tmp, parser->ctx);
    fmpz_mod_poly_clear(rhs, parser->ctx);
}

static void lp_univar_expression(lp_univar_parser_t *parser, fmpz_mod_poly_t result)
{
    fmpz_mod_poly_t rhs, tmp;

    fmpz_mod_poly_init(rhs, parser->ctx);
    fmpz_mod_poly_init(tmp, parser->ctx);
    lp_univar_term(parser, result);
    lp_skip_ws(&parser->cursor);

    while (parser->ok && (*parser->cursor == '+' || *parser->cursor == '-')) {
        char op = *parser->cursor++;
        lp_univar_term(parser, rhs);
        if (!parser->ok)
            break;
        if (op == '+')
            fmpz_mod_poly_add(tmp, result, rhs, parser->ctx);
        else
            fmpz_mod_poly_sub(tmp, result, rhs, parser->ctx);
        fmpz_mod_poly_swap(result, tmp, parser->ctx);
        lp_skip_ws(&parser->cursor);
    }

    fmpz_mod_poly_clear(tmp, parser->ctx);
    fmpz_mod_poly_clear(rhs, parser->ctx);
}

static int lp_resultant_string_to_poly(const char *resultant_str,
                                       const char *var_name,
                                       const fmpz_t prime,
                                       fmpz_mod_poly_t poly,
                                       const fmpz_mod_ctx_t modctx)
{
    lp_univar_parser_t parser;

    (void) prime;
    fmpz_mod_poly_zero(poly, modctx);
    parser.cursor = resultant_str;
    parser.var_name = var_name;
    parser.ctx = modctx;
    parser.ok = 1;

    lp_univar_expression(&parser, poly);
    lp_skip_ws(&parser.cursor);
    return parser.ok && *parser.cursor == '\0';
}

static int lp_resultant_roots(const char *resultant_str,
                              const char *var_name,
                              const fmpz_t prime,
                              lp_root_list_t *roots,
                              slong *degree_out)
{
    fmpz_mod_ctx_t ctx;
    fmpz_mod_poly_t poly, factor_poly;
    fmpz_mod_poly_factor_t factors;
    fmpz_t c0, c1, root, tmp;
    slong degree;
    int ok = 0;

    *degree_out = -1;
    fmpz_mod_ctx_init(ctx, prime);
    fmpz_mod_poly_init(poly, ctx);

    if (!lp_resultant_string_to_poly(resultant_str, var_name, prime, poly, ctx))
        goto cleanup_poly;

    degree = fmpz_mod_poly_degree(poly, ctx);
    *degree_out = degree;
    if (degree <= 0) {
        ok = 1;
        goto cleanup_poly;
    }

    fmpz_mod_poly_factor_init(factors, ctx);
    fmpz_mod_poly_roots(factors, poly, 1, ctx);
    fmpz_mod_poly_init(factor_poly, ctx);
    fmpz_init(c0);
    fmpz_init(c1);
    fmpz_init(root);
    fmpz_init(tmp);

    for (slong i = 0; i < factors->num; i++) {
        fmpz_mod_poly_factor_get_poly(factor_poly, factors, i, ctx);
        if (fmpz_mod_poly_degree(factor_poly, ctx) != 1)
            continue;

        fmpz_mod_poly_get_coeff_fmpz(c0, factor_poly, 0, ctx);
        fmpz_mod_poly_get_coeff_fmpz(c1, factor_poly, 1, ctx);
        if (!fmpz_mod_is_invertible(c1, ctx))
            continue;

        fmpz_mod_inv(tmp, c1, ctx);
        fmpz_mod_neg(root, c0, ctx);
        fmpz_mod_mul(root, root, tmp, ctx);

        if (!lp_root_list_append(roots, root, factors->exp[i]))
            goto cleanup_factors;
    }

    ok = 1;

cleanup_factors:
    fmpz_clear(tmp);
    fmpz_clear(root);
    fmpz_clear(c1);
    fmpz_clear(c0);
    fmpz_mod_poly_clear(factor_poly, ctx);
    fmpz_mod_poly_factor_clear(factors, ctx);
cleanup_poly:
    fmpz_mod_poly_clear(poly, ctx);
    fmpz_mod_ctx_clear(ctx);
    return ok;
}

static int lp_verify_full_assignment(lp_solver_state_t *state, fmpz_t *current_values)
{
    fmpz_t eval;
    int ok = 1;

    fmpz_init(eval);
    for (slong i = 0; i < state->num_original_polys; i++) {
        if (!lp_evaluate_polynomial_mod(state->original_polys[i],
                                        state->original_vars,
                                        state->num_original_vars,
                                        current_values,
                                        state->sols->ctx,
                                        eval)) {
            ok = 0;
            break;
        }
        if (!fmpz_is_zero(eval))
            ok = 0;
        if (!ok)
            break;
    }
    fmpz_clear(eval);
    return ok;
}

static int lp_append_current_solution(lp_solver_state_t *state, fmpz_t *current_values)
{
    large_prime_solutions_t *sols = state->sols;
    fmpz_t **new_rows;
    slong idx;

    if (!lp_verify_full_assignment(state, current_values))
        return 1;

    for (slong s = 0; s < sols->num_solution_sets; s++) {
        int same = 1;
        for (slong v = 0; v < sols->num_variables; v++) {
            if (!fmpz_equal(sols->solution_values[s][v], current_values[v])) {
                same = 0;
                break;
            }
        }
        if (same)
            return 1;
    }

    new_rows = (fmpz_t **) realloc(sols->solution_values,
                                   (size_t) (sols->num_solution_sets + 1) * sizeof(fmpz_t *));
    if (!new_rows)
        return 0;
    sols->solution_values = new_rows;

    idx = sols->num_solution_sets++;
    sols->solution_values[idx] = (fmpz_t *) malloc((size_t) sols->num_variables * sizeof(fmpz_t));
    if (!sols->solution_values[idx])
        return 0;

    for (slong v = 0; v < sols->num_variables; v++) {
        fmpz_init(sols->solution_values[idx][v]);
        fmpz_set(sols->solution_values[idx][v], current_values[v]);
    }

    return 1;
}

static int lp_solve_recursive(lp_solver_state_t *state,
                              char **poly_strings, slong num_polys,
                              rational_variable_info_t *vars, slong num_vars,
                              fmpz_t *current_values)
{
    char **active_polys = NULL;
    slong active_count = 0;
    int status = LP_SOLVE_NONE;

    lp_progress("large-prime recursion: vars=%ld polys=%ld", num_vars, num_polys);

    if (!lp_duplicate_nonzero_equations(state, poly_strings, num_polys,
                                        vars, num_vars, current_values,
                                        &active_polys, &active_count)) {
        lp_set_error(state->sols, "Failed to analyze a large-prime equation");
        return LP_SOLVE_ERROR;
    }

    if (active_count == 0) {
        lp_free_string_array(active_polys, active_count);
        if (num_vars == 0) {
            if (!lp_append_current_solution(state, current_values)) {
                lp_set_error(state->sols, "Failed to append a large-prime solution");
                return LP_SOLVE_ERROR;
            }
            return LP_SOLVE_FOUND;
        }
        return LP_SOLVE_DIMPOS;
    }

    if (num_vars == 0) {
        lp_free_string_array(active_polys, active_count);
        return LP_SOLVE_NONE;
    }

    if (active_count < num_vars) {
        lp_free_string_array(active_polys, active_count);
        return LP_SOLVE_DIMPOS;
    }

    {
        rational_equation_combination_t *combinations = NULL;
        slong num_combinations = 0;
        slong target_equations = num_vars;
        int found_solution = 0;
        int found_dimpos = 0;

        if (active_count == target_equations) {
            combinations = (rational_equation_combination_t *) malloc(sizeof(rational_equation_combination_t));
            if (!combinations) {
                lp_free_string_array(active_polys, active_count);
                lp_set_error(state->sols, "Out of memory while preparing equation combinations");
                return LP_SOLVE_ERROR;
            }
            combinations[0].equation_indices = (slong *) malloc((size_t) target_equations * sizeof(slong));
            if (!combinations[0].equation_indices) {
                free(combinations);
                lp_free_string_array(active_polys, active_count);
                lp_set_error(state->sols, "Out of memory while preparing equation combinations");
                return LP_SOLVE_ERROR;
            }
            combinations[0].num_equations = target_equations;
            combinations[0].tried = 0;
            combinations[0].success = 0;
            for (slong i = 0; i < target_equations; i++)
                combinations[0].equation_indices[i] = i;
            num_combinations = 1;
        } else {
            rational_generate_equation_combinations(active_count, target_equations,
                                                    &combinations, &num_combinations);
        }

        for (slong comb_idx = 0; comb_idx < num_combinations && !found_solution; comb_idx++) {
            char *selected_polys = NULL;
            char *elim_vars = NULL;
            char *resultant = NULL;
            lp_root_list_t roots;
            slong resultant_degree = -1;

            if (!lp_join_strings_by_indices(active_polys,
                                            combinations[comb_idx].equation_indices,
                                            target_equations,
                                            &selected_polys) ||
                (num_vars > 1 && !lp_join_elimination_vars(vars, num_vars, &elim_vars))) {
                free(selected_polys);
                free(elim_vars);
                rational_free_equation_combinations(combinations, num_combinations);
                lp_free_string_array(active_polys, active_count);
                lp_set_error(state->sols, "Out of memory while building a Dixon elimination task");
                return LP_SOLVE_ERROR;
            }

            if (num_vars == 1)
                resultant = strdup(selected_polys);
            else
                resultant = dixon_str_large_prime(selected_polys, elim_vars, state->sols->prime);
            free(selected_polys);
            free(elim_vars);

            if (!resultant) {
                continue;
            }

            if (strcmp(resultant, "0") == 0) {
                found_dimpos = 1;
                free(resultant);
                continue;
            }

            lp_root_list_init(&roots);
            if (!lp_resultant_roots(resultant, vars[0].name, state->sols->prime, &roots, &resultant_degree)) {
                free(resultant);
                lp_root_list_clear(&roots);
                rational_free_equation_combinations(combinations, num_combinations);
                lp_free_string_array(active_polys, active_count);
                lp_set_error(state->sols, "Failed to parse a large-prime resultant or extract its roots");
                return LP_SOLVE_ERROR;
            }
            free(resultant);

            if (resultant_degree <= 0 || roots.count == 0) {
                lp_root_list_clear(&roots);
                continue;
            }

            for (slong r = 0; r < roots.count; r++) {
                char **reduced_polys = NULL;
                slong reduced_count = 0;
                int sub_status;

                fmpz_set(current_values[vars[0].index], roots.values[r]);

                if (!lp_reduce_system_by_root(state, active_polys, active_count, vars, num_vars,
                                             roots.values[r], current_values,
                                             &reduced_polys, &reduced_count)) {
                    lp_root_list_clear(&roots);
                    rational_free_equation_combinations(combinations, num_combinations);
                    lp_free_string_array(active_polys, active_count);
                    lp_set_error(state->sols, "Failed to substitute a large-prime root back into the system");
                    return LP_SOLVE_ERROR;
                }

                if (reduced_count == 0) {
                    if (num_vars - 1 == 0) {
                        if (!lp_append_current_solution(state, current_values)) {
                            lp_free_string_array(reduced_polys, reduced_count);
                            lp_root_list_clear(&roots);
                            rational_free_equation_combinations(combinations, num_combinations);
                            lp_free_string_array(active_polys, active_count);
                            lp_set_error(state->sols, "Failed to append a large-prime solution");
                            return LP_SOLVE_ERROR;
                        }
                        found_solution = 1;
                    } else {
                        found_dimpos = 1;
                    }
                    lp_free_string_array(reduced_polys, reduced_count);
                    continue;
                }

                sub_status = lp_solve_recursive(state, reduced_polys, reduced_count,
                                                vars + 1, num_vars - 1, current_values);
                lp_free_string_array(reduced_polys, reduced_count);

                if (sub_status == LP_SOLVE_ERROR) {
                    lp_root_list_clear(&roots);
                    rational_free_equation_combinations(combinations, num_combinations);
                    lp_free_string_array(active_polys, active_count);
                    return LP_SOLVE_ERROR;
                }
                if (sub_status == LP_SOLVE_FOUND)
                    found_solution = 1;
                else if (sub_status == LP_SOLVE_DIMPOS)
                    found_dimpos = 1;
            }

            lp_root_list_clear(&roots);
        }

        rational_free_equation_combinations(combinations, num_combinations);
        lp_free_string_array(active_polys, active_count);

        if (found_solution)
            status = LP_SOLVE_FOUND;
        else if (found_dimpos)
            status = LP_SOLVE_DIMPOS;
    }

    return status;
}

static void lp_print_solution_row(FILE *fp, const large_prime_solutions_t *sols, slong row)
{
    for (slong v = 0; v < sols->num_variables; v++) {
        if (v > 0)
            fprintf(fp, ", ");
        fprintf(fp, "%s = ", sols->variable_names[v]);
        fmpz_fprint(fp, sols->solution_values[row][v]);
    }
}

void large_prime_solutions_clear(large_prime_solutions_t *sols)
{
    if (!sols)
        return;

    if (sols->variable_names) {
        for (slong i = 0; i < sols->num_variables; i++)
            free(sols->variable_names[i]);
        free(sols->variable_names);
    }

    if (sols->solution_values) {
        for (slong s = 0; s < sols->num_solution_sets; s++) {
            if (sols->solution_values[s]) {
                for (slong v = 0; v < sols->num_variables; v++)
                    fmpz_clear(sols->solution_values[s][v]);
                free(sols->solution_values[s]);
            }
        }
        free(sols->solution_values);
    }

    if (sols->error_message)
        free(sols->error_message);

    fmpz_clear(sols->prime);
    fmpz_mod_ctx_clear(sols->ctx);
}

void print_large_prime_solutions(const large_prime_solutions_t *sols)
{
    if (!sols) {
        printf("Error: no large-prime solver result available.\n");
        return;
    }

    if (!sols->is_valid) {
        printf("Error: %s\n", sols->error_message ? sols->error_message : "large-prime solving failed");
        return;
    }

    if (sols->has_no_solutions == -1) {
        printf("System appears to have dimension > 0 over F_");
        fmpz_print(sols->prime);
        printf(".\n");
        return;
    }

    if (sols->has_no_solutions || sols->num_solution_sets == 0) {
        printf("No solutions found over F_");
        fmpz_print(sols->prime);
        printf(".\n");
        return;
    }

    printf("Found %ld solution set(s) over F_", sols->num_solution_sets);
    fmpz_print(sols->prime);
    printf(":\n");

    for (slong s = 0; s < sols->num_solution_sets; s++) {
        printf("  Solution %ld: ", s + 1);
        lp_print_solution_row(stdout, sols, s);
        printf("\n");
    }
}

void save_large_prime_solver_result_to_file(const char *filename,
                                            const char *polys_str,
                                            const fmpz_t prime,
                                            const large_prime_solutions_t *sols,
                                            double cpu_time,
                                            double wall_time,
                                            int threads_num)
{
    FILE *fp = fopen(filename, "w");
    if (!fp)
        return;

    fprintf(fp, "Large-prime polynomial system solver result\n");
    fprintf(fp, "==========================================\n");
    fprintf(fp, "Field: F_");
    fmpz_fprint(fp, prime);
    fprintf(fp, "\n");
    fprintf(fp, "Polynomials: %s\n", polys_str);
    fprintf(fp, "CPU time: %.3f seconds\n", cpu_time);
    fprintf(fp, "Wall time: %.3f seconds\n", wall_time);
    fprintf(fp, "Threads: %d\n\n", threads_num);

    if (!sols) {
        fprintf(fp, "Error: no solver result available.\n");
    } else if (!sols->is_valid) {
        fprintf(fp, "Error: %s\n", sols->error_message ? sols->error_message : "large-prime solving failed");
    } else if (sols->has_no_solutions == -1) {
        fprintf(fp, "System appears to have dimension > 0 over this field.\n");
    } else if (sols->has_no_solutions || sols->num_solution_sets == 0) {
        fprintf(fp, "No solutions found.\n");
    } else {
        fprintf(fp, "Found %ld solution set(s):\n", sols->num_solution_sets);
        for (slong s = 0; s < sols->num_solution_sets; s++) {
            fprintf(fp, "  Solution %ld: ", s + 1);
            lp_print_solution_row(fp, sols, s);
            fprintf(fp, "\n");
        }
    }

    fclose(fp);
}

static char **lp_compute_remaining_vars(const char *polys_str, const char *elim_str, slong *num_remaining_out)
{
    slong num_polys = 0;
    slong num_all_vars = 0;
    slong num_elim = 0;
    char **poly_array = NULL;
    char **all_vars = NULL;
    char **elim_vars = NULL;
    char **remaining = NULL;
    slong num_remaining = 0;

    *num_remaining_out = 0;

    poly_array = split_string(polys_str, &num_polys);
    if (!poly_array || num_polys == 0) {
        free_split_strings(poly_array, num_polys);
        return NULL;
    }

    all_vars = rational_extract_variables_improved(poly_array, num_polys, &num_all_vars);
    free_split_strings(poly_array, num_polys);
    if (!all_vars)
        return NULL;

    elim_vars = split_string(elim_str, &num_elim);
    remaining = (char **) calloc((size_t) (num_all_vars > 0 ? num_all_vars : 1), sizeof(char *));
    if (!remaining) {
        free_split_strings(all_vars, num_all_vars);
        free_split_strings(elim_vars, num_elim);
        return NULL;
    }

    for (slong i = 0; i < num_all_vars; i++) {
        if (!lp_is_name_in_list(all_vars[i], elim_vars, num_elim)) {
            remaining[num_remaining++] = strdup(all_vars[i]);
            if (!remaining[num_remaining - 1]) {
                lp_free_string_array(remaining, num_remaining - 1);
                free_split_strings(all_vars, num_all_vars);
                free_split_strings(elim_vars, num_elim);
                return NULL;
            }
        }
    }

    free_split_strings(all_vars, num_all_vars);
    free_split_strings(elim_vars, num_elim);
    *num_remaining_out = num_remaining;
    return remaining;
}

void large_prime_print_roots_from_resultant_string(const char *resultant_str,
                                                   const char *polys_str,
                                                   const char *vars_str,
                                                   const fmpz_t prime,
                                                   FILE *fp,
                                                   int print_to_stdout)
{
    char **remaining_vars = NULL;
    slong num_remaining = 0;
    lp_root_list_t roots;
    slong degree = -1;

    if (!resultant_str || strcmp(resultant_str, "0") == 0)
        return;

    remaining_vars = lp_compute_remaining_vars(polys_str, vars_str, &num_remaining);
    if (!remaining_vars || num_remaining != 1) {
        lp_free_string_array(remaining_vars, num_remaining);
        return;
    }

    lp_root_list_init(&roots);
    if (!lp_resultant_roots(resultant_str, remaining_vars[0], prime, &roots, &degree)) {
        lp_root_list_clear(&roots);
        lp_free_string_array(remaining_vars, num_remaining);
        return;
    }

    if (degree > 0 && roots.count > 0) {
        if (print_to_stdout) {
            printf("\nRoots in %s (degree %ld):\n", remaining_vars[0], degree);
            printf("Find %ld roots:\n", roots.count);
            for (slong i = 0; i < roots.count; i++) {
                printf("  Root %ld: ", i + 1);
                fmpz_print(roots.values[i]);
                printf(" (Multiplicity: %ld)\n", roots.mults[i]);
            }
        }
        if (fp) {
            fprintf(fp, "\nRoots in %s (degree %ld):\n", remaining_vars[0], degree);
            fprintf(fp, "Find %ld roots:\n", roots.count);
            for (slong i = 0; i < roots.count; i++) {
                fprintf(fp, "  Root %ld: ", i + 1);
                fmpz_fprint(fp, roots.values[i]);
                fprintf(fp, " (Multiplicity: %ld)\n", roots.mults[i]);
            }
        }
    }

    lp_root_list_clear(&roots);
    lp_free_string_array(remaining_vars, num_remaining);
}

large_prime_solutions_t *solve_large_prime_polynomial_system_string(const char *poly_string,
                                                                    const fmpz_t prime)
{
    slong num_polys = 0;
    char **poly_strings = split_string(poly_string, &num_polys);
    rational_variable_info_t *vars = NULL;
    slong num_vars = 0;
    large_prime_solutions_t *sols;
    lp_solver_state_t state;
    fmpz_t *current_values = NULL;
    int status;

    if (!poly_strings || num_polys == 0)
        return NULL;

    if (!rational_extract_and_sort_variables(poly_strings, num_polys, &vars, &num_vars)) {
        free_split_strings(poly_strings, num_polys);
        return NULL;
    }

    sols = (large_prime_solutions_t *) malloc(sizeof(large_prime_solutions_t));
    if (!sols) {
        lp_free_rational_vars(vars, num_vars);
        free_split_strings(poly_strings, num_polys);
        return NULL;
    }

    lp_solutions_init(sols, num_vars, prime);
    sols->num_equations = num_polys;

    for (slong i = 0; i < num_vars; i++) {
        sols->variable_names[i] = strdup(vars[i].name);
        if (!sols->variable_names[i]) {
            lp_set_error(sols, "Out of memory while storing large-prime variable names");
            free_split_strings(poly_strings, num_polys);
            lp_free_rational_vars(vars, num_vars);
            return sols;
        }
    }

    current_values = (fmpz_t *) malloc((size_t) (num_vars > 0 ? num_vars : 1) * sizeof(fmpz_t));
    if (!current_values) {
        lp_set_error(sols, "Out of memory while preparing large-prime current values");
        free_split_strings(poly_strings, num_polys);
        lp_free_rational_vars(vars, num_vars);
        return sols;
    }
    for (slong i = 0; i < num_vars; i++)
        fmpz_init(current_values[i]);

    state.sols = sols;
    state.original_polys = poly_strings;
    state.num_original_polys = num_polys;
    state.original_vars = vars;
    state.num_original_vars = num_vars;

    status = lp_solve_recursive(&state, poly_strings, num_polys, vars, num_vars, current_values);
    if (status == LP_SOLVE_ERROR) {
        sols->is_valid = 0;
    } else {
        sols->is_valid = 1;
        if (status == LP_SOLVE_FOUND) {
            sols->has_no_solutions = 0;
        } else if (status == LP_SOLVE_DIMPOS) {
            sols->has_no_solutions = -1;
        } else {
            sols->has_no_solutions = 1;
        }
    }

    for (slong i = 0; i < num_vars; i++)
        fmpz_clear(current_values[i]);
    free(current_values);

    lp_free_rational_vars(vars, num_vars);
    free_split_strings(poly_strings, num_polys);
    return sols;
}
