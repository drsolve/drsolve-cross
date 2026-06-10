#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "complex_solver.h"
#include "dixon_interface_flint.h"
#include "rational_system_solver.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

static int complex_solver_roots_reserve(complex_solver_roots_t *roots, slong need)
{
    if (!roots || need <= roots->alloc) {
        return 1;
    }

    slong new_alloc = roots->alloc ? roots->alloc : 4;
    while (new_alloc < need) {
        new_alloc *= 2;
    }

    acb_t *new_roots = (acb_t *) realloc(roots->roots, (size_t) new_alloc * sizeof(acb_t));
    slong *new_mults = (slong *) realloc(roots->multiplicities, (size_t) new_alloc * sizeof(slong));
    if (!new_roots || !new_mults) {
        if (new_roots) roots->roots = new_roots;
        if (new_mults) roots->multiplicities = new_mults;
        return 0;
    }

    roots->roots = new_roots;
    roots->multiplicities = new_mults;
    for (slong i = roots->alloc; i < new_alloc; i++) {
        acb_init(roots->roots[i]);
        roots->multiplicities[i] = 0;
    }
    roots->alloc = new_alloc;
    return 1;
}

void complex_solver_roots_init(complex_solver_roots_t *roots)
{
    if (!roots) {
        return;
    }

    roots->roots = NULL;
    roots->multiplicities = NULL;
    roots->num_roots = 0;
    roots->alloc = 0;
    roots->prec = 0;
    roots->source_poly = NULL;
    roots->variable_name = NULL;
}

void complex_solver_roots_clear(complex_solver_roots_t *roots)
{
    if (!roots) {
        return;
    }

    if (roots->roots) {
        for (slong i = 0; i < roots->alloc; i++) {
            acb_clear(roots->roots[i]);
        }
        free(roots->roots);
    }
    free(roots->multiplicities);
    free(roots->source_poly);
    free(roots->variable_name);

    roots->roots = NULL;
    roots->multiplicities = NULL;
    roots->num_roots = 0;
    roots->alloc = 0;
    roots->prec = 0;
    roots->source_poly = NULL;
    roots->variable_name = NULL;
}

void complex_solver_roots_print(const complex_solver_roots_t *roots, FILE *fp, slong digits)
{
    FILE *out = fp ? fp : stdout;

    if (!roots) {
        fprintf(out, "No complex root data available.\n");
        return;
    }

    fprintf(out, "Found %ld approximate complex root(s)", roots->num_roots);
    if (roots->variable_name) {
        fprintf(out, " for %s", roots->variable_name);
    }
    fprintf(out, ":\n");

    for (slong i = 0; i < roots->num_roots; i++) {
        fprintf(out, "  Root %ld: ", i + 1);
        acb_fprintd(out, roots->roots[i], digits);
        if (roots->multiplicities) {
            fprintf(out, " (Multiplicity: %ld)", roots->multiplicities[i]);
        }
        fprintf(out, "\n");
    }
}

static int complex_solver_solution_list_reserve(complex_solver_solution_list_t *sols, slong need)
{
    if (!sols || need <= sols->alloc) {
        return 1;
    }

    slong new_alloc = sols->alloc ? sols->alloc : 4;
    while (new_alloc < need) {
        new_alloc *= 2;
    }

    complex_solver_solution_pair_t *new_items =
        (complex_solver_solution_pair_t *) realloc(sols->items,
                                                   (size_t) new_alloc * sizeof(complex_solver_solution_pair_t));
    if (!new_items) {
        return 0;
    }

    sols->items = new_items;
    for (slong i = sols->alloc; i < new_alloc; i++) {
        acb_init(sols->items[i].x);
        acb_init(sols->items[i].y);
        acb_init(sols->items[i].residual1);
        acb_init(sols->items[i].residual2);
    }
    sols->alloc = new_alloc;
    return 1;
}

void complex_solver_solution_list_init(complex_solver_solution_list_t *sols)
{
    if (!sols) {
        return;
    }

    sols->items = NULL;
    sols->num_solutions = 0;
    sols->alloc = 0;
    sols->prec = 0;
    sols->source_system = NULL;
    sols->x_name = NULL;
    sols->y_name = NULL;
}

void complex_solver_solution_list_clear(complex_solver_solution_list_t *sols)
{
    if (!sols) {
        return;
    }

    if (sols->items) {
        for (slong i = 0; i < sols->alloc; i++) {
            acb_clear(sols->items[i].x);
            acb_clear(sols->items[i].y);
            acb_clear(sols->items[i].residual1);
            acb_clear(sols->items[i].residual2);
        }
        free(sols->items);
    }
    free(sols->source_system);
    free(sols->x_name);
    free(sols->y_name);

    sols->items = NULL;
    sols->num_solutions = 0;
    sols->alloc = 0;
    sols->prec = 0;
    sols->source_system = NULL;
    sols->x_name = NULL;
    sols->y_name = NULL;
}

void complex_solver_solution_list_print(const complex_solver_solution_list_t *sols, FILE *fp, slong digits)
{
    FILE *out = fp ? fp : stdout;
    const char *x_name = (sols && sols->x_name) ? sols->x_name : "x";
    const char *y_name = (sols && sols->y_name) ? sols->y_name : "y";

    if (!sols) {
        fprintf(out, "No complex 2x2 solution data available.\n");
        return;
    }

    fprintf(out, "Found %ld approximate complex solution set(s):\n", sols->num_solutions);
    for (slong i = 0; i < sols->num_solutions; i++) {
        fprintf(out, "  [%ld] %s = ", i + 1, x_name);
        acb_fprintd(out, sols->items[i].x, digits);
        fprintf(out, ", %s = ", y_name);
        acb_fprintd(out, sols->items[i].y, digits);
        fprintf(out, "\n");
        fprintf(out, "      residuals: f1 = ");
        acb_fprintd(out, sols->items[i].residual1, digits);
        fprintf(out, ", f2 = ");
        acb_fprintd(out, sols->items[i].residual2, digits);
        fprintf(out, "\n");
    }
}

static int complex_solver_set_root_metadata(complex_solver_roots_t *roots,
                                            const char *poly_str,
                                            const char *var_name,
                                            slong prec)
{
    char *poly_copy = NULL;
    char *var_copy = NULL;

    if (!roots) {
        return 0;
    }

    if (poly_str) {
        poly_copy = strdup(poly_str);
        if (!poly_copy) {
            return 0;
        }
    }

    if (var_name) {
        var_copy = strdup(var_name);
        if (!var_copy) {
            free(poly_copy);
            return 0;
        }
    }

    free(roots->source_poly);
    free(roots->variable_name);
    roots->source_poly = poly_copy;
    roots->variable_name = var_copy;
    roots->prec = prec;
    return 1;
}

int complex_solver_parse_univariate_fmpq_poly(const char *poly_str,
                                              const char *var_name,
                                              fmpq_poly_t poly)
{
    const char *ptr;
    size_t var_len;

    if (!poly_str || !var_name) {
        return 0;
    }

    fmpq_poly_init(poly);
    fmpq_poly_zero(poly);

    if (strcmp(poly_str, "0") == 0 || strcmp(poly_str, "") == 0) {
        return 1;
    }

    ptr = poly_str;
    var_len = strlen(var_name);

    while (*ptr) {
        int sign = 1;
        slong term_degree = 0;
        int parsed_factor = 0;
        fmpq_t term_coeff;

        while (isspace((unsigned char) *ptr)) ptr++;
        while (*ptr == ')') {
            ptr++;
            while (isspace((unsigned char) *ptr)) ptr++;
        }
        if (*ptr == '\0') break;

        if (*ptr == '+') {
            ptr++;
        } else if (*ptr == '-') {
            sign = -1;
            ptr++;
        }
        while (isspace((unsigned char) *ptr)) ptr++;

        fmpq_init(term_coeff);
        fmpq_one(term_coeff);

        while (*ptr) {
            while (isspace((unsigned char) *ptr)) ptr++;

            if (*ptr == ')') {
                ptr++;
                continue;
            }

            if (*ptr == '(') {
                const char *factor_ptr = ptr + 1;
                int factor_sign = 1;

                while (isspace((unsigned char) *factor_ptr)) factor_ptr++;
                if (*factor_ptr == '+') {
                    factor_ptr++;
                    while (isspace((unsigned char) *factor_ptr)) factor_ptr++;
                } else if (*factor_ptr == '-') {
                    factor_sign = -1;
                    factor_ptr++;
                    while (isspace((unsigned char) *factor_ptr)) factor_ptr++;
                }

                if (strncmp(factor_ptr, var_name, var_len) == 0 &&
                    !isalnum((unsigned char) factor_ptr[var_len])) {
                    slong factor_degree = 1;
                    factor_ptr += var_len;
                    while (isspace((unsigned char) *factor_ptr)) factor_ptr++;
                    if (*factor_ptr == '^') {
                        factor_ptr++;
                        while (isspace((unsigned char) *factor_ptr)) factor_ptr++;
                        if (isdigit((unsigned char) *factor_ptr)) {
                            factor_degree = strtol(factor_ptr, (char **) &factor_ptr, 10);
                        }
                    }
                    while (isspace((unsigned char) *factor_ptr)) factor_ptr++;
                    if (*factor_ptr == ')') {
                        factor_ptr++;
                        term_degree += factor_degree;
                        if (factor_sign < 0) {
                            fmpq_neg(term_coeff, term_coeff);
                        }
                        ptr = factor_ptr;
                        parsed_factor = 1;
                    } else {
                        break;
                    }
                } else if (isdigit((unsigned char) *factor_ptr) || *factor_ptr == '.') {
                    const char *num_start = factor_ptr;
                    size_t num_len;
                    char *num_str;
                    char *slash;
                    fmpq_t factor_coeff;

                    while (isdigit((unsigned char) *factor_ptr) || *factor_ptr == '.' || *factor_ptr == '/') {
                        factor_ptr++;
                    }
                    while (isspace((unsigned char) *factor_ptr)) factor_ptr++;
                    if (*factor_ptr != ')') {
                        break;
                    }

                    num_len = (size_t) (factor_ptr - num_start);
                    while (num_len > 0 && isspace((unsigned char) num_start[num_len - 1])) num_len--;
                    num_str = (char *) malloc(num_len + 1);
                    if (!num_str) {
                        fmpq_clear(term_coeff);
                        return 0;
                    }
                    memcpy(num_str, num_start, num_len);
                    num_str[num_len] = '\0';

                    fmpq_init(factor_coeff);
                    slash = strchr(num_str, '/');
                    if (slash) {
                        fmpz_t num, den;
                        *slash = '\0';
                        fmpz_init(num);
                        fmpz_init(den);
                        fmpz_set_str(num, num_str, 10);
                        fmpz_set_str(den, slash + 1, 10);
                        fmpq_set_fmpz_frac(factor_coeff, num, den);
                        fmpz_clear(num);
                        fmpz_clear(den);
                    } else {
                        fmpz_t num;
                        fmpz_init(num);
                        fmpz_set_str(num, num_str, 10);
                        fmpq_set_fmpz(factor_coeff, num);
                        fmpz_clear(num);
                    }

                    if (factor_sign < 0) {
                        fmpq_neg(factor_coeff, factor_coeff);
                    }
                    fmpq_mul(term_coeff, term_coeff, factor_coeff);
                    fmpq_clear(factor_coeff);
                    free(num_str);
                    ptr = factor_ptr + 1;
                    parsed_factor = 1;
                } else {
                    break;
                }
            } else if (isdigit((unsigned char) *ptr) || *ptr == '.') {
                const char *num_start = ptr;
                size_t num_len;
                char *num_str;
                char *slash;
                fmpq_t factor_coeff;

                while (isdigit((unsigned char) *ptr) || *ptr == '.' || *ptr == '/') ptr++;
                num_len = (size_t) (ptr - num_start);
                num_str = (char *) malloc(num_len + 1);
                if (!num_str) {
                    fmpq_clear(term_coeff);
                    return 0;
                }
                memcpy(num_str, num_start, num_len);
                num_str[num_len] = '\0';

                fmpq_init(factor_coeff);
                slash = strchr(num_str, '/');
                if (slash) {
                    fmpz_t num, den;
                    *slash = '\0';
                    fmpz_init(num);
                    fmpz_init(den);
                    fmpz_set_str(num, num_str, 10);
                    fmpz_set_str(den, slash + 1, 10);
                    fmpq_set_fmpz_frac(factor_coeff, num, den);
                    fmpz_clear(num);
                    fmpz_clear(den);
                } else {
                    fmpz_t num;
                    fmpz_init(num);
                    fmpz_set_str(num, num_str, 10);
                    fmpq_set_fmpz(factor_coeff, num);
                    fmpz_clear(num);
                }
                fmpq_mul(term_coeff, term_coeff, factor_coeff);
                fmpq_clear(factor_coeff);
                free(num_str);
                parsed_factor = 1;
            } else if (strncmp(ptr, var_name, var_len) == 0 &&
                       !isalnum((unsigned char) ptr[var_len])) {
                slong factor_degree = 1;
                ptr += var_len;
                while (isspace((unsigned char) *ptr)) ptr++;
                if (*ptr == '^') {
                    ptr++;
                    while (isspace((unsigned char) *ptr)) ptr++;
                    if (isdigit((unsigned char) *ptr)) {
                        factor_degree = strtol(ptr, (char **) &ptr, 10);
                    }
                }
                term_degree += factor_degree;
                parsed_factor = 1;
            } else {
                break;
            }

            while (isspace((unsigned char) *ptr)) ptr++;
            if (*ptr == '*') {
                ptr++;
                continue;
            }
            break;
        }

        if (!parsed_factor) {
            fmpq_clear(term_coeff);
            if (*ptr) {
                ptr++;
                continue;
            }
            break;
        }

        if (sign < 0) {
            fmpq_neg(term_coeff, term_coeff);
        }

        {
            fmpq_t existing_coeff;
            fmpq_init(existing_coeff);
            fmpq_poly_get_coeff_fmpq(existing_coeff, poly, term_degree);
            fmpq_add(existing_coeff, existing_coeff, term_coeff);
            fmpq_poly_set_coeff_fmpq(poly, term_degree, existing_coeff);
            fmpq_clear(existing_coeff);
        }

        fmpq_clear(term_coeff);
    }

    return 1;
}

int complex_solver_univariate_complex_roots_from_string(const char *poly_str,
                                                        const char *var_name,
                                                        slong prec,
                                                        complex_solver_roots_t *roots_out)
{
    fmpq_poly_t poly;
    acb_roots_t acb_roots;

    if (!poly_str || !var_name || !roots_out || prec <= 0) {
        return 0;
    }

    complex_solver_roots_clear(roots_out);
    complex_solver_roots_init(roots_out);

    if (!complex_solver_parse_univariate_fmpq_poly(poly_str, var_name, poly)) {
        return 0;
    }

    acb_roots_init(&acb_roots);
    fmpq_poly_acb_roots(&acb_roots, poly, prec);

    if (!complex_solver_roots_reserve(roots_out, acb_roots.num_roots)) {
        acb_roots_clear(&acb_roots);
        fmpq_poly_clear(poly);
        return 0;
    }

    for (slong i = 0; i < acb_roots.num_roots; i++) {
        acb_set(roots_out->roots[i], acb_roots.roots[i]);
        roots_out->multiplicities[i] = acb_roots.multiplicities ? acb_roots.multiplicities[i] : 1;
    }
    roots_out->num_roots = acb_roots.num_roots;

    if (!complex_solver_set_root_metadata(roots_out, poly_str, var_name, prec)) {
        acb_roots_clear(&acb_roots);
        fmpq_poly_clear(poly);
        return 0;
    }

    acb_roots_clear(&acb_roots);
    fmpq_poly_clear(poly);
    return 1;
}

static int complex_solver_split_polynomials(const char *poly_string,
                                            char **poly1_out,
                                            char **poly2_out)
{
    slong count = 0;
    char **parts = NULL;
    int ok = 0;

    if (poly1_out) *poly1_out = NULL;
    if (poly2_out) *poly2_out = NULL;
    if (!poly_string || !poly1_out || !poly2_out) {
        return 0;
    }

    parts = split_string(poly_string, &count);
    if (!parts || count != 2) {
        goto cleanup;
    }

    *poly1_out = strdup(parts[0]);
    *poly2_out = strdup(parts[1]);
    if (!*poly1_out || !*poly2_out) {
        free(*poly1_out);
        free(*poly2_out);
        *poly1_out = NULL;
        *poly2_out = NULL;
        goto cleanup;
    }

    ok = 1;

cleanup:
    if (parts) {
        free_split_strings(parts, count);
    }
    return ok;
}

static void complex_solver_skip_ws(const char **ptr)
{
    while (ptr && *ptr && isspace((unsigned char) **ptr)) {
        (*ptr)++;
    }
}

static int complex_solver_match_var_factor(const char **ptr,
                                           const char *var_name,
                                           slong *degree_out)
{
    const char *p;
    size_t var_len;
    slong degree = 1;

    if (!ptr || !*ptr || !var_name || !degree_out) {
        return 0;
    }

    p = *ptr;
    var_len = strlen(var_name);
    if (strncmp(p, var_name, var_len) != 0 || isalnum((unsigned char) p[var_len])) {
        return 0;
    }

    p += var_len;
    complex_solver_skip_ws(&p);
    if (*p == '^') {
        p++;
        complex_solver_skip_ws(&p);
        if (!isdigit((unsigned char) *p)) {
            return 0;
        }
        degree = strtol(p, (char **) &p, 10);
    }

    *degree_out = degree;
    *ptr = p;
    return 1;
}

static int complex_solver_match_numeric_factor(const char **ptr, fmpq_t value_out)
{
    const char *start;
    const char *p;
    size_t len;
    char *tmp;
    char *slash;

    if (!ptr || !*ptr || !value_out) {
        return 0;
    }

    p = *ptr;
    if (!(isdigit((unsigned char) *p) || *p == '.')) {
        return 0;
    }

    start = p;
    while (isdigit((unsigned char) *p) || *p == '.' || *p == '/') {
        p++;
    }

    len = (size_t) (p - start);
    tmp = (char *) malloc(len + 1);
    if (!tmp) {
        return 0;
    }
    memcpy(tmp, start, len);
    tmp[len] = '\0';

    slash = strchr(tmp, '/');
    if (slash) {
        fmpz_t num, den;
        *slash = '\0';
        fmpz_init(num);
        fmpz_init(den);
        fmpz_set_str(num, tmp, 10);
        fmpz_set_str(den, slash + 1, 10);
        fmpq_set_fmpz_frac(value_out, num, den);
        fmpz_clear(num);
        fmpz_clear(den);
    } else {
        fmpz_t num;
        fmpz_init(num);
        fmpz_set_str(num, tmp, 10);
        fmpq_set_fmpz(value_out, num);
        fmpz_clear(num);
    }

    free(tmp);
    *ptr = p;
    return 1;
}

static int complex_solver_iterate_bivariate_terms(const char *poly_str,
                                                  const char *x_name,
                                                  const char *y_name,
                                                  int (*cb)(void *userdata,
                                                            const fmpq_t coeff,
                                                            slong x_deg,
                                                            slong y_deg),
                                                  void *userdata)
{
    const char *ptr = poly_str;

    if (!poly_str || !x_name || !y_name || !cb) {
        return 0;
    }

    while (*ptr) {
        int sign = 1;
        int parsed_factor = 0;
        slong x_deg = 0;
        slong y_deg = 0;
        fmpq_t coeff;

        complex_solver_skip_ws(&ptr);
        while (*ptr == ')') {
            ptr++;
            complex_solver_skip_ws(&ptr);
        }
        if (*ptr == '\0') break;

        if (*ptr == '+') {
            ptr++;
        } else if (*ptr == '-') {
            sign = -1;
            ptr++;
        }
        complex_solver_skip_ws(&ptr);

        fmpq_init(coeff);
        fmpq_one(coeff);

        while (*ptr) {
            complex_solver_skip_ws(&ptr);

            if (*ptr == ')') {
                ptr++;
                continue;
            }

            if (*ptr == '(') {
                const char *factor_ptr = ptr + 1;
                int factor_sign = 1;
                slong deg = 0;
                fmpq_t factor_coeff;

                complex_solver_skip_ws(&factor_ptr);
                if (*factor_ptr == '+') {
                    factor_ptr++;
                    complex_solver_skip_ws(&factor_ptr);
                } else if (*factor_ptr == '-') {
                    factor_sign = -1;
                    factor_ptr++;
                    complex_solver_skip_ws(&factor_ptr);
                }

                if (complex_solver_match_var_factor(&factor_ptr, x_name, &deg)) {
                    x_deg += deg;
                    if (factor_sign < 0) fmpq_neg(coeff, coeff);
                    complex_solver_skip_ws(&factor_ptr);
                    if (*factor_ptr != ')') {
                        fmpq_clear(coeff);
                        return 0;
                    }
                    ptr = factor_ptr + 1;
                    parsed_factor = 1;
                } else if (complex_solver_match_var_factor(&factor_ptr, y_name, &deg)) {
                    y_deg += deg;
                    if (factor_sign < 0) fmpq_neg(coeff, coeff);
                    complex_solver_skip_ws(&factor_ptr);
                    if (*factor_ptr != ')') {
                        fmpq_clear(coeff);
                        return 0;
                    }
                    ptr = factor_ptr + 1;
                    parsed_factor = 1;
                } else {
                    fmpq_init(factor_coeff);
                    if (!complex_solver_match_numeric_factor(&factor_ptr, factor_coeff)) {
                        fmpq_clear(factor_coeff);
                        fmpq_clear(coeff);
                        return 0;
                    }
                    if (factor_sign < 0) fmpq_neg(factor_coeff, factor_coeff);
                    fmpq_mul(coeff, coeff, factor_coeff);
                    complex_solver_skip_ws(&factor_ptr);
                    fmpq_clear(factor_coeff);
                    if (*factor_ptr != ')') {
                        fmpq_clear(coeff);
                        return 0;
                    }
                    ptr = factor_ptr + 1;
                    parsed_factor = 1;
                }
            } else {
                slong deg = 0;
                fmpq_t factor_coeff;

                if (complex_solver_match_var_factor(&ptr, x_name, &deg)) {
                    x_deg += deg;
                    parsed_factor = 1;
                } else if (complex_solver_match_var_factor(&ptr, y_name, &deg)) {
                    y_deg += deg;
                    parsed_factor = 1;
                } else {
                    fmpq_init(factor_coeff);
                    if (!complex_solver_match_numeric_factor(&ptr, factor_coeff)) {
                        fmpq_clear(factor_coeff);
                        break;
                    }
                    fmpq_mul(coeff, coeff, factor_coeff);
                    fmpq_clear(factor_coeff);
                    parsed_factor = 1;
                }
            }

            complex_solver_skip_ws(&ptr);
            if (*ptr == '*') {
                ptr++;
                continue;
            }
            break;
        }

        if (!parsed_factor) {
            fmpq_clear(coeff);
            if (*ptr) {
                ptr++;
                continue;
            }
            break;
        }

        if (sign < 0) {
            fmpq_neg(coeff, coeff);
        }

        if (!cb(userdata, coeff, x_deg, y_deg)) {
            fmpq_clear(coeff);
            return 0;
        }
        fmpq_clear(coeff);
    }

    return 1;
}

typedef struct {
    acb_poly_struct *poly;
    acb_srcptr fixed_value;
    slong solve_var;
    slong prec;
} complex_solver_poly_builder_t;

static int complex_solver_build_poly_term_cb(void *userdata,
                                             const fmpq_t coeff,
                                             slong x_deg,
                                             slong y_deg)
{
    complex_solver_poly_builder_t *builder = (complex_solver_poly_builder_t *) userdata;
    slong target_deg = (builder->solve_var == 0) ? x_deg : y_deg;
    slong fixed_deg = (builder->solve_var == 0) ? y_deg : x_deg;
    acb_t term_coeff;
    acb_t existing;

    acb_init(term_coeff);
    acb_init(existing);
    acb_set_fmpq(term_coeff, coeff, builder->prec);
    if (fixed_deg > 0) {
        acb_t tmp;
        acb_init(tmp);
        acb_pow_ui(tmp, builder->fixed_value, (ulong) fixed_deg, builder->prec);
        acb_mul(term_coeff, term_coeff, tmp, builder->prec);
        acb_clear(tmp);
    }

    acb_poly_get_coeff_acb(existing, builder->poly, target_deg);
    acb_add(existing, existing, term_coeff, builder->prec);
    acb_poly_set_coeff_acb(builder->poly, target_deg, existing);

    acb_clear(term_coeff);
    acb_clear(existing);
    return 1;
}

static int complex_solver_build_acb_poly_from_bivariate_string(acb_poly_t out,
                                                               const char *poly_str,
                                                               const char *x_name,
                                                               const char *y_name,
                                                               slong solve_var,
                                                               const acb_t fixed_value,
                                                               slong prec)
{
    complex_solver_poly_builder_t builder;

    builder.poly = out;
    builder.fixed_value = fixed_value;
    builder.solve_var = solve_var;
    builder.prec = prec;

    acb_poly_zero(out);
    return complex_solver_iterate_bivariate_terms(poly_str, x_name, y_name,
                                                  complex_solver_build_poly_term_cb,
                                                  &builder);
}

static int complex_solver_pair_present(const complex_solver_solution_list_t *sols,
                                       const acb_t x,
                                       const acb_t y)
{
    if (!sols) {
        return 0;
    }

    for (slong i = 0; i < sols->num_solutions; i++) {
        if (acb_overlaps(sols->items[i].x, x) && acb_overlaps(sols->items[i].y, y)) {
            return 1;
        }
    }
    return 0;
}

static int complex_solver_solution_list_append_unique(complex_solver_solution_list_t *sols,
                                                      const acb_t x,
                                                      const acb_t y,
                                                      const acb_t residual1,
                                                      const acb_t residual2)
{
    if (!sols || !x || !y) {
        return 0;
    }
    if (complex_solver_pair_present(sols, x, y)) {
        return 1;
    }
    if (!complex_solver_solution_list_reserve(sols, sols->num_solutions + 1)) {
        return 0;
    }

    acb_set(sols->items[sols->num_solutions].x, x);
    acb_set(sols->items[sols->num_solutions].y, y);
    acb_set(sols->items[sols->num_solutions].residual1, residual1);
    acb_set(sols->items[sols->num_solutions].residual2, residual2);
    sols->num_solutions++;
    return 1;
}

typedef struct {
    acb_t acc;
    acb_srcptr x;
    acb_srcptr y;
    slong prec;
} complex_solver_eval_acc_t;

static int complex_solver_eval_term_cb(void *userdata,
                                       const fmpq_t coeff,
                                       slong x_deg,
                                       slong y_deg)
{
    complex_solver_eval_acc_t *acc = (complex_solver_eval_acc_t *) userdata;
    acb_t term;

    acb_init(term);
    acb_set_fmpq(term, coeff, acc->prec);
    if (x_deg > 0) {
        acb_t tmp;
        acb_init(tmp);
        acb_pow_ui(tmp, acc->x, (ulong) x_deg, acc->prec);
        acb_mul(term, term, tmp, acc->prec);
        acb_clear(tmp);
    }
    if (y_deg > 0) {
        acb_t tmp;
        acb_init(tmp);
        acb_pow_ui(tmp, acc->y, (ulong) y_deg, acc->prec);
        acb_mul(term, term, tmp, acc->prec);
        acb_clear(tmp);
    }
    acb_add(acc->acc, acc->acc, term, acc->prec);
    acb_clear(term);
    return 1;
}

static int complex_solver_evaluate_bivariate_string(acb_t out,
                                                    const char *poly_str,
                                                    const char *x_name,
                                                    const char *y_name,
                                                    const acb_t x,
                                                    const acb_t y,
                                                    slong prec)
{
    complex_solver_eval_acc_t acc;

    acb_init(acc.acc);
    acb_zero(acc.acc);
    acc.x = x;
    acc.y = y;
    acc.prec = prec;

    if (!complex_solver_iterate_bivariate_terms(poly_str, x_name, y_name,
                                                complex_solver_eval_term_cb, &acc)) {
        acb_clear(acc.acc);
        return 0;
    }

    acb_set(out, acc.acc);
    acb_clear(acc.acc);
    return 1;
}

static int complex_solver_verify_solution_pair_from_strings(const char *poly1_str,
                                                            const char *poly2_str,
                                                            const char *x_name,
                                                            const char *y_name,
                                                            const acb_t x,
                                                            const acb_t y,
                                                            acb_t residual1_out,
                                                            acb_t residual2_out,
                                                            slong prec)
{
    acb_t res1;
    acb_t res2;
    int ok;

    acb_init(res1);
    acb_init(res2);
    if (!complex_solver_evaluate_bivariate_string(res1, poly1_str, x_name, y_name, x, y, prec) ||
        !complex_solver_evaluate_bivariate_string(res2, poly2_str, x_name, y_name, x, y, prec)) {
        acb_clear(res1);
        acb_clear(res2);
        return 0;
    }
    if (residual1_out) {
        acb_set(residual1_out, res1);
    }
    if (residual2_out) {
        acb_set(residual2_out, res2);
    }
    ok = acb_contains_zero(res1) && acb_contains_zero(res2);
    acb_clear(res1);
    acb_clear(res2);
    return ok;
}

int complex_solver_solve_bivariate_2x2_from_string(const char *poly_string,
                                                   slong prec,
                                                   complex_solver_solution_list_t *sols_out)
{
    char *poly1_str = NULL;
    char *poly2_str = NULL;
    char **var_names = NULL;
    slong num_vars = 0;
    char *resultant = NULL;
    complex_solver_roots_t y_roots;
    int ok = 0;

    if (!poly_string || !sols_out || prec <= 0) {
        return 0;
    }

    complex_solver_solution_list_clear(sols_out);
    complex_solver_solution_list_init(sols_out);
    complex_solver_roots_init(&y_roots);

    if (!complex_solver_split_polynomials(poly_string, &poly1_str, &poly2_str)) {
        goto cleanup;
    }

    {
        char *poly_array[2];
        poly_array[0] = poly1_str;
        poly_array[1] = poly2_str;
        var_names = rational_extract_variables_improved(poly_array, 2, &num_vars);
    }
    if (!var_names || num_vars != 2) {
        goto cleanup;
    }

    resultant = dixon_str_rational(poly_string, var_names[0]);
    if (!resultant) {
        goto cleanup;
    }

    if (!complex_solver_univariate_complex_roots_from_string(resultant, var_names[1], prec, &y_roots)) {
        goto cleanup;
    }

    for (slong i = 0; i < y_roots.num_roots; i++) {
        acb_poly_t x_poly;
        acb_ptr x_root_vec = NULL;
        slong x_root_count = 0;

        acb_poly_init(x_poly);

        if (complex_solver_build_acb_poly_from_bivariate_string(x_poly, poly1_str,
                                                                var_names[0], var_names[1],
                                                                0, y_roots.roots[i], prec) &&
            acb_poly_degree(x_poly) > 0) {
            x_root_count = acb_poly_degree(x_poly);
            x_root_vec = _acb_vec_init(x_root_count);
            acb_poly_find_roots(x_root_vec, x_poly, NULL, 0, prec);

            for (slong j = 0; j < x_root_count; j++) {
                acb_t residual1;
                acb_t residual2;
                acb_init(residual1);
                acb_init(residual2);
                if (complex_solver_verify_solution_pair_from_strings(poly1_str, poly2_str,
                                                                     var_names[0], var_names[1],
                                                                     x_root_vec + j, y_roots.roots[i],
                                                                     residual1, residual2,
                                                                     prec)) {
                    complex_solver_solution_list_append_unique(sols_out,
                                                               x_root_vec + j,
                                                               y_roots.roots[i],
                                                               residual1,
                                                               residual2);
                }
                acb_clear(residual1);
                acb_clear(residual2);
            }
        }

        if (x_root_vec) {
            _acb_vec_clear(x_root_vec, x_root_count);
        }
        acb_poly_clear(x_poly);
    }

    sols_out->prec = prec;
    sols_out->source_system = strdup(poly_string);
    sols_out->x_name = strdup(var_names[0]);
    sols_out->y_name = strdup(var_names[1]);
    ok = (sols_out->source_system && sols_out->x_name && sols_out->y_name);
cleanup:
    complex_solver_roots_clear(&y_roots);
    free(resultant);
    free(poly1_str);
    free(poly2_str);
    if (var_names) {
        for (slong i = 0; i < num_vars; i++) {
            free(var_names[i]);
        }
        free(var_names);
    }
    if (!ok) {
        complex_solver_solution_list_clear(sols_out);
    }
    return ok;
}
