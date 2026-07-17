#include "fmpq_acb_roots.h"
#include "dixon_flint.h"
#include <flint/arith.h>
#include <flint/acb_poly.h>
#include <flint/arb.h>
#include <stdarg.h>
#include <time.h>

extern int g_dixon_verbose_level;
extern rational_root_scan_mode_t g_rational_root_scan_mode;

static int root_trace_enabled(void)
{
    return g_dixon_verbose_level >= 3;
}

static int root_info_enabled(void)
{
    return g_dixon_verbose_level >= 1;
}

static void root_trace_log(const char *fmt, ...)
{
    va_list args;

    if (!root_trace_enabled() || !fmt) return;

    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
}

void fmpq_roots_init(fmpq_roots_t *roots) {
    roots->roots = NULL;
    roots->multiplicities = NULL;
    roots->num_roots = 0;
    roots->alloc = 0;
}

void fmpq_roots_clear(fmpq_roots_t *roots) {
    if (roots->roots) {
        for (slong i = 0; i < roots->num_roots; i++) {
            fmpq_clear(roots->roots[i]);
        }
        free(roots->roots);
    }
    if (roots->multiplicities) {
        free(roots->multiplicities);
    }
    roots->roots = NULL;
    roots->multiplicities = NULL;
    roots->num_roots = 0;
    roots->alloc = 0;
}

static void fmpq_roots_add_root(fmpq_roots_t *roots, const fmpq_t root, slong mult) {
    if (roots->num_roots >= roots->alloc) {
        slong new_alloc = roots->alloc ? 2 * roots->alloc : 4;
        roots->roots = (fmpq_t *)realloc(roots->roots, new_alloc * sizeof(fmpq_t));
        roots->multiplicities = (slong *)realloc(roots->multiplicities, new_alloc * sizeof(slong));
        for (slong i = roots->alloc; i < new_alloc; i++) {
            fmpq_init(roots->roots[i]);
        }
        roots->alloc = new_alloc;
    }
    fmpq_set(roots->roots[roots->num_roots], root);
    roots->multiplicities[roots->num_roots] = mult;
    roots->num_roots++;
}

slong fmpq_poly_roots(fmpq_roots_t *roots, const fmpq_poly_t poly, int with_multiplicity) {
    slong degree = fmpq_poly_degree(poly);
    clock_t start_time = clock();
     
    if (degree <= 0) {
        return 0;
    }
    
    if (degree > FMPQ_ROOT_SEARCH_MAX_DEGREE) {
        return 0;
    }
    
    fmpz_poly_t int_poly, prim_poly, num_divs, den_divs;
    fmpz_t common_den, content, abs_const, abs_lead, coeff, gcd_nd;
    slong zero_mult = 0;
    
    fmpz_poly_init(int_poly);
    fmpz_poly_init(prim_poly);
    fmpz_poly_init(num_divs);
    fmpz_poly_init(den_divs);
    fmpz_init(common_den);
    fmpz_init(content);
    fmpz_init(abs_const);
    fmpz_init(abs_lead);
    fmpz_init(coeff);
    fmpz_init(gcd_nd);
    
    fmpq_poly_get_numerator(int_poly, poly);
    fmpq_poly_get_denominator(common_den, poly);
    
    fmpz_poly_content(content, int_poly);
    if (fmpz_is_zero(content)) {
        goto cleanup;
    }
    
    fmpz_poly_scalar_divexact_fmpz(prim_poly, int_poly, content);
    
    fmpz_poly_get_coeff_fmpz(coeff, prim_poly, fmpz_poly_degree(prim_poly));
    if (fmpz_sgn(coeff) < 0) {
        fmpz_poly_neg(prim_poly, prim_poly);
    }
    
    while (zero_mult <= fmpz_poly_degree(prim_poly)) {
        fmpz_poly_get_coeff_fmpz(coeff, prim_poly, zero_mult);
        if (!fmpz_is_zero(coeff)) break;
        zero_mult++;
    }
    
    if (zero_mult > 0) {
        fmpq_t zero_root;
        fmpq_init(zero_root);
        fmpq_zero(zero_root);
        fmpq_roots_add_root(roots, zero_root, with_multiplicity ? zero_mult : 1);
        fmpq_clear(zero_root);
    }
    
    if (zero_mult <= fmpz_poly_degree(prim_poly)) {
        fmpq_t candidate, value;
        
        fmpq_init(candidate);
        fmpq_init(value);
        
        fmpz_poly_get_coeff_fmpz(abs_const, prim_poly, zero_mult);
        fmpz_abs(abs_const, abs_const);
        fmpz_poly_get_coeff_fmpz(abs_lead, prim_poly, fmpz_poly_degree(prim_poly));
        fmpz_abs(abs_lead, abs_lead);
        
        arith_divisors(num_divs, abs_const);
        arith_divisors(den_divs, abs_lead);
         
        slong candidate_count = 2 * fmpz_poly_length(num_divs) * fmpz_poly_length(den_divs);
        root_trace_log("[root-debug:v3] Rational-root test: degree=%ld, numerator divisors=%ld, denominator divisors=%ld, candidates=%ld.\n",
                       degree, fmpz_poly_length(num_divs), fmpz_poly_length(den_divs), candidate_count);
        if (candidate_count > FMPQ_ROOT_SEARCH_HARD_MAX_CANDIDATES) {
            root_trace_log("[root-debug:v3] Skipping rational-root scan: %ld candidates exceed hard cap %ld.\n",
                           candidate_count, (slong) FMPQ_ROOT_SEARCH_HARD_MAX_CANDIDATES);
            fmpq_clear(candidate);
            fmpq_clear(value);
            goto cleanup;
        }
        if (g_rational_root_scan_mode == RATIONAL_ROOT_SCAN_OFF) {
            root_trace_log("[root-debug:v3] Rational-root scan disabled by CLI option.\n");
            fmpq_clear(candidate);
            fmpq_clear(value);
            goto cleanup;
        }
        if (g_rational_root_scan_mode == RATIONAL_ROOT_SCAN_AUTO &&
            candidate_count > FMPQ_ROOT_SEARCH_AUTO_MAX_CANDIDATES) {
            if (root_info_enabled()) {
                printf("Skipping rational-root scan: %ld candidate(s) exceed auto threshold %ld.\n",
                       candidate_count, (slong) FMPQ_ROOT_SEARCH_AUTO_MAX_CANDIDATES);
            }
            root_trace_log("[root-debug:v3] Skipping exhaustive rational-root scan: %ld candidates exceed auto threshold %ld.\n",
                           candidate_count, (slong) FMPQ_ROOT_SEARCH_AUTO_MAX_CANDIDATES);
            fmpq_clear(candidate);
            fmpq_clear(value);
            goto cleanup;
        }
        if (g_rational_root_scan_mode == RATIONAL_ROOT_SCAN_FORCE &&
            candidate_count > FMPQ_ROOT_SEARCH_AUTO_MAX_CANDIDATES) {
            root_trace_log("[root-debug:v3] Forcing exhaustive rational-root scan despite %ld candidates (auto threshold %ld).\n",
                           candidate_count, (slong) FMPQ_ROOT_SEARCH_AUTO_MAX_CANDIDATES);
        }
         
        for (slong i = 0; i < fmpz_poly_length(num_divs); i++) {
            fmpz_t num;
            fmpz_init(num);
            fmpz_poly_get_coeff_fmpz(num, num_divs, i);
            
            for (slong j = 0; j < fmpz_poly_length(den_divs); j++) {
                fmpz_t den;
                fmpz_init(den);
                fmpz_poly_get_coeff_fmpz(den, den_divs, j);
                
                fmpz_gcd(gcd_nd, num, den);
                if (fmpz_is_one(gcd_nd)) {
                    for (int sign = -1; sign <= 1; sign += 2) {
                        if (sign < 0) fmpz_neg(num, num);
                        fmpq_set_fmpz_frac(candidate, num, den);
                        fmpq_canonicalise(candidate);
                        fmpq_poly_evaluate_fmpq(value, poly, candidate);
                        
                        if (fmpq_is_zero(value)) {
                            int already_found = 0;
                            for (slong k = 0; k < roots->num_roots; k++) {
                                if (fmpq_equal(roots->roots[k], candidate)) {
                                    already_found = 1;
                                    break;
                                }
                            }
                            if (!already_found) {
                                fmpq_roots_add_root(roots, candidate, 1);
                            }
                        }
                        if (sign < 0) fmpz_neg(num, num);
                    }
                }
                
                fmpz_clear(den);
            }
            
            fmpz_clear(num);
        }
        
        fmpq_clear(candidate);
        fmpq_clear(value);
    }
    
cleanup:
    fmpz_poly_clear(int_poly);
    fmpz_poly_clear(prim_poly);
    fmpz_poly_clear(num_divs);
    fmpz_poly_clear(den_divs);
    fmpz_clear(common_den);
    fmpz_clear(content);
    fmpz_clear(abs_const);
    fmpz_clear(abs_lead);
    fmpz_clear(coeff);
    fmpz_clear(gcd_nd);

    root_trace_log("[root-debug:v3] Rational-root test finished in %.3f s, found %ld root(s).\n",
                   (double) (clock() - start_time) / CLOCKS_PER_SEC,
                   roots->num_roots);
     
    return roots->num_roots;
}

void fmpq_roots_print(const fmpq_roots_t *roots) {
    printf("Found %ld roots:\n", roots->num_roots);
    for (slong i = 0; i < roots->num_roots; i++) {
        printf("  Root %ld: ", i + 1);
        fmpq_print(roots->roots[i]);
        printf(" (Multiplicity: %ld)\n", roots->multiplicities[i]);
    }
}

char* fmpq_roots_to_string(const fmpq_roots_t *roots) {
    if (roots->num_roots == 0) {
        return strdup("No rational roots found");
    }
    
    size_t total_len = 1;
    for (slong i = 0; i < roots->num_roots; i++) {
        char *root_str = fmpq_get_str(NULL, 10, roots->roots[i]);
        total_len += strlen(root_str) + 10;
        flint_free(root_str);
    }
    
    char *result = (char *)malloc(total_len);
    result[0] = '\0';
    
    for (slong i = 0; i < roots->num_roots; i++) {
        char *root_str = fmpq_get_str(NULL, 10, roots->roots[i]);
        if (i > 0) {
            strcat(result, ", ");
        }
        strcat(result, root_str);
        flint_free(root_str);
    }
    
    return result;
}

void acb_roots_init(acb_roots_t *roots) {
    roots->roots = NULL;
    roots->multiplicities = NULL;
    roots->num_roots = 0;
    roots->alloc = 0;
}

void acb_roots_clear(acb_roots_t *roots) {
    if (roots->roots) {
        for (slong i = 0; i < roots->num_roots; i++) {
            acb_clear(roots->roots[i]);
        }
        free(roots->roots);
    }
    if (roots->multiplicities) {
        free(roots->multiplicities);
    }
    roots->roots = NULL;
    roots->multiplicities = NULL;
    roots->num_roots = 0;
    roots->alloc = 0;
}

static void acb_roots_add_root(acb_roots_t *roots, const acb_t root, slong mult) {
    if (roots->num_roots >= roots->alloc) {
        slong new_alloc = roots->alloc ? 2 * roots->alloc : 4;
        roots->roots = (acb_t *)realloc(roots->roots, new_alloc * sizeof(acb_t));
        roots->multiplicities = (slong *)realloc(roots->multiplicities, new_alloc * sizeof(slong));
        for (slong i = roots->alloc; i < new_alloc; i++) {
            acb_init(roots->roots[i]);
        }
        roots->alloc = new_alloc;
    }
    acb_set(roots->roots[roots->num_roots], root);
    roots->multiplicities[roots->num_roots] = mult;
    roots->num_roots++;
}

slong acb_poly_roots(acb_roots_t *roots, const acb_poly_t poly, slong prec) {
    slong degree = acb_poly_degree(poly);
    clock_t start_time = clock();
     
    if (degree <= 0) {
        return 0;
    }
    
    acb_ptr roots_array = _acb_vec_init(degree);
    
    acb_poly_find_roots(roots_array, poly, NULL, 0, prec);
    
    for (slong i = 0; i < degree; i++) {
        acb_roots_add_root(roots, roots_array + i, 1);
    }
    
    _acb_vec_clear(roots_array, degree);

    root_trace_log("[root-debug:v3] acb_poly_find_roots finished in %.3f s for degree %ld.\n",
                   (double) (clock() - start_time) / CLOCKS_PER_SEC,
                   degree);
     
    return roots->num_roots;
}

slong fmpq_poly_acb_roots(acb_roots_t *roots, const fmpq_poly_t poly, slong prec) {
    slong degree = fmpq_poly_degree(poly);
    clock_t convert_start = clock();
    clock_t solve_start;
     
    if (degree <= 0) {
        return 0;
    }
    
    acb_poly_t acb_poly;
    acb_poly_init(acb_poly);
    acb_poly_set_fmpq_poly(acb_poly, poly, prec);
    root_trace_log("[root-debug:v3] Converted fmpq polynomial to acb polynomial in %.3f s (degree %ld, prec=%ld).\n",
                   (double) (clock() - convert_start) / CLOCKS_PER_SEC,
                   degree, prec);
     
    solve_start = clock();
    slong result = acb_poly_roots(roots, acb_poly, prec);
    root_trace_log("[root-debug:v3] Approximate complex-root phase total %.3f s.\n",
                   (double) (clock() - solve_start) / CLOCKS_PER_SEC);
     
    acb_poly_clear(acb_poly);
     
    return result;
}

static void approx_root_report_one(FILE *fp, slong index, const acb_t root,
                                   const acb_t residual, int complex_root,
                                   slong digits)
{
    if (!fp) return;

    fprintf(fp, "  Root %ld interval: ", index);
    if (complex_root) {
        acb_fprintd(fp, root, digits);
    } else {
        arb_fprintd(fp, acb_realref(root), digits);
    }
    fprintf(fp, "\n    Residual interval: ");
    acb_fprintd(fp, residual, digits);
    fprintf(fp, "\n");
}

static void approx_real_root_report_one(FILE *fp, slong index, const arb_t root,
                                        const arb_t normalized_residual,
                                        slong multiplicity,
                                        slong digits)
{
    if (!fp) return;

    fprintf(fp, "  Root %ld interval: ", index);
    arb_fprintd(fp, root, digits);
    fprintf(fp, " (Multiplicity: %ld)", multiplicity);
    fprintf(fp, "\n    Normalized residual interval: ");
    arb_fprintd(fp, normalized_residual, digits);
    fprintf(fp, "\n");
}

static void fmpq_poly_get_primitive_integer_poly(fmpz_poly_t integer_poly,
                                                 const fmpq_poly_t rational_poly)
{
    fmpz_poly_t numerator;
    fmpz_poly_init(numerator);
    fmpq_poly_get_numerator(numerator, rational_poly);
    fmpz_poly_primitive_part(integer_poly, numerator);
    fmpz_poly_clear(numerator);
}

static void arb_roots_add_root(arb_roots_t *roots, const arb_t root, slong mult);

static slong fmpz_poly_isolate_real_roots_with_multiplicity(
    arb_roots_t *roots, const fmpz_poly_t integer_poly, slong prec)
{
    fmpz_poly_factor_t factors;

    fmpz_poly_factor_init(factors);
    fmpz_poly_factor_squarefree(factors, integer_poly);

    for (slong factor_idx = 0; factor_idx < factors->num; factor_idx++) {
        slong degree = fmpz_poly_degree(factors->p + factor_idx);
        arb_ptr factor_roots;
        slong count;

        if (degree <= 0) continue;
        factor_roots = _arb_vec_init(degree);
        count = arb_fmpz_poly_real_roots(factor_roots, factors->p + factor_idx,
                                         0, prec);
        for (slong i = 0; i < count; i++) {
            if (arb_is_finite(factor_roots + i)) {
                arb_roots_add_root(roots, factor_roots + i,
                                   factors->exp[factor_idx]);
            }
        }
        _arb_vec_clear(factor_roots, degree);
    }

    fmpz_poly_factor_clear(factors);
    return roots->num_roots;
}

static void fmpz_poly_normalized_residual(arb_t normalized,
                                          const fmpz_poly_t poly,
                                          const arb_t root, slong prec)
{
    arb_t residual, abs_root, scale_x, power, denominator, term, one;
    fmpz_t coeff;
    slong length = fmpz_poly_length(poly);

    arb_init(residual);
    arb_init(abs_root);
    arb_init(scale_x);
    arb_init(power);
    arb_init(denominator);
    arb_init(term);
    arb_init(one);
    fmpz_init(coeff);

    arb_fmpz_poly_evaluate_arb(residual, poly, root, prec);
    arb_abs(abs_root, root);
    arb_one(one);
    arb_max(scale_x, abs_root, one, prec);
    arb_one(power);
    arb_zero(denominator);

    for (slong i = 0; i < length; i++) {
        fmpz_poly_get_coeff_fmpz(coeff, poly, i);
        fmpz_abs(coeff, coeff);
        arb_set_fmpz(term, coeff);
        arb_mul(term, term, power, prec);
        arb_add(denominator, denominator, term, prec);
        arb_mul(power, power, scale_x, prec);
    }

    if (arb_contains_zero(denominator)) {
        arb_indeterminate(normalized);
    } else {
        arb_div(normalized, residual, denominator, prec);
    }

    fmpz_clear(coeff);
    arb_clear(residual);
    arb_clear(abs_root);
    arb_clear(scale_x);
    arb_clear(power);
    arb_clear(denominator);
    arb_clear(term);
    arb_clear(one);
}

slong fmpq_poly_approx_roots_report(const fmpq_poly_t poly, slong prec,
                                    int complex_roots, FILE *fp_file,
                                    int print_to_stdout)
{
    slong degree = fmpq_poly_degree(poly);
    slong reported = 0;
    slong digits = FLINT_MAX(15, (slong) (0.30103 * (double) prec));

    if (degree <= 0 || prec < 2) return 0;

    if (!complex_roots) {
        fmpz_poly_t integer_poly;
        arb_roots_t real_roots;

        fmpz_poly_init(integer_poly);
        fmpq_poly_get_primitive_integer_poly(integer_poly, poly);
        arb_roots_init(&real_roots);
        fmpz_poly_isolate_real_roots_with_multiplicity(&real_roots, integer_poly, prec);

        if (print_to_stdout) {
            printf("Approximate coefficient real-root analysis (%ld-bit Arb precision):\n", prec);
        }
        if (fp_file) {
            fprintf(fp_file, "Approximate coefficient real-root analysis (%ld-bit Arb precision):\n", prec);
        }

        for (slong i = 0; i < real_roots.num_roots; i++) {
            arb_t normalized_residual;
            arb_init(normalized_residual);
            fmpz_poly_normalized_residual(normalized_residual, integer_poly,
                                          real_roots.roots[i], prec);
            reported++;
            if (print_to_stdout) {
                approx_real_root_report_one(stdout, reported, real_roots.roots[i], normalized_residual,
                                            real_roots.multiplicities[i], digits);
            }
            if (fp_file) {
                approx_real_root_report_one(fp_file, reported, real_roots.roots[i], normalized_residual,
                                            real_roots.multiplicities[i], digits);
            }
            arb_clear(normalized_residual);
        }

        if (reported == 0) {
            if (print_to_stdout) printf("No isolated approximate real roots found.\n");
            if (fp_file) fprintf(fp_file, "No isolated approximate real roots found.\n");
        }

        arb_roots_clear(&real_roots);
        fmpz_poly_clear(integer_poly);
        return reported;
    }

    acb_poly_t approx_poly;
    acb_ptr roots_array;

    acb_poly_init(approx_poly);
    acb_poly_set_fmpq_poly(approx_poly, poly, prec);
    roots_array = _acb_vec_init(degree);
    acb_poly_find_roots(roots_array, approx_poly, NULL, 0, prec);

    if (print_to_stdout) {
        printf("Approximate coefficient root analysis (%ld-bit Arb precision):\n", prec);
    }
    if (fp_file) {
        fprintf(fp_file, "Approximate coefficient root analysis (%ld-bit Arb precision):\n", prec);
    }

    for (slong i = 0; i < degree; i++) {
        acb_t residual;

        acb_init(residual);
        acb_poly_evaluate(residual, approx_poly, roots_array + i, prec);
        reported++;
        if (print_to_stdout) {
            approx_root_report_one(stdout, reported, roots_array + i, residual,
                                   complex_roots, digits);
        }
        if (fp_file) {
            approx_root_report_one(fp_file, reported, roots_array + i, residual,
                                   complex_roots, digits);
        }
        acb_clear(residual);
    }

    if (reported == 0) {
        if (print_to_stdout) printf("No approximate %s roots found.\n", complex_roots ? "complex" : "real");
        if (fp_file) fprintf(fp_file, "No approximate %s roots found.\n", complex_roots ? "complex" : "real");
    }

    _acb_vec_clear(roots_array, degree);
    acb_poly_clear(approx_poly);
    return reported;
}

void acb_roots_print(const acb_roots_t *roots) {
    printf("Found %ld approximate roots:\n", roots->num_roots);
    for (slong i = 0; i < roots->num_roots; i++) {
        printf("  Root %ld: ", i + 1);
        acb_printd(roots->roots[i], 15);
        printf(" (Multiplicity: %ld)\n", roots->multiplicities[i]);
    }
}

void fmpq_acb_roots_init(fmpq_acb_roots_t *roots) {
    fmpq_roots_init(&roots->rational_roots);
    acb_roots_init(&roots->approximate_roots);
    arb_roots_init(&roots->real_roots);
}

void fmpq_acb_roots_clear(fmpq_acb_roots_t *roots) {
    fmpq_roots_clear(&roots->rational_roots);
    acb_roots_clear(&roots->approximate_roots);
    arb_roots_clear(&roots->real_roots);
}

slong fmpq_poly_all_roots(fmpq_acb_roots_t *roots, const fmpq_poly_t poly, slong prec) {
    slong degree = fmpq_poly_degree(poly);
    
    if (degree <= 0) {
        return 0;
    }
    
    fmpq_poly_roots(&roots->rational_roots, poly, 0);
    
    fmpq_poly_acb_roots(&roots->approximate_roots, poly, prec);
    
    return roots->rational_roots.num_roots + roots->approximate_roots.num_roots;
}

void fmpq_acb_roots_print(const fmpq_acb_roots_t *roots) {
    if (roots->rational_roots.num_roots > 0) {
        printf("\nRational roots:\n");
        fmpq_roots_print(&roots->rational_roots);
    }
    
    if (roots->approximate_roots.num_roots > 0) {
        printf("\nApproximate complex roots:\n");
        acb_roots_print(&roots->approximate_roots);
    }
    
    if (roots->rational_roots.num_roots == 0 && roots->approximate_roots.num_roots == 0) {
        printf("No roots found.\n");
    }
}

void arb_roots_init(arb_roots_t *roots) {
    roots->roots = NULL;
    roots->multiplicities = NULL;
    roots->num_roots = 0;
    roots->alloc = 0;
}

void arb_roots_clear(arb_roots_t *roots) {
    if (roots->roots) {
        for (slong i = 0; i < roots->num_roots; i++) {
            arb_clear(roots->roots[i]);
        }
        free(roots->roots);
    }
    if (roots->multiplicities) {
        free(roots->multiplicities);
    }
    roots->roots = NULL;
    roots->multiplicities = NULL;
    roots->num_roots = 0;
    roots->alloc = 0;
}

static void arb_roots_add_root(arb_roots_t *roots, const arb_t root, slong mult) {
    if (roots->num_roots >= roots->alloc) {
        slong new_alloc = roots->alloc ? 2 * roots->alloc : 4;
        roots->roots = (arb_t *)realloc(roots->roots, new_alloc * sizeof(arb_t));
        roots->multiplicities = (slong *)realloc(roots->multiplicities, new_alloc * sizeof(slong));
        for (slong i = roots->alloc; i < new_alloc; i++) {
            arb_init(roots->roots[i]);
        }
        roots->alloc = new_alloc;
    }
    arb_set(roots->roots[roots->num_roots], root);
    roots->multiplicities[roots->num_roots] = mult;
    roots->num_roots++;
}

slong acb_roots_to_real(arb_roots_t *real_roots, const acb_roots_t *complex_roots, slong prec) {
    clock_t start_time = clock();
    (void) prec;
    for (slong i = 0; i < complex_roots->num_roots; i++) {
        arb_t real_part, imag_part;
        arb_init(real_part);
        arb_init(imag_part);
        
        acb_get_real(real_part, complex_roots->roots[i]);
        acb_get_imag(imag_part, complex_roots->roots[i]);
        
        int is_real = 0;
        if (arb_contains_zero(imag_part)) {
            is_real = 1;
        } else {
            mag_t imag_mag;
            mag_init(imag_mag);
            arb_get_mag(imag_mag, imag_part);
            
            mag_t threshold;
            mag_init(threshold);
            mag_set_ui_2exp_si(threshold, 1, -30);
            
            if (mag_cmp(imag_mag, threshold) <= 0) {
                is_real = 1;
            }
            
            mag_clear(imag_mag);
            mag_clear(threshold);
        }
        
        if (is_real) {
            int already_found = 0;
            for (slong j = 0; j < real_roots->num_roots; j++) {
                if (arb_overlaps(real_roots->roots[j], real_part)) {
                    already_found = 1;
                    break;
                }
            }
            
            if (!already_found) {
                arb_roots_add_root(real_roots, real_part, complex_roots->multiplicities[i]);
            }
        }
        
        arb_clear(real_part);
        arb_clear(imag_part);
    }

    root_trace_log("[root-debug:v3] Filtered %ld complex root(s) down to %ld real root interval(s) in %.3f s.\n",
                   complex_roots->num_roots, real_roots->num_roots,
                   (double) (clock() - start_time) / CLOCKS_PER_SEC);
     
    return real_roots->num_roots;
}

void arb_roots_print(const arb_roots_t *roots) {
    printf("Found %ld approximate real roots:\n", roots->num_roots);
    for (slong i = 0; i < roots->num_roots; i++) {
        printf("  Root %ld: ", i + 1);
        arb_printd(roots->roots[i], 15);
        printf(" (Multiplicity: %ld)\n", roots->multiplicities[i]);
    }
}

slong fmpq_poly_real_roots(fmpq_acb_roots_t *roots, const fmpq_poly_t poly, slong prec) {
    slong degree = fmpq_poly_degree(poly);
    clock_t start_time = clock();
     
    if (degree <= 0) {
        return 0;
    }
     
    fmpz_poly_t integer_poly;
    slong num_real_roots;

    fmpz_poly_init(integer_poly);
    fmpq_poly_get_primitive_integer_poly(integer_poly, poly);
    num_real_roots = fmpz_poly_isolate_real_roots_with_multiplicity(
        &roots->real_roots, integer_poly, prec);
    fmpz_poly_clear(integer_poly);

    root_trace_log("[root-debug:v3] Dedicated real-root pipeline summary: degree=%ld, isolated=%ld, finite=%ld, total=%.3f s.\n",
                   degree,
                   num_real_roots,
                   roots->real_roots.num_roots,
                   (double) (clock() - start_time) / CLOCKS_PER_SEC);
     
    return roots->real_roots.num_roots;
}

void fmpq_acb_roots_print_real(const fmpq_acb_roots_t *roots) {
    if (roots->rational_roots.num_roots > 0) {
        printf("\nRational roots:\n");
        fmpq_roots_print(&roots->rational_roots);
    }
    
    if (roots->real_roots.num_roots > 0) {
        printf("\nApproximate real roots:\n");
        arb_roots_print(&roots->real_roots);
    }
    
    if (roots->rational_roots.num_roots == 0 && roots->real_roots.num_roots == 0) {
        printf("No roots found.\n");
    }
}

static void fmpq_poly_real_roots_report_to_file(FILE *fp,
                                                const fmpz_poly_t integer_poly,
                                                const arb_roots_t *roots,
                                                slong prec, slong digits)
{
    if (!fp) return;

    fprintf(fp, "Found %ld certified real root interval(s):\n", roots->num_roots);
    for (slong i = 0; i < roots->num_roots; i++) {
        arb_t normalized_residual;
        arb_init(normalized_residual);
        fmpz_poly_normalized_residual(normalized_residual, integer_poly,
                                      roots->roots[i], prec);
        fprintf(fp, "  Root %ld: ", i + 1);
        arb_fprintd(fp, roots->roots[i], digits);
        fprintf(fp, " (Multiplicity: %ld)\n", roots->multiplicities[i]);
        fprintf(fp, "    Normalized residual interval: ");
        arb_fprintd(fp, normalized_residual, digits);
        fprintf(fp, "\n");
        arb_clear(normalized_residual);
    }
}

void fmpq_poly_real_roots_report(const fmpq_poly_t poly,
                                 const arb_roots_t *roots, slong prec,
                                 FILE *fp_file, int print_to_stdout)
{
    fmpz_poly_t integer_poly;
    slong digits = FLINT_MAX(15, (slong) (0.30103 * (double) prec));

    fmpz_poly_init(integer_poly);
    fmpq_poly_get_primitive_integer_poly(integer_poly, poly);
    if (print_to_stdout) {
        fmpq_poly_real_roots_report_to_file(stdout, integer_poly, roots, prec, digits);
    }
    if (fp_file) {
        fmpq_poly_real_roots_report_to_file(fp_file, integer_poly, roots, prec, digits);
    }
    fmpz_poly_clear(integer_poly);
}

void fmpq_roots_print_to_file(FILE *fp, const fmpq_roots_t *roots) {
    fprintf(fp, "Found %ld rational roots:\n", roots->num_roots);
    for (slong i = 0; i < roots->num_roots; i++) {
        fprintf(fp, "  Root %ld: ", i + 1);
        fmpq_fprint(fp, roots->roots[i]);
        fprintf(fp, " (Multiplicity: %ld)\n", roots->multiplicities[i]);
    }
}

void acb_roots_print_to_file(FILE *fp, const acb_roots_t *roots) {
    fprintf(fp, "Found %ld approximate complex roots:\n", roots->num_roots);
    for (slong i = 0; i < roots->num_roots; i++) {
        fprintf(fp, "  Root %ld: ", i + 1);
        acb_fprintd(fp, roots->roots[i], 15);
        fprintf(fp, " (Multiplicity: %ld)\n", roots->multiplicities[i]);
    }
}

void arb_roots_print_to_file(FILE *fp, const arb_roots_t *roots) {
    fprintf(fp, "Found %ld approximate real roots:\n", roots->num_roots);
    for (slong i = 0; i < roots->num_roots; i++) {
        fprintf(fp, "  Root %ld: ", i + 1);
        arb_fprintd(fp, roots->roots[i], 15);
        fprintf(fp, " (Multiplicity: %ld)\n", roots->multiplicities[i]);
    }
}

void fmpq_acb_roots_print_all_to_file(FILE *fp, const fmpq_acb_roots_t *roots) {
    if (roots->rational_roots.num_roots > 0) {
        fprintf(fp, "\nRational roots:\n");
        fmpq_roots_print_to_file(fp, &roots->rational_roots);
    }
    
    if (roots->approximate_roots.num_roots > 0) {
        fprintf(fp, "\nApproximate complex roots:\n");
        acb_roots_print_to_file(fp, &roots->approximate_roots);
    }
    
    if (roots->rational_roots.num_roots == 0 && roots->approximate_roots.num_roots == 0) {
        fprintf(fp, "No roots found.\n");
    }
}
