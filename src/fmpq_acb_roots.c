#include "fmpq_acb_roots.h"
#include <flint/arith.h>
#include <flint/acb_poly.h>
#include <flint/arb.h>

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
        if (candidate_count > FMPQ_ROOT_SEARCH_MAX_CANDIDATES) {
            fmpq_clear(candidate);
            fmpq_clear(value);
            goto cleanup;
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
    
    if (degree <= 0) {
        return 0;
    }
    
    acb_ptr roots_array = _acb_vec_init(degree);
    
    acb_poly_find_roots(roots_array, poly, NULL, 0, prec);
    
    for (slong i = 0; i < degree; i++) {
        acb_roots_add_root(roots, roots_array + i, 1);
    }
    
    _acb_vec_clear(roots_array, degree);
    
    return roots->num_roots;
}

slong fmpq_poly_acb_roots(acb_roots_t *roots, const fmpq_poly_t poly, slong prec) {
    slong degree = fmpq_poly_degree(poly);
    
    if (degree <= 0) {
        return 0;
    }
    
    acb_poly_t acb_poly;
    acb_poly_init(acb_poly);
    acb_poly_set_fmpq_poly(acb_poly, poly, prec);
    
    slong result = acb_poly_roots(roots, acb_poly, prec);
    
    acb_poly_clear(acb_poly);
    
    return result;
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
    
    if (degree <= 0) {
        return 0;
    }
    
    fmpq_poly_roots(&roots->rational_roots, poly, 0);
    
    fmpq_poly_acb_roots(&roots->approximate_roots, poly, prec);
    
    acb_roots_to_real(&roots->real_roots, &roots->approximate_roots, prec);
    
    return roots->rational_roots.num_roots + roots->real_roots.num_roots;
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
