#include "fmpq_acb_roots.h"
#include <flint/arith.h>

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
