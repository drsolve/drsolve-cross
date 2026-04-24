#include "fq_mvpoly.h"

/* ============================================================================
 * Field Equation Reduction Mode
 * ============================================================================ */

/* When enabled, polynomial multiplication reduces each variable modulo x^q - x,
 * i.e. any exponent >= q is replaced using x^q = x (Frobenius: a^q = a in F_q).
 * Enable via fq_mvpoly_set_field_equation_reduction(1) before computation. */
int g_field_equation_reduction = 0;

void fq_mvpoly_set_field_equation_reduction(int enable)
{
    g_field_equation_reduction = enable;
}

/* Reduce a single exponent e modulo x^q - x.
 * The relation x^q = x means x^(q+k) = x^(1+k), so for e >= q we reduce:
 *   e -> ((e - 1) % (q - 1)) + 1   when e >= 1
 * For e == 0 the term is a constant and needs no reduction. */
static slong reduce_exp_field(slong e, slong q)
{
    if (e == 0 || e < q) return e;
    /* x^q = x^1, so the exponent cycles with period q-1 starting at 1 */
    return ((e - 1) % (q - 1)) + 1;
}

/* Reduce all variable (and parameter) exponents of poly in-place using x^q = x.
 * Terms that become identical after reduction are combined. */
void fq_mvpoly_reduce_field_equation(fq_mvpoly_t *poly)
{
    if (!poly || poly->nterms == 0) return;

    /* Compute q = p^d from the field context */
    mp_limb_t p = fq_nmod_ctx_prime(poly->ctx);
    slong     d = fq_nmod_ctx_degree(poly->ctx);
    slong     q = 1;
    for (slong i = 0; i < d; i++) {
        if (q > WORD_MAX / (slong)p) {
            q = WORD_MAX;
            break;
        }
        q *= (slong)p;
    }

    if (q <= 1) return; /* nothing to reduce */

    /* Apply the exponent reduction to every term in-place */
    for (slong i = 0; i < poly->nterms; i++) {
        if (poly->terms[i].var_exp) {
            for (slong k = 0; k < poly->nvars; k++)
                poly->terms[i].var_exp[k] =
                    reduce_exp_field(poly->terms[i].var_exp[k], q);
        }
        if (poly->terms[i].par_exp) {
            for (slong k = 0; k < poly->npars; k++)
                poly->terms[i].par_exp[k] =
                    reduce_exp_field(poly->terms[i].par_exp[k], q);
        }
    }

    /* Re-combine like terms that may now share the same monomial.
     * Build a fresh polynomial by re-inserting every term through add_term
     * (which handles duplicate detection and coefficient accumulation). */
    fq_mvpoly_t combined;
    fq_mvpoly_init(&combined, poly->nvars, poly->npars, poly->ctx);
    for (slong i = 0; i < poly->nterms; i++) {
        fq_mvpoly_add_term(&combined,
                           poly->terms[i].var_exp,
                           poly->terms[i].par_exp,
                           poly->terms[i].coeff);
    }
    fq_mvpoly_clear(poly);
    *poly = combined;
}

/* ============================================================================
 * Basic Operations Implementation
 * ============================================================================ */

void fq_mvpoly_init(fq_mvpoly_t *p, slong nvars, slong npars, const fq_nmod_ctx_t ctx) {
    p->nvars = nvars;
    p->npars = npars;
    p->nterms = 0;
    p->alloc = 16;
    p->terms = (fq_monomial_t*) flint_malloc(p->alloc * sizeof(fq_monomial_t));
    
    memset(p->terms, 0, p->alloc * sizeof(fq_monomial_t));
    
    p->ctx = ctx; 
}

void fq_mvpoly_clear(fq_mvpoly_t *p) {
    if (p->terms) {
        for (slong i = 0; i < p->nterms; i++) {
            fq_nmod_clear(p->terms[i].coeff, p->ctx);
            if (p->terms[i].var_exp) flint_free(p->terms[i].var_exp);
            if (p->terms[i].par_exp) flint_free(p->terms[i].par_exp);
        }

        flint_free(p->terms);
        p->terms = NULL;
    }
    p->nterms = 0;
    p->alloc = 0;
    p->nvars = 0;
    p->npars = 0;
}

void fq_mvpoly_clear_safe(fq_mvpoly_t *p) {
    /* Only clear if structure appears to be properly initialized */
    if (p->terms && p->nterms >= 0 && p->ctx) {
        fq_mvpoly_clear(p);
    } else {
        /* Reset to safe state */
        memset(p, 0, sizeof(fq_mvpoly_t));
    }
}

/* ============================================================================
 * Comparison Function Implementation
 * ============================================================================ */

int compare_fq_degrees(const void *a, const void *b) {
    fq_index_degree_pair *pa = (fq_index_degree_pair *)a;
    fq_index_degree_pair *pb = (fq_index_degree_pair *)b;
    
    /* Sort by degree (ascending) */
    if (pa->degree < pb->degree) return -1;
    if (pa->degree > pb->degree) return 1;
    
    /* Same degree, sort by index for stability */
    if (pa->index < pb->index) return -1;
    if (pa->index > pb->index) return 1;
    return 0;
}

/* ============================================================================
 * Term Management Implementation
 * ============================================================================ */

void fq_mvpoly_add_term_fast(fq_mvpoly_t *p, const slong *var_exp, const slong *par_exp, const fq_nmod_t coeff) {
    if (fq_nmod_is_zero(coeff, p->ctx)) return;
    
    /* Add new term without checking for duplicates */
    if (p->nterms >= p->alloc) {
        p->alloc = p->alloc ? p->alloc * 2 : 16;
        p->terms = (fq_monomial_t*) flint_realloc(p->terms, p->alloc * sizeof(fq_monomial_t));
    }
    
    /* Initialize coefficient */
    fq_nmod_init(p->terms[p->nterms].coeff, p->ctx);
    fq_nmod_set(p->terms[p->nterms].coeff, coeff, p->ctx);
    
    /* Allocate and copy exponents */
    if (p->nvars > 0) {
        p->terms[p->nterms].var_exp = (slong*) flint_calloc(p->nvars, sizeof(slong));
        if (var_exp) {
            memcpy(p->terms[p->nterms].var_exp, var_exp, p->nvars * sizeof(slong));
        }
    } else {
        p->terms[p->nterms].var_exp = NULL;
    }
    
    if (p->npars > 0) {
        p->terms[p->nterms].par_exp = (slong*) flint_calloc(p->npars, sizeof(slong));
        if (par_exp) {
            memcpy(p->terms[p->nterms].par_exp, par_exp, p->npars * sizeof(slong));
        }
    } else {
        p->terms[p->nterms].par_exp = NULL;
    }
    
    p->nterms++;
}

void fq_mvpoly_copy(fq_mvpoly_t *dest, const fq_mvpoly_t *src) {
    if (dest == src) {
        return;
    }
    
    fq_mvpoly_t temp;
    fq_mvpoly_init(&temp, src->nvars, src->npars, src->ctx);
    
    for (slong i = 0; i < src->nterms; i++) {
        fq_mvpoly_add_term_fast(&temp, src->terms[i].var_exp, src->terms[i].par_exp, src->terms[i].coeff);
    }
    
    *dest = temp;
}

void fq_mvpoly_add_term(fq_mvpoly_t *p, const slong *var_exp, const slong *par_exp, const fq_nmod_t coeff) {
    if (fq_nmod_is_zero(coeff, p->ctx)) return;
    
    /* Check if term already exists */
    for (slong i = 0; i < p->nterms; i++) {
        int same = 1;
        
        /* Check variable exponents */
        if (p->nvars > 0) {
            for (slong j = 0; j < p->nvars; j++) {
                slong e1 = var_exp ? var_exp[j] : 0;
                slong e2 = p->terms[i].var_exp ? p->terms[i].var_exp[j] : 0;
                if (e1 != e2) {
                    same = 0;
                    break;
                }
            }
        }
        
        /* Check parameter exponents */
        if (same && p->npars > 0) {
            for (slong j = 0; j < p->npars; j++) {
                slong e1 = par_exp ? par_exp[j] : 0;
                slong e2 = p->terms[i].par_exp ? p->terms[i].par_exp[j] : 0;
                if (e1 != e2) {
                    same = 0;
                    break;
                }
            }
        }
        
        if (same) {
            fq_nmod_add(p->terms[i].coeff, p->terms[i].coeff, coeff, p->ctx);
            if (fq_nmod_is_zero(p->terms[i].coeff, p->ctx)) {
                /* Remove zero term */
                fq_nmod_clear(p->terms[i].coeff, p->ctx);
                if (p->terms[i].var_exp) flint_free(p->terms[i].var_exp);
                if (p->terms[i].par_exp) flint_free(p->terms[i].par_exp);
                for (slong j = i; j < p->nterms - 1; j++) {
                    p->terms[j] = p->terms[j + 1];
                }
                p->nterms--;
            }
            return;
        }
    }
    
    /* Add new term */
    fq_mvpoly_add_term_fast(p, var_exp, par_exp, coeff);
}

// Output functions
void fq_nmod_print_pretty_enhanced(const fq_nmod_t a, const fq_nmod_ctx_t ctx) {
    if (fq_nmod_is_zero(a, ctx)) {
        printf("0");
        return;
    }
    
    slong degree = fq_nmod_ctx_degree(ctx);
    
    if (degree == 1) {
        nmod_poly_t poly;
        nmod_poly_init(poly, fq_nmod_ctx_prime(ctx));
        fq_nmod_get_nmod_poly(poly, a, ctx);
        
        if (nmod_poly_degree(poly) >= 0) {
            printf("%lu", nmod_poly_get_coeff_ui(poly, 0));
        } else {
            printf("0");
        }
        nmod_poly_clear(poly);
    } else {
        nmod_poly_t poly;
        nmod_poly_init(poly, fq_nmod_ctx_prime(ctx));
        fq_nmod_get_nmod_poly(poly, a, ctx);
        
        slong deg = nmod_poly_degree(poly);
        int first_term = 1;
        
        for (slong i = deg; i >= 0; i--) {
            mp_limb_t coeff = nmod_poly_get_coeff_ui(poly, i);
            if (coeff != 0) {
                if (!first_term) {
                    printf(" + ");
                }
                first_term = 0;
                
                if (i == 0) {
                    printf("%lu", coeff);
                } else if (i == 1) {
                    if (coeff == 1) {
                        printf("t");
                    } else {
                        printf("%lu*t", coeff);
                    }
                } else {
                    if (coeff == 1) {
                        printf("t^%ld", i);
                    } else {
                        printf("%lu*t^%ld", coeff, i);
                    }
                }
            }
        }
        
        if (first_term) {
            printf("0");
        }
        
        nmod_poly_clear(poly);
    }
}


/* ============================================================================
 * Display Functions Implementation
 * ============================================================================ */

  // Enhanced coefficient printing with proper parentheses for extension fields
void fq_nmod_print_with_parentheses(const fq_nmod_t a, const fq_nmod_ctx_t ctx, const char *gen_name) {
    if (fq_nmod_is_zero(a, ctx)) {
        printf("0");
        return;
    }
    
    slong degree = fq_nmod_ctx_degree(ctx);
    
    if (degree == 1) {
        // Prime field - no parentheses needed
        nmod_poly_t poly;
        nmod_poly_init(poly, fq_nmod_ctx_prime(ctx));
        fq_nmod_get_nmod_poly(poly, a, ctx);
        
        if (nmod_poly_degree(poly) >= 0) {
            printf("%lu", nmod_poly_get_coeff_ui(poly, 0));
        } else {
            printf("0");
        }
        nmod_poly_clear(poly);
    } else {
        // Extension field - ALWAYS use parentheses for multi-term expressions
        nmod_poly_t poly;
        nmod_poly_init(poly, fq_nmod_ctx_prime(ctx));
        fq_nmod_get_nmod_poly(poly, a, ctx);
        
        slong deg = nmod_poly_degree(poly);
        
        // Count non-zero terms to decide on parentheses
        int term_count = 0;
        for (slong i = deg; i >= 0; i--) {
            if (nmod_poly_get_coeff_ui(poly, i) != 0) {
                term_count++;
            }
        }
        
        // Use parentheses for multi-term expressions or when explicitly requested
        int use_parens = (term_count > 1);
        
        if (use_parens) printf("(");
        
        int first_term = 1;
        for (slong i = deg; i >= 0; i--) {
            mp_limb_t coeff = nmod_poly_get_coeff_ui(poly, i);
            if (coeff != 0) {
                if (!first_term) {
                    printf(" + ");
                }
                first_term = 0;
                
                if (i == 0) {
                    printf("%lu", coeff);
                } else if (i == 1) {
                    if (coeff == 1) {
                        printf("%s", gen_name ? gen_name : "t");
                    } else {
                        printf("%lu*%s", coeff, gen_name ? gen_name : "t");
                    }
                } else {
                    if (coeff == 1) {
                        printf("%s^%ld", gen_name ? gen_name : "t", i);
                    } else {
                        printf("%lu*%s^%ld", coeff, gen_name ? gen_name : "t", i);
                    }
                }
            }
        }
        
        if (first_term) {
            printf("0");
        }
        
        if (use_parens) printf(")");
        
        nmod_poly_clear(poly);
    }
}

// Main consolidated print function with proper name support
void fq_mvpoly_print_with_names(const fq_mvpoly_t *poly, const char *poly_name,
                               char **var_names, char **par_names, 
                               const char *gen_name, int expanded_format) {
    
    // Print polynomial name if provided
    if (poly_name && strlen(poly_name) > 0) {
        printf("%s = ", poly_name);
    }
    
    if (poly->nterms == 0) {
        printf("0\n");
        return;
    }
    
    // Default names if not provided
    char default_var_names[] = {'x', 'y', 'z', 'w', 'v', 'u'};
    char default_par_names[] = {'a', 'b', 'c', 'd'};
    
    for (slong i = 0; i < poly->nterms; i++) {
        if (i > 0) printf(" + ");
        
        // Check if we have variables or parameters for this term
        int has_vars_or_pars = 0;
        
        // Check variables
        if (poly->nvars > 0 && poly->terms[i].var_exp) {
            for (slong j = 0; j < poly->nvars; j++) {
                if (poly->terms[i].var_exp[j] > 0) {
                    has_vars_or_pars = 1;
                    break;
                }
            }
        }
        
        // Check parameters
        if (!has_vars_or_pars && poly->npars > 0 && poly->terms[i].par_exp) {
            for (slong j = 0; j < poly->npars; j++) {
                if (poly->terms[i].par_exp[j] > 0) {
                    has_vars_or_pars = 1;
                    break;
                }
            }
        }
        
        // Handle coefficient
        fq_nmod_t one;
        fq_nmod_init(one, poly->ctx);
        fq_nmod_one(one, poly->ctx);
        
        if (fq_nmod_is_one(poly->terms[i].coeff, poly->ctx) && has_vars_or_pars) {
            // Coefficient is 1 and we have variables/parameters - don't print coefficient
        } else {
            // Print coefficient with proper parentheses
            fq_nmod_print_with_parentheses(poly->terms[i].coeff, poly->ctx, gen_name);
            
            if (has_vars_or_pars) {
                printf("*");
            }
        }
        
        fq_nmod_clear(one, poly->ctx);
        
        // Track whether we've printed anything for this term yet
        int term_printed = 0;
        
        // Print variables
        // Fix for fq_mvpoly_print_with_names function
        // Replace the variable printing sections:
        
        // Print variables
        if (expanded_format && poly->nvars > 0 && poly->nvars % 2 == 0) {
            // Expanded format with dual variables (for Dixon polynomials)
            slong actual_nvars = poly->nvars / 2;
            
            // Regular variables
            for (slong j = 0; j < actual_nvars; j++) {
                if (poly->terms[i].var_exp && poly->terms[i].var_exp[j] > 0) {
                    // FIXED: Only add * if something was already printed for this term
                    if (term_printed) {
                        printf("*");
                    }
                    
                    // Use provided name or default
                    if (var_names && var_names[j]) {
                        printf("%s", var_names[j]);
                    } else if (j < 6) {
                        printf("%c", default_var_names[j]);
                    } else {
                        printf("x_%ld", j);
                    }
                    
                    if (poly->terms[i].var_exp[j] > 1) {
                        printf("^%ld", poly->terms[i].var_exp[j]);
                    }
                    term_printed = 1;
                }
            }
            
            // Dual variables with tilde
            for (slong j = actual_nvars; j < poly->nvars; j++) {
                if (poly->terms[i].var_exp && poly->terms[i].var_exp[j] > 0) {
                    // FIXED: Only add * if something was already printed for this term
                    if (term_printed) {
                        printf("*");
                    }
                    slong orig_idx = j - actual_nvars;
                    
                    if (var_names && var_names[orig_idx]) {
                        printf("~%s", var_names[orig_idx]);
                    } else if (orig_idx < 6) {
                        printf("~%c", default_var_names[orig_idx]);
                    } else {
                        printf("~x_%ld", orig_idx);
                    }
                    
                    if (poly->terms[i].var_exp[j] > 1) {
                        printf("^%ld", poly->terms[i].var_exp[j]);
                    }
                    term_printed = 1;
                }
            }
        } else {
            // Normal variable format
            for (slong j = 0; j < poly->nvars; j++) {
                if (poly->terms[i].var_exp && poly->terms[i].var_exp[j] > 0) {
                    // FIXED: Only add * if something was already printed for this term
                    if (term_printed) {
                        printf("*");
                    }
                    
                    if (var_names && var_names[j]) {
                        printf("%s", var_names[j]);
                    } else if (j < 6) {
                        printf("%c", default_var_names[j]);
                    } else {
                        printf("x_%ld", j);
                    }
                    
                    if (poly->terms[i].var_exp[j] > 1) {
                        printf("^%ld", poly->terms[i].var_exp[j]);
                    }
                    term_printed = 1;
                }
            }
        }
        
        // Print parameters with actual names
        for (slong j = 0; j < poly->npars; j++) {
            if (poly->terms[i].par_exp && poly->terms[i].par_exp[j] > 0) {
                // FIXED: Only add * if something was already printed for this term
                if (term_printed) {
                    printf("*");
                }
                
                // Use provided parameter name or default
                if (par_names && par_names[j]) {
                    printf("%s", par_names[j]);
                } else if (j < 4) {
                    printf("%c", default_par_names[j]);
                } else {
                    printf("p_%ld", j);
                }
                
                if (poly->terms[i].par_exp[j] > 1) {
                    printf("^%ld", poly->terms[i].par_exp[j]);
                }
                term_printed = 1;
            }
        }
    }
    printf("\n");
}


// Convenience wrapper functions
void fq_mvpoly_print(const fq_mvpoly_t *poly, const char *name) {
    fq_mvpoly_print_with_names(poly, name, NULL, NULL, NULL, 0);
}

void fq_mvpoly_print_expanded(const fq_mvpoly_t *poly, const char *name, int use_dual) {
    fq_mvpoly_print_with_names(poly, name, NULL, NULL, NULL, use_dual);
}

// Enhanced version that accepts name arrays
void fq_mvpoly_print_enhanced(const fq_mvpoly_t *poly, const char *name,
                            char **var_names, char **par_names, const char *gen_name) {
    fq_mvpoly_print_with_names(poly, name, var_names, par_names, gen_name, 0);
}


/* ============================================================================
 * Arithmetic Operations Implementation
 * ============================================================================ */

void fq_mvpoly_mul(fq_mvpoly_t *result, const fq_mvpoly_t *a, const fq_mvpoly_t *b) {
    fq_mvpoly_t temp;
    fq_mvpoly_init(&temp, a->nvars, a->npars, a->ctx);
    
    for (slong i = 0; i < a->nterms; i++) {
        for (slong j = 0; j < b->nterms; j++) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, a->ctx);
            fq_nmod_mul(coeff, a->terms[i].coeff, b->terms[j].coeff, a->ctx);
            
            slong *var_exp = NULL;
            slong *par_exp = NULL;
            
            if (a->nvars > 0) {
                var_exp = (slong*) flint_calloc(a->nvars, sizeof(slong));
                for (slong k = 0; k < a->nvars; k++) {
                    var_exp[k] = (a->terms[i].var_exp ? a->terms[i].var_exp[k] : 0) +
                                 (b->terms[j].var_exp ? b->terms[j].var_exp[k] : 0);
                }
            }
            
            if (a->npars > 0) {
                par_exp = (slong*) flint_calloc(a->npars, sizeof(slong));
                for (slong k = 0; k < a->npars; k++) {
                    par_exp[k] = (a->terms[i].par_exp ? a->terms[i].par_exp[k] : 0) +
                                 (b->terms[j].par_exp ? b->terms[j].par_exp[k] : 0);
                }
            }
            
            fq_mvpoly_add_term(&temp, var_exp, par_exp, coeff);
            
            fq_nmod_clear(coeff, a->ctx);
            if (var_exp) flint_free(var_exp);
            if (par_exp) flint_free(par_exp);
        }
    }
    
    /* Handle case where result might already be initialized */
    if (result == a || result == b) {
        fq_mvpoly_clear(result);
    }
    *result = temp;

    if (g_field_equation_reduction)
        fq_mvpoly_reduce_field_equation(result);
}

void fq_mvpoly_pow(fq_mvpoly_t *result, const fq_mvpoly_t *base, slong power) {
    fq_mvpoly_t temp;
    fq_mvpoly_init(&temp, base->nvars, base->npars, base->ctx);
    
    fq_nmod_t one;
    fq_nmod_init(one, base->ctx);
    fq_nmod_one(one, base->ctx);
    fq_mvpoly_add_term(&temp, NULL, NULL, one); // temp = 1
    fq_nmod_clear(one, base->ctx);
    
    if (power == 0) {
        *result = temp;
        return;
    }
    
    fq_mvpoly_t base_copy, temp2;
    fq_mvpoly_copy(&base_copy, base);
    
    while (power > 0) {
        if (power & 1) {
            fq_mvpoly_mul(&temp2, &temp, &base_copy);
            fq_mvpoly_clear(&temp);
            temp = temp2;  // struct assignment
        }
        
        if (power > 1) {
            fq_mvpoly_mul(&temp2, &base_copy, &base_copy);
            fq_mvpoly_clear(&base_copy);
            base_copy = temp2;  // struct assignment
        }
        
        power >>= 1;
    }
    
    fq_mvpoly_clear(&base_copy);
    
    /* Handle case where result might already be initialized */
    if (result == base) {
        fq_mvpoly_clear(result);
    }
    *result = temp;
}

void fq_mvpoly_scalar_mul(fq_mvpoly_t *result, const fq_mvpoly_t *p, const fq_nmod_t scalar) {
    fq_mvpoly_t temp;
    fq_mvpoly_init(&temp, p->nvars, p->npars, p->ctx);
    
    for (slong i = 0; i < p->nterms; i++) {
        fq_nmod_t new_coeff;
        fq_nmod_init(new_coeff, p->ctx);
        fq_nmod_mul(new_coeff, p->terms[i].coeff, scalar, p->ctx);
        fq_mvpoly_add_term_fast(&temp, p->terms[i].var_exp, p->terms[i].par_exp, new_coeff);
        fq_nmod_clear(new_coeff, p->ctx);
    }
    
    /* Handle case where result might already be initialized */
    if (result == p) {
        fq_mvpoly_clear(result);
    }
    *result = temp;
}

/* Optimized fq_mvpoly_sub function - uses fq_nmod_mpoly for efficient subtraction */
void fq_mvpoly_sub_optimized(fq_mvpoly_t *result, const fq_mvpoly_t *a, const fq_mvpoly_t *b) {
    /* Handle special cases */
    if (a->nterms == 0 && b->nterms == 0) {
        if (result->terms != NULL) {
            fq_mvpoly_clear(result);
        }
        fq_mvpoly_init(result, a->nvars, a->npars, a->ctx);
        return;
    }
    
    if (a->nterms == 0) {
        /* result = -b */
        fq_mvpoly_t temp;
        fq_mvpoly_init(&temp, b->nvars, b->npars, b->ctx);
        
        for (slong i = 0; i < b->nterms; i++) {
            fq_nmod_t neg_coeff;
            fq_nmod_init(neg_coeff, b->ctx);
            fq_nmod_neg(neg_coeff, b->terms[i].coeff, b->ctx);
            fq_mvpoly_add_term_fast(&temp, b->terms[i].var_exp, b->terms[i].par_exp, neg_coeff);
            fq_nmod_clear(neg_coeff, b->ctx);
        }
        
        if (result->terms != NULL) {
            fq_mvpoly_clear(result);
        }
        *result = temp;
        return;
    }
    
    if (b->nterms == 0) {
        /* result = a */
        if (result == a) {
            return; // already correct result
        }
        
        if (result->terms != NULL) {
            fq_mvpoly_clear(result);
        }
        fq_mvpoly_copy(result, a);
        return;
    }
    
    /* Create mpoly context */
    slong total_vars = a->nvars + a->npars;
    fq_nmod_mpoly_ctx_t mpoly_ctx;
    fq_nmod_mpoly_ctx_init(mpoly_ctx, total_vars, ORD_LEX, a->ctx);
    
    /* Create mpoly variables */
    fq_nmod_mpoly_t mpoly_a, mpoly_b, mpoly_result;
    fq_nmod_mpoly_init(mpoly_a, mpoly_ctx);
    fq_nmod_mpoly_init(mpoly_b, mpoly_ctx);
    fq_nmod_mpoly_init(mpoly_result, mpoly_ctx);
    
    /* Convert to mpoly format */
    fq_mvpoly_to_fq_nmod_mpoly(mpoly_a, a, mpoly_ctx);
    fq_mvpoly_to_fq_nmod_mpoly(mpoly_b, b, mpoly_ctx);
    
    /* Execute efficient subtraction */
    fq_nmod_mpoly_sub(mpoly_result, mpoly_a, mpoly_b, mpoly_ctx);
    
    /* Convert back to fq_mvpoly format */
    fq_mvpoly_t temp;
    fq_nmod_mpoly_to_fq_mvpoly(&temp, mpoly_result, a->nvars, a->npars, mpoly_ctx, a->ctx);
    
    /* Handle result */
    if (result == a || result == b) {
        fq_mvpoly_clear(result);
    } else if (result->terms != NULL) {
        fq_mvpoly_clear(result);
    }
    *result = temp;
    
    /* Cleanup */
    fq_nmod_mpoly_clear(mpoly_a, mpoly_ctx);
    fq_nmod_mpoly_clear(mpoly_b, mpoly_ctx);
    fq_nmod_mpoly_clear(mpoly_result, mpoly_ctx);
    fq_nmod_mpoly_ctx_clear(mpoly_ctx);
}

/* Also optimize fq_mvpoly_add function */
void fq_mvpoly_add(fq_mvpoly_t *result, const fq_mvpoly_t *a, const fq_mvpoly_t *b) {
    if (a->nterms == 0 && b->nterms == 0) {
        if (result->terms != NULL) {
            fq_mvpoly_clear(result);
        }
        fq_mvpoly_init(result, a->nvars, a->npars, a->ctx);
        return;
    }
    
    if (a->nterms == 0) {
        if (result == b) return;
        if (result->terms != NULL) {
            fq_mvpoly_clear(result);
        }
        fq_mvpoly_copy(result, b);
        return;
    }
    
    if (b->nterms == 0) {
        if (result == a) return;
        if (result->terms != NULL) {
            fq_mvpoly_clear(result);
        }
        fq_mvpoly_copy(result, a);
        return;
    }
    
    slong total_vars = a->nvars + a->npars;
    fq_nmod_mpoly_ctx_t mpoly_ctx;
    fq_nmod_mpoly_ctx_init(mpoly_ctx, total_vars, ORD_LEX, a->ctx);
    
    fq_nmod_mpoly_t mpoly_a, mpoly_b, mpoly_result;
    fq_nmod_mpoly_init(mpoly_a, mpoly_ctx);
    fq_nmod_mpoly_init(mpoly_b, mpoly_ctx);
    fq_nmod_mpoly_init(mpoly_result, mpoly_ctx);
    
    fq_mvpoly_to_fq_nmod_mpoly(mpoly_a, a, mpoly_ctx);
    fq_mvpoly_to_fq_nmod_mpoly(mpoly_b, b, mpoly_ctx);
    
    fq_nmod_mpoly_add(mpoly_result, mpoly_a, mpoly_b, mpoly_ctx);
    
    fq_mvpoly_t temp;
    fq_nmod_mpoly_to_fq_mvpoly(&temp, mpoly_result, a->nvars, a->npars, mpoly_ctx, a->ctx);

    if (result == a || result == b) {
        fq_mvpoly_clear(result);
    }
    *result = temp;
    
    fq_nmod_mpoly_clear(mpoly_a, mpoly_ctx);
    fq_nmod_mpoly_clear(mpoly_b, mpoly_ctx);
    fq_nmod_mpoly_clear(mpoly_result, mpoly_ctx);
    fq_nmod_mpoly_ctx_clear(mpoly_ctx);
}

void fq_mvpoly_sub(fq_mvpoly_t *result, const fq_mvpoly_t *a, const fq_mvpoly_t *b) {

    if (a->nterms == 0 && b->nterms == 0) {
        if (result->terms != NULL) {
            fq_mvpoly_clear(result);
        }
        fq_mvpoly_init(result, a->nvars, a->npars, a->ctx);
        return;
    }
    
    if (a->nterms == 0) {
        fq_mvpoly_t temp;
        fq_mvpoly_init(&temp, b->nvars, b->npars, b->ctx);
        
        for (slong i = 0; i < b->nterms; i++) {
            fq_nmod_t neg_coeff;
            fq_nmod_init(neg_coeff, b->ctx);
            fq_nmod_neg(neg_coeff, b->terms[i].coeff, b->ctx);
            fq_mvpoly_add_term_fast(&temp, b->terms[i].var_exp, b->terms[i].par_exp, neg_coeff);
            fq_nmod_clear(neg_coeff, b->ctx);
        }
        
        if (result->terms != NULL) {
            fq_mvpoly_clear(result);
        }
        *result = temp;
        return;
    }
    
    if (b->nterms == 0) {
        if (result == a) return;
        
        if (result->terms != NULL) {
            fq_mvpoly_clear(result);
        }
        fq_mvpoly_copy(result, a);
        return;
    }

    slong total_vars = a->nvars + a->npars;
    fq_nmod_mpoly_ctx_t mpoly_ctx;
    fq_nmod_mpoly_ctx_init(mpoly_ctx, total_vars, ORD_LEX, a->ctx);
    
    fq_nmod_mpoly_t mpoly_a, mpoly_b, mpoly_result;
    fq_nmod_mpoly_init(mpoly_a, mpoly_ctx);
    fq_nmod_mpoly_init(mpoly_b, mpoly_ctx);
    fq_nmod_mpoly_init(mpoly_result, mpoly_ctx);
    
    fq_mvpoly_to_fq_nmod_mpoly(mpoly_a, a, mpoly_ctx);
    fq_mvpoly_to_fq_nmod_mpoly(mpoly_b, b, mpoly_ctx);
    
    fq_nmod_mpoly_sub(mpoly_result, mpoly_a, mpoly_b, mpoly_ctx);

    fq_mvpoly_t temp;
    fq_nmod_mpoly_to_fq_mvpoly(&temp, mpoly_result, a->nvars, a->npars, mpoly_ctx, a->ctx);

    if (result == a || result == b) {
        fq_mvpoly_clear(result);
    }
    *result = temp;

    fq_nmod_mpoly_clear(mpoly_a, mpoly_ctx);
    fq_nmod_mpoly_clear(mpoly_b, mpoly_ctx);
    fq_nmod_mpoly_clear(mpoly_result, mpoly_ctx);
    fq_nmod_mpoly_ctx_clear(mpoly_ctx);
}

void print_fq_matrix_mvpoly(fq_mvpoly_t **matrix, slong nrows, slong ncols, 
                            const char *matrix_name, int show_details) {
    printf("\n=== %s Matrix (%ld x %ld) ===\n", matrix_name, nrows, ncols);
    
    if (show_details) {
        for (slong i = 0; i < nrows; i++) {
            for (slong j = 0; j < ncols; j++) {
                printf("M[%ld][%ld]: ", i, j);
                if (matrix[i][j].nterms == 0) {
                    printf("0\n");
                } else if (matrix[i][j].nterms <= 10) {
                    fq_mvpoly_print_expanded(&matrix[i][j], "", 1);
                } else {
                    printf("%ld terms", matrix[i][j].nterms);
                    
                    slong max_var_deg = 0, max_par_deg = 0;
                    for (slong t = 0; t < matrix[i][j].nterms; t++) {
                        if (matrix[i][j].terms[t].var_exp) {
                            for (slong k = 0; k < matrix[i][j].nvars; k++) {
                                if (matrix[i][j].terms[t].var_exp[k] > max_var_deg) {
                                    max_var_deg = matrix[i][j].terms[t].var_exp[k];
                                }
                            }
                        }
                        if (matrix[i][j].terms[t].par_exp && matrix[i][j].npars > 0) {
                            for (slong k = 0; k < matrix[i][j].npars; k++) {
                                if (matrix[i][j].terms[t].par_exp[k] > max_par_deg) {
                                    max_par_deg = matrix[i][j].terms[t].par_exp[k];
                                }
                            }
                        }
                    }
                    printf(" (max var deg: %ld, max par deg: %ld)\n", max_var_deg, max_par_deg);
                }
            }
        }
    } else {
        printf("Matrix term counts:\n");
        for (slong i = 0; i < nrows; i++) {
            printf("  Row %ld: ", i);
            for (slong j = 0; j < ncols; j++) {
                printf("%3ld", matrix[i][j].nterms);
                if (j < ncols - 1) printf(" ");
            }
            printf("\n");
        }
    }
    printf("=== End %s Matrix ===\n\n", matrix_name);
}

void analyze_fq_matrix_mvpoly(fq_mvpoly_t **matrix, slong nrows, slong ncols, 
                              const char *matrix_name) {
    printf("\n--- Analysis of %s Matrix ---\n", matrix_name);
    
    slong total_terms = 0;
    slong zero_entries = 0;
    slong max_terms = 0;
    slong max_var_deg = 0;
    slong max_par_deg = 0;
    
    for (slong i = 0; i < nrows; i++) {
        for (slong j = 0; j < ncols; j++) {
            total_terms += matrix[i][j].nterms;
            if (matrix[i][j].nterms == 0) {
                zero_entries++;
            }
            if (matrix[i][j].nterms > max_terms) {
                max_terms = matrix[i][j].nterms;
            }
            
            for (slong t = 0; t < matrix[i][j].nterms; t++) {
                if (matrix[i][j].terms[t].var_exp) {
                    for (slong k = 0; k < matrix[i][j].nvars; k++) {
                        if (matrix[i][j].terms[t].var_exp[k] > max_var_deg) {
                            max_var_deg = matrix[i][j].terms[t].var_exp[k];
                        }
                    }
                }
                if (matrix[i][j].terms[t].par_exp && matrix[i][j].npars > 0) {
                    for (slong k = 0; k < matrix[i][j].npars; k++) {
                        if (matrix[i][j].terms[t].par_exp[k] > max_par_deg) {
                            max_par_deg = matrix[i][j].terms[t].par_exp[k];
                        }
                    }
                }
            }
        }
    }
    
    printf("  Total entries: %ld\n", nrows * ncols);
    printf("  Zero entries: %ld (%.1f%%)\n", zero_entries, 
           100.0 * zero_entries / (nrows * ncols));
    printf("  Non-zero entries: %ld\n", nrows * ncols - zero_entries);
    printf("  Total terms: %ld\n", total_terms);
    printf("  Average terms per entry: %.2f\n", (double)total_terms / (nrows * ncols));
    printf("  Maximum terms in single entry: %ld\n", max_terms);
    printf("  Maximum variable degree: %ld\n", max_var_deg);
    printf("  Maximum parameter degree: %ld\n", max_par_deg);
    printf("--- End Analysis ---\n\n");
}


// ============ Kronecker+HNF helpers ============

slong exp_to_kronecker_index(const slong *exp, const slong *degs, slong n) {
    slong index = 0;
    slong stride = 1;
    
    for (slong i = 0; i < n; i++) {
        index += exp[i] * stride;
        stride *= (degs[i] + 1);
    }
    
    return index;
}

void kronecker_index_to_exp(slong index, slong *exp, const slong *degs, slong n) {
    for (slong i = 0; i < n; i++) {
        exp[i] = index % (degs[i] + 1);
        index /= (degs[i] + 1);
    }
}

// Convert fq_mvpoly to univariate via Kronecker
void fq_mvpoly_to_kronecker_full(fq_nmod_poly_t out, const fq_mvpoly_t *p, 
                                const slong *var_degs, const slong *par_degs) {
    // Calculate total degree
    slong total_vars = p->nvars + p->npars;
    slong *all_degs = (slong*) flint_malloc(total_vars * sizeof(slong));
    
    // Copy degree bounds
    for (slong i = 0; i < p->nvars; i++) {
        all_degs[i] = var_degs[i];
    }
    for (slong i = 0; i < p->npars; i++) {
        all_degs[p->nvars + i] = par_degs[i];
    }
    
    // Calculate maximum degree
    slong max_deg = 1;
    for (slong i = 0; i < total_vars; i++) {
        max_deg *= (all_degs[i] + 1);
    }
    max_deg--;
    
    fq_nmod_poly_fit_length(out, max_deg + 1, p->ctx);
    fq_nmod_poly_zero(out, p->ctx);
    
    // Convert each term
    slong *combined_exp = (slong*) flint_calloc(total_vars, sizeof(slong));
    
    for (slong i = 0; i < p->nterms; i++) {
        // Combine variable and parameter exponents
        for (slong j = 0; j < p->nvars; j++) {
            combined_exp[j] = p->terms[i].var_exp ? p->terms[i].var_exp[j] : 0;
        }
        for (slong j = 0; j < p->npars; j++) {
            combined_exp[p->nvars + j] = p->terms[i].par_exp ? p->terms[i].par_exp[j] : 0;
        }
        
        slong idx = exp_to_kronecker_index(combined_exp, all_degs, total_vars);
        
        fq_nmod_t old_coeff, new_coeff;
        fq_nmod_init(old_coeff, p->ctx);
        fq_nmod_init(new_coeff, p->ctx);
        
        fq_nmod_poly_get_coeff(old_coeff, out, idx, p->ctx);
        fq_nmod_add(new_coeff, old_coeff, p->terms[i].coeff, p->ctx);
        fq_nmod_poly_set_coeff(out, idx, new_coeff, p->ctx);
        
        fq_nmod_clear(old_coeff, p->ctx);
        fq_nmod_clear(new_coeff, p->ctx);
    }
    
    _fq_nmod_poly_normalise(out, p->ctx);
    
    flint_free(combined_exp);
    flint_free(all_degs);
}

void kronecker_to_fq_mvpoly_full(fq_mvpoly_t *out, const fq_nmod_poly_t in, 
                                const slong *var_degs, slong nvars,
                                const slong *par_degs, slong npars, const fq_nmod_ctx_t ctx) {
    fq_mvpoly_init(out, nvars, npars, ctx);
    
    slong deg = fq_nmod_poly_degree(in, ctx);
    if (deg < 0) return;
    
    // Pre-calculate non-zero terms
    slong nterms = 0;
    fq_nmod_t coeff;
    fq_nmod_init(coeff, ctx);
    
    for (slong i = 0; i <= deg; i++) {
        fq_nmod_poly_get_coeff(coeff, in, i, ctx);
        if (!fq_nmod_is_zero(coeff, ctx)) {
            nterms++;
        }
    }
    
    if (nterms == 0) {
        fq_nmod_clear(coeff, ctx);
        return;
    }
    
    // Pre-allocate space
    if (out->alloc < nterms) {
        out->alloc = nterms;
        out->terms = (fq_monomial_t*) flint_realloc(out->terms, out->alloc * sizeof(fq_monomial_t));
    }
    
    // Prepare work arrays
    slong total_vars = nvars + npars;
    slong *all_degs = (slong*) flint_malloc(total_vars * sizeof(slong));
    
    for (slong i = 0; i < nvars; i++) {
        all_degs[i] = var_degs[i];
    }
    for (slong i = 0; i < npars; i++) {
        all_degs[nvars + i] = par_degs[i];
    }
    
    slong *combined_exp = (slong*) flint_calloc(total_vars, sizeof(slong));
    
    // Fill terms
    slong term_idx = 0;
    for (slong i = 0; i <= deg; i++) {
        fq_nmod_poly_get_coeff(coeff, in, i, ctx);
        if (!fq_nmod_is_zero(coeff, ctx)) {
            // Decode exponents
            slong idx = i;
            for (slong j = 0; j < total_vars; j++) {
                combined_exp[j] = idx % (all_degs[j] + 1);
                idx /= (all_degs[j] + 1);
            }
            
            // Initialize and set coefficient
            fq_nmod_init(out->terms[term_idx].coeff, ctx);
            fq_nmod_set(out->terms[term_idx].coeff, coeff, ctx);
            
            if (nvars > 0) {
                out->terms[term_idx].var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
                for (slong j = 0; j < nvars; j++) {
                    out->terms[term_idx].var_exp[j] = combined_exp[j];
                }
            } else {
                out->terms[term_idx].var_exp = NULL;
            }
            
            if (npars > 0) {
                out->terms[term_idx].par_exp = (slong*) flint_calloc(npars, sizeof(slong));
                for (slong j = 0; j < npars; j++) {
                    out->terms[term_idx].par_exp[j] = combined_exp[nvars + j];
                }
            } else {
                out->terms[term_idx].par_exp = NULL;
            }
            
            term_idx++;
        }
    }
    
    out->nterms = term_idx;
    
    fq_nmod_clear(coeff, ctx);
    flint_free(combined_exp);
    flint_free(all_degs);
}


// Improved division by linear factor (x_i - ~x_i) - direct method
void divide_by_fq_linear_factor_improved(fq_mvpoly_t *quotient, const fq_mvpoly_t *dividend,
                                        slong var_idx, slong nvars, slong npars) {
    fq_mvpoly_t temp;
    fq_mvpoly_init(&temp, nvars, npars, dividend->ctx);
    
    if (dividend->nterms == 0) {
        *quotient = temp;
        return;
    }
    
    slong actual_nvars = nvars / 2;  // Since we have dual variables
    
    // Process each term in the dividend
    for (slong i = 0; i < dividend->nterms; i++) {
        slong deg_x = dividend->terms[i].var_exp ? dividend->terms[i].var_exp[var_idx] : 0;
        slong deg_dual = dividend->terms[i].var_exp ? dividend->terms[i].var_exp[var_idx + actual_nvars] : 0;
        
        if (deg_x == deg_dual) {
            // This term is divisible by (x_i - ~x_i), the result has the same exponents
            fq_mvpoly_add_term(&temp, dividend->terms[i].var_exp, dividend->terms[i].par_exp, 
                              dividend->terms[i].coeff);
        } else if (deg_x > deg_dual) {
            // Reduce x degree by 1
            slong *new_exp = (slong*) flint_calloc(nvars, sizeof(slong));
            if (dividend->terms[i].var_exp) {
                memcpy(new_exp, dividend->terms[i].var_exp, nvars * sizeof(slong));
            }
            new_exp[var_idx]--;
            
            fq_mvpoly_add_term(&temp, new_exp, dividend->terms[i].par_exp, 
                              dividend->terms[i].coeff);
            flint_free(new_exp);
        } else { // deg_dual > deg_x
            // Reduce ~x degree by 1, with negative coefficient
            slong *new_exp = (slong*) flint_calloc(nvars, sizeof(slong));
            if (dividend->terms[i].var_exp) {
                memcpy(new_exp, dividend->terms[i].var_exp, nvars * sizeof(slong));
            }
            new_exp[var_idx + actual_nvars]--;
            
            fq_nmod_t neg_coeff;
            fq_nmod_init(neg_coeff, dividend->ctx);
            fq_nmod_neg(neg_coeff, dividend->terms[i].coeff, dividend->ctx);
            
            fq_mvpoly_add_term(&temp, new_exp, dividend->terms[i].par_exp, neg_coeff);
            
            fq_nmod_clear(neg_coeff, dividend->ctx);
            flint_free(new_exp);
        }
    }
    
    //if (quotient == dividend) {
        //fq_mvpoly_clear(quotient);
    //}
    if (quotient->terms != NULL) {
        fq_mvpoly_clear(quotient);
    }
    *quotient = temp;
}


// Division using FLINT's built-in functions (alternative approach)
void divide_by_fq_linear_factor_flint(fq_mvpoly_t *quotient, const fq_mvpoly_t *dividend,
                                     slong var_idx, slong nvars, slong npars) {
    fq_mvpoly_init(quotient, nvars, npars, dividend->ctx);
    
    if (dividend->nterms == 0) return;
    
    slong total_vars = nvars + npars;
    fq_nmod_mpoly_ctx_t mpoly_ctx;
    fq_nmod_mpoly_ctx_init(mpoly_ctx, total_vars, ORD_LEX, dividend->ctx);
    
    fq_nmod_mpoly_t dividend_mpoly;
    fq_nmod_mpoly_init(dividend_mpoly, mpoly_ctx);
    fq_mvpoly_to_fq_nmod_mpoly(dividend_mpoly, dividend, mpoly_ctx);
    
    fq_nmod_mpoly_t divisor;
    fq_nmod_mpoly_init(divisor, mpoly_ctx);
    
    slong actual_nvars = nvars / 2;
    
    fq_nmod_t one, neg_one;
    fq_nmod_init(one, dividend->ctx);
    fq_nmod_init(neg_one, dividend->ctx);
    fq_nmod_one(one, dividend->ctx);
    fq_nmod_neg(neg_one, one, dividend->ctx);
    
    ulong *exps = (ulong*) flint_calloc(total_vars, sizeof(ulong));
    
    memset(exps, 0, total_vars * sizeof(ulong));
    exps[var_idx] = 1;
    fq_nmod_mpoly_push_term_fq_nmod_ui(divisor, one, exps, mpoly_ctx);
    
    memset(exps, 0, total_vars * sizeof(ulong));
    exps[var_idx + actual_nvars] = 1;
    fq_nmod_mpoly_push_term_fq_nmod_ui(divisor, neg_one, exps, mpoly_ctx);
    
    fq_nmod_mpoly_sort_terms(divisor, mpoly_ctx);
    fq_nmod_mpoly_combine_like_terms(divisor, mpoly_ctx);
    
    fq_nmod_mpoly_t quotient_mpoly;
    fq_nmod_mpoly_init(quotient_mpoly, mpoly_ctx);
    
    int divides = fq_nmod_mpoly_divides(quotient_mpoly, dividend_mpoly, divisor, mpoly_ctx);
    
    if (divides) {

        fq_mvpoly_clear(quotient);
        fq_nmod_mpoly_to_fq_mvpoly(quotient, quotient_mpoly, nvars, npars, mpoly_ctx, dividend->ctx);
    } else {
        fq_nmod_mpoly_t remainder;
        fq_nmod_mpoly_init(remainder, mpoly_ctx);
        
        fq_nmod_mpoly_divrem(quotient_mpoly, remainder, dividend_mpoly, divisor, mpoly_ctx);
        
        if (fq_nmod_mpoly_is_zero(remainder, mpoly_ctx)) {
            fq_mvpoly_clear(quotient);
            fq_nmod_mpoly_to_fq_mvpoly(quotient, quotient_mpoly, nvars, npars, mpoly_ctx, dividend->ctx);
        } else {
            printf("    Warning: FLINT division not exact, falling back to direct method\n");
            fq_mvpoly_clear(quotient);
            divide_by_fq_linear_factor_improved(quotient, dividend, var_idx, nvars, npars);
        }
        
        fq_nmod_mpoly_clear(remainder, mpoly_ctx);
    }
    
    flint_free(exps);
    fq_nmod_clear(one, dividend->ctx);
    fq_nmod_clear(neg_one, dividend->ctx);
    fq_nmod_mpoly_clear(dividend_mpoly, mpoly_ctx);
    fq_nmod_mpoly_clear(divisor, mpoly_ctx);
    fq_nmod_mpoly_clear(quotient_mpoly, mpoly_ctx);
    fq_nmod_mpoly_ctx_clear(mpoly_ctx);
}


// ============ Evaluation functions ============

void evaluate_fq_mvpoly_at_params(fq_nmod_t result, const fq_mvpoly_t *poly, const fq_nmod_t *param_vals) {
    // Create unified context with auto-detection
    field_ctx_t unified_ctx;
    field_ctx_init(&unified_ctx, poly->ctx);
    
    // Get context pointer
    void *ctx_ptr = NULL;
    switch (unified_ctx.field_id) {
        case FIELD_ID_NMOD:
            ctx_ptr = (void*)&unified_ctx.ctx.nmod_ctx;
            break;
        case FIELD_ID_FQ_ZECH:
            ctx_ptr = (void*)unified_ctx.ctx.zech_ctx;
            break;
        default:
            ctx_ptr = (void*)unified_ctx.ctx.fq_ctx;
            break;
    }
    
    // Initialize unified result to zero
    field_elem_u unified_result;
    field_init_elem(&unified_result, unified_ctx.field_id, ctx_ptr);
    field_set_zero(&unified_result, unified_ctx.field_id, ctx_ptr);
    
    // Convert parameters to unified format
    field_elem_u *unified_params = NULL;
    if (poly->npars > 0) {
        unified_params = (field_elem_u*)malloc(poly->npars * sizeof(field_elem_u));
        for (slong i = 0; i < poly->npars; i++) {
            field_init_elem(&unified_params[i], unified_ctx.field_id, ctx_ptr);
            fq_nmod_to_field_elem(&unified_params[i], param_vals[i], &unified_ctx);
        }
    }
    
    // Same evaluation logic, but with unified operations
    for (slong i = 0; i < poly->nterms; i++) {
        field_elem_u term_val;
        field_init_elem(&term_val, unified_ctx.field_id, ctx_ptr);
        
        // Convert coefficient to unified format
        fq_nmod_to_field_elem(&term_val, poly->terms[i].coeff, &unified_ctx);
        
        // Multiply by parameter powers
        if (poly->npars > 0 && poly->terms[i].par_exp) {
            for (slong j = 0; j < poly->npars; j++) {
                for (slong k = 0; k < poly->terms[i].par_exp[j]; k++) {
                    field_mul(&term_val, &term_val, &unified_params[j], 
                             unified_ctx.field_id, ctx_ptr);
                }
            }
        }
        
        // Add to result
        field_add(&unified_result, &unified_result, &term_val, 
                 unified_ctx.field_id, ctx_ptr);
        field_clear_elem(&term_val, unified_ctx.field_id, ctx_ptr);
    }
    
    // Convert result back to fq_nmod
    field_elem_to_fq_nmod(result, &unified_result, &unified_ctx);
    
    // Cleanup
    field_clear_elem(&unified_result, unified_ctx.field_id, ctx_ptr);
    if (unified_params) {
        for (slong i = 0; i < poly->npars; i++) {
            field_clear_elem(&unified_params[i], unified_ctx.field_id, ctx_ptr);
        }
        free(unified_params);
    }
    field_ctx_clear(&unified_ctx);
}

void evaluate_fq_mvpoly_safe(fq_nmod_t result, fq_mvpoly_t ***matrix, 
                             slong i, slong j, const fq_nmod_t *param_vals,
                             slong npars, const fq_nmod_ctx_t ctx) {
    if (matrix[i][j] == NULL) {
        fq_nmod_zero(result, ctx);
    } else {
        evaluate_fq_mvpoly_at_params(result, matrix[i][j], param_vals);
    }
}

void fq_mvpoly_make_monic(fq_mvpoly_t *poly) {
    if (poly->nterms == 0) {
        return;  // Empty polynomial, nothing to do
    }
    
    // Find the leading term (highest total degree in parameters)
    slong leading_idx = 0;
    slong max_deg = -1;
    
    for (slong i = 0; i < poly->nterms; i++) {
        slong total_deg = 0;
        
        // Calculate total degree in parameters
        if (poly->terms[i].par_exp && poly->npars > 0) {
            for (slong j = 0; j < poly->npars; j++) {
                total_deg += poly->terms[i].par_exp[j];
            }
        }
        
        // Update leading term if this has higher degree
        if (total_deg > max_deg) {
            max_deg = total_deg;
            leading_idx = i;
        }
    }
    
    // Check if already monic
    if (fq_nmod_is_one(poly->terms[leading_idx].coeff, poly->ctx)) {
        return;  // Already monic
    }
    
    // Check if leading coefficient is zero (shouldn't happen for valid polynomials)
    if (fq_nmod_is_zero(poly->terms[leading_idx].coeff, poly->ctx)) {
        printf("Warning: leading coefficient is zero in fq_mvpoly_make_monic\n");
        return;
    }
    
    // Compute inverse of leading coefficient
    fq_nmod_t inv_leading;
    fq_nmod_init(inv_leading, poly->ctx);
    fq_nmod_inv(inv_leading, poly->terms[leading_idx].coeff, poly->ctx);
    
    // Multiply entire polynomial by inverse of leading coefficient
    fq_mvpoly_t normalized;
    fq_mvpoly_scalar_mul(&normalized, poly, inv_leading);
    
    // Replace original polynomial with normalized version
    fq_mvpoly_clear(poly);
    *poly = normalized;
    
    // Cleanup
    fq_nmod_clear(inv_leading, poly->ctx);
}
