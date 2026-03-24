/**
 * fq_multivariate_interpolation.c - Implementation of multivariate polynomial interpolation
 * 
 * Minimally modified fq_nmod interpolation - just add parallelization to core loop
 * Based on the original optimized version
 */

#include "fq_multivariate_interpolation.h"
extern int g_field_equation_reduction;

// Global control for parallelization
int USE_PARALLEL = 1;

static ulong interp_field_size_q(const fq_nmod_ctx_t ctx) {
    mp_limb_t p = fq_nmod_ctx_prime(ctx);
    slong d = fq_nmod_ctx_degree(ctx);
    ulong q = 1;
    for (slong i = 0; i < d; i++) {
        if (q > WORD_MAX / p) return WORD_MAX;
        q *= p;
    }
    return q;
}

static void fq_nmod_poly_reduce_field_equation_interp(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx) {
    if (!g_field_equation_reduction) return;
    slong deg = fq_nmod_poly_degree(poly, ctx);
    if (deg <= 0) return;
    ulong q = interp_field_size_q(ctx);
    if (q <= 1 || q == WORD_MAX) return;
    if ((ulong)deg < q) return;

    fq_nmod_poly_t reduced;
    fq_nmod_poly_init(reduced, ctx);
    fq_nmod_poly_zero(reduced, ctx);
    fq_nmod_t coeff, acc;
    fq_nmod_init(coeff, ctx);
    fq_nmod_init(acc, ctx);

    for (slong i = 0; i <= deg; i++) {
        fq_nmod_poly_get_coeff(coeff, poly, i, ctx);
        if (fq_nmod_is_zero(coeff, ctx)) continue;
        slong tgt = (i == 0 || (ulong)i < q) ? i : (slong)(((ulong)(i - 1) % (q - 1)) + 1);
        fq_nmod_poly_get_coeff(acc, reduced, tgt, ctx);
        fq_nmod_add(acc, acc, coeff, ctx);
        fq_nmod_poly_set_coeff(reduced, tgt, acc, ctx);
    }

    fq_nmod_poly_set(poly, reduced, ctx);
    fq_nmod_clear(coeff, ctx);
    fq_nmod_clear(acc, ctx);
    fq_nmod_poly_clear(reduced, ctx);
}


typedef struct {
    field_elem_u *coeffs;
    slong degree;
} interp_unified_poly_t;

static void interp_build_unified_product_tree(interp_unified_poly_t *result,
                                              const field_elem_u *nodes_arr,
                                              slong start,
                                              slong end,
                                              field_id_t field_id,
                                              void *ctx_ptr) {
    if (start == end) {
        result->degree = 1;
        result->coeffs = (field_elem_u*)malloc(2 * sizeof(field_elem_u));

        field_init_elem(&result->coeffs[0], field_id, ctx_ptr);
        field_init_elem(&result->coeffs[1], field_id, ctx_ptr);

        field_set_one(&result->coeffs[1], field_id, ctx_ptr);
        field_neg(&result->coeffs[0], &nodes_arr[start], field_id, ctx_ptr);
        return;
    }

    slong mid = (start + end) / 2;

    interp_unified_poly_t left, right;
    interp_build_unified_product_tree(&left, nodes_arr, start, mid, field_id, ctx_ptr);
    interp_build_unified_product_tree(&right, nodes_arr, mid + 1, end, field_id, ctx_ptr);

    result->degree = left.degree + right.degree;
    result->coeffs = (field_elem_u*)malloc((result->degree + 1) * sizeof(field_elem_u));

    for (slong i = 0; i <= result->degree; i++) {
        field_init_elem(&result->coeffs[i], field_id, ctx_ptr);
        field_set_zero(&result->coeffs[i], field_id, ctx_ptr);
    }

    field_elem_u temp;
    field_init_elem(&temp, field_id, ctx_ptr);

    for (slong i = 0; i <= left.degree; i++) {
        for (slong j = 0; j <= right.degree; j++) {
            field_mul(&temp, &left.coeffs[i], &right.coeffs[j], field_id, ctx_ptr);
            field_add(&result->coeffs[i + j], &result->coeffs[i + j], &temp, field_id, ctx_ptr);
        }
    }

    field_clear_elem(&temp, field_id, ctx_ptr);

    for (slong i = 0; i <= left.degree; i++) {
        field_clear_elem(&left.coeffs[i], field_id, ctx_ptr);
    }
    for (slong i = 0; i <= right.degree; i++) {
        field_clear_elem(&right.coeffs[i], field_id, ctx_ptr);
    }
    free(left.coeffs);
    free(right.coeffs);
}

static void interp_divide_by_linear(interp_unified_poly_t *quotient,
                                    const interp_unified_poly_t *dividend,
                                    const field_elem_u *node,
                                    field_id_t field_id,
                                    void *ctx_ptr) {
    quotient->degree = dividend->degree - 1;
    quotient->coeffs = (field_elem_u*)malloc((quotient->degree + 1) * sizeof(field_elem_u));

    field_elem_u temp;
    field_init_elem(&temp, field_id, ctx_ptr);

    for (slong i = 0; i <= quotient->degree; i++) {
        field_init_elem(&quotient->coeffs[i], field_id, ctx_ptr);
    }

    field_set_elem(&quotient->coeffs[quotient->degree],
                   &dividend->coeffs[dividend->degree],
                   field_id, ctx_ptr);

    for (slong i = quotient->degree - 1; i >= 0; i--) {
        field_mul(&temp, &quotient->coeffs[i + 1], node, field_id, ctx_ptr);
        field_add(&quotient->coeffs[i], &dividend->coeffs[i + 1], &temp, field_id, ctx_ptr);
    }

    field_clear_elem(&temp, field_id, ctx_ptr);
}

void fq_interpolation_set_parallel(int use_parallel) {
    USE_PARALLEL = use_parallel;
    //printf("Parallelization %s\n", use_parallel ? "enabled" : "disabled");
}

void fq_interpolation_use_half_threads(void) {
    #ifdef _OPENMP
    int max_threads = omp_get_max_threads();
    int half_threads = (max_threads + 1) / 2;
    omp_set_num_threads(half_threads);
    //printf("Using %d threads (half of %d available)\n", half_threads, max_threads);
    #endif
}

// ============= Core optimization: batch matrix evaluation =============

/**
 * Batch evaluation with optional parallelization - minimal change to original
 */
void fq_evaluate_matrix_at_point_batch(fq_nmod_mat_t result_mat,
                                               fq_mvpoly_t **poly_matrix,
                                               slong size,
                                               const fq_nmod_t *var_vals,
                                               const fq_nmod_t *param_vals,
                                               const fq_nmod_ctx_t ctx) {
    // Quick check for empty matrix
    if (size == 0) return;
    
    // Check if we're in GF(2^n) and can use optimized operations
    mp_limb_t prime = fq_nmod_ctx_prime(ctx);
    slong field_degree = fq_nmod_ctx_degree(ctx);
    
    // Initialize unified context
    field_ctx_t unified_ctx;
    field_ctx_init(&unified_ctx, ctx);
    void *ctx_ptr = (unified_ctx.field_id == FIELD_ID_NMOD) ? 
                   (void*)&unified_ctx.ctx.nmod_ctx : 
                   (void*)unified_ctx.ctx.fq_ctx;
    
    //printf("Using %s for matrix evaluation\n", unified_ctx.description);
    
    // Find maximum exponents
    slong max_var_exp = 0, max_par_exp = 0;
    slong max_nvars = 0, max_npars = 0;
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_mvpoly_t *poly = &poly_matrix[i][j];
            
            if (poly->nvars > max_nvars) max_nvars = poly->nvars;
            if (poly->npars > max_npars) max_npars = poly->npars;
            
            for (slong t = 0; t < poly->nterms; t++) {
                if (poly->terms[t].var_exp) {
                    for (slong v = 0; v < poly->nvars; v++) {
                        if (poly->terms[t].var_exp[v] > max_var_exp) {
                            max_var_exp = poly->terms[t].var_exp[v];
                        }
                    }
                }
                if (poly->terms[t].par_exp) {
                    for (slong p = 0; p < poly->npars; p++) {
                        if (poly->terms[t].par_exp[p] > max_par_exp) {
                            max_par_exp = poly->terms[t].par_exp[p];
                        }
                    }
                }
            }
        }
    }
    
    // Precompute powers using unified interface
    field_elem_u **var_powers_unified = NULL;
    field_elem_u **par_powers_unified = NULL;
    
    if (max_nvars > 0 && var_vals && max_var_exp > 0) {
        var_powers_unified = (field_elem_u**)malloc(max_nvars * sizeof(field_elem_u*));
        
        for (slong j = 0; j < max_nvars; j++) {
            var_powers_unified[j] = (field_elem_u*)malloc((max_var_exp + 1) * sizeof(field_elem_u));
            
            // Convert base value
            field_elem_u base;
            field_init_elem(&base, unified_ctx.field_id, ctx_ptr);
            fq_nmod_to_field_elem(&base, var_vals[j], &unified_ctx);
            
            // Compute powers
            field_init_elem(&var_powers_unified[j][0], unified_ctx.field_id, ctx_ptr);
            field_set_one(&var_powers_unified[j][0], unified_ctx.field_id, ctx_ptr);
            
            for (slong k = 1; k <= max_var_exp; k++) {
                field_init_elem(&var_powers_unified[j][k], unified_ctx.field_id, ctx_ptr);
                field_mul(&var_powers_unified[j][k], &var_powers_unified[j][k-1], 
                         &base, unified_ctx.field_id, ctx_ptr);
            }
            
            field_clear_elem(&base, unified_ctx.field_id, ctx_ptr);
        }
    }
    
    if (max_npars > 0 && param_vals && max_par_exp > 0) {
        par_powers_unified = (field_elem_u**)malloc(max_npars * sizeof(field_elem_u*));
        
        for (slong j = 0; j < max_npars; j++) {
            par_powers_unified[j] = (field_elem_u*)malloc((max_par_exp + 1) * sizeof(field_elem_u));
            
            // Convert base value
            field_elem_u base;
            field_init_elem(&base, unified_ctx.field_id, ctx_ptr);
            fq_nmod_to_field_elem(&base, param_vals[j], &unified_ctx);
            
            // Compute powers
            field_init_elem(&par_powers_unified[j][0], unified_ctx.field_id, ctx_ptr);
            field_set_one(&par_powers_unified[j][0], unified_ctx.field_id, ctx_ptr);
            
            for (slong k = 1; k <= max_par_exp; k++) {
                field_init_elem(&par_powers_unified[j][k], unified_ctx.field_id, ctx_ptr);
                field_mul(&par_powers_unified[j][k], &par_powers_unified[j][k-1], 
                         &base, unified_ctx.field_id, ctx_ptr);
            }
            
            field_clear_elem(&base, unified_ctx.field_id, ctx_ptr);
        }
    }
    
    // Process matrix using unified operations
    field_elem_u term_val, result_elem;
    field_init_elem(&term_val, unified_ctx.field_id, ctx_ptr);
    field_init_elem(&result_elem, unified_ctx.field_id, ctx_ptr);
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_mvpoly_t *poly = &poly_matrix[i][j];
            
            field_set_zero(&result_elem, unified_ctx.field_id, ctx_ptr);
            
            // Evaluate polynomial
            for (slong t = 0; t < poly->nterms; t++) {
                // Convert coefficient
                fq_nmod_to_field_elem(&term_val, poly->terms[t].coeff, &unified_ctx);
                
                // Multiply by variable powers
                if (var_powers_unified && poly->terms[t].var_exp) {
                    for (slong v = 0; v < poly->nvars; v++) {
                        slong exp = poly->terms[t].var_exp[v];
                        if (exp > 0) {
                            field_mul(&term_val, &term_val, &var_powers_unified[v][exp], 
                                     unified_ctx.field_id, ctx_ptr);
                        }
                    }
                }
                
                // Multiply by parameter powers
                if (par_powers_unified && poly->terms[t].par_exp) {
                    for (slong p = 0; p < poly->npars; p++) {
                        slong exp = poly->terms[t].par_exp[p];
                        if (exp > 0) {
                            field_mul(&term_val, &term_val, &par_powers_unified[p][exp], 
                                     unified_ctx.field_id, ctx_ptr);
                        }
                    }
                }
                
                // Add to result
                field_add(&result_elem, &result_elem, &term_val, 
                         unified_ctx.field_id, ctx_ptr);
            }
            
            // Convert back to fq_nmod and store
            field_elem_to_fq_nmod(fq_nmod_mat_entry(result_mat, i, j), 
                                 &result_elem, &unified_ctx);
        }
    }
    
    // Cleanup
    field_clear_elem(&term_val, unified_ctx.field_id, ctx_ptr);
    field_clear_elem(&result_elem, unified_ctx.field_id, ctx_ptr);
    
    if (var_powers_unified) {
        for (slong j = 0; j < max_nvars; j++) {
            for (slong k = 0; k <= max_var_exp; k++) {
                field_clear_elem(&var_powers_unified[j][k], unified_ctx.field_id, ctx_ptr);
            }
            free(var_powers_unified[j]);
        }
        free(var_powers_unified);
    }
    
    if (par_powers_unified) {
        for (slong j = 0; j < max_npars; j++) {
            for (slong k = 0; k <= max_par_exp; k++) {
                field_clear_elem(&par_powers_unified[j][k], unified_ctx.field_id, ctx_ptr);
            }
            free(par_powers_unified[j]);
        }
        free(par_powers_unified);
    }
}
// ============= Keep all other functions unchanged =============

// Modified interpolation functions with detailed timing

// Global timing statistics
static InterpolationStats g_stats = {0};

void reset_interpolation_stats(void) {
    memset(&g_stats, 0, sizeof(InterpolationStats));
}

void print_interpolation_stats(void) {
    printf("\n=== Detailed Interpolation Statistics ===\n");
    printf("Total Lagrange interpolation time:     %.3f s (%ld calls, %.3f ms/call)\n", 
           g_stats.lagrange_time, g_stats.lagrange_calls, 
           g_stats.lagrange_calls > 0 ? 1000.0 * g_stats.lagrange_time / g_stats.lagrange_calls : 0);
    printf("Total monomial collection time:        %.3f s\n", g_stats.monomial_collection_time);
    printf("Total monomial interpolation time:     %.3f s\n", g_stats.monomial_interpolation_time);
    printf("Total result construction time:        %.3f s\n", g_stats.result_construction_time);
    printf("Total memory management time:          %.3f s\n", g_stats.memory_time);
    printf("Other time:                           %.3f s\n", g_stats.other_time);
    printf("\nStatistics:\n");
    printf("  Recursive calls:                     %ld\n", g_stats.recursive_calls);
    printf("  Total unique monomials:              %ld\n", g_stats.total_monomials);
    printf("  Total terms processed:               %ld\n", g_stats.total_terms_processed);
    printf("=========================================\n");
}

// Timer utility
double get_time(void) {
    #ifdef _OPENMP
    return omp_get_wtime();
    #else
    return (double)clock() / CLOCKS_PER_SEC;
    #endif
}


// Modified Lagrange interpolation with detailed timing

// Divide-and-conquer approach for building product polynomials
void fq_build_product_tree(fq_nmod_poly_t result, 
                          const fq_nmod_t *nodes, 
                          slong start, slong end,
                          const fq_nmod_ctx_t ctx) {
    if (start == end) {
        // Base case: (x - nodes[start])
        fq_nmod_t neg_node, one;
        fq_nmod_init(neg_node, ctx);
        fq_nmod_init(one, ctx);
        fq_nmod_one(one, ctx);
        fq_nmod_neg(neg_node, nodes[start], ctx);
        
        fq_nmod_poly_zero(result, ctx);
        fq_nmod_poly_set_coeff(result, 1, one, ctx);
        fq_nmod_poly_set_coeff(result, 0, neg_node, ctx);
        
        fq_nmod_clear(neg_node, ctx);
        fq_nmod_clear(one, ctx);
        return;
    }
    
    // Divide
    slong mid = (start + end) / 2;
    
    fq_nmod_poly_t left, right;
    fq_nmod_poly_init(left, ctx);
    fq_nmod_poly_init(right, ctx);
    
    // Conquer
    fq_build_product_tree(left, nodes, start, mid, ctx);
    fq_build_product_tree(right, nodes, mid + 1, end, ctx);
    
    // Combine
    fq_nmod_poly_mul(result, left, right, ctx);
    fq_nmod_poly_reduce_field_equation_interp(result, ctx);
    
    fq_nmod_poly_clear(left, ctx);
    fq_nmod_poly_clear(right, ctx);
}

// Automatically selects the best implementation based on field type

void fq_lagrange_interpolation_optimized(fq_nmod_poly_t result, 
                                             const fq_nmod_t *nodes, 
                                             const fq_nmod_t *values, 
                                             slong k, 
                                             const fq_nmod_ctx_t ctx) {
    double start_time = get_time();
    
    // Track this call
    g_stats.lagrange_calls++;
    
    fq_nmod_poly_zero(result, ctx);
    
    if (k == 0) {
        g_stats.lagrange_time += get_time() - start_time;
        return;
    }
    if (k == 1) {
        fq_nmod_poly_set_coeff(result, 0, values[0], ctx);
        g_stats.lagrange_time += get_time() - start_time;
        return;
    }
    
    // Initialize unified field context for optimal field operations
    field_ctx_t unified_ctx;
    field_ctx_init(&unified_ctx, ctx);
    
    void *ctx_ptr = (unified_ctx.field_id == FIELD_ID_NMOD) ? 
                   (void*)&unified_ctx.ctx.nmod_ctx : 
                   (void*)unified_ctx.ctx.fq_ctx;
    
    //printf("Using optimized %s for Lagrange interpolation\n", unified_ctx.description);
    
    // Convert input nodes and values to unified format
    field_elem_u *unified_nodes = (field_elem_u*)malloc(k * sizeof(field_elem_u));
    field_elem_u *unified_values = (field_elem_u*)malloc(k * sizeof(field_elem_u));
    
    for (slong i = 0; i < k; i++) {
        field_init_elem(&unified_nodes[i], unified_ctx.field_id, ctx_ptr);
        field_init_elem(&unified_values[i], unified_ctx.field_id, ctx_ptr);
        fq_nmod_to_field_elem(&unified_nodes[i], nodes[i], &unified_ctx);
        fq_nmod_to_field_elem(&unified_values[i], values[i], &unified_ctx);
    }
    
    // Step 1: Build the full product polynomial ∏(x - x_i) using divide-and-conquer
    double tree_start = get_time();
    
    // Use unified operations for building product tree
    interp_unified_poly_t full_product;
    interp_build_unified_product_tree(&full_product, unified_nodes, 0, k - 1,
                                      unified_ctx.field_id, ctx_ptr);
    
    double tree_time = get_time() - tree_start;
    
    // Step 2: Compute denominators using unified operations and batch inversion
    double batch_start = get_time();
    field_elem_u *denominators = (field_elem_u*)malloc(k * sizeof(field_elem_u));
    field_elem_u *den_inverses = (field_elem_u*)malloc(k * sizeof(field_elem_u));
    
    for (slong j = 0; j < k; j++) {
        field_init_elem(&denominators[j], unified_ctx.field_id, ctx_ptr);
        field_init_elem(&den_inverses[j], unified_ctx.field_id, ctx_ptr);
        field_set_one(&denominators[j], unified_ctx.field_id, ctx_ptr);
        
        field_elem_u diff;
        field_init_elem(&diff, unified_ctx.field_id, ctx_ptr);
        
        for (slong m = 0; m < k; m++) {
            if (m != j) {
                field_sub(&diff, &unified_nodes[j], &unified_nodes[m], 
                         unified_ctx.field_id, ctx_ptr);
                field_mul(&denominators[j], &denominators[j], &diff, 
                         unified_ctx.field_id, ctx_ptr);
            }
        }
        
        field_clear_elem(&diff, unified_ctx.field_id, ctx_ptr);
    }
    
    // Check for zero denominators
    for (slong j = 0; j < k; j++) {
        if (field_is_zero(&denominators[j], unified_ctx.field_id, ctx_ptr)) {
            printf("ERROR: Zero denominator at index %ld in Lagrange interpolation\n", j);
            // Handle error appropriately
            goto cleanup;
        }
    }
    
    // Batch inversion using unified operations
    field_elem_u *products = (field_elem_u*)malloc(k * sizeof(field_elem_u));
    for (slong i = 0; i < k; i++) {
        field_init_elem(&products[i], unified_ctx.field_id, ctx_ptr);
    }
    
    field_set_elem(&products[0], &denominators[0], unified_ctx.field_id, ctx_ptr);
    for (slong i = 1; i < k; i++) {
        field_mul(&products[i], &products[i-1], &denominators[i], 
                 unified_ctx.field_id, ctx_ptr);
    }
    
    field_elem_u inv;
    field_init_elem(&inv, unified_ctx.field_id, ctx_ptr);
    field_inv(&inv, &products[k-1], unified_ctx.field_id, ctx_ptr);
    
    for (slong i = k - 1; i >= 0; i--) {
        if (i > 0) {
            field_mul(&den_inverses[i], &inv, &products[i-1], 
                     unified_ctx.field_id, ctx_ptr);
            field_mul(&inv, &inv, &denominators[i], unified_ctx.field_id, ctx_ptr);
        } else {
            field_set_elem(&den_inverses[0], &inv, unified_ctx.field_id, ctx_ptr);
        }
    }
    
    field_clear_elem(&inv, unified_ctx.field_id, ctx_ptr);
    for (slong i = 0; i < k; i++) {
        field_clear_elem(&products[i], unified_ctx.field_id, ctx_ptr);
    }
    free(products);
    double batch_time = get_time() - batch_start;
    
    // Step 3: Use polynomial division to get numerators efficiently
    double interp_start = get_time();
    
    // Build result polynomial
    interp_unified_poly_t result_poly;
    result_poly.degree = k - 1;
    result_poly.coeffs = (field_elem_u*)malloc(k * sizeof(field_elem_u));
    
    for (slong i = 0; i < k; i++) {
        field_init_elem(&result_poly.coeffs[i], unified_ctx.field_id, ctx_ptr);
        field_set_zero(&result_poly.coeffs[i], unified_ctx.field_id, ctx_ptr);
    }
    
    field_elem_u coeff, temp_coeff;
    field_init_elem(&coeff, unified_ctx.field_id, ctx_ptr);
    field_init_elem(&temp_coeff, unified_ctx.field_id, ctx_ptr);
    
    for (slong j = 0; j < k; j++) {
        if (field_is_zero(&unified_values[j], unified_ctx.field_id, ctx_ptr)) continue;
        
        // Divide full_product by (x - nodes[j]) to get numerator
        interp_unified_poly_t quotient;
        interp_divide_by_linear(&quotient, &full_product, &unified_nodes[j],
                                unified_ctx.field_id, ctx_ptr);
        
        // Compute coefficient: values[j] * den_inverses[j]
        field_mul(&coeff, &unified_values[j], &den_inverses[j], 
                 unified_ctx.field_id, ctx_ptr);
        
        // Add quotient * coeff to result
        for (slong i = 0; i <= quotient.degree; i++) {
            field_mul(&temp_coeff, &quotient.coeffs[i], &coeff, 
                     unified_ctx.field_id, ctx_ptr);
            field_add(&result_poly.coeffs[i], &result_poly.coeffs[i], &temp_coeff, 
                     unified_ctx.field_id, ctx_ptr);
        }
        
        // Clean up quotient
        for (slong i = 0; i <= quotient.degree; i++) {
            field_clear_elem(&quotient.coeffs[i], unified_ctx.field_id, ctx_ptr);
        }
        free(quotient.coeffs);
    }
    
    double interp_time = get_time() - interp_start;
    
    // Convert result back to fq_nmod_poly_t
    fq_nmod_t temp_fq;
    fq_nmod_init(temp_fq, ctx);
    
    for (slong i = 0; i <= result_poly.degree; i++) {
        field_elem_to_fq_nmod(temp_fq, &result_poly.coeffs[i], &unified_ctx);
        fq_nmod_poly_set_coeff(result, i, temp_fq, ctx);
    }
    
cleanup:
    // Cleanup
    field_clear_elem(&coeff, unified_ctx.field_id, ctx_ptr);
    field_clear_elem(&temp_coeff, unified_ctx.field_id, ctx_ptr);
    fq_nmod_clear(temp_fq, ctx);
    
    // Clean up result polynomial
    for (slong i = 0; i <= result_poly.degree; i++) {
        field_clear_elem(&result_poly.coeffs[i], unified_ctx.field_id, ctx_ptr);
    }
    free(result_poly.coeffs);
    
    // Clean up full product
    for (slong i = 0; i <= full_product.degree; i++) {
        field_clear_elem(&full_product.coeffs[i], unified_ctx.field_id, ctx_ptr);
    }
    free(full_product.coeffs);
    
    // Clean up denominators and inverses
    for (slong j = 0; j < k; j++) {
        field_clear_elem(&denominators[j], unified_ctx.field_id, ctx_ptr);
        field_clear_elem(&den_inverses[j], unified_ctx.field_id, ctx_ptr);
    }
    free(denominators);
    free(den_inverses);
    
    // Clean up unified nodes and values
    for (slong i = 0; i < k; i++) {
        field_clear_elem(&unified_nodes[i], unified_ctx.field_id, ctx_ptr);
        field_clear_elem(&unified_values[i], unified_ctx.field_id, ctx_ptr);
    }
    free(unified_nodes);
    free(unified_values);
    
    // Calculate and update timing
    double total_time = get_time() - start_time;
    g_stats.lagrange_time += total_time;
    
    /*
    printf("\nUnified Lagrange interpolation optimization results:\n");
    printf("  Field type:               %s\n", unified_ctx.description);
    printf("  Product tree construction: %.3f s\n", tree_time);
    printf("  Batch inversion:          %.3f s\n", batch_time);
    printf("  Interpolation:            %.3f s\n", interp_time);
    printf("  Total time:               %.3f s\n", total_time);
    */
}

// Copy tensor interpolation functions unchanged
// Modified version of fq_tensor_interpolation_recursive_optimized with fixed timing

void fq_tensor_interpolation_recursive_optimized(fq_mvpoly_t *result,
                                                slong current_dim,
                                                const fq_nmod_t **grids,
                                                const slong *grid_sizes,
                                                const fq_nmod_t *flat_values,
                                                slong *value_offset,
                                                slong total_dims,
                                                const fq_nmod_ctx_t ctx) {
    double func_start = get_time();
    g_stats.recursive_calls++;
    
    // Local time tracking variables
    double local_mem_time = 0.0;
    double local_lagrange_time = 0.0;  // Track it locally!
    double local_result_time = 0.0;
    double local_monomial_collect_time = 0.0;
    double local_monomial_interp_time = 0.0;
    double local_recursive_time = 0.0;
    
    FQ_INTERP_PRINT("Recursive interpolation: dim=%ld, total=%ld\n", current_dim, total_dims);
    
    // Progress indicator for top-level calls
    static slong top_level_progress = 0;
    static slong top_level_total = 0;
    if (current_dim == total_dims - 1) {
        top_level_total = 1;
        for (slong i = 0; i < total_dims; i++) {
            top_level_total *= grid_sizes[i];
        }
        top_level_progress = 0;
    }
    
    if (current_dim == 0) {
        // Base case
        double base_start = get_time();

        //printf("DEBUG: Base case interpolation with %ld points\n", grid_sizes[0]);
        //printf("DEBUG: Grid points for dimension 0:\n");

        fq_nmod_t *current_values = (fq_nmod_t*) flint_malloc(grid_sizes[0] * sizeof(fq_nmod_t));
        for (slong i = 0; i < grid_sizes[0]; i++) {
            fq_nmod_init(current_values[i], ctx);
            fq_nmod_set(current_values[i], flat_values[(*value_offset)++], ctx);
        }
        
        local_mem_time = get_time() - base_start;
        
        fq_nmod_poly_t uni_result;
        fq_nmod_poly_init(uni_result, ctx);
        
        // Track Lagrange interpolation time locally
        double lagrange_start = get_time();
        
        // Get the current number of Lagrange calls before our call
        slong lagrange_calls_before = g_stats.lagrange_calls;
        double lagrange_time_before = g_stats.lagrange_time;
        
        fq_lagrange_interpolation_optimized(uni_result, grids[0], current_values, 
                                          grid_sizes[0], ctx);
        
        // Calculate how much time was added by this specific call
        local_lagrange_time = g_stats.lagrange_time - lagrange_time_before;
        
        double result_start = get_time();
        fq_mvpoly_init(result, total_dims, 0, ctx);
        
        slong deg = fq_nmod_poly_degree(uni_result, ctx);
        for (slong i = 0; i <= deg; i++) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_nmod_poly_get_coeff(coeff, uni_result, i, ctx);
            
            if (!fq_nmod_is_zero(coeff, ctx)) {
                slong *var_exp = (slong*) flint_calloc(total_dims, sizeof(slong));
                var_exp[0] = i;
                fq_mvpoly_add_term(result, var_exp, NULL, coeff);
                flint_free(var_exp);
            }
            fq_nmod_clear(coeff, ctx);
        }
        local_result_time = get_time() - result_start;
        
        double cleanup_start = get_time();
        for (slong i = 0; i < grid_sizes[0]; i++) {
            fq_nmod_clear(current_values[i], ctx);
        }
        flint_free(current_values);
        fq_nmod_poly_clear(uni_result, ctx);
        local_mem_time += get_time() - cleanup_start;
        
        // Update global stats with LOCAL times (excluding lagrange which is already tracked)
        g_stats.memory_time += local_mem_time;
        g_stats.result_construction_time += local_result_time;
        
        // Calculate other time for this call only
        double total_local_time = get_time() - func_start;
        double tracked_local_time = local_mem_time + local_result_time + local_lagrange_time;
        g_stats.other_time += total_local_time - tracked_local_time;
        
        // Update progress
        if (current_dim == 0 && total_dims > 1) {
            top_level_progress += grid_sizes[0];
            if (top_level_progress % 100 == 0 || top_level_progress == top_level_total) {
                printf("\rInterpolation progress: %ld/%ld", 
                       top_level_progress, top_level_total);
                fflush(stdout);
            }
        }
        
        return;
    }
    
    // Recursive case
    double init_start = get_time();
    fq_mvpoly_init(result, total_dims, 0, ctx);
    
    slong block_size = 1;
    for (slong i = 0; i < current_dim; i++) {
        block_size *= grid_sizes[i];
    }
    
    FQ_INTERP_PRINT("Block size for dim %ld: %ld\n", current_dim, block_size);
    
    fq_mvpoly_t *H_polys = (fq_mvpoly_t*) flint_malloc(grid_sizes[current_dim] * sizeof(fq_mvpoly_t));
    local_mem_time = get_time() - init_start;
    
    // Recursive calls
    double recursive_start = get_time();
    for (slong j = 0; j < grid_sizes[current_dim]; j++) {
        FQ_INTERP_PRINT("Processing grid point %ld in dim %ld\n", j, current_dim);
        fq_tensor_interpolation_recursive_optimized(&H_polys[j], current_dim - 1, grids, grid_sizes,
                                                   flat_values, value_offset, total_dims, ctx);
    }
    local_recursive_time = get_time() - recursive_start;
    FQ_INTERP_PRINT("Collecting unique monomials\n");
    // Collect unique monomials
    double monomial_start = get_time();
    typedef struct {
        slong *exp;
        slong index;
    } monomial_info_t;
    
    monomial_info_t *monomials = NULL;
    slong n_monomials = 0;
    slong alloc_monomials = 0;
    
    for (slong j = 0; j < grid_sizes[current_dim]; j++) {
        g_stats.total_terms_processed += H_polys[j].nterms;
        
        for (slong t = 0; t < H_polys[j].nterms; t++) {
            int found = 0;
            for (slong m = 0; m < n_monomials; m++) {
                int same = 1;
                for (slong k = 0; k < total_dims; k++) {
                    if (k == current_dim) continue;
                    slong exp1 = H_polys[j].terms[t].var_exp ? H_polys[j].terms[t].var_exp[k] : 0;
                    slong exp2 = monomials[m].exp[k];
                    if (exp1 != exp2) {
                        same = 0;
                        break;
                    }
                }
                if (same) {
                    found = 1;
                    break;
                }
            }
            
            if (!found) {
                if (n_monomials >= alloc_monomials) {
                    alloc_monomials = alloc_monomials ? alloc_monomials * 2 : 16;
                    monomials = (monomial_info_t*) realloc(monomials, 
                                                          alloc_monomials * sizeof(monomial_info_t));
                }
                
                monomials[n_monomials].exp = (slong*) flint_calloc(total_dims, sizeof(slong));
                if (H_polys[j].terms[t].var_exp) {
                    for (slong k = 0; k < total_dims; k++) {
                        if (k != current_dim) {
                            monomials[n_monomials].exp[k] = H_polys[j].terms[t].var_exp[k];
                        }
                    }
                }
                monomials[n_monomials].index = n_monomials;
                n_monomials++;
            }
        }
    }
    
    g_stats.total_monomials += n_monomials;
    local_monomial_collect_time = get_time() - monomial_start;
    
    FQ_INTERP_PRINT("Found %ld unique monomial patterns in dim %ld\n", n_monomials, current_dim);
    
    // Interpolate each monomial
    double interp_start = get_time();
    fq_nmod_t *interp_values = (fq_nmod_t*) flint_malloc(grid_sizes[current_dim] * sizeof(fq_nmod_t));
    for (slong i = 0; i < grid_sizes[current_dim]; i++) {
        fq_nmod_init(interp_values[i], ctx);
    }
    
    // Track lagrange time before the loop
    double lagrange_time_before_loop = g_stats.lagrange_time;
    
    field_ctx_t unified_ctx;
    field_ctx_init(&unified_ctx, ctx);
    void *ctx_ptr = (unified_ctx.field_id == FIELD_ID_NMOD) ? 
                   (void*)&unified_ctx.ctx.nmod_ctx : 
                   (void*)unified_ctx.ctx.fq_ctx;
    
    typedef struct term_node {
        slong term_index;
        struct term_node *next;
    } term_node_t;
    
    term_node_t ***monomial_lists = (term_node_t***)malloc(grid_sizes[current_dim] * sizeof(term_node_t**));
    
    for (slong j = 0; j < grid_sizes[current_dim]; j++) {
        monomial_lists[j] = (term_node_t**)calloc(n_monomials, sizeof(term_node_t*));
        
        for (slong t = 0; t < H_polys[j].nterms; t++) {
            for (slong m = 0; m < n_monomials; m++) {
                int same = 1;
                for (slong k = 0; k < total_dims; k++) {
                    if (k == current_dim) continue;
                    slong exp1 = H_polys[j].terms[t].var_exp ? H_polys[j].terms[t].var_exp[k] : 0;
                    slong exp2 = monomials[m].exp[k];
                    if (exp1 != exp2) {
                        same = 0;
                        break;
                    }
                }
                if (same) {
                    term_node_t *node = (term_node_t*)malloc(sizeof(term_node_t));
                    node->term_index = t;
                    node->next = monomial_lists[j][m];
                    monomial_lists[j][m] = node;
                }
            }
        }
    }
    
    slong *final_exp = (slong*)flint_calloc(total_dims, sizeof(slong));
    field_elem_u *unified_interp_values = (field_elem_u*)malloc(grid_sizes[current_dim] * sizeof(field_elem_u));
    field_elem_u sum_elem, term_elem, result_coeff;
    
    for (slong j = 0; j < grid_sizes[current_dim]; j++) {
        field_init_elem(&unified_interp_values[j], unified_ctx.field_id, ctx_ptr);
    }
    field_init_elem(&sum_elem, unified_ctx.field_id, ctx_ptr);
    field_init_elem(&term_elem, unified_ctx.field_id, ctx_ptr);
    field_init_elem(&result_coeff, unified_ctx.field_id, ctx_ptr);
    
    for (slong m = 0; m < n_monomials; m++) {
        /*
        if (n_monomials > 100 && (m % (n_monomials / 10) == 0 || m == n_monomials - 1)) {
            printf("\rMonomial interpolation progress: %ld/%ld (%.1f%%)", 
                   m + 1, n_monomials, 100.0 * (m + 1) / n_monomials);
            //fflush(stdout);
        }
        */
        for (slong j = 0; j < grid_sizes[current_dim]; j++) {
            field_set_zero(&sum_elem, unified_ctx.field_id, ctx_ptr);
            
            term_node_t *node = monomial_lists[j][m];
            while (node) {
                slong t = node->term_index;
                fq_nmod_to_field_elem(&term_elem, H_polys[j].terms[t].coeff, &unified_ctx);
                field_add(&sum_elem, &sum_elem, &term_elem, unified_ctx.field_id, ctx_ptr);
                node = node->next;
            }
            
            field_set_elem(&unified_interp_values[j], &sum_elem, unified_ctx.field_id, ctx_ptr);
            
            field_elem_to_fq_nmod(interp_values[j], &unified_interp_values[j], &unified_ctx);
        }
        
        fq_nmod_poly_t coeff_poly;
        fq_nmod_poly_init(coeff_poly, ctx);
        fq_lagrange_interpolation_optimized(coeff_poly, grids[current_dim], interp_values, 
                                          grid_sizes[current_dim], ctx);
        
        slong deg = fq_nmod_poly_degree(coeff_poly, ctx);
        for (slong d = 0; d <= deg; d++) {
            fq_nmod_t temp_coeff;
            fq_nmod_init(temp_coeff, ctx);
            fq_nmod_poly_get_coeff(temp_coeff, coeff_poly, d, ctx);
            
            fq_nmod_to_field_elem(&result_coeff, temp_coeff, &unified_ctx);
            
            if (!field_is_zero(&result_coeff, unified_ctx.field_id, ctx_ptr)) {
                memset(final_exp, 0, total_dims * sizeof(slong));
                final_exp[current_dim] = d;
                for (slong k = 0; k < total_dims; k++) {
                    if (k != current_dim) {
                        final_exp[k] = monomials[m].exp[k];
                    }
                }
                fq_mvpoly_add_term_fast(result, final_exp, NULL, temp_coeff);
            }
            fq_nmod_clear(temp_coeff, ctx);
        }
        fq_nmod_poly_clear(coeff_poly, ctx);
    }
    
    for (slong j = 0; j < grid_sizes[current_dim]; j++) {
        field_clear_elem(&unified_interp_values[j], unified_ctx.field_id, ctx_ptr);
    }
    field_clear_elem(&sum_elem, unified_ctx.field_id, ctx_ptr);
    field_clear_elem(&term_elem, unified_ctx.field_id, ctx_ptr);
    field_clear_elem(&result_coeff, unified_ctx.field_id, ctx_ptr);
    
    for (slong j = 0; j < grid_sizes[current_dim]; j++) {
        for (slong m = 0; m < n_monomials; m++) {
            term_node_t *node = monomial_lists[j][m];
            while (node) {
                term_node_t *next = node->next;
                free(node);
                node = next;
            }
        }
        free(monomial_lists[j]);
    }
    free(monomial_lists);
    free(unified_interp_values);
    flint_free(final_exp);
    
    // Calculate lagrange time used in this loop
    local_lagrange_time = g_stats.lagrange_time - lagrange_time_before_loop;
    local_monomial_interp_time = get_time() - interp_start - local_lagrange_time;
    
    // Cleanup
    double cleanup_start = get_time();
    for (slong i = 0; i < grid_sizes[current_dim]; i++) {
        fq_nmod_clear(interp_values[i], ctx);
        fq_mvpoly_clear(&H_polys[i]);
    }
    flint_free(interp_values);
    flint_free(H_polys);
    
    for (slong i = 0; i < n_monomials; i++) {
        flint_free(monomials[i].exp);
    }
    if (monomials) free(monomials);
    local_mem_time += get_time() - cleanup_start;
    
    FQ_INTERP_PRINT("Completed dim %ld with %ld result terms\n", current_dim, result->nterms);
    
    // Update global stats with LOCAL times
    g_stats.memory_time += local_mem_time;
    g_stats.monomial_collection_time += local_monomial_collect_time;
    g_stats.monomial_interpolation_time += local_monomial_interp_time;
    
    // Calculate other time for this call - now includes all tracked time
    double total_local_time = get_time() - func_start;
    double tracked_local_time = local_mem_time + local_recursive_time + 
                               local_monomial_collect_time + local_monomial_interp_time + 
                               local_lagrange_time;
    g_stats.other_time += total_local_time - tracked_local_time;
}

void fq_tensor_interpolation_all_vars_optimized(fq_mvpoly_t *result,
                                               const fq_nmod_t **grids,
                                               const fq_nmod_t *values,
                                               const slong *grid_sizes,
                                               slong nvars,
                                               slong npars,
                                               const fq_nmod_ctx_t ctx) {
    slong total_dims = nvars + npars;
    
    FQ_INTERP_PRINT("Starting tensor interpolation: nvars=%ld, npars=%ld, total=%ld\n", 
                    nvars, npars, total_dims);
    
    // Reset statistics
    reset_interpolation_stats();
    
    if (total_dims == 0) {
        fq_mvpoly_init(result, 0, 0, ctx);
        if (!fq_nmod_is_zero(values[0], ctx)) {
            fq_mvpoly_add_term(result, NULL, NULL, values[0]);
        }
        return;
    }
    
    slong value_offset = 0;
    double start_time = get_time();
    
    fq_tensor_interpolation_recursive_optimized(result, total_dims - 1, grids, grid_sizes,
                                               values, &value_offset, total_dims, ctx);
    
    printf("\n");
    
    if (result->nvars != nvars || result->npars != npars) {
        FQ_INTERP_PRINT("Fixing result structure: %ld->%ld vars, %ld->%ld pars\n",
                        result->nvars, nvars, result->npars, npars);
        
        double fix_start = get_time();
        fq_mvpoly_t fixed_result;
        fq_mvpoly_init(&fixed_result, nvars, npars, ctx);
        
        for (slong t = 0; t < result->nterms; t++) {
            slong *var_exp = NULL;
            slong *par_exp = NULL;
            
            if (nvars > 0) {
                var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
                if (result->terms[t].var_exp) {
                    for (slong i = 0; i < FLINT_MIN(nvars, total_dims); i++) {
                        var_exp[i] = result->terms[t].var_exp[i];
                    }
                }
            }
            
            if (npars > 0) {
                par_exp = (slong*) flint_calloc(npars, sizeof(slong));
                if (result->terms[t].var_exp && nvars < total_dims) {
                    for (slong i = 0; i < FLINT_MIN(npars, total_dims - nvars); i++) {
                        par_exp[i] = result->terms[t].var_exp[nvars + i];
                    }
                }
            }
            
            fq_mvpoly_add_term_fast(&fixed_result, var_exp, par_exp, result->terms[t].coeff);
            
            if (var_exp) flint_free(var_exp);
            if (par_exp) flint_free(par_exp);
        }
        
        fq_mvpoly_clear(result);
        *result = fixed_result;
        g_stats.result_construction_time += get_time() - fix_start;
    }
    //printf("Interpolation Over\n");
    // Print detailed statistics
    // print_interpolation_stats();
}

void fq_generate_evaluation_points_optimized(fq_nmod_t **grids, slong *grid_sizes, 
                                            slong total_vars, slong *degrees, 
                                            const fq_nmod_ctx_t ctx) {
    slong extra_points = 1;
    
    slong field_degree = fq_nmod_ctx_degree(ctx);
    mp_limb_t prime = fq_nmod_ctx_prime(ctx);
    
    // Calculate field size
    slong field_size;
    if (field_degree == 1) {
        field_size = prime;
    } else {
        // For extension fields, field_size = prime^field_degree
        // Be careful with overflow for large fields
        field_size = 1;
        for (slong i = 0; i < field_degree; i++) {
            if (field_size > WORD_MAX / prime) {
                // Field is too large to represent as slong
                field_size = WORD_MAX;
                break;
            }
            field_size *= prime;
        }
    }
    
    FQ_INTERP_PRINT("Field size: %ld (prime=%lu, degree=%ld)\n", 
                    field_size, prime, field_degree);
    
    // Initialize unified field context for ALL field types
    field_ctx_t unified_ctx;
    field_ctx_init(&unified_ctx, ctx);
    
    // Initialize unified elements
    void *ctx_ptr = (unified_ctx.field_id == FIELD_ID_NMOD) ? 
                   (void*)&unified_ctx.ctx.nmod_ctx : 
                   (void*)unified_ctx.ctx.fq_ctx;
    
    field_elem_u generator_unified, temp_unified, one_unified;
    field_init_elem(&generator_unified, unified_ctx.field_id, ctx_ptr);
    field_init_elem(&temp_unified, unified_ctx.field_id, ctx_ptr);
    field_init_elem(&one_unified, unified_ctx.field_id, ctx_ptr);
    field_set_one(&one_unified, unified_ctx.field_id, ctx_ptr);
    
    // Regular FLINT elements for compatibility
    fq_nmod_t generator, temp, one;
    fq_nmod_init(generator, ctx);
    fq_nmod_init(temp, ctx);
    fq_nmod_init(one, ctx);
    fq_nmod_one(one, ctx);
    
    fq_nmod_gen(generator, ctx);
    // Convert generator to unified format
    fq_nmod_to_field_elem(&generator_unified, generator, &unified_ctx);
    
    FQ_INTERP_PRINT("Field generator = ");
    
    // Generate evaluation points for each variable
    for (slong i = 0; i < total_vars; i++) {
        slong requested_degree = degrees[i];
        slong effective_degree = requested_degree;
        if (effective_degree < 0) {
            effective_degree = 0;
        }
        if (g_field_equation_reduction && field_size > 1 && effective_degree >= field_size) {
            effective_degree = field_size - 1;
        }
        grid_sizes[i] = effective_degree + extra_points;
        
        // IMPORTANT: Cap grid size at field size
        if (grid_sizes[i] > field_size && field_size > 0) {
            printf("Variable %ld: capping grid size from %ld to field size %ld\n", 
                   i, grid_sizes[i], field_size);
            grid_sizes[i] = field_size;
        }
        
        grids[i] = (fq_nmod_t*) flint_malloc(grid_sizes[i] * sizeof(fq_nmod_t));
        
        FQ_INTERP_PRINT("Generating grid for variable %ld, requested_degree=%ld, effective_degree=%ld, grid_size=%ld\n", 
                       i, requested_degree, effective_degree, grid_sizes[i]);
        
        // Initialize all grid points
        for (slong j = 0; j < grid_sizes[i]; j++) {
            fq_nmod_init(grids[i][j], ctx);
        }
        
        if (effective_degree == 0 && extra_points == 1) {
            // Special case: degree 0 only needs one point
            FQ_INTERP_PRINT("Case: degree 0, using single point\n");
            fq_nmod_set(grids[i][0], one, ctx);
        } else {
            FQ_INTERP_PRINT("Generating %ld distinct points for variable %ld\n", grid_sizes[i], i);
            
            slong points_generated = 0;
            
            // ========== PRIME FIELD FAST PATH ==========
            // For prime fields GF(p), use consecutive integers 0, 1, 2, ..., n-1
            // This is simple, efficient, and guaranteed to give distinct points
            if (field_degree == 1 && grid_sizes[i] <= prime) {
                FQ_INTERP_PRINT("  Prime field detected, using consecutive integers\n");
                
                for (slong j = 0; j < grid_sizes[i]; j++) {
                    // Set grid[i][j] = j (as field element)
                    fq_nmod_set_ui(grids[i][j], j, ctx);
                }
                points_generated = grid_sizes[i];
                
                FQ_INTERP_PRINT("  Generated %ld points: 0, 1, 2, ..., %ld\n", 
                               points_generated, points_generated - 1);
                
            } else {
                // ========== EXTENSION FIELD OR COMPLEX CASE ==========
                // For extension fields or when we need more sophisticated point generation
                
                // Point 1: Always use 0
                fq_nmod_zero(grids[i][0], ctx);
                points_generated = 1;
                
                if (grid_sizes[i] > 1) {
                    // Point 2: Always use 1
                    fq_nmod_one(grids[i][1], ctx);
                    points_generated = 2;
                    
                    if (grid_sizes[i] > 2) {
                        // Phase 1: Try powers of generator (with duplicate detection)
                        field_elem_u current_power;
                        field_init_elem(&current_power, unified_ctx.field_id, ctx_ptr);
                        field_set_elem(&current_power, &generator_unified, unified_ctx.field_id, ctx_ptr);
                        
                        int cycle_detected = 0;
                        for (slong j = 2; j < grid_sizes[i] && !cycle_detected; j++) {
                            // Convert current power to fq_nmod
                            field_elem_to_fq_nmod(grids[i][j], &current_power, &unified_ctx);
                            
                            // Check for duplicates with all previous points
                            int is_duplicate = 0;
                            for (slong k = 0; k < points_generated; k++) {
                                if (fq_nmod_equal(grids[i][j], grids[i][k], ctx)) {
                                    is_duplicate = 1;
                                    cycle_detected = 1;  // Generator cycle completed
                                    break;
                                }
                            }
                            
                            if (!is_duplicate) {
                                points_generated++;
                                // Compute next power: current = current * generator
                                field_mul(&current_power, &current_power, &generator_unified, 
                                         unified_ctx.field_id, ctx_ptr);
                            } else {
                                break;
                            }
                        }
                        
                        field_clear_elem(&current_power, unified_ctx.field_id, ctx_ptr);
                        
                        FQ_INTERP_PRINT("  Generated %ld points from powers of generator\n", points_generated);
                        
                        // Phase 2: If we still need more points, use linear combinations
                        if (points_generated < grid_sizes[i]) {
                            FQ_INTERP_PRINT("  Need %ld more points, trying linear combinations...\n", 
                                           grid_sizes[i] - points_generated);
                            
                            field_elem_u combo, gen_mult, const_part;
                            field_init_elem(&combo, unified_ctx.field_id, ctx_ptr);
                            field_init_elem(&gen_mult, unified_ctx.field_id, ctx_ptr);
                            field_init_elem(&const_part, unified_ctx.field_id, ctx_ptr);
                            
                            // Try combinations: a*g + b for small a, b
                            int found_enough = 0;
                            for (slong a = 1; a <= 10 && !found_enough; a++) {
                                for (slong b = 2; b <= 10 && !found_enough; b++) {  // Start b at 2 (already have 0,1)
                                    // Compute a*generator
                                    field_set_zero(&gen_mult, unified_ctx.field_id, ctx_ptr);
                                    for (slong m = 0; m < a; m++) {
                                        field_add(&gen_mult, &gen_mult, &generator_unified,
                                                 unified_ctx.field_id, ctx_ptr);
                                    }
                                    
                                    // Compute b*1
                                    field_set_zero(&const_part, unified_ctx.field_id, ctx_ptr);
                                    for (slong m = 0; m < b; m++) {
                                        field_add(&const_part, &const_part, &one_unified,
                                                 unified_ctx.field_id, ctx_ptr);
                                    }
                                    
                                    // combo = a*g + b
                                    field_add(&combo, &gen_mult, &const_part,
                                             unified_ctx.field_id, ctx_ptr);
                                    
                                    // Convert to fq_nmod
                                    fq_nmod_t candidate;
                                    fq_nmod_init(candidate, ctx);
                                    field_elem_to_fq_nmod(candidate, &combo, &unified_ctx);
                                    
                                    // Check if unique
                                    int is_duplicate = 0;
                                    for (slong k = 0; k < points_generated; k++) {
                                        if (fq_nmod_equal(candidate, grids[i][k], ctx)) {
                                            is_duplicate = 1;
                                            break;
                                        }
                                    }
                                    
                                    if (!is_duplicate) {
                                        fq_nmod_set(grids[i][points_generated], candidate, ctx);
                                        points_generated++;
                                        
                                        if (points_generated >= grid_sizes[i]) {
                                            found_enough = 1;
                                        }
                                    }
                                    
                                    fq_nmod_clear(candidate, ctx);
                                }
                            }
                            
                            field_clear_elem(&combo, unified_ctx.field_id, ctx_ptr);
                            field_clear_elem(&gen_mult, unified_ctx.field_id, ctx_ptr);
                            field_clear_elem(&const_part, unified_ctx.field_id, ctx_ptr);
                            
                            FQ_INTERP_PRINT("  After linear combinations: %ld total points\n", points_generated);
                        }
                    }
                }
            }
            
            // ========== CRITICAL CHECK ==========
            // Verify we generated enough distinct points
            if (points_generated < grid_sizes[i]) {
                printf("\n!!! CRITICAL WARNING !!!\n");
                printf("Field GF(%lu^%ld) is too small!\n", prime, field_degree);
                printf("  Field size: %ld elements\n", field_size);
                printf("  Requested: %ld distinct points for variable %ld\n", 
                       grid_sizes[i], i);
                printf("  Generated: %ld distinct points\n", points_generated);
                printf("  Degree bound: %ld (needs %ld+1 points)\n", 
                       effective_degree, effective_degree);
                
                if (points_generated < effective_degree + 1) {
                    printf("\n!!! FATAL ERROR !!!\n");
                    printf("Cannot perform interpolation: need at least %ld points for degree %ld\n",
                           effective_degree + 1, effective_degree);
                    printf("Possible solutions:\n");
                    printf("  1. Use a larger field (e.g., increase field degree)\n");
                    printf("  2. Reduce degree bounds\n");
                    printf("  3. Use modular methods or different interpolation approach\n");
                    
                    // Cleanup before returning
                    for (slong j = 0; j < i; j++) {
                        for (slong k = 0; k < grid_sizes[j]; k++) {
                            fq_nmod_clear(grids[j][k], ctx);
                        }
                        free(grids[j]);
                    }
                    for (slong j = 0; j < grid_sizes[i]; j++) {
                        fq_nmod_clear(grids[i][j], ctx);
                    }
                    free(grids[i]);
                    
                    fq_nmod_clear(generator, ctx);
                    fq_nmod_clear(temp, ctx);
                    fq_nmod_clear(one, ctx);
                    field_clear_elem(&generator_unified, unified_ctx.field_id, ctx_ptr);
                    field_clear_elem(&temp_unified, unified_ctx.field_id, ctx_ptr);
                    field_clear_elem(&one_unified, unified_ctx.field_id, ctx_ptr);
                    
                    return;  // Fatal error - cannot continue
                }
                
                printf("Adjusting grid size from %ld to %ld for variable %ld\n",
                       grid_sizes[i], points_generated, i);
                grid_sizes[i] = points_generated;
            }
            
            // ========== FINAL VERIFICATION ==========
            // Verify all points are actually distinct (sanity check)
            FQ_INTERP_PRINT("  Verifying all %ld points are distinct...\n", grid_sizes[i]);
            for (slong j = 0; j < grid_sizes[i]; j++) {
                for (slong k = j + 1; k < grid_sizes[i]; k++) {
                    if (fq_nmod_equal(grids[i][j], grids[i][k], ctx)) {
                        printf("\n!!! FATAL BUG !!!\n");
                        printf("Duplicate points detected after generation!\n");
                        printf("  Variable %ld: point[%ld] == point[%ld]\n", i, j, k);
                        printf("  Point value: ");
                        fq_nmod_print_pretty(grids[i][j], ctx);
                        printf("\n");
                        printf("This is a bug in the point generation code!\n");
                        
                        // Cleanup and abort
                        for (slong v = 0; v <= i; v++) {
                            for (slong p = 0; p < grid_sizes[v]; p++) {
                                fq_nmod_clear(grids[v][p], ctx);
                            }
                            free(grids[v]);
                        }
                        
                        fq_nmod_clear(generator, ctx);
                        fq_nmod_clear(temp, ctx);
                        fq_nmod_clear(one, ctx);
                        field_clear_elem(&generator_unified, unified_ctx.field_id, ctx_ptr);
                        field_clear_elem(&temp_unified, unified_ctx.field_id, ctx_ptr);
                        field_clear_elem(&one_unified, unified_ctx.field_id, ctx_ptr);
                        
                        return;  // Fatal bug
                    }
                }
            }
            FQ_INTERP_PRINT("  ✓ All points verified distinct\n");
        }
        
        // Debug: print first few points
        FQ_INTERP_PRINT("Grid[%ld]: ", i);
    }
    
    // Cleanup
    fq_nmod_clear(generator, ctx);
    fq_nmod_clear(temp, ctx);
    fq_nmod_clear(one, ctx);
    
    field_clear_elem(&generator_unified, unified_ctx.field_id, ctx_ptr);
    field_clear_elem(&temp_unified, unified_ctx.field_id, ctx_ptr);
    field_clear_elem(&one_unified, unified_ctx.field_id, ctx_ptr);
}
    

void fq_compute_det_degree_bounds_optimized(slong *bounds, fq_mvpoly_t **matrix, 
                                           slong size, slong total_vars) {
    FQ_INTERP_PRINT("Computing degree bounds for %ld variables\n", total_vars);
    
    for (slong var = 0; var < total_vars; var++) {
        bounds[var] = 0;
        
        for (slong row = 0; row < size; row++) {
            slong row_max_deg = 0;
            
            for (slong col = 0; col < size; col++) {
                for (slong t = 0; t < matrix[row][col].nterms; t++) {
                    slong deg = 0;
                    if (matrix[row][col].terms[t].var_exp && var < matrix[row][col].nvars) {
                        deg = matrix[row][col].terms[t].var_exp[var];
                    } else if (matrix[row][col].terms[t].par_exp && 
                              var >= matrix[row][col].nvars && 
                              var - matrix[row][col].nvars < matrix[row][col].npars) {
                        deg = matrix[row][col].terms[t].par_exp[var - matrix[row][col].nvars];
                    }
                    
                    if (deg > row_max_deg) {
                        row_max_deg = deg;
                    }
                }
            }
            bounds[var] += row_max_deg;
        }
        
        if (bounds[var] < 0) bounds[var] = 0;
        
        FQ_INTERP_PRINT("Variable %ld degree bound: %ld\n", var, bounds[var]);
    }
}

// Main interpolation function with PARALLELIZATION ON POINTS
// Main interpolation function with PARALLELIZATION ON POINTS
void fq_compute_det_by_interpolation_optimized(fq_mvpoly_t *result,
                                              fq_mvpoly_t **matrix,
                                              slong size,
                                              slong nvars,
                                              slong npars,
                                              const fq_nmod_ctx_t ctx,
                                              slong *degree_bounds) {
    FQ_INTERP_PRINT("\n=== Optimized FQ Determinant Interpolation ===\n");
    #ifdef _OPENMP
    printf("OpenMP available: %d threads max\n", omp_get_max_threads());
    //printf("Parallelization: %s\n", USE_PARALLEL ? "enabled" : "disabled");
    #else
    printf("OpenMP not available - running sequential\n");
    #endif
    
    // Time tracking - FIX: Use consistent timing method
    #ifdef _OPENMP
    double start_total_wtime = omp_get_wtime();  // Add this declaration!
    #else
    clock_t start_total = clock();
    #endif
    
    double time_matrix_eval = 0.0;
    double time_det_computation = 0.0;
    double time_interpolation = 0.0;
    double time_preprocessing = 0.0;
    
    // Setup phase (unchanged)
    #ifdef _OPENMP
    double temp_time = omp_get_wtime();
    #else
    clock_t temp_time = clock();
    #endif
    
    slong actual_nvars = 0;
    slong actual_npars = npars;
    
    if (size > 0) {
        actual_nvars = matrix[0][0].nvars;
        actual_npars = matrix[0][0].npars;
        FQ_INTERP_PRINT("Matrix structure: nvars=%ld, npars=%ld\n", actual_nvars, actual_npars);
    }
    
    slong total_vars = actual_nvars + actual_npars;
    
    if (total_vars == 0) {
        // Constant case (unchanged)
        fq_nmod_mat_t const_mat;
        fq_nmod_mat_init(const_mat, size, size, ctx);
        
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                if (matrix[i][j].nterms > 0) {
                    fq_nmod_set(fq_nmod_mat_entry(const_mat, i, j), 
                               matrix[i][j].terms[0].coeff, ctx);
                } else {
                    fq_nmod_zero(fq_nmod_mat_entry(const_mat, i, j), ctx);
                }
            }
        }
        
        fq_nmod_t det;
        fq_nmod_init(det, ctx);
        fq_nmod_mat_det(det, const_mat, ctx);
        
        fq_mvpoly_init(result, actual_nvars, actual_npars, ctx);
        if (!fq_nmod_is_zero(det, ctx)) {
            fq_mvpoly_add_term(result, NULL, NULL, det);
        }
        
        fq_nmod_clear(det, ctx);
        fq_nmod_mat_clear(const_mat, ctx);
        return;
    }
    
    // Compute degree bounds
    slong *computed_bounds = NULL;
    if (!degree_bounds) {
        computed_bounds = (slong*) malloc(total_vars * sizeof(slong));
        fq_compute_det_degree_bounds_optimized(computed_bounds, matrix, size, total_vars);
        degree_bounds = computed_bounds;
    }
    
    // Generate interpolation grids
    fq_nmod_t **grids = (fq_nmod_t**) malloc(total_vars * sizeof(fq_nmod_t*));
    slong *grid_sizes = (slong*) malloc(total_vars * sizeof(slong));
    
    fq_generate_evaluation_points_optimized(grids, grid_sizes, total_vars, degree_bounds, ctx);
    
    // Calculate total points
    slong total_points = 1;
    for (slong i = 0; i < total_vars; i++) {
        total_points *= grid_sizes[i];
    }
    
    FQ_INTERP_PRINT("Total interpolation points: %ld\n", total_points);
    
    #ifdef _OPENMP
    time_preprocessing = omp_get_wtime() - temp_time;
    #else
    time_preprocessing = (double)(clock() - temp_time) / CLOCKS_PER_SEC;
    #endif
    
    if (total_points > 100000) {
        printf("Warning: Very large interpolation problem (%ld points)\n", total_points);
    }
    
    // Allocate values array
    fq_nmod_t *values = (fq_nmod_t*) malloc(total_points * sizeof(fq_nmod_t));
    for (slong i = 0; i < total_points; i++) {
        fq_nmod_init(values[i], ctx);
    }
    
    FQ_INTERP_PRINT("Computing interpolation points...\n");
    
    // *** MAIN CHANGE: Parallelize over interpolation points ***
    #ifdef _OPENMP
    if (USE_PARALLEL && total_points > 10) {
        // Parallel version
        double wall_eval_time = 0.0;
        double wall_det_time = 0.0;
        double cpu_eval_time = 0.0;
        double cpu_det_time = 0.0;
        slong par_eval_count = 0;
        slong par_det_count = 0;
        
        double parallel_start = omp_get_wtime();
        
        #pragma omp parallel
        {
            // Thread-local evaluation matrix
            fq_nmod_mat_t eval_mat;
            fq_nmod_mat_init(eval_mat, size, size, ctx);
            
            // Thread-local timing
            double thread_eval_time = 0.0;
            double thread_det_time = 0.0;
            slong thread_eval_count = 0;
            slong thread_det_count = 0;
            
            #pragma omp for schedule(dynamic, 1)
            for (slong point = 0; point < total_points; point++) {
                // Calculate indices for this point
                slong temp = point;
                slong *indices = (slong*) calloc(total_vars, sizeof(slong));
                for (slong i = 0; i < total_vars; i++) {
                    indices[i] = temp % grid_sizes[i];
                    temp /= grid_sizes[i];
                }
                
                // Calculate current point coordinates
                fq_nmod_t *coords = (fq_nmod_t*) malloc(total_vars * sizeof(fq_nmod_t));
                for (slong i = 0; i < total_vars; i++) {
                    fq_nmod_init(coords[i], ctx);
                    fq_nmod_set(coords[i], grids[i][indices[i]], ctx);
                }
                
                // Evaluate matrix at this point
                double matrix_eval_start = omp_get_wtime();
                
                fq_nmod_t *var_vals = NULL;
                fq_nmod_t *param_vals = NULL;
                
                if (actual_nvars > 0) {
                    var_vals = coords;
                }
                if (actual_npars > 0) {
                    param_vals = coords + actual_nvars;
                }
                
                // Use the NON-PARALLEL batch evaluation
                fq_evaluate_matrix_at_point_batch(eval_mat, matrix, size, 
                                                 var_vals, param_vals, ctx);
                
                thread_eval_time += omp_get_wtime() - matrix_eval_start;
                thread_eval_count += size * size;
                
                // Compute determinant
                double matrix_det_start = omp_get_wtime();
                fq_nmod_mat_det(values[point], eval_mat, ctx);
                thread_det_time += omp_get_wtime() - matrix_det_start;
                thread_det_count++;
                
                // Cleanup
                for (slong i = 0; i < total_vars; i++) {
                    fq_nmod_clear(coords[i], ctx);
                }
                free(coords);
                free(indices);
                
                // Progress report (only from thread 0)
                if (omp_get_thread_num() == 0 && (point % 100 == 0)) {
                    // Use atomic read to get a better estimate
                    slong completed = point;  // This is just an approximation
                    printf("\rProgress: %ld/%ld", completed, total_points);
                    fflush(stdout);
                }
            }
            
            // Accumulate thread timings
            #pragma omp critical
            {
                cpu_eval_time  += thread_eval_time;
                cpu_det_time += thread_det_time;
                par_eval_count += thread_eval_count;
                par_det_count += thread_det_count;
            }
            
            fq_nmod_mat_clear(eval_mat, ctx);
            
        }
        double parallel_end = omp_get_wtime();
        wall_eval_time = (parallel_end - parallel_start) * 
                         (cpu_eval_time / (cpu_eval_time + cpu_det_time));
        wall_det_time = (parallel_end - parallel_start) * 
                        (cpu_det_time / (cpu_eval_time + cpu_det_time));
        
        time_matrix_eval = wall_eval_time;
        time_det_computation = wall_det_time;

        // Store CPU times for detailed report
        double cpu_time_matrix_eval = cpu_eval_time;
        double cpu_time_det_computation = cpu_det_time;
        /*
        printf("\rProgress: 100.0%%\n");

        // Print both wall time and CPU time statistics
        printf("\n=== Parallel Execution Statistics ===\n");
        printf("Wall time (real time):\n");
        printf("  Matrix evaluation:    %.3f s\n", wall_eval_time);
        printf("  Determinant comp:     %.3f s\n", wall_det_time);
        printf("  Total parallel:       %.3f s\n", parallel_end - parallel_start);
        printf("\nCPU time (sum across threads):\n");
        printf("  Matrix evaluation:    %.3f s\n", cpu_time_matrix_eval);
        printf("  Determinant comp:     %.3f s\n", cpu_time_det_computation);
        printf("  Total CPU:            %.3f s\n", cpu_time_matrix_eval + cpu_time_det_computation);
        printf("\nEfficiency:\n");
        printf("  Parallel speedup:     %.2fx\n", 
               (cpu_time_matrix_eval + cpu_time_det_computation) / (parallel_end - parallel_start));
        printf("  Thread utilization:   %.1f%%\n", 
               100.0 * (cpu_time_matrix_eval + cpu_time_det_computation) / 
               ((parallel_end - parallel_start) * omp_get_max_threads()));
        printf("===================================\n");
        */
    } else
    #endif
    {
        // Sequential version (original) - also fix timing here
        #ifdef _OPENMP
        double seq_start = omp_get_wtime();
        #endif
        
        slong *indices = (slong*) calloc(total_vars, sizeof(slong));
        slong eval_count = 0;
        slong det_count = 0;
        
        for (slong point = 0; point < total_points; point++) {
            #ifdef _OPENMP
            double point_start = omp_get_wtime();
            #else
            clock_t point_start = clock();
            #endif
            
            // Calculate current point coordinates
            fq_nmod_t *coords = (fq_nmod_t*) malloc(total_vars * sizeof(fq_nmod_t));
            for (slong i = 0; i < total_vars; i++) {
                fq_nmod_init(coords[i], ctx);
                fq_nmod_set(coords[i], grids[i][indices[i]], ctx);
            }
            
            // Create evaluation matrix
            fq_nmod_mat_t eval_mat;
            fq_nmod_mat_init(eval_mat, size, size, ctx);
            
            // Batch evaluation
            #ifdef _OPENMP
            double matrix_eval_start = omp_get_wtime();
            #else
            clock_t matrix_eval_start = clock();
            #endif
            
            fq_nmod_t *var_vals = NULL;
            fq_nmod_t *param_vals = NULL;
            
            if (actual_nvars > 0) {
                var_vals = coords;
            }
            if (actual_npars > 0) {
                param_vals = coords + actual_nvars;
            }
            
            fq_evaluate_matrix_at_point_batch(eval_mat, matrix, size, 
                                             var_vals, param_vals, ctx);
            
            #ifdef _OPENMP
            double matrix_eval_time = omp_get_wtime() - matrix_eval_start;
            #else
            double matrix_eval_time = (double)(clock() - matrix_eval_start) / CLOCKS_PER_SEC;
            #endif
            time_matrix_eval += matrix_eval_time;
            eval_count += size * size;
            
            // Compute determinant
            #ifdef _OPENMP
            double det_start = omp_get_wtime();
            #else
            clock_t det_start = clock();
            #endif
            
            fq_nmod_mat_det(values[point], eval_mat, ctx);
            
            #ifdef _OPENMP
            double det_time = omp_get_wtime() - det_start;
            #else
            double det_time = (double)(clock() - det_start) / CLOCKS_PER_SEC;
            #endif
            time_det_computation += det_time;
            det_count++;
            
            // Cleanup
            fq_nmod_mat_clear(eval_mat, ctx);
            for (slong i = 0; i < total_vars; i++) {
                fq_nmod_clear(coords[i], ctx);
            }
            free(coords);
            
            // Update indices
            for (slong i = 0; i < total_vars; i++) {
                indices[i]++;
                if (indices[i] < grid_sizes[i]) break;
                indices[i] = 0;
            }
            
            // Progress report
            if (point % 1 == 0 || point == total_points - 1) {
                printf("\rProgress: %ld/%ld", 
                       point + 1, total_points);
                fflush(stdout);
            }
        }
        
        free(indices);
    }
    /*
    printf("DEBUG: Generated grids:\n");
    for (slong var = 0; var < total_vars && var < 3; var++) {
        printf("Variable %ld: ", var);
        for (slong j = 0; j < grid_sizes[var] && j < 5; j++) {
            fq_nmod_print_pretty(grids[var][j], ctx);
            printf(" ");
        }
        printf("\n");
    }
    */
    printf("\n");
    
    FQ_INTERP_PRINT("Starting tensor interpolation...\n");
    
    // Tensor product interpolation (unchanged)
    #ifdef _OPENMP
    double interp_start = omp_get_wtime();
    #else
    clock_t interp_start = clock();
    #endif
    
    fq_tensor_interpolation_all_vars_optimized(result, (const fq_nmod_t**)grids, values, 
                                              grid_sizes, actual_nvars, actual_npars, ctx);
    
    #ifdef _OPENMP
    time_interpolation = omp_get_wtime() - interp_start;
    #else
    time_interpolation = (double)(clock() - interp_start) / CLOCKS_PER_SEC;
    #endif
    
    FQ_INTERP_PRINT("Interpolation complete: %ld terms\n", result->nterms);
    
    // Calculate total time - FIX: Use consistent timing
    #ifdef _OPENMP
    double total_time = omp_get_wtime() - start_total_wtime;
    #else
    double total_time = (double)(clock() - start_total) / CLOCKS_PER_SEC;
    #endif
    
    double time_other = total_time - time_preprocessing - time_matrix_eval - time_det_computation - time_interpolation;
    /*
    // Time statistics report
    printf("\n=== Optimized FQ Interpolation Time Statistics ===\n");
    printf("Total computation time: %.3f seconds\n", total_time);
    #ifdef _OPENMP
    if (USE_PARALLEL && total_points > 10) {
        printf("Parallel execution with %d threads\n", omp_get_max_threads());
    }
    #endif
    
    printf("\nBreakdown by operation:\n");
    printf("  1. Preprocessing (bounds, grids):    %.3f s (%.1f%%)\n", 
           time_preprocessing, 100.0 * time_preprocessing / total_time);
    printf("  2. Matrix evaluation (OPTIMIZED):    %.3f s (%.1f%%)\n", 
           time_matrix_eval, 100.0 * time_matrix_eval / total_time);
    printf("  3. Determinant computation:          %.3f s (%.1f%%)\n", 
           time_det_computation, 100.0 * time_det_computation / total_time);
    printf("  4. Tensor interpolation:             %.3f s (%.1f%%)\n", 
           time_interpolation, 100.0 * time_interpolation / total_time);
    printf("  5. Other overhead:                   %.3f s (%.1f%%)\n", 
           time_other, 100.0 * time_other / total_time);
    
    printf("\nOptimization effects:\n");
    printf("  - Total interpolation points:        %ld\n", total_points);
    printf("  - Matrix size:                       %ld x %ld\n", size, size);
    
    double efficiency_ratio = time_matrix_eval / time_det_computation;
    printf("  - Matrix eval / Det ratio:           %.2f\n", efficiency_ratio);
    
    if (efficiency_ratio < 1.0) {
        printf("  → SUCCESS: Matrix evaluation optimized effectively!\n");
    } else if (efficiency_ratio < 2.0) {
        printf("  → GOOD: Matrix evaluation well optimized\n");
    } else {
        printf("  → Room for more optimization in matrix evaluation\n");
    }
    
    printf("=== End Optimized Time Statistics ===\n\n");
    */
    // Cleanup
    for (slong i = 0; i < total_vars; i++) {
        for (slong j = 0; j < grid_sizes[i]; j++) {
            fq_nmod_clear(grids[i][j], ctx);
        }
        free(grids[i]);
    }
    free(grids);
    free(grid_sizes);
    
    for (slong i = 0; i < total_points; i++) {
        fq_nmod_clear(values[i], ctx);
    }
    free(values);
    
    if (computed_bounds) {
        free(computed_bounds);
    }
}
// Compatible interface wrapper
void fq_compute_det_by_interpolation(fq_mvpoly_t *result,
                                     fq_mvpoly_t **matrix,
                                     slong size,
                                     slong nvars,
                                     slong npars,
                                     const fq_nmod_ctx_t ctx,
                                     slong uniform_bound) {
    FQ_INTERP_PRINT("=== Using Optimized FQ Interpolation Algorithm ===\n");
    FQ_INTERP_PRINT("Input parameters: size=%ld, nvars=%ld, npars=%ld\n", size, nvars, npars);
    
    slong actual_nvars = 0;
    slong actual_npars = npars;
    
    if (size > 0) {
        actual_nvars = matrix[0][0].nvars;
        actual_npars = matrix[0][0].npars;
    }
    
    slong total_vars = actual_nvars + actual_npars;
    
    slong *degree_bounds = (slong*) malloc(total_vars * sizeof(slong));
    
    for (slong i = 0; i < total_vars; i++) {
        degree_bounds[i] = uniform_bound;
    }
    
    fq_compute_det_by_interpolation_optimized(result, matrix, size, nvars, npars, 
                                             ctx, degree_bounds);
    
    free(degree_bounds);
}
