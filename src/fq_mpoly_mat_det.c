/*
 * Implementation of optimized polynomial matrix determinant computation
 * Contains all algorithm implementations for various determinant methods
 */

#include "fq_mpoly_mat_det.h"
extern int g_field_equation_reduction;
static inline ulong reduce_exp_field_ui(ulong e, ulong q);
static ulong field_size_q_from_fq_ctx(const fq_nmod_ctx_t ctx);
static void fq_nmod_poly_reduce_field_equation_inplace(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx);
static void fq_nmod_mpoly_reduce_field_equation_inplace(fq_nmod_mpoly_t poly, const fq_nmod_mpoly_ctx_t ctx);
// ============= Timing Utilities Implementation =============
double get_cpu_time(void) {
    return ((double)clock()) / CLOCKS_PER_SEC;
}

timing_info_t start_timing(void) {
    timing_info_t t;
    t.wall_time = get_wall_time();
    t.cpu_time = get_cpu_time();
    return t;
}

timing_info_t end_timing(timing_info_t start) {
    timing_info_t elapsed;
    elapsed.wall_time = get_wall_time() - start.wall_time;
    elapsed.cpu_time = get_cpu_time() - start.cpu_time;
    return elapsed;
}

void print_timing(const char* label, timing_info_t elapsed) {
    DET_PRINT("%s: Wall time: %.6f s, CPU time: %.6f s", 
              label, elapsed.wall_time, elapsed.cpu_time);
    if (elapsed.wall_time > 0) {
        DET_PRINT(" (CPU efficiency: %.1f%%)\n", 
                  (elapsed.cpu_time / elapsed.wall_time) * 100.0);
    } else {
        DET_PRINT("\n");
    }
}

// ============= Polynomial Operation Optimizations Implementation =============

// Custom multiplication for dense polynomials
static inline void poly_mul_dense_optimized(fq_nmod_mpoly_t c, 
                                           const fq_nmod_mpoly_t a,
                                           const fq_nmod_mpoly_t b,
                                           const fq_nmod_mpoly_ctx_t ctx) {
    slong alen = fq_nmod_mpoly_length(a, ctx);
    slong blen = fq_nmod_mpoly_length(b, ctx);
    
    // For very dense polynomials, try different multiplication algorithms
    if (alen > 100 && blen > 100) {
        // Try Johnson's multiplication for dense polynomials
        fq_nmod_mpoly_mul_johnson(c, a, b, ctx);
    } else {
        // Default multiplication
        fq_nmod_mpoly_mul(c, a, b, ctx);
    }
    fq_nmod_mpoly_reduce_field_equation_inplace(c, ctx);
}

// ============= Conversion Functions for Polynomial Recursive Implementation =============

// Convert fq_mvpoly to fq_nmod_poly for a specific variable
void mvpoly_to_fq_nmod_poly(fq_nmod_poly_t poly, const fq_mvpoly_t *mvpoly, 
                           slong var_index, const fq_nmod_ctx_t ctx) {
    fq_nmod_poly_zero(poly, ctx);
    
    slong nvars = mvpoly->nvars;
    
    for (slong t = 0; t < mvpoly->nterms; t++) {
        slong degree = 0;
        int other_vars_zero = 1;
        
        // Check if this term has non-zero exponents in other variables
        if (mvpoly->terms[t].var_exp) {
            for (slong v = 0; v < nvars; v++) {
                if (v == var_index && var_index < nvars) {
                    degree = mvpoly->terms[t].var_exp[v];
                } else if (mvpoly->terms[t].var_exp[v] > 0) {
                    other_vars_zero = 0;
                    break;
                }
            }
        }
        
        if (mvpoly->terms[t].par_exp && other_vars_zero) {
            for (slong p = 0; p < mvpoly->npars; p++) {
                slong idx = nvars + p;
                if (idx == var_index) {
                    degree = mvpoly->terms[t].par_exp[p];
                } else if (mvpoly->terms[t].par_exp[p] > 0) {
                    other_vars_zero = 0;
                    break;
                }
            }
        }
        
        // Only include terms where all other variables have zero exponent
        if (other_vars_zero) {
            fq_nmod_t existing;
            fq_nmod_init(existing, ctx);
            fq_nmod_poly_get_coeff(existing, poly, degree, ctx);
            fq_nmod_add(existing, existing, mvpoly->terms[t].coeff, ctx);
            fq_nmod_poly_set_coeff(poly, degree, existing, ctx);
            fq_nmod_clear(existing, ctx);
        }
    }
}

// Convert fq_nmod_poly back to fq_mvpoly
void fq_nmod_poly_to_mvpoly(fq_mvpoly_t *mvpoly, const fq_nmod_poly_t poly,
                           slong var_index, slong nvars, slong npars,
                           const fq_nmod_ctx_t ctx) {
    fq_mvpoly_init(mvpoly, nvars, npars, ctx);
    
    slong degree = fq_nmod_poly_degree(poly, ctx);
    if (degree < 0) return;
    
    for (slong d = 0; d <= degree; d++) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, ctx);
        fq_nmod_poly_get_coeff(coeff, poly, d, ctx);
        
        if (!fq_nmod_is_zero(coeff, ctx)) {
            slong *var_exp = NULL;
            slong *par_exp = NULL;
            
            if (var_index < nvars && nvars > 0) {
                var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
                var_exp[var_index] = d;
            } else if (var_index >= nvars && npars > 0) {
                par_exp = (slong*) flint_calloc(npars, sizeof(slong));
                par_exp[var_index - nvars] = d;
            }
            
            fq_mvpoly_add_term(mvpoly, var_exp, par_exp, coeff);
            
            if (var_exp) flint_free(var_exp);
            if (par_exp) flint_free(par_exp);
        }
        
        fq_nmod_clear(coeff, ctx);
    }
}

// ============= Polynomial Recursive Determinant Implementation =============

// Recursive determinant computation using fq_nmod_poly operations
void compute_det_poly_recursive_helper(fq_nmod_poly_t det, 
                                      fq_nmod_poly_t **matrix,
                                      slong size, 
                                      const fq_nmod_ctx_t ctx) {
    if (size == 0) {
        fq_nmod_poly_one(det, ctx);
        return;
    }
    
    if (size == 1) {
        fq_nmod_poly_set(det, matrix[0][0], ctx);
        return;
    }
    
    if (size == 2) {
        fq_nmod_poly_t ad, bc;
        fq_nmod_poly_init(ad, ctx);
        fq_nmod_poly_init(bc, ctx);
        
        fq_nmod_poly_mul(ad, matrix[0][0], matrix[1][1], ctx);
        fq_nmod_poly_reduce_field_equation_inplace(ad, ctx);
        fq_nmod_poly_mul(bc, matrix[0][1], matrix[1][0], ctx);
        fq_nmod_poly_reduce_field_equation_inplace(bc, ctx);
        fq_nmod_poly_sub(det, ad, bc, ctx);
        
        fq_nmod_poly_clear(ad, ctx);
        fq_nmod_poly_clear(bc, ctx);
        return;
    }
    
    // General case: Laplace expansion
    fq_nmod_poly_zero(det, ctx);
    
    fq_nmod_poly_t cofactor, subdet;
    fq_nmod_poly_init(cofactor, ctx);
    fq_nmod_poly_init(subdet, ctx);
    
    // Allocate submatrix
    fq_nmod_poly_t **submatrix = (fq_nmod_poly_t**) malloc((size-1) * sizeof(fq_nmod_poly_t*));
    for (slong i = 0; i < size-1; i++) {
        submatrix[i] = (fq_nmod_poly_t*) malloc((size-1) * sizeof(fq_nmod_poly_t));
        for (slong j = 0; j < size-1; j++) {
            fq_nmod_poly_init(submatrix[i][j], ctx);
        }
    }
    
    for (slong col = 0; col < size; col++) {
        if (fq_nmod_poly_is_zero(matrix[0][col], ctx)) {
            continue;
        }
        
        // Build submatrix
        for (slong i = 1; i < size; i++) {
            slong sub_j = 0;
            for (slong j = 0; j < size; j++) {
                if (j != col) {
                    fq_nmod_poly_set(submatrix[i-1][sub_j], matrix[i][j], ctx);
                    sub_j++;
                }
            }
        }
        
        // Recursive call
        compute_det_poly_recursive_helper(subdet, submatrix, size-1, ctx);
        
        // Multiply by matrix element
        fq_nmod_poly_mul(cofactor, matrix[0][col], subdet, ctx);
        fq_nmod_poly_reduce_field_equation_inplace(cofactor, ctx);
        
        // Add or subtract based on sign
        if (col % 2 == 0) {
            fq_nmod_poly_add(det, det, cofactor, ctx);
        } else {
            fq_nmod_poly_sub(det, det, cofactor, ctx);
        }
    }
    
    // Cleanup
    for (slong i = 0; i < size-1; i++) {
        for (slong j = 0; j < size-1; j++) {
            fq_nmod_poly_clear(submatrix[i][j], ctx);
        }
        free(submatrix[i]);
    }
    free(submatrix);
    
    fq_nmod_poly_clear(cofactor, ctx);
    fq_nmod_poly_clear(subdet, ctx);
}

// Main function for polynomial recursive algorithm
void compute_fq_det_poly_recursive(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    if (size <= 0) {
        fq_mvpoly_init(result, matrix[0][0].nvars, matrix[0][0].npars, matrix[0][0].ctx);
        return;
    }
    
    timing_info_t total_start = start_timing();
    
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    slong total_vars = nvars + npars;
    
    DET_PRINT("Computing %ldx%ld determinant via polynomial recursive with Kronecker\n", size, size);
    DET_PRINT("Variables: %ld, Parameters: %ld\n", nvars, npars);
    
    // Special case: if already univariate, no Kronecker needed
    if (total_vars == 1) {
        DET_PRINT("Already univariate, using direct polynomial recursive\n");
        
        // Convert to fq_nmod_poly format
        fq_nmod_poly_t **poly_matrix = (fq_nmod_poly_t**) malloc(size * sizeof(fq_nmod_poly_t*));
        for (slong i = 0; i < size; i++) {
            poly_matrix[i] = (fq_nmod_poly_t*) malloc(size * sizeof(fq_nmod_poly_t));
            for (slong j = 0; j < size; j++) {
                fq_nmod_poly_init(poly_matrix[i][j], ctx);
                
                // Convert mvpoly to poly (direct conversion for univariate)
                for (slong t = 0; t < matrix[i][j].nterms; t++) {
                    slong degree = 0;
                    if (matrix[i][j].terms[t].var_exp && nvars > 0) {
                        degree = matrix[i][j].terms[t].var_exp[0];
                    } else if (matrix[i][j].terms[t].par_exp && npars > 0) {
                        degree = matrix[i][j].terms[t].par_exp[0];
                    }
                    fq_nmod_poly_set_coeff(poly_matrix[i][j], degree, 
                                          matrix[i][j].terms[t].coeff, ctx);
                }
            }
        }
        
        // Compute determinant
        fq_nmod_poly_t det_poly;
        fq_nmod_poly_init(det_poly, ctx);
        compute_det_poly_recursive_helper(det_poly, poly_matrix, size, ctx);
        
        // Convert back to mvpoly
        fq_mvpoly_init(result, nvars, npars, ctx);
        slong degree = fq_nmod_poly_degree(det_poly, ctx);
        for (slong d = 0; d <= degree; d++) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_nmod_poly_get_coeff(coeff, det_poly, d, ctx);
            
            if (!fq_nmod_is_zero(coeff, ctx)) {
                if (nvars > 0) {
                    slong *var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
                    var_exp[0] = d;
                    fq_mvpoly_add_term(result, var_exp, NULL, coeff);
                    flint_free(var_exp);
                } else {
                    slong *par_exp = (slong*) flint_calloc(npars, sizeof(slong));
                    par_exp[0] = d;
                    fq_mvpoly_add_term(result, NULL, par_exp, coeff);
                    flint_free(par_exp);
                }
            }
            
            fq_nmod_clear(coeff, ctx);
        }
        
        // Cleanup
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                fq_nmod_poly_clear(poly_matrix[i][j], ctx);
            }
            free(poly_matrix[i]);
        }
        free(poly_matrix);
        fq_nmod_poly_clear(det_poly, ctx);
        
        timing_info_t total_elapsed = end_timing(total_start);
        print_timing("Total poly recursive (univariate)", total_elapsed);
        return;
    }
    
    // Multivariate case: use Kronecker+HNF
    
    // Step 1: Compute variable bounds (same as in Kronecker algorithm)
    timing_info_t bounds_start = start_timing();
    slong *var_bounds = (slong*) malloc(total_vars * sizeof(slong));
    compute_kronecker_bounds(var_bounds, matrix, size, nvars, npars);
    timing_info_t bounds_elapsed = end_timing(bounds_start);
    
    // Step 2: Compute substitution powers
    slong *substitution_powers = (slong*) malloc(total_vars * sizeof(slong));
    substitution_powers[0] = 1;
    for (slong v = 1; v < total_vars; v++) {
        substitution_powers[v] = substitution_powers[v-1] * var_bounds[v-1];
    }
    
    DET_PRINT("Substitution powers: ");
    for (slong v = 0; v < total_vars; v++) {
        DET_PRINT("%ld ", substitution_powers[v]);
    }
    DET_PRINT("\n");
    
    // Step 3: Convert matrix to univariate using Kronecker
    timing_info_t convert_start = start_timing();
    fq_nmod_poly_t **poly_matrix = (fq_nmod_poly_t**) malloc(size * sizeof(fq_nmod_poly_t*));
    for (slong i = 0; i < size; i++) {
        poly_matrix[i] = (fq_nmod_poly_t*) malloc(size * sizeof(fq_nmod_poly_t));
        for (slong j = 0; j < size; j++) {
            fq_nmod_poly_init(poly_matrix[i][j], ctx);
            mvpoly_to_univariate_kronecker(poly_matrix[i][j], &matrix[i][j], 
                                          substitution_powers, ctx);
        }
    }
    timing_info_t convert_elapsed = end_timing(convert_start);
    
    // Step 4: Compute determinant using recursive algorithm
    timing_info_t det_start = start_timing();
    fq_nmod_poly_t det_poly;
    fq_nmod_poly_init(det_poly, ctx);
    
    compute_det_poly_recursive_helper(det_poly, poly_matrix, size, ctx);
    
    timing_info_t det_elapsed = end_timing(det_start);
    DET_PRINT("Univariate determinant degree: %ld\n", fq_nmod_poly_degree(det_poly, ctx));
    
    // Step 5: Convert back to multivariate
    timing_info_t back_start = start_timing();
    univariate_to_mvpoly_kronecker(result, det_poly, substitution_powers, 
                                  var_bounds, nvars, npars, ctx);
    timing_info_t back_elapsed = end_timing(back_start);
    
    // Cleanup
    free(var_bounds);
    free(substitution_powers);
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_nmod_poly_clear(poly_matrix[i][j], ctx);
        }
        free(poly_matrix[i]);
    }
    free(poly_matrix);
    fq_nmod_poly_clear(det_poly, ctx);
    
    timing_info_t total_elapsed = end_timing(total_start);
    
    printf("\n=== Polynomial Recursive Time Statistics ===\n");
    print_timing("Compute bounds", bounds_elapsed);
    print_timing("Convert to univariate", convert_elapsed);
    print_timing("Recursive determinant", det_elapsed);
    print_timing("Convert back", back_elapsed);
    print_timing("Total poly recursive", total_elapsed);
    printf("Final result: %ld terms\n", result->nterms);
    printf("============================================\n");
}

// ============= Kronecker+HNF Implementation =============

// Compute the Kronecker bound for a multivariate polynomial matrix
void compute_kronecker_bounds(slong *var_bounds, fq_mvpoly_t **matrix, 
                             slong size, slong nvars, slong npars) {
    slong total_vars = nvars + npars;
    
    // Initialize bounds
    for (slong v = 0; v < total_vars; v++) {
        var_bounds[v] = 0;
    }
    
    // Find maximum degree in each variable across all matrix entries
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_mvpoly_t *poly = &matrix[i][j];
            
            for (slong t = 0; t < poly->nterms; t++) {
                // Check variable degrees
                if (poly->terms[t].var_exp && nvars > 0) {
                    for (slong v = 0; v < nvars && v < poly->nvars; v++) {
                        if (poly->terms[t].var_exp[v] > var_bounds[v]) {
                            var_bounds[v] = poly->terms[t].var_exp[v];
                        }
                    }
                }
                
                // Check parameter degrees
                if (poly->terms[t].par_exp && npars > 0) {
                    for (slong p = 0; p < npars && p < poly->npars; p++) {
                        if (poly->terms[t].par_exp[p] > var_bounds[nvars + p]) {
                            var_bounds[nvars + p] = poly->terms[t].par_exp[p];
                        }
                    }
                }
            }
        }
    }
    
    // Compute degree bound for determinant (sum of row maximums)
    for (slong v = 0; v < total_vars; v++) {
        slong det_bound = 0;
        
        for (slong row = 0; row < size; row++) {
            slong row_max = 0;
            
            for (slong col = 0; col < size; col++) {
                fq_mvpoly_t *poly = &matrix[row][col];
                
                for (slong t = 0; t < poly->nterms; t++) {
                    slong deg = 0;
                    
                    if (v < nvars && poly->terms[t].var_exp && v < poly->nvars) {
                        deg = poly->terms[t].var_exp[v];
                    } else if (v >= nvars && poly->terms[t].par_exp && 
                              v - nvars < poly->npars) {
                        deg = poly->terms[t].par_exp[v - nvars];
                    }
                    
                    if (deg > row_max) row_max = deg;
                }
            }
            det_bound += row_max;
        }
        
        var_bounds[v] = det_bound + 1;  // Add 1 for safety
    }
    for (slong v = 0; v < nvars/2; v++) {
        var_bounds[v] = var_bounds[v + nvars/2];
    }
}

// Convert multivariate polynomial to univariate using Kronecker+HNF
void mvpoly_to_univariate_kronecker(fq_nmod_poly_t uni_poly,
                                   const fq_mvpoly_t *mv_poly,
                                   const slong *substitution_powers,
                                   const fq_nmod_ctx_t ctx) {
    fq_nmod_poly_zero(uni_poly, ctx);
    
    if (mv_poly->nterms == 0) return;
    
    slong total_vars = mv_poly->nvars + mv_poly->npars;
    
    for (slong t = 0; t < mv_poly->nterms; t++) {
        slong uni_exp = 0;
        
        // Compute univariate exponent: sum of var_exp[i] * substitution_powers[i]
        if (mv_poly->terms[t].var_exp) {
            for (slong v = 0; v < mv_poly->nvars; v++) {
                uni_exp += mv_poly->terms[t].var_exp[v] * substitution_powers[v];
            }
        }
        
        if (mv_poly->terms[t].par_exp) {
            for (slong p = 0; p < mv_poly->npars; p++) {
                uni_exp += mv_poly->terms[t].par_exp[p] * 
                          substitution_powers[mv_poly->nvars + p];
            }
        }
        
        // Add coefficient at computed degree
        fq_nmod_t existing;
        fq_nmod_init(existing, ctx);
        fq_nmod_poly_get_coeff(existing, uni_poly, uni_exp, ctx);
        fq_nmod_add(existing, existing, mv_poly->terms[t].coeff, ctx);
        fq_nmod_poly_set_coeff(uni_poly, uni_exp, existing, ctx);
        fq_nmod_clear(existing, ctx);
    }
}

// Convert univariate polynomial back to multivariate
void univariate_to_mvpoly_kronecker(fq_mvpoly_t *mv_poly,
                                   const fq_nmod_poly_t uni_poly,
                                   const slong *substitution_powers,
                                   const slong *var_bounds,
                                   slong nvars, slong npars,
                                   const fq_nmod_ctx_t ctx) {
    fq_mvpoly_init(mv_poly, nvars, npars, ctx);
    
    slong degree = fq_nmod_poly_degree(uni_poly, ctx);
    if (degree < 0) return;
    
    slong total_vars = nvars + npars;
    
    for (slong d = 0; d <= degree; d++) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, ctx);
        fq_nmod_poly_get_coeff(coeff, uni_poly, d, ctx);
        
        if (!fq_nmod_is_zero(coeff, ctx)) {
            // Decompose d into multivariate exponents
            slong *var_exp = NULL;
            slong *par_exp = NULL;
            
            if (nvars > 0) {
                var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
            }
            if (npars > 0) {
                par_exp = (slong*) flint_calloc(npars, sizeof(slong));
            }
            
            slong remaining = d;
            
            // Extract exponents in reverse order (largest substitution power first)
            for (slong v = total_vars - 1; v >= 0; v--) {
                slong exp = remaining / substitution_powers[v];
                remaining = remaining % substitution_powers[v];
                
                if (v < nvars && var_exp) {
                    var_exp[v] = exp;
                } else if (v >= nvars && par_exp) {
                    par_exp[v - nvars] = exp;
                }
            }
            
            fq_mvpoly_add_term(mv_poly, var_exp, par_exp, coeff);
            
            if (var_exp) flint_free(var_exp);
            if (par_exp) flint_free(par_exp);
        }
        
        fq_nmod_clear(coeff, ctx);
    }
}

// Compute determinant using Kronecker+HNF
void compute_fq_det_kronecker(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    if (size <= 0) {
        fq_mvpoly_init(result, matrix[0][0].nvars, matrix[0][0].npars, matrix[0][0].ctx);
        return;
    }
    
    timing_info_t total_start = start_timing();
    
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    slong total_vars = nvars + npars;
    
    DET_PRINT("Computing %ldx%ld determinant via Kronecker+HNF\n", size, size);
    DET_PRINT("Variables: %ld, Parameters: %ld\n", nvars, npars);
    
    // Step 1: Compute variable bounds
    timing_info_t bounds_start = start_timing();
    slong *var_bounds = (slong*) malloc(total_vars * sizeof(slong));
    compute_kronecker_bounds(var_bounds, matrix, size, nvars, npars);
    timing_info_t bounds_elapsed = end_timing(bounds_start);
    
    DET_PRINT("Variable bounds: ");
    for (slong v = 0; v < total_vars; v++) {
        DET_PRINT("%ld ", var_bounds[v]);
    }
    DET_PRINT("\n");
    
    // Step 2: Compute substitution powers
    slong *substitution_powers = (slong*) malloc(total_vars * sizeof(slong));
    substitution_powers[0] = 1;
    for (slong v = 1; v < total_vars; v++) {
        substitution_powers[v] = substitution_powers[v-1] * var_bounds[v-1];
    }
    
    DET_PRINT("Substitution powers: ");
    for (slong v = 0; v < total_vars; v++) {
        DET_PRINT("%ld ", substitution_powers[v]);
    }
    DET_PRINT("\n");
    
    // Step 3: Convert matrix to univariate
    timing_info_t convert_start = start_timing();
    fq_nmod_poly_mat_t uni_mat;
    fq_nmod_poly_mat_init(uni_mat, size, size, ctx);
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            mvpoly_to_univariate_kronecker(fq_nmod_poly_mat_entry(uni_mat, i, j),
                                          &matrix[i][j], substitution_powers, ctx);
        }
    }
    timing_info_t convert_elapsed = end_timing(convert_start);
    
    // Check maximum degree
    slong max_degree = 0;
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            slong deg = fq_nmod_poly_degree(fq_nmod_poly_mat_entry(uni_mat, i, j), ctx);
            if (deg > max_degree) max_degree = deg;
        }
    }
    DET_PRINT("Maximum univariate degree after conversion: %ld\n", max_degree);
    
    // Step 4: Compute univariate determinant
    timing_info_t det_start = start_timing();
    fq_nmod_poly_t det_poly;
    fq_nmod_poly_init(det_poly, ctx);
    
    // Use the optimized univariate determinant function
    fq_nmod_poly_mat_det_iter(det_poly, uni_mat, ctx);
    
    timing_info_t det_elapsed = end_timing(det_start);
    DET_PRINT("Univariate determinant degree: %ld\n", fq_nmod_poly_degree(det_poly, ctx));
    
    // Step 5: Convert back to multivariate
    timing_info_t back_start = start_timing();
    univariate_to_mvpoly_kronecker(result, det_poly, substitution_powers, 
                                  var_bounds, nvars, npars, ctx);
    timing_info_t back_elapsed = end_timing(back_start);
    
    // Cleanup
    free(var_bounds);
    free(substitution_powers);
    fq_nmod_poly_clear(det_poly, ctx);
    fq_nmod_poly_mat_clear(uni_mat, ctx);
    
    timing_info_t total_elapsed = end_timing(total_start);
    
    printf("\n=== Kronecker+HNF Time Statistics ===\n");
    print_timing("Compute bounds", bounds_elapsed);
    print_timing("Convert to univariate", convert_elapsed);
    print_timing("Univariate determinant", det_elapsed);
    print_timing("Convert back", back_elapsed);
    print_timing("Total Kronecker", total_elapsed);
    printf("Final result: %ld terms\n", result->nterms);
    printf("==============================================\n");
}

// ============= Prime Field Conversion Functions Implementation =============

void fq_mvpoly_to_nmod_mpoly(nmod_mpoly_t mpoly, const fq_mvpoly_t *poly, 
                             nmod_mpoly_ctx_t mpoly_ctx) {
    nmod_mpoly_zero(mpoly, mpoly_ctx);
    
    if (poly->nterms == 0) return;
    
    slong total_vars = poly->nvars + poly->npars;
    
    // Pre-allocate space for better performance
    nmod_mpoly_fit_length(mpoly, poly->nterms, mpoly_ctx);
    
    for (slong i = 0; i < poly->nterms; i++) {
        ulong *exps = (ulong*) flint_calloc(total_vars, sizeof(ulong));
        
        if (poly->terms[i].var_exp && poly->nvars > 0) {
            for (slong j = 0; j < poly->nvars; j++) {
                exps[j] = (ulong)poly->terms[i].var_exp[j];
            }
        }
        
        if (poly->terms[i].par_exp && poly->npars > 0) {
            for (slong j = 0; j < poly->npars; j++) {
                exps[poly->nvars + j] = (ulong)poly->terms[i].par_exp[j];
            }
        }
        
        // For prime fields, extract the coefficient as ulong
        ulong coeff_val = nmod_poly_get_coeff_ui(poly->terms[i].coeff, 0);
        nmod_mpoly_push_term_ui_ui(mpoly, coeff_val, exps, mpoly_ctx);
        flint_free(exps);
    }
    
    nmod_mpoly_sort_terms(mpoly, mpoly_ctx);
    nmod_mpoly_combine_like_terms(mpoly, mpoly_ctx);
}

void nmod_mpoly_to_fq_mvpoly(fq_mvpoly_t *result, const nmod_mpoly_t poly,
                             slong nvars, slong npars,
                             const nmod_mpoly_ctx_t mpoly_ctx,
                             const fq_nmod_ctx_t field_ctx) {
    fq_mvpoly_init(result, nvars, npars, field_ctx);
    
    slong nterms = nmod_mpoly_length(poly, mpoly_ctx);
    if (nterms == 0) return;
    
    slong total_vars = nmod_mpoly_ctx_nvars(mpoly_ctx);
    
    if (result->alloc < nterms) {
        result->alloc = nterms;
        result->terms = (fq_monomial_t*) flint_realloc(result->terms, 
                                                        result->alloc * sizeof(fq_monomial_t));
    }
    
    // Batch allocate exponent buffer
    ulong *exp_buffer = (ulong*) flint_malloc(total_vars * sizeof(ulong));
    
    // Batch convert with proper initialization
    for (slong i = 0; i < nterms; i++) {
        // Get coefficient
        mp_limb_t coeff_limb = nmod_mpoly_get_term_coeff_ui(poly, i, mpoly_ctx);
        
        // Initialize coefficient (important!)
        fq_nmod_init(result->terms[i].coeff, field_ctx);
        fq_nmod_set_ui(result->terms[i].coeff, coeff_limb, field_ctx);
        
        // Get exponents
        nmod_mpoly_get_term_exp_ui(exp_buffer, poly, i, mpoly_ctx);
        
        // Allocate and set variable exponents
        if (nvars > 0) {
            result->terms[i].var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
            for (slong j = 0; j < nvars && j < total_vars; j++) {
                result->terms[i].var_exp[j] = (slong)exp_buffer[j];
            }
        } else {
            result->terms[i].var_exp = NULL;
        }
        
        // Allocate and set parameter exponents
        if (npars > 0 && total_vars > nvars) {
            result->terms[i].par_exp = (slong*) flint_calloc(npars, sizeof(slong));
            for (slong j = 0; j < npars && (nvars + j) < total_vars; j++) {
                result->terms[i].par_exp[j] = (slong)exp_buffer[nvars + j];
            }
        } else {
            result->terms[i].par_exp = NULL;
        }
    }
    
    // Set term count
    result->nterms = nterms;
    
    // Cleanup
    flint_free(exp_buffer);
}

void fq_matrix_mvpoly_to_nmod_mpoly(nmod_mpoly_t **mpoly_matrix, 
                                   fq_mvpoly_t **mvpoly_matrix, 
                                   slong size, 
                                   nmod_mpoly_ctx_t mpoly_ctx) {
    DET_PRINT("Converting %ld x %ld matrix to nmod_mpoly\n", size, size);
    
    timing_info_t start = start_timing();
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            nmod_mpoly_init(mpoly_matrix[i][j], mpoly_ctx);
            fq_mvpoly_to_nmod_mpoly(mpoly_matrix[i][j], &mvpoly_matrix[i][j], mpoly_ctx);
        }
    }
    
    timing_info_t elapsed = end_timing(start);
    print_timing("Matrix conversion to nmod_mpoly", elapsed);
}

// ============= Prime Field Determinant Computation Implementation =============

// Hand-optimized 3x3 determinant for nmod_mpoly
void compute_det_3x3_nmod_optimized(nmod_mpoly_t det, 
                                   nmod_mpoly_t **m,
                                   nmod_mpoly_ctx_t ctx) {
    nmod_mpoly_t t1, t2, t3, t4, t5, t6, sum;
    
    // Initialize temporaries
    nmod_mpoly_init(t1, ctx);
    nmod_mpoly_init(t2, ctx);
    nmod_mpoly_init(t3, ctx);
    nmod_mpoly_init(t4, ctx);
    nmod_mpoly_init(t5, ctx);
    nmod_mpoly_init(t6, ctx);
    nmod_mpoly_init(sum, ctx);
    
    // Compute 6 products in parallel if beneficial
    #pragma omp parallel sections if(omp_get_max_threads() > 2)
    {
        #pragma omp section
        {
            nmod_mpoly_mul(t1, m[1][1], m[2][2], ctx);
            nmod_mpoly_mul(t1, m[0][0], t1, ctx);
        }
        #pragma omp section
        {
            nmod_mpoly_mul(t2, m[1][2], m[2][0], ctx);
            nmod_mpoly_mul(t2, m[0][1], t2, ctx);
        }
        #pragma omp section
        {
            nmod_mpoly_mul(t3, m[1][0], m[2][1], ctx);
            nmod_mpoly_mul(t3, m[0][2], t3, ctx);
        }
        #pragma omp section
        {
            nmod_mpoly_mul(t4, m[1][0], m[2][2], ctx);
            nmod_mpoly_mul(t4, m[0][1], t4, ctx);
        }
        #pragma omp section
        {
            nmod_mpoly_mul(t5, m[1][1], m[2][0], ctx);
            nmod_mpoly_mul(t5, m[0][2], t5, ctx);
        }
        #pragma omp section
        {
            nmod_mpoly_mul(t6, m[1][2], m[2][1], ctx);
            nmod_mpoly_mul(t6, m[0][0], t6, ctx);
        }
    }
    
    // Sum with signs
    nmod_mpoly_add(sum, t1, t2, ctx);
    nmod_mpoly_add(sum, sum, t3, ctx);
    nmod_mpoly_sub(sum, sum, t4, ctx);
    nmod_mpoly_sub(sum, sum, t5, ctx);
    nmod_mpoly_sub(det, sum, t6, ctx);
    
    // Cleanup
    nmod_mpoly_clear(t1, ctx);
    nmod_mpoly_clear(t2, ctx);
    nmod_mpoly_clear(t3, ctx);
    nmod_mpoly_clear(t4, ctx);
    nmod_mpoly_clear(t5, ctx);
    nmod_mpoly_clear(t6, ctx);
    nmod_mpoly_clear(sum, ctx);
}

// Recursive determinant for nmod_mpoly
void compute_nmod_mpoly_det_recursive(nmod_mpoly_t det_result, 
                                     nmod_mpoly_t **mpoly_matrix, 
                                     slong size, 
                                     nmod_mpoly_ctx_t mpoly_ctx) {
    if (size <= 0) {
        nmod_mpoly_one(det_result, mpoly_ctx);
        return;
    }
    
    if (size == 1) {
        nmod_mpoly_set(det_result, mpoly_matrix[0][0], mpoly_ctx);
        return;
    }
    
    if (size == 2) {
        nmod_mpoly_t ad, bc;
        nmod_mpoly_init(ad, mpoly_ctx);
        nmod_mpoly_init(bc, mpoly_ctx);
        
        nmod_mpoly_mul(ad, mpoly_matrix[0][0], mpoly_matrix[1][1], mpoly_ctx);
        nmod_mpoly_mul(bc, mpoly_matrix[0][1], mpoly_matrix[1][0], mpoly_ctx);
        nmod_mpoly_sub(det_result, ad, bc, mpoly_ctx);
        
        nmod_mpoly_clear(ad, mpoly_ctx);
        nmod_mpoly_clear(bc, mpoly_ctx);
        return;
    }
    
    if (size == 3) {
        compute_det_3x3_nmod_optimized(det_result, mpoly_matrix, mpoly_ctx);
        return;
    }
    
    // General case: Laplace expansion
    nmod_mpoly_zero(det_result, mpoly_ctx);
    
    nmod_mpoly_t temp_result, cofactor, subdet;
    nmod_mpoly_init(temp_result, mpoly_ctx);
    nmod_mpoly_init(cofactor, mpoly_ctx);
    nmod_mpoly_init(subdet, mpoly_ctx);
    
    for (slong col = 0; col < size; col++) {
        if (nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
            continue;
        }
        
        // Create submatrix
        nmod_mpoly_t **submatrix = (nmod_mpoly_t**) flint_malloc((size-1) * sizeof(nmod_mpoly_t*));
        for (slong i = 0; i < size-1; i++) {
            submatrix[i] = (nmod_mpoly_t*) flint_malloc((size-1) * sizeof(nmod_mpoly_t));
            for (slong j = 0; j < size-1; j++) {
                nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
            }
        }
        
        // Fill submatrix
        for (slong i = 1; i < size; i++) {
            slong sub_j = 0;
            for (slong j = 0; j < size; j++) {
                if (j != col) {
                    nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                    sub_j++;
                }
            }
        }
        
        // Recursive computation
        compute_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
        
        // Compute cofactor
        nmod_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
        
        // Add/subtract to result
        if (col % 2 == 0) {
            nmod_mpoly_add(temp_result, det_result, cofactor, mpoly_ctx);
        } else {
            nmod_mpoly_sub(temp_result, det_result, cofactor, mpoly_ctx);
        }
        nmod_mpoly_set(det_result, temp_result, mpoly_ctx);
        
        // Cleanup submatrix
        for (slong i = 0; i < size-1; i++) {
            for (slong j = 0; j < size-1; j++) {
                nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
            }
            flint_free(submatrix[i]);
        }
        flint_free(submatrix);
    }
    
    nmod_mpoly_clear(temp_result, mpoly_ctx);
    nmod_mpoly_clear(cofactor, mpoly_ctx);
    nmod_mpoly_clear(subdet, mpoly_ctx);
}

// Parallel determinant computation for nmod_mpoly with proper nested parallelism
void compute_nmod_mpoly_det_parallel_optimized(nmod_mpoly_t det_result, 
                                              nmod_mpoly_t **mpoly_matrix, 
                                              slong size, 
                                              nmod_mpoly_ctx_t mpoly_ctx,
                                              slong depth) {
    // For deep recursion or small matrices, use sequential
    if (size < PARALLEL_THRESHOLD || depth >= MAX_PARALLEL_DEPTH) {
        compute_nmod_mpoly_det_recursive(det_result, mpoly_matrix, size, mpoly_ctx);
        return;
    }
    
    DET_PRINT("Parallel nmod computation for %ld x %ld matrix (depth %ld)\n", size, size, depth);
    
    if (size <= 3) {
        compute_nmod_mpoly_det_recursive(det_result, mpoly_matrix, size, mpoly_ctx);
        return;
    }
    
    nmod_mpoly_zero(det_result, mpoly_ctx);
    
    // Count non-zero entries in first row
    slong nonzero_count = 0;
    for (slong col = 0; col < size; col++) {
        if (!nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
            nonzero_count++;
        }
    }
    
    if (nonzero_count < 2) {
        compute_nmod_mpoly_det_recursive(det_result, mpoly_matrix, size, mpoly_ctx);
        return;
    }
    
    // Allocate space for partial results
    nmod_mpoly_t *partial_results = (nmod_mpoly_t*) flint_malloc(size * sizeof(nmod_mpoly_t));
    for (slong i = 0; i < size; i++) {
        nmod_mpoly_init(partial_results[i], mpoly_ctx);
        nmod_mpoly_zero(partial_results[i], mpoly_ctx);
    }
    
    // Determine parallelism strategy based on depth
    if (depth == 0) {
        // First level: use parallel for with nested parallelism enabled
        #pragma omp parallel for schedule(static) num_threads(FLINT_MIN(nonzero_count, omp_get_max_threads()))
        for (slong col = 0; col < size; col++) {
            if (nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
                continue;
            }
            
            nmod_mpoly_t cofactor, subdet;
            nmod_mpoly_init(cofactor, mpoly_ctx);
            nmod_mpoly_init(subdet, mpoly_ctx);
            
            // Create submatrix
            nmod_mpoly_t **submatrix = (nmod_mpoly_t**) flint_malloc((size-1) * sizeof(nmod_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (nmod_mpoly_t*) flint_malloc((size-1) * sizeof(nmod_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
                }
            }
            
            // Fill submatrix
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                        sub_j++;
                    }
                }
            }
            
            // Recursive call - this will use nested parallelism at depth 1
            compute_nmod_mpoly_det_parallel_optimized(subdet, submatrix, size-1, mpoly_ctx, depth+1);
            
            // Compute cofactor
            nmod_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
            
            // Store with sign
            if (col % 2 == 0) {
                nmod_mpoly_set(partial_results[col], cofactor, mpoly_ctx);
            } else {
                nmod_mpoly_neg(partial_results[col], cofactor, mpoly_ctx);
            }
            
            // Cleanup
            for (slong i = 0; i < size-1; i++) {
                for (slong j = 0; j < size-1; j++) {
                    nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
                }
                flint_free(submatrix[i]);
            }
            flint_free(submatrix);
            
            nmod_mpoly_clear(cofactor, mpoly_ctx);
            nmod_mpoly_clear(subdet, mpoly_ctx);
        }
    } else if (depth == 1 && size >= PARALLEL_THRESHOLD) {
        // Second level: also use parallel for, but with fewer threads
        slong max_threads_level2 = FLINT_MAX(1, omp_get_max_threads() / size);
        
        #pragma omp parallel for schedule(static) num_threads(FLINT_MIN(nonzero_count, max_threads_level2)) if(nonzero_count >= 3)
        for (slong col = 0; col < size; col++) {
            if (nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
                continue;
            }
            
            nmod_mpoly_t cofactor, subdet;
            nmod_mpoly_init(cofactor, mpoly_ctx);
            nmod_mpoly_init(subdet, mpoly_ctx);
            
            // Create submatrix
            nmod_mpoly_t **submatrix = (nmod_mpoly_t**) flint_malloc((size-1) * sizeof(nmod_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (nmod_mpoly_t*) flint_malloc((size-1) * sizeof(nmod_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
                }
            }
            
            // Fill submatrix
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                        sub_j++;
                    }
                }
            }
            
            // At depth > 1, use sequential computation
            compute_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
            
            // Compute cofactor
            nmod_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
            
            // Store with sign
            if (col % 2 == 0) {
                nmod_mpoly_set(partial_results[col], cofactor, mpoly_ctx);
            } else {
                nmod_mpoly_neg(partial_results[col], cofactor, mpoly_ctx);
            }
            
            // Cleanup
            for (slong i = 0; i < size-1; i++) {
                for (slong j = 0; j < size-1; j++) {
                    nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
                }
                flint_free(submatrix[i]);
            }
            flint_free(submatrix);
            
            nmod_mpoly_clear(cofactor, mpoly_ctx);
            nmod_mpoly_clear(subdet, mpoly_ctx);
        }
    } else {
        // Sequential fallback for deeper levels
        for (slong col = 0; col < size; col++) {
            if (nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
                continue;
            }
            
            nmod_mpoly_t cofactor, subdet;
            nmod_mpoly_init(cofactor, mpoly_ctx);
            nmod_mpoly_init(subdet, mpoly_ctx);
            
            // Create and fill submatrix
            nmod_mpoly_t **submatrix = (nmod_mpoly_t**) flint_malloc((size-1) * sizeof(nmod_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (nmod_mpoly_t*) flint_malloc((size-1) * sizeof(nmod_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
                }
            }
            
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                        sub_j++;
                    }
                }
            }
            
            // Sequential computation
            compute_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
            
            // Compute cofactor and store
            nmod_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
            if (col % 2 == 0) {
                nmod_mpoly_set(partial_results[col], cofactor, mpoly_ctx);
            } else {
                nmod_mpoly_neg(partial_results[col], cofactor, mpoly_ctx);
            }
            
            // Cleanup
            for (slong i = 0; i < size-1; i++) {
                for (slong j = 0; j < size-1; j++) {
                    nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
                }
                flint_free(submatrix[i]);
            }
            flint_free(submatrix);
            
            nmod_mpoly_clear(cofactor, mpoly_ctx);
            nmod_mpoly_clear(subdet, mpoly_ctx);
        }
    }
    
    // Sum results (sequential to avoid race conditions)
    nmod_mpoly_t temp_sum;
    nmod_mpoly_init(temp_sum, mpoly_ctx);
    
    for (slong col = 0; col < size; col++) {
        if (!nmod_mpoly_is_zero(partial_results[col], mpoly_ctx)) {
            nmod_mpoly_add(temp_sum, det_result, partial_results[col], mpoly_ctx);
            nmod_mpoly_set(det_result, temp_sum, mpoly_ctx);
        }
        nmod_mpoly_clear(partial_results[col], mpoly_ctx);
    }
    
    nmod_mpoly_clear(temp_sum, mpoly_ctx);
    flint_free(partial_results);
}

// ============= Univariate Optimization Implementation =============

int is_univariate_matrix(fq_mvpoly_t **matrix, slong size) {
    if (size == 0) return 0;
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    return (nvars == 1 && npars == 0);
}

void compute_fq_det_univariate_optimized(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    if (size <= 0) {
        fq_mvpoly_init(result, matrix[0][0].nvars, matrix[0][0].npars, matrix[0][0].ctx);
        return;
    }
    
    DET_PRINT("Using univariate polynomial matrix optimization for %ldx%ld matrix\n", size, size);
    
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    fq_mvpoly_init(result, 1, 0, ctx);
    
    timing_info_t start = start_timing();
    
    fq_nmod_poly_mat_t poly_mat;
    fq_nmod_poly_mat_init(poly_mat, size, size, ctx);
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_nmod_poly_struct *entry = fq_nmod_poly_mat_entry(poly_mat, i, j);
            fq_nmod_poly_zero(entry, ctx);
            
            for (slong k = 0; k < matrix[i][j].nterms; k++) {
                fq_monomial_t *term = &matrix[i][j].terms[k];
                slong degree = 0;
                if (term->var_exp && matrix[i][j].nvars > 0) {
                    degree = term->var_exp[0];
                }
                fq_nmod_poly_set_coeff(entry, degree, term->coeff, ctx);
            }
        }
    }
    
    fq_nmod_poly_t det_poly;
    fq_nmod_poly_init(det_poly, ctx);
    
    fq_nmod_poly_mat_det_iter(det_poly, poly_mat, ctx);
    
    timing_info_t conv_elapsed = end_timing(start);
    print_timing("Univariate matrix determinant", conv_elapsed);
    
    slong degree = fq_nmod_poly_degree(det_poly, ctx);
    if (degree >= 0) {
        for (slong d = 0; d <= degree; d++) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_nmod_poly_get_coeff(coeff, det_poly, d, ctx);
            
            if (!fq_nmod_is_zero(coeff, ctx)) {
                slong *var_exp = (slong*) flint_calloc(1, sizeof(slong));
                var_exp[0] = d;
                fq_mvpoly_add_term(result, var_exp, NULL, coeff);
                flint_free(var_exp);
            }
            
            fq_nmod_clear(coeff, ctx);
        }
    }
    
    fq_nmod_poly_clear(det_poly, ctx);
    fq_nmod_poly_mat_clear(poly_mat, ctx);
}

// ============= Conversion Functions Implementation =============

void fq_mvpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t mpoly, const fq_mvpoly_t *poly, 
                               fq_nmod_mpoly_ctx_t mpoly_ctx) {
    fq_nmod_mpoly_zero(mpoly, mpoly_ctx);
    
    if (poly->nterms == 0) return;
    
    slong total_vars = poly->nvars + poly->npars;
    
    // Pre-allocate space for better performance
    fq_nmod_mpoly_fit_length(mpoly, poly->nterms, mpoly_ctx);
    
    for (slong i = 0; i < poly->nterms; i++) {
        ulong *exps = (ulong*) flint_calloc(total_vars, sizeof(ulong));
        
        if (poly->terms[i].var_exp && poly->nvars > 0) {
            for (slong j = 0; j < poly->nvars; j++) {
                exps[j] = (ulong)poly->terms[i].var_exp[j];
            }
        }
        
        if (poly->terms[i].par_exp && poly->npars > 0) {
            for (slong j = 0; j < poly->npars; j++) {
                exps[poly->nvars + j] = (ulong)poly->terms[i].par_exp[j];
            }
        }
        
        fq_nmod_mpoly_push_term_fq_nmod_ui(mpoly, poly->terms[i].coeff, exps, mpoly_ctx);
        flint_free(exps);
    }
    
    fq_nmod_mpoly_sort_terms(mpoly, mpoly_ctx);
    fq_nmod_mpoly_combine_like_terms(mpoly, mpoly_ctx);
}

void fq_nmod_mpoly_to_fq_mvpoly(fq_mvpoly_t *poly, const fq_nmod_mpoly_t mpoly,
                               slong nvars, slong npars, 
                               fq_nmod_mpoly_ctx_t mpoly_ctx, const fq_nmod_ctx_t ctx) {
    fq_mvpoly_init(poly, nvars, npars, ctx);
    slong nterms = fq_nmod_mpoly_length(mpoly, mpoly_ctx);
    if (nterms == 0) return;
    
    slong total_vars = fq_nmod_mpoly_ctx_nvars(mpoly_ctx);
    
    // Pre-allocate the terms array
    if (poly->alloc < nterms) {
        // First, clear any existing terms
        for (slong i = 0; i < poly->nterms; i++) {
            fq_nmod_clear(poly->terms[i].coeff, ctx);
            if (poly->terms[i].var_exp) flint_free(poly->terms[i].var_exp);
            if (poly->terms[i].par_exp) flint_free(poly->terms[i].par_exp);
        }
        
        poly->alloc = nterms;
        poly->terms = (fq_monomial_t*) flint_realloc(poly->terms, 
                                                      poly->alloc * sizeof(fq_monomial_t));
        
        // Initialize all term structures
        for (slong i = 0; i < poly->alloc; i++) {
            // Zero out the structure first
            memset(&poly->terms[i], 0, sizeof(fq_monomial_t));
        }
    }
   
    // Allocate a single exponent buffer for reading
    ulong *exp_buffer = (ulong*) flint_malloc(total_vars * sizeof(ulong));
    
    // Process all terms - but allocate individually for compatibility
    for (slong i = 0; i < nterms; i++) {
        // Initialize coefficient
        fq_nmod_init(poly->terms[i].coeff, ctx);
        fq_nmod_mpoly_get_term_coeff_fq_nmod(poly->terms[i].coeff, mpoly, i, mpoly_ctx);
        // Get exponents for this term
        fq_nmod_mpoly_get_term_exp_ui(exp_buffer, mpoly, i, mpoly_ctx);
       
        // Allocate and set variable exponents
        if (nvars > 0) {
            poly->terms[i].var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
            for (slong j = 0; j < nvars && j < total_vars; j++) {
                poly->terms[i].var_exp[j] = (slong)exp_buffer[j];
            }
        } else {
            poly->terms[i].var_exp = NULL;
        }
        // Allocate and set parameter exponents
        if (npars > 0 && total_vars > nvars) {
            poly->terms[i].par_exp = (slong*) flint_calloc(npars, sizeof(slong));
            for (slong j = 0; j < npars && (nvars + j) < total_vars; j++) {
                poly->terms[i].par_exp[j] = (slong)exp_buffer[nvars + j];
            }
        } else {
            poly->terms[i].par_exp = NULL;
        }
    }

// Set the number of terms
    poly->nterms = nterms;
    if (g_field_equation_reduction) {
        fq_mvpoly_reduce_field_equation(poly);
    }
    
    // Cleanup
    flint_free(exp_buffer);
}

void fq_matrix_mvpoly_to_mpoly(fq_nmod_mpoly_t **mpoly_matrix, 
                              fq_mvpoly_t **mvpoly_matrix, 
                              slong size, 
                              fq_nmod_mpoly_ctx_t mpoly_ctx) {
    DET_PRINT("Converting %ld x %ld matrix\n", size, size);
    
    timing_info_t start = start_timing();
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_nmod_mpoly_init(mpoly_matrix[i][j], mpoly_ctx);
            fq_mvpoly_to_fq_nmod_mpoly(mpoly_matrix[i][j], &mvpoly_matrix[i][j], mpoly_ctx);
        }
    }
    
    timing_info_t elapsed = end_timing(start);
    print_timing("Matrix conversion", elapsed);
}

// ============= Optimized Determinant Computation Implementation =============

// Hand-optimized 3x3 determinant
void compute_det_3x3_optimized(fq_nmod_mpoly_t det, 
                              fq_nmod_mpoly_t **m,
                              fq_nmod_mpoly_ctx_t ctx) {
    fq_nmod_mpoly_t t1, t2, t3, t4, t5, t6, sum;
    
    // Initialize temporaries
    fq_nmod_mpoly_init(t1, ctx);
    fq_nmod_mpoly_init(t2, ctx);
    fq_nmod_mpoly_init(t3, ctx);
    fq_nmod_mpoly_init(t4, ctx);
    fq_nmod_mpoly_init(t5, ctx);
    fq_nmod_mpoly_init(t6, ctx);
    fq_nmod_mpoly_init(sum, ctx);
    
    // Compute 6 products in parallel if beneficial
    #pragma omp parallel sections if(omp_get_max_threads() > 2)
    {
        #pragma omp section
        {
            fq_nmod_mpoly_mul(t1, m[1][1], m[2][2], ctx);
            poly_mul_dense_optimized(t1, m[0][0], t1, ctx);
        }
        #pragma omp section
        {
            fq_nmod_mpoly_mul(t2, m[1][2], m[2][0], ctx);
            poly_mul_dense_optimized(t2, m[0][1], t2, ctx);
        }
        #pragma omp section
        {
            fq_nmod_mpoly_mul(t3, m[1][0], m[2][1], ctx);
            poly_mul_dense_optimized(t3, m[0][2], t3, ctx);
        }
        #pragma omp section
        {
            fq_nmod_mpoly_mul(t4, m[1][0], m[2][2], ctx);
            poly_mul_dense_optimized(t4, m[0][1], t4, ctx);
        }
        #pragma omp section
        {
            fq_nmod_mpoly_mul(t5, m[1][1], m[2][0], ctx);
            poly_mul_dense_optimized(t5, m[0][2], t5, ctx);
        }
        #pragma omp section
        {
            fq_nmod_mpoly_mul(t6, m[1][2], m[2][1], ctx);
            poly_mul_dense_optimized(t6, m[0][0], t6, ctx);
        }
    }
    
    // Sum with signs
    fq_nmod_mpoly_add(sum, t1, t2, ctx);
    fq_nmod_mpoly_add(sum, sum, t3, ctx);
    fq_nmod_mpoly_sub(sum, sum, t4, ctx);
    fq_nmod_mpoly_sub(sum, sum, t5, ctx);
    fq_nmod_mpoly_sub(det, sum, t6, ctx);
    
    // Cleanup
    fq_nmod_mpoly_clear(t1, ctx);
    fq_nmod_mpoly_clear(t2, ctx);
    fq_nmod_mpoly_clear(t3, ctx);
    fq_nmod_mpoly_clear(t4, ctx);
    fq_nmod_mpoly_clear(t5, ctx);
    fq_nmod_mpoly_clear(t6, ctx);
    fq_nmod_mpoly_clear(sum, ctx);
}

// Recursive determinant with optimizations
void compute_fq_nmod_mpoly_det_recursive(fq_nmod_mpoly_t det_result, 
                                        fq_nmod_mpoly_t **mpoly_matrix, 
                                        slong size, 
                                        fq_nmod_mpoly_ctx_t mpoly_ctx) {
    if (size <= 0) {
        fq_nmod_mpoly_one(det_result, mpoly_ctx);
        return;
    }
    
    if (size == 1) {
        fq_nmod_mpoly_set(det_result, mpoly_matrix[0][0], mpoly_ctx);
        return;
    }
    
    if (size == 2) {
        fq_nmod_mpoly_t ad, bc;
        fq_nmod_mpoly_init(ad, mpoly_ctx);
        fq_nmod_mpoly_init(bc, mpoly_ctx);
        
        poly_mul_dense_optimized(ad, mpoly_matrix[0][0], mpoly_matrix[1][1], mpoly_ctx);
        poly_mul_dense_optimized(bc, mpoly_matrix[0][1], mpoly_matrix[1][0], mpoly_ctx);
        fq_nmod_mpoly_sub(det_result, ad, bc, mpoly_ctx);
        
        fq_nmod_mpoly_clear(ad, mpoly_ctx);
        fq_nmod_mpoly_clear(bc, mpoly_ctx);
        return;
    }
    
    if (size == 3) {
        compute_det_3x3_optimized(det_result, mpoly_matrix, mpoly_ctx);
        return;
    }
    
    // General case: Laplace expansion
    fq_nmod_mpoly_zero(det_result, mpoly_ctx);
    
    fq_nmod_mpoly_t temp_result, cofactor, subdet;
    fq_nmod_mpoly_init(temp_result, mpoly_ctx);
    fq_nmod_mpoly_init(cofactor, mpoly_ctx);
    fq_nmod_mpoly_init(subdet, mpoly_ctx);
    
    for (slong col = 0; col < size; col++) {
        if (fq_nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
            continue;
        }
        
        // Create submatrix
        fq_nmod_mpoly_t **submatrix = (fq_nmod_mpoly_t**) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t*));
        for (slong i = 0; i < size-1; i++) {
            submatrix[i] = (fq_nmod_mpoly_t*) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t));
            for (slong j = 0; j < size-1; j++) {
                fq_nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
            }
        }
        
        // Fill submatrix
        for (slong i = 1; i < size; i++) {
            slong sub_j = 0;
            for (slong j = 0; j < size; j++) {
                if (j != col) {
                    fq_nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                    sub_j++;
                }
            }
        }
        
        // Recursive computation
        compute_fq_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
        
        // Compute cofactor
        poly_mul_dense_optimized(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
        
        // Add/subtract to result
        if (col % 2 == 0) {
            fq_nmod_mpoly_add(temp_result, det_result, cofactor, mpoly_ctx);
        } else {
            fq_nmod_mpoly_sub(temp_result, det_result, cofactor, mpoly_ctx);
        }
        fq_nmod_mpoly_set(det_result, temp_result, mpoly_ctx);
        
        // Cleanup submatrix
        for (slong i = 0; i < size-1; i++) {
            for (slong j = 0; j < size-1; j++) {
                fq_nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
            }
            flint_free(submatrix[i]);
        }
        flint_free(submatrix);
    }
    
    fq_nmod_mpoly_clear(temp_result, mpoly_ctx);
    fq_nmod_mpoly_clear(cofactor, mpoly_ctx);
    fq_nmod_mpoly_clear(subdet, mpoly_ctx);
}

void compute_fq_nmod_mpoly_det_parallel_optimized(fq_nmod_mpoly_t det_result, 
                                                  fq_nmod_mpoly_t **mpoly_matrix, 
                                                  slong size, 
                                                  fq_nmod_mpoly_ctx_t mpoly_ctx,
                                                  slong depth) {
    if (size < PARALLEL_THRESHOLD || depth >= MAX_PARALLEL_DEPTH) {
        compute_fq_nmod_mpoly_det_recursive(det_result, mpoly_matrix, size, mpoly_ctx);
        return;
    }
    
    // DET_PRINT("Parallel computation for %ld x %ld matrix (depth %ld)\n", size, size, depth);
    
    if (size <= 3) {
        compute_fq_nmod_mpoly_det_recursive(det_result, mpoly_matrix, size, mpoly_ctx);
        return;
    }
    
    fq_nmod_mpoly_zero(det_result, mpoly_ctx);
    
    // Count non-zero entries in first row
    slong nonzero_count = 0;
    for (slong col = 0; col < size; col++) {
        if (!fq_nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
            nonzero_count++;
        }
    }
    
    if (nonzero_count < 2) {
        compute_fq_nmod_mpoly_det_recursive(det_result, mpoly_matrix, size, mpoly_ctx);
        return;
    }
    
    // Allocate space for partial results
    fq_nmod_mpoly_t *partial_results = (fq_nmod_mpoly_t*) flint_malloc(size * sizeof(fq_nmod_mpoly_t));
    for (slong i = 0; i < size; i++) {
        fq_nmod_mpoly_init(partial_results[i], mpoly_ctx);
        fq_nmod_mpoly_zero(partial_results[i], mpoly_ctx);
    }
    
    // Determine parallelism strategy based on depth
    if (depth == 0) {
        // First level: use parallel for with nested parallelism enabled
        #pragma omp parallel for schedule(static) num_threads(FLINT_MIN(nonzero_count, omp_get_max_threads()))
        for (slong col = 0; col < size; col++) {
            if (fq_nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
                continue;
            }
            
            fq_nmod_mpoly_t cofactor, subdet;
            fq_nmod_mpoly_init(cofactor, mpoly_ctx);
            fq_nmod_mpoly_init(subdet, mpoly_ctx);
            
            // Create submatrix
            fq_nmod_mpoly_t **submatrix = (fq_nmod_mpoly_t**) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (fq_nmod_mpoly_t*) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    fq_nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
                }
            }
            
            // Fill submatrix
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        fq_nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                        sub_j++;
                    }
                }
            }
            
            // Recursive call - this will use nested parallelism at depth 1
            compute_fq_nmod_mpoly_det_parallel_optimized(subdet, submatrix, size-1, mpoly_ctx, depth+1);
            
            // Compute cofactor
            poly_mul_dense_optimized(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
            
            // Store with sign
            if (col % 2 == 0) {
                fq_nmod_mpoly_set(partial_results[col], cofactor, mpoly_ctx);
            } else {
                fq_nmod_mpoly_neg(partial_results[col], cofactor, mpoly_ctx);
            }
            
            // Cleanup
            for (slong i = 0; i < size-1; i++) {
                for (slong j = 0; j < size-1; j++) {
                    fq_nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
                }
                flint_free(submatrix[i]);
            }
            flint_free(submatrix);
            
            fq_nmod_mpoly_clear(cofactor, mpoly_ctx);
            fq_nmod_mpoly_clear(subdet, mpoly_ctx);
        }
    } else if (depth == 1 && size >= PARALLEL_THRESHOLD) {
        // Second level: also use parallel for, but with fewer threads
        slong max_threads_level2 = FLINT_MAX(1, omp_get_max_threads() / size);
        
        #pragma omp parallel for schedule(static) num_threads(FLINT_MIN(nonzero_count, max_threads_level2)) if(nonzero_count >= 3)
        for (slong col = 0; col < size; col++) {
            if (fq_nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
                continue;
            }
            
            fq_nmod_mpoly_t cofactor, subdet;
            fq_nmod_mpoly_init(cofactor, mpoly_ctx);
            fq_nmod_mpoly_init(subdet, mpoly_ctx);
            
            // Create submatrix
            fq_nmod_mpoly_t **submatrix = (fq_nmod_mpoly_t**) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (fq_nmod_mpoly_t*) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    fq_nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
                }
            }
            
            // Fill submatrix
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        fq_nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                        sub_j++;
                    }
                }
            }
            
            // At depth > 1, use sequential computation
            compute_fq_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
            
            // Compute cofactor
            poly_mul_dense_optimized(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
            
            // Store with sign
            if (col % 2 == 0) {
                fq_nmod_mpoly_set(partial_results[col], cofactor, mpoly_ctx);
            } else {
                fq_nmod_mpoly_neg(partial_results[col], cofactor, mpoly_ctx);
            }
            
            // Cleanup
            for (slong i = 0; i < size-1; i++) {
                for (slong j = 0; j < size-1; j++) {
                    fq_nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
                }
                flint_free(submatrix[i]);
            }
            flint_free(submatrix);
            
            fq_nmod_mpoly_clear(cofactor, mpoly_ctx);
            fq_nmod_mpoly_clear(subdet, mpoly_ctx);
        }
    } else {
        // Sequential fallback for deeper levels or small matrices
        for (slong col = 0; col < size; col++) {
            if (fq_nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
                continue;
            }
            
            fq_nmod_mpoly_t cofactor, subdet;
            fq_nmod_mpoly_init(cofactor, mpoly_ctx);
            fq_nmod_mpoly_init(subdet, mpoly_ctx);
            
            // Create and fill submatrix
            fq_nmod_mpoly_t **submatrix = (fq_nmod_mpoly_t**) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (fq_nmod_mpoly_t*) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    fq_nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
                }
            }
            
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        fq_nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                        sub_j++;
                    }
                }
            }
            
            // Sequential computation
            compute_fq_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
            
            // Compute cofactor and store
            poly_mul_dense_optimized(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
            if (col % 2 == 0) {
                fq_nmod_mpoly_set(partial_results[col], cofactor, mpoly_ctx);
            } else {
                fq_nmod_mpoly_neg(partial_results[col], cofactor, mpoly_ctx);
            }
            
            // Cleanup
            for (slong i = 0; i < size-1; i++) {
                for (slong j = 0; j < size-1; j++) {
                    fq_nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
                }
                flint_free(submatrix[i]);
            }
            flint_free(submatrix);
            
            fq_nmod_mpoly_clear(cofactor, mpoly_ctx);
            fq_nmod_mpoly_clear(subdet, mpoly_ctx);
        }
    }
    
    // Sum results (sequential to avoid race conditions)
    fq_nmod_mpoly_t temp_sum;
    fq_nmod_mpoly_init(temp_sum, mpoly_ctx);
    
    for (slong col = 0; col < size; col++) {
        if (!fq_nmod_mpoly_is_zero(partial_results[col], mpoly_ctx)) {
            fq_nmod_mpoly_add(temp_sum, det_result, partial_results[col], mpoly_ctx);
            fq_nmod_mpoly_set(det_result, temp_sum, mpoly_ctx);
        }
        fq_nmod_mpoly_clear(partial_results[col], mpoly_ctx);
    }
    
    fq_nmod_mpoly_clear(temp_sum, mpoly_ctx);
    flint_free(partial_results);
}

void compute_fq_det_huang_interpolation(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    if (size <= 0) {
        fq_mvpoly_init(result, matrix[0][0].nvars, matrix[0][0].npars, matrix[0][0].ctx);
        return;
    }
    
    timing_info_t total_start = start_timing();
    
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    /*
    // Check if we're in a prime field
    if (!is_prime_field(ctx)) {
        printf("ERROR: sparse interpolation requires prime field\n");
        compute_fq_det_recursive(result, matrix, size);
        return;
    }
    */
    DET_PRINT("Computing %ldx%ld determinant via sparse interpolation\n", size, size);
    DET_PRINT("Variables: %ld, Parameters: %ld\n", nvars, npars);
    
    // Get the prime
    mp_limb_t p = fq_nmod_ctx_modulus(ctx)->mod.n;
    
    // Create nmod_mpoly context - FIX: Add modulus parameter
    nmod_mpoly_ctx_t nmod_ctx;
    nmod_mpoly_ctx_init(nmod_ctx, nvars + npars, ORD_LEX, p);  // Added p as 4th parameter
    
    // Convert matrix to poly_mat_t format for huang.h
    poly_mat_t huang_mat;
    poly_mat_init(&huang_mat, size, size, nmod_ctx);
    
    // Convert each entry
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            nmod_mpoly_t temp;
            nmod_mpoly_init(temp, nmod_ctx);
            
            // Convert fq_mvpoly to nmod_mpoly
            fq_mvpoly_to_nmod_mpoly(temp, &matrix[i][j], nmod_ctx);
            poly_mat_entry_set(&huang_mat, i, j, temp, nmod_ctx);
            
            nmod_mpoly_clear(temp, nmod_ctx);
        }
    }
    
    // Call Sparse.pdf-style determinant probing
    nmod_mpoly_t det_nmod;
    nmod_mpoly_init(det_nmod, nmod_ctx);
    
    timing_info_t huang_start = start_timing();
    ComputePolyMatrixDet(det_nmod, &huang_mat, nvars + npars, p, nmod_ctx);
    timing_info_t huang_elapsed = end_timing(huang_start);
    print_timing("sparse interpolation", huang_elapsed);
    
    // Convert result back to fq_mvpoly
    nmod_mpoly_to_fq_mvpoly(result, det_nmod, nvars, npars, nmod_ctx, ctx);
    
    DET_PRINT("Final result: %ld terms\n", result->nterms);
    
    // Cleanup
    poly_mat_clear(&huang_mat, nmod_ctx);
    nmod_mpoly_clear(det_nmod, nmod_ctx);
    nmod_mpoly_ctx_clear(nmod_ctx);
    
    timing_info_t total_elapsed = end_timing(total_start);
    print_timing("Total sparse interpolation method", total_elapsed);
}

/* Main function extracted from the #else branch */
void compute_fq_det_unified_interface(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    if (size <= 0) {
        fq_mvpoly_init(result, matrix[0][0].nvars, matrix[0][0].npars, matrix[0][0].ctx);
        return;
    }

    timing_info_t total_start = start_timing();
    
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    slong total_vars = nvars + npars;
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    slong max_threads = omp_get_max_threads();
    
    DET_PRINT("Computing %ldx%ld determinant (OpenMP: %ld threads available)\n", 
              size, size, max_threads);

    fq_mvpoly_init(result, nvars, npars, ctx);

    // Check for univariate optimization
    if (is_univariate_matrix(matrix, size) && size >= UNIVARIATE_THRESHOLD) {
        DET_PRINT("Detected univariate matrix, using specialized optimization\n");
        compute_fq_det_univariate_optimized(result, matrix, size);
        
        timing_info_t total_elapsed = end_timing(total_start);
        //print_timing("Total univariate computation", total_elapsed);
        return;
    }

    // ===== USE UNIFIED INTERFACE =====
    DET_PRINT("Using unified multivariate polynomial interface\n");
    
    // DEBUG: Print original matrix
    //debug_print_fq_mvpoly_matrix(matrix, size, "ORIGINAL");
    // Step 1: Create field context wrapper
    field_ctx_t field_ctx;
    field_ctx_init(&field_ctx, ctx);  // Properly initialize the field context

    // Verify the field type detection
    mp_limb_t p = fq_nmod_ctx_modulus(ctx)->mod.n;
    slong degree = fq_nmod_ctx_degree(ctx);
    DET_PRINT("Field: p=%lu, degree=%ld, detected type=%d\n", p, degree, field_ctx.field_id);
    
    // Debug: print field context details
    //printf("Field context details:\n");
    //printf("  field_id: %d\n", field_ctx.field_id);
    //printf("  modulus: ");
    //fq_nmod_ctx_modulus_print_pretty(ctx, "t");
    //printf("\n");
    //printf("  degree: %ld\n", degree);
    //printf("  characteristic: %lu\n", p);
    // Step 2: Create unified multivariate context
    unified_mpoly_ctx_t unified_ctx = unified_mpoly_ctx_init(total_vars, ORD_LEX, &field_ctx);
    if (!unified_ctx) {
        printf("ERROR: Failed to create unified context\n");
        return;
    }
    // Step 3: Allocate unified polynomial matrix
    unified_mpoly_t **unified_matrix = (unified_mpoly_t**) malloc(size * sizeof(unified_mpoly_t*));
    for (slong i = 0; i < size; i++) {
        unified_matrix[i] = (unified_mpoly_t*) malloc(size * sizeof(unified_mpoly_t));
        for (slong j = 0; j < size; j++) {
            unified_matrix[i][j] = unified_mpoly_init(unified_ctx);
            if (!unified_matrix[i][j]) {
                printf("ERROR: Failed to initialize unified polynomial at (%ld,%ld)\n", i, j);
                return;
            }
        }
    }
    // Step 4: Convert fq_mvpoly matrix to unified_mpoly matrix - Fast batch version
    timing_info_t convert_start = start_timing();
    DET_PRINT("Converting %ldx%ld matrix to unified format\n", size, size);
    
    // Create mpoly context for intermediate conversion
    fq_nmod_mpoly_ctx_t mpoly_ctx;
    fq_nmod_mpoly_ctx_init(mpoly_ctx, total_vars, ORD_LEX, ctx);
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            unified_mpoly_zero(unified_matrix[i][j]);
            
            fq_mvpoly_t *src_poly = &matrix[i][j];
            
            if (src_poly->nterms == 0) {
                continue; // Skip empty polynomials
            }
            
            // Method 1: Use push_term interface for batch processing
            fq_nmod_mpoly_t temp_poly;
            fq_nmod_mpoly_init(temp_poly, mpoly_ctx);
            
            // Push all terms at once
            for (slong t = 0; t < src_poly->nterms; t++) {
                if (fq_nmod_is_zero(src_poly->terms[t].coeff, ctx)) {
                    continue;
                }
                
                // Build combined exponent vector
                ulong *exp = (ulong*) flint_calloc(total_vars, sizeof(ulong));
                
                // Copy variable exponents
                if (src_poly->terms[t].var_exp && nvars > 0) {
                    for (slong v = 0; v < nvars; v++) {
                        exp[v] = (ulong)src_poly->terms[t].var_exp[v];
                    }
                }
                
                // Copy parameter exponents  
                if (src_poly->terms[t].par_exp && npars > 0) {
                    for (slong p = 0; p < npars; p++) {
                        exp[nvars + p] = (ulong)src_poly->terms[t].par_exp[p];
                    }
                }
                
                // Push term directly - much faster than repeated add
                fq_nmod_mpoly_push_term_fq_nmod_ui(temp_poly, src_poly->terms[t].coeff, exp, mpoly_ctx);
                
                flint_free(exp);
            }
            
            // Sort and combine like terms once at the end
            fq_nmod_mpoly_sort_terms(temp_poly, mpoly_ctx);
            fq_nmod_mpoly_combine_like_terms(temp_poly, mpoly_ctx);
            
            // Convert from fq_nmod_mpoly to unified_mpoly
            slong temp_length = fq_nmod_mpoly_length(temp_poly, mpoly_ctx);
            
            for (slong k = 0; k < temp_length; k++) {
                // Get coefficient
                fq_nmod_t coeff;
                fq_nmod_init(coeff, ctx);
                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, temp_poly, k, mpoly_ctx);
                
                // Get exponent vector
                ulong *exp = (ulong*) flint_calloc(total_vars, sizeof(ulong));
                fq_nmod_mpoly_get_term_exp_ui(exp, temp_poly, k, mpoly_ctx);
                
                // Convert coefficient to field element
                field_elem_u field_coeff;
                fq_nmod_to_field_elem(&field_coeff, coeff, &field_ctx);
                
                // Set in unified polynomial
                unified_mpoly_set_coeff_ui(unified_matrix[i][j], &field_coeff, exp);
                
                // Cleanup
                fq_nmod_clear(coeff, ctx);
                flint_free(exp);
            }
            
            fq_nmod_mpoly_clear(temp_poly, mpoly_ctx);
        }
    }
    
    fq_nmod_mpoly_ctx_clear(mpoly_ctx);
    timing_info_t convert_elapsed = end_timing(convert_start);
    // DEBUG: Print converted matrix
    //debug_print_unified_matrix(unified_matrix, size, "CONVERTED");
    // Step 5: Enable optimizations if applicable
    /*
    if (field_ctx.field_id == FIELD_ID_GF28) {
        unified_mpoly_enable_optimizations(FIELD_ID_GF28, 1);
        DET_PRINT("Enabled GF(2^8) optimizations\n");
    } else if (field_ctx.field_id == FIELD_ID_GF2128) {
        unified_mpoly_enable_optimizations(FIELD_ID_GF2128, 1);
        DET_PRINT("Enabled GF(2^128) optimizations\n");
    }
    */
    // Step 6: Compute determinant using unified interface
    unified_mpoly_t det_unified = unified_mpoly_init(unified_ctx);
    if (!det_unified) {
        printf("ERROR: Failed to initialize determinant polynomial\n");
        return;
    }

    timing_info_t det_start = start_timing();
    int use_parallel = (size >= PARALLEL_THRESHOLD && max_threads > 1);
    compute_unified_mpoly_det(det_unified, unified_matrix, size, unified_ctx, use_parallel);
    timing_info_t det_elapsed = end_timing(det_start);
    //print_timing("Determinant computation (unified)", det_elapsed);

    DET_PRINT("Unified determinant has %ld terms\n", unified_mpoly_length(det_unified));
    // Step 7: Convert result back to fq_mvpoly (with debugging)
    timing_info_t result_start = start_timing();
    fq_mvpoly_clear(result);  // Clear the initialization from the beginning
    fq_mvpoly_init(result, nvars, npars, ctx);
    
   // printf("\n--- Converting result back to fq_mvpoly ---\n");
    //printf("Field type for result conversion: %d\n", field_ctx.field_id);

    // Convert based on field type
    if (field_ctx.field_id == FIELD_ID_NMOD || is_prime_field(ctx)) {
        //printf("Using prime field result conversion\n");
        // For prime fields, convert from nmod_mpoly
        nmod_mpoly_struct *nmod_poly = GET_NMOD_POLY(det_unified);
        nmod_mpoly_ctx_struct *nmod_ctx = &(unified_ctx->ctx.nmod_ctx);

        // Use the existing conversion function
        slong nterms = nmod_mpoly_length(nmod_poly, nmod_ctx);
        DET_PRINT("Converting nmod_mpoly with %ld terms\n", nterms);

        // Pre-allocate the result polynomial to avoid reallocations
        if (result->alloc < nterms) {
            result->alloc = nterms + nterms/10; // Add 10% extra space
            result->terms = (fq_monomial_t*) flint_realloc(result->terms, 
                                                            result->alloc * sizeof(fq_monomial_t));
        }

        // Allocate a temporary exponent array once
        ulong *exp_ui = (ulong*) flint_malloc(total_vars * sizeof(ulong));

        // Batch convert all terms
        result->nterms = 0;
        for (slong i = 0; i < nterms; i++) {
            // Get coefficient
            mp_limb_t coeff_ui = nmod_mpoly_get_term_coeff_ui(nmod_poly, i, nmod_ctx);
            //printf("Result term %ld: coeff_ui = %lu\n", i, coeff_ui);

            // Get exponents
            nmod_mpoly_get_term_exp_ui(exp_ui, nmod_poly, i, nmod_ctx);
            /*
            printf("  exponents: ");
            for (slong k = 0; k < total_vars; k++) {
                printf("%lu ", exp_ui[k]);
            }
            printf("\n");
            */
            // Directly set the term without using fq_mvpoly_add_term
            fq_nmod_init(result->terms[result->nterms].coeff, ctx);
            fq_nmod_set_ui(result->terms[result->nterms].coeff, coeff_ui, ctx);

            // Split exponents
            if (nvars > 0) {
                result->terms[result->nterms].var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
                for (slong v = 0; v < nvars; v++) {
                    result->terms[result->nterms].var_exp[v] = (slong)exp_ui[v];
                }
            } else {
                result->terms[result->nterms].var_exp = NULL;
            }

            if (npars > 0) {
                result->terms[result->nterms].par_exp = (slong*) flint_calloc(npars, sizeof(slong));
                for (slong p = 0; p < npars; p++) {
                    result->terms[result->nterms].par_exp[p] = (slong)exp_ui[nvars + p];
                }
            } else {
                result->terms[result->nterms].par_exp = NULL;
            }

            result->nterms++;

            // Progress indicator for large conversions
            if (i > 0 && i % 10000 == 0) {
                DET_PRINT("Converted %ld/%ld terms...\n", i, nterms);
            }
        }

        flint_free(exp_ui);
        DET_PRINT("Conversion complete: %ld terms\n", result->nterms);

    } else if (field_ctx.field_id == FIELD_ID_FQ_ZECH) {
        //printf("Using Zech field result conversion\n");
        // For Zech logarithm fields, convert from fq_zech_mpoly
        fq_zech_mpoly_struct *zech_poly = GET_ZECH_POLY(det_unified);
        fq_zech_mpoly_ctx_struct *zech_ctx = &(unified_ctx->ctx.zech_ctx);

        slong nterms = fq_zech_mpoly_length(zech_poly, zech_ctx);
        DET_PRINT("Converting fq_zech_mpoly with %ld terms\n", nterms);

        if (nterms > 0) {
            // Pre-allocate the result polynomial
            if (result->alloc < nterms) {
                result->alloc = nterms + nterms/10; // Add 10% extra space
                result->terms = (fq_monomial_t*) flint_realloc(result->terms, 
                                                                result->alloc * sizeof(fq_monomial_t));
            }

            // Allocate temporary storage
            ulong *exp_ui = (ulong*) flint_malloc(total_vars * sizeof(ulong));
            fq_zech_t zech_coeff;
            fq_zech_init(zech_coeff, field_ctx.ctx.zech_ctx);

            // Convert each term
            result->nterms = 0;
            for (slong i = 0; i < nterms; i++) {
                // Get coefficient from Zech representation
                fq_zech_mpoly_get_term_coeff_fq_zech(zech_coeff, zech_poly, i, zech_ctx);

                // Get exponents
                fq_zech_mpoly_get_term_exp_ui(exp_ui, zech_poly, i, zech_ctx);

                // Initialize result coefficient
                fq_nmod_init(result->terms[result->nterms].coeff, ctx);

                // Convert Zech coefficient to fq_nmod
                fq_zech_get_fq_nmod(result->terms[result->nterms].coeff, zech_coeff, field_ctx.ctx.zech_ctx);

                // Split exponents
                if (nvars > 0) {
                    result->terms[result->nterms].var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
                    for (slong v = 0; v < nvars; v++) {
                        result->terms[result->nterms].var_exp[v] = (slong)exp_ui[v];
                    }
                } else {
                    result->terms[result->nterms].var_exp = NULL;
                }

                if (npars > 0) {
                    result->terms[result->nterms].par_exp = (slong*) flint_calloc(npars, sizeof(slong));
                    for (slong p = 0; p < npars; p++) {
                        result->terms[result->nterms].par_exp[p] = (slong)exp_ui[nvars + p];
                    }
                } else {
                    result->terms[result->nterms].par_exp = NULL;
                }

                result->nterms++;

                // Progress indicator for large conversions
                if (i > 0 && i % 10000 == 0) {
                    DET_PRINT("Converted %ld/%ld terms from Zech...\n", i, nterms);
                }
            }

            // Cleanup
            flint_free(exp_ui);
            fq_zech_clear(zech_coeff, field_ctx.ctx.zech_ctx);
            DET_PRINT("Zech conversion complete: %ld terms\n", result->nterms);
        }

    } else {
        //printf("Using extension field result conversion\n");
        // For extension fields, convert from fq_nmod_mpoly
        fq_nmod_mpoly_struct *fq_poly = GET_FQ_POLY(det_unified);
        fq_nmod_mpoly_ctx_struct *fq_ctx = &(unified_ctx->ctx.fq_ctx);

        // Use the existing conversion function if available
        slong nterms = fq_nmod_mpoly_length(fq_poly, fq_ctx);
        DET_PRINT("Converting fq_nmod_mpoly with %ld terms\n", nterms);

        if (nterms > 0) {
            // Use existing conversion function
            fq_nmod_mpoly_to_fq_mvpoly(result, fq_poly, nvars, npars, fq_ctx, ctx);
        }
    }
    
    timing_info_t result_elapsed = end_timing(result_start);
    //print_timing("Result conversion from unified format", result_elapsed);

    DET_PRINT("Final result: %ld terms\n", result->nterms);

    // DEBUG: Print final result
    //debug_print_fq_mvpoly_matrix(&result, 1, "FINAL RESULT");
    // Step 8: Cleanup
/*
    // Disable optimizations
    if (field_ctx.field_id == FIELD_ID_GF28) {
        unified_mpoly_enable_optimizations(FIELD_ID_GF28, 0);
    } else if (field_ctx.field_id == FIELD_ID_GF2128) {
        unified_mpoly_enable_optimizations(FIELD_ID_GF2128, 0);
    }
*/
    // Free unified matrix
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            unified_mpoly_clear(unified_matrix[i][j]);
        }
        free(unified_matrix[i]);
    }
    free(unified_matrix);

    // Clear unified polynomial and context
    unified_mpoly_clear(det_unified);
    unified_mpoly_ctx_clear(unified_ctx);
    
    timing_info_t total_elapsed = end_timing(total_start);
    (void) total_elapsed;
}
// ============= Main Interface with Algorithm Selection Implementation =============

void compute_fq_det_recursive_flint(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    if (size <= 0) {
        fq_mvpoly_init(result, matrix[0][0].nvars, matrix[0][0].npars, matrix[0][0].ctx);
        return;
    }
    
    timing_info_t total_start = start_timing();
    
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    
    // Choose algorithm based on configuration
    #if DET_ALGORITHM == DET_ALGORITHM_INTERPOLATION
    {
        printf("Using multivariate interpolation algorithm\n");
        
        // Include the interpolation header if not already included
        #ifndef FQ_NMOD_INTERPOLATION_OPTIMIZED_H
        #include "fq_multivariate_interpolation.h"
        #endif
        slong total_vars = nvars + npars;
        slong *var_bounds = (slong*) malloc(total_vars * sizeof(slong));
        compute_kronecker_bounds(var_bounds, matrix, size, nvars, npars);
        // Use interpolation algorithm
        fq_compute_det_by_interpolation_optimized(result, matrix, size, 
                                                 nvars, npars, ctx, var_bounds);
        return;
    }
    #elif DET_ALGORITHM == DET_ALGORITHM_KRONECKER
    {
        printf("Using Kronecker+HNF algorithm\n");
        compute_fq_det_kronecker(result, matrix, size);
        return;
    }
    #elif DET_ALGORITHM == DET_ALGORITHM_POLY_RECURSIVE
    {
        printf("Using polynomial recursive algorithm\n");
        compute_fq_det_poly_recursive(result, matrix, size);
        return;
    }
    #elif DET_ALGORITHM == DET_ALGORITHM_HUANG
    {
        // Check if we're in a prime field
        if (!is_prime_field(ctx)) {
            // Fall back to recursive algorithm (use code below)
        } else {
            compute_fq_det_huang_interpolation(result, matrix, size);
            return;
        }
    }
    #else
    {
        compute_fq_det_unified_interface(result, matrix, size);
    }
    #endif
}

// Compatibility interface
void compute_fq_det_recursive(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    compute_fq_det_recursive_flint(result, matrix, size);
}

static inline ulong reduce_exp_field_ui(ulong e, ulong q) {
    if (e == 0 || e < q) return e;
    return ((e - 1) % (q - 1)) + 1;
}

static ulong field_size_q_from_fq_ctx(const fq_nmod_ctx_t ctx) {
    mp_limb_t p = fq_nmod_ctx_prime(ctx);
    slong d = fq_nmod_ctx_degree(ctx);
    ulong q = 1;
    for (slong i = 0; i < d; i++) {
        if (q > WORD_MAX / p) return WORD_MAX;
        q *= p;
    }
    return q;
}

static void fq_nmod_poly_reduce_field_equation_inplace(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx) {
    if (!g_field_equation_reduction) return;
    slong deg = fq_nmod_poly_degree(poly, ctx);
    if (deg <= 0) return;
    ulong q = field_size_q_from_fq_ctx(ctx);
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

static void fq_nmod_mpoly_reduce_field_equation_inplace(fq_nmod_mpoly_t poly, const fq_nmod_mpoly_ctx_t ctx) {
    if (!g_field_equation_reduction) return;
    slong nterms = fq_nmod_mpoly_length(poly, ctx);
    if (nterms <= 0) return;
    ulong q = field_size_q_from_fq_ctx(ctx->fqctx);
    if (q <= 1 || q == WORD_MAX) return;

    slong nvars = fq_nmod_mpoly_ctx_nvars(ctx);
    ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
    fq_nmod_t coeff;
    fq_nmod_init(coeff, ctx->fqctx);
    fq_nmod_mpoly_t reduced;
    fq_nmod_mpoly_init(reduced, ctx);

    for (slong i = 0; i < nterms; i++) {
        fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, poly, i, ctx);
        fq_nmod_mpoly_get_term_exp_ui(exp, poly, i, ctx);
        for (slong k = 0; k < nvars; k++) exp[k] = reduce_exp_field_ui(exp[k], q);
        fq_nmod_mpoly_push_term_fq_nmod_ui(reduced, coeff, exp, ctx);
    }

    fq_nmod_mpoly_sort_terms(reduced, ctx);
    fq_nmod_mpoly_combine_like_terms(reduced, ctx);
    fq_nmod_mpoly_set(poly, reduced, ctx);
    fq_nmod_mpoly_clear(reduced, ctx);
    fq_nmod_clear(coeff, ctx->fqctx);
    flint_free(exp);
}
