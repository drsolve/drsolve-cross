/*
 * dixon_flint.c - Dixon Resultant Implementation for Finite Extension Fields
 *
 * This file contains the complete implementation of Dixon resultant computation
 * over finite fields using FLINT library for polynomial arithmetic and matrix operations.
 */

#include "dixon_flint.h"
#include <sys/time.h>

// Global method selection variable definitions
det_method_t dixon_global_method_step1 = -1;
det_method_t dixon_global_method_step4 = -1;
det_method_t dixon_global_method = -1; // deprecated compatibility alias

static const char *dixon_det_method_name(det_method_t method)
{
    switch (method) {
        case DET_METHOD_RECURSIVE:
            return "recursive expansion";
        case DET_METHOD_KRONECKER:
            return "Kronecker+HNF";
        case DET_METHOD_INTERPOLATION:
            return "interpolation";
        case DET_METHOD_HUANG:
            return "sparse interpolation";
        default:
            return "default";
    }
}

static void init_evaluation_parameters(fq_nmod_t *param_vals, slong npars,
                                      const fq_nmod_ctx_t ctx,
                                      slong attempt)
{
    for (slong i = 0; i < npars; i++) {
        fq_nmod_init(param_vals[i], ctx);
        fq_nmod_set_si(param_vals[i], 7 * (attempt + 1) * (i + 1) + 13, ctx);
        if (fq_nmod_is_zero(param_vals[i], ctx)) {
            fq_nmod_one(param_vals[i], ctx);
        }
    }
}

static void clear_evaluation_parameters(fq_nmod_t *param_vals, slong npars,
                                       const fq_nmod_ctx_t ctx)
{
    for (slong i = 0; i < npars; i++) {
        fq_nmod_clear(param_vals[i], ctx);
    }
    flint_free(param_vals);
}

static void init_extension_evaluation_parameters(fq_nmod_t *param_vals,
                                                slong npars,
                                                const fq_nmod_ctx_t ext_ctx,
                                                slong attempt)
{
    fq_nmod_t gen, constant;
    fq_nmod_init(gen, ext_ctx);
    fq_nmod_init(constant, ext_ctx);
    fq_nmod_gen(gen, ext_ctx);

    for (slong i = 0; i < npars; i++) {
        fq_nmod_init(param_vals[i], ext_ctx);
        fq_nmod_set(param_vals[i], gen, ext_ctx);
        fq_nmod_set_ui(constant, 7 * (attempt + 1) * (i + 1) + 13, ext_ctx);
        fq_nmod_add(param_vals[i], param_vals[i], constant, ext_ctx);
        if (fq_nmod_is_zero(param_vals[i], ext_ctx)) {
            fq_nmod_set(param_vals[i], gen, ext_ctx);
        }
    }

    fq_nmod_clear(gen, ext_ctx);
    fq_nmod_clear(constant, ext_ctx);
}

static void evaluate_fq_mvpoly_at_extension_params(fq_nmod_t result,
                                                  const fq_mvpoly_t *poly,
                                                  const fq_nmod_t *param_vals,
                                                  const fq_nmod_ctx_t ext_ctx)
{
    fq_nmod_t term_val, coeff_val;
    fq_nmod_init(term_val, ext_ctx);
    fq_nmod_init(coeff_val, ext_ctx);
    fq_nmod_zero(result, ext_ctx);

    if (poly == NULL || poly->nterms == 0) {
        fq_nmod_clear(term_val, ext_ctx);
        fq_nmod_clear(coeff_val, ext_ctx);
        return;
    }

    for (slong i = 0; i < poly->nterms; i++) {
        fq_nmod_set_ui(coeff_val, nmod_poly_get_coeff_ui(poly->terms[i].coeff, 0), ext_ctx);
        fq_nmod_set(term_val, coeff_val, ext_ctx);

        if (poly->terms[i].par_exp && poly->npars > 0) {
            for (slong j = 0; j < poly->npars; j++) {
                slong exp = poly->terms[i].par_exp[j];
                if (exp > 0) {
                    fq_nmod_t pow_val;
                    fq_nmod_init(pow_val, ext_ctx);
                    fq_nmod_pow_ui(pow_val, param_vals[j], exp, ext_ctx);
                    fq_nmod_mul(term_val, term_val, pow_val, ext_ctx);
                    fq_nmod_clear(pow_val, ext_ctx);
                }
            }
        }

        fq_nmod_add(result, result, term_val, ext_ctx);
    }

    fq_nmod_clear(term_val, ext_ctx);
    fq_nmod_clear(coeff_val, ext_ctx);
}

static slong evaluate_selected_submatrix_rank(fq_mvpoly_t ***full_matrix,
                                             const slong *row_indices,
                                             const slong *col_indices,
                                             slong size,
                                             const fq_nmod_t *param_vals,
                                             const fq_nmod_ctx_t ctx)
{
    fq_nmod_mat_t mat;
    fq_nmod_t value;
    fq_nmod_mat_init(mat, size, size, ctx);
    fq_nmod_init(value, ctx);

    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_mvpoly_t *entry = full_matrix[row_indices[i]][col_indices[j]];
            if (entry != NULL) {
                evaluate_fq_mvpoly_at_params(value, entry, param_vals);
                fq_nmod_set(fq_nmod_mat_entry(mat, i, j), value, ctx);
            } else {
                fq_nmod_zero(fq_nmod_mat_entry(mat, i, j), ctx);
            }
        }
    }

    slong rank = fq_nmod_mat_rank(mat, ctx);
    fq_nmod_clear(value, ctx);
    fq_nmod_mat_clear(mat, ctx);
    return rank;
}

static slong evaluate_selected_submatrix_rank_extension(fq_mvpoly_t ***full_matrix,
                                                       const slong *row_indices,
                                                       const slong *col_indices,
                                                       slong size,
                                                       const fq_nmod_t *param_vals,
                                                       const fq_nmod_ctx_t ext_ctx)
{
    fq_nmod_mat_t mat;
    fq_nmod_t value;
    fq_nmod_mat_init(mat, size, size, ext_ctx);
    fq_nmod_init(value, ext_ctx);

    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_mvpoly_t *entry = full_matrix[row_indices[i]][col_indices[j]];
            if (entry != NULL) {
                evaluate_fq_mvpoly_at_extension_params(value, entry, param_vals, ext_ctx);
                fq_nmod_set(fq_nmod_mat_entry(mat, i, j), value, ext_ctx);
            } else {
                fq_nmod_zero(fq_nmod_mat_entry(mat, i, j), ext_ctx);
            }
        }
    }

    slong rank = fq_nmod_mat_rank(mat, ext_ctx);
    fq_nmod_clear(value, ext_ctx);
    fq_nmod_mat_clear(mat, ext_ctx);
    return rank;
}

// Build cancellation matrix in multivariate form
void build_fq_cancellation_matrix_mvpoly(fq_mvpoly_t ***M, fq_mvpoly_t *polys, 
                                        slong nvars, slong npars) {
    slong n = nvars + 1;
    
    // Allocate matrix space
    *M = (fq_mvpoly_t**) flint_malloc(n * sizeof(fq_mvpoly_t*));
    for (slong i = 0; i < n; i++) {
        (*M)[i] = (fq_mvpoly_t*) flint_malloc(n * sizeof(fq_mvpoly_t));
    }
    
    // Build matrix entries
    for (slong i = 0; i < n; i++) {
        for (slong j = 0; j < n; j++) {
            fq_mvpoly_init(&(*M)[i][j], 2 * nvars, npars, polys[0].ctx);
            
            // Substitute variables according to row index
            for (slong t = 0; t < polys[j].nterms; t++) {
                slong *new_var_exp = (slong*) flint_calloc(2 * nvars, sizeof(slong));
                
                for (slong k = 0; k < nvars; k++) {
                    slong orig_exp = polys[j].terms[t].var_exp ? polys[j].terms[t].var_exp[k] : 0;
                    
                    if (k < i) {
                        // Use dual variable ~x_k
                        new_var_exp[nvars + k] = orig_exp;
                    } else {
                        // Use original variable x_k
                        new_var_exp[k] = orig_exp;
                    }
                }
                
                fq_mvpoly_add_term_fast(&(*M)[i][j], new_var_exp, polys[j].terms[t].par_exp, 
                                  polys[j].terms[t].coeff);
                flint_free(new_var_exp);
            }
        }
    }
}

void perform_fq_matrix_row_operations_mvpoly(fq_mvpoly_t ***new_matrix, fq_mvpoly_t ***original_matrix,
                                                   slong nvars, slong npars) {
    slong n = nvars + 1;
    
    *new_matrix = (fq_mvpoly_t**) flint_malloc(n * sizeof(fq_mvpoly_t*));
    for (slong i = 0; i < n; i++) {
        (*new_matrix)[i] = (fq_mvpoly_t*) flint_malloc(n * sizeof(fq_mvpoly_t));
    }
    
    for (slong j = 0; j < n; j++) {
        fq_mvpoly_copy(&(*new_matrix)[0][j], &(*original_matrix)[0][j]);
    }
    
    for (slong i = 1; i < n; i++) {
        for (slong j = 0; j < n; j++) {
            fq_mvpoly_t diff;
            fq_mvpoly_sub(&diff, &(*original_matrix)[i][j], &(*original_matrix)[i-1][j]);
      
            if (diff.nterms > 0) {

                divide_by_fq_linear_factor_flint(&(*new_matrix)[i][j], &diff, 
                                               i-1, 2*nvars, npars);
            } else {
                fq_mvpoly_init(&(*new_matrix)[i][j], 2*nvars, npars, diff.ctx);
            }
            
            fq_mvpoly_clear(&diff);
        }
    }
}

// Compute Dixon resultant degree bound
slong compute_fq_dixon_resultant_degree_bound(fq_mvpoly_t *polys, slong npolys, slong nvars, slong npars) {
    slong degree_product = 1;
    
    for (slong i = 0; i < npolys; i++) {
        slong max_total_deg = 0;
        
        // Find maximum total degree of polynomial i
        for (slong t = 0; t < polys[i].nterms; t++) {
            slong total_deg = 0;
            
            // Sum variable degrees
            if (polys[i].terms[t].var_exp) {
                for (slong j = 0; j < nvars; j++) {
                    total_deg += polys[i].terms[t].var_exp[j];
                }
            }
            
            // Sum parameter degrees
            if (polys[i].terms[t].par_exp && npars > 0) {
                for (slong j = 0; j < npars; j++) {
                    total_deg += polys[i].terms[t].par_exp[j];
                }
            }
            
            if (total_deg > max_total_deg) {
                max_total_deg = total_deg;
            }
        }
        
        degree_product *= max_total_deg;
    }
    
    return degree_product + 1;
}

void compute_fq_coefficient_matrix_det(fq_mvpoly_t *result, fq_mvpoly_t **coeff_matrix,
                                       slong size, slong npars, const fq_nmod_ctx_t ctx,
                                       det_method_t method, slong res_deg_bound) {
    if (size == 0) {
        fq_mvpoly_init(result, 0, npars, ctx);
        return;
    }
    
    if (npars == 0) {
        fq_mvpoly_init(result, 0, npars, ctx);
        
        fq_nmod_mat_t scalar_mat;
        fq_nmod_mat_init(scalar_mat, size, size, ctx);
        
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                if (coeff_matrix[i][j].nterms > 0) {
                    fq_nmod_set(fq_nmod_mat_entry(scalar_mat, i, j), 
                                coeff_matrix[i][j].terms[0].coeff, ctx);
                } else {
                    fq_nmod_zero(fq_nmod_mat_entry(scalar_mat, i, j), ctx);
                }
            }
        }
        
        printf("\nComputing Resultant\n");
        clock_t cpu_start = clock();
        double wall_start = get_wall_time();
        
        fq_nmod_t det;
        fq_nmod_init(det, ctx);
        fq_nmod_mat_det(det, scalar_mat, ctx);
        
        clock_t cpu_end = clock();
        double wall_end = get_wall_time();
        double cpu_elapsed = (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;
        double wall_elapsed = wall_end - wall_start;
        
        int threads = 1;
        #ifdef _OPENMP
        threads = omp_get_max_threads();
        #endif
        //printf("CPU time: %.3f seconds | Wall time: %.3f seconds | Threads: %d\n", cpu_elapsed, wall_elapsed, threads);
        
        if (!fq_nmod_is_zero(det, ctx)) {
            fq_mvpoly_add_term_fast(result, NULL, NULL, det);
        }
        
        fq_nmod_clear(det, ctx);
        fq_nmod_mat_clear(scalar_mat, ctx);
        
    } else if (npars == 1) {
        clock_t cpu_start = clock();
        double wall_start = get_wall_time();
        
        if (method == DET_METHOD_INTERPOLATION) {
            printf("Method: interpolation\n");
            
            fq_compute_det_by_interpolation(result, coeff_matrix, size,
                                           0, npars, ctx, res_deg_bound);
        } else {
            fq_mvpoly_init(result, 0, npars, ctx);
            
            fq_nmod_poly_mat_t poly_mat;
            fq_nmod_poly_mat_init(poly_mat, size, size, ctx);
            
            for (slong i = 0; i < size; i++) {
                for (slong j = 0; j < size; j++) {
                    fq_nmod_poly_struct *entry = fq_nmod_poly_mat_entry(poly_mat, i, j);
                    fq_nmod_poly_zero(entry, ctx);
                    
                    for (slong t = 0; t < coeff_matrix[i][j].nterms; t++) {
                        slong deg = coeff_matrix[i][j].terms[t].par_exp ? 
                                   coeff_matrix[i][j].terms[t].par_exp[0] : 0;
                        fq_nmod_poly_set_coeff(entry, deg, 
                                              coeff_matrix[i][j].terms[t].coeff, ctx);
                    }
                }
            }
            
            fq_nmod_poly_t det_poly;
            fq_nmod_poly_init(det_poly, ctx);
            
            printf("Method: HNF\n");
            
            fq_nmod_poly_mat_det_iter(det_poly, poly_mat, ctx);
            
            slong det_deg = fq_nmod_poly_degree(det_poly, ctx);
            if (det_deg >= 0) {
                for (slong i = 0; i <= det_deg; i++) {
                    fq_nmod_t coeff;
                    fq_nmod_init(coeff, ctx);
                    fq_nmod_poly_get_coeff(coeff, det_poly, i, ctx);
                    if (!fq_nmod_is_zero(coeff, ctx)) {
                        slong par_exp[1] = {i};
                        fq_mvpoly_add_term_fast(result, NULL, par_exp, coeff);
                    }
                    fq_nmod_clear(coeff, ctx);
                }
            }
            
            fq_nmod_poly_clear(det_poly, ctx);
            fq_nmod_poly_mat_clear(poly_mat, ctx);
        }
        
        clock_t cpu_end = clock();
        double wall_end = get_wall_time();
        double cpu_elapsed = (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;
        double wall_elapsed = wall_end - wall_start;
        
        int threads = 1;
        #ifdef _OPENMP
        threads = omp_get_max_threads();
        #endif
        // printf("CPU time: %.3f seconds | Wall time: %.3f seconds | Threads: %d\n", cpu_elapsed, wall_elapsed, threads);
        
    } else {
        clock_t cpu_start = clock();
        double wall_start = get_wall_time();
        switch (method) {
            case DET_METHOD_INTERPOLATION:
                printf("Method: interpolation\n");
                
                fq_compute_det_by_interpolation(result, coeff_matrix, size,
                                               0, npars, ctx, res_deg_bound);
                break;
                
            case DET_METHOD_RECURSIVE:
                printf("Method: recursive expansion\n");
                
                compute_fq_det_recursive(result, coeff_matrix, size);
                break;
                
            case DET_METHOD_KRONECKER:
                printf("Method: Kronecker+HNF\n");
                
                compute_fq_det_kronecker(result, coeff_matrix, size);
                break;

            case DET_METHOD_HUANG:
                printf("Method: sparse interpolation\n");
                
                compute_fq_det_huang_interpolation(result, coeff_matrix, size);
                break;
                
            default:
                printf("Method: interpolation (default)\n");
                fq_compute_det_by_interpolation(result, coeff_matrix, size,
                                               0, npars, ctx, res_deg_bound);
                break;
        }
        clock_t cpu_end = clock();
        double wall_end = get_wall_time();
        double cpu_elapsed = (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;
        double wall_elapsed = wall_end - wall_start;
        
        int threads = 1;
        #ifdef _OPENMP
        threads = omp_get_max_threads();
        #endif
        // printf("CPU time: %.3f seconds | Wall time: %.3f seconds | Threads: %d\n", cpu_elapsed, wall_elapsed, threads);
    }

    if (g_field_equation_reduction) {
        fq_mvpoly_reduce_field_equation(result);
    }
}

// Extended tracker structure with pre-allocated workspace
// Initialize optimized tracker
static void unified_row_basis_tracker_init(unified_row_basis_tracker_t *tracker, 
                                            slong max_size, slong ncols, 
                                            field_ctx_t *ctx) {
    tracker->max_size = max_size;
    tracker->ncols = ncols;
    tracker->ctx = ctx;
    tracker->current_rank = 0;
    tracker->initialized = 1;
    tracker->workspace_initialized = 0;
    
    void *ctx_ptr = (ctx->field_id == FIELD_ID_NMOD) ? 
                   (void*)&ctx->ctx.nmod_ctx : 
                   (void*)ctx->ctx.fq_ctx;
    
    // Allocate main storage
    tracker->reduced_rows = (field_elem_u*) flint_calloc(max_size * ncols, sizeof(field_elem_u));
    tracker->pivot_cols = (slong*) flint_calloc(max_size, sizeof(slong));
    tracker->selected_indices = (slong*) flint_calloc(max_size, sizeof(slong));
    
    // Pre-allocate workspace
    tracker->work_row = (field_elem_u*) flint_malloc(ncols * sizeof(field_elem_u));
    tracker->temp_vars = (field_elem_u*) flint_malloc(4 * sizeof(field_elem_u)); // factor, temp, pivot_val, neg_temp
    
    // Initialize all field elements
    for (slong i = 0; i < max_size * ncols; i++) {
        field_init_elem(&tracker->reduced_rows[i], ctx->field_id, ctx_ptr);
        field_set_zero(&tracker->reduced_rows[i], ctx->field_id, ctx_ptr);
    }
    
    // Initialize workspace
    for (slong j = 0; j < ncols; j++) {
        field_init_elem(&tracker->work_row[j], ctx->field_id, ctx_ptr);
    }
    for (slong i = 0; i < 4; i++) {
        field_init_elem(&tracker->temp_vars[i], ctx->field_id, ctx_ptr);
    }
    tracker->workspace_initialized = 1;
    
    // Initialize pivot columns to -1
    for (slong i = 0; i < max_size; i++) {
        tracker->pivot_cols[i] = -1;
        tracker->selected_indices[i] = -1;
    }
}

// Clear optimized tracker
static void unified_row_basis_tracker_clear(unified_row_basis_tracker_t *tracker) {
    if (!tracker->initialized) return;
    
    void *ctx_ptr = (tracker->ctx->field_id == FIELD_ID_NMOD) ? 
                   (void*)&tracker->ctx->ctx.nmod_ctx : 
                   (void*)tracker->ctx->ctx.fq_ctx;
    
    // Clear main storage
    if (tracker->reduced_rows) {
        for (slong i = 0; i < tracker->max_size * tracker->ncols; i++) {
            field_clear_elem(&tracker->reduced_rows[i], tracker->ctx->field_id, ctx_ptr);
        }
        flint_free(tracker->reduced_rows);
    }
    
    // Clear workspace
    if (tracker->workspace_initialized) {
        for (slong j = 0; j < tracker->ncols; j++) {
            field_clear_elem(&tracker->work_row[j], tracker->ctx->field_id, ctx_ptr);
        }
        for (slong i = 0; i < 4; i++) {
            field_clear_elem(&tracker->temp_vars[i], tracker->ctx->field_id, ctx_ptr);
        }
        flint_free(tracker->work_row);
        flint_free(tracker->temp_vars);
    }
    
    if (tracker->pivot_cols) flint_free(tracker->pivot_cols);
    if (tracker->selected_indices) flint_free(tracker->selected_indices);
    
    tracker->initialized = 0;
    tracker->workspace_initialized = 0;
}

// Optimized version of adding row to basis - mathematical logic completely unchanged
static int unified_try_add_row_to_basis(unified_row_basis_tracker_t *tracker, 
                                         const field_elem_u *unified_mat,
                                         slong new_row_idx, slong ncols) {
    if (!tracker->initialized || tracker->current_rank >= tracker->max_size) {
        return 0;
    }
    
    void *ctx_ptr = (tracker->ctx->field_id == FIELD_ID_NMOD) ? 
                   (void*)&tracker->ctx->ctx.nmod_ctx : 
                   (void*)tracker->ctx->ctx.fq_ctx;
    
    // Use pre-allocated work row to avoid allocation each time
    field_elem_u *work_row = tracker->work_row;
    
    // Copy input row to work row (reuse allocated space)
    for (slong j = 0; j < ncols; j++) {
        field_set_elem(&work_row[j], &unified_mat[new_row_idx * ncols + j], 
                      tracker->ctx->field_id, ctx_ptr);
    }
    
    // Use pre-allocated temporary variables
    field_elem_u *factor = &tracker->temp_vars[0];
    field_elem_u *temp = &tracker->temp_vars[1];
    field_elem_u *pivot_val = &tracker->temp_vars[2];
    field_elem_u *neg_temp = &tracker->temp_vars[3];
    
    // Perform elimination for each existing basis vector (mathematical logic unchanged)
    for (slong i = 0; i < tracker->current_rank; i++) {
        slong pivot_col = tracker->pivot_cols[i];
        
        // Validate pivot column index
        if (pivot_col < 0 || pivot_col >= ncols) continue;
        
        // If work row is non-zero at pivot position, perform elimination
        if (!field_is_zero(&work_row[pivot_col], tracker->ctx->field_id, ctx_ptr)) {
            // Get pivot value of basis vector
            slong base_idx = i * ncols + pivot_col;
            field_set_elem(pivot_val, &tracker->reduced_rows[base_idx], 
                          tracker->ctx->field_id, ctx_ptr);
            
            // Calculate elimination factor = work_row[pivot_col] / pivot_val
            field_inv(temp, pivot_val, tracker->ctx->field_id, ctx_ptr);
            field_mul(factor, &work_row[pivot_col], temp, tracker->ctx->field_id, ctx_ptr);
            
            // Perform elimination: work_row -= factor * basis_row[i]
            // Optimization: pre-determine if neg_temp is needed to avoid conditional branches in inner loop
            int use_direct_add = (tracker->ctx->field_id >= FIELD_ID_GF28 && 
                                 tracker->ctx->field_id <= FIELD_ID_GF2128);
            
            for (slong j = 0; j < ncols; j++) {
                slong idx = i * ncols + j;
                field_mul(temp, factor, &tracker->reduced_rows[idx], 
                         tracker->ctx->field_id, ctx_ptr);
                
                if (use_direct_add) {
                    // For GF(2^n), subtraction equals addition
                    field_add(&work_row[j], &work_row[j], temp, 
                             tracker->ctx->field_id, ctx_ptr);
                } else {
                    // For other fields, use actual subtraction
                    field_neg(neg_temp, temp, tracker->ctx->field_id, ctx_ptr);
                    field_add(&work_row[j], &work_row[j], neg_temp, 
                             tracker->ctx->field_id, ctx_ptr);
                }
            }
        }
    }
    
    // Find first non-zero position (mathematical logic unchanged)
    slong first_nonzero = -1;
    for (slong j = 0; j < ncols; j++) {
        if (!field_is_zero(&work_row[j], tracker->ctx->field_id, ctx_ptr)) {
            first_nonzero = j;
            break;
        }
    }
    
    // If all zeros, then linearly dependent
    if (first_nonzero == -1) {
        return 0;  // Work row doesn't need cleanup because it's pre-allocated
    }
    
    // Normalize work row (make first non-zero element 1)
    field_inv(temp, &work_row[first_nonzero], tracker->ctx->field_id, ctx_ptr);
    
    // Store normalized row to basis
    slong base_row_start = tracker->current_rank * ncols;
    for (slong j = 0; j < ncols; j++) {
        slong idx = base_row_start + j;
        if (j < first_nonzero) {
            // Positions before pivot should be 0
            field_set_zero(&tracker->reduced_rows[idx], tracker->ctx->field_id, ctx_ptr);
        } else {
            // Normalize and store
            field_mul(&tracker->reduced_rows[idx], &work_row[j], temp, 
                     tracker->ctx->field_id, ctx_ptr);
        }
    }
    
    // Update tracking information (mathematical logic unchanged)
    tracker->pivot_cols[tracker->current_rank] = first_nonzero;
    tracker->selected_indices[tracker->current_rank] = new_row_idx;
    tracker->current_rank++;
    
    return 1;
}

// Helper function to compute maximum degree (not total degree) of a polynomial
static slong compute_fq_polynomial_total_degree(fq_mvpoly_t *poly, slong npars) {
    if (poly == NULL || poly->nterms == 0) {
        return 0;
    }
    
    slong max_degree = 0;
    
    for (slong t = 0; t < poly->nterms; t++) {
        // Check variable degrees
        if (poly->terms[t].var_exp) {
            for (slong v = 0; v < poly->nvars; v++) {
                if (poly->terms[t].var_exp[v] > max_degree) {
                    max_degree = poly->terms[t].var_exp[v];
                }
            }
        }
        
        // Check parameter degrees
        if (poly->terms[t].par_exp && npars > 0) {
            for (slong p = 0; p < npars; p++) {
                if (poly->terms[t].par_exp[p] > max_degree) {
                    max_degree = poly->terms[t].par_exp[p];
                }
            }
        }
    }
    
    return max_degree;
}

// Compute row maximum total degree
static slong compute_fq_row_max_total_degree(fq_mvpoly_t **matrix_row, slong ncols, slong npars) {
    slong max_degree = -1;
    
    for (slong j = 0; j < ncols; j++) {
        slong poly_deg = compute_fq_polynomial_total_degree(matrix_row[j], npars);
        if (poly_deg > max_degree) {
            max_degree = poly_deg;
        }
    }
    
    return max_degree;
}

// Compute column maximum total degree
static slong compute_fq_col_max_total_degree(fq_mvpoly_t ***matrix, slong col_idx, slong nrows, slong npars) {
    slong max_degree = -1;
    
    for (slong i = 0; i < nrows; i++) {
        slong poly_deg = compute_fq_polynomial_total_degree(matrix[i][col_idx], npars);
        if (poly_deg > max_degree) {
            max_degree = poly_deg;
        }
    }
    
    return max_degree;
}


// Compute maximum total degree of a column in selected rows submatrix
static slong compute_fq_selected_rows_col_max_total_degree(fq_mvpoly_t ***full_matrix, 
                                                           slong *selected_rows, 
                                                           slong num_selected_rows,
                                                           slong col_idx, 
                                                           slong npars) {
    slong max_degree = -1;
    
    for (slong i = 0; i < num_selected_rows; i++) {
        slong row_idx = selected_rows[i];
        slong poly_deg = compute_fq_polynomial_total_degree(full_matrix[row_idx][col_idx], npars);
        if (poly_deg > max_degree) {
            max_degree = poly_deg;
        }
    }
    
    return max_degree;
}

// Compute maximum total degree of a row in selected columns submatrix
static slong compute_fq_selected_cols_row_max_total_degree(fq_mvpoly_t ***full_matrix, 
                                                           slong row_idx,
                                                           slong *selected_cols, 
                                                           slong num_selected_cols,
                                                           slong npars) {
    slong max_degree = -1;
    
    for (slong j = 0; j < num_selected_cols; j++) {
        slong col_idx = selected_cols[j];
        slong poly_deg = compute_fq_polynomial_total_degree(full_matrix[row_idx][col_idx], npars);
        if (poly_deg > max_degree) {
            max_degree = poly_deg;
        }
    }
    
    return max_degree;
}

// Check if two index arrays are identical
static int indices_equal(slong *indices1, slong *indices2, slong size) {
    for (slong i = 0; i < size; i++) {
        if (indices1[i] != indices2[i]) {
            return 0;
        }
    }
    return 1;
}
// Find pivot rows for maximal rank submatrix - ensure linear independence
// More efficient version: using incremental rank checking
void find_pivot_rows_nmod_fixed(slong **selected_rows_out, slong *num_selected,
                                const nmod_mat_t mat) {
    
    slong nrows = nmod_mat_nrows(mat);
    slong ncols = nmod_mat_ncols(mat);
    nmod_t mod = mat->mod;
    
    if (nrows == 0 || ncols == 0) {
        *selected_rows_out = NULL;
        *num_selected = 0;
        return;
    }
    
    slong min_dim = FLINT_MIN(nrows, ncols);
    slong *selected_rows = (slong*) flint_malloc(min_dim * sizeof(slong));
    slong rank = 0;
    
    // Create working copy
    nmod_mat_t A;
    nmod_mat_init(A, nrows, ncols, mod.n);
    nmod_mat_set(A, mat);
    
    // Row permutation tracking
    slong *P = (slong*) flint_malloc(nrows * sizeof(slong));
    for (slong i = 0; i < nrows; i++) {
        P[i] = i;
    }
    
    // PROPER LU DECOMPOSITION using FLINT-style algorithm
    slong current_row = 0;
    for (slong col = 0; col < ncols && current_row < nrows; col++) {
        // printf("%d %d\n",col,ncols);
        // Find pivot in current column
        slong pivot_row = -1;
        for (slong i = current_row; i < nrows; i++) {
            if (nmod_mat_entry(A, i, col) != 0) {
                pivot_row = i;
                break;
            }
        }
        
        // No pivot found in this column
        if (pivot_row == -1) {
            continue;
        }
        
        // Record this pivot row
        selected_rows[rank] = P[pivot_row];
        rank++;
        
        // Swap rows if needed - use FLINT's efficient method
        if (pivot_row != current_row) {
            // Swap permutation
            slong temp_idx = P[current_row];
            P[current_row] = P[pivot_row];
            P[pivot_row] = temp_idx;
            
            // Swap matrix rows efficiently
            //nn_ptr temp_row = A->rows[current_row];
            //A->rows[current_row] = A->rows[pivot_row];
            //A->rows[pivot_row] = temp_row;
            // Swap matrix rows using compatible API
			
            for (slong j = 0; j < ncols; j++) {
                mp_limb_t temp_val = nmod_mat_entry(A, current_row, j);
                nmod_mat_entry(A, current_row, j) = nmod_mat_entry(A, pivot_row, j);
                nmod_mat_entry(A, pivot_row, j) = temp_val;
            }
        }
        
        // Eliminate below pivot using FLINT vectorized operations
        mp_limb_t pivot = nmod_mat_entry(A, current_row, col);
        mp_limb_t pivot_inv = n_invmod(pivot, mod.n);
        
        for (slong i = current_row + 1; i < nrows; i++) {
            mp_limb_t factor = nmod_mat_entry(A, i, col);
            if (factor == 0) continue;
            
            factor = n_mulmod2_preinv(factor, pivot_inv, mod.n, mod.ninv);
            mp_limb_t neg_factor = nmod_neg(factor, mod);
            
            // Use FLINT's vectorized subtraction
            /*
            slong remaining_cols = ncols - col;
            if (remaining_cols > 0) {
                _nmod_vec_scalar_addmul_nmod(A->rows[i] + col,
                                           A->rows[current_row] + col,
                                           remaining_cols, neg_factor, mod);
            }
            */
            mp_limb_t *row_i = &nmod_mat_entry(A, i, 0);
            mp_limb_t *row_current = &nmod_mat_entry(A, current_row, 0);
            if (ncols - col > 0) {
                _nmod_vec_scalar_addmul_nmod(row_i + col, 
                                           row_current + col, 
                                           ncols - col, 
                                           neg_factor, mod);
            }
        }
        
        current_row++;
    }
    
    // Set output
    *num_selected = rank;
    if (rank > 0) {
        *selected_rows_out = (slong*) flint_realloc(selected_rows, rank * sizeof(slong));
    } else {
        flint_free(selected_rows);
        *selected_rows_out = NULL;
    }
    
    // Cleanup
    nmod_mat_clear(A);
    flint_free(P);
}
// ============================================================================
// Optimized find_pivot_rows_simple - directly calls nmod version
// ============================================================================
void find_pivot_rows_simple(slong **selected_rows_out, slong *num_selected,
                                        const field_elem_u *unified_mat, 
                                        slong nrows, slong ncols,
                                        field_ctx_t *ctx) {
    
    // ============================================================================
    // Prime field fast path: directly use nmod_fixed implementation (optimal performance)
    // ============================================================================
    if (ctx->field_id == FIELD_ID_NMOD) {
        //printf("Using direct nmod_fixed implementation for prime field (optimal performance)\n");
        clock_t start = clock();
        
        // Convert to nmod_mat format
        nmod_mat_t nmod_mat;
        nmod_mat_init(nmod_mat, nrows, ncols, ctx->ctx.nmod_ctx.n);
        
        // Efficient data copy (direct access to nmod field)
        for (slong i = 0; i < nrows; i++) {
            for (slong j = 0; j < ncols; j++) {
                nmod_mat_entry(nmod_mat, i, j) = unified_mat[i * ncols + j].nmod;
            }
        }
        
        // Directly use nmod_fixed version (avoid adaptive layer selection overhead)
        find_pivot_rows_nmod_fixed(selected_rows_out, num_selected, nmod_mat);
        
        clock_t end = clock();
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        
        nmod_mat_clear(nmod_mat);
        //printf("Prime field computation completed in %.4f seconds using direct nmod path\n", elapsed);
        return;
    }
    
    // ============================================================================
    // Non-prime field: use unified interface implementation
    // ============================================================================
    //printf("Using unified interface for non-prime field computation\n");
    
    void *ctx_ptr = (ctx->field_id == FIELD_ID_FQ_ZECH) ? 
                   (void*)ctx->ctx.zech_ctx : 
                   (void*)ctx->ctx.fq_ctx;
    
    slong max_rank = FLINT_MIN(nrows, ncols);
    slong *selected_rows = (slong*) flint_malloc(max_rank * sizeof(slong));
    slong rank = 0;
    
    // Lightweight linear independence checker
    typedef struct {
        field_elem_u *rows;     
        slong *pivot_positions; 
        slong count;            
        slong ncols;
    } basis_tracker_t;
    
    basis_tracker_t tracker;
    tracker.rows = (field_elem_u*) flint_calloc(max_rank * ncols, sizeof(field_elem_u));
    tracker.pivot_positions = (slong*) flint_malloc(max_rank * sizeof(slong));
    tracker.count = 0;
    tracker.ncols = ncols;
    
    // Initialize all field elements
    for (slong i = 0; i < max_rank * ncols; i++) {
        field_init_elem(&tracker.rows[i], ctx->field_id, ctx_ptr);
    }
    
    clock_t start = clock();
    
    // Process row by row
    for (slong row = 0; row < nrows && rank < max_rank; row++) {
        // Create test vector
        field_elem_u *test_vec = (field_elem_u*) flint_malloc(ncols * sizeof(field_elem_u));
        for (slong j = 0; j < ncols; j++) {
            field_init_elem(&test_vec[j], ctx->field_id, ctx_ptr);
            field_set_elem(&test_vec[j], &unified_mat[row * ncols + j], ctx->field_id, ctx_ptr);
        }
        
        // Eliminate using existing basis vectors
        for (slong i = 0; i < tracker.count; i++) {
            slong pivot_col = tracker.pivot_positions[i];
            if (!field_is_zero(&test_vec[pivot_col], ctx->field_id, ctx_ptr)) {
                field_elem_u factor;
                field_init_elem(&factor, ctx->field_id, ctx_ptr);
                field_set_elem(&factor, &test_vec[pivot_col], ctx->field_id, ctx_ptr);
                
                for (slong j = 0; j < ncols; j++) {
                    field_elem_u temp;
                    field_init_elem(&temp, ctx->field_id, ctx_ptr);
                    field_mul(&temp, &factor, &tracker.rows[i * ncols + j], ctx->field_id, ctx_ptr);
                    field_sub(&test_vec[j], &test_vec[j], &temp, ctx->field_id, ctx_ptr);
                    field_clear_elem(&temp, ctx->field_id, ctx_ptr);
                }
                
                field_clear_elem(&factor, ctx->field_id, ctx_ptr);
            }
        }
        
        // Find first non-zero position
        slong pivot_pos = -1;
        for (slong j = 0; j < ncols; j++) {
            if (!field_is_zero(&test_vec[j], ctx->field_id, ctx_ptr)) {
                pivot_pos = j;
                break;
            }
        }
        
        if (pivot_pos >= 0) {
            // Linearly independent, add to basis
            selected_rows[rank] = row;
            tracker.pivot_positions[tracker.count] = pivot_pos;
            
            // Normalize and store
            field_elem_u inv;
            field_init_elem(&inv, ctx->field_id, ctx_ptr);
            field_inv(&inv, &test_vec[pivot_pos], ctx->field_id, ctx_ptr);
            
            for (slong j = 0; j < ncols; j++) {
                field_mul(&tracker.rows[tracker.count * ncols + j], &test_vec[j], &inv, ctx->field_id, ctx_ptr);
            }
            
            field_clear_elem(&inv, ctx->field_id, ctx_ptr);
            tracker.count++;
            rank++;
        }
        
        // Clean up test vector
        for (slong j = 0; j < ncols; j++) {
            field_clear_elem(&test_vec[j], ctx->field_id, ctx_ptr);
        }
        flint_free(test_vec);
    }
    
    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    
    // Set output
    if (rank > 0) {
        *selected_rows_out = (slong*) flint_malloc(rank * sizeof(slong));
        memcpy(*selected_rows_out, selected_rows, rank * sizeof(slong));
    } else {
        *selected_rows_out = NULL;
    }
    *num_selected = rank;
    
    // Cleanup
    for (slong i = 0; i < max_rank * ncols; i++) {
        field_clear_elem(&tracker.rows[i], ctx->field_id, ctx_ptr);
    }
    flint_free(tracker.rows);
    flint_free(tracker.pivot_positions);
    flint_free(selected_rows);
    
    //printf("Non-prime field computation completed in %.4f seconds using unified interface\n", elapsed);
}

int g_matrix_transpose_threshold = 1000;

// Optimized find_fq_optimal_maximal_rank_submatrix
void find_fq_optimal_maximal_rank_submatrix(fq_mvpoly_t ***full_matrix, 
                                           slong nrows, slong ncols,
                                           slong **row_indices_out, 
                                           slong **col_indices_out,
                                           slong *num_rows, slong *num_cols,
                                           slong npars) {
    // Get context
    const fq_nmod_ctx_struct *ctx = NULL;
    for (slong i = 0; i < nrows && !ctx; i++) {
        for (slong j = 0; j < ncols && !ctx; j++) {
            if (full_matrix[i][j] != NULL) {
                ctx = full_matrix[i][j]->ctx;
                break;
            }
        }
    }
    
    const int use_extension_specialization = (fq_nmod_ctx_degree(ctx) == 1 && fq_nmod_ctx_prime(ctx) <= 100 && npars > 0);
    const slong MAX_SELECTION_ATTEMPTS = use_extension_specialization ? 3 : 2;

    field_ctx_t selection_ctx;
    fq_nmod_ctx_t extension_eval_ctx;
    if (use_extension_specialization) {
        fmpz_t selection_prime;
        fmpz_init_set_ui(selection_prime, fq_nmod_ctx_prime(ctx));
        fq_nmod_ctx_init(extension_eval_ctx, selection_prime, 2, "u");
        fmpz_clear(selection_prime);
        field_ctx_init(&selection_ctx, extension_eval_ctx);
    } else {
        field_ctx_init(&selection_ctx, ctx);
    }

    void *ctx_ptr = (selection_ctx.field_id == FIELD_ID_NMOD) ?
                   (void*)&selection_ctx.ctx.nmod_ctx :
                   (void*)selection_ctx.ctx.fq_ctx;

    slong accepted_size = 0;
    slong *accepted_rows = NULL;
    slong *accepted_cols = NULL;
    slong accepted_score = -1;
    int selection_stable = 0;

    for (slong selection_attempt = 0;
         selection_attempt < MAX_SELECTION_ATTEMPTS && !selection_stable;
         selection_attempt++) {
        fq_nmod_t *param_vals = NULL;
        if (use_extension_specialization) {
            param_vals = (fq_nmod_t*) flint_malloc(npars * sizeof(fq_nmod_t));
            init_extension_evaluation_parameters(param_vals, npars,
                                                 extension_eval_ctx,
                                                 selection_attempt);
        } else {
            param_vals = (fq_nmod_t*) flint_malloc(npars * sizeof(fq_nmod_t));
            init_evaluation_parameters(param_vals, npars, ctx, selection_attempt);
        }

        field_elem_u *unified_mat = (field_elem_u*) flint_malloc(nrows * ncols * sizeof(field_elem_u));
        slong *current_row_indices = NULL;
        slong *current_col_indices = NULL;
        slong *prev_row_indices = NULL;
        slong *prev_col_indices = NULL;
        slong current_size = 0;
        slong prev_size = 0;
        const slong MAX_ITERATIONS = 10;
        slong iteration = 0;
        int converged = 0;

        for (slong i = 0; i < nrows; i++) {
            for (slong j = 0; j < ncols; j++) {
                slong idx = i * ncols + j;
                field_init_elem(&unified_mat[idx], selection_ctx.field_id, ctx_ptr);

                if (use_extension_specialization) {
                    fq_nmod_t val;
                    fq_nmod_init(val, extension_eval_ctx);
                    if (full_matrix[i][j] == NULL) {
                        fq_nmod_zero(val, extension_eval_ctx);
                    } else {
                        evaluate_fq_mvpoly_at_extension_params(val, full_matrix[i][j],
                                                               param_vals, extension_eval_ctx);
                    }
                    fq_nmod_to_field_elem(&unified_mat[idx], val, &selection_ctx);
                    fq_nmod_clear(val, extension_eval_ctx);
                } else {
                    fq_nmod_t val;
                    fq_nmod_init(val, ctx);
                    if (full_matrix[i][j] == NULL) {
                        fq_nmod_zero(val, ctx);
                    } else {
                        evaluate_fq_mvpoly_at_params(val, full_matrix[i][j], param_vals);
                    }
                    fq_nmod_to_field_elem(&unified_mat[idx], val, &selection_ctx);
                    fq_nmod_clear(val, ctx);
                }
            }
        }

        clock_t iter_start = clock();

        while (iteration < MAX_ITERATIONS && !converged) {
            if (iteration > 0) {
                prev_row_indices = (slong*) flint_malloc(current_size * sizeof(slong));
                prev_col_indices = (slong*) flint_malloc(current_size * sizeof(slong));
                memcpy(prev_row_indices, current_row_indices, current_size * sizeof(slong));
                memcpy(prev_col_indices, current_col_indices, current_size * sizeof(slong));
                prev_size = current_size;
            }

            slong row_rank_selected = 0;

            if (iteration == 0) {
                slong *selected_rows = NULL;
                slong num_selected = 0;

                find_pivot_rows_simple(&selected_rows, &num_selected,
                                      unified_mat, nrows, ncols, &selection_ctx);

                current_size = num_selected;
                row_rank_selected = num_selected;
                current_row_indices = (slong*) flint_malloc(current_size * sizeof(slong));
                memcpy(current_row_indices, selected_rows, current_size * sizeof(slong));

                if (current_size > g_matrix_transpose_threshold) {
                    field_elem_u *transposed_mat = (field_elem_u*) flint_malloc(ncols * current_size * sizeof(field_elem_u));

                    for (slong i = 0; i < ncols * current_size; i++) {
                        field_init_elem(&transposed_mat[i], selection_ctx.field_id, ctx_ptr);
                    }

                    for (slong i = 0; i < current_size; i++) {
                        slong orig_row = current_row_indices[i];
                        for (slong j = 0; j < ncols; j++) {
                            slong src_idx = orig_row * ncols + j;
                            slong dst_idx = j * current_size + i;
                            field_set_elem(&transposed_mat[dst_idx], &unified_mat[src_idx],
                                          selection_ctx.field_id, ctx_ptr);
                        }
                    }

                    slong *selected_cols = NULL;
                    slong num_selected_cols = 0;

                    find_pivot_rows_simple(&selected_cols, &num_selected_cols,
                                          transposed_mat, ncols, current_size, &selection_ctx);

                    current_col_indices = (slong*) flint_malloc(num_selected_cols * sizeof(slong));
                    memcpy(current_col_indices, selected_cols, num_selected_cols * sizeof(slong));
                    current_size = FLINT_MIN(current_size, num_selected_cols);

                    for (slong i = 0; i < ncols * num_selected; i++) {
                        field_clear_elem(&transposed_mat[i], selection_ctx.field_id, ctx_ptr);
                    }
                    flint_free(transposed_mat);
                    flint_free(selected_rows);
                    flint_free(selected_cols);

                    converged = 1;
                    break;
                } else {
                    flint_free(selected_rows);
                }
            } else {
                fq_index_degree_pair *row_degrees = (fq_index_degree_pair*) flint_malloc(nrows * sizeof(fq_index_degree_pair));

                for (slong i = 0; i < nrows; i++) {
                    row_degrees[i].index = i;
                    row_degrees[i].degree = compute_fq_selected_cols_row_max_total_degree(
                        full_matrix, i, current_col_indices, current_size, npars);
                }

                qsort(row_degrees, nrows, sizeof(fq_index_degree_pair), compare_fq_degrees);

                field_elem_u *col_submat = (field_elem_u*) flint_malloc(nrows * current_size * sizeof(field_elem_u));
                slong col_submat_cols = current_size;
                for (slong i = 0; i < nrows * current_size; i++) {
                    field_init_elem(&col_submat[i], selection_ctx.field_id, ctx_ptr);
                }

                for (slong i = 0; i < nrows; i++) {
                    for (slong j = 0; j < current_size; j++) {
                        slong src_idx = i * ncols + current_col_indices[j];
                        slong dst_idx = i * current_size + j;
                        field_set_elem(&col_submat[dst_idx], &unified_mat[src_idx],
                                      selection_ctx.field_id, ctx_ptr);
                    }
                }

                unified_row_basis_tracker_t row_tracker;
                unified_row_basis_tracker_init(&row_tracker, current_size, current_size, &selection_ctx);

                for (slong i = 0; i < nrows && row_tracker.current_rank < current_size; i++) {
                    slong row_idx = row_degrees[i].index;
                    unified_try_add_row_to_basis(&row_tracker, col_submat, row_idx, current_size);
                }

                flint_free(current_row_indices);
                current_row_indices = (slong*) flint_malloc(row_tracker.current_rank * sizeof(slong));
                memcpy(current_row_indices, row_tracker.selected_indices, row_tracker.current_rank * sizeof(slong));
                row_rank_selected = row_tracker.current_rank;
                current_size = row_tracker.current_rank;

                for (slong i = 0; i < nrows * col_submat_cols; i++) {
                    field_clear_elem(&col_submat[i], selection_ctx.field_id, ctx_ptr);
                }
                flint_free(col_submat);
                unified_row_basis_tracker_clear(&row_tracker);
                flint_free(row_degrees);
            }

            fq_index_degree_pair *col_degrees = (fq_index_degree_pair*) flint_malloc(ncols * sizeof(fq_index_degree_pair));

            for (slong j = 0; j < ncols; j++) {
                col_degrees[j].index = j;
                col_degrees[j].degree = compute_fq_selected_rows_col_max_total_degree(
                    full_matrix, current_row_indices, current_size, j, npars);
            }

            qsort(col_degrees, ncols, sizeof(fq_index_degree_pair), compare_fq_degrees);

            field_elem_u *transposed = (field_elem_u*) flint_malloc(ncols * current_size * sizeof(field_elem_u));
            slong transposed_rows = current_size;
            for (slong i = 0; i < ncols * current_size; i++) {
                field_init_elem(&transposed[i], selection_ctx.field_id, ctx_ptr);
            }

            for (slong i = 0; i < current_size; i++) {
                for (slong j = 0; j < ncols; j++) {
                    slong src_idx = current_row_indices[i] * ncols + j;
                    slong dst_idx = j * current_size + i;
                    field_set_elem(&transposed[dst_idx], &unified_mat[src_idx],
                                  selection_ctx.field_id, ctx_ptr);
                }
            }

            unified_row_basis_tracker_t col_tracker;
            unified_row_basis_tracker_init(&col_tracker, ncols, current_size, &selection_ctx);

            for (slong j = 0; j < ncols && col_tracker.current_rank < current_size; j++) {
                slong col_idx = col_degrees[j].index;
                unified_try_add_row_to_basis(&col_tracker, transposed, col_idx, current_size);
            }

            slong col_rank_selected = col_tracker.current_rank;

            if (iteration == 0) {
                current_col_indices = (slong*) flint_malloc(col_tracker.current_rank * sizeof(slong));
            } else {
                flint_free(current_col_indices);
                current_col_indices = (slong*) flint_malloc(col_tracker.current_rank * sizeof(slong));
            }
            memcpy(current_col_indices, col_tracker.selected_indices, col_tracker.current_rank * sizeof(slong));
            current_size = FLINT_MIN(current_size, col_tracker.current_rank);

            for (slong i = 0; i < ncols * transposed_rows; i++) {
                field_clear_elem(&transposed[i], selection_ctx.field_id, ctx_ptr);
            }
            flint_free(transposed);
            unified_row_basis_tracker_clear(&col_tracker);
            flint_free(col_degrees);

            if (iteration > 0) {
                if (current_size == prev_size &&
                    current_size == FLINT_MIN(col_rank_selected, row_rank_selected) &&
                    indices_equal(current_row_indices, prev_row_indices, current_size) &&
                    indices_equal(current_col_indices, prev_col_indices, current_size)) {
                    converged = 1;
                }

                flint_free(prev_row_indices);
                flint_free(prev_col_indices);
                prev_row_indices = NULL;
                prev_col_indices = NULL;
            }

            iteration++;
        }

        clock_t iter_end = clock();
        (void) iter_start;
        (void) iter_end;

        slong final_rank;
        slong verify_rank;
        if (use_extension_specialization) {
            final_rank = evaluate_selected_submatrix_rank_extension(
                full_matrix, current_row_indices, current_col_indices, current_size,
                param_vals, extension_eval_ctx);

            fq_nmod_t *verify_vals = (fq_nmod_t*) flint_malloc(npars * sizeof(fq_nmod_t));
            init_extension_evaluation_parameters(verify_vals, npars,
                                                 extension_eval_ctx,
                                                 selection_attempt + 1);
            verify_rank = evaluate_selected_submatrix_rank_extension(
                full_matrix, current_row_indices, current_col_indices, current_size,
                verify_vals, extension_eval_ctx);
            clear_evaluation_parameters(verify_vals, npars, extension_eval_ctx);
        } else {
            final_rank = evaluate_selected_submatrix_rank(
                full_matrix, current_row_indices, current_col_indices, current_size,
                param_vals, ctx);

            fq_nmod_t *verify_vals = (fq_nmod_t*) flint_malloc(npars * sizeof(fq_nmod_t));
            init_evaluation_parameters(verify_vals, npars, ctx, selection_attempt + 1);
            verify_rank = evaluate_selected_submatrix_rank(
                full_matrix, current_row_indices, current_col_indices, current_size,
                verify_vals, ctx);
            clear_evaluation_parameters(verify_vals, npars, ctx);
        }

        slong candidate_score = 0;
        if (final_rank == current_size) candidate_score++;
        if (verify_rank == current_size) candidate_score++;

        if (current_size > accepted_size ||
            (current_size == accepted_size && candidate_score > accepted_score)) {
            if (accepted_rows) flint_free(accepted_rows);
            if (accepted_cols) flint_free(accepted_cols);
            accepted_rows = (slong*) flint_malloc(current_size * sizeof(slong));
            accepted_cols = (slong*) flint_malloc(current_size * sizeof(slong));
            memcpy(accepted_rows, current_row_indices, current_size * sizeof(slong));
            memcpy(accepted_cols, current_col_indices, current_size * sizeof(slong));
            accepted_size = current_size;
            accepted_score = candidate_score;
        }

        if (candidate_score == 2) {
            selection_stable = 1;
        } else {
            printf("Submatrix verification failed; retrying with a new specialization.\n");
        }

        if (prev_row_indices) flint_free(prev_row_indices);
        if (prev_col_indices) flint_free(prev_col_indices);
        flint_free(current_row_indices);
        flint_free(current_col_indices);
        if (param_vals) {
            clear_evaluation_parameters(param_vals, npars,
                                       use_extension_specialization ? extension_eval_ctx : ctx);
        }
        for (slong i = 0; i < nrows * ncols; i++) {
            field_clear_elem(&unified_mat[i], selection_ctx.field_id, ctx_ptr);
        }
        flint_free(unified_mat);
    }

    if (!accepted_rows || !accepted_cols) {
        *row_indices_out = NULL;
        *col_indices_out = NULL;
        *num_rows = 0;
        *num_cols = 0;
    } else {
        if (!selection_stable) {
            printf("Warning: using best available submatrix after %ld attempts.\n",
                   MAX_SELECTION_ATTEMPTS);
        }
        *row_indices_out = accepted_rows;
        *col_indices_out = accepted_cols;
        *num_rows = accepted_size;
        *num_cols = accepted_size;
    }

    field_ctx_clear(&selection_ctx);
    if (use_extension_specialization) {
        fq_nmod_ctx_clear(extension_eval_ctx);
    }
}

// Optimized monomial collection function - replaces the original O(n²) loop
void collect_unique_monomials(
    monom_t **x_monoms_out, slong *nx_monoms_out,
    monom_t **dual_monoms_out, slong *ndual_monoms_out,
    const fq_mvpoly_t *dixon_poly, 
    const slong *d0, const slong *d1, slong nvars) {
    
    if (dixon_poly->nterms == 0) {
        *x_monoms_out = NULL; *nx_monoms_out = 0;
        *dual_monoms_out = NULL; *ndual_monoms_out = 0;
        return;
    }
    
    /* Simple hash table size - power of 2 for fast modulo */
    slong hash_size = 1024;
    while (hash_size < dixon_poly->nterms) hash_size <<= 1;
    
    /* Hash tables for tracking unique monomials */
    hash_entry_t **x_buckets = (hash_entry_t**) flint_calloc(hash_size, sizeof(hash_entry_t*));
    hash_entry_t **dual_buckets = (hash_entry_t**) flint_calloc(hash_size, sizeof(hash_entry_t*));
    
    /* * COALESCED ALLOCATION STRATEGY:
     * We allocate a single block for both the monom_t array and the exponent data.
     * This ensures that when the caller frees the returned pointer, all associated 
     * memory is released without needing extra arguments.
     */
    slong max_terms = dixon_poly->nterms;
    
    /* Calculate sizes for X-monomials */
    slong x_structs_size = max_terms * sizeof(monom_t);
    slong x_data_size = max_terms * nvars * sizeof(slong);
    char *x_combined = (char *) flint_malloc(x_structs_size + x_data_size);
    
    monom_t *x_monoms = (monom_t *) x_combined;
    slong *x_exp_storage = (slong *) (x_combined + x_structs_size);
    
    /* Calculate sizes for Dual-monomials */
    slong dual_structs_size = max_terms * sizeof(monom_t);
    slong dual_data_size = max_terms * nvars * sizeof(slong);
    char *dual_combined = (char *) flint_malloc(dual_structs_size + dual_data_size);
    
    monom_t *dual_monoms = (monom_t *) dual_combined;
    slong *dual_exp_storage = (slong *) (dual_combined + dual_structs_size);

    slong nx_monoms = 0, ndual_monoms = 0;
    
    /* Process each term in the Dixon polynomial */
    for (slong i = 0; i < dixon_poly->nterms; i++) {
        const slong *var_exp = dixon_poly->terms[i].var_exp;
        if (!var_exp) continue;
        
        /* Check degree bounds d0 and d1 */
        int valid = 1;
        for (slong k = 0; k < nvars && valid; k++) {
            if (var_exp[k] >= d0[k] || var_exp[nvars + k] >= d1[k]) {
                valid = 0;
            }
        }
        if (!valid) continue;
        
        /* 1. Process x-monomial (first nvars components) */
        ulong x_hash = 0;
        for (slong k = 0; k < nvars; k++) {
            x_hash = x_hash * 31 + var_exp[k];
        }
        x_hash &= (hash_size - 1);
        
        hash_entry_t *entry = x_buckets[x_hash];
        int found = 0;
        while (entry && !found) {
            if (memcmp(entry->exp, var_exp, nvars * sizeof(slong)) == 0) {
                found = 1;
            } else {
                entry = entry->next;
            }
        }
        
        if (!found) {
            /* Map the .exp pointer to the correct offset in the storage block */
            slong *x_exp = &x_exp_storage[nx_monoms * nvars];
            memcpy(x_exp, var_exp, nvars * sizeof(slong));
            
            x_monoms[nx_monoms].exp = x_exp;
            x_monoms[nx_monoms].idx = nx_monoms;
            
            /* Add to hash table for future lookup */
            hash_entry_t *new_entry = (hash_entry_t*) flint_malloc(sizeof(hash_entry_t));
            new_entry->exp = x_exp;
            new_entry->idx = nx_monoms;
            new_entry->next = x_buckets[x_hash];
            x_buckets[x_hash] = new_entry;
            
            nx_monoms++;
        }
        
        /* 2. Process dual monomial (next nvars components) */
        const slong *dual_exp_src = &var_exp[nvars];
        ulong dual_hash = 0;
        for (slong k = 0; k < nvars; k++) {
            dual_hash = dual_hash * 31 + dual_exp_src[k];
        }
        dual_hash &= (hash_size - 1);
        
        entry = dual_buckets[dual_hash];
        found = 0;
        while (entry && !found) {
            if (memcmp(entry->exp, dual_exp_src, nvars * sizeof(slong)) == 0) {
                found = 1;
            } else {
                entry = entry->next;
            }
        }
        
        if (!found) {
            slong *dual_exp = &dual_exp_storage[ndual_monoms * nvars];
            memcpy(dual_exp, dual_exp_src, nvars * sizeof(slong));
            
            dual_monoms[ndual_monoms].exp = dual_exp;
            dual_monoms[ndual_monoms].idx = ndual_monoms;
            
            hash_entry_t *new_entry = (hash_entry_t*) flint_malloc(sizeof(hash_entry_t));
            new_entry->exp = dual_exp;
            new_entry->idx = ndual_monoms;
            new_entry->next = dual_buckets[dual_hash];
            dual_buckets[dual_hash] = new_entry;
            
            ndual_monoms++;
        }
    }
    
    /* Clean up hash table metadata (actual monom data is in combined blocks) */
    for (slong i = 0; i < hash_size; i++) {
        hash_entry_t *curr;
        curr = x_buckets[i];
        while (curr) {
            hash_entry_t *next = curr->next;
            flint_free(curr);
            curr = next;
        }
        curr = dual_buckets[i];
        while (curr) {
            hash_entry_t *next = curr->next;
            flint_free(curr);
            curr = next;
        }
    }
    flint_free(x_buckets);
    flint_free(dual_buckets);
    
    /* Handle empty results */
    if (nx_monoms == 0) {
        flint_free(x_combined);
        x_monoms = NULL;
    }
    if (ndual_monoms == 0) {
        flint_free(dual_combined);
        dual_monoms = NULL;
    }
    
    /* Final output assignment */
    *x_monoms_out = x_monoms;
    *nx_monoms_out = nx_monoms;
    *dual_monoms_out = dual_monoms;
    *ndual_monoms_out = ndual_monoms;
    
}

// Allocate single element on demand
fq_mvpoly_t* get_matrix_entry_lazy(fq_mvpoly_t ***matrix, slong i, slong j,
                                  slong npars, const fq_nmod_ctx_t ctx) {
    if (!matrix[i][j]) {
        matrix[i][j] = (fq_mvpoly_t*) flint_malloc(sizeof(fq_mvpoly_t));
        fq_mvpoly_init(matrix[i][j], 0, npars, ctx);
    }
    return matrix[i][j];
}
void fill_coefficient_matrix_optimized(fq_mvpoly_t ***full_matrix,
                                      monom_t *x_monoms, slong nx_monoms,
                                      monom_t *dual_monoms, slong ndual_monoms,
                                      const fq_mvpoly_t *dixon_poly,
                                      const slong *d0, const slong *d1, 
                                      slong nvars, slong npars) {
    
    // Find corresponding matrix positions for each Dixon polynomial term
    for (slong t = 0; t < dixon_poly->nterms; t++) {
        if (!dixon_poly->terms[t].var_exp) continue;
        
        // Check degree bounds
        int valid = 1;
        for (slong k = 0; k < nvars; k++) {
            if (dixon_poly->terms[t].var_exp[k] >= d0[k] || 
                dixon_poly->terms[t].var_exp[nvars + k] >= d1[k]) {
                valid = 0;
                break;
            }
        }
        if (!valid) continue;
        
        // Find x-monomial row index
        slong row = -1;
        for (slong i = 0; i < nx_monoms; i++) {
            if (memcmp(x_monoms[i].exp, dixon_poly->terms[t].var_exp, 
                      nvars * sizeof(slong)) == 0) {
                row = i;
                break;
            }
        }
        
        // Find dual-monomial column index
        slong col = -1;
        for (slong j = 0; j < ndual_monoms; j++) {
            if (memcmp(dual_monoms[j].exp, &dixon_poly->terms[t].var_exp[nvars], 
                      nvars * sizeof(slong)) == 0) {
                col = j;
                break;
            }
        }
        
        // Only allocate memory for positions that are actually needed
        if (row >= 0 && col >= 0) {
            fq_mvpoly_t *entry = get_matrix_entry_lazy(full_matrix, row, col, 
                                                      npars, dixon_poly->ctx);
            fq_mvpoly_add_term_fast(entry, NULL, dixon_poly->terms[t].par_exp, 
                                   dixon_poly->terms[t].coeff);
        }
    }
}
// Optimized version of find_fq_optimal_maximal_rank_submatrix
// ============ Extract coefficient matrix ============

void extract_fq_coefficient_matrix_from_dixon(fq_mvpoly_t ***coeff_matrix,
                                              slong *row_indices, slong *col_indices,
                                              slong *matrix_size,
                                              const fq_mvpoly_t *dixon_poly,
                                              slong nvars, slong npars) {
    printf("\nStep 2: Construct Dixon matrix\n");
    
    slong *d0 = (slong*) flint_calloc(nvars, sizeof(slong));
    slong *d1 = (slong*) flint_calloc(nvars, sizeof(slong));
    
    for (slong i = 0; i < dixon_poly->nterms; i++) {
        if (dixon_poly->terms[i].var_exp) {
            for (slong j = 0; j < nvars; j++) {
                if (dixon_poly->terms[i].var_exp[j] > d0[j]) {
                    d0[j] = dixon_poly->terms[i].var_exp[j];
                }
            }
            for (slong j = 0; j < nvars; j++) {
                if (dixon_poly->terms[i].var_exp[nvars + j] > d1[j]) {
                    d1[j] = dixon_poly->terms[i].var_exp[nvars + j];
                }
            }
        }
    }
    
    for (slong i = 0; i < nvars; i++) {
        d0[i]++;
        d1[i]++;
    }
    
    slong expected_rows = 1;
    slong expected_cols = 1;
    for (slong i = 0; i < nvars; i++) {
        expected_rows *= d0[i];
        expected_cols *= d1[i];
    }

    monom_t *x_monoms = NULL;
    monom_t *dual_monoms = NULL;
    slong nx_monoms = 0, ndual_monoms = 0;
    
    collect_unique_monomials(&x_monoms, &nx_monoms,
                        &dual_monoms, &ndual_monoms,
                        dixon_poly, d0, d1, nvars);

    printf("Dixon matrix size: %ld x %ld\n", nx_monoms, ndual_monoms);
    
    if (nx_monoms == 0 || ndual_monoms == 0) {
        printf("Warning: Empty coefficient matrix\n");
        *matrix_size = 0;
        flint_free(d0);
        flint_free(d1);
        if (x_monoms) flint_free(x_monoms);
        if (dual_monoms) flint_free(dual_monoms);
        return;
    }

    fq_mvpoly_t ***full_matrix = (fq_mvpoly_t***) flint_malloc(nx_monoms * sizeof(fq_mvpoly_t**));
    for (slong i = 0; i < nx_monoms; i++) {
        full_matrix[i] = (fq_mvpoly_t**) flint_calloc(ndual_monoms, sizeof(fq_mvpoly_t*));
    }

    fill_coefficient_matrix_optimized(full_matrix, x_monoms, nx_monoms, 
                                     dual_monoms, ndual_monoms, dixon_poly, 
                                     d0, d1, nvars, npars);

    slong *row_idx_array = NULL;
    slong *col_idx_array = NULL;
    slong num_rows, num_cols;
    
    if (npars == 0) {
        fq_nmod_mat_t eval_mat;
        fq_nmod_mat_init(eval_mat, nx_monoms, ndual_monoms, dixon_poly->ctx);
        for (slong i = 0; i < nx_monoms; i++) {
            for (slong j = 0; j < ndual_monoms; j++) {
                if (full_matrix[i][j] != NULL && full_matrix[i][j]->nterms > 0) {
                    fq_nmod_set(fq_nmod_mat_entry(eval_mat, i, j), 
                               full_matrix[i][j]->terms[0].coeff, dixon_poly->ctx);
                } else {
                    fq_nmod_zero(fq_nmod_mat_entry(eval_mat, i, j), dixon_poly->ctx);
                }
            }
        }
        slong rank = fq_nmod_mat_rank(eval_mat, dixon_poly->ctx);
        slong min_size = FLINT_MIN(nx_monoms, ndual_monoms);
        slong actual_size = FLINT_MIN(rank, min_size);
        
        row_idx_array = (slong*) flint_malloc(actual_size * sizeof(slong));
        col_idx_array = (slong*) flint_malloc(actual_size * sizeof(slong));
        
        for (slong i = 0; i < actual_size; i++) {
            row_idx_array[i] = i;
            col_idx_array[i] = i;
        }
        num_rows = actual_size;
        num_cols = actual_size;
        
        fq_nmod_mat_clear(eval_mat, dixon_poly->ctx);
    } else {
        slong small_size = 1;
        if (nx_monoms < small_size && ndual_monoms < small_size && 
            expected_rows < small_size && expected_cols < small_size) {
            slong min_size = FLINT_MIN(nx_monoms, ndual_monoms);
            row_idx_array = (slong*) flint_malloc(min_size * sizeof(slong));
            col_idx_array = (slong*) flint_malloc(min_size * sizeof(slong));
            
            for (slong i = 0; i < min_size; i++) {
                row_idx_array[i] = i;
                col_idx_array[i] = i;
            }
            num_rows = min_size;
            num_cols = min_size;
        } else {
            find_fq_optimal_maximal_rank_submatrix(full_matrix, nx_monoms, ndual_monoms,
                                                  &row_idx_array, &col_idx_array, 
                                                  &num_rows, &num_cols,
                                                  npars);
        }
    }

    slong submat_rank = FLINT_MIN(num_rows, num_cols);
    printf("\nStep 3: Extract maximal-rank submatrix\n");
    
    if (submat_rank == 0) {
        printf("Warning: Matrix has rank 0\n");
        *matrix_size = 0;
        
        if (full_matrix) {
            for (slong i = 0; i < nx_monoms; i++) {
                if (full_matrix[i]) {
                    for (slong j = 0; j < ndual_monoms; j++) {
                        if (full_matrix[i][j] != NULL) {
                            fq_mvpoly_clear(full_matrix[i][j]);
                            flint_free(full_matrix[i][j]);
                        }
                    }
                    flint_free(full_matrix[i]);
                }
            }
            flint_free(full_matrix);
        }
        
        if (row_idx_array) flint_free(row_idx_array);
        if (col_idx_array) flint_free(col_idx_array);
        if (x_monoms) flint_free(x_monoms);
        if (dual_monoms) flint_free(dual_monoms);
        if (d0) flint_free(d0);
        if (d1) flint_free(d1);
        return;
    }

    printf("Submatrix size: %ld x %ld\n", submat_rank, submat_rank);
    
    *coeff_matrix = (fq_mvpoly_t**) flint_malloc(submat_rank * sizeof(fq_mvpoly_t*));
    for (slong i = 0; i < submat_rank; i++) {
        (*coeff_matrix)[i] = (fq_mvpoly_t*) flint_malloc(submat_rank * sizeof(fq_mvpoly_t));
        for (slong j = 0; j < submat_rank; j++) {
            fq_mvpoly_t *source = full_matrix[row_idx_array[i]][col_idx_array[j]];
            if (source != NULL) {
                fq_mvpoly_copy(&(*coeff_matrix)[i][j], source);
            } else {
                fq_mvpoly_init(&(*coeff_matrix)[i][j], 0, npars, dixon_poly->ctx);
            }
        }
    }
    
    for (slong i = 0; i < submat_rank; i++) {
        row_indices[i] = row_idx_array[i];
        col_indices[i] = col_idx_array[i];
    }
    *matrix_size = submat_rank;
    
    if (full_matrix) {
        for (slong i = 0; i < nx_monoms; i++) {
            if (full_matrix[i]) {
                for (slong j = 0; j < ndual_monoms; j++) {
                    if (full_matrix[i][j]) {
                        fq_mvpoly_clear(full_matrix[i][j]);
                        flint_free(full_matrix[i][j]);
                    }
                }
                flint_free(full_matrix[i]);
            }
        }
        flint_free(full_matrix);
    }
    
    if (row_idx_array) flint_free(row_idx_array);
    if (col_idx_array) flint_free(col_idx_array);
    if (x_monoms) flint_free(x_monoms);
    if (dual_monoms) flint_free(dual_monoms);
    if (d0) flint_free(d0);
    if (d1) flint_free(d1);
}

// Compute determinant of cancellation matrix
void compute_fq_cancel_matrix_det(fq_mvpoly_t *result, fq_mvpoly_t **modified_M_mvpoly,
                                  slong nvars, slong npars, det_method_t method) {
    clock_t start = clock();
    switch (method) {
        case DET_METHOD_INTERPOLATION:
            fq_compute_det_by_interpolation_optimized(result, modified_M_mvpoly,
                                                      nvars + 1, nvars, npars,
                                                      modified_M_mvpoly[0][0].ctx, NULL);
            break;
        case DET_METHOD_KRONECKER:
            compute_fq_det_kronecker(result, modified_M_mvpoly, nvars + 1);
            break;
        case DET_METHOD_HUANG:
            compute_fq_det_huang_interpolation(result, modified_M_mvpoly, nvars + 1);
            break;
        case DET_METHOD_RECURSIVE:
        default:
            compute_fq_det_recursive(result, modified_M_mvpoly, nvars + 1);
            break;
    }
    clock_t end = clock();
    (void) start;
    (void) end;
}

// ============ Get size of Dixon matrix ============

static slong dixon_binomial(slong n, slong k) {
    if (k > n || k < 0) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;
    
    slong result = 1;
    for (slong i = 0; i < k; i++) {
        result = result * (n - i) / (i + 1);
    }
    return result;
}

static void dixon_add_generic_monomials(fq_mvpoly_t *poly,
                                        slong nvars,
                                        slong *exp,
                                        slong pos,
                                        slong remaining,
                                        flint_rand_t state,
                                        const fq_nmod_ctx_t ctx) {
    if (pos == nvars) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, ctx);
        do {
            fq_nmod_randtest(coeff, state, ctx);
        } while (fq_nmod_is_zero(coeff, ctx));

        slong *var_exp = (slong*) malloc(nvars * sizeof(slong));
        memcpy(var_exp, exp, nvars * sizeof(slong));
        fq_mvpoly_add_term(poly, var_exp, NULL, coeff);

        fq_nmod_clear(coeff, ctx);
        free(var_exp);
        return;
    }

    for (slong d = 0; d <= remaining; d++) {
        exp[pos] = d;
        dixon_add_generic_monomials(poly, nvars, exp, pos + 1, remaining - d, state, ctx);
    }
}

// Calculate actual Dixon matrix size by building the system and extracting matrix
slong dixon_matrix_size(slong nvars, slong degree, ulong prime, slong field_degree) {
    printf("\nCalculating Dixon matrix size for n=%ld, d=%ld\n", nvars, degree);
    
    // Initialize field context
    fq_nmod_ctx_t ctx;
    fmpz_t p;
    fmpz_init_set_ui(p, prime);
    fq_nmod_ctx_init(ctx, p, field_degree, "t");
    fmpz_clear(p);
    
    // Generate generic polynomial system
    nvars--;
    slong npolys = nvars + 1;
    fq_mvpoly_t *polys = (fq_mvpoly_t*) malloc(npolys * sizeof(fq_mvpoly_t));
    
    // Generate dense generic polynomials of degree d
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, 12345, 67890);
    
    for (slong i = 0; i < npolys; i++) {
        fq_mvpoly_init(&polys[i], nvars, 0, ctx);  // No parameters
                
        // Generate all monomials of degree <= d
        slong total_monomials = dixon_binomial(nvars + degree, degree);
        
        // Enumerate all monomials recursively
        slong *temp_exp = (slong*) calloc(nvars, sizeof(slong));
        for (slong d = 0; d <= degree; d++) {
            dixon_add_generic_monomials(&polys[i], nvars, temp_exp, 0, d, state, ctx);
        }
        free(temp_exp);
    }
    
    // Build cancellation matrix
    fq_mvpoly_t **M_mvpoly;
    build_fq_cancellation_matrix_mvpoly(&M_mvpoly, polys, nvars, 0);
    
    // Perform row operations
    fq_mvpoly_t **modified_M_mvpoly;
    perform_fq_matrix_row_operations_mvpoly(&modified_M_mvpoly, &M_mvpoly, nvars, 0);
    
    // Compute determinant to get Dixon polynomial
    fq_mvpoly_t d_poly;
    compute_fq_cancel_matrix_det(&d_poly, modified_M_mvpoly, nvars, 0, DET_METHOD_RECURSIVE);
    
    // Calculate degree bounds
    slong *d0 = (slong*) calloc(nvars, sizeof(slong));
    slong *d1 = (slong*) calloc(nvars, sizeof(slong));
    
    for (slong i = 0; i < d_poly.nterms; i++) {
        if (d_poly.terms[i].var_exp) {
            for (slong j = 0; j < nvars; j++) {
                if (d_poly.terms[i].var_exp[j] > d0[j]) {
                    d0[j] = d_poly.terms[i].var_exp[j];
                }
            }
            for (slong j = 0; j < nvars; j++) {
                if (d_poly.terms[i].var_exp[nvars + j] > d1[j]) {
                    d1[j] = d_poly.terms[i].var_exp[nvars + j];
                }
            }
        }
    }
    
    for (slong i = 0; i < nvars; i++) {
        d0[i]++;
        d1[i]++;
    }
    
    // Collect unique monomials
    monom_t *x_monoms = NULL;
    monom_t *dual_monoms = NULL;
    slong nx_monoms = 0, ndual_monoms = 0;
    
    collect_unique_monomials(&x_monoms, &nx_monoms,
                            &dual_monoms, &ndual_monoms,
                            &d_poly, d0, d1, nvars);
    
    // The actual matrix size (should be square, so take minimum)
    slong matrix_size = FLINT_MAX(nx_monoms, ndual_monoms);
    
    printf("Dixon matrix actual size: %ld x %ld\n", nx_monoms, ndual_monoms);
    
    // Cleanup
    if (x_monoms) flint_free(x_monoms);
    if (dual_monoms) flint_free(dual_monoms);
    flint_free(d0);
    flint_free(d1);
    
    fq_mvpoly_clear(&d_poly);
    
    for (slong i = 0; i <= nvars; i++) {
        for (slong j = 0; j <= nvars; j++) {
            fq_mvpoly_clear(&M_mvpoly[i][j]);
            fq_mvpoly_clear(&modified_M_mvpoly[i][j]);
        }
        flint_free(M_mvpoly[i]);
        flint_free(modified_M_mvpoly[i]);
    }
    flint_free(M_mvpoly);
    flint_free(modified_M_mvpoly);
    
    for (slong i = 0; i < npolys; i++) {
        fq_mvpoly_clear(&polys[i]);
    }
    free(polys);
    
    flint_rand_clear(state);
    fq_nmod_ctx_clear(ctx);
    
    return matrix_size;
}

// ============ Main Dixon resultant function ============
void fq_dixon_resultant(fq_mvpoly_t *result, fq_mvpoly_t *polys, 
                       slong nvars, slong npars) {
    cleanup_unified_workspace();
    printf("\nStep 1: Build Dixon polynomial\n");
    clock_t step1_start = clock();
    fq_mvpoly_t **M_mvpoly;
    printf("Build Cancellation Matrix\n");
    build_fq_cancellation_matrix_mvpoly(&M_mvpoly, polys, nvars, npars);
    
    // Display analysis of original matrix
    //analyze_fq_matrix_mvpoly(M_mvpoly, nvars + 1, nvars + 1, "Original Cancellation");
    
    fq_mvpoly_t **modified_M_mvpoly;
    printf("Perform Matrix Row Operations\n");
    perform_fq_matrix_row_operations_mvpoly(&modified_M_mvpoly, &M_mvpoly, nvars, npars);
    
    fq_mvpoly_t d_poly;
    det_method_t step1_method = DET_METHOD_RECURSIVE;
    if (dixon_global_method_step1 != -1) {
        step1_method = dixon_global_method_step1;
        printf("Step 1 method override active: %d (%s)\n",
               dixon_global_method_step1, dixon_det_method_name(dixon_global_method_step1));
    }
    printf("Computing cancellation matrix determinant using %s...\n",
           dixon_det_method_name(step1_method));
    compute_fq_cancel_matrix_det(&d_poly, modified_M_mvpoly, nvars, npars, step1_method);
    
    if (d_poly.nterms <= 100) {
        printf("Dixon polynomial: %ld terms\n", d_poly.nterms);
        fq_mvpoly_print_expanded(&d_poly, "DixonPoly", 1);
    } else {
        printf("Dixon polynomial: %ld terms (not shown)\n", d_poly.nterms);
    }
    printf("Time: %.3f seconds\n", (double)(clock() - step1_start) / CLOCKS_PER_SEC);

    for (slong i = 0; i <= nvars; i++) {
        for (slong j = 0; j <= nvars; j++) {
            fq_mvpoly_clear(&M_mvpoly[i][j]);
            fq_mvpoly_clear(&modified_M_mvpoly[i][j]);
        }
        flint_free(M_mvpoly[i]);
        flint_free(modified_M_mvpoly[i]);
    }
    flint_free(M_mvpoly);
    flint_free(modified_M_mvpoly);
    
    fq_mvpoly_t **coeff_matrix;
    slong *row_indices = (slong*) flint_malloc(d_poly.nterms * sizeof(slong));
    slong *col_indices = (slong*) flint_malloc(d_poly.nterms * sizeof(slong));
    slong matrix_size;
    
    extract_fq_coefficient_matrix_from_dixon(&coeff_matrix, row_indices, col_indices,
                                            &matrix_size, &d_poly, nvars, npars);
    
    if (matrix_size > 0) {
        printf("\nStep 4: Compute resultant\n");
        
        slong res_deg_bound = compute_fq_dixon_resultant_degree_bound(polys, nvars+1, nvars, npars);
        printf("Degree bound: %ld\n", res_deg_bound);
        
        ulong field_size = 1;
        for (slong i = 0; i < fq_nmod_ctx_degree(polys[0].ctx); i++) {
            field_size *= fq_nmod_ctx_prime(polys[0].ctx);
        }
        
        det_method_t coeff_method; // DET_METHOD_RECURSIVE DET_METHOD_KRONECKER DET_METHOD_INTERPOLATION DET_METHOD_HUANG
        #ifdef _OPENMP
        if (npars > 1) {
            coeff_method = DET_METHOD_INTERPOLATION;
        } else 
        #endif
        if (matrix_size < 9) {
            coeff_method = DET_METHOD_RECURSIVE;
        } else {
            coeff_method = DET_METHOD_KRONECKER;
        }
        //coeff_method = DET_METHOD_INTERPOLATION;
        if (dixon_global_method_step4 != -1) {
            coeff_method = dixon_global_method_step4;
            printf("Step 4 method override active: %d (%s)\n",
                   dixon_global_method_step4, dixon_det_method_name(dixon_global_method_step4));
        }
        
        compute_fq_coefficient_matrix_det(result, coeff_matrix, matrix_size,
                                         npars, polys[0].ctx, coeff_method, res_deg_bound);
        
        if (result->nterms <= 100) {
            fq_mvpoly_print(result, "Final Resultant");
        } else {
            printf("Final resultant too large to display (%ld terms)\n", result->nterms);
        }
        fq_mvpoly_make_monic(result);
        // Cleanup coefficient matrix
        for (slong i = 0; i < matrix_size; i++) {
            for (slong j = 0; j < matrix_size; j++) {
                fq_mvpoly_clear(&coeff_matrix[i][j]);
            }
            flint_free(coeff_matrix[i]);
        }
        flint_free(coeff_matrix);
    } else {
        fq_mvpoly_init(result, 0, npars, polys[0].ctx);
        printf("Warning: Empty coefficient matrix, resultant is 0\n");
    }
    
    // Cleanup
    flint_free(row_indices);
    flint_free(col_indices);
    
    fq_mvpoly_clear(&d_poly);

    printf("\n=== Dixon Resultant Computation Complete ===\n");
}

void fq_dixon_resultant_with_names(fq_mvpoly_t *result, fq_mvpoly_t *polys, 
                                  slong nvars, slong npars,
                                  char **var_names, char **par_names, 
                                  const char *gen_name) {
    cleanup_unified_workspace();
    
    printf("\nStep 1: Build Dixon polynomial\n");
    clock_t step1_cpu_start = clock();
    double step1_wall_start = get_wall_time();
    fq_mvpoly_t **M_mvpoly;
    printf("Build Cancellation Matrix\n");
    build_fq_cancellation_matrix_mvpoly(&M_mvpoly, polys, nvars, npars);
    
    fq_mvpoly_t **modified_M_mvpoly;
    printf("Perform Matrix Row Operations\n");
    perform_fq_matrix_row_operations_mvpoly(&modified_M_mvpoly, &M_mvpoly, nvars, npars);

    fq_mvpoly_t d_poly;
    det_method_t step1_method = DET_METHOD_RECURSIVE;
    if (dixon_global_method_step1 != -1) {
        step1_method = dixon_global_method_step1;
        printf("Step 1 method override active: %d (%s)\n",
               dixon_global_method_step1, dixon_det_method_name(dixon_global_method_step1));
    }
    printf("Computing cancellation matrix determinant using %s...\n",
           dixon_det_method_name(step1_method));
    compute_fq_cancel_matrix_det(&d_poly, modified_M_mvpoly, nvars, npars, step1_method);
    
    if (d_poly.nterms <= 100) {
        printf("Dixon polynomial: %ld terms\n", d_poly.nterms);
        fq_mvpoly_print_with_names(&d_poly, "DixonPoly", var_names, par_names, gen_name, 1);
    } else {
        printf("Dixon polynomial: %ld terms (not shown)\n", d_poly.nterms);
    }
    clock_t step1_cpu_end = clock();
    double step1_wall_end = get_wall_time();
    double step1_cpu_elapsed = (double)(step1_cpu_end - step1_cpu_start) / CLOCKS_PER_SEC;
    double step1_wall_elapsed = step1_wall_end - step1_wall_start;
    int threads = 1;
    #ifdef _OPENMP
    threads = omp_get_max_threads();
    #endif
    // printf("CPU time: %.3f seconds | Wall time: %.3f seconds | Threads: %d\n", step1_cpu_elapsed, step1_wall_elapsed, threads);
    
    fq_mvpoly_t **coeff_matrix;
    slong max_indices = d_poly.nterms > 0 ? d_poly.nterms : 1;
    slong *row_indices = (slong*) flint_malloc(max_indices * sizeof(slong));
    slong *col_indices = (slong*) flint_malloc(max_indices * sizeof(slong));
    slong matrix_size;
    
    extract_fq_coefficient_matrix_from_dixon(&coeff_matrix, row_indices, col_indices,
                                            &matrix_size, &d_poly, nvars, npars);
    
    if (matrix_size > 0) {
        printf("\nStep 4: Compute resultant\n");
        
        slong res_deg_bound = compute_fq_dixon_resultant_degree_bound(polys, nvars+1, nvars, npars);
        printf("Degree bound: %ld\n", res_deg_bound);
        
        det_method_t coeff_method;
        #ifdef _OPENMP
        if (npars > 1) {
            coeff_method = DET_METHOD_INTERPOLATION;
        } else 
        #endif
        if (matrix_size < 9) {
            coeff_method = DET_METHOD_RECURSIVE;
        } else {
            coeff_method = DET_METHOD_KRONECKER;
        }
        if (dixon_global_method_step4 != -1) {
            coeff_method = dixon_global_method_step4;
            printf("Step 4 method override active: %d (%s)\n",
                   dixon_global_method_step4, dixon_det_method_name(dixon_global_method_step4));
        }
        
        compute_fq_coefficient_matrix_det(result, coeff_matrix, matrix_size,
                                         npars, polys[0].ctx, coeff_method, res_deg_bound);
        
        if (result->nterms <= 100) {
            fq_mvpoly_print_with_names(result, "Final Resultant", NULL, par_names, gen_name, 0);
        } else {
            printf("Final resultant too large to display (%ld terms)\n", result->nterms);
        }
        fq_mvpoly_make_monic(result);
        
        for (slong i = 0; i < matrix_size; i++) {
            for (slong j = 0; j < matrix_size; j++) {
                fq_mvpoly_clear(&coeff_matrix[i][j]);
            }
            flint_free(coeff_matrix[i]);
        }
        flint_free(coeff_matrix);
    } else {
        fq_mvpoly_init(result, 0, npars, polys[0].ctx);
        printf("Warning: Empty coefficient matrix, resultant is 0\n");
    }
    
    flint_free(row_indices);
    flint_free(col_indices);
    
    fq_mvpoly_clear(&d_poly);
    
    for (slong i = 0; i <= nvars; i++) {
        for (slong j = 0; j <= nvars; j++) {
            fq_mvpoly_clear(&M_mvpoly[i][j]);
            fq_mvpoly_clear(&modified_M_mvpoly[i][j]);
        }
        flint_free(M_mvpoly[i]);
        flint_free(modified_M_mvpoly[i]);
    }
    flint_free(M_mvpoly);
    flint_free(modified_M_mvpoly);

    printf("\n=== Dixon Resultant Computation Complete ===\n");
}

