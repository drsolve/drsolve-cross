/* unified_mpoly_det.c - Implementation of unified polynomial matrix determinant computation */

#include "unified_mpoly_det.h"
#include "dixon_flint.h"

static void debug_print_unified_mpoly_generic(const unified_mpoly_t poly)
{
    slong nvars = (poly != NULL && poly->ctx_ptr != NULL) ? poly->ctx_ptr->nvars : 0;
    if (nvars <= 0) {
        const char *fallback[] = {"x"};
        unified_mpoly_print_pretty(poly, fallback);
        return;
    }

    char **vars = (char **) malloc((size_t) nvars * sizeof(char *));
    if (vars == NULL) {
        printf("<debug-print allocation failed>");
        return;
    }
    for (slong i = 0; i < nvars; i++) {
        vars[i] = (char *) malloc(24);
        if (vars[i] == NULL) {
            for (slong j = 0; j < i; j++) free(vars[j]);
            free(vars);
            printf("<debug-print allocation failed>");
            return;
        }
        snprintf(vars[i], 24, "x%ld", i + 1);
    }

    unified_mpoly_print_pretty(poly, (const char **) vars);

    for (slong i = 0; i < nvars; i++) free(vars[i]);
    free(vars);
}

static int debug_perm_sign(const slong *perm, slong size)
{
    int sign = 1;
    for (slong i = 0; i < size; i++) {
        for (slong j = i + 1; j < size; j++) {
            if (perm[i] > perm[j]) sign = -sign;
        }
    }
    return sign;
}

static slong debug_small_factorial(slong n)
{
    slong r = 1;
    for (slong i = 2; i <= n; i++) r *= i;
    return r;
}

static void debug_dump_leibniz_paths_dfs(unified_mpoly_t **mpoly_matrix,
                                         slong size,
                                         unified_mpoly_ctx_t ctx,
                                         slong row,
                                         slong *perm,
                                         int *used,
                                         slong *path_index)
{
    if (row == size) {
        unified_mpoly_t prefix = unified_mpoly_init(ctx);
        unified_mpoly_t full = unified_mpoly_init(ctx);
        slong prefix_terms;
        slong last_terms;
        int sign = debug_perm_sign(perm, size);

        unified_mpoly_one(prefix);
        for (slong i = 0; i < size - 1; i++) {
            unified_mpoly_mul(prefix, prefix, mpoly_matrix[i][perm[i]]);
        }
        prefix_terms = unified_mpoly_length(prefix);
        last_terms = unified_mpoly_length(mpoly_matrix[size - 1][perm[size - 1]]);

        unified_mpoly_mul(full, prefix, mpoly_matrix[size - 1][perm[size - 1]]);
        if (sign < 0) {
            unified_mpoly_neg(full, full);
        }

        printf("[direct-det path %ld/%ld] perm=[",
               *path_index + 1, debug_small_factorial(size));
        for (slong i = 0; i < size; i++) {
            printf("%ld%s", perm[i], (i + 1 < size) ? "," : "");
        }
        printf("] sign=%c\n", sign > 0 ? '+' : '-');
        printf("  final multiply operand term counts: left=%ld, right=%ld\n",
               prefix_terms, last_terms);
        printf("  path result before summation: terms=%ld\n",
               unified_mpoly_length(full));

        unified_mpoly_clear(prefix);
        unified_mpoly_clear(full);
        (*path_index)++;
        return;
    }

    for (slong col = 0; col < size; col++) {
        if (used[col]) continue;
        used[col] = 1;
        perm[row] = col;
        debug_dump_leibniz_paths_dfs(mpoly_matrix, size, ctx, row + 1, perm, used, path_index);
        used[col] = 0;
    }
}

static void debug_dump_leibniz_paths(unified_mpoly_t **mpoly_matrix,
                                     slong size,
                                     unified_mpoly_ctx_t ctx)
{
    slong *perm = NULL;
    int *used = NULL;
    slong path_index = 0;

    if (g_dixon_verbose_level < 3 || size > 5 || size <= 0) {
        return;
    }

    perm = (slong *) calloc((size_t) size, sizeof(slong));
    used = (int *) calloc((size_t) size, sizeof(int));
    if (perm == NULL || used == NULL) {
        free(perm);
        free(used);
        return;
    }

    printf("\n=== Direct determinant Leibniz path dump (size=%ld) ===\n", size);
    debug_dump_leibniz_paths_dfs(mpoly_matrix, size, ctx, 0, perm, used, &path_index);
    printf("=== End direct determinant Leibniz path dump ===\n");

    free(perm);
    free(used);
}

/* ============================================================================
   RECURSIVE DETERMINANT COMPUTATION WITH UNIFIED INTERFACE
   ============================================================================ */

/* Hand-optimized 3x3 determinant for unified_mpoly */
static void compute_det_3x3_unified(unified_mpoly_t det, 
                                   unified_mpoly_t **m,
                                   unified_mpoly_ctx_t ctx) {
    unified_mpoly_t t1, t2, t3, t4, t5, t6, sum;
    /*
    const char *vars[] = {"x", "y", "z", "u", "v", "w"};
    for (int i = 0; i < 3; i++) {
        printf("[ ");
        for (int j = 0; j < 3; j++) {
            unified_mpoly_print_pretty(m[i][j], vars);
            if (j < 2) printf(",\n  ");
        }
        printf(" ]");
        if (i < 2) printf(",\n");
        printf("\n");
    }
    */
    
    /* Initialize temporaries */
    t1 = unified_mpoly_init(ctx);
    t2 = unified_mpoly_init(ctx);
    t3 = unified_mpoly_init(ctx);
    t4 = unified_mpoly_init(ctx);
    t5 = unified_mpoly_init(ctx);
    t6 = unified_mpoly_init(ctx);
    sum = unified_mpoly_init(ctx);
    
    /* Compute 6 products in parallel if beneficial */
    #ifdef _OPENMP
    #pragma omp parallel sections if(omp_get_max_threads() > 2)
    {
        #pragma omp section
        {
            unified_mpoly_mul(t1, m[1][1], m[2][2]);
            unified_mpoly_mul(t1, m[0][0], t1);
        }
        #pragma omp section
        {
            unified_mpoly_mul(t2, m[1][2], m[2][0]);
            unified_mpoly_mul(t2, m[0][1], t2);
        }
        #pragma omp section
        {
            unified_mpoly_mul(t3, m[1][0], m[2][1]);
            unified_mpoly_mul(t3, m[0][2], t3);
        }
        #pragma omp section
        {
            unified_mpoly_mul(t4, m[1][0], m[2][2]);
            unified_mpoly_mul(t4, m[0][1], t4);
        }
        #pragma omp section
        {
            unified_mpoly_mul(t5, m[1][1], m[2][0]);
            unified_mpoly_mul(t5, m[0][2], t5);
        }
        #pragma omp section
        {
            unified_mpoly_mul(t6, m[1][2], m[2][1]);
            unified_mpoly_mul(t6, m[0][0], t6);
        }
    }
    #else
    /* Sequential computation when OpenMP is not available */
    unified_mpoly_mul(t1, m[1][1], m[2][2]);
    unified_mpoly_mul(t1, m[0][0], t1);
    
    unified_mpoly_mul(t2, m[1][2], m[2][0]);
    unified_mpoly_mul(t2, m[0][1], t2);
    
    unified_mpoly_mul(t3, m[1][0], m[2][1]);
    unified_mpoly_mul(t3, m[0][2], t3);
    
    unified_mpoly_mul(t4, m[1][0], m[2][2]);
    unified_mpoly_mul(t4, m[0][1], t4);
    
    unified_mpoly_mul(t5, m[1][1], m[2][0]);
    unified_mpoly_mul(t5, m[0][2], t5);
    
    unified_mpoly_mul(t6, m[1][2], m[2][1]);
    unified_mpoly_mul(t6, m[0][0], t6);
    #endif
    
    /* Sum with signs: det = t1 + t2 + t3 - t4 - t5 - t6 */
    unified_mpoly_add(sum, t1, t2);
    unified_mpoly_add(sum, sum, t3);
    unified_mpoly_sub(sum, sum, t4);
    unified_mpoly_sub(sum, sum, t5);
    unified_mpoly_sub(det, sum, t6);
    
    /* Cleanup */
    unified_mpoly_clear(t1);
    unified_mpoly_clear(t2);
    unified_mpoly_clear(t3);
    unified_mpoly_clear(t4);
    unified_mpoly_clear(t5);
    unified_mpoly_clear(t6);
    unified_mpoly_clear(sum);
}

static int unified_mpoly_find_nonzero_pivot(unified_mpoly_t **mat,
                                            slong size,
                                            slong start)
{
    for (slong i = start; i < size; i++) {
        for (slong j = start; j < size; j++) {
            if (!unified_mpoly_is_zero(mat[i][j])) {
                return (int) i;
            }
        }
    }
    return -1;
}

static int unified_mpoly_find_nonzero_pivot_col(unified_mpoly_t **mat,
                                                slong size,
                                                slong row,
                                                slong start_col)
{
    for (slong j = start_col; j < size; j++) {
        if (!unified_mpoly_is_zero(mat[row][j])) {
            return (int) j;
        }
    }
    return -1;
}

int compute_unified_mpoly_det_bareiss(unified_mpoly_t det_result,
                                     unified_mpoly_t **mpoly_matrix,
                                     slong size,
                                     unified_mpoly_ctx_t ctx) {
    if (size <= 0) {
        unified_mpoly_one(det_result);
        return 1;
    }

    if (size == 1) {
        unified_mpoly_set(det_result, mpoly_matrix[0][0]);
        return 1;
    }

    if (size == 2) {
        unified_mpoly_t ad, bc;
        ad = unified_mpoly_init(ctx);
        bc = unified_mpoly_init(ctx);
        unified_mpoly_mul(ad, mpoly_matrix[0][0], mpoly_matrix[1][1]);
        unified_mpoly_mul(bc, mpoly_matrix[0][1], mpoly_matrix[1][0]);
        unified_mpoly_sub(det_result, ad, bc);
        unified_mpoly_clear(ad);
        unified_mpoly_clear(bc);
        return 1;
    }

    unified_mpoly_t **work = unified_mpoly_mat_init(size, size, ctx);
    unified_mpoly_t prev_pivot = unified_mpoly_init(ctx);
    unified_mpoly_t pivot = unified_mpoly_init(ctx);
    unified_mpoly_t zero_poly = unified_mpoly_init(ctx);
    int sign = 1;
    int ok = 1;

    unified_mpoly_zero(zero_poly);
    unified_mpoly_one(prev_pivot);

    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            unified_mpoly_set(work[i][j], mpoly_matrix[i][j]);
        }
    }

    for (slong k = 0; k < size - 1; k++) {
        slong pivot_row = k;
        slong pivot_col = k;

        if (unified_mpoly_is_zero(work[k][k])) {
            int row_found = unified_mpoly_find_nonzero_pivot(work, size, k);
            if (row_found < 0) {
                unified_mpoly_zero(det_result);
                goto cleanup;
            }
            pivot_row = (slong) row_found;
            pivot_col = (slong) unified_mpoly_find_nonzero_pivot_col(work, size, pivot_row, k);
            if (pivot_col < 0) {
                unified_mpoly_zero(det_result);
                goto cleanup;
            }
            if (pivot_row != k) {
                unified_mpoly_t *tmp_row = work[k];
                work[k] = work[pivot_row];
                work[pivot_row] = tmp_row;
                sign = -sign;
            }
            if (pivot_col != k) {
                for (slong i = 0; i < size; i++) {
                    unified_mpoly_swap(work[i][k], work[i][pivot_col]);
                }
                sign = -sign;
            }
        }

        unified_mpoly_set(pivot, work[k][k]);
        if (unified_mpoly_is_zero(pivot)) {
            unified_mpoly_zero(det_result);
            goto cleanup;
        }

        if (unified_mpoly_is_zero(prev_pivot) && k > 0) {
            unified_mpoly_zero(det_result);
            ok = 0;
            goto cleanup;
        }

        #ifdef _OPENMP
        if ((size - (k + 1)) >= 2 && omp_get_max_threads() > 1) {
            volatile int parallel_failed = 0;
            #pragma omp parallel for collapse(2) schedule(static)
            for (slong i = k + 1; i < size; i++) {
                for (slong j = k + 1; j < size; j++) {
                    unified_mpoly_t t1 = unified_mpoly_init(ctx);
                    unified_mpoly_t t2 = unified_mpoly_init(ctx);
                    unified_mpoly_t numer = unified_mpoly_init(ctx);
                    unified_mpoly_t quotient = unified_mpoly_init(ctx);

                    unified_mpoly_mul(t1, work[i][j], work[k][k]);
                    unified_mpoly_mul(t2, work[i][k], work[k][j]);
                    unified_mpoly_sub(numer, t1, t2);

                    if (k == 0) {
                        unified_mpoly_set(work[i][j], numer);
                    } else {
                        if (!unified_mpoly_divexact(quotient, numer, prev_pivot)) {
                            parallel_failed = 1;
                        } else {
                            unified_mpoly_set(work[i][j], quotient);
                        }
                    }

                    unified_mpoly_clear(t1);
                    unified_mpoly_clear(t2);
                    unified_mpoly_clear(numer);
                    unified_mpoly_clear(quotient);
                }
            }
            if (parallel_failed) {
                ok = 0;
                goto cleanup;
            }
        } else
        #endif
        {
            unified_mpoly_t t1 = unified_mpoly_init(ctx);
            unified_mpoly_t t2 = unified_mpoly_init(ctx);
            unified_mpoly_t numer = unified_mpoly_init(ctx);
            unified_mpoly_t quotient = unified_mpoly_init(ctx);

            for (slong i = k + 1; i < size; i++) {
                for (slong j = k + 1; j < size; j++) {
                    unified_mpoly_mul(t1, work[i][j], work[k][k]);
                    unified_mpoly_mul(t2, work[i][k], work[k][j]);
                    unified_mpoly_sub(numer, t1, t2);

                    if (k == 0) {
                        unified_mpoly_set(work[i][j], numer);
                    } else {
                        if (!unified_mpoly_divexact(quotient, numer, prev_pivot)) {
                            ok = 0;
                            unified_mpoly_clear(t1);
                            unified_mpoly_clear(t2);
                            unified_mpoly_clear(numer);
                            unified_mpoly_clear(quotient);
                            goto cleanup;
                        }
                        unified_mpoly_set(work[i][j], quotient);
                    }
                }
            }

            unified_mpoly_clear(t1);
            unified_mpoly_clear(t2);
            unified_mpoly_clear(numer);
            unified_mpoly_clear(quotient);
        }

        for (slong i = k + 1; i < size; i++) {
            unified_mpoly_zero(work[i][k]);
            unified_mpoly_zero(work[k][i]);
        }

        unified_mpoly_set(prev_pivot, pivot);
    }

    unified_mpoly_set(det_result, work[size - 1][size - 1]);
    if (sign < 0) {
        unified_mpoly_sub(det_result, zero_poly, det_result);
    }

cleanup:
    unified_mpoly_clear(prev_pivot);
    unified_mpoly_clear(pivot);
    unified_mpoly_clear(zero_poly);
    unified_mpoly_mat_clear(work, size, size);
    return ok;
}

/* Recursive determinant computation using unified_mpoly interface */
void compute_unified_mpoly_det_recursive(unified_mpoly_t det_result, 
                                        unified_mpoly_t **mpoly_matrix, 
                                        slong size, 
                                        unified_mpoly_ctx_t ctx) {
    
    static int recursion_depth = 0;
    recursion_depth++;
    
    if (size <= 0) {
        unified_mpoly_one(det_result);
        recursion_depth--;
        return;
    }
    
    if (size == 1) {
        unified_mpoly_set(det_result, mpoly_matrix[0][0]);
        recursion_depth--;
        return;
    }
    
    if (size == 2) {
        unified_mpoly_t ad, bc;
        ad = unified_mpoly_init(ctx);
        bc = unified_mpoly_init(ctx);
        
        /* det = a*d - b*c */
        unified_mpoly_mul(ad, mpoly_matrix[0][0], mpoly_matrix[1][1]);
        unified_mpoly_mul(bc, mpoly_matrix[0][1], mpoly_matrix[1][0]);
        unified_mpoly_sub(det_result, ad, bc);
        
        unified_mpoly_clear(ad);
        unified_mpoly_clear(bc);
        recursion_depth--;
        return;
    }
    
    if (size == 3) {
        compute_det_3x3_unified(det_result, mpoly_matrix, ctx);
        recursion_depth--;
        return;
    }
    
    /* General case: Laplace expansion along first row */
    unified_mpoly_zero(det_result);
    
    unified_mpoly_t temp_result, cofactor, subdet;
    temp_result = unified_mpoly_init(ctx);
    cofactor = unified_mpoly_init(ctx);
    subdet = unified_mpoly_init(ctx);
    
    /* Allocate submatrix */
    unified_mpoly_t **submatrix = (unified_mpoly_t**) malloc((size-1) * sizeof(unified_mpoly_t*));
    for (slong i = 0; i < size-1; i++) {
        submatrix[i] = (unified_mpoly_t*) malloc((size-1) * sizeof(unified_mpoly_t));
        for (slong j = 0; j < size-1; j++) {
            submatrix[i][j] = unified_mpoly_init(ctx);
        }
    }
    
    /* Laplace expansion along the first row */
    for (slong col = 0; col < size; col++) {
        /* Skip if the element is zero */
        int is_zero = unified_mpoly_is_zero(mpoly_matrix[0][col]);
        
        if (is_zero) {
            continue;
        }
        
        /* Build submatrix by removing row 0 and column col */
        for (slong i = 1; i < size; i++) {
            slong sub_j = 0;
            for (slong j = 0; j < size; j++) {
                if (j != col) {
                    unified_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j]);
                    sub_j++;
                }
            }
        }
        
        /* Recursive computation of minor */
        compute_unified_mpoly_det_recursive(subdet, submatrix, size-1, ctx);
        
        /* Compute cofactor = (-1)^col * matrix[0][col] * subdet */
        unified_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet);
        
        /* Add or subtract based on sign */
        if (col % 2 == 0) {
            unified_mpoly_add(temp_result, det_result, cofactor);
        } else {
            unified_mpoly_sub(temp_result, det_result, cofactor);
        }
        unified_mpoly_set(det_result, temp_result);
    }
    
    /* Cleanup submatrix */
    for (slong i = 0; i < size-1; i++) {
        for (slong j = 0; j < size-1; j++) {
            unified_mpoly_clear(submatrix[i][j]);
        }
        free(submatrix[i]);
    }
    free(submatrix);
    
    unified_mpoly_clear(temp_result);
    unified_mpoly_clear(cofactor);
    unified_mpoly_clear(subdet);

    recursion_depth--;
}

/* ============================================================================
   PARALLEL DETERMINANT COMPUTATION WITH NESTED PARALLELISM
   ============================================================================ */

/* Initialize nested parallelism */
void init_nested_parallelism(void) {
    #ifdef _OPENMP
    omp_set_nested(1);  /* Enable nested parallelism */
    omp_set_max_active_levels(2);  /* Allow 2 levels of parallelism */
    #endif
}

/* Parallel determinant computation with proper nested parallelism */
void compute_unified_mpoly_det_parallel(unified_mpoly_t det_result, 
                                       unified_mpoly_t **mpoly_matrix, 
                                       slong size, 
                                       unified_mpoly_ctx_t ctx,
                                       slong depth) {
    /* For deep recursion or small matrices, use sequential */
    if (size < PARALLEL_THRESHOLD || depth >= MAX_PARALLEL_DEPTH) {
        compute_unified_mpoly_det_recursive(det_result, mpoly_matrix, size, ctx);
        return;
    }
    
    if (size <= 3) {
        compute_unified_mpoly_det_recursive(det_result, mpoly_matrix, size, ctx);
        return;
    }
    
    unified_mpoly_zero(det_result);
    
    /* Count non-zero entries in first row */
    slong nonzero_count = 0;
    for (slong col = 0; col < size; col++) {
        if (!unified_mpoly_is_zero(mpoly_matrix[0][col])) {
            nonzero_count++;
        }
    }
    
    if (nonzero_count < 2) {
        compute_unified_mpoly_det_recursive(det_result, mpoly_matrix, size, ctx);
        return;
    }
    
    /* Allocate space for partial results */
    unified_mpoly_t *partial_results = (unified_mpoly_t*) malloc(size * sizeof(unified_mpoly_t));
    for (slong i = 0; i < size; i++) {
        partial_results[i] = unified_mpoly_init(ctx);
        unified_mpoly_zero(partial_results[i]);
    }
    
    /* Determine parallelism strategy based on depth */
    #ifdef _OPENMP
    if (depth == 0) {
        /* First level: use parallel for with nested parallelism enabled */
        #pragma omp parallel for schedule(static) num_threads(FLINT_MIN(nonzero_count, omp_get_max_threads()))
        for (slong col = 0; col < size; col++) {
            if (unified_mpoly_is_zero(mpoly_matrix[0][col])) {
                continue;
            }
            
            unified_mpoly_t cofactor = unified_mpoly_init(ctx);
            unified_mpoly_t subdet = unified_mpoly_init(ctx);
            
            /* Create submatrix */
            unified_mpoly_t **submatrix = (unified_mpoly_t**) malloc((size-1) * sizeof(unified_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (unified_mpoly_t*) malloc((size-1) * sizeof(unified_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    submatrix[i][j] = unified_mpoly_init(ctx);
                }
            }
            
            /* Fill submatrix */
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        unified_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j]);
                        sub_j++;
                    }
                }
            }
            
            /* Recursive call - this will use nested parallelism at depth 1 */
            compute_unified_mpoly_det_parallel(subdet, submatrix, size-1, ctx, depth+1);
            
            /* Compute cofactor */
            unified_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet);
            
            /* Store with sign */
            if (col % 2 == 0) {
                unified_mpoly_set(partial_results[col], cofactor);
            } else {
                unified_mpoly_neg(partial_results[col], cofactor);
            }
            
            /* Cleanup */
            for (slong i = 0; i < size-1; i++) {
                for (slong j = 0; j < size-1; j++) {
                    unified_mpoly_clear(submatrix[i][j]);
                }
                free(submatrix[i]);
            }
            free(submatrix);
            
            unified_mpoly_clear(cofactor);
            unified_mpoly_clear(subdet);
        }
    } else if (depth == 1 && size >= PARALLEL_THRESHOLD) {
        /* Second level: also use parallel for, but with fewer threads */
        slong max_threads_level2 = FLINT_MAX(1, omp_get_max_threads() / size);
        
        #pragma omp parallel for schedule(static) num_threads(FLINT_MIN(nonzero_count, max_threads_level2)) if(nonzero_count >= 3)
        for (slong col = 0; col < size; col++) {
            if (unified_mpoly_is_zero(mpoly_matrix[0][col])) {
                continue;
            }
            
            unified_mpoly_t cofactor = unified_mpoly_init(ctx);
            unified_mpoly_t subdet = unified_mpoly_init(ctx);
            
            /* Create submatrix */
            unified_mpoly_t **submatrix = (unified_mpoly_t**) malloc((size-1) * sizeof(unified_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (unified_mpoly_t*) malloc((size-1) * sizeof(unified_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    submatrix[i][j] = unified_mpoly_init(ctx);
                }
            }
            
            /* Fill submatrix */
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        unified_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j]);
                        sub_j++;
                    }
                }
            }
            
            /* At depth > 1, use sequential computation */
            compute_unified_mpoly_det_recursive(subdet, submatrix, size-1, ctx);
            
            /* Compute cofactor */
            unified_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet);
            
            /* Store with sign */
            if (col % 2 == 0) {
                unified_mpoly_set(partial_results[col], cofactor);
            } else {
                unified_mpoly_neg(partial_results[col], cofactor);
            }
            
            /* Cleanup */
            for (slong i = 0; i < size-1; i++) {
                for (slong j = 0; j < size-1; j++) {
                    unified_mpoly_clear(submatrix[i][j]);
                }
                free(submatrix[i]);
            }
            free(submatrix);
            
            unified_mpoly_clear(cofactor);
            unified_mpoly_clear(subdet);
        }
    } else
    #endif
    {
        /* Sequential fallback for deeper levels or when OpenMP is not available */
        for (slong col = 0; col < size; col++) {
            if (unified_mpoly_is_zero(mpoly_matrix[0][col])) {
                continue;
            }
            
            unified_mpoly_t cofactor = unified_mpoly_init(ctx);
            unified_mpoly_t subdet = unified_mpoly_init(ctx);
            
            /* Create and fill submatrix */
            unified_mpoly_t **submatrix = (unified_mpoly_t**) malloc((size-1) * sizeof(unified_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (unified_mpoly_t*) malloc((size-1) * sizeof(unified_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    submatrix[i][j] = unified_mpoly_init(ctx);
                }
            }
            
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        unified_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j]);
                        sub_j++;
                    }
                }
            }
            
            /* Sequential computation */
            compute_unified_mpoly_det_recursive(subdet, submatrix, size-1, ctx);
            
            /* Compute cofactor and store */
            unified_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet);
            if (col % 2 == 0) {
                unified_mpoly_set(partial_results[col], cofactor);
            } else {
                unified_mpoly_neg(partial_results[col], cofactor);
            }
            
            /* Cleanup */
            for (slong i = 0; i < size-1; i++) {
                for (slong j = 0; j < size-1; j++) {
                    unified_mpoly_clear(submatrix[i][j]);
                }
                free(submatrix[i]);
            }
            free(submatrix);
            
            unified_mpoly_clear(cofactor);
            unified_mpoly_clear(subdet);
        }
    }
    
    /* Sum results (sequential to avoid race conditions) */
    unified_mpoly_t temp_sum = unified_mpoly_init(ctx);
    
    for (slong col = 0; col < size; col++) {
        if (!unified_mpoly_is_zero(partial_results[col])) {
            unified_mpoly_add(temp_sum, det_result, partial_results[col]);
            unified_mpoly_set(det_result, temp_sum);
        }
        unified_mpoly_clear(partial_results[col]);
    }
    
    unified_mpoly_clear(temp_sum);
    free(partial_results);
}

/* ============================================================================
   CONVENIENCE WRAPPER FUNCTIONS
   ============================================================================ */

/* Main entry point for determinant computation */
void compute_unified_mpoly_det(unified_mpoly_t det_result,
                              unified_mpoly_t **mpoly_matrix,
                              slong size,
                              unified_mpoly_ctx_t ctx,
                              int use_parallel) {
    compute_unified_mpoly_det_with_method(det_result, mpoly_matrix, size, ctx,
                                          use_parallel, 0);
}

void compute_unified_mpoly_det_with_method(unified_mpoly_t det_result,
                                          unified_mpoly_t **mpoly_matrix,
                                          slong size,
                                          unified_mpoly_ctx_t ctx,
                                          int use_parallel,
                                          int prefer_bareiss) {
    if (!prefer_bareiss && g_dixon_verbose_level >= 3 && size <= 5) {
        debug_dump_leibniz_paths(mpoly_matrix, size, ctx);
    }

    if (prefer_bareiss) {
        if (!compute_unified_mpoly_det_bareiss(det_result, mpoly_matrix, size, ctx)) {
            compute_unified_mpoly_det_recursive(det_result, mpoly_matrix, size, ctx);
        }
        return;
    }
    
    /* Special case: 3x3 matrix with parallel computation available */
    #ifdef _OPENMP
    if (size == 3 && use_parallel && omp_get_max_threads() > 1) {
        compute_det_3x3_unified(det_result, mpoly_matrix, ctx);
        return;
    }
    #endif
    
    if (use_parallel && size >= PARALLEL_THRESHOLD) {
        #ifdef _OPENMP
        static int nested_init = 0;
        if (!nested_init) {
            init_nested_parallelism();
            nested_init = 1;
        }
        compute_unified_mpoly_det_parallel(det_result, mpoly_matrix, size, ctx, 0);
        #else
        compute_unified_mpoly_det_recursive(det_result, mpoly_matrix, size, ctx);
        #endif
    } else {
        compute_unified_mpoly_det_recursive(det_result, mpoly_matrix, size, ctx);
    }
}

/* ============================================================================
   HELPER FUNCTIONS FOR MATRIX ALLOCATION/DEALLOCATION
   ============================================================================ */

/* Allocate a unified_mpoly matrix */
unified_mpoly_t** unified_mpoly_mat_init(slong rows, slong cols, unified_mpoly_ctx_t ctx) {
    unified_mpoly_t **mat = (unified_mpoly_t**) malloc(rows * sizeof(unified_mpoly_t*));
    for (slong i = 0; i < rows; i++) {
        mat[i] = (unified_mpoly_t*) malloc(cols * sizeof(unified_mpoly_t));
        for (slong j = 0; j < cols; j++) {
            mat[i][j] = unified_mpoly_init(ctx);
        }
    }
    return mat;
}

/* Free a unified_mpoly matrix */
void unified_mpoly_mat_clear(unified_mpoly_t **mat, slong rows, slong cols) {
    for (slong i = 0; i < rows; i++) {
        for (slong j = 0; j < cols; j++) {
            unified_mpoly_clear(mat[i][j]);
        }
        free(mat[i]);
    }
    free(mat);
}

/* Print a unified_mpoly matrix */
void unified_mpoly_mat_print_pretty(unified_mpoly_t **mat, slong rows, slong cols,
                                   const char **vars) {
    printf("[\n");
    for (slong i = 0; i < rows; i++) {
        printf("  [");
        for (slong j = 0; j < cols; j++) {
            unified_mpoly_print_pretty(mat[i][j], vars);
            if (j < cols - 1) printf(", ");
        }
        printf("]");
        if (i < rows - 1) printf(",");
        printf("\n");
    }
    printf("]\n");
}

/* ============================================================================
   EXAMPLE/TEST FUNCTIONS
   ============================================================================ */

/* Test the determinant computation */
void test_unified_mpoly_det(void) {
    printf("\n=== Testing Unified Polynomial Matrix Determinant ===\n");
    
    /* Initialize field context for GF(2^8) */
    fq_nmod_ctx_t fq_ctx;
    nmod_poly_t mod;
    nmod_poly_init(mod, 2);
    
    /* Set modulus for GF(2^8): x^8 + x^4 + x^3 + x^2 + 1 */
    nmod_poly_set_coeff_ui(mod, 0, 1);
    nmod_poly_set_coeff_ui(mod, 2, 1);
    nmod_poly_set_coeff_ui(mod, 3, 1);
    nmod_poly_set_coeff_ui(mod, 4, 1);
    nmod_poly_set_coeff_ui(mod, 8, 1);
    
    fq_nmod_ctx_init_modulus(fq_ctx, mod, "a");
    nmod_poly_clear(mod);
    
    /* Create field context */
    field_ctx_t field_ctx;
    field_ctx_init(&field_ctx, fq_ctx);
    
    /* Create multivariate context for 2 variables */
    unified_mpoly_ctx_t ctx = unified_mpoly_ctx_init(2, ORD_LEX, &field_ctx);
    
    /* Test 3x3 determinant */
    slong size = 3;
    unified_mpoly_t **mat = unified_mpoly_mat_init(size, size, ctx);
    
    /* Set up a simple test matrix:
     * [x+1,  y,   1]
     * [y,    x,   0]
     * [1,    0,   x+y]
     */
    field_elem_u one;
    field_set_one(&one, field_ctx.field_id, (void*)field_ctx.ctx.fq_ctx);
    
    ulong exp[2];
    
    /* mat[0][0] = x + 1 */
    exp[0] = 1; exp[1] = 0;  /* x */
    unified_mpoly_set_coeff_ui(mat[0][0], &one, exp);
    exp[0] = 0; exp[1] = 0;  /* 1 */
    unified_mpoly_set_coeff_ui(mat[0][0], &one, exp);
    
    /* mat[0][1] = y */
    exp[0] = 0; exp[1] = 1;  /* y */
    unified_mpoly_set_coeff_ui(mat[0][1], &one, exp);
    
    /* mat[0][2] = 1 */
    exp[0] = 0; exp[1] = 0;  /* 1 */
    unified_mpoly_set_coeff_ui(mat[0][2], &one, exp);
    
    /* mat[1][0] = y */
    exp[0] = 0; exp[1] = 1;  /* y */
    unified_mpoly_set_coeff_ui(mat[1][0], &one, exp);
    
    /* mat[1][1] = x */
    exp[0] = 1; exp[1] = 0;  /* x */
    unified_mpoly_set_coeff_ui(mat[1][1], &one, exp);
    
    /* mat[1][2] = 0 */
    unified_mpoly_zero(mat[1][2]);
    
    /* mat[2][0] = 1 */
    exp[0] = 0; exp[1] = 0;  /* 1 */
    unified_mpoly_set_coeff_ui(mat[2][0], &one, exp);
    
    /* mat[2][1] = 0 */
    unified_mpoly_zero(mat[2][1]);
    
    /* mat[2][2] = x + y */
    exp[0] = 1; exp[1] = 0;  /* x */
    unified_mpoly_set_coeff_ui(mat[2][2], &one, exp);
    exp[0] = 0; exp[1] = 1;  /* y */
    unified_mpoly_set_coeff_ui(mat[2][2], &one, exp);
    
    /* Print the matrix */
    const char *vars[] = {"x", "y"};
    printf("Test matrix:\n");
    unified_mpoly_mat_print_pretty(mat, size, size, vars);
    
    /* Compute determinant */
    unified_mpoly_t det = unified_mpoly_init(ctx);
    
    /* Test sequential computation */
    printf("\nSequential computation:\n");
    clock_t start = clock();
    compute_unified_mpoly_det(det, mat, size, ctx, 0);
    clock_t end = clock();
    double seq_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    printf("Determinant = ");
    unified_mpoly_print_pretty(det, vars);
    printf("\n");
    printf("Time: %.6f seconds\n", seq_time);
    
    /* Test parallel computation (if available) */
    #ifdef _OPENMP
    printf("\nParallel computation:\n");
    unified_mpoly_t det_parallel = unified_mpoly_init(ctx);
    
    start = clock();
    compute_unified_mpoly_det(det_parallel, mat, size, ctx, 1);
    end = clock();
    double par_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    printf("Determinant = ");
    unified_mpoly_print_pretty(det_parallel, vars);
    printf("\n");
    printf("Time: %.6f seconds\n", par_time);
    
    /* Verify results match */
    unified_mpoly_sub(det_parallel, det_parallel, det);
    if (unified_mpoly_is_zero(det_parallel)) {
        printf("Results match!\n");
    } else {
        printf("ERROR: Results don't match!\n");
    }
    
    unified_mpoly_clear(det_parallel);
    #endif
    
    /* Cleanup */
    unified_mpoly_clear(det);
    unified_mpoly_mat_clear(mat, size, size);
    unified_mpoly_ctx_clear(ctx);
    fq_nmod_ctx_clear(fq_ctx);
    
    /* Test with larger matrix */
    printf("\n=== Testing with larger 5x5 matrix ===\n");
    
    /* Reinitialize for prime field test */
    fq_nmod_ctx_init_ui(fq_ctx, 101, 1, "x");
    field_ctx_init(&field_ctx, fq_ctx);
    ctx = unified_mpoly_ctx_init(2, ORD_LEX, &field_ctx);
    
    size = 5;
    mat = unified_mpoly_mat_init(size, size, ctx);
    
    /* Create a random-ish matrix */
    flint_rand_t state;
    flint_rand_init(state);
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            /* Add some random terms */
            field_elem_u coeff;
            for (int k = 0; k < 3; k++) {
                exp[0] = n_randint(state, 3);
                exp[1] = n_randint(state, 3);
                coeff.nmod = n_randint(state, 100) + 1;
                unified_mpoly_set_coeff_ui(mat[i][j], &coeff, exp);
            }
        }
    }
    
    printf("Computing 5x5 determinant...\n");
    det = unified_mpoly_init(ctx);
    
    start = clock();
    compute_unified_mpoly_det(det, mat, size, ctx, 1);  /* Use parallel if available */
    end = clock();
    double time_5x5 = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    printf("Determinant computed in %.6f seconds\n", time_5x5);
    printf("Result has %ld terms\n", unified_mpoly_length(det));
    
    /* Cleanup */
    flint_rand_clear(state);
    unified_mpoly_clear(det);
    unified_mpoly_mat_clear(mat, size, size);
    unified_mpoly_ctx_clear(ctx);
    fq_nmod_ctx_clear(fq_ctx);
}
