/* fq_poly_mat_det.c - Fixed Matrix determinant implementation */

#include "fq_poly_mat_det.h"

/* ============================================================================
   GLOBAL VARIABLES FOR PROGRESS TRACKING
   ============================================================================ */

int _is_first_step = 0;
int _show_progress = 0;
int g_show_progress = 0;

/* ============================================================================
   UNIFIED POLYNOMIAL MATRIX OPERATIONS - IMPLEMENTATIONS
   ============================================================================ */

/* Swap rows in unified poly mat */
void unified_poly_mat_swap_rows_ptr(unified_poly_mat_struct *mat, slong r, slong s) {
    if (r != s) {
        unified_poly_struct *tmp = mat->rows[r];
        mat->rows[r] = mat->rows[s];
        mat->rows[s] = tmp;
    }
}

/* Rotate rows upward */
void unified_poly_mat_rotate_rows_upward_ptr(unified_poly_mat_struct *mat, slong i, slong j) {
    if (i != j) {
        unified_poly_struct *tmp_mat = mat->rows[i];
        for (slong ii = i; ii < j; ii++) {
            mat->rows[ii] = mat->rows[ii+1];
        }
        mat->rows[j] = tmp_mat;
    }
}

/* Permute rows */
void unified_poly_mat_permute_rows_ptr(unified_poly_mat_struct *mat, const slong *perm) {
    unified_poly_struct **new_rows = (unified_poly_struct **)malloc(mat->r * sizeof(unified_poly_struct *));
    
    for (slong i = 0; i < mat->r; i++) {
        new_rows[i] = mat->rows[perm[i]];
    }
    
    memcpy(mat->rows, new_rows, mat->r * sizeof(unified_poly_struct *));
    free(new_rows);
}

/* Window resize columns */
void unified_poly_mat_window_resize_columns_ptr(unified_poly_mat_struct *mat, slong c) {
    if (c >= 0)
        mat->c = c;
    else
        mat->c = mat->c + c + 1;
}

/* ============================================================================
   PIVOT COMPUTATION - IMPLEMENTATIONS
   ============================================================================ */

void unified_poly_vec_pivot_profile(slong *pivind, slong *pivdeg,
                                   unified_poly_struct *vec,
                                   slong len) {
    slong max = -1;
    slong piv = -1;
    
    for (slong j = 0; j < len; j++) {
        slong d = unified_poly_degree(&vec[j]);
        if (d >= 0 && d > max) {
            max = d;
            piv = j;
        }
    }
    
    *pivdeg = max;
    *pivind = (piv == -1) ? len : piv;
}

/* ============================================================================
   COLLISION SOLVING - FULLY UNIFIED WITH COMPILE-TIME DISPATCH
   ============================================================================ */

void atomic_solve_pivot_collision_unified_optimized(unified_poly_mat_struct *mat,
                                                   unified_poly_mat_struct *other,
                                                   slong pi1, slong pi2, slong j) {
    unified_poly_struct *entry1 = mat->rows[pi1] + j;
    unified_poly_struct *entry2 = mat->rows[pi2] + j;
    
    slong deg1 = unified_poly_degree_fast(entry1);
    slong deg2 = unified_poly_degree_fast(entry2);
    slong exp = deg1 - deg2;
    
    if (exp < 0) return;
    
    /* Ensure workspace is initialized for THIS field context */
    ensure_workspace_initialized(mat->ctx);
    
    /* Get leading coefficients directly */
    field_elem_u *lc1_ptr = unified_poly_get_coeff_ptr(entry1, deg1);
    field_elem_u *lc2_ptr = unified_poly_get_coeff_ptr(entry2, deg2);
    
    if (!lc1_ptr || !lc2_ptr) return;
    
    /* Use workspace elements */
    field_elem_u *cst = &g_unified_workspace.cst;
    
    /* Get context pointer */
    void *ctx_ptr = (mat->ctx->field_id == FIELD_ID_NMOD) ? 
                   (void*)&mat->ctx->ctx.nmod_ctx : 
                   (void*)mat->ctx->ctx.fq_ctx;
    
    /* Compute the constant: cst = lc1/lc2 for all fields */
    /* No special case needed - the unified interface handles field differences */
    if (field_is_one(lc2_ptr, mat->ctx->field_id, ctx_ptr)) {
        field_set_elem(cst, lc1_ptr, mat->ctx->field_id, ctx_ptr);
    } else {
        field_inv(&g_unified_workspace.inv, lc2_ptr, mat->ctx->field_id, ctx_ptr);
        field_mul(cst, lc1_ptr, &g_unified_workspace.inv, mat->ctx->field_id, ctx_ptr);
    }
    
    field_neg(cst, cst, mat->ctx->field_id, ctx_ptr);
    
    /* Process columns */
    for (slong jj = 0; jj < mat->c; jj++) {
        unified_poly_struct *e1 = mat->rows[pi1] + jj;
        unified_poly_struct *e2 = mat->rows[pi2] + jj;
        
        if (unified_poly_is_zero(e2)) continue;
        
        /* Always use multiply and add approach for correctness */
        unified_poly_scalar_mul(&g_unified_workspace.tmp, e2, cst);
        unified_poly_shift_left_add_inplace(e1, &g_unified_workspace.tmp, exp, mat->ctx);
        unified_poly_normalise(e1);
    }
    
    /* Process other matrix if provided */
    if (other) {
        for (slong jj = 0; jj < other->c; jj++) {
            unified_poly_struct *e1 = other->rows[pi1] + jj;
            unified_poly_struct *e2 = other->rows[pi2] + jj;
            
            if (unified_poly_is_zero(e2)) continue;
            
            unified_poly_scalar_mul(&g_unified_workspace.tmp, e2, cst);
            unified_poly_shift_left_add_inplace(e1, &g_unified_workspace.tmp, exp, other->ctx);
            unified_poly_normalise(e1);
        }
    }
}

/* ============================================================================
   WEAK POPOV FORM COMPUTATION - UNIFIED VERSION
   ============================================================================ */

slong unified_poly_mat_weak_popov_iter_submat_rowbyrow_optimized(
    unified_poly_mat_struct *mat,
    unified_poly_mat_struct *tsf,
    slong *det,
    slong *pivind,
    slong *rrp,
    slong rstart,
    slong cstart,
    slong rdim,
    slong cdim,
    slong early_exit_zr) {
    
    if (rdim == 0 || cdim == 0)
        return 0;
    
    const slong zpiv = cdim;
    slong rk = 0;
    slong zr = 0;
    
    /* Pre-allocate pivot tracking */
    slong *pivot_row = (slong *)calloc(cdim, sizeof(slong));
    for (slong i = 0; i < cdim; i++) {
        pivot_row[i] = -1L;
    }
    
    /* Don't use degree cache for general finite fields - it can become stale */
    /* The overhead of recomputing degrees is worth the correctness */
    
    slong pivdeg;
    clock_t row_start_time = clock();
    
    while (rk + zr < rdim && zr < early_exit_zr) {
        clock_t current_time = clock();
        double elapsed = ((double)(current_time - row_start_time)) / CLOCKS_PER_SEC;
        
        /* Progress reporting 
        if (_show_progress && _is_first_step && rdim >= 10) {
            printf("\r  Row %ld/%ld | Rank: %ld | Zero rows: %ld | Time: %.1fs | Rate: %.1f rows/s",
                   rk + zr + 1, rdim, rk, zr, elapsed,
                   (rk + zr) > 0 ? (rk + zr) / elapsed : 0.0);
            //fflush(stdout);
        }
        */
        /* Always recompute pivot - no caching for correctness */
        slong max_deg = -1;
        slong max_idx = zpiv;
        
        /* Ensure all polynomials are normalized before degree computation */
        for (slong j = 0; j < cdim; j++) {
            unified_poly_normalise(mat->rows[rstart + rk] + cstart + j);
        }
        
        unified_poly_vec_pivot_profile(&max_idx, &max_deg,
                                      mat->rows[rstart + rk] + cstart, cdim);
        
        pivind[rk] = max_idx;
        pivdeg = max_deg;
        
        if (pivind[rk] == zpiv) {
            /* Zero row */
            unified_poly_mat_rotate_rows_upward_ptr(mat, rstart + rk, rstart + rdim - 1);
            if (tsf) {
                unified_poly_mat_rotate_rows_upward_ptr(tsf, rstart + rk, rstart + rdim - 1);
            }
            
            if (det && (rdim - 1 - rk) % 2 == 1) {
                *det = -*det;
            }
            zr++;
        } else if (pivot_row[pivind[rk]] == -1) {
            /* New pivot */
            pivot_row[pivind[rk]] = rk;
            if (rrp) {
                rrp[rk] = rk + zr;
            }
            rk++;
        } else {
            /* Pivot collision */
            slong pi = pivot_row[pivind[rk]];
            
            /* Ensure pivot polynomial is normalized */
            unified_poly_normalise(mat->rows[rstart + pi] + cstart + pivind[rk]);
            unified_poly_struct *entry_pi = mat->rows[rstart + pi] + cstart + pivind[rk];
            slong deg_pi = unified_poly_degree_fast(entry_pi);
            
            if (pivdeg < deg_pi) {
                unified_poly_mat_swap_rows_ptr(mat, rstart + pi, rstart + rk);
                if (tsf) {
                    unified_poly_mat_swap_rows_ptr(tsf, rstart + pi, rstart + rk);
                }
                
                if (det && pi != rk) {
                    *det = -*det;
                }
            }
            
            /* Use collision solver */
            atomic_solve_pivot_collision_unified_optimized(mat, tsf, rstart + rk, rstart + pi, 
                                                          cstart + pivind[rk]);
            
            /* Ensure all modified polynomials are normalized */
            for (slong j = 0; j < mat->c; j++) {
                unified_poly_normalise(mat->rows[rstart + rk] + j);
            }
        }
    }
    
    free(pivot_row);
    
    if (zr >= early_exit_zr)
        return -rk;
    else
        return rk;
}

/* ============================================================================
   HELPER FUNCTIONS - IMPLEMENTATIONS
   ============================================================================ */

/* Permutation functions */
int slong_pair_compare(const void *a, const void *b) {
    slong_pair aa = *(const slong_pair *)a;
    slong_pair bb = *(const slong_pair *)b;
    if (aa.value == bb.value)
        return (aa.index < bb.index) ? -1 : ((aa.index > bb.index) ? 1 : 0);
    else
        return (aa.value < bb.value) ? -1 : 1;
}

void vec_sort_permutation(slong *perm, slong *sorted_vec,
                         const slong *vec, slong n, slong_pair *pair_tmp) {
    for (slong i = 0; i < n; i++) {
        pair_tmp[i].value = vec[i];
        pair_tmp[i].index = i;
    }
    
    qsort(pair_tmp, n, sizeof(slong_pair), slong_pair_compare);
    
    for (slong i = 0; i < n; i++)
        perm[i] = pair_tmp[i].index;
    if (sorted_vec)
        for (slong i = 0; i < n; i++)
            sorted_vec[i] = pair_tmp[i].value;
}

void unified_poly_mat_permute_rows_by_sorting_vec_ptr(unified_poly_mat_struct *mat,
                                                     slong r, slong *vec, slong *perm) {
    slong_pair *tmp = (slong_pair *)malloc(r * sizeof(slong_pair));
    vec_sort_permutation(perm, vec, vec, r, tmp);
    for (slong i = r; i < mat->r; i++)
        perm[i] = i;
    free(tmp);
    unified_poly_mat_permute_rows_ptr(mat, perm);
}

/* Check permutation parity */
int perm_parity(const slong *perm, slong n) {
    int parity = 0;
    int *visited = (int *)calloc(n, sizeof(int));
    
    for (slong i = 0; i < n; i++) {
        if (!visited[i]) {
            slong j = i;
            slong cycle_len = 0;
            
            while (!visited[j]) {
                visited[j] = 1;
                j = perm[j];
                cycle_len++;
            }
            
            if (cycle_len % 2 == 0)
                parity ^= 1;
        }
    }
    
    free(visited);
    return parity;
}

#ifdef HAVE_PML
static ulong dixonres_nmod_poly_mat_eval_det_at(const nmod_poly_mat_t mat, ulong point) {
    nmod_mat_t eval_mat;
    ulong det_eval;

    nmod_mat_init(eval_mat, mat->r, mat->c, mat->modulus);
    nmod_poly_mat_evaluate_nmod(eval_mat, mat, point);
    det_eval = nmod_mat_det(eval_mat);
    nmod_mat_clear(eval_mat);

    return det_eval;
}

static int dixonres_nmod_poly_mat_det_generic_exact(nmod_poly_t det,
                                                    const nmod_poly_mat_t mat) {
    if (!nmod_poly_mat_det_generic(det, mat))
        return 0;

    if (nmod_poly_is_zero(det))
        return 0;

    {
        const ulong modulus = mat->modulus;
        const ulong max_trials = (modulus < UWORD(32)) ? modulus : UWORD(32);
        ulong scale = UWORD(1);
        int have_scale = 0;
        slong good_points = 0;

        for (ulong point = 0; point < max_trials; point++) {
            const ulong det_poly = nmod_poly_evaluate_nmod(det, point);
            const ulong det_eval = dixonres_nmod_poly_mat_eval_det_at(mat, point);

            if (det_poly == UWORD(0)) {
                if (det_eval != UWORD(0))
                    return 0;
                continue;
            }

            if (det_eval == UWORD(0))
                return 0;

            if (!have_scale) {
                scale = nmod_mul(det_eval,
                                 n_invmod(det_poly, modulus),
                                 det->mod);
                have_scale = 1;
            } else if (nmod_mul(det_poly, scale, det->mod) != det_eval) {
                return 0;
            }

            ++good_points;
            if (good_points >= 2)
                break;
        }

        if (!have_scale || good_points < 2)
            return 0;

        if (scale != UWORD(1))
            nmod_poly_scalar_mul_nmod(det, det, scale);
    }

    return 1;
}
#endif

/* ============================================================================
   MAIN ENTRY POINT - IMPLEMENTATION
   ============================================================================ */

void fq_nmod_poly_mat_det_iter(fq_nmod_poly_t det,
                              fq_nmod_poly_mat_t mat,
                              const fq_nmod_ctx_t ctx) {
    if (mat->r == 0) {
        fq_nmod_poly_one(det, ctx);
        return;
    }
    
    if (mat->r != mat->c) {
        fq_nmod_poly_zero(det, ctx);
        return;
    }

#ifdef HAVE_PML
    /* Check if we're in a prime field (degree 1) */
    if (fq_nmod_ctx_degree(ctx) == 1) {
        /* Convert to nmod_poly_mat and use the optimized prime field version */
        printf("Using fast determinant backend from PML library\n");

        nmod_poly_mat_t nmod_mat;
        nmod_poly_t nmod_det;
        
        ulong p = fq_nmod_ctx_prime(ctx);
        nmod_poly_mat_init(nmod_mat, mat->r, mat->c, p);
        nmod_poly_init(nmod_det, p);
        
        /* Convert fq_nmod_poly_mat to nmod_poly_mat */
        for (slong i = 0; i < mat->r; i++) {
            for (slong j = 0; j < mat->c; j++) {
                fq_nmod_poly_struct *src = fq_nmod_poly_mat_entry(mat, i, j);
                nmod_poly_struct *dst = nmod_poly_mat_entry(nmod_mat, i, j);
                
                /* Convert each polynomial coefficient */
                slong len = fq_nmod_poly_length(src, ctx);
                nmod_poly_fit_length(dst, len);
                
                for (slong k = 0; k < len; k++) {
                    fq_nmod_t coeff;
                    fq_nmod_init(coeff, ctx);
                    fq_nmod_poly_get_coeff(coeff, src, k, ctx);
                    
                    /* In prime field, fq_nmod element is just a polynomial of degree 0 */
                    ulong val = nmod_poly_get_coeff_ui(coeff, 0);
                    nmod_poly_set_coeff_ui(dst, k, val);
                    
                    fq_nmod_clear(coeff, ctx);
                }
                _nmod_poly_set_length(dst, len);
                _nmod_poly_normalise(dst);
            }
        }
        
        /* Call the optimized nmod version */
        clock_t nmod_start = clock();
        if (!dixonres_nmod_poly_mat_det_generic_exact(nmod_det, nmod_mat)) {
            nmod_poly_mat_det_iter(nmod_det, nmod_mat);
        }
        clock_t nmod_end = clock();
        
        /* Convert result back to fq_nmod_poly */
        fq_nmod_poly_zero(det, ctx);
        slong det_len = nmod_poly_length(nmod_det);
        
        for (slong i = 0; i < det_len; i++) {
            ulong coeff_val = nmod_poly_get_coeff_ui(nmod_det, i);
            if (coeff_val != 0) {
                fq_nmod_t coeff;
                fq_nmod_init(coeff, ctx);
                nmod_poly_set_coeff_ui(coeff, 0, coeff_val);
                fq_nmod_poly_set_coeff(det, i, coeff, ctx);
                fq_nmod_clear(coeff, ctx);
            }
        }
        
        /* Cleanup */
        nmod_poly_mat_clear(nmod_mat);
        nmod_poly_clear(nmod_det);
        
        double nmod_time = ((double)(nmod_end - nmod_start)) / CLOCKS_PER_SEC;
        (void) nmod_time;
        
        return;
    }

#else
    /* No PML optimization available */
    if (fq_nmod_ctx_degree(ctx) == 1) {
        /* fall through to the general algorithm */
    }
#endif
    
    /* Initialize field context with compile-time dispatch */
    field_ctx_t field_ctx;
    field_ctx_init(&field_ctx, ctx);

    /* Initialize workspace for the field context */
    ensure_workspace_initialized(&field_ctx);
    /*
    // Print field information 
    printf("\n=== Field Information ===\n");
    printf("Field size: GF(%lu^%ld)\n", fq_nmod_ctx_prime(ctx), fq_nmod_ctx_degree(ctx));
    printf("Matrix size: %ld x %ld\n", mat->r, mat->c);
    printf("Using compile-time dispatch: %s\n", field_ctx.description);
    */
    clock_t start_time = clock();
    
    /* Convert entire matrix to unified format ONCE */
    //printf("Converting matrix to unified format...\n");
    unified_poly_mat_t unified_mat;
    unified_poly_mat_init(unified_mat, mat->r, mat->c, &field_ctx);
    fq_nmod_poly_mat_to_unified(unified_mat, mat, ctx, &field_ctx);
    
    /* Initialize progress tracking */
    _show_progress = (mat->r >= 50);
    g_show_progress = (mat->r >= 20);
    
    if (_show_progress)
        printf("Computing determinant for %ldx%ld matrix...\n", mat->r, mat->r);
    
    slong udet = 1;
    slong rk;
    slong *pivind = (slong *)malloc(mat->r * sizeof(slong));
    slong *perm = (slong *)malloc(mat->r * sizeof(slong));
    
    /* Create window view as a copy of the matrix struct */
    unified_poly_mat_struct view;
    view.entries = NULL;  /* Indicates this is a window/view */
    view.r = unified_mat->r;
    view.c = unified_mat->c;
    view.rows = unified_mat->rows;
    view.ctx = unified_mat->ctx;
    
    /* Process each dimension from largest to smallest */
    for (slong i = mat->r - 1; i >= 1; i--) {
        _is_first_step = (i == mat->r - 1);
        /*
        if (!_is_first_step && g_show_progress) {
            printf("\rProcessing dimension %ld x %ld... ", i + 1, i);
            //fflush(stdout);
        }
        */
        /* Call weak popov on unified matrix */
        rk = unified_poly_mat_weak_popov_iter_submat_rowbyrow_optimized(
            &view, NULL, &udet, pivind, NULL, 0, 0, i+1, i, 2);
        
        /* Early exit if rank-deficient */
        if (rk < i || unified_poly_is_zero(view.rows[i] + i)) {
            unified_poly_mat_clear(unified_mat);
            free(pivind);
            free(perm);
            fq_nmod_poly_zero(det, ctx);
            
            if (_show_progress || g_show_progress) {
                printf("\r                                                                          \r");
                printf("Matrix is singular at dimension %ld\n", i + 1);
            }
            return;
        }
        
        /* Sort and permute rows */
        unified_poly_mat_permute_rows_by_sorting_vec_ptr(&view, rk, pivind, perm);
        unified_poly_mat_window_resize_columns_ptr(&view, -1);
        
        if (perm_parity(perm, rk))
            udet = -udet;
    }
    
    /* Retrieve determinant as product of diagonal entries */
    unified_poly_mat_window_resize_columns_ptr(&view, mat->r - 1);
    
    /* Compute determinant in unified format */
    unified_poly_t det_unified;
    unified_poly_init(det_unified, &field_ctx);
    unified_poly_set(det_unified, view.rows[0] + 0);
    
    if (unified_poly_is_zero(det_unified)) {
        unified_poly_clear(det_unified);
        unified_poly_mat_clear(unified_mat);
        free(pivind);
        free(perm);
        fq_nmod_poly_zero(det, ctx);
        return;
    }
    
    /* Get context pointer */
    void *ctx_ptr = (field_ctx.field_id == FIELD_ID_NMOD) ? 
                   (void*)&field_ctx.ctx.nmod_ctx : 
                   (void*)field_ctx.ctx.fq_ctx;
    
    if (udet == -1) {
        /* Negate by multiplying by -1 */
        field_elem_u neg_one;
        field_init_elem(&neg_one, field_ctx.field_id, ctx_ptr);
        field_set_one(&neg_one, field_ctx.field_id, ctx_ptr);
        field_neg(&neg_one, &neg_one, field_ctx.field_id, ctx_ptr);
        unified_poly_scalar_mul(det_unified, det_unified, &neg_one);
        field_clear_elem(&neg_one, field_ctx.field_id, ctx_ptr);
    }
    
    /* Multiply by remaining diagonal elements */
    unified_poly_t tmp;
    unified_poly_init(tmp, &field_ctx);
    
    for (slong i = 1; i < view.r; i++) {
        unified_poly_mul(tmp, det_unified, view.rows[i] + i);
        unified_poly_set(det_unified, tmp);
    }
    
    unified_poly_clear(tmp);
    
    /* Convert result back to fq_nmod_poly */
    unified_to_fq_nmod_poly(det, det_unified, ctx, &field_ctx);
    
    /* Cleanup */
    unified_poly_clear(det_unified);
    unified_poly_mat_clear(unified_mat);
    free(pivind);
    free(perm);
    
    clock_t end_time = clock();
    double total_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    
    if (_show_progress || g_show_progress) {
        printf("\r                                                                          \r");
        printf("  Time: %.2f seconds\n", total_time);
    }
}
