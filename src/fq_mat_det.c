/* fq_mat_det.c - Unified Matrix Determinant Implementation */

#include "fq_mat_det.h"

/* ============================================================================
   UNIFIED MATRIX FUNCTIONS - IMPLEMENTATIONS
   ============================================================================ */

void unified_mat_init(unified_mat_t mat, slong rows, slong cols, field_ctx_t *ctx) {
    mat->entries = (field_elem_u *)calloc(rows * cols, sizeof(field_elem_u));
    mat->r = rows;
    mat->c = cols;
    mat->ctx = ctx;
    
    /* Initialize all entries */
    void *ctx_ptr = (ctx->field_id == FIELD_ID_NMOD) ? 
                   (void*)&ctx->ctx.nmod_ctx : 
                   (void*)ctx->ctx.fq_ctx;
    
    for (slong i = 0; i < rows * cols; i++) {
        field_init_elem(&mat->entries[i], ctx->field_id, ctx_ptr);
    }
}

void unified_mat_clear(unified_mat_t mat) {
    void *ctx_ptr = (mat->ctx->field_id == FIELD_ID_NMOD) ? 
                   (void*)&mat->ctx->ctx.nmod_ctx : 
                   (void*)mat->ctx->ctx.fq_ctx;
    
    for (slong i = 0; i < mat->r * mat->c; i++) {
        field_clear_elem(&mat->entries[i], mat->ctx->field_id, ctx_ptr);
    }
    free(mat->entries);
}

field_elem_u* unified_mat_entry(unified_mat_t mat, slong i, slong j) {
    return &mat->entries[i * mat->c + j];
}

const field_elem_u* unified_mat_entry_const(const unified_mat_t mat, slong i, slong j) {
    return &mat->entries[i * mat->c + j];
}

void unified_mat_swap_rows(unified_mat_t mat, slong r, slong s) {
    if (r != s) {
        void *ctx_ptr = (mat->ctx->field_id == FIELD_ID_NMOD) ? 
                       (void*)&mat->ctx->ctx.nmod_ctx : 
                       (void*)mat->ctx->ctx.fq_ctx;
        
        for (slong j = 0; j < mat->c; j++) {
            field_elem_u tmp;
            field_init_elem(&tmp, mat->ctx->field_id, ctx_ptr);
            
            /* tmp = mat[r][j] */
            field_set_elem(&tmp, unified_mat_entry(mat, r, j), mat->ctx->field_id, ctx_ptr);
            
            /* mat[r][j] = mat[s][j] */
            field_set_elem(unified_mat_entry(mat, r, j), 
                          unified_mat_entry(mat, s, j), 
                          mat->ctx->field_id, ctx_ptr);
            
            /* mat[s][j] = tmp */
            field_set_elem(unified_mat_entry(mat, s, j), &tmp, mat->ctx->field_id, ctx_ptr);
            
            field_clear_elem(&tmp, mat->ctx->field_id, ctx_ptr);
        }
    }
}

/* ============================================================================
   UNIFIED LU DETERMINANT COMPUTATION - IMPLEMENTATION
   
   This single implementation works for all finite fields through the unified
   interface, providing optimal performance for each field type:
   - Direct native operations for prime fields
   - Lookup tables for GF(2^8) and GF(2^16)
   - PCLMUL instructions for GF(2^32), GF(2^64), GF(2^128) when available
   - General FLINT operations for other finite fields
   ============================================================================ */

void unified_mat_det_lu(field_elem_u *det, const unified_mat_t mat) {
    slong n = mat->r;
    field_ctx_t *ctx = mat->ctx;
    void *ctx_ptr = (ctx->field_id == FIELD_ID_NMOD) ? 
                   (void*)&ctx->ctx.nmod_ctx : 
                   (void*)ctx->ctx.fq_ctx;
    
    if (n != mat->c) {
        flint_throw(FLINT_ERROR, "Matrix must be square for determinant calculation\n");
    }
    
    if (n == 0) {
        field_set_one(det, ctx->field_id, ctx_ptr);
        return;
    }
    
    if (n == 1) {
        field_set_elem(det, unified_mat_entry_const(mat, 0, 0), ctx->field_id, ctx_ptr);
        return;
    }
    
    /* Create working copy */
    unified_mat_t A;
    unified_mat_init(A, n, n, ctx);
    
    /* Copy matrix entries */
    for (slong i = 0; i < n; i++) {
        for (slong j = 0; j < n; j++) {
            field_set_elem(unified_mat_entry(A, i, j), 
                          unified_mat_entry_const(mat, i, j), 
                          ctx->field_id, ctx_ptr);
        }
    }
    
    /* Initialize det = 1 */
    field_set_one(det, ctx->field_id, ctx_ptr);
    
    /* Track row swaps for permutation sign */
    slong row_swaps = 0;
    
    /* LU decomposition with partial pivoting */
    for (slong k = 0; k < n - 1; k++) {
        /* Find pivot */
        slong pivot_row = k;
        const field_elem_u *pivot_val = unified_mat_entry_const(A, k, k);
        
        /* Only find other pivot when diagonal element is zero */
        if (field_is_zero(pivot_val, ctx->field_id, ctx_ptr)) {
            for (slong i = k + 1; i < n; i++) {
                const field_elem_u *val = unified_mat_entry_const(A, i, k);
                if (!field_is_zero(val, ctx->field_id, ctx_ptr)) {
                    pivot_row = i;
                    pivot_val = val;
                    break;
                }
            }
        }
        
        /* Check if matrix is singular */
        if (field_is_zero(pivot_val, ctx->field_id, ctx_ptr)) {
            field_set_zero(det, ctx->field_id, ctx_ptr);
            unified_mat_clear(A);
            return;
        }
        
        /* Swap rows if needed */
        if (pivot_row != k) {
            unified_mat_swap_rows(A, k, pivot_row);
            row_swaps++;
        }
        
        /* Get pivot and its inverse */
        field_elem_u pivot, pivot_inv;
        field_init_elem(&pivot, ctx->field_id, ctx_ptr);
        field_init_elem(&pivot_inv, ctx->field_id, ctx_ptr);
        
        field_set_elem(&pivot, unified_mat_entry(A, k, k), ctx->field_id, ctx_ptr);
        field_inv(&pivot_inv, &pivot, ctx->field_id, ctx_ptr);
        
        /* Eliminate column k below diagonal */
        for (slong i = k + 1; i < n; i++) {
            field_elem_u factor;
            field_init_elem(&factor, ctx->field_id, ctx_ptr);
            
            /* factor = A[i][k] / A[k][k] */
            field_mul(&factor, unified_mat_entry(A, i, k), &pivot_inv, ctx->field_id, ctx_ptr);
            
            /* Store L factor */
            field_set_elem(unified_mat_entry(A, i, k), &factor, ctx->field_id, ctx_ptr);
            
            /* Update row i */
            for (slong j = k + 1; j < n; j++) {
                field_elem_u prod, sum;
                field_init_elem(&prod, ctx->field_id, ctx_ptr);
                field_init_elem(&sum, ctx->field_id, ctx_ptr);
                
                /* A[i][j] = A[i][j] - factor * A[k][j] */
                field_mul(&prod, &factor, unified_mat_entry(A, k, j), ctx->field_id, ctx_ptr);
                field_sub(&sum, unified_mat_entry(A, i, j), &prod, ctx->field_id, ctx_ptr);
                field_set_elem(unified_mat_entry(A, i, j), &sum, ctx->field_id, ctx_ptr);
                
                field_clear_elem(&prod, ctx->field_id, ctx_ptr);
                field_clear_elem(&sum, ctx->field_id, ctx_ptr);
            }
            
            field_clear_elem(&factor, ctx->field_id, ctx_ptr);
        }
        
        field_clear_elem(&pivot, ctx->field_id, ctx_ptr);
        field_clear_elem(&pivot_inv, ctx->field_id, ctx_ptr);
    }
    
    /* Compute determinant as product of diagonal elements */
    field_elem_u temp;
    field_init_elem(&temp, ctx->field_id, ctx_ptr);
    
    for (slong i = 0; i < n; i++) {
        field_set_elem(&temp, det, ctx->field_id, ctx_ptr);
        field_mul(det, &temp, unified_mat_entry(A, i, i), ctx->field_id, ctx_ptr);
    }
    
    field_clear_elem(&temp, ctx->field_id, ctx_ptr);
    
    /* Handle permutation sign for non-binary fields */
    /* For fields with characteristic != 2, odd number of swaps means negate */
    if (ctx->field_id == FIELD_ID_NMOD) {
        /* Prime field case */
        if (ctx->ctx.nmod_ctx.n != 2 && (row_swaps % 2) == 1) {
            field_neg(det, det, ctx->field_id, ctx_ptr);
        }
    } else if (ctx->field_id == FIELD_ID_FQ) {
        /* General finite field case */
        if (fq_nmod_ctx_prime(ctx->ctx.fq_ctx) != 2 && (row_swaps % 2) == 1) {
            field_neg(det, det, ctx->field_id, ctx_ptr);
        }
    }
    /* For GF(2^n) fields, no negation needed since -1 = 1 */
    
    unified_mat_clear(A);
}

/* ============================================================================
   MAIN DETERMINANT FUNCTIONS - IMPLEMENTATIONS
   ============================================================================ */

/* Main determinant function - uses unified interface for ALL fields */
void fq_nmod_mat_det(fq_nmod_t det, const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx) {
    slong n = fq_nmod_mat_nrows(mat, ctx);
    
    if (n != fq_nmod_mat_ncols(mat, ctx)) {
        flint_throw(FLINT_ERROR, "Matrix must be square for determinant calculation\n");
    }
    
    if (n == 0) {
        fq_nmod_one(det, ctx);
        return;
    }
    
    if (n == 1) {
        fq_nmod_set(det, fq_nmod_mat_entry(mat, 0, 0), ctx);
        return;
    }
    
    /* Always use unified interface - it supports all field types:
     * - FIELD_ID_NMOD for prime fields
     * - FIELD_ID_GF28/16/32/64/128 for optimized binary extension fields  
     * - FIELD_ID_FQ for general finite fields
     */
    
    /* Initialize unified context */
    field_ctx_t unified_ctx;
    field_ctx_init(&unified_ctx, ctx);
    void *ctx_ptr = (unified_ctx.field_id == FIELD_ID_NMOD) ? 
                   (void*)&unified_ctx.ctx.nmod_ctx : 
                   (void*)unified_ctx.ctx.fq_ctx;
    
    /* Create unified matrix */
    unified_mat_t mat_unified;
    unified_mat_init(mat_unified, n, n, &unified_ctx);
    
    /* Convert matrix to unified format */
    for (slong i = 0; i < n; i++) {
        for (slong j = 0; j < n; j++) {
            fq_nmod_to_field_elem(unified_mat_entry(mat_unified, i, j), 
                                 fq_nmod_mat_entry(mat, i, j), 
                                 &unified_ctx);
        }
    }
    
    /* Compute determinant using unified implementation */
    field_elem_u det_unified;
    field_init_elem(&det_unified, unified_ctx.field_id, ctx_ptr);
    
    unified_mat_det_lu(&det_unified, mat_unified);
    
    /* Convert result back to fq_nmod */
    field_elem_to_fq_nmod(det, &det_unified, &unified_ctx);
    
    /* Cleanup */
    field_clear_elem(&det_unified, unified_ctx.field_id, ctx_ptr);
    unified_mat_clear(mat_unified);
}

/* Alternative: using characteristic polynomial method */
void fq_nmod_mat_det_charpoly(fq_nmod_t det, const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx) {
    slong n = fq_nmod_mat_nrows(mat, ctx);
    
    if (n != fq_nmod_mat_ncols(mat, ctx)) {
        flint_throw(FLINT_ERROR, "Matrix must be square for determinant calculation\n");
    }
    
    if (n == 0) {
        fq_nmod_one(det, ctx);
        return;
    }
    
    if (n == 1) {
        fq_nmod_set(det, fq_nmod_mat_entry(mat, 0, 0), ctx);
        return;
    }
    
    /* Compute characteristic polynomial */
    fq_nmod_poly_t charpoly;
    fq_nmod_poly_init(charpoly, ctx);
    fq_nmod_mat_charpoly(charpoly, mat, ctx);
    
    /* Determinant is the constant term times (-1)^n */
    fq_nmod_poly_get_coeff(det, charpoly, 0, ctx);
    
    if (n % 2 == 1 && fq_nmod_ctx_prime(ctx) != 2) {
        fq_nmod_neg(det, det, ctx);
    }
    
    fq_nmod_poly_clear(charpoly, ctx);
}