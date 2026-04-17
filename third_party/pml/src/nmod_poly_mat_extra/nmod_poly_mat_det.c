/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_vec.h>
#include "nmod_poly_mat_utils.h" // for permute_rows_by_sorting_vec
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_kernel.h"

#ifndef NMOD_POLY_MAT_DET_GENERIC_BASECASE
#define NMOD_POLY_MAT_DET_GENERIC_BASECASE 8
#endif

static slong
_nmod_poly_mat_det_permutation_sign_from_block_order(const slong * first,
                                                     slong first_len,
                                                     const slong * second,
                                                     slong second_len)
{
    slong inv_parity = 0;

    for (slong i = 0; i < first_len; i++)
        for (slong j = 0; j < second_len; j++)
            inv_parity ^= (second[j] < first[i]);

    return inv_parity ? -1 : 1;
}

static void
_nmod_poly_mat_det_extract_rows(nmod_poly_mat_t sub,
                                const nmod_poly_mat_t mat,
                                const slong * rows,
                                slong nb_rows)
{
    for (slong i = 0; i < nb_rows; i++)
        for (slong j = 0; j < mat->c; j++)
            nmod_poly_set(nmod_poly_mat_entry(sub, i, j),
                          nmod_poly_mat_entry(mat, rows[i], j));
}

static void
_nmod_poly_mat_det_extract_columns(nmod_poly_mat_t sub,
                                   const nmod_poly_mat_t mat,
                                   const slong * cols,
                                   slong nb_cols)
{
    for (slong i = 0; i < mat->r; i++)
        for (slong j = 0; j < nb_cols; j++)
            nmod_poly_set(nmod_poly_mat_entry(sub, i, j),
                          nmod_poly_mat_entry(mat, i, cols[j]));
}

static void
_nmod_poly_mat_det_transpose_kernel_basis_rows(nmod_poly_mat_t kerbas,
                                               const nmod_poly_mat_t kercols,
                                               slong ker_dim)
{
    for (slong i = 0; i < ker_dim; i++)
        for (slong j = 0; j < kercols->r; j++)
            nmod_poly_set(nmod_poly_mat_entry(kerbas, i, j),
                          nmod_poly_mat_entry(kercols, j, i));
}

static int
_nmod_poly_mat_det_generic_recursive(nmod_poly_t det,
                                     const nmod_poly_mat_t pmat)
{
    const slong dim = pmat->r;

    if (pmat->r != pmat->c)
        return 0;

    if (dim == 0)
    {
        nmod_poly_one(det);
        return 1;
    }

    if (dim <= NMOD_POLY_MAT_DET_GENERIC_BASECASE)
    {
        nmod_poly_mat_det(det, pmat);
        return 1;
    }

    const slong cdim1 = dim >> 1;
    const slong cdim2 = dim - cdim1;
    int ok = 0;

    nmod_poly_mat_t pmat_l;
    nmod_poly_mat_t pmat_r;
    nmod_poly_mat_t pmat_l_t;
    nmod_poly_mat_t kercols;
    nmod_poly_mat_t kerbas;
    nmod_poly_mat_t leading_block;
    nmod_poly_mat_t schur_input;
    nmod_poly_mat_t schur;

    nmod_poly_mat_init(pmat_l, dim, cdim1, pmat->modulus);
    nmod_poly_mat_init(pmat_r, dim, cdim2, pmat->modulus);
    nmod_poly_mat_init(pmat_l_t, cdim1, dim, pmat->modulus);
    nmod_poly_mat_init(kercols, dim, dim, pmat->modulus);
    nmod_poly_mat_init(kerbas, cdim2, dim, pmat->modulus);
    nmod_poly_mat_init(leading_block, cdim1, cdim1, pmat->modulus);
    nmod_poly_mat_init(schur_input, cdim2, cdim2, pmat->modulus);
    nmod_poly_mat_init(schur, cdim2, cdim2, pmat->modulus);

    slong * kerdeg = FLINT_ARRAY_ALLOC(dim, slong);
    slong * pivind = FLINT_ARRAY_ALLOC(cdim2, slong);
    slong * comp_rows = FLINT_ARRAY_ALLOC(cdim1, slong);
    slong * shift = FLINT_ARRAY_ALLOC(dim, slong);
    char * is_pivot_row = FLINT_ARRAY_ALLOC(dim, char);

    nmod_poly_t det_leading, det_schur, det_kernel_pivots, num, rem;
    nmod_poly_init(det_leading, pmat->modulus);
    nmod_poly_init(det_schur, pmat->modulus);
    nmod_poly_init(det_kernel_pivots, pmat->modulus);
    nmod_poly_init(num, pmat->modulus);
    nmod_poly_init(rem, pmat->modulus);

    for (slong i = 0; i < dim; i++)
    {
        for (slong j = 0; j < cdim1; j++)
            nmod_poly_set(nmod_poly_mat_entry(pmat_l, i, j),
                          nmod_poly_mat_entry(pmat, i, j));
        for (slong j = 0; j < cdim2; j++)
            nmod_poly_set(nmod_poly_mat_entry(pmat_r, i, j),
                          nmod_poly_mat_entry(pmat, i, j + cdim1));
        is_pivot_row[i] = 0;
    }

    nmod_poly_mat_row_degree(shift, pmat_l, NULL);
    for (slong i = 0; i < dim; i++)
        if (shift[i] < 0)
            shift[i] = 0;

    nmod_poly_mat_transpose(pmat_l_t, pmat_l);
    slong ker_dim = nmod_poly_mat_kernel_zls(kercols, kerdeg, pmat_l_t, shift, 2.0);
    if (ker_dim != cdim2)
        goto cleanup;

    _nmod_poly_mat_det_transpose_kernel_basis_rows(kerbas, kercols, cdim2);
    if (nmod_poly_mat_ordered_weak_popov_iter(kerbas, shift, NULL, pivind, NULL, ROW_LOWER) != cdim2)
        goto cleanup;

    for (slong i = 0; i < cdim2; i++)
    {
        if (pivind[i] < 0 || pivind[i] >= dim || is_pivot_row[pivind[i]])
            goto cleanup;
        is_pivot_row[pivind[i]] = 1;
    }

    {
        slong idx = 0;
        for (slong i = 0; i < dim; i++)
            if (!is_pivot_row[i])
                comp_rows[idx++] = i;
        if (idx != cdim1)
            goto cleanup;
    }

    _nmod_poly_mat_det_extract_rows(leading_block, pmat_l, comp_rows, cdim1);
    _nmod_poly_mat_det_extract_columns(schur_input, kerbas, pivind, cdim2);
    nmod_poly_mat_mul(schur, kerbas, pmat_r);

    if (!_nmod_poly_mat_det_generic_recursive(det_leading, leading_block))
        goto cleanup;
    if (!_nmod_poly_mat_det_generic_recursive(det_schur, schur))
        goto cleanup;
    if (!_nmod_poly_mat_det_generic_recursive(det_kernel_pivots, schur_input))
        goto cleanup;
    if (nmod_poly_is_zero(det_kernel_pivots))
        goto cleanup;

    nmod_poly_mul(num, det_leading, det_schur);
    if (_nmod_poly_mat_det_permutation_sign_from_block_order(comp_rows, cdim1, pivind, cdim2) < 0)
        nmod_poly_neg(num, num);

    nmod_poly_divrem(det, rem, num, det_kernel_pivots);
    ok = nmod_poly_is_zero(rem);

cleanup:
    nmod_poly_clear(rem);
    nmod_poly_clear(num);
    nmod_poly_clear(det_kernel_pivots);
    nmod_poly_clear(det_schur);
    nmod_poly_clear(det_leading);

    flint_free(is_pivot_row);
    flint_free(shift);
    flint_free(comp_rows);
    flint_free(pivind);
    flint_free(kerdeg);

    nmod_poly_mat_clear(schur);
    nmod_poly_mat_clear(schur_input);
    nmod_poly_mat_clear(leading_block);
    nmod_poly_mat_clear(kerbas);
    nmod_poly_mat_clear(kercols);
    nmod_poly_mat_clear(pmat_l_t);
    nmod_poly_mat_clear(pmat_r);
    nmod_poly_mat_clear(pmat_l);

    return ok;
}

// Mulders-Storjohann determinant algorithm
// matrix must be square, not tested
// this modifies mat, if all goes well (input matrix nonsingular) at the end it
// is upper triangular (up to permutation)
void nmod_poly_mat_det_iter(nmod_poly_t det, nmod_poly_mat_t mat)
{
    if (mat->r == 0) { nmod_poly_one(det); return; }

    // determinant (+1 or -1), pivot index, and rank for weak Popov form computations
    slong udet = 1; 
    slong rk;
    slong * pivind = flint_malloc(mat->r * sizeof(slong));
    // permutation for putting into ordered weak Popov
    slong * perm = flint_malloc(mat->r * sizeof(slong));

    // window of original mat
    nmod_poly_mat_t view;
    nmod_poly_mat_window_init(view, mat, 0, 0, mat->r, mat->r);
    // TODO the following is necessary for Flint <= v3.0.0
    //view->modulus = mat->modulus; 

    for (slong i = mat->r -1; i >= 1; i--)
    {
        //            [ V  * ]                 [ V ]
        // mat is now [ 0  D ] and view is now [ 0 ] of size mat->r x (i+1),
        // with det = det(D) and V of size (i+1) x (i+1)
        //                 [ V'  *  * ]
        // -> transform to [ 0   d  * ] via weak Popov form of V[:,:i]
        //                 [ 0   0  D ]

        // this applies weak Popov form on view[:i+1,:i],
        // with transformations applied to the whole rows view[:i+1,:] = V
        // with early exit if detecting rank < i (in which case rk < 0)
        // with update of the determinant of unimodular transformation (+1 or -1)
        rk = _nmod_poly_mat_weak_popov_iter_submat_rowbyrow(view, NULL, NULL, &udet, pivind, NULL, 0, 0, i+1, i, 2, ROW_UPPER);

        // early exit if rank-deficient
        if (rk < i || nmod_poly_is_zero(nmod_poly_mat_entry(view, i, i)))
        {
            nmod_poly_zero(det);
            return;
        }

        // permute into ordered weak Popov form
        _nmod_poly_mat_permute_rows_by_sorting_vec(view, rk, pivind, perm);
        _nmod_poly_mat_window_resize_columns(view, -1);
        if (_perm_parity(perm, rk)) // odd permutation, negate udet
            udet = -udet;
    }
    flint_free(pivind);
    flint_free(perm);

    // retrieve determinant as product of diagonal entries
    // use view rather than mat, since mat has not been row-permuted
    // so it is only triangular up to permutation
    _nmod_poly_mat_window_resize_columns(view, mat->r -1); // reset view->c to mat->c
    nmod_poly_set(det, nmod_poly_mat_entry(view, 0, 0)); // recall here mat->r == mat->c > 0
    if (nmod_poly_is_zero(det))
        return; // rank deficient early exit, [0,0] had not been tested yet
    if (udet == -1)
        _nmod_vec_neg(det->coeffs, det->coeffs, det->length, det->mod);
    for (slong i = 1; i < view->r; i++)
        nmod_poly_mul(det, det, nmod_poly_mat_entry(view, i, i));
    nmod_poly_mat_window_clear(view);
}

int nmod_poly_mat_det_generic(nmod_poly_t det,
                              const nmod_poly_mat_t mat)
{
    return _nmod_poly_mat_det_generic_recursive(det, mat);
}

int nmod_poly_mat_det_generic_knowing_degree(nmod_poly_t det,
                                             const nmod_poly_mat_t mat,
                                             slong degree)
{
    const int ok = _nmod_poly_mat_det_generic_recursive(det, mat);
    return ok && (degree == nmod_poly_degree(det));
}
