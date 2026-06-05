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
#include <time.h>
#include <stdio.h>
#include "nmod_poly_mat_utils.h" // for permute_rows_by_sorting_vec
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_kernel.h"

#ifndef NMOD_POLY_MAT_DET_HNF_BASECASE
#define NMOD_POLY_MAT_DET_HNF_BASECASE 8
#endif

extern int g_dixon_verbose_level;
extern int g_dixon_debug_mode;

typedef struct
{
    slong calls;
    slong basecase_calls;
    slong failures;
    slong max_depth;
    slong max_dim;
    slong total_dim;
    slong total_degree;
    slong max_input_degree;
    slong max_schur_degree;
    slong max_kernel_pivot_degree;
    slong exact_divisions;
    double total_time;
    double split_copy_time;
    double degree_shift_time;
    double transpose_time;
    double kernel_time;
    double popov_time;
    double pivot_time;
    double leading_extract_time;
    double schur_input_extract_time;
    double schur_mul_time;
    double recurse_time;
    double division_time;
} nmod_poly_mat_det_hnf_profile_t;

typedef struct
{
    slong calls;
    slong total_dim;
    slong max_dim;
    slong loop_iterations;
    slong rank_defects;
    double total_time;
    double weak_popov_time;
    double permute_time;
    double diagonal_product_time;
} nmod_poly_mat_det_iter_profile_t;

static nmod_poly_mat_det_hnf_profile_t g_nmod_det_hnf_profile;
static nmod_poly_mat_det_iter_profile_t g_nmod_det_iter_profile;
static slong g_nmod_det_hnf_depth = 0;
static slong g_nmod_det_iter_depth = 0;

static double
_nmod_det_now_seconds(void)
{
    return ((double) clock()) / CLOCKS_PER_SEC;
}

static int
_nmod_det_profile_enabled(void)
{
    return g_dixon_verbose_level >= 3 && g_dixon_debug_mode;
}

static int
_nmod_det_profile_heavy_enabled(void)
{
    return g_dixon_verbose_level >= 3 && g_dixon_debug_mode;
}

static void
_nmod_det_hnf_profile_reset(void)
{
    g_nmod_det_hnf_profile = (nmod_poly_mat_det_hnf_profile_t) {0};
}

static void
_nmod_det_iter_profile_reset(void)
{
    g_nmod_det_iter_profile = (nmod_poly_mat_det_iter_profile_t) {0};
}

static void
_nmod_det_hnf_profile_print(void)
{
    if (!_nmod_det_profile_enabled())
        return;

    const nmod_poly_mat_det_hnf_profile_t *p = &g_nmod_det_hnf_profile;
    printf("\n  [PML HNF det profile]\n");
    printf("    calls=%ld basecase=%ld failures=%ld max_depth=%ld\n",
           p->calls, p->basecase_calls, p->failures, p->max_depth);
    printf("    max_dim=%ld avg_dim=%.2f avg_input_deg=%.2f max_input_deg=%ld max_schur_deg=%ld max_kernel_pivot_deg=%ld\n",
           p->max_dim,
           (p->calls > 0) ? ((double) p->total_dim / (double) p->calls) : 0.0,
           (p->calls > 0) ? ((double) p->total_degree / (double) p->calls) : 0.0,
           p->max_input_degree, p->max_schur_degree, p->max_kernel_pivot_degree);
    printf("    exact_divisions=%ld total=%.6fs recurse=%.6fs kernel=%.6fs schur_mul=%.6fs popov=%.6fs division=%.6fs\n",
           p->exact_divisions, p->total_time, p->recurse_time, p->kernel_time,
           p->schur_mul_time, p->popov_time, p->division_time);
    printf("    split_copy=%.6fs degree_shift=%.6fs transpose=%.6fs pivot=%.6fs lead_extract=%.6fs schur_input_extract=%.6fs\n",
           p->split_copy_time, p->degree_shift_time, p->transpose_time,
           p->pivot_time, p->leading_extract_time, p->schur_input_extract_time);
    nmod_poly_mat_kernel_zls_profile_print();
}

static void
_nmod_det_iter_profile_print(void)
{
    if (!_nmod_det_profile_enabled())
        return;

    const nmod_poly_mat_det_iter_profile_t *p = &g_nmod_det_iter_profile;
    printf("\n  [PML iter det profile]\n");
    printf("    calls=%ld loop_iterations=%ld rank_defects=%ld max_dim=%ld avg_dim=%.2f\n",
           p->calls, p->loop_iterations, p->rank_defects, p->max_dim,
           (p->calls > 0) ? ((double) p->total_dim / (double) p->calls) : 0.0);
    printf("    total=%.6fs weak_popov=%.6fs permute=%.6fs diagonal_product=%.6fs\n",
           p->total_time, p->weak_popov_time, p->permute_time, p->diagonal_product_time);
}

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
_nmod_poly_mat_det_hnf_recursive(nmod_poly_t det,
                                 const nmod_poly_mat_t pmat)
{
    const slong dim = pmat->r;
    const double call_start = _nmod_det_now_seconds();
    const int collect_profile = _nmod_det_profile_enabled();
    const int heavy_profile = _nmod_det_profile_heavy_enabled();
    const slong depth = g_nmod_det_hnf_depth++;
    slong input_degree = -1;

    if (pmat->r != pmat->c)
    {
        if (collect_profile && depth == 0)
            _nmod_det_hnf_profile_print();
        g_nmod_det_hnf_depth--;
        return 0;
    }

    if (collect_profile)
    {
        nmod_poly_mat_det_hnf_profile_t *profile = &g_nmod_det_hnf_profile;
        profile->calls++;
        profile->total_dim += dim;
        if (dim > profile->max_dim)
            profile->max_dim = dim;
        if (depth > profile->max_depth)
            profile->max_depth = depth;
        input_degree = nmod_poly_mat_degree(pmat);
        profile->total_degree += (input_degree >= 0) ? input_degree : 0;
        if (input_degree > profile->max_input_degree)
            profile->max_input_degree = input_degree;
    }

    if (dim == 0)
    {
        nmod_poly_one(det);
        if (collect_profile)
            g_nmod_det_hnf_profile.basecase_calls++;
        g_nmod_det_hnf_profile.total_time += _nmod_det_now_seconds() - call_start;
        if (collect_profile && depth == 0)
            _nmod_det_hnf_profile_print();
        g_nmod_det_hnf_depth--;
        return 1;
    }

    if (dim <= NMOD_POLY_MAT_DET_HNF_BASECASE)
    {
        nmod_poly_mat_det(det, pmat);
        if (collect_profile)
            g_nmod_det_hnf_profile.basecase_calls++;
        g_nmod_det_hnf_profile.total_time += _nmod_det_now_seconds() - call_start;
        if (collect_profile && heavy_profile && 0)
            printf("    [PML HNF node] depth=%ld dim=%ld input_deg=%ld basecase=1\n",
                   depth, dim, input_degree);
        if (collect_profile && depth == 0)
            _nmod_det_hnf_profile_print();
        g_nmod_det_hnf_depth--;
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
    if (collect_profile)
        g_nmod_det_hnf_profile.split_copy_time += _nmod_det_now_seconds() - call_start;

    {
        double t0 = _nmod_det_now_seconds();
    nmod_poly_mat_row_degree(shift, pmat_l, NULL);
    for (slong i = 0; i < dim; i++)
        if (shift[i] < 0)
            shift[i] = 0;
        if (collect_profile)
            g_nmod_det_hnf_profile.degree_shift_time += _nmod_det_now_seconds() - t0;
    }

    {
        double t0 = _nmod_det_now_seconds();
    nmod_poly_mat_transpose(pmat_l_t, pmat_l);
        if (collect_profile)
            g_nmod_det_hnf_profile.transpose_time += _nmod_det_now_seconds() - t0;
    }
    {
        double t0 = _nmod_det_now_seconds();
    slong ker_dim = nmod_poly_mat_kernel_zls(kercols, kerdeg, pmat_l_t, shift, 2.0);
        if (collect_profile)
            g_nmod_det_hnf_profile.kernel_time += _nmod_det_now_seconds() - t0;
    if (ker_dim != cdim2)
        goto cleanup;
    }

    _nmod_poly_mat_det_transpose_kernel_basis_rows(kerbas, kercols, cdim2);
    {
        double t0 = _nmod_det_now_seconds();
    if (nmod_poly_mat_ordered_weak_popov_iter(kerbas, shift, NULL, pivind, NULL, ROW_LOWER) != cdim2)
        goto cleanup;
        if (collect_profile)
            g_nmod_det_hnf_profile.popov_time += _nmod_det_now_seconds() - t0;
    }

    {
        double t0 = _nmod_det_now_seconds();
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
        if (collect_profile)
            g_nmod_det_hnf_profile.pivot_time += _nmod_det_now_seconds() - t0;
    }

    {
        double t0 = _nmod_det_now_seconds();
    _nmod_poly_mat_det_extract_rows(leading_block, pmat_l, comp_rows, cdim1);
        if (collect_profile)
            g_nmod_det_hnf_profile.leading_extract_time += _nmod_det_now_seconds() - t0;
    }
    {
        double t0 = _nmod_det_now_seconds();
    _nmod_poly_mat_det_extract_columns(schur_input, kerbas, pivind, cdim2);
        if (collect_profile)
            g_nmod_det_hnf_profile.schur_input_extract_time += _nmod_det_now_seconds() - t0;
    }
    {
        double t0 = _nmod_det_now_seconds();
    nmod_poly_mat_mul(schur, kerbas, pmat_r);
        if (collect_profile)
            g_nmod_det_hnf_profile.schur_mul_time += _nmod_det_now_seconds() - t0;
    }
    if (collect_profile)
    {
        slong schur_deg = nmod_poly_mat_degree(schur);
        slong schur_input_deg = nmod_poly_mat_degree(schur_input);
        if (schur_deg > g_nmod_det_hnf_profile.max_schur_degree)
            g_nmod_det_hnf_profile.max_schur_degree = schur_deg;
        if (schur_input_deg > g_nmod_det_hnf_profile.max_kernel_pivot_degree)
            g_nmod_det_hnf_profile.max_kernel_pivot_degree = schur_input_deg;
        if (heavy_profile && 0)
            printf("    [PML HNF node] depth=%ld dim=%ld split=(%ld,%ld) input_deg=%ld schur_deg=%ld ker_pivot_deg=%ld\n",
                   depth, dim, cdim1, cdim2, input_degree, schur_deg, schur_input_deg);
    }

    {
        double t0 = _nmod_det_now_seconds();
    if (!_nmod_poly_mat_det_hnf_recursive(det_leading, leading_block))
        goto cleanup;
    if (!_nmod_poly_mat_det_hnf_recursive(det_schur, schur))
        goto cleanup;
    if (!_nmod_poly_mat_det_hnf_recursive(det_kernel_pivots, schur_input))
        goto cleanup;
        if (collect_profile)
            g_nmod_det_hnf_profile.recurse_time += _nmod_det_now_seconds() - t0;
    }
    if (nmod_poly_is_zero(det_kernel_pivots))
        goto cleanup;

    nmod_poly_mul(num, det_leading, det_schur);
    if (_nmod_poly_mat_det_permutation_sign_from_block_order(comp_rows, cdim1, pivind, cdim2) < 0)
        nmod_poly_neg(num, num);

    {
        double t0 = _nmod_det_now_seconds();
    nmod_poly_divrem(det, rem, num, det_kernel_pivots);
        if (collect_profile)
        {
            g_nmod_det_hnf_profile.division_time += _nmod_det_now_seconds() - t0;
            g_nmod_det_hnf_profile.exact_divisions++;
        }
    }
    ok = nmod_poly_is_zero(rem);

cleanup:
    if (collect_profile && !ok)
        g_nmod_det_hnf_profile.failures++;
    if (collect_profile)
        g_nmod_det_hnf_profile.total_time += _nmod_det_now_seconds() - call_start;

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

    if (collect_profile && depth == 0)
        _nmod_det_hnf_profile_print();
    g_nmod_det_hnf_depth--;
    return ok;
}

// Mulders-Storjohann determinant algorithm
// matrix must be square, not tested
// this modifies mat, if all goes well (input matrix nonsingular) at the end it
// is upper triangular (up to permutation)
void nmod_poly_mat_det_iter(nmod_poly_t det, nmod_poly_mat_t mat)
{
    const double call_start = _nmod_det_now_seconds();
    const int collect_profile = _nmod_det_profile_enabled();
    const slong depth = g_nmod_det_iter_depth++;

    if (mat->r == 0) {
        nmod_poly_one(det);
        if (collect_profile) {
            g_nmod_det_iter_profile.calls++;
            g_nmod_det_iter_profile.total_time += _nmod_det_now_seconds() - call_start;
            if (depth == 0)
                _nmod_det_iter_profile_print();
        }
        g_nmod_det_iter_depth--;
        return;
    }

    if (collect_profile)
    {
        g_nmod_det_iter_profile.calls++;
        g_nmod_det_iter_profile.total_dim += mat->r;
        if (mat->r > g_nmod_det_iter_profile.max_dim)
            g_nmod_det_iter_profile.max_dim = mat->r;
    }

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
        double t0 = _nmod_det_now_seconds();
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
        if (collect_profile)
        {
            g_nmod_det_iter_profile.loop_iterations++;
            g_nmod_det_iter_profile.weak_popov_time += _nmod_det_now_seconds() - t0;
        }

        // early exit if rank-deficient
        if (rk < i || nmod_poly_is_zero(nmod_poly_mat_entry(view, i, i)))
        {
            nmod_poly_zero(det);
            if (collect_profile)
                g_nmod_det_iter_profile.rank_defects++;
            goto cleanup;
        }

        // permute into ordered weak Popov form
        t0 = _nmod_det_now_seconds();
        _nmod_poly_mat_permute_rows_by_sorting_vec(view, rk, pivind, perm);
        _nmod_poly_mat_window_resize_columns(view, -1);
        if (_perm_parity(perm, rk)) // odd permutation, negate udet
            udet = -udet;
        if (collect_profile)
            g_nmod_det_iter_profile.permute_time += _nmod_det_now_seconds() - t0;
    }
    // retrieve determinant as product of diagonal entries
    // use view rather than mat, since mat has not been row-permuted
    // so it is only triangular up to permutation
    _nmod_poly_mat_window_resize_columns(view, mat->r -1); // reset view->c to mat->c
    nmod_poly_set(det, nmod_poly_mat_entry(view, 0, 0)); // recall here mat->r == mat->c > 0
    if (nmod_poly_is_zero(det))
        goto cleanup; // rank deficient early exit, [0,0] had not been tested yet
    if (udet == -1)
        _nmod_vec_neg(det->coeffs, det->coeffs, det->length, det->mod);
    {
        double t0 = _nmod_det_now_seconds();
    for (slong i = 1; i < view->r; i++)
        nmod_poly_mul(det, det, nmod_poly_mat_entry(view, i, i));
        if (collect_profile)
            g_nmod_det_iter_profile.diagonal_product_time += _nmod_det_now_seconds() - t0;
    }

cleanup:
    nmod_poly_mat_window_clear(view);
    flint_free(pivind);
    flint_free(perm);
    if (collect_profile)
        g_nmod_det_iter_profile.total_time += _nmod_det_now_seconds() - call_start;
    if (collect_profile && depth == 0)
        _nmod_det_iter_profile_print();
    g_nmod_det_iter_depth--;
}

int nmod_poly_mat_det_hnf(nmod_poly_t det,
                          const nmod_poly_mat_t mat)
{
    if (_nmod_det_profile_enabled()) {
        _nmod_det_hnf_profile_reset();
        nmod_poly_mat_kernel_zls_profile_reset();
    }
    return _nmod_poly_mat_det_hnf_recursive(det, mat);
}

int nmod_poly_mat_det_hnf_knowing_degree(nmod_poly_t det,
                                         const nmod_poly_mat_t mat,
                                         slong degree)
{
    if (_nmod_det_profile_enabled()) {
        _nmod_det_hnf_profile_reset();
        nmod_poly_mat_kernel_zls_profile_reset();
    }
    const int ok = _nmod_poly_mat_det_hnf_recursive(det, mat);
    return ok && (degree == nmod_poly_degree(det));
}
