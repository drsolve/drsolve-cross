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
#include <flint/nmod_mat.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
    slong constant_fastpath;
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
    double constant_time;
    /* counters for the LNZ 2-branch split (Algorithm 2 of
     * Labahn-Neiger-Zhou, J. Complexity 2017) */
    slong lnz_splits;
    slong lnz_fast_nodes;
    slong lnz_full_nodes;
    slong lnz_kappa_retries;
    slong lnz_node_failures;
    slong verify_checks;
    slong verify_failures;
    double vu_kernel_time;
    double k_kernel_time;
    double completion_time;
    double stack_time;
    double verify_time;
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
    printf("    calls=%ld basecase=%ld constant_fastpath=%ld failures=%ld max_depth=%ld\n",
           p->calls, p->basecase_calls, p->constant_fastpath, p->failures, p->max_depth);
    printf("    max_dim=%ld avg_dim=%.2f avg_input_deg=%.2f max_input_deg=%ld max_schur_deg=%ld max_kernel_pivot_deg=%ld\n",
           p->max_dim,
           (p->calls > 0) ? ((double) p->total_dim / (double) p->calls) : 0.0,
           (p->calls > 0) ? ((double) p->total_degree / (double) p->calls) : 0.0,
           p->max_input_degree, p->max_schur_degree, p->max_kernel_pivot_degree);
    printf("    exact_divisions=%ld total=%.6fs recurse=%.6fs kernel=%.6fs schur_mul=%.6fs popov=%.6fs division=%.6fs constant=%.6fs\n",
           p->exact_divisions, p->total_time, p->recurse_time, p->kernel_time,
           p->schur_mul_time, p->popov_time, p->division_time, p->constant_time);
    printf("    split_copy=%.6fs degree_shift=%.6fs transpose=%.6fs pivot=%.6fs lead_extract=%.6fs schur_input_extract=%.6fs\n",
           p->split_copy_time, p->degree_shift_time, p->transpose_time,
           p->pivot_time, p->leading_extract_time, p->schur_input_extract_time);
    printf("    lnz_splits=%ld lnz_fast_nodes=%ld lnz_full_nodes=%ld lnz_kappa_retries=%ld lnz_node_failures=%ld verify_checks=%ld verify_failures=%ld\n",
           p->lnz_splits, p->lnz_fast_nodes, p->lnz_full_nodes,
           p->lnz_kappa_retries, p->lnz_node_failures,
           p->verify_checks, p->verify_failures);
    printf("    lnz_vu_kernel=%.6fs lnz_k_kernel=%.6fs lnz_completion=%.6fs lnz_stack=%.6fs verify=%.6fs\n",
           p->vu_kernel_time, p->k_kernel_time, p->completion_time,
           p->stack_time, p->verify_time);
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

/* -------------------------------------------------------------------- */
/* algorithm selection and top-level verification                        */
/* -------------------------------------------------------------------- */

#define NMOD_DET_ALGO_LNZ 0
#define NMOD_DET_ALGO_LEGACY 1

/* set to 1 while re-running with the legacy split after a failed
 * verification of the LNZ result (see nmod_poly_mat_det_hnf) */
static int g_nmod_det_force_legacy = 0;

/* DRSOLVE_PML_DET_ALGO: "lnz" (default) = 2-branch split of
 * Labahn-Neiger-Zhou 2017 (Algorithm 2, via kernel bases and the
 * Lemma 4.4 constant-determinant trick); "legacy" = previous 3-branch
 * column split with exact polynomial division. */
static int
_nmod_det_algo(void)
{
    static int cached = -1;
    if (g_nmod_det_force_legacy)
        return NMOD_DET_ALGO_LEGACY;
    if (cached < 0)
    {
        const char * env = getenv("DRSOLVE_PML_DET_ALGO");
        cached = (env != NULL && strcmp(env, "legacy") == 0)
                 ? NMOD_DET_ALGO_LEGACY : NMOD_DET_ALGO_LNZ;
    }
    return cached;
}

/* DRSOLVE_PML_DET_KAPPA: initial order factor for the ZLS kernel calls
 * (default 2.0, as in the legacy code; per-node retry with 3.0). */
static double
_nmod_det_kappa(void)
{
    static double cached = -1.0;
    if (cached < 0.0)
    {
        const char * env = getenv("DRSOLVE_PML_DET_KAPPA");
        cached = 2.0;
        if (env != NULL)
        {
            const double v = atof(env);
            if (v >= 2.0)
                cached = v;
        }
    }
    return cached;
}

/* DRSOLVE_PML_DET_VERIFY: "0" disables the final randomized-evaluation
 * check det(A)(alpha) == det(A(alpha)) (enabled by default: it costs
 * O(dim^2 deg + dim^3) base field operations, negligible compared to the
 * determinant computation itself, and it replaces the exact-division
 * self-check of the legacy algorithm). */
static int
_nmod_det_verify_enabled(void)
{
    static int cached = -1;
    if (cached < 0)
    {
        const char * env = getenv("DRSOLVE_PML_DET_VERIFY");
        cached = (env == NULL || strcmp(env, "0") != 0);
    }
    return cached;
}

/* Schwartz-Zippel style check at a pseudo-random point: a wrong result is
 * accepted with probability at most deg(det)/modulus. A mismatch proves
 * the result is wrong (det(A)(alpha) = det(A(alpha)) is an identity). */
static int
_nmod_det_verify_random_eval(const nmod_poly_t det,
                             const nmod_poly_mat_t mat)
{
    static ulong state = UWORD(0x9E3779B97F4A7C15);
    const slong dim = mat->r;
    ulong alpha, lhs, rhs;
    nmod_mat_t cmat;

    if (dim == 0)
        return nmod_poly_is_one(det);

    /* splitmix64 step: cheap, deterministic across runs, no dependency on
     * the FLINT random state API (which changed names across versions) */
    state += UWORD(0x9E3779B97F4A7C15);
    alpha = state;
    alpha = (alpha ^ (alpha >> 30)) * UWORD(0xBF58476D1CE4E5B9);
    alpha = (alpha ^ (alpha >> 27)) * UWORD(0x94D049BB133111EB);
    alpha = alpha ^ (alpha >> 31);
    alpha = alpha % mat->modulus;

    lhs = nmod_poly_evaluate_nmod(det, alpha);

    nmod_mat_init(cmat, dim, dim, mat->modulus);
    for (slong i = 0; i < dim; i++)
        for (slong j = 0; j < dim; j++)
            nmod_mat_entry(cmat, i, j) =
                nmod_poly_evaluate_nmod(nmod_poly_mat_entry(mat, i, j), alpha);
    rhs = nmod_mat_det(cmat);
    nmod_mat_clear(cmat);

    return lhs == rhs;
}

static int _nmod_poly_mat_det_hnf_recursive(nmod_poly_t det,
                                            const nmod_poly_mat_t pmat);
static int _nmod_poly_mat_det_hnf_legacy_split(nmod_poly_t det,
                                               const nmod_poly_mat_t pmat);
static int _nmod_poly_mat_det_lnz_split(nmod_poly_t det,
                                        const nmod_poly_mat_t pmat,
                                        double kappa);

static int
_nmod_poly_mat_det_hnf_recursive(nmod_poly_t det,
                                 const nmod_poly_mat_t pmat)
{
    const slong dim = pmat->r;
    const int collect_profile = _nmod_det_profile_enabled();
    const double call_start = collect_profile ? _nmod_det_now_seconds() : 0.0;
    const slong depth = g_nmod_det_hnf_depth++;
    slong input_degree;

    if (pmat->r != pmat->c)
    {
        g_nmod_det_hnf_depth--;
        return 0;
    }

    /* degree scan is O(dim^2) pointer reads: always cheap compared to the
     * rest, and it enables the constant-matrix fast path below */
    input_degree = nmod_poly_mat_degree(pmat);

    if (collect_profile)
    {
        nmod_poly_mat_det_hnf_profile_t *profile = &g_nmod_det_hnf_profile;
        profile->calls++;
        profile->total_dim += dim;
        if (dim > profile->max_dim)
            profile->max_dim = dim;
        if (depth > profile->max_depth)
            profile->max_depth = depth;
        profile->total_degree += (input_degree >= 0) ? input_degree : 0;
        if (input_degree > profile->max_input_degree)
            profile->max_input_degree = input_degree;
    }

    if (dim == 0)
    {
        nmod_poly_one(det);
        if (collect_profile)
        {
            g_nmod_det_hnf_profile.basecase_calls++;
            g_nmod_det_hnf_profile.total_time += _nmod_det_now_seconds() - call_start;
        }
        g_nmod_det_hnf_depth--;
        return 1;
    }

    /* Constant matrix (degree <= 0, including the zero matrix): the
     * determinant lives in the base field, compute it with a single dense
     * LU over Z/p instead of running the whole polynomial machinery
     * (kernel basis, weak Popov, 3-way recursion). This is a frequent case
     * deep in the recursion: e.g. in the second benchmark of time.txt the
     * average input degree over all recursive calls is 0.34, i.e. most
     * nodes are constant or nearly so. */
    if (input_degree <= 0)
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        nmod_mat_t cmat;
        nmod_mat_init(cmat, dim, dim, pmat->modulus);
        for (slong i = 0; i < dim; i++)
            for (slong j = 0; j < dim; j++)
                nmod_mat_entry(cmat, i, j) =
                    nmod_poly_get_coeff_ui(nmod_poly_mat_entry(pmat, i, j), 0);
        const ulong cdet = nmod_mat_det(cmat);
        nmod_mat_clear(cmat);
        nmod_poly_zero(det);
        nmod_poly_set_coeff_ui(det, 0, cdet);
        if (collect_profile)
        {
            g_nmod_det_hnf_profile.constant_fastpath++;
            g_nmod_det_hnf_profile.constant_time += _nmod_det_now_seconds() - t0;
            g_nmod_det_hnf_profile.total_time += _nmod_det_now_seconds() - call_start;
        }
        g_nmod_det_hnf_depth--;
        return 1;
    }

    if (dim <= NMOD_POLY_MAT_DET_HNF_BASECASE)
    {
        nmod_poly_mat_det(det, pmat);
        if (collect_profile)
        {
            g_nmod_det_hnf_profile.basecase_calls++;
            g_nmod_det_hnf_profile.total_time += _nmod_det_now_seconds() - call_start;
        }
        g_nmod_det_hnf_depth--;
        return 1;
    }

    /* ---- non-trivial dimension: dispatch on the selected split ---- */
    {
        int ok;
        if (_nmod_det_algo() == NMOD_DET_ALGO_LEGACY)
            ok = _nmod_poly_mat_det_hnf_legacy_split(det, pmat);
        else
        {
            const double kappa = _nmod_det_kappa();
            ok = _nmod_poly_mat_det_lnz_split(det, pmat, kappa);
            if (!ok && kappa < 3.0)
            {
                /* the usual ZLS failure mode is an approximant order that
                 * is too small: retry this node once with a larger order
                 * factor before giving up on the fast split */
                if (collect_profile)
                    g_nmod_det_hnf_profile.lnz_kappa_retries++;
                ok = _nmod_poly_mat_det_lnz_split(det, pmat, 3.0);
            }
            if (!ok)
            {
                /* structural fallback for this node only: the legacy
                 * 3-branch split with its exact-division self-check */
                if (collect_profile)
                    g_nmod_det_hnf_profile.lnz_node_failures++;
                ok = _nmod_poly_mat_det_hnf_legacy_split(det, pmat);
            }
        }

        if (collect_profile)
        {
            if (!ok)
                g_nmod_det_hnf_profile.failures++;
            g_nmod_det_hnf_profile.total_time += _nmod_det_now_seconds() - call_start;
        }
        g_nmod_det_hnf_depth--;
        return ok;
    }
}

/* ==================================================================== */
/* legacy 3-branch column split (previous algorithm):                   */
/*   det(A) = +- det(A_l[rows]) det(K A_r) / det(K[:,piv])              */
/* with K a left kernel basis of A_l. The third determinant and the     */
/* final exact polynomial division are what the LNZ split removes.      */
/* ==================================================================== */
static int
_nmod_poly_mat_det_hnf_legacy_split(nmod_poly_t det,
                                    const nmod_poly_mat_t pmat)
{
    const slong dim = pmat->r;
    const int collect_profile = _nmod_det_profile_enabled();
    const slong cdim1 = dim >> 1;
    const slong cdim2 = dim - cdim1;
    int ok = 0;

    /* Read-only column windows into pmat: no deep copy of the two halves. */
    nmod_poly_mat_t pmat_l;
    nmod_poly_mat_t pmat_r;
    nmod_poly_mat_window_init(pmat_l, pmat, 0, 0, dim, cdim1);
    nmod_poly_mat_window_init(pmat_r, pmat, 0, cdim1, dim, dim);
#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR == 0)
    /* windows do not carry the modulus field in old FLINT versions */
    pmat_l->modulus = pmat->modulus;
    pmat_r->modulus = pmat->modulus;
#endif

    nmod_poly_mat_t pmat_l_t;
    nmod_poly_mat_t kercols;
    nmod_poly_mat_t kerbas;
    nmod_poly_mat_t leading_block;
    nmod_poly_mat_t schur_input;
    nmod_poly_mat_t schur;

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
        is_pivot_row[i] = 0;

    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        nmod_poly_mat_row_degree(shift, pmat_l, NULL);
        for (slong i = 0; i < dim; i++)
            if (shift[i] < 0)
                shift[i] = 0;
        if (collect_profile)
            g_nmod_det_hnf_profile.degree_shift_time += _nmod_det_now_seconds() - t0;
    }

    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        nmod_poly_mat_transpose(pmat_l_t, pmat_l);
        if (collect_profile)
            g_nmod_det_hnf_profile.transpose_time += _nmod_det_now_seconds() - t0;
    }
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        slong ker_dim = nmod_poly_mat_kernel_zls(kercols, kerdeg, pmat_l_t, shift, 2.0);
        if (collect_profile)
            g_nmod_det_hnf_profile.kernel_time += _nmod_det_now_seconds() - t0;
        if (ker_dim != cdim2)
            goto cleanup;
    }

    _nmod_poly_mat_det_transpose_kernel_basis_rows(kerbas, kercols, cdim2);

    /* NOTE: this weak Popov step is NOT only pivot selection. It reduces the
     * degrees of kerbas (the ZLS output with kappa=2 is generally not fully
     * reduced), and both `schur = kerbas * pmat_r` and
     * `schur_input = kerbas[:, pivind]` inherit these smaller degrees: the
     * three recursive subproblems (and in particular the pmbasis orders
     * inside their kernel computations) get cheaper. Replacing it by a
     * constant leading-matrix pivot selection makes this step itself faster
     * but roughly doubles the pmbasis time further down -- measured overall
     * regression. Its own cost is negligible here (see popov= in the
     * profiles). Do not remove. */
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        if (nmod_poly_mat_ordered_weak_popov_iter(kerbas, shift, NULL, pivind, NULL, ROW_LOWER) != cdim2)
        {
            if (collect_profile)
                g_nmod_det_hnf_profile.popov_time += _nmod_det_now_seconds() - t0;
            goto cleanup;
        }
        if (collect_profile)
            g_nmod_det_hnf_profile.popov_time += _nmod_det_now_seconds() - t0;
    }

    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
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
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        _nmod_poly_mat_det_extract_rows(leading_block, pmat_l, comp_rows, cdim1);
        if (collect_profile)
            g_nmod_det_hnf_profile.leading_extract_time += _nmod_det_now_seconds() - t0;
    }
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        _nmod_poly_mat_det_extract_columns(schur_input, kerbas, pivind, cdim2);
        if (collect_profile)
            g_nmod_det_hnf_profile.schur_input_extract_time += _nmod_det_now_seconds() - t0;
    }
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
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
    }

    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
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
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        nmod_poly_divrem(det, rem, num, det_kernel_pivots);
        if (collect_profile)
        {
            g_nmod_det_hnf_profile.division_time += _nmod_det_now_seconds() - t0;
            g_nmod_det_hnf_profile.exact_divisions++;
        }
    }
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
    nmod_poly_mat_window_clear(pmat_r);
    nmod_poly_mat_window_clear(pmat_l);

    return ok;
}

/* ==================================================================== */
/* LNZ 2-branch split: Algorithm 2 of                                   */
/*   Labahn, Neiger, Zhou, "Fast, deterministic computation of the      */
/*   Hermite normal form and determinant of a polynomial matrix",       */
/*   Journal of Complexity 42 (2017) 44-71,                             */
/* implemented with kernel bases only.                                  */
/*                                                                      */
/* Split A = [Au; Ad] into its top m = ceil(dim/2) and bottom           */
/* k = dim - m rows. Then:                                              */
/*                                                                      */
/*  1. N  := minimal right kernel basis of Au   (dim x k, Au N = 0).    */
/*     N is a basis of a saturated module, so it can be completed to a  */
/*     unimodular U = [Ul | N]; in particular N(0) has full column      */
/*     rank.                                                            */
/*  2. Vu := minimal left kernel basis of N     (m x dim, Vu N = 0).    */
/*     Its row space is { v : v N = 0 }, which contains the rows of     */
/*     Au, so Au = B1 Vu for a unique nonsingular B1; moreover Vu is    */
/*     surjective (also a saturated basis), hence B1 is a column basis  */
/*     of Au, and Vu is the top block of U^{-1} (Lemma 4.2).            */
/*  3. B2 := Ad N                               (k x k).                */
/*     A U = [[B1, 0], [Ad Ul, B2]], so                                 */
/*        det(A) = det(B1) det(B2) / det(U)     with det(U) in K*.      */
/*  4. det(U) via Lemma 4.4: with constant matrices N0 = N mod x,       */
/*     Vu0 = Vu mod x, and any constant Ul* such that [Ul* | N0] is     */
/*     nonsingular,                                                     */
/*        1/det(U) = det(V) = det(Vu0 Ul*) / det([Ul* | N0]).           */
/*     We take Ul* = columns of the identity indexed by the complement  */
/*     of a row rank profile of N0, which makes [Ul* | N0] nonsingular  */
/*     by construction and turns Vu0 Ul* into a column selection.       */
/*  5. det(B1) without computing B1: let K = [K1 | K2] be a minimal     */
/*     left kernel basis of the stacked matrix [Vu; Au] (2m rows).      */
/*     Every kernel row (c, g) satisfies c Vu + g Au = 0, i.e.          */
/*     c = -g B1, so K1 = -K2 B1 with K2 unimodular (the g-components   */
/*     of a kernel basis form a basis of K[x]^{1 x m}). Hence           */
/*        det(B1) = (-1)^m det(K1) / det(K2(0)),                        */
/*     and we recurse on K1 (m x m) instead of B1.                      */
/*                                                                      */
/* Total: two recursive calls (K1 and B2), three kernel basis calls,    */
/* one polynomial matrix product, and a handful of constant             */
/* determinants. No polynomial division. The number of recursive        */
/* branches drops from 3 to 2, which is what makes the total cost       */
/* softly O(dim^omega * average degree) for any omega > 2.              */
/* ==================================================================== */
static int
_nmod_poly_mat_det_lnz_split(nmod_poly_t det,
                             const nmod_poly_mat_t pmat,
                             double kappa)
{
    const slong dim = pmat->r;
    const slong m = (dim + 1) / 2;     /* rows of Au (paper: ceil(n/2)) */
    const slong k = dim - m;           /* rows of Ad */
    const int collect_profile = _nmod_det_profile_enabled();
    const ulong mod = pmat->modulus;
    int ok = 0;

    nmod_t md;
    nmod_init(&md, mod);

    /* row windows */
    nmod_poly_mat_t Au, Ad;
    nmod_poly_mat_window_init(Au, pmat, 0, 0, m, dim);
    nmod_poly_mat_window_init(Ad, pmat, m, 0, dim, dim);
#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR == 0)
    Au->modulus = mod;
    Ad->modulus = mod;
#endif

    nmod_poly_mat_t N;      /* right kernel of Au, columns 0..k-1 used */
    nmod_poly_mat_t Nt;     /* its transpose (k x dim), weak Popov reduced */
    nmod_poly_mat_t Nred;   /* reduced kernel basis, back to dim x k */
    nmod_poly_mat_t VuC;    /* right kernel of Nt, columns 0..m-1 used */
    nmod_poly_mat_t Vu;     /* left kernel basis of Nred (m x dim) */
    nmod_poly_mat_t St;     /* transpose of the stack [Vu; Au] (dim x 2m) */
    nmod_poly_mat_t KC;     /* right kernel of St, columns 0..m-1 used */
    nmod_poly_mat_t K;      /* left kernel basis of [Vu; Au] (m x 2m) */
    nmod_poly_mat_t B2;     /* Ad * Nred (k x k) */

    nmod_poly_mat_init(N, dim, dim, mod);
    nmod_poly_mat_init(Nt, k, dim, mod);
    nmod_poly_mat_init(Nred, dim, k, mod);
    nmod_poly_mat_init(VuC, dim, dim, mod);
    nmod_poly_mat_init(Vu, m, dim, mod);
    nmod_poly_mat_init(St, dim, 2 * m, mod);
    nmod_poly_mat_init(KC, 2 * m, 2 * m, mod);
    nmod_poly_mat_init(K, m, 2 * m, mod);
    nmod_poly_mat_init(B2, k, k, mod);

    nmod_mat_t N0, Ustar, VuC0, K20;
    nmod_mat_init(N0, dim, k, mod);
    nmod_mat_init(Ustar, dim, dim, mod);
    nmod_mat_init(VuC0, m, m, mod);
    nmod_mat_init(K20, m, m, mod);

    slong * kerdeg = FLINT_ARRAY_ALLOC(2 * dim, slong);
    slong * shiftA = FLINT_ARRAY_ALLOC(dim, slong);
    slong * shiftN = FLINT_ARRAY_ALLOC(dim, slong);
    slong * shiftS = FLINT_ARRAY_ALLOC(2 * m, slong);
    slong * pivots = FLINT_ARRAY_ALLOC(2 * m, slong);
    slong * luperm = FLINT_ARRAY_ALLOC(dim, slong);
    slong * comp = FLINT_ARRAY_ALLOC(m, slong);
    char * is_pivot_row = FLINT_ARRAY_ALLOC(dim, char);

    nmod_poly_t detK1, detB2;
    nmod_poly_init(detK1, mod);
    nmod_poly_init(detB2, mod);

    if (collect_profile)
        g_nmod_det_hnf_profile.lnz_splits++;

    /* ---- 1) N: minimal right kernel basis of Au ---- */
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        nmod_poly_mat_column_degree(shiftA, Au, NULL);
        for (slong j = 0; j < dim; j++)
            if (shiftA[j] < 0)
                shiftA[j] = 0;
        if (collect_profile)
            g_nmod_det_hnf_profile.degree_shift_time += _nmod_det_now_seconds() - t0;
    }
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        const slong nz = nmod_poly_mat_kernel_zls(N, kerdeg, Au, shiftA, kappa);
        if (collect_profile)
            g_nmod_det_hnf_profile.kernel_time += _nmod_det_now_seconds() - t0;
        if (nz != k)
            goto cleanup;
    }
    _nmod_poly_mat_det_transpose_kernel_basis_rows(Nt, N, k);
    /* Weak Popov reduction: as in the legacy split, this reduces the
     * degrees of the (kappa=2) ZLS output, which every later step
     * inherits. Correctness only needs *a* kernel basis, but the degree
     * reduction is measurably important. */
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        const slong rk = nmod_poly_mat_ordered_weak_popov_iter(Nt, shiftA, NULL, pivots, NULL, ROW_LOWER);
        if (collect_profile)
            g_nmod_det_hnf_profile.popov_time += _nmod_det_now_seconds() - t0;
        if (rk != k)
            goto cleanup;
    }
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        nmod_poly_mat_transpose(Nred, Nt);
        if (collect_profile)
            g_nmod_det_hnf_profile.transpose_time += _nmod_det_now_seconds() - t0;
    }

    /* ---- 2) constants shared by both paths: N0, its row rank profile,
     *         the completion Ustar = [Ul* | N0], and du = det(Ustar) ---- */
    ulong du;
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        slong rank;
        nmod_mat_t lu;

        for (slong i = 0; i < dim; i++)
            for (slong j = 0; j < k; j++)
                nmod_mat_entry(N0, i, j) =
                    nmod_poly_get_coeff_ui(nmod_poly_mat_entry(Nred, i, j), 0);

        /* row rank profile of N0 via LU with row pivoting */
        nmod_mat_init_set(lu, N0);
        rank = nmod_mat_lu(luperm, lu, 0);
        nmod_mat_clear(lu);
        if (rank != k)
        {
            /* cannot happen for a genuine kernel basis (N is a basis of a
             * saturated module, hence completable to a unimodular matrix,
             * hence N(0) has full column rank); treat as a ZLS failure */
            if (collect_profile)
                g_nmod_det_hnf_profile.completion_time += _nmod_det_now_seconds() - t0;
            goto cleanup;
        }

        for (slong i = 0; i < dim; i++)
            is_pivot_row[i] = 0;
        for (slong i = 0; i < k; i++)
            is_pivot_row[luperm[i]] = 1;
        {
            slong t = 0;
            for (slong i = 0; i < dim; i++)
                if (!is_pivot_row[i])
                    comp[t++] = i;
        }

        /* Ustar = [e_{comp[0]} ... e_{comp[m-1]} | N0]: nonsingular by
         * construction (kill the unit columns: what remains is the pivot
         * submatrix N0[luperm[0..k), :], which is invertible) */
        nmod_mat_zero(Ustar);
        for (slong j = 0; j < m; j++)
            nmod_mat_entry(Ustar, comp[j], j) = 1;
        for (slong i = 0; i < dim; i++)
            for (slong j = 0; j < k; j++)
                nmod_mat_entry(Ustar, i, m + j) = nmod_mat_entry(N0, i, j);
        du = nmod_mat_det(Ustar);
        if (du == 0)
        {
            if (collect_profile)
                g_nmod_det_hnf_profile.completion_time += _nmod_det_now_seconds() - t0;
            goto cleanup;
        }
        if (collect_profile)
            g_nmod_det_hnf_profile.completion_time += _nmod_det_now_seconds() - t0;
    }

    /* ---- 3) B2 := Ad * Nred, and its determinant (needed by both
     *         paths: this is the only recursive branch of the fast path,
     *         and one of the two branches of the full path) ---- */
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        nmod_poly_mat_mul(B2, Ad, Nred);
        if (collect_profile)
            g_nmod_det_hnf_profile.schur_mul_time += _nmod_det_now_seconds() - t0;
    }
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        const int rec_ok = _nmod_poly_mat_det_hnf_recursive(detB2, B2);
        if (collect_profile)
            g_nmod_det_hnf_profile.recurse_time += _nmod_det_now_seconds() - t0;
        if (!rec_ok)
            goto cleanup;
    }

    /* ---- 4) FAST PATH (generic case): take Vu = Au itself. ----
     * If the row space of Au is saturated -- which holds generically,
     * since it fails only when the m x m minors of Au have a nontrivial
     * common factor -- then Au itself is a basis of the saturation, so it
     * is a valid choice of Vu with B1 = I. The whole determinant then
     * collapses to
     *     det(A) = dV * det(B2),   dV = det(Au(0)[:, comp]) / du,
     * with a single recursive branch and a single kernel basis call in
     * this node. Saturation is not tested directly: the resulting value
     * is checked at one random point below and on mismatch we fall through
     * to the deterministically justified full path. This is intentionally
     * permissive: even over small prime fields we still try the shortcut,
     * accepting the Monte Carlo risk that a wrong value may occasionally
     * pass the single-point check. */
    {
        const int fast_path_allowed = (dim > 0);

        if (fast_path_allowed)
        {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        ulong dnum_f;
        for (slong i = 0; i < m; i++)
            for (slong j = 0; j < m; j++)
                nmod_mat_entry(VuC0, i, j) =
                    nmod_poly_get_coeff_ui(nmod_poly_mat_entry(Au, i, comp[j]), 0);
        dnum_f = nmod_mat_det(VuC0);
        if (collect_profile)
            g_nmod_det_hnf_profile.completion_time += _nmod_det_now_seconds() - t0;

        if (dnum_f != 0)
        {
            const ulong dV_f = nmod_mul(dnum_f, nmod_inv(du, md), md);
            nmod_poly_scalar_mul_nmod(det, detB2, dV_f);

            t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
            const int pass = _nmod_det_verify_random_eval(det, pmat);
            if (collect_profile)
            {
                g_nmod_det_hnf_profile.verify_time += _nmod_det_now_seconds() - t0;
                if (pass)
                    g_nmod_det_hnf_profile.lnz_fast_nodes++;
            }
            if (pass)
            {
                ok = 1;
                goto cleanup;
            }
        }
        }
    }
    if (collect_profile)
        g_nmod_det_hnf_profile.lnz_full_nodes++;

    /* ---- 5) FULL PATH: Vu := minimal left kernel basis of Nred ----
     * (right kernel of Nt = Nred transposed, then transpose back).
     * This is the unconditional construction: the row space of Vu is
     * { v : v N = 0 }, the saturation of the row space of Au, whatever
     * the input. */
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        nmod_poly_mat_column_degree(shiftN, Nt, NULL);
        for (slong j = 0; j < dim; j++)
            if (shiftN[j] < 0)
                shiftN[j] = 0;
        if (collect_profile)
            g_nmod_det_hnf_profile.degree_shift_time += _nmod_det_now_seconds() - t0;
    }
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        const slong nz = nmod_poly_mat_kernel_zls(VuC, kerdeg, Nt, shiftN, kappa);
        if (collect_profile)
            g_nmod_det_hnf_profile.vu_kernel_time += _nmod_det_now_seconds() - t0;
        if (nz != m)
            goto cleanup;
    }
    _nmod_poly_mat_det_transpose_kernel_basis_rows(Vu, VuC, m);
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        const slong rk = nmod_poly_mat_ordered_weak_popov_iter(Vu, shiftN, NULL, pivots, NULL, ROW_LOWER);
        if (collect_profile)
            g_nmod_det_hnf_profile.popov_time += _nmod_det_now_seconds() - t0;
        if (rk != m)
            goto cleanup;
    }

    /* dV = det(Vu(0)[:, comp]) / du  (Lemma 4.4; Vu0 * Ul* is just the
     * column selection Vu(0)[:, comp]) */
    ulong dV;
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        ulong dnum;
        for (slong i = 0; i < m; i++)
            for (slong j = 0; j < m; j++)
                nmod_mat_entry(VuC0, i, j) =
                    nmod_poly_get_coeff_ui(nmod_poly_mat_entry(Vu, i, comp[j]), 0);
        dnum = nmod_mat_det(VuC0);
        if (collect_profile)
            g_nmod_det_hnf_profile.completion_time += _nmod_det_now_seconds() - t0;
        if (dnum == 0)
        {
            /* Lemma 4.3 guarantees Vu U*l unimodular, so dnum = 0
             * indicates a defective kernel basis */
            goto cleanup;
        }
        dV = nmod_mul(dnum, nmod_inv(du, md), md);
    }

    /* ---- 6) K = [K1 | K2]: minimal left kernel basis of [Vu; Au] ---- */
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        for (slong i = 0; i < m; i++)
            for (slong j = 0; j < dim; j++)
            {
                nmod_poly_set(nmod_poly_mat_entry(St, j, i),
                              nmod_poly_mat_entry(Vu, i, j));
                nmod_poly_set(nmod_poly_mat_entry(St, j, m + i),
                              nmod_poly_mat_entry(Au, i, j));
            }
        nmod_poly_mat_column_degree(shiftS, St, NULL);
        for (slong j = 0; j < 2 * m; j++)
            if (shiftS[j] < 0)
                shiftS[j] = 0;
        if (collect_profile)
            g_nmod_det_hnf_profile.stack_time += _nmod_det_now_seconds() - t0;
    }
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        const slong nz = nmod_poly_mat_kernel_zls(KC, kerdeg, St, shiftS, kappa);
        if (collect_profile)
            g_nmod_det_hnf_profile.k_kernel_time += _nmod_det_now_seconds() - t0;
        if (nz != m)
            goto cleanup;
    }
    _nmod_poly_mat_det_transpose_kernel_basis_rows(K, KC, m);
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        const slong rk = nmod_poly_mat_ordered_weak_popov_iter(K, shiftS, NULL, pivots, NULL, ROW_LOWER);
        if (collect_profile)
            g_nmod_det_hnf_profile.popov_time += _nmod_det_now_seconds() - t0;
        if (rk != m)
            goto cleanup;
    }

    /* c2 = det(K2(0)): K2 is unimodular for a genuine kernel basis, so
     * its constant term must be invertible -- a strong sanity check */
    ulong c2;
    {
        for (slong i = 0; i < m; i++)
            for (slong j = 0; j < m; j++)
                nmod_mat_entry(K20, i, j) =
                    nmod_poly_get_coeff_ui(nmod_poly_mat_entry(K, i, m + j), 0);
        c2 = nmod_mat_det(K20);
        if (c2 == 0)
            goto cleanup;
    }

    /* ---- 7) recurse on K1 (m x m) ---- */
    {
        nmod_poly_mat_t K1;
        nmod_poly_mat_window_init(K1, K, 0, 0, m, m);
#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR == 0)
        K1->modulus = mod;
#endif
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        const int rec_ok = _nmod_poly_mat_det_hnf_recursive(detK1, K1);
        if (collect_profile)
            g_nmod_det_hnf_profile.recurse_time += _nmod_det_now_seconds() - t0;
        nmod_poly_mat_window_clear(K1);
        if (!rec_ok)
            goto cleanup;
    }

    /* ---- 8) det(A) = dV * (-1)^m * det(K1) det(B2) / c2 ---- */
    {
        ulong scalar = nmod_mul(dV, nmod_inv(c2, md), md);
        if (m & 1)
            scalar = nmod_neg(scalar, md);
        nmod_poly_mul(det, detK1, detB2);
        nmod_poly_scalar_mul_nmod(det, det, scalar);
    }
    /* same random-point check as the fast path: catches a defective
     * kernel basis that slipped past the structural checks above, and
     * lets the caller retry with a larger kappa or the legacy split */
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        ok = _nmod_det_verify_random_eval(det, pmat);
        if (collect_profile)
            g_nmod_det_hnf_profile.verify_time += _nmod_det_now_seconds() - t0;
    }

cleanup:
    nmod_poly_clear(detB2);
    nmod_poly_clear(detK1);

    flint_free(is_pivot_row);
    flint_free(comp);
    flint_free(luperm);
    flint_free(pivots);
    flint_free(shiftS);
    flint_free(shiftN);
    flint_free(shiftA);
    flint_free(kerdeg);

    nmod_mat_clear(K20);
    nmod_mat_clear(VuC0);
    nmod_mat_clear(Ustar);
    nmod_mat_clear(N0);

    nmod_poly_mat_clear(B2);
    nmod_poly_mat_clear(K);
    nmod_poly_mat_clear(KC);
    nmod_poly_mat_clear(St);
    nmod_poly_mat_clear(Vu);
    nmod_poly_mat_clear(VuC);
    nmod_poly_mat_clear(Nred);
    nmod_poly_mat_clear(Nt);
    nmod_poly_mat_clear(N);

    nmod_poly_mat_window_clear(Ad);
    nmod_poly_mat_window_clear(Au);

    return ok;
}

// Mulders-Storjohann determinant algorithm
// matrix must be square, not tested
// this modifies mat, if all goes well (input matrix nonsingular) at the end it
// is upper triangular (up to permutation)
void nmod_poly_mat_det_iter(nmod_poly_t det, nmod_poly_mat_t mat)
{
    const int collect_profile = _nmod_det_profile_enabled();
    const double call_start = collect_profile ? _nmod_det_now_seconds() : 0.0;
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
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
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
        t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
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
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
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
    {
        g_nmod_det_iter_profile.total_time += _nmod_det_now_seconds() - call_start;
        if (depth == 0)
            _nmod_det_iter_profile_print();
    }
    g_nmod_det_iter_depth--;
}

/* Run the recursion, then (unless disabled) verify the result at a random
 * point. If the LNZ algorithm produced a value that fails verification --
 * which the per-node sanity checks should already prevent -- recompute
 * once with the legacy 3-branch algorithm, whose exact-division check is
 * a different failure detector, and verify again. */
static int
_nmod_poly_mat_det_hnf_toplevel(nmod_poly_t det,
                                const nmod_poly_mat_t mat)
{
    const int collect_profile = _nmod_det_profile_enabled();
    if (collect_profile) {
        _nmod_det_hnf_profile_reset();
        nmod_poly_mat_kernel_zls_profile_reset();
    }

    int ok = _nmod_poly_mat_det_hnf_recursive(det, mat);

    if (ok && _nmod_det_verify_enabled())
    {
        double t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
        int pass = _nmod_det_verify_random_eval(det, mat);
        if (collect_profile)
        {
            g_nmod_det_hnf_profile.verify_checks++;
            if (!pass)
                g_nmod_det_hnf_profile.verify_failures++;
            g_nmod_det_hnf_profile.verify_time += _nmod_det_now_seconds() - t0;
        }
        if (!pass && _nmod_det_algo() != NMOD_DET_ALGO_LEGACY)
        {
            g_nmod_det_force_legacy = 1;
            ok = _nmod_poly_mat_det_hnf_recursive(det, mat);
            g_nmod_det_force_legacy = 0;
            if (ok)
            {
                t0 = collect_profile ? _nmod_det_now_seconds() : 0.0;
                pass = _nmod_det_verify_random_eval(det, mat);
                if (collect_profile)
                {
                    g_nmod_det_hnf_profile.verify_checks++;
                    if (!pass)
                        g_nmod_det_hnf_profile.verify_failures++;
                    g_nmod_det_hnf_profile.verify_time += _nmod_det_now_seconds() - t0;
                }
            }
        }
        ok = ok && pass;
    }

    if (collect_profile)
        _nmod_det_hnf_profile_print();

    return ok;
}

int nmod_poly_mat_det_hnf(nmod_poly_t det,
                          const nmod_poly_mat_t mat)
{
    return _nmod_poly_mat_det_hnf_toplevel(det, mat);
}

int nmod_poly_mat_det_hnf_knowing_degree(nmod_poly_t det,
                                         const nmod_poly_mat_t mat,
                                         slong degree)
{
    const int ok = _nmod_poly_mat_det_hnf_toplevel(det, mat);
    return ok && (degree == nmod_poly_degree(det));
}
