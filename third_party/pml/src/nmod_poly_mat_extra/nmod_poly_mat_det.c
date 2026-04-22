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
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "nmod_poly_mat_utils.h" // for permute_rows_by_sorting_vec
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_kernel.h"

#ifndef NMOD_POLY_MAT_DET_GENERIC_BASECASE
#define NMOD_POLY_MAT_DET_GENERIC_BASECASE 8
#endif

#ifndef NMOD_POLY_MAT_DET_GENERIC_PARALLEL_THRESHOLD
#define NMOD_POLY_MAT_DET_GENERIC_PARALLEL_THRESHOLD 16
#endif

#ifndef NMOD_POLY_MAT_DET_GENERIC_MAX_PARALLEL_DEPTH
#define NMOD_POLY_MAT_DET_GENERIC_MAX_PARALLEL_DEPTH 1
#endif

enum
{
    NMOD_DET_PROFILE_COPY_SPLIT = 0,
    NMOD_DET_PROFILE_SHIFT,
    NMOD_DET_PROFILE_TRANSPOSE,
    NMOD_DET_PROFILE_KERNEL,
    NMOD_DET_PROFILE_WEAK_POPOV,
    NMOD_DET_PROFILE_PIVOT_SELECT,
    NMOD_DET_PROFILE_EXTRACT,
    NMOD_DET_PROFILE_SCHUR_MUL,
    NMOD_DET_PROFILE_RECURSE_LEADING,
    NMOD_DET_PROFILE_RECURSE_SCHUR,
    NMOD_DET_PROFILE_RECURSE_KERNEL,
    NMOD_DET_PROFILE_COMBINE,
    NMOD_DET_PROFILE_BASECASE,
    NMOD_DET_PROFILE_PHASE_COUNT
};

enum
{
    NMOD_DET_BRANCH_NONE = -1,
    NMOD_DET_BRANCH_LEADING = 0,
    NMOD_DET_BRANCH_SCHUR,
    NMOD_DET_BRANCH_KERNEL,
    NMOD_DET_BRANCH_COUNT
};

#ifndef NMOD_DET_PROFILE_MAX_TRACK_DEPTH
#define NMOD_DET_PROFILE_MAX_TRACK_DEPTH 16
#endif

typedef struct
{
    slong root_dim;
    int use_parallel_requested;
    int top_level_parallel_used;
    int detail_depth_limit;
    slong node_count;
    slong basecase_count;
    slong max_depth;
    double wall_time;
    double accum[NMOD_DET_PROFILE_PHASE_COUNT];
    double top_level[NMOD_DET_PROFILE_PHASE_COUNT];
    slong nodes_by_depth[NMOD_DET_PROFILE_MAX_TRACK_DEPTH];
    double by_depth[NMOD_DET_PROFILE_MAX_TRACK_DEPTH][NMOD_DET_PROFILE_PHASE_COUNT];
    double branch_accum[NMOD_DET_BRANCH_COUNT][NMOD_DET_PROFILE_PHASE_COUNT];
    double branch_wall[NMOD_DET_BRANCH_COUNT];
} nmod_poly_mat_det_profile_t;

static double
_nmod_poly_mat_det_profile_now(void)
{
#ifdef _OPENMP
    return omp_get_wtime();
#else
    return ((double) clock()) / CLOCKS_PER_SEC;
#endif
}

static int
_nmod_poly_mat_det_profile_enabled(void)
{
    static int initialized = 0;
    static int enabled = 0;

    if (!initialized)
    {
        const char * env = getenv("DIXON_NMOD_DET_PROFILE");
        enabled = (env != NULL && env[0] != '\0' && env[0] != '0');
        initialized = 1;
    }

    return enabled;
}

static int
_nmod_poly_mat_det_profile_detail_depth(void)
{
    static int initialized = 0;
    static int depth_limit = 2;

    if (!initialized)
    {
        const char * env = getenv("DIXON_NMOD_DET_PROFILE_DEPTH");
        if (env != NULL && env[0] != '\0')
        {
            const long parsed = strtol(env, NULL, 10);
            if (parsed >= 0)
                depth_limit = (int) parsed;
        }
        initialized = 1;
    }

    return depth_limit;
}

static void
_nmod_poly_mat_det_profile_add(nmod_poly_mat_det_profile_t * profile,
                               int phase,
                               double elapsed,
                               slong depth,
                               int top_branch)
{
    if (profile == NULL || elapsed <= 0.0)
        return;

#ifdef _OPENMP
    #pragma omp atomic update
    profile->accum[phase] += elapsed;
#else
    profile->accum[phase] += elapsed;
#endif

    if (depth == 0)
    {
#ifdef _OPENMP
        #pragma omp atomic update
        profile->top_level[phase] += elapsed;
#else
        profile->top_level[phase] += elapsed;
 #endif
    }

    if (depth < NMOD_DET_PROFILE_MAX_TRACK_DEPTH)
    {
#ifdef _OPENMP
        #pragma omp atomic update
        profile->by_depth[depth][phase] += elapsed;
#else
        profile->by_depth[depth][phase] += elapsed;
#endif
    }

    if (top_branch >= 0 && top_branch < NMOD_DET_BRANCH_COUNT)
    {
#ifdef _OPENMP
        #pragma omp atomic update
        profile->branch_accum[top_branch][phase] += elapsed;
#else
        profile->branch_accum[top_branch][phase] += elapsed;
#endif
    }
}

static void
_nmod_poly_mat_det_profile_note_call(nmod_poly_mat_det_profile_t * profile,
                                     slong depth)
{
    if (profile == NULL)
        return;

#ifdef _OPENMP
    #pragma omp atomic update
    profile->node_count += 1;
#else
    profile->node_count += 1;
#endif

    if (depth < NMOD_DET_PROFILE_MAX_TRACK_DEPTH)
    {
#ifdef _OPENMP
        #pragma omp atomic update
        profile->nodes_by_depth[depth] += 1;
#else
        profile->nodes_by_depth[depth] += 1;
#endif
    }

    if (depth > profile->max_depth)
    {
#ifdef _OPENMP
        #pragma omp critical(nmod_det_profile_max_depth)
        {
            if (depth > profile->max_depth)
                profile->max_depth = depth;
        }
#else
        profile->max_depth = depth;
#endif
    }
}

static void
_nmod_poly_mat_det_profile_note_basecase(nmod_poly_mat_det_profile_t * profile)
{
    if (profile == NULL)
        return;

#ifdef _OPENMP
    #pragma omp atomic update
    profile->basecase_count += 1;
#else
    profile->basecase_count += 1;
#endif
}

static const char *
_nmod_poly_mat_det_profile_phase_name(int phase)
{
    static const char * names[NMOD_DET_PROFILE_PHASE_COUNT] =
    {
        "copy_split",
        "shift",
        "transpose",
        "kernel",
        "weak_popov",
        "pivot_select",
        "extract",
        "schur_mul",
        "recurse_leading",
        "recurse_schur",
        "recurse_kernel",
        "combine",
        "basecase"
    };

    return names[phase];
}

static const char *
_nmod_poly_mat_det_profile_branch_name(int branch)
{
    static const char * names[NMOD_DET_BRANCH_COUNT] =
    {
        "leading",
        "schur",
        "kernel_pivots"
    };

    return names[branch];
}

static void
_nmod_poly_mat_det_profile_print(const nmod_poly_mat_det_profile_t * profile,
                                 int ok)
{
    int top_phase = 0;
    double top_phase_time = profile->top_level[0];
    double top_prefix = 0.0;
    double top_children_max = 0.0;
    double top_children_sum = 0.0;

    for (int i = 0; i < NMOD_DET_PROFILE_PHASE_COUNT; i++)
    {
        if (profile->top_level[i] > top_phase_time)
        {
            top_phase_time = profile->top_level[i];
            top_phase = i;
        }
    }

    for (int i = NMOD_DET_PROFILE_COPY_SPLIT; i <= NMOD_DET_PROFILE_SCHUR_MUL; i++)
        top_prefix += profile->top_level[i];

    for (int i = NMOD_DET_PROFILE_RECURSE_LEADING; i <= NMOD_DET_PROFILE_RECURSE_KERNEL; i++)
    {
        top_children_sum += profile->top_level[i];
        if (profile->top_level[i] > top_children_max)
            top_children_max = profile->top_level[i];
    }

    fprintf(stderr,
            "[nmod_det_profile] dim=%ld ok=%d use_parallel=%d top_parallel=%d wall=%.6f nodes=%ld basecases=%ld max_depth=%ld\n",
            (long) profile->root_dim,
            ok,
            profile->use_parallel_requested,
            profile->top_level_parallel_used,
            profile->wall_time,
            (long) profile->node_count,
            (long) profile->basecase_count,
            (long) profile->max_depth);

    fprintf(stderr,
            "[nmod_det_profile] top_level prefix: copy_split=%.6f shift=%.6f transpose=%.6f kernel=%.6f weak_popov=%.6f pivot_select=%.6f extract=%.6f schur_mul=%.6f\n",
            profile->top_level[NMOD_DET_PROFILE_COPY_SPLIT],
            profile->top_level[NMOD_DET_PROFILE_SHIFT],
            profile->top_level[NMOD_DET_PROFILE_TRANSPOSE],
            profile->top_level[NMOD_DET_PROFILE_KERNEL],
            profile->top_level[NMOD_DET_PROFILE_WEAK_POPOV],
            profile->top_level[NMOD_DET_PROFILE_PIVOT_SELECT],
            profile->top_level[NMOD_DET_PROFILE_EXTRACT],
            profile->top_level[NMOD_DET_PROFILE_SCHUR_MUL]);

    fprintf(stderr,
            "[nmod_det_profile] top_level recurse: leading=%.6f schur=%.6f kernel_pivots=%.6f combine=%.6f\n",
            profile->top_level[NMOD_DET_PROFILE_RECURSE_LEADING],
            profile->top_level[NMOD_DET_PROFILE_RECURSE_SCHUR],
            profile->top_level[NMOD_DET_PROFILE_RECURSE_KERNEL],
            profile->top_level[NMOD_DET_PROFILE_COMBINE]);

    fprintf(stderr,
            "[nmod_det_profile] top_level summary: prefix_total=%.6f recurse_max=%.6f recurse_sum=%.6f hottest_phase=%s(%.6f)\n",
            top_prefix,
            top_children_max,
            top_children_sum,
            _nmod_poly_mat_det_profile_phase_name(top_phase),
            top_phase_time);

    for (int branch = 0; branch < NMOD_DET_BRANCH_COUNT; branch++)
    {
        int branch_phase = 0;
        double branch_phase_time = profile->branch_accum[branch][0];
        double branch_prefix = 0.0;

        for (int i = 0; i < NMOD_DET_PROFILE_PHASE_COUNT; i++)
            if (profile->branch_accum[branch][i] > branch_phase_time)
            {
                branch_phase_time = profile->branch_accum[branch][i];
                branch_phase = i;
            }

        for (int i = NMOD_DET_PROFILE_COPY_SPLIT; i <= NMOD_DET_PROFILE_SCHUR_MUL; i++)
            branch_prefix += profile->branch_accum[branch][i];

        fprintf(stderr,
                "[nmod_det_profile] branch=%s wall=%.6f prefix_total=%.6f kernel=%.6f schur_mul=%.6f recurse_leading=%.6f recurse_schur=%.6f recurse_kernel=%.6f basecase=%.6f hottest_phase=%s(%.6f)\n",
                _nmod_poly_mat_det_profile_branch_name(branch),
                profile->branch_wall[branch],
                branch_prefix,
                profile->branch_accum[branch][NMOD_DET_PROFILE_KERNEL],
                profile->branch_accum[branch][NMOD_DET_PROFILE_SCHUR_MUL],
                profile->branch_accum[branch][NMOD_DET_PROFILE_RECURSE_LEADING],
                profile->branch_accum[branch][NMOD_DET_PROFILE_RECURSE_SCHUR],
                profile->branch_accum[branch][NMOD_DET_PROFILE_RECURSE_KERNEL],
                profile->branch_accum[branch][NMOD_DET_PROFILE_BASECASE],
                _nmod_poly_mat_det_profile_phase_name(branch_phase),
                branch_phase_time);
    }

    {
        const slong print_depth = FLINT_MIN(profile->max_depth,
                                            (slong) FLINT_MIN(profile->detail_depth_limit,
                                                              NMOD_DET_PROFILE_MAX_TRACK_DEPTH - 1));
        for (slong depth = 0; depth <= print_depth; depth++)
        {
            int depth_phase = 0;
            double depth_phase_time = profile->by_depth[depth][0];
            double depth_prefix = 0.0;

            for (int i = 0; i < NMOD_DET_PROFILE_PHASE_COUNT; i++)
                if (profile->by_depth[depth][i] > depth_phase_time)
                {
                    depth_phase_time = profile->by_depth[depth][i];
                    depth_phase = i;
                }

            for (int i = NMOD_DET_PROFILE_COPY_SPLIT; i <= NMOD_DET_PROFILE_SCHUR_MUL; i++)
                depth_prefix += profile->by_depth[depth][i];

            fprintf(stderr,
                    "[nmod_det_profile] depth=%ld nodes=%ld prefix_total=%.6f kernel=%.6f schur_mul=%.6f recurse_leading=%.6f recurse_schur=%.6f recurse_kernel=%.6f basecase=%.6f hottest_phase=%s(%.6f)\n",
                    (long) depth,
                    (long) profile->nodes_by_depth[depth],
                    depth_prefix,
                    profile->by_depth[depth][NMOD_DET_PROFILE_KERNEL],
                    profile->by_depth[depth][NMOD_DET_PROFILE_SCHUR_MUL],
                    profile->by_depth[depth][NMOD_DET_PROFILE_RECURSE_LEADING],
                    profile->by_depth[depth][NMOD_DET_PROFILE_RECURSE_SCHUR],
                    profile->by_depth[depth][NMOD_DET_PROFILE_RECURSE_KERNEL],
                    profile->by_depth[depth][NMOD_DET_PROFILE_BASECASE],
                    _nmod_poly_mat_det_profile_phase_name(depth_phase),
                    depth_phase_time);
        }
    }

    fprintf(stderr, "[nmod_det_profile] accumulated_work:");
    for (int i = 0; i < NMOD_DET_PROFILE_PHASE_COUNT; i++)
        fprintf(stderr, " %s=%.6f", _nmod_poly_mat_det_profile_phase_name(i), profile->accum[i]);
    fprintf(stderr, "\n");
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
_nmod_poly_mat_det_generic_recursive(nmod_poly_t det,
                                     const nmod_poly_mat_t pmat,
                                     int use_parallel,
                                     slong depth,
                                     nmod_poly_mat_det_profile_t * profile,
                                     int top_branch)
{
    const slong dim = pmat->r;
    double t0, t1;

    _nmod_poly_mat_det_profile_note_call(profile, depth);

    if (pmat->r != pmat->c)
        return 0;

    if (dim == 0)
    {
        nmod_poly_one(det);
        return 1;
    }

    if (dim <= NMOD_POLY_MAT_DET_GENERIC_BASECASE)
    {
        t0 = _nmod_poly_mat_det_profile_now();
        nmod_poly_mat_det(det, pmat);
        t1 = _nmod_poly_mat_det_profile_now();
        _nmod_poly_mat_det_profile_note_basecase(profile);
        _nmod_poly_mat_det_profile_add(profile, NMOD_DET_PROFILE_BASECASE, t1 - t0, depth, top_branch);
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

    t0 = _nmod_poly_mat_det_profile_now();
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
    t1 = _nmod_poly_mat_det_profile_now();
    _nmod_poly_mat_det_profile_add(profile, NMOD_DET_PROFILE_COPY_SPLIT, t1 - t0, depth, top_branch);

    t0 = _nmod_poly_mat_det_profile_now();
    nmod_poly_mat_row_degree(shift, pmat_l, NULL);
    for (slong i = 0; i < dim; i++)
        if (shift[i] < 0)
            shift[i] = 0;
    t1 = _nmod_poly_mat_det_profile_now();
    _nmod_poly_mat_det_profile_add(profile, NMOD_DET_PROFILE_SHIFT, t1 - t0, depth, top_branch);

    t0 = _nmod_poly_mat_det_profile_now();
    nmod_poly_mat_transpose(pmat_l_t, pmat_l);
    t1 = _nmod_poly_mat_det_profile_now();
    _nmod_poly_mat_det_profile_add(profile, NMOD_DET_PROFILE_TRANSPOSE, t1 - t0, depth, top_branch);

    t0 = _nmod_poly_mat_det_profile_now();
    slong ker_dim = nmod_poly_mat_kernel_zls(kercols, kerdeg, pmat_l_t, shift, 2.0);
    t1 = _nmod_poly_mat_det_profile_now();
    _nmod_poly_mat_det_profile_add(profile, NMOD_DET_PROFILE_KERNEL, t1 - t0, depth, top_branch);
    if (ker_dim != cdim2)
        goto cleanup;

    t0 = _nmod_poly_mat_det_profile_now();
    _nmod_poly_mat_det_transpose_kernel_basis_rows(kerbas, kercols, cdim2);
    if (nmod_poly_mat_ordered_weak_popov_iter(kerbas, shift, NULL, pivind, NULL, ROW_LOWER) != cdim2)
    {
        t1 = _nmod_poly_mat_det_profile_now();
        _nmod_poly_mat_det_profile_add(profile, NMOD_DET_PROFILE_WEAK_POPOV, t1 - t0, depth, top_branch);
        goto cleanup;
    }
    t1 = _nmod_poly_mat_det_profile_now();
    _nmod_poly_mat_det_profile_add(profile, NMOD_DET_PROFILE_WEAK_POPOV, t1 - t0, depth, top_branch);

    t0 = _nmod_poly_mat_det_profile_now();
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
    t1 = _nmod_poly_mat_det_profile_now();
    _nmod_poly_mat_det_profile_add(profile, NMOD_DET_PROFILE_PIVOT_SELECT, t1 - t0, depth, top_branch);

    t0 = _nmod_poly_mat_det_profile_now();
    _nmod_poly_mat_det_extract_rows(leading_block, pmat_l, comp_rows, cdim1);
    _nmod_poly_mat_det_extract_columns(schur_input, kerbas, pivind, cdim2);
    t1 = _nmod_poly_mat_det_profile_now();
    _nmod_poly_mat_det_profile_add(profile, NMOD_DET_PROFILE_EXTRACT, t1 - t0, depth, top_branch);

    t0 = _nmod_poly_mat_det_profile_now();
    nmod_poly_mat_mul(schur, kerbas, pmat_r);
    t1 = _nmod_poly_mat_det_profile_now();
    _nmod_poly_mat_det_profile_add(profile, NMOD_DET_PROFILE_SCHUR_MUL, t1 - t0, depth, top_branch);

    {
        int ok_leading = 1;
        int ok_schur = 1;
        int ok_kernel_pivots = 1;
        double t_leading0, t_leading1, t_schur0, t_schur1, t_kernel0, t_kernel1;

#ifdef _OPENMP
        if (use_parallel
            && depth < NMOD_POLY_MAT_DET_GENERIC_MAX_PARALLEL_DEPTH
            && dim >= NMOD_POLY_MAT_DET_GENERIC_PARALLEL_THRESHOLD
            && omp_get_max_threads() > 1)
        {
            if (depth == 0 && profile != NULL)
                profile->top_level_parallel_used = 1;
            #pragma omp parallel sections num_threads(FLINT_MIN((slong) 3, (slong) omp_get_max_threads()))
            {
                #pragma omp section
                {
                    t_leading0 = _nmod_poly_mat_det_profile_now();
                    ok_leading = _nmod_poly_mat_det_generic_recursive(det_leading,
                                                                      leading_block,
                                                                      use_parallel,
                                                                      depth + 1,
                                                                      profile,
                                                                      (depth == 0) ? NMOD_DET_BRANCH_LEADING : top_branch);
                    t_leading1 = _nmod_poly_mat_det_profile_now();
                    _nmod_poly_mat_det_profile_add(profile, NMOD_DET_PROFILE_RECURSE_LEADING, t_leading1 - t_leading0, depth, top_branch);
                    if (depth == 0 && profile != NULL)
                        profile->branch_wall[NMOD_DET_BRANCH_LEADING] = t_leading1 - t_leading0;
                }
                #pragma omp section
                {
                    t_schur0 = _nmod_poly_mat_det_profile_now();
                    ok_schur = _nmod_poly_mat_det_generic_recursive(det_schur,
                                                                    schur,
                                                                    use_parallel,
                                                                    depth + 1,
                                                                    profile,
                                                                    (depth == 0) ? NMOD_DET_BRANCH_SCHUR : top_branch);
                    t_schur1 = _nmod_poly_mat_det_profile_now();
                    _nmod_poly_mat_det_profile_add(profile, NMOD_DET_PROFILE_RECURSE_SCHUR, t_schur1 - t_schur0, depth, top_branch);
                    if (depth == 0 && profile != NULL)
                        profile->branch_wall[NMOD_DET_BRANCH_SCHUR] = t_schur1 - t_schur0;
                }
                #pragma omp section
                {
                    t_kernel0 = _nmod_poly_mat_det_profile_now();
                    ok_kernel_pivots = _nmod_poly_mat_det_generic_recursive(det_kernel_pivots,
                                                                            schur_input,
                                                                            use_parallel,
                                                                            depth + 1,
                                                                            profile,
                                                                            (depth == 0) ? NMOD_DET_BRANCH_KERNEL : top_branch);
                    t_kernel1 = _nmod_poly_mat_det_profile_now();
                    _nmod_poly_mat_det_profile_add(profile, NMOD_DET_PROFILE_RECURSE_KERNEL, t_kernel1 - t_kernel0, depth, top_branch);
                    if (depth == 0 && profile != NULL)
                        profile->branch_wall[NMOD_DET_BRANCH_KERNEL] = t_kernel1 - t_kernel0;
                }
            }
        }
        else
#endif
        {
            t_leading0 = _nmod_poly_mat_det_profile_now();
            ok_leading = _nmod_poly_mat_det_generic_recursive(det_leading,
                                                              leading_block,
                                                              use_parallel,
                                                              depth + 1,
                                                              profile,
                                                              (depth == 0) ? NMOD_DET_BRANCH_LEADING : top_branch);
            t_leading1 = _nmod_poly_mat_det_profile_now();
            _nmod_poly_mat_det_profile_add(profile, NMOD_DET_PROFILE_RECURSE_LEADING, t_leading1 - t_leading0, depth, top_branch);
            if (depth == 0 && profile != NULL)
                profile->branch_wall[NMOD_DET_BRANCH_LEADING] = t_leading1 - t_leading0;

            t_schur0 = _nmod_poly_mat_det_profile_now();
            ok_schur = _nmod_poly_mat_det_generic_recursive(det_schur,
                                                            schur,
                                                            use_parallel,
                                                            depth + 1,
                                                            profile,
                                                            (depth == 0) ? NMOD_DET_BRANCH_SCHUR : top_branch);
            t_schur1 = _nmod_poly_mat_det_profile_now();
            _nmod_poly_mat_det_profile_add(profile, NMOD_DET_PROFILE_RECURSE_SCHUR, t_schur1 - t_schur0, depth, top_branch);
            if (depth == 0 && profile != NULL)
                profile->branch_wall[NMOD_DET_BRANCH_SCHUR] = t_schur1 - t_schur0;

            t_kernel0 = _nmod_poly_mat_det_profile_now();
            ok_kernel_pivots = _nmod_poly_mat_det_generic_recursive(det_kernel_pivots,
                                                                    schur_input,
                                                                    use_parallel,
                                                                    depth + 1,
                                                                    profile,
                                                                    (depth == 0) ? NMOD_DET_BRANCH_KERNEL : top_branch);
            t_kernel1 = _nmod_poly_mat_det_profile_now();
            _nmod_poly_mat_det_profile_add(profile, NMOD_DET_PROFILE_RECURSE_KERNEL, t_kernel1 - t_kernel0, depth, top_branch);
            if (depth == 0 && profile != NULL)
                profile->branch_wall[NMOD_DET_BRANCH_KERNEL] = t_kernel1 - t_kernel0;
        }

        if (!ok_leading || !ok_schur || !ok_kernel_pivots)
            goto cleanup;
    }
    if (nmod_poly_is_zero(det_kernel_pivots))
        goto cleanup;

    t0 = _nmod_poly_mat_det_profile_now();
    nmod_poly_mul(num, det_leading, det_schur);
    if (_nmod_poly_mat_det_permutation_sign_from_block_order(comp_rows, cdim1, pivind, cdim2) < 0)
        nmod_poly_neg(num, num);

    nmod_poly_divrem(det, rem, num, det_kernel_pivots);
    ok = nmod_poly_is_zero(rem);
    t1 = _nmod_poly_mat_det_profile_now();
    _nmod_poly_mat_det_profile_add(profile, NMOD_DET_PROFILE_COMBINE, t1 - t0, depth, top_branch);

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

int nmod_poly_mat_det_generic_with_opts(nmod_poly_t det,
                                        const nmod_poly_mat_t mat,
                                        int use_parallel)
{
    nmod_poly_mat_det_profile_t profile;
    nmod_poly_mat_det_profile_t * profile_ptr = NULL;
    double start_time = 0.0;
    int ok;

    if (_nmod_poly_mat_det_profile_enabled())
    {
        for (int i = 0; i < NMOD_DET_PROFILE_PHASE_COUNT; i++)
        {
            profile.accum[i] = 0.0;
            profile.top_level[i] = 0.0;
        }
        for (int depth = 0; depth < NMOD_DET_PROFILE_MAX_TRACK_DEPTH; depth++)
        {
            profile.nodes_by_depth[depth] = 0;
            for (int i = 0; i < NMOD_DET_PROFILE_PHASE_COUNT; i++)
                profile.by_depth[depth][i] = 0.0;
        }
        for (int branch = 0; branch < NMOD_DET_BRANCH_COUNT; branch++)
        {
            profile.branch_wall[branch] = 0.0;
            for (int i = 0; i < NMOD_DET_PROFILE_PHASE_COUNT; i++)
                profile.branch_accum[branch][i] = 0.0;
        }
        profile.root_dim = mat->r;
        profile.use_parallel_requested = use_parallel;
        profile.top_level_parallel_used = 0;
        profile.detail_depth_limit = _nmod_poly_mat_det_profile_detail_depth();
        profile.node_count = 0;
        profile.basecase_count = 0;
        profile.max_depth = 0;
        profile.wall_time = 0.0;
        profile_ptr = &profile;
        start_time = _nmod_poly_mat_det_profile_now();
    }

    ok = _nmod_poly_mat_det_generic_recursive(det, mat, use_parallel, 0, profile_ptr, NMOD_DET_BRANCH_NONE);

    if (profile_ptr != NULL)
    {
        profile.wall_time = _nmod_poly_mat_det_profile_now() - start_time;
        _nmod_poly_mat_det_profile_print(&profile, ok);
    }

    return ok;
}

int nmod_poly_mat_det_generic(nmod_poly_t det,
                              const nmod_poly_mat_t mat)
{
    return nmod_poly_mat_det_generic_with_opts(det, mat, 0);
}

int nmod_poly_mat_det_generic_knowing_degree_with_opts(nmod_poly_t det,
                                                       const nmod_poly_mat_t mat,
                                                       slong degree,
                                                       int use_parallel)
{
    const int ok = nmod_poly_mat_det_generic_with_opts(det, mat, use_parallel);
    return ok && (degree == nmod_poly_degree(det));
}

int nmod_poly_mat_det_generic_knowing_degree(nmod_poly_t det,
                                             const nmod_poly_mat_t mat,
                                             slong degree)
{
    return nmod_poly_mat_det_generic_knowing_degree_with_opts(det, mat, degree, 0);
}