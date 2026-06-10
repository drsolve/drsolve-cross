#include <flint/nmod_vec.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "nmod_mat_poly.h"
#include "nmod_mat_extra.h"

extern int g_dixon_verbose_level;
extern int g_dixon_debug_mode;

typedef struct
{
    slong value;
    slong index;
} slong_pair;

typedef struct
{
    slong calls;
    slong zero_input_returns;
    slong nullity_zero_breaks;
    slong partial_nullity_updates;
    slong full_nullity_skips;
    slong max_order;
    slong max_rows;
    slong max_cols;
    slong max_rank;
    slong max_nullity;
    slong gemm_calls;
    double total_time;
    double init_time;
    double residual_init_time;
    double permute_time;
    double nullspace_time;
    double shift_time;
    double appbas_update_time;
    double appbas_update_mul_time;
    double residual_update_time;
    double residual_update_mul_time;
    double row_shift_time;
    double cleanup_time;
} nmod_mat_poly_mbasis_resupdate_profile_t;

static nmod_mat_poly_mbasis_resupdate_profile_t g_nmod_mat_poly_mbasis_resupdate_profile;

static double
_nmod_resupdate_now_seconds(void)
{
    return ((double) clock()) / CLOCKS_PER_SEC;
}

static int
_nmod_resupdate_profile_enabled(void)
{
    return g_dixon_verbose_level >= 2 || g_dixon_debug_mode;
}

void
nmod_mat_poly_mbasis_resupdate_profile_reset(void)
{
    g_nmod_mat_poly_mbasis_resupdate_profile =
        (nmod_mat_poly_mbasis_resupdate_profile_t) {0};
}

void
nmod_mat_poly_mbasis_resupdate_profile_print(void)
{
    if (!_nmod_resupdate_profile_enabled())
        return;

    const nmod_mat_poly_mbasis_resupdate_profile_t *p =
        &g_nmod_mat_poly_mbasis_resupdate_profile;

    printf("        [mbasis_resupdate profile]\n");
    printf("          calls=%ld zero_input=%ld nullity_zero_breaks=%ld partial_updates=%ld full_nullity_skips=%ld max_order=%ld max_shape=%ldx%ld max_rank=%ld max_nullity=%ld\n",
           p->calls, p->zero_input_returns, p->nullity_zero_breaks,
           p->partial_nullity_updates, p->full_nullity_skips,
           p->max_order, p->max_rows, p->max_cols, p->max_rank, p->max_nullity);
    printf("          total=%.6fs init=%.6fs residual_init=%.6fs permute=%.6fs nullspace=%.6fs shift=%.6fs appbas_update=%.6fs residual_update=%.6fs row_shift=%.6fs cleanup=%.6fs\n",
           p->total_time, p->init_time, p->residual_init_time, p->permute_time,
           p->nullspace_time, p->shift_time, p->appbas_update_time,
           p->residual_update_time, p->row_shift_time, p->cleanup_time);
    printf("          update mul: appbas=%.6fs residual=%.6fs gemm_calls=%ld\n",
           p->appbas_update_mul_time, p->residual_update_mul_time,
           p->gemm_calls);
}

static int
_slong_pair_compare_resupdate(const void * a, const void * b)
{
    slong_pair aa = *(const slong_pair *) a;
    slong_pair bb = *(const slong_pair *) b;

    if (aa.value == bb.value)
    {
        if (aa.index < bb.index)
            return -1;
        if (aa.index > bb.index)
            return 1;
        return 0;
    }

    return (aa.value < bb.value) ? -1 : 1;
}

static inline void
_find_shift_permutation_resupdate(slong * perm,
                                  const slong * shift,
                                  slong n,
                                  slong_pair * pair_tmp)
{
    for (slong i = 0; i < n; i++)
    {
        pair_tmp[i].value = shift[i];
        pair_tmp[i].index = i;
    }

    qsort(pair_tmp, n, sizeof(slong_pair), _slong_pair_compare_resupdate);

    for (slong i = 0; i < n; i++)
        perm[i] = pair_tmp[i].index;
}

static void
_nmod_mat_poly_resupdate_add_compact_kernel_rows(nmod_mat_t mat,
                                                 const nmod_mat_t nsbas,
                                                 slong rank,
                                                 slong nullity,
                                                 nmod_mat_t scratch_mul)
{
    nmod_mat_t mat_win;
    nmod_mat_window_init(mat_win, mat, 0, 0, rank, mat->c);
    nmod_mat_mul(scratch_mul, nsbas, mat_win);

    for (slong i = 0; i < nullity; ++i)
    {
        _nmod_vec_add(nmod_mat_entry_ptr(mat, rank + i, 0),
                      nmod_mat_entry_ptr(mat, rank + i, 0),
                      nmod_mat_entry_ptr(scratch_mul, i, 0),
                      mat->c,
                      mat->mod);
    }

    nmod_mat_window_clear(mat_win);
}

static void
_nmod_mat_poly_resupdate_shift_top_rows(nmod_mat_poly_t matp,
                                        slong src_lo,
                                        slong rank)
{
    if (rank == 0 || src_lo >= matp->length)
        return;

    for (slong d = matp->length - 1; d > src_lo; --d)
        for (slong i = 0; i < rank; ++i)
            _nmod_vec_set(nmod_mat_entry_ptr(matp->coeffs + d, i, 0),
                          nmod_mat_entry_ptr(matp->coeffs + (d - 1), i, 0),
                          matp->c);

    for (slong i = 0; i < rank; ++i)
        _nmod_vec_zero(nmod_mat_entry_ptr(matp->coeffs + src_lo, i, 0), matp->c);
}

void
nmod_mat_poly_mbasis_resupdate(nmod_mat_poly_t appbas,
                               slong * shift,
                               const nmod_mat_poly_t matp,
                               slong order)
{
    const double call_start = _nmod_resupdate_now_seconds();
    const int collect_profile = _nmod_resupdate_profile_enabled();
    double t0 = 0.0;
    const slong m = matp->r;
    const slong n = matp->c;

    if (collect_profile)
    {
        g_nmod_mat_poly_mbasis_resupdate_profile.calls++;
        if (order > g_nmod_mat_poly_mbasis_resupdate_profile.max_order)
            g_nmod_mat_poly_mbasis_resupdate_profile.max_order = order;
        if (m > g_nmod_mat_poly_mbasis_resupdate_profile.max_rows)
            g_nmod_mat_poly_mbasis_resupdate_profile.max_rows = m;
        if (n > g_nmod_mat_poly_mbasis_resupdate_profile.max_cols)
            g_nmod_mat_poly_mbasis_resupdate_profile.max_cols = n;
    }

    nmod_mat_poly_one(appbas);

    if (nmod_mat_poly_is_zero(matp) || order <= 0)
    {
        if (collect_profile)
        {
            g_nmod_mat_poly_mbasis_resupdate_profile.zero_input_returns++;
            g_nmod_mat_poly_mbasis_resupdate_profile.total_time +=
                _nmod_resupdate_now_seconds() - call_start;
        }
        return;
    }

    t0 = _nmod_resupdate_now_seconds();
    nmod_mat_poly_t residuals;
    nmod_mat_poly_init2(residuals, m, n, matp->mod.n, order);
    _nmod_mat_poly_set_length(residuals, order);
    for (slong d = 0; d < order; ++d)
    {
        if (d < matp->length)
            nmod_mat_set(residuals->coeffs + d, matp->coeffs + d);
        else
            nmod_mat_zero(residuals->coeffs + d);
    }
    if (collect_profile)
        g_nmod_mat_poly_mbasis_resupdate_profile.residual_init_time +=
            _nmod_resupdate_now_seconds() - t0;

    t0 = _nmod_resupdate_now_seconds();
    nmod_mat_t res;
    nmod_mat_init(res, m, n, matp->mod.n);

    slong * perm = _perm_init(m);
    slong * pivots = (slong *) flint_malloc(m * sizeof(slong));
    slong * pivots_inv = (slong *) flint_malloc(m * sizeof(slong));
    slong_pair * pair_tmp = (slong_pair *) flint_malloc(m * sizeof(slong_pair));
    nmod_mat_t nsbas;
    nmod_mat_t scratch_mul;
    int nsbas_initialized = 0;
    int scratch_mul_initialized = 0;
    if (collect_profile)
        g_nmod_mat_poly_mbasis_resupdate_profile.init_time +=
            _nmod_resupdate_now_seconds() - t0;

    for (slong ord = 0; ord < order; ++ord)
    {
        slong nullity, rank;

        t0 = _nmod_resupdate_now_seconds();
        _find_shift_permutation_resupdate(perm, shift, m, pair_tmp);
        nmod_mat_set(res, residuals->coeffs + ord);
        nmod_mat_permute_rows(res, perm, NULL);
        if (collect_profile)
            g_nmod_mat_poly_mbasis_resupdate_profile.permute_time +=
                _nmod_resupdate_now_seconds() - t0;

        t0 = _nmod_resupdate_now_seconds();
        if (nsbas_initialized)
            nmod_mat_clear(nsbas);
        nullity = nmod_mat_left_nullspace_compact(nsbas, pivots, res);
        nsbas_initialized = 1;
        rank = m - nullity;
        if (collect_profile)
        {
            if (rank > g_nmod_mat_poly_mbasis_resupdate_profile.max_rank)
                g_nmod_mat_poly_mbasis_resupdate_profile.max_rank = rank;
            if (nullity > g_nmod_mat_poly_mbasis_resupdate_profile.max_nullity)
                g_nmod_mat_poly_mbasis_resupdate_profile.max_nullity = nullity;
            g_nmod_mat_poly_mbasis_resupdate_profile.nullspace_time +=
                _nmod_resupdate_now_seconds() - t0;
        }

        if (nullity == 0)
        {
            if (collect_profile)
                g_nmod_mat_poly_mbasis_resupdate_profile.nullity_zero_breaks++;
            nmod_mat_poly_shift_left(appbas, appbas, order - ord);
            for (slong i = 0; i < m; ++i)
                shift[i] += order - ord;
            break;
        }

        if (nullity == m)
        {
            if (collect_profile)
                g_nmod_mat_poly_mbasis_resupdate_profile.full_nullity_skips++;
            continue;
        }

        if (collect_profile)
            g_nmod_mat_poly_mbasis_resupdate_profile.partial_nullity_updates++;

        if (!scratch_mul_initialized ||
            scratch_mul->r != nullity ||
            scratch_mul->c != appbas->c)
        {
            if (scratch_mul_initialized)
                nmod_mat_clear(scratch_mul);
            nmod_mat_init(scratch_mul, nullity, appbas->c, matp->mod.n);
            scratch_mul_initialized = 1;
        }

        t0 = _nmod_resupdate_now_seconds();
        _perm_compose(pivots, perm, pivots, m);
        _perm_inv(pivots_inv, pivots, m);
        for (slong i = 0; i < rank; ++i)
            shift[pivots[i]] += 1;
        if (collect_profile)
            g_nmod_mat_poly_mbasis_resupdate_profile.shift_time +=
                _nmod_resupdate_now_seconds() - t0;

        t0 = _nmod_resupdate_now_seconds();
        nmod_mat_poly_permute_rows(appbas, pivots, NULL);
        if (collect_profile)
            g_nmod_mat_poly_mbasis_resupdate_profile.permute_time +=
                _nmod_resupdate_now_seconds() - t0;

        t0 = _nmod_resupdate_now_seconds();
        for (slong d = 0; d < appbas->length; ++d)
        {
            double tm = _nmod_resupdate_now_seconds();
            _nmod_mat_poly_resupdate_add_compact_kernel_rows(appbas->coeffs + d,
                                                             nsbas,
                                                             rank,
                                                             nullity,
                                                             scratch_mul);
            if (collect_profile)
            {
                g_nmod_mat_poly_mbasis_resupdate_profile.appbas_update_mul_time +=
                    _nmod_resupdate_now_seconds() - tm;
                g_nmod_mat_poly_mbasis_resupdate_profile.gemm_calls++;
            }
        }
        if (collect_profile)
            g_nmod_mat_poly_mbasis_resupdate_profile.appbas_update_time +=
                _nmod_resupdate_now_seconds() - t0;

        t0 = _nmod_resupdate_now_seconds();
        for (slong i = 0; i < rank; ++i)
        {
            if (!_nmod_vec_is_zero(nmod_mat_poly_entry_ptr(appbas, appbas->length - 1, i, 0), m))
            {
                nmod_mat_poly_fit_length(appbas, appbas->length + 1);
                _nmod_mat_poly_set_length(appbas, appbas->length + 1);
                break;
            }
        }

        for (slong d = appbas->length - 1; d > 0; --d)
            for (slong i = 0; i < rank; ++i)
                _nmod_vec_set(nmod_mat_poly_entry_ptr(appbas, d, i, 0),
                              nmod_mat_poly_entry_ptr(appbas, d - 1, i, 0),
                              m);
        for (slong i = 0; i < rank; ++i)
            _nmod_vec_zero(nmod_mat_poly_entry_ptr(appbas, 0, i, 0), m);
        if (collect_profile)
            g_nmod_mat_poly_mbasis_resupdate_profile.row_shift_time +=
                _nmod_resupdate_now_seconds() - t0;

        t0 = _nmod_resupdate_now_seconds();
        nmod_mat_poly_permute_rows(appbas, pivots_inv, NULL);
        if (collect_profile)
            g_nmod_mat_poly_mbasis_resupdate_profile.permute_time +=
                _nmod_resupdate_now_seconds() - t0;

        t0 = _nmod_resupdate_now_seconds();
        for (slong d = ord + 1; d < order; ++d)
        {
            double tm;
            if (!scratch_mul_initialized ||
                scratch_mul->r != nullity ||
                scratch_mul->c != residuals->c)
            {
                if (scratch_mul_initialized)
                    nmod_mat_clear(scratch_mul);
                nmod_mat_init(scratch_mul, nullity, residuals->c, matp->mod.n);
                scratch_mul_initialized = 1;
            }
            nmod_mat_permute_rows(residuals->coeffs + d, pivots, NULL);
            tm = _nmod_resupdate_now_seconds();
            _nmod_mat_poly_resupdate_add_compact_kernel_rows(residuals->coeffs + d,
                                                             nsbas,
                                                             rank,
                                                             nullity,
                                                             scratch_mul);
            if (collect_profile)
            {
                g_nmod_mat_poly_mbasis_resupdate_profile.residual_update_mul_time +=
                    _nmod_resupdate_now_seconds() - tm;
                g_nmod_mat_poly_mbasis_resupdate_profile.gemm_calls++;
            }
        }
        nmod_mat_permute_rows(residuals->coeffs + ord, pivots, NULL);
        if (collect_profile)
            g_nmod_mat_poly_mbasis_resupdate_profile.residual_update_time +=
                _nmod_resupdate_now_seconds() - t0;

        t0 = _nmod_resupdate_now_seconds();
        _nmod_mat_poly_resupdate_shift_top_rows(residuals, ord, rank);
        for (slong d = ord; d < order; ++d)
            nmod_mat_permute_rows(residuals->coeffs + d, pivots_inv, NULL);
        if (collect_profile)
        {
            g_nmod_mat_poly_mbasis_resupdate_profile.row_shift_time +=
                _nmod_resupdate_now_seconds() - t0;
            g_nmod_mat_poly_mbasis_resupdate_profile.permute_time += 0.0;
        }
    }

    t0 = _nmod_resupdate_now_seconds();
    nmod_mat_clear(res);
    if (nsbas_initialized)
        nmod_mat_clear(nsbas);
    if (scratch_mul_initialized)
        nmod_mat_clear(scratch_mul);
    nmod_mat_poly_clear(residuals);
    _perm_clear(perm);
    flint_free(pivots);
    flint_free(pivots_inv);
    flint_free(pair_tmp);
    if (collect_profile)
    {
        g_nmod_mat_poly_mbasis_resupdate_profile.cleanup_time +=
            _nmod_resupdate_now_seconds() - t0;
        g_nmod_mat_poly_mbasis_resupdate_profile.total_time +=
            _nmod_resupdate_now_seconds() - call_start;
    }
}
