#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>
#include <flint/perm.h>
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

/* Kernel-row update, applied to all matrix coefficients coeffs[d_lo..d_hi-1]
 * at once. In terms of the row indices pivots[] (pivots[0..rank-1] = pivot
 * rows, pivots[rank..rank+nullity-1] = kernel rows), it performs, for each d:
 *
 *     coeffs[d][pivots[rank+i], :] += sum_j nsbas[i][j] * coeffs[d][pivots[j], :]
 *
 * Instead of one small nmod_mat_mul per coefficient (as previously done,
 * with the extra cost of permuting each full coefficient back and forth),
 * the pivot rows of many consecutive coefficients are gathered side by side
 * into a single rank x (chunk*ncols) matrix, multiplied ONCE by nsbas, and
 * the product is scattered-added back into the kernel rows. This turns
 * (d_hi-d_lo) small products into a handful of large ones -- exactly the same
 * field operations, but far less per-call overhead in FLINT's nmod_mat_mul
 * and much better use of its blocking/vectorized kernels. The gathered block
 * is capped at ~8MB so memory stays bounded for large shapes.
 */
static void
_resupdate_batched_kernel_update(nmod_mat_struct * coeffs,
                                 slong d_lo,
                                 slong d_hi,
                                 const nmod_mat_t nsbas,
                                 const slong * pivots,
                                 slong rank,
                                 slong nullity,
                                 slong ncols,
                                 ulong modn,
                                 nmod_t mod,
                                 slong * gemm_counter)
{
    const slong nb = d_hi - d_lo;
    if (nb <= 0 || rank <= 0 || nullity <= 0 || ncols <= 0)
        return;

    /* number of coefficients gathered per product: pivot block <= 2^20 entries */
    slong chunk = ((slong) 1 << 20) / (rank * ncols);
    if (chunk < 1)
        chunk = 1;
    if (chunk > nb)
        chunk = nb;

    for (slong start = d_lo; start < d_hi; start += chunk)
    {
        const slong len = FLINT_MIN(chunk, d_hi - start);

        nmod_mat_t gath, prod;
        nmod_mat_init(gath, rank, len * ncols, modn);
        nmod_mat_init(prod, nullity, len * ncols, modn);

        for (slong dd = 0; dd < len; dd++)
            for (slong j = 0; j < rank; j++)
                _nmod_vec_set(nmod_mat_entry_ptr(gath, j, dd * ncols),
                              nmod_mat_entry_ptr(coeffs + start + dd, pivots[j], 0),
                              ncols);

        nmod_mat_mul(prod, nsbas, gath);
        if (gemm_counter)
            (*gemm_counter)++;

        for (slong dd = 0; dd < len; dd++)
            for (slong i = 0; i < nullity; i++)
                _nmod_vec_add(nmod_mat_entry_ptr(coeffs + start + dd, pivots[rank + i], 0),
                              nmod_mat_entry_ptr(coeffs + start + dd, pivots[rank + i], 0),
                              nmod_mat_entry_ptr(prod, i, dd * ncols),
                              ncols,
                              mod);

        nmod_mat_clear(prod);
        nmod_mat_clear(gath);
    }
}

void
nmod_mat_poly_mbasis_resupdate(nmod_mat_poly_t appbas,
                               slong * shift,
                               const nmod_mat_poly_t matp,
                               slong order)
{
    const int collect_profile = _nmod_resupdate_profile_enabled();
    const double call_start = collect_profile ? _nmod_resupdate_now_seconds() : 0.0;
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

    if (collect_profile)
        t0 = _nmod_resupdate_now_seconds();
    nmod_mat_poly_t residuals;
    nmod_mat_poly_init2(residuals, m, n, matp->mod.n, order);
    _nmod_mat_poly_set_length(residuals, order);
    for (slong d = 0; d < order; ++d)
    {
        if (d < matp->length)
            nmod_mat_set(residuals->coeffs + d, matp->coeffs + d);
        /* _nmod_mat_poly_set_length zero-initializes the coefficients
         * beyond matp->length, nothing to do for them */
    }
    if (collect_profile)
        g_nmod_mat_poly_mbasis_resupdate_profile.residual_init_time +=
            _nmod_resupdate_now_seconds() - t0;

    if (collect_profile)
        t0 = _nmod_resupdate_now_seconds();
    nmod_mat_t res;
    nmod_mat_init(res, m, n, matp->mod.n);

    slong * perm = _perm_init(m);
    slong * pivots = (slong *) flint_malloc(m * sizeof(slong));
    slong_pair * pair_tmp = (slong_pair *) flint_malloc(m * sizeof(slong_pair));
    nmod_mat_t nsbas;
    int nsbas_initialized = 0;
    if (collect_profile)
        g_nmod_mat_poly_mbasis_resupdate_profile.init_time +=
            _nmod_resupdate_now_seconds() - t0;

    for (slong ord = 0; ord < order; ++ord)
    {
        slong nullity, rank;

        /* res <- rows of residuals[ord] sorted by nondecreasing shift.
         * Direct row gather: replaces the former full copy + full row
         * permutation of the coefficient. */
        if (collect_profile)
            t0 = _nmod_resupdate_now_seconds();
        _find_shift_permutation_resupdate(perm, shift, m, pair_tmp);
        for (slong i = 0; i < m; ++i)
            _nmod_vec_set(nmod_mat_entry_ptr(res, i, 0),
                          nmod_mat_entry_ptr(residuals->coeffs + ord, perm[i], 0),
                          n);
        if (collect_profile)
            g_nmod_mat_poly_mbasis_resupdate_profile.permute_time +=
                _nmod_resupdate_now_seconds() - t0;

        if (collect_profile)
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

        /* pivots[] in original row indices: pivots[0..rank-1] are the pivot
         * rows (row rank profile of res, pulled back through perm),
         * pivots[rank..m-1] the kernel rows. All updates below address rows
         * through this list directly -- no coefficient of appbas or of the
         * residual is ever permuted. */
        if (collect_profile)
            t0 = _nmod_resupdate_now_seconds();
        _perm_compose(pivots, perm, pivots, m);
        for (slong i = 0; i < rank; ++i)
            shift[pivots[i]] += 1;
        if (collect_profile)
            g_nmod_mat_poly_mbasis_resupdate_profile.shift_time +=
                _nmod_resupdate_now_seconds() - t0;

        /* appbas: kernel-row update on all coefficients, one batched product */
        if (collect_profile)
            t0 = _nmod_resupdate_now_seconds();
        _resupdate_batched_kernel_update(appbas->coeffs, 0, appbas->length,
                                         nsbas, pivots, rank, nullity,
                                         appbas->c, matp->mod.n, matp->mod,
                                         collect_profile ?
                                         &g_nmod_mat_poly_mbasis_resupdate_profile.gemm_calls
                                         : NULL);
        if (collect_profile)
        {
            const double dt = _nmod_resupdate_now_seconds() - t0;
            g_nmod_mat_poly_mbasis_resupdate_profile.appbas_update_time += dt;
            g_nmod_mat_poly_mbasis_resupdate_profile.appbas_update_mul_time += dt;
        }

        /* multiply pivot rows of appbas by x, growing the length if needed */
        if (collect_profile)
            t0 = _nmod_resupdate_now_seconds();
        for (slong i = 0; i < rank; ++i)
        {
            if (!_nmod_vec_is_zero(nmod_mat_poly_entry_ptr(appbas, appbas->length - 1, pivots[i], 0), appbas->c))
            {
                nmod_mat_poly_fit_length(appbas, appbas->length + 1);
                _nmod_mat_poly_set_length(appbas, appbas->length + 1);
                break;
            }
        }

        for (slong d = appbas->length - 1; d > 0; --d)
            for (slong i = 0; i < rank; ++i)
                _nmod_vec_set(nmod_mat_poly_entry_ptr(appbas, d, pivots[i], 0),
                              nmod_mat_poly_entry_ptr(appbas, d - 1, pivots[i], 0),
                              appbas->c);
        for (slong i = 0; i < rank; ++i)
            _nmod_vec_zero(nmod_mat_poly_entry_ptr(appbas, 0, pivots[i], 0), appbas->c);
        if (collect_profile)
            g_nmod_mat_poly_mbasis_resupdate_profile.row_shift_time +=
                _nmod_resupdate_now_seconds() - t0;

        /* residual: same kernel-row update on the remaining coefficients,
         * batched as well */
        if (collect_profile)
            t0 = _nmod_resupdate_now_seconds();
        _resupdate_batched_kernel_update(residuals->coeffs, ord + 1, order,
                                         nsbas, pivots, rank, nullity,
                                         n, matp->mod.n, matp->mod,
                                         collect_profile ?
                                         &g_nmod_mat_poly_mbasis_resupdate_profile.gemm_calls
                                         : NULL);
        if (collect_profile)
        {
            const double dt = _nmod_resupdate_now_seconds() - t0;
            g_nmod_mat_poly_mbasis_resupdate_profile.residual_update_time += dt;
            g_nmod_mat_poly_mbasis_resupdate_profile.residual_update_mul_time += dt;
        }

        /* multiply pivot rows of the residual by x: shift them one
         * coefficient up, starting from coefficient ord (whose pivot rows
         * are consumed here and never read again) */
        if (collect_profile)
            t0 = _nmod_resupdate_now_seconds();
        for (slong d = order - 1; d > ord; --d)
            for (slong i = 0; i < rank; ++i)
                _nmod_vec_set(nmod_mat_entry_ptr(residuals->coeffs + d, pivots[i], 0),
                              nmod_mat_entry_ptr(residuals->coeffs + (d - 1), pivots[i], 0),
                              n);
        if (collect_profile)
            g_nmod_mat_poly_mbasis_resupdate_profile.row_shift_time +=
                _nmod_resupdate_now_seconds() - t0;
    }

    if (collect_profile)
        t0 = _nmod_resupdate_now_seconds();
    nmod_mat_clear(res);
    if (nsbas_initialized)
        nmod_mat_clear(nsbas);
    nmod_mat_poly_clear(residuals);
    _perm_clear(perm);
    flint_free(pivots);
    flint_free(pair_tmp);
    if (collect_profile)
    {
        g_nmod_mat_poly_mbasis_resupdate_profile.cleanup_time +=
            _nmod_resupdate_now_seconds() - t0;
        g_nmod_mat_poly_mbasis_resupdate_profile.total_time +=
            _nmod_resupdate_now_seconds() - call_start;
    }
}