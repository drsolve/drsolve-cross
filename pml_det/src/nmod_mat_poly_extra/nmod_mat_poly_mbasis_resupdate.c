#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nmod_mat_poly.h"
#include "nmod_mat_extra.h"

extern int g_dixon_verbose_level;
extern int g_dixon_debug_mode;

typedef struct
{
    slong value;
    slong index;
} slong_pair;

static inline int
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
                                                 slong nullity)
{
    nmod_mat_t ns_mul, mat_win;
    nmod_mat_init(ns_mul, nullity, mat->c, mat->mod.n);
    nmod_mat_window_init(mat_win, mat, 0, 0, rank, mat->c);
    nmod_mat_mul(ns_mul, nsbas, mat_win);

    for (slong i = 0; i < nullity; ++i)
    {
        _nmod_vec_add(nmod_mat_entry_ptr(mat, rank + i, 0),
                      nmod_mat_entry_ptr(mat, rank + i, 0),
                      nmod_mat_entry_ptr(ns_mul, i, 0),
                      mat->c,
                      mat->mod);
    }

    nmod_mat_window_clear(mat_win);
    nmod_mat_clear(ns_mul);
}

static void
_nmod_mat_poly_resupdate_shift_top_rows(nmod_mat_poly_t matp,
                                        slong src_lo,
                                        slong rank)
{
    if (rank == 0 || src_lo >= matp->length)
        return;

    for (slong d = matp->length - 1; d > src_lo; --d)
    {
        for (slong i = 0; i < rank; ++i)
        {
            _nmod_vec_set(nmod_mat_entry_ptr(matp->coeffs + d, i, 0),
                          nmod_mat_entry_ptr(matp->coeffs + (d - 1), i, 0),
                          matp->c);
        }
    }

    for (slong i = 0; i < rank; ++i)
        _nmod_vec_zero(nmod_mat_entry_ptr(matp->coeffs + src_lo, i, 0), matp->c);
}

void
nmod_mat_poly_mbasis_resupdate(nmod_mat_poly_t appbas,
                               slong * shift,
                               const nmod_mat_poly_t matp,
                               slong order)
{
    const slong m = matp->r;
    const slong n = matp->c;

    nmod_mat_poly_one(appbas);

    if (nmod_mat_poly_is_zero(matp) || order <= 0)
        return;

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

    nmod_mat_t res;
    nmod_mat_init(res, m, n, matp->mod.n);

    slong * perm = _perm_init(m);
    slong * pivots = (slong *) flint_malloc(m * sizeof(slong));
    slong * pivots_inv = (slong *) flint_malloc(m * sizeof(slong));
    slong_pair * pair_tmp = (slong_pair *) flint_malloc(m * sizeof(slong_pair));

    nmod_mat_t nsbas;

    for (slong ord = 0; ord < order; ++ord)
    {
        slong nullity, rank;

        _find_shift_permutation_resupdate(perm, shift, m, pair_tmp);
        nmod_mat_set(res, residuals->coeffs + ord);
        nmod_mat_permute_rows(res, perm, NULL);

        if (ord > 0)
            nmod_mat_clear(nsbas);
        nullity = nmod_mat_left_nullspace_compact(nsbas, pivots, res);
        rank = m - nullity;

        if (nullity == 0)
        {
            nmod_mat_poly_shift_left(appbas, appbas, order - ord);
            for (slong i = 0; i < m; ++i)
                shift[i] += order - ord;
            break;
        }

        if (nullity == m)
            continue;

        _perm_compose(pivots, perm, pivots, m);
        _perm_inv(pivots_inv, pivots, m);

        for (slong i = 0; i < rank; ++i)
            shift[pivots[i]] += 1;

        nmod_mat_poly_permute_rows(appbas, pivots, NULL);

        for (slong d = 0; d < appbas->length; ++d)
            _nmod_mat_poly_resupdate_add_compact_kernel_rows(appbas->coeffs + d, nsbas, rank, nullity);

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

        nmod_mat_poly_permute_rows(appbas, pivots_inv, NULL);

        for (slong d = ord + 1; d < order; ++d)
        {
            nmod_mat_permute_rows(residuals->coeffs + d, pivots, NULL);
            _nmod_mat_poly_resupdate_add_compact_kernel_rows(residuals->coeffs + d, nsbas, rank, nullity);
        }

        nmod_mat_permute_rows(residuals->coeffs + ord, pivots, NULL);

        _nmod_mat_poly_resupdate_shift_top_rows(residuals, ord, rank);

        for (slong d = ord; d < order; ++d)
            nmod_mat_permute_rows(residuals->coeffs + d, pivots_inv, NULL);
    }

    nmod_mat_clear(res);
    nmod_mat_clear(nsbas);
    nmod_mat_poly_clear(residuals);
    _perm_clear(perm);
    flint_free(pivots);
    flint_free(pivots_inv);
    flint_free(pair_tmp);

    if (0 && g_dixon_debug_mode && g_dixon_verbose_level >= 3)
    {
        fprintf(stderr,
                "  [PML mbasis_resupdate] experimental FLINT residual-update path used.\n");
    }
}
