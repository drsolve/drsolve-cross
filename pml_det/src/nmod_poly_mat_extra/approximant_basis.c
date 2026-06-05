/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "nmod_poly_mat_multiply.h"  // for middle product
#include "nmod_poly_mat_approximant.h"
#include "nmod_mat_poly.h"
#include "nmod_poly_mat_utils.h"

extern int g_dixon_verbose_level;
extern int g_dixon_debug_mode;

typedef struct
{
    slong mbasis_calls;
    slong pmbasis_calls;
    slong pmbasis_basecase_calls;
    slong max_depth;
    slong max_order;
    slong max_rows;
    slong max_cols;
    double mbasis_total_time;
    double pmbasis_total_time;
    double pmbasis_basecase_time;
    double pmbasis_first_half_time;
    double pmbasis_middle_product_time;
    double pmbasis_second_half_time;
    double pmbasis_final_mul_time;
} nmod_poly_mat_pmbasis_profile_t;

static nmod_poly_mat_pmbasis_profile_t g_nmod_pmbasis_profile;
static slong g_nmod_pmbasis_depth = 0;

static double
_nmod_pmbasis_now_seconds(void)
{
    return ((double) clock()) / CLOCKS_PER_SEC;
}

static int
_nmod_pmbasis_profile_enabled(void)
{
    return g_dixon_verbose_level >= 3 && g_dixon_debug_mode;
}

void
nmod_poly_mat_pmbasis_profile_reset(void)
{
    g_nmod_pmbasis_profile = (nmod_poly_mat_pmbasis_profile_t) {0};
    nmod_mat_poly_mbasis_profile_reset();
}

void
nmod_poly_mat_pmbasis_profile_print(void)
{
    if (!_nmod_pmbasis_profile_enabled())
        return;

    const nmod_poly_mat_pmbasis_profile_t *p = &g_nmod_pmbasis_profile;
    printf("      [pmbasis profile]\n");
    printf("        pmbasis_calls=%ld basecase=%ld mbasis_calls=%ld max_depth=%ld max_order=%ld max_shape=%ldx%ld\n",
           p->pmbasis_calls, p->pmbasis_basecase_calls, p->mbasis_calls,
           p->max_depth, p->max_order, p->max_rows, p->max_cols);
    printf("        total=%.6fs basecase=%.6fs first_half=%.6fs middle_product=%.6fs second_half=%.6fs final_mul=%.6fs mbasis_total=%.6fs\n",
           p->pmbasis_total_time, p->pmbasis_basecase_time,
           p->pmbasis_first_half_time, p->pmbasis_middle_product_time,
           p->pmbasis_second_half_time, p->pmbasis_final_mul_time,
           p->mbasis_total_time);
    nmod_mat_poly_mbasis_profile_print();
}

void nmod_poly_mat_mbasis(nmod_poly_mat_t appbas,
                          slong * shift,
                          const nmod_poly_mat_t pmat,
                          ulong order)
{
    const double call_start = _nmod_pmbasis_now_seconds();
    const char *mbasis_mode = getenv("DRSOLVE_PML_MBASIS_MODE");
    nmod_mat_poly_t app, matp;
    if (_nmod_pmbasis_profile_enabled()) {
        g_nmod_pmbasis_profile.mbasis_calls++;
        if ((slong) order > g_nmod_pmbasis_profile.max_order)
            g_nmod_pmbasis_profile.max_order = (slong) order;
        if (pmat->r > g_nmod_pmbasis_profile.max_rows)
            g_nmod_pmbasis_profile.max_rows = pmat->r;
        if (pmat->c > g_nmod_pmbasis_profile.max_cols)
            g_nmod_pmbasis_profile.max_cols = pmat->c;
    }
    nmod_mat_poly_init(matp, pmat->r, pmat->c, pmat->modulus);
    nmod_mat_poly_set_trunc_from_poly_mat(matp, pmat, order);
    nmod_mat_poly_init(app, pmat->r, pmat->r, pmat->modulus);
    if (mbasis_mode != NULL &&
        (strcmp(mbasis_mode, "legacy") == 0 ||
         strcmp(mbasis_mode, "classic") == 0 ||
         strcmp(mbasis_mode, "rescomp") == 0))
        nmod_mat_poly_mbasis(app, shift, matp, order);
    else
        nmod_mat_poly_mbasis_resupdate(app, shift, matp, order);
    nmod_poly_mat_set_from_mat_poly(appbas, app);
    nmod_mat_poly_clear(matp);
    nmod_mat_poly_clear(app);
    if (_nmod_pmbasis_profile_enabled())
        g_nmod_pmbasis_profile.mbasis_total_time += _nmod_pmbasis_now_seconds() - call_start;
}

void nmod_poly_mat_pmbasis(nmod_poly_mat_t appbas,
                           slong * shift,
                           const nmod_poly_mat_t pmat,
                           slong order)
{
    const double call_start = _nmod_pmbasis_now_seconds();
    const int collect_profile = _nmod_pmbasis_profile_enabled();
    const slong depth = g_nmod_pmbasis_depth++;

    if (collect_profile) {
        g_nmod_pmbasis_profile.pmbasis_calls++;
        if (depth > g_nmod_pmbasis_profile.max_depth)
            g_nmod_pmbasis_profile.max_depth = depth;
        if (order > g_nmod_pmbasis_profile.max_order)
            g_nmod_pmbasis_profile.max_order = order;
        if (pmat->r > g_nmod_pmbasis_profile.max_rows)
            g_nmod_pmbasis_profile.max_rows = pmat->r;
        if (pmat->c > g_nmod_pmbasis_profile.max_cols)
            g_nmod_pmbasis_profile.max_cols = pmat->c;
    }

    if (order <= PMBASIS_THRES)
    {
        if (collect_profile)
            g_nmod_pmbasis_profile.pmbasis_basecase_calls++;
        {
            double t0 = _nmod_pmbasis_now_seconds();
        nmod_poly_mat_mbasis(appbas, shift, pmat, order);
            if (collect_profile)
                g_nmod_pmbasis_profile.pmbasis_basecase_time += _nmod_pmbasis_now_seconds() - t0;
        }
        if (collect_profile)
            g_nmod_pmbasis_profile.pmbasis_total_time += _nmod_pmbasis_now_seconds() - call_start;
        g_nmod_pmbasis_depth--;
        return;
    }

    const long order1 = order>>1;
    const long order2 = order - order1;
    nmod_poly_mat_t appbas2, residual;

    nmod_poly_mat_init(appbas2, pmat->r, pmat->r, pmat->modulus);
    nmod_poly_mat_init(residual, pmat->r, pmat->c, pmat->modulus);

    {
        double t0 = _nmod_pmbasis_now_seconds();
    nmod_poly_mat_pmbasis(appbas, shift, pmat, order1);
        if (collect_profile)
            g_nmod_pmbasis_profile.pmbasis_first_half_time += _nmod_pmbasis_now_seconds() - t0;
    }

    {
        double t0 = _nmod_pmbasis_now_seconds();
    nmod_poly_mat_middle_product_naive(residual, appbas, pmat, order1, order2-1);
        if (collect_profile)
            g_nmod_pmbasis_profile.pmbasis_middle_product_time += _nmod_pmbasis_now_seconds() - t0;
    }

    {
        double t0 = _nmod_pmbasis_now_seconds();
    nmod_poly_mat_pmbasis(appbas2, shift, residual, order2);
        if (collect_profile)
            g_nmod_pmbasis_profile.pmbasis_second_half_time += _nmod_pmbasis_now_seconds() - t0;
    }

    {
        double t0 = _nmod_pmbasis_now_seconds();
    nmod_poly_mat_mul(appbas, appbas2, appbas);
        if (collect_profile)
            g_nmod_pmbasis_profile.pmbasis_final_mul_time += _nmod_pmbasis_now_seconds() - t0;
    }

    nmod_poly_mat_clear(appbas2);
    nmod_poly_mat_clear(residual);
    if (collect_profile)
        g_nmod_pmbasis_profile.pmbasis_total_time += _nmod_pmbasis_now_seconds() - call_start;
    g_nmod_pmbasis_depth--;
}
