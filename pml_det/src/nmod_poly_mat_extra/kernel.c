/*
    Copyright (C) 2025 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_kernel.h"

extern int g_dixon_verbose_level;
extern int g_dixon_debug_mode;

typedef struct
{
    slong kernel_calls;
    slong zls_calls;
    slong zls_zero_residue_returns;
    slong zls_m1_returns;
    slong zls_recursive_fallback_returns;
    slong zls_recursive_success_returns;
    slong max_depth;
    slong max_m;
    slong max_n;
    slong max_order;
    double kernel_total_time;
    double kernel_copy_time;
    double kernel_sort_time;
    double kernel_zls_time;
    double kernel_unpermute_time;
    double zls_total_time;
    double zls_shift_prep_time;
    double zls_transpose_time;
    double zls_pmbasis_time;
    double zls_residual_mul_time;
    double zls_zero_split_time;
    double zls_p2_sort_time;
    double zls_shift_right_split_time;
    double zls_recursive_time;
    double zls_middle_mul_time;
    double zls_combine_time;
} nmod_poly_mat_kernel_zls_profile_t;

static nmod_poly_mat_kernel_zls_profile_t g_nmod_kernel_zls_profile;
static slong g_nmod_kernel_zls_depth = 0;

static double
_nmod_kernel_now_seconds(void)
{
    return ((double) clock()) / CLOCKS_PER_SEC;
}

static int
_nmod_kernel_profile_enabled(void)
{
    return g_dixon_verbose_level >= 3 && g_dixon_debug_mode;
}

void
nmod_poly_mat_kernel_zls_profile_reset(void)
{
    g_nmod_kernel_zls_profile = (nmod_poly_mat_kernel_zls_profile_t) {0};
    nmod_poly_mat_pmbasis_profile_reset();
}

void
nmod_poly_mat_kernel_zls_profile_print(void)
{
    if (!_nmod_kernel_profile_enabled())
        return;

    const nmod_poly_mat_kernel_zls_profile_t *p = &g_nmod_kernel_zls_profile;
    printf("    [kernel_zls profile]\n");
    printf("      kernel_calls=%ld zls_calls=%ld max_depth=%ld max_shape=%ldx%ld max_order=%ld\n",
           p->kernel_calls, p->zls_calls, p->max_depth, p->max_m, p->max_n, p->max_order);
    printf("      returns: zero_residue=%ld m1=%ld recursive_fallback=%ld recursive_success=%ld\n",
           p->zls_zero_residue_returns, p->zls_m1_returns,
           p->zls_recursive_fallback_returns, p->zls_recursive_success_returns);
    printf("      kernel_total=%.6fs copy=%.6fs sort=%.6fs zls=%.6fs unpermute=%.6fs\n",
           p->kernel_total_time, p->kernel_copy_time, p->kernel_sort_time,
           p->kernel_zls_time, p->kernel_unpermute_time);
    printf("      zls_total=%.6fs shift_prep=%.6fs transpose=%.6fs pmbasis=%.6fs residual_mul=%.6fs zero_split=%.6fs\n",
           p->zls_total_time, p->zls_shift_prep_time, p->zls_transpose_time,
           p->zls_pmbasis_time, p->zls_residual_mul_time, p->zls_zero_split_time);
    printf("      p2_sort=%.6fs shift_right_split=%.6fs recurse=%.6fs middle_mul=%.6fs combine=%.6fs\n",
           p->zls_p2_sort_time, p->zls_shift_right_split_time, p->zls_recursive_time,
           p->zls_middle_mul_time, p->zls_combine_time);
    nmod_poly_mat_pmbasis_profile_print();
}

/**
 *
 * The content of the file has not yet been tested very extensively
 *
 */

/**
 *
 * In place, M is m x n, shifted column degrees
 *   permute the columns in nondecreasing order
 *
 * Todo/to see: zero is set to degree 0 (not -1), for specific use
 *     in ZLS algorithm for the kernel
 *
 * Input:
 *  M  m x n
 *  ishift[m]
 *
 * Output:
 *  M is modified in place
 *  perm[n],initialized outside, the carried out permutation
 *  sdeg[n] initialized outside
 *
 */

// Dummy function for qsort
 int cmp(const void *a, const void *b)
{
    const slong *x = a;
    const slong *y = b;

    if (*x < *y) return -1;
    if (*x > *y) return 1;
    return 0;
}

void _nmod_poly_mat_sort_permute_columns_zls(nmod_poly_mat_t M, slong *sdeg, \
                                            slong *perm, const slong *ishift)
{

    slong n = M->c;
    slong j;

    nmod_poly_mat_column_degree(sdeg, M, ishift);

    for (j=0; j<n; j++)
        if (sdeg[j] < 0) sdeg[j]=0;

    _nmod_poly_mat_permute_columns_by_sorting_vec(M, n, sdeg, perm);

}

/**
 *
 *  Internal function
 *
 *  Right shifted kernel of a polynomial matrix, assuming that the columns has been
 *     sorted by shifted degree
 *
 *  Algorithm of Wei Zhou, George Labahn, and Arne Storjohann
 *   "Computing Minimal Nullspace Bases"
 *    ISSAC 2012, https://dl.acm.org/doi/abs/10.1145/2442829.2442881
 *
 *  TODO/TO SEE:
 *
 *  Input:
 *    A in m x n
 *    ishift[n], NULL (the degrees are computed) or initialized outside,
 *      the shift for the kernel
 *      values should be at least 0 (even for zero columns in A
 *      "with entries arranged in non-decreasing order and bounding the
 *       corresponding column degrees of A."
 *    kappa, a double >= 2, for the order of the order bases
 *              kappa * s instead of 3 *s in ZLS
 *
 *  Output:
 *    returns the dimension w of the kernel, which may be zero
 *    N, n x w, is initialized here only when w > 0,
 *           should be freed outside in that case only,
 *            gives a minimal basis of the right kernel
 *    degN[n], initialized outside, its first w entries are concerned,
 *        they are the ishift shift degrees of the kernel basis
 *
 */

int nmod_poly_mat_zls_sorted(nmod_poly_mat_t N, slong *degN, const nmod_poly_mat_t A, \
                                 const slong *ishift, const double kappa)
{
    const double call_start = _nmod_kernel_now_seconds();
    const int collect_profile = _nmod_kernel_profile_enabled();
    const slong depth = g_nmod_kernel_zls_depth++;

    slong i,j,k;

    slong m = A->r;
    slong n = A->c;


    slong min_mn;
    slong rho=0;

    if (collect_profile) {
        g_nmod_kernel_zls_profile.zls_calls++;
        if (depth > g_nmod_kernel_zls_profile.max_depth)
            g_nmod_kernel_zls_profile.max_depth = depth;
        if (m > g_nmod_kernel_zls_profile.max_m)
            g_nmod_kernel_zls_profile.max_m = m;
        if (n > g_nmod_kernel_zls_profile.max_n)
            g_nmod_kernel_zls_profile.max_n = n;
    }

    // Tuning the order of the subsequent approximant
    // ----------------------------------------------

    /* Temporary shift workspace.
       Use slong, not plain long: on Windows/LLP64 FLINT defines slong as
       long long while long stays 32-bit, and qsort(..., sizeof(slong), ...)
       would otherwise corrupt memory. */
    slong * shift = FLINT_ARRAY_ALLOC(n, slong);


    // In order to sort for computing the order of the approximant, for the test m=1 to be ok
    {
        double t0 = _nmod_kernel_now_seconds();
        for (i=0; i<n; i++)
        {
            shift[i]=ishift[i];
        }
        qsort(shift, n, sizeof(slong), cmp);
        if (collect_profile)
            g_nmod_kernel_zls_profile.zls_shift_prep_time += _nmod_kernel_now_seconds() - t0;
    }

    if (m <= n)
    {
        min_mn=m;
        for (i=n-m; i<n; i++)
            rho+=shift[i];
    }
    else
    {   // No sort needed, we consider all degrees
        min_mn=n;
        for (i=0; i<n; i++)
            rho+=ishift[i];
        // printf("\n m > n  %ld  %ld",m,n);   // To see, rectangular case
    }


    slong s;
    s = ceil((double) rho/min_mn);

    // Transposed A for the approximant PT on the left
    // PT will then be use without transposing

    nmod_poly_mat_t AT;
    nmod_poly_mat_init(AT, n, m, A->modulus);
    {
        double t0 = _nmod_kernel_now_seconds();
        nmod_poly_mat_transpose(AT,A);
        if (collect_profile)
            g_nmod_kernel_zls_profile.zls_transpose_time += _nmod_kernel_now_seconds() - t0;
    }

    nmod_poly_mat_t PT;
    nmod_poly_mat_init(PT, n, n, A->modulus);

    // Approximant PT
    // shift is modified in place
    // --------------------------

    for (i=0; i<n; i++)
    {
        shift[i]=ishift[i];
    }

    slong ks;
    ks=ceil((double) kappa*s);
    if (collect_profile && ks + 1 > g_nmod_kernel_zls_profile.max_order)
        g_nmod_kernel_zls_profile.max_order = ks + 1;
    {
        double t0 = _nmod_kernel_now_seconds();
        nmod_poly_mat_pmbasis(PT, shift, AT, ks+1);
        if (collect_profile)
            g_nmod_kernel_zls_profile.zls_pmbasis_time += _nmod_kernel_now_seconds() - t0;
    }

    // Looking for zero residues and non zero residues
    //  the global residue is reused later
    // -----------------------------------------------

    nmod_poly_mat_t RT;
    nmod_poly_mat_init(RT, n, m, A->modulus);
    {
        double t0 = _nmod_kernel_now_seconds();
        nmod_poly_mat_mul(RT,PT,AT);
        if (collect_profile)
            g_nmod_kernel_zls_profile.zls_residual_mul_time += _nmod_kernel_now_seconds() - t0;
    }

    slong cdeg[n];
    nmod_poly_mat_row_degree(cdeg, RT, NULL);

    nmod_poly_mat_clear(AT);

    // Zero residue matrix P1  n x n1
    // nonzero residues matrix P2  n x n2
    //-----------------------------------

    slong n1=0;
    slong n2;

    {
        double t0 = _nmod_kernel_now_seconds();
        for (j=0; j<n; j++) {
            if (cdeg[j]<0)
                n1+=1;
        }
        if (collect_profile)
            g_nmod_kernel_zls_profile.zls_zero_split_time += _nmod_kernel_now_seconds() - t0;
    }

    nmod_poly_mat_t P1;

    if (n1>0) {
        nmod_poly_mat_init(P1, n, n1, A->modulus);

        k=0;
        for (j = 0; j < n; j++)
        {
            if (cdeg[j]<0) {

                for (i = 0; i < n; i++)
                    nmod_poly_set(nmod_poly_mat_entry(P1, i, k), nmod_poly_mat_entry(PT, j, i));
                k+=1;
            }
        }
    }

    n2=n-n1;

    // the kernel is found
    if (n2==0) {
        nmod_poly_mat_init_set(N,P1);
        nmod_poly_mat_column_degree(degN, P1, ishift);
        if (collect_profile) {
            g_nmod_kernel_zls_profile.zls_zero_residue_returns++;
            g_nmod_kernel_zls_profile.zls_total_time += _nmod_kernel_now_seconds() - call_start;
        }

        nmod_poly_mat_clear(P1);
        nmod_poly_mat_clear(RT);
        nmod_poly_mat_clear(PT);
        flint_free(shift);
        g_nmod_kernel_zls_depth--;
        return n1;
    }

    nmod_poly_mat_t P2;
    nmod_poly_mat_init(P2, n, n2, A->modulus);

    k=0;
    for (j = 0; j < n; j++)
    {
        if (cdeg[j]>=0) {

            for (i = 0; i < n; i++)
                nmod_poly_set(nmod_poly_mat_entry(P2, i, k), nmod_poly_mat_entry(PT, j, i));
            k+=1;
        }

    }

    //nmod_poly_mat_clear(PT);

    //  Special case m=1
    //  before the divide and conquer on m
    //  the kernel is found
    // -----------------------------------

    if (m==1){

        if (n1==0) {
            if (collect_profile) {
                g_nmod_kernel_zls_profile.zls_m1_returns++;
                g_nmod_kernel_zls_profile.zls_total_time += _nmod_kernel_now_seconds() - call_start;
            }
            nmod_poly_mat_clear(P2);
            nmod_poly_mat_clear(RT);
            nmod_poly_mat_clear(PT);
            flint_free(shift);
            g_nmod_kernel_zls_depth--;
            return 0;
        }
        else {

            nmod_poly_mat_init_set(N,P1);
            nmod_poly_mat_column_degree(degN, P1, ishift);
            if (collect_profile) {
                g_nmod_kernel_zls_profile.zls_m1_returns++;
                g_nmod_kernel_zls_profile.zls_total_time += _nmod_kernel_now_seconds() - call_start;
            }

            nmod_poly_mat_clear(P1);
            nmod_poly_mat_clear(P2);
            nmod_poly_mat_clear(RT);
            nmod_poly_mat_clear(PT);
            flint_free(shift);
            g_nmod_kernel_zls_depth--;
            return n1;
        }
    }

    //  Now n2 <> 0 and m > 1
    //    one can proceed to the divide and conquer from the rows
    //    of the nonzero residue P2, and of A.P2
    // --------------------------------------------------------

    slong * perm = flint_malloc(n2 * sizeof(slong));


    slong degP2[n2];

    // perm will be used also later for G from RT
    {
        double t0 = _nmod_kernel_now_seconds();
        _nmod_poly_mat_sort_permute_columns_zls(P2,degP2,perm,ishift);
        if (collect_profile)
            g_nmod_kernel_zls_profile.zls_p2_sort_time += _nmod_kernel_now_seconds() - t0;
    }

    for (i = 0; i < n2; i++) {
        degP2[i]=degP2[i]-ks;  // used below for the recursive calls
        if (degP2[i] < 0)
            degP2[i] = 0;
    }


    //+++++++++ OLD

    // nmod_poly_mat_t G;
    // nmod_poly_mat_init(G, m, n2, A->modulus);

    // nmod_poly_mat_mul(G, A, P2);

    // nmod_poly_mat_t TT;
    // nmod_poly_mat_init(TT, m, n2, A->modulus);


    // nmod_poly_mat_shift_right(TT,G,ks);


    // slong new_m=floor((double) m/2);

    // nmod_poly_mat_t G1;
    // nmod_poly_mat_init(G1, new_m, n2, A->modulus);

    // for (i = 0; i < new_m; i++){
    //     for (j = 0; j < n2; j++) {
    //         nmod_poly_set(nmod_poly_mat_entry(G1, i, j), nmod_poly_mat_entry(TT, i, j));
    //     }
    // }

    // nmod_poly_mat_t G2;
    // nmod_poly_mat_init(G2, m-new_m, n2, A->modulus);

    // for (i = 0; i < m-new_m; i++) {
    //     for (j = 0; j < n2; j++) {
    //         nmod_poly_set(nmod_poly_mat_entry(G2, i, j), nmod_poly_mat_entry(TT, i+new_m, j));
    //     }
    // }


    //++++++++  END OLD


    //++++++++++ NEW ++++++++++++++
    nmod_poly_mat_t TT;
    nmod_poly_mat_init(TT, m, n2, A->modulus);

    // We extract G (temporary TT) from RT as we have been extracting P2 from PT
    //   and apply above ordering of P2
    {
        double t0 = _nmod_kernel_now_seconds();
        k=0;
        for (j = 0; j < n; j++)
        {
            if (cdeg[j]>=0) {

                for (i = 0; i < m; i++)
                    nmod_poly_set(nmod_poly_mat_entry(TT, i, k), nmod_poly_mat_entry(RT, j, i));
                k+=1;
            }

        }
        nmod_poly_mat_clear(RT);

        nmod_poly_mat_permute_columns(TT, perm, NULL);

        nmod_poly_mat_t G;
        nmod_poly_mat_init(G, m, n2, A->modulus);

        nmod_poly_mat_shift_right(G,TT,ks);


    // We split G for the recursive call
    // ---------------------------------

        slong new_m=floor((double) m/2);

        nmod_poly_mat_t G1;
        nmod_poly_mat_init(G1, new_m, n2, A->modulus);

        for (i = 0; i < new_m; i++){
            for (j = 0; j < n2; j++) {
                nmod_poly_set(nmod_poly_mat_entry(G1, i, j), nmod_poly_mat_entry(G, i, j));
            }
        }

        nmod_poly_mat_t G2;
        nmod_poly_mat_init(G2, m-new_m, n2, A->modulus);

        for (i = 0; i < m-new_m; i++) {
            for (j = 0; j < n2; j++) {
                nmod_poly_set(nmod_poly_mat_entry(G2, i, j), nmod_poly_mat_entry(G, i+new_m, j));
            }
        }
        if (collect_profile)
            g_nmod_kernel_zls_profile.zls_shift_right_split_time += _nmod_kernel_now_seconds() - t0;

        nmod_poly_mat_clear(TT);
        nmod_poly_mat_clear(G);
    //++++++++++++++  END NEW ++++++++++++++++++

    // Recursive calls
    // ---------------

        nmod_poly_mat_t N1;
        nmod_poly_mat_t N2;

        slong c1=0;
        slong c2=0;


        {
            double t1 = _nmod_kernel_now_seconds();
            c1=nmod_poly_mat_zls_sorted(N1, degN, G1, degP2, kappa);
            if (collect_profile)
                g_nmod_kernel_zls_profile.zls_recursive_time += _nmod_kernel_now_seconds() - t1;
        }


        if (c1 != 0) {

            nmod_poly_mat_t G3;
            nmod_poly_mat_init(G3, m-new_m, c1, A->modulus);

            {
                double t1 = _nmod_kernel_now_seconds();
                nmod_poly_mat_mul(G3, G2, N1);
                if (collect_profile)
                    g_nmod_kernel_zls_profile.zls_middle_mul_time += _nmod_kernel_now_seconds() - t1;
            }

            for (i=0; i<c1; i++) {
                shift[i]=degN[i];
            }

            {
                double t1 = _nmod_kernel_now_seconds();
                c2=nmod_poly_mat_zls_sorted(N2, degN, G3, shift, kappa);
                if (collect_profile)
                    g_nmod_kernel_zls_profile.zls_recursive_time += _nmod_kernel_now_seconds() - t1;
            }
            nmod_poly_mat_clear(G3);

        }


        nmod_poly_mat_clear(G1);
        nmod_poly_mat_clear(G2);

    // the recursive calls did not provide anything more

        if ((c1==0) || (c2==0)){

            if (c1 != 0)
                nmod_poly_mat_clear(N1);
            if (c2 != 0)
                nmod_poly_mat_clear(N2);
            flint_free(perm);
            nmod_poly_mat_clear(PT);

            if (collect_profile)
                g_nmod_kernel_zls_profile.zls_recursive_fallback_returns++;

            if (n1==0) {
                nmod_poly_mat_clear(P2);
                if (n1 > 0)
                    nmod_poly_mat_clear(P1);
                flint_free(shift);
                if (collect_profile)
                    g_nmod_kernel_zls_profile.zls_total_time += _nmod_kernel_now_seconds() - call_start;
                g_nmod_kernel_zls_depth--;
                return 0;
            }
            else {
                nmod_poly_mat_init_set(N,P1);
                nmod_poly_mat_column_degree(degN, P1, ishift);

                nmod_poly_mat_clear(P1);
                nmod_poly_mat_clear(P2);
                flint_free(shift);
                if (collect_profile)
                    g_nmod_kernel_zls_profile.zls_total_time += _nmod_kernel_now_seconds() - call_start;
                g_nmod_kernel_zls_depth--;
                return n1;
            }
        }

    // A new part of the kernel has been found by the recursive calls, Q
    // -----------------------------------------------------------------


        nmod_poly_mat_t Q1;
        nmod_poly_mat_init(Q1, n, c1, A->modulus);

        {
            double t1 = _nmod_kernel_now_seconds();
            nmod_poly_mat_mul(Q1, P2, N1);
            if (collect_profile)
                g_nmod_kernel_zls_profile.zls_middle_mul_time += _nmod_kernel_now_seconds() - t1;
        }

        nmod_poly_mat_t Q;
        nmod_poly_mat_init(Q, n, c2, A->modulus);

        {
            double t1 = _nmod_kernel_now_seconds();
            nmod_poly_mat_mul(Q, Q1, N2);
            if (collect_profile)
                g_nmod_kernel_zls_profile.zls_combine_time += _nmod_kernel_now_seconds() - t1;
        }

        nmod_poly_mat_clear(N1);
        nmod_poly_mat_clear(N2);
        nmod_poly_mat_clear(Q1);
        nmod_poly_mat_clear(P2);
        flint_free(perm);
        nmod_poly_mat_clear(PT);

        if (collect_profile)
            g_nmod_kernel_zls_profile.zls_recursive_success_returns++;

        if (n1 ==0) {

            nmod_poly_mat_init_set(N,Q); // We should not need to copy
            nmod_poly_mat_column_degree(degN, Q, ishift);

            nmod_poly_mat_clear(Q);
            if (n1 > 0)
                nmod_poly_mat_clear(P1);
            flint_free(shift);
            if (collect_profile)
                g_nmod_kernel_zls_profile.zls_total_time += _nmod_kernel_now_seconds() - call_start;
            g_nmod_kernel_zls_depth--;
            return c2;

        }
        else {
            nmod_poly_mat_init(N, n, n1+c2, A->modulus);

            {
                double t1 = _nmod_kernel_now_seconds();
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n1; j++) {
                        nmod_poly_set(nmod_poly_mat_entry(N, i, j), nmod_poly_mat_entry(P1,i,j));
                    }
                }

                for (i = 0; i < n; i++) {
                    for (j = 0; j < c2; j++) {
                        nmod_poly_set(nmod_poly_mat_entry(N, i, j+n1), nmod_poly_mat_entry(Q, i, j));
                    }
                }
                if (collect_profile)
                    g_nmod_kernel_zls_profile.zls_combine_time += _nmod_kernel_now_seconds() - t1;
            }

            slong odeg[n1+c2];

            nmod_poly_mat_column_degree(odeg, P1, ishift);

            for (i=0; i<n1; i++) {
                degN[i]=odeg[i];
            }

            nmod_poly_mat_column_degree(odeg, Q, ishift);

            for (i=0; i<c2; i++) {
                degN[i+n1]=odeg[i];
            }

            nmod_poly_mat_clear(P1);
            nmod_poly_mat_clear(Q);

            flint_free(shift);
            if (collect_profile)
                g_nmod_kernel_zls_profile.zls_total_time += _nmod_kernel_now_seconds() - call_start;
            g_nmod_kernel_zls_depth--;
            return n1+c2;
        }

    }

    flint_free(shift);
    if (collect_profile)
        g_nmod_kernel_zls_profile.zls_total_time += _nmod_kernel_now_seconds() - call_start;
    g_nmod_kernel_zls_depth--;
    return 0;
}


/**
 *
 *  Right shifted kernel of a polynomial matrix, assuming that the columns has been
 *     sorted by shifted degree
 *
 *  Algorithm of Wei Zhou, George Labahn, and Arne Storjohann
 *   "Computing Minimal Nullspace Bases"
 *    ISSAC 2012, https://dl.acm.org/doi/abs/10.1145/2442829.2442881
 *
 *  Calls nmod_poly_mat_zls_sorted after an initial sorting
 *
 *  TODO/TO SEE:
 *
 * Input:
 *    iA in m x n
 *     ishift[n], NULL (the degrees are computed) or initialized outside,
 *      the shift for the kernel
 *      values should be at least 0 (even for zero columns in A
 *      "with entries arranged in non-decreasing order and bounding the
 *       corresponding column degrees of A."
 *    kappa, a double >= 2, for the order of the order bases
 *              kappa * s instead of 3 *s in ZLS
 *
 *  Output:
 *    returns the dimension w of the kernel, which may be zero
 *    N, is initialized  n x n outside
 *       its first w columns give a minimal basis of the kernel
 *    degN[n], initialized outside, its first w entries are concerned,
 *        they are the ishift shifted degrees of the kernel basis
 *
 */

/**
 *  Calls nmod_poly_mat_zls_sorted after an initial sorting
 */

int nmod_poly_mat_kernel_zls(nmod_poly_mat_t N, slong *degN, const nmod_poly_mat_t iA, \
                         const slong *ishift, const double kappa)
{
    const double call_start = _nmod_kernel_now_seconds();
    const int collect_profile = _nmod_kernel_profile_enabled();

    slong j,k;

    slong n = iA->c;

    if (collect_profile)
        g_nmod_kernel_zls_profile.kernel_calls++;


    nmod_poly_mat_t A;
    {
        double t0 = _nmod_kernel_now_seconds();
        nmod_poly_mat_init_set(A, iA);
        if (collect_profile)
            g_nmod_kernel_zls_profile.kernel_copy_time += _nmod_kernel_now_seconds() - t0;
    }

    slong * perm = flint_malloc(n * sizeof(slong));
    slong sdeg[n];

    // No input shift, simply the column degrees, then ordered
    //   a permutation is carried out in place
    // -------------------------------------------------------

    if (ishift == NULL) {
        double t0 = _nmod_kernel_now_seconds();
        nmod_poly_mat_column_degree(sdeg, A, NULL);

        for (j=0; j<n; j++)
            if (sdeg[j] < 0) sdeg[j]=0;

        _nmod_poly_mat_permute_columns_by_sorting_vec(A, n, sdeg, perm);
        if (collect_profile)
            g_nmod_kernel_zls_profile.kernel_sort_time += _nmod_kernel_now_seconds() - t0;
    }
    // Input shift, we sort
    // --------------------
    else {
        double t0 = _nmod_kernel_now_seconds();

        for (j=0; j<n; j++)
            sdeg[j]=ishift[j];

        _nmod_poly_mat_permute_columns_by_sorting_vec(A, n, sdeg, perm);
        if (collect_profile)
            g_nmod_kernel_zls_profile.kernel_sort_time += _nmod_kernel_now_seconds() - t0;
    }

    // Call to ZLS
    // -----------

    slong nz;
    slong tdeg[n];

    nmod_poly_mat_t NT;
    {
        double t0 = _nmod_kernel_now_seconds();
        nz=nmod_poly_mat_zls_sorted(NT, tdeg, A, sdeg, kappa);
        if (collect_profile)
            g_nmod_kernel_zls_profile.kernel_zls_time += _nmod_kernel_now_seconds() - t0;
    }


    // Undo the permutation for the kernel of the input matrix
    // -------------------------------------------------------

    if (nz !=0) {
        double t0 = _nmod_kernel_now_seconds();

        for (k = 0; k < n; k++) {
            for (j = 0; j < nz; j++){
                nmod_poly_set(nmod_poly_mat_entry(N, perm[k], j), nmod_poly_mat_entry(NT,k,j));
                degN[perm[k]]=tdeg[k];
            }
        }
        if (collect_profile)
            g_nmod_kernel_zls_profile.kernel_unpermute_time += _nmod_kernel_now_seconds() - t0;

        nmod_poly_mat_clear(NT); // Since nz > 0
        flint_free(perm);
        nmod_poly_mat_clear(A);
        if (collect_profile)
            g_nmod_kernel_zls_profile.kernel_total_time += _nmod_kernel_now_seconds() - call_start;

        return nz;
    }

    flint_free(perm);
    nmod_poly_mat_clear(A);
    if (collect_profile)
        g_nmod_kernel_zls_profile.kernel_total_time += _nmod_kernel_now_seconds() - call_start;

    return 0;

}


/**
 * Experimental, should not be really considered
 *
 */

int nmod_poly_mat_approximant_kernel(nmod_poly_mat_t N, slong *degN, const nmod_poly_mat_t A, \
                                 const slong *ishift)
{

    slong i,j,k;

    slong m = A->r;
    slong n = A->c;

    slong deg[n];

    nmod_poly_mat_column_degree(deg, A, NULL);

    slong degmax=deg[0];
    for (i=1; i<n; i++)
        if (degmax < deg[i])
            degmax=deg[i];

    slong min_nm=n;
    if (n >m)
        min_nm=m;

    nmod_poly_mat_t AT;
    nmod_poly_mat_init(AT, n, m, A->modulus);
    nmod_poly_mat_transpose(AT,A);

    nmod_poly_mat_t PT;
    nmod_poly_mat_init(PT, n, n, A->modulus);


    // Approximant PT
    // shift is modified in place
    // --------------------------

    slong shift[n];
    for (i=0; i<n; i++)
    {
        shift[i]=ishift[i];
    }

    nmod_poly_mat_pmbasis(PT, shift, AT, degmax*(min_nm +1)+1);


    // Looking for zero residues and non zero residues
    //  the global residue is reused later
    // -----------------------------------------------

    nmod_poly_mat_t RT;
    nmod_poly_mat_init(RT, n, m, A->modulus);
    nmod_poly_mat_mul(RT,PT,AT);

    slong cdeg[n];
    nmod_poly_mat_row_degree(cdeg, RT, NULL);

    nmod_poly_mat_clear(AT);
    nmod_poly_mat_clear(RT);


    // Zero residue matrix P1
    //-----------------------

    slong n1=0;

    for (j=0; j<n; j++) {
        if (cdeg[j]<0) {
            n1+=1;
        }
    }


    if (n1 > 0) {

        k=0;

        for (j=0; j<n; j++) {

            if (cdeg[j]<0) {
                for (i = 0; i < n; i++)
                    nmod_poly_set(nmod_poly_mat_entry(N, i, k), nmod_poly_mat_entry(PT, j, i));
                k+=1;
            }
        }

        nmod_poly_mat_clear(PT);
        return n1;
    }

    nmod_poly_mat_clear(PT);
    return 0;
}
