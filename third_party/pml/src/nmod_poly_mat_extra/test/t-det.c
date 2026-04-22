/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_poly_mat.h>
#include <flint/nmod_types.h>
#include <flint/profiler.h>
#include <flint/test_helpers.h>

#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_extra.h"  // det_iter currently in _extra.h

// test one given input
int core_test_determinant(const nmod_poly_mat_t mat)
{
    nmod_poly_t det;
    nmod_poly_init(det, mat->modulus);
    nmod_poly_t det_correct;
    nmod_poly_init(det_correct, mat->modulus);
    // init copy of mat
    nmod_poly_mat_t copy_mat;
    nmod_poly_mat_init(copy_mat, mat->r, mat->c, mat->modulus);

    // verification of determinant
    int verif_det;

    nmod_poly_mat_det(det_correct, mat);

    { // Mulders and Storjohann's algorithm, row by row variant
        nmod_poly_mat_set(copy_mat, mat);
        nmod_poly_mat_det_iter(det, copy_mat);
        verif_det = nmod_poly_equal(det_correct, det);
    }

    nmod_poly_mat_clear(copy_mat);
    nmod_poly_clear(det_correct);
    nmod_poly_clear(det);

    return verif_det;
}

int core_test_generic_determinant_parallel(const nmod_poly_mat_t mat)
{
    int verif_det = 1;
    nmod_poly_t det_correct, det_seq, det_par, det_legacy;
    nmod_poly_init(det_correct, mat->modulus);
    nmod_poly_init(det_seq, mat->modulus);
    nmod_poly_init(det_par, mat->modulus);
    nmod_poly_init(det_legacy, mat->modulus);

    nmod_poly_mat_det(det_correct, mat);

    {
        const int ok_seq = nmod_poly_mat_det_generic_with_opts(det_seq, mat, 0);
        const int ok_par = nmod_poly_mat_det_generic_with_opts(det_par, mat, 1);
        const int ok_legacy = nmod_poly_mat_det_generic(det_legacy, mat);

        verif_det = (ok_seq == ok_par) && (ok_seq == ok_legacy);

        if (verif_det && ok_seq)
        {
            verif_det = nmod_poly_equal(det_seq, det_par)
                     && nmod_poly_equal(det_seq, det_legacy);

            if (verif_det && !nmod_poly_is_zero(det_correct))
            {
                const slong degree = nmod_poly_degree(det_correct);
                nmod_poly_t det_seq_deg, det_par_deg, det_legacy_deg;
                nmod_poly_init(det_seq_deg, mat->modulus);
                nmod_poly_init(det_par_deg, mat->modulus);
                nmod_poly_init(det_legacy_deg, mat->modulus);

                verif_det =
                    nmod_poly_mat_det_generic_knowing_degree_with_opts(det_seq_deg, mat, degree, 0)
                    && nmod_poly_mat_det_generic_knowing_degree_with_opts(det_par_deg, mat, degree, 1)
                    && nmod_poly_mat_det_generic_knowing_degree(det_legacy_deg, mat, degree)
                    && nmod_poly_equal(det_seq, det_seq_deg)
                    && nmod_poly_equal(det_seq, det_par_deg)
                    && nmod_poly_equal(det_seq, det_legacy_deg);

                nmod_poly_clear(det_legacy_deg);
                nmod_poly_clear(det_par_deg);
                nmod_poly_clear(det_seq_deg);
            }
        }
    }

    nmod_poly_clear(det_legacy);
    nmod_poly_clear(det_par);
    nmod_poly_clear(det_seq);
    nmod_poly_clear(det_correct);

    return verif_det;
}

TEST_FUNCTION_START(nmod_poly_mat_det, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong nbits = 2 + n_randint(state, 62);
        ulong rdim = 1 + n_randint(state, 20);
        ulong order = n_randint(state, 40);

        ulong prime = n_randprime(state, nbits, 1);

        nmod_poly_mat_t mat;
        nmod_poly_mat_init(mat, rdim, rdim, prime);

        if (i < 30)
            nmod_poly_mat_randtest_sparse(mat, state, order, 0.2);
        else if (i < 60)
            nmod_poly_mat_rand(mat, state, order);
        else
            nmod_poly_mat_randtest(mat, state, order);

        result = core_test_determinant(mat);

        nmod_poly_mat_clear(mat);

        if (!result)
            TEST_FUNCTION_FAIL("dim = %wu, degree = %wu, p = %wu\n",
                    rdim, order, prime);
    }

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        ulong nbits = 2 + n_randint(state, 20);
        ulong rdim = 12 + n_randint(state, 8);
        ulong order = 1 + n_randint(state, 8);

        ulong prime = n_randprime(state, nbits, 1);

        nmod_poly_mat_t mat;
        nmod_poly_mat_init(mat, rdim, rdim, prime);

        if (i % 2 == 0)
            nmod_poly_mat_randtest_sparse(mat, state, order, 0.3);
        else
            nmod_poly_mat_randtest(mat, state, order);

        result = core_test_generic_determinant_parallel(mat);

        nmod_poly_mat_clear(mat);

        if (!result)
            TEST_FUNCTION_FAIL("generic det mismatch: dim = %wu, degree = %wu, p = %wu\n",
                    rdim, order, prime);
    }

    TEST_FUNCTION_END(state);
}