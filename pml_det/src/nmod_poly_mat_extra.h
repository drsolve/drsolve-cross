/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_MAT_EXTRA_H
#define NMOD_POLY_MAT_EXTRA_H

/**
 * \file nmod_poly_mat_extra.h
 * \brief Convenience umbrella header for the drsolve-bundled PML determinant subset
 *
 * This is a trimmed copy of the upstream PML umbrella header. In this bundled
 * subset it only re-exports the internal determinant-related headers that are
 * still present under `pml_det/src/`.
 */

#include <flint/nmod_poly_mat.h>

#include "pml.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_kernel.h"

NMOD_POLY_MAT_INLINE void
apply_perm_to_vector(slong *res, const slong *initial_vect,
                     const slong *perm, slong length)
{
    for (slong i = 0; i < length; i++)
        res[perm[i]] = initial_vect[i];
}

void nmod_poly_mat_det_iter(nmod_poly_t det, nmod_poly_mat_t mat);

/* Kernel-basis / HNF-style determinant recursion completed for FLINT
   from the original PML ntl-extras determinant_generic_knowing_degree idea. */
int nmod_poly_mat_det_hnf(nmod_poly_t det,
                          const nmod_poly_mat_t mat);

int nmod_poly_mat_det_hnf_knowing_degree(nmod_poly_t det,
                                         const nmod_poly_mat_t mat,
                                         slong degree);

#endif /* NMOD_POLY_MAT_EXTRA_H */
