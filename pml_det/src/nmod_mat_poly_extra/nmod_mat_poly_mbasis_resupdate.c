#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nmod_mat_poly.h"

extern int g_dixon_verbose_level;
extern int g_dixon_debug_mode;

/*
 * Experimental scaffold for a FLINT port of PML's mbasis_resupdate path.
 *
 * The intended design is:
 *  - maintain a truncated residual polynomial matrix across iterations,
 *  - update residual coefficients incrementally after each kernel step,
 *  - avoid repeated nmod_mat_poly_mul_coeff(...) in the low-order large-matrix regime.
 *
 * This file currently provides the integration point and controlled dispatch.
 * The actual resupdate loop is not ported yet; for correctness it falls back
 * to the existing mbasis implementation.
 */
void
nmod_mat_poly_mbasis_resupdate(nmod_mat_poly_t appbas,
                               slong * shift,
                               const nmod_mat_poly_t matp,
                               slong order)
{
    if (g_dixon_debug_mode || g_dixon_verbose_level >= 3)
    {
        fprintf(stderr,
                "  [PML mbasis_resupdate] experimental scaffold selected; "
                "falling back to current mbasis implementation.\n");
    }

    nmod_mat_poly_mbasis(appbas, shift, matp, order);
}
