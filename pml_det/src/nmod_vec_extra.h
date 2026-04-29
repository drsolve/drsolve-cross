/*
 * Minimal nmod_vec_extra interface retained for drsolve sparse interpolation.
 *
 * This is a trimmed local header derived from PML's nmod_vec_extra.h and
 * intentionally exposes only the unbalanced modular dot-product routine that
 * drsolve still uses in fq_sparse_interpolation.c.
 */
#ifndef DRSOLVE_NMOD_VEC_EXTRA_H
#define DRSOLVE_NMOD_VEC_EXTRA_H

#include <flint/flint.h>
#include <flint/nmod_types.h>

#ifdef __cplusplus
extern "C" {
#endif

ulong nmod_vec_dot_product_unbalanced(nn_srcptr v1, nn_srcptr v2,
                                      ulong len, ulong max_bits1, ulong max_bits2,
                                      nmod_t mod);

#ifdef __cplusplus
}
#endif

#endif /* DRSOLVE_NMOD_VEC_EXTRA_H */
