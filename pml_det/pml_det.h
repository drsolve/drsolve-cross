/*
 * Minimal public header for the drsolve-bundled PML determinant subset.
 *
 * This header intentionally exposes only the polynomial-matrix determinant
 * entry points used by drsolve.
 */
#ifndef DRSOLVE_PML_DET_H
#define DRSOLVE_PML_DET_H

#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#ifdef __cplusplus
extern "C" {
#endif

void nmod_poly_mat_det_iter(nmod_poly_t det, nmod_poly_mat_t mat);

int nmod_poly_mat_det_hnf(nmod_poly_t det,
                          const nmod_poly_mat_t mat);

int nmod_poly_mat_det_hnf_knowing_degree(nmod_poly_t det,
                                         const nmod_poly_mat_t mat,
                                         slong degree);

#ifdef __cplusplus
}
#endif

#endif /* DRSOLVE_PML_DET_H */
