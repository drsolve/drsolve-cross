/* fq_poly_mat_det.h - Fixed Matrix determinant with compile-time dispatch */
#ifndef FQ_POLY_MAT_DET_H
#define FQ_POLY_MAT_DET_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <flint/flint.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/nmod_vec.h>

/* Conditionally include PML headers if available */
#ifdef HAVE_PML
#include <nmod_poly_mat_utils.h>
#include <nmod_poly_mat_extra.h>
#endif

/* Include specialized implementations */
#include "gf2n_field.h"
#include "fq_unified_interface.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
   GLOBAL VARIABLES FOR PROGRESS TRACKING
   ============================================================================ */

extern int _is_first_step;
extern int _show_progress;
extern int g_show_progress;

/* ============================================================================
   UNIFIED POLYNOMIAL MATRIX OPERATIONS - DECLARATIONS
   ============================================================================ */

/* Swap rows in unified poly mat */
void unified_poly_mat_swap_rows_ptr(unified_poly_mat_struct *mat, slong r, slong s);

/* Rotate rows upward */
void unified_poly_mat_rotate_rows_upward_ptr(unified_poly_mat_struct *mat, slong i, slong j);

/* Permute rows */
void unified_poly_mat_permute_rows_ptr(unified_poly_mat_struct *mat, const slong *perm);

/* Window resize columns */
void unified_poly_mat_window_resize_columns_ptr(unified_poly_mat_struct *mat, slong c);

/* ============================================================================
   PIVOT COMPUTATION - DECLARATIONS
   ============================================================================ */

void unified_poly_vec_pivot_profile(slong *pivind, slong *pivdeg,
                                   unified_poly_struct *vec,
                                   slong len);

/* ============================================================================
   COLLISION SOLVING - DECLARATIONS
   ============================================================================ */

void atomic_solve_pivot_collision_unified_optimized(unified_poly_mat_struct *mat,
                                                   unified_poly_mat_struct *other,
                                                   slong pi1, slong pi2, slong j);

/* ============================================================================
   WEAK POPOV FORM COMPUTATION - DECLARATIONS
   ============================================================================ */

slong unified_poly_mat_weak_popov_iter_submat_rowbyrow_optimized(
    unified_poly_mat_struct *mat,
    unified_poly_mat_struct *tsf,
    slong *det,
    slong *pivind,
    slong *rrp,
    slong rstart,
    slong cstart,
    slong rdim,
    slong cdim,
    slong early_exit_zr);

/* ============================================================================
   HELPER FUNCTIONS - DECLARATIONS
   ============================================================================ */

/* Permutation functions */
typedef struct {
    slong value;
    slong index;
} slong_pair;

int slong_pair_compare(const void *a, const void *b);

void vec_sort_permutation(slong *perm, slong *sorted_vec,
                         const slong *vec, slong n, slong_pair *pair_tmp);

void unified_poly_mat_permute_rows_by_sorting_vec_ptr(unified_poly_mat_struct *mat,
                                                     slong r, slong *vec, slong *perm);

/* Check permutation parity */
int perm_parity(const slong *perm, slong n);

/* ============================================================================
   MAIN ENTRY POINT - DECLARATION
   ============================================================================ */

void fq_nmod_poly_mat_det_iter(fq_nmod_poly_t det,
                              fq_nmod_poly_mat_t mat,
                              const fq_nmod_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif /* FQ_POLY_MAT_DET_H */
