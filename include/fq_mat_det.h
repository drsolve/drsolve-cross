/* fq_mat_det.h - Unified Matrix Determinant Header */
#ifndef FQ_MAT_DET_H
#define FQ_MAT_DET_H

#include <flint/fq_nmod_mpoly.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mat.h>
#include <flint/fq_nmod_poly.h>
#include "gf2n_field.h"
#include "fq_unified_interface.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
   Unified Matrix Determinant Implementation
   
   Uses the unified field interface for ALL finite fields:
   - Prime fields (FIELD_ID_NMOD)
   - Optimized binary extension fields GF(2^8), GF(2^16), GF(2^32), GF(2^64), GF(2^128)
   - General finite fields (FIELD_ID_FQ)
   
   This provides a single, consistent implementation with optimal performance
   for each field type through the unified interface.
   ============================================================================ */

/* ============================================================================
   UNIFIED MATRIX STRUCTURE - DECLARATIONS
   ============================================================================ */

typedef struct {
    field_elem_u *entries;
    slong r;
    slong c;
    field_ctx_t *ctx;
} unified_mat_struct;
typedef unified_mat_struct unified_mat_t[1];

/* ============================================================================
   UNIFIED MATRIX FUNCTIONS - DECLARATIONS
   ============================================================================ */

/* Matrix initialization and cleanup */
void unified_mat_init(unified_mat_t mat, slong rows, slong cols, field_ctx_t *ctx);
void unified_mat_clear(unified_mat_t mat);

/* Matrix element access */
field_elem_u* unified_mat_entry(unified_mat_t mat, slong i, slong j);
const field_elem_u* unified_mat_entry_const(const unified_mat_t mat, slong i, slong j);

/* Matrix operations */
void unified_mat_swap_rows(unified_mat_t mat, slong r, slong s);

/* ============================================================================
   UNIFIED LU DETERMINANT COMPUTATION - DECLARATION
   
   This single implementation works for all finite fields through the unified
   interface, providing optimal performance for each field type:
   - Direct native operations for prime fields
   - Lookup tables for GF(2^8) and GF(2^16)
   - PCLMUL instructions for GF(2^32), GF(2^64), GF(2^128) when available
   - General FLINT operations for other finite fields
   ============================================================================ */

void unified_mat_det_lu(field_elem_u *det, const unified_mat_t mat);

/* ============================================================================
   MAIN DETERMINANT FUNCTIONS - DECLARATIONS
   ============================================================================ */

/* Main determinant function - uses unified interface for ALL fields */
void fq_nmod_mat_det(fq_nmod_t det, const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx);

/* Alternative: using characteristic polynomial method */
void fq_nmod_mat_det_charpoly(fq_nmod_t det, const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif /* FQ_MAT_DET_H */