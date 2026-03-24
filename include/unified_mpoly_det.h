/* unified_mpoly_det.h - Unified polynomial matrix determinant computation */

#ifndef UNIFIED_MPOLY_DET_H
#define UNIFIED_MPOLY_DET_H

#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "unified_mpoly_interface.h"

/* Configuration */
#define PARALLEL_THRESHOLD 3
#define MAX_PARALLEL_DEPTH 2

/* ============================================================================
   FUNCTION DECLARATIONS
   ============================================================================ */

/* Parallel initialization */
void init_nested_parallelism(void);

/* Core determinant computation functions */
void compute_unified_mpoly_det_recursive(unified_mpoly_t det_result, 
                                        unified_mpoly_t **mpoly_matrix, 
                                        slong size, 
                                        unified_mpoly_ctx_t ctx);

void compute_unified_mpoly_det_parallel(unified_mpoly_t det_result, 
                                       unified_mpoly_t **mpoly_matrix, 
                                       slong size, 
                                       unified_mpoly_ctx_t ctx,
                                       slong depth);

/* Main entry point for determinant computation */
void compute_unified_mpoly_det(unified_mpoly_t det_result,
                              unified_mpoly_t **mpoly_matrix,
                              slong size,
                              unified_mpoly_ctx_t ctx,
                              int use_parallel);

/* Matrix helper functions */
unified_mpoly_t** unified_mpoly_mat_init(slong rows, slong cols, unified_mpoly_ctx_t ctx);
void unified_mpoly_mat_clear(unified_mpoly_t **mat, slong rows, slong cols);
void unified_mpoly_mat_print_pretty(unified_mpoly_t **mat, slong rows, slong cols,
                                   const char **vars);

/* Test function */
void test_unified_mpoly_det(void);

#endif /* UNIFIED_MPOLY_DET_H */