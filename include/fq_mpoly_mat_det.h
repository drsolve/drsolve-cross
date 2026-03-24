/*
 * Optimized polynomial matrix determinant computation for small matrices
 * Supports multiple algorithms including recursive, interpolation, and Kronecker substitution
 * Enhanced with prime field optimization using nmod_mpoly
 */

#ifndef FQ_MPOLY_MAT_DET_H
#define FQ_MPOLY_MAT_DET_H

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#endif
#include <flint/fq_nmod_mpoly.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mat.h>
#include <flint/fq_nmod_poly.h>
#include <flint/nmod_mpoly.h>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include "fq_poly_mat_det.h"
#include "fq_multivariate_interpolation.h"
#include "fq_sparse_interpolation.h"
#include "unified_mpoly_det.h"

// Algorithm selection options
#define DET_ALGORITHM_RECURSIVE       0  // Original recursive expansion algorithm
#define DET_ALGORITHM_INTERPOLATION   1  // Multivariate interpolation algorithm
#define DET_ALGORITHM_KRONECKER       2  // Kronecker substitution to univariate
#define DET_ALGORITHM_POLY_RECURSIVE  3  // Convert to fq_nmod_poly and use recursive
#define DET_ALGORITHM_HUANG           4  // Huang's sparse interpolation (only for SPARSE polynomials over PRIME fields)

// Default algorithm selection 
#ifndef DET_ALGORITHM
#define DET_ALGORITHM DET_ALGORITHM_RECURSIVE
#endif

// Configuration constants
#define PARALLEL_THRESHOLD 3
#define MAX_PARALLEL_DEPTH 2
#define UNIVARIATE_THRESHOLD 3

// Debug control
#define DEBUG_FQ_DET 0

#if DEBUG_FQ_DET
#define DET_PRINT(fmt, ...) printf("[FQ_DET] " fmt, ##__VA_ARGS__)
#else
#define DET_PRINT(fmt, ...)
#endif

// ============= Timing Utilities =============

typedef struct {
    double wall_time;
    double cpu_time;
} timing_info_t;

// Get current CPU time
double get_cpu_time(void);

// Start timing measurement
timing_info_t start_timing(void);

// End timing measurement and return elapsed time
timing_info_t end_timing(timing_info_t start);

// Print timing information with label
void print_timing(const char* label, timing_info_t elapsed);

// ============= Prime Field Detection =============

// Check if context represents a prime field (degree 1)
static inline int is_prime_field(const fq_nmod_ctx_t ctx) {
    return fq_nmod_ctx_degree(ctx) == 1;
}

// ============= Polynomial Operation Optimizations =============

// ============= Conversion Functions for Polynomial Recursive =============

// Convert fq_mvpoly to fq_nmod_poly for a specific variable
void mvpoly_to_fq_nmod_poly(fq_nmod_poly_t poly, const fq_mvpoly_t *mvpoly, 
                           slong var_index, const fq_nmod_ctx_t ctx);

// Convert fq_nmod_poly back to fq_mvpoly
void fq_nmod_poly_to_mvpoly(fq_mvpoly_t *mvpoly, const fq_nmod_poly_t poly,
                           slong var_index, slong nvars, slong npars,
                           const fq_nmod_ctx_t ctx);

// ============= Polynomial Recursive Determinant =============

// Recursive determinant computation using fq_nmod_poly operations
void compute_det_poly_recursive_helper(fq_nmod_poly_t det, 
                                      fq_nmod_poly_t **matrix,
                                      slong size, 
                                      const fq_nmod_ctx_t ctx);

// Main function for polynomial recursive algorithm
void compute_fq_det_poly_recursive(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size);

// ============= Kronecker Substitution Implementation =============

// Compute bounds for Kronecker substitution
void compute_kronecker_bounds(slong *var_bounds, fq_mvpoly_t **matrix, 
                             slong size, slong nvars, slong npars);

// Convert multivariate polynomial to univariate using Kronecker substitution
void mvpoly_to_univariate_kronecker(fq_nmod_poly_t uni_poly,
                                   const fq_mvpoly_t *mv_poly,
                                   const slong *substitution_powers,
                                   const fq_nmod_ctx_t ctx);

// Convert univariate polynomial back to multivariate
void univariate_to_mvpoly_kronecker(fq_mvpoly_t *mv_poly,
                                   const fq_nmod_poly_t uni_poly,
                                   const slong *substitution_powers,
                                   const slong *var_bounds,
                                   slong nvars, slong npars,
                                   const fq_nmod_ctx_t ctx);

// Compute determinant using Kronecker substitution
void compute_fq_det_kronecker(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size);

// ============= Prime Field Conversion Functions =============

// Convert fq_mvpoly to nmod_mpoly for prime fields
void fq_mvpoly_to_nmod_mpoly(nmod_mpoly_t mpoly, const fq_mvpoly_t *poly, 
                            nmod_mpoly_ctx_t mpoly_ctx);

// Convert nmod_mpoly to fq_mvpoly
void nmod_mpoly_to_fq_mvpoly(fq_mvpoly_t *result, const nmod_mpoly_t poly,
                            slong nvars, slong npars,
                            const nmod_mpoly_ctx_t mpoly_ctx,
                            const fq_nmod_ctx_t field_ctx);

// Convert matrix from fq_mvpoly to nmod_mpoly format
void fq_matrix_mvpoly_to_nmod_mpoly(nmod_mpoly_t **mpoly_matrix, 
                                   fq_mvpoly_t **mvpoly_matrix, 
                                   slong size, 
                                   nmod_mpoly_ctx_t mpoly_ctx);

// ============= Prime Field Determinant Computation =============

// Optimized 3x3 determinant for nmod_mpoly
void compute_det_3x3_nmod_optimized(nmod_mpoly_t det, 
                                   nmod_mpoly_t **m,
                                   nmod_mpoly_ctx_t ctx);

// Recursive determinant for nmod_mpoly
void compute_nmod_mpoly_det_recursive(nmod_mpoly_t det_result, 
                                     nmod_mpoly_t **mpoly_matrix, 
                                     slong size, 
                                     nmod_mpoly_ctx_t mpoly_ctx);

// Parallel determinant computation for nmod_mpoly
void compute_nmod_mpoly_det_parallel_optimized(nmod_mpoly_t det_result, 
                                              nmod_mpoly_t **mpoly_matrix, 
                                              slong size, 
                                              nmod_mpoly_ctx_t mpoly_ctx,
                                              slong depth);

// ============= Univariate Optimization =============

// Check if matrix contains only univariate polynomials
int is_univariate_matrix(fq_mvpoly_t **matrix, slong size);

// Compute determinant for univariate polynomial matrices
void compute_fq_det_univariate_optimized(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size);

// ============= Conversion Functions =============

// Convert fq_mvpoly to fq_nmod_mpoly
void fq_mvpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t mpoly, const fq_mvpoly_t *poly, 
                               fq_nmod_mpoly_ctx_t mpoly_ctx);

// Convert fq_nmod_mpoly to fq_mvpoly
void fq_nmod_mpoly_to_fq_mvpoly(fq_mvpoly_t *poly, const fq_nmod_mpoly_t mpoly,
                               slong nvars, slong npars, 
                               fq_nmod_mpoly_ctx_t mpoly_ctx, const fq_nmod_ctx_t ctx);

// Convert matrix from fq_mvpoly to fq_nmod_mpoly format
void fq_matrix_mvpoly_to_mpoly(fq_nmod_mpoly_t **mpoly_matrix, 
                              fq_mvpoly_t **mvpoly_matrix, 
                              slong size, 
                              fq_nmod_mpoly_ctx_t mpoly_ctx);

// ============= Optimized Determinant Computation =============

// Optimized 3x3 determinant for fq_nmod_mpoly
void compute_det_3x3_optimized(fq_nmod_mpoly_t det, 
                              fq_nmod_mpoly_t **m,
                              fq_nmod_mpoly_ctx_t ctx);

// Recursive determinant with optimizations for fq_nmod_mpoly
void compute_fq_nmod_mpoly_det_recursive(fq_nmod_mpoly_t det_result, 
                                        fq_nmod_mpoly_t **mpoly_matrix, 
                                        slong size, 
                                        fq_nmod_mpoly_ctx_t mpoly_ctx);

// Parallel determinant computation for fq_nmod_mpoly
void compute_fq_nmod_mpoly_det_parallel_optimized(fq_nmod_mpoly_t det_result, 
                                                  fq_nmod_mpoly_t **mpoly_matrix, 
                                                  slong size, 
                                                  fq_nmod_mpoly_ctx_t mpoly_ctx,
                                                  slong depth);

// Huang sparse interpolation determinant
void compute_fq_det_huang_interpolation(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size);

// Direct recursive algorithm
void compute_fq_det_unified_interface(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size);

// ============= Main Interface with Algorithm Selection =============

// Main determinant computation function with algorithm selection
void compute_fq_det_recursive_flint(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size);

// Compatibility interface
void compute_fq_det_recursive(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size);

#endif // FQ_MPOLY_MAT_DET_H
