/**
 * fq_multivariate_interpolation.h - Header file for multivariate polynomial interpolation
 * 
 * Minimally modified fq_nmod interpolation - just add parallelization to core loop
 * Based on the original optimized version
 */

#ifndef FQ_MULTIVARIATE_INTERPOLATION_H
#define FQ_MULTIVARIATE_INTERPOLATION_H

#include <flint/fq_nmod_poly.h>
#include <flint/fq_nmod_mat.h>
#include <flint/fq_nmod.h>
#include <string.h>
#include <time.h>
#include "fq_mat_det.h"
#include "fq_mvpoly.h"

// OpenMP support (optional)
#ifdef _OPENMP
#include <omp.h>
#define USE_OPENMP 1
#else
#define USE_OPENMP 0
#define omp_get_max_threads() 1
#define omp_set_num_threads(n)
#endif

// Global control for parallelization
extern int USE_PARALLEL;

// Debug switch
#define DEBUG_FQ_INTERPOLATION 0

#if DEBUG_FQ_INTERPOLATION
#define FQ_INTERP_PRINT(fmt, ...) printf("[FQ_INTERP] " fmt, ##__VA_ARGS__)
#else
#define FQ_INTERP_PRINT(fmt, ...)
#endif

// Forward declaration
//typedef struct fq_mvpoly_struct fq_mvpoly_t;

// Global timing statistics
typedef struct {
    double lagrange_time;
    double monomial_collection_time;
    double monomial_interpolation_time;
    double result_construction_time;
    double memory_time;
    double other_time;
    slong lagrange_calls;
    slong recursive_calls;
    slong total_monomials;
    slong total_terms_processed;
} InterpolationStats;

// Forward declarations
void fq_mvpoly_init(fq_mvpoly_t *p, slong nvars, slong npars, const fq_nmod_ctx_t ctx);
void fq_mvpoly_clear(fq_mvpoly_t *p);
void fq_mvpoly_add_term(fq_mvpoly_t *p, const slong *var_exp, const slong *par_exp, const fq_nmod_t coeff);

// Parallelization control functions
void fq_interpolation_set_parallel(int use_parallel);
void fq_interpolation_use_half_threads(void);

// Statistics functions
void reset_interpolation_stats(void);
void print_interpolation_stats(void);

// Timer utility
double get_time(void);

// Core optimization: batch matrix evaluation
void fq_evaluate_matrix_at_point_batch(fq_nmod_mat_t result_mat,
                                               fq_mvpoly_t **poly_matrix,
                                               slong size,
                                               const fq_nmod_t *var_vals,
                                               const fq_nmod_t *param_vals,
                                               const fq_nmod_ctx_t ctx);

// Divide-and-conquer approach for building product polynomials
void fq_build_product_tree(fq_nmod_poly_t result, 
                          const fq_nmod_t *nodes, 
                          slong start, slong end,
                          const fq_nmod_ctx_t ctx);

// Optimized Lagrange interpolation using divide-and-conquer
void fq_lagrange_interpolation_optimized(fq_nmod_poly_t result, 
                                             const fq_nmod_t *nodes, 
                                             const fq_nmod_t *values, 
                                             slong k, 
                                             const fq_nmod_ctx_t ctx);

// Modified version of fq_tensor_interpolation_recursive_optimized with fixed timing
void fq_tensor_interpolation_recursive_optimized(fq_mvpoly_t *result,
                                                slong current_dim,
                                                const fq_nmod_t **grids,
                                                const slong *grid_sizes,
                                                const fq_nmod_t *flat_values,
                                                slong *value_offset,
                                                slong total_dims,
                                                const fq_nmod_ctx_t ctx);

void fq_tensor_interpolation_all_vars_optimized(fq_mvpoly_t *result,
                                               const fq_nmod_t **grids,
                                               const fq_nmod_t *values,
                                               const slong *grid_sizes,
                                               slong nvars,
                                               slong npars,
                                               const fq_nmod_ctx_t ctx);

void fq_generate_evaluation_points_optimized(fq_nmod_t **grids, slong *grid_sizes, 
                                            slong total_vars, slong *degrees, 
                                            const fq_nmod_ctx_t ctx);

void fq_compute_det_degree_bounds_optimized(slong *bounds, fq_mvpoly_t **matrix, 
                                           slong size, slong total_vars);

// Main interpolation function with PARALLELIZATION ON POINTS
void fq_compute_det_by_interpolation_optimized(fq_mvpoly_t *result,
                                              fq_mvpoly_t **matrix,
                                              slong size,
                                              slong nvars,
                                              slong npars,
                                              const fq_nmod_ctx_t ctx,
                                              slong *degree_bounds);

// Compatible interface wrapper
void fq_compute_det_by_interpolation(fq_mvpoly_t *result,
                                     fq_mvpoly_t **matrix,
                                     slong size,
                                     slong nvars,
                                     slong npars,
                                     const fq_nmod_ctx_t ctx,
                                     slong uniform_bound);

#endif // FQ_MULTIVARIATE_INTERPOLATION_H