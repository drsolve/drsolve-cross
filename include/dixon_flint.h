#ifndef DIXON_FLINT_H
#define DIXON_FLINT_H

/*
 * dixon_flint.h - Complete Dixon Resultant Implementation for Finite Extension Fields
 *
 * This header provides functions for computing Dixon resultants over finite fields
 * using FLINT library for polynomial arithmetic and matrix operations.
 *
 * Compile with: gcc -O3 -march=native -o dixon_flint dixon_flint.c -lflint -lmpfr -lgmp -lpthread
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include <unistd.h>
#include <ctype.h>

#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>
#include <flint/fq_nmod_mat.h>
#include <flint/fq_nmod_mpoly.h>
#include <flint/nmod_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/profiler.h>

#include "fq_mvpoly.h"
#include "fq_mpoly_mat_det.h"
#include "fq_unified_interface.h"
#include "dixon_interface_flint.h"

// Debug output control - set to 0 to disable all output
#define DEBUG_OUTPUT_D 0

#if DEBUG_OUTPUT_D
    #define DEBUG_PRINT_D(...) printf(__VA_ARGS__)
#else
    #define DEBUG_PRINT_D(...) ((void)0)
#endif


// Method enumeration for determinant computation
typedef enum {
    DET_METHOD_RECURSIVE = 0,     // Recursive expansion method
    DET_METHOD_KRONECKER = 1,     // Kronecker substitution method
    DET_METHOD_INTERPOLATION = 2, // Interpolation method
    DET_METHOD_HUANG = 3          // Huang's interpolation method
} det_method_t;

// Global method selection variable
extern det_method_t dixon_global_method;
extern int g_matrix_transpose_threshold;

// Matrix operations
void build_fq_cancellation_matrix_mvpoly(fq_mvpoly_t ***M, fq_mvpoly_t *polys, 
                                        slong nvars, slong npars);

void perform_fq_matrix_row_operations_mvpoly(fq_mvpoly_t ***new_matrix, fq_mvpoly_t ***original_matrix,
                                           slong nvars, slong npars);

// Degree bound computation
slong compute_fq_dixon_resultant_degree_bound(fq_mvpoly_t *polys, slong npolys, 
                                             slong nvars, slong npars);

// Coefficient matrix determinant computation
void compute_fq_coefficient_matrix_det(fq_mvpoly_t *result, fq_mvpoly_t **coeff_matrix,
                                      slong size, slong npars, const fq_nmod_ctx_t ctx,
                                      det_method_t method, slong res_deg_bound);

// Row basis tracker structure for linear independence checking
typedef struct {
    field_elem_u *reduced_rows;    // Reduced row vectors
    slong *pivot_cols;             // Pivot column positions for each row
    slong *selected_indices;       // Selected original row indices
    slong current_rank;            // Current rank
    slong max_size;
    slong ncols;
    field_ctx_t *ctx;              // Unified field context
    int initialized;               // Initialization flag
    
    // Pre-allocated workspace to avoid repeated allocation
    field_elem_u *work_row;        // Working row
    field_elem_u *temp_vars;       // Temporary variable pool: [factor, temp, pivot_val, neg_temp]
    int workspace_initialized;     // Workspace initialization flag
} unified_row_basis_tracker_t;

// Pivot row finding functions
void find_pivot_rows_nmod_fixed(slong **selected_rows_out, slong *num_selected,
                               const nmod_mat_t mat);

void find_pivot_rows_simple(slong **selected_rows_out, slong *num_selected,
                           const field_elem_u *unified_mat, 
                           slong nrows, slong ncols,
                           field_ctx_t *ctx);

// Maximal rank submatrix finding
void find_fq_optimal_maximal_rank_submatrix(fq_mvpoly_t ***full_matrix, 
                                           slong nrows, slong ncols,
                                           slong **row_indices_out, 
                                           slong **col_indices_out,
                                           slong *num_rows, slong *num_cols,
                                           slong npars);

// Monomial collection structures and functions
typedef struct {
    slong *exp;
    slong idx;
} monom_t;

typedef struct hash_entry {
    slong *exp;
    slong idx;
    struct hash_entry *next;
} hash_entry_t;

// Optimized monomial collection with hash table
void collect_unique_monomials(
    monom_t **x_monoms_out, slong *nx_monoms_out,
    monom_t **dual_monoms_out, slong *ndual_monoms_out,
    const fq_mvpoly_t *dixon_poly, 
    const slong *d0, const slong *d1, slong nvars);

// Lazy matrix entry allocation
fq_mvpoly_t* get_matrix_entry_lazy(fq_mvpoly_t ***matrix, slong i, slong j,
                                  slong npars, const fq_nmod_ctx_t ctx);

// Optimized coefficient matrix filling
void fill_coefficient_matrix_optimized(fq_mvpoly_t ***full_matrix,
                                      monom_t *x_monoms, slong nx_monoms,
                                      monom_t *dual_monoms, slong ndual_monoms,
                                      const fq_mvpoly_t *dixon_poly,
                                      const slong *d0, const slong *d1, 
                                      slong nvars, slong npars);

// Extract coefficient matrix from Dixon polynomial
void extract_fq_coefficient_matrix_from_dixon(fq_mvpoly_t ***coeff_matrix,
                                             slong *row_indices, slong *col_indices,
                                             slong *matrix_size,
                                             const fq_mvpoly_t *dixon_poly,
                                             slong nvars, slong npars);

// Compute determinant of cancellation matrix
void compute_fq_cancel_matrix_det(fq_mvpoly_t *result, fq_mvpoly_t **modified_M_mvpoly,
                                 slong nvars, slong npars, det_method_t method);

slong dixon_matrix_size(slong nvars, slong degree, ulong prime, slong field_degree);

// Main Dixon resultant computation function
void fq_dixon_resultant(fq_mvpoly_t *result, fq_mvpoly_t *polys, 
                       slong nvars, slong npars);

void fq_dixon_resultant_with_names(fq_mvpoly_t *result, fq_mvpoly_t *polys, 
                                  slong nvars, slong npars,
                                  char **var_names, char **par_names, 
                                  const char *gen_name);
#endif
