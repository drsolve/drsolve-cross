#ifndef FQ_MVPOLY_H
#define FQ_MVPOLY_H

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

#include "fq_unified_interface.h"

/* ============================================================================
 * Data Structures for Multivariate Polynomials over Finite Extension Fields
 * ============================================================================ */

// Monomial structure for finite extension fields
// Represents a single monomial term with separate exponent vectors for 
// variables and parameters, along with a coefficient in F_{p^d}
typedef struct {
    slong *var_exp;      // Exponent vector for variables (including duals)
    slong *par_exp;      // Exponent vector for parameters
    fq_nmod_t coeff;     // Coefficient in F_{p^d}
} fq_monomial_t;

// Multivariate polynomial with parameters over F_{p^d}
// Represents a multivariate polynomial where variables can include dual variables
// for Dixon resultant computation, and coefficients can depend on parameters
typedef struct {
    slong nvars;         // Number of variables (not including duals)
    slong npars;         // Number of parameters
    slong nterms;        // Number of terms
    slong alloc;         // Allocated space
    fq_monomial_t *terms; // Array of terms
    const fq_nmod_ctx_struct *ctx;   // Finite field context (stored as pointer to struct)
} fq_mvpoly_t;

// Comparison structure for sorting degrees
typedef struct {
    slong index;
    slong degree;
} fq_index_degree_pair;

/* ============================================================================
 * Field Equation Reduction
 * ============================================================================ */

/* Global flag: when non-zero, fq_mvpoly_mul reduces each variable modulo
 * x^q - x after every multiplication (working in F_q[vars] / <x^q - x, ...>). */
extern int g_field_equation_reduction;

/* Enable or disable field-equation reduction mode. */
void fq_mvpoly_set_field_equation_reduction(int enable);

/* Reduce all variable exponents of poly in-place using x^q = x. */
void fq_mvpoly_reduce_field_equation(fq_mvpoly_t *poly);

/* ============================================================================
 * Basic Operations
 * ============================================================================ */

// Initialize a multivariate polynomial
void fq_mvpoly_init(fq_mvpoly_t *p, slong nvars, slong npars, const fq_nmod_ctx_t ctx);

// Clear a multivariate polynomial and free memory
void fq_mvpoly_clear(fq_mvpoly_t *p);

// Safe clear function that checks for proper initialization
void fq_mvpoly_clear_safe(fq_mvpoly_t *p);

// Copy one polynomial to another
void fq_mvpoly_copy(fq_mvpoly_t *dest, const fq_mvpoly_t *src);

/* ============================================================================
 * Term Management
 * ============================================================================ */

// Add a term to the polynomial with duplicate checking
void fq_mvpoly_add_term(fq_mvpoly_t *p, const slong *var_exp, const slong *par_exp, const fq_nmod_t coeff);

// Add a term to the polynomial without duplicate checking (faster)
void fq_mvpoly_add_term_fast(fq_mvpoly_t *p, const slong *var_exp, const slong *par_exp, const fq_nmod_t coeff);

/* ============================================================================
 * Display Functions
 * ============================================================================ */

// Print polynomial in standard format
void fq_mvpoly_print(const fq_mvpoly_t *p, const char *name);

// Print polynomial with expanded dual variable notation
void fq_mvpoly_print_expanded(const fq_mvpoly_t *p, const char *name, int use_dual);
void fq_mvpoly_print_with_names(const fq_mvpoly_t *poly, const char *poly_name,
                               char **var_names, char **par_names, 
                               const char *gen_name, int expanded_format);
/* ============================================================================
 * Arithmetic Operations
 * ============================================================================ */

// Multiply two multivariate polynomials
void fq_mvpoly_mul(fq_mvpoly_t *result, const fq_mvpoly_t *a, const fq_mvpoly_t *b);

// Compute polynomial to a power
void fq_mvpoly_pow(fq_mvpoly_t *result, const fq_mvpoly_t *base, slong power);

// Multiply polynomial by scalar
void fq_mvpoly_scalar_mul(fq_mvpoly_t *result, const fq_mvpoly_t *p, const fq_nmod_t scalar);

// Add two polynomials
void fq_mvpoly_add(fq_mvpoly_t *result, const fq_mvpoly_t *a, const fq_mvpoly_t *b);

// Subtract two polynomials
void fq_mvpoly_sub(fq_mvpoly_t *result, const fq_mvpoly_t *a, const fq_mvpoly_t *b);

// Make polynomial monic (leading coefficient = 1)
void fq_mvpoly_make_monic(fq_mvpoly_t *poly);

/* ============================================================================
 * Conversion Functions
 * ============================================================================ */

// Convert fq_mvpoly to FLINT's fq_nmod_mpoly format
void fq_mvpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t mpoly, const fq_mvpoly_t *poly, 
                               fq_nmod_mpoly_ctx_t mpoly_ctx);

// Convert FLINT's fq_nmod_mpoly to fq_mvpoly format
void fq_nmod_mpoly_to_fq_mvpoly(fq_mvpoly_t *poly, const fq_nmod_mpoly_t mpoly,
                                slong nvars, slong npars, 
                                fq_nmod_mpoly_ctx_t mpoly_ctx, const fq_nmod_ctx_t ctx);

/* ============================================================================
 * Kronecker Substitution Helpers
 * ============================================================================ */

// Convert exponent vector to Kronecker index
slong exp_to_kronecker_index(const slong *exp, const slong *degs, slong n);

// Convert Kronecker index to exponent vector
void kronecker_index_to_exp(slong index, slong *exp, const slong *degs, slong n);

// Convert fq_mvpoly to univariate via Kronecker substitution
void fq_mvpoly_to_kronecker_full(fq_nmod_poly_t out, const fq_mvpoly_t *p, 
                                const slong *var_degs, const slong *par_degs);

// Convert univariate polynomial back to fq_mvpoly
void kronecker_to_fq_mvpoly_full(fq_mvpoly_t *out, const fq_nmod_poly_t in, 
                                const slong *var_degs, slong nvars,
                                const slong *par_degs, slong npars, const fq_nmod_ctx_t ctx);

/* ============================================================================
 * Division Operations
 * ============================================================================ */

// Divide by linear factor (x_i - ~x_i) using improved direct method
void divide_by_fq_linear_factor_improved(fq_mvpoly_t *quotient, const fq_mvpoly_t *dividend,
                                        slong var_idx, slong nvars, slong npars);

// Divide by linear factor using FLINT's built-in functions
void divide_by_fq_linear_factor_flint(fq_mvpoly_t *quotient, const fq_mvpoly_t *dividend,
                                     slong var_idx, slong nvars, slong npars);

/* ============================================================================
 * Evaluation Functions
 * ============================================================================ */

// Evaluate polynomial at given parameter values
void evaluate_fq_mvpoly_at_params(fq_nmod_t result, const fq_mvpoly_t *poly, const fq_nmod_t *param_vals);

// Safe evaluation of matrix entry with NULL checking
void evaluate_fq_mvpoly_safe(fq_nmod_t result, fq_mvpoly_t ***matrix, 
                             slong i, slong j, const fq_nmod_t *param_vals,
                             slong npars, const fq_nmod_ctx_t ctx);

/* ============================================================================
 * Matrix Display and Analysis Functions
 * ============================================================================ */

// Print matrix of polynomials with optional details
void print_fq_matrix_mvpoly(fq_mvpoly_t **matrix, slong nrows, slong ncols, 
                            const char *matrix_name, int show_details);

// Analyze and display statistics about a matrix of polynomials
void analyze_fq_matrix_mvpoly(fq_mvpoly_t **matrix, slong nrows, slong ncols, 
                              const char *matrix_name);

/* ============================================================================
 * Comparison Function
 * ============================================================================ */

// Comparison function for sorting degrees
int compare_fq_degrees(const void *a, const void *b);

#endif /* FQ_MVPOLY_H */
