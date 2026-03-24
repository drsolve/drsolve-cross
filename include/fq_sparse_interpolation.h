/**
 * fq_sparse_interpolation.h - Header file for sparse polynomial interpolation
 * 
 * This code is for testing only. Experiments show this method is much less 
 * efficient than usual for dense polynomials from the AO algorithm, and is 
 * applicable only in very sparse cases.
 */

#ifndef FQ_SPARSE_INTERPOLATION_H
#define FQ_SPARSE_INTERPOLATION_H

#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_mpoly.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_factor.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <gmp.h>
#include "gf2n_field.h"

// Debug flag
#define DEBUG 0
#define TIMING 1
#define DETAILED_TIMING 0

// Timing structure
typedef struct {
    clock_t start;
    double elapsed;
} my_timer_t;

// Global random state
extern flint_rand_t global_state;

// Polynomial matrix structure
typedef struct {
    nmod_mpoly_struct** entries;
    slong rows;
    slong cols;
} poly_mat_t;

// Term pair structure
typedef struct {
    mp_limb_t coeff;
    mp_limb_t root;
} term_pair;

typedef struct {
    term_pair* pairs;
    slong length;
} term_list;

// Hash table structure for omega powers
typedef struct {
    mp_limb_t* keys;
    slong* values;
    slong size;
    slong capacity;
} omega_hash_table;

// Timer functions
void timer_start(my_timer_t* t);
void timer_stop(my_timer_t* t);
void timer_print(my_timer_t* t, const char* label);

// Polynomial matrix functions
void poly_mat_init(poly_mat_t* mat, slong rows, slong cols, const nmod_mpoly_ctx_t ctx);
void poly_mat_clear(poly_mat_t* mat, const nmod_mpoly_ctx_t ctx);
void poly_mat_entry_set(poly_mat_t* mat, slong i, slong j, const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx);
void poly_mat_det_2x2(nmod_mpoly_t det, const poly_mat_t* mat, const nmod_mpoly_ctx_t ctx);
void poly_mat_det(nmod_mpoly_t det, const poly_mat_t* mat, const nmod_mpoly_ctx_t ctx);

// Polynomial utility functions
slong poly_max_total_degree(const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx);
void poly_max_degrees_per_var(slong* max_degs, const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx);

// Berlekamp-Massey Algorithm (BM)
void BM(nmod_poly_t C, const mp_limb_t* s, slong N, nmod_t mod);

// Vinvert - compute coefficients from roots using partial fraction decomposition
void Vinvert(mp_limb_t* c1, const nmod_poly_t c, const mp_limb_t* v, 
             const mp_limb_t* a, slong n, nmod_t mod);

// Compare function for sorting by root value
int compare_pairs_by_root(const void* a, const void* b);

// Use FLINT's built-in functions for root finding
void find_roots_flint(nmod_poly_t poly, mp_limb_t* roots, 
                     slong* num_roots, nmod_t mod);

// MC - Monomials and Coefficients (Algorithm 3.1)
term_list MC(const mp_limb_t* a, slong N, nmod_t mod);

// Compute discrete logarithm with bounded search
slong mylog(const mp_limb_t omega, const mp_limb_t B, slong d, nmod_t mod);

// Myeval with detailed timing
void Myeval(nmod_mat_t c, const nmod_mpoly_t f, slong n, 
            const mp_limb_t* alpha, nmod_t mod, const nmod_mpoly_ctx_t mctx);

// TotalMyeval with detailed timing
void TotalMyeval(nmod_mat_t M, const nmod_mpoly_t f, slong n,
                 const mp_limb_t* alpha, const mp_limb_t omega, 
                 nmod_t mod, const nmod_mpoly_ctx_t mctx);

// Diversification transformation
void Mydiver(nmod_mpoly_t g, const nmod_mpoly_t f, slong n,
             const mp_limb_t* zeta, nmod_t mod, 
             const nmod_mpoly_ctx_t mctx);

// Hash table functions for omega powers
void omega_hash_init(omega_hash_table* ht, slong max_power);
void omega_hash_clear(omega_hash_table* ht);
void omega_hash_insert(omega_hash_table* ht, const mp_limb_t key, slong value);
slong omega_hash_find(omega_hash_table* ht, const mp_limb_t key);

// MBOT function
void MBOT(nmod_mpoly_t result, const nmod_mat_t M, slong n, slong T,
          const mp_limb_t omega, const mp_limb_t* zeta, nmod_t mod,
          slong* max_degs_per_var, const nmod_mpoly_ctx_t mctx);

// Generate random polynomial with specified parameters
void myrandpoly(nmod_mpoly_t f, slong n, slong T, slong D, 
                nmod_t mod, const nmod_mpoly_ctx_t mctx);

// Check if field size constraint is satisfied
int check_field_constraint(mp_limb_t p, slong n, slong T, slong d);

// Find primitive root of prime p
void find_primitive_root(mp_limb_t* omega, mp_limb_t p);

// Main interpolation algorithm with constraint checking
void ComputePolyMatrixDet(nmod_mpoly_t det_poly, nmod_mpoly_t actual_det,
                         const poly_mat_t* A, slong n, mp_limb_t p,
                         const nmod_mpoly_ctx_t mctx);

// Test function for random polynomial
void test_random_polynomial(void);

// Test function
int huang_test(void);

#endif
