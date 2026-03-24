#ifndef FQ_NMOD_ROOTS_H
#define FQ_NMOD_ROOTS_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>
#include <flint/fq_nmod_poly_factor.h>
#include <flint/fmpz.h>
#include <gmp.h>
#include <flint/ulong_extras.h>
#include "gf2n_field.h"

// Prime constants for testing
#define PRIME 9223372036854775783ULL
#define SMALL_PRIME 1073741827ULL

// Root storage structure for nmod_poly version
typedef struct {
    mp_limb_t *roots;
    slong *mult;
    slong num;
    slong alloc;
} nmod_roots_struct;
typedef nmod_roots_struct nmod_roots_t[1];

// Root storage structure for fq_nmod_poly version
typedef struct {
    fq_nmod_struct *roots;
    slong *mult;
    slong num;
    slong alloc;
    mp_limb_t p;
    slong d;
} fq_nmod_roots_struct;
typedef fq_nmod_roots_struct fq_nmod_roots_t[1];

// Function declarations

// nmod_roots utility functions
void nmod_roots_init(nmod_roots_t roots);
void nmod_roots_clear(nmod_roots_t roots);
void nmod_roots_fit_length(nmod_roots_t roots, slong len);
void nmod_roots_add(nmod_roots_t roots, mp_limb_t root, slong mult);

// fq_nmod_roots utility functions
void fq_nmod_roots_init(fq_nmod_roots_t roots, const fq_nmod_ctx_t ctx);
void fq_nmod_roots_clear(fq_nmod_roots_t roots, const fq_nmod_ctx_t ctx);
void fq_nmod_roots_fit_length(fq_nmod_roots_t roots, slong len, const fq_nmod_ctx_t ctx);
void fq_nmod_roots_add(fq_nmod_roots_t roots, const fq_nmod_t root, slong mult, const fq_nmod_ctx_t ctx);

// Utility functions
double get_time_roots(void);

// nmod_poly version implementations
void nmod_simple_power_x_mod(nmod_poly_t result, mp_limb_t exp, const nmod_poly_t modulus);
void nmod_extract_linear_factors(nmod_roots_t roots, const nmod_poly_t poly, flint_rand_t state);
slong our_nmod_poly_roots(nmod_roots_t roots, const nmod_poly_t poly, int with_multiplicity);

// fq_nmod_poly version implementations
void fq_nmod_simple_power_x_mod(fq_nmod_poly_t result, mp_limb_t exp, const fq_nmod_poly_t modulus, const fq_nmod_ctx_t ctx);
void fq_nmod_extract_linear_factors(fq_nmod_roots_t roots, const fq_nmod_poly_t poly, flint_rand_t state, const fq_nmod_ctx_t ctx);
slong our_fq_nmod_poly_roots(fq_nmod_roots_t roots, const fq_nmod_poly_t poly, int with_multiplicity, const fq_nmod_ctx_t ctx);

// Polynomial generation functions for testing
void generate_nmod_poly(nmod_poly_t poly, flint_rand_t state, slong degree, mp_limb_t p);
void generate_fq_nmod_poly(fq_nmod_poly_t poly, flint_rand_t state, slong degree, const fq_nmod_ctx_t ctx);

// Benchmark functions
double benchmark_nmod_roots(slong degree, int num_tests);
double benchmark_fq_nmod_roots(slong degree, slong extension, int num_tests);
double benchmark_flint_fq_nmod_factor(slong degree, slong extension, int num_tests);

// Test and verification functions
void test_fq_nmod_correctness(void);
void run_unified_comparison(void);

#endif // FQ_NMOD_ROOTS_H
