#ifndef DIXON_TEST_H
#define DIXON_TEST_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "dixon_flint.h"
#include "dixon_interface_flint.h"
#include "fq_mvpoly.h"
#include "fq_unified_interface.h"
#include "unified_mpoly_resultant.h"
#include "dixon_with_ideal_reduction.h"
#include "resultant_with_ideal_reduction.h"
#include "dixon_complexity.h"
#include "polynomial_system_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

// ============= Macro Definitions =============

#define UNWRAP(...) __VA_ARGS__
#define COUNT_ARGS(...) (sizeof((const char*[]){__VA_ARGS__}) / sizeof(const char*))

#define DIXON(polys, vars) \
    dixon((const char*[]){UNWRAP polys}, COUNT_ARGS polys, \
          (const char*[]){UNWRAP vars}, COUNT_ARGS vars, ctx)

#define DIXON_WITH_IDEAL(polys, vars) \
    dixon_with_ideal((const char*[]){UNWRAP polys}, COUNT_ARGS polys, \
          (const char*[]){UNWRAP vars}, COUNT_ARGS vars, ideal, ctx)

#define GET_ARG1(a, b) a
#define GET_ARG2(a, b) b
#define APPLY(macro, args) macro args
#define RESULTANT(polys, var) \
    bivariate_resultant(APPLY(GET_ARG1, (UNWRAP polys)), \
                       APPLY(GET_ARG2, (UNWRAP polys)), \
                       UNWRAP(UNWRAP var), ctx)

#define DIXON_COMPLEXITY(polys, vars) \
    dixon_complexity_auto((const char*[]){UNWRAP polys}, COUNT_ARGS polys, \
          (const char*[]){UNWRAP vars}, COUNT_ARGS vars, ctx)

#define MAX_COMPLEXITY(...) \
    extract_max_complexity((const char*[]){__VA_ARGS__}, COUNT_ARGS(__VA_ARGS__))

// ============= Type Definitions =============

typedef struct {
    slong *exponents;      // Combined exponents for all indeterminates
    slong total_degree;    // Total degree of the monomial
} monomial_t;

// ============= Function Declarations =============

// Math Utility Functions
slong binomial_coefficient(slong n, slong k);
slong count_possible_monomials(slong nvars, slong npars, slong max_degree);

// Polynomial Generation Functions
void enumerate_all_monomials(monomial_t **monomials, slong *count, 
                            slong total_indeterminates, slong max_degree);

void generate_random_polynomial(fq_mvpoly_t *poly, slong nvars, slong npars,
                               slong max_degree, double density_ratio,
                               const fq_nmod_ctx_t ctx, flint_rand_t state);

void generate_polynomial_system(fq_mvpoly_t **polys, slong nvars, slong npolys, 
                               slong npars, const slong *degrees,
                               double density_ratio,
                               const fq_nmod_ctx_t ctx, flint_rand_t state);

// Test Functions
void test_dixon_system(const char *test_name, slong nvars, slong npars,
                      ulong p, slong field_degree, const slong *degrees,
                      slong npolys, double density_ratio, flint_rand_t state,
                      int test_case_index);

void test_xhash(void);
int test_dixon_resultant(void);
void test_dixon(int test_mode);
    
#ifdef __cplusplus
}
#endif

#endif /* DIXON_TEST_H */