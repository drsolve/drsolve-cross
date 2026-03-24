#ifndef RESULTANT_WITH_IDEAL_REDUCTION_H
#define RESULTANT_WITH_IDEAL_REDUCTION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_mpoly.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mpoly.h>
#include "dixon_with_ideal_reduction.h"

// Bivariate resultant with ideal reduction
char* resultant_with_ideal_reduction(const char *poly1_string, const char *poly2_string,
                                   const char *elim_var, 
                                   const fq_nmod_ctx_t ctx,
                                   unified_triangular_ideal_t *ideal);

// String interface for resultant with ideal reduction
char* resultant_with_ideal_reduction_str(const char *poly1_string,
                                       const char *poly2_string, 
                                       const char *elim_var_string,
                                       const char *ideal_gens_string,
                                       const fq_nmod_ctx_t ctx);

// Smart elimination that chooses between Dixon and resultant methods based on input
char* elimination_with_ideal_reduction_str(const char *poly_string,
                                         const char *elim_vars_string,
                                         const char *ideal_gens_string,
                                         const fq_nmod_ctx_t ctx);

// Elimination with ideal - general interface for multiple polynomials
char* elimination_with_ideal(const char **poly_strings,
                           slong num_polys,
                           const char **elim_vars,
                           slong num_elim_vars,
                           const char *ideal_string,
                           const fq_nmod_ctx_t ctx);

// Unified elimination - automatically chooses Dixon or resultant method
char* elimination_with_ideal_reduction(const char **poly_strings, slong num_polys,
                                      const char **elim_vars, slong num_elim_vars,
                                      const fq_nmod_ctx_t ctx,
                                      unified_triangular_ideal_t *ideal);

// Test function for iterative elimination
void test_iterative_elimination_str2(void);

#endif