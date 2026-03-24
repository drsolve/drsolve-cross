/* unified_mpoly_resultant.h - Resultant computation using generic ring interface with Zech support */

#ifndef UNIFIED_MPOLY_RESULTANT_H
#define UNIFIED_MPOLY_RESULTANT_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mpoly.h>
#include <flint/fq_zech.h>
#include <flint/fq_zech_mpoly.h>
#include <flint/mpoly.h>
#include "unified_mpoly_interface.h"

/* ============================================================================
   MISSING FQ_ZECH_MPOLY FUNCTIONS DECLARATIONS
   ============================================================================ */

/* Implementation of missing fq_zech_mpoly_get_term_exp_ui function */
void fq_zech_mpoly_get_term_exp_ui(ulong * exp, const fq_zech_mpoly_t A, 
                                   slong i, const fq_zech_mpoly_ctx_t ctx);

void fq_zech_mpoly_get_term_coeff_fq_zech(fq_zech_t c, 
                                          const fq_zech_mpoly_t A, 
                                          slong i, 
                                          const fq_zech_mpoly_ctx_t ctx);

/* ============================================================================
   UNIVARIATE POLYNOMIAL TYPE FOR UNIFIED MPOLY
   ============================================================================ */

typedef struct {
    unified_mpoly_t *coeffs;
    fmpz *exps;
    slong alloc;
    slong length;
} unified_mpoly_univar_struct;

typedef unified_mpoly_univar_struct unified_mpoly_univar_t[1];

/* ============================================================================
   MISSING MPOLY_UNIVAR FUNCTIONS DECLARATIONS
   ============================================================================ */

/* Zero out a univariate polynomial */
void mpoly_univar_zero(mpoly_univar_t A, mpoly_void_ring_t R);

/* ============================================================================
   UNIVARIATE POLYNOMIAL OPERATIONS DECLARATIONS
   ============================================================================ */

void unified_mpoly_univar_init(unified_mpoly_univar_t A, unified_mpoly_ctx_t ctx);
void unified_mpoly_univar_clear(unified_mpoly_univar_t A, unified_mpoly_ctx_t ctx);
void unified_mpoly_univar_fit_length(unified_mpoly_univar_t A, slong len, 
                                    unified_mpoly_ctx_t ctx);

/* Debug print for univariate polynomial */
void unified_mpoly_univar_print(const unified_mpoly_univar_t A, const char *var,
                               unified_mpoly_ctx_t ctx);

/* ============================================================================
   HELPER FUNCTIONS DECLARATIONS
   ============================================================================ */

/* Get term exponent - handles Zech */
void unified_mpoly_get_term_exp_ui(ulong *exp, const unified_mpoly_t poly,
                                  slong i, unified_mpoly_ctx_t ctx);

/* Get term coefficient - handles Zech */
void unified_mpoly_get_term_coeff_ui(field_elem_u *coeff, const unified_mpoly_t poly,
                                    slong i, unified_mpoly_ctx_t ctx);

/* Convert multivariate polynomial to univariate */
void unified_mpoly_to_univar(unified_mpoly_univar_t A, const unified_mpoly_t B,
                            slong var, unified_mpoly_ctx_t ctx);

/* Convert unified_mpoly_univar to mpoly_univar for use with generic ring operations */
void unified_mpoly_univar_to_mpoly_univar(mpoly_univar_t A, 
                                         unified_mpoly_univar_t B,
                                         mpoly_void_ring_t R);

/* ============================================================================
   MAIN RESULTANT FUNCTION DECLARATION
   ============================================================================ */

int unified_mpoly_resultant(unified_mpoly_t R, const unified_mpoly_t A,
                           const unified_mpoly_t B, slong var,
                           unified_mpoly_ctx_t ctx);

/* ============================================================================
   RANDOM POLYNOMIAL GENERATION FOR TESTING
   ============================================================================ */

void fq_nmod_mpoly_randtest_custom(fq_nmod_mpoly_t poly, flint_rand_t state,
                                  slong length, slong exp_bound,
                                  const fq_nmod_mpoly_ctx_t ctx);

#endif /* UNIFIED_MPOLY_RESULTANT_H */