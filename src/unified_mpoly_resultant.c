/* unified_mpoly_resultant.c - Implementation of resultant computation using generic ring interface with Zech support */

#include "unified_mpoly_resultant.h"
extern int g_field_equation_reduction;

/* ============================================================================
   MISSING FQ_ZECH_MPOLY FUNCTIONS IMPLEMENTATION
   ============================================================================ */

/* Implementation of missing fq_zech_mpoly_get_term_exp_ui function */
void fq_zech_mpoly_get_term_exp_ui(ulong * exp, const fq_zech_mpoly_t A, 
                                   slong i, const fq_zech_mpoly_ctx_t ctx)
{
    slong N;

    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "index out of range in fq_zech_mpoly_get_term_exp_ui");
    }

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    mpoly_get_monomial_ui(exp, A->exps + N*i, A->bits, ctx->minfo);
}

void fq_zech_mpoly_get_term_coeff_fq_zech(fq_zech_t c, 
                                          const fq_zech_mpoly_t A, 
                                          slong i, 
                                          const fq_zech_mpoly_ctx_t ctx)
{
    if (i >= (slong) A->length)
    {
        flint_throw(FLINT_ERROR, "fq_zech_mpoly_get_term_coeff_fq_zech: index out of range");
    }

    fq_zech_set(c, A->coeffs + i, ctx->fqctx);
}

/* ============================================================================
   MISSING MPOLY_UNIVAR FUNCTIONS IMPLEMENTATION
   ============================================================================ */

/* Zero out a univariate polynomial */
void mpoly_univar_zero(mpoly_univar_t A, mpoly_void_ring_t R) {
    /* Clear existing coefficients */
    for (slong i = 0; i < A->length; i++) {
        void *coeff_ptr = (char *)A->coeffs + i * R->elem_size;
        R->zero(coeff_ptr, R->ctx);
    }
    A->length = 0;
}

/* ============================================================================
   UNIVARIATE POLYNOMIAL OPERATIONS IMPLEMENTATION
   ============================================================================ */

void unified_mpoly_univar_init(unified_mpoly_univar_t A, unified_mpoly_ctx_t ctx) {
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
}

void unified_mpoly_univar_clear(unified_mpoly_univar_t A, unified_mpoly_ctx_t ctx) {
    slong i;
    for (i = 0; i < A->alloc; i++) {
        if (A->coeffs && A->coeffs[i]) {
            unified_mpoly_clear(A->coeffs[i]);
            A->coeffs[i] = NULL;
        }
        if (A->exps) {
            fmpz_clear(A->exps + i);
        }
    }
    if (A->coeffs) {
        free(A->coeffs);
        A->coeffs = NULL;
    }
    if (A->exps) {
        free(A->exps);
        A->exps = NULL;
    }
    A->length = 0;
    A->alloc = 0;
}

void unified_mpoly_univar_fit_length(unified_mpoly_univar_t A, slong len, 
                                    unified_mpoly_ctx_t ctx) {
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(len, 2*A->alloc);
    
    if (len <= old_alloc) {
        return;
    }
    
    if (new_alloc < 1) new_alloc = 1;
    
    A->exps = (fmpz *)realloc(A->exps, new_alloc * sizeof(fmpz));
    A->coeffs = (unified_mpoly_t *)realloc(A->coeffs, new_alloc * sizeof(unified_mpoly_t));
    
    for (i = old_alloc; i < new_alloc; i++) {
        fmpz_init(A->exps + i);
        A->coeffs[i] = unified_mpoly_init(ctx);
    }
    
    A->alloc = new_alloc;
}

/* Debug print for univariate polynomial */
void unified_mpoly_univar_print(const unified_mpoly_univar_t A, const char *var,
                               unified_mpoly_ctx_t ctx) {
    slong i;
    printf("Univariate polynomial in %s:\n", var);
    for (i = 0; i < A->length; i++) {
        printf("  Coefficient of %s^", var);
        fmpz_print(A->exps + i);
        printf(": ");
        unified_mpoly_print_pretty(A->coeffs[i], (const char *[]){"x", "y", "z"});
        printf("\n");
    }
}

/* ============================================================================
   HELPER FUNCTIONS WITH ZECH SUPPORT IMPLEMENTATION
   ============================================================================ */

/* Get term exponent - handles Zech */
void unified_mpoly_get_term_exp_ui(ulong *exp, const unified_mpoly_t poly,
                                  slong i, unified_mpoly_ctx_t ctx) {
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_get_term_exp_ui(exp, GET_NMOD_POLY(poly), i, GET_NMOD_CTX(ctx));
            break;
            
        case FIELD_ID_FQ_ZECH:
            fq_zech_mpoly_get_term_exp_ui(exp, GET_ZECH_POLY(poly), i, GET_ZECH_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_get_term_exp_ui(exp, GET_FQ_POLY(poly), i, GET_FQ_CTX(ctx));
            break;
    }
}

/* Get term coefficient - handles Zech */
void unified_mpoly_get_term_coeff_ui(field_elem_u *coeff, const unified_mpoly_t poly,
                                    slong i, unified_mpoly_ctx_t ctx) {
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            coeff->nmod = nmod_mpoly_get_term_coeff_ui(GET_NMOD_POLY(poly), i, 
                                                       GET_NMOD_CTX(ctx));
            break;
            
        case FIELD_ID_FQ_ZECH:
            {
                /* Now we can properly get the coefficient using our new function */
                fq_zech_init(&coeff->fq_zech, ctx->field_ctx->ctx.zech_ctx);
                fq_zech_mpoly_get_term_coeff_fq_zech(&coeff->fq_zech, 
                                                    GET_ZECH_POLY(poly), i, 
                                                    GET_ZECH_CTX(ctx));
            }
            break;
            
        default:
            {
                fq_nmod_t temp;
                fq_nmod_init(temp, ctx->field_ctx->ctx.fq_ctx);
                fq_nmod_mpoly_get_term_coeff_fq_nmod(temp, GET_FQ_POLY(poly), i, 
                                                    GET_FQ_CTX(ctx));
                fq_nmod_to_field_elem(coeff, temp, ctx->field_ctx);
                fq_nmod_clear(temp, ctx->field_ctx->ctx.fq_ctx);
            }
            break;
    }
}

/* Convert multivariate polynomial to univariate */
void unified_mpoly_to_univar(unified_mpoly_univar_t A, const unified_mpoly_t B,
                            slong var, unified_mpoly_ctx_t ctx) {
    slong i, j;
    slong nvars = ctx->nvars;
    slong Blen = unified_mpoly_length(B);
    ulong *exp = (ulong *)calloc(nvars, sizeof(ulong));
    ulong *exp_copy = (ulong *)calloc(nvars, sizeof(ulong));
    field_elem_u coeff;
    
    /* Clear A */
    A->length = 0;
    
    if (Blen == 0) {
        free(exp);
        free(exp_copy);
        return;
    }
    
    void *ctx_ptr = NULL;
    switch (ctx->field_ctx->field_id) {
        case FIELD_ID_NMOD:
            ctx_ptr = (void*)&ctx->field_ctx->ctx.nmod_ctx;
            break;
        case FIELD_ID_FQ_ZECH:
            ctx_ptr = (void*)ctx->field_ctx->ctx.zech_ctx;
            break;
        default:
            ctx_ptr = (void*)ctx->field_ctx->ctx.fq_ctx;
            break;
    }
    field_init_elem(&coeff, ctx->field_ctx->field_id, ctx_ptr);
    
    /* Process each term of B */
    for (i = 0; i < Blen; i++) {
        /* Get exponent and coefficient for term i */
        unified_mpoly_get_term_exp_ui(exp, B, i, ctx);
        unified_mpoly_get_term_coeff_ui(&coeff, B, i, ctx);
        
        /* Find or create entry for this power of var */
        ulong var_exp = exp[var];
        slong pos = -1;
        
        for (j = 0; j < A->length; j++) {
            if (fmpz_equal_ui(A->exps + j, var_exp)) {
                pos = j;
                break;
            }
        }
        
        if (pos == -1) {
            /* New power - add it */
            unified_mpoly_univar_fit_length(A, A->length + 1, ctx);
            pos = A->length;
            A->length++;
            fmpz_set_ui(A->exps + pos, var_exp);
            unified_mpoly_zero(A->coeffs[pos]);
        }
        
        /* Add monomial to coefficient */
        memcpy(exp_copy, exp, nvars * sizeof(ulong));
        exp_copy[var] = 0;  /* Remove power of var */
        
        unified_mpoly_t temp = unified_mpoly_init(ctx);
        unified_mpoly_set_coeff_ui(temp, &coeff, exp_copy);
        unified_mpoly_add(A->coeffs[pos], A->coeffs[pos], temp);
        unified_mpoly_clear(temp);
    }
    
    field_clear_elem(&coeff, ctx->field_ctx->field_id, ctx_ptr);
    free(exp);
    free(exp_copy);
    
    /* Sort by decreasing powers */
    for (i = 0; i < A->length - 1; i++) {
        for (j = i + 1; j < A->length; j++) {
            if (fmpz_cmp(A->exps + i, A->exps + j) < 0) {
                /* Swap */
                fmpz_swap(A->exps + i, A->exps + j);
                unified_mpoly_t temp = A->coeffs[i];
                A->coeffs[i] = A->coeffs[j];
                A->coeffs[j] = temp;
            }
        }
    }
}

/* Convert unified_mpoly_univar to mpoly_univar for use with generic ring operations */
void unified_mpoly_univar_to_mpoly_univar(mpoly_univar_t A, 
                                         unified_mpoly_univar_t B,
                                         mpoly_void_ring_t R) {
    slong i;
    
    /* Clear and initialize A */
    mpoly_univar_zero(A, R);

    mpoly_univar_fit_length(A, B->length, R);
    A->length = B->length;
    
    for (i = 0; i < B->length; i++) {
        fmpz_set(A->exps + i, B->exps + i);
        
        void *coeff_ptr = (char *)A->coeffs + i * R->elem_size;
        R->set(coeff_ptr, &(B->coeffs[i]), R->ctx);
    }
}
/* ============================================================================
   MAIN RESULTANT FUNCTION USING GENERIC RING INTERFACE
   ============================================================================ */

int unified_mpoly_resultant(unified_mpoly_t R, const unified_mpoly_t A,
                           const unified_mpoly_t B, slong var,
                           unified_mpoly_ctx_t ctx) {
    int success;
    int saved_field_eq_mode = g_field_equation_reduction;
    g_field_equation_reduction = 0;
    
    /* Generic implementation using ring interface */
    unified_mpoly_univar_t Ax, Bx;
    mpoly_univar_t Ax_generic, Bx_generic;
    mpoly_void_ring_t ring;
    
    /* Initialize our ring structure */
    unified_mpoly_ring_init(ring, ctx);
    
    /* Initialize univariate polynomials */
    unified_mpoly_univar_init(Ax, ctx);
    unified_mpoly_univar_init(Bx, ctx);
    //printf("Convert A to univariate...\n");
    /* Convert to univariate */
    unified_mpoly_to_univar(Ax, A, var, ctx);

    //printf("Convert B to univariate...\n");
    unified_mpoly_to_univar(Bx, B, var, ctx);
    
    /* Debug print (uncommented for debugging)
    printf("\nDebug: A as univariate:\n");
    unified_mpoly_univar_print(Ax, "X", ctx);
    printf("\nDebug: B as univariate:\n");
    unified_mpoly_univar_print(Bx, "X", ctx);
    */
    
    /* Initialize generic univariate polynomials */
    mpoly_univar_init(Ax_generic, ring);
    mpoly_univar_init(Bx_generic, ring);
    
    /* Convert to generic format - this creates new polynomials owned by Ax_generic/Bx_generic */
    //printf("Convert A to generic format...\n");
    unified_mpoly_univar_to_mpoly_univar(Ax_generic, Ax, ring);
    //printf("Convert B to generic format...\n");
    unified_mpoly_univar_to_mpoly_univar(Bx_generic, Bx, ring);
    
    /* Initialize result properly */
    unified_mpoly_zero(R);
    unified_mpoly_t *R_ptr = &R;
    
    /* Compute resultant */
    success = mpoly_univar_resultant(R_ptr, Ax_generic, Bx_generic, ring);
    //printf("end compute...\n");
    
    /* Cleanup - mpoly_univar_clear will clean up the coefficients it owns */
    mpoly_univar_clear(Ax_generic, ring);
    mpoly_univar_clear(Bx_generic, ring);
    
    /* Clear our univariate structures */
    unified_mpoly_univar_clear(Ax, ctx);
    unified_mpoly_univar_clear(Bx, ctx);
    g_field_equation_reduction = saved_field_eq_mode;

    return success;
}

/* ============================================================================
   RANDOM POLYNOMIAL GENERATION FOR TESTING
   ============================================================================ */

void fq_nmod_mpoly_randtest_custom(fq_nmod_mpoly_t poly, flint_rand_t state,
                                  slong length, slong exp_bound,
                                  const fq_nmod_mpoly_ctx_t ctx) {
    fq_nmod_mpoly_zero(poly, ctx);
    
    slong nvars = ctx->minfo->nvars;
    ulong *exp = (ulong *)calloc(nvars, sizeof(ulong));
    fq_nmod_t coeff;
    fq_nmod_init(coeff, ctx->fqctx);
    
    for (slong i = 0; i < length; i++) {
        /* Random exponents */
        for (slong j = 0; j < nvars; j++) {
            exp[j] = n_randint(state, exp_bound);
        }
        
        /* Random non-zero coefficient */
        fq_nmod_randtest_not_zero(coeff, state, ctx->fqctx);
        fq_nmod_mpoly_set_coeff_fq_nmod_ui(poly, coeff, exp, ctx);
    }
    
    free(exp);
    fq_nmod_clear(coeff, ctx->fqctx);
}
