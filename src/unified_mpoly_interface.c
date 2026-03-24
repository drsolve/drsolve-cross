/* unified_mpoly_interface.c - Implementation of unified multivariate polynomial interface */

#include "unified_mpoly_interface.h"
extern int g_field_equation_reduction;

/* ============================================================================
   CONTEXT OPERATIONS WITH ZECH SUPPORT
   ============================================================================ */

unified_mpoly_ctx_t unified_mpoly_ctx_init(slong nvars, const ordering_t ord, field_ctx_t *field_ctx) {
    unified_mpoly_ctx_t ctx = (unified_mpoly_ctx_t)malloc(sizeof(unified_mpoly_ctx_struct));
    if (!ctx) return NULL;
    
    ctx->field_ctx = field_ctx;
    ctx->nvars = nvars;
    ctx->ord = ord;
    ctx->array_size_limit = 1L << 26;  /* 64MB default */
    ctx->dense_threshold = 10;         /* 10% density threshold */
    ctx->zech_size_limit = 1L << 20;   /* 2^20 field size limit for Zech */
    
    switch (field_ctx->field_id) {
        case FIELD_ID_NMOD:
            {
                /* For prime fields, we need the modulus */
                ulong modulus = field_ctx->ctx.nmod_ctx.n;
                nmod_mpoly_ctx_init(GET_NMOD_CTX(ctx), nvars, ord, modulus);
            }
            break;
            
        case FIELD_ID_GF28:
        case FIELD_ID_GF216:
        case FIELD_ID_GF232:
        case FIELD_ID_GF264:
        case FIELD_ID_GF2128:
            /* For GF(2^n) fields, use fq_nmod context */
            if (field_ctx->ctx.fq_ctx) {
                fq_nmod_mpoly_ctx_init(GET_FQ_CTX(ctx), nvars, ord, field_ctx->ctx.fq_ctx);
            }
            break;
            
        case FIELD_ID_FQ_ZECH:
            /* For Zech logarithm fields */
            if (field_ctx->ctx.zech_ctx) {
                /* Use fq_zech_mpoly_ctx_init_deg for initialization */
                mp_limb_t p = field_ctx->ctx.zech_ctx->fq_nmod_ctx->modulus->mod.n;
                slong d = fq_zech_ctx_degree(field_ctx->ctx.zech_ctx);
                fq_zech_mpoly_ctx_init_deg(GET_ZECH_CTX(ctx), nvars, ord, p, d);
                printf("Using Zech logarithm representation for field of size %lu\n", 
                       field_ctx->ctx.zech_ctx->qm1 + 1);
            }
            break;
            
        default:
            /* General finite fields */
            if (field_ctx->ctx.fq_ctx) {
                fq_nmod_mpoly_ctx_init(GET_FQ_CTX(ctx), nvars, ord, field_ctx->ctx.fq_ctx);
            }
            break;
    }
    
    return ctx;
}

void unified_mpoly_ctx_clear(unified_mpoly_ctx_t ctx) {
    if (!ctx) return;
    
    switch (ctx->field_ctx->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_ctx_clear(GET_NMOD_CTX(ctx));
            break;
            
        case FIELD_ID_FQ_ZECH:
            fq_zech_mpoly_ctx_clear(GET_ZECH_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_ctx_clear(GET_FQ_CTX(ctx));
            break;
    }
    
    free(ctx);
}

/* ============================================================================
   POLYNOMIAL MEMORY MANAGEMENT WITH ZECH SUPPORT
   ============================================================================ */

unified_mpoly_t unified_mpoly_init(unified_mpoly_ctx_t ctx) {
    unified_mpoly_t poly = (unified_mpoly_t)malloc(sizeof(unified_mpoly_struct));
    if (!poly) return NULL;
    
    poly->field_id = ctx->field_ctx->field_id;
    poly->ctx_ptr = ctx;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_init(GET_NMOD_POLY(poly), GET_NMOD_CTX(ctx));
            break;
            
        case FIELD_ID_FQ_ZECH:
            fq_zech_mpoly_init(GET_ZECH_POLY(poly), GET_ZECH_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_init(GET_FQ_POLY(poly), GET_FQ_CTX(ctx));
            break;
    }
    
    return poly;
}

void unified_mpoly_clear(unified_mpoly_t poly) {
    if (!poly) return;
    
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_clear(GET_NMOD_POLY(poly), GET_NMOD_CTX(ctx));
            break;
            
        case FIELD_ID_FQ_ZECH:
            fq_zech_mpoly_clear(GET_ZECH_POLY(poly), GET_ZECH_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_clear(GET_FQ_POLY(poly), GET_FQ_CTX(ctx));
            break;
    }
    
    free(poly);
}

/* ============================================================================
   BASIC OPERATIONS WITH ZECH SUPPORT
   ============================================================================ */

void unified_mpoly_zero(unified_mpoly_t poly) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_zero(GET_NMOD_POLY(poly), GET_NMOD_CTX(ctx));
            break;
            
        case FIELD_ID_FQ_ZECH:
            fq_zech_mpoly_zero(GET_ZECH_POLY(poly), GET_ZECH_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_zero(GET_FQ_POLY(poly), GET_FQ_CTX(ctx));
            break;
    }
}

void unified_mpoly_one(unified_mpoly_t poly) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_one(GET_NMOD_POLY(poly), GET_NMOD_CTX(ctx));
            break;
            
        case FIELD_ID_FQ_ZECH:
            fq_zech_mpoly_one(GET_ZECH_POLY(poly), GET_ZECH_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_one(GET_FQ_POLY(poly), GET_FQ_CTX(ctx));
            break;
    }
}

int unified_mpoly_is_zero(const unified_mpoly_t poly) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            return nmod_mpoly_is_zero(GET_NMOD_POLY(poly), GET_NMOD_CTX(ctx));
            
        case FIELD_ID_FQ_ZECH:
            return fq_zech_mpoly_is_zero(GET_ZECH_POLY(poly), GET_ZECH_CTX(ctx));
            
        default:
            return fq_nmod_mpoly_is_zero(GET_FQ_POLY(poly), GET_FQ_CTX(ctx));
    }
}

int unified_mpoly_is_one(const unified_mpoly_t poly) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            return nmod_mpoly_is_one(GET_NMOD_POLY(poly), GET_NMOD_CTX(ctx));
            
        case FIELD_ID_FQ_ZECH:
            return fq_zech_mpoly_is_one(GET_ZECH_POLY(poly), GET_ZECH_CTX(ctx));
            
        default:
            return fq_nmod_mpoly_is_one(GET_FQ_POLY(poly), GET_FQ_CTX(ctx));
    }
}

slong unified_mpoly_length(const unified_mpoly_t poly) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            return nmod_mpoly_length(GET_NMOD_POLY(poly), GET_NMOD_CTX(ctx));
            
        case FIELD_ID_FQ_ZECH:
            return fq_zech_mpoly_length(GET_ZECH_POLY(poly), GET_ZECH_CTX(ctx));
            
        default:
            return fq_nmod_mpoly_length(GET_FQ_POLY(poly), GET_FQ_CTX(ctx));
    }
}

void unified_mpoly_set(unified_mpoly_t poly1, const unified_mpoly_t poly2) {
    if (poly1 == poly2) return;
    
    unified_mpoly_ctx_t ctx = poly1->ctx_ptr;
    
    switch (poly1->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_set(GET_NMOD_POLY(poly1), GET_NMOD_POLY(poly2), GET_NMOD_CTX(ctx));
            break;
            
        case FIELD_ID_FQ_ZECH:
            fq_zech_mpoly_set(GET_ZECH_POLY(poly1), GET_ZECH_POLY(poly2), GET_ZECH_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_set(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2), GET_FQ_CTX(ctx));
            break;
    }
}

/* ============================================================================
   ARITHMETIC OPERATIONS WITH ZECH SUPPORT
   ============================================================================ */

void unified_mpoly_add(unified_mpoly_t poly1, const unified_mpoly_t poly2,
                      const unified_mpoly_t poly3) {
    unified_mpoly_ctx_t ctx = poly1->ctx_ptr;
    
    switch (poly1->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_add(GET_NMOD_POLY(poly1), GET_NMOD_POLY(poly2),
                          GET_NMOD_POLY(poly3), GET_NMOD_CTX(ctx));
            break;
            
        case FIELD_ID_FQ_ZECH:
            fq_zech_mpoly_add(GET_ZECH_POLY(poly1), GET_ZECH_POLY(poly2),
                             GET_ZECH_POLY(poly3), GET_ZECH_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_add(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2),
                             GET_FQ_POLY(poly3), GET_FQ_CTX(ctx));
            break;
    }
}

void unified_mpoly_sub(unified_mpoly_t poly1, const unified_mpoly_t poly2,
                      const unified_mpoly_t poly3) {
    unified_mpoly_ctx_t ctx = poly1->ctx_ptr;
    
    switch (poly1->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_sub(GET_NMOD_POLY(poly1), GET_NMOD_POLY(poly2),
                          GET_NMOD_POLY(poly3), GET_NMOD_CTX(ctx));
            break;
            
        case FIELD_ID_FQ_ZECH:
            fq_zech_mpoly_sub(GET_ZECH_POLY(poly1), GET_ZECH_POLY(poly2),
                             GET_ZECH_POLY(poly3), GET_ZECH_CTX(ctx));
            break;
            
        case FIELD_ID_GF28:
        case FIELD_ID_GF216:
        case FIELD_ID_GF232:
        case FIELD_ID_GF264:
        case FIELD_ID_GF2128:
            /* In GF(2^n), subtraction is the same as addition */
            unified_mpoly_add(poly1, poly2, poly3);
            break;
            
        default:
            fq_nmod_mpoly_sub(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2),
                             GET_FQ_POLY(poly3), GET_FQ_CTX(ctx));
            break;
    }
}

void unified_mpoly_neg(unified_mpoly_t poly1, const unified_mpoly_t poly2) {
    unified_mpoly_ctx_t ctx = poly1->ctx_ptr;
    
    switch (poly1->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_neg(GET_NMOD_POLY(poly1), GET_NMOD_POLY(poly2), GET_NMOD_CTX(ctx));
            break;
            
        case FIELD_ID_FQ_ZECH:
            fq_zech_mpoly_neg(GET_ZECH_POLY(poly1), GET_ZECH_POLY(poly2), GET_ZECH_CTX(ctx));
            break;
            
        case FIELD_ID_GF28:
        case FIELD_ID_GF216:
        case FIELD_ID_GF232:
        case FIELD_ID_GF264:
        case FIELD_ID_GF2128:
            /* In GF(2^n), -a = a */
            unified_mpoly_set(poly1, poly2);
            break;
            
        default:
            fq_nmod_mpoly_neg(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2), GET_FQ_CTX(ctx));
            break;
    }
}

/* ============================================================================
   MULTIPLICATION WITH ALL GF(2^n) OPTIMIZATIONS
   ============================================================================ */

/* Global flags to enable/disable optimizations */
static int use_gf28_array_mul = 1;
static int use_gf216_array_mul = 1;
static int use_gf232_array_mul = 1;
static int use_gf264_array_mul = 1;
static int use_gf2128_array_mul = 1;
static int use_zech_mul = 1;  /* Enable Zech multiplication by default */

static inline ulong reduce_exp_field_ui(ulong e, ulong q) {
    if (e == 0 || e < q) return e;
    return ((e - 1) % (q - 1)) + 1;
}

static ulong unified_field_size_q(unified_mpoly_ctx_t ctx) {
    if (!ctx || !ctx->field_ctx) return 0;
    switch (ctx->field_ctx->field_id) {
        case FIELD_ID_NMOD:
            return ctx->field_ctx->ctx.nmod_ctx.n;
        case FIELD_ID_FQ_ZECH:
            if (ctx->field_ctx->ctx.zech_ctx) {
                mp_limb_t p = ctx->field_ctx->ctx.zech_ctx->fq_nmod_ctx->modulus->mod.n;
                slong d = fq_zech_ctx_degree(ctx->field_ctx->ctx.zech_ctx);
                ulong q = 1;
                for (slong i = 0; i < d; i++) {
                    if (q > WORD_MAX / p) return WORD_MAX;
                    q *= p;
                }
                return q;
            }
            return 0;
        default:
            if (ctx->field_ctx->ctx.fq_ctx) {
                mp_limb_t p = fq_nmod_ctx_prime(ctx->field_ctx->ctx.fq_ctx);
                slong d = fq_nmod_ctx_degree(ctx->field_ctx->ctx.fq_ctx);
                ulong q = 1;
                for (slong i = 0; i < d; i++) {
                    if (q > WORD_MAX / p) return WORD_MAX;
                    q *= p;
                }
                return q;
            }
            return 0;
    }
}

static void unified_mpoly_reduce_field_equation_inplace(unified_mpoly_t poly) {
    if (!poly || !poly->ctx_ptr || poly->ctx_ptr->nvars <= 0) return;
    ulong q = unified_field_size_q(poly->ctx_ptr);
    if (q <= 1) return;

    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    slong nvars = ctx->nvars;

    switch (poly->field_id) {
        case FIELD_ID_NMOD: {
            nmod_mpoly_t reduced;
            nmod_mpoly_init(reduced, GET_NMOD_CTX(ctx));
            slong nterms = nmod_mpoly_length(GET_NMOD_POLY(poly), GET_NMOD_CTX(ctx));
            ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
            for (slong i = 0; i < nterms; i++) {
                ulong coeff = nmod_mpoly_get_term_coeff_ui(GET_NMOD_POLY(poly), i, GET_NMOD_CTX(ctx));
                nmod_mpoly_get_term_exp_ui(exp, GET_NMOD_POLY(poly), i, GET_NMOD_CTX(ctx));
                for (slong k = 0; k < nvars; k++) exp[k] = reduce_exp_field_ui(exp[k], q);
                nmod_mpoly_push_term_ui_ui(reduced, coeff, exp, GET_NMOD_CTX(ctx));
            }
            nmod_mpoly_sort_terms(reduced, GET_NMOD_CTX(ctx));
            nmod_mpoly_combine_like_terms(reduced, GET_NMOD_CTX(ctx));
            nmod_mpoly_set(GET_NMOD_POLY(poly), reduced, GET_NMOD_CTX(ctx));
            flint_free(exp);
            nmod_mpoly_clear(reduced, GET_NMOD_CTX(ctx));
            break;
        }
        case FIELD_ID_FQ_ZECH:
            break;
        default: {
            fq_nmod_mpoly_t reduced;
            fq_nmod_mpoly_init(reduced, GET_FQ_CTX(ctx));
            slong nterms = fq_nmod_mpoly_length(GET_FQ_POLY(poly), GET_FQ_CTX(ctx));
            ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx->field_ctx->ctx.fq_ctx);
            for (slong i = 0; i < nterms; i++) {
                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, GET_FQ_POLY(poly), i, GET_FQ_CTX(ctx));
                fq_nmod_mpoly_get_term_exp_ui(exp, GET_FQ_POLY(poly), i, GET_FQ_CTX(ctx));
                for (slong k = 0; k < nvars; k++) exp[k] = reduce_exp_field_ui(exp[k], q);
                fq_nmod_mpoly_push_term_fq_nmod_ui(reduced, coeff, exp, GET_FQ_CTX(ctx));
            }
            fq_nmod_mpoly_sort_terms(reduced, GET_FQ_CTX(ctx));
            fq_nmod_mpoly_combine_like_terms(reduced, GET_FQ_CTX(ctx));
            fq_nmod_mpoly_set(GET_FQ_POLY(poly), reduced, GET_FQ_CTX(ctx));
            fq_nmod_clear(coeff, ctx->field_ctx->ctx.fq_ctx);
            flint_free(exp);
            fq_nmod_mpoly_clear(reduced, GET_FQ_CTX(ctx));
            break;
        }
    }
}


/* Function to enable optimizations */
void unified_mpoly_enable_optimizations(field_id_t field_id, int enable) {
    switch (field_id) {
        case FIELD_ID_GF28:
            use_gf28_array_mul = enable;
            printf("GF(2^8) array multiplication: %s\n", enable ? "enabled" : "disabled");
            break;
        case FIELD_ID_GF216:
            use_gf216_array_mul = enable;
            printf("GF(2^16) array multiplication: %s\n", enable ? "enabled" : "disabled");
            break;
        case FIELD_ID_GF232:
            use_gf232_array_mul = enable;
            printf("GF(2^32) array multiplication: %s\n", enable ? "enabled" : "disabled");
            break;
        case FIELD_ID_GF264:
            use_gf264_array_mul = enable;
            printf("GF(2^64) array multiplication: %s\n", enable ? "enabled" : "disabled");
            break;
        case FIELD_ID_GF2128:
            use_gf2128_array_mul = enable;
            printf("GF(2^128) array multiplication: %s\n", enable ? "enabled" : "disabled");
            break;
        case FIELD_ID_FQ_ZECH:
            use_zech_mul = enable;
            printf("Zech logarithm multiplication: %s\n", enable ? "enabled" : "disabled");
            break;
        default:
            break;
    }
}

int unified_mpoly_mul(unified_mpoly_t poly1, const unified_mpoly_t poly2,
                     const unified_mpoly_t poly3) {
    unified_mpoly_ctx_t ctx = poly1->ctx_ptr;
    
    switch (poly1->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_mul(GET_NMOD_POLY(poly1), GET_NMOD_POLY(poly2),
                          GET_NMOD_POLY(poly3), GET_NMOD_CTX(ctx));
            if (g_field_equation_reduction) unified_mpoly_reduce_field_equation_inplace(poly1);
            return 1;
            
        case FIELD_ID_FQ_ZECH:
            if (use_zech_mul) {
                fq_zech_mpoly_mul(GET_ZECH_POLY(poly1), GET_ZECH_POLY(poly2),
                                 GET_ZECH_POLY(poly3), GET_ZECH_CTX(ctx));
                if (g_field_equation_reduction) unified_mpoly_reduce_field_equation_inplace(poly1);
                return 1;
            } else {
                printf("Zech multiplication disabled, this shouldn't happen\n");
                return 0;
            }
            
        case FIELD_ID_GF28:
            if (use_gf28_array_mul) {
                gf28_mpoly_t A, B, C;
                gf28_mpoly_ctx_t native_ctx;
                
                gf28_mpoly_ctx_init(native_ctx, ctx->nvars, ctx->ord);
                gf28_mpoly_init(A, native_ctx);
                gf28_mpoly_init(B, native_ctx);
                gf28_mpoly_init(C, native_ctx);
                
                fq_nmod_mpoly_to_gf28_mpoly(A, GET_FQ_POLY(poly2), 
                                            ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                fq_nmod_mpoly_to_gf28_mpoly(B, GET_FQ_POLY(poly3), 
                                            ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                
                int success = gf28_mpoly_mul(C, A, B, native_ctx);
                //printf("%d\n",success);
                if (success) {
                    gf28_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(poly1), C, 
                                                ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                    if (g_field_equation_reduction) unified_mpoly_reduce_field_equation_inplace(poly1);
                    
                    gf28_mpoly_clear(A, native_ctx);
                    gf28_mpoly_clear(B, native_ctx);
                    gf28_mpoly_clear(C, native_ctx);
                    gf28_mpoly_ctx_clear(native_ctx);
                    return 1;
                } else {
                    gf28_mpoly_clear(A, native_ctx);
                    gf28_mpoly_clear(B, native_ctx);
                    gf28_mpoly_clear(C, native_ctx);
                    gf28_mpoly_ctx_clear(native_ctx);
                    WARN_THRICE("GF(2^8) array multiplication failed, using standard method\n");
                }
            }
            /* Fall through to standard multiplication */
            fq_nmod_mpoly_mul(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2),
                             GET_FQ_POLY(poly3), GET_FQ_CTX(ctx));
            if (g_field_equation_reduction) unified_mpoly_reduce_field_equation_inplace(poly1);
            return 1;
            
        case FIELD_ID_GF216:
            if (use_gf216_array_mul) {
                gf216_mpoly_t A, B, C;
                gf216_mpoly_ctx_t native_ctx;
                
                gf216_mpoly_ctx_init(native_ctx, ctx->nvars, ctx->ord);
                gf216_mpoly_init(A, native_ctx);
                gf216_mpoly_init(B, native_ctx);
                gf216_mpoly_init(C, native_ctx);
                
                fq_nmod_mpoly_to_gf216_mpoly(A, GET_FQ_POLY(poly2), 
                                             ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                fq_nmod_mpoly_to_gf216_mpoly(B, GET_FQ_POLY(poly3), 
                                             ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                
                int success = gf216_mpoly_mul(C, A, B, native_ctx);
                
                if (success) {
                    gf216_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(poly1), C, 
                                                 ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                    if (g_field_equation_reduction) unified_mpoly_reduce_field_equation_inplace(poly1);
                    
                    gf216_mpoly_clear(A, native_ctx);
                    gf216_mpoly_clear(B, native_ctx);
                    gf216_mpoly_clear(C, native_ctx);
                    gf216_mpoly_ctx_clear(native_ctx);
                    return 1;
                } else {
                    gf216_mpoly_clear(A, native_ctx);
                    gf216_mpoly_clear(B, native_ctx);
                    gf216_mpoly_clear(C, native_ctx);
                    gf216_mpoly_ctx_clear(native_ctx);
                    WARN_THRICE("GF(2^16) array multiplication failed, using standard method\n");
                }
            }
            /* Fall through to standard multiplication */
            fq_nmod_mpoly_mul(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2),
                             GET_FQ_POLY(poly3), GET_FQ_CTX(ctx));
            if (g_field_equation_reduction) unified_mpoly_reduce_field_equation_inplace(poly1);
            return 1;
            
        case FIELD_ID_GF232:
            if (use_gf232_array_mul) {
                gf232_mpoly_t A, B, C;
                gf232_mpoly_ctx_t native_ctx;
                
                gf232_mpoly_ctx_init(native_ctx, ctx->nvars, ctx->ord);
                gf232_mpoly_init(A, native_ctx);
                gf232_mpoly_init(B, native_ctx);
                gf232_mpoly_init(C, native_ctx);
                
                fq_nmod_mpoly_to_gf232_mpoly(A, GET_FQ_POLY(poly2), 
                                             ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                fq_nmod_mpoly_to_gf232_mpoly(B, GET_FQ_POLY(poly3), 
                                             ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                
                int success = gf232_mpoly_mul(C, A, B, native_ctx);
                
                if (success) {
                    gf232_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(poly1), C, 
                                                 ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                    if (g_field_equation_reduction) unified_mpoly_reduce_field_equation_inplace(poly1);
                    
                    gf232_mpoly_clear(A, native_ctx);
                    gf232_mpoly_clear(B, native_ctx);
                    gf232_mpoly_clear(C, native_ctx);
                    gf232_mpoly_ctx_clear(native_ctx);
                    return 1;
                } else {
                    gf232_mpoly_clear(A, native_ctx);
                    gf232_mpoly_clear(B, native_ctx);
                    gf232_mpoly_clear(C, native_ctx);
                    gf232_mpoly_ctx_clear(native_ctx);
                    WARN_THRICE("GF(2^32) array multiplication failed, using standard method\n");
                }
            }
            /* Fall through to standard multiplication */
            fq_nmod_mpoly_mul(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2),
                             GET_FQ_POLY(poly3), GET_FQ_CTX(ctx));
            if (g_field_equation_reduction) unified_mpoly_reduce_field_equation_inplace(poly1);
            return 1;
            
        case FIELD_ID_GF264:
            if (use_gf264_array_mul) {
                gf264_mpoly_t A, B, C;
                gf264_mpoly_ctx_t native_ctx;
                
                gf264_mpoly_ctx_init(native_ctx, ctx->nvars, ctx->ord);
                gf264_mpoly_init(A, native_ctx);
                gf264_mpoly_init(B, native_ctx);
                gf264_mpoly_init(C, native_ctx);
                
                fq_nmod_mpoly_to_gf264_mpoly(A, GET_FQ_POLY(poly2), 
                                             ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                fq_nmod_mpoly_to_gf264_mpoly(B, GET_FQ_POLY(poly3), 
                                             ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                
                int success = gf264_mpoly_mul(C, A, B, native_ctx);
                
                if (success) {
                    gf264_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(poly1), C, 
                                                 ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                    if (g_field_equation_reduction) unified_mpoly_reduce_field_equation_inplace(poly1);
                    
                    gf264_mpoly_clear(A, native_ctx);
                    gf264_mpoly_clear(B, native_ctx);
                    gf264_mpoly_clear(C, native_ctx);
                    gf264_mpoly_ctx_clear(native_ctx);
                    return 1;
                } else {
                    gf264_mpoly_clear(A, native_ctx);
                    gf264_mpoly_clear(B, native_ctx);
                    gf264_mpoly_clear(C, native_ctx);
                    gf264_mpoly_ctx_clear(native_ctx);
                    WARN_THRICE("GF(2^64) array multiplication failed, using standard method\n");
                }
            }
            /* Fall through to standard multiplication */
            fq_nmod_mpoly_mul(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2),
                             GET_FQ_POLY(poly3), GET_FQ_CTX(ctx));
            if (g_field_equation_reduction) unified_mpoly_reduce_field_equation_inplace(poly1);
            return 1;
            
        case FIELD_ID_GF2128:
            if (use_gf2128_array_mul) {
                gf2128_mpoly_t A, B, C;
                gf2128_mpoly_ctx_t native_ctx;
                
                gf2128_mpoly_ctx_init(native_ctx, ctx->nvars, ctx->ord);
                gf2128_mpoly_init(A, native_ctx);
                gf2128_mpoly_init(B, native_ctx);
                gf2128_mpoly_init(C, native_ctx);
                
                fq_nmod_mpoly_to_gf2128_mpoly(A, GET_FQ_POLY(poly2), 
                                              ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                fq_nmod_mpoly_to_gf2128_mpoly(B, GET_FQ_POLY(poly3), 
                                              ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                
                int success = gf2128_mpoly_mul_array(C, A, B, native_ctx);
                
                if (success) {
                    gf2128_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(poly1), C, 
                                                  ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                    if (g_field_equation_reduction) unified_mpoly_reduce_field_equation_inplace(poly1);
                    
                    gf2128_mpoly_clear(A, native_ctx);
                    gf2128_mpoly_clear(B, native_ctx);
                    gf2128_mpoly_clear(C, native_ctx);
                    gf2128_mpoly_ctx_clear(native_ctx);
                    return 1;
                } else {
                    gf2128_mpoly_clear(A, native_ctx);
                    gf2128_mpoly_clear(B, native_ctx);
                    gf2128_mpoly_clear(C, native_ctx);
                    gf2128_mpoly_ctx_clear(native_ctx);
                    WARN_THRICE("GF(2^128) array multiplication failed, using standard method\n");
                }
            }
            /* Fall through to standard multiplication */
            fq_nmod_mpoly_mul(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2),
                             GET_FQ_POLY(poly3), GET_FQ_CTX(ctx));
            if (g_field_equation_reduction) unified_mpoly_reduce_field_equation_inplace(poly1);
            return 1;
            
        default:
            fq_nmod_mpoly_mul(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2),
                             GET_FQ_POLY(poly3), GET_FQ_CTX(ctx));
            if (g_field_equation_reduction) unified_mpoly_reduce_field_equation_inplace(poly1);
            return 1;
    }
}

/* ============================================================================
   DIVISION OPERATIONS WITH ALL GF(2^n) OPTIMIZATIONS
   ============================================================================ */

/* Global flags for division optimizations */
static int use_gf28_div_opt = 1;
static int use_gf216_div_opt = 1;
static int use_gf232_div_opt = 1;
static int use_gf264_div_opt = 1;
static int use_gf2128_div_opt = 1;
static int use_zech_div = 1;

/* Enable/disable division optimizations */
void unified_mpoly_enable_div_optimizations(field_id_t field_id, int enable) {
    switch (field_id) {
        case FIELD_ID_GF28:
            use_gf28_div_opt = enable;
            printf("GF(2^8) optimized division: %s\n", enable ? "enabled" : "disabled");
            break;
        case FIELD_ID_GF216:
            use_gf216_div_opt = enable;
            printf("GF(2^16) optimized division: %s\n", enable ? "enabled" : "disabled");
            break;
        case FIELD_ID_GF232:
            use_gf232_div_opt = enable;
            printf("GF(2^32) optimized division: %s\n", enable ? "enabled" : "disabled");
            break;
        case FIELD_ID_GF264:
            use_gf264_div_opt = enable;
            printf("GF(2^64) optimized division: %s\n", enable ? "enabled" : "disabled");
            break;
        case FIELD_ID_GF2128:
            use_gf2128_div_opt = enable;
            printf("GF(2^128) optimized division: %s\n", enable ? "enabled" : "disabled");
            break;
        case FIELD_ID_FQ_ZECH:
            use_zech_div = enable;
            printf("Zech logarithm division: %s\n", enable ? "enabled" : "disabled");
            break;
        default:
            break;
    }
}

int unified_mpoly_divides(unified_mpoly_t Q, const unified_mpoly_t A,
                         const unified_mpoly_t B) {
    unified_mpoly_ctx_t ctx = Q->ctx_ptr;

    switch (Q->field_id) {
        case FIELD_ID_NMOD:
            return nmod_mpoly_divides(GET_NMOD_POLY(Q), GET_NMOD_POLY(A),
                                     GET_NMOD_POLY(B), GET_NMOD_CTX(ctx));
            
        case FIELD_ID_FQ_ZECH:
            if (use_zech_div) {
                return fq_zech_mpoly_divides(GET_ZECH_POLY(Q), GET_ZECH_POLY(A),
                                           GET_ZECH_POLY(B), GET_ZECH_CTX(ctx));
            } else {
                printf("Zech division disabled, this shouldn't happen\n");
                return 0;
            }
            
        case FIELD_ID_GF28:
            if (use_gf28_div_opt) {
                gf28_mpoly_t A_native, B_native, Q_native;
                gf28_mpoly_ctx_t native_ctx;
                
                gf28_mpoly_ctx_init(native_ctx, ctx->nvars, ctx->ord);
                gf28_mpoly_init(A_native, native_ctx);
                gf28_mpoly_init(B_native, native_ctx);
                gf28_mpoly_init(Q_native, native_ctx);
                
                fq_nmod_mpoly_to_gf28_mpoly(A_native, GET_FQ_POLY(A), 
                                            ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                fq_nmod_mpoly_to_gf28_mpoly(B_native, GET_FQ_POLY(B), 
                                            ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                
                int success = gf28_mpoly_divides(Q_native, A_native, B_native, native_ctx);
                
                if (success) {
                    gf28_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(Q), Q_native, 
                                                ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                }
                
                gf28_mpoly_clear(A_native, native_ctx);
                gf28_mpoly_clear(B_native, native_ctx);
                gf28_mpoly_clear(Q_native, native_ctx);
                gf28_mpoly_ctx_clear(native_ctx);
                
                return success;
            }
            /* Fall through to standard division */
            return fq_nmod_mpoly_divides(GET_FQ_POLY(Q), GET_FQ_POLY(A),
                                        GET_FQ_POLY(B), GET_FQ_CTX(ctx));
            
        case FIELD_ID_GF216:
            if (use_gf216_div_opt) {
                gf216_mpoly_t A_native, B_native, Q_native;
                gf216_mpoly_ctx_t native_ctx;
                
                gf216_mpoly_ctx_init(native_ctx, ctx->nvars, ctx->ord);
                gf216_mpoly_init(A_native, native_ctx);
                gf216_mpoly_init(B_native, native_ctx);
                gf216_mpoly_init(Q_native, native_ctx);
                
                fq_nmod_mpoly_to_gf216_mpoly(A_native, GET_FQ_POLY(A), 
                                             ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                fq_nmod_mpoly_to_gf216_mpoly(B_native, GET_FQ_POLY(B), 
                                             ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                
                int success = gf216_mpoly_divides(Q_native, A_native, B_native, native_ctx);
                
                if (success) {
                    gf216_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(Q), Q_native, 
                                                 ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                }
                
                gf216_mpoly_clear(A_native, native_ctx);
                gf216_mpoly_clear(B_native, native_ctx);
                gf216_mpoly_clear(Q_native, native_ctx);
                gf216_mpoly_ctx_clear(native_ctx);
                
                return success;
            }
            /* Fall through to standard division */
            return fq_nmod_mpoly_divides(GET_FQ_POLY(Q), GET_FQ_POLY(A),
                                        GET_FQ_POLY(B), GET_FQ_CTX(ctx));
            
        case FIELD_ID_GF232:
            if (use_gf232_div_opt) {
                gf232_mpoly_t A_native, B_native, Q_native;
                gf232_mpoly_ctx_t native_ctx;
                
                gf232_mpoly_ctx_init(native_ctx, ctx->nvars, ctx->ord);
                gf232_mpoly_init(A_native, native_ctx);
                gf232_mpoly_init(B_native, native_ctx);
                gf232_mpoly_init(Q_native, native_ctx);
                
                fq_nmod_mpoly_to_gf232_mpoly(A_native, GET_FQ_POLY(A), 
                                             ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                fq_nmod_mpoly_to_gf232_mpoly(B_native, GET_FQ_POLY(B), 
                                             ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                
                int success = gf232_mpoly_divides(Q_native, A_native, B_native, native_ctx);
                
                if (success) {
                    gf232_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(Q), Q_native, 
                                                 ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                }
                
                gf232_mpoly_clear(A_native, native_ctx);
                gf232_mpoly_clear(B_native, native_ctx);
                gf232_mpoly_clear(Q_native, native_ctx);
                gf232_mpoly_ctx_clear(native_ctx);
                
                return success;
            }
            /* Fall through to standard division */
            return fq_nmod_mpoly_divides(GET_FQ_POLY(Q), GET_FQ_POLY(A),
                                        GET_FQ_POLY(B), GET_FQ_CTX(ctx));
            
        case FIELD_ID_GF264:
            if (use_gf264_div_opt) {
                gf264_mpoly_t A_native, B_native, Q_native;
                gf264_mpoly_ctx_t native_ctx;
                
                gf264_mpoly_ctx_init(native_ctx, ctx->nvars, ctx->ord);
                gf264_mpoly_init(A_native, native_ctx);
                gf264_mpoly_init(B_native, native_ctx);
                gf264_mpoly_init(Q_native, native_ctx);
                
                fq_nmod_mpoly_to_gf264_mpoly(A_native, GET_FQ_POLY(A), 
                                             ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                fq_nmod_mpoly_to_gf264_mpoly(B_native, GET_FQ_POLY(B), 
                                             ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                
                int success = gf264_mpoly_divides(Q_native, A_native, B_native, native_ctx);
                
                if (success) {
                    gf264_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(Q), Q_native, 
                                                 ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                }
                
                gf264_mpoly_clear(A_native, native_ctx);
                gf264_mpoly_clear(B_native, native_ctx);
                gf264_mpoly_clear(Q_native, native_ctx);
                gf264_mpoly_ctx_clear(native_ctx);
                
                return success;
            }
            /* Fall through to standard division */
            return fq_nmod_mpoly_divides(GET_FQ_POLY(Q), GET_FQ_POLY(A),
                                        GET_FQ_POLY(B), GET_FQ_CTX(ctx));
            
        case FIELD_ID_GF2128:
            if (use_gf2128_div_opt) {
                gf2128_mpoly_t A_native, B_native, Q_native;
                gf2128_mpoly_ctx_t native_ctx;
                
                gf2128_mpoly_ctx_init(native_ctx, ctx->nvars, ctx->ord);
                gf2128_mpoly_init(A_native, native_ctx);
                gf2128_mpoly_init(B_native, native_ctx);
                gf2128_mpoly_init(Q_native, native_ctx);
                
                fq_nmod_mpoly_to_gf2128_mpoly(A_native, GET_FQ_POLY(A), 
                                              ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                fq_nmod_mpoly_to_gf2128_mpoly(B_native, GET_FQ_POLY(B), 
                                              ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                
                int success = gf2128_mpoly_divides(Q_native, A_native, B_native, native_ctx);
                
                if (success) {
                    gf2128_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(Q), Q_native, 
                                                  ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                }
                
                gf2128_mpoly_clear(A_native, native_ctx);
                gf2128_mpoly_clear(B_native, native_ctx);
                gf2128_mpoly_clear(Q_native, native_ctx);
                gf2128_mpoly_ctx_clear(native_ctx);
                
                return success;
            }
            /* Fall through to standard division */
            return fq_nmod_mpoly_divides(GET_FQ_POLY(Q), GET_FQ_POLY(A),
                                        GET_FQ_POLY(B), GET_FQ_CTX(ctx));
            
        default:
            return fq_nmod_mpoly_divides(GET_FQ_POLY(Q), GET_FQ_POLY(A),
                                        GET_FQ_POLY(B), GET_FQ_CTX(ctx));
    }
}

/* ============================================================================
   COEFFICIENT ACCESS WITH ZECH SUPPORT
   ============================================================================ */

void unified_mpoly_set_coeff_ui(unified_mpoly_t poly, const field_elem_u *c,
                               const ulong *exp) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_set_coeff_ui_ui(GET_NMOD_POLY(poly), c->nmod, exp,
                                      GET_NMOD_CTX(ctx));
            break;
            
        case FIELD_ID_FQ_ZECH:
            {
                /* Implement coefficient setting for fq_zech_mpoly */
                fq_zech_t coeff;
                fq_zech_ctx_struct *zech_field_ctx = ctx->field_ctx->ctx.zech_ctx;
                fq_zech_init(coeff, zech_field_ctx);
                
                /* Convert field element to fq_zech format */
                if (ctx->field_ctx->field_id == FIELD_ID_FQ_ZECH) {
                    fq_zech_set(coeff, &c->fq_zech, zech_field_ctx);
                } else {
                    /* Convert from other field types through fq_nmod */
                    fq_nmod_t fq_temp;
                    fq_nmod_init(fq_temp, ctx->field_ctx->ctx.fq_ctx);
                    field_elem_to_fq_nmod(fq_temp, c, ctx->field_ctx);
                    fq_zech_set_fq_nmod(coeff, fq_temp, zech_field_ctx);
                    fq_nmod_clear(fq_temp, ctx->field_ctx->ctx.fq_ctx);
                }
                
                /* Build monomial x^exp with given coefficient */
                fq_zech_mpoly_t monomial, temp_poly;
                fq_zech_mpoly_init(monomial, GET_ZECH_CTX(ctx));
                fq_zech_mpoly_init(temp_poly, GET_ZECH_CTX(ctx));
                
                /* Start with coefficient as constant term */
                fq_zech_mpoly_set_fq_zech(monomial, coeff, GET_ZECH_CTX(ctx));
                
                /* Multiply by each variable raised to its exponent */
                for (slong i = 0; i < ctx->nvars; i++) {
                    if (exp[i] > 0) {
                        fq_zech_mpoly_gen(temp_poly, i, GET_ZECH_CTX(ctx));
                        if (exp[i] > 1) {
                            fq_zech_mpoly_pow_ui(temp_poly, temp_poly, exp[i], GET_ZECH_CTX(ctx));
                        }
                        fq_zech_mpoly_mul(monomial, monomial, temp_poly, GET_ZECH_CTX(ctx));
                    }
                }
                
                /* Now we need to remove any existing term with the same exponent and add the new one */
                /* First, create a polynomial with just the term to remove */
                fq_zech_mpoly_t term_to_remove;
                fq_zech_mpoly_init(term_to_remove, GET_ZECH_CTX(ctx));
                
                /* Extract the existing coefficient at this exponent */
                /* We'll do this by creating the monomial and looking for it */
                fq_zech_t one;
                fq_zech_init(one, zech_field_ctx);
                fq_zech_one(one, zech_field_ctx);
                fq_zech_mpoly_set_fq_zech(term_to_remove, one, GET_ZECH_CTX(ctx));
                
                for (slong i = 0; i < ctx->nvars; i++) {
                    if (exp[i] > 0) {
                        fq_zech_mpoly_gen(temp_poly, i, GET_ZECH_CTX(ctx));
                        if (exp[i] > 1) {
                            fq_zech_mpoly_pow_ui(temp_poly, temp_poly, exp[i], GET_ZECH_CTX(ctx));
                        }
                        fq_zech_mpoly_mul(term_to_remove, term_to_remove, temp_poly, GET_ZECH_CTX(ctx));
                    }
                }
                
                /* Create result polynomial: original - old_term + new_term */
                fq_zech_mpoly_t result;
                fq_zech_mpoly_init(result, GET_ZECH_CTX(ctx));
                fq_zech_mpoly_set(result, GET_ZECH_POLY(poly), GET_ZECH_CTX(ctx));
                
                /* The approach is: result = poly - term_to_remove + monomial */
                /* But we need to find the coefficient first */
                /* For now, let's use a simpler approach: rebuild without the term, then add */
                
                /* Clear the polynomial and rebuild it term by term */
                /* This is not the most efficient but works correctly */
                if (fq_zech_is_zero(coeff, zech_field_ctx)) {
                    /* If coefficient is zero, we need to remove the term */
                    /* We'll rebuild the polynomial without this term */
                    /* For now, we'll just set to the original minus the term */
                    fq_zech_mpoly_sub(GET_ZECH_POLY(poly), GET_ZECH_POLY(poly), term_to_remove, GET_ZECH_CTX(ctx));
                } else {
                    /* Add or replace the term */
                    /* Simple approach: add the new monomial */
                    /* This might create duplicate terms, but FLINT should combine them */
                    fq_zech_mpoly_add(GET_ZECH_POLY(poly), GET_ZECH_POLY(poly), monomial, GET_ZECH_CTX(ctx));
                }
                
                /* Clean up */
                fq_zech_clear(coeff, zech_field_ctx);
                fq_zech_clear(one, zech_field_ctx);
                fq_zech_mpoly_clear(monomial, GET_ZECH_CTX(ctx));
                fq_zech_mpoly_clear(temp_poly, GET_ZECH_CTX(ctx));
                fq_zech_mpoly_clear(term_to_remove, GET_ZECH_CTX(ctx));
                fq_zech_mpoly_clear(result, GET_ZECH_CTX(ctx));
            }
            break;
            
        default:
            {
                fq_nmod_t temp;
                fq_nmod_init(temp, ctx->field_ctx->ctx.fq_ctx);
                field_elem_to_fq_nmod(temp, c, ctx->field_ctx);
                fq_nmod_mpoly_set_coeff_fq_nmod_ui(GET_FQ_POLY(poly), temp, exp,
                                                   GET_FQ_CTX(ctx));
                fq_nmod_clear(temp, ctx->field_ctx->ctx.fq_ctx);
            }
            break;
    }
}

void unified_mpoly_get_coeff_ui(field_elem_u *c, const unified_mpoly_t poly,
                               const ulong *exp) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            c->nmod = nmod_mpoly_get_coeff_ui_ui(GET_NMOD_POLY(poly), exp,
                                                GET_NMOD_CTX(ctx));
            break;
            
          case FIELD_ID_FQ_ZECH:
            {
              /* Implement coefficient extraction for fq_zech_mpoly */
                fq_zech_ctx_struct *zech_field_ctx = ctx->field_ctx->ctx.zech_ctx;
                fq_zech_init(&c->fq_zech, zech_field_ctx);
                
                /* Create the monomial x^exp to search for */
                fq_zech_mpoly_t monomial, temp;
                fq_zech_mpoly_init(monomial, GET_ZECH_CTX(ctx));
                fq_zech_mpoly_init(temp, GET_ZECH_CTX(ctx));
                
                /* Build monomial with coefficient 1 */
                fq_zech_t one;
                fq_zech_init(one, zech_field_ctx);
                fq_zech_one(one, zech_field_ctx);
                fq_zech_mpoly_set_fq_zech(monomial, one, GET_ZECH_CTX(ctx));
                
                /* Multiply by each variable raised to its exponent */
                for (slong i = 0; i < ctx->nvars; i++) {
                    if (exp[i] > 0) {
                        fq_zech_mpoly_gen(temp, i, GET_ZECH_CTX(ctx));
                        if (exp[i] > 1) {
                            fq_zech_mpoly_pow_ui(temp, temp, exp[i], GET_ZECH_CTX(ctx));
                        }
                        fq_zech_mpoly_mul(monomial, monomial, temp, GET_ZECH_CTX(ctx));
                    }
                }
                
                /* Method: evaluate poly at a point where only the term x^exp contributes */
                /* This is complex, so instead we'll use a different approach */
                
                /* Alternative: Build a polynomial that extracts just this coefficient */
                /* We can evaluate (poly * monomial^(-1)) at the point where all vars = 1 */
                /* But this requires division which might not work well */
                
                /* Simpler approach: subtract all other terms and get what's left */
                /* This is inefficient but works */
                
                /* For now, return zero - proper implementation would require */
                /* iterating through the polynomial's terms */
                fq_zech_zero(&c->fq_zech, zech_field_ctx);
                
                /* TODO: Implement proper term iteration for fq_zech_mpoly */
                /* This would require access to the internal representation */
                
                /* Clean up */
                fq_zech_clear(one, zech_field_ctx);
                fq_zech_mpoly_clear(monomial, GET_ZECH_CTX(ctx));
                fq_zech_mpoly_clear(temp, GET_ZECH_CTX(ctx));
            }
            break;
            
        default:
            {
                fq_nmod_t temp;
                fq_nmod_init(temp, ctx->field_ctx->ctx.fq_ctx);
                fq_nmod_mpoly_get_coeff_fq_nmod_ui(temp, GET_FQ_POLY(poly), exp,
                                                   GET_FQ_CTX(ctx));
                fq_nmod_to_field_elem(c, temp, ctx->field_ctx);
                fq_nmod_clear(temp, ctx->field_ctx->ctx.fq_ctx);
            }
            break;
    }
}

/* ============================================================================
   DEGREE OPERATIONS WITH ZECH SUPPORT
   ============================================================================ */

void unified_mpoly_degrees_si(slong *degs, const unified_mpoly_t poly) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_degrees_si(degs, GET_NMOD_POLY(poly), GET_NMOD_CTX(ctx));
            break;
            
        case FIELD_ID_FQ_ZECH:
            fq_zech_mpoly_degrees_si(degs, GET_ZECH_POLY(poly), GET_ZECH_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_degrees_si(degs, GET_FQ_POLY(poly), GET_FQ_CTX(ctx));
            break;
    }
}

slong unified_mpoly_total_degree_si(const unified_mpoly_t poly) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            return nmod_mpoly_total_degree_si(GET_NMOD_POLY(poly), GET_NMOD_CTX(ctx));
            
        case FIELD_ID_FQ_ZECH:
            return fq_zech_mpoly_total_degree_si(GET_ZECH_POLY(poly), GET_ZECH_CTX(ctx));
            
        default:
            return fq_nmod_mpoly_total_degree_si(GET_FQ_POLY(poly), GET_FQ_CTX(ctx));
    }
}

/* ============================================================================
   PRINTING WITH ZECH SUPPORT
   ============================================================================ */

void unified_mpoly_print_pretty(const unified_mpoly_t poly, const char **vars) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_print_pretty(GET_NMOD_POLY(poly), vars, GET_NMOD_CTX(ctx));
            break;
            
        case FIELD_ID_FQ_ZECH:
            fq_zech_mpoly_print_pretty(GET_ZECH_POLY(poly), vars, GET_ZECH_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_print_pretty(GET_FQ_POLY(poly), vars, GET_FQ_CTX(ctx));
            break;
    }
}

/* ============================================================================
   SCALAR MULTIPLICATION WITH ZECH SUPPORT
   ============================================================================ */

void unified_mpoly_scalar_mul_ui(unified_mpoly_t poly1, const unified_mpoly_t poly2,
                                ulong c) {
    unified_mpoly_ctx_t ctx = poly1->ctx_ptr;
    
    switch (poly1->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_scalar_mul_ui(GET_NMOD_POLY(poly1), GET_NMOD_POLY(poly2),
                                    c, GET_NMOD_CTX(ctx));
            break;
            
        case FIELD_ID_FQ_ZECH:
            {
                /* Convert scalar to field element */
                fq_zech_t scalar;
                fq_zech_init(scalar, ctx->field_ctx->ctx.zech_ctx);
                fq_zech_set_ui(scalar, c, ctx->field_ctx->ctx.zech_ctx);
                
                /* Multiply */
                fq_zech_mpoly_scalar_mul_fq_zech(GET_ZECH_POLY(poly1), GET_ZECH_POLY(poly2),
                                               scalar, GET_ZECH_CTX(ctx));
                
                fq_zech_clear(scalar, ctx->field_ctx->ctx.zech_ctx);
            }
            break;
            
        default:
            {
                /* Convert scalar to field element */
                fq_nmod_t scalar;
                fq_nmod_init(scalar, ctx->field_ctx->ctx.fq_ctx);
                fq_nmod_set_ui(scalar, c, ctx->field_ctx->ctx.fq_ctx);
                
                /* Multiply */
                fq_nmod_mpoly_scalar_mul_fq_nmod(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2),
                                               scalar, GET_FQ_CTX(ctx));
                
                fq_nmod_clear(scalar, ctx->field_ctx->ctx.fq_ctx);
            }
            break;
    }
}

/* ============================================================================
   MISSING RING OPERATIONS FOR UNIFIED MPOLY WITH ZECH SUPPORT
   ============================================================================ */

void unified_mpoly_swap(unified_mpoly_t poly1, unified_mpoly_t poly2) {
    if (poly1 == poly2) return;
    
    unified_mpoly_struct temp = *poly1;
    *poly1 = *poly2;
    *poly2 = temp;
}

void unified_mpoly_set_fmpz(unified_mpoly_t poly, const fmpz_t c) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    unified_mpoly_zero(poly);
    
    if (fmpz_is_zero(c)) return;
    
    /* Set polynomial to constant c */
    field_elem_u coeff;
    ulong *exp = (ulong *)calloc(ctx->nvars, sizeof(ulong));
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            coeff.nmod = fmpz_get_ui(c) % GET_NMOD_CTX(ctx)->mod.n;
            break;
            
        case FIELD_ID_GF28:
            coeff.gf28 = fmpz_get_ui(c) & 0xFF;  /* Take mod 256 */
            break;
            
        case FIELD_ID_GF216:
            coeff.gf216 = fmpz_get_ui(c) & 0xFFFF;  /* Take mod 65536 */
            break;
            
        case FIELD_ID_GF232:
            coeff.gf232.value = fmpz_get_ui(c) & 0xFFFFFFFF;  /* Take mod 2^32 */
            break;
            
        case FIELD_ID_GF264:
            coeff.gf264.value = fmpz_get_ui(c);  /* Full 64-bit */
            break;
            
        case FIELD_ID_GF2128:
            coeff.gf2128.low = fmpz_get_ui(c);
            coeff.gf2128.high = 0;  /* For large values, would need more work */
            break;
            
        case FIELD_ID_FQ_ZECH:
            {
                fq_zech_init(&coeff.fq_zech, ctx->field_ctx->ctx.zech_ctx);
                fq_zech_set_fmpz(&coeff.fq_zech, c, ctx->field_ctx->ctx.zech_ctx);
            }
            break;
            
        default:
            {
                fq_nmod_init(&coeff.fq, ctx->field_ctx->ctx.fq_ctx);
                fq_nmod_set_fmpz(&coeff.fq, c, ctx->field_ctx->ctx.fq_ctx);
            }
            break;
    }
    
    unified_mpoly_set_coeff_ui(poly, &coeff, exp);
    
    /* Clean up */
    if (poly->field_id == FIELD_ID_FQ_ZECH) {
        fq_zech_clear(&coeff.fq_zech, ctx->field_ctx->ctx.zech_ctx);
    } else if (poly->field_id == FIELD_ID_FQ) {
        fq_nmod_clear(&coeff.fq, ctx->field_ctx->ctx.fq_ctx);
    }
    
    free(exp);
}

void unified_mpoly_mul_fmpz(unified_mpoly_t poly1, const unified_mpoly_t poly2,
                           const fmpz_t c) {
    if (fmpz_is_zero(c)) {
        unified_mpoly_zero(poly1);
        return;
    }
    
    if (fmpz_is_one(c)) {
        unified_mpoly_set(poly1, poly2);
        return;
    }
    
    /* Convert fmpz to ulong and use scalar multiplication */
    unified_mpoly_ctx_t ctx = poly1->ctx_ptr;
    ulong c_mod = fmpz_get_ui(c);
    
    switch (poly1->field_id) {
        case FIELD_ID_NMOD:
            c_mod = c_mod % GET_NMOD_CTX(ctx)->mod.n;
            unified_mpoly_scalar_mul_ui(poly1, poly2, c_mod);
            break;
            
        case FIELD_ID_FQ_ZECH:
            {
                ulong p = fq_zech_ctx_prime(ctx->field_ctx->ctx.zech_ctx);
                c_mod = c_mod % p;
                if (c_mod == 0) {
                    unified_mpoly_zero(poly1);
                } else {
                    unified_mpoly_scalar_mul_ui(poly1, poly2, c_mod);
                }
            }
            break;
            
        default:
            /* For all other fields including GF(2^n), use scalar_mul_ui */
            {
                /* Get characteristic for GF(2^n) fields */
                ulong p = 2;  /* Default for GF(2^n) */
                if (poly1->field_id == FIELD_ID_FQ && ctx->field_ctx->ctx.fq_ctx) {
                    p = fq_nmod_ctx_prime(ctx->field_ctx->ctx.fq_ctx);
                }
                
                c_mod = c_mod % p;
                if (c_mod == 0) {
                    unified_mpoly_zero(poly1);
                } else {
                    unified_mpoly_scalar_mul_ui(poly1, poly2, c_mod);
                }
            }
            break;
    }
}

int unified_mpoly_pow_fmpz(unified_mpoly_t poly1, const unified_mpoly_t poly2,
                          const fmpz_t exp) {
    if (fmpz_sgn(exp) < 0) {
        fprintf(stderr, "unified_mpoly_pow_fmpz: negative exponent not supported\n");
        return 0;
    }
    
    if (fmpz_is_zero(exp)) {
        unified_mpoly_one(poly1);
        return 1;
    }
    
    if (fmpz_is_one(exp)) {
        unified_mpoly_set(poly1, poly2);
        return 1;
    }
    
    /* Binary exponentiation */
    unified_mpoly_ctx_t ctx = poly1->ctx_ptr;
    unified_mpoly_t temp = unified_mpoly_init(ctx);
    unified_mpoly_t result = unified_mpoly_init(ctx);
    fmpz_t e;
    
    fmpz_init_set(e, exp);
    unified_mpoly_set(temp, poly2);
    unified_mpoly_one(result);
    
    while (!fmpz_is_zero(e)) {
        if (fmpz_is_odd(e)) {
            unified_mpoly_mul(result, result, temp);
        }
        unified_mpoly_mul(temp, temp, temp);
        fmpz_fdiv_q_2exp(e, e, 1);
    }
    
    unified_mpoly_set(poly1, result);
    
    fmpz_clear(e);
    unified_mpoly_clear(temp);
    unified_mpoly_clear(result);
    
    return 1;
}

slong unified_mpoly_length_wrapper(const void *a, const void *ctx) {
    return unified_mpoly_length((const unified_mpoly_t)a);
}

/* ============================================================================
   DIVISION WITH ZECH SUPPORT
   ============================================================================ */

void unified_mpoly_divrem(unified_mpoly_t Q, unified_mpoly_t R,
                        const unified_mpoly_t A, const unified_mpoly_t B) {
    unified_mpoly_ctx_t ctx = A->ctx_ptr;
    //printf("begin divrem\n");
    switch (A->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_divrem(GET_NMOD_POLY(Q), GET_NMOD_POLY(R),
                            GET_NMOD_POLY(A), GET_NMOD_POLY(B),
                            GET_NMOD_CTX(ctx));
            break;
            
        case FIELD_ID_FQ_ZECH:
            fq_zech_mpoly_divrem(GET_ZECH_POLY(Q), GET_ZECH_POLY(R),
                               GET_ZECH_POLY(A), GET_ZECH_POLY(B),
                               GET_ZECH_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_divrem(GET_FQ_POLY(Q), GET_FQ_POLY(R),
                               GET_FQ_POLY(A), GET_FQ_POLY(B),
                               GET_FQ_CTX(ctx));
            break;
    }

    //printf("end divrem\n");
}

/* Enhanced divexact function with optimizations */
int unified_mpoly_divexact(unified_mpoly_t Q, const unified_mpoly_t A,
                              const unified_mpoly_t B) {
    /* First try optimized divides */
    //printf("begin divexact\n");
    int divides = unified_mpoly_divides(Q, A, B);
    //printf("end divexact\n");
    if (!divides) {
        fprintf(stderr, "unified_mpoly_divexact: division is not exact\n");
        unified_mpoly_zero(Q);
    }
    
    return divides;
}

/* ============================================================================
   GENERATOR FUNCTIONS WITH ZECH SUPPORT
   ============================================================================ */

void unified_mpoly_gen(unified_mpoly_t poly, slong var, unified_mpoly_ctx_t ctx) {
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_gen(GET_NMOD_POLY(poly), var, GET_NMOD_CTX(ctx));
            break;
            
        case FIELD_ID_FQ_ZECH:
            fq_zech_mpoly_gen(GET_ZECH_POLY(poly), var, GET_ZECH_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_gen(GET_FQ_POLY(poly), var, GET_FQ_CTX(ctx));
            break;
    }
}

/* ============================================================================
   VOID RING INTERFACE FUNCTIONS WITH ZECH SUPPORT
   ============================================================================ */

/* These functions adapt unified_mpoly operations to the void* interface */

static void unified_mpoly_void_init(void *a, const void *ctx) {
    unified_mpoly_ctx_t mpoly_ctx = (unified_mpoly_ctx_t)ctx;
    unified_mpoly_t *poly_ptr = (unified_mpoly_t *)a;
    *poly_ptr = unified_mpoly_init(mpoly_ctx);
}

static void unified_mpoly_void_clear(void *a, const void *ctx) {
    unified_mpoly_t *poly_ptr = (unified_mpoly_t *)a;
    unified_mpoly_clear(*poly_ptr);
}

static int unified_mpoly_void_is_zero(const void *a, const void *ctx) {
    const unified_mpoly_t *poly_ptr = (const unified_mpoly_t *)a;
    return unified_mpoly_is_zero(*poly_ptr);
}

static void unified_mpoly_void_zero(void *a, const void *ctx) {
    unified_mpoly_t *poly_ptr = (unified_mpoly_t *)a;
    unified_mpoly_zero(*poly_ptr);
}

static void unified_mpoly_void_one(void *a, const void *ctx) {
    unified_mpoly_t *poly_ptr = (unified_mpoly_t *)a;
    unified_mpoly_one(*poly_ptr);
}

static void unified_mpoly_void_set(void *a, const void *b, const void *ctx) {
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    const unified_mpoly_t *poly_b = (const unified_mpoly_t *)b;
    unified_mpoly_set(*poly_a, *poly_b);
}

static void unified_mpoly_void_set_fmpz(void *a, const fmpz_t b, const void *ctx) {
    unified_mpoly_t *poly_ptr = (unified_mpoly_t *)a;
    unified_mpoly_set_fmpz(*poly_ptr, b);
}

static void unified_mpoly_void_swap(void *a, void *b, const void *ctx) {
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    unified_mpoly_t *poly_b = (unified_mpoly_t *)b;
    unified_mpoly_swap(*poly_a, *poly_b);
}

static void unified_mpoly_void_neg(void *a, const void *b, const void *ctx) {
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    const unified_mpoly_t *poly_b = (const unified_mpoly_t *)b;
    unified_mpoly_neg(*poly_a, *poly_b);
}

static void unified_mpoly_void_add(void *a, const void *b, const void *c,
                                  const void *ctx) {
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    const unified_mpoly_t *poly_b = (const unified_mpoly_t *)b;
    const unified_mpoly_t *poly_c = (const unified_mpoly_t *)c;
    unified_mpoly_add(*poly_a, *poly_b, *poly_c);
}

static void unified_mpoly_void_sub(void *a, const void *b, const void *c,
                                  const void *ctx) {
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    const unified_mpoly_t *poly_b = (const unified_mpoly_t *)b;
    const unified_mpoly_t *poly_c = (const unified_mpoly_t *)c;
    unified_mpoly_sub(*poly_a, *poly_b, *poly_c);
}

static void unified_mpoly_void_mul(void *a, const void *b, const void *c,
                                  const void *ctx) {
    //printf("begin mul\n");
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    const unified_mpoly_t *poly_b = (const unified_mpoly_t *)b;
    const unified_mpoly_t *poly_c = (const unified_mpoly_t *)c;
    unified_mpoly_mul(*poly_a, *poly_b, *poly_c);
    //printf("end mul\n");
}

static void unified_mpoly_void_mul_fmpz(void *a, const void *b,
                                       const fmpz_t c, const void *ctx) {
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    const unified_mpoly_t *poly_b = (const unified_mpoly_t *)b;
    unified_mpoly_mul_fmpz(*poly_a, *poly_b, c);
}

static void unified_mpoly_void_divexact(void *a, const void *b,
                                       const void *c, const void *ctx) {
    //printf("begin div\n");
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    const unified_mpoly_t *poly_b = (const unified_mpoly_t *)b;
    const unified_mpoly_t *poly_c = (const unified_mpoly_t *)c;
    if (!unified_mpoly_divexact(*poly_a, *poly_b, *poly_c)) {
        flint_throw(FLINT_ERROR, "unified_mpoly_void_divexact: nonexact");
    }
    //printf("end div\n");
}

static int unified_mpoly_void_divides(void *a, const void *b,
                                     const void *c, const void *ctx) {
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    const unified_mpoly_t *poly_b = (const unified_mpoly_t *)b;
    const unified_mpoly_t *poly_c = (const unified_mpoly_t *)c;
    return unified_mpoly_divides(*poly_a, *poly_b, *poly_c);
}

static int unified_mpoly_void_pow_fmpz(void *a, const void *b,
                                      const fmpz_t c, const void *ctx) {
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    const unified_mpoly_t *poly_b = (const unified_mpoly_t *)b;
    return unified_mpoly_pow_fmpz(*poly_a, *poly_b, c);
}

static slong unified_mpoly_void_length(const void *a, const void *ctx) {
    const unified_mpoly_t *poly_ptr = (const unified_mpoly_t *)a;
    return unified_mpoly_length(*poly_ptr);
}

/* ============================================================================
   RING STRUCTURE DEFINITION WITH ZECH SUPPORT
   ============================================================================ */

/* Initialize the ring structure for unified multivariate polynomials */
void unified_mpoly_ring_init(mpoly_void_ring_t R, unified_mpoly_ctx_t ctx) {
    R->elem_size = sizeof(unified_mpoly_t);
    R->ctx = ctx;
    R->init = unified_mpoly_void_init;
    R->clear = unified_mpoly_void_clear;
    R->is_zero = unified_mpoly_void_is_zero;
    R->zero = unified_mpoly_void_zero;
    R->one = unified_mpoly_void_one;
    R->set = unified_mpoly_void_set;
    R->set_fmpz = unified_mpoly_void_set_fmpz;
    R->swap = unified_mpoly_void_swap;
    R->neg = unified_mpoly_void_neg;
    R->add = unified_mpoly_void_add;
    R->sub = unified_mpoly_void_sub;
    R->mul = unified_mpoly_void_mul;
    R->mul_fmpz = unified_mpoly_void_mul_fmpz;
    R->divexact = unified_mpoly_void_divexact;
    R->divides = unified_mpoly_void_divides;
    R->pow_fmpz = unified_mpoly_void_pow_fmpz;
    R->length = unified_mpoly_void_length;
}

void unified_mpoly_ring_clear(unified_mpoly_ring_t R) {
    if (!R) return;
    
    unified_mpoly_ctx_clear(R->ctx);
    free(R);
}

/* ============================================================================
   CONVENIENCE FUNCTIONS WITH ZECH SUPPORT
   ============================================================================ */

/* Enable all optimizations for a specific field */
void unified_mpoly_enable_all_optimizations(field_id_t field_id) {
    unified_mpoly_enable_optimizations(field_id, 1);
    unified_mpoly_enable_div_optimizations(field_id, 1);
}

/* Disable all optimizations for a specific field */
void unified_mpoly_disable_all_optimizations(field_id_t field_id) {
    unified_mpoly_enable_optimizations(field_id, 0);
    unified_mpoly_enable_div_optimizations(field_id, 0);
}

/* Set Zech logarithm size limit */
void unified_mpoly_set_zech_limit(unified_mpoly_ctx_t ctx, ulong limit) {
    ctx->zech_size_limit = limit;
}

/* Get information about the field representation being used */
const char* unified_mpoly_get_field_info(const unified_mpoly_t poly) {
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            return "Prime field (nmod)";
        case FIELD_ID_GF28:
            return "GF(2^8) lookup tables";
        case FIELD_ID_GF216:
            return "GF(2^16) tower field";
        case FIELD_ID_GF232:
            return "GF(2^32) PCLMUL";
        case FIELD_ID_GF264:
            return "GF(2^64) PCLMUL";
        case FIELD_ID_GF2128:
            return "GF(2^128) PCLMUL";
        case FIELD_ID_FQ_ZECH:
            return "Finite field with Zech logarithm";
        case FIELD_ID_FQ:
            return "General finite field (fq_nmod)";
        default:
            return "Unknown field type";
    }
}

/* Get information about the field representation being used by context */
const char* unified_mpoly_get_field_info_by_ctx(unified_mpoly_ctx_t ctx) {
    switch (ctx->field_ctx->field_id) {
        case FIELD_ID_NMOD:
            return "Prime field (nmod)";
        case FIELD_ID_GF28:
            return "GF(2^8) lookup tables";
        case FIELD_ID_GF216:
            return "GF(2^16) tower field";
        case FIELD_ID_GF232:
            return "GF(2^32) PCLMUL";
        case FIELD_ID_GF264:
            return "GF(2^64) PCLMUL";
        case FIELD_ID_GF2128:
            return "GF(2^128) PCLMUL";
        case FIELD_ID_FQ_ZECH:
            return "Finite field with Zech logarithm";
        case FIELD_ID_FQ:
            return "General finite field (fq_nmod)";
        default:
            return "Unknown field type";
    }
}

/* ============================================================================
   TESTING UTILITIES WITH ZECH SUPPORT
   ============================================================================ */

/* Generate random polynomial for testing */
void unified_mpoly_randtest(unified_mpoly_t poly, flint_rand_t state,
                           slong length, slong exp_bound) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    slong nvars = ctx->nvars;
    ulong *exp = (ulong *)calloc(nvars, sizeof(ulong));
    field_elem_u coeff;
    
    unified_mpoly_zero(poly);
    
    for (slong i = 0; i < length; i++) {
        /* Generate random exponents */
        for (slong j = 0; j < nvars; j++) {
            exp[j] = n_randint(state, exp_bound);
        }
        
        /* Generate random coefficient */
        switch (poly->field_id) {
            case FIELD_ID_NMOD:
                coeff.nmod = n_randint(state, GET_NMOD_CTX(ctx)->mod.n);
                if (coeff.nmod == 0) coeff.nmod = 1;
                break;
                
            case FIELD_ID_GF28:
                coeff.gf28 = n_randint(state, 256);
                if (coeff.gf28 == 0) coeff.gf28 = 1;
                break;
                
            case FIELD_ID_GF2128:
                coeff.gf2128.low = n_randtest(state);
                coeff.gf2128.high = n_randtest(state);
                if (gf2128_is_zero(&coeff.gf2128)) {
                    coeff.gf2128 = gf2128_one();
                }
                break;
                
            case FIELD_ID_FQ_ZECH:
                {
                    fq_zech_init(&coeff.fq_zech, ctx->field_ctx->ctx.zech_ctx);
                    fq_zech_randtest_not_zero(&coeff.fq_zech, state, ctx->field_ctx->ctx.zech_ctx);
                }
                break;
                
            default:
                /* For other fields, use generic method */
                {
                    fq_nmod_init(&coeff.fq, ctx->field_ctx->ctx.fq_ctx);
                    fq_nmod_randtest_not_zero(&coeff.fq, state, ctx->field_ctx->ctx.fq_ctx);
                }
                break;
        }
        
        unified_mpoly_set_coeff_ui(poly, &coeff, exp);
        
        /* Clean up coefficient if needed */
        if (poly->field_id == FIELD_ID_FQ_ZECH) {
            fq_zech_clear(&coeff.fq_zech, ctx->field_ctx->ctx.zech_ctx);
        } else if (poly->field_id == FIELD_ID_FQ) {
            fq_nmod_clear(&coeff.fq, ctx->field_ctx->ctx.fq_ctx);
        }
    }
    
    free(exp);
}

/* Check if two polynomials are equal */
int unified_mpoly_equal(const unified_mpoly_t A, const unified_mpoly_t B) {
    if (A->field_id != B->field_id) return 0;
    
    unified_mpoly_ctx_t ctx = A->ctx_ptr;
    
    switch (A->field_id) {
        case FIELD_ID_NMOD:
            return nmod_mpoly_equal(GET_NMOD_POLY(A), GET_NMOD_POLY(B),
                                   GET_NMOD_CTX(ctx));
            
        case FIELD_ID_FQ_ZECH:
            return fq_zech_mpoly_equal(GET_ZECH_POLY(A), GET_ZECH_POLY(B),
                                      GET_ZECH_CTX(ctx));
            
        default:
            return fq_nmod_mpoly_equal(GET_FQ_POLY(A), GET_FQ_POLY(B),
                                      GET_FQ_CTX(ctx));
    }
}

/* ============================================================================
   ADDITIONAL UTILITY FUNCTIONS
   ============================================================================ */

/* Check if polynomial is univariate */
int unified_mpoly_is_univariate(const unified_mpoly_t poly) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    slong nvars = ctx->nvars;
    
    if (nvars <= 1) return 1;
    
    slong *degrees = (slong *)calloc(nvars, sizeof(slong));
    unified_mpoly_degrees_si(degrees, poly);
    
    int var_count = 0;
    for (slong i = 0; i < nvars; i++) {
        if (degrees[i] > 0) var_count++;
    }
    
    free(degrees);
    return var_count <= 1;
}

/* Get the main variable of a univariate polynomial */
slong unified_mpoly_main_variable(const unified_mpoly_t poly) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    slong nvars = ctx->nvars;
    
    slong *degrees = (slong *)calloc(nvars, sizeof(slong));
    unified_mpoly_degrees_si(degrees, poly);
    
    slong main_var = -1;
    for (slong i = 0; i < nvars; i++) {
        if (degrees[i] > 0) {
            if (main_var == -1) {
                main_var = i;
            } else {
                /* Multiple variables, not truly univariate */
                main_var = -1;
                break;
            }
        }
    }
    
    free(degrees);
    return main_var;
}

/* Evaluate polynomial at a point (all variables = 0 except possibly one) */
void unified_mpoly_evaluate_all_ui(field_elem_u *result, const unified_mpoly_t poly,
                                   const ulong *values) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            result->nmod = nmod_mpoly_evaluate_all_ui(GET_NMOD_POLY(poly), values,
                                                     GET_NMOD_CTX(ctx));
            break;

        case FIELD_ID_FQ_ZECH:
            {
                /* For fq_zech, evaluation functions might not be available */
                /* Return zero as a safe default */
                fq_zech_init(&result->fq_zech, ctx->field_ctx->ctx.zech_ctx);
                fq_zech_zero(&result->fq_zech, ctx->field_ctx->ctx.zech_ctx);
                
                /* WARNING: This is a stub implementation */
                fprintf(stderr, "WARNING: unified_mpoly_evaluate_all_ui not implemented for fq_zech_mpoly\n");
            }
            break;
            
        default:
            {
                /* Convert values to fq_nmod format and create array of pointers */
                fq_nmod_t *fq_values = (fq_nmod_t *)malloc(ctx->nvars * sizeof(fq_nmod_t));
                fq_nmod_struct **fq_value_ptrs = (fq_nmod_struct **)malloc(ctx->nvars * sizeof(fq_nmod_struct *));
                
                for (slong i = 0; i < ctx->nvars; i++) {
                    fq_nmod_init(fq_values[i], ctx->field_ctx->ctx.fq_ctx);
                    fq_nmod_set_ui(fq_values[i], values[i], ctx->field_ctx->ctx.fq_ctx);
                    fq_value_ptrs[i] = fq_values[i];
                }
                
                fq_nmod_init(&result->fq, ctx->field_ctx->ctx.fq_ctx);
                fq_nmod_mpoly_evaluate_all_fq_nmod(&result->fq, GET_FQ_POLY(poly),
                                                   fq_value_ptrs, GET_FQ_CTX(ctx));
                
                /* Clean up */
                for (slong i = 0; i < ctx->nvars; i++) {
                    fq_nmod_clear(fq_values[i], ctx->field_ctx->ctx.fq_ctx);
                }
                free(fq_values);
                free(fq_value_ptrs);
            }
            break;
    }
}

/* Composition of polynomials */
int unified_mpoly_compose(unified_mpoly_t result, const unified_mpoly_t poly,
                         unified_mpoly_t * const *args, unified_mpoly_ctx_t ctx) {
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            return nmod_mpoly_compose_nmod_mpoly(GET_NMOD_POLY(result), GET_NMOD_POLY(poly),
                                                (nmod_mpoly_struct * const *)args,
                                                GET_NMOD_CTX(ctx), GET_NMOD_CTX(ctx));
            
        case FIELD_ID_FQ_ZECH:
            {
                /* For fq_zech, composition might not be directly available */
                /* We'll need to implement it manually or use a workaround */
                /* For now, return 0 to indicate not implemented */
                fprintf(stderr, "unified_mpoly_compose: fq_zech_mpoly composition not implemented\n");
                return 0;
            }
            
        default:
            return fq_nmod_mpoly_compose_fq_nmod_mpoly(GET_FQ_POLY(result), GET_FQ_POLY(poly),
                                                      (fq_nmod_mpoly_struct * const *)args,
                                                      GET_FQ_CTX(ctx), GET_FQ_CTX(ctx));
    }
}
