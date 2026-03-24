/* fq_unified_interface.c - Unified field operations implementation */

#include "fq_unified_interface.h"
extern int g_field_equation_reduction;

/* ============================================================================
   GLOBAL WORKSPACE
   ============================================================================ */

/* Per-thread workspace */
__thread unified_workspace_t g_unified_workspace = {0};

/* ============================================================================
   CONTEXT MANAGEMENT IMPLEMENTATIONS
   ============================================================================ */

/* Try to create Zech context from fq_nmod context */
int try_create_zech_context(fq_zech_ctx_struct *zech_ctx, const fq_nmod_ctx_t fq_ctx) {
    const nmod_poly_struct *modulus = fq_nmod_ctx_modulus(fq_ctx);
    
    /* Try to initialize - fq_zech_ctx_init_modulus returns void, but we can check validity */
    fq_zech_ctx_init_modulus(zech_ctx, modulus, "x");
    
    /* Check if initialization was successful by verifying the context */
    /* If qm1 is 0, initialization likely failed */
    return (zech_ctx->qm1 > 0);
}

/* Enhanced field context initialization with Zech support */
void field_ctx_init_enhanced(field_ctx_t *ctx, const fq_nmod_ctx_t fq_ctx, ulong zech_limit) {
    slong degree = fq_nmod_ctx_degree(fq_ctx);
    ulong prime = fq_nmod_ctx_prime(fq_ctx);
    
    ctx->zech_size_limit = zech_limit;
    
    if (degree == 1) {
        /* Prime field */
        ctx->field_id = FIELD_ID_NMOD;
        nmod_init(&ctx->ctx.nmod_ctx, prime);
        ctx->elem_size = sizeof(ulong);
        ctx->description = "Prime field (nmod)";
    } else if (prime == 2) {
        /* Check for optimized GF(2^n) implementations */
        switch (degree) {
            case 8:
                ctx->field_id = FIELD_ID_GF28;
                ctx->ctx.fq_ctx = fq_ctx;
                ctx->elem_size = sizeof(uint8_t);
                init_gf28_standard();
                init_gf28_conversion(fq_ctx);
                ctx->description = "GF(2^8) lookup tables";
                break;
            case 16:
                ctx->field_id = FIELD_ID_GF216;
                ctx->ctx.fq_ctx = fq_ctx;
                ctx->elem_size = sizeof(uint16_t);
                init_gf216_standard();
                init_gf216_conversion(fq_ctx);
                ctx->description = "GF(2^16) tower field";
                break;
            case 32:
                ctx->field_id = FIELD_ID_GF232;
                ctx->ctx.fq_ctx = fq_ctx;
                ctx->elem_size = sizeof(gf232_t);
                init_gf232();
                init_gf232_conversion(fq_ctx);
                ctx->description = "GF(2^32) PCLMUL";
                break;
            case 64:
                ctx->field_id = FIELD_ID_GF264;
                ctx->ctx.fq_ctx = fq_ctx;
                ctx->elem_size = sizeof(gf264_t);
                init_gf264();
                init_gf264_conversion(fq_ctx);
                ctx->description = "GF(2^64) PCLMUL";
                break;
            case 128:
                ctx->field_id = FIELD_ID_GF2128;
                ctx->ctx.fq_ctx = fq_ctx;
                ctx->elem_size = sizeof(gf2128_t);
                init_gf2128();
                init_gf2128_conversion(fq_ctx);
                ctx->description = "GF(2^128) PCLMUL";
                break;
            default:
                goto check_zech;
        }
    } else {
check_zech:
        /* Check if Zech logarithm is suitable */
        if (should_use_zech_logarithm(fq_ctx, zech_limit)) {
            fq_zech_ctx_struct *zech_ctx = (fq_zech_ctx_struct*)malloc(sizeof(fq_zech_ctx_struct));
            
            if (zech_ctx && try_create_zech_context(zech_ctx, fq_ctx)) {
                /* Successfully created Zech context */
                ctx->field_id = FIELD_ID_FQ_ZECH;
                ctx->ctx.zech_ctx = zech_ctx;
                ctx->elem_size = sizeof(fq_zech_struct);
                ctx->description = "Finite field with Zech logarithm";
                
                //printf("Auto-selected Zech logarithm for field of size %lu\n",
                //       calculate_field_size(prime, degree));
            } else {
                /* Failed to create Zech context, use fq_nmod */
                if (zech_ctx) {
                    /* Try to clear if partially initialized */
                    if (zech_ctx->qm1 > 0) {
                        fq_zech_ctx_clear(zech_ctx);
                    }
                    free(zech_ctx);
                }
                goto use_fq_nmod;
            }
        } else {
use_fq_nmod:
            /* Use general fq_nmod representation */
            ctx->field_id = FIELD_ID_FQ;
            ctx->ctx.fq_ctx = fq_ctx;
            ctx->elem_size = sizeof(fq_nmod_struct);
            ctx->description = "General finite field";
        }
    }
}

/* Standard field context initialization (backward compatibility) */
void field_ctx_init(field_ctx_t *ctx, const fq_nmod_ctx_t fq_ctx) {
    field_ctx_init_enhanced(ctx, fq_ctx, 1L << 20);  /* Default Zech limit: 2^20 */
}

void field_ctx_set_zech_limit(field_ctx_t *ctx, ulong limit) {
    ctx->zech_size_limit = limit;
}

const char* field_ctx_get_description(const field_ctx_t *ctx) {
    return ctx->description;
}

int field_ctx_is_using_zech(const field_ctx_t *ctx) {
    return ctx->field_id == FIELD_ID_FQ_ZECH;
}

ulong field_ctx_get_size(const field_ctx_t *ctx) {
    switch (ctx->field_id) {
        case FIELD_ID_NMOD:
            return ctx->ctx.nmod_ctx.n;
        case FIELD_ID_GF28:
            return 256;
        case FIELD_ID_GF216:
            return 65536;
        case FIELD_ID_GF232:
            return 1UL << 32;
        case FIELD_ID_GF264:
            return 0; /* Too large for ulong */
        case FIELD_ID_GF2128:
            return 0; /* Too large for ulong */
        case FIELD_ID_FQ_ZECH:
            return ctx->ctx.zech_ctx->qm1 + 1;
        case FIELD_ID_FQ:
            {
                ulong p = fq_nmod_ctx_prime(ctx->ctx.fq_ctx);
                slong d = fq_nmod_ctx_degree(ctx->ctx.fq_ctx);
                return calculate_field_size(p, d);
            }
        default:
            return 0;
    }
}

void field_ctx_clear(field_ctx_t *ctx) {
    if (ctx->field_id == FIELD_ID_FQ_ZECH && ctx->ctx.zech_ctx) {
        fq_zech_ctx_clear(ctx->ctx.zech_ctx);
        free(ctx->ctx.zech_ctx);
        ctx->ctx.zech_ctx = NULL;
    }
    /* Other field types don't own their contexts */
}

/* ============================================================================
   FIELD OPERATIONS IMPLEMENTATION
   ============================================================================ */

void field_neg(field_elem_u *res, const field_elem_u *a,
               field_id_t field_id, const void *ctx) {
    switch (field_id) {
        case FIELD_ID_GF28:
        case FIELD_ID_GF216:
        case FIELD_ID_GF232:
        case FIELD_ID_GF264:
        case FIELD_ID_GF2128:
            *res = *a;  /* In GF(2^n), -a = a */
            break;
        case FIELD_ID_NMOD:
            res->nmod = nmod_neg(a->nmod, *(const nmod_t*)ctx);
            break;
        case FIELD_ID_FQ_ZECH:
            if (res != a) {
                fq_zech_init(&res->fq_zech, (const fq_zech_ctx_struct *)ctx);
            }
            fq_zech_neg(&res->fq_zech, &a->fq_zech, (const fq_zech_ctx_struct *)ctx);
            break;
        case FIELD_ID_FQ:
            fq_nmod_neg(&res->fq, &a->fq, (const fq_nmod_ctx_struct *)ctx);
            break;
    }
}

void field_inv(field_elem_u *res, const field_elem_u *a,
               field_id_t field_id, const void *ctx) {
    switch (field_id) {
        case FIELD_ID_GF28:
            res->gf28 = gf28_inv(a->gf28);
            break;
        case FIELD_ID_GF216:
            res->gf216 = gf216_inv(a->gf216);
            break;
        case FIELD_ID_GF232:
            res->gf232 = gf232_inv(&a->gf232);
            break;
        case FIELD_ID_GF264:
            res->gf264 = gf264_inv(&a->gf264);
            break;
        case FIELD_ID_GF2128:
            res->gf2128 = gf2128_inv(&a->gf2128);
            break;
        case FIELD_ID_NMOD:
            res->nmod = n_invmod(a->nmod, ((const nmod_t*)ctx)->n);
            break;
        case FIELD_ID_FQ_ZECH:
            if (res != a) {
                fq_zech_init(&res->fq_zech, (const fq_zech_ctx_struct *)ctx);
            }
            fq_zech_inv(&res->fq_zech, &a->fq_zech, (const fq_zech_ctx_struct *)ctx);
            break;
        case FIELD_ID_FQ:
            fq_nmod_inv(&res->fq, &a->fq, (const fq_nmod_ctx_struct *)ctx);
            break;
    }
}

void field_set_zero(field_elem_u *res, field_id_t field_id, const void *ctx) {
    switch (field_id) {
        case FIELD_ID_GF28:
            res->gf28 = 0;
            break;
        case FIELD_ID_GF216:
            res->gf216 = 0;
            break;
        case FIELD_ID_GF232:
            res->gf232 = gf232_zero();
            break;
        case FIELD_ID_GF264:
            res->gf264 = gf264_zero();
            break;
        case FIELD_ID_GF2128:
            res->gf2128 = gf2128_zero();
            break;
        case FIELD_ID_NMOD:
            res->nmod = 0;
            break;
        case FIELD_ID_FQ_ZECH:
            fq_zech_zero(&res->fq_zech, (const fq_zech_ctx_struct *)ctx);
            break;
        case FIELD_ID_FQ:
            fq_nmod_zero(&res->fq, (const fq_nmod_ctx_struct *)ctx);
            break;
    }
}

void field_set_one(field_elem_u *res, field_id_t field_id, const void *ctx) {
    switch (field_id) {
        case FIELD_ID_GF28:
            res->gf28 = 1;
            break;
        case FIELD_ID_GF216:
            res->gf216 = 1;
            break;
        case FIELD_ID_GF232:
            res->gf232 = gf232_one();
            break;
        case FIELD_ID_GF264:
            res->gf264 = gf264_one();
            break;
        case FIELD_ID_GF2128:
            res->gf2128 = gf2128_one();
            break;
        case FIELD_ID_NMOD:
            res->nmod = 1;
            break;
        case FIELD_ID_FQ_ZECH:
            fq_zech_one(&res->fq_zech, (const fq_zech_ctx_struct *)ctx);
            break;
        case FIELD_ID_FQ:
            fq_nmod_one(&res->fq, (const fq_nmod_ctx_struct *)ctx);
            break;
    }
}

int field_equal(const field_elem_u *a, const field_elem_u *b, 
                field_id_t field_id, const void *ctx) {
    switch (field_id) {
        case FIELD_ID_GF28:
            return a->gf28 == b->gf28;
        case FIELD_ID_GF216:
            return a->gf216 == b->gf216;
        case FIELD_ID_GF232:
            return a->gf232.value == b->gf232.value;
        case FIELD_ID_GF264:
            return a->gf264.value == b->gf264.value;
        case FIELD_ID_GF2128:
            return a->gf2128.low == b->gf2128.low && a->gf2128.high == b->gf2128.high;
        case FIELD_ID_NMOD:
            return a->nmod == b->nmod;
        case FIELD_ID_FQ_ZECH:
            return fq_zech_equal(&a->fq_zech, &b->fq_zech, (const fq_zech_ctx_struct *)ctx);
        case FIELD_ID_FQ:
            return fq_nmod_equal(&a->fq, &b->fq, (const fq_nmod_ctx_struct *)ctx);
    }
    return 0;
}

void field_init_elem(field_elem_u *elem, field_id_t field_id, const void *ctx) {
    switch (field_id) {
        case FIELD_ID_GF28:
            elem->gf28 = 0;
            break;
        case FIELD_ID_GF216:
            elem->gf216 = 0;
            break;
        case FIELD_ID_GF232:
            elem->gf232 = gf232_zero();
            break;
        case FIELD_ID_GF264:
            elem->gf264 = gf264_zero();
            break;
        case FIELD_ID_GF2128:
            elem->gf2128 = gf2128_zero();
            break;
        case FIELD_ID_NMOD:
            elem->nmod = 0;
            break;
        case FIELD_ID_FQ_ZECH:
            fq_zech_init(&elem->fq_zech, (const fq_zech_ctx_struct *)ctx);
            break;
        case FIELD_ID_FQ:
            fq_nmod_init(&elem->fq, (const fq_nmod_ctx_struct *)ctx);
            break;
    }
}

void field_clear_elem(field_elem_u *elem, field_id_t field_id, const void *ctx) {
    switch (field_id) {
        case FIELD_ID_FQ_ZECH:
            fq_zech_clear(&elem->fq_zech, (const fq_zech_ctx_struct *)ctx);
            break;
        case FIELD_ID_FQ:
            fq_nmod_clear(&elem->fq, (const fq_nmod_ctx_struct *)ctx);
            break;
        default:
            /* No cleanup needed for other field types */
            break;
    }
}

void field_set_elem(field_elem_u *res, const field_elem_u *a,
                    field_id_t field_id, const void *ctx) {
    switch (field_id) {
        case FIELD_ID_GF28:
            res->gf28 = a->gf28;
            break;
        case FIELD_ID_GF216:
            res->gf216 = a->gf216;
            break;
        case FIELD_ID_GF232:
            res->gf232 = a->gf232;
            break;
        case FIELD_ID_GF264:
            res->gf264 = a->gf264;
            break;
        case FIELD_ID_GF2128:
            res->gf2128 = a->gf2128;
            break;
        case FIELD_ID_NMOD:
            res->nmod = a->nmod;
            break;
        case FIELD_ID_FQ_ZECH:
            fq_zech_set(&res->fq_zech, &a->fq_zech, (const fq_zech_ctx_struct *)ctx);
            break;
        case FIELD_ID_FQ:
            fq_nmod_set(&res->fq, &a->fq, (const fq_nmod_ctx_struct *)ctx);
            break;
    }
}

void field_sub(field_elem_u *res, const field_elem_u *a, const field_elem_u *b,
               field_id_t field_id, const void *ctx) {
    switch (field_id) {
        case FIELD_ID_GF28:
        case FIELD_ID_GF216:
        case FIELD_ID_GF232:
        case FIELD_ID_GF264:
        case FIELD_ID_GF2128:
            /* In GF(2^n), subtraction equals addition */
            field_add(res, a, b, field_id, ctx);
            break;
        case FIELD_ID_NMOD:
            res->nmod = nmod_sub(a->nmod, b->nmod, *(const nmod_t*)ctx);
            break;
        case FIELD_ID_FQ_ZECH:
            if (res != a && res != b) {
                fq_zech_init(&res->fq_zech, (const fq_zech_ctx_struct *)ctx);
            }
            fq_zech_sub(&res->fq_zech, &a->fq_zech, &b->fq_zech, (const fq_zech_ctx_struct *)ctx);
            break;
        case FIELD_ID_FQ:
            if (res != a && res != b) {
                fq_nmod_init(&res->fq, (const fq_nmod_ctx_struct *)ctx);
            }
            fq_nmod_sub(&res->fq, &a->fq, &b->fq, (const fq_nmod_ctx_struct *)ctx);
            break;
    }
}

/* ============================================================================
   CONVERSION FUNCTIONS WITH ZECH SUPPORT
   ============================================================================ */

void fq_nmod_to_field_elem(field_elem_u *res, const fq_nmod_t elem, 
                          const field_ctx_t *ctx) {
    switch (ctx->field_id) {
        case FIELD_ID_NMOD:
            res->nmod = nmod_poly_get_coeff_ui(elem, 0);
            break;
        case FIELD_ID_GF28:
            res->gf28 = fq_nmod_to_gf28_elem(elem, ctx->ctx.fq_ctx);
            break;
        case FIELD_ID_GF216:
            res->gf216 = fq_nmod_to_gf216_elem(elem, ctx->ctx.fq_ctx);
            break;
        case FIELD_ID_GF232:
            res->gf232 = fq_nmod_to_gf232(elem, ctx->ctx.fq_ctx);
            break;
        case FIELD_ID_GF264:
            res->gf264 = fq_nmod_to_gf264(elem, ctx->ctx.fq_ctx);
            break;
        case FIELD_ID_GF2128:
            res->gf2128 = fq_nmod_to_gf2128(elem, ctx->ctx.fq_ctx);
            break;
        case FIELD_ID_FQ_ZECH:
            fq_zech_init(&res->fq_zech, ctx->ctx.zech_ctx);
            fq_zech_set_fq_nmod(&res->fq_zech, elem, ctx->ctx.zech_ctx);
            break;
        case FIELD_ID_FQ:
            fq_nmod_init(&res->fq, ctx->ctx.fq_ctx);
            fq_nmod_set(&res->fq, elem, ctx->ctx.fq_ctx);
            break;
    }
}

void field_elem_to_fq_nmod(fq_nmod_t res, const field_elem_u *elem,
                          const field_ctx_t *ctx) {
    switch (ctx->field_id) {
        case FIELD_ID_NMOD:
            fq_nmod_zero(res, ctx->ctx.fq_ctx);
            nmod_poly_set_coeff_ui(res, 0, elem->nmod);
            break;
        case FIELD_ID_GF28:
            gf28_elem_to_fq_nmod(res, elem->gf28, ctx->ctx.fq_ctx);
            break;
        case FIELD_ID_GF216:
            gf216_elem_to_fq_nmod(res, elem->gf216, ctx->ctx.fq_ctx);
            break;
        case FIELD_ID_GF232:
            gf232_to_fq_nmod(res, &elem->gf232, ctx->ctx.fq_ctx);
            break;
        case FIELD_ID_GF264:
            gf264_to_fq_nmod(res, &elem->gf264, ctx->ctx.fq_ctx);
            break;
        case FIELD_ID_GF2128:
            gf2128_to_fq_nmod(res, &elem->gf2128, ctx->ctx.fq_ctx);
            break;
        case FIELD_ID_FQ_ZECH:
            fq_zech_get_fq_nmod(res, &elem->fq_zech, ctx->ctx.zech_ctx);
            break;
        case FIELD_ID_FQ:
            fq_nmod_set(res, &elem->fq, ctx->ctx.fq_ctx);
            break;
    }
}

/* ============================================================================
   POLYNOMIAL OPERATIONS IMPLEMENTATION
   ============================================================================ */

void unified_poly_init(unified_poly_t poly, field_ctx_t *ctx) {
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
    poly->ctx = ctx;
}

void unified_poly_clear(unified_poly_t poly) {
    if (poly->coeffs) {
        void *ctx_ptr = NULL;
        switch (poly->ctx->field_id) {
            case FIELD_ID_NMOD:
                ctx_ptr = (void*)&poly->ctx->ctx.nmod_ctx;
                break;
            case FIELD_ID_FQ_ZECH:
                ctx_ptr = (void*)poly->ctx->ctx.zech_ctx;
                break;
            default:
                ctx_ptr = (void*)poly->ctx->ctx.fq_ctx;
                break;
        }
        
        if (poly->ctx->field_id == FIELD_ID_FQ || poly->ctx->field_id == FIELD_ID_FQ_ZECH) {
            for (slong i = 0; i < poly->alloc; i++) {
                field_clear_elem(&poly->coeffs[i], poly->ctx->field_id, ctx_ptr);
            }
        }
        free(poly->coeffs);
        poly->coeffs = NULL;
    }
    poly->length = 0;
    poly->alloc = 0;
}

void unified_poly_fit_length(unified_poly_t poly, slong len) {
    if (len > poly->alloc) {
        slong new_alloc = poly->alloc;
        if (new_alloc == 0) new_alloc = 16;
        while (new_alloc < len) new_alloc *= 2;
        
        void *ctx_ptr = NULL;
        switch (poly->ctx->field_id) {
            case FIELD_ID_NMOD:
                ctx_ptr = (void*)&poly->ctx->ctx.nmod_ctx;
                break;
            case FIELD_ID_FQ_ZECH:
                ctx_ptr = (void*)poly->ctx->ctx.zech_ctx;
                break;
            default:
                ctx_ptr = (void*)poly->ctx->ctx.fq_ctx;
                break;
        }
        
        /* For FIELD_ID_FQ and FIELD_ID_FQ_ZECH, we need to be careful about initialization */
        if (poly->ctx->field_id == FIELD_ID_FQ || poly->ctx->field_id == FIELD_ID_FQ_ZECH) {
            /* Allocate new array */
            field_elem_u *new_coeffs = (field_elem_u *)malloc(new_alloc * sizeof(field_elem_u));
            
            /* Copy existing coefficients */
            for (slong i = 0; i < poly->length; i++) {
                /* Initialize new location */
                field_init_elem(&new_coeffs[i], poly->ctx->field_id, ctx_ptr);
                /* Copy value */
                field_set_elem(&new_coeffs[i], &poly->coeffs[i], poly->ctx->field_id, ctx_ptr);
            }
            
            /* Initialize remaining elements */
            for (slong i = poly->length; i < new_alloc; i++) {
                field_init_elem(&new_coeffs[i], poly->ctx->field_id, ctx_ptr);
            }
            
            /* Clear old coefficients */
            if (poly->coeffs) {
                for (slong i = 0; i < poly->alloc; i++) {
                    field_clear_elem(&poly->coeffs[i], poly->ctx->field_id, ctx_ptr);
                }
                free(poly->coeffs);
            }
            
            poly->coeffs = new_coeffs;
        } else {
            /* For other field types, use realloc */
            field_elem_u *new_coeffs = (field_elem_u *)realloc(poly->coeffs, 
                                                               new_alloc * sizeof(field_elem_u));
            if (!new_coeffs) {
                printf("Memory allocation failed in unified_poly_fit_length\n");
                return;
            }
            poly->coeffs = new_coeffs;
            
            /* Initialize new elements */
            for (slong i = poly->alloc; i < new_alloc; i++) {
                field_init_elem(&poly->coeffs[i], poly->ctx->field_id, ctx_ptr);
            }
        }
        
        poly->alloc = new_alloc;
    }
}

void unified_poly_normalise(unified_poly_t poly) {
    void *ctx_ptr = NULL;
    switch (poly->ctx->field_id) {
        case FIELD_ID_NMOD:
            ctx_ptr = (void*)&poly->ctx->ctx.nmod_ctx;
            break;
        case FIELD_ID_FQ_ZECH:
            ctx_ptr = (void*)poly->ctx->ctx.zech_ctx;
            break;
        default:
            ctx_ptr = (void*)poly->ctx->ctx.fq_ctx;
            break;
    }
    
    /* Remove all trailing zeros */
    while (poly->length > 0 && 
           field_is_zero(&poly->coeffs[poly->length - 1], poly->ctx->field_id, ctx_ptr)) {
        poly->length--;
    }
}

void unified_poly_zero(unified_poly_t poly) {
    poly->length = 0;
}

int unified_poly_is_zero(const unified_poly_t poly) {
    return poly->length == 0;
}

slong unified_poly_degree(const unified_poly_t poly) {
    return poly->length - 1;
}

void unified_poly_set(unified_poly_t res, const unified_poly_t poly) {
    if (res == poly) return;
    
    unified_poly_fit_length(res, poly->length);
    res->length = poly->length;
    
    void *ctx_ptr = NULL;
    switch (res->ctx->field_id) {
        case FIELD_ID_NMOD:
            ctx_ptr = (void*)&res->ctx->ctx.nmod_ctx;
            break;
        case FIELD_ID_FQ_ZECH:
            ctx_ptr = (void*)res->ctx->ctx.zech_ctx;
            break;
        default:
            ctx_ptr = (void*)res->ctx->ctx.fq_ctx;
            break;
    }
    
    /* Use field_set_elem for all copies */
    for (slong i = 0; i < poly->length; i++) {
        field_set_elem(&res->coeffs[i], &poly->coeffs[i], res->ctx->field_id, ctx_ptr);
    }
}

void unified_poly_add(unified_poly_t res, const unified_poly_t a, const unified_poly_t b) {
    slong max_len = FLINT_MAX(a->length, b->length);
    slong min_len = FLINT_MIN(a->length, b->length);
    
    if (max_len == 0) {
        unified_poly_zero(res);
        return;
    }
    
    unified_poly_fit_length(res, max_len);
    
    void *ctx_ptr = NULL;
    switch (res->ctx->field_id) {
        case FIELD_ID_NMOD:
            ctx_ptr = (void*)&res->ctx->ctx.nmod_ctx;
            break;
        case FIELD_ID_FQ_ZECH:
            ctx_ptr = (void*)res->ctx->ctx.zech_ctx;
            break;
        default:
            ctx_ptr = (void*)res->ctx->ctx.fq_ctx;
            break;
    }
    
    /* Add common coefficients */
    for (slong i = 0; i < min_len; i++) {
        field_add(&res->coeffs[i], &a->coeffs[i], &b->coeffs[i], 
                 res->ctx->field_id, ctx_ptr);
    }
    
    /* Copy remaining coefficients */
    if (a->length > b->length) {
        for (slong i = min_len; i < a->length; i++) {
            field_set_elem(&res->coeffs[i], &a->coeffs[i], res->ctx->field_id, ctx_ptr);
        }
        res->length = a->length;
    } else if (b->length > a->length) {
        for (slong i = min_len; i < b->length; i++) {
            field_set_elem(&res->coeffs[i], &b->coeffs[i], res->ctx->field_id, ctx_ptr);
        }
        res->length = b->length;
    } else {
        res->length = min_len;
    }
    
    unified_poly_normalise(res);
}

void unified_poly_scalar_mul(unified_poly_t res, const unified_poly_t poly, 
                            const field_elem_u *c) {
    void *ctx_ptr = NULL;
    switch (poly->ctx->field_id) {
        case FIELD_ID_NMOD:
            ctx_ptr = (void*)&poly->ctx->ctx.nmod_ctx;
            break;
        case FIELD_ID_FQ_ZECH:
            ctx_ptr = (void*)poly->ctx->ctx.zech_ctx;
            break;
        default:
            ctx_ptr = (void*)poly->ctx->ctx.fq_ctx;
            break;
    }
    
    if (field_is_zero(c, poly->ctx->field_id, ctx_ptr)) {
        unified_poly_zero(res);
        return;
    }
    
    if (field_is_one(c, poly->ctx->field_id, ctx_ptr)) {
        unified_poly_set(res, poly);
        return;
    }
    
    unified_poly_fit_length(res, poly->length);
    res->length = poly->length;
    
    for (slong i = 0; i < poly->length; i++) {
        field_mul(&res->coeffs[i], &poly->coeffs[i], c, poly->ctx->field_id, ctx_ptr);
    }
    
    unified_poly_normalise(res);
}

static inline slong unified_reduce_degree_field(slong e, ulong q) {
    if (e == 0 || (ulong)e < q) return e;
    return (slong)(((ulong)(e - 1) % (q - 1)) + 1);
}

static ulong unified_poly_field_size_q(const field_ctx_t *ctx) {
    if (!ctx) return 0;
    switch (ctx->field_id) {
        case FIELD_ID_NMOD:
            return ctx->ctx.nmod_ctx.n;
        case FIELD_ID_FQ_ZECH:
            if (ctx->ctx.zech_ctx) {
                mp_limb_t p = ctx->ctx.zech_ctx->fq_nmod_ctx->modulus->mod.n;
                slong d = fq_zech_ctx_degree(ctx->ctx.zech_ctx);
                ulong q = 1;
                for (slong i = 0; i < d; i++) {
                    if (q > WORD_MAX / p) return WORD_MAX;
                    q *= p;
                }
                return q;
            }
            return 0;
        default:
            if (ctx->ctx.fq_ctx) {
                mp_limb_t p = fq_nmod_ctx_prime(ctx->ctx.fq_ctx);
                slong d = fq_nmod_ctx_degree(ctx->ctx.fq_ctx);
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

static void unified_poly_reduce_field_equation_inplace(unified_poly_t poly) {
    if (!poly || poly->length <= 1 || !poly->ctx) return;
    ulong q = unified_poly_field_size_q(poly->ctx);
    if (q <= 1 || q == WORD_MAX) return;
    if ((ulong)(poly->length - 1) < q) return;

    slong max_keep = (q > (ulong)WORD_MAX) ? poly->length : FLINT_MIN(poly->length, (slong)q);
    if (max_keep <= 0) return;

    unified_poly_struct reduced;
    unified_poly_init(&reduced, poly->ctx);
    unified_poly_fit_length(&reduced, max_keep);

    void *ctx_ptr = NULL;
    switch (poly->ctx->field_id) {
        case FIELD_ID_NMOD:
            ctx_ptr = (void*)&poly->ctx->ctx.nmod_ctx;
            break;
        case FIELD_ID_FQ_ZECH:
            ctx_ptr = (void*)poly->ctx->ctx.zech_ctx;
            break;
        default:
            ctx_ptr = (void*)poly->ctx->ctx.fq_ctx;
            break;
    }

    for (slong i = 0; i < max_keep; i++) {
        field_set_zero(&reduced.coeffs[i], poly->ctx->field_id, ctx_ptr);
    }

    slong max_used = 0;
    for (slong i = 0; i < poly->length; i++) {
        slong tgt = unified_reduce_degree_field(i, q);
        field_add(&reduced.coeffs[tgt], &reduced.coeffs[tgt], &poly->coeffs[i],
                  poly->ctx->field_id, ctx_ptr);
        if (tgt > max_used) max_used = tgt;
    }

    reduced.length = max_used + 1;
    unified_poly_normalise(&reduced);
    unified_poly_set(poly, &reduced);
    unified_poly_clear(&reduced);
}

void unified_poly_mul(unified_poly_t res, const unified_poly_t a, const unified_poly_t b) {
    if (unified_poly_is_zero(a) || unified_poly_is_zero(b)) {
        unified_poly_zero(res);
        return;
    }
    
    slong rlen = a->length + b->length - 1;
    
    /* Always use temporary for FIELD_ID_FQ and FIELD_ID_FQ_ZECH to avoid aliasing issues */
    if (res->ctx->field_id == FIELD_ID_FQ || res->ctx->field_id == FIELD_ID_FQ_ZECH || res == a || res == b) {
        unified_poly_struct temp;
        unified_poly_init(&temp, res->ctx);
        unified_poly_fit_length(&temp, rlen);
        
        void *ctx_ptr = NULL;
        switch (res->ctx->field_id) {
            case FIELD_ID_NMOD:
                ctx_ptr = (void*)&res->ctx->ctx.nmod_ctx;
                break;
            case FIELD_ID_FQ_ZECH:
                ctx_ptr = (void*)res->ctx->ctx.zech_ctx;
                break;
            default:
                ctx_ptr = (void*)res->ctx->ctx.fq_ctx;
                break;
        }
        
        /* Initialize result coefficients to zero */
        for (slong i = 0; i < rlen; i++) {
            field_set_zero(&temp.coeffs[i], res->ctx->field_id, ctx_ptr);
        }
        
        /* Multiply */
        field_elem_u prod;
        field_init_elem(&prod, res->ctx->field_id, ctx_ptr);
        
        for (slong i = 0; i < a->length; i++) {
            for (slong j = 0; j < b->length; j++) {
                field_mul(&prod, &a->coeffs[i], &b->coeffs[j], res->ctx->field_id, ctx_ptr);
                field_add(&temp.coeffs[i + j], &temp.coeffs[i + j], &prod, 
                         res->ctx->field_id, ctx_ptr);
            }
        }
        
        field_clear_elem(&prod, res->ctx->field_id, ctx_ptr);
        
        temp.length = rlen;
        unified_poly_normalise(&temp);
        if (g_field_equation_reduction) {
            unified_poly_reduce_field_equation_inplace(&temp);
        }
        
        /* Copy result back */
        unified_poly_set(res, &temp);
        unified_poly_clear(&temp);
    } else {
        /* Original implementation for non-FQ fields */
        unified_poly_fit_length(res, rlen);
        
        void *ctx_ptr = NULL;
        switch (res->ctx->field_id) {
            case FIELD_ID_NMOD:
                ctx_ptr = (void*)&res->ctx->ctx.nmod_ctx;
                break;
            case FIELD_ID_FQ_ZECH:
                ctx_ptr = (void*)res->ctx->ctx.zech_ctx;
                break;
            default:
                ctx_ptr = (void*)res->ctx->ctx.fq_ctx;
                break;
        }
        
        /* Initialize result coefficients to zero */
        for (slong i = 0; i < rlen; i++) {
            field_set_zero(&res->coeffs[i], res->ctx->field_id, ctx_ptr);
        }
        
        /* Multiply */
        field_elem_u prod;
        field_init_elem(&prod, res->ctx->field_id, ctx_ptr);
        
        for (slong i = 0; i < a->length; i++) {
            for (slong j = 0; j < b->length; j++) {
                field_mul(&prod, &a->coeffs[i], &b->coeffs[j], res->ctx->field_id, ctx_ptr);
                field_add(&res->coeffs[i + j], &res->coeffs[i + j], &prod, 
                         res->ctx->field_id, ctx_ptr);
            }
        }
        
        field_clear_elem(&prod, res->ctx->field_id, ctx_ptr);
        
        res->length = rlen;
        unified_poly_normalise(res);
        if (g_field_equation_reduction) {
            unified_poly_reduce_field_equation_inplace(res);
        }
    }
}

void unified_poly_get_coeff(field_elem_u *coeff, const unified_poly_t poly, slong i) {
    void *ctx_ptr = NULL;
    switch (poly->ctx->field_id) {
        case FIELD_ID_NMOD:
            ctx_ptr = (void*)&poly->ctx->ctx.nmod_ctx;
            break;
        case FIELD_ID_FQ_ZECH:
            ctx_ptr = (void*)poly->ctx->ctx.zech_ctx;
            break;
        default:
            ctx_ptr = (void*)poly->ctx->ctx.fq_ctx;
            break;
    }
    
    if (i < poly->length) {
        field_set_elem(coeff, &poly->coeffs[i], poly->ctx->field_id, ctx_ptr);
    } else {
        field_set_zero(coeff, poly->ctx->field_id, ctx_ptr);
    }
}

void unified_poly_set_coeff(unified_poly_t poly, slong i, const field_elem_u *coeff) {
    unified_poly_fit_length(poly, i + 1);
    
    void *ctx_ptr = NULL;
    switch (poly->ctx->field_id) {
        case FIELD_ID_NMOD:
            ctx_ptr = (void*)&poly->ctx->ctx.nmod_ctx;
            break;
        case FIELD_ID_FQ_ZECH:
            ctx_ptr = (void*)poly->ctx->ctx.zech_ctx;
            break;
        default:
            ctx_ptr = (void*)poly->ctx->ctx.fq_ctx;
            break;
    }
    
    if (i >= poly->length) {
        /* Zero out intermediate coefficients */
        for (slong j = poly->length; j < i; j++) {
            field_set_zero(&poly->coeffs[j], poly->ctx->field_id, ctx_ptr);
        }
        poly->length = i + 1;
    }
    
    field_set_elem(&poly->coeffs[i], coeff, poly->ctx->field_id, ctx_ptr);
    
    /* Update length if setting leading coefficient to zero */
    if (i == poly->length - 1) {
        unified_poly_normalise(poly);
    }
}

void unified_poly_shift_left(unified_poly_t res, const unified_poly_t poly, slong n) {
    if (n == 0) {
        unified_poly_set(res, poly);
        return;
    }
    
    if (unified_poly_is_zero(poly)) {
        unified_poly_zero(res);
        return;
    }
    
    slong new_len = poly->length + n;
    unified_poly_fit_length(res, new_len);
    
    void *ctx_ptr = NULL;
    switch (res->ctx->field_id) {
        case FIELD_ID_NMOD:
            ctx_ptr = (void*)&res->ctx->ctx.nmod_ctx;
            break;
        case FIELD_ID_FQ_ZECH:
            ctx_ptr = (void*)res->ctx->ctx.zech_ctx;
            break;
        default:
            ctx_ptr = (void*)res->ctx->ctx.fq_ctx;
            break;
    }
    
    /* Zero lower coefficients */
    for (slong i = 0; i < n; i++) {
        field_set_zero(&res->coeffs[i], res->ctx->field_id, ctx_ptr);
    }
    
    /* Copy coefficients properly */
    if (res == poly) {
        /* In-place: shift from right to left */
        for (slong i = new_len - 1; i >= n; i--) {
            field_set_elem(&res->coeffs[i], &res->coeffs[i - n], res->ctx->field_id, ctx_ptr);
        }
    } else {
        /* Copy shifted */
        for (slong i = 0; i < poly->length; i++) {
            field_set_elem(&res->coeffs[i + n], &poly->coeffs[i], res->ctx->field_id, ctx_ptr);
        }
    }
    
    res->length = new_len;
}

/* ============================================================================
   CONVERSION FUNCTIONS
   ============================================================================ */

void fq_nmod_poly_to_unified(unified_poly_t res, const fq_nmod_poly_t poly,
                            const fq_nmod_ctx_t ctx, field_ctx_t *field_ctx) {
    slong len = fq_nmod_poly_length(poly, ctx);
    
    if (len == 0) {
        unified_poly_zero(res);
        return;
    }
    
    unified_poly_fit_length(res, len);
    res->length = len;
    
    for (slong i = 0; i < len; i++) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, ctx);
        fq_nmod_poly_get_coeff(coeff, poly, i, ctx);
        fq_nmod_to_field_elem(&res->coeffs[i], coeff, field_ctx);
        fq_nmod_clear(coeff, ctx);
    }
    
    unified_poly_normalise(res);
}

void unified_to_fq_nmod_poly(fq_nmod_poly_t res, const unified_poly_t poly,
                            const fq_nmod_ctx_t ctx, field_ctx_t *field_ctx) {
    fq_nmod_poly_zero(res, ctx);
    
    void *ctx_ptr = NULL;
    switch (field_ctx->field_id) {
        case FIELD_ID_NMOD:
            ctx_ptr = (void*)&field_ctx->ctx.nmod_ctx;
            break;
        case FIELD_ID_FQ_ZECH:
            ctx_ptr = (void*)field_ctx->ctx.zech_ctx;
            break;
        default:
            ctx_ptr = (void*)field_ctx->ctx.fq_ctx;
            break;
    }
    
    for (slong i = 0; i < poly->length; i++) {
        if (!field_is_zero(&poly->coeffs[i], field_ctx->field_id, ctx_ptr)) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            field_elem_to_fq_nmod(coeff, &poly->coeffs[i], field_ctx);
            fq_nmod_poly_set_coeff(res, i, coeff, ctx);
            fq_nmod_clear(coeff, ctx);
        }
    }
}

/* ============================================================================
   MATRIX OPERATIONS IMPLEMENTATION
   ============================================================================ */

void unified_poly_mat_init(unified_poly_mat_t mat, slong rows, slong cols,
                          field_ctx_t *ctx) {
    mat->entries = NULL;
    mat->rows = NULL;
    mat->ctx = ctx;
    
    if (rows > 0 && cols > 0) {
        mat->entries = (unified_poly_struct *)malloc(rows * cols * sizeof(unified_poly_struct));
        mat->rows = (unified_poly_struct **)malloc(rows * sizeof(unified_poly_struct *));
        
        for (slong i = 0; i < rows * cols; i++) {
            unified_poly_init(mat->entries + i, ctx);
        }
        
        for (slong i = 0; i < rows; i++) {
            mat->rows[i] = mat->entries + i * cols;
        }
    }
    
    mat->r = rows;
    mat->c = cols;
}

void unified_poly_mat_clear(unified_poly_mat_t mat) {
    if (mat->entries != NULL) {
        for (slong i = 0; i < mat->r * mat->c; i++) {
            unified_poly_clear(mat->entries + i);
        }
        free(mat->entries);
        free(mat->rows);
    }
}

unified_poly_struct *unified_poly_mat_entry(unified_poly_mat_t mat, slong i, slong j) {
    return mat->rows[i] + j;
}

const unified_poly_struct *unified_poly_mat_entry_const(const unified_poly_mat_t mat, slong i, slong j) {
    return mat->rows[i] + j;
}

void fq_nmod_poly_mat_to_unified(unified_poly_mat_t res, const fq_nmod_poly_mat_t mat,
                                const fq_nmod_ctx_t ctx, field_ctx_t *field_ctx) {
    for (slong i = 0; i < mat->r; i++) {
        for (slong j = 0; j < mat->c; j++) {
            fq_nmod_poly_to_unified(unified_poly_mat_entry(res, i, j),
                                   fq_nmod_poly_mat_entry(mat, i, j),
                                   ctx, field_ctx);
        }
    }
}

void unified_to_fq_nmod_poly_mat(fq_nmod_poly_mat_t res, const unified_poly_mat_t mat,
                                const fq_nmod_ctx_t ctx, field_ctx_t *field_ctx) {
    for (slong i = 0; i < mat->r; i++) {
        for (slong j = 0; j < mat->c; j++) {
            unified_to_fq_nmod_poly(fq_nmod_poly_mat_entry(res, i, j),
                                   unified_poly_mat_entry_const(mat, i, j),
                                   ctx, field_ctx);
        }
    }
}

/* ============================================================================
   WORKSPACE MANAGEMENT IMPLEMENTATION
   ============================================================================ */

/* Clear workspace when switching fields */
static void clear_workspace(unified_workspace_t *ws, field_ctx_t *ctx) {
    if (ws->initialized && ws->field_ctx) {
        void *ctx_ptr = ws->field_ctx;
        
        field_clear_elem(&ws->lc1, ws->field_id, ctx_ptr);
        field_clear_elem(&ws->lc2, ws->field_id, ctx_ptr);
        field_clear_elem(&ws->cst, ws->field_id, ctx_ptr);
        field_clear_elem(&ws->inv, ws->field_id, ctx_ptr);
        unified_poly_clear(&ws->tmp);
        unified_poly_clear(&ws->tmp2);
        
        ws->initialized = 0;
        ws->field_id = 0;
        ws->field_ctx = NULL;
    }
}

void cleanup_unified_workspace(void) {
    if (g_unified_workspace.initialized) {
        void *ctx_ptr = g_unified_workspace.field_ctx;
        
        if (ctx_ptr) {
            field_clear_elem(&g_unified_workspace.lc1, g_unified_workspace.field_id, ctx_ptr);
            field_clear_elem(&g_unified_workspace.lc2, g_unified_workspace.field_id, ctx_ptr);
            field_clear_elem(&g_unified_workspace.cst, g_unified_workspace.field_id, ctx_ptr);
            field_clear_elem(&g_unified_workspace.inv, g_unified_workspace.field_id, ctx_ptr);
        }
        
        if (g_unified_workspace.tmp.coeffs) {
            unified_poly_clear(&g_unified_workspace.tmp);
        }
        if (g_unified_workspace.tmp2.coeffs) {
            unified_poly_clear(&g_unified_workspace.tmp2);
        }
        
        g_unified_workspace.initialized = 0;
        g_unified_workspace.field_id = 0;
        g_unified_workspace.field_ctx = NULL;
    }
}

/* Ensure workspace is initialized for the current field */
void ensure_workspace_initialized(field_ctx_t *ctx) {
    void *ctx_ptr = NULL;
    switch (ctx->field_id) {
        case FIELD_ID_NMOD:
            ctx_ptr = (void*)&ctx->ctx.nmod_ctx;
            break;
        case FIELD_ID_FQ_ZECH:
            ctx_ptr = (void*)ctx->ctx.zech_ctx;
            break;
        default:
            ctx_ptr = (void*)ctx->ctx.fq_ctx;
            break;
    }
    
    /* Always reinitialize if field context changes OR if it's FQ/FQ_ZECH field */
    /* For FQ/FQ_ZECH fields, we should be more conservative about reusing workspace */
    if (g_unified_workspace.initialized) {
        int need_reinit = 0;
        
        /* Check if field type changed */
        if (g_unified_workspace.field_id != ctx->field_id) {
            need_reinit = 1;
        }
        /* For FQ/FQ_ZECH fields, check if the field parameters match */
        else if (ctx->field_id == FIELD_ID_FQ || ctx->field_id == FIELD_ID_FQ_ZECH) {
            /* Compare field degree and characteristic */
            if (ctx->field_id == FIELD_ID_FQ) {
                const fq_nmod_ctx_struct *old_ctx = (const fq_nmod_ctx_struct *)g_unified_workspace.field_ctx;
                const fq_nmod_ctx_struct *new_ctx = (const fq_nmod_ctx_struct *)ctx_ptr;
                
                if (!old_ctx || !new_ctx || 
                    fq_nmod_ctx_degree(old_ctx) != fq_nmod_ctx_degree(new_ctx) ||
                    fq_nmod_ctx_prime(old_ctx) != fq_nmod_ctx_prime(new_ctx)) {
                    need_reinit = 1;
                }
            } else if (ctx->field_id == FIELD_ID_FQ_ZECH) {
                const fq_zech_ctx_struct *old_ctx = (const fq_zech_ctx_struct *)g_unified_workspace.field_ctx;
                const fq_zech_ctx_struct *new_ctx = (const fq_zech_ctx_struct *)ctx_ptr;
                
                if (!old_ctx || !new_ctx || old_ctx->qm1 != new_ctx->qm1) {
                    need_reinit = 1;
                }
            }
        }
        /* For other fields, check context pointer */
        else if (g_unified_workspace.field_ctx != ctx_ptr) {
            need_reinit = 1;
        }
        
        if (need_reinit) {
            /* Clear old workspace */
            cleanup_unified_workspace();
        } else {
            /* Workspace already initialized for this field */
            return;
        }
    }
    
    /* Initialize workspace for current field */
    g_unified_workspace.field_id = ctx->field_id;
    g_unified_workspace.field_ctx = ctx_ptr;
    
    field_init_elem(&g_unified_workspace.lc1, ctx->field_id, ctx_ptr);
    field_init_elem(&g_unified_workspace.lc2, ctx->field_id, ctx_ptr);
    field_init_elem(&g_unified_workspace.cst, ctx->field_id, ctx_ptr);
    field_init_elem(&g_unified_workspace.inv, ctx->field_id, ctx_ptr);
    unified_poly_init(&g_unified_workspace.tmp, ctx);
    unified_poly_init(&g_unified_workspace.tmp2, ctx);
    
    g_unified_workspace.initialized = 1;
}

/* Also add cleanup for conversion tables in gf2n_field.h */
void reset_all_gf2n_conversions(void) {
    /* Reset GF(2^8) conversion */
    if (g_gf28_conversion) {
        g_gf28_conversion->initialized = 0;
    }
    
    /* Reset GF(2^16) conversion */
    if (g_gf216_conversion) {
        g_gf216_conversion->initialized = 0;
    }
    
    /* Reset GF(2^32) conversion */
    if (g_gf232_conversion) {
        g_gf232_conversion->initialized = 0;
    }
    
    /* Reset GF(2^64) conversion */
    if (g_gf264_conversion) {
        g_gf264_conversion->initialized = 0;
    }
    
    /* Reset GF(2^128) conversion */
    if (g_gf2128_conversion) {
        g_gf2128_conversion->initialized = 0;
    }
}

/* ============================================================================
   FIELD INITIALIZATION AND CLEANUP
   ============================================================================ */

field_ctx_t* field_init(field_id_t field_id, ulong characteristic) {
    field_ctx_t *ctx = (field_ctx_t*)malloc(sizeof(field_ctx_t));
    if (!ctx) return NULL;
    
    ctx->field_id = field_id;
    ctx->zech_size_limit = 1L << 20;  /* Default limit */
    
    switch (field_id) {
        case FIELD_ID_NMOD:
            /* Prime field */
            nmod_init(&ctx->ctx.nmod_ctx, characteristic);
            ctx->elem_size = sizeof(ulong);
            ctx->description = "Prime field (nmod)";
            break;
            
        case FIELD_ID_GF28:
            {
                /* GF(2^8) */
                fq_nmod_ctx_struct *fq_ctx = (fq_nmod_ctx_struct*)malloc(sizeof(fq_nmod_ctx_struct));
                
                /* Initialize with standard irreducible polynomial x^8 + x^4 + x^3 + x + 1 */
                nmod_poly_t modulus;
                nmod_poly_init(modulus, 2);
                nmod_poly_set_coeff_ui(modulus, 0, 1);
                nmod_poly_set_coeff_ui(modulus, 2, 1);
                nmod_poly_set_coeff_ui(modulus, 3, 1);
                nmod_poly_set_coeff_ui(modulus, 4, 1);
                nmod_poly_set_coeff_ui(modulus, 8, 1);
                
                fq_nmod_ctx_init_modulus(fq_ctx, modulus, "x");
                nmod_poly_clear(modulus);
                
                ctx->ctx.fq_ctx = fq_ctx;
                ctx->elem_size = sizeof(uint8_t);
                ctx->description = "GF(2^8) lookup tables";
                
                /* Initialize GF(2^8) tables and conversion */
                init_gf28_standard();
                init_gf28_conversion(fq_ctx);
            }
            break;
            
        case FIELD_ID_GF216:
            {
                /* GF(2^16) */
                fq_nmod_ctx_struct *fq_ctx = (fq_nmod_ctx_struct*)malloc(sizeof(fq_nmod_ctx_struct));
                
                /* Initialize with standard irreducible polynomial */
                nmod_poly_t modulus;
                nmod_poly_init(modulus, 2);
                nmod_poly_set_coeff_ui(modulus, 0, 1);
                nmod_poly_set_coeff_ui(modulus, 5, 1);
                nmod_poly_set_coeff_ui(modulus, 3, 1);
                nmod_poly_set_coeff_ui(modulus, 1, 1);
                nmod_poly_set_coeff_ui(modulus, 16, 1);
                
                fq_nmod_ctx_init_modulus(fq_ctx, modulus, "x");
                nmod_poly_clear(modulus);
                
                ctx->ctx.fq_ctx = fq_ctx;
                ctx->elem_size = sizeof(uint16_t);
                ctx->description = "GF(2^16) tower field";
                
                /* Initialize GF(2^16) conversion */
                init_gf216_standard();
                init_gf216_conversion(fq_ctx);
            }
            break;
            
        case FIELD_ID_GF232:
            {
                /* GF(2^32) */
                fq_nmod_ctx_struct *fq_ctx = (fq_nmod_ctx_struct*)malloc(sizeof(fq_nmod_ctx_struct));
                
                /* Initialize with standard irreducible polynomial */
                nmod_poly_t modulus;
                nmod_poly_init(modulus, 2);
                nmod_poly_set_coeff_ui(modulus, 0, 1);
                nmod_poly_set_coeff_ui(modulus, 1, 1);
                nmod_poly_set_coeff_ui(modulus, 2, 1);
                nmod_poly_set_coeff_ui(modulus, 22, 1);
                nmod_poly_set_coeff_ui(modulus, 32, 1);
                
                fq_nmod_ctx_init_modulus(fq_ctx, modulus, "x");
                nmod_poly_clear(modulus);
                
                ctx->ctx.fq_ctx = fq_ctx;
                ctx->elem_size = sizeof(gf232_t);
                ctx->description = "GF(2^32) PCLMUL";
                
                /* Initialize GF(2^32) */
                init_gf232();
                init_gf232_conversion(fq_ctx);
            }
            break;
            
        case FIELD_ID_GF264:
            {
                /* GF(2^64) */
                fq_nmod_ctx_struct *fq_ctx = (fq_nmod_ctx_struct*)malloc(sizeof(fq_nmod_ctx_struct));
                
                /* Initialize with standard irreducible polynomial */
                nmod_poly_t modulus;
                nmod_poly_init(modulus, 2);
                nmod_poly_set_coeff_ui(modulus, 0, 1);
                nmod_poly_set_coeff_ui(modulus, 1, 1);
                nmod_poly_set_coeff_ui(modulus, 3, 1);
                nmod_poly_set_coeff_ui(modulus, 4, 1);
                nmod_poly_set_coeff_ui(modulus, 64, 1);
                
                fq_nmod_ctx_init_modulus(fq_ctx, modulus, "x");
                nmod_poly_clear(modulus);
                
                ctx->ctx.fq_ctx = fq_ctx;
                ctx->elem_size = sizeof(gf264_t);
                ctx->description = "GF(2^64) PCLMUL";
                
                /* Initialize GF(2^64) */
                init_gf264();
                init_gf264_conversion(fq_ctx);
            }
            break;
            
        case FIELD_ID_GF2128:
            {
                /* GF(2^128) */
                fq_nmod_ctx_struct *fq_ctx = (fq_nmod_ctx_struct*)malloc(sizeof(fq_nmod_ctx_struct));
                
                /* Initialize with standard irreducible polynomial x^128 + x^7 + x^2 + x + 1 */
                nmod_poly_t modulus;
                nmod_poly_init(modulus, 2);
                nmod_poly_set_coeff_ui(modulus, 0, 1);
                nmod_poly_set_coeff_ui(modulus, 1, 1);
                nmod_poly_set_coeff_ui(modulus, 2, 1);
                nmod_poly_set_coeff_ui(modulus, 7, 1);
                nmod_poly_set_coeff_ui(modulus, 128, 1);
                
                fq_nmod_ctx_init_modulus(fq_ctx, modulus, "x");
                nmod_poly_clear(modulus);
                
                ctx->ctx.fq_ctx = fq_ctx;
                ctx->elem_size = sizeof(gf2128_t);
                ctx->description = "GF(2^128) PCLMUL";
                
                /* Initialize GF(2^128) */
                init_gf2128();
                init_gf2128_conversion(fq_ctx);
            }
            break;
            
        default:
            /* General finite field */
            free(ctx);
            return NULL;
    }
    
    return ctx;
}

void field_clear(field_ctx_t *ctx) {
    if (!ctx) return;
    
    switch (ctx->field_id) {
        case FIELD_ID_NMOD:
            /* Nothing to clear for nmod */
            break;
            
        case FIELD_ID_GF28:
        case FIELD_ID_GF216:
        case FIELD_ID_GF232:
        case FIELD_ID_GF264:
        case FIELD_ID_GF2128:
        case FIELD_ID_FQ:
            if (ctx->ctx.fq_ctx) {
                fq_nmod_ctx_clear((fq_nmod_ctx_struct *)ctx->ctx.fq_ctx);
                free((void*)ctx->ctx.fq_ctx);
            }
            break;
            
        case FIELD_ID_FQ_ZECH:
            if (ctx->ctx.zech_ctx) {
                fq_zech_ctx_clear(ctx->ctx.zech_ctx);
                free(ctx->ctx.zech_ctx);
            }
            break;
    }
    
    free(ctx);
}
