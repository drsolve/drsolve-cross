/* fq_unified_interface.h - Unified field operations header with Zech logarithm support */
#ifndef FQ_UNIFIED_INTERFACE_H
#define FQ_UNIFIED_INTERFACE_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <flint/flint.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>
#include <flint/fq_zech.h>  /* Add Zech support */
#include <flint/nmod.h>
#include "gf2n_field.h"
#include "gf2n_poly.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
   UNIFIED FIELD ELEMENT TYPE WITH ZECH SUPPORT
   ============================================================================ */

/* Extended field element type enumeration */
typedef enum {
    FIELD_ID_NMOD   = 0,  /* Prime field Z/pZ */
    FIELD_ID_GF28   = 1,  /* GF(2^8) */
    FIELD_ID_GF216  = 2,  /* GF(2^16) */
    FIELD_ID_GF232  = 3,  /* GF(2^32) */
    FIELD_ID_GF264  = 4,  /* GF(2^64) */
    FIELD_ID_GF2128 = 5,  /* GF(2^128) */
    FIELD_ID_FQ_ZECH = 6, /* Small finite fields with Zech logarithm */
    FIELD_ID_FQ     = 7   /* General finite field */
} field_id_t;

/* Unified field element using union with Zech support */
typedef union {
    ulong nmod;              /* For prime fields */
    uint8_t gf28;            /* For GF(2^8) */
    uint16_t gf216;          /* For GF(2^16) */
    gf232_t gf232;           /* For GF(2^32) */
    gf264_t gf264;           /* For GF(2^64) */
    gf2128_t gf2128;         /* For GF(2^128) */
    fq_zech_struct fq_zech;  /* For Zech logarithm fields */
    fq_nmod_struct fq;       /* For general fields */
} field_elem_u;

/* ============================================================================
   FIELD CONTEXT STRUCTURE WITH ZECH SUPPORT
   ============================================================================ */

typedef struct {
    field_id_t field_id;
    size_t elem_size;
    ulong zech_size_limit;   /* Size limit for using Zech logarithms */
    
    /* Context data */
    union {
        nmod_t nmod_ctx;                      /* For prime fields */
        const fq_nmod_ctx_struct *fq_ctx;     /* For general fields */
        fq_zech_ctx_struct *zech_ctx;         /* For Zech logarithm fields - owned */
    } ctx;
    
    /* Optimization info */
    const char *description;
} field_ctx_t;

/* ============================================================================
   UNIFIED POLYNOMIAL TYPE
   ============================================================================ */

typedef struct {
    field_elem_u *coeffs;     /* Array of coefficients */
    slong length;             /* Current length */
    slong alloc;              /* Allocated size */
    field_ctx_t *ctx;         /* Field context */
} unified_poly_struct;

typedef unified_poly_struct unified_poly_t[1];

/* ============================================================================
   POLYNOMIAL MATRIX TYPE
   ============================================================================ */

typedef struct {
    unified_poly_struct *entries;
    slong r, c;
    unified_poly_struct **rows;
    field_ctx_t *ctx;
} unified_poly_mat_struct;

typedef unified_poly_mat_struct unified_poly_mat_t[1];

/* ============================================================================
   ZECH LOGARITHM DETECTION AND UTILITIES
   ============================================================================ */

#ifndef CALCULATE_FIELD_SIZE_DEFINED
#define CALCULATE_FIELD_SIZE_DEFINED

/* Calculate field size p^d with overflow protection */
static inline ulong calculate_field_size(ulong p, slong d) {
    if (d <= 0) return 0;
    if (d > 30) return UWORD_MAX;  /* Avoid overflow */
    
    ulong result = 1;
    for (slong i = 0; i < d; i++) {
        if (result > UWORD_MAX / p) return UWORD_MAX;  /* Overflow */
        result *= p;
    }
    return result;
}

/* Check if field should use Zech logarithm representation */
static inline int should_use_zech_logarithm(const fq_nmod_ctx_t fq_ctx, ulong size_limit) {
    ulong p = fq_nmod_ctx_prime(fq_ctx);
    slong d = fq_nmod_ctx_degree(fq_ctx);
    
    /* Don't use Zech for prime fields */
    if (d == 1) return 0;
    
    /* Don't use Zech for already optimized GF(2^n) fields */
    if (p == 2) {
        switch (d) {
            case 8:   /* GF(2^8) has custom optimization */
            case 16:  /* GF(2^16) has custom optimization */
            case 32:  /* GF(2^32) has custom optimization */
            case 64:  /* GF(2^64) has custom optimization */
            case 128: /* GF(2^128) has custom optimization */
                return 0;
            default:
                break;
        }
    }
    
    /* Check field size */
    ulong field_size = calculate_field_size(p, d);
    if (field_size == UWORD_MAX || field_size > size_limit) {
        return 0;
    }
    
    /* Additional checks could be added here:
     * - Verify polynomial is primitive
     * - Check memory constraints
     * - Consider performance characteristics
     */
    return 1;
}

#endif /* CALCULATE_FIELD_SIZE_DEFINED */

/* ============================================================================
   WORKSPACE FOR HOT PATHS
   ============================================================================ */

typedef struct {
    field_elem_u lc1, lc2, cst, inv;
    unified_poly_struct tmp, tmp2;
    field_id_t field_id;  /* Track which field this workspace is for */
    void *field_ctx;      /* Store the context pointer */
    int initialized;
} unified_workspace_t;

/* Per-thread workspace */
extern __thread unified_workspace_t g_unified_workspace;

/* ============================================================================
   FUNCTION DECLARATIONS - CONTEXT MANAGEMENT
   ============================================================================ */

/* Try to create Zech context from fq_nmod context */
int try_create_zech_context(fq_zech_ctx_struct *zech_ctx, const fq_nmod_ctx_t fq_ctx);

/* Enhanced field context initialization with Zech support */
void field_ctx_init_enhanced(field_ctx_t *ctx, const fq_nmod_ctx_t fq_ctx, ulong zech_limit);

/* Standard field context initialization (backward compatibility) */
void field_ctx_init(field_ctx_t *ctx, const fq_nmod_ctx_t fq_ctx);

/* Context management functions */
void field_ctx_set_zech_limit(field_ctx_t *ctx, ulong limit);
const char* field_ctx_get_description(const field_ctx_t *ctx);
int field_ctx_is_using_zech(const field_ctx_t *ctx);
ulong field_ctx_get_size(const field_ctx_t *ctx);
void field_ctx_clear(field_ctx_t *ctx);

/* ============================================================================
   FUNCTION DECLARATIONS - FIELD OPERATIONS
   ============================================================================ */

/* Field operations - non-inline versions */
void field_neg(field_elem_u *res, const field_elem_u *a, field_id_t field_id, const void *ctx);
void field_inv(field_elem_u *res, const field_elem_u *a, field_id_t field_id, const void *ctx);
void field_set_zero(field_elem_u *res, field_id_t field_id, const void *ctx);
void field_set_one(field_elem_u *res, field_id_t field_id, const void *ctx);
int field_equal(const field_elem_u *a, const field_elem_u *b, field_id_t field_id, const void *ctx);
void field_init_elem(field_elem_u *elem, field_id_t field_id, const void *ctx);
void field_clear_elem(field_elem_u *elem, field_id_t field_id, const void *ctx);
void field_set_elem(field_elem_u *res, const field_elem_u *a, field_id_t field_id, const void *ctx);
void field_sub(field_elem_u *res, const field_elem_u *a, const field_elem_u *b,
               field_id_t field_id, const void *ctx);

/* Conversion functions */
void fq_nmod_to_field_elem(field_elem_u *res, const fq_nmod_t elem, const field_ctx_t *ctx);
void field_elem_to_fq_nmod(fq_nmod_t res, const field_elem_u *elem, const field_ctx_t *ctx);

/* ============================================================================
   FUNCTION DECLARATIONS - POLYNOMIAL OPERATIONS
   ============================================================================ */

/* Polynomial operations */
void unified_poly_init(unified_poly_t poly, field_ctx_t *ctx);
void unified_poly_clear(unified_poly_t poly);
void unified_poly_fit_length(unified_poly_t poly, slong len);
void unified_poly_normalise(unified_poly_t poly);
void unified_poly_zero(unified_poly_t poly);
int unified_poly_is_zero(const unified_poly_t poly);
slong unified_poly_degree(const unified_poly_t poly);
void unified_poly_set(unified_poly_t res, const unified_poly_t poly);
void unified_poly_add(unified_poly_t res, const unified_poly_t a, const unified_poly_t b);
void unified_poly_scalar_mul(unified_poly_t res, const unified_poly_t poly, const field_elem_u *c);
void unified_poly_mul(unified_poly_t res, const unified_poly_t a, const unified_poly_t b);
void unified_poly_get_coeff(field_elem_u *coeff, const unified_poly_t poly, slong i);
void unified_poly_set_coeff(unified_poly_t poly, slong i, const field_elem_u *coeff);
void unified_poly_shift_left(unified_poly_t res, const unified_poly_t poly, slong n);

/* Conversion functions */
void fq_nmod_poly_to_unified(unified_poly_t res, const fq_nmod_poly_t poly,
                            const fq_nmod_ctx_t ctx, field_ctx_t *field_ctx);
void unified_to_fq_nmod_poly(fq_nmod_poly_t res, const unified_poly_t poly,
                            const fq_nmod_ctx_t ctx, field_ctx_t *field_ctx);

/* ============================================================================
   FUNCTION DECLARATIONS - MATRIX OPERATIONS
   ============================================================================ */

/* Matrix operations */
void unified_poly_mat_init(unified_poly_mat_t mat, slong rows, slong cols, field_ctx_t *ctx);
void unified_poly_mat_clear(unified_poly_mat_t mat);
unified_poly_struct *unified_poly_mat_entry(unified_poly_mat_t mat, slong i, slong j);
const unified_poly_struct *unified_poly_mat_entry_const(const unified_poly_mat_t mat, slong i, slong j);
void fq_nmod_poly_mat_to_unified(unified_poly_mat_t res, const fq_nmod_poly_mat_t mat,
                                const fq_nmod_ctx_t ctx, field_ctx_t *field_ctx);
void unified_to_fq_nmod_poly_mat(fq_nmod_poly_mat_t res, const unified_poly_mat_t mat,
                                const fq_nmod_ctx_t ctx, field_ctx_t *field_ctx);

/* ============================================================================
   FUNCTION DECLARATIONS - WORKSPACE MANAGEMENT
   ============================================================================ */

/* Workspace management */
void ensure_workspace_initialized(field_ctx_t *ctx);
void cleanup_unified_workspace(void);
void reset_all_gf2n_conversions(void);

/* Field creation and cleanup */
field_ctx_t* field_init(field_id_t field_id, ulong characteristic);
void field_clear(field_ctx_t *ctx);

/* ============================================================================
   INLINE IMPLEMENTATIONS - ONLY THE MOST CRITICAL HOT PATH FUNCTIONS
   ============================================================================ */

/* Core multiplication - inline for hot path */
static inline void field_mul(field_elem_u *res, const field_elem_u *a, const field_elem_u *b, 
                            field_id_t field_id, const void *ctx) {
    switch (field_id) {
        case FIELD_ID_GF28:
            res->gf28 = gf28_mul(a->gf28, b->gf28);
            break;
        case FIELD_ID_GF216:
            res->gf216 = gf216_mul(a->gf216, b->gf216);
            break;
        case FIELD_ID_GF232:
            res->gf232 = gf232_mul(&a->gf232, &b->gf232);
            break;
        case FIELD_ID_GF264:
            res->gf264 = gf264_mul(&a->gf264, &b->gf264);
            break;
        case FIELD_ID_GF2128:
            res->gf2128 = gf2128_mul(&a->gf2128, &b->gf2128);
            break;
        case FIELD_ID_NMOD:
            res->nmod = nmod_mul(a->nmod, b->nmod, *(const nmod_t*)ctx);
            break;
        case FIELD_ID_FQ_ZECH:
            /* Zech logarithm multiplication - very fast! */
            if (res != a && res != b) {
                fq_zech_init(&res->fq_zech, (const fq_zech_ctx_struct *)ctx);
            }
            fq_zech_mul(&res->fq_zech, &a->fq_zech, &b->fq_zech, (const fq_zech_ctx_struct *)ctx);
            break;
        case FIELD_ID_FQ:
            /* Ensure result is initialized before operation */
            if (res != a && res != b) {
                fq_nmod_init(&res->fq, (const fq_nmod_ctx_struct *)ctx);
            }
            fq_nmod_mul(&res->fq, &a->fq, &b->fq, (const fq_nmod_ctx_struct *)ctx);
            break;
    }
}

/* Core addition - inline for hot path */
static inline void field_add(field_elem_u *res, const field_elem_u *a, const field_elem_u *b,
                            field_id_t field_id, const void *ctx) {
    switch (field_id) {
        case FIELD_ID_GF28:
            res->gf28 = a->gf28 ^ b->gf28;
            break;
        case FIELD_ID_GF216:
            res->gf216 = a->gf216 ^ b->gf216;
            break;
        case FIELD_ID_GF232:
            res->gf232.value = a->gf232.value ^ b->gf232.value;
            break;
        case FIELD_ID_GF264:
            res->gf264.value = a->gf264.value ^ b->gf264.value;
            break;
        case FIELD_ID_GF2128:
            res->gf2128.low = a->gf2128.low ^ b->gf2128.low;
            res->gf2128.high = a->gf2128.high ^ b->gf2128.high;
            break;
        case FIELD_ID_NMOD:
            res->nmod = nmod_add(a->nmod, b->nmod, *(const nmod_t*)ctx);
            break;
        case FIELD_ID_FQ_ZECH:
            /* Zech logarithm addition */
            if (res != a && res != b) {
                fq_zech_init(&res->fq_zech, (const fq_zech_ctx_struct *)ctx);
            }
            fq_zech_add(&res->fq_zech, &a->fq_zech, &b->fq_zech, (const fq_zech_ctx_struct *)ctx);
            break;
        case FIELD_ID_FQ:
            /* Ensure result is initialized before operation */
            if (res != a && res != b) {
                fq_nmod_init(&res->fq, (const fq_nmod_ctx_struct *)ctx);
            }
            fq_nmod_add(&res->fq, &a->fq, &b->fq, (const fq_nmod_ctx_struct *)ctx);
            break;
    }
}

/* Check if zero - inline for hot path */
static inline int field_is_zero(const field_elem_u *a, field_id_t field_id, const void *ctx) {
    switch (field_id) {
        case FIELD_ID_GF28:
            return a->gf28 == 0;
        case FIELD_ID_GF216:
            return a->gf216 == 0;
        case FIELD_ID_GF232:
            return a->gf232.value == 0;
        case FIELD_ID_GF264:
            return a->gf264.value == 0;
        case FIELD_ID_GF2128:
            return a->gf2128.low == 0 && a->gf2128.high == 0;
        case FIELD_ID_NMOD:
            return a->nmod == 0;
        case FIELD_ID_FQ_ZECH:
            return fq_zech_is_zero(&a->fq_zech, (const fq_zech_ctx_struct *)ctx);
        case FIELD_ID_FQ:
            return fq_nmod_is_zero(&a->fq, (const fq_nmod_ctx_struct *)ctx);
    }
    return 0;
}

/* Check if one - inline for hot path */
static inline int field_is_one(const field_elem_u *a, field_id_t field_id, const void *ctx) {
    switch (field_id) {
        case FIELD_ID_GF28:
            return a->gf28 == 1;
        case FIELD_ID_GF216:
            return a->gf216 == 1;
        case FIELD_ID_GF232:
            return a->gf232.value == 1;
        case FIELD_ID_GF264:
            return a->gf264.value == 1;
        case FIELD_ID_GF2128:
            return a->gf2128.low == 1 && a->gf2128.high == 0;
        case FIELD_ID_NMOD:
            return a->nmod == 1;
        case FIELD_ID_FQ_ZECH:
            return fq_zech_is_one(&a->fq_zech, (const fq_zech_ctx_struct *)ctx);
        case FIELD_ID_FQ:
            return fq_nmod_is_one(&a->fq, (const fq_nmod_ctx_struct *)ctx);
    }
    return 0;
}

/* Fast degree computation - inline for hot path */
static inline slong unified_poly_degree_fast(const unified_poly_struct *poly) {
    slong len = poly->length;
    if (len == 0) return -1;
    
    /* Check for trailing zeros for ALL field types, not just GF(2^n) */
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
    
    /* Remove trailing zeros */
    while (len > 0 && field_is_zero(&poly->coeffs[len-1], poly->ctx->field_id, ctx_ptr)) {
        len--;
    }
    
    return len - 1;
}

/* Get coefficient pointer directly - inline for hot path */
static inline field_elem_u* unified_poly_get_coeff_ptr(unified_poly_struct *poly, slong i) {
    if (i < poly->length) {
        return &poly->coeffs[i];
    }
    return NULL;
}

/* In-place polynomial addition - inline for hot path */
static inline void unified_poly_add_inplace(unified_poly_struct *res, 
                                           const unified_poly_struct *a,
                                           field_ctx_t *ctx) {
    slong min_len = FLINT_MIN(res->length, a->length);
    
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
    
    /* Add coefficients */
    for (slong i = 0; i < min_len; i++) {
        field_add(&res->coeffs[i], &res->coeffs[i], &a->coeffs[i], ctx->field_id, ctx_ptr);
    }
    
    /* Handle remaining coefficients if a is longer */
    if (a->length > res->length) {
        unified_poly_fit_length(res, a->length);
        for (slong i = res->length; i < a->length; i++) {
            field_set_elem(&res->coeffs[i], &a->coeffs[i], ctx->field_id, ctx_ptr);
        }
        res->length = a->length;
    }
}

/* Combined shift-left and add operation - inline for hot path */
static inline void unified_poly_shift_left_add_inplace(unified_poly_struct *res,
                                                      const unified_poly_struct *a,
                                                      slong shift,
                                                      field_ctx_t *ctx) {
    if (shift == 0) {
        unified_poly_add_inplace(res, a, ctx);
        return;
    }
    
    slong new_len = a->length + shift;
    unified_poly_fit_length(res, new_len);
    
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
    
    /* Ensure we have zeros up to shift position */
    if (res->length < shift) {
        for (slong i = res->length; i < shift; i++) {
            field_set_zero(&res->coeffs[i], ctx->field_id, ctx_ptr);
        }
    }
    
    /* Add shifted coefficients */
    for (slong i = 0; i < a->length; i++) {
        if (i + shift < res->length) {
            field_add(&res->coeffs[i + shift], &res->coeffs[i + shift], 
                     &a->coeffs[i], ctx->field_id, ctx_ptr);
        } else {
            field_set_elem(&res->coeffs[i + shift], &a->coeffs[i], ctx->field_id, ctx_ptr);
        }
    }
    
    if (new_len > res->length) {
        res->length = new_len;
    }
}

/* ============================================================================
   FLINT VERSION COMPATIBILITY
   ============================================================================ */

#if __FLINT_VERSION >= 3
#define FLINT_RAND_INIT(state) flint_rand_init(state)
#define FLINT_RAND_CLEAR(state) flint_rand_clear(state)
#define FLINT_RAND_SEED(state, seed) flint_rand_set_seed(state, seed, seed + 1)
#else
#define FLINT_RAND_INIT(state) flint_randinit(state)
#define FLINT_RAND_CLEAR(state) flint_randclear(state)
#define FLINT_RAND_SEED(state, seed) flint_randseed(state, seed, seed)
#endif

#ifdef __cplusplus
}
#endif

#endif /* FQ_UNIFIED_INTERFACE_H */