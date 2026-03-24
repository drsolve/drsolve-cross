/* unified_mpoly_interface.h - Complete single-file implementation with Zech logarithm support */

/*
 * NOTE: This implementation includes support for fq_zech_mpoly (Zech logarithm representation).
 * However, some fq_zech_mpoly functions may not be available in all FLINT versions:
 * - Direct coefficient access functions are not available, so we use workarounds
 * - evaluate_all function is implemented using sequential single-variable evaluation
 * - compose function is not implemented for fq_zech_mpoly
 * 
 * For maximum compatibility, consider using fq_nmod_mpoly instead of fq_zech_mpoly
 * unless the Zech logarithm optimization is specifically needed.
 */

#ifndef UNIFIED_MPOLY_INTERFACE_H
#define UNIFIED_MPOLY_INTERFACE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <flint/flint.h>
#include <flint/nmod_mpoly.h>
#include <flint/fq_nmod_mpoly.h>
#include <flint/fq_zech_mpoly.h>
#include <flint/fmpz.h>
#include <flint/mpoly.h>
#include <flint/longlong.h>
#include "fq_unified_interface.h"
#include "gf2n_mpoly.h"


#define GET_NMOD_POLY(poly) (&(poly)->data.nmod_poly)
#define GET_FQ_POLY(poly) (&(poly)->data.fq_poly)
#define GET_ZECH_POLY(poly) (&(poly)->data.zech_poly)

#define GET_NMOD_CTX(ctx) (&(ctx)->ctx.nmod_ctx)
#define GET_FQ_CTX(ctx) (&(ctx)->ctx.fq_ctx)
#define GET_ZECH_CTX(ctx) (&(ctx)->ctx.zech_ctx)


#define WARN_ONCE_N(n, fmt, ...) \
    do { \
        static int _warn_count = 0; \
        if (_warn_count < (n)) { \
            _warn_count++; \
            fprintf(stderr, fmt, ##__VA_ARGS__); \
            if (_warn_count == (n)) \
                fprintf(stderr, "(further identical warnings suppressed)\n"); \
        } \
    } while (0)

#define WARN_ONCE(fmt, ...)    WARN_ONCE_N(1, fmt, ##__VA_ARGS__)
#define WARN_THRICE(fmt, ...)  WARN_ONCE_N(3, fmt, ##__VA_ARGS__)

/* ============================================================================
   UNIFIED MULTIVARIATE POLYNOMIAL TYPES AND STRUCTURES WITH ZECH SUPPORT
   ============================================================================ */

/* Context for multivariate polynomials with Zech support */
typedef struct {
    field_ctx_t *field_ctx;
    slong nvars;
    ordering_t ord;
    
    /* Different contexts based on field type */
    union {
        nmod_mpoly_ctx_struct nmod_ctx;
        fq_nmod_mpoly_ctx_struct fq_ctx;
        fq_zech_mpoly_ctx_struct zech_ctx;
    } ctx;
    
    /* Method selection parameters */
    slong array_size_limit;
    slong dense_threshold;
    ulong zech_size_limit;  /* Zech logarithm table size limit (default 2^20) */
} unified_mpoly_ctx_struct;

typedef unified_mpoly_ctx_struct *unified_mpoly_ctx_t;

/* Unified multivariate polynomial type with Zech support */
typedef struct {
    /* Field type identifier */
    field_id_t field_id;
    
    /* Storage for different polynomial types */
    union {
        nmod_mpoly_struct nmod_poly;
        fq_nmod_mpoly_struct fq_poly;
        fq_zech_mpoly_struct zech_poly;
    } data;
    
    /* Context pointer */
    unified_mpoly_ctx_t ctx_ptr;
} unified_mpoly_struct;

typedef unified_mpoly_struct *unified_mpoly_t;

/* ============================================================================
   HELPER MACROS WITH ZECH SUPPORT
   ============================================================================ */

#define GET_NMOD_POLY(poly) (&(poly)->data.nmod_poly)
#define GET_FQ_POLY(poly) (&(poly)->data.fq_poly)
#define GET_ZECH_POLY(poly) (&(poly)->data.zech_poly)

#define GET_NMOD_CTX(ctx) (&(ctx)->ctx.nmod_ctx)
#define GET_FQ_CTX(ctx) (&(ctx)->ctx.fq_ctx)
#define GET_ZECH_CTX(ctx) (&(ctx)->ctx.zech_ctx)

/* ============================================================================
   FUNCTION DECLARATIONS
   ============================================================================ */

/* Context operations */
unified_mpoly_ctx_t unified_mpoly_ctx_init(slong nvars, const ordering_t ord, field_ctx_t *field_ctx);
void unified_mpoly_ctx_clear(unified_mpoly_ctx_t ctx);

/* Polynomial memory management */
unified_mpoly_t unified_mpoly_init(unified_mpoly_ctx_t ctx);
void unified_mpoly_clear(unified_mpoly_t poly);

/* Basic operations */
void unified_mpoly_zero(unified_mpoly_t poly);
void unified_mpoly_one(unified_mpoly_t poly);
int unified_mpoly_is_zero(const unified_mpoly_t poly);
int unified_mpoly_is_one(const unified_mpoly_t poly);
slong unified_mpoly_length(const unified_mpoly_t poly);
void unified_mpoly_set(unified_mpoly_t poly1, const unified_mpoly_t poly2);

/* Arithmetic operations */
void unified_mpoly_add(unified_mpoly_t poly1, const unified_mpoly_t poly2, const unified_mpoly_t poly3);
void unified_mpoly_sub(unified_mpoly_t poly1, const unified_mpoly_t poly2, const unified_mpoly_t poly3);
void unified_mpoly_neg(unified_mpoly_t poly1, const unified_mpoly_t poly2);

/* Multiplication with optimization support */
void unified_mpoly_enable_optimizations(field_id_t field_id, int enable);
int unified_mpoly_mul(unified_mpoly_t poly1, const unified_mpoly_t poly2, const unified_mpoly_t poly3);

/* Division operations */
void unified_mpoly_enable_div_optimizations(field_id_t field_id, int enable);
int unified_mpoly_divides(unified_mpoly_t Q, const unified_mpoly_t A, const unified_mpoly_t B);
void unified_mpoly_divrem(unified_mpoly_t Q, unified_mpoly_t R, const unified_mpoly_t A, const unified_mpoly_t B);
int unified_mpoly_divexact(unified_mpoly_t Q, const unified_mpoly_t A, const unified_mpoly_t B);

/* Coefficient access */
void unified_mpoly_set_coeff_ui(unified_mpoly_t poly, const field_elem_u *c, const ulong *exp);
void unified_mpoly_get_coeff_ui(field_elem_u *c, const unified_mpoly_t poly, const ulong *exp);

/* Degree operations */
void unified_mpoly_degrees_si(slong *degs, const unified_mpoly_t poly);
slong unified_mpoly_total_degree_si(const unified_mpoly_t poly);

/* Printing */
void unified_mpoly_print_pretty(const unified_mpoly_t poly, const char **vars);

/* Scalar multiplication */
void unified_mpoly_scalar_mul_ui(unified_mpoly_t poly1, const unified_mpoly_t poly2, ulong c);

/* Ring operations */
void unified_mpoly_swap(unified_mpoly_t poly1, unified_mpoly_t poly2);
void unified_mpoly_set_fmpz(unified_mpoly_t poly, const fmpz_t c);
void unified_mpoly_mul_fmpz(unified_mpoly_t poly1, const unified_mpoly_t poly2, const fmpz_t c);
int unified_mpoly_pow_fmpz(unified_mpoly_t poly1, const unified_mpoly_t poly2, const fmpz_t exp);
slong unified_mpoly_length_wrapper(const void *a, const void *ctx);

/* Generator functions */
void unified_mpoly_gen(unified_mpoly_t poly, slong var, unified_mpoly_ctx_t ctx);

/* Ring structure */
typedef struct {
    mpoly_void_ring_t ring;
    unified_mpoly_ctx_t ctx;
} unified_mpoly_ring_struct;

typedef unified_mpoly_ring_struct *unified_mpoly_ring_t;

void unified_mpoly_ring_init(mpoly_void_ring_t R, unified_mpoly_ctx_t ctx);
void unified_mpoly_ring_clear(unified_mpoly_ring_t R);

/* Convenience functions */
void unified_mpoly_enable_all_optimizations(field_id_t field_id);
void unified_mpoly_disable_all_optimizations(field_id_t field_id);
void unified_mpoly_set_zech_limit(unified_mpoly_ctx_t ctx, ulong limit);
const char* unified_mpoly_get_field_info(const unified_mpoly_t poly);
const char* unified_mpoly_get_field_info_by_ctx(unified_mpoly_ctx_t ctx);

/* Testing utilities */
void unified_mpoly_randtest(unified_mpoly_t poly, flint_rand_t state, slong length, slong exp_bound);
int unified_mpoly_equal(const unified_mpoly_t A, const unified_mpoly_t B);

/* Additional utility functions */
int unified_mpoly_is_univariate(const unified_mpoly_t poly);
slong unified_mpoly_main_variable(const unified_mpoly_t poly);
void unified_mpoly_evaluate_all_ui(field_elem_u *result, const unified_mpoly_t poly, const ulong *values);
int unified_mpoly_compose(unified_mpoly_t result, const unified_mpoly_t poly, unified_mpoly_t * const *args, unified_mpoly_ctx_t ctx);

#endif /* UNIFIED_MPOLY_INTERFACE_H */
