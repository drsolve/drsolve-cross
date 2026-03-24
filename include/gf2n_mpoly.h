/* gf2n_mpoly.h - Unified GF(2^n) Multivariate Polynomial Header */
#ifndef GF2N_MPOLY_H
#define GF2N_MPOLY_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>

#include <flint/flint.h>
#include <flint/mpoly.h>
#include <flint/fmpz.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mpoly.h>
#include <flint/longlong.h>
#include "gf2n_field.h"
#include "gf2n_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Debug flag */
#define DEBUG_DIVISION 0

/* Block size for array operations */
#define BLOCK 256

/* ============================================================================
   MACRO DEFINITIONS FOR GENERIC OPERATIONS
   ============================================================================ */

/* Generic field operations macros */
#define GF_ADD(field, a, b) gf##field##_add(a, b)
#define GF_MUL(field, a, b) gf##field##_mul(a, b)
#define GF_INV(field, a) gf##field##_inv(a)
#define GF_IS_ZERO(field, a) gf##field##_is_zero(a)
#define GF_ZERO(field) gf##field##_zero()

/* ============================================================================
   GF(2^8) MULTIVARIATE POLYNOMIAL STRUCTURES
   ============================================================================ */

typedef struct {
    uint8_t *coeffs;
    ulong *exps;
    slong length;
    slong coeffs_alloc;
    slong exps_alloc;
    flint_bitcnt_t bits;
} gf28_mpoly_struct;

typedef gf28_mpoly_struct gf28_mpoly_t[1];

typedef struct {
    mpoly_ctx_t minfo;
} gf28_mpoly_ctx_struct;

typedef gf28_mpoly_ctx_struct gf28_mpoly_ctx_t[1];

/* Dense polynomial structures for division */
typedef struct {
    uint8_t *coeffs;
    slong alloc;
    slong *deg_bounds;
    slong nvars;
} gf28_mpolyd_struct;

typedef gf28_mpolyd_struct gf28_mpolyd_t[1];

/* Dynamic array structure */
typedef struct {
    uint8_t *data;
    slong size;
} gf28_dynamic_array_t;

/* ============================================================================
   GF(2^16) MULTIVARIATE POLYNOMIAL STRUCTURES
   ============================================================================ */

typedef struct {
    uint16_t *coeffs;
    ulong *exps;
    slong length;
    slong coeffs_alloc;
    slong exps_alloc;
    flint_bitcnt_t bits;
} gf216_mpoly_struct;

typedef gf216_mpoly_struct gf216_mpoly_t[1];

typedef struct {
    mpoly_ctx_t minfo;
} gf216_mpoly_ctx_struct;

typedef gf216_mpoly_ctx_struct gf216_mpoly_ctx_t[1];

typedef struct {
    uint16_t *coeffs;
    slong alloc;
    slong *deg_bounds;
    slong nvars;
} gf216_mpolyd_struct;

typedef gf216_mpolyd_struct gf216_mpolyd_t[1];

typedef struct {
    uint16_t *data;
    slong size;
} gf216_dynamic_array_t;

/* ============================================================================
   GF(2^32) MULTIVARIATE POLYNOMIAL STRUCTURES
   ============================================================================ */

typedef struct {
    gf232_t *coeffs;
    ulong *exps;
    slong length;
    slong coeffs_alloc;
    slong exps_alloc;
    flint_bitcnt_t bits;
} gf232_mpoly_struct;

typedef gf232_mpoly_struct gf232_mpoly_t[1];

typedef struct {
    mpoly_ctx_t minfo;
} gf232_mpoly_ctx_struct;

typedef gf232_mpoly_ctx_struct gf232_mpoly_ctx_t[1];

typedef struct {
    gf232_t *coeffs;
    slong alloc;
    slong *deg_bounds;
    slong nvars;
} gf232_mpolyd_struct;

typedef gf232_mpolyd_struct gf232_mpolyd_t[1];

typedef struct {
    gf232_t *data;
    slong size;
} gf232_dynamic_array_t;

/* ============================================================================
   GF(2^64) MULTIVARIATE POLYNOMIAL STRUCTURES
   ============================================================================ */

typedef struct {
    gf264_t *coeffs;
    ulong *exps;
    slong length;
    slong coeffs_alloc;
    slong exps_alloc;
    flint_bitcnt_t bits;
} gf264_mpoly_struct;

typedef gf264_mpoly_struct gf264_mpoly_t[1];

typedef struct {
    mpoly_ctx_t minfo;
} gf264_mpoly_ctx_struct;

typedef gf264_mpoly_ctx_struct gf264_mpoly_ctx_t[1];

typedef struct {
    gf264_t *coeffs;
    slong alloc;
    slong *deg_bounds;
    slong nvars;
} gf264_mpolyd_struct;

typedef gf264_mpolyd_struct gf264_mpolyd_t[1];

typedef struct {
    gf264_t *data;
    slong size;
} gf264_dynamic_array_t;

/* ============================================================================
   GF(2^128) MULTIVARIATE POLYNOMIAL STRUCTURES
   ============================================================================ */

typedef struct {
    gf2128_t *coeffs;
    ulong *exps;
    slong length;
    slong coeffs_alloc;
    slong exps_alloc;
    flint_bitcnt_t bits;
} gf2128_mpoly_struct;

typedef gf2128_mpoly_struct gf2128_mpoly_t[1];

typedef struct {
    mpoly_ctx_t minfo;
} gf2128_mpoly_ctx_struct;

typedef gf2128_mpoly_ctx_struct gf2128_mpoly_ctx_t[1];

typedef struct {
    gf2128_t *coeffs;
    slong alloc;
    slong *deg_bounds;
    slong nvars;
} gf2128_mpolyd_struct;

typedef gf2128_mpolyd_struct gf2128_mpolyd_t[1];

typedef struct {
    gf2128_t *data;
    slong size;
    slong chunk_size;
    int use_chunks;
} gf2128_dynamic_array_t;

/* ============================================================================
   FUNCTION DECLARATIONS - GF(2^8)
   ============================================================================ */

/* Basic operations */
void gf28_mpoly_init(gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx);
void gf28_mpoly_clear(gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx);
void gf28_mpoly_ctx_init(gf28_mpoly_ctx_t ctx, slong nvars, const ordering_t ord);
void gf28_mpoly_ctx_clear(gf28_mpoly_ctx_t ctx);
void gf28_mpoly_zero(gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx);
void gf28_mpoly_set(gf28_mpoly_t res, const gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx);
/* Coefficient access */
void gf28_mpoly_set_coeff_ui_ui(gf28_mpoly_t poly, uint8_t c, 
                                const ulong *exp, const gf28_mpoly_ctx_t ctx);

/* Multiplication */
int gf28_mpoly_mul(gf28_mpoly_t res, const gf28_mpoly_t a, const gf28_mpoly_t b, 
                   const gf28_mpoly_ctx_t ctx);
int gf28_mpoly_mul_array(gf28_mpoly_t A, const gf28_mpoly_t B,
                         const gf28_mpoly_t C, const gf28_mpoly_ctx_t ctx);

/* Division */
int gf28_mpoly_divides(gf28_mpoly_t Q, const gf28_mpoly_t A, 
                       const gf28_mpoly_t B, const gf28_mpoly_ctx_t ctx);

/* Utility */
void gf28_mpoly_print(const gf28_mpoly_t poly, const char **vars, 
                      const gf28_mpoly_ctx_t ctx);
int gf28_mpoly_equal(const gf28_mpoly_t A, const gf28_mpoly_t B, 
                     const gf28_mpoly_ctx_t ctx);
void gf28_mpoly_randtest(gf28_mpoly_t poly, flint_rand_t state,
                         slong length, slong exp_bound, 
                         const gf28_mpoly_ctx_t ctx);

/* Conversion */
void fq_nmod_mpoly_to_gf28_mpoly(gf28_mpoly_t res, const fq_nmod_mpoly_t poly,
                                 const fq_nmod_ctx_t fqctx, 
                                 const fq_nmod_mpoly_ctx_t fq_mpoly_ctx);
void gf28_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf28_mpoly_t poly,
                                 const fq_nmod_ctx_t fqctx, 
                                 const fq_nmod_mpoly_ctx_t fq_mpoly_ctx);

/* ============================================================================
   FUNCTION DECLARATIONS - GF(2^16)
   ============================================================================ */

void gf216_mpoly_init(gf216_mpoly_t poly, const gf216_mpoly_ctx_t ctx);
void gf216_mpoly_clear(gf216_mpoly_t poly, const gf216_mpoly_ctx_t ctx);
void gf216_mpoly_ctx_init(gf216_mpoly_ctx_t ctx, slong nvars, const ordering_t ord);
void gf216_mpoly_ctx_clear(gf216_mpoly_ctx_t ctx);
void gf216_mpoly_zero(gf216_mpoly_t poly, const gf216_mpoly_ctx_t ctx);
void gf216_mpoly_set(gf216_mpoly_t res, const gf216_mpoly_t poly, const gf216_mpoly_ctx_t ctx);

void gf216_mpoly_set_coeff_ui_ui(gf216_mpoly_t poly, uint16_t c, 
                                 const ulong *exp, const gf216_mpoly_ctx_t ctx);

int gf216_mpoly_mul(gf216_mpoly_t res, const gf216_mpoly_t a, const gf216_mpoly_t b, 
                    const gf216_mpoly_ctx_t ctx);
int gf216_mpoly_mul_array(gf216_mpoly_t A, const gf216_mpoly_t B,
                          const gf216_mpoly_t C, const gf216_mpoly_ctx_t ctx);

int gf216_mpoly_divides(gf216_mpoly_t Q, const gf216_mpoly_t A, 
                        const gf216_mpoly_t B, const gf216_mpoly_ctx_t ctx);

void gf216_mpoly_print(const gf216_mpoly_t poly, const char **vars, 
                       const gf216_mpoly_ctx_t ctx);
int gf216_mpoly_equal(const gf216_mpoly_t A, const gf216_mpoly_t B, 
                      const gf216_mpoly_ctx_t ctx);
void gf216_mpoly_randtest(gf216_mpoly_t poly, flint_rand_t state,
                          slong length, slong exp_bound, 
                          const gf216_mpoly_ctx_t ctx);

void fq_nmod_mpoly_to_gf216_mpoly(gf216_mpoly_t res, const fq_nmod_mpoly_t poly,
                                  const fq_nmod_ctx_t fqctx, 
                                  const fq_nmod_mpoly_ctx_t fq_mpoly_ctx);
void gf216_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf216_mpoly_t poly,
                                  const fq_nmod_ctx_t fqctx, 
                                  const fq_nmod_mpoly_ctx_t fq_mpoly_ctx);

/* ============================================================================
   FUNCTION DECLARATIONS - GF(2^32)
   ============================================================================ */

void gf232_mpoly_init(gf232_mpoly_t poly, const gf232_mpoly_ctx_t ctx);
void gf232_mpoly_clear(gf232_mpoly_t poly, const gf232_mpoly_ctx_t ctx);
void gf232_mpoly_ctx_init(gf232_mpoly_ctx_t ctx, slong nvars, const ordering_t ord);
void gf232_mpoly_ctx_clear(gf232_mpoly_ctx_t ctx);
void gf232_mpoly_zero(gf232_mpoly_t poly, const gf232_mpoly_ctx_t ctx);
void gf232_mpoly_set(gf232_mpoly_t res, const gf232_mpoly_t poly, const gf232_mpoly_ctx_t ctx);

void gf232_mpoly_set_coeff_ui_ui(gf232_mpoly_t poly, const gf232_t *c, 
                                 const ulong *exp, const gf232_mpoly_ctx_t ctx);

int gf232_mpoly_mul(gf232_mpoly_t res, const gf232_mpoly_t a, const gf232_mpoly_t b, 
                    const gf232_mpoly_ctx_t ctx);
int gf232_mpoly_mul_array(gf232_mpoly_t A, const gf232_mpoly_t B,
                          const gf232_mpoly_t C, const gf232_mpoly_ctx_t ctx);

int gf232_mpoly_divides(gf232_mpoly_t Q, const gf232_mpoly_t A, 
                        const gf232_mpoly_t B, const gf232_mpoly_ctx_t ctx);

void gf232_mpoly_print(const gf232_mpoly_t poly, const char **vars, 
                       const gf232_mpoly_ctx_t ctx);
int gf232_mpoly_equal(const gf232_mpoly_t A, const gf232_mpoly_t B, 
                      const gf232_mpoly_ctx_t ctx);
void gf232_mpoly_randtest(gf232_mpoly_t poly, flint_rand_t state,
                          slong length, slong exp_bound, 
                          const gf232_mpoly_ctx_t ctx);

void fq_nmod_mpoly_to_gf232_mpoly(gf232_mpoly_t res, const fq_nmod_mpoly_t poly,
                                  const fq_nmod_ctx_t fqctx, 
                                  const fq_nmod_mpoly_ctx_t fq_mpoly_ctx);
void gf232_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf232_mpoly_t poly,
                                  const fq_nmod_ctx_t fqctx, 
                                  const fq_nmod_mpoly_ctx_t fq_mpoly_ctx);

/* ============================================================================
   FUNCTION DECLARATIONS - GF(2^64)
   ============================================================================ */

void gf264_mpoly_init(gf264_mpoly_t poly, const gf264_mpoly_ctx_t ctx);
void gf264_mpoly_clear(gf264_mpoly_t poly, const gf264_mpoly_ctx_t ctx);
void gf264_mpoly_ctx_init(gf264_mpoly_ctx_t ctx, slong nvars, const ordering_t ord);
void gf264_mpoly_ctx_clear(gf264_mpoly_ctx_t ctx);
void gf264_mpoly_zero(gf264_mpoly_t poly, const gf264_mpoly_ctx_t ctx);
void gf264_mpoly_set(gf264_mpoly_t res, const gf264_mpoly_t poly, const gf264_mpoly_ctx_t ctx);

void gf264_mpoly_set_coeff_ui_ui(gf264_mpoly_t poly, const gf264_t *c, 
                                 const ulong *exp, const gf264_mpoly_ctx_t ctx);

int gf264_mpoly_mul(gf264_mpoly_t res, const gf264_mpoly_t a, const gf264_mpoly_t b, 
                    const gf264_mpoly_ctx_t ctx);
int gf264_mpoly_mul_array(gf264_mpoly_t A, const gf264_mpoly_t B,
                          const gf264_mpoly_t C, const gf264_mpoly_ctx_t ctx);

int gf264_mpoly_divides(gf264_mpoly_t Q, const gf264_mpoly_t A, 
                        const gf264_mpoly_t B, const gf264_mpoly_ctx_t ctx);

void gf264_mpoly_print(const gf264_mpoly_t poly, const char **vars, 
                       const gf264_mpoly_ctx_t ctx);
int gf264_mpoly_equal(const gf264_mpoly_t A, const gf264_mpoly_t B, 
                      const gf264_mpoly_ctx_t ctx);
void gf264_mpoly_randtest(gf264_mpoly_t poly, flint_rand_t state,
                          slong length, slong exp_bound, 
                          const gf264_mpoly_ctx_t ctx);

void fq_nmod_mpoly_to_gf264_mpoly(gf264_mpoly_t res, const fq_nmod_mpoly_t poly,
                                  const fq_nmod_ctx_t fqctx, 
                                  const fq_nmod_mpoly_ctx_t fq_mpoly_ctx);
void gf264_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf264_mpoly_t poly,
                                  const fq_nmod_ctx_t fqctx, 
                                  const fq_nmod_mpoly_ctx_t fq_mpoly_ctx);

/* ============================================================================
   FUNCTION DECLARATIONS - GF(2^128)
   ============================================================================ */

void gf2128_mpoly_init(gf2128_mpoly_t poly, const gf2128_mpoly_ctx_t ctx);
void gf2128_mpoly_clear(gf2128_mpoly_t poly, const gf2128_mpoly_ctx_t ctx);
void gf2128_mpoly_ctx_init(gf2128_mpoly_ctx_t ctx, slong nvars, const ordering_t ord);
void gf2128_mpoly_ctx_clear(gf2128_mpoly_ctx_t ctx);
void gf2128_mpoly_zero(gf2128_mpoly_t poly, const gf2128_mpoly_ctx_t ctx);
void gf2128_mpoly_set(gf2128_mpoly_t res, const gf2128_mpoly_t poly, const gf2128_mpoly_ctx_t ctx);

void gf2128_mpoly_set_coeff_ui_ui(gf2128_mpoly_t poly, const gf2128_t *c, 
                                  const ulong *exp, const gf2128_mpoly_ctx_t ctx);

int gf2128_mpoly_mul(gf2128_mpoly_t res, const gf2128_mpoly_t a, const gf2128_mpoly_t b, 
                     const gf2128_mpoly_ctx_t ctx);
int gf2128_mpoly_mul_array(gf2128_mpoly_t A, const gf2128_mpoly_t B,
                           const gf2128_mpoly_t C, const gf2128_mpoly_ctx_t ctx);

int gf2128_mpoly_divides(gf2128_mpoly_t Q, const gf2128_mpoly_t A, 
                         const gf2128_mpoly_t B, const gf2128_mpoly_ctx_t ctx);

void gf2128_mpoly_print(const gf2128_mpoly_t poly, const char **vars, 
                        const gf2128_mpoly_ctx_t ctx);
int gf2128_mpoly_equal(const gf2128_mpoly_t A, const gf2128_mpoly_t B, 
                       const gf2128_mpoly_ctx_t ctx);
void gf2128_mpoly_randtest(gf2128_mpoly_t poly, flint_rand_t state,
                           slong length, slong exp_bound, 
                           const gf2128_mpoly_ctx_t ctx);

void fq_nmod_mpoly_to_gf2128_mpoly(gf2128_mpoly_t res, const fq_nmod_mpoly_t poly,
                                   const fq_nmod_ctx_t fqctx, 
                                   const fq_nmod_mpoly_ctx_t fq_mpoly_ctx);
void gf2128_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf2128_mpoly_t poly,
                                   const fq_nmod_ctx_t fqctx, 
                                   const fq_nmod_mpoly_ctx_t fq_mpoly_ctx);

/* ============================================================================
   TIMING UTILITIES
   ============================================================================ */

double get_wall_time(void);
int gf28_mpoly_test_mul_array(
    const slong *nvars_list, slong nvars_count,
    slong exp_bound, slong nterms, slong ntrials, int verbose);
int gf28_mpoly_test_mul_vs_flint(
    const slong *nvars_list, slong nvars_count,
    slong exp_bound, slong nterms, slong ntrials, int verbose);
#ifdef __cplusplus
}
#endif

#endif /* GF2N_MPOLY_H */
