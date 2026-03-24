/* gf2n_poly.h - Polynomial operations over GF(2^n) fields */
#ifndef GF2N_POLY_H
#define GF2N_POLY_H

#include "gf2n_field.h"
#include <flint/fq_nmod_poly.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
   POLYNOMIAL STRUCTURES
   ============================================================================ */

/* GF(2^8) polynomial structure */
typedef struct {
    uint8_t *coeffs;
    slong length;
    slong alloc;
} gf28_poly_struct;
typedef gf28_poly_struct gf28_poly_t[1];

/* GF(2^16) polynomial structure */
typedef struct {
    uint16_t *coeffs;
    slong length;
    slong alloc;
} gf216_poly_struct;
typedef gf216_poly_struct gf216_poly_t[1];

/* GF(2^32) polynomial structure */
typedef struct {
    gf232_t *coeffs;
    slong length;
    slong alloc;
} gf232_poly_struct;
typedef gf232_poly_struct gf232_poly_t[1];

/* GF(2^64) polynomial structure */
typedef struct {
    gf264_t *coeffs;
    slong length;
    slong alloc;
} gf264_poly_struct;
typedef gf264_poly_struct gf264_poly_t[1];

/* GF(2^128) polynomial structure */
typedef struct {
    gf2128_t *coeffs;
    slong length;
    slong alloc;
} gf2128_poly_struct;
typedef gf2128_poly_struct gf2128_poly_t[1];

/* ============================================================================
   MATRIX STRUCTURES
   ============================================================================ */

/* GF(2^8) polynomial matrix structure */
typedef struct {
    gf28_poly_struct *entries;
    slong r, c;
    gf28_poly_struct **rows;
} gf28_poly_mat_struct;
typedef gf28_poly_mat_struct gf28_poly_mat_t[1];

/* GF(2^16) polynomial matrix structure */
typedef struct {
    gf216_poly_struct *entries;
    slong r, c;
    gf216_poly_struct **rows;
} gf216_poly_mat_struct;
typedef gf216_poly_mat_struct gf216_poly_mat_t[1];

/* GF(2^32) polynomial matrix structure */
typedef struct {
    gf232_poly_struct *entries;
    slong r, c;
    gf232_poly_struct **rows;
} gf232_poly_mat_struct;
typedef gf232_poly_mat_struct gf232_poly_mat_t[1];

/* GF(2^64) polynomial matrix structure */
typedef struct {
    gf264_poly_struct *entries;
    slong r, c;
    gf264_poly_struct **rows;
} gf264_poly_mat_struct;
typedef gf264_poly_mat_struct gf264_poly_mat_t[1];

/* GF(2^128) polynomial matrix structure */
typedef struct {
    gf2128_poly_struct *entries;
    slong r, c;
    gf2128_poly_struct **rows;
} gf2128_poly_mat_struct;
typedef gf2128_poly_mat_struct gf2128_poly_mat_t[1];

/* FLINT polynomial matrix structure */
typedef struct {
    fq_nmod_poly_struct *entries;
    slong r;
    slong c;
    fq_nmod_poly_struct **rows;
    fq_nmod_ctx_struct *ctx;
} fq_nmod_poly_mat_struct;
typedef fq_nmod_poly_mat_struct fq_nmod_poly_mat_t[1];

/* ============================================================================
   GF(2^8) POLYNOMIAL OPERATIONS
   ============================================================================ */

/* Basic operations */
void gf28_poly_init(gf28_poly_t poly);
void gf28_poly_clear(gf28_poly_t poly);
void gf28_poly_fit_length(gf28_poly_t poly, slong len);
void gf28_poly_normalise(gf28_poly_t poly);

/* Polynomial properties */
void gf28_poly_zero(gf28_poly_t poly);
int gf28_poly_is_zero(const gf28_poly_t poly);
slong gf28_poly_degree(const gf28_poly_t poly);

/* Coefficient access */
void gf28_poly_set(gf28_poly_t res, const gf28_poly_t poly);
uint8_t gf28_poly_get_coeff(const gf28_poly_t poly, slong i);
void gf28_poly_set_coeff(gf28_poly_t poly, slong i, uint8_t c);

/* Arithmetic operations */
void gf28_poly_add(gf28_poly_t res, const gf28_poly_t a, const gf28_poly_t b);
void gf28_poly_scalar_mul(gf28_poly_t res, const gf28_poly_t poly, uint8_t c);
void gf28_poly_shift_left(gf28_poly_t res, const gf28_poly_t poly, slong n);
void gf28_poly_mul_schoolbook(gf28_poly_t res, const gf28_poly_t a, const gf28_poly_t b);
void gf28_poly_mul_karatsuba(gf28_poly_t res, const gf28_poly_t a, const gf28_poly_t b);
void gf28_poly_mul(gf28_poly_t res, const gf28_poly_t a, const gf28_poly_t b);

/* ============================================================================
   GF(2^16) POLYNOMIAL OPERATIONS
   ============================================================================ */

/* Basic operations */
void gf216_poly_init(gf216_poly_t poly);
void gf216_poly_clear(gf216_poly_t poly);
void gf216_poly_fit_length(gf216_poly_t poly, slong len);
void gf216_poly_normalise(gf216_poly_t poly);

/* Polynomial properties */
void gf216_poly_zero(gf216_poly_t poly);
int gf216_poly_is_zero(const gf216_poly_t poly);
slong gf216_poly_degree(const gf216_poly_t poly);

/* Coefficient access */
void gf216_poly_set(gf216_poly_t res, const gf216_poly_t poly);
void gf216_poly_set_coeff(gf216_poly_t poly, slong i, uint16_t c);
uint16_t gf216_poly_get_coeff(const gf216_poly_t poly, slong i);

/* Arithmetic operations */
void gf216_poly_add(gf216_poly_t res, const gf216_poly_t a, const gf216_poly_t b);
void gf216_poly_scalar_mul(gf216_poly_t res, const gf216_poly_t poly, uint16_t c);
void gf216_poly_shift_left(gf216_poly_t res, const gf216_poly_t poly, slong n);
void gf216_poly_mul(gf216_poly_t res, const gf216_poly_t a, const gf216_poly_t b);

/* ============================================================================
   GF(2^32) POLYNOMIAL OPERATIONS
   ============================================================================ */

/* Basic operations */
void gf232_poly_init(gf232_poly_t poly);
void gf232_poly_clear(gf232_poly_t poly);
void gf232_poly_fit_length(gf232_poly_t poly, slong len);
void gf232_poly_normalise(gf232_poly_t poly);

/* Polynomial properties */
void gf232_poly_zero(gf232_poly_t poly);
int gf232_poly_is_zero(const gf232_poly_t poly);
slong gf232_poly_degree(const gf232_poly_t poly);

/* Coefficient access */
void gf232_poly_set(gf232_poly_t res, const gf232_poly_t poly);
gf232_t gf232_poly_get_coeff(const gf232_poly_t poly, slong i);
void gf232_poly_set_coeff(gf232_poly_t poly, slong i, const gf232_t *c);

/* Arithmetic operations */
void gf232_poly_add(gf232_poly_t res, const gf232_poly_t a, const gf232_poly_t b);
void gf232_poly_scalar_mul(gf232_poly_t res, const gf232_poly_t poly, const gf232_t *c);
void gf232_poly_shift_left(gf232_poly_t res, const gf232_poly_t poly, slong n);
void gf232_poly_mul_schoolbook(gf232_poly_t res, const gf232_poly_t a, const gf232_poly_t b);
void gf232_poly_mul_karatsuba(gf232_poly_t res, const gf232_poly_t a, const gf232_poly_t b);
void gf232_poly_mul(gf232_poly_t res, const gf232_poly_t a, const gf232_poly_t b);

/* ============================================================================
   GF(2^64) POLYNOMIAL OPERATIONS
   ============================================================================ */

/* Basic operations */
void gf264_poly_init(gf264_poly_t poly);
void gf264_poly_clear(gf264_poly_t poly);
void gf264_poly_fit_length(gf264_poly_t poly, slong len);
void gf264_poly_normalise(gf264_poly_t poly);

/* Polynomial properties */
void gf264_poly_zero(gf264_poly_t poly);
int gf264_poly_is_zero(const gf264_poly_t poly);
slong gf264_poly_degree(const gf264_poly_t poly);

/* Coefficient access */
void gf264_poly_set(gf264_poly_t res, const gf264_poly_t poly);
gf264_t gf264_poly_get_coeff(const gf264_poly_t poly, slong i);
void gf264_poly_set_coeff(gf264_poly_t poly, slong i, const gf264_t *c);

/* Arithmetic operations */
void gf264_poly_add(gf264_poly_t res, const gf264_poly_t a, const gf264_poly_t b);
void gf264_poly_scalar_mul(gf264_poly_t res, const gf264_poly_t poly, const gf264_t *c);
void gf264_poly_mul(gf264_poly_t res, const gf264_poly_t a, const gf264_poly_t b);

/* ============================================================================
   GF(2^128) POLYNOMIAL OPERATIONS
   ============================================================================ */

/* Basic operations */
void gf2128_poly_init(gf2128_poly_t poly);
void gf2128_poly_clear(gf2128_poly_t poly);
void gf2128_poly_fit_length(gf2128_poly_t poly, slong len);
void gf2128_poly_normalise(gf2128_poly_t poly);

/* Polynomial properties */
void gf2128_poly_zero(gf2128_poly_t poly);
int gf2128_poly_is_zero(const gf2128_poly_t poly);
slong gf2128_poly_degree(const gf2128_poly_t poly);

/* Coefficient access */
void gf2128_poly_set(gf2128_poly_t res, const gf2128_poly_t poly);
gf2128_t gf2128_poly_get_coeff(const gf2128_poly_t poly, slong i);
void gf2128_poly_set_coeff(gf2128_poly_t poly, slong i, const gf2128_t *c);

/* Arithmetic operations */
void gf2128_poly_add(gf2128_poly_t res, const gf2128_poly_t a, const gf2128_poly_t b);
void gf2128_poly_scalar_mul(gf2128_poly_t res, const gf2128_poly_t poly, const gf2128_t *c);
void gf2128_poly_shift_left(gf2128_poly_t res, const gf2128_poly_t poly, slong n);
void gf2128_poly_mul_schoolbook(gf2128_poly_t res, const gf2128_poly_t a, const gf2128_poly_t b);
void gf2128_poly_mul_karatsuba(gf2128_poly_t res, const gf2128_poly_t a, const gf2128_poly_t b);
void gf2128_poly_mul(gf2128_poly_t res, const gf2128_poly_t a, const gf2128_poly_t b);

/* ============================================================================
   MATRIX OPERATIONS
   ============================================================================ */

/* GF(2^8) matrix operations */
void gf28_poly_mat_init(gf28_poly_mat_t mat, slong rows, slong cols);
void gf28_poly_mat_clear(gf28_poly_mat_t mat);
gf28_poly_struct *gf28_poly_mat_entry(gf28_poly_mat_t mat, slong i, slong j);
void gf28_poly_mat_swap_rows(gf28_poly_mat_t mat, slong r, slong s);

/* GF(2^16) matrix operations */
void gf216_poly_mat_init(gf216_poly_mat_t mat, slong rows, slong cols);
void gf216_poly_mat_clear(gf216_poly_mat_t mat);
gf216_poly_struct *gf216_poly_mat_entry(gf216_poly_mat_t mat, slong i, slong j);
void gf216_poly_mat_swap_rows(gf216_poly_mat_t mat, slong r, slong s);
void gf216_poly_mat_permute_rows(gf216_poly_mat_t mat, const slong *perm);

/* GF(2^32) matrix operations */
void gf232_poly_mat_init(gf232_poly_mat_t mat, slong rows, slong cols);
void gf232_poly_mat_clear(gf232_poly_mat_t mat);
gf232_poly_struct *gf232_poly_mat_entry(gf232_poly_mat_t mat, slong i, slong j);

/* GF(2^64) matrix operations */
void gf264_poly_mat_init(gf264_poly_mat_t mat, slong rows, slong cols);
void gf264_poly_mat_clear(gf264_poly_mat_t mat);
gf264_poly_struct *gf264_poly_mat_entry(gf264_poly_mat_t mat, slong i, slong j);

/* GF(2^128) matrix operations */
void gf2128_poly_mat_init(gf2128_poly_mat_t mat, slong rows, slong cols);
void gf2128_poly_mat_clear(gf2128_poly_mat_t mat);
gf2128_poly_struct *gf2128_poly_mat_entry(gf2128_poly_mat_t mat, slong i, slong j);
const gf2128_poly_struct *gf2128_poly_mat_entry_const(const gf2128_poly_mat_t mat, slong i, slong j);
void gf2128_poly_mat_swap_rows(gf2128_poly_mat_t mat, slong r, slong s);

/* FLINT polynomial matrix operations */
void fq_nmod_poly_mat_init(fq_nmod_poly_mat_t mat, slong rows, slong cols, const fq_nmod_ctx_t ctx);
void fq_nmod_poly_mat_clear(fq_nmod_poly_mat_t mat, const fq_nmod_ctx_t ctx);
fq_nmod_poly_struct* fq_nmod_poly_mat_entry(const fq_nmod_poly_mat_t mat, slong i, slong j);

/* ============================================================================
   CONVERSION FUNCTIONS BETWEEN POLYNOMIALS
   ============================================================================ */

/* Convert between FLINT and optimized polynomial formats */
void fq_nmod_poly_to_gf28_poly(gf28_poly_t res, const fq_nmod_poly_struct *poly, const fq_nmod_ctx_t ctx);
void gf28_poly_to_fq_nmod_poly(fq_nmod_poly_struct *res, const gf28_poly_t poly, const fq_nmod_ctx_t ctx);

void fq_nmod_poly_to_gf216_poly(gf216_poly_t res, const fq_nmod_poly_struct *poly, const fq_nmod_ctx_t ctx);
void gf216_poly_to_fq_nmod_poly(fq_nmod_poly_struct *res, const gf216_poly_t poly, const fq_nmod_ctx_t ctx);

void fq_nmod_poly_to_gf232_poly(gf232_poly_t res, const fq_nmod_poly_struct *poly, const fq_nmod_ctx_t ctx);
void gf232_poly_to_fq_nmod_poly(fq_nmod_poly_struct *res, const gf232_poly_t poly, const fq_nmod_ctx_t ctx);

void fq_nmod_poly_to_gf264_poly(gf264_poly_t res, const fq_nmod_poly_struct *poly, const fq_nmod_ctx_t ctx);
void gf264_poly_to_fq_nmod_poly(fq_nmod_poly_struct *res, const gf264_poly_t poly, const fq_nmod_ctx_t ctx);

void fq_nmod_poly_to_gf2128_poly(gf2128_poly_t res, const fq_nmod_poly_struct *poly, const fq_nmod_ctx_t ctx);
void gf2128_poly_to_fq_nmod_poly(fq_nmod_poly_struct *res, const gf2128_poly_t poly, const fq_nmod_ctx_t ctx);

/* Convert between FLINT and optimized matrix formats */
void fq_nmod_poly_mat_to_gf28(gf28_poly_mat_t res, const fq_nmod_poly_mat_t mat, const fq_nmod_ctx_t ctx);
void fq_nmod_poly_mat_to_gf216(gf216_poly_mat_t res, const fq_nmod_poly_mat_t mat, const fq_nmod_ctx_t ctx);
void fq_nmod_poly_mat_to_gf2128(gf2128_poly_mat_t res, const fq_nmod_poly_mat_t mat, const fq_nmod_ctx_t ctx);
void gf2128_poly_mat_to_fq_nmod(fq_nmod_poly_mat_t res, const gf2128_poly_mat_t mat, const fq_nmod_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif /* GF2N_POLY_H */
