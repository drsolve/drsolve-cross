/* gf2n_field.h - GF(2^n) field operations with native PCLMUL support */
#ifndef GF2N_FIELD_H
#define GF2N_FIELD_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#if defined(__x86_64__) || defined(__i386__)
#define DIXON_X86_SIMD 1
#include <immintrin.h>
#include <cpuid.h>
#else
#define DIXON_X86_SIMD 0
#endif
#include <flint/flint.h>
#include <flint/fq_nmod.h>
#include <flint/nmod_poly.h>
#include <gmp.h>
#include <time.h>
#include <assert.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 1)
  #ifdef fq_nmod_ctx_prime
  #  undef fq_nmod_ctx_prime
     static inline ulong fq_nmod_ctx_prime(const fq_nmod_ctx_t ctx) {
         return fmpz_get_ui(&((ctx)->p));
     }
  #endif
#endif

#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 2)
  #define flint_rand_init     flint_randinit
  #define flint_rand_set_seed flint_randseed
  #define flint_rand_clear    flint_randclear
#endif

/* ============================================================================
   FORWARD DECLARATIONS AND TYPE DEFINITIONS
   ============================================================================ */

/* GF(2^8) type definition */
typedef uint8_t gf28_t;

/* GF(2^16) type definition */
typedef uint16_t gf216_t;

/* GF(2^32) type definition */
typedef struct {
    uint32_t value;
} gf232_t;

/* GF(2^64) type definition */
typedef struct {
    uint64_t value;
} gf264_t;

/* GF(2^128) type definition */
typedef struct {
    uint64_t low;
    uint64_t high;
} gf2128_t;

/* ============================================================================
   LOOKUP TABLE STRUCTURES
   ============================================================================ */

/* GF(2^8) lookup tables */
typedef struct {
    uint8_t log_table[256];
    uint16_t exp_table[512];
    uint8_t inv_table[256];
    uint8_t mul_table[256][256];
    uint8_t generator; 
} gf28_tables_t;

typedef struct {
    uint8_t mul_table[256][256];
    uint8_t sqr_table[256];
    uint8_t inv_table[256];
    uint8_t double_table[256];
    int initialized;
} gf28_complete_tables_t;

/* GF(2^16) lookup tables */
typedef struct {
    uint16_t log_table[65536];
    uint32_t exp_table[131072];
    uint16_t inv_table[65536];
    uint16_t generator;
} gf216_tables_t;

typedef struct {
    uint16_t *mul_table;
    uint16_t sqr_table[65536];
    uint16_t inv_table[65536];
    uint16_t double_table[65536];
    int initialized;
} gf216_complete_tables_t;

/* Conversion table structures */
typedef struct {
   uint8_t flint_to_gf28[256];
   uint8_t gf28_to_flint[256];
   int initialized;
} gf28_conversion_t;

typedef struct {
   uint16_t *flint_to_gf216;
   uint16_t *gf216_to_flint;
   int initialized;
} gf216_conversion_t;

typedef struct {
    uint32_t flint_to_gf232[256];
    uint32_t gf232_to_flint[256];
    uint32_t flint_poly;
    uint32_t our_poly;
    int initialized;
} gf232_conversion_t;

typedef struct {
    uint64_t *flint_to_gf264_low;
    uint64_t *flint_to_gf264_high;
    uint64_t *gf264_to_flint_low;
    uint64_t *gf264_to_flint_high;
    uint64_t flint_poly_low;
    uint64_t flint_poly_high;
    int initialized;
} gf264_conversion_t;

typedef struct {
    int initialized;
    uint64_t flint_poly_low;
    uint64_t flint_poly_high;
} gf2128_conversion_t;

/* ============================================================================
   FUNCTION POINTER TYPES
   ============================================================================ */

/* Function pointers for multiplication */
typedef gf232_t (*gf232_mul_func)(const gf232_t*, const gf232_t*);
typedef gf264_t (*gf264_mul_func)(const gf264_t*, const gf264_t*);
typedef gf2128_t (*gf2128_mul_func)(const gf2128_t*, const gf2128_t*);

/* ============================================================================
   CPU FEATURE DETECTION
   ============================================================================ */

/* Check CPU features */
int has_pclmulqdq(void);

/* ============================================================================
   GF(2^8) OPERATIONS
   ============================================================================ */

/* Basic operations */
gf28_t gf28_add(gf28_t a, gf28_t b);
gf28_t gf28_sub(gf28_t a, gf28_t b);
gf28_t gf28_mul(gf28_t a, gf28_t b);
gf28_t gf28_sqr(gf28_t a);
gf28_t gf28_inv(gf28_t a);
gf28_t gf28_div(gf28_t a, gf28_t b);

/* Utility functions */
const uint8_t* gf28_get_scalar_row(uint8_t scalar);

/* Initialization and cleanup */
void init_gf28_tables(uint8_t irred_poly);
void init_gf28_complete_tables(void);
void init_gf28_standard(void);
void cleanup_gf28_tables(void);

/* ============================================================================
   GF(2^16) OPERATIONS
   ============================================================================ */

/* Basic operations */
gf216_t gf216_add(gf216_t a, gf216_t b);
gf216_t gf216_mul(gf216_t a, gf216_t b);
gf216_t gf216_sqr(gf216_t a);
gf216_t gf216_inv(gf216_t a);
gf216_t gf216_div(gf216_t a, gf216_t b);

/* Initialization and cleanup */
void init_gf216_tables(uint16_t irred_poly);
void init_gf216_complete_tables(void);
void init_gf216_standard(void);
void cleanup_gf216_tables(void);
void print_gf216_memory_usage(void);

/* ============================================================================
   GF(2^32) OPERATIONS
   ============================================================================ */

/* Basic operations */
gf232_t gf232_create(uint32_t val);
gf232_t gf232_zero(void);
gf232_t gf232_one(void);
int gf232_is_zero(const gf232_t *a);
int gf232_equal(const gf232_t *a, const gf232_t *b);
gf232_t gf232_add(const gf232_t *a, const gf232_t *b);
gf232_t gf232_mul_software(const gf232_t *a, const gf232_t *b);
gf232_t gf232_mul_pclmul(const gf232_t *a, const gf232_t *b);
gf232_t gf232_sqr(const gf232_t *a);
gf232_t gf232_inv(const gf232_t *a);
gf232_t gf232_div(const gf232_t *a, const gf232_t *b);
void gf232_print(const gf232_t *a);

/* Initialization */
void init_gf232(void);

/* ============================================================================
   GF(2^64) OPERATIONS
   ============================================================================ */

/* Basic operations */
gf264_t gf264_create(uint64_t val);
gf264_t gf264_zero(void);
gf264_t gf264_one(void);
int gf264_is_zero(const gf264_t *a);
int gf264_equal(const gf264_t *a, const gf264_t *b);
gf264_t gf264_add(const gf264_t *a, const gf264_t *b);
gf264_t gf264_mul_software(const gf264_t *a, const gf264_t *b);
gf264_t gf264_mul_pclmul(const gf264_t *a, const gf264_t *b);
gf264_t gf264_sqr(const gf264_t *a);
gf264_t gf264_inv(const gf264_t *a);
gf264_t gf264_div(const gf264_t *a, const gf264_t *b);
void gf264_print(const gf264_t *a);

/* Initialization */
void init_gf264(void);

/* ============================================================================
   GF(2^128) OPERATIONS
   ============================================================================ */

/* Basic operations */
gf2128_t gf2128_create(uint64_t low, uint64_t high);
gf2128_t gf2128_zero(void);
gf2128_t gf2128_one(void);
int gf2128_is_zero(const gf2128_t *a);
int gf2128_equal(const gf2128_t *a, const gf2128_t *b);
gf2128_t gf2128_add(const gf2128_t *a, const gf2128_t *b);
gf2128_t gf2128_mul_software(const gf2128_t *a, const gf2128_t *b);
gf2128_t gf2128_mul_clmul(const gf2128_t *a, const gf2128_t *b);
gf2128_t gf2128_sqr(const gf2128_t *a);
gf2128_t gf2128_inv(const gf2128_t *a);
gf2128_t gf2128_div(const gf2128_t *a, const gf2128_t *b);
void gf2128_print(const gf2128_t *a);

/* Initialization */
void init_gf2128(void);

/* ============================================================================
   CONVERSION FUNCTIONS
   ============================================================================ */

/* Helper functions */
int _fq_nmod_ctx_is_gf2n(const fq_nmod_ctx_t ctx);
uint64_t extract_irred_poly(const fq_nmod_ctx_t ctx);
void extract_gf264_poly(const fq_nmod_ctx_t ctx, uint64_t *poly_low, uint64_t *poly_high);

/* GF(2^8) conversions */
void init_gf28_conversion(const fq_nmod_ctx_t ctx);
void cleanup_gf28_conversion(void);
uint8_t fq_nmod_to_gf28_elem(const fq_nmod_t elem, const fq_nmod_ctx_t ctx);
void gf28_elem_to_fq_nmod(fq_nmod_t res, uint8_t elem, const fq_nmod_ctx_t ctx);

/* GF(2^16) conversions */
void init_gf216_conversion(const fq_nmod_ctx_t ctx);
void cleanup_gf216_conversion(void);
uint16_t fq_nmod_to_gf216_elem(const fq_nmod_t elem, const fq_nmod_ctx_t ctx);
void gf216_elem_to_fq_nmod(fq_nmod_t res, uint16_t elem, const fq_nmod_ctx_t ctx);

/* GF(2^32) conversions */
void init_gf232_conversion(const fq_nmod_ctx_t ctx);
void cleanup_gf232_conversion(void);
gf232_t fq_nmod_to_gf232(const fq_nmod_t elem, const fq_nmod_ctx_t ctx);
void gf232_to_fq_nmod(fq_nmod_t res, const gf232_t *elem, const fq_nmod_ctx_t ctx);

/* GF(2^64) conversions */
void init_gf264_conversion(const fq_nmod_ctx_t ctx);
void cleanup_gf264_conversion(void);
gf264_t fq_nmod_to_gf264(const fq_nmod_t elem, const fq_nmod_ctx_t ctx);
void gf264_to_fq_nmod(fq_nmod_t res, const gf264_t *elem, const fq_nmod_ctx_t ctx);

/* GF(2^128) conversions */
void init_gf2128_conversion(const fq_nmod_ctx_t ctx);
gf2128_t fq_nmod_to_gf2128(const fq_nmod_t elem, const fq_nmod_ctx_t ctx);
void gf2128_to_fq_nmod(fq_nmod_t res, const gf2128_t *elem, const fq_nmod_ctx_t ctx);

/* ============================================================================
   UTILITY FUNCTIONS
   ============================================================================ */

void gf2128_to_hex(const gf2128_t *a, char *buf);
int gf2128_from_hex(gf2128_t *a, const char *hex);

/* Global initialization and cleanup */
void init_all_gf2n_fields(void);
void cleanup_all_gf2n_fields(void);

/* Global multiplication function pointers */
extern gf232_mul_func gf232_mul;
extern gf264_mul_func gf264_mul;
extern gf2128_mul_func gf2128_mul;

/* ============================================================================
   GLOBAL VARIABLES
   ============================================================================ */

/* Global lookup tables */
static gf28_tables_t *g_gf28_tables = NULL;
static gf28_complete_tables_t g_gf28_complete_tables = {0};
static gf216_tables_t *g_gf216_tables = NULL;
static gf216_complete_tables_t g_gf216_complete_tables = {0};

/* Conversion tables */
static gf28_conversion_t *g_gf28_conversion = NULL;
static gf216_conversion_t *g_gf216_conversion = NULL;
static gf232_conversion_t *g_gf232_conversion = NULL;
static gf264_conversion_t *g_gf264_conversion = NULL;
static gf2128_conversion_t *g_gf2128_conversion = NULL;
    
#ifdef __cplusplus
}
#endif

#endif /* GF2N_FIELD_H */
