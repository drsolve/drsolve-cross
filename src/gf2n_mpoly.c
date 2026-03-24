/* gf2n_mpoly.c - Unified GF(2^n) Multivariate Polynomial Implementation
 * Up to ~9+ variables with 8-bit exponent packing under LEX ordering
 */

#include "gf2n_mpoly.h"
#if DIXON_X86_SIMD
#include <emmintrin.h>  /* SSE2 */
#endif

/* ============================================================================
   HELPER: multi-word packed exponent field access (N <= 2)
   ============================================================================ */

FLINT_INLINE ulong
_mpoly_get_field(const ulong *exps, flint_bitcnt_t bits, slong k)
{
    slong fbit = (slong)bits * k;
    slong word = fbit / FLINT_BITS;
    slong boff = fbit % FLINT_BITS;
    ulong mask = (bits < (flint_bitcnt_t)FLINT_BITS)
               ? (((ulong)1 << bits) - 1)
               : ~(ulong)0;
    return (exps[word] >> boff) & mask;
}

FLINT_INLINE void
_mpoly_set_field2(ulong exp[2], flint_bitcnt_t bits, slong k, ulong val)
{
    slong fbit = (slong)bits * k;
    slong word = fbit / FLINT_BITS;
    slong boff = fbit % FLINT_BITS;
    if (word < 2)
        exp[word] |= val << boff;
}

FLINT_INLINE void _sub2(ulong dst[2], const ulong dec[2])
{
    ulong old = dst[0];
    dst[0] -= dec[0];
    if (dst[0] > old)
        dst[1]--;
    dst[1] -= dec[1];
}

/* ============================================================================
   MAIN-VARIABLE SPLIT FOR N = 1 OR 2
   ============================================================================ */
/*
 * FLINT LEX field convention (critical for correctness):
 *
 *   For ORD_LEX with nvars variables x_0 > x_1 > ... > x_{n-1}:
 *     field k  =  x_{nvars-1-k}   (reversed!)
 *     field 0  =  x_{n-1}  at bits [0 .. bits-1]       (least significant in LEX)
 *     field n-1=  x_0      at bits [(n-1)*bits .. n*bits-1] (most significant in LEX)
 *
 *   num = nfields-1 = nvars-1  =>  main variable is FIELD num = x_0 (most significant).
 *   Sub-variables covered by the flat array: fields 0..num-1 = x_{n-1}..x_1.
 *
 *   mpoly_max_fields_fmpz gives maxBf[k] = max of field k = max of x_{n-1-k}.
 *   So mults[k] = 1 + maxBf[k]_B + maxBf[k]_C  refers to x_{n-1-k}.
 *
 *   The flat index for a term is:
 *     flat = sum_{j=0}^{num-1}  field_j * stride_j
 *          = x_{n-1}*1 + x_{n-2}*mults[0] + ... + x_1*mults[0]*...*mults[n-3]
 *
 *   This matches how mults are ordered by _gf28_mpoly_mul_array_LEX (mults[0]
 *   is the range for field 0 = x_{n-1}, mults[n-2] for field n-2 = x_1).
 *
 *   Terms in LEX descending order => main_val = field num = x_0 exponent
 *   is non-increasing as i increases. The monotone boundary-fill is valid.
 */
static void
_mpoly_main_variable_split_LEX_N(
    slong          *Amain,
    ulong          *Apexp,
    const ulong    *Aexp,
    slong           Alength,
    slong           Atotal,
    const ulong    *mults,
    slong           num,
    flint_bitcnt_t  Abits,
    slong           N,
    const mpoly_ctx_t minfo)
{
    slong nvars = minfo->nvars;
    ulong *expv = (ulong *)flint_malloc(nvars * sizeof(ulong));

    for (slong i = 0; i < Atotal; i++) {
        mpoly_get_monomial_ui(expv, Aexp + N * i, Abits, minfo);
        ulong flat = 0, stride = 1;
        for (slong j = 0; j < num; j++) {
            flat   += expv[nvars - 1 - j] * stride;
            stride *= mults[j];
        }
        Apexp[i] = flat;
    }

    for (slong k = 0; k <= Alength; k++)
        Amain[k] = Atotal;
    
    slong cur_k = 0;
    for (slong i = 0; i < Atotal; i++) {
        mpoly_get_monomial_ui(expv, Aexp + N * i, Abits, minfo);
        slong main_val = (slong)expv[0]; 
        slong group_k  = Alength - 1 - main_val;
        while (cur_k <= group_k) {
            Amain[cur_k] = i;
            cur_k++;
        }
    }

    flint_free(expv);
}

/* ============================================================================
   MACRO TEMPLATE: CORE ARRAY MULTIPLICATION (supports N = 1 and N = 2)
   ============================================================================ */

#define DEFINE_MPOLY_ARRAY_MUL(FIELD, COEFF_T, IS_ZERO, ADD_OP, MUL_OP, ARRAY_LIMIT) \
\
static void _##FIELD##_mpoly_addmul_array1_safe(                               \
    COEFF_T *poly1, slong array_size,                                          \
    const COEFF_T *poly2, const ulong *exp2, slong len2,                       \
    const COEFF_T *poly3, const ulong *exp3, slong len3)                       \
{                                                                              \
    for (slong ii = 0; ii < len2 + BLOCK; ii += BLOCK) {                      \
        for (slong jj = 0; jj < len3 + BLOCK; jj += BLOCK) {                 \
            for (slong i = ii; i < FLINT_MIN(ii + BLOCK, len2); i++) {        \
                slong off2 = (slong)exp2[i];                                   \
                if (off2 >= array_size) continue;                              \
                COEFF_T *c2 = poly1 + off2;                                    \
                if (!(IS_ZERO(poly2[i]))) {                                    \
                    for (slong j = jj; j < FLINT_MIN(jj + BLOCK, len3); j++) {\
                        slong off3 = (slong)exp3[j];                           \
                        if (off2 + off3 >= array_size) continue;               \
                        COEFF_T prod = MUL_OP(poly2[i], poly3[j]);             \
                        c2[off3] = ADD_OP(c2[off3], prod);                     \
                    }                                                          \
                }                                                              \
            }                                                                  \
        }                                                                      \
    }                                                                          \
}                                                                              \
\
static void _##FIELD##_mpoly_repack_exps(                                      \
    FIELD##_mpoly_t poly, flint_bitcnt_t new_bits,                             \
    const FIELD##_mpoly_ctx_t ctx)                                             \
{                                                                              \
    if (poly->bits == new_bits) return;                                        \
    if (poly->length == 0) { poly->bits = new_bits; return; }                 \
    slong nvars = ctx->minfo->nvars;                                           \
    slong N_old = mpoly_words_per_exp(poly->bits, ctx->minfo);                \
    slong N_new = mpoly_words_per_exp(new_bits,   ctx->minfo);                \
    ulong *ne   = (ulong *)flint_malloc(N_new * poly->length * sizeof(ulong));\
    ulong *expv = (ulong *)flint_malloc(nvars * sizeof(ulong));                \
    for (slong ii = 0; ii < poly->length; ii++) {                              \
        mpoly_get_monomial_ui(expv, poly->exps + N_old*ii, poly->bits, ctx->minfo); \
        mpoly_set_monomial_ui(ne   + N_new*ii, expv, new_bits, ctx->minfo);   \
    }                                                                          \
    flint_free(expv);                                                          \
    flint_free(poly->exps);                                                    \
    poly->exps       = ne;                                                     \
    poly->exps_alloc = N_new * poly->length;                                   \
    poly->bits       = new_bits;                                               \
}                                                                              \
\
static slong FIELD##_mpoly_append_array_LEX_safe(                             \
    FIELD##_mpoly_t P, slong Plen, COEFF_T *coeff_array,                      \
    const ulong *mults, slong num, slong array_size,                          \
    slong top, const FIELD##_mpoly_ctx_t ctx)                                 \
{                                                                              \
    slong N     = mpoly_words_per_exp(P->bits, ctx->minfo);                   \
    slong nvars = ctx->minfo->nvars;                                           \
    ulong *expv = (ulong *)flint_malloc(nvars * sizeof(ulong));                \
    ulong *pe   = (ulong *)flint_malloc(N     * sizeof(ulong));                \
    for (slong off = array_size - 1; off >= 0; off--) {                       \
        if (IS_ZERO(coeff_array[off])) continue;                              \
        COEFF_T coeff    = coeff_array[off];                                   \
        slong d = off;                                                         \
        for (slong j = 0; j < num; j++) {                                     \
            expv[nvars - 1 - j] = (ulong)(d % mults[j]);                      \
            d /= mults[j];                                                     \
        }                                                                      \
        expv[0] = (ulong)top;                                                  \
        mpoly_set_monomial_ui(pe, expv, P->bits, ctx->minfo);                 \
        if (Plen >= P->coeffs_alloc) {                                        \
            slong na = FLINT_MAX(Plen + 1, 2 * P->coeffs_alloc);              \
            P->coeffs = (COEFF_T *)flint_realloc(P->coeffs, na*sizeof(COEFF_T));\
            P->coeffs_alloc = na;                                              \
        }                                                                      \
        if (N * (Plen + 1) > P->exps_alloc) {                                \
            slong na = FLINT_MAX(N*(Plen+1), 2*P->exps_alloc);                \
            P->exps = (ulong *)flint_realloc(P->exps, na*sizeof(ulong));      \
            P->exps_alloc = na;                                                \
        }                                                                      \
        memcpy(P->exps + N*Plen, pe, N*sizeof(ulong));                        \
        P->coeffs[Plen] = coeff;                                               \
        Plen++;                                                                \
    }                                                                          \
    flint_free(expv); flint_free(pe);                                          \
    return Plen;                                                               \
}                                                                              \
\
static void _##FIELD##_mpoly_mul_array_chunked_LEX(                           \
    FIELD##_mpoly_t P,                                                         \
    const FIELD##_mpoly_t A, const FIELD##_mpoly_t B,                         \
    const ulong *mults, const FIELD##_mpoly_ctx_t ctx)                        \
{                                                                              \
    slong num   = ctx->minfo->nfields - 1;                                     \
    slong nvars = ctx->minfo->nvars;                                           \
    slong NA    = mpoly_words_per_exp(A->bits, ctx->minfo);                   \
    slong NB    = mpoly_words_per_exp(B->bits, ctx->minfo);                   \
    slong array_size = 1;                                                      \
    for (slong i = 0; i < num; i++) array_size *= mults[i];                  \
    ulong *_tv = (ulong *)flint_malloc(nvars * sizeof(ulong));                 \
    mpoly_get_monomial_ui(_tv, A->exps, A->bits, ctx->minfo);                 \
    slong Al = 1 + (slong)_tv[0];                                              \
    mpoly_get_monomial_ui(_tv, B->exps, B->bits, ctx->minfo);                 \
    slong Bl = 1 + (slong)_tv[0];                                              \
    flint_free(_tv);                                                           \
    slong *Amain = (slong *)flint_malloc((Al + 1) * sizeof(slong));           \
    slong *Bmain = (slong *)flint_malloc((Bl + 1) * sizeof(slong));           \
    ulong *Apexp = (ulong *)flint_malloc(A->length * sizeof(ulong));           \
    ulong *Bpexp = (ulong *)flint_malloc(B->length * sizeof(ulong));           \
    _mpoly_main_variable_split_LEX_N(Amain, Apexp, A->exps, Al,              \
        A->length, mults, num, A->bits, NA, ctx->minfo);                       \
    _mpoly_main_variable_split_LEX_N(Bmain, Bpexp, B->exps, Bl,              \
        B->length, mults, num, B->bits, NB, ctx->minfo);                       \
    slong Pl = Al + Bl - 1, Plen = 0;                                         \
    COEFF_T *coeff_array = (COEFF_T *)calloc(array_size, sizeof(COEFF_T));    \
    if (!coeff_array) goto done;                                               \
    for (slong Pi = 0; Pi < Pl; Pi++) {                                       \
        memset(coeff_array, 0, array_size * sizeof(COEFF_T));                 \
        for (slong i = 0, j = Pi; i < Al && j >= 0; i++, j--) {              \
            if (j < Bl)                                                        \
                _##FIELD##_mpoly_addmul_array1_safe(coeff_array, array_size,  \
                    A->coeffs + Amain[i], Apexp + Amain[i],                   \
                    Amain[i+1] - Amain[i],                                    \
                    B->coeffs + Bmain[j], Bpexp + Bmain[j],                   \
                    Bmain[j+1] - Bmain[j]);                                   \
        }                                                                      \
        Plen = FIELD##_mpoly_append_array_LEX_safe(P, Plen, coeff_array,     \
            mults, num, array_size, Pl - Pi - 1, ctx);                        \
    }                                                                          \
    P->length = Plen;                                                          \
    free(coeff_array);                                                         \
done:                                                                          \
    flint_free(Amain); flint_free(Bmain);                                      \
    flint_free(Apexp); flint_free(Bpexp);                                      \
}                                                                              \
\
static int _##FIELD##_mpoly_mul_array_LEX(                                    \
    FIELD##_mpoly_t A,                                                         \
    const FIELD##_mpoly_t B, fmpz *maxBfields,                                \
    const FIELD##_mpoly_t C, fmpz *maxCfields,                                \
    const FIELD##_mpoly_ctx_t ctx)                                             \
{                                                                              \
    FLINT_ASSERT(ctx->minfo->nvars > 0);                                       \
    FLINT_ASSERT(B->length != 0 && C->length != 0);                           \
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);                                 \
    slong i = ctx->minfo->nfields - 1;                                         \
    ulong *mults = (ulong *)flint_malloc(ctx->minfo->nfields * sizeof(ulong));\
    mults[i] = 1 + fmpz_get_ui(maxBfields+i) + fmpz_get_ui(maxCfields+i);    \
    ulong max = mults[i];                                                      \
    if (((slong)mults[i]) <= 0) { flint_free(mults); return 0; }              \
    slong array_size = WORD(1);                                                \
    for (i--; i >= 0; i--) {                                                  \
        mults[i] = 1 + fmpz_get_ui(maxBfields+i) + fmpz_get_ui(maxCfields+i);\
        max |= mults[i];                                                       \
        ulong hi;                                                              \
        umul_ppmm(hi, array_size, array_size, mults[i]);                      \
        if (hi != 0 || (slong)mults[i] <= 0 || array_size <= 0)              \
            { flint_free(mults); return 0; }                                   \
    }                                                                          \
    if (array_size > ARRAY_LIMIT) { flint_free(mults); return 0; }            \
    flint_bitcnt_t exp_bits = FLINT_MAX(MPOLY_MIN_BITS,                       \
        FLINT_BIT_COUNT(max) + 1);                                             \
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);                          \
    if (mpoly_words_per_exp(exp_bits, ctx->minfo) > 2)                        \
        { flint_free(mults); return 0; }                                       \
    FIELD##_mpoly_t Btmp, Ctmp;                                                \
    FIELD##_mpoly_init(Btmp, ctx);                                             \
    FIELD##_mpoly_init(Ctmp, ctx);                                             \
    FIELD##_mpoly_set(Btmp, B, ctx);                                           \
    FIELD##_mpoly_set(Ctmp, C, ctx);                                           \
    _##FIELD##_mpoly_repack_exps(Btmp, exp_bits, ctx);                        \
    _##FIELD##_mpoly_repack_exps(Ctmp, exp_bits, ctx);                        \
    if (A == B || A == C) {                                                    \
        FIELD##_mpoly_t T;                                                     \
        FIELD##_mpoly_init(T, ctx);                                            \
        T->bits = exp_bits;                                                    \
        _##FIELD##_mpoly_mul_array_chunked_LEX(T, Ctmp, Btmp, mults, ctx);    \
        FIELD##_mpoly_struct tmp = *T; *T = *A; *A = tmp;                     \
        FIELD##_mpoly_clear(T, ctx);                                           \
    } else {                                                                   \
        A->bits = exp_bits;                                                    \
        _##FIELD##_mpoly_mul_array_chunked_LEX(A, Ctmp, Btmp, mults, ctx);    \
    }                                                                          \
    FIELD##_mpoly_clear(Btmp, ctx);                                            \
    FIELD##_mpoly_clear(Ctmp, ctx);                                            \
    flint_free(mults);                                                         \
    return 1;                                                                  \
}

/* ============================================================================
   INSTANTIATE ARRAY MUL FOR EACH FIELD
   ============================================================================ */

#define GF28_IS_ZERO(x)   ((x) == 0)
#define GF28_ADD(a, b)    ((a) ^ (b))
#define GF28_MUL(a, b)    (gf28_mul((a), (b)))

#define GF216_IS_ZERO(x)  ((x) == 0)
#define GF216_ADD(a, b)   ((a) ^ (b))
#define GF216_MUL(a, b)   (gf216_mul((a), (b)))
DEFINE_MPOLY_ARRAY_MUL(gf216, uint16_t, GF216_IS_ZERO, GF216_ADD, GF216_MUL, (1L << 27))

#define GF232_IS_ZERO(x)  (gf232_is_zero(&(x)))
#define GF232_ADD(a, b)   (gf232_add(&(a), &(b)))
#define GF232_MUL(a, b)   (gf232_mul(&(a), &(b)))
DEFINE_MPOLY_ARRAY_MUL(gf232, gf232_t, GF232_IS_ZERO, GF232_ADD, GF232_MUL, (1L << 27))

#define GF264_IS_ZERO(x)  (gf264_is_zero(&(x)))
#define GF264_ADD(a, b)   (gf264_add(&(a), &(b)))
#define GF264_MUL(a, b)   (gf264_mul(&(a), &(b)))
DEFINE_MPOLY_ARRAY_MUL(gf264, gf264_t, GF264_IS_ZERO, GF264_ADD, GF264_MUL, (1L << 26))

#define GF2128_IS_ZERO(x) (gf2128_is_zero(&(x)))
#define GF2128_ADD(a, b)  (gf2128_add(&(a), &(b)))
#define GF2128_MUL(a, b)  (gf2128_mul(&(a), &(b)))
DEFINE_MPOLY_ARRAY_MUL(gf2128, gf2128_t, GF2128_IS_ZERO, GF2128_ADD, GF2128_MUL, (1L << 26))

/* ============================================================================
   DIVISION IMPLEMENTATION
   ============================================================================ */

#define DEFINE_MPOLY_DIVISION_IMPL(FIELD, COEFF_T, IS_ZERO, ADD_OP, INV_OP, MUL_OP, SET_COEFF_CALL) \
\
static void FIELD##_mpolyd_init(FIELD##_mpolyd_t A, slong nvars) {            \
    A->coeffs = NULL; A->alloc = 0;                                            \
    A->deg_bounds = (slong *)calloc(nvars, sizeof(slong));                     \
    A->nvars = nvars;                                                          \
}                                                                              \
static void FIELD##_mpolyd_clear(FIELD##_mpolyd_t A) {                        \
    if (A->coeffs)     free(A->coeffs);                                        \
    if (A->deg_bounds) free(A->deg_bounds);                                    \
}                                                                              \
static slong FIELD##_mpolyd_offset(const FIELD##_mpolyd_t A, const ulong *exp) { \
    slong off = 0, stride = 1;                                                 \
    for (slong i = 0; i < A->nvars; i++) {                                    \
        if (exp[i] >= (ulong)A->deg_bounds[i]) return -1;                     \
        off += (slong)exp[i] * stride;                                         \
        stride *= A->deg_bounds[i];                                            \
    }                                                                          \
    return off;                                                                \
}                                                                              \
static int FIELD##_mpolyd_set_degbounds(FIELD##_mpolyd_t A, const slong *bounds) { \
    slong size = 1;                                                            \
    for (slong i = 0; i < A->nvars; i++) {                                    \
        A->deg_bounds[i] = bounds[i];                                          \
        if (bounds[i] <= 0) return 0;                                          \
        if (size > WORD_MAX / bounds[i]) return 0;                             \
        size *= bounds[i];                                                     \
    }                                                                          \
    if (size > (1L << 26)) return 0;                                           \
    A->alloc = size;                                                           \
    A->coeffs = (COEFF_T *)calloc(size, sizeof(COEFF_T));                     \
    return A->coeffs != NULL;                                                  \
}                                                                              \
static void FIELD##_mpoly_to_mpolyd(FIELD##_mpolyd_t A,                       \
    const FIELD##_mpoly_t B, const FIELD##_mpoly_ctx_t ctx)                   \
{                                                                              \
    slong nvars = ctx->minfo->nvars;                                           \
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);                       \
    ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));                 \
    memset(A->coeffs, 0, A->alloc * sizeof(COEFF_T));                         \
    for (slong i = 0; i < B->length; i++) {                                   \
        mpoly_get_monomial_ui(exp, B->exps + N*i, B->bits, ctx->minfo);       \
        slong off = FIELD##_mpolyd_offset(A, exp);                            \
        if (off >= 0 && off < A->alloc) A->coeffs[off] = B->coeffs[i];       \
    }                                                                          \
    flint_free(exp);                                                           \
}                                                                              \
static void FIELD##_mpolyd_to_mpoly(FIELD##_mpoly_t A,                        \
    const FIELD##_mpolyd_t B, const FIELD##_mpoly_ctx_t ctx)                  \
{                                                                              \
    slong nvars = ctx->minfo->nvars;                                           \
    ulong *exp = (ulong *)calloc(nvars, sizeof(ulong));                        \
    flint_bitcnt_t bits_needed = MPOLY_MIN_BITS;                               \
    for (slong i = 0; i < nvars; i++) {                                        \
        if (B->deg_bounds[i] > 0)                                              \
            bits_needed = FLINT_MAX(bits_needed,                               \
                (flint_bitcnt_t)FLINT_BIT_COUNT(B->deg_bounds[i] - 1));       \
    }                                                                          \
    bits_needed = mpoly_fix_bits(bits_needed, ctx->minfo);                    \
    FIELD##_mpoly_zero(A, ctx); A->bits = bits_needed;                         \
    for (slong off = B->alloc - 1; off >= 0; off--) {                         \
        if (IS_ZERO(B->coeffs[off])) continue;                                 \
        slong temp = off;                                                      \
        for (slong i = 0; i < nvars; i++) {                                    \
            exp[i] = (ulong)(temp % B->deg_bounds[i]);                         \
            temp  /= B->deg_bounds[i];                                         \
        }                                                                      \
        SET_COEFF_CALL(A, B->coeffs[off], exp, ctx);                          \
    }                                                                          \
    free(exp);                                                                 \
}                                                                              \
static void FIELD##_mpolyd_divrem_univar(FIELD##_mpolyd_t Q,                  \
    FIELD##_mpolyd_t R, const FIELD##_mpolyd_t A, const FIELD##_mpolyd_t B)   \
{                                                                              \
    slong n = A->alloc;                                                        \
    memcpy(R->coeffs, A->coeffs, n * sizeof(COEFF_T));                        \
    memset(Q->coeffs, 0, Q->alloc * sizeof(COEFF_T));                         \
    slong degB = -1;                                                           \
    for (slong i = n - 1; i >= 0; i--)                                        \
        if (!(IS_ZERO(B->coeffs[i]))) { degB = i; break; }                    \
    if (degB < 0) return;                                                      \
    COEFF_T lc_B_inv = INV_OP(B->coeffs[degB]);                               \
    for (slong i = n - 1; i >= degB; i--) {                                   \
        if (!(IS_ZERO(R->coeffs[i]))) {                                        \
            COEFF_T q = MUL_OP(R->coeffs[i], lc_B_inv);                       \
            if (i - degB < Q->alloc) Q->coeffs[i - degB] = q;                \
            for (slong j = 0; j <= degB && i - degB + j < n; j++) {           \
                COEFF_T prod = MUL_OP(q, B->coeffs[j]);                       \
                R->coeffs[i-degB+j] = ADD_OP(R->coeffs[i-degB+j], prod);     \
            }                                                                  \
        }                                                                      \
    }                                                                          \
}                                                                              \
static int FIELD##_mpolyd_is_zero(const FIELD##_mpolyd_t A) {                 \
    for (slong i = 0; i < A->alloc; i++)                                       \
        if (!(IS_ZERO(A->coeffs[i]))) return 0;                                \
    return 1;                                                                  \
}                                                                              \
static int FIELD##_mpoly_divides_dense(FIELD##_mpoly_t Q,                     \
    const FIELD##_mpoly_t A, const FIELD##_mpoly_t B,                         \
    const FIELD##_mpoly_ctx_t ctx)                                             \
{                                                                              \
    slong nvars = ctx->minfo->nvars;                                           \
    slong *degs_A  = (slong *)flint_malloc(nvars * sizeof(slong));             \
    slong *degs_B  = (slong *)flint_malloc(nvars * sizeof(slong));             \
    slong *bounds  = (slong *)flint_malloc(nvars * sizeof(slong));             \
    mpoly_degrees_si(degs_A, A->exps, A->length, A->bits, ctx->minfo);        \
    mpoly_degrees_si(degs_B, B->exps, B->length, B->bits, ctx->minfo);        \
    for (slong i = 0; i < nvars; i++) {                                        \
        if (degs_A[i] < degs_B[i]) {                                           \
            flint_free(degs_A); flint_free(degs_B); flint_free(bounds);        \
            FIELD##_mpoly_zero(Q, ctx); return 0;                              \
        }                                                                      \
        bounds[i] = degs_A[i] + 1;                                             \
    }                                                                          \
    FIELD##_mpolyd_t Ad, Bd, Qd, Rd;                                          \
    FIELD##_mpolyd_init(Ad, nvars); FIELD##_mpolyd_init(Bd, nvars);           \
    FIELD##_mpolyd_init(Qd, nvars); FIELD##_mpolyd_init(Rd, nvars);           \
    int success = 0;                                                           \
    if (!FIELD##_mpolyd_set_degbounds(Ad, bounds) ||                           \
        !FIELD##_mpolyd_set_degbounds(Bd, bounds) ||                           \
        !FIELD##_mpolyd_set_degbounds(Qd, bounds) ||                           \
        !FIELD##_mpolyd_set_degbounds(Rd, bounds)) goto cleanup;               \
    FIELD##_mpoly_to_mpolyd(Ad, A, ctx);                                       \
    FIELD##_mpoly_to_mpolyd(Bd, B, ctx);                                       \
    FIELD##_mpolyd_divrem_univar(Qd, Rd, Ad, Bd);                             \
    if (FIELD##_mpolyd_is_zero(Rd)) {                                          \
        FIELD##_mpolyd_to_mpoly(Q, Qd, ctx); success = 1;                     \
    } else {                                                                   \
        FIELD##_mpoly_zero(Q, ctx); success = 0;                               \
    }                                                                          \
cleanup:                                                                       \
    FIELD##_mpolyd_clear(Ad); FIELD##_mpolyd_clear(Bd);                       \
    FIELD##_mpolyd_clear(Qd); FIELD##_mpolyd_clear(Rd);                       \
    flint_free(degs_A); flint_free(degs_B); flint_free(bounds);                \
    return success;                                                            \
}

DEFINE_MPOLY_DIVISION_IMPL(gf216, uint16_t, GF216_IS_ZERO, GF216_ADD,
    gf216_inv, GF216_MUL, gf216_mpoly_set_coeff_ui_ui)

#define GF232_INV(x)          (gf232_inv(&(x)))
#define GF232_SET_COEFF(A,c,e,ctx) gf232_mpoly_set_coeff_ui_ui(A,&(c),e,ctx)
DEFINE_MPOLY_DIVISION_IMPL(gf232, gf232_t, GF232_IS_ZERO, GF232_ADD,
    GF232_INV, GF232_MUL, GF232_SET_COEFF)

#define GF264_INV(x)          (gf264_inv(&(x)))
#define GF264_SET_COEFF(A,c,e,ctx) gf264_mpoly_set_coeff_ui_ui(A,&(c),e,ctx)
DEFINE_MPOLY_DIVISION_IMPL(gf264, gf264_t, GF264_IS_ZERO, GF264_ADD,
    GF264_INV, GF264_MUL, GF264_SET_COEFF)

#define GF2128_INV(x)         (gf2128_inv(&(x)))
#define GF2128_SET_COEFF(A,c,e,ctx) gf2128_mpoly_set_coeff_ui_ui(A,&(c),e,ctx)
DEFINE_MPOLY_DIVISION_IMPL(gf2128, gf2128_t, GF2128_IS_ZERO, GF2128_ADD,
    GF2128_INV, GF2128_MUL, GF2128_SET_COEFF)

/* ============================================================================
   PUBLIC API - BASIC OPERATIONS
   ============================================================================ */

double get_wall_time(void) {
    struct timeval time;
    if (gettimeofday(&time, NULL)) return 0;
    return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
}

#define DEFINE_BASIC_OPS(FIELD, COEFF_T, IS_ZERO_CHECK)                        \
void FIELD##_mpoly_init(FIELD##_mpoly_t poly, const FIELD##_mpoly_ctx_t ctx) { \
    poly->coeffs = NULL; poly->exps = NULL; poly->length = 0;                  \
    poly->coeffs_alloc = 0; poly->exps_alloc = 0; poly->bits = MPOLY_MIN_BITS; \
}                                                                               \
void FIELD##_mpoly_clear(FIELD##_mpoly_t poly, const FIELD##_mpoly_ctx_t ctx) {\
    if (poly->coeffs) flint_free(poly->coeffs);                                 \
    if (poly->exps)   flint_free(poly->exps);                                   \
}                                                                               \
void FIELD##_mpoly_ctx_init(FIELD##_mpoly_ctx_t ctx, slong nvars,              \
    const ordering_t ord) { mpoly_ctx_init(ctx->minfo, nvars, ord); }          \
void FIELD##_mpoly_ctx_clear(FIELD##_mpoly_ctx_t ctx) {                        \
    mpoly_ctx_clear(ctx->minfo); }                                              \
void FIELD##_mpoly_zero(FIELD##_mpoly_t poly, const FIELD##_mpoly_ctx_t ctx) { \
    poly->length = 0; }                                                         \
void FIELD##_mpoly_set(FIELD##_mpoly_t res, const FIELD##_mpoly_t poly,        \
    const FIELD##_mpoly_ctx_t ctx)                                              \
{                                                                               \
    if (res == poly) return;                                                    \
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);                     \
    if (res->coeffs_alloc < poly->length) {                                     \
        res->coeffs = (COEFF_T *)flint_realloc(res->coeffs,                    \
            poly->length * sizeof(COEFF_T));                                    \
        res->coeffs_alloc = poly->length;                                       \
    }                                                                           \
    if (res->exps_alloc < N * poly->length) {                                   \
        res->exps = (ulong *)flint_realloc(res->exps,                           \
            N * poly->length * sizeof(ulong));                                  \
        res->exps_alloc = N * poly->length;                                     \
    }                                                                           \
    memcpy(res->coeffs, poly->coeffs, poly->length * sizeof(COEFF_T));         \
    memcpy(res->exps, poly->exps, N * poly->length * sizeof(ulong));           \
    res->length = poly->length; res->bits = poly->bits;                        \
}                                                                               \
int FIELD##_mpoly_equal(const FIELD##_mpoly_t A, const FIELD##_mpoly_t B,      \
    const FIELD##_mpoly_ctx_t ctx)                                              \
{                                                                               \
    if (A->length != B->length) return 0;                                       \
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);                        \
    for (slong i = 0; i < A->length; i++) {                                     \
        if (IS_ZERO_CHECK(A->coeffs[i], B->coeffs[i])) return 0;               \
        if (!mpoly_monomial_equal(A->exps + N*i, B->exps + N*i, N)) return 0;  \
    }                                                                           \
    return 1;                                                                   \
}                                                                               \
int FIELD##_mpoly_mul_array(FIELD##_mpoly_t A, const FIELD##_mpoly_t B,        \
    const FIELD##_mpoly_t C, const FIELD##_mpoly_ctx_t ctx)                    \
{                                                                               \
    if (B->length == 0 || C->length == 0) { FIELD##_mpoly_zero(A,ctx); return 1; } \
    fmpz *maxBf = (fmpz *)flint_malloc(ctx->minfo->nfields * sizeof(fmpz));    \
    fmpz *maxCf = (fmpz *)flint_malloc(ctx->minfo->nfields * sizeof(fmpz));    \
    for (slong i = 0; i < ctx->minfo->nfields; i++) {                          \
        fmpz_init(maxBf+i); fmpz_init(maxCf+i); }                              \
    mpoly_max_fields_fmpz(maxBf, B->exps, B->length, B->bits, ctx->minfo);    \
    mpoly_max_fields_fmpz(maxCf, C->exps, C->length, C->bits, ctx->minfo);    \
    int success = 0;                                                            \
    if (ctx->minfo->ord == ORD_LEX)                                             \
        success = _##FIELD##_mpoly_mul_array_LEX(A,B,maxBf,C,maxCf,ctx);      \
    for (slong i = 0; i < ctx->minfo->nfields; i++) {                          \
        fmpz_clear(maxBf+i); fmpz_clear(maxCf+i); }                            \
    flint_free(maxBf); flint_free(maxCf);                                       \
    return success;                                                             \
}                                                                               \
int FIELD##_mpoly_mul(FIELD##_mpoly_t res, const FIELD##_mpoly_t a,            \
    const FIELD##_mpoly_t b, const FIELD##_mpoly_ctx_t ctx)                    \
{   return FIELD##_mpoly_mul_array(res, a, b, ctx); }                          \
int FIELD##_mpoly_divides(FIELD##_mpoly_t Q, const FIELD##_mpoly_t A,          \
    const FIELD##_mpoly_t B, const FIELD##_mpoly_ctx_t ctx)                    \
{                                                                               \
    if (B->length == 0) return A->length == 0;                                  \
    if (A->length == 0) { FIELD##_mpoly_zero(Q,ctx); return 1; }               \
    return FIELD##_mpoly_divides_dense(Q, A, B, ctx);                          \
}

#define GF216_EQ_CHECK(a,b)   ((a)!=(b))
#define GF232_EQ_CHECK(a,b)   (!gf232_equal(&(a),&(b)))
#define GF264_EQ_CHECK(a,b)   (!gf264_equal(&(a),&(b)))
#define GF2128_EQ_CHECK(a,b)  (!gf2128_equal(&(a),&(b)))

DEFINE_BASIC_OPS(gf216,  uint16_t, GF216_EQ_CHECK)
DEFINE_BASIC_OPS(gf232,  gf232_t,  GF232_EQ_CHECK)
DEFINE_BASIC_OPS(gf264,  gf264_t,  GF264_EQ_CHECK)
DEFINE_BASIC_OPS(gf2128, gf2128_t, GF2128_EQ_CHECK)

/* ============================================================================
   FIT LENGTH AND RESET BITS
   ============================================================================ */

#define DEFINE_FIT_LENGTH_RESET_BITS(FIELD, COEFF_T)                           \
void FIELD##_mpoly_fit_length_reset_bits(FIELD##_mpoly_t poly, slong len,      \
    flint_bitcnt_t bits, const FIELD##_mpoly_ctx_t ctx)                        \
{                                                                               \
    slong N = mpoly_words_per_exp(bits, ctx->minfo);                           \
    if (len > poly->coeffs_alloc) {                                             \
        slong na = FLINT_MAX(len, 2 * poly->coeffs_alloc);                     \
        poly->coeffs = (COEFF_T *)flint_realloc(poly->coeffs,na*sizeof(COEFF_T));\
        poly->coeffs_alloc = na;                                                \
    }                                                                           \
    if (N * len > poly->exps_alloc) {                                           \
        slong na = FLINT_MAX(N * len, 2 * poly->exps_alloc);                   \
        poly->exps = (ulong *)flint_realloc(poly->exps, na * sizeof(ulong));   \
        poly->exps_alloc = na;                                                  \
    }                                                                           \
    poly->bits = bits;                                                          \
}

DEFINE_FIT_LENGTH_RESET_BITS(gf216,  uint16_t)
DEFINE_FIT_LENGTH_RESET_BITS(gf232,  gf232_t)
DEFINE_FIT_LENGTH_RESET_BITS(gf264,  gf264_t)
DEFINE_FIT_LENGTH_RESET_BITS(gf2128, gf2128_t)

/* ============================================================================
   COEFFICIENT SETTING
   ============================================================================ */

void gf216_mpoly_set_coeff_ui_ui(gf216_mpoly_t poly, uint16_t c,
    const ulong *exp, const gf216_mpoly_ctx_t ctx)
{
    if (poly->bits == 0) poly->bits = MPOLY_MIN_BITS;
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    ulong *pe = (ulong *)flint_malloc(N * sizeof(ulong));
    mpoly_set_monomial_ui(pe, exp, poly->bits, ctx->minfo);
    for (slong pos = 0; pos < poly->length; pos++) {
        if (mpoly_monomial_equal(poly->exps + N*pos, pe, N)) {
            if (c == 0) {
                for (slong i = pos; i < poly->length - 1; i++) {
                    poly->coeffs[i] = poly->coeffs[i+1];
                    mpoly_monomial_set(poly->exps+N*i, poly->exps+N*(i+1), N);
                }
                poly->length--;
            } else { poly->coeffs[pos] = c; }
            flint_free(pe); return;
        }
    }
    if (c != 0) {
        if (poly->length >= poly->coeffs_alloc) {
            slong na = FLINT_MAX(poly->length+1, 2*poly->coeffs_alloc);
            poly->coeffs = (uint16_t *)flint_realloc(poly->coeffs, na*sizeof(uint16_t));
            poly->coeffs_alloc = na;
        }
        if (N*poly->length >= poly->exps_alloc) {
            slong na = FLINT_MAX(N*(poly->length+1), 2*poly->exps_alloc);
            poly->exps = (ulong *)flint_realloc(poly->exps, na*sizeof(ulong));
            poly->exps_alloc = na;
        }
        poly->coeffs[poly->length] = c;
        mpoly_monomial_set(poly->exps + N*poly->length, pe, N);
        poly->length++;
    }
    flint_free(pe);
}

#define DEFINE_SET_COEFF_STRUCT(FIELD, COEFF_T)                                \
void FIELD##_mpoly_set_coeff_ui_ui(FIELD##_mpoly_t poly, const COEFF_T *c,    \
    const ulong *exp, const FIELD##_mpoly_ctx_t ctx)                           \
{                                                                               \
    if (poly->bits == 0) poly->bits = MPOLY_MIN_BITS;                          \
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);                     \
    ulong *pe = (ulong *)flint_malloc(N * sizeof(ulong));                      \
    mpoly_set_monomial_ui(pe, exp, poly->bits, ctx->minfo);                    \
    for (slong pos = 0; pos < poly->length; pos++) {                            \
        if (mpoly_monomial_equal(poly->exps + N*pos, pe, N)) {                 \
            if (FIELD##_is_zero(c)) {                                           \
                for (slong i = pos; i < poly->length-1; i++) {                  \
                    poly->coeffs[i] = poly->coeffs[i+1];                        \
                    mpoly_monomial_set(poly->exps+N*i, poly->exps+N*(i+1), N);  \
                }                                                               \
                poly->length--;                                                  \
            } else { poly->coeffs[pos] = *c; }                                  \
            flint_free(pe); return;                                              \
        }                                                                       \
    }                                                                           \
    if (!FIELD##_is_zero(c)) {                                                  \
        if (poly->length >= poly->coeffs_alloc) {                               \
            slong na = FLINT_MAX(poly->length+1, 2*poly->coeffs_alloc);         \
            poly->coeffs = (COEFF_T *)flint_realloc(poly->coeffs,na*sizeof(COEFF_T));\
            poly->coeffs_alloc = na;                                             \
        }                                                                       \
        if (N*poly->length >= poly->exps_alloc) {                               \
            slong na = FLINT_MAX(N*(poly->length+1), 2*poly->exps_alloc);       \
            poly->exps = (ulong *)flint_realloc(poly->exps, na*sizeof(ulong));  \
            poly->exps_alloc = na;                                               \
        }                                                                       \
        poly->coeffs[poly->length] = *c;                                         \
        mpoly_monomial_set(poly->exps + N*poly->length, pe, N);                 \
        poly->length++;                                                          \
    }                                                                           \
    flint_free(pe);                                                              \
}

DEFINE_SET_COEFF_STRUCT(gf232,  gf232_t)
DEFINE_SET_COEFF_STRUCT(gf264,  gf264_t)
DEFINE_SET_COEFF_STRUCT(gf2128, gf2128_t)

/* ============================================================================
   PRINTING
   ============================================================================ */

void gf216_mpoly_print(const gf216_mpoly_t poly, const char **vars,
    const gf216_mpoly_ctx_t ctx)
{
    if (poly->length == 0) { printf("0"); return; }
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
    for (slong i = 0; i < poly->length; i++) {
        if (i > 0) printf(" + ");
        printf("0x%04x", poly->coeffs[i]);
        mpoly_get_monomial_ui(exp, poly->exps + N*i, poly->bits, ctx->minfo);
        for (slong j = 0; j < nvars; j++)
            if (exp[j] > 0) { printf("*%s", vars[j]); if (exp[j]>1) printf("^%lu",exp[j]); }
    }
    flint_free(exp);
}

#define DEFINE_PRINT_STRUCT(FIELD, PRINT_FUNC)                                  \
void FIELD##_mpoly_print(const FIELD##_mpoly_t poly, const char **vars,        \
    const FIELD##_mpoly_ctx_t ctx)                                              \
{                                                                               \
    if (poly->length == 0) { printf("0"); return; }                             \
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);                     \
    slong nvars = ctx->minfo->nvars;                                            \
    ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));                  \
    for (slong i = 0; i < poly->length; i++) {                                  \
        if (i > 0) printf(" + ");                                               \
        PRINT_FUNC(&poly->coeffs[i]);                                            \
        mpoly_get_monomial_ui(exp, poly->exps + N*i, poly->bits, ctx->minfo);  \
        for (slong j = 0; j < nvars; j++)                                       \
            if (exp[j] > 0) { printf("*%s",vars[j]); if(exp[j]>1) printf("^%lu",exp[j]); } \
    }                                                                           \
    flint_free(exp);                                                             \
}
DEFINE_PRINT_STRUCT(gf232,  gf232_print)
DEFINE_PRINT_STRUCT(gf264,  gf264_print)
DEFINE_PRINT_STRUCT(gf2128, gf2128_print)

/* ============================================================================
   RANDOM TESTING
   ============================================================================ */

void gf28_mpoly_randtest(gf28_mpoly_t poly, flint_rand_t state,
    slong length, slong exp_bound, const gf28_mpoly_ctx_t ctx)
{
    slong nvars = ctx->minfo->nvars;
    ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
    gf28_mpoly_zero(poly, ctx);
    for (slong i = 0; i < length; i++) {
        for (slong j = 0; j < nvars; j++) exp[j] = n_randint(state, exp_bound);
        uint8_t c = n_randint(state, 255) + 1;
        gf28_mpoly_set_coeff_ui_ui(poly, c, exp, ctx);
    }
    flint_free(exp);
}

void gf216_mpoly_randtest(gf216_mpoly_t poly, flint_rand_t state,
    slong length, slong exp_bound, const gf216_mpoly_ctx_t ctx)
{
    slong nvars = ctx->minfo->nvars;
    ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
    gf216_mpoly_zero(poly, ctx);
    for (slong i = 0; i < length; i++) {
        for (slong j = 0; j < nvars; j++) exp[j] = n_randint(state, exp_bound);
        uint16_t c = n_randint(state, 65535) + 1;
        gf216_mpoly_set_coeff_ui_ui(poly, c, exp, ctx);
    }
    flint_free(exp);
}

#define DEFINE_RANDTEST_STRUCT(FIELD, COEFF_T, CREATE_RANDOM)                  \
void FIELD##_mpoly_randtest(FIELD##_mpoly_t poly, flint_rand_t state,          \
    slong length, slong exp_bound, const FIELD##_mpoly_ctx_t ctx)              \
{                                                                               \
    slong nvars = ctx->minfo->nvars;                                            \
    ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));                  \
    FIELD##_mpoly_zero(poly, ctx);                                              \
    for (slong i = 0; i < length; i++) {                                        \
        for (slong j = 0; j < nvars; j++) exp[j] = n_randint(state, exp_bound);\
        COEFF_T c = CREATE_RANDOM;                                              \
        if (!FIELD##_is_zero(&c)) FIELD##_mpoly_set_coeff_ui_ui(poly,&c,exp,ctx);\
    }                                                                           \
    flint_free(exp);                                                             \
}
DEFINE_RANDTEST_STRUCT(gf232,  gf232_t,  gf232_create(n_randtest(state)))
DEFINE_RANDTEST_STRUCT(gf264,  gf264_t,  gf264_create(n_randtest(state)))

void gf2128_mpoly_randtest(gf2128_mpoly_t poly, flint_rand_t state,
    slong length, slong exp_bound, const gf2128_mpoly_ctx_t ctx)
{
    slong nvars = ctx->minfo->nvars;
    ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
    gf2128_mpoly_zero(poly, ctx);
    for (slong i = 0; i < length; i++) {
        for (slong j = 0; j < nvars; j++) exp[j] = n_randint(state, exp_bound);
        gf2128_t c; c.low = n_randtest(state); c.high = n_randtest(state);
        if (!gf2128_is_zero(&c)) gf2128_mpoly_set_coeff_ui_ui(poly, &c, exp, ctx);
    }
    flint_free(exp);
}

/* ============================================================================
   CONVERSION FUNCTIONS
   ============================================================================ */

void fq_nmod_mpoly_to_gf216_mpoly(gf216_mpoly_t res, const fq_nmod_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx)
{
    gf216_mpoly_ctx_t ctx;
    gf216_mpoly_ctx_init(ctx, fq_mpoly_ctx->minfo->nvars, fq_mpoly_ctx->minfo->ord);
    gf216_mpoly_zero(res, ctx);
    slong len = fq_nmod_mpoly_length(poly, fq_mpoly_ctx);
    if (len == 0) { gf216_mpoly_ctx_clear(ctx); return; }
    if (!g_gf216_conversion || !g_gf216_conversion->initialized) init_gf216_conversion(fqctx);
    flint_bitcnt_t bits = FLINT_MAX(poly->bits, MPOLY_MIN_BITS);
    gf216_mpoly_fit_length_reset_bits(res, len, bits, ctx);
    res->length = 0;
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong actual = 0;
    for (slong i = 0; i < len; i++) {
        fq_nmod_t coeff; fq_nmod_init(coeff, fqctx);
        fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, poly, i, fq_mpoly_ctx);
        uint16_t c = fq_nmod_to_gf216_elem(coeff, fqctx);
        if (c != 0) {
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
            fq_nmod_mpoly_get_term_exp_ui(exp, poly, i, fq_mpoly_ctx);
            mpoly_set_monomial_ui(res->exps + N*actual, exp, bits, ctx->minfo);
            res->coeffs[actual++] = c;
            flint_free(exp);
        }
        fq_nmod_clear(coeff, fqctx);
    }
    res->length = actual;
    gf216_mpoly_ctx_clear(ctx);
}

void gf216_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf216_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx)
{
    if (!g_gf216_conversion || !g_gf216_conversion->initialized) init_gf216_conversion(fqctx);
    fq_nmod_mpoly_zero(res, fq_mpoly_ctx);
    if (poly->length == 0) return;
    slong N = mpoly_words_per_exp(poly->bits, fq_mpoly_ctx->minfo);
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    for (slong i = 0; i < poly->length; i++) {
        if (poly->coeffs[i] != 0) {
            fq_nmod_t coeff; fq_nmod_init(coeff, fqctx);
            gf216_elem_to_fq_nmod(coeff, poly->coeffs[i], fqctx);
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
            mpoly_get_monomial_ui(exp, poly->exps+N*i, poly->bits, fq_mpoly_ctx->minfo);
            fq_nmod_mpoly_set_coeff_fq_nmod_ui(res, coeff, exp, fq_mpoly_ctx);
            flint_free(exp); fq_nmod_clear(coeff, fqctx);
        }
    }
}

void fq_nmod_mpoly_to_gf232_mpoly(gf232_mpoly_t res, const fq_nmod_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx)
{
    gf232_mpoly_ctx_t ctx;
    gf232_mpoly_ctx_init(ctx, fq_mpoly_ctx->minfo->nvars, fq_mpoly_ctx->minfo->ord);
    gf232_mpoly_zero(res, ctx);
    slong len = fq_nmod_mpoly_length(poly, fq_mpoly_ctx);
    if (len == 0) { gf232_mpoly_ctx_clear(ctx); return; }
    if (!g_gf232_conversion || !g_gf232_conversion->initialized) init_gf232_conversion(fqctx);
    flint_bitcnt_t bits = FLINT_MAX(poly->bits, MPOLY_MIN_BITS);
    gf232_mpoly_fit_length_reset_bits(res, len, bits, ctx);
    res->length = 0;
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong actual = 0;
    for (slong i = 0; i < len; i++) {
        fq_nmod_t coeff; fq_nmod_init(coeff, fqctx);
        fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, poly, i, fq_mpoly_ctx);
        gf232_t c = fq_nmod_to_gf232(coeff, fqctx);
        if (!gf232_is_zero(&c)) {
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
            fq_nmod_mpoly_get_term_exp_ui(exp, poly, i, fq_mpoly_ctx);
            mpoly_set_monomial_ui(res->exps + N*actual, exp, bits, ctx->minfo);
            res->coeffs[actual++] = c;
            flint_free(exp);
        }
        fq_nmod_clear(coeff, fqctx);
    }
    res->length = actual;
    gf232_mpoly_ctx_clear(ctx);
}

void gf232_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf232_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx)
{
    if (!g_gf232_conversion || !g_gf232_conversion->initialized) init_gf232_conversion(fqctx);
    fq_nmod_mpoly_zero(res, fq_mpoly_ctx);
    if (poly->length == 0) return;
    slong N = mpoly_words_per_exp(poly->bits, fq_mpoly_ctx->minfo);
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    for (slong i = 0; i < poly->length; i++) {
        if (!gf232_is_zero(&poly->coeffs[i])) {
            fq_nmod_t coeff; fq_nmod_init(coeff, fqctx);
            gf232_to_fq_nmod(coeff, &poly->coeffs[i], fqctx);
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
            mpoly_get_monomial_ui(exp, poly->exps+N*i, poly->bits, fq_mpoly_ctx->minfo);
            fq_nmod_mpoly_set_coeff_fq_nmod_ui(res, coeff, exp, fq_mpoly_ctx);
            flint_free(exp); fq_nmod_clear(coeff, fqctx);
        }
    }
}

void fq_nmod_mpoly_to_gf264_mpoly(gf264_mpoly_t res, const fq_nmod_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx)
{
    gf264_mpoly_ctx_t ctx;
    gf264_mpoly_ctx_init(ctx, fq_mpoly_ctx->minfo->nvars, fq_mpoly_ctx->minfo->ord);
    gf264_mpoly_zero(res, ctx);
    slong len = fq_nmod_mpoly_length(poly, fq_mpoly_ctx);
    if (len == 0) { gf264_mpoly_ctx_clear(ctx); return; }
    if (!g_gf264_conversion || !g_gf264_conversion->initialized) init_gf264_conversion(fqctx);
    flint_bitcnt_t bits = FLINT_MAX(poly->bits, MPOLY_MIN_BITS);
    gf264_mpoly_fit_length_reset_bits(res, len, bits, ctx);
    res->length = 0;
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong actual = 0;
    for (slong i = 0; i < len; i++) {
        fq_nmod_t coeff; fq_nmod_init(coeff, fqctx);
        fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, poly, i, fq_mpoly_ctx);
        gf264_t c = fq_nmod_to_gf264(coeff, fqctx);
        if (!gf264_is_zero(&c)) {
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
            fq_nmod_mpoly_get_term_exp_ui(exp, poly, i, fq_mpoly_ctx);
            mpoly_set_monomial_ui(res->exps + N*actual, exp, bits, ctx->minfo);
            res->coeffs[actual++] = c;
            flint_free(exp);
        }
        fq_nmod_clear(coeff, fqctx);
    }
    res->length = actual;
    gf264_mpoly_ctx_clear(ctx);
}

void gf264_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf264_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx)
{
    if (!g_gf264_conversion || !g_gf264_conversion->initialized) init_gf264_conversion(fqctx);
    fq_nmod_mpoly_zero(res, fq_mpoly_ctx);
    if (poly->length == 0) return;
    slong N = mpoly_words_per_exp(poly->bits, fq_mpoly_ctx->minfo);
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    for (slong i = 0; i < poly->length; i++) {
        if (!gf264_is_zero(&poly->coeffs[i])) {
            fq_nmod_t coeff; fq_nmod_init(coeff, fqctx);
            gf264_to_fq_nmod(coeff, &poly->coeffs[i], fqctx);
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
            mpoly_get_monomial_ui(exp, poly->exps+N*i, poly->bits, fq_mpoly_ctx->minfo);
            fq_nmod_mpoly_set_coeff_fq_nmod_ui(res, coeff, exp, fq_mpoly_ctx);
            flint_free(exp); fq_nmod_clear(coeff, fqctx);
        }
    }
}

void fq_nmod_mpoly_to_gf2128_mpoly(gf2128_mpoly_t res, const fq_nmod_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx)
{
    gf2128_mpoly_ctx_t ctx;
    gf2128_mpoly_ctx_init(ctx, fq_mpoly_ctx->minfo->nvars, fq_mpoly_ctx->minfo->ord);
    gf2128_mpoly_zero(res, ctx);
    slong len = fq_nmod_mpoly_length(poly, fq_mpoly_ctx);
    if (len == 0) { gf2128_mpoly_ctx_clear(ctx); return; }
    if (!g_gf2128_conversion || !g_gf2128_conversion->initialized) init_gf2128_conversion(fqctx);
    flint_bitcnt_t bits = FLINT_MAX(poly->bits, MPOLY_MIN_BITS);
    gf2128_mpoly_fit_length_reset_bits(res, len, bits, ctx);
    res->length = 0;
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong actual = 0;
    for (slong i = 0; i < len; i++) {
        fq_nmod_t coeff; fq_nmod_init(coeff, fqctx);
        fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, poly, i, fq_mpoly_ctx);
        gf2128_t c = fq_nmod_to_gf2128(coeff, fqctx);
        if (!gf2128_is_zero(&c)) {
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
            fq_nmod_mpoly_get_term_exp_ui(exp, poly, i, fq_mpoly_ctx);
            mpoly_set_monomial_ui(res->exps + N*actual, exp, bits, ctx->minfo);
            res->coeffs[actual++] = c;
            flint_free(exp);
        }
        fq_nmod_clear(coeff, fqctx);
    }
    res->length = actual;
    gf2128_mpoly_ctx_clear(ctx);
}

void gf2128_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf2128_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx)
{
    if (!g_gf2128_conversion || !g_gf2128_conversion->initialized) init_gf2128_conversion(fqctx);
    fq_nmod_mpoly_zero(res, fq_mpoly_ctx);
    if (poly->length == 0) return;
    slong N = mpoly_words_per_exp(poly->bits, fq_mpoly_ctx->minfo);
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    for (slong i = 0; i < poly->length; i++) {
        if (!gf2128_is_zero(&poly->coeffs[i])) {
            fq_nmod_t coeff; fq_nmod_init(coeff, fqctx);
            gf2128_to_fq_nmod(coeff, &poly->coeffs[i], fqctx);
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
            mpoly_get_monomial_ui(exp, poly->exps+N*i, poly->bits, fq_mpoly_ctx->minfo);
            fq_nmod_mpoly_set_coeff_fq_nmod_ui(res, coeff, exp, fq_mpoly_ctx);
            flint_free(exp); fq_nmod_clear(coeff, fqctx);
        }
    }
}

/* ============================================================================
   GF(2^8) — Optimized implementation (supports N = 1 and N = 2)
   ============================================================================ */

void gf28_mpoly_init(gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx) {
    poly->coeffs = NULL; poly->exps = NULL; poly->length = 0;
    poly->coeffs_alloc = 0; poly->exps_alloc = 0; poly->bits = MPOLY_MIN_BITS;
}
void gf28_mpoly_clear(gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx) {
    if (poly->coeffs) flint_free(poly->coeffs);
    if (poly->exps)   flint_free(poly->exps);
}
void gf28_mpoly_ctx_init(gf28_mpoly_ctx_t ctx, slong nvars, const ordering_t ord) {
    mpoly_ctx_init(ctx->minfo, nvars, ord);
}
void gf28_mpoly_ctx_clear(gf28_mpoly_ctx_t ctx) { mpoly_ctx_clear(ctx->minfo); }
void gf28_mpoly_zero(gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx) { poly->length = 0; }

void _gf28_mpoly_fit_length(uint8_t **coeffs, slong *coeffs_alloc,
                             ulong **exps, slong *exps_alloc, slong N, slong length) {
    if (length > *coeffs_alloc) {
        slong na = FLINT_MAX(length, 2 * (*coeffs_alloc));
        *coeffs = (uint8_t *)flint_realloc(*coeffs, na * sizeof(uint8_t));
        *coeffs_alloc = na;
    }
    if (N * length > *exps_alloc) {
        slong na = FLINT_MAX(N * length, 2 * (*exps_alloc));
        *exps = (ulong *)flint_realloc(*exps, na * sizeof(ulong));
        *exps_alloc = na;
    }
}

void gf28_mpoly_fit_length_reset_bits(gf28_mpoly_t poly, slong len,
                                       flint_bitcnt_t bits, const gf28_mpoly_ctx_t ctx) {
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    _gf28_mpoly_fit_length(&poly->coeffs, &poly->coeffs_alloc,
                            &poly->exps,   &poly->exps_alloc, N, len);
    poly->bits = bits;
}

void gf28_mpoly_init3(gf28_mpoly_t poly, slong alloc, flint_bitcnt_t bits,
                       const gf28_mpoly_ctx_t ctx) {
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    poly->coeffs = (uint8_t *)flint_malloc(alloc * sizeof(uint8_t));
    poly->exps   = (ulong *)flint_malloc(N * alloc * sizeof(ulong));
    poly->coeffs_alloc = alloc;
    poly->exps_alloc   = N * alloc;
    poly->length = 0; poly->bits = bits;
}

void gf28_mpoly_swap(gf28_mpoly_t p1, gf28_mpoly_t p2, const gf28_mpoly_ctx_t ctx) {
    gf28_mpoly_struct t = *p1; *p1 = *p2; *p2 = t;
}
void _gf28_mpoly_set_length(gf28_mpoly_t poly, slong len, const gf28_mpoly_ctx_t ctx) {
    poly->length = len;
}

void gf28_mpoly_set(gf28_mpoly_t res, const gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx) {
    if (res == poly) return;
    gf28_mpoly_fit_length_reset_bits(res, poly->length, poly->bits, ctx);
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    memcpy(res->coeffs, poly->coeffs, poly->length * sizeof(uint8_t));
    memcpy(res->exps,   poly->exps,   N * poly->length * sizeof(ulong));
    res->length = poly->length;
}

static void _gf28_mpoly_repack_exps(
    gf28_mpoly_t poly,
    flint_bitcnt_t new_bits,
    const gf28_mpoly_ctx_t ctx)
{
    if (poly->bits == new_bits) return;
    if (poly->length == 0) { poly->bits = new_bits; return; }

    slong nvars  = ctx->minfo->nvars;
    slong N_old  = mpoly_words_per_exp(poly->bits, ctx->minfo);
    slong N_new  = mpoly_words_per_exp(new_bits,   ctx->minfo);
    ulong *new_exps = (ulong *)flint_malloc(N_new * poly->length * sizeof(ulong));
    ulong *expv     = (ulong *)flint_malloc(nvars  * sizeof(ulong));

    for (slong i = 0; i < poly->length; i++) {
        mpoly_get_monomial_ui(expv,
            poly->exps + N_old * i, poly->bits, ctx->minfo);
        mpoly_set_monomial_ui(new_exps + N_new * i,
            expv, new_bits, ctx->minfo);
    }

    flint_free(expv);
    flint_free(poly->exps);
    poly->exps       = new_exps;
    poly->exps_alloc = N_new * poly->length;
    poly->bits       = new_bits;
}

void gf28_mpoly_set_coeff_ui_ui(gf28_mpoly_t poly, uint8_t c,
                                 const ulong *exp, const gf28_mpoly_ctx_t ctx) {
    if (poly->bits == 0) poly->bits = MPOLY_MIN_BITS;
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    ulong *cmpmask  = (ulong *)flint_malloc(N * FLINT_BITS * sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, poly->bits, ctx->minfo);
    ulong *pe = (ulong *)flint_malloc(N * sizeof(ulong));
    mpoly_set_monomial_ui(pe, exp, poly->bits, ctx->minfo);
    slong pos;
    for (pos = 0; pos < poly->length; pos++) {
        if (mpoly_monomial_equal(poly->exps + N*pos, pe, N)) {
            if (c == 0) {
                for (slong i = pos; i < poly->length-1; i++) {
                    poly->coeffs[i] = poly->coeffs[i+1];
                    mpoly_monomial_set(poly->exps+N*i, poly->exps+N*(i+1), N);
                }
                poly->length--;
            } else { poly->coeffs[pos] = c; }
            flint_free(cmpmask); flint_free(pe); return;
        }
        if (mpoly_monomial_lt(poly->exps + N*pos, pe, N, cmpmask)) break;
    }
    if (c != 0) {
        _gf28_mpoly_fit_length(&poly->coeffs, &poly->coeffs_alloc,
                                &poly->exps,   &poly->exps_alloc, N, poly->length+1);
        for (slong i = poly->length; i > pos; i--) {
            poly->coeffs[i] = poly->coeffs[i-1];
            mpoly_monomial_set(poly->exps+N*i, poly->exps+N*(i-1), N);
        }
        poly->coeffs[pos] = c;
        mpoly_monomial_set(poly->exps+N*pos, pe, N);
        poly->length++;
    }
    flint_free(cmpmask); flint_free(pe);
}

void gf28_mpoly_print(const gf28_mpoly_t poly, const char **vars,
                       const gf28_mpoly_ctx_t ctx) {
    if (poly->length == 0) { printf("0"); return; }
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    for (slong i = 0; i < poly->length; i++) {
        if (i > 0) printf(" + ");
        printf("0x%02x", poly->coeffs[i]);
        ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
        mpoly_get_monomial_ui(exp, poly->exps + N*i, poly->bits, ctx->minfo);
        int fv = 1;
        for (slong j = 0; j < nvars; j++) {
            if (exp[j] > 0) {
                printf("%s%s", fv?"*":"*", vars[j]);
                if (exp[j] > 1) printf("^%lu", exp[j]);
                fv = 0;
            }
        }
        flint_free(exp);
    }
}

/* ============================================================================
   GF(2^8) ARRAY MULTIPLICATION (supports N = 1 and N = 2)
   ============================================================================ */

void _gf28_mpoly_addmul_array1_safe(uint8_t *poly1, slong array_size,
    const uint8_t *poly2, const ulong *exp2, slong len2,
    const uint8_t *poly3, const ulong *exp3, slong len3)
{
    for (slong ii = 0; ii < len2 + BLOCK; ii += BLOCK) {
        for (slong jj = 0; jj < len3 + BLOCK; jj += BLOCK) {
            for (slong i = ii; i < FLINT_MIN(ii + BLOCK, len2); i++) {
                slong off2 = (slong)exp2[i];
                if (off2 >= array_size) continue;
                uint8_t *c2 = poly1 + off2;
                if (poly2[i] != 0) {
                    const uint8_t *mul_row = gf28_get_scalar_row(poly2[i]);
                    for (slong j = jj; j < FLINT_MIN(jj + BLOCK, len3); j++) {
                        slong off3 = (slong)exp3[j];
                        if (off2 + off3 >= array_size) continue;
                        c2[off3] ^= mul_row[poly3[j]];
                    }
                }
            }
        }
    }
}

slong gf28_mpoly_append_array_LEX_safe(gf28_mpoly_t P, slong Plen,
    uint8_t *coeff_array, const ulong *mults, slong num, slong array_size,
    slong top, const gf28_mpoly_ctx_t ctx)
{
    slong N     = mpoly_words_per_exp(P->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong *expv = (ulong *)flint_malloc(nvars * sizeof(ulong));
    ulong *pe   = (ulong *)flint_malloc(N     * sizeof(ulong));

    for (slong off = array_size - 1; off >= 0; off--) {
        if (coeff_array[off] == 0) continue;

        uint8_t coeff    = coeff_array[off];
        coeff_array[off] = 0; 
        slong d = off;
        for (slong j = 0; j < num; j++) {
            expv[nvars - 1 - j] = (ulong)(d % mults[j]);
            d /= mults[j];
        }
        expv[0] = (ulong)top; 

        mpoly_set_monomial_ui(pe, expv, P->bits, ctx->minfo);

        _gf28_mpoly_fit_length(&P->coeffs, &P->coeffs_alloc,
                                &P->exps,   &P->exps_alloc, N, Plen + 1);
        memcpy(P->exps + N * Plen, pe, N * sizeof(ulong));
        P->coeffs[Plen] = coeff;
        Plen++;
    }

    flint_free(expv);
    flint_free(pe);
    return Plen;
}

void gf28_dynamic_array_init(gf28_dynamic_array_t *arr, slong size) {
    arr->size = size;
    arr->data = (uint8_t *)calloc(size, sizeof(uint8_t));
    if (!arr->data && size > 0) { flint_printf("OOM: %ld bytes\n", size); flint_abort(); }
}
void gf28_dynamic_array_clear(gf28_dynamic_array_t *arr) {
    if (arr->data) { free(arr->data); arr->data = NULL; }
}

void _gf28_mpoly_mul_array_chunked_LEX(gf28_mpoly_t P,
    const gf28_mpoly_t A, const gf28_mpoly_t B,
    const ulong *mults, const gf28_mpoly_ctx_t ctx)
{
    slong num   = ctx->minfo->nfields - 1;
    slong nvars = ctx->minfo->nvars;
    slong NA    = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong NB    = mpoly_words_per_exp(B->bits, ctx->minfo);

    slong array_size = 1;
    for (slong i = 0; i < num; i++) array_size *= mults[i];

    ulong *tmpexp = (ulong *)flint_malloc(nvars * sizeof(ulong));
    mpoly_get_monomial_ui(tmpexp, A->exps, A->bits, ctx->minfo);
    slong Al = 1 + (slong)tmpexp[0];
    mpoly_get_monomial_ui(tmpexp, B->exps, B->bits, ctx->minfo);
    slong Bl = 1 + (slong)tmpexp[0];
    flint_free(tmpexp);

    TMP_INIT; TMP_START;
    slong *Amain = (slong *)TMP_ALLOC((Al + 1) * sizeof(slong));
    slong *Bmain = (slong *)TMP_ALLOC((Bl + 1) * sizeof(slong));
    ulong *Apexp = (ulong *)flint_malloc(A->length * sizeof(ulong));
    ulong *Bpexp = (ulong *)flint_malloc(B->length * sizeof(ulong));

    _mpoly_main_variable_split_LEX_N(Amain, Apexp, A->exps, Al,
        A->length, mults, num, A->bits, NA, ctx->minfo);
    _mpoly_main_variable_split_LEX_N(Bmain, Bpexp, B->exps, Bl,
        B->length, mults, num, B->bits, NB, ctx->minfo);

    slong Pl = Al + Bl - 1, Plen = 0;

    gf28_dynamic_array_t ca;
    gf28_dynamic_array_init(&ca, array_size);

    for (slong Pi = 0; Pi < Pl; Pi++) {
        memset(ca.data, 0, array_size);
        for (slong i = 0, j = Pi; i < Al && j >= 0; i++, j--) {
            if (j < Bl)
                _gf28_mpoly_addmul_array1_safe(ca.data, array_size,
                    A->coeffs + Amain[i], Apexp + Amain[i], Amain[i+1] - Amain[i],
                    B->coeffs + Bmain[j], Bpexp + Bmain[j], Bmain[j+1] - Bmain[j]);
        }
        Plen = gf28_mpoly_append_array_LEX_safe(P, Plen, ca.data,
            mults, num, array_size, Pl - Pi - 1, ctx);
    }
    _gf28_mpoly_set_length(P, Plen, ctx);
    gf28_dynamic_array_clear(&ca);
    flint_free(Apexp); flint_free(Bpexp);
    TMP_END;
}

int _gf28_mpoly_mul_array_LEX(gf28_mpoly_t A,
    const gf28_mpoly_t B, fmpz *maxBfields,
    const gf28_mpoly_t C, fmpz *maxCfields,
    const gf28_mpoly_ctx_t ctx)
{
    slong exp_bits;
    ulong array_size, max, *mults;
    int success;
    TMP_INIT; TMP_START;

    FLINT_ASSERT(ctx->minfo->nvars > 0);
    FLINT_ASSERT(B->length != 0 && C->length != 0);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    mults = (ulong *)TMP_ALLOC(ctx->minfo->nfields * sizeof(ulong));

    slong i = ctx->minfo->nfields - 1;
    mults[i] = 1 + fmpz_get_ui(maxBfields + i) + fmpz_get_ui(maxCfields + i);
    max = mults[i];
    if (((slong)mults[i]) <= 0) { success = 0; goto cleanup; }

    array_size = WORD(1);
    for (i--; i >= 0; i--) {
        ulong hi;
        mults[i] = 1 + fmpz_get_ui(maxBfields + i) + fmpz_get_ui(maxCfields + i);
        max |= mults[i];
        umul_ppmm(hi, array_size, array_size, mults[i]);
        if (hi != 0 || (slong)mults[i] <= 0 || array_size <= 0)
            { success = 0; goto cleanup; }
    }
    if (array_size > (1L << 28)) { success = 0; goto cleanup; }

    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, FLINT_BIT_COUNT(max) + 1);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    if (mpoly_words_per_exp(exp_bits, ctx->minfo) > 2)
        { success = 0; goto cleanup; }

    {
        gf28_mpoly_t Btmp, Ctmp;
        gf28_mpoly_init(Btmp, ctx);
        gf28_mpoly_init(Ctmp, ctx);
        gf28_mpoly_set(Btmp, B, ctx);
        gf28_mpoly_set(Ctmp, C, ctx);
        _gf28_mpoly_repack_exps(Btmp, exp_bits, ctx); 
        _gf28_mpoly_repack_exps(Ctmp, exp_bits, ctx);

        if (A == B || A == C) {
            gf28_mpoly_t T;
            gf28_mpoly_init3(T, B->length + C->length - 1, exp_bits, ctx);
            _gf28_mpoly_mul_array_chunked_LEX(T, Ctmp, Btmp, mults, ctx);
            gf28_mpoly_swap(T, A, ctx);
            gf28_mpoly_clear(T, ctx);
        } else {
            gf28_mpoly_fit_length_reset_bits(A, B->length + C->length - 1,
                                              exp_bits, ctx);
            _gf28_mpoly_mul_array_chunked_LEX(A, Ctmp, Btmp, mults, ctx);
        }

        gf28_mpoly_clear(Btmp, ctx);
        gf28_mpoly_clear(Ctmp, ctx);
    }
    success = 1;

cleanup:
    TMP_END;
    return success;
}

int gf28_mpoly_mul_array(gf28_mpoly_t A, const gf28_mpoly_t B,
                          const gf28_mpoly_t C, const gf28_mpoly_ctx_t ctx)
{
    slong i;
    int   success;
    fmpz *maxBf, *maxCf;
    TMP_INIT;

    if (!g_gf28_complete_tables.initialized)
        init_gf28_standard();

    if (B->length == 0 || C->length == 0) { gf28_mpoly_zero(A, ctx); return 1; }
    if (B->bits == 0 || C->bits == 0)      return 0;

    if (ctx->minfo->nvars < 1 ||
        mpoly_words_per_exp(B->bits, ctx->minfo) > 2 ||
        mpoly_words_per_exp(C->bits, ctx->minfo) > 2)
        return 0;

    TMP_START;
    maxBf = (fmpz *)TMP_ALLOC(ctx->minfo->nfields * sizeof(fmpz));
    maxCf = (fmpz *)TMP_ALLOC(ctx->minfo->nfields * sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++) {
        fmpz_init(maxBf + i);
        fmpz_init(maxCf + i);
    }
    mpoly_max_fields_fmpz(maxBf, B->exps, B->length, B->bits, ctx->minfo);
    mpoly_max_fields_fmpz(maxCf, C->exps, C->length, C->bits, ctx->minfo);

    switch (ctx->minfo->ord) {
        case ORD_LEX:
            success = _gf28_mpoly_mul_array_LEX(A, B, maxBf, C, maxCf, ctx);
            break;
        default:
            success = 0;
    }

    for (i = 0; i < ctx->minfo->nfields; i++) {
        fmpz_clear(maxBf + i);
        fmpz_clear(maxCf + i);
    }
    TMP_END;
    return success;
}

double get_wall_time_8(void) {
    struct timeval time;
    if (gettimeofday(&time, NULL)) return 0;
    return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
}

/* ============================================================================
   GF(2^8) DIVISION
   ============================================================================ */

void gf28_mpolyd_init(gf28_mpolyd_t A, slong nvars) {
    A->coeffs = NULL; A->alloc = 0;
    A->deg_bounds = (slong *)calloc(nvars, sizeof(slong)); A->nvars = nvars;
}
void gf28_mpolyd_clear(gf28_mpolyd_t A) {
    if (A->coeffs)     free(A->coeffs);
    if (A->deg_bounds) free(A->deg_bounds);
}
slong gf28_mpolyd_offset(const gf28_mpolyd_t A, const ulong *exp) {
    slong off = 0, stride = 1;
    for (slong i = 0; i < A->nvars; i++) {
        if (exp[i] >= (ulong)A->deg_bounds[i]) return -1;
        off += (slong)exp[i] * stride; stride *= A->deg_bounds[i];
    }
    return off;
}
int gf28_mpolyd_set_degbounds(gf28_mpolyd_t A, const slong *bounds) {
    slong size = 1;
    for (slong i = 0; i < A->nvars; i++) {
        A->deg_bounds[i] = bounds[i];
        if (bounds[i] <= 0) return 0;
        if (size > WORD_MAX / bounds[i]) return 0;
        size *= bounds[i];
    }
    if (size > (1L << 28)) return 0;
    A->alloc = size;
    A->coeffs = (uint8_t *)calloc(size, sizeof(uint8_t));
    return A->coeffs != NULL;
}

void batch_xor_sse2(uint8_t *dst, const uint8_t *src, slong len) {
    slong i = 0;
#if DIXON_X86_SIMD
    for (; i + 16 <= len; i += 16) {
        __m128i a = _mm_loadu_si128((__m128i *)(dst + i));
        __m128i b = _mm_loadu_si128((__m128i *)(src + i));
        _mm_storeu_si128((__m128i *)(dst + i), _mm_xor_si128(a, b));
    }
#endif
    for (; i < len; i++) dst[i] ^= src[i];
}

void gf28_mpolyd_divrem_univar(gf28_mpolyd_t Q, gf28_mpolyd_t R,
    const gf28_mpolyd_t A, const gf28_mpolyd_t B)
{
    slong n = A->alloc;
    memcpy(R->coeffs, A->coeffs, n * sizeof(uint8_t));
    memset(Q->coeffs, 0, Q->alloc * sizeof(uint8_t));
    slong degB = -1;
    for (slong i = n-1; i >= 0; i--)
        if (B->coeffs[i] != 0) { degB = i; break; }
    if (degB < 0) return;
    uint8_t lc_B_inv = gf28_inv(B->coeffs[degB]);
    uint8_t (*Bm)[degB + 1] = malloc(256 * (degB + 1) * sizeof(uint8_t));
    if (!Bm) return;
    for (int k = 0; k < 256; k++) {
        const uint8_t *row = gf28_get_scalar_row(k);
        for (slong j = 0; j <= degB; j++) Bm[k][j] = row[B->coeffs[j]];
    }
    for (slong i = n-1; i >= degB; i--) {
        if (R->coeffs[i] != 0) {
            uint8_t q = gf28_mul(R->coeffs[i], lc_B_inv);
            if (i - degB < Q->alloc) Q->coeffs[i - degB] = q;
            const uint8_t *Bq = Bm[q];
            slong len = FLINT_MIN(degB + 1, n - (i - degB));
            for (slong j = 0; j < len; j++) R->coeffs[i-degB+j] ^= Bq[j];
        }
    }
    free(Bm);
}

int gf28_mpolyd_is_zero(const gf28_mpolyd_t A) {
    for (slong i = 0; i < A->alloc; i++) if (A->coeffs[i] != 0) return 0;
    return 1;
}
int gf28_mpoly_is_monomial(const gf28_mpoly_t poly)       { return poly->length == 1; }
uint8_t gf28_mpoly_get_monomial_coeff(const gf28_mpoly_t poly) {
    return poly->length != 1 ? 0 : poly->coeffs[0]; }
void gf28_mpoly_get_monomial_exp(ulong *exp, const gf28_mpoly_t poly,
    const gf28_mpoly_ctx_t ctx) {
    if (poly->length != 1) return;
    mpoly_get_monomial_ui(exp, poly->exps, poly->bits, ctx->minfo);
}

void gf28_mpoly_to_mpolyd(gf28_mpolyd_t A, const gf28_mpoly_t B,
    const gf28_mpoly_ctx_t ctx)
{
    slong nvars = ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
    ulong *exp = (ulong *)malloc(nvars * sizeof(ulong));
    memset(A->coeffs, 0, A->alloc * sizeof(uint8_t));
    for (slong i = 0; i < B->length; i++) {
        mpoly_get_monomial_ui(exp, B->exps + N*i, B->bits, ctx->minfo);
        slong off = gf28_mpolyd_offset(A, exp);
        if (off >= 0 && off < A->alloc) A->coeffs[off] = B->coeffs[i];
    }
    free(exp);
}

void gf28_mpolyd_to_mpoly_fast(gf28_mpoly_t A, const gf28_mpolyd_t B,
    const gf28_mpoly_ctx_t ctx)
{
    slong nvars = ctx->minfo->nvars;
    flint_bitcnt_t bits_needed = MPOLY_MIN_BITS;
    for (slong i = 0; i < nvars; i++)
        if (B->deg_bounds[i] > 0)
            bits_needed = FLINT_MAX(bits_needed,
                (flint_bitcnt_t)FLINT_BIT_COUNT(B->deg_bounds[i]-1));
    if (bits_needed < 16) bits_needed = 16;
    bits_needed = mpoly_fix_bits(bits_needed, ctx->minfo);
    slong N = mpoly_words_per_exp(bits_needed, ctx->minfo);
    slong nz = 0;
    for (slong off = 0; off < B->alloc; off++) if (B->coeffs[off] != 0) nz++;
    if (nz == 0) { gf28_mpoly_zero(A, ctx); return; }
    gf28_mpoly_fit_length_reset_bits(A, nz, bits_needed, ctx);
    A->length = 0;
    ulong *eb = (ulong *)malloc(nvars * sizeof(ulong));
    ulong *pb = (ulong *)malloc(N * sizeof(ulong));
    for (slong off = B->alloc - 1; off >= 0; off--) {
        if (B->coeffs[off] != 0) {
            slong temp = off;
            for (slong i = 0; i < nvars; i++) {
                eb[i] = (ulong)(temp % B->deg_bounds[i]);
                temp /= B->deg_bounds[i];
            }
            mpoly_set_monomial_ui(pb, eb, bits_needed, ctx->minfo);
            A->coeffs[A->length] = B->coeffs[off];
            mpoly_monomial_set(A->exps + N*A->length, pb, N);
            A->length++;
        }
    }
    free(eb); free(pb);
}

void gf28_mpolyd_to_mpoly_univariate(gf28_mpoly_t A, const gf28_mpolyd_t B,
    const gf28_mpoly_ctx_t ctx)
{
    if (ctx->minfo->nvars != 1) { gf28_mpolyd_to_mpoly_fast(A, B, ctx); return; }
    slong db = B->deg_bounds[0], nz = 0, hi = -1;
    for (slong i = 0; i < db; i++) if (B->coeffs[i] != 0) { nz++; hi = i; }
    if (nz == 0) { gf28_mpoly_zero(A, ctx); return; }
    flint_bitcnt_t bn = FLINT_MAX((flint_bitcnt_t)FLINT_BIT_COUNT(hi)+1, MPOLY_MIN_BITS);
    bn = mpoly_fix_bits(bn, ctx->minfo);
    gf28_mpoly_fit_length_reset_bits(A, nz, bn, ctx);
    A->length = 0;
    slong N = mpoly_words_per_exp(bn, ctx->minfo);
    for (slong deg = hi; deg >= 0; deg--) {
        if (B->coeffs[deg] != 0) {
            A->coeffs[A->length] = B->coeffs[deg];
            if (N == 1) {
                A->exps[A->length] = (ulong)deg;
            } else {
                ulong e[1] = {(ulong)deg};
                mpoly_set_monomial_ui(A->exps + A->length*N, e, bn, ctx->minfo);
            }
            A->length++;
        }
    }
}

void gf28_mpolyd_to_mpoly(gf28_mpoly_t A, const gf28_mpolyd_t B,
    const gf28_mpoly_ctx_t ctx) {
    if (ctx->minfo->nvars == 1) gf28_mpolyd_to_mpoly_univariate(A, B, ctx);
    else                        gf28_mpolyd_to_mpoly_fast(A, B, ctx);
}

int gf28_mpoly_is_univariate(const gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx,
    slong *main_var)
{
    if (poly->length == 0) { *main_var = 0; return 1; }
    slong nvars = ctx->minfo->nvars;
    if (nvars == 1) { *main_var = 0; return 1; }
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    ulong *exp = (ulong *)malloc(nvars * sizeof(ulong));
    int fv = -1;
    for (slong i = 0; i < poly->length; i++) {
        mpoly_get_monomial_ui(exp, poly->exps+N*i, poly->bits, ctx->minfo);
        for (slong j = 0; j < nvars; j++) {
            if (exp[j] > 0) {
                if (fv == -1) fv = j;
                else if (fv != j) { free(exp); return 0; }
            }
        }
    }
    free(exp); *main_var = fv >= 0 ? fv : 0; return 1;
}

slong gf28_mpoly_univariate_degree(const gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx,
    slong main_var)
{
    if (poly->length == 0) return -1;
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong *exp = (ulong *)malloc(nvars * sizeof(ulong));
    slong md = 0;
    for (slong i = 0; i < poly->length; i++) {
        mpoly_get_monomial_ui(exp, poly->exps+N*i, poly->bits, ctx->minfo);
        if ((slong)exp[main_var] > md) md = (slong)exp[main_var];
    }
    free(exp); return md;
}

void gf28_mpoly_to_gf28_poly_univar(gf28_poly_t res, const gf28_mpoly_t poly,
    const gf28_mpoly_ctx_t ctx, slong main_var)
{
    gf28_poly_zero(res);
    if (poly->length == 0) return;
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong *exp = (ulong *)malloc(nvars * sizeof(ulong));
    for (slong i = 0; i < poly->length; i++) {
        mpoly_get_monomial_ui(exp, poly->exps+N*i, poly->bits, ctx->minfo);
        gf28_poly_set_coeff(res, exp[main_var], poly->coeffs[i]);
    }
    free(exp);
}

void gf28_poly_to_gf28_mpoly_univar(gf28_mpoly_t res, const gf28_poly_t poly,
    const gf28_mpoly_ctx_t ctx, slong main_var)
{
    gf28_mpoly_zero(res, ctx);
    if (poly->length == 0) return;
    slong nz = 0;
    for (slong i = 0; i < poly->length; i++) if (poly->coeffs[i] != 0) nz++;
    if (nz == 0) return;
    slong md = poly->length - 1;
    flint_bitcnt_t rb = FLINT_MAX((flint_bitcnt_t)FLINT_BIT_COUNT(md)+1, MPOLY_MIN_BITS);
    rb = mpoly_fix_bits(rb, ctx->minfo);
    gf28_mpoly_fit_length_reset_bits(res, nz, rb, ctx);
    slong nvars = ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(rb, ctx->minfo);
    slong ti = 0;
    for (slong i = poly->length - 1; i >= 0; i--) {
        if (poly->coeffs[i] != 0) {
            res->coeffs[ti] = poly->coeffs[i];
            ulong *te = (ulong *)calloc(nvars, sizeof(ulong));
            te[main_var] = i;
            mpoly_set_monomial_ui(res->exps + N*ti, te, rb, ctx->minfo);
            free(te); ti++;
        }
    }
    res->length = nz;
}

void gf28_mpoly_mul_flint_univar(gf28_mpoly_t res, const gf28_mpoly_t a,
    const gf28_mpoly_t b, const gf28_mpoly_ctx_t ctx, slong main_var)
{
    nmod_poly_t modulus; fq_nmod_ctx_t fq_ctx;
    nmod_poly_init(modulus, 2);
    nmod_poly_set_coeff_ui(modulus, 0, 1);
    nmod_poly_set_coeff_ui(modulus, 2, 1);
    nmod_poly_set_coeff_ui(modulus, 3, 1);
    nmod_poly_set_coeff_ui(modulus, 4, 1);
    nmod_poly_set_coeff_ui(modulus, 8, 1);
    fq_nmod_ctx_init_modulus(fq_ctx, modulus, "x");
    gf28_poly_t pa, pb, pr;
    gf28_poly_init(pa); gf28_poly_init(pb); gf28_poly_init(pr);
    gf28_mpoly_to_gf28_poly_univar(pa, a, ctx, main_var);
    gf28_mpoly_to_gf28_poly_univar(pb, b, ctx, main_var);
    fq_nmod_poly_t fa, fb, fr;
    fq_nmod_poly_init(fa,fq_ctx); fq_nmod_poly_init(fb,fq_ctx); fq_nmod_poly_init(fr,fq_ctx);
    gf28_poly_to_fq_nmod_poly(fa, pa, fq_ctx);
    gf28_poly_to_fq_nmod_poly(fb, pb, fq_ctx);
    fq_nmod_poly_mul(fr, fa, fb, fq_ctx);
    fq_nmod_poly_to_gf28_poly(pr, fr, fq_ctx);
    gf28_poly_to_gf28_mpoly_univar(res, pr, ctx, main_var);
    gf28_poly_clear(pa); gf28_poly_clear(pb); gf28_poly_clear(pr);
    fq_nmod_poly_clear(fa,fq_ctx); fq_nmod_poly_clear(fb,fq_ctx); fq_nmod_poly_clear(fr,fq_ctx);
    fq_nmod_ctx_clear(fq_ctx); nmod_poly_clear(modulus);
}

int gf28_mpoly_divides_monomial(gf28_mpoly_t Q, const gf28_mpoly_t A,
    const gf28_mpoly_t B, const gf28_mpoly_ctx_t ctx)
{
    slong nvars = ctx->minfo->nvars;
    ulong *eb = (ulong *)malloc(nvars*sizeof(ulong));
    ulong *ea = (ulong *)malloc(nvars*sizeof(ulong));
    ulong *eq = (ulong *)malloc(nvars*sizeof(ulong));
    uint8_t cb = gf28_mpoly_get_monomial_coeff(B);
    if (cb == 0) { free(eb);free(ea);free(eq); return 0; }
    uint8_t cbi = gf28_inv(cb);
    gf28_mpoly_get_monomial_exp(eb, B, ctx);
    flint_bitcnt_t bits = FLINT_MAX(A->bits, B->bits);
    if (bits < 16) bits = 16;
    gf28_mpoly_fit_length_reset_bits(Q, A->length, bits, ctx);
    Q->length = 0;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    for (slong i = 0; i < A->length; i++) {
        mpoly_get_monomial_ui(ea, A->exps+N*i, A->bits, ctx->minfo);
        int ok = 1;
        for (slong j = 0; j < nvars; j++) {
            if (ea[j] < eb[j]) { ok = 0; break; }
            eq[j] = ea[j] - eb[j];
        }
        if (!ok) { free(eb);free(ea);free(eq); gf28_mpoly_zero(Q,ctx); return 0; }
        uint8_t cq = gf28_mul(A->coeffs[i], cbi);
        if (cq != 0) {
            mpoly_set_monomial_ui(Q->exps+N*Q->length, eq, bits, ctx->minfo);
            Q->coeffs[Q->length] = cq; Q->length++;
        }
    }
    free(eb);free(ea);free(eq); return 1;
}

int gf28_mpoly_divides_dense(gf28_mpoly_t Q, const gf28_mpoly_t A,
    const gf28_mpoly_t B, const gf28_mpoly_ctx_t ctx)
{
    slong nvars = ctx->minfo->nvars; int success = 0;
    slong *dA = (slong *)malloc(nvars*sizeof(slong));
    slong *dB = (slong *)malloc(nvars*sizeof(slong));
    slong *bo = (slong *)malloc(nvars*sizeof(slong));
    mpoly_degrees_si(dA, A->exps, A->length, A->bits, ctx->minfo);
    mpoly_degrees_si(dB, B->exps, B->length, B->bits, ctx->minfo);
    slong mva, mvb;
    int ua = gf28_mpoly_is_univariate(A, ctx, &mva);
    int ub = gf28_mpoly_is_univariate(B, ctx, &mvb);
    for (slong i = 0; i < nvars; i++) {
        if (dA[i] < dB[i]) { free(dA);free(dB);free(bo); gf28_mpoly_zero(Q,ctx); return 0; }
        bo[i] = dA[i] + 1;
    }
    gf28_mpolyd_t Ad, Bd, Qd, Rd;
    gf28_mpolyd_init(Ad,nvars); gf28_mpolyd_init(Bd,nvars);
    gf28_mpolyd_init(Qd,nvars); gf28_mpolyd_init(Rd,nvars);
    if (!gf28_mpolyd_set_degbounds(Ad,bo)||!gf28_mpolyd_set_degbounds(Bd,bo)||
        !gf28_mpolyd_set_degbounds(Qd,bo)||!gf28_mpolyd_set_degbounds(Rd,bo))
        goto cleanup;
    gf28_mpoly_to_mpolyd(Ad,A,ctx);
    gf28_mpoly_to_mpolyd(Bd,B,ctx);
    gf28_mpolyd_divrem_univar(Qd,Rd,Ad,Bd);
    if (gf28_mpolyd_is_zero(Rd)) {
        if (ua && ub && mva == mvb) gf28_mpolyd_to_mpoly_univariate(Q,Qd,ctx);
        else                        gf28_mpolyd_to_mpoly_fast(Q,Qd,ctx);
        success = 1;
    } else { gf28_mpoly_zero(Q,ctx); success = 0; }
cleanup:
    gf28_mpolyd_clear(Ad);gf28_mpolyd_clear(Bd);
    gf28_mpolyd_clear(Qd);gf28_mpolyd_clear(Rd);
    free(dA);free(dB);free(bo); return success;
}

int gf28_mpoly_mul(gf28_mpoly_t res, const gf28_mpoly_t a,
    const gf28_mpoly_t b, const gf28_mpoly_ctx_t ctx)
{
    if (a->length == 0 || b->length == 0) { gf28_mpoly_zero(res,ctx); return 1; }
    slong mva, mvb;
    int ua = gf28_mpoly_is_univariate(a, ctx, &mva);
    int ub = gf28_mpoly_is_univariate(b, ctx, &mvb);
    if (ua && ub && mva == mvb) {
        slong da = gf28_mpoly_univariate_degree(a,ctx,mva);
        slong db = gf28_mpoly_univariate_degree(b,ctx,mvb);
        slong mx = FLINT_MAX(da,db);
        gf28_poly_t pa,pb,pr;
        gf28_poly_init(pa); gf28_poly_init(pb); gf28_poly_init(pr);
        gf28_mpoly_to_gf28_poly_univar(pa,a,ctx,mva);
        gf28_mpoly_to_gf28_poly_univar(pb,b,ctx,mvb);
        if (mx < 10)          gf28_poly_mul_schoolbook(pr,pa,pb);
        else if (mx <= 10000) gf28_poly_mul_karatsuba(pr,pa,pb);
        else {
            gf28_poly_clear(pa);gf28_poly_clear(pb);gf28_poly_clear(pr);
            gf28_mpoly_mul_flint_univar(res,a,b,ctx,mva); return 1;
        }
        gf28_poly_to_gf28_mpoly_univar(res,pr,ctx,mva);
        gf28_poly_clear(pa);gf28_poly_clear(pb);gf28_poly_clear(pr);
        return 1;
    }
    return gf28_mpoly_mul_array(res, a, b, ctx);
}

int gf28_mpoly_divides(gf28_mpoly_t Q, const gf28_mpoly_t A,
    const gf28_mpoly_t B, const gf28_mpoly_ctx_t ctx)
{
    if (B->length == 0) return A->length == 0;
    if (A->length == 0) { gf28_mpoly_zero(Q,ctx); return 1; }
    if (gf28_mpoly_is_monomial(B)) return gf28_mpoly_divides_monomial(Q,A,B,ctx);
    return gf28_mpoly_divides_dense(Q,A,B,ctx);
}

void fq_nmod_mpoly_to_gf28_mpoly(gf28_mpoly_t res, const fq_nmod_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx)
{
    gf28_mpoly_ctx_t ctx;
    gf28_mpoly_ctx_init(ctx, fq_mpoly_ctx->minfo->nvars, fq_mpoly_ctx->minfo->ord);
    gf28_mpoly_zero(res, ctx);
    slong len = fq_nmod_mpoly_length(poly, fq_mpoly_ctx);
    if (len == 0) { gf28_mpoly_ctx_clear(ctx); return; }
    flint_bitcnt_t bits = FLINT_MAX(poly->bits, MPOLY_MIN_BITS);
    gf28_mpoly_fit_length_reset_bits(res, len, bits, ctx);
    res->length = 0;
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong actual = 0;
    for (slong i = 0; i < len; i++) {
        fq_nmod_t coeff; fq_nmod_init(coeff, fqctx);
        fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, poly, i, fq_mpoly_ctx);
        uint8_t c = fq_nmod_to_gf28_elem(coeff, fqctx);
        if (c != 0) {
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
            fq_nmod_mpoly_get_term_exp_ui(exp, poly, i, fq_mpoly_ctx);
            mpoly_set_monomial_ui(res->exps + N*actual, exp, bits, ctx->minfo);
            res->coeffs[actual++] = c;
            flint_free(exp);
        }
        fq_nmod_clear(coeff, fqctx);
    }
    res->length = actual;
    gf28_mpoly_ctx_clear(ctx);
}

void gf28_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf28_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx)
{
    fq_nmod_mpoly_zero(res, fq_mpoly_ctx);
    if (poly->length == 0) return;
    slong N = mpoly_words_per_exp(poly->bits, fq_mpoly_ctx->minfo);
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    for (slong i = 0; i < poly->length; i++) {
        if (poly->coeffs[i] != 0) {
            fq_nmod_t coeff; fq_nmod_init(coeff, fqctx);
            gf28_elem_to_fq_nmod(coeff, poly->coeffs[i], fqctx);
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
            mpoly_get_monomial_ui(exp, poly->exps+N*i, poly->bits, fq_mpoly_ctx->minfo);
            fq_nmod_mpoly_set_coeff_fq_nmod_ui(res, coeff, exp, fq_mpoly_ctx);
            flint_free(exp); fq_nmod_clear(coeff, fqctx);
        }
    }
}

/* ============================================================================
   HELPER: set up the GF(2^8) context used throughout this file.
   Irreducible poly: x^8 + x^4 + x^3 + x^2 + 1  (same as gf28_mpoly_mul_flint_univar)
   ============================================================================ */
static void _gf28_fq_ctx_init(fq_nmod_ctx_t fq_ctx)
{
    nmod_poly_t modulus;
    nmod_poly_init(modulus, 2);
    nmod_poly_set_coeff_ui(modulus, 0, 1);
    nmod_poly_set_coeff_ui(modulus, 2, 1);
    nmod_poly_set_coeff_ui(modulus, 3, 1);
    nmod_poly_set_coeff_ui(modulus, 4, 1);
    nmod_poly_set_coeff_ui(modulus, 8, 1);
    fq_nmod_ctx_init_modulus(fq_ctx, modulus, "z");
    nmod_poly_clear(modulus);
}

/* ============================================================================
   HELPER: estimate worst-case array_size for a given nvars and exp_bound.
   Each mults[i] <= 1 + 2*(exp_bound-1) = 2*exp_bound - 1.
   Returns the product, or SLONG_MAX on overflow.
   ============================================================================ */
static slong _estimate_array_size(slong nvars, slong exp_bound)
{
    slong max_mult = 2 * exp_bound - 1;   /* upper bound on each mults[i] */
    if (max_mult <= 0) return 1;
    double sz = 1.0;
    for (slong i = 0; i < nvars; i++) {
        sz *= (double)max_mult;
        if (sz > (double)(1L << 30)) return LONG_MAX;
    }
    return (slong)sz;
}

/* Choose the largest exp_bound such that _estimate_array_size(nvars, eb) < limit. */
static slong _safe_exp_bound(slong nvars, slong requested_eb, slong array_limit)
{
    if (_estimate_array_size(nvars, requested_eb) < array_limit)
        return requested_eb;
    /* Binary-search downward */
    slong lo = 2, hi = requested_eb;
    while (lo < hi) {
        slong mid = (lo + hi + 1) / 2;
        if (_estimate_array_size(nvars, mid) < array_limit) lo = mid;
        else hi = mid - 1;
    }
    return lo;
}

static uint8_t gf28_mpoly_get_coeff_at(
    const gf28_mpoly_t poly,
    const ulong       *exp,
    const gf28_mpoly_ctx_t ctx)
{
    if (poly->length == 0) return 0;
    slong N  = mpoly_words_per_exp(poly->bits, ctx->minfo);
    ulong *pe = (ulong *)flint_malloc(N * sizeof(ulong));
    mpoly_set_monomial_ui(pe, exp, poly->bits, ctx->minfo);
    uint8_t result = 0;
    for (slong pos = 0; pos < poly->length; pos++) {
        if (mpoly_monomial_equal(poly->exps + N * pos, pe, N)) {
            result = poly->coeffs[pos];
            break;
        }
    }
    flint_free(pe);
    return result;
}

/* ============================================================================
   SCHOOLBOOK (naive) reference multiplication — O(n^2 * nvars), correct by
   construction.  Used as ground truth in the vs-flint test when a simpler
   reference is needed and to cross-check the conversion path.
   ============================================================================ */
static void _gf28_mpoly_mul_schoolbook(
    gf28_mpoly_t       res,
    const gf28_mpoly_t A,
    const gf28_mpoly_t B,
    const gf28_mpoly_ctx_t ctx)
{
    if (!g_gf28_complete_tables.initialized)
        init_gf28_standard();

    slong nvars = ctx->minfo->nvars;
    gf28_mpoly_zero(res, ctx);
    if (A->length == 0 || B->length == 0) return;

    slong NA = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    ulong *ea = (ulong *)flint_malloc(nvars * sizeof(ulong));
    ulong *eb = (ulong *)flint_malloc(nvars * sizeof(ulong));
    ulong *ec = (ulong *)flint_malloc(nvars * sizeof(ulong));

    for (slong i = 0; i < A->length; i++) {
        mpoly_get_monomial_ui(ea, A->exps + NA * i, A->bits, ctx->minfo);
        uint8_t ca = A->coeffs[i];
        for (slong j = 0; j < B->length; j++) {
            mpoly_get_monomial_ui(eb, B->exps + NB * j, B->bits, ctx->minfo);
            uint8_t cprod = gf28_mul(ca, B->coeffs[j]);
            if (cprod == 0) continue;
            for (slong k = 0; k < nvars; k++) ec[k] = ea[k] + eb[k];
            uint8_t existing  = gf28_mpoly_get_coeff_at(res, ec, ctx);
            uint8_t new_coeff = existing ^ cprod; 
            gf28_mpoly_set_coeff_ui_ui(res, new_coeff, ec, ctx);
        }
    }
    flint_free(ea);
    flint_free(eb);
    flint_free(ec);
}

