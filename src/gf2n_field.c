/* gf2n_field.c - Implementation of GF(2^n) field operations */
#include "gf2n_field.h"


/* ============================================================================
   CPU FEATURE DETECTION
   ============================================================================ */

/* Check CPU features */
int has_pclmulqdq(void) {
#if DIXON_X86_SIMD
    unsigned int eax, ebx, ecx, edx;
    if (__get_cpuid(1, &eax, &ebx, &ecx, &edx)) {
        return (ecx & bit_PCLMUL) != 0;
    }
#endif
    return 0;
}

/* Global multiplication function pointers */
gf232_mul_func gf232_mul = NULL;
gf264_mul_func gf264_mul = NULL;
gf2128_mul_func gf2128_mul = NULL;

/* ============================================================================
   GF(2^8) OPERATIONS IMPLEMENTATION
   ============================================================================ */

/* Initialize GF(2^8) lookup tables */
void init_gf28_tables(uint8_t irred_poly) {
    if (g_gf28_tables) return;
    
    g_gf28_tables = (gf28_tables_t *)malloc(sizeof(gf28_tables_t));
    
    memset(g_gf28_tables->exp_table, 0, sizeof(g_gf28_tables->exp_table));
    memset(g_gf28_tables->log_table, 0, sizeof(g_gf28_tables->log_table));
    memset(g_gf28_tables->inv_table, 0, sizeof(g_gf28_tables->inv_table));
    
    /* Find a primitive element (generator) */
    uint8_t generator = 0;
    for (uint8_t g = 2; g < 256; g++) {
        uint16_t alpha = 1;
        int is_primitive = 1;
        
        for (int i = 0; i < 255; i++) {
            uint16_t new_alpha = 0;
            uint16_t a = alpha;
            uint8_t b = g;
            
            while (b) {
                if (b & 1) new_alpha ^= a;
                a <<= 1;
                if (a & 0x100) a ^= (0x100 | irred_poly);
                b >>= 1;
            }
            
            alpha = new_alpha & 0xFF;
            
            if (i < 254 && alpha == 1) {
                is_primitive = 0;
                break;
            }
        }
        
        if (is_primitive && alpha == 1) {
            generator = g;
            break;
        }
    }
    
    g_gf28_tables->generator = generator;
    
    /* Generate lookup tables using the primitive element */
    uint16_t alpha = 1;
    for (int i = 0; i < 255; i++) {
        g_gf28_tables->exp_table[i] = (uint8_t)alpha;
        g_gf28_tables->log_table[(uint8_t)alpha] = i;
        
        uint16_t new_alpha = 0;
        uint16_t a = alpha;
        uint8_t b = generator;
        
        while (b) {
            if (b & 1) new_alpha ^= a;
            a <<= 1;
            if (a & 0x100) a ^= (0x100 | irred_poly);
            b >>= 1;
        }
        
        alpha = new_alpha & 0xFF;
    }
    
    for (int i = 0; i < 255; i++) {
        g_gf28_tables->exp_table[i + 255] = g_gf28_tables->exp_table[i];
    }
    
    g_gf28_tables->exp_table[255] = 1;
    g_gf28_tables->log_table[0] = 0;
    
    g_gf28_tables->inv_table[0] = 0;
    for (int i = 1; i < 256; i++) {
        int log_inv = 255 - g_gf28_tables->log_table[i];
        g_gf28_tables->inv_table[i] = g_gf28_tables->exp_table[log_inv];
    }
}

/* Initialize complete lookup tables */
void init_gf28_complete_tables(void) {
    if (g_gf28_complete_tables.initialized) return;
    
    if (!g_gf28_tables) {
        init_gf28_tables(0x1D);
    }
    
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            if (i == 0 || j == 0) {
                g_gf28_complete_tables.mul_table[i][j] = 0;
            } else {
                int log_sum = g_gf28_tables->log_table[i] + g_gf28_tables->log_table[j];
                g_gf28_complete_tables.mul_table[i][j] = g_gf28_tables->exp_table[log_sum];
            }
        }
        
        g_gf28_complete_tables.sqr_table[i] = g_gf28_complete_tables.mul_table[i][i];
        g_gf28_complete_tables.inv_table[i] = g_gf28_tables->inv_table[i];
        g_gf28_complete_tables.double_table[i] = g_gf28_complete_tables.mul_table[2][i];
    }
    
    g_gf28_complete_tables.initialized = 1;
}

void init_gf28_standard(void) {
    if (g_gf28_tables) return;
    init_gf28_tables(0x1D);
    init_gf28_complete_tables();
}

void cleanup_gf28_tables(void) {
    if (g_gf28_tables) {
        free(g_gf28_tables);
        g_gf28_tables = NULL;
    }
}

/* GF(2^8) arithmetic operations */
inline gf28_t gf28_add(gf28_t a, gf28_t b) {
    return a ^ b;
}

inline gf28_t gf28_sub(gf28_t a, gf28_t b) {
    return a ^ b;  // In GF(2^8), subtraction is the same as addition (XOR)
}

inline gf28_t gf28_mul(gf28_t a, gf28_t b) {
    return g_gf28_complete_tables.mul_table[a][b];
}

inline const uint8_t* gf28_get_scalar_row(uint8_t scalar) {
    return g_gf28_complete_tables.mul_table[scalar];
}

inline gf28_t gf28_sqr(gf28_t a) {
    return g_gf28_complete_tables.sqr_table[a];
}

inline gf28_t gf28_inv(gf28_t a) {
    return g_gf28_complete_tables.inv_table[a];
}

inline gf28_t gf28_div(gf28_t a, gf28_t b) {
    if (b == 0) return 0;
    return gf28_mul(a, gf28_inv(b));
}

/* ============================================================================
   GF(2^16) OPERATIONS IMPLEMENTATION
   ============================================================================ */

void init_gf216_tables(uint16_t irred_poly) {
    if (g_gf216_tables) return;
    
    g_gf216_tables = (gf216_tables_t *)malloc(sizeof(gf216_tables_t));
    
    memset(g_gf216_tables->exp_table, 0, sizeof(g_gf216_tables->exp_table));
    memset(g_gf216_tables->log_table, 0, sizeof(g_gf216_tables->log_table));
    memset(g_gf216_tables->inv_table, 0, sizeof(g_gf216_tables->inv_table));
    
    uint16_t generator = 0;
    for (uint16_t g = 2; g < 65536; g++) {
        uint32_t alpha = 1;
        int is_primitive = 1;
        
        for (int i = 0; i < 65535; i++) {
            uint32_t new_alpha = 0;
            uint32_t a = alpha;
            uint16_t b = g;
            
            while (b) {
                if (b & 1) new_alpha ^= a;
                a <<= 1;
                if (a & 0x10000) a ^= (0x10000 | irred_poly);
                b >>= 1;
            }
            
            alpha = new_alpha & 0xFFFF;
            
            if (i < 65534 && alpha == 1) {
                is_primitive = 0;
                break;
            }
        }
        
        if (is_primitive && alpha == 1) {
            generator = g;
            break;
        }
    }
    
    g_gf216_tables->generator = generator;
    
    uint32_t alpha = 1;
    for (int i = 0; i < 65535; i++) {
        g_gf216_tables->exp_table[i] = (uint16_t)alpha;
        g_gf216_tables->log_table[(uint16_t)alpha] = i;
        
        uint32_t new_alpha = 0;
        uint32_t a = alpha;
        uint16_t b = generator;
        
        while (b) {
            if (b & 1) new_alpha ^= a;
            a <<= 1;
            if (a & 0x10000) a ^= (0x10000 | irred_poly);
            b >>= 1;
        }
        
        alpha = new_alpha & 0xFFFF;
    }
    
    for (int i = 0; i < 65535; i++) {
        g_gf216_tables->exp_table[i + 65535] = g_gf216_tables->exp_table[i];
    }
    
    g_gf216_tables->exp_table[65535] = 1;
    g_gf216_tables->log_table[0] = 0;
    
    g_gf216_tables->inv_table[0] = 0;
    for (int i = 1; i < 65536; i++) {
        int log_inv = 65535 - g_gf216_tables->log_table[i];
        g_gf216_tables->inv_table[i] = g_gf216_tables->exp_table[log_inv];
    }
}

void init_gf216_complete_tables(void) {
    if (g_gf216_complete_tables.initialized) return;
    
    if (!g_gf216_tables) {
        init_gf216_tables(0x002D);
    }
    
    for (int i = 0; i < 65536; i++) {
        if (i == 0) {
            g_gf216_complete_tables.sqr_table[i] = 0;
        } else {
            int log_sq = (2 * g_gf216_tables->log_table[i]) % 65535;
            g_gf216_complete_tables.sqr_table[i] = g_gf216_tables->exp_table[log_sq];
        }
        
        g_gf216_complete_tables.inv_table[i] = g_gf216_tables->inv_table[i];
        
        if (i == 0) {
            g_gf216_complete_tables.double_table[i] = 0;
        } else {
            int log_2 = g_gf216_tables->log_table[2];
            int log_prod = (log_2 + g_gf216_tables->log_table[i]) % 65535;
            g_gf216_complete_tables.double_table[i] = g_gf216_tables->exp_table[log_prod];
        }
    }
    
    g_gf216_complete_tables.initialized = 1;
}

void init_gf216_standard(void) {
    if (g_gf216_tables) return;
    init_gf216_tables(0x002D);
    init_gf216_complete_tables();
}

void cleanup_gf216_tables(void) {
    if (g_gf216_tables) {
        free(g_gf216_tables);
        g_gf216_tables = NULL;
    }
    if (g_gf216_complete_tables.mul_table) {
        free(g_gf216_complete_tables.mul_table);
        g_gf216_complete_tables.mul_table = NULL;
    }
}

void print_gf216_memory_usage(void) {
    size_t total = 0;
    
    printf("=== GF(2^16) Memory Usage ===\n");
    
    size_t basic = sizeof(gf216_tables_t);
    printf("Basic tables (log/exp/inv): %zu bytes (%.1f MB)\n", basic, basic/1024.0/1024.0);
    total += basic;
    
    size_t other = sizeof(g_gf216_complete_tables) - sizeof(uint16_t*);
    printf("Other tables (sqr/inv/double): %zu bytes (%.1f KB)\n", other, other/1024.0);
    total += other;
    
    printf("Total memory: %zu bytes (%.1f MB)\n", total, total/1024.0/1024.0);
}

inline gf216_t gf216_add(gf216_t a, gf216_t b) {
    return a ^ b;
}

inline gf216_t gf216_mul(gf216_t a, gf216_t b) {
    if (a == 0 || b == 0) return 0;
    int log_sum = g_gf216_tables->log_table[a] + g_gf216_tables->log_table[b];
    return g_gf216_tables->exp_table[log_sum];
}

inline gf216_t gf216_sqr(gf216_t a) {
    return g_gf216_complete_tables.sqr_table[a];
}

inline gf216_t gf216_inv(gf216_t a) {
    return g_gf216_complete_tables.inv_table[a];
}

inline gf216_t gf216_div(gf216_t a, gf216_t b) {
    if (b == 0) return 0;
    return gf216_mul(a, gf216_inv(b));
}

/* ============================================================================
   GF(2^32) OPERATIONS IMPLEMENTATION
   ============================================================================ */

/* GF(2^32): x^32 + x^15 + x^9 + x^7 + x^4 + x^3 + 1 = 0x8299 */
#define GF232_MODULUS 0x8299

/* Basic operations */
inline gf232_t gf232_create(uint32_t val) {
    gf232_t result = {val};
    return result;
}

inline gf232_t gf232_zero(void) {
    return gf232_create(0);
}

inline gf232_t gf232_one(void) {
    return gf232_create(1);
}

inline int gf232_is_zero(const gf232_t *a) {
    return a->value == 0;
}

inline int gf232_equal(const gf232_t *a, const gf232_t *b) {
    return a->value == b->value;
}

inline gf232_t gf232_add(const gf232_t *a, const gf232_t *b) {
    gf232_t result;
    result.value = a->value ^ b->value;
    return result;
}

void gf232_print(const gf232_t *a) {
    printf("%08x", a->value);
}



/* Software multiplication for GF(2^32) with new polynomial */
gf232_t gf232_mul_software(const gf232_t *a, const gf232_t *b) {
    uint64_t result = 0;
    uint64_t temp = a->value;
    uint32_t mask = b->value;
    
    while (mask) {
        if (mask & 1) {
            result ^= temp;
        }
        temp <<= 1;
        mask >>= 1;
    }
    
    /* Reduction modulo x^32 + x^15 + x^9 + x^7 + x^4 + x^3 + 1 */
    for (int i = 63; i >= 32; i--) {
        if ((result >> i) & 1) {
            result ^= ((uint64_t)GF232_MODULUS << (i - 32));
        }
    }
    
    return gf232_create((uint32_t)result);
}

/* Optimized PCLMUL-based multiplication for GF(2^32) with new polynomial */
gf232_t gf232_mul_pclmul(const gf232_t *a, const gf232_t *b) {
#if DIXON_X86_SIMD
    __m128i x = _mm_set_epi64x(0, a->value);
    __m128i y = _mm_set_epi64x(0, b->value);
    
    /* 32x32 -> 64 bit multiplication */
    __m128i prod = _mm_clmulepi64_si128(x, y, 0x00);
    
    uint64_t result = _mm_extract_epi64(prod, 0);
    
    /* Optimized reduction for x^32 + x^15 + x^9 + x^7 + x^4 + x^3 + 1 */
    uint64_t hi = result >> 32;  // Upper 32 bits
    uint64_t lo = result & 0xFFFFFFFF;  // Lower 32 bits
    
    /* For GF(2^32) with polynomial x^32 + x^15 + x^9 + x^7 + x^4 + x^3 + 1:
     * x^32 ≡ x^15 + x^9 + x^7 + x^4 + x^3 + 1 = 0x8299
     */
    uint64_t fold = 0;
    uint32_t h = (uint32_t)hi;
    
    /* Compute h * 0x8299 using shifts and XORs */
    fold = h;
    fold ^= (uint64_t)h << 3;
    fold ^= (uint64_t)h << 4;
    fold ^= (uint64_t)h << 7;
    fold ^= (uint64_t)h << 9;
    fold ^= (uint64_t)h << 15;
    
    result = lo ^ (fold & 0xFFFFFFFF);
    
    uint32_t fold_hi = fold >> 32;
    if (fold_hi) {
        uint32_t extra = fold_hi;
        extra ^= fold_hi << 3;
        extra ^= fold_hi << 4;
        extra ^= fold_hi << 7;
        extra ^= fold_hi << 9;
        extra ^= fold_hi << 15;
        result ^= extra;
    }
    
    return gf232_create((uint32_t)result);
#else
    return gf232_mul_software(a, b);
#endif
}

/* Fast squaring in GF(2^32) with new polynomial */
gf232_t gf232_sqr(const gf232_t *a) {
    uint64_t result = 0;
    uint32_t val = a->value;
    
    /* Expand bits */
    for (int i = 0; i < 16; i++) {
        if ((val >> i) & 1) {
            result |= (uint64_t)1 << (2 * i);
        }
    }
    for (int i = 16; i < 32; i++) {
        if ((val >> i) & 1) {
            result |= (uint64_t)1 << (2 * i);
        }
    }
    
    /* Reduction modulo x^32 + x^15 + x^9 + x^7 + x^4 + x^3 + 1 */
    for (int i = 63; i >= 32; i--) {
        if ((result >> i) & 1) {
            result ^= ((uint64_t)GF232_MODULUS << (i - 32));
        }
    }
    
    return gf232_create((uint32_t)result);
}
/* Initialize multiplication function */
void init_gf232(void) {
    if (!gf232_mul) {
        if (has_pclmulqdq()) {
            gf232_mul = gf232_mul_pclmul;
            printf("Using PCLMULQDQ for GF(2^32) multiplication\n");
        } else {
            gf232_mul = gf232_mul_software;
            printf("Using software implementation for GF(2^32) multiplication\n");
        }
    }
}


/* Inversion using repeated squaring */
gf232_t gf232_inv(const gf232_t *a) {
    if (gf232_is_zero(a)) {
        return gf232_zero();
    }
    
    // In GF(2^32), a^(2^32-1) = 1, so a^(2^32-2) = a^(-1)
    // Using square-and-multiply algorithm
    
    gf232_t a_pow = *a;
    
    // Start from a^2
    a_pow = gf232_sqr(&a_pow);
    gf232_t result = a_pow;
    
    // Continue with a^4, a^8, ..., a^(2^31)
    for (int i = 2; i < 32; i++) {
        a_pow = gf232_sqr(&a_pow);
        result = gf232_mul(&result, &a_pow);
    }
    
    return result;
}

/* Division in GF(2^32) */
inline gf232_t gf232_div(const gf232_t *a, const gf232_t *b) {
    gf232_t b_inv = gf232_inv(b);
    return gf232_mul(a, &b_inv);
}

/* ============================================================================
   GF(2^64) OPERATIONS IMPLEMENTATION
   ============================================================================ */

/* GF(2^64): x^64 + x^33 + x^30 + x^26 + x^25 + x^24 + x^23 + x^22 + x^21 + 
             x^20 + x^18 + x^13 + x^12 + x^11 + x^10 + x^7 + x^5 + x^4 + 
             x^2 + x + 1 = 0x247F43CB7ULL */
#define GF264_MODULUS 0x247F43CB7ULL


/* Basic operations */
inline gf264_t gf264_create(uint64_t val) {
    gf264_t result = {val};
    return result;
}

inline gf264_t gf264_zero(void) {
    return gf264_create(0);
}

inline gf264_t gf264_one(void) {
    return gf264_create(1);
}

inline int gf264_is_zero(const gf264_t *a) {
    return a->value == 0;
}

inline int gf264_equal(const gf264_t *a, const gf264_t *b) {
    return a->value == b->value;
}

inline gf264_t gf264_add(const gf264_t *a, const gf264_t *b) {
    gf264_t result;
    result.value = a->value ^ b->value;
    return result;
}

void gf264_print(const gf264_t *a) {
    printf("%016llx", (unsigned long long)a->value);
}

/* Helper function to reduce a 128-bit value modulo the GF(2^64) polynomial 
 * This implements Barrett/Montgomery-style reduction for GF(2^n)
 */
static void gf264_reduce_128(uint64_t *lo, uint64_t *hi) {
    /* The modulus is: x^64 + 0x247F43CB7
     * So x^64 ≡ 0x247F43CB7 (the lower 64-bit terms)
     * 
     * Strategy: Process hi bits from high to low
     * For each bit k in hi: x^(64+k) ≡ x^k * 0x247F43CB7
     */
    
    uint64_t h = *hi;
    uint64_t l = *lo;
    
    /* Process all 64 bits of hi
     * We go from bit 63 down to bit 0 */
    for (int k = 63; k >= 0; k--) {
        if ((h >> k) & 1ULL) {
            /* Clear this bit - we're reducing x^(64+k) */
            h ^= (1ULL << k);
            
            /* Replace with x^k * 0x247F43CB7 
             * Need to compute: 0x247F43CB7 << k and split across 128 bits */
            
            if (k == 0) {
                /* x^64 becomes 0x247F43CB7 */
                l ^= GF264_MODULUS;
            } else if (k < 31) {
                /* GF264_MODULUS << k fits entirely in lower 64 bits */
                l ^= (GF264_MODULUS << k);
            } else if (k == 31) {
                /* GF264_MODULUS << 31: bottom 31 bits go to l, top bits to h */
                l ^= (GF264_MODULUS << 31);
                h ^= (GF264_MODULUS >> 33);  /* top 3 bits */
            } else {
                /* k > 31: result spans both parts */
                l ^= (GF264_MODULUS << k);
                h ^= (GF264_MODULUS >> (64 - k));
            }
        }
    }
    
    *lo = l;
    *hi = h;
}

/* Software multiplication for GF(2^64) */
gf264_t gf264_mul_software(const gf264_t *a, const gf264_t *b) {
    uint64_t result_low = 0, result_high = 0;
    uint64_t a_val = a->value;
    uint64_t b_val = b->value;
    
    /* Step 1: 64x64 -> 128 bit carryless multiplication */
    for (int i = 0; i < 64; i++) {
        if ((b_val >> i) & 1) {
            result_low ^= (a_val << i);
            if (i > 0) {
                result_high ^= (a_val >> (64 - i));
            }
        }
    }
    
    /* Step 2: Reduce modulo the polynomial */
    gf264_reduce_128(&result_low, &result_high);
    
    return gf264_create(result_low);
}

/* PCLMUL-based multiplication for GF(2^64) */
gf264_t gf264_mul_pclmul(const gf264_t *a, const gf264_t *b) {
#if DIXON_X86_SIMD
    __m128i x = _mm_set_epi64x(0, a->value);
    __m128i y = _mm_set_epi64x(0, b->value);
    
    /* 64x64 -> 128 bit multiplication */
    __m128i prod = _mm_clmulepi64_si128(x, y, 0x00);
    
    uint64_t lo = _mm_extract_epi64(prod, 0);
    uint64_t hi = _mm_extract_epi64(prod, 1);
    
    /* Reduce using the helper function */
    gf264_reduce_128(&lo, &hi);
    
    return gf264_create(lo);
#else
    return gf264_mul_software(a, b);
#endif
}

/* Fast squaring in GF(2^64) */
gf264_t gf264_sqr(const gf264_t *a) {
    uint64_t result_low = 0, result_high = 0;
    uint64_t val = a->value;
    
    /* Expand bits: a^2 spreads bits apart
     * bit i goes to bit 2i */
    for (int i = 0; i < 32; i++) {
        if ((val >> i) & 1) {
            result_low |= (1ULL << (2 * i));
        }
    }
    for (int i = 32; i < 64; i++) {
        if ((val >> i) & 1) {
            result_high |= (1ULL << (2 * (i - 32)));
        }
    }
    
    /* Reduce modulo the polynomial */
    gf264_reduce_128(&result_low, &result_high);
    
    return gf264_create(result_low);
}



/* Initialize multiplication function */
void init_gf264(void) {
    if (!gf264_mul) {
        if (has_pclmulqdq()) {
            gf264_mul = gf264_mul_pclmul;
            printf("Using PCLMULQDQ for GF(2^64) multiplication\n");
        } else {
            gf264_mul = gf264_mul_software;
            printf("Using software implementation for GF(2^64) multiplication\n");
        }
    }
}

/* Inversion using addition chain optimized for GF(2^64) */
gf264_t gf264_inv(const gf264_t *a) {
    if (gf264_is_zero(a)) {
        return gf264_zero();
    }
    
    /* In GF(2^64), we need to compute a^(2^64-2) = a^(0xFFFFFFFFFFFFFFFE) */
    /* a^(2^64-2) = a^2 * a^4 * a^8 * ... * a^(2^63) */
    
    gf264_t result = *a;
    gf264_t a_power = *a;
    
    /* Start with a^2 */
    a_power = gf264_sqr(&a_power);
    result = a_power;
    
    /* Multiply by a^4, a^8, ..., a^(2^63) */
    for (int i = 2; i <= 63; i++) {
        a_power = gf264_sqr(&a_power);  /* a_power = a^(2^i) */
        result = gf264_mul(&result, &a_power);
    }
    
    return result;
}

/* Division in GF(2^64) */
inline gf264_t gf264_div(const gf264_t *a, const gf264_t *b) {
    gf264_t b_inv = gf264_inv(b);
    return gf264_mul(a, &b_inv);
}

/* ============================================================================
   GF(2^128) OPERATIONS IMPLEMENTATION
   ============================================================================ */

/* Basic operations */
inline gf2128_t gf2128_create(uint64_t low, uint64_t high) {
    gf2128_t result = {low, high};
    return result;
}

inline gf2128_t gf2128_zero(void) {
    return gf2128_create(0, 0);
}

inline gf2128_t gf2128_one(void) {
    return gf2128_create(1, 0);
}

inline int gf2128_is_zero(const gf2128_t *a) {
    return (a->low == 0 && a->high == 0);
}

inline int gf2128_equal(const gf2128_t *a, const gf2128_t *b) {
    return (a->low == b->low && a->high == b->high);
}

inline gf2128_t gf2128_add(const gf2128_t *a, const gf2128_t *b) {
    gf2128_t result;
    result.low = a->low ^ b->low;
    result.high = a->high ^ b->high;
    return result;
}

void gf2128_print(const gf2128_t *a) {
    printf("%016lx%016lx", a->high, a->low);
}

/* Software multiplication for GF(2^128) */
gf2128_t gf2128_mul_software(const gf2128_t *a, const gf2128_t *b) {
    uint64_t a0 = a->low, a1 = a->high;
    uint64_t b0 = b->low, b1 = b->high;
    uint64_t z0 = 0, z1 = 0, z2 = 0, z3 = 0;
    
    /* 128x128 bit multiplication using schoolbook method */
    for (int i = 0; i < 64; i++) {
        if ((b0 >> i) & 1) {
            z0 ^= a0 << i;
            if (i > 0) {
                z1 ^= a0 >> (64 - i);
            }
        }
    }
    
    for (int i = 0; i < 64; i++) {
        if ((b0 >> i) & 1) {
            z1 ^= a1 << i;
            if (i > 0) {
                z2 ^= a1 >> (64 - i);
            }
        }
    }
    
    for (int i = 0; i < 64; i++) {
        if ((b1 >> i) & 1) {
            z1 ^= a0 << i;
            if (i > 0) {
                z2 ^= a0 >> (64 - i);
            }
        }
    }
    
    for (int i = 0; i < 64; i++) {
        if ((b1 >> i) & 1) {
            z2 ^= a1 << i;
            if (i > 0) {
                z3 ^= a1 >> (64 - i);
            }
        }
    }
    
    /* Reduction modulo x^128 + x^7 + x^2 + x + 1 */
    /* Reduce z3 (bits 192-255) */
    for (int i = 63; i >= 0; i--) {
        if ((z3 >> i) & 1) {
            int pos = 64 + i;
            z3 ^= (1ULL << i);
            
            if (pos + 7 < 128) {
                z1 ^= (1ULL << (pos + 7 - 64));
            } else {
                z2 ^= (1ULL << (pos + 7 - 128));
            }
            
            if (pos + 2 < 128) {
                z1 ^= (1ULL << (pos + 2 - 64));
            } else {
                z2 ^= (1ULL << (pos + 2 - 128));
            }
            
            if (pos + 1 < 128) {
                z1 ^= (1ULL << (pos + 1 - 64));
            } else {
                z2 ^= (1ULL << (pos + 1 - 128));
            }
            
            if (pos < 128) {
                z1 ^= (1ULL << (pos - 64));
            } else {
                z2 ^= (1ULL << (pos - 128));
            }
        }
    }
    
    /* Reduce z2 (bits 128-191) */
    for (int i = 63; i >= 0; i--) {
        if ((z2 >> i) & 1) {
            int pos = i;
            z2 ^= (1ULL << i);
            
            if (pos + 7 < 64) {
                z0 ^= (1ULL << (pos + 7));
            } else {
                z1 ^= (1ULL << (pos + 7 - 64));
            }
            
            if (pos + 2 < 64) {
                z0 ^= (1ULL << (pos + 2));
            } else {
                z1 ^= (1ULL << (pos + 2 - 64));
            }
            
            if (pos + 1 < 64) {
                z0 ^= (1ULL << (pos + 1));
            } else {
                z1 ^= (1ULL << (pos + 1 - 64));
            }
            
            if (pos < 64) {
                z0 ^= (1ULL << pos);
            } else {
                z1 ^= (1ULL << (pos - 64));
            }
        }
    }
    
    return gf2128_create(z0, z1);
}

/* CLMUL-based multiplication for GF(2^128) */
inline gf2128_t gf2128_mul_clmul(const gf2128_t *a, const gf2128_t *b) {
#if DIXON_X86_SIMD
    __m128i x = _mm_set_epi64x(a->high, a->low);
    __m128i y = _mm_set_epi64x(b->high, b->low);
    
    /* Step 1: 128x128 -> 256 bit multiplication */
    __m128i t0, t1, t2;
    
    t0 = _mm_clmulepi64_si128(x, y, 0x00);
    t1 = _mm_clmulepi64_si128(x, y, 0x11);
    t2 = _mm_xor_si128(
        _mm_clmulepi64_si128(x, y, 0x10),
        _mm_clmulepi64_si128(x, y, 0x01)
    );
    
    __m128i lo = _mm_xor_si128(t0, _mm_slli_si128(t2, 8));
    __m128i hi = _mm_xor_si128(t1, _mm_srli_si128(t2, 8));
    
    const uint64_t g = 0x87;
    
    __m128i tmp = _mm_clmulepi64_si128(hi, _mm_set_epi64x(0, g), 0x01);
    hi = _mm_xor_si128(hi, _mm_srli_si128(tmp, 8));
    lo = _mm_xor_si128(lo, _mm_slli_si128(tmp, 8));
    
    tmp = _mm_clmulepi64_si128(hi, _mm_set_epi64x(0, g), 0x00);
    lo = _mm_xor_si128(lo, tmp);
    
    gf2128_t res;
    res.low = _mm_extract_epi64(lo, 0);
    res.high = _mm_extract_epi64(lo, 1);
    
    return res;
#else
    return gf2128_mul_software(a, b);
#endif
}

/* Initialize multiplication function */
void init_gf2128(void) {
    if (!gf2128_mul) {
        if (has_pclmulqdq()) {
            gf2128_mul = gf2128_mul_clmul;
            printf("Using PCLMULQDQ for GF(2^128) multiplication\n");
        } else {
            gf2128_mul = gf2128_mul_software;
            printf("Using software implementation for GF(2^128) multiplication\n");
        }
    }
}

/* Fast squaring in GF(2^128) */
gf2128_t gf2128_sqr(const gf2128_t *a) {
    uint64_t r0 = 0, r1 = 0, r2 = 0, r3 = 0;
    
    /* Expand low 64 bits */
    for (int i = 0; i < 32; i++) {
        if ((a->low >> i) & 1) {
            r0 |= 1ULL << (2 * i);
        }
    }
    for (int i = 32; i < 64; i++) {
        if ((a->low >> i) & 1) {
            r1 |= 1ULL << (2 * (i - 32));
        }
    }
    
    /* Expand high 64 bits */
    for (int i = 0; i < 32; i++) {
        if ((a->high >> i) & 1) {
            r2 |= 1ULL << (2 * i);
        }
    }
    for (int i = 32; i < 64; i++) {
        if ((a->high >> i) & 1) {
            r3 |= 1ULL << (2 * (i - 32));
        }
    }
    
    /* Reduction modulo x^128 + x^7 + x^2 + x + 1 */
    for (int i = 63; i >= 0; i--) {
        if ((r3 >> i) & 1) {
            int pos = 64 + i;
            r3 ^= (1ULL << i);
            
            if (pos + 7 < 128) {
                r1 ^= (1ULL << (pos + 7 - 64));
            } else {
                r2 ^= (1ULL << (pos + 7 - 128));
            }
            
            if (pos + 2 < 128) {
                r1 ^= (1ULL << (pos + 2 - 64));
            } else {
                r2 ^= (1ULL << (pos + 2 - 128));
            }
            
            if (pos + 1 < 128) {
                r1 ^= (1ULL << (pos + 1 - 64));
            } else {
                r2 ^= (1ULL << (pos + 1 - 128));
            }
            
            if (pos < 128) {
                r1 ^= (1ULL << (pos - 64));
            } else {
                r2 ^= (1ULL << (pos - 128));
            }
        }
    }
    
    for (int i = 63; i >= 0; i--) {
        if ((r2 >> i) & 1) {
            int pos = i;
            r2 ^= (1ULL << i);
            
            if (pos + 7 < 64) {
                r0 ^= (1ULL << (pos + 7));
            } else {
                r1 ^= (1ULL << (pos + 7 - 64));
            }
            
            if (pos + 2 < 64) {
                r0 ^= (1ULL << (pos + 2));
            } else {
                r1 ^= (1ULL << (pos + 2 - 64));
            }
            
            if (pos + 1 < 64) {
                r0 ^= (1ULL << (pos + 1));
            } else {
                r1 ^= (1ULL << (pos + 1 - 64));
            }
            
            if (pos < 64) {
                r0 ^= (1ULL << pos);
            } else {
                r1 ^= (1ULL << (pos - 64));
            }
        }
    }
    
    return gf2128_create(r0, r1);
}

/* Forward declaration for conversion function */
gf2128_t fq_nmod_to_gf2128(const fq_nmod_t elem, const fq_nmod_ctx_t ctx);

/* Inversion using repeated squaring - optimized for GF(2^128) */
gf2128_t gf2128_inv(const gf2128_t *a) {
    if (gf2128_is_zero(a)) {
        return gf2128_zero();
    }
    
    /* For debugging: use FLINT to compute the correct inverse */
    /* Create GF(2^128) context */
    fq_nmod_ctx_t ctx;
    nmod_poly_t mod;
    nmod_poly_init(mod, 2);
    
    /* Set modulus to x^128 + x^7 + x^2 + x + 1 */
    nmod_poly_set_coeff_ui(mod, 0, 1);
    nmod_poly_set_coeff_ui(mod, 1, 1);
    nmod_poly_set_coeff_ui(mod, 2, 1);
    nmod_poly_set_coeff_ui(mod, 7, 1);
    nmod_poly_set_coeff_ui(mod, 128, 1);
    
    fq_nmod_ctx_init_modulus(ctx, mod, "t");
    
    /* Convert to FLINT format */
    fq_nmod_t a_flint, inv_flint;
    fq_nmod_init(a_flint, ctx);
    fq_nmod_init(inv_flint, ctx);
    
    gf2128_to_fq_nmod(a_flint, a, ctx);
    
    /* Compute inverse using FLINT */
    fq_nmod_inv(inv_flint, a_flint, ctx);
    
    /* Convert back */
    gf2128_t result = fq_nmod_to_gf2128(inv_flint, ctx);
    
    /* Cleanup */
    fq_nmod_clear(a_flint, ctx);
    fq_nmod_clear(inv_flint, ctx);
    fq_nmod_ctx_clear(ctx);
    nmod_poly_clear(mod);
    
    return result;
}

/* Division in GF(2^128) */
inline gf2128_t gf2128_div(const gf2128_t *a, const gf2128_t *b) {
    gf2128_t b_inv = gf2128_inv(b);
    return gf2128_mul(a, &b_inv);
}

/* ============================================================================
   UTILITY FUNCTIONS IMPLEMENTATION
   ============================================================================ */

inline int _fq_nmod_ctx_is_gf2n(const fq_nmod_ctx_t ctx) {
   const nmod_poly_struct *mod = fq_nmod_ctx_modulus(ctx);
   return (mod->mod.n == 2);
}

uint64_t extract_irred_poly(const fq_nmod_ctx_t ctx) {
   const nmod_poly_struct *mod = fq_nmod_ctx_modulus(ctx);
   uint64_t irred = 0;
   
   for (slong i = 0; i <= nmod_poly_degree(mod); i++) {
       if (nmod_poly_get_coeff_ui(mod, i))
           irred |= (1ULL << i);
   }
   
   return irred;
}

void extract_gf264_poly(const fq_nmod_ctx_t ctx, uint64_t *poly_low, uint64_t *poly_high) {
    const nmod_poly_struct *mod = fq_nmod_ctx_modulus(ctx);
    *poly_low = 0;
    *poly_high = 0;
    
    for (slong i = 0; i < 64 && i <= nmod_poly_degree(mod); i++) {
        if (nmod_poly_get_coeff_ui(mod, i)) {
            *poly_low |= (1ULL << i);
        }
    }
    
    for (slong i = 64; i <= nmod_poly_degree(mod); i++) {
        if (nmod_poly_get_coeff_ui(mod, i)) {
            *poly_high |= (1ULL << (i - 64));
        }
    }
}

void gf2128_to_hex(const gf2128_t *a, char *buf) {
   sprintf(buf, "%016llx%016llx", 
           (unsigned long long)a->high, 
           (unsigned long long)a->low);
}

int gf2128_from_hex(gf2128_t *a, const char *hex) {
   if (strlen(hex) > 32) return -1;
   
   char high_str[17] = {0};
   char low_str[17] = {0};
   
   size_t len = strlen(hex);
   if (len > 16) {
       strncpy(high_str, hex, len - 16);
       strncpy(low_str, hex + len - 16, 16);
   } else {
       strncpy(low_str, hex, len);
   }
   
   uint64_t high = strtoull(high_str, NULL, 16);
   uint64_t low = strtoull(low_str, NULL, 16);
   
   *a = gf2128_create(low, high);
   
   return 0;
}

/* Initialize all GF(2^n) fields */
void init_all_gf2n_fields(void) {
    init_gf28_standard();
    init_gf216_standard();
    init_gf232();
    init_gf264();
    init_gf2128();
}

void cleanup_all_gf2n_fields(void) {
    cleanup_gf28_tables();
    cleanup_gf28_conversion();
    cleanup_gf216_tables();
    cleanup_gf216_conversion();
    cleanup_gf232_conversion();
    cleanup_gf264_conversion();
}

/* ============================================================================
   CONVERSION FUNCTIONS IMPLEMENTATION
   ============================================================================ */

/* Initialize GF(2^128) conversion - stub implementation */
void init_gf2128_conversion(const fq_nmod_ctx_t ctx) {
    if (g_gf2128_conversion && g_gf2128_conversion->initialized) {
        return;
    }
    
    if (!g_gf2128_conversion) {
        g_gf2128_conversion = (gf2128_conversion_t *)calloc(1, sizeof(gf2128_conversion_t));
    }
    
    // For now, just mark as initialized
    g_gf2128_conversion->initialized = 1;
}

/* Convert between fq_nmod and native GF(2^128) - IMPLEMENTATION */
gf2128_t fq_nmod_to_gf2128(const fq_nmod_t elem, const fq_nmod_ctx_t ctx) {
   if (fq_nmod_is_zero(elem, ctx)) {
       return gf2128_zero();
   }
   
   uint64_t poly_low = 0, poly_high = 0;
   
   for (int i = 0; i < 64; i++) {
       if (nmod_poly_get_coeff_ui(elem, i)) {
           poly_low |= (1ULL << i);
       }
   }
   for (int i = 64; i < 128; i++) {
       if (nmod_poly_get_coeff_ui(elem, i)) {
           poly_high |= (1ULL << (i - 64));
       }
   }
   
   return gf2128_create(poly_low, poly_high);
}

void gf2128_to_fq_nmod(fq_nmod_t res, const gf2128_t *elem, const fq_nmod_ctx_t ctx) {
   fq_nmod_zero(res, ctx);
   
   if (gf2128_is_zero(elem)) {
       return;
   }
   
   for (int i = 0; i < 64; i++) {
       if ((elem->low >> i) & 1) {
           nmod_poly_set_coeff_ui(res, i, 1);
       }
   }
   for (int i = 0; i < 64; i++) {
       if ((elem->high >> i) & 1) {
           nmod_poly_set_coeff_ui(res, i + 64, 1);
       }
   }
}

/* Initialize GF(2^8) conversion tables */
void init_gf28_conversion(const fq_nmod_ctx_t ctx) {
   if (g_gf28_conversion && g_gf28_conversion->initialized) {
       return;
   }
   
   if (!g_gf28_conversion) {
       g_gf28_conversion = (gf28_conversion_t *)calloc(1, sizeof(gf28_conversion_t));
   }
   
   if (!g_gf28_tables) {
       init_gf28_tables(0x1D);
   }
   
   uint64_t flint_poly = extract_irred_poly(ctx);
   
   if (flint_poly == 0x11D) {
       for (int i = 0; i < 256; i++) {
           g_gf28_conversion->flint_to_gf28[i] = i;
           g_gf28_conversion->gf28_to_flint[i] = i;
       }
   } else {
       g_gf28_conversion->flint_to_gf28[0] = 0;
       g_gf28_conversion->gf28_to_flint[0] = 0;
       
       fq_nmod_t gen, elem;
       fq_nmod_init(gen, ctx);
       fq_nmod_init(elem, ctx);
       
       fq_nmod_gen(gen, ctx);
       
       uint8_t flint_log[256];
       uint8_t flint_exp[256];
       memset(flint_log, 0xFF, 256);
       
       fq_nmod_one(elem, ctx);
       for (int i = 0; i < 255; i++) {
           uint8_t poly_val = 0;
           for (int j = 0; j < 8; j++) {
               if (nmod_poly_get_coeff_ui(elem, j)) {
                   poly_val |= (1 << j);
               }
           }
           
           flint_exp[i] = poly_val;
           flint_log[poly_val] = i;
           
           fq_nmod_mul(elem, elem, gen, ctx);
       }
       
       for (int i = 0; i < 255; i++) {
           uint8_t flint_elem = flint_exp[i];
           uint8_t our_elem = g_gf28_tables->exp_table[i];
           
           g_gf28_conversion->flint_to_gf28[flint_elem] = our_elem;
           g_gf28_conversion->gf28_to_flint[our_elem] = flint_elem;
       }
       
       fq_nmod_clear(gen, ctx);
       fq_nmod_clear(elem, ctx);
   }
   
   g_gf28_conversion->initialized = 1;
}

/* Initialize GF(2^16) conversion tables */
void init_gf216_conversion(const fq_nmod_ctx_t ctx) {
    if (g_gf216_conversion && g_gf216_conversion->initialized) {
        return;
    }
    
    if (!g_gf216_conversion) {
        g_gf216_conversion = (gf216_conversion_t *)calloc(1, sizeof(gf216_conversion_t));
        g_gf216_conversion->flint_to_gf216 = (uint16_t *)calloc(65536, sizeof(uint16_t));
        g_gf216_conversion->gf216_to_flint = (uint16_t *)calloc(65536, sizeof(uint16_t));
    }
    
    init_gf216_standard();
    
    // Check if polynomials match
    uint64_t flint_poly = extract_irred_poly(ctx);
    
    if (flint_poly == 0x1002D) {
        // Same polynomial, use identity mapping
        for (int i = 0; i < 65536; i++) {
            g_gf216_conversion->flint_to_gf216[i] = i;
            g_gf216_conversion->gf216_to_flint[i] = i;
        }
    } else {
        // Different polynomials, build conversion table
        fq_nmod_t elem, gen;
        fq_nmod_init(elem, ctx);
        fq_nmod_init(gen, ctx);
        
        fq_nmod_gen(gen, ctx);
        
        // Map zero
        g_gf216_conversion->flint_to_gf216[0] = 0;
        g_gf216_conversion->gf216_to_flint[0] = 0;
        
        // Build logarithm tables for both representations
        uint16_t flint_log[65536];
        uint16_t flint_exp[65536];
        memset(flint_log, 0xFF, sizeof(flint_log));
        
        fq_nmod_one(elem, ctx);
        for (int i = 0; i < 65535; i++) {
            uint16_t flint_rep = 0;
            for (int j = 0; j < 16; j++) {
                if (nmod_poly_get_coeff_ui(elem, j)) {
                    flint_rep |= (1 << j);
                }
            }
            
            flint_exp[i] = flint_rep;
            flint_log[flint_rep] = i;
            
            fq_nmod_mul(elem, elem, gen, ctx);
        }
        
        // Map between representations using logarithms
        for (int i = 0; i < 65535; i++) {
            uint16_t flint_elem = flint_exp[i];
            uint16_t our_elem = g_gf216_tables->exp_table[i];
            
            g_gf216_conversion->flint_to_gf216[flint_elem] = our_elem;
            g_gf216_conversion->gf216_to_flint[our_elem] = flint_elem;
        }
        
        fq_nmod_clear(elem, ctx);
        fq_nmod_clear(gen, ctx);
    }
    
    g_gf216_conversion->initialized = 1;
}

void cleanup_gf28_conversion(void) {
   if (g_gf28_conversion) {
       free(g_gf28_conversion);
       g_gf28_conversion = NULL;
   }
}

void cleanup_gf216_conversion(void) {
   if (g_gf216_conversion) {
       if (g_gf216_conversion->flint_to_gf216) {
           free(g_gf216_conversion->flint_to_gf216);
       }
       if (g_gf216_conversion->gf216_to_flint) {
           free(g_gf216_conversion->gf216_to_flint);
       }
       free(g_gf216_conversion);
       g_gf216_conversion = NULL;
   }
}

uint8_t fq_nmod_to_gf28_elem(const fq_nmod_t elem, const fq_nmod_ctx_t ctx) {
   if (!g_gf28_conversion || !g_gf28_conversion->initialized) {
       init_gf28_conversion(ctx);
   }
   
   if (fq_nmod_is_zero(elem, ctx)) {
       return 0;
   }
   
   uint8_t flint_rep = 0;
   for (int j = 0; j < 8; j++) {
       if (nmod_poly_get_coeff_ui(elem, j)) {
           flint_rep |= (1 << j);
       }
   }
   
   return g_gf28_conversion->flint_to_gf28[flint_rep];
}

void gf28_elem_to_fq_nmod(fq_nmod_t res, uint8_t elem, const fq_nmod_ctx_t ctx) {
   if (!g_gf28_conversion || !g_gf28_conversion->initialized) {
       init_gf28_conversion(ctx);
   }
   
   fq_nmod_zero(res, ctx);
   
   if (elem == 0) {
       return;
   }
   
   uint8_t flint_rep = g_gf28_conversion->gf28_to_flint[elem];
   
   for (int j = 0; j < 8; j++) {
       if (flint_rep & (1 << j)) {
           nmod_poly_set_coeff_ui(res, j, 1);
       }
   }
}

uint16_t fq_nmod_to_gf216_elem(const fq_nmod_t elem, const fq_nmod_ctx_t ctx) {
   if (!g_gf216_conversion || !g_gf216_conversion->initialized) {
       init_gf216_conversion(ctx);
   }
   
   if (fq_nmod_is_zero(elem, ctx)) {
       return 0;
   }
   
   uint16_t flint_rep = 0;
   for (int j = 0; j < 16; j++) {
       if (nmod_poly_get_coeff_ui(elem, j)) {
           flint_rep |= (1 << j);
       }
   }
   
   return g_gf216_conversion->flint_to_gf216[flint_rep];
}

void gf216_elem_to_fq_nmod(fq_nmod_t res, uint16_t elem, const fq_nmod_ctx_t ctx) {
   if (!g_gf216_conversion || !g_gf216_conversion->initialized) {
       init_gf216_conversion(ctx);
   }
   
   fq_nmod_zero(res, ctx);
   
   if (elem == 0) {
       return;
   }
   
   uint16_t flint_rep = g_gf216_conversion->gf216_to_flint[elem];
   
   for (int j = 0; j < 16; j++) {
       if (flint_rep & (1 << j)) {
           nmod_poly_set_coeff_ui(res, j, 1);
       }
   }
}

/* Initialize GF(2^32) conversion tables */
void init_gf232_conversion(const fq_nmod_ctx_t ctx) {
    if (g_gf232_conversion && g_gf232_conversion->initialized) {
        return;
    }
    
    if (!g_gf232_conversion) {
        g_gf232_conversion = (gf232_conversion_t *)calloc(1, sizeof(gf232_conversion_t));
    }
    
    uint64_t flint_poly = extract_irred_poly(ctx);
    g_gf232_conversion->flint_poly = (uint32_t)(flint_poly & 0xFFFFFFFF);
    g_gf232_conversion->our_poly = 0x8299; // x^32 + x^15 + x^9 + x^7 + x^4 + x^3 + 1
    
    if ((flint_poly >> 32) == 1 && (flint_poly & 0xFFFFFFFF) == g_gf232_conversion->our_poly) {
        g_gf232_conversion->initialized = 1;
        return;
    }
    
    g_gf232_conversion->initialized = 1;
}


/* Enhanced conversion functions for GF(2^32) */
gf232_t fq_nmod_to_gf232(const fq_nmod_t elem, const fq_nmod_ctx_t ctx) {
    if (!g_gf232_conversion || !g_gf232_conversion->initialized) {
        init_gf232_conversion(ctx);
    }
    
    if (fq_nmod_is_zero(elem, ctx)) {
        return gf232_zero();
    }
    
    // If polynomials are the same, direct conversion
    if (g_gf232_conversion->flint_poly == g_gf232_conversion->our_poly) {
        uint32_t poly_val = 0;
        for (int i = 0; i < 32; i++) {
            if (nmod_poly_get_coeff_ui(elem, i)) {
                poly_val |= (1UL << i);
            }
        }
        return gf232_create(poly_val);
    }
    
    // Otherwise, we need to convert between different field representations
    // This is more complex for GF(2^32), so we use FLINT for conversion
    fq_nmod_ctx_t our_ctx;
    nmod_poly_t our_mod;
    nmod_poly_init(our_mod, 2);
    
    // Set modulus: x^32 + x^15 + x^9 + x^7 + x^4 + x^3 + 1
    nmod_poly_set_coeff_ui(our_mod, 0, 1);
    nmod_poly_set_coeff_ui(our_mod, 3, 1);
    nmod_poly_set_coeff_ui(our_mod, 4, 1);
    nmod_poly_set_coeff_ui(our_mod, 7, 1);
    nmod_poly_set_coeff_ui(our_mod, 9, 1);
    nmod_poly_set_coeff_ui(our_mod, 15, 1);
    nmod_poly_set_coeff_ui(our_mod, 32, 1);
    
    fq_nmod_ctx_init_modulus(our_ctx, our_mod, "t");
    
    // Convert element representation
    fq_nmod_t converted;
    fq_nmod_init(converted, our_ctx);
    
    // Copy polynomial coefficients
    for (int i = 0; i < 32; i++) {
        if (nmod_poly_get_coeff_ui(elem, i)) {
            nmod_poly_set_coeff_ui(converted, i, 1);
        }
    }
    
    // Reduce in our field
    fq_nmod_reduce(converted, our_ctx);
    
    // Extract result
    uint32_t result = 0;
    for (int i = 0; i < 32; i++) {
        if (nmod_poly_get_coeff_ui(converted, i)) {
            result |= (1UL << i);
        }
    }
    
    fq_nmod_clear(converted, our_ctx);
    fq_nmod_ctx_clear(our_ctx);
    nmod_poly_clear(our_mod);
    
    return gf232_create(result);
}

void gf232_to_fq_nmod(fq_nmod_t res, const gf232_t *elem, const fq_nmod_ctx_t ctx) {
    if (!g_gf232_conversion || !g_gf232_conversion->initialized) {
        init_gf232_conversion(ctx);
    }
    
    fq_nmod_zero(res, ctx);
    
    if (gf232_is_zero(elem)) {
        return;
    }
    
    if (g_gf232_conversion->flint_poly == g_gf232_conversion->our_poly) {
        for (int i = 0; i < 32; i++) {
            if ((elem->value >> i) & 1) {
                nmod_poly_set_coeff_ui(res, i, 1);
            }
        }
        return;
    }
    
    // Different polynomials - convert
    fq_nmod_ctx_t our_ctx;
    nmod_poly_t our_mod;
    nmod_poly_init(our_mod, 2);
    
    // Set modulus: x^32 + x^15 + x^9 + x^7 + x^4 + x^3 + 1
    nmod_poly_set_coeff_ui(our_mod, 0, 1);
    nmod_poly_set_coeff_ui(our_mod, 3, 1);
    nmod_poly_set_coeff_ui(our_mod, 4, 1);
    nmod_poly_set_coeff_ui(our_mod, 7, 1);
    nmod_poly_set_coeff_ui(our_mod, 9, 1);
    nmod_poly_set_coeff_ui(our_mod, 15, 1);
    nmod_poly_set_coeff_ui(our_mod, 32, 1);
    
    fq_nmod_ctx_init_modulus(our_ctx, our_mod, "t");
    
    fq_nmod_t temp;
    fq_nmod_init(temp, our_ctx);
    
    for (int i = 0; i < 32; i++) {
        if ((elem->value >> i) & 1) {
            nmod_poly_set_coeff_ui(temp, i, 1);
        }
    }
    
    for (int i = 0; i < 32; i++) {
        if ((elem->value >> i) & 1) {
            nmod_poly_set_coeff_ui(res, i, 1);
        }
    }
    fq_nmod_reduce(res, ctx);
    
    fq_nmod_clear(temp, our_ctx);
    fq_nmod_ctx_clear(our_ctx);
    nmod_poly_clear(our_mod);
}

/* Initialize GF(2^64) conversion with full polynomial matching */
void init_gf264_conversion(const fq_nmod_ctx_t ctx) {
    if (g_gf264_conversion && g_gf264_conversion->initialized) {
        return;
    }
    
    if (!g_gf264_conversion) {
        g_gf264_conversion = (gf264_conversion_t *)calloc(1, sizeof(gf264_conversion_t));
    }
    
    extract_gf264_poly(ctx, &g_gf264_conversion->flint_poly_low, 
                            &g_gf264_conversion->flint_poly_high);
    
    /* Check if polynomials match: high should be 0x2, low should be 0x47F43CB7 */
    if (g_gf264_conversion->flint_poly_high == 0x2 &&
        g_gf264_conversion->flint_poly_low == 0x47F43CB7) {
        g_gf264_conversion->initialized = 1;
        return;
    }
    
    /* Different polynomials - need conversion tables */
    g_gf264_conversion->flint_to_gf264_low = (uint64_t *)calloc(65536, sizeof(uint64_t));
    g_gf264_conversion->gf264_to_flint_low = (uint64_t *)calloc(65536, sizeof(uint64_t));
    
    fq_nmod_t gen_flint, elem;
    fq_nmod_init(gen_flint, ctx);
    fq_nmod_init(elem, ctx);
    
    fq_nmod_gen(gen_flint, ctx);
    
    fq_nmod_one(elem, ctx);
    for (int i = 0; i < 65536; i++) {
        uint64_t flint_val = 0;
        for (int j = 0; j < 64; j++) {
            if (nmod_poly_get_coeff_ui(elem, j)) {
                flint_val |= (1ULL << j);
            }
        }
        
        g_gf264_conversion->flint_to_gf264_low[i] = flint_val;
        if (flint_val < 65536) {
            g_gf264_conversion->gf264_to_flint_low[flint_val] = i;
        }
        
        fq_nmod_mul(elem, elem, gen_flint, ctx);
    }
    
    fq_nmod_clear(gen_flint, ctx);
    fq_nmod_clear(elem, ctx);
    
    g_gf264_conversion->initialized = 1;
}

/* Enhanced conversion functions for GF(2^64) */
gf264_t fq_nmod_to_gf264(const fq_nmod_t elem, const fq_nmod_ctx_t ctx) {
    if (!g_gf264_conversion || !g_gf264_conversion->initialized) {
        init_gf264_conversion(ctx);
    }
    
    if (fq_nmod_is_zero(elem, ctx)) {
        return gf264_zero();
    }
    
    // Check if same polynomial
    if (g_gf264_conversion->flint_poly_high == 0x2 && 
        g_gf264_conversion->flint_poly_low == 0x47F43CB7) {
        // Direct conversion
        uint64_t poly_val = 0;
        for (int i = 0; i < 64; i++) {
            if (nmod_poly_get_coeff_ui(elem, i)) {
                poly_val |= (1ULL << i);
            }
        }
        return gf264_create(poly_val);
    }
    
    // Different polynomials - need field isomorphism
    // For now, we'll use FLINT to handle the conversion
    
    // Create our field context
    fq_nmod_ctx_t our_ctx;
    nmod_poly_t our_mod;
    nmod_poly_init(our_mod, 2);
    
    nmod_poly_set_coeff_ui(our_mod, 0, 1);
    nmod_poly_set_coeff_ui(our_mod, 1, 1);
    nmod_poly_set_coeff_ui(our_mod, 2, 1);
    nmod_poly_set_coeff_ui(our_mod, 4, 1);
    nmod_poly_set_coeff_ui(our_mod, 5, 1);
    nmod_poly_set_coeff_ui(our_mod, 7, 1);
    nmod_poly_set_coeff_ui(our_mod, 10, 1);
    nmod_poly_set_coeff_ui(our_mod, 11, 1);
    nmod_poly_set_coeff_ui(our_mod, 12, 1);
    nmod_poly_set_coeff_ui(our_mod, 13, 1);
    nmod_poly_set_coeff_ui(our_mod, 18, 1);
    nmod_poly_set_coeff_ui(our_mod, 20, 1);
    nmod_poly_set_coeff_ui(our_mod, 21, 1);
    nmod_poly_set_coeff_ui(our_mod, 22, 1);
    nmod_poly_set_coeff_ui(our_mod, 23, 1);
    nmod_poly_set_coeff_ui(our_mod, 24, 1);
    nmod_poly_set_coeff_ui(our_mod, 25, 1);
    nmod_poly_set_coeff_ui(our_mod, 26, 1);
    nmod_poly_set_coeff_ui(our_mod, 30, 1);
    nmod_poly_set_coeff_ui(our_mod, 33, 1);
    nmod_poly_set_coeff_ui(our_mod, 64, 1);
    fq_nmod_ctx_init_modulus(our_ctx, our_mod, "t");
    
    // Convert element by polynomial evaluation
    fq_nmod_t result;
    fq_nmod_init(result, our_ctx);
    
    // Copy polynomial representation
    for (int i = 0; i < 64; i++) {
        if (nmod_poly_get_coeff_ui(elem, i)) {
            nmod_poly_set_coeff_ui(result, i, 1);
        }
    }
    
    // Reduce in our field
    fq_nmod_reduce(result, our_ctx);
    
    // Extract result
    uint64_t val = 0;
    for (int i = 0; i < 64; i++) {
        if (nmod_poly_get_coeff_ui(result, i)) {
            val |= (1ULL << i);
        }
    }
    
    fq_nmod_clear(result, our_ctx);
    fq_nmod_ctx_clear(our_ctx);
    nmod_poly_clear(our_mod);
    
    return gf264_create(val);
}

void gf264_to_fq_nmod(fq_nmod_t res, const gf264_t *elem, const fq_nmod_ctx_t ctx) {
    if (!g_gf264_conversion || !g_gf264_conversion->initialized) {
        init_gf264_conversion(ctx);
    }
    
    fq_nmod_zero(res, ctx);
    
    if (gf264_is_zero(elem)) {
        return;
    }
    
    if (g_gf264_conversion->flint_poly_high == 0x2 && 
        g_gf264_conversion->flint_poly_low == 0x47F43CB7) {
        for (int i = 0; i < 64; i++) {
            if ((elem->value >> i) & 1) {
                nmod_poly_set_coeff_ui(res, i, 1);
            }
        }
        return;
    }
    
    // Different polynomials - convert
    fq_nmod_ctx_t our_ctx;
    nmod_poly_t our_mod;
    nmod_poly_init(our_mod, 2);
    
    // Set modulus with all terms
    nmod_poly_set_coeff_ui(our_mod, 0, 1);
    nmod_poly_set_coeff_ui(our_mod, 1, 1);
    nmod_poly_set_coeff_ui(our_mod, 2, 1);
    nmod_poly_set_coeff_ui(our_mod, 4, 1);
    nmod_poly_set_coeff_ui(our_mod, 5, 1);
    nmod_poly_set_coeff_ui(our_mod, 7, 1);
    nmod_poly_set_coeff_ui(our_mod, 10, 1);
    nmod_poly_set_coeff_ui(our_mod, 11, 1);
    nmod_poly_set_coeff_ui(our_mod, 12, 1);
    nmod_poly_set_coeff_ui(our_mod, 13, 1);
    nmod_poly_set_coeff_ui(our_mod, 18, 1);
    nmod_poly_set_coeff_ui(our_mod, 20, 1);
    nmod_poly_set_coeff_ui(our_mod, 21, 1);
    nmod_poly_set_coeff_ui(our_mod, 22, 1);
    nmod_poly_set_coeff_ui(our_mod, 23, 1);
    nmod_poly_set_coeff_ui(our_mod, 24, 1);
    nmod_poly_set_coeff_ui(our_mod, 25, 1);
    nmod_poly_set_coeff_ui(our_mod, 26, 1);
    nmod_poly_set_coeff_ui(our_mod, 30, 1);
    nmod_poly_set_coeff_ui(our_mod, 33, 1);
    nmod_poly_set_coeff_ui(our_mod, 64, 1);
    
    fq_nmod_ctx_init_modulus(our_ctx, our_mod, "t");
    
    fq_nmod_t temp;
    fq_nmod_init(temp, our_ctx);
    
    for (int i = 0; i < 64; i++) {
        if ((elem->value >> i) & 1) {
            nmod_poly_set_coeff_ui(temp, i, 1);
        }
    }
    
    for (int i = 0; i < 64; i++) {
        if ((elem->value >> i) & 1) {
            nmod_poly_set_coeff_ui(res, i, 1);
        }
    }
    
    fq_nmod_reduce(res, ctx);
    
    fq_nmod_clear(temp, our_ctx);
    fq_nmod_ctx_clear(our_ctx);
    nmod_poly_clear(our_mod);
}

void cleanup_gf232_conversion(void) {
    if (g_gf232_conversion) {
        free(g_gf232_conversion);
        g_gf232_conversion = NULL;
    }
}

void cleanup_gf264_conversion(void) {
    if (g_gf264_conversion) {
        free(g_gf264_conversion);
        g_gf264_conversion = NULL;
    }
}