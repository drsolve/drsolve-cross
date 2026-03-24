/* gf2n_poly.c - Implementation of polynomial operations over GF(2^n) fields */
#include "gf2n_poly.h"

/* ============================================================================
   GF(2^8) POLYNOMIAL OPERATIONS IMPLEMENTATION
   ============================================================================ */

/* Basic operations */
void gf28_poly_init(gf28_poly_t poly) {
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

void gf28_poly_clear(gf28_poly_t poly) {
    if (poly->coeffs) {
        free(poly->coeffs);
        poly->coeffs = NULL;
    }
    poly->length = 0;
    poly->alloc = 0;
}

void gf28_poly_fit_length(gf28_poly_t poly, slong len) {
    if (len > poly->alloc) {
        slong new_alloc = FLINT_MAX(len, poly->alloc * 2);
        poly->coeffs = (uint8_t *)realloc(poly->coeffs, new_alloc * sizeof(uint8_t));
        memset(poly->coeffs + poly->alloc, 0, (new_alloc - poly->alloc) * sizeof(uint8_t));
        poly->alloc = new_alloc;
    }
}

void gf28_poly_normalise(gf28_poly_t poly) {
    while (poly->length > 0 && poly->coeffs[poly->length - 1] == 0) {
        poly->length--;
    }
}

inline void gf28_poly_zero(gf28_poly_t poly) {
    poly->length = 0;
}

inline int gf28_poly_is_zero(const gf28_poly_t poly) {
    return poly->length == 0;
}

inline slong gf28_poly_degree(const gf28_poly_t poly) {
    return poly->length - 1;
}

void gf28_poly_set(gf28_poly_t res, const gf28_poly_t poly) {
    if (res == poly) return;
    gf28_poly_fit_length(res, poly->length);
    memcpy(res->coeffs, poly->coeffs, poly->length * sizeof(uint8_t));
    res->length = poly->length;
}

uint8_t gf28_poly_get_coeff(const gf28_poly_t poly, slong i) {
    return (i < poly->length) ? poly->coeffs[i] : 0;
}

void gf28_poly_set_coeff(gf28_poly_t poly, slong i, uint8_t c) {
    if (i >= poly->alloc) {
        slong new_alloc = FLINT_MAX(i + 1, 2 * poly->alloc);
        poly->coeffs = (uint8_t *)flint_realloc(poly->coeffs, new_alloc * sizeof(uint8_t));
        poly->alloc = new_alloc;
    }

    if (i >= poly->length) {
        for (slong j = poly->length; j < i; j++) {
            poly->coeffs[j] = 0; // GF(2^8) zero element is 0
        }
        poly->length = i + 1;
    }

    poly->coeffs[i] = c;

    if (c == 0 && i == poly->length - 1) {
        while (poly->length > 0 && poly->coeffs[poly->length - 1] == 0) {
            poly->length--;
        }
    }
}

void gf28_poly_add(gf28_poly_t res, const gf28_poly_t a, const gf28_poly_t b) {
    slong max_len = FLINT_MAX(a->length, b->length);
    if (max_len == 0) {
        gf28_poly_zero(res);
        return;
    }
    
    gf28_poly_fit_length(res, max_len);
    slong min_len = FLINT_MIN(a->length, b->length);
    
    for (slong i = 0; i < min_len; i++) {
        res->coeffs[i] = gf28_add(a->coeffs[i], b->coeffs[i]);
    }
    
    if (a->length > b->length) {
        memcpy(res->coeffs + min_len, a->coeffs + min_len, 
               (a->length - min_len) * sizeof(uint8_t));
        res->length = a->length;
    } else if (b->length > a->length) {
        memcpy(res->coeffs + min_len, b->coeffs + min_len, 
               (b->length - min_len) * sizeof(uint8_t));
        res->length = b->length;
    } else {
        res->length = min_len;
    }
    
    gf28_poly_normalise(res);
}

void gf28_poly_scalar_mul(gf28_poly_t res, const gf28_poly_t poly, uint8_t c) {
    if (c == 0) {
        gf28_poly_zero(res);
        return;
    }
    
    if (c == 1) {
        gf28_poly_set(res, poly);
        return;
    }
    
    gf28_poly_fit_length(res, poly->length);
    res->length = poly->length;
    
    for (slong i = 0; i < poly->length; i++) {
        res->coeffs[i] = gf28_mul(poly->coeffs[i], c);
    }
    
    gf28_poly_normalise(res);
}

void gf28_poly_shift_left(gf28_poly_t res, const gf28_poly_t poly, slong n) {
    if (n == 0) {
        gf28_poly_set(res, poly);
        return;
    }
    
    if (gf28_poly_is_zero(poly)) {
        gf28_poly_zero(res);
        return;
    }
    
    slong new_len = poly->length + n;
    gf28_poly_fit_length(res, new_len);
    
    memmove(res->coeffs + n, poly->coeffs, poly->length * sizeof(uint8_t));
    memset(res->coeffs, 0, n * sizeof(uint8_t));
    res->length = new_len;
}

void gf28_poly_mul_schoolbook(gf28_poly_t res, const gf28_poly_t a, const gf28_poly_t b) {
    if (gf28_poly_is_zero(a) || gf28_poly_is_zero(b)) {
        gf28_poly_zero(res);
        return;
    }
    
    slong rlen = a->length + b->length - 1;
    uint8_t *temp = (uint8_t *)calloc(rlen, sizeof(uint8_t));
    
    for (slong i = 0; i < a->length; i++) {
        if (a->coeffs[i] == 0) continue;
        for (slong j = 0; j < b->length; j++) {
            if (b->coeffs[j] == 0) continue;
            temp[i + j] = gf28_add(temp[i + j], gf28_mul(a->coeffs[i], b->coeffs[j]));
        }
    }
    
    gf28_poly_fit_length(res, rlen);
    memcpy(res->coeffs, temp, rlen * sizeof(uint8_t));
    res->length = rlen;
    gf28_poly_normalise(res);
    
    free(temp);
}

/* Fast polynomial multiplication for GF(2^8) using Karatsuba algorithm */
void gf28_poly_mul_karatsuba(gf28_poly_t res, const gf28_poly_t a, const gf28_poly_t b) {
    if (gf28_poly_is_zero(a) || gf28_poly_is_zero(b)) {
        gf28_poly_zero(res);
        return;
    }
    
    slong alen = a->length;
    slong blen = b->length;
    
    /* For small polynomials, use schoolbook multiplication */
    if (alen < 32 || blen < 32) {
        slong rlen = alen + blen - 1;
        uint8_t *temp = (uint8_t *)calloc(rlen, sizeof(uint8_t));
        
        for (slong i = 0; i < alen; i++) {
            if (a->coeffs[i] == 0) continue;
            const uint8_t* row = gf28_get_scalar_row(a->coeffs[i]);
            for (slong j = 0; j < blen; j++) {
                temp[i + j] ^= row[b->coeffs[j]];
            }
        }
        
        gf28_poly_fit_length(res, rlen);
        memcpy(res->coeffs, temp, rlen * sizeof(uint8_t));
        res->length = rlen;
        gf28_poly_normalise(res);
        free(temp);
        return;
    }
    
    /* Karatsuba algorithm */
    slong split = FLINT_MAX(alen, blen) / 2;
    
    gf28_poly_t a0, a1, b0, b1;
    gf28_poly_t a0b0, a1b1, amid, bmid, mid;
    
    gf28_poly_init(a0);
    gf28_poly_init(a1);
    gf28_poly_init(b0);
    gf28_poly_init(b1);
    gf28_poly_init(a0b0);
    gf28_poly_init(a1b1);
    gf28_poly_init(amid);
    gf28_poly_init(bmid);
    gf28_poly_init(mid);
    
    /* Split polynomials: a = a0 + a1*x^split */
    gf28_poly_fit_length(a0, FLINT_MIN(split, alen));
    a0->length = FLINT_MIN(split, alen);
    memcpy(a0->coeffs, a->coeffs, a0->length * sizeof(uint8_t));
    gf28_poly_normalise(a0);
    
    if (alen > split) {
        gf28_poly_fit_length(a1, alen - split);
        a1->length = alen - split;
        memcpy(a1->coeffs, a->coeffs + split, a1->length * sizeof(uint8_t));
        gf28_poly_normalise(a1);
    }
    
    /* Split b similarly */
    gf28_poly_fit_length(b0, FLINT_MIN(split, blen));
    b0->length = FLINT_MIN(split, blen);
    memcpy(b0->coeffs, b->coeffs, b0->length * sizeof(uint8_t));
    gf28_poly_normalise(b0);
    
    if (blen > split) {
        gf28_poly_fit_length(b1, blen - split);
        b1->length = blen - split;
        memcpy(b1->coeffs, b->coeffs + split, b1->length * sizeof(uint8_t));
        gf28_poly_normalise(b1);
    }
    
    /* Compute three products recursively */
    gf28_poly_mul_karatsuba(a0b0, a0, b0);
    gf28_poly_mul_karatsuba(a1b1, a1, b1);
    
    gf28_poly_add(amid, a0, a1);
    gf28_poly_add(bmid, b0, b1);
    gf28_poly_mul_karatsuba(mid, amid, bmid);
    
    /* Combine results: res = a0b0 + ((a0+a1)(b0+b1) - a0b0 - a1b1)*x^split + a1b1*x^(2*split) */
    gf28_poly_t temp1, temp2;
    gf28_poly_init(temp1);
    gf28_poly_init(temp2);
    
    /* temp1 = mid - a0b0 - a1b1 */
    gf28_poly_add(temp1, mid, a0b0);
    gf28_poly_add(temp1, temp1, a1b1);
    
    /* Build result */
    slong rlen = alen + blen - 1;
    gf28_poly_fit_length(res, rlen);
    memset(res->coeffs, 0, rlen * sizeof(uint8_t));
    
    /* Add a0b0 */
    for (slong i = 0; i < a0b0->length; i++) {
        res->coeffs[i] ^= a0b0->coeffs[i];
    }
    
    /* Add temp1 * x^split */
    for (slong i = 0; i < temp1->length; i++) {
        res->coeffs[i + split] ^= temp1->coeffs[i];
    }
    
    /* Add a1b1 * x^(2*split) */
    for (slong i = 0; i < a1b1->length; i++) {
        res->coeffs[i + 2*split] ^= a1b1->coeffs[i];
    }
    
    res->length = rlen;
    gf28_poly_normalise(res);
    
    /* Cleanup */
    gf28_poly_clear(a0);
    gf28_poly_clear(a1);
    gf28_poly_clear(b0);
    gf28_poly_clear(b1);
    gf28_poly_clear(a0b0);
    gf28_poly_clear(a1b1);
    gf28_poly_clear(amid);
    gf28_poly_clear(bmid);
    gf28_poly_clear(mid);
    gf28_poly_clear(temp1);
    gf28_poly_clear(temp2);
}

/* Replace the original gf28_poly_mul with the fast version */
void gf28_poly_mul(gf28_poly_t res, const gf28_poly_t a, const gf28_poly_t b) {
    gf28_poly_mul_schoolbook(res, a, b); // gf28_poly_mul_karatsuba gf28_poly_mul_schoolbook
}

/* ============================================================================
   GF(2^16) POLYNOMIAL OPERATIONS IMPLEMENTATION
   ============================================================================ */

/* Basic operations */
void gf216_poly_init(gf216_poly_t poly) {
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

void gf216_poly_clear(gf216_poly_t poly) {
    if (poly->coeffs) {
        free(poly->coeffs);
        poly->coeffs = NULL;
    }
    poly->length = 0;
    poly->alloc = 0;
}

void gf216_poly_fit_length(gf216_poly_t poly, slong len) {
    if (len > poly->alloc) {
        slong new_alloc = FLINT_MAX(len, poly->alloc * 2);
        poly->coeffs = (uint16_t *)realloc(poly->coeffs, new_alloc * sizeof(uint16_t));
        memset(poly->coeffs + poly->alloc, 0, (new_alloc - poly->alloc) * sizeof(uint16_t));
        poly->alloc = new_alloc;
    }
}

void gf216_poly_normalise(gf216_poly_t poly) {
    while (poly->length > 0 && poly->coeffs[poly->length - 1] == 0) {
        poly->length--;
    }
}

inline void gf216_poly_zero(gf216_poly_t poly) {
    poly->length = 0;
}

inline int gf216_poly_is_zero(const gf216_poly_t poly) {
    return poly->length == 0;
}

inline slong gf216_poly_degree(const gf216_poly_t poly) {
    return poly->length - 1;
}

void gf216_poly_set(gf216_poly_t res, const gf216_poly_t poly) {
    if (res == poly) return;
    gf216_poly_fit_length(res, poly->length);
    memcpy(res->coeffs, poly->coeffs, poly->length * sizeof(uint16_t));
    res->length = poly->length;
}

void gf216_poly_set_coeff(gf216_poly_t poly, slong i, uint16_t c) {
    if (i >= poly->alloc) {
        slong new_alloc = FLINT_MAX(i + 1, 2 * poly->alloc);
        poly->coeffs = (uint16_t *)flint_realloc(poly->coeffs, new_alloc * sizeof(uint16_t));
        poly->alloc = new_alloc;
    }

    if (i >= poly->length) {
        for (slong j = poly->length; j < i; j++) {
            poly->coeffs[j] = 0; // GF(2^16) zero element is 0
        }
        poly->length = i + 1;
    }

    poly->coeffs[i] = c;

    if (c == 0 && i == poly->length - 1) {
        while (poly->length > 0 && poly->coeffs[poly->length - 1] == 0) {
            poly->length--;
        }
    }
}

uint16_t gf216_poly_get_coeff(const gf216_poly_t poly, slong i) {
    return (i < poly->length) ? poly->coeffs[i] : 0;
}

void gf216_poly_add(gf216_poly_t res, const gf216_poly_t a, const gf216_poly_t b) {
    slong max_len = FLINT_MAX(a->length, b->length);
    if (max_len == 0) {
        gf216_poly_zero(res);
        return;
    }
    
    gf216_poly_fit_length(res, max_len);
    slong min_len = FLINT_MIN(a->length, b->length);
    
    for (slong i = 0; i < min_len; i++) {
        res->coeffs[i] = a->coeffs[i] ^ b->coeffs[i];
    }
    
    if (a->length > b->length) {
        memcpy(res->coeffs + min_len, a->coeffs + min_len, 
               (a->length - min_len) * sizeof(uint16_t));
        res->length = a->length;
    } else if (b->length > a->length) {
        memcpy(res->coeffs + min_len, b->coeffs + min_len, 
               (b->length - min_len) * sizeof(uint16_t));
        res->length = b->length;
    } else {
        res->length = min_len;
    }
    
    gf216_poly_normalise(res);
}

void gf216_poly_scalar_mul(gf216_poly_t res, const gf216_poly_t poly, uint16_t c) {
    if (c == 0) {
        gf216_poly_zero(res);
        return;
    }
    
    if (c == 1) {
        gf216_poly_set(res, poly);
        return;
    }
    
    gf216_poly_fit_length(res, poly->length);
    res->length = poly->length;
    
    for (slong i = 0; i < poly->length; i++) {
        res->coeffs[i] = gf216_mul(poly->coeffs[i], c);
    }
    
    gf216_poly_normalise(res);
}

void gf216_poly_shift_left(gf216_poly_t res, const gf216_poly_t poly, slong n) {
    if (n == 0) {
        gf216_poly_set(res, poly);
        return;
    }
    
    if (gf216_poly_is_zero(poly)) {
        gf216_poly_zero(res);
        return;
    }
    
    slong new_len = poly->length + n;
    gf216_poly_fit_length(res, new_len);
    
    memmove(res->coeffs + n, poly->coeffs, poly->length * sizeof(uint16_t));
    memset(res->coeffs, 0, n * sizeof(uint16_t));
    res->length = new_len;
}

void gf216_poly_mul(gf216_poly_t res, const gf216_poly_t a, const gf216_poly_t b) {
    if (gf216_poly_is_zero(a) || gf216_poly_is_zero(b)) {
        gf216_poly_zero(res);
        return;
    }
    
    slong rlen = a->length + b->length - 1;
    uint16_t *temp = (uint16_t *)calloc(rlen, sizeof(uint16_t));
    
    for (slong i = 0; i < a->length; i++) {
        if (a->coeffs[i] == 0) continue;
        for (slong j = 0; j < b->length; j++) {
            if (b->coeffs[j] == 0) continue;
            temp[i + j] ^= gf216_mul(a->coeffs[i], b->coeffs[j]);
        }
    }
    
    gf216_poly_fit_length(res, rlen);
    memcpy(res->coeffs, temp, rlen * sizeof(uint16_t));
    res->length = rlen;
    gf216_poly_normalise(res);
    
    free(temp);
}

/* ============================================================================
   GF(2^32) POLYNOMIAL OPERATIONS IMPLEMENTATION
   ============================================================================ */

void gf232_poly_init(gf232_poly_t poly) {
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

void gf232_poly_clear(gf232_poly_t poly) {
    if (poly->coeffs) {
        free(poly->coeffs);
        poly->coeffs = NULL;
    }
    poly->length = 0;
    poly->alloc = 0;
}

void gf232_poly_fit_length(gf232_poly_t poly, slong len) {
    if (len > poly->alloc) {
        slong new_alloc = FLINT_MAX(len, poly->alloc * 2);
        poly->coeffs = (gf232_t *)realloc(poly->coeffs, new_alloc * sizeof(gf232_t));
        for (slong i = poly->alloc; i < new_alloc; i++) {
            poly->coeffs[i] = gf232_zero();
        }
        poly->alloc = new_alloc;
    }
}

void gf232_poly_normalise(gf232_poly_t poly) {
    while (poly->length > 0 && gf232_is_zero(&poly->coeffs[poly->length - 1])) {
        poly->length--;
    }
}

inline void gf232_poly_zero(gf232_poly_t poly) {
    poly->length = 0;
}

inline int gf232_poly_is_zero(const gf232_poly_t poly) {
    return poly->length == 0;
}

inline slong gf232_poly_degree(const gf232_poly_t poly) {
    return poly->length - 1;
}

void gf232_poly_set(gf232_poly_t res, const gf232_poly_t poly) {
    if (res == poly) return;
    gf232_poly_fit_length(res, poly->length);
    memcpy(res->coeffs, poly->coeffs, poly->length * sizeof(gf232_t));
    res->length = poly->length;
}

gf232_t gf232_poly_get_coeff(const gf232_poly_t poly, slong i) {
    if (i < poly->length) {
        return poly->coeffs[i];
    } else {
        return gf232_zero();
    }
}

void gf232_poly_set_coeff(gf232_poly_t poly, slong i, const gf232_t *c) {
    gf232_poly_fit_length(poly, i + 1);
    if (i >= poly->length) {
        for (slong j = poly->length; j < i; j++) {
            poly->coeffs[j] = gf232_zero();
        }
        poly->length = i + 1;
    }
    poly->coeffs[i] = *c;
    
    if (gf232_is_zero(c) && i == poly->length - 1) {
        gf232_poly_normalise(poly);
    }
}

void gf232_poly_add(gf232_poly_t res, const gf232_poly_t a, const gf232_poly_t b) {
    slong max_len = FLINT_MAX(a->length, b->length);
    slong min_len = FLINT_MIN(a->length, b->length);
    
    if (max_len == 0) {
        gf232_poly_zero(res);
        return;
    }
    
    gf232_poly_fit_length(res, max_len);
    
    for (slong i = 0; i < min_len; i++) {
        res->coeffs[i] = gf232_add(&a->coeffs[i], &b->coeffs[i]);
    }
    
    if (a->length > b->length) {
        memcpy(res->coeffs + min_len, a->coeffs + min_len, 
               (a->length - min_len) * sizeof(gf232_t));
        res->length = a->length;
    } else if (b->length > a->length) {
        memcpy(res->coeffs + min_len, b->coeffs + min_len, 
               (b->length - min_len) * sizeof(gf232_t));
        res->length = b->length;
    } else {
        res->length = min_len;
    }
    
    gf232_poly_normalise(res);
}

void gf232_poly_scalar_mul(gf232_poly_t res, const gf232_poly_t poly, const gf232_t *c) {
    if (gf232_is_zero(c)) {
        gf232_poly_zero(res);
        return;
    }
    
    gf232_t one = gf232_one();
    if (gf232_equal(c, &one)) {
        gf232_poly_set(res, poly);
        return;
    }
    
    gf232_poly_fit_length(res, poly->length);
    res->length = poly->length;
    
    for (slong i = 0; i < poly->length; i++) {
        res->coeffs[i] = gf232_mul(&poly->coeffs[i], c);
    }
    
    gf232_poly_normalise(res);
}

void gf232_poly_shift_left(gf232_poly_t res, const gf232_poly_t poly, slong n) {
    if (n == 0) {
        gf232_poly_set(res, poly);
        return;
    }
    
    if (gf232_poly_is_zero(poly)) {
        gf232_poly_zero(res);
        return;
    }
    
    slong new_len = poly->length + n;
    gf232_poly_fit_length(res, new_len);
    
    memmove(res->coeffs + n, poly->coeffs, poly->length * sizeof(gf232_t));
    
    for (slong i = 0; i < n; i++) {
        res->coeffs[i] = gf232_zero();
    }
    
    res->length = new_len;
}

/* Schoolbook multiplication for small polynomials */
void gf232_poly_mul_schoolbook(gf232_poly_t res, const gf232_poly_t a, const gf232_poly_t b) {
    if (gf232_poly_is_zero(a) || gf232_poly_is_zero(b)) {
        gf232_poly_zero(res);
        return;
    }
    
    slong rlen = a->length + b->length - 1;
    gf232_t *temp = (gf232_t *)calloc(rlen, sizeof(gf232_t));
    
    for (slong i = 0; i < a->length; i++) {
        for (slong j = 0; j < b->length; j++) {
            gf232_t prod = gf232_mul(&a->coeffs[i], &b->coeffs[j]);
            temp[i + j] = gf232_add(&temp[i + j], &prod);
        }
    }
    
    gf232_poly_fit_length(res, rlen);
    memcpy(res->coeffs, temp, rlen * sizeof(gf232_t));
    res->length = rlen;
    gf232_poly_normalise(res);
    
    free(temp);
}

/* Fast polynomial multiplication for GF(2^32) using Karatsuba algorithm */
void gf232_poly_mul_karatsuba(gf232_poly_t res, const gf232_poly_t a, const gf232_poly_t b) {
    if (gf232_poly_is_zero(a) || gf232_poly_is_zero(b)) {
        gf232_poly_zero(res);
        return;
    }
    
    slong alen = a->length;
    slong blen = b->length;
    
    /* For small polynomials, use schoolbook multiplication */
    if (alen <= 16 || blen <= 16) {
        gf232_poly_mul_schoolbook(res, a, b);
        return;
    }
    
    /* Handle extremely unbalanced polynomials */
    if (alen > 8 * blen || blen > 8 * alen) {
        gf232_poly_mul_schoolbook(res, a, b);
        return;
    }
    
    /* Karatsuba algorithm - use balanced split */
    slong split = (FLINT_MAX(alen, blen) + 1) / 2;
    
    gf232_poly_t a0, a1, b0, b1;
    gf232_poly_t z0, z1, z2, temp1, temp2;
    
    /* Initialize all polynomials */
    gf232_poly_init(a0);
    gf232_poly_init(a1);
    gf232_poly_init(b0);
    gf232_poly_init(b1);
    gf232_poly_init(z0);
    gf232_poly_init(z1);
    gf232_poly_init(z2);
    gf232_poly_init(temp1);
    gf232_poly_init(temp2);
    
    /* Split a: a = a0 + a1*x^split */
    for (slong i = 0; i < FLINT_MIN(split, alen); i++) {
        gf232_poly_set_coeff(a0, i, &a->coeffs[i]);
    }
    
    for (slong i = split; i < alen; i++) {
        gf232_poly_set_coeff(a1, i - split, &a->coeffs[i]);
    }
    
    /* Split b: b = b0 + b1*x^split */
    for (slong i = 0; i < FLINT_MIN(split, blen); i++) {
        gf232_poly_set_coeff(b0, i, &b->coeffs[i]);
    }
    
    for (slong i = split; i < blen; i++) {
        gf232_poly_set_coeff(b1, i - split, &b->coeffs[i]);
    }
    
    /* Compute z0 = a0 * b0 */
    gf232_poly_mul_karatsuba(z0, a0, b0);
    
    /* Compute z2 = a1 * b1 */
    gf232_poly_mul_karatsuba(z2, a1, b1);
    
    /* Compute z1 = (a0 + a1) * (b0 + b1) - z0 - z2 */
    gf232_poly_add(temp1, a0, a1);  /* temp1 = a0 + a1 */
    gf232_poly_add(temp2, b0, b1);  /* temp2 = b0 + b1 */
    gf232_poly_mul_karatsuba(z1, temp1, temp2);  /* z1 = (a0 + a1) * (b0 + b1) */
    
    /* z1 = z1 - z0 - z2 = z1 + z0 + z2 (in GF(2^n)) */
    gf232_poly_add(z1, z1, z0);
    gf232_poly_add(z1, z1, z2);
    
    /* Construct result: res = z0 + z1*x^split + z2*x^(2*split) */
    slong result_len = alen + blen - 1;
    gf232_poly_fit_length(res, result_len);
    
    /* Initialize result to zero */
    for (slong i = 0; i < result_len; i++) {
        res->coeffs[i] = gf232_zero();
    }
    res->length = result_len;
    
    /* Add z0 */
    for (slong i = 0; i < z0->length; i++) {
        res->coeffs[i] = gf232_add(&res->coeffs[i], &z0->coeffs[i]);
    }
    
    /* Add z1 * x^split */
    for (slong i = 0; i < z1->length; i++) {
        if (i + split < result_len) {
            res->coeffs[i + split] = gf232_add(&res->coeffs[i + split], &z1->coeffs[i]);
        }
    }
    
    /* Add z2 * x^(2*split) */
    for (slong i = 0; i < z2->length; i++) {
        if (i + 2*split < result_len) {
            res->coeffs[i + 2*split] = gf232_add(&res->coeffs[i + 2*split], &z2->coeffs[i]);
        }
    }
    
    /* Normalize result */
    gf232_poly_normalise(res);
    
    /* Cleanup */
    gf232_poly_clear(a0);
    gf232_poly_clear(a1);
    gf232_poly_clear(b0);
    gf232_poly_clear(b1);
    gf232_poly_clear(z0);
    gf232_poly_clear(z1);
    gf232_poly_clear(z2);
    gf232_poly_clear(temp1);
    gf232_poly_clear(temp2);
}

/* Replace the original gf232_poly_mul with the Karatsuba version */
void gf232_poly_mul(gf232_poly_t res, const gf232_poly_t a, const gf232_poly_t b) {
    gf232_poly_mul_karatsuba(res, a, b);
}

/* ============================================================================
   GF(2^64) POLYNOMIAL OPERATIONS IMPLEMENTATION
   ============================================================================ */

void gf264_poly_init(gf264_poly_t poly) {
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

void gf264_poly_clear(gf264_poly_t poly) {
    if (poly->coeffs) {
        free(poly->coeffs);
        poly->coeffs = NULL;
    }
    poly->length = 0;
    poly->alloc = 0;
}

void gf264_poly_fit_length(gf264_poly_t poly, slong len) {
    if (len > poly->alloc) {
        slong new_alloc = FLINT_MAX(len, poly->alloc * 2);
        poly->coeffs = (gf264_t *)realloc(poly->coeffs, new_alloc * sizeof(gf264_t));
        for (slong i = poly->alloc; i < new_alloc; i++) {
            poly->coeffs[i] = gf264_zero();
        }
        poly->alloc = new_alloc;
    }
}

void gf264_poly_normalise(gf264_poly_t poly) {
    while (poly->length > 0 && gf264_is_zero(&poly->coeffs[poly->length - 1])) {
        poly->length--;
    }
}

inline void gf264_poly_zero(gf264_poly_t poly) {
    poly->length = 0;
}

inline int gf264_poly_is_zero(const gf264_poly_t poly) {
    return poly->length == 0;
}

inline slong gf264_poly_degree(const gf264_poly_t poly) {
    return poly->length - 1;
}

void gf264_poly_set(gf264_poly_t res, const gf264_poly_t poly) {
    if (res == poly) return;
    gf264_poly_fit_length(res, poly->length);
    memcpy(res->coeffs, poly->coeffs, poly->length * sizeof(gf264_t));
    res->length = poly->length;
}

gf264_t gf264_poly_get_coeff(const gf264_poly_t poly, slong i) {
    if (i < poly->length) {
        return poly->coeffs[i];
    } else {
        return gf264_zero();
    }
}

void gf264_poly_set_coeff(gf264_poly_t poly, slong i, const gf264_t *c) {
    gf264_poly_fit_length(poly, i + 1);
    if (i >= poly->length) {
        for (slong j = poly->length; j < i; j++) {
            poly->coeffs[j] = gf264_zero();
        }
        poly->length = i + 1;
    }
    poly->coeffs[i] = *c;
    
    if (gf264_is_zero(c) && i == poly->length - 1) {
        gf264_poly_normalise(poly);
    }
}

void gf264_poly_add(gf264_poly_t res, const gf264_poly_t a, const gf264_poly_t b) {
    slong max_len = FLINT_MAX(a->length, b->length);
    slong min_len = FLINT_MIN(a->length, b->length);
    
    if (max_len == 0) {
        gf264_poly_zero(res);
        return;
    }
    
    gf264_poly_fit_length(res, max_len);
    
    for (slong i = 0; i < min_len; i++) {
        res->coeffs[i] = gf264_add(&a->coeffs[i], &b->coeffs[i]);
    }
    
    if (a->length > b->length) {
        memcpy(res->coeffs + min_len, a->coeffs + min_len, 
               (a->length - min_len) * sizeof(gf264_t));
        res->length = a->length;
    } else if (b->length > a->length) {
        memcpy(res->coeffs + min_len, b->coeffs + min_len, 
               (b->length - min_len) * sizeof(gf264_t));
        res->length = b->length;
    } else {
        res->length = min_len;
    }
    
    gf264_poly_normalise(res);
}

void gf264_poly_scalar_mul(gf264_poly_t res, const gf264_poly_t poly, const gf264_t *c) {
    if (gf264_is_zero(c)) {
        gf264_poly_zero(res);
        return;
    }
    
    gf264_t one = gf264_one();
    if (gf264_equal(c, &one)) {
        gf264_poly_set(res, poly);
        return;
    }
    
    gf264_poly_fit_length(res, poly->length);
    res->length = poly->length;
    
    for (slong i = 0; i < poly->length; i++) {
        res->coeffs[i] = gf264_mul(&poly->coeffs[i], c);
    }
    
    gf264_poly_normalise(res);
}

void gf264_poly_mul(gf264_poly_t res, const gf264_poly_t a, const gf264_poly_t b) {
    if (gf264_poly_is_zero(a) || gf264_poly_is_zero(b)) {
        gf264_poly_zero(res);
        return;
    }
    
    slong rlen = a->length + b->length - 1;
    gf264_t *temp = (gf264_t *)calloc(rlen, sizeof(gf264_t));
    
    for (slong i = 0; i < a->length; i++) {
        for (slong j = 0; j < b->length; j++) {
            gf264_t prod = gf264_mul(&a->coeffs[i], &b->coeffs[j]);
            temp[i + j] = gf264_add(&temp[i + j], &prod);
        }
    }
    
    gf264_poly_fit_length(res, rlen);
    memcpy(res->coeffs, temp, rlen * sizeof(gf264_t));
    res->length = rlen;
    gf264_poly_normalise(res);
    
    free(temp);
}

/* ============================================================================
   GF(2^128) POLYNOMIAL OPERATIONS IMPLEMENTATION
   ============================================================================ */

void gf2128_poly_init(gf2128_poly_t poly) {
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

void gf2128_poly_clear(gf2128_poly_t poly) {
    if (poly->coeffs) {
        free(poly->coeffs);
        poly->coeffs = NULL;
    }
    poly->length = 0;
    poly->alloc = 0;
}

void gf2128_poly_fit_length(gf2128_poly_t poly, slong len) {
    if (len > poly->alloc) {
        slong new_alloc = FLINT_MAX(len, poly->alloc * 2);
        poly->coeffs = (gf2128_t *)realloc(poly->coeffs, new_alloc * sizeof(gf2128_t));
        for (slong i = poly->alloc; i < new_alloc; i++) {
            poly->coeffs[i] = gf2128_zero();
        }
        poly->alloc = new_alloc;
    }
}

void gf2128_poly_normalise(gf2128_poly_t poly) {
    while (poly->length > 0 && gf2128_is_zero(&poly->coeffs[poly->length - 1])) {
        poly->length--;
    }
}

/* Helper function to set polynomial coefficient */
void gf2128_poly_set_coeff(gf2128_poly_t poly, slong i, const gf2128_t *c) {
    gf2128_poly_fit_length(poly, i + 1);
    if (i >= poly->length) {
        // Zero out coefficients between old length and i
        for (slong j = poly->length; j < i; j++) {
            poly->coeffs[j] = gf2128_zero();
        }
        poly->length = i + 1;
    }
    poly->coeffs[i] = *c;
    
    // Update length if we're setting the coefficient to zero
    if (gf2128_is_zero(c) && i == poly->length - 1) {
        gf2128_poly_normalise(poly);
    }
}

inline void gf2128_poly_zero(gf2128_poly_t poly) {
    poly->length = 0;
}

inline int gf2128_poly_is_zero(const gf2128_poly_t poly) {
    return poly->length == 0;
}

inline slong gf2128_poly_degree(const gf2128_poly_t poly) {
    return poly->length - 1;
}

void gf2128_poly_set(gf2128_poly_t res, const gf2128_poly_t poly) {
    if (res == poly) return;
    gf2128_poly_fit_length(res, poly->length);
    memcpy(res->coeffs, poly->coeffs, poly->length * sizeof(gf2128_t));
    res->length = poly->length;
}

gf2128_t gf2128_poly_get_coeff(const gf2128_poly_t poly, slong i) {
    if (i < poly->length) {
        return poly->coeffs[i];
    } else {
        return gf2128_zero();
    }
}

void gf2128_poly_add(gf2128_poly_t res, const gf2128_poly_t a, const gf2128_poly_t b) {
    slong max_len = FLINT_MAX(a->length, b->length);
    slong min_len = FLINT_MIN(a->length, b->length);
    
    if (max_len == 0) {
        gf2128_poly_zero(res);
        return;
    }
    
    gf2128_poly_fit_length(res, max_len);
    
    for (slong i = 0; i < min_len; i++) {
        res->coeffs[i] = gf2128_add(&a->coeffs[i], &b->coeffs[i]);
    }
    
    if (a->length > b->length) {
        memcpy(res->coeffs + min_len, a->coeffs + min_len, 
               (a->length - min_len) * sizeof(gf2128_t));
        res->length = a->length;
    } else if (b->length > a->length) {
        memcpy(res->coeffs + min_len, b->coeffs + min_len, 
               (b->length - min_len) * sizeof(gf2128_t));
        res->length = b->length;
    } else {
        res->length = min_len;
    }
    
    gf2128_poly_normalise(res);
}

void gf2128_poly_scalar_mul(gf2128_poly_t res, const gf2128_poly_t poly, const gf2128_t *c) {
    if (gf2128_is_zero(c)) {
        gf2128_poly_zero(res);
        return;
    }
    
    gf2128_t one = gf2128_one();
    if (gf2128_equal(c, &one)) {
        gf2128_poly_set(res, poly);
        return;
    }
    
    gf2128_poly_fit_length(res, poly->length);
    res->length = poly->length;
    
    for (slong i = 0; i < poly->length; i++) {
        res->coeffs[i] = gf2128_mul(&poly->coeffs[i], c);
    }
    
    gf2128_poly_normalise(res);
}

void gf2128_poly_shift_left(gf2128_poly_t res, const gf2128_poly_t poly, slong n) {
    if (n == 0) {
        gf2128_poly_set(res, poly);
        return;
    }
    
    if (gf2128_poly_is_zero(poly)) {
        gf2128_poly_zero(res);
        return;
    }
    
    slong new_len = poly->length + n;
    gf2128_poly_fit_length(res, new_len);
    
    memmove(res->coeffs + n, poly->coeffs, poly->length * sizeof(gf2128_t));
    
    for (slong i = 0; i < n; i++) {
        res->coeffs[i] = gf2128_zero();
    }
    
    res->length = new_len;
}

void gf2128_poly_mul_schoolbook(gf2128_poly_t res, const gf2128_poly_t a, const gf2128_poly_t b) {
    if (gf2128_poly_is_zero(a) || gf2128_poly_is_zero(b)) {
        gf2128_poly_zero(res);
        return;
    }
    
    slong rlen = a->length + b->length - 1;
    gf2128_t *temp = (gf2128_t *)calloc(rlen, sizeof(gf2128_t));
    
    for (slong i = 0; i < a->length; i++) {
        for (slong j = 0; j < b->length; j++) {
            gf2128_t prod = gf2128_mul(&a->coeffs[i], &b->coeffs[j]);
            temp[i + j] = gf2128_add(&temp[i + j], &prod);
        }
    }
    
    gf2128_poly_fit_length(res, rlen);
    memcpy(res->coeffs, temp, rlen * sizeof(gf2128_t));
    res->length = rlen;
    gf2128_poly_normalise(res);
    
    free(temp);
}

/* Fast polynomial multiplication for GF(2^128) using Karatsuba algorithm */
void gf2128_poly_mul_karatsuba(gf2128_poly_t res, const gf2128_poly_t a, const gf2128_poly_t b) {
    if (gf2128_poly_is_zero(a) || gf2128_poly_is_zero(b)) {
        gf2128_poly_zero(res);
        return;
    }
    
    slong alen = a->length;
    slong blen = b->length;
    
    /* For small polynomials, use schoolbook multiplication */
    if (alen <= 8 || blen <= 8) {
        gf2128_poly_mul_schoolbook(res, a, b);
        return;
    }
    
    /* Handle extremely unbalanced polynomials */
    if (alen > 8 * blen || blen > 8 * alen) {
        gf2128_poly_mul_schoolbook(res, a, b);
        return;
    }
    
    /* Karatsuba algorithm - use balanced split */
    slong split = (FLINT_MAX(alen, blen) + 1) / 2;
    
    gf2128_poly_t a0, a1, b0, b1;
    gf2128_poly_t z0, z1, z2, temp1, temp2;
    
    /* Initialize all polynomials */
    gf2128_poly_init(a0);
    gf2128_poly_init(a1);
    gf2128_poly_init(b0);
    gf2128_poly_init(b1);
    gf2128_poly_init(z0);
    gf2128_poly_init(z1);
    gf2128_poly_init(z2);
    gf2128_poly_init(temp1);
    gf2128_poly_init(temp2);
    
    /* Split a: a = a0 + a1*x^split */
    for (slong i = 0; i < FLINT_MIN(split, alen); i++) {
        gf2128_poly_set_coeff(a0, i, &a->coeffs[i]);
    }
    
    for (slong i = split; i < alen; i++) {
        gf2128_poly_set_coeff(a1, i - split, &a->coeffs[i]);
    }
    
    /* Split b: b = b0 + b1*x^split */
    for (slong i = 0; i < FLINT_MIN(split, blen); i++) {
        gf2128_poly_set_coeff(b0, i, &b->coeffs[i]);
    }
    
    for (slong i = split; i < blen; i++) {
        gf2128_poly_set_coeff(b1, i - split, &b->coeffs[i]);
    }
    
    /* Compute z0 = a0 * b0 */
    gf2128_poly_mul_karatsuba(z0, a0, b0);
    
    /* Compute z2 = a1 * b1 */
    gf2128_poly_mul_karatsuba(z2, a1, b1);
    
    /* Compute z1 = (a0 + a1) * (b0 + b1) - z0 - z2 */
    gf2128_poly_add(temp1, a0, a1);  /* temp1 = a0 + a1 */
    gf2128_poly_add(temp2, b0, b1);  /* temp2 = b0 + b1 */
    gf2128_poly_mul_karatsuba(z1, temp1, temp2);  /* z1 = (a0 + a1) * (b0 + b1) */
    
    /* z1 = z1 - z0 - z2 = z1 + z0 + z2 (in GF(2^n)) */
    gf2128_poly_add(z1, z1, z0);
    gf2128_poly_add(z1, z1, z2);
    
    /* Construct result: res = z0 + z1*x^split + z2*x^(2*split) */
    slong result_len = alen + blen - 1;
    gf2128_poly_fit_length(res, result_len);
    
    /* Initialize result to zero */
    for (slong i = 0; i < result_len; i++) {
        res->coeffs[i] = gf2128_zero();
    }
    res->length = result_len;
    
    /* Add z0 */
    for (slong i = 0; i < z0->length; i++) {
        res->coeffs[i] = gf2128_add(&res->coeffs[i], &z0->coeffs[i]);
    }
    
    /* Add z1 * x^split */
    for (slong i = 0; i < z1->length; i++) {
        if (i + split < result_len) {
            res->coeffs[i + split] = gf2128_add(&res->coeffs[i + split], &z1->coeffs[i]);
        }
    }
    
    /* Add z2 * x^(2*split) */
    for (slong i = 0; i < z2->length; i++) {
        if (i + 2*split < result_len) {
            res->coeffs[i + 2*split] = gf2128_add(&res->coeffs[i + 2*split], &z2->coeffs[i]);
        }
    }
    
    /* Normalize result */
    gf2128_poly_normalise(res);
    
    /* Cleanup */
    gf2128_poly_clear(a0);
    gf2128_poly_clear(a1);
    gf2128_poly_clear(b0);
    gf2128_poly_clear(b1);
    gf2128_poly_clear(z0);
    gf2128_poly_clear(z1);
    gf2128_poly_clear(z2);
    gf2128_poly_clear(temp1);
    gf2128_poly_clear(temp2);
}

void gf2128_poly_mul(gf2128_poly_t res, const gf2128_poly_t a, const gf2128_poly_t b){
    gf2128_poly_mul_karatsuba(res, a, b);
}

/* ============================================================================
   MATRIX OPERATIONS IMPLEMENTATION
   ============================================================================ */

/* GF(2^8) matrix operations */
void gf28_poly_mat_init(gf28_poly_mat_t mat, slong rows, slong cols) {
    mat->entries = NULL;
    mat->rows = NULL;
    
    if (rows > 0 && cols > 0) {
        mat->entries = (gf28_poly_struct *)malloc(rows * cols * sizeof(gf28_poly_struct));
        mat->rows = (gf28_poly_struct **)malloc(rows * sizeof(gf28_poly_struct *));
        
        for (slong i = 0; i < rows * cols; i++) {
            gf28_poly_init(mat->entries + i);
        }
        
        for (slong i = 0; i < rows; i++) {
            mat->rows[i] = mat->entries + i * cols;
        }
    }
    
    mat->r = rows;
    mat->c = cols;
}

void gf28_poly_mat_clear(gf28_poly_mat_t mat) {
    if (mat->entries != NULL) {
        for (slong i = 0; i < mat->r * mat->c; i++) {
            gf28_poly_clear(mat->entries + i);
        }
        free(mat->entries);
        free(mat->rows);
    }
}

inline gf28_poly_struct *gf28_poly_mat_entry(gf28_poly_mat_t mat, slong i, slong j) {
    return mat->rows[i] + j;
}

void gf28_poly_mat_swap_rows(gf28_poly_mat_t mat, slong r, slong s) {
    if (r != s) {
        gf28_poly_struct *tmp = mat->rows[r];
        mat->rows[r] = mat->rows[s];
        mat->rows[s] = tmp;
    }
}

/* GF(2^16) matrix operations */
void gf216_poly_mat_init(gf216_poly_mat_t mat, slong rows, slong cols) {
    mat->entries = NULL;
    mat->rows = NULL;
    
    if (rows > 0 && cols > 0) {
        mat->entries = (gf216_poly_struct *)malloc(rows * cols * sizeof(gf216_poly_struct));
        mat->rows = (gf216_poly_struct **)malloc(rows * sizeof(gf216_poly_struct *));
        
        for (slong i = 0; i < rows * cols; i++) {
            gf216_poly_init(mat->entries + i);
        }
        
        for (slong i = 0; i < rows; i++) {
            mat->rows[i] = mat->entries + i * cols;
        }
    }
    
    mat->r = rows;
    mat->c = cols;
}

void gf216_poly_mat_clear(gf216_poly_mat_t mat) {
    if (mat->entries != NULL) {
        for (slong i = 0; i < mat->r * mat->c; i++) {
            gf216_poly_clear(mat->entries + i);
        }
        free(mat->entries);
        free(mat->rows);
    }
}

inline gf216_poly_struct *gf216_poly_mat_entry(gf216_poly_mat_t mat, slong i, slong j) {
    return mat->rows[i] + j;
}

void gf216_poly_mat_swap_rows(gf216_poly_mat_t mat, slong r, slong s) {
    if (r != s) {
        gf216_poly_struct *tmp = mat->rows[r];
        mat->rows[r] = mat->rows[s];
        mat->rows[s] = tmp;
    }
}

void gf216_poly_mat_permute_rows(gf216_poly_mat_t mat, const slong *perm) {
    gf216_poly_struct **new_rows = (gf216_poly_struct **)malloc(mat->r * sizeof(gf216_poly_struct *));
    
    for (slong i = 0; i < mat->r; i++) {
        new_rows[i] = mat->rows[perm[i]];
    }
    
    memcpy(mat->rows, new_rows, mat->r * sizeof(gf216_poly_struct *));
    free(new_rows);
}

/* GF(2^32) matrix operations */
void gf232_poly_mat_init(gf232_poly_mat_t mat, slong rows, slong cols) {
    mat->entries = NULL;
    mat->rows = NULL;
    
    if (rows > 0 && cols > 0) {
        mat->entries = (gf232_poly_struct *)malloc(rows * cols * sizeof(gf232_poly_struct));
        mat->rows = (gf232_poly_struct **)malloc(rows * sizeof(gf232_poly_struct *));
        
        for (slong i = 0; i < rows * cols; i++) {
            gf232_poly_init(mat->entries + i);
        }
        
        for (slong i = 0; i < rows; i++) {
            mat->rows[i] = mat->entries + i * cols;
        }
    }
    
    mat->r = rows;
    mat->c = cols;
}

void gf232_poly_mat_clear(gf232_poly_mat_t mat) {
    if (mat->entries != NULL) {
        for (slong i = 0; i < mat->r * mat->c; i++) {
            gf232_poly_clear(mat->entries + i);
        }
        free(mat->entries);
        free(mat->rows);
    }
}

inline gf232_poly_struct *gf232_poly_mat_entry(gf232_poly_mat_t mat, slong i, slong j) {
    return mat->rows[i] + j;
}

/* GF(2^64) matrix operations */
void gf264_poly_mat_init(gf264_poly_mat_t mat, slong rows, slong cols) {
    mat->entries = NULL;
    mat->rows = NULL;
    
    if (rows > 0 && cols > 0) {
        mat->entries = (gf264_poly_struct *)malloc(rows * cols * sizeof(gf264_poly_struct));
        mat->rows = (gf264_poly_struct **)malloc(rows * sizeof(gf264_poly_struct *));
        
        for (slong i = 0; i < rows * cols; i++) {
            gf264_poly_init(mat->entries + i);
        }
        
        for (slong i = 0; i < rows; i++) {
            mat->rows[i] = mat->entries + i * cols;
        }
    }
    
    mat->r = rows;
    mat->c = cols;
}

void gf264_poly_mat_clear(gf264_poly_mat_t mat) {
    if (mat->entries != NULL) {
        for (slong i = 0; i < mat->r * mat->c; i++) {
            gf264_poly_clear(mat->entries + i);
        }
        free(mat->entries);
        free(mat->rows);
    }
}

inline gf264_poly_struct *gf264_poly_mat_entry(gf264_poly_mat_t mat, slong i, slong j) {
    return mat->rows[i] + j;
}

/* GF(2^128) matrix operations */
void gf2128_poly_mat_init(gf2128_poly_mat_t mat, slong rows, slong cols) {
    mat->entries = NULL;
    mat->rows = NULL;
    
    if (rows > 0 && cols > 0) {
        mat->entries = (gf2128_poly_struct *)malloc(rows * cols * sizeof(gf2128_poly_struct));
        mat->rows = (gf2128_poly_struct **)malloc(rows * sizeof(gf2128_poly_struct *));
        
        for (slong i = 0; i < rows * cols; i++) {
           gf2128_poly_init(mat->entries + i);
       }
       
       for (slong i = 0; i < rows; i++) {
           mat->rows[i] = mat->entries + i * cols;
       }
   }
   
   mat->r = rows;
   mat->c = cols;
}

void gf2128_poly_mat_clear(gf2128_poly_mat_t mat) {
   if (mat->entries != NULL) {
       for (slong i = 0; i < mat->r * mat->c; i++) {
           gf2128_poly_clear(mat->entries + i);
       }
       free(mat->entries);
       free(mat->rows);
   }
}

inline gf2128_poly_struct *gf2128_poly_mat_entry(gf2128_poly_mat_t mat, slong i, slong j) {
   return mat->rows[i] + j;
}

inline const gf2128_poly_struct *gf2128_poly_mat_entry_const(const gf2128_poly_mat_t mat, slong i, slong j) {
   return mat->rows[i] + j;
}

void gf2128_poly_mat_swap_rows(gf2128_poly_mat_t mat, slong r, slong s) {
   if (r != s) {
       gf2128_poly_struct *tmp = mat->rows[r];
       mat->rows[r] = mat->rows[s];
       mat->rows[s] = tmp;
   }
}

/* FLINT polynomial matrix operations */
void fq_nmod_poly_mat_init(fq_nmod_poly_mat_t mat, slong rows, slong cols,
                          const fq_nmod_ctx_t ctx) {
    mat->entries = NULL;
    mat->rows = NULL;
    
    if (rows > 0 && cols > 0) {
        mat->entries = (fq_nmod_poly_struct *)flint_malloc(rows * cols * sizeof(fq_nmod_poly_struct));
        mat->rows = (fq_nmod_poly_struct **)flint_malloc(rows * sizeof(fq_nmod_poly_struct *));
        
        for (slong i = 0; i < rows * cols; i++)
            fq_nmod_poly_init(mat->entries + i, ctx);
        
        for (slong i = 0; i < rows; i++)
            mat->rows[i] = mat->entries + i * cols;
    }
    
    mat->r = rows;
    mat->c = cols;
    mat->ctx = (fq_nmod_ctx_struct *)ctx;
}

void fq_nmod_poly_mat_clear(fq_nmod_poly_mat_t mat, const fq_nmod_ctx_t ctx) {
    if (mat->entries != NULL) {
        for (slong i = 0; i < mat->r * mat->c; i++)
            fq_nmod_poly_clear(mat->entries + i, ctx);
        
        flint_free(mat->entries);
        flint_free(mat->rows);
    }
}

inline fq_nmod_poly_struct* fq_nmod_poly_mat_entry(const fq_nmod_poly_mat_t mat, 
                                                          slong i, slong j) {
    return mat->rows[i] + j;
}

/* ============================================================================
   CONVERSION FUNCTIONS BETWEEN POLYNOMIALS IMPLEMENTATION
   ============================================================================ */

/* Convert polynomials */
void fq_nmod_poly_to_gf2128_poly(gf2128_poly_t res,
                                        const fq_nmod_poly_struct *poly,
                                        const fq_nmod_ctx_t ctx) {
   slong len = fq_nmod_poly_length(poly, ctx);
   if (len == 0) {
       gf2128_poly_zero(res);
       return;
   }
   
   gf2128_poly_fit_length(res, len);
   res->length = len;
   
   for (slong i = 0; i < len; i++) {
       fq_nmod_t coeff;
       fq_nmod_init(coeff, ctx);
       fq_nmod_poly_get_coeff(coeff, poly, i, ctx);
       
       res->coeffs[i] = fq_nmod_to_gf2128(coeff, ctx);
       
       fq_nmod_clear(coeff, ctx);
   }
   
   gf2128_poly_normalise(res);
}

void gf2128_poly_to_fq_nmod_poly(fq_nmod_poly_struct *res,
                                        const gf2128_poly_t poly,
                                        const fq_nmod_ctx_t ctx) {
   fq_nmod_poly_zero(res, ctx);
   
   for (slong i = 0; i < poly->length; i++) {
       if (!gf2128_is_zero(&poly->coeffs[i])) {
           fq_nmod_t coeff;
           fq_nmod_init(coeff, ctx);
           
           gf2128_to_fq_nmod(coeff, &poly->coeffs[i], ctx);
           fq_nmod_poly_set_coeff(res, i, coeff, ctx);
           
           fq_nmod_clear(coeff, ctx);
       }
   }
}

void fq_nmod_poly_to_gf28_poly(gf28_poly_t res,
                                     const fq_nmod_poly_struct *poly,
                                     const fq_nmod_ctx_t ctx) {
   if (!g_gf28_conversion || !g_gf28_conversion->initialized) {
       init_gf28_conversion(ctx);
   }
   
   slong len = fq_nmod_poly_length(poly, ctx);
   if (len == 0) {
       gf28_poly_zero(res);
       return;
   }
   
   gf28_poly_fit_length(res, len);
   res->length = len;
   
   for (slong i = 0; i < len; i++) {
       fq_nmod_t coeff;
       fq_nmod_init(coeff, ctx);
       fq_nmod_poly_get_coeff(coeff, poly, i, ctx);
       
       res->coeffs[i] = fq_nmod_to_gf28_elem(coeff, ctx);
       
       fq_nmod_clear(coeff, ctx);
   }
   
   gf28_poly_normalise(res);
}

void gf28_poly_to_fq_nmod_poly(fq_nmod_poly_struct *res,
                                     const gf28_poly_t poly,
                                     const fq_nmod_ctx_t ctx) {
   if (!g_gf28_conversion || !g_gf28_conversion->initialized) {
       init_gf28_conversion(ctx);
   }
   
   fq_nmod_poly_zero(res, ctx);
   
   for (slong i = 0; i < poly->length; i++) {
       if (poly->coeffs[i] != 0) {
           fq_nmod_t coeff;
           fq_nmod_init(coeff, ctx);
           
           gf28_elem_to_fq_nmod(coeff, poly->coeffs[i], ctx);
           fq_nmod_poly_set_coeff(res, i, coeff, ctx);
           
           fq_nmod_clear(coeff, ctx);
       }
   }
}

void fq_nmod_poly_to_gf216_poly(gf216_poly_t res,
                                     const fq_nmod_poly_struct *poly,
                                     const fq_nmod_ctx_t ctx) {
   if (!g_gf216_conversion || !g_gf216_conversion->initialized) {
       init_gf216_conversion(ctx);
   }
   
   slong len = fq_nmod_poly_length(poly, ctx);
   if (len == 0) {
       gf216_poly_zero(res);
       return;
   }
   
   gf216_poly_fit_length(res, len);
   res->length = len;
   
   for (slong i = 0; i < len; i++) {
       fq_nmod_t coeff;
       fq_nmod_init(coeff, ctx);
       fq_nmod_poly_get_coeff(coeff, poly, i, ctx);
       
       res->coeffs[i] = fq_nmod_to_gf216_elem(coeff, ctx);
       
       fq_nmod_clear(coeff, ctx);
   }
   
   gf216_poly_normalise(res);
}

void gf216_poly_to_fq_nmod_poly(fq_nmod_poly_struct *res,
                                     const gf216_poly_t poly,
                                     const fq_nmod_ctx_t ctx) {
   if (!g_gf216_conversion || !g_gf216_conversion->initialized) {
       init_gf216_conversion(ctx);
   }
   
   fq_nmod_poly_zero(res, ctx);
   
   for (slong i = 0; i < poly->length; i++) {
       if (poly->coeffs[i] != 0) {
           fq_nmod_t coeff;
           fq_nmod_init(coeff, ctx);
           
           gf216_elem_to_fq_nmod(coeff, poly->coeffs[i], ctx);
           fq_nmod_poly_set_coeff(res, i, coeff, ctx);
           
           fq_nmod_clear(coeff, ctx);
       }
   }
}

/* Convert polynomials for GF(2^32) and GF(2^64) */
void fq_nmod_poly_to_gf232_poly(gf232_poly_t res,
                                       const fq_nmod_poly_struct *poly,
                                       const fq_nmod_ctx_t ctx) {
    slong len = fq_nmod_poly_length(poly, ctx);
    if (len == 0) {
        gf232_poly_zero(res);
        return;
    }
    
    gf232_poly_fit_length(res, len);
    res->length = len;
    
    for (slong i = 0; i < len; i++) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, ctx);
        fq_nmod_poly_get_coeff(coeff, poly, i, ctx);
        
        res->coeffs[i] = fq_nmod_to_gf232(coeff, ctx);
        
        fq_nmod_clear(coeff, ctx);
    }
    
    gf232_poly_normalise(res);
}

void gf232_poly_to_fq_nmod_poly(fq_nmod_poly_struct *res,
                                       const gf232_poly_t poly,
                                       const fq_nmod_ctx_t ctx) {
    fq_nmod_poly_zero(res, ctx);
    
    for (slong i = 0; i < poly->length; i++) {
        if (!gf232_is_zero(&poly->coeffs[i])) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            
            gf232_to_fq_nmod(coeff, &poly->coeffs[i], ctx);
            fq_nmod_poly_set_coeff(res, i, coeff, ctx);
            
            fq_nmod_clear(coeff, ctx);
        }
    }
}

void fq_nmod_poly_to_gf264_poly(gf264_poly_t res,
                                       const fq_nmod_poly_struct *poly,
                                       const fq_nmod_ctx_t ctx) {
    slong len = fq_nmod_poly_length(poly, ctx);
    if (len == 0) {
        gf264_poly_zero(res);
        return;
    }
    
    gf264_poly_fit_length(res, len);
    res->length = len;
    
    for (slong i = 0; i < len; i++) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, ctx);
        fq_nmod_poly_get_coeff(coeff, poly, i, ctx);
        
        res->coeffs[i] = fq_nmod_to_gf264(coeff, ctx);
        
        fq_nmod_clear(coeff, ctx);
    }
    
    gf264_poly_normalise(res);
}

void gf264_poly_to_fq_nmod_poly(fq_nmod_poly_struct *res,
                                       const gf264_poly_t poly,
                                       const fq_nmod_ctx_t ctx) {
    fq_nmod_poly_zero(res, ctx);
    
    for (slong i = 0; i < poly->length; i++) {
        if (!gf264_is_zero(&poly->coeffs[i])) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            
            gf264_to_fq_nmod(coeff, &poly->coeffs[i], ctx);
            fq_nmod_poly_set_coeff(res, i, coeff, ctx);
            
            fq_nmod_clear(coeff, ctx);
        }
    }
}

/* Convert matrices */
void fq_nmod_poly_mat_to_gf2128(gf2128_poly_mat_t res,
                                       const fq_nmod_poly_mat_t mat,
                                       const fq_nmod_ctx_t ctx) {
   for (slong i = 0; i < mat->r; i++) {
       for (slong j = 0; j < mat->c; j++) {
           fq_nmod_poly_to_gf2128_poly(
               gf2128_poly_mat_entry(res, i, j),
               fq_nmod_poly_mat_entry(mat, i, j),
               ctx);
       }
   }
}

void gf2128_poly_mat_to_fq_nmod(fq_nmod_poly_mat_t res,
                                       const gf2128_poly_mat_t mat,
                                       const fq_nmod_ctx_t ctx) {
   for (slong i = 0; i < mat->r; i++) {
       for (slong j = 0; j < mat->c; j++) {
           gf2128_poly_to_fq_nmod_poly(
               fq_nmod_poly_mat_entry(res, i, j),
               gf2128_poly_mat_entry_const(mat, i, j),
               ctx);
       }
   }
}

void fq_nmod_poly_mat_to_gf28(gf28_poly_mat_t res,
                                    const fq_nmod_poly_mat_t mat,
                                    const fq_nmod_ctx_t ctx) {
   for (slong i = 0; i < mat->r; i++) {
       for (slong j = 0; j < mat->c; j++) {
           fq_nmod_poly_to_gf28_poly(
               gf28_poly_mat_entry(res, i, j),
               fq_nmod_poly_mat_entry(mat, i, j),
               ctx);
       }
   }
}

void fq_nmod_poly_mat_to_gf216(gf216_poly_mat_t res,
                                    const fq_nmod_poly_mat_t mat,
                                    const fq_nmod_ctx_t ctx) {
   for (slong i = 0; i < mat->r; i++) {
       for (slong j = 0; j < mat->c; j++) {
           fq_nmod_poly_to_gf216_poly(
               gf216_poly_mat_entry(res, i, j),
               fq_nmod_poly_mat_entry(mat, i, j),
               ctx);
       }
   }
}