#include <flint/longlong.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"

static inline ulong
_nmod_vec_dot_bound_limbs_unbalanced(ulong len, ulong max_bits1, ulong max_bits2)
{
    const ulong a1 = (max_bits1 == FLINT_BITS) ? (UWORD_MAX) : (UWORD(1) << max_bits1) - 1;
    const ulong a2 = (max_bits2 == FLINT_BITS) ? (UWORD_MAX) : (UWORD(1) << max_bits2) - 1;

    ulong t2, t1, t0, u1, u0;
    umul_ppmm(t1, t0, a1, a2);
    umul_ppmm(t2, t1, t1, len);
    umul_ppmm(u1, u0, t0, len);
    add_sssaaaaaa(t2, t1, t0,  t2, t1, UWORD(0),  UWORD(0), u1, u0);

    if (t2 != 0)
        return 3;
    if (t1 != 0)
        return 2;
    return (t0 != 0);
}

static inline ulong
_nmod_vec_dot_product_1(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod)
{
    ulong res;
    slong i;
    _NMOD_VEC_DOT1(res, i, len, v1[i], v2[i], mod);
    return res;
}

static inline ulong
_nmod_vec_dot_product_2(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod)
{
    ulong res;
    slong i;
    _NMOD_VEC_DOT2(res, i, len, v1[i], v2[i], mod);
    return res;
}

static inline ulong
_nmod_vec_dot_product_3(nn_srcptr v1, nn_srcptr v2,
                        ulong len, ulong max_bits1, ulong max_bits2,
                        nmod_t mod)
{
    const ulong log_nt = 2 * FLINT_BITS - (max_bits1 + max_bits2);
    const ulong num_terms = (log_nt < FLINT_BITS) ? (UWORD(1) << log_nt) : (UWORD_MAX);

    ulong s0, s1, u0, u1;
    ulong t2 = UWORD(0);
    ulong t1 = UWORD(0);
    ulong t0 = UWORD(0);

    ulong i = 0;
    if (num_terms >= 8)
        for (; i + 7 < len; i += 8)
        {
            umul_ppmm(u1, u0, v1[i+0], v2[i+0]);
            umul_ppmm(s1, s0, v1[i+1], v2[i+1]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[i+2], v2[i+2]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[i+3], v2[i+3]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[i+4], v2[i+4]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[i+5], v2[i+5]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[i+6], v2[i+6]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[i+7], v2[i+7]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), u1, u0);
        }
    else
        for (; i + num_terms < len; i += num_terms)
        {
            umul_ppmm(u1, u0, v1[i], v2[i]);
            for (ulong ii = 1; ii < num_terms; ii++)
            {
                umul_ppmm(s1, s0, v1[i+ii], v2[i+ii]);
                add_ssaaaa(u1, u0, u1, u0, s1, s0);
            }
            add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), u1, u0);
        }

    u0 = UWORD(0);
    u1 = UWORD(0);
    for (; i < len; i++)
    {
        umul_ppmm(s1, s0, v1[i], v2[i]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
    }

    add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), u1, u0);
    NMOD_RED(t2, t2, mod);

    ulong res;
    NMOD_RED3(res, t2, t1, t0, mod);
    return res;
}

ulong nmod_vec_dot_product_unbalanced(nn_srcptr v1, nn_srcptr v2,
                                      ulong len, ulong max_bits1, ulong max_bits2,
                                      nmod_t mod)
{
    const ulong n_limbs = _nmod_vec_dot_bound_limbs_unbalanced(len, max_bits1, max_bits2);

    if (n_limbs == 2)
        return _nmod_vec_dot_product_2(v1, v2, len, mod);
    if (n_limbs == 3)
        return _nmod_vec_dot_product_3(v1, v2, len, max_bits1, max_bits2, mod);
    if (n_limbs == 1)
        return _nmod_vec_dot_product_1(v1, v2, len, mod);
    return UWORD(0);
}
