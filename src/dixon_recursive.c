/*
 * dixon_recursive.c - Recursive block-matrix Dixon construction.
 *
 * This implementation follows the block recursion from Zhao/Fu (2005) and
 * the Maple reference in paper/dixonmaple/快速计算.mpl. At each elimination
 * level x_N, we:
 *
 *   1. split every polynomial by the degree of x_N;
 *   2. recursively build the lower-level Dixon blocks DD[k][j];
 *   3. build the Sylvester-style blocks P[i][k];
 *   4. combine them through
 *        D_{i,j} = sum_k (-1)^k P[i][k] * DD[k][j]
 *      followed by the block recurrence
 *        D_{i,j} <- D_{i,j} + D_{i-1,j+1}.
 *
 * The result is the Dixon matrix directly, without first constructing the
 * Dixon polynomial.
 */

#include "dixon_recursive.h"

#include <stdarg.h>
#include <string.h>

#include "dixon_flint.h"

#ifdef _OPENMP
#include <omp.h>
#endif

double get_wall_time(void);

typedef struct {
    slong rows;
    slong cols;
    slong npars;
    const fq_nmod_ctx_struct *ctx;
    fq_mvpoly_t *entries;
} fast_dixon_matrix_t;

typedef struct {
    slong rows;
    slong total_nonzero;
    slong *row_counts;
    slong *row_offsets;
    slong *col_indices;
} fast_dixon_matrix_support_t;

#define FAST_DIXON_ENTRY(mat, i, j) ((mat)->entries[(i) * (mat)->cols + (j)])

#define FAST_DIXON_MAX_LEVELS 32

typedef struct {
    slong calls;
    slong base_case_calls;
    slong current_vars;
    slong current_degree;
    slong lower_size;
    slong output_rows;
    slong output_cols;
    slong s_blocks;
    slong f_blocks;
    slong d_block_products;
    slong f_recursive_nodes;
    slong f_leaf_calls;
    slong f_valid_leaves;
    slong f_zero_skips;
    slong f_pruned_branches;
    slong cache_hits;
    slong cache_misses;
    double total_time;
    double split_time;
    double s_build_time;
    double f_build_time;
    double f_child_build_time;
    double mul_time;
    double recurrence_time;
    double assemble_time;
    double base_case_time;
} fast_dixon_level_profile_t;

typedef struct {
    slong total_nvars;
    slong recursive_calls;
    slong base_case_calls;
    slong matrix_inits;
    unsigned long long matrix_entries_allocated;
    unsigned long long poly_mul_ops;
    unsigned long long poly_mul_generated_terms;
    unsigned long long poly_add_ops;
    unsigned long long f_recursive_nodes;
    unsigned long long f_leaf_calls;
    unsigned long long f_valid_leaves;
    unsigned long long f_zero_skips;
    unsigned long long f_pruned_branches;
    unsigned long long cache_hits;
    unsigned long long cache_misses;
    fast_dixon_level_profile_t levels[FAST_DIXON_MAX_LEVELS];
} fast_dixon_profile_t;

static fast_dixon_profile_t g_fast_dixon_profile;

typedef struct {
    slong key_len;
    slong count;
    slong alloc;
    slong slot_count;
    ulong *hashes;
    slong *slots;
    fq_mvpoly_t **keys;
    fast_dixon_matrix_t *values;
    fast_dixon_matrix_support_t *supports;
} fast_dixon_subproblem_cache_t;

static void fast_dixon_matrix_clear(fast_dixon_matrix_t *mat);
static void fast_dixon_matrix_support_init(fast_dixon_matrix_support_t *support,
                                           const fast_dixon_matrix_t *mat);
static void fast_dixon_matrix_support_clear(fast_dixon_matrix_support_t *support);
static void fast_dixon_poly_accumulate_inplace(fq_mvpoly_t *dest,
                                               const fq_mvpoly_t *src,
                                               int subtract);

static void fast_dixon_build_matrix(fast_dixon_matrix_t *out,
                                    const fq_mvpoly_t *const *polys,
                                    const slong *degrees,
                                    slong total_nvars,
                                    slong pos,
                                    slong npars,
                                    const fq_nmod_ctx_t ctx);

static fast_dixon_subproblem_cache_t *g_fast_dixon_subproblem_cache_pool = NULL;
static int g_fast_dixon_subproblem_cache_threads = 0;
static slong g_fast_dixon_trace_logs_by_level[FAST_DIXON_MAX_LEVELS];
static slong g_fast_dixon_trace_total_logs = 0;

static int fast_dixon_profile_enabled(void)
{
    return g_dixon_verbose_level >= 2;
}

static int fast_dixon_profile_heavy_enabled(void)
{
    return g_dixon_verbose_level >= 3;
}

static void fast_dixon_profile_reset(slong total_nvars)
{
    memset(&g_fast_dixon_profile, 0, sizeof(g_fast_dixon_profile));
    g_fast_dixon_profile.total_nvars = total_nvars;
    memset(g_fast_dixon_trace_logs_by_level, 0, sizeof(g_fast_dixon_trace_logs_by_level));
    g_fast_dixon_trace_total_logs = 0;
}

static fast_dixon_level_profile_t *fast_dixon_get_level_profile(slong pos)
{
    if (pos < 0 || pos >= FAST_DIXON_MAX_LEVELS) {
        return NULL;
    }
    return &g_fast_dixon_profile.levels[pos];
}

static void fast_dixon_trace_log(slong depth, const char *fmt, ...)
{
    va_list args;

    if (g_dixon_verbose_level < 3) {
        return;
    }
    if (depth < 0 || depth >= FAST_DIXON_MAX_LEVELS) {
        return;
    }
    if (g_fast_dixon_trace_total_logs >= 24) {
        return;
    }
    if (depth > 0 && g_fast_dixon_trace_logs_by_level[depth] >= 3) {
        return;
    }
    g_fast_dixon_trace_logs_by_level[depth]++;
    g_fast_dixon_trace_total_logs++;

    for (slong i = 0; i < depth * 2; i++) {
        putchar(' ');
    }

    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
}

static void fast_dixon_poly_init_zero(fq_mvpoly_t *poly,
                                      slong nvars,
                                      slong npars,
                                      const fq_nmod_ctx_t ctx)
{
    poly->nvars = nvars;
    poly->npars = npars;
    poly->nterms = 0;
    poly->alloc = 0;
    poly->terms = NULL;
    poly->ctx = ctx;
}

static void fast_dixon_subproblem_cache_init(fast_dixon_subproblem_cache_t *cache,
                                             slong key_len)
{
    cache->key_len = key_len;
    cache->count = 0;
    cache->alloc = 0;
    cache->slot_count = 0;
    cache->hashes = NULL;
    cache->slots = NULL;
    cache->keys = NULL;
    cache->values = NULL;
    cache->supports = NULL;
}

static void fast_dixon_subproblem_cache_clear(fast_dixon_subproblem_cache_t *cache)
{
    if (cache->values != NULL) {
        for (slong i = 0; i < cache->count; i++) {
            fast_dixon_matrix_clear(&cache->values[i]);
        }
        flint_free(cache->values);
    }
    if (cache->supports != NULL) {
        for (slong i = 0; i < cache->count; i++) {
            fast_dixon_matrix_support_clear(&cache->supports[i]);
        }
        flint_free(cache->supports);
    }
    if (cache->keys != NULL) {
        for (slong i = 0; i < cache->count; i++) {
            if (cache->keys[i] != NULL) {
                for (slong j = 0; j < cache->key_len; j++) {
                    fq_mvpoly_clear(&cache->keys[i][j]);
                }
                flint_free(cache->keys[i]);
            }
        }
        flint_free(cache->keys);
    }
    if (cache->hashes != NULL) {
        flint_free(cache->hashes);
    }
    if (cache->slots != NULL) {
        flint_free(cache->slots);
    }

    cache->key_len = 0;
    cache->count = 0;
    cache->alloc = 0;
    cache->slot_count = 0;
    cache->hashes = NULL;
    cache->slots = NULL;
    cache->keys = NULL;
    cache->values = NULL;
    cache->supports = NULL;
}

static ulong fast_dixon_hash_mix(ulong h, ulong x)
{
    h ^= x + 0x9e3779b97f4a7c15UL + (h << 6) + (h >> 2);
    return h;
}

static ulong fast_dixon_fq_nmod_get_ui_prime_field(const fq_nmod_t coeff,
                                                   const fq_nmod_ctx_t ctx)
{
    nmod_poly_t poly;
    ulong value = 0;

    nmod_poly_init(poly, fq_nmod_ctx_prime(ctx));
    fq_nmod_get_nmod_poly(poly, coeff, ctx);
    if (nmod_poly_degree(poly) >= 0) {
        value = nmod_poly_get_coeff_ui(poly, 0);
    }
    nmod_poly_clear(poly);
    return value;
}

static ulong fast_dixon_poly_hash_prime_field(const fq_mvpoly_t *poly)
{
    ulong h = 1469598103934665603UL;

    h = fast_dixon_hash_mix(h, (ulong) poly->nvars);
    h = fast_dixon_hash_mix(h, (ulong) poly->npars);
    h = fast_dixon_hash_mix(h, (ulong) poly->nterms);

    for (slong t = 0; t < poly->nterms; t++) {
        for (slong v = 0; v < poly->nvars; v++) {
            ulong e = (ulong) (poly->terms[t].var_exp ? poly->terms[t].var_exp[v] + 1 : 0);
            h = fast_dixon_hash_mix(h, e);
        }
        for (slong p = 0; p < poly->npars; p++) {
            ulong e = (ulong) (poly->terms[t].par_exp ? poly->terms[t].par_exp[p] + 1 : 0);
            h = fast_dixon_hash_mix(h, e);
        }
        h = fast_dixon_hash_mix(h,
                                fast_dixon_fq_nmod_get_ui_prime_field(poly->terms[t].coeff,
                                                                      poly->ctx));
    }

    return h;
}

static int fast_dixon_poly_equal_prime_field(const fq_mvpoly_t *lhs,
                                             const fq_mvpoly_t *rhs)
{
    if (lhs->nvars != rhs->nvars || lhs->npars != rhs->npars || lhs->nterms != rhs->nterms) {
        return 0;
    }

    for (slong t = 0; t < lhs->nterms; t++) {
        for (slong v = 0; v < lhs->nvars; v++) {
            slong lhs_exp = lhs->terms[t].var_exp ? lhs->terms[t].var_exp[v] : 0;
            slong rhs_exp = rhs->terms[t].var_exp ? rhs->terms[t].var_exp[v] : 0;
            if (lhs_exp != rhs_exp) {
                return 0;
            }
        }
        for (slong p = 0; p < lhs->npars; p++) {
            slong lhs_exp = lhs->terms[t].par_exp ? lhs->terms[t].par_exp[p] : 0;
            slong rhs_exp = rhs->terms[t].par_exp ? rhs->terms[t].par_exp[p] : 0;
            if (lhs_exp != rhs_exp) {
                return 0;
            }
        }
        if (!fq_nmod_equal(lhs->terms[t].coeff, rhs->terms[t].coeff, lhs->ctx)) {
            return 0;
        }
    }

    return 1;
}

static ulong fast_dixon_subproblem_cache_key_hash(const fq_mvpoly_t *const *key,
                                                  slong key_len)
{
    ulong h = 0xcbf29ce484222325UL;

    for (slong i = 0; i < key_len; i++) {
        h = fast_dixon_hash_mix(h, fast_dixon_poly_hash_prime_field(key[i]));
    }

    return h;
}

static int fast_dixon_subproblem_cache_key_equal(const fq_mvpoly_t *const *lhs,
                                                 const fq_mvpoly_t *rhs,
                                                 slong key_len)
{
    for (slong i = 0; i < key_len; i++) {
        if (!fast_dixon_poly_equal_prime_field(lhs[i], &rhs[i])) {
            return 0;
        }
    }
    return 1;
}

static ulong fast_dixon_subproblem_cache_key_hash_from_parts(const ulong *parts,
                                                             slong key_len)
{
    ulong h = 0xcbf29ce484222325UL;

    for (slong i = 0; i < key_len; i++) {
        h = fast_dixon_hash_mix(h, parts[i]);
    }

    return h;
}

static void fast_dixon_subproblem_cache_reset_all(void)
{
    g_fast_dixon_subproblem_cache_pool = NULL;
    g_fast_dixon_subproblem_cache_threads = 0;
}

static void fast_dixon_subproblem_cache_clear_all(void)
{
    if (g_fast_dixon_subproblem_cache_pool == NULL) {
        return;
    }

    for (slong i = 0; i < (slong) g_fast_dixon_subproblem_cache_threads * FAST_DIXON_MAX_LEVELS; i++) {
        fast_dixon_subproblem_cache_clear(&g_fast_dixon_subproblem_cache_pool[i]);
    }

    flint_free(g_fast_dixon_subproblem_cache_pool);
    g_fast_dixon_subproblem_cache_pool = NULL;
    g_fast_dixon_subproblem_cache_threads = 0;
}

static void fast_dixon_subproblem_cache_prepare_pool(int thread_count)
{
    if (thread_count < 1) {
        thread_count = 1;
    }

    g_fast_dixon_subproblem_cache_pool = (fast_dixon_subproblem_cache_t *)
        flint_calloc((size_t) thread_count * FAST_DIXON_MAX_LEVELS,
                     sizeof(fast_dixon_subproblem_cache_t));
    g_fast_dixon_subproblem_cache_threads = thread_count;
}

static void fast_dixon_subproblem_cache_rebuild_slots(fast_dixon_subproblem_cache_t *cache,
                                                      slong slot_count)
{
    cache->slots = (slong *) flint_realloc(cache->slots, (size_t) slot_count * sizeof(slong));
    memset(cache->slots, 0, (size_t) slot_count * sizeof(slong));
    cache->slot_count = slot_count;

    for (slong i = 0; i < cache->count; i++) {
        slong slot = (slong) (cache->hashes[i] & (ulong) (cache->slot_count - 1));
        while (cache->slots[slot] != 0) {
            slot = (slot + 1) & (cache->slot_count - 1);
        }
        cache->slots[slot] = i + 1;
    }
}

static void fast_dixon_subproblem_cache_reserve(fast_dixon_subproblem_cache_t *cache,
                                                slong new_alloc)
{
    cache->hashes = (ulong *) flint_realloc(cache->hashes, (size_t) new_alloc * sizeof(ulong));
    cache->keys = (fq_mvpoly_t **) flint_realloc(cache->keys, (size_t) new_alloc * sizeof(fq_mvpoly_t *));
    cache->values = (fast_dixon_matrix_t *) flint_realloc(cache->values,
                                                          (size_t) new_alloc * sizeof(fast_dixon_matrix_t));
    cache->supports = (fast_dixon_matrix_support_t *) flint_realloc(cache->supports,
                                                                    (size_t) new_alloc * sizeof(fast_dixon_matrix_support_t));
    for (slong i = cache->alloc; i < new_alloc; i++) {
        cache->keys[i] = NULL;
        cache->supports[i].rows = 0;
        cache->supports[i].total_nonzero = 0;
        cache->supports[i].row_counts = NULL;
        cache->supports[i].row_offsets = NULL;
        cache->supports[i].col_indices = NULL;
    }
    cache->alloc = new_alloc;

    if (cache->slot_count < 2 * cache->alloc) {
        slong slot_count = 16;
        while (slot_count < 2 * cache->alloc) {
            slot_count <<= 1;
        }
        fast_dixon_subproblem_cache_rebuild_slots(cache, slot_count);
    }
}

static slong fast_dixon_subproblem_cache_lookup(const fast_dixon_subproblem_cache_t *cache,
                                                const fq_mvpoly_t *const *selected,
                                                ulong key_hash)
{
    slong slot;

    if (cache == NULL || cache->slot_count == 0) {
        return -1;
    }

    slot = (slong) (key_hash & (ulong) (cache->slot_count - 1));
    while (cache->slots[slot] != 0) {
        slong idx = cache->slots[slot] - 1;
        if (cache->hashes[idx] == key_hash &&
            fast_dixon_subproblem_cache_key_equal(selected, cache->keys[idx], cache->key_len)) {
            return idx;
        }
        slot = (slot + 1) & (cache->slot_count - 1);
    }

    return -1;
}

static slong fast_dixon_subproblem_cache_insert(fast_dixon_subproblem_cache_t *cache,
                                                const fq_mvpoly_t *const *selected,
                                                ulong key_hash)
{
    slong new_index;
    slong slot;

    if (cache->count >= cache->alloc) {
        fast_dixon_subproblem_cache_reserve(cache, cache->alloc ? 2 * cache->alloc : 16);
    }

    new_index = cache->count++;
    cache->hashes[new_index] = key_hash;
    cache->keys[new_index] =
        (fq_mvpoly_t *) flint_malloc((size_t) cache->key_len * sizeof(fq_mvpoly_t));
    for (slong i = 0; i < cache->key_len; i++) {
        fast_dixon_poly_init_zero(&cache->keys[new_index][i],
                                  selected[i]->nvars,
                                  selected[i]->npars,
                                  selected[i]->ctx);
        fq_mvpoly_copy(&cache->keys[new_index][i], selected[i]);
    }

    if (4 * cache->count >= 3 * cache->slot_count) {
        slong slot_count = cache->slot_count ? cache->slot_count : 16;
        while (4 * cache->count >= 3 * slot_count) {
            slot_count <<= 1;
        }
        fast_dixon_subproblem_cache_rebuild_slots(cache, slot_count);
    }

    slot = (slong) (key_hash & (ulong) (cache->slot_count - 1));
    while (cache->slots[slot] != 0) {
        slot = (slot + 1) & (cache->slot_count - 1);
    }
    cache->slots[slot] = new_index + 1;

    return new_index;
}

static fast_dixon_subproblem_cache_t *fast_dixon_get_shared_subproblem_cache(slong child_pos,
                                                                             slong total_nvars,
                                                                             const fq_nmod_ctx_t ctx)
{
    fast_dixon_subproblem_cache_t *cache;
    slong key_len;
    int tid = 0;

    if (fq_nmod_ctx_degree(ctx) != 1) {
        return NULL;
    }
    if (child_pos != total_nvars - 1) {
        return NULL;
    }
    if (child_pos < 0 || child_pos >= FAST_DIXON_MAX_LEVELS) {
        return NULL;
    }
    if (g_fast_dixon_subproblem_cache_pool == NULL || g_fast_dixon_subproblem_cache_threads <= 0) {
        return NULL;
    }

#ifdef _OPENMP
    if (omp_in_parallel()) {
        return NULL;
    }
#endif
    if (tid < 0 || tid >= g_fast_dixon_subproblem_cache_threads) {
        tid = 0;
    }

    key_len = total_nvars - child_pos + 1;
    cache = &g_fast_dixon_subproblem_cache_pool[tid * FAST_DIXON_MAX_LEVELS + child_pos];
    if (cache->key_len == 0 && cache->count == 0 && cache->alloc == 0 &&
        cache->slot_count == 0 && cache->hashes == NULL && cache->slots == NULL &&
        cache->keys == NULL && cache->values == NULL) {
        fast_dixon_subproblem_cache_init(cache, key_len);
    }

    return cache;
}

static int fast_dixon_parallel_threads_for_level(slong pos,
                                                 slong current_degree,
                                                 slong num_f_blocks)
{
#ifdef _OPENMP
    slong work_items = current_degree * num_f_blocks;
    int max_threads;
    int capped_threads;

    if (g_dixon_verbose_level >= 3) {
        return 0;
    }

    if (pos != 0) {
        return 0;
    }

    if (omp_in_parallel()) {
        return 0;
    }

    if (work_items < 8) {
        return 0;
    }

    max_threads = omp_get_max_threads();
    if (max_threads <= 1) {
        return 0;
    }

    capped_threads = max_threads;
    if (work_items >= 256) {
        if (capped_threads > 16) {
            capped_threads = 16;
        }
    } else if (capped_threads > 12) {
        capped_threads = 12;
    }
    if (capped_threads > (int) work_items) {
        capped_threads = (int) work_items;
    }

    return capped_threads > 1 ? capped_threads : 0;
#else
    (void) pos;
    (void) current_degree;
    (void) num_f_blocks;
    return 0;
#endif
}

static const char *fast_dixon_det_method_name(det_method_t method)
{
    switch (method) {
        case DET_METHOD_RECURSIVE:
            return "recursive expansion";
        case DET_METHOD_KRONECKER:
            return "HNF";
        case DET_METHOD_INTERPOLATION:
            return "interpolation";
        case DET_METHOD_HUANG:
            return "sparse interpolation";
        case DET_METHOD_KRONECKER_NMOD:
            return "Bareiss";
        default:
            return "default";
    }
}

static void fast_dixon_info_log(const char *fmt, ...)
{
    va_list args;

    if (g_dixon_verbose_level < 1) {
        return;
    }

    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
}

static slong fast_dixon_factorial(slong n)
{
    slong value = 1;

    for (slong i = 2; i <= n; i++) {
        value *= i;
    }

    return value;
}

static slong fast_dixon_level_size(const slong *degrees,
                                   slong total_nvars,
                                   slong pos)
{
    slong remaining = total_nvars - pos;
    slong size = 1;

    if (remaining <= 0) {
        return 1;
    }

    size = fast_dixon_factorial(remaining);
    for (slong i = pos; i < total_nvars; i++) {
        size *= degrees[i];
    }

    return size;
}

static void fast_dixon_matrix_init(fast_dixon_matrix_t *mat,
                                   slong rows,
                                   slong cols,
                                   slong npars,
                                   const fq_nmod_ctx_t ctx)
{
    mat->rows = rows;
    mat->cols = cols;
    mat->npars = npars;
    mat->ctx = ctx;
    mat->entries = (fq_mvpoly_t *) flint_malloc((size_t) (rows * cols) * sizeof(fq_mvpoly_t));

    for (slong i = 0; i < rows * cols; i++) {
        fast_dixon_poly_init_zero(&mat->entries[i], 0, npars, ctx);
    }

    if (fast_dixon_profile_heavy_enabled()) {
        g_fast_dixon_profile.matrix_inits++;
        g_fast_dixon_profile.matrix_entries_allocated +=
            (unsigned long long) rows * (unsigned long long) cols;
    }
}

static void fast_dixon_matrix_clear(fast_dixon_matrix_t *mat)
{
    if (mat->entries == NULL) {
        return;
    }

    for (slong i = 0; i < mat->rows * mat->cols; i++) {
        fq_mvpoly_clear(&mat->entries[i]);
    }

    flint_free(mat->entries);
    mat->entries = NULL;
    mat->rows = 0;
    mat->cols = 0;
    mat->npars = 0;
    mat->ctx = NULL;
}

static void fast_dixon_matrix_add_inplace(fast_dixon_matrix_t *dest,
                                          const fast_dixon_matrix_t *src,
                                          int subtract)
{
    for (slong i = 0; i < dest->rows; i++) {
        for (slong j = 0; j < dest->cols; j++) {
            fq_mvpoly_t *dest_entry = &FAST_DIXON_ENTRY(dest, i, j);
            const fq_mvpoly_t *src_entry = &FAST_DIXON_ENTRY(src, i, j);

            if (src_entry->nterms == 0) {
                continue;
            }

            fast_dixon_poly_accumulate_inplace(dest_entry, src_entry, subtract);

            if (fast_dixon_profile_heavy_enabled()) {
                g_fast_dixon_profile.poly_add_ops++;
            }
        }
    }
}

static void fast_dixon_poly_assign_copy(fq_mvpoly_t *dest, const fq_mvpoly_t *src)
{
    if (dest->terms != NULL) {
        fq_mvpoly_clear(dest);
    }
    fq_mvpoly_copy(dest, src);
}

static void fast_dixon_poly_assign_neg_copy(fq_mvpoly_t *dest, const fq_mvpoly_t *src)
{
    if (dest->terms != NULL) {
        fq_mvpoly_clear(dest);
    }

    fq_mvpoly_init(dest, src->nvars, src->npars, src->ctx);
    for (slong t = 0; t < src->nterms; t++) {
        fq_nmod_t neg_coeff;
        fq_nmod_init(neg_coeff, src->ctx);
        fq_nmod_neg(neg_coeff, src->terms[t].coeff, src->ctx);
        fq_mvpoly_add_term_fast(dest, src->terms[t].var_exp, src->terms[t].par_exp, neg_coeff);
        fq_nmod_clear(neg_coeff, src->ctx);
    }
}

static int fast_dixon_poly_use_direct_accumulate(const fq_mvpoly_t *dest,
                                                 const fq_mvpoly_t *src)
{
    slong sum_terms;
    slong prod_terms;

    if (dest->ctx != src->ctx ||
        dest->nvars != src->nvars ||
        dest->npars != src->npars ||
        fq_nmod_ctx_degree(dest->ctx) != 1) {
        return 0;
    }

    sum_terms = dest->nterms + src->nterms;
    prod_terms = dest->nterms * src->nterms;
    return (sum_terms <= 48 || prod_terms <= 256);
}

static void fast_dixon_poly_accumulate_inplace(fq_mvpoly_t *dest,
                                               const fq_mvpoly_t *src,
                                               int subtract)
{
    if (src->nterms == 0) {
        return;
    }

    if (dest->nterms == 0) {
        if (subtract) {
            fast_dixon_poly_assign_neg_copy(dest, src);
        } else {
            fast_dixon_poly_assign_copy(dest, src);
        }
        return;
    }

    if (fast_dixon_poly_use_direct_accumulate(dest, src)) {
        for (slong t = 0; t < src->nterms; t++) {
            if (subtract) {
                fq_nmod_t neg_coeff;
                fq_nmod_init(neg_coeff, src->ctx);
                fq_nmod_neg(neg_coeff, src->terms[t].coeff, src->ctx);
                fq_mvpoly_add_term(dest,
                                   src->terms[t].var_exp,
                                   src->terms[t].par_exp,
                                   neg_coeff);
                fq_nmod_clear(neg_coeff, src->ctx);
            } else {
                fq_mvpoly_add_term(dest,
                                   src->terms[t].var_exp,
                                   src->terms[t].par_exp,
                                   src->terms[t].coeff);
            }
        }
        return;
    }

    if (subtract) {
        fq_mvpoly_sub(dest, dest, src);
    } else {
        fq_mvpoly_add(dest, dest, src);
    }
}

static void fast_dixon_matrix_support_init(fast_dixon_matrix_support_t *support,
                                           const fast_dixon_matrix_t *mat)
{
    slong total_nonzero = 0;

    support->rows = mat->rows;
    support->total_nonzero = 0;
    support->row_counts = (slong *) flint_calloc((size_t) mat->rows, sizeof(slong));
    support->row_offsets = (slong *) flint_calloc((size_t) (mat->rows + 1), sizeof(slong));
    support->col_indices = NULL;

    for (slong i = 0; i < mat->rows; i++) {
        slong count = 0;
        for (slong j = 0; j < mat->cols; j++) {
            if (FAST_DIXON_ENTRY(mat, i, j).nterms > 0) {
                count++;
            }
        }
        support->row_counts[i] = count;
        support->row_offsets[i] = total_nonzero;
        total_nonzero += count;
    }
    support->row_offsets[mat->rows] = total_nonzero;
    support->total_nonzero = total_nonzero;

    if (total_nonzero > 0) {
        support->col_indices = (slong *) flint_malloc((size_t) total_nonzero * sizeof(slong));
        total_nonzero = 0;
        for (slong i = 0; i < mat->rows; i++) {
            for (slong j = 0; j < mat->cols; j++) {
                if (FAST_DIXON_ENTRY(mat, i, j).nterms > 0) {
                    support->col_indices[total_nonzero++] = j;
                }
            }
        }
    }
}

static void fast_dixon_matrix_support_clear(fast_dixon_matrix_support_t *support)
{
    if (support->row_counts != NULL) {
        flint_free(support->row_counts);
    }
    if (support->row_offsets != NULL) {
        flint_free(support->row_offsets);
    }
    if (support->col_indices != NULL) {
        flint_free(support->col_indices);
    }

    support->rows = 0;
    support->total_nonzero = 0;
    support->row_counts = NULL;
    support->row_offsets = NULL;
    support->col_indices = NULL;
}

static void fast_dixon_matrix_mul_accumulate(fast_dixon_matrix_t *dest,
                                             const fast_dixon_matrix_t *left,
                                             const fast_dixon_matrix_t *right,
                                             int subtract)
{
    for (slong i = 0; i < left->rows; i++) {
        for (slong k = 0; k < left->cols; k++) {
            const fq_mvpoly_t *left_entry = &FAST_DIXON_ENTRY(left, i, k);

            if (left_entry->nterms == 0) {
                continue;
            }

            for (slong j = 0; j < right->cols; j++) {
                const fq_mvpoly_t *right_entry = &FAST_DIXON_ENTRY(right, k, j);
                fq_mvpoly_t product;

                if (right_entry->nterms == 0) {
                    continue;
                }

                if (fast_dixon_profile_heavy_enabled()) {
                    g_fast_dixon_profile.poly_mul_ops++;
                }
                fq_mvpoly_mul(&product, left_entry, right_entry);
                if (fast_dixon_profile_heavy_enabled()) {
                    g_fast_dixon_profile.poly_mul_generated_terms +=
                        (unsigned long long) product.nterms;
                }
                if (subtract) {
                    fast_dixon_poly_accumulate_inplace(&FAST_DIXON_ENTRY(dest, i, j),
                                                       &product,
                                                       1);
                } else {
                    fast_dixon_poly_accumulate_inplace(&FAST_DIXON_ENTRY(dest, i, j),
                                                       &product,
                                                       0);
                }
                if (fast_dixon_profile_heavy_enabled()) {
                    g_fast_dixon_profile.poly_add_ops++;
                }
                fq_mvpoly_clear(&product);
            }
        }
    }
}

static void fast_dixon_matrix_mul_accumulate_supported(fast_dixon_matrix_t *dest,
                                                       const fast_dixon_matrix_t *left,
                                                       const fast_dixon_matrix_support_t *left_support,
                                                       const fast_dixon_matrix_t *right,
                                                       const fast_dixon_matrix_support_t *right_support,
                                                       int subtract)
{
    if (left_support == NULL || right_support == NULL) {
        fast_dixon_matrix_mul_accumulate(dest, left, right, subtract);
        return;
    }

    for (slong i = 0; i < left->rows; i++) {
        slong left_count = left_support->row_counts[i];
        slong left_offset = left_support->row_offsets[i];

        for (slong left_idx = 0; left_idx < left_count; left_idx++) {
            slong k = left_support->col_indices[left_offset + left_idx];
            const fq_mvpoly_t *left_entry = &FAST_DIXON_ENTRY(left, i, k);
            slong right_count = right_support->row_counts[k];
            slong right_offset = right_support->row_offsets[k];

            for (slong right_idx = 0; right_idx < right_count; right_idx++) {
                slong j = right_support->col_indices[right_offset + right_idx];
                const fq_mvpoly_t *right_entry = &FAST_DIXON_ENTRY(right, k, j);
                fq_mvpoly_t product;

                if (fast_dixon_profile_heavy_enabled()) {
                    g_fast_dixon_profile.poly_mul_ops++;
                }
                fq_mvpoly_mul(&product, left_entry, right_entry);
                if (fast_dixon_profile_heavy_enabled()) {
                    g_fast_dixon_profile.poly_mul_generated_terms +=
                        (unsigned long long) product.nterms;
                }
                if (subtract) {
                    fast_dixon_poly_accumulate_inplace(&FAST_DIXON_ENTRY(dest, i, j),
                                                       &product,
                                                       1);
                } else {
                    fast_dixon_poly_accumulate_inplace(&FAST_DIXON_ENTRY(dest, i, j),
                                                       &product,
                                                       0);
                }
                if (fast_dixon_profile_heavy_enabled()) {
                    g_fast_dixon_profile.poly_add_ops++;
                }
                fq_mvpoly_clear(&product);
            }
        }
    }
}

static void fast_dixon_matrix_copy_block(fast_dixon_matrix_t *dest,
                                         slong row_offset,
                                         slong col_offset,
                                         const fast_dixon_matrix_t *src)
{
    for (slong i = 0; i < src->rows; i++) {
        for (slong j = 0; j < src->cols; j++) {
            if (FAST_DIXON_ENTRY(src, i, j).nterms == 0) {
                continue;
            }
            fq_mvpoly_copy(&FAST_DIXON_ENTRY(dest, row_offset + i, col_offset + j),
                           &FAST_DIXON_ENTRY(src, i, j));
        }
    }
}

static void fast_dixon_compute_degree_bounds(slong *degrees,
                                             const fq_mvpoly_t *const *polys,
                                             slong npolys,
                                             slong nvars)
{
    for (slong v = 0; v < nvars; v++) {
        slong max_deg = 0;

        for (slong p = 0; p < npolys; p++) {
            for (slong t = 0; t < polys[p]->nterms; t++) {
                slong deg = polys[p]->terms[t].var_exp ? polys[p]->terms[t].var_exp[v] : 0;
                if (deg > max_deg) {
                    max_deg = deg;
                }
            }
        }

        degrees[v] = max_deg;
    }
}

static fq_mvpoly_t **fast_dixon_split_polys_by_degree(const fq_mvpoly_t *const *polys,
                                                      slong npolys,
                                                      slong current_degree)
{
    fq_mvpoly_t **coeffs = (fq_mvpoly_t **) flint_malloc((size_t) npolys * sizeof(fq_mvpoly_t *));

    for (slong p = 0; p < npolys; p++) {
        slong remaining_vars = polys[p]->nvars > 0 ? polys[p]->nvars - 1 : 0;

        coeffs[p] = (fq_mvpoly_t *) flint_malloc((size_t) (current_degree + 1) * sizeof(fq_mvpoly_t));
        for (slong d = 0; d <= current_degree; d++) {
            fq_mvpoly_init(&coeffs[p][d], remaining_vars, polys[p]->npars, polys[p]->ctx);
        }

        for (slong t = 0; t < polys[p]->nterms; t++) {
            slong deg = polys[p]->terms[t].var_exp ? polys[p]->terms[t].var_exp[0] : 0;

            if (deg < 0 || deg > current_degree) {
                continue;
            }

            if (remaining_vars > 0) {
                slong *new_var_exp = (slong *) flint_calloc((size_t) remaining_vars, sizeof(slong));
                if (polys[p]->terms[t].var_exp != NULL) {
                    memcpy(new_var_exp,
                           polys[p]->terms[t].var_exp + 1,
                           (size_t) remaining_vars * sizeof(slong));
                }
                fq_mvpoly_add_term(&coeffs[p][deg],
                                   new_var_exp,
                                   polys[p]->terms[t].par_exp,
                                   polys[p]->terms[t].coeff);
                flint_free(new_var_exp);
            } else {
                fq_mvpoly_add_term(&coeffs[p][deg],
                                   NULL,
                                   polys[p]->terms[t].par_exp,
                                   polys[p]->terms[t].coeff);
            }
        }
    }

    return coeffs;
}

static ulong **fast_dixon_build_split_hashes(fq_mvpoly_t **coeffs,
                                             slong npolys,
                                             slong current_degree,
                                             const fq_nmod_ctx_t ctx)
{
    ulong **hashes;

    if (fq_nmod_ctx_degree(ctx) != 1) {
        return NULL;
    }

    hashes = (ulong **) flint_malloc((size_t) npolys * sizeof(ulong *));
    for (slong p = 0; p < npolys; p++) {
        hashes[p] = (ulong *) flint_malloc((size_t) (current_degree + 1) * sizeof(ulong));
        for (slong d = 0; d <= current_degree; d++) {
            hashes[p][d] = fast_dixon_poly_hash_prime_field(&coeffs[p][d]);
        }
    }

    return hashes;
}

static void fast_dixon_clear_split_coeffs(fq_mvpoly_t **coeffs,
                                          slong npolys,
                                          slong current_degree)
{
    if (coeffs == NULL) {
        return;
    }

    for (slong p = 0; p < npolys; p++) {
        for (slong d = 0; d <= current_degree; d++) {
            fq_mvpoly_clear(&coeffs[p][d]);
        }
        flint_free(coeffs[p]);
    }

    flint_free(coeffs);
}

static void fast_dixon_clear_split_hashes(ulong **hashes,
                                          slong npolys)
{
    if (hashes == NULL) {
        return;
    }

    for (slong p = 0; p < npolys; p++) {
        flint_free(hashes[p]);
    }

    flint_free(hashes);
}

static void fast_dixon_decode_rectangular_index(slong index,
                                                const slong *counts,
                                                slong nvars,
                                                slong *exp_out)
{
    for (slong i = 0; i < nvars; i++) {
        exp_out[i] = 0;
    }

    /*
     * Keep the monomial order compatible with the standard Dixon code and the
     * Maple reference: x_{pos+1} > x_{pos+2} > ... > x_{n}, so the last
     * remaining variable varies fastest.
     */
    for (slong i = nvars - 1; i >= 0; i--) {
        exp_out[i] = counts[i] > 0 ? index % counts[i] : 0;
        if (counts[i] > 0) {
            index /= counts[i];
        }
    }
}

static slong *fast_dixon_build_multiplier_table(const slong *degrees,
                                                slong total_nvars,
                                                slong pos,
                                                slong *count_out,
                                                slong *nvars_out)
{
    slong rem_vars = total_nvars - pos - 1;
    slong count = fast_dixon_level_size(degrees, total_nvars, pos + 1);
    slong *table;
    slong *counts;

    *count_out = count;
    *nvars_out = rem_vars;

    if (rem_vars <= 0) {
        table = (slong *) flint_malloc(sizeof(slong));
        table[0] = 0;
        return table;
    }

    counts = (slong *) flint_malloc((size_t) rem_vars * sizeof(slong));
    for (slong i = 0; i < rem_vars; i++) {
        counts[i] = (i + 1) * degrees[pos + 1 + i];
    }

    table = (slong *) flint_calloc((size_t) count * rem_vars, sizeof(slong));
    for (slong idx = 0; idx < count; idx++) {
        fast_dixon_decode_rectangular_index(idx,
                                            counts,
                                            rem_vars,
                                            table + idx * rem_vars);
    }

    flint_free(counts);
    return table;
}

static void fast_dixon_build_p_block(fast_dixon_matrix_t *out,
                                     const fq_mvpoly_t *coeff_poly,
                                     const slong *degrees,
                                     slong total_nvars,
                                     slong pos,
                                     const slong *mult_table,
                                     const slong *row_counts,
                                     slong lower_size,
                                     slong rem_vars)
{
    slong current_vars = total_nvars - pos;
    slong syl_rows = current_vars * lower_size;

    fast_dixon_matrix_init(out, syl_rows, lower_size, coeff_poly->npars, coeff_poly->ctx);
    if (coeff_poly->nterms == 0) {
        return;
    }

    for (slong col = 0; col < lower_size; col++) {
        const slong *mult_exp = rem_vars > 0 ? mult_table + col * rem_vars : NULL;

        for (slong t = 0; t < coeff_poly->nterms; t++) {
            slong row_idx = 0;
            slong stride = 1;

            for (slong v = rem_vars - 1; v >= 0; v--) {
                slong term_exp = coeff_poly->terms[t].var_exp ? coeff_poly->terms[t].var_exp[v] : 0;
                row_idx += (mult_exp[v] + term_exp) * stride;
                stride *= row_counts[v];
            }

            fq_mvpoly_add_term(&FAST_DIXON_ENTRY(out, row_idx, col),
                               NULL,
                               coeff_poly->terms[t].par_exp,
                               coeff_poly->terms[t].coeff);
        }
    }
}

static void fast_dixon_build_base_case(fast_dixon_matrix_t *out,
                                       fq_mvpoly_t **coeffs,
                                       slong current_degree,
                                       slong npars,
                                       const fq_nmod_ctx_t ctx)
{
    fast_dixon_matrix_init(out, current_degree, current_degree, npars, ctx);

    for (slong high = 1; high <= current_degree; high++) {
        for (slong low = 0; low < high; low++) {
            fq_mvpoly_t a_high_b_low;
            fq_mvpoly_t b_high_a_low;
            fq_mvpoly_t coeff;
            slong span = high - low;

            coeff.terms = NULL;
            coeff.nterms = 0;
            coeff.alloc = 0;
            coeff.nvars = 0;
            coeff.npars = 0;
            coeff.ctx = ctx;

            fq_mvpoly_mul(&a_high_b_low, &coeffs[0][high], &coeffs[1][low]);
            fq_mvpoly_mul(&b_high_a_low, &coeffs[1][high], &coeffs[0][low]);
            fq_mvpoly_sub(&coeff, &a_high_b_low, &b_high_a_low);

            fq_mvpoly_clear(&a_high_b_low);
            fq_mvpoly_clear(&b_high_a_low);

            if (coeff.nterms > 0) {
                for (slong t = 0; t < span; t++) {
                    fq_mvpoly_add(&FAST_DIXON_ENTRY(out, high - 1 - t, low + t),
                                  &FAST_DIXON_ENTRY(out, high - 1 - t, low + t),
                                  &coeff);
                }
            }

            fq_mvpoly_clear(&coeff);
        }
    }
}

static void fast_dixon_matrix_add_interleaved_rows(fast_dixon_matrix_t *dest,
                                                   const fast_dixon_matrix_t *src,
                                                   const fast_dixon_matrix_support_t *src_support,
                                                   slong interleave_index,
                                                   slong interleave_width,
                                                   int subtract)
{
    if (src_support != NULL) {
        for (slong i = 0; i < src->rows; i++) {
            slong dest_row = i * interleave_width + interleave_index;
            slong row_count = src_support->row_counts[i];
            slong row_offset = src_support->row_offsets[i];

            for (slong idx = 0; idx < row_count; idx++) {
                slong j = src_support->col_indices[row_offset + idx];
                fq_mvpoly_t *dest_entry = &FAST_DIXON_ENTRY(dest, dest_row, j);
                const fq_mvpoly_t *src_entry = &FAST_DIXON_ENTRY(src, i, j);

                fast_dixon_poly_accumulate_inplace(dest_entry, src_entry, subtract);
                if (fast_dixon_profile_heavy_enabled()) {
                    g_fast_dixon_profile.poly_add_ops++;
                }
            }
        }
        return;
    }

    for (slong i = 0; i < src->rows; i++) {
        slong dest_row = i * interleave_width + interleave_index;

        for (slong j = 0; j < src->cols; j++) {
            fq_mvpoly_t *dest_entry = &FAST_DIXON_ENTRY(dest, dest_row, j);
            const fq_mvpoly_t *src_entry = &FAST_DIXON_ENTRY(src, i, j);

            if (src_entry->nterms == 0) {
                continue;
            }

            fast_dixon_poly_accumulate_inplace(dest_entry, src_entry, subtract);
            if (fast_dixon_profile_heavy_enabled()) {
                g_fast_dixon_profile.poly_add_ops++;
            }
        }
    }
}

static void fast_dixon_matrix_scatter_columns(fast_dixon_matrix_t *dest,
                                              const fast_dixon_matrix_t *src,
                                              slong scatter_index,
                                              slong scatter_width)
{
    for (slong i = 0; i < src->rows; i++) {
        for (slong j = 0; j < src->cols; j++) {
            if (FAST_DIXON_ENTRY(src, i, j).nterms == 0) {
                continue;
            }
            fq_mvpoly_copy(&FAST_DIXON_ENTRY(dest, i, j * scatter_width + scatter_index),
                           &FAST_DIXON_ENTRY(src, i, j));
        }
    }
}

static void fast_dixon_accumulate_f_recursive(fast_dixon_matrix_t *f_block,
                                              fq_mvpoly_t **coeffs,
                                              ulong **coeff_hashes,
                                              const slong *poly_indices,
                                              const slong *degrees,
                                              slong total_nvars,
                                              slong pos,
                                              slong npars,
                                              const fq_nmod_ctx_t ctx,
                                              slong depth,
                                              slong tuple_len,
                                              slong target_sum,
                                              slong current_sum,
                                              int has_zero,
                                              const fq_mvpoly_t **selected,
                                              ulong *selected_hashes,
                                              fast_dixon_subproblem_cache_t *cache,
                                              slong omit_idx,
                                              slong npolys)
{
    slong current_degree = degrees[pos];
    fast_dixon_level_profile_t *level_profile =
        fast_dixon_profile_heavy_enabled() ? fast_dixon_get_level_profile(pos) : NULL;

    if (fast_dixon_profile_heavy_enabled()) {
        g_fast_dixon_profile.f_recursive_nodes++;
        if (level_profile != NULL) {
            level_profile->f_recursive_nodes++;
        }
    }

    if (depth == tuple_len) {
        const fast_dixon_matrix_t *submatrix;
        const fast_dixon_matrix_support_t *submatrix_support = NULL;
        fast_dixon_matrix_t uncached_submatrix;
        fast_dixon_matrix_support_t uncached_support;
        double child_start;
        double child_elapsed;
        slong cached_index = -1;
        ulong key_hash = 0;

        if (fast_dixon_profile_heavy_enabled()) {
            g_fast_dixon_profile.f_leaf_calls++;
            if (level_profile != NULL) {
                level_profile->f_leaf_calls++;
            }
        }

        if (has_zero) {
            if (fast_dixon_profile_heavy_enabled()) {
                g_fast_dixon_profile.f_zero_skips++;
                if (level_profile != NULL) {
                    level_profile->f_zero_skips++;
                }
            }
            return;
        }

        if (current_sum != target_sum) {
            return;
        }

        if (fast_dixon_profile_heavy_enabled()) {
            g_fast_dixon_profile.f_valid_leaves++;
            if (level_profile != NULL) {
                level_profile->f_valid_leaves++;
            }
        }

        submatrix = NULL;
        if (cache != NULL) {
            key_hash = fast_dixon_subproblem_cache_key_hash_from_parts(selected_hashes,
                                                                       cache->key_len);
            cached_index = fast_dixon_subproblem_cache_lookup(cache, selected, key_hash);
            if (cached_index >= 0) {
                submatrix = &cache->values[cached_index];
                submatrix_support = &cache->supports[cached_index];
                if (fast_dixon_profile_heavy_enabled()) {
                    g_fast_dixon_profile.cache_hits++;
                    if (level_profile != NULL) {
                        level_profile->cache_hits++;
                    }
                }
            } else {
                slong new_index = fast_dixon_subproblem_cache_insert(cache, selected, key_hash);

                child_start = get_wall_time();
                fast_dixon_build_matrix(&cache->values[new_index],
                                        selected,
                                        degrees,
                                        total_nvars,
                                        pos + 1,
                                        npars,
                                        ctx);
                child_elapsed = get_wall_time() - child_start;
                submatrix = &cache->values[new_index];
                fast_dixon_matrix_support_init(&cache->supports[new_index], submatrix);
                submatrix_support = &cache->supports[new_index];
                if (fast_dixon_profile_heavy_enabled()) {
                    g_fast_dixon_profile.cache_misses++;
                    if (level_profile != NULL) {
                        level_profile->cache_misses++;
                        level_profile->f_child_build_time += child_elapsed;
                    }
                }
            }
        }

        if (submatrix == NULL) {
            uncached_support.rows = 0;
            uncached_support.total_nonzero = 0;
            uncached_support.row_counts = NULL;
            uncached_support.row_offsets = NULL;
            uncached_support.col_indices = NULL;
            child_start = get_wall_time();
            fast_dixon_build_matrix(&uncached_submatrix,
                                    selected,
                                    degrees,
                                    total_nvars,
                                    pos + 1,
                                    npars,
                                    ctx);
            submatrix = &uncached_submatrix;
            fast_dixon_matrix_support_init(&uncached_support, submatrix);
            submatrix_support = &uncached_support;
            child_elapsed = get_wall_time() - child_start;
            if (fast_dixon_profile_heavy_enabled()) {
                g_fast_dixon_profile.cache_misses++;
                if (level_profile != NULL) {
                    level_profile->cache_misses++;
                    level_profile->f_child_build_time += child_elapsed;
                }
            }
            fast_dixon_matrix_add_interleaved_rows(f_block,
                                                   submatrix,
                                                   submatrix_support,
                                                   omit_idx,
                                                   npolys,
                                                   (omit_idx & 1) != 0);
            fast_dixon_matrix_support_clear(&uncached_support);
            fast_dixon_matrix_clear(&uncached_submatrix);
            return;
        }

        fast_dixon_matrix_add_interleaved_rows(f_block,
                                               submatrix,
                                               submatrix_support,
                                               omit_idx,
                                               npolys,
                                               (omit_idx & 1) != 0);
        return;
    }

    for (slong deg = 0; deg <= current_degree; deg++) {
        slong remaining_slots = tuple_len - depth - 1;
        slong new_sum = current_sum + deg;
        slong min_possible = new_sum;
        slong max_possible = new_sum + remaining_slots * current_degree;

        if (min_possible > target_sum) {
            if (fast_dixon_profile_heavy_enabled()) {
                g_fast_dixon_profile.f_pruned_branches++;
                if (level_profile != NULL) {
                    level_profile->f_pruned_branches++;
                }
            }
            break;
        }
        if (max_possible < target_sum) {
            if (fast_dixon_profile_heavy_enabled()) {
                g_fast_dixon_profile.f_pruned_branches++;
                if (level_profile != NULL) {
                    level_profile->f_pruned_branches++;
                }
            }
            continue;
        }

        selected[depth] = &coeffs[poly_indices[depth]][deg];
        if (selected_hashes != NULL && coeff_hashes != NULL) {
            selected_hashes[depth] = coeff_hashes[poly_indices[depth]][deg];
        }
        fast_dixon_accumulate_f_recursive(f_block,
                                          coeffs,
                                          coeff_hashes,
                                          poly_indices,
                                          degrees,
                                          total_nvars,
                                          pos,
                                          npars,
                                          ctx,
                                          depth + 1,
                                          tuple_len,
                                          target_sum,
                                          new_sum,
                                          has_zero || (selected[depth]->nterms == 0),
                                          selected,
                                          selected_hashes,
                                          cache,
                                          omit_idx,
                                          npolys);
    }
}

static void fast_dixon_build_f_block(fast_dixon_matrix_t *out,
                                     fq_mvpoly_t **coeffs,
                                     ulong **coeff_hashes,
                                     const slong *degrees,
                                     slong total_nvars,
                                     slong pos,
                                     slong npars,
                                     const fq_nmod_ctx_t ctx,
                                     slong block_index,
                                     fast_dixon_subproblem_cache_t *cache)
{
    slong current_vars = total_nvars - pos;
    slong npolys = current_vars + 1;
    slong tuple_len = current_vars;
    slong lower_size = fast_dixon_level_size(degrees, total_nvars, pos + 1);
    slong *poly_indices = (slong *) flint_malloc((size_t) tuple_len * sizeof(slong));
    const fq_mvpoly_t **selected = (const fq_mvpoly_t **) flint_malloc((size_t) tuple_len * sizeof(fq_mvpoly_t *));
    ulong *selected_hashes = coeff_hashes != NULL
        ? (ulong *) flint_malloc((size_t) tuple_len * sizeof(ulong))
        : NULL;
    fast_dixon_subproblem_cache_t local_cache;
    int use_local_cache = 0;

    if (cache == NULL && coeff_hashes != NULL) {
        fast_dixon_subproblem_cache_init(&local_cache, tuple_len);
        cache = &local_cache;
        use_local_cache = 1;
    }

    fast_dixon_matrix_init(out, npolys * lower_size, lower_size, npars, ctx);

    for (slong omit_idx = 0; omit_idx < npolys; omit_idx++) {
        slong cursor = 0;

        for (slong p = 0; p < npolys; p++) {
            if (p == omit_idx) {
                continue;
            }
            poly_indices[cursor++] = p;
        }

        fast_dixon_accumulate_f_recursive(out,
                                          coeffs,
                                          coeff_hashes,
                                          poly_indices,
                                          degrees,
                                          total_nvars,
                                          pos,
                                          npars,
                                          ctx,
                                          0,
                                          tuple_len,
                                          block_index + 1,
                                          0,
                                          0,
                                          selected,
                                          selected_hashes,
                                          cache,
                                          omit_idx,
                                          npolys);
    }

    flint_free(poly_indices);
    flint_free(selected);
    if (selected_hashes != NULL) {
        flint_free(selected_hashes);
    }
    if (use_local_cache) {
        fast_dixon_subproblem_cache_clear(&local_cache);
    }
}

static void fast_dixon_build_s_block(fast_dixon_matrix_t *out,
                                     fq_mvpoly_t **coeffs,
                                     const slong *degrees,
                                     slong total_nvars,
                                     slong pos,
                                     slong npars,
                                     const fq_nmod_ctx_t ctx,
                                     slong coeff_degree,
                                     const slong *mult_table,
                                     const slong *row_counts,
                                     slong lower_size,
                                     slong rem_vars)
{
    slong current_vars = total_nvars - pos;
    slong npolys = current_vars + 1;

    fast_dixon_matrix_init(out,
                           current_vars * lower_size,
                           npolys * lower_size,
                           npars,
                           ctx);

    for (slong poly_idx = 0; poly_idx < npolys; poly_idx++) {
        fast_dixon_matrix_t poly_block;

        fast_dixon_build_p_block(&poly_block,
                                 &coeffs[poly_idx][coeff_degree],
                                 degrees,
                                 total_nvars,
                                 pos,
                                 mult_table,
                                 row_counts,
                                 lower_size,
                                 rem_vars);
        fast_dixon_matrix_scatter_columns(out, &poly_block, poly_idx, npolys);
        fast_dixon_matrix_clear(&poly_block);
    }
}

static void fast_dixon_build_matrix(fast_dixon_matrix_t *out,
                                    const fq_mvpoly_t *const *polys,
                                    const slong *degrees,
                                    slong total_nvars,
                                    slong pos,
                                    slong npars,
                                    const fq_nmod_ctx_t ctx)
{
    slong current_vars = total_nvars - pos;
    slong current_degree = degrees[pos];
    slong npolys = current_vars + 1;
    fq_mvpoly_t **coeffs;
    ulong **coeff_hashes = NULL;
    int collect_level_profile = fast_dixon_profile_enabled() &&
                                (pos == 0 || fast_dixon_profile_heavy_enabled());
    fast_dixon_level_profile_t *level_profile =
        collect_level_profile ? fast_dixon_get_level_profile(pos) : NULL;
    double call_start = get_wall_time();
    double phase_start;

    if (fast_dixon_profile_heavy_enabled()) {
        g_fast_dixon_profile.recursive_calls++;
        if (level_profile != NULL) {
            level_profile->calls++;
            level_profile->current_vars = current_vars;
            level_profile->current_degree = current_degree;
            level_profile->lower_size = fast_dixon_level_size(degrees, total_nvars, pos + 1);
        }
    }

    fast_dixon_trace_log(pos,
                         "[fast-dixon] enter level=%ld vars_left=%ld degree=%ld\n",
                         pos, current_vars, current_degree);

    phase_start = get_wall_time();
    coeffs = fast_dixon_split_polys_by_degree(polys, npolys, current_degree);
    coeff_hashes = fast_dixon_build_split_hashes(coeffs, npolys, current_degree, ctx);
    if (level_profile != NULL) {
        level_profile->split_time += get_wall_time() - phase_start;
    }

    if (current_vars == 1) {
        phase_start = get_wall_time();
        fast_dixon_build_base_case(out, coeffs, current_degree, npars, ctx);
        if (fast_dixon_profile_heavy_enabled()) {
            double base_elapsed = get_wall_time() - phase_start;
            g_fast_dixon_profile.base_case_calls++;
            if (level_profile != NULL) {
                level_profile->base_case_calls++;
                level_profile->base_case_time += base_elapsed;
                level_profile->output_rows = out->rows;
                level_profile->output_cols = out->cols;
                level_profile->total_time += get_wall_time() - call_start;
            }
        }
        fast_dixon_trace_log(pos,
                             "[fast-dixon] base level=%ld output=%ldx%ld time=%.3fs\n",
                             pos, out->rows, out->cols, get_wall_time() - call_start);
        fast_dixon_clear_split_hashes(coeff_hashes, npolys);
        fast_dixon_clear_split_coeffs(coeffs, npolys, current_degree);
        return;
    }

    {
        slong lower_size = fast_dixon_level_size(degrees, total_nvars, pos + 1);
        slong d_block_rows = current_vars * lower_size;
        slong num_f_blocks = current_vars * current_degree;
        slong rem_vars = total_nvars - pos - 1;
        slong mult_count = 0;
        slong mult_nvars = 0;
        slong *mult_table = fast_dixon_build_multiplier_table(degrees, total_nvars, pos,
                                                              &mult_count, &mult_nvars);
        slong *row_counts = rem_vars > 0
            ? (slong *) flint_malloc((size_t) rem_vars * sizeof(slong))
            : NULL;
        int parallel_threads = fast_dixon_parallel_threads_for_level(pos, current_degree, num_f_blocks);
        fast_dixon_matrix_t *s_blocks =
            (fast_dixon_matrix_t *) flint_malloc((size_t) current_degree * sizeof(fast_dixon_matrix_t));
        fast_dixon_matrix_support_t *s_supports =
            (fast_dixon_matrix_support_t *) flint_malloc((size_t) current_degree * sizeof(fast_dixon_matrix_support_t));
        fast_dixon_matrix_t *f_blocks =
            (fast_dixon_matrix_t *) flint_malloc((size_t) num_f_blocks * sizeof(fast_dixon_matrix_t));
        fast_dixon_matrix_support_t *f_supports =
            (fast_dixon_matrix_support_t *) flint_malloc((size_t) num_f_blocks * sizeof(fast_dixon_matrix_support_t));
        fast_dixon_matrix_t *d_blocks =
            (fast_dixon_matrix_t *) flint_malloc((size_t) current_degree * num_f_blocks * sizeof(fast_dixon_matrix_t));

        for (slong i = 0; i < current_degree; i++) {
            s_supports[i].rows = 0;
            s_supports[i].total_nonzero = 0;
            s_supports[i].row_counts = NULL;
            s_supports[i].row_offsets = NULL;
            s_supports[i].col_indices = NULL;
        }
        for (slong i = 0; i < num_f_blocks; i++) {
            f_supports[i].rows = 0;
            f_supports[i].total_nonzero = 0;
            f_supports[i].row_counts = NULL;
            f_supports[i].row_offsets = NULL;
            f_supports[i].col_indices = NULL;
        }

        for (slong i = 0; i < rem_vars; i++) {
            row_counts[i] = (i + 2) * degrees[pos + 1 + i];
        }

        if (level_profile != NULL) {
            level_profile->s_blocks += current_degree;
            level_profile->f_blocks += num_f_blocks;
            level_profile->d_block_products += current_degree * num_f_blocks;
            level_profile->output_rows = current_degree * d_block_rows;
            level_profile->output_cols = num_f_blocks * lower_size;
        }

        phase_start = get_wall_time();
#ifdef _OPENMP
        #pragma omp parallel for if(parallel_threads > 1) num_threads(parallel_threads) schedule(dynamic)
#endif
        for (slong i = 0; i < current_degree; i++) {
            fast_dixon_build_s_block(&s_blocks[i],
                                     coeffs,
                                     degrees,
                                     total_nvars,
                                     pos,
                                     npars,
                                     ctx,
                                     i,
                                     mult_table,
                                     row_counts,
                                     lower_size,
                                     rem_vars);
            fast_dixon_matrix_support_init(&s_supports[i], &s_blocks[i]);
        }
        if (level_profile != NULL) {
            level_profile->s_build_time += get_wall_time() - phase_start;
        }

        phase_start = get_wall_time();
#ifdef _OPENMP
        #pragma omp parallel for if(parallel_threads > 1) num_threads(parallel_threads) schedule(dynamic)
#endif
        for (slong j = 0; j < num_f_blocks; j++) {
            fast_dixon_subproblem_cache_t *child_cache =
                fast_dixon_get_shared_subproblem_cache(pos + 1, total_nvars, ctx);
            fast_dixon_build_f_block(&f_blocks[j],
                                     coeffs,
                                     coeff_hashes,
                                     degrees,
                                     total_nvars,
                                     pos,
                                     npars,
                                     ctx,
                                     j,
                                     child_cache);
            fast_dixon_matrix_support_init(&f_supports[j], &f_blocks[j]);
        }
        if (level_profile != NULL) {
            level_profile->f_build_time += get_wall_time() - phase_start;
        }

        phase_start = get_wall_time();
        {
            slong total_products = current_degree * num_f_blocks;
#ifdef _OPENMP
            #pragma omp parallel for if(parallel_threads > 1) num_threads(parallel_threads) schedule(dynamic)
#endif
            for (slong block_idx = 0; block_idx < total_products; block_idx++) {
                slong i = block_idx / num_f_blocks;
                slong j = block_idx % num_f_blocks;

                fast_dixon_matrix_init(&d_blocks[block_idx],
                                       d_block_rows, lower_size, npars, ctx);
                fast_dixon_matrix_mul_accumulate_supported(&d_blocks[block_idx],
                                                           &s_blocks[i],
                                                           &s_supports[i],
                                                           &f_blocks[j],
                                                           &f_supports[j],
                                                           0);
            }
        }
        if (level_profile != NULL) {
            level_profile->mul_time += get_wall_time() - phase_start;
        }

        phase_start = get_wall_time();
        for (slong i = 1; i < current_degree; i++) {
            for (slong j = 0; j < num_f_blocks - 1; j++) {
                fast_dixon_matrix_add_inplace(&d_blocks[i * num_f_blocks + j],
                                              &d_blocks[(i - 1) * num_f_blocks + (j + 1)],
                                              0);
            }
        }
        if (level_profile != NULL) {
            level_profile->recurrence_time += get_wall_time() - phase_start;
        }

        phase_start = get_wall_time();
        fast_dixon_matrix_init(out,
                               current_degree * d_block_rows,
                               num_f_blocks * lower_size,
                               npars,
                               ctx);

        {
            slong total_products = current_degree * num_f_blocks;
#ifdef _OPENMP
            #pragma omp parallel for if(parallel_threads > 1) num_threads(parallel_threads) schedule(static)
#endif
            for (slong block_idx = 0; block_idx < total_products; block_idx++) {
                slong i = block_idx / num_f_blocks;
                slong j = block_idx % num_f_blocks;

                fast_dixon_matrix_copy_block(out,
                                             i * d_block_rows,
                                             j * lower_size,
                                             &d_blocks[block_idx]);
            }
        }
        if (level_profile != NULL) {
            level_profile->assemble_time += get_wall_time() - phase_start;
        }

        for (slong i = 0; i < current_degree * num_f_blocks; i++) {
            fast_dixon_matrix_clear(&d_blocks[i]);
        }
        for (slong i = 0; i < current_degree; i++) {
            fast_dixon_matrix_support_clear(&s_supports[i]);
            fast_dixon_matrix_clear(&s_blocks[i]);
        }
        for (slong i = 0; i < num_f_blocks; i++) {
            fast_dixon_matrix_support_clear(&f_supports[i]);
            fast_dixon_matrix_clear(&f_blocks[i]);
        }

        flint_free(d_blocks);
        flint_free(s_blocks);
        flint_free(s_supports);
        flint_free(f_blocks);
        flint_free(f_supports);
        if (row_counts != NULL) {
            flint_free(row_counts);
        }
        flint_free(mult_table);
        (void) mult_count;
        (void) mult_nvars;
    }

    if (level_profile != NULL) {
        level_profile->total_time += get_wall_time() - call_start;
    }

    fast_dixon_trace_log(pos,
                         "[fast-dixon] exit  level=%ld output=%ldx%ld time=%.3fs\n",
                         pos, out->rows, out->cols, get_wall_time() - call_start);

    fast_dixon_clear_split_hashes(coeff_hashes, npolys);
    fast_dixon_clear_split_coeffs(coeffs, npolys, current_degree);
}

static det_method_t choose_fast_dixon_det_method(slong matrix_size, slong npars)
{
    det_method_t coeff_method;

#ifdef _OPENMP
    if (npars > 1) {
        coeff_method = DET_METHOD_INTERPOLATION;
    } else
#endif
    if (matrix_size < 9) {
        coeff_method = DET_METHOD_RECURSIVE;
    } else {
        coeff_method = DET_METHOD_KRONECKER;
    }

    if (dixon_global_method_step4 != -1) {
        coeff_method = dixon_global_method_step4;
    }
    fast_dixon_info_log("  Determinant method: %s\n", fast_dixon_det_method_name(coeff_method));

    return coeff_method;
}

static void clear_fast_dixon_coeff_matrix(fq_mvpoly_t **coeff_matrix, slong matrix_size)
{
    if (coeff_matrix == NULL) {
        return;
    }

    for (slong i = 0; i < matrix_size; i++) {
        for (slong j = 0; j < matrix_size; j++) {
            fq_mvpoly_clear(&coeff_matrix[i][j]);
        }
        flint_free(coeff_matrix[i]);
    }
    flint_free(coeff_matrix);
}

static fq_mvpoly_t ***fast_dixon_build_pointer_grid(const fast_dixon_matrix_t *matrix)
{
    fq_mvpoly_t ***grid = (fq_mvpoly_t ***) flint_malloc((size_t) matrix->rows * sizeof(fq_mvpoly_t **));

    for (slong i = 0; i < matrix->rows; i++) {
        grid[i] = (fq_mvpoly_t **) flint_malloc((size_t) matrix->cols * sizeof(fq_mvpoly_t *));
        for (slong j = 0; j < matrix->cols; j++) {
            fq_mvpoly_t *entry = &FAST_DIXON_ENTRY(matrix, i, j);
            grid[i][j] = entry->nterms > 0 ? entry : NULL;
        }
    }

    return grid;
}

static void fast_dixon_free_pointer_grid(fq_mvpoly_t ***grid, slong rows)
{
    if (grid == NULL) {
        return;
    }

    for (slong i = 0; i < rows; i++) {
        flint_free(grid[i]);
    }
    flint_free(grid);
}

static slong fast_dixon_matrix_count_nonzero(const fast_dixon_matrix_t *matrix)
{
    slong nnz = 0;

    for (slong i = 0; i < matrix->rows * matrix->cols; i++) {
        if (matrix->entries[i].nterms > 0) {
            nnz++;
        }
    }

    return nnz;
}

static void fast_dixon_print_bottleneck_hint(const fast_dixon_level_profile_t *top)
{
    double phases[4];
    const char *names[4];
    slong best = 0;

    if (!fast_dixon_profile_enabled() || top == NULL) {
        return;
    }

    phases[0] = top->f_build_time;
    phases[1] = top->mul_time;
    phases[2] = top->s_build_time;
    phases[3] = top->assemble_time + top->recurrence_time;

    names[0] = "F-block recursive construction";
    names[1] = "S_i * F_j block multiplication";
    names[2] = "S-block construction";
    names[3] = "block recurrence / assembly";

    for (slong i = 1; i < 4; i++) {
        if (phases[i] > phases[best]) {
            best = i;
        }
    }

    fast_dixon_info_log("  Bottleneck candidate: %s (top-level %.3f s)\n",
                        names[best], phases[best]);
    if (best == 0) {
        fast_dixon_info_log("    -> likely dominated by repeatedly constructing lower-level Dixon blocks for many valid F_j tuples.\n");
    } else if (best == 1) {
        fast_dixon_info_log("    -> likely dominated by repeated fq_mvpoly_mul/fq_mvpoly_add during all S_i * F_j products.\n");
    } else if (best == 2) {
        fast_dixon_info_log("    -> likely dominated by scattering coefficients into Sylvester-style S_i blocks.\n");
    } else {
        fast_dixon_info_log("    -> likely dominated by copying/accumulating block matrices after products are computed.\n");
    }
}

static double fast_dixon_report_time(double value)
{
    return value < 0.0 ? 0.0 : value;
}

static void fast_dixon_print_profile_report(const fast_dixon_matrix_t *full_matrix)
{
    fast_dixon_level_profile_t *top = fast_dixon_get_level_profile(0);
    slong nnz;

    if (!fast_dixon_profile_enabled() || top == NULL) {
        return;
    }

    fast_dixon_info_log("\n  Fast Dixon construction profile:\n");
    if (fast_dixon_profile_heavy_enabled()) {
        nnz = fast_dixon_matrix_count_nonzero(full_matrix);
        fast_dixon_info_log("    Recursive calls: %ld total (%ld base cases)\n",
                            g_fast_dixon_profile.recursive_calls,
                            g_fast_dixon_profile.base_case_calls);
        fast_dixon_info_log("    Matrix allocations: %ld blocks, %llu total entry slots\n",
                            g_fast_dixon_profile.matrix_inits,
                            g_fast_dixon_profile.matrix_entries_allocated);
        fast_dixon_info_log("    Full matrix support: %ld / %ld non-zero entries (%.2f%% density)\n",
                            nnz,
                            full_matrix->rows * full_matrix->cols,
                            full_matrix->rows * full_matrix->cols > 0
                                ? (100.0 * (double) nnz /
                                   (double) (full_matrix->rows * full_matrix->cols))
                                : 0.0);
        fast_dixon_info_log("    Polynomial kernel ops: %llu muls, %llu generated terms, %llu adds/subs\n",
                            g_fast_dixon_profile.poly_mul_ops,
                            g_fast_dixon_profile.poly_mul_generated_terms,
                            g_fast_dixon_profile.poly_add_ops);
        fast_dixon_info_log("    F-block tuple search: %llu recursion nodes, %llu leaves, %llu valid tuples, %llu zero skips, %llu pruned branches\n",
                            g_fast_dixon_profile.f_recursive_nodes,
                            g_fast_dixon_profile.f_leaf_calls,
                            g_fast_dixon_profile.f_valid_leaves,
                            g_fast_dixon_profile.f_zero_skips,
                            g_fast_dixon_profile.f_pruned_branches);
        fast_dixon_info_log("    Recursive subproblem cache: %llu hits, %llu misses (hit rate %.2f%%)\n",
                            g_fast_dixon_profile.cache_hits,
                            g_fast_dixon_profile.cache_misses,
                            (g_fast_dixon_profile.cache_hits + g_fast_dixon_profile.cache_misses) > 0
                                ? (100.0 * (double) g_fast_dixon_profile.cache_hits /
                                   (double) (g_fast_dixon_profile.cache_hits + g_fast_dixon_profile.cache_misses))
                                : 0.0);
    } else {
        fast_dixon_info_log("    Detailed recursive counters are disabled at -v 2 to avoid profiling overhead.\n");
    }

    fast_dixon_info_log("    Top-level phase times:\n");
    fast_dixon_info_log("      split coefficients      : %.3f s\n", top->split_time);
    fast_dixon_info_log("      build S_i blocks        : %.3f s\n", top->s_build_time);
    fast_dixon_info_log("      build F_j blocks        : %.3f s (child lower-level builds %.3f s)\n",
                        top->f_build_time, top->f_child_build_time);
    fast_dixon_info_log("      multiply S_i * F_j      : %.3f s\n", top->mul_time);
    fast_dixon_info_log("      recurrence accumulation : %.3f s\n", top->recurrence_time);
    fast_dixon_info_log("      final block assembly    : %.3f s\n", top->assemble_time);

    if (fast_dixon_profile_heavy_enabled()) {
        fast_dixon_info_log("    Per-level inclusive summary:\n");
        fast_dixon_info_log("      (lower-level phase totals are aggregated over many recursive calls and may overlap with child-build time; do not sum them directly)\n");
        for (slong pos = 0; pos < FAST_DIXON_MAX_LEVELS; pos++) {
            fast_dixon_level_profile_t *level = &g_fast_dixon_profile.levels[pos];

            if (level->calls == 0) {
                continue;
            }

            fast_dixon_info_log("      level %ld: calls=%ld vars_left=%ld degree=%ld lower=%ld output=%ldx%ld total=%.3f s\n",
                                pos,
                                level->calls,
                                level->current_vars,
                                level->current_degree,
                                level->lower_size,
                                level->output_rows,
                                level->output_cols,
                                level->total_time);
            fast_dixon_info_log("               S=%ld F=%ld products=%ld | split=%.3f Sbuild=%.3f Fbuild=%.3f child=%.3f mul=%.3f recur=%.3f asm=%.3f base=%.3f\n",
                                level->s_blocks,
                                level->f_blocks,
                                level->d_block_products,
                                fast_dixon_report_time(level->split_time),
                                fast_dixon_report_time(level->s_build_time),
                                fast_dixon_report_time(level->f_build_time),
                                fast_dixon_report_time(level->f_child_build_time),
                                fast_dixon_report_time(level->mul_time),
                                fast_dixon_report_time(level->recurrence_time),
                                fast_dixon_report_time(level->assemble_time),
                                fast_dixon_report_time(level->base_case_time));
            fast_dixon_info_log("               F-search nodes=%ld leaves=%ld valid=%ld zero=%ld pruned=%ld cache=%ld/%ld\n",
                                level->f_recursive_nodes,
                                level->f_leaf_calls,
                                level->f_valid_leaves,
                                level->f_zero_skips,
                                level->f_pruned_branches,
                                level->cache_hits,
                                level->cache_misses);
        }
    }

    fast_dixon_print_bottleneck_hint(top);
}

static void fast_dixon_extract_square_submatrix(fq_mvpoly_t ***coeff_matrix_out,
                                                slong *matrix_size_out,
                                                const fast_dixon_matrix_t *full_matrix,
                                                slong npars)
{
    fq_mvpoly_t ***grid = fast_dixon_build_pointer_grid(full_matrix);
    fq_mvpoly_t ***trimmed_grid = NULL;
    slong *trimmed_rows = NULL;
    slong *trimmed_cols = NULL;
    slong *row_idx_array = NULL;
    slong *col_idx_array = NULL;
    slong num_rows = 0;
    slong num_cols = 0;
    slong submat_rank = 0;
    slong trimmed_nrows = 0;
    slong trimmed_ncols = 0;
    double support_scan_start;
    double support_scan_elapsed;
    double trim_grid_start;
    double trim_grid_elapsed;
    double step3_wall_start;
    double rank_select_start;
    double rank_select_elapsed;
    double copy_start;
    double copy_elapsed;

    if (fast_dixon_profile_enabled()) {
        fast_dixon_info_log("  Raw Dixon matrix size: %ld x %ld\n",
                            full_matrix->rows, full_matrix->cols);
    }

    support_scan_start = get_wall_time();
    trimmed_rows = (slong *) flint_malloc((size_t) full_matrix->rows * sizeof(slong));
    trimmed_cols = (slong *) flint_malloc((size_t) full_matrix->cols * sizeof(slong));

    for (slong i = 0; i < full_matrix->rows; i++) {
        int nonzero = 0;
        for (slong j = 0; j < full_matrix->cols; j++) {
            if (grid[i][j] != NULL && grid[i][j]->nterms > 0) {
                nonzero = 1;
                break;
            }
        }
        if (nonzero) {
            trimmed_rows[trimmed_nrows++] = i;
        }
    }

    for (slong j = 0; j < full_matrix->cols; j++) {
        int nonzero = 0;
        for (slong i = 0; i < full_matrix->rows; i++) {
            if (grid[i][j] != NULL && grid[i][j]->nterms > 0) {
                nonzero = 1;
                break;
            }
        }
        if (nonzero) {
            trimmed_cols[trimmed_ncols++] = j;
        }
    }
    support_scan_elapsed = get_wall_time() - support_scan_start;

    if (trimmed_nrows == 0 || trimmed_ncols == 0) {
        *coeff_matrix_out = NULL;
        *matrix_size_out = 0;
        fast_dixon_free_pointer_grid(grid, full_matrix->rows);
        flint_free(trimmed_rows);
        flint_free(trimmed_cols);
        fast_dixon_info_log("Warning: Fast Dixon matrix has no non-zero support\n");
        return;
    }

    trim_grid_start = get_wall_time();
    trimmed_grid = (fq_mvpoly_t ***) flint_malloc((size_t) trimmed_nrows * sizeof(fq_mvpoly_t **));
    for (slong i = 0; i < trimmed_nrows; i++) {
        trimmed_grid[i] = (fq_mvpoly_t **) flint_malloc((size_t) trimmed_ncols * sizeof(fq_mvpoly_t *));
        for (slong j = 0; j < trimmed_ncols; j++) {
            trimmed_grid[i][j] = grid[trimmed_rows[i]][trimmed_cols[j]];
        }
    }
    trim_grid_elapsed = get_wall_time() - trim_grid_start;

    fast_dixon_info_log("  Dixon matrix size: %ld x %ld\n", trimmed_nrows, trimmed_ncols);
    if (fast_dixon_profile_heavy_enabled()) {
        fast_dixon_info_log("  Support scan: %.3f s | trim grid build: %.3f s\n",
                            support_scan_elapsed, trim_grid_elapsed);
    }

    fast_dixon_info_log("\nStep 3: Extract maximal-rank submatrix\n");
    step3_wall_start = get_wall_time();
    rank_select_start = get_wall_time();
    if (npars == 0) {
        fq_nmod_mat_t eval_mat;
        fq_nmod_mat_init(eval_mat, trimmed_nrows, trimmed_ncols, full_matrix->ctx);

        for (slong i = 0; i < trimmed_nrows; i++) {
            for (slong j = 0; j < trimmed_ncols; j++) {
                fq_mvpoly_t *entry = trimmed_grid[i][j];
                if (entry != NULL && entry->nterms > 0) {
                    fq_nmod_set(fq_nmod_mat_entry(eval_mat, i, j),
                                entry->terms[0].coeff,
                                full_matrix->ctx);
                } else {
                    fq_nmod_zero(fq_nmod_mat_entry(eval_mat, i, j), full_matrix->ctx);
                }
            }
        }

        submat_rank = fq_nmod_mat_rank(eval_mat, full_matrix->ctx);
        submat_rank = FLINT_MIN(submat_rank, FLINT_MIN(trimmed_nrows, trimmed_ncols));

        if (submat_rank > 0) {
            row_idx_array = (slong *) flint_malloc((size_t) submat_rank * sizeof(slong));
            col_idx_array = (slong *) flint_malloc((size_t) submat_rank * sizeof(slong));
            for (slong i = 0; i < submat_rank; i++) {
                row_idx_array[i] = i;
                col_idx_array[i] = i;
            }
        }

        fq_nmod_mat_clear(eval_mat, full_matrix->ctx);
        num_rows = submat_rank;
        num_cols = submat_rank;
    } else {
        find_fq_optimal_maximal_rank_submatrix(trimmed_grid,
                                               trimmed_nrows,
                                               trimmed_ncols,
                                               &row_idx_array,
                                               &col_idx_array,
                                               &num_rows,
                                               &num_cols,
                                               npars,
                                               g_dixon_fast_use_ksy_precondition
                                                   ? g_dixon_fast_ksy_constant_col
                                                   : -1);
        submat_rank = FLINT_MIN(num_rows, num_cols);
    }
    rank_select_elapsed = get_wall_time() - rank_select_start;

    if (submat_rank <= 0) {
        *coeff_matrix_out = NULL;
        *matrix_size_out = 0;
        fast_dixon_free_pointer_grid(trimmed_grid, trimmed_nrows);
        fast_dixon_free_pointer_grid(grid, full_matrix->rows);
        flint_free(trimmed_rows);
        flint_free(trimmed_cols);
        if (row_idx_array != NULL) flint_free(row_idx_array);
        if (col_idx_array != NULL) flint_free(col_idx_array);
        fast_dixon_info_log("Warning: Fast Dixon matrix has rank 0\n");
        dixon_maybe_print_step_time("Step 3", get_wall_time() - step3_wall_start);
        return;
    }

    copy_start = get_wall_time();
    *coeff_matrix_out = (fq_mvpoly_t **) flint_malloc((size_t) submat_rank * sizeof(fq_mvpoly_t *));
    for (slong i = 0; i < submat_rank; i++) {
        (*coeff_matrix_out)[i] = (fq_mvpoly_t *) flint_malloc((size_t) submat_rank * sizeof(fq_mvpoly_t));
        for (slong j = 0; j < submat_rank; j++) {
            fq_mvpoly_t *source = trimmed_grid[row_idx_array[i]][col_idx_array[j]];
            if (source != NULL) {
                fq_mvpoly_copy(&(*coeff_matrix_out)[i][j], source);
            } else {
                fast_dixon_poly_init_zero(&(*coeff_matrix_out)[i][j], 0, npars, full_matrix->ctx);
            }
        }
    }
    copy_elapsed = get_wall_time() - copy_start;

    *matrix_size_out = submat_rank;
    fast_dixon_info_log("  Submatrix size: %ld x %ld\n", submat_rank, submat_rank);
    if (fast_dixon_profile_heavy_enabled()) {
        fast_dixon_info_log("  Rank/select: %.3f s | copy selected submatrix: %.3f s\n",
                            rank_select_elapsed, copy_elapsed);
    }
    dixon_maybe_print_step_time("Step 3", get_wall_time() - step3_wall_start);

    fast_dixon_free_pointer_grid(trimmed_grid, trimmed_nrows);
    fast_dixon_free_pointer_grid(grid, full_matrix->rows);
    flint_free(trimmed_rows);
    flint_free(trimmed_cols);
    if (row_idx_array != NULL) flint_free(row_idx_array);
    if (col_idx_array != NULL) flint_free(col_idx_array);
}

static void fast_dixon_print_small_matrix(const fast_dixon_matrix_t *matrix,
                                          const char *name)
{
    fq_mvpoly_t **rows;

    if (g_dixon_verbose_level < 3 || matrix->rows > 10 || matrix->cols > 10) {
        return;
    }

    rows = (fq_mvpoly_t **) flint_malloc((size_t) matrix->rows * sizeof(fq_mvpoly_t *));
    for (slong i = 0; i < matrix->rows; i++) {
        rows[i] = matrix->entries + i * matrix->cols;
    }

    print_fq_matrix_mvpoly(rows, matrix->rows, matrix->cols, name, 1);
    flint_free(rows);
}

static void fq_dixon_fast_resultant_common(fq_mvpoly_t *result, fq_mvpoly_t *polys,
                                           slong nvars, slong npars,
                                           char **var_names, char **par_names,
                                           const char *gen_name)
{
    const fq_mvpoly_t **poly_ptrs;
    slong *degrees;
    fast_dixon_matrix_t full_matrix;
    fq_mvpoly_t **coeff_matrix = NULL;
    slong matrix_size = 0;

    full_matrix.rows = 0;
    full_matrix.cols = 0;
    full_matrix.npars = 0;
    full_matrix.ctx = NULL;
    full_matrix.entries = NULL;

    cleanup_unified_workspace();
    fast_dixon_profile_reset(nvars);
    fast_dixon_subproblem_cache_reset_all();

    poly_ptrs = (const fq_mvpoly_t **) flint_malloc((size_t) (nvars + 1) * sizeof(fq_mvpoly_t *));
    for (slong i = 0; i < nvars + 1; i++) {
        poly_ptrs[i] = &polys[i];
    }

    degrees = (slong *) flint_malloc((size_t) nvars * sizeof(slong));
    fast_dixon_compute_degree_bounds(degrees, poly_ptrs, nvars + 1, nvars);

#ifdef _OPENMP
    {
        int cache_threads = omp_get_max_threads();
        if (cache_threads > 16) {
            cache_threads = 16;
        }
        fast_dixon_subproblem_cache_prepare_pool(cache_threads);
    }
#else
    fast_dixon_subproblem_cache_prepare_pool(1);
#endif

    for (slong i = 0; i < nvars; i++) {
        if (degrees[i] <= 0) {
            fast_dixon_info_log("\nFast Dixon fallback: variable %ld has degree %ld, using standard Dixon pipeline.\n",
                                i, degrees[i]);
            fast_dixon_subproblem_cache_clear_all();
            flint_free(degrees);
            flint_free(poly_ptrs);
            fq_dixon_resultant_with_names(result, polys, nvars, npars,
                                          var_names, par_names, gen_name);
            return;
        }
    }

    fast_dixon_info_log("\nStep 1: Skip explicit Dixon polynomial construction\n");
    fast_dixon_info_log("  Recursive block construction builds the Dixon matrix directly.\n");
    fast_dixon_info_log("\nStep 2: Construct Dixon matrix (recursive block construction)\n");
    {
        clock_t step2_cpu_start = clock();
        double step2_wall_start = get_wall_time();

        fast_dixon_build_matrix(&full_matrix, poly_ptrs, degrees, nvars, 0, npars, polys[0].ctx);

        fast_dixon_print_small_matrix(&full_matrix, "Dixon");
        fast_dixon_print_profile_report(&full_matrix);
        fast_dixon_subproblem_cache_clear_all();

        dixon_maybe_print_parallel_step_time("Step 2",
                                             ((double) (clock() - step2_cpu_start) / CLOCKS_PER_SEC),
                                             get_wall_time() - step2_wall_start);
    }

    fast_dixon_extract_square_submatrix(&coeff_matrix, &matrix_size, &full_matrix, npars);

    if (matrix_size > 0) {
        det_method_t coeff_method;
        slong res_deg_bound;
        clock_t step4_cpu_start;
        double step4_wall_start;

        fast_dixon_info_log("\nStep 4: Compute resultant\n");
        res_deg_bound = compute_fq_dixon_resultant_degree_bound(polys, nvars + 1, nvars, npars);
        coeff_method = choose_fast_dixon_det_method(matrix_size, npars);
        step4_cpu_start = clock();
        step4_wall_start = get_wall_time();

        compute_fq_coefficient_matrix_det(result, coeff_matrix, matrix_size,
                                          npars, polys[0].ctx, coeff_method, res_deg_bound);
        dixon_maybe_print_step_method_time("Step 4",
                                           coeff_method,
                                           ((double) (clock() - step4_cpu_start) / CLOCKS_PER_SEC),
                                           get_wall_time() - step4_wall_start);
        fq_mvpoly_make_monic(result);
        print_resultant_summary(result, par_names, npars);

        if (result->nterms <= 100) {
            if (var_names || par_names || gen_name) {
                fq_mvpoly_print_with_names(result, "  Final Resultant",
                                           NULL, par_names, gen_name, 0);
            } else {
                fq_mvpoly_print(result, "  Final Resultant");
            }
        } else {
            fast_dixon_info_log("  Final resultant too large to display (%ld terms)\n",
                                result->nterms);
        }
    } else {
        fq_mvpoly_init(result, 0, npars, polys[0].ctx);
        fast_dixon_info_log("Warning: Empty fast Dixon coefficient matrix, resultant is 0\n");
    }

    clear_fast_dixon_coeff_matrix(coeff_matrix, matrix_size);
    fast_dixon_matrix_clear(&full_matrix);
    flint_free(degrees);
    flint_free(poly_ptrs);

    fast_dixon_info_log("\n=== Fast Dixon Resultant Computation Complete ===\n");
}

void fq_dixon_fast_resultant(fq_mvpoly_t *result, fq_mvpoly_t *polys,
                             slong nvars, slong npars)
{
    fq_dixon_fast_resultant_common(result, polys, nvars, npars,
                                   NULL, NULL, NULL);
}

void fq_dixon_fast_resultant_with_names(fq_mvpoly_t *result, fq_mvpoly_t *polys,
                                        slong nvars, slong npars,
                                        char **var_names, char **par_names,
                                        const char *gen_name)
{
    fq_dixon_fast_resultant_common(result, polys, nvars, npars,
                                   var_names, par_names, gen_name);
}
