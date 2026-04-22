#include "fq_sparse_interpolation.h"
#include "nmod_vec_extra.h"
#include <flint/nmod_vec.h>

#include <limits.h>
#include <stdarg.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define SPARSE_BASE_RECOVERY_ATTEMPTS 24
#define SPARSE_DEFAULT_START_BUDGET 32
#define SPARSE_ESTIMATE_CAP ((slong) (LONG_MAX / 4))
#define SPARSE_FAST_ROOT_THRESHOLD 64
#define SPARSE_FAST_ROOT_SPLIT_ATTEMPTS 12
#define SPARSE_BM_HIST_BINS 11

typedef int (*sparse_probe_fn_t)(mp_limb_t* a,
                                 mp_limb_t* h,
                                 const mp_limb_t* alpha,
                                 slong nvars,
                                 slong Tbound,
                                 void* user_ctx,
                                 nmod_t mod);

typedef struct {
    double total;
    double prepare;
    double eval_entries;
    double det;
    double inv;
    double direct_trace;
    double adjugate_trace;
    slong point_count;
    slong entry_term_count;
    slong inverse_count;
    slong singular_count;
} sparse_probe_timing_t;

typedef struct {
    const nmod_mpoly_struct* f;
    slong nvars;
    const nmod_mpoly_ctx_struct* mctx;
    sparse_probe_timing_t* timing;
} explicit_poly_probe_ctx_t;

typedef struct {
    const poly_mat_t* A;
    slong matrix_size;
    slong nvars;
    const nmod_mpoly_ctx_struct* mctx;
    sparse_probe_timing_t* timing;
} determinant_probe_ctx_t;

typedef struct {
    slong len;
    mp_limb_t* coeffs;
    mp_limb_t* roots;
    mp_limb_t* powers;
    ulong* exps;
} geometric_poly_cache_t;

typedef struct {
    slong npoints;
    slong leaf_count;
    slong node_count;
    mp_limb_t* points;
    nmod_poly_struct* prod;
    nmod_poly_struct* rem;
} sparse_multipoint_tree_t;

typedef struct {
    double total;
    double probe;
    double bm;
    double bm_discrepancy;
    double bm_update;
    double bm_finalize;
    double bm_update_prepare;
    double bm_update_axpy;
    double bm_update_trim;
    double bm_update_promote;
    slong bm_update_count;
    slong bm_promote_count;
    slong bm_update_len_sum;
    slong bm_update_hist[SPARSE_BM_HIST_BINS];
    double root_factor;
    double root_squarefree;
    double root_frobenius;
    double root_split;
    double root_sort;
    double root_fallback_factor;
    double vinvert;
    double vinvert_denom;
    double vinvert_tree_build;
    double vinvert_build;
    double vinvert_mul;
    double vinvert_extract;
    double vinvert_eval;
    double reconstruct;
} sparse_interpolation_timing_t;

typedef struct {
    int initialized;
    slong capacity_N;
    slong processed_terms;
    slong L;
    slong m;
    slong c_len;
    slong b_len;
    mp_limb_t b_inv;
    mp_limb_t* c;
    mp_limb_t* bpoly;
    mp_limb_t* temp;
    mp_limb_t* s_rev;
} sparse_bm_state_t;

typedef struct {
    int valid;
    slong nvars;
    mp_limb_t* alpha;
    sparse_bm_state_t bm_state;
} sparse_scout_reuse_t;

typedef struct {
    int used;
    int likely_underbudget;
    slong attempts;
    slong attempt_limit;
    slong budget;
    slong last_relation_degree;
    slong last_reconstructed_terms;
    size_t workspace_bytes;
    sparse_probe_timing_t probe_timing;
    sparse_interpolation_timing_t interpolation_timing;
    char last_reason[256];
} sparse_recovery_info_t;

flint_rand_t global_state;
static int global_state_initialized = 0;

static void ensure_global_state(void)
{
    if (!global_state_initialized) {
        flint_rand_init(global_state);
        global_state_initialized = 1;
    }
}

static void clear_global_state_if_needed(void)
{
    if (global_state_initialized) {
        flint_rand_clear(global_state);
        global_state_initialized = 0;
    }
}

static double timer_seconds_since(clock_t start);
static void sparse_multipoint_tree_init_empty(sparse_multipoint_tree_t* tree);
static int sparse_bm_inner_timing_enabled(void);
static int sparse_bm_hist_enabled(void);
static inline void sparse_bm_axpy(mp_limb_t* res, const mp_limb_t* vec, slong len, mp_limb_t negcoef, nmod_t mod);
static int sparse_bm_hist_bin(slong len);
static void sparse_print_bm_histogram(const sparse_interpolation_timing_t* timing, const char* prefix);
static void sparse_bm_state_init_empty(sparse_bm_state_t* state);
static void sparse_bm_state_clear(sparse_bm_state_t* state);
static int sparse_bm_state_ensure_capacity(sparse_bm_state_t* state, slong N);
static void sparse_bm_state_reset_sequence(sparse_bm_state_t* state);
static void BM_timed_incremental(nmod_poly_t C, const mp_limb_t* s, slong N, nmod_t mod,
                                 sparse_bm_state_t* state, sparse_interpolation_timing_t* timing);
static void sparse_scout_reuse_init_empty(sparse_scout_reuse_t* reuse);
static void sparse_scout_reuse_clear(sparse_scout_reuse_t* reuse);
static int sparse_scout_reuse_prepare(sparse_scout_reuse_t* reuse, slong nvars, slong budget);

void timer_start(my_timer_t* t)
{
    t->start = clock();
}

static inline void sparse_bm_axpy(mp_limb_t* res, const mp_limb_t* vec, slong len, mp_limb_t negcoef, nmod_t mod)
{
    if (len <= 0 || negcoef == 0) {
        return;
    }

    if (len <= 8) {
        for (slong i = 0; i < len; i++) {
            res[i] = nmod_add(res[i], nmod_mul(vec[i], negcoef, mod), mod);
        }
        return;
    }

    _nmod_vec_scalar_addmul_nmod(res, vec, len, negcoef, mod);
}

static int sparse_bm_hist_bin(slong len)
{
    if (len <= 8) return 0;
    if (len <= 16) return 1;
    if (len <= 32) return 2;
    if (len <= 64) return 3;
    if (len <= 128) return 4;
    if (len <= 256) return 5;
    if (len <= 512) return 6;
    if (len <= 1024) return 7;
    if (len <= 2048) return 8;
    if (len <= 4096) return 9;
    return 10;
}

static void sparse_print_bm_histogram(const sparse_interpolation_timing_t* timing, const char* prefix)
{
    static const char* labels[SPARSE_BM_HIST_BINS] = {
        "<=8", "9-16", "17-32", "33-64", "65-128",
        "129-256", "257-512", "513-1024", "1025-2048",
        "2049-4096", ">4096"
    };

    if (timing == NULL || timing->bm_update_count <= 0) {
        return;
    }

    printf("%s[Sparse BM b_len] updates=%ld promotes=%ld avg_len=%.2f\n",
           prefix,
           timing->bm_update_count,
           timing->bm_promote_count,
           timing->bm_update_count > 0
               ? ((double) timing->bm_update_len_sum) / ((double) timing->bm_update_count)
               : 0.0);

    printf("%s[Sparse BM b_len hist]", prefix);
    for (int i = 0; i < SPARSE_BM_HIST_BINS; i++) {
        if (timing->bm_update_hist[i] > 0) {
            printf(" %s:%ld", labels[i], timing->bm_update_hist[i]);
        }
    }
    printf("\n");
}

static void sparse_bm_state_init_empty(sparse_bm_state_t* state)
{
    if (state == NULL) return;
    memset(state, 0, sizeof(*state));
    state->m = 1;
    state->c_len = 1;
    state->b_len = 1;
    state->b_inv = 1;
}

static void sparse_bm_state_clear(sparse_bm_state_t* state)
{
    if (state == NULL) return;
    free(state->c);
    free(state->bpoly);
    free(state->temp);
    free(state->s_rev);
    sparse_bm_state_init_empty(state);
}

static int sparse_bm_state_ensure_capacity(sparse_bm_state_t* state, slong N)
{
    if (state == NULL || N <= 0) return 0;
    if (state->capacity_N >= N &&
        state->c != NULL && state->bpoly != NULL &&
        state->temp != NULL && state->s_rev != NULL) {
        return 1;
    }

    slong oldN = state->capacity_N;
    mp_limb_t* new_c = (mp_limb_t*) realloc(state->c, (size_t) (N + 1) * sizeof(mp_limb_t));
    mp_limb_t* new_bpoly = (mp_limb_t*) realloc(state->bpoly, (size_t) (N + 1) * sizeof(mp_limb_t));
    mp_limb_t* new_temp = (mp_limb_t*) realloc(state->temp, (size_t) (N + 1) * sizeof(mp_limb_t));
    mp_limb_t* new_s_rev = (mp_limb_t*) realloc(state->s_rev, (size_t) (2 * N) * sizeof(mp_limb_t));
    if (new_c == NULL || new_bpoly == NULL || new_temp == NULL || new_s_rev == NULL) {
        free(new_c);
        free(new_bpoly);
        free(new_temp);
        free(new_s_rev);
        sparse_bm_state_clear(state);
        return 0;
    }

    state->c = new_c;
    state->bpoly = new_bpoly;
    state->temp = new_temp;
    state->s_rev = new_s_rev;
    if (oldN < N) {
        memset(state->c + oldN + 1, 0, (size_t) (N - oldN) * sizeof(mp_limb_t));
        memset(state->bpoly + oldN + 1, 0, (size_t) (N - oldN) * sizeof(mp_limb_t));
        memset(state->temp + oldN + 1, 0, (size_t) (N - oldN) * sizeof(mp_limb_t));
    }
    state->capacity_N = N;
    if (!state->initialized) {
        state->initialized = 1;
        state->processed_terms = 0;
        state->L = 0;
        state->m = 1;
        state->c_len = 1;
        state->b_len = 1;
        state->b_inv = 1;
        memset(state->c, 0, (size_t) (N + 1) * sizeof(mp_limb_t));
        memset(state->bpoly, 0, (size_t) (N + 1) * sizeof(mp_limb_t));
        state->c[0] = 1;
        state->bpoly[0] = 1;
    }
    return 1;
}

static void sparse_bm_state_reset_sequence(sparse_bm_state_t* state)
{
    if (state == NULL || !state->initialized || state->capacity_N <= 0) {
        return;
    }
    state->processed_terms = 0;
    state->L = 0;
    state->m = 1;
    state->c_len = 1;
    state->b_len = 1;
    state->b_inv = 1;
    memset(state->c, 0, (size_t) (state->capacity_N + 1) * sizeof(mp_limb_t));
    memset(state->bpoly, 0, (size_t) (state->capacity_N + 1) * sizeof(mp_limb_t));
    state->c[0] = 1;
    state->bpoly[0] = 1;
}

static void sparse_scout_reuse_init_empty(sparse_scout_reuse_t* reuse)
{
    if (reuse == NULL) return;
    reuse->valid = 0;
    reuse->nvars = 0;
    reuse->alpha = NULL;
    sparse_bm_state_init_empty(&reuse->bm_state);
}

static void sparse_scout_reuse_clear(sparse_scout_reuse_t* reuse)
{
    if (reuse == NULL) return;
    free(reuse->alpha);
    reuse->alpha = NULL;
    sparse_bm_state_clear(&reuse->bm_state);
    reuse->valid = 0;
    reuse->nvars = 0;
}

static int sparse_scout_reuse_prepare(sparse_scout_reuse_t* reuse, slong nvars, slong budget)
{
    if (reuse == NULL) return 0;
    if (reuse->nvars != nvars || reuse->alpha == NULL) {
        mp_limb_t* alpha = (mp_limb_t*) realloc(reuse->alpha, (size_t) nvars * sizeof(mp_limb_t));
        if (alpha == NULL) {
            sparse_scout_reuse_clear(reuse);
            return 0;
        }
        reuse->alpha = alpha;
        reuse->nvars = nvars;
    }
    if (!sparse_bm_state_ensure_capacity(&reuse->bm_state, budget)) {
        sparse_scout_reuse_clear(reuse);
        return 0;
    }
    return 1;
}

void timer_stop(my_timer_t* t)
{
    t->elapsed = ((double)(clock() - t->start)) / CLOCKS_PER_SEC;
}

void timer_print(my_timer_t* t, const char* label)
{
    if (TIMING) {
        printf("  [TIMING] %s: %.6f seconds\n", label, t->elapsed);
    }
}

void poly_mat_init(poly_mat_t* mat, slong rows, slong cols, const nmod_mpoly_ctx_t ctx)
{
    mat->rows = rows;
    mat->cols = cols;
    mat->entries = (nmod_mpoly_struct**) malloc((size_t) rows * sizeof(nmod_mpoly_struct*));
    for (slong i = 0; i < rows; i++) {
        mat->entries[i] = (nmod_mpoly_struct*) malloc((size_t) cols * sizeof(nmod_mpoly_struct));
        for (slong j = 0; j < cols; j++) {
            nmod_mpoly_init(mat->entries[i] + j, ctx);
        }
    }
}

void poly_mat_clear(poly_mat_t* mat, const nmod_mpoly_ctx_t ctx)
{
    if (mat == NULL || mat->entries == NULL) {
        return;
    }

    for (slong i = 0; i < mat->rows; i++) {
        for (slong j = 0; j < mat->cols; j++) {
            nmod_mpoly_clear(mat->entries[i] + j, ctx);
        }
        free(mat->entries[i]);
    }
    free(mat->entries);
    mat->entries = NULL;
    mat->rows = 0;
    mat->cols = 0;
}

void poly_mat_entry_set(poly_mat_t* mat, slong i, slong j, const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_set(mat->entries[i] + j, poly, ctx);
}

void poly_mat_det_2x2(nmod_mpoly_t det, const poly_mat_t* mat, const nmod_mpoly_ctx_t ctx)
{
    if (mat->rows != 2 || mat->cols != 2) {
        nmod_mpoly_zero(det, ctx);
        return;
    }

    nmod_mpoly_t temp1, temp2;
    nmod_mpoly_init(temp1, ctx);
    nmod_mpoly_init(temp2, ctx);

    nmod_mpoly_mul(temp1, mat->entries[0] + 0, mat->entries[1] + 1, ctx);
    nmod_mpoly_mul(temp2, mat->entries[0] + 1, mat->entries[1] + 0, ctx);
    nmod_mpoly_sub(det, temp1, temp2, ctx);

    nmod_mpoly_clear(temp1, ctx);
    nmod_mpoly_clear(temp2, ctx);
}

void poly_mat_det(nmod_mpoly_t det, const poly_mat_t* mat, const nmod_mpoly_ctx_t ctx)
{
    if (mat->rows != mat->cols) {
        nmod_mpoly_zero(det, ctx);
        return;
    }

    slong n = mat->rows;
    if (n == 0) {
        nmod_mpoly_zero(det, ctx);
        return;
    }
    if (n == 1) {
        nmod_mpoly_set(det, mat->entries[0] + 0, ctx);
        return;
    }
    if (n == 2) {
        poly_mat_det_2x2(det, mat, ctx);
        return;
    }

    nmod_mpoly_zero(det, ctx);

    nmod_mpoly_t minor_det, temp;
    nmod_mpoly_init(minor_det, ctx);
    nmod_mpoly_init(temp, ctx);

    for (slong j = 0; j < n; j++) {
        poly_mat_t minor;
        poly_mat_init(&minor, n - 1, n - 1, ctx);

        for (slong i = 1; i < n; i++) {
            slong col_idx = 0;
            for (slong k = 0; k < n; k++) {
                if (k != j) {
                    nmod_mpoly_set(minor.entries[i - 1] + col_idx, mat->entries[i] + k, ctx);
                    col_idx++;
                }
            }
        }

        poly_mat_det(minor_det, &minor, ctx);
        nmod_mpoly_mul(temp, mat->entries[0] + j, minor_det, ctx);

        if ((j & 1) == 0) {
            nmod_mpoly_add(det, det, temp, ctx);
        } else {
            nmod_mpoly_sub(det, det, temp, ctx);
        }

        poly_mat_clear(&minor, ctx);
    }

    nmod_mpoly_clear(minor_det, ctx);
    nmod_mpoly_clear(temp, ctx);
}

slong poly_max_total_degree(const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)
{
    slong max_deg = 0;
    slong nvars = nmod_mpoly_ctx_nvars(ctx);
    slong len = nmod_mpoly_length(poly, ctx);
    ulong* exp = (ulong*) malloc((size_t) nvars * sizeof(ulong));

    for (slong i = 0; i < len; i++) {
        slong total_deg = 0;
        nmod_mpoly_get_term_exp_ui(exp, poly, i, ctx);
        for (slong j = 0; j < nvars; j++) {
            total_deg += (slong) exp[j];
        }
        if (total_deg > max_deg) {
            max_deg = total_deg;
        }
    }

    free(exp);
    return max_deg;
}

void poly_max_degrees_per_var(slong* max_degs, const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)
{
    slong nvars = nmod_mpoly_ctx_nvars(ctx);
    slong len = nmod_mpoly_length(poly, ctx);
    ulong* exp = (ulong*) malloc((size_t) nvars * sizeof(ulong));

    for (slong i = 0; i < nvars; i++) {
        max_degs[i] = 0;
    }

    for (slong i = 0; i < len; i++) {
        nmod_mpoly_get_term_exp_ui(exp, poly, i, ctx);
        for (slong j = 0; j < nvars; j++) {
            if ((slong) exp[j] > max_degs[j]) {
                max_degs[j] = (slong) exp[j];
            }
        }
    }

    free(exp);
}

static void BM_reference(nmod_poly_t C, const mp_limb_t* s, slong N, nmod_t mod)
{
    my_timer_t timer;
    if (DETAILED_TIMING) timer_start(&timer);

    nmod_poly_t B, T;
    mp_limb_t d, b, temp, temp2;
    slong L = 0, k = 1;

    nmod_poly_init(B, mod.n);
    nmod_poly_init(T, mod.n);

    nmod_poly_set_coeff_ui(B, 0, 1);
    nmod_poly_set_coeff_ui(C, 0, 1);
    b = 1;

    for (slong n = 0; n < 2 * N; n++) {
        d = s[n];
        for (slong i = 1; i <= L && i <= n; i++) {
            temp = nmod_poly_get_coeff_ui(C, i);
            temp2 = nmod_mul(temp, s[n - i], mod);
            d = nmod_add(d, temp2, mod);
        }

        if (d == 0) {
            k++;
        } else if (n < 2 * L) {
            nmod_poly_t xkB;
            nmod_poly_init(xkB, mod.n);
            nmod_poly_shift_left(xkB, B, k);
            temp = nmod_inv(b, mod);
            temp = nmod_mul(d, temp, mod);
            nmod_poly_scalar_mul_nmod(xkB, xkB, temp);
            nmod_poly_sub(C, C, xkB);
            nmod_poly_clear(xkB);
            k++;
        } else {
            nmod_poly_set(T, C);

            nmod_poly_t xkB;
            nmod_poly_init(xkB, mod.n);
            nmod_poly_shift_left(xkB, B, k);
            temp = nmod_inv(b, mod);
            temp = nmod_mul(d, temp, mod);
            nmod_poly_scalar_mul_nmod(xkB, xkB, temp);
            nmod_poly_sub(C, C, xkB);
            nmod_poly_clear(xkB);

            nmod_poly_set(B, T);
            L = n + 1 - L;
            k = 1;
            b = d;
        }
    }

    nmod_poly_clear(B);
    nmod_poly_clear(T);

    if (DETAILED_TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "    BM algorithm");
    }
}

static void BM_timed_incremental(nmod_poly_t C,
                                 const mp_limb_t* s,
                                 slong N,
                                 nmod_t mod,
                                 sparse_bm_state_t* state,
                                 sparse_interpolation_timing_t* timing)
{
    my_timer_t timer;
    if (DETAILED_TIMING) timer_start(&timer);

    if (N <= 0) {
        nmod_poly_zero(C);
        nmod_poly_set_coeff_ui(C, 0, 1);
        return;
    }

    if (state == NULL || !sparse_bm_state_ensure_capacity(state, N)) {
        BM_reference(C, s, N, mod);
        return;
    }

    for (slong i = 0; i < 2 * N; i++) {
        state->s_rev[i] = s[2 * N - 1 - i];
    }

    ulong mod_bits = FLINT_BIT_COUNT(mod.n);
    int inner_timing = (timing != NULL && sparse_bm_inner_timing_enabled());
    int hist_enabled = (timing != NULL && sparse_bm_hist_enabled());

    for (slong n = state->processed_terms; n < 2 * N; n++) {
        clock_t stage_start = 0;
        mp_limb_t d = s[n];
        if (inner_timing) {
            stage_start = clock();
        }
        if (state->L > 0) {
            d = nmod_add(d,
                         nmod_vec_dot_product_unbalanced(state->c + 1,
                                                         state->s_rev + (2 * N - n),
                                                         (ulong) state->L,
                                                         mod_bits,
                                                         mod_bits,
                                                         mod),
                         mod);
        }
        if (inner_timing) {
            timing->bm_discrepancy += timer_seconds_since(stage_start);
        }

        if (d == 0) {
            state->m++;
            continue;
        }

        if (hist_enabled) {
            timing->bm_update_count++;
            timing->bm_update_len_sum += state->b_len;
            timing->bm_update_hist[sparse_bm_hist_bin(state->b_len)]++;
        }

        clock_t update_total_start = 0;
        if (inner_timing) {
            update_total_start = clock();
            stage_start = update_total_start;
        }
        mp_limb_t coef = nmod_mul(d, state->b_inv, mod);
        slong old_c_len = state->c_len;
        int promote = (n >= 2 * state->L);

        if (promote) {
            memcpy(state->temp, state->c, (size_t) old_c_len * sizeof(mp_limb_t));
        }
        if (inner_timing) {
            timing->bm_update_prepare += timer_seconds_since(stage_start);
            stage_start = clock();
        }

        slong update_len = state->b_len + state->m;
        if (update_len > state->c_len) {
            state->c_len = update_len;
        }
        sparse_bm_axpy(state->c + state->m, state->bpoly, state->b_len, nmod_neg(coef, mod), mod);
        if (inner_timing) {
            timing->bm_update_axpy += timer_seconds_since(stage_start);
            stage_start = clock();
        }
        if (update_len >= old_c_len) {
            while (state->c_len > 1 && state->c[state->c_len - 1] == 0) {
                state->c_len--;
            }
        }
        if (inner_timing) {
            timing->bm_update_trim += timer_seconds_since(stage_start);
            stage_start = clock();
        }

        if (promote) {
            if (hist_enabled) {
                timing->bm_promote_count++;
            }
            mp_limb_t* swap = state->bpoly;
            state->bpoly = state->temp;
            state->temp = swap;
            state->b_len = old_c_len;
            state->L = n + 1 - state->L;
            state->m = 1;
            state->b_inv = nmod_inv(d, mod);
        } else {
            state->m++;
        }

        if (inner_timing) {
            timing->bm_update_promote += timer_seconds_since(stage_start);
        }
        if (inner_timing) {
            timing->bm_update += timer_seconds_since(update_total_start);
        }
    }
    state->processed_terms = 2 * N;

    clock_t stage_start = 0;
    if (inner_timing) {
        stage_start = clock();
    }
    nmod_poly_zero(C);
    for (slong i = 0; i < state->c_len; i++) {
        if (state->c[i] != 0) {
            nmod_poly_set_coeff_ui(C, i, state->c[i]);
        }
    }
    if (inner_timing) {
        timing->bm_finalize += timer_seconds_since(stage_start);
    }

    if (DETAILED_TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "    BM algorithm");
    }
}

static void BM_timed(nmod_poly_t C,
                     const mp_limb_t* s,
                     slong N,
                     nmod_t mod,
                     sparse_interpolation_timing_t* timing)
{
    sparse_bm_state_t state;
    sparse_bm_state_init_empty(&state);
    BM_timed_incremental(C, s, N, mod, &state, timing);
    sparse_bm_state_clear(&state);
}

void BM(nmod_poly_t C, const mp_limb_t* s, slong N, nmod_t mod)
{
    BM_timed(C, s, N, mod, NULL);
}

void Vinvert(mp_limb_t* c1, const nmod_poly_t c, const mp_limb_t* v,
             const mp_limb_t* a, slong n, nmod_t mod)
{
    my_timer_t timer;
    if (DETAILED_TIMING) timer_start(&timer);

    nmod_poly_t d, q, q1, q2;
    mp_limb_t temp;

    nmod_poly_init(d, mod.n);
    nmod_poly_init(q, mod.n);
    nmod_poly_init(q1, mod.n);
    nmod_poly_init(q2, mod.n);

    nmod_poly_zero(d);
    for (slong i = 0; i < n; i++) {
        nmod_poly_set_coeff_ui(d, n - i, a[i]);
    }

    nmod_poly_mul(q, c, d);

    nmod_poly_zero(q1);
    for (slong i = 0; i < n; i++) {
        temp = nmod_poly_get_coeff_ui(q, 2 * n - i);
        nmod_poly_set_coeff_ui(q1, n - 1 - i, temp);
    }

    nmod_poly_derivative(q2, c);

    for (slong i = 0; i < n; i++) {
        mp_limb_t q1_val = nmod_poly_evaluate_nmod(q1, v[i]);
        mp_limb_t q2_val = nmod_poly_evaluate_nmod(q2, v[i]);

        if (q2_val != 0) {
            temp = nmod_inv(q2_val, mod);
            c1[i] = nmod_mul(q1_val, temp, mod);
        } else {
            c1[i] = a[i];
        }
    }

    nmod_poly_clear(d);
    nmod_poly_clear(q);
    nmod_poly_clear(q1);
    nmod_poly_clear(q2);

    if (DETAILED_TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "    Vinvert");
    }
}

static int precompute_locator_inverse_denominators(mp_limb_t* denom_inv,
                                                   const nmod_poly_t c,
                                                   const mp_limb_t* v,
                                                   slong n,
                                                   nmod_t mod)
{
    nmod_poly_t deriv;
    nmod_poly_init(deriv, mod.n);
    nmod_poly_derivative(deriv, c);

    for (slong i = 0; i < n; i++) {
        mp_limb_t q2_val = nmod_poly_evaluate_nmod(deriv, v[i]);
        if (q2_val == 0) {
            nmod_poly_clear(deriv);
            return 0;
        }
        denom_inv[i] = nmod_inv(q2_val, mod);
    }

    nmod_poly_clear(deriv);
    return 1;
}

static int sparse_multipoint_tree_init(sparse_multipoint_tree_t* tree,
                                       const mp_limb_t* points,
                                       slong npoints,
                                       nmod_t mod)
{
    slong leaf_count = 1;
    sparse_multipoint_tree_init_empty(tree);

    if (npoints <= 0) {
        return 1;
    }

    while (leaf_count < npoints) {
        leaf_count <<= 1;
    }

    tree->npoints = npoints;
    tree->leaf_count = leaf_count;
    tree->node_count = 2 * leaf_count;
    tree->points = (mp_limb_t*) malloc((size_t) npoints * sizeof(mp_limb_t));
    tree->prod = (nmod_poly_struct*) malloc((size_t) tree->node_count * sizeof(nmod_poly_struct));
    tree->rem = (nmod_poly_struct*) malloc((size_t) tree->node_count * sizeof(nmod_poly_struct));
    if (tree->points == NULL || tree->prod == NULL || tree->rem == NULL) {
        free(tree->points);
        free(tree->prod);
        free(tree->rem);
        sparse_multipoint_tree_init_empty(tree);
        return 0;
    }

    memcpy(tree->points, points, (size_t) npoints * sizeof(mp_limb_t));
    for (slong i = 1; i < tree->node_count; i++) {
        nmod_poly_init(tree->prod + i, mod.n);
        nmod_poly_init(tree->rem + i, mod.n);
    }

    for (slong i = 0; i < leaf_count; i++) {
        slong leaf = leaf_count + i;
        if (i < npoints) {
            nmod_poly_zero(tree->prod + leaf);
            nmod_poly_set_coeff_ui(tree->prod + leaf, 1, 1);
            nmod_poly_set_coeff_ui(tree->prod + leaf, 0, nmod_neg(points[i], mod));
        } else {
            nmod_poly_one(tree->prod + leaf);
        }
    }

    for (slong node = leaf_count - 1; node >= 1; node--) {
        nmod_poly_mul(tree->prod + node,
                      tree->prod + (2 * node),
                      tree->prod + (2 * node + 1));
    }

    return 1;
}

static void sparse_multipoint_tree_clear(sparse_multipoint_tree_t* tree)
{
    if (tree == NULL) return;
    if (tree->prod != NULL && tree->rem != NULL) {
        for (slong i = 1; i < tree->node_count; i++) {
            nmod_poly_clear(tree->prod + i);
            nmod_poly_clear(tree->rem + i);
        }
    }
    free(tree->points);
    free(tree->prod);
    free(tree->rem);
    sparse_multipoint_tree_init_empty(tree);
}

static void sparse_multipoint_tree_evaluate(mp_limb_t* values,
                                            sparse_multipoint_tree_t* tree,
                                            const nmod_poly_t poly,
                                            nmod_t mod)
{
    nmod_poly_t quotient;
    nmod_poly_init(quotient, mod.n);

    nmod_poly_set(tree->rem + 1, poly);
    for (slong node = 1; node < tree->leaf_count; node++) {
        slong left = 2 * node;
        slong right = left + 1;
        slong rem_deg = nmod_poly_degree(tree->rem + node);

        if (nmod_poly_degree(tree->prod + left) > 0) {
            if (rem_deg < nmod_poly_degree(tree->prod + left)) {
                nmod_poly_set(tree->rem + left, tree->rem + node);
            } else {
                nmod_poly_divrem(quotient, tree->rem + left, tree->rem + node, tree->prod + left);
            }
        }

        if (nmod_poly_degree(tree->prod + right) > 0) {
            if (rem_deg < nmod_poly_degree(tree->prod + right)) {
                nmod_poly_set(tree->rem + right, tree->rem + node);
            } else {
                nmod_poly_divrem(quotient, tree->rem + right, tree->rem + node, tree->prod + right);
            }
        }
    }

    for (slong i = 0; i < tree->npoints; i++) {
        slong leaf = tree->leaf_count + i;
        values[i] = nmod_poly_degree(tree->rem + leaf) >= 0
                    ? nmod_poly_get_coeff_ui(tree->rem + leaf, 0)
                    : 0;
    }

    nmod_poly_clear(quotient);
}

static void Vinvert_preconditioned(mp_limb_t* c1,
                                   const nmod_poly_t c,
                                   const mp_limb_t* v,
                                   sparse_multipoint_tree_t* eval_tree,
                                   const mp_limb_t* denom_inv,
                                   const mp_limb_t* a,
                                   slong n,
                                   nmod_t mod,
                                   sparse_interpolation_timing_t* timing)
{
    nmod_poly_t d, q, q1;
    clock_t stage_start;

    nmod_poly_init(d, mod.n);
    nmod_poly_init(q, mod.n);
    nmod_poly_init(q1, mod.n);

    stage_start = clock();
    nmod_poly_zero(d);
    for (slong i = 0; i < n; i++) {
        nmod_poly_set_coeff_ui(d, n - i, a[i]);
    }
    if (timing != NULL) {
        timing->vinvert_build += timer_seconds_since(stage_start);
    }

    stage_start = clock();
    nmod_poly_mul(q, c, d);
    if (timing != NULL) {
        timing->vinvert_mul += timer_seconds_since(stage_start);
    }

    stage_start = clock();
    nmod_poly_zero(q1);
    for (slong i = 0; i < n; i++) {
        nmod_poly_set_coeff_ui(q1, n - 1 - i, nmod_poly_get_coeff_ui(q, 2 * n - i));
    }
    if (timing != NULL) {
        timing->vinvert_extract += timer_seconds_since(stage_start);
    }

    stage_start = clock();
    if (eval_tree != NULL && eval_tree->npoints == n) {
        sparse_multipoint_tree_evaluate(c1, eval_tree, q1, mod);
    } else {
        for (slong i = 0; i < n; i++) {
            c1[i] = nmod_poly_evaluate_nmod(q1, v[i]);
        }
    }
    for (slong i = 0; i < n; i++) {
        c1[i] = nmod_mul(c1[i], denom_inv[i], mod);
    }
    if (timing != NULL) {
        timing->vinvert_eval += timer_seconds_since(stage_start);
    }

    nmod_poly_clear(d);
    nmod_poly_clear(q);
    nmod_poly_clear(q1);
}

static int compare_limb_values(const void* a, const void* b)
{
    const mp_limb_t va = *(const mp_limb_t*) a;
    const mp_limb_t vb = *(const mp_limb_t*) b;
    if (va < vb) return -1;
    if (va > vb) return 1;
    return 0;
}

static void sparse_recovery_info_reset(sparse_recovery_info_t* info)
{
    if (info == NULL) {
        return;
    }

    info->used = 1;
    info->likely_underbudget = 0;
    info->attempts = 0;
    info->attempt_limit = 0;
    info->budget = 0;
    info->last_relation_degree = -1;
    info->last_reconstructed_terms = 0;
    info->workspace_bytes = 0;
    memset(&info->probe_timing, 0, sizeof(info->probe_timing));
    memset(&info->interpolation_timing, 0, sizeof(info->interpolation_timing));
    info->last_reason[0] = '\0';
}

static void sparse_recovery_info_set_reason(sparse_recovery_info_t* info, const char* fmt, ...)
{
    if (info == NULL) {
        return;
    }

    va_list args;
    va_start(args, fmt);
    vsnprintf(info->last_reason, sizeof(info->last_reason), fmt, args);
    va_end(args);
}

static slong sparse_get_env_slong(const char* name, slong default_value, slong min_value)
{
    const char* value = getenv(name);
    if (value == NULL || *value == '\0') {
        return default_value;
    }

    char* endptr = NULL;
    long long parsed = strtoll(value, &endptr, 10);
    if (endptr == value || *endptr != '\0') {
        return default_value;
    }
    if (parsed < (long long) min_value) {
        return min_value;
    }
    if (parsed > (long long) LONG_MAX) {
        return LONG_MAX;
    }
    return (slong) parsed;
}

static int sparse_debug_level(void)
{
    static int initialized = 0;
    static int level = 0;
    if (!initialized) {
        level = (int) sparse_get_env_slong("DIXON_SPARSE_DEBUG", 0, 0);
        initialized = 1;
    }
    return level;
}

static int sparse_bm_inner_timing_enabled(void)
{
    static int initialized = 0;
    static int enabled = 0;
    if (!initialized) {
        enabled = (int) sparse_get_env_slong("DIXON_SPARSE_BM_INNER_TIMING", 0, 0);
        initialized = 1;
    }
    return enabled;
}

static int sparse_bm_hist_enabled(void)
{
    static int initialized = 0;
    static int enabled = 0;
    if (!initialized) {
        enabled = (int) sparse_get_env_slong("DIXON_SPARSE_BM_HIST", 0, 0);
        initialized = 1;
    }
    return enabled;
}

static double timer_seconds_since(clock_t start)
{
    return ((double) (clock() - start)) / CLOCKS_PER_SEC;
}

static void sparse_probe_timing_add(sparse_probe_timing_t* dst, const sparse_probe_timing_t* src)
{
    if (dst == NULL || src == NULL) return;
    dst->total += src->total;
    dst->prepare += src->prepare;
    dst->eval_entries += src->eval_entries;
    dst->det += src->det;
    dst->inv += src->inv;
    dst->direct_trace += src->direct_trace;
    dst->adjugate_trace += src->adjugate_trace;
    dst->point_count += src->point_count;
    dst->entry_term_count += src->entry_term_count;
    dst->inverse_count += src->inverse_count;
    dst->singular_count += src->singular_count;
}

static void sparse_interpolation_timing_add(sparse_interpolation_timing_t* dst,
                                            const sparse_interpolation_timing_t* src)
{
    if (dst == NULL || src == NULL) return;
    dst->total += src->total;
    dst->probe += src->probe;
    dst->bm += src->bm;
    dst->bm_discrepancy += src->bm_discrepancy;
    dst->bm_update += src->bm_update;
    dst->bm_finalize += src->bm_finalize;
    dst->bm_update_prepare += src->bm_update_prepare;
    dst->bm_update_axpy += src->bm_update_axpy;
    dst->bm_update_trim += src->bm_update_trim;
    dst->bm_update_promote += src->bm_update_promote;
    dst->bm_update_count += src->bm_update_count;
    dst->bm_promote_count += src->bm_promote_count;
    dst->bm_update_len_sum += src->bm_update_len_sum;
    for (int i = 0; i < SPARSE_BM_HIST_BINS; i++) {
        dst->bm_update_hist[i] += src->bm_update_hist[i];
    }
    dst->root_factor += src->root_factor;
    dst->root_squarefree += src->root_squarefree;
    dst->root_frobenius += src->root_frobenius;
    dst->root_split += src->root_split;
    dst->root_sort += src->root_sort;
    dst->root_fallback_factor += src->root_fallback_factor;
    dst->vinvert += src->vinvert;
    dst->vinvert_denom += src->vinvert_denom;
    dst->vinvert_tree_build += src->vinvert_tree_build;
    dst->vinvert_build += src->vinvert_build;
    dst->vinvert_mul += src->vinvert_mul;
    dst->vinvert_extract += src->vinvert_extract;
    dst->vinvert_eval += src->vinvert_eval;
    dst->reconstruct += src->reconstruct;
}

static slong sparse_root_split_attempt_limit(slong degree)
{
    slong env_attempts = sparse_get_env_slong("DIXON_SPARSE_ROOT_SPLIT_ATTEMPTS", -1, 1);
    if (env_attempts > 0) {
        return env_attempts;
    }

    slong attempts = SPARSE_FAST_ROOT_SPLIT_ATTEMPTS;
    slong deg = degree;
    while (deg > 1) {
        deg >>= 1;
        attempts++;
    }
    if (attempts > 32) attempts = 32;
    return attempts;
}

static void sparse_multipoint_tree_init_empty(sparse_multipoint_tree_t* tree)
{
    if (tree == NULL) return;
    tree->npoints = 0;
    tree->leaf_count = 0;
    tree->node_count = 0;
    tree->points = NULL;
    tree->prod = NULL;
    tree->rem = NULL;
}

static slong sparse_get_attempt_limit(slong budget)
{
    slong env_attempts = sparse_get_env_slong("DIXON_SPARSE_ATTEMPTS", -1, 1);
    if (env_attempts > 0) {
        return env_attempts;
    }

    slong attempts = SPARSE_BASE_RECOVERY_ATTEMPTS;
    if (budget >= 4096) attempts += 8;
    if (budget >= 16384) attempts += 8;
    if (budget >= 65536) attempts += 8;
    return attempts;
}

static slong sparse_get_scout_attempt_limit(slong budget, slong upper_bound)
{
    slong env_scout = sparse_get_env_slong("DIXON_SPARSE_SCOUT_ATTEMPTS", 2, 1);
    slong full = sparse_get_attempt_limit(budget);
    if (budget >= upper_bound) {
        return full;
    }
    if (budget >= upper_bound / 4) {
        return FLINT_MIN(full, env_scout + 2);
    }
    return FLINT_MIN(full, env_scout);
}

static slong sparse_round_budget_up_pow2(slong target, slong upper_bound)
{
    slong budget = 1;
    if (target < 1) target = 1;
    while (budget < target && budget < upper_bound) {
        budget <<= 1;
    }
    if (budget > upper_bound) {
        budget = upper_bound;
    }
    return budget;
}

static slong sparse_get_start_budget(slong upper_bound,
                                     const slong* max_degs_per_var,
                                     slong nvars)
{
    slong default_start = SPARSE_DEFAULT_START_BUDGET;
    if (upper_bound > SPARSE_DEFAULT_START_BUDGET) {
        slong active_vars = 0;
        slong degree_sum = 0;
        slong shift_target;
        slong shape_target;
        slong target;

        if (max_degs_per_var != NULL) {
            for (slong i = 0; i < nvars; i++) {
                if (max_degs_per_var[i] > 0) {
                    active_vars++;
                }
                if (degree_sum < upper_bound) {
                    slong add = max_degs_per_var[i];
                    if (add > upper_bound - degree_sum) {
                        degree_sum = upper_bound;
                    } else {
                        degree_sum += add;
                    }
                }
            }
        } else {
            active_vars = nvars;
        }

        shift_target = upper_bound;
        {
            slong shift = 3 + active_vars / 2;
            if (shift > 12) shift = 12;
            while (shift > 0 && shift_target > SPARSE_DEFAULT_START_BUDGET) {
                shift_target >>= 1;
                shift--;
            }
        }

        shape_target = degree_sum + 1;
        if (active_vars > 0) {
            slong scale_shift = active_vars;
            slong scale;
            if (scale_shift > 5) scale_shift = 5;
            scale = ((slong) 1) << scale_shift;
            if (shape_target > upper_bound / scale) {
                shape_target = upper_bound;
            } else {
                shape_target *= scale;
            }
        }

        if (active_vars <= 3) {
            target = shift_target;
        } else if (active_vars <= 6) {
            target = FLINT_MAX(shape_target, shift_target >> 1);
        } else {
            target = shape_target;
        }

        if (target < SPARSE_DEFAULT_START_BUDGET) target = SPARSE_DEFAULT_START_BUDGET;
        default_start = sparse_round_budget_up_pow2(target, upper_bound);
    }
    slong start = sparse_get_env_slong("DIXON_SPARSE_START_BUDGET", default_start, 1);
    if (start > upper_bound) {
        start = upper_bound;
    }
    return start;
}

static size_t sparse_workspace_bytes(slong nvars, slong Tbound)
{
    if (Tbound <= 0) {
        return 0;
    }

    size_t limb_bytes = sizeof(mp_limb_t);
    size_t ulong_bytes = sizeof(ulong);
    size_t total = 0;
    total += (size_t) nvars * limb_bytes;
    total += (size_t) (2 * Tbound) * limb_bytes;
    total += (size_t) nvars * (size_t) (2 * Tbound) * limb_bytes;
    total += (size_t) Tbound * limb_bytes;
    total += (size_t) Tbound * limb_bytes;
    total += (size_t) Tbound * limb_bytes;
    total += (size_t) nvars * (size_t) Tbound * limb_bytes;
    total += (size_t) nvars * ulong_bytes;
    return total;
}

static mp_limb_t nmod_rand_nonzero(nmod_t mod)
{
    mp_limb_t x;
    do {
        x = n_randint(global_state, mod.n);
    } while (x == 0);
    return x;
}

static void build_locator_from_relation(nmod_poly_t locator, const nmod_poly_t relation, slong degree)
{
    nmod_poly_zero(locator);
    nmod_poly_set_coeff_ui(locator, degree, 1);
    for (slong i = 1; i <= degree; i++) {
        nmod_poly_set_coeff_ui(locator, degree - i, nmod_poly_get_coeff_ui(relation, i));
    }
}

static void build_locator_from_roots(nmod_poly_t locator, const mp_limb_t* roots, slong t, nmod_t mod)
{
    nmod_poly_one(locator);
    for (slong i = 0; i < t; i++) {
        nmod_poly_t factor;
        nmod_poly_init(factor, mod.n);
        nmod_poly_set_coeff_ui(factor, 1, 1);
        nmod_poly_set_coeff_ui(factor, 0, nmod_neg(roots[i], mod));
        nmod_poly_mul(locator, locator, factor);
        nmod_poly_clear(factor);
    }
}

static int sparse_extract_linear_root(mp_limb_t* root, const nmod_poly_t poly, nmod_t mod)
{
    mp_limb_t a, b;

    if (nmod_poly_degree(poly) != 1) {
        return 0;
    }

    a = nmod_poly_get_coeff_ui(poly, 1);
    if (a == 0) {
        return 0;
    }

    b = nmod_poly_get_coeff_ui(poly, 0);
    *root = nmod_neg(nmod_mul(b, nmod_inv(a, mod), mod), mod);
    return 1;
}

static int sparse_collect_linear_roots_recursive(mp_limb_t* roots,
                                                 slong* count,
                                                 const nmod_poly_t poly,
                                                 nmod_t mod)
{
    slong deg = nmod_poly_degree(poly);
    slong base_count = *count;
    if (deg <= 0) {
        return 1;
    }

    if (deg == 1) {
        return sparse_extract_linear_root(roots + (*count)++, poly, mod);
    }

    nmod_poly_t g, h, split, quotient;
    int ok = 0;

    nmod_poly_init(g, mod.n);
    nmod_poly_init(h, mod.n);
    nmod_poly_init(split, mod.n);
    nmod_poly_init(quotient, mod.n);

    for (slong attempt = 0; attempt < sparse_root_split_attempt_limit(deg); attempt++) {
        *count = base_count;
        nmod_poly_zero(g);
        nmod_poly_set_coeff_ui(g, 1, 1);
        nmod_poly_set_coeff_ui(g, 0, n_randint(global_state, mod.n));
        nmod_poly_powmod_ui_binexp(h, g, (mod.n - 1) / 2, poly);
        nmod_poly_set_coeff_ui(h, 0, n_submod(nmod_poly_get_coeff_ui(h, 0), 1, mod.n));
        nmod_poly_gcd(split, h, poly);

        slong split_deg = nmod_poly_degree(split);
        if (split_deg <= 0 || split_deg >= deg) {
            continue;
        }

        nmod_poly_div(quotient, poly, split);
        ok = sparse_collect_linear_roots_recursive(roots, count, split, mod) &&
             sparse_collect_linear_roots_recursive(roots, count, quotient, mod);
        if (ok) {
            break;
        }
    }

    if (!ok) {
        *count = base_count;
    }

    nmod_poly_clear(g);
    nmod_poly_clear(h);
    nmod_poly_clear(split);
    nmod_poly_clear(quotient);
    return ok;
}

static int get_simple_roots_fast(mp_limb_t* roots,
                                 const nmod_poly_t poly,
                                 slong expected_degree,
                                 nmod_t mod,
                                 sparse_interpolation_timing_t* timing)
{
    nmod_poly_t deriv, square_gcd, x, x_to_p, frobenius, root_poly;
    slong count = 0;
    int ok = 0;
    clock_t stage_start;

    if (nmod_poly_degree(poly) != expected_degree) {
        return 0;
    }

    ensure_global_state();

    nmod_poly_init(deriv, mod.n);
    nmod_poly_init(square_gcd, mod.n);
    nmod_poly_init(x, mod.n);
    nmod_poly_init(x_to_p, mod.n);
    nmod_poly_init(frobenius, mod.n);
    nmod_poly_init(root_poly, mod.n);

    stage_start = clock();
    nmod_poly_derivative(deriv, poly);
    nmod_poly_gcd(square_gcd, poly, deriv);
    if (timing != NULL) {
        timing->root_squarefree += timer_seconds_since(stage_start);
    }
    if (nmod_poly_degree(square_gcd) > 0) {
        goto cleanup;
    }

    stage_start = clock();
    nmod_poly_set_coeff_ui(x, 1, 1);
    nmod_poly_powmod_ui_binexp(x_to_p, x, mod.n, poly);
    nmod_poly_sub(frobenius, x_to_p, x);
    nmod_poly_gcd(root_poly, poly, frobenius);
    if (timing != NULL) {
        timing->root_frobenius += timer_seconds_since(stage_start);
    }
    if (nmod_poly_degree(root_poly) != expected_degree) {
        goto cleanup;
    }

    stage_start = clock();
    ok = sparse_collect_linear_roots_recursive(roots, &count, root_poly, mod);
    if (timing != NULL) {
        timing->root_split += timer_seconds_since(stage_start);
    }
    if (!ok || count != expected_degree) {
        ok = 0;
        goto cleanup;
    }

    stage_start = clock();
    qsort(roots, (size_t) count, sizeof(mp_limb_t), compare_limb_values);
    for (slong i = 1; i < count; i++) {
        if (roots[i] == roots[i - 1]) {
            ok = 0;
            goto cleanup;
        }
    }
    if (timing != NULL) {
        timing->root_sort += timer_seconds_since(stage_start);
    }

cleanup:
    nmod_poly_clear(deriv);
    nmod_poly_clear(square_gcd);
    nmod_poly_clear(x);
    nmod_poly_clear(x_to_p);
    nmod_poly_clear(frobenius);
    nmod_poly_clear(root_poly);
    return ok;
}

static int get_simple_roots_strict(mp_limb_t* roots,
                                   const nmod_poly_t poly,
                                   slong expected_degree,
                                   nmod_t mod,
                                   sparse_interpolation_timing_t* timing)
{
    if (expected_degree == 0) {
        return 1;
    }

    if (expected_degree >= SPARSE_FAST_ROOT_THRESHOLD &&
        get_simple_roots_fast(roots, poly, expected_degree, mod, timing)) {
        return 1;
    }

    nmod_poly_factor_t fac;
    nmod_poly_factor_init(fac);
    clock_t stage_start = clock();
    nmod_poly_factor(fac, poly);

    slong count = 0;
    int ok = 1;

    for (slong i = 0; i < fac->num && ok; i++) {
        if (fac->exp[i] != 1 || nmod_poly_degree(fac->p + i) != 1) {
            ok = 0;
            break;
        }

        mp_limb_t a = nmod_poly_get_coeff_ui(fac->p + i, 1);
        mp_limb_t b = nmod_poly_get_coeff_ui(fac->p + i, 0);
        if (a == 0) {
            ok = 0;
            break;
        }

        roots[count++] = nmod_neg(nmod_mul(b, nmod_inv(a, mod), mod), mod);
    }

    if (count != expected_degree) {
        ok = 0;
    }

    if (ok) {
        qsort(roots, (size_t) count, sizeof(mp_limb_t), compare_limb_values);
        for (slong i = 1; i < count; i++) {
            if (roots[i] == roots[i - 1]) {
                ok = 0;
                break;
            }
        }
    }

    if (timing != NULL) {
        timing->root_fallback_factor += timer_seconds_since(stage_start);
    }

    nmod_poly_factor_clear(fac);
    return ok;
}

static slong sat_add_slong(slong a, slong b, slong cap)
{
    if (a >= cap || b >= cap) return cap;
    if (b > cap - a) return cap;
    return a + b;
}

static slong sat_mul_slong(slong a, slong b, slong cap)
{
    if (a == 0 || b == 0) return 0;
    if (a >= cap || b >= cap) return cap;
    if (a > cap / b) return cap;
    return a * b;
}

static mp_limb_t evaluate_root_from_alpha(const mp_limb_t* alpha,
                                          const ulong* exp,
                                          slong nvars,
                                          nmod_t mod)
{
    mp_limb_t value = 1;
    for (slong i = 0; i < nvars; i++) {
        ulong e = exp[i];
        if (e == 0) continue;
        if (e == 1) {
            value = nmod_mul(value, alpha[i], mod);
        } else {
            value = nmod_mul(value, nmod_pow_ui(alpha[i], e, mod), mod);
        }
    }
    return value;
}

static int geometric_poly_cache_init(geometric_poly_cache_t* cache,
                                     const nmod_mpoly_t poly,
                                     const mp_limb_t* alpha,
                                     slong nvars,
                                     nmod_t mod,
                                     ulong* exp_scratch,
                                     const nmod_mpoly_ctx_t mctx)
{
    cache->len = nmod_mpoly_length(poly, mctx);
    cache->coeffs = NULL;
    cache->roots = NULL;
    cache->powers = NULL;
    cache->exps = NULL;

    if (cache->len <= 0) {
        return 1;
    }

    cache->coeffs = (mp_limb_t*) malloc((size_t) cache->len * sizeof(mp_limb_t));
    cache->roots = (mp_limb_t*) malloc((size_t) cache->len * sizeof(mp_limb_t));
    cache->powers = (mp_limb_t*) malloc((size_t) cache->len * sizeof(mp_limb_t));
    cache->exps = (ulong*) malloc((size_t) cache->len * (size_t) nvars * sizeof(ulong));
    if (cache->coeffs == NULL || cache->roots == NULL || cache->powers == NULL || cache->exps == NULL) {
        return 0;
    }

    for (slong term = 0; term < cache->len; term++) {
        nmod_mpoly_get_term_exp_ui(exp_scratch, poly, term, mctx);
        cache->coeffs[term] = nmod_mpoly_get_term_coeff_ui(poly, term, mctx);
        cache->roots[term] = evaluate_root_from_alpha(alpha, exp_scratch, nvars, mod);
        cache->powers[term] = 1;
        memcpy(cache->exps + (size_t) term * (size_t) nvars, exp_scratch, (size_t) nvars * sizeof(ulong));
    }

    return 1;
}

static void geometric_poly_cache_clear(geometric_poly_cache_t* cache)
{
    if (cache == NULL) return;
    free(cache->coeffs);
    free(cache->roots);
    free(cache->powers);
    free(cache->exps);
    cache->coeffs = NULL;
    cache->roots = NULL;
    cache->powers = NULL;
    cache->exps = NULL;
    cache->len = 0;
}

static void geometric_poly_cache_accumulate_and_advance(mp_limb_t* value,
                                                        mp_limb_t* scaled_derivs,
                                                        geometric_poly_cache_t* cache,
                                                        slong nvars,
                                                        nmod_t mod)
{
    *value = 0;
    for (slong i = 0; i < nvars; i++) {
        scaled_derivs[i] = 0;
    }

    for (slong term = 0; term < cache->len; term++) {
        mp_limb_t term_value = nmod_mul(cache->coeffs[term], cache->powers[term], mod);
        *value = nmod_add(*value, term_value, mod);
        ulong* exp = cache->exps + (size_t) term * (size_t) nvars;
        for (slong i = 0; i < nvars; i++) {
            ulong e = exp[i];
            if (e == 0) continue;
            scaled_derivs[i] = nmod_add(scaled_derivs[i],
                                        nmod_mul(term_value, e % mod.n, mod),
                                        mod);
        }
        cache->powers[term] = nmod_mul(cache->powers[term], cache->roots[term], mod);
    }
}

static mp_limb_t nmod_mat_minor_det(const nmod_mat_t A, slong remove_row, slong remove_col, nmod_t mod)
{
    slong n = nmod_mat_nrows(A);
    if (n == 1) {
        return 1;
    }

    nmod_mat_t minor;
    nmod_mat_init(minor, n - 1, n - 1, mod.n);

    slong rr = 0;
    for (slong i = 0; i < n; i++) {
        if (i == remove_row) continue;
        slong cc = 0;
        for (slong j = 0; j < n; j++) {
            if (j == remove_col) continue;
            nmod_mat_entry(minor, rr, cc) = nmod_mat_entry(A, i, j);
            cc++;
        }
        rr++;
    }

    mp_limb_t det = nmod_mat_det(minor);
    nmod_mat_clear(minor);
    return det;
}

static void nmod_mat_adjugate_classical(nmod_mat_t adj, const nmod_mat_t A, nmod_t mod)
{
    slong n = nmod_mat_nrows(A);
    if (n == 0) {
        return;
    }
    if (n == 1) {
        nmod_mat_entry(adj, 0, 0) = 1;
        return;
    }

    for (slong r = 0; r < n; r++) {
        for (slong c = 0; c < n; c++) {
            mp_limb_t cofactor = nmod_mat_minor_det(A, c, r, mod);
            if (((r + c) & 1) != 0) {
                cofactor = nmod_neg(cofactor, mod);
            }
            nmod_mat_entry(adj, r, c) = cofactor;
        }
    }
}

static int explicit_poly_probe_ctx_init(explicit_poly_probe_ctx_t* ctx,
                                        const nmod_mpoly_t f,
                                        slong nvars,
                                        const nmod_mpoly_ctx_t mctx)
{
    ctx->f = f;
    ctx->nvars = nvars;
    ctx->mctx = mctx;
    ctx->timing = NULL;
    return 1;
}

static void explicit_poly_probe_ctx_clear(explicit_poly_probe_ctx_t* ctx)
{
    (void) ctx;
}

static int explicit_poly_probe_sequences(mp_limb_t* a,
                                         mp_limb_t* h,
                                         const mp_limb_t* alpha,
                                         slong nvars,
                                         slong Tbound,
                                         void* user_ctx,
                                         nmod_t mod)
{
    explicit_poly_probe_ctx_t* ctx = (explicit_poly_probe_ctx_t*) user_ctx;
    slong seq_len = 2 * Tbound;
    mp_limb_t* scaled_derivs = (mp_limb_t*) malloc((size_t) nvars * sizeof(mp_limb_t));
    ulong* exp_scratch = (ulong*) malloc((size_t) nvars * sizeof(ulong));
    geometric_poly_cache_t cache;
    clock_t total_start = clock();
    if (scaled_derivs == NULL || exp_scratch == NULL) {
        free(scaled_derivs);
        free(exp_scratch);
        return 0;
    }
    if (!geometric_poly_cache_init(&cache, ctx->f, alpha, nvars, mod, exp_scratch, ctx->mctx)) {
        free(scaled_derivs);
        free(exp_scratch);
        geometric_poly_cache_clear(&cache);
        return 0;
    }
    if (ctx->timing != NULL) {
        ctx->timing->prepare += timer_seconds_since(total_start);
        ctx->timing->entry_term_count += cache.len;
    }

    for (slong j = 0; j < seq_len; j++) {
        clock_t eval_start = clock();
        geometric_poly_cache_accumulate_and_advance(a + j, scaled_derivs, &cache, nvars, mod);
        if (ctx->timing != NULL) {
            ctx->timing->eval_entries += timer_seconds_since(eval_start);
            ctx->timing->point_count++;
        }
        for (slong i = 0; i < nvars; i++) {
            h[i * seq_len + j] = scaled_derivs[i];
        }
    }

    free(scaled_derivs);
    free(exp_scratch);
    geometric_poly_cache_clear(&cache);
    if (ctx->timing != NULL) {
        ctx->timing->total += timer_seconds_since(total_start);
    }
    return 1;
}

static int determinant_probe_ctx_init(determinant_probe_ctx_t* ctx,
                                      const poly_mat_t* A,
                                      slong nvars,
                                      const nmod_mpoly_ctx_t mctx)
{
    ctx->A = A;
    ctx->nvars = nvars;
    ctx->matrix_size = A->rows;
    ctx->mctx = mctx;
    ctx->timing = NULL;
    return 1;
}

static void determinant_probe_ctx_clear(determinant_probe_ctx_t* ctx)
{
    (void) ctx;
}

static int determinant_probe_sequences(mp_limb_t* a,
                                       mp_limb_t* h,
                                       const mp_limb_t* alpha,
                                       slong nvars,
                                       slong Tbound,
                                       void* user_ctx,
                                       nmod_t mod)
{
    determinant_probe_ctx_t* ctx = (determinant_probe_ctx_t*) user_ctx;
    slong seq_len = 2 * Tbound;
    mp_limb_t* entry_scaled = (mp_limb_t*) malloc((size_t) ctx->matrix_size * (size_t) ctx->matrix_size * (size_t) nvars * sizeof(mp_limb_t));
    mp_limb_t* scaled_derivs = (mp_limb_t*) malloc((size_t) nvars * sizeof(mp_limb_t));
    ulong* exp_scratch = (ulong*) malloc((size_t) nvars * sizeof(ulong));
    geometric_poly_cache_t* caches = (geometric_poly_cache_t*) malloc((size_t) ctx->matrix_size * (size_t) ctx->matrix_size * sizeof(geometric_poly_cache_t));
    clock_t total_start = clock();
    if (entry_scaled == NULL || scaled_derivs == NULL || exp_scratch == NULL || caches == NULL) {
        free(entry_scaled);
        free(scaled_derivs);
        free(exp_scratch);
        free(caches);
        return 0;
    }

    nmod_mat_t eval_mat, inv_mat, adj_mat;
    nmod_mat_init(eval_mat, ctx->matrix_size, ctx->matrix_size, mod.n);
    nmod_mat_init(inv_mat, ctx->matrix_size, ctx->matrix_size, mod.n);
    nmod_mat_init(adj_mat, ctx->matrix_size, ctx->matrix_size, mod.n);

    clock_t prepare_start = clock();
    for (slong r = 0; r < ctx->matrix_size; r++) {
        for (slong c = 0; c < ctx->matrix_size; c++) {
            size_t idx = (size_t) r * (size_t) ctx->matrix_size + (size_t) c;
            if (!geometric_poly_cache_init(caches + idx,
                                           ctx->A->entries[r] + c,
                                           alpha,
                                           nvars,
                                           mod,
                                           exp_scratch,
                                           ctx->mctx)) {
                for (size_t k = 0; k <= idx; k++) {
                    geometric_poly_cache_clear(caches + k);
                }
                nmod_mat_clear(eval_mat);
                nmod_mat_clear(inv_mat);
                nmod_mat_clear(adj_mat);
                free(entry_scaled);
                free(scaled_derivs);
                free(exp_scratch);
                free(caches);
                return 0;
            }
            if (ctx->timing != NULL) {
                ctx->timing->entry_term_count += caches[idx].len;
            }
        }
    }
    if (ctx->timing != NULL) {
        ctx->timing->prepare += timer_seconds_since(prepare_start);
    }

    for (slong j = 0; j < seq_len; j++) {
        clock_t eval_start = clock();
        for (slong r = 0; r < ctx->matrix_size; r++) {
            for (slong c = 0; c < ctx->matrix_size; c++) {
                size_t base = ((size_t) r * (size_t) ctx->matrix_size + (size_t) c) * (size_t) nvars;
                geometric_poly_cache_accumulate_and_advance(&nmod_mat_entry(eval_mat, r, c),
                                                            scaled_derivs,
                                                            caches + ((size_t) r * (size_t) ctx->matrix_size + (size_t) c),
                                                            nvars,
                                                            mod);
                for (slong i = 0; i < nvars; i++) {
                    entry_scaled[base + (size_t) i] = scaled_derivs[i];
                }
            }
        }
        if (ctx->timing != NULL) {
            ctx->timing->eval_entries += timer_seconds_since(eval_start);
            ctx->timing->point_count++;
        }

        clock_t det_start = clock();
        a[j] = nmod_mat_det(eval_mat);
        if (ctx->timing != NULL) {
            ctx->timing->det += timer_seconds_since(det_start);
        }

        int have_inverse = 0;
        if (ctx->matrix_size > 0 && a[j] != 0) {
            clock_t inv_start = clock();
            have_inverse = nmod_mat_inv(inv_mat, eval_mat);
            if (ctx->timing != NULL) {
                ctx->timing->inv += timer_seconds_since(inv_start);
                if (have_inverse) {
                    ctx->timing->inverse_count++;
                }
            }
        }

        if (have_inverse) {
            clock_t trace_start = clock();
            for (slong var = 0; var < nvars; var++) {
                h[var * seq_len + j] = 0;
            }
            for (slong r = 0; r < ctx->matrix_size; r++) {
                for (slong c = 0; c < ctx->matrix_size; c++) {
                    mp_limb_t weight = nmod_mat_entry(inv_mat, c, r);
                    if (weight == 0) continue;
                    size_t base = ((size_t) r * (size_t) ctx->matrix_size + (size_t) c) * (size_t) nvars;
                    for (slong var = 0; var < nvars; var++) {
                        h[var * seq_len + j] = nmod_add(h[var * seq_len + j],
                                                        nmod_mul(weight, entry_scaled[base + (size_t) var], mod),
                                                        mod);
                    }
                }
            }
            for (slong var = 0; var < nvars; var++) {
                h[var * seq_len + j] = nmod_mul(a[j], h[var * seq_len + j], mod);
            }
            if (ctx->timing != NULL) {
                ctx->timing->direct_trace += timer_seconds_since(trace_start);
            }
        } else {
            clock_t adj_start = clock();
            nmod_mat_adjugate_classical(adj_mat, eval_mat, mod);
            for (slong var = 0; var < nvars; var++) {
                h[var * seq_len + j] = 0;
            }
            for (slong r = 0; r < ctx->matrix_size; r++) {
                for (slong c = 0; c < ctx->matrix_size; c++) {
                    mp_limb_t weight = nmod_mat_entry(adj_mat, c, r);
                    if (weight == 0) continue;
                    size_t base = ((size_t) r * (size_t) ctx->matrix_size + (size_t) c) * (size_t) nvars;
                    for (slong var = 0; var < nvars; var++) {
                        h[var * seq_len + j] = nmod_add(h[var * seq_len + j],
                                                        nmod_mul(weight, entry_scaled[base + (size_t) var], mod),
                                                        mod);
                    }
                }
            }
            if (ctx->timing != NULL) {
                ctx->timing->adjugate_trace += timer_seconds_since(adj_start);
                ctx->timing->singular_count++;
            }
        }
    }

    nmod_mat_clear(eval_mat);
    nmod_mat_clear(inv_mat);
    nmod_mat_clear(adj_mat);
    for (slong r = 0; r < ctx->matrix_size; r++) {
        for (slong c = 0; c < ctx->matrix_size; c++) {
            geometric_poly_cache_clear(caches + ((size_t) r * (size_t) ctx->matrix_size + (size_t) c));
        }
    }
    free(entry_scaled);
    free(scaled_derivs);
    free(exp_scratch);
    free(caches);
    if (ctx->timing != NULL) {
        ctx->timing->total += timer_seconds_since(total_start);
    }
    return 1;
}

static void sparse_probe_attach_timing(sparse_probe_fn_t probe_fn, void* probe_ctx, sparse_probe_timing_t* timing)
{
    if (probe_ctx == NULL) {
        return;
    }
    if (probe_fn == explicit_poly_probe_sequences) {
        ((explicit_poly_probe_ctx_t*) probe_ctx)->timing = timing;
    } else if (probe_fn == determinant_probe_sequences) {
        ((determinant_probe_ctx_t*) probe_ctx)->timing = timing;
    }
}

static slong estimate_term_bound_by_products(const poly_mat_t* A, slong cap, int by_rows, const nmod_mpoly_ctx_t mctx)
{
    slong outer = by_rows ? A->rows : A->cols;
    slong inner = by_rows ? A->cols : A->rows;
    slong bound = 1;

    for (slong i = 0; i < outer; i++) {
        slong sum = 0;
        for (slong j = 0; j < inner; j++) {
            const nmod_mpoly_struct* entry = by_rows ? (A->entries[i] + j) : (A->entries[j] + i);
            sum = sat_add_slong(sum, nmod_mpoly_length(entry, mctx), cap);
        }
        bound = sat_mul_slong(bound, sum, cap);
    }

    return bound;
}

static slong estimate_term_bound_permanent(const poly_mat_t* A, slong cap, const nmod_mpoly_ctx_t mctx)
{
    slong size = A->rows;
    if (size == 0) return 1;
    if (size > 20) {
        slong row_bound = estimate_term_bound_by_products(A, cap, 1, mctx);
        slong col_bound = estimate_term_bound_by_products(A, cap, 0, mctx);
        return row_bound < col_bound ? row_bound : col_bound;
    }

    size_t subsets = ((size_t) 1) << size;
    slong* dp = (slong*) calloc(subsets, sizeof(slong));
    slong* next = (slong*) calloc(subsets, sizeof(slong));
    if (dp == NULL || next == NULL) {
        free(dp);
        free(next);
        slong row_bound = estimate_term_bound_by_products(A, cap, 1, mctx);
        slong col_bound = estimate_term_bound_by_products(A, cap, 0, mctx);
        return row_bound < col_bound ? row_bound : col_bound;
    }

    dp[0] = 1;
    for (slong row = 0; row < size; row++) {
        memset(next, 0, subsets * sizeof(slong));
        for (size_t mask = 0; mask < subsets; mask++) {
            if (dp[mask] == 0) continue;
            for (slong col = 0; col < size; col++) {
                size_t bit = ((size_t) 1) << col;
                if ((mask & bit) == 0) {
                    slong len = nmod_mpoly_length(A->entries[row] + col, mctx);
                    slong contrib = sat_mul_slong(dp[mask], len, cap);
                    next[mask | bit] = sat_add_slong(next[mask | bit], contrib, cap);
                }
            }
        }
        memcpy(dp, next, subsets * sizeof(slong));
    }

    slong result = dp[subsets - 1];
    free(dp);
    free(next);
    return result;
}

static slong estimate_term_bound_degree_box(const slong* max_degs, slong nvars, slong cap)
{
    slong bound = 1;
    for (slong i = 0; i < nvars; i++) {
        bound = sat_mul_slong(bound, max_degs[i] + 1, cap);
    }
    return bound;
}

static void estimate_matrix_det_bounds(slong* max_degs,
                                       slong* term_bound,
                                       const poly_mat_t* A,
                                       slong nvars,
                                       const nmod_mpoly_ctx_t mctx)
{
    slong* row_max = (slong*) calloc((size_t) A->rows * nvars, sizeof(slong));
    slong* col_max = (slong*) calloc((size_t) A->cols * nvars, sizeof(slong));
    slong* entry_degs = (slong*) malloc((size_t) nvars * sizeof(slong));

    for (slong i = 0; i < nvars; i++) {
        max_degs[i] = 0;
    }

    if (row_max != NULL && col_max != NULL && entry_degs != NULL) {
        for (slong r = 0; r < A->rows; r++) {
            for (slong c = 0; c < A->cols; c++) {
                poly_max_degrees_per_var(entry_degs, A->entries[r] + c, mctx);
                for (slong v = 0; v < nvars; v++) {
                    slong deg = entry_degs[v];
                    if (deg > row_max[r * nvars + v]) row_max[r * nvars + v] = deg;
                    if (deg > col_max[c * nvars + v]) col_max[c * nvars + v] = deg;
                }
            }
        }

        for (slong v = 0; v < nvars; v++) {
            slong row_bound = 0;
            slong col_bound = 0;
            for (slong r = 0; r < A->rows; r++) row_bound += row_max[r * nvars + v];
            for (slong c = 0; c < A->cols; c++) col_bound += col_max[c * nvars + v];
            max_degs[v] = row_bound < col_bound ? row_bound : col_bound;
        }
    }

    free(row_max);
    free(col_max);
    free(entry_degs);

    slong permanent_bound = estimate_term_bound_permanent(A, SPARSE_ESTIMATE_CAP, mctx);
    slong degree_box_bound = estimate_term_bound_degree_box(max_degs, nvars, SPARSE_ESTIMATE_CAP);
    *term_bound = permanent_bound < degree_box_bound ? permanent_bound : degree_box_bound;
}

static int SparseInterpolateFromProbes(nmod_mpoly_t result,
                                       slong nvars,
                                       slong Tbound,
                                       int scout_mode,
                                       nmod_t mod,
                                       const slong* max_degs_per_var,
                                       sparse_probe_fn_t probe_fn,
                                       void* probe_ctx,
                                       slong attempt_limit,
                                       sparse_recovery_info_t* info,
                                       const nmod_mpoly_ctx_t mctx)
{
    nmod_mpoly_zero(result, mctx);
    sparse_recovery_info_reset(info);
    sparse_probe_attach_timing(probe_fn, probe_ctx, info != NULL ? &info->probe_timing : NULL);
    clock_t total_start = clock();
    if (Tbound <= 0) {
        if (info != NULL) {
            info->budget = Tbound;
            info->last_reconstructed_terms = 0;
            sparse_recovery_info_set_reason(info, "zero polynomial");
        }
        return 1;
    }

    slong seq_len = 2 * Tbound;
    mp_limb_t* alpha = (mp_limb_t*) malloc((size_t) nvars * sizeof(mp_limb_t));
    mp_limb_t* a = (mp_limb_t*) malloc((size_t) seq_len * sizeof(mp_limb_t));
    mp_limb_t* h = (mp_limb_t*) malloc((size_t) nvars * seq_len * sizeof(mp_limb_t));
    mp_limb_t* roots = (mp_limb_t*) malloc((size_t) Tbound * sizeof(mp_limb_t));
    mp_limb_t* denom_inv = (mp_limb_t*) malloc((size_t) Tbound * sizeof(mp_limb_t));
    mp_limb_t* coeffs = (mp_limb_t*) malloc((size_t) Tbound * sizeof(mp_limb_t));
    mp_limb_t* ce = (mp_limb_t*) malloc((size_t) nvars * Tbound * sizeof(mp_limb_t));
    ulong* recovered_exp = (ulong*) malloc((size_t) nvars * sizeof(ulong));
    sparse_multipoint_tree_t eval_tree;
    int success = 0;

    sparse_multipoint_tree_init_empty(&eval_tree);

    if (info != NULL) {
        info->budget = Tbound;
        info->attempt_limit = attempt_limit;
        info->workspace_bytes = sparse_workspace_bytes(nvars, Tbound);
    }

    if (!alpha || !a || !h || !roots || !denom_inv || !coeffs || !ce || !recovered_exp) {
        sparse_recovery_info_set_reason(info,
                                        "workspace allocation failed for T=%ld (~%.2f MiB)",
                                        Tbound,
                                        info != NULL ? ((double) info->workspace_bytes) / (1024.0 * 1024.0) : 0.0);
        goto cleanup;
    }

    for (slong attempt = 0; attempt < attempt_limit && !success; attempt++) {
        if (info != NULL) {
            info->attempts = attempt + 1;
        }
        if (sparse_debug_level() >= 2) {
            printf("  [Sparse probe] budget=%ld attempt=%ld/%ld\n", Tbound, attempt + 1, attempt_limit);
        }
        for (slong i = 0; i < nvars; i++) {
            alpha[i] = nmod_rand_nonzero(mod);
        }

        clock_t probe_start = clock();
        if (!probe_fn(a, h, alpha, nvars, Tbound, probe_ctx, mod)) {
            sparse_recovery_info_set_reason(info, "probe callback failed on attempt %ld", attempt + 1);
            continue;
        }
        if (info != NULL) {
            info->interpolation_timing.probe += timer_seconds_since(probe_start);
        }

        nmod_poly_t relation, locator;
        nmod_poly_init(relation, mod.n);
        nmod_poly_init(locator, mod.n);

        clock_t bm_start = clock();
        BM_timed(relation, a, Tbound, mod,
                 info != NULL ? &info->interpolation_timing : NULL);
        if (info != NULL) {
            info->interpolation_timing.bm += timer_seconds_since(bm_start);
        }
        slong t = nmod_poly_degree(relation);
        if (t < 0) t = 0;
        if (info != NULL) {
            info->last_relation_degree = t;
        }

        if (t == 0) {
            success = 1;
            sparse_recovery_info_set_reason(info, "recovered zero polynomial");
            nmod_poly_clear(relation);
            nmod_poly_clear(locator);
            break;
        }
        if (t > Tbound) {
            if (info != NULL) {
                info->likely_underbudget = 1;
            }
            sparse_recovery_info_set_reason(info,
                                            "BM returned relation degree %ld > budget %ld",
                                            t,
                                            Tbound);
            nmod_poly_clear(relation);
            nmod_poly_clear(locator);
            continue;
        }

        if (scout_mode && t == Tbound) {
            if (info != NULL) {
                info->likely_underbudget = 1;
            }
            sparse_recovery_info_set_reason(info,
                                            "BM saturated scout budget %ld before root recovery",
                                            Tbound);
            nmod_poly_clear(relation);
            nmod_poly_clear(locator);
            break;
        }

        build_locator_from_relation(locator, relation, t);
        clock_t root_start = clock();
        if (!get_simple_roots_strict(roots, locator, t, mod,
                                     info != NULL ? &info->interpolation_timing : NULL)) {
            if (info != NULL && t == Tbound) {
                info->likely_underbudget = 1;
            }
            sparse_recovery_info_set_reason(info,
                                            "locator did not split into %ld distinct simple roots",
                                            t);
            if (info != NULL) {
                info->interpolation_timing.root_factor += timer_seconds_since(root_start);
            }
            nmod_poly_clear(relation);
            nmod_poly_clear(locator);
            continue;
        }
        if (info != NULL) {
            info->interpolation_timing.root_factor += timer_seconds_since(root_start);
        }

        if (t >= 512) {
            clock_t tree_start = clock();
            if (!sparse_multipoint_tree_init(&eval_tree, roots, t, mod)) {
                sparse_recovery_info_set_reason(info,
                                                "failed to allocate fixed-roots multipoint tree for %ld roots",
                                                t);
                nmod_poly_clear(relation);
                nmod_poly_clear(locator);
                continue;
            }
            if (info != NULL) {
                info->interpolation_timing.vinvert_tree_build += timer_seconds_since(tree_start);
            }
        }

        clock_t vinvert_start = clock();
        clock_t denom_start = clock();
        if (!precompute_locator_inverse_denominators(denom_inv, locator, roots, t, mod)) {
            sparse_recovery_info_set_reason(info,
                                            "failed to precompute locator derivative inverses for %ld roots",
                                            t);
            sparse_multipoint_tree_clear(&eval_tree);
            nmod_poly_clear(relation);
            nmod_poly_clear(locator);
            continue;
        }
        if (info != NULL) {
            info->interpolation_timing.vinvert_denom += timer_seconds_since(denom_start);
        }
        Vinvert_preconditioned(coeffs, locator, roots, &eval_tree, denom_inv, a, t, mod,
                              info != NULL ? &info->interpolation_timing : NULL);
        for (slong i = 0; i < nvars; i++) {
            Vinvert_preconditioned(ce + i * Tbound, locator, roots, &eval_tree, denom_inv, h + i * seq_len, t, mod,
                                  info != NULL ? &info->interpolation_timing : NULL);
        }
        if (info != NULL) {
            info->interpolation_timing.vinvert += timer_seconds_since(vinvert_start);
        }
        sparse_multipoint_tree_clear(&eval_tree);

        nmod_mpoly_zero(result, mctx);
        int bad = 0;
        clock_t reconstruct_start = clock();
        for (slong k = 0; k < t && !bad; k++) {
            mp_limb_t coeff = coeffs[k];
            if (coeff == 0) {
                bad = 1;
                sparse_recovery_info_set_reason(info,
                                                "recovered zero coefficient at term %ld/%ld",
                                                k + 1,
                                                t);
                break;
            }

            mp_limb_t coeff_inv = nmod_inv(coeff, mod);
            for (slong i = 0; i < nvars; i++) {
                mp_limb_t exp_mod = nmod_mul(ce[i * Tbound + k], coeff_inv, mod);
                if (max_degs_per_var != NULL && (slong) exp_mod > max_degs_per_var[i]) {
                    bad = 1;
                    sparse_recovery_info_set_reason(info,
                                                    "recovered exponent %lu exceeds bound %ld for variable %ld",
                                                    (ulong) exp_mod,
                                                    max_degs_per_var[i],
                                                    i);
                    break;
                }
                recovered_exp[i] = (ulong) exp_mod;
            }

            if (!bad) {
                nmod_mpoly_push_term_ui_ui(result, coeff, recovered_exp, mctx);
            }
        }

        nmod_mpoly_combine_like_terms(result, mctx);
        if (info != NULL) {
            info->last_reconstructed_terms = nmod_mpoly_length(result, mctx);
            info->interpolation_timing.reconstruct += timer_seconds_since(reconstruct_start);
        }
        if (!bad && nmod_mpoly_length(result, mctx) == t) {
            success = 1;
            sparse_recovery_info_set_reason(info,
                                            "success with %ld recovered terms",
                                            t);
        } else if (!bad) {
            if (info != NULL && t == Tbound) {
                info->likely_underbudget = 1;
            }
            sparse_recovery_info_set_reason(info,
                                            "reconstructed %ld terms after combine_like_terms, expected %ld",
                                            nmod_mpoly_length(result, mctx),
                                            t);
        }

        if (!success) {
            nmod_mpoly_zero(result, mctx);
        }

        nmod_poly_clear(relation);
        nmod_poly_clear(locator);
    }

cleanup:
    free(alpha);
    free(a);
    free(h);
    free(roots);
    free(denom_inv);
    sparse_multipoint_tree_clear(&eval_tree);
    free(coeffs);
    free(ce);
    free(recovered_exp);

    if (!success) {
        if (info != NULL && info->last_reason[0] == '\0') {
            sparse_recovery_info_set_reason(info, "probe recovery failed after %ld attempts", info->attempts);
        }
        nmod_mpoly_zero(result, mctx);
    }
    if (info != NULL) {
        info->interpolation_timing.total += timer_seconds_since(total_start);
    }
    return success;
}

static int SparseInterpolateScoutBudgetReuse(nmod_mpoly_t result,
                                             slong nvars,
                                             slong Tbound,
                                             nmod_t mod,
                                             const slong* max_degs_per_var,
                                             sparse_probe_fn_t probe_fn,
                                             void* probe_ctx,
                                             sparse_scout_reuse_t* reuse,
                                             sparse_recovery_info_t* info,
                                             const nmod_mpoly_ctx_t mctx)
{
    nmod_mpoly_zero(result, mctx);
    sparse_recovery_info_reset(info);
    sparse_probe_attach_timing(probe_fn, probe_ctx, info != NULL ? &info->probe_timing : NULL);
    clock_t total_start = clock();
    slong seq_len = 2 * Tbound;
    mp_limb_t* a = (mp_limb_t*) malloc((size_t) seq_len * sizeof(mp_limb_t));
    mp_limb_t* h = (mp_limb_t*) malloc((size_t) nvars * seq_len * sizeof(mp_limb_t));
    mp_limb_t* roots = (mp_limb_t*) malloc((size_t) Tbound * sizeof(mp_limb_t));
    mp_limb_t* denom_inv = (mp_limb_t*) malloc((size_t) Tbound * sizeof(mp_limb_t));
    mp_limb_t* coeffs = (mp_limb_t*) malloc((size_t) Tbound * sizeof(mp_limb_t));
    mp_limb_t* ce = (mp_limb_t*) malloc((size_t) nvars * Tbound * sizeof(mp_limb_t));
    ulong* recovered_exp = (ulong*) malloc((size_t) nvars * sizeof(ulong));
    sparse_multipoint_tree_t eval_tree;
    int success = 0;
    int fresh_alpha = 0;

    sparse_multipoint_tree_init_empty(&eval_tree);

    if (info != NULL) {
        info->budget = Tbound;
        info->attempt_limit = 1;
        info->attempts = 1;
        info->workspace_bytes = sparse_workspace_bytes(nvars, Tbound);
    }

    if (!sparse_scout_reuse_prepare(reuse, nvars, Tbound) ||
        a == NULL || h == NULL || roots == NULL || denom_inv == NULL ||
        coeffs == NULL || ce == NULL || recovered_exp == NULL) {
        sparse_recovery_info_set_reason(info,
                                        "workspace allocation failed for scout reuse T=%ld (~%.2f MiB)",
                                        Tbound,
                                        info != NULL ? ((double) info->workspace_bytes) / (1024.0 * 1024.0) : 0.0);
        goto cleanup;
    }

    if (!reuse->valid) {
        for (slong i = 0; i < nvars; i++) {
            reuse->alpha[i] = nmod_rand_nonzero(mod);
        }
        sparse_bm_state_reset_sequence(&reuse->bm_state);
        reuse->valid = 1;
        fresh_alpha = 1;
    }

    if (sparse_debug_level() >= 2) {
        printf("  [Sparse probe] budget=%ld attempt=1/1%s\n",
               Tbound, fresh_alpha ? " (fresh scout alpha)" : " (reused scout alpha)");
    }

    clock_t probe_start = clock();
    if (!probe_fn(a, h, reuse->alpha, nvars, Tbound, probe_ctx, mod)) {
        sparse_recovery_info_set_reason(info, "probe callback failed in scout reuse");
        goto cleanup;
    }
    if (info != NULL) {
        info->interpolation_timing.probe += timer_seconds_since(probe_start);
    }

    nmod_poly_t relation, locator;
    nmod_poly_init(relation, mod.n);
    nmod_poly_init(locator, mod.n);

    clock_t bm_start = clock();
    BM_timed_incremental(relation, a, Tbound, mod, &reuse->bm_state,
                         info != NULL ? &info->interpolation_timing : NULL);
    if (info != NULL) {
        info->interpolation_timing.bm += timer_seconds_since(bm_start);
    }
    slong t = nmod_poly_degree(relation);
    if (t < 0) t = 0;
    if (info != NULL) {
        info->last_relation_degree = t;
    }

    if (t == 0) {
        success = 1;
        sparse_recovery_info_set_reason(info, "recovered zero polynomial");
        nmod_poly_clear(relation);
        nmod_poly_clear(locator);
        goto cleanup;
    }

    if (t > Tbound || t == Tbound) {
        if (info != NULL) {
            info->likely_underbudget = 1;
        }
        if (t > Tbound) {
            sparse_recovery_info_set_reason(info,
                                            "BM returned relation degree %ld > budget %ld",
                                            t, Tbound);
        } else {
            sparse_recovery_info_set_reason(info,
                                            "BM saturated scout budget %ld before root recovery",
                                            Tbound);
        }
        nmod_poly_clear(relation);
        nmod_poly_clear(locator);
        goto cleanup;
    }

    build_locator_from_relation(locator, relation, t);
    clock_t root_start = clock();
    if (!get_simple_roots_strict(roots, locator, t, mod,
                                 info != NULL ? &info->interpolation_timing : NULL)) {
        sparse_recovery_info_set_reason(info,
                                        "locator did not split into %ld distinct simple roots",
                                        t);
        if (info != NULL) {
            info->interpolation_timing.root_factor += timer_seconds_since(root_start);
        }
        nmod_poly_clear(relation);
        nmod_poly_clear(locator);
        goto cleanup;
    }
    if (info != NULL) {
        info->interpolation_timing.root_factor += timer_seconds_since(root_start);
    }

    if (t >= 512) {
        clock_t tree_start = clock();
        if (!sparse_multipoint_tree_init(&eval_tree, roots, t, mod)) {
            sparse_recovery_info_set_reason(info,
                                            "failed to allocate fixed-roots multipoint tree for %ld roots",
                                            t);
            nmod_poly_clear(relation);
            nmod_poly_clear(locator);
            goto cleanup;
        }
        if (info != NULL) {
            info->interpolation_timing.vinvert_tree_build += timer_seconds_since(tree_start);
        }
    }

    clock_t vinvert_start = clock();
    clock_t denom_start = clock();
    if (!precompute_locator_inverse_denominators(denom_inv, locator, roots, t, mod)) {
        sparse_recovery_info_set_reason(info,
                                        "failed to precompute locator derivative inverses for %ld roots",
                                        t);
        nmod_poly_clear(relation);
        nmod_poly_clear(locator);
        goto cleanup;
    }
    if (info != NULL) {
        info->interpolation_timing.vinvert_denom += timer_seconds_since(denom_start);
    }
    Vinvert_preconditioned(coeffs, locator, roots, &eval_tree, denom_inv, a, t, mod,
                          info != NULL ? &info->interpolation_timing : NULL);
    for (slong i = 0; i < nvars; i++) {
        Vinvert_preconditioned(ce + i * Tbound, locator, roots, &eval_tree, denom_inv, h + i * seq_len, t, mod,
                              info != NULL ? &info->interpolation_timing : NULL);
    }
    if (info != NULL) {
        info->interpolation_timing.vinvert += timer_seconds_since(vinvert_start);
    }
    sparse_multipoint_tree_clear(&eval_tree);

    nmod_mpoly_zero(result, mctx);
    int bad = 0;
    clock_t reconstruct_start = clock();
    for (slong k = 0; k < t && !bad; k++) {
        mp_limb_t coeff = coeffs[k];
        if (coeff == 0) {
            bad = 1;
            sparse_recovery_info_set_reason(info,
                                            "recovered zero coefficient at term %ld/%ld",
                                            k + 1, t);
            break;
        }
        mp_limb_t coeff_inv = nmod_inv(coeff, mod);
        for (slong i = 0; i < nvars; i++) {
            mp_limb_t exp_mod = nmod_mul(ce[i * Tbound + k], coeff_inv, mod);
            if (max_degs_per_var != NULL && (slong) exp_mod > max_degs_per_var[i]) {
                bad = 1;
                sparse_recovery_info_set_reason(info,
                                                "recovered exponent %lu exceeds bound %ld for variable %ld",
                                                (ulong) exp_mod, max_degs_per_var[i], i);
                break;
            }
            recovered_exp[i] = (ulong) exp_mod;
        }
        if (!bad) {
            nmod_mpoly_push_term_ui_ui(result, coeff, recovered_exp, mctx);
        }
    }
    nmod_mpoly_combine_like_terms(result, mctx);
    if (info != NULL) {
        info->last_reconstructed_terms = nmod_mpoly_length(result, mctx);
        info->interpolation_timing.reconstruct += timer_seconds_since(reconstruct_start);
    }
    if (!bad && nmod_mpoly_length(result, mctx) == t) {
        success = 1;
        sparse_recovery_info_set_reason(info, "success with %ld recovered terms", t);
    } else if (!bad) {
        sparse_recovery_info_set_reason(info,
                                        "reconstructed %ld terms after combine_like_terms, expected %ld",
                                        nmod_mpoly_length(result, mctx), t);
    }
    if (!success) {
        nmod_mpoly_zero(result, mctx);
    }

    nmod_poly_clear(relation);
    nmod_poly_clear(locator);

cleanup:
    free(a);
    free(h);
    free(roots);
    free(denom_inv);
    sparse_multipoint_tree_clear(&eval_tree);
    free(coeffs);
    free(ce);
    free(recovered_exp);

    if (!success && info != NULL && info->last_reason[0] == '\0') {
        sparse_recovery_info_set_reason(info, "scout reuse attempt failed");
    }
    if (!success && info != NULL) {
        nmod_mpoly_zero(result, mctx);
    }
    if (info != NULL) {
        info->interpolation_timing.total += timer_seconds_since(total_start);
    }
    return success;
}

static int SparseInterpolateFromProbesAdaptive(nmod_mpoly_t result,
                                               slong nvars,
                                               slong Tupper,
                                               nmod_t mod,
                                               const slong* max_degs_per_var,
                                               sparse_probe_fn_t probe_fn,
                                               void* probe_ctx,
                                               sparse_recovery_info_t* info,
                                               const nmod_mpoly_ctx_t mctx)
{
    sparse_recovery_info_t last_info;
    sparse_recovery_info_reset(&last_info);
    sparse_scout_reuse_t scout_reuse;
    sparse_scout_reuse_init_empty(&scout_reuse);
    sparse_probe_timing_t total_probe_timing;
    sparse_interpolation_timing_t total_interp_timing;
    memset(&total_probe_timing, 0, sizeof(total_probe_timing));
    memset(&total_interp_timing, 0, sizeof(total_interp_timing));
    slong total_attempts = 0;
    slong budgets_tried = 0;

    if (Tupper <= 0) {
        nmod_mpoly_zero(result, mctx);
        if (info != NULL) {
            *info = last_info;
            info->budget = 0;
            sparse_recovery_info_set_reason(info, "zero polynomial");
        }
        return 1;
    }

    slong budget = sparse_get_start_budget(Tupper, max_degs_per_var, nvars);
    if (sparse_debug_level() >= 1 && max_degs_per_var != NULL) {
        slong active_vars = 0;
        slong degree_sum = 0;
        for (slong i = 0; i < nvars; i++) {
            if (max_degs_per_var[i] > 0) active_vars++;
            degree_sum += max_degs_per_var[i];
        }
        printf("  [Sparse adaptive] start budget %ld chosen from degree profile (active_vars=%ld, degree_sum=%ld, upper=%ld)\n",
               budget, active_vars, degree_sum, Tupper);
    }
    while (1) {
        sparse_recovery_info_t current_info;
        slong attempt_limit = sparse_get_scout_attempt_limit(budget, Tupper);
        if (sparse_debug_level() >= 1) {
            printf("  [Sparse adaptive] trying budget %ld / upper bound %ld with %ld attempts\n",
                   budget, Tupper, attempt_limit);
        }
        int success;
        if (budget < Tupper) {
            if (sparse_debug_level() >= 1 && scout_reuse.valid) {
                printf("    [Sparse adaptive] reusing scout alpha/BM state from previous saturated budget.\n");
            }
            success = SparseInterpolateScoutBudgetReuse(result,
                                                        nvars,
                                                        budget,
                                                        mod,
                                                        max_degs_per_var,
                                                        probe_fn,
                                                        probe_ctx,
                                                        &scout_reuse,
                                                        &current_info,
                                                        mctx);
            if (!success && !current_info.likely_underbudget) {
                if (sparse_debug_level() >= 1) {
                    printf("    [Sparse adaptive] scout reuse path failed without clear underbudget evidence; retrying with regular attempts.\n");
                }
                sparse_scout_reuse_clear(&scout_reuse);
                success = SparseInterpolateFromProbes(result,
                                                      nvars,
                                                      budget,
                                                      1,
                                                      mod,
                                                      max_degs_per_var,
                                                      probe_fn,
                                                      probe_ctx,
                                                      attempt_limit,
                                                      &current_info,
                                                      mctx);
            }
            if (!(current_info.likely_underbudget && current_info.last_relation_degree >= budget)) {
                sparse_scout_reuse_clear(&scout_reuse);
            }
        } else {
            sparse_scout_reuse_clear(&scout_reuse);
            success = SparseInterpolateFromProbes(result,
                                                  nvars,
                                                  budget,
                                                  0,
                                                  mod,
                                                  max_degs_per_var,
                                                  probe_fn,
                                                  probe_ctx,
                                                  attempt_limit,
                                                  &current_info,
                                                  mctx);
        }
        last_info = current_info;
        budgets_tried++;
        total_attempts += current_info.attempts;
        sparse_probe_timing_add(&total_probe_timing, &current_info.probe_timing);
        sparse_interpolation_timing_add(&total_interp_timing, &current_info.interpolation_timing);
        if (sparse_debug_level() >= 1) {
            printf("    [Sparse budget result] attempts=%ld/%ld relation_deg=%ld terms=%ld underbudget=%s reason=%s\n",
                   current_info.attempts,
                   current_info.attempt_limit,
                   current_info.last_relation_degree,
                   current_info.last_reconstructed_terms,
                   current_info.likely_underbudget ? "yes" : "no",
                   current_info.last_reason[0] ? current_info.last_reason : "n/a");
            printf("    [Sparse budget timing] probe=%.6f BM=%.6f roots=%.6f Vinvert=%.6f reconstruct=%.6f total=%.6f\n",
                   current_info.interpolation_timing.probe,
                   current_info.interpolation_timing.bm,
                   current_info.interpolation_timing.root_factor,
                   current_info.interpolation_timing.vinvert,
                   current_info.interpolation_timing.reconstruct,
                   current_info.interpolation_timing.total);
            if (current_info.interpolation_timing.bm_discrepancy > 0.0 ||
                current_info.interpolation_timing.bm_update > 0.0 ||
                current_info.interpolation_timing.bm_finalize > 0.0) {
                printf("      [Sparse BM timing] discrepancy=%.6f update=%.6f finalize=%.6f\n",
                       current_info.interpolation_timing.bm_discrepancy,
                       current_info.interpolation_timing.bm_update,
                       current_info.interpolation_timing.bm_finalize);
            }
            if (current_info.interpolation_timing.bm_update_prepare > 0.0 ||
                current_info.interpolation_timing.bm_update_axpy > 0.0 ||
                current_info.interpolation_timing.bm_update_trim > 0.0 ||
                current_info.interpolation_timing.bm_update_promote > 0.0) {
                printf("        [Sparse BM update timing] prepare=%.6f axpy=%.6f trim=%.6f promote=%.6f\n",
                       current_info.interpolation_timing.bm_update_prepare,
                       current_info.interpolation_timing.bm_update_axpy,
                       current_info.interpolation_timing.bm_update_trim,
                       current_info.interpolation_timing.bm_update_promote);
            }
            if (current_info.interpolation_timing.bm_update_count > 0) {
                sparse_print_bm_histogram(&current_info.interpolation_timing, "        ");
            }
            printf("      [Sparse root timing] squarefree=%.6f frobenius=%.6f split=%.6f sort=%.6f fallback_factor=%.6f\n",
                   current_info.interpolation_timing.root_squarefree,
                   current_info.interpolation_timing.root_frobenius,
                   current_info.interpolation_timing.root_split,
                   current_info.interpolation_timing.root_sort,
                   current_info.interpolation_timing.root_fallback_factor);
            printf("      [Sparse Vinvert timing] denom=%.6f tree_build=%.6f build=%.6f mul=%.6f extract=%.6f eval=%.6f\n",
                   current_info.interpolation_timing.vinvert_denom,
                   current_info.interpolation_timing.vinvert_tree_build,
                   current_info.interpolation_timing.vinvert_build,
                   current_info.interpolation_timing.vinvert_mul,
                   current_info.interpolation_timing.vinvert_extract,
                   current_info.interpolation_timing.vinvert_eval);
        }
        if (success) {
            if (info != NULL) {
                *info = current_info;
                info->attempts = total_attempts;
                info->probe_timing = total_probe_timing;
                info->interpolation_timing = total_interp_timing;
            }
            if (sparse_debug_level() >= 1) {
                printf("  [Sparse adaptive summary] budgets=%ld total_attempts=%ld success_budget=%ld\n",
                       budgets_tried, total_attempts, budget);
            }
            sparse_scout_reuse_clear(&scout_reuse);
            return 1;
        }

        if (budget >= Tupper) {
            break;
        }

        if (current_info.likely_underbudget && sparse_debug_level() >= 1) {
            printf("    [Sparse adaptive] escalating because budget %ld is likely below the true sparsity.\n", budget);
        }

        if (budget > Tupper / 2) {
            budget = Tupper;
        } else {
            budget *= 2;
        }
    }

    if (info != NULL) {
        *info = last_info;
        info->attempts = total_attempts;
        info->probe_timing = total_probe_timing;
        info->interpolation_timing = total_interp_timing;
    }
    sparse_scout_reuse_clear(&scout_reuse);
    nmod_mpoly_zero(result, mctx);
    return 0;
}

int DerivativeSparseInterpolatePolynomial(nmod_mpoly_t result,
                                         const nmod_mpoly_t f,
                                         slong n,
                                         nmod_t mod,
                                         const slong* max_degs_per_var,
                                         const nmod_mpoly_ctx_t mctx)
{
    slong t = nmod_mpoly_length(f, mctx);
    slong* local_max_degs = NULL;
    const slong* bounds = max_degs_per_var;
    explicit_poly_probe_ctx_t probe_ctx;

    if (bounds == NULL) {
        local_max_degs = (slong*) malloc((size_t) n * sizeof(slong));
        if (local_max_degs == NULL) {
            nmod_mpoly_zero(result, mctx);
            return 0;
        }
        poly_max_degrees_per_var(local_max_degs, f, mctx);
        bounds = local_max_degs;
    }

    for (slong i = 0; i < n; i++) {
        if ((mp_limb_t) bounds[i] >= mod.n) {
            if (TIMING) {
                printf("DerivativeSparseInterpolatePolynomial failed: partial degree bound reaches characteristic.\n");
            }
            free(local_max_degs);
            nmod_mpoly_zero(result, mctx);
            return 0;
        }
    }

    if (!explicit_poly_probe_ctx_init(&probe_ctx, f, n, mctx)) {
        free(local_max_degs);
        nmod_mpoly_zero(result, mctx);
        return 0;
    }

    sparse_recovery_info_t recovery_info;
    sparse_probe_attach_timing(explicit_poly_probe_sequences, &probe_ctx, &recovery_info.probe_timing);
    int success = SparseInterpolateFromProbes(result, n, t, 0, mod, bounds,
                                              explicit_poly_probe_sequences,
                                              &probe_ctx,
                                              sparse_get_attempt_limit(t),
                                              &recovery_info,
                                              mctx);

    if (TIMING && sparse_debug_level() >= 1) {
        printf("  [Sparse polynomial timing] probe=%.6f BM=%.6f roots=%.6f Vinvert=%.6f reconstruct=%.6f total=%.6f\n",
               recovery_info.interpolation_timing.probe,
               recovery_info.interpolation_timing.bm,
               recovery_info.interpolation_timing.root_factor,
               recovery_info.interpolation_timing.vinvert,
               recovery_info.interpolation_timing.reconstruct,
               recovery_info.interpolation_timing.total);
        if (recovery_info.interpolation_timing.bm_discrepancy > 0.0 ||
            recovery_info.interpolation_timing.bm_update > 0.0 ||
            recovery_info.interpolation_timing.bm_finalize > 0.0) {
            printf("    [Sparse BM timing] discrepancy=%.6f update=%.6f finalize=%.6f\n",
                   recovery_info.interpolation_timing.bm_discrepancy,
                   recovery_info.interpolation_timing.bm_update,
                   recovery_info.interpolation_timing.bm_finalize);
        }
        if (recovery_info.interpolation_timing.bm_update_prepare > 0.0 ||
            recovery_info.interpolation_timing.bm_update_axpy > 0.0 ||
            recovery_info.interpolation_timing.bm_update_trim > 0.0 ||
            recovery_info.interpolation_timing.bm_update_promote > 0.0) {
            printf("      [Sparse BM update timing] prepare=%.6f axpy=%.6f trim=%.6f promote=%.6f\n",
                   recovery_info.interpolation_timing.bm_update_prepare,
                   recovery_info.interpolation_timing.bm_update_axpy,
                   recovery_info.interpolation_timing.bm_update_trim,
                   recovery_info.interpolation_timing.bm_update_promote);
        }
        if (recovery_info.interpolation_timing.bm_update_count > 0) {
            sparse_print_bm_histogram(&recovery_info.interpolation_timing, "      ");
        }
        printf("    [Sparse root timing] squarefree=%.6f frobenius=%.6f split=%.6f sort=%.6f fallback_factor=%.6f\n",
               recovery_info.interpolation_timing.root_squarefree,
               recovery_info.interpolation_timing.root_frobenius,
               recovery_info.interpolation_timing.root_split,
               recovery_info.interpolation_timing.root_sort,
               recovery_info.interpolation_timing.root_fallback_factor);
        printf("    [Sparse Vinvert timing] denom=%.6f tree_build=%.6f build=%.6f mul=%.6f extract=%.6f eval=%.6f\n",
               recovery_info.interpolation_timing.vinvert_denom,
               recovery_info.interpolation_timing.vinvert_tree_build,
               recovery_info.interpolation_timing.vinvert_build,
               recovery_info.interpolation_timing.vinvert_mul,
               recovery_info.interpolation_timing.vinvert_extract,
               recovery_info.interpolation_timing.vinvert_eval);
    }

    if (!success && TIMING) {
        printf("DerivativeSparseInterpolatePolynomial failed: %s\n",
               recovery_info.last_reason[0] ? recovery_info.last_reason : "unknown reason");
    }

    explicit_poly_probe_ctx_clear(&probe_ctx);
    free(local_max_degs);
    return success;
}

void ComputePolyMatrixDet(nmod_mpoly_t det_poly,
                         const poly_mat_t* A,
                         slong nvars,
                         mp_limb_t p,
                         const nmod_mpoly_ctx_t mctx)
{
    my_timer_t timer;
    ensure_global_state();
    nmod_mpoly_zero(det_poly, mctx);

    if (A->rows != A->cols) {
        return;
    }

    nmod_t mod;
    nmod_init(&mod, p);

    slong* max_degs = (slong*) malloc((size_t) nvars * sizeof(slong));
    sparse_recovery_info_t recovery_info;
    sparse_recovery_info_reset(&recovery_info);
    int fallback_direct = 0;
    slong Tbound = 0;
    slong user_hard_limit = sparse_get_env_slong("DIXON_SPARSE_T_HARD_LIMIT", 0, 0);
    char fallback_reason[512];
    fallback_reason[0] = '\0';

    if (max_degs == NULL) {
        fallback_direct = 1;
        snprintf(fallback_reason,
                 sizeof(fallback_reason),
                 "failed to allocate max degree bounds");
    } else {
        estimate_matrix_det_bounds(max_degs, &Tbound, A, nvars, mctx);
        for (slong i = 0; i < nvars; i++) {
            if ((mp_limb_t) max_degs[i] >= mod.n) {
                fallback_direct = 1;
                snprintf(fallback_reason,
                         sizeof(fallback_reason),
                         "theoretical restriction: partial degree bound for variable %ld is %ld >= characteristic %lu",
                         i,
                         max_degs[i],
                         mod.n);
                break;
            }
        }
        if (!fallback_direct && user_hard_limit > 0 && Tbound > user_hard_limit) {
            fallback_direct = 1;
            snprintf(fallback_reason,
                     sizeof(fallback_reason),
                     "user hard limit DIXON_SPARSE_T_HARD_LIMIT=%ld is below estimated T upper bound %ld",
                     user_hard_limit,
                     Tbound);
        }
    }

    if (TIMING && sparse_debug_level() >= 1) {
        printf("Sparse.pdf probe bounds: T <= %ld\n", Tbound);
        printf("Sparse.pdf partial degree bounds: ");
        for (slong i = 0; i < nvars; i++) {
            printf("%ld ", max_degs ? max_degs[i] : -1L);
        }
        printf("\n");
        if (user_hard_limit > 0) {
            printf("Sparse.pdf hard T limit from environment: %ld\n", user_hard_limit);
        } else {
            printf("Sparse.pdf hard T limit from environment: disabled\n");
        }
    }

    if (Tbound == 0 && !fallback_direct) {
        free(max_degs);
        nmod_mpoly_zero(det_poly, mctx);
        return;
    }

    if (!fallback_direct) {
        determinant_probe_ctx_t probe_ctx;

        if (!determinant_probe_ctx_init(&probe_ctx, A, nvars, mctx)) {
            fallback_direct = 1;
            snprintf(fallback_reason,
                     sizeof(fallback_reason),
                     "failed to initialize determinant derivative probe context");
        } else {
            if (TIMING) timer_start(&timer);
            int success = SparseInterpolateFromProbesAdaptive(det_poly, nvars, Tbound, mod, max_degs,
                                                              determinant_probe_sequences,
                                                              &probe_ctx,
                                                              &recovery_info,
                                                              mctx);
            if (TIMING && sparse_debug_level() >= 1) {
                timer_stop(&timer);
                timer_print(&timer, "    Sparse.pdf determinant probing");
                printf("    [Sparse determinant timing] prepare=%.6f eval=%.6f det=%.6f inv=%.6f direct-trace=%.6f adjugate=%.6f points=%ld terms=%ld invertible=%ld singular=%ld\n",
                       recovery_info.probe_timing.prepare,
                       recovery_info.probe_timing.eval_entries,
                       recovery_info.probe_timing.det,
                       recovery_info.probe_timing.inv,
                       recovery_info.probe_timing.direct_trace,
                       recovery_info.probe_timing.adjugate_trace,
                       recovery_info.probe_timing.point_count,
                       recovery_info.probe_timing.entry_term_count,
                       recovery_info.probe_timing.inverse_count,
                       recovery_info.probe_timing.singular_count);
                printf("    [Sparse interpolation timing] probe=%.6f BM=%.6f roots=%.6f Vinvert=%.6f reconstruct=%.6f total=%.6f\n",
                       recovery_info.interpolation_timing.probe,
                       recovery_info.interpolation_timing.bm,
                       recovery_info.interpolation_timing.root_factor,
                       recovery_info.interpolation_timing.vinvert,
                       recovery_info.interpolation_timing.reconstruct,
                       recovery_info.interpolation_timing.total);
                if (recovery_info.interpolation_timing.bm_discrepancy > 0.0 ||
                    recovery_info.interpolation_timing.bm_update > 0.0 ||
                    recovery_info.interpolation_timing.bm_finalize > 0.0) {
                    printf("      [Sparse BM timing] discrepancy=%.6f update=%.6f finalize=%.6f\n",
                           recovery_info.interpolation_timing.bm_discrepancy,
                           recovery_info.interpolation_timing.bm_update,
                           recovery_info.interpolation_timing.bm_finalize);
                }
                if (recovery_info.interpolation_timing.bm_update_prepare > 0.0 ||
                    recovery_info.interpolation_timing.bm_update_axpy > 0.0 ||
                    recovery_info.interpolation_timing.bm_update_trim > 0.0 ||
                    recovery_info.interpolation_timing.bm_update_promote > 0.0) {
                    printf("        [Sparse BM update timing] prepare=%.6f axpy=%.6f trim=%.6f promote=%.6f\n",
                           recovery_info.interpolation_timing.bm_update_prepare,
                           recovery_info.interpolation_timing.bm_update_axpy,
                           recovery_info.interpolation_timing.bm_update_trim,
                           recovery_info.interpolation_timing.bm_update_promote);
                }
                if (recovery_info.interpolation_timing.bm_update_count > 0) {
                    sparse_print_bm_histogram(&recovery_info.interpolation_timing, "        ");
                }
                printf("      [Sparse root timing] squarefree=%.6f frobenius=%.6f split=%.6f sort=%.6f fallback_factor=%.6f\n",
                       recovery_info.interpolation_timing.root_squarefree,
                       recovery_info.interpolation_timing.root_frobenius,
                       recovery_info.interpolation_timing.root_split,
                       recovery_info.interpolation_timing.root_sort,
                       recovery_info.interpolation_timing.root_fallback_factor);
                printf("      [Sparse Vinvert timing] denom=%.6f tree_build=%.6f build=%.6f mul=%.6f extract=%.6f eval=%.6f\n",
                       recovery_info.interpolation_timing.vinvert_denom,
                       recovery_info.interpolation_timing.vinvert_tree_build,
                       recovery_info.interpolation_timing.vinvert_build,
                       recovery_info.interpolation_timing.vinvert_mul,
                       recovery_info.interpolation_timing.vinvert_extract,
                       recovery_info.interpolation_timing.vinvert_eval);
            }
            if (!success) {
                fallback_direct = 1;
                snprintf(fallback_reason,
                         sizeof(fallback_reason),
                         "probe recovery failed after budget %ld, %ld/%ld attempts, workspace ~%.2f MiB; last reason: %s",
                         recovery_info.budget,
                         recovery_info.attempts,
                         recovery_info.attempt_limit,
                         ((double) recovery_info.workspace_bytes) / (1024.0 * 1024.0),
                         recovery_info.last_reason[0] ? recovery_info.last_reason : "unknown");
            }
            determinant_probe_ctx_clear(&probe_ctx);
        }
    }

    if (fallback_direct) {
        if (TIMING) {
            printf("WARNING: falling back to direct determinant.\n");
            if (fallback_reason[0] != '\0') {
                printf("         reason: %s\n", fallback_reason);
            }
            if (recovery_info.used) {
                printf("         last probe budget: %ld, attempts: %ld/%ld, workspace: %.2f MiB\n",
                       recovery_info.budget,
                       recovery_info.attempts,
                       recovery_info.attempt_limit,
                       ((double) recovery_info.workspace_bytes) / (1024.0 * 1024.0));
                if (recovery_info.last_reason[0] != '\0') {
                    printf("         last probe detail: %s\n", recovery_info.last_reason);
                }
            }
            timer_start(&timer);
        }
        poly_mat_det(det_poly, A, mctx);
        if (TIMING) {
            timer_stop(&timer);
            timer_print(&timer, "    Direct determinant fallback");
        }
    }

    free(max_degs);
}

void myrandpoly(nmod_mpoly_t f, slong n, slong T, slong D,
                nmod_t mod, const nmod_mpoly_ctx_t mctx)
{
    ensure_global_state();
    nmod_mpoly_zero(f, mctx);

    ulong** used_exps = (ulong**) malloc((size_t) T * sizeof(ulong*));
    slong num_terms = 0;

    while (num_terms < T) {
        ulong* exp = (ulong*) calloc((size_t) n, sizeof(ulong));
        for (slong j = 0; j < n; j++) {
            exp[j] = n_randint(global_state, D + 1);
        }

        int duplicate = 0;
        for (slong k = 0; k < num_terms; k++) {
            int same = 1;
            for (slong j = 0; j < n; j++) {
                if (exp[j] != used_exps[k][j]) {
                    same = 0;
                    break;
                }
            }
            if (same) {
                duplicate = 1;
                break;
            }
        }

        if (!duplicate) {
            mp_limb_t coeff;
            do {
                coeff = n_randint(global_state, mod.n);
            } while (coeff == 0);
            nmod_mpoly_push_term_ui_ui(f, coeff, exp, mctx);
            used_exps[num_terms++] = exp;
        } else {
            free(exp);
        }
    }

    for (slong i = 0; i < num_terms; i++) {
        free(used_exps[i]);
    }
    free(used_exps);

    nmod_mpoly_combine_like_terms(f, mctx);
}

void test_random_polynomial(void)
{
    printf("\n========== Testing Sparse.pdf Derivative Interpolation ==========\n");

    slong n = 3;
    slong T = 8000;
    slong d = 1000;
    mp_limb_t p = n_nextprime(2 * (n + 2) * T * T * d, 1);

    printf("Parameters: n=%ld, T=%ld, d=%ld\n", n, T, d);
    printf("Using prime: %lu\n", p);

    nmod_t mod;
    nmod_init(&mod, p);

    char** vars = (char**) malloc((size_t) n * sizeof(char*));
    for (slong i = 0; i < n; i++) {
        vars[i] = (char*) malloc(10);
        sprintf(vars[i], "x%ld", i);
    }

    nmod_mpoly_ctx_t mctx;
    nmod_mpoly_ctx_init(mctx, n, ORD_LEX, p);

    nmod_mpoly_t f, g;
    nmod_mpoly_init(f, mctx);
    nmod_mpoly_init(g, mctx);

    myrandpoly(f, n, T, d, mod, mctx);

    slong* max_degs = (slong*) malloc((size_t) n * sizeof(slong));
    poly_max_degrees_per_var(max_degs, f, mctx);

    if (!DerivativeSparseInterpolatePolynomial(g, f, n, mod, max_degs, mctx)) {
        printf("Derivative sparse interpolation failed.\n");
    }

    nmod_mpoly_sort_terms(f, mctx);
    nmod_mpoly_sort_terms(g, mctx);

    printf("Recovered polynomial with %ld terms\n", nmod_mpoly_length(g, mctx));
    if (nmod_mpoly_equal(f, g, mctx)) {
        printf("Success: Polynomials match!\n");
    } else {
        printf("Error: Polynomials do not match!\n");
    }

    free(max_degs);
    nmod_mpoly_clear(f, mctx);
    nmod_mpoly_clear(g, mctx);
    nmod_mpoly_ctx_clear(mctx);

    for (slong i = 0; i < n; i++) {
        free(vars[i]);
    }
    free(vars);
}

int huang_test(void)
{
    printf("Code version: probe-based Sparse.pdf interpolation without old method\n");
    ensure_global_state();

    printf("========== Testing Polynomial Matrix Determinant ==========\n");

    mp_limb_t p = 65537;
    slong n = 3;
    slong k = 2;

    char** vars = (char**) malloc((size_t) n * sizeof(char*));
    for (slong i = 0; i < n; i++) {
        vars[i] = (char*) malloc(10);
        sprintf(vars[i], "x%ld", i);
    }

    nmod_mpoly_ctx_t mctx;
    nmod_mpoly_ctx_init(mctx, n, ORD_LEX, p);

    poly_mat_t A;
    poly_mat_init(&A, k, k, mctx);

    nmod_mpoly_t entry;
    nmod_mpoly_init(entry, mctx);
    ulong* exp = (ulong*) calloc((size_t) n, sizeof(ulong));

    nmod_mpoly_zero(entry, mctx);
    exp[0] = 0; exp[1] = 1; exp[2] = 0;
    nmod_mpoly_push_term_ui_ui(entry, 3, exp, mctx);
    exp[0] = 3; exp[1] = 1; exp[2] = 0;
    nmod_mpoly_push_term_ui_ui(entry, 22, exp, mctx);
    poly_mat_entry_set(&A, 0, 0, entry, mctx);

    nmod_mpoly_zero(entry, mctx);
    exp[0] = 1; exp[1] = 2; exp[2] = 0;
    nmod_mpoly_push_term_ui_ui(entry, 64, exp, mctx);
    poly_mat_entry_set(&A, 0, 1, entry, mctx);

    nmod_mpoly_zero(entry, mctx);
    exp[0] = 2; exp[1] = 0; exp[2] = 0;
    nmod_mpoly_push_term_ui_ui(entry, 61, exp, mctx);
    exp[0] = 2; exp[1] = 1; exp[2] = 0;
    nmod_mpoly_push_term_ui_ui(entry, 91, exp, mctx);
    poly_mat_entry_set(&A, 1, 0, entry, mctx);

    nmod_mpoly_zero(entry, mctx);
    exp[0] = 1; exp[1] = 3; exp[2] = 0;
    nmod_mpoly_push_term_ui_ui(entry, 87, exp, mctx);
    exp[0] = 0; exp[1] = 2; exp[2] = 4;
    nmod_mpoly_push_term_ui_ui(entry, 26, exp, mctx);
    exp[0] = 0; exp[1] = 3; exp[2] = 2;
    nmod_mpoly_push_term_ui_ui(entry, 89, exp, mctx);
    poly_mat_entry_set(&A, 1, 1, entry, mctx);

    free(exp);
    nmod_mpoly_clear(entry, mctx);

    nmod_mpoly_t det_poly, actual_det;
    nmod_mpoly_init(det_poly, mctx);
    nmod_mpoly_init(actual_det, mctx);

    poly_mat_det(actual_det, &A, mctx);
    ComputePolyMatrixDet(det_poly, &A, n, p, mctx);

    nmod_mpoly_sort_terms(det_poly, mctx);
    nmod_mpoly_sort_terms(actual_det, mctx);

    printf("\n========== FINAL RESULTS ==========\n");
    printf("Computed determinant:\n");
    nmod_mpoly_print_pretty(det_poly, (const char**) vars, mctx);
    printf("\n\nActual determinant:\n");
    nmod_mpoly_print_pretty(actual_det, (const char**) vars, mctx);
    printf("\n");

    if (nmod_mpoly_equal(det_poly, actual_det, mctx)) {
        printf("\nSuccess: Determinants match!\n");
    } else {
        printf("\nError: Determinants do not match!\n");
    }

    nmod_mpoly_clear(det_poly, mctx);
    nmod_mpoly_clear(actual_det, mctx);
    poly_mat_clear(&A, mctx);
    nmod_mpoly_ctx_clear(mctx);

    for (slong i = 0; i < n; i++) {
        free(vars[i]);
    }
    free(vars);

    test_random_polynomial();

    clear_global_state_if_needed();
    flint_cleanup();
    return 0;
}
