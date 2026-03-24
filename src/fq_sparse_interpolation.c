/**
 * fq_sparse_interpolation.c - Implementation of sparse polynomial interpolation
 * 
 * This code is for testing only. Experiments show this method is much less 
 * efficient than usual for dense polynomials from the AO algorithm, and is 
 * applicable only in very sparse cases.
 */

#include "fq_sparse_interpolation.h"

/* Global random state */
flint_rand_t global_state;

/* Timer functions */
void timer_start(my_timer_t* t) {
    t->start = clock();
}

void timer_stop(my_timer_t* t) {
    t->elapsed = ((double)(clock() - t->start)) / CLOCKS_PER_SEC;
}

void timer_print(my_timer_t* t, const char* label) {
    if (TIMING) {
        printf("  [TIMING] %s: %.6f seconds\n", label, t->elapsed);
    }
}

/* Polynomial matrix functions */
void poly_mat_init(poly_mat_t* mat, slong rows, slong cols, const nmod_mpoly_ctx_t ctx) {
    mat->rows = rows;
    mat->cols = cols;
    mat->entries = (nmod_mpoly_struct**)malloc(rows * sizeof(nmod_mpoly_struct*));
    for (slong i = 0; i < rows; i++) {
        mat->entries[i] = (nmod_mpoly_struct*)malloc(cols * sizeof(nmod_mpoly_struct));
        for (slong j = 0; j < cols; j++) {
            nmod_mpoly_init(mat->entries[i] + j, ctx);
        }
    }
}

void poly_mat_clear(poly_mat_t* mat, const nmod_mpoly_ctx_t ctx) {
    for (slong i = 0; i < mat->rows; i++) {
        for (slong j = 0; j < mat->cols; j++) {
            nmod_mpoly_clear(mat->entries[i] + j, ctx);
        }
        free(mat->entries[i]);
    }
    free(mat->entries);
}

void poly_mat_entry_set(poly_mat_t* mat, slong i, slong j, const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx) {
    nmod_mpoly_set(mat->entries[i] + j, poly, ctx);
}

/* Compute 2x2 determinant for testing */
void poly_mat_det_2x2(nmod_mpoly_t det, const poly_mat_t* mat, 
                      const nmod_mpoly_ctx_t ctx) {
    if (mat->rows != 2 || mat->cols != 2) {
        nmod_mpoly_zero(det, ctx);
        return;
    }
    
    nmod_mpoly_t temp1, temp2;
    nmod_mpoly_init(temp1, ctx);
    nmod_mpoly_init(temp2, ctx);
    
    /* det = a00*a11 - a01*a10 */
    nmod_mpoly_mul(temp1, mat->entries[0] + 0, mat->entries[1] + 1, ctx);
    nmod_mpoly_mul(temp2, mat->entries[0] + 1, mat->entries[1] + 0, ctx);
    nmod_mpoly_sub(det, temp1, temp2, ctx);
    
    nmod_mpoly_clear(temp1, ctx);
    nmod_mpoly_clear(temp2, ctx);
}

/* General determinant computation using recursive expansion */
void poly_mat_det(nmod_mpoly_t det, const poly_mat_t* mat, 
                  const nmod_mpoly_ctx_t ctx) {

    if (mat->rows != mat->cols) {
        nmod_mpoly_zero(det, ctx);
        return;
    }

    slong n = mat->rows;
    
    if (n == 1) {
        nmod_mpoly_set(det, mat->entries[0] + 0, ctx);
        return;
    }
    
    if (n == 2) {
        poly_mat_det_2x2(det, mat, ctx);
        return;
    }

    /* For larger matrices, expand along first row */
    nmod_mpoly_zero(det, ctx);
    
    nmod_mpoly_t minor_det, temp;
    nmod_mpoly_init(minor_det, ctx);
    nmod_mpoly_init(temp, ctx);
  
    for (slong j = 0; j < n; j++) {
        /* Build minor matrix */
        poly_mat_t minor;
        poly_mat_init(&minor, n-1, n-1, ctx);
        
        for (slong i = 1; i < n; i++) {
            slong col_idx = 0;
            for (slong k = 0; k < n; k++) {
                if (k != j) {
                    nmod_mpoly_set(minor.entries[i-1] + col_idx, 
                                      mat->entries[i] + k, ctx);
                    col_idx++;
                }
            }
        }
       
        /* Compute determinant of minor */
        poly_mat_det(minor_det, &minor, ctx);
        
        /* Multiply by corresponding element */
        nmod_mpoly_mul(temp, mat->entries[0] + j, minor_det, ctx);
  
        /* Add with sign */
        if (j % 2 == 0) {
            nmod_mpoly_add(det, det, temp, ctx);
                         
        } else {
            nmod_mpoly_sub(det, det, temp, ctx);
        }
        poly_mat_clear(&minor, ctx);
    }
    fflush(stdout);
    nmod_mpoly_clear(minor_det, ctx);
    nmod_mpoly_clear(temp, ctx);
}

/* Get the maximum total degree of a polynomial */
slong poly_max_total_degree(const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx) {
    slong max_deg = 0;
    slong nvars = nmod_mpoly_ctx_nvars(ctx);
    slong len = nmod_mpoly_length(poly, ctx);
    
    for (slong i = 0; i < len; i++) {
        ulong* exp = (ulong*)malloc(nvars * sizeof(ulong));
        nmod_mpoly_get_term_exp_ui(exp, poly, i, ctx);
        
        slong total_deg = 0;
        for (slong j = 0; j < nvars; j++) {
            total_deg += exp[j];
        }
        
        if (total_deg > max_deg) {
            max_deg = total_deg;
        }
        
        free(exp);
    }
    
    return max_deg;
}

/* Get the maximum degree in each variable of a polynomial */
void poly_max_degrees_per_var(slong* max_degs, const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx) {
    slong nvars = nmod_mpoly_ctx_nvars(ctx);
    slong len = nmod_mpoly_length(poly, ctx);
    
    /* Initialize to 0 */
    for (slong i = 0; i < nvars; i++) {
        max_degs[i] = 0;
    }
    
    /* Find max degree in each variable */
    for (slong i = 0; i < len; i++) {
        ulong* exp = (ulong*)malloc(nvars * sizeof(ulong));
        nmod_mpoly_get_term_exp_ui(exp, poly, i, ctx);
        
        for (slong j = 0; j < nvars; j++) {
            if (exp[j] > max_degs[j]) {
                max_degs[j] = exp[j];
            }
        }
        
        free(exp);
    }
}

/* Berlekamp-Massey Algorithm (BM) */
void BM(nmod_poly_t C, const mp_limb_t* s, slong N, nmod_t mod) {
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
    
    for (slong n = 0; n < 2*N; n++) {
        /* Calculate d = s[n] + sum(C[i] * s[n-i]) for i=1 to L */
        d = s[n];
        for (slong i = 1; i <= L && i <= n; i++) {
            temp = nmod_poly_get_coeff_ui(C, i);
            temp2 = nmod_mul(temp, s[n - i], mod);
            d = nmod_add(d, temp2, mod);
        }
        
        if (d == 0) {
            k++;
        } else if (n < 2*L) {
            /* C = C - (d/b) * x^k * B */
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
            /* Save old C */
            nmod_poly_set(T, C);
            
            /* C = C - (d/b) * x^k * B */
            nmod_poly_t xkB;
            nmod_poly_init(xkB, mod.n);
            nmod_poly_shift_left(xkB, B, k);
            temp = nmod_inv(b, mod);
            temp = nmod_mul(d, temp, mod);
            nmod_poly_scalar_mul_nmod(xkB, xkB, temp);
            nmod_poly_sub(C, C, xkB);
            nmod_poly_clear(xkB);
            
            /* Update B, L, k, b */
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

/* Vinvert - compute coefficients from roots using partial fraction decomposition */
void Vinvert(mp_limb_t* c1, const nmod_poly_t c, const mp_limb_t* v, 
             const mp_limb_t* a, slong n, nmod_t mod) {
    
    my_timer_t timer;
    if (DETAILED_TIMING) timer_start(&timer);
    
    nmod_poly_t d, q, q1, q2;
    mp_limb_t temp;
    
    nmod_poly_init(d, mod.n);
    nmod_poly_init(q, mod.n);
    nmod_poly_init(q1, mod.n);
    nmod_poly_init(q2, mod.n);
    
    /* Build polynomial d = sum(a[i] * z^(n+1-i)) */
    nmod_poly_zero(d);
    for (slong i = 0; i < n; i++) {
        nmod_poly_set_coeff_ui(d, n - i, a[i]);
    }
    
    /* q = c * d */
    nmod_poly_mul(q, c, d);
    
    /* Extract coefficients for q1 */
    nmod_poly_zero(q1);
    for (slong i = 0; i < n; i++) {
        temp = nmod_poly_get_coeff_ui(q, 2*n - i);
        nmod_poly_set_coeff_ui(q1, n - 1 - i, temp);
    }
    
    /* q2 = derivative of c */
    nmod_poly_derivative(q2, c);
    
    /* Evaluate and compute coefficients */
    for (slong i = 0; i < n; i++) {
        mp_limb_t num, den, q1_val, q2_val;
        
        /* Evaluate q1 at v[i] */
        q1_val = nmod_poly_evaluate_nmod(q1, v[i]);
        
        /* Evaluate q2 at v[i] */
        q2_val = nmod_poly_evaluate_nmod(q2, v[i]);
        
        /* Compute num / den */
        if (q2_val != 0) {
            temp = nmod_inv(q2_val, mod);
            c1[i] = nmod_mul(q1_val, temp, mod);
        } else {
            /* Fallback: just use a[i] if derivative is zero */
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

/* Compare function for sorting by root value */
int compare_pairs_by_root(const void* a, const void* b) {
    const term_pair* pa = (const term_pair*)a;
    const term_pair* pb = (const term_pair*)b;
    if (pa->root < pb->root) return -1;
    if (pa->root > pb->root) return 1;
    return 0;
}

/* Use FLINT's built-in functions for root finding */
void find_roots_flint(nmod_poly_t poly, mp_limb_t* roots, 
                     slong* num_roots, nmod_t mod) {
    slong degree = nmod_poly_degree(poly);
    *num_roots = 0;
    
    if (degree == 0) return;
    
    if (degree == 1) {
        /* Linear polynomial: ax + b = 0 */
        mp_limb_t a, b;
        a = nmod_poly_get_coeff_ui(poly, 1);
        b = nmod_poly_get_coeff_ui(poly, 0);
        
        if (a != 0) {
            mp_limb_t a_inv = nmod_inv(a, mod);
            roots[0] = nmod_mul(b, a_inv, mod);
            roots[0] = nmod_neg(roots[0], mod);
            *num_roots = 1;
        }
        return;
    }
    
    /* For higher degree, use factorization */
    nmod_poly_factor_t fac;
    nmod_poly_factor_init(fac);
    nmod_poly_factor(fac, poly);
    
    for (slong i = 0; i < fac->num; i++) {
        if (nmod_poly_degree(fac->p + i) == 1) {
            mp_limb_t a, b;
            a = nmod_poly_get_coeff_ui(fac->p + i, 1);
            b = nmod_poly_get_coeff_ui(fac->p + i, 0);
            
            if (a != 0) {
                mp_limb_t a_inv = nmod_inv(a, mod);
                mp_limb_t root = nmod_mul(b, a_inv, mod);
                root = nmod_neg(root, mod);
                
                /* Add root with multiplicity */
                for (slong j = 0; j < fac->exp[i] && *num_roots < degree; j++) {
                    roots[*num_roots] = root;
                    (*num_roots)++;
                }
            }
        }
    }
    
    nmod_poly_factor_clear(fac);
}

/* MC - Monomials and Coefficients (Algorithm 3.1) */
term_list MC(const mp_limb_t* a, slong N, nmod_t mod) {
    my_timer_t timer, timer_inner;
    if (DETAILED_TIMING) timer_start(&timer);
    
    term_list result = {NULL, 0};
    nmod_poly_t f, Lambda;
    nmod_poly_factor_t fac;
    mp_limb_t* v;
    mp_limb_t* c;
    mp_limb_t* q;
    slong m;
    
    nmod_poly_init(f, mod.n);
    nmod_poly_init(Lambda, mod.n);
    nmod_poly_factor_init(fac);
    
    /* Find minimal polynomial using BM */
    BM(f, a, N, mod);
    m = nmod_poly_degree(f);
    
    if (m == 0) {
        nmod_poly_clear(f);
        nmod_poly_clear(Lambda);
        nmod_poly_factor_clear(fac);
        return result;
    }
    
    /* Build Lambda(z) = z^m + sum(f[i] * z^(m-i)) */
    nmod_poly_set_coeff_ui(Lambda, m, 1);
    for (slong i = 1; i <= m; i++) {
        mp_limb_t coeff = nmod_poly_get_coeff_ui(f, i);
        nmod_poly_set_coeff_ui(Lambda, m - i, coeff);
    }
    
    /* Allocate arrays */
    v = (mp_limb_t*)malloc(m * sizeof(mp_limb_t));
    c = (mp_limb_t*)malloc(m * sizeof(mp_limb_t));
    q = (mp_limb_t*)malloc(m * sizeof(mp_limb_t));
    for (slong i = 0; i < m; i++) {
        q[i] = a[i];  /* Get first m values */
    }
    
    /* Use FLINT's root finding */
    if (DETAILED_TIMING) timer_start(&timer_inner);
    
    slong root_count = 0;
    find_roots_flint(Lambda, v, &root_count, mod);
    
    if (DETAILED_TIMING) {
        timer_stop(&timer_inner);
        timer_print(&timer_inner, "    Root finding (FLINT)");
    }
    
    /* Verify roots in debug mode */
    if (DEBUG) {
        printf("Verifying roots...\n");
        for (slong i = 0; i < root_count; i++) {
            mp_limb_t eval = nmod_poly_evaluate_nmod(Lambda, v[i]);
            if (eval != 0) {
                printf("WARNING: Root %ld is not valid!\n", i);
            }
        }
    }
    
    /* If not enough roots found, use default values */
    while (root_count < m) {
        if (root_count > 0) {
            v[root_count] = v[0];  /* Duplicate first root */
        } else {
            v[root_count] = 1;  /* Use 1 as default */
        }
        root_count++;
    }
    
    /* Compute coefficients using Vinvert */
    if (m > 1) {
        Vinvert(c, Lambda, v, q, m, mod);
    } else {
        c[0] = a[0];
    }
    
    /* Build result */
    result.length = m;
    result.pairs = (term_pair*)malloc(m * sizeof(term_pair));
    
    for (slong i = 0; i < m; i++) {
        result.pairs[i].coeff = c[i];
        result.pairs[i].root = v[i];
    }
    
    /* Sort by root value */
    qsort(result.pairs, result.length, sizeof(term_pair), compare_pairs_by_root);
    
    /* Cleanup */
    free(v);
    free(c);
    free(q);
    
    nmod_poly_clear(f);
    nmod_poly_clear(Lambda);
    nmod_poly_factor_clear(fac);
    
    if (DETAILED_TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "  MC (total)");
    }
    
    return result;
}

/* Compute discrete logarithm with bounded search */
slong mylog(const mp_limb_t omega, const mp_limb_t B, slong d, nmod_t mod) {
    if (B == 1) {
        return 0;
    }
    
    /* Try all exponents from 0 to d */
    mp_limb_t omega_pow = 1;
    
    for (slong k = 0; k <= d; k++) {
        if (omega_pow == B) {
            return k;
        }
        omega_pow = nmod_mul(omega_pow, omega, mod);
    }
    
    /* If not found, print warning in debug mode */
    if (DEBUG) {
        printf("WARNING: mylog failed to find discrete log of %lu (searched up to %ld)\n", B, d);
    }
    
    return -1; /* Not found */
}

/* Myeval with detailed timing */
void Myeval(nmod_mat_t c, const nmod_mpoly_t f, slong n, 
            const mp_limb_t* alpha, nmod_t mod, const nmod_mpoly_ctx_t mctx) {
    
    my_timer_t timer_total, timer_step;
    if (DETAILED_TIMING) timer_start(&timer_total);
    
    slong t = nmod_mpoly_length(f, mctx);
    
    if (t == 0) {
        nmod_mat_zero(c);
        return;
    }
    
    /* Step 1: Get monomials and coefficients */
    if (DETAILED_TIMING) timer_start(&timer_step);
    mp_limb_t* coeffs = (mp_limb_t*)malloc(t * sizeof(mp_limb_t));
    ulong** exps = (ulong**)malloc(t * sizeof(ulong*));
    
    for (slong i = 0; i < t; i++) {
        exps[i] = (ulong*)malloc(n * sizeof(ulong));
        coeffs[i] = nmod_mpoly_get_term_coeff_ui(f, i, mctx);
        nmod_mpoly_get_term_exp_ui(exps[i], f, i, mctx);
    }
    if (DETAILED_TIMING) {
        timer_stop(&timer_step);
        timer_print(&timer_step, "      Get coeffs/exps");
    }
    
    /* Step 2: Evaluate monomials at alpha */
    if (DETAILED_TIMING) timer_start(&timer_step);
    nmod_mat_t V;
    nmod_mat_init(V, t, 1, mod.n);
    
    mp_limb_t u, temp;
    
    for (slong i = 0; i < t; i++) {
        u = 1;
        for (slong j = 0; j < n; j++) {
            if (exps[i][j] > 0) {
                temp = nmod_pow_ui(alpha[j], exps[i][j], mod);
                u = nmod_mul(u, temp, mod);
            }
        }
        nmod_mat_entry(V, i, 0) = u;
    }
    if (DETAILED_TIMING) {
        timer_stop(&timer_step);
        timer_print(&timer_step, "      Evaluate monomials");
    }
    
    /* Step 3: Build Vandermonde matrix */
    if (DETAILED_TIMING) timer_start(&timer_step);
    nmod_mat_t M;
    nmod_mat_init(M, t, t, mod.n);
    
    for (slong i = 0; i < t; i++) {
        u = 1;
        for (slong j = 0; j < t; j++) {
            nmod_mat_entry(M, j, i) = u;
            u = nmod_mul(u, nmod_mat_entry(V, i, 0), mod);
        }
    }
    if (DETAILED_TIMING) {
        timer_stop(&timer_step);
        timer_print(&timer_step, "      Build Vandermonde");
    }
    
    /* Step 4: Coefficient vector */
    if (DETAILED_TIMING) timer_start(&timer_step);
    nmod_mat_t C1;
    nmod_mat_init(C1, t, 1, mod.n);
    for (slong i = 0; i < t; i++) {
        nmod_mat_entry(C1, i, 0) = coeffs[i];
    }
    if (DETAILED_TIMING) {
        timer_stop(&timer_step);
        timer_print(&timer_step, "      Setup coeff vector");
    }
    
    /* Step 5: First matrix multiplication a = M * C1 */
    if (DETAILED_TIMING) timer_start(&timer_step);
    nmod_mat_t a;
    nmod_mat_init(a, t, 1, mod.n);
    nmod_mat_mul(a, M, C1);
    if (DETAILED_TIMING) {
        timer_stop(&timer_step);
        timer_print(&timer_step, "      First mat mul");
    }
    
    /* Step 6: Scale for second evaluation */
    if (DETAILED_TIMING) timer_start(&timer_step);
    for (slong i = 0; i < t; i++) {
        temp = nmod_mat_entry(V, i, 0);
        u = nmod_pow_ui(temp, t, mod);
        for (slong j = 0; j < t; j++) {
            temp = nmod_mat_entry(M, j, i);
            temp = nmod_mul(temp, u, mod);
            nmod_mat_entry(M, j, i) = temp;
        }
    }
    if (DETAILED_TIMING) {
        timer_stop(&timer_step);
        timer_print(&timer_step, "      Scale matrix");
    }
    
    /* Step 7: Second matrix multiplication b = M * C1 */
    if (DETAILED_TIMING) timer_start(&timer_step);
    nmod_mat_t b;
    nmod_mat_init(b, t, 1, mod.n);
    nmod_mat_mul(b, M, C1);
    if (DETAILED_TIMING) {
        timer_stop(&timer_step);
        timer_print(&timer_step, "      Second mat mul");
    }
    
    /* Step 8: Combine results */
    if (DETAILED_TIMING) timer_start(&timer_step);
    for (slong i = 0; i < t; i++) {
        nmod_mat_entry(c, 0, i) = nmod_mat_entry(a, i, 0);
        nmod_mat_entry(c, 0, i + t) = nmod_mat_entry(b, i, 0);
    }
    if (DETAILED_TIMING) {
        timer_stop(&timer_step);
        timer_print(&timer_step, "      Combine results");
    }
    
    /* Cleanup */
    for (slong i = 0; i < t; i++) {
        free(exps[i]);
    }
    free(coeffs);
    free(exps);
    
    nmod_mat_clear(V);
    nmod_mat_clear(M);
    nmod_mat_clear(C1);
    nmod_mat_clear(a);
    nmod_mat_clear(b);
    
    if (DETAILED_TIMING) {
        timer_stop(&timer_total);
        timer_print(&timer_total, "    Myeval total");
    }
}

/* TotalMyeval with detailed timing */
void TotalMyeval(nmod_mat_t M, const nmod_mpoly_t f, slong n,
                 const mp_limb_t* alpha, const mp_limb_t omega, 
                 nmod_t mod, const nmod_mpoly_ctx_t mctx) {
    
    my_timer_t timer, timer_step;
    if (TIMING) timer_start(&timer);
    
    slong T = nmod_mpoly_length(f, mctx);
    
    if (T == 0) {
        nmod_mat_zero(M);
        return;
    }
    
    printf("\n  TotalMyeval details (T=%ld, n=%ld):\n", T, n);
    
    /* First row - evaluate at alpha */
    if (DETAILED_TIMING) timer_start(&timer_step);
    nmod_mat_t row;
    nmod_mat_init(row, 1, 2*T, mod.n);
    
    printf("  Row 0 (base evaluation):\n");
    Myeval(row, f, n, alpha, mod, mctx);
    
    for (slong i = 0; i < 2*T; i++) {
        nmod_mat_entry(M, 0, i) = nmod_mat_entry(row, 0, i);
    }
    nmod_mat_clear(row);
    if (DETAILED_TIMING) {
        timer_stop(&timer_step);
        timer_print(&timer_step, "    First row total");
    }
    
    /* Remaining rows - modify one variable at a time */
    if (DETAILED_TIMING) timer_start(&timer_step);
    mp_limb_t* beta = (mp_limb_t*)malloc(n * sizeof(mp_limb_t));
    if (DETAILED_TIMING) {
        timer_stop(&timer_step);
        timer_print(&timer_step, "    Beta allocation");
    }
    
    for (slong k = 0; k < n; k++) {
        my_timer_t timer_row;
        if (DETAILED_TIMING) timer_start(&timer_row);
        
        printf("  Row %ld (modify var %ld):\n", k+1, k);
        
        /* Copy alpha to beta */
        if (DETAILED_TIMING) timer_start(&timer_step);
        for (slong i = 0; i < n; i++) {
            beta[i] = alpha[i];
        }
        
        /* Modify k-th variable: beta[k] = alpha[k] * omega */
        beta[k] = nmod_mul(alpha[k], omega, mod);
        if (DETAILED_TIMING) {
            timer_stop(&timer_step);
            timer_print(&timer_step, "    Setup beta");
        }
        
        /* Evaluate */
        nmod_mat_init(row, 1, 2*T, mod.n);
        Myeval(row, f, n, beta, mod, mctx);
        
        if (DETAILED_TIMING) timer_start(&timer_step);
        for (slong i = 0; i < 2*T; i++) {
            nmod_mat_entry(M, k + 1, i) = nmod_mat_entry(row, 0, i);
        }
        nmod_mat_clear(row);
        if (DETAILED_TIMING) {
            timer_stop(&timer_step);
            timer_print(&timer_step, "    Copy to result");
        }
        
        if (DETAILED_TIMING) {
            timer_stop(&timer_row);
            timer_print(&timer_row, "    Row total");
        }
    }
    
    /* Cleanup */
    free(beta);
    
    if (TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "TotalMyeval");
    }
}

/* Diversification transformation */
void Mydiver(nmod_mpoly_t g, const nmod_mpoly_t f, slong n,
             const mp_limb_t* zeta, nmod_t mod, 
             const nmod_mpoly_ctx_t mctx) {
    
    my_timer_t timer;
    if (TIMING) timer_start(&timer);
    
    slong t = nmod_mpoly_length(f, mctx);
    
    if (t == 0) {
        nmod_mpoly_zero(g, mctx);
        return;
    }
    
    nmod_mpoly_zero(g, mctx);
    
    mp_limb_t coeff, u, temp;
    
    /* Process each term */
    for (slong i = 0; i < t; i++) {
        /* Get coefficient and exponents */
        coeff = nmod_mpoly_get_term_coeff_ui(f, i, mctx);
        ulong* exp = (ulong*)malloc(n * sizeof(ulong));
        nmod_mpoly_get_term_exp_ui(exp, f, i, mctx);
        
        /* Multiply coefficient by zeta powers */
        u = 1;
        for (slong j = 0; j < n; j++) {
            if (exp[j] > 0) {
                temp = nmod_pow_ui(zeta[j], exp[j], mod);
                u = nmod_mul(u, temp, mod);
            }
        }
        
        /* Add transformed term */
        temp = nmod_mul(u, coeff, mod);
        nmod_mpoly_push_term_ui_ui(g, temp, exp, mctx);
        
        free(exp);
    }
    
    if (TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "Mydiver");
    }
}

/* Hash function for mp_limb_t */
static inline ulong limb_hash(const mp_limb_t x, ulong capacity) {
    return (x * 2654435761UL) % capacity;
}

void omega_hash_init(omega_hash_table* ht, slong max_power) {
    ht->capacity = max_power * 2 + 100;
    ht->size = 0;
    ht->keys = (mp_limb_t*)malloc(ht->capacity * sizeof(mp_limb_t));
    ht->values = (slong*)malloc(ht->capacity * sizeof(slong));
    
    for (slong i = 0; i < ht->capacity; i++) {
        ht->keys[i] = 0;
        ht->values[i] = -1;
    }
}

void omega_hash_clear(omega_hash_table* ht) {
    free(ht->keys);
    free(ht->values);
}

void omega_hash_insert(omega_hash_table* ht, const mp_limb_t key, slong value) {
    ulong idx = limb_hash(key, ht->capacity);
    
    while (ht->values[idx] != -1) {
        if (ht->keys[idx] == key) {
            return;
        }
        idx = (idx + 1) % ht->capacity;
    }
    
    ht->keys[idx] = key;
    ht->values[idx] = value;
    ht->size++;
}

slong omega_hash_find(omega_hash_table* ht, const mp_limb_t key) {
    ulong idx = limb_hash(key, ht->capacity);
    
    while (ht->values[idx] != -1) {
        if (ht->keys[idx] == key) {
            return ht->values[idx];
        }
        idx = (idx + 1) % ht->capacity;
        
        static slong max_probes = 0;
        if (max_probes == 0) {
            max_probes = ht->capacity;
        }
        slong probes = 0;
        if (++probes > max_probes) {
            break;
        }
    }
    
    return -1;
}

/* MBOT function */
void MBOT(nmod_mpoly_t result, const nmod_mat_t M, slong n, slong T,
          const mp_limb_t omega, const mp_limb_t* zeta, nmod_t mod,
          slong* max_degs_per_var, const nmod_mpoly_ctx_t mctx) {
    
    my_timer_t timer, timer_total;
    timer_start(&timer_total);
    
    printf("\n=== MBOT Performance Analysis ===\n");
    printf("Parameters: n=%ld, T=%ld\n", n, T);
    
    /* Process first row to get reference terms */
    if (TIMING) timer_start(&timer);
    mp_limb_t* first_row = (mp_limb_t*)malloc(2 * T * sizeof(mp_limb_t));
    for (slong j = 0; j < 2 * T; j++) {
        first_row[j] = nmod_mat_entry(M, 0, j);
    }
    if (TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "Extract first row");
    }
    
    /* MC for first row */
    printf("\nProcessing first row:\n");
    term_list a = MC(first_row, T, mod);
    
    if (a.length == 0) {
        nmod_mpoly_zero(result, mctx);
        goto cleanup_first;
    }
    
    slong t = a.length;
    
    if (t != T) {
        printf("WARNING: MC found %ld terms, expected %ld\n", t, T);
    }
    
    /* Extract coefficients and roots from first row */
    mp_limb_t* c = (mp_limb_t*)malloc(t * sizeof(mp_limb_t));
    mp_limb_t* first_roots = (mp_limb_t*)malloc(t * sizeof(mp_limb_t));
    for (slong i = 0; i < t; i++) {
        c[i] = a.pairs[i].coeff;
        first_roots[i] = a.pairs[i].root;
    }
    
    /* Precompute omega powers */
    if (TIMING) timer_start(&timer);
    slong max_power = 0;
    for (slong i = 0; i < n; i++) {
        if (max_degs_per_var[i] > max_power) {
            max_power = max_degs_per_var[i];
        }
    }
    max_power = max_power * 2 + 100;
    
    if (max_power > 1000000) {
        printf("ERROR: Max power %ld is too large, limiting to 1000000\n", max_power);
        max_power = 1000000;
    }
    
    printf("\nPrecomputing omega powers up to %ld\n", max_power);
    mp_limb_t* omega_powers = (mp_limb_t*)malloc((max_power + 1) * sizeof(mp_limb_t));
    if (omega_powers == NULL) {
        printf("ERROR: Failed to allocate memory for omega powers\n");
        nmod_mpoly_zero(result, mctx);
        goto cleanup_first;
    }
    
    for (slong k = 0; k <= max_power; k++) {
        omega_powers[k] = nmod_pow_ui(omega, k, mod);
    }
    if (TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "Precompute omega powers");
    }
    
    /* Build hash table */
    printf("\nBuilding omega power hash table...\n");
    if (TIMING) timer_start(&timer);
    
    omega_hash_table omega_table;
    omega_hash_init(&omega_table, max_power);
    
    for (slong k = 0; k <= max_power; k++) {
        omega_hash_insert(&omega_table, omega_powers[k], k);
    }
    
    if (TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "Build omega hash table");
    }
    
    /* Initialize root matrix N */
    nmod_mat_t N;
    nmod_mat_init(N, n + 1, t, mod.n);
    
    /* Set first row with roots */
    for (slong j = 0; j < t; j++) {
        nmod_mat_entry(N, 0, j) = first_roots[j];
    }
    
    /* Process each subsequent row */
    printf("\nProcessing remaining %ld rows:\n", n);
    if (TIMING) timer_start(&timer);
    
    int all_matched = 1;
    for (slong i = 1; i <= n; i++) {
        my_timer_t timer_row;
        if (DETAILED_TIMING) timer_start(&timer_row);
        
        printf("  Row %ld:\n", i);
        
        /* Extract current row */
        mp_limb_t* current_row = (mp_limb_t*)malloc(2 * T * sizeof(mp_limb_t));
        for (slong j = 0; j < 2 * T; j++) {
            current_row[j] = nmod_mat_entry(M, i, j);
        }
        
        term_list current = MC(current_row, T, mod);
        
        if (current.length != t) {
            printf("    WARNING: Row %ld has %ld terms, expected %ld\n", i, current.length, t);
            all_matched = 0;
        }
        
        /* Root matching section */
        my_timer_t timer_match;
        if (DETAILED_TIMING) timer_start(&timer_match);
        
        /* Create a permutation array to match roots */
        slong* perm = (slong*)malloc(t * sizeof(slong));
        int* used = (int*)calloc(current.length, sizeof(int));
        
        /* Initialize permutation to -1 (not matched) */
        for (slong j = 0; j < t; j++) {
            perm[j] = -1;
        }
        
        /* Try to match roots based on ratio being a power of omega */
        slong comparisons = 0;
        mp_limb_t ratio, inv;
        
        for (slong j = 0; j < t; j++) {
            if (first_roots[j] == 0) continue;
            
            inv = nmod_inv(first_roots[j], mod);
            
            /* Find best match */
            for (slong k = 0; k < current.length && k < t; k++) {
                if (used[k]) continue;
                comparisons++;
                
                /* Compute ratio = current_root[k] / first_root[j] */
                ratio = nmod_mul(current.pairs[k].root, inv, mod);
                
                /* Use hash table lookup */
                slong exp = omega_hash_find(&omega_table, ratio);
                if (exp >= 0) {
                    perm[j] = k;
                    used[k] = 1;
                    break;
                }
            }
        }
        
        if (DETAILED_TIMING) {
            timer_stop(&timer_match);
            printf("    Root matching (%ld comparisons): %.6f seconds\n", comparisons, timer_match.elapsed);
        }
        
        /* Fill unmatched positions with remaining roots */
        slong next_unused = 0;
        for (slong j = 0; j < t; j++) {
            if (perm[j] == -1) {
                while (next_unused < current.length && used[next_unused]) {
                    next_unused++;
                }
                if (next_unused < current.length) {
                    perm[j] = next_unused;
                    used[next_unused] = 1;
                } else {
                    /* No more roots available, use first root as fallback */
                    if (current.length > 0) {
                        perm[j] = 0;
                    }
                }
            }
        }
        
        /* Set matched roots in matrix N */
        for (slong j = 0; j < t; j++) {
            if (perm[j] >= 0 && perm[j] < current.length) {
                nmod_mat_entry(N, i, j) = current.pairs[perm[j]].root;
            } else {
                /* Use first root as default */
                if (j < t) {
                    nmod_mat_entry(N, i, j) = first_roots[j];
                }
            }
        }
        
        /* Cleanup row data */
        free(perm);
        free(used);
        
        if (current.pairs) {
            free(current.pairs);
        }
        free(current_row);
        
        if (DETAILED_TIMING) {
            timer_stop(&timer_row);
            printf("  Total row %ld: %.6f seconds\n", i, timer_row.elapsed);
        }
    }
    
    if (TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "Process all rows");
    }
    
    if (!all_matched) {
        printf("WARNING: Not all rows had matching term counts\n");
    }
    
    /* Compute exponents */
    printf("\nComputing exponents:\n");
    if (TIMING) timer_start(&timer);
    
    slong** E = (slong**)malloc(n * sizeof(slong*));
    for (slong i = 0; i < n; i++) {
        E[i] = (slong*)calloc(t, sizeof(slong));
    }
    
    mp_limb_t ratio;
    
    slong total_log_searches = 0;
    for (slong i = 0; i < n; i++) {
        for (slong j = 0; j < t; j++) {
            mp_limb_t n0j, nij, inv;
            
            n0j = nmod_mat_entry(N, 0, j);
            nij = nmod_mat_entry(N, i + 1, j);
            
            if (n0j != 0) {
                /* ratio = N[i+1,j] / N[0,j] */
                inv = nmod_inv(n0j, mod);
                ratio = nmod_mul(nij, inv, mod);
                
                /* Use hash table lookup */
                total_log_searches++;
                E[i][j] = omega_hash_find(&omega_table, ratio);
                
                if (E[i][j] < 0) {
                    E[i][j] = 0; /* Default to 0 if not found */
                }
            } else {
                E[i][j] = 0;
            }
        }
    }
    
    if (TIMING) {
        timer_stop(&timer);
        printf("  Total discrete log searches: %ld\n", total_log_searches);
        timer_print(&timer, "Compute exponents");
    }
    
    /* Build polynomial */
    printf("\nBuilding final polynomial:\n");
    if (TIMING) timer_start(&timer);
    
    nmod_mpoly_zero(result, mctx);
    mp_limb_t coeff, zeta_pow, inv;
    
    for (slong j = 0; j < t; j++) {
        /* Start with coefficient from first row */
        coeff = c[j];
        
        /* Create monomial */
        ulong* exp_vec = (ulong*)calloc(n, sizeof(ulong));
        
        /* Process each variable and adjust coefficient */
        for (slong i = 0; i < n; i++) {
            exp_vec[i] = (E[i][j] >= 0) ? E[i][j] : 0;
            
            /* Divide coefficient by zeta[i]^E[i][j] */
            if (E[i][j] > 0) {
                zeta_pow = nmod_pow_ui(zeta[i], E[i][j], mod);
                inv = nmod_inv(zeta_pow, mod);
                coeff = nmod_mul(coeff, inv, mod);
            }
        }
        
        /* Add term to polynomial */
        nmod_mpoly_push_term_ui_ui(result, coeff, exp_vec, mctx);
        
        free(exp_vec);
    }
    
    /* Combine any like terms that may have been created */
    nmod_mpoly_combine_like_terms(result, mctx);
    
    if (TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "Build polynomial");
    }
    
    /* Cleanup */
    omega_hash_clear(&omega_table);
    
    free(omega_powers);
    
    for (slong i = 0; i < n; i++) {
        free(E[i]);
    }
    free(E);
    
    free(c);
    free(first_roots);
    
    nmod_mat_clear(N);
    
    if (a.pairs) {
        free(a.pairs);
    }
    
cleanup_first:
    free(first_row);
    
    timer_stop(&timer_total);
    printf("\n=== MBOT Total Time: %.6f seconds ===\n", timer_total.elapsed);
}

/* Generate random polynomial with specified parameters */
void myrandpoly(nmod_mpoly_t f, slong n, slong T, slong D, 
                nmod_t mod, const nmod_mpoly_ctx_t mctx) {
    
    nmod_mpoly_zero(f, mctx);
    
    mp_limb_t coeff;
    
    /* Generate T distinct monomials to avoid duplicates */
    ulong** used_exps = (ulong**)malloc(T * sizeof(ulong*));
    slong num_terms = 0;
    
    while (num_terms < T) {
        /* Generate random exponent vector */
        ulong* exp = (ulong*)calloc(n, sizeof(ulong));
        for (slong j = 0; j < n; j++) {
            exp[j] = n_randint(global_state, D + 1);
        }
        
        /* Check if this monomial already exists */
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
            /* Random non-zero coefficient */
            do {
                coeff = n_randint(global_state, mod.n);
            } while (coeff == 0);
            
            /* Add term */
            nmod_mpoly_push_term_ui_ui(f, coeff, exp, mctx);
            
            /* Store exponent */
            used_exps[num_terms] = exp;
            num_terms++;
        } else {
            free(exp);
        }
    }
    
    /* Cleanup */
    for (slong i = 0; i < num_terms; i++) {
        free(used_exps[i]);
    }
    free(used_exps);
    
    /* Combine like terms if any */
    nmod_mpoly_combine_like_terms(f, mctx);
}

/* Check if field size constraint is satisfied */
int check_field_constraint(mp_limb_t p, slong n, slong T, slong d) {
    /* Check if p > 2(n+2)*T^2*d */
    mp_limb_t bound = 2 * (n + 2) * T * T * d;
    
    int result = p > bound;
    
    if (!result) {
        printf("\n*** WARNING: Field size constraint violated! ***\n");
        printf("Field size p = %lu\n", p);
        printf("Required: p > 2(n+2)*T^2*d = %lu\n", bound);
        printf("where n=%ld, T=%ld, d=%ld\n", n, T, d);
    }
    
    return result;
}

/* Find primitive root of prime p */
void find_primitive_root(mp_limb_t* omega, mp_limb_t p) {
    /* For small known primes, use known primitive roots */
    if (p == 101) {
        *omega = 2;
        return;
    }
    if (p == 65537) {
        *omega = 3;
        return;
    }
    
    /* Otherwise, use FLINT's primitive root finder */
    nmod_t mod;
    nmod_init(&mod, p);
    *omega = n_primitive_root_prime(p);
}

/* Main interpolation algorithm with constraint checking */
void ComputePolyMatrixDet(nmod_mpoly_t det_poly, nmod_mpoly_t actual_det,
                         const poly_mat_t* A, slong n, mp_limb_t p,
                         const nmod_mpoly_ctx_t mctx) {
    my_timer_t timer_det;
    if (DEBUG) printf("\n========== Computing Polynomial Matrix Determinant ==========\n");
    
    /* Set up modulus */
    nmod_t mod;
    nmod_init(&mod, p);
    
    /* Find primitive root */
    mp_limb_t omega;
    find_primitive_root(&omega, p);
    
    if (DEBUG) {
        printf("Prime p = %lu, Primitive root omega = %lu\n", p, omega);
    }
    
    /* Generate random values */
    mp_limb_t* zeta = (mp_limb_t*)malloc(n * sizeof(mp_limb_t));
    mp_limb_t* alpha = (mp_limb_t*)malloc(n * sizeof(mp_limb_t));
    
    for (slong i = 0; i < n; i++) {
        zeta[i] = n_randint(global_state, p);
        alpha[i] = n_randint(global_state, p);
        /* Ensure non-zero */
        if (zeta[i] == 0) zeta[i] = 1;
        if (alpha[i] == 0) alpha[i] = 1;
    }
    
    if (DEBUG) {
        printf("Random zeta = (");
        for (slong i = 0; i < n; i++) {
            printf("%lu", zeta[i]);
            if (i < n-1) printf(",");
        }
        printf(")\n");
        printf("Random alpha = (");
        for (slong i = 0; i < n; i++) {
            printf("%lu", alpha[i]);
            if (i < n-1) printf(",");
        }
        printf(")\n");
    }
    
    if (TIMING) timer_start(&timer_det);
    /* Compute actual determinant */
    poly_mat_det(actual_det, A, mctx);
    if (TIMING) {
        timer_stop(&timer_det);
        timer_print(&timer_det, "    Direct determinant Computing");
    }
    
    slong T = nmod_mpoly_length(actual_det, mctx);
    if (T == 0) {
        nmod_mpoly_zero(det_poly, mctx);
        goto cleanup;
    }
    
    /* Get maximum degrees per variable */
    slong* max_degs = (slong*)malloc(n * sizeof(slong));
    poly_max_degrees_per_var(max_degs, actual_det, mctx);
    
    /* Get maximum total degree for field size check */
    slong d = 0;
    for (slong i = 0; i < n; i++) {
        d += max_degs[i];
    }
    
    if (DEBUG) {
        printf("\nActual determinant has %ld terms\n", T);
        printf("Max degrees per variable: ");
        for (slong i = 0; i < n; i++) {
            printf("%ld ", max_degs[i]);
        }
        printf("\nMax total degree: %ld\n", d);
    }

    if (TIMING) timer_start(&timer_det);
    
    /* Check field size constraint */
    check_field_constraint(p, n, T, d);
    
    /* Apply diversification transformation */
    nmod_mpoly_t det_transformed;
    nmod_mpoly_init(det_transformed, mctx);
    Mydiver(det_transformed, actual_det, n, zeta, mod, mctx);
    
    /* Build evaluation matrix */
    nmod_mat_t M;
    nmod_mat_init(M, n + 1, 2*T, p);
    TotalMyeval(M, det_transformed, n, alpha, omega, mod, mctx);
    
    if (DEBUG) {
        printf("\nEvaluation matrix M:\n");
        for (slong i = 0; i <= n; i++) {
            printf("Row %ld: ", i);
            for (slong j = 0; j < 2*T && j < 10; j++) {
                printf("%lu ", nmod_mat_entry(M, i, j));
            }
            if (2*T > 10) printf("...");
            printf("\n");
        }
    }
    
    /* Apply MBOT algorithm with per-variable degree bounds */
    MBOT(det_poly, M, n, T, omega, zeta, mod, max_degs, mctx);

    if (TIMING) {
        timer_stop(&timer_det);
        timer_print(&timer_det, "    MBOT determinant Computing");
    }
    
    /* Cleanup */
    free(max_degs);
    nmod_mpoly_clear(det_transformed, mctx);
    nmod_mat_clear(M);
    
cleanup:
    free(zeta);
    free(alpha);
}

/* Test function for random polynomial */
void test_random_polynomial(void) {
    printf("\n========== Testing Random Polynomial Interpolation ==========\n");
    
    /* Parameters */
    slong n = 3;      /* Number of variables */
    slong T = 8000;   /* Number of terms */
    slong d = 1000;   /* Maximum degree */
    
    /* Calculate required prime size */
    mp_limb_t p = 2 * (n + 2) * T * T * d;
    
    /* Find next prime */
    p = n_nextprime(p, 1);
    
    printf("Parameters: n=%ld, T=%ld, d=%ld\n", n, T, d);
    printf("Minimum prime required: %lu\n", 2 * (n + 2) * T * T * d);
    printf("Using prime: %lu\n", p);
    
    /* Initialize contexts */
    nmod_t mod;
    nmod_init(&mod, p);
    
    char** vars = (char**)malloc(n * sizeof(char*));
    for (slong i = 0; i < n; i++) {
        vars[i] = (char*)malloc(10);
        sprintf(vars[i], "x%ld", i);
    }
    
    nmod_mpoly_ctx_t mctx;
    nmod_mpoly_ctx_init(mctx, n, ORD_LEX, p);
    
    /* Generate random polynomial */
    nmod_mpoly_t f;
    nmod_mpoly_init(f, mctx);
    myrandpoly(f, n, T, d, mod, mctx);
    
    slong actual_terms = nmod_mpoly_length(f, mctx);
    printf("\nGenerated polynomial with %ld terms (requested %ld):\n", actual_terms, T);
    if (actual_terms <= 20) {
        nmod_mpoly_print_pretty(f, (const char**)vars, mctx);
        printf("\n");
    } else {
        printf("(too large to display)\n");
    }
    
    /* Find primitive root */
    mp_limb_t omega;
    find_primitive_root(&omega, p);
    
    /* Generate random transformation values */
    mp_limb_t* zeta = (mp_limb_t*)malloc(n * sizeof(mp_limb_t));
    mp_limb_t* alpha = (mp_limb_t*)malloc(n * sizeof(mp_limb_t));
    
    for (slong i = 0; i < n; i++) {
        zeta[i] = n_randint(global_state, p);
        alpha[i] = n_randint(global_state, p);
        if (zeta[i] == 0) zeta[i] = 1;
        if (alpha[i] == 0) alpha[i] = 1;
    }
    
    /* Apply diversification */
    nmod_mpoly_t f1;
    nmod_mpoly_init(f1, mctx);
    Mydiver(f1, f, n, zeta, mod, mctx);
    
    /* Build evaluation matrix */
    nmod_mat_t M;
    nmod_mat_init(M, n + 1, 2*actual_terms, p);
    TotalMyeval(M, f1, n, alpha, omega, mod, mctx);
    
    /* Apply MBOT to recover polynomial */
    nmod_mpoly_t g;
    nmod_mpoly_init(g, mctx);
    
    /* Get max degrees per variable */
    slong* max_degs_per_var = (slong*)malloc(n * sizeof(slong));
    poly_max_degrees_per_var(max_degs_per_var, f, mctx);
    
    printf("Max degrees per variable: ");
    for (slong i = 0; i < n; i++) {
        printf("%ld ", max_degs_per_var[i]);
    }
    printf("\n");
    
    MBOT(g, M, n, actual_terms, omega, zeta, mod, max_degs_per_var, mctx);
    
    /* Check if recovery was successful */
    printf("\nRecovered polynomial with %ld terms\n", nmod_mpoly_length(g, mctx));
    
    /* Sort both polynomials to compare */
    nmod_mpoly_sort_terms(f, mctx);
    nmod_mpoly_sort_terms(g, mctx);
    
    if (nmod_mpoly_equal(f, g, mctx)) {
        printf("Success: Polynomials match!\n");
    } else {
        printf("Error: Polynomials don't match!\n");
        if (nmod_mpoly_length(g, mctx) <= 20) {
            printf("Recovered: ");
            nmod_mpoly_print_pretty(g, (const char**)vars, mctx);
            printf("\n");
        }
        
        /* Debug: Check for differences */
        nmod_mpoly_t diff;
        nmod_mpoly_init(diff, mctx);
        nmod_mpoly_sub(diff, f, g, mctx);
        printf("Difference has %ld terms\n", nmod_mpoly_length(diff, mctx));
        if (nmod_mpoly_length(diff, mctx) <= 10) {
            printf("Missing/extra terms: ");
            nmod_mpoly_print_pretty(diff, (const char**)vars, mctx);
            printf("\n");
        }
        nmod_mpoly_clear(diff, mctx);
    }
    
    /* Cleanup */
    nmod_mpoly_clear(f, mctx);
    nmod_mpoly_clear(f1, mctx);
    nmod_mpoly_clear(g, mctx);
    nmod_mat_clear(M);
    
    free(zeta);
    free(alpha);
    free(max_degs_per_var);
    
    nmod_mpoly_ctx_clear(mctx);
    
    for (slong i = 0; i < n; i++) {
        free(vars[i]);
    }
    free(vars);
}

/* Test function */
int huang_test(void) {
    /* Version check */
    printf("Code version: nmod version with FLINT root finding\n");
    
    /* Initialize random state */
    flint_rand_init(global_state);
    
    /* First run the original test */
    printf("========== Testing Polynomial Matrix Determinant ==========\n");
    
    /* Set prime and dimensions */
    mp_limb_t p = 65537; /* Use a moderate prime */
    slong n = 3; /* Number of variables */
    slong k = 2; /* Matrix size (2x2) */
    
    /* Initialize contexts */
    char** vars = (char**)malloc(n * sizeof(char*));
    for (slong i = 0; i < n; i++) {
        vars[i] = (char*)malloc(10);
        sprintf(vars[i], "x%ld", i);
    }
    
    nmod_mpoly_ctx_t mctx;
    nmod_mpoly_ctx_init(mctx, n, ORD_LEX, p);
    
    /* Create a test polynomial matrix */
    poly_mat_t A;
    poly_mat_init(&A, k, k, mctx);
    
    /* Fill with specific polynomials for testing */
    nmod_mpoly_t entry;
    nmod_mpoly_init(entry, mctx);
    
    /* A[0,0] = 3*x1 + 22*x0^3*x1 */
    nmod_mpoly_zero(entry, mctx);
    ulong* exp = (ulong*)calloc(n, sizeof(ulong));
    
    exp[0] = 0; exp[1] = 1; exp[2] = 0; /* x1 */
    nmod_mpoly_push_term_ui_ui(entry, 3, exp, mctx);
    
    exp[0] = 3; exp[1] = 1; exp[2] = 0; /* x0^3*x1 */
    nmod_mpoly_push_term_ui_ui(entry, 22, exp, mctx);
    poly_mat_entry_set(&A, 0, 0, entry, mctx);
    
    /* A[0,1] = 64*x0*x1^2 */
    nmod_mpoly_zero(entry, mctx);
    exp[0] = 1; exp[1] = 2; exp[2] = 0; /* x0*x1^2 */
    nmod_mpoly_push_term_ui_ui(entry, 64, exp, mctx);
    poly_mat_entry_set(&A, 0, 1, entry, mctx);
    
    /* A[1,0] = 61*x0^2 + 91*x0^2*x1 */
    nmod_mpoly_zero(entry, mctx);
    exp[0] = 2; exp[1] = 0; exp[2] = 0; /* x0^2 */
    nmod_mpoly_push_term_ui_ui(entry, 61, exp, mctx);
    
    exp[0] = 2; exp[1] = 1; exp[2] = 0; /* x0^2*x1 */
    nmod_mpoly_push_term_ui_ui(entry, 91, exp, mctx);
    poly_mat_entry_set(&A, 1, 0, entry, mctx);
    
    /* A[1,1] = 87*x0*x1^3 + 26*x1^2*x2^4 + 89*x1^3*x2^2 */
    nmod_mpoly_zero(entry, mctx);
    exp[0] = 1; exp[1] = 3; exp[2] = 0; /* x0*x1^3 */
    nmod_mpoly_push_term_ui_ui(entry, 87, exp, mctx);
    
    exp[0] = 0; exp[1] = 2; exp[2] = 4; /* x1^2*x2^4 */
    nmod_mpoly_push_term_ui_ui(entry, 26, exp, mctx);
    
    exp[0] = 0; exp[1] = 3; exp[2] = 2; /* x1^3*x2^2 */
    nmod_mpoly_push_term_ui_ui(entry, 89, exp, mctx);
    poly_mat_entry_set(&A, 1, 1, entry, mctx);
    
    free(exp);
    nmod_mpoly_clear(entry, mctx);
    
    /* Compute determinant */
    nmod_mpoly_t det_poly, actual_det;
    nmod_mpoly_init(det_poly, mctx);
    nmod_mpoly_init(actual_det, mctx);
    
    ComputePolyMatrixDet(det_poly, actual_det, &A, n, p, mctx);
    
    /* Sort before comparison */
    nmod_mpoly_sort_terms(det_poly, mctx);
    nmod_mpoly_sort_terms(actual_det, mctx);
    
    /* Print results */
    printf("\n========== FINAL RESULTS ==========\n");
    printf("Computed determinant:\n");
    nmod_mpoly_print_pretty(det_poly, (const char**)vars, mctx);
    printf("\n\nActual determinant:\n");
    nmod_mpoly_print_pretty(actual_det, (const char**)vars, mctx);
    printf("\n");
    
    /* Check if they're equal */
    if (nmod_mpoly_equal(det_poly, actual_det, mctx)) {
        printf("\nSuccess: Determinants match!\n");
    } else {
        printf("\nError: Determinants don't match!\n");
    }
    
    /* Cleanup */
    nmod_mpoly_clear(det_poly, mctx);
    nmod_mpoly_clear(actual_det, mctx);
    poly_mat_clear(&A, mctx);
    nmod_mpoly_ctx_clear(mctx);
    
    for (slong i = 0; i < n; i++) {
        free(vars[i]);
    }
    free(vars);
    
    /* Now run the random polynomial test */
    test_random_polynomial();
    
    flint_rand_clear(global_state);
    flint_cleanup();
    
    return 0;
}