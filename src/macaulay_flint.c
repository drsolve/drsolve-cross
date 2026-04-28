/*
 * macaulay_flint.c - Macaulay resultant construction for finite extension fields
 */

#include "macaulay_flint.h"
#include "dixon_flint.h"

static const char *macaulay_det_method_name(det_method_t method)
{
    switch (method) {
        case DET_METHOD_RECURSIVE:
            return "recursive expansion";
        case DET_METHOD_KRONECKER:
            return "Kronecker+HNF";
        case DET_METHOD_INTERPOLATION:
            return "interpolation";
        case DET_METHOD_HUANG:
            return "sparse interpolation";
        default:
            return "default";
    }
}

static slong macaulay_binomial(slong n, slong k)
{
    if (k > n || k < 0) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;

    slong result = 1;
    for (slong i = 0; i < k; i++) {
        result = result * (n - i) / (i + 1);
    }
    return result;
}

static slong compute_fq_polynomial_elimination_degree(const fq_mvpoly_t *poly, slong nvars)
{
    slong max_degree = 0;

    if (poly == NULL) {
        return 0;
    }

    for (slong t = 0; t < poly->nterms; t++) {
        slong degree = 0;
        if (poly->terms[t].var_exp) {
            for (slong v = 0; v < nvars && v < poly->nvars; v++) {
                degree += poly->terms[t].var_exp[v];
            }
        }
        if (degree > max_degree) {
            max_degree = degree;
        }
    }

    return max_degree;
}

static void compute_macaulay_dimensions(const fq_mvpoly_t *polys,
                                        slong npolys,
                                        slong nvars,
                                        slong *elim_degrees,
                                        slong *macaulay_degree_out,
                                        slong *nrows_out,
                                        slong *ncols_out)
{
    slong macaulay_degree = 0;
    slong nrows = 0;
    slong ncols = 0;

    for (slong i = 0; i < npolys; i++) {
        elim_degrees[i] = compute_fq_polynomial_elimination_degree(&polys[i], nvars);
        macaulay_degree += elim_degrees[i];
    }
    macaulay_degree -= nvars;
    if (macaulay_degree < 0) {
        macaulay_degree = 0;
    }

    ncols = macaulay_binomial(nvars + macaulay_degree, macaulay_degree);
    for (slong i = 0; i < npolys; i++) {
        slong multiplier_degree = macaulay_degree - elim_degrees[i];
        nrows += macaulay_binomial(nvars + multiplier_degree, multiplier_degree);
    }

    *macaulay_degree_out = macaulay_degree;
    *nrows_out = nrows;
    *ncols_out = ncols;
}

static ulong hash_monomial_exp(const slong *exp, slong nvars)
{
    ulong h = 1469598103934665603ULL;
    for (slong i = 0; i < nvars; i++) {
        h ^= (ulong) (exp[i] + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
        h *= 1099511628211ULL;
    }
    return h;
}

static void macaulay_fill_monomials_leq_recursive(monom_t *monoms,
                                                  slong *count,
                                                  slong *storage,
                                                  slong nvars,
                                                  slong pos,
                                                  slong remaining,
                                                  slong *current_exp)
{
    if (pos == nvars) {
        slong *dest = &storage[(*count) * nvars];
        memcpy(dest, current_exp, (size_t) nvars * sizeof(slong));
        monoms[*count].exp = dest;
        monoms[*count].idx = *count;
        (*count)++;
        return;
    }

    for (slong d = 0; d <= remaining; d++) {
        current_exp[pos] = d;
        macaulay_fill_monomials_leq_recursive(monoms, count, storage,
                                              nvars, pos + 1, remaining - d,
                                              current_exp);
    }
}

static void enumerate_monomials_degree_leq(monom_t **monoms_out,
                                           slong *count_out,
                                           slong nvars,
                                           slong degree_bound)
{
    slong total_monomials;
    slong structs_size;
    slong exps_size;
    char *combined;
    monom_t *monoms;
    slong *storage;
    slong *current_exp;
    slong count = 0;

    if (degree_bound < 0) {
        *monoms_out = NULL;
        *count_out = 0;
        return;
    }

    total_monomials = macaulay_binomial(nvars + degree_bound, degree_bound);
    structs_size = total_monomials * sizeof(monom_t);
    exps_size = total_monomials * nvars * sizeof(slong);
    combined = (char *) flint_malloc((size_t) structs_size + (size_t) exps_size);
    monoms = (monom_t *) combined;
    storage = (slong *) (combined + structs_size);
    current_exp = (slong *) flint_calloc((size_t) (nvars > 0 ? nvars : 1), sizeof(slong));

    macaulay_fill_monomials_leq_recursive(monoms, &count, storage,
                                          nvars, 0, degree_bound, current_exp);

    flint_free(current_exp);
    *monoms_out = monoms;
    *count_out = count;
}

static hash_entry_t **build_monomial_hash_table(const monom_t *monoms,
                                                slong count,
                                                slong nvars,
                                                slong *hash_size_out)
{
    slong hash_size = 1024;
    hash_entry_t **buckets;

    while (hash_size < 2 * FLINT_MAX(1, count)) {
        hash_size <<= 1;
    }

    buckets = (hash_entry_t **) flint_calloc((size_t) hash_size, sizeof(hash_entry_t *));
    for (slong i = 0; i < count; i++) {
        ulong hash = hash_monomial_exp(monoms[i].exp, nvars) & (ulong) (hash_size - 1);
        hash_entry_t *entry = (hash_entry_t *) flint_malloc(sizeof(hash_entry_t));
        entry->exp = monoms[i].exp;
        entry->idx = monoms[i].idx;
        entry->next = buckets[hash];
        buckets[hash] = entry;
    }

    *hash_size_out = hash_size;
    return buckets;
}

static void free_hash_table_entries(hash_entry_t **buckets, slong hash_size)
{
    if (buckets == NULL) {
        return;
    }
    for (slong i = 0; i < hash_size; i++) {
        hash_entry_t *entry = buckets[i];
        while (entry) {
            hash_entry_t *next = entry->next;
            flint_free(entry);
            entry = next;
        }
    }
    flint_free(buckets);
}

static slong lookup_monomial_index(hash_entry_t **buckets,
                                   slong hash_size,
                                   const slong *exp,
                                   slong nvars)
{
    ulong hash = hash_monomial_exp(exp, nvars) & (ulong) (hash_size - 1);
    hash_entry_t *entry = buckets[hash];

    while (entry) {
        if (memcmp(entry->exp, exp, (size_t) nvars * sizeof(slong)) == 0) {
            return entry->idx;
        }
        entry = entry->next;
    }

    return -1;
}

static void extract_fq_coefficient_matrix_from_macaulay(fq_mvpoly_t ***coeff_matrix,
                                                        slong *matrix_size,
                                                        const fq_mvpoly_t *polys,
                                                        slong npolys,
                                                        slong nvars,
                                                        slong npars,
                                                        slong *full_rows_out,
                                                        slong *full_cols_out,
                                                        slong *macaulay_degree_out)
{
    const fq_nmod_ctx_struct *ctx = polys[0].ctx;
    double step2_wall_start;
    double step3_wall_start;
    slong *elim_degrees = (slong *) flint_calloc((size_t) npolys, sizeof(slong));
    slong macaulay_degree = 0;
    fq_mvpoly_t ***full_matrix = NULL;
    monom_t *column_monoms = NULL;
    slong ncols = 0;
    slong nrows = 0;
    hash_entry_t **col_hash = NULL;
    slong col_hash_size = 0;
    slong row_counter = 0;

    compute_macaulay_dimensions(polys, npolys, nvars, elim_degrees,
                                &macaulay_degree, &nrows, &ncols);

    enumerate_monomials_degree_leq(&column_monoms, &ncols, nvars, macaulay_degree);

    printf("\nStep 2: Construct Macaulay matrix\n");
    step2_wall_start = get_wall_time();
    printf("Macaulay degree: %ld\n", macaulay_degree);
    printf("Macaulay matrix size: %ld x %ld\n", nrows, ncols);

    *full_rows_out = nrows;
    *full_cols_out = ncols;
    *macaulay_degree_out = macaulay_degree;

    if (nrows == 0 || ncols == 0) {
        printf("Warning: Empty Macaulay coefficient matrix\n");
        *coeff_matrix = NULL;
        *matrix_size = 0;
        dixon_maybe_print_step_time("Step 2", get_wall_time() - step2_wall_start);
        flint_free(elim_degrees);
        if (column_monoms) flint_free(column_monoms);
        return;
    }

    full_matrix = (fq_mvpoly_t ***) flint_malloc((size_t) nrows * sizeof(fq_mvpoly_t **));
    for (slong i = 0; i < nrows; i++) {
        full_matrix[i] = (fq_mvpoly_t **) flint_calloc((size_t) ncols, sizeof(fq_mvpoly_t *));
    }

    col_hash = build_monomial_hash_table(column_monoms, ncols, nvars, &col_hash_size);

    for (slong poly_idx = 0; poly_idx < npolys; poly_idx++) {
        monom_t *multipliers = NULL;
        slong num_multipliers = 0;
        slong multiplier_degree = macaulay_degree - elim_degrees[poly_idx];

        enumerate_monomials_degree_leq(&multipliers, &num_multipliers, nvars, multiplier_degree);

        for (slong mi = 0; mi < num_multipliers; mi++) {
            for (slong t = 0; t < polys[poly_idx].nterms; t++) {
                slong *col_exp = (slong *) flint_calloc((size_t) (nvars > 0 ? nvars : 1), sizeof(slong));
                slong col = -1;

                for (slong v = 0; v < nvars; v++) {
                    slong term_exp = (polys[poly_idx].terms[t].var_exp && v < polys[poly_idx].nvars)
                                   ? polys[poly_idx].terms[t].var_exp[v]
                                   : 0;
                    col_exp[v] = multipliers[mi].exp[v] + term_exp;
                }

                col = lookup_monomial_index(col_hash, col_hash_size, col_exp, nvars);
                if (col >= 0) {
                    fq_mvpoly_t *entry = get_matrix_entry_lazy(full_matrix, row_counter, col,
                                                               npars, ctx);
                    fq_mvpoly_add_term(entry, NULL,
                                       polys[poly_idx].terms[t].par_exp,
                                       polys[poly_idx].terms[t].coeff);
                }

                flint_free(col_exp);
            }
            row_counter++;
        }

        if (multipliers) flint_free(multipliers);
    }
    dixon_maybe_print_step_time("Step 2", get_wall_time() - step2_wall_start);

    printf("\nStep 3: Extract maximal-rank submatrix\n");
    step3_wall_start = get_wall_time();
    {
        slong *row_idx_array = NULL;
        slong *col_idx_array = NULL;
        slong submat_rows = 0;
        slong submat_cols = 0;
        slong submat_rank = 0;

        find_fq_optimal_maximal_rank_submatrix(full_matrix, nrows, ncols,
                                               &row_idx_array, &col_idx_array,
                                               &submat_rows, &submat_cols, npars);

        submat_rank = FLINT_MIN(submat_rows, submat_cols);

        printf("Submatrix size: %ld x %ld\n", submat_rank, submat_rank);

        if (submat_rank <= 0) {
            *coeff_matrix = NULL;
            *matrix_size = 0;
        } else {
            *coeff_matrix = (fq_mvpoly_t **) flint_malloc((size_t) submat_rank * sizeof(fq_mvpoly_t *));
            for (slong i = 0; i < submat_rank; i++) {
                (*coeff_matrix)[i] = (fq_mvpoly_t *) flint_malloc((size_t) submat_rank * sizeof(fq_mvpoly_t));
                for (slong j = 0; j < submat_rank; j++) {
                    fq_mvpoly_t *source = full_matrix[row_idx_array[i]][col_idx_array[j]];
                    if (source != NULL) {
                        fq_mvpoly_copy(&(*coeff_matrix)[i][j], source);
                    } else {
                        fq_mvpoly_init(&(*coeff_matrix)[i][j], 0, npars, ctx);
                    }
                }
            }
            *matrix_size = submat_rank;
        }

        if (row_idx_array) flint_free(row_idx_array);
        if (col_idx_array) flint_free(col_idx_array);
    }
    dixon_maybe_print_step_time("Step 3", get_wall_time() - step3_wall_start);

    for (slong i = 0; i < nrows; i++) {
        if (full_matrix[i]) {
            for (slong j = 0; j < ncols; j++) {
                if (full_matrix[i][j]) {
                    fq_mvpoly_clear(full_matrix[i][j]);
                    flint_free(full_matrix[i][j]);
                }
            }
            flint_free(full_matrix[i]);
        }
    }
    flint_free(full_matrix);
    free_hash_table_entries(col_hash, col_hash_size);
    flint_free(elim_degrees);
    if (column_monoms) flint_free(column_monoms);
}

static void clear_macaulay_coeff_matrix(fq_mvpoly_t **coeff_matrix, slong matrix_size)
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

static det_method_t choose_macaulay_det_method(slong matrix_size, slong npars)
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
        printf("Macaulay determinant method override active: %d (%s)\n",
               dixon_global_method_step4, macaulay_det_method_name(dixon_global_method_step4));
    } else if (dixon_global_method_step1 != -1) {
        coeff_method = dixon_global_method_step1;
        printf("Macaulay determinant method override active: %d (%s)\n",
               dixon_global_method_step1, macaulay_det_method_name(dixon_global_method_step1));
    }

    return coeff_method;
}

void fq_macaulay_resultant(fq_mvpoly_t *result, fq_mvpoly_t *polys,
                           slong nvars, slong npars)
{
    cleanup_unified_workspace();

    printf("\nStep 1: Build Macaulay coefficient matrix\n");
    double step1_wall_start = get_wall_time();

    fq_mvpoly_t **coeff_matrix = NULL;
    slong macaulay_degree = 0;
    slong full_rows = 0, full_cols = 0;
    slong matrix_size = 0;

    extract_fq_coefficient_matrix_from_macaulay(&coeff_matrix, &matrix_size,
                                                polys, nvars + 1, nvars, npars,
                                                &full_rows, &full_cols, &macaulay_degree);
    dixon_maybe_print_step_time("Step 1", get_wall_time() - step1_wall_start);

    if (matrix_size > 0) {
        det_method_t coeff_method;
        slong res_deg_bound = compute_fq_dixon_resultant_degree_bound(polys, nvars + 1, nvars, npars);
        clock_t step4_cpu_start;
        double step4_wall_start;

        printf("\nStep 4: Compute Macaulay resultant\n");
        printf("Resultant degree bound: %ld\n", res_deg_bound);

        coeff_method = choose_macaulay_det_method(matrix_size, npars);
        step4_cpu_start = clock();
        step4_wall_start = get_wall_time();
        compute_fq_coefficient_matrix_det(result, coeff_matrix, matrix_size,
                                          npars, polys[0].ctx, coeff_method, res_deg_bound);
        dixon_maybe_print_step_method_time("Step 4",
                                           coeff_method,
                                           ((double)(clock() - step4_cpu_start) / CLOCKS_PER_SEC),
                                           get_wall_time() - step4_wall_start);

        if (result->nterms <= 100) {
            fq_mvpoly_print(result, "Final Resultant");
        } else {
            printf("Final resultant too large to display (%ld terms)\n", result->nterms);
        }
        fq_mvpoly_make_monic(result);
        clear_macaulay_coeff_matrix(coeff_matrix, matrix_size);
    } else {
        fq_mvpoly_init(result, 0, npars, polys[0].ctx);
        printf("Warning: Empty Macaulay coefficient matrix, resultant is 0\n");
    }

    printf("\n=== Macaulay Resultant Computation Complete ===\n");
}

void fq_macaulay_resultant_with_names(fq_mvpoly_t *result, fq_mvpoly_t *polys,
                                      slong nvars, slong npars,
                                      char **var_names, char **par_names,
                                      const char *gen_name)
{
    (void) var_names;
    cleanup_unified_workspace();

    printf("\nStep 1: Build Macaulay coefficient matrix\n");
    double step1_wall_start = get_wall_time();

    fq_mvpoly_t **coeff_matrix = NULL;
    slong macaulay_degree = 0;
    slong full_rows = 0, full_cols = 0;
    slong matrix_size = 0;

    extract_fq_coefficient_matrix_from_macaulay(&coeff_matrix, &matrix_size,
                                                polys, nvars + 1, nvars, npars,
                                                &full_rows, &full_cols, &macaulay_degree);
    dixon_maybe_print_step_time("Step 1", get_wall_time() - step1_wall_start);

    if (matrix_size > 0) {
        det_method_t coeff_method;
        slong res_deg_bound = compute_fq_dixon_resultant_degree_bound(polys, nvars + 1, nvars, npars);
        clock_t step4_cpu_start;
        double step4_wall_start;

        printf("\nStep 4: Compute Macaulay resultant\n");
        printf("Resultant degree bound: %ld\n", res_deg_bound);

        coeff_method = choose_macaulay_det_method(matrix_size, npars);
        step4_cpu_start = clock();
        step4_wall_start = get_wall_time();
        compute_fq_coefficient_matrix_det(result, coeff_matrix, matrix_size,
                                          npars, polys[0].ctx, coeff_method, res_deg_bound);
        dixon_maybe_print_step_method_time("Step 4",
                                           coeff_method,
                                           ((double)(clock() - step4_cpu_start) / CLOCKS_PER_SEC),
                                           get_wall_time() - step4_wall_start);

        if (result->nterms <= 100) {
            fq_mvpoly_print_with_names(result, "Final Resultant", NULL, par_names, gen_name, 0);
        } else {
            printf("Final resultant too large to display (%ld terms)\n", result->nterms);
        }
        fq_mvpoly_make_monic(result);
        clear_macaulay_coeff_matrix(coeff_matrix, matrix_size);
    } else {
        fq_mvpoly_init(result, 0, npars, polys[0].ctx);
        printf("Warning: Empty Macaulay coefficient matrix, resultant is 0\n");
    }

    printf("\n=== Macaulay Resultant Computation Complete ===\n");
}
