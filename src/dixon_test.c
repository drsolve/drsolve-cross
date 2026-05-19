// dixon_test.c

#include "dixon_test.h"

#include "large_prime_system_solver.h"

#include <fcntl.h>
#include <time.h>
#include <unistd.h>

// Global threshold variable (extern declaration should be in header)
extern int g_matrix_transpose_threshold;

// ============= DEGREE CHECKING STRUCTURES =============

typedef struct {
    char test_name[256];
    slong nvars;
    slong npars;
    slong npolys;
    slong *input_degrees;
    slong resultant_degree;
    slong expected_degree_product;
    int degree_satisfies_bound;
    int algorithm_type; // 0 = traditional, 1 = degree-optimal
    int test_case_index; // Index of test case (0-4 for the 5 test cases)
} degree_check_result_t;

// Global storage for degree check results
static degree_check_result_t *degree_results = NULL;
static slong num_degree_results = 0;
static slong degree_results_capacity = 0;
static int g_bezout_test_quiet_mode = 0;

typedef struct {
    const char *label;
    slong nvars;
    slong npolys;
    slong degrees[8];
} solver_degree_case_t;

typedef struct {
    const char *label;
    int rational_mode;
    int large_prime_mode;
    ulong prime_ui;
    ulong power;
    const char *field_poly;
    const char *generator;
    const char *prime_str;
} solver_field_case_t;

typedef struct {
    slong passed_cases;
    slong no_solution_cases;
    slong failed_cases;
} solver_correctness_stats_t;

static int dtest_append_text(char **buffer, size_t *capacity, size_t *length, const char *text)
{
    size_t text_len;
    char *grown;

    if (!buffer || !capacity || !length || !text)
        return 0;

    text_len = strlen(text);
    if (*buffer == NULL || *capacity == 0) {
        *capacity = text_len + 32;
        *buffer = (char *) malloc(*capacity);
        if (!*buffer)
            return 0;
        (*buffer)[0] = '\0';
        *length = 0;
    }

    if (*length + text_len + 1 > *capacity) {
        while (*length + text_len + 1 > *capacity) {
            *capacity *= 2;
        }
        grown = (char *) realloc(*buffer, *capacity);
        if (!grown)
            return 0;
        *buffer = grown;
    }

    memcpy(*buffer + *length, text, text_len + 1);
    *length += text_len;
    return 1;
}

static void dtest_redirect_fd_to_devnull(int fd, int *saved_fd)
{
    int devnull;

    if (!saved_fd)
        return;

    *saved_fd = dup(fd);
    if (*saved_fd == -1)
        return;

    devnull = open("/dev/null", O_WRONLY);
    if (devnull == -1) {
        close(*saved_fd);
        *saved_fd = -1;
        return;
    }

    dup2(devnull, fd);
    close(devnull);
}

static void dtest_restore_fd(int fd, int saved_fd)
{
    if (saved_fd == -1)
        return;

    fflush(fd == STDOUT_FILENO ? stdout : stderr);
    dup2(saved_fd, fd);
    close(saved_fd);
}

static void dtest_free_enumerated_monomials(monomial_t *monomials, slong count)
{
    if (!monomials)
        return;
    for (slong i = 0; i < count; i++)
        free(monomials[i].exponents);
    free(monomials);
}

static int dtest_build_random_system_strings(slong nvars,
                                             slong num_elim_vars,
                                             char **elim_vars_out,
                                             char **all_vars_out)
{
    char *elim_vars = NULL;
    char *all_vars = NULL;
    size_t elim_cap = 0, elim_len = 0;
    size_t all_cap = 0, all_len = 0;

    if (!elim_vars_out || !all_vars_out || nvars <= 0 || num_elim_vars < 0 || num_elim_vars > nvars)
        return 0;

    for (slong i = 0; i < nvars; i++) {
        char name[32];
        snprintf(name, sizeof(name), "x%ld", i);
        if (i > 0 && !dtest_append_text(&all_vars, &all_cap, &all_len, ","))
            goto fail;
        if (!dtest_append_text(&all_vars, &all_cap, &all_len, name))
            goto fail;

        if (i < num_elim_vars) {
            if (i > 0 && !dtest_append_text(&elim_vars, &elim_cap, &elim_len, ","))
                goto fail;
            if (!dtest_append_text(&elim_vars, &elim_cap, &elim_len, name))
                goto fail;
        }
    }

    if (!elim_vars)
        elim_vars = strdup("");
    if (!all_vars)
        all_vars = strdup("");
    if (!elim_vars || !all_vars)
        goto fail;

    *elim_vars_out = elim_vars;
    *all_vars_out = all_vars;
    return 1;

fail:
    free(elim_vars);
    free(all_vars);
    return 0;
}

static int dtest_append_fmpz_monomial_text(char **buffer,
                                           size_t *capacity,
                                           size_t *length,
                                           const fmpz_t coeff,
                                           const slong *exp,
                                           slong nvars,
                                           int *first_term)
{
    int has_var = 0;
    char *coeff_str = NULL;

    if (!buffer || !capacity || !length || !first_term)
        return 0;
    if (fmpz_is_zero(coeff))
        return 0;

    for (slong j = 0; j < nvars; j++) {
        if (exp[j] != 0) {
            has_var = 1;
            break;
        }
    }

    if (!*first_term) {
        if (!dtest_append_text(buffer, capacity, length, " + "))
            return 0;
    }

    if (!fmpz_is_one(coeff) || !has_var) {
        coeff_str = fmpz_get_str(NULL, 10, coeff);
        if (!coeff_str)
            return 0;
        if (!dtest_append_text(buffer, capacity, length, coeff_str)) {
            flint_free(coeff_str);
            return 0;
        }
        flint_free(coeff_str);
    }

    if (has_var) {
        int wrote_any_var = 0;
        for (slong j = 0; j < nvars; j++) {
            char piece[64];
            if (exp[j] == 0)
                continue;
            if ((!fmpz_is_one(coeff)) || wrote_any_var) {
                if (!dtest_append_text(buffer, capacity, length, "*"))
                    return 0;
            }
            snprintf(piece, sizeof(piece), "x%ld", j);
            if (!dtest_append_text(buffer, capacity, length, piece))
                return 0;
            if (exp[j] != 1) {
                snprintf(piece, sizeof(piece), "^%ld", exp[j]);
                if (!dtest_append_text(buffer, capacity, length, piece))
                    return 0;
            }
            wrote_any_var = 1;
        }
    }

    *first_term = 0;
    return 1;
}

static int dtest_generate_random_large_prime_system_string(const slong *degrees,
                                                           slong npolys,
                                                           slong nvars,
                                                           double density_ratio,
                                                           ulong seed,
                                                           const fmpz_t prime,
                                                           char **polys_str_out)
{
    char *polys_str = NULL;
    size_t polys_cap = 0;
    size_t polys_len = 0;
    flint_rand_t rstate;

    if (!degrees || !polys_str_out || npolys <= 0 || nvars <= 0)
        return 0;

    flint_rand_init(rstate);
    flint_rand_set_seed(rstate, seed, seed + 1);

    for (slong i = 0; i < npolys; i++) {
        monomial_t *monomials = NULL;
        slong monomial_count = 0;
        slong degree = degrees[i] > 0 ? degrees[i] : 1;
        slong target_terms;
        slong *indices = NULL;
        slong leading_idx = -1;
        char *poly_buf = NULL;
        size_t poly_cap = 0;
        size_t poly_len = 0;
        int first_term = 1;
        slong selected = 0;
        fmpz_t coeff;

        fmpz_init(coeff);
        enumerate_all_monomials(&monomials, &monomial_count, nvars, degree);
        if (!monomials || monomial_count <= 0) {
            fmpz_clear(coeff);
            dtest_free_enumerated_monomials(monomials, monomial_count);
            flint_rand_clear(rstate);
            free(polys_str);
            return 0;
        }

        target_terms = (slong) (density_ratio * monomial_count);
        if (target_terms < 1) target_terms = 1;
        if (target_terms > monomial_count) target_terms = monomial_count;

        indices = (slong *) malloc((size_t) monomial_count * sizeof(slong));
        if (!indices) {
            fmpz_clear(coeff);
            dtest_free_enumerated_monomials(monomials, monomial_count);
            flint_rand_clear(rstate);
            free(polys_str);
            return 0;
        }

        for (slong j = 0; j < monomial_count; j++) indices[j] = j;
        for (slong j = monomial_count - 1; j > 0; j--) {
            slong k = n_randint(rstate, j + 1);
            slong tmp = indices[j];
            indices[j] = indices[k];
            indices[k] = tmp;
        }

        for (slong j = 0; j < monomial_count; j++) {
            if (monomials[indices[j]].total_degree == degree) {
                leading_idx = indices[j];
                break;
            }
        }
        if (leading_idx < 0) leading_idx = indices[0];

        do {
            fmpz_randm(coeff, rstate, prime);
        } while (fmpz_is_zero(coeff));
        if (!dtest_append_fmpz_monomial_text(&poly_buf, &poly_cap, &poly_len,
                                             coeff, monomials[leading_idx].exponents,
                                             nvars, &first_term)) {
            fmpz_clear(coeff);
            free(indices);
            free(poly_buf);
            dtest_free_enumerated_monomials(monomials, monomial_count);
            flint_rand_clear(rstate);
            free(polys_str);
            return 0;
        }
        selected++;

        for (slong j = 0; j < monomial_count && selected < target_terms; j++) {
            slong idx = indices[j];
            if (idx == leading_idx) continue;

            do {
                fmpz_randm(coeff, rstate, prime);
            } while (fmpz_is_zero(coeff));
            if (!dtest_append_fmpz_monomial_text(&poly_buf, &poly_cap, &poly_len,
                                                 coeff, monomials[idx].exponents,
                                                 nvars, &first_term)) {
                fmpz_clear(coeff);
                free(indices);
                free(poly_buf);
                dtest_free_enumerated_monomials(monomials, monomial_count);
                flint_rand_clear(rstate);
                free(polys_str);
                return 0;
            }
            selected++;
        }

        if (i > 0 && !dtest_append_text(&polys_str, &polys_cap, &polys_len, ", ")) {
            fmpz_clear(coeff);
            free(indices);
            free(poly_buf);
            dtest_free_enumerated_monomials(monomials, monomial_count);
            flint_rand_clear(rstate);
            free(polys_str);
            return 0;
        }
        if (!dtest_append_text(&polys_str, &polys_cap, &polys_len, poly_buf ? poly_buf : "0")) {
            fmpz_clear(coeff);
            free(indices);
            free(poly_buf);
            dtest_free_enumerated_monomials(monomials, monomial_count);
            flint_rand_clear(rstate);
            free(polys_str);
            return 0;
        }

        fmpz_clear(coeff);
        free(indices);
        free(poly_buf);
        dtest_free_enumerated_monomials(monomials, monomial_count);
    }

    flint_rand_clear(rstate);
    *polys_str_out = polys_str;
    return 1;
}

static int dtest_generate_random_system_string(const solver_field_case_t *field,
                                               const slong *degrees,
                                               slong npolys,
                                               slong nvars,
                                               double density_ratio,
                                               ulong seed,
                                               char **polys_str_out)
{
    fq_nmod_ctx_t ctx;
    fmpz_t p;
    char *polys_str = NULL;
    char **poly_strs = NULL;
    fq_mvpoly_t *polys = NULL;
    char *gen_name = NULL;
    size_t total_len = 4;
    flint_rand_t rstate;
    slong *deg_copy = NULL;
    int ok = 0;
    int saved_stdout = -1;

    if (!field || !degrees || !polys_str_out)
        return 0;

    if (field->large_prime_mode) {
        fmpz_init(p);
        if (fmpz_set_str(p, field->prime_str, 10) != 0) {
            fmpz_clear(p);
            return 0;
        }
        ok = dtest_generate_random_large_prime_system_string(degrees, npolys, nvars,
                                                             density_ratio, seed, p,
                                                             &polys_str);
        fmpz_clear(p);
        if (!ok)
            return 0;
        *polys_str_out = polys_str;
        return 1;
    }

    fmpz_init_set_ui(p, field->prime_ui);
    fq_nmod_ctx_init(ctx, p, field->power, field->generator ? field->generator : "t");
    fmpz_clear(p);

    deg_copy = (slong *) malloc((size_t) npolys * sizeof(slong));
    if (!deg_copy) {
        fq_nmod_ctx_clear(ctx);
        return 0;
    }
    for (slong i = 0; i < npolys; i++)
        deg_copy[i] = degrees[i];

    flint_rand_init(rstate);
    flint_rand_set_seed(rstate, seed, seed + 1);

    dtest_redirect_fd_to_devnull(STDOUT_FILENO, &saved_stdout);
    generate_polynomial_system(&polys, nvars, npolys, 0,
                               deg_copy, density_ratio, ctx, rstate);
    dtest_restore_fd(STDOUT_FILENO, saved_stdout);

    free(deg_copy);
    deg_copy = NULL;
    flint_rand_clear(rstate);

    if (!polys) {
        fq_nmod_ctx_clear(ctx);
        return 0;
    }

    poly_strs = (char **) malloc((size_t) npolys * sizeof(char *));
    if (!poly_strs) {
        for (slong i = 0; i < npolys; i++) fq_mvpoly_clear(&polys[i]);
        free(polys);
        fq_nmod_ctx_clear(ctx);
        return 0;
    }

    gen_name = get_generator_name(ctx);
    for (slong i = 0; i < npolys; i++) {
        char *s = fq_mvpoly_to_string(&polys[i], NULL, gen_name);
        poly_strs[i] = (s && strlen(s) > 0) ? s : (free(s), strdup("0"));
        if (!poly_strs[i]) {
            for (slong j = 0; j <= i; j++) free(poly_strs[j]);
            free(poly_strs);
            for (slong j = 0; j < npolys; j++) fq_mvpoly_clear(&polys[j]);
            free(polys);
            if (gen_name) free(gen_name);
            fq_nmod_ctx_clear(ctx);
            return 0;
        }
        total_len += strlen(poly_strs[i]) + 3;
    }

    polys_str = (char *) malloc(total_len);
    if (!polys_str) {
        for (slong i = 0; i < npolys; i++) free(poly_strs[i]);
        free(poly_strs);
        for (slong i = 0; i < npolys; i++) fq_mvpoly_clear(&polys[i]);
        free(polys);
        if (gen_name) free(gen_name);
        fq_nmod_ctx_clear(ctx);
        return 0;
    }
    polys_str[0] = '\0';
    for (slong i = 0; i < npolys; i++) {
        if (i > 0) strcat(polys_str, ", ");
        strcat(polys_str, poly_strs[i]);
    }

    for (slong i = 0; i < npolys; i++) free(poly_strs[i]);
    free(poly_strs);
    for (slong i = 0; i < npolys; i++) fq_mvpoly_clear(&polys[i]);
    free(polys);
    if (gen_name) free(gen_name);
    fq_nmod_ctx_clear(ctx);

    *polys_str_out = polys_str;
    return 1;
}

static char *dtest_compute_subres_resultant_str(const char *polys_str,
                                                const char *vars_str,
                                                const fq_nmod_ctx_t ctx)
{
    slong num_polys = 0, num_vars = 0;
    char **poly_array = split_string(polys_str, &num_polys);
    char **vars_array = split_string(vars_str, &num_vars);
    char *result = NULL;

    if (num_polys == 2 && num_vars == 1) {
        result = bivariate_resultant(poly_array[0], poly_array[1], vars_array[0], ctx);
    }

    free_split_strings(poly_array, num_polys);
    free_split_strings(vars_array, num_vars);
    return result;
}

static int dtest_verify_polynomial_solver_solutions(const char *polys_str,
                                                    polynomial_solutions_t *sols)
{
    (void) polys_str;

    if (!sols || !sols->is_valid)
        return 0;
    if (sols->has_no_solutions == 1)
        return 1;
    if (sols->has_no_solutions == -1)
        return 0;

    return sols->num_solution_sets > 0 &&
           sols->checked_solution_sets == sols->num_solution_sets &&
           sols->verified_solution_sets == sols->num_solution_sets;
}

static int dtest_verify_large_prime_solutions(const char *polys_str,
                                              large_prime_solutions_t *sols)
{
    (void) polys_str;

    if (!sols || !sols->is_valid)
        return 0;
    if (sols->has_no_solutions == 1)
        return 1;
    if (sols->has_no_solutions == -1)
        return 0;

    if (sols->num_solution_sets <= 0)
        return 0;

    if (!sols->solution_verified)
        return 0;

    for (slong set = 0; set < sols->num_solution_sets; set++) {
        if (!sols->solution_verified[set])
            return 0;
    }
    return 1;
}

// ============= DEGREE CHECKING FUNCTIONS =============

// Get maximum total degree of a polynomial
slong fq_mvpoly_max_total_degree(const fq_mvpoly_t *poly) {
    slong max_deg = 0;
    
    for (slong i = 0; i < poly->nterms; i++) {
        slong total_deg = 0;
        
        // Sum variable exponents
        if (poly->terms[i].var_exp) {
            for (slong j = 0; j < poly->nvars; j++) {
                total_deg += poly->terms[i].var_exp[j];
            }
        }
        
        // Sum parameter exponents
        if (poly->terms[i].par_exp) {
            for (slong j = 0; j < poly->npars; j++) {
                total_deg += poly->terms[i].par_exp[j];
            }
        }
        
        if (total_deg > max_deg) {
            max_deg = total_deg;
        }
    }
    
    return max_deg;
}

// Add degree check result to global storage
void add_degree_check_result(const char *test_name, slong nvars, slong npars, slong npolys,
                             const slong *input_degrees, slong resultant_degree, 
                             int algorithm_type, int test_case_index) {
    // Ensure capacity
    if (num_degree_results >= degree_results_capacity) {
        degree_results_capacity = degree_results_capacity == 0 ? 10 : degree_results_capacity * 2;
        degree_results = (degree_check_result_t*) realloc(degree_results, 
                                                          degree_results_capacity * sizeof(degree_check_result_t));
    }
    
    degree_check_result_t *result = &degree_results[num_degree_results];
    
    // Copy test name
    strncpy(result->test_name, test_name, 255);
    result->test_name[255] = '\0';
    
    result->nvars = nvars;
    result->npars = npars;
    result->npolys = npolys;
    result->algorithm_type = algorithm_type;
    result->test_case_index = test_case_index;
    
    // Copy input degrees
    result->input_degrees = (slong*) malloc(npolys * sizeof(slong));
    memcpy(result->input_degrees, input_degrees, npolys * sizeof(slong));
    
    result->resultant_degree = resultant_degree;
    
    // Calculate expected degree product
    result->expected_degree_product = 1;
    for (slong i = 0; i < npolys; i++) {
        result->expected_degree_product *= input_degrees[i];
    }
    
    // Check if degree bound is satisfied
    result->degree_satisfies_bound = (resultant_degree <= result->expected_degree_product);
    
    num_degree_results++;
}

// Print comprehensive summary of all degree checks with algorithm comparison
void print_degree_check_summary(void) {
    printf("\n");
    printf("╔════════════════════════════════════════════════════════════════════════════╗\n");
    printf("║         DIXON RESULTANT DEGREE ANALYSIS - ALGORITHM COMPARISON             ║\n");
    printf("╚════════════════════════════════════════════════════════════════════════════╝\n");
    printf("\n");
    
    if (num_degree_results == 0) {
        printf("No degree check results collected.\n");
        return;
    }
    
    // Find number of unique test cases
    int max_test_case_index = -1;
    for (slong i = 0; i < num_degree_results; i++) {
        if (degree_results[i].test_case_index > max_test_case_index) {
            max_test_case_index = degree_results[i].test_case_index;
        }
    }
    int num_test_cases = max_test_case_index + 1;
    
    const char *algorithm_names[2] = {
        "Traditional Max-Rank Submatrix Extraction",
        "Degree-Optimal Submatrix Selection"
    };
    
    // Print results for each test case
    for (int tc = 0; tc < num_test_cases; tc++) {
        printf("\n");
        printf("═════════════════════════════════════════════════════════════════════════════\n");
        printf("TEST CASE %d\n", tc + 1);
        
        // Get test case name from first result of this test case
        for (slong i = 0; i < num_degree_results; i++) {
            if (degree_results[i].test_case_index == tc) {
                printf("Name: %s\n", degree_results[i].test_name);
                printf("System: %ld vars, %ld params, %ld polynomials\n", 
                       degree_results[i].nvars, degree_results[i].npars, degree_results[i].npolys);
                printf("Degrees: [");
                for (slong j = 0; j < degree_results[i].npolys; j++) {
                    printf("%ld", degree_results[i].input_degrees[j]);
                    if (j < degree_results[i].npolys - 1) printf(", ");
                }
                printf("]\n");
                printf("Expected bound: %ld\n", degree_results[i].expected_degree_product);
                break;
            }
        }
        printf("═════════════════════════════════════════════════════════════════════════════\n");
        
        // For each algorithm, calculate statistics for this test case
        for (int algo = 0; algo <= 1; algo++) {
            printf("\n");
            printf("─────────────────────────────────────────────────────────────────────────────\n");
            printf("Algorithm: %s\n", algorithm_names[algo]);
            printf("─────────────────────────────────────────────────────────────────────────────\n");
            
            slong count = 0;
            slong satisfied = 0;
            slong violated = 0;
            slong min_deg = LONG_MAX;
            slong max_deg = 0;
            double avg_deg = 0.0;
            double avg_ratio = 0.0;
            
            // Collect statistics
            for (slong i = 0; i < num_degree_results; i++) {
                if (degree_results[i].test_case_index != tc || 
                    degree_results[i].algorithm_type != algo) {
                    continue;
                }
                
                count++;
                if (degree_results[i].degree_satisfies_bound) {
                    satisfied++;
                } else {
                    violated++;
                }
                
                slong deg = degree_results[i].resultant_degree;
                if (deg < min_deg) min_deg = deg;
                if (deg > max_deg) max_deg = deg;
                avg_deg += deg;
                
                if (degree_results[i].expected_degree_product > 0) {
                    avg_ratio += (double)deg / degree_results[i].expected_degree_product;
                }
            }
            
            if (count > 0) {
                avg_deg /= count;
                avg_ratio /= count;
            }
            
            printf("  Repetitions: %ld\n", count);
            printf("  Bound satisfied: %ld (%.1f%%)\n", 
                   satisfied, count > 0 ? 100.0 * satisfied / count : 0.0);
            printf("  Bound violated: %ld (%.1f%%)\n", 
                   violated, count > 0 ? 100.0 * violated / count : 0.0);
            printf("\n");
            printf("  Resultant Degree Statistics:\n");
            printf("    Min: %ld\n", min_deg == LONG_MAX ? 0 : min_deg);
            printf("    Max: %ld\n", max_deg);
            printf("    Avg: %.2f\n", avg_deg);
            printf("    Avg ratio (deg/bound): %.3f\n", avg_ratio);
        }
        
        // Comparison for this test case
        printf("\n");
        printf("─────────────────────────────────────────────────────────────────────────────\n");
        printf("Comparison Summary for Test Case %d:\n", tc + 1);
        printf("─────────────────────────────────────────────────────────────────────────────\n");
        
        slong trad_count = 0, opt_count = 0;
        slong trad_satisfied = 0, opt_satisfied = 0;
        double trad_avg_ratio = 0.0, opt_avg_ratio = 0.0;
        
        for (slong i = 0; i < num_degree_results; i++) {
            if (degree_results[i].test_case_index != tc) continue;
            
            if (degree_results[i].algorithm_type == 0) {
                trad_count++;
                if (degree_results[i].degree_satisfies_bound) trad_satisfied++;
                if (degree_results[i].expected_degree_product > 0) {
                    trad_avg_ratio += (double)degree_results[i].resultant_degree / 
                                     degree_results[i].expected_degree_product;
                }
            } else {
                opt_count++;
                if (degree_results[i].degree_satisfies_bound) opt_satisfied++;
                if (degree_results[i].expected_degree_product > 0) {
                    opt_avg_ratio += (double)degree_results[i].resultant_degree / 
                                    degree_results[i].expected_degree_product;
                }
            }
        }
        
        if (trad_count > 0) trad_avg_ratio /= trad_count;
        if (opt_count > 0) opt_avg_ratio /= opt_count;
        
        printf("\n");
        printf("┌─────────────────────────────────────┬──────────────┬──────────────┐\n");
        printf("│ Metric                              │ Traditional  │ Deg-Optimal  │\n");
        printf("├─────────────────────────────────────┼──────────────┼──────────────┤\n");
        printf("│ Repetitions                         │ %12ld │ %12ld │\n", 
               trad_count, opt_count);
        printf("│ Bound satisfied                     │ %12ld │ %12ld │\n", 
               trad_satisfied, opt_satisfied);
        printf("│ Success rate (%%)                    │ %11.1f%% │ %11.1f%% │\n",
               trad_count > 0 ? 100.0 * trad_satisfied / trad_count : 0.0,
               opt_count > 0 ? 100.0 * opt_satisfied / opt_count : 0.0);
        printf("│ Avg ratio (deg/bound)               │ %12.3f │ %12.3f │\n",
               trad_avg_ratio, opt_avg_ratio);
        printf("└─────────────────────────────────────┴──────────────┴──────────────┘\n");
    }
    
    // Print overall summary
    printf("\n");
    printf("═════════════════════════════════════════════════════════════════════════════\n");
    printf("OVERALL COMPARATIVE ANALYSIS (All Test Cases)\n");
    printf("═════════════════════════════════════════════════════════════════════════════\n");
    
    slong traditional_count = 0, optimal_count = 0;
    slong traditional_satisfied = 0, optimal_satisfied = 0;
    slong traditional_violated = 0, optimal_violated = 0;
    
    for (slong i = 0; i < num_degree_results; i++) {
        if (degree_results[i].algorithm_type == 0) {
            traditional_count++;
            if (degree_results[i].degree_satisfies_bound) {
                traditional_satisfied++;
            } else {
                traditional_violated++;
            }
        } else {
            optimal_count++;
            if (degree_results[i].degree_satisfies_bound) {
                optimal_satisfied++;
            } else {
                optimal_violated++;
            }
        }
    }
    
    printf("\n");
    printf("┌─────────────────────────────────────┬──────────────┬──────────────┐\n");
    printf("│ Metric                              │ Traditional  │ Deg-Optimal  │\n");
    printf("├─────────────────────────────────────┼──────────────┼──────────────┤\n");
    printf("│ Total tests                         │ %12ld │ %12ld │\n", 
           traditional_count, optimal_count);
    printf("│ Bound satisfied                     │ %12ld │ %12ld │\n", 
           traditional_satisfied, optimal_satisfied);
    printf("│ Bound violated                      │ %12ld │ %12ld │\n", 
           traditional_violated, optimal_violated);
    printf("│ Success rate (%%)                    │ %11.1f%% │ %11.1f%% │\n",
           traditional_count > 0 ? 100.0 * traditional_satisfied / traditional_count : 0.0,
           optimal_count > 0 ? 100.0 * optimal_satisfied / optimal_count : 0.0);
    printf("└─────────────────────────────────────┴──────────────┴──────────────┘\n");
    
    printf("\n");
    printf("═════════════════════════════════════════════════════════════════════════════\n");
}

// Cleanup degree check results
void cleanup_degree_check_results(void) {
    if (degree_results) {
        for (slong i = 0; i < num_degree_results; i++) {
            if (degree_results[i].input_degrees) {
                free(degree_results[i].input_degrees);
            }
        }
        free(degree_results);
        degree_results = NULL;
    }
    num_degree_results = 0;
    degree_results_capacity = 0;
}


// ============= Math Utility Functions =============

// Calculate binomial coefficient C(n, k)
slong binomial_coefficient(slong n, slong k) {
    if (k > n || k < 0) return 0;
    if (k == 0 || k == n) return 1;
    
    // Use symmetry to optimize calculation
    if (k > n - k) k = n - k;
    
    slong result = 1;
    for (slong i = 0; i < k; i++) {
        result = result * (n - i) / (i + 1);
    }
    return result;
}

// Calculate total number of possible monomials for given variables, parameters and degree
slong count_possible_monomials(slong nvars, slong npars, slong max_degree) {
    slong total_indeterminates = nvars + npars;
    return binomial_coefficient(total_indeterminates + max_degree, max_degree);
}


static void enumerate_all_monomials_recursive(monomial_t *monomials,
                                              slong *count,
                                              slong total_indeterminates,
                                              slong *current_exp,
                                              slong pos,
                                              slong remaining_degree) {
    if (pos == total_indeterminates) {
        if (remaining_degree == 0) {
            monomials[*count].exponents = (slong*) malloc(total_indeterminates * sizeof(slong));
            memcpy(monomials[*count].exponents, current_exp, total_indeterminates * sizeof(slong));

            slong total_deg = 0;
            for (slong i = 0; i < total_indeterminates; i++) {
                total_deg += current_exp[i];
            }
            monomials[*count].total_degree = total_deg;
            (*count)++;
        }
        return;
    }

    for (slong deg = 0; deg <= remaining_degree; deg++) {
        current_exp[pos] = deg;
        enumerate_all_monomials_recursive(monomials, count, total_indeterminates,
                                          current_exp, pos + 1, remaining_degree - deg);
    }
}

// ============= Polynomial Generation =============
// Generate all possible monomials with degree <= max_degree
void enumerate_all_monomials(monomial_t **monomials, slong *count, 
                            slong total_indeterminates, slong max_degree) {
    // Calculate maximum possible monomials
    slong max_possible = 1;
    for (slong d = 0; d <= max_degree; d++) {
        max_possible += binomial_coefficient(total_indeterminates + d - 1, d);
    }
    
    *monomials = (monomial_t*) malloc(max_possible * sizeof(monomial_t));
    *count = 0;
    
    // Generate monomials for each total degree
    slong *temp_exp = (slong*) calloc(total_indeterminates, sizeof(slong));
    for (slong degree = 0; degree <= max_degree; degree++) {
        enumerate_all_monomials_recursive(*monomials, count, total_indeterminates,
                                          temp_exp, 0, degree);
    }
    free(temp_exp);
}

// Improved polynomial generation with better density control
void generate_random_polynomial(fq_mvpoly_t *poly, slong nvars, slong npars,
                               slong max_degree, double density_ratio,
                               const fq_nmod_ctx_t ctx, flint_rand_t state) {
    
    fq_mvpoly_init(poly, nvars, npars, ctx);
    
    slong total_indeterminates = nvars + npars;
    
    // Generate all possible monomials
    monomial_t *all_monomials;
    slong total_monomials;
    enumerate_all_monomials(&all_monomials, &total_monomials, total_indeterminates, max_degree);
    
    if (!g_bezout_test_quiet_mode) {
        printf("    Total possible monomials: %ld\n", total_monomials);
    }
    
    // Calculate target number of terms
    slong target_terms = (slong)(density_ratio * total_monomials);
    if (target_terms < 1) target_terms = 1;
    if (target_terms > total_monomials) target_terms = total_monomials;
    
    if (!g_bezout_test_quiet_mode) {
        printf("    Target terms: %ld (%.1f%% density)\n", target_terms,
               (double) target_terms / total_monomials * 100);
    }
    
    // Shuffle the monomial array to randomize selection
    for (slong i = total_monomials - 1; i > 0; i--) {
        slong j = n_randint(state, i + 1);
        // Swap monomials[i] and monomials[j]
        monomial_t temp = all_monomials[i];
        all_monomials[i] = all_monomials[j];
        all_monomials[j] = temp;
    }
    
    // Select first target_terms monomials and add them to polynomial
    for (slong i = 0; i < target_terms; i++) {
        monomial_t *mon = &all_monomials[i];
        
        // Split exponents into variable and parameter parts
        slong *var_exp = NULL;
        slong *par_exp = NULL;
        
        // Extract variable exponents (first nvars positions)
        if (nvars > 0) {
            int has_var = 0;
            for (slong j = 0; j < nvars; j++) {
                if (mon->exponents[j] > 0) {
                    has_var = 1;
                    break;
                }
            }
            if (has_var || (nvars > 0 && npars == 0)) {  // Include even if all zeros for pure variable case
                var_exp = (slong*) malloc(nvars * sizeof(slong));
                memcpy(var_exp, mon->exponents, nvars * sizeof(slong));
            }
        }
        
        // Extract parameter exponents (last npars positions)
        if (npars > 0) {
            int has_par = 0;
            for (slong j = nvars; j < total_indeterminates; j++) {
                if (mon->exponents[j] > 0) {
                    has_par = 1;
                    break;
                }
            }
            if (has_par || (npars > 0 && nvars == 0)) {  // Include even if all zeros for pure parameter case
                par_exp = (slong*) malloc(npars * sizeof(slong));
                memcpy(par_exp, mon->exponents + nvars, npars * sizeof(slong));
            }
        }
        
        // Generate random non-zero coefficient
        fq_nmod_t coeff;
        fq_nmod_init(coeff, ctx);
        do {
            fq_nmod_randtest(coeff, state, ctx);
        } while (fq_nmod_is_zero(coeff, ctx));
        
        // Add the term to polynomial
        fq_mvpoly_add_term(poly, var_exp, par_exp, coeff);
        
        // Cleanup
        fq_nmod_clear(coeff, ctx);
        if (var_exp) free(var_exp);
        if (par_exp) free(par_exp);
    }
    
    // Cleanup monomial list
    for (slong i = 0; i < total_monomials; i++) {
        free(all_monomials[i].exponents);
    }
    free(all_monomials);
    
    if (!g_bezout_test_quiet_mode) {
        printf("    Generated polynomial with %ld terms (target: %ld, achieved density: %.1f%%)\n",
               poly->nterms, target_terms,
               (double) poly->nterms / total_monomials * 100);
    }
}

// Generate polynomial system with specified degrees and density
void generate_polynomial_system(fq_mvpoly_t **polys, slong nvars, slong npolys, 
                               slong npars, const slong *degrees,
                               double density_ratio,
                               const fq_nmod_ctx_t ctx, flint_rand_t state) {
    
    if (degrees == NULL) {
        printf("Error: degrees array cannot be NULL\n");
        return;
    }
    
    *polys = (fq_mvpoly_t*) malloc(npolys * sizeof(fq_mvpoly_t));
    if (!*polys) {
        printf("[ERROR] Failed to allocate polynomial array\n");
        return;
    }
    
    for (slong i = 0; i < npolys; i++) {
        
        slong max_degree = degrees[i];
        if (max_degree <= 0) {
            max_degree = 2;  // Default value
        }
        
        // First few polynomials include parameters
        slong poly_npars = npars;
        
        generate_random_polynomial(&(*polys)[i], nvars, poly_npars, 
                                  max_degree, density_ratio, ctx, state);
        
    }
    
}

// ============= Test Functions =============

void test_dixon_system(const char *test_name, slong nvars, slong npars,
                      ulong p, slong field_degree, const slong *degrees,
                      slong npolys, double density_ratio, flint_rand_t state,
                      int test_case_index) {
    printf("\n=== %s ===\n", test_name);
    printf("Algorithm: %s\n", g_matrix_transpose_threshold == 0 ? 
           "Traditional Max-Rank" : "Degree-Optimal Selection");
    printf("Field: GF(%lu^%ld), Variables: %ld, Parameters: %ld, Density: %.1f%%\n", 
           p, field_degree, nvars, npars, density_ratio * 100);
    
    // Print degree information
    printf("Polynomial degrees: [");
    for (slong i = 0; i < npolys; i++) {
        printf("%ld", degrees[i]);
        if (i < npolys - 1) printf(", ");
    }
    printf("]\n");
    
    // Calculate theoretical complexity
    slong total_theoretical_terms = 0;
    for (slong i = 0; i < npolys; i++) {
        slong possible = count_possible_monomials(nvars, npars, degrees[i]);
        printf("Polynomial %ld: max possible terms = %ld\n", i, possible);
        total_theoretical_terms += possible;
    }
    printf("System theoretical total terms: %ld\n", total_theoretical_terms);
    
    // Initialize field
    fq_nmod_ctx_t ctx;
    fmpz_t p_fmpz;
    fmpz_init(p_fmpz);
    fmpz_set_ui(p_fmpz, p);
    fq_nmod_ctx_init(ctx, p_fmpz, field_degree, "t");
    fmpz_clear(p_fmpz);
    
    // Generate polynomial system
    fq_mvpoly_t *polys;
    generate_polynomial_system(&polys, nvars, npolys, npars, 
                              degrees, density_ratio, ctx, state);
    
    // Print system information
    printf("\nPolynomial system:\n");
    slong actual_total_terms = 0;
    for (slong i = 0; i < npolys; i++) {
        printf("  p%ld (degree %ld): %ld terms", i, degrees[i], polys[i].nterms);
        if (polys[i].npars > 0) {
            printf(" (with %ld parameters)", polys[i].npars);
        }
        printf("\n");
        actual_total_terms += polys[i].nterms;
    }
    printf("System actual total terms: %ld (density: %.1f%%)\n", 
           actual_total_terms, (double)actual_total_terms / total_theoretical_terms * 100);
    
    // Compute Dixon resultant
    fq_mvpoly_t result;
    
    clock_t start = clock();
    fq_dixon_resultant(&result, polys, nvars, npars);
    clock_t end = clock();
    
    
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Computation time: %.3f seconds\n", elapsed);
    printf("Resultant: %ld terms", result.nterms);
    if (result.npars > 0) {
        printf(" with %ld parameters", result.npars);
    }
    printf("\n");

    // ========== DEGREE CHECKING CODE ==========
    
    // Calculate resultant degree
    slong resultant_degree = fq_mvpoly_max_total_degree(&result);
    printf("Resultant total degree: %ld\n", resultant_degree);
    
    // Calculate expected degree bound (product of input degrees)
    slong expected_product = 1;
    for (slong i = 0; i < npolys; i++) {
        expected_product *= degrees[i];
    }
    printf("Expected degree bound (product): %ld\n", expected_product);
    
    // Check if bound is satisfied
    if (resultant_degree <= expected_product) {
        printf("✓ Degree check PASSED: deg(resultant) ≤ ∏deg(inputs)\n");
    } else {
        printf("✗ Degree check FAILED: deg(resultant) > ∏deg(inputs)\n");
    }
    
    // Store result for summary (with algorithm type and test case index)
    int algo_type = (g_matrix_transpose_threshold == 0) ? 0 : 1;
    add_degree_check_result(test_name, nvars, npars, npolys, degrees, 
                           resultant_degree, algo_type, test_case_index);
    
    // ========== END DEGREE CHECKING CODE ==========
    
    if (result.nterms > 0 && result.nterms <= 10) {
        printf("  ");
        fq_mvpoly_print(&result, "R");
        printf("\n");
    }
    
    // Cleanup
    fq_mvpoly_clear(&result);
    for (slong i = 0; i < npolys; i++) {
        fq_mvpoly_clear(&polys[i]);
    }
    free(polys);
    fq_nmod_ctx_clear(ctx);
}

static void test_dixon_system_quiet(const char *test_name, slong nvars, slong npars,
                                    ulong p, slong field_degree, const slong *degrees,
                                    slong npolys, double density_ratio, flint_rand_t state,
                                    int test_case_index) {
    fq_nmod_ctx_t ctx;
    fmpz_t p_fmpz;
    fq_mvpoly_t *polys;
    fq_mvpoly_t result;
    slong resultant_degree;
    int algo_type;
    int saved_verbose_level = g_dixon_verbose_level;
    int saved_debug_mode = g_dixon_debug_mode;

    fmpz_init(p_fmpz);
    fmpz_set_ui(p_fmpz, p);
    fq_nmod_ctx_init(ctx, p_fmpz, field_degree, "t");
    fmpz_clear(p_fmpz);

    generate_polynomial_system(&polys, nvars, npolys, npars,
                               degrees, density_ratio, ctx, state);

    g_dixon_verbose_level = 0;
    g_dixon_debug_mode = 0;
    fq_dixon_resultant(&result, polys, nvars, npars);
    g_dixon_verbose_level = saved_verbose_level;
    g_dixon_debug_mode = saved_debug_mode;

    resultant_degree = fq_mvpoly_max_total_degree(&result);
    algo_type = (g_matrix_transpose_threshold == 0) ? 0 : 1;
    add_degree_check_result(test_name, nvars, npars, npolys, degrees,
                            resultant_degree, algo_type, test_case_index);

    fq_mvpoly_clear(&result);
    for (slong i = 0; i < npolys; i++) {
        fq_mvpoly_clear(&polys[i]);
    }
    free(polys);
    fq_nmod_ctx_clear(ctx);
}

void test_xhash() {
    printf("=== Prime Field Polynomial System Dixon Resultant Implementation ===\n\n");    
    const char *ideal = "x7^3 = -19561*x4^3 + 2061*x4^2*x5 + 25073*x4^2*x6 + 19787*x4^2 + 779*x4*x5^2 + 31516*x4*x5*x6 - 17049*x4*x5 + 9065*x4*x6^2 - 14413*x4*x6 - 9964*x4 + 24021*x5^3 - 31858*x5^2*x6 - 20667*x5^2 - 19224*x5*x6^2 - 15430*x5*x6 + 19731*x5 + 31617*x6^3 + 6653*x6^2 - 13781*x6 - 15782,\
    x8^3 = 31617*x4^3 + 9065*x4^2*x5 - 19224*x4^2*x6 + 6653*x4^2 + 25073*x4*x5^2 + 31516*x4*x5*x6 - 14413*x4*x5 - 31858*x4*x6^2 - 15430*x4*x6 - 13781*x4 - 19561*x5^3 + 2061*x5^2*x6 + 19787*x5^2 + 779*x5*x6^2 - 17049*x5*x6 - 9964*x5 + 24021*x6^3 - 20667*x6^2 + 19731*x6 - 15782,\
    x9^3 = 24021*x4^3 - 31858*x4^2*x5 + 779*x4^2*x6 - 20667*x4^2 - 19224*x4*x5^2 + 31516*x4*x5*x6 - 15430*x4*x5 + 2061*x4*x6^2 - 17049*x4*x6 + 19731*x4 + 31617*x5^3 + 9065*x5^2*x6 + 6653*x5^2 + 25073*x5*x6^2 - 14413*x5*x6 - 13781*x5 - 19561*x6^3 + 19787*x6^2 - 9964*x6 - 15782";
    
    const char* f1 = "6653*x0^2 - 13781*x0 - x2^2 - 15782";
    const char* f2 = "20667*x0^2 + 19731*x0 - x3^2 - 15782";
    const char* f3 = "3*x1^3 + 16*x1^2*x2 - 24563*x1^2 - 15*x1*x2^2 + 3*x1*x2*x3 + 8202*x1*x2 + 16*x1*x3^2 - 8170*x1*x3 - 6753*x1 - 16*x2^3 + 8*x2^2*x3 - 24637*x2^2 - 9*x2*x3^2 + 27861*x2*x3 - 19866*x2 - 13*x3^3 + 26199*x3^2 - 26963*x3 - x4 - 32551";
    const char* f4 = "-4*x1^3 - 4*x1^2*x2 - 15*x1^2*x3 - 11483*x1^2 + 15*x1*x2^2 - 12*x1*x2*x3 + 16410*x1*x2 - x1*x3^2 + 19631*x1*x3 + 28842*x1 + 11*x2^3 + 6*x2^2*x3 - 3245*x2^2 - 4*x2*x3^2 - 26213*x2*x3 - 5476*x2 + 12*x3^3 + 11475*x3^2 + 12887*x3 - x5 - 8584";
    const char* f5 = "16*x1^3 - 11*x1^2*x2 + 15*x1^2*x3 + 6597*x1^2 + 4*x1*x2^2 - 15*x1*x2*x3 - 21332*x1*x2 + 14*x1*x3^2 + 18051*x1*x3 - 22387*x1 + 2*x2^3 - 5*x2^2*x3 + 8180*x2^2 - 8*x2*x3^2 - 27885*x2*x3 - 15205*x2 + 15*x3^3 - 21257*x3^2 - 3092*x3 - x6 - 25478";
    const char* f6 = "-19561*x4^3 + 2061*x4^2*x5 + 25073*x4^2*x6 + 19787*x4^2 + 779*x4*x5^2 + 31516*x4*x5*x6 - 17049*x4*x5 + 9065*x4*x6^2 - 14413*x4*x6 - 9964*x4 + 24021*x5^3 - 31858*x5^2*x6 - 20667*x5^2 - 19224*x5*x6^2 - 15430*x5*x6 + 19731*x5 + 31617*x6^3 + 6653*x6^2 - 13781*x6 - x7^3 - 15782";
    const char* f7 = "31617*x4^3 + 9065*x4^2*x5 - 19224*x4^2*x6 + 6653*x4^2 + 25073*x4*x5^2 + 31516*x4*x5*x6 - 14413*x4*x5 - 31858*x4*x6^2 - 15430*x4*x6 - 13781*x4 - 19561*x5^3 + 2061*x5^2*x6 + 19787*x5^2 + 779*x5*x6^2 - 17049*x5*x6 - 9964*x5 + 24021*x6^3 - 20667*x6^2 + 19731*x6 - x8^3 - 15782";
    const char* f8 = "24021*x4^3 - 31858*x4^2*x5 + 779*x4^2*x6 - 20667*x4^2 - 19224*x4*x5^2 + 31516*x4*x5*x6 - 15430*x4*x5 + 2061*x4*x6^2 - 17049*x4*x6 + 19731*x4 + 31617*x5^3 + 9065*x5^2*x6 + 6653*x5^2 + 25073*x5*x6^2 - 14413*x5*x6 - 13781*x5 - 19561*x6^3 + 19787*x6^2 - 9964*x6 - x9^3 - 15782";
    const char* f9 = "3*x7^3 + 16*x7^2*x8 - 24563*x7^2 - 15*x7*x8^2 + 3*x7*x8*x9 + 8202*x7*x8 + 16*x7*x9^2 - 8170*x7*x9 - 6753*x7 - 16*x8^3 + 8*x8^2*x9 - 24637*x8^2 - 9*x8*x9^2 + 27861*x8*x9 - 19866*x8 - 13*x9^3 + 26199*x9^2 - 26963*x9 - 32551";

    fq_nmod_ctx_t ctx;
    mp_limb_t prime = 65537;
    fmpz_t p;
    fmpz_init_set_ui(p, prime);
    fq_nmod_ctx_init(ctx, p, 1, "t");    
    
    printf("Prime field: GF(p) where p = 65537\n");
    printf("\n");
    
    clock_t total_start = clock();
    
    char* r1 = RESULTANT((f1, f2), ("x0"));
    char* r2 = DIXON((r1, f3, f4, f5), ("x1", "x2", "x3"));
    char* r3 = DIXON_WITH_IDEAL((f9, f8), ("x9"));
    char* r4 = DIXON_WITH_IDEAL((r3, f7), ("x8"));
    char* r5 = DIXON_WITH_IDEAL((r4, f6), ("x7"));
    char* r6 = DIXON_COMPLEXITY((r2, r5), ("x4"));

    printf("%s",r6);
    clock_t total_end = clock();
    double total_time = (double)(total_end - total_start) / CLOCKS_PER_SEC;
    
    printf("Total computation time: %.3f seconds\n", total_time);
    
    free(r1);
    free(r2);
    free(r3);
    free(r4);
    free(r5);
    free(r6);
    
    fq_nmod_ctx_clear(ctx);
    printf("\n=== Computation Complete ===\n");
} 

int test_polynomial_solver_correctness(void)
{
    static const solver_degree_case_t degree_cases[] = {
        {"[5]*2", 2, 2, {5, 5, 0, 0}},
        {"[3]*3", 3, 3, {3, 3, 3, 0}},
        {"[2]*4", 4, 4, {2, 2, 2, 2}},
        {"[2,3,4]", 3, 3, {2, 3, 4, 0}}
    };
    static const solver_field_case_t field_cases[] = {
        {"17", 0, 0, 17, 1, NULL, "t", NULL},
        {"65537", 0, 0, 65537, 1, NULL, "t", NULL},
        {"2^8", 0, 0, 2, 8, NULL, "t", NULL},
        {"18446744073709551629", 0, 1, 0, 1, NULL, "t", "18446744073709551629"}
    };
    const slong repetitions = 8;
    const double density_ratio = 1.0;
    solver_correctness_stats_t stats[4][4];
    slong total_cases = 0;
    slong total_passed = 0;
    slong total_no_solution = 0;
    slong total_failed = 0;

    memset(stats, 0, sizeof(stats));

    printf("=== Polynomial Solver Correctness Test ===\n\n");
    printf("Random systems per combination: %ld\n", repetitions);
    printf("Degree sets: 4\n");
    printf("Fields: 4\n");
    printf("Total cases: %ld\n\n", (slong) (4 * 4 * repetitions));

    for (slong d = 0; d < 4; d++) {
        const solver_degree_case_t *deg_case = &degree_cases[d];
        for (slong f = 0; f < 4; f++) {
            const solver_field_case_t *field = &field_cases[f];

            printf("[%ld/%ld] %s over %s\n",
                   (slong) (d * 4 + f + 1), (slong) 16,
                   deg_case->label, field->label);

            for (slong rep = 0; rep < repetitions; rep++) {
                char *polys_str = NULL;
                ulong seed = 1000003UL + (ulong) (d * 1000 + f * 100 + rep);

                printf("  Progress %ld/%ld\r", total_cases + 1, (slong) (4 * 4 * repetitions));
                fflush(stdout);

                total_cases++;
                if (!dtest_generate_random_system_string(field,
                                                        deg_case->degrees,
                                                        deg_case->npolys,
                                                        deg_case->nvars,
                                                        density_ratio,
                                                        seed,
                                                        &polys_str)) {
                    stats[d][f].failed_cases++;
                    total_failed++;
                    continue;
                }

                if (field->large_prime_mode) {
                    large_prime_solutions_t *sols;
                    fmpz_t prime;
                    int ok;
                    int saved_stdout = -1;
                    int saved_stderr = -1;

                    fmpz_init(prime);
                    if (fmpz_set_str(prime, field->prime_str, 10) != 0) {
                        fmpz_clear(prime);
                        free(polys_str);
                        stats[d][f].failed_cases++;
                        total_failed++;
                        continue;
                    }

                    large_prime_solver_set_realtime_progress(0);
                    dtest_redirect_fd_to_devnull(STDOUT_FILENO, &saved_stdout);
                    dtest_redirect_fd_to_devnull(STDERR_FILENO, &saved_stderr);
                    sols = solve_large_prime_polynomial_system_string(polys_str, prime);
                    dtest_restore_fd(STDOUT_FILENO, saved_stdout);
                    dtest_restore_fd(STDERR_FILENO, saved_stderr);
                    fmpz_clear(prime);

                    ok = dtest_verify_large_prime_solutions(polys_str, sols);
                    if (!sols || !sols->is_valid || sols->has_no_solutions == -1 || !ok) {
                        stats[d][f].failed_cases++;
                        total_failed++;
                    } else if (sols->has_no_solutions || sols->num_solution_sets == 0) {
                        stats[d][f].no_solution_cases++;
                        total_no_solution++;
                    } else {
                        stats[d][f].passed_cases++;
                        total_passed++;
                    }

                    if (sols) {
                        large_prime_solutions_clear(sols);
                        free(sols);
                    }
                } else {
                    polynomial_solutions_t *sols;
                    fq_nmod_ctx_t ctx;
                    fmpz_t p;
                    int ok;
                    int saved_stdout = -1;
                    int saved_stderr = -1;

                    fmpz_init_set_ui(p, field->prime_ui);
                    fq_nmod_ctx_init(ctx, p, field->power, field->generator ? field->generator : "t");
                    fmpz_clear(p);

                    polynomial_solver_set_realtime_progress(0);
                    polynomial_solver_set_internal_trace(0);
                    dtest_redirect_fd_to_devnull(STDOUT_FILENO, &saved_stdout);
                    dtest_redirect_fd_to_devnull(STDERR_FILENO, &saved_stderr);
                    sols = solve_polynomial_system_string(polys_str, ctx);
                    dtest_restore_fd(STDOUT_FILENO, saved_stdout);
                    dtest_restore_fd(STDERR_FILENO, saved_stderr);

                    ok = dtest_verify_polynomial_solver_solutions(polys_str, sols);
                    if (!sols || !sols->is_valid || sols->has_no_solutions == -1 || !ok) {
                        stats[d][f].failed_cases++;
                        total_failed++;
                    } else if (sols->has_no_solutions || sols->num_solution_sets == 0) {
                        stats[d][f].no_solution_cases++;
                        total_no_solution++;
                    } else {
                        stats[d][f].passed_cases++;
                        total_passed++;
                    }

                    if (sols) {
                        polynomial_solutions_clear(sols);
                        free(sols);
                    }
                    fq_nmod_ctx_clear(ctx);
                }

                free(polys_str);
            }

            printf("  Result: passed=%ld, failed=%ld, no-solution=%ld\n",
                   stats[d][f].passed_cases,
                   stats[d][f].failed_cases,
                   stats[d][f].no_solution_cases);
        }
    }

    printf("Summary by degree/field combination:\n");
    for (slong d = 0; d < 4; d++) {
        for (slong f = 0; f < 4; f++) {
            printf("  %s over %s: passed=%ld, failed=%ld, no-solution=%ld\n",
                   degree_cases[d].label,
                   field_cases[f].label,
                   stats[d][f].passed_cases,
                   stats[d][f].failed_cases,
                   stats[d][f].no_solution_cases);
        }
    }

    printf("\nOverall: passed=%ld, failed=%ld, no-solution=%ld, total=%ld\n",
           total_passed, total_failed, total_no_solution, total_cases);

    return (total_failed == 0) ? 0 : 1;
}

int test_polynomial_solver_performance(void)
{
    static const solver_degree_case_t perf_cases[] = {
        {"2 vars deg 8", 2, 2, {8, 8, 0, 0, 0, 0, 0, 0}},
        {"2 vars deg 16", 2, 2, {16, 16, 0, 0, 0, 0, 0, 0}},
        {"2 vars deg 32", 2, 2, {32, 32, 0, 0, 0, 0, 0, 0}},
        {"2 vars deg 64", 2, 2, {64, 64, 0, 0, 0, 0, 0, 0}},
        {"2 vars deg 128", 2, 2, {128, 128, 0, 0, 0, 0, 0, 0}},
        {"3 vars deg 4", 3, 3, {4, 4, 4, 0, 0, 0, 0, 0}},
        {"3 vars deg 8", 3, 3, {8, 8, 8, 0, 0, 0, 0, 0}},
        {"3 vars deg 12", 3, 3, {12, 12, 12, 0, 0, 0, 0, 0}},
        {"3 vars deg 16", 3, 3, {16, 16, 16, 0, 0, 0, 0, 0}},
        {"4 vars deg 3", 4, 4, {3, 3, 3, 3, 0, 0, 0, 0}},
        {"4 vars deg 4", 4, 4, {4, 4, 4, 4, 0, 0, 0, 0}},
        {"4 vars deg 5", 4, 4, {5, 5, 5, 5, 0, 0, 0, 0}},
        {"4 vars deg 6", 4, 4, {6, 6, 6, 6, 0, 0, 0, 0}},
        {"5 vars deg 2", 5, 5, {2, 2, 2, 2, 2, 0, 0, 0}},
        {"5 vars deg 3", 5, 5, {3, 3, 3, 3, 3, 0, 0, 0}},
        {"5 vars deg 4", 5, 5, {4, 4, 4, 4, 4, 0, 0, 0}},
        {"6 vars deg 2", 6, 6, {2, 2, 2, 2, 2, 2, 0, 0}},
        {"6 vars deg 3", 6, 6, {3, 3, 3, 3, 3, 3, 0, 0}},
        {"7 vars deg 2", 7, 7, {2, 2, 2, 2, 2, 2, 2, 0}},
        {"8 vars deg 2", 8, 8, {2, 2, 2, 2, 2, 2, 2, 2}}
    };
    const solver_field_case_t field = {"65537", 0, 0, 65537, 1, NULL, "t", NULL};
    const slong total_cases = (slong) (sizeof(perf_cases) / sizeof(perf_cases[0]));
    double total_seconds = 0.0;

    printf("=== Resultant Performance Test ===\n\n");
    printf("Field: 65537\n");
    printf("Total cases: %ld\n\n", total_cases);

    for (slong i = 0; i < total_cases; i++) {
        const solver_degree_case_t *test = &perf_cases[i];
        char *polys_str = NULL;
        char *elim_vars = NULL;
        char *all_vars = NULL;
        char *result = NULL;
        fq_nmod_ctx_t ctx;
        fmpz_t p;
        int saved_stdout = -1;
        int saved_stderr = -1;
        struct timespec start;
        struct timespec end;
        double elapsed;
        ulong seed = 2000003UL + (ulong) i;
        const char *status = "failed";

        printf("[%ld/%ld] %s ... ", i + 1, total_cases, test->label);
        fflush(stdout);

        if (!dtest_generate_random_system_string(&field,
                                                test->degrees,
                                                test->npolys,
                                                test->nvars,
                                                1.0,
                                                seed,
                                                &polys_str)) {
            printf("generation failed\n");
            continue;
        }

        if (!dtest_build_random_system_strings(test->nvars, test->npolys - 1,
                                               &elim_vars, &all_vars)) {
            free(polys_str);
            printf("elim-vars failed\n");
            continue;
        }

        fmpz_init_set_ui(p, field.prime_ui);
        fq_nmod_ctx_init(ctx, p, field.power, field.generator);
        fmpz_clear(p);

        dtest_redirect_fd_to_devnull(STDOUT_FILENO, &saved_stdout);
        dtest_redirect_fd_to_devnull(STDERR_FILENO, &saved_stderr);
        clock_gettime(CLOCK_MONOTONIC, &start);
        if (test->npolys == 2) {
            result = dtest_compute_subres_resultant_str(polys_str, elim_vars, ctx);
        } else {
            resultant_method_t saved_method = g_resultant_method;
            if (test->npolys == 3) {
                g_resultant_method = RESULTANT_METHOD_DIXON_RECURSIVE;
            } else {
                g_resultant_method = RESULTANT_METHOD_DIXON;
            }
            result = dixon_str(polys_str, elim_vars, ctx);
            g_resultant_method = saved_method;
        }
        clock_gettime(CLOCK_MONOTONIC, &end);
        dtest_restore_fd(STDOUT_FILENO, saved_stdout);
        dtest_restore_fd(STDERR_FILENO, saved_stderr);

        elapsed = (double) (end.tv_sec - start.tv_sec) +
                  (double) (end.tv_nsec - start.tv_nsec) / 1000000000.0;
        total_seconds += elapsed;

        if (result) {
            status = "success";
        }

        printf("%.3f s [%s]\n", elapsed, status);

        free(result);
        free(elim_vars);
        free(all_vars);
        fq_nmod_ctx_clear(ctx);
        free(polys_str);
    }

    printf("\nTotal time: %.3f s\n", total_seconds);
    return 0;
}

// Main test function with algorithm comparison
int test_bezout_bound() {
    // Initialize random state
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL) + getpid(), time(NULL) * getpid());
    
    printf("=================================================\n");
    printf("Bezout Bound Test Suite - Algorithm Comparison\n");
    printf("=================================================\n");
    
    // Test configurations
    const int num_reps = 200;
    
    // Define test cases
    typedef struct {
        const char *name;
        slong nvars;
        slong npars;
        slong npolys;
        slong degrees[10];
        double density;
    } test_case_t;
    
    test_case_t test_cases[] = {
        {"Mixed degree (3,4,5) system", 2, 1, 3, {3, 4, 5}, 0.9},
        {"Quintic-sextic (5,6,6) system", 2, 1, 3, {5, 6, 6}, 1.0},
        {"Uniform quintic (5,5,5) system", 2, 1, 3, {5, 5, 5}, 0.8},
        {"Six quadratic equations (2,2,2,2,2,2)", 5, 1, 6, {2, 2, 2, 2, 2, 2}, 0.7},
        {"Mixed quadratic-cubic with parameters (2,2,3,3)", 3, 2, 4, {2, 2, 3, 3}, 0.8}
    };
    int num_test_cases = 5;

    {
        int total_runs = 2 * num_reps * num_test_cases;
        int run_index = 0;

        g_bezout_test_quiet_mode = 1;

        for (int algo = 0; algo <= 1; algo++) {
            if (algo == 0) {
                g_matrix_transpose_threshold = 0;
                printf("\n[Algorithm 1/2] Traditional Max-Rank Submatrix Extraction\n");
            } else {
                g_matrix_transpose_threshold = 1000;
                printf("\n[Algorithm 2/2] Degree-Optimal Submatrix Selection\n");
            }

            for (int rep = 0; rep < num_reps; rep++) {
                printf("  Repetition %d/%d\n", rep + 1, num_reps);

                for (int tc = 0; tc < num_test_cases; tc++) {
                    test_case_t *test = &test_cases[tc];
                    run_index++;
                    printf("    Progress %d/%d: %s\n", run_index, total_runs, test->name);
                    test_dixon_system_quiet(test->name, test->nvars, test->npars,
                                           65537, 1, test->degrees, test->npolys,
                                           test->density, state, tc);
                }
            }
        }

        g_bezout_test_quiet_mode = 0;
        printf("\nAll test runs completed. Final summary:\n");
        print_degree_check_summary();
        cleanup_degree_check_results();
        flint_rand_clear(state);
        flint_cleanup();
        return 0;
    }
    
    // Run tests for both algorithms
    for (int algo = 0; algo <= 1; algo++) {
        // Set algorithm threshold
        if (algo == 0) {
            g_matrix_transpose_threshold = 0;  // Traditional
            printf("\n╔════════════════════════════════════════════════════════════════╗\n");
            printf("║  ALGORITHM: Traditional Max-Rank Submatrix Extraction         ║\n");
            printf("╚════════════════════════════════════════════════════════════════╝\n");
        } else {
            g_matrix_transpose_threshold = 1000;  // Degree-Optimal
            printf("\n╔════════════════════════════════════════════════════════════════╗\n");
            printf("║  ALGORITHM: Degree-Optimal Submatrix Selection                 ║\n");
            printf("╚════════════════════════════════════════════════════════════════╝\n");
        }
        
        // Run all test cases with multiple repetitions
        for (int rep = 0; rep < num_reps; rep++) {
            printf("\n--- Repetition %d/%d ---\n", rep + 1, num_reps);
            
            for (int tc = 0; tc < num_test_cases; tc++) {
                test_case_t *test = &test_cases[tc];
                test_dixon_system(test->name, test->nvars, test->npars, 
                                65537, 1, test->degrees, test->npolys, 
                                test->density, state, tc);  // Pass test case index
            }
        }
    }
    
    // Print comprehensive comparison summary
    print_degree_check_summary();

    // Cleanup
    cleanup_degree_check_results();
    flint_rand_clear(state);
    flint_cleanup();
    
    return 0;
}

// ============================================================================
// Dixon Matrix Size Test Function
// ============================================================================

// Test data structure
typedef struct {
    slong n;
    slong d;
    slong fc_bound;
    slong unrestricted_bound;
    slong ks96_bound;
    slong variable_bound;
    slong true_size;
} dixon_size_test_data_t;

// ============================================================================
// Dixon Matrix Size Test Function
// ============================================================================

void test_dixon_matrix_sizes() {
    printf("\n");
    printf("╔════════════════════════════════════════════════════════════════════════════╗\n");
    printf("║                    DIXON MATRIX SIZE VERIFICATION TEST                     ║\n");
    printf("╚════════════════════════════════════════════════════════════════════════════╝\n");
    printf("\n");
    
    // Test parameters
    struct {
        slong n;
        slong d;
        slong expected;
    } test_cases[] = {
        {3, 2, 5},
        {3, 3, 12},
        {3, 4, 22},
        {3, 5, 35},
        {4, 2, 14},
        {4, 3, 55},
        {4, 4, 140},
        {4, 5, 285},
        {5, 2, 42},
        {5, 3, 273},
        {5, 4, 969},
        {5, 5, 2530},
        {6, 2, 132},
        {6, 3, 1428},
        {6, 4, 7084},
        {6, 5, 23751},
        {7, 2, 429},
        {7, 3, 7752},
        //{7, 4, 53820},
        {8, 2, 1430},
        //{8, 3, 43263},
        {9, 2, 4862},
        {10, 2, 16796}
    };
    
    int num_tests = sizeof(test_cases) / sizeof(test_cases[0]);
    
    printf("Testing Dixon matrix sizes...\n");
    printf("Field: GF(65537)\n");
    printf("\n");
    
    // Print table header
    printf("┌────┬────┬──────────┬──────────┬────────┐\n");
    printf("│ n  │ d  │ Expected │ Computed │ Match? │\n");
    printf("├────┼────┼──────────┼──────────┼────────┤\n");
    
    int passed = 0;
    int failed = 0;
    
    clock_t total_start = clock();
    
    for (int i = 0; i < num_tests; i++) {
        slong n = test_cases[i].n;
        slong d = test_cases[i].d;
        slong expected = test_cases[i].expected;
        
        // Compute actual matrix size
        slong computed = dixon_matrix_size(n, d, 65537, 1);
        
        int match = (computed == expected);
        if (match) {
            passed++;
        } else {
            failed++;
        }
        
        printf("│ %2ld │ %2ld │ %8ld │ %8ld │ %6s │\n",
               n, d, expected, computed, match ? "✓" : "✗");
    }
    
    printf("└────┴────┴──────────┴──────────┴────────┘\n");
    
    clock_t total_end = clock();
    double total_time = (double)(total_end - total_start) / CLOCKS_PER_SEC;
    
    printf("\n");
    printf("═════════════════════════════════════════════════════════════════════════════\n");
    printf("Test Summary:\n");
    printf("─────────────────────────────────────────────────────────────────────────────\n");
    printf("  Total tests: %d\n", num_tests);
    printf("  Passed:      %d\n", passed);
    printf("  Failed:      %d\n", failed);
    printf("  Success rate: %.1f%%\n", 
           num_tests > 0 ? 100.0 * passed / num_tests : 0.0);
    printf("  Total time:  %.3f seconds\n", total_time);
    printf("═════════════════════════════════════════════════════════════════════════════\n");
}


// Print test menu
void print_test_menu() {
    printf("\n");
    printf("╔════════════════════════════════════════════════════════════════════════════╗\n");
    printf("║                       DIXON RESULTANT TEST SUITE                           ║\n");
    printf("╚════════════════════════════════════════════════════════════════════════════╝\n");
    printf("\n");
    printf("Available Test Modes:\n");
    printf("─────────────────────────────────────────────────────────────────────────────\n");
    printf("\n");
    printf("  0 - Display this help menu\n");
    printf("      Shows available test options and descriptions\n");
    printf("\n");
    printf("  1 - Test Dixon Matrix Sizes\n");
    printf("\n");
    printf("  2 - Bezout Bound Verification\n");
    printf("\n");
    printf("  3 - Polynomial Solver Correctness\n");
    printf("\n");
    printf("  4 - Resultant Performance\n");
    printf("\n");
    printf("  5 - XHash Example System\n");
    printf("\n");
    printf("─────────────────────────────────────────────────────────────────────────────\n");
}


void test_dixon(int test_mode){
    switch(test_mode) {
        case 0:
            // Display help menu
            print_test_menu();
            break;
            
        case 1:
            printf("\n>>> Test Dixon Matrix Sizes <<<\n");
            test_dixon_matrix_sizes();
            break;            

        case 2:
            printf("\n>>> Running Bezout Bound Test Suite <<<\n");
            test_bezout_bound();
            break;
            
        case 3:
            printf("\n>>> Running Polynomial Solver Correctness Test <<<\n");
            test_polynomial_solver_correctness();
            break;

        case 4:
            printf("\n>>> Running Resultant Performance Test <<<\n");
            test_polynomial_solver_performance();
            break;

        case 5:
            printf("\n>>> Running XHash Example System <<<\n");
            test_xhash();
            break;
            
        default:
            printf("\n>>> Error: Invalid test mode %d <<<\n", test_mode);
            printf("Valid modes: 0 (help), 1 (Test Dixon Matrix Sizes), 2 (Bezout bound), 3 (Polynomial solver correctness), 4 (Resultant performance), 5 (XHash)\n");
            printf("Run with mode 0 to see detailed descriptions.\n");
            break;
    }
    
    if (test_mode != 0) {
        printf("\n=== Test Complete ===\n");
    }
}
