// dixon_complexity.c - Modified to use Hessenberg method

#include "dixon_complexity.h"

// Comparison function for descending order sorting
int compare_desc(const void *a, const void *b) {
    long long_a = *(const long*)a;
    long long_b = *(const long*)b;
    return (long_b > long_a) - (long_b < long_a);
}

/**
 * Unified boundary height calculation function
 * 
 * @param d_list: degree list
 * @param n: list length
 * @param shift: slope correction (Dixon typically uses -1)
 * @param include_prefix_zero: whether to add a[0] = 0 at the beginning
 *        - 1:  a = [0, s[0], s[0]+s[1], ...] (for determinant/Hessenberg)
 *        - 0:  a = [s[0], s[0]+s[1], ...]     (for DP)
 * @param a_out: output boundary array
 * @param slopes_out: output slopes array (can be NULL)
 */
void get_boundary_heights(long *d_list, long n, long shift, 
                          int include_prefix_zero,
                          long **a_out, long **slopes_out) {
    // 1. Calculate slopes and sort in descending order
    long *slopes = (long *)malloc(n * sizeof(long));
    for (long i = 0; i < n; i++) {
        slopes[i] = d_list[i] + shift;
    }
    qsort(slopes, n, sizeof(long), compare_desc);
    
    // 2. Calculate upper bound a
    long *a;
    if (include_prefix_zero) {
        // Determinant/Hessenberg mode: a = [0, s[0], s[0]+s[1], ...]
        a = (long *)calloc(n, sizeof(long));
        a[0] = 0;
        long current_h = 0;
        for (long i = 0; i < n - 1; i++) {
            current_h += slopes[i];
            a[i + 1] = current_h;
        }
    } else {
        // DP mode: a = [s[0], s[0]+s[1], ...]
        a = (long *)malloc(n * sizeof(long));
        long current_h = 0;
        for (long i = 0; i < n; i++) {
            current_h += slopes[i];
            a[i] = current_h;
        }
    }
    
    *a_out = a;
    if (slopes_out != NULL) {
        *slopes_out = slopes;
    } else {
        free(slopes);
    }
}

/**
 * Calculate Dixon matrix size using Hessenberg recurrence - O(n^2)
 * This is the main method replacing the old inequality system approach
 */
void dixon_size(fmpz_t result, const long *d_list, int len, int show_details) {
    fmpz_set_ui(result, 0);
    
    if (len == 0) {
        return;
    }
    
    // For Dixon resultant, we use shift = -1
    long shift = -1;
    long n = len;
    
    // Get boundary heights
    long *a, *slopes;
    get_boundary_heights((long*)d_list, n, shift, 1, &a, &slopes);
    
    if (show_details) {
        printf("Input degrees: [");
        for (int i = 0; i < len; i++) {
            printf("%ld%s", d_list[i], i < len - 1 ? ", " : "");
        }
        printf("]\n");
        
        printf("Slopes (sorted): [");
        for (long i = 0; i < n; i++) {
            printf("%ld%s", slopes[i], i < n - 1 ? ", " : "");
        }
        printf("]\n");
        
        printf("Boundary a_i: [");
        for (long i = 0; i < n; i++) {
            printf("%ld%s", a[i], i < n - 1 ? ", " : "");
        }
        printf("]\n");
    }
    
    // D[k] stores the determinant of the first k x k submatrix
    fmpz *D = (fmpz *)malloc((n + 1) * sizeof(fmpz));
    for (long i = 0; i <= n; i++) {
        fmpz_init(D + i);
    }
    fmpz_one(D + 0);  // D[0] = 1
    
    fmpz_t m_val, term, sum_val;
    fmpz_init(m_val);
    fmpz_init(term);
    fmpz_init(sum_val);
    
    // Hessenberg recurrence
    for (long k = 1; k <= n; k++) {
        fmpz_zero(sum_val);
        
        for (long i = 1; i <= k; i++) {
            long row_idx = i - 1;
            long col_idx = k - 1;
            
            // Binomial coefficient: binom(a[row_idx] + 1, col_idx - row_idx + 1)
            long upper = a[row_idx] + 1;
            long lower = col_idx - row_idx + 1;
            
            if (lower >= 0 && lower <= upper) {
                fmpz_bin_uiui(m_val, (unsigned long)upper, (unsigned long)lower);
            } else {
                fmpz_zero(m_val);
            }
            
            // Sign: (-1)^(k-i)
            fmpz_mul(term, m_val, D + (i - 1));
            if ((k - i) % 2 == 1) {
                fmpz_neg(term, term);
            }
            
            fmpz_add(sum_val, sum_val, term);
        }
        
        fmpz_set(D + k, sum_val);
    }
    
    fmpz_set(result, D + n);
    
    if (show_details) {
        printf("Dixon matrix size (Hessenberg): ");
        fmpz_print(result);
        printf("\n");
    }
    
    // Cleanup
    for (long i = 0; i <= n; i++) {
        fmpz_clear(D + i);
    }
    free(D);
    fmpz_clear(m_val);
    fmpz_clear(term);
    fmpz_clear(sum_val);
    free(a);
    free(slopes);
}

// Calculate Dixon complexity
double dixon_complexity(const long *a_values, int len, int n, double omega) {
    fmpz_t size;
    fmpz_init(size);
    
    dixon_size(size, a_values, len, 0);
    
    // Check if size is zero to avoid computation issues
    if (fmpz_is_zero(size)) {
        fmpz_clear(size);
        return 0.0;
    }
    
    long d = 0;
    int m = len;
    
    if (m == n + 1) {
        d = 1;
    } else if (m == n) {
        for (int i = 0; i < len; i++) {
            d += a_values[i];
        }
    } else if (m < n) {
        for (int i = 0; i < len; i++) {
            d += a_values[i];
        }
        // Calculate (d+1)^(n-m+1)
        long exponent = n - m + 1;
        long base = d + 1;
        
        // Use double to avoid integer overflow
        double d_double = pow((double)base, (double)exponent);
        
        if (!isfinite(d_double) || d_double <= 0) {
            fmpz_clear(size);
            return INFINITY;
        }
        
        // Calculate log₂(d * size^omega)
        double size_double = fmpz_get_d(size);
        
        if (size_double <= 0 || !isfinite(size_double)) {
            fmpz_clear(size);
            return 0.0;
        }
        
        double powered_size = pow(size_double, omega);
        
        if (!isfinite(powered_size) || powered_size <= 0) {
            fmpz_clear(size);
            return INFINITY;
        }
        
        double product = d_double * powered_size;
        
        if (!isfinite(product) || product <= 0) {
            fmpz_clear(size);
            return INFINITY;
        }
        
        double result = log2(product);
        fmpz_clear(size);
        return result;
        
    } else {
        fmpz_clear(size);
        return 0.0;
    }
    
    // Handle m >= n cases
    if (d <= 0) {
        fmpz_clear(size);
        return 0.0;
    }
    
    double size_double = fmpz_get_d(size);
    
    if (size_double <= 0 || !isfinite(size_double)) {
        fmpz_clear(size);
        return 0.0;
    }
    
    double powered_size = pow(size_double, omega);
    
    if (!isfinite(powered_size) || powered_size <= 0) {
        fmpz_clear(size);
        return INFINITY;
    }
    
    double product = (double)d * powered_size;
    
    if (!isfinite(product) || product <= 0) {
        fmpz_clear(size);
        return INFINITY;
    }
    
    double result = log2(product);
    
    fmpz_clear(size);
    return result;
}

static double log2_add_exp(double a, double b) {
    if (!isfinite(a)) return b;
    if (!isfinite(b)) return a;
    if (a < b) {
        double tmp = a;
        a = b;
        b = tmp;
    }
    return a + log2(1.0 + pow(2.0, b - a));
}

static double log2_factorial_slong(slong n) {
    double acc = 0.0;
    for (slong i = 2; i <= n; i++) {
        acc += log2((double) i);
    }
    return acc;
}

static double log2_binomial_upper(slong n, slong k) {
    if (k < 0 || k > n) return -INFINITY;
    long double ln_binom =
        lgammal((long double) n + 1.0L) -
        lgammal((long double) k + 1.0L) -
        lgammal((long double) (n - k) + 1.0L);
    return (double) (ln_binom / logl(2.0L));
}

static slong saturated_binomial_slong(slong n, slong k) {
    slong result = 1;

    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;

    for (slong i = 1; i <= k; i++) {
        if (result > WORD_MAX / (n - k + i)) {
            return WORD_MAX;
        }
        result = result * (n - k + i) / i;
    }

    return result;
}

static double log2_soft_fft_multiply_from_degree_log2(double log2_degree_bound) {
    if (!(log2_degree_bound > 0.0) || !isfinite(log2_degree_bound)) {
        return 0.0;
    }

    double result = log2_degree_bound;
    if (log2_degree_bound > 1.0) {
        result += log2(log2_degree_bound);
    }
    if (log2_degree_bound > 4.0) {
        result += log2(log2(log2_degree_bound));
    }
    return result;
}

static int compare_long_desc_qsort(const void *a, const void *b) {
    long av = *(const long *) a;
    long bv = *(const long *) b;
    if (av < bv) return 1;
    if (av > bv) return -1;
    return 0;
}

static void build_sorted_degree_prefix(long **sorted_out,
                                       long **prefix_out,
                                       const long *degrees,
                                       slong num_polys) {
    long *sorted = NULL;
    long *prefix = NULL;

    if (num_polys > 0) {
        sorted = (long *) malloc((size_t) num_polys * sizeof(long));
        prefix = (long *) calloc((size_t) (num_polys + 1), sizeof(long));
    }

    for (slong i = 0; i < num_polys; i++) {
        sorted[i] = degrees[i] > 0 ? degrees[i] : 0;
    }
    qsort(sorted, (size_t) num_polys, sizeof(long), compare_long_desc_qsort);

    for (slong i = 0; i < num_polys; i++) {
        prefix[i + 1] = prefix[i] + sorted[i];
    }

    *sorted_out = sorted;
    *prefix_out = prefix;
}

static long prefix_degree_bound(const long *prefix, slong count) {
    long bound;

    if (count <= 0) {
        return 0;
    }

    bound = prefix[count] - (long) count;
    return bound > 0 ? bound : 0;
}

static double log2_kronecker_degree_from_bounds(const long *var_degree_bounds,
                                                slong num_vars) {
    double log2_weight = 0.0;
    double log2_degree = -INFINITY;

    for (slong i = 0; i < num_vars; i++) {
        long bound = var_degree_bounds[i];
        if (bound > 0) {
            log2_degree = log2_add_exp(log2_degree,
                                       log2((double) bound) + log2_weight);
        }
        log2_weight += log2((double) bound + 1.0);
    }

    if (!isfinite(log2_degree) || log2_degree < 0.0) {
        return 0.0;
    }

    return log2_degree;
}

static int step1_row_highest_active_var_index(slong row_idx,
                                              slong num_elim_vars,
                                              slong num_parameter_vars) {
    slong total_vars = 2 * num_elim_vars + num_parameter_vars;

    if (total_vars <= 0) {
        return -1;
    }

    if (num_parameter_vars > 0) {
        return (int) (total_vars - 1);
    }

    if (row_idx <= 0) {
        return num_elim_vars > 0 ? (int) (num_elim_vars - 1) : -1;
    }

    if (row_idx <= num_elim_vars) {
        return (int) (num_elim_vars + row_idx - 1);
    }

    return -1;
}

static double step1_entry_degree_log2(slong row_idx,
                                      long col_degree,
                                      const double *var_weight_log2,
                                      slong num_elim_vars,
                                      slong num_parameter_vars) {
    int highest_idx = step1_row_highest_active_var_index(row_idx,
                                                         num_elim_vars,
                                                         num_parameter_vars);
    long effective_degree = col_degree;

    if (row_idx > 0) {
        effective_degree -= 1;
    }

    if (highest_idx < 0 || effective_degree <= 0) {
        return -INFINITY;
    }

    return log2((double) effective_degree) + var_weight_log2[highest_idx];
}

void dixon_complexity_report_from_degrees(dixon_complexity_report_t *report,
                                          const long *degrees,
                                          slong num_polys,
                                          slong num_all_vars,
                                          slong num_elim_vars,
                                          slong num_parameter_vars,
                                          double omega) {
    memset(report, 0, sizeof(*report));

    report->det_size = num_polys;
    report->num_all_vars = num_all_vars;
    report->num_elim_vars = num_elim_vars;
    report->num_parameter_vars = num_parameter_vars;
    report->step1_var_count = 2 * num_elim_vars + num_parameter_vars;
    if (report->step1_var_count < 0) {
        report->step1_var_count = 0;
    }

    {
        int standard_step1_shape =
            (num_polys > 0 &&
             num_elim_vars >= 0 &&
             num_polys == num_elim_vars + 1 &&
             report->step1_var_count >= 0);
        long *sorted = NULL;
        long *prefix = NULL;
        long *var_degree_bounds = NULL;
        double *var_weight_log2 = NULL;
        double log2_row_sum = -INFINITY;
        double log2_col_sum = -INFINITY;

        for (slong i = 0; i < num_polys; i++) {
            long di = degrees[i] > 0 ? degrees[i] : 0;
            report->step1_det_total_degree += di;
        }

        report->step1_det_total_degree -= (num_elim_vars > 0 ? num_elim_vars : 0);
        if (report->step1_det_total_degree < 0) {
            report->step1_det_total_degree = 0;
        }

        if (!standard_step1_shape || report->step1_var_count <= 0) {
            double log2_det_uni_degree = 0.0;

            if (report->step1_det_total_degree > 0 && report->step1_var_count > 0) {
                log2_det_uni_degree =
                    ((double) report->step1_var_count) *
                    log2((double) report->step1_det_total_degree + 1.0);
            }
            report->step1_kronecker_degree_log2 = log2_det_uni_degree;
            report->step1_direct_log2 =
                log2_factorial_slong(report->det_size) +
                log2_soft_fft_multiply_from_degree_log2(log2_det_uni_degree);

            {
                double log2_col_avg = log2_det_uni_degree;
                double log2_row_avg = log2_det_uni_degree;
                double log2_s = log2_col_avg;
                if (isfinite(log2_row_avg) && (!isfinite(log2_s) || log2_row_avg < log2_s)) {
                    log2_s = log2_row_avg;
                }
                if (!isfinite(log2_s) || log2_s < 0.0) {
                    log2_s = 0.0;
                }
                report->step1_hnf_log2 =
                    (report->det_size > 1 ? omega * log2((double) report->det_size) : 0.0) +
                    log2_s;
            }
        } else {
            build_sorted_degree_prefix(&sorted, &prefix, degrees, num_polys);

            if (report->step1_var_count > 0) {
                var_degree_bounds = (long *) calloc((size_t) report->step1_var_count,
                                                    sizeof(long));
                var_weight_log2 = (double *) calloc((size_t) report->step1_var_count,
                                                    sizeof(double));
            }

            for (slong i = 0; i < num_elim_vars; i++) {
                slong count_orig = i + 1;
                slong count_dual = num_elim_vars - i;
                long orig_bound = prefix_degree_bound(prefix, count_orig);
                long dual_bound = prefix_degree_bound(prefix, count_dual);

                if (orig_bound < 0) orig_bound = 0;
                if (dual_bound < 0) dual_bound = 0;

                var_degree_bounds[i] = orig_bound;
                var_degree_bounds[num_elim_vars + i] = dual_bound;
            }

            for (slong i = 0; i < num_parameter_vars; i++) {
                long param_bound = report->step1_det_total_degree;
                if (param_bound < 0) param_bound = 0;
                var_degree_bounds[2 * num_elim_vars + i] = param_bound;
            }

            if (report->step1_var_count > 0) {
                double log2_weight = 0.0;
                for (slong i = 0; i < report->step1_var_count; i++) {
                    var_weight_log2[i] = log2_weight;
                    log2_weight += log2((double) var_degree_bounds[i] + 1.0);
                }
            }

            report->step1_kronecker_degree_log2 =
                log2_kronecker_degree_from_bounds(var_degree_bounds,
                                                  report->step1_var_count);
            report->step1_direct_log2 =
                log2_factorial_slong(report->det_size) +
                log2_soft_fft_multiply_from_degree_log2(report->step1_kronecker_degree_log2);

            for (slong row = 0; row < num_polys; row++) {
                long row_max_degree = 0;
                double row_log2_degree = -INFINITY;
                if (row == 0) {
                    for (slong j = 0; j < num_polys; j++) {
                        long dj = degrees[j] > 0 ? degrees[j] : 0;
                        if (dj > row_max_degree) {
                            row_max_degree = dj;
                        }
                    }
                } else {
                    for (slong j = 0; j < num_polys; j++) {
                        long dj = degrees[j] > 0 ? degrees[j] : 0;
                        if (dj > 0 && dj - 1 > row_max_degree) {
                            row_max_degree = dj - 1;
                        }
                    }
                }

                if (row_max_degree > 0 && var_weight_log2 != NULL) {
                    long effective_col_degree = (row == 0) ? row_max_degree : (row_max_degree + 1);
                    row_log2_degree = step1_entry_degree_log2(row,
                                                              effective_col_degree,
                                                              var_weight_log2,
                                                              num_elim_vars,
                                                              num_parameter_vars);
                }
                log2_row_sum = log2_add_exp(log2_row_sum, row_log2_degree);
            }

            for (slong col = 0; col < num_polys; col++) {
                long dj = degrees[col] > 0 ? degrees[col] : 0;
                double col_log2_degree = -INFINITY;
                for (slong row = 0; row < num_polys; row++) {
                    double entry_log2 = step1_entry_degree_log2(row,
                                                                dj,
                                                                var_weight_log2,
                                                                num_elim_vars,
                                                                num_parameter_vars);
                    if (!isfinite(col_log2_degree) || entry_log2 > col_log2_degree) {
                        col_log2_degree = entry_log2;
                    }
                }
                log2_col_sum = log2_add_exp(log2_col_sum, col_log2_degree);
            }

            {
                double log2_row_avg = -INFINITY;
                double log2_col_avg = -INFINITY;
                double log2_s = -INFINITY;

                if (num_polys > 0 && isfinite(log2_row_sum)) {
                    log2_row_avg = log2_row_sum - log2((double) num_polys);
                }
                if (num_polys > 0 && isfinite(log2_col_sum)) {
                    log2_col_avg = log2_col_sum - log2((double) num_polys);
                }

                log2_s = log2_col_avg;
                if (isfinite(log2_row_avg) && (!isfinite(log2_s) || log2_row_avg < log2_s)) {
                    log2_s = log2_row_avg;
                }
                if (!isfinite(log2_s) || log2_s < 0.0) {
                    log2_s = 0.0;
                }

                report->step1_hnf_log2 =
                    (report->det_size > 1 ? omega * log2((double) report->det_size) : 0.0) +
                    log2_s;
            }

            free(sorted);
            free(prefix);
            free(var_degree_bounds);
            free(var_weight_log2);
        }
    }

    {
        double total_degree_term_log2 =
            log2_binomial_upper(report->step1_det_total_degree + report->step1_var_count,
                                report->step1_var_count);
        double box_term_log2 = 0.0;
        double log2_probe_unit = 0.0;
        double log2_det_inv = (report->det_size > 1)
            ? omega * log2((double) report->det_size)
            : 0.0;
        double log2_trace = (report->det_size > 0 && report->step1_var_count > 0)
            ? log2((double) report->step1_var_count) +
              2.0 * log2((double) report->det_size)
            : 0.0;

        if (report->step1_var_count > 0) {
            int standard_step1_shape =
                (num_polys > 0 &&
                 num_elim_vars >= 0 &&
                 num_polys == num_elim_vars + 1);
            if (standard_step1_shape) {
                long *sorted = NULL;
                long *prefix = NULL;
                build_sorted_degree_prefix(&sorted, &prefix, degrees, num_polys);
                for (slong i = 0; i < num_elim_vars; i++) {
                    long orig_bound = prefix_degree_bound(prefix, i + 1);
                    long dual_bound = prefix_degree_bound(prefix, num_elim_vars - i);
                    box_term_log2 += log2((double) orig_bound + 1.0);
                    box_term_log2 += log2((double) dual_bound + 1.0);
                }
                for (slong i = 0; i < num_parameter_vars; i++) {
                    box_term_log2 += log2((double) report->step1_det_total_degree + 1.0);
                }
                free(sorted);
                free(prefix);
            } else {
                box_term_log2 =
                    ((double) report->step1_var_count) *
                    log2((double) report->step1_det_total_degree + 1.0);
            }
        }

        report->step1_sparse_term_bound_log2 = total_degree_term_log2;
        if (isfinite(box_term_log2) &&
            (!isfinite(report->step1_sparse_term_bound_log2) ||
             box_term_log2 < report->step1_sparse_term_bound_log2)) {
            report->step1_sparse_term_bound_log2 = box_term_log2;
        }
        if (!isfinite(report->step1_sparse_term_bound_log2) ||
            report->step1_sparse_term_bound_log2 < 0.0) {
            report->step1_sparse_term_bound_log2 = 0.0;
        }

        log2_probe_unit = log2_add_exp(log2_det_inv, log2_trace);
        report->step1_sparse_log2 =
            report->step1_sparse_term_bound_log2 +
            log2_probe_unit;
    }

    {
        long max_degree = 0;
        double entry_param_log2 = 0.0;
        for (slong i = 0; i < num_polys; i++) {
            long di = degrees[i] > 0 ? degrees[i] : 0;
            if (di > max_degree) {
                max_degree = di;
            }
            report->macaulay_degree += di;
        }
        report->macaulay_degree -= num_elim_vars;
        if (report->macaulay_degree < 0) {
            report->macaulay_degree = 0;
        }

        if (num_elim_vars > 0) {
            report->macaulay_cols =
                saturated_binomial_slong(num_elim_vars + report->macaulay_degree,
                                         num_elim_vars);
            report->macaulay_rows = 0;
            for (slong i = 0; i < num_polys; i++) {
                slong multiplier_degree =
                    report->macaulay_degree - (degrees[i] > 0 ? degrees[i] : 0);
                if (multiplier_degree < 0) multiplier_degree = 0;
                report->macaulay_rows +=
                    saturated_binomial_slong(num_elim_vars + multiplier_degree,
                                             num_elim_vars);
            }
            report->macaulay_square_size = FLINT_MIN(report->macaulay_rows,
                                                     report->macaulay_cols);
        }

        if (num_parameter_vars > 0 && max_degree > 0) {
            entry_param_log2 =
                ((double) num_parameter_vars) * log2((double) max_degree + 1.0);
        }
        report->macaulay_log2 =
            (report->macaulay_square_size > 1
                ? omega * log2((double) report->macaulay_square_size)
                : 0.0) +
            entry_param_log2;

        {
            slong gb_n = num_all_vars > 0 ? num_all_vars : num_elim_vars;
            if (gb_n < 1) gb_n = 1;
            report->grobner_dreg = gb_n * (max_degree > 0 ? (max_degree - 1) : 0) + 1;
            report->grobner_log2 =
                omega * log2_binomial_upper(gb_n + report->grobner_dreg, gb_n);
            if (!isfinite(report->grobner_log2) || report->grobner_log2 < 0.0) {
                report->grobner_log2 = 0.0;
            }
        }
    }

    report->step4_log2 = dixon_complexity(degrees, (int) num_polys, (int) num_all_vars, omega);
    report->total_direct_log2 = log2_add_exp(report->step1_direct_log2, report->step4_log2);
    report->total_hnf_log2 = log2_add_exp(report->step1_hnf_log2, report->step4_log2);
    report->total_sparse_log2 = log2_add_exp(report->step1_sparse_log2, report->step4_log2);
}

static int complexity_is_elimination_var(const char *var_name,
                                         char *const *elim_vars,
                                         slong num_elim_vars) {
    for (slong i = 0; i < num_elim_vars; i++) {
        if (strcmp(var_name, elim_vars[i]) == 0) {
            return 1;
        }
    }
    return 0;
}

static void dixon_complexity_print_field_label(FILE *out,
                                               const fmpz_t prime,
                                               ulong power) {
    if (fmpz_is_zero(prime)) {
        fprintf(out, "Q");
        return;
    }

    fprintf(out, "F_");
    fmpz_fprint(out, prime);
    if (power > 1) {
        fprintf(out, "^%lu", power);
        if (fmpz_abs_fits_ui(prime)) {
            mp_limb_t prime_ui = fmpz_get_ui(prime);
            mp_limb_t field_size = 1;
            int overflow = 0;

            for (ulong i = 0; i < power; i++) {
                if (prime_ui != 0 && field_size > UWORD_MAX / prime_ui) {
                    overflow = 1;
                    break;
                }
                field_size *= prime_ui;
            }

            if (!overflow) {
                fprintf(out, " (size %lu)", field_size);
            }
        }
    }
}

static void dixon_complexity_fprint_name_list(FILE *fp,
                                              char *const *names,
                                              slong count) {
    if (count <= 0) {
        fprintf(fp, "(none)");
        return;
    }

    for (slong i = 0; i < count; i++) {
        if (i > 0) {
            fprintf(fp, ", ");
        }
        fprintf(fp, "%s", names[i]);
    }
}

static void dixon_complexity_fprint_parameter_vars(FILE *fp,
                                                   char *const *all_vars,
                                                   slong num_all_vars,
                                                   char *const *elim_vars,
                                                   slong num_elim_vars) {
    int first = 1;

    for (slong i = 0; i < num_all_vars; i++) {
        if (!complexity_is_elimination_var(all_vars[i], elim_vars, num_elim_vars)) {
            if (!first) {
                fprintf(fp, ", ");
            }
            fprintf(fp, "%s", all_vars[i]);
            first = 0;
        }
    }

    if (first) {
        fprintf(fp, "(none)");
    }
}

static double dixon_complexity_best_total_log2(const dixon_complexity_report_t *report) {
    double best = report->total_direct_log2;
    if (report->total_hnf_log2 < best) {
        best = report->total_hnf_log2;
    }
    if (report->total_sparse_log2 < best) {
        best = report->total_sparse_log2;
    }
    return best;
}

static void dixon_complexity_write_report_body(
        FILE *fp,
        slong num_polys,
        slong num_all_vars,
        char **all_vars,
        char **elim_var_list,
        slong num_elim,
        slong num_parameter_vars,
        const long *degrees,
        const fmpz_t matrix_size,
        long bezout_bound,
        const dixon_complexity_report_t *report,
        double omega) {
    fprintf(fp, "--- Raw parameters ---\n");
    fprintf(fp, "Equations: %ld\n", num_polys);
    fprintf(fp, "All vars (%ld): ", num_all_vars);
    dixon_complexity_fprint_name_list(fp, all_vars, num_all_vars);
    fprintf(fp, "\n");
    fprintf(fp, "Elimination vars (%ld): ", num_elim);
    dixon_complexity_fprint_name_list(fp, elim_var_list, num_elim);
    fprintf(fp, "\n");
    fprintf(fp, "Parameter vars (%ld): ", num_parameter_vars);
    dixon_complexity_fprint_parameter_vars(fp, all_vars, num_all_vars,
                                           elim_var_list, num_elim);
    fprintf(fp, "\n");
    fprintf(fp, "Degree sequence: [");
    for (slong i = 0; i < num_polys; i++) {
        if (i > 0) {
            fprintf(fp, ", ");
        }
        fprintf(fp, "%ld", degrees[i]);
    }
    fprintf(fp, "]\n");
    fprintf(fp, "Bezout bound (degree product): %ld\n", bezout_bound);

    fprintf(fp, "\n--- Step 1 ---\n");
    fprintf(fp, "Cancellation matrix size: %ld x %ld\n",
            report->det_size, report->det_size);
    fprintf(fp, "Step 1 indeterminates (2*elim + params): %ld = 2*%ld + %ld\n",
            report->step1_var_count, report->num_elim_vars, report->num_parameter_vars);
    fprintf(fp, "Step 1 determinant total degree upper bound: %ld\n",
            report->step1_det_total_degree);
    fprintf(fp, "Step 1 Kronecker univariate degree upper bound (log2): %.6f\n",
            report->step1_kronecker_degree_log2);
    fprintf(fp, "Step 1 sparse term upper bound (log2 T): %.6f\n",
            report->step1_sparse_term_bound_log2);
    fprintf(fp, "Step 1 direct Leibniz + Kronecker/FFT (log2): %.6f\n",
            report->step1_direct_log2);
    fprintf(fp, "Step 1 Kronecker + HNF (log2): %.6f\n",
            report->step1_hnf_log2);
    fprintf(fp, "Step 1 derivative sparse interpolation (log2): %.6f\n",
            report->step1_sparse_log2);
    if (num_polys != num_elim + 1) {
        fprintf(fp,
                "Step 1 note: standard Dixon resultant shape expects #polys = #elim + 1; current input is %ld vs %ld + 1.\n",
                num_polys, num_elim);
    }

    fprintf(fp, "\n--- Step 4 ---\n");
    fprintf(fp, "Step 4 ambient variable count n: %ld\n", report->num_all_vars);
    fprintf(fp, "Dixon matrix size: ");
    fmpz_fprint(fp, matrix_size);
    fprintf(fp, "\n");
    fprintf(fp, "Step 4 resultant degree estimate (Bezout): %ld\n", bezout_bound);
    fprintf(fp, "Step 4 Dixon complexity (log2, omega=%.4g): %.6f\n",
            omega, report->step4_log2);

    fprintf(fp, "\n--- Overall ---\n");
    fprintf(fp, "Total with direct step 1 (log2): %.6f\n",
            report->total_direct_log2);
    fprintf(fp, "Total with HNF step 1 (log2): %.6f\n",
            report->total_hnf_log2);
    fprintf(fp, "Total with sparse step 1 (log2): %.6f\n",
            report->total_sparse_log2);
    fprintf(fp, "Best overall estimate used for encoding (log2): %.6f\n",
            dixon_complexity_best_total_log2(report));

    fprintf(fp, "\n--- Comparison: Macaulay / Groebner ---\n");
    fprintf(fp, "Macaulay degree: %ld\n", report->macaulay_degree);
    fprintf(fp, "Macaulay matrix size upper bound: %ld x %ld (square <= %ld)\n",
            report->macaulay_rows, report->macaulay_cols, report->macaulay_square_size);
    fprintf(fp, "Macaulay resultant + HNF estimate (log2): %.6f\n",
            report->macaulay_log2);
    fprintf(fp, "Groebner degree of regularity estimate: %ld\n",
            report->grobner_dreg);
    fprintf(fp, "Groebner basis estimate (log2): %.6f\n",
            report->grobner_log2);
}

static void save_comp_result_to_file(
        const char   *filename,
        const char   *polys_str,
        const char   *vars_str,
        const fmpz_t  prime,
        ulong         power,
        slong         num_polys,
        slong         num_all_vars,
        char        **all_vars,
        char        **elim_var_list,
        slong         num_elim,
        slong         num_parameter_vars,
        const long   *degrees,
        const fmpz_t  matrix_size,
        long          bezout_bound,
        const dixon_complexity_report_t *report,
        double        omega,
        double        comp_time) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Warning: Cannot create output file '%s'\n", filename);
        return;
    }

    fprintf(fp, "Dixon Complexity Analysis\n");
    fprintf(fp, "=========================\n");
    fprintf(fp, "Field: ");
    dixon_complexity_print_field_label(fp, prime, power);
    fprintf(fp, "\n");
    fprintf(fp, "Polynomials: %s\n", polys_str);
    fprintf(fp, "Eliminate:   %s\n", vars_str);
    fprintf(fp, "Computation time: %.3f seconds\n\n", comp_time);

    dixon_complexity_write_report_body(fp,
                                       num_polys,
                                       num_all_vars,
                                       all_vars,
                                       elim_var_list,
                                       num_elim,
                                       num_parameter_vars,
                                       degrees,
                                       matrix_size,
                                       bezout_bound,
                                       report,
                                       omega);
    fclose(fp);
}

void run_complexity_analysis(
        const char      *polys_str,
        const char      *vars_str,
        const fmpz_t     prime,
        ulong            power,
        const fq_nmod_ctx_t ctx,
        const char      *output_filename,
        int              silent_mode,
        double           comp_time,
        double           omega) {
    slong num_polys = 0;
    slong num_elim = 0;
    char **poly_arr = split_string(polys_str, &num_polys);
    char **elim_arr = split_string(vars_str, &num_elim);
    char *gen_name = (ctx == NULL) ? NULL : get_generator_name(ctx);
    char **all_vars = NULL;
    slong num_all_vars = 0;
    slong num_parameter_vars = 0;
    long *degrees = NULL;
    long bezout = 1;
    fmpz_t matrix_size;
    dixon_complexity_report_t report;

    collect_variables((const char **) poly_arr, num_polys,
                      gen_name, &all_vars, &num_all_vars);

    for (slong i = 0; i < num_all_vars; i++) {
        if (!complexity_is_elimination_var(all_vars[i], elim_arr, num_elim)) {
            num_parameter_vars++;
        }
    }

    if (num_polys <= 0) {
        if (!silent_mode) {
            fprintf(stderr, "Error: no polynomials to analyze\n");
        }
        free_split_strings(poly_arr, num_polys);
        free_split_strings(elim_arr, num_elim);
        if (gen_name) {
            free(gen_name);
        }
        for (slong i = 0; i < num_all_vars; i++) {
            free(all_vars[i]);
        }
        free(all_vars);
        return;
    }

    degrees = (long *) calloc((size_t) num_polys, sizeof(long));
    for (slong i = 0; i < num_polys; i++) {
        degrees[i] = get_poly_total_degree(poly_arr[i], gen_name);
        bezout *= degrees[i];
    }

    fmpz_init(matrix_size);
    dixon_size(matrix_size, degrees, (int) num_polys, 0);

    dixon_complexity_report_from_degrees(&report,
                                         degrees,
                                         num_polys,
                                         num_all_vars,
                                         num_elim,
                                         num_parameter_vars,
                                         omega);

    if (!silent_mode) {
        printf("\n=== Complexity Analysis ===\n");
        dixon_complexity_write_report_body(stdout,
                                           num_polys,
                                           num_all_vars,
                                           all_vars,
                                           elim_arr,
                                           num_elim,
                                           num_parameter_vars,
                                           degrees,
                                           matrix_size,
                                           bezout,
                                           &report,
                                           omega);
        if (output_filename) {
            printf("\nReport saved to: %s\n", output_filename);
        }
        printf("===========================\n");
    }

    if (output_filename) {
        save_comp_result_to_file(output_filename, polys_str, vars_str,
                                 prime, power,
                                 num_polys, num_all_vars, all_vars,
                                 elim_arr, num_elim, num_parameter_vars,
                                 degrees, matrix_size, bezout,
                                 &report, omega, comp_time);
    }

    fmpz_clear(matrix_size);
    free(degrees);
    for (slong i = 0; i < num_all_vars; i++) {
        free(all_vars[i]);
    }
    free(all_vars);
    if (gen_name) {
        free(gen_name);
    }
    free_split_strings(poly_arr, num_polys);
    free_split_strings(elim_arr, num_elim);
}

// Initialize polynomial analysis structure
static void poly_analysis_init(poly_analysis_t *analysis, slong num_polys, const fq_nmod_ctx_t ctx) {
    analysis->max_vars = 16;
    analysis->all_vars = (char**) malloc(analysis->max_vars * sizeof(char*));
    analysis->num_all_vars = 0;
    analysis->degrees = (long*) calloc(num_polys, sizeof(long));
    analysis->num_polys = num_polys;
    analysis->ctx = ctx;
}

// Clear polynomial analysis structure
static void poly_analysis_clear(poly_analysis_t *analysis) {
    for (slong i = 0; i < analysis->num_all_vars; i++) {
        free(analysis->all_vars[i]);
    }
    free(analysis->all_vars);
    free(analysis->degrees);
}

// Simple hash function for variable names
static slong hash_string(const char *str, slong bucket_count) {
    slong hash = 5381;
    while (*str) {
        hash = ((hash << 5) + hash) + *str++;
    }
    return hash < 0 ? (-hash) % bucket_count : hash % bucket_count;
}

// Initialize variable hash table
static void var_hash_init(var_hash_table_t *table, slong initial_buckets) {
    table->bucket_count = initial_buckets;
    table->buckets = (var_entry_t**) calloc(initial_buckets, sizeof(var_entry_t*));
    table->count = 0;
}

// Clear variable hash table
static void var_hash_clear(var_hash_table_t *table) {
    for (slong i = 0; i < table->bucket_count; i++) {
        var_entry_t *entry = table->buckets[i];
        while (entry) {
            var_entry_t *next = entry->next;
            free(entry->name);
            free(entry);
            entry = next;
        }
    }
    free(table->buckets);
    table->buckets = NULL;
    table->bucket_count = 0;
    table->count = 0;
}

// Find variable in hash table (returns index or -1 if not found)
static slong var_hash_find(var_hash_table_t *table, const char *name) {
    if (!table->buckets) return -1;
    
    slong bucket = hash_string(name, table->bucket_count);
    var_entry_t *entry = table->buckets[bucket];
    
    while (entry) {
        if (strcmp(entry->name, name) == 0) {
            return entry->index;
        }
        entry = entry->next;
    }
    return -1;
}

// Add variable to hash table (returns index)
static slong var_hash_add(var_hash_table_t *table, const char *name) {
    // Check if already exists
    slong existing = var_hash_find(table, name);
    if (existing >= 0) return existing;
    
    // Add new entry
    slong bucket = hash_string(name, table->bucket_count);
    var_entry_t *entry = (var_entry_t*) malloc(sizeof(var_entry_t));
    if (!entry) return -1;
    
    entry->name = strdup(name);
    if (!entry->name) {
        free(entry);
        return -1;
    }
    entry->index = table->count++;
    entry->next = table->buckets[bucket];
    table->buckets[bucket] = entry;
    
    return entry->index;
}

static slong find_variable_optimized(poly_analysis_t *analysis, const char *var_name) {
    // For small numbers of variables, linear search is still fast
    if (analysis->num_all_vars < 16) {
        for (slong i = 0; i < analysis->num_all_vars; i++) {
            if (strcmp(analysis->all_vars[i], var_name) == 0) {
                return i;
            }
        }
        return -1;
    }
    
    // For larger numbers, linear search with optimizations
    for (slong i = 0; i < analysis->num_all_vars; i++) {
        const char *existing = analysis->all_vars[i];
        // Quick length check first
        if (strlen(existing) == strlen(var_name) && 
            strcmp(existing, var_name) == 0) {
            return i;
        }
    }
    return -1;
}

static int add_variable_optimized(poly_analysis_t *analysis, const char *var_name) {
    // Check if already exists using optimized search
    if (find_variable_optimized(analysis, var_name) >= 0) {
        return 1; // Already exists
    }
    
    // Expand array if needed
    if (analysis->num_all_vars >= analysis->max_vars) {
        slong new_max = analysis->max_vars + (analysis->max_vars >> 1) + 8;
        char **new_vars = (char**) realloc(analysis->all_vars, 
                                           new_max * sizeof(char*));
        if (!new_vars) {
            return 0; // Memory allocation failed
        }
        analysis->all_vars = new_vars;
        analysis->max_vars = new_max;
    }
    
    // Add the variable
    analysis->all_vars[analysis->num_all_vars] = strdup(var_name);
    if (!analysis->all_vars[analysis->num_all_vars]) {
        return 0; // strdup failed
    }
    
    analysis->num_all_vars++;
    return 1;
}

// Fast degree-only parsing without full polynomial construction

static int parse_and_extract_degree(lightweight_parser_t *parser) {
    parser->max_degree_found = 0;
    var_hash_init(&parser->var_table, 16);

    long current_term_degree = 0;   /* accumulated total degree of current monomial */
    int  in_term = 0;               /* are we inside a monomial? */

    while (parser->pos < parser->len) {
        char c = parser->input[parser->pos];

        /* ── whitespace / explicit multiply: stay in current monomial ── */
        if (isspace(c) || c == '*' || c == '(' || c == ')') {
            if (c == '(' || c == ')') {
                /* parentheses: flush current term first */
                if (in_term && current_term_degree > parser->max_degree_found)
                    parser->max_degree_found = current_term_degree;
                current_term_degree = 0;
                in_term = 0;
            }
            parser->pos++;
            continue;
        }

        /* ── additive operator: end of current monomial ── */
        if (c == '+' || c == '-') {
            if (in_term && current_term_degree > parser->max_degree_found)
                parser->max_degree_found = current_term_degree;
            current_term_degree = 0;
            in_term = 0;
            parser->pos++;
            continue;
        }

        /* ── numeric literal (coefficient): skip digits, don't add to degree ── */
        if (isdigit(c)) {
            while (parser->pos < parser->len &&
                   isdigit(parser->input[parser->pos]))
                parser->pos++;
            /* a bare number with no variable following is a constant monomial
               of degree 0 — mark that we are inside a term so a subsequent
               '+'/'-' triggers a flush */
            in_term = 1;
            continue;
        }

        /* ── identifier: variable or field generator ── */
        if (isalpha(c) || c == '_') {
            size_t start = parser->pos;
            while (parser->pos < parser->len &&
                   (isalnum(parser->input[parser->pos]) ||
                    parser->input[parser->pos] == '_'))
                parser->pos++;

            size_t id_len  = parser->pos - start;
            char  *var_name = malloc(id_len + 1);
            if (!var_name) return 0;
            strncpy(var_name, parser->input + start, id_len);
            var_name[id_len] = '\0';

            /* read optional exponent (applies to this identifier) */
            long exponent = 1;
            if (parser->pos < parser->len &&
                parser->input[parser->pos] == '^') {
                parser->pos++;          /* skip '^' */
                exponent = 0;
                while (parser->pos < parser->len &&
                       isdigit(parser->input[parser->pos])) {
                    exponent = exponent * 10 +
                               (parser->input[parser->pos] - '0');
                    parser->pos++;
                }
                if (exponent == 0) exponent = 1;
            }

            /* skip field generator — it contributes no "polynomial" degree */
            if (parser->generator_name &&
                strcmp(var_name, parser->generator_name) == 0) {
                free(var_name);
                in_term = 1;   /* still inside the monomial */
                continue;
            }

            /* register variable (for the caller's variable-discovery side-effect) */
            var_hash_add(&parser->var_table, var_name);
            free(var_name);

            /* accumulate total degree of this monomial */
            current_term_degree += exponent;
            in_term = 1;
            continue;
        }

        /* unknown character — skip */
        parser->pos++;
    }

    /* flush the last monomial */
    if (in_term && current_term_degree > parser->max_degree_found)
        parser->max_degree_found = current_term_degree;

    return 1;
}


/* =========================================================================
 * Complexity analysis helpers (internal to main)
 * ========================================================================= */

/*
 * Return the total degree of a polynomial string, ignoring the field
 * generator (gen_name).  Each monomial's contribution is the sum of
 * exponents of all non-generator identifiers.
 */
long get_poly_total_degree(const char *poly_str, const char *gen_name)
{
    if (!poly_str) return 0;

    long max_deg = 0, cur_deg = 0;
    int  in_term = 0;
    const char *p = poly_str;

    while (*p) {
        /* whitespace / multiply / parentheses */
        if (isspace((unsigned char)*p) || *p == '*' ||
            *p == '(' || *p == ')') {
            p++;
            continue;
        }

        /* additive operator → end current monomial */
        if (*p == '+' || *p == '-') {
            if (in_term && cur_deg > max_deg) max_deg = cur_deg;
            cur_deg = 0;
            in_term = 0;
            p++;
            continue;
        }

        /* numeric literal (coefficient) – skip digits */
        if (isdigit((unsigned char)*p)) {
            while (*p && isdigit((unsigned char)*p)) p++;
            in_term = 1;
            continue;
        }

        /* identifier */
        if (isalpha((unsigned char)*p) || *p == '_') {
            char name[64];
            int  ni = 0;
            while (*p && (isalnum((unsigned char)*p) || *p == '_') && ni < 63)
                name[ni++] = *p++;
            name[ni] = '\0';

            /* read exponent */
            long exp = 1;
            if (*p == '^') {
                p++;
                exp = 0;
                while (*p && isdigit((unsigned char)*p)) {
                    exp = exp * 10 + (*p - '0');
                    p++;
                }
                if (exp == 0) exp = 1;
            }

            /* skip field generator */
            if (gen_name && strcmp(name, gen_name) == 0) {
                in_term = 1;
                continue;
            }

            cur_deg += exp;
            in_term = 1;
            continue;
        }

        p++;
    }

    if (in_term && cur_deg > max_deg) max_deg = cur_deg;
    return max_deg;
}

/*
 * Collect unique variable names (excluding the field generator) that appear
 * across all polynomial strings.  Returns malloc'd array; caller frees.
 */
void collect_variables(const char **polys, slong npolys,
                               const char *gen_name,
                               char ***vars_out, slong *nvars_out)
{
    slong  cap  = 16;
    char **vars = malloc(cap * sizeof(char *));
    slong  nv   = 0;

    for (slong i = 0; i < npolys; i++) {
        const char *p = polys[i];
        while (*p) {
            if (isalpha((unsigned char)*p) || *p == '_') {
                char name[64];
                int  ni = 0;
                while (*p && (isalnum((unsigned char)*p) || *p == '_') && ni < 63)
                    name[ni++] = *p++;
                name[ni] = '\0';

                if (gen_name && strcmp(name, gen_name) == 0) continue;

                /* already recorded? */
                int found = 0;
                for (slong j = 0; j < nv; j++) {
                    if (strcmp(vars[j], name) == 0) { found = 1; break; }
                }
                if (!found) {
                    if (nv >= cap) {
                        cap *= 2;
                        vars = realloc(vars, cap * sizeof(char *));
                    }
                    vars[nv++] = strdup(name);
                }
            } else {
                p++;
            }
        }
    }

    *vars_out  = vars;
    *nvars_out = nv;
}


static void analyze_single_polynomial(poly_analysis_t *analysis, slong poly_idx, 
                                     const char *poly_str) {
    // Input validation
    if (!analysis || !poly_str || poly_idx >= analysis->num_polys) {
        if (analysis && poly_idx < analysis->num_polys) {
            analysis->degrees[poly_idx] = 0;
        }
        return;
    }
    
    // Early exit for empty polynomial
    size_t poly_len = strlen(poly_str);
    if (poly_len == 0) {
        analysis->degrees[poly_idx] = 0;
        return;
    }
    
    // Get generator name
    char *gen_name = get_generator_name(analysis->ctx);
    if (!gen_name) {
        analysis->degrees[poly_idx] = 0;
        return;
    }
    
    // Use lightweight parser for degree calculation
    lightweight_parser_t light_parser;
    light_parser.input = poly_str;
    light_parser.pos = 0;
    light_parser.len = poly_len;
    light_parser.ctx = analysis->ctx;
    light_parser.generator_name = gen_name;
    
    int parse_success = parse_and_extract_degree(&light_parser);
    
    if (parse_success) {
        // Extract variables from the hash table and add to global analysis
        for (slong i = 0; i < light_parser.var_table.bucket_count; i++) {
            var_entry_t *entry = light_parser.var_table.buckets[i];
            while (entry) {
                if (!add_variable_optimized(analysis, entry->name)) {
                    break;
                }
                entry = entry->next;
            }
        }
        
        analysis->degrees[poly_idx] = light_parser.max_degree_found;
    } else {
        analysis->degrees[poly_idx] = 0;
    }
    
    // Cleanup
    var_hash_clear(&light_parser.var_table);
    free(gen_name);
}

// Check if a variable is in the elimination list
static int is_elimination_var(const char *var_name, const char **elim_vars, slong num_elim_vars) {
    for (slong i = 0; i < num_elim_vars; i++) {
        if (strcmp(var_name, elim_vars[i]) == 0) {
            return 1;
        }
    }
    return 0;
}

// Main function: dixon_complexity_auto with complexity encoding
char* dixon_complexity_auto(const char **poly_strings, slong num_polys,
                           const char **elim_vars, slong num_elim_vars,
                           const fq_nmod_ctx_t ctx) {
    
    printf("\n=== Dixon Complexity Analysis ===\n");
    
    // Initialize analysis structure
    poly_analysis_t analysis;
    poly_analysis_init(&analysis, num_polys, ctx);
    
    // Analyze each polynomial
    for (slong i = 0; i < num_polys; i++) {
        analyze_single_polynomial(&analysis, i, poly_strings[i]);
    }
    
    // Identify remaining variables
    char **remaining_vars = (char**) malloc(analysis.num_all_vars * sizeof(char*));
    slong num_remaining = 0;
    
    for (slong i = 0; i < analysis.num_all_vars; i++) {
        if (!is_elimination_var(analysis.all_vars[i], elim_vars, num_elim_vars)) {
            remaining_vars[num_remaining] = strdup(analysis.all_vars[i]);
            num_remaining++;
        }
    }
    
    // Calculate total degree product
    long td = 1;
    for (slong i = 0; i < num_polys; i++) {
        td *= analysis.degrees[i];
    }

    dixon_complexity_report_t report;
    dixon_complexity_report_from_degrees(&report,
                                         analysis.degrees,
                                         num_polys,
                                         analysis.num_all_vars,
                                         num_elim_vars,
                                         num_remaining,
                                         DIXON_OMEGA);
    
    fmpz_t matrix_size;
    fmpz_init(matrix_size);
    dixon_size(matrix_size, analysis.degrees, num_polys, 0);

    dixon_complexity_write_report_body(stdout,
                                       num_polys,
                                       analysis.num_all_vars,
                                       analysis.all_vars,
                                       (char **) elim_vars,
                                       num_elim_vars,
                                       num_remaining,
                                       analysis.degrees,
                                       matrix_size,
                                       td,
                                       &report,
                                       DIXON_OMEGA);
    
    // Encode a single compatibility value for callers that expect one scalar.
    double complexity = dixon_complexity_best_total_log2(&report);
    long complexity_encoded = (long)(complexity * 1000.0 + 0.5);
    
    // Generate evaluation polynomial
    printf("\n=== Generating Evaluation Polynomial ===\n");
    
    size_t eval_poly_size = 1024;
    char *eval_poly = (char*) malloc(eval_poly_size);
    eval_poly[0] = '\0';
    
    if (num_remaining == 0) {
        sprintf(eval_poly, "%ld", complexity_encoded);
        printf("No remaining variables, evaluation polynomial: %s\n", eval_poly);
    } else {
        for (slong i = 0; i < num_remaining; i++) {
            if (i > 0) {
                strcat(eval_poly, " + ");
            }
            
            size_t needed = strlen(eval_poly) + strlen(remaining_vars[i]) + 50;
            if (needed >= eval_poly_size) {
                eval_poly_size = needed * 2;
                eval_poly = (char*) realloc(eval_poly, eval_poly_size);
            }
            
            strcat(eval_poly, remaining_vars[i]);
            if (td > 1) {
                char exp_str[32];
                sprintf(exp_str, "^%ld", td);
                strcat(eval_poly, exp_str);
            }
        }
        
        char complexity_str[64];
        sprintf(complexity_str, " + %ld", complexity_encoded);
        
        size_t needed = strlen(eval_poly) + strlen(complexity_str) + 1;
        if (needed >= eval_poly_size) {
            eval_poly_size = needed * 2;
            eval_poly = (char*) realloc(eval_poly, eval_poly_size);
        }
        
        strcat(eval_poly, complexity_str);
        printf("Evaluation polynomial: %s\n", eval_poly);
    }
    
    // Cleanup
    for (slong i = 0; i < num_remaining; i++) {
        free(remaining_vars[i]);
    }
    free(remaining_vars);
    poly_analysis_clear(&analysis);
    fmpz_clear(matrix_size);
    
    printf("=== Analysis Complete ===\n\n");
    
    return eval_poly;
}

// String interface version
char* dixon_complexity_auto_str(const char *poly_string,
                                const char *vars_string,
                                const fq_nmod_ctx_t ctx) {
    
    // Split input strings
    slong num_polys, num_vars;
    char **poly_array = split_string(poly_string, &num_polys);
    char **vars_array = split_string(vars_string, &num_vars);
    
    // Convert to const char**
    const char **poly_strings = (const char**) malloc(num_polys * sizeof(char*));
    const char **elim_vars = (const char**) malloc(num_vars * sizeof(char*));
    
    for (slong i = 0; i < num_polys; i++) {
        poly_strings[i] = poly_array[i];
    }
    for (slong i = 0; i < num_vars; i++) {
        elim_vars[i] = vars_array[i];
    }
    
    // Call main function
    char *result = dixon_complexity_auto(poly_strings, num_polys, elim_vars, num_vars, ctx);
    
    // Cleanup
    free_split_strings(poly_array, num_polys);
    free_split_strings(vars_array, num_vars);
    free(poly_strings);
    free(elim_vars);
    
    return result;
}

// Extract constant term from polynomial string
static long extract_constant_term(const char *poly_str) {
    if (!poly_str || strlen(poly_str) == 0) {
        return 0;
    }
    
    // Simple case: entire string is a number
    char *endptr;
    long value = strtol(poly_str, &endptr, 10);
    if (*endptr == '\0') {
        return value;
    }
    
    // Complex case: find constant terms
    long max_constant = 0;
    const char *ptr = poly_str;
    
    while (*ptr) {
        while (*ptr && isspace(*ptr)) ptr++;
        if (!*ptr) break;
        
        if (*ptr == '+') {
            ptr++;
            while (*ptr && isspace(*ptr)) ptr++;
        }
        
        if (isdigit(*ptr) || (*ptr == '-' && isdigit(*(ptr+1)))) {
            char *term_end;
            long term_value = strtol(ptr, &term_end, 10);
            
            const char *check_ptr = term_end;
            while (*check_ptr && isspace(*check_ptr)) check_ptr++;
            
            if (!*check_ptr || *check_ptr == '+' || *check_ptr == '-') {
                if (term_value > max_constant) {
                    max_constant = term_value;
                }
            }
            
            ptr = term_end;
        } else {
            while (*ptr && *ptr != '+' && *ptr != '-') ptr++;
            if (*ptr == '-') ptr++;
        }
    }
    
    return max_constant;
}

// Extract maximum complexity from multiple polynomials
double extract_max_complexity(const char **poly_strings, slong num_polys) {
    if (!poly_strings || num_polys == 0) {
        return 0.0;
    }
    
    long max_encoded = 0;
    
    printf("\n=== Extracting Complexity Values ===\n");
    
    for (slong i = 0; i < num_polys; i++) {
        long constant = extract_constant_term(poly_strings[i]);
        printf("Polynomial %ld: constant term = %ld\n", i + 1, constant);
        
        if (constant > max_encoded) {
            max_encoded = constant;
        }
    }
    
    double complexity = (double)max_encoded / 1000.0;
    printf("Maximum encoded value: %ld\n", max_encoded);
    printf("Decoded complexity: %.6f\n", complexity);
    printf("=== Extraction Complete ===\n\n");
    
    return complexity;
}

// String interface for complexity extraction
double extract_max_complexity_str(const char *poly_string) {
    if (!poly_string || strlen(poly_string) == 0) {
        return 0.0;
    }
    
    slong num_polys;
    char **poly_array = split_string(poly_string, &num_polys);
    
    if (!poly_array || num_polys == 0) {
        return 0.0;
    }
    
    const char **poly_strings = (const char**) malloc(num_polys * sizeof(char*));
    for (slong i = 0; i < num_polys; i++) {
        poly_strings[i] = poly_array[i];
    }
    
    double complexity = extract_max_complexity(poly_strings, num_polys);
    
    free_split_strings(poly_array, num_polys);
    free(poly_strings);
    
    return complexity;
}

int test_dixon_complexity() {
    // Test data
    long a1[] = {1000, 1000, 1000, 1001, 1002, 1003};
    int len = sizeof(a1) / sizeof(a1[0]);
    double omega = DIXON_OMEGA;
    
    printf("Dixon Complexity Results (Hessenberg Method):\n");
    for (int n = 5; n < 10; n++) {
        double complexity = dixon_complexity(a1, len, n, omega);
        printf("n=%d: %.6f\n", n, complexity);
    }
    
    // Test dixon_size function
    printf("\nTesting dixon_size with Hessenberg method:\n");
    fmpz_t test_result;
    fmpz_init(test_result);
    dixon_size(test_result, a1, len, 1);
    fmpz_clear(test_result);
    
    return 0;
}
