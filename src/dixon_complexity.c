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
    
    printf("\n=== Dixon Complexity Analysis (Hessenberg Method) ===\n");
    
    // Initialize analysis structure
    poly_analysis_t analysis;
    poly_analysis_init(&analysis, num_polys, ctx);
    
    // Analyze each polynomial
    printf("Analyzing %ld polynomials...\n", num_polys);
    for (slong i = 0; i < num_polys; i++) {
        printf("  Polynomial %ld: ", i + 1);
        analyze_single_polynomial(&analysis, i, poly_strings[i]);
        printf("degree = %ld\n", analysis.degrees[i]);
    }
    
    // Print all discovered variables
    printf("\nDiscovered variables (%ld): ", analysis.num_all_vars);
    for (slong i = 0; i < analysis.num_all_vars; i++) {
        if (i > 0) printf(", ");
        printf("%s", analysis.all_vars[i]);
    }
    printf("\n");
    
    printf("Elimination variables (%ld): ", num_elim_vars);
    for (slong i = 0; i < num_elim_vars; i++) {
        if (i > 0) printf(", ");
        printf("%s", elim_vars[i]);
    }
    printf("\n");
    
    // Identify remaining variables
    char **remaining_vars = (char**) malloc(analysis.num_all_vars * sizeof(char*));
    slong num_remaining = 0;
    
    for (slong i = 0; i < analysis.num_all_vars; i++) {
        if (!is_elimination_var(analysis.all_vars[i], elim_vars, num_elim_vars)) {
            remaining_vars[num_remaining] = strdup(analysis.all_vars[i]);
            num_remaining++;
        }
    }
    
    printf("Remaining variables (%ld): ", num_remaining);
    for (slong i = 0; i < num_remaining; i++) {
        if (i > 0) printf(", ");
        printf("%s", remaining_vars[i]);
    }
    printf("\n");
    
    // Calculate total degree product
    long td = 1;
    printf("\nDegree calculation: ");
    for (slong i = 0; i < num_polys; i++) {
        if (i > 0) printf(" * ");
        printf("%ld", analysis.degrees[i]);
        td *= analysis.degrees[i];
    }
    printf(" = %ld\n", td);
    
    // Compute Dixon matrix size and complexity using Hessenberg method
    printf("\n=== Dixon Matrix Analysis (Hessenberg Method) ===\n");
    
    fmpz_t matrix_size;
    fmpz_init(matrix_size);
    dixon_size(matrix_size, analysis.degrees, num_polys, 1);
    
    printf("Dixon matrix size: ");
    fmpz_print(matrix_size);
    printf("\n");
    
    // Calculate complexity
    double complexity = dixon_complexity(analysis.degrees, num_polys, 
                                       analysis.num_all_vars, DIXON_OMEGA);
    printf("Dixon complexity (log_2): %.6f\n", complexity);
    printf("Omega parameter: %.1f\n", DIXON_OMEGA);
    
    // Encode complexity as integer
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
    double omega = 2.3;
    
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
