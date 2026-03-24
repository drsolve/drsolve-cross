#include "polynomial_system_solver.h"
#include <ctype.h>
#include <stdarg.h>

static int solver_realtime_progress_enabled = 0;

void polynomial_solver_set_realtime_progress(int enabled) {
    solver_realtime_progress_enabled = enabled;
}

static void solver_progress(const char *fmt, ...) {
    if (!solver_realtime_progress_enabled) {
        return;
    }

    va_list args;
    va_start(args, fmt);
    fprintf(stderr, "Progress: ");
    vfprintf(stderr, fmt, args);
    fprintf(stderr, "\n");
    fflush(stderr);
    va_end(args);
}

static char *join_variable_names(variable_info_t *vars, slong num_vars) {
    if (!vars || num_vars <= 0) {
        return NULL;
    }

    size_t total_len = 1;
    for (slong i = 0; i < num_vars; i++) {
        total_len += strlen(vars[i].name) + 2;
    }

    char *result = (char *) malloc(total_len);
    if (!result) {
        return NULL;
    }

    result[0] = '\0';
    for (slong i = 0; i < num_vars; i++) {
        if (i > 0) {
            strcat(result, ", ");
        }
        strcat(result, vars[i].name);
    }

    return result;
}

static char *format_equation_index_list(const slong *indices, slong count) {
    if (!indices || count <= 0) {
        return NULL;
    }

    size_t total_len = 1;
    for (slong i = 0; i < count; i++) {
        total_len += 24;
    }

    char *result = (char *) malloc(total_len);
    if (!result) {
        return NULL;
    }

    result[0] = '\0';
    for (slong i = 0; i < count; i++) {
        char piece[32];
        if (i > 0) {
            strcat(result, ",");
        }
        snprintf(piece, sizeof(piece), "%ld", indices[i] + 1);
        strcat(result, piece);
    }

    return result;
}

static void append_string_entry(char ***items, slong *count, slong *cap, char *entry) {
    if (!entry) {
        return;
    }

    if (*count >= *cap) {
        slong new_cap = (*cap == 0) ? 8 : (*cap * 2);
        char **new_items = (char **) realloc(*items, new_cap * sizeof(char *));
        if (!new_items) {
            free(entry);
            return;
        }
        *items = new_items;
        *cap = new_cap;
    }

    (*items)[*count] = entry;
    (*count)++;
}

static void store_resultant_step(polynomial_solutions_t *sols, char *buffer) {
    if (!sols || !buffer) {
        free(buffer);
        return;
    }

    append_string_entry(&sols->resultant_steps,
                        &sols->num_resultant_steps,
                        &sols->resultant_steps_cap,
                        buffer);
}

static void add_resultant_step(polynomial_solutions_t *sols, const char *fmt, ...) {
    if (!sols) {
        return;
    }

    va_list args;
    va_list args_copy;
    va_start(args, fmt);
    va_copy(args_copy, args);
    int needed = vsnprintf(NULL, 0, fmt, args_copy);
    va_end(args_copy);

    if (needed < 0) {
        va_end(args);
        return;
    }

    char *buffer = (char *) malloc((size_t) needed + 1);
    if (!buffer) {
        va_end(args);
        return;
    }

    vsnprintf(buffer, (size_t) needed + 1, fmt, args);
    va_end(args);

    store_resultant_step(sols, buffer);
    solver_progress("%s", buffer);
}

static char *format_solution_set_line_from_arrays(char **variable_names, slong num_variables,
                                                  fq_nmod_t **solution_set, slong *solutions_per_var,
                                                  const fq_nmod_ctx_t ctx) {
    if (!variable_names || !solution_set || !solutions_per_var || num_variables <= 0) {
        return NULL;
    }

    size_t total_len = 1;
    for (slong var = 0; var < num_variables; var++) {
        total_len += strlen(variable_names[var]) + 16;
        slong num_sols = solutions_per_var[var];
        if (num_sols == 0) {
            total_len += strlen("no solution");
        } else {
            for (slong sol = 0; sol < num_sols; sol++) {
                char *sol_str = fq_nmod_get_str_pretty(solution_set[var][sol], ctx);
                if (sol_str) {
                    total_len += strlen(sol_str) + 2;
                    free(sol_str);
                }
            }
        }
    }

    char *line = (char *) malloc(total_len);
    if (!line) {
        return NULL;
    }

    line[0] = '\0';
    for (slong var = 0; var < num_variables; var++) {
        slong num_sols = solutions_per_var[var];

        if (var > 0) {
            strcat(line, ", ");
        }
        strcat(line, variable_names[var]);
        strcat(line, " = ");

        if (num_sols == 0) {
            strcat(line, "no solution");
        } else if (num_sols == 1) {
            char *sol_str = fq_nmod_get_str_pretty(solution_set[var][0], ctx);
            if (sol_str) {
                strcat(line, sol_str);
                free(sol_str);
            }
        } else {
            strcat(line, "{");
            for (slong sol = 0; sol < num_sols; sol++) {
                char *sol_str = fq_nmod_get_str_pretty(solution_set[var][sol], ctx);
                if (sol > 0) {
                    strcat(line, ", ");
                }
                if (sol_str) {
                    strcat(line, sol_str);
                    free(sol_str);
                }
            }
            strcat(line, "}");
        }
    }

    return line;
}

static void add_candidate_solution(polynomial_solutions_t *sols, char *line, int pass) {
    if (!sols || !line) {
        return;
    }

    if (sols->num_candidate_solution_lines >= sols->candidate_solution_lines_cap) {
        slong new_cap = (sols->candidate_solution_lines_cap == 0)
                        ? 8 : (sols->candidate_solution_lines_cap * 2);
        char **new_lines = (char **) realloc(sols->candidate_solution_lines,
                                             new_cap * sizeof(char *));
        int *new_pass = (int *) realloc(sols->candidate_solution_pass,
                                        new_cap * sizeof(int));
        if (!new_lines || !new_pass) {
            free(line);
            if (new_lines) sols->candidate_solution_lines = new_lines;
            if (new_pass) sols->candidate_solution_pass = new_pass;
            return;
        }
        sols->candidate_solution_lines = new_lines;
        sols->candidate_solution_pass = new_pass;
        sols->candidate_solution_lines_cap = new_cap;
    }

    sols->candidate_solution_lines[sols->num_candidate_solution_lines] = line;
    sols->candidate_solution_pass[sols->num_candidate_solution_lines] = pass;
    sols->num_candidate_solution_lines++;
}

static void merge_solver_logs(polynomial_solutions_t *dst, const polynomial_solutions_t *src) {
    if (!dst || !src) {
        return;
    }

    for (slong i = 0; i < src->num_resultant_steps; i++) {
        store_resultant_step(dst, strdup(src->resultant_steps[i]));
    }
}

// ============= BASIC UTILITY FUNCTIONS =============

// Comparison function: sort by maximum degree in ascending order
static int compare_variables_by_degree(const void *a, const void *b) {
    const variable_info_t *va = (const variable_info_t*)a;
    const variable_info_t *vb = (const variable_info_t*)b;
    
    if (va->max_degree < vb->max_degree) return -1;
    if (va->max_degree > vb->max_degree) return 1;
    return strcmp(va->name, vb->name); // Sort by name if degrees are equal
}

// Check if a polynomial string contains a specific variable
int contains_variable(const char *poly_str, const char *var_name) {
    const char *ptr = poly_str;
    size_t var_len = strlen(var_name);
    
    while ((ptr = strstr(ptr, var_name)) != NULL) {
        // Check if this is a complete variable name, not part of another identifier
        if ((ptr == poly_str || !isalnum(*(ptr-1))) && 
            !isalnum(*(ptr + var_len))) {
            return 1;
        }
        ptr++;
    }
    return 0;
}

// Get extension field generator name (reference DIXON_INTERFACE_FLINT_H)
char* get_generator_name_for_solver(const fq_nmod_ctx_t ctx) {
    if (fq_nmod_ctx_degree(ctx) > 1) {
        if (ctx->var && strlen(ctx->var) > 0) {
            return strdup(ctx->var);
        } else {
            return strdup("t");  // Default extension generator name
        }
    } else {
        return NULL;  // Prime field doesn't need generator
    }
}

// Check if string contains specific identifier (improved version considering extension generator)
int contains_identifier_improved(const char *poly_str, const char *identifier, const fq_nmod_ctx_t ctx) {
    const char *ptr = poly_str;
    size_t id_len = strlen(identifier);
    
    // Get extension generator name
    char *gen_name = get_generator_name_for_solver(ctx);
    
    while ((ptr = strstr(ptr, identifier)) != NULL) {
        // Check if this is a complete identifier name, not part of another identifier
        if ((ptr == poly_str || !isalnum(*(ptr-1))) && 
            !isalnum(*(ptr + id_len))) {
            
            // Special handling for extension generator
            if (gen_name && strcmp(identifier, gen_name) == 0) {
                free(gen_name);
                return 1;  // Found extension generator
            }
            
            free(gen_name);
            return 1;  // Found regular variable
        }
        ptr++;
    }
    
    if (gen_name) free(gen_name);
    return 0;
}

// Improved variable extraction function
char** extract_variables_improved(char **poly_strings, slong num_polys, slong *num_vars_out, const fq_nmod_ctx_t ctx) {
    // Temporary storage for discovered variables
    char **temp_vars = (char**) malloc(100 * sizeof(char*));  // Maximum 100 variables
    slong temp_count = 0;
    
    // Get extension generator name
    char *gen_name = get_generator_name_for_solver(ctx);
    
    printf("=== Improved Variable Extraction ===\n");
    if (gen_name) {
        printf("Extension field generator: %s\n", gen_name);
    } else {
        printf("Prime field - no generator\n");
    }
    
    // Iterate through all polynomial strings
    for (slong poly_idx = 0; poly_idx < num_polys; poly_idx++) {
        const char *poly = poly_strings[poly_idx];
        size_t len = strlen(poly);
        
        for (size_t i = 0; i < len; i++) {
            if (isalpha(poly[i]) || poly[i] == '_') {
                // Extract complete identifier
                size_t start = i;
                while (i < len && (isalnum(poly[i]) || poly[i] == '_')) {
                    i++;
                }
                
                // Create identifier string
                size_t id_len = i - start;
                char *identifier = (char*) malloc(id_len + 1);
                strncpy(identifier, poly + start, id_len);
                identifier[id_len] = '\0';
                
                // Skip numbers (not variables)
                if (isdigit(identifier[0])) {
                    free(identifier);
                    continue;
                }
                
                // Check if this is extension generator
                int is_generator = 0;
                if (gen_name && strcmp(identifier, gen_name) == 0) {
                    is_generator = 1;
                    //printf("Found generator in polynomial %ld: %s\n", poly_idx, identifier);
                }
                
                // Check if already in list (but extension generator doesn't count as variable)
                int already_exists = 0;
                if (!is_generator) {
                    for (slong j = 0; j < temp_count; j++) {
                        if (strcmp(temp_vars[j], identifier) == 0) {
                            already_exists = 1;
                            break;
                        }
                    }
                    
                    // If it's a new variable, add to list
                    if (!already_exists) {
                        temp_vars[temp_count] = strdup(identifier);
                        temp_count++;
                        printf("Found variable: %s\n", identifier);
                    }
                }
                
                free(identifier);
                i--; // Because for loop will i++ again
            }
        }
    }
    
    if (gen_name) free(gen_name);
    
    // Create final variable array
    char **final_vars = NULL;
    if (temp_count > 0) {
        final_vars = (char**) malloc(temp_count * sizeof(char*));
        for (slong i = 0; i < temp_count; i++) {
            final_vars[i] = temp_vars[i];
        }
    }
    
    free(temp_vars);
    *num_vars_out = temp_count;
    
    printf("Total variables found: %ld\n", temp_count);
    for (slong i = 0; i < temp_count; i++) {
        printf("  Variable %ld: %s\n", i, final_vars[i]);
    }
    
    return final_vars;
}

// Improved variable degree calculation function
slong get_variable_max_degree_in_polynomial(const char *poly_str, const char *var_name, const fq_nmod_ctx_t ctx) {
    const char *ptr = poly_str;
    slong max_degree = 0;
    size_t var_len = strlen(var_name);
    
    // Get extension generator name to ensure no confusion with extension generator
    char *gen_name = get_generator_name_for_solver(ctx);
    int is_generator = (gen_name && strcmp(var_name, gen_name) == 0);
    
    while ((ptr = strstr(ptr, var_name)) != NULL) {
        // Ensure this is a complete variable name
        if ((ptr == poly_str || !isalnum(*(ptr-1))) && 
            !isalnum(*(ptr + var_len))) {
            
            const char *exp_ptr = ptr + var_len;
            while (isspace(*exp_ptr)) exp_ptr++;
            
            slong degree = 1;  // Default degree is 1
            if (*exp_ptr == '^') {
                exp_ptr++;
                while (isspace(*exp_ptr)) exp_ptr++;
                if (isdigit(*exp_ptr)) {
                    degree = strtol(exp_ptr, NULL, 10);
                }
            }
            
            if (degree > max_degree) {
                max_degree = degree;
            }
        }
        ptr++;
    }
    
    if (gen_name) free(gen_name);
    return max_degree;
}

// ============= EQUATION COMBINATION MANAGEMENT =============

// Generate all possible combinations of n equations taken k at a time
void generate_equation_combinations(slong num_equations, slong target_equations, 
                                   equation_combination_t **combinations, slong *num_combinations) {
    // Calculate number of combinations: C(n,k)
    slong total_combinations = 1;
    for (slong i = 0; i < target_equations; i++) {
        total_combinations = total_combinations * (num_equations - i) / (i + 1);
    }
    
    *num_combinations = total_combinations;
    *combinations = (equation_combination_t*) malloc(total_combinations * sizeof(equation_combination_t));
    
    // Generate combinations using lexicographic order
    slong comb_idx = 0;
    slong *current_combination = (slong*) malloc(target_equations * sizeof(slong));
    
    // Initialize first combination: [0, 1, 2, ..., target_equations-1]
    for (slong i = 0; i < target_equations; i++) {
        current_combination[i] = i;
    }
    
    while (comb_idx < total_combinations) {
        // Store current combination
        (*combinations)[comb_idx].equation_indices = (slong*) malloc(target_equations * sizeof(slong));
        (*combinations)[comb_idx].num_equations = target_equations;
        (*combinations)[comb_idx].tried = 0;
        (*combinations)[comb_idx].success = 0;
        
        for (slong i = 0; i < target_equations; i++) {
            (*combinations)[comb_idx].equation_indices[i] = current_combination[i];
        }
        comb_idx++;
        
        // Generate next combination
        slong k = target_equations - 1;
        while (k >= 0 && current_combination[k] == num_equations - target_equations + k) {
            k--;
        }
        
        if (k < 0) break; // All combinations generated
        
        current_combination[k]++;
        for (slong i = k + 1; i < target_equations; i++) {
            current_combination[i] = current_combination[i-1] + 1;
        }
    }
    
    free(current_combination);
}

// Free equation combinations
void free_equation_combinations(equation_combination_t *combinations, slong num_combinations) {
    if (combinations) {
        for (slong i = 0; i < num_combinations; i++) {
            if (combinations[i].equation_indices) {
                free(combinations[i].equation_indices);
            }
        }
        free(combinations);
    }
}

// ============= SOLUTION STRUCTURE MANAGEMENT =============

// Initialize the enhanced solution structure
void polynomial_solutions_init(polynomial_solutions_t *sols, slong num_vars, 
                               const fq_nmod_ctx_t ctx) {
    sols->num_variables = num_vars;
    sols->ctx = ctx;
    sols->is_valid = 0;
    sols->has_no_solutions = 0;
    sols->error_message = NULL;
    sols->num_solution_sets = 0;
    sols->num_equations = 0;
    sols->variable_order = NULL;
    sols->elimination_summary = NULL;
    sols->total_combinations = 0;
    sols->successful_combinations = 0;
    sols->num_base_solutions = 0;
    sols->checked_solution_sets = -1;
    sols->verified_solution_sets = -1;
    sols->resultant_steps = NULL;
    sols->num_resultant_steps = 0;
    sols->resultant_steps_cap = 0;
    sols->candidate_solution_lines = NULL;
    sols->candidate_solution_pass = NULL;
    sols->num_candidate_solution_lines = 0;
    sols->candidate_solution_lines_cap = 0;
    
    if (num_vars > 0) {
        sols->variable_names = (char**) calloc(num_vars, sizeof(char*));
        sols->solution_sets = NULL;  // Will be allocated when we know the number of solution sets
        sols->solutions_per_var = NULL;
    } else {
        sols->variable_names = NULL;
        sols->solution_sets = NULL;
        sols->solutions_per_var = NULL;
    }
}

// Clear the enhanced solution structure
void polynomial_solutions_clear(polynomial_solutions_t *sols) {
    if (sols->variable_names) {
        for (slong i = 0; i < sols->num_variables; i++) {
            if (sols->variable_names[i]) {
                free(sols->variable_names[i]);
            }
        }
        free(sols->variable_names);
    }
    
    if (sols->solution_sets) {
        for (slong set = 0; set < sols->num_solution_sets; set++) {
            if (sols->solution_sets[set]) {
                for (slong var = 0; var < sols->num_variables; var++) {
                    if (sols->solution_sets[set][var]) {
                        for (slong sol = 0; sol < sols->solutions_per_var[set * sols->num_variables + var]; sol++) {
                            fq_nmod_clear(sols->solution_sets[set][var][sol], sols->ctx);
                        }
                        free(sols->solution_sets[set][var]);
                    }
                }
                free(sols->solution_sets[set]);
            }
        }
        free(sols->solution_sets);
    }
    
    if (sols->solutions_per_var) {
        free(sols->solutions_per_var);
    }
    
    if (sols->error_message) {
        free(sols->error_message);
    }

    if (sols->variable_order) {
        free(sols->variable_order);
    }

    if (sols->elimination_summary) {
        free(sols->elimination_summary);
    }

    if (sols->resultant_steps) {
        for (slong i = 0; i < sols->num_resultant_steps; i++) {
            free(sols->resultant_steps[i]);
        }
        free(sols->resultant_steps);
    }

    if (sols->candidate_solution_lines) {
        for (slong i = 0; i < sols->num_candidate_solution_lines; i++) {
            free(sols->candidate_solution_lines[i]);
        }
        free(sols->candidate_solution_lines);
    }

    if (sols->candidate_solution_pass) {
        free(sols->candidate_solution_pass);
    }
    
    memset(sols, 0, sizeof(polynomial_solutions_t));
}

// Enhanced print solutions function
void print_polynomial_solutions(const polynomial_solutions_t *sols) {
    if (!sols) {
        printf("No solution data available.\n");
        return;
    }

    if (sols->num_resultant_steps > 0) {
        for (slong i = 0; i < sols->num_resultant_steps; i++) {
            printf("Progress: %s\n", sols->resultant_steps[i]);
        }
        printf("\n");
    }

    if (!sols->is_valid) {
        printf("Solve failed");
        if (sols->error_message) {
            printf(": %s", sols->error_message);
        }
        printf("\n");
        return;
    }

    if (sols->has_no_solutions == -1) {
        printf("System has positive dimension; please use Dixon resultant elimination.\n");
        return;
    }

    if (sols->has_no_solutions == 1) {
        printf("No solutions over the finite field.\n");
        return;
    }

    if (sols->num_variables == 0) {
        printf("No variables.\n");
        return;
    }

    if (sols->num_solution_sets == 0) {
        printf("No solutions found.\n");
        return;
    }

    if (sols->num_candidate_solution_lines > 0) {
        printf("Candidate sets:\n");
        for (slong i = 0; i < sols->num_candidate_solution_lines; i++) {
            printf("  [%ld] %s [%s]\n",
                   i + 1,
                   sols->candidate_solution_lines[i],
                   sols->candidate_solution_pass[i] ? "PASS" : "FAIL");
        }
    }
    printf("Found %ld solution set(s):\n", sols->num_solution_sets);

    for (slong set = 0; set < sols->num_solution_sets; set++) {
        printf("  [%ld] ", set + 1);
        for (slong var = 0; var < sols->num_variables; var++) {
            slong num_sols = sols->solutions_per_var[set * sols->num_variables + var];

            if (var > 0) {
                printf(", ");
            }

            printf("%s = ", sols->variable_names[var]);
            if (num_sols == 0) {
                printf("no solution");
            } else if (num_sols == 1) {
                fq_nmod_print_pretty(sols->solution_sets[set][var][0], sols->ctx);
            } else {
                printf("{");
                for (slong sol = 0; sol < num_sols; sol++) {
                    if (sol > 0) printf(", ");
                    fq_nmod_print_pretty(sols->solution_sets[set][var][sol], sols->ctx);
                }
                printf("}");
            }
        }
        printf("\n");
    }
}

// ============= VARIABLE ANALYSIS FUNCTIONS =============

// Extract variables and sort by degree
int extract_and_sort_variables(char **poly_strings, slong num_polys,
                                       variable_info_t **sorted_vars_out, slong *num_vars_out,
                                       const fq_nmod_ctx_t ctx) {
    
    printf("=== Analyzing polynomial variables (improved) ===\n");
    
    // Use improved variable extraction
    char **all_vars = extract_variables_improved(poly_strings, num_polys, num_vars_out, ctx);
    
    if (*num_vars_out == 0) {
        *sorted_vars_out = NULL;
        return 1;  // No variables is also successful
    }
    
    // Analyze degree of each polynomial
    for (slong i = 0; i < num_polys; i++) {
        slong total_degree = 0;
        
        // Calculate total degree of this polynomial
        for (slong j = 0; j < *num_vars_out; j++) {
            slong var_degree = get_variable_max_degree_in_polynomial(poly_strings[i], all_vars[j], ctx);
            if (var_degree > total_degree) {
                total_degree = var_degree;  // Simplified: use highest variable degree
            }
        }
        
        printf("Polynomial %ld: estimated degree = %ld\n", i + 1, total_degree);
    }
    
    printf("Discovered variables (%ld): ", *num_vars_out);
    for (slong i = 0; i < *num_vars_out; i++) {
        if (i > 0) printf(", ");
        printf("%s", all_vars[i]);
    }
    printf("\n");
    
    // Create variable info array
    variable_info_t *var_info = (variable_info_t*) malloc(*num_vars_out * sizeof(variable_info_t));
    
    // Calculate maximum degree of each variable across all polynomials
    for (slong i = 0; i < *num_vars_out; i++) {
        var_info[i].name = strdup(all_vars[i]);
        var_info[i].index = i;
        var_info[i].max_degree = 0;
        
        // Find maximum degree of this variable across all polynomials
        for (slong j = 0; j < num_polys; j++) {
            slong degree = get_variable_max_degree_in_polynomial(poly_strings[j], all_vars[i], ctx);
            if (degree > var_info[i].max_degree) {
                var_info[i].max_degree = degree;
            }
        }
        
        free(all_vars[i]);  // Free temporary string
    }
    
    free(all_vars);  // Free temporary array
    
    // Sort by degree
    qsort(var_info, *num_vars_out, sizeof(variable_info_t), compare_variables_by_degree);
    
    printf("\nVariables sorted by maximum degree:\n");
    for (slong i = 0; i < *num_vars_out; i++) {
        printf("  %s: max degree %ld\n", var_info[i].name, var_info[i].max_degree);
    }
    
    *sorted_vars_out = var_info;
    return 1;
}

// ============= POLYNOMIAL MANIPULATION FUNCTIONS =============

// Solve univariate equation and return all roots
fq_nmod_t* solve_univariate_equation_all_roots(const char *poly_str, const char *var_name,
                                                       slong *num_roots_out, const fq_nmod_ctx_t ctx) {
    
    // Use parsing mechanism from DIXON_INTERFACE_FLINT_H
    char *gen_name = get_generator_name_for_solver(ctx);
    
    parser_state_t state = {0};
    state.var_names = (char**) malloc(sizeof(char*));
    state.var_names[0] = strdup(var_name);
    state.nvars = 1;
    state.npars = 0;
    state.max_pars = 16;
    state.par_names = (char**) malloc(state.max_pars * sizeof(char*));
    state.ctx = ctx;
    state.current.str = NULL;
    fq_nmod_init(state.current.value, ctx);
    state.generator_name = gen_name;
    
    // Parse polynomial
    fq_mvpoly_t poly;
    fq_mvpoly_init(&poly, 1, state.max_pars, ctx);
    
    state.input = poly_str;
    state.pos = 0;
    state.len = strlen(poly_str);
    next_token(&state);
    
    parse_expression(&state, &poly);
    
    // Convert to univariate polynomial
    fq_nmod_poly_t univ_poly;
    fq_nmod_poly_init(univ_poly, ctx);
    
    // Convert fq_mvpoly to fq_nmod_poly (assuming only one variable)
    for (slong i = 0; i < poly.nterms; i++) {
        slong degree = 0;
        if (poly.terms[i].var_exp && poly.terms[i].var_exp[0] > 0) {
            degree = poly.terms[i].var_exp[0];
        }
        
        // Check if there are parameter terms (final univariate polynomial shouldn't have any)
        if (poly.terms[i].par_exp) {
            int has_params = 0;
            for (slong j = 0; j < state.npars; j++) {
                if (poly.terms[i].par_exp[j] > 0) {
                    has_params = 1;
                    break;
                }
            }
            if (has_params) {
                printf("Warning: univariate equation still contains parameters or generator\n");
            }
        }
        
        fq_nmod_poly_set_coeff(univ_poly, degree, poly.terms[i].coeff, ctx);
    }
    
    // Find roots
    fq_nmod_roots_t roots;
    fq_nmod_roots_init(roots, ctx);
    
    slong num_roots = our_fq_nmod_poly_roots(roots, univ_poly, 1, ctx);
    
    fq_nmod_t *result_roots = NULL;
    if (num_roots > 0) {
        result_roots = (fq_nmod_t*) malloc(roots->num * sizeof(fq_nmod_t));
        for (slong i = 0; i < roots->num; i++) {
            fq_nmod_init(result_roots[i], ctx);
            fq_nmod_set(result_roots[i], &roots->roots[i], ctx);
        }
        *num_roots_out = roots->num;
    } else {
        *num_roots_out = 0;
    }
    
    // Cleanup
    fq_nmod_roots_clear(roots, ctx);
    fq_nmod_poly_clear(univ_poly, ctx);
    fq_mvpoly_clear(&poly);
    
    for (slong i = 0; i < state.npars; i++) {
        free(state.par_names[i]);
    }
    free(state.par_names);
    free(state.var_names[0]);
    free(state.var_names);
    if (state.generator_name) free(state.generator_name);
    fq_nmod_clear(state.current.value, ctx);
    if (state.current.str) free(state.current.str);
    
    return result_roots;
}

// Substitute a specific variable with its value in a polynomial string
char* substitute_variable_in_polynomial(const char *poly_str, const char *var_name,
                                                const fq_nmod_t value, const fq_nmod_ctx_t ctx) {
    // Use parsing mechanism from DIXON_INTERFACE_FLINT_H
    char *gen_name = get_generator_name_for_solver(ctx);
    
    parser_state_t state = {0};
    state.var_names = (char**) malloc(sizeof(char*));
    state.var_names[0] = strdup(var_name);
    state.nvars = 1;
    state.npars = 0;
    state.max_pars = 16;
    state.par_names = (char**) malloc(state.max_pars * sizeof(char*));
    state.ctx = ctx;
    state.current.str = NULL;
    fq_nmod_init(state.current.value, ctx);
    state.generator_name = gen_name;
    
    // Parse polynomial
    fq_mvpoly_t poly;
    fq_mvpoly_init(&poly, 1, state.max_pars, ctx);
    
    state.input = poly_str;
    state.pos = 0;
    state.len = strlen(poly_str);
    next_token(&state);
    
    parse_expression(&state, &poly);
    
    // Create result polynomial by substituting variable value
    fq_mvpoly_t result_poly;
    fq_mvpoly_init(&result_poly, 0, state.npars, ctx);
    
    // Substitute variable for each term
    for (slong i = 0; i < poly.nterms; i++) {
        fq_nmod_t term_coeff;
        fq_nmod_init(term_coeff, ctx);
        fq_nmod_set(term_coeff, poly.terms[i].coeff, ctx);
        
        // If this term contains the variable to substitute
        if (poly.terms[i].var_exp && poly.terms[i].var_exp[0] > 0) {
            slong power = poly.terms[i].var_exp[0];
            fq_nmod_t var_power;
            fq_nmod_init(var_power, ctx);
            fq_nmod_pow_ui(var_power, value, power, ctx);
            fq_nmod_mul(term_coeff, term_coeff, var_power, ctx);
            fq_nmod_clear(var_power, ctx);
        }
        
        // Add processed term to result (preserve parameter exponents)
        fq_mvpoly_add_term(&result_poly, NULL, poly.terms[i].par_exp, term_coeff);
        fq_nmod_clear(term_coeff, ctx);
    }
    
    // Convert result to string
    char *result_str = fq_mvpoly_to_string(&result_poly, state.par_names, gen_name);
    
    // Cleanup
    fq_mvpoly_clear(&poly);
    fq_mvpoly_clear(&result_poly);
    
    for (slong i = 0; i < state.npars; i++) {
        free(state.par_names[i]);
    }
    free(state.par_names);
    free(state.var_names[0]);
    free(state.var_names);
    if (state.generator_name) free(state.generator_name);
    fq_nmod_clear(state.current.value, ctx);
    if (state.current.str) free(state.current.str);
    
    return result_str;
}

// ENHANCED: Use Dixon resultant to eliminate variable with equation selection
char* eliminate_variable_dixon_with_selection(char **poly_strings, slong num_polys, 
                                              const char *elim_var, const fq_nmod_ctx_t ctx,
                                              equation_combination_t *combination) {
    
    if (combination->num_equations == 2) {
        // Bivariate case, use bivariate resultant directly
        char *poly1 = poly_strings[combination->equation_indices[0]];
        char *poly2 = poly_strings[combination->equation_indices[1]];
        
        printf("    Using equations %ld and %ld for elimination\n", 
               combination->equation_indices[0], combination->equation_indices[1]);
        
        return bivariate_resultant(poly1, poly2, elim_var, ctx);
    } else if (combination->num_equations > 2) {
        // Multi-polynomial case, use Dixon resultant
        const char **poly_const = (const char**) malloc(combination->num_equations * sizeof(char*));
        const char **vars_array = (const char**) malloc(sizeof(char*));
        
        printf("    Using equations: ");
        for (slong i = 0; i < combination->num_equations; i++) {
            poly_const[i] = poly_strings[combination->equation_indices[i]];
            if (i > 0) printf(", ");
            printf("%ld", combination->equation_indices[i]);
        }
        printf(" for elimination\n");
        
        vars_array[0] = elim_var;
        
        char *result = dixon(poly_const, combination->num_equations, vars_array, 1, ctx);
        
        free(poly_const);
        free(vars_array);
        
        return result;
    } else {
        // Single polynomial, cannot eliminate
        return strdup("0");
    }
}

// ============= SOLUTION VERIFICATION =============

// Verify a single solution by substituting into original equations
int verify_solution_set(char **original_polys, slong num_polys,
                       variable_info_t *sorted_vars, slong num_vars,
                       fq_nmod_t **solution_values, const fq_nmod_ctx_t ctx) {
    
    printf("    Verifying solution: ");
    for (slong v = 0; v < num_vars; v++) {
        if (v > 0) printf(", ");
        printf("%s=", sorted_vars[v].name);
        fq_nmod_print_pretty(solution_values[v][0], ctx);
    }
    printf("\n");
    
    // Substitute all variables into each original equation
    for (slong poly_idx = 0; poly_idx < num_polys; poly_idx++) {
        char *current_poly = strdup(original_polys[poly_idx]);
        
        // Substitute each variable in turn
        for (slong var_idx = 0; var_idx < num_vars; var_idx++) {
            char *next_poly = substitute_variable_in_polynomial(current_poly,
                                                               sorted_vars[var_idx].name,
                                                               solution_values[var_idx][0],
                                                               ctx);
            free(current_poly);
            current_poly = next_poly;
        }
        
        // Check if the result is zero (or close to zero)
        // Parse the result to check if it's zero
        char *gen_name = get_generator_name(ctx);
        parser_state_t state = {0};
        state.var_names = NULL;
        state.nvars = 0;
        state.npars = 0;
        state.max_pars = 16;
        state.par_names = (char**) malloc(state.max_pars * sizeof(char*));
        state.ctx = ctx;
        state.current.str = NULL;
        fq_nmod_init(state.current.value, ctx);
        state.generator_name = strdup(gen_name);
        
        fq_mvpoly_t result_poly;
        fq_mvpoly_init(&result_poly, 0, state.max_pars, ctx);
        
        state.input = current_poly;
        state.pos = 0;
        state.len = strlen(current_poly);
        next_token(&state);
        
        parse_expression(&state, &result_poly);
        
        // Check if polynomial is zero
        int is_zero = (result_poly.nterms == 0);
        if (!is_zero && result_poly.nterms == 1) {
            // Check if the single term is zero coefficient
            is_zero = fq_nmod_is_zero(result_poly.terms[0].coeff, ctx);
        }
        
        if (!is_zero) {
            printf("    Verification FAILED for equation %ld: %s != 0\n", poly_idx + 1, current_poly);
            
            // Cleanup
            fq_mvpoly_clear(&result_poly);
            for (slong i = 0; i < state.npars; i++) {
                free(state.par_names[i]);
            }
            free(state.par_names);
            if (state.generator_name) free(state.generator_name);
            fq_nmod_clear(state.current.value, ctx);
            if (state.current.str) free(state.current.str);
            free(gen_name);
            free(current_poly);
            return 0;
        }
        
        // Cleanup
        fq_mvpoly_clear(&result_poly);
        for (slong i = 0; i < state.npars; i++) {
            free(state.par_names[i]);
        }
        free(state.par_names);
        if (state.generator_name) free(state.generator_name);
        fq_nmod_clear(state.current.value, ctx);
        if (state.current.str) free(state.current.str);
        free(gen_name);
        free(current_poly);
    }
    
    printf("    Verification PASSED\n");
    return 1;
}

// Filter solution sets by verification
void filter_solutions_by_verification(polynomial_solutions_t *sols,
                                     char **original_polys, slong num_polys,
                                     variable_info_t *sorted_vars) {
    if (!sols->is_valid || sols->has_no_solutions || sols->num_solution_sets == 0) {
        return;
    }
    
    printf("\n=== Verifying Solution Sets ===\n");
    
    slong verified_count = 0;
    int *valid_sets = (int*) calloc(sols->num_solution_sets, sizeof(int));
    
    // Check each solution set
    for (slong set = 0; set < sols->num_solution_sets; set++) {
        printf("  Checking solution set %ld...\n", set + 1);
        int pass = 0;
        char *candidate_line = format_solution_set_line_from_arrays(
            sols->variable_names,
            sols->num_variables,
            sols->solution_sets[set],
            &sols->solutions_per_var[set * sols->num_variables],
            sols->ctx
        );
        
        // Check if this solution set is complete
        int is_complete = 1;
        for (slong var = 0; var < sols->num_variables; var++) {
            if (sols->solutions_per_var[set * sols->num_variables + var] == 0) {
                printf("    Solution set %ld incomplete (missing value for %s)\n", 
                       set + 1, sols->variable_names[var]);
                is_complete = 0;
                break;
            }
        }
        
        if (!is_complete) {
            add_candidate_solution(sols, candidate_line, 0);
            continue;
        }
        
        // Verify this solution
        if (verify_solution_set(original_polys, num_polys, sorted_vars, sols->num_variables,
                               sols->solution_sets[set], sols->ctx)) {
            valid_sets[set] = 1;
            verified_count++;
            pass = 1;
        }
        add_candidate_solution(sols, candidate_line, pass);
    }
    
    printf("  Verification complete: %ld out of %ld solution sets are valid\n", 
           verified_count, sols->num_solution_sets);
    sols->checked_solution_sets = sols->num_solution_sets;
    sols->verified_solution_sets = verified_count;
    
    if (verified_count == 0) {
        // No valid solutions - mark as having no solutions
        sols->has_no_solutions = 1;
        
        // Clear existing solution sets
        if (sols->solution_sets) {
            for (slong set = 0; set < sols->num_solution_sets; set++) {
                if (sols->solution_sets[set]) {
                    for (slong var = 0; var < sols->num_variables; var++) {
                        if (sols->solution_sets[set][var]) {
                            for (slong sol = 0; sol < sols->solutions_per_var[set * sols->num_variables + var]; sol++) {
                                fq_nmod_clear(sols->solution_sets[set][var][sol], sols->ctx);
                            }
                            free(sols->solution_sets[set][var]);
                        }
                    }
                    free(sols->solution_sets[set]);
                }
            }
            free(sols->solution_sets);
            sols->solution_sets = NULL;
        }
        
        if (sols->solutions_per_var) {
            free(sols->solutions_per_var);
            sols->solutions_per_var = NULL;
        }
        
        sols->num_solution_sets = 0;
        
    } else if (verified_count < sols->num_solution_sets) {
        // Some solutions invalid - create new filtered structure
        fq_nmod_t ***new_solution_sets = (fq_nmod_t***) malloc(verified_count * sizeof(fq_nmod_t**));
        slong *new_solutions_per_var = (slong*) calloc(verified_count * sols->num_variables, sizeof(slong));
        
        slong new_idx = 0;
        for (slong set = 0; set < sols->num_solution_sets; set++) {
            if (valid_sets[set]) {
                // Copy this valid solution set
                new_solution_sets[new_idx] = sols->solution_sets[set];
                for (slong var = 0; var < sols->num_variables; var++) {
                    new_solutions_per_var[new_idx * sols->num_variables + var] = 
                        sols->solutions_per_var[set * sols->num_variables + var];
                }
                sols->solution_sets[set] = NULL; // Prevent double-free
                new_idx++;
            } else {
                // Free invalid solution set
                if (sols->solution_sets[set]) {
                    for (slong var = 0; var < sols->num_variables; var++) {
                        if (sols->solution_sets[set][var]) {
                            for (slong sol = 0; sol < sols->solutions_per_var[set * sols->num_variables + var]; sol++) {
                                fq_nmod_clear(sols->solution_sets[set][var][sol], sols->ctx);
                            }
                            free(sols->solution_sets[set][var]);
                        }
                    }
                    free(sols->solution_sets[set]);
                }
            }
        }
        
        // Replace old structure with new one
        free(sols->solution_sets);
        free(sols->solutions_per_var);
        sols->solution_sets = new_solution_sets;
        sols->solutions_per_var = new_solutions_per_var;
        sols->num_solution_sets = verified_count;
    }
    
    free(valid_sets);
    printf("=== Verification Complete ===\n\n");
}

// ============= RECURSIVE BACK SUBSTITUTION =============

// ENHANCED: Back substitution with robust equation selection
int solve_by_back_substitution_recursive_enhanced(char **original_polys, slong num_polys,
                                                  variable_info_t *sorted_vars, slong num_vars,
                                                  fq_nmod_t *base_solutions, slong num_base_solutions,
                                                  polynomial_solutions_t *sols) {
    
    // Skip if already dimension > 0
    if (sols->has_no_solutions == -1) {
        printf("Back substitution SKIPPED: system already has dimension > 0 status\n");
        return 1;
    }
    
    if (num_vars <= 1) {
        return 1; 
    }
    
    printf("Starting back substitution process...\n");
    printf("Base variable %s has %ld solutions\n", sorted_vars[0].name, num_base_solutions);
    
    // Collect solution sets
    fq_nmod_t ***all_solution_sets = NULL;
    slong *all_solutions_per_var = NULL;
    slong total_solution_sets = 0;
    
    for (slong base_idx = 0; base_idx < num_base_solutions; base_idx++) {
        printf("\n--- Processing base solution %ld: ", base_idx + 1);
        fq_nmod_print_pretty(base_solutions[base_idx], sols->ctx);
        printf(" ---\n");
        // Substitute base solution
        char **reduced_polys = (char**) malloc(num_polys * sizeof(char*));
        slong num_nonzero_polys = 0;
        
        for (slong poly_idx = 0; poly_idx < num_polys; poly_idx++) {
            char *substituted = substitute_variable_in_polynomial(original_polys[poly_idx],
                                                                  sorted_vars[0].name,
                                                                  base_solutions[base_idx],
                                                                  sols->ctx);
            printf("Equation %ld after substituting %s: %s\n", 
                   poly_idx, sorted_vars[0].name, substituted);
            
            if (strcmp(substituted, "0") != 0) {
                reduced_polys[num_nonzero_polys] = substituted;
                num_nonzero_polys++;
            } else {
                printf("Equation %ld became zero (satisfied)\n", poly_idx);
                free(substituted);
            }
        }
        
        if (num_nonzero_polys == 0) {
            printf("All equations satisfied by base solution\n");
            // All equations satisfied - valid solution
            if (total_solution_sets == 0) {
                all_solution_sets = (fq_nmod_t***) malloc(sizeof(fq_nmod_t**));
                all_solutions_per_var = (slong*) calloc(num_vars, sizeof(slong));
            } else {
                all_solution_sets = (fq_nmod_t***) realloc(all_solution_sets, 
                                                          (total_solution_sets + 1) * sizeof(fq_nmod_t**));
                all_solutions_per_var = (slong*) realloc(all_solutions_per_var,
                                                         (total_solution_sets + 1) * num_vars * sizeof(slong));
            }
            
            all_solution_sets[total_solution_sets] = (fq_nmod_t**) calloc(num_vars, sizeof(fq_nmod_t*));
            all_solution_sets[total_solution_sets][0] = (fq_nmod_t*) malloc(sizeof(fq_nmod_t));
            fq_nmod_init(all_solution_sets[total_solution_sets][0][0], sols->ctx);
            fq_nmod_set(all_solution_sets[total_solution_sets][0][0], base_solutions[base_idx], sols->ctx);
            all_solutions_per_var[total_solution_sets * num_vars + 0] = 1;
            
            for (slong v = 1; v < num_vars; v++) {
                all_solutions_per_var[total_solution_sets * num_vars + v] = 0;
            }
            
            total_solution_sets++;
            free(reduced_polys);
            continue;
        }
        
        // Try to solve remaining system
        slong target_equations = num_vars - 1;
        
        if (num_nonzero_polys >= target_equations) {
            equation_combination_t *combinations = NULL;
            slong num_combinations = 0;
            
            if (num_nonzero_polys == target_equations) {
                num_combinations = 1;
                combinations = (equation_combination_t*) malloc(sizeof(equation_combination_t));
                combinations[0].equation_indices = (slong*) malloc(target_equations * sizeof(slong));
                combinations[0].num_equations = target_equations;
                for (slong i = 0; i < target_equations; i++) {
                    combinations[0].equation_indices[i] = i;
                }
            } else {
                generate_equation_combinations(num_nonzero_polys, target_equations, 
                                             &combinations, &num_combinations);
            }
            
            printf("Generated %ld equation combinations to try\n", num_combinations);
            
            // CRITICAL: Track the types of failures
            int found_finite_solution = 0;
            int found_no_solution = 0;
            int found_dimension_gt_zero = 0;
            
            for (slong comb_idx = 0; comb_idx < num_combinations; comb_idx++) {
                printf("  Trying combination %ld: equations ", comb_idx + 1);
                for (slong i = 0; i < target_equations; i++) {
                    if (i > 0) printf(", ");
                    printf("%ld", combinations[comb_idx].equation_indices[i]);
                }
                printf("\n");
                
                char **final_polys = (char**) malloc(target_equations * sizeof(char*));
                for (slong i = 0; i < target_equations; i++) {
                    final_polys[i] = strdup(reduced_polys[combinations[comb_idx].equation_indices[i]]);
                }

                {
                    char *base_str = fq_nmod_get_str_pretty(base_solutions[base_idx], sols->ctx);
                    char *eq_list = format_equation_index_list(combinations[comb_idx].equation_indices,
                                                               target_equations);
                    if (target_equations == 1) {
                        add_resultant_step(sols,
                                           "Base %s = %s -> solve reduced eq(%s)",
                                           sorted_vars[0].name,
                                           base_str ? base_str : "?",
                                           eq_list ? eq_list : "?");
                    } else {
                        add_resultant_step(sols,
                                           "Base %s = %s -> compute reduced resultants from eq(%s)",
                                           sorted_vars[0].name,
                                           base_str ? base_str : "?",
                                           eq_list ? eq_list : "?");
                    }
                    if (eq_list) free(eq_list);
                    if (base_str) free(base_str);
                }
                
                variable_info_t *remaining_vars = (variable_info_t*) malloc((num_vars - 1) * sizeof(variable_info_t));
                for (slong i = 1; i < num_vars; i++) {
                    remaining_vars[i - 1] = sorted_vars[i];
                }
                
                // Recursively solve
                polynomial_solutions_t *test_sols = solve_polynomial_system_array_with_vars(final_polys, 
                                                                                           target_equations,
                                                                                           remaining_vars,
                                                                                           num_vars - 1,
                                                                                           sols->ctx);
                
                if (test_sols && test_sols->is_valid) {
                    merge_solver_logs(sols, test_sols);
                    if (test_sols->has_no_solutions == -1) {
                        printf("  Combination %ld: dimension > 0\n", comb_idx + 1);
                        found_dimension_gt_zero++;
                    } else if (test_sols->has_no_solutions == 1) {
                        printf("  Combination %ld: no solutions\n", comb_idx + 1);
                        found_no_solution++;
                    } else {
                        printf("  ✓ Combination %ld succeeded with finite solutions!\n", comb_idx + 1);
                        found_finite_solution++;
                        
                        // Process successful solutions (existing code)
                        if (test_sols->num_solution_sets > 0) {
                            slong old_total = total_solution_sets;
                            total_solution_sets += test_sols->num_solution_sets;
                            
                            if (old_total == 0) {
                                all_solution_sets = (fq_nmod_t***) malloc(total_solution_sets * sizeof(fq_nmod_t**));
                                all_solutions_per_var = (slong*) calloc(total_solution_sets * num_vars, sizeof(slong));
                            } else {
                                all_solution_sets = (fq_nmod_t***) realloc(all_solution_sets, 
                                                                          total_solution_sets * sizeof(fq_nmod_t**));
                                all_solutions_per_var = (slong*) realloc(all_solutions_per_var,
                                                                         total_solution_sets * num_vars * sizeof(slong));
                            }
                            
                            // Combine solutions
                            for (slong rec_set = 0; rec_set < test_sols->num_solution_sets; rec_set++) {
                                slong combined_idx = old_total + rec_set;
                                
                                all_solution_sets[combined_idx] = (fq_nmod_t**) calloc(num_vars, sizeof(fq_nmod_t*));
                                
                                all_solution_sets[combined_idx][0] = (fq_nmod_t*) malloc(sizeof(fq_nmod_t));
                                fq_nmod_init(all_solution_sets[combined_idx][0][0], sols->ctx);
                                fq_nmod_set(all_solution_sets[combined_idx][0][0], base_solutions[base_idx], sols->ctx);
                                all_solutions_per_var[combined_idx * num_vars + 0] = 1;
                                
                                for (slong var = 1; var < num_vars; var++) {
                                    slong rec_var_idx = var - 1;
                                    slong num_rec_sols = test_sols->solutions_per_var[rec_set * test_sols->num_variables + rec_var_idx];
                                    
                                    if (num_rec_sols > 0) {
                                        all_solution_sets[combined_idx][var] = (fq_nmod_t*) malloc(sizeof(fq_nmod_t));
                                        fq_nmod_init(all_solution_sets[combined_idx][var][0], sols->ctx);
                                        fq_nmod_set(all_solution_sets[combined_idx][var][0],
                                                   test_sols->solution_sets[rec_set][rec_var_idx][0], sols->ctx);
                                        all_solutions_per_var[combined_idx * num_vars + var] = 1;
                                    } else {
                                        all_solutions_per_var[combined_idx * num_vars + var] = 0;
                                    }
                                }
                            }
                        }
                        
                        // Found working combination, can break
                        if (test_sols) {
                            polynomial_solutions_clear(test_sols);
                            free(test_sols);
                        }
                        for (slong i = 0; i < target_equations; i++) {
                            free(final_polys[i]);
                        }
                        free(final_polys);
                        free(remaining_vars);
                        break; // Found working solution, stop trying more combinations
                    }
                } else {
                    printf("  ✗ Combination %ld failed to solve\n", comb_idx + 1);
                }
                
                if (test_sols) {
                    polynomial_solutions_clear(test_sols);
                    free(test_sols);
                }
                for (slong i = 0; i < target_equations; i++) {
                    free(final_polys[i]);
                }
                free(final_polys);
                free(remaining_vars);
            }
            
            // CRITICAL: Analyze the types of failures
            if (found_finite_solution == 0) {
                // No finite solutions found
                if (found_dimension_gt_zero > 0) {
                    printf("All equation combinations for base solution %ld resulted in dimension > 0\n", base_idx + 1);
                    printf("PROPAGATING DIMENSION > 0 to main system\n");
                    
                    // CRITICAL: If all recursive calls have dimension > 0, the main system has dimension > 0
                    sols->has_no_solutions = -1;
                    
                    // Cleanup and return immediately
                    free_equation_combinations(combinations, num_combinations);
                    for (slong i = 0; i < num_nonzero_polys; i++) {
                        free(reduced_polys[i]);
                    }
                    free(reduced_polys);
                    return 1;  // Early return with dimension > 0 status
                } else {
                    printf("All equation combinations failed (no solutions) for base solution %ld\n", base_idx + 1);
                }
            }
            
            free_equation_combinations(combinations, num_combinations);
        } else {
            printf("Not enough equations (%ld) for remaining variables (%ld)\n", 
                   num_nonzero_polys, num_vars - 1);
        }
        
        // Cleanup
        for (slong i = 0; i < num_nonzero_polys; i++) {
            free(reduced_polys[i]);
        }
        free(reduced_polys);
    }
    
    // Update solution structure - but only if we haven't detected dimension > 0
    if (sols->has_no_solutions != -1) {
        if (total_solution_sets > 0) {
            sols->num_solution_sets = total_solution_sets;
            sols->solution_sets = all_solution_sets;
            sols->solutions_per_var = all_solutions_per_var;
            printf("Back substitution completed: %ld total solution sets\n", total_solution_sets);
        } else {
            printf("Back substitution completed: system has no solutions\n");
            sols->has_no_solutions = 1;
        }
    } else {
        printf("Back substitution completed: dimension > 0 status maintained\n");
    }
    
    return 1;
}

// ============= MAIN SOLVING ALGORITHMS =============

// ENHANCED: Step-by-step elimination solving with robust equation selection
int solve_by_elimination_enhanced(char **poly_strings, slong num_polys,
                                 variable_info_t *sorted_vars, slong num_vars,
                                 polynomial_solutions_t *sols) {
    
    printf("=== Starting enhanced step-by-step elimination solving ===\n");
    
    if (num_vars == 1) {
        // Single variable case (unchanged)
        char *final_var = sorted_vars[0].name;
        slong num_base_solutions;
        sols->total_combinations = 1;
        fq_nmod_t *base_solutions = solve_univariate_equation_all_roots(poly_strings[0], final_var, 
                                                                       &num_base_solutions, sols->ctx);
        sols->num_base_solutions = num_base_solutions;
        if (!sols->elimination_summary) {
            size_t len = strlen(final_var) + 32;
            sols->elimination_summary = (char *) malloc(len);
            if (sols->elimination_summary) {
                snprintf(sols->elimination_summary, len, "Directly solve for %s", final_var);
            }
        }
        
        if (num_base_solutions == 0) {
            sols->is_valid = 1;
            sols->has_no_solutions = 1;
            return 1;
        }
        
        if (!base_solutions) {
            return 0;
        }
        
        sols->num_solution_sets = num_base_solutions;
        sols->solution_sets = (fq_nmod_t***) malloc(num_base_solutions * sizeof(fq_nmod_t**));
        sols->solutions_per_var = (slong*) calloc(num_base_solutions * num_vars, sizeof(slong));
        
        for (slong i = 0; i < num_base_solutions; i++) {
            sols->solution_sets[i] = (fq_nmod_t**) calloc(num_vars, sizeof(fq_nmod_t*));
            sols->solution_sets[i][0] = (fq_nmod_t*) malloc(sizeof(fq_nmod_t));
            fq_nmod_init(sols->solution_sets[i][0][0], sols->ctx);
            fq_nmod_set(sols->solution_sets[i][0][0], base_solutions[i], sols->ctx);
            sols->solutions_per_var[i * num_vars + 0] = 1;
        }
        sols->successful_combinations = (num_base_solutions > 0) ? 1 : 0;
        
        for (slong i = 0; i < num_base_solutions; i++) {
            fq_nmod_clear(base_solutions[i], sols->ctx);
        }
        free(base_solutions);
        return 1;
    }
    
    // Multi-variable elimination (existing logic unchanged)
    printf("%ld*%ld system detected\n", num_vars, num_polys);
    
    char *elim_var = sorted_vars[num_vars-1].name;
    char *keep_var = sorted_vars[0].name;

    if (sols->elimination_summary) {
        free(sols->elimination_summary);
        sols->elimination_summary = NULL;
    }
    {
        size_t len = strlen(elim_var) + strlen(keep_var) + 32;
        sols->elimination_summary = (char *) malloc(len);
        if (sols->elimination_summary) {
            snprintf(sols->elimination_summary, len, "Eliminate %s, keep %s", elim_var, keep_var);
        }
    }
    
    printf("--- Attempting elimination: eliminate %s, keep %s ---\n", elim_var, keep_var);
    
    // Try ALL elimination combinations and count successes
    slong successful_combinations = 0;
    slong total_combinations = 0;
    char *working_resultant = NULL;
    
    // [Existing elimination logic - try all combinations, count successes]
    // This part was correct in the previous version
    
    if (num_vars == 2) {
        if (num_polys == 2) {
            total_combinations = 1;
            add_resultant_step(sols,
                               "Compute resultant of eq(1,2) eliminating %s",
                               elim_var);
            char *resultant = bivariate_resultant(poly_strings[0], poly_strings[1], elim_var, sols->ctx);
            if (resultant && strcmp(resultant, "0") != 0) {
                working_resultant = resultant;
                successful_combinations = 1;
                add_resultant_step(sols,
                                   "eq(1,2) -> non-zero resultant in %s",
                                   keep_var);
            } else {
                printf("Bivariate resultant is zero\n");
                add_resultant_step(sols,
                                   "eq(1,2) -> zero resultant");
                if (resultant) free(resultant);
            }
        } else if (num_polys > 2) {
            equation_combination_t *combinations = NULL;
            slong num_combinations = 0;
            generate_equation_combinations(num_polys, 2, &combinations, &num_combinations);
            
            total_combinations = num_combinations;
            
            for (slong comb_idx = 0; comb_idx < num_combinations; comb_idx++) {
                printf("  Trying combination %ld: equations %ld and %ld\n", 
                       comb_idx + 1, 
                       combinations[comb_idx].equation_indices[0],
                       combinations[comb_idx].equation_indices[1]);
                add_resultant_step(sols,
                                   "Compute resultant of eq(%ld,%ld) eliminating %s",
                                   combinations[comb_idx].equation_indices[0] + 1,
                                   combinations[comb_idx].equation_indices[1] + 1,
                                   elim_var);
                
                char *test_resultant = eliminate_variable_dixon_with_selection(poly_strings, num_polys, 
                                                                              elim_var, sols->ctx,
                                                                              &combinations[comb_idx]);
                
                if (test_resultant && strcmp(test_resultant, "0") != 0) {
                    printf("  ✓ Combination %ld succeeded!\n", comb_idx + 1);
                    add_resultant_step(sols,
                                       "eq(%ld,%ld) -> non-zero resultant in %s",
                                       combinations[comb_idx].equation_indices[0] + 1,
                                       combinations[comb_idx].equation_indices[1] + 1,
                                       keep_var);
                    if (!working_resultant) {
                        working_resultant = test_resultant;
                    } else {
                        free(test_resultant);
                    }
                    successful_combinations++;
                } else {
                    printf("  ✗ Combination %ld failed (resultant = 0)\n", comb_idx + 1);
                    add_resultant_step(sols,
                                       "eq(%ld,%ld) -> zero resultant",
                                       combinations[comb_idx].equation_indices[0] + 1,
                                       combinations[comb_idx].equation_indices[1] + 1);
                    if (test_resultant) free(test_resultant);
                }
            }
            
            free_equation_combinations(combinations, num_combinations);
        }
    } else {
        // Multivariate case (similar logic)
        slong num_elim_vars = num_vars - 1;
        const char **elim_vars = (const char**) malloc(num_elim_vars * sizeof(char*));
        for (slong i = 0; i < num_elim_vars; i++) {
            elim_vars[i] = sorted_vars[i + 1].name;
        }
        
        if (num_polys == num_vars) {
            total_combinations = 1;
            const char **poly_const = (const char**) malloc(num_polys * sizeof(char*));
            for (slong i = 0; i < num_polys; i++) {
                poly_const[i] = poly_strings[i];
            }
            add_resultant_step(sols,
                               "Compute Dixon resultant of eq(1..%ld) eliminating %ld variable(s)",
                               num_polys, num_elim_vars);
            
            char *resultant = dixon(poly_const, num_polys, elim_vars, num_elim_vars, sols->ctx);
            
            if (resultant && strcmp(resultant, "0") != 0) {
                working_resultant = resultant;
                successful_combinations = 1;
                add_resultant_step(sols,
                                   "eq(1..%ld) -> non-zero resultant in %s",
                                   num_polys, keep_var);
            } else {
                add_resultant_step(sols,
                                   "eq(1..%ld) -> zero resultant",
                                   num_polys);
                if (resultant) free(resultant);
            }
            
            free(poly_const);
        } else if (num_polys > num_vars) {
            equation_combination_t *combinations = NULL;
            slong num_combinations = 0;
            generate_equation_combinations(num_polys, num_vars, &combinations, &num_combinations);
            
            total_combinations = num_combinations;
            
            for (slong comb_idx = 0; comb_idx < num_combinations; comb_idx++) {
                printf("  Trying combination %ld: equations ", comb_idx + 1);
                for (slong i = 0; i < num_vars; i++) {
                    if (i > 0) printf(", ");
                    printf("%ld", combinations[comb_idx].equation_indices[i]);
                }
                printf("\n");
                char *eq_list = format_equation_index_list(combinations[comb_idx].equation_indices,
                                                           num_vars);
                add_resultant_step(sols,
                                   "Compute Dixon resultant of eq(%s) eliminating %ld variable(s)",
                                   eq_list ? eq_list : "?",
                                   num_elim_vars);
                
                const char **poly_const = (const char**) malloc(num_vars * sizeof(char*));
                for (slong i = 0; i < num_vars; i++) {
                    poly_const[i] = poly_strings[combinations[comb_idx].equation_indices[i]];
                }
                
                char *test_resultant = dixon(poly_const, num_vars, elim_vars, num_elim_vars, sols->ctx);
                
                free(poly_const);
                
                if (test_resultant && strcmp(test_resultant, "0") != 0) {
                    printf("  ✓ Combination %ld succeeded!\n", comb_idx + 1);
                    add_resultant_step(sols,
                                       "eq(%s) -> non-zero resultant in %s",
                                       eq_list ? eq_list : "?",
                                       keep_var);
                    if (!working_resultant) {
                        working_resultant = test_resultant;
                    } else {
                        free(test_resultant);
                    }
                    successful_combinations++;
                } else {
                    printf("  ✗ Combination %ld failed (resultant = 0)\n", comb_idx + 1);
                    add_resultant_step(sols,
                                       "eq(%s) -> zero resultant",
                                       eq_list ? eq_list : "?");
                    if (test_resultant) free(test_resultant);
                }
                if (eq_list) free(eq_list);
            }
            
            free_equation_combinations(combinations, num_combinations);
        }
        
        free(elim_vars);
    }
    
    // CRITICAL DECISION: Only report dimension > 0 if ALL combinations failed
    if (successful_combinations == 0) {
        sols->total_combinations = total_combinations;
        sols->successful_combinations = 0;
        printf("All %ld equation combinations resulted in zero resultant\n", total_combinations);
        printf("System has dimension greater than zero\n");
        sols->is_valid = 1;
        sols->has_no_solutions = -1;
        printf("EARLY RETURN: Dimension > 0 detected, skipping all further processing\n");
        return 1;
    }

    sols->total_combinations = total_combinations;
    sols->successful_combinations = successful_combinations;
    
    printf("Found %ld non-zero resultant(s) out of %ld combinations, proceeding with solving...\n", 
           successful_combinations, total_combinations);
    
    // Solve univariate resultant
    slong num_base_solutions;
    fq_nmod_t *base_solutions = solve_univariate_equation_all_roots(working_resultant, keep_var, 
                                                                   &num_base_solutions, sols->ctx);
    sols->num_base_solutions = num_base_solutions;
    
    if (num_base_solutions == 0) {
        printf("Univariate resultant has no solutions\n");
        sols->is_valid = 1;
        sols->has_no_solutions = 1;
        free(working_resultant);
        return 1;
    }
    
    if (!base_solutions) {
        printf("Failed to solve univariate resultant\n");
        free(working_resultant);
        return 0;
    }
    
    printf("Found %ld value(s) for %s, starting back substitution...\n", num_base_solutions, keep_var);
    free(working_resultant);
    
    // Do back substitution
    int backtrack_success = solve_by_back_substitution_recursive_enhanced(poly_strings, num_polys, 
                                                                         sorted_vars, num_vars, 
                                                                         base_solutions, num_base_solutions, sols);
    
    // Cleanup
    for (slong i = 0; i < num_base_solutions; i++) {
        fq_nmod_clear(base_solutions[i], sols->ctx);
    }
    free(base_solutions);
    
    return backtrack_success;
}

// ============= MAIN SOLVER INTERFACES =============

// ENHANCED: Main solver function with robust equation selection
polynomial_solutions_t* solve_polynomial_system_array_with_vars(char **poly_strings, slong num_polys,
                                                               variable_info_t *original_vars, slong num_original_vars,
                                                               const fq_nmod_ctx_t ctx) {
    printf("\n=== Enhanced Polynomial System Solver (With Variable Info) ===\n");
    printf("Input equation count: %ld\n", num_polys);
    printf("Original variable count: %ld\n", num_original_vars);
    
    // Initialize solution structure
    polynomial_solutions_t *sols = (polynomial_solutions_t*) malloc(sizeof(polynomial_solutions_t));
    polynomial_solutions_init(sols, num_original_vars, ctx);
    sols->num_equations = num_polys;
    sols->variable_order = join_variable_names(original_vars, num_original_vars);
    
    // Copy variable names
    for (slong i = 0; i < num_original_vars; i++) {
        sols->variable_names[i] = strdup(original_vars[i].name);
    }
    
    // Execute enhanced step-by-step elimination solving
    int success = solve_by_elimination_enhanced(poly_strings, num_polys, original_vars, num_original_vars, sols);
    
    if (success) {
        sols->is_valid = 1;
        
        if (sols->error_message) {
            printf("Enhanced solving failed: %s\n", sols->error_message);
            sols->is_valid = 0; 
        } else if (sols->has_no_solutions == -1) {
            printf("Enhanced solving completed: polynomial system dimension greater than zero\n");
        } else if (sols->has_no_solutions == 1) {
            printf("Enhanced solving successful: system has no solutions over the finite field\n");
        } else {
            printf("Enhanced solving successful\n");
        }
    } else {
        printf("Enhanced solving failed\n");
        if (!sols->error_message) {
            sols->error_message = strdup("Enhanced elimination process failed");
        }
    }
    
    return sols;
}

// ENHANCED: Main solver function - array interface with robust equation selection and verification
polynomial_solutions_t* solve_polynomial_system_array(char **poly_strings, slong num_polys,
                                                      const fq_nmod_ctx_t ctx) {
    printf("\n=== Enhanced Polynomial System Solver ===\n");
    
    // Extract variables (existing code)
    variable_info_t *sorted_vars = NULL;
    slong num_vars = 0;
    
    if (!extract_and_sort_variables(poly_strings, num_polys, &sorted_vars, &num_vars, ctx)) {
        polynomial_solutions_t *sols = (polynomial_solutions_t*) malloc(sizeof(polynomial_solutions_t));
        polynomial_solutions_init(sols, 0, ctx);
        sols->error_message = strdup("Variable extraction failed");
        return sols;
    }
    
    polynomial_solutions_t *sols = (polynomial_solutions_t*) malloc(sizeof(polynomial_solutions_t));
    polynomial_solutions_init(sols, num_vars, ctx);
    sols->num_equations = num_polys;
    sols->variable_order = join_variable_names(sorted_vars, num_vars);
    
    for (slong i = 0; i < num_vars; i++) {
        sols->variable_names[i] = strdup(sorted_vars[i].name);
    }
    
    // Solve the system
    int success = solve_by_elimination_enhanced(poly_strings, num_polys, sorted_vars, num_vars, sols);
    
    if (success) {
        sols->is_valid = 1;
        
        if (sols->has_no_solutions == -1) {
            printf("FINAL STATUS: polynomial system dimension greater than zero\n");
        } else if (sols->has_no_solutions == 1) {
            printf("FINAL STATUS: system has no solutions over the finite field\n");
        } else {
            printf("FINAL STATUS: found finite solutions\n");
            if (sols->has_no_solutions == 0) {
                filter_solutions_by_verification(sols, poly_strings, num_polys, sorted_vars);
            }
        }
    } else {
        printf("Solving failed\n");
        sols->is_valid = 0;
        if (!sols->error_message) {
            sols->error_message = strdup("Elimination process failed");
        }
    }
    // Cleanup
    for (slong i = 0; i < num_vars; i++) {
        free(sorted_vars[i].name);
    }
    free(sorted_vars);
    return sols;
}

// String interface
polynomial_solutions_t* solve_polynomial_system_string(const char *poly_string, 
                                                       const fq_nmod_ctx_t ctx) {
    if (!poly_string || strlen(poly_string) == 0) {
        polynomial_solutions_t *sols = (polynomial_solutions_t*) malloc(sizeof(polynomial_solutions_t));
        polynomial_solutions_init(sols, 0, ctx);
        sols->error_message = strdup("Input string is empty");
        return sols;
    }
    
    // Split string
    slong num_polys;
    char **poly_array = split_string(poly_string, &num_polys);
    
    if (!poly_array || num_polys == 0) {
        polynomial_solutions_t *sols = (polynomial_solutions_t*) malloc(sizeof(polynomial_solutions_t));
        polynomial_solutions_init(sols, 0, ctx);
        sols->error_message = strdup("Polynomial string parsing failed");
        return sols;
    }
    
    // Call array interface
    polynomial_solutions_t *result = solve_polynomial_system_array(poly_array, num_polys, ctx);
    // Cleanup
    free_split_strings(poly_array, num_polys);
    return result;
}

// ============= TEST FUNCTION =============

// Test function
void test_polynomial_solver(void) {
    printf("\n=== Testing Enhanced Polynomial System Solver ===\n");
    
    fq_nmod_ctx_t ctx;
    mp_limb_t prime = 7;
    fmpz_t p;
    fmpz_init_set_ui(p, prime);
    fq_nmod_ctx_init(ctx, p, 1, "t");  
    
    // Test 1: Simple linear system
    printf("\n--- Test 1: Simple linear system ---\n");
    char *linear_polys[] = {
        "x + y - 3",
        "2*x - y - 0"
    };
    
    polynomial_solutions_t *sols1 = solve_polynomial_system_array(linear_polys, 2, ctx);
    print_polynomial_solutions(sols1);
    polynomial_solutions_clear(sols1);
    free(sols1);
    
    // Test 2: Example system from user (should demonstrate robust equation selection)
    printf("\n--- Test 2: User example system (x^2*y^9*z+x*y, x*y+z, x^4+z+1) ---\n");
    const char *user_system = "x^2*y^9*z+x*y, x*y+z, x^4+z+1";
    
    polynomial_solutions_t *sols2 = solve_polynomial_system_string(user_system, ctx);
    print_polynomial_solutions(sols2);
    polynomial_solutions_clear(sols2);
    free(sols2);
    
    // Test 3: System that should trigger "dimension > 0" error
    printf("\n--- Test 3: System with dimension > 0 ---\n");
    const char *high_dim_system = "x*y, x*y, x^2 + 1";
    
    polynomial_solutions_t *sols3 = solve_polynomial_system_string(high_dim_system, ctx);
    print_polynomial_solutions(sols3);
    polynomial_solutions_clear(sols3);
    free(sols3);
    
    printf("=== Enhanced Testing Complete ===\n");
}
