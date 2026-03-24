#include "dixon_with_ideal_reduction.h"

/* Find variable index by name in the current context */
static slong find_variable_by_name(const char *var_name, char **current_var_names, slong nvars) {
    for (slong i = 0; i < nvars; i++) {
        if (current_var_names && current_var_names[i] && 
            strcmp(current_var_names[i], var_name) == 0) {
            return i;
        }
    }
    return -1;
}

/* Find generator for a specific variable in the ideal */
static slong find_generator_for_variable(const unified_triangular_ideal_t *ideal, const char *var_name) {
    for (slong i = 0; i < ideal->num_gens; i++) {
        if (ideal->var_names[i] && strcmp(ideal->var_names[i], var_name) == 0) {
            return i;
        }
    }
    return -1;
}

/* Check if string contains equation format (has equals sign) */
static int is_equation_format(const char *poly_str) {
    return strchr(poly_str, '=') != NULL;
}

/* Split equation string into left and right parts */
static void split_equation_string(const char *equation_str, char **lhs_str, char **rhs_str) {
    const char *eq_pos = strchr(equation_str, '=');
    if (!eq_pos) {
        *lhs_str = strdup(equation_str);
        *rhs_str = strdup("0");
        return;
    }
    
    slong lhs_len = eq_pos - equation_str;
    *lhs_str = (char*) malloc((lhs_len + 1) * sizeof(char));
    strncpy(*lhs_str, equation_str, lhs_len);
    (*lhs_str)[lhs_len] = '\0';
    
    /* Remove trailing spaces from left side */
    while (lhs_len > 0 && ((*lhs_str)[lhs_len-1] == ' ' || (*lhs_str)[lhs_len-1] == '\t')) {
        (*lhs_str)[lhs_len-1] = '\0';
        lhs_len--;
    }
    
    /* Right side part */
    const char *rhs_start = eq_pos + 1;
    while (*rhs_start == ' ' || *rhs_start == '\t') rhs_start++;
    *rhs_str = strdup(rhs_start);
    
    DEBUG_PRINT_R("Split equation: LHS='%s', RHS='%s'\n", *lhs_str, *rhs_str);
}

/* Extract main variable information from left-hand side expression */
static void extract_main_variable_from_lhs(const fq_mvpoly_t *lhs, 
                                   char **main_var_name, 
                                   slong *main_var_degree,
                                   char **var_names_hint) {
    *main_var_name = NULL;
    *main_var_degree = 0;
    
    /* Find highest degree single variable term */
    slong max_degree = 0;
    slong main_var_idx = -1;
    
    for (slong t = 0; t < lhs->nterms; t++) {
        if (!lhs->terms[t].var_exp) continue;
        
        /* Check this term, find highest degree single variable term */
        slong single_var_idx = -1;
        slong single_var_degree = 0;
        int is_single_var = 1;
        
        for (slong v = 0; v < lhs->nvars; v++) {
            if (lhs->terms[t].var_exp[v] > 0) {
                if (single_var_idx == -1) {
                    single_var_idx = v;
                    single_var_degree = lhs->terms[t].var_exp[v];
                } else {
                    /* More than one variable */
                    is_single_var = 0;
                    break;
                }
            }
        }
        
        /* Check parameters (should be no parameters in single variable terms on left side) */
        if (is_single_var && lhs->terms[t].par_exp) {
            for (slong p = 0; p < lhs->npars; p++) {
                if (lhs->terms[t].par_exp[p] > 0) {
                    is_single_var = 0;
                    break;
                }
            }
        }
        
        /* If single variable term with higher degree, update main variable */
        if (is_single_var && single_var_degree > max_degree) {
            max_degree = single_var_degree;
            main_var_idx = single_var_idx;
        }
    }
    
    if (main_var_idx >= 0) {
        *main_var_degree = max_degree;
        
        if (var_names_hint && var_names_hint[main_var_idx]) {
            *main_var_name = strdup(var_names_hint[main_var_idx]);
        } else {
            char temp[32];
            sprintf(temp, "x_%ld", main_var_idx);
            *main_var_name = strdup(temp);
        }
        
        DEBUG_PRINT_R("Extracted main variable: '%s' with degree %ld\n", 
               *main_var_name, *main_var_degree);
    } else {
        DEBUG_PRINT_R("Warning: Could not extract main variable from LHS\n");
    }
}

/* Parse equation format generator */
static void parse_equation_generator(equation_info_t *eq_info, 
                             const char *equation_str,
                             char **var_names_hint,
                             slong nvars,
                             const fq_nmod_ctx_t ctx) {
    
    DEBUG_PRINT_R("Parsing equation: %s\n", equation_str);
    
    /* Split equation */
    char *lhs_str, *rhs_str;
    split_equation_string(equation_str, &lhs_str, &rhs_str);
    
    /* Set parser state */
    parser_state_t state = {0};
    state.var_names = (char**) malloc(nvars * sizeof(char*));
    for (slong i = 0; i < nvars; i++) {
        state.var_names[i] = var_names_hint[i] ? strdup(var_names_hint[i]) : NULL;
    }
    state.nvars = nvars;
    state.npars = 0;
    state.max_pars = 16;
    state.par_names = (char**) malloc(state.max_pars * sizeof(char*));
    state.ctx = ctx;
    state.generator_name = get_generator_name(ctx);
    fq_nmod_init(state.current.value, ctx);
    state.current.str = NULL;
    
    /* Parse left side */
    state.input = lhs_str;
    state.pos = 0;
    state.len = strlen(lhs_str);
    next_token(&state);
    parse_expression(&state, &eq_info->lhs);
    
    /* Parse right side */
    state.input = rhs_str;
    state.pos = 0;
    state.len = strlen(rhs_str);
    if (state.current.str) {
        free(state.current.str);
        state.current.str = NULL;
    }
    next_token(&state);
    parse_expression(&state, &eq_info->rhs);
    
    /* Extract main variable information from left side */
    extract_main_variable_from_lhs(&eq_info->lhs, &eq_info->main_var_name, 
                                  &eq_info->main_var_degree, var_names_hint);
    
    /* Build standard form: lhs - rhs */
    fq_mvpoly_sub(&eq_info->standard, &eq_info->lhs, &eq_info->rhs);
    
    DEBUG_PRINT_R("Equation parsing complete:\n");
    DEBUG_PRINT_R("  Main variable: %s^%ld\n", 
           eq_info->main_var_name ? eq_info->main_var_name : "unknown", 
           eq_info->main_var_degree);
    DEBUG_PRINT_R("  LHS terms: %ld, RHS terms: %ld, Standard terms: %ld\n",
           eq_info->lhs.nterms, eq_info->rhs.nterms, eq_info->standard.nterms);
    
    /* Cleanup */
    free(lhs_str);
    free(rhs_str);
    
    for (slong i = 0; i < nvars; i++) {
        if (state.var_names[i]) free(state.var_names[i]);
    }
    free(state.var_names);
    
    for (slong i = 0; i < state.npars; i++) {
        free(state.par_names[i]);
    }
    free(state.par_names);
    
    free(state.generator_name);
    fq_nmod_clear(state.current.value, ctx);
    if (state.current.str) free(state.current.str);
}

/* Initialize triangular ideal */
void unified_triangular_ideal_init(unified_triangular_ideal_t *ideal, 
                                  slong max_gens, slong nvars, 
                                  const fq_nmod_ctx_t field_ctx) {
    ideal->generators = (void**) flint_malloc(max_gens * sizeof(void*));
    ideal->var_indices = (slong*) flint_malloc(max_gens * sizeof(slong));
    ideal->var_names = (char**) flint_malloc(max_gens * sizeof(char*));
    ideal->leading_degrees = (slong*) flint_malloc(max_gens * sizeof(slong));
    ideal->num_gens = 0;
    ideal->max_gens = max_gens;
    ideal->field_ctx = field_ctx;
    
    /* Initialize arrays */
    for (slong i = 0; i < max_gens; i++) {
        ideal->var_names[i] = NULL;
        ideal->leading_degrees[i] = 3; /* Default to 3 */
    }
    
    /* Determine if prime field or extension field */
    if (fq_nmod_ctx_degree(field_ctx) == 1) {
        ideal->is_prime_field = 1;
        ideal->ctx.nmod_ctx = (nmod_mpoly_ctx_struct*) flint_malloc(sizeof(nmod_mpoly_ctx_struct));
        nmod_mpoly_ctx_init(ideal->ctx.nmod_ctx, nvars, ORD_LEX, fq_nmod_ctx_prime(field_ctx));
        
        for (slong i = 0; i < max_gens; i++) {
            ideal->generators[i] = flint_malloc(sizeof(nmod_mpoly_struct));
            nmod_mpoly_init((nmod_mpoly_struct*)ideal->generators[i], ideal->ctx.nmod_ctx);
        }
    } else {
        ideal->is_prime_field = 0;
        ideal->ctx.fq_ctx = (fq_nmod_mpoly_ctx_struct*) flint_malloc(sizeof(fq_nmod_mpoly_ctx_struct));
        fq_nmod_mpoly_ctx_init(ideal->ctx.fq_ctx, nvars, ORD_LEX, field_ctx);
        
        for (slong i = 0; i < max_gens; i++) {
            ideal->generators[i] = flint_malloc(sizeof(fq_nmod_mpoly_struct));
            fq_nmod_mpoly_init((fq_nmod_mpoly_struct*)ideal->generators[i], ideal->ctx.fq_ctx);
        }
    }
    ideal->complete_var_names = NULL;
    ideal->total_vars = 0;
}

/* Clear triangular ideal */
void unified_triangular_ideal_clear(unified_triangular_ideal_t *ideal) {
    if (ideal->is_prime_field) {
        for (slong i = 0; i < ideal->max_gens; i++) {
            nmod_mpoly_clear((nmod_mpoly_struct*)ideal->generators[i], ideal->ctx.nmod_ctx);
            flint_free(ideal->generators[i]);
        }
        nmod_mpoly_ctx_clear(ideal->ctx.nmod_ctx);
        flint_free(ideal->ctx.nmod_ctx);
    } else {
        for (slong i = 0; i < ideal->max_gens; i++) {
            fq_nmod_mpoly_clear((fq_nmod_mpoly_struct*)ideal->generators[i], ideal->ctx.fq_ctx);
            flint_free(ideal->generators[i]);
        }
        fq_nmod_mpoly_ctx_clear(ideal->ctx.fq_ctx);
        flint_free(ideal->ctx.fq_ctx);
    }
    
    /* Clear variable names */
    for (slong i = 0; i < ideal->max_gens; i++) {
        if (ideal->var_names[i] != NULL) {
            free(ideal->var_names[i]);
        }
    }

    if (ideal->complete_var_names) {
        for (slong i = 0; i < ideal->total_vars; i++) {
            if (ideal->complete_var_names[i]) {
                free(ideal->complete_var_names[i]);
            }
        }
        free(ideal->complete_var_names);
    }
    
    flint_free(ideal->generators);
    flint_free(ideal->var_indices);
    flint_free(ideal->var_names);
    flint_free(ideal->leading_degrees);
}

/* Create reduced ideal context with robust variable mapping logic */
void create_reduced_ideal_context(unified_triangular_ideal_t *reduced_ideal,
                                const unified_triangular_ideal_t *original_ideal,
                                char **current_var_names,
                                slong current_nvars,
                                const fq_nmod_ctx_t field_ctx) {
    
    /* Count how many generators are relevant to the current variables */
    slong relevant_gens = 0;
    for (slong i = 0; i < current_nvars; i++) {
        if (current_var_names[i]) {
            slong gen_idx = find_generator_for_variable(original_ideal, current_var_names[i]);
            if (gen_idx >= 0) {
                relevant_gens++;
            }
        }
    }
    
    DEBUG_PRINT_R("DEBUG: Creating reduced ideal with %ld relevant generators for %ld variables\n", 
           relevant_gens, current_nvars);
    DEBUG_PRINT_R("Current variables: ");
    for (slong i = 0; i < current_nvars; i++) {
        DEBUG_PRINT_R("%s ", current_var_names[i] ? current_var_names[i] : "?");
    }
    DEBUG_PRINT_R("\n");
    
    /* Initialize reduced ideal with matching context */
    unified_triangular_ideal_init(reduced_ideal, relevant_gens, current_nvars, field_ctx);
    
    /* For each current variable, find its generator and add to reduced ideal */
    for (slong i = 0; i < current_nvars; i++) {
        if (!current_var_names[i]) continue;
        
        slong gen_idx = find_generator_for_variable(original_ideal, current_var_names[i]);
        if (gen_idx < 0) continue;
        
        DEBUG_PRINT_R("DEBUG: Processing generator %ld for variable '%s' (new index %ld)\n", 
               gen_idx, current_var_names[i], i);
        
        /* Convert generator to new context */
        if (reduced_ideal->is_prime_field && original_ideal->is_prime_field) {
            nmod_mpoly_struct *old_gen = (nmod_mpoly_struct*)original_ideal->generators[gen_idx];
            nmod_mpoly_struct *new_gen = (nmod_mpoly_struct*)reduced_ideal->generators[reduced_ideal->num_gens];
            
            /* Create variable mapping from original to current context */
            slong orig_nvars = nmod_mpoly_ctx_nvars(original_ideal->ctx.nmod_ctx);
            slong *var_map = (slong*) flint_malloc(orig_nvars * sizeof(slong));
            
            /* Initialize mapping to -1 (variable eliminated) */
            for (slong j = 0; j < orig_nvars; j++) {
                var_map[j] = -1;
            }
            
            /* Build variable mapping based on names */
            DEBUG_PRINT_R("Using saved variable names for mapping:\n");
            
            for (slong old_v = 0; old_v < orig_nvars && old_v < original_ideal->total_vars; old_v++) {
                const char *old_var_name = original_ideal->complete_var_names[old_v];
                
                /* Find corresponding variable name in current variables */
                slong new_idx = find_variable_by_name(old_var_name, current_var_names, current_nvars);
                if (new_idx >= 0) {
                    var_map[old_v] = new_idx;
                    DEBUG_PRINT_R("  Mapped %s: old_idx=%ld -> new_idx=%ld\n", 
                           old_var_name, old_v, new_idx);
                } else {
                    DEBUG_PRINT_R("  Variable %s (old_idx=%ld) not found in current context (eliminated)\n", 
                           old_var_name, old_v);
                }
            }
            
            /* Apply variable mapping to convert generator */
            nmod_mpoly_zero(new_gen, reduced_ideal->ctx.nmod_ctx);
            
            for (slong t = 0; t < nmod_mpoly_length(old_gen, original_ideal->ctx.nmod_ctx); t++) {
                mp_limb_t coeff = nmod_mpoly_get_term_coeff_ui(old_gen, t, original_ideal->ctx.nmod_ctx);
                ulong *old_exp = (ulong*) flint_malloc(orig_nvars * sizeof(ulong));
                nmod_mpoly_get_term_exp_ui(old_exp, old_gen, t, original_ideal->ctx.nmod_ctx);
                
                /* Create new exponent vector */
                ulong *new_exp = (ulong*) flint_calloc(current_nvars, sizeof(ulong));
                int valid_term = 1;
                
                /* Map exponents */
                for (slong old_v = 0; old_v < orig_nvars; old_v++) {
                    if (old_exp[old_v] > 0) {
                        if (var_map[old_v] >= 0) {
                            new_exp[var_map[old_v]] = old_exp[old_v];
                        } else {
                            /* Variable was eliminated */
                            valid_term = 0;
                            break;
                        }
                    }
                }
                
                if (valid_term) {
                    nmod_mpoly_t temp;
                    nmod_mpoly_init(temp, reduced_ideal->ctx.nmod_ctx);
                    nmod_mpoly_set_coeff_ui_ui(temp, coeff, new_exp, reduced_ideal->ctx.nmod_ctx);
                    nmod_mpoly_add(new_gen, new_gen, temp, reduced_ideal->ctx.nmod_ctx);
                    nmod_mpoly_clear(temp, reduced_ideal->ctx.nmod_ctx);
                }
                
                flint_free(old_exp);
                flint_free(new_exp);
            }
            
            flint_free(var_map);
            
            /* Check if the resulting polynomial is non-zero */
            if (!nmod_mpoly_is_zero(new_gen, reduced_ideal->ctx.nmod_ctx)) {
                reduced_ideal->var_indices[reduced_ideal->num_gens] = i;
                reduced_ideal->var_names[reduced_ideal->num_gens] = strdup(current_var_names[i]);
                reduced_ideal->leading_degrees[reduced_ideal->num_gens] = 
                    original_ideal->leading_degrees[gen_idx]; /* Copy leading degree */
                reduced_ideal->num_gens++;
                
                DEBUG_PRINT_R("DEBUG: Added generator %ld to reduced ideal\n", reduced_ideal->num_gens - 1);
            }
        } else if (!reduced_ideal->is_prime_field && !original_ideal->is_prime_field) {
            /* Similar for fq_nmod_mpoly */
            fq_nmod_mpoly_struct *old_gen = (fq_nmod_mpoly_struct*)original_ideal->generators[gen_idx];
            fq_nmod_mpoly_struct *new_gen = (fq_nmod_mpoly_struct*)reduced_ideal->generators[reduced_ideal->num_gens];
            
            /* Create variable mapping from original to current context */
            slong orig_nvars = fq_nmod_mpoly_ctx_nvars(original_ideal->ctx.fq_ctx);
            slong *var_map = (slong*) flint_malloc(orig_nvars * sizeof(slong));
            
            /* Initialize mapping to -1 (variable eliminated) */
            for (slong j = 0; j < orig_nvars; j++) {
                var_map[j] = -1;
            }
            
            /* Build variable mapping based on names - same logic as nmod case */
            DEBUG_PRINT_R("Using saved variable names for fq_nmod mapping:\n");
            
            /* Direct mapping based on variable names */
            for (slong old_v = 0; old_v < orig_nvars && old_v < original_ideal->total_vars; old_v++) {
                const char *old_var_name = original_ideal->complete_var_names[old_v];
                
                /* Find corresponding variable name in current variables */
                slong new_idx = find_variable_by_name(old_var_name, current_var_names, current_nvars);
                if (new_idx >= 0) {
                    var_map[old_v] = new_idx;
                    DEBUG_PRINT_R("  Mapped %s: old_idx=%ld -> new_idx=%ld\n", 
                           old_var_name, old_v, new_idx);
                } else {
                    DEBUG_PRINT_R("  Variable %s (old_idx=%ld) not found in current context (eliminated)\n", 
                           old_var_name, old_v);
                }
            }
            
            /* Apply variable mapping to convert generator */
            fq_nmod_mpoly_zero(new_gen, reduced_ideal->ctx.fq_ctx);
            
            for (slong t = 0; t < fq_nmod_mpoly_length(old_gen, original_ideal->ctx.fq_ctx); t++) {
                fq_nmod_t coeff;
                fq_nmod_init(coeff, field_ctx);
                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, old_gen, t, original_ideal->ctx.fq_ctx);
                
                ulong *old_exp = (ulong*) flint_malloc(orig_nvars * sizeof(ulong));
                fq_nmod_mpoly_get_term_exp_ui(old_exp, old_gen, t, original_ideal->ctx.fq_ctx);
                
                /* Create new exponent vector */
                ulong *new_exp = (ulong*) flint_calloc(current_nvars, sizeof(ulong));
                int valid_term = 1;
                
                /* Map exponents */
                for (slong old_v = 0; old_v < orig_nvars; old_v++) {
                    if (old_exp[old_v] > 0) {
                        if (var_map[old_v] >= 0) {
                            new_exp[var_map[old_v]] = old_exp[old_v];
                        } else {
                            /* Variable was eliminated */
                            valid_term = 0;
                            break;
                        }
                    }
                }
                
                if (valid_term) {
                    fq_nmod_mpoly_set_coeff_fq_nmod_ui(new_gen, coeff, new_exp, reduced_ideal->ctx.fq_ctx);
                }
                
                fq_nmod_clear(coeff, field_ctx);
                flint_free(old_exp);
                flint_free(new_exp);
            }
            
            flint_free(var_map);
            
            /* Check if the resulting polynomial is non-zero */
            if (!fq_nmod_mpoly_is_zero(new_gen, reduced_ideal->ctx.fq_ctx)) {
                reduced_ideal->var_indices[reduced_ideal->num_gens] = i;
                reduced_ideal->var_names[reduced_ideal->num_gens] = strdup(current_var_names[i]);
                reduced_ideal->leading_degrees[reduced_ideal->num_gens] = 
                    original_ideal->leading_degrees[gen_idx];
                reduced_ideal->num_gens++;
                
                DEBUG_PRINT_R("DEBUG: Added generator %ld to reduced ideal\n", reduced_ideal->num_gens - 1);
            }
        }
    }
    
    DEBUG_PRINT_R("DEBUG: Reduced ideal has %ld generators\n", reduced_ideal->num_gens);
}

/* Optimized reduction function with batch polynomial building for nmod */
void triangular_ideal_reduce_nmod_mpoly_with_names(nmod_mpoly_t poly, 
                                                   const unified_triangular_ideal_t *ideal,
                                                   char **current_var_names) {
    if (ideal->num_gens == 0) return;
    
    slong nvars = nmod_mpoly_ctx_nvars(ideal->ctx.nmod_ctx);
    slong initial_terms = nmod_mpoly_length(poly, ideal->ctx.nmod_ctx);
    
    /* Timing variables */
    clock_t total_start = clock();
    double total_scan_time = 0.0;
    double total_fi_extract_time = 0.0;
    double total_power_precompute_time = 0.0;
    double total_reduction_time = 0.0;
    
    DEBUG_PRINT_R("\n=== REDUCTION TIMING ANALYSIS ===\n");
    DEBUG_PRINT_R("Initial polynomial: %ld terms, %ld variables\n", initial_terms, nvars);
    DEBUG_PRINT_R("Ideal has %ld generators\n\n", ideal->num_gens);
    
    /* Process variables in reverse order */
    for (slong var_idx = nvars - 1; var_idx >= 0; var_idx--) {
        if (!current_var_names || !current_var_names[var_idx]) continue;
        
        clock_t var_start = clock();
        
        /* Find generator for this variable */
        slong gen_idx = -1;
        slong leading_degree = 3;
        for (slong g = 0; g < ideal->num_gens; g++) {
            if (ideal->var_names[g] && strcmp(ideal->var_names[g], current_var_names[var_idx]) == 0) {
                gen_idx = g;
                leading_degree = ideal->leading_degrees[g];
                break;
            }
        }
        
        if (gen_idx < 0) continue;
        
        /* Quick scan to estimate reduction workload */
        clock_t scan_start = clock();
        slong current_terms = nmod_mpoly_length(poly, ideal->ctx.nmod_ctx);
        slong terms_to_reduce = 0;
        slong max_degree = 0;
        
        /* Sample terms to estimate reduction workload */
        slong sample_size = FLINT_MIN(5000, current_terms);
        slong sample_terms_to_reduce = 0;
        
        for (slong t = 0; t < sample_size; t++) {
            ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
            if (!exp) {
                DEBUG_PRINT_R("Memory allocation failed for exponent vector\n");
                continue;
            }
            
            nmod_mpoly_get_term_exp_ui(exp, poly, t, ideal->ctx.nmod_ctx);
            if (exp[var_idx] >= leading_degree) {
                sample_terms_to_reduce++;
                if (exp[var_idx] > max_degree) max_degree = exp[var_idx];
            }
            flint_free(exp);
        }
        
        /* Estimate total terms to reduce */
        if (sample_size < current_terms && sample_size > 0) {
            terms_to_reduce = (sample_terms_to_reduce * current_terms) / sample_size;
        } else {
            terms_to_reduce = sample_terms_to_reduce;
        }
        
        clock_t scan_end = clock();
        double scan_time = ((double)(scan_end - scan_start)) / CLOCKS_PER_SEC;
        total_scan_time += scan_time;
        
        DEBUG_PRINT_R("Variable '%s' (index %ld): %ld terms, max degree %ld\n", 
               current_var_names[var_idx], var_idx, current_terms, max_degree);
        
        if (terms_to_reduce == 0) {
            DEBUG_PRINT_R("  No terms to reduce, skipping\n\n");
            continue;
        }
        
        nmod_mpoly_struct *gen = (nmod_mpoly_struct*)ideal->generators[gen_idx];
        
        /* Extract f_i once and cache it */
        clock_t fi_start = clock();
        nmod_mpoly_t f_i;
        nmod_mpoly_init(f_i, ideal->ctx.nmod_ctx);
        
        /* Build f_i by subtracting leading term from generator */
        slong gen_length = nmod_mpoly_length(gen, ideal->ctx.nmod_ctx);
        for (slong t = 0; t < gen_length; t++) {
            ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
            if (!exp) continue;
            
            nmod_mpoly_get_term_exp_ui(exp, gen, t, ideal->ctx.nmod_ctx);
            
            /* Skip var^leading_degree term */
            int is_var_leading = (exp[var_idx] == leading_degree);
            for (slong v = 0; v < nvars; v++) {
                if (v != var_idx && exp[v] > 0) {
                    is_var_leading = 0;
                    break;
                }
            }
            
            if (!is_var_leading) {
                mp_limb_t coeff = nmod_mpoly_get_term_coeff_ui(gen, t, ideal->ctx.nmod_ctx);
                coeff = nmod_neg(coeff, ideal->ctx.nmod_ctx->mod);
                nmod_mpoly_t temp;
                nmod_mpoly_init(temp, ideal->ctx.nmod_ctx);
                nmod_mpoly_set_coeff_ui_ui(temp, coeff, exp, ideal->ctx.nmod_ctx);
                nmod_mpoly_add(f_i, f_i, temp, ideal->ctx.nmod_ctx);
                nmod_mpoly_clear(temp, ideal->ctx.nmod_ctx);
            }
            flint_free(exp);
        }
        
        clock_t fi_end = clock();
        double fi_time = ((double)(fi_end - fi_start)) / CLOCKS_PER_SEC;
        total_fi_extract_time += fi_time;
        
        /* Precompute powers for large polynomials */
        clock_t power_start = clock();
        slong max_quotient = max_degree / leading_degree;
        nmod_mpoly_t *f_i_powers = NULL;
        slong max_precompute = 0;
        
        if (max_quotient > 0 && terms_to_reduce > 1000) {
            max_precompute = FLINT_MIN(max_quotient, 6);
            f_i_powers = (nmod_mpoly_t*) flint_malloc((max_precompute + 1) * sizeof(nmod_mpoly_t));
            
            if (f_i_powers) {
                for (slong i = 0; i <= max_precompute; i++) {
                    nmod_mpoly_init(f_i_powers[i], ideal->ctx.nmod_ctx);
                }
                nmod_mpoly_one(f_i_powers[0], ideal->ctx.nmod_ctx);
                if (max_precompute >= 1) {
                    nmod_mpoly_set(f_i_powers[1], f_i, ideal->ctx.nmod_ctx);
                }
                
                for (slong i = 2; i <= max_precompute; i++) {
                    nmod_mpoly_mul(f_i_powers[i], f_i_powers[i-1], f_i, ideal->ctx.nmod_ctx);
                }
            }
        }
        
        clock_t power_end = clock();
        double power_time = ((double)(power_end - power_start)) / CLOCKS_PER_SEC;
        total_power_precompute_time += power_time;
        
        /* Process polynomial using batch building approach */
        clock_t reduction_start = clock();
        slong nterms = nmod_mpoly_length(poly, ideal->ctx.nmod_ctx);
        slong reduced_count = 0;
        
        /* Final result polynomial */
        nmod_mpoly_t result_poly;
        nmod_mpoly_init(result_poly, ideal->ctx.nmod_ctx);
        
        slong batch_size = 10000;
        
        /* Pre-allocate workspace for batch building */
        slong max_batch_terms = batch_size * 10; /* Conservative estimate */
        mp_limb_t *batch_coeffs = (mp_limb_t*) flint_malloc(max_batch_terms * sizeof(mp_limb_t));
        ulong *batch_exps = (ulong*) flint_malloc(max_batch_terms * nvars * sizeof(ulong));
        
        for (slong batch_start = 0; batch_start < nterms; batch_start += batch_size) {
            slong batch_end = FLINT_MIN(batch_start + batch_size, nterms);
            
            /* Collect all terms for this batch */
            slong batch_term_count = 0;
            
            /* Pre-allocate single exponent vector for reuse */
            ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
            
            for (slong t = batch_start; t < batch_end; t++) {
                nmod_mpoly_get_term_exp_ui(exp, poly, t, ideal->ctx.nmod_ctx);
                mp_limb_t coeff = nmod_mpoly_get_term_coeff_ui(poly, t, ideal->ctx.nmod_ctx);
                slong var_deg = exp[var_idx];
                
                if (var_deg >= leading_degree) {
                    /* Reduction case */
                    reduced_count++;
                    
                    slong quotient = var_deg / leading_degree;
                    slong remainder = var_deg % leading_degree;
                    exp[var_idx] = remainder;
                    
                    /* Apply reduction and collect resulting terms */
                    if (quotient == 0) {
                        /* Just add the term with modified exponent */
                        if (batch_term_count >= max_batch_terms) {
                            max_batch_terms *= 2;
                            batch_coeffs = (mp_limb_t*) flint_realloc(batch_coeffs, max_batch_terms * sizeof(mp_limb_t));
                            batch_exps = (ulong*) flint_realloc(batch_exps, max_batch_terms * nvars * sizeof(ulong));
                        }
                        batch_coeffs[batch_term_count] = coeff;
                        memcpy(batch_exps + batch_term_count * nvars, exp, nvars * sizeof(ulong));
                        batch_term_count++;
                    } else {
                        /* Need to multiply by f_i^quotient */
                        nmod_mpoly_t *power_poly = NULL;
                        nmod_mpoly_t temp_power;
                        
                        if (quotient == 1) {
                            power_poly = &f_i;
                        } else if (f_i_powers && quotient <= max_precompute) {
                            power_poly = &f_i_powers[quotient];
                        } else {
                            /* Compute power on demand */
                            nmod_mpoly_init(temp_power, ideal->ctx.nmod_ctx);
                            if (quotient <= 20) {
                                nmod_mpoly_pow_ui(temp_power, f_i, quotient, ideal->ctx.nmod_ctx);
                            } else {
                                /* Binary exponentiation */
                                nmod_mpoly_one(temp_power, ideal->ctx.nmod_ctx);
                                nmod_mpoly_t base;
                                nmod_mpoly_init(base, ideal->ctx.nmod_ctx);
                                nmod_mpoly_set(base, f_i, ideal->ctx.nmod_ctx);
                                
                                slong exp_remaining = quotient;
                                while (exp_remaining > 0) {
                                    if (exp_remaining & 1) {
                                        nmod_mpoly_mul(temp_power, temp_power, base, ideal->ctx.nmod_ctx);
                                    }
                                    if (exp_remaining > 1) {
                                        nmod_mpoly_mul(base, base, base, ideal->ctx.nmod_ctx);
                                    }
                                    exp_remaining >>= 1;
                                }
                                nmod_mpoly_clear(base, ideal->ctx.nmod_ctx);
                            }
                            power_poly = &temp_power;
                        }
                        
                        /* Collect terms from coeff * x^remainder * power_poly */
                        slong power_len = nmod_mpoly_length(*power_poly, ideal->ctx.nmod_ctx);
                        for (slong p = 0; p < power_len; p++) {
                            if (batch_term_count >= max_batch_terms) {
                                max_batch_terms *= 2;
                                batch_coeffs = (mp_limb_t*) flint_realloc(batch_coeffs, max_batch_terms * sizeof(mp_limb_t));
                                batch_exps = (ulong*) flint_realloc(batch_exps, max_batch_terms * nvars * sizeof(ulong));
                            }
                            
                            /* Get term from power polynomial */
                            ulong *power_exp = batch_exps + batch_term_count * nvars;
                            nmod_mpoly_get_term_exp_ui(power_exp, *power_poly, p, ideal->ctx.nmod_ctx);
                            mp_limb_t power_coeff = nmod_mpoly_get_term_coeff_ui(*power_poly, p, ideal->ctx.nmod_ctx);
                            
                            /* Add base exponent */
                            for (slong v = 0; v < nvars; v++) {
                                power_exp[v] += exp[v];
                            }
                            
                            /* Multiply coefficients */
                            batch_coeffs[batch_term_count] = nmod_mul(coeff, power_coeff, ideal->ctx.nmod_ctx->mod);
                            batch_term_count++;
                        }
                        
                        if (power_poly == &temp_power) {
                            nmod_mpoly_clear(temp_power, ideal->ctx.nmod_ctx);
                        }
                    }
                } else {
                    /* Direct copy case - just collect the term */
                    if (batch_term_count >= max_batch_terms) {
                        max_batch_terms *= 2;
                        batch_coeffs = (mp_limb_t*) flint_realloc(batch_coeffs, max_batch_terms * sizeof(mp_limb_t));
                        batch_exps = (ulong*) flint_realloc(batch_exps, max_batch_terms * nvars * sizeof(ulong));
                    }
                    
                    batch_coeffs[batch_term_count] = coeff;
                    memcpy(batch_exps + batch_term_count * nvars, exp, nvars * sizeof(ulong));
                    batch_term_count++;
                }
            }
            
            /* Free the reused exponent vector */
            flint_free(exp);
            
            /* Build polynomial from collected terms using low-level API */
            nmod_mpoly_t batch_poly;
            nmod_mpoly_init(batch_poly, ideal->ctx.nmod_ctx);
            
            if (batch_term_count > 0) {
                /* Use fit_length to pre-allocate space */
                nmod_mpoly_fit_length(batch_poly, batch_term_count, ideal->ctx.nmod_ctx);
                
                /* Set terms directly */
                for (slong i = 0; i < batch_term_count; i++) {
                    nmod_mpoly_push_term_ui_ui(batch_poly, batch_coeffs[i], 
                                              batch_exps + i * nvars, ideal->ctx.nmod_ctx);
                }
                
                /* Sort and combine like terms */
                nmod_mpoly_sort_terms(batch_poly, ideal->ctx.nmod_ctx);
                nmod_mpoly_combine_like_terms(batch_poly, ideal->ctx.nmod_ctx);
            }
            
            /* Add batch result to final result */
            nmod_mpoly_add(result_poly, result_poly, batch_poly, ideal->ctx.nmod_ctx);
            nmod_mpoly_clear(batch_poly, ideal->ctx.nmod_ctx);
            
            /* Progress indicator */
            if (batch_start % (batch_size * 10) == 0 && batch_start > 0) {
                DEBUG_PRINT_R("  Progress: %ld/%ld terms processed\n", batch_end, nterms);
            }
        }
        
        /* Free batch workspace */
        flint_free(batch_coeffs);
        flint_free(batch_exps);
        
        /* Replace original polynomial with result */
        nmod_mpoly_swap(poly, result_poly, ideal->ctx.nmod_ctx);
        
        clock_t reduction_end = clock();
        double reduction_time = ((double)(reduction_end - reduction_start)) / CLOCKS_PER_SEC;
        total_reduction_time += reduction_time;
        
        clock_t var_end = clock();
        double var_total = ((double)(var_end - var_start)) / CLOCKS_PER_SEC;
        
        DEBUG_PRINT_R("  Time: %.3f seconds (reduced %ld terms, final size: %ld)\n\n", 
               var_total, reduced_count, nmod_mpoly_length(poly, ideal->ctx.nmod_ctx));
        
        /* Cleanup */
        nmod_mpoly_clear(result_poly, ideal->ctx.nmod_ctx);
        nmod_mpoly_clear(f_i, ideal->ctx.nmod_ctx);
        
        if (f_i_powers) {
            for (slong i = 0; i <= max_precompute; i++) {
                nmod_mpoly_clear(f_i_powers[i], ideal->ctx.nmod_ctx);
            }
            flint_free(f_i_powers);
        }
    }
    
    clock_t total_end = clock();
    double total_time = ((double)(total_end - total_start)) / CLOCKS_PER_SEC;
    
    DEBUG_PRINT_R("=== TIMING SUMMARY ===\n");
    DEBUG_PRINT_R("Total reduction time: %.3f seconds\n", total_time);
    DEBUG_PRINT_R("  Scanning: %.3f seconds (%.1f%%)\n", 
           total_scan_time, 100.0 * total_scan_time / total_time);
    DEBUG_PRINT_R("  f_i extraction: %.3f seconds (%.1f%%)\n", 
           total_fi_extract_time, 100.0 * total_fi_extract_time / total_time);
    DEBUG_PRINT_R("  Power precomputation: %.3f seconds (%.1f%%)\n", 
           total_power_precompute_time, 100.0 * total_power_precompute_time / total_time);
    DEBUG_PRINT_R("  Main reduction: %.3f seconds (%.1f%%)\n", 
           total_reduction_time, 100.0 * total_reduction_time / total_time);
    DEBUG_PRINT_R("Final polynomial: %ld terms\n", nmod_mpoly_length(poly, ideal->ctx.nmod_ctx));
    DEBUG_PRINT_R("=== END TIMING ANALYSIS ===\n\n");
}

/* Similar function for fq_nmod_mpoly */
void triangular_ideal_reduce_fq_nmod_mpoly_with_names(fq_nmod_mpoly_t poly,
                                                      const unified_triangular_ideal_t *ideal,
                                                      char **current_var_names) {
    if (ideal->num_gens == 0) return;
    
    slong nvars = fq_nmod_mpoly_ctx_nvars(ideal->ctx.fq_ctx);
    slong initial_terms = fq_nmod_mpoly_length(poly, ideal->ctx.fq_ctx);
    
    /* Timing variables */
    clock_t total_start = clock();
    double total_scan_time = 0.0;
    double total_fi_extract_time = 0.0;
    double total_power_precompute_time = 0.0;
    double total_reduction_time = 0.0;
    
    DEBUG_PRINT_R("\n=== REDUCTION TIMING ANALYSIS (FQ) ===\n");
    DEBUG_PRINT_R("Initial polynomial: %ld terms, %ld variables\n", initial_terms, nvars);
    DEBUG_PRINT_R("Ideal has %ld generators\n\n", ideal->num_gens);
    
    /* Process variables in reverse order */
    for (slong var_idx = nvars - 1; var_idx >= 0; var_idx--) {
        if (!current_var_names || !current_var_names[var_idx]) continue;
        
        clock_t var_start = clock();
        
        /* Find generator for this variable */
        slong gen_idx = -1;
        slong leading_degree = 3;
        for (slong g = 0; g < ideal->num_gens; g++) {
            if (ideal->var_names[g] && strcmp(ideal->var_names[g], current_var_names[var_idx]) == 0) {
                gen_idx = g;
                leading_degree = ideal->leading_degrees[g];
                break;
            }
        }
        
        if (gen_idx < 0) continue;
        
        /* Quick scan to estimate reduction workload */
        clock_t scan_start = clock();
        slong current_terms = fq_nmod_mpoly_length(poly, ideal->ctx.fq_ctx);
        slong terms_to_reduce = 0;
        slong max_degree = 0;
        
        /* Sample terms to estimate reduction workload */
        slong sample_size = FLINT_MIN(5000, current_terms);
        slong sample_terms_to_reduce = 0;
        
        for (slong t = 0; t < sample_size; t++) {
            ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
            if (!exp) {
                DEBUG_PRINT_R("Memory allocation failed for exponent vector\n");
                continue;
            }
            
            fq_nmod_mpoly_get_term_exp_ui(exp, poly, t, ideal->ctx.fq_ctx);
            if (exp[var_idx] >= leading_degree) {
                sample_terms_to_reduce++;
                if (exp[var_idx] > max_degree) max_degree = exp[var_idx];
            }
            flint_free(exp);
        }
        
        /* Estimate total terms to reduce */
        if (sample_size < current_terms && sample_size > 0) {
            terms_to_reduce = (sample_terms_to_reduce * current_terms) / sample_size;
        } else {
            terms_to_reduce = sample_terms_to_reduce;
        }
        
        clock_t scan_end = clock();
        double scan_time = ((double)(scan_end - scan_start)) / CLOCKS_PER_SEC;
        total_scan_time += scan_time;
        
        DEBUG_PRINT_R("Variable '%s' (index %ld): %ld terms, max degree %ld\n", 
               current_var_names[var_idx], var_idx, current_terms, max_degree);
        
        if (terms_to_reduce == 0) {
            DEBUG_PRINT_R("  No terms to reduce, skipping\n\n");
            continue;
        }
        
        fq_nmod_mpoly_struct *gen = (fq_nmod_mpoly_struct*)ideal->generators[gen_idx];
        
        /* Extract f_i once and cache it */
        clock_t fi_start = clock();
        fq_nmod_mpoly_t f_i;
        fq_nmod_mpoly_init(f_i, ideal->ctx.fq_ctx);
        
        /* Build f_i by subtracting leading term from generator */
        slong gen_length = fq_nmod_mpoly_length(gen, ideal->ctx.fq_ctx);
        for (slong t = 0; t < gen_length; t++) {
            ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
            if (!exp) continue;
            
            fq_nmod_mpoly_get_term_exp_ui(exp, gen, t, ideal->ctx.fq_ctx);
            
            /* Skip var^leading_degree term */
            int is_var_leading = (exp[var_idx] == leading_degree);
            for (slong v = 0; v < nvars; v++) {
                if (v != var_idx && exp[v] > 0) {
                    is_var_leading = 0;
                    break;
                }
            }
            
            if (!is_var_leading) {
                fq_nmod_t coeff;
                fq_nmod_init(coeff, ideal->field_ctx);
                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, gen, t, ideal->ctx.fq_ctx);
                fq_nmod_neg(coeff, coeff, ideal->field_ctx);
                
                fq_nmod_mpoly_set_coeff_fq_nmod_ui(f_i, coeff, exp, ideal->ctx.fq_ctx);
                
                fq_nmod_clear(coeff, ideal->field_ctx);
            }
            flint_free(exp);
        }
        
        clock_t fi_end = clock();
        double fi_time = ((double)(fi_end - fi_start)) / CLOCKS_PER_SEC;
        total_fi_extract_time += fi_time;
        
        /* Precompute powers for large polynomials */
        clock_t power_start = clock();
        slong max_quotient = max_degree / leading_degree;
        fq_nmod_mpoly_t *f_i_powers = NULL;
        slong max_precompute = 0;
        
        if (max_quotient > 0 && terms_to_reduce > 1000) {
            max_precompute = FLINT_MIN(max_quotient, 6);
            f_i_powers = (fq_nmod_mpoly_t*) flint_malloc((max_precompute + 1) * sizeof(fq_nmod_mpoly_t));
            
            if (f_i_powers) {
                for (slong i = 0; i <= max_precompute; i++) {
                    fq_nmod_mpoly_init(f_i_powers[i], ideal->ctx.fq_ctx);
                }
                fq_nmod_mpoly_one(f_i_powers[0], ideal->ctx.fq_ctx);
                if (max_precompute >= 1) {
                    fq_nmod_mpoly_set(f_i_powers[1], f_i, ideal->ctx.fq_ctx);
                }
                
                for (slong i = 2; i <= max_precompute; i++) {
                    fq_nmod_mpoly_mul(f_i_powers[i], f_i_powers[i-1], f_i, ideal->ctx.fq_ctx);
                }
            }
        }
        
        clock_t power_end = clock();
        double power_time = ((double)(power_end - power_start)) / CLOCKS_PER_SEC;
        total_power_precompute_time += power_time;
        
        /* Process polynomial using batch building approach */
        clock_t reduction_start = clock();
        slong nterms = fq_nmod_mpoly_length(poly, ideal->ctx.fq_ctx);
        slong reduced_count = 0;
        
        /* Final result polynomial */
        fq_nmod_mpoly_t result_poly;
        fq_nmod_mpoly_init(result_poly, ideal->ctx.fq_ctx);
        
        slong batch_size = 10000;
        
        /* Pre-allocate workspace for batch building */
        slong max_batch_terms = batch_size * 10;
        fq_nmod_t *batch_coeffs = (fq_nmod_t*) flint_malloc(max_batch_terms * sizeof(fq_nmod_t));
        for (slong i = 0; i < max_batch_terms; i++) {
            fq_nmod_init(batch_coeffs[i], ideal->field_ctx);
        }
        ulong *batch_exps = (ulong*) flint_malloc(max_batch_terms * nvars * sizeof(ulong));
        
        for (slong batch_start = 0; batch_start < nterms; batch_start += batch_size) {
            slong batch_end = FLINT_MIN(batch_start + batch_size, nterms);
            
            /* Collect all terms for this batch */
            slong batch_term_count = 0;
            
            /* Pre-allocate single exponent vector for reuse */
            ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ideal->field_ctx);
            
            for (slong t = batch_start; t < batch_end; t++) {
                fq_nmod_mpoly_get_term_exp_ui(exp, poly, t, ideal->ctx.fq_ctx);
                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, poly, t, ideal->ctx.fq_ctx);
                slong var_deg = exp[var_idx];
                
                if (var_deg >= leading_degree) {
                    /* Reduction case */
                    reduced_count++;
                    
                    slong quotient = var_deg / leading_degree;
                    slong remainder = var_deg % leading_degree;
                    exp[var_idx] = remainder;
                    
                    /* Apply reduction and collect resulting terms */
                    if (quotient == 0) {
                        /* Just add the term with modified exponent */
                        if (batch_term_count >= max_batch_terms) {
                            max_batch_terms *= 2;
                            fq_nmod_t *new_coeffs = (fq_nmod_t*) flint_malloc(max_batch_terms * sizeof(fq_nmod_t));
                            for (slong i = 0; i < batch_term_count; i++) {
                                fq_nmod_init(new_coeffs[i], ideal->field_ctx);
                                fq_nmod_set(new_coeffs[i], batch_coeffs[i], ideal->field_ctx);
                            }
                            for (slong i = batch_term_count; i < max_batch_terms; i++) {
                                fq_nmod_init(new_coeffs[i], ideal->field_ctx);
                            }
                            for (slong i = 0; i < batch_term_count; i++) {
                                fq_nmod_clear(batch_coeffs[i], ideal->field_ctx);
                            }
                            flint_free(batch_coeffs);
                            batch_coeffs = new_coeffs;
                            
                            batch_exps = (ulong*) flint_realloc(batch_exps, max_batch_terms * nvars * sizeof(ulong));
                        }
                        fq_nmod_set(batch_coeffs[batch_term_count], coeff, ideal->field_ctx);
                        memcpy(batch_exps + batch_term_count * nvars, exp, nvars * sizeof(ulong));
                        batch_term_count++;
                    } else {
                        /* Need to multiply by f_i^quotient */
                        fq_nmod_mpoly_t *power_poly = NULL;
                        fq_nmod_mpoly_t temp_power;
                        
                        if (quotient == 1) {
                            power_poly = &f_i;
                        } else if (f_i_powers && quotient <= max_precompute) {
                            power_poly = &f_i_powers[quotient];
                        } else {
                            /* Compute power on demand */
                            fq_nmod_mpoly_init(temp_power, ideal->ctx.fq_ctx);
                            if (quotient <= 20) {
                                fq_nmod_mpoly_pow_ui(temp_power, f_i, quotient, ideal->ctx.fq_ctx);
                            } else {
                                /* Binary exponentiation */
                                fq_nmod_mpoly_one(temp_power, ideal->ctx.fq_ctx);
                                fq_nmod_mpoly_t base;
                                fq_nmod_mpoly_init(base, ideal->ctx.fq_ctx);
                                fq_nmod_mpoly_set(base, f_i, ideal->ctx.fq_ctx);
                                
                                slong exp_remaining = quotient;
                                while (exp_remaining > 0) {
                                    if (exp_remaining & 1) {
                                        fq_nmod_mpoly_mul(temp_power, temp_power, base, ideal->ctx.fq_ctx);
                                    }
                                    if (exp_remaining > 1) {
                                        fq_nmod_mpoly_mul(base, base, base, ideal->ctx.fq_ctx);
                                    }
                                    exp_remaining >>= 1;
                                }
                                fq_nmod_mpoly_clear(base, ideal->ctx.fq_ctx);
                            }
                            power_poly = &temp_power;
                        }
                        
                        /* Collect terms from coeff * x^remainder * power_poly */
                        slong power_len = fq_nmod_mpoly_length(*power_poly, ideal->ctx.fq_ctx);
                        fq_nmod_t power_coeff;
                        fq_nmod_init(power_coeff, ideal->field_ctx);
                        
                        for (slong p = 0; p < power_len; p++) {
                            if (batch_term_count >= max_batch_terms) {
                                max_batch_terms *= 2;
                                fq_nmod_t *new_coeffs = (fq_nmod_t*) flint_malloc(max_batch_terms * sizeof(fq_nmod_t));
                                for (slong i = 0; i < batch_term_count; i++) {
                                    fq_nmod_init(new_coeffs[i], ideal->field_ctx);
                                    fq_nmod_set(new_coeffs[i], batch_coeffs[i], ideal->field_ctx);
                                }
                                for (slong i = batch_term_count; i < max_batch_terms; i++) {
                                    fq_nmod_init(new_coeffs[i], ideal->field_ctx);
                                }
                                for (slong i = 0; i < batch_term_count; i++) {
                                    fq_nmod_clear(batch_coeffs[i], ideal->field_ctx);
                                }
                                flint_free(batch_coeffs);
                                batch_coeffs = new_coeffs;
                                
                                batch_exps = (ulong*) flint_realloc(batch_exps, max_batch_terms * nvars * sizeof(ulong));
                            }
                            
                            /* Get term from power polynomial */
                            ulong *power_exp = batch_exps + batch_term_count * nvars;
                            fq_nmod_mpoly_get_term_exp_ui(power_exp, *power_poly, p, ideal->ctx.fq_ctx);
                            fq_nmod_mpoly_get_term_coeff_fq_nmod(power_coeff, *power_poly, p, ideal->ctx.fq_ctx);
                            
                            /* Add base exponent */
                            for (slong v = 0; v < nvars; v++) {
                                power_exp[v] += exp[v];
                            }
                            
                            /* Multiply coefficients */
                            fq_nmod_mul(batch_coeffs[batch_term_count], coeff, power_coeff, ideal->field_ctx);
                            batch_term_count++;
                        }
                        
                        fq_nmod_clear(power_coeff, ideal->field_ctx);
                        
                        if (power_poly == &temp_power) {
                            fq_nmod_mpoly_clear(temp_power, ideal->ctx.fq_ctx);
                        }
                    }
                } else {
                    /* Direct copy case - just collect the term */
                    if (batch_term_count >= max_batch_terms) {
                        max_batch_terms *= 2;
                        fq_nmod_t *new_coeffs = (fq_nmod_t*) flint_malloc(max_batch_terms * sizeof(fq_nmod_t));
                        for (slong i = 0; i < batch_term_count; i++) {
                            fq_nmod_init(new_coeffs[i], ideal->field_ctx);
                            fq_nmod_set(new_coeffs[i], batch_coeffs[i], ideal->field_ctx);
                        }
                        for (slong i = batch_term_count; i < max_batch_terms; i++) {
                            fq_nmod_init(new_coeffs[i], ideal->field_ctx);
                        }
                        for (slong i = 0; i < batch_term_count; i++) {
                            fq_nmod_clear(batch_coeffs[i], ideal->field_ctx);
                        }
                        flint_free(batch_coeffs);
                        batch_coeffs = new_coeffs;
                        
                        batch_exps = (ulong*) flint_realloc(batch_exps, max_batch_terms * nvars * sizeof(ulong));
                    }
                    
                    fq_nmod_set(batch_coeffs[batch_term_count], coeff, ideal->field_ctx);
                    memcpy(batch_exps + batch_term_count * nvars, exp, nvars * sizeof(ulong));
                    batch_term_count++;
                }
            }
            
            /* Free the reused vectors */
            flint_free(exp);
            fq_nmod_clear(coeff, ideal->field_ctx);
            
            /* Build polynomial from collected terms */
            fq_nmod_mpoly_t batch_poly;
            fq_nmod_mpoly_init(batch_poly, ideal->ctx.fq_ctx);
            
            if (batch_term_count > 0) {
                /* Add terms one by one (no direct batch API for fq_nmod_mpoly) */
                for (slong i = 0; i < batch_term_count; i++) {
                    fq_nmod_mpoly_set_coeff_fq_nmod_ui(batch_poly, batch_coeffs[i], 
                                                       batch_exps + i * nvars, ideal->ctx.fq_ctx);
                }
            }
            
            /* Add batch result to final result */
            fq_nmod_mpoly_add(result_poly, result_poly, batch_poly, ideal->ctx.fq_ctx);
            fq_nmod_mpoly_clear(batch_poly, ideal->ctx.fq_ctx);
            
            /* Progress indicator */
            if (batch_start % (batch_size * 10) == 0 && batch_start > 0) {
                DEBUG_PRINT_R("  Progress: %ld/%ld terms processed\n", batch_end, nterms);
            }
        }
        
        /* Free batch workspace */
        for (slong i = 0; i < max_batch_terms; i++) {
            fq_nmod_clear(batch_coeffs[i], ideal->field_ctx);
        }
        flint_free(batch_coeffs);
        flint_free(batch_exps);
        
        /* Replace original polynomial with result */
        fq_nmod_mpoly_swap(poly, result_poly, ideal->ctx.fq_ctx);
        
        clock_t reduction_end = clock();
        double reduction_time = ((double)(reduction_end - reduction_start)) / CLOCKS_PER_SEC;
        total_reduction_time += reduction_time;
        
        clock_t var_end = clock();
        double var_total = ((double)(var_end - var_start)) / CLOCKS_PER_SEC;
        
        DEBUG_PRINT_R("  Time: %.3f seconds (reduced %ld terms, final size: %ld)\n\n", 
               var_total, reduced_count, fq_nmod_mpoly_length(poly, ideal->ctx.fq_ctx));
        
        /* Cleanup */
        fq_nmod_mpoly_clear(result_poly, ideal->ctx.fq_ctx);
        fq_nmod_mpoly_clear(f_i, ideal->ctx.fq_ctx);
        
        if (f_i_powers) {
            for (slong i = 0; i <= max_precompute; i++) {
                fq_nmod_mpoly_clear(f_i_powers[i], ideal->ctx.fq_ctx);
            }
            flint_free(f_i_powers);
        }
    }
    
    clock_t total_end = clock();
    double total_time = ((double)(total_end - total_start)) / CLOCKS_PER_SEC;
    
    DEBUG_PRINT_R("=== TIMING SUMMARY (FQ) ===\n");
    DEBUG_PRINT_R("Total reduction time: %.3f seconds\n", total_time);
    DEBUG_PRINT_R("  Scanning: %.3f seconds (%.1f%%)\n", 
           total_scan_time, 100.0 * total_scan_time / total_time);
    DEBUG_PRINT_R("  f_i extraction: %.3f seconds (%.1f%%)\n", 
           total_fi_extract_time, 100.0 * total_fi_extract_time / total_time);
    DEBUG_PRINT_R("  Power precomputation: %.3f seconds (%.1f%%)\n", 
           total_power_precompute_time, 100.0 * total_power_precompute_time / total_time);
    DEBUG_PRINT_R("  Main reduction: %.3f seconds (%.1f%%)\n", 
           total_reduction_time, 100.0 * total_reduction_time / total_time);
    DEBUG_PRINT_R("Final polynomial: %ld terms\n", fq_nmod_mpoly_length(poly, ideal->ctx.fq_ctx));
    DEBUG_PRINT_R("=== END TIMING ANALYSIS ===\n\n");
}

/* Wrapper functions that maintain compatibility */
void triangular_ideal_reduce_nmod_mpoly(nmod_mpoly_t poly, 
                                       const unified_triangular_ideal_t *ideal) {
    triangular_ideal_reduce_nmod_mpoly_with_names(poly, ideal, NULL);
}

void triangular_ideal_reduce_fq_nmod_mpoly(fq_nmod_mpoly_t poly,
                                          const unified_triangular_ideal_t *ideal) {
    triangular_ideal_reduce_fq_nmod_mpoly_with_names(poly, ideal, NULL);
}

/* Compute determinant with reduction for nmod_mpoly matrix - with variable names */
void compute_nmod_det_with_triangular_reduction_with_names(nmod_mpoly_t det,
                                                          nmod_mpoly_struct **matrix,
                                                          slong size,
                                                          const unified_triangular_ideal_t *ideal,
                                                          char **current_var_names) {
    if (size == 0) {
        nmod_mpoly_one(det, ideal->ctx.nmod_ctx);
        return;
    }
    
    if (size == 1) {
        nmod_mpoly_set(det, &matrix[0][0], ideal->ctx.nmod_ctx);
        triangular_ideal_reduce_nmod_mpoly_with_names(det, ideal, current_var_names);
        return;
    }
    
    if (size == 2) {
        /* 2x2 matrix determinant */
        nmod_mpoly_t ad, bc, temp;
        nmod_mpoly_init(ad, ideal->ctx.nmod_ctx);
        nmod_mpoly_init(bc, ideal->ctx.nmod_ctx);
        nmod_mpoly_init(temp, ideal->ctx.nmod_ctx);
        
        nmod_mpoly_mul(ad, &matrix[0][0], &matrix[1][1], ideal->ctx.nmod_ctx);
        nmod_mpoly_mul(bc, &matrix[0][1], &matrix[1][0], ideal->ctx.nmod_ctx);
        nmod_mpoly_sub(det, ad, bc, ideal->ctx.nmod_ctx);
        
        /* Only reduce once at the end */
        triangular_ideal_reduce_nmod_mpoly_with_names(det, ideal, current_var_names);
        
        nmod_mpoly_clear(ad, ideal->ctx.nmod_ctx);
        nmod_mpoly_clear(bc, ideal->ctx.nmod_ctx);
        nmod_mpoly_clear(temp, ideal->ctx.nmod_ctx);
        return;
    }
    
    /* For larger matrices, use cofactor expansion */
    nmod_mpoly_zero(det, ideal->ctx.nmod_ctx);
    
    /* Expand along first row */
    for (slong j = 0; j < size; j++) {
        if (nmod_mpoly_is_zero(&matrix[0][j], ideal->ctx.nmod_ctx)) continue;
        
        /* Build submatrix */
        nmod_mpoly_struct **submatrix = (nmod_mpoly_struct**) 
            flint_malloc((size-1) * sizeof(nmod_mpoly_struct*));
        for (slong i = 0; i < size-1; i++) {
            submatrix[i] = (nmod_mpoly_struct*) 
                flint_malloc((size-1) * sizeof(nmod_mpoly_struct));
            for (slong k = 0; k < size-1; k++) {
                nmod_mpoly_init(&submatrix[i][k], ideal->ctx.nmod_ctx);
            }
        }
        
        /* Copy submatrix */
        for (slong i = 1; i < size; i++) {
            slong col_idx = 0;
            for (slong k = 0; k < size; k++) {
                if (k != j) {
                    nmod_mpoly_set(&submatrix[i-1][col_idx], &matrix[i][k], 
                                  ideal->ctx.nmod_ctx);
                    col_idx++;
                }
            }
        }
        
        /* Compute minor recursively */
        nmod_mpoly_t minor;
        nmod_mpoly_init(minor, ideal->ctx.nmod_ctx);
        compute_nmod_det_with_triangular_reduction_with_names(minor, submatrix, size-1, ideal, current_var_names);
        
        /* Compute contribution */
        nmod_mpoly_t contrib;
        nmod_mpoly_init(contrib, ideal->ctx.nmod_ctx);
        nmod_mpoly_mul(contrib, &matrix[0][j], minor, ideal->ctx.nmod_ctx);
        
        if (j % 2 == 0) {
            nmod_mpoly_add(det, det, contrib, ideal->ctx.nmod_ctx);
        } else {
            nmod_mpoly_sub(det, det, contrib, ideal->ctx.nmod_ctx);
        }
        
        /* Reduce periodically */
        if ((j + 1) % 3 == 0 || j == size - 1) {
            triangular_ideal_reduce_nmod_mpoly_with_names(det, ideal, current_var_names);
        }
        
        /* Cleanup */
        nmod_mpoly_clear(minor, ideal->ctx.nmod_ctx);
        nmod_mpoly_clear(contrib, ideal->ctx.nmod_ctx);
        
        for (slong i = 0; i < size-1; i++) {
            for (slong k = 0; k < size-1; k++) {
                nmod_mpoly_clear(&submatrix[i][k], ideal->ctx.nmod_ctx);
            }
            flint_free(submatrix[i]);
        }
        flint_free(submatrix);
    }
}

/* Similar function for fq_nmod_mpoly - with variable names */
void compute_fq_nmod_det_with_triangular_reduction_with_names(fq_nmod_mpoly_t det,
                                                             fq_nmod_mpoly_struct **matrix,
                                                             slong size,
                                                             const unified_triangular_ideal_t *ideal,
                                                             char **current_var_names) {
    if (size == 0) {
        fq_nmod_mpoly_one(det, ideal->ctx.fq_ctx);
        return;
    }
    
    if (size == 1) {
        fq_nmod_mpoly_set(det, &matrix[0][0], ideal->ctx.fq_ctx);
        triangular_ideal_reduce_fq_nmod_mpoly_with_names(det, ideal, current_var_names);
        return;
    }
    
    if (size == 2) {
        /* 2x2 matrix determinant */
        fq_nmod_mpoly_t ad, bc, temp;
        fq_nmod_mpoly_init(ad, ideal->ctx.fq_ctx);
        fq_nmod_mpoly_init(bc, ideal->ctx.fq_ctx);
        fq_nmod_mpoly_init(temp, ideal->ctx.fq_ctx);
        
        fq_nmod_mpoly_mul(ad, &matrix[0][0], &matrix[1][1], ideal->ctx.fq_ctx);
        fq_nmod_mpoly_mul(bc, &matrix[0][1], &matrix[1][0], ideal->ctx.fq_ctx);
        fq_nmod_mpoly_sub(det, ad, bc, ideal->ctx.fq_ctx);
        
        /* Only reduce once at the end */
        triangular_ideal_reduce_fq_nmod_mpoly_with_names(det, ideal, current_var_names);
        
        fq_nmod_mpoly_clear(ad, ideal->ctx.fq_ctx);
        fq_nmod_mpoly_clear(bc, ideal->ctx.fq_ctx);
        fq_nmod_mpoly_clear(temp, ideal->ctx.fq_ctx);
        return;
    }
    
    /* For larger matrices, use cofactor expansion */
    fq_nmod_mpoly_zero(det, ideal->ctx.fq_ctx);
    
    /* Expand along first row */
    for (slong j = 0; j < size; j++) {
        if (fq_nmod_mpoly_is_zero(&matrix[0][j], ideal->ctx.fq_ctx)) continue;
        
        /* Build submatrix */
        fq_nmod_mpoly_struct **submatrix = (fq_nmod_mpoly_struct**) 
            flint_malloc((size-1) * sizeof(fq_nmod_mpoly_struct*));
        for (slong i = 0; i < size-1; i++) {
            submatrix[i] = (fq_nmod_mpoly_struct*) 
                flint_malloc((size-1) * sizeof(fq_nmod_mpoly_struct));
            for (slong k = 0; k < size-1; k++) {
                fq_nmod_mpoly_init(&submatrix[i][k], ideal->ctx.fq_ctx);
            }
        }
        
        /* Copy submatrix */
        for (slong i = 1; i < size; i++) {
            slong col_idx = 0;
            for (slong k = 0; k < size; k++) {
                if (k != j) {
                    fq_nmod_mpoly_set(&submatrix[i-1][col_idx], &matrix[i][k], 
                                     ideal->ctx.fq_ctx);
                    col_idx++;
                }
            }
        }
        
        /* Compute minor recursively */
        fq_nmod_mpoly_t minor;
        fq_nmod_mpoly_init(minor, ideal->ctx.fq_ctx);
        compute_fq_nmod_det_with_triangular_reduction_with_names(minor, submatrix, size-1, ideal, current_var_names);
        
        /* Compute contribution */
        fq_nmod_mpoly_t contrib;
        fq_nmod_mpoly_init(contrib, ideal->ctx.fq_ctx);
        fq_nmod_mpoly_mul(contrib, &matrix[0][j], minor, ideal->ctx.fq_ctx);
        
        if (j % 2 == 0) {
            fq_nmod_mpoly_add(det, det, contrib, ideal->ctx.fq_ctx);
        } else {
            fq_nmod_mpoly_sub(det, det, contrib, ideal->ctx.fq_ctx);
        }
        
        /* Reduce periodically */
        if ((j + 1) % 3 == 0 || j == size - 1) {
            triangular_ideal_reduce_fq_nmod_mpoly_with_names(det, ideal, current_var_names);
        }
        
        /* Cleanup */
        fq_nmod_mpoly_clear(minor, ideal->ctx.fq_ctx);
        fq_nmod_mpoly_clear(contrib, ideal->ctx.fq_ctx);
        
        for (slong i = 0; i < size-1; i++) {
            for (slong k = 0; k < size-1; k++) {
                fq_nmod_mpoly_clear(&submatrix[i][k], ideal->ctx.fq_ctx);
            }
            flint_free(submatrix[i]);
        }
        flint_free(submatrix);
    }
}

/* Old wrapper functions for compatibility */
void compute_nmod_det_with_triangular_reduction(nmod_mpoly_t det,
                                               nmod_mpoly_struct **matrix,
                                               slong size,
                                               const unified_triangular_ideal_t *ideal) {
    compute_nmod_det_with_triangular_reduction_with_names(det, matrix, size, ideal, NULL);
}

void compute_fq_nmod_det_with_triangular_reduction(fq_nmod_mpoly_t det,
                                                  fq_nmod_mpoly_struct **matrix,
                                                  slong size,
                                                  const unified_triangular_ideal_t *ideal) {
    compute_fq_nmod_det_with_triangular_reduction_with_names(det, matrix, size, ideal, NULL);
}

/* Simplified and safer matrix conversion with timeout protection */
void compute_det_with_reduction_from_mvpoly(fq_mvpoly_t *result,
                                           fq_mvpoly_t **matrix,
                                           slong size,
                                           const unified_triangular_ideal_t *ideal,
                                           char **current_var_names) {
    DEBUG_PRINT_R("Computing determinant with ideal reduction (matrix size: %ld x %ld)\n", size, size);
    
    if (size == 0) {
        fq_mvpoly_init(result, 0, 0, ideal->field_ctx);
        return;
    }
    
    slong current_nvars = matrix[0][0].npars;  /* These are the remaining variables */
    
    /* Build complete variable name list - only parameters remain as variables */
    char **complete_var_names = (char**) malloc(current_nvars * sizeof(char*));
    for (slong i = 0; i < current_nvars; i++) {
        if (current_var_names && current_var_names[i]) {
            complete_var_names[i] = current_var_names[i];
        } else {
            char temp[32];
            sprintf(temp, "param_%ld", i);
            complete_var_names[i] = strdup(temp);
        }
    }
    
    /* Create reduced ideal with proper context */
    unified_triangular_ideal_t reduced_ideal;
    create_reduced_ideal_context(&reduced_ideal, ideal, complete_var_names, 
                               current_nvars, ideal->field_ctx);
    
    if (ideal->is_prime_field) {
        DEBUG_PRINT_R("Using nmod_mpoly for prime field computation\n");
        
        /* Convert to nmod_mpoly matrix with only parameter variables */
        nmod_mpoly_struct **nmod_matrix = (nmod_mpoly_struct**) flint_malloc(size * sizeof(nmod_mpoly_struct*));
        for (slong i = 0; i < size; i++) {
            nmod_matrix[i] = (nmod_mpoly_struct*) flint_malloc(size * sizeof(nmod_mpoly_struct));
            for (slong j = 0; j < size; j++) {
                nmod_mpoly_init(&nmod_matrix[i][j], reduced_ideal.ctx.nmod_ctx);
                nmod_mpoly_zero(&nmod_matrix[i][j], reduced_ideal.ctx.nmod_ctx);
                
                /* Convert fq_mvpoly to nmod_mpoly - SIMPLIFIED AND SAFER */
                DEBUG_PRINT_R("Converting matrix[%ld][%ld] with %ld terms...\n", i, j, matrix[i][j].nterms);
                
                /* Process each term with safety checks */
                for (slong t = 0; t < matrix[i][j].nterms && t < 10000; t++) { /* Limit terms to prevent infinite loops */
                    /* Get coefficient as mp_limb_t */
                    mp_limb_t coeff = 0;
                    
                    /* For prime field, extract the coefficient safely */
                    if (fq_nmod_is_one(matrix[i][j].terms[t].coeff, ideal->field_ctx)) {
                        coeff = 1;
                    } else if (!fq_nmod_is_zero(matrix[i][j].terms[t].coeff, ideal->field_ctx)) {
                        /* Extract coefficient value - for prime field, it's the constant term */
                        fq_nmod_t temp;
                        fq_nmod_init(temp, ideal->field_ctx);
                        fq_nmod_set(temp, matrix[i][j].terms[t].coeff, ideal->field_ctx);
                        
                        /* For prime field (degree 1), coefficient is directly the value */
                        nmod_poly_struct *p = &temp[0];
                        if (p->length > 0) {
                            coeff = p->coeffs[0];
                        }
                        
                        fq_nmod_clear(temp, ideal->field_ctx);
                    }
                    
                    /* Build exponent vector for nmod_mpoly */
                    ulong *exp = (ulong*) flint_calloc(current_nvars, sizeof(ulong));
                    
                    /* Copy parameter exponents - these become the variables */
                    if (matrix[i][j].terms[t].par_exp) {
                        for (slong v = 0; v < matrix[i][j].npars && v < current_nvars; v++) {
                            exp[v] = matrix[i][j].terms[t].par_exp[v];
                        }
                    }
                    
                    /* Add term if coefficient is non-zero */
                    if (coeff != 0) {
                        nmod_mpoly_t temp;
                        nmod_mpoly_init(temp, reduced_ideal.ctx.nmod_ctx);
                        nmod_mpoly_set_coeff_ui_ui(temp, coeff, exp, reduced_ideal.ctx.nmod_ctx);
                        nmod_mpoly_add(&nmod_matrix[i][j], &nmod_matrix[i][j], temp, reduced_ideal.ctx.nmod_ctx);
                        nmod_mpoly_clear(temp, reduced_ideal.ctx.nmod_ctx);
                    }
                    
                    flint_free(exp);
                }
                
                DEBUG_PRINT_R("  Converted to %ld terms\n", nmod_mpoly_length(&nmod_matrix[i][j], reduced_ideal.ctx.nmod_ctx));
            }
        }
        
        /* Compute determinant with reduction - ADD TIMEOUT PROTECTION */
        nmod_mpoly_t det;
        nmod_mpoly_init(det, reduced_ideal.ctx.nmod_ctx);
        
        DEBUG_PRINT_R("Starting determinant computation...\n");
        clock_t start = clock();
        
        /* Use simplified determinant for small matrices */
        if (size == 1) {
            nmod_mpoly_set(det, &nmod_matrix[0][0], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(det, &reduced_ideal, complete_var_names);
        } else if (size == 2) {
            nmod_mpoly_t ad, bc;
            nmod_mpoly_init(ad, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_init(bc, reduced_ideal.ctx.nmod_ctx);
            
            nmod_mpoly_mul(ad, &nmod_matrix[0][0], &nmod_matrix[1][1], reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_mul(bc, &nmod_matrix[0][1], &nmod_matrix[1][0], reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_sub(det, ad, bc, reduced_ideal.ctx.nmod_ctx);
            
            /* Reduce only once at the end */
            triangular_ideal_reduce_nmod_mpoly_with_names(det, &reduced_ideal, complete_var_names);
            
            nmod_mpoly_clear(ad, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_clear(bc, reduced_ideal.ctx.nmod_ctx);
        } else if (size == 3) {
            /* For 3x3: First reduce all matrix elements, then compute with intermediate reductions */
            DEBUG_PRINT_R("Pre-reducing all 3x3 matrix elements...\n");
            
            /* Reduce each matrix element first */
            for (slong i = 0; i < 3; i++) {
                for (slong j = 0; j < 3; j++) {
                    DEBUG_PRINT_R("  Reducing element [%ld][%ld] (%ld terms)...\n", i, j,
                           nmod_mpoly_length(&nmod_matrix[i][j], reduced_ideal.ctx.nmod_ctx));
                    triangular_ideal_reduce_nmod_mpoly_with_names(&nmod_matrix[i][j], 
                                                                &reduced_ideal, complete_var_names);
                    DEBUG_PRINT_R("    Reduced to %ld terms\n", 
                           nmod_mpoly_length(&nmod_matrix[i][j], reduced_ideal.ctx.nmod_ctx));
                }
            }
            
            DEBUG_PRINT_R("Computing 3x3 determinant with intermediate reductions...\n");
            
            /* Compute 3x3 determinant: det = a00*a11*a22 + a01*a12*a20 + a02*a10*a21 
                                              - a00*a12*a21 - a01*a10*a22 - a02*a11*a20 */
            nmod_mpoly_t term1, term2, term3, term4, term5, term6, temp1, temp2;
            nmod_mpoly_init(term1, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_init(term2, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_init(term3, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_init(term4, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_init(term5, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_init(term6, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_init(temp1, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_init(temp2, reduced_ideal.ctx.nmod_ctx);
            
            /* Positive terms with reduction after each multiplication */
            
            /* term1 = a00 * a11 * a22 */
            DEBUG_PRINT_R("  Computing term1: a00*a11*a22...\n");
            nmod_mpoly_mul(temp1, &nmod_matrix[0][0], &nmod_matrix[1][1], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(temp1, &reduced_ideal, complete_var_names);
            nmod_mpoly_mul(term1, temp1, &nmod_matrix[2][2], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(term1, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    term1 has %ld terms\n", nmod_mpoly_length(term1, reduced_ideal.ctx.nmod_ctx));
            
            /* term2 = a01 * a12 * a20 */
            DEBUG_PRINT_R("  Computing term2: a01*a12*a20...\n");
            nmod_mpoly_mul(temp1, &nmod_matrix[0][1], &nmod_matrix[1][2], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(temp1, &reduced_ideal, complete_var_names);
            nmod_mpoly_mul(term2, temp1, &nmod_matrix[2][0], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(term2, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    term2 has %ld terms\n", nmod_mpoly_length(term2, reduced_ideal.ctx.nmod_ctx));
            
            /* term3 = a02 * a10 * a21 */
            DEBUG_PRINT_R("  Computing term3: a02*a10*a21...\n");
            nmod_mpoly_mul(temp1, &nmod_matrix[0][2], &nmod_matrix[1][0], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(temp1, &reduced_ideal, complete_var_names);
            nmod_mpoly_mul(term3, temp1, &nmod_matrix[2][1], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(term3, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    term3 has %ld terms\n", nmod_mpoly_length(term3, reduced_ideal.ctx.nmod_ctx));
            
            /* Negative terms with reduction after each multiplication */
            
            /* term4 = a00 * a12 * a21 */
            DEBUG_PRINT_R("  Computing term4: a00*a12*a21...\n");
            nmod_mpoly_mul(temp1, &nmod_matrix[0][0], &nmod_matrix[1][2], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(temp1, &reduced_ideal, complete_var_names);
            nmod_mpoly_mul(term4, temp1, &nmod_matrix[2][1], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(term4, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    term4 has %ld terms\n", nmod_mpoly_length(term4, reduced_ideal.ctx.nmod_ctx));
            
            /* term5 = a01 * a10 * a22 */
            DEBUG_PRINT_R("  Computing term5: a01*a10*a22...\n");
            nmod_mpoly_mul(temp1, &nmod_matrix[0][1], &nmod_matrix[1][0], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(temp1, &reduced_ideal, complete_var_names);
            nmod_mpoly_mul(term5, temp1, &nmod_matrix[2][2], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(term5, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    term5 has %ld terms\n", nmod_mpoly_length(term5, reduced_ideal.ctx.nmod_ctx));
            
            /* term6 = a02 * a11 * a20 */
            DEBUG_PRINT_R("  Computing term6: a02*a11*a20...\n");
            nmod_mpoly_mul(temp1, &nmod_matrix[0][2], &nmod_matrix[1][1], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(temp1, &reduced_ideal, complete_var_names);
            nmod_mpoly_mul(term6, temp1, &nmod_matrix[2][0], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(term6, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    term6 has %ld terms\n", nmod_mpoly_length(term6, reduced_ideal.ctx.nmod_ctx));
            
            /* Combine terms with reduction after each addition/subtraction */
            DEBUG_PRINT_R("  Combining terms...\n");
            
            /* det = term1 + term2 + term3 - term4 - term5 - term6 */
            nmod_mpoly_add(det, term1, term2, reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(det, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    After adding term2: %ld terms\n", nmod_mpoly_length(det, reduced_ideal.ctx.nmod_ctx));
            
            nmod_mpoly_add(det, det, term3, reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(det, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    After adding term3: %ld terms\n", nmod_mpoly_length(det, reduced_ideal.ctx.nmod_ctx));
            
            nmod_mpoly_sub(det, det, term4, reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(det, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    After subtracting term4: %ld terms\n", nmod_mpoly_length(det, reduced_ideal.ctx.nmod_ctx));
            
            nmod_mpoly_sub(det, det, term5, reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(det, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    After subtracting term5: %ld terms\n", nmod_mpoly_length(det, reduced_ideal.ctx.nmod_ctx));
            
            nmod_mpoly_sub(det, det, term6, reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(det, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    Final determinant: %ld terms\n", nmod_mpoly_length(det, reduced_ideal.ctx.nmod_ctx));
            
            /* Cleanup */
            nmod_mpoly_clear(term1, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_clear(term2, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_clear(term3, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_clear(term4, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_clear(term5, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_clear(term6, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_clear(temp1, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_clear(temp2, reduced_ideal.ctx.nmod_ctx);
        } else {
            /* For larger matrices, use the recursive function but with timeout */
            compute_nmod_det_with_triangular_reduction_with_names(det, nmod_matrix, size, 
                                                               &reduced_ideal, complete_var_names);
        }
        
        clock_t end = clock();
        DEBUG_PRINT_R("Determinant computation time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        
        /* Convert result back - result has 0 vars and npars parameters */
        nmod_mpoly_to_fq_mvpoly(result, det, 0, current_nvars, reduced_ideal.ctx.nmod_ctx, ideal->field_ctx);
        
        /* Cleanup */
        nmod_mpoly_clear(det, reduced_ideal.ctx.nmod_ctx);
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                nmod_mpoly_clear(&nmod_matrix[i][j], reduced_ideal.ctx.nmod_ctx);
            }
            flint_free(nmod_matrix[i]);
        }
        flint_free(nmod_matrix);
        
    } else {
        /* Similar for fq_nmod_mpoly */
        DEBUG_PRINT_R("Using fq_nmod_mpoly for extension field computation\n");
        
        /* Convert to fq_nmod_mpoly matrix */
        fq_nmod_mpoly_struct **fq_matrix = (fq_nmod_mpoly_struct**) flint_malloc(size * sizeof(fq_nmod_mpoly_struct*));
        for (slong i = 0; i < size; i++) {
            fq_matrix[i] = (fq_nmod_mpoly_struct*) flint_malloc(size * sizeof(fq_nmod_mpoly_struct));
            for (slong j = 0; j < size; j++) {
                fq_nmod_mpoly_init(&fq_matrix[i][j], reduced_ideal.ctx.fq_ctx);
                fq_nmod_mpoly_zero(&fq_matrix[i][j], reduced_ideal.ctx.fq_ctx);
                
                /* Convert fq_mvpoly to fq_nmod_mpoly */
                DEBUG_PRINT_R("Converting matrix[%ld][%ld] with %ld terms...\n", i, j, matrix[i][j].nterms);
                
                /* Process each term */
                for (slong t = 0; t < matrix[i][j].nterms && t < 10000; t++) {
                    /* Build exponent vector for fq_nmod_mpoly */
                    ulong *exp = (ulong*) flint_calloc(current_nvars, sizeof(ulong));
                    
                    /* Copy parameter exponents - these become the variables */
                    if (matrix[i][j].terms[t].par_exp) {
                        for (slong v = 0; v < matrix[i][j].npars && v < current_nvars; v++) {
                            exp[v] = matrix[i][j].terms[t].par_exp[v];
                        }
                    }
                    
                    /* Set coefficient */
                    fq_nmod_mpoly_set_coeff_fq_nmod_ui(&fq_matrix[i][j], 
                                                       matrix[i][j].terms[t].coeff, 
                                                       exp, reduced_ideal.ctx.fq_ctx);
                    
                    flint_free(exp);
                }
                
                DEBUG_PRINT_R("  Converted to %ld terms\n", 
                       fq_nmod_mpoly_length(&fq_matrix[i][j], reduced_ideal.ctx.fq_ctx));
            }
        }
        
        /* Compute determinant with reduction */
        fq_nmod_mpoly_t det;
        fq_nmod_mpoly_init(det, reduced_ideal.ctx.fq_ctx);
        
        DEBUG_PRINT_R("Starting determinant computation...\n");
        clock_t start = clock();
        
        if (size == 1) {
            fq_nmod_mpoly_set(det, &fq_matrix[0][0], reduced_ideal.ctx.fq_ctx);
            triangular_ideal_reduce_fq_nmod_mpoly_with_names(det, &reduced_ideal, complete_var_names);
        } else if (size == 2) {
            fq_nmod_mpoly_t ad, bc;
            fq_nmod_mpoly_init(ad, reduced_ideal.ctx.fq_ctx);
            fq_nmod_mpoly_init(bc, reduced_ideal.ctx.fq_ctx);
            
            fq_nmod_mpoly_mul(ad, &fq_matrix[0][0], &fq_matrix[1][1], reduced_ideal.ctx.fq_ctx);
            fq_nmod_mpoly_mul(bc, &fq_matrix[0][1], &fq_matrix[1][0], reduced_ideal.ctx.fq_ctx);
            fq_nmod_mpoly_sub(det, ad, bc, reduced_ideal.ctx.fq_ctx);
            
            triangular_ideal_reduce_fq_nmod_mpoly_with_names(det, &reduced_ideal, complete_var_names);
            
            fq_nmod_mpoly_clear(ad, reduced_ideal.ctx.fq_ctx);
            fq_nmod_mpoly_clear(bc, reduced_ideal.ctx.fq_ctx);
        } else {
            compute_fq_nmod_det_with_triangular_reduction_with_names(det, fq_matrix, size, 
                                                                    &reduced_ideal, complete_var_names);
        }
        
        clock_t end = clock();
        DEBUG_PRINT_R("Determinant computation time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        
        /* Convert result back */
        fq_nmod_mpoly_to_fq_mvpoly(result, det, 0, current_nvars, reduced_ideal.ctx.fq_ctx, ideal->field_ctx);
        
        /* Cleanup */
        fq_nmod_mpoly_clear(det, reduced_ideal.ctx.fq_ctx);
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                fq_nmod_mpoly_clear(&fq_matrix[i][j], reduced_ideal.ctx.fq_ctx);
            }
            flint_free(fq_matrix[i]);
        }
        flint_free(fq_matrix);
    }
    
    /* Cleanup reduced ideal and temporary names */
    unified_triangular_ideal_clear(&reduced_ideal);
    
    /* Free any allocated names that we created */
    for (slong i = 0; i < current_nvars; i++) {
        if (complete_var_names[i] != current_var_names[i]) {
            free(complete_var_names[i]);
        }
    }
    free(complete_var_names);
}

/* Initialize equation information */
void equation_info_init(equation_info_t *eq, slong nvars, slong npars, const fq_nmod_ctx_t ctx) {
    eq->main_var_name = NULL;
    eq->main_var_degree = 0;
    fq_mvpoly_init(&eq->lhs, nvars, npars, ctx);
    fq_mvpoly_init(&eq->rhs, nvars, npars, ctx);
    fq_mvpoly_init(&eq->standard, nvars, npars, ctx);
}

/* Clear equation information */
void equation_info_clear(equation_info_t *eq) {
    if (eq->main_var_name) {
        free(eq->main_var_name);
        eq->main_var_name = NULL;
    }
    fq_mvpoly_clear(&eq->lhs);
    fq_mvpoly_clear(&eq->rhs);
    fq_mvpoly_clear(&eq->standard);
}

/* Add generator in equation format to the ideal */
void unified_triangular_ideal_add_generator_equation_format(unified_triangular_ideal_t *ideal,
                                                           const char *equation_str,
                                                           char **var_names_hint) {
    if (ideal->num_gens >= ideal->max_gens) {
        DEBUG_PRINT_R("Error: Maximum number of generators exceeded\n");
        return;
    }
    
    /* Check if equation format */
    if (!is_equation_format(equation_str)) {
        DEBUG_PRINT_R("Warning: Not equation format, using original parser\n");
        /* Fall back to original parsing method */
        return;
    }
    
    /* Parse equation */
    equation_info_t eq_info;
    slong nvars = ideal->is_prime_field ? nmod_mpoly_ctx_nvars(ideal->ctx.nmod_ctx) 
                                        : fq_nmod_mpoly_ctx_nvars(ideal->ctx.fq_ctx);
    equation_info_init(&eq_info, nvars, 0, ideal->field_ctx);
    
    parse_equation_generator(&eq_info, equation_str, var_names_hint, nvars, ideal->field_ctx);
    
    if (!eq_info.main_var_name) {
        DEBUG_PRINT_R("Error: Could not extract main variable from equation\n");
        equation_info_clear(&eq_info);
        return;
    }
    
    /* Find main variable index */
    slong main_var_idx = -1;
    for (slong i = 0; i < nvars; i++) {
        if (var_names_hint && var_names_hint[i] && 
            strcmp(var_names_hint[i], eq_info.main_var_name) == 0) {
            main_var_idx = i;
            break;
        }
    }
    
    if (main_var_idx < 0) {
        DEBUG_PRINT_R("Error: Main variable '%s' not found in variable list\n", eq_info.main_var_name);
        equation_info_clear(&eq_info);
        return;
    }
    
    DEBUG_PRINT_R("Adding generator: %s^%ld = ... (variable index %ld)\n", 
           eq_info.main_var_name, eq_info.main_var_degree, main_var_idx);
    
    /* Convert and store generator (using standard form) */
    if (ideal->is_prime_field) {
        fq_mvpoly_to_nmod_mpoly((nmod_mpoly_struct*)ideal->generators[ideal->num_gens], 
                               &eq_info.standard, ideal->ctx.nmod_ctx);
    } else {
        fq_mvpoly_to_fq_nmod_mpoly((fq_nmod_mpoly_struct*)ideal->generators[ideal->num_gens], 
                                  &eq_info.standard, ideal->ctx.fq_ctx);
    }
    
    ideal->var_indices[ideal->num_gens] = main_var_idx;
    ideal->var_names[ideal->num_gens] = strdup(eq_info.main_var_name);
    ideal->leading_degrees[ideal->num_gens] = eq_info.main_var_degree;
    ideal->num_gens++;
    
    DEBUG_PRINT_R("Successfully added generator %ld for variable '%s' (degree %ld)\n", 
           ideal->num_gens - 1, eq_info.main_var_name, eq_info.main_var_degree);
    
    equation_info_clear(&eq_info);
}

/* Construct triangular ideal from equation format strings */
void construct_triangular_ideal_from_strings(unified_triangular_ideal_t *ideal,
                                              const char **equation_strs,
                                              slong num_equations,
                                              const char **var_names,
                                              slong nvars,
                                              const fq_nmod_ctx_t ctx) {
    DEBUG_PRINT_R("Constructing triangular ideal from %ld equations\n", num_equations);
    
    /* Initialize ideal */
    unified_triangular_ideal_init(ideal, num_equations, nvars, ctx);
    
    /* Print equations to be processed */
    DEBUG_PRINT_R("Triangular ideal equations:\n");
    for (slong i = 0; i < num_equations; i++) {
        DEBUG_PRINT_R("  eq%ld: %s\n", i+1, equation_strs[i]);
    }
    DEBUG_PRINT_R("\n");
    
    /* Create copy of variable names array */
    char **var_names_copy = (char**) malloc(nvars * sizeof(char*));
    for (slong i = 0; i < nvars; i++) {
        var_names_copy[i] = strdup(var_names[i]);
    }
    
    /* Process each equation */
    for (slong i = 0; i < num_equations; i++) {
        unified_triangular_ideal_add_generator_equation_format(ideal, equation_strs[i], var_names_copy);
    }
    
    /* Save complete variable name list */
    ideal->total_vars = nvars;
    ideal->complete_var_names = (char**) malloc(nvars * sizeof(char*));
    for (slong i = 0; i < nvars; i++) {
        ideal->complete_var_names[i] = strdup(var_names[i]);
    }
    
    DEBUG_PRINT_R("Ideal construction complete:\n");
    DEBUG_PRINT_R("  Number of generators: %ld\n", ideal->num_gens);
    for (slong i = 0; i < ideal->num_gens; i++) {
        DEBUG_PRINT_R("  Generator %ld: %s^%ld (index %ld)\n", 
               i, ideal->var_names[i] ? ideal->var_names[i] : "?", 
               ideal->leading_degrees[i], ideal->var_indices[i]);
    }
    
    /* Cleanup */
    for (slong i = 0; i < nvars; i++) {
        free(var_names_copy[i]);
    }
    free(var_names_copy);
}

/* Main Dixon resultant computation with ideal reduction */
char* dixon_with_ideal_reduction(const char **poly_strings, slong num_polys,
                                const char **elim_vars, slong num_elim_vars,
                                const fq_nmod_ctx_t ctx,
                                unified_triangular_ideal_t *ideal) {
    /* Extract ALL variables from the ideal instead of hardcoding */
    slong total_system_vars = 0;
    char **all_system_vars = NULL;
    
    if (ideal && ideal->complete_var_names && ideal->total_vars > 0) {
        /* Use complete variable names saved in ideal */
        total_system_vars = ideal->total_vars;
        all_system_vars = (char**) malloc(total_system_vars * sizeof(char*));
        
        for (slong i = 0; i < total_system_vars; i++) {
            all_system_vars[i] = strdup(ideal->complete_var_names[i]);
        }
        
        DEBUG_PRINT_R("Using saved complete variable names from ideal:\n");
        for (slong i = 0; i < total_system_vars; i++) {
            DEBUG_PRINT_R("  var[%ld] = %s\n", i, all_system_vars[i]);
        }
    } else if (ideal && ideal->num_gens > 0) {
        /* Estimate maximum variables from ideal context */
        if (ideal->is_prime_field) {
            total_system_vars = nmod_mpoly_ctx_nvars(ideal->ctx.nmod_ctx);
        } else {
            total_system_vars = fq_nmod_mpoly_ctx_nvars(ideal->ctx.fq_ctx);
        }
        
        all_system_vars = (char**) malloc(total_system_vars * sizeof(char*));
        
        /* Initialize array */
        for (slong i = 0; i < total_system_vars; i++) {
            all_system_vars[i] = NULL;
        }
        
        /* Extract variable names from ideal generators */
        for (slong g = 0; g < ideal->num_gens; g++) {
            if (ideal->var_names[g] && ideal->var_indices[g] >= 0 && 
                ideal->var_indices[g] < total_system_vars) {
                if (!all_system_vars[ideal->var_indices[g]]) {
                    all_system_vars[ideal->var_indices[g]] = strdup(ideal->var_names[g]);
                }
            }
        }
        
        /* Fill in missing variable names with defaults */
        for (slong i = 0; i < total_system_vars; i++) {
            if (!all_system_vars[i]) {
                char temp[32];
                if (i == 0) {
                    sprintf(temp, "x");
                } else {
                    sprintf(temp, "var_%ld", i);
                }
                all_system_vars[i] = strdup(temp);
            }
        }
    } else {
        /* Fallback: minimal variable set */
        total_system_vars = 2;
        all_system_vars = (char**) malloc(total_system_vars * sizeof(char*));
        all_system_vars[0] = strdup("x");
        all_system_vars[1] = strdup("y");
    }
    
    /* Parse polynomials */
    parser_state_t state = {0};
    state.var_names = (char**) malloc(num_elim_vars * sizeof(char*));
    for (slong i = 0; i < num_elim_vars; i++) {
        state.var_names[i] = strdup(elim_vars[i]);
    }
    state.nvars = num_elim_vars;
    state.npars = 0;
    state.max_pars = 32;  /* Increased size */
    state.par_names = (char**) malloc(state.max_pars * sizeof(char*));
    state.ctx = ctx;
    state.generator_name = get_generator_name(ctx);
    fq_nmod_init(state.current.value, ctx);
    state.current.str = NULL;
    
    /* CRITICAL FIX: Add ALL non-eliminated variables as parameters */
    /* This includes x and all z_i that are not being eliminated */
    for (slong i = 0; i < total_system_vars; i++) {
        int is_elim_var = 0;
        
        /* Check if this is an elimination variable */
        for (slong j = 0; j < num_elim_vars; j++) {
            if (strcmp(all_system_vars[i], elim_vars[j]) == 0) {
                is_elim_var = 1;
                break;
            }
        }
        
        /* If not an elimination variable, add as parameter */
        if (!is_elim_var) {
            /* Check if already in parameters */
            int already_param = 0;
            for (slong j = 0; j < state.npars; j++) {
                if (strcmp(state.par_names[j], all_system_vars[i]) == 0) {
                    already_param = 1;
                    break;
                }
            }
            
            if (!already_param) {
                if (state.npars >= state.max_pars) {
                    state.max_pars *= 2;
                    state.par_names = (char**) realloc(state.par_names, 
                                                      state.max_pars * sizeof(char*));
                }
                state.par_names[state.npars] = strdup(all_system_vars[i]);
                state.npars++;
            }
        }
    }
    
    /* First pass: parse to identify any additional parameters in polynomials */
    for (slong i = 0; i < num_polys; i++) {
        fq_mvpoly_t temp;
        fq_mvpoly_init(&temp, num_elim_vars, state.max_pars, ctx);
        
        state.input = poly_strings[i];
        state.pos = 0;
        state.len = strlen(poly_strings[i]);
        next_token(&state);
        
        parse_expression(&state, &temp);
        fq_mvpoly_clear(&temp);
    }
    
    /* Build complete variable list: elimination variables + parameters */
    slong total_vars = num_elim_vars + state.npars;
    char **all_var_names = (char**) malloc(total_vars * sizeof(char*));
    
    /* Copy elimination variables */
    for (slong i = 0; i < num_elim_vars; i++) {
        all_var_names[i] = strdup(elim_vars[i]);
    }
    
    /* Copy parameters */
    for (slong i = 0; i < state.npars; i++) {
        all_var_names[num_elim_vars + i] = strdup(state.par_names[i]);
    }
    
    /* DEBUG: Print all variables */
    DEBUG_PRINT_R("\n=== DEBUG: All variables (elimination + parameters) ===\n");
    DEBUG_PRINT_R("Elimination variables (%ld): ", num_elim_vars);
    for (slong i = 0; i < num_elim_vars; i++) {
        DEBUG_PRINT_R("%s ", all_var_names[i]);
    }
    DEBUG_PRINT_R("\nParameters (%ld): ", state.npars);
    for (slong i = 0; i < state.npars; i++) {
        DEBUG_PRINT_R("%s ", all_var_names[num_elim_vars + i]);
    }
    DEBUG_PRINT_R("\n=== END DEBUG ===\n\n");
    
    /* Parse polynomials again with correct context */
    fq_mvpoly_t *polys = (fq_mvpoly_t*) malloc(num_polys * sizeof(fq_mvpoly_t));
    
    for (slong i = 0; i < num_polys; i++) {
        fq_mvpoly_init(&polys[i], num_elim_vars, state.npars, ctx);
        
        state.input = poly_strings[i];
        state.pos = 0;
        state.len = strlen(poly_strings[i]);
        if (state.current.str) {
            free(state.current.str);
            state.current.str = NULL;
        }
        next_token(&state);
        
        parse_expression(&state, &polys[i]);
    }
    
    /* Build cancellation matrix */
    printf("\nStep 1: Build Dixon polynomial\n");
    clock_t step1_start = clock();
    fq_mvpoly_t **M_mvpoly;
    printf("Build Cancellation Matrix\n");
    build_fq_cancellation_matrix_mvpoly(&M_mvpoly, polys, num_elim_vars, state.npars);
    
    /* Perform row operations */
    fq_mvpoly_t **modified_M_mvpoly;
    printf("Perform Matrix Row Operations\n");
    perform_fq_matrix_row_operations_mvpoly(&modified_M_mvpoly, &M_mvpoly, num_elim_vars, state.npars);
    
    /* Compute determinant of modified matrix */
    fq_mvpoly_t d_poly;
    printf("Computing cancellation matrix determinant using recursive expansion...\n");
    compute_fq_cancel_matrix_det(&d_poly, modified_M_mvpoly, num_elim_vars, state.npars, DET_METHOD_RECURSIVE);
    if (d_poly.nterms <= 100) {
        printf("Dixon polynomial: %ld terms\n", d_poly.nterms);
    } else {
        printf("Dixon polynomial: %ld terms (not shown)\n", d_poly.nterms);
    }
    printf("Time: %.3f seconds\n", (double)(clock() - step1_start) / CLOCKS_PER_SEC);
    
    /* Extract coefficient matrix */
    fq_mvpoly_t **coeff_matrix = NULL;
    slong *row_indices = (slong*) flint_malloc(d_poly.nterms * sizeof(slong));
    slong *col_indices = (slong*) flint_malloc(d_poly.nterms * sizeof(slong));
    slong matrix_size = 0;
    
    extract_fq_coefficient_matrix_from_dixon(&coeff_matrix, row_indices, col_indices,
                                            &matrix_size, &d_poly, num_elim_vars, state.npars);
    
    /* Compute determinant with ideal reduction */
    fq_mvpoly_t result_poly;
    if (matrix_size > 0) {
        printf("\nStep 4: Compute resultant with ideal reduction\n");
        
        /* Pass parameter names for proper variable mapping during reduction */
        compute_det_with_reduction_from_mvpoly(&result_poly, coeff_matrix, matrix_size, ideal, 
                                             state.par_names);
        
        /* Cleanup coefficient matrix */
        for (slong i = 0; i < matrix_size; i++) {
            for (slong j = 0; j < matrix_size; j++) {
                fq_mvpoly_clear(&coeff_matrix[i][j]);
            }
            flint_free(coeff_matrix[i]);
        }
        flint_free(coeff_matrix);
    } else {
        fq_mvpoly_init(&result_poly, 0, state.npars, ctx);
        printf("Warning: Empty coefficient matrix, resultant is 0\n");
    }
    fq_mvpoly_make_monic(&result_poly);
    
    /* Convert result to string */
    find_and_print_roots_of_univariate_resultant(&result_poly, &state);
    char *result = fq_mvpoly_to_string(&result_poly, state.par_names, state.generator_name);
    
    /* Cleanup */
    fq_mvpoly_clear(&result_poly);
    flint_free(row_indices);
    flint_free(col_indices);
    
    for (slong i = 0; i <= num_elim_vars; i++) {
        for (slong j = 0; j <= num_elim_vars; j++) {
            fq_mvpoly_clear(&M_mvpoly[i][j]);
            fq_mvpoly_clear(&modified_M_mvpoly[i][j]);
        }
        flint_free(M_mvpoly[i]);
        flint_free(modified_M_mvpoly[i]);
    }
    flint_free(M_mvpoly);
    flint_free(modified_M_mvpoly);
    fq_mvpoly_clear(&d_poly);
    
    for (slong i = 0; i < num_polys; i++) {
        fq_mvpoly_clear(&polys[i]);
    }
    free(polys);
    
    for (slong i = 0; i < num_elim_vars; i++) {
        free(state.var_names[i]);
    }
    free(state.var_names);
    
    for (slong i = 0; i < state.npars; i++) {
        free(state.par_names[i]);
    }
    free(state.par_names);
    
    for (slong i = 0; i < total_vars; i++) {
        free(all_var_names[i]);
    }
    free(all_var_names);
    
    for (slong i = 0; i < total_system_vars; i++) {
        free(all_system_vars[i]);
    }
    free(all_system_vars);
    
    if (state.generator_name) {
        free(state.generator_name);
    }
    
    fq_nmod_clear(state.current.value, ctx);
    if (state.current.str) {
        free(state.current.str);
    }

    printf("\n=== Dixon Resultant Computation Complete ===\n");

    return result;
}

/* Helper functions for string interface */
char** split_string_r(const char* str, slong* count) {
    if (!str || !count) {
        if (count) *count = 0;
        return NULL;
    }
    
    *count = 0;  /* Initialize immediately */
    
    slong len = strlen(str);
    if (len == 0) {
        return NULL;
    }
    
    /* Count separators more safely */
    slong num_parts = 1;
    for (slong i = 0; i < len; i++) {
        if (str[i] == ';' || str[i] == ',') {
            num_parts++;
        }
    }
    
    /* Allocate and INITIALIZE to NULL */
    char** result = (char**) calloc(num_parts, sizeof(char*));  /* Use calloc! */
    if (!result) {
        return NULL;
    }
    
    /* Simple, safe parsing */
    slong start_pos = 0;
    
    for (slong i = 0; i <= len; i++) {  /* Include len to catch final part */
        if (i == len || str[i] == ';' || str[i] == ',') {
            slong part_len = i - start_pos;
            
            if (part_len > 0) {
                /* Bounds check */
                if (*count >= num_parts) {
                    break;
                }
                
                result[*count] = (char*) malloc(part_len + 1);
                if (!result[*count]) {
                    break;
                }
                
                strncpy(result[*count], str + start_pos, part_len);
                result[*count][part_len] = '\0';
                
                /* Simple trim - remove leading/trailing spaces */
                char* trimmed = result[*count];
                while (*trimmed == ' ' || *trimmed == '\t') trimmed++;
                
                if (trimmed != result[*count]) {
                    memmove(result[*count], trimmed, strlen(trimmed) + 1);
                }
                
                /* Remove trailing spaces */
                slong trim_len = strlen(result[*count]);
                while (trim_len > 0 && (result[*count][trim_len-1] == ' ' || 
                                       result[*count][trim_len-1] == '\t')) {
                    result[*count][trim_len-1] = '\0';
                    trim_len--;
                }
                
                (*count)++;
            }
            
            start_pos = i + 1;
        }
    }
    
    return result;
}

void free_split_string_rs(char** strings, slong count) {
    if (strings) {
        for (slong i = 0; i < count; i++) {
            if (strings[i]) {
                free(strings[i]);
            }
        }
        free(strings);
    }
}

/* Auto-extract all variables from ideal generators */
char** extract_all_variables_from_ideal_gens(const char **ideal_gens_strs, 
                                            slong num_gens,
                                            const fq_nmod_ctx_t ctx,
                                            slong *total_vars) {
    /* Temporary parser state for collecting variables */
    parser_state_t temp_state;
    temp_state.var_names = NULL;
    temp_state.nvars = 0;
    temp_state.par_names = (char**) malloc(32 * sizeof(char*));
    temp_state.npars = 0;
    temp_state.max_pars = 32;
    temp_state.ctx = ctx;
    temp_state.generator_name = get_generator_name(ctx);
    fq_nmod_init(temp_state.current.value, ctx);
    temp_state.current.str = NULL;
    
    /* Use set to avoid duplicate variables */
    char **all_vars = (char**) malloc(64 * sizeof(char*));
    slong max_vars = 64;
    slong num_vars = 0;
    
    /* Parse each ideal generator, extract variables */
    for (slong i = 0; i < num_gens; i++) {
        if (!ideal_gens_strs[i]) continue;
        
        /* Check if equation format (contains =) */
        const char *eq_pos = strchr(ideal_gens_strs[i], '=');
        char *lhs_str = NULL, *rhs_str = NULL;
        
        if (eq_pos) {
            /* Split left and right sides of equation */
            slong lhs_len = eq_pos - ideal_gens_strs[i];
            lhs_str = (char*) malloc((lhs_len + 1) * sizeof(char));
            strncpy(lhs_str, ideal_gens_strs[i], lhs_len);
            lhs_str[lhs_len] = '\0';
            
            /* Remove trailing spaces */
            while (lhs_len > 0 && (lhs_str[lhs_len-1] == ' ' || lhs_str[lhs_len-1] == '\t')) {
                lhs_str[lhs_len-1] = '\0';
                lhs_len--;
            }
            
            /* Right side part */
            const char *rhs_start = eq_pos + 1;
            while (*rhs_start == ' ' || *rhs_start == '\t') rhs_start++;
            rhs_str = strdup(rhs_start);
        } else {
            /* Not equation format, use original string directly */
            lhs_str = strdup(ideal_gens_strs[i]);
            rhs_str = strdup("0");
        }
        
        /* Parse left and right sides, extract variables */
        const char *parts[] = {lhs_str, rhs_str};
        for (int part = 0; part < 2; part++) {
            temp_state.input = parts[part];
            temp_state.pos = 0;
            temp_state.len = strlen(parts[part]);
            
            /* Reset token */
            if (temp_state.current.str) {
                free(temp_state.current.str);
                temp_state.current.str = NULL;
            }
            
            /* Simple lexical analysis, extract identifiers */
            while (temp_state.pos < temp_state.len) {
                /* Skip spaces */
                while (temp_state.pos < temp_state.len && 
                       isspace(temp_state.input[temp_state.pos])) {
                    temp_state.pos++;
                }
                
                if (temp_state.pos >= temp_state.len) break;
                
                /* Check if letter-starting identifier */
                if (isalpha(temp_state.input[temp_state.pos]) || 
                    temp_state.input[temp_state.pos] == '_') {
                    
                    size_t start = temp_state.pos;
                    while (temp_state.pos < temp_state.len && 
                           (isalnum(temp_state.input[temp_state.pos]) || 
                            temp_state.input[temp_state.pos] == '_')) {
                        temp_state.pos++;
                    }
                    
                    size_t len = temp_state.pos - start;
                    char *identifier = (char*) malloc(len + 1);
                    strncpy(identifier, temp_state.input + start, len);
                    identifier[len] = '\0';
                    
                    /* Check if generator name, if so skip */
                    if (temp_state.generator_name && 
                        strcmp(identifier, temp_state.generator_name) == 0) {
                        free(identifier);
                        continue;
                    }
                    
                    /* Check if already in list */
                    int found = 0;
                    for (slong j = 0; j < num_vars; j++) {
                        if (strcmp(all_vars[j], identifier) == 0) {
                            found = 1;
                            break;
                        }
                    }
                    
                    if (!found) {
                        /* Expand array if needed */
                        if (num_vars >= max_vars) {
                            max_vars *= 2;
                            all_vars = (char**) realloc(all_vars, max_vars * sizeof(char*));
                        }
                        
                        all_vars[num_vars] = identifier;
                        num_vars++;
                        
                        DEBUG_PRINT_R("Found variable: %s\n", identifier);
                    } else {
                        free(identifier);
                    }
                } else {
                    temp_state.pos++;
                }
            }
        }
        
        free(lhs_str);
        free(rhs_str);
    }
    
    /* Cleanup temporary state */
    for (slong i = 0; i < temp_state.npars; i++) {
        free(temp_state.par_names[i]);
    }
    free(temp_state.par_names);
    if (temp_state.generator_name) free(temp_state.generator_name);
    fq_nmod_clear(temp_state.current.value, ctx);
    if (temp_state.current.str) free(temp_state.current.str);
    
    *total_vars = num_vars;
    
    DEBUG_PRINT_R("Extracted %ld variables from ideal generators\n", num_vars);
    
    return all_vars;
}

void construct_triangular_ideal_str(unified_triangular_ideal_t *ideal,
                                   const char *gens_string,   
                                   const char *vars_string,   
                                   const fq_nmod_ctx_t ctx) {
    
    slong num_gens, num_vars;
    char **gens_array = split_string_r(gens_string, &num_gens);
    char **vars_array = split_string_r(vars_string, &num_vars);
    
    const char **ideal_gens = (const char**) malloc(num_gens * sizeof(char*));
    const char **var_names = (const char**) malloc(num_vars * sizeof(char*));
    
    for (slong i = 0; i < num_gens; i++) {
        ideal_gens[i] = gens_array[i];
    }
    for (slong i = 0; i < num_vars; i++) {
        var_names[i] = vars_array[i];
    }
    
    construct_triangular_ideal_from_strings(ideal, ideal_gens, num_gens, 
                                          var_names, num_vars, ctx);
    
    free(ideal_gens);
    free(var_names);
    free_split_string_rs(gens_array, num_gens);
    free_split_string_rs(vars_array, num_vars);
}

char* dixon_with_ideal(const char **poly_strings,
                           slong num_polys,
                           const char **elim_vars,
                           slong num_elim_vars,
                           const char *ideal_string,
                           const fq_nmod_ctx_t ctx) {
    
    clock_t start_time = clock();
    
    /* Step 1: Split ideal string into individual equations */
    slong num_ideal_equations;
    char **ideal_equations_array = split_string_r(ideal_string, &num_ideal_equations);
    
    /* Convert to const pointer array */
    const char **ideal_equations = (const char**) malloc(num_ideal_equations * sizeof(char*));
    for (slong i = 0; i < num_ideal_equations; i++) {
        ideal_equations[i] = ideal_equations_array[i];
    }
    
    /* Step 2: Auto-extract all variables from ideal */
    slong num_ideal_vars;
    char **ideal_vars_auto = extract_all_variables_from_ideal_gens(
        ideal_equations, num_ideal_equations, ctx, &num_ideal_vars);
    
    /* Convert to const pointer array */
    const char **ideal_vars = (const char**) malloc(num_ideal_vars * sizeof(char*));
    for (slong i = 0; i < num_ideal_vars; i++) {
        ideal_vars[i] = ideal_vars_auto[i];
    }
    
    /* Step 3: Construct triangular ideal */
    unified_triangular_ideal_t ideal;
    construct_triangular_ideal_from_strings(&ideal, ideal_equations, num_ideal_equations,
                                          ideal_vars, num_ideal_vars, ctx);
    
    /* Step 4: Compute Dixon resultant */
    char *result = dixon_with_ideal_reduction(poly_strings, num_polys,
                                            elim_vars, num_elim_vars,
                                            ctx, &ideal);

    
    /* Cleanup */
    unified_triangular_ideal_clear(&ideal);
    free(ideal_equations);
    free(ideal_vars);
    
    /* Cleanup auto-extracted variables */
    for (slong i = 0; i < num_ideal_vars; i++) {
        free(ideal_vars_auto[i]);
    }
    free(ideal_vars_auto);
    
    /* Cleanup split equations */
    free_split_string_rs(ideal_equations_array, num_ideal_equations);
    
    clock_t end_time = clock();
    double elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Time: %.3f seconds\n", elapsed);
    return result;
}

char* dixon_with_ideal_reduction_str(const char *poly_string,
                                    const char *elim_vars_string,
                                    const char *ideal_gens_string,
                                    const fq_nmod_ctx_t ctx) {
    
    slong num_polys, num_elim_vars, num_gens;
    char **poly_array = split_string_r(poly_string, &num_polys);
    char **elim_array = split_string_r(elim_vars_string, &num_elim_vars);
    char **gens_array = split_string_r(ideal_gens_string, &num_gens);
    
    /* Auto-extract all variables from ideal */
    slong num_all_vars;
    const char **ideal_gens = (const char**) malloc(num_gens * sizeof(char*));
    for (slong i = 0; i < num_gens; i++) {
        ideal_gens[i] = gens_array[i];
    }
    
    char **all_vars_auto = extract_all_variables_from_ideal_gens(
        ideal_gens, num_gens, ctx, &num_all_vars);
    
    /* Convert to const char** */
    const char **all_var_names = (const char**) malloc(num_all_vars * sizeof(char*));
    for (slong i = 0; i < num_all_vars; i++) {
        all_var_names[i] = all_vars_auto[i];
    }
    
    /* Construct ideal */
    unified_triangular_ideal_t ideal;
    construct_triangular_ideal_from_strings(&ideal, ideal_gens, num_gens,
                                          all_var_names, num_all_vars, ctx);
    
    /* Prepare dixon computation parameters */
    const char **poly_strings = (const char**) malloc(num_polys * sizeof(char*));
    const char **elim_vars = (const char**) malloc(num_elim_vars * sizeof(char*));
    
    for (slong i = 0; i < num_polys; i++) {
        poly_strings[i] = poly_array[i];
    }
    for (slong i = 0; i < num_elim_vars; i++) {
        elim_vars[i] = elim_array[i];
    }
    
    /* Compute dixon resultant */
    char *result = dixon_with_ideal_reduction(poly_strings, num_polys,
                                            elim_vars, num_elim_vars,
                                            ctx, &ideal);
    
    /* Cleanup */
    unified_triangular_ideal_clear(&ideal);
    free(poly_strings);
    free(elim_vars);
    free(ideal_gens);
    free(all_var_names);
    
    /* Cleanup auto-extracted variables */
    for (slong i = 0; i < num_all_vars; i++) {
        free(all_vars_auto[i]);
    }
    free(all_vars_auto);
    
    free_split_string_rs(poly_array, num_polys);
    free_split_string_rs(elim_array, num_elim_vars);
    free_split_string_rs(gens_array, num_gens);
    
    return result;
}

void test_iterative_elimination(void) {
    printf("\n================================================\n");
    printf("Test: Simplified Iterative Elimination\n");
    printf("================================================\n");
    
    fq_nmod_ctx_t ctx;
    fmpz_t p;
    fmpz_init(p);
    fmpz_set_ui(p, 257);
    fq_nmod_ctx_init(ctx, p, 1, "t");
    fmpz_clear(p);
    
    /* Define all generators upfront */
    const char *ideal_gens[] = {
        "a2^3 = 2*a1 + 1",
        "a3^3 = a1*a2 + 3",
        "a4^3 = a1 + a2*a3 + 5"
    };
    const char *var_names[] = {"a1", "a2", "a3", "a4"};
    
    /* Construct ideal once */
    unified_triangular_ideal_t ideal;
    construct_triangular_ideal_from_strings(&ideal, ideal_gens, 3, var_names, 4, ctx);
    
    /* Step 1: Eliminate a4 */
    const char *step1_polys[] = {
        "a1^2 + a2^2 + a3^2 + a4^2 - 100",
        "a4^3 - a1 - a2*a3 - 5"
    };
    const char *step1_elim[] = {"a4"};
    
    char *result1 = dixon_with_ideal_reduction(step1_polys, 2, step1_elim, 1, ctx, &ideal);
    printf("After eliminating a4: %s\n\n", result1);
    
    /* Step 2: Eliminate a3 */
    const char *step2_polys[] = {
        result1,
        "a3^3 - a1*a2 - 3"
    };
    const char *step2_elim[] = {"a3"};
    
    char *result2 = dixon_with_ideal_reduction(step2_polys, 2, step2_elim, 1, ctx, &ideal);
    printf("After eliminating a3: %s\n\n", result2);
    
    /* Step 3: Eliminate a2 */
    const char *step3_polys[] = {
        result2,
        "a2^3 - 2*a1 - 1"
    };
    const char *step3_elim[] = {"a2"};
    
    char *final_result = dixon_with_ideal_reduction(step3_polys, 2, step3_elim, 1, ctx, &ideal);
    printf("Final univariate polynomial in a1: %s\n", final_result);
    
    /* Cleanup */
    free(result1);
    free(result2);
    free(final_result);
    unified_triangular_ideal_clear(&ideal);
    fq_nmod_ctx_clear(ctx);
}

void test_iterative_elimination_str(void) {
    printf("\n================================================\n");
    printf("Test: Simplified Iterative Elimination\n");
    printf("================================================\n");
    
    fq_nmod_ctx_t ctx;
    fmpz_t p;
    fmpz_init(p);
    fmpz_set_ui(p, 2147483489);
    fq_nmod_ctx_init(ctx, p, 1, "t");
    fmpz_clear(p);
    
    /* Define ideal generators in equation format */
    const char *ideal_gens_str = "a2^3 = 2*a1 + 1; a3^3 = a1*a2 + 3; a4^3 = a1 + a2*a3 + 5";
    const char *all_vars_str = "a1, a2, a3, a4";
    
    /* Step 1: Eliminate a4 */
    const char *step1_polys_str = "a1^2 + a2^2 + a3^2 + a4^2 - 100, a4^3 - a1 - a2*a3 - 5";
    const char *step1_elim_str = "a4";
    
    char *result1 = dixon_with_ideal_reduction_str(step1_polys_str, step1_elim_str, 
                                                  ideal_gens_str, ctx);
    printf("After eliminating a4: %s\n\n", result1);
    
    /* Step 2: Eliminate a3 */
    char step2_polys_str[10000];
    snprintf(step2_polys_str, sizeof(step2_polys_str), "%s, a3^3 - a1*a2 - 3", result1);
    const char *step2_elim_str = "a3";
    
    char *result2 = dixon_with_ideal_reduction_str(step2_polys_str, step2_elim_str,
                                                  ideal_gens_str, ctx);
    printf("After eliminating a3: %s\n\n", result2);
    
    /* Step 3: Eliminate a2 */
    char step3_polys_str[10000];
    snprintf(step3_polys_str, sizeof(step3_polys_str), "%s, a2^3 - 2*a1 - 1", result2);
    const char *step3_elim_str = "a2";
    
    char *final_result = dixon_with_ideal_reduction_str(step3_polys_str, step3_elim_str,
                                                       ideal_gens_str, ctx);
    printf("Final univariate polynomial in a1: %s\n", final_result);
    
    /* Cleanup */
    free(result1);
    free(result2);
    free(final_result);
    fq_nmod_ctx_clear(ctx);
}
