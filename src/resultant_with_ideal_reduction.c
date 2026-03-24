#include "resultant_with_ideal_reduction.h"

// Bivariate resultant with ideal reduction - completely rewritten
char* resultant_with_ideal_reduction(const char *poly1_string, const char *poly2_string,
                                   const char *elim_var, 
                                   const fq_nmod_ctx_t ctx,
                                   unified_triangular_ideal_t *ideal) {
    printf("\n=== Resultant with Ideal Reduction ===\n");
    printf("Eliminating variable: %s\n", elim_var);
    
    DEBUG_PRINT_R("Input polynomials:\n");
    DEBUG_PRINT_R("  p1: %s\n", poly1_string);
    DEBUG_PRINT_R("  p2: %s\n", poly2_string);
    
    // Extract complete variable list from ideal or construct default
    char **all_system_vars = NULL;
    slong total_system_vars = 0;
    
    if (ideal && ideal->complete_var_names && ideal->total_vars > 0) {
        // Use saved complete variable names from ideal
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
        // Extract variables from ideal context and generators
        if (ideal->is_prime_field) {
            total_system_vars = nmod_mpoly_ctx_nvars(ideal->ctx.nmod_ctx);
        } else {
            total_system_vars = fq_nmod_mpoly_ctx_nvars(ideal->ctx.fq_ctx);
        }
        
        all_system_vars = (char**) malloc(total_system_vars * sizeof(char*));
        
        // Initialize to NULL
        for (slong i = 0; i < total_system_vars; i++) {
            all_system_vars[i] = NULL;
        }
        
        // Extract variable names from ideal generators
        for (slong g = 0; g < ideal->num_gens; g++) {
            if (ideal->var_names[g] && ideal->var_indices[g] >= 0 && 
                ideal->var_indices[g] < total_system_vars) {
                if (!all_system_vars[ideal->var_indices[g]]) {
                    all_system_vars[ideal->var_indices[g]] = strdup(ideal->var_names[g]);
                }
            }
        }
        
        // Fill missing variables with defaults
        for (slong i = 0; i < total_system_vars; i++) {
            if (!all_system_vars[i]) {
                char temp[32];
                sprintf(temp, "var_%ld", i);
                all_system_vars[i] = strdup(temp);
            }
        }
        
        DEBUG_PRINT_R("Extracted variables from ideal context:\n");
        for (slong i = 0; i < total_system_vars; i++) {
            DEBUG_PRINT_R("  var[%ld] = %s\n", i, all_system_vars[i]);
        }
    } else {
        // Create minimal default variable set
        total_system_vars = 0;
        all_system_vars = NULL;
        
        DEBUG_PRINT_R("No ideal provided, will extract variables from polynomials\n");
    }
    
    // Initialize parser state for variable discovery
    parser_state_t discovery_state;
    discovery_state.var_names = (char**) malloc(1 * sizeof(char*));
    discovery_state.var_names[0] = strdup(elim_var);
    discovery_state.nvars = 1;
    discovery_state.npars = 0;
    discovery_state.max_pars = 32;
    discovery_state.par_names = (char**) malloc(discovery_state.max_pars * sizeof(char*));
    discovery_state.ctx = ctx;
    discovery_state.generator_name = get_generator_name(ctx);
    fq_nmod_init(discovery_state.current.value, ctx);
    discovery_state.current.str = NULL;
    
    // First pass: parse polynomials to discover all variables
    DEBUG_PRINT_R("First pass: discovering variables from polynomials\n");
    
    fq_mvpoly_t temp1, temp2;
    fq_mvpoly_init(&temp1, 1, discovery_state.max_pars, ctx);
    fq_mvpoly_init(&temp2, 1, discovery_state.max_pars, ctx);
    
    // Parse first polynomial
    discovery_state.input = poly1_string;
    discovery_state.pos = 0;
    discovery_state.len = strlen(poly1_string);
    next_token(&discovery_state);
    parse_expression(&discovery_state, &temp1);
    
    // Parse second polynomial
    discovery_state.input = poly2_string;
    discovery_state.pos = 0;
    discovery_state.len = strlen(poly2_string);
    if (discovery_state.current.str) {
        free(discovery_state.current.str);
        discovery_state.current.str = NULL;
    }
    next_token(&discovery_state);
    parse_expression(&discovery_state, &temp2);
    
    fq_mvpoly_clear(&temp1);
    fq_mvpoly_clear(&temp2);
    
    DEBUG_PRINT_R("Discovered parameters from polynomials:\n");
    for (slong i = 0; i < discovery_state.npars; i++) {
        DEBUG_PRINT_R("  param[%ld] = %s\n", i, discovery_state.par_names[i]);
    }
    
    // Merge discovered variables with ideal variables if needed
    if (total_system_vars == 0) {
        // No ideal variables, use discovered parameters as system variables
        total_system_vars = discovery_state.npars + 1; // +1 for elimination variable
        all_system_vars = (char**) malloc(total_system_vars * sizeof(char*));
        
        // Add elimination variable first
        all_system_vars[0] = strdup(elim_var);
        
        // Add discovered parameters
        for (slong i = 0; i < discovery_state.npars; i++) {
            all_system_vars[1 + i] = strdup(discovery_state.par_names[i]);
        }
        
        DEBUG_PRINT_R("Created system variables from discovered parameters:\n");
        for (slong i = 0; i < total_system_vars; i++) {
            DEBUG_PRINT_R("  sys_var[%ld] = %s\n", i, all_system_vars[i]);
        }
    } else {
        // Check if all discovered variables are in the system variables
        for (slong i = 0; i < discovery_state.npars; i++) {
            int found = 0;
            for (slong j = 0; j < total_system_vars; j++) {
                if (strcmp(discovery_state.par_names[i], all_system_vars[j]) == 0) {
                    found = 1;
                    break;
                }
            }
            if (!found) {
                // Add missing variable
                total_system_vars++;
                all_system_vars = (char**) realloc(all_system_vars, total_system_vars * sizeof(char*));
                all_system_vars[total_system_vars - 1] = strdup(discovery_state.par_names[i]);
                
                DEBUG_PRINT_R("Added missing variable to system: %s\n", discovery_state.par_names[i]);
            }
        }
        
        // Ensure elimination variable is in the system
        int elim_found = 0;
        for (slong i = 0; i < total_system_vars; i++) {
            if (strcmp(elim_var, all_system_vars[i]) == 0) {
                elim_found = 1;
                break;
            }
        }
        if (!elim_found) {
            total_system_vars++;
            all_system_vars = (char**) realloc(all_system_vars, total_system_vars * sizeof(char*));
            all_system_vars[total_system_vars - 1] = strdup(elim_var);
        }
    }
    
    // Now determine which variables are parameters (non-elimination variables)
    slong num_params = 0;
    char **param_names = (char**) malloc(total_system_vars * sizeof(char*));
    
    for (slong i = 0; i < total_system_vars; i++) {
        if (strcmp(all_system_vars[i], elim_var) != 0) {
            param_names[num_params] = strdup(all_system_vars[i]);
            num_params++;
        }
    }
    
    DEBUG_PRINT_R("Final parameter list (%ld parameters):\n", num_params);
    for (slong i = 0; i < num_params; i++) {
        DEBUG_PRINT_R("  param[%ld] = %s\n", i, param_names[i]);
    }
    
    // Clean up discovery state
    for (slong i = 0; i < discovery_state.nvars; i++) {
        free(discovery_state.var_names[i]);
    }
    free(discovery_state.var_names);
    for (slong i = 0; i < discovery_state.npars; i++) {
        free(discovery_state.par_names[i]);
    }
    free(discovery_state.par_names);
    if (discovery_state.generator_name) free(discovery_state.generator_name);
    fq_nmod_clear(discovery_state.current.value, ctx);
    if (discovery_state.current.str) free(discovery_state.current.str);
    
    // Initialize final parser state with correct variable and parameter setup
    parser_state_t state = {0};
    state.var_names = (char**) malloc(1 * sizeof(char*));
    state.var_names[0] = strdup(elim_var);
    state.nvars = 1;
    state.npars = num_params;
    state.max_pars = num_params;
    state.par_names = param_names; // Use the prepared parameter names
    state.ctx = ctx;
    state.generator_name = get_generator_name(ctx);
    fq_nmod_init(state.current.value, ctx);
    state.current.str = NULL;
    
    // Parse polynomials with correct context
    fq_mvpoly_t poly1, poly2;
    fq_mvpoly_init(&poly1, 1, state.npars, ctx);
    fq_mvpoly_init(&poly2, 1, state.npars, ctx);
    
    DEBUG_PRINT_R("Second pass: parsing polynomials with correct context\n");
    
    state.input = poly1_string;
    state.pos = 0;
    state.len = strlen(poly1_string);
    next_token(&state);
    parse_expression(&state, &poly1);
    
    state.input = poly2_string;
    state.pos = 0;
    state.len = strlen(poly2_string);
    if (state.current.str) {
        free(state.current.str);
        state.current.str = NULL;
    }
    next_token(&state);
    parse_expression(&state, &poly2);
    
    DEBUG_PRINT_R("Parsed polynomials:\n");
    DEBUG_PRINT_R("  poly1: %ld terms\n", poly1.nterms);
    DEBUG_PRINT_R("  poly2: %ld terms\n", poly2.nterms);
    
    // Create unified field context
    field_ctx_t field_ctx;
    field_ctx_init(&field_ctx, ctx);
    
    // Create unified mpoly context
    slong total_vars = 1 + state.npars; // elimination variable + parameters
    unified_mpoly_ctx_t unified_ctx = unified_mpoly_ctx_init(total_vars, ORD_LEX, &field_ctx);
    
    // Convert to unified mpoly format
    unified_mpoly_t A = unified_mpoly_init(unified_ctx);
    unified_mpoly_t B = unified_mpoly_init(unified_ctx);
    
    DEBUG_PRINT_R("Converting poly1 to unified format\n");
    // Convert poly1
    for (slong i = 0; i < poly1.nterms; i++) {
        field_elem_u coeff;
        ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
        
        fq_nmod_to_field_elem(&coeff, poly1.terms[i].coeff, &field_ctx);
        
        // Set elimination variable exponent
        if (poly1.terms[i].var_exp) {
            exp[0] = poly1.terms[i].var_exp[0];
        }
        
        // Set parameter exponents
        if (poly1.terms[i].par_exp) {
            for (slong j = 0; j < state.npars && j < poly1.npars; j++) {
                exp[1 + j] = poly1.terms[i].par_exp[j];
            }
        }
        
        unified_mpoly_set_coeff_ui(A, &coeff, exp);
        free(exp);
    }
    
    DEBUG_PRINT_R("Converting poly2 to unified format\n");
    // Convert poly2
    for (slong i = 0; i < poly2.nterms; i++) {
        field_elem_u coeff;
        ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
        
        fq_nmod_to_field_elem(&coeff, poly2.terms[i].coeff, &field_ctx);
        
        // Set elimination variable exponent
        if (poly2.terms[i].var_exp) {
            exp[0] = poly2.terms[i].var_exp[0];
        }
        
        // Set parameter exponents
        if (poly2.terms[i].par_exp) {
            for (slong j = 0; j < state.npars && j < poly2.npars; j++) {
                exp[1 + j] = poly2.terms[i].par_exp[j];
            }
        }
        
        unified_mpoly_set_coeff_ui(B, &coeff, exp);
        free(exp);
    }
    
    printf("\nStep 1: Compute resultant w.r.t. %s\n", elim_var);
    printf("  Poly1 has %ld terms\n", unified_mpoly_length(A));
    printf("  Poly2 has %ld terms\n", unified_mpoly_length(B));
    
    // Compute resultant
    unified_mpoly_t R = unified_mpoly_init(unified_ctx);
    clock_t start = clock();
    int success = unified_mpoly_resultant(R, A, B, 0, unified_ctx);
    clock_t end = clock();
    
    if (!success) {
        printf("Resultant computation failed!\n");
        unified_mpoly_clear(A);
        unified_mpoly_clear(B);
        unified_mpoly_clear(R);
        unified_mpoly_ctx_clear(unified_ctx);
        
        // Cleanup
        for (slong i = 0; i < state.nvars; i++) {
            free(state.var_names[i]);
        }
        free(state.var_names);
        for (slong i = 0; i < state.npars; i++) {
            free(state.par_names[i]);
        }
        free(state.par_names);
        for (slong i = 0; i < total_system_vars; i++) {
            free(all_system_vars[i]);
        }
        free(all_system_vars);
        if (state.generator_name) free(state.generator_name);
        fq_nmod_clear(state.current.value, ctx);
        if (state.current.str) free(state.current.str);
        fq_mvpoly_clear(&poly1);
        fq_mvpoly_clear(&poly2);
        
        return strdup("0");
    }
    
    printf("Resultant computation time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("Resultant has %ld terms\n", unified_mpoly_length(R));
    
    // Convert result back to fq_mvpoly first to inspect it
    fq_mvpoly_t debug_mvpoly;
    fq_mvpoly_init(&debug_mvpoly, 0, state.npars, ctx);
    
    DEBUG_PRINT_R("Converting result to fq_mvpoly for inspection\n");
    
    // Convert from unified format for debugging
    if (field_ctx.field_id == FIELD_ID_NMOD) {
        nmod_mpoly_struct *nmod_res = GET_NMOD_POLY(R);
        nmod_mpoly_ctx_struct *nmod_ctx = &(unified_ctx->ctx.nmod_ctx);
        
        for (slong i = 0; i < nmod_mpoly_length(nmod_res, nmod_ctx); i++) {
            ulong coeff_ui = nmod_mpoly_get_term_coeff_ui(nmod_res, i, nmod_ctx);
            ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
            nmod_mpoly_get_term_exp_ui(exp, nmod_res, i, nmod_ctx);
            
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_nmod_set_ui(coeff, coeff_ui, ctx);
            
            slong *par_exp = NULL;
            if (state.npars > 0) {
                par_exp = (slong*) calloc(state.npars, sizeof(slong));
                for (slong j = 0; j < state.npars; j++) {
                    par_exp[j] = exp[1 + j]; // Skip elimination variable at index 0
                }
            }
            
            fq_mvpoly_add_term_fast(&debug_mvpoly, NULL, par_exp, coeff);
            
            fq_nmod_clear(coeff, ctx);
            free(exp);
            if (par_exp) free(par_exp);
        }
    } else {
        fq_nmod_mpoly_struct *fq_res = GET_FQ_POLY(R);
        fq_nmod_mpoly_ctx_struct *fq_ctx = &(unified_ctx->ctx.fq_ctx);
        
        for (slong i = 0; i < fq_nmod_mpoly_length(fq_res, fq_ctx); i++) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, fq_res, i, fq_ctx);
            
            ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
            fq_nmod_mpoly_get_term_exp_ui(exp, fq_res, i, fq_ctx);
            
            slong *par_exp = NULL;
            if (state.npars > 0) {
                par_exp = (slong*) calloc(state.npars, sizeof(slong));
                for (slong j = 0; j < state.npars; j++) {
                    par_exp[j] = exp[1 + j];
                }
            }
            
            fq_mvpoly_add_term_fast(&debug_mvpoly, NULL, par_exp, coeff);
            
            fq_nmod_clear(coeff, ctx);
            free(exp);
            if (par_exp) free(par_exp);
        }
    }
    
    // Print the polynomial before reduction
    printf("\nBefore ideal reduction, resultant polynomial:\n");
    char *before_reduction = fq_mvpoly_to_string(&debug_mvpoly, state.par_names, state.generator_name);
    printf("  %s\n", before_reduction);
    free(before_reduction);
    
    // Apply ideal reduction if ideal is provided
    if (ideal && ideal->num_gens > 0) {
        printf("\nStep 2: Apply ideal reduction\n");
        
        // Create reduced ideal context for current parameters
        unified_triangular_ideal_t reduced_ideal;
        create_reduced_ideal_context(&reduced_ideal, ideal, state.par_names, 
                                   state.npars, ctx);
        
        printf("Reduced ideal has %ld generators for %ld parameters\n", 
               reduced_ideal.num_gens, state.npars);
        
        // Convert the resultant polynomial from unified context to parameter-only context for reduction
        if (field_ctx.field_id == FIELD_ID_NMOD && reduced_ideal.is_prime_field) {
            nmod_mpoly_struct *nmod_res = GET_NMOD_POLY(R);
            nmod_mpoly_ctx_struct *nmod_ctx = &(unified_ctx->ctx.nmod_ctx);
            
            printf("Converting from unified context (%ld vars) to parameter context (%ld vars)\n",
                   nmod_mpoly_ctx_nvars(nmod_ctx), nmod_mpoly_ctx_nvars(reduced_ideal.ctx.nmod_ctx));
            
            // Create polynomial in parameter-only context
            nmod_mpoly_t param_poly;
            nmod_mpoly_init(param_poly, reduced_ideal.ctx.nmod_ctx);
            
            // Convert terms from unified context to parameter context
            for (slong i = 0; i < nmod_mpoly_length(nmod_res, nmod_ctx); i++) {
                ulong coeff_ui = nmod_mpoly_get_term_coeff_ui(nmod_res, i, nmod_ctx);
                ulong *unified_exp = (ulong*) calloc(total_vars, sizeof(ulong));
                nmod_mpoly_get_term_exp_ui(unified_exp, nmod_res, i, nmod_ctx);
                
                // Extract parameter exponents (skip elimination variable at index 0)
                ulong *param_exp = (ulong*) calloc(state.npars, sizeof(ulong));
                for (slong j = 0; j < state.npars; j++) {
                    param_exp[j] = unified_exp[1 + j];
                }
                
                // Add term to parameter polynomial
                nmod_mpoly_t temp;
                nmod_mpoly_init(temp, reduced_ideal.ctx.nmod_ctx);
                nmod_mpoly_set_coeff_ui_ui(temp, coeff_ui, param_exp, reduced_ideal.ctx.nmod_ctx);
                nmod_mpoly_add(param_poly, param_poly, temp, reduced_ideal.ctx.nmod_ctx);
                nmod_mpoly_clear(temp, reduced_ideal.ctx.nmod_ctx);
                
                free(unified_exp);
                free(param_exp);
            }
            
            printf("Applying nmod reduction...\n");
            printf("Before reduction: %ld terms\n", nmod_mpoly_length(param_poly, reduced_ideal.ctx.nmod_ctx));
            
            // Apply reduction on parameter polynomial
            triangular_ideal_reduce_nmod_mpoly_with_names(param_poly, &reduced_ideal, state.par_names);
            
            printf("After reduction: %ld terms\n", nmod_mpoly_length(param_poly, reduced_ideal.ctx.nmod_ctx));
            
            // Convert back to unified format (copy reduced polynomial back)
            nmod_mpoly_zero(nmod_res, nmod_ctx);
            for (slong i = 0; i < nmod_mpoly_length(param_poly, reduced_ideal.ctx.nmod_ctx); i++) {
                ulong coeff_ui = nmod_mpoly_get_term_coeff_ui(param_poly, i, reduced_ideal.ctx.nmod_ctx);
                ulong *param_exp = (ulong*) calloc(state.npars, sizeof(ulong));
                nmod_mpoly_get_term_exp_ui(param_exp, param_poly, i, reduced_ideal.ctx.nmod_ctx);
                
                // Create unified exponent (elimination var = 0, then parameters)
                ulong *unified_exp = (ulong*) calloc(total_vars, sizeof(ulong));
                unified_exp[0] = 0; // elimination variable eliminated
                for (slong j = 0; j < state.npars; j++) {
                    unified_exp[1 + j] = param_exp[j];
                }
                
                nmod_mpoly_t temp;
                nmod_mpoly_init(temp, nmod_ctx);
                nmod_mpoly_set_coeff_ui_ui(temp, coeff_ui, unified_exp, nmod_ctx);
                nmod_mpoly_add(nmod_res, nmod_res, temp, nmod_ctx);
                nmod_mpoly_clear(temp, nmod_ctx);
                
                free(param_exp);
                free(unified_exp);
            }
            
            nmod_mpoly_clear(param_poly, reduced_ideal.ctx.nmod_ctx);
            
        } else if (field_ctx.field_id != FIELD_ID_NMOD && !reduced_ideal.is_prime_field) {
            fq_nmod_mpoly_struct *fq_res = GET_FQ_POLY(R);
            fq_nmod_mpoly_ctx_struct *fq_ctx = &(unified_ctx->ctx.fq_ctx);
            
            printf("Converting from unified context (%ld vars) to parameter context (%ld vars)\n",
                   fq_nmod_mpoly_ctx_nvars(fq_ctx), fq_nmod_mpoly_ctx_nvars(reduced_ideal.ctx.fq_ctx));
            
            // Create polynomial in parameter-only context
            fq_nmod_mpoly_t param_poly;
            fq_nmod_mpoly_init(param_poly, reduced_ideal.ctx.fq_ctx);
            
            // Convert terms from unified context to parameter context
            for (slong i = 0; i < fq_nmod_mpoly_length(fq_res, fq_ctx); i++) {
                fq_nmod_t coeff;
                fq_nmod_init(coeff, ctx);
                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, fq_res, i, fq_ctx);
                
                ulong *unified_exp = (ulong*) calloc(total_vars, sizeof(ulong));
                fq_nmod_mpoly_get_term_exp_ui(unified_exp, fq_res, i, fq_ctx);
                
                // Extract parameter exponents (skip elimination variable at index 0)
                ulong *param_exp = (ulong*) calloc(state.npars, sizeof(ulong));
                for (slong j = 0; j < state.npars; j++) {
                    param_exp[j] = unified_exp[1 + j];
                }
                
                // Add term to parameter polynomial
                fq_nmod_mpoly_set_coeff_fq_nmod_ui(param_poly, coeff, param_exp, reduced_ideal.ctx.fq_ctx);
                
                fq_nmod_clear(coeff, ctx);
                free(unified_exp);
                free(param_exp);
            }
            
            printf("Applying fq_nmod reduction...\n");
            printf("Before reduction: %ld terms\n", fq_nmod_mpoly_length(param_poly, reduced_ideal.ctx.fq_ctx));
            
            // Apply reduction on parameter polynomial
            triangular_ideal_reduce_fq_nmod_mpoly_with_names(param_poly, &reduced_ideal, state.par_names);
            
            printf("After reduction: %ld terms\n", fq_nmod_mpoly_length(param_poly, reduced_ideal.ctx.fq_ctx));
            
            // Convert back to unified format
            fq_nmod_mpoly_zero(fq_res, fq_ctx);
            for (slong i = 0; i < fq_nmod_mpoly_length(param_poly, reduced_ideal.ctx.fq_ctx); i++) {
                fq_nmod_t coeff;
                fq_nmod_init(coeff, ctx);
                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, param_poly, i, reduced_ideal.ctx.fq_ctx);
                
                ulong *param_exp = (ulong*) calloc(state.npars, sizeof(ulong));
                fq_nmod_mpoly_get_term_exp_ui(param_exp, param_poly, i, reduced_ideal.ctx.fq_ctx);
                
                // Create unified exponent
                ulong *unified_exp = (ulong*) calloc(total_vars, sizeof(ulong));
                unified_exp[0] = 0; // elimination variable eliminated
                for (slong j = 0; j < state.npars; j++) {
                    unified_exp[1 + j] = param_exp[j];
                }
                
                fq_nmod_mpoly_set_coeff_fq_nmod_ui(fq_res, coeff, unified_exp, fq_ctx);
                
                fq_nmod_clear(coeff, ctx);
                free(param_exp);
                free(unified_exp);
            }
            
            fq_nmod_mpoly_clear(param_poly, reduced_ideal.ctx.fq_ctx);
            
        } else {
            printf("Warning: Field context mismatch, skipping reduction\n");
        }
        
        unified_triangular_ideal_clear(&reduced_ideal);
    } else {
        printf("\nStep 2: Skipping ideal reduction (no ideal provided)\n");
    }
    
    fq_mvpoly_clear(&debug_mvpoly);
    
    // Convert result back to fq_mvpoly
    fq_mvpoly_t result_mvpoly;
    fq_mvpoly_init(&result_mvpoly, 0, state.npars, ctx);
    
    DEBUG_PRINT_R("Converting result back to fq_mvpoly format\n");
    
    // Convert from unified format based on field type
    if (field_ctx.field_id == FIELD_ID_NMOD) {
        nmod_mpoly_struct *nmod_res = GET_NMOD_POLY(R);
        nmod_mpoly_ctx_struct *nmod_ctx = &(unified_ctx->ctx.nmod_ctx);
        
        for (slong i = 0; i < nmod_mpoly_length(nmod_res, nmod_ctx); i++) {
            ulong coeff_ui = nmod_mpoly_get_term_coeff_ui(nmod_res, i, nmod_ctx);
            ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
            nmod_mpoly_get_term_exp_ui(exp, nmod_res, i, nmod_ctx);
            
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_nmod_set_ui(coeff, coeff_ui, ctx);
            
            slong *par_exp = NULL;
            if (state.npars > 0) {
                par_exp = (slong*) calloc(state.npars, sizeof(slong));
                for (slong j = 0; j < state.npars; j++) {
                    par_exp[j] = exp[1 + j]; // Skip elimination variable at index 0
                }
            }
            
            fq_mvpoly_add_term_fast(&result_mvpoly, NULL, par_exp, coeff);
            
            fq_nmod_clear(coeff, ctx);
            free(exp);
            if (par_exp) free(par_exp);
        }
    } else if (field_ctx.field_id == FIELD_ID_FQ_ZECH) {
        fq_zech_mpoly_struct *zech_res = GET_ZECH_POLY(R);
        fq_zech_mpoly_ctx_struct *zech_ctx = &(unified_ctx->ctx.zech_ctx);
        
        for (slong i = 0; i < fq_zech_mpoly_length(zech_res, zech_ctx); i++) {
            fq_zech_t zech_coeff;
            fq_zech_init(zech_coeff, zech_ctx->fqctx);
            fq_zech_mpoly_get_term_coeff_fq_zech(zech_coeff, zech_res, i, zech_ctx);
            
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_zech_get_fq_nmod(coeff, zech_coeff, zech_ctx->fqctx);
            
            ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
            fq_zech_mpoly_get_term_exp_ui(exp, zech_res, i, zech_ctx);
            
            slong *par_exp = NULL;
            if (state.npars > 0) {
                par_exp = (slong*) calloc(state.npars, sizeof(slong));
                for (slong j = 0; j < state.npars; j++) {
                    par_exp[j] = exp[1 + j];
                }
            }
            
            fq_mvpoly_add_term_fast(&result_mvpoly, NULL, par_exp, coeff);
            
            fq_nmod_clear(coeff, ctx);
            fq_zech_clear(zech_coeff, zech_ctx->fqctx);
            free(exp);
            if (par_exp) free(par_exp);
        }
    } else {
        fq_nmod_mpoly_struct *fq_res = GET_FQ_POLY(R);
        fq_nmod_mpoly_ctx_struct *fq_ctx = &(unified_ctx->ctx.fq_ctx);
        
        for (slong i = 0; i < fq_nmod_mpoly_length(fq_res, fq_ctx); i++) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, fq_res, i, fq_ctx);
            
            ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
            fq_nmod_mpoly_get_term_exp_ui(exp, fq_res, i, fq_ctx);
            
            slong *par_exp = NULL;
            if (state.npars > 0) {
                par_exp = (slong*) calloc(state.npars, sizeof(slong));
                for (slong j = 0; j < state.npars; j++) {
                    par_exp[j] = exp[1 + j];
                }
            }
            
            fq_mvpoly_add_term_fast(&result_mvpoly, NULL, par_exp, coeff);
            
            fq_nmod_clear(coeff, ctx);
            free(exp);
            if (par_exp) free(par_exp);
        }
    }
    
    // Convert result to string
    char *result = fq_mvpoly_to_string(&result_mvpoly, state.par_names, state.generator_name);
    
    // Print remaining variables
    printf("Remaining variables: ");
    if (state.npars == 0) {
        printf("none");
    } else {
        for (slong i = 0; i < state.npars; i++) {
            if (i > 0) printf(", ");
            printf("%s", state.par_names[i]);
        }
    }
    printf("\n");
    
    // Cleanup
    unified_mpoly_clear(A);
    unified_mpoly_clear(B);
    unified_mpoly_clear(R);
    unified_mpoly_ctx_clear(unified_ctx);
    fq_mvpoly_clear(&poly1);
    fq_mvpoly_clear(&poly2);
    fq_mvpoly_clear(&result_mvpoly);
    
    for (slong i = 0; i < state.nvars; i++) {
        free(state.var_names[i]);
    }
    free(state.var_names);
    
    for (slong i = 0; i < state.npars; i++) {
        free(state.par_names[i]);
    }
    free(state.par_names);
    
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
    
    printf("\n=== Resultant with Ideal Reduction Complete ===\n");
    
    return result;
}

// String interface for resultant with ideal reduction
char* resultant_with_ideal_reduction_str(const char *poly1_string,
                                       const char *poly2_string, 
                                       const char *elim_var_string,
                                       const char *ideal_gens_string,
                                       const fq_nmod_ctx_t ctx) {
    
    slong num_gens, num_all_vars;
    char **gens_array = split_string_r(ideal_gens_string, &num_gens);
    
    unified_triangular_ideal_t ideal;
    const char **ideal_gens = (const char**) malloc(num_gens * sizeof(char*));
    const char **all_var_names = (const char**) malloc(num_all_vars * sizeof(char*));

    for (slong i = 0; i < num_gens; i++) {
        ideal_gens[i] = gens_array[i];
    }
    
    char **all_vars_auto = extract_all_variables_from_ideal_gens(
        ideal_gens, num_gens, ctx, &num_all_vars);
    for (slong i = 0; i < num_all_vars; i++) {
        all_var_names[i] = all_vars_auto[i];
    }
    
    construct_triangular_ideal_from_strings(&ideal, ideal_gens, num_gens,
                                          all_var_names, num_all_vars, ctx);
    
    char *result = resultant_with_ideal_reduction(poly1_string, poly2_string,
                                                elim_var_string, ctx, &ideal);
    
    unified_triangular_ideal_clear(&ideal);
    free(ideal_gens);
    free(all_var_names);
    free_split_string_rs(gens_array, num_gens);
    
    return result;
}

char* elimination_with_ideal_reduction_str(const char *poly_string,
                                               const char *elim_vars_string,
                                               const char *ideal_gens_string,
                                               const fq_nmod_ctx_t ctx) {
    
    printf("\n=== Smart Elimination with Ideal Reduction ===\n");
    
    // Parse input polynomials
    slong num_polys;
    char **poly_array = split_string_r(poly_string, &num_polys);
    
    printf("Detected %ld polynomial(s)\n", num_polys);
    
    char *result = NULL;
    
    if (num_polys == 2) {
        // Use resultant method
        printf("Using bivariate resultant method\n");
        
        // Parse elimination variables - resultant only supports single variable elimination
        slong num_elim_vars;
        char **elim_array = split_string_r(elim_vars_string, &num_elim_vars);
        
        if (num_elim_vars != 1) {
            printf("Warning: Resultant method requires exactly one elimination variable, "
                   "using first variable: %s\n", elim_array[0]);
        }
        
        const char *elim_var = (num_elim_vars > 0) ? elim_array[0] : "x";
        
        // Call resultant function
        result = resultant_with_ideal_reduction_str(poly_array[0], poly_array[1], 
                                                  elim_var, ideal_gens_string, ctx);
        
        // Cleanup
        free_split_string_rs(elim_array, num_elim_vars);
        
    } else if (num_polys > 2) {
        // Use Dixon method
        printf("Using Dixon resultant method for %ld polynomials\n", num_polys);
        
        // Call dixon_with_ideal_reduction_str
        result = dixon_with_ideal_reduction_str(poly_string, elim_vars_string,
                                              ideal_gens_string, ctx);
        
    } else {
        // No polynomials
        printf("No polynomials provided\n");
        result = strdup("1");
    }
    
    // Cleanup polynomial array
    free_split_string_rs(poly_array, num_polys);
    
    printf("\n=== Smart Elimination Complete ===\n");
    
    return result;
}

char* elimination_with_ideal(const char **poly_strings,
                           slong num_polys,
                           const char **elim_vars,
                           slong num_elim_vars,
                           const char *ideal_string,
                           const fq_nmod_ctx_t ctx) {
    
    clock_t start_time = clock();
    char *result;
    if (num_polys == 2) {
        // Use resultant method
        printf("Using bivariate resultant method\n");
        // Call resultant function
        result = resultant_with_ideal_reduction_str(poly_strings[0], poly_strings[1], 
                                                  elim_vars[0], ideal_string, ctx);
    }
    else {
        printf("\n=== Dixon with Ideal Reduction ===\n");
        printf("Processing %ld polynomials, eliminating %ld variables\n", num_polys, num_elim_vars);
        printf("Ideal string: %s\n", ideal_string);
        
        // Print elimination variables
        printf("Eliminating variables: ");
        for (slong i = 0; i < num_elim_vars; i++) {
            if (i > 0) printf(", ");
            printf("%s", elim_vars[i]);
        }
        printf("\n");
        
        // Print polynomials
        printf("Input polynomials:\n");
        for (slong i = 0; i < num_polys; i++) {
            printf("  p%ld: %s\n", i, poly_strings[i]);
        }
        
        // Step 1: Split ideal string into individual equations
        slong num_ideal_equations;
        char **ideal_equations_array = split_string_r(ideal_string, &num_ideal_equations);
        
        printf("\nExtracted %ld equations from ideal:\n", num_ideal_equations);
        for (slong i = 0; i < num_ideal_equations; i++) {
            printf("  eq%ld: %s\n", i, ideal_equations_array[i]);
        }
        
        // Convert to const pointer array
        const char **ideal_equations = (const char**) malloc(num_ideal_equations * sizeof(char*));
        for (slong i = 0; i < num_ideal_equations; i++) {
            ideal_equations[i] = ideal_equations_array[i];
        }
        
        // Step 2: Auto-extract all variables from ideal
        slong num_ideal_vars;
        char **ideal_vars_auto = extract_all_variables_from_ideal_gens(
            ideal_equations, num_ideal_equations, ctx, &num_ideal_vars);
        
        printf("\nAuto-extracted %ld variables from ideal: ", num_ideal_vars);
        for (slong i = 0; i < num_ideal_vars; i++) {
            if (i > 0) printf(", ");
            printf("%s", ideal_vars_auto[i]);
        }
        printf("\n");
        
        // Convert to const pointer array
        const char **ideal_vars = (const char**) malloc(num_ideal_vars * sizeof(char*));
        for (slong i = 0; i < num_ideal_vars; i++) {
            ideal_vars[i] = ideal_vars_auto[i];
        }
        
        // Step 3: Construct triangular ideal
        unified_triangular_ideal_t ideal;
        construct_triangular_ideal_from_strings(&ideal, ideal_equations, num_ideal_equations,
                                              ideal_vars, num_ideal_vars, ctx);
        
        printf("\nIdeal construction complete:\n");
        printf("  Number of generators: %ld\n", ideal.num_gens);
        for (slong i = 0; i < ideal.num_gens; i++) {
            printf("  Generator %ld: %s^%ld (index %ld)\n", 
                   i, ideal.var_names[i] ? ideal.var_names[i] : "?", 
                   ideal.leading_degrees[i], ideal.var_indices[i]);
        }
        
        // Step 4: Compute Dixon resultant
        printf("\nComputing Dixon resultant with ideal reduction...\n");
        result = dixon_with_ideal_reduction(poly_strings, num_polys,
                                                elim_vars, num_elim_vars,
                                                ctx, &ideal);
        
        // Cleanup
        unified_triangular_ideal_clear(&ideal);
        free(ideal_equations);
        free(ideal_vars);
        
        // Cleanup auto-extracted variables
        for (slong i = 0; i < num_ideal_vars; i++) {
            free(ideal_vars_auto[i]);
        }
        free(ideal_vars_auto);
        
        // Cleanup split equations
        free_split_string_rs(ideal_equations_array, num_ideal_equations);
        
        printf("=== Dixon Computation Complete ===\n\n");
    }
    clock_t end_time = clock();
    double elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Time: %.3f seconds\n", elapsed);
    return result;
}

// Unified elimination with ideal reduction - automatically chooses Dixon or resultant
char* elimination_with_ideal_reduction(const char **poly_strings, slong num_polys,
                                      const char **elim_vars, slong num_elim_vars,
                                      const fq_nmod_ctx_t ctx,
                                      unified_triangular_ideal_t *ideal) {
    printf("\n=== Unified Elimination with Ideal Reduction ===\n");
    printf("Eliminating variables: ");
    for (slong i = 0; i < num_elim_vars; i++) {
        if (i > 0) printf(", ");
        printf("%s", elim_vars[i]);
    }
    printf("\n");
    
    printf("Input: %ld polynomial(s), %ld elimination variable(s)\n", 
           num_polys, num_elim_vars);
    
    char *result = NULL;
    
    // Case 1: Two polynomials and one elimination variable -> use bivariate resultant
    if (num_polys == 2 && num_elim_vars == 1) {
        printf("Strategy: Using bivariate resultant method\n");
        result = resultant_with_ideal_reduction(poly_strings[0], poly_strings[1],
                                              elim_vars[0], ctx, ideal);
    }
    // Case 2: Multiple polynomials or multiple variables -> use Dixon method
    else if (num_polys > 2 || num_elim_vars > 1) {
        printf("Strategy: Using Dixon resultant method\n");
        result = dixon_with_ideal_reduction(poly_strings, num_polys,
                                          elim_vars, num_elim_vars,
                                          ctx, ideal);
    }
    // Case 3: Single polynomial -> apply ideal reduction directly
    else if (num_polys == 1) {
        printf("Strategy: Single polynomial - applying ideal reduction only\n");

        result = strdup(poly_strings[0]); // Simplified - should apply reduction
        
        printf("Note: Single polynomial ideal reduction not fully implemented\n");
    }
    // Case 4: No polynomials
    else {
        printf("Strategy: No polynomials provided - returning unit element\n");
        result = strdup("1");
    }
    
    printf("\n=== Unified Elimination Complete ===\n");
    
    return result;
}

void test_iterative_elimination_str2(void) {
    printf("\n================================================\n");
    printf("Test: Simplified Iterative Elimination\n");
    printf("================================================\n");
    
    fq_nmod_ctx_t ctx;
    fmpz_t p;
    fmpz_init(p);
    fmpz_set_ui(p, 2147483489);
    fq_nmod_ctx_init(ctx, p, 1, "t");
    fmpz_clear(p);
    
    // Define ideal generators in equation format
    const char *ideal_gens_str = "a2^3 = 2*a1 + 1; a3^3 = a1*a2 + 3; a4^3 = a1 + a2*a3 + 5";
    const char *all_vars_str = "a1, a2, a3, a4";
    
    // Step 1: Eliminate a4
    const char *step1_polys_str = "a1^2 + a2^2 + a3^2 + a4^2 - 100, a4^3 - a1 - a2*a3 - 5";
    const char *step1_elim_str = "a4";
    
    char *result1 = elimination_with_ideal_reduction_str(step1_polys_str, step1_elim_str, 
                                                  ideal_gens_str, ctx);
    printf("After eliminating a4: %s\n\n", result1);
    
    // Step 2: Eliminate a3
    char step2_polys_str[10000];
    snprintf(step2_polys_str, sizeof(step2_polys_str), "%s, a3^3 - a1*a2 - 3", result1);
    const char *step2_elim_str = "a3";
    
    char *result2 = elimination_with_ideal_reduction_str(step2_polys_str, step2_elim_str,
                                                  ideal_gens_str, ctx);
    printf("After eliminating a3: %s\n\n", result2);
    
    // Step 3: Eliminate a2
    char step3_polys_str[10000];
    snprintf(step3_polys_str, sizeof(step3_polys_str), "%s, a2^3 - 2*a1 - 1", result2);
    const char *step3_elim_str = "a2";
    
    char *final_result = elimination_with_ideal_reduction_str(step3_polys_str, step3_elim_str,
                                                       ideal_gens_str, ctx);
    printf("Final univariate polynomial in a1: %s\n", final_result);
    
    // Cleanup
    free(result1);
    free(result2);
    free(final_result);
    fq_nmod_ctx_clear(ctx);
}
