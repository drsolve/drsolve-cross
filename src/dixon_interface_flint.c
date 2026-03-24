// Complete fixed Dixon resultant string interface implementation
#include "dixon_interface_flint.h"
#include "dixon_complexity.h"
#include <flint/arith.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly.h>
#include <fcntl.h>
#include <unistd.h>

// Fixed string parser implementation

static int g_suppress_univariate_root_reporting = 0;
static const slong QQ_ROOT_SEARCH_MAX_DEGREE = 1000;
static const slong QQ_ROOT_SEARCH_MAX_CANDIDATES = 1000000;

static int at_end(parser_state_t *state) {
    return state->pos >= state->len;
}

static char peek(parser_state_t *state) {
    if (at_end(state)) return '\0';
    return state->input[state->pos];
}

static char advance(parser_state_t *state) {
    if (at_end(state)) return '\0';
    return state->input[state->pos++];
}

static void skip_whitespace(parser_state_t *state) {
    while (!at_end(state) && isspace(peek(state))) {
        advance(state);
    }
}

static void parse_number(parser_state_t *state) {
    size_t start = state->pos;
    
    while (!at_end(state) && isdigit(peek(state))) {
        advance(state);
    }
    
    size_t len = state->pos - start;
    state->current.str = (char*) malloc(len + 1);
    strncpy(state->current.str, state->input + start, len);
    state->current.str[len] = '\0';

    fmpz_t parsed_value;
    fmpz_init(parsed_value);

    if (fmpz_set_str(parsed_value, state->current.str, 10) == 0) {
        fq_nmod_set_fmpz(state->current.value, parsed_value, state->ctx);
        if (fmpz_fits_si(parsed_value)) {
            state->current.int_value = fmpz_get_si(parsed_value);
            state->current.int_value_valid = 1;
        } else {
            state->current.int_value = 0;
            state->current.int_value_valid = 0;
        }
    } else {
        state->current.int_value = 0;
        state->current.int_value_valid = 0;
        fq_nmod_zero(state->current.value, state->ctx);
    }

    fmpz_clear(parsed_value);
    state->current.type = TOK_NUMBER;
    
    DEBUG_PRINT("Number: %s (%ld)%s\n", state->current.str, state->current.int_value,
                state->current.int_value_valid ? "" : " [no slong fit]");
}

static void parse_identifier(parser_state_t *state) {
    size_t start = state->pos;
    
    while (!at_end(state) && (isalnum(peek(state)) || peek(state) == '_')) {
        advance(state);
    }
    
    size_t len = state->pos - start;
    state->current.str = (char*) malloc(len + 1);
    strncpy(state->current.str, state->input + start, len);
    state->current.str[len] = '\0';
    
    if (state->generator_name && strcmp(state->current.str, state->generator_name) == 0) {
        state->current.type = TOK_GENERATOR;
        fq_nmod_gen(state->current.value, state->ctx);
        DEBUG_PRINT("Generator: %s\n", state->current.str);
    } else {
        state->current.type = TOK_VARIABLE;
        DEBUG_PRINT("Identifier: %s\n", state->current.str);
    }
}

void next_token(parser_state_t *state) {
    skip_whitespace(state);
    state->current.int_value = 0;
    state->current.int_value_valid = 0;
    
    if (state->current.str) {
        free(state->current.str);
        state->current.str = NULL;
    }
    
    if (at_end(state)) {
        state->current.type = TOK_EOF;
        state->current.int_value_valid = 0;
        return;
    }
    
    char c = peek(state);
    
    switch (c) {
        case '+': advance(state); state->current.type = TOK_PLUS; return;
        case '-': advance(state); state->current.type = TOK_MINUS; return;
        case '*': advance(state); state->current.type = TOK_MULT; return;
        case '^': advance(state); state->current.type = TOK_POWER; return;
        case '(': advance(state); state->current.type = TOK_LPAREN; return;
        case ')': advance(state); state->current.type = TOK_RPAREN; return;
    }
    
    if (isdigit(c)) {
        parse_number(state);
    } else if (isalpha(c) || c == '_') {
        parse_identifier(state);
    } else {
        advance(state);
        next_token(state);
    }
}

static slong find_or_add_parameter(parser_state_t *state, const char *name) {
    // Check if it's a variable
    for (slong i = 0; i < state->nvars; i++) {
        if (strcmp(state->var_names[i], name) == 0) {
            return -1;
        }
    }
    
    // Check existing parameters
    for (slong i = 0; i < state->npars; i++) {
        if (strcmp(state->par_names[i], name) == 0) {
            return i;
        }
    }
    
    // Add new parameter
    if (state->npars >= state->max_pars) {
        state->max_pars *= 2;
        state->par_names = (char**) realloc(state->par_names, 
                                           state->max_pars * sizeof(char*));
    }
    
    state->par_names[state->npars] = strdup(name);
    DEBUG_PRINT("New parameter: %s (index %ld)\n", name, state->npars);
    return state->npars++;
}

static slong get_variable_index(parser_state_t *state, const char *name) {
    for (slong i = 0; i < state->nvars; i++) {
        if (strcmp(state->var_names[i], name) == 0) {
            return i;
        }
    }
    return -1;
}

void parse_primary(parser_state_t *state, fq_mvpoly_t *poly) {
    if (state->current.type == TOK_NUMBER) {
        fq_mvpoly_add_term(poly, NULL, NULL, state->current.value);
        next_token(state);
        
    } else if (state->current.type == TOK_GENERATOR) {
        fq_mvpoly_add_term(poly, NULL, NULL, state->current.value);
        next_token(state);
        
    } else if (state->current.type == TOK_VARIABLE) {
        char *name = strdup(state->current.str);
        next_token(state);
        
        slong var_idx = get_variable_index(state, name);
        if (var_idx >= 0) {
            slong *var_exp = (slong*) calloc(state->nvars, sizeof(slong));
            var_exp[var_idx] = 1;
            
            fq_nmod_t one;
            fq_nmod_init(one, state->ctx);
            fq_nmod_one(one, state->ctx);
            fq_mvpoly_add_term(poly, var_exp, NULL, one);
            fq_nmod_clear(one, state->ctx);
            free(var_exp);
        } else {
            slong par_idx = find_or_add_parameter(state, name);
            if (par_idx >= 0) {
                slong *par_exp = (slong*) calloc(state->max_pars, sizeof(slong));
                par_exp[par_idx] = 1;
                
                fq_nmod_t one;
                fq_nmod_init(one, state->ctx);
                fq_nmod_one(one, state->ctx);
                fq_mvpoly_add_term(poly, NULL, par_exp, one);
                fq_nmod_clear(one, state->ctx);
                free(par_exp);
            }
        }
        free(name);
        
    } else if (state->current.type == TOK_LPAREN) {
        next_token(state);
        parse_expression(state, poly);
        if (state->current.type == TOK_RPAREN) {
            next_token(state);
        }
        
    } else if (state->current.type == TOK_MINUS) {
        next_token(state);
        fq_mvpoly_t temp;
        fq_mvpoly_init(&temp, state->nvars, state->max_pars, state->ctx);
        parse_primary(state, &temp);
        
        for (slong i = 0; i < temp.nterms; i++) {
            fq_nmod_t neg_coeff;
            fq_nmod_init(neg_coeff, state->ctx);
            fq_nmod_neg(neg_coeff, temp.terms[i].coeff, state->ctx);
            fq_mvpoly_add_term(poly, temp.terms[i].var_exp, temp.terms[i].par_exp, neg_coeff);
            fq_nmod_clear(neg_coeff, state->ctx);
        }
        fq_mvpoly_clear(&temp);
    }
}

void parse_factor(parser_state_t *state, fq_mvpoly_t *poly) {
    fq_mvpoly_t base;
    fq_mvpoly_init(&base, state->nvars, state->max_pars, state->ctx);
    parse_primary(state, &base);
    
    if (state->current.type == TOK_POWER) {
        next_token(state);
        if (state->current.type == TOK_NUMBER) {
            if (!state->current.int_value_valid) {
                fprintf(stderr, "Error: exponent '%s' exceeds supported slong range\n",
                        state->current.str ? state->current.str : "<unknown>");
                fq_mvpoly_clear(&base);
                return;
            }
            slong exp = state->current.int_value;
            next_token(state);
            
            fq_mvpoly_t result;
            fq_mvpoly_pow(&result, &base, exp);
            
            for (slong i = 0; i < result.nterms; i++) {
                fq_mvpoly_add_term(poly, result.terms[i].var_exp, result.terms[i].par_exp, result.terms[i].coeff);
            }
            fq_mvpoly_clear(&result);
        }
    } else {
        for (slong i = 0; i < base.nterms; i++) {
            fq_mvpoly_add_term(poly, base.terms[i].var_exp, base.terms[i].par_exp, base.terms[i].coeff);
        }
    }
    
    fq_mvpoly_clear(&base);
}

void parse_term(parser_state_t *state, fq_mvpoly_t *poly) {
    fq_mvpoly_t result;
    fq_mvpoly_init(&result, state->nvars, state->max_pars, state->ctx);
    
    parse_factor(state, &result);
    
    while (state->current.type == TOK_MULT) {
        next_token(state);
        
        fq_mvpoly_t factor;
        fq_mvpoly_init(&factor, state->nvars, state->max_pars, state->ctx);
        parse_factor(state, &factor);
        
        fq_mvpoly_t temp;
        fq_mvpoly_mul(&temp, &result, &factor);
        
        fq_mvpoly_clear(&result);
        fq_mvpoly_clear(&factor);
        fq_mvpoly_copy(&result, &temp);
        fq_mvpoly_clear(&temp);
    }
    
    for (slong i = 0; i < result.nterms; i++) {
        fq_mvpoly_add_term_fast(poly, result.terms[i].var_exp, result.terms[i].par_exp, result.terms[i].coeff);
    }
    fq_mvpoly_clear(&result);
}

void parse_expression(parser_state_t *state, fq_mvpoly_t *poly) {
    int negate = 0;
    if (state->current.type == TOK_MINUS) {
        negate = 1;
        next_token(state);
    } else if (state->current.type == TOK_PLUS) {
        next_token(state);
    }
    
    fq_mvpoly_t first_term;
    fq_mvpoly_init(&first_term, state->nvars, state->max_pars, state->ctx);
    parse_term(state, &first_term);
    
    if (negate) {
        for (slong i = 0; i < first_term.nterms; i++) {
            fq_nmod_t neg_coeff;
            fq_nmod_init(neg_coeff, state->ctx);
            fq_nmod_neg(neg_coeff, first_term.terms[i].coeff, state->ctx);
            fq_mvpoly_add_term_fast(poly, first_term.terms[i].var_exp, first_term.terms[i].par_exp, neg_coeff);
            fq_nmod_clear(neg_coeff, state->ctx);
        }
    } else {
        for (slong i = 0; i < first_term.nterms; i++) {
            fq_mvpoly_add_term_fast(poly, first_term.terms[i].var_exp, first_term.terms[i].par_exp, 
                           first_term.terms[i].coeff);
        }
    }
    fq_mvpoly_clear(&first_term);
    
    while (state->current.type == TOK_PLUS || state->current.type == TOK_MINUS) {
        int subtract = (state->current.type == TOK_MINUS);
        next_token(state);
        
        fq_mvpoly_t term;
        fq_mvpoly_init(&term, state->nvars, state->max_pars, state->ctx);
        parse_term(state, &term);
        
        for (slong i = 0; i < term.nterms; i++) {
            if (subtract) {
                fq_nmod_t neg_coeff;
                fq_nmod_init(neg_coeff, state->ctx);
                fq_nmod_neg(neg_coeff, term.terms[i].coeff, state->ctx);
                fq_mvpoly_add_term_fast(poly, term.terms[i].var_exp, term.terms[i].par_exp, neg_coeff);
                fq_nmod_clear(neg_coeff, state->ctx);
            } else {
                fq_mvpoly_add_term_fast(poly, term.terms[i].var_exp, term.terms[i].par_exp, term.terms[i].coeff);
            }
        }
        fq_mvpoly_clear(&term);
    }
}

void find_and_print_roots_of_univariate_resultant(const fq_mvpoly_t *result, parser_state_t *state) {
    if (g_suppress_univariate_root_reporting) {
        return;
    }

    if (result->nvars != 0) {
        return;
    }
    
    int *par_used = (int*) calloc(result->npars, sizeof(int));
    slong actual_par_count = 0;
    slong main_par_idx = -1;
    
    for (slong t = 0; t < result->nterms; t++) {
        if (result->terms[t].par_exp) {
            for (slong p = 0; p < result->npars; p++) {
                if (result->terms[t].par_exp[p] > 0 && !par_used[p]) {
                    par_used[p] = 1;
                    main_par_idx = p;
                    actual_par_count++;
                }
            }
        }
    }
    
    if (actual_par_count > 1) {
        free(par_used);
        return;
    }
    
    if (actual_par_count == 0) {
        free(par_used);
        return;
    }
    
    fq_nmod_poly_t poly;
    fq_nmod_poly_init(poly, result->ctx);
    
    for (slong i = 0; i < result->nterms; i++) {
        slong degree = 0;
        if (result->terms[i].par_exp && result->terms[i].par_exp[main_par_idx] > 0) {
            degree = result->terms[i].par_exp[main_par_idx];
        }
        fq_nmod_poly_set_coeff(poly, degree, result->terms[i].coeff, result->ctx);
    }
    
    const char *var_name = "unknown";
    if (state->par_names && state->par_names[main_par_idx]) {
        var_name = state->par_names[main_par_idx];
    }
    
    slong degree = fq_nmod_poly_degree(poly, result->ctx);
    if (degree <= 0) {
        fq_nmod_poly_clear(poly, result->ctx);
        free(par_used);
        return;
    }

    slong field_degree = fq_nmod_ctx_degree(result->ctx);
    if (field_degree == 1) {
        mp_limb_t prime = fq_nmod_ctx_prime(result->ctx);
        
        nmod_poly_t nmod_poly;
        nmod_poly_init(nmod_poly, prime);
        
        for (slong i = 0; i <= degree; i++) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, result->ctx);
            fq_nmod_poly_get_coeff(coeff, poly, i, result->ctx);
            
            nmod_poly_t temp_poly;
            nmod_poly_init(temp_poly, prime);
            fq_nmod_get_nmod_poly(temp_poly, coeff, result->ctx);
            
            mp_limb_t coeff_ui = 0;
            if (nmod_poly_degree(temp_poly) >= 0) {
                coeff_ui = nmod_poly_get_coeff_ui(temp_poly, 0);
            }
            
            nmod_poly_set_coeff_ui(nmod_poly, i, coeff_ui);
            
            fq_nmod_clear(coeff, result->ctx);
            nmod_poly_clear(temp_poly);
        }
        
        nmod_roots_t nmod_roots;
        nmod_roots_init(nmod_roots);
        slong num_roots = our_nmod_poly_roots(nmod_roots, nmod_poly, 1);

        printf("\nRoots in %s (degree %ld):\n", var_name, degree);
        printf("Find %ld roots:\n", num_roots);
        for (slong i = 0; i < nmod_roots->num; i++) {
            printf("  Root %ld: %lu (Multiplicity: %ld)\n", i + 1, 
                   nmod_roots->roots[i], nmod_roots->mult[i]);
        }
        
        nmod_roots_clear(nmod_roots);
        nmod_poly_clear(nmod_poly);
        
    } else {
        fq_nmod_roots_t roots;
        fq_nmod_roots_init(roots, result->ctx);
        slong num_roots = our_fq_nmod_poly_roots(roots, poly, 1, result->ctx);

        printf("\nRoots in %s (degree %ld):\n", var_name, degree);
        printf("Find %ld roots:\n", num_roots);
        for (slong i = 0; i < roots->num; i++) {
            printf("  Root %ld: ", i + 1);
            fq_nmod_print_pretty(roots->roots + i, result->ctx);
            printf(" (Multiplicity: %ld)\n", roots->mult[i]);
        }
        
        fq_nmod_roots_clear(roots, result->ctx);
    }

    fq_nmod_poly_clear(poly, result->ctx);
    free(par_used);
}

// Helper functions

// Get field generator name
char* get_generator_name(const fq_nmod_ctx_t ctx) {
    return strdup(ctx->var);
}

// Convert fq_nmod_t to string with field generator
char* fq_nmod_to_string_with_gen(const fq_nmod_t elem, const fq_nmod_ctx_t ctx, const char *gen_name) {
    // Initially allocate a large buffer
    size_t capacity = 256;
    char *buffer = (char*) malloc(capacity);
    buffer[0] = '\0';
    
    if (fq_nmod_is_zero(elem, ctx)) {
        strcpy(buffer, "0");
        return buffer;
    }
    
    if (fq_nmod_ctx_degree(ctx) == 1) {
        // Prime field element
        nmod_poly_t poly;
        nmod_poly_init(poly, fq_nmod_ctx_prime(ctx));
        fq_nmod_get_nmod_poly(poly, elem, ctx);
        
        if (nmod_poly_degree(poly) >= 0) {
            sprintf(buffer, "%lu", nmod_poly_get_coeff_ui(poly, 0));
        } else {
            strcpy(buffer, "0");
        }
        nmod_poly_clear(poly);
    } else {
        // Extension field element
        nmod_poly_t poly;
        nmod_poly_init(poly, fq_nmod_ctx_prime(ctx));
        fq_nmod_get_nmod_poly(poly, elem, ctx);
        
        slong deg = nmod_poly_degree(poly);
        int first = 1;
        size_t len = 0;
        
        // Dynamically build string
        strcat(buffer, "(");
        len = 1;
        
        for (slong i = deg; i >= 0; i--) {
            mp_limb_t coeff = nmod_poly_get_coeff_ui(poly, i);
            if (coeff != 0) {
                char temp[128];
                
                if (!first) {
                    sprintf(temp, " + ");
                } else {
                    temp[0] = '\0';
                }
                first = 0;
                
                if (i == 0) {
                    sprintf(temp + strlen(temp), "%lu", coeff);
                } else if (i == 1) {
                    if (coeff == 1) {
                        sprintf(temp + strlen(temp), "%s", gen_name);
                    } else {
                        sprintf(temp + strlen(temp), "%lu*%s", coeff, gen_name);
                    }
                } else {
                    if (coeff == 1) {
                        sprintf(temp + strlen(temp), "%s^%ld", gen_name, i);
                    } else {
                        sprintf(temp + strlen(temp), "%lu*%s^%ld", coeff, gen_name, i);
                    }
                }
                
                // Check if buffer expansion is needed
                size_t temp_len = strlen(temp);
                if (len + temp_len + 2 >= capacity) {
                    capacity = capacity * 2 + temp_len + 100;
                    char *new_buffer = realloc(buffer, capacity);
                    if (!new_buffer) {
                        free(buffer);
                        nmod_poly_clear(poly);
                        return NULL;
                    }
                    buffer = new_buffer;
                }
                
                strcat(buffer, temp);
                len += temp_len;
            }
        }
        
        strcat(buffer, ")");
        
        if (first) {
            strcpy(buffer, "(0)");
        }
        
        nmod_poly_clear(poly);
    }
    
    return buffer;
}

// String builder helper functions
static void sb_init(string_builder_t *sb, size_t initial_capacity) {
    sb->capacity = initial_capacity < 1024 ? 1024 : initial_capacity;
    sb->buffer = (char*) malloc(sb->capacity);
    if (!sb->buffer) {
        sb->capacity = 0;
        sb->length = 0;
        return;
    }
    sb->buffer[0] = '\0';
    sb->length = 0;
}

static int sb_ensure_capacity(string_builder_t *sb, size_t additional) {
    size_t required = sb->length + additional + 1;
    
    if (required > sb->capacity) {
        size_t new_capacity = sb->capacity * 2;
        while (new_capacity < required) {
            new_capacity *= 2;
        }
        
        char *new_buffer = (char*) realloc(sb->buffer, new_capacity);
        if (!new_buffer) {
            return 0;
        }
        
        sb->buffer = new_buffer;
        sb->capacity = new_capacity;
    }
    return 1;
}

static void sb_append(string_builder_t *sb, const char *str) {
    size_t len = strlen(str);
    if (sb_ensure_capacity(sb, len)) {
        memcpy(sb->buffer + sb->length, str, len + 1);
        sb->length += len;
    }
}

static void sb_append_char(string_builder_t *sb, char c) {
    if (sb_ensure_capacity(sb, 1)) {
        sb->buffer[sb->length++] = c;
        sb->buffer[sb->length] = '\0';
    }
}

static void sb_append_long(string_builder_t *sb, long value) {
    char temp[32];
    sprintf(temp, "%ld", value);
    sb_append(sb, temp);
}

static void sb_append_ulong(string_builder_t *sb, unsigned long value) {
    char temp[32];
    sprintf(temp, "%lu", value);
    sb_append(sb, temp);
}

static char* sb_finalize(string_builder_t *sb) {
    char *result = sb->buffer;
    sb->buffer = NULL;
    sb->capacity = 0;
    sb->length = 0;
    return result;
}

// Optimized coefficient to string function
// Updated string builder version for consistent formatting
void fq_nmod_to_string_builder(string_builder_t *sb, const fq_nmod_t elem, 
                                    const fq_nmod_ctx_t ctx, const char *gen_name) {
    if (fq_nmod_is_zero(elem, ctx)) {
        sb_append(sb, "0");
        return;
    }
    
    if (fq_nmod_ctx_degree(ctx) == 1) {
        // Prime field element - no parentheses needed
        nmod_poly_t poly;
        nmod_poly_init(poly, fq_nmod_ctx_prime(ctx));
        fq_nmod_get_nmod_poly(poly, elem, ctx);
        
        if (nmod_poly_degree(poly) >= 0) {
            sb_append_ulong(sb, nmod_poly_get_coeff_ui(poly, 0));
        } else {
            sb_append(sb, "0");
        }
        nmod_poly_clear(poly);
    } else {
        // Extension field element
        nmod_poly_t poly;
        nmod_poly_init(poly, fq_nmod_ctx_prime(ctx));
        fq_nmod_get_nmod_poly(poly, elem, ctx);
        
        slong deg = nmod_poly_degree(poly);
        
        // Count non-zero terms
        int term_count = 0;
        for (slong i = deg; i >= 0; i--) {
            if (nmod_poly_get_coeff_ui(poly, i) != 0) {
                term_count++;
            }
        }
        
        // Use parentheses for multi-term expressions
        int use_parens = (term_count > 1);
        
        if (use_parens) sb_append_char(sb, '(');
        
        int first = 1;
        for (slong i = deg; i >= 0; i--) {
            mp_limb_t coeff = nmod_poly_get_coeff_ui(poly, i);
            if (coeff != 0) {
                if (!first) {
                    sb_append(sb, " + ");
                }
                first = 0;
                
                if (i == 0) {
                    sb_append_ulong(sb, coeff);
                } else if (i == 1) {
                    if (coeff == 1) {
                        sb_append(sb, gen_name ? gen_name : "t");
                    } else {
                        sb_append_ulong(sb, coeff);
                        sb_append_char(sb, '*');
                        sb_append(sb, gen_name ? gen_name : "t");
                    }
                } else {
                    if (coeff == 1) {
                        sb_append(sb, gen_name ? gen_name : "t");
                        sb_append_char(sb, '^');
                        sb_append_long(sb, i);
                    } else {
                        sb_append_ulong(sb, coeff);
                        sb_append_char(sb, '*');
                        sb_append(sb, gen_name ? gen_name : "t");
                        sb_append_char(sb, '^');
                        sb_append_long(sb, i);
                    }
                }
            }
        }
        
        if (first) {
            sb_append(sb, "0");
        }
        
        if (use_parens) sb_append_char(sb, ')');
        
        nmod_poly_clear(poly);
    }
}

// New optimized version function - directly replace original fq_mvpoly_to_string
char* fq_mvpoly_to_string(const fq_mvpoly_t *poly, char **par_names, const char *gen_name) {
    if (poly->nterms == 0) {
        return strdup("0");
    }

    size_t estimated_size = poly->nterms * (50 + 10 * (poly->nvars + poly->npars));
    if (estimated_size < 1024) estimated_size = 1024;

    string_builder_t sb;
    sb_init(&sb, estimated_size);

    for (slong i = 0; i < poly->nterms; i++) {
        if (i > 0) {
            sb_append(&sb, " + ");
        }

        int has_vars_or_pars = 0;

        if (poly->nvars > 0 && poly->terms[i].var_exp) {
            for (slong j = 0; j < poly->nvars; j++) {
                if (poly->terms[i].var_exp[j] > 0) { has_vars_or_pars = 1; break; }
            }
        }
        if (!has_vars_or_pars && poly->npars > 0 && poly->terms[i].par_exp) {
            for (slong j = 0; j < poly->npars; j++) {
                if (poly->terms[i].par_exp[j] > 0) { has_vars_or_pars = 1; break; }
            }
        }

        fq_nmod_t one;
        fq_nmod_init(one, poly->ctx);
        fq_nmod_one(one, poly->ctx);

        if (fq_nmod_is_one(poly->terms[i].coeff, poly->ctx) && has_vars_or_pars) {
        } else {
            fq_nmod_to_string_builder(&sb, poly->terms[i].coeff, poly->ctx, gen_name);
            if (has_vars_or_pars) {
                sb_append_char(&sb, '*');
            }
        }
        fq_nmod_clear(one, poly->ctx);

        int term_has_content = 0;

        if (poly->nvars > 0 && poly->terms[i].var_exp) {
            for (slong j = 0; j < poly->nvars; j++) {
                if (poly->terms[i].var_exp[j] > 0) {
                    if (term_has_content) sb_append_char(&sb, '*');
                    sb_append_char(&sb, 'x');
                    sb_append_long(&sb, j);
                    if (poly->terms[i].var_exp[j] > 1) {
                        sb_append_char(&sb, '^');
                        sb_append_long(&sb, poly->terms[i].var_exp[j]);
                    }
                    term_has_content = 1;
                }
            }
        }

        if (poly->npars > 0 && poly->terms[i].par_exp) {
            for (slong j = 0; j < poly->npars; j++) {
                if (poly->terms[i].par_exp[j] > 0) {
                    if (term_has_content) sb_append_char(&sb, '*');
                    if (par_names && par_names[j]) {
                        sb_append(&sb, par_names[j]);
                    } else {
                        sb_append_char(&sb, 'p');
                        sb_append_long(&sb, j);
                    }
                    if (poly->terms[i].par_exp[j] > 1) {
                        sb_append_char(&sb, '^');
                        sb_append_long(&sb, poly->terms[i].par_exp[j]);
                    }
                    term_has_content = 1;
                }
            }
        }
    }

    return sb_finalize(&sb);
}

// Core Dixon Function

// Internal computation function
char* compute_dixon_internal(const char **poly_strings, slong npoly_strings,
                           const char **var_names, slong nvars,
                           const fq_nmod_ctx_t ctx,
                           char ***remaining_vars, slong *num_remaining) {
    
    if (npoly_strings != nvars + 1) {
        fprintf(stderr, "Error: Need exactly %ld polynomials for %ld variables\n",
                nvars + 1, nvars);
        *remaining_vars = NULL;
        *num_remaining = 0;
        return strdup("0");
    }
    
    // Get generator name
    char *gen_name = get_generator_name(ctx);
    
    // Initialize parser state with original variable names
    parser_state_t state = {0};
    state.var_names = (char**) malloc(nvars * sizeof(char*));
    for (slong i = 0; i < nvars; i++) {
        state.var_names[i] = strdup(var_names[i]);
    }
    state.nvars = nvars;
    state.npars = 0;
    state.max_pars = 16;
    state.par_names = (char**) malloc(state.max_pars * sizeof(char*));
    state.ctx = ctx;
    state.current.str = NULL;
    fq_nmod_init(state.current.value, ctx);
    state.generator_name = strdup(gen_name);
    
    // First pass: identify parameters
    for (slong i = 0; i < npoly_strings; i++) {
        fq_mvpoly_t temp;
        fq_mvpoly_init(&temp, nvars, state.max_pars, ctx);
        
        state.input = poly_strings[i];
        state.pos = 0;
        state.len = strlen(poly_strings[i]);
        next_token(&state);
        
        parse_expression(&state, &temp);
        fq_mvpoly_clear(&temp);
    }
    
    // Save parameters as remaining variables
    *num_remaining = state.npars;
    if (state.npars > 0) {
        *remaining_vars = (char**) malloc(state.npars * sizeof(char*));
        for (slong i = 0; i < state.npars; i++) {
            (*remaining_vars)[i] = strdup(state.par_names[i]);
        }
    } else {
        *remaining_vars = NULL;
    }
    
    // Second pass: parse polynomials
    fq_mvpoly_t *polys = (fq_mvpoly_t*) malloc(npoly_strings * sizeof(fq_mvpoly_t));
    
    for (slong i = 0; i < npoly_strings; i++) {
        fq_mvpoly_init(&polys[i], nvars, state.npars, ctx);
        
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
    // Compute Dixon resultant with original names
    fq_mvpoly_t dixon_result_poly;
    fq_dixon_resultant_with_names(&dixon_result_poly, polys, nvars, state.npars,
                                 state.var_names, state.par_names, gen_name);

    // Find roots with proper parameter names
    find_and_print_roots_of_univariate_resultant(&dixon_result_poly, &state);
    
    // Convert result to string with original parameter names
    char *result_string;
    if (dixon_result_poly.nterms == 0) {
        result_string = strdup("0");
    } else {
        result_string = fq_mvpoly_to_string(&dixon_result_poly, state.par_names, gen_name);
    }
    
    // Cleanup
    fq_mvpoly_clear(&dixon_result_poly);
    for (slong i = 0; i < npoly_strings; i++) {
        fq_mvpoly_clear(&polys[i]);
    }
    free(polys);
    
    for (slong i = 0; i < nvars; i++) {
        free(state.var_names[i]);
    }
    free(state.var_names);
    
    for (slong i = 0; i < state.npars; i++) {
        free(state.par_names[i]);
    }
    free(state.par_names);
    
    if (state.generator_name) {
        free(state.generator_name);
    }
    
    fq_nmod_clear(state.current.value, ctx);
    if (state.current.str) {
        free(state.current.str);
    }
    
    free(gen_name);
    
    return result_string;
}

// Helper function: remove whitespace from both ends of string
static char* trim_whitespace(char *str) {
    if (!str) return NULL;
    
    // Remove leading whitespace
    while (isspace(*str)) str++;
    
    if (*str == '\0') return str;
    
    // Remove trailing whitespace
    char *end = str + strlen(str) - 1;
    while (end > str && isspace(*end)) end--;
    *(end + 1) = '\0';
    
    return str;
}

// Helper function: split string (by comma)
char** split_string(const char *input, slong *count) {
    if (!input || strlen(input) == 0) {
        *count = 0;
        return NULL;
    }
    
    size_t input_len = strlen(input);
    
    char *work_str = (char*) malloc(input_len + 1);
    if (!work_str) {
        *count = 0;
        return NULL;
    }
    memcpy(work_str, input, input_len + 1);
    
    slong num_elements = 1;
    int in_parentheses = 0;
    for (size_t i = 0; i < input_len; i++) {
        if (input[i] == '(') in_parentheses++;
        else if (input[i] == ')') in_parentheses--;
        else if (input[i] == ',' && in_parentheses == 0) {
            num_elements++;
        }
    }
    
    char **result = (char**) malloc(num_elements * sizeof(char*));
    if (!result) {
        free(work_str);
        *count = 0;
        return NULL;
    }
    
    slong idx = 0;
    size_t start = 0;
    in_parentheses = 0;
    
    for (size_t i = 0; i <= input_len; i++) {
        if (i < input_len) {
            if (input[i] == '(') in_parentheses++;
            else if (input[i] == ')') in_parentheses--;
        }
        
        if ((i == input_len || (input[i] == ',' && in_parentheses == 0)) && i > start) {
            size_t poly_len = i - start;
            char *poly = (char*) malloc(poly_len + 1);
            memcpy(poly, input + start, poly_len);
            poly[poly_len] = '\0';
            
            char *trimmed = trim_whitespace(poly);
            result[idx++] = strdup(trimmed);

            free(poly);
            
            start = i + 1;
        }
    }
    
    *count = idx;
    free(work_str);
    
    return result;
}

// Free split string array
void free_split_strings(char **strings, slong count) {
    if (!strings) return;
    for (slong i = 0; i < count; i++) {
        if (strings[i]) free(strings[i]);
    }
    free(strings);
}

typedef struct {
    slong *par_exp;
    fmpz_t residue;
} qq_term_t;

typedef struct {
    qq_term_t *terms;
    slong nterms;
    slong alloc;
    slong npars;
    char **par_names;
    fmpz_t modulus;
} qq_poly_recon_t;

static int qq_same_monomial(const slong *a, const slong *b, slong npars) {
    for (slong i = 0; i < npars; i++) {
        if (a[i] != b[i]) return 0;
    }
    return 1;
}

static void qq_poly_recon_init(qq_poly_recon_t *poly, char **par_names, slong npars) {
    poly->terms = NULL;
    poly->nterms = 0;
    poly->alloc = 0;
    poly->npars = npars;
    poly->par_names = (char**) malloc((size_t) npars * sizeof(char*));
    for (slong i = 0; i < npars; i++) {
        poly->par_names[i] = strdup(par_names[i]);
    }
    fmpz_init(poly->modulus);
    fmpz_one(poly->modulus);
}

static void qq_poly_recon_clear_terms(qq_poly_recon_t *poly) {
    if (!poly || !poly->terms) {
        if (poly) {
            poly->terms = NULL;
            poly->nterms = 0;
            poly->alloc = 0;
        }
        return;
    }

    for (slong i = 0; i < poly->nterms; i++) {
        if (poly->terms[i].par_exp) free(poly->terms[i].par_exp);
        fmpz_clear(poly->terms[i].residue);
    }
    free(poly->terms);
    poly->terms = NULL;
    poly->nterms = 0;
    poly->alloc = 0;
}

static void qq_poly_recon_clear(qq_poly_recon_t *poly) {
    qq_poly_recon_clear_terms(poly);
    if (poly->par_names) {
        for (slong i = 0; i < poly->npars; i++) {
            free(poly->par_names[i]);
        }
        free(poly->par_names);
    }
    fmpz_clear(poly->modulus);
}

static void qq_poly_recon_copy(qq_poly_recon_t *dst, const qq_poly_recon_t *src) {
    if (!dst || !src || dst == src) return;

    qq_poly_recon_clear_terms(dst);
    fmpz_set(dst->modulus, src->modulus);

    if (src->nterms == 0) {
        return;
    }

    dst->terms = (qq_term_t*) malloc((size_t) src->nterms * sizeof(qq_term_t));
    dst->nterms = src->nterms;
    dst->alloc = src->nterms;

    for (slong i = 0; i < src->nterms; i++) {
        dst->terms[i].par_exp = (slong*) calloc((size_t) dst->npars, sizeof(slong));
        if (dst->npars > 0 && src->terms[i].par_exp) {
            memcpy(dst->terms[i].par_exp, src->terms[i].par_exp,
                   (size_t) dst->npars * sizeof(slong));
        }
        fmpz_init(dst->terms[i].residue);
        fmpz_set(dst->terms[i].residue, src->terms[i].residue);
    }
}

static slong qq_poly_recon_find_term(const qq_poly_recon_t *poly, const slong *par_exp) {
    for (slong i = 0; i < poly->nterms; i++) {
        if (qq_same_monomial(poly->terms[i].par_exp, par_exp, poly->npars)) {
            return i;
        }
    }
    return -1;
}

static slong qq_poly_recon_add_term(qq_poly_recon_t *poly, const slong *par_exp) {
    if (poly->nterms >= poly->alloc) {
        poly->alloc = poly->alloc ? 2 * poly->alloc : 16;
        poly->terms = (qq_term_t*) realloc(poly->terms, (size_t) poly->alloc * sizeof(qq_term_t));
    }

    slong idx = poly->nterms++;
    poly->terms[idx].par_exp = (slong*) calloc((size_t) poly->npars, sizeof(slong));
    if (par_exp && poly->npars > 0) {
        memcpy(poly->terms[idx].par_exp, par_exp, (size_t) poly->npars * sizeof(slong));
    }
    fmpz_init(poly->terms[idx].residue);
    fmpz_zero(poly->terms[idx].residue);
    return idx;
}

static ulong fq_nmod_get_ui_prime_field(const fq_nmod_t coeff, const fq_nmod_ctx_t ctx) {
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

static int parse_result_string_fixed_params(const char *result_str,
                                           char **par_names,
                                           slong npars,
                                           const fq_nmod_ctx_t ctx,
                                           fq_mvpoly_t *poly) {
    parser_state_t state = {0};

    state.input = result_str;
    state.pos = 0;
    state.len = strlen(result_str);
    state.var_names = NULL;
    state.nvars = 0;
    state.npars = npars;
    state.max_pars = FLINT_MAX(npars, 16);
    state.par_names = (char**) malloc((size_t) state.max_pars * sizeof(char*));
    for (slong i = 0; i < npars; i++) {
        state.par_names[i] = strdup(par_names[i]);
    }
    state.ctx = ctx;
    state.current.str = NULL;
    fq_nmod_init(state.current.value, ctx);
    state.generator_name = NULL;

    fq_mvpoly_init(poly, 0, npars, ctx);
    next_token(&state);
    parse_expression(&state, poly);

    int ok = (state.npars == npars && state.current.type == TOK_EOF);

    for (slong i = 0; i < state.npars; i++) {
        free(state.par_names[i]);
    }
    free(state.par_names);
    fq_nmod_clear(state.current.value, ctx);
    if (state.current.str) free(state.current.str);

    if (!ok) {
        fq_mvpoly_clear(poly);
        memset(poly, 0, sizeof(*poly));
    }

    return ok;
}

static char **compute_remaining_vars_from_input(const char *poly_string,
                                                const char *vars_string,
                                                slong *num_remaining_out) {
    slong num_polys, num_elim, num_all;
    char **poly_array = split_string(poly_string, &num_polys);
    char **elim_array = split_string(vars_string, &num_elim);
    char **all_vars = NULL;
    char **remaining = NULL;
    slong num_remaining = 0;

    collect_variables((const char **) poly_array, num_polys, NULL, &all_vars, &num_all);

    remaining = (char**) malloc((size_t) num_all * sizeof(char*));
    for (slong i = 0; i < num_all; i++) {
        int is_elim = 0;
        for (slong j = 0; j < num_elim; j++) {
            if (strcmp(all_vars[i], elim_array[j]) == 0) {
                is_elim = 1;
                break;
            }
        }
        if (!is_elim) {
            remaining[num_remaining++] = strdup(all_vars[i]);
        }
    }

    for (slong i = 0; i < num_all; i++) free(all_vars[i]);
    free(all_vars);
    free_split_strings(poly_array, num_polys);
    free_split_strings(elim_array, num_elim);

    *num_remaining_out = num_remaining;
    return remaining;
}

static void qq_poly_recon_absorb_modular_result(qq_poly_recon_t *acc,
                                                const fq_mvpoly_t *mod_poly,
                                                const fq_nmod_ctx_t ctx) {
    ulong prime = fq_nmod_ctx_prime(ctx);
    int *matched = (int*) calloc((size_t) acc->nterms, sizeof(int));
    fmpz_t prev_mod, next_mod;

    fmpz_init(prev_mod);
    fmpz_init(next_mod);
    fmpz_set(prev_mod, acc->modulus);
    fmpz_mul_ui(next_mod, prev_mod, prime);

    for (slong i = 0; i < mod_poly->nterms; i++) {
        ulong residue_ui = fq_nmod_get_ui_prime_field(mod_poly->terms[i].coeff, ctx);
        slong idx = qq_poly_recon_find_term(acc, mod_poly->terms[i].par_exp);
        if (idx < 0) {
            idx = qq_poly_recon_add_term(acc, mod_poly->terms[i].par_exp);
            if (!fmpz_is_one(prev_mod)) {
                fmpz_t tmp_zero;
                fmpz_init(tmp_zero);
                fmpz_zero(tmp_zero);
                fmpz_CRT_ui(acc->terms[idx].residue, tmp_zero, prev_mod, residue_ui, prime, 0);
                fmpz_clear(tmp_zero);
            } else {
                fmpz_set_ui(acc->terms[idx].residue, residue_ui);
            }
            matched = (int*) realloc(matched, (size_t) acc->nterms * sizeof(int));
            matched[idx] = 1;
        } else {
            matched[idx] = 1;
            fmpz_CRT_ui(acc->terms[idx].residue, acc->terms[idx].residue, prev_mod, residue_ui, prime, 0);
        }
    }

    for (slong i = 0; i < acc->nterms; i++) {
        if (!matched[i]) {
            fmpz_CRT_ui(acc->terms[i].residue, acc->terms[i].residue, prev_mod, 0, prime, 0);
        }
    }

    fmpz_set(acc->modulus, next_mod);
    fmpz_clear(prev_mod);
    fmpz_clear(next_mod);
    free(matched);
}

static void append_text(string_builder_t *sb, const char *text) {
    if (!text) return;
    size_t len = strlen(text);
    if (sb->length + len + 1 > sb->capacity) {
        while (sb->length + len + 1 > sb->capacity) {
            sb->capacity = sb->capacity ? 2 * sb->capacity : 128;
        }
        sb->buffer = (char*) realloc(sb->buffer, sb->capacity);
    }
    memcpy(sb->buffer + sb->length, text, len + 1);
    sb->length += len;
}

static void append_char(string_builder_t *sb, char c) {
    if (sb->length + 2 > sb->capacity) {
        sb->capacity = sb->capacity ? 2 * sb->capacity : 128;
        sb->buffer = (char*) realloc(sb->buffer, sb->capacity);
    }
    sb->buffer[sb->length++] = c;
    sb->buffer[sb->length] = '\0';
}

static int qq_poly_recon_to_string(char **out, const qq_poly_recon_t *acc) {
    string_builder_t sb = {NULL, 0, 0};
    int first_term = 1;

    for (slong i = 0; i < acc->nterms; i++) {
        fmpq_t coeff;
        char *coeff_str;
        int has_pars = 0;
        int coeff_sign;
        int coeff_is_pm1;

        fmpq_init(coeff);
        if (!fmpz_is_zero(acc->terms[i].residue)) {
            if (!_fmpq_reconstruct_fmpz(fmpq_numref(coeff), fmpq_denref(coeff),
                                       acc->terms[i].residue, acc->modulus)) {
                fmpq_clear(coeff);
                free(sb.buffer);
                return 0;
            }
            _fmpq_canonicalise(fmpq_numref(coeff), fmpq_denref(coeff));
        } else {
            fmpq_zero(coeff);
        }

        if (fmpq_is_zero(coeff)) {
            fmpq_clear(coeff);
            continue;
        }

        for (slong j = 0; j < acc->npars; j++) {
            if (acc->terms[i].par_exp[j] != 0) {
                has_pars = 1;
                break;
            }
        }

        coeff_sign = fmpq_sgn(coeff);
        coeff_is_pm1 = fmpq_is_pm1(coeff);

        if (!first_term) {
            append_text(&sb, coeff_sign < 0 ? " - " : " + ");
        } else if (coeff_sign < 0) {
            append_char(&sb, '-');
        }

        if (!(has_pars && coeff_is_pm1)) {
            fmpq_t abs_coeff;
            fmpq_init(abs_coeff);
            if (coeff_sign < 0) fmpq_neg(abs_coeff, coeff);
            else fmpq_set(abs_coeff, coeff);
            coeff_str = fmpq_get_str(NULL, 10, abs_coeff);
            append_text(&sb, coeff_str);
            flint_free(coeff_str);
            fmpq_clear(abs_coeff);
            if (has_pars) append_char(&sb, '*');
        }

        if (has_pars) {
            int first_factor = 1;
            for (slong j = 0; j < acc->npars; j++) {
                slong exp = acc->terms[i].par_exp[j];
                if (exp == 0) continue;
                if (!first_factor) append_char(&sb, '*');
                append_text(&sb, acc->par_names[j]);
                if (exp != 1) {
                    char exp_buf[64];
                    snprintf(exp_buf, sizeof(exp_buf), "^%ld", exp);
                    append_text(&sb, exp_buf);
                }
                first_factor = 0;
            }
        }

        first_term = 0;
        fmpq_clear(coeff);
    }

    if (first_term) {
        *out = strdup("0");
        free(sb.buffer);
        return 1;
    }

    *out = sb.buffer;
    return 1;
}

static int qq_poly_recon_reconstruct_coeff(fmpq_t coeff, const qq_poly_recon_t *acc, slong term_idx) {
    if (fmpz_is_zero(acc->terms[term_idx].residue)) {
        fmpq_zero(coeff);
        return 1;
    }

    if (!_fmpq_reconstruct_fmpz(fmpq_numref(coeff), fmpq_denref(coeff),
                                acc->terms[term_idx].residue, acc->modulus)) {
        return 0;
    }

    _fmpq_canonicalise(fmpq_numref(coeff), fmpq_denref(coeff));
    return 1;
}

static int qq_poly_recon_to_mod_prime_string(char **out,
                                             const qq_poly_recon_t *acc,
                                             const fmpz_t prime) {
    string_builder_t sb = {NULL, 0, 0};
    int first_term = 1;
    fmpq_t coeff;
    fmpz_t coeff_mod, num_mod, den_mod, den_inv;

    fmpq_init(coeff);
    fmpz_init(coeff_mod);
    fmpz_init(num_mod);
    fmpz_init(den_mod);
    fmpz_init(den_inv);

    for (slong i = 0; i < acc->nterms; i++) {
        int has_pars = 0;

        if (!qq_poly_recon_reconstruct_coeff(coeff, acc, i)) {
            goto fail;
        }

        if (fmpq_is_zero(coeff)) {
            continue;
        }

        fmpz_mod(num_mod, fmpq_numref(coeff), prime);
        fmpz_mod(den_mod, fmpq_denref(coeff), prime);
        if (fmpz_is_zero(den_mod) || !fmpz_invmod(den_inv, den_mod, prime)) {
            goto fail;
        }

        fmpz_mul(coeff_mod, num_mod, den_inv);
        fmpz_mod(coeff_mod, coeff_mod, prime);
        if (fmpz_is_zero(coeff_mod)) {
            continue;
        }

        for (slong j = 0; j < acc->npars; j++) {
            if (acc->terms[i].par_exp[j] != 0) {
                has_pars = 1;
                break;
            }
        }

        if (!first_term) {
            append_text(&sb, " + ");
        }

        if (!(has_pars && fmpz_is_one(coeff_mod))) {
            char *coeff_str = fmpz_get_str(NULL, 10, coeff_mod);
            append_text(&sb, coeff_str);
            flint_free(coeff_str);
            if (has_pars) append_char(&sb, '*');
        }

        if (has_pars) {
            int first_factor = 1;
            for (slong j = 0; j < acc->npars; j++) {
                slong exp = acc->terms[i].par_exp[j];
                if (exp == 0) continue;
                if (!first_factor) append_char(&sb, '*');
                append_text(&sb, acc->par_names[j]);
                if (exp != 1) {
                    char exp_buf[64];
                    snprintf(exp_buf, sizeof(exp_buf), "^%ld", exp);
                    append_text(&sb, exp_buf);
                }
                first_factor = 0;
            }
        }

        first_term = 0;
    }

    if (first_term) {
        *out = strdup("0");
    } else {
        *out = sb.buffer;
        sb.buffer = NULL;
    }

    free(sb.buffer);
    fmpq_clear(coeff);
    fmpz_clear(coeff_mod);
    fmpz_clear(num_mod);
    fmpz_clear(den_mod);
    fmpz_clear(den_inv);
    return 1;

fail:
    free(sb.buffer);
    fmpq_clear(coeff);
    fmpz_clear(coeff_mod);
    fmpz_clear(num_mod);
    fmpz_clear(den_mod);
    fmpz_clear(den_inv);
    *out = NULL;
    return 0;
}

static void qq_select_reconstruction_primes(ulong *primes, slong *num_primes_out, slong max_primes) {
    ulong candidate;
    slong num_primes = 0;

    candidate = (FLINT_BITS >= 64) ? ((((ulong) 1) << 62) - 4096)
                                   : ((((ulong) 1) << 30) - 4096);
    for (slong i = 0; i < max_primes; i++) {
        candidate = n_nextprime(candidate, 0);
        primes[num_primes++] = candidate;
        candidate += 1024;
    }

    *num_primes_out = num_primes;
}

static char *qq_reconstruct_from_modular_dixon(const char *poly_string,
                                               const char *vars_string,
                                               qq_poly_recon_t *best_recon_out,
                                               int *have_best_recon_out) {
    const slong max_primes = 24;
    const slong min_primes_before_stopping = 12;
    const slong stable_rounds_before_stopping = 8;
    ulong primes[max_primes];
    slong num_primes = 0;
    slong num_remaining;
    char **remaining_vars;
    qq_poly_recon_t recon;
    qq_poly_recon_t best_recon;
    char *best_result = NULL;
    int have_best_recon = 0;
    int stable_count = 0;

    qq_select_reconstruction_primes(primes, &num_primes, max_primes);

    remaining_vars = compute_remaining_vars_from_input(poly_string, vars_string, &num_remaining);
    qq_poly_recon_init(&recon, remaining_vars, num_remaining);
    qq_poly_recon_init(&best_recon, remaining_vars, num_remaining);

    for (slong i = 0; i < num_remaining; i++) free(remaining_vars[i]);
    free(remaining_vars);

    printf("Reconstructing over Q from modular Dixon resultants...\n");

    g_suppress_univariate_root_reporting = 1;
    for (slong i = 0; i < num_primes; i++) {
        fq_nmod_ctx_t ctx;
        fmpz_t p;
        char *mod_result;
        fq_mvpoly_t mod_poly;
        char *candidate_result = NULL;
        int saved_stdout = -1;
        int devnull = -1;

        fmpz_init(p);
        fmpz_set_ui(p, primes[i]);
        fq_nmod_ctx_init(ctx, p, 1, "t");

        printf("Prime %ld/%ld: ", i + 1, num_primes);
        printf("Computing Dixon resultant modulo %lu...\n", primes[i]);

        fflush(stdout);
        saved_stdout = dup(STDOUT_FILENO);
        devnull = open("/dev/null", O_WRONLY);
        if (saved_stdout != -1 && devnull != -1) {
            dup2(devnull, STDOUT_FILENO);
        }

        mod_result = dixon_str(poly_string, vars_string, ctx);

        fflush(stdout);
        if (saved_stdout != -1) {
            dup2(saved_stdout, STDOUT_FILENO);
            close(saved_stdout);
        }
        if (devnull != -1) {
            close(devnull);
        }

        if (!parse_result_string_fixed_params(mod_result, recon.par_names, recon.npars, ctx, &mod_poly)) {
            printf("Skipping p = %lu: failed to parse modular resultant.\n", primes[i]);
            free(mod_result);
            fq_nmod_ctx_clear(ctx);
            fmpz_clear(p);
            continue;
        }

        qq_poly_recon_absorb_modular_result(&recon, &mod_poly, ctx);
        if (qq_poly_recon_to_string(&candidate_result, &recon)) {
            qq_poly_recon_copy(&best_recon, &recon);
            have_best_recon = 1;
            if (best_result && strcmp(best_result, candidate_result) == 0) {
                stable_count++;
                //printf("Reconstruction unchanged after p = %lu.\n", primes[i]);
            } else {
                stable_count = 0;
                free(best_result);
                best_result = strdup(candidate_result);
                printf("Reconstruction updated after p = %lu.\n", primes[i]);
            }
            free(candidate_result);
            if (i + 1 >= min_primes_before_stopping &&
                stable_count >= stable_rounds_before_stopping) {
                printf("Reconstruction stabilized after %ld prime(s).\n", i + 1);
                fq_mvpoly_clear(&mod_poly);
                free(mod_result);
                fq_nmod_ctx_clear(ctx);
                fmpz_clear(p);
                break;
            }
        }

        fq_mvpoly_clear(&mod_poly);
        free(mod_result);
        fq_nmod_ctx_clear(ctx);
        fmpz_clear(p);
    }
    g_suppress_univariate_root_reporting = 0;

    if (!best_result) {
        best_result = strdup("0");
    }

    if (best_recon_out) {
        *best_recon_out = best_recon;
    } else {
        qq_poly_recon_clear(&best_recon);
    }

    if (have_best_recon_out) {
        *have_best_recon_out = have_best_recon;
    }

    qq_poly_recon_clear(&recon);
    return best_result;
}

static void find_and_print_rational_roots_of_univariate_resultant(const qq_poly_recon_t *acc) {
    slong actual_par_count = 0;
    slong main_par_idx = -1;
    int *par_used;

    if (!acc) return;

    par_used = (int*) calloc((size_t) acc->npars, sizeof(int));

    for (slong t = 0; t < acc->nterms; t++) {
        if (acc->terms[t].par_exp) {
            for (slong p = 0; p < acc->npars; p++) {
                if (acc->terms[t].par_exp[p] > 0 && !par_used[p]) {
                    par_used[p] = 1;
                    main_par_idx = p;
                    actual_par_count++;
                }
            }
        }
    }

    if (actual_par_count > 1) {
        free(par_used);
        return;
    }

    if (actual_par_count == 0) {
        free(par_used);
        return;
    }

    fmpq_poly_t rat_poly;
    fmpq_poly_init(rat_poly);

    for (slong i = 0; i < acc->nterms; i++) {
        slong degree = 0;
        fmpq_t coeff;

        if (acc->terms[i].par_exp && acc->terms[i].par_exp[main_par_idx] > 0) {
            degree = acc->terms[i].par_exp[main_par_idx];
        }

        fmpq_init(coeff);
        if (!qq_poly_recon_reconstruct_coeff(coeff, acc, i)) {
            printf("Failed to reconstruct rational coefficient, skipping root search.\n");
            fmpq_clear(coeff);
            fmpq_poly_clear(rat_poly);
            free(par_used);
            return;
        }

        fmpq_poly_set_coeff_fmpq(rat_poly, degree, coeff);
        fmpq_clear(coeff);
    }

    fmpq_poly_canonicalise(rat_poly);

    const char *var_name = acc->par_names[main_par_idx];
    slong degree = fmpq_poly_degree(rat_poly);

    if (degree <= 0) {
        fmpq_poly_clear(rat_poly);
        free(par_used);
        return;
    }

    if (degree > QQ_ROOT_SEARCH_MAX_DEGREE) {
        printf("\nSkipping rational root search in %s: degree %ld exceeds limit %ld.\n",
               var_name, degree, QQ_ROOT_SEARCH_MAX_DEGREE);
        fmpq_poly_clear(rat_poly);
        free(par_used);
        return;
    }

    fmpz_poly_t int_poly, prim_poly, num_divs, den_divs;
    fmpz_t common_den, content, abs_const, abs_lead, coeff, gcd_nd;
    slong zero_mult = 0;
    slong num_roots = 0;
    char **root_strings = NULL;
    slong roots_alloc = 0;

    fmpz_poly_init(int_poly);
    fmpz_poly_init(prim_poly);
    fmpz_poly_init(num_divs);
    fmpz_poly_init(den_divs);
    fmpz_init(common_den);
    fmpz_init(content);
    fmpz_init(abs_const);
    fmpz_init(abs_lead);
    fmpz_init(coeff);
    fmpz_init(gcd_nd);

    fmpq_poly_get_numerator(int_poly, rat_poly);
    fmpq_poly_get_denominator(common_den, rat_poly);

    fmpz_poly_content(content, int_poly);
    if (fmpz_is_zero(content)) {
        goto rational_cleanup;
    }

    fmpz_poly_scalar_divexact_fmpz(prim_poly, int_poly, content);

    fmpz_poly_get_coeff_fmpz(coeff, prim_poly, fmpz_poly_degree(prim_poly));
    if (fmpz_sgn(coeff) < 0) {
        fmpz_poly_neg(prim_poly, prim_poly);
    }

    while (zero_mult <= fmpz_poly_degree(prim_poly)) {
        fmpz_poly_get_coeff_fmpz(coeff, prim_poly, zero_mult);
        if (!fmpz_is_zero(coeff)) break;
        zero_mult++;
    }

    if (zero_mult > 0) {
        if (num_roots >= roots_alloc) {
            roots_alloc = roots_alloc ? 2 * roots_alloc : 4;
            root_strings = (char**) realloc(root_strings, (size_t) roots_alloc * sizeof(char*));
        }
        root_strings[num_roots++] = strdup("0");
    }

    if (zero_mult <= fmpz_poly_degree(prim_poly)) {
        fmpq_t candidate, value;

        fmpq_init(candidate);
        fmpq_init(value);

        fmpz_poly_get_coeff_fmpz(abs_const, prim_poly, zero_mult);
        fmpz_abs(abs_const, abs_const);
        fmpz_poly_get_coeff_fmpz(abs_lead, prim_poly, fmpz_poly_degree(prim_poly));
        fmpz_abs(abs_lead, abs_lead);

        arith_divisors(num_divs, abs_const);
        arith_divisors(den_divs, abs_lead);

        slong candidate_count = 2 * fmpz_poly_length(num_divs) * fmpz_poly_length(den_divs);
        if (candidate_count > QQ_ROOT_SEARCH_MAX_CANDIDATES) {
            printf("\nSkipping rational root search in %s: %ld candidates exceed limit %ld.\n",
                   var_name, candidate_count, QQ_ROOT_SEARCH_MAX_CANDIDATES);
            fmpq_clear(candidate);
            fmpq_clear(value);
            goto rational_cleanup;
        }

        for (slong i = 0; i < fmpz_poly_length(num_divs); i++) {
            fmpz_t num;
            fmpz_init(num);
            fmpz_poly_get_coeff_fmpz(num, num_divs, i);

            for (slong j = 0; j < fmpz_poly_length(den_divs); j++) {
                fmpz_t den;
                fmpz_init(den);
                fmpz_poly_get_coeff_fmpz(den, den_divs, j);

                fmpz_gcd(gcd_nd, num, den);
                if (fmpz_is_one(gcd_nd)) {
                    for (int sign = -1; sign <= 1; sign += 2) {
                        char *root_str;

                        if (sign < 0) fmpz_neg(num, num);
                        fmpq_set_fmpz_frac(candidate, num, den);
                        fmpq_canonicalise(candidate);
                        fmpq_poly_evaluate_fmpq(value, rat_poly, candidate);

                        if (fmpq_is_zero(value)) {
                            if (num_roots >= roots_alloc) {
                                roots_alloc = roots_alloc ? 2 * roots_alloc : 4;
                                root_strings = (char**) realloc(root_strings, (size_t) roots_alloc * sizeof(char*));
                            }
                            root_str = fmpq_get_str(NULL, 10, candidate);
                            root_strings[num_roots++] = strdup(root_str);
                            flint_free(root_str);
                        }
                        if (sign < 0) fmpz_neg(num, num);
                    }
                }

                fmpz_clear(den);
            }

            fmpz_clear(num);
        }

        fmpq_clear(candidate);
        fmpq_clear(value);
    }

    printf("\nRoots in %s (degree %ld):\n", var_name, degree);
    printf("Find %ld roots:\n", num_roots);
    for (slong i = 0; i < num_roots; i++) {
        printf("  Root %ld: %s\n", i + 1, root_strings[i]);
    }

rational_cleanup:
    if (root_strings) {
        for (slong i = 0; i < num_roots; i++) free(root_strings[i]);
        free(root_strings);
    }
    fmpz_poly_clear(int_poly);
    fmpz_poly_clear(prim_poly);
    fmpz_poly_clear(num_divs);
    fmpz_poly_clear(den_divs);
    fmpz_clear(common_den);
    fmpz_clear(content);
    fmpz_clear(abs_const);
    fmpz_clear(abs_lead);
    fmpz_clear(coeff);
    fmpz_clear(gcd_nd);
    fmpq_poly_clear(rat_poly);
    free(par_used);
}

char* dixon_str_rational(const char *poly_string,
                         const char *vars_string) {
    qq_poly_recon_t best_recon;
    char *best_result;
    int have_best_recon = 0;

    best_result = qq_reconstruct_from_modular_dixon(poly_string, vars_string,
                                                    &best_recon, &have_best_recon);
    
    printf("Final reconstruction over Q completed.\n");
    
    if (have_best_recon && strcmp(best_result, "0") != 0) {
        find_and_print_rational_roots_of_univariate_resultant(&best_recon);
    }
    qq_poly_recon_clear(&best_recon);
    return best_result;
}

char* dixon_str_large_prime(const char *poly_string,
                            const char *vars_string,
                            const fmpz_t prime) {
    qq_poly_recon_t best_recon;
    char *best_result;
    char *mod_result = NULL;
    char *prime_str;
    int have_best_recon = 0;

    best_result = qq_reconstruct_from_modular_dixon(poly_string, vars_string,
                                                    &best_recon, &have_best_recon);

    prime_str = fmpz_get_str(NULL, 10, prime);
    printf("Final reconstruction over Q completed.\n");
    printf("Reducing reconstructed resultant modulo p = %s...\n", prime_str);
    flint_free(prime_str);

    if (!have_best_recon || strcmp(best_result, "0") == 0) {
        mod_result = strdup(best_result);
    } else if (!qq_poly_recon_to_mod_prime_string(&mod_result, &best_recon, prime)) {
        printf("Failed to reduce rational reconstruction modulo the target prime.\n");
        mod_result = NULL;
    }

    free(best_result);
    qq_poly_recon_clear(&best_recon);
    return mod_result;
}

// Main Dixon Interface

// Unified Dixon function - accepts array of polynomial strings
// Returns result as string, optionally outputs remaining variables
char* dixon(const char **poly_strings, slong num_polys, 
            const char **elim_vars, slong num_elim_vars,
            const fq_nmod_ctx_t ctx) {

    clock_t start = clock();
    /*
    for (slong i = 0; i < num_polys; i++) {
        printf("  p%ld: %s\n", i, poly_strings[i]);
    }
    */
    // Compute Dixon resultant
    char **remaining_vars = NULL;
    slong num_remaining = 0;
    
    char *result = compute_dixon_internal(poly_strings, num_polys, 
                                         elim_vars, num_elim_vars, ctx,
                                         &remaining_vars, &num_remaining);
    
    // Cleanup remaining vars
    if (remaining_vars) {
        for (slong i = 0; i < num_remaining; i++) {
            free(remaining_vars[i]);
        }
        free(remaining_vars);
    }

    clock_t end = clock();
    (void) end;
    (void) start;

    return result;
}

// Implementation using unified_mpoly_resultant
char* bivariate_resultant(const char *poly1_str, const char *poly2_str,
                                         const char *elim_var, const fq_nmod_ctx_t ctx) {
    clock_t total_start = clock();
    // Get generator name
    char *gen_name = get_generator_name(ctx);
    printf("Use FLINT's built-in resultant to compute two polynomials\n");
    char **remaining_vars = NULL;
    slong num_remaining = 0;
    
    // Initialize parser state
    parser_state_t state = {0};
    state.var_names = (char**) malloc(1 * sizeof(char*));
    state.var_names[0] = strdup(elim_var);
    state.nvars = 1;
    state.npars = 0;
    state.max_pars = 16;
    state.par_names = (char**) malloc(state.max_pars * sizeof(char*));
    state.ctx = ctx;
    state.current.str = NULL;
    fq_nmod_init(state.current.value, ctx);
    state.generator_name = strdup(gen_name);
    
    // First pass: parse to identify parameters
    fq_mvpoly_t temp1, temp2;
    fq_mvpoly_init(&temp1, 1, state.max_pars, ctx);
    fq_mvpoly_init(&temp2, 1, state.max_pars, ctx);

    state.input = poly1_str;
    state.pos = 0;
    state.len = strlen(poly1_str);
    next_token(&state);
    parse_expression(&state, &temp1);

    state.input = poly2_str;
    state.pos = 0;
    state.len = strlen(poly2_str);
    if (state.current.str) {
        free(state.current.str);
        state.current.str = NULL;
    }
    next_token(&state);
    parse_expression(&state, &temp2);
    
    fq_mvpoly_clear(&temp1);
    fq_mvpoly_clear(&temp2);
    
    // Save parameters
    num_remaining = state.npars;
    if (state.npars > 0) {
        remaining_vars = (char**) malloc(state.npars * sizeof(char*));
        for (slong i = 0; i < state.npars; i++) {
            (remaining_vars)[i] = strdup(state.par_names[i]);
        }
    } else {
        remaining_vars = NULL;
    }

    // Second pass: formal parsing
    fq_mvpoly_t poly1, poly2;
    fq_mvpoly_init(&poly1, 1, state.npars, ctx);
    fq_mvpoly_init(&poly2, 1, state.npars, ctx);
    
    state.input = poly1_str;
    state.pos = 0;
    state.len = strlen(poly1_str);
    if (state.current.str) {
        free(state.current.str);
        state.current.str = NULL;
    }
    next_token(&state);
    parse_expression(&state, &poly1);
    
    state.input = poly2_str;
    state.pos = 0;
    state.len = strlen(poly2_str);
    if (state.current.str) {
        free(state.current.str);
        state.current.str = NULL;
    }
    next_token(&state);
    parse_expression(&state, &poly2);
    
    // Initialize unified field context
    field_ctx_t field_ctx;
    field_ctx_init(&field_ctx, ctx);
    
    // Create unified multivariate polynomial context
    slong total_vars = 1 + state.npars;  // elimination variable + parameters
    unified_mpoly_ctx_t unified_ctx = unified_mpoly_ctx_init(total_vars, ORD_LEX, &field_ctx);
    
    // Initialize unified polynomials
    unified_mpoly_t A = unified_mpoly_init(unified_ctx);
    unified_mpoly_t B = unified_mpoly_init(unified_ctx);
    unified_mpoly_t R = unified_mpoly_init(unified_ctx);

    // Convert first polynomial
    for (slong i = 0; i < poly1.nterms; i++) {
        field_elem_u coeff;
        ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
        
        // Convert coefficient
        fq_nmod_to_field_elem(&coeff, poly1.terms[i].coeff, &field_ctx);
        
        // Set exponents
        if (poly1.terms[i].var_exp) {
            exp[0] = poly1.terms[i].var_exp[0];
        }
        if (poly1.terms[i].par_exp) {
            for (slong j = 0; j < state.npars; j++) {
                exp[1 + j] = poly1.terms[i].par_exp[j];
            }
        }
        
        unified_mpoly_set_coeff_ui(A, &coeff, exp);
        free(exp);
    }

    // Convert second polynomial
    for (slong i = 0; i < poly2.nterms; i++) {
        field_elem_u coeff;
        ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
        
        // Convert coefficient
        fq_nmod_to_field_elem(&coeff, poly2.terms[i].coeff, &field_ctx);
        
        // Set exponents
        if (poly2.terms[i].var_exp) {
            exp[0] = poly2.terms[i].var_exp[0];
        }
        if (poly2.terms[i].par_exp) {
            for (slong j = 0; j < state.npars; j++) {
                exp[1 + j] = poly2.terms[i].par_exp[j];
            }
        }
        
        unified_mpoly_set_coeff_ui(B, &coeff, exp);
        free(exp);
    }
    
    // Compute resultant
    clock_t start = clock();
    int success = unified_mpoly_resultant(R, A, B, 0, unified_ctx);
    clock_t end = clock();
    
    if (!success) {
        printf("Unified resultant computation failed!\n");
        unified_mpoly_clear(A);
        unified_mpoly_clear(B);
        unified_mpoly_clear(R);
        unified_mpoly_ctx_clear(unified_ctx);
        
        // Clean up other resources
        for (slong i = 0; i < state.nvars; i++) {
            free(state.var_names[i]);
        }
        free(state.var_names);
        for (slong i = 0; i < state.npars; i++) {
            free(state.par_names[i]);
        }
        free(state.par_names);
        if (state.generator_name) free(state.generator_name);
        fq_nmod_clear(state.current.value, ctx);
        if (state.current.str) free(state.current.str);
        free(gen_name);
        fq_mvpoly_clear(&poly1);
        fq_mvpoly_clear(&poly2);
        
        return strdup("0");
    }
    
    // Convert result back to fq_mvpoly format
    fq_mvpoly_t result_mvpoly;
    fq_mvpoly_init(&result_mvpoly, 0, state.npars, ctx);
    
    // Convert from unified format back
    slong result_len = unified_mpoly_length(R);
    if (result_len > 0) {
        // Need to iterate through all terms of the result
        // This requires accessing the internal structure of unified_mpoly
        // Since unified_mpoly is a wrapper for nmod_mpoly or fq_nmod_mpoly
        // We need to handle based on field_id
        
        if (field_ctx.field_id == FIELD_ID_NMOD) {
            // Handle nmod case
            nmod_mpoly_struct *nmod_res = GET_NMOD_POLY(R);
            nmod_mpoly_ctx_struct *nmod_ctx = &(unified_ctx->ctx.nmod_ctx);
            
            for (slong i = 0; i < nmod_mpoly_length(nmod_res, nmod_ctx); i++) {
                ulong coeff_ui = nmod_mpoly_get_term_coeff_ui(nmod_res, i, nmod_ctx);
                ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
                nmod_mpoly_get_term_exp_ui(exp, nmod_res, i, nmod_ctx);
                
                // Convert coefficient
                fq_nmod_t coeff;
                fq_nmod_init(coeff, ctx);
                fq_nmod_set_ui(coeff, coeff_ui, ctx);
                
                // Extract parameter exponents
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
        } else if (field_ctx.field_id == FIELD_ID_FQ_ZECH) {
        // Handle fq_zech case
        fq_zech_mpoly_struct *zech_res = GET_ZECH_POLY(R);
        fq_zech_mpoly_ctx_struct *zech_ctx = &(unified_ctx->ctx.zech_ctx);
        
        for (slong i = 0; i < fq_zech_mpoly_length(zech_res, zech_ctx); i++) {
            fq_zech_t zech_coeff;
            fq_zech_init(zech_coeff, zech_ctx->fqctx);
            fq_zech_mpoly_get_term_coeff_fq_zech(zech_coeff, zech_res, i, zech_ctx);
            
            // Convert to fq_nmod
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_zech_get_fq_nmod(coeff, zech_coeff, zech_ctx->fqctx);
            
            ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
            fq_zech_mpoly_get_term_exp_ui(exp, zech_res, i, zech_ctx);
            
            // Extract parameter exponents
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
            // Handle fq_nmod case
            fq_nmod_mpoly_struct *fq_res = GET_FQ_POLY(R);
            fq_nmod_mpoly_ctx_struct *fq_ctx = &(unified_ctx->ctx.fq_ctx);
            for (slong i = 0; i < fq_nmod_mpoly_length(fq_res, fq_ctx); i++) {
                fq_nmod_t coeff;
                fq_nmod_init(coeff, ctx);
                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, fq_res, i, fq_ctx);
                
                ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
                fq_nmod_mpoly_get_term_exp_ui(exp, fq_res, i, fq_ctx);
                
                // Extract parameter exponents
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
    }
    if (g_field_equation_reduction) {
        fq_mvpoly_reduce_field_equation(&result_mvpoly);
    }
    fq_mvpoly_make_monic(&result_mvpoly);

    // If it's a univariate polynomial, try to find roots
    find_and_print_roots_of_univariate_resultant(&result_mvpoly, &state);
    //printf("fq_mvpoly_to_string\n");
    // Convert result to string
    char *result_string;
    if (result_mvpoly.nterms == 0) {
        result_string = strdup("0");
    } else {
        result_string = fq_mvpoly_to_string(&result_mvpoly, state.par_names, gen_name);
    }
    //printf("%s", result_string);
    //printf("clean up\n");

    // Output remaining variable info
    if (num_remaining > 0) {
        for (slong i = 0; i < num_remaining; i++) {
            free(remaining_vars[i]);
        }
        free(remaining_vars);
    }
    // Clean up
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
    if (state.generator_name) free(state.generator_name);
    fq_nmod_clear(state.current.value, ctx);
    if (state.current.str) free(state.current.str);
    free(gen_name);

    clock_t total_end = clock();
    (void) total_end;
    (void) total_start;

    return result_string;
}

char* dixon_str(const char *poly_string,
                const char *vars_string,
                const fq_nmod_ctx_t ctx) {

    slong num_polys, num_vars;
    char **poly_array = split_string(poly_string, &num_polys);
    char **vars_array = split_string(vars_string, &num_vars);
    
    char *result = NULL;
    
    if (num_polys == 2 && num_vars == 1) {
        result = bivariate_resultant(poly_array[0], poly_array[1], 
                                     vars_array[0], ctx);        
    } else {
        const char **poly_strings = (const char**) malloc(num_polys * sizeof(char*));
        const char **elim_vars   = (const char**) malloc(num_vars  * sizeof(char*));
        
        for (slong i = 0; i < num_polys; i++) poly_strings[i] = poly_array[i];
        for (slong i = 0; i < num_vars;  i++) elim_vars[i]    = vars_array[i];
        
        result = dixon(poly_strings, num_polys, elim_vars, num_vars, ctx);
        
        free(poly_strings);
        free(elim_vars);
    }
    
    free_split_strings(poly_array, num_polys);
    free_split_strings(vars_array, num_vars);
    
    return result;
}
