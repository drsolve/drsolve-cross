// Complete fixed Dixon resultant string interface header
#ifndef DIXON_INTERFACE_FLINT_H
#define DIXON_INTERFACE_FLINT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <flint/flint.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>
#include <flint/fq_nmod_mat.h>
#include <flint/fmpz.h>
#include <flint/fq_nmod_poly_factor.h>
#include "unified_mpoly_resultant.h"
#include "fq_nmod_roots.h"
#include "fmpq_acb_roots.h"
#include "fq_mvpoly.h"
#include "dixon_flint.h"
// Debug switch
#define DEBUG_PARSER 0

#if DEBUG_PARSER
#define DEBUG_PRINT(fmt, ...) printf("[PARSER] " fmt, ##__VA_ARGS__)
#else
#define DEBUG_PRINT(fmt, ...)
#endif

// Dixon resultant computation function
void fq_dixon_resultant(fq_mvpoly_t *result, fq_mvpoly_t *polys, 
                       slong nvars, slong npars);

// Token type definitions
typedef enum {
    TOK_NUMBER,      // number
    TOK_VARIABLE,    // variable (x, y, z etc)
    TOK_PARAMETER,   // parameter (a, b, c etc)
    TOK_GENERATOR,   // field generator (t)
    TOK_PLUS,        // +
    TOK_MINUS,       // -
    TOK_MULT,        // *
    TOK_POWER,       // ^
    TOK_LPAREN,      // (
    TOK_RPAREN,      // )
    TOK_EOF          // end
} token_type_t;

typedef struct {
    token_type_t type;
    char *str;
    fq_nmod_t value;
    slong int_value;
    int int_value_valid;
    const fq_nmod_ctx_struct *ctx;
} token_t;

typedef struct {
    const char *input;
    size_t pos;
    size_t len;
    token_t current;
    
    char **var_names;
    slong nvars;
    char **par_names;
    slong npars;
    slong max_pars;
    
    const fq_nmod_ctx_struct *ctx;
    char *generator_name;
} parser_state_t;

// String builder structure for efficient string construction
typedef struct {
    char *buffer;
    size_t length;      // current string length
    size_t capacity;    // total allocated capacity
} string_builder_t;

void next_token(parser_state_t *state);
void parse_expression(parser_state_t *state, fq_mvpoly_t *poly);
void parse_term(parser_state_t *state, fq_mvpoly_t *poly);
void parse_factor(parser_state_t *state, fq_mvpoly_t *poly);
void parse_primary(parser_state_t *state, fq_mvpoly_t *poly);
char** split_string(const char *input, slong *count);
void free_split_strings(char **strings, slong count);

// Helper functions
char* get_generator_name(const fq_nmod_ctx_t ctx);
char* fq_nmod_to_string_with_gen(const fq_nmod_t elem, const fq_nmod_ctx_t ctx, const char *gen_name);
char* fq_mvpoly_to_string(const fq_mvpoly_t *poly, char **var_names, const char *gen_name);

// Output functions
void fq_nmod_print_pretty_enhanced(const fq_nmod_t a, const fq_nmod_ctx_t ctx);
void fq_mvpoly_print_enhanced(const fq_mvpoly_t *p, const char *name);
void find_and_print_roots_of_univariate_resultant(const fq_mvpoly_t *result, parser_state_t *state);
void find_and_print_roots_of_univariate_resultant_with_file(const fq_mvpoly_t *result, parser_state_t *state, FILE *fp_file, int print_to_stdout);

// Main Dixon interface functions
char* dixon(const char **poly_strings, slong num_polys, 
            const char **elim_vars, slong num_elim_vars,
            const fq_nmod_ctx_t ctx);

char* bivariate_resultant(const char *poly1_str, const char *poly2_str,
                         const char *elim_var, const fq_nmod_ctx_t ctx);

char* dixon_str(const char *poly_string,    // comma-separated polynomials
                const char *vars_string,     // comma-separated variables
                const fq_nmod_ctx_t ctx);

char* dixon_str_with_file(const char *poly_string,    // comma-separated polynomials
                          const char *vars_string,     // comma-separated variables
                          const fq_nmod_ctx_t ctx,
                          FILE *fp_file,
                          int print_to_stdout);

char* dixon_str_rational_with_file(const char *poly_string,
                                    const char *vars_string,
                                    const char *output_filename);

char* dixon_str_rational(const char *poly_string,
                         const char *vars_string);

char* dixon_str_large_prime(const char *poly_string,
                            const char *vars_string,
                            const fmpz_t prime);

void append_roots_to_file_from_result(const char *result_str,
                                       const char *polys_str,
                                       const char *vars_str,
                                       const fq_nmod_ctx_t ctx,
                                       FILE *fp_file);

#endif // DIXON_INTERFACE_FLINT_H
