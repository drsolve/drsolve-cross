#ifndef DIXON_WITH_IDEAL_REDUCTION_H
#define DIXON_WITH_IDEAL_REDUCTION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_mpoly.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mpoly.h>
#include "dixon_flint.h"

/* Debug output control macro - set to 0 to disable all output */
#define DEBUG_OUTPUT_R 0

#if DEBUG_OUTPUT_R
    #define DEBUG_PRINT_R(...) printf(__VA_ARGS__)
#else
    #define DEBUG_PRINT_R(...) ((void)0)
#endif

/* Unified triangular ideal structure supporting both prime and extension fields */
typedef struct {
    void **generators;            /* Array of polynomial generators */
    slong num_gens;              /* Number of generators */
    slong *var_indices;          /* Which variable each generator eliminates */
    char **var_names;            /* Names of variables being reduced */
    slong *leading_degrees;      /* Leading degree of each generator for optimization */
    slong max_gens;              /* Maximum number of generators */
    int is_prime_field;          /* 1 for prime field, 0 for extension field */
    union {
        nmod_mpoly_ctx_struct *nmod_ctx;
        fq_nmod_mpoly_ctx_struct *fq_ctx;
    } ctx;
    const fq_nmod_ctx_struct *field_ctx;
    char **complete_var_names;   /* Complete list of variable names */
    slong total_vars;            /* Total number of variables */
} unified_triangular_ideal_t;

/* Equation information structure for parsing equation format generators */
typedef struct {
    char *main_var_name;      /* Main variable name */
    slong main_var_degree;    /* Main variable degree */
    fq_mvpoly_t lhs;         /* Left-hand side expression */
    fq_mvpoly_t rhs;         /* Right-hand side expression */
    fq_mvpoly_t standard;    /* Standard form: lhs - rhs */
} equation_info_t;

/* Core triangular ideal functions */
void unified_triangular_ideal_init(unified_triangular_ideal_t *ideal, 
                                  slong max_gens, slong nvars, 
                                  const fq_nmod_ctx_t field_ctx);

void unified_triangular_ideal_clear(unified_triangular_ideal_t *ideal);

void create_reduced_ideal_context(unified_triangular_ideal_t *reduced_ideal,
                                const unified_triangular_ideal_t *original_ideal,
                                char **current_var_names,
                                slong current_nvars,
                                const fq_nmod_ctx_t field_ctx);

/* Polynomial reduction functions with ideal */
void triangular_ideal_reduce_nmod_mpoly_with_names(nmod_mpoly_t poly, 
                                                   const unified_triangular_ideal_t *ideal,
                                                   char **current_var_names);

void triangular_ideal_reduce_fq_nmod_mpoly_with_names(fq_nmod_mpoly_t poly,
                                                      const unified_triangular_ideal_t *ideal,
                                                      char **current_var_names);

void triangular_ideal_reduce_nmod_mpoly(nmod_mpoly_t poly, 
                                       const unified_triangular_ideal_t *ideal);

void triangular_ideal_reduce_fq_nmod_mpoly(fq_nmod_mpoly_t poly,
                                          const unified_triangular_ideal_t *ideal);

/* Matrix determinant computation with ideal reduction */
void compute_nmod_det_with_triangular_reduction_with_names(nmod_mpoly_t det,
                                                          nmod_mpoly_struct **matrix,
                                                          slong size,
                                                          const unified_triangular_ideal_t *ideal,
                                                          char **current_var_names);

void compute_fq_nmod_det_with_triangular_reduction_with_names(fq_nmod_mpoly_t det,
                                                             fq_nmod_mpoly_struct **matrix,
                                                             slong size,
                                                             const unified_triangular_ideal_t *ideal,
                                                             char **current_var_names);

void compute_nmod_det_with_triangular_reduction(nmod_mpoly_t det,
                                               nmod_mpoly_struct **matrix,
                                               slong size,
                                               const unified_triangular_ideal_t *ideal);

void compute_fq_nmod_det_with_triangular_reduction(fq_nmod_mpoly_t det,
                                                  fq_nmod_mpoly_struct **matrix,
                                                  slong size,
                                                  const unified_triangular_ideal_t *ideal);

void compute_det_with_reduction_from_mvpoly(fq_mvpoly_t *result,
                                           fq_mvpoly_t **matrix,
                                           slong size,
                                           const unified_triangular_ideal_t *ideal,
                                           char **current_var_names);

/* Equation information functions */
void equation_info_init(equation_info_t *eq, slong nvars, slong npars, const fq_nmod_ctx_t ctx);

void equation_info_clear(equation_info_t *eq);

/* Ideal construction from equation format strings */
void unified_triangular_ideal_add_generator_equation_format(unified_triangular_ideal_t *ideal,
                                                           const char *equation_str,
                                                           char **var_names_hint);

void construct_triangular_ideal_from_strings(unified_triangular_ideal_t *ideal,
                                              const char **equation_strs,
                                              slong num_equations,
                                              const char **var_names,
                                              slong nvars,
                                              const fq_nmod_ctx_t ctx);

void construct_triangular_ideal_str(unified_triangular_ideal_t *ideal,
                                   const char *gens_string,   
                                   const char *vars_string,   
                                   const fq_nmod_ctx_t ctx);

/* Main Dixon resultant computation with ideal reduction */
char* dixon_with_ideal_reduction(const char **poly_strings, slong num_polys,
                                const char **elim_vars, slong num_elim_vars,
                                const fq_nmod_ctx_t ctx,
                                unified_triangular_ideal_t *ideal);

char* dixon_with_ideal(const char **poly_strings,
                           slong num_polys,
                           const char **elim_vars,
                           slong num_elim_vars,
                           const char *ideal_string,
                           const fq_nmod_ctx_t ctx);

char* dixon_with_ideal_reduction_str(const char *poly_string,
                                    const char *elim_vars_string,
                                    const char *ideal_gens_string,
                                    const fq_nmod_ctx_t ctx);

/* String manipulation utilities */
char** split_string_r(const char* str, slong* count);

void free_split_string_rs(char** strings, slong count);

char** extract_all_variables_from_ideal_gens(const char **ideal_gens_strs, 
                                            slong num_gens,
                                            const fq_nmod_ctx_t ctx,
                                            slong *total_vars);

/* Test functions */
void test_iterative_elimination(void);

void test_iterative_elimination_str(void);

#endif /* DIXON_WITH_IDEAL_REDUCTION_H */