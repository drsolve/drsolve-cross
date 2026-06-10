#ifndef COMPLEX_SOLVER_H
#define COMPLEX_SOLVER_H

#include <stdio.h>

#include "fmpq_acb_roots.h"
#include "flint/fmpq_poly.h"

/*
 * Initial complex-solver module.
 *
 * The first target is Q-univariate complex roots and a small self-test hook.
 * The API is intentionally small so the later 2x2 complex backsolve path can
 * extend this module without coupling it to rational_system_solver.c internals.
 */

typedef struct {
    acb_t *roots;
    slong *multiplicities;
    slong num_roots;
    slong alloc;
    slong prec;
    char *source_poly;
    char *variable_name;
} complex_solver_roots_t;

typedef struct {
    acb_t x;
    acb_t y;
    acb_t residual1;
    acb_t residual2;
} complex_solver_solution_pair_t;

typedef struct {
    complex_solver_solution_pair_t *items;
    slong num_solutions;
    slong alloc;
    slong prec;
    char *source_system;
    char *x_name;
    char *y_name;
} complex_solver_solution_list_t;

void complex_solver_roots_init(complex_solver_roots_t *roots);
void complex_solver_roots_clear(complex_solver_roots_t *roots);
void complex_solver_roots_print(const complex_solver_roots_t *roots, FILE *fp, slong digits);

void complex_solver_solution_list_init(complex_solver_solution_list_t *sols);
void complex_solver_solution_list_clear(complex_solver_solution_list_t *sols);
void complex_solver_solution_list_print(const complex_solver_solution_list_t *sols, FILE *fp, slong digits);

int complex_solver_parse_univariate_fmpq_poly(const char *poly_str,
                                              const char *var_name,
                                              fmpq_poly_t poly);

int complex_solver_univariate_complex_roots_from_string(const char *poly_str,
                                                        const char *var_name,
                                                        slong prec,
                                                        complex_solver_roots_t *roots_out);

int complex_solver_solve_bivariate_2x2_from_string(const char *poly_string,
                                                   slong prec,
                                                   complex_solver_solution_list_t *sols_out);

#endif
