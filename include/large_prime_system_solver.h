#ifndef LARGE_PRIME_SYSTEM_SOLVER_H
#define LARGE_PRIME_SYSTEM_SOLVER_H

#include <stdio.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    slong num_variables;
    slong num_solution_sets;
    slong num_equations;
    char **variable_names;
    fmpz_t **solution_values;
    int is_valid;
    int has_no_solutions;
    char *error_message;
    fmpz_t prime;
    fmpz_mod_ctx_t ctx;
} large_prime_solutions_t;

void large_prime_solver_set_realtime_progress(int enabled);
large_prime_solutions_t* solve_large_prime_polynomial_system_string(const char *poly_string,
                                                                   const fmpz_t prime);
void large_prime_solutions_clear(large_prime_solutions_t *sols);
void print_large_prime_solutions(const large_prime_solutions_t *sols);
void save_large_prime_solver_result_to_file(const char *filename,
                                            const char *polys_str,
                                            const fmpz_t prime,
                                            const large_prime_solutions_t *sols,
                                            double cpu_time,
                                            double wall_time,
                                            int threads_num);
void large_prime_print_roots_from_resultant_string(const char *resultant_str,
                                                   const char *polys_str,
                                                   const char *vars_str,
                                                   const fmpz_t prime,
                                                   FILE *fp,
                                                   int print_to_stdout);

#ifdef __cplusplus
}
#endif

#endif