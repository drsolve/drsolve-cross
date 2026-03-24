#ifndef POLYNOMIAL_SYSTEM_SOLVER_H
#define POLYNOMIAL_SYSTEM_SOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// Include the existing headers
#include "dixon_flint.h"
#include "fq_mvpoly.h"
#include "dixon_interface_flint.h"
#include "dixon_complexity.h"

// Forward declarations
typedef struct polynomial_solutions polynomial_solutions_t;
typedef struct variable_info variable_info_t;
typedef struct equation_combination equation_combination_t;

// ============= DATA STRUCTURES =============

// Enhanced solution structure to handle multiple complete solution sets
typedef struct polynomial_solutions {
    char **variable_names;      // Array of variable names
    slong num_variables;        // Number of variables
    fq_nmod_t ***solution_sets; // 3D array: [solution_set][variable][solution_index]
    slong *solutions_per_var;   // Number of solutions for each variable in each set
    slong num_solution_sets;    // Number of complete solution sets
    slong num_equations;        // Number of input equations
    int is_valid;               // Whether the solution is valid
    int has_no_solutions;       // 1 = no solutions, 0 = has solutions, -1 = dimension > 0
    char *error_message;        // Error message
    char *variable_order;       // Variable order used by solver
    char *elimination_summary;  // Summary of elimination strategy
    slong total_combinations;   // Total equation combinations tried
    slong successful_combinations; // Non-zero resultant combinations
    slong num_base_solutions;   // Number of roots of final univariate polynomial
    slong checked_solution_sets;    // Candidate sets checked by verification
    slong verified_solution_sets;   // Candidate sets that passed verification
    char **resultant_steps;     // Logged resultant/elimination steps
    slong num_resultant_steps;  // Number of logged steps
    slong resultant_steps_cap;  // Capacity of logged steps
    char **candidate_solution_lines; // Human-readable candidate solutions
    int *candidate_solution_pass;    // Verification result for each candidate
    slong num_candidate_solution_lines; // Number of candidate solution lines
    slong candidate_solution_lines_cap; // Capacity of candidate solution lines
    const fq_nmod_ctx_struct *ctx;  // Finite field context
} polynomial_solutions_t;

// Variable information structure
typedef struct variable_info {
    char *name;
    slong max_degree;
    slong index;
} variable_info_t;

// Equation combination structure for tracking different equation selections
typedef struct equation_combination {
    slong *equation_indices;    // Indices of selected equations
    slong num_equations;        // Number of equations in this combination
    int tried;                  // Whether this combination has been tried
    int success;                // Whether this combination succeeded
} equation_combination_t;

// ============= FUNCTION DECLARATIONS =============

// Basic utility functions
int contains_variable(const char *poly_str, const char *var_name);
char* get_generator_name_for_solver(const fq_nmod_ctx_t ctx);
int contains_identifier_improved(const char *poly_str, const char *identifier, const fq_nmod_ctx_t ctx);
char** extract_variables_improved(char **poly_strings, slong num_polys, slong *num_vars_out, const fq_nmod_ctx_t ctx);
slong get_variable_max_degree_in_polynomial(const char *poly_str, const char *var_name, const fq_nmod_ctx_t ctx);

// Equation combination management
void generate_equation_combinations(slong num_equations, slong target_equations, 
                                   equation_combination_t **combinations, slong *num_combinations);
void free_equation_combinations(equation_combination_t *combinations, slong num_combinations);

// Solution structure management
void polynomial_solver_set_realtime_progress(int enabled);
void polynomial_solutions_init(polynomial_solutions_t *sols, slong num_vars, 
                               const fq_nmod_ctx_t ctx);
void polynomial_solutions_clear(polynomial_solutions_t *sols);
void print_polynomial_solutions(const polynomial_solutions_t *sols);

// Variable analysis functions
int extract_and_sort_variables(char **poly_strings, slong num_polys,
                              variable_info_t **sorted_vars_out, slong *num_vars_out,
                              const fq_nmod_ctx_t ctx);

// Polynomial manipulation functions
fq_nmod_t* solve_univariate_equation_all_roots(const char *poly_str, const char *var_name,
                                               slong *num_roots_out, const fq_nmod_ctx_t ctx);
char* substitute_variable_in_polynomial(const char *poly_str, const char *var_name,
                                        const fq_nmod_t value, const fq_nmod_ctx_t ctx);
char* eliminate_variable_dixon_with_selection(char **poly_strings, slong num_polys, 
                                              const char *elim_var, const fq_nmod_ctx_t ctx,
                                              equation_combination_t *combination);

// Solution verification
int verify_solution_set(char **original_polys, slong num_polys,
                       variable_info_t *sorted_vars, slong num_vars,
                       fq_nmod_t **solution_values, const fq_nmod_ctx_t ctx);
void filter_solutions_by_verification(polynomial_solutions_t *sols,
                                     char **original_polys, slong num_polys,
                                     variable_info_t *sorted_vars);

// Recursive back substitution
int solve_by_back_substitution_recursive_enhanced(char **original_polys, slong num_polys,
                                                  variable_info_t *sorted_vars, slong num_vars,
                                                  fq_nmod_t *base_solutions, slong num_base_solutions,
                                                  polynomial_solutions_t *sols);

// Main solving algorithms
int solve_by_elimination_enhanced(char **poly_strings, slong num_polys,
                                 variable_info_t *sorted_vars, slong num_vars,
                                 polynomial_solutions_t *sols);

// Main solver interfaces
polynomial_solutions_t* solve_polynomial_system_array_with_vars(char **poly_strings, slong num_polys,
                                                               variable_info_t *original_vars, slong num_original_vars,
                                                               const fq_nmod_ctx_t ctx);
polynomial_solutions_t* solve_polynomial_system_array(char **poly_strings, slong num_polys,
                                                      const fq_nmod_ctx_t ctx);
polynomial_solutions_t* solve_polynomial_system_string(const char *poly_string, 
                                                       const fq_nmod_ctx_t ctx);

// Test function
void test_polynomial_solver(void);

#endif // POLYNOMIAL_SYSTEM_SOLVER_H
