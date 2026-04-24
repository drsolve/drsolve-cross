#ifndef RATIONAL_SYSTEM_SOLVER_H
#define RATIONAL_SYSTEM_SOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "fmpq_acb_roots.h"
#include "dixon_interface_flint.h"
#include "flint/fmpq.h"
#include "flint/fmpq_poly.h"
#include "flint/fmpq_mpoly.h"

// Forward declarations
typedef struct rational_solutions rational_solutions_t;
typedef struct rational_variable_info rational_variable_info_t;
typedef struct rational_equation_combination rational_equation_combination_t;

// Enhanced solution structure to handle multiple complete solution sets
typedef struct rational_solutions {
    char **variable_names;      
    slong num_variables;        
    fmpq_t ***solution_sets;    
    slong *solutions_per_var;   
    slong num_solution_sets;    
    arb_t **solution_residuals;
    arb_t ***real_solution_sets;
    slong *real_solutions_per_var;
    slong num_real_solution_sets;
    arb_t **real_solution_residuals;
    slong real_solution_precision;
    char *real_root_summary;
    slong num_equations;        
    int is_valid;               
    int has_no_solutions;       
    char *error_message;        
    char *variable_order;       
    char *elimination_summary;  
    slong total_combinations;   
    slong successful_combinations; 
    slong num_base_solutions;   
    slong checked_solution_sets;    
    slong verified_solution_sets;   
    char **resultant_steps;     
    slong num_resultant_steps;  
    slong resultant_steps_cap;  
    char **candidate_solution_lines; 
    int *candidate_solution_pass;    
    slong num_candidate_solution_lines; 
    slong candidate_solution_lines_cap; 
} rational_solutions_t;

typedef struct rational_variable_info {
    char *name;
    slong max_degree;
    slong index;
} rational_variable_info_t;

typedef struct rational_equation_combination {
    slong *equation_indices;    
    slong num_equations;        
    int tried;                  
    int success;                
} rational_equation_combination_t;

void rational_solver_set_realtime_progress(int enabled);
void rational_solutions_init(rational_solutions_t *sols, slong num_vars);
void rational_solutions_clear(rational_solutions_t *sols);
void print_rational_solutions(const rational_solutions_t *sols);

int rational_contains_variable(const char *poly_str, const char *var_name);
int rational_contains_identifier_improved(const char *poly_str, const char *identifier);
char** rational_extract_variables_improved(char **poly_strings, slong num_polys, slong *num_vars_out);
slong rational_get_variable_max_degree_in_polynomial(const char *poly_str, const char *var_name);

void rational_generate_equation_combinations(slong num_equations, slong target_equations, 
                                             rational_equation_combination_t **combinations, slong *num_combinations);
void rational_free_equation_combinations(rational_equation_combination_t *combinations, slong num_combinations);

int rational_extract_and_sort_variables(char **poly_strings, slong num_polys,
                                         rational_variable_info_t **sorted_vars_out, slong *num_vars_out);

fmpq_t* rational_solve_univariate_equation_all_roots(const char *poly_str, const char *var_name,
                                                      slong *num_roots_out);
arb_t* rational_solve_univariate_equation_all_real_roots(const char *poly_str, const char *var_name,
                                                          slong *num_roots_out, slong prec);
char* rational_substitute_variable_in_polynomial(const char *poly_str, const char *var_name,
                                                  const fmpq_t value);
char* rational_eliminate_variable_dixon_with_selection(char **poly_strings, slong num_polys, 
                                                        const char **elim_vars, slong num_elim_vars,
                                                        rational_equation_combination_t *combination);

int rational_verify_solution_set(char **original_polys, slong num_polys,
                                 rational_variable_info_t *sorted_vars, slong num_vars,
                                 fmpq_t **solution_values);
void rational_filter_solutions_by_verification(rational_solutions_t *sols,
                                                char **original_polys, slong num_polys,
                                                rational_variable_info_t *sorted_vars);

int rational_solve_by_back_substitution_recursive_enhanced(char **original_polys, slong num_polys,
                                                             rational_variable_info_t *sorted_vars, slong num_vars,
                                                             fmpq_t *base_solutions, slong num_base_solutions,
                                                             rational_solutions_t *sols);

int rational_solve_by_elimination_enhanced(char **poly_strings, slong num_polys,
                                             rational_variable_info_t *sorted_vars, slong num_vars,
                                             rational_solutions_t *sols);

rational_solutions_t* solve_rational_polynomial_system_array_with_vars(char **poly_strings, slong num_polys,
                                                                       rational_variable_info_t *original_vars, slong num_original_vars);
rational_solutions_t* solve_rational_polynomial_system_array(char **poly_strings, slong num_polys);
rational_solutions_t* solve_rational_polynomial_system_string(const char *poly_string);

void test_rational_polynomial_solver(void);

#endif
