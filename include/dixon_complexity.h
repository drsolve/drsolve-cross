#ifndef DIXON_COMPLEXITY_H
#define DIXON_COMPLEXITY_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include "dixon_interface_flint.h"

// Define omega parameter for complexity calculation
#define DIXON_OMEGA 2.81

// Helper structure to collect polynomial information
typedef struct {
    char **all_vars;           // All variables found in polynomials
    slong num_all_vars;        // Total number of unique variables
    slong max_vars;            // Allocated space for variables
    long *degrees;             // Degree of each polynomial
    slong num_polys;           // Number of polynomials
    const fq_nmod_ctx_struct *ctx;
} poly_analysis_t;

// Variable hash table entry for fast variable lookups
typedef struct var_entry {
    char *name;
    slong index;
    struct var_entry *next;
} var_entry_t;

// Hash table for variable management
typedef struct {
    var_entry_t **buckets;
    slong bucket_count;
    slong count;
} var_hash_table_t;

// Lightweight parser state for degree calculation only
typedef struct {
    const char *input;
    size_t pos;
    size_t len;
    var_hash_table_t var_table;
    long max_degree_found;
    const fq_nmod_ctx_struct *ctx;
    const char *generator_name;
} lightweight_parser_t;

typedef struct {
    slong det_size;
    double det_size_log2;
    double det_factorial_log2;
    slong num_all_vars;
    slong num_elim_vars;
    slong num_parameter_vars;
    slong step1_var_count;
    long step1_det_total_degree;
    double step1_kronecker_degree_log2;
    double step1_direct_log2;
    double step1_direct_factorial_log2;
    double step1_direct_fft_log2;
    double step1_hnf_log2;
    double step1_hnf_linear_algebra_log2;
    double step1_hnf_degree_density_log2;
    double step1_sparse_log2;
    double step1_sparse_term_bound_log2;
    double step1_sparse_slp_length_log2;
    long step1_sparse_param_degree_bound;
    long step1_sparse_partial_degree_bound;
    double step1_sparse_success_prob_lb;
    double step1_sparse_retry_factor;
    double step1_sparse_expected_log2;
    double step1_sparse_q_for_three_quarters_log2;
    slong macaulay_degree;
    slong macaulay_rows;
    slong macaulay_cols;
    slong macaulay_square_size;
    double macaulay_log2;
    slong grobner_dreg;
    double grobner_log2;
    double step4_log2;
    double total_direct_log2;
    double total_hnf_log2;
    double total_sparse_log2;
} dixon_complexity_report_t;

// Function declarations

// Basic utility functions
int compare_desc(const void *a, const void *b);

// Core Dixon complexity calculations
void dixon_size(fmpz_t result, const long *a_values, int len, int show_details);
double dixon_complexity(const long *a_values, int len, int n, double omega);

long get_poly_total_degree(const char *poly_str, const char *gen_name);
void collect_variables(const char **polys, slong npolys,
                               const char *gen_name,
                               char ***vars_out, slong *nvars_out);

// Main Dixon complexity analysis functions
char* dixon_complexity_auto(const char **poly_strings, slong num_polys,
                           const char **elim_vars, slong num_elim_vars,
                           const fq_nmod_ctx_t ctx);
char* dixon_complexity_auto_str(const char *poly_string, const char *vars_string, const fq_nmod_ctx_t ctx);

void run_complexity_analysis(const char *polys_str,
                             const char *vars_str,
                             const fmpz_t prime,
                             ulong power,
                             const fq_nmod_ctx_t ctx,
                             const char *output_filename,
                             int silent_mode,
                             double comp_time,
                             double omega);

void run_complexity_analysis_from_degrees(const long *degrees,
                                          slong num_polys,
                                          slong num_all_vars,
                                          slong num_elim_vars,
                                          const fmpz_t prime,
                                          ulong power,
                                          const fq_nmod_ctx_t ctx,
                                          const char *output_filename,
                                          int silent_mode,
                                          double comp_time,
                                          double omega,
                                          const char *system_spec);

void dixon_complexity_report_from_degrees(dixon_complexity_report_t *report,
                                          const long *degrees,
                                          slong num_polys,
                                          slong num_all_vars,
                                          slong num_elim_vars,
                                          slong num_parameter_vars,
                                          const fmpz_t field_order,
                                          double omega);

// Complexity extraction functions
double extract_max_complexity(const char **poly_strings, slong num_polys);
double extract_max_complexity_str(const char *poly_string);

// Test function
int test_dixon_complexity(void);

#endif // DIXON_COMPLEXITY_H
