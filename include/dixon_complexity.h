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
#define DIXON_OMEGA 2.3

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

// Complexity extraction functions
double extract_max_complexity(const char **poly_strings, slong num_polys);
double extract_max_complexity_str(const char *poly_string);

// Test function
int test_dixon_complexity(void);

#endif // DIXON_COMPLEXITY_H
