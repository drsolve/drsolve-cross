#ifndef FMPQ_ACB_ROOTS_H
#define FMPQ_ACB_ROOTS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly.h>
#include <flint/acb.h>
#include <flint/acb_poly.h>

#define FMPQ_ROOT_SEARCH_MAX_DEGREE 1000
#define FMPQ_ROOT_SEARCH_MAX_CANDIDATES 1000000

typedef struct {
    fmpq_t *roots;
    slong *multiplicities;
    slong num_roots;
    slong alloc;
} fmpq_roots_t;

typedef fmpq_roots_t *fmpq_roots_ptr;
typedef const fmpq_roots_t *fmpq_roots_srcptr;

typedef struct {
    acb_t *roots;
    slong *multiplicities;
    slong num_roots;
    slong alloc;
} acb_roots_t;

typedef acb_roots_t *acb_roots_ptr;
typedef const acb_roots_t *acb_roots_srcptr;

typedef struct {
    arb_t *roots;
    slong *multiplicities;
    slong num_roots;
    slong alloc;
} arb_roots_t;

typedef arb_roots_t *arb_roots_ptr;
typedef const arb_roots_t *arb_roots_srcptr;

typedef struct {
    fmpq_roots_t rational_roots;
    acb_roots_t approximate_roots;
    arb_roots_t real_roots;
} fmpq_acb_roots_t;

void fmpq_roots_init(fmpq_roots_t *roots);
void fmpq_roots_clear(fmpq_roots_t *roots);

slong fmpq_poly_roots(fmpq_roots_t *roots, const fmpq_poly_t poly, int with_multiplicity);

void fmpq_roots_print(const fmpq_roots_t *roots);
char* fmpq_roots_to_string(const fmpq_roots_t *roots);

void acb_roots_init(acb_roots_t *roots);
void acb_roots_clear(acb_roots_t *roots);

slong acb_poly_roots(acb_roots_t *roots, const acb_poly_t poly, slong prec);
slong fmpq_poly_acb_roots(acb_roots_t *roots, const fmpq_poly_t poly, slong prec);

void acb_roots_print(const acb_roots_t *roots);

void arb_roots_init(arb_roots_t *roots);
void arb_roots_clear(arb_roots_t *roots);

slong acb_roots_to_real(arb_roots_t *real_roots, const acb_roots_t *complex_roots, slong prec);

void arb_roots_print(const arb_roots_t *roots);

void fmpq_acb_roots_init(fmpq_acb_roots_t *roots);
void fmpq_acb_roots_clear(fmpq_acb_roots_t *roots);

slong fmpq_poly_all_roots(fmpq_acb_roots_t *roots, const fmpq_poly_t poly, slong prec);
slong fmpq_poly_real_roots(fmpq_acb_roots_t *roots, const fmpq_poly_t poly, slong prec);

void fmpq_acb_roots_print(const fmpq_acb_roots_t *roots);
void fmpq_acb_roots_print_real(const fmpq_acb_roots_t *roots);

void fmpq_roots_print_to_file(FILE *fp, const fmpq_roots_t *roots);
void acb_roots_print_to_file(FILE *fp, const acb_roots_t *roots);
void arb_roots_print_to_file(FILE *fp, const arb_roots_t *roots);
void fmpq_acb_roots_print_all_to_file(FILE *fp, const fmpq_acb_roots_t *roots);

#endif
