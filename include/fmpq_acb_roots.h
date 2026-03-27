#ifndef FMPQ_ACB_ROOTS_H
#define FMPQ_ACB_ROOTS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly.h>

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

void fmpq_roots_init(fmpq_roots_t *roots);
void fmpq_roots_clear(fmpq_roots_t *roots);

slong fmpq_poly_roots(fmpq_roots_t *roots, const fmpq_poly_t poly, int with_multiplicity);

void fmpq_roots_print(const fmpq_roots_t *roots);
char* fmpq_roots_to_string(const fmpq_roots_t *roots);

#endif
