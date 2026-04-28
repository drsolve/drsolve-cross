#ifndef MACAULAY_FLINT_H
#define MACAULAY_FLINT_H

#include <flint/flint.h>
#include "fq_mvpoly.h"

void fq_macaulay_resultant(fq_mvpoly_t *result, fq_mvpoly_t *polys,
                           slong nvars, slong npars);

void fq_macaulay_resultant_with_names(fq_mvpoly_t *result, fq_mvpoly_t *polys,
                                      slong nvars, slong npars,
                                      char **var_names, char **par_names,
                                      const char *gen_name);

#endif
