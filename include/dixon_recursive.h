#ifndef DIXON_RECURSIVE_H
#define DIXON_RECURSIVE_H

#include <flint/flint.h>
#include "fq_mvpoly.h"

void fq_dixon_fast_resultant(fq_mvpoly_t *result, fq_mvpoly_t *polys,
                             slong nvars, slong npars);

void fq_dixon_fast_resultant_with_names(fq_mvpoly_t *result, fq_mvpoly_t *polys,
                                        slong nvars, slong npars,
                                        char **var_names, char **par_names,
                                        const char *gen_name);

#endif
