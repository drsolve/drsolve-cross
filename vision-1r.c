#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <flint/flint.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>
#include <flint/fq_nmod_mat.h>
#include <flint/fmpz.h>

#include "dixon_test.h"

const char* p0 = "z8*x0*x1 + (z8^7 + z8^5 + z8^3 + z8)*x1 + 1";

const char* p1 = "(z8^2 + z8)*x0*x2 + (z8^6 + z8^5 + z8^4 + z8^3 + z8 + 1)*x2 + 1";

const char* p2 = "(z8^6 + z8^5 + z8^3)*x1^4 + (z8^6 + z8^4 + z8^3 + z8^2)*x2^4 + (z8^6 + z8^4 + z8^3 + z8^2)*x1^2 + (z8^6 + z8^5 + z8^4 + z8)*x2^2 + (z8^7 + z8^6 + z8^4)*x1 + (z8^7 + z8^5 + z8^4 + z8^3)*x2 + (z8^5 + z8^4 + z8^3 + 1)";

int main() {
    printf("=== Vision Attack using DIXON Syntax ===\n\n");
    
    fq_nmod_ctx_t ctx;
    nmod_poly_t modulus;
    nmod_poly_init(modulus, 2);
    nmod_poly_set_coeff_ui(modulus, 0, 1); 
    nmod_poly_set_coeff_ui(modulus, 2, 1);
    nmod_poly_set_coeff_ui(modulus, 3, 1); 
    nmod_poly_set_coeff_ui(modulus, 4, 1); 
    nmod_poly_set_coeff_ui(modulus, 8, 1);
    fq_nmod_ctx_init_modulus(ctx, modulus, "z8");
    nmod_poly_clear(modulus);
    
    char* r1 = RESULTANT((p0, p1), ("x0"));
    char* d = RESULTANT((r1, p2), ("x1"));
    
    fq_nmod_ctx_clear(ctx);
    printf("\n=== Computation Complete ===\n");
    return 0;
}