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

const char* p0 = "(z8^7 + z8^6 + z8^5 + 1)*x1^4*x3^4 + (z8^4 + z8^3 + z8^2 + z8 + 1)*x2^4*x3^4 + (z8^6 + 1)*x1^4*x3^2 + (z8^7 + z8^6 + z8^5 + z8^3 + z8^2 + z8 + 1)*x2^4*x3^2 + (z8^6 + 1)*x1^2*x3^4 + (z8^7 + z8^6 + z8^5 + z8^3 + z8^2 + z8 + 1)*x2^2*x3^4 + (z8^7 + z8^4 + z8^3 + z8^2 + 1)*x1^4*x3 + (z8^6 + z8^4 + z8^3 + z8^2 + 1)*x2^4*x3 + (z8^7 + z8^4 + z8^3 + z8^2 + 1)*x1*x3^4 + (z8^6 + z8^4 + z8^3 + z8^2 + 1)*x2*x3^4 + (z8^7 + z8^4 + z8^2)*x1^2*x3^2 + (z8^7 + z8^6 + z8^4 + z8^3 + z8^2 + z8)*x2^2*x3^2 + (z8^4 + z8^3 + 1)*x3^4 + (z8^7 + z8^6 + z8^4 + z8^3 + z8)*x1^2*x3 + (z8^7 + z8^5 + z8^4 + z8^2 + z8 + 1)*x2^2*x3 + (z8^7 + z8^6 + z8^4 + z8^3 + z8)*x1*x3^2 + (z8^7 + z8^5 + z8^4 + z8^2 + z8 + 1)*x2*x3^2 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^2)*x1*x3 + (z8^7 + z8^3 + z8^2 + z8)*x2*x3 + (z8^6 + z8^5 + z8^3 + z8)*x3^2 + (z8^6 + z8^3 + z8^2 + z8 + 1)*x3 + 1";

const char* p1 = " (z8^5 + z8^4 + z8^3 + z8^2 + z8)*x1^4*x4^4 + (z8^7 + z8^6)*x2^4*x4^4 + (z8^7 + z8^6 + z8 + 1)*x1^4*x4^2 + (z8^6 + z8^5 + z8^3 + z8^2 + 1)*x2^4*x4^2 + (z8^7 + z8^6 + z8 + 1)*x1^2*x4^4 + (z8^6 + z8^5 + z8^3 + z8^2 + 1)*x2^2*x4^4 + (z8^7 + z8^5 + z8^4 + z8^3 + z8)*x1^4*x4 + (z8^6 + z8^5 + z8^4 + z8^3 + z8)*x2^4*x4 + (z8^7 + z8^5 + z8^4 + z8^3 + z8)*x1*x4^4 + (z8^6 + z8^5 + z8^4 + z8^3 + z8)*x2*x4^4 + (z8^7 + z8^5 + 1)*x1^2*x4^2 + (z8^7 + z8^6 + z8^5 + z8^3 + z8 + 1)*x2^2*x4^2 + (z8^6 + z8^5 + z8 + 1)*x4^4 + (z8^6 + z8^5 + z8^4 + z8 + 1)*x1^2*x4 + (z8^4 + z8^3 + z8^2 + z8)*x2^2*x4 + (z8^6 + z8^5 + z8^4 + z8 + 1)*x1*x4^2 + (z8^4 + z8^3 + z8^2 + z8)*x2*x4^2 + x1*x4 + (z8^6 + z8^5 + z8^4 + z8^3 + z8 + 1)*x2*x4 + (z8^6 + z8^5 + z8^4 + z8^2)*x4^2 + (z8^5 + z8^4 + z8^2)*x4 + 1";

const char* p2 = "z8*x0*x1 + z8*x1 + 1";

const char* p3 = "(z8^2 + z8)*x0*x2 + (z8^7 + z8^5 + z8^3 + z8 + 1)*x2 + 1";

const char* p4 = "z8*x3*x5 + (z8 + 1)*x4*x5 + (z8^5 + z8^4 + z8^2)*x5 + 1";

const char* p5 = "(z8^2 + z8)*x3*x6 + (z8^2 + z8 + 1)*x4*x6 + (z8^6 + z8^3 + z8^2 + z8)*x6 + 1";

const char* p6 = "(z8^7 + z8^6 + 1)*x5^4 + (z8^5 + z8^3 + z8^2 + z8 + 1)*x6^4 + (z8^5 + z8^3 + 1)*x5^2 + (z8^7 + z8^5 + z8^4 + z8 + 1)*x6^2 + (z8^6 + z8^5 + z8^4 + z8^3 + z8^2 + 1)*x5 + (z8^7 + z8^6 + z8^3 + z8^2 + 1)*x6 + (z8^6 + z8^4 + z8)";

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
    
    char* r1 = RESULTANT((p2, p3), ("x0"));
    char* r2 = DIXON((r1, p0, p1), ("x1", "x2"));
    char* s1 = RESULTANT((p6, p5), ("x6"));
    char* s2 = RESULTANT((s1, p4), ("x5"));
    char* d = RESULTANT((r2, s2), ("x3"));    
    
    fq_nmod_ctx_clear(ctx);
    
    printf("\n=== Computation Complete ===\n");
    
    return 0;
}