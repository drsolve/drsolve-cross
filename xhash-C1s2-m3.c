#include "dixon_test.h"

const char *ideal = "x1^3 = -7*x0^3 + 8*x0^2 + 8*x0 + 4, \
x2^3 = -x0^3 + 3*x0^2 - 3*x0 - 5, \
x3^3 = 8*x0^3 + 5*x0^2 - 8*x0 - 7, \
x4 = -3*x1^2*x2 + 4*x1^2*x3 - 4*x1^2 + 5*x1*x2^2 + 5*x1*x2*x3 - 8*x1*x2 + x1*x3^2 + 2*x1*x3 - 5*x1 + 5*x2^3 + 6*x2^2*x3 + 6*x2^2 - 7*x2*x3^2 - 5*x2*x3 + 7*x2 - 7*x3^3 + 5*x3^2 - 7*x3 + 7, \
x5 = -x1^3 - 3*x1^2*x2 - 2*x1^2*x3 + 7*x1^2 - 4*x1*x2^2 + x1*x2*x3 - 6*x1*x2 - 3*x1*x3^2 + 4*x1*x3 + 8*x1 - 3*x2^3 - x2^2*x3 - 2*x2^2 + 7*x2*x3^2 - x2*x3 + 8*x2 + 7*x3^3 + 2*x3^2 - 2*x3, \
x6 = -7*x1^3 - 7*x1^2*x2 - 2*x1^2*x3 - 7*x1^2 - 2*x1*x2^2 - 3*x1*x2 - x1*x3^2 - 7*x1*x3 - 4*x1 + 3*x2^3 - 3*x2^2*x3 - 4*x2^2 + x2*x3^2 + 3*x2*x3 - 8*x2 + 5*x3^3 - 5*x3 + 5";

const char* f0 = "-7*x0^3 + 8*x0^2 + 8*x0 - x1^3 + 4";
const char* f1 = "-x0^3 + 3*x0^2 - 3*x0 - x2^3 - 5";
const char* f2 = "8*x0^3 + 5*x0^2 - 8*x0 - x3^3 - 7";
const char* f3 = "-3*x1^2*x2 + 4*x1^2*x3 - 4*x1^2 + 5*x1*x2^2 + 5*x1*x2*x3 - 8*x1*x2 + x1*x3^2 + 2*x1*x3 - 5*x1 + 5*x2^3 + 6*x2^2*x3 + 6*x2^2 - 7*x2*x3^2 - 5*x2*x3 + 7*x2 - 7*x3^3 + 5*x3^2 - 7*x3 - x4 + 7";
const char* f4 = "-x1^3 - 3*x1^2*x2 - 2*x1^2*x3 + 7*x1^2 - 4*x1*x2^2 + x1*x2*x3 - 6*x1*x2 - 3*x1*x3^2 + 4*x1*x3 + 8*x1 - 3*x2^3 - x2^2*x3 - 2*x2^2 + 7*x2*x3^2 - x2*x3 + 8*x2 + 7*x3^3 + 2*x3^2 - 2*x3 - x5";
const char* f5 = "-7*x1^3 - 7*x1^2*x2 - 2*x1^2*x3 - 7*x1^2 - 2*x1*x2^2 - 3*x1*x2 - x1*x3^2 - 7*x1*x3 - 4*x1 + 3*x2^3 - 3*x2^2*x3 - 4*x2^2 + x2*x3^2 + 3*x2*x3 - 8*x2 + 5*x3^3 - 5*x3 - x6 + 5";
const char* f6 = "-7*x4^3 - 7*x4^2*x5 + 2*x4^2*x6 - 7*x4^2 + 7*x4*x5^2 + 6*x4*x5*x6 + 2*x4*x5 + 8*x4*x6^2 + 5*x4*x6 - x4 + 8*x5^3 - 7*x5^2*x6 - 8*x5^2 - 8*x5*x6^2 - 5*x5*x6 - x6^3 + x6^2 - 7*x6 - 6";

int main() {    
    fq_nmod_ctx_t ctx;
    mp_limb_t prime = 4611686018427388039; 
    fmpz_t p;
    fmpz_init_set_ui(p, prime);
    fq_nmod_ctx_init(ctx, p, 1, "t");   
    
    clock_t total_start = clock();
    
    char* r1 = DIXON_WITH_IDEAL((f6, f5), ("x6"));
    char* r2 = DIXON_WITH_IDEAL((r1, f4), ("x5"));
    char* r3 = DIXON_WITH_IDEAL((r2, f3), ("x4"));
    char* r4 = DIXON_WITH_IDEAL((r3, f2), ("x3"));
    char* r5 = DIXON_WITH_IDEAL((r4, f1), ("x2"));
    char* r6 = DIXON_WITH_IDEAL((r5, f0), ("x1"));
    
    clock_t total_end = clock();
    double total_time = (double)(total_end - total_start) / CLOCKS_PER_SEC;
    
    printf("Total computation time: %.3f seconds\n", total_time);
    
    
    // Clear field context
    fq_nmod_ctx_clear(ctx);
    
    printf("\n=== Computation Complete ===\n");
    
    return 0;
} 