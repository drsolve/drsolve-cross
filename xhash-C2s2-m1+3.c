#include "dixon_test.h"

const char *ideal = "x2^3 = -3*x0^3 + 4*x0^2*x1 + 8*x0^2 - 3*x0*x1^2 - 8*x0*x1 - x0 + 2*x1^3 - x1^2 + 3*x1, \
x3^3 = 2*x0^3 - x0^2*x1 + x0^2 + 2*x0*x1^2 - 2*x0*x1 - 8*x0 - 3*x1^3 - x1^2 - 6*x1 - 2, \
x4^3 = 7*x0^3 - x0^2*x1 - 7*x0^2 + 6*x0*x1^2 - 4*x0*x1 + 8*x0 + 2*x1^3 - 6*x1^2 - 7*x1 - 1, \
x5^3 = 2*x0^3 + 4*x0^2*x1 + x0*x1^2 + 7*x0*x1 - 2*x0 + 7*x1^3 - 6*x1^2 + 2*x1 + 1";

const char* f0 = "-3*x0^3 + 4*x0^2*x1 + 8*x0^2 - 3*x0*x1^2 - 8*x0*x1 - x0 + 2*x1^3 - x1^2 + 3*x1 - x2^3";
const char* f1 = "2*x0^3 - x0^2*x1 + x0^2 + 2*x0*x1^2 - 2*x0*x1 - 8*x0 - 3*x1^3 - x1^2 - 6*x1 - x3^3 - 2";
const char* f2 = "7*x0^3 - x0^2*x1 - 7*x0^2 + 6*x0*x1^2 - 4*x0*x1 + 8*x0 + 2*x1^3 - 6*x1^2 - 7*x1 - x4^3 - 1";
const char* f3 = "2*x0^3 + 4*x0^2*x1 + x0*x1^2 + 7*x0*x1 - 2*x0 + 7*x1^3 - 6*x1^2 + 2*x1 - x5^3 + 1";
const char* f4 = "4*x2^3 + 3*x2^2*x3 + 4*x2^2 + 6*x2*x3^2 + 2*x2*x3 - 6*x2 + 7*x3^3 - 6*x3^2 + 5*x3 + 8*x4 + 6*x5 - 4";
const char* f5 = "4*x2^3 - 3*x2^2*x3 + x2^2 - 8*x2*x3^2 + 3*x2*x3 - 2*x2 - 3*x3^3 + 3*x3^2 - 3*x3 - 8*x4 + 8*x5 - 1";

int main() {    
    fq_nmod_ctx_t ctx;
    mp_limb_t prime = 4611686018427388039; 
    fmpz_t p;
    fmpz_init_set_ui(p, prime);
    fq_nmod_ctx_init(ctx, p, 1, "t");   
    
    clock_t total_start = clock();
    
    char* r1 = DIXON_WITH_IDEAL((f5, f3), ("x5"));
    char* r2 = DIXON_WITH_IDEAL((r1, f2), ("x4"));
    char* r3 = DIXON_WITH_IDEAL((r2, f1), ("x3"));
    char* r4 = DIXON_WITH_IDEAL((r3, f0), ("x2"));
    
    char* l1 = DIXON_WITH_IDEAL((f4, f3), ("x5"));
    char* l2 = DIXON_WITH_IDEAL((l1, f2), ("x4"));
    char* l3 = DIXON_WITH_IDEAL((l2, f1), ("x3"));
    char* l4 = DIXON_WITH_IDEAL((l3, f0), ("x2"));
    
    char* r = RESULTANT((r4,l4),("x1"));
    clock_t total_end = clock();
    double total_time = (double)(total_end - total_start) / CLOCKS_PER_SEC;
    
    printf("Total computation time: %.3f seconds\n", total_time);
    
    
    // Clear field context
    fq_nmod_ctx_clear(ctx);
    
    printf("\n=== Computation Complete ===\n");
    
    return 0;
} 