#include <stdio.h>
#include "measure.h"
#include "base.h"

int main() {
    size_t num_iter = 5;
    size_t N = 10;
    size_t M = 10;
    size_t T = 10;
    
    unsigned long cycles = perf(baum_welch_base, N, M, T, num_iter);    
    
    printf("%lu\n", cycles);
    return 0;
}
