#include <stdio.h>
#include "measure.h"
#include "base.h"

INLINE void measured(){
    int i;
    volatile double x = 0;
    for(i = 0; i < 1e7; ++i){
        x += 1;
    }
}

int main() {
    size_t num_iter = 5;
    size_t N = 100;
    size_t M = 200;
    size_t T = 500;
    
    unsigned long cycles = perf(baum_welch_base, N, M, T, num_iter);    
    
    printf("%lu\n", cycles);
    return 0;
}
