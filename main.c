#include <stdio.h>
#include "measure.h"

INLINE void measured(){
    int i;
    volatile double x = 0;
    for(i = 0; i < 1e7; ++i){
        x += 1;
    }
}

int main() {
    unsigned long cycles = perf(measured);
    printf("%lu\n", cycles);
    return 0;
}
