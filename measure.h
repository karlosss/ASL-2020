#ifndef ASL_PROJECT_MEASURE_H
#define ASL_PROJECT_MEASURE_H

#include <stdlib.h>
#include "common.h"

#define REP 100

// Returns the number of cycles required to run the function func
unsigned long perf(fn func, size_t N, size_t M, size_t T, size_t num_iter) {
    int* O; // observation
    double* PI; // Initial probabilities
    double* A; // Transition probabilities
    double* B; // Emission probabilities

    generate_observation(&O, T, M);
    generate_m(&PI, 1, N);
    generate_m(&A, N, N);
    generate_m(&B, N, M);

    double* FW = static_cast<double *>(malloc(N*T * sizeof(double)));
    double* BW = static_cast<double *>(malloc(N*T * sizeof(double)));
    double* C = static_cast<double *>(malloc(T * sizeof(double)));

    unsigned int cycles_low, cycles_high, cycles_low1, cycles_high1;
    unsigned long start, end;
    int i;
    unsigned long time[REP];

    for (i = 0; i < REP; i++) {
        asm volatile ("CPUID\n\t"
                      "RDTSC\n\t"
                      "mov %%edx, %0\n\t"
                      "mov %%eax, %1\n\t": "=r" (cycles_high), "=r"(cycles_low)::"%rax", "%rbx", "%rcx", "%rdx");
        func(PI, A, B, O, FW, BW, C, N, M, T, num_iter);
        asm volatile("RDTSCP\n\t"
                     "mov %%edx, %0\n\t"
                     "mov %%eax, %1\n\t"
                     "CPUID\n\t": "=r" (cycles_high1), "=r"(cycles_low1)::"%rax", "%rbx", "%rcx", "%rdx");

        start = ( ((unsigned long)cycles_high << 32) | cycles_low );
        end = ( ((unsigned long)cycles_high1 << 32) | cycles_low1 );

        time[i] = end-start;
    }

    free(O); free(PI); free(A); free(B);
    free(FW); free(BW); free(C);

    qsort(time, REP, sizeof(unsigned long), double_cmp);
    return time[REP/2];
}

#endif //ASL_PROJECT_MEASURE_H
