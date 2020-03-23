#ifndef ASL_PROJECT_MEASURE_H
#define ASL_PROJECT_MEASURE_H

#include <stdlib.h>
#include "common.h"

#define REP 100

unsigned long perf(fn func) {
    unsigned int cycles_low, cycles_high, cycles_low1, cycles_high1;
    unsigned long start, end;
    int i;
    unsigned long time[REP];

    asm volatile("CPUID\n\t"
                 "RDTSC\n\t"
                 "mov %%edx, %0\n\t"
                 "mov %%eax, %1\n\t": "=r" (cycles_high), "=r" (cycles_low)::"%rax", "%rbx", "%rcx", "%rdx");
    asm volatile("RDTSCP\n\t"
                 "mov %%edx, %0\n\t"
                 "mov %%eax, %1\n\t"
                 "CPUID\n\t": "=r" (cycles_high1), "=r" (cycles_low1)::"%rax", "%rbx", "%rcx", "%rdx");
    asm volatile("CPUID\n\t"
                 "RDTSC\n\t"
                 "mov %%edx, %0\n\t"
                 "mov %%eax, %1\n\t": "=r" (cycles_high), "=r" (cycles_low)::"%rax", "%rbx", "%rcx", "%rdx");
    asm volatile("RDTSCP\n\t"
                 "mov %%edx, %0\n\t"
                 "mov %%eax, %1\n\t"
                 "CPUID\n\t": "=r" (cycles_high1), "=r" (cycles_low1)::"%rax", "%rbx", "%rcx", "%rdx");


    for (i = 0; i < REP; i++) {
        asm volatile ("CPUID\n\t"
                      "RDTSC\n\t"
                      "mov %%edx, %0\n\t"
                      "mov %%eax, %1\n\t": "=r" (cycles_high), "=r"(cycles_low)::"%rax", "%rbx", "%rcx", "%rdx");
        func();
        asm volatile("RDTSCP\n\t"
                     "mov %%edx, %0\n\t"
                     "mov %%eax, %1\n\t"
                     "CPUID\n\t": "=r" (cycles_high1), "=r"(cycles_low1)::"%rax", "%rbx", "%rcx", "%rdx");

        start = ( ((unsigned long)cycles_high << 32) | cycles_low );
        end = ( ((unsigned long)cycles_high1 << 32) | cycles_low1 );

        time[i] = end-start;
    }

    qsort(time, REP, sizeof(unsigned long), double_cmp);
    return time[REP/2];
}

#endif //ASL_PROJECT_MEASURE_H
