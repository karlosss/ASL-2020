#ifndef ASL_PROJECT_CLOCK_H
#define ASL_PROJECT_CLOCK_H

#define START_CLOCK(hi, lo) asm volatile ("CPUID\n\t" \
                                          "RDTSC\n\t" \
                                          "mov %%edx, %0\n\t" \
                                          "mov %%eax, %1\n\t": "=r" (hi), "=r"(lo)::"%rax", "%rbx", "%rcx", "%rdx")

#define STOP_CLOCK(hi, lo) asm volatile("RDTSCP\n\t" \
                                        "mov %%edx, %0\n\t" \
                                        "mov %%eax, %1\n\t" \
                                        "CPUID\n\t": "=r" (hi), "=r"(lo)::"%rax", "%rbx", "%rcx", "%rdx")

#endif //ASL_PROJECT_CLOCK_H
