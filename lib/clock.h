#ifndef ASL_PROJECT_CLOCK_H
#define ASL_PROJECT_CLOCK_H

#define COUNTER_LO(a) ((a).int32.lo)
#define COUNTER_HI(a) ((a).int32.hi)
#define COUNTER_VAL(a) ((a).int64)

#define START_CLOCK(c) asm volatile ("CPUID\n\t" \
                                    "RDTSC\n\t" \
                                    "mov %%edx, %0\n\t" \
                                    "mov %%eax, %1\n\t": "=r" (COUNTER_HI(c)), "=r"(COUNTER_LO(c))::"%rax", "%rbx", "%rcx", "%rdx")

#define STOP_CLOCK(c) asm volatile( "RDTSCP\n\t" \
                                    "mov %%edx, %0\n\t" \
                                    "mov %%eax, %1\n\t" \
                                    "CPUID\n\t": "=r" (COUNTER_HI(c)), "=r"(COUNTER_LO(c))::"%rax", "%rbx", "%rcx", "%rdx")

typedef union{
    size_t int64;
	struct {unsigned int lo, hi;} int32;
} tsc_counter;

static size_t start_clock() {
    tsc_counter start;
    START_CLOCK(start);
    return COUNTER_VAL(start);
}

static size_t cycles_elapsed(size_t start) {
	tsc_counter end;
	STOP_CLOCK(end);
	return COUNTER_VAL(end) - start;
}

#endif //ASL_PROJECT_CLOCK_H
