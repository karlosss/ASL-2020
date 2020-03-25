#include <bits/stdc++.h>
#include "common.h"
#include "clock.h"

#define REP 100

static int size_t_cmp(const void *a, const void *b) {
    if (*(size_t *) a > *(size_t *) b) return 1;
    else if (*(size_t *) a < *(size_t *) b) return -1;
    else return 0;
}

size_t perf(size_t N, size_t M, size_t T, size_t num_iter) {
    int *O; // observation
    double *PI; // Initial probabilities
    double *A; // Transition probabilities
    double *B; // Emission probabilities
    double *FW, *BW, *C;
    generate_input(N, M, T, &O, &PI, &A, &B, &FW, &BW, &C);

    unsigned int cycles_low, cycles_high, cycles_low1, cycles_high1;
    size_t start, end;
    int i;
    size_t time[REP];

    for (i = 0; i < REP; i++) {
        START_CLOCK(cycles_high, cycles_low);
        baum_welch(PI, A, B, O, FW, BW, C, N, M, T, num_iter);
        STOP_CLOCK(cycles_high1, cycles_low1);

        start = (((size_t) cycles_high << 32) | cycles_low);
        end = (((size_t) cycles_high1 << 32) | cycles_low1);

        time[i] = end - start;
    }

    free(O);
    free(PI);
    free(A);
    free(B);
    free(FW);
    free(BW);
    free(C);

    qsort(time, REP, sizeof(size_t), size_t_cmp);
    return time[REP / 2];
}
