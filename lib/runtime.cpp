#include <bits/stdc++.h>
#include "common.h"
#include "clock.h"

#define REP 100
#define CYCLES_REQUIRED 1e8


size_t runtime(int N, int M, int T, int num_iter) {
    int *O; // observation
    double *PI; // Initial probabilities
    double *A; // Transition probabilities
    double *B; // Emission probabilities
    double *FW, *BW, *C;
    generate_input(N, M, T, &O, &PI, &A, &B, &FW, &BW, &C);

    size_t start, cycles;
    size_t num_runs = 10;
    double multiplier = 1;
    size_t cycle_measurements[REP];
    do {
        num_runs = num_runs * multiplier;
        start = start_clock();
        for (size_t i = 0; i < num_runs; i++) {
            baum_welch(PI, A, B, O, FW, BW, C, N, M, T, num_iter);
        }
        cycles = cycles_elapsed(start);

        multiplier = (CYCLES_REQUIRED) / (double) cycles;

    } while (multiplier > 2);

    for (size_t j = 0; j < REP; j++) {

        start = start_clock();
        for (size_t i = 0; i < num_runs; ++i) {
            baum_welch(PI, A, B, O, FW, BW, C, N, M, T, num_iter);
        }
        cycles = cycles_elapsed(start);

        cycle_measurements[j] = cycles / num_runs;
    }

    free(O);
    free(PI);
    free(A);
    free(B);
    free(FW);
    free(BW);
    free(C);

    std::sort(cycle_measurements, cycle_measurements + REP);
    return cycle_measurements[REP / 2];
}
