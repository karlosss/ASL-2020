#include <iostream>
#include <random>
#include <cassert>

#include "base.h"
#include "common.h"

int main() {
    size_t num_iter = 5;
    size_t limit = 50;
    
    
    for(size_t N = 2; N < limit; N++) {
    for(size_t M = 2; M < limit; M++) {
    for(size_t T = 2; T < limit; T++) {
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
        baum_welch_test(PI, A, B, O, FW, BW, C, N, M, T, num_iter);

        free(O); free(PI); free(A); free(B);
        free(FW); free(BW); free(C);
    }}}
    std::cout << "Testing completed. No errors found!" << '\n';
    return 0;
}