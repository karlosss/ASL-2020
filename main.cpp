#include <iostream>
#include <random>
#include <cassert>

#include "base.h"

#define EPS (1e-6)
size_t N, M, T;

void generate_observation(int **PO, size_t len)
{   
    *PO = static_cast<int *>(malloc(len * sizeof(int)));
    int* O = *PO;

    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_int_distribution<> dist(0, M-1);
    
    for (size_t i = 0; i < len; ++i)  
        O[i] = dist(gen);
}

void generate_m(double ** PM, size_t row, size_t col)
{   
    *PM = static_cast<double *>(malloc(row*col * sizeof(double)));
    double* M = *PM;

    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<double> dist(1., 2.);

    for (size_t i = 0; i < row; ++i) {
        double sum = 0.;
        for(size_t j = 0; j < col; ++j) {
            M[i*row + j] = dist(gen);
            sum += M[i*row + j];
        }
        assert(sum > 0.);
        for(size_t j = 0; j < col; ++j) {
            M[i*row + j] /= sum;
        }
        double sum_prob = 0.;
        for(size_t j = 0; j < col; ++j) {
            sum_prob += M[i*row + j];
        }
        assert(abs(sum_prob - 1.) < EPS);
    }   
}

void print_m(double* M, size_t rows, size_t cols){
    for(size_t i = 0; i < rows; i++) {
        for(size_t j = 0; j < cols; j++) {
            std::cout << M[i*rows + j] << " ";
        }
        std::cout << '\n';
    }
}

void test_base() {
    N = 2; M = 2; T = 3;

    int *O; // observation
    double *PI, *A, *B; // Initial, Transition and Emission probabilities

    generate_observation(&O, T);
    generate_m(&PI, 1, N);
    generate_m(&A, N, N);
    generate_m(&B, N, M);
    /*std::cout << "Printing PI:" << '\n';
    print_m(PI, 1, N);
    std::cout << "Printing A:" << '\n';
    print_m(A, N, N);
    std::cout << "Printing B:" << '\n';
    print_m(B, N, M);*/

    double* FW = static_cast<double *>(malloc(N*T * sizeof(double)));
    double* BW = static_cast<double *>(malloc(N*T * sizeof(double)));
    double* C = static_cast<double *>(malloc(T * sizeof(double)));
    baum_welch_test(PI, A, B, O, FW, BW, C, N, M, T, 1);

    free(O); free(PI); free(A); free(B);

    return;
}

int main() {
    test_base();
}
