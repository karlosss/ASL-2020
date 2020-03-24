#include <math.h>
#include <random>
#include <cassert>
#include <iostream>

void generate_observation(int **PO, size_t len, size_t M) {   
    *PO = static_cast<int *>(malloc(len * sizeof(int)));
    int* O = *PO;

    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_int_distribution<> dist(0, M-1);
    
    for (size_t i = 0; i < len; ++i)  
        O[i] = dist(gen);
}

void generate_m(double ** PM, size_t row, size_t col) {   
    *PM = static_cast<double *>(malloc(row*col * sizeof(double)));
    double* M = *PM;

    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<double> dist(0., 1.);

    for (size_t i = 0; i < row; ++i) {
        double sum = 0.;
        for(size_t j = 0; j < col; ++j) {
            M[i*col + j] = dist(gen);
            sum += M[i*col + j];
        }
        assert(sum > 0.);
        for(size_t j = 0; j < col; ++j) {
            M[i*col + j] /= sum;
        }
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