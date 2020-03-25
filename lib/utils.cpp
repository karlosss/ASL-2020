#include <bits/stdc++.h>

void generate_observation(int** O, size_t len, size_t M) {
    *O = static_cast<int *>(malloc(len * sizeof(int)));

    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_int_distribution<> dist(0, M - 1);

    for (size_t i = 0; i < len; ++i)
        (*O)[i] = dist(gen);
}

void generate_m(double **M, size_t row, size_t col) {
    *M = static_cast<double *>(malloc(row * col * sizeof(double)));

    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<double> dist(0., 1.);

    for (size_t i = 0; i < row; ++i) {
        double sum = 0.;
        for (size_t j = 0; j < col; ++j) {
            (*M)[i * col + j] = dist(gen);
            sum += (*M)[i * col + j];
        }
        assert(sum > 0.);
        for (size_t j = 0; j < col; ++j) {
            (*M)[i * col + j] /= sum;
        }
    }
}

void generate_input(size_t N, size_t M, size_t T, int** O, double** PI, double** A, double** B, double** FW, double** BW, double** C){
    generate_observation(O, T, M);
    generate_m(PI, 1, N);
    generate_m(A, N, N);
    generate_m(B, N, M);
    *FW = static_cast<double *>(malloc(N * T * sizeof(double)));
    *BW = static_cast<double *>(malloc(N * T * sizeof(double)));
    *C = static_cast<double *>(malloc(T * sizeof(double)));
}


void copy_array(const double* src, double* dest, size_t s){
    for(size_t i = 0; i < s; ++i){
        dest[i] = src[i];
    }
}


void copy_array(const int* src, int* dest, size_t s){
    for(size_t i = 0; i < s; ++i){
        dest[i] = src[i];
    }
}


bool compare_arrays(const double* a, const double* b, size_t s){
    for(size_t i = 0; i < s; ++i){
        if(a[i] != b[i]) return false;
    }
    return true;
}

bool compare_arrays(const int* a, const int* b, size_t s){
    for(size_t i = 0; i < s; ++i){
        if(a[i] != b[i]) return false;
    }
    return true;
}


void clone_input(size_t N, size_t M, size_t T, int* O, double* PI, double* A, double* B,
                 int** O2, double** PI2, double** A2, double** B2, double** FW2, double** BW2, double** C2){

    generate_input(N, M, T, O2, PI2, A2, B2, FW2, BW2, C2);
    copy_array(O, *O2, T);
    copy_array(PI, *PI2, N);
    copy_array(A, *A2, N*N);
    copy_array(B, *B2, N*M);
}

bool compare(size_t N, size_t M, size_t T, int* O, double* PI, double* A, double* B, double* FW, double* BW, double* C,
             int* O2, double* PI2, double* A2, double* B2, double* FW2, double* BW2, double* C2){
    return
    compare_arrays(O, O2, T) &&
    compare_arrays(PI, PI2, N) &&
    compare_arrays(A, A2, N*N) &&
    compare_arrays(B, B2, N*M) &&
    compare_arrays(FW, FW2, N*T) &&
    compare_arrays(BW, BW2, N*T) &&
    compare_arrays(C, C2, T);
}
