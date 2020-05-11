#include <math.h>
#include <bits/stdc++.h>
#include <papi.h>
#define EPS (1e-3)

void generate_observation(int** O, size_t len, size_t M) {
    //*O = static_cast<int *>(malloc(len * sizeof(int)));

    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_int_distribution<> dist(0, M - 1);

    for (size_t i = 0; i < len; ++i)
        (*O)[i] = dist(gen);
}

void generate_m(double** M, size_t row, size_t col) {
    //*M = static_cast<double *>(malloc(row * col * sizeof(double)));

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
    char* arr = (char*) malloc(
            T*sizeof(int) + // generate_observation(O, T, M);
            1*N*sizeof(double) + // generate_m(PI, 1, N);
            N*N*sizeof(double) + // generate_m(A, N, N);
            N*M*sizeof(double) + // generate_m(B, N, M);
            N*T*sizeof(double) + // *FW = static_cast<double *>(malloc(N * T * sizeof(double)));
            N*T*sizeof(double) + // *BW = static_cast<double *>(malloc(N * T * sizeof(double)));
            T*sizeof(double)     // *C = static_cast<double *>(malloc(T * sizeof(double)));
    );

    size_t ptr = 0;

    *O = (int*) arr;
    ptr += T*sizeof(int);

    *PI = (double*) (arr + ptr);
    ptr += 1*N*sizeof(double);

    *A = (double*) (arr + ptr);
    ptr += N*N*sizeof(double);

    *B = (double*) (arr + ptr);
    ptr += N*M*sizeof(double);

    *FW = (double*) (arr + ptr);
    ptr += N*T*sizeof(double);

    *BW = (double*) (arr + ptr);
    ptr += N*T*sizeof(double);

    *C = (double*) (arr + ptr);

    generate_observation(O, T, M);
    generate_m(PI, 1, N);
    generate_m(A, N, N);
    generate_m(B, N, M);
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


double nrm_sqr_diff(double *x, double *y, int n) {
    double nrm_sqr = 0.0;
    for(int i = 0; i < n; i++) {
        nrm_sqr += (x[i] - y[i]) * (x[i] - y[i]);
    }
    
    if (isnan(nrm_sqr)) {
      nrm_sqr = INFINITY;
    }
    
    return nrm_sqr;
}


void clone_input(size_t N, size_t M, size_t T, int* O, double* PI, double* A, double* B,
                 int** O2, double** PI2, double** A2, double** B2, double** FW2, double** BW2, double** C2){

    generate_input(N, M, T, O2, PI2, A2, B2, FW2, BW2, C2);
    copy_array(O, *O2, T);
    copy_array(PI, *PI2, N);
    copy_array(A, *A2, N*N);
    copy_array(B, *B2, N*M);
}

bool compare_outputs(size_t N, size_t M, double* PI, double* A, double* B,
                     double* PI2, double* A2, double* B2){
    return
    nrm_sqr_diff(PI, PI2, N) < EPS &&
    nrm_sqr_diff(A, A2, N*N) < EPS &&
    nrm_sqr_diff(B, B2, N*M) < EPS;
}

void transpose(size_t N, size_t M, double* X) {
    double* Y = (double*) malloc(N*M*sizeof(double));
    memcpy(Y, X, N*M*sizeof(double));

    for(int i = 0; i < N; i++) {
        for(int j = 0; j < M; j++) {
            X[j*N+i] = Y[i*M + j];
        }
    }
}

void handle_error (int retval){
    std::cout << "PAPI error " <<  retval << ": " << PAPI_strerror(retval) << std::endl;
    exit(1);
}
