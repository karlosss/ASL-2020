#ifndef ASL_PROJECT_COMMON_H
#define ASL_PROJECT_COMMON_H

#include <bits/stdc++.h>
#include <papi.h>

#define REGION_BEGIN(name) for (int papi_retval = PAPI_hl_region_begin(#name); papi_retval != PAPI_OK; ) {handle_error(papi_retval); break;}
#define REGION_END(name) for (int papi_retval = PAPI_hl_region_end(#name); papi_retval != PAPI_OK; ) {handle_error(papi_retval); break;}

void generate_input(size_t N, size_t M, size_t T, int** O, double** PI, double** A, double** B, double** FW, double** BW, double** C);

void clone_input(size_t N, size_t M, size_t T, int* O, double* PI, double* A, double* B,
        int** O2, double** PI2, double** A2, double** B2, double** FW2, double** BW2, double** C2);

bool compare_outputs(size_t N, size_t M, double* PI, double* A, double* B, 
                                                 double* PI2, double* A2, double* B2);

void baum_welch(double* PI, double* A, double* B, int* O, double* FW, double* BW, double* C, int N, int M, int T, int n_iter);

void baum_welch_base(double* PI, double* A, double* B, int* O, double* FW, double* BW, double* C, int N, int M, int T, int n_iter);

size_t flop_count(int N, int M, int T, int n_iter);

void handle_error (int retval);

#endif //ASL_PROJECT_COMMON_H
