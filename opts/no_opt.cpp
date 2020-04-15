#include <bits/stdc++.h>
#include <papi.h>
#include "common.h"

size_t flop_count(int N, int M, int T, int n_iter){
    size_t add = 4*(T-1)*N*N + T*N*M + 2*T*N;
    size_t mult = 8*(T-1)*N*N + T*N*M + 3*T*N + N;
    size_t div = (T-1)*N + N*N + T*N*M + T*N + N*M + N + T;
    return n_iter*(add + mult + div);
}

void baum_welch(double* PI, double* A, double* B, int* O, double* FW, double* BW, double* C, int N, int M, int T, int n_iter) {
    int retval;
    retval = PAPI_hl_region_begin("BaumWelch");
    if ( retval != PAPI_OK ) handle_error(retval);

    for(int it = 0; it < n_iter; it++) {

        // Calculate the Forward trellis (scaled)
        retval = PAPI_hl_region_begin("forwardVars");
        if ( retval != PAPI_OK ) handle_error(retval);
        double scale = 0.;
        for(int i = 0; i < N; i++) {
            FW[i*T + 0] = PI[i]*B[i*M + O[0]];
            scale += FW[i*T + 0];
        }
        // Scale timestep 0
        C[0] = 1./scale;
        for(int i = 0; i < N; i++) {
            FW[i*T + 0]*= C[0];
        }

        for(int t = 1; t < T; t++) {
            scale = 0.;
            for(int i = 0; i < N; i++) {
                FW[i*T + t] = 0.;
                for(int j = 0; j < N; j++) {
                    FW[i*T + t] += FW[j*T + t-1]*A[j*N + i]*B[i*M + O[t]];
                }
                scale += FW[i*T + t];
            }
            C[t] = 1./scale;

            for(int i = 0; i < N; i++) {
                FW[i*T + t]*= C[t];
            }
        }

        retval = PAPI_hl_region_end("forwardVars");
        if ( retval != PAPI_OK ) handle_error(retval);
        
        retval = PAPI_hl_region_begin("backwardVars");
        if ( retval != PAPI_OK ) handle_error(retval);
        // Calculate the Backward trellis (scaled)
        for(int i = 0; i < N; i++) {
            BW[i*T + T-1] = C[T-1];
        }

        for(int t = T-2; t >= 0; t--) {
            for(int i = 0; i < N; i++) {
                BW[i*T + t] = 0.;
                for(int j = 0; j < N; j++) {
                    BW[i*T + t] += BW[j*T + t+1]*A[i*N+j]*B[j*M + O[t+1]];
                }
                BW[i*T + t]*= C[t];
            }
        }
        retval = PAPI_hl_region_end("backwardVars");
        if ( retval != PAPI_OK ) handle_error(retval);

        retval = PAPI_hl_region_begin("updateInitial");
        if ( retval != PAPI_OK ) handle_error(retval);
        // update the Initial State probabilities
        for(int i = 0; i < N; i++) {
            PI[i] = FW[i*T + 0]*BW[i*T+0]/C[0];
        }
        retval = PAPI_hl_region_end("updateInitial");
        if ( retval != PAPI_OK ) handle_error(retval);

        retval = PAPI_hl_region_begin("updateTransition");
        if ( retval != PAPI_OK ) handle_error(retval);
        // update the State Transition probabilities
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                double num = 0.;
                double denom = 0.;
                for(int t = 0; t < T-1; t++) {
                    num += FW[i*T + t]*A[i*N + j]*B[j*M + O[t+1]]*BW[j*T + t+1];
                    denom += FW[i*T + t]*BW[i*T + t]/C[t];
                }
                A[i*N + j] = num/denom;
            }
        }
        retval = PAPI_hl_region_end("updateTransition");
        if ( retval != PAPI_OK ) handle_error(retval);

        retval = PAPI_hl_region_begin("updateEmission");
        if ( retval != PAPI_OK ) handle_error(retval);
        // update the State Emission probabilities
        for(int o = 0; o < M; o++) {
            for(int j = 0; j < N; j++) {
                double sum_o = 0.;
                double denom = 0.;
                for(int t = 0; t < T; t++) {
                    denom += FW[j*T + t]*BW[j*T + t]/C[t];
                    if(O[t] == o) sum_o += FW[j*T + t]*BW[j*T + t]/C[t];;
                }
                B[j*M + o] = sum_o/denom;
            }
        }

        retval = PAPI_hl_region_end("updateEmission");
        if ( retval != PAPI_OK ) handle_error(retval);
    }

    retval = PAPI_hl_region_end("BaumWelch");
    if ( retval != PAPI_OK ) handle_error(retval);
}
