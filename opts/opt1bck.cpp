#include <bits/stdc++.h>
#include "common.h"

size_t flop_count(int N, int M, int T, int n_iter){
    size_t add = 0;
    size_t mul = 0;
    size_t div = 0;

    // forward vars
    mul += N; // FW[i*T + 0] = PI[i]*B[i*M + O[0]];
    add += N; // scale += FW[i*T + 0];

    div += 1; // C[0] = 1./scale;
    mul += N; // FW[i*T + 0]*= C[0];

    add += (T-1)*N*N; // FW[i*T + t] += FW[j*T + t-1]*A[j*N + i]*B[i*M + O[t]];
    mul += 2*(T-1)*N*N;
    add += N*(T-1); // scale += FW[i*T + t];
    div += T-1; // C[t] = 1./scale;
    mul += (T-1)*N; // FW[i*T + t]*= C[t];

    // backward vars
    add += (T-1)*N*N; // BW[i*T + t] += BW[j*T + t+1]*A[i*N+j]*B[j*M + O[t+1]];
    mul += 2*(T-1)*N*N;
    mul += (T-1)*N; // BW[i*T + t]*= C[t];

    // update initial
    mul += N; // PI[i] = FW[i*T + 0]*BW[i*T+0]/C[0];
    div += N;

    // update transition
    add += (T-1)*N*N; // num += FW[i*T + t]*A[i*N + j]*B[j*M + O[t+1]]*BW[j*T + t+1];
    mul += 3*(T-1)*N*N;
    add += (T-1)*N*N; // denom += FW[i*T + t]*BW[i*T + t]/C[t];
    mul += (T-1)*N*N;
    div += (T-1)*N*N;
    div += N*N; // A[i*N + j] = num/denom;

    // update emission
    add += N*M*T; // denom += FW[j*T + t]*BW[j*T + t]/C[t];
    mul += N*M*T;
    div += N*M*T;

    add += N*T; // if(O[t] == o) sum_o += FW[j*T + t]*BW[j*T + t]/C[t];
    mul += N*T;
    div += N*T;

    div += N*M; // B[j*M + o] = sum_o/denom;
    return n_iter*(add + mul + div);
}

void baum_welch(double* PI, double* A, double* B, int* O, double* FW, double* BW, double* C, int N, int M, int T, int n_iter) {
    double* scales = C + T;

    REGION_BEGIN(baum_welch)

    for(int it = 0; it < n_iter; it++) {

        // Calculate the Forward trellis (scaled)
        REGION_BEGIN(forward_vars)

        double scale = 0.;
        int o0 = O[0];
        for(int i = 0; i < N; i++) {
            double fwit = PI[i]*B[i*M + o0];
            FW[i*T] = fwit;
            scale += fwit;
        }

        double c0 = 1./scale;
        C[0] = c0;
        double scales0 = scale;

        for(int i = 0; i < N; i++) {
            FW[i*T] *= c0;
        }

        for(int t = 1; t < T; t++) {
            scale = 0.;
            int ot = O[t];
            for(int i = 0; i < N; i++) {
                double fwitt = 0.;
                double bimot = B[i*M + ot];

                for(int j = 0; j < N; j++) {
                    fwitt += FW[j*T + t-1]*A[j*N + i];
                }
                fwitt *= bimot;

                FW[i*T + t] = fwitt;
                scale += fwitt;
            }

            double ct = 1./scale;
            C[t] = ct;
            scales[t] = scale;

            for(int i = 0; i < N; i++) {
                FW[i*T + t] *= ct;
            }
        }

        REGION_END(forward_vars)

        REGION_BEGIN(backward_vars)
        // Calculate the Backward trellis (scaled)
        double ct1 = C[T-1];
        for(int i = 0; i < N; i++) {
            BW[i*T + T-1] = ct1;
        }

        for(int t = T-2; t >= 0; t--) {
            int ot1 = O[t+1];
            for(int i = 0; i < N; i++) {
                double bwitt = 0.;
                for(int j = 0; j < N; j++) {
                    bwitt += BW[j*T + t+1]*A[i*N+j]*B[j*M + ot1];
                }
                bwitt *= C[t];
                BW[i*T + t] = bwitt;
            }
        }
        REGION_END(backward_vars)

        REGION_BEGIN(update_initial)
        REGION_END(update_initial)

        REGION_BEGIN(update_transition)
        // update the State Transition probabilities

        for(int i = 0; i < N; i++) {

            PI[i] = FW[i*T]*BW[i*T]*scales0;

            for(int j = 0; j < N; j++) {

                double num = 0.;
                double denom = 0.;

                double ainj = A[i*N + j];

                for(int t = 0; t < T-1; t++) {
                    double fwitt = FW[i*T + t];
                    num += fwitt * ainj * B[j*M + O[t+1]] * BW[j*T + t+1];
                    denom += fwitt * BW[i*T + t] * scales[t];
                }

                A[i*N + j] = num/denom;
            }
        }
        REGION_END(update_transition)

        REGION_BEGIN(update_emission)
        // update the State Emission probabilities
        for(int o = 0; o < M; o++) {
            for(int j = 0; j < N; j++) {
                double sum_o = 0.;
                double denom = 0.;
                for(int t = 0; t < T; t++) {
                    denom += FW[j*T + t]*BW[j*T + t]*scales[t];
                    if(O[t] == o) sum_o += FW[j*T + t]*BW[j*T + t]*scales[t];
                }
                B[j*M + o] = sum_o/denom;
            }
        }

        REGION_END(update_emission)
    }

    REGION_END(baum_welch)
}
