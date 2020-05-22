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
    double* scales = (double*) malloc(T*sizeof(double));
    
    REGION_BEGIN(baum_welch)

    for(int it = 0; it < n_iter; it++) {

        // Calculate the Forward trellis (scaled)
        REGION_BEGIN(forward_vars)
        double scale = 0.;
        for(int i = 0; i < N; i++) {
            FW[i*T + 0] = PI[i]*B[i*M + O[0]];
            scale += FW[i*T + 0];
        }
        // Scale timestep 0
        C[0] = 1./scale;
        scales[0] = scale;
        for(int i = 0; i < N; i++) {
            FW[i*T + 0]*= C[0];
        }

        for(int t = 1; t < T; t++) {
            scale = 0.;
            int obs = O[t];
            for(int i = 0; i < N-3; i+=4) {
                double accum0 = 0.;
                double accum1 = 0.;
                double accum2 = 0.;
                double accum3 = 0.;

                double accum4 = 0.;
                double accum5 = 0.;
                double accum6 = 0.;
                double accum7 = 0.;

                

                for(int j = 0; j < N-1; j+=2) {
                    double fw_j = FW[j*T + t-1];
                    accum0 += fw_j*A[j*N + i];
                    accum1 += fw_j*A[j*N + i+1];
                    accum2 += fw_j*A[j*N + i+2];
                    accum3 += fw_j*A[j*N + i+3];

                    double fw_j1 = FW[(j+1)*T + t-1];
                    accum4 += fw_j1*A[(j+1)*N + i];
                    accum5 += fw_j1*A[(j+1)*N + i+1];
                    accum6 += fw_j1*A[(j+1)*N + i+2];
                    accum7 += fw_j1*A[(j+1)*N + i+3];
                }
                FW[i*T + t] = B[i*M + obs]*(accum0+accum4);
                FW[(i+1)*T + t] = B[(i+1)*M + obs]*(accum1+accum5);
                FW[(i+2)*T + t] = B[(i+2)*M + obs]*(accum2+accum6);
                FW[(i+3)*T + t] = B[(i+3)*M + obs]*(accum3+accum7);
                scale += FW[i*T + t] + FW[(i+1)*T + t] + FW[(i+2)*T + t] + FW[(i+3)*T + t];
            }
            C[t] = 1./scale;
            scales[t] = scale;
            double c_t = C[t];
            for(int i = 0; i < N; i++) {
                FW[i*T + t]*= c_t;
            }
        }

        REGION_END(forward_vars)

        REGION_BEGIN(backward_vars)
        // Calculate the Backward trellis (scaled)
        double bw_c_t_1 = C[T-1];
        for(int i = 0; i < N; i++) {
            BW[i*T + T-1] = bw_c_t_1;
        }

        for(int t = T-2; t >= 0; t--) {
            int obs = O[t+1];
            double c_t = C[t];
            for(int i = 0; i < N-3; i+=4) {
                double accum0 = 0.;
                double accum1 = 0.;
                double accum2 = 0.;
                double accum3 = 0.;
                
                for(int j = 0; j < N; j++) {
                    double bw_j = BW[j*T + t+1];
                    double b_j = B[j*M + obs];
                    accum0 += bw_j*A[i*N+j]*b_j;
                    accum1 += bw_j*A[(i+1)*N+j]*b_j;
                    accum2 += bw_j*A[(i+2)*N+j]*b_j;
                    accum3 += bw_j*A[(i+3)*N+j]*b_j;

                }

                BW[i*T + t] = c_t*accum0;
                BW[(i+1)*T + t] = c_t*accum1;
                BW[(i+2)*T + t] = c_t*accum2;
                BW[(i+3)*T + t] = c_t*accum3;
            }
        }
        REGION_END(backward_vars)

        REGION_BEGIN(update_initial)
        // update the Initial State probabilities
        for(int i = 0; i < N; i++) {
            PI[i] = FW[i*T + 0]*BW[i*T+0]*scales[0];
        }
        REGION_END(update_initial)

        REGION_BEGIN(update_transition)
        // update the State Transition probabilities
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
            
                double n_accum0 = 0.;
                double n_accum1 = 0.;
                double n_accum2 = 0.;
                double n_accum3 = 0.;

                double d_accum0 = 0.;
                double d_accum1 = 0.;
                double d_accum2 = 0.;
                double d_accum3 = 0.;

                int t = 0;
                for(; t < T-4; t+=4) {
                    n_accum0 += FW[i*T + t]*A[i*N + j]*B[j*M + O[t+1]]*BW[j*T + t+1];
                    n_accum1 += FW[i*T + t+1]*A[i*N + j]*B[j*M + O[t+2]]*BW[j*T + t+2];
                    n_accum2 += FW[i*T + t+2]*A[i*N + j]*B[j*M + O[t+3]]*BW[j*T + t+3];
                    n_accum3 += FW[i*T + t+3]*A[i*N + j]*B[j*M + O[t+4]]*BW[j*T + t+4];

                    d_accum0 += FW[i*T + t]*BW[i*T + t]*scales[t];
                    d_accum1 += FW[i*T + t+1]*BW[i*T + t+1]*scales[t+1];
                    d_accum2 += FW[i*T + t+2]*BW[i*T + t+2]*scales[t+2];
                    d_accum3 += FW[i*T + t+3]*BW[i*T + t+3]*scales[t+3];
                }
                for(; t < T-1; t++) {
                    n_accum0 += FW[i*T + t]*A[i*N + j]*B[j*M + O[t+1]]*BW[j*T + t+1];
                    d_accum0 += FW[i*T + t]*BW[i*T + t]*scales[t]; 
                }
                A[i*N + j] = (n_accum0 + n_accum1 + n_accum2 + n_accum3)/(d_accum0 + d_accum1 + d_accum2 + d_accum3);
            }
        }
        REGION_END(update_transition)

        REGION_BEGIN(update_emission)
        // update the State Emission probabilities
        for(int o = 0; o < M; o++) {
            for(int j = 0; j < N; j++) {
                double n_accum0 = 0.;
                double n_accum1 = 0.;
                double n_accum2 = 0.;
                double n_accum3 = 0.;

                double d_accum0 = 0.;
                double d_accum1 = 0.;
                double d_accum2 = 0.;
                double d_accum3 = 0.;
                for(int t = 0; t < T-3; t+=4) {
                    double mul0 = FW[j*T + t]*BW[j*T + t]*scales[t];
                    d_accum0 += mul0;
                    double mul1 = FW[j*T + t+1]*BW[j*T + t+1]*scales[t+1];
                    d_accum1 += mul1;
                    double mul2 = FW[j*T + t+2]*BW[j*T + t+2]*scales[t+2];
                    d_accum2 += mul2;
                    double mul3 = FW[j*T + t+3]*BW[j*T + t+3]*scales[t+3];
                    d_accum3 += mul3;

                    if(O[t] == o) n_accum0 += mul0;
                    if(O[t+1] == o) n_accum1 += mul1;
                    if(O[t+2] == o) n_accum2 += mul2;
                    if(O[t+3] == o) n_accum3 += mul3;
                }
                B[j*M + o] = (n_accum0+n_accum1+n_accum2+n_accum3)/(d_accum0+d_accum1+d_accum2+d_accum3);
            }
        }

        REGION_END(update_emission)
    }

    REGION_END(baum_welch)
}
