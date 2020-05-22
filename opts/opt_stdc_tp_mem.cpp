#include <bits/stdc++.h>
#include "common.h"
#include <immintrin.h>

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
    double* sum_os = scales + T;
    double* denoms = sum_os + M*N;
    init_zero(sum_os, M*N);

    transpose(N, M, B);

    REGION_BEGIN(baum_welch)

    for(int it = 0; it < n_iter; it++) {
        
        // Calculate the Forward trellis (scaled)
        REGION_BEGIN(forward_vars)
        double scale = 0.;
        int o0 = O[0];
        for(int i = 0; i < N; i++) {
            FW[0*N + i] = PI[i]*B[o0*N + i];
            scale += FW[0*N + i];
        }
        // Scale timestep 0
        C[0] = 1./scale;
        scales[0] = scale;
        double scales0 = scale;
        for(int i = 0; i < N; i++) {
            FW[0*N + i]*= C[0];
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
                    double fw_j = FW[(t-1)*N + j];
                    accum0 += fw_j*A[j*N + i];
                    accum1 += fw_j*A[j*N + i+1];
                    accum2 += fw_j*A[j*N + i+2];
                    accum3 += fw_j*A[j*N + i+3];

                    double fw_j1 = FW[(t-1)*N + j+1];
                    accum4 += fw_j1*A[(j+1)*N + i];
                    accum5 += fw_j1*A[(j+1)*N + i+1];
                    accum6 += fw_j1*A[(j+1)*N + i+2];
                    accum7 += fw_j1*A[(j+1)*N + i+3];
                }
                FW[t*N + i] = B[obs*N + i]*(accum0+accum4);
                FW[t*N + i+1] = B[obs*N + i+1]*(accum1+accum5);
                FW[t*N + i+2] = B[obs*N + i+2]*(accum2+accum6);
                FW[t*N + i+3] = B[obs*N + i+3]*(accum3+accum7);
                scale += FW[t*N + i] + FW[t*N + i+1] + FW[t*N + i+2] + FW[t*N + i+3];
            }
            C[t] = 1./scale;
            scales[t] = scale;
            double c_t = C[t];
            for(int i = 0; i < N; i++) {
                FW[t*N + i]*= c_t;
            }
        }

        REGION_END(forward_vars)

        REGION_BEGIN(backward_vars)
        // Calculate the Backward trellis (scaled)
        double bw_c_t_1 = C[T-1];
        for(int i = 0; i < N; i++) {
            BW[(T-1)*N + i] = bw_c_t_1;
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
                    double bw_j = BW[(t+1)*N + j];
                    double b_j = B[obs*N + j];
                    accum0 += bw_j*A[i*N+j]*b_j;
                    accum1 += bw_j*A[(i+1)*N+j]*b_j;
                    accum2 += bw_j*A[(i+2)*N+j]*b_j;
                    accum3 += bw_j*A[(i+3)*N+j]*b_j;

                }

                BW[t*N + i] = c_t*accum0;
                BW[t*N+i+1] = c_t*accum1;
                BW[t*N+i+2] = c_t*accum2;
                BW[t*N+i+3] = c_t*accum3;
            }
        }
        REGION_END(backward_vars)

        REGION_BEGIN(update_initial)
        REGION_END(update_initial)

        REGION_BEGIN(update_transition)
        // update the State Transition probabilities

        double scalest1 = scales[T-1];
        int ot1 = O[T-1];

        for(int i = 0; i < N; i+=4) {

            double pi0 = FW[i+0] * BW[i+0] * scales0;
            double pi1 = FW[i+1] * BW[i+1] * scales0;
            double pi2 = FW[i+2] * BW[i+2] * scales0;
            double pi3 = FW[i+3] * BW[i+3] * scales0;

            sum_os[(i+0) + o0*N] = pi0;
            sum_os[(i+1) + o0*N] = pi1;
            sum_os[(i+2) + o0*N] = pi2;
            sum_os[(i+3) + o0*N] = pi3;

            PI[i+0] = pi0;
            PI[i+1] = pi1;
            PI[i+2] = pi2;
            PI[i+3] = pi3;

            double denomi0 = 0.0;
            double denomi1 = 0.0;
            double denomi2 = 0.0;
            double denomi3 = 0.0;
            for(int t = 1; t < T-1; t+=1) {
                double toaddi0 = FW[(i+0) + t*N] * BW[(i+0) + t*N] * scales[t];
                double toaddi1 = FW[(i+1) + t*N] * BW[(i+1) + t*N] * scales[t];
                double toaddi2 = FW[(i+2) + t*N] * BW[(i+2) + t*N] * scales[t];
                double toaddi3 = FW[(i+3) + t*N] * BW[(i+3) + t*N] * scales[t];

                denomi0 += toaddi0;
                denomi1 += toaddi1;
                denomi2 += toaddi2;
                denomi3 += toaddi3;

                sum_os[(i+0) + O[t]*N] += toaddi0;
                sum_os[(i+1) + O[t]*N] += toaddi1;
                sum_os[(i+2) + O[t]*N] += toaddi2;
                sum_os[(i+3) + O[t]*N] += toaddi3;
            }
            denomi0 += pi0;
            denomi1 += pi1;
            denomi2 += pi2;
            denomi3 += pi3;

            double lastaddi0 = FW[(i+0) + (T-1)*N] * BW[(i+0) + (T-1)*N] * scalest1;
            double lastaddi1 = FW[(i+1) + (T-1)*N] * BW[(i+1) + (T-1)*N] * scalest1;
            double lastaddi2 = FW[(i+2) + (T-1)*N] * BW[(i+2) + (T-1)*N] * scalest1;
            double lastaddi3 = FW[(i+3) + (T-1)*N] * BW[(i+3) + (T-1)*N] * scalest1;

            denoms[i+0] = denomi0 + lastaddi0;
            denoms[i+1] = denomi1 + lastaddi1;
            denoms[i+2] = denomi2 + lastaddi2;
            denoms[i+3] = denomi3 + lastaddi3;

            sum_os[(i+0) + ot1*N] += lastaddi0;
            sum_os[(i+1) + ot1*N] += lastaddi1;
            sum_os[(i+2) + ot1*N] += lastaddi2;
            sum_os[(i+3) + ot1*N] += lastaddi3;

            for(int j = 0; j < N; j += 4) {
                double numi0j0 = 0.0;
                double numi1j0 = 0.0;
                double numi2j0 = 0.0;
                double numi3j0 = 0.0;

                double numi0j1 = 0.0;
                double numi1j1 = 0.0;
                double numi2j1 = 0.0;
                double numi3j1 = 0.0;

                double numi0j2 = 0.0;
                double numi1j2 = 0.0;
                double numi2j2 = 0.0;
                double numi3j2 = 0.0;

                double numi0j3 = 0.0;
                double numi1j3 = 0.0;
                double numi2j3 = 0.0;
                double numi3j3 = 0.0;

                for(int t = 0; t < T-1; t+=1) {
                    numi0j0 += FW[(i+0) + (t  )*N] * B[(j+0) + O[t+1]*N] * BW[(j+0) + (t+1)*N];
                    numi1j0 += FW[(i+1) + (t  )*N] * B[(j+0) + O[t+1]*N] * BW[(j+0) + (t+1)*N];
                    numi2j0 += FW[(i+2) + (t  )*N] * B[(j+0) + O[t+1]*N] * BW[(j+0) + (t+1)*N];
                    numi3j0 += FW[(i+3) + (t  )*N] * B[(j+0) + O[t+1]*N] * BW[(j+0) + (t+1)*N];

                    numi0j1 += FW[(i+0) + (t  )*N] * B[(j+1) + O[t+1]*N] * BW[(j+1) + (t+1)*N];
                    numi1j1 += FW[(i+1) + (t  )*N] * B[(j+1) + O[t+1]*N] * BW[(j+1) + (t+1)*N];
                    numi2j1 += FW[(i+2) + (t  )*N] * B[(j+1) + O[t+1]*N] * BW[(j+1) + (t+1)*N];
                    numi3j1 += FW[(i+3) + (t  )*N] * B[(j+1) + O[t+1]*N] * BW[(j+1) + (t+1)*N];

                    numi0j2 += FW[(i+0) + (t  )*N] * B[(j+2) + O[t+1]*N] * BW[(j+2) + (t+1)*N];
                    numi1j2 += FW[(i+1) + (t  )*N] * B[(j+2) + O[t+1]*N] * BW[(j+2) + (t+1)*N];
                    numi2j2 += FW[(i+2) + (t  )*N] * B[(j+2) + O[t+1]*N] * BW[(j+2) + (t+1)*N];
                    numi3j2 += FW[(i+3) + (t  )*N] * B[(j+2) + O[t+1]*N] * BW[(j+2) + (t+1)*N];

                    numi0j3 += FW[(i+0) + (t  )*N] * B[(j+3) + O[t+1]*N] * BW[(j+3) + (t+1)*N];
                    numi1j3 += FW[(i+1) + (t  )*N] * B[(j+3) + O[t+1]*N] * BW[(j+3) + (t+1)*N];
                    numi2j3 += FW[(i+2) + (t  )*N] * B[(j+3) + O[t+1]*N] * BW[(j+3) + (t+1)*N];
                    numi3j3 += FW[(i+3) + (t  )*N] * B[(j+3) + O[t+1]*N] * BW[(j+3) + (t+1)*N];
                }
                A[(i+0)*N + (j+0)] = numi0j0*A[(i+0)*N + (j+0)]/denomi0;
                A[(i+1)*N + (j+0)] = numi1j0*A[(i+1)*N + (j+0)]/denomi1;
                A[(i+2)*N + (j+0)] = numi2j0*A[(i+2)*N + (j+0)]/denomi2;
                A[(i+3)*N + (j+0)] = numi3j0*A[(i+3)*N + (j+0)]/denomi3;

                A[(i+0)*N + (j+1)] = numi0j1*A[(i+0)*N + (j+1)]/denomi0;
                A[(i+1)*N + (j+1)] = numi1j1*A[(i+1)*N + (j+1)]/denomi1;
                A[(i+2)*N + (j+1)] = numi2j1*A[(i+2)*N + (j+1)]/denomi2;
                A[(i+3)*N + (j+1)] = numi3j1*A[(i+3)*N + (j+1)]/denomi3;

                A[(i+0)*N + (j+2)] = numi0j2*A[(i+0)*N + (j+2)]/denomi0;
                A[(i+1)*N + (j+2)] = numi1j2*A[(i+1)*N + (j+2)]/denomi1;
                A[(i+2)*N + (j+2)] = numi2j2*A[(i+2)*N + (j+2)]/denomi2;
                A[(i+3)*N + (j+2)] = numi3j2*A[(i+3)*N + (j+2)]/denomi3;

                A[(i+0)*N + (j+3)] = numi0j3*A[(i+0)*N + (j+3)]/denomi0;
                A[(i+1)*N + (j+3)] = numi1j3*A[(i+1)*N + (j+3)]/denomi1;
                A[(i+2)*N + (j+3)] = numi2j3*A[(i+2)*N + (j+3)]/denomi2;
                A[(i+3)*N + (j+3)] = numi3j3*A[(i+3)*N + (j+3)]/denomi3;
            }
        }
        REGION_END(update_transition)

        REGION_BEGIN(update_emission)

        for(int o = 0; o < M; o += 1){
            for(int i = 0; i < N-3; i+=4) {
                double denomsi = denoms[i];
                double denomsi1 = denoms[i+1];
                double denomsi2 = denoms[i+2];
                double denomsi3 = denoms[i+3];
                B[o*N + i] = sum_os[o*N + i]/denomsi;
                B[o*N + i+1] = sum_os[o*N + i+1]/denomsi1;
                B[o*N + i+2] = sum_os[o*N + i+2]/denomsi2;
                B[o*N + i+3] = sum_os[o*N + i+3]/denomsi3;
                sum_os[o*N + i] = 0;
                sum_os[o*N + i+1] = 0;
                sum_os[o*N + i+2] = 0;
                sum_os[o*N + i+3] = 0;
            }
        }
        REGION_END(update_emission)
    }

    REGION_END(baum_welch)
    transpose(M, N, B);
}
