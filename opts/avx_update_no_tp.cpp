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
    double* sum_os = scales + T;
    double* denoms = sum_os + M*N;
    init_zero(sum_os, M*N);

    transpose(N,M,B);

    REGION_BEGIN(baum_welch)

    for(int it = 0; it < n_iter; it++) {

        // Calculate the Forward trellis (scaled)
        REGION_BEGIN(forward_vars)

        double scale0 = 0.;
        double scale1 = 0.;
        double scale2 = 0.;
        double scale3 = 0.;

        int o0 = O[0];
        for(int i = 0; i < N; i+=4) {
            double fwit0 = PI[i]*B[i*M + o0];
            double fwit1 = PI[i+1]*B[(i+1)*M + o0];
            double fwit2 = PI[i+2]*B[(i+2)*M + o0];
            double fwit3 = PI[i+3]*B[(i+3)*M + o0];

            FW[i*T] = fwit0;
            FW[(i+1)*T] = fwit1;
            FW[(i+2)*T] = fwit2;
            FW[(i+3)*T] = fwit3;

            scale0 += fwit0;
            scale1 += fwit1;
            scale2 += fwit2;
            scale3 += fwit3;
        }

        double scale = scale0 + scale1 + scale2 + scale3;

        double c0 = 1./scale;
        C[0] = c0;
        double scales0 = scale;
        scales[0] = scale;

        for(int i = 0; i < N; i += 4) {
            FW[i*T] *= c0;
            FW[(i+1)*T] *= c0;
            FW[(i+2)*T] *= c0;
            FW[(i+3)*T] *= c0;
        }

        for(int t = 1; t < T; t++) {
            scale = 0.;
            int ot = O[t];
            for(int i = 0; i < N; i++) {
                double fwitt0 = 0.;
                double fwitt1 = 0.;
                double fwitt2 = 0.;
                double fwitt3 = 0.;

                for(int j = 0; j < N; j += 4) {
                    fwitt0 += FW[j*T + t-1]*A[j*N + i];
                    fwitt1 += FW[(j+1)*T + t-1]*A[(j+1)*N + i];
                    fwitt2 += FW[(j+2)*T + t-1]*A[(j+2)*N + i];
                    fwitt3 += FW[(j+3)*T + t-1]*A[(j+3)*N + i];
                }

                double fwitt = (fwitt0+fwitt1+fwitt2+fwitt3) * B[i*M + ot];

                FW[i*T + t] = fwitt;
                scale += fwitt;
            }

            double ct = 1./scale;
            C[t] = ct;
            scales[t] = scale;

            for(int i = 0; i < N; i += 4) {
                FW[i*T + t] *= ct;
                FW[(i+1)*T + t] *= ct;
                FW[(i+2)*T + t] *= ct;
                FW[(i+3)*T + t] *= ct;
            }
        }

        REGION_END(forward_vars)

        REGION_BEGIN(backward_vars)
        // Calculate the Backward trellis (scaled)
        double ct1 = C[T-1];

        for(int i = 0; i < N; i+=4) {
            BW[i*T + T-1] = ct1;
            BW[(i+1)*T + T-1] = ct1;
            BW[(i+2)*T + T-1] = ct1;
            BW[(i+3)*T + T-1] = ct1;
        }

        for(int t = T-2; t >= 0; t--) {
            int ot1 = O[t+1];
            double ct = C[t];

            for(int i = 0; i < N; i++) {
                double bwitt0 = 0.;
                double bwitt1 = 0.;
                double bwitt2 = 0.;
                double bwitt3 = 0.;
                for(int j = 0; j < N; j+=4) {
                    bwitt0 += BW[j*T + t+1]*A[i*N+j]*B[j*M + ot1];
                    bwitt1 += BW[(j+1)*T + t+1]*A[i*N+j+1]*B[(j+1)*M + ot1];
                    bwitt2 += BW[(j+2)*T + t+1]*A[i*N+j+2]*B[(j+2)*M + ot1];
                    bwitt3 += BW[(j+3)*T + t+1]*A[i*N+j+3]*B[(j+3)*M + ot1];
                }
                double bwitt = (bwitt0+bwitt1+bwitt2+bwitt3)*ct;
                BW[i*T + t] = bwitt;
            }
        }
        REGION_END(backward_vars)

        REGION_BEGIN(update_initial)
        REGION_END(update_initial)

        REGION_BEGIN(update_transition)
        // update the State Transition probabilities

        double scalest1 = scales[T-1];
        int ot1 = O[T-1];
        int ot2 = O[T-2];
        int ot3 = O[T-3];

        for(int i = 0; i < N; i++) {

            double pi = FW[i*T] * BW[i*T] * scales0;

            sum_os[i*M + o0] = pi;
            PI[i] = pi;

            double denom0 = 0;
            double denom1 = 0;
            double denom2 = 0;
            double denom3 = 0;

            for(int t = 1; t < T-3; t+=4) {
                double toadd0 = FW[i*T + t] * BW[i*T + t] * scales[t];
                double toadd1 = FW[i*T + t+1] * BW[i*T + t+1] * scales[t+1];
                double toadd2 = FW[i*T + t+2] * BW[i*T + t+2] * scales[t+2];
                double toadd3 = FW[i*T + t+3] * BW[i*T + t+3] * scales[t+3];

                denom0 += toadd0;
                denom1 += toadd1;
                denom2 += toadd2;
                denom3 += toadd3;

                sum_os[i*M + O[t]] += toadd0;
                sum_os[i*M + O[t+1]] += toadd1;
                sum_os[i*M + O[t+2]] += toadd2;
                sum_os[i*M + O[t+3]] += toadd3;
            }

            // leftover 2 iterations
            double toadd0 = FW[i*T + T-3] * BW[i*T + T-3] * scales[T-3];
            double toadd1 = FW[i*T + T-2] * BW[i*T + T-2] * scales[T-2];

            denom0 += toadd0;
            denom1 += toadd1;

            sum_os[i*M + ot3] += toadd0;
            sum_os[i*M + ot2] += toadd1;

            double denom = denom0+denom1+denom2+denom3+pi;

            double lastadd = FW[i*T + T-1] * BW[i*T + T-1] * scalest1;
            denoms[i] = denom + lastadd;
            sum_os[i*M + ot1] += lastadd;

            for(int j = 0; j < N; j++) {
                double num0 = 0.;
                double num1 = 0.;
                double num2 = 0.;
                double num3 = 0.;

                for(int t = 0; t < T-4; t+=4) {
                    num0 += FW[i*T + t] * B[j*M + O[t+1]] * BW[j*T + t+1];
                    num1 += FW[i*T + t+1] * B[j*M + O[t+2]] * BW[j*T + t+2];
                    num2 += FW[i*T + t+2] * B[j*M + O[t+3]] * BW[j*T + t+3];
                    num3 += FW[i*T + t+3] * B[j*M + O[t+4]] * BW[j*T + t+4];
                }

                num0 += FW[i*T + T-4] * B[j*M + ot3] * BW[j*T + T-3];
                num1 += FW[i*T + T-3] * B[j*M + ot2] * BW[j*T + T-2];
                num2 += FW[i*T + T-2] * B[j*M + ot1] * BW[j*T + T-1];

                double num = (num0+num1+num2+num3)*A[i*N + j];
                A[i*N + j] = num/denom;
            }
        }
        REGION_END(update_transition)

        REGION_BEGIN(update_emission)
        // update the State Emission probabilities
        for(int i = 0; i < N; i++) {
            double denomsi = denoms[i];
            for(int o = 0; o < M; o += 4){
                B[i*M + o] = sum_os[i*M + o]/denomsi;
                B[i*M + o+1] = sum_os[i*M + o+1]/denomsi;
                B[i*M + o+2] = sum_os[i*M + o+2]/denomsi;
                B[i*M + o+3] = sum_os[i*M + o+3]/denomsi;

                sum_os[i*M + o] = 0;
                sum_os[i*M + o+1] = 0;
                sum_os[i*M + o+2] = 0;
                sum_os[i*M + o+3] = 0;
            }
        }
        REGION_END(update_emission)
    }

    REGION_END(baum_welch)

    transpose(M,N,B);
}
