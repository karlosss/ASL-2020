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
        __m256d scale_vec0 = _mm256_setzero_pd();
        int o0 = O[0];
        for(int i = 0; i < N-3; i+=4) {
            __m256d pi_vec = _mm256_loadu_pd(PI + i);
            __m256d b_vec = _mm256_loadu_pd(B + O[0]*N + i);
            __m256d result = _mm256_mul_pd(pi_vec, b_vec);
            _mm256_storeu_pd(FW + 0*N + i, result);
            scale_vec0 = _mm256_add_pd(scale_vec0, result);
        }
        __m128d scale_low0 = _mm256_castpd256_pd128(scale_vec0);
        __m128d scale_high0 = _mm256_extractf128_pd(scale_vec0, 1);
        __m128d scale_sum0 = _mm_add_pd(scale_low0, scale_high0);
        __m128d get_high0 = _mm_unpackhi_pd(scale_sum0, scale_sum0);
        double scale0 = _mm_cvtsd_f64(_mm_add_pd(scale_sum0, get_high0));
        // Scale timestep 0

        C[0] = 1./scale0;
        double c0 = 1./scale0;
        scales[0] = scale0;
        double scales0 = scale0;

        __m256d c_0 = _mm256_set1_pd(C[0]);
        for(int i = 0; i < N-3; i+=4) {
            __m256d fw_vec = _mm256_loadu_pd(FW + 0*N + i); 
            __m256d result = _mm256_mul_pd(fw_vec, c_0);
            _mm256_storeu_pd(FW + 0*N + i, result);
        }

        for(int t = 1; t < T; t++) {
            __m256d scale_vec = _mm256_setzero_pd();
            int obs = O[t];

            for(int i = 0; i < N-3; i+=4) {
                __m256d accum0 = _mm256_setzero_pd();

                for(int j = 0; j < N; j++) {

                    __m256d fw_vec0 = _mm256_set1_pd(FW[(t-1)*N + j]);
                    __m256d a_vec0 = _mm256_loadu_pd(A + j*N + i);
                    accum0 = _mm256_fmadd_pd(fw_vec0, a_vec0, accum0);

                }
                __m256d b_vec = _mm256_loadu_pd(B + obs*N + i);
                __m256d result = _mm256_mul_pd(accum0, b_vec);
                scale_vec = _mm256_add_pd(result, scale_vec);
                _mm256_storeu_pd(FW + t*N + i, result);
            }

            __m128d scale_low = _mm256_castpd256_pd128(scale_vec);
            __m128d scale_high = _mm256_extractf128_pd(scale_vec, 1);
            __m128d scale_sum = _mm_add_pd(scale_low, scale_high);
            __m128d get_high = _mm_unpackhi_pd(scale_sum, scale_sum);
            double scale = _mm_cvtsd_f64(_mm_add_pd(scale_sum, get_high));

            C[t] = 1./scale;
            scales[t] = scale;

            for(int i = 0; i < N; i++) {
                FW[t*N + i]*= C[t];
            }
        }

        REGION_END(forward_vars)

        REGION_BEGIN(backward_vars)
        // Calculate the Backward trellis (scaled)
        __m256d bw_c_t_1 = _mm256_set1_pd(C[T-1]);
        for(int i = 0; i < N-3; i+= 4) {
            _mm256_storeu_pd(BW + (T-1)*N + i, bw_c_t_1);
        }

        for(int t = T-2; t >= 0; t--) {
            int obs = O[t+1];
            int limit = N-3;
            for(int i = 0; i < limit; i+=4) {
                __m256d accum0 = _mm256_setzero_pd();
                __m256d accum1 = _mm256_setzero_pd();
                __m256d accum2 = _mm256_setzero_pd();
                __m256d accum3 = _mm256_setzero_pd();

                for(int j = 0; j < limit; j+=4) {
                    __m256d A_i0 = _mm256_loadu_pd(A + i*N + j);
                    __m256d A_i1 = _mm256_loadu_pd(A + (i+1)*N + j);
                    __m256d A_i2 = _mm256_loadu_pd(A + (i+2)*N + j);
                    __m256d A_i3 = _mm256_loadu_pd(A + (i+3)*N + j);

                    __m256d BW_t = _mm256_loadu_pd(BW + (t+1)*N + j);
                    __m256d B_t = _mm256_loadu_pd(B + obs*N + j); 

                    __m256d mul = _mm256_mul_pd(BW_t, B_t);

                    accum0 = _mm256_fmadd_pd(A_i0, mul, accum0);
                    accum1 = _mm256_fmadd_pd(A_i1, mul, accum1);
                    accum2 = _mm256_fmadd_pd(A_i2, mul, accum2);
                    accum3 = _mm256_fmadd_pd(A_i3, mul, accum3);
                }
                __m256d hadd0 = _mm256_hadd_pd(accum0, accum1);
                __m256d hadd1 = _mm256_hadd_pd(accum2, accum3);
                
                __m256d blend = _mm256_blend_pd(hadd0, hadd1, 0xC);
                __m256d permute = _mm256_permute2f128_pd(hadd0, hadd1, 0x21); 
                __m256d sum = _mm256_add_pd(blend, permute);
                __m256d c_t = _mm256_set1_pd(C[t]);
                __m256d result = _mm256_mul_pd(c_t,sum);
                _mm256_storeu_pd(BW + t*N + i, result);
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

            double pi = FW[i] * BW[i] * scales0;

            sum_os[i*M + o0] = pi;
            PI[i] = pi;

            double denom0 = 0;
            double denom1 = 0;
            double denom2 = 0;
            double denom3 = 0;

            for(int t = 1; t < T-3; t+=4) {
                double toadd0 = FW[i + (t  )*N] * BW[i + (t  )*N] * scales[t];
                double toadd1 = FW[i + (t+1)*N] * BW[i + (t+1)*N] * scales[t+1];
                double toadd2 = FW[i + (t+2)*N] * BW[i + (t+2)*N] * scales[t+2];
                double toadd3 = FW[i + (t+3)*N] * BW[i + (t+3)*N] * scales[t+3];

                denom0 += toadd0;
                denom1 += toadd1;
                denom2 += toadd2;
                denom3 += toadd3;

                sum_os[i*M + O[t  ]] += toadd0;
                sum_os[i*M + O[t+1]] += toadd1;
                sum_os[i*M + O[t+2]] += toadd2;
                sum_os[i*M + O[t+3]] += toadd3;
            }

            // leftover 2 iterations
            double toadd0 = FW[i + (T-3)*N] * BW[i + (T-3)*N] * scales[T-3];
            double toadd1 = FW[i + (T-2)*N] * BW[i + (T-2)*N] * scales[T-2];

            denom0 += toadd0;
            denom1 += toadd1;

            sum_os[i*M + ot3] += toadd0;
            sum_os[i*M + ot2] += toadd1;

            double denom = denom0+denom1+denom2+denom3+pi;

            double lastadd = FW[i + (T-1)*N] * BW[i + (T-1)*N] * scalest1;
            denoms[i] = denom + lastadd;
            sum_os[i*M + ot1] += lastadd;

            for(int j = 0; j < N; j++) {
                double num0 = 0.;
                double num1 = 0.;
                double num2 = 0.;
                double num3 = 0.;

                for(int t = 0; t < T-4; t+=4) {
                    num0 += FW[i + (t  )*N] * B[j + O[t+1]*N] * BW[j + (t+1)*N];
                    num1 += FW[i + (t+1)*N] * B[j + O[t+2]*N] * BW[j + (t+2)*N];
                    num2 += FW[i + (t+2)*N] * B[j + O[t+3]*N] * BW[j + (t+3)*N];
                    num3 += FW[i + (t+3)*N] * B[j + O[t+4]*N] * BW[j + (t+4)*N];
                }

                num0 += FW[i + (T-4)*N] * B[j + ot3*N] * BW[j + (T-3)*N];
                num1 += FW[i + (T-3)*N] * B[j + ot2*N] * BW[j + (T-2)*N];
                num2 += FW[i + (T-2)*N] * B[j + ot1*N] * BW[j + (T-1)*N];

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
                B[i + (o  )*N] = sum_os[i*M + o  ]/denomsi;
                B[i + (o+1)*N] = sum_os[i*M + o+1]/denomsi;
                B[i + (o+2)*N] = sum_os[i*M + o+2]/denomsi;
                B[i + (o+3)*N] = sum_os[i*M + o+3]/denomsi;

                sum_os[i*M + o] = 0;
                sum_os[i*M + o+1] = 0;
                sum_os[i*M + o+2] = 0;
                sum_os[i*M + o+3] = 0;
            }
        }
        REGION_END(update_emission)
    }

    REGION_END(baum_welch)
    transpose(M, N, B);
}
