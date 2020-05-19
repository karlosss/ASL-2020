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
    double* scales = (double*) malloc(T*sizeof(double));

    transpose(N, M, B);

    REGION_BEGIN(baum_welch)

    for(int it = 0; it < n_iter; it++) {

        // Calculate the Forward trellis (scaled)
        REGION_BEGIN(forward_vars)
        __m256d scale_vec0 = _mm256_setzero_pd();
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
        scales[0] = scale0;

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
        // update the Initial State probabilities
        for(int i = 0; i < N; i++) {
            PI[i] = FW[0*N + i]*BW[0*N + i]*scales[0];
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
                    n_accum0 += FW[t*N + i]*A[i*N + j]*B[O[t+1]*N + j]*BW[(t+1)*N + j];
                    n_accum1 += FW[(t+1)*N + i]*A[i*N + j]*B[O[t+2]*N + j]*BW[(t+2)*N + j];
                    n_accum2 += FW[(t+2)*N + i]*A[i*N + j]*B[O[t+3]*N + j]*BW[(t+3)*N + j];
                    n_accum3 += FW[(t+3)*N + i]*A[i*N + j]*B[O[t+4]*N + j]*BW[(t+4)*N + j];

                    d_accum0 += FW[t*N + i]*BW[t*N + i]*scales[t];
                    d_accum1 += FW[(t+1)*N + i]*BW[(t+1)*N + i]*scales[t+1];
                    d_accum2 += FW[(t+2)*N + i]*BW[(t+2)*N + i]*scales[t+2];
                    d_accum3 += FW[(t+3)*N + i]*BW[(t+3)*N + i]*scales[t+3];
                }
                for(; t < T-1; t++) {
                    n_accum0 += FW[t*N + i]*A[i*N + j]*B[O[t+1]*N + j]*BW[(t+1)*N + j];
                    d_accum0 += FW[t*N + i]*BW[t*N + i]*scales[t]; 
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
                    d_accum0 += FW[t*N + j]*BW[t*N + j]*scales[t];
                    d_accum1 += FW[(t+1)*N + j]*BW[(t+1)*N + j]*scales[t+1];
                    d_accum2 += FW[(t+2)*N + j]*BW[(t+2)*N + j]*scales[t+2];
                    d_accum3 += FW[(t+3)*N + j]*BW[(t+3)*N + j]*scales[t+3];

                    if(O[t] == o) n_accum0 += FW[t*N + j]*BW[t*N + j]*scales[t];
                    if(O[t+1] == o) n_accum1 += FW[(t+1)*N + j]*BW[(t+1)*N + j]*scales[t+1];
                    if(O[t+2] == o) n_accum2 += FW[(t+2)*N + j]*BW[(t+2)*N + j]*scales[t+2];
                    if(O[t+3] == o) n_accum3 += FW[(t+3)*N + j]*BW[(t+3)*N + j]*scales[t+3];
                }
                B[o*N + j] = (n_accum0+n_accum1+n_accum2+n_accum3)/(d_accum0+d_accum1+d_accum2+d_accum3);
            }
        }

        REGION_END(update_emission)
    }

    REGION_END(baum_welch)
    transpose(M, N, B);
}
