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
        __m256d scale_vec1 = _mm256_setzero_pd();
        for(int i = 0; i < N-7; i+=8) {
            __m256d pi_vec = _mm256_loadu_pd(PI + i);
            __m256d pi_vec1 = _mm256_loadu_pd(PI + i+4);
            __m256d b_vec = _mm256_loadu_pd(B + O[0]*N + i);
            __m256d b_vec1 = _mm256_loadu_pd(B + O[0]*N + i+4);
            __m256d result = _mm256_mul_pd(pi_vec, b_vec);
            __m256d result1 = _mm256_mul_pd(pi_vec1, b_vec1);
            _mm256_storeu_pd(FW + 0*N + i, result);
            _mm256_storeu_pd(FW + 0*N + i + 4, result1);
            scale_vec0 = _mm256_add_pd(scale_vec0, result);
            scale_vec1 = _mm256_add_pd(scale_vec1, result1);
        }
        scale_vec0 = _mm256_add_pd(scale_vec0, scale_vec1);
        scale_vec0 = _mm256_hadd_pd(scale_vec0, _mm256_permute2f128_pd(scale_vec0, scale_vec0, 1));
        scale_vec0 = _mm256_hadd_pd(scale_vec0, scale_vec0);
        double scale0 = _mm_cvtsd_f64(_mm256_castpd256_pd128(scale_vec0));
        // Scale timestep 0

        C[0] = 1./scale0;
        scales[0] = scale0;

        __m256d c_0 = _mm256_broadcast_sd(&C[0]);
        for(int i = 0; i < N-7; i+=8) {
            __m256d fw_vec = _mm256_loadu_pd(FW + 0*N + i);
            __m256d fw_vec1 = _mm256_loadu_pd(FW + 0*N + i+4);
            _mm256_storeu_pd(FW + 0*N + i, _mm256_mul_pd(fw_vec, c_0));
            _mm256_storeu_pd(FW + 0*N + i+4, _mm256_mul_pd(fw_vec1, c_0));
        }

        for(int t = 1; t < T; t++) {
            __m256d scale_vec = _mm256_setzero_pd();
            int obs = O[t];

            for(int i = 0; i < N-3; i+=4) {
                __m256d accum0 = _mm256_setzero_pd();
                __m256d accum1 = _mm256_setzero_pd();
                __m256d accum2 = _mm256_setzero_pd();
                __m256d accum3 = _mm256_setzero_pd();

                for(int j = 0; j < N-3; j+=4) {

                    __m256d fw_vec0 = _mm256_broadcast_sd(&FW[(t-1)*N + j]);
                    __m256d a_vec0 = _mm256_loadu_pd(A + j*N + i);
                    accum0 = _mm256_fmadd_pd(fw_vec0, a_vec0, accum0);

                    __m256d fw_vec1 = _mm256_broadcast_sd(&FW[(t-1)*N + j+1]);
                    __m256d a_vec1 = _mm256_loadu_pd(A + (j+1)*N + i);
                    accum1 = _mm256_fmadd_pd(fw_vec1, a_vec1, accum1);

                    __m256d fw_vec2 = _mm256_broadcast_sd(&FW[(t-1)*N + j+2]);
                    __m256d a_vec2 = _mm256_loadu_pd(A + (j+2)*N + i);
                    accum2 = _mm256_fmadd_pd(fw_vec2, a_vec2, accum2);

                    __m256d fw_vec3 = _mm256_broadcast_sd(&FW[(t-1)*N + j+3]);
                    __m256d a_vec3 = _mm256_loadu_pd(A + (j+3)*N + i);
                    accum3 = _mm256_fmadd_pd(fw_vec3, a_vec3, accum3);

                }
                __m256d sum0 = _mm256_add_pd(accum0, accum1);
                __m256d sum1 = _mm256_add_pd(accum2, accum3);
                __m256d sum2 = _mm256_add_pd(sum0, sum1);
                __m256d b_vec = _mm256_loadu_pd(B + obs*N + i);
                __m256d result = _mm256_mul_pd(sum2, b_vec);
                scale_vec = _mm256_add_pd(result, scale_vec);
                _mm256_storeu_pd(FW + t*N + i, result);
            }

            scale_vec = _mm256_hadd_pd(scale_vec, _mm256_permute2f128_pd(scale_vec, scale_vec, 1));
            scale_vec = _mm256_hadd_pd(scale_vec, scale_vec);
            double scale = _mm_cvtsd_f64(_mm256_castpd256_pd128(scale_vec));

            C[t] = 1./scale;
            scales[t] = scale;
            __m256d c_t = _mm256_broadcast_sd(&C[t]);
            for(int i = 0; i < N-7; i+=8) {
                __m256d fw_vec = _mm256_loadu_pd(FW + t*N + i); 
                __m256d fw_vec1 = _mm256_loadu_pd(FW + t*N + i+4);
                _mm256_storeu_pd(FW + t*N + i, _mm256_mul_pd(fw_vec, c_t));
                _mm256_storeu_pd(FW + t*N + i+4, _mm256_mul_pd(fw_vec1, c_t));
            }
        }

        REGION_END(forward_vars)

        REGION_BEGIN(backward_vars)
        // Calculate the Backward trellis (scaled)
        __m256d bw_c_t_1 = _mm256_broadcast_sd(&C[T-1]);
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
                __m256d c_t = _mm256_broadcast_sd(&C[t]);
                __m256d result = _mm256_mul_pd(c_t,sum);
                _mm256_storeu_pd(BW + t*N + i, result);
            }
        }
        REGION_END(backward_vars)

        REGION_BEGIN(update_initial)
        // update the Initial State probabilities
        __m256d scale_0 = _mm256_broadcast_sd(&scales[0]);
        for(int i = 0; i < N-3; i+=4) {
            __m256d fw_vec0 = _mm256_loadu_pd(FW + 0*N + i);
            __m256d fw_sc_mul = _mm256_mul_pd(fw_vec0,scale_0);
            __m256d bw_vec0 = _mm256_loadu_pd(BW + 0*N + i);
            __m256d pi_res = _mm256_mul_pd(fw_sc_mul, bw_vec0);
            _mm256_storeu_pd(PI + i, pi_res);
        }
        REGION_END(update_initial)

        REGION_BEGIN(update_transition)
        // update the State Transition probabilities
        for(int i = 0; i < N; i++) {
            
            for(int j = 0; j < N-7; j+=8) {
            
                __m256d n_accum0 = _mm256_setzero_pd();
                __m256d n_accum1 = _mm256_setzero_pd();

                __m256d d_accum0 = _mm256_setzero_pd();

                __m256d a_vec0 = _mm256_loadu_pd(A + i*N + j);
                __m256d a_vec1 = _mm256_loadu_pd(A + i*N + j+4);

                for(int t = 0; t < T-1; t++) {
                    //n_accum0 += FW[t*N + i]*A[i*N + j]*B[O[t+1]*N + j]*BW[(t+1)*N + j];
                    __m256d fw_vec0 = _mm256_broadcast_sd(&FW[t*N + i]);

                    __m256d b_vec1 = _mm256_loadu_pd(B + O[t+1]*N + j);
                    __m256d bwj_vec1 = _mm256_loadu_pd(BW + (t+1)*N + j);
                    __m256d fw_a_mul0 = _mm256_mul_pd(fw_vec0, a_vec0);
                    __m256d b_bw_mul0 = _mm256_mul_pd(b_vec1, bwj_vec1);
                    n_accum0 = _mm256_fmadd_pd(fw_a_mul0, b_bw_mul0, n_accum0);

                    __m256d b_vec2 = _mm256_loadu_pd(B + O[t+1]*N + j+4);
                    __m256d bwj_vec2 = _mm256_loadu_pd(BW + (t+1)*N + j+4);
                    __m256d fw_a_mul1 = _mm256_mul_pd(fw_vec0, a_vec1);
                    __m256d b_bw_mul1 = _mm256_mul_pd(b_vec2, bwj_vec2);
                    n_accum1 = _mm256_fmadd_pd(fw_a_mul1, b_bw_mul1, n_accum1);

                    //d_accum0 += FW[t*N + i]*BW[t*N + i]*scales[t];
                    __m256d bwi_vec0 = _mm256_broadcast_sd(&BW[t*N + i]);
                    __m256d scale_t0 = _mm256_broadcast_sd(&scales[t]);
                    __m256d fw_bwi_vec0 = _mm256_mul_pd(fw_vec0, bwi_vec0);
                    d_accum0 = _mm256_fmadd_pd(fw_bwi_vec0, scale_t0, d_accum0);
                }

                __m256d result0 = _mm256_div_pd(n_accum0, d_accum0);
                __m256d result1 = _mm256_div_pd(n_accum1, d_accum0);
                _mm256_storeu_pd(A + i*N + j, result0);
                _mm256_storeu_pd(A + i*N + j+4, result1);
            }
        }
        REGION_END(update_transition)

        REGION_BEGIN(update_emission)
        // update the State Emission probabilities
        for(int o = 0; o < M; o++) {
            for(int j = 0; j < N-7; j+=8) {
                __m256d n_accum0 = _mm256_setzero_pd();
                __m256d n_accum1 = _mm256_setzero_pd();

                __m256d d_accum0 = _mm256_setzero_pd();
                __m256d d_accum1 = _mm256_setzero_pd();

                for(int t = 0; t < T; t++) {
                    //d_accum0 += FW[t*N + j]*BW[t*N + j]*scales[t];
                    __m256d fw_vec0 = _mm256_loadu_pd(FW + t*N + j);
                    __m256d fw_vec1 = _mm256_loadu_pd(FW + t*N + j + 4);
                    __m256d bw_vec0 = _mm256_loadu_pd(BW + t*N + j);
                    __m256d bw_vec1 = _mm256_loadu_pd(BW + t*N + j+4);
                    __m256d scale_t0 = _mm256_broadcast_sd(&scales[t]);
                    __m256d fw_bw_mul0 = _mm256_mul_pd(fw_vec0, bw_vec0);
                    __m256d fw_bw_mul1 = _mm256_mul_pd(fw_vec1, bw_vec1);
                    d_accum0 = _mm256_fmadd_pd(fw_bw_mul0, scale_t0, d_accum0);
                    d_accum1 = _mm256_fmadd_pd(fw_bw_mul1, scale_t0, d_accum1);
                    if(O[t] == o) {
                        //n_accum0 += FW[t*N + j]*BW[t*N + j]*scales[t];
                        n_accum0 = _mm256_fmadd_pd(fw_bw_mul0, scale_t0, n_accum0);
                        n_accum1 = _mm256_fmadd_pd(fw_bw_mul1, scale_t0, n_accum1);
                    }
                }
                __m256d result0 = _mm256_div_pd(n_accum0, d_accum0);
                __m256d result1 = _mm256_div_pd(n_accum1, d_accum1);
                _mm256_storeu_pd(B + o*N + j, result0);
                _mm256_storeu_pd(B + o*N + j+4, result1);
            }
        }

        REGION_END(update_emission)
    }

    REGION_END(baum_welch)
    transpose(M, N, B);
}
