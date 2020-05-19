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
    
    REGION_BEGIN(baum_welch)
    //unroll factor and clean up var
    int unroll_n_outer = 16;
    int unroll_n_inner = 4;

    for(int it = 0; it < n_iter; it++) {

        // Calculate the Forward trellis (scaled)
        REGION_BEGIN(forward_vars)
        int o0 = O[0];
        //constant offset vectors
        __m256d zero_vec = _mm256_setzero_pd();
        //__m256i gather_N = _mm256_set_epi64x(3*N, 2*N, N, 0);
        __m256i gather_T = _mm256_set_epi64x(3*T, 2*T, T, 0);        
        __m256i gather_M = _mm256_set_epi64x(3*M, 2*M, M, 0);
        __m256i B_offset = _mm256_set_epi64x(4*M, 4*M, 4*M, 4*M);
        __m256d scale_vec = _mm256_setzero_pd();

        double scale = 0.;

        __m256d forward_vec, pi_vec, b_vec;
        for(int i = 0; i < N; i+=unroll_n_inner) {

            pi_vec = _mm256_loadu_pd(PI+i);
            b_vec = _mm256_i64gather_pd(B+i*M+O[0], gather_M, 8);
            forward_vec = _mm256_mul_pd(pi_vec, b_vec);
            scale_vec = _mm256_add_pd(scale_vec, forward_vec);

            __m128d first_half = _mm256_extractf128_pd(forward_vec, 0);
            __m128d second_half = _mm256_extractf128_pd(forward_vec, 1);
            _mm_storel_pd(FW + i*T, first_half);
            _mm_storeh_pd(FW + (i+1)*T, first_half);
            _mm_storel_pd(FW + (i+2)*T, second_half);
            _mm_storeh_pd(FW + (i+3)*T, second_half);
        }

        __m256d sum0 = _mm256_unpacklo_pd(scale_vec, zero_vec);
        __m256d sum1 = _mm256_unpackhi_pd(scale_vec, zero_vec);
        sum0 = _mm256_add_pd(sum0, sum1); 
        __m128d top = _mm256_extractf128_pd(sum0, 1);
        sum1 = _mm256_insertf128_pd(zero_vec, top, 0);
        sum0 = _mm256_add_pd(sum0, sum1);
        
        scale = ((double *)&sum0)[0];
        
        // Scale timestep 0
        C[0] = 1./scale;
        scales[0] = scale;
        double c0 = 1./scale;
        double scales0 = scale;

        //scalar replace C[0] (decide wheter to load or calculate again)
        //double C0 = C[0];
        __m256d c_vec = _mm256_set1_pd(C[0]);
        for(int i = 0; i < N; i+=unroll_n_inner) {
            forward_vec = _mm256_i64gather_pd(FW+i*T, gather_T, 8);
            forward_vec = _mm256_mul_pd(forward_vec, c_vec);

            __m128d first_half = _mm256_extractf128_pd(forward_vec, 0);
            __m128d second_half = _mm256_extractf128_pd(forward_vec, 1);
            _mm_storel_pd(FW + i*T, first_half);
            _mm_storeh_pd(FW + (i+1)*T, first_half);
            _mm_storel_pd(FW + (i+2)*T, second_half);
            _mm_storeh_pd(FW + (i+3)*T, second_half);
        }

        for(int t = 1; t < T; t++) {

            //preload O[t]
            int Ot = O[t];

            scale = 0.;
            __m256d scale0 = _mm256_setzero_pd();
            __m256d scale1 = _mm256_setzero_pd();
            for(int i = 0; i < N; i+=unroll_n_outer) {
                
                //4 times due to outer unroll
                //scalar replace FW[i*T+t]
                __m256d a_vec0, a_vec1, a_vec2, a_vec3;

                __m256d acc0 = _mm256_setzero_pd();
                __m256d acc1 = _mm256_setzero_pd();
                __m256d acc2 = _mm256_setzero_pd();
                __m256d acc3 = _mm256_setzero_pd();

                for(int j = 0; j < N; j+=unroll_n_inner) {

                    //strength reduce
                    int jTtm = j*T+t-1;
                    //scalar replace
                    __m256d forward_0 = _mm256_set1_pd(FW[jTtm]);
                    __m256d forward_1 = _mm256_set1_pd(FW[jTtm+T]);
                    __m256d forward_2 = _mm256_set1_pd(FW[jTtm+2*T]);
                    __m256d forward_3 = _mm256_set1_pd(FW[jTtm+3*T]);
                    //strength reduce
                    int jNi0 = j*N+i;
                    int jNi1 = j*N+i+N;
                    int jNi2 = j*N+i+2*N;
                    int jNi3 = j*N+i+3*N;

                    //////////////////////////////////////////////
                    a_vec0 = _mm256_loadu_pd(A+jNi0);
                    acc0 = _mm256_fmadd_pd(forward_0, a_vec0, acc0);
                    a_vec1 = _mm256_loadu_pd(A+jNi0+4);
                    acc1 = _mm256_fmadd_pd(forward_0, a_vec1, acc1);
                    a_vec2 = _mm256_loadu_pd(A+jNi0+8);
                    acc2 = _mm256_fmadd_pd(forward_0, a_vec2, acc2);
                    a_vec3 = _mm256_loadu_pd(A+jNi0+12);
                    acc3 = _mm256_fmadd_pd(forward_0, a_vec3, acc3);
                    //////////////////////////////////////////////

                    //////////////////////////////////////////////
                    a_vec0 = _mm256_loadu_pd(A+jNi1);
                    acc0 = _mm256_fmadd_pd(forward_1, a_vec0, acc0);
                    a_vec1 = _mm256_loadu_pd(A+jNi1+4);
                    acc1 = _mm256_fmadd_pd(forward_1, a_vec1, acc1);
                    a_vec2 = _mm256_loadu_pd(A+jNi1+8);
                    acc2 = _mm256_fmadd_pd(forward_1, a_vec2, acc2);
                    a_vec3 = _mm256_loadu_pd(A+jNi1+12);
                    acc3 = _mm256_fmadd_pd(forward_1, a_vec3, acc3);
                    //////////////////////////////////////////////

                    //////////////////////////////////////////////
                    a_vec0 = _mm256_loadu_pd(A+jNi2);
                    acc0 = _mm256_fmadd_pd(forward_2, a_vec0, acc0);
                    a_vec1 = _mm256_loadu_pd(A+jNi2+4);
                    acc1 = _mm256_fmadd_pd(forward_2, a_vec1, acc1);
                    a_vec2 = _mm256_loadu_pd(A+jNi2+8);
                    acc2 = _mm256_fmadd_pd(forward_2, a_vec2, acc2);
                    a_vec3 = _mm256_loadu_pd(A+jNi2+12);
                    acc3 = _mm256_fmadd_pd(forward_2, a_vec3, acc3);
                    //////////////////////////////////////////////

                    //////////////////////////////////////////////
                    a_vec0 = _mm256_loadu_pd(A+jNi3);
                    acc0 = _mm256_fmadd_pd(forward_3, a_vec0, acc0);
                    a_vec1 = _mm256_loadu_pd(A+jNi3+4);
                    acc1 = _mm256_fmadd_pd(forward_3, a_vec1, acc1);
                    a_vec2 = _mm256_loadu_pd(A+jNi3+8);
                    acc2 = _mm256_fmadd_pd(forward_3, a_vec2, acc2);
                    a_vec3 = _mm256_loadu_pd(A+jNi3+12);
                    acc3 = _mm256_fmadd_pd(forward_3, a_vec3, acc3);
                    //////////////////////////////////////////////
   
                }
                
                //move multiplication post loop
                //strength reduce
                int iMOt = i*M+Ot;

                __m256i b_gather = gather_M;
                b_vec = _mm256_i64gather_pd(B+iMOt, b_gather, 8);
                acc0 = _mm256_mul_pd(acc0, b_vec);
                b_gather = _mm256_add_epi64(b_gather, B_offset);
                b_vec = _mm256_i64gather_pd(B+iMOt, b_gather, 8);
                acc1 = _mm256_mul_pd(acc1, b_vec);
                b_gather = _mm256_add_epi64(b_gather, B_offset);
                b_vec = _mm256_i64gather_pd(B+iMOt, b_gather, 8);
                acc2 = _mm256_mul_pd(acc2, b_vec);
                b_gather = _mm256_add_epi64(b_gather, B_offset);
                b_vec = _mm256_i64gather_pd(B+iMOt, b_gather, 8);
                acc3 = _mm256_mul_pd(acc3, b_vec);


                scale0 = _mm256_add_pd(acc0, acc1);
                scale1 = _mm256_add_pd(acc2, acc3);
                scale0 = _mm256_add_pd(scale0, scale1);

                __m256d sum0 = _mm256_unpacklo_pd(scale0, zero_vec);
                __m256d sum1 = _mm256_unpackhi_pd(scale0, zero_vec);
                sum0 = _mm256_add_pd(sum0, sum1); 
                __m128d top = _mm256_extractf128_pd(sum0, 1);
                sum1 = _mm256_insertf128_pd(zero_vec, top, 0);
                sum0 = _mm256_add_pd(sum0, sum1);
        
                scale += ((double *)&sum0)[0];

                
                //revert scalar replace
                //strenght reduce
                int iTt = i*T+t;

                __m128d first_half = _mm256_extractf128_pd(acc0, 0);
                __m128d second_half = _mm256_extractf128_pd(acc0, 1);
                _mm_storel_pd(FW + iTt, first_half);
                _mm_storeh_pd(FW + iTt+T, first_half);
                _mm_storel_pd(FW + iTt+2*T, second_half);
                _mm_storeh_pd(FW + iTt+3*T, second_half);

                first_half = _mm256_extractf128_pd(acc1, 0);
                second_half = _mm256_extractf128_pd(acc1, 1);
                _mm_storel_pd(FW + iTt+4*T, first_half);
                _mm_storeh_pd(FW + iTt+5*T, first_half);
                _mm_storel_pd(FW + iTt+6*T, second_half);
                _mm_storeh_pd(FW + iTt+7*T, second_half);

                first_half = _mm256_extractf128_pd(acc2, 0);
                second_half = _mm256_extractf128_pd(acc2, 1);
                _mm_storel_pd(FW + iTt+8*T, first_half);
                _mm_storeh_pd(FW + iTt+9*T, first_half);
                _mm_storel_pd(FW + iTt+10*T, second_half);
                _mm_storeh_pd(FW + iTt+11*T, second_half);

                first_half = _mm256_extractf128_pd(acc3, 0);
                second_half = _mm256_extractf128_pd(acc3, 1);
                _mm_storel_pd(FW + iTt+12*T, first_half);
                _mm_storeh_pd(FW + iTt+13*T, first_half);
                _mm_storel_pd(FW + iTt+14*T, second_half);
                _mm_storeh_pd(FW + iTt+15*T, second_half);
      
            }
            
            C[t] = 1./scale;
            scales[t] = scale;
            c_vec = _mm256_set1_pd(C[t]);

            for(int i = 0; i < N; i+=unroll_n_inner) {

                forward_vec = _mm256_i64gather_pd(FW+i*T+t, gather_T, 8);
                forward_vec = _mm256_mul_pd(forward_vec, c_vec);

                __m128d first_half = _mm256_extractf128_pd(forward_vec, 0);
                __m128d second_half = _mm256_extractf128_pd(forward_vec, 1);
                _mm_storel_pd(FW + i*T+t, first_half);
                _mm_storeh_pd(FW + i*T+t+T, first_half);
                _mm_storel_pd(FW + i*T+t+2*T, second_half);
                _mm_storeh_pd(FW + i*T+t+3*T, second_half);
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
            int limit = N-3;
            for(int i = 0; i < limit; i+=4) {
                __m256d accum0 = _mm256_setzero_pd();
                __m256d accum1 = _mm256_setzero_pd();
                __m256d accum2 = _mm256_setzero_pd();
                __m256d accum3 = _mm256_setzero_pd();
                
                for(int j = 0; j < limit; j+=4) {
                    __m256i idx = _mm256_set_epi64x(3*N,2*N, N, 0);

                    __m256d BW0 = _mm256_set1_pd(BW[j*T + t+1]);
                    __m256d B0 = _mm256_set1_pd(B[j*M + obs]);
                    __m256d A0 = _mm256_i64gather_pd(A + i*N + j, idx, 8);
                    __m256d mul0 = _mm256_mul_pd(BW0, B0);
                    accum0 = _mm256_fmadd_pd(A0, mul0, accum0);
                    
                    __m256d BW1 = _mm256_set1_pd(BW[(j+1)*T + t+1]);
                    __m256d B1 = _mm256_set1_pd(B[(j+1)*M + obs]);
                    __m256d A1 = _mm256_i64gather_pd(A + i*N + j+1, idx, 8);
                    __m256d mul1 = _mm256_mul_pd(BW1, B1);
                    accum1 = _mm256_fmadd_pd(A1, mul1, accum1);

                    __m256d BW2 = _mm256_set1_pd(BW[(j+2)*T + t+1]);
                    __m256d B2 = _mm256_set1_pd(B[(j+2)*M + obs]);
                    __m256d A2 = _mm256_i64gather_pd(A + i*N + j+2, idx, 8);
                    __m256d mul2 = _mm256_mul_pd(BW2, B2);
                    accum2 = _mm256_fmadd_pd(A2, mul2, accum2);

                    __m256d BW3 = _mm256_set1_pd(BW[(j+3)*T + t+1]);
                    __m256d B3 = _mm256_set1_pd(B[(j+3)*M + obs]);
                    __m256d A3 = _mm256_i64gather_pd(A + i*N + j+3, idx, 8);
                    __m256d mul3 = _mm256_mul_pd(BW3, B3);
                    accum3 = _mm256_fmadd_pd(A3, mul3, accum3);
                }
                __m256d c_t = _mm256_set1_pd(C[t]);
                __m256d res0 = _mm256_add_pd(accum0, accum1);
                __m256d res1 = _mm256_add_pd(accum2, accum3);
                __m256d res2 = _mm256_add_pd(res0, res1);
                __m256d result = _mm256_mul_pd(c_t, res2);
                
                __m128d first_half = _mm256_extractf128_pd(result, 0);
                __m128d second_half = _mm256_extractf128_pd(result, 1);
                _mm_storel_pd(BW + i*T + t, first_half);
                _mm_storeh_pd(BW + (i+1)*T + t, first_half);
                _mm_storel_pd(BW + (i+2)*T + t, second_half);
                _mm_storeh_pd(BW + (i+3)*T + t, second_half);
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
}