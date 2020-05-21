#include <bits/stdc++.h>
#include "common.h"
#include <immintrin.h>

bool double_eq(double a, double b) {
    return fabs(a-b) < 0.001;
}

void compvec(__m256d vec, double a, double b, double c, double d){
    if(!double_eq(((double*) &vec)[0], a) || !double_eq(((double*) &vec)[1], b) || !double_eq(((double*) &vec)[2], c) || !double_eq(((double*) &vec)[3], d)){
        std::cout << "Expected: " << a << " " << b << " " << c << " " << d << ", got " << ((double*) &vec)[0] << " " <<
                  ((double*) &vec)[1] << " " << ((double*) &vec)[2] << " " << ((double*) &vec)[3] << "\n";
    }
}

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

        __m256d vec_scales0 = _mm256_set1_pd(scales0);

        for(int i = 0; i < N; i+=8) {
            __m256d vec_fw0 = _mm256_loadu_pd(FW+i);
            __m256d vec_fw1 = _mm256_loadu_pd(FW+i+4);
            __m256d vec_bw0 = _mm256_loadu_pd(BW+i);
            __m256d vec_bw1 = _mm256_loadu_pd(BW+i+4);

            __m256d vec_pi0 = _mm256_mul_pd(vec_fw0, vec_bw0);
            __m256d vec_pi1 = _mm256_mul_pd(vec_fw1, vec_bw1);

            vec_pi0 = _mm256_mul_pd(vec_pi0, vec_scales0);
            vec_pi1 = _mm256_mul_pd(vec_pi1, vec_scales0);

            _mm256_storeu_pd(sum_os + o0*N+i, vec_pi0);
            _mm256_storeu_pd(sum_os + o0*N+i+4, vec_pi1);

            _mm256_storeu_pd(PI + i, vec_pi0);
            _mm256_storeu_pd(PI + i+4, vec_pi1);

            __m256d vec_denomi0 = _mm256_setzero_pd();
            __m256d vec_denomi1 = _mm256_setzero_pd();

            for(int t = 1; t < T-1; t+=1) {
                __m256d vec_scalest = _mm256_broadcast_sd(scales+t);

                __m256d vec_fwt0 = _mm256_loadu_pd(FW+t*N+i);
                __m256d vec_fwt1 = _mm256_loadu_pd(FW+t*N+i+4);
                __m256d vec_bwt0 = _mm256_loadu_pd(BW+t*N+i);
                __m256d vec_bwt1 = _mm256_loadu_pd(BW+t*N+i+4);

                __m256d vec_toaddi0 = _mm256_mul_pd(vec_bwt0, vec_fwt0);
                __m256d vec_toaddi1 = _mm256_mul_pd(vec_bwt1, vec_fwt1);

                vec_denomi0 = _mm256_fmadd_pd(vec_toaddi0, vec_scalest, vec_denomi0);
                vec_denomi1 = _mm256_fmadd_pd(vec_toaddi1, vec_scalest, vec_denomi1);

                __m256d vec_sumos0 = _mm256_loadu_pd(sum_os + O[t]*N + i);
                __m256d vec_sumos1 = _mm256_loadu_pd(sum_os + O[t]*N + i + 4);

                _mm256_storeu_pd(sum_os + O[t]*N + i, _mm256_fmadd_pd(vec_toaddi0, vec_scalest, vec_sumos0));
                _mm256_storeu_pd(sum_os + O[t]*N + i+4, _mm256_fmadd_pd(vec_toaddi1, vec_scalest, vec_sumos1));
            }

            vec_denomi0 = _mm256_add_pd(vec_denomi0, vec_pi0);
            vec_denomi1 = _mm256_add_pd(vec_denomi1, vec_pi1);

            __m256d vec_scalest1 = _mm256_set1_pd(scalest1);

            __m256d vec_fwl0 = _mm256_loadu_pd(FW+(T-1)*N+i);
            __m256d vec_fwl1 = _mm256_loadu_pd(FW+(T-1)*N+i+4);

            __m256d vec_bwl0 = _mm256_loadu_pd(BW+(T-1)*N+i);
            __m256d vec_bwl1 = _mm256_loadu_pd(BW+(T-1)*N+i+4);

            __m256d vec_lastadd0 = _mm256_mul_pd(vec_fwl0, vec_bwl0);
            __m256d vec_lastadd1 = _mm256_mul_pd(vec_fwl1, vec_bwl1);

            _mm256_storeu_pd(denoms+i, _mm256_fmadd_pd(vec_lastadd0, vec_scalest1, vec_denomi0));
            _mm256_storeu_pd(denoms+i+4, _mm256_fmadd_pd(vec_lastadd1, vec_scalest1, vec_denomi1));

            __m256d vec_sumos0 = _mm256_loadu_pd(sum_os + ot1*N + i);
            __m256d vec_sumos1 = _mm256_loadu_pd(sum_os + ot1*N + i + 4);

            _mm256_storeu_pd(sum_os + ot1*N + i, _mm256_fmadd_pd(vec_lastadd0, vec_scalest1, vec_sumos0));
            _mm256_storeu_pd(sum_os + ot1*N + i+4, _mm256_fmadd_pd(vec_lastadd1, vec_scalest1, vec_sumos1));

            for(int j = 0; j < N; j += 4) {
                __m256d vec_numi00 = _mm256_setzero_pd();
                __m256d vec_numi01 = _mm256_setzero_pd();
                __m256d vec_numi10 = _mm256_setzero_pd();
                __m256d vec_numi11 = _mm256_setzero_pd();
                __m256d vec_numi20 = _mm256_setzero_pd();
                __m256d vec_numi21 = _mm256_setzero_pd();
                __m256d vec_numi30 = _mm256_setzero_pd();
                __m256d vec_numi31 = _mm256_setzero_pd();

                for(int t = 0; t < T-1; t+=1) {
                    __m256d vec_bwt0 = _mm256_broadcast_sd(BW+(t+1)*N+j);
                    __m256d vec_bwt1 = _mm256_broadcast_sd(BW+(t+1)*N+j+1);
                    __m256d vec_bwt2 = _mm256_broadcast_sd(BW+(t+1)*N+j+2);
                    __m256d vec_bwt3 = _mm256_broadcast_sd(BW+(t+1)*N+j+3);

                    __m256d vec_bt0 = _mm256_broadcast_sd(B+O[t+1]*N+j);
                    __m256d vec_bt1 = _mm256_broadcast_sd(B+O[t+1]*N+j+1);
                    __m256d vec_bt2 = _mm256_broadcast_sd(B+O[t+1]*N+j+2);
                    __m256d vec_bt3 = _mm256_broadcast_sd(B+O[t+1]*N+j+3);

                    __m256d vec_fwt0 = _mm256_loadu_pd(FW+t*N+i);
                    __m256d vec_fwt1 = _mm256_loadu_pd(FW+t*N+i+4);

                    __m256d vec_mul_0 = _mm256_mul_pd(vec_bt0, vec_bwt0);
                    __m256d vec_mul_1 = _mm256_mul_pd(vec_bt1, vec_bwt1);
                    __m256d vec_mul_2 = _mm256_mul_pd(vec_bt2, vec_bwt2);
                    __m256d vec_mul_3 = _mm256_mul_pd(vec_bt3, vec_bwt3);

                    vec_numi00 = _mm256_fmadd_pd(vec_fwt0, vec_mul_0, vec_numi00);
                    vec_numi01 = _mm256_fmadd_pd(vec_fwt1, vec_mul_0, vec_numi01);

                    vec_numi10 = _mm256_fmadd_pd(vec_fwt0, vec_mul_1, vec_numi10);
                    vec_numi11 = _mm256_fmadd_pd(vec_fwt1, vec_mul_1, vec_numi11);

                    vec_numi20 = _mm256_fmadd_pd(vec_fwt0, vec_mul_2, vec_numi20);
                    vec_numi21 = _mm256_fmadd_pd(vec_fwt1, vec_mul_2, vec_numi21);

                    vec_numi30 = _mm256_fmadd_pd(vec_fwt0, vec_mul_3, vec_numi30);
                    vec_numi31 = _mm256_fmadd_pd(vec_fwt1, vec_mul_3, vec_numi31);
                }

                __m256d a0 = _mm256_loadu_pd(A+i*N+j);
                __m256d a1 = _mm256_loadu_pd(A+(i+1)*N+j);
                __m256d a2 = _mm256_loadu_pd(A+(i+2)*N+j);
                __m256d a3 = _mm256_loadu_pd(A+(i+3)*N+j);
                __m256d a4 = _mm256_loadu_pd(A+(i+4)*N+j);
                __m256d a5 = _mm256_loadu_pd(A+(i+5)*N+j);
                __m256d a6 = _mm256_loadu_pd(A+(i+6)*N+j);
                __m256d a7 = _mm256_loadu_pd(A+(i+7)*N+j);

                __m256d numi0 = _mm256_insertf128_pd(_mm256_shuffle_pd(vec_numi00,vec_numi10, 0b0000),_mm256_castpd256_pd128(_mm256_shuffle_pd(vec_numi20,vec_numi30, 0b0000)), 1);
                __m256d numi1 = _mm256_insertf128_pd(_mm256_shuffle_pd(vec_numi00,vec_numi10, 0b0011),_mm256_castpd256_pd128(_mm256_shuffle_pd(vec_numi20,vec_numi30, 0b0011)), 1);
                __m256d numi2 = _mm256_unpacklo_pd(_mm256_permute2f128_pd(vec_numi00, vec_numi20, 0x31), _mm256_permute2f128_pd(vec_numi10, vec_numi30, 0x31));
                __m256d numi3 = _mm256_unpackhi_pd(_mm256_permute2f128_pd(vec_numi00, vec_numi20, 0x31), _mm256_permute2f128_pd(vec_numi10, vec_numi30, 0x31));
                __m256d numi4 = _mm256_insertf128_pd(_mm256_shuffle_pd(vec_numi01,vec_numi11, 0b0000),_mm256_castpd256_pd128(_mm256_shuffle_pd(vec_numi21,vec_numi31, 0b0000)), 1);
                __m256d numi5 = _mm256_insertf128_pd(_mm256_shuffle_pd(vec_numi01,vec_numi11, 0b0011),_mm256_castpd256_pd128(_mm256_shuffle_pd(vec_numi21,vec_numi31, 0b0011)), 1);
                __m256d numi6 = _mm256_unpacklo_pd(_mm256_permute2f128_pd(vec_numi01, vec_numi21, 0x31), _mm256_permute2f128_pd(vec_numi11, vec_numi31, 0x31));
                __m256d numi7 = _mm256_unpackhi_pd(_mm256_permute2f128_pd(vec_numi01, vec_numi21, 0x31), _mm256_permute2f128_pd(vec_numi11, vec_numi31, 0x31));

                __m256d t1 = _mm256_permute2f128_pd(vec_denomi0, vec_denomi0, 0x0);
                __m256d t2 = _mm256_permute2f128_pd(vec_denomi0, vec_denomi0, 0x11);
                __m256d bc_denomi0 = _mm256_permute_pd(t1,0);
                __m256d bc_denomi1 = _mm256_permute_pd(t1,0xf);
                __m256d bc_denomi2 = _mm256_permute_pd(t2,0);
                __m256d bc_denomi3 = _mm256_permute_pd(t2,0xf);

                t1 = _mm256_permute2f128_pd(vec_denomi1, vec_denomi1, 0x0);
                t2 = _mm256_permute2f128_pd(vec_denomi1, vec_denomi1, 0x11);
                __m256d bc_denomi4 = _mm256_permute_pd(t1,0);
                __m256d bc_denomi5 = _mm256_permute_pd(t1,0xf);
                __m256d bc_denomi6 = _mm256_permute_pd(t2,0);
                __m256d bc_denomi7 = _mm256_permute_pd(t2,0xf);

                _mm256_storeu_pd(A+i*N+j, _mm256_div_pd(_mm256_mul_pd(numi0, a0), bc_denomi0));
                _mm256_storeu_pd(A+(i+1)*N+j, _mm256_div_pd(_mm256_mul_pd(numi1, a1), bc_denomi1));
                _mm256_storeu_pd(A+(i+2)*N+j, _mm256_div_pd(_mm256_mul_pd(numi2, a2), bc_denomi2));
                _mm256_storeu_pd(A+(i+3)*N+j, _mm256_div_pd(_mm256_mul_pd(numi3, a3), bc_denomi3));
                _mm256_storeu_pd(A+(i+4)*N+j, _mm256_div_pd(_mm256_mul_pd(numi4, a4), bc_denomi4));
                _mm256_storeu_pd(A+(i+5)*N+j, _mm256_div_pd(_mm256_mul_pd(numi5, a5), bc_denomi5));
                _mm256_storeu_pd(A+(i+6)*N+j, _mm256_div_pd(_mm256_mul_pd(numi6, a6), bc_denomi6));
                _mm256_storeu_pd(A+(i+7)*N+j, _mm256_div_pd(_mm256_mul_pd(numi7, a7), bc_denomi7));
            }
        }
        REGION_END(update_transition)

        REGION_BEGIN(update_emission)
        for(int o = 0; o < M; o += 1){
            for(int i = 0; i < N; i+=8) {
                __m256d denoms_i0 = _mm256_loadu_pd(denoms+i);
                __m256d denoms_i1 = _mm256_loadu_pd(denoms+i+4);

                __m256d sum_os_io0 = _mm256_loadu_pd(sum_os+o*N+i);
                __m256d sum_os_io1 = _mm256_loadu_pd(sum_os+o*N+i+4);

                __m256d res0 = _mm256_div_pd(sum_os_io0, denoms_i0);
                __m256d res1 = _mm256_div_pd(sum_os_io1, denoms_i1);

                _mm256_storeu_pd(B + o*N + i, res0);
                _mm256_storeu_pd(B + o*N + i + 4, res1);

                _mm256_storeu_pd(sum_os + i + o*N, _mm256_setzero_pd());
                _mm256_storeu_pd(sum_os + i + o*N + 4, _mm256_setzero_pd());
            }
        }
        REGION_END(update_emission)

    }

    REGION_END(baum_welch)
    transpose(M, N, B);
}
