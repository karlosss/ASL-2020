#include <bits/stdc++.h>
#include "common.h"

void debug_mode() {
    size_t N, M, T, num_iter;
    std::cin >> N >> M >> T >> num_iter;
    if(std::cin.fail() || std::cin.eof()){
        std::cout << "I/O error!";
        exit(2);
    }

    int* O;
    double *PI, *A, *B, *FW, *BW, *C;

    generate_input(N, M, T, &O, &PI, &A, &B, &FW, &BW, &C);

    baum_welch(PI, A, B, O, FW, BW, C, N, M, T, num_iter);
}

void test_mode(){
    size_t N = 20;
    size_t M = 20;
    size_t T = 100;
    size_t num_iter = 100;

    int* O;
    double *PI, *A, *B, *FW, *BW, *C;
    int* O2;
    double *PI2, *A2, *B2, *FW2, *BW2, *C2;

    generate_input(N, M, T, &O, &PI, &A, &B, &FW, &BW, &C);
    clone_input(N, M, T, O, PI, A, B, &O2, &PI2, &A2, &B2, &FW2, &BW2, &C2);

    baum_welch_base(PI, A, B, O, FW, BW, C, N, M, T, num_iter);
    baum_welch(PI2, A2, B2, O2, FW2, BW2, C2, N, M, T, num_iter);

    if(compare_outputs(N, M, PI, A, B, PI2, A2, B2)) exit(0);
    exit(1);
}

void runtime_mode(){
    size_t N, M, T, num_iter;
    std::cin >> N >> M >> T >> num_iter;
    if(std::cin.fail() || std::cin.eof()){
        std::cout << "I/O error!";
        exit(2);
    }

    std::cout << runtime(N, M, T, num_iter) << std::endl;

    exit(0);
}

void flops_mode(){
    size_t N, M, T, num_iter;
    std::cin >> N >> M >> T >> num_iter;
    if(std::cin.fail() || std::cin.eof()){
        std::cout << "I/O error!";
        exit(2);
    }

    std::cout << flop_count(N, M, T, num_iter) << std::endl;

    exit(0);
}
