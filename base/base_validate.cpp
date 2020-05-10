#include <bits/stdc++.h>
#include "./common.h"

int main() {
    int N = 2; int M = 2; int T = 4; int num_iter = 1;
    
    double PI[2] = {0.2, 0.8};
    double A[4] = {1./2., 1./2., 1./4., 3./4.};
    double B[4] = {1./4., 3./4., 4./5., 1./5.};
    int O[4] = {1, 0, 0, 1};
    double FW[8] = {0.};
    double BW[8] = {0.};
    double C[4] = {0.};

    baum_welch_base(PI, A, B, O, FW, BW, C, N, M, T, num_iter);

    double PI2[2] = {112485.0/267581.0, 155096.0/267581.0};
    double A2[4] = {1205.0/3433.0, 2228.0/3433.0, 21025.0/76741.0, 55716.0/76741.0};
    double B2[4] = {7633.0/34696.0, 27063.0/34696.0, 114708.0/180841.0, 66133.0/180841.0};
    if(compare_outputs(N, M, PI, A, B, PI2, A2, B2)) exit(0);
    exit(1);
}