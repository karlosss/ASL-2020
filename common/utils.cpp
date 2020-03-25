#include <bits/stdc++.h>

void print_m(double* M, size_t rows, size_t cols){
    for(size_t i = 0; i < rows; i++) {
        for(size_t j = 0; j < cols; j++) {
            std::cout << M[i*rows + j] << " ";
        }
        std::cout << '\n';
    }
}
