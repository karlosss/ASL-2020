#ifndef ASL_PROJECT_COMMON_H
#define ASL_PROJECT_COMMON_H

#include <bits/stdc++.h>
#include "inline.h"

#define EPS (1e-6)

void generate_observation(int **PO, size_t len, size_t M);

void generate_m(double **PM, size_t row, size_t col);

INLINE void baum_welch(double*, double*, double*, int*, double*, double*, double*, int, int, int, int) {}

size_t perf(size_t N, size_t M, size_t T, size_t num_iter);

#endif //ASL_PROJECT_COMMON_H
