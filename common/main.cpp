#include <bits/stdc++.h>
#include "common.h"

int main() {
    size_t N, M, T, num_iter;
    std::cin >> N >> M >> T >> num_iter;

    size_t cycles = perf(N, M, T, num_iter);

    std::cout << cycles << std::endl;

    return 0;
}
