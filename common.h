#ifndef ASL_PROJECT_COMMON_H
#define ASL_PROJECT_COMMON_H

#define INLINE static inline __attribute__((always_inline))

typedef void (*fn)(double*, double*, double*, int*, double*, double*, double*, int, int,  int, int);

static int double_cmp (const void * a, const void * b) {
    if (*(double*)a > *(double*)b) return 1;
    else if (*(double*)a < *(double*)b) return -1;
    else return 0;
}

void generate_observation(int **PO, size_t len, size_t M);
void generate_m(double ** PM, size_t row, size_t col);

#endif //ASL_PROJECT_COMMON_H
