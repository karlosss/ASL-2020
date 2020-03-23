#ifndef ASL_PROJECT_COMMON_H
#define ASL_PROJECT_COMMON_H

#define INLINE static inline __attribute__((always_inline))

typedef void (*fn)();

static int double_cmp (const void * a, const void * b) {
    if (*(double*)a > *(double*)b) return 1;
    else if (*(double*)a < *(double*)b) return -1;
    else return 0;
}

#endif //ASL_PROJECT_COMMON_H
