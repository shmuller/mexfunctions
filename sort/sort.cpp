#include "sort.h"

#include <algorithm>
#include "timsort.hpp"

void stdsort(int *a, int n) {
    std::sort(a, a + n);
}

void timsort(int *a, int n) {
    gfx::timsort(a, a + n);
}

