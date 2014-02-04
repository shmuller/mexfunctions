#include <algorithm>

#include "timsort.hpp"

extern "C" {

void stdsort(int *a, int n) {
    std::sort(a, a + n);
}

void timsort(int *a, int n) {
    gfx::timsort(a, a + n);
}

}
