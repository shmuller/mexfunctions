#include <algorithm>

extern "C" {

void stdsort(int *a, int n) {
    std::sort(a, a + n);
}

}
