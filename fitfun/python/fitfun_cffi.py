from cffi import FFI
ffi = FFI()
ffi.cdef("""\
    void meth_e2(size_t P, size_t x, size_t y, int n, int m);
""")

ff = ffi.verify("""\
#include <fitfun.h>

data D;

void parse_args(data *D, size_t P, size_t x, size_t y, int n, int m) {
    D->P = (double*) P;
    D->x = (double*) x;
    D->y = (double*) y;
    D->n = n;
    D->m = m;
}

void meth_e2(size_t P, size_t x, size_t y, int n, int m) {
    parse_args(&D, P, x, y, n, m);
    e2(&D);
}
""", libraries=['fitfun'])

def get_ptr_array(x):
    return x.__array_interface__['data'][0]

try:
    from accel import _get_ptr as get_ptr
except ImportError:
    get_ptr = get_ptr_array

try:
    import numpypy as np
except ImportError:
    import numpy as np


def e2(P, x, y):
    ff.meth_e2(get_ptr(P), get_ptr(x), get_ptr(y), P.size, x.size)
    return y

P = np.array([1., 1.])
x = np.arange(4.)
y = np.empty(x.shape, x.dtype)

print e2(P, x, y)

print P[0]*np.exp(-P[1]*x)


