from ctypes import *
from ctypes.util import find_library

odefun_t = CFUNCTYPE(None, POINTER(c_int), POINTER(c_double), 
                     POINTER(c_double), POINTER(c_double))

oderoot_t = CFUNCTYPE(None, POINTER(c_int), POINTER(c_double), 
                      POINTER(c_double), POINTER(c_int), POINTER(c_double))

class struct_data(Structure):
    __slots__ = ['f', 'neq', 'n', 't', 'y', 'points_done', 'ng', 'g', 'ibbox', 'bbox']

    _fields_ = [
        ('f', odefun_t),
        ('neq', c_int),
        ('n', c_int),
        ('t', c_void_p),
        ('y', c_void_p),
        ('points_done', c_int),
        ('ng', c_int),
        ('g', oderoot_t),
        ('ibbox', c_void_p),
        ('bbox', c_void_p)]

D = struct_data()

def get_ptr_ctypes(x):
    return x.ctypes._data

def get_ptr_array(x):
    return x.__array_interface__['data'][0]

try:
    from accel import _get_ptr as get_ptr
except ImportError:
    get_ptr = get_ptr_array


libname = find_library('odepack')
codepack = CDLL(libname)
codepack.odesolve.argtypes = [POINTER(struct_data)]

def odesolve(f, y, t, ng=0, g=None, ibbox=None, bbox=None):
    D.f = odefun_t(f)
    D.n, D.neq = y.shape
    D.y, D.t = get_ptr(y), get_ptr(t)

    if ng != 0:
        D.ng = ng
        if g is None:
            D.g = oderoot_t()
            D.ibbox, D.bbox = get_ptr(ibbox), get_ptr(bbox)
        else:
            D.g = oderoot_t(g)

    codepack.odesolve(byref(D))
    return D.points_done


if __name__ == "__main__":
    try:
        import numpy as np
    except ImportError:
        import numpypy as np

    neq = 2
    n = 6

    def f(neq, t, y, ydot):
        ydot[0] = -y[0]
        ydot[1] = -y[1]

    y = np.ones((n, neq))
    t = np.arange(n, dtype=np.float64)

    points_done = odesolve(f, y, t)
    t1, y1 = t[:points_done], y[:points_done]

    print y1[:,0]
    print np.exp(-t1)


    def g(neq, t, y, ng, gout):
        gout[0] = y[0] - 0.2

    points_done = odesolve(f, y, t, 1, g)
    t2, y2 = t[:points_done], y[:points_done]

    print y2[:,0]
    print np.exp(-t2)


    ibbox = np.array([0], 'i')
    bbox = np.array([0.2])

    points_done = odesolve(f, y, t, 1, None, ibbox, bbox)
    t3, y3 = t[:points_done], y[:points_done]

    print y3[:,0]
    print np.exp(-t3)

