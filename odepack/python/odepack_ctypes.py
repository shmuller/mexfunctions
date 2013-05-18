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

if __name__ == "__main__":
    try:
        import numpy as np
    except ImportError:
        import numpypy as np

    neq = 2
    n = 7
    ng = 1

    def f(neq, t, y, ydot):
        ydot[0] = -y[0]
        ydot[1] = -y[1]

    def g(neq, t, y, ng, gout):
        gout[0] = y[0] - 0.2

    y = np.ones((n, neq))

    t = np.arange(n, dtype=np.float64)
    
    D.f = odefun_t(f)
    D.neq = neq
    D.n = n
    D.t = get_ptr(t)
    D.y = get_ptr(y)
    
    #D.ng = 1
    #D.g = oderoot_t(g)

    #"""
    ibbox = np.array([0], 'i')
    bbox = np.array([0.2])

    D.g = oderoot_t()
    D.ng = 1
    D.ibbox = get_ptr(ibbox)
    D.bbox = get_ptr(bbox)
    #"""
    codepack.odesolve(byref(D))

    points_done = D.points_done

    t = t[:points_done]
    y = y[:points_done]

    print y[:, 0]
    print y[:, 1]
    print np.exp(-t)


