from odepack_test import test_odesolve

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

def _cast_odefun_t(f):
    def wrapper(neq, t, y, ydot):
        f(y, t, ydot)
    return odefun_t(wrapper)

def _cast_oderoot_t(g):
    def wrapper(neq, t, y, ng, gout):
        g(y, t, gout)
    return oderoot_t(wrapper)


def odesolve(f, y, t, ng=0, g=None, ibbox=None, bbox=None):
    D.f = _cast_odefun_t(f)
    D.n, D.neq = y.shape
    D.y, D.t = get_ptr(y), get_ptr(t)

    if ng != 0:
        D.ng = ng
        if g is None:
            D.g = oderoot_t()
            D.ibbox, D.bbox = get_ptr(ibbox), get_ptr(bbox)
        else:
            D.g = _cast_oderoot_t(g)

    codepack.odesolve(byref(D))
    return D.points_done


if __name__ == "__main__":
    test_odesolve(odesolve)

