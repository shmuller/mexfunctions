from ctypes import c_int, c_double, c_void_p, Structure, byref, CDLL
from ctypes.util import find_library

p_t = p_int_t = c_void_p

class struct_data(Structure):
    __slots__ = ['n', 'm', 'P', 'do_var', 'x', 'y', 'ydata', 'w', 'a']

    _fields_ = [
        ('n', c_int),
        ('m', c_int),
        ('P', p_t),
        ('do_var', p_int_t),
        ('x', p_t),
        ('y', p_t),
        ('ydata', p_t),
        ('w', p_t),
        ('a', p_t)]

D = struct_data()

def get_ptr_ctypes(x):
    return x.ctypes._data

def get_ptr_array(x):
    return x.__array_interface__['data'][0]

try:
    from accel import _get_ptr as get_ptr
except ImportError:
    get_ptr = get_ptr_array


def parse_args(P, x, y, a=None):
    D.P = get_ptr(P)
    D.x = get_ptr(x)
    D.y = get_ptr(y)
    D.n = P.size
    D.m = x.size
    if a is not None:
        D.a = get_ptr(a)

def parse_args_ydata(P, x, y, ydata, a=None):
    parse_args(P, x, y, a)
    D.ydata = get_ptr(ydata)


libname = find_library('fitfun')

cfitfun = CDLL(libname)

def _fun_factory(name):
    getattr(cfitfun, name + '_rms').restype = c_double

    def fun(P, x, y, a=None):
        parse_args(P, x, y, a)
        getattr(cfitfun, name)(byref(D))
        return y

    def fun_diff(P, x, y, ydata, a=None):
        parse_args_ydata(P, x, y, ydata, a)
        getattr(cfitfun, name + '_diff')(byref(D))
        return y

    def fun_rms(P, x, y, a=None):
        parse_args(P, x, y, a)
        return getattr(cfitfun, name + '_rms')(byref(D))

    def fun_fit(P, x, y, ydata, a=None, do_var=None):
        parse_args_ydata(P, x, y, ydata, a)
        if do_var is not None:
            D.do_var = get_ptr(do_var)
        getattr(cfitfun, name + '_fit')(byref(D))
        return P

    return fun, fun_diff, fun_rms, fun_fit

e2, e2_diff, e2_rms, e2_fit  = _fun_factory('e2')
IV3, IV3_diff, IV3_rms, IV3_fit = _fun_factory('IV3')
IV4, IV4_diff, IV4_rms, IV4_fit = _fun_factory('IV4')
IV5, IV5_diff, IV5_rms, IV5_fit = _fun_factory('IV5')
IV6, IV6_diff, IV6_rms, IV6_fit = _fun_factory('IV6')
IVdbl, IVdbl_diff, IVdbl_rms, IVdbl_fit = _fun_factory('IVdbl')
IVdbl2, IVdbl2_diff, IVdbl2_rms, IVdbl2_fit = _fun_factory('IVdbl2')

