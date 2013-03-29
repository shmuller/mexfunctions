import ctypes as C
import ctypes.util

p_t = C.POINTER(C.c_double)
p_int_t = C.POINTER(C.c_int)

class struct_data(C.Structure):
    __slots__ = ['n', 'm', 'P', 'do_var', 'x', 'y', 'ydata', 'w', 'a']

    _fields_ = [
        ('n', C.c_int),
        ('m', C.c_int),
        ('P', p_t),
        ('do_var', p_int_t),
        ('x', p_t),
        ('y', p_t),
        ('ydata', p_t),
        ('w', p_t),
        ('a', p_t)]

D = struct_data()

def get_ptr_ctypes(x, t=p_t):
    return x.ctypes.data_as(t)

def get_ptr_array(x, t=p_t):
    return C.cast(x.__array_interface__['data'][0], t)

try:
    from accel import _get_ptr
    def get_ptr(x, t=p_t):
        return C.cast(_get_ptr(x), t)
except ImportError:
    get_ptr = get_ptr_array


def parse_args(P, x, y, a=None):
    D.P = get_ptr(P)
    D.x = get_ptr(x)
    D.y = get_ptr(y)
    D.n = C.c_int(P.size)
    D.m = C.c_int(x.size)
    if a is not None:
        D.a = get_ptr(a)

def parse_args_ydata(P, x, y, ydata, a=None):
    parse_args(P, x, y, a)
    D.ydata = get_ptr(ydata)


libname = C.util.find_library('fitfun')

cfitfun = C.CDLL(libname)

def _fun_factory(name):
    for suffix in ('', '_diff', '_rms', '_fit'):
        getattr(cfitfun, name + suffix).argtypes = [C.POINTER(struct_data)]

    getattr(cfitfun, name + '_rms').restype = C.c_double

    def fun(P, x, y, a=None):
        parse_args(P, x, y, a)
        getattr(cfitfun, name)(C.byref(D))
        return y

    def fun_diff(P, x, y, ydata, a=None):
        parse_args_ydata(P, x, y, ydata, a)
        getattr(cfitfun, name + '_diff')(C.byref(D))
        return y

    def fun_rms(P, x, y, a=None):
        parse_args(P, x, y, a)
        return getattr(cfitfun, name + '_rms')(C.byref(D))

    def fun_fit(P, x, y, ydata, a=None, do_var=None):
        parse_args_ydata(P, x, y, ydata, a)
        if do_var is not None:
            D.do_var = get_ptr(do_var, p_int_t)
        getattr(cfitfun, name + '_fit')(C.byref(D))
        return P

    return fun, fun_diff, fun_rms, fun_fit

e2, e2_diff, e2_rms, e2_fit  = _fun_factory('e2')
IV3, IV3_diff, IV3_rms, IV3_fit = _fun_factory('IV3')
IV4, IV4_diff, IV4_rms, IV4_fit = _fun_factory('IV4')
IV5, IV5_diff, IV5_rms, IV5_fit = _fun_factory('IV5')
IV6, IV6_diff, IV6_rms, IV6_fit = _fun_factory('IV6')
IVdbl, IVdbl_diff, IVdbl_rms, IVdbl_fit = _fun_factory('IVdbl')
IVdbl2, IVdbl2_diff, IVdbl2_rms, IVdbl2_fit = _fun_factory('IVdbl2')

