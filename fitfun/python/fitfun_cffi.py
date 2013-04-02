def get_ptr_array(x):
    return x.__array_interface__['data'][0]

try:
    from accel import _get_ptr as ptr
except ImportError:
    ptr = get_ptr_array


typedef = """\
typedef size_t p_t;

"""

header = """\
#include <fitfun.h>

data D;

void parse_args(int n, p_t P, int m, p_t x, p_t y) {
    D.n = n;
    D.P = (double*) P;
    D.m = m;
    D.x = (double*) x;
    D.y = (double*) y;
}

void parse_args_a(int n, p_t P, int m, p_t x, p_t y, p_t a) {
    parse_args(n, P, m, x, y);
    D.a = (double*) a;
}

void parse_args_ydata(int n, p_t P, int m, p_t x, p_t y, p_t ydata) {
    parse_args(n, P, m, x, y);
    D.ydata = (double*) ydata;
}

void parse_args_a_ydata(int n, p_t P, int m, p_t x, p_t y, p_t ydata, p_t a) {
    parse_args_a(n, P, m, x, y, a);
    D.ydata = (double*) ydata;
}


"""

template_header = """\
void _{fun}(int n, p_t P, int m, p_t x, p_t y);
void _{fun}_diff(int n, p_t P, int m, p_t x, p_t y, p_t ydata);
double _{fun}_rms(int n, p_t P, int m, p_t x, p_t y);
void _{fun}_fit(int n, p_t P, int m, p_t x, p_t y, p_t ydata, p_t do_var);

"""

template_header_a = """\
void _{fun}(int n, p_t P, int m, p_t x, p_t y, p_t a);
void _{fun}_diff(int n, p_t P, int m, p_t x, p_t y, p_t ydata, p_t a);
double _{fun}_rms(int n, p_t P, int m, p_t x, p_t y, p_t a);
void _{fun}_fit(int n, p_t P, int m, p_t x, p_t y, p_t ydata, p_t a, p_t do_var);

"""


template = """\
void _{fun}(int n, p_t P, int m, p_t x, p_t y) {left}
    parse_args(n, P, m, x, y);
    {fun}(&D);
{right}

void _{fun}_diff(int n, p_t P, int m, p_t x, p_t y, p_t ydata) {left}
    parse_args_ydata(n, P, m, x, y, ydata);
    {fun}_diff(&D);
{right}

double _{fun}_rms(int n, p_t P, int m, p_t x, p_t y) {left}
    parse_args(n, P, m, x, y);
    return {fun}_rms(&D);
{right}

void _{fun}_fit(int n, p_t P, int m, p_t x, p_t y, p_t ydata, p_t do_var) {left}
    parse_args_ydata(n, P, m, x, y, ydata);
    if (do_var) D.do_var = (int*) do_var;
    {fun}_fit(&D);
{right}


"""

template_a = """\
void _{fun}(int n, p_t P, int m, p_t x, p_t y, p_t a) {left}
    parse_args_a(n, P, m, x, y, a);
    {fun}(&D);
{right}

void _{fun}_diff(int n, p_t P, int m, p_t x, p_t y, p_t ydata, p_t a) {left}
    parse_args_a_ydata(n, P, m, x, y, ydata, a);
    {fun}_diff(&D);
{right}

double _{fun}_rms(int n, p_t P, int m, p_t x, p_t y, p_t a) {left}
    parse_args_a(n, P, m, x, y, a);
    return {fun}_rms(&D);
{right}

void _{fun}_fit(int n, p_t P, int m, p_t x, p_t y, p_t ydata, p_t a, p_t do_var) {left}
    parse_args_a_ydata(n, P, m, x, y, ydata, a);
    if (do_var) D.do_var = (int*) do_var;
    {fun}_fit(&D);
{right}


"""

funs = ('e2', 'IV3', 'IVdbl')
funs_a = ('IV4', 'IV5', 'IV6', 'IVdbl2')

cdef = typedef
verify = typedef + header

subst = dict(left='{', right='}')

for fun in funs:
    cdef += template_header.format(fun=fun, **subst)
    verify += template.format(fun=fun, **subst)

for fun in funs_a:
    cdef += template_header_a.format(fun=fun, **subst)
    verify += template_a.format(fun=fun, **subst)


from cffi import FFI
ffi = FFI()
ffi.cdef(cdef)

ff = ffi.verify(verify, libraries=['fitfun'])


def _fun_factory(name):
    name_diff, name_rms, name_fit = name + '_diff', name + '_rms', name + '_fit'

    def fun(P, x, y):
        getattr(ff, name)(P.size, ptr(P), 
                          x.size, ptr(x), ptr(y))
        return y

    def fun_diff(P, x, y, ydata):
        getattr(ff, name_diff)(P.size, ptr(P), 
                               x.size, ptr(x), ptr(y), ptr(ydata))
        return y

    def fun_rms(P, x, y):
        return getattr(ff, name_rms)(P.size, ptr(P), 
                                     x.size, ptr(x), ptr(y))

    def fun_fit(P, x, y, ydata, do_var=0):
        if do_var is not 0:
            do_var = ptr(do_var)
        getattr(ff, name_fit)(P.size, ptr(P), 
                              x.size, ptr(x), ptr(y), ptr(ydata), do_var)
        return P
    return fun, fun_diff, fun_rms, fun_fit


def _fun_a_factory(name):
    name_diff, name_rms, name_fit = name + '_diff', name + '_rms', name + '_fit'

    def fun(P, x, y, a):
        getattr(ff, name)(P.size, ptr(P), 
                          x.size, ptr(x), ptr(y), ptr(a))
        return y

    def fun_diff(P, x, y, ydata, a):
        getattr(ff, name_diff)(P.size, ptr(P), 
                               x.size, ptr(x), ptr(y), ptr(ydata), ptr(a))
        return y

    def fun_rms(P, x, y, a):
        return getattr(ff, name_rms)(P.size, ptr(P), 
                                     x.size, ptr(x), ptr(y), ptr(a))

    def fun_fit(P, x, y, ydata, a, do_var=0):
        if do_var is not 0:
            do_var = ptr(do_var)
        getattr(ff, name_fit)(P.size, ptr(P), 
                              x.size, ptr(x), ptr(y), ptr(ydata), ptr(a), do_var)
        return P
    return fun, fun_diff, fun_rms, fun_fit


e2, e2_diff, e2_rms, e2_fit  = _fun_factory('_e2')
IV3, IV3_diff, IV3_rms, IV3_fit = _fun_factory('_IV3')
IVdbl, IVdbl_diff, IVdbl_rms, IVdbl_fit = _fun_factory('_IVdbl')

IV4, IV4_diff, IV4_rms, IV4_fit = _fun_a_factory('_IV4')
IV5, IV5_diff, IV5_rms, IV5_fit = _fun_a_factory('_IV5')
IV6, IV6_diff, IV6_rms, IV6_fit = _fun_a_factory('_IV6')
IVdbl2, IVdbl2_diff, IVdbl2_rms, IVdbl2_fit = _fun_a_factory('_IVdbl2')


if __name__ == "__main__":
    try:
        import numpypy as np
    except ImportError:
        import numpy as np

    P = np.array([1., 1.])
    x = np.arange(4.)
    y = np.empty(x.shape, x.dtype)

    print e2(P, x, y)

    print P[0]*np.exp(-P[1]*x)


