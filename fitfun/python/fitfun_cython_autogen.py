template_header_pxd = """\
cdef extern from "{header}":
    ctypedef struct data:
        int n
        int m
        {dtype}* P
        int* do_var
        {dtype}* x
        {dtype}* y
        {dtype}* ydata
        {dtype}* w
        {dtype}* a

"""

template_pxd = """\
    void {fun}(data* D)
    void {fun}_diff(data* D)
    {dtype} {fun}_rms(data* D)
    void {fun}_fit(data *D)

"""

template_header = """\
import cython
from numpy cimport ndarray, {dtype}_t, int32_t
ctypedef {dtype}_t DTYPE_t

cimport {cfitfun}

cdef {cfitfun}.data D

cdef parse_args({p_t} P, {p_t} x, {p_t} y):
    D.n = P.size
    D.m = x.size
    D.P = <{dtype}*> P.data
    D.x = <{dtype}*> x.data
    D.y = <{dtype}*> y.data

cdef parse_args_a({p_t} P, {p_t} x, {p_t} y, {p_t} a):
    parse_args(P, x, y)
    D.a = <{dtype}*> a.data

cdef parse_args_ydata({p_t} P, {p_t} x, {p_t} y, {p_t} ydata):
    parse_args(P, x, y)
    D.ydata = <{dtype}*> ydata.data

cdef parse_args_a_ydata({p_t} P, {p_t} x, {p_t} y, {p_t} ydata, {p_t} a):
    parse_args_a(P, x, y, a)
    D.ydata = <{dtype}*> ydata.data


"""

template = """\
cpdef {fun}({p_t} P, {p_t} x, {p_t} y):
    parse_args(P, x, y)
    {cfitfun}.{fun}(&D)
    return y

cpdef {fun}_diff({p_t} P, {p_t} x, {p_t} y, {p_t} ydata):
    parse_args_ydata(P, x, y, ydata)
    {cfitfun}.{fun}_diff(&D)
    return y

cpdef {fun}_rms({p_t} P, {p_t} x, {p_t} y):
    parse_args(P, x, y)
    return {cfitfun}.{fun}_rms(&D)

cpdef {fun}_fit({p_t} P, {p_t} x, {p_t} y, {p_t} ydata, ndarray do_var=None):
    parse_args_ydata(P, x, y, ydata)
    if do_var is not None:
        D.do_var = <int*> do_var.data
    {cfitfun}.{fun}_fit(&D)
    return P


"""

template_a = """\
cpdef {fun}({p_t} P, {p_t} x, {p_t} y, {p_t} a):
    parse_args_a(P, x, y, a)
    {cfitfun}.{fun}(&D)
    return y

cpdef {fun}_diff({p_t} P, {p_t} x, {p_t} y, {p_t} ydata, {p_t} a):
    parse_args_a_ydata(P, x, y, ydata, a)
    {cfitfun}.{fun}_diff(&D)
    return y

cpdef {fun}_rms({p_t} P, {p_t} x, {p_t} y, {p_t} a):
    parse_args_a(P, x, y, a)
    return {cfitfun}.{fun}_rms(&D)

cpdef {fun}_fit({p_t} P, {p_t} x, {p_t} y, {p_t} ydata, {p_t} a, ndarray do_var=None):
    parse_args_a_ydata(P, x, y, ydata, a)
    if do_var is not None:
        D.do_var = <int*> do_var.data
    {cfitfun}.{fun}_fit(&D)
    return P


"""

class AutoGen:
    def __init__(self, header='fitfun.h', pyx_file='fitfun_cython.pyx', 
            funs=(), funs_a=()):
        self.cfitfun = 'c' + header.rpartition('.h')[0]

        self.pxd_file = self.cfitfun + '.pxd'
        self.pyx_file = pyx_file

        self.funs = funs
        self.funs_a = funs_a
        self.funs_all = funs + funs_a

        self.subst = dict(
            header = header,
            cfitfun = self.cfitfun,
            p_t = "ndarray[DTYPE_t]",
            dtype = "double")
        

    def pxd_gen(self):
        f = open(self.pxd_file, "w")
        f.write(template_header_pxd.format(**self.subst))
        for fun in self.funs_all:
            f.write(template_pxd.format(fun=fun, **self.subst))
        f.close()

    def pyx_gen(self):
        f = open(self.pyx_file, "w")
        f.write(template_header.format(**self.subst))
        for fun in self.funs:
            f.write(template.format(fun=fun, **self.subst))
        for fun in self.funs_a:
            f.write(template_a.format(fun=fun, **self.subst))
        f.close()


if __name__ == "__main__":
    funs = ('e2', 'IV3', 'IVdbl')
    funs_a = ('IV4', 'IV5', 'IV6', 'IVdbl2')

    autogen = AutoGen(funs=funs, funs_a=funs_a)
    autogen.pxd_gen()
    autogen.pyx_gen()


