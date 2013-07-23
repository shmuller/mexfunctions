import py

from pypy.module.cpyext.pyobject import PyObject
from pypy.module.cpyext.test.test_api import BaseApiTest
from rpython.rtyper.lltypesystem import rffi, lltype

from pypy.module.micronumpy.interp_numarray import W_NDimArray
from pypy.module.micronumpy.interp_dtype import get_dtype_cache
from pypy.objspace.std.listobject import W_ListObject

def scalar(space):
    dtype = get_dtype_cache(space).w_float64dtype
    return W_NDimArray.new_scalar(space, dtype, space.wrap(10.))

def array(space, shape):
    dtype = get_dtype_cache(space).w_float64dtype
    return W_NDimArray.from_shape(shape, dtype, order='C')

def iarray(space, shape):
    dtype = get_dtype_cache(space).w_int64dtype
    return W_NDimArray.from_shape(shape, dtype, order='C')


NULL = lltype.nullptr(rffi.VOIDP.TO)

class TestNDArrayObject(BaseApiTest):

    def test_NDIM(self, space, api):
        a = array(space, [10, 5, 3])
        assert api.PyArray_NDIM(a) == 3

    def test_DIM(self, space, api):
        a = array(space, [10, 5, 3])
        assert api.PyArray_DIM(a, 1) == 5

    def test_SIZE(self, space, api):
        a = array(space, [10, 5, 3])
        assert api.PyArray_SIZE(a) == 150

    def test_ITEMSIZE(self, space, api):
        a = array(space, [10, 5, 3])
        assert api.PyArray_ITEMSIZE(a) == 8
    
    def test_NBYTES(self, space, api):
        a = array(space, [10, 5, 3])
        assert api.PyArray_NBYTES(a) == 1200

    def test_TYPE(self, space, api):
        a = array(space, [10, 5, 3])
        assert api.PyArray_TYPE(a) == 12

    def test_DATA(self, space, api):
        a = array(space, [10, 5, 3])
        addr = api.PyArray_DATA(a)
        addr2 = rffi.cast(rffi.VOIDP, a.implementation.storage)
        assert addr == addr2

    def test_FromAny_scalar(self, space, api):
        a0 = scalar(space)
        assert a0.implementation.get_scalar_value().value == 10.

        a = api.PyArray_FromAny(a0, NULL, 0, 0, 0, NULL)
        assert api.PyArray_NDIM(a) == 0
        
        ptr = rffi.cast(rffi.DOUBLEP, api.PyArray_DATA(a))
        assert ptr[0] == 10.

    def test_FromAny(self, space, api):
        a = array(space, [10, 5, 3])
        assert api.PyArray_FromAny(a, NULL, 0, 0, 0, NULL) is a

    def test_list_from_fixedptr(self, space, api):
        A = lltype.GcArray(lltype.Float)
        ptr = lltype.malloc(A, 3)
        assert isinstance(ptr, lltype._ptr)
        ptr[0] = 10.
        ptr[1] = 5.
        ptr[2] = 3.
        l = list(ptr)
        assert l == [10., 5., 3.]

    def test_list_from_openptr(self, space, api):
        nd = 3
        a = array(space, [nd])
        ptr = rffi.cast(rffi.DOUBLEP, api.PyArray_DATA(a))
        ptr[0] = 10.
        ptr[1] = 5.
        ptr[2] = 3.
        l = []
        for i in range(nd):
            l.append(ptr[i])
        assert l == [10., 5., 3.]

    def test_SimpleNew_scalar(self, space, api):
        ptr_s = lltype.nullptr(rffi.LONGP.TO)
        a = api.PyArray_SimpleNew(0, ptr_s, 12)
        
        dtype = get_dtype_cache(space).w_float64dtype
        
        a.set_scalar_value(dtype.itemtype.box(10.))
        assert a.get_scalar_value().value == 10.

    def test_SimpleNewFromData_scalar(self, space, api):
        a = array(space, [1])
        num = api.PyArray_TYPE(a)
        ptr_a = api.PyArray_DATA(a)

        x = rffi.cast(rffi.DOUBLEP, ptr_a)
        x[0] = float(10.)

        ptr_s = lltype.nullptr(rffi.LONGP.TO)
        
        res = api.PyArray_SimpleNewFromData(0, ptr_s, num, ptr_a)
        assert res.is_scalar()
        assert res.get_scalar_value().value == 10.

    def test_SimpleNew(self, space, api):
        shape = [10, 5, 3]
        nd = len(shape)

        s = iarray(space, [nd])
        ptr_s = rffi.cast(rffi.LONGP, api.PyArray_DATA(s))
        ptr_s[0] = 10
        ptr_s[1] = 5
        ptr_s[2] = 3
        
        a = api.PyArray_SimpleNew(nd, ptr_s, 12)
        
        #assert list(api.PyArray_DIMS(a))[:3] == shape

        ptr_a = api.PyArray_DATA(a)

        x = rffi.cast(rffi.DOUBLEP, ptr_a)
        for i in range(150):
            x[i] = float(i)

        for i in range(150):
            assert x[i] == float(i)

    def test_SimpleNewFromData(self, space, api):
        shape = [10, 5, 3]
        nd = len(shape)

        s = iarray(space, [nd])
        ptr_s = rffi.cast(rffi.LONGP, api.PyArray_DATA(s))
        ptr_s[0] = 10
        ptr_s[1] = 5
        ptr_s[2] = 3

        a = array(space, shape)
        num = api.PyArray_TYPE(a)
        ptr_a = api.PyArray_DATA(a)

        x = rffi.cast(rffi.DOUBLEP, ptr_a)
        for i in range(150):
            x[i] = float(i)

        res = api.PyArray_SimpleNewFromData(nd, ptr_s, num, ptr_a)
        assert api.PyArray_TYPE(res) == num
        assert api.PyArray_DATA(res) == ptr_a
        for i in range(nd):
            assert api.PyArray_DIM(res, i) == shape[i]
        ptr_r = rffi.cast(rffi.DOUBLEP, api.PyArray_DATA(res))
        for i in range(150):
            assert ptr_r[i] == float(i)

    def test_SimpleNewFromData_complex(self, space, api):
        a = array(space, [2])
        ptr_a = api.PyArray_DATA(a)
        
        x = rffi.cast(rffi.DOUBLEP, ptr_a)
        x[0] = 3.
        x[1] = 4.
                
        ptr_s = lltype.nullptr(rffi.LONGP.TO)
        
        res = api.PyArray_SimpleNewFromData(0, ptr_s, 15, ptr_a)
        assert res.get_scalar_value().real == 3.
        assert res.get_scalar_value().imag == 4.
