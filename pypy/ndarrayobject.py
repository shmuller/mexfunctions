from pypy.interpreter.error import OperationError
from rpython.rtyper.lltypesystem import rffi, lltype
from pypy.module.cpyext.api import cpython_api, Py_ssize_t, CANNOT_FAIL
from pypy.module.cpyext.pyerrors import PyErr_BadInternalCall
from pypy.module.cpyext.pyobject import PyObject
from pypy.module.micronumpy.interp_numarray import W_NDimArray, convert_to_array
from pypy.module.micronumpy.interp_dtype import get_dtype_cache
from pypy.module.micronumpy.arrayimpl.scalar import Scalar
from rpython.rlib.rawstorage import RAW_STORAGE_PTR, raw_storage_getitem
from pypy.objspace.std.tupleobject import W_TupleObject
from pypy.module.cpyext.pyobject import Py_DecRef

# the asserts are needed, otherwise the translation fails

@cpython_api([PyObject], Py_ssize_t, error=CANNOT_FAIL)
def PyArray_NDIM(space, w_array):
    assert isinstance(w_array, W_NDimArray)
    return len(w_array.get_shape())

@cpython_api([PyObject, Py_ssize_t], Py_ssize_t, error=CANNOT_FAIL)
def PyArray_DIM(space, w_array, n):
    assert isinstance(w_array, W_NDimArray)
    return w_array.get_shape()[n]

@cpython_api([PyObject], Py_ssize_t, error=CANNOT_FAIL)
def PyArray_SIZE(space, w_array):
    assert isinstance(w_array, W_NDimArray)
    return w_array.get_size()

@cpython_api([PyObject], Py_ssize_t, error=CANNOT_FAIL)
def PyArray_ITEMSIZE(space, w_array):
    assert isinstance(w_array, W_NDimArray)
    return w_array.get_dtype().get_size()

@cpython_api([PyObject], Py_ssize_t, error=CANNOT_FAIL)
def PyArray_NBYTES(space, w_array):
    assert isinstance(w_array, W_NDimArray)
    return w_array.get_size() * w_array.get_dtype().get_size()

@cpython_api([PyObject], Py_ssize_t, error=CANNOT_FAIL)
def PyArray_TYPE(space, w_array):
    assert isinstance(w_array, W_NDimArray)
    return w_array.get_dtype().num


@cpython_api([PyObject], rffi.VOIDP, error=CANNOT_FAIL)
def PyArray_DATA(space, w_array):
    assert isinstance(w_array, W_NDimArray)
    return rffi.cast(rffi.VOIDP, w_array.implementation.storage)


@cpython_api([PyObject, rffi.VOIDP, Py_ssize_t, Py_ssize_t, Py_ssize_t, rffi.VOIDP], 
             PyObject, error=CANNOT_FAIL)
def PyArray_FromAny(space, w_obj, dtype, min_depth, max_depth, requirements, context):
    w_array = convert_to_array(space, w_obj)
    if w_array.is_scalar():
        impl = w_array.implementation
        w_array = W_NDimArray.from_shape([1], impl.dtype)
        w_array.implementation.setitem(0, impl.value)
        w_array.implementation.shape = []
    return w_array


@cpython_api([Py_ssize_t, rffi.LONGP, Py_ssize_t], PyObject, error=CANNOT_FAIL)
def PyArray_SimpleNew(space, nd, dims, typenum):
    dtype = get_dtype_cache(space).dtypes_by_num[typenum]
    shape = []
    for i in range(nd):
        # back-and-forth wrapping needed to translate
        shape.append(space.int_w(space.wrap(dims[i])))

    return W_NDimArray.from_shape(shape, dtype)


def simple_new_from_data(space, nd, dims, typenum, data, owning):
    dtype = get_dtype_cache(space).dtypes_by_num[typenum]
    storage = rffi.cast(RAW_STORAGE_PTR, data)
    if nd == 0:
        w_val = dtype.itemtype.box_raw_data(storage)
        return W_NDimArray(Scalar(dtype, w_val))
    else:
        shape = []
        for i in range(nd):
            # back-and-forth wrapping needed to translate
            shape.append(space.int_w(space.wrap(dims[i])))
        
        return W_NDimArray.from_shape_and_storage(shape, storage, dtype, owning=owning)

@cpython_api([Py_ssize_t, rffi.LONGP, Py_ssize_t, rffi.VOIDP], PyObject, error=CANNOT_FAIL)
def PyArray_SimpleNewFromData(space, nd, dims, typenum, data):
    return simple_new_from_data(space, nd, dims, typenum, data, owning=False)

@cpython_api([Py_ssize_t, rffi.LONGP, Py_ssize_t, rffi.VOIDP], PyObject, error=CANNOT_FAIL)
def PyArray_SimpleNewFromDataOwning(space, nd, dims, typenum, data):
    return simple_new_from_data(space, nd, dims, typenum, data, owning=True)

"""
#shape_array = lltype.Array(rffi.LONG)
#shape_ptr = lltype.malloc(shape_array, 10, immortal=True)

@cpython_api([PyObject], rffi.LONGP, error=CANNOT_FAIL)
def PyArray_DIMS(space, w_array):
    assert isinstance(w_array, W_NDimArray)
    if w_array.is_scalar():
        return lltype.nullptr(rffi.LONGP.TO)
    else:
        shape = w_array.get_shape()
        dtype = get_dtype_cache(space).w_int64dtype
        w_shape = W_NDimArray.from_shape([len(shape)], dtype)
        shape_ptr = rffi.cast(rffi.LONGP, w_shape.implementation.storage)
        for i in range(len(shape)):
            shape_ptr[i] = shape[i]
        return shape_ptr
"""

