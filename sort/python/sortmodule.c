/* sort wrapper for Python. Compile with:
 * 
 * python setup.py install
 *
 * S. H. Muller, 2014/02/04
 */

#include <Python.h>
#include <numpy/arrayobject.h>

#include "../sort.h"

static PyObject* sort(PyObject *self, PyArrayObject *arr)
{
    int ndims = PyArray_NDIM(arr);
    npy_intp N = PyArray_DIM(arr, ndims-1);
    int typenum = PyArray_TYPE(arr);
    void *x = PyArray_DATA(arr);

    quick_sort(x, N);

    Py_RETURN_NONE;
}

static PyObject* qsort_meth(PyObject *self, PyArrayObject *arr)
{
    int ndims = PyArray_NDIM(arr);
    npy_intp N = PyArray_DIM(arr, ndims-1);
    int typenum = PyArray_TYPE(arr);
    void *x = PyArray_DATA(arr);

    qsort_wrap(x, N);

    Py_RETURN_NONE;
}

static PyObject* quicksort_meth(PyObject *self, PyArrayObject *arr)
{
    int ndims = PyArray_NDIM(arr);
    npy_intp N = PyArray_DIM(arr, ndims-1);
    int typenum = PyArray_TYPE(arr);
    void *x = PyArray_DATA(arr);

    quicksort(x, N);

    Py_RETURN_NONE;
}

static PyObject* stdsort_meth(PyObject *self, PyArrayObject *arr)
{
    int ndims = PyArray_NDIM(arr);
    npy_intp N = PyArray_DIM(arr, ndims-1);
    int typenum = PyArray_TYPE(arr);
    void *x = PyArray_DATA(arr);

    stdsort(x, N);

    Py_RETURN_NONE;
}


static PyMethodDef methods[] = {
    {"qsort", (PyCFunction)qsort_meth, METH_O, "qsort"},
    {"sort", (PyCFunction)sort, METH_O, "Sort"},
    {"quicksort", (PyCFunction)quicksort_meth, METH_O, "Quicksort"},
    {"stdsort", (PyCFunction)stdsort_meth, METH_O, "std::sort"},
    {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
initsort(void)
{
    Py_InitModule("sort", methods);
    import_array();
}

