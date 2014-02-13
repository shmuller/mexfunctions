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

    quick_sort((int*)x, N);

    Py_RETURN_NONE;
}

static PyObject* qsort_meth(PyObject *self, PyArrayObject *arr)
{
    int ndims = PyArray_NDIM(arr);
    npy_intp N = PyArray_DIM(arr, ndims-1);
    int typenum = PyArray_TYPE(arr);
    void *x = PyArray_DATA(arr);

    qsort_wrap((int*)x, N);

    Py_RETURN_NONE;
}

static PyObject* numpy_quicksort_meth(PyObject *self, PyArrayObject *arr)
{
    int ndims = PyArray_NDIM(arr);
    npy_intp N = PyArray_DIM(arr, ndims-1);
    int typenum = PyArray_TYPE(arr);
    void *x = PyArray_DATA(arr);

    numpy_quicksort((int*)x, N);

    Py_RETURN_NONE;
}

static PyObject* numpy_mergesort_meth(PyObject *self, PyArrayObject *arr)
{
    int ndims = PyArray_NDIM(arr);
    npy_intp N = PyArray_DIM(arr, ndims-1);
    int typenum = PyArray_TYPE(arr);
    void *x = PyArray_DATA(arr);

    numpy_mergesort((int*)x, N);

    Py_RETURN_NONE;
}

static PyObject* insert_sort_meth(PyObject *self, PyArrayObject *arr)
{
    int ndims = PyArray_NDIM(arr);
    npy_intp N = PyArray_DIM(arr, ndims-1);
    int typenum = PyArray_TYPE(arr);
    void *x = PyArray_DATA(arr);

    insert_sort((int*)x, N);

    Py_RETURN_NONE;
}

static PyObject* select_sort_meth(PyObject *self, PyArrayObject *arr)
{
    int ndims = PyArray_NDIM(arr);
    npy_intp N = PyArray_DIM(arr, ndims-1);
    int typenum = PyArray_TYPE(arr);
    void *x = PyArray_DATA(arr);

    select_sort((int*)x, N);

    Py_RETURN_NONE;
}

static PyObject* merge_sort_meth(PyObject *self, PyArrayObject *arr)
{
    int ndims = PyArray_NDIM(arr);
    npy_intp N = PyArray_DIM(arr, ndims-1);
    int typenum = PyArray_TYPE(arr);
    void *x = PyArray_DATA(arr);

    merge_sort((int*)x, N);

    Py_RETURN_NONE;
}

#include <algorithm>

static PyObject* stdsort_meth(PyObject *self, PyArrayObject *arr)
{
    int ndims = PyArray_NDIM(arr);
    npy_intp N = PyArray_DIM(arr, ndims-1);
    int typenum = PyArray_TYPE(arr);
    void *x = PyArray_DATA(arr);

    std::sort((int*)x, (int*)x + N);

    Py_RETURN_NONE;
}

#include "../timsort.hpp"

static PyObject* timsort_meth(PyObject *self, PyArrayObject *arr)
{
    int ndims = PyArray_NDIM(arr);
    npy_intp N = PyArray_DIM(arr, ndims-1);
    int typenum = PyArray_TYPE(arr);
    void *x = PyArray_DATA(arr);

    gfx::timsort((int*)x, (int*)x + N);

    Py_RETURN_NONE;
}

#define SORT_NAME swenson
#define SORT_TYPE int

#include "swenson_sort.h"

static PyObject* sw_timsort_meth(PyObject *self, PyArrayObject *arr)
{
    int ndims = PyArray_NDIM(arr);
    npy_intp N = PyArray_DIM(arr, ndims-1);
    int typenum = PyArray_TYPE(arr);
    void *x = PyArray_DATA(arr);

    swenson_tim_sort((int*)x, N);

    Py_RETURN_NONE;
}

static PyObject* sw_quicksort_meth(PyObject *self, PyArrayObject *arr)
{
    int ndims = PyArray_NDIM(arr);
    npy_intp N = PyArray_DIM(arr, ndims-1);
    int typenum = PyArray_TYPE(arr);
    void *x = PyArray_DATA(arr);

    swenson_quick_sort((int*)x, N);

    Py_RETURN_NONE;
}

static PyObject* bitmap_sort_meth(PyObject *self, PyArrayObject *arr)
{
    int ndims = PyArray_NDIM(arr);
    npy_intp N = PyArray_DIM(arr, ndims-1);
    int typenum = PyArray_TYPE(arr);
    void *x = PyArray_DATA(arr);

    int ncount = bitmap_sort((unsigned int*)x, N, NULL);

    return Py_BuildValue("i", ncount);
}


static PyMethodDef methods[] = {
    {"qsort", (PyCFunction)qsort_meth, METH_O, "qsort"},
    {"sort", (PyCFunction)sort, METH_O, "Sort"},
    {"numpy_quicksort", (PyCFunction)numpy_quicksort_meth, METH_O, "Numpy quicksort"},
    {"numpy_mergesort", (PyCFunction)numpy_mergesort_meth, METH_O, "Numpy mergesort"},
    {"insert_sort", (PyCFunction)insert_sort_meth, METH_O, "Insertion sort"},
    {"select_sort", (PyCFunction)select_sort_meth, METH_O, "Selection sort"},
    {"merge_sort", (PyCFunction)merge_sort_meth, METH_O, "Mergesort"},
    {"stdsort", (PyCFunction)stdsort_meth, METH_O, "std::sort"},
    {"timsort", (PyCFunction)timsort_meth, METH_O, "gfx::timsort"},
    {"sw_timsort", (PyCFunction)sw_timsort_meth, METH_O, "Swenson timsort"},
    {"sw_quicksort", (PyCFunction)sw_quicksort_meth, METH_O, "Swenson quicksort"},
    {"bitmap_sort", (PyCFunction)bitmap_sort_meth, METH_O, "Bitmap sort"},
    {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
initsort(void)
{
    Py_InitModule("sort", methods);
    import_array();
}

