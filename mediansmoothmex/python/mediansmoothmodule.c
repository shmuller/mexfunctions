/* mediansmooth wrapper for Python. Compile with:
 * 
 * python setup.py install
 *
 * S. H. Muller, 2012/09/11
 */

#include <Python.h>
#include <numpy/arrayobject.h>

#include "../mediansmooth.h"

static PyObject* mediansmooth(PyObject *self, PyObject *args)
{
    PyObject *in;
    int w, bdry = 0;
    if (!PyArg_ParseTuple(args, "Oi|i", &in, &w, &bdry)) {
        PyErr_SetString(PyExc_TypeError, "x, w, [bdry] expected");
        return NULL;
    }
    PyObject *arr = PyArray_FromAny(in, NULL, 0, 0, 0, NULL);
    
    int ndims = PyArray_NDIM(arr);
    npy_intp N = PyArray_DIM(arr, ndims-1);
    int typenum = PyArray_TYPE(arr);
    void *x = PyArray_DATA(arr);

    switch (typenum) {
        case NPY_DOUBLE:
            median_filt_pqueue_double(x, N, w, bdry);
            break;
        case NPY_FLOAT:
            median_filt_pqueue_float(x, N, w, bdry);
            break;
        case NPY_LONG:
            median_filt_pqueue_int64(x, N, w, bdry);
            break;
        case NPY_INT:
            median_filt_pqueue_int32(x, N, w, bdry);
            break;
    }
    Py_RETURN_NONE;
}

static PyObject* mediansmooth_old(PyObject *self, PyObject *args)
{
    PyObject *in;
    int w;
    if (!PyArg_ParseTuple(args, "Oi", &in, &w)) {
        PyErr_SetString(PyExc_TypeError, "x, w expected");
        return NULL;
    }
    PyObject *arr = PyArray_FromAny(in, NULL, 0, 0, 0, NULL);
    
    int ndims = PyArray_NDIM(arr);
    npy_intp N = PyArray_DIM(arr, ndims-1);
    int typenum = PyArray_TYPE(arr);
    void *x = PyArray_DATA(arr);

    PyObject *retval = PyArray_SimpleNew(1, &N, NPY_INT);
    int *ind = PyArray_DATA(retval);

    switch (typenum) {
        case NPY_DOUBLE:
            median_filt_double(x, N, w, ind);
            break;
        case NPY_FLOAT:
            median_filt_float(x, N, w, ind);
            break;
    }
    return retval;
}


static PyMethodDef methods[] = {
    {"mediansmooth", mediansmooth, METH_VARARGS, "Sliding median filter"},
    {"mediansmooth_old", mediansmooth_old, METH_VARARGS, "Sliding median filter - old"},
    {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
initmediansmooth(void)
{
    Py_InitModule("mediansmooth", methods);
    import_array();
}

PyMODINIT_FUNC
initmediansmooth2(void)
{
    Py_InitModule("mediansmooth2", methods);
    import_array();
}


