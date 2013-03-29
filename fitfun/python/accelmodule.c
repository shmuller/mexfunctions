#include <Python.h>
#include <numpy/arrayobject.h>

PyObject *_get_ptr(PyObject *self, PyObject *obj) {
    return PyLong_FromVoidPtr(PyArray_DATA(obj));
}

static PyMethodDef methods[] = {
    {"_get_ptr", _get_ptr, METH_O, "Wrapper to PyArray_DATA()"},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initaccel(void)
{
    Py_InitModule("accel", methods);
}

