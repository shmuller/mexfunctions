#include <Python.h>
#include <numpy/arrayobject.h>

PyObject *get_ptr(PyObject *self, PyObject *obj) {
    return Py_BuildValue("N", PyLong_FromVoidPtr(PyArray_DATA(obj)));
}

static PyMethodDef methods[] = {
    {"get_ptr", get_ptr, METH_O, "Wrapper to PyArray_DATA()"},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initaccel(void)
{
    Py_InitModule("accel", methods);
}

