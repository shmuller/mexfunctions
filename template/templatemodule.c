#include <Python.h>
#include <numpy/arrayobject.h>

static PyObject* test(PyObject *self, PyObject *args)
{
    Py_RETURN_NONE;
}


static PyMethodDef methods[] = {
    {"test", test, METH_VARARGS, "Test function"},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
inittemplate(void)
{
    Py_InitModule("template", methods);
    import_array();
}


