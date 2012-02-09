#include <Python.h>
#include "numpy/noprefix.h"

#include "atomic.h"

static PyObject* atomic(PyObject *self, PyObject *args)
{
	int i;
	PyObject *pyx;
	char *target=NULL, *type=NULL, *model=NULL;
	
	if (!PyArg_ParseTuple(args, "Ozz|z", &pyx, &target, &type, &model)) {
        return NULL;
    }
    
    int n = PyArray_DIM(pyx, 0);
	double *x = PyArray_DATA(pyx);
	
	atomic_desc D = get_atomic_desc(target, type, model);
	
    sigmav_vec(n, x, &D);
    
    Py_RETURN_NONE;
}



static PyMethodDef AtomicMethods[] = {
    {"atomic", atomic, METH_VARARGS, "Atomic properties"},
    {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
initatomic(void)
{
    Py_InitModule("atomic", AtomicMethods);
}
