#include <Python.h>
#include "numpy/noprefix.h"

#include "atomic.h"

static PyObject* atomic(PyObject *self, PyObject *args)
{
	int i;
	PyObject *pyx;
	
	if (!PyArg_ParseTuple(args, "O", &pyx)) {
        return NULL;
    }
    
    int n = PyArray_DIM(pyx, 0);
	double *x = PyArray_DATA(pyx);
	
	atomic_desc D = get_atomic_desc("D", "ion", "BEB");
	
    sigmav_vec(n, x, &D);
    
	
	return Py_BuildValue("d", D.mu_2q);
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
