#include <Python.h>
#include "numpy/arrayobject.h"

#include "../specfun.h"

static PyObject* specfun(PyObject *self, PyObject *args)
{
	int i;
	PyObject *pyx;
	char *name=NULL;
	
	if (!PyArg_ParseTuple(args, "zO", &name, &pyx)) {
        return NULL;
    }
    
    int n = PyArray_DIM(pyx, 0);
	double *x = PyArray_DATA(pyx);
	
	func* fun = kv_select(KV_LEN(kv_specfun), kv_specfun, name);
	if (fun==NULL) {
		return NULL;
	}
	
	for (i=n; i--; )
	    *x++ = fun(x);
    
    Py_RETURN_NONE;
}

static PyMethodDef SpecfunMethods[] = {
    {"specfun", specfun, METH_VARARGS, "Special functions"},
    {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
initspecfun(void)
{
    Py_InitModule("specfun", SpecfunMethods);
}
