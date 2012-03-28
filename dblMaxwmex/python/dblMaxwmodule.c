/* triblock(A,B,IPIV,C,D)
 * Tri-block-diagnonal solve step. Compile with:
 *
 * python setup.py install
 *
 * S. H. Muller, 2012/03/04
 */

#include <Python.h>
#include <numpy/arrayobject.h>

#include "dblMaxw.h"

static PyObject *integrate(PyObject *self, PyObject *args)
{
    /*
    int is_half = 0;
    PyObject *a[5];
    if (!PyArg_ParseTuple(args, "OOOOO|i", a, a+1, a+2, a+3, a+4, &is_half)) {
        PyErr_SetString(PyExc_TypeError, "D1, U2, IPIV, L1, D2 expected");
        return NULL;
    }

	double *D1 = PyArray_DATA(a[0]);
	double *U2 = PyArray_DATA(a[1]);
	int *IPIV  = PyArray_DATA(a[2]);
	double *L1 = PyArray_DATA(a[3]);
	double *D2 = PyArray_DATA(a[4]);

	int N = PyArray_DIM(a[1],1), NRHS = PyArray_DIM(a[1],0);
    int INFO;

    _triblock(D1, U2, IPIV, L1, D2, is_half, N, NRHS, &INFO);

    */
    Py_RETURN_NONE;
}


static PyMethodDef methods[] = {
    {"integrate", integrate, METH_VARARGS, "Double-Maxwellian integral"},
    {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
initdblMaxw(void)
{
    Py_InitModule("dblMaxw", methods);
}

