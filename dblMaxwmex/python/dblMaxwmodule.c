/* dblMaxw.integrate("f_r",vt,v0,ut,u0,IJ,(VU))
 * Double-Maxwellian integrals. Compile with:
 *
 * python setup.py install
 *
 * S. H. Muller, 2012/03/29
 */

#include <Python.h>
#include <numpy/arrayobject.h>

#include "dblMaxw.h"

static PyObject *integrate(PyObject *self, PyObject *args)
{
    int i;
    char *f_r, *nrm;
    double vt, ut;
    PyObject *a[3], *VU_args=NULL, *VU_arg, *L;
    if (!PyArg_ParseTuple(args, "zdOdOOz|O", &f_r, &vt, a, &ut, a+1, a+2, &nrm, &VU_args)) {
        PyErr_SetString(PyExc_TypeError, "f_r, vt, v0, ut, u0, IJ, nrm, (VU) expected");
        return NULL;
    }

    double *v0 = PyArray_DATA(a[0]);
    double *u0 = PyArray_DATA(a[1]);
    int *IJ    = PyArray_DATA(a[2]);
    
    int mV = 3-IJ[0], mU = 3-IJ[1], m = mV+mU;

    if (((VU_args) ? PyTuple_Size(VU_args) : 0) != m) {
        PyErr_SetString(PyExc_TypeError, "Incorrect size for tuple VU");
        return NULL;
    }
    
    double **VU = (m==0) ? NULL : malloc(m*(sizeof(double*)+sizeof(int)+sizeof(npy_int)));
    double **V = VU, **U=VU+mV;
    int *DI = VU+m;
    npy_int *dims = DI+m;

    for(i=0; i<m; i++) {
        VU_arg = PyTuple_GetItem(VU_args,i);
        VU[i] = PyArray_DATA(VU_arg);
        dims[m-1-i] = DI[i] = PyArray_SIZE(VU_arg);
    }
    
    if (m == 0) {
        npy_int ones[] = {1};
        L = PyArray_FromDims(1, ones, NPY_DOUBLE);
    } else {
        L = PyArray_FromDims(m, dims, NPY_DOUBLE);
    }
    double *Y = PyArray_DATA(L);
    
    dblMaxw(f_r, &vt, v0, &ut, u0, IJ, nrm, DI, V, U, Y);

    if (VU) free(VU);

    return L;
}


static PyMethodDef methods[] = {
    {"integrate", integrate, METH_VARARGS, "Double-Maxwellian integral"},
    {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
initdblMaxw(void)
{
    Py_InitModule("dblMaxw", methods);
    import_array();
}

