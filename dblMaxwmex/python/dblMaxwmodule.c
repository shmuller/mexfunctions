/* dblMaxw.integrate("f_r",vt,v0,ut,u0,IJ,(VU))
 * Double-Maxwellian integrals. Compile with:
 *
 * python setup.py install
 *
 * S. H. Muller, 2012/03/29
 */

#include <Python.h>
#include <numpy/arrayobject.h>

#include "../dblMaxw.h"

int dims_ok(PyObject *in, int m, npy_int *dims)
{
    int i;
    if (PyArray_NDIM(in) != m) return 0;

    npy_intp *dv = PyArray_DIMS(in);
    for(i=0; i<m; i++) if (dv[i] != dims[i]) return 0;
    return 1;
}

int size_ok(PyObject *in, int m, npy_int *dims)
{
    int i, N=1;
    for(i=0; i<m; i++) N *= dims[i];
    
    return PyArray_SIZE(in) == N;
}


static PyObject *integrate(PyObject *self, PyObject *args)
{
    int i;
    char *f_r, *nrm;
    double vt, ut, r0;
    PyObject *a[2], *fM_args=NULL, *VU_args=NULL, *VU_arg, *out=NULL;
    if (!PyArg_ParseTuple(args, "zdOOz|OO", 
            &f_r, &r0, &fM_args, a, &nrm, &VU_args, &out)) {
        PyErr_SetString(PyExc_TypeError, "f_r, r0, (fM), IJ, nrm, [(VU), out] expected");
        return NULL;
    }
    int *IJ = PyArray_DATA(a[0]);

    if (!PyArg_ParseTuple(fM_args, "dOdO", &vt, a, &ut, a+1)) {
        PyErr_SetString(PyExc_TypeError, "fM = (vt, v0, ut, u0) expected");
        return NULL;
    }
    double *v0 = PyArray_DATA(a[0]);
    double *u0 = PyArray_DATA(a[1]);
        
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
    
    if (out == NULL) {
        if (m == 0) {
            npy_int ones[] = {1};
            out = PyArray_FromDims(1, ones, NPY_DOUBLE);
        } else {
            out = PyArray_FromDims(m, dims, NPY_DOUBLE);
        }
    } else if (!size_ok(out, m, dims)) {
        PyErr_SetString(PyExc_TypeError, "Incorrect size of output array");
        if (VU) free(VU);
        return NULL;
    } else {
        Py_INCREF(out);
    }
    double *Y = PyArray_DATA(out);
    
    dblMaxw(f_r, &r0, &vt, v0, &ut, u0, IJ, nrm, DI, V, U, Y);

    if (VU) free(VU);
    return out;
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

