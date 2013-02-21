#include <Python.h>
#include <numpy/arrayobject.h>

#include <minpack.h>

typedef struct {
    PyObject *func;
    PyObject *args;
    double **P;
    double **y;
} container2;

container2 CONTAINER2;


void fun(data *D)
{
    container2 *C = &CONTAINER2;

    // temporarily override data sections of Python arguments
    *C->P = D->P;
    *C->y = D->y;

    // call back to Python for function evaluation
    PyEval_CallObject(C->func, C->args);
}


static PyObject *meth_leastsq(PyObject *self, PyObject *args)
{
    data DATA = {0};
    data *D = &DATA;
    container2 *C = &CONTAINER2;
    double *P_save, *y_save;

    // expect 2 arguments, the function and its arguments (p0, x, y, ...)
    C->func = PyTuple_GET_ITEM(args, 0);
    C->args = PyTuple_GET_ITEM(args, 1);

    // from the arguments, we need to understand more about p0 and y
    PyObject *P = PyTuple_GET_ITEM(C->args, 0);
    PyObject *y = PyTuple_GET_ITEM(C->args, 2);

    // store the address of the data block of p0 and y
    C->P = (double**) &((PyArrayObject*)P)->data;
    C->y = (double**) &((PyArrayObject*)y)->data;

    // fill D for leastsq, also save address of original data blocks
    D->P = P_save = PyArray_DATA(P);
    D->y = y_save = PyArray_DATA(y);
    D->n = PyArray_DIM(P, 0);
    D->m = PyArray_DIM(y, 0);

    // call fitting using wrapper function fun
    leastsq(fun, D);

    // reattach original data blocks to Python arguments
    *C->P = P_save;
    *C->y = y_save;

    // return reference to overwritten p0
    Py_INCREF(P);
    return P;
}


static PyMethodDef methods[] = {
    {"leastsq", meth_leastsq, METH_VARARGS, "Least-square fit"},
    {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
initminpack(void)
{
    Py_InitModule("minpack", methods);
}

