#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>

#include <odepack.h>

data DATA = {0};
data *D = &DATA;

PyObject *odefun=NULL, *odeargs=NULL, *odeterm=NULL;


void python_wrapper(int *neq, double *t, double *y, double *ydot)
{
    // temporarily override data sections of Python arguments
    *(double**)D->y = y;
    *(double**)D->t = t;
    *(double**)D->ydot = ydot;

    // call back to Python for function evaluation
    PyEval_CallObject(odefun, odeargs);
}


int python_term(data *D) {
    return PyObject_IsTrue(PyEval_CallObject(odeterm, odeargs));
}


static PyObject *meth_odesolve(PyObject *self, PyObject *args)
{
    int i, nargs=PyTuple_GET_SIZE(args);
    PyObject *res, *time, *t, *y, *ydot;
    void *y_save, *t_save, *ydot_save;

    D->wrapper = python_wrapper;

    // get callback function, results and time vectors
    odefun = PyTuple_GET_ITEM(args, 0);
    res  = PyTuple_GET_ITEM(args, 1);
    time = PyTuple_GET_ITEM(args, 2);

    D->res  = PyArray_DATA(res);
    D->time = PyArray_DATA(time);

    D->n = PyArray_DIM(time, 0);
    for (i=PyArray_NDIM(res), D->neq=1.; i--; ) {
        D->neq *= PyArray_DIM(res, i);
    }
    D->neq /= D->n;
    
    // get arguments to be passed to odefun
    odeargs = PyTuple_GET_ITEM(args, 3);

    if (nargs > 4) {
        odeterm = PyTuple_GET_ITEM(args, 4);
        D->term = python_term;
    }

    // parse odeargs for y, t and ydot
    y    = PyTuple_GET_ITEM(odeargs, 0);
    t    = PyTuple_GET_ITEM(odeargs, 1);
    ydot = PyTuple_GET_ITEM(odeargs, 2);

    // store the addresses of the data blocks for simple substitution
    D->y = (double*) &((PyArrayObject*)y)->data;
    D->t = (double*) &((PyArrayObject*)t)->data;
    D->ydot = (double*) &((PyArrayObject*)ydot)->data;
    
    // save the original data vectors
    y_save = *(void**)D->y;
    t_save = *(void**)D->t;
    ydot_save = *(void**)D->ydot;

    // call solver
    odesolve(D);

    // reattach original data vectors to Python arguments
    *(void**)D->y = y_save;
    *(void**)D->t = t_save;
    *(void**)D->ydot = ydot_save;

    // return points_done
    return Py_BuildValue("i", D->points_done);
}


static PyMethodDef methods[] = {
    {"odesolve", meth_odesolve, METH_VARARGS, "ODE solver"},
    {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
initodepack(void)
{
    Py_InitModule("odepack", methods);
}

