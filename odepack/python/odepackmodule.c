#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>

#include <odepack.h>

data DATA;
data *D = &DATA;

PyObject *odefun=NULL, *odeargs=NULL, *oderoot=NULL;
void **p_y=NULL, **p_t=NULL, **p_ydot=NULL, **p_gout=NULL;


void python_f(int *neq, double *t, double *y, double *ydot) {
    // temporarily override data sections of Python arguments
    *p_y = y;
    *p_t = t;
    *p_ydot = ydot;

    // call back to Python for function evaluation
    PyEval_CallObject(odefun, odeargs);
}

void python_g(int *neq, double *t, double *y, int *ng, double *gout) {
    // temporarily override data sections of Python arguments
    *p_y = y;
    *p_t = t;
    *p_gout = gout;

    // call back to Python for constraint evaluation
    PyEval_CallObject(oderoot, odeargs);
}


static PyObject *meth_odesolve(PyObject *self, PyObject *args)
{
    int i, nargs=PyTuple_GET_SIZE(args);
    PyObject *res, *time, *y, *t, *ydot, *gout;
    void *y_save, *t_save, *ydot_save, *gout_save;

    // reset all structure fields to 0 between calls
    DATA = empty_data;

    D->f = python_f;

    // get callback function, results and time vectors
    odefun = PyTuple_GET_ITEM(args, 0);
    res  = PyTuple_GET_ITEM(args, 1);
    time = PyTuple_GET_ITEM(args, 2);

    D->y = PyArray_DATA(res);
    D->t = PyArray_DATA(time);

    D->n = PyArray_DIM(time, 0);
    for (i=PyArray_NDIM(res), D->neq=1.; i--; ) {
        D->neq *= PyArray_DIM(res, i);
    }
    D->neq /= D->n;
    
    // get arguments to be passed to odefun
    odeargs = PyTuple_GET_ITEM(args, 3);

    // parse odeargs for y, t and ydot
    y    = PyTuple_GET_ITEM(odeargs, 0);
    t    = PyTuple_GET_ITEM(odeargs, 1);
    ydot = PyTuple_GET_ITEM(odeargs, 2);

    // store the addresses of the data blocks for simple substitution
    p_y = (void**) &((PyArrayObject*)y)->data;
    p_t = (void**) &((PyArrayObject*)t)->data;
    p_ydot = (void**) &((PyArrayObject*)ydot)->data;
    
    // save the original data vectors
    y_save = *p_y;
    t_save = *p_t;
    ydot_save = *p_ydot;

    if (nargs < 5) {
        // call solver without termination conditions
        odesolve(D);
    } else {
        D->g = python_g;
        oderoot = PyTuple_GET_ITEM(args, 4);
        gout = PyTuple_GET_ITEM(odeargs, 3);
        D->ng = PyArray_DIM(gout, 0);
        p_gout = (void**) &((PyArrayObject*)gout)->data;
        gout_save = *p_gout;

        // call solver with termination conditions
        odesolve(D);

        *p_gout = gout_save;
    }
        
    // reattach original data vectors to Python arguments
    *p_y = y_save;
    *p_t = t_save;
    *p_ydot = ydot_save;

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

