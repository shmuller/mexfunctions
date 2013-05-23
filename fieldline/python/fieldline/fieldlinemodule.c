#include <Python.h>
#include <numpy/arrayobject.h>

#include <fieldline.h>

#include <odepack.h>
#include <dierckx.h>

#include <stdio.h>

void get_dierckx_data(dierckx_data *spl, PyObject *arg) {
    PyObject *obj_tx, *obj_ty, *obj_c;

    PyArg_ParseTuple(arg, "OOOii", &obj_tx, &obj_ty, &obj_c, &spl->kx, &spl->ky);

    spl->tx = PyArray_DATA(obj_tx);
    spl->nx = PyArray_DIM(obj_tx, 0);

    spl->ty = PyArray_DATA(obj_ty);
    spl->ny = PyArray_DIM(obj_ty, 0);

    spl->c = PyArray_DATA(obj_c);

    spl->mx = spl->my = 1;
}


static PyObject *meth_solve_bdry(PyObject *self, PyObject *args)
{
    dierckx_data SPLR={0}, SPLz={0};
    data DATA={0};
    data *D = &DATA;

    fieldline_data FIELDLINE = {&SPLR, &SPLz, &DATA};
    fieldline_data *F = &FIELDLINE;

    get_dierckx_data(F->splR, PyTuple_GET_ITEM(args, 0));
    get_dierckx_data(F->splz, PyTuple_GET_ITEM(args, 1));

    PyObject *obj;

    obj = PyTuple_GET_ITEM(args, 2);
    D->y = PyArray_DATA(obj);
    
    obj = PyTuple_GET_ITEM(args, 3);
    D->t = PyArray_DATA(obj);
    D->n = PyArray_DIM(obj, 0);

    obj = PyTuple_GET_ITEM(args, 4);
    F->bdry = PyArray_DATA(obj);
    F->n_bdry = PyArray_DIM(obj, 0) / 2;

    solve_bdry(F);

    return Py_BuildValue("i", D->points_done);
}


static PyMethodDef methods[] = {
    {"solve_bdry", meth_solve_bdry, METH_VARARGS, "Field line solver"},
    {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
init_fieldline(void)
{
    Py_InitModule("_fieldline", methods);
}

