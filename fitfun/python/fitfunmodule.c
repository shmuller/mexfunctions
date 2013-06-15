#include <Python.h>
#include <numpy/arrayobject.h>

#include <math.h>

#include "../fitfun.h"

#define get_arr(args, i) PyArray_DATA(PyTuple_GET_ITEM(args, i))

void parse_args(PyObject *args, data *D)
{
    PyObject *obj;
    obj = PyTuple_GET_ITEM(args, 0);
    D->P = PyArray_DATA(obj);
    D->n = PyArray_DIM(obj, 0);

    obj = PyTuple_GET_ITEM(args, 1);
    D->x = PyArray_DATA(obj);
    D->m = PyArray_DIM(obj, 0);

    obj = PyTuple_GET_ITEM(args, 2);
    D->y = PyArray_DATA(obj);
}

void parse_args_a(PyObject *args, data *D)
{
    parse_args(args, D);
    D->a = get_arr(args, 3);
}

void parse_args_ydata(PyObject *args, data *D)
{
    parse_args(args, D);
    D->ydata = get_arr(args, 3);
}

void parse_args_a_ydata(PyObject *args, data *D)
{
    parse_args_ydata(args, D);
    D->a = get_arr(args, 4);
}


#define meth_template_passthru(fun, parser, i)                \
static PyObject* meth_##fun(PyObject *self, PyObject *args) { \
    data D = {0};                                             \
    parser(args, &D);                                         \
    fun(&D);                                                  \
    PyObject *obj = PyTuple_GET_ITEM(args, i);                \
    Py_INCREF(obj);                                           \
    return obj;                                               \
}

#define meth_template_double(fun, parser)                     \
static PyObject* meth_##fun(PyObject *self, PyObject *args) { \
    data D = {0};                                             \
    parser(args, &D);                                         \
    double d = fun(&D);                                       \
    return Py_BuildValue("d", d);                             \
}

#define meth_template_passthru_fit(fun, parser, i)                          \
static PyObject* meth_##fun(PyObject *self, PyObject *args, PyObject *kw) { \
    PyObject *obj;                                                          \
    data D = {0};                                                           \
    parser(args, &D);                                                       \
    if ((obj = PyDict_GetItemString(kw, "do_var")))                         \
        D.do_var = PyArray_DATA(obj);                                       \
    fun(&D);                                                                \
    obj = PyTuple_GET_ITEM(args, i);                                        \
    Py_INCREF(obj);                                                         \
    return obj;                                                             \
}


#define meth_template_all(fun, parser)                   \
meth_template_passthru(fun, parser, 2)                   \
meth_template_passthru(fun##_diff, parser##_ydata, 2)    \
meth_template_passthru_fit(fun##_fit, parser##_ydata, 0) \
meth_template_double(fun##_rms, parser)

meth_template_all(e2, parse_args)
meth_template_all(IV3, parse_args)
meth_template_all(IV4, parse_args_a)
meth_template_all(IV5, parse_args_a)
meth_template_all(IV6, parse_args_a)
meth_template_all(IVdbl, parse_args)
meth_template_all(IVdbl2, parse_args_a)


#define method_table_entries(fun, comment)                               \
{#fun, meth_##fun, METH_VARARGS, comment},                               \
{#fun"_diff", meth_##fun##_diff, METH_VARARGS, "Difference to "comment}, \
{#fun"_rms", meth_##fun##_rms, METH_VARARGS, "rms for "comment},         \
{#fun"_fit", (PyCFunction) meth_##fun##_fit, METH_VARARGS|METH_KEYWORDS, "Fit with "comment}

static PyMethodDef methods[] = {
    method_table_entries(e2, "Exp with 2 parameters"),
    method_table_entries(IV3, "IV curve with 3 parameters"),
    method_table_entries(IV4, "IV curve with 4 parameters"),
    method_table_entries(IV5, "IV curve with 5 parameters"),
    method_table_entries(IV6, "IV curve with 6 parameters"),
    method_table_entries(IVdbl, "IV curve for double probe"),
    method_table_entries(IVdbl2, "IV curve for double probe with linear variations"),
    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC
initfitfun(void)
{
    Py_InitModule("fitfun", methods);
    import_array();
}

