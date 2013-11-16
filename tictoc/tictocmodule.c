#include <Python.h>
#include <mach/mach_time.h>
#define MAXDEPTH 100

uint64_t start[MAXDEPTH];
int lvl=0;

static PyObject* tic(PyObject *self, PyObject *args)
{
    if (lvl == MAXDEPTH) {
        PyErr_SetString(PyExc_Exception, "Maximum depth exceeded in tic()");
        return NULL;
    }
    start[lvl++] = mach_absolute_time();
    Py_RETURN_NONE;
}

static PyObject* toc(PyObject *self, PyObject *args)
{
    if (lvl == 0) {
        PyErr_SetString(PyExc_Exception, "Called toc() without tic()");
        return NULL;
    }
    return PyFloat_FromDouble(
            (double)(mach_absolute_time() - start[--lvl]) / 1000000000L);
}

static PyObject* res(PyObject *self, PyObject *args)
{
    return tic(NULL, NULL), toc(NULL, NULL);
}


static PyMethodDef methods[] = {
    {"tic", tic, METH_NOARGS, "Start timer"},
    {"toc", toc, METH_NOARGS, "Stop timer"},
    {"res", res, METH_NOARGS, "Test timer resolution"},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
inittictoc(void)
{
    //...to convince yourself that absolute time is actually in nanoseconds...
    //static mach_timebase_info_data_t sTimebaseInfo;
    //mach_timebase_info(&sTimebaseInfo);
    //printf("multiplier %u / %u\n", sTimebaseInfo.numer, sTimebaseInfo.denom);

    Py_InitModule("tictoc", methods);
}


