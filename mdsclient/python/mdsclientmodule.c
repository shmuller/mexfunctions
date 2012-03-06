/* Mdsplus client wrapper for Python. Compile with:
 * 
 * python setup.py install
 *
 * S. H. Muller, 2012/03/06
 */

#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <Python.h>
#include <numpy/arrayobject.h>

#include "../mdsclient.h"

#define ERROR(t,m) PyErr_SetString(t,m); return NULL

void mds2py_type(w_dtype_t w_dtype, int *typenum)
{
    switch (w_dtype)
    {
        case w_dtype_CSTRING         :  *typenum = NPY_STRING;  break;
        case w_dtype_UCHAR           :  *typenum = NPY_UINT8;   break;
        case w_dtype_CHAR            :  *typenum = NPY_INT8;    break;
        case w_dtype_USHORT          :  *typenum = NPY_UINT16;  break;
        case w_dtype_SHORT           :  *typenum = NPY_INT16;   break;
        case w_dtype_ULONG           :  *typenum = NPY_UINT32;  break;
        case w_dtype_LONG            :  *typenum = NPY_INT32;   break;
        case w_dtype_ULONGLONG       :  *typenum = NPY_UINT64;  break;
        case w_dtype_LONGLONG        :  *typenum = NPY_INT64;   break;
        case w_dtype_FLOAT           :  *typenum = NPY_FLOAT;   break;
        case w_dtype_DOUBLE          :  *typenum = NPY_DOUBLE;  break;
        case w_dtype_COMPLEX         :  *typenum = NPY_CFLOAT;  break;
        case w_dtype_COMPLEX_DOUBLE  :  *typenum = NPY_CDOUBLE; break;
        default                      :  *typenum = NPY_NOTYPE;  break;
    }
}

void py2mds_type(int typenum, w_dtype_t* w_dtype)
{
    switch (typenum)
    {
        case NPY_STRING   :  *w_dtype = w_dtype_CSTRING;        break;
        case NPY_UINT8    :  *w_dtype = w_dtype_UCHAR;          break;
        case NPY_INT8     :  *w_dtype = w_dtype_CHAR;           break;
        case NPY_UINT16   :  *w_dtype = w_dtype_USHORT;         break;
        case NPY_INT16    :  *w_dtype = w_dtype_SHORT;          break;
        case NPY_UINT32   :  *w_dtype = w_dtype_ULONG;          break;
        case NPY_INT32    :  *w_dtype = w_dtype_LONG;           break;
        case NPY_UINT64   :  *w_dtype = w_dtype_ULONGLONG;      break;
        case NPY_INT64    :  *w_dtype = w_dtype_LONGLONG;       break;
        case NPY_FLOAT    :  *w_dtype = w_dtype_FLOAT;          break;
        case NPY_DOUBLE   :  *w_dtype = w_dtype_DOUBLE;         break;
        case NPY_CFLOAT   :  *w_dtype = w_dtype_COMPLEX;        break;
        case NPY_CDOUBLE  :  *w_dtype = w_dtype_COMPLEX_DOUBLE; break;
        default           :  *w_dtype = w_dtype_UNKNOWN;        break;
    }
}

void py2mds_dims(Descrip *D, const PyObject *in)
{
    int i, num, siz;
    int ndims = PyArray_NDIM(in);
    npy_intp *dv = PyArray_DIMS(in);
    
    int *dims = (ndims==0) ? NULL : (int*) malloc(ndims*sizeof(int));
    for(i=0,num=1; i<ndims; i++) num *= dims[i] = dv[ndims-1-i];
    siz = (num==0) ? 0 : PyArray_ITEMSIZE(in);

    mkDescrip_dims(D, ndims, dims, num, siz);
}

void py2mds(Descrip *D, const PyObject *in)
{
    if (PyString_Check(in)) {
        D->siz = PyString_Size(in);
        D->ptr = PyString_AsString(in);
        mkDescrip(D, w_dtype_CSTRING, 0, NULL, 0, D->siz, D->ptr);
        return;
    }

    PyObject *arr = PyArray_FromAny(in, NULL, 0, 0, 0, NULL);
    int typenum = PyArray_TYPE(arr);
    py2mds_type(typenum, &D->w_dtype);
    py2mds_dims(D, arr);

    D->ptr = PyArray_DATA(arr);
}


void mds2py(PyObject **out, const Descrip *D)
{
    int i, tmp, typenum;
    mds2py_type(D->w_dtype, &typenum);

    npy_int *dims = (D->ndims==0) ? NULL : malloc(D->ndims*sizeof(npy_int));
    for(i=0; i<D->ndims; i++) dims[i] = D->dims[D->ndims-1-i];

    if (D->w_dtype == w_dtype_UNKNOWN) {
        *out = Py_BuildValue("");
    } else if (D->w_dtype == w_dtype_CSTRING) {
        *out = Py_BuildValue("s#", D->ptr, D->num);
    } else {
        *out = PyArray_FromDims(D->ndims, dims, typenum);
        memcpy(PyArray_DATA(*out), D->ptr, D->num*D->siz);
    }
    if (dims) free(dims);
}


static PyObject* mdsconnect(PyObject *self, PyObject *args)
{
    int sock;
    char *host=NULL;
    if (!PyArg_ParseTuple(args, "s", &host)) {
        ERROR(PyExc_TypeError, "Host string expected. E.g.: 'localhost:8000'");
    }
    switch (sock=sm_mdsconnect(host)) {
        case -1:
        case -2:
        case -3:
        case -4: ERROR(PyExc_Exception, "Could not connect to server");
        case -5: ERROR(PyExc_Exception, "Could not authenticate user");
    }
    return Py_BuildValue("i", sock);
}

static PyObject* mdsdisconnect(PyObject *self, PyObject *args)
{
    int sock;
    if (!PyArg_ParseTuple(args, "i", &sock)) {
        ERROR(PyExc_TypeError, "Socket argument expected.");
    }
    sm_mdsdisconnect(sock);
    Py_RETURN_NONE;
}

static PyObject* mdsvalue(PyObject *self, PyObject *args)
{
    int i, sock;
    Descrip l, *R;
    int nR = PyTuple_Size(args);
    
    R = (Descrip*) malloc(nR*sizeof(Descrip));
    
    for(i=0; i<nR; i++) {
        py2mds(&R[i], PyTuple_GetItem(args,i));
    }

    void *mem;
    sock = *((int*)R[0].ptr);
    sm_mdsvalue(sock, &l, nR-1, R+1, &mem);

    PyObject *retval;
    mds2py(&retval, &l);
    if (mem) free(mem);

    for(i=0; i<nR; i++) if (R[i].dims) free(R[i].dims);
    free(R);

    return retval;    
}


static PyObject* insert_arg(PyObject *args, int idx, PyObject *arg)
{
    int i, len=PyTuple_GET_SIZE(args);
    PyObject *new_args = PyTuple_New(len+1);
    for(i=0; i<idx; i++) {
        PyTuple_SET_ITEM(new_args, i, PyTuple_GET_ITEM(args,i));
    }
    PyTuple_SET_ITEM(new_args, idx, arg);
    for(; i<len; i++) {
        PyTuple_SET_ITEM(new_args, i+1, PyTuple_GET_ITEM(args,i));
    }
    return new_args;
}

static PyObject* mdsopen(PyObject *self, PyObject *args)
{
    PyObject *arg = PyString_FromString("TreeOpen($,$)");
    PyObject *new_args = insert_arg(args, 1, arg);
    return mdsvalue(self, new_args);
}

static PyObject* mdsclose(PyObject *self, PyObject *args)
{
    PyObject *arg = PyString_FromString("TreeClose()");
    PyObject *new_args = insert_arg(args, 1, arg);
    return mdsvalue(self, new_args);
}


static PyMethodDef methods[] = {
    {"mdsconnect", mdsconnect, METH_VARARGS, "Connect to mdsplus server"},
    {"mdsdisconnect", mdsdisconnect, METH_VARARGS, "Disconnect from mdsplus server"},
    {"mdsvalue", mdsvalue, METH_VARARGS, "Evaluate mdsplus expression"},
    {"mdsopen", mdsopen, METH_VARARGS, "Open mdsplus tree"},
    {"mdsclose", mdsclose, METH_VARARGS, "Close mdsplus tree"},
    {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
initmdsclient(void)
{
    Py_InitModule("mdsclient", methods);
    import_array();
}


