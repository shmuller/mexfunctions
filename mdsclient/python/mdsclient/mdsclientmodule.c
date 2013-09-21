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

#ifndef PyArray_SimpleNewFromDataOwning
PyObject *PyArray_SimpleNewFromDataOwning(
    int nd, npy_intp* dims, int typenum, void* data) 
{
    PyObject *arr = PyArray_SimpleNewFromData(nd, dims, typenum, data);

#ifdef NPY_ARRAY_OWNDATA
    PyArray_ENABLEFLAGS((PyArrayObject*)arr, NPY_ARRAY_OWNDATA);
#else
    ((PyArrayObject*)arr)->flags |= NPY_OWNDATA;
#endif
    return arr;
}
#endif

#include <mdsclient.h>

#define ERROR(t,m) PyErr_SetString(t,m); return NULL

void mds2py_type(w_dtype_t w_dtype, int *typenum)
{
    switch (w_dtype)
    {
        case w_dtype_CSTRING         :  *typenum = NPY_STRING;  break;
        case w_dtype_UCHAR           :  *typenum = NPY_UBYTE;   break;
        case w_dtype_CHAR            :  *typenum = NPY_BYTE;    break;
        case w_dtype_USHORT          :  *typenum = NPY_USHORT;  break;
        case w_dtype_SHORT           :  *typenum = NPY_SHORT;   break;
        case w_dtype_ULONG           :  *typenum = NPY_UINT;    break;
        case w_dtype_LONG            :  *typenum = NPY_INT;     break;
        case w_dtype_ULONGLONG       :  *typenum = NPY_ULONG;   break;
        case w_dtype_LONGLONG        :  *typenum = NPY_LONG;    break;
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
        case NPY_UBYTE    :  *w_dtype = w_dtype_UCHAR;          break;
        case NPY_BYTE     :  *w_dtype = w_dtype_CHAR;           break;
        case NPY_USHORT   :  *w_dtype = w_dtype_USHORT;         break;
        case NPY_SHORT    :  *w_dtype = w_dtype_SHORT;          break;
        case NPY_UINT     :  *w_dtype = w_dtype_ULONG;          break;
        case NPY_INT      :  *w_dtype = w_dtype_LONG;           break;
        case NPY_ULONG    :  *w_dtype = w_dtype_ULONGLONG;      break;
        case NPY_LONG     :  *w_dtype = w_dtype_LONGLONG;       break;
        case NPY_FLOAT    :  *w_dtype = w_dtype_FLOAT;          break;
        case NPY_DOUBLE   :  *w_dtype = w_dtype_DOUBLE;         break;
        case NPY_CFLOAT   :  *w_dtype = w_dtype_COMPLEX;        break;
        case NPY_CDOUBLE  :  *w_dtype = w_dtype_COMPLEX_DOUBLE; break;
        default           :  *w_dtype = w_dtype_UNKNOWN;        break;
    }
}

void py2mds_dims(Descrip *D, PyObject *in)
{
    int i, num, siz;
    int ndims = PyArray_NDIM(in);
    
    int *dims = (ndims==0) ? NULL : (int*) malloc(ndims*sizeof(int));
    for(i=0,num=1; i<ndims; i++) num *= dims[i] = PyArray_DIM(in, ndims-1-i);
    siz = (num==0) ? 0 : PyArray_ITEMSIZE(in);

    mkDescrip_dims(D, ndims, dims, num, siz);
}

void py2mds(Descrip *D, PyObject *in)
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


void mds2py(PyObject **out, Descrip *D)
{
    int i, typenum;
    mds2py_type(D->w_dtype, &typenum);

    npy_intp *dims = (D->ndims==0) ? NULL : malloc(D->ndims*sizeof(npy_intp));
    for(i=0; i<D->ndims; i++) dims[i] = D->dims[D->ndims-1-i];

    if (D->w_dtype == w_dtype_UNKNOWN) {
        *out = Py_BuildValue("");
    } else if (D->w_dtype == w_dtype_CSTRING) {
        *out = Py_BuildValue("s#", D->ptr, D->num);
    } else {
        /* Mdsplus allocates memory in GetMdsMsg() in mdstcpip/GetMdsMsg.c. However,
         * D->ptr doesn't point to the beginning of the allocated memory, so Python 
         * would be unable to free it. Hence we must make a copy of the data section.
         */
        //*out = PyArray_SimpleNew(D->ndims, dims, typenum);
        //memcpy(PyArray_DATA(*out), D->ptr, D->num*D->siz);

        /* SHM 2013/05/27: The sources relevant for Mdsplus client functionality were
         * extracted into mdsiputil.c and modified in such a way that a pointer to the 
         * beginning of the data section is returned. The memory can thus be used
         * directly by Numpy.
         */
        *out = PyArray_SimpleNewFromDataOwning(D->ndims, dims, typenum, D->ptr);
        D->ptr = NULL; // flag to avoid freeing
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

    sock = *((int*)R[0].ptr);
    sm_mdsvalue(sock, &l, nR-1, R+1);

    PyObject *retval;
    mds2py(&retval, &l);
    if (l.ptr) free(l.ptr);

    for(i=0; i<nR; i++) if (R[i].dims) free(R[i].dims);
    free(R);

    return retval;    
}


#define USAGE_MDSCONNECT    "sock = mdsconnect('hostname:port')"
#define USAGE_MDSDISCONNECT "mdsdisconnect(sock)"
#define USAGE_MDSVALUE      "x = mdsvalue(sock, 'mdsexpr', arg1, arg2, ...)"

static PyMethodDef methods[] = {
    {"mdsconnect", mdsconnect, METH_VARARGS, USAGE_MDSCONNECT},
    {"mdsdisconnect", mdsdisconnect, METH_VARARGS, USAGE_MDSDISCONNECT},
    {"mdsvalue", mdsvalue, METH_VARARGS, USAGE_MDSVALUE},
    {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
init_mdsclient(void)
{
    Py_InitModule("_mdsclient", methods);
    import_array();
}


