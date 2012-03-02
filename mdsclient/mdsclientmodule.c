#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <Python.h>
#include <numpy/noprefix.h>
//#include <numpy/arrayobject.h>

#include "mdsclient.h"

#define ERROR(x) fprintf(stderr,"%s\n",x); return NULL

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


/*
void *octGetData(const octave_value &in)
{
    if (in.is_complex_type() && !in.is_scalar_type()) {
        // handle complex data types separately, but only if not scalar!
        if (in.is_double_type()) {
            const ComplexNDArray t = in.complex_array_value();
            return (void*) t.data();
	} else if (in.is_single_type()) {
            const FloatComplexNDArray t = in.float_complex_array_value();
            return (void*) t.data();
        } else {
            error("Data type not implemented.");
            return NULL;
	}
    } else {
        // handle bulk of data types with mex_get_data()
        return in.mex_get_data();
    }
}


void oct2mds_dtype(Descrip *D, const octave_value &in)
{
    if (in.is_string()) {
        D->w_dtype = w_dtype_CSTRING;
    } else if (in.is_real_type()) {
        if (in.is_double_type()) {
            D->w_dtype = w_dtype_DOUBLE;
        } else if (in.is_single_type()) {
            D->w_dtype = w_dtype_FLOAT;
        } else if (in.is_uint8_type()) {
            D->w_dtype = w_dtype_UCHAR;
        } else if (in.is_int8_type()) {
            D->w_dtype = w_dtype_CHAR;
        } else if (in.is_uint16_type()) {
            D->w_dtype = w_dtype_USHORT;
        } else if (in.is_int16_type()) {
            D->w_dtype = w_dtype_SHORT;
        } else if (in.is_uint32_type()) {
            D->w_dtype = w_dtype_ULONG;
        } else if (in.is_int32_type()) {
            D->w_dtype = w_dtype_LONG;
        } else if (in.is_uint64_type()) {
            D->w_dtype = w_dtype_ULONGLONG;
        } else if (in.is_int64_type()) {
            D->w_dtype = w_dtype_LONGLONG;
        }
    } else if (in.is_double_type()) {
        D->w_dtype = w_dtype_COMPLEX_DOUBLE;
    } else if (in.is_single_type()) {
        D->w_dtype = w_dtype_COMPLEX;
    } else {
        error("Unknown data type.");
    }
}

void oct2mds_dims(Descrip *D, const octave_value &in)
{
    int i, num, siz;
    int ndims = in.ndims();
    dim_vector dv = in.dims();

    // remove singleton dimensions
    for(i=ndims-1; i>=0; i--) if (dv(i)==1) ndims--; else break;

    int *dims = (ndims==0) ? NULL : (int*) malloc(ndims*sizeof(int));
    for(i=0,num=1; i<ndims; i++) num *= dims[i] = dv(i);
    siz = (num==0) ? 0 : in.byte_size()/num;

    mkDescrip_dims(D, ndims, dims, num, siz);
}

void oct2mds(Descrip *D, const octave_value &in)
{
    oct2mds_dims(D, in);
    oct2mds_dtype(D, in);

    D->ptr = octGetData(in);

    if (in.is_string()) {
        void *ptr = calloc(D->num+1,sizeof(char));
        memcpy(ptr,D->ptr,D->num);
	mkDescrip(D, D->w_dtype, 0, NULL, 0, D->siz, ptr);
    }
}


void mds2oct(octave_value &out, const Descrip *D)
{
    if (D->w_dtype == w_dtype_UNKNOWN) {
        out = NDArray();
	return;
    }
    int i, numbytes = D->num*D->siz;
    int ndims = (D->ndims > 2) ? D->ndims : 2;

    dim_vector dv;
    dv.resize(ndims);
    for(i=0; i<D->ndims; i++) dv(i) = D->dims[i];
    for(; i<ndims; i++) dv(i) = 1;

    if (D->w_dtype == w_dtype_CSTRING) dv(1) = D->num;

    switch (D->w_dtype) {
        case w_dtype_CSTRING:        out = charNDArray(dv);         break;
        case w_dtype_UCHAR:          out = uint8NDArray(dv);        break;
        case w_dtype_CHAR:           out = int8NDArray(dv);         break;
        case w_dtype_USHORT:         out = uint16NDArray(dv);       break;
        case w_dtype_SHORT:          out = int16NDArray(dv);        break;
        case w_dtype_ULONG:          out = uint32NDArray(dv);       break;
        case w_dtype_LONG:           out = int32NDArray(dv);        break;
        case w_dtype_ULONGLONG:      out = uint64NDArray(dv);       break;
        case w_dtype_LONGLONG:       out = int64NDArray(dv);        break;
        case w_dtype_FLOAT:          out = FloatNDArray(dv);        break;
        case w_dtype_COMPLEX:        out = FloatComplexNDArray(dv); break;
        case w_dtype_DOUBLE:         out = NDArray(dv);             break;
        case w_dtype_COMPLEX_DOUBLE: out = ComplexNDArray(dv);      break;
    }
    memcpy(octGetData(out), D->ptr, numbytes);
}



DEFUN_DLD(mdsclientmex, args, nargout, "MDSplus client")
{
    int i, sock;
    octave_value_list retval;

    Descrip l, *R;    
    int nR = args.length();

    R = (Descrip*) malloc(nR*sizeof(Descrip));
    
    for(i=0; i<nR; i++) {
        oct2mds(&R[i], args(i));
    }

    char *cmd = (char*) R[0].ptr;

    if (strcmp(cmd,"mdsconnect")==0) 
    {
        char *host = (char*) R[1].ptr;
        switch (sock=sm_mdsconnect(host)) {
            case -1:
            case -2:
            case -3:
            case -4: ERROR("Could not connect to server");
            case -5: ERROR("Could not authenticate user");
        }
	dim_vector dv(1);
	int32NDArray t(dv,sock);
	retval(0) = t;
    } 
    else if (strcmp(cmd,"mdsvalue")==0)
    {
        void *mem;
        sock = *((int*)R[1].ptr);
        sm_mdsvalue(sock, &l, nR-2, R+2, &mem);

        mds2oct(retval(0), &l);
        if (mem) free(mem);
    }
    else if (strcmp(cmd,"mdsdisconnect")==0)
    {
        sock = *((int*)R[1].ptr);
        sm_mdsdisconnect(sock);
    }
 
    for(i=0; i<nR; i++) if (R[i].dims) free(R[i].dims);
    free(R);
    
    return retval;
}

*/


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
    free(dims);
}


static PyObject* mdsconnect(PyObject *self, PyObject *args)
{
    int sock;
    char *host=NULL;
    if (!PyArg_ParseTuple(args, "z", &host)) {
        return NULL;
    }
    switch (sock=sm_mdsconnect(host)) {
        case -1:
        case -2:
        case -3:
        case -4: ERROR("Could not connect to server");
        case -5: ERROR("Could not authenticate user");
    }
    return Py_BuildValue("i", sock);

    Py_RETURN_NONE;
}

static PyObject* mdsdisconnect(PyObject *self, PyObject *args)
{
    int sock;
    if (!PyArg_ParseTuple(args, "i", &sock)) {
        return NULL;
    }
    sm_mdsdisconnect(sock);
    Py_RETURN_NONE;
}


static PyObject* mdsvalue(PyObject *self, PyObject *args)
{
    int sock;
    char *cmd=NULL;
    if (!PyArg_ParseTuple(args, "iz", &sock, &cmd)) {
        return NULL;
    }
    Descrip l, R;
    mkDescrip(&R, w_dtype_CSTRING, 0, NULL, 0, strlen(cmd), cmd);

    void *mem;
    sm_mdsvalue(sock, &l, 1, &R, &mem);

    PyObject *retval;
    mds2py(&retval, &l);
    if (mem) free(mem);

    return retval;
}


static PyObject* mdsopen(PyObject *self, PyObject *args)
{
    Py_RETURN_NONE;
}

static PyObject* mdsclose(PyObject *self, PyObject *args)
{
    Py_RETURN_NONE;
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


