#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <octave/oct.h>

#include "mdsclient.h"

#define ERROR(x) error(x); return octave_value_list()

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
    //for(i=ndims-1; i>=0; i--) if (dv(i)==1) ndims--; else break;

    int *dims = (int*) malloc(ndims*sizeof(int));
    for(i=0,num=1; i<ndims; i++) num *= dims[i] = dv(i);
    siz = in.byte_size()/num;

    mkDescrip_dims(D, ndims, dims, num, siz);
}

void oct2mds(Descrip *D, const octave_value &in)
{
    oct2mds_dims(D, in);
    oct2mds_dtype(D, in);

    if (in.is_string()) 
    {
        void *ptr = calloc(D->num+1,sizeof(char));
        memcpy(ptr,in.mex_get_data(),D->num);
	mkDescrip(D, D->w_dtype, 0, NULL, 0, D->siz, ptr);
    } 
    else if (in.is_real_type() || in.is_scalar_type()) 
    {
        // real or scalar complex
        D->ptr = in.mex_get_data();
    } 
    else if (in.is_double_type()) 
    {
        const ComplexNDArray t = in.complex_array_value();
        D->ptr = (void*)t.data();
    } 
    else if (in.is_single_type()) 
    {
        const FloatComplexNDArray t = in.float_complex_array_value();
        D->ptr = (void*)t.data();
    }
}


void mds2oct(octave_value &out, const Descrip *D)
{
    int i, numbytes = D->num*D->siz;
    int ndims = (D->ndims > 2) ? D->ndims : 2;

    dim_vector dv;
    dv.resize(ndims);
    for(i=0; i<D->ndims; i++) dv(i) = D->dims[i];
    for(; i<ndims; i++) dv(i) = 1;

    switch (D->w_dtype) {
       case w_dtype_CSTRING: {
           dv(1) = D->num;
           charNDArray t(dv);
	   memcpy((void*)t.data(), D->ptr, numbytes);
           out = t; break;
       }
       case w_dtype_UCHAR: {
           uint8NDArray t(dv);
           memcpy((void*)t.data(), D->ptr, numbytes);
           out = t; break;
       }
       case w_dtype_CHAR: {
           int8NDArray t(dv);
           memcpy((void*)t.data(), D->ptr, numbytes);
           out = t; break;
       }
       case w_dtype_USHORT: {
           uint16NDArray t(dv);
           memcpy((void*)t.data(), D->ptr, numbytes);
           out = t; break;
       }
       case w_dtype_SHORT: {
           int16NDArray t(dv);
           memcpy((void*)t.data(), D->ptr, numbytes);
           out = t; break;
       }
       case w_dtype_ULONG: {
           uint32NDArray t(dv);
           memcpy((void*)t.data(), D->ptr, numbytes);
           out = t; break;
       }
       case w_dtype_LONG: {
           int32NDArray t(dv);
           memcpy((void*)t.data(), D->ptr, numbytes);
           out = t; break;
       }
       case w_dtype_ULONGLONG: {
           uint64NDArray t(dv);
           memcpy((void*)t.data(), D->ptr, numbytes);
           out = t; break;
       }
       case w_dtype_LONGLONG: {
           int64NDArray t(dv);
           memcpy((void*)t.data(), D->ptr, numbytes);
           out = t; break;
       }
       case w_dtype_FLOAT: {
           FloatNDArray t(dv);
           memcpy((void*)t.data(), D->ptr, numbytes);
           out = t; break;
       }
       case w_dtype_COMPLEX: {
           FloatComplexNDArray t(dv);
           memcpy((void*)t.data(), D->ptr, numbytes);
           out = t; break;
       }
       case w_dtype_DOUBLE: {
           NDArray t(dv);
           memcpy((void*)t.data(), D->ptr, numbytes);
           out = t; break;
       }
       case w_dtype_COMPLEX_DOUBLE: {
           ComplexNDArray t(dv);
           memcpy((void*)t.data(), D->ptr, numbytes);
           out = t; break;
       }
    }
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
        if ((sock=sm_mdsconnect(host)) < 0) {
            ERROR("Could not connect.");
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
 
    for(i=0; i<nR; i++) free(R[i].dims);
    free(R);
    
    return retval;
}

