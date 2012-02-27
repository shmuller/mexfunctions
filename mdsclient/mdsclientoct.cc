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

    if (in.is_string()) {
        void *ptr = calloc(D->num+1,sizeof(char));
        memcpy(ptr,in.mex_get_data(),D->num);
	mkDescrip(D, D->w_dtype, 0, NULL, 0, D->siz, ptr);
    } else if (in.is_real_type()) {
        D->ptr = in.mex_get_data();
    } else if (in.is_double_type()) {
        const ComplexNDArray tmp = in.complex_array_value();
        D->ptr = (void*)tmp.data();
    } else if (in.is_single_type()) {
        const FloatComplexNDArray tmp = in.float_complex_array_value();
        D->ptr = (void*)tmp.data();
    }
    
    /*
    if (in.is_string()) {
        const charNDArray tmp = in.char_array_value();
        void *out = calloc(D->num+1,sizeof(char));
        memcpy(out,tmp.data(),D->num);
        mkDescrip(D, w_dtype_CSTRING, 0, NULL, 0, D->siz, out);
    } else if (in.is_real_type()) {
        if (in.is_double_type()) {
  	    const NDArray tmp = in.array_value();
            mkDescrip_data(D, w_dtype_DOUBLE, (void*)tmp.data());
        } else if (in.is_single_type()) {
            const FloatNDArray tmp = in.float_array_value();
            mkDescrip_data(D, w_dtype_FLOAT, (void*)tmp.data());
        } else if (in.is_uint8_type()) {
            const uint8NDArray tmp = in.uint8_array_value();
            mkDescrip_data(D, w_dtype_UCHAR, (void*)tmp.data());
        } else if (in.is_int8_type()) {
            const int8NDArray tmp = in.int8_array_value();
            mkDescrip_data(D, w_dtype_CHAR, (void*)tmp.data());
        } else if (in.is_uint16_type()) {
            const uint16NDArray tmp = in.uint16_array_value();
            mkDescrip_data(D, w_dtype_USHORT, (void*)tmp.data());
        } else if (in.is_int16_type()) {
            const int16NDArray tmp = in.int16_array_value();
            mkDescrip_data(D, w_dtype_SHORT, (void*)tmp.data());
        } else if (in.is_uint32_type()) {
            const uint32NDArray tmp = in.uint32_array_value();
            mkDescrip_data(D, w_dtype_ULONG, (void*)tmp.data());
        } else if (in.is_int32_type()) {
            const int32NDArray tmp = in.int32_array_value();
            mkDescrip_data(D, w_dtype_LONG, (void*)tmp.data());
        } else if (in.is_uint64_type()) {
            const uint64NDArray tmp = in.uint64_array_value();
            mkDescrip_data(D, w_dtype_ULONGLONG, (void*)tmp.data());
        } else if (in.is_int64_type()) {
            const int64NDArray tmp = in.int64_array_value();
            mkDescrip_data(D, w_dtype_LONGLONG, (void*)tmp.data());
        }
    } else if (in.is_double_type()) {
        const ComplexNDArray tmp = in.complex_array_value();
        mkDescrip_data(D, w_dtype_COMPLEX_DOUBLE, (void*)tmp.data());
    } else if (in.is_single_type()) {
        const FloatComplexNDArray tmp = in.float_complex_array_value();
        mkDescrip_data(D, w_dtype_COMPLEX, (void*)tmp.data());
    }
    */
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
           charNDArray tmp(dv);
	   memcpy((void*)tmp.data(), D->ptr, numbytes);
           out = tmp; break;
       }
       case w_dtype_UCHAR: {
           uint8NDArray tmp(dv);
           memcpy((void*)tmp.data(), D->ptr, numbytes);
           out = tmp; break;
       }
       case w_dtype_CHAR: {
           int8NDArray tmp(dv);
           memcpy((void*)tmp.data(), D->ptr, numbytes);
           out = tmp; break;
       }
       case w_dtype_USHORT: {
           uint16NDArray tmp(dv);
           memcpy((void*)tmp.data(), D->ptr, numbytes);
           out = tmp; break;
       }
       case w_dtype_SHORT: {
           int16NDArray tmp(dv);
           memcpy((void*)tmp.data(), D->ptr, numbytes);
           out = tmp; break;
       }
       case w_dtype_ULONG: {
           uint32NDArray tmp(dv);
           memcpy((void*)tmp.data(), D->ptr, numbytes);
           out = tmp; break;
       }
       case w_dtype_LONG: {
           int32NDArray tmp(dv);
           memcpy((void*)tmp.data(), D->ptr, numbytes);
           out = tmp; break;
       }
       case w_dtype_ULONGLONG: {
           uint64NDArray tmp(dv);
           memcpy((void*)tmp.data(), D->ptr, numbytes);
           out = tmp; break;
       }
       case w_dtype_LONGLONG: {
           int64NDArray tmp(dv);
           memcpy((void*)tmp.data(), D->ptr, numbytes);
           out = tmp; break;
       }
       case w_dtype_FLOAT: {
           FloatNDArray tmp(dv);
           memcpy((void*)tmp.data(), D->ptr, numbytes);
           out = tmp; break;
       }
       case w_dtype_COMPLEX: {
           FloatComplexNDArray tmp(dv);
           memcpy((void*)tmp.data(), D->ptr, numbytes);
           out = tmp; break;
       }
       case w_dtype_DOUBLE: {
           NDArray tmp(dv);
           memcpy((void*)tmp.data(), D->ptr, numbytes);
           out = tmp; break;
       }
       case w_dtype_COMPLEX_DOUBLE: {
           ComplexNDArray tmp(dv);
           memcpy((void*)tmp.data(), D->ptr, numbytes);
           out = tmp; break;
       }
    }
}


DEFUN_DLD(mdsclientmex, args, nargout, "MDSplus client")
{
    int i;
    octave_value_list retval;

    char host[] = "localhost:8001";
    int sock;

    if ((sock=sm_mdsconnect(host)) < 0) {
        ERROR("Could not connect.");
    }

    printf("sock = %d\n", sock);

    //sm_mdsopen(sock, "rcp", 132777);

    Descrip l, *R;
    
    int nR = args.length();

    R = (Descrip*) malloc(nR*sizeof(Descrip));
    
    for(i=0; i<nR; i++) {
        oct2mds(&R[i], args(i));
    }

    void *mem;
    sm_mdsvalue(sock, &l, nR, R, &mem);

    mds2oct(retval(0), &l);

    if (mem) free(mem);
 
    for(i=0; i<nR; i++) free(R[i].dims);
    free(R);

    //sm_mdsclose(sock);

    sm_mdsdisconnect(sock);

    return retval;
}

