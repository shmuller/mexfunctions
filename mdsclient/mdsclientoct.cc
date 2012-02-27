#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <octave/oct.h>

#include "mdsclient.h"

int sm_error(const char *errstr) 
{
    fprintf(stderr, "%s\n", errstr);
    exit(0);
}

int oct2mds(Descrip *D, const octave_value &in)
{
    int i, num, siz;
    int ndims = in.ndims();
    dim_vector dv = in.dims();
    for(i=ndims-1; i>=0; i--) if (dv(i)==1) ndims--; else break;

    int *dims = (int*) malloc(ndims*sizeof(int));
    for(i=0,num=1; i<ndims; i++) num *= dims[i] = dv(i);
    siz = in.byte_size()/num;

    if (in.is_string()) {
        const charNDArray tmp = in.char_array_value();
        void *out = calloc(num+1,sizeof(char));
        memcpy(out,tmp.data(),num);
        mkDescrip(D, w_dtype_CSTRING, 0, NULL, 0, siz, out);
    } else if (in.is_real_type()) {
        if (in.is_double_type()) {
  	    const NDArray tmp = in.array_value();
            mkDescrip(D, w_dtype_DOUBLE, ndims, dims, num, siz, (void*)tmp.data());
        } else if (in.is_single_type()) {
            const FloatNDArray tmp = in.float_array_value();
            mkDescrip(D, w_dtype_FLOAT, ndims, dims, num, siz, (void*)tmp.data());
        } else if (in.is_uint8_type()) {
            const uint8NDArray tmp = in.uint8_array_value();
            mkDescrip(D, w_dtype_UCHAR, ndims, dims, num, siz, (void*)tmp.data());
        } else if (in.is_int8_type()) {
            const int8NDArray tmp = in.int8_array_value();
            mkDescrip(D, w_dtype_CHAR, ndims, dims, num, siz, (void*)tmp.data());
        } else if (in.is_uint16_type()) {
            const uint16NDArray tmp = in.uint16_array_value();
            mkDescrip(D, w_dtype_USHORT, ndims, dims, num, siz, (void*)tmp.data());
        } else if (in.is_int16_type()) {
            const int16NDArray tmp = in.int16_array_value();
            mkDescrip(D, w_dtype_SHORT, ndims, dims, num, siz, (void*)tmp.data());
        } else if (in.is_uint32_type()) {
            const uint32NDArray tmp = in.uint32_array_value();
            mkDescrip(D, w_dtype_ULONG, ndims, dims, num, siz, (void*)tmp.data());
        } else if (in.is_int32_type()) {
            const int32NDArray tmp = in.int32_array_value();
            mkDescrip(D, w_dtype_LONG, ndims, dims, num, siz, (void*)tmp.data());
        } else if (in.is_uint64_type()) {
            const uint64NDArray tmp = in.uint64_array_value();
            mkDescrip(D, w_dtype_ULONGLONG, ndims, dims, num, siz, (void*)tmp.data());
        } else if (in.is_int64_type()) {
            const int64NDArray tmp = in.int64_array_value();
            mkDescrip(D, w_dtype_LONGLONG, ndims, dims, num, siz, (void*)tmp.data());
        }
    } else if (in.is_double_type()) {
        const ComplexNDArray tmp = in.complex_array_value();
        mkDescrip(D, w_dtype_COMPLEX_DOUBLE, ndims, dims, num, siz, (void*)tmp.data());
    } else if (in.is_single_type()) {
        const FloatComplexNDArray tmp = in.float_complex_array_value();
        mkDescrip(D, w_dtype_COMPLEX, ndims, dims, num, siz, (void*)tmp.data());
    }
    return 0;
}

int mds2oct(octave_value &out, const Descrip *D)
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
    return 0;
}


DEFUN_DLD(mdsclientmex, args, nargout, "MDSplus client")
{
    int i;
    octave_value_list retval;

    char host[] = "localhost:8010";
    int sock;

    if ((sock=sm_mdsconnect(host)) < 0) {
        sm_error("Could not connect.\n");
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
