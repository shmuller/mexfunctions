#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include "mex.h"

#include "mdsclient.h"

#define ERROR(x) mexErrMsgTxt(x)

void mex2mds_type(mxClassID mxID, mxComplexity mxCo, w_dtype_t *w_dtype)
{
    switch (mxID)
    {
        case mxCHAR_CLASS    :  *w_dtype = w_dtype_CSTRING;    break;
        case mxUINT8_CLASS   :  *w_dtype = w_dtype_UCHAR;      break;
        case mxINT8_CLASS    :  *w_dtype = w_dtype_CHAR;       break;
        case mxUINT16_CLASS  :  *w_dtype = w_dtype_USHORT;     break;
        case mxINT16_CLASS   :  *w_dtype = w_dtype_SHORT;      break;
        case mxUINT32_CLASS  :  *w_dtype = w_dtype_ULONG;      break;
        case mxINT32_CLASS   :  *w_dtype = w_dtype_LONG;       break;
        case mxUINT64_CLASS  :  *w_dtype = w_dtype_ULONGLONG;  break;
        case mxINT64_CLASS   :  *w_dtype = w_dtype_LONGLONG;   break;
        case mxSINGLE_CLASS  :  *w_dtype = (mxCo==mxREAL) ? w_dtype_FLOAT  : w_dtype_COMPLEX;        break;
        case mxDOUBLE_CLASS  :  *w_dtype = (mxCo==mxREAL) ? w_dtype_DOUBLE : w_dtype_COMPLEX_DOUBLE; break;
        default              :  *w_dtype = w_dtype_UNKNOWN;    break;
    }
}

void mds2mex_type(w_dtype_t w_dtype, mxClassID *mxID, mxComplexity *mxCo)
{
    switch (w_dtype)
    {
        case w_dtype_CSTRING         :  *mxID = mxCHAR_CLASS;    *mxCo = mxREAL;     break;
        case w_dtype_UCHAR           :  *mxID = mxUINT8_CLASS;   *mxCo = mxREAL;     break;
        case w_dtype_CHAR            :  *mxID = mxINT8_CLASS;    *mxCo = mxREAL;     break;
        case w_dtype_USHORT          :  *mxID = mxUINT16_CLASS;  *mxCo = mxREAL;     break;
        case w_dtype_SHORT           :  *mxID = mxINT16_CLASS;   *mxCo = mxREAL;     break;
        case w_dtype_ULONG           :  *mxID = mxUINT32_CLASS;  *mxCo = mxREAL;     break;
        case w_dtype_LONG            :  *mxID = mxINT32_CLASS;   *mxCo = mxREAL;     break;
        case w_dtype_ULONGLONG       :  *mxID = mxUINT64_CLASS;  *mxCo = mxREAL;     break;
        case w_dtype_LONGLONG        :  *mxID = mxINT64_CLASS;   *mxCo = mxREAL;     break;
        case w_dtype_FLOAT           :  *mxID = mxSINGLE_CLASS;  *mxCo = mxREAL;     break;
        case w_dtype_DOUBLE          :  *mxID = mxDOUBLE_CLASS;  *mxCo = mxREAL;     break;
        case w_dtype_COMPLEX         :  *mxID = mxSINGLE_CLASS;  *mxCo = mxCOMPLEX;  break;
        case w_dtype_COMPLEX_DOUBLE  :  *mxID = mxDOUBLE_CLASS;  *mxCo = mxCOMPLEX;  break;
        default                      :  *mxID = mxUNKNOWN_CLASS; *mxCo = mxREAL;     break;
    }
}

void mex2mds_cmplx(void **buf, void *pr, void *pi, int num, int siz)
{
    int i, s = siz/2;
    void *b;
    *buf = malloc(num*siz);

	for(i=0,b=*buf; i<num; i++,b+=siz,pr+=s,pi+=s) {
		memcpy(b  ,pr,s);
		memcpy(b+s,pi,s);
	}
}

void mds2mex_cmplx(void *buf, void *pr, void *pi, int num, int siz) 
{
	int i, s = siz/2;
	void *b;

	for(i=0,b=buf; i<num; i++,b+=siz,pr+=s,pi+=s) {
		memcpy(pr,b  ,s);
		memcpy(pi,b+s,s);
	}
}


void mex2mds_dtype(Descrip *D, const mxArray *in)
{
    mxClassID mxID = mxGetClassID(in); 
    mxComplexity mxCo = (mxIsComplex(in)) ? mxCOMPLEX : mxREAL;

    mex2mds_type(mxID, mxCo, &D->w_dtype);
}

void mex2mds_dims(Descrip *D, const mxArray *in)
{
    int i, num, siz;
    int ndims = mxGetNumberOfDimensions(in);
    const mwSize *dv = mxGetDimensions(in);

    /* remove singleton dimensions */
    for(i=ndims-1; i>=0; i--) if (dv[i]==1) ndims--; else break;

    int *dims = (ndims==0) ? NULL : (int*) malloc(ndims*sizeof(int));
    for(i=0,num=1; i<ndims; i++) num *= dims[i] = dv[i];
    siz = (num==0) ? 0 : mxGetElementSize(in);
    if (mxIsComplex(in)) siz *= 2;

    mkDescrip_dims(D, ndims, dims, num, siz);
}

void mex2mds(Descrip *D, const mxArray *in)
{
    mex2mds_dims(D, in);
    mex2mds_dtype(D, in);
        
    if (mxIsChar(in)) {
        D->ptr = mxArrayToString(in);
        mkDescrip_dims(D, 0, NULL, 0, D->siz);
    } else if (!mxIsComplex(in)) {
        D->ptr = mxGetData(in);
    } else {
        mex2mds_cmplx(&D->ptr, mxGetData(in), mxGetImagData(in), D->num, D->siz);
    }
}


void mds2mex(mxArray **out, const Descrip *D)
{
    if (D->w_dtype == w_dtype_UNKNOWN) {
        *out = mxCreateNumericMatrix(0,0,mxDOUBLE_CLASS,mxREAL);
        return;
    } else if (D->w_dtype == w_dtype_CSTRING) {
        void *ptr = calloc(D->num+1,sizeof(char));
        memcpy(ptr, D->ptr, D->num*sizeof(char));
        *out = mxCreateString(ptr);
        free(ptr);
        return;
    }
    int i, numbytes = D->num*D->siz;
    int ndims = (D->ndims > 2) ? D->ndims : 2;

    mwSize *dv = (mwSize*) malloc(ndims*sizeof(mwSize));
    for(i=0; i<D->ndims; i++) dv[i] = D->dims[i];
    for(; i<ndims; i++) dv[i] = 1;

    if (D->w_dtype == w_dtype_CSTRING) dv[1] = D->num;
    
    mxClassID mxID;
    mxComplexity mxCo;
    mds2mex_type(D->w_dtype, &mxID, &mxCo);
    *out = mxCreateNumericArray(ndims, dv, mxID, mxCo);
    free(dv);
    
    if (mxCo == mxREAL) {
        memcpy(mxGetData(*out), D->ptr, numbytes);
    } else {
        mds2mex_cmplx(D->ptr, mxGetData(*out), mxGetImagData(*out), D->num, D->siz);
    }
}


void mexFunction(int nL, mxArray *retval[], int nR, const mxArray *args[])
{
    int i, sock;

    Descrip l, *R;

    R = (Descrip*) malloc(nR*sizeof(Descrip));
    
    for(i=0; i<nR; i++) {
        mex2mds(&R[i], args[i]);
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
        retval[0] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
        *((int*)mxGetData(retval[0])) = sock;
    } 
    else if (strcmp(cmd,"mdsvalue")==0)
    {
        void *mem;
        sock = *((int*)R[1].ptr);
        sm_mdsvalue(sock, &l, nR-2, R+2, &mem);

        mds2mex(retval, &l);
        if (mem) free(mem);
    }
    else if (strcmp(cmd,"mdsdisconnect")==0)
    {
        sock = *((int*)R[1].ptr);
        sm_mdsdisconnect(sock);
    }
 
    for(i=0; i<nR; i++) if (R[i].dims) free(R[i].dims);
    free(R);
}

