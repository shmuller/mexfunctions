/* mdsclient
 *
 * Compile on Linux: 
 * mex -v mdsclient.c COPTIMFLAGS=-O3 LDOPTIMFLAGS=-O3 -I/usr/local/mdsplus/include -L/usr/local/mdsplus/lib -lMdsIpShr
 * 
 * Compile on Windows:
 * mex -v mdsclient.c OPTIMFLAGS=-O3 -I"C:\PROGRA~1\MDSplus\DEVTOOLS\include" -L"C:\PROGRA~1\MDSplus\DEVTOOLS\lib" -lMdsIpShr
 *
 * S. H. Muller, 2008/02/05
 */

#include <stdio.h>
#include <stdlib.h>

#include <mex.h>
#include <matrix.h>
#include <ipdesc.h>
#include <string.h>

#ifndef min
#define min(a,b) ((a)<(b) ? (a) : (b))
#endif
#ifndef max
#define max(a,b) ((a)>(b) ? (a) : (b))
#endif

#define status_ok(status) (((status) & 1) == 1)

#ifndef INVALID_SOCKET
#define INVALID_SOCKET -1
#endif

#ifdef COMPAT_V70
typedef int mwSize;
#endif

typedef struct {
	char dtype;
	char ndims;
	int length;
	int *dims; 
    void *ptr;
} mdsDescrip;

typedef mxArray wrap_Array;
typedef mxClassID wrap_ClassID;
typedef mxComplexity wrap_Complexity;

typedef struct {
    mxClassID ID;
    mxComplexity Co;
} wrap_dtype;

int wrap_dtype_isequal(const wrap_dtype *a, const wrap_dtype *b) {
    return (a->ID == b->ID && a->Co == b->Co);
}

typedef struct {
    wrap_dtype w_dtype;
    char dtype;
} dtype_item;

static const dtype_item dtype_table[] = {
    {{mxCHAR_CLASS  , mxREAL   }, DTYPE_CSTRING       },
    {{mxUINT8_CLASS , mxREAL   }, DTYPE_UCHAR         },
    {{mxINT8_CLASS  , mxREAL   }, DTYPE_CHAR          },
    {{mxUINT16_CLASS, mxREAL   }, DTYPE_USHORT        },
    {{mxINT16_CLASS , mxREAL   }, DTYPE_SHORT         },
    {{mxUINT32_CLASS, mxREAL   }, DTYPE_ULONG         },
    {{mxINT32_CLASS , mxREAL   }, DTYPE_LONG          },
    {{mxUINT64_CLASS, mxREAL   }, DTYPE_ULONGLONG     },
    {{mxINT64_CLASS , mxREAL   }, DTYPE_LONGLONG      },
    {{mxSINGLE_CLASS, mxREAL   }, DTYPE_FLOAT         },
    {{mxSINGLE_CLASS, mxCOMPLEX}, DTYPE_COMPLEX       },
    {{mxDOUBLE_CLASS, mxREAL   }, DTYPE_DOUBLE        },
    {{mxDOUBLE_CLASS, mxCOMPLEX}, DTYPE_COMPLEX_DOUBLE}
};

#define LEN(x) (sizeof(x)/sizeof((x)[0]))

static const int dtype_table_len = LEN(dtype_table);

wrap_dtype wrap_getdtype(const wrap_Array *r) {
    wrap_dtype w_dtype; 
    w_dtype.ID = mxGetClassID(r);
    w_dtype.Co = mxIsComplex(r) ? mxCOMPLEX : mxREAL;
    return w_dtype;
}

typedef struct {
    wrap_dtype w_dtype;
    int siz;
    char ndims;
	int num;
	mwSize *dims; 
    void *pr;
    void *pi;
} wrap_Descrip;

wrap_Descrip wrap_getDescrip(const wrap_Array *r) {
    wrap_Descrip w_D;
    w_D.w_dtype = wrap_getdtype(r);
    w_D.siz     = mxGetElementSize(r);
    w_D.ndims   = mxGetNumberOfDimensions(r);
    w_D.num     = mxGetNumberOfElements(r);
    w_D.dims    = (mwSize*) mxGetDimensions(r);
    
    if (w_D.w_dtype.ID != mxCHAR_CLASS) {
        w_D.pr = mxGetData(r);
        w_D.pi = mxGetImagData(r);
    } else {
        w_D.pr = malloc((w_D.num+1)*sizeof(char));
	    mxGetString(r, w_D.pr, w_D.num+1);
    }
    return w_D;
}


void wrap_error(char *err) {
    mexErrMsgTxt(err);
}

int wrap_getNumberOfElements(const wrap_Array *r) {
    return mxGetNumberOfElements(r);
}

int wrap_getElementSize(const wrap_Array *r) {
    return mxGetElementSize(r);
}

int wrap_getNumberOfDimensions(const wrap_Array *r) {
    return mxGetNumberOfDimensions(r);
}


/* completely wrapped code */

const wrap_dtype *mds2wrap_dtype(const char *dtype) {
    int i;
    const dtype_item *t;
    if (dtype==NULL) return NULL;
    for(i=dtype_table_len,t=dtype_table; i--; t++) {
        if (t->dtype==*dtype) return &t->w_dtype;
    }
    return NULL;
}

const char *wrap2mds_dtype(const wrap_dtype *w_dtype) {
    int i;
    const dtype_item *t;
    if (w_dtype==NULL) return NULL;
    for(i=dtype_table_len,t=dtype_table; i--; t++) {
        if (wrap_dtype_isequal(&t->w_dtype,w_dtype)) return &t->dtype;
    }
    return NULL;
}



void *getnumarg(const wrap_Array *r, wrap_ClassID id) {
	if (wrap_getdtype(r).ID == id) { 
		return mxGetData(r);
	} else {
		wrap_error("Wrong argument type");
	}
}

void *getstringarg(const wrap_Array *r) {
	short len = wrap_getNumberOfElements(r);
	void *str = malloc((len+1)*sizeof(char));

	mxGetString(r,str,len+1);
	return str;
}


int mds2mex_dims(struct descrip *d, char *ndims, mwSize **dims)
{
	int i;
	*ndims = max(d->ndims,2);
	*dims = malloc(*ndims*sizeof(mwSize));

	for(i=0; i<d->ndims; i++) (*dims)[i] = d->dims[i];
	for(; i<*ndims; i++) (*dims)[i] = 1;
	return(1);
}

void *mex2mds_cmplx(const wrap_Descrip *w_D) {
	size_t i,num,siz;
	void *pr,*pi,*buf,*b;
    
    num = w_D->num;
    siz = w_D->siz;
    pr  = w_D->pr;
    pi  = w_D->pi;
    
	buf = malloc(2*num*siz);

	for(i=0,b=buf; i<num; i++,b+=2*siz,pr+=siz,pi+=siz) {
		memcpy(b    ,pr,siz);
		memcpy(b+siz,pi,siz);
	}
	return buf;
}

void mds2mex_cmplx(const wrap_Descrip *w_D, void *buf) {
	size_t i,num,siz;
	void *pr,*pi,*b;
    
    num = w_D->num;
    siz = w_D->siz;
    pr  = w_D->pr;
    pi  = w_D->pi;

	for(i=0,b=buf; i<num; i++,b+=2*siz,pr+=siz,pi+=siz) {
		memcpy(pr,b    ,siz);
		memcpy(pi,b+siz,siz);
	}
}

void getmdsDescrip(const wrap_Descrip *w_D, mdsDescrip *D) {
	int i;
	mwSize *dimsR = w_D->dims;
	D->ndims = w_D->ndims;
    
    wrap_dtype w_dtype = w_D->w_dtype;
    D->dtype = *wrap2mds_dtype(&w_dtype);
    
    wrap_ClassID ID = w_dtype.ID;
    wrap_Complexity Co = w_dtype.Co;
    
	if (ID != mxCHAR_CLASS) {
		/* remove singleton dimensions */
		for(i=D->ndims-1; i>=0; i--) if (dimsR[i]==1) (D->ndims)--; else break;

		D->dims = calloc(D->ndims,sizeof(int));
		for(i=0; i<D->ndims; i++) (D->dims)[i]=dimsR[i];

		if (Co==mxREAL) {
			D->ptr = w_D->pr;
		} else {
			if (ID != mxDOUBLE_CLASS && ID != mxSINGLE_CLASS) {
				wrap_error("Complex data must be single or double");
			}
			D->ptr = mex2mds_cmplx(w_D);
		}
	} else {
		D->ndims = 0; D->dims = NULL;
		D->ptr = w_D->pr;
	}
}

#include "tcp.c"

int sm_mdsconnect(int nL, wrap_Array *L[], int nR, const wrap_Array *R[]) {
	char *host = getstringarg(R[1]);
    int sock;
    
    switch (sock=tcpconnect(host)) {
        case -1:
        case -2:
        case -3:
        case -4: wrap_error("Could not connect to server");
        case -5: wrap_error("Could not authenticate user");
    }
    
	L[0] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
	*((int*)mxGetData(L[0])) = sock;
	return(1);
}

int sm_mdsdisconnect(int nL, wrap_Array *L[], int nR, const wrap_Array *R[]) {
	int sock = *((int*)getnumarg(R[1],mxINT32_CLASS));
	tcpdisconnect(sock);
	return(1);
}

int sm_mdsopen(int nL, wrap_Array *L[], int nR, const wrap_Array *R[]) {
	SOCKET sock = *((int*)getnumarg(R[1],mxINT32_CLASS));
	char *tree = getstringarg(R[2]);
	int shot = *((int*)getnumarg(R[3],mxINT32_CLASS));
	int stat = MdsOpen(sock,tree,shot);
	if (!status_ok(stat)) {
		wrap_error("Could not open tree");
	}
	return(1);
}

int sm_mdsclose(int nL, wrap_Array *L[], int nR, const wrap_Array *R[]) {
	SOCKET sock = *((int*)getnumarg(R[1],mxINT32_CLASS));
	int stat = MdsClose(sock);
	if (!status_ok(stat)) {
		wrap_error("Could not close tree");
	}
	return(1);
}

int sm_mdsvalue(int nL, wrap_Array *L[], int nR, const wrap_Array *R[]) {
	SOCKET sock = *((int*)getnumarg(R[1],mxINT32_CLASS));
	mdsDescrip *D = malloc(sizeof(mdsDescrip));
	
	struct descrip exparg, *arg;
	int i = 0, numbytes = 0, stat = 0;
	void *mem = 0, *out = 0;

    wrap_Descrip w_D;
	char ndims;
	int *dims;
	mwSize *dimsL;
	wrap_ClassID ID;
	wrap_Complexity Co;

	for(i=2; i<nR; i++) {
        w_D = wrap_getDescrip(R[i]);
		getmdsDescrip(&w_D,D);
		arg = MakeDescrip(&exparg,D->dtype,D->ndims,D->dims,D->ptr);
		stat = SendArg(sock, i-2, arg->dtype, nR-2, ArgLen(arg), arg->ndims, arg->dims, arg->ptr);
		if (D->dims) free(D->dims);
	}
	
	stat = GetAnswerInfoTS(sock, &arg->dtype, &arg->length, &arg->ndims, arg->dims, &numbytes, &arg->ptr, &mem);
    
    const wrap_dtype *w_dtype = mds2wrap_dtype(&arg->dtype);
    ID = w_dtype->ID;
    Co = w_dtype->Co;
    
	if (ID == mxUNKNOWN_CLASS) {
		L[0] = mxCreateNumericArray(0,0,mxDOUBLE_CLASS,mxREAL);
	} else if (ID == mxCHAR_CLASS) {
		out = calloc(numbytes+1,sizeof(char));
		memcpy(out,arg->ptr,numbytes);
		L[0] = mxCreateString((char*)out);
	} else {
		stat = mds2mex_dims(arg,&ndims,&dimsL);
		L[0] = mxCreateNumericArray(ndims,dimsL,ID,Co);
        w_D = wrap_getDescrip(L[0]);
		if (Co==mxREAL) {
			memcpy(w_D.pr,arg->ptr,numbytes);
		} else {
			mds2mex_cmplx(&w_D,arg->ptr);
		}   
	}
	if (mem) free(mem);
	free(D);
	return(1);
}



void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
	char *cmd = getstringarg(R[0]), errstr[256];
	int stat;

	if (strcmp(cmd,"mdsvalue")==0) {
		stat = sm_mdsvalue(nL,L,nR,R);
	} else if (strcmp(cmd,"mdsconnect")==0) {
		stat = sm_mdsconnect(nL,L,nR,R);
	} else if (strcmp(cmd,"mdsopen")==0) {
		stat = sm_mdsopen(nL,L,nR,R);
	} else if (strcmp(cmd,"mdsclose")==0) {
		stat = sm_mdsclose(nL,L,nR,R);
	} else if (strcmp(cmd,"mdsdisconnect")==0) {
		stat = sm_mdsdisconnect(nL,L,nR,R);
	} else {
		wrap_error("Unknown command");
	}

	if (stat < 0) {
		sprintf(errstr,"Untrapped error occurred: %d",stat);
		wrap_error(errstr);
	}
}

