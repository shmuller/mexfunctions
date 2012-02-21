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

typedef struct mdsdescrip {
	char dtype;
	char ndims;
	int length;
	int *dims; 
    void *ptr;
} mdsDescrip;

typedef mxArray wrap_Array;
typedef mxClassID wrap_ClassID;
typedef mxComplexity wrap_Complexity;

wrap_ClassID wrap_getClassID(const wrap_Array *r) {
    return mxGetClassID(r);
}

wrap_Complexity wrap_getComplexity(const wrap_Array *r) {
    return mxIsComplex(r) ? mxCOMPLEX : mxREAL;
}

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

const wrap_dtype *mds2wrap(const char *dtype) {
    int i;
    const dtype_item *t;
    if (dtype==NULL) return NULL;
    for(i=dtype_table_len,t=dtype_table; i--; t++) {
        if (t->dtype==*dtype) return &t->w_dtype;
    }
    return NULL;
}

const char *wrap2mds(const wrap_dtype *w_dtype) {
    int i;
    const dtype_item *t;
    if (w_dtype==NULL) return NULL;
    for(i=dtype_table_len,t=dtype_table; i--; t++) {
        if (wrap_dtype_isequal(&t->w_dtype,w_dtype)) return &t->dtype;
    }
    return NULL;
}




void *getnumarg(const wrap_Array *r, wrap_ClassID id) {
	if (wrap_getClassID(r) == id) { 
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

/*
int mex2mds_type(wrap_ClassID ID, wrap_Complexity Co, char *dtype)
{
  switch (ID)
  {
	case mxCHAR_CLASS    :  *dtype = DTYPE_CSTRING;    break;
	case mxUINT8_CLASS   :  *dtype = DTYPE_UCHAR;      break;
	case mxINT8_CLASS    :  *dtype = DTYPE_CHAR;       break;
	case mxUINT16_CLASS  :  *dtype = DTYPE_USHORT;     break;
	case mxINT16_CLASS   :  *dtype = DTYPE_SHORT;      break;
	case mxUINT32_CLASS  :  *dtype = DTYPE_ULONG;      break;
	case mxINT32_CLASS   :  *dtype = DTYPE_LONG;       break;
	case mxUINT64_CLASS  :  *dtype = DTYPE_ULONGLONG;  break;
	case mxINT64_CLASS   :  *dtype = DTYPE_LONGLONG;   break;
	case mxSINGLE_CLASS  :  *dtype = (Co==mxREAL) ? DTYPE_FLOAT  : DTYPE_COMPLEX;        break;
	case mxDOUBLE_CLASS  :  *dtype = (Co==mxREAL) ? DTYPE_DOUBLE : DTYPE_COMPLEX_DOUBLE; break;
  }
  return(1);
}

int mds2mex_type(struct descrip *d, wrap_ClassID *ID, wrap_Complexity *Co)
{
  switch (d->dtype)
  {
	case DTYPE_CSTRING         :  *ID = mxCHAR_CLASS;    *Co = mxREAL;     break;
	case DTYPE_UCHAR           :  *ID = mxUINT8_CLASS;   *Co = mxREAL;     break;
	case DTYPE_CHAR            :  *ID = mxINT8_CLASS;    *Co = mxREAL;     break;
	case DTYPE_USHORT          :  *ID = mxUINT16_CLASS;  *Co = mxREAL;     break;
	case DTYPE_SHORT           :  *ID = mxINT16_CLASS;   *Co = mxREAL;     break;
	case DTYPE_ULONG           :  *ID = mxUINT32_CLASS;  *Co = mxREAL;     break;
	case DTYPE_LONG            :  *ID = mxINT32_CLASS;   *Co = mxREAL;     break;
	case DTYPE_ULONGLONG       :  *ID = mxUINT64_CLASS;  *Co = mxREAL;     break;
	case DTYPE_LONGLONG        :  *ID = mxINT64_CLASS;   *Co = mxREAL;     break;
	case DTYPE_FLOAT           :  *ID = mxSINGLE_CLASS;  *Co = mxREAL;     break;
	case DTYPE_DOUBLE          :  *ID = mxDOUBLE_CLASS;  *Co = mxREAL;     break;
	case DTYPE_COMPLEX         :  *ID = mxSINGLE_CLASS;  *Co = mxCOMPLEX;  break;
	case DTYPE_COMPLEX_DOUBLE  :  *ID = mxDOUBLE_CLASS;  *Co = mxCOMPLEX;  break;
	default                    :  *ID = mxUNKNOWN_CLASS; *Co = mxREAL;     break;
  }
  return(1);
}
*/

int mds2mex_dims(struct descrip *d, char *ndims, mwSize **dims)
{
	int i;
	*ndims = max(d->ndims,2);
	*dims = malloc(*ndims*sizeof(mwSize));

	for(i=0; i<d->ndims; i++) (*dims)[i] = d->dims[i];
	for(; i<*ndims; i++) (*dims)[i] = 1;
	return(1);
}

void *mex2mds_cmplx(const wrap_Array *r) {
	size_t i,num,siz;
	void *pr,*pi,*buf,*b;

	num = wrap_getNumberOfElements(r);
	siz = wrap_getElementSize(r);
	pr  = mxGetData(r);
	pi  = mxGetImagData(r);
	buf = malloc(2*num*siz);

	for(i=0,b=buf; i<num; i++,b+=2*siz,pr+=siz,pi+=siz) {
		memcpy(b    ,pr,siz);
		memcpy(b+siz,pi,siz);
	}
	return buf;
}

void mds2mex_cmplx(const wrap_Array *r, void *buf) {
	size_t i,num,siz;
	void *pr,*pi,*b;

	num = wrap_getNumberOfElements(r);
	siz = wrap_getElementSize(r);
	pr  = mxGetData(r);
	pi  = mxGetImagData(r);

	for(i=0,b=buf; i<num; i++,b+=2*siz,pr+=siz,pi+=siz) {
		memcpy(pr,b    ,siz);
		memcpy(pi,b+siz,siz);
	}
}

void getmdsDescrip(const wrap_Array *r, mdsDescrip *D) {
	int i;
	const mwSize *dimsR = mxGetDimensions(r);
	D->ndims = wrap_getNumberOfDimensions(r);
    /*
	wrap_ClassID ID = wrap_getClassID(r); 
	wrap_Complexity Co = wrap_getComplexity(r);
	mex2mds_type(ID,Co,&D->dtype);
	*/
    
    wrap_dtype w_dtype = wrap_getdtype(r);
    D->dtype = *wrap2mds(&w_dtype);
    
    wrap_ClassID ID = w_dtype.ID;
    wrap_Complexity Co = w_dtype.Co;
    
	if (ID != mxCHAR_CLASS) {
		/* remove singleton dimensions */
		for(i=D->ndims-1; i>=0; i--) if (dimsR[i]==1) (D->ndims)--; else break;

		D->dims = calloc(D->ndims,sizeof(int));
		for(i=0; i<D->ndims; i++) (D->dims)[i]=dimsR[i];

		if (Co==mxREAL) {
			D->ptr = mxGetData(r);
		} else {
			if (ID != mxDOUBLE_CLASS && ID != mxSINGLE_CLASS) {
				wrap_error("Complex data must be single or double");
			}
			D->ptr = mex2mds_cmplx(r);
		}
	} else {
		D->ndims = 0; D->dims = NULL;
		D->ptr = getstringarg(r);
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

	char ndims;
	int *dims;
	mwSize *dimsL;
	wrap_ClassID ID;
	wrap_Complexity Co;

	for(i=2; i<nR; i++) {
		getmdsDescrip(R[i],D);
		arg = MakeDescrip(&exparg,D->dtype,D->ndims,D->dims,D->ptr);
		stat = SendArg(sock, i-2, arg->dtype, nR-2, ArgLen(arg), arg->ndims, arg->dims, arg->ptr);
		if (D->dims) free(D->dims);
	}
	
	stat = GetAnswerInfoTS(sock, &arg->dtype, &arg->length, &arg->ndims, arg->dims, &numbytes, &arg->ptr, &mem);
    /*
    stat = mds2mex_type(arg,&ID,&Co);
    */
    const wrap_dtype *w_dtype = mds2wrap(&arg->dtype);
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
		if (Co==mxREAL) {
			memcpy(mxGetData(L[0]),arg->ptr,numbytes);
		} else {
			mds2mex_cmplx(L[0],arg->ptr);
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

