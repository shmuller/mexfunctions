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

void *getnumarg(const mxArray *r, mxClassID id) {
	if (mxGetClassID(r) == id) { 
		return mxGetData(r);
	} else {
		mexErrMsgTxt("Wrong argument type");
	}
}

void *getstringarg(const mxArray *r) {
	short len = mxGetNumberOfElements(r);
	void *str = malloc((len+1)*sizeof(char));

	mxGetString(r,str,len+1);
	return str;
}

int mex2mds_type(mxClassID mxID, mxComplexity mxCo, char *dtype)
{
  switch (mxID)
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
	case mxSINGLE_CLASS  :  *dtype = (mxCo==mxREAL) ? DTYPE_FLOAT  : DTYPE_COMPLEX;        break;
	case mxDOUBLE_CLASS  :  *dtype = (mxCo==mxREAL) ? DTYPE_DOUBLE : DTYPE_COMPLEX_DOUBLE; break;
  }
  return(1);
}

int mds2mex_type(struct descrip *d, mxClassID *mxID, mxComplexity *mxCo)
{
  switch (d->dtype)
  {
	case DTYPE_CSTRING         :  *mxID = mxCHAR_CLASS;    *mxCo = mxREAL;     break;
	case DTYPE_UCHAR           :  *mxID = mxUINT8_CLASS;   *mxCo = mxREAL;     break;
	case DTYPE_CHAR            :  *mxID = mxINT8_CLASS;    *mxCo = mxREAL;     break;
	case DTYPE_USHORT          :  *mxID = mxUINT16_CLASS;  *mxCo = mxREAL;     break;
	case DTYPE_SHORT           :  *mxID = mxINT16_CLASS;   *mxCo = mxREAL;     break;
	case DTYPE_ULONG           :  *mxID = mxUINT32_CLASS;  *mxCo = mxREAL;     break;
	case DTYPE_LONG            :  *mxID = mxINT32_CLASS;   *mxCo = mxREAL;     break;
	case DTYPE_ULONGLONG       :  *mxID = mxUINT64_CLASS;  *mxCo = mxREAL;     break;
	case DTYPE_LONGLONG        :  *mxID = mxINT64_CLASS;   *mxCo = mxREAL;     break;
	case DTYPE_FLOAT           :  *mxID = mxSINGLE_CLASS;  *mxCo = mxREAL;     break;
	case DTYPE_DOUBLE          :  *mxID = mxDOUBLE_CLASS;  *mxCo = mxREAL;     break;
	case DTYPE_COMPLEX         :  *mxID = mxSINGLE_CLASS;  *mxCo = mxCOMPLEX;  break;
	case DTYPE_COMPLEX_DOUBLE  :  *mxID = mxDOUBLE_CLASS;  *mxCo = mxCOMPLEX;  break;
	default                    :  *mxID = mxUNKNOWN_CLASS; *mxCo = mxREAL;     break;
  }
  return(1);
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

void *mex2mds_cmplx(const mxArray *r) {
	size_t i,num,siz;
	void *pr,*pi,*buf,*b;

	num = mxGetNumberOfElements(r);
	siz = mxGetElementSize(r);
	pr  = mxGetData(r);
	pi  = mxGetImagData(r);
	buf = malloc(2*num*siz);

	for(i=0,b=buf; i<num; i++,b+=2*siz,pr+=siz,pi+=siz) {
		memcpy(b    ,pr,siz);
		memcpy(b+siz,pi,siz);
	}
	return buf;
}

void mds2mex_cmplx(const mxArray *r, void *buf) {
	size_t i,num,siz;
	void *pr,*pi,*b;

	num = mxGetNumberOfElements(r);
	siz = mxGetElementSize(r);
	pr  = mxGetData(r);
	pi  = mxGetImagData(r);

	for(i=0,b=buf; i<num; i++,b+=2*siz,pr+=siz,pi+=siz) {
		memcpy(pr,b    ,siz);
		memcpy(pi,b+siz,siz);
	}
}

void getmdsDescrip(const mxArray *r, mdsDescrip *D) {
	int i;
	const mwSize *dimsR = mxGetDimensions(r);
	D->ndims = mxGetNumberOfDimensions(r);
	mxClassID mxID = mxGetClassID(r); 
	mxComplexity mxCo = (mxIsComplex(r)) ? mxCOMPLEX : mxREAL;
	mex2mds_type(mxID,mxCo,&D->dtype);
	
	if (mxID != mxCHAR_CLASS) {
		/* remove singleton dimensions */
		for(i=D->ndims-1; i>=0; i--) if (dimsR[i]==1) (D->ndims)--; else break;

		D->dims = calloc(D->ndims,sizeof(int));
		for(i=0; i<D->ndims; i++) (D->dims)[i]=dimsR[i];

		if (mxCo==mxREAL) {
			D->ptr = mxGetData(r);
		} else {
			if (mxID != mxDOUBLE_CLASS && mxID != mxSINGLE_CLASS) {
				mexErrMsgTxt("Complex data must be single or double");
			}
			D->ptr = mex2mds_cmplx(r);
		}
	} else {
		D->ndims = 0; D->dims = NULL;
		D->ptr = getstringarg(r);
	}
}

#include "tcp.c"

int sm_mdsconnect(int nL, mxArray *L[], int nR, const mxArray *R[]) {
	char *host = getstringarg(R[1]);
    int sock;
    
    switch (sock=tcpconnect(host)) {
        case -1:
        case -2:
        case -3:
        case -4: mexErrMsgTxt("Could not connect to server");
        case -5: mexErrMsgTxt("Could not authenticate user");
    }
    
	L[0] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
	*((int*)mxGetData(L[0])) = sock;
	return(1);
}

int sm_mdsdisconnect(int nL, mxArray *L[], int nR, const mxArray *R[]) {
	int sock = *((int*)getnumarg(R[1],mxINT32_CLASS));
	tcpdisconnect(sock);
	return(1);
}

int sm_mdsopen(int nL, mxArray *L[], int nR, const mxArray *R[]) {
	SOCKET sock = *((int*)getnumarg(R[1],mxINT32_CLASS));
	char *tree = getstringarg(R[2]);
	int shot = *((int*)getnumarg(R[3],mxINT32_CLASS));
	int stat = MdsOpen(sock,tree,shot);
	if (!status_ok(stat)) {
		mexErrMsgTxt("Could not open tree");
	}
	return(1);
}

int sm_mdsclose(int nL, mxArray *L[], int nR, const mxArray *R[]) {
	SOCKET sock = *((int*)getnumarg(R[1],mxINT32_CLASS));
	int stat = MdsClose(sock);
	if (!status_ok(stat)) {
		mexErrMsgTxt("Could not close tree");
	}
	return(1);
}

int sm_mdsvalue(int nL, mxArray *L[], int nR, const mxArray *R[]) {
	SOCKET sock = *((int*)getnumarg(R[1],mxINT32_CLASS));
	mdsDescrip *D = malloc(sizeof(mdsDescrip));
	
	struct descrip exparg, *arg;
	int i = 0, numbytes = 0, stat = 0;
	void *mem = 0, *out = 0;

	char ndims;
	int *dims;
	mwSize *dimsL;
	mxClassID mxID;
	mxComplexity mxCo;

	for(i=2; i<nR; i++) {
		getmdsDescrip(R[i],D);
		arg = MakeDescrip(&exparg,D->dtype,D->ndims,D->dims,D->ptr);
		stat = SendArg(sock, i-2, arg->dtype, nR-2, ArgLen(arg), arg->ndims, arg->dims, arg->ptr);
		if (D->dims) free(D->dims);
	}
	
	stat = GetAnswerInfoTS(sock, &arg->dtype, &arg->length, &arg->ndims, arg->dims, &numbytes, &arg->ptr, &mem);
	stat = mds2mex_type(arg,&mxID,&mxCo);

	if (mxID == mxUNKNOWN_CLASS) {
		L[0] = mxCreateNumericArray(0,0,mxDOUBLE_CLASS,mxREAL);
	} else if (mxID == mxCHAR_CLASS) {
		out = calloc(numbytes+1,sizeof(char));
		memcpy(out,arg->ptr,numbytes);
		L[0] = mxCreateString((char*)out);
	} else {
		stat = mds2mex_dims(arg,&ndims,&dimsL);
		L[0] = mxCreateNumericArray(ndims,dimsL,mxID,mxCo);
		if (mxCo==mxREAL) {
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
		mexErrMsgTxt("Unknown command");
	}

	if (stat < 0) {
		sprintf(errstr,"Untrapped error occurred: %d",stat);
		mexErrMsgTxt(errstr);
	}
}

