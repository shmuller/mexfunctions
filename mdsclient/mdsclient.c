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

#include "stdio.h"
#include "stdlib.h"

#include "mex.h"
#include "matrix.h"
#include "ipdesc.h"
#include "string.h"

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

void* getnumarg(const mxArray *r, mxClassID id) {
	if (mxGetClassID(r) == id) { 
		return((void*)mxGetPr(r));
	} else {
		mexErrMsgTxt("Wrong argument type");
	}
}

char *getstringarg(const mxArray *r) {
	short len = mxGetNumberOfElements(r);
	char *str = malloc((len+1)*sizeof(char));

	mxGetString(r,str,len+1);
	return(str);
}

void *getarg(const mxArray *r, mxClassID *mxID, char *ndims, const mwSize **dims) {
    *mxID = mxGetClassID(r);
    *ndims = mxGetNumberOfDimensions(r);
    *dims = mxGetDimensions(r);
    if (*mxID != mxCHAR_CLASS) {
        return((void*)mxGetPr(r));
    } else {
        return((void*)getstringarg(r));
    }
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

	for(i=0; i<d->ndims; i++) (*dims)[i] = max(d->dims[i],1);
	for(; i<*ndims; i++) (*dims)[i] = 1;
	return(1);
}


int sm_mdsconnect(int nL, mxArray *L[], int nR, const mxArray *R[]) {
	char *serv = getstringarg(R[1]);
	SOCKET sock = ConnectToMds(serv);
	if (sock == INVALID_SOCKET) {
		mexErrMsgTxt("Could not connect to server");
	}
	L[0] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
	*((int*)mxGetPr(L[0])) = sock;
	return(1);
}

int sm_mdsopen(int nL, mxArray *L[], int nR, const mxArray *R[]) {
	SOCKET sock = *((int*)getnumarg(R[1],mxINT32_CLASS));
	char *tree = getstringarg(R[2]);
	int shot = *((int*)getnumarg(R[3],mxINT32_CLASS));
	int stat = MdsOpen(sock,tree,shot);
/*	if (!status_ok(stat)) {
		mexErrMsgTxt("Could not open tree");
	}
*/	return(stat);
}

int sm_mdsclose(int nL, mxArray *L[], int nR, const mxArray *R[]) {
	SOCKET sock = *((int*)getnumarg(R[1],mxINT32_CLASS));
	int stat = MdsClose(sock);
/*	if (!status_ok(stat)) {
		mexErrMsgTxt("Could not close tree");
	}
*/	return(stat);
}

int sm_mdsdisconnect(int nL, mxArray *L[], int nR, const mxArray *R[]) {
	SOCKET sock = *((int*)getnumarg(R[1],mxINT32_CLASS));
	int stat = DisconnectFromMds(sock);
/*	if (!status_ok(stat)) {
		mexErrMsgTxt("Could not disconnect from server");
	}
*/	return(stat);
}

int sm_mdsvalue(int nL, mxArray *L[], int nR, const mxArray *R[]) {
	SOCKET sock = *((int*)getnumarg(R[1],mxINT32_CLASS));
	void *mxArg;
    mxClassID mxID;
	mxComplexity mxCo;
	
	struct descrip exparg, *arg;
	int idx = 0, nargs = 1, numbytes = 0, stat = 0;
	void *mem = 0, *out = 0;

	char ndims;
	const mwSize *dimsR;
    mwSize *dimsL;

    mxArg = getarg(R[2],&mxID,&ndims,&dimsR);
    
	arg = MakeDescrip(&exparg,DTYPE_CSTRING,0,0,mxArg);
	stat = SendArg(sock, idx, arg->dtype, nargs, ArgLen(arg), arg->ndims, arg->dims, arg->ptr);

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
		out = (void*) mxGetPr(L[0]);
		memcpy(out,arg->ptr,numbytes);
	}
	if (mem) free(mem);
	return(1);
}



void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
	char *cmd = getstringarg(R[0]);
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

	if (!status_ok(stat)) {
		mexErrMsgTxt("Untrapped error occurred");
	}
}

