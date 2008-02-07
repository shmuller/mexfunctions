/* mdsclient
 *
 * Compile on Linux: 
 * mex -v COPTIMFLAGS=-O3 LDOPTIMFLAGS=-O3 mdsclient.c
 * 
 * Compile on Windows
 * mex -v OPTIMFLAGS=-O3 mdsclient.c
 *
 * S. H. Muller, 2008/02/05
 */

#include "stdio.h"
#include "stdlib.h"

#include "mex.h"
#include "matrix.h"
#include "ipdesc.h"
#include "string.h"

#define min(a,b) ((a)<(b) ? (a) : (b))
#define max(a,b) ((a)>(b) ? (a) : (b))

#define status_ok(status) (((status) & 1) == 1)

#ifndef INVALID_SOCKET
#define INVALID_SOCKET -1
#endif



int mds2mex_type(struct descrip *d, mxClassID *mxID, mxComplexity *mxCo)
{
  switch (d->dtype)
  {
	case DTYPE_CSTRING		:  *mxID = -1;             *mxCo = mxREAL; break;
	case DTYPE_UCHAR		:  *mxID = mxUINT8_CLASS;  *mxCo = mxREAL; break;
	case DTYPE_CHAR			:  *mxID = mxINT8_CLASS;   *mxCo = mxREAL; break;
	case DTYPE_USHORT		:  *mxID = mxUINT16_CLASS; *mxCo = mxREAL; break;
	case DTYPE_SHORT		:  *mxID = mxINT16_CLASS;  *mxCo = mxREAL; break;
	case DTYPE_ULONG		:  *mxID = mxUINT32_CLASS; *mxCo = mxREAL; break;
	case DTYPE_LONG			:  *mxID = mxINT32_CLASS;  *mxCo = mxREAL; break;
	case DTYPE_ULONGLONG		:  *mxID = mxUINT64_CLASS; *mxCo = mxREAL; break;
	case DTYPE_LONGLONG		:  *mxID = mxINT64_CLASS;  *mxCo = mxREAL; break;
	case DTYPE_FLOAT		:  *mxID = mxSINGLE_CLASS; *mxCo = mxREAL; break;
	case DTYPE_DOUBLE		:  *mxID = mxDOUBLE_CLASS; *mxCo = mxREAL; break;
	case DTYPE_COMPLEX		:  *mxID = mxSINGLE_CLASS; *mxCo = mxCOMPLEX; break;
	case DTYPE_COMPLEX_DOUBLE	:  *mxID = mxDOUBLE_CLASS; *mxCo = mxCOMPLEX; break;
	default				:  *mxID = 0;              *mxCo = mxREAL; break;
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


void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
	SOCKET sock;
	int stat,shot=0,zero=0;
	char ndimsR = mxGetNumberOfDimensions(R[0]);
	const mwSize *dimsR = mxGetDimensions(R[0]);

	short len = mxGetNumberOfElements(R[0]);
	char *expression = malloc((len+1)*sizeof(char));
	mxGetString(R[0],expression,len+1);

	mxClassID mxID;
	mxComplexity mxCo;
	char ndims;
	mwSize *dims;
    
/*	sock = ConnectToMds("plaspc03.ucsd.edu");*/
	sock = ConnectToMds("localhost:8001");

	if (sock == INVALID_SOCKET) {
		mexErrMsgTxt("Could not connect to server");
	}

	stat = MdsOpen(sock,"csdx",shot);
	
	struct descrip exparg, *arg=&exparg;

	int idx=0, nargs=1;
    	int numbytes = 0;
    	void *mem = 0;
	void *out = 0;

	arg = MakeDescrip(&exparg,DTYPE_CSTRING,0,0,expression);
	stat = SendArg(sock, idx, arg->dtype, nargs, ArgLen(arg), arg->ndims, arg->dims, arg->ptr);
	
	stat = GetAnswerInfoTS(sock, &arg->dtype, &arg->length, &arg->ndims, arg->dims, &numbytes, &arg->ptr, &mem);

	stat = mds2mex_type(arg,&mxID,&mxCo);

	if (mxID == 0) {
		L[0] = mxCreateNumericArray(0,0,mxDOUBLE_CLASS,mxREAL);
	} else if (mxID == -1) {
		out = calloc(numbytes+1,sizeof(char));
		memcpy(out,arg->ptr,numbytes);
		L[0] = mxCreateString((char*)out);
	} else {
		stat = mds2mex_dims(arg,&ndims,&dims);
		L[0] = mxCreateNumericArray(ndims,dims,mxID,mxCo);
		out = (void*) mxGetPr(L[0]);
		memcpy(out,arg->ptr,numbytes);
	}

	stat = MdsClose(sock);

	stat = DisconnectFromMds(sock);
}

