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
#include "math.h"
#include "mdslib.h"

#define min(a,b) ((a)<(b) ? (a) : (b))
#define max(a,b) ((a)>(b) ? (a) : (b))

#define status_ok(status) (((status) & 1) == 1)


void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
	SOCKET sock;
	int stat,shot=0,zero=0;
	int dtype_long = DTYPE_LONG;
	char dtype=DTYPE_CSTRING;
	int nbytes=0;
	void **dptr;

  	int idesc = descr2(&dtype_long, &zero);
    
	register int i;
	unsigned char nargs = 0;
	char ndims = mxGetNumberOfDimensions(R[0]);
	const mwSize *Dims = mxGetDimensions(R[0]);
	short len = mxGetNumberOfElements(R[0]);
	char *ptr = (char*) mxGetPr(R[0]);
    
	int *dims = malloc(ndims*sizeof(int));
	for(i=0; i<ndims; i++) dims[i] = (int)Dims[i];

	int *out;
    
	L[0] = mxCreateNumericArray(ndims,Dims,mxINT32_CLASS,mxREAL);
	out = (int*) mxGetPr(L[0]);
    
	sock = MdsConnect("plaspc03.ucsd.edu");
	if (sock == INVALID_SOCKET) {
		mxErrMsgTxt("Could not connect to server");
	}

	stat = MdsOpen("csdx",&shot);
/*	
	printf("dtype=%d, len=%d, ndims=%d, dims[0]=%d, nbytes=%d\n",dtype,len,ndims,dims[0],nbytes);

	stat = SendArg(sock, 0, dtype, nargs, len, ndims, dims, ptr);
	if (!status_ok(stat)) {
		mxErrMsgTxt("SendArg error");
	}

	stat = GetAnswerInfo(sock, &dtype, &len, &ndims, dims, &nbytes, dptr);
	if (!status_ok(stat)) {
		mxErrMsgTxt("GetAnswerInfo error");
	}

	printf("dtype=%d, len=%d, ndims=%d, dims[0]=%d, nbytes=%d\n",dtype,len,ndims,dims[0],nbytes);
*/
	stat = MdsValue2("$shot",&idesc,out,&zero);
	
	stat = MdsClose("csdx",&shot);

	stat = DisconnectFromMds(sock);


}
