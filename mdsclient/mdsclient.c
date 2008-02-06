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
#include "ipdesc.h"

#define min(a,b) ((a)<(b) ? (a) : (b))
#define max(a,b) ((a)>(b) ? (a) : (b))

#define status_ok(status) (((status) & 1) == 1)

#ifndef INVALID_SOCKET
#define INVALID_SOCKET -1
#endif

void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
	SOCKET sock;
	int stat,shot=0,zero=0;    

	char ndims = mxGetNumberOfDimensions(R[0]);
	const mwSize *dims = mxGetDimensions(R[0]);
	short len = mxGetNumberOfElements(R[0]);
	char *expression = malloc((len+1)*sizeof(char));

	mxGetString(R[0],expression,len+1);
	printf("%s\n",expression);    

	L[0] = mxCreateNumericArray(1,dims,mxINT32_CLASS,mxREAL);
	void *out = (void*) mxGetPr(L[0]);
    
	sock = ConnectToMds("plaspc03.ucsd.edu");
	if (sock == INVALID_SOCKET) {
		mexErrMsgTxt("Could not connect to server");
	}

	stat = MdsOpen(sock,"csdx",shot);

	int idx=0, nargs=0;
	struct descrip exparg, *arg=&exparg;
/*
	arg = MakeDescrip((struct descrip *)&exparg,DTYPE_CSTRING,0,0,expression);
	status = SendArg(sock, idx, arg->dtype, nargs, ArgLen(arg), arg->ndims, arg->dims, arg->ptr);
*/	
	stat = MdsValue(sock, expression, arg);
	*((int*)out) = *((int*)arg->ptr);

	stat = MdsClose(sock);

	stat = DisconnectFromMds(sock);

}
