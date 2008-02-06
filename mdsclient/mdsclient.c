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
	int dtype_long = DTYPE_LONG;
	char dtype=DTYPE_CSTRING;
	int nbytes=0;

/*	int idesc = descr2(&dtype_long, &zero);*/
    
	register int i;
	unsigned char nargs = 0;
	char ndims = mxGetNumberOfDimensions(R[0]);
	const mwSize *Dims = mxGetDimensions(R[0]);
	short len = mxGetNumberOfElements(R[0]);
	char *ptr = malloc((len+1)*sizeof(char));

	mxGetString(R[0],ptr,len+1);
	printf("%s\n",ptr);    

	int *dims = malloc(ndims*sizeof(int));
	for(i=0; i<ndims; i++) dims[i] = (int)Dims[i];

	int *out;
    
	L[0] = mxCreateNumericArray(1,Dims,mxINT32_CLASS,mxREAL);
	out = (int*) mxGetPr(L[0]);
    
	sock = ConnectToMds("plaspc03.ucsd.edu");
	if (sock == INVALID_SOCKET) {
		mxErrMsgTxt("Could not connect to server");
	}

	stat = MdsOpen(sock,"csdx",shot);

	struct descrip ans_arg;
	
	stat = MdsValue(sock, ptr, &ans_arg);

	*out = *((int*)ans_arg.ptr);


/*
	struct descrip exparg, *arg, *ans_arg;
	ans_arg = arg = MakeDescrip(&exparg,DTYPE_CSTRING,0,0,ptr);

	printf("dtype=%d, len=%d, ndims=%d, dims[0]=%d, numbytes=%d\n",arg->dtype,ArgLen(arg),arg->ndims,arg->dims[0],0);
	
	int idx=0;
    	stat = SendArg(sock, idx, arg->dtype, nargs, ArgLen(arg), arg->ndims, arg->dims, arg->ptr);

    	int numbytes;
    	void *dptr;
    	void *mem = 0;
    	stat = GetAnswerInfoTS(sock, &ans_arg->dtype, &len, &ans_arg->ndims, ans_arg->dims, &numbytes, &dptr, &mem);
	ans_arg->length = len;

	printf("dtype=%d, len=%d, ndims=%d, dims[0]=%d, numbytes=%d\n",ans_arg->dtype,ans_arg->length,ans_arg->ndims,ans_arg->dims[0],numbytes);
	
	*out = *((int*)dptr);
*/

/*	stat = MdsValue2(ptr,&idesc,out,&zero);*/
	
	stat = MdsClose(sock);

	stat = DisconnectFromMds(sock);


}
