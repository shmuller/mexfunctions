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
#include "string.h"

#define min(a,b) ((a)<(b) ? (a) : (b))
#define max(a,b) ((a)>(b) ? (a) : (b))

#define status_ok(status) (((status) & 1) == 1)

#ifndef INVALID_SOCKET
#define INVALID_SOCKET -1
#endif


int myMdsValue(SOCKET sock, char *expression, struct descrip *ans_arg);


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
	
	struct descrip exparg, *arg=&exparg;
/*
	int idx=0, nargs=1;
    	int numbytes = 0;
    	void *dptr;
    	void *mem = 0;

	arg = MakeDescrip(&exparg,DTYPE_CSTRING,0,0,expression);
	stat = SendArg(sock, idx, arg->dtype, nargs, ArgLen(arg), arg->ndims, arg->dims, arg->ptr);
	
	printf("dtype=%d, len=%d, ndims=%d, dims[0]=%d, numbytes=%d\n",arg->dtype,arg->length,arg->ndims,arg->dims[0],numbytes);

	stat = GetAnswerInfoTS(sock, &arg->dtype, &len, &arg->ndims, arg->dims, &numbytes, &dptr, &mem);
    	arg->length = len;

	printf("dtype=%d, len=%d, ndims=%d, dims[0]=%d, numbytes=%d\n",arg->dtype,arg->length,arg->ndims,arg->dims[0],numbytes);

	*((int*)out) = *((int*)dptr);
*/
	
	stat = myMdsValue(sock, expression, arg);
	*((int*)out) = *((int*)arg->ptr);

	stat = MdsClose(sock);

	stat = DisconnectFromMds(sock);
}



int myMdsValue(SOCKET sock, char *expression, struct descrip *ans_arg)
{
  int i;
  unsigned char nargs = 1;
  unsigned char idx = 0;
  int status = 1;
  struct descrip exparg;
  struct descrip *arg = &exparg;
  
  arg = MakeDescrip(&exparg,DTYPE_CSTRING,0,0,expression);

  status = SendArg(sock, idx, arg->dtype, nargs, ArgLen(arg), arg->ndims, arg->dims, arg->ptr);
  
  if (status & 1)
  {
    short len;
    int numbytes;
    void *dptr;
    void *mem = 0;
    status = GetAnswerInfoTS(sock, &ans_arg->dtype, &len, &ans_arg->ndims, ans_arg->dims, &numbytes, &dptr, &mem);
    ans_arg->length = len;
    if (numbytes)
    {
      if (ans_arg->dtype == DTYPE_CSTRING)
      {
        ans_arg->ptr = malloc(numbytes+1);
        ((char *)ans_arg->ptr)[numbytes] = 0;
      }
      else if (numbytes > 0)
        ans_arg->ptr = malloc(numbytes);
      if (numbytes > 0)
        memcpy(ans_arg->ptr,dptr,numbytes);
    }
    else
      ans_arg->ptr = NULL;
    if (mem) free(mem);
  }
  else
    ans_arg->ptr = NULL;
  return status;
}


