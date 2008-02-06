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
	int stat,shot=0,zero=0,size;
	int dtype_long = DTYPE_LONG;
  	int idesc = descr2(&dtype_long, &zero);
    
	register int i;
	const mwSize ndims = mxGetNumberOfDimensions(R[0]);
	const mwSize *dims = mxGetDimensions(R[0]);
	const int npts = mxGetNumberOfElements(R[0]);
	register double *x=mxGetPr(R[0]), *y=mxGetPr(R[1]);
    
	const int npol = mxGetNumberOfElements(R[2]);
	const double *xp=mxGetPr(R[2]), *yp=mxGetPr(R[3]);
	const double *px=xp, *py=yp;
    
	int *out;
    
	L[0] = mxCreateNumericArray(ndims,dims,mxINT32_CLASS,mxREAL);
	out = (int*) mxGetPr(L[0]);
    
	sock = MdsConnect("plaspc03.ucsd.edu");
	if (sock == INVALID_SOCKET) {
		mxErrMsgTxt("Could not connect to server");
	}

	stat = MdsOpen("csdx",&shot);

	stat = MdsValue2("\\NXPIX",&idesc,out,&zero);
	
	stat = MdsClose("csdx",&shot);

	stat = DisconnectFromMds(sock);


}
