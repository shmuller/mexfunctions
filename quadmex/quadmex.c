/* y = quadmex(fun, x)
 * Use gauss_legendre() to calculate integral over Matlab function fun.
 *
 * Compile mex file with (using gnumex):
 *
 * mex -v quadmex.c gauss_legendre.o
 *
 * S. H. Muller, 2012/01/12
 */

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "string.h"

#include "gauss_legendre.h"

mxArray *gauss_legendre_matlab(int n, int d, int nR, mxArray **R)
{
    static mxArray *Rd = NULL;
    mxArray *L, *Y, *AB = R[d];
    
	double* x = NULL;
	double* w = NULL;
	double A,B,Ax,s;
	int i, j, dtbl, m, M;

	m = (n+1)>>1;

    double *ab = mxGetData(AB);
    A = 0.5*(ab[1]-ab[0]);
	B = 0.5*(ab[1]+ab[0]);
    
    double *p, *q, *y, *X = malloc(n*sizeof(double));
    
    if (Rd == NULL) {
        Rd = mxCreateNumericMatrix(1,0,mxDOUBLE_CLASS,mxREAL);
        mexMakeArrayPersistent(Rd);
    }
    mxSetN(Rd, n);
    mxSetData(Rd, X);
    
    dtbl = gauss_legendre_load_tbl(n, &x, &w);

    i = 0;
    p = X; q = p+n-1;
    if(n&1) {
        *p++ = B;
        i = 1;
    }
    for(;i<m;i++) {
        Ax = A*x[i];
        *p++ = B+Ax;
        *q-- = B-Ax;
    }
    
    R[d] = Rd;
    mexCallMATLAB(1, &Y, nR, R, "feval");
    R[d] = AB;
    
    M = mxGetM(Y);
    y = mxGetData(Y);
    
    for(i=n&1,p=y+i,q=y+n-1; i<m; i++) {
        *p++ += *q--;
    }
    
    L = mxCreateNumericMatrix(M,1,mxDOUBLE_CLASS,mxREAL);
    p = mxGetData(L);
    
    *p = 0.;
    for(i=0,q=y; i<m; i++) {
        *p += w[i]*(*q++);
    }
    *p *= A;
    
    mxDestroyArray(Y);
    free(X);
    
	if (dtbl) {
		free(x);
		free(w);
	}
	return L;
}


void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{   
    register int i;
    const mwSize ndims = mxGetNumberOfDimensions(R[1]);
    const mwSize *dims = mxGetDimensions(R[1]);
    const int npts = mxGetNumberOfElements(R[1]);
    
    int n = *((int*)mxGetData(R[0]));
    int d = *((int*)mxGetData(R[1]));
    
    L[0] = gauss_legendre_matlab(n, d, nR-2, (mxArray**)(R+2));
    
}

