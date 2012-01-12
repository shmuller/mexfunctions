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

double gauss_legendre_matlab(int n, mxArray *fun, double a, double b)
{
    static mxArray *R[] = {NULL, NULL};
    
	double* x = NULL;
	double* w = NULL;
	double A,B,Ax,s=0.;
	int i, dtbl, m;

	m = (n+1)>>1;

    double *X = malloc(n*sizeof(double));
    double *p;
    
    if (R[1] == NULL) {
        R[1] = mxCreateNumericArray(0,NULL,mxDOUBLE_CLASS,mxREAL);
        mexMakeArrayPersistent(R[1]);
    }
    
    R[0] = fun;
    
    mwSize dims[] = {n};
    mxSetDimensions(R[1], dims, 1);
    mxSetData(R[1], X);
    mxArray *L;
    
    dtbl = gauss_legendre_load_tbl(n, &x, &w);
    
	A = 0.5*(b-a);
	B = 0.5*(b+a);

    i = 0;
    p = X;
    if(n&1) {
        *p++ = B;
        i = 1;
    }
    for(;i<m;i++) {
        Ax = A*x[i];
        *p++ = B+Ax;
        *p++ = B-Ax;
    }
    
    mexCallMATLAB(1, &L, 2, R, "feval");
    
    i = 0;
    p = mxGetData(L);
    if(n&1) {
        s = w[0]*(*p++);
        i = 1;
    }
    for(;i<m;i++) {
        s += w[i]*(p[0] + p[1]);
        p += 2;
    }
    
    mxDestroyArray(L);
    free(X);
    
	if (dtbl) {
		free(x);
		free(w);
	}
	return A*s;
}


void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{   
    register int i;
    const mwSize ndims = mxGetNumberOfDimensions(R[1]);
    const mwSize *dims = mxGetDimensions(R[1]);
    const int npts = mxGetNumberOfElements(R[1]);
    
    double *x=mxGetData(R[1]), *y;
        
    int n = 256;
        
    L[0] = mxCreateDoubleScalar(0.);
    y = mxGetData(L[0]);
    
    *y = gauss_legendre_matlab(n, (mxArray*)R[0], x[0], x[1]);
    
}

