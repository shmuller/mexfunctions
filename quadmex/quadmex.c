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

double gauss_legendre_matlab(int n, int d, int nR, mxArray **R)
{
    static mxArray *Rd = NULL;
    mxArray *L, *AB = R[d];
    
	double* x = NULL;
	double* w = NULL;
	double A,B,Ax,s=0.;
	int i, dtbl, m;

	m = (n+1)>>1;

    double *ab = mxGetData(AB);
    A = 0.5*(ab[1]-ab[0]);
	B = 0.5*(ab[1]+ab[0]);
    
    double *p, *X = malloc(n*sizeof(double));
    
    if (Rd == NULL) {
        Rd = mxCreateNumericArray(0,NULL,mxDOUBLE_CLASS,mxREAL);
        mexMakeArrayPersistent(Rd);
    }
    
    mwSize dims[] = {n};
    mxSetDimensions(Rd, dims, 1);
    mxSetData(Rd, X);
    
    dtbl = gauss_legendre_load_tbl(n, &x, &w);

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
    
    R[d] = Rd;
    mexCallMATLAB(1, &L, nR, R, "feval");
    R[d] = AB;
    
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
    
    double *y;
    
    int n = *((int*)mxGetData(R[0]));
    int d = *((int*)mxGetData(R[1]));    
    
    L[0] = mxCreateDoubleScalar(0.);
    y = mxGetData(L[0]);
    
    *y = gauss_legendre_matlab(n, d, nR-2, (mxArray**)(R+2));
    
}

