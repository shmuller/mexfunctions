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
    mxArray *L, *Y, *Rd = R[d];
    
	double* x = NULL;
	double* w = NULL;
	double A,B,Ax,wi;
	int i, j, dtbl, o, m, M, M2, Mo, Mm, nn, mm;

    o = n&1;
	m = (n+1)>>1;

    double *ab = mxGetData(Rd);
    A = 0.5*(ab[1]-ab[0]);
	B = 0.5*(ab[1]+ab[0]);
    
    double *p, *q, *y, *X = malloc(n*sizeof(double));
    
    // load abscissae and weight table
    dtbl = gauss_legendre_load_tbl(n, &x, &w);

    // populate X with scaled abscissae
    p = X+m; q = p-1;
    if(o) {
        *q-- = B;
    }
    for(i=o; i<m; i++) {
        Ax = A*x[i];
        *p++ = B+Ax;
        *q-- = B-Ax;
    }
    
    // attach abscissa vector X to Rd
    mm = mxGetM(Rd);
    nn = mxGetN(Rd);
    mxSetM(Rd, 1);
    mxSetN(Rd, n);
    mxSetData(Rd, X);
    
    // evaluate function at abscissa vector X
    mexCallMATLAB(1, &Y, nR, R, "feval");
    
    // re-attach ab to Rd
    mxSetM(Rd, mm);
    mxSetN(Rd, nn);
    mxSetData(Rd, ab);
    
    M = mxGetM(Y);
    y = mxGetData(Y);
    
    M2 = M*2;
    Mo = M*o;
    Mm = M*m;
    
    // add symmetric pairs
    for(i=o,p=y+Mm,q=p-M-Mo; i<m; i++,q-=M2) {
        for(j=0; j<M; j++) {
            *p++ += *q++;
        }
    }
    
    // create output matrix
    L = mxCreateNumericMatrix(M,1,mxDOUBLE_CLASS,mxREAL);
    p = mxGetData(L);
    
    // populate output matrix with weighted sums
    for(i=0,q=y+Mm-Mo; i<m; i++,p-=M) {
        for(j=0,wi=w[i]; j<M; j++) {
            *p++ += wi*(*q++);
        }
    }
    // apply final scaling
    for(j=0; j<M; j++) {
        *p++ *= A;
    }
    
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

