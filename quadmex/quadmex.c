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

typedef double real;

typedef real (fun1)(real);
typedef real (fun2)(real, void *);

real wrap(real x, void *fun)
{
    static mxArray *R[2] = {NULL,NULL};
    static double *p = NULL;
    double y;
    mxArray *L;
    
    if (p == NULL) {
        R[1] = mxCreateDoubleScalar(0.);
        mexMakeArrayPersistent(R[1]);
        p = mxGetData(R[1]);
    }
    
    R[0] = (mxArray*) fun;
    *p = x;
            
    mexCallMATLAB(1, &L, 2, R, "feval");
    y = *((double*)mxGetData(L));
    
    mxDestroyArray(L);
    return y;
    
    //return (*(fun1*)fun)(x);
}

void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{   
    register int i;
    const mwSize ndims = mxGetNumberOfDimensions(R[1]);
    const mwSize *dims = mxGetDimensions(R[1]);
    const int npts = mxGetNumberOfElements(R[1]);
    
    double *x=mxGetData(R[1]), *y;
    
    int n = 32;
    
    //L[0] = mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
    L[0] = mxCreateDoubleScalar(0.);
    y = mxGetData(L[0]);
    
    *y = gauss_legendre(n, wrap, (void*)R[0], x[0], x[1]);
    
    //*y = wrap(x[0], (void*)R[0]);
    
    //mexCallMATLAB(1, L, 2, R, "feval");
    
}

