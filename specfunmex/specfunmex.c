/* y = specfunmex(nu,x,scale)
 *
 * S. H. Muller, 2012/01/09
 */

#include "mex.h"
#include "matrix.h"
#include "math.h"

#define besi0 besi0_
#define besei0 besei0_

extern double besi0(double *);
extern double besei0(double *);

void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{   
    register int i;
    const mwSize ndims = mxGetNumberOfDimensions(R[1]);
    const mwSize *dims = mxGetDimensions(R[1]);
    const int npts = mxGetNumberOfElements(R[1]);
    double *nu=mxGetData(R[0]), *x=mxGetData(R[1]), *scale=mxGetData(R[2]), *y;
    
    double (*bess)(double *) = (*scale > 0.) ? besei0 : besi0;
    
    L[0] = mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
    y = mxGetData(L[0]);
    
    for(i=npts; i--; )
        *y++ = bess(x++);
}

