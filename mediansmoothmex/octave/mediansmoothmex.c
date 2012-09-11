/* y = mediansmoothmex(double(x),int32(w))
 *
 * Compile on Linux: 
 * mex -v COPTIMFLAGS=-O3 LDOPTIMFLAGS=-O3 mediansmoothmex.c
 * 
 * Compile on Windows
 * mex -v OPTIMFLAGS=-O3 mediansmoothmex.c
 *
 * S. H. Muller, 2010/02/22
 */


#include <mex.h>

#include "../mediansmooth.h"

void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
    const mwSize ndims = mxGetNumberOfDimensions(R[0]);
    const mwSize *dims = mxGetDimensions(R[0]);
    const mxClassID mxID = mxGetClassID(R[0]);
    const int N = mxGetNumberOfElements(R[0]);
    const int bytes = mxGetElementSize(R[0]);
    const void *x = mxGetData(R[0]);
    
    const int *w = mxGetData(R[1]);
    
    L[0] = mxCreateNumericArray(ndims,dims,mxINT32_CLASS,mxREAL);
    int *ind = mxGetData(L[0]);
    
    switch (mxID) {
        case mxDOUBLE_CLASS:
            median_filt_double(x,N,*w,ind);
            break;
        case mxSINGLE_CLASS:
            median_filt_float(x,N,*w,ind);
            break;
        default:
            break;
    }
}
