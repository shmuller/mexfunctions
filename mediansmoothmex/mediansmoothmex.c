/* y = mediansmoothmex(x,w)
 *
 * Compile on Linux: 
 * mex -v COPTIMFLAGS=-O3 LDOPTIMFLAGS=-O3 mediansmoothmex.c
 * 
 * Compile on Windows
 * mex -v OPTIMFLAGS=-O3 mediansmoothmex.c
 *
 * S. H. Muller, 2010/02/22
 */

#include <stdlib.h>

#include "mex.h"
#include "matrix.h"
#include "math.h"

#define torben median_double
#define elem_type double
#include "torben.c"
#undef torben
#undef elem_type

#define torben median_float
#define elem_type float
#include "torben.c"
#undef torben
#undef elem_type

#define torben median_int
#define elem_type int
#include "torben.c"
#undef torben
#undef elem_type

typedef int (*medianfcn) (const void*, const int);

void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
    register int i;
    
    const mwSize ndims = mxGetNumberOfDimensions(R[0]);
    const mwSize *dims = mxGetDimensions(R[0]);
    const mxClassID mxID = mxGetClassID(R[0]);
    const size_t N = mxGetNumberOfElements(R[0]);
    const size_t bytes = mxGetElementSize(R[0]);
    const void *x = mxGetData(R[0]);
    
    const int *w = mxGetData(R[1]);
    
    medianfcn median;
    
    switch (mxID) {
        case mxDOUBLE_CLASS: median = median_double; break;
        case mxSINGLE_CLASS: median = median_float;  break;
        case mxINT32_CLASS:  median = median_int;    break;
        default:
            mexErrMsgTxt("Unsupported data type");
    }
    
    L[0] = mxCreateNumericArray(ndims,dims,mxINT32_CLASS,mxREAL);
    int *ind = mxGetData(L[0]);
    
    if (*w+1 >= N) {
        int t = median(x,N);
        for (i=0; i<N; i++) *ind++ = t;
    } else {
        for (i=0; i<*w; i++) {
            *ind++ = median(x, i+1+*w);
        }
        for (i=0; i<N-2*(*w); i++) {
            *ind++ = i + median(x+i*bytes, 2*(*w)+1);
        }
        for (; i<N-*w; i++) {
            *ind++ = i + median(x+i*bytes, N-i);
        }
    }
    
}
