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

typedef int (*comparfcn) (const void*, const void*);

int compar_int(const void *a, const void *b)
{
    return *(const int*)b - *(const int*)a;
}

int compar_float(const void *a, const void *b)
{
    float d = *(const float*)b - *(const float*)a;
    return (d<0)-(d>0);
}

int compar_double(const void *a, const void *b)
{
    double d = *(const double*)b - *(const double*)a;
    return (d<0)-(d>0);
}

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
    
    comparfcn compar;
    
    switch (mxID) {
        case mxDOUBLE_CLASS: compar = compar_double; break;
        case mxSINGLE_CLASS: compar = compar_float;  break;
        case mxINT32_CLASS:  compar = compar_int;    break;
        default:
            mexErrMsgTxt("Unsupported data type");
    }
    
    L[0] = mxCreateNumericArray(ndims,dims,mxID,mxREAL);
    void *y = mxGetData(L[0]);
    
    memcpy(y,x,N*bytes);
    
    qsort(y,N,bytes,compar);
    
    switch (mxID) {
        case mxDOUBLE_CLASS:
            *(double*)y = median_double(y,N);
        case mxSINGLE_CLASS:
            *(float*)y = median_float(y,N);
    }
    
}
