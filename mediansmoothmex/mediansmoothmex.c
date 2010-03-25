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
#include <matrix.h>

#include <float.h>

#define TYPE double
#define TYPE_MIN (-DBL_MAX)
#define TYPE_MAX DBL_MAX
#define elem_type elem_type_double
#define elem_ptr elem_ptr_double
#define print print_double
#define compar compar_double
#define median_init median_init_double
#define find_spot find_spot_double
#define median_add median_add_double
#define median_remove median_remove_double
#define median_replace median_replace_double
#define median_filt median_filt_double

#include "mediansmooth.c"

#define TYPE float
#define TYPE_MIN (-FLT_MAX)
#define TYPE_MAX FLT_MAX
#define elem_type elem_type_float
#define elem_ptr elem_ptr_float
#define print print_float
#define compar compar_float
#define median_init median_init_float
#define find_spot find_spot_float
#define median_add median_add_float
#define median_remove median_remove_float
#define median_replace median_replace_float
#define median_filt median_filt_float

#include "mediansmooth.c"

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
