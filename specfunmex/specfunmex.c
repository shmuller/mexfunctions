/* y = specfunmex(name, x)
 * Matlab interface for specfun.f90. Compile object file with:
 *
 * gfortran -c -O3 -fno-underscoring specfun.f90
 *
 * Compile mex file with (using gnumex with gfortran as linker):
 *
 * mex -v specfunmex.c specfun.o ../common/common.o
 *
 * S. H. Muller, 2012/01/09
 */

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "string.h"

#include "../common/common.h"
#include "specfun.h"

#define STRLEN 1024

void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{   
    register int i;
    char name[STRLEN];
    const mwSize ndims = mxGetNumberOfDimensions(R[1]);
    const mwSize *dims = mxGetDimensions(R[1]);
    const int npts = mxGetNumberOfElements(R[1]);
    
    double *x=mxGetData(R[1]), *y;
    func *fun;
    
    mxGetString(R[0],name,STRLEN);
    
    fun = kv_select(KV_LEN(kv_specfun), kv_specfun, name);
    if (fun == NULL) {
        mexErrMsgTxt("specfunmex: Unknown function name");
    }
    
    L[0] = mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
    y = mxGetData(L[0]);
    
    for(i=npts; i--; )
        *y++ = fun(x++);
    
}

