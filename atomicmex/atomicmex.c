/* sigv = atomicmex(w,'target','model')
 * sigma*v for ionization for different atoms and molecules
 *
 * Compile mex file with (using gnumex):
 *
 * mex atomicmex.c atomic.o ../common/common.o
 *
 * S. H. Muller, 2012/01/21
 */

#include "mex.h"

#include "atomic.h"

#define STRLEN 1024



void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
    int i;
    char target[STRLEN];
    char type[STRLEN];
    char model[STRLEN];
    
    const int n = mxGetNumberOfElements(R[0]);
    double *x = mxGetData(R[0]), *w;
    
    mxGetString(R[1],target,STRLEN);
    mxGetString(R[2],type,STRLEN);
    mxGetString(R[3],model,STRLEN);
    
    atomic_desc D = get_atomic_desc(target, type, model);
    
    L[0] = mxCreateNumericMatrix(n,1,mxDOUBLE_CLASS,mxREAL);
    w = mxGetData(L[0]);
    
    memcpy(w,x,n*sizeof(double));
    
    sigmav_vec(n, w, &D);
    
}
