/* y = quadmex(int32([di,dj]),int32([ni,nj]),@fun,p1,...,li,...,lj,...,pn)
 * Use gauss_legendre() to calculate integral over Matlab function fun.
 *
 * Compile mex file with (using gnumex):
 *
 * mex -v quadmex.c quad.o gauss_legendre.o
 *
 * S. H. Muller, 2012/01/12
 */

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "string.h"

#include "quadfun.h"


double Fun(double *par)
{
    double x=par[0], y=par[1], z=par[2];
    return x*y*y*z*z*z;
}


mxArray *gauss_legendre_fun(int nI, const int *d, const int *n, int nR, mxArray **R)
{
    int i, j, N, Np, Ns;
    mxArray *res;
    double *y;
    func *fun = Fun;
    
    // prepare input permutations
    link *li, *lp, *ls, *LI = malloc(nR*sizeof(link));
    
    int *isI = malloc(nR*sizeof(int));
    for(i=0; i<nR; i++) isI[i] = -1;
    for(i=0; i<nI; i++) isI[d[i]-1] = i;
    
    for(i=0,lp=LI+nI,ls=LI+nR-1; i<nR; i++) {
        N = mxGetNumberOfElements(R[i]);
        li = (isI[i] < 0) ? ((N > 1) ? lp++ : ls--) : LI+isI[i];
        li->x = mxGetData(R[i]);
        li->N = N;
        li->o = i;
    }
    Np = lp-LI-nI;
    Ns = LI+nR-1-ls;
    free(isI);
    
    
    N = (nI==nR) ? 1 : (LI+nI)->N;
    res = mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL);
    y = mxGetData(res);
    
    quadfun(fun, LI, nI, Np, Ns, N, y, n);
    
    free(LI);
    
    return res;
}


void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
    const int nI = mxGetNumberOfElements(R[0]);
    const int *d = mxGetData(R[0]);
    const int *n = mxGetData(R[1]);
    
    L[0] = gauss_legendre_fun(nI, d, n, nR-2, (mxArray**)(R+2));
}
