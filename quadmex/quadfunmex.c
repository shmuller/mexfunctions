/* y = quadfunmex(int32([di,dj]),int32([ni,nj]),'fun',p1,...,li,...,lj,...,pn)
 * Use gauss_legendre() to calculate integral over function 'fun' in
 * functions.c
 *
 * Compile mex file with (using gnumex):
 *
 * mex quadfunmex.c quadfun.o gauss_legendre.o ../specfunmex/specfun.o
 *
 * S. H. Muller, 2012/01/18
 */

#include "mex.h"
#include "math.h"
#include "string.h"

#include "../common/common.h"

#include "quadfun.h"

#include "functions.c"

#define STRLEN 1024

link *mk_link(int nI, const int *d, const int *n, int nR, mxArray **R, 
    int *Np, int *Ns)
{
    int i, N;
    link *li, *lp, *ls, *LI = malloc(nR*sizeof(link));
    
    int *isI = malloc(nR*sizeof(int));
    for(i=0; i<nR; i++) isI[i] = -1;
    for(i=0; i<nI; i++) isI[d[i]-1] = i;
    
    for(i=0,lp=LI+nI,ls=LI+nR-1; i<nR; i++) {
        N = mxGetNumberOfElements(R[i]);
        if (isI[i] < 0) {
            li = (N > 1) ? lp++ : ls--;
        } else {
            li = LI+isI[i];
            N = n[isI[i]];
        }
        li->x = mxGetData(R[i]);
        li->N = N;
        li->o = i;
    }
    *Np = lp-LI-nI;
    *Ns = LI+nR-1-ls;
    
    free(isI);
    return LI;
}

mxArray *gauss_legendre_fun(func *fun, int nI, const int *d, const int *n, 
    int nR, mxArray **R)
{
    int i, j, N, Np, Ns;
    mxArray *res;
    double *y;
    
    // prepare input permutations
    link *LI = mk_link(nI, d, n, nR, R, &Np, &Ns);
    
    N = (nI==nR) ? 1 : (LI+nI)->N;
    res = mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL);
    y = mxGetData(res);
    
    quadfun(fun, LI, nI, Np, Ns, N, y);
    
    free(LI);
    
    return res;
}


void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
    int i;
    char name[STRLEN];
    const int nI = mxGetNumberOfElements(R[0]);
    const int *d = mxGetData(R[0]);
    const int *n = mxGetData(R[1]);

    mxGetString(R[2],name,STRLEN);

    func *fun = kv_select(KV_LEN(functions), functions, name);
    if (fun == NULL) {
        mexErrMsgTxt("quadfunmex: Unknown function name");
    }
     
    L[0] = gauss_legendre_fun(fun, nI, d, n, nR-3, (mxArray**)(R+3));
}
