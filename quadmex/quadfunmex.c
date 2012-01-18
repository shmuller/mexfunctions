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

#include "gauss_legendre.h"

#define max(a,b) ((b) > (a) ? (b) : (a))

typedef double (func)(double *);

typedef struct {
    func *fun;
    int n;
    double *ab;
    double *par;
} intpar;

double Fun(double *par)
{
    double x=par[0], y=par[1], z=par[2];
    return x*y*y*z*z*z;
}

double integrate(double x, void *data)
{
    intpar *IP = data;
    *IP->par = x;
    if (IP->fun != NULL) {
        return IP->fun((double*)(IP+1));
    } else {
        return gauss_legendre(IP->n, integrate, IP+1, IP->ab[0], IP->ab[1]);
    }
}

typedef struct {
    double *x;
    int N;
    int o;
} link;

mxArray *gauss_legendre_fun(int nI, const int *d, const int *n, int nR, mxArray **R)
{
    int i, j, N, Np, Ns;
    mxArray *res;
    double *y;
    
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
    
    
    N = (LI+nI)->N;
    res = mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL);
    y = mxGetData(res);
    
    
    intpar *ip, *IP = malloc(nI*sizeof(intpar)+nR*sizeof(double));
    double *ab, *par = (double*)(IP+nI);
    
    // get integration variables
    ab = LI->x;
    for(i=1,ip=IP,li=LI; i<nI; i++,ip++,li++) {
        ip->fun = NULL;
        ip->n = n[i];
        ip->ab = (li+1)->x;
        ip->par = par + li->o;
    }
    ip->fun = Fun;
    ip->par = par + li->o;
    
    
    // write singleton parameters
    for(i=0,ls=LI+nR-1; i<Ns; i++,ls--) {
        par[ls->o] = *ls->x;
    }
    
    // perform integrations
    for(i=0; i<N; i++) {
        for(j=0,lp=LI+nI; j<Np; j++,lp++) {
            par[lp->o] = lp->x[i];
        }
        *y++ = gauss_legendre(n[0], integrate, IP, ab[0], ab[1]);
    }
    
    free(IP);
    
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
