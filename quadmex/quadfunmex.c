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

typedef double (func)(void*);

typedef struct {
    func *fun;
    int n;
    double *ab;
    double *par;
} intpar;

double Fun(void *data)
{
    double *par = data;
    double x=par[0], y=par[1], z=par[2];
    return x*y*y*z*z*z;
}

double integrate(double x, void *data)
{
    intpar *IP = data;
    *IP->par = x;
    if (IP->fun != NULL) {
        return IP->fun(IP+1);
    } else {
        return gauss_legendre(IP->n, integrate, IP+1, IP->ab[0], IP->ab[1]);
    }
}


mxArray *gauss_legendre_fun(int nI, const int *d, const int *n, int nR, mxArray **R)
{
    int i;
    mxArray *res;
    double y;
        
    intpar *ip, *IP = malloc(nI*sizeof(intpar)+nR*sizeof(double));
    double *ab, *par = (double*)(IP+nI);
    
    for(i=1,ip=IP; i<nI; i++,ip++) {
        ip->fun = NULL;
        ip->n = n[i];
        ip->ab = mxGetData(R[i]);
        ip->par = par+i-1;
    }
    ip->fun = Fun;
    ip->par = par+i-1;
    
    ab = mxGetData(R[0]);
    y = gauss_legendre(n[0], integrate, IP, ab[0], ab[1]);
    
    res = mxCreateDoubleScalar(y);
    
    free(IP);
    
    return res;
}


void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
    const int nI = mxGetNumberOfElements(R[0]);
    const int *d = mxGetData(R[0]);
    const int *n = mxGetData(R[1]);
    
    L[0] = gauss_legendre_fun(nI, d, n, nR-2, (mxArray**)(R+2));
}
