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

typedef double (func)(double, void*);

typedef struct {
    func *fun;
    int n;
    double *ab;
    double *par;
} intpar;

double Fun(double z, void *data)
{
    double *par = data;
    double x=par[0], y=par[1];
    return x*y*y*z*z*z;
}

double integrate(double x, void *data)
{
    intpar *IP = data;
    *IP->par = x;
    return gauss_legendre(IP->n, IP->fun, IP+1, IP->ab[0], IP->ab[1]);
}


mxArray *gauss_legendre_fun(int nI, const int *d, const int *n, int nR, mxArray **R)
{
    int i, j, m = nR-1;
    mxArray *res;
    double y;
    
    intpar *ip, *IP = malloc(m*(sizeof(intpar)+sizeof(double)));
    double *ab, *par = (double*)(IP+m);
    
    for(i=1,ip=IP; i<=m; i++,ip++) {
        ip->fun = (i < m) ? integrate : Fun;
        ip->n = n[i];
        ip->ab = mxGetData(R[i]);
        ip->par = par++;
    }
    
    ab = mxGetData(R[0]);
    y = gauss_legendre(n[0], integrate, IP, ab[0], ab[1]);
    
    /*
    double *ab1 = mxGetData(R[0]);
    double *ab2 = mxGetData(R[1]);
    double *ab3 = mxGetData(R[2]);
    
    intpar IP1 = {Fun, par, *n, ab1, par};
    intpar IP2 = {integrate, &IP1, *n, ab2, par+1};
    y = gauss_legendre(*n, integrate, &IP2, ab3[0], ab3[1]);
    */
    
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
