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
    int n;
    func *fun;
    void *data;
    double *ab;
    int pos;
    double *par;
} intpar;

double Fun(double x, void *data)
{
    double *par = data;
    double y=par[0], z=par[1];
    return x*y*y*z*z*z;
}

double integrate(double x, void *data)
{
    intpar *IP = (intpar*)data;
    IP->par[IP->pos] = x;
    return gauss_legendre(IP->n, IP->fun, IP->data, IP->ab[0], IP->ab[1]);
}


mxArray *gauss_legendre_fun(int nI, const int *d, const int *n, int nR, mxArray **R)
{
    mxArray *res;
    double y;
    double par[2];
    
    double *ab1 = mxGetData(R[0]);
    double *ab2 = mxGetData(R[1]);
    double *ab3 = mxGetData(R[2]);
    
    intpar IP = {*n, Fun, par, ab1, 0, par};
    intpar IP2 = {*n, integrate, &IP, ab2, 1, par};
    
    y = gauss_legendre(*n, integrate, &IP2, ab3[0], ab3[1]);
    
    res = mxCreateDoubleScalar(y);
    
    return res;
}


void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
    const int nI = mxGetNumberOfElements(R[0]);
    const int *d = mxGetData(R[0]);
    const int *n = mxGetData(R[1]);
    
    L[0] = gauss_legendre_fun(nI, d, n, nR-2, (mxArray**)(R+2));
}
