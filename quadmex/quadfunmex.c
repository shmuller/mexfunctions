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


mxArray *gauss_legendre_fun(int nI, const int *d, const int *n, int nR, mxArray **R)
{
    int i, j, l, L=1, c, C, nP = nR-nI;
    mxArray *res;
    double *y;
    double **P, **p;
    
    intpar *ip, *IP = malloc(nI*sizeof(intpar)+nR*sizeof(double));
    double *ab, *par = (double*)(IP+nI);
    
    // get parameters
    if (nP > 0) {
        P = malloc(2*nP*sizeof(double*));
        p = P+nP;
        for(j=0,C=0; j<nR; j++) {
            for(i=0; i<nI; i++) if (j==d[i]-1) break;
            if (i == nI) {
                P[C] = mxGetData(R[j]);
                l = mxGetNumberOfElements(R[j]);
                if (l == 1) {
                    par[j] = *P[C];  // write singleton parameters
                } else {
                    L = max(L,l);
                    p[C] = par+j;
                    C++;
                }
            }
        }
    }
    
    res = mxCreateNumericMatrix(L,1,mxDOUBLE_CLASS,mxREAL);
    y = mxGetData(res);
    
    // get integration variables
    for(i=1,ip=IP; i<nI; i++,ip++) {
        ip->fun = NULL;
        ip->n = n[i];
        ip->ab = mxGetData(R[d[i]-1]);
        ip->par = par+d[i-1]-1;
    }
    ip->fun = Fun;
    ip->par = par+d[i-1]-1;
    ab = mxGetData(R[d[0]-1]);
    
    // perform integrations
    for(i=0; i<L; i++) {
        for(c=0; c<C; c++) {
            *p[c] = P[c][i];
        }
        *y++ = gauss_legendre(n[0], integrate, IP, ab[0], ab[1]);
    }
    
    if (nP > 0) free(P);
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
