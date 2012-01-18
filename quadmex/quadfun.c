/* y = quadmex(int32([di,dj]),int32([ni,nj]),@fun,p1,...,li,...,lj,...,pn)
 * Use gauss_legendre() to calculate integral over Matlab function fun.
 *
 * Compile mex file with (using gnumex):
 *
 * mex -v quadmex.c quad.o gauss_legendre.o
 *
 * S. H. Muller, 2012/01/12
 */

#include "string.h"

#include "gauss_legendre.h"

#include "quadfun.h"

typedef struct {
    func *fun;
    int n;
    double *ab;
    double *par;
} intpar;


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

void quadfun(func *fun, link *LI, int nI, int Np, int Ns, int N, double *y, const int *n)
{
    int i, j, nR=nI+Np+Ns;
    link *li, *lp, *ls;
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
    ip->fun = fun;
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
}

