/* quadfun.c
 *
 * S. H. Muller, 2012/01/18
 */

#include "stdlib.h"

#include "gauss_legendre.h"

#include "quadfun.h"

typedef struct {
    func *fun;
    int n;
    SM_REAL *ab;
    SM_REAL *par;
} intpar;


SM_REAL integrate(SM_REAL x, void *data)
{
    intpar *IP = data;
    *IP->par = x;
    if (IP->fun != NULL) {
        return IP->fun((SM_REAL*)(IP+1));  // note contiguous memory!
    } else {
        return gauss_legendre(IP->n, integrate, IP+1, IP->ab[0], IP->ab[1]);
    }
}

void quadfun(func *fun, link *LI, int nI, int Np, int Ns, int N, SM_REAL *y)
{
    int i, j, nR=nI+Np+Ns, n;
    link *li, *lp, *ls;
    
    // allocate contiguous memory
    intpar *ip, *IP = malloc(nI*sizeof(intpar)+nR*sizeof(SM_REAL));
    SM_REAL *ab, *par = (SM_REAL*)(IP+nI);
    
    // write singleton parameters
    for(i=0,ls=LI+nR-1; i<Ns; i++,ls--) {
        par[ls->o] = *ls->x;
    }
    
    if (nI == 0) {
        // don't integrate, just evaluate
        for(i=0; i<N; i++) {
            for(j=0,lp=LI+nI; j<Np; j++,lp++) {
                par[lp->o] = lp->x[i];
            }
            *y++ = fun(par);
        }
    } else {
        // get integration variables
        n = LI->N;
        ab = LI->x;
        for(i=1,ip=IP,li=LI; i<nI; i++,ip++,li++) {
            ip->n = (li+1)->N;
            ip->ab = (li+1)->x;
            ip->fun = NULL;
            ip->par = par + li->o;
        }
        ip->fun = fun;
        ip->par = par + li->o;
        
        // perform integrations
        for(i=0; i<N; i++) {
            for(j=0,lp=LI+nI; j<Np; j++,lp++) {
                par[lp->o] = lp->x[i];
            }
            *y++ = gauss_legendre(n, integrate, IP, ab[0], ab[1]);
        }
    }
    
    free(IP);
}

