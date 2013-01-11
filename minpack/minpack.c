#include <stdlib.h>
#include "minpack.h"

#define chkder chkder_
#define hybrd hybrd_
#define hybrj hybrj_
#define lmdif lmdif_
#define lmder lmder_
#define lmstr lmstr_

extern void chkder_(int*,int*,double*,double*,double*,int*,double*,double*,int*,double*);
extern void hybrd_(void*,int*,double*,double*,double*,int*,int*,int*,double*,double*,int*,double*,int*,int*,int*,double*,int*,double*,int*,double*,double*,double*,double*,double*);
extern void hybrj_(void*,int*,double*,double*,double*,int*,double*,int*,double*,int*,double*,int*,int*,int*,int*,double*,int*,double*,double*,double*,double*,double*);
extern void lmdif_(void*,int*,int*,double*,double*,double*,double*,double*,int*,double*,double*,int*,double*,int*,int*,int*,double*,int*,int*,double*,double*,double*,double*,double*);
extern void lmder_(void*,int*,int*,double*,double*,double*,int*,double*,double*,double*,int*,double*,int*,double*,int*,int*,int*,int*,int*,double*,double*,double*,double*,double*);
extern void lmstr_(void*,int*,int*,double*,double*,double*,int*,double*,double*,double*,int*,double*,int*,double*,int*,int*,int*,int*,int*,double*,double*,double*,double*,double*);


typedef struct {
    func *f;
    data *D;
} container;

container CONTAINER;

static int wrapper(int *m, int *n, double *x, double *fvec, int *iflag)
{
    func *f = CONTAINER.f;
    data *D = CONTAINER.D;
    f(D);
    return 0;
}

int leastsq(func *f, data *D)
{
    double xtol = 1.49012e-8, ftol = 1.49012e-8;
    double gtol = 0.0, epsfcn = 0.0, factor = 1.0e2;
    int mode = 1, nprint = 0, info, nfev, ldfjac;
    int n = D->n, m = D->m, maxfev = 200*(n+1);

    void *mem = malloc((n+n*m+n+3*n+m)*sizeof(double)+n*sizeof(int));
    double *diag = mem;
    double *fjac = diag + n;
    double *qtf = fjac + n*m;
    double *wa1 = qtf + n;
    double *wa2 = wa1 + n;
    double *wa3 = wa2 + n;
    double *wa4 = wa3 + n;
    int *ipvt = (int*)(wa4 + m);
    
    CONTAINER.f = f;
    CONTAINER.D = D;

    lmdif(wrapper, &m, &n, D->P, D->y, &ftol, &xtol, &gtol, &maxfev, &epsfcn, 
          diag, &mode, &factor, &nprint, &info, &nfev, fjac, &ldfjac, ipvt, qtf, 
          wa1, wa2, wa3, wa4);
    
    free(mem);
    return 0;
}


