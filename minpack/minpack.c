#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "minpack.h"

extern void lmdif_(void *fcn, int *m, int *n, double *x, double *fvec,
    double *ftol, double *xtol, double *gtol, int *maxfev, double *epsfcn,
    double *diag, int *mode, double *factor, int *nprint, int *info, int *nfev,
    double *fjac, int *ldfjac, int *ipvt, double *qtf, 
    double *wa1, double *wa2, double *wa3, double *wa4);

extern void lmdif1_(void *fcn, int *m, int *n, double *x, double *fvec,
    double *tol, int *info, int *iwa, double *wa, int *lwa);

typedef struct {
    func *f;
    data *D;
} container;

container CONTAINER;

static int wrapper(int *m, int *n, double *x, double *fvec, int *iflag)
{
    func *f = CONTAINER.f;
    data *D = CONTAINER.D;

    D->P = x;
    D->y = fvec;
    f(D);
    return 0;
}

int leastsq(func *f, data *D)
{
    int m = D->m, n = D->n;
    double *x = D->P, *fvec = D->y;
    
    double ftol = 1.49012e-8, xtol = 1.49012e-8, gtol = 0.0, epsfcn = 0.0, factor = 1.0e2;
    int maxfev = 200*(n+1), mode = 1, nprint = 0, info = 99, nfev = -1, ldfjac = m;

    void *mem = malloc((m*n + 5*n + m)*sizeof(double) + n*sizeof(int));
    double *diag = mem;
    double *fjac = diag + n;
    double *qtf = fjac + m*n;
    double *wa1 = qtf + n;
    double *wa2 = wa1 + n;
    double *wa3 = wa2 + n;
    double *wa4 = wa3 + n;
    int *ipvt = (int*)(wa4 + m);
    
    CONTAINER.f = f;
    CONTAINER.D = D;

    lmdif_(wrapper, &m, &n, x, fvec, &ftol, &xtol, &gtol, &maxfev, &epsfcn, 
            diag, &mode, &factor, &nprint, &info, &nfev, fjac, &ldfjac, ipvt, qtf, 
            wa1, wa2, wa3, wa4);

    D->P = x;
    D->y = fvec;
    
    free(mem);
    return 0;
}

int leastsq1(func *f, data *D)
{
    double tol = 1.49012e-8;
    int m = D->m, n = D->n, lwa = m*n + 5*n + m, info = 99;

    double *x = D->P, *fvec = D->y;

    void *mem = malloc(lwa*sizeof(double) + n*sizeof(int));
    double *wa = mem;
    int *iwa = (int*)(wa + lwa);

    CONTAINER.f = f;
    CONTAINER.D = D;

    lmdif1_(wrapper, &m, &n, x, fvec, &tol, &info, iwa, wa, &lwa);

    D->P = x;
    D->y = fvec;
   
    free(mem);
    return 0;
}


