#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "minpack.h"

#define lmdif lmdif_

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

    printf("P = [%f, %f], [%f, %f]\n",
            D->P[0], D->P[1], x[0], x[1]);
    fflush(stdout);

    memcpy(D->P, x, (*n)*sizeof(double));
    f(D);
    memcpy(fvec, D->y, (*m)*sizeof(double));
    return 0;
}

int leastsq(func *f, data *D)
{
    double xtol = 1.49012e-8, ftol = 1.49012e-8;
    double gtol = 0.0, epsfcn = 0.0, factor = 1.0e2;
    int mode = 1, nprint = 0, info = 99, nfev, ldfjac;
    int n = D->n, m = D->m, maxfev = 200*(n+1);

    void *mem = calloc((n+n*m+n+3*n+m)*sizeof(double)+n*sizeof(int), 1);
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

    printf("info = %d\n", info);
    
    free(mem);
    return 0;
}

int leastsq1(func *f, data *D)
{
    double tol = 1.49012e-8;
    int m = D->m, n = D->m, lwa = n*m + 5*n + m, info = 99;

    void *mem = calloc(lwa*sizeof(double) + n*sizeof(int), 1);
    double *wa = mem;
    int *iwa = (int*)(wa + lwa);

    CONTAINER.f = f;
    CONTAINER.D = D;

    double *fvec = malloc(m*sizeof(double));

    double *x = malloc(n*sizeof(double));
    memcpy(x, D->P, n*sizeof(double));

    lmdif1_(wrapper, &m, &n, x, fvec, &tol, &info, iwa, wa, &lwa);

    printf("info = %d\n", info);
    
    free(x);
    free(fvec);
    free(mem);
    return 0;
}


