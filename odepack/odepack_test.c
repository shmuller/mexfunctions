#include "odepack.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define max(x, y) (((x) > (y)) ? (x) : (y))

extern dlsoda_(void *f, int *neq, double *y, double *t, double *tout, int *itol,
               double *rtol, double *atol, int *itask, int *istate, int *iopt, 
               double *rwork, int *lrw, int *iwork, int *liw, void *jac, int *jt);

void fun(int *neq, double *t, double *y, double *ydot) {
    *ydot = -(*y);
}

int main() {
    int neq=1, itol=1, itask=1, istate=1, iopt=0, jt=2;
    double t=0., tout=1.;
    double y=1.;
    double rtol=1.49012e-8, atol=1.49012e-8;
    
    int lrw = 22 + neq*max(16, neq+9);
    int liw = 20 + neq;

    void *mem = malloc(lrw*sizeof(double) + liw*sizeof(int));
    double *rwork = mem;
    int *iwork = (int*)(rwork + lrw);

    void *jac = NULL;

    dlsoda_(fun, &neq, &y, &t, &tout, &itol, &rtol, &atol, &itask, &istate, 
            &iopt, rwork, &lrw, iwork, &liw, jac, &jt);

    printf("%e, %e\n", y, exp(-tout));

    free(mem);

    return 0;
}


