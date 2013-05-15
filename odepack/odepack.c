#include <stdlib.h>

#include "odepack.h"

#define max(x, y) (((x) > (y)) ? (x) : (y))

extern dlsoda_(void *f, int *neq, double *y, double *t, double *tout, int *itol,
               double *rtol, double *atol, int *itask, int *istate, int *iopt, 
               double *rwork, int *lrw, int *iwork, int *liw, void *jac, int *jt);

data *DD;

void wrapper(int *neq, double *t, double *y, double *ydot) {
    DD->t = t;
    DD->y = y;
    DD->ydot = ydot;
    DD->dy_dt(DD);
}

int odeint(data *D) {
    int neq=D->neq, itol=1, itask=1, istate=1, iopt=0, jt=2;
    double rtol=1.49012e-8, atol=1.49012e-8;
    
    int lrw = 22 + neq*max(16, neq+9);
    int liw = 20 + neq;

    void *mem = malloc(lrw*sizeof(double) + liw*sizeof(int));
    double *rwork = mem;
    int *iwork = (int*)(rwork + lrw);

    void *jac = NULL;
    
    DD = D;

    int i;
    double *t = D->time;
    double *y = D->res;
    double tout;

    tout = t[1];
    t[1] = t[0]; ++t;
    for (i=0; i<neq; ++i,++y) y[neq] = y[0];

    dlsoda_(wrapper, &neq, y, t, &tout, &itol, &rtol, &atol, 
            &itask, &istate, &iopt, rwork, &lrw, iwork, &liw, jac, &jt);

    free(mem);
}


