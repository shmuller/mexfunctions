#include <stdlib.h>

#include "odepack.h"

#define max(x, y) (((x) > (y)) ? (x) : (y))

extern dlsoda_(void *f, int *neq, double *y, double *t, double *tout, int *itol,
               double *rtol, double *atol, int *itask, int *istate, int *iopt, 
               double *rwork, int *lrw, int *iwork, int *liw, void *jac, int *jt);

extern dlsodar_(void *f, int *neq, double *y, double *t, double *tout, int *itol,
                double *rtol, double *atol, int *itask, int *istate, int *iopt, 
                double *rwork, int *lrw, int *iwork, int *liw, void *jac, int *jt,
                void *g, int *ng, int *jroot);

data *DD;

int odesolve(data *D) {
    int neq=D->neq, itol=1, itask=1, istate=1, iopt=0, jt=2, ng=D->ng;
    double rtol=1.49012e-8, atol=1.49012e-8;
    
    int lrw = 22 + neq*max(16, neq+9) + 3*ng;
    int liw = 20 + neq;

    void *mem = malloc(lrw*sizeof(double) + liw*sizeof(int));
    double *rwork = mem;
    int *iwork = (int*)(rwork + lrw);

    void *jac = NULL;

    int i, j;
    double *t = D->t;
    double *y = D->y;
    double t0;

    odefunf_t *f = D->f;
    odefun_t *term = D->term;

    DD = D;

    // setup initial conditions
    t0 = *t++;
    for (j=neq; j--; ++y) y[neq] = y[0];
    
    if (term) {
        for (i=D->n; --i; y+=neq,++t) {
            dlsoda_(f, &neq, y, &t0, t, &itol, &rtol, &atol, 
                    &itask, &istate, &iopt, rwork, &lrw, iwork, &liw, jac, &jt);

            if (term(D)) {
                --i;
                break;
            }
        }
    } else {
        for (i=D->n; --i; y+=neq,++t) {
            dlsoda_(f, &neq, y, &t0, t, &itol, &rtol, &atol, 
                    &itask, &istate, &iopt, rwork, &lrw, iwork, &liw, jac, &jt);
        }
    }
    D->points_done = D->n - i;
    free(mem);
}


