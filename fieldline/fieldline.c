#include "fieldline.h"

#include <odepack.h>
#include <dierckx.h>
#include <math.h>

double twopi=2.*M_PI;
double BR_Bphi, Bz_Bphi;
double *X, *Y;
int n_bdry;

dierckx_data *splR, *splz;


void f(int *neq, double *t, double *y, double *ydot) {
    double fact = *y*twopi;
    splR->y = splz->y = y;
    splR->x = splz->x = y + 1;
    dierckx_bispev(splR);
    dierckx_bispev(splz);

    ydot[0] = fact*BR_Bphi;
    ydot[1] = fact*Bz_Bphi;
    ydot[2] = fact*sqrt(1. + BR_Bphi*BR_Bphi + Bz_Bphi*Bz_Bphi);
}


int pnpoly(double x, double y) {
    int i, j, c=0;
    for (i=0, j=n_bdry-1; i<n_bdry; j=i++) {
        if (((Y[i] > y) != (Y[j] > y)) &&
	        (x < (X[j] - X[i]) * (y - Y[i]) / (Y[j] - Y[i]) + X[i]))
            c = !c;
    }
    return c;
}


void g(int *neq, double *t, double *y, int *ng, double *gout) {
    gout[0] = (pnpoly(y[0], y[1])) ? 1.0 : -1.0;
}


void solve_bdry(fieldline_data *F) {
    splR = F->splR;
    splz = F->splz;
    splR->z = &BR_Bphi;
    splz->z = &Bz_Bphi;    

    void *memR = dierckx_make_workspace(splR);
    void *memz = dierckx_make_workspace(splz);

    n_bdry = F->n_bdry;
    X = F->bdry;
    Y = X + n_bdry;

    data *D = F->D;
    D->f = f;
    D->neq = 3;
    D->g = g;
    D->ng = 1;

    odesolve(D);

    free(memR);
    free(memz);
}

