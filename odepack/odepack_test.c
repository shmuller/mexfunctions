#include "odepack.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


void f(int *neq, double *t, double *y, double *ydot) {
    ydot[0] = -y[0];
    ydot[1] = -y[1];
}

void g(int *neq, double *t, double *y, int *ng, double *gout) {
    gout[0] = *t - 3.;
}

#define SIZE(x) (sizeof(x)/sizeof(x[0]))

int main() {

    int i;
    double t[] = {0., 1., 2., 3., 4., 5., 6.};
    double y[SIZE(t)*2] = {1., 1.};

    data DATA = {0};
    data *D = &DATA;

    D->f = f;
    D->neq = 2;
    D->n = SIZE(t);
    D->t = t;
    D->y = y;

    D->ng = 1;
    D->g = g;

    odesolve(D);

    for (i=0; i<D->points_done; i++) {
        printf("%e, %e, %e\n", y[2*i], y[2*i+1], exp(-t[i]));
    }

    return 0;
}


