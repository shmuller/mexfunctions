#include "odepack.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


void f(int *neq, double *t, double *y, double *ydot) {
    ydot[0] = -y[0];
    ydot[1] = -y[1];
}

/*
int term(data *D) {
    printf("t = %.10e, y[0] = %.10e, y[1] = %.10e\n", D->t[0], D->y[0], D->y[1]);
    return *D->t > 3.;
}
*/

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
    //D->term = term;

    odesolve(D);

    for (i=0; i<D->points_done; i++) {
        printf("%e, %e, %e\n", y[2*i], y[2*i+1], exp(-t[i]));
    }

    return 0;
}


