#include "odepack.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int fun(data *D) {
    D->ydot[0] = -D->y[0];
    D->ydot[1] = -D->y[1];
    return 0;
}

int term(data *D) {
    printf("t = %.10e, y[0] = %.10e, y[1] = %.10e\n", D->t[0], D->y[0], D->y[1]);
    return *D->t > 3.;
}

#define SIZE(x) (sizeof(x)/sizeof(x[0]))

int main() {

    int i;
    double time[] = {0., 1., 2., 3., 4., 5., 6.};
    double res[SIZE(time)*2] = {1., 1.};

    data DATA = {0};
    data *D = &DATA;

    D->dy_dt = fun;
    D->neq = 2;
    D->n = SIZE(time);
    D->time = time;
    D->res = res;
    D->term = term;

    odesolve(D);

    for (i=0; i<D->points_done; i++) {
        printf("%e, %e, %e\n", D->res[2*i], D->res[2*i+1], exp(-D->time[i]));
    }

    return 0;
}


