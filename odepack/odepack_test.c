#include "odepack.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void fun(data *D) {
    D->ydot[0] = -D->y[0];
    D->ydot[1] = -D->y[1];
}

#define SIZE(x) (sizeof(x)/sizeof(x[0]))

int main() {

    int i;
    double time[] = {0., 1., 2., 3.};
    double res[SIZE(time)*2] = {1., 1.};

    data DATA = {0};
    data *D = &DATA;

    D->dy_dt = fun;
    D->neq = 2;
    D->n = SIZE(time);
    D->time = time;
    D->res = res;

    odesolve(D);

    for (i=0; i<D->n; i++) {
        printf("%e, %e, %e\n", D->res[2*i], D->res[2*i+1], exp(-D->time[i]));
    }

    return 0;
}


