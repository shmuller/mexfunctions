#include "odepack.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void fun(data *D) {
    *D->ydot = -(*D->y);
}

int main() {

    double time[] = {0., 1.};
    double res[2] = {1.};

    data DATA;
    data *D = &DATA;

    D->dy_dt = fun;
    D->neq = 1;
    D->time = time;
    D->res = res;

    odeint(D);

    printf("%e, %e\n", D->res[1], exp(-D->time[1]));

    return 0;
}


