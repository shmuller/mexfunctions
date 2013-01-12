#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <minpack.h>

static void test(data *D) {
    int i;
    double *P = D->P, *x = D->x, *y = D->y, *ydata = D->ydata;
    for (i=D->m; i--; ) {
        *y++ = P[0]*exp(-P[1] * *x++) - *ydata++;
    }
}

int main() {
    data D;

    D.m = 10;
    D.n = 2;
    
    void *mem = malloc((3*D.m + D.n)*sizeof(double));
    D.x = mem;
    D.y = D.x + D.m;
    D.ydata = D.y + D.m;
    D.P = D.ydata + D.m;
    
    double noise[] = {0.01, -0.03, -0.01, 0.05, 0.02, -0.01, 0.03, -0.02, 0.04, 0.01};

    int i;
    for(i=0; i < D.m; ++i) {
        D.x[i] = i;
        D.ydata[i] = exp(-D.x[i]) + noise[i];
    }
    
    D.P[0] = 0.8;
    D.P[1] = 1.2;

    leastsq1(test, &D);

    printf("P = [%f, %f]\n", D.P[0], D.P[1]);

    D.P[0] = 0.8;
    D.P[1] = 1.2;

    leastsq(test, &D);

    printf("P = [%f, %f]\n", D.P[0], D.P[1]);

    free(mem);
    return 0;
}

