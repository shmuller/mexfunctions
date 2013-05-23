#include <stdlib.h>

#include "dierckx.h"

void bispev_(double *tx, int *nx, double *ty, int *ny, double *c,
        int *kx, int *ky, double *x, int *mx, double *y, int *my, double *z,
        double *wrk, int *lwrk, int *iwrk, int *kwrk, int *ier);


void *dierckx_make_workspace(dierckx_data *D) {
    D->lwrk = D->mx*(D->kx+1)+D->my*(D->ky+1);
    D->kwrk = D->mx+D->my;

    void *mem = malloc(D->lwrk*sizeof(double) + D->kwrk*sizeof(int));
    D->wrk = mem;
    D->iwrk = (int*)(D->wrk + D->lwrk);
    return mem;
}

void dierckx_bispev(dierckx_data *D) {
    bispev_(D->tx, &D->nx, D->ty, &D->ny, D->c,
        &D->kx, &D->ky, D->x, &D->mx, D->y, &D->my, D->z,
        D->wrk, &D->lwrk, D->iwrk, &D->kwrk, &D->ier);
}

