#include <stdio.h>

#include "../quadmex/gauss_legendre.h"

#include "../common/common.h"

#include "angle_int.c"

typedef struct {
    func *fun;
    SM_REAL *par;
} DATA;


SM_REAL intgrd(SM_REAL r, void *data)
{
    DATA *D = (DATA*) data;
    *D->par = r;
    return D->fun(D->par);
}


int main()
{
    int i, j;
    SM_REAL res;

    /* r, R0, Rt, z0, zt */
    SM_REAL par[][5] = {
        {0., 0., 1. , 0., 1. },
        {0., 1., 1. , 1., 1. },
        {0., 1., 0.5, 1., 1. },
        {0., 1., 1. , 1., 0.5},
        {0., 0., 0.5, 1., 1. },
        {0., 0., 1. , 1., 0.5}
    };

    func *fun[] = {AngleInt, AngleInt2};

    int Ni = KV_LEN(par), Nj = KV_LEN(fun);

    for(i=0; i<Ni; i++) {
        for(j=0; j<Nj; j++) {
            DATA D = {fun[j], par[i]};

            res = gauss_legendre(32, intgrd, (void*) &D, 0., 10.);

            printf("%.16e ",res);
        }
        printf("\n");
    }

    return 0;
}
