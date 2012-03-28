#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "dblMaxw.h"

#include "../common/common.h"
#include "../atomicmex/atomic.h"

#include "../quadmex/gauss_legendre.h"

#include "angle_int.h"


#define SQRT2PI 2.50662827463100050241576528481105

typedef SM_REAL (func_r)(SM_REAL, void *);

typedef struct {
    SM_REAL x0;
    SM_REAL r0;
    func_r *f_r;
    void *data_r;
    func *f_angle;
    SM_REAL *data_angle;
} intgrd_data;

SM_REAL one(SM_REAL r, void *data)
{
    return 1.;
}

SM_REAL vrel(SM_REAL r, void *data)
{
    return r;
}

SM_REAL fM(SM_REAL x, SM_REAL xt)
{
    x /= xt;
    return exp(-0.5*x*x)/(SQRT2PI*xt);
}

SM_REAL intgrd(SM_REAL x, void *data)
{
    intgrd_data *ID = data;
    SM_REAL r = hypot(x,ID->x0)*ID->r0;
    *ID->data_angle = x;
    return ID->f_r(r,ID->data_r)*ID->f_angle(ID->data_angle);
}

typedef int (compare_t)(void *, const void *, const void *);

int compare_r(SM_REAL *DATA, const unsigned *A, const unsigned *B)
{
    SM_REAL a = DATA[*A];
    SM_REAL b = DATA[*B];

    return (a < b) ? -1 : (a > b);
}


#define R0 data_angle[1]
#define Rt data_angle[2]
#define z0 data_angle[3]
#define zt data_angle[4]

#define Nr 32

void dblMaxw(char *f_r, double *vt, double *v0, double *ut, double *u0,
    int *IJ, int *DI, double **V, double **U, double *Y)
{
    int i, j, k, l;

    double *y;

    double wt = hypot(*vt,*ut);
    double rM = 10.*wt;

    double w0[] = {u0[0]-v0[0], u0[1]-v0[1], u0[2]-v0[2]};

    double data_angle[5];
    intgrd_data ID = {0., 1., NULL, NULL, AngleInt3_, data_angle};

    int m = 6-IJ[0]-IJ[1], N=1;
    for(i=0; i<m; i++) N *= DI[i];

    atomic_desc D;

    if (strcmp(f_r,"ion")==0) {
        D = get_atomic_desc("D", "ion", "BEB");
        ID.f_r = (func_r*) sigmav;
        ID.data_r = &D;
    } else if (strcmp(f_r,"CX")==0) {
        D = get_atomic_desc("D", "CX", "");
        ID.f_r = (func_r*) sigmav;
        ID.data_r = &D;
    } else if (strcmp(f_r,"vrel")==0) {
        ID.f_r = vrel;
    } else {
        ID.f_r = one;
    }

    if (IJ[1] == 3) {
        if (IJ[0] == 3) {
            R0 = 0.;
            Rt = wt;
            z0 = sqrt(w0[0]*w0[0]+w0[1]*w0[1]+w0[2]*w0[2]);
            zt = wt;
                
            *Y = gauss_legendre(Nr, intgrd, &ID, 0., rM);
        
        } else if (IJ[0] == 2) {   
            double *v3 = V[0];

            R0 = hypot(w0[0],w0[1]);
            Rt = wt;
            zt = *ut;

            for(i=0,y=Y; i<N; i++) {
                z0 = u0[2]-v3[i];
                *y++ = gauss_legendre(Nr, intgrd, &ID, 0., rM);
            }

            for(i=0,y=Y; i<N; i++) {
                *y++ *= fM(v3[i]-v0[2],*vt);
            }

        } else if (IJ[0] == 1) {
            double *v1 = V[0], *v2 = V[1];
            double dx, dy, fM1, fM2;

            Rt = *ut;
            z0 = w0[2];
            zt = wt;

            for(j=0,y=Y; j<DI[1]; j++) {
                dy = u0[1]-v2[j];
                for(i=0; i<DI[0]; i++) {
                    dx = u0[0]-v1[i];
                    R0 = hypot(dx,dy);
                    *y++ = gauss_legendre(Nr, intgrd, &ID, 0., rM);
                }
            }

            for(j=0,y=Y; j<DI[1]; j++) {
                fM2 = fM(v2[j]-v0[1],*vt);
                for(i=0; i<DI[0]; i++) {
                    fM1 = fM2*fM(v1[i]-v0[0],*vt);
                    *y++ *= fM1;
                }
            }

        }

    } else if (IJ[1] == 2 && IJ[0] == 2) {
        double *v3 = V[0], *u3 = U[0];
        double fM1, fM2;
        
        ID.f_angle = AngleInt2;
        R0 = hypot(w0[0],w0[1]);
        Rt = wt;

        for(j=0,y=Y; j<DI[1]; j++) for(i=0; i<DI[0]; i++) {
            ID.x0 = u3[j]-v3[i];
            *y++ = gauss_legendre(Nr, intgrd, &ID, 0., rM);
        }

        for(j=0,y=Y; j<DI[1]; j++) {
            fM2 = fM(u3[j]-u0[2],*ut);
            for(i=0; i<DI[0]; i++) {
                fM1 = fM2*fM(v3[i]-v0[2],*vt);
                *y++ *= fM1;
            }
        }

    } else if (IJ[1] == 1 && IJ[0] == 1) {
        double *v1 = V[0], *v2 = V[1], *u1 = U[0], *u2 = U[1];
        double val, *ys, fM1, fM2, fM3, fM4;

        unsigned *seq = malloc(N*sizeof(unsigned));
        unsigned idx, *s;

        ID.f_angle = AngleInt1;
        R0 = w0[2];
        Rt = wt;

        for(l=0,y=Y,s=seq,idx=0; l<DI[3]; l++) for(k=0; k<DI[2]; k++) {
            for(j=0; j<DI[1]; j++) for(i=0; i<DI[0]; i++) {
                *y++ = hypot(u1[k]-v1[i],u2[l]-v2[j]);
                *s++ = idx++;
            }
        }

        /*
        for(i=0,y=Y; i<N; i++) {
            ID.x0 = *y;
            *y++ = gauss_legendre(Nr, intgrd, &ID, 0., rM);
        }
        */

        qsort_r(seq, N, sizeof(unsigned), Y, compare_r);

        for(i=N,s=seq,ID.x0=-1.; i--; ) {
            ys = Y + *s++;
            if (ID.x0 != *ys) {
                ID.x0 = *ys;
                val = gauss_legendre(Nr, intgrd, &ID, 0., rM);
            }
            *ys = val;
        }

        for(l=0,y=Y; l<DI[3]; l++) {
            fM4 = fM(u2[l]-u0[1],*ut);
            for(k=0; k<DI[2]; k++) {
                fM3 = fM4*fM(u1[k]-u0[0],*ut);
                for(j=0; j<DI[1]; j++) {
                    fM2 = fM3*fM(v2[j]-v0[1],*vt);
                    for(i=0; i<DI[0]; i++) {
                        fM1 = fM2*fM(v1[i]-v0[0],*vt);
                        *y++ *= fM1;
                    }
                }
            }
        }

        free(seq);
    }
    
}

