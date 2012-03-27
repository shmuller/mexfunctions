#include "mex.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "../common/common.h"
#include "../atomicmex/atomic.h"

#include "../quadmex/gauss_legendre.h"

#include "angle_int.c"

#define STRLEN 1024

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


typedef struct {
    SM_REAL val;
    unsigned idx;
} PAIR;


int compare(const void *A, const void *B)
{
    SM_REAL a = ((PAIR*)A)->val;
    SM_REAL b = ((PAIR*)B)->val;
    return (a > b) - (a < b);
}

#define R0 data_angle[1]
#define Rt data_angle[2]
#define z0 data_angle[3]
#define zt data_angle[4]

#define Nr 32

void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
    int i, j, k, l;

    char f_r[STRLEN];
    mxGetString(*R++,f_r,STRLEN);

    double *vt = mxGetData(*R++);
    double *v0 = mxGetData(*R++);
    double *ut = mxGetData(*R++);
    double *u0 = mxGetData(*R++);
    
    double wt = hypot(*vt,*ut);
    double rM = 10.*wt;

    int *IJ = mxGetData(*R++);
    int m = 6-IJ[0]-IJ[1];

    double w0[] = {u0[0]-v0[0], u0[1]-v0[1], u0[2]-v0[2]};
    double data_angle[5], *y;
    
    intgrd_data ID = {0., 1., NULL, NULL, AngleInt3_, data_angle};

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
            L[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
            y = mxGetData(L[0]);

            R0 = 0.;
            Rt = wt;
            z0 = sqrt(w0[0]*w0[0]+w0[1]*w0[1]+w0[2]*w0[2]);
            zt = wt;
                
            *y = gauss_legendre(Nr, intgrd, &ID, 0., rM);
        
        } else if (IJ[0] == 2) {
               
            double *v3 = mxGetData(R[0]);
            int ndims = mxGetNumberOfDimensions(R[0]);
            mwSize *dims = mxGetDimensions(R[0]);
            int N = mxGetNumberOfElements(R[0]);

            L[0] = mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
            y = mxGetData(L[0]);

            R0 = hypot(w0[0],w0[1]);
            Rt = wt;
            zt = *ut;

            for(i=0; i<N; i++) {
                z0 = u0[2]-v3[i];
                *y++ = gauss_legendre(Nr, intgrd, &ID, 0., rM);
                //*y++ = fM(v3[i]-v0[2],*vt)*gauss_legendre(Nr, intgrd, &ID, 0., rM);
            }

        } else if (IJ[0] == 1) {
            double *v1 = mxGetData(R[0]);
            double *v2 = mxGetData(R[1]);
            int N1 = mxGetNumberOfElements(R[0]);
            int N2 = mxGetNumberOfElements(R[1]);

            mwSize dims[] = {N1,N2};
            L[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
            y = mxGetData(L[0]);

            double dx, dy;

            Rt = *ut;
            z0 = w0[2];
            zt = wt;

            for(j=0; j<N2; j++) {
                dy = u0[1]-v2[j];
                for(i=0; i<N1; i++) {
                    dx = u0[0]-v1[i];
                    R0 = hypot(dx,dy);
                    *y++ = gauss_legendre(Nr, intgrd, &ID, 0., rM);
                }
            }
        }

    } else if (IJ[1] == 2 && IJ[0] == 2) {
        double *v3 = mxGetData(R[0]);
        double *u3 = mxGetData(R[1]);
        int N1 = mxGetNumberOfElements(R[0]);
        int N2 = mxGetNumberOfElements(R[1]);

        mwSize dims[] = {N1,N2};
        L[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
        y = mxGetData(L[0]);

        ID.f_angle = AngleInt2;
        R0 = hypot(w0[0],w0[1]);
        Rt = wt;

        for(j=0; j<N2; j++) {
            for(i=0; i<N1; i++) {
                ID.x0 = u3[j]-v3[i];
                *y++ = gauss_legendre(Nr, intgrd, &ID, 0., rM);
            }
        }

    } else if (IJ[1] == 1 && IJ[0] == 1) {
        double *v1 = mxGetData(R[0]);
        double *u1 = mxGetData(R[1]);
        double *v2 = mxGetData(R[2]);
        double *u2 = mxGetData(R[3]);
        int N1 = mxGetNumberOfElements(R[0]);
        int N2 = mxGetNumberOfElements(R[1]);
        int N3 = mxGetNumberOfElements(R[2]);
        int N4 = mxGetNumberOfElements(R[3]);
        int N = N1*N2*N3*N4;

        mwSize dims[] = {N1,N2,N3,N4};
        L[0] = mxCreateNumericArray(4,dims,mxDOUBLE_CLASS,mxREAL);
        y = mxGetData(L[0]);

        PAIR *pair = malloc(N*sizeof(PAIR));
        PAIR *p = pair;
        unsigned idx = 0;

        ID.f_angle = AngleInt1;
        R0 = w0[2];
        Rt = wt;

        double w1,w2;
        for(l=0; l<N4; l++) {
            for(k=0; k<N3; k++) {
                w2 = u2[l]-v2[k];
                for(j=0; j<N2; j++) {
                    for(i=0; i<N1; i++) {
                        w1 = u1[j]-v1[i];
                        ID.x0 = hypot(w1,w2);

                        p->val = ID.x0;
                        p->idx = idx++;
                        p++;
                        
                        //*y++ = gauss_legendre(Nr, intgrd, &ID, 0., rM);
                    }
                }
            }
        }

        qsort(pair,N,sizeof(PAIR),compare);

        idx = 0;
        double val;
        for(i=0,p=pair,ID.x0=-1.; i<N; i++,p++) {
            if (ID.x0 != p->val) {
                ID.x0 = p->val;
                val = gauss_legendre(Nr, intgrd, &ID, 0., rM);
                idx++;
            }
            y[p->idx] = val;
        }

        printf("idx = %lu\n", idx);
    }
    
}

