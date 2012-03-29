#include "mex.h"

#include "dblMaxw.h"

#define STRLEN 1024

void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
    int i;
    char f_r[STRLEN], nrm[STRLEN];
    mxGetString(*R++,f_r,STRLEN);

    double *vt = mxGetData(*R++);
    double *v0 = mxGetData(*R++);
    double *ut = mxGetData(*R++);
    double *u0 = mxGetData(*R++);
    
    int *IJ = mxGetData(*R++);
    int mV = 3-IJ[0], mU = 3-IJ[1], m = mV+mU;
    
    mxGetString(*R++,nrm,STRLEN);

    double **VU = (m==0) ? NULL : malloc(m*(sizeof(double*)+sizeof(int)));
    double **V = VU, **U=VU+mV;
    int *DI = VU+m;

    for(i=0; i<m; i++) {
        VU[i] = mxGetData(R[i]);
        DI[i] = mxGetNumberOfElements(R[i]);
    }

    if (m == 0) {
        L[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
    } else {
        L[0] = mxCreateNumericArray(m,DI,mxDOUBLE_CLASS,mxREAL);
    }
    double *Y = mxGetData(L[0]);

    dblMaxw(f_r, vt, v0, ut, u0, IJ, nrm, DI, V, U, Y);

    if (VU) free(VU);
}



