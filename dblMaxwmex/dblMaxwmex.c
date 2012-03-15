#include "mex.h"
#include "math.h"
#include "string.h"

#include "../common/common.h"

#include "../quadmex/gauss_legendre.h"

#include "../quadmex/functions.c"

#include "../atomicmex/atomic.h"

#define STRLEN 1024

#define SQRT2PI 2.50662827463100050241576528481105

typedef real (func_r)(real, void *);

typedef struct {
    real x0;
    real r0;
    func_r *f_r;
    void *data_r;
    func *f_angle;
    real *data_angle;
} intgrd_data;

real one(real r, void *data)
{
    return 1.;
}

real vrel(real r, void *data)
{
    return r;
}

real Maxwellian(real x, real x0, real xt)
{
    x -= x0;
    x /= xt;
    return exp(-0.5*x*x)/(SQRT2PI*xt);
}

real intgrd(real x, void *data)
{
    intgrd_data *ID = data;
    real r = hypot(x,ID->x0)*ID->r0;
    *ID->data_angle = x;
    return ID->f_r(r,ID->data_r)*ID->f_angle(ID->data_angle);
}


void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
    int i;
    
    double *vt = mxGetData(R[0]);
    double *v0 = mxGetData(R[1]);
    double *ut = mxGetData(R[2]);
    double *u0 = mxGetData(R[3]);
    
    double wt = hypot(*vt,*ut);
    
    int *IJ = mxGetData(R[4]);
    int m = 6-IJ[0]-IJ[1];
    
    double w0[] = {u0[0]-v0[0], u0[1]-v0[1], u0[2]-v0[2]};
    double R0, z0, angle_data[4], *y;
    
    intgrd_data ID = {0., 1., NULL, NULL, NULL, angle_data};
    
    ID.f_r = vrel;
    
    //atomic_desc D = get_atomic_desc("D", "CX", "");
    //ID.f_r = (func_r*) sigmav;
    //ID.data_r = &D;
    
    if (IJ[0] == 3 && IJ[1] == 3) {
        R0 = hypot(w0[0],w0[1]); 
        z0 = w0[2];
        ID.r0 = wt;
        ID.f_angle = angle_int;
        angle_data[1] = 1.;
        angle_data[2] = z0/wt;
        angle_data[3] = R0/wt;
        
        L[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        y = mxGetData(L[0]);
        
        *y = gauss_legendre(256, intgrd, &ID, 0., 10.);
        
    } else if (IJ[0] == 2 && IJ[1] == 3) {
        
    }
    
}
