#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "string.h"

#include "../common/common.h"

#include "../quadmex/gauss_legendre.h"

#include "../quadmex/functions.c"

#include "../atomicmex/atomic.h"

#define STRLEN 1024

typedef real (func_r)(real, void *);

typedef struct {
    func_r *f_r;
    void *data_r;
    real *data_angle;
} fun_data;

real one(real x, void *data)
{
    return 1.;
}

real vrel(real x, void *data)
{
    double *v = data;
    return x*(*v);
}

real intgrd(real x, void *data)
{
    fun_data *f = data;
    *f->data_angle = x;
    return f->f_r(x,f->data_r)*angle_int(f->data_angle);
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
    int n = 6-IJ[0]-IJ[1];
    
    double w0[3];
    for(i=0; i<3; i++) w0[i]=u0[i]-v0[i];
    
    double R0=hypot(w0[0],w0[1]), z0=w0[2];
    
    double angle_data[] = {0., 1., z0/wt, R0/wt};
 
    //atomic_desc D = get_atomic_desc("D", "CX", "");
    
    //fun_data f = {(func_r*) sigmav, &D, x};
    
    fun_data f = {vrel, &wt, angle_data};
    
    double res = gauss_legendre(256, intgrd, &f, 0., 10.);
    
    L[0] = mxCreateDoubleScalar(res);
    
}
