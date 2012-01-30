#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "string.h"

#include "../common/common.h"

#include "../quadmex/gauss_legendre.h"

#include "../quadmex/functions.c"

#define STRLEN 1024

typedef struct {
    func *fun;
    real *data;
} fun_data;

real one(real *x)
{
    return 1.;
}

real intgrd(real x, void *data)
{
    fun_data *f = data;
    *f->data = x;
    return f->fun(&x)*angle_int(f->data);
}


void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
    int i;
    double *x = mxGetData(R[0]); 
    
    fun_data f = {one, x};
    
    //double res = intgrd(x[0], &f);
    
    double res = gauss_legendre(32, intgrd, &f, 0., 20.);
    
    L[0] = mxCreateDoubleScalar(res);
    
}
