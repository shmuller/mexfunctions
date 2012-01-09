/* y = bessmex(nu,x,scale)
 *
 * S. H. Muller, 2012/01/09
 */

#include "mex.h"
#include "matrix.h"
#include "math.h"

typedef float real;

real bessi0(real x)
{
	real ax,ans;
	double y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
			+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
	} else {
		y=3.75/ax;
		ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2))))))));
	}
	return ans;
}

real bessi0_scaled(real x)
{
	real ax,ans;
	double y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=(1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
			+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2))))))*exp(-ax);
	} else {
		y=3.75/ax;
		ans=(0.39894228+y*(0.1328592e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2))))))))/sqrt(ax);
	}
	return ans;
}


void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{   
    register int i;
    const mwSize ndims = mxGetNumberOfDimensions(R[1]);
    const mwSize *dims = mxGetDimensions(R[1]);
    const int npts = mxGetNumberOfElements(R[1]);
    double *nu=mxGetData(R[0]), *x=mxGetData(R[1]), *scale=mxGetData(R[2]), *y;
    
    real (*bess)(real) = (*scale > 0.) ? bessi0_scaled : bessi0;
    
    L[0] = mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
    y = mxGetData(L[0]);
    
    for(i=npts; i--; )
        *y++ = bess(*x++);
}

