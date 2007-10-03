/* pnpolymex.c
 *
 * Algorithm: 
 * http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
 *
 * Compile on Linux: 
 * mex -v COPTIMFLAGS=-O3 LDOPTIMFLAGS=-O3 pnpolymex.c
 * 
 * Compile on Windows
 * mex -v OPTIMFLAGS=-O3 pnpolymex.c
 *
 * S. H. Muller, 2007/10/03
 */

#include "mex.h"
#include "matrix.h"

char pnpoly(int npol, const double *xp, const double *yp, double x, double y)
{
    char c = 0;
    int i, j;
    for (i = 0, j = npol-1; i < npol; j = i++) {
        if ((((yp[i]<=y) && (y<yp[j])) ||
            ((yp[j]<=y) && (y<yp[i]))) &&
            (x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))
        c = !c;
    }
    return c;
}

void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
    register int i;
    
    const mwSize ndims = mxGetNumberOfDimensions(R[0]);
    const mwSize *dims = mxGetDimensions(R[0]);
    const int npts = mxGetNumberOfElements(R[0]);
    register double *x=mxGetPr(R[0]), *y=mxGetPr(R[1]);
    
    const int npol = mxGetNumberOfElements(R[2]);
    const double *xp=mxGetPr(R[2]), *yp=mxGetPr(R[3]);
    
    char *in;
    
    L[0] = mxCreateNumericArray(ndims,dims,mxLOGICAL_CLASS,mxREAL);
    in = (char*) mxGetPr(L[0]);
    
    for (i=npts; i--; ) {
        *in++ = pnpoly(npol,xp,yp,*x++,*y++);
    }
}
