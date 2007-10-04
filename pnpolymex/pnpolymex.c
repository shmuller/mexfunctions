/* in = pnpolymex(x,y,xp,yp)
 * Finds points included in closed polygon. The result for points on the
 * polygon is random.
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
#include "math.h"

#define min(a,b) ((a)<(b) ? (a) : (b))
#define max(a,b) ((a)>(b) ? (a) : (b))

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
    const double *px=xp, *py=yp;
    double xm=INFINITY, xM=-INFINITY, ym=INFINITY, yM=-INFINITY; 
    
    char *in;
    
    L[0] = mxCreateNumericArray(ndims,dims,mxLOGICAL_CLASS,mxREAL);
    in = (char*) mxGetPr(L[0]);
    
    /* calculate bounding box */
    for (i=npol; i--; px++,py++) {
        xm=min(xm,*px); xM=max(xM,*px);
        ym=min(ym,*py); yM=max(yM,*py);
    }

    /* test each point first against bounding box then against polygon */
    for (i=npts; i--; x++,y++) {
        *in++ = (xm<=*x && *x<=xM && ym<=*y && *y<=yM) && pnpoly(npol,xp,yp,*x,*y);
    }
}
