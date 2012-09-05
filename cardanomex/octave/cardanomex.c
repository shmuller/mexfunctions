/* x = cardanomex(coefs)
 * Real roots of cubic polynomial a*x^3+b*x^2+c*x+d=0 using Cardano's
 * formula. The coefficients are given as coefs = [a,b,c,d].'.
 *
 * Compile on Linux: 
 * mex -v COPTIMFLAGS=-O3 LDOPTIMFLAGS=-O3 cardanomex.c
 * 
 * Compile on Windows
 * mex -v OPTIMFLAGS=-O3 cardanomex.c
 *
 * S. H. Muller, 2011/07/31
 */

#include "mex.h"
#include "../cardano.h"

void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
    register int i; 
    const int npts = mxGetNumberOfElements(R[0])/4;
    const mwSize dims[] = {3,npts};
    
    const double *coefs = mxGetPr(R[0]);
    
    L[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
    double *x = mxGetPr(L[0]);
   
    for (i=npts; i--; x+=3,coefs+=4) {
        cardano(coefs,x);
    }
}

