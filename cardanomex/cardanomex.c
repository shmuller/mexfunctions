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
#include "math.h"

#define SQRT3 1.73205080756887729

int cardano(const double *coefs, double *x)
{
    const double a3 = coefs[1]/(3.*coefs[0]);
    const double B = (coefs[2]-a3*coefs[1])/coefs[0];
    const double C = coefs[3]/coefs[0]-a3*(B+a3*a3);
    const double B3 = -B/3.;
    const double C2 = -C/2.;
    const double D2 = C2*C2-B3*B3*B3;
    double D, S, s, phi3, s_cos, s_sin;
    
    if (coefs[1] == 0. && coefs[2] == 0. && coefs[3] == 0.) {
        x[0] = x[1] = x[2] = 0.;
    } else if (D2 >= 0) {
        /* 1 real root */
        D = sqrt(D2);
        S = (C2 > 0) ? C2+D : C2-D;
        s = (S > 0) ? pow(S,1./3.) : -pow(-S,1./3.);  /* faster than s = cbrt(S); */
        x[0] = s + B3/s - a3;
        x[1] = x[2] = NAN;
    } else {
        /* 3 real roots */
        s = sqrt(B3);
        phi3 = acos(C2/(s*s*s))/3.;
        s_cos = s*cos(phi3);
        s_sin = SQRT3*s*sin(phi3);
        x[0] = 2.*s_cos - a3;
        x[1] = x[2] = -s_cos - a3;
        x[1] += s_sin;
        x[2] -= s_sin;
    }
    return 0;
}

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
