#include "math.h"

#include "../specfunmex/specfun.h"

#define ISQRT2PI 0.39894228040143267793994605993438

real Fun(real *par)
{
    real x=par[0], y=par[1], z=par[2];
    return x*y*y*z*z*z;
}

real dim1(real *par)
{
    real x=par[0], a=par[1], r, ra;
    
    if (x != 0.) {
        r = 1./x-1; ra = r-a;
        return ISQRT2PI*exp(-0.5*ra*ra)*(1+exp(-2*a*r))/(x*x);
    } else {
        return 0.;
    }
}

real dim2(real *par)
{
    real x=par[0], a=par[1], r, ra, ar;
    
    if (x != 0.) {
        r = 1./x-1; ra = r-a; ar = a*r;
        return r*exp(-0.5*ra*ra)*besei0(&ar)/(x*x);
    } else {
        return 0.;
    }
}

static const pair functions[] = {
    "Fun", Fun,
    "dim1", dim1,
    "dim2", dim2
};
