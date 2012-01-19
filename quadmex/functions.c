#include "math.h"

#include "../specfunmex/specfun.h"


#define ISQRT2   0.70710678118654752440084436210485
#define ISQRTPI  0.56418958354775628694807945156077
#define SQRT2_PI 0.79788456080286535587989211986876
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




real _I_a_eq_1(real r, real z0)
{
    real rz = r-z0, zr = 2.*z0*r, sinhcs = (zr==0) ? 1. : (1-exp(-zr))/zr;
    
    return SQRT2_PI*r*r*exp(-0.5*rz*rz)*sinhcs;
}

real _I_a_lt_1(real r, real z0, real a2)
{
    return 0.;
}

real _I_a_gt_1(real r, real z0, real a2)
{
    return 0.;
}

real angle_int_driftZ(real *par)
{
    real r=par[0], w=par[1], z0=fabs(par[2]), a2=1./(w*w);

    if (a2 == 1.) {
        return _I_a_eq_1(r,z0);
    } else if (a2 < 1) {
        return _I_a_lt_1(r,z0,a2);
    } else {
        return _I_a_gt_1(r,z0,a2);
    }
}



static const pair functions[] = {
    "Fun", Fun,
    "dim1", dim1,
    "dim2", dim2,
    "angle_int_driftZ", angle_int_driftZ
};
