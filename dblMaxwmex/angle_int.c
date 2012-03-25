#include "math.h"

#include "../common/common.h"
#include "../specfunmex/specfun.h"

#include "../quadmex/gauss_legendre.h"

#define TWOPI    6.28318530717958647692528676655901
#define I2PI     0.15915494309189533576888376337251
#define ISQRT2   0.70710678118654752440084436210485
#define ISQRTPI  0.56418958354775628694807945156077
#define SQRT2_PI 0.79788456080286535587989211986876
#define SQRT2PI  2.50662827463100050241576528481105
#define ISQRT2PI 0.39894228040143267793994605993438


/* private functions */
SM_REAL _I_a_eq_1(SM_REAL r, SM_REAL z0)
{
    SM_REAL rz=r-z0, zr=2.*z0*r, sinhcs=(zr==0) ? 1. : (1-exp(-zr))/zr;
    return SQRT2_PI*r*r*exp(-0.5*rz*rz)*sinhcs;
}

SM_REAL _L_a_gt_1(SM_REAL r, SM_REAL z0, SM_REAL b)
{
    SM_REAL rz=r-z0, x=ISQRT2*(z0/b+r*b);
    return exp(-0.5*rz*rz)*daw(&x);
}

SM_REAL _L_a_lt_1(SM_REAL r, SM_REAL z0, SM_REAL b)
{
    SM_REAL rz=r+z0, z0_b=z0/b, x=ISQRT2*(z0_b+r*b);
    
    if (x < 0.) {
        return exp(-0.5*(1-b*b)*(r*r-z0_b*z0_b))*r8_erfc(&x);
    } else {
        return exp(-0.5*rz*rz)*r8_erfcx(&x);
    }
}

SM_REAL _I_a_gt_1(SM_REAL r, SM_REAL z0, SM_REAL a2)
{
    SM_REAL b=sqrt(a2-1);
    return ISQRTPI*(a2/b)*r*(_L_a_gt_1(r,z0,b) - _L_a_gt_1(-r,z0,b));
}

SM_REAL _I_a_lt_1(SM_REAL r, SM_REAL z0, SM_REAL a2)
{
    SM_REAL b=sqrt(1-a2);
    return 0.5*(a2/b)*r*(_L_a_lt_1(-r,z0,b) - _L_a_lt_1(r,z0,b));
}

SM_REAL _cos_th_int(SM_REAL c, void *data)
{
    SM_REAL *par = data;
    SM_REAL r=par[0], ar=par[1], z0=par[2], aR0=par[3];
    SM_REAL s=sqrt(1-c*c), aR=s*ar, daR=aR-aR0, dz=c*r-z0, x=aR0*aR;
    
    return exp(-0.5*(daR*daR+dz*dz))*besei0(&x);
}


/* public functions */
SM_REAL angle_int_iso(SM_REAL *par)
{
    return _I_a_eq_1(par[0], par[1]);
}

SM_REAL angle_int_Z(SM_REAL *par)
{
    SM_REAL r=par[0], a=par[1], z0=fabs(par[2]), a2=a*a;

    if (a2 == 1.) {
        return _I_a_eq_1(r,z0);
    } else if (a2 > 1) {
        return _I_a_gt_1(r,z0,a2);
    } else {
        return _I_a_lt_1(r,z0,a2);
    }
}

SM_REAL angle_int_RZ(SM_REAL *par)
{
    SM_REAL r=par[0], a=par[1], z0=par[2], R0=par[3], ar=a*r, aR0=a*R0;
    SM_REAL par2[] = {r,ar,z0,aR0};
    
    return ISQRT2PI*ar*ar*gauss_legendre(32, _cos_th_int, par2, -1., 1.);
}

/* r, alpha, z0, R0 */
SM_REAL angle_int(SM_REAL *par)
{
    if (par[1] == 1.) {
        return _I_a_eq_1(par[0], hypot(par[2],par[3]));
    } else if (par[3] == 0) {
        return angle_int_Z(par);
    } else {
        return angle_int_RZ(par);
    }
}

/* r, R0, Rt, z0, zt */
SM_REAL AngleInt(SM_REAL *par)
{
    SM_REAL r=par[0], R0=par[1], Rt=par[2], z0=par[3], zt=par[4];
    SM_REAL f=1./zt;
    SM_REAL par2[] = {f*r,zt/Rt,f*z0,f*R0};

    return f*angle_int(par2);
}


/**********************************************************/

SM_REAL _L_zt_ne_1(SM_REAL r, SM_REAL z0, SM_REAL izt, SM_REAL b)
{
    SM_REAL t, rz=(r+z0)*izt, z0_b=z0/b;
    if (izt < 1.) {
        t = ISQRT2*(z0_b-r*b)*izt;
        return 2.*ISQRTPI*exp(-0.5*rz*rz)*daw(&t);
    } else {
        t = ISQRT2*(z0_b+r*b)*izt;
        if (t > 0.) {
            return exp(-0.5*rz*rz)*r8_erfcx(&t);
        } else {
            return exp(-0.5*(r*r-z0_b*z0_b))*r8_erfc(&t);
        }
    }
}

SM_REAL angle_int_Z2(SM_REAL *par)
{
    SM_REAL r=par[0], izt=par[1], z0=fabs(par[2]), b;

    if (izt == 1.) {
        return _I_a_eq_1(r,z0);
    } else {
        b = sqrt(fabs(1.-1./(izt*izt)));
        return 0.5*r/b*(_L_zt_ne_1(-r,z0,izt,b)-_L_zt_ne_1(r,z0,izt,b));
    }
}


SM_REAL _cos_th_int2(SM_REAL c, void *data)
{
    SM_REAL *par = data;
    SM_REAL r=par[0], izt=par[1], z0=par[2], R0=par[3];
    SM_REAL s=sqrt(1-c*c), R=s*r, dR=R-R0, dz=(c*r-z0)*izt, x=R*R0;
    
    return exp(-0.5*(dR*dR+dz*dz))*besei0(&x);
}

SM_REAL angle_int_RZ2(SM_REAL *par)
{
    SM_REAL r=par[0], izt=par[1];
    return ISQRT2PI*r*r*izt*gauss_legendre(32, _cos_th_int2, par, -1., 1.);
}


/* r, izt, z0, R0 */
SM_REAL angle_int2(SM_REAL *par)
{
    if (par[1] == 1.) {
        return _I_a_eq_1(par[0], hypot(par[2],par[3]));
    } else if (par[3] == 0.) {
        return angle_int_Z2(par);
    } else {
        return angle_int_RZ2(par);
    }
}

/* r, R0, Rt, z0, zt */
SM_REAL AngleInt2(SM_REAL *par)
{
    SM_REAL r=par[0], R0=par[1], Rt=par[2], z0=par[3], zt=par[4];
    SM_REAL f=1./Rt;
    SM_REAL par2[] = {f*r,1./(f*zt),f*z0,f*R0};

    return f*angle_int2(par2);
}



