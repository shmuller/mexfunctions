#include "math.h"

#include "../common/common.h"
#include "../specfunmex/specfun.h"

#include "gauss_legendre.h"

#define TWOPI    6.28318530717958647692528676655901
#define I2PI     0.15915494309189533576888376337251
#define ISQRT2   0.70710678118654752440084436210485
#define ISQRTPI  0.56418958354775628694807945156077
#define SQRT2_PI 0.79788456080286535587989211986876
#define SQRT2PI  2.50662827463100050241576528481105
#define ISQRT2PI 0.39894228040143267793994605993438

real dim1(real *par)
{
    real x=par[0], a=par[1], r, ra;
    
    if (x != 0.) {
        r=1./x-1; ra=r-a;
        return ISQRT2PI*exp(-0.5*ra*ra)*(1+exp(-2*a*r))/(x*x);
    } else {
        return 0.;
    }
}

real dim2(real *par)
{
    real x=par[0], a=par[1], r, ra, ar;
    
    if (x != 0.) {
        r=1./x-1; ra=r-a; ar=a*r;
        return r*exp(-0.5*ra*ra)*besei0(&ar)/(x*x);
    } else {
        return 0.;
    }
}



/* private functions */
real _I_a_eq_1(real r, real z0)
{
    real rz=r-z0, zr=2.*z0*r, sinhcs=(zr==0) ? 1. : (1-exp(-zr))/zr;
    return SQRT2_PI*r*r*exp(-0.5*rz*rz)*sinhcs;
}

real _L_a_gt_1(real r, real z0, real b)
{
    real rz=r-z0, x=ISQRT2*(z0/b+r*b);
    return exp(-0.5*rz*rz)*daw(&x);
}

real _L_a_lt_1(real r, real z0, real b)
{
    real rz=r+z0, z0_b=z0/b, x=ISQRT2*(z0_b+r*b);
    
    if (x < 0.) {
        return exp(-0.5*(1-b*b)*(r*r-z0_b*z0_b))*r8_erfc(&x);
    } else {
        return exp(-0.5*rz*rz)*r8_erfcx(&x);
    }
}

real _I_a_gt_1(real r, real z0, real a2)
{
    real b=sqrt(a2-1);
    return ISQRTPI*(a2/b)*r*(_L_a_gt_1(r,z0,b) - _L_a_gt_1(-r,z0,b));
}

real _I_a_lt_1(real r, real z0, real a2)
{
    real b=sqrt(1-a2);
    return 0.5*(a2/b)*r*(_L_a_lt_1(-r,z0,b) - _L_a_lt_1(r,z0,b));
}

real _cos_th_int(real c, void *data)
{
    real *par = data;
    real r=par[0], ar=par[1], z0=par[2], aR0=par[3];
    real s=sqrt(1-c*c), aR=s*ar, daR=aR-aR0, dz=c*r-z0, x=aR0*aR;
    
    return exp(-0.5*(daR*daR+dz*dz))*besei0(&x);
}


/* public functions */
real angle_int_iso(real *par)
{
    return _I_a_eq_1(par[0], par[1]);
}

real angle_int_Z(real *par)
{
    real r=par[0], a=par[1], z0=fabs(par[2]), a2=a*a;

    if (a2 == 1.) {
        return _I_a_eq_1(r,z0);
    } else if (a2 > 1) {
        return _I_a_gt_1(r,z0,a2);
    } else {
        return _I_a_lt_1(r,z0,a2);
    }
}

real angle_int_RZ(real *par)
{
    real r=par[0], a=par[1], z0=par[2], R0=par[3], ar=a*r, aR0=a*R0;
    real par2[] = {r,ar,z0,aR0};
    
    return ISQRT2PI*ar*ar*gauss_legendre(32, _cos_th_int, par2, -1., 1.);
}

real angle_int(real *par)
{
    if (par[1] == 1.) {
        return _I_a_eq_1(par[0], hypot(par[2],par[3]));
    } else if (par[3] == 0) {
        return angle_int_Z(par);
    } else {
        return angle_int_RZ(par);
    }
}


/* tests */
real Maxw(real *par)
{
    real x=par[0], y=par[1], z=par[2];
    real a=par[3], b=par[4], c=par[5], v=par[6];
    real r2=x*x+y*y+z*z, f, M;
    
    x -= a; y -= b; z -= c;
    
    M = exp(-0.5*(x*x+y*y+z*z)/(v*v));
    f = SQRT2PI*v;
    M /= f*f*f;
    return sqrt(r2)*M;
}

real Maxw_r(real *par)
{
    real r=par[0];
    real a=par[1], b=par[2], c=par[3], v=par[4];
    real R0=hypot(a,b), z0=c;
    
    real par2[] = {r/v, 1., z0/v, R0/v};
    
    return r*angle_int(par2)/v;
}


real dblMaxw(real *par)
{
    real x=par[ 0], y=par[ 1], z=par[ 2];
    real a=par[ 3], b=par[ 4], c=par[ 5], v=par[ 6];
   
    real X=par[ 7], Y=par[ 8], Z=par[ 9];
    real A=par[10], B=par[11], C=par[12], V=par[13];
    
    real dx, dy, dz, vrel, f, M;
    
    dx = x-X; dy = y-Y; dz = z-Z;
    vrel = sqrt(dx*dx+dy*dy+dz*dz);
    
    x -= a; y -= b; z -= c;
    X -= A; Y -= B; Z -= C;
    
    M = exp(-0.5*((x*x+y*y+z*z)/(v*v)+(X*X+Y*Y+Z*Z)/(V*V)));
    f = TWOPI*v*V;
    M /= f*f*f;
    
    return vrel*M;
}

static const keyval functions[] = {
    "dim1", dim1,
    "dim2", dim2,
    "angle_int_iso", angle_int_iso,
    "angle_int_Z", angle_int_Z,
    "angle_int_RZ", angle_int_RZ,
    "angle_int", angle_int,
    "Maxw", Maxw,
    "Maxw_r", Maxw_r,
    "dblMaxw", dblMaxw
};

