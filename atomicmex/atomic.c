/*
function AP = atomic_params(target)
%AP = atomic_params(target)
%   Obtain atomic physics parameters for target atom.
%
%   S. H. Muller, 2009/03/29
%
%   B ....... binding energy
%   U ....... orbital kinetic energy
%   N ....... electron occupation number
%   C ....... series expansion coefficients for df_dw (differential oscillator strength)
%   Ni ...... int df_dw dw
%   Mi2 ..... R/B int 1/(w+1) df_dw dw
%   K ....... 2-Ni/N 
%   Q ....... 2*B*Mi2/(N*R) (dipole constant)

*/

#include "../common/common.h"
#include "atomic.h"
#include "atomic_param.h"

#define ME_Q 5.685557358631881E-012
#define MP_Q 1.043939583073274E-008

/* differential oscillator strengths */
typedef void (diff_osc_fun)(double, const atomic_param *, double *, double *);

void diff_osc_BEQ(double t, const atomic_param *AP, double *D0, double *K)
{
    *D0 = AP->Q*0.5*(1.-1./(t*t));
    *K = 2.-AP->Q;
}

void diff_osc_BEB(double t, const atomic_param *AP, double *D0, double *K)
{
    *D0 = 0.5*(1.-1./(t*t));
    *K = 1.;
}

const keyval kv_diff_osc[] = {
    "BEQ", diff_osc_BEQ,
    "BEB", diff_osc_BEB
};

diff_osc_fun *get_diff_osc(const char *model)
{
    return kv_select(KV_LEN(kv_diff_osc), kv_diff_osc, model);  
}


/* total integrated cross sections */
typedef double (TICS_fun)(double, const atomic_desc *);

double TICS_ion_shell(double T, const atomic_param *AP, diff_osc_fun *diff_osc)
{
    double t=T/AP->B, u=AP->U/AP->B, A=AP->fact*AP->S/(t+(u+1)/AP->n);
    double D0, K;
    
    diff_osc(t, AP, &D0, &K);
    
    return (t < 1.) ? 0. : A*((D0-K/(t+1))*log(t)+K*(1-1./t));
}

double TICS_ion(double T, const atomic_desc *D)
{
    int i;
    double s = 0.;
    const atomic_param *AP = D->AP;
    diff_osc_fun *diff_osc = D->diff_osc;
    
    for(i=AP->shells; i--; ) {
        s += TICS_ion_shell(T, AP++, diff_osc);
    }
    return s;
}

double TICS_CX(double T, const atomic_desc *D)
{
        static const double a = 7.042E-10, b = 0.414E-10; 
        double sqrt_s = a;
        const atomic_param *AP = D->AP;
        
        T /= AP->A;
        if (T > 1.) {
            sqrt_s -= b*log(T);
        }
        return sqrt_s*sqrt_s;
}

const keyval kv_TICS[] = {
    "ion", TICS_ion,
    "CX", TICS_CX
};

TICS_fun *get_TICS(const char *type)
{
    return kv_select(KV_LEN(kv_TICS), kv_TICS, type);  
}


/* public functions */
atomic_desc get_atomic_desc(const char *target, const char *type, const char *model)
{
    atomic_desc D;
    const atomic_param *AP = get_AP(target);
    
    D.AP = AP;
    D.TICS = get_TICS(type);
    D.diff_osc = get_diff_osc(model);
    
    if (strcmp(type,"ion")==0) {
        D.mu_2q = 0.5*ME_Q;
    } else {
        D.mu_2q = 0.25*MP_Q*AP->A;
    }
    return D;
}


double sigmav(double w, const atomic_desc *D)
{
    double T = D->mu_2q*(w*w);
    TICS_fun *TICS = D->TICS;
    return w*TICS(T, D);
}


void sigmav_vec(int n, double *w, const atomic_desc *D)
{
    int i, N;
    double T, mu_2q = D->mu_2q;
    TICS_fun *TICS = D->TICS;
    
    for(i=n; i--; ) {
        T = mu_2q*(*w)*(*w);
        *w++ *= TICS(T, D);
    }
}

