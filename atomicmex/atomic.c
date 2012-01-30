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

#define MQ_E 5.685557358631881E-012
#define MQ_P 1.043939583073274E-008


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
    return select(LENGTH(kv_diff_osc), kv_diff_osc, model);  
}


double TICS_ion_shell(double T, const atomic_param *AP, diff_osc_fun *diff_osc)
{
    double t=T/AP->B, u=AP->U/AP->B, A=AP->fact*AP->S/(t+(u+1)/AP->n);
    double D0, K;
    
    diff_osc(t, AP, &D0, &K);
    
    return (t < 1.) ? 0. : A*((D0-K/(t+1))*log(t)+K*(1-1./t));
}

double TICS_ion(double T, int shells, const atomic_param *AP, diff_osc_fun *diff_osc)
{
    int i;
    double s = 0.;
    
    for(i=0; i<shells; i++) {
        s += TICS_ion_shell(T, AP++, diff_osc);
    }
    return s;
}

void sigmav_ion(int n, double *w, const char *target, const char *model)
{
    int i, N;
    double T;
    const atomic_param *AP = get_AP(target);
    diff_osc_fun *diff_osc = get_diff_osc(model);
    
    for(i=0; i<n; i++) {
        T = 0.5*MQ_E*(*w)*(*w);
        *w++ *= TICS_ion(T, AP->shells, AP, diff_osc);
    }
}

