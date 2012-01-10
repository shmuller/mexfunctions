/* y = specfunmex(name, x)
 * Matlab interface for specfun.f90. Compile object file with:
 *
 * gfortran -c -O3 -fno-underscoring specfun.f90
 *
 * Compile mex file with (using gnumex with gfortran as linker):
 *
 * mex -v specfunmex.c specfun.o
 *
 * S. H. Muller, 2012/01/09
 */

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "string.h"

#define STRLEN 1024

typedef double real;

typedef real (func)(real *);

typedef struct {
    char *name;
    func *fun;
} pair;

extern func besei0;   // exponentially scaled Bessel I0(X) function
extern func besei1;   // exponentially scaled Bessel I1(X) function
extern func besek0;   // exponentially scaled Bessel K0(X) function
extern func besek1;   // exponentially scaled Bessel K1(X) function
extern func besi0;    // Bessel I0(X) function
extern func besi1;    // Bessel I1(X) function
extern func besj0;    // Bessel J0(X) function
extern func besj1;    // Bessel J1(X) function
extern func besk0;    // Bessel K0(X) function
extern func besk1;    // Bessel K1(X) function
extern func besy0;    // Bessel Y0(X) function
extern func besy1;    // Bessel Y1(X) function
extern func daw;      // Dawson's integral function
extern func dlgama;   // log ( Gamma ( X ) ) for a real argument
extern func ei;       // exponential integral Ei(X)
extern func eone;     // exponential integral E1(X)
extern func expei;    // scaled exponential integral exp(-X) * Ei(X)
extern func r8_erf;   // error function
extern func r8_erfc;  // complementary error function
extern func r8_erfcx; // exponentially scaled complementary error function
extern func r8_gamma; // Gamma(X) for a real argument
extern func r8_psi;   // Psi(X)

static const pair P[] = {
    "besei0", besei0,
    "besei1", besei1,
    "besek0", besek0,
    "besek1", besek1,
    "besi0",  besi0,
    "besi1",  besi1,
    "besj0",  besj0,
    "besj1",  besj1,
    "besk0",  besk0,
    "besk1",  besk1,
    "besy0",  besy0,
    "besy1",  besy1,
    "daw",    daw,
    "dlgama", dlgama,
    "ei",     ei,
    "eone",   eone,
    "expei",  expei,
    "erf",    r8_erf,
    "erfc",   r8_erfc,
    "erfcx",  r8_erfcx,
    "gamma",  r8_gamma,
    "psi",    r8_psi
};

void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{   
    register int i;
    char name[STRLEN];
    const mwSize ndims = mxGetNumberOfDimensions(R[1]);
    const mwSize *dims = mxGetDimensions(R[1]);
    const int npts = mxGetNumberOfElements(R[1]);
    
    double *x=mxGetData(R[1]), *y;
    func *fun = NULL;
    const pair *p;
    
    mxGetString(R[0],name,STRLEN);
    
    for(i=sizeof(P)/sizeof(pair),p=P; i--; p++) {
        if (strcmp(name,p->name)==0) {
            fun = p->fun;
            break;
        }
    }
    if (fun == NULL) {
        mexErrMsgTxt("specfunmex: Unknown function name");
    }
        
    L[0] = mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
    y = mxGetData(L[0]);
    
    for(i=npts; i--; )
        *y++ = fun(x++);
    
}

