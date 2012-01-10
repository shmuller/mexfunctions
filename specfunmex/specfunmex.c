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

typedef struct {
    char *name;
    real (*fun)(real *);
} pair;

extern real besei0(real *);   // exponentially scaled Bessel I0(X) function
extern real besei1(real *);   // exponentially scaled Bessel I1(X) function
extern real besek0(real *);   // exponentially scaled Bessel K0(X) function
extern real besek1(real *);   // exponentially scaled Bessel K1(X) function
extern real besi0(real *);    // Bessel I0(X) function
extern real besi1(real *);    // Bessel I1(X) function
extern real besj0(real *);    // Bessel J0(X) function
extern real besj1(real *);    // Bessel J1(X) function
extern real besk0(real *);    // Bessel K0(X) function
extern real besk1(real *);    // Bessel K1(X) function
extern real besy0(real *);    // Bessel Y0(X) function
extern real besy1(real *);    // Bessel Y1(X) function
extern real daw(real *);      // Dawson's integral function
extern real dlgama(real *);   // log ( Gamma ( X ) ) for a real argument
extern real ei(real *);       // exponential integral Ei(X)
extern real eone(real *);     // exponential integral E1(X)
extern real expei(real *);    // scaled exponential integral exp(-X) * Ei(X)
extern real r8_erf(real *);   // error function
extern real r8_erfc(real *);  // complementary error function
extern real r8_erfcx(real *); // exponentially scaled complementary error function
extern real r8_gamma(real *); // Gamma(X) for a real argument
extern real r8_psi(real *);   // Psi(X)

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
    real (*f)(real *) = NULL;
    const pair *p;
    
    mxGetString(R[0],name,STRLEN);
    
    for(i=sizeof(P)/sizeof(pair),p=P; i--; p++) {
        if (strcmp(name,p->name)==0) {
            f = p->fun;
            break;
        }
    }
    if (f == NULL) {
        mexErrMsgTxt("specfunmex: Unknown function name");
    }
        
    L[0] = mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
    y = mxGetData(L[0]);
    
    for(i=npts; i--; )
        *y++ = f(x++);
    
}

