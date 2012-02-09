/* Interface definition for specfun.f90.
 *
 * S. H. Muller, 2012/01/12
 */

/*
#ifndef __COMMON__
#define __COMMON__
typedef double real;

typedef real (func)(real*);

typedef struct {
    char *name;
    void *fun;
} pair;
#endif
*/

#include "../common/common.h"

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

static const keyval kv_specfun[] = {
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
