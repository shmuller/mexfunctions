/* Interface definition for specfun.f90.
 *
 * S. H. Muller, 2012/01/12
 */

#ifdef __cplusplus
extern "C"
{
#endif

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
    "besei0", (const void*) besei0,
    "besei1", (const void*) besei1,
    "besek0", (const void*) besek0,
    "besek1", (const void*) besek1,
    "besi0",  (const void*) besi0,
    "besi1",  (const void*) besi1,
    "besj0",  (const void*) besj0,
    "besj1",  (const void*) besj1,
    "besk0",  (const void*) besk0,
    "besk1",  (const void*) besk1,
    "besy0",  (const void*) besy0,
    "besy1",  (const void*) besy1,
    "daw",    (const void*) daw,
    "dlgama", (const void*) dlgama,
    "ei",     (const void*) ei,
    "eone",   (const void*) eone,
    "expei",  (const void*) expei,
    "erf",    (const void*) r8_erf,
    "erfc",   (const void*) r8_erfc,
    "erfcx",  (const void*) r8_erfcx,
    "gamma",  (const void*) r8_gamma,
    "psi",    (const void*) r8_psi
};

#ifdef __cplusplus
}  /* end extern "C" */
#endif


