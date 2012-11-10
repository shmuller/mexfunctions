/* Interface definition for specfun.f90.
 *
 * S. H. Muller, 2012/01/12
 */

#ifdef __cplusplus
extern "C"
{
#endif

#include "../common/common.h"

#define besei0 besei0_
#define besei1 besei1_
#define besek0 besek0_
#define besek1 besek1_
#define besi0 besi0_
#define besi1 besi1_
#define besj0 besj0_
#define besj1 besj1_
#define besk0 besk0_
#define besk1 besk1_
#define besy0 besy0_
#define besy1 besy1_
#define daw daw_
#define dlgama dlgama_
#define ei ei_
#define eone eone_
#define expei expei_
#define r8_erf r8_erf_
#define r8_erfc r8_erfc_
#define r8_erfcx r8_erfcx_
#define r8_gamma r8_gamma_
#define r8_psi r8_psi_

extern func besei0_;   // exponentially scaled Bessel I0(X) function
extern func besei1_;   // exponentially scaled Bessel I1(X) function
extern func besek0_;   // exponentially scaled Bessel K0(X) function
extern func besek1_;   // exponentially scaled Bessel K1(X) function
extern func besi0_;    // Bessel I0(X) function
extern func besi1_;    // Bessel I1(X) function
extern func besj0_;    // Bessel J0(X) function
extern func besj1_;    // Bessel J1(X) function
extern func besk0_;    // Bessel K0(X) function
extern func besk1_;    // Bessel K1(X) function
extern func besy0_;    // Bessel Y0(X) function
extern func besy1_;    // Bessel Y1(X) function
extern func daw_;      // Dawson's integral function
extern func dlgama_;   // log ( Gamma ( X ) ) for a real argument
extern func ei_;       // exponential integral Ei(X)
extern func eone_;     // exponential integral E1(X)
extern func expei_;    // scaled exponential integral exp(-X) * Ei(X)
extern func r8_erf_;   // error function
extern func r8_erfc_;  // complementary error function
extern func r8_erfcx_; // exponentially scaled complementary error function
extern func r8_gamma_; // Gamma(X) for a real argument
extern func r8_psi_;   // Psi(X)

static const keyval kv_specfun[] = {
    "besei0", (const void*) besei0_,
    "besei1", (const void*) besei1_,
    "besek0", (const void*) besek0_,
    "besek1", (const void*) besek1_,
    "besi0",  (const void*) besi0_,
    "besi1",  (const void*) besi1_,
    "besj0",  (const void*) besj0_,
    "besj1",  (const void*) besj1_,
    "besk0",  (const void*) besk0_,
    "besk1",  (const void*) besk1_,
    "besy0",  (const void*) besy0_,
    "besy1",  (const void*) besy1_,
    "daw",    (const void*) daw_,
    "dlgama", (const void*) dlgama_,
    "ei",     (const void*) ei_,
    "eone",   (const void*) eone_,
    "expei",  (const void*) expei_,
    "erf",    (const void*) r8_erf_,
    "erfc",   (const void*) r8_erfc_,
    "erfcx",  (const void*) r8_erfcx_,
    "gamma",  (const void*) r8_gamma_,
    "psi",    (const void*) r8_psi_
};

#ifdef __cplusplus
}  /* end extern "C" */
#endif


