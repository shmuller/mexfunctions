/* minpack.h
 *
 * Wrap minpack Fortran subroutines lmdif and lmdif1 to the 
 * (func*, data*) interface.
 *
 * S. H. Muller, 2012
 */

typedef struct {
    int n;          // length of P, do_var
    int m;          // length of x, y, ydata, w, a
    double *P;      // parameter vector to be optimized
    int *do_var;    // bool array to switch off optimization for some parameters
    double *x;      // independent coordinate of data sites
    double *y;      // vector to hold the result of function evaluations
    double *ydata;  // data values
    double *w;      // weights (not implemented)
    double *a;      // parameter for fits allowing linear variations of P
} data;

typedef void (func)(data *D);

int leastsq(func *f, data *D);

int leastsq1(func *f, data *D);


