/* odepack.h
 *
 * Wrap odepack Fortran subroutine dlsode to the 
 * (func*, data*) interface.
 *
 * S. H. Muller, 2013/05/15
 */


typedef struct _data data;

typedef int (odefun_t)(data *D);
typedef void (odefunf_t)(int *neq, double *t, double *y, double *ydot);
typedef void (rootfunf_t)(int *neq, double *t, double *y, int *ng, double *gout);

struct _data {
    odefunf_t *f;             // wrapper with Fortran call signature
    int neq;                  // number of equations (length of y, ydot)
    int n;                    // number of time points (length of time, res)
    double *t;                // full time vector
    double *y;                // full results vector
    void *jac;                // function calculating the jacobian
    odefun_t *term;           // function checking termination conditions
    int points_done;          // number of time points done
    int ng;                   // number of constraints in g
    rootfunf_t *g;            // constraint function
};

static const data empty_data = {0};

int odesolve(data *D);


