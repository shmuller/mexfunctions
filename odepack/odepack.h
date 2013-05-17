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

struct _data {
    odefunf_t *wrapper;       // wrapper with Fortran call signature
    odefun_t *dy_dt;          // function calculating the RHS of the ode
    void *args;               // any additional arguments
    int neq;                  // number of equations (length of y, ydot)
    int n;                    // number of time points (length of time, res)
    double *time;             // full time vector
    double *res;              // full results vector
    double *t;                // current time
    double *y;                // current value
    double *ydot;             // vector in which to return current RHS
    void *jac;                // function calculating the jacobian
    odefun_t *term;           // function checking termination conditions
    int points_done;          // number of time points done
};

static const data empty_data = {0};

int odesolve(data *D);


