/* odepack.h
 *
 * Wrap odepack Fortran subroutine dlsode to the 
 * (func*, data*) interface.
 *
 * S. H. Muller, 2013/05/15
 */


typedef struct _data data;

struct _data {
    void (*dy_dt)(data*);   // function calculating the RHS of the ode
    int neq;                // number of equations (length of y, ydot)
    double *time;           // full time vector
    double *res;            // full results vector
    double *t;              // current time
    double *y;              // current time
    double *ydot;           // vector in which to return current RHS
    void (*jac)(data*);     // function calculating the jacobian
    void (*term)(data*);    // function checking termination conditions
};

int odeint(data *D);


