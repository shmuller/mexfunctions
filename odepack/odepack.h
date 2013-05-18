/* odepack.h
 *
 * Wrap odepack Fortran subroutine dlsode to the 
 * (func*, data*) interface.
 *
 * S. H. Muller, 2013/05/15
 */


typedef struct _data data;

typedef void (odefun_t)(int *neq, double *t, double *y, double *ydot);
typedef void (oderoot_t)(int *neq, double *t, double *y, int *ng, double *gout);

struct _data {
    odefun_t *f;              // wrapper with Fortran call signature
    int neq;                  // number of equations (length of y, ydot)
    int n;                    // number of time points (length of time, res)
    double *t;                // full time vector
    double *y;                // full results vector
    int points_done;          // number of time points done
    int ng;                   // number of constraints in g
    oderoot_t *g;             // constraint function
    int *ibbox;               // indices of y for which bbox is checked (length ng)
    double *bbox;             // bounding box (length ng)
};

static const data empty_data = {0};

oderoot_t bbox_g;

void odesolve(data *D);


