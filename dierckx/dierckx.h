typedef struct {
    double *tx;
    int nx;
    double *ty;
    int ny;
    double *c;
    int kx;
    int ky;
    double *x;
    int mx;
    double *y;
    int my;
    double *z;
    double *wrk;
    int lwrk;
    int *iwrk;
    int kwrk;
    int ier;
} dierckx_data;

void *dierckx_make_workspace(dierckx_data *D);
void dierckx_bispev(dierckx_data *D);


