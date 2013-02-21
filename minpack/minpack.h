
typedef struct {
    int n;
    int m;
    double *P;
    int *do_var;
    double *x;
    double *y;
    double *ydata;
    double *w;
    double *a;
} data;

typedef void (func)(data *D);

int leastsq(func *f, data *D);

int leastsq1(func *f, data *D);


