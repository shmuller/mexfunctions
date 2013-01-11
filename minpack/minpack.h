
typedef struct {
    int n;
    int m;
    double *P;
    double *x;
    double *y;
    double *ydata;
    double *w;
    double *a;
} data;

typedef void (func)(data *D);

int leastsq(func *f, data *D);


