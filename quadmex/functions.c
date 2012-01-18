
typedef struct {
    char *name;
    func *fun;
} pair;


double Fun(double *par)
{
    double x=par[0], y=par[1], z=par[2];
    return x*y*y*z*z*z;
}


static const pair P[] = {
    "Fun", Fun
};
