typedef struct {
    void *splR;
    void *splz;
    void *D;
    double *bdry;
    int n_bdry;
} fieldline_data;

void solve_bdry(fieldline_data *F);

