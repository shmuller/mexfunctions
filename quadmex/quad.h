#ifndef __QUAD_H__
#define __QUAD_H__

#ifdef __cplusplus
extern "C"
{
#endif

void looppar(const int rank, const int *N, int *s, int *no, int *ni, const int dim);
void filldim(int rank, const int *N, double *X, const double *x, int dim);
void filldim2(int rank, const int *N, double *X, const double *x, int dim);
void weighted_sum(double *y, const double *w, int n, int ni);
double *scaleX(int n, double *x, double A, double B);
double getScaledX(int N, int *L, double *X, double *ab, int *dtbl, double **x, double **w, int d);

#ifdef __cplusplus
}
#endif

#endif /* __QUAD_H__ */
