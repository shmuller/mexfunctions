#ifndef __QUADFUN_H__
#define __QUADFUN_H__

#ifdef __cplusplus
extern "C"
{
#endif

typedef double (func)(double *);

typedef struct {
    double *x;
    int N;
    int o;
} link;

void quadfun(func *fun, link *LI, int nI, int Np, int Ns, int N, double *y);

#ifdef __cplusplus
}
#endif

#endif /* __QUADFUN_H__ */
