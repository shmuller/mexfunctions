#ifndef __QUADFUN_H__
#define __QUADFUN_H__

#ifdef __cplusplus
extern "C"
{
#endif

#ifndef __COMMON__
#define __COMMON__
typedef double real;

typedef real (func)(real*);

typedef struct {
    char *name;
    void *fun;
} pair;
#endif /* __COMMON__ */

typedef struct {
    real *x;
    int N;
    int o;
} link;

void quadfun(func *fun, link *LI, int nI, int Np, int Ns, int N, real *y);

#ifdef __cplusplus
}
#endif

#endif /* __QUADFUN_H__ */
