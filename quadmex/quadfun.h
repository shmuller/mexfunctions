#ifndef __QUADFUN_H__
#define __QUADFUN_H__

#ifdef __cplusplus
extern "C"
{
#endif

#include "../common/common.h"

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
